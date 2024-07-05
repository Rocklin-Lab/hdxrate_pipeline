import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from hxrate.backexchange import calc_back_exchange
from scipy import odr
import argparse
from hxrate.methods import normalize_mass_distribution_array, gauss_fit_to_isotope_dist_array, \
    correct_centroids_using_backexchange
from hxrate.hxdata import load_data_from_hdx_ms_dist, load_tp_dependent_dict


def gen_list_of_backexchange(list_of_dist_files,
                             list_of_seq,
                             d2o_frac,
                             d2o_purity):
    """

    :param list_of_dist_files:
    :param list_of_seq:
    :param d2o_frac:
    :param d2o_purity:
    :return:
    """

    backexchange_list = []

    for ind, (fpath, seq) in enumerate(zip(list_of_dist_files, list_of_seq)):

        data_dict = load_data_from_hdx_ms_dist(fpath=fpath)
        tp = data_dict['tp']
        dist = data_dict['mass_dist']
        norm_dist = normalize_mass_distribution_array(mass_dist_array=dist)

        bkexch_obj = calc_back_exchange(sequence=seq,
                                        experimental_isotope_dist=norm_dist[-1],
                                        timepoints_array=tp,
                                        d2o_fraction=d2o_frac,
                                        d2o_purity=d2o_purity)
        backexchange_list.append(bkexch_obj.backexchange_value)

    return backexchange_list


def plot_backexchange_correlation(low_ph_corr_bkexch,
                                  high_ph_corr_bkexch,
                                  low_ph_nocorr_backexch,
                                  high_ph_nocorr_backexch,
                                  corr_slope,
                                  corr_intercept,
                                  corr_r,
                                  output_path):
    """
    plot backexchange correlation with low and high ph data
    :param low_ph_corr_bkexch:
    :param high_ph_corr_bkexch:
    :param low_ph_saturate_backexch:
    :param high_ph_saturate_backexch:
    :param low_ph_others:
    :param high_ph_others:
    :param corr_slope:
    :param corr_intercept:
    :param corr_r:
    :param output_path:
    :return:
    """

    x_ = np.linspace(start=min(low_ph_corr_bkexch) - 0.05, stop=max(low_ph_corr_bkexch) + 0.05, num=20)
    y_ = corr_slope * x_ + corr_intercept

    title = 'Slope: %.4f | Intercept: %.4f | R: %.4f' % (corr_slope, corr_intercept, corr_r)

    plt.scatter(low_ph_nocorr_backexch, high_ph_nocorr_backexch, color='aquamarine', label='unsatisfied')
    plt.scatter(low_ph_corr_bkexch, high_ph_corr_bkexch, color='darkcyan', label='satisfied')
    plt.plot(x_, y_, ls='--', color='darkcyan')
    plt.xlabel('low ph d2o saturation level')
    plt.ylabel('high ph d2o saturation level')
    plt.grid()
    plt.title(title)
    plt.savefig(output_path)
    plt.close()


# def check_for_exchange_saturation(centroid_list, mass_rate_threshold=0.03):
#     """

#     :param dist_list:
#     :param mass_rate_threshold:
#     :return:
#     """

#     mass_rate = abs((centroid_list[1] - centroid_list[0])/centroid_list[0])
#     if mass_rate > mass_rate_threshold:
#         return False
#     else:
#         return True

def check_for_exchange_saturation(centroid_list, mass_rate_threshold=1):
    """

    :param dist_list:
    :param mass_rate_threshold:
    :return:
    """

    mass_rate = abs(centroid_list[1] - centroid_list[0])
    if mass_rate > mass_rate_threshold:
        return False
    else:
        return True


def check_saturation_from_file(fpath,
                               sequence,
                               d2o_frac,
                               d2o_pur,
                               mass_rate_threshold=0.03,
                               dist_indices=[-1, -5],
                               bkexch_corr_fpath=None):
    """

    :param fpath:
    :param mass_rate_threshold:
    :param dist_indices:
    :return:
    """

    data_dict = load_data_from_hdx_ms_dist(fpath=fpath)
    tp = data_dict['tp']
    dist_list = data_dict['mass_dist']
    norm_dist_list = normalize_mass_distribution_array(mass_dist_array=dist_list)

    # use backexchange correction to correct the centroids and then check for saturation if bkexch corr filepath given
    if bkexch_corr_fpath is not None:
        bkexch_corr_dict = load_tp_dependent_dict(bkexch_corr_fpath)
        bkexchange_object = calc_back_exchange(sequence=sequence,
                                               experimental_isotope_dist=norm_dist_list[-1],
                                               timepoints_array=tp,
                                               d2o_fraction=d2o_frac,
                                               d2o_purity=d2o_pur,
                                               backexchange_corr_dict=bkexch_corr_dict)
    else:
        bkexchange_object = calc_back_exchange(sequence=sequence,
                                               experimental_isotope_dist=norm_dist_list[-1],
                                               timepoints_array=tp,
                                               d2o_fraction=d2o_frac,
                                               d2o_purity=d2o_pur)
    if len(norm_dist_list) > 2:
        dist_to_check = norm_dist_list[dist_indices]
        gauss_fit_list = gauss_fit_to_isotope_dist_array(isotope_dist=dist_to_check)
        centroid_list = [x.centroid for x in gauss_fit_list]
        backexchange_for_check = bkexchange_object.backexchange_array[dist_indices]
        corr_centroid_list = correct_centroids_using_backexchange(centroids=np.array(centroid_list),
                                                                          backexchange_array=backexchange_for_check)

        satur_bool = check_for_exchange_saturation(centroid_list=corr_centroid_list,
                                                   mass_rate_threshold=mass_rate_threshold)
    else:
        satur_bool = False

    outdict = dict()
    outdict['backexchange_value'] = bkexchange_object.backexchange_value
    outdict['saturation_bool'] = satur_bool

    return outdict


def gen_list_of_bkexch_and_saturation_bool(filepath_list,
                                           sequence_list,
                                           d2o_fraction,
                                           d2o_purity,
                                           mass_rate_threshold,
                                           dist_indices,
                                           bkexch_corr_fpath):

    # gen backexchange and check for saturation
    backexchange_arr = np.zeros(len(filepath_list))
    check_saturation_bool_list = []
    for ind, (fpath, seq) in enumerate(zip(filepath_list, sequence_list)):
        outdict = check_saturation_from_file(fpath=fpath,
                                             sequence=seq,
                                             d2o_frac=d2o_fraction,
                                             d2o_pur=d2o_purity,
                                             mass_rate_threshold=mass_rate_threshold,
                                             dist_indices=dist_indices,
                                             bkexch_corr_fpath=bkexch_corr_fpath)
        backexchange_arr[ind] = outdict['backexchange_value']
        check_saturation_bool_list.append(outdict['saturation_bool'])

    out_dict = dict()
    out_dict['backexchange_array'] = backexchange_arr
    out_dict['saturation_bool_list'] = check_saturation_bool_list

    return out_dict

#
# def gen_list_of_bkexch_and_saturation_bool_from_list_of_bkexch_files(bkexch_file_list,
#                                                                      )


def write_low_high_backexchange_array(low_ph_protein_name,
                                      high_ph_protein_name,
                                      sequence_array,
                                      low_ph_backexchange_array,
                                      high_ph_backexchange_array,
                                      low_ph_saturation,
                                      high_ph_saturation,
                                      corr_include_indices,
                                      low_ph_backexchange_new,
                                      low_ph_backexchange_all_corr,
                                      output_path):
    """

    :param low_ph_protein_name:
    :param high_ph_protein_name:
    :param low_ph_backexchange_array:
    :param high_ph_backexchange_array:
    :param corr_include_indices:
    :param output_path:
    :return:
    """

    corr_include_arr = np.zeros(len(low_ph_backexchange_array))
    corr_include_arr[corr_include_indices] = 1

    header = 'low_ph_protein_name,high_ph_protein_name,low_high_name,sequence,low_ph_backexchange,high_ph_backexchange,low_ph_saturation,high_ph_saturation,corr_include,low_ph_backexchange_new,low_ph_backexchange_all_corr\n'
    data_string = ''

    for num in range(len(corr_include_arr)):

        data_string += '{},{},{},{},{},{},{},{},{},{},{}\n'.format(low_ph_protein_name[num],
                                                                   high_ph_protein_name[num],
                                                                   low_ph_protein_name[num] + '_' + high_ph_protein_name[num],
                                                                   sequence_array[num],
                                                                   low_ph_backexchange_array[num],
                                                                   high_ph_backexchange_array[num],
                                                                   low_ph_saturation[num],
                                                                   high_ph_saturation[num],
                                                                   corr_include_arr[num],
                                                                   low_ph_backexchange_new[num],
                                                                   low_ph_backexchange_all_corr[num]
                                                                   )

    with open(output_path, 'w') as outfile:
        outfile.write(header + data_string)
        outfile.close()


def write_low_high_backexchange_correlation(slope,
                                            intercept,
                                            rvalue,
                                            std_slope_error,
                                            std_intercept_error,
                                            output_path):

    header = 'slope,intercept,rvalue,std_slope_error,std_intercept_error\n'
    data_string = '{},{},{},{},{}\n'.format(slope,
                                               intercept,
                                               rvalue,
                                               std_slope_error,
                                               std_intercept_error)

    with open(output_path, 'w') as outfile:
        outfile.write(header + data_string)
        outfile.close()


def gen_low_high_ph_bkexchange_(merge_sample_list_fpath: str,
                                low_ph_d2o_frac: float,
                                low_ph_d2o_pur: float,
                                high_ph_d2o_frac: float,
                                high_ph_d2o_pur: float,
                                bkexchange_ub: float,
                                bkexchange_lb: float,
                                saturation_mass_rate_threshold: float,
                                lowph_backexchange_correction_filepath: str,
                                highph_backexchange_correction_filepath: str,
                                plot_output_path: str = None,
                                bkexchange_corr_output_path: str = None,
                                bkexchange_output_path: str = None,
                                return_flag: bool = False):
    """

    :param merge_sample_list_fpath:
    :param low_ph_d2o_frac:
    :param low_ph_d2o_pur:
    :param high_ph_d2o_frac:
    :param high_ph_d2o_pur:
    :param bkexchange_ub:
    :param bkexchange_lb:
    :param plot_output_path:
    :param bkexchange_corr_output_path:
    :param bkexchange_output_path:
    :param return_flag:
    :return:
    """

    merge_sample_df = pd.read_csv(merge_sample_list_fpath)

    seq_array = merge_sample_df['sequence'].values

    high_ph_fpath_list = merge_sample_df['high_ph_data_fpath'].values
    low_ph_fpath_list = merge_sample_df['low_ph_data_fpath'].values

    low_ph_bkexch_satur = gen_list_of_bkexch_and_saturation_bool(filepath_list=low_ph_fpath_list,
                                                                 sequence_list=seq_array,
                                                                 d2o_fraction=low_ph_d2o_frac,
                                                                 d2o_purity=low_ph_d2o_pur,
                                                                 mass_rate_threshold=saturation_mass_rate_threshold,
                                                                 dist_indices=[-1, -5],
                                                                 bkexch_corr_fpath=lowph_backexchange_correction_filepath)

    high_ph_bkexch_satur = gen_list_of_bkexch_and_saturation_bool(filepath_list=high_ph_fpath_list,
                                                                  sequence_list=seq_array,
                                                                  d2o_fraction=high_ph_d2o_frac,
                                                                  d2o_purity=high_ph_d2o_pur,
                                                                  mass_rate_threshold=saturation_mass_rate_threshold,
                                                                  dist_indices=[-1, -5],
                                                                  bkexch_corr_fpath=highph_backexchange_correction_filepath)

    check_saturation_bool = []
    for bool_low, bool_high in zip(low_ph_bkexch_satur['saturation_bool_list'], high_ph_bkexch_satur['saturation_bool_list']):
        if bool_low:
            if bool_high:
                check_saturation_bool.append(True)
            else:
                check_saturation_bool.append(False)
        else:
            check_saturation_bool.append(False)

    corr_include_indices = []
    no_corr_indices = []
    for ind, (low_ph_bkexch_, high_ph_bkexch_) in enumerate(zip(low_ph_bkexch_satur['backexchange_array'],
                                                                high_ph_bkexch_satur['backexchange_array'])):
        if check_saturation_bool[ind]:
            if bkexchange_lb <= low_ph_bkexch_ <= bkexchange_ub:
                if bkexchange_lb <= high_ph_bkexch_ <= bkexchange_ub:
                    corr_include_indices.append(ind)
                else:
                    no_corr_indices.append(ind)
            else:
                no_corr_indices.append(ind)
        else:
            no_corr_indices.append(ind)

    corr_include_indices = np.array(corr_include_indices)
    no_corr_indices = np.array(no_corr_indices)

    low_ph_bkexch_corr = low_ph_bkexch_satur['backexchange_array'][corr_include_indices]
    high_ph_bkexch_corr = high_ph_bkexch_satur['backexchange_array'][corr_include_indices]

    low_ph_bkexch_nocorr = low_ph_bkexch_satur['backexchange_array'][no_corr_indices]
    high_ph_bkexch_nocorr = high_ph_bkexch_satur['backexchange_array'][no_corr_indices]

    # use orthogonal distance regression for low and high ph correlation
    data = odr.Data(x=low_ph_bkexch_corr, y=high_ph_bkexch_corr)
    odr_object = odr.ODR(data=data, model=odr.unilinear)
    odr_output = odr_object.run()

    corr_coef = np.corrcoef(x=low_ph_bkexch_corr,
                            y=high_ph_bkexch_corr)

    # generate low ph backexchange if saturation not achieved
    low_ph_bkexch_new_arr = np.zeros(len(low_ph_bkexch_satur['backexchange_array']))
    low_ph_bkexch_all_corr_arr = np.zeros(len(low_ph_bkexch_satur['backexchange_array']))

    for ind, (low_ph_bkexch, high_ph_bkexch, low_ph_satur) in enumerate(zip(low_ph_bkexch_satur['backexchange_array'],
                                                                            high_ph_bkexch_satur['backexchange_array'],
                                                                            low_ph_bkexch_satur['saturation_bool_list'])):

        corr_low_ph_bkexch = (high_ph_bkexch - odr_output.beta[1])/odr_output.beta[0]
        low_ph_bkexch_all_corr_arr[ind] = corr_low_ph_bkexch
        if low_ph_satur:
            low_ph_bkexch_new_arr[ind] = low_ph_bkexch
        else:
            low_ph_bkexch_new_arr[ind] = corr_low_ph_bkexch

    plot_backexchange_correlation(low_ph_bkexch_corr,
                                  high_ph_bkexch_corr,
                                  low_ph_nocorr_backexch=low_ph_bkexch_nocorr,
                                  high_ph_nocorr_backexch=high_ph_bkexch_nocorr,
                                  corr_slope=odr_output.beta[0],
                                  corr_intercept=odr_output.beta[1],
                                  corr_r=corr_coef[0, 1],
                                  output_path=plot_output_path)

    write_low_high_backexchange_correlation(slope=odr_output.beta[0],
                                            intercept=odr_output.beta[1],
                                            rvalue=corr_coef[0, 1],
                                            std_slope_error=odr_output.sd_beta[0],
                                            std_intercept_error=odr_output.sd_beta[1],
                                            output_path=bkexchange_corr_output_path)

    write_low_high_backexchange_array(low_ph_protein_name=merge_sample_df['protein_name_low_ph'].values,
                                      high_ph_protein_name=merge_sample_df['protein_name_high_ph'].values,
                                      sequence_array=seq_array,
                                      low_ph_backexchange_array=low_ph_bkexch_satur['backexchange_array'],
                                      high_ph_backexchange_array=high_ph_bkexch_satur['backexchange_array'],
                                      low_ph_saturation=low_ph_bkexch_satur['saturation_bool_list'],
                                      high_ph_saturation=high_ph_bkexch_satur['saturation_bool_list'],
                                      corr_include_indices=corr_include_indices,
                                      low_ph_backexchange_new=low_ph_bkexch_new_arr,
                                      low_ph_backexchange_all_corr=low_ph_bkexch_all_corr_arr,
                                      output_path=bkexchange_output_path)

    if return_flag:
        out_dict = dict()
        out_dict['low_ph_backexchange_array'] = low_ph_bkexch_satur['backexchange_array']
        out_dict['high_ph_backexchange_array'] = high_ph_bkexch_satur['backexchange_array']
        out_dict['low_ph_saturation_bool_list'] = low_ph_bkexch_satur['saturation_bool_list']
        out_dict['high_ph_saturation_bool_list'] = high_ph_bkexch_satur['saturation_bool_list']
        return out_dict


def gen_low_high_ph_bkexchange_from_file_list(list_of_low_ph_hxms_files: list,
                                              list_of_high_ph_hxms_files: list,
                                              low_ph_protname_list: list,
                                              high_ph_protname_list: list,
                                              sequence_list: list,
                                              low_ph_d2o_frac: float,
                                              low_ph_d2o_pur: float,
                                              high_ph_d2o_frac: float,
                                              high_ph_d2o_pur: float,
                                              bkexchange_ub: float,
                                              bkexchange_lb: float,
                                              saturation_mass_rate_threshold: float,
                                              lowph_backexchange_correction_filepath: str,
                                              highph_backexchange_correction_filepath: str,
                                              plot_output_path: str = None,
                                              bkexchange_corr_output_path: str = None,
                                              bkexchange_output_path: str = None,
                                              return_flag: bool = False):
    """

    :param merge_sample_list_fpath:
    :param low_ph_d2o_frac:
    :param low_ph_d2o_pur:
    :param high_ph_d2o_frac:
    :param high_ph_d2o_pur:
    :param bkexchange_ub:
    :param bkexchange_lb:
    :param plot_output_path:
    :param bkexchange_corr_output_path:
    :param bkexchange_output_path:
    :param return_flag:
    :return:
    """

    low_ph_bkexch_satur = gen_list_of_bkexch_and_saturation_bool(filepath_list=list_of_low_ph_hxms_files,
                                                                 sequence_list=sequence_list,
                                                                 d2o_fraction=low_ph_d2o_frac,
                                                                 d2o_purity=low_ph_d2o_pur,
                                                                 mass_rate_threshold=saturation_mass_rate_threshold,
                                                                 dist_indices=[-1, -5],
                                                                 bkexch_corr_fpath=lowph_backexchange_correction_filepath)

    high_ph_bkexch_satur = gen_list_of_bkexch_and_saturation_bool(filepath_list=list_of_high_ph_hxms_files,
                                                                  sequence_list=sequence_list,
                                                                  d2o_fraction=high_ph_d2o_frac,
                                                                  d2o_purity=high_ph_d2o_pur,
                                                                  mass_rate_threshold=saturation_mass_rate_threshold,
                                                                  dist_indices=[-1, -5],
                                                                  bkexch_corr_fpath=highph_backexchange_correction_filepath)

    check_saturation_bool = []
    for bool_low, bool_high in zip(low_ph_bkexch_satur['saturation_bool_list'], high_ph_bkexch_satur['saturation_bool_list']):
        if bool_low:
            if bool_high:
                check_saturation_bool.append(True)
            else:
                check_saturation_bool.append(False)
        else:
            check_saturation_bool.append(False)

    corr_include_indices = []
    no_corr_indices = []
    for ind, (low_ph_bkexch_, high_ph_bkexch_) in enumerate(zip(low_ph_bkexch_satur['backexchange_array'],
                                                                high_ph_bkexch_satur['backexchange_array'])):
        if check_saturation_bool[ind]:
            if bkexchange_lb <= low_ph_bkexch_ <= bkexchange_ub:
                if bkexchange_lb <= high_ph_bkexch_ <= bkexchange_ub:
                    corr_include_indices.append(ind)
                else:
                    no_corr_indices.append(ind)
            else:
                no_corr_indices.append(ind)
        else:
            no_corr_indices.append(ind)

    corr_include_indices = np.array(corr_include_indices)
    no_corr_indices = np.array(no_corr_indices)

    low_ph_bkexch_corr = low_ph_bkexch_satur['backexchange_array'][corr_include_indices]
    high_ph_bkexch_corr = high_ph_bkexch_satur['backexchange_array'][corr_include_indices]

    low_ph_bkexch_nocorr = low_ph_bkexch_satur['backexchange_array'][no_corr_indices]
    high_ph_bkexch_nocorr = high_ph_bkexch_satur['backexchange_array'][no_corr_indices]

    # use orthogonal distance regression for low and high ph correlation
    data = odr.Data(x=low_ph_bkexch_corr, y=high_ph_bkexch_corr)
    odr_object = odr.ODR(data=data, model=odr.unilinear)
    odr_output = odr_object.run()

    corr_coef = np.corrcoef(x=low_ph_bkexch_corr,
                            y=high_ph_bkexch_corr)

    # generate low ph backexchange if saturation not achieved
    low_ph_bkexch_new_arr = np.zeros(len(low_ph_bkexch_satur['backexchange_array']))
    low_ph_bkexch_corr_arr = np.zeros(len(low_ph_bkexch_satur['backexchange_array']))
    for ind, (low_ph_bkexch, high_ph_bkexch, low_ph_satur) in enumerate(zip(low_ph_bkexch_satur['backexchange_array'],
                                                                            high_ph_bkexch_satur['backexchange_array'],
                                                                            low_ph_bkexch_satur['saturation_bool_list'])):

        low_ph_corr_bkexch = (high_ph_bkexch - odr_output.beta[1])/odr_output.beta[0]
        low_ph_bkexch_corr_arr[ind] = low_ph_corr_bkexch
        if low_ph_satur:
            low_ph_bkexch_new_arr[ind] = low_ph_bkexch
        else:
            low_ph_bkexch_new_arr[ind] = low_ph_corr_bkexch

    plot_backexchange_correlation(low_ph_bkexch_corr,
                                  high_ph_bkexch_corr,
                                  low_ph_nocorr_backexch=low_ph_bkexch_nocorr,
                                  high_ph_nocorr_backexch=high_ph_bkexch_nocorr,
                                  corr_slope=odr_output.beta[0],
                                  corr_intercept=odr_output.beta[1],
                                  corr_r=corr_coef[0, 1],
                                  output_path=plot_output_path)

    write_low_high_backexchange_correlation(slope=odr_output.beta[0],
                                            intercept=odr_output.beta[1],
                                            rvalue=corr_coef[0, 1],
                                            std_slope_error=odr_output.sd_beta[0],
                                            std_intercept_error=odr_output.sd_beta[1],
                                            output_path=bkexchange_corr_output_path)

    write_low_high_backexchange_array(low_ph_protein_name=low_ph_protname_list,
                                      high_ph_protein_name=high_ph_protname_list,
                                      sequence_array=sequence_list,
                                      low_ph_backexchange_array=low_ph_bkexch_satur['backexchange_array'],
                                      high_ph_backexchange_array=high_ph_bkexch_satur['backexchange_array'],
                                      low_ph_saturation=low_ph_bkexch_satur['saturation_bool_list'],
                                      high_ph_saturation=high_ph_bkexch_satur['saturation_bool_list'],
                                      corr_include_indices=corr_include_indices,
                                      low_ph_backexchange_new=low_ph_bkexch_new_arr,
                                      low_ph_backexchange_all_corr=low_ph_bkexch_corr_arr,
                                      output_path=bkexchange_output_path)

    if return_flag:
        out_dict = dict()
        out_dict['low_ph_backexchange_array'] = low_ph_bkexch_satur['backexchange_array']
        out_dict['high_ph_backexchange_array'] = high_ph_bkexch_satur['backexchange_array']
        out_dict['low_ph_saturation_bool_list'] = low_ph_bkexch_satur['saturation_bool_list']
        out_dict['high_ph_saturation_bool_list'] = high_ph_bkexch_satur['saturation_bool_list']
        return out_dict
#
#
# def gen_low_high_ph_bkexchange_from_bkexch_file_list(list_of_low_ph_bkexch_files: list,
#                                                      list_of_high_ph_bkexch_files: list,
#                                                      low_ph_protname_list: list,
#                                                      high_ph_protname_list: list,
#                                                      sequence_list: list,
#                                                      low_ph_d2o_frac: float,
#                                                      low_ph_d2o_pur: float,
#                                                      high_ph_d2o_frac: float,
#                                                      high_ph_d2o_pur: float,
#                                                      bkexchange_ub: float,
#                                                      bkexchange_lb: float,
#                                                      saturation_mass_rate_threshold: float,
#                                                      lowph_backexchange_correction_filepath: str,
#                                                      highph_backexchange_correction_filepath: str,
#                                                      plot_output_path: str = None,
#                                                      bkexchange_corr_output_path: str = None,
#                                                      bkexchange_output_path: str = None,
#                                                      return_flag: bool = False):
#     """
#
#     :param merge_sample_list_fpath:
#     :param low_ph_d2o_frac:
#     :param low_ph_d2o_pur:
#     :param high_ph_d2o_frac:
#     :param high_ph_d2o_pur:
#     :param bkexchange_ub:
#     :param bkexchange_lb:
#     :param plot_output_path:
#     :param bkexchange_corr_output_path:
#     :param bkexchange_output_path:
#     :param return_flag:
#     :return:
#     """
#
#     low_ph_bkexch_satur = gen_list_of_bkexch_and_saturation_bool(filepath_list=list_of_low_ph_hxms_files,
#                                                                  sequence_list=sequence_list,
#                                                                  d2o_fraction=low_ph_d2o_frac,
#                                                                  d2o_purity=low_ph_d2o_pur,
#                                                                  mass_rate_threshold=saturation_mass_rate_threshold,
#                                                                  dist_indices=[-1, -3],
#                                                                  bkexch_corr_fpath=lowph_backexchange_correction_filepath)
#
#     high_ph_bkexch_satur = gen_list_of_bkexch_and_saturation_bool(filepath_list=list_of_high_ph_hxms_files,
#                                                                   sequence_list=sequence_list,
#                                                                   d2o_fraction=high_ph_d2o_frac,
#                                                                   d2o_purity=high_ph_d2o_pur,
#                                                                   mass_rate_threshold=saturation_mass_rate_threshold,
#                                                                   dist_indices=[-1, -3],
#                                                                   bkexch_corr_fpath=highph_backexchange_correction_filepath)
#
#     check_saturation_bool = []
#     for bool_low, bool_high in zip(low_ph_bkexch_satur['saturation_bool_list'], high_ph_bkexch_satur['saturation_bool_list']):
#         if bool_low:
#             if bool_high:
#                 check_saturation_bool.append(True)
#             else:
#                 check_saturation_bool.append(False)
#         else:
#             check_saturation_bool.append(False)
#
#     corr_include_indices = []
#     no_corr_indices = []
#     for ind, (low_ph_bkexch_, high_ph_bkexch_) in enumerate(zip(low_ph_bkexch_satur['backexchange_array'],
#                                                                 high_ph_bkexch_satur['backexchange_array'])):
#         if check_saturation_bool[ind]:
#             if bkexchange_lb <= low_ph_bkexch_ <= bkexchange_ub:
#                 if bkexchange_lb <= high_ph_bkexch_ <= bkexchange_ub:
#                     corr_include_indices.append(ind)
#                 else:
#                     no_corr_indices.append(ind)
#             else:
#                 no_corr_indices.append(ind)
#         else:
#             no_corr_indices.append(ind)
#
#     corr_include_indices = np.array(corr_include_indices)
#     no_corr_indices = np.array(no_corr_indices)
#
#     low_ph_bkexch_corr = low_ph_bkexch_satur['backexchange_array'][corr_include_indices]
#     high_ph_bkexch_corr = high_ph_bkexch_satur['backexchange_array'][corr_include_indices]
#
#     low_ph_bkexch_nocorr = low_ph_bkexch_satur['backexchange_array'][no_corr_indices]
#     high_ph_bkexch_nocorr = high_ph_bkexch_satur['backexchange_array'][no_corr_indices]
#
#     # use orthogonal distance regression for low and high ph correlation
#     data = odr.Data(x=low_ph_bkexch_corr, y=high_ph_bkexch_corr)
#     odr_object = odr.ODR(data=data, model=odr.unilinear)
#     odr_output = odr_object.run()
#
#     corr_coef = np.corrcoef(x=low_ph_bkexch_corr,
#                             y=high_ph_bkexch_corr)
#
#     # generate low ph backexchange if saturation not achieved
#     low_ph_bkexch_new_arr = np.zeros(len(low_ph_bkexch_satur['backexchange_array']))
#     for ind, (low_ph_bkexch, high_ph_bkexch, low_ph_satur) in enumerate(zip(low_ph_bkexch_satur['backexchange_array'],
#                                                                             high_ph_bkexch_satur['backexchange_array'],
#                                                                             low_ph_bkexch_satur['saturation_bool_list'])):
#         if low_ph_satur:
#             low_ph_bkexch_new_arr[ind] = low_ph_bkexch
#         else:
#             new_low_ph_bkexch = (high_ph_bkexch - odr_output.beta[1])/odr_output.beta[0]
#             low_ph_bkexch_new_arr[ind] = new_low_ph_bkexch
#
#     plot_backexchange_correlation(low_ph_bkexch_corr,
#                                   high_ph_bkexch_corr,
#                                   low_ph_nocorr_backexch=low_ph_bkexch_nocorr,
#                                   high_ph_nocorr_backexch=high_ph_bkexch_nocorr,
#                                   corr_slope=odr_output.beta[0],
#                                   corr_intercept=odr_output.beta[1],
#                                   corr_r=corr_coef[0, 1],
#                                   output_path=plot_output_path)
#
#     write_low_high_backexchange_correlation(slope=odr_output.beta[0],
#                                             intercept=odr_output.beta[1],
#                                             rvalue=corr_coef[0, 1],
#                                             std_slope_error=odr_output.sd_beta[0],
#                                             std_intercept_error=odr_output.sd_beta[1],
#                                             output_path=bkexchange_corr_output_path)
#
#     write_low_high_backexchange_array(low_ph_protein_name=low_ph_protname_list,
#                                       high_ph_protein_name=high_ph_protname_list,
#                                       sequence_array=sequence_list,
#                                       low_ph_backexchange_array=low_ph_bkexch_satur['backexchange_array'],
#                                       high_ph_backexchange_array=high_ph_bkexch_satur['backexchange_array'],
#                                       low_ph_saturation=low_ph_bkexch_satur['saturation_bool_list'],
#                                       high_ph_saturation=high_ph_bkexch_satur['saturation_bool_list'],
#                                       corr_include_indices=corr_include_indices,
#                                       low_ph_backexchange_new=low_ph_bkexch_new_arr,
#                                       output_path=bkexchange_output_path)
#
#     if return_flag:
#         out_dict = dict()
#         out_dict['low_ph_backexchange_array'] = low_ph_bkexch_satur['backexchange_array']
#         out_dict['high_ph_backexchange_array'] = high_ph_bkexch_satur['backexchange_array']
#         out_dict['low_ph_saturation_bool_list'] = low_ph_bkexch_satur['saturation_bool_list']
#         out_dict['high_ph_saturation_bool_list'] = high_ph_bkexch_satur['saturation_bool_list']
#         return out_dict
#
#

def gen_parser_args():
    """
    generate parser
    :return: parser
    """

    parser = argparse.ArgumentParser(prog='GEN BACKEXCHANGE HIGH LOW', description='Generate backexchange for high and low phs')
    parser.add_argument('-s', '--samplepath', help='merge sample list file path .csv')
    parser.add_argument('-ldf', '--lowdfrac', help='low d2o fraction')
    parser.add_argument('-ldp', '--lowdpur', help='low d2o purity')
    parser.add_argument('-hdf', '--highdfrac', help='high d2o fraction')
    parser.add_argument('-hdp', '--highdpur', help='high d2o purity')
    parser.add_argument('-hbc', '--highbkexcorr', help='high ph backexchange correction fielpath')
    parser.add_argument('-lbc', '--lowbkexcorr', help='low ph backexchange correction filepath')
    parser.add_argument('-smr', '--saturmassrate', help='saturation mass rate threshold')
    parser.add_argument('-blb', '--bklowbound', help='backexchange low bound for correlation')
    parser.add_argument('-bub', '--bkupbound', help='backexchange upper bound for correlation')
    parser.add_argument('-bco', '--bkcorrout', help='backexchange correlation csv output path')
    parser.add_argument('-bcp', '--bkplotout', help='backexchange correlation plot path')
    parser.add_argument('-bko', '--bkoutput', help='backexchange high low csv output path')
    return parser


def gen_parser_args_v2():
    """
    generate parser
    :return: parser
    """

    parser = argparse.ArgumentParser(prog='GEN BACKEXCHANGE HIGH LOW', description='Generate backexchange for high and low phs')
    parser.add_argument('-lpfl', '--lowphfilelist', nargs='+', help='low ph file path list')
    parser.add_argument('-hpfl', '--highphfilelist', nargs='+', help='high ph file path list')
    parser.add_argument('-lpnl', '--lowphnamelist', nargs='+', help='low ph name list')
    parser.add_argument('-hpnl', '--highphnamelist', nargs='+', help='high ph name list')
    parser.add_argument('-sl', '--sequencelist', nargs='+', help='sequence list')
    parser.add_argument('-ldf', '--lowdfrac', help='low d2o fraction')
    parser.add_argument('-ldp', '--lowdpur', help='low d2o purity')
    parser.add_argument('-hdf', '--highdfrac', help='high d2o fraction')
    parser.add_argument('-hdp', '--highdpur', help='high d2o purity')
    parser.add_argument('-hbc', '--highbkexcorr', help='high ph backexchange correction fielpath')
    parser.add_argument('-lbc', '--lowbkexcorr', help='low ph backexchange correction filepath')
    parser.add_argument('-smr', '--saturmassrate', help='saturation mass rate threshold')
    parser.add_argument('-blb', '--bklowbound', help='backexchange low bound for correlation')
    parser.add_argument('-bub', '--bkupbound', help='backexchange upper bound for correlation')
    parser.add_argument('-bco', '--bkcorrout', help='backexchange correlation csv output path')
    parser.add_argument('-bcp', '--bkplotout', help='backexchange correlation plot path')
    parser.add_argument('-bko', '--bkoutput', help='backexchange high low csv output path')
    return parser


def run_from_parser():

    parser = gen_parser_args()

    options = parser.parse_args()

    gen_low_high_ph_bkexchange_(merge_sample_list_fpath=options.samplepath,
                                low_ph_d2o_frac=float(options.lowdfrac),
                                low_ph_d2o_pur=float(options.lowdpur),
                                high_ph_d2o_frac=float(options.highdfrac),
                                high_ph_d2o_pur=float(options.highdpur),
                                lowph_backexchange_correction_filepath=options.lowbkexcorr,
                                highph_backexchange_correction_filepath=options.highbkexcorr,
                                bkexchange_ub=float(options.bkupbound),
                                bkexchange_lb=float(options.bklowbound),
                                saturation_mass_rate_threshold=float(options.saturmassrate),
                                plot_output_path=options.bkplotout,
                                bkexchange_corr_output_path=options.bkcorrout,
                                bkexchange_output_path=options.bkoutput,
                                return_flag=False)


def run_from_parser_v2():

    parser = gen_parser_args_v2()

    options = parser.parse_args()

    gen_low_high_ph_bkexchange_from_file_list(list_of_low_ph_hxms_files=options.lowphfilelist,
                                              list_of_high_ph_hxms_files=options.highphfilelist,
                                              low_ph_protname_list=options.lowphnamelist,
                                              high_ph_protname_list=options.highphnamelist,
                                              sequence_list=options.sequencelist,
                                              low_ph_d2o_frac=float(options.lowdfrac),
                                              low_ph_d2o_pur=float(options.lowdpur),
                                              high_ph_d2o_frac=float(options.highdfrac),
                                              high_ph_d2o_pur=float(options.highdpur),
                                              lowph_backexchange_correction_filepath=options.lowbkexcorr,
                                              highph_backexchange_correction_filepath=options.highbkexcorr,
                                              bkexchange_ub=float(options.bkupbound),
                                              bkexchange_lb=float(options.bklowbound),
                                              saturation_mass_rate_threshold=float(options.saturmassrate),
                                              plot_output_path=options.bkplotout,
                                              bkexchange_corr_output_path=options.bkcorrout,
                                              bkexchange_output_path=options.bkoutput,
                                              return_flag=False)


if __name__ == '__main__':

    run_from_parser_v2()
