import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
from scipy.optimize import minimize
from scipy import interpolate
import matplotlib.pyplot as plt
import methods
import hx_rate_fit
import hxdata


def interpol_and_calc_mse(log_high_tp_refac,
                          log_low_tp,
                          high_centroids,
                          low_centroids,
                          slice_index_low,
                          slice_index_high):

    # generate a new time axis that's common for both high and low ph data

    new_axis = log_high_tp_refac[slice_index_low:slice_index_high]
    high_centroids_comp = high_centroids[slice_index_low:slice_index_high]

    # create an interpolation function for low data
    low_interp_func = interpolate.interp1d(x=log_low_tp, y=low_centroids)

    # generate interpolated centroids on the new axis
    low_centroids_comp = low_interp_func(x=new_axis)

    try:
        mse = mean_squared_error(y_true=high_centroids_comp, y_pred=low_centroids_comp)
    except ValueError:
        mse = 100

    return mse


def gen_mse_with_factor_shift(factor: float,
                              low_tp: np.ndarray,
                              low_centroids: np.ndarray,
                              high_tp: np.ndarray,
                              high_centroids: np.ndarray) -> float:
    """
    gen comparison arrays when the high tp is shifted by some factor
    :param factor: factor by which to shift the high timepoints
    :param low_tp: low ph timepoints not including zero
    :param low_centroids: low ph centroinds not including for zero
    :param high_tp: high ph timepoints not including zero
    :param high_centroids: high ph centroids not including for zero
    :return: mse
    """

    # shift the tp by some factor
    # factor = 5
    high_tp_refac = high_tp * factor

    log_low_tp = np.log10(low_tp)
    log_high_tp_refac = np.log10(high_tp_refac)

    # generate low and high indices where to slice the timepoint for comparison that is common on both high and low data
    tp_diff_high = np.subtract(log_low_tp[-1], log_high_tp_refac)
    tp_diff_low = np.subtract(log_low_tp[0], log_high_tp_refac)
    # tp_diff = np.subtract(log_low_tp[-1], log_high_tp_refac)
    high_ind = 0
    low_ind = 0
    for ind, diff in enumerate(tp_diff_high):
        if diff < 0:
            high_ind = ind
            break

    for ind2, diff in enumerate(tp_diff_low):
        if diff < 0:
            low_ind = ind2
            break

    # use the high and low indices to slice the data
    # make statements for if high_ind == len(tpdiff) -1 or low_ind == len(tpdifflow) -1

    mse = 100
    if high_ind == 0 and low_ind == 0:
        mse = 100
    elif high_ind <= low_ind:
        mse = 100
    elif low_ind >= high_ind:
        mse = 100
    elif high_ind > 0:
        if high_ind == len(tp_diff_high) -1:
            mse = interpol_and_calc_mse(log_high_tp_refac=log_high_tp_refac,
                                        log_low_tp=log_low_tp,
                                        high_centroids=high_centroids,
                                        low_centroids=low_centroids,
                                        slice_index_low=low_ind,
                                        slice_index_high=-1)
        else:
            mse = interpol_and_calc_mse(log_high_tp_refac=log_high_tp_refac,
                                        log_low_tp=log_low_tp,
                                        high_centroids=high_centroids,
                                        low_centroids=low_centroids,
                                        slice_index_low=low_ind,
                                        slice_index_high=high_ind)


    # elif high_ind == 0 and low_ind > 0:
    #     if low_ind >= high_ind:
    #         mse = 100
    #     else:
    #         mse = interpol_and_calc_mse(log_high_tp_refac=log_high_tp_refac,
    #                                     log_low_tp=log_low_tp,
    #                                     high_centroids=high_centroids,
    #                                     low_centroids=low_centroids,
    #                                     slice_index_low=low_ind,
    #                                     slice_index_high=high_ind)
    # elif high_ind > 0 and low_ind == 0:
    #     if high_ind == len(tp_diff_high) -1:
    #         mse = interpol_and_calc_mse(log_high_tp_refac=log_high_tp_refac,
    #                                     log_low_tp=log_low_tp,
    #                                     high_centroids=high_centroids,
    #                                     low_centroids=low_centroids,
    #                                     slice_index_low=low_ind,
    #                                     slice_index_high=None)
    # elif high_ind > 0 and low_ind > 0:
    #     mse = 100
    #
    #
    #
    #
    # # make sure high end is not zero. That means it didn't find any common timepoints for comparison. in that case
    # # return the full centroids for comparison
    # if high_ind == 0:
    #     print('no common timepoints for comparison. Returning the entire array as it is for computing mse.')
    #     high_centroids_comp = high_centroids
    #     low_centroids_comp = low_centroids
    # # if high_ind == len(log_high_tp_refac)-1:
    # #     new_axis = log_high_tp_refac
    # #     high_centroids_comp
    # elif high_ind == len(tp_diff) - 1:
    #     new_axis = log_high_tp_refac
    #     high_centroids_comp = high_centroids
    #     low_interp_func = interpolate.interp1d(x=log_low_tp, y=low_centroids)
    #     low_centroids_comp = low_interp_func(x=new_axis)
    # else:
    #     # generate a new time axis that's common for both high and low ph data
    #     new_axis = log_high_tp_refac[:high_ind]
    #     high_centroids_comp = high_centroids[:high_ind]
    #
    #     # create an interpolation function for low data
    #     low_interp_func = interpolate.interp1d(x=log_low_tp, y=low_centroids)
    #
    #     # generate interpolated centroids on the new axis
    #     low_centroids_comp = low_interp_func(x=new_axis)
    #
    # # computing mse
    # try:
    #     mse = mean_squared_error(y_true=high_centroids_comp, y_pred=low_centroids_comp)
    # except ValueError:
    #     mse = 100

    return mse


def optimize_factor_for_alignment(low_tp: np.ndarray,
                                  low_centroids: np.ndarray,
                                  high_tp: np.ndarray,
                                  high_centroids: np.ndarray,
                                  opt_init_guess: float = 1.0) -> dict:
    """

    :param low_tp: non zero low ph timepoints
    :param low_centroids: non zero low ph centroids
    :param high_tp: non zero high ph timepoints
    :param high_centroids: non zero high ph centroids
    :param opt_method: optimize method
    :param opt_init_guess: initial guess for the factor
    :return: optimize dictionary
    """

    # optimize for factor that best aligns the centroids from low to high ph data

    opt = minimize(fun=gen_mse_with_factor_shift,
                   x0=np.array([opt_init_guess]),
                   args=(low_tp, low_centroids, high_tp, high_centroids))

    # store important info in a dictionary for output
    out_dict = dict()
    # out_dict['opt_method'] = opt_method
    out_dict['opt_init_guess'] = opt_init_guess
    out_dict['opt_mse'] = opt.fun
    out_dict['opt_nfev'] = opt.nfev
    out_dict['opt_nit'] = opt.nit
    out_dict['opt_success'] = opt.success
    out_dict['opt_message'] = opt.message
    out_dict['opt_x'] = opt.x[0]

    return out_dict


def gen_diff_array(array):

    diff_arr = np.zeros(len(array)-1)
    for num in range(len(diff_arr)):
        diff_arr[num] = array[num+1] - array[num]

    return diff_arr


def check_diff_arr(arr):

    unique = True
    num = 0
    for ind, item in enumerate(arr):
        if item <= 1e-1:
            unique = False
            num = ind
            break

    return unique, num


def modify_to_add_uniq_values(array, add_value=0.4):

    diff_arr = gen_diff_array(array)

    uniq, num = check_diff_arr(diff_arr)

    while not uniq:

        array[num+1] = array[num+1] + add_value
        diff_arr = gen_diff_array(array)
        uniq, num = check_diff_arr(diff_arr)

    return array


def gen_merged_dstirbution(sequence: str,
                           low_tp: np.ndarray,
                           low_dist: np.ndarray,
                           low_d2o_frac: float,
                           low_d2o_purity: float,
                           low_user_backexchange: float,
                           low_backexchange_correction_dict: dict,
                           high_tp: np.ndarray,
                           high_dist: np.ndarray,
                           high_d2o_frac: float,
                           high_d2o_purity: float,
                           high_user_backexchange: float,
                           high_backexchange_correction_dict: dict,
                           factor_init_list: list = [1, 10, 100, 1000, 10000]) -> dict:
    """

    :param sequence:
    :param low_tp:
    :param low_dist:
    :param low_d2o_frac:
    :param low_d2o_purity:
    :param low_user_backexchange:
    :param low_backexchange_correction_arr:
    :param high_tp:
    :param high_dist:
    :param high_d2o_frac:
    :param high_d2o_purity:
    :param high_user_backexchange:
    :param high_backexchange_correction_arr:
    :return:
    """

    # fit gaussians to both distributions
    low_gauss_fit_list = methods.gauss_fit_to_isotope_dist_array(low_dist)
    high_gauss_fit_list = methods.gauss_fit_to_isotope_dist_array(high_dist)

    # get backexchange array

    low_backexchange_obj = hx_rate_fit.calc_back_exchange(sequence=sequence,
                                                          experimental_isotope_dist=low_dist[-1],
                                                          timepoints_array=low_tp,
                                                          d2o_fraction=low_d2o_frac,
                                                          d2o_purity=low_d2o_purity,
                                                          usr_backexchange=low_user_backexchange,
                                                          backexchange_corr_dict=low_backexchange_correction_dict)

    high_backexchange_obj = hx_rate_fit.calc_back_exchange(sequence=sequence,
                                                           experimental_isotope_dist=high_dist[-1],
                                                           timepoints_array=high_tp,
                                                           d2o_fraction=high_d2o_frac,
                                                           d2o_purity=high_d2o_purity,
                                                           usr_backexchange=high_user_backexchange,
                                                           backexchange_corr_dict=high_backexchange_correction_dict)

    # generate centroids array for both distributions
    low_centroids = np.array([x.centroid for x in low_gauss_fit_list])
    high_centroids = np.array([x.centroid for x in high_gauss_fit_list])

    # generate centroids with backexhcnage corrections
    low_centroids_corr = methods.correct_centroids_using_backexchange(centroids=low_centroids,
                                                                      backexchange_array=low_backexchange_obj.backexchange_array,
                                                                      include_zero_dist=True)
    high_centroids_corr = methods.correct_centroids_using_backexchange(centroids=high_centroids,
                                                                       backexchange_array=high_backexchange_obj.backexchange_array,
                                                                       include_zero_dist=True)
    opt_list = []
    for factor_init_value in factor_init_list:
        opt_dict = optimize_factor_for_alignment(low_tp=low_tp[1:],
                                                 low_centroids=low_centroids_corr[1:],
                                                 high_tp=high_tp[1:],
                                                 high_centroids=high_centroids_corr[1:],
                                                 opt_init_guess=factor_init_value)
        opt_list.append(opt_dict)

    mse_arr = np.array([x['opt_mse'] for x in opt_list])
    min_mse_ind = np.argmin(mse_arr)

    opt_ = opt_list[min_mse_ind]

    # if the factor value is less than 0
    if opt_['opt_x'] < 0:
        opt_['opt_x'] = 1.1
        opt_['opt_success'] = False
        opt_['opt_message'] = 'Negative factor value. Setting the value to 1.1'

    # generate the high ph timepoints with factor applied
    high_tp_factor = high_tp * opt_['opt_x']

    # merge the low and high ph timepoints together
    merge_tp = np.concatenate([low_tp, high_tp_factor[1:]])  # don't inlcude the zero timepoint from high ph data

    # generate a sorting key
    tp_sort_index = np.argsort(merge_tp)

    # generate a sorted timepoint array
    sorted_merged_tp = merge_tp[tp_sort_index]

    # if timepoints are less that 0.1 apart, make it at least 0.1 apart
    mod_sorted_merge_tp = modify_to_add_uniq_values(array=sorted_merged_tp, add_value=0.4)
    # generate index for sorting
    tp_sort_index_2 = np.argsort(mod_sorted_merge_tp)
    sorted_merged_tp = mod_sorted_merge_tp[tp_sort_index_2]

    # ind = 0
    # while ind < len(sorted_merged_tp)-1:
    #     t_diff = sorted_merged_tp[ind+1] - sorted_merged_tp[ind]
    #     if t_diff <= 0.1:
    #         sorted_merged_tp[ind+1] = sorted_merged_tp[ind+1] + 0.1
    #     # if t_diff == 0:
    #     #     sorted_merged_tp[ind+1] = sorted_merged_tp[ind+1] + 0.5
    #     ind += 1

    # generate a sorted merged distributions
    merge_dist = np.concatenate([low_dist, high_dist[1:, :]], axis=0)
    sorted_merged_dist = merge_dist[tp_sort_index, :]
    sorted_merged_dist = sorted_merged_dist[tp_sort_index_2, :]

    # generate a sorted merged backexchange array
    merged_backexchange = np.concatenate([low_backexchange_obj.backexchange_array, high_backexchange_obj.backexchange_array[1:]])
    sorted_merged_backexchange = merged_backexchange[tp_sort_index]
    sorted_merged_backexchange = sorted_merged_backexchange[tp_sort_index_2]

    merged_backexchange_correction_array = methods.gen_backexchange_correction_from_backexchange_array(backexchange_array=sorted_merged_backexchange)

    out_dict = dict()
    out_dict['low_tp'] = low_tp
    out_dict['low_centroids'] = low_centroids
    out_dict['low_centroids_corr'] = low_centroids_corr
    out_dict['high_tp'] = high_tp
    out_dict['high_tp_factor'] = high_tp_factor
    out_dict['high_centroids'] = high_centroids
    out_dict['high_centroids_corr'] = high_centroids_corr
    out_dict['factor_optimization'] = opt_
    out_dict['merged_timepoints'] = sorted_merged_tp
    out_dict['merged_dist'] = sorted_merged_dist
    out_dict['merged_backexchange'] = sorted_merged_backexchange
    out_dict['merged_backexchange_correction_array'] = merged_backexchange_correction_array

    return out_dict


def get_backexchange_corr_arr(bakcexchange_corr_fpath):
    """
    given backexchange correction filepath
    :param bakcexchange_corr_fpath: backexchange correction filepath
    :return: backexchange correction
    """
    df = pd.read_csv(bakcexchange_corr_fpath)
    backexchange_corr_arr = df.iloc[:, 1].values
    return backexchange_corr_arr


def plot_merge_centroids(low_tp, low_centroids, high_tp_refactored, high_centroids, factor, mse, output_path):
    """

    :param low_tp:
    :param low_centroids:
    :param high_tp_refactored:
    :param high_centroids:
    :param output_path:
    :return:
    """
    plt.scatter(np.log(low_tp), low_centroids, marker='o', color='blue', edgecolor='black', label='low_ph', s=60)
    plt.scatter(np.log(high_tp_refactored), high_centroids, marker='o', color='orange', edgecolor='black',
                label='high_ph', s=60)
    plt.xlabel('log time (seconds)')
    plt.ylabel('Center of mass')
    plt.grid()
    plt.legend(loc='best')
    plt.title('Merge low and high ph data (Factor: %.4f | MSE: %.4f)' % (factor, mse))
    plt.savefig(output_path)
    plt.close()


# def gen_high_low_merged_from_to_file(sequence: str,
#                                      low_ph_data_fpath: str,
#                                      low_d2o_frac: float,
#                                      low_d2o_purity: float,
#                                      low_user_backexchange: float or None,
#                                      low_backexchange_corr_fpath: str,
#                                      high_ph_data_fpath: str,
#                                      high_d2o_frac: float,
#                                      high_d2o_purity: float,
#                                      high_user_backexchange: float or None,
#                                      high_backexchange_corr_fpath: str,
#                                      merged_backexchange_fpath: str or None,
#                                      merged_data_fpath: str or None,
#                                      merged_backexchange_correction_fpath: str or None,
#                                      factor_fpath: str or None,
#                                      merge_plot_fpath: str or None,
#                                      return_flag: bool = False):
#
#     # read the ph data
#     low_tp, low_dists = hxdata.load_data_from_hdx_ms_dist_(low_ph_data_fpath)
#     high_tp, high_dists = hxdata.load_data_from_hdx_ms_dist_(high_ph_data_fpath)
#
#     # get the backexchange correction array
#     if low_backexchange_corr_fpath is not None:
#         low_backexchange_corr_arr = get_backexchange_corr_arr(low_backexchange_corr_fpath)
#     else:
#         low_backexchange_corr_arr = None
#
#     if high_backexchange_corr_fpath is not None:
#         high_backexchange_corr_arr = get_backexchange_corr_arr(high_backexchange_corr_fpath)
#     else:
#         high_backexchange_corr_arr = None
#
#     merged_data_dict = gen_merged_dstirbution(sequence=sequence,
#                                               low_tp=low_tp,
#                                               low_dist=low_dists,
#                                               low_d2o_frac=low_d2o_frac,
#                                               low_d2o_purity=low_d2o_purity,
#                                               low_user_backexchange=low_user_backexchange,
#                                               low_backexchange_correction_arr=low_backexchange_corr_arr,
#                                               high_tp=high_tp,
#                                               high_dist=high_dists,
#                                               high_d2o_frac=high_d2o_frac,
#                                               high_d2o_purity=high_d2o_purity,
#                                               high_user_backexchange=high_user_backexchange,
#                                               high_backexchange_correction_arr=high_backexchange_corr_arr)
#
#     # write the merged timepoints and isotope distribution
#     if merged_data_fpath is not None:
#         hxdata.write_isotope_dist_timepoints(timepoints=merged_data_dict['merged_timepoints'],
#                                              isotope_dist_array=merged_data_dict['merged_dist'],
#                                              output_path=merged_data_fpath)
#
#     # write the merged backexchange
#     if merged_backexchange_fpath is not None:
#         hxdata.write_backexchange_array(timepoints=merged_data_dict['merged_timepoints'],
#                                         backexchange_array=merged_data_dict['merged_backexchange'],
#                                         output_path=merged_backexchange_fpath)
#
#     # write the merged backexchange correction
#     if merged_backexchange_correction_fpath is not None:
#         hxdata.write_backexchange_correction_array(timepoints=merged_data_dict['merged_timepoints'],
#                                                    backexchange_correction_array=merged_data_dict['merged_backexchange_correction_array'],
#                                                    output_path=merged_backexchange_correction_fpath)
#
#     # write factor and optimization info to a file
#     if factor_fpath is not None:
#         hxdata.write_merge_factor(merge_factor=merged_data_dict['factor_optimization']['opt_x'],
#                                   opt_mse=merged_data_dict['factor_optimization']['opt_mse'],
#                                   opt_nfev=merged_data_dict['factor_optimization']['opt_nfev'],
#                                   opt_nit=merged_data_dict['factor_optimization']['opt_nit'],
#                                   opt_success=merged_data_dict['factor_optimization']['opt_success'],
#                                   opt_message=merged_data_dict['factor_optimization']['opt_message'],
#                                   output_path=factor_fpath)
#
#     if merge_plot_fpath is not None:
#         plot_merge_centroids(low_tp=merged_data_dict['low_tp'][1:],
#                              low_centroids=merged_data_dict['low_centroids_corr'][1:],
#                              high_tp_refactored=merged_data_dict['high_tp_factor'][1:],
#                              high_centroids=merged_data_dict['high_centroids_corr'][1:],
#                              factor=merged_data_dict['factor_optimization']['opt_x'],
#                              mse=merged_data_dict['factor_optimization']['opt_mse'],
#                              output_path=merge_plot_fpath)
#
#     if return_flag:
#         return merged_data_dict


def gen_high_low_merged_from_to_file_v2(sequence: str,
                                        low_ph_data_fpath: str,
                                        low_ph_prot_name: str,
                                        low_d2o_frac: float,
                                        low_d2o_purity: float,
                                        low_user_backexchange: float or None,
                                        low_backexchange_corr_fpath: str,
                                        high_ph_data_fpath: str,
                                        high_ph_prot_name: str,
                                        high_d2o_frac: float,
                                        high_d2o_purity: float,
                                        high_user_backexchange: float or None,
                                        high_backexchange_corr_fpath: str,
                                        low_high_backexchange_list_fpath: str or None,
                                        merged_backexchange_fpath: str or None,
                                        merged_data_fpath: str or None,
                                        merged_backexchange_correction_fpath: str or None,
                                        factor_fpath: str or None,
                                        merge_plot_fpath: str or None,
                                        return_flag: bool = False):

    # read the ph data
    low_tp, low_dists = hxdata.load_data_from_hdx_ms_dist_(low_ph_data_fpath)
    high_tp, high_dists = hxdata.load_data_from_hdx_ms_dist_(high_ph_data_fpath)

    # get the backexchange correction array
    if low_backexchange_corr_fpath is not None:
        low_backexchange_corr_dict = hxdata.load_tp_dependent_dict(low_backexchange_corr_fpath)
    else:
        low_backexchange_corr_dict = None

    if high_backexchange_corr_fpath is not None:
        high_backexchange_corr_dict = hxdata.load_tp_dependent_dict(high_backexchange_corr_fpath)
    else:
        high_backexchange_corr_dict = None

    if low_high_backexchange_list_fpath is not None:
        low_high_bkexch_df = pd.read_csv(low_high_backexchange_list_fpath)
        low_high_name = np.array([x+'_'+y for x, y in zip(low_high_bkexch_df['low_ph_protein_name'].values, low_high_bkexch_df['high_ph_protein_name'].values)])
        low_high_bkexch_df['low_high_name'] = low_high_name
        low_high_bkexch_prot_df = low_high_bkexch_df[low_high_bkexch_df['low_high_name'] == low_ph_prot_name + '_' + high_ph_prot_name]
        low_user_backexchange = low_high_bkexch_prot_df['low_ph_backexchange_new'].values[0]
        high_user_backexchange = low_high_bkexch_prot_df['high_ph_backexchange'].values[0]
        print('low_user_backexchange: ', low_user_backexchange)
        print('high_user_backexchange: ', high_user_backexchange)

    merged_data_dict = gen_merged_dstirbution(sequence=sequence,
                                              low_tp=low_tp,
                                              low_dist=low_dists,
                                              low_d2o_frac=low_d2o_frac,
                                              low_d2o_purity=low_d2o_purity,
                                              low_user_backexchange=low_user_backexchange,
                                              low_backexchange_correction_dict=low_backexchange_corr_dict,
                                              high_tp=high_tp,
                                              high_dist=high_dists,
                                              high_d2o_frac=high_d2o_frac,
                                              high_d2o_purity=high_d2o_purity,
                                              high_user_backexchange=high_user_backexchange,
                                              high_backexchange_correction_dict=high_backexchange_corr_dict)

    # write the merged timepoints and isotope distribution
    if merged_data_fpath is not None:
        hxdata.write_isotope_dist_timepoints(timepoints=merged_data_dict['merged_timepoints'],
                                             isotope_dist_array=merged_data_dict['merged_dist'],
                                             output_path=merged_data_fpath)

    # write the merged backexchange
    if merged_backexchange_fpath is not None:
        hxdata.write_backexchange_array(timepoints=merged_data_dict['merged_timepoints'],
                                        backexchange_array=merged_data_dict['merged_backexchange'],
                                        output_path=merged_backexchange_fpath)

    # write the merged backexchange correction
    if merged_backexchange_correction_fpath is not None:
        hxdata.write_backexchange_correction_array(timepoints=merged_data_dict['merged_timepoints'],
                                                   backexchange_correction_array=merged_data_dict['merged_backexchange_correction_array'],
                                                   output_path=merged_backexchange_correction_fpath)

    # write factor and optimization info to a file
    if factor_fpath is not None:
        hxdata.write_merge_factor(merge_factor=merged_data_dict['factor_optimization']['opt_x'],
                                  opt_mse=merged_data_dict['factor_optimization']['opt_mse'],
                                  opt_nfev=merged_data_dict['factor_optimization']['opt_nfev'],
                                  opt_nit=merged_data_dict['factor_optimization']['opt_nit'],
                                  opt_success=merged_data_dict['factor_optimization']['opt_success'],
                                  opt_message=merged_data_dict['factor_optimization']['opt_message'],
                                  output_path=factor_fpath)

    if merge_plot_fpath is not None:
        plot_merge_centroids(low_tp=merged_data_dict['low_tp'][1:],
                             low_centroids=merged_data_dict['low_centroids_corr'][1:],
                             high_tp_refactored=merged_data_dict['high_tp_factor'][1:],
                             high_centroids=merged_data_dict['high_centroids_corr'][1:],
                             factor=merged_data_dict['factor_optimization']['opt_x'],
                             mse=merged_data_dict['factor_optimization']['opt_mse'],
                             output_path=merge_plot_fpath)

    if return_flag:
        return merged_data_dict


# def gen_parse_args():
#
#     parser = argparse.ArgumentParser()
#
#     parser.add_argument('-seq', '--sequence', help='protein sequence one letter code')
#     parser.add_argument('-ldata', '--lowphdata', help='low ph data filepath')
#     parser.add_argument('-ldf', '--lowdfrac', help='low d2o fraction')
#     parser.add_argument('-ldp', '--lowdpur', help='low d2o purity')
#     parser.add_argument('-lbk', '--lowbkex', help='low ph backexchange user value')
#     parser.add_argument('-lbkc', '--lowbkexcorr', help='low ph backexchange correction filepath')
#     parser.add_argument('-hdata', '--highphdata', help='high ph data filepath')
#     parser.add_argument('-hdf', '--highdfrac', help='high d2o fraction')
#     parser.add_argument('-hdp', '--highdpur', help='high d2o purity')
#     parser.add_argument('-hbk', '--highbkex', help='high ph backexchange user value')
#     parser.add_argument('-hbkc', '--highbkexcorr', help='high ph backexchange correction filepath')
#     parser.add_argument('-mbk', '--mergebkex', help='merge backexchange output filepath .csv')
#     parser.add_argument('-mbkc', '--mergebkexcorr', help='merge backexchange correction output filepath .csv')
#     parser.add_argument('-mdp', '--mergedatapath', help='merge distribution output filepath .csv')
#     parser.add_argument('-mpp', '--mergeplotpath', help='merge plot output filepath .pdf')
#     parser.add_argument('-mfp', '--mergefactorpath', help='merge factor output filepath .csv')
#
#     return parser


def gen_parse_args_v2():

    parser = argparse.ArgumentParser()

    parser.add_argument('-seq', '--sequence', help='protein sequence one letter code')
    parser.add_argument('-ldata', '--lowphdata', help='low ph data filepath')
    parser.add_argument('-lpn', '--lowprotname', help='low ph protein name')
    parser.add_argument('-hpn', '--highprotname', help='high ph protein name')
    parser.add_argument('-ldf', '--lowdfrac', help='low d2o fraction')
    parser.add_argument('-ldp', '--lowdpur', help='low d2o purity')
    parser.add_argument('-lbk', '--lowbkex', help='low ph backexchange user value')
    parser.add_argument('-lbkc', '--lowbkexcorr', help='low ph backexchange correction filepath')
    parser.add_argument('-hdata', '--highphdata', help='high ph data filepath')
    parser.add_argument('-hdf', '--highdfrac', help='high d2o fraction')
    parser.add_argument('-hdp', '--highdpur', help='high d2o purity')
    parser.add_argument('-hbk', '--highbkex', help='high ph backexchange user value')
    parser.add_argument('-hbkc', '--highbkexcorr', help='high ph backexchange correction filepath')
    parser.add_argument('-lhblf', '--lowhighbkexchlistfile', help='low high backexchange list filepath')
    parser.add_argument('-mbk', '--mergebkex', help='merge backexchange output filepath .csv')
    parser.add_argument('-mbkc', '--mergebkexcorr', help='merge backexchange correction output filepath .csv')
    parser.add_argument('-mdp', '--mergedatapath', help='merge distribution output filepath .csv')
    parser.add_argument('-mpp', '--mergeplotpath', help='merge plot output filepath .pdf')
    parser.add_argument('-mfp', '--mergefactorpath', help='merge factor output filepath .csv')

    return parser


# def run_from_parser():
#
#     parser_ = gen_parse_args()
#     options = parser_.parse_args()
#
#     gen_high_low_merged_from_to_file(sequence=options.sequence,
#                                      low_ph_data_fpath=options.lowphdata,
#                                      low_d2o_frac=float(options.lowdfrac),
#                                      low_d2o_purity=float(options.lowdpur),
#                                      low_user_backexchange=float(options.lowbkex),
#                                      low_backexchange_corr_fpath=options.lowbkexcorr,
#                                      high_ph_data_fpath=options.highphdata,
#                                      high_d2o_frac=float(options.highdfrac),
#                                      high_d2o_purity=float(options.highdpur),
#                                      high_user_backexchange=float(options.highbkex),
#                                      high_backexchange_corr_fpath=options.highbkexcorr,
#                                      merged_backexchange_fpath=options.mergebkex,
#                                      merged_backexchange_correction_fpath=options.mergebkexcorr,
#                                      merged_data_fpath=options.mergedatapath,
#                                      factor_fpath=options.mergefactorpath,
#                                      merge_plot_fpath=options.mergeplotpath,
#                                      return_flag=False)


def run_from_parser_v2():

    parser_ = gen_parse_args_v2()
    options = parser_.parse_args()

    gen_high_low_merged_from_to_file_v2(sequence=options.sequence,
                                        low_ph_data_fpath=options.lowphdata,
                                        low_ph_prot_name=options.lowprotname,
                                        low_d2o_frac=float(options.lowdfrac),
                                        low_d2o_purity=float(options.lowdpur),
                                        low_user_backexchange=None,
                                        low_backexchange_corr_fpath=options.lowbkexcorr,
                                        high_ph_data_fpath=options.highphdata,
                                        high_ph_prot_name=options.highprotname,
                                        high_d2o_frac=float(options.highdfrac),
                                        high_d2o_purity=float(options.highdpur),
                                        high_user_backexchange=None,
                                        high_backexchange_corr_fpath=options.highbkexcorr,
                                        low_high_backexchange_list_fpath=options.lowhighbkexchlistfile,
                                        merged_backexchange_fpath=options.mergebkex,
                                        merged_backexchange_correction_fpath=options.mergebkexcorr,
                                        merged_data_fpath=options.mergedatapath,
                                        factor_fpath=options.mergefactorpath,
                                        merge_plot_fpath=options.mergeplotpath,
                                        return_flag=False)


if __name__ == '__main__':

    run_from_parser_v2()

    # sequence = 'HMYTVKEGETLAQIAKNYKTSLEEIIKANKIEDPSLLSIGQRLFIP'
    # # #
    # low_ph_data_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/hx_rates_library/test_repo/merge_test/ph6/A0A1J4VH79.1_187-230_17.85463_winner_multibody.cpickle.zlib.csv'
    # high_ph_data_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/hx_rates_library/test_repo/merge_test/ph9/A0A1J4VH79.1_187-230_17.85472_winner_multibody.cpickle.zlib.csv'
    # # #
    # d2o_frac = 0.95
    # d2o_pur = 0.95
    # # #
    # low_bkexch_corr_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/hx_rates_library/test_repo/merge_test/backexchange/low_ph_bkexch_corr.csv'
    # high_bkexch_corr_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/hx_rates_library/test_repo/merge_test/backexchange/high_ph_bkexch_corr.csv'
    # # #
    # high_low_bkexch_list = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/test/Lib11/backexchange/high_low_backexchange_list.csv'
    # # #
    # output_path = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/hx_rates_library/test_repo/merge_test/test_output'
    # # #
    # gen_high_low_merged_from_to_file_v2(sequence=sequence,
    #                                     low_ph_data_fpath=low_ph_data_fpath,
    #                                     low_ph_prot_name='A0A1M4T304.1_578-633_14.88156',
    #                                     low_d2o_frac=d2o_frac,
    #                                     low_d2o_purity=d2o_pur,
    #                                     low_user_backexchange=None,
    #                                     low_backexchange_corr_fpath=low_bkexch_corr_fpath,
    #                                     high_ph_data_fpath=high_ph_data_fpath,
    #                                     high_ph_prot_name='A0A1M4T304.1_578-633_14.67927',
    #                                     high_d2o_frac=d2o_frac,
    #                                     high_d2o_purity=d2o_pur,
    #                                     high_user_backexchange=None,
    #                                     high_backexchange_corr_fpath=high_bkexch_corr_fpath,
    #                                     low_high_backexchange_list_fpath=high_low_bkexch_list,
    #                                     merged_backexchange_fpath=output_path + '/merge_bkexch_.csv',
    #                                     merged_data_fpath=output_path + '/merge_data_.csv',
    #                                     merged_backexchange_correction_fpath=output_path + '/merge_bkexch_corr_.csv',
    #                                     factor_fpath=output_path + '/merge_bkexch_factor_.csv',
    #                                     merge_plot_fpath=output_path + '/merge_data_plot.pdf',
    #                                     return_flag=False)

    # common_backexchange_filepath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bkexch_corr_output/common_2.csv_bkexchange.csv'
    # common_backexchange_df = pd.read_csv(common_backexchange_filepath)
    # # common_backexchange_df['comp_bool'] = np.array([1 for x in range(len(common_backexchange_df['sequence'].values))])
    #
    # low_backexchange_corr_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/lib15_ph6_sample.csv_backexchange_correction.csv'
    # high_backexchange_corr_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/lib15_ph7_sample.csv_backexchange_correction.csv'
    #
    # comp_df = common_backexchange_df[common_backexchange_df['comp_bool'] == 1]
    #
    # low_fpath_dir = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/lib15_ph6'
    # high_fpath_dir = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/lib15_ph7'
    #
    # d2o_frac = 0.95
    # d2o_pur = 0.95
    #
    # output_top_dir = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/merged_data_ph6_ph7'
    #
    # import glob
    # import os
    #
    # for ind, (low_name, high_name, seq) in enumerate(zip(comp_df['ph6'].values, comp_df['ph7'].values, comp_df['sequence'].values)):
    #
    #     low_phdata_fpath = glob.glob(low_fpath_dir + '/' + low_name + '*.csv')[0]
    #     high_phdata_fpath = glob.glob(high_fpath_dir + '/' + high_name + '*.csv')[0]
    #
    #     output_dir = hxdata.make_new_dir(output_top_dir + '/' + low_name + '_' + high_name)
    #     merged_data_path = os.path.join(output_dir, low_name + '_' + high_name + '_merged_data.csv')
    #     merged_backexchange_path = os.path.join(output_dir, low_name + '_' + high_name + '_merged_backexchange.csv')
    #     merged_backexchange_corr_path = os.path.join(output_dir, low_name + '_' + high_name + '_merged_backexchange_correction.csv')
    #     factor_path = os.path.join(output_dir, low_name + '_' + high_name + '_factor.csv')
    #     factor_plot_path = os.path.join(output_dir, low_name + '_' + high_name + '_factor_plot.pdf')
    #
    #     gen_high_low_merged_from_to_file(sequence=seq,
    #                                      low_ph_data_fpath=low_phdata_fpath,
    #                                      low_d2o_frac=d2o_frac,
    #                                      low_d2o_purity=d2o_pur,
    #                                      low_user_backexchange=None,
    #                                      low_backexchange_corr_fpath=low_backexchange_corr_fpath,
    #                                      high_ph_data_fpath=high_phdata_fpath,
    #                                      high_d2o_frac=d2o_frac,
    #                                      high_d2o_purity=d2o_pur,
    #                                      high_user_backexchange=None,
    #                                      high_backexchange_corr_fpath=high_backexchange_corr_fpath,
    #                                      merged_backexchange_fpath=merged_backexchange_path,
    #                                      merged_backexchange_correction_fpath=merged_backexchange_corr_path,
    #                                      merged_data_fpath=merged_data_path,
    #                                      factor_fpath=factor_path,
    #                                      merge_plot_fpath=factor_plot_path)
