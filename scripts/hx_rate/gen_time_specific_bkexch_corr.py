import os.path
from collections import OrderedDict
from hxrate.hxdata import load_data_from_hdx_ms_dist, write_backexchange_correction_array
from hxrate.methods import gauss_fit_to_isotope_dist_array, normalize_mass_distribution_array
import numpy as np
import argparse
from dataclasses import dataclass
import matplotlib.pyplot as plt


@dataclass
class MassRateProtein(object):
    protein_name: str = None
    timepoints: np.ndarray = None
    mass_diff_rate_array: np.ndarray = None
    num_timepoints: int = None


@dataclass
class MassRateTimepoint(object):
    protein_name: str = None
    mass_rate_value: float = None


def calc_mass_diff_rate(ref_mass, mass_array):

    diff_rate = np.zeros(len(mass_array))
    for ind, mass_val in enumerate(mass_array):
        rate = (mass_val - ref_mass) / ref_mass
        diff_rate[ind] = rate
    return diff_rate


def gen_mass_rate_protein(protein_name: str,
                          isotope_distribution: np.ndarray,
                          timepoints: np.ndarray):

    # fit gaussian to the isotope distributions
    gauss_fit = gauss_fit_to_isotope_dist_array(isotope_dist=isotope_distribution)

    # generate a centroid array
    centroid_array = np.array([x.centroid for x in gauss_fit])

    mass_rate_arr_ = calc_mass_diff_rate(ref_mass=centroid_array[-1], mass_array=centroid_array)

    mass_rate_protein = MassRateProtein(protein_name=protein_name,
                                        timepoints=timepoints,
                                        mass_diff_rate_array=mass_rate_arr_,
                                        num_timepoints=len(timepoints))

    return mass_rate_protein


def gen_list_of_tp_diff_rates(list_of_mass_rate_proteins):

    tp_mass_rates_dict = OrderedDict()

    for mass_rate_prots in list_of_mass_rate_proteins:

        for ind, (tp, mass_rate_val) in enumerate(zip(mass_rate_prots.timepoints, mass_rate_prots.mass_diff_rate_array)):

            if tp not in tp_mass_rates_dict.keys():

                tp_mass_rates_dict[tp] = []
                tp_mass_rates_dict[tp].append(MassRateTimepoint(protein_name=mass_rate_prots.protein_name,
                                                                mass_rate_value=mass_rate_val))

            else:
                tp_mass_rates_dict[tp].append(MassRateTimepoint(protein_name=mass_rate_prots.protein_name,
                                                                mass_rate_value=mass_rate_val))

    return tp_mass_rates_dict


def gen_filtered_rates_dict(tp_mass_rates_dict, rate_threshold, ch_rate_threshold, min_number_points):

    filter_mass_rates_dict = OrderedDict()

    for keys, values in tp_mass_rates_dict.items():

        if keys != 0:

            if len(values) > 0:

                ch_rate_threshold_flag = False

                while not ch_rate_threshold_flag:

                    filter_mass_rates_dict[keys] = []

                    for mass_rate_tp in values:

                        abs_mass_rate_value = abs(mass_rate_tp.mass_rate_value)
                        if abs_mass_rate_value <= rate_threshold:
                            filter_mass_rates_dict[keys].append(mass_rate_tp)

                    if len(filter_mass_rates_dict[keys]) < min_number_points:
                        rate_threshold = rate_threshold * (1+ch_rate_threshold)
                    else:
                        ch_rate_threshold_flag = True

            else:
                
                filter_mass_rates_dict[keys] = []

    return filter_mass_rates_dict


def gen_correction_array(filtered_tp_mass_rates_dict):

    timepoint_list = []
    avg_mass_rate_list = []

    for key, values in filtered_tp_mass_rates_dict.items():

        mass_rate_value_arr = np.array([x.mass_rate_value for x in values])
        avg_mass_rate_value = np.median(mass_rate_value_arr)
        avg_mass_rate_list.append(avg_mass_rate_value)
        timepoint_list.append(key)

    timepoint_list = [0] + timepoint_list
    avg_mass_rate_list = [-0.9] + avg_mass_rate_list

    timepoints_arr = np.array(timepoint_list)
    avg_mass_rate_arr = np.array(avg_mass_rate_list)

    sort_index = np.argsort(timepoints_arr)

    return timepoints_arr[sort_index], avg_mass_rate_arr[sort_index]


def gen_list_of_mass_rate_proteins(list_of_hx_ms_files,
                                   hxms_dist_fpath_delim_str):

    list_of_mass_rate_proteins = []

    for ind, hx_ms_fpath in enumerate(list_of_hx_ms_files):

        hxms_dist_fname = os.path.split(hx_ms_fpath)[1]
        prot_name = hxms_dist_fname.split(hxms_dist_fpath_delim_str)[0]

        # tp, iso_dist_array = load_data_from_hdx_ms_dist_(hx_ms_fpath)
        data_dict = load_data_from_hdx_ms_dist(fpath=hx_ms_fpath)
        tp = data_dict['tp']
        iso_dist_array = data_dict['mass_dist']
        norm_dist = normalize_mass_distribution_array(mass_dist_array=iso_dist_array)
        gauss_fit = gauss_fit_to_isotope_dist_array(isotope_dist=norm_dist)
        centroid_array = np.array([x.centroid for x in gauss_fit])
        mass_rate_arr_ = calc_mass_diff_rate(ref_mass=centroid_array[-1], mass_array=centroid_array)
        mass_rate_prot = MassRateProtein(protein_name=prot_name,
                                         timepoints=tp,
                                         mass_diff_rate_array=mass_rate_arr_,
                                         num_timepoints=len(tp))
        list_of_mass_rate_proteins.append(mass_rate_prot)

    return list_of_mass_rate_proteins


def plot_mass_rate_all(list_of_mass_rate_prot_objs, correction_tp_arr, correction_bk_arr, output_path):
    # todo: add param description

    for mass_rate_prot in list_of_mass_rate_prot_objs:

        plt.plot(mass_rate_prot.timepoints, mass_rate_prot.mass_diff_rate_array, color='grey', alpha=0.1)

    plt.plot(correction_tp_arr, correction_bk_arr, color='red', ls='--')

    plt.xscale('log')
    plt.ylim((-1, 1))

    plt.grid(axis='x', alpha=0.25)
    plt.grid(axis='y', alpha=0.25)

    plt.xlabel('Timepoint')
    plt.ylabel('dM rate')

    plt.savefig(output_path)
    plt.close()


def gen_bkexch_correction_(list_of_hxms_files,
                           hxms_dist_fpath_delim_str,
                           rate_threshold,
                           ch_rate_threshold,
                           min_num_points,
                           bkexch_corr_csv_output_path=None,
                           bkexch_corr_plot_output_path=None,
                           return_flag=False):

    list_of_mass_rate_prot_objs = gen_list_of_mass_rate_proteins(list_of_hx_ms_files=list_of_hxms_files,
                                                                 hxms_dist_fpath_delim_str=hxms_dist_fpath_delim_str)

    list_of_tp_rate_dict = gen_list_of_tp_diff_rates(list_of_mass_rate_proteins=list_of_mass_rate_prot_objs)

    filter_tp_rate_dict = gen_filtered_rates_dict(tp_mass_rates_dict=list_of_tp_rate_dict,
                                                  rate_threshold=rate_threshold,
                                                  ch_rate_threshold=ch_rate_threshold,
                                                  min_number_points=min_num_points)

    timepoint_arr, bkexch_corr_arr = gen_correction_array(filtered_tp_mass_rates_dict=filter_tp_rate_dict)

    if bkexch_corr_csv_output_path is not None:
        write_backexchange_correction_array(timepoints=timepoint_arr,
                                            backexchange_correction_array=bkexch_corr_arr,
                                            output_path=bkexch_corr_csv_output_path)

    if bkexch_corr_plot_output_path is not None:
        plot_mass_rate_all(list_of_mass_rate_prot_objs=list_of_mass_rate_prot_objs,
                           correction_tp_arr=timepoint_arr,
                           correction_bk_arr=bkexch_corr_arr,
                           output_path=bkexch_corr_plot_output_path)

    if return_flag:
        correction_dict = dict()
        for ind, (tp_, corr_) in enumerate(zip(timepoint_arr, bkexch_corr_arr)):
            correction_dict[tp_] = corr_

        return correction_dict


def gen_parser_arguments():
    """
    generate commandline arguements to generate backexchange correction file
    :return:parser
    """
    parser = argparse.ArgumentParser(prog='Backexchange correction',
                                     description='Generate backexchange correction path')
    parser.add_argument('-l', '--list_of_hxms_files', nargs='+', help='list of hx ms dist files')
    parser.add_argument('-s', '--delim_string', help='delimiter string')
    parser.add_argument('-r', '--rate_tol', help='rate tolerance', default=0.08)
    parser.add_argument('-c', '--change_rate_threshold', help='change rate threshold', default=0.02)
    parser.add_argument('-m', '--min_number_points', help='min number of points', default=3)
    parser.add_argument('-p', '--plot_path', help='plot file path', default='../../workfolder/bkexch_corr.pdf')
    parser.add_argument('-o', '--output_path', help='correction file path', default='../../workfolder/bkexch_corr.csv')
    parser.add_argument('-g', '--return_flag', help='return the value bool', default=False)

    return parser


def run_from_parser():

    parser_ = gen_parser_arguments()
    options = parser_.parse_args()

    gen_bkexch_correction_(list_of_hxms_files=options.list_of_hxms_files,
                           hxms_dist_fpath_delim_str=options.delim_string,
                           rate_threshold=float(options.rate_tol),
                           ch_rate_threshold=float(options.change_rate_threshold),
                           min_num_points=int(options.min_number_points),
                           bkexch_corr_csv_output_path=options.output_path,
                           bkexch_corr_plot_output_path=options.plot_path)


if __name__ == '__main__':

    run_from_parser()
