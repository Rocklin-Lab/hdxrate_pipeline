# hx_rate_fitting code
# authors: Gabe Rocklin and Suggie

import argparse
import pandas as pd
import numpy as np
from hxrate.ratefit import fit_rate
from hxrate.methods import normalize_mass_distribution_array
from hxrate.hxdata import load_data_from_hdx_ms_dist, load_tp_dependent_dict


def fit_rate_from_to_file(prot_name: str,
                          sequence: str,
                          hx_ms_dist_fpath: str or list,
                          d2o_fraction: float,
                          d2o_purity: float,
                          prot_rt_name: str = 'PROTEIN_RT',
                          merge_exp: bool = False,
                          exp_label: str or list = None,
                          usr_backexchange: float or list = None,
                          low_high_backexchange_list_fpath: str = None,
                          low_high_bkexch_corr: bool = True,
                          backexchange_corr_fpath: str or list = None,
                          backexchange_array_fpath: str or list = None,
                          adjust_backexchange: bool = True,
                          num_chains: int = 4,
                          num_warmups: int = 100,
                          num_samples: int = 500,
                          sample_backexchange: bool = False,
                          save_posterior_samples: bool = True,
                          hx_rate_output_path: str = None,
                          hx_rate_csv_output_path: str = None,
                          hx_isotope_dist_output_path: str = None,
                          merge_hx_ms_dist_output_path: str = None,
                          hx_rate_plot_path: str = None,
                          posterior_plot_path: str = None,
                          return_flag: bool = False) -> object:
    # todo: add param descriptions for the function here

    # make it accept list

    if merge_exp:
        timepoints, mass_dist_list, timepoint_label = [], [], []
        for hxmsdist_fp in hx_ms_dist_fpath:
            data_dict = load_data_from_hdx_ms_dist(fpath=hxmsdist_fp)
            timepoints.append(data_dict['tp'])
            mass_dist_list.append(data_dict['mass_dist'])
            if 'tp_ind' in list(data_dict.keys()):
                timepoint_label.append(data_dict['tp_ind'])
            else:
                timepoint_label = None

        norm_dist = [normalize_mass_distribution_array(mass_dist_array=x) for x in mass_dist_list]

        bkexch_corr_dict = None
        if backexchange_corr_fpath is not None:
            bkexch_corr_dict = []
            for bkexch_corr_fpath_ in backexchange_corr_fpath:
                bkexch_corr_dict.append(load_tp_dependent_dict(bkexch_corr_fpath_))

        backexchange_array = None
        if backexchange_array_fpath is not None:
            backexchange_array = []
            for ind, bkexch_arr_fpath_ in enumerate(backexchange_array_fpath):
                bkexch_arr_dict = load_tp_dependent_dict(filepath=bkexch_arr_fpath_)
                backexchange_array_tp = np.zeros(len(timepoints[ind]))
                for ind_, tp_ in enumerate(timepoints[ind]):
                    backexchange_array_tp[ind_] = bkexch_arr_dict[tp_]
                backexchange_array.append(backexchange_array_tp)

        # get backexchange value from user list or backexchange list file
        backexchange_value = None

        # if backexchange list filepath is present
        if low_high_backexchange_list_fpath is not None:
            backexchange_value = []
            low_high_bkexch_df = pd.read_csv(low_high_backexchange_list_fpath)
            # todo: delete the following low_high_name df assignment after making sure that is already present in
            #  the csv file
            low_high_bkexch_df['low_high_name'] = np.array([x+'_'+y for x, y in zip(low_high_bkexch_df['low_ph_protein_name'].values, low_high_bkexch_df['high_ph_protein_name'].values)])
            low_high_bkexch_prot_df = low_high_bkexch_df[low_high_bkexch_df['low_high_name'] == prot_rt_name]
            if low_high_bkexch_corr:
                low_user_backexchange = low_high_bkexch_prot_df['low_ph_backexchange_all_corr'].values
            else:
                low_user_backexchange = low_high_bkexch_prot_df['low_ph_backexchange'].values
            high_user_backexchange = low_high_bkexch_prot_df['high_ph_backexchange'].values

            # low_high_name should be unique so we expect a single value on the backexchange array
            # put the backexchange values in a list
            if len(low_user_backexchange) == 1:
                backexchange_value.append(low_user_backexchange[0])
            else:
                backexchange_value.append(None)

            if len(high_user_backexchange) == 1:
                backexchange_value.append(high_user_backexchange[0])
            else:
                backexchange_value.append(None)

        # if backexchange list filepath is not present, see if user backexchange has values and if not get a
        # list of Nones
        else:
            if usr_backexchange is None:
                backexchange_value = [None for _ in range(len(hx_ms_dist_fpath))]
            else:
                backexchange_value = usr_backexchange
    else:

        if type(hx_ms_dist_fpath) == list:
            hx_ms_dist_fpath = hx_ms_dist_fpath[0]

        # timepoints, mass_dist = load_data_from_hdx_ms_dist_(fpath=hx_ms_dist_fpath)
        data_dict = load_data_from_hdx_ms_dist(fpath=hx_ms_dist_fpath)
        timepoints = data_dict['tp']
        mass_dist = data_dict['mass_dist']
        if 'tp_ind' in list(data_dict.keys()):
            timepoint_label = data_dict['tp_ind']
        else:
            timepoint_label = None

        norm_dist = normalize_mass_distribution_array(mass_dist_array=mass_dist)

        # get backexchange correction
        bkexch_corr_dict = None
        if backexchange_corr_fpath is not None:
            if type(backexchange_corr_fpath) == list:
                backexchange_corr_fpath = backexchange_corr_fpath[0]
            bkexch_corr_dict = load_tp_dependent_dict(filepath=backexchange_corr_fpath)

        # get backexchange array
        backexchange_array = None
        if backexchange_array_fpath is not None:
            if type(backexchange_array_fpath) == list:
                backexchange_array_fpath = backexchange_array_fpath[0]
            backexchange_array_dict = load_tp_dependent_dict(filepath=backexchange_array_fpath)
            backexchange_array = np.zeros(len(timepoints))
            for ind, tp in enumerate(timepoints):
                backexchange_array[ind] = backexchange_array_dict[tp]

        if type(usr_backexchange) == list:
            usr_backexchange = float(usr_backexchange[0])
        elif type(usr_backexchange) == float:
            usr_backexchange = usr_backexchange
        else:
            usr_backexchange = None
        backexchange_value = usr_backexchange

    # make sure exp label gets passed as a single string if its only one exp
    if type(exp_label) == list:
        if len(exp_label) > 1:
            exp_label = exp_label
        else:
            exp_label = exp_label[0]

    # fit rate
    hxrate_out = fit_rate(prot_name=prot_name,
                          sequence=sequence,
                          time_points=timepoints,
                          timepoint_label=timepoint_label,
                          exp_label=exp_label,
                          norm_mass_distribution_array=norm_dist,
                          d2o_fraction=d2o_fraction,
                          d2o_purity=d2o_purity,
                          num_chains=num_chains,
                          num_warmups=num_warmups,
                          num_samples=num_samples,
                          prot_rt_name=prot_rt_name,
                          merge_exp=merge_exp,
                          sample_backexchange=sample_backexchange,
                          adj_backexchange=adjust_backexchange,
                          backexchange_value=backexchange_value,
                          backexchange_correction_dict=bkexch_corr_dict,
                          backexchange_array=backexchange_array)

    # write hxrate_pipeline as csv file
    if hx_rate_csv_output_path is not None:
        hxrate_out.write_rates_to_csv(hx_rate_csv_output_path)

    # write hxrate_pipeline pred distributions as a csv file
    if hx_isotope_dist_output_path is not None:
        hxrate_out.pred_dist_to_csv(hx_isotope_dist_output_path)

    # write merge distribution
    if merge_exp:
        if merge_hx_ms_dist_output_path is not None:
            hxrate_out.exp_merge_dist_to_csv(merge_hx_ms_dist_output_path)

    # save the ratefit output dictionary to pickle
    if hx_rate_output_path is not None:
        hxrate_out.output_to_pickle(hx_rate_output_path, save_posterior_samples=save_posterior_samples)

    # plot distribution
    if hx_rate_plot_path is not None:
        hxrate_out.plot_hxrate_output(hx_rate_plot_path)

    # plot posterior distribution
    if posterior_plot_path is not None:
        hxrate_out.plot_bayes_samples(posterior_plot_path)

    if return_flag:
        return hxrate_out.output


def gen_parser_arguments():
    """
    generate commandline arguements to run the hx rate fitting algorithm
    :return:parser
    """

    parser = argparse.ArgumentParser(prog='HX_RATE_FIT', description='Run HX rate fitting algorithm')
    parser.add_argument('-p', '--protname', help='protein name', default='PROTEIN')
    parser.add_argument('-pr', '--protrtname', help='protein rt name', default='PROTEIN_RT')
    parser.add_argument('-s', '--sequence', help='protein sequence one letter amino acid', default='PROTEIN')
    parser.add_argument('-i', '--hxdist', nargs='+', help='hx mass distribution input file .csv')
    parser.add_argument('-df', '--d2o_frac', help='d2o fracion', type=float, default=0.95)
    parser.add_argument('-dp', '--d2o_pur', help='d2o purity', type=float, default=0.95)
    parser.add_argument('-ub', '--user_bkexchange', help='user defined backexchange', nargs='+', type=float or list, default=None)
    parser.add_argument('-expl', '--exp_label', help='Exp label in str or [str, str,..]', nargs='+', type=str or list, default=None)
    parser.add_argument('-bcl', '--backexchange_low_high_list_fpath', type=str, help='backexchange low high ph list filepath .csv', default=None)
    parser.add_argument('-bcf', '--bkexchange_corr_fpath', nargs='+', help='backexchange correction filepath .csv')
    parser.add_argument('-baf', '--bkexchange_array_fpath', nargs='+', help='backexchange array filepath .csv')
    parser.add_argument('--adjust_backexchange', help='adjust backexchange boolean', default=True,
                        action=argparse.BooleanOptionalAction)
    parser.add_argument('--lohibkexchcorr', help='correct backexchange for low ph data using high ph data', default=True,
                        action=argparse.BooleanOptionalAction)
    parser.add_argument('-nc', '--num_chains', help='number of independent markov chains for MCMC', type=int, default=4)
    parser.add_argument('-nw', '--num_warmups', help='number of warmups for MCMC', type=int, default=100)
    parser.add_argument('-ns', '--num_samples', help='number of samples for MCMC', type=int, default=500)
    parser.add_argument('--merge', help='merge distributions', default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument('--sample_backexchange', help='sample backexchange for MCMC', default=False,
                        action=argparse.BooleanOptionalAction)
    parser.add_argument('--save_posterior_samples', help='save posterior samples', default=True,
                        action=argparse.BooleanOptionalAction)
    parser.add_argument('-o', '--output_pickle_file', help='output pickle filepath', default=None)
    parser.add_argument('-or', '--output_rate_csv', help='output rates csv filepath', default=None)
    parser.add_argument('-op', '--output_rate_plot', help='output_rate_plot filepath', default=None)
    parser.add_argument('-opp', '--output_posterior_plot', help='output_posterior_plot filepath', default=None)
    parser.add_argument('-od', '--output_iso_dist', help='output isotope distribution filepath', default=None)
    parser.add_argument('-md', '--merge_dist_output', help='merge distribution output .csv filepath', default=None)

    return parser


def hx_rate_fitting_from_parser(parser):
    """
    from the parser arguments, generate essential arguments for hx rate fitting function and run the function
    :param parser: parser
    :return:
    """

    options = parser.parse_args()

    fit_rate_from_to_file(prot_name=options.protname,
                          prot_rt_name=options.protrtname,
                          sequence=options.sequence,
                          hx_ms_dist_fpath=options.hxdist,
                          d2o_fraction=options.d2o_frac,
                          d2o_purity=options.d2o_pur,
                          usr_backexchange=options.user_bkexchange,
                          merge_exp=options.merge,
                          exp_label=options.exp_label,
                          backexchange_corr_fpath=options.bkexchange_corr_fpath,
                          backexchange_array_fpath=options.bkexchange_array_fpath,
                          low_high_backexchange_list_fpath=options.backexchange_low_high_list_fpath,
                          adjust_backexchange=options.adjust_backexchange,
                          num_chains=options.num_chains,
                          num_warmups=options.num_warmups,
                          num_samples=options.num_samples,
                          sample_backexchange=options.sample_backexchange,
                          hx_rate_output_path=options.output_pickle_file,
                          hx_rate_csv_output_path=options.output_rate_csv,
                          hx_isotope_dist_output_path=options.output_iso_dist,
                          hx_rate_plot_path=options.output_rate_plot,
                          merge_hx_ms_dist_output_path=options.merge_dist_output,
                          posterior_plot_path=options.output_posterior_plot)


if __name__ == '__main__':
    pass

    parser_ = gen_parser_arguments()
    hx_rate_fitting_from_parser(parser_)
