# hx_rate_fitting code
# authors: Gabe Rocklin and Suggie

import multiprocessing as mp
import numpy as np
from dataclasses import dataclass

import pandas as pd
from scipy.special import expit
from scipy.optimize import fmin_powell
from methods import isotope_dist_from_PoiBin, gen_temp_rates, gen_theoretical_isotope_dist_for_all_timepoints, \
    normalize_mass_distribution_array, hx_rate_fitting_optimization, compute_rmse_exp_thr_iso_dist, \
    gauss_fit_to_isotope_dist_array, convert_hxrate_object_to_dict, plot_hx_rate_fitting_, gen_corr_backexchange
from hxdata import load_data_from_hdx_ms_dist_, write_pickle_object, write_hx_rate_output, write_isotope_dist_timepoints, load_backexhange_correction_
import time


@dataclass
class ExpData(object):
    """
    class container to store exp data
    """
    protein_name: str = None
    protein_sequence: str = None
    timepoints: np.ndarray = None
    exp_isotope_dist_array: np.ndarray = None
    gauss_fit: list = None


@dataclass
class HXRate(object):
    """
    class container to store hx rate fitting data
    """
    exp_data: object = None
    back_exchange: object = None
    back_exchange_res_subtract: int = None
    optimization_cost: float = None
    optimization_func_evals: int = None
    optimization_init_rate_guess: np.ndarray = None
    hx_rates: np.ndarray = None
    thr_isotope_dist_array: np.ndarray = None
    thr_isotope_dist_gauss_fit: list = None
    fit_rmse_each_timepoint: np.ndarray = None
    total_fit_rmse: float = None


@dataclass
class BackExchange(object):
    """
    class container to store backexchange data
    """
    backexchange_value: float = None
    fit_rmse: float = None
    theoretical_isotope_dist: np.ndarray = None
    backexchange_correction_dict: dict = None
    backexchange_array: np.ndarray = None

    def gen_tp_backexchange_array(self, timepoints):
        """
        generate timepoint specific backexchange array
        :param timepoints:
        :return:
        """
        if self.backexchange_correction_dict is None:
            self.backexchange_array = np.array([self.backexchange_value for _ in range(len(timepoints))])
        else:
            tp_backexchange_corr_array = np.zeros(len(timepoints))
            for ind, tp in enumerate(timepoints):
                tp_backexchange_corr_array[ind] = self.backexchange_correction_dict[tp]
            self.backexchange_array = gen_corr_backexchange(mass_rate_array=tp_backexchange_corr_array,
                                                            fix_backexchange_value=self.backexchange_value)


def calc_back_exchange(sequence: str,
                       experimental_isotope_dist: np.ndarray,
                       d2o_fraction: float,
                       d2o_purity: float,
                       timepoint: float = 1e9,
                       usr_backexchange: float = None) -> object:
    """
    calculate back exchange from the experimental isotope distribution
    :param sequence: protein sequence
    :param experimental_isotope_dist: experimental isotope distribution to be used for backexchange calculation
    :param d2o_fraction: d2o fraction
    :param d2o_purity: d2o purity
    :param timepoint: hdx timepoint in seconds. default to 1e9 for full deuteration
    :param usr_backexchange: if usr backexchange provided, will generate backexchange object based on that backexchange
    value
    :return: backexchange data object
    """
    # set high rates for calculating back exchange
    rates = gen_temp_rates(sequence=sequence,
                           rate_value=1e4)

    # set the number of bins for isotope distribtuion
    num_bins = len(experimental_isotope_dist)

    if usr_backexchange is None:
        print('\nCALCULATING BACK EXCHANGE ... ')

        opt = fmin_powell(lambda x: compute_rmse_exp_thr_iso_dist(exp_isotope_dist=experimental_isotope_dist,
                                                                  thr_isotope_dist=isotope_dist_from_PoiBin(sequence=sequence,
                                                                                timepoint=timepoint,
                                                                                inv_backexchange=expit(x),
                                                                                rates=rates,
                                                                                d2o_fraction=d2o_fraction,
                                                                                d2o_purity=d2o_purity,
                                                                                num_bins=num_bins),
                                                                  squared=False),
                          x0=0.05,
                          disp=True)

        back_exchange = 1 - expit(opt)[0]

    else:
        print('\nSETTING USER BACK EXCHANGE ... ')
        back_exchange = usr_backexchange

    thr_iso_dist_full_deut = isotope_dist_from_PoiBin(sequence=sequence,
                                                      timepoint=timepoint,
                                                      inv_backexchange=1 - back_exchange,
                                                      rates=rates,
                                                      d2o_fraction=d2o_fraction,
                                                      d2o_purity=d2o_purity,
                                                      num_bins=num_bins)

    fit_rmse = compute_rmse_exp_thr_iso_dist(exp_isotope_dist=experimental_isotope_dist,
                                             thr_isotope_dist=thr_iso_dist_full_deut,
                                             squared=False)

    backexchange_obj = BackExchange(backexchange_value=back_exchange,
                                    fit_rmse=fit_rmse,
                                    theoretical_isotope_dist=thr_iso_dist_full_deut)

    return backexchange_obj


def fit_rate_without_backexchange_adjust(prot_name: str,
                                         sequence: str,
                                         time_points: np.ndarray,
                                         norm_mass_distribution_array: np.ndarray,
                                         d2o_fraction: float,
                                         d2o_purity: float,
                                         opt_iter: int = 500,
                                         opt_temp: float = 0.00003,
                                         opt_step_size: float = 0.02,
                                         multi_proc: bool = True,
                                         number_of_cores: int = 6,
                                         rate_init_list: list = None,
                                         free_energy_values: np.ndarray = None,
                                         temperature: float = None,
                                         backexchange_value: float = None,
                                         backexchange_correction_dict: dict = None) -> object:

    # todo: add param descritpions

    # initialize hxrate data object
    hxrate = HXRate()

    # store exp data in the object
    expdata = ExpData(protein_name=prot_name,
                      protein_sequence=sequence,
                      timepoints=time_points,
                      exp_isotope_dist_array=norm_mass_distribution_array)

    # fit gaussian to exp data
    expdata.gauss_fit = gauss_fit_to_isotope_dist_array(norm_mass_distribution_array)

    # store expdata object in hxrate data object
    hxrate.exp_data = expdata

    # calculate back exchange first
    if backexchange_value is None:
        back_exchange = calc_back_exchange(sequence=sequence,
                                           experimental_isotope_dist=norm_mass_distribution_array[-1],
                                           d2o_fraction=d2o_fraction,
                                           d2o_purity=d2o_purity)
    else:
        back_exchange = calc_back_exchange(sequence=sequence,
                                           experimental_isotope_dist=norm_mass_distribution_array[-1],
                                           d2o_fraction=d2o_fraction,
                                           d2o_purity=d2o_purity,
                                           usr_backexchange=backexchange_value)

    # enter backexchange correction dictionary to class
    back_exchange.backexchange_correction_dict = backexchange_correction_dict

    # generate timepoint specific backexchange array
    back_exchange.gen_tp_backexchange_array(timepoints=time_points)

    # if backexchange_array is None:
    #     if backexchange_correction_array is None:
    #         # generate an array of backexchange with same backexchange value to an array of length of timepoints
    #         back_exchange.backexchange_array = np.array([back_exchange.backexchange_value for x in time_points])
    #     else:
    #         back_exchange.backexchange_array = gen_corr_backexchange(mass_rate_array=backexchange_correction_array,
    #                                                                  fix_backexchange_value=back_exchange.backexchange_value)
    # else:
    #     back_exchange.backexchange_array = backexchange_array

    # if backexchange_correction_array is None:
    # # generate an array of backexchange with same backexchange value to an array of length of timepoints
    #     back_exchange.backexchange_array = np.array([back_exchange.backexchange_value for x in time_points])
    # else:
    #     back_exchange.backexchange_array = gen_corr_backexchange(mass_rate_array=backexchange_correction_array,
    #                                                              fix_backexchange_value=back_exchange.backexchange_value)

    # store backexchange object in hxrate object
    hxrate.back_exchange = back_exchange

    inv_back_exchange_array = np.subtract(1, back_exchange.backexchange_array)

    # generate a temporary rates to determine what residues can't exchange
    temp_rates = gen_temp_rates(sequence=sequence, rate_value=1)

    # get the indices of residues that don't exchange
    zero_indices = np.where(temp_rates == 0)[0]

    # calculate the number of residues that can exchange
    num_rates = len(temp_rates) - len(zero_indices)

    num_bins_ = len(norm_mass_distribution_array[0])

    print('\nHX RATE FITTING ... ')

    # get init rate list as initial guesses for rate fit optimization
    if rate_init_list is not None:
        init_rates_list = rate_init_list
    else:
        init_rates_list = [np.array([-2 for x in range(num_rates)]),
                           np.array([-5 for x in range(num_rates)]),
                           np.linspace(-7, 0, num_rates),
                           np.linspace(-7, -2, num_rates),
                           np.linspace(-8, -4, num_rates),
                           np.linspace(1, -12, num_rates)]

    # create a list to store the opt object
    store_opt_object = [0]

    if multi_proc:

        # set the max number of cores to be equal to the number of initial rate guesses
        if number_of_cores > len(init_rates_list):
            number_of_cores = len(init_rates_list)

        pool = mp.Pool(processes=number_of_cores)

        results = []

        for init_rate in init_rates_list:

            kwds_pool = {'init_rate_guess': init_rate,
                         'exp_isotope_dist_array': norm_mass_distribution_array,
                         'sequence': sequence,
                         'timepoints': time_points,
                         'inv_backexchange_array': inv_back_exchange_array,
                         'd2o_fraction': d2o_fraction,
                         'd2o_purity': d2o_purity,
                         'num_bins': num_bins_,
                         'free_energy_values': free_energy_values,
                         'temperature': temperature,
                         'opt_iter': opt_iter,
                         'opt_temp': opt_temp,
                         'opt_step_size': opt_step_size,
                         'return_tuple': True}

            pool_result = pool.apply_async(hx_rate_fitting_optimization, kwds=kwds_pool)
            results.append(pool_result)

        pool.close()

        tuple_results_list = [p.get() for p in results]
        opt_object_list = [tup[0] for tup in tuple_results_list]
        opt_obj_fun_list = np.array([x.fun for x in opt_object_list])
        min_fun_ind = np.argmin(opt_obj_fun_list)
        store_opt_object[0] = opt_object_list[min_fun_ind]
        init_rate_used = tuple_results_list[min_fun_ind][1]

    else:

        opt_cost = 10
        init_rate_final_ind = -1

        for ind, init_rate in enumerate(init_rates_list[:1]):

            opt_ = hx_rate_fitting_optimization(exp_isotope_dist_array=norm_mass_distribution_array,
                                                sequence=sequence,
                                                timepoints=time_points,
                                                inv_backexchange_array=inv_back_exchange_array,
                                                d2o_fraction=d2o_fraction,
                                                d2o_purity=d2o_purity,
                                                num_bins=num_bins_,
                                                free_energy_values=free_energy_values,
                                                temperature=temperature,
                                                init_rate_guess=init_rate,
                                                opt_iter=opt_iter,
                                                opt_temp=opt_temp,
                                                opt_step_size=opt_step_size)

            new_opt_cost = opt_.fun
            if new_opt_cost < opt_cost:
                init_rate_final_ind = ind
                store_opt_object[0] = opt_
            opt_cost = new_opt_cost

        init_rate_used = init_rates_list[init_rate_final_ind]

    opt_object = store_opt_object[0]
    hxrate.optimization_cost = opt_object.fun
    hxrate.optimization_func_evals = opt_object.nfev
    hxrate.optimization_init_rate_guess = init_rate_used
    hxrate.hx_rates = opt_object.x

    # generate theoretical isotope distribution array
    hxrate.thr_isotope_dist_array = gen_theoretical_isotope_dist_for_all_timepoints(sequence=sequence,
                                                                                    timepoints=time_points,
                                                                                    rates=np.exp(hxrate.hx_rates),
                                                                                    inv_backexchange_array=inv_back_exchange_array,
                                                                                    d2o_fraction=d2o_fraction,
                                                                                    d2o_purity=d2o_purity,
                                                                                    num_bins=num_bins_,
                                                                                    free_energy_values=free_energy_values,
                                                                                    temperature=temperature)

    # compute the fit mse for each distribution
    hxrate.fit_rmse_each_timepoint = np.zeros(len(time_points))
    for ind, (exp_dist, thr_dist) in enumerate(zip(norm_mass_distribution_array, hxrate.thr_isotope_dist_array)):
        hxrate.fit_rmse_each_timepoint[ind] = compute_rmse_exp_thr_iso_dist(exp_isotope_dist=exp_dist,
                                                                            thr_isotope_dist=thr_dist,
                                                                            squared=False)

    # compute the total fit mse distribution
    exp_isotope_dist_concat = np.concatenate(norm_mass_distribution_array)
    thr_isotope_dist_concat = np.concatenate(hxrate.thr_isotope_dist_array)
    thr_isotope_dist_concat[np.isnan(thr_isotope_dist_concat)] = 0
    hxrate.total_fit_rmse = compute_rmse_exp_thr_iso_dist(exp_isotope_dist=exp_isotope_dist_concat,
                                                          thr_isotope_dist=thr_isotope_dist_concat,
                                                          squared=False)

    # fit gaussian to the thr isotope dist
    hxrate.thr_isotope_dist_gauss_fit = gauss_fit_to_isotope_dist_array(hxrate.thr_isotope_dist_array)

    return hxrate


def fit_rate_with_backexchange_adjust(prot_name: str,
                                      sequence: str,
                                      time_points: np.ndarray,
                                      norm_mass_distribution_array: np.ndarray,
                                      d2o_fraction: float,
                                      d2o_purity: float,
                                      opt_iter: int = 500,
                                      opt_temp: float = 0.00003,
                                      opt_step_size: float = 0.02,
                                      multi_proc: bool = True,
                                      number_of_cores: int = 6,
                                      rate_init_list: list = None,
                                      free_energy_values: np.ndarray = None,
                                      temperature: float = None,
                                      backexchange_value: float = None,
                                      backexchange_correction_dict: dict = None,
                                      max_res_subtract_for_backexchange: int = 3,
                                      slow_rates_max_diff: float = 1.6) -> object:

    # todo: add param descritpions

    # initialize hxrate data object
    hxrate = HXRate()

    # store exp data in the object
    expdata = ExpData(protein_name=prot_name,
                      protein_sequence=sequence,
                      timepoints=time_points,
                      exp_isotope_dist_array=norm_mass_distribution_array)

    # fit gaussian to exp data
    expdata.gauss_fit = gauss_fit_to_isotope_dist_array(norm_mass_distribution_array)

    # store expdata object in hxrate data object
    hxrate.exp_data = expdata

    # set a flag for backexchange adjusting
    backexchange_adjust = False
    hxrate.back_exchange_res_subtract = 0

    # calculate backexchange first
    if backexchange_value is None:
        back_exchange = calc_back_exchange(sequence=sequence,
                                           experimental_isotope_dist=norm_mass_distribution_array[-1],
                                           d2o_fraction=d2o_fraction,
                                           d2o_purity=d2o_purity)
    else:
        back_exchange = calc_back_exchange(sequence=sequence,
                                           experimental_isotope_dist=norm_mass_distribution_array[-1],
                                           d2o_fraction=d2o_fraction,
                                           d2o_purity=d2o_purity,
                                           usr_backexchange=backexchange_value)

    # enter correction dict into the object
    back_exchange.backexchange_correction_dict = backexchange_correction_dict

    # generate timepoint specific backexhange array
    back_exchange.gen_tp_backexchange_array(timepoints=time_points)

    # set the original backexchange for reference
    original_backexchange = back_exchange.backexchange_value

    while backexchange_adjust is False:

        # if backexchange_correction_array is None:
        #     # generate an array of backexchange with same backexchange value to an array of length of timepoints
        #     back_exchange.backexchange_array = np.array([back_exchange.backexchange_value for x in time_points])
        # else:
        #     back_exchange.backexchange_array = gen_corr_backexchange(mass_rate_array=backexchange_correction_array,
        #                                                              fix_backexchange_value=back_exchange.backexchange_value)

        # store backexchange object in hxrate object
        hxrate.back_exchange = back_exchange

        inv_back_exchange_array = np.subtract(1, back_exchange.backexchange_array)

        # generate a temporary rates to determine what residues can't exchange
        temp_rates = gen_temp_rates(sequence=sequence, rate_value=1)

        # get the indices of residues that don't exchange
        zero_indices = np.where(temp_rates == 0)[0]

        # calculate the number of residues that can exchange
        num_rates = len(temp_rates) - len(zero_indices)

        num_bins_ = len(norm_mass_distribution_array[0])

        print('\nHX RATE FITTING ... ')

        # get init rate list as initial guesses for rate fit optimization
        if rate_init_list is not None:
            init_rates_list = rate_init_list
        else:
            init_rates_list = [np.array([-2 for x in range(num_rates)]),
                               np.array([-5 for x in range(num_rates)]),
                               np.linspace(-7, 0, num_rates),
                               np.linspace(-7, -2, num_rates),
                               np.linspace(-8, -4, num_rates),
                               np.linspace(1, -12, num_rates)]

        # create a list to store the opt object
        store_opt_object = [0]

        if multi_proc:

            # set the max number of cores to be equal to the number of initial rate guesses
            if number_of_cores > len(init_rates_list):
                number_of_cores = len(init_rates_list)

            pool = mp.Pool(processes=number_of_cores)

            results = []

            for init_rate in init_rates_list:

                kwds_pool = {'init_rate_guess': init_rate,
                             'exp_isotope_dist_array': norm_mass_distribution_array,
                             'sequence': sequence,
                             'timepoints': time_points,
                             'inv_backexchange_array': inv_back_exchange_array,
                             'd2o_fraction': d2o_fraction,
                             'd2o_purity': d2o_purity,
                             'num_bins': num_bins_,
                             'free_energy_values': free_energy_values,
                             'temperature': temperature,
                             'opt_iter': opt_iter,
                             'opt_temp': opt_temp,
                             'opt_step_size': opt_step_size,
                             'return_tuple': True}

                pool_result = pool.apply_async(hx_rate_fitting_optimization, kwds=kwds_pool)
                results.append(pool_result)

            pool.close()

            tuple_results_list = [p.get() for p in results]
            opt_object_list = [tup[0] for tup in tuple_results_list]
            opt_obj_fun_list = np.array([x.fun for x in opt_object_list])
            min_fun_ind = np.argmin(opt_obj_fun_list)
            store_opt_object[0] = opt_object_list[min_fun_ind]
            init_rate_used = tuple_results_list[min_fun_ind][1]

        else:

            opt_cost = 10
            init_rate_final_ind = -1

            for ind, init_rate in enumerate(init_rates_list[:1]):

                opt_ = hx_rate_fitting_optimization(exp_isotope_dist_array=norm_mass_distribution_array,
                                                    sequence=sequence,
                                                    timepoints=time_points,
                                                    inv_backexchange_array=inv_back_exchange_array,
                                                    d2o_fraction=d2o_fraction,
                                                    d2o_purity=d2o_purity,
                                                    num_bins=num_bins_,
                                                    free_energy_values=free_energy_values,
                                                    temperature=temperature,
                                                    init_rate_guess=init_rate,
                                                    opt_iter=opt_iter,
                                                    opt_temp=opt_temp,
                                                    opt_step_size=opt_step_size)

                new_opt_cost = opt_.fun
                if new_opt_cost < opt_cost:
                    init_rate_final_ind = ind
                    store_opt_object[0] = opt_
                opt_cost = new_opt_cost

            init_rate_used = init_rates_list[init_rate_final_ind]

        opt_object = store_opt_object[0]
        hxrate.optimization_cost = opt_object.fun
        hxrate.optimization_func_evals = opt_object.nfev
        hxrate.optimization_init_rate_guess = init_rate_used
        hxrate.hx_rates = opt_object.x

        # sort the rates
        sort_rates = sorted(hxrate.hx_rates)

        # get the difference of the slowest rate and 3rd slowest.
        rate_diff = sort_rates[2] - sort_rates[0]

        # if the rate difference is smaller than 1.6, hxrate optimization ends
        if rate_diff < slow_rates_max_diff:
            backexchange_adjust = True

        # else, backexchange is adjusted by reducing its value in proportion to having less residues.
        # then redo the hx rate fitting
        else:
            hxrate.back_exchange_res_subtract += 1
            # if the backexchange res subtract exceeds the max res subtract, terminate the fitting routine
            if hxrate.back_exchange_res_subtract > max_res_subtract_for_backexchange:
                backexchange_adjust = True
            else:
                print('adjusting backexchange value ... \n')
                backexchange_value = original_backexchange * ((len(sort_rates) + hxrate.back_exchange_res_subtract)/len(sort_rates))
                back_exchange = calc_back_exchange(sequence=sequence,
                                                   experimental_isotope_dist=norm_mass_distribution_array[-1],
                                                   d2o_fraction=d2o_fraction,
                                                   d2o_purity=d2o_purity,
                                                   usr_backexchange=backexchange_value)

                # enter correction dict into the object
                back_exchange.backexchange_correction_dict = backexchange_correction_dict

                # generate timepoint specific backexchange array
                back_exchange.gen_tp_backexchange_array(timepoints=time_points)

    # generate theoretical isotope distribution array
    hxrate.thr_isotope_dist_array = gen_theoretical_isotope_dist_for_all_timepoints(sequence=sequence,
                                                                                    timepoints=time_points,
                                                                                    rates=np.exp(hxrate.hx_rates),
                                                                                    inv_backexchange_array=inv_back_exchange_array,
                                                                                    d2o_fraction=d2o_fraction,
                                                                                    d2o_purity=d2o_purity,
                                                                                    num_bins=num_bins_,
                                                                                    free_energy_values=free_energy_values,
                                                                                    temperature=temperature)

    # compute the fit mse for each distribution
    hxrate.fit_rmse_each_timepoint = np.zeros(len(time_points))
    for ind, (exp_dist, thr_dist) in enumerate(zip(norm_mass_distribution_array, hxrate.thr_isotope_dist_array)):
        hxrate.fit_rmse_each_timepoint[ind] = compute_rmse_exp_thr_iso_dist(exp_isotope_dist=exp_dist,
                                                                            thr_isotope_dist=thr_dist,
                                                                            squared=False)

    # compute the total fit mse distribution
    exp_isotope_dist_concat = np.concatenate(norm_mass_distribution_array)
    thr_isotope_dist_concat = np.concatenate(hxrate.thr_isotope_dist_array)
    thr_isotope_dist_concat[np.isnan(thr_isotope_dist_concat)] = 0
    hxrate.total_fit_rmse = compute_rmse_exp_thr_iso_dist(exp_isotope_dist=exp_isotope_dist_concat,
                                                          thr_isotope_dist=thr_isotope_dist_concat,
                                                          squared=False)

    # fit gaussian to the thr isotope dist
    hxrate.thr_isotope_dist_gauss_fit = gauss_fit_to_isotope_dist_array(hxrate.thr_isotope_dist_array)

    return hxrate


def fit_rate_from_to_file(prot_name: str,
                          sequence: str,
                          hx_ms_dist_fpath: str,
                          d2o_fraction: float,
                          d2o_purity: float,
                          opt_iter: int,
                          opt_temp: float,
                          opt_step_size: float,
                          multi_proc: bool = False,
                          number_of_cores: int = 6,
                          free_energy_values: np.ndarray = None,
                          temperature: float = None,
                          usr_backexchange: float = None,
                          backexchange_corr_fpath: str = None,
                          adjust_backexchange: bool = True,
                          hx_rate_output_path: str = None,
                          hx_rate_csv_output_path: str = None,
                          hx_isotope_dist_output_path: str = None,
                          hx_rate_plot_path: str = None,
                          return_flag: bool = False) -> object:
    # todo: add param descriptions for the function here

    # get the mass distribution data at each hdx timepoint
    timepoints, mass_dist = load_data_from_hdx_ms_dist_(hx_ms_dist_fpath)

    # normalize the mass distribution
    norm_dist = normalize_mass_distribution_array(mass_dist_array=mass_dist)

    # use backexchange correction dictionary
    bkexch_corr_dict = None
    if backexchange_corr_fpath is not None:
        bkexch_corr_dict = load_backexhange_correction_(correction_fpath=backexchange_corr_fpath)

    # fit rate
    if adjust_backexchange:
        hxrate_object = fit_rate_with_backexchange_adjust(prot_name=prot_name,
                                                          sequence=sequence,
                                                          time_points=timepoints,
                                                          norm_mass_distribution_array=norm_dist,
                                                          d2o_fraction=d2o_fraction,
                                                          d2o_purity=d2o_purity,
                                                          opt_iter=opt_iter,
                                                          opt_temp=opt_temp,
                                                          opt_step_size=opt_step_size,
                                                          multi_proc=multi_proc,
                                                          number_of_cores=number_of_cores,
                                                          free_energy_values=free_energy_values,
                                                          temperature=temperature,
                                                          backexchange_value=usr_backexchange,
                                                          backexchange_correction_dict=bkexch_corr_dict)
    else:
        hxrate_object = fit_rate_without_backexchange_adjust(prot_name=prot_name,
                                                             sequence=sequence,
                                                             time_points=timepoints,
                                                             norm_mass_distribution_array=norm_dist,
                                                             d2o_fraction=d2o_fraction,
                                                             d2o_purity=d2o_purity,
                                                             opt_iter=opt_iter,
                                                             opt_temp=opt_temp,
                                                             opt_step_size=opt_step_size,
                                                             multi_proc=multi_proc,
                                                             number_of_cores=number_of_cores,
                                                             free_energy_values=free_energy_values,
                                                             temperature=temperature,
                                                             backexchange_value=usr_backexchange,
                                                             backexchange_correction_dict=bkexch_corr_dict)

    # convert hxrate object to dict and save as a pickle file

    # write hxrate as a csv file
    if hx_rate_csv_output_path is not None:
        write_hx_rate_output(hxrate_object.hx_rates, hx_rate_csv_output_path)

    # write hxrate distributions as a csv file
    if hx_isotope_dist_output_path is not None:
        write_isotope_dist_timepoints(timepoints=timepoints,
                                      isotope_dist_array=hxrate_object.thr_isotope_dist_array,
                                      output_path=hx_isotope_dist_output_path)

    # plot hxrate output as .pdf file
    if hx_rate_plot_path is not None:

        exp_centroid_arr = np.array([x.centroid for x in hxrate_object.exp_data.gauss_fit])
        exp_width_arr = np.array([x.width for x in hxrate_object.exp_data.gauss_fit])

        thr_centroid_arr = np.array([x.centroid for x in hxrate_object.thr_isotope_dist_gauss_fit])
        thr_width_arr = np.array([x.width for x in hxrate_object.thr_isotope_dist_gauss_fit])

        plot_hx_rate_fitting_(prot_name=prot_name,
                              hx_rates=hxrate_object.hx_rates,
                              exp_isotope_dist=norm_dist,
                              thr_isotope_dist=hxrate_object.thr_isotope_dist_array,
                              exp_isotope_centroid_array=exp_centroid_arr,
                              thr_isotope_centroid_array=thr_centroid_arr,
                              exp_isotope_width_array=exp_width_arr,
                              thr_isotope_width_array=thr_width_arr,
                              timepoints=timepoints,
                              fit_rmse_timepoints=hxrate_object.fit_rmse_each_timepoint,
                              fit_rmse_total=hxrate_object.total_fit_rmse,
                              backexchange=hxrate_object.back_exchange.backexchange_value,
                              backexchange_array=hxrate_object.back_exchange.backexchange_array,
                              d2o_fraction=d2o_fraction,
                              d2o_purity=d2o_purity,
                              output_path=hx_rate_plot_path)

    # save hxrate object into a pickle object
    if hx_rate_output_path is not None:
        hxrate_dict = convert_hxrate_object_to_dict(hxrate_object=hxrate_object)
        write_pickle_object(hxrate_dict, hx_rate_output_path)

    if return_flag:
        return hxrate_object


if __name__ == '__main__':
    pass

    # prot_name = 'PDB1Z96_12.30242'
    # prot_sequence = 'HMDPGLNSKIAQLVSMGFDPLEAAQALDAANGDLDVAASFLL'
    # hx_dist_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/lib15_ph6/PDB1Z96_12.30242_winner.cpickle.zlib.csv'
    # bk_corr_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/merged_data_ph6_ph7/EEHEE_rd4_0642.pdb_9.047_EEHEE_rd4_0642.pdb_9.11541/EEHEE_rd4_0642.pdb_9.047_EEHEE_rd4_0642.pdb_9.11541_merged_backexchange_correction.csv'
    #
    # fit_rate_from_to_file(prot_name=prot_name,
    #                       sequence=prot_sequence,
    #                       hx_ms_dist_fpath=hx_dist_fpath,
    #                       d2o_fraction=0.95,
    #                       d2o_purity=0.95,
    #                       opt_iter=20,
    #                       opt_temp=0.0003,
    #                       opt_step_size=0.02,
    #                       multi_proc=True,
    #                       number_of_cores=6,
    #                       backexchange_corr_fpath=bk_corr_fpath,
    #                       hx_rate_output_path=hx_dist_fpath + '_rate.pickle',
    #                       hx_rate_csv_output_path=hx_dist_fpath + '_rate.csv',
    #                       hx_rate_plot_path=hx_dist_fpath + '_rate.pdf',
    #                       hx_isotope_dist_output_path=hx_dist_fpath + '_rate_iso_dist.csv')
