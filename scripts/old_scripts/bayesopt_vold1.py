import copy
import numpy as np
from sklearn.metrics import mean_squared_error
import numpyro
from jax import random
from numpyro.infer import NUTS, MCMC
import jax.numpy as jnp
import molmass
from scipy.ndimage import center_of_mass
from scipy.optimize import curve_fit
from scipy.stats import linregress
from dataclasses import dataclass
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from itertools import cycle
import pickle

# global variables
r_constant = 0.0019872036


@dataclass
class GaussFit(object):
    """
    class container to store gauss fit data
    """
    fit_success: bool
    gauss_fit_dist: np.ndarray
    y_baseline: float
    amplitude: float
    centroid: float
    width: float
    r_sq: float
    rmse: float


@dataclass
class BayesOutputSample(object):
    posterior_samples_by_chain: np.ndarray = None
    mean: np.ndarray or float = None
    median: np.ndarray or float = None
    std: np.ndarray or float = None
    ci_5: np.ndarray or float = None
    ci_95: np.ndarray or float = None
    n_eff: np.ndarray or float = None
    rhat: np.ndarray or float = None


class ExpDataRateFit(object):

    """
    Prepare experiment data for rate fitting
    """

    def __init__(self,
                 sequence: str,
                 prot_name: str,
                 prot_rt_name: str,
                 timepoints: list or np.ndarray,
                 timepoint_index_list: list or None,
                 exp_label: str or list or None,
                 exp_distribution: list or np.ndarray,
                 backexchange: list or np.ndarray,
                 merge_exp: bool,
                 d2o_purity: float,
                 d2o_fraction: float):

        """
        init params
        :param sequence: protein sequence
        :param timepoints: timepoints array
        :param exp_distribution: ms exp distribution
        :param backexchange: backexchange array timepoints
        :param merge_exp: bool. whether to merge timepoints or not
        :param d2o_purity: d2o purity
        :param d2o_fraction: d2o fraction
        """

        self.protein_name = prot_name
        self.protein_rt_name = prot_rt_name
        self.sequence = sequence
        self.d2o_purity = d2o_purity
        self.d2o_fraction = d2o_fraction
        self.merge_exp = merge_exp

        self.timepoints = timepoints
        self.timepoint_label = gen_timepoint_label(timepoints=timepoints,
                                                   merge_exp=merge_exp,
                                                   tp_ind_label=timepoint_index_list,
                                                   exp_label=exp_label)
        self.exp_distribution = exp_distribution
        self.backexchange = backexchange

        self.num_rates = self.gen_num_exchange_rates()

        self.num_bins_ms = self.num_bins_ms_data()

        self.num_merge_facs = 0

        # update data if merging data
        # delete the first timepoint zero for list elements in 1:
        if self.merge_exp:

            self.num_merge_facs = len(self.timepoints) - 1
            self.update_backexchange_and_expdist_for_additional_exp()

            # save the input exp_dist, bkexchange, and tp list
            self.exp_distribution_list = self.exp_distribution
            self.backexchange_list = self.backexchange
            self.timepoints_list = self.timepoints

            # concat backexchange list
            self.backexchange = np.concatenate(self.backexchange)

            # concat timepoints
            # self.timepoints = np.concatenate(self.timepoints_list)

            # update exp distribution to get same number of bins
            self.check_num_bins_and_update_exp_dist(num_bins_ref=self.num_bins_ms)

            # stack exp distribution
            self.exp_distribution_stack = np.vstack(self.exp_distribution_list)

        self.flat_nonzero_exp_dist, self.nonzero_exp_dist_indices = self.gen_flat_nonzero_exp_dist()

    def gen_num_exchange_rates(self):

        temp_rates = gen_temp_rates(sequence=self.sequence)
        zero_indices = np.where(temp_rates == 0)[0]
        num_rates = len(temp_rates) - len(zero_indices)
        return num_rates

    def num_bins_ms_data(self):

        if self.merge_exp:
            num_bins = len(self.exp_distribution[0][0])
        else:
            num_bins = len(self.exp_distribution[0])

        return num_bins

    def check_num_bins_and_update_exp_dist(self, num_bins_ref: int):

        for ind, exp_dist_arr in enumerate(self.exp_distribution):

            for ind2, dist_arr in enumerate(exp_dist_arr):

                num_bins = len(exp_dist_arr[0])

                if num_bins != num_bins_ref:

                    print('Num bins unmatch for exp_distribution[%s][%s]' % (ind, ind2))

                    if num_bins < num_bins_ref:
                        pad_width = num_bins_ref - num_bins
                        print('Padding with %s 0s' % pad_width)
                        update_arr = np.pad(dist_arr, pad_width=pad_width)
                    else:
                        print('Trimming to match num bins')
                        update_arr = dist_arr[:num_bins_ref]

                    self.exp_distribution[ind][ind2] = update_arr

    def gen_flat_nonzero_exp_dist(self):

        if self.merge_exp:
            concat_arr_list = []
            for exp_arr in self.exp_distribution:
                concat_exp_arr = np.concatenate(exp_arr)
                concat_arr_list.append(concat_exp_arr)
            flat_exp_dist = np.concatenate(concat_arr_list)
        else:
            flat_exp_dist = np.concatenate(self.exp_distribution)

        non_zero_exp_dist_indices = np.nonzero(flat_exp_dist)[0]
        flat_exp_dist_non_zero = flat_exp_dist[non_zero_exp_dist_indices]

        return flat_exp_dist_non_zero, non_zero_exp_dist_indices

    def update_backexchange_and_expdist_for_additional_exp(self):

        if self.merge_exp:

            backexchange_list = []
            exp_dist_list = []
            update_tp_list = []

            for ind, (tp_arr, bkexch_arr, exp_dist_arr) in enumerate(zip(self.timepoints,
                                                                         self.backexchange,
                                                                         self.exp_distribution)):

                if ind == 0:
                    exp_dist_list.append(exp_dist_arr)
                    backexchange_list.append(bkexch_arr)
                    update_tp_list.append(tp_arr)

                else:

                    update_bkexch = []
                    update_exp_dist = []
                    update_tp_ = []

                    for ind_, tp_ in enumerate(tp_arr):
                        if tp_ != 0:
                            update_bkexch.append(bkexch_arr[ind_])
                            update_exp_dist.append(exp_dist_arr[ind_])
                            update_tp_.append(tp_)

                    update_bkexch_arr = np.array(update_bkexch)
                    update_exp_dist_arr = np.array(update_exp_dist)
                    update_tp_arr = np.array(update_tp_)

                    backexchange_list.append(update_bkexch_arr)
                    exp_dist_list.append(update_exp_dist_arr)
                    update_tp_list.append(update_tp_arr)


            # update backexchange and exp dist list
            self.backexchange = backexchange_list
            self.exp_distribution = exp_dist_list
            self.timepoints = update_tp_list

        else:
            print('del_zero_tp_backexchange: No operations to be done')


@dataclass
class RateChainDiagnostics(object):

    init_num_chains: int = None
    rmse_tol: float = 1e-2
    overall_rmse: float = None
    chain_rmse_list: list = None
    min_chain_rmse: float = None
    chain_pass_list: list = None
    discard_chain_indices: list = None
    rerun_opt: bool = False
    num_rerun_opt: int = None

def gen_timepoint_label(timepoints, merge_exp=False, tp_ind_label=None, exp_label=None):

    tp_unk_label_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

    if merge_exp:
        tp_label = [[] for _ in range(len(timepoints))]
        if tp_ind_label is None:
            if exp_label is None:
                for ind, tp_arr in enumerate(timepoints):
                    tp_label[ind] = [str(tp_unk_label_list[ind]) + '_' + str(x) for x in range(len(tp_arr))]
            else:
                for ind, tp_arr in enumerate(timepoints):
                    tp_label[ind] = [str(exp_label[ind]) + '_' + str(x) for x in range(len(tp_arr))]
        else:
            if exp_label is None:
                for ind, tp_ind_l in enumerate(tp_ind_label):
                    tp_label[ind] = [str(tp_unk_label_list[ind]) + '_' + str(x) for x in tp_ind_l]
            else:
                for ind, tp_ind_l in enumerate(tp_ind_label):
                    tp_label[ind] = [str(exp_label[ind]) + '_' + str(x) for x in tp_ind_l]
    else:
        tp_label = []
        if tp_ind_label is None:
            if exp_label is None:
                tp_label = [str(tp_unk_label_list[0]) + '_' + str(x) for x in range(len(timepoints))]
            else:
                tp_label = [str(exp_label) + '_' + str(x) for x in range(len(timepoints))]
        else:
            if exp_label is None:
                tp_label = [str(tp_unk_label_list[0]) + '_' + str(x) for x in tp_ind_label]
            else:
                tp_label = [str(exp_label) + '_' + str(x) for x in tp_ind_label]
    return tp_label


def diagnose_posterior_sample_rate_among_chains(posterior_samples_by_chain_dict,
                                                protein_sequence,
                                                timepoints,
                                                backexchange_array,
                                                d2o_fraction,
                                                d2o_purity,
                                                num_bins,
                                                exp_distribution,
                                                init_num_chains=4,
                                                rmse_tol=1e-2):

    chain_diag = RateChainDiagnostics(init_num_chains=init_num_chains,
                                      rmse_tol=rmse_tol)

    chain_diag.chain_pass_list = []
    chain_diag.chain_rmse_list = []
    chain_diag.discard_chain_indices = []

    summary = numpyro.diagnostics.summary(samples=posterior_samples_by_chain_dict)

    mean_rates_ = np.exp(summary['rate']['mean'])

    thr_iso_dist_all_chains_comb = np.array(gen_theoretical_isotope_dist_for_all_timepoints(sequence=protein_sequence,
                                                                                            timepoints=timepoints,
                                                                                            rates=np.exp(summary['rate']['mean']),
                                                                                            inv_backexchange_array=np.subtract(1, backexchange_array),
                                                                                            d2o_purity=d2o_purity,
                                                                                            d2o_fraction=d2o_fraction,
                                                                                            num_bins=num_bins))

    chain_diag.overall_rmse = compute_rmse_exp_thr_iso_dist(exp_isotope_dist=np.concatenate(exp_distribution),
                                                            thr_isotope_dist=np.concatenate(thr_iso_dist_all_chains_comb),
                                                            squared=False)

    for ind, posterior_rate_samples in enumerate(posterior_samples_by_chain_dict['rate']):

        mean_rates = np.mean(posterior_rate_samples, axis=0)

        thr_iso_dist = np.array(gen_theoretical_isotope_dist_for_all_timepoints(sequence=protein_sequence,
                                                                                timepoints=timepoints,
                                                                                rates=np.exp(mean_rates),
                                                                                inv_backexchange_array=np.subtract(1, backexchange_array),
                                                                                d2o_purity=d2o_purity,
                                                                                d2o_fraction=d2o_fraction,
                                                                                num_bins=num_bins))

        chain_rmse = compute_rmse_exp_thr_iso_dist(exp_isotope_dist=np.concatenate(exp_distribution),
                                                   thr_isotope_dist=np.concatenate(thr_iso_dist),
                                                   squared=False)

        chain_diag.chain_rmse_list.append(chain_rmse)

    chain_diag.min_chain_rmse = min(chain_diag.chain_rmse_list)

    for ind, chain_rmse_ in enumerate(chain_diag.chain_rmse_list):

        if abs(chain_rmse_ - chain_diag.min_chain_rmse) > rmse_tol:
            chain_diag.chain_pass_list.append(False)
            chain_diag.discard_chain_indices.append(ind)
        else:
            chain_diag.chain_pass_list.append(True)

    if len(chain_diag.discard_chain_indices) == init_num_chains:
        chain_diag.rerun_opt = True

    if len(chain_diag.discard_chain_indices) == 0:
        chain_diag.discard_chain_indices = None

    return chain_diag


def discard_chain_from_posterior_samples(posterior_samples_by_chain_dict,
                                         discard_chain_indices):

    new_dict = dict()

    for key_, posterior_samples_by_chain in posterior_samples_by_chain_dict.items():

        store_list = []

        for ind, posterior_samples_chain in enumerate(posterior_samples_by_chain):

            if ind not in discard_chain_indices:
                posterior_samples = np.array(posterior_samples_chain)
                store_list.append(posterior_samples)

        new_dict[key_] = np.array(store_list)

    return new_dict


class BayesRateFit(object):
    """
    Bayes Rate Fit class
    """

    def __init__(self,
                 num_chains: int = 4,
                 num_warmups: int = 100,
                 num_samples: int = 500,
                 sample_backexchange: bool = False):
        """
        initialize the class with mcmc key parameters
        :param num_chains: number of chains
        :param num_warmups: number of warmup or burn ins
        :param num_samples: number of samples for posterior
        :param return_posterior_distributions: bool. If True, returns all the posterior distributions
        :param sample_backexchange: If True, will sample backexchange for rate fitting
        """

        self.num_chains = num_chains
        self.num_warmups = num_warmups
        self.num_samples = num_samples
        self.sample_backexchange = sample_backexchange
        self.output = None

    def fit_rate(self, exp_data_object):

        # nuts kernel

        # set the number of cores to be the number of chains
        numpyro.set_host_device_count(n=self.num_chains)

        if exp_data_object.merge_exp:
            nuts_kernel = NUTS(model=rate_fit_model_norm_priors_with_merge)
        else:
            nuts_kernel = NUTS(model=rate_fit_model_norm_priors)

        # initialize MCMC
        mcmc = MCMC(nuts_kernel, num_warmup=self.num_warmups, num_samples=self.num_samples, num_chains=self.num_chains,
                    chain_method='parallel')

        # gen random key
        rng_key = random.PRNGKey(0)

        rerun_opt = True
        num_rerun_opt = -1

        while rerun_opt:

            # run mcmc

            if exp_data_object.merge_exp:
                mcmc.run(rng_key=rng_key,
                         num_rates=exp_data_object.num_rates,
                         sequence=exp_data_object.sequence,
                         timepoints_array_list=exp_data_object.timepoints,
                         num_merge_facs=len(exp_data_object.timepoints)-1,
                         backexchange_array=jnp.asarray(exp_data_object.backexchange),
                         d2o_fraction=exp_data_object.d2o_fraction,
                         d2o_purity=exp_data_object.d2o_purity,
                         num_bins=exp_data_object.num_bins_ms,
                         obs_dist_nonzero_flat=jnp.asarray(exp_data_object.flat_nonzero_exp_dist),
                         nonzero_indices=exp_data_object.nonzero_exp_dist_indices,
                         extra_fields=('potential_energy',))

            else:
                mcmc.run(rng_key=rng_key,
                         num_rates=exp_data_object.num_rates,
                         sequence=exp_data_object.sequence,
                         timepoints=jnp.asarray(exp_data_object.timepoints),
                         backexchange_array=jnp.asarray(exp_data_object.backexchange),
                         d2o_fraction=exp_data_object.d2o_fraction,
                         d2o_purity=exp_data_object.d2o_purity,
                         num_bins=exp_data_object.num_bins_ms,
                         obs_dist_nonzero_flat=jnp.asarray(exp_data_object.flat_nonzero_exp_dist),
                         nonzero_indices=exp_data_object.nonzero_exp_dist_indices,
                         extra_fields=('potential_energy',))

            # get posterior samples by chain
            posterior_samples_by_chain = sort_posterior_rates_in_samples(
                posterior_samples_by_chain=mcmc.get_samples(group_by_chain=True))

            # generate summary from the posterior samples
            summary_ = numpyro.diagnostics.summary(samples=posterior_samples_by_chain)

            if exp_data_object.merge_exp:
                if not type(summary_['merge_fac']['mean']) == np.ndarray:
                    merge_fac_mean = np.array(summary_['merge_fac']['mean'])
                else:
                    merge_fac_mean = summary_['merge_fac']['mean']

                timepoints_for_diag = np.array(recalc_timepoints_with_merge(timepoints_array_list=exp_data_object.timepoints,
                                                                            merge_facs_=merge_fac_mean))
            else:
                timepoints_for_diag = exp_data_object.timepoints

            # diagnose chains for convergence
            if exp_data_object.merge_exp:
                exp_dist_stack = exp_data_object.exp_distribution_stack
            else:
                exp_dist_stack = exp_data_object.exp_distribution

            chain_diagnostics = diagnose_posterior_sample_rate_among_chains(posterior_samples_by_chain_dict=posterior_samples_by_chain,
                                                                            protein_sequence=exp_data_object.sequence,
                                                                            timepoints=timepoints_for_diag,
                                                                            backexchange_array=exp_data_object.backexchange,
                                                                            d2o_purity=exp_data_object.d2o_purity,
                                                                            d2o_fraction=exp_data_object.d2o_fraction,
                                                                            num_bins=exp_data_object.num_bins_ms,
                                                                            exp_distribution=exp_dist_stack,
                                                                            init_num_chains=self.num_chains)

            rerun_opt = chain_diagnostics.rerun_opt

            num_rerun_opt += 1
            chain_diagnostics.num_rerun_opt = num_rerun_opt

        if chain_diagnostics.discard_chain_indices is not None:

            posterior_samples_by_chain_with_discard = discard_chain_from_posterior_samples(posterior_samples_by_chain_dict=posterior_samples_by_chain,
                                                                                           discard_chain_indices=chain_diagnostics.discard_chain_indices)

            # generate summary from the posterior samples
            summary_ = numpyro.diagnostics.summary(samples=posterior_samples_by_chain_with_discard)



        # save bayes output data to objects

        self.output = dict()

        # save some key data
        self.output['protein_name'] = exp_data_object.protein_name
        self.output['protein_rt_name'] = exp_data_object.protein_rt_name
        self.output['sequence'] = exp_data_object.sequence
        self.output['d2o_purity'] = exp_data_object.d2o_purity
        self.output['d2o_fraction'] = exp_data_object.d2o_fraction
        self.output['merge_exp'] = exp_data_object.merge_exp
        self.output['num_merge_facs'] = exp_data_object.num_merge_facs
        self.output['chain_diagnostics'] = vars(chain_diagnostics)

        # get keys from summary
        bayes_out_keys = list(summary_.keys())

        self.output['bayes_sample'] = dict()

        for ind, output_keys in enumerate(bayes_out_keys):

            self.output['bayes_sample'][output_keys] = vars(BayesOutputSample(posterior_samples_by_chain=np.array(posterior_samples_by_chain[output_keys]),
                                                                              mean=summary_[output_keys]['mean'],
                                                                              median=summary_[output_keys]['median'],
                                                                              std=summary_[output_keys]['std'],
                                                                              ci_5=summary_[output_keys]['5.0%'],
                                                                              ci_95=summary_[output_keys]['95.0%'],
                                                                              n_eff=summary_[output_keys]['n_eff'],
                                                                              rhat=summary_[output_keys]['r_hat']))

        # save exp data

        if exp_data_object.merge_exp:

            if not type(summary_['merge_fac']['mean']) == np.ndarray:
                merge_fac_mean = np.array(summary_['merge_fac']['mean'])
            else:
                merge_fac_mean = summary_['merge_fac']['mean']

            self.output['timepoints'] = np.array(recalc_timepoints_with_merge(timepoints_array_list=exp_data_object.timepoints,
                                                                              merge_facs_=merge_fac_mean))

            self.output['exp_distribution'] = exp_data_object.exp_distribution_stack
            self.output['num_merge_facs'] = exp_data_object.num_merge_facs

            self.output['tp_ind_label'] = np.concatenate(exp_data_object.timepoint_label)
            tp_lab_list = [exp_data_object.timepoint_label[0]]
            for ind, tp_lab in enumerate(exp_data_object.timepoint_label):
                if ind > 0:
                    tp_lab_list.append(tp_lab[1:])
            self.output['tp_ind_label'] = np.concatenate(tp_lab_list)

        else:
            self.output['timepoints'] = exp_data_object.timepoints
            self.output['exp_distribution'] = exp_data_object.exp_distribution
            self.output['tp_ind_label'] = np.array(exp_data_object.timepoint_label)

        self.output['timepoints_sort_indices'] = np.argsort(self.output['timepoints'])

        self.output['backexchange'] = exp_data_object.backexchange

        # calculate theoretical distribution and rmse
        self.output['pred_distribution'] = np.array(
            gen_theoretical_isotope_dist_for_all_timepoints(sequence=exp_data_object.sequence,
                                                            timepoints=jnp.asarray(self.output['timepoints']),
                                                            rates=jnp.exp(summary_['rate']['mean']),
                                                            inv_backexchange_array=jnp.subtract(1, exp_data_object.backexchange),
                                                            d2o_fraction=exp_data_object.d2o_fraction,
                                                            d2o_purity=exp_data_object.d2o_purity,
                                                            num_bins=exp_data_object.num_bins_ms))

        flat_thr_dist = np.concatenate(self.output['pred_distribution'])
        flat_thr_dist_non_zero = flat_thr_dist[exp_data_object.nonzero_exp_dist_indices]

        # calculate rmse
        self.output['rmse'] = dict()
        self.output['rmse']['total'] = compute_rmse_exp_thr_iso_dist(exp_isotope_dist=exp_data_object.flat_nonzero_exp_dist,
                                                                     thr_isotope_dist=flat_thr_dist_non_zero,
                                                                     squared=False)

        self.output['rmse']['per_timepoint'] = np.zeros(len(self.output['timepoints']))

        for ind, (exp_dist, thr_dist) in enumerate(zip(self.output['exp_distribution'], self.output['pred_distribution'])):

            self.output['rmse']['per_timepoint'][ind] = compute_rmse_exp_thr_iso_dist(exp_isotope_dist=exp_dist,
                                                                                      thr_isotope_dist=thr_dist,
                                                                                      squared=False)

        # fit gaussian to exp and pred data
        self.output['exp_dist_gauss_fit'] = [vars(x) for x in gauss_fit_to_isotope_dist_array(isotope_dist=self.output['exp_distribution'])]
        self.output['pred_dist_guass_fit'] = [vars(x) for x in gauss_fit_to_isotope_dist_array(isotope_dist=self.output['pred_distribution'])]

    def exp_merge_dist_to_csv(self, output_path: str):

        if self.output is not None:

            sort_tp = self.output['timepoints'][self.output['timepoints_sort_indices']]
            sort_expdata = self.output['exp_distribution'][self.output['timepoints_sort_indices']]

            write_isotope_dist_timepoints(timepoints=sort_tp,
                                          isotope_dist_array=sort_expdata,
                                          output_path=output_path)
        else:
            print('Output is None')

    def pred_dist_to_csv(self, output_path: str):

        if self.output is not None:

            sort_tp = self.output['timepoints'][self.output['timepoints_sort_indices']]
            sort_tp_labels = self.output['tp_ind_label'][self.output['timepoints_sort_indices']]
            sort_tp_ind = [x.split('_')[-1] for x in sort_tp_labels]
            sort_pred_dist = self.output['pred_distribution'][self.output['timepoints_sort_indices']]

            write_isotope_dist_timepoints(timepoints=sort_tp,
                                          isotope_dist_array=sort_pred_dist,
                                          output_path=output_path,
                                          timepoint_label=sort_tp_ind)

        else:
            print('Output is None')

    def write_rates_to_csv(self, output_path: str):

        if self.output is not None:

            write_hx_rate_output_csv(hxrate_mean_array=self.output['bayes_sample']['rate']['mean'],
                                     hxrate_median_array=self.output['bayes_sample']['rate']['median'],
                                     hxrate_std_array=self.output['bayes_sample']['rate']['std'],
                                     hxrate_5percent_array=self.output['bayes_sample']['rate']['ci_5'],
                                     hxrate_95percent_array=self.output['bayes_sample']['rate']['ci_95'],
                                     neff_array=self.output['bayes_sample']['rate']['n_eff'],
                                     r_hat_array=self.output['bayes_sample']['rate']['rhat'],
                                     output_path=output_path)
        else:
            print('Output is None')

    def plot_hxrate_output(self, output_path: str):

        if self.output is not None:

            # do some operations for correct input

            # calculate hxrate_pipeline error based on CI
            hxrate_error = np.zeros((2, len(self.output['bayes_sample']['rate']['mean'])))
            hxrate_error[0] = np.subtract(self.output['bayes_sample']['rate']['mean'], self.output['bayes_sample']['rate']['ci_5'])
            hxrate_error[1] = np.subtract(self.output['bayes_sample']['rate']['ci_95'], self.output['bayes_sample']['rate']['mean'])

            # sort dist with timepoints
            sort_tp = self.output['timepoints'][self.output['timepoints_sort_indices']]
            sort_tp_label = self.output['tp_ind_label'][self.output['timepoints_sort_indices']]
            sort_exp_dist = self.output['exp_distribution'][self.output['timepoints_sort_indices']]
            sort_pred_dist = self.output['pred_distribution'][self.output['timepoints_sort_indices']]

            exp_centroid_arr_ = np.array([x['centroid'] for x in self.output['exp_dist_gauss_fit']])[self.output['timepoints_sort_indices']]
            pred_centroid_arr_ = np.array([x['centroid'] for x in self.output['pred_dist_guass_fit']])[self.output['timepoints_sort_indices']]

            exp_width_arr_ = np.array([x['width'] for x in self.output['exp_dist_gauss_fit']])[self.output['timepoints_sort_indices']]
            pred_width_arr_ = np.array([x['width'] for x in self.output['pred_dist_guass_fit']])[self.output['timepoints_sort_indices']]

            rmse_per_timepoint = self.output['rmse']['per_timepoint'][self.output['timepoints_sort_indices']]

            backexchange_array = self.output['backexchange'][self.output['timepoints_sort_indices']]

            if self.output['merge_exp']:
                num_merge_facs = self.output['num_merge_facs']
                merge_fac = self.output['bayes_sample']['merge_fac']['mean']
            else:
                num_merge_facs = None
                merge_fac = None

            plot_hx_rate_fitting_bayes(prot_name=self.output['protein_name'],
                                       hx_rates=self.output['bayes_sample']['rate']['mean'],
                                       hx_rates_error=hxrate_error,
                                       timepoints=sort_tp,
                                       timepoint_label=sort_tp_label,
                                       exp_isotope_dist=sort_exp_dist,
                                       thr_isotope_dist=sort_pred_dist,
                                       exp_isotope_centroid_array=exp_centroid_arr_,
                                       thr_isotope_centroid_array=pred_centroid_arr_,
                                       exp_isotope_width_array=exp_width_arr_,
                                       thr_isotope_width_array=pred_width_arr_,
                                       fit_rmse_timepoints=rmse_per_timepoint,
                                       fit_rmse_total=self.output['rmse']['total'],
                                       backexchange=self.output['backexchange'][-1],
                                       backexchange_array=backexchange_array,
                                       d2o_fraction=self.output['d2o_fraction'],
                                       d2o_purity=self.output['d2o_purity'],
                                       merge_exp=self.output['merge_exp'],
                                       num_merge_factors=num_merge_facs,
                                       merge_factor=merge_fac,
                                       output_path=output_path)
        else:
            print('Output is None')

    def plot_bayes_samples(self, output_path: str):

        if self.output is not None:

            plot_posteriors(bayesfit_sample_dict=self.output['bayes_sample'],
                            discard_chain_bool_list=self.output['chain_diagnostics']['chain_pass_list'],
                            output_path=output_path)

        else:
            print('Output is None')

    def output_to_pickle(self, output_path: str, save_posterior_samples: bool):

        output_dict = copy.deepcopy(self.output)

        if not save_posterior_samples:
            # delete posterior sample and save the dictionary
            for dicts_ in output_dict['bayes_sample'].values():
                dicts_.pop('posterior_samples_by_chain', None)

        write_pickle_object(obj=output_dict,
                            filepath=output_path)


def gen_temp_rates(sequence):

    # set high rate value as 1e2
    high_rate_value = 1e2
    rates = np.array([high_rate_value] * len(sequence), dtype=float)

    # set the rates for the first two residues as 0
    rates[:2] = 0

    # set the rate for proline to be 0
    if 'P' in sequence:
        amino_acid_list = [x for x in sequence]
        for ind, amino_acid in enumerate(amino_acid_list):
            if amino_acid == 'P':
                rates[ind] = 0

    return rates


def write_pickle_object(obj, filepath):
    """
    write an object to a pickle file
    :param obj: object
    :param filepath: pickle file path
    :return: None
    """
    with open(filepath, 'wb') as outfile:
        pickle.dump(obj, outfile)


def reshape_posterior_samples(posterior_samples):
    """
    reshape the posterior samples nd array
    :param posterior_samples:
    :return:
    """

    posterior_samples_shape = posterior_samples.shape

    reshape_array = np.zeros((posterior_samples_shape[-1], posterior_samples_shape[0], posterior_samples_shape[1]))

    for chain_ind, chain_arrs in enumerate(posterior_samples):
        chain_arr_transpose = chain_arrs.T
        for posterior_ind, posterior_arr in enumerate(chain_arr_transpose):
            for sample_ind, sample_value in enumerate(posterior_arr):
                reshape_array[posterior_ind][chain_ind][sample_ind] = sample_value

    return reshape_array


def plot_posteriors(bayesfit_sample_dict, discard_chain_bool_list=[False, False, False, False], output_path=None):

    num_fig_grids = 0

    for dict_items in bayesfit_sample_dict.values():
        mean = dict_items['mean']
        if type(mean) == np.float32:
            num_fig_grids += 1
        else:
            num_fig_grids += len(mean)

    # number of columns 2
    num_columns = 2

    num_rows = num_fig_grids/num_columns

    if num_fig_grids % 2 != 0:
        num_rows = divmod(num_fig_grids, num_columns)[0] + 1

    num_rows = int(num_rows)

    font_size = 8

    fig_size = (20, 2.5 * num_rows)

    fig = plt.figure(figsize=fig_size)

    plt.rcParams.update({'font.size': font_size})

    gs0 = gridspec.GridSpec(nrows=num_rows, ncols=num_columns)

    # do one when there's only one mean, if more than one then do another action

    counter = 0

    for keys_, items_ in bayesfit_sample_dict.items():

        if type(items_['mean']) == np.float32:

            plot_posteriors_grid(fig_obj=fig,
                                 sample=items_['posterior_samples_by_chain'],
                                 sample_mean=items_['mean'],
                                 sample_std=items_['std'],
                                 sample_5percent=items_['ci_5'],
                                 sample_95percent=items_['ci_95'],
                                 sample_rhat=items_['rhat'],
                                 sample_label=keys_,
                                 gridspec_obj=gs0,
                                 gridspec_index=counter,
                                 discard_chain_bool_list=discard_chain_bool_list)

            counter += 1

        else:

            num_sample_per_key = len(items_['mean'])

            reshape_array = reshape_posterior_samples(posterior_samples=items_['posterior_samples_by_chain'])

            for num in range(num_sample_per_key):

                plot_posteriors_grid(fig_obj=fig,
                                     sample=reshape_array[num],
                                     sample_mean=items_['mean'][num],
                                     sample_std=items_['std'][num],
                                     sample_5percent=items_['ci_5'][num],
                                     sample_95percent=items_['ci_95'][num],
                                     sample_rhat=items_['rhat'][num],
                                     sample_label='%s_%s' % (keys_, num),
                                     gridspec_obj=gs0,
                                     gridspec_index=counter,
                                     discard_chain_bool_list=discard_chain_bool_list)

                counter += 1

    plt.subplots_adjust(hspace=1.2, wspace=0.12, top=0.95)

    plt.savefig(output_path, bbox_inches="tight")
    plt.close()


def plot_posteriors_grid(fig_obj,
                         sample,
                         sample_mean,
                         sample_std,
                         sample_5percent,
                         sample_95percent,
                         sample_rhat,
                         sample_label,
                         gridspec_obj,
                         gridspec_index,
                         discard_chain_bool_list=[False,False,False,False]):
    """

    :param fig_obj:
    :param sample:
    :param sample_mean:
    :param sample_std:
    :param sample_5percent:
    :param sample_95percent:
    :param sample_rhat:
    :param sample_label:
    :param gridspec_obj:
    :param gridspec_index:
    :return:
    """

    gs00 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gridspec_obj[gridspec_index])

    ax00 = fig_obj.add_subplot(gs00[0, :-1])

    chain_colors = cycle(['indianred', 'mediumaquamarine', 'deepskyblue', 'mediumpurple', 'palevioletred'])

    discard_chain_color = 'gray'

    ind_arr = np.arange(0, len(sample[0]))
    for ind, (sample_per_chain, discard_chain, color_) in enumerate(zip(sample, discard_chain_bool_list, chain_colors)):
        if discard_chain:
            color_ = discard_chain_color
        if ind == 0:
            ax00.plot(ind_arr, sample_per_chain, color=color_, linewidth=0.5)
        else:
            last_ind = ind_arr[-1]
            ind_arr = np.arange(last_ind+1, last_ind + 1 + len(sample_per_chain))
            ax00.plot(ind_arr, sample_per_chain, color=color_, linewidth=0.5)

    ax00.hlines(y=sample_mean, xmin=0, xmax=ind_arr[-1], ls='--', colors='red', linewidth=0.5)
    ax00.hlines(y=sample_5percent, xmin=0, xmax=ind_arr[-1], ls='--', colors='black', linewidth=0.5)
    ax00.hlines(y=sample_95percent, xmin=0, xmax=ind_arr[-1], ls='--', colors='black', linewidth=0.5)
    ax00.spines['right'].set_visible(False)
    ax00.spines['top'].set_visible(False)
    plt.grid(alpha=0.25)
    ax00.tick_params(length=3, pad=3)

    # put stats on the right side
    plt.text(0.95,
             1.2,
             "mean = %.4f\nstd = %.4f" % (sample_mean, sample_std),
             horizontalalignment='right',
             verticalalignment='top',
             transform=ax00.transAxes)

    # put the label on the left side
    plt.text(0.01, 1.2, sample_label,
             horizontalalignment="left",
             verticalalignment="top",
             transform=ax00.transAxes)

    ax01 = fig_obj.add_subplot(gs00[0, -1])
    ax01.scatter(0, sample_rhat, color='black')
    ax01.hlines(y=1+0.05, xmin=-1, xmax=1, ls='--', colors='black', linewidth=1)
    ax01.hlines(y=1-0.05, xmin=-1, xmax=1, ls='--', colors='black', linewidth=1)
    ax01.spines['right'].set_visible(False)
    ax01.spines['top'].set_visible(False)
    ax01.spines['bottom'].set_visible(False)
    plt.xticks([])
    ax01.tick_params(length=3, pad=3)

    # put rhat on the right side
    plt.text(0.95,
             1.2,
             "rhat = %.4f" % sample_rhat,
             horizontalalignment='right',
             verticalalignment='top',
             transform=ax01.transAxes)


def correct_centroids_using_backexchange(centroids: np.ndarray,
                                         backexchange_array: np.ndarray,
                                         include_zero_dist: bool = False) -> np.ndarray:
    """

    :param centroids: uncorrected centroids
    :param backexchange_array: backexchange array for each timepoint
    :return: corrected centroids
    """

    # generate corr centroids zero arrays and fill with corrections. timepoint 0 doesn't need any correction
    corr_centroids = np.zeros(len(centroids))

    # apply correction to centroids based on backexchange
    for ind, (centr, bkexch) in enumerate(zip(centroids, backexchange_array)):
        corr_centroids[ind] = centr/(1-bkexch)

    if include_zero_dist:
        corr_centroids[0] = centroids[0]

    return corr_centroids


def plot_hx_rate_fitting_bayes(prot_name: str,
                               hx_rates: np.ndarray,
                               hx_rates_error: np.ndarray,
                               exp_isotope_dist: np.ndarray,
                               thr_isotope_dist: np.ndarray,
                               exp_isotope_centroid_array: np.ndarray,
                               thr_isotope_centroid_array: np.ndarray,
                               exp_isotope_width_array: np.ndarray,
                               thr_isotope_width_array: np.ndarray,
                               timepoints: np.ndarray,
                               timepoint_label: np.ndarray,
                               fit_rmse_timepoints: np.ndarray,
                               fit_rmse_total: float,
                               backexchange: float,
                               backexchange_array: np.ndarray,
                               d2o_fraction: float,
                               d2o_purity: float,
                               output_path: str,
                               merge_exp: bool = False,
                               num_merge_factors: int = None,
                               merge_factor: np.ndarray=None):
    """
    generate several plots for visualizing the hx rate fitting output
    :param prot_name: protein name
    :param hx_rates: in ln scale
    :param exp_isotope_dist: exp isotope dist array
    :param thr_isotope_dist: thr isotope dist array from hx rates
    :param exp_isotope_centroid_array: exp isotope centroid in an array
    :param thr_isotope_centroid_array: thr isotope centroid in an array
    :param exp_isotope_width_array: exp isotope width in an array
    :param thr_isotope_width_array: thr isotope width in an array
    :param timepoints: time points
    :param fit_rmse_timepoints: fit rmse for each timepoint
    :param fit_rmse_total: total fit rmse
    :param backexchange: backexchange value
    :param backexchange_array: backexchange array
    :param d2o_fraction: d2o fraction
    :param d2o_purity: d2o purity
    :param output_path: plot saveing output path
    :return:
    """

    # define figure size
    num_columns = 2

    # this is a constant
    num_plots_second_row = 8

    num_rows = len(timepoints)
    row_width = 1.0
    if len(timepoints) < num_plots_second_row:
        num_rows = num_plots_second_row
        row_width = 2.5
    fig_size = (25, row_width * num_rows)

    font_size = 10

    fig = plt.figure(figsize=fig_size)
    gs = fig.add_gridspec(nrows=num_rows, ncols=num_columns)

    plt.rcParams.update({'font.size': font_size})

    if fit_rmse_timepoints is None:
        fit_rmse_tp = np.zeros(len(timepoints))
    else:
        fit_rmse_tp = fit_rmse_timepoints

    # def the color for different experiments
    num_exp = 1
    if merge_exp:
        num_exp = num_merge_factors + 1

    color_exp = ['dodgerblue', 'forestgreen', 'fuchsia', 'teal', 'gold', 'sienna']

    exp_labels = [x.split('_')[0] for x in timepoint_label]
    # use this in title for indicating the exp
    unique_exp_labels = np.unique(exp_labels)
    unique_exp_color_dict = dict()
    for ind, uniq_exp in enumerate(unique_exp_labels):
        unique_exp_color_dict[uniq_exp] = color_exp[ind]

    color_exp_list = [unique_exp_color_dict[x] for x in exp_labels]

    #######################################################
    #######################################################
    # start plotting the exp and thr isotope dist
    for num, (timepoint, exp_dist, thr_dist, exp_centroid, bkexch, tp_label, exp_label) in enumerate(zip(timepoints, exp_isotope_dist, thr_isotope_dist, exp_isotope_centroid_array, backexchange_array, timepoint_label, exp_labels)):

        if fit_rmse_timepoints is None:
            rmse_tp = compute_rmse_exp_thr_iso_dist(exp_isotope_dist=exp_dist,
                                                    thr_isotope_dist=thr_dist,
                                                    squared=False)
            fit_rmse_tp[num] = rmse_tp
        else:
            rmse_tp = fit_rmse_tp[num]

        # plot exp and thr isotope dist
        ax = fig.add_subplot(gs[num, 0])
        plt.plot(exp_dist, color=unique_exp_color_dict[exp_label], marker='o', ls='-', markersize=3)
        plt.plot(thr_dist, color='red')
        ax.set_yticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # ax.spines['left'].set_visible(False)
        plt.xticks(range(0, len(exp_dist) + 5, 5))
        ax.set_xticklabels(range(0, len(exp_dist) + 5, 5))
        plt.grid(axis='x', alpha=0.25)
        ax.tick_params(length=3, pad=3)

        # put the rmse on the right side of the plot and delta centroid
        if num == 0:
            delta_centroid = 0.0
        else:
            prev_centroid = exp_isotope_centroid_array[num - 1]
            delta_centroid = exp_centroid - prev_centroid

        # delta_centroid_text = 'dmz=%.2f' % delta_centroid

        plt.text(1.0, 1.2, "fit rmse = %.4f\nd_mz = %.2f\nbkexch = %.2f" % (rmse_tp, delta_centroid, bkexch*100),
                 horizontalalignment="right",
                 verticalalignment="top",
                 transform=ax.transAxes)

        # put timepoint on  the left side of the plot
        plt.text(0.01, 1.2, "%s\n%s" % (tp_label, timepoint),
                 horizontalalignment="left",
                 verticalalignment="top",
                 transform=ax.transAxes)

        # put the centroid information by the peak max
        plt.text(
            exp_centroid,
            1.1,
            "%.1f" % exp_centroid,
            horizontalalignment="center",
            verticalalignment="bottom",
            fontsize=8)
    #######################################################
    #######################################################

    # 8 plots on the second row
    second_plot_row_thickness = int(len(timepoints)/num_plots_second_row)
    if len(timepoints) < num_plots_second_row:
        second_plot_row_thickness = 1
    second_plot_indices = [(num*second_plot_row_thickness) for num in range(num_plots_second_row)]

    #######################################################
    #######################################################
    # plot timepoint specific backexchange
    ax0 = fig.add_subplot(gs[second_plot_indices[0]: second_plot_indices[1], 1])

    plt.scatter(x=np.arange(len(timepoints))[1:], y=backexchange_array[1:]*100, color=color_exp_list[1:])
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    plt.xticks(range(-1, len(timepoints) + 1, 1))
    ax0.set_xticklabels(range(-1, len(timepoints) + 1, 1))
    plt.grid(axis='x', alpha=0.25)
    plt.grid(axis='y', alpha=0.25)
    plt.xlabel('Timepoint index')
    plt.ylabel('Back Exchange (%)')
    ax0.tick_params(length=3, pad=3)

    #######################################################
    #######################################################
    # plot fit rmse

    ax1 = fig.add_subplot(gs[second_plot_indices[1]: second_plot_indices[2], 1])
    plt.scatter(np.arange(len(timepoints)), fit_rmse_tp, color=color_exp_list)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    plt.xticks(range(0, len(timepoints) + 1, 1))
    ax1.set_xticklabels(range(0, len(timepoints) + 1, 1))
    if max(fit_rmse_tp) <= 0.15:
        y_ticks = np.round(np.linspace(0, 0.15, num=16), 2)
        plt.yticks(y_ticks)
        ax1.set_yticklabels(y_ticks)
    else:
        plt.axhline(y=0.15, ls='--', color='black')
    plt.grid(axis='x', alpha=0.25)
    plt.grid(axis='y', alpha=0.25)
    plt.xlabel('Timepoint index')
    plt.ylabel('Fit RMSE')
    ax1.tick_params(length=3, pad=3)

    #######################################################
    #######################################################

    #######################################################
    #######################################################
    # plot center of mass exp and thr

    timepoints_v2 = np.array([x for x in timepoints])
    timepoints_v2[0] = timepoints_v2[2] - timepoints_v2[1]

    ax2 = fig.add_subplot(gs[second_plot_indices[2]: second_plot_indices[3], 1])
    ax2.plot(timepoints_v2, exp_isotope_centroid_array, ls='--', color='gray')
    ax2.plot(timepoints_v2, thr_isotope_centroid_array, ls='--', color='gray')
    plt.scatter(timepoints_v2, exp_isotope_centroid_array, marker='o', ls='-', color=color_exp_list)
    plt.scatter(timepoints_v2, thr_isotope_centroid_array, marker='o', ls='-', color='red')
    ax2.set_xscale('log')
    ax2.set_xticks(timepoints_v2)
    ax2.set_xticklabels([])
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    plt.grid(axis='x', alpha=0.25)
    plt.grid(axis='y', alpha=0.25)
    plt.xlabel('log(timepoint)')
    plt.ylabel('Centroid')
    #######################################################
    #######################################################

    #######################################################
    #######################################################
    # plot center of mass exp and thr corrected using backexchange

    exp_isotope_centroid_array_corr = correct_centroids_using_backexchange(centroids=exp_isotope_centroid_array,
                                                                           backexchange_array=backexchange_array,
                                                                           include_zero_dist=True)
    thr_isotope_centroid_array_corr = correct_centroids_using_backexchange(centroids=thr_isotope_centroid_array,
                                                                           backexchange_array=backexchange_array,
                                                                           include_zero_dist=True)

    ax2 = fig.add_subplot(gs[second_plot_indices[3]: second_plot_indices[4], 1])
    ax2.plot(timepoints_v2, exp_isotope_centroid_array_corr, ls='--', color='gray')
    ax2.plot(timepoints_v2, thr_isotope_centroid_array_corr, ls='--', color='gray')
    plt.scatter(timepoints_v2, exp_isotope_centroid_array_corr, marker='o', ls='-', color=color_exp_list)
    plt.scatter(timepoints_v2, thr_isotope_centroid_array_corr, marker='o', ls='-', color='red')

    ax2.set_xscale('log')
    ax2.set_xticks(timepoints_v2)
    ax2.set_xticklabels([])
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    plt.grid(axis='x', alpha=0.25)
    plt.grid(axis='y', alpha=0.25)
    plt.xlabel('log(timepoint)')
    plt.ylabel('Corrected Centroid')
    #######################################################
    #######################################################

    #######################################################
    #######################################################
    # plot the error in center of mass between exp and thr distributions

    com_difference = np.subtract(exp_isotope_centroid_array_corr, thr_isotope_centroid_array_corr)

    ax3 = fig.add_subplot(gs[second_plot_indices[4]: second_plot_indices[5], 1])
    ax3.scatter(np.arange(len(timepoints)), com_difference, color=color_exp_list)
    plt.axhline(y=0, ls='--', color='black', alpha=0.50)
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    plt.xticks(range(0, len(timepoints) + 1, 1))
    ax3.set_xticklabels(range(0, len(timepoints) + 1, 1))
    plt.grid(axis='x', alpha=0.25)
    plt.grid(axis='y', alpha=0.25)
    plt.xlabel('Timepoint index')
    plt.ylabel('Centroid difference (E-T)')
    #######################################################
    #######################################################

    #######################################################
    #######################################################
    # plot the width of exp and thr distributions

    ax4 = fig.add_subplot(gs[second_plot_indices[5]: second_plot_indices[6], 1])
    ax4.scatter(np.arange(len(timepoints)), exp_isotope_width_array, color=color_exp_list)
    ax4.scatter(np.arange(len(timepoints)), thr_isotope_width_array, color='red')
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    plt.xticks(range(0, len(timepoints) + 1, 1))
    ax4.set_xticklabels(range(0, len(timepoints) + 1, 1))
    plt.grid(axis='x', alpha=0.25)
    plt.grid(axis='y', alpha=0.25)
    plt.xlabel('Timepoint index')
    plt.ylabel('Width')
    #######################################################
    #######################################################

    #######################################################
    #######################################################
    # plot the error in center of mass between exp and thr distributions

    width_difference = np.subtract(exp_isotope_width_array, thr_isotope_width_array)

    ax5 = fig.add_subplot(gs[second_plot_indices[6]: second_plot_indices[7], 1])
    ax5.scatter(np.arange(len(timepoints)), width_difference, color=color_exp_list)
    plt.axhline(y=0, ls='--', color='black', alpha=0.50)
    ax5.spines['right'].set_visible(False)
    ax5.spines['top'].set_visible(False)
    plt.xticks(range(0, len(timepoints) + 1, 1))
    ax5.set_xticklabels(range(0, len(timepoints) + 1, 1))
    plt.grid(axis='x', alpha=0.25)
    plt.grid(axis='y', alpha=0.25)
    plt.xlabel('Timepoint index')
    plt.ylabel('Width difference (E-T)')
    #######################################################
    #######################################################

    #######################################################
    #######################################################
    # plot the rates in log10 scale
    hx_rates_log10 = np.log10(np.exp(hx_rates))
    hx_rates_error_log10 = np.log10(np.exp(hx_rates_error))

    ax6 = fig.add_subplot(gs[second_plot_indices[7]:, 1])
    plt.errorbar(x=np.arange(len(hx_rates_log10)), y=hx_rates_log10, yerr=hx_rates_error_log10, marker='o', ls='-',
                 color='red', markerfacecolor='red', markeredgecolor='black')
    # plt.plot(np.arange(len(hx_rates_log10)), hx_rates_log10[sort_ind], marker='o', ls='-', color='red',
    #          markerfacecolor='red', markeredgecolor='black')
    plt.xticks(range(0, len(hx_rates_log10) + 2, 2))
    ax6.set_xticklabels(range(0, len(hx_rates_log10) + 2, 2))
    ax6.spines['right'].set_visible(False)
    ax6.spines['top'].set_visible(False)
    plt.grid(axis='x', alpha=0.25)
    plt.grid(axis='y', alpha=0.25)
    plt.xlabel('Residues (Ranked from slowest to fastest exchanging)')
    plt.ylabel('Rate: log k (1/s)')
    #######################################################
    #######################################################

    # adjust some plot properties and add title
    plt.subplots_adjust(hspace=1.2, wspace=0.1, top=0.96)

    merge_factor_str = 'One Exp'

    if merge_exp:

        merge_fac_str_list = ['merge_fac_' + str(ind) + ': ' + str(np.round(x, 4)) for ind, x in enumerate(merge_factor)]

        if len(merge_fac_str_list) > 1:
            merge_factor_str = ' | '.join([x for x in merge_fac_str_list])
        else:
            merge_factor_str = merge_fac_str_list[0]

    title_1 = 'Fit RMSE: %.4f | %s | Backexchange: %.2f %% | D2O Purity: %.1f %% | D2O_Fraction: %.1f %%' %(fit_rmse_total,
                                                                                                            merge_factor_str,
                                                                                                            backexchange*100,
                                                                                                            d2o_purity*100,
                                                                                                            d2o_fraction*100)

    plot_title = prot_name + ' (' + title_1 + ')'

    plt.suptitle(plot_title)

    # title set according to unique exp color dict

    incr_x_ind = 0.020
    curr_x_loc = 0.495

    curr_y_loc = 0.968

    for ind, (key, values) in enumerate(unique_exp_color_dict.items()):
        curr_x_loc = curr_x_loc + (ind * incr_x_ind)
        plt.figtext(curr_x_loc, curr_y_loc, key, color=values, ha='right', fontsize=8)

    plt.figtext(curr_x_loc + incr_x_ind, curr_y_loc, "FIT", color='red', ha='right', fontsize=8)

    # plt.figtext(0.498, 0.968, "EXP", color='blue', ha='right', fontsize=8)
    # plt.figtext(0.502, 0.968, "FIT", color='red', ha='left', fontsize=8)

    plt.savefig(output_path, bbox_inches="tight")
    plt.close()


def write_isotope_dist_timepoints(timepoints, isotope_dist_array, output_path, timepoint_label=None):

    if timepoint_label is None:
        timepoint_label = np.arange(0, len(timepoints))

    tp_label_str = ','.join([str(x) for x in timepoint_label])
    header1 = '#tp_ind,' + tp_label_str + '\n'

    timepoint_str = ','.join(['%.4f' % x for x in timepoints])
    header2 = '#tp,' + timepoint_str + '\n'
    data_string = ''
    for ind, arr in enumerate(isotope_dist_array.T):
        arr_str = ','.join([str(x) for x in arr])
        data_string += '{},{}\n'.format(ind, arr_str)

    with open(output_path, 'w') as outfile:
        outfile.write(header1 + header2 + data_string)
        outfile.close()


def write_hx_rate_output_csv(hxrate_mean_array: np.ndarray,
                             hxrate_median_array: np.ndarray,
                             hxrate_std_array: np.ndarray,
                             hxrate_5percent_array: np.ndarray,
                             hxrate_95percent_array: np.ndarray,
                             neff_array: np.ndarray,
                             r_hat_array: np.ndarray,
                             output_path: str):

    header = 'ind,rate_mean,rate_median,rate_std,rate_5%,rate_95%,n_eff,r_hat\n'
    data_string = ''

    for num in range(len(hxrate_mean_array)):

        data_string += '{},{},{},{},{},{},{},{}\n'.format(num,
                                                          hxrate_mean_array[num],
                                                          hxrate_median_array[num],
                                                          hxrate_std_array[num],
                                                          hxrate_95percent_array[num],
                                                          hxrate_5percent_array[num],
                                                          neff_array[num],
                                                          r_hat_array[num])

    with open(output_path, 'w') as outfile:
        outfile.write(header + data_string)
        outfile.close()


def estimate_gauss_param(ydata: np.ndarray,
                         xdata: np.ndarray,
                         baseline: float = 0.0,
                         width_ht: float = 0.7) -> list:
    ymax = np.max(ydata)
    maxindex = np.nonzero(ydata == ymax)[0]
    peakmax_x = xdata[maxindex][0]
    norm_arr = ydata/max(ydata)
    bins_for_width = norm_arr[norm_arr > width_ht]
    width_bin = len(bins_for_width)
    init_guess = [baseline, ymax, peakmax_x, width_bin]
    return init_guess


def gauss_func(x, y0, A, xc, w):
    """
    gaussian function with baseline
    :param x: xdata
    :param y0: baseline
    :param A: amplitude
    :param xc: centroid
    :param w: width
    :return: gauss(x)
    """
    rxc = ((x - xc) ** 2) / (2 * (w ** 2))
    y = y0 + A * (np.exp(-rxc))
    return y


def fit_gaussian(data: np.ndarray) -> object:
    """
    fit gaussian to data
    :param data: xdata to fit gaussian
    :return: gauss fit object
    """
    xdata = np.arange(len(data))
    guess_params = estimate_gauss_param(ydata=data,
                                        xdata=xdata)

    # initialize gauss fit object with fit success as false
    mean = sum(xdata * data) / sum(data)
    sigma = np.sqrt(sum(data * (xdata - mean) ** 2) / sum(data))
    gaussfit = GaussFit(fit_success=False,
                        gauss_fit_dist=data,
                        y_baseline=guess_params[0],
                        amplitude=guess_params[1],
                        centroid=center_of_mass_(data_array=data),
                        width=sigma,
                        r_sq=0.00,
                        rmse=100.0)

    if np.all(data == data[0]):
        return gaussfit
    else:
        try:

            # fit gaussian
            popt, pcov = curve_fit(gauss_func, xdata, data, p0=guess_params, maxfev=100000)

            # if the centroid is smaller than 0, return the false gaussfit object
            if popt[2] < 0.0:
                return gaussfit

            # if the width is smaller than 0, return the false gauss fit object
            if popt[3] < 0.0 or popt[3] > len(data):
                return gaussfit

            # for successful gaussian fit
            else:
                gaussfit.fit_success = True
                gaussfit.y_baseline = popt[0]
                gaussfit.amplitude = popt[1]
                gaussfit.centroid = popt[2]
                gaussfit.width = popt[3]
                gaussfit.gauss_fit_dist = gauss_func(xdata, *popt)
                gaussfit.rmse = compute_rmse_exp_thr_iso_dist(exp_isotope_dist=data,
                                                              thr_isotope_dist=gaussfit.gauss_fit_dist,
                                                              squared=True)
                slope, intercept, rvalue, pvalue, stderr = linregress(data, gaussfit.gauss_fit_dist)
                gaussfit.r_sq = rvalue**2

                return gaussfit

        except RuntimeError:
            return gaussfit


def gauss_fit_to_isotope_dist_array(isotope_dist: np.ndarray) -> list:
    """
    fit gaussian to each isotope dist array and store the output in a list
    :param isotope_dist: 2d array of isotope dist
    :return: gauss fit list
    """

    gauss_fit_list = []

    for dist in isotope_dist:

        gauss_fit = fit_gaussian(data=dist)
        gauss_fit_list.append(gauss_fit)

    return gauss_fit_list


def center_of_mass_(data_array):
    com = center_of_mass(data_array)[0]
    return com


def replace_nans_with_zeros_1d_array(array):
    """
    replace nans with 0s in 1 d array
    :param array:
    :return:
    """
    # check for nans and replace by 0s
    nan_inds = np.argwhere(np.isnan(array))
    if len(nan_inds) > 0:
        array[nan_inds] = 0
    return array


def compute_rmse_exp_thr_iso_dist(exp_isotope_dist: np.ndarray,
                                  thr_isotope_dist: np.ndarray,
                                  squared: bool = False):
    """
    compute the mean squared error between exp and thr isotope distribution only with values of exp_dist that are higher
    than 0
    :param exp_isotope_dist:
    :param thr_isotope_dist:
    :param squared:
    :return:
    """

    # replace nans with 0s
    exp_isotope_dist_nonans = replace_nans_with_zeros_1d_array(exp_isotope_dist)
    thr_isotope_dist_nonans = replace_nans_with_zeros_1d_array(thr_isotope_dist)

    # get the array with only values greater than 0 for the exp data
    exp_isotope_dist_comp = exp_isotope_dist_nonans[exp_isotope_dist_nonans > 0]
    thr_isotope_dist_comp = thr_isotope_dist_nonans[exp_isotope_dist_nonans > 0]

    if len(exp_isotope_dist_comp) > 0 and len(thr_isotope_dist_comp) > 0:
        rmse = mean_squared_error(exp_isotope_dist_comp, thr_isotope_dist_comp, squared=squared)
    else:
        rmse = np.nan
    return rmse


def PoiBin(success_probabilities):
    """
    poisson binomial probability distribution function
    :param success_probabilities:
    :return: probabiliy distribution
    """
    number_trials = success_probabilities.size

    omega = 2 * jnp.pi / (number_trials + 1)

    chi = jnp.empty(number_trials + 1, dtype=complex)
    chi = chi.at[0].set(1)
    # chi[0] = 1
    half_number_trials = int(
        number_trials / 2 + number_trials % 2)
    # set first half of chis:

    # idx_array = np.arange(1, half_number_trials + 1)
    exp_value = jnp.exp(omega * jnp.arange(1, half_number_trials + 1) * 1j)
    xy = 1 - success_probabilities + success_probabilities * exp_value[:, jnp.newaxis]
    # sum over the principal values of the arguments of z:
    argz_sum = jnp.arctan2(xy.imag, xy.real).sum(axis=1)
    # get d value:
    # exparg = np.log(np.abs(xy)).sum(axis=1)
    d_value = jnp.exp(jnp.log(jnp.abs(xy)).sum(axis=1))
    # get chi values:
    chi = chi.at[1:half_number_trials + 1].set(d_value * jnp.exp(argz_sum * 1j))
    # chi[1:half_number_trials + 1] = d_value * jnp.exp(argz_sum * 1j)

    # set second half of chis:

    chi = chi.at[half_number_trials + 1: number_trials + 1].set(jnp.conjugate(
        chi[1:number_trials - half_number_trials + 1][::-1]))

    # chi[half_number_trials + 1:number_trials + 1] = jnp.conjugate(
    #     chi[1:number_trials - half_number_trials + 1][::-1])
    chi /= number_trials + 1
    xi = jnp.fft.fft(chi)
    return xi.real


def theoretical_isotope_dist(sequence, num_isotopes=None):
    """
    calculate theoretical isotope distribtuion from a given one letter sequence of protein chain
    :param sequence: protein sequence in one letter code
    :param num_isotopes: number of isotopes to include. If none, includes all
    :return: isotope distribution
    """
    seq_formula = molmass.Formula(sequence)
    isotope_dist = jnp.array([x[1] for x in seq_formula.spectrum().values()])
    isotope_dist = isotope_dist / jnp.max(isotope_dist)
    if num_isotopes:
        if num_isotopes < len(isotope_dist):
            isotope_dist = isotope_dist[:num_isotopes]
        else:
            fill_arr = jnp.zeros(num_isotopes - len(isotope_dist))
            isotope_dist = jnp.append(isotope_dist, fill_arr)
    return isotope_dist


def calc_hx_prob(timepoint,
                 rate_constant,
                 inv_back_exchange,
                 d2o_purity,
                 d2o_fraction):
    """
    calculate the exchange probability for each residue at the timepoint given the rate constant, backexchange, d2o purity and fraction
    :param timepoint: timepoint in seconds
    :param rate_constant: array of rate_constant
    :param inv_back_exchange: 1 - backexchange value
    :param d2o_purity: d2o purity
    :param d2o_fraction: d2o fraction
    :return: array of probabilities
    """

    prob = (1.0 - jnp.exp(-rate_constant * timepoint)) * (d2o_fraction * d2o_purity * inv_back_exchange)
    return prob


def hx_rates_probability_distribution(timepoint,
                                      rates,
                                      inv_backexchange,
                                      d2o_fraction,
                                      d2o_purity,
                                      free_energy_values=None,
                                      temperature=None):
    """
    generate rate of hx probabilities for all residues
    :param timepoint: hdx timepoint in seconds
    :param rates: rates
    :param inv_backexchange: 1 - backexchange
    :param d2o_fraction: d2o fractyion
    :param d2o_purity: d2o purity
    :param free_energy_values: free energy values
    :param temperature: temperature
    :return: hx probablities of all resdiues given the hdx rates
    """

    fractions = jnp.array([1 for x in rates])
    if free_energy_values is not None:
        if temperature is None:
            raise ValueError('You need to specify temperature (K) in order to use free energy values')
        else:
            fractions = jnp.exp(-free_energy_values / (r_constant * temperature)) / (
                        1.0 + jnp.exp(-free_energy_values / (r_constant * temperature)))

    rate_constant_values = rates * fractions

    hx_probabs = calc_hx_prob(timepoint=timepoint,
                              rate_constant=rate_constant_values,
                              inv_back_exchange=inv_backexchange,
                              d2o_purity=d2o_purity,
                              d2o_fraction=d2o_fraction)

    return hx_probabs


def isotope_dist_from_PoiBin(sequence,
                             timepoint,
                             inv_backexchange,
                             rates: np.ndarray,
                             d2o_fraction,
                             d2o_purity,
                             num_bins,
                             free_energy_values=None,
                             temperature=None):
    """
    generate theoretical isotope distribution based on hdx rates, timepoint, and other conditions
    :param sequence: protein sequence str
    :param timepoint: timepoint in seconds
    :param inv_backexchange: 1 - backexchange
    :param rates: hdx rates
    :param d2o_fraction: d2o fraction
    :param d2o_purity: d2o purity
    :param num_bins: number of bins to include for isotope distribution
    :param free_energy_values: free energy values to be used for calculating isotope distribution
    :param temperature: temperature to be used when free energy is not None
    :return: isotope distribution normalized
    """

    hx_probs = hx_rates_probability_distribution(timepoint=timepoint,
                                                 rates=rates,
                                                 inv_backexchange=inv_backexchange,
                                                 d2o_fraction=d2o_fraction,
                                                 d2o_purity=d2o_purity,
                                                 free_energy_values=free_energy_values,
                                                 temperature=temperature)

    pmf_hx_probs = PoiBin(hx_probs)

    seq_isotope_dist = theoretical_isotope_dist(sequence=sequence, num_isotopes=num_bins)

    isotope_dist_poibin = jnp.convolve(pmf_hx_probs, seq_isotope_dist)[:num_bins]
    isotope_dist_poibin_norm = isotope_dist_poibin / jnp.max(isotope_dist_poibin)

    return isotope_dist_poibin_norm


def gen_theoretical_isotope_dist_for_all_timepoints(sequence,
                                                    timepoints,
                                                    rates,
                                                    inv_backexchange_array,
                                                    d2o_fraction,
                                                    d2o_purity,
                                                    num_bins,
                                                    free_energy_values=None,
                                                    temperature=None):
    """

    :param sequence: protein sequence
    :param timepoints: array of hdx timepoints
    :param rates: rates
    :param inv_backexchange_array: inv backexchange array with length equals to length of timepoints
    :param d2o_fraction: d2o fraction
    :param d2o_purity: d2o purity
    :param num_bins: number of bins for isotope distribution
    :param free_energy_values: free energy values array (optional)
    :param temperature: temperature in Kelvin (optional)
    :return: array of theoretical isotope distributions for each timepoint
    """

    out_array = jnp.zeros((len(timepoints), num_bins))

    for ind, (tp, inv_backexch) in enumerate(zip(timepoints, inv_backexchange_array)):
        isotope_dist = isotope_dist_from_PoiBin(sequence=sequence,
                                                timepoint=tp,
                                                inv_backexchange=inv_backexch,
                                                rates=rates,
                                                d2o_fraction=d2o_fraction,
                                                d2o_purity=d2o_purity,
                                                num_bins=num_bins,
                                                free_energy_values=free_energy_values,
                                                temperature=temperature)
        out_array = out_array.at[ind].set(isotope_dist)

    return out_array


def rate_fit_model(num_rates,
                   sequence,
                   timepoints,
                   backexchange_array,
                   d2o_fraction,
                   d2o_purity,
                   num_bins,
                   obs_dist_nonzero_flat,
                   nonzero_indices):
    """
    rate fit model for opt
    :param num_rates: number of rates
    :param sequence: protein sequence
    :param timepoints: timepoints array
    :param inv_backexchange_array:  1- backexchange array
    :param d2o_fraction: d2o fraction
    :param d2o_purity: d2o purity
    :param num_bins: number of bins in integrated mz data
    :param obs_dist_nonzero_flat: observed experimental distribution flattened and with non zero values
    :param nonzero_indices: indices in exp distribution with non zero values
    :return: numpyro sample object
    """

    with numpyro.plate(name='rates', size=num_rates):
        rates_ = numpyro.sample(name='rate',
                                fn=numpyro.distributions.Uniform(low=-15, high=5))

    thr_dists = gen_theoretical_isotope_dist_for_all_timepoints(sequence=sequence,
                                                                timepoints=timepoints,
                                                                rates=jnp.exp(rates_),
                                                                inv_backexchange_array=jnp.subtract(1,
                                                                                                    backexchange_array),
                                                                d2o_fraction=d2o_fraction,
                                                                d2o_purity=d2o_purity,
                                                                num_bins=num_bins)

    flat_thr_dist = jnp.concatenate(thr_dists)
    flat_thr_dist_non_zero = flat_thr_dist[nonzero_indices]

    sigma = numpyro.sample(name='sigma',
                           fn=numpyro.distributions.Normal(loc=0.5, scale=0.5))

    with numpyro.plate(name='bins', size=len(flat_thr_dist_non_zero)):
        return numpyro.sample(name='bin_preds',
                              fn=numpyro.distributions.Normal(loc=flat_thr_dist_non_zero, scale=sigma),
                              obs=obs_dist_nonzero_flat)


# def rate_fit_model_norm_priors_with_backexchange_sampling(num_rates,
#                       sequence,
#                       timepoints,
#                       backexchange_array,
#                       d2o_fraction,
#                       d2o_purity,
#                       num_bins,
#                       obs_dist_nonzero_flat,
#                       nonzero_indices):
#     """
#     rate fit model for opt
#     :param num_rates: number of rates
#     :param sequence: protein sequence
#     :param timepoints: timepoints array
#     :param inv_backexchange_array:  1- backexchange array
#     :param d2o_fraction: d2o fraction
#     :param d2o_purity: d2o purity
#     :param num_bins: number of bins in integrated mz data
#     :param obs_dist_nonzero_flat: observed experimental distribution flattened and with non zero values
#     :param nonzero_indices: indices in exp distribution with non zero values
#     :return: numpyro sample object
#     """
#
#     rate_center = np.linspace(start=-15, stop=5, num=num_rates)
#     rate_sigma = 2.5
#     with numpyro.plate(name='rates', size=num_rates):
#         rates_ = numpyro.sample(name='rate',
#                                 fn=numpyro.distributions.Normal(loc=rate_center, scale=rate_sigma))
#
#     # todo: need to re evaluate the prior distribution params
#     backexchange_sigma = 0.01
#
#     with numpyro.plate(name='backexchange_values', size=len(backexchange_array)):
#         backexchange = numpyro.sample(name='backexchange',
#                                       fn=numpyro.distributions.Normal(loc=backexchange_array,
#                                                                       scale=backexchange_sigma))
#
#     thr_dists = gen_theoretical_isotope_dist_for_all_timepoints(sequence=sequence,
#                                                                 timepoints=timepoints,
#                                                                 rates=jnp.exp(rates_),
#                                                                 inv_backexchange_array=jnp.subtract(1,
#                                                                                                     backexchange),
#                                                                 d2o_fraction=d2o_fraction,
#                                                                 d2o_purity=d2o_purity,
#                                                                 num_bins=num_bins)
#
#     flat_thr_dist = jnp.concatenate(thr_dists)
#     flat_thr_dist_non_zero = flat_thr_dist[nonzero_indices]
#
#     sigma = numpyro.sample(name='sigma',
#                            fn=numpyro.distributions.Normal(loc=0.5, scale=0.5))
#
#     with numpyro.plate(name='bins', size=len(flat_thr_dist_non_zero)):
#         return numpyro.sample(name='bin_preds',
#                               fn=numpyro.distributions.Normal(loc=flat_thr_dist_non_zero, scale=sigma),
#                               obs=obs_dist_nonzero_flat)


def rate_fit_model_norm_priors(num_rates,
                      sequence,
                      timepoints,
                      backexchange_array,
                      d2o_fraction,
                      d2o_purity,
                      num_bins,
                      obs_dist_nonzero_flat,
                      nonzero_indices):
    """
    rate fit model for opt
    :param num_rates: number of rates
    :param sequence: protein sequence
    :param timepoints: timepoints array
    :param inv_backexchange_array:  1- backexchange array
    :param d2o_fraction: d2o fraction
    :param d2o_purity: d2o purity
    :param num_bins: number of bins in integrated mz data
    :param obs_dist_nonzero_flat: observed experimental distribution flattened and with non zero values
    :param nonzero_indices: indices in exp distribution with non zero values
    :return: numpyro sample object
    """

    rate_center = np.linspace(start=-15, stop=5, num=num_rates)
    rate_sigma = 2.5
    with numpyro.plate(name='rates', size=num_rates):
        rates_ = numpyro.sample(name='rate',
                                fn=numpyro.distributions.Normal(loc=rate_center, scale=rate_sigma))

    thr_dists = gen_theoretical_isotope_dist_for_all_timepoints(sequence=sequence,
                                                                timepoints=timepoints,
                                                                rates=jnp.exp(rates_),
                                                                inv_backexchange_array=jnp.subtract(1,
                                                                                                    backexchange_array),
                                                                d2o_fraction=d2o_fraction,
                                                                d2o_purity=d2o_purity,
                                                                num_bins=num_bins)

    flat_thr_dist = jnp.concatenate(thr_dists)
    flat_thr_dist_non_zero = flat_thr_dist[nonzero_indices]

    sigma = numpyro.sample(name='sigma',
                           fn=numpyro.distributions.Normal(loc=0.5, scale=0.5))

    with numpyro.plate(name='bins', size=len(flat_thr_dist_non_zero)):
        return numpyro.sample(name='bin_preds',
                              fn=numpyro.distributions.Normal(loc=flat_thr_dist_non_zero, scale=sigma),
                              obs=obs_dist_nonzero_flat)


def recalc_timepoints_with_merge(timepoints_array_list, merge_facs_):

    multiply_tps_list = [timepoints_array_list[0]]

    for ind, merge_facs_sample in enumerate(merge_facs_):
        mult_tparr = jnp.multiply(10**merge_facs_sample, timepoints_array_list[ind+1])
        multiply_tps_list.append(mult_tparr)

    concat_tp_array = jnp.concatenate(multiply_tps_list)

    return concat_tp_array


def rate_fit_model_norm_priors_with_merge(num_rates,
                                          sequence,
                                          timepoints_array_list,
                                          num_merge_facs,
                                          backexchange_array,
                                          d2o_fraction,
                                          d2o_purity,
                                          num_bins,
                                          obs_dist_nonzero_flat,
                                          nonzero_indices):
    """
    rate fit model for opt
    :param num_rates: number of rates
    :param sequence: protein sequence
    :param timepoints: timepoints array
    :param inv_backexchange_array:  1- backexchange array
    :param d2o_fraction: d2o fraction
    :param d2o_purity: d2o purity
    :param num_bins: number of bins in integrated mz data
    :param obs_dist_nonzero_flat: observed experimental distribution flattened and with non zero values
    :param nonzero_indices: indices in exp distribution with non zero values
    :return: numpyro sample object
    """

    rate_center = np.linspace(start=-15, stop=5, num=num_rates)
    rate_sigma = 2.5
    with numpyro.plate(name='rates', size=num_rates):
        rates_ = numpyro.sample(name='rate',
                                fn=numpyro.distributions.Normal(loc=rate_center, scale=rate_sigma))

    # gen merge priors
    log_merge_prior_sigma = 1.0
    log_merge_prior_center = 3.0

    with numpyro.plate(name='merge_facs', size=num_merge_facs):
        merge_facs_ = numpyro.sample(name='merge_fac', fn=numpyro.distributions.TruncatedNormal(loc=log_merge_prior_center,
                                                                                                scale=log_merge_prior_sigma,
                                                                                                low=0.0))

    concat_tp_arr_ = recalc_timepoints_with_merge(timepoints_array_list=timepoints_array_list,
                                                  merge_facs_=merge_facs_)

    thr_dists = gen_theoretical_isotope_dist_for_all_timepoints(sequence=sequence,
                                                                timepoints=concat_tp_arr_,
                                                                rates=jnp.exp(rates_),
                                                                inv_backexchange_array=jnp.subtract(1,
                                                                                                    backexchange_array),
                                                                d2o_fraction=d2o_fraction,
                                                                d2o_purity=d2o_purity,
                                                                num_bins=num_bins)

    flat_thr_dist = jnp.concatenate(thr_dists)
    flat_thr_dist_non_zero = flat_thr_dist[nonzero_indices]

    sigma = numpyro.sample(name='sigma',
                           fn=numpyro.distributions.Normal(loc=0.5, scale=0.5))

    with numpyro.plate(name='bins', size=len(flat_thr_dist_non_zero)):
        return numpyro.sample(name='bin_preds',
                              fn=numpyro.distributions.Normal(loc=flat_thr_dist_non_zero, scale=sigma),
                              obs=obs_dist_nonzero_flat)

# def rate_fit_model_v2(num_rates,
#                       sequence,
#                       timepoints,
#                       backexchange_array,
#                       d2o_fraction,
#                       d2o_purity,
#                       num_bins,
#                       obs_dist_nonzero_flat,
#                       nonzero_indices):
#     """
#     rate fit model for opt
#     :param num_rates: number of rates
#     :param sequence: protein sequence
#     :param timepoints: timepoints array
#     :param inv_backexchange_array:  1- backexchange array
#     :param d2o_fraction: d2o fraction
#     :param d2o_purity: d2o purity
#     :param num_bins: number of bins in integrated mz data
#     :param obs_dist_nonzero_flat: observed experimental distribution flattened and with non zero values
#     :param nonzero_indices: indices in exp distribution with non zero values
#     :return: numpyro sample object
#     """
#
#     with numpyro.plate(name='rates', size=num_rates):
#         rates_ = numpyro.sample(name='rate',
#                                 fn=numpyro.distributions.Uniform(low=-15, high=5))
#
#     # todo: need to re evaluate the prior distribution params
#     backexchange_sigma = numpyro.sample(name='backexchange_sigma',
#                                         fn=numpyro.distributions.HalfNormal(scale=0.001))
#
#     with numpyro.plate(name='backexchange_values', size=len(backexchange_array)):
#         backexchange = numpyro.sample(name='backexchange',
#                                       fn=numpyro.distributions.Normal(loc=backexchange_array,
#                                                                       scale=backexchange_sigma))
#
#     thr_dists = gen_theoretical_isotope_dist_for_all_timepoints(sequence=sequence,
#                                                                 timepoints=timepoints,
#                                                                 rates=jnp.exp(rates_),
#                                                                 inv_backexchange_array=jnp.subtract(1, backexchange),
#                                                                 d2o_fraction=d2o_fraction,
#                                                                 d2o_purity=d2o_purity,
#                                                                 num_bins=num_bins)
#
#     flat_thr_dist = jnp.concatenate(thr_dists)
#     flat_thr_dist_non_zero = flat_thr_dist[nonzero_indices]
#
#     sigma = numpyro.sample(name='sigma',
#                            fn=numpyro.distributions.Normal(loc=0.5, scale=0.5))
#
#     with numpyro.plate(name='bins', size=len(flat_thr_dist_non_zero)):
#         return numpyro.sample(name='bin_preds',
#                               fn=numpyro.distributions.Normal(loc=flat_thr_dist_non_zero, scale=sigma),
#                               obs=obs_dist_nonzero_flat)


def sort_posterior_rates_in_samples(posterior_samples_by_chain):
    """
    sort the posterior rates in each chain
    :param posterior_samples_by_chain: posterior samples by chain provided by mcmc
    :return: posterior samples by chain
    """

    for key, values in posterior_samples_by_chain.items():

        if key == 'rate':

            store_chain_values = []

            for ind, chain_values in enumerate(values):
                new_chain_values = jnp.sort(chain_values, axis=1)
                store_chain_values.append(new_chain_values)

            store_chain_array = jnp.array(store_chain_values)

            posterior_samples_by_chain['rate'] = store_chain_array

    return posterior_samples_by_chain


if __name__ == '__main__':
    pass

    # import numpy as np
    # from methods import normalize_mass_distribution_array, gauss_fit_to_isotope_dist_array, plot_hx_rate_fitting_bayes
    # from hx_rate_fit import calc_back_exchange
    # from hxdata import load_tp_dependent_dict, load_data_from_hdx_ms_dist_
    #
    #
    #
    # temp = 295.0
    #
    # d2o_frac = 0.95
    # d2o_pur = 0.95
    #
    # low_ph = 6.0
    # high_ph = 9.0
    #
    # sequence = 'HMVAVPQLIGSTVKEARAKAEKAGLKIDAGDAKSNDRVLVQNPLPGFSAERDSVITVKTV'
    #
    # low_ph_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/dG_May2023/debug_ratefits/Lib08/hdxlim/ph6/C4LI33.1_661-718_9.77_winner_multibody.cpickle.zlib.csv'
    # high_ph_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/dG_May2023/debug_ratefits/Lib08/hdxlim/ph9/C4LI33.1_661-718_9.54_winner_multibody.cpickle.zlib.csv'
    #
    # prot_name = 'C4LI33.1_661-718'
    # prot_rt_name_ = 'C4LI33.1_661-718_9.77_C4LI33.1_661-718_9.54'
    #
    # low_ph_bkexch_corr_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/dG_May2023/debug_ratefits/Lib08/backexchange/low_ph_bkexch_corr.csv'
    # high_ph_bkexch_corr_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/dG_May2023/debug_ratefits/Lib08/backexchange/high_ph_bkexch_corr.csv'
    #
    # low_ph_bkexch_corr_dict = load_tp_dependent_dict(filepath=low_ph_bkexch_corr_fpath)
    # high_ph_bkexch_corr_dict = load_tp_dependent_dict(filepath=high_ph_bkexch_corr_fpath)
    #
    # lowph_hxmsdata_dict = load_data_from_hdx_ms_dist_(fpath=low_ph_fpath)
    # low_ph_timepoints = lowph_hxmsdata_dict['tp']
    # low_ph_ms_dist = lowph_hxmsdata_dict['mass_dist']
    #
    # highph_hxmsdata_dict = load_data_from_hdx_ms_dist_(fpath=high_ph_fpath)
    # high_ph_timepoints = highph_hxmsdata_dict['tp']
    # high_ph_ms_dist = highph_hxmsdata_dict['mass_dist']
    #
    # low_ph_ms_norm_dist = normalize_mass_distribution_array(mass_dist_array=low_ph_ms_dist)
    # high_ph_ms_norm_dist = normalize_mass_distribution_array(mass_dist_array=high_ph_ms_dist)
    #
    # # try for single ph rate fitting
    #
    # backexchange_obj_lowph = calc_back_exchange(sequence=sequence,
    #                                             experimental_isotope_dist=low_ph_ms_norm_dist[-1],
    #                                             timepoints_array=low_ph_timepoints,
    #                                             d2o_fraction=d2o_frac,
    #                                             d2o_purity=d2o_pur,
    #                                             backexchange_corr_dict=low_ph_bkexch_corr_dict)
    #
    # backexchange_obj_highph = calc_back_exchange(sequence=sequence,
    #                                              experimental_isotope_dist=high_ph_ms_norm_dist[-1],
    #                                              timepoints_array=high_ph_timepoints,
    #                                              d2o_fraction=d2o_frac,
    #                                              d2o_purity=d2o_pur,
    #                                              backexchange_corr_dict=high_ph_bkexch_corr_dict)
    #
    # expdata_obj = ExpDataRateFit(sequence=sequence,
    #                              prot_name=prot_name,
    #                              prot_rt_name=prot_rt_name_,
    #                              timepoints=[low_ph_timepoints, high_ph_timepoints],
    #                              timepoint_index_list=None,
    #                              exp_distribution=[low_ph_ms_norm_dist, high_ph_ms_norm_dist],
    #                              exp_label=['ph6', 'ph9'],
    #                              backexchange=[backexchange_obj_lowph.backexchange_array, backexchange_obj_highph.backexchange_array],
    #                              merge_exp=True,
    #                              d2o_purity=d2o_pur,
    #                              d2o_fraction=d2o_frac)
    #
    # print('heho')
    #
    # bayesopt = BayesRateFit(num_chains=4,
    #                         num_warmups=5,
    #                         num_samples=5,
    #                         sample_backexchange=False)
    #
    # bayesopt.fit_rate(exp_data_object=expdata_obj)
