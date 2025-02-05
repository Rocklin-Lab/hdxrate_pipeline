import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error as mse
import pickle
from scipy.signal import savgol_filter
from scipy.stats import percentileofscore

import pdb


def process_ex1_metrics(df_original):
    # Prepare initial dataframe
    df = df_original[["name", "library", "sequence", "name_rt-group_pH6", "name_rt-group_pH9", "file_path"]].copy()
    out = []

    # Load and process pickle files
    for i, line in df.iterrows():
        fp = line["file_path"]
        row = line.tolist()
        with open(fp, 'rb') as file:
            fe_pickle = pickle.load(file)
        row.extend([
            fe_pickle['hxrate_output']['exp_distribution'],
            fe_pickle['hxrate_output']['pred_distribution'],
            [x['width'] for x in fe_pickle['hxrate_output']['pred_dist_guass_fit']],
            [x['width'] for x in fe_pickle['hxrate_output']['exp_dist_gauss_fit']],
            [x['rmse'] for x in fe_pickle['hxrate_output']['exp_dist_gauss_fit']],
            [x['centroid'] for x in fe_pickle['hxrate_output']['pred_dist_guass_fit']],
            fe_pickle['hxrate_output']['backexchange'],
#            fe_pickle['hxrate_output']['timepoints'],
            fe_pickle['hxrate_output']['tp_ind_label'],
#            fe_pickle['hxrate_output']['timepoints_sort_indices'],
            fe_pickle['hxrate_output']['merge_exp'],
        ])
        out.append(row)

    # Build intermediate dataframe
    rdf = pd.DataFrame(out, columns=[
        'name', 'library', 'sequence', 'name_rt-group_pH6', 'name_rt-group_pH9', 'file_path',
        'exp_distribution', 'pred_distribution', 'fit_width', 'exp_width', 'exp_rmse',
        'fit_center', 'backexchange', # 'timepoints',
        'tp_ind_label',
        #'timepoints_sort_indices',
        'merge_exp'
    ])

    # Add computed metrics
    rdf['savgol'] = [savgol_filter(expdata, window_length=5, polyorder=3) for expdata in rdf['exp_width']]
    rdf['backexchange'] = [np.array([0] + list(x[1:])) for x in rdf['backexchange']]
    rdf['fit_center'] = [np.array(x) for x in rdf['fit_center']]
    rdf['true_center'] = [fit_center / (1 - backexchange) for fit_center, backexchange in zip(rdf['fit_center'], rdf['backexchange'])]
    rdf['fit_ratio'] = [np.array(savgol) / np.array(fit) for fit, savgol in zip(rdf['fit_width'], rdf['savgol'])]
    rdf['fit_fraction_end'] = [true_center[np.argmax(fit_ratio)] / true_center[-1] for true_center, fit_ratio in zip(rdf['true_center'], rdf['fit_ratio'])]
    rdf['max_ratio'] = [max(x) for x in rdf['fit_ratio']]

    # Merge processed data with the original dataframe
    df = pd.merge(df_original, rdf, on=['name', 'library', 'sequence', 'name_rt-group_pH6', 'name_rt-group_pH9', 'file_path'], how="left")

    #pdb.set_trace()

    # Define helper function for true time calculation
    def true_time(timepoints, tp_ind_label, merge_factor_mean):
        if pd.isna(merge_factor_mean):
            return timepoints
        scale = np.array(['pH9' in x for x in tp_ind_label])
        scalefac = np.ones(len(timepoints))
        scalefac[scale] = 10 ** merge_factor_mean
        return timepoints / scalefac

    # Add more computed columns
    df['true_time'] = [true_time(tp, lbl, mf) for tp, lbl, mf in zip(df['timepoints'], df['tp_ind_label'], df['merge_exp'])]
    df['log_true_time'] = [[np.log10(a) if a != 0 else 0 for a in x] for x in df['true_time']]

    df['pH6']=[np.array(['pH6' in a for a in x]) for x in df['tp_ind_label'].values]

    df['avg_fit_ratio_pH6'] = [np.average(x[pH6]) for x, pH6 in zip(df['fit_ratio'], df['pH6'])]
    df['median_fit_ratio_pH6'] = [np.median(x[pH6]) for x, pH6 in zip(df['fit_ratio'], df['pH6'])]
    df['max_fit_ratio_pH6'] = [max(x[pH6]) for x, pH6 in zip(df['fit_ratio'], df['pH6'])]
    df['max_over_median_ratio_pH6'] = [max(x[pH6]) / np.median(x[pH6]) for x, pH6 in zip(df['fit_ratio'], df['pH6'])]

    df['avg_fit_ratio_pH9'] = [np.average(x[~pH6]) for x, pH6 in zip(df['fit_ratio'], df['pH6'])]
    df['median_fit_ratio_pH9'] = [np.median(x[~pH6]) for x, pH6 in zip(df['fit_ratio'], df['pH6'])]
    df['max_fit_ratio_pH9'] = [max(x[~pH6], default=np.nan) for x, pH6 in zip(df['fit_ratio'], df['pH6'])]
    df['max_over_median_ratio_pH9'] = [max(x[~pH6], default=np.nan) / np.median(x[~pH6]) for x, pH6 in zip(df['fit_ratio'], df['pH6'])]

    df['progress'] = [(fit_center / (1.0 - backexchange)) / (fit_center[-1] / (1.0 - backexchange[-1])) for fit_center, backexchange in zip(df['fit_center'], df['backexchange'])]
    df['progress_pH6'] = [progress[pH6] for progress, pH6 in zip(df['progress'], df['pH6'])]
    df['progress_pH9'] = [progress[~pH6] for progress, pH6 in zip(df['progress'], df['pH6'])]

    df['ratio_window_pH6'] = [np.average(fit_ratio[pH6][(progress_pH6 > 0.80) & (progress_pH6 < 0.90)]) for fit_ratio, pH6, progress_pH6 in zip(df['fit_ratio'], df['pH6'], df['progress_pH6'])]
    df['ratio_window_pH9'] = [np.average(fit_ratio[~pH6][(progress_pH9 > 0.80) & (progress_pH9 < 0.90)]) for fit_ratio, pH6, progress_pH9 in zip(df['fit_ratio'], df['pH6'], df['progress_pH9'])]

    df['progress_window2'] = [(progress > 0.65) & (progress < 0.95) for progress in df['progress']]

    # EXP VS SMOOTH EXPERIMENTAL DISTRIBUTION

    df['obs_max_pH6'] = [np.argmax(np.array(exp_distribution)[pH6 & progress_window2],axis=1) for exp_distribution, pH6, progress_window2 in zip(df['exp_distribution'], df['pH6'], df['progress_window2'])]
    df['obs_max_pH9'] = [np.argmax(np.array(exp_distribution)[~pH6 & progress_window2],axis=1) for exp_distribution, pH6, progress_window2 in zip(df['exp_distribution'], df['pH6'], df['progress_window2'])]
    df['fit_center_pH6'] = [fit_center[pH6 & progress_window2] for fit_center, pH6, progress_window2 in zip(df['fit_center'], df['pH6'], df['progress_window2'])]
    df['fit_center_pH9'] = [fit_center[~pH6 & progress_window2] for fit_center, pH6, progress_window2 in zip(df['fit_center'], df['pH6'], df['progress_window2'])]

    df['fit_minus_obs_max_pH6'] = [np.max(fit_center - obs_max, initial=0) for fit_center, obs_max in zip(df['fit_center_pH6'], df['obs_max_pH6'])]
    df['fit_minus_obs_min_pH6'] = [np.min(fit_center - obs_max, initial=0) for fit_center, obs_max in zip(df['fit_center_pH6'], df['obs_max_pH6'])]
    df['fit_minus_obs_range_pH6'] = [np.max(fit_center - obs_max, initial=0) - np.min(fit_center - obs_max, initial=0) for fit_center, obs_max in zip(df['fit_center_pH6'], df['obs_max_pH6'])]
    df['fit_minus_obs_max_index_pH6'] = [np.argmax(list(fit_center - obs_max) + [0]) for fit_center, obs_max in zip(df['fit_center_pH6'], df['obs_max_pH6'])]
    df['fit_minus_obs_min_index_pH6'] = [np.argmin(list(fit_center - obs_max) + [0]) for fit_center, obs_max in zip(df['fit_center_pH6'], df['obs_max_pH6'])]
    df['fit_ahead_progress_pH6'] = [np.array(list(progress_pH6) + [0])[fit_minus_obs_max_index_pH6] for progress_pH6, fit_minus_obs_max_index_pH6 in zip(df['progress_pH6'], df['fit_minus_obs_max_index_pH6'])]


    df['fit_minus_obs_max_pH9'] = [np.max(fit_center - obs_max, initial=0) for fit_center, obs_max in zip(df['fit_center_pH9'], df['obs_max_pH9'])]
    df['fit_minus_obs_min_pH9'] = [np.min(fit_center - obs_max, initial=0) for fit_center, obs_max in zip(df['fit_center_pH9'], df['obs_max_pH9'])]
    df['fit_minus_obs_range_pH9'] = [np.max(fit_center - obs_max, initial=0) - np.min(fit_center - obs_max, initial=0) for fit_center, obs_max in zip(df['fit_center_pH9'], df['obs_max_pH9'])]
    df['fit_minus_obs_max_index_pH9'] = [np.argmax(list(fit_center - obs_max) + [0]) for fit_center, obs_max in zip(df['fit_center_pH9'], df['obs_max_pH9'])]
    df['fit_minus_obs_min_index_pH9'] = [np.argmin(list(fit_center - obs_max) + [0]) for fit_center, obs_max in zip(df['fit_center_pH9'], df['obs_max_pH9'])]
    df['fit_ahead_progress_pH9'] = [np.array(list(progress_pH9) + [0])[fit_minus_obs_max_index_pH9] for progress_pH9, fit_minus_obs_max_index_pH9 in zip(df['progress_pH9'], df['fit_minus_obs_max_index_pH9'])]

    df['rate_fitting_rmse_percentile'] = [percentileofscore(df['rate_fitting_rmse_total'].values, x) for x in df['rate_fitting_rmse_total']]

    # Add additional metrics as required for ex1 classification
    ex1_width_criteria_pH6=df.query('(ratio_window_pH6 > 1.25 & max_over_median_ratio_pH6 > 1.35)')
    ex1_width_criteria_pH9=df.query('(ratio_window_pH9 > 1.25 & max_over_median_ratio_pH9 > 1.35)')

    # ex1_width_criteria_pH6=df.query('( max_over_median_ratio_pH6 > 1.4)')
    # ex1_width_criteria_pH9=df.query('( max_over_median_ratio_pH9 > 1.4)')

    ex1_jump_criteria_pH6=df.query('( ( (fit_minus_obs_range_pH6 > 5) | (fit_minus_obs_min_pH6 < -2 & fit_minus_obs_max_pH6 > 2) ) & fit_minus_obs_max_index_pH6 < fit_minus_obs_min_index_pH6 & fit_ahead_progress_pH6 < 0.9)')
    ex1_jump_criteria_pH9=df.query('( ( (fit_minus_obs_range_pH9 > 5) | (fit_minus_obs_min_pH9 < -2 & fit_minus_obs_max_pH9 > 2) ) & fit_minus_obs_max_index_pH9 < fit_minus_obs_min_index_pH9 & fit_ahead_progress_pH9 < 0.9)')

    df['ex1_width_criteria_pH6'] = [((x in ex1_width_criteria_pH6['name_rt-group_pH6'].values) and (y is not None)) for (x, y) in df[['name_rt-group_pH6','pH']].values]
    df['ex1_width_criteria_pH9'] = [((x in ex1_width_criteria_pH9['name_rt-group_pH9'].values) and (y is None)) for (x, y) in df[['name_rt-group_pH9','pH']].values]
    df['ex1_jump_criteria_pH6'] = [((x in ex1_jump_criteria_pH6['name_rt-group_pH6'].values) and (y is not None)) for (x, y) in df[['name_rt-group_pH6','pH']].values]
    df['ex1_jump_criteria_pH9'] = [((x in ex1_jump_criteria_pH9['name_rt-group_pH9'].values) and (y is None)) for (x, y) in df[['name_rt-group_pH9','pH']].values]
    df['ex1_pH6'] = df['ex1_width_criteria_pH6'] | df['ex1_jump_criteria_pH6']
    df['ex1_pH9'] = df['ex1_width_criteria_pH9'] | df['ex1_jump_criteria_pH9']
    df['ex1_either'] = df['ex1_pH6'] | df['ex1_pH9']

    return df


def get_reg_coeff(df, col="centroids", n_last=5):
    df["saturation_ang_coeff"] = df.apply(lambda x: np.polyfit(np.arange(n_last), x[col][-n_last:], 1)[0], axis=1)
    df["saturation_lin_coeff"] = df.apply(lambda x: np.polyfit(np.arange(n_last), x[col][-n_last:], 1)[1], axis=1)
    df["saturation_pred_values"] = df.apply(lambda x: np.polyval(
        [x["saturation_ang_coeff"], x["saturation_lin_coeff"]], np.arange(n_last)
    ), axis=1)
    df["saturation_mse"] = df.apply(lambda x: mse(
        x[col][-n_last:],
        x["saturation_pred_values"]
    ), axis=1)
    df["saturation_abs_diff"] = df.apply(
        lambda x: np.abs(x["saturation_pred_values"][-1] - x["saturation_pred_values"][-n_last]), axis=1)

def get_concatenated_centroids(df):
    df["centroids"] = df.apply(lambda x: x["centroids_pH6"] + x["centroids_pH9"][1:], axis=1)

def get_corrected_centroids(df):
#    pdb.set_trace()
    df["centroids_corr"] = df.apply(lambda x: [x["centroids"][0]] + list(np.array(x["centroids"][1:]) / (1 - np.array(x["backexchange_per_tp"][1:]))), axis=1)

def get_ph_indexes(df):
    df["pH_index"] = df.apply(lambda x: ["pH6" if "pH6" in i else "pH9" for i in x["backexchange_ind_label"]], axis=1)

def get_bx_range(df):
    df["backexchange_range"] = df.apply(
        lambda x: np.max(x["backexchange_per_tp"][1:]) - np.min(x["backexchange_per_tp"][1:]), axis=1)

def remove_data_low_masses(df, centroid_threshold=10):
    return df[~(df.apply(lambda x: x["centroids"][-1] < centroid_threshold, axis=1))]

def get_highest_rate_rmse(df):
    df["rate_fitting_rmse_max"] = df.apply(lambda x: np.max(x["rate_fitting_rmse_per_timepoint"]), axis=1)

def get_overlap_count(df):
    df["tp_overlap_count"] = df.apply(lambda x: np.sum(
        np.array(x["pH_index"])[x["timepoints_sort_indices"]][1:] != np.array(x["pH_index"])[
                                                                         x["timepoints_sort_indices"]][:-1]), axis=1)

def find_neighbor_indexes(lst):
    ph6_index = [i for (i, j) in enumerate(lst) if 'pH6' in j][-1]
    ph9_index = [i for (i, j) in enumerate(lst) if 'pH9' in j][0]
    return ph6_index, ph9_index

def find_mass_diff(row):
    index_1, index_2 = find_neighbor_indexes(row["pH_index"])
    x1, x2 = np.array(row["centroids_corr"])[[index_1, index_2]]
    return x2 - x1

def get_ph6_ph9_mass_diff(df):
    df["overlap_mass_diff"] = df.apply(lambda x: find_mass_diff(x), axis=1)

def fill_with_important_data_merged(df):
    print(f"Getting concatenated centroids - new column: centroids")
    get_concatenated_centroids(df)
    print(f"Getting backexchange corrected centroids - new column: centroids_corr")
#    pdb.set_trace()
    get_corrected_centroids(df)
    print(f"Getting linear reg coeff based on dM/dt for the last 5 timepoints - new columns: saturation_ang_coeff, saturation_lin_coeff, and saturation_mse")
    get_reg_coeff(df, col="centroids_corr")
    print(f"Getting indexes for centroids based on pH6 or pH9 data - new column: pH_index")
    get_ph_indexes(df)
    print(f"Getting range of observed backexchanges - before correction? - new column: backexchange_range")
    get_bx_range(df)
    print(f"Getting max rate rmse - new columns: rate_fitting_rmse_max")
    get_highest_rate_rmse(df)
    print(f"Getting pH6/pH9 overlap count - new column: tp_overlap_count")
    get_overlap_count(df)
    print(f"Getting mass diff between last pH6 tp and first pH9 tp - new column: overlap_mass_diff")
    get_ph6_ph9_mass_diff(df)
    print(f"Removing unreasonable low masses based on last timepoint")
    df = remove_data_low_masses(df)
    print(f"Reseting indexes")
    df.reset_index(drop=True, inplace=True)
    return df

def fill_with_important_data_unmerged(df):

    print(f"Getting backexchange corrected centroids - new column: centroids_corr")
    get_corrected_centroids(df)

    print(
        f"Getting linear reg coeff based on dM/dt for the last 5 timepoints - new columns: saturation_ang_coeff, saturation_lin_coeff, and saturation_mse")
    get_reg_coeff(df, col="centroids_corr")

    print(f"Getting indexes for centroids based on pH6 or pH9 data - new column: pH_index")
    get_ph_indexes(df)

    print(f"Getting range of observed backexchanges - before correction? - new column: backexchange_range")
    get_bx_range(df)

    print(f"Getting max rate rmse - new columns: rate_fitting_rmse_max")
    get_highest_rate_rmse(df)

    print(f"Removing unreasonable low masses based on last timepoint")
    df = remove_data_low_masses(df)

    print(f"Reseting indexes")
    df.reset_index(drop=True, inplace=True)

    return df

def get_tp_boundaries(k, timepoints=None, lower_bound=0.15, upper_bound=0.85, ti=25, exp_factor=1.325, n_timepoints=30, merge_factor=0):
    if timepoints is None:
        timepoints = 10 ** merge_factor * (ti * exp_factor ** np.arange(n_timepoints))
    m = 1 - np.exp(-np.exp(k) * timepoints)
    mask = (m >= lower_bound) & (m <= upper_bound)
    return timepoints[mask]

def get_mean_rmse_most_stable(x, n_res=5):
    k_mean = np.mean(x["rates_mean"][:n_res])
    timepoints = np.array(x["timepoints"])
    tp_boundaries = get_tp_boundaries(k_mean, timepoints=timepoints, lower_bound=0.1, upper_bound=0.9)
    if len(tp_boundaries) < 2:
        return 0.5
    else:
        tp_min, tp_max = tp_boundaries[0], tp_boundaries[-1]
    mask = ((timepoints >= tp_min) & (timepoints <= tp_max))
    rmse_mean = np.mean(np.array(x["rate_fitting_rmse_per_timepoint"])[mask])
    if np.isnan(rmse_mean):
        return 0.5
    else:
        return rmse_mean

def get_mean_rmse_before_full_deuteration(x, merge_factor=0):
    k_mean = np.mean(x["rates_mean"][:5])
    timepoints = np.array(x["timepoints"])
    tp_boundaries = get_tp_boundaries(k_mean, timepoints=timepoints)
    if len(tp_boundaries) < 2:
        return 0.5
    else:
        tp_min, tp_max = tp_boundaries[0], tp_boundaries[-1]
    mask = (timepoints <= tp_max)
    rmse_mean = np.mean(np.array(x["rate_fitting_rmse_per_timepoint"])[mask])
    if np.isnan(rmse_mean):
        return 0.5
    else:
        return rmse_mean

def get_rate_fitting_relative_error_max(x, fastest_rate_threshold=-4):
    mask = np.array(x["rates_mean"]) < fastest_rate_threshold
    errors_subset = np.array(x["rate_fitting_relative_error"])[mask][:5]
    if len(errors_subset) == 0:
        return 0
    else:
        return np.max(errors_subset)

def get_delta_mass(x, lower_threshold=5e2):
    timepoints = np.array(x["timepoints"])[x["timepoints_sort_indices"]]
    centroids = np.array(x["centroids_corr"])[x["timepoints_sort_indices"]]
    mask = timepoints > lower_threshold
    return (centroids[mask][1:] - centroids[mask][:-1]) / np.max(centroids[mask])

def get_delta_mass_avg(x, window=2):
    arg_max_ = np.argmax(np.abs(x["delta_mass"]))
    if arg_max_ < window:
        return np.mean(x["delta_mass"][:arg_max_ + window])
    else:
        return np.mean(x["delta_mass"][arg_max_ - window:arg_max_ + window])

def process_dataframe_unmerged(df):
    df = fill_with_important_data_unmerged(df)
    print(
        f"Filling dataframe with rmse percentiles -  new columns: rate_fitting_rmse_percentile_95, rate_fitting_rmse_percentile_90, rate_fitting_rmse_most_stable, rate_fitting_rmse_before_full_deuteration")
    df["rate_fitting_rmse_percentile_95"] = df.apply(lambda x: np.percentile(x["rate_fitting_rmse_per_timepoint"], 95),
                                                     axis=1)
    df["rate_fitting_rmse_percentile_90"] = df.apply(lambda x: np.percentile(x["rate_fitting_rmse_per_timepoint"], 90),
                                                     axis=1)
    df["rate_fitting_rmse_most_stable"] = df.apply(lambda x: get_mean_rmse_most_stable(x), axis=1)
    df["rate_fitting_rmse_before_full_deuteration"] = df.apply(lambda x: get_mean_rmse_before_full_deuteration(x),
                                                               axis=1)

    print(f"Filling dataframe with slowest_rate")
    df["slowest_rate"] = df.apply(lambda x: x["rates_mean"][0], axis=1)

    print("Get number of experimental timepoints the descrive kinetics of slowest exchangeble residues")
    df["n_tp_stable_res"] = df.apply(
        lambda x: len(get_tp_boundaries(np.array(x["slowest_rate"]), timepoints=np.array(x["timepoints"]))), axis=1)
    print("Get relative error of modelled rates as evaluated by their confidence interval")
    df["rate_fitting_relative_error"] = df.apply(
        lambda x: np.abs((np.array(x["rates_ci_95"]) - np.array(x["rates_ci_5"])) / np.array(x["rates_mean"])), axis=1)
    print("Get max relative error for rates slower than e^-4")
    df["rate_fitting_relative_error_max"] = df.apply(lambda x: get_rate_fitting_relative_error_max(x), axis=1)

    print(f"Compute delta mass between adjacent timepoints above 500s")
    df["delta_mass"] = df.apply(lambda x: get_delta_mass(x), axis=1)
    df["delta_mass_max"] = df.apply(lambda x: np.max(x["delta_mass"]), axis=1)
    df["mass_max"] = df.apply(lambda x: np.max(x["centroids_corr"]), axis=1)
    df["delta_mass_max_abs"] = df.apply(lambda x: np.max(np.abs(x["delta_mass"])), axis=1)
    df["delta_mass_avg"] = df.apply(lambda x: get_delta_mass_avg(x, window=2), axis=1)

    print(f"Compute rate different between slowest rate and fifth slowest rate")
    df["rates_mean_gap_3"] = df.apply(lambda x: x["rates_mean"][2] - x["rates_mean"][0], axis=1)
    df["rates_mean_gap_5"] = df.apply(lambda x: x["rates_mean"][4] - x["rates_mean"][0], axis=1)

    df["RT"] = [float(i.split("_")[-1]) for i in df["name_rt-group"]]

    df["PO_total_score_monobody"] = df[
        [i for i in df.keys() if "PO" in i and "delta" not in i and not "monobody" in i][1:-1]].sum(axis=1)


    df["saturation_abs_diff_pH6"] = df.apply(lambda x: np.abs(x["centroids_corr"][-1] - x["centroids_corr"][-5]), axis=1)

    return df


def process_dataframe_merged(df):
    df = fill_with_important_data_merged(df)
    print(f"Filling dataframe with rmse percentiles -  new columns: rate_fitting_rmse_percentile_95, rate_fitting_rmse_percentile_90, rate_fitting_rmse_most_stable, rate_fitting_rmse_before_full_deuteration")
    df["rate_fitting_rmse_percentile_95"] = df.apply(lambda x: np.percentile(x["rate_fitting_rmse_per_timepoint"], 95), axis=1)
    df["rate_fitting_rmse_percentile_90"] = df.apply(lambda x: np.percentile(x["rate_fitting_rmse_per_timepoint"], 90), axis=1)
    df["rate_fitting_rmse_most_stable"] = df.apply(lambda x: get_mean_rmse_most_stable(x), axis=1)
    df["rate_fitting_rmse_before_full_deuteration"] = df.apply(lambda x: get_mean_rmse_before_full_deuteration(x), axis=1)
    print(f"Filling dataframe with slowest_rate")
    df["slowest_rate"] = df.apply(lambda x: x["rates_mean"][0], axis=1)
    print("Get number of experimental timepoints the describe kinetics of slowest exchangeable residues")
    df["n_tp_stable_res"] = df.apply(lambda x: len(get_tp_boundaries(np.array(x["slowest_rate"]), timepoints=np.array(x["timepoints"]))), axis=1)
    print("Get relative error of modelled rates as evaluated by their confidence interval")
    df["rate_fitting_relative_error"] = df.apply(lambda x: np.abs((np.array(x["rates_ci_95"]) - np.array(x["rates_ci_5"])) / np.array(x["rates_mean"])), axis=1)
    print("Get max relative error for rates slower than e^-4")
    df["rate_fitting_relative_error_max"] = df.apply(lambda x: get_rate_fitting_relative_error_max(x), axis=1)
    print(f"Compute delta mass between adjacent timepoints above 500s")
    df["delta_mass"] = df.apply(lambda x: get_delta_mass(x), axis=1)
    df["delta_mass_max"] = df.apply(lambda x: np.max(x["delta_mass"]), axis=1)
    df["mass_max"] = df.apply(lambda x: np.max(x["centroids_corr"]), axis=1)
    df["delta_mass_max_abs"] = df.apply(lambda x: np.max(np.abs(x["delta_mass"])), axis=1)
    df["delta_mass_avg"] = df.apply(lambda x: get_delta_mass_avg(x, window=2), axis=1)
    df["delta_mass_avg_abs"] = np.abs(df["delta_mass_avg"])
    df["PO_total_score"] = df.apply(lambda x: x["PO_total_score_pH6"] + x["PO_total_score_pH9"], axis=1)
    df["merge_factor_ci_range"] = df.apply(lambda x: x["merge_factor_ci_95"] - x["merge_factor_ci_5"], axis=1)
    df["saturation_abs_diff_pH6"] = df.apply(lambda x: np.abs(x["centroids_corr"][len(x["centroids_pH6"])] - x["centroids_corr"][len(x["centroids_pH6"]) - 5]), axis=1)
    print(f"Compute rate different between slowest rate and fifth slowest rate")
    df["rates_mean_gap_3"] = df.apply(lambda x: x["rates_mean"][2] - x["rates_mean"][0], axis=1)
    df["rates_mean_gap_5"] = df.apply(lambda x: x["rates_mean"][4] - x["rates_mean"][0], axis=1)
    df["RT"] = [np.round(np.mean([float(i.split("_")[-1]), float(j.split("_")[-1])]), 2) for (i, j) in zip(df["name_rt-group_pH6"].values, df["name_rt-group_pH9"].values)]
    df["PO_total_score_monobody"] = df[[i for i in df.keys() if "PO" in i and ("delta" not in i) and ("winner" not in i) and ("total" not in i)]].sum(axis=1) / 2
    df['saturation_abs_diff_pH9'] = df.apply(lambda x: np.abs(x['centroids_corr'][-1] - x['centroids_corr'][-5]), axis=1)

    return df
