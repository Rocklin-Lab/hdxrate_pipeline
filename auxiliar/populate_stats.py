import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error as mse

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
    return df