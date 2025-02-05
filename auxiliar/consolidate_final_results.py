import pandas as pd
import numpy as np
import argparse
import json
import os
import pdb

from populate_stats import process_ex1_metrics

def max_value(rows):
    columns = ["idotp", "idotp_pH6", "idotp_pH9"]
    max_values = []

    for col in columns:
        if col in rows and rows[col] is not np.nan:
            max_values.append(max(rows[col].values()))

    if max_values:
        return np.min(max_values)
    else:
        return np.nan

def main(args):
    # Apply prefix to file paths if prefix is provided
    single_pH_file = args.prefix + "PO_and_RateFittingResults_nomatches.json"
    two_pHs_file = args.prefix + "PO_and_RateFittingResults_matches.json"
    output_unfiltered = args.prefix + 'unfiltered.json'
    output_filtered = args.prefix + 'filtered.json'
    output_deduplicated = args.prefix + 'deduplicated.json'


    # Load data
    if os.path.exists(single_pH_file):
        df_unmerged = pd.read_json(single_pH_file)
    else:
        df_unmerged = pd.DataFrame()

    if os.path.exists(two_pHs_file):
        df_merged = pd.read_json(two_pHs_file)
    else:
        df_merged = pd.DataFrame()


    if not df_merged.empty:
        df_merged["max_idotp"] = df_merged.apply(lambda x: max_value(x), axis=1)
    if not df_unmerged.empty:
        df_unmerged["max_idotp"] = df_unmerged.apply(lambda x: max_value(x), axis=1)

    #pdb.set_trace()

    items = ['PO_baseline_peak_error', 'PO_dt_ground_rmse', 'PO_rt_ground_rmse', 'PO_dt_ground_fit', 'PO_rt_ground_fit', 'PO_rmses_sum']

    for item in items:
        if not df_merged.empty:
            df_merged[item] = (df_merged[item + "_pH6"] + df_merged[item + "_pH9"]) / 2

    if not df_merged.empty:
        df_merged["PO_total_score"] = df_merged["PO_total_score"] / 2

    # Concatenate dataframes
    df_unfiltered = pd.concat([df_unmerged, df_merged]).reset_index(drop=True)

    # Populate with metrics
    df_unfiltered["n_exch_res"] = df_unfiltered["free_energy"].apply(lambda x: len(x))

    df_unfiltered["n_measurable_obs_rates"] = df_unfiltered.apply(
        lambda x: np.max(np.where(
            (np.exp(x['rates_ci_95']) - np.exp(x['rates_ci_5'])) / np.abs(np.exp(x['rates_mean'])) < 1
        )[0]) + 1 if np.any(
            (np.exp(x['rates_ci_95']) - np.exp(x['rates_ci_5'])) / np.abs(np.exp(x['rates_mean'])) < 1
        ) else 0, axis=1
    )


    df_unfiltered["pos_max_measurable_rate"] = df_unfiltered.apply(
        lambda x: np.max(np.where(
            (np.exp(x['rates_ci_95']) - np.exp(x['rates_ci_5'])) / np.abs(np.exp(x['rates_mean'])) < 1
        )[0]) if np.any(
            (np.exp(x['rates_ci_95']) - np.exp(x['rates_ci_5'])) / np.abs(np.exp(x['rates_mean'])) < 1
        ) else 0, axis=1
    )

    df_unfiltered['n_measurable_intrinsic_rates'] = df_unfiltered.apply(
        lambda x: np.sum(
            (np.array(x['intrinsic_rates']) < 1) & (np.array(x['intrinsic_rates']) > 0)
        ), axis=1)



    df_unfiltered['tp1_back_centroid'] = ((np.array([i[1] - i[0] for i in df_unfiltered.centroids_corr]) < 0) | (np.array([i[2] - i[0] for i in df_unfiltered.centroids_corr]) < 0))


    df_unfiltered["ratio_n_measurable_obs_rates_n_exch"] = df_unfiltered["n_measurable_obs_rates"] / df_unfiltered["n_exch_res"]
    df_unfiltered['ratio_n_measurable_intrinsic_rates_n_exch'] = df_unfiltered['n_measurable_intrinsic_rates']/df_unfiltered['n_exch_res']

    df_unfiltered = process_ex1_metrics(df_unfiltered)

    #### finish populating with metrics

    df_unfiltered.to_json(output_unfiltered)

    # Define filters
    base_filters = "PO_baseline_peak_error < 4 & max_idotp > 0.99 & ~tp1_back_centroid" # max_qvalue < 0.024

    filters_unmerged_too_fast = "saturation_abs_diff_pH6 < 2 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & delta_mass_avg < 0.07 & (slowest_rate > -4.5 | ratio_n_measurable_obs_rates_n_exch <= 0.20 | dg_mean < 2)"

    filters_unmerged = "saturation_abs_diff_pH6 < 2 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & n_tp_stable_res > 4 & rate_fitting_rmse_most_stable < 0.20 & rate_fitting_rmse_before_full_deuteration < 0.2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & delta_mass_avg < 0.07 & ratio_n_measurable_obs_rates_n_exch > 0.20 & slowest_rate < -4.5 & slowest_rate > -11.5 & dg_mean > 2"

    filters_merged_measurable = "saturation_abs_diff_pH9 < 2 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & slowest_rate < -4.5 & n_tp_stable_res > 4 & rate_fitting_rmse_most_stable < 0.20 & rate_fitting_rmse_before_full_deuteration < 0.2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & delta_mass_avg < 0.07 & tp_overlap_count > 1 & overlap_mass_diff < 2 & ratio_n_measurable_obs_rates_n_exch > 0.20 & dg_mean > 2"

    filters_merged_not_fully_deuterated = "saturation_abs_diff_pH9 > 2 & tp_overlap_count > 1 & overlap_mass_diff < 2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & rate_fitting_rmse_percentile_90 < 0.3 & delta_mass_max_abs < 0.07 & slowest_rate < -17 & ratio_n_measurable_obs_rates_n_exch > 0.20"


    # Filter and print results
    print("Too fast unmerged")
    set_0 = set(df_unfiltered.query(f"{base_filters} & {filters_unmerged_too_fast}").sequence)
    print(len(df_unfiltered.query(f"{base_filters}  & {filters_unmerged_too_fast}")), len(set_0))

    print("Measurable unmerged")
    set_1 = set(df_unfiltered.query(f"{base_filters} & {filters_unmerged}").sequence)
    print(len(df_unfiltered.query(f"{base_filters} & {filters_unmerged}")), len(set_1))

    print("Measurable merged")
    set_2 = set(df_unfiltered.query(f"{base_filters} & {filters_merged_measurable}").sequence)
    print(len(df_unfiltered.query(f"{base_filters} & {filters_merged_measurable}")), len(set_2))

    print("Not fully deuterated merged")
    set_3 = set(df_unfiltered.query(f"{base_filters} & {filters_merged_not_fully_deuterated}").sequence)
    print(len(df_unfiltered.query(f"{base_filters} & {filters_merged_not_fully_deuterated}")), len(set_3))

    print("Total unique sequences")
    print(len(set_0.union(set_1).union(set_2).union(set_3)))

    df_0 = df_unfiltered.query(f"{base_filters} & {filters_unmerged_too_fast}").reset_index(drop=True).copy()
    df_0["group"] = "group_0: unmeasurable unmerged"

    df_1 = df_unfiltered.query(f"{base_filters}  & {filters_unmerged}").reset_index(drop=True).copy()
    df_1["group"] = "group_1: measurable unmerged"

    df_2 = df_unfiltered.query(f"{base_filters} & {filters_merged_measurable}").reset_index(drop=True).copy()
    df_2["group"] = "group_2: measurable merged"

    df_3 = df_unfiltered.query(f"{base_filters} & {filters_merged_not_fully_deuterated}").reset_index(drop=True).copy()
    df_3["group"] = "group_3: not fully deuterated merged"

    df_filtered = pd.concat([df_0, df_1, df_2, df_3]).reset_index(drop=True)
    print(f"Total rows in filtered data: {len(df_filtered)} \n # of unique sequences: {len(set(df_filtered.sequence))} ")
    df_filtered.to_json(output_filtered)

    df_deduplicated = df_filtered.sort_values('PO_total_score').drop_duplicates('sequence').query("~ex1_either").reset_index(drop=True)
    df_deduplicated.to_json(output_deduplicated)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process JSON files and filter data.')
    parser.add_argument('-p', '--prefix', type=str, default=None, help='Prefix to preserve mapping between input and output files')

    args = parser.parse_args()
    main(args)
