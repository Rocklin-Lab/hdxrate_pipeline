import pandas as pd
import numpy as np
import argparse
import json


def max_value(rows):
    l = [max(dictionary.values()) for dictionary in rows[["idotp", "idotp_pH6", "idotp_pH9"]] if
         dictionary is not np.nan]
    return np.min(l)

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
    # Define file paths based on prefix
    prefix = args.prefix if args.prefix else ""
    single_pH_file = f"{prefix}PO_and_RateFittingResults_nomatches.json"
    two_pHs_file = f"{prefix}PO_and_RateFittingResults_matches.json"
    output_unfiltered = f"{prefix}df_unfiltered_data.json"
    output_filtered = f"{prefix}df_filtered.json"
    output_deduplicated = f"{prefix}df_deduplicated.json"

    # Load data
    df_unmerged = pd.read_json(single_pH_file) if args.single_pH_file else pd.DataFrame()
    df_merged = pd.read_json(two_pHs_file) if args.two_pHs_file else pd.DataFrame()

    if not df_merged.empty:
        df_merged["max_idotp"] = df_merged.apply(lambda x: max_value(x), axis=1)
    if not df_unmerged.empty:
        df_unmerged["max_idotp"] = df_unmerged.apply(lambda x: max_value(x), axis=1)

    items = ['PO_baseline_peak_error', 'PO_dt_ground_rmse', 'PO_rt_ground_rmse', 'PO_dt_ground_fit', 'PO_rt_ground_fit',
             'PO_rmses_sum']
    
    for item in items:
        if not df_merged.empty:
            df_merged[item] = (df_merged[item + "_pH6"] + df_merged[item + "_pH9"]) / 2

    if not df_merged.empty:
        df_merged["PO_total_score"] = df_merged["PO_total_score"] / 2

    # Concatenate dataframes
    df_unfiltered = pd.concat([df_unmerged, df_merged]).reset_index(drop=True)

    df_unfiltered.to_json(args.output_unfiltered)

    # Define filters
    # Define filters

    base_filters = "PO_baseline_peak_error < 4 & max_idotp > 0.99" # TODO: This needs to be addressed later --> max_qvalue < 0.024
    filters_unmerged_too_fast = "saturation_abs_diff_pH6 < 2 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & delta_mass_avg < 0.07 & (slowest_rate > -4.5 | ratio_n_measurable_obs_rates_n_exch <= 0.20 | dg_mean < 2)"
    filters_unmerged = "saturation_abs_diff_pH6 < 2 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & n_tp_stable_res > 4 & rate_fitting_rmse_most_stable < 0.20 & rate_fitting_rmse_before_full_deuteration < 0.2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & delta_mass_avg < 0.07 & ratio_n_measurable_obs_rates_n_exch > 0.20 & slowest_rate < -4.5 & slowest_rate > -11.5 & dg_mean > 2"
    filters_merged_measurable = "saturation_abs_diff_pH9 < 2 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & slowest_rate < -4.5 & n_tp_stable_res > 4 & rate_fitting_rmse_most_stable < 0.20 & rate_fitting_rmse_before_full_deuteration < 0.2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & delta_mass_avg < 0.07 & tp_overlap_count > 1 & overlap_mass_diff < 2 & ratio_n_measurable_obs_rates_n_exch > 0.20 & dg_mean > 2"
    filters_merged_not_fully_deuterated = "saturation_abs_diff_pH9 > 2 & tp_overlap_count > 1 & overlap_mass_diff < 2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & rate_fitting_rmse_percentile_90 < 0.3 & delta_mass_max_abs < 0.07 & slowest_rate < -17 & ratio_n_measurable_obs_rates_n_exch > 0.20"
    

    # Filter and print results
    print("Too fast unmerged")
    set_0 = set(df_unmerged.query(f"{base_filters} & {filters_unmerged_too_fast}").sequence)
    print(len(df_unmerged.query(f"{base_filters}  & {filters_unmerged_too_fast}")), len(set_0))

    print("Measurable unmerged")
    set_1 = set(df_unmerged.query(f"{base_filters} & {filters_unmerged}").sequence)
    print(len(df_unmerged.query(f"{base_filters} & {filters_unmerged}")), len(set_1))

    print("Measurable merged")
    set_2 = set(df_merged.query(f"{base_filters} & {filters_merged_measurable}").sequence)
    print(len(df_merged.query(f"{base_filters} & {filters_merged_measurable}")), len(set_2))

    print("Not fully deuterated merged")
    set_3 = set(df_merged.query(f"{base_filters} & {filters_merged_not_fully_deuterated}").sequence)
    print(len(df_merged.query(f"{base_filters} & {filters_merged_not_fully_deuterated}")), len(set_3))

    print("Total unique sequences")
    print(len(set_0.union(set_1).union(set_2).union(set_3)))

    df_0 = df_unmerged.query(f"{base_filters} & {filters_unmerged_too_fast}").reset_index(drop=True).copy()
    df_0["group"] = "group_0: unmeasurable unmerged"

    df_1 = df_unmerged.query(f"{base_filters}  & {filters_unmerged}").reset_index(drop=True).copy()
    df_1["group"] = "group_1: measurable unmerged"

    df_2 = df_merged.query(f"{base_filters} & {filters_merged_measurable}").reset_index(drop=True).copy()
    df_2["group"] = "group_2: measurable merged"

    df_3 = df_merged.query(f"{base_filters} & {filters_merged_not_fully_deuterated}").reset_index(drop=True).copy()
    df_3["group"] = "group_3: not fully deuterated merged"

    df_filtered = pd.concat([df_0, df_1, df_2, df_3]).reset_index(drop=True)
    print(f"Total rows in filtered data: {len(df_filtered)} \n # of unique sequences: {len(set(df_filtered.sequence))} ")
    df_filtered.to_json(args.output_filtered)
    
    df_deduplicated = df_filtered.sort_values('PO_total_score').drop_duplicates('sequence').reset_index(drop=True)
    df_deduplicated.to_json(args.output_deduplicated)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process JSON files and filter data.')
    parser.add_argument('-s', '--single_pH_file', type=str, default='PO_and_RateFittingResults_nomatches.json',
                        help='Path to the merged JSON file')
    parser.add_argument('-t', '--two_pHs_file', type=str, default='PO_and_RateFittingResults_matches.json',
                        help='Path to the unmerged JSON file')
    parser.add_argument('-ou', '--output_unfiltered', type=str, default='df_unfiltered_data.json',
                        help='Output path for unfiltered data')
    parser.add_argument('-of','--output_filtered', type=str, default='df_filtered.json', help='Output path for filtered data')
    parser.add_argument('-od', '--output_deduplicated', type=str, default='df_deduplicated.json', help='Output path for filtered data')
    parser.add_argument('-p', '--prefix', type=str, default='', help='Prefix to preserve mapping between input and output files)

    args = parser.parse_args()
    main(args)
