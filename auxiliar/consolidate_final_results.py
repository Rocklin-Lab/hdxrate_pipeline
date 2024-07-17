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
    # Load data
    df_unmerged = pd.read_json(args.single_pH_file) if args.unmerged_file else pd.DataFrame()
    df_merged = pd.read_json(args.two_pHs_file) if args.merged_file else pd.DataFrame()
    if not df_merged.empty:
        df_merged["max_idotp"] = df_merged.apply(lambda x: max_value(x), axis=1)
    if not df_unmerged.empty:
        df_unmerged["max_idotp"] = df_unmerged.apply(lambda x: max_value(x), axis=1)

#    df_merged['max_idotp'] = 1
#    df_unmerged['max_idotp'] = 1

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
    base_filters = "PO_baseline_peak_error < 4 & max_idotp > 0.99"
    filters_unmerged_too_fast = "saturation_abs_diff < 1 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & slowest_rate > -4.5 & delta_mass_avg < 0.07"
    filters_unmerged = "saturation_abs_diff < 1 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & slowest_rate < -4.5 & slowest_rate > -11.5 & n_tp_stable_res > 4 & rate_fitting_rmse_most_stable < 0.20 & rate_fitting_rmse_before_full_deuteration < 0.2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & delta_mass_avg < 0.07"
    filters_merged_measurable = "saturation_abs_diff < 1 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & slowest_rate < -4.5 & n_tp_stable_res > 4 & rate_fitting_rmse_most_stable < 0.20 & rate_fitting_rmse_before_full_deuteration < 0.2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & delta_mass_avg < 0.07 & tp_overlap_count > 1 & overlap_mass_diff < 2"
    filters_merged_not_fully_deuterated = "saturation_abs_diff > 1 & tp_overlap_count > 1 & overlap_mass_diff < 2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & rate_fitting_rmse_percentile_90 < 0.3 & delta_mass_max_abs < 0.07"

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
    parser.add_argument('-s', '--single_pH_file', type=str, default='consolidated_po_and_rates_unmerged.json',
                        help='Path to the merged JSON file')
    parser.add_argument('-t', '--two_pHs_file', type=str, default='consolidated_po_and_rates_merged.json',
                        help='Path to the unmerged JSON file')
    parser.add_argument('-ou', '--output_unfiltered', type=str, default='df_unfiltered_data.json',
                        help='Output path for unfiltered data')
    parser.add_argument('-of','--output_filtered', type=str, default='df_filtered.json', help='Output path for filtered data')
    parser.add_argument('-od', '--output_deduplicated', type=str, default='df_deduplicated.json', help='Output path for filtered data')

    args = parser.parse_args()
    main(args)
