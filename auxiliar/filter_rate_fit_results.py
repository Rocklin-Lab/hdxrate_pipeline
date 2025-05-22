import pandas as pd
import numpy as np
import argparse
import os
from populate_stats import process_ex1_metrics
import ast

def max_value(row):
    columns = ["idotp", "idotp_pH6", "idotp_pH9"]
    max_values = [max(row[col].values()) for col in columns if col in row and pd.notna(row[col])]
    return np.min(max_values) if max_values else np.nan

def load_data(file_path):
    return pd.read_json(file_path) if os.path.exists(file_path) else pd.DataFrame()

def calculate_metrics(df, ph_condition):
    df["n_exch_res"] = df["free_energy"].apply(len)

    measurable_mask = lambda x: (np.exp(x['rates_ci_95']) - np.exp(x['rates_ci_5'])) / np.abs(np.exp(x['rates_mean'])) < 1

    df['n_measurable_obs_rates'] = df.apply(
        lambda x: np.max(np.where(measurable_mask(x) & (np.array(x['rates_mean']) < -3.45))[0]) + 1
        if np.any(measurable_mask(x) & (np.array(x['rates_mean']) < -3.45)) else 0, axis=1)

    df["pos_max_measurable_rate"] = df.apply(
        lambda x: np.max(np.where(measurable_mask(x))[0]) if np.any(measurable_mask(x)) else 0, axis=1)

    df['n_measurable_intrinsic_rates'] = df.apply(
        lambda x: np.sum((np.array(x['intrinsic_rates']) < 1) & (np.array(x['intrinsic_rates']) > 0)), axis=1)

    df['tp1_back_centroid'] = df['centroids_corr'].apply(
        lambda centroids: (centroids[1] - centroids[0] < 0) or (centroids[2] - centroids[0] < 0))

    df["ratio_n_measurable_obs_rates_n_exch"] = df["n_measurable_obs_rates"] / df["n_exch_res"]
    df['ratio_n_measurable_intrinsic_rates_n_exch'] = df['n_measurable_intrinsic_rates'] / df['n_exch_res']

    return process_ex1_metrics(df, ph_condition)

def apply_filters(df, filters):
    return df.query(filters).reset_index(drop=True)

def extract_column_names_from_query(expr):
    class ColumnVisitor(ast.NodeVisitor):
        def __init__(self):
            self.names = set()

        def visit_Name(self, node):
            self.names.add(node.id)

    try:
        tree = ast.parse(expr, mode='eval')
    except SyntaxError:
        return set()

    visitor = ColumnVisitor()
    visitor.visit(tree)
    return visitor.names

def process_and_label(df, base_filters, specific_filters, label):
    full_filter_expr = f"{base_filters} & {specific_filters}"
    required_cols = extract_column_names_from_query(full_filter_expr)

    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        print(f"Skipping '{label}' due to missing columns: {missing_cols}")
        return pd.DataFrame()

    filtered_df = apply_filters(df, full_filter_expr)
    filtered_df["group"] = label
    return filtered_df


def main(args):
    df_unmerged = load_data(args.single_pH_file)
    if args.two_pHs_file is None:
        df_merged = pd.DataFrame()
    else:
        df_merged = load_data(args.two_pHs_file)

    items = ['PO_baseline_peak_error', 'PO_dt_ground_rmse', 'PO_rt_ground_rmse', 'PO_dt_ground_fit', 'PO_rt_ground_fit', 'PO_rmses_sum']

    if not df_merged.empty:
        df_merged["max_idotp"] = df_merged.apply(max_value, axis=1)
        for item in items:
            df_merged[item] = (df_merged[f"{item}_pH6"] + df_merged[f"{item}_pH9"]) / 2
        df_merged["PO_total_score"] /= 2

    if not df_unmerged.empty:
        df_unmerged["max_idotp"] = df_unmerged.apply(max_value, axis=1)
        for item in items:
            df_unmerged[item] = df_unmerged[f"{item}_pH6"]
        df_unmerged["PO_total_score"] = df_unmerged["PO_total_score_pH6"]

    df_unfiltered = pd.concat([df_unmerged, df_merged], ignore_index=True)

    if args.two_pHs_file is None:
        df_unfiltered = calculate_metrics(df_unfiltered, ph_condition='pH6')
    else:
        df_unfiltered = calculate_metrics(df_unfiltered, ph_condition='both')

    df_unfiltered.to_json(args.output_unfiltered, orient='records', indent=4)

    base_filters = "PO_baseline_peak_error < 4 & max_idotp > 0.99 & ~tp1_back_centroid"

    filter_conditions = {
        "group_0: unmeasurable unmerged": "saturation_abs_diff_pH6 < 2 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & delta_mass_avg < 0.07 & (slowest_rate > -4.5 | ratio_n_measurable_obs_rates_n_exch <= 0.20 | dg_mean < 2)",
        "group_1: measurable unmerged": "saturation_abs_diff_pH6 < 2 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & n_tp_stable_res > 4 & rate_fitting_rmse_most_stable < 0.20 & rate_fitting_rmse_before_full_deuteration < 0.2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & delta_mass_avg < 0.07 & ratio_n_measurable_obs_rates_n_exch > 0.20 & slowest_rate < -4.5 & slowest_rate > -11.5 & dg_mean > 2",
        "group_2: measurable merged": "saturation_abs_diff_pH9 < 2 & rate_fitting_rmse_percentile_90 < 0.3 & backexchange_range < 0.1 & backexchange_final_tp < 0.45 & backexchange_final_tp > 0.1 & slowest_rate < -4.5 & n_tp_stable_res > 4 & rate_fitting_rmse_most_stable < 0.20 & rate_fitting_rmse_before_full_deuteration < 0.2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & delta_mass_avg < 0.07 & tp_overlap_count > 1 & overlap_mass_diff < 2 & ratio_n_measurable_obs_rates_n_exch > 0.20 & dg_mean > 2",
        "group_3: not fully deuterated merged": "saturation_abs_diff_pH9 > 2 & tp_overlap_count > 1 & overlap_mass_diff < 2 & rates_mean_gap_3 < 2 & rates_mean_gap_5 < 3 & rate_fitting_rmse_percentile_90 < 0.3 & delta_mass_max_abs < 0.07 & slowest_rate < -17 & ratio_n_measurable_obs_rates_n_exch > 0.20"
    }

    df_filtered = pd.concat(
        [process_and_label(df_unfiltered, base_filters, condition, label) for label, condition in filter_conditions.items()],
        ignore_index=True
    )

    df_filtered.to_json(args.output_filtered, orient='records', indent=4)

    df_deduplicated = df_filtered.sort_values('PO_total_score').drop_duplicates('sequence').query("~ex1_either").reset_index(drop=True)
    df_deduplicated.to_json(args.output_deduplicated, orient='records', indent=4)

    print(f"Total rows in filtered data: {len(df_filtered)}\nNumber of unique sequences: {df_filtered['sequence'].nunique()}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process JSON files and filter data.')
    parser.add_argument('-s', '--single_pH_file', type=str, required=True, help='Path to the unmerged JSON file')
    parser.add_argument('-t', '--two_pHs_file', type=str, default=None, help='Path to the merged JSON file')
    parser.add_argument('-ou', '--output_unfiltered', type=str, default='df_unfiltered_data.json', help='Output path for unfiltered data')
    parser.add_argument('-of', '--output_filtered', type=str, default='df_filtered.json', help='Output path for filtered data')
    parser.add_argument('-od', '--output_deduplicated', type=str, default='df_deduplicated.json', help='Output path for deduplicated data')

    args = parser.parse_args()
    main(args)
