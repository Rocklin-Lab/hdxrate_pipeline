import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd

from populate_stats import process_dataframe_unmerged, process_dataframe_merged

import pdb

def replace_nan_with_zeros(array):
    return np.where(np.isnan(array) | (array < 0), 0, array)

def extract_info_from_dg_pickle_file(f, n_highest, library):
    d = pd.read_pickle(f)
    pH = "pH6" if 'pH6' in f else ("pH9" if 'pH9' in f else None)
    protein_name = d["hxrate_output"]["protein_name"]
    protein_rt_name = d["hxrate_output"]["protein_rt_name"]
    sequence = d["hxrate_output"]["sequence"]
    intrinsic_rates, intrinsic_rates_median, reg_intrinsic_rates, netcharge = d["intrinsic_rates"], d[
        "intrinsic_rates_median"], d["reg_intrinsic_rates"], d["netcharge"]
    chain_diag = d["hxrate_output"]["chain_diagnostics"]
    chain_rmse_list = d["hxrate_output"]["chain_diagnostics"]["chain_rmse_list"]
    chain_pass_list = d["hxrate_output"]["chain_diagnostics"]["chain_pass_list"]
    bayes_sample = d["hxrate_output"]["bayes_sample"]
    rates_mean = bayes_sample["rate"]["mean"]
    rates_std = bayes_sample["rate"]["std"]
    rates_ci_5 = bayes_sample["rate"]["ci_5"]
    rates_ci_95 = bayes_sample["rate"]["ci_95"]
    rate_fitting_rmse_total = d["hxrate_output"]["rmse"]["total"]
    rate_fitting_rmse_per_timepoint = d["hxrate_output"]["rmse"]["per_timepoint"]
    backexchange_per_tp = d["hxrate_output"]["backexchange"]
    backexchange_final_tp = backexchange_per_tp[-1]
    backexchange_n_res_subtract = d["hxrate_output"]["back_exchange_res_subtract"]
    backexchange_ind_label = d["hxrate_output"]["tp_ind_label"]
    timepoints = d["hxrate_output"]["timepoints"]
    timepoints_sort_indices = d["hxrate_output"]["timepoints_sort_indices"]
    free_energy_ = d["free_energy"]
    free_energy = replace_nan_with_zeros(free_energy_)
    dg_max = free_energy[0]
    dg_mean = np.mean(free_energy[:n_highest])
    dg_median = np.median(free_energy[:n_highest])
    file_path = os.path.abspath(f)

    return [protein_name, library, pH, protein_rt_name, sequence,
            intrinsic_rates, intrinsic_rates_median, reg_intrinsic_rates, netcharge, chain_rmse_list,
            chain_pass_list,
            rates_mean, rates_std, rates_ci_5, rates_ci_95, rate_fitting_rmse_total,
            rate_fitting_rmse_per_timepoint,
            backexchange_per_tp, backexchange_final_tp, backexchange_n_res_subtract, backexchange_ind_label,
            timepoints, timepoints_sort_indices,
            free_energy_, free_energy, dg_max, dg_mean, dg_median, file_path]

def extract_info_from_dg_pickle_file_with_merge(f, n_highest, library):
    d = pd.read_pickle(f)
    protein_name = d["hxrate_output"]["protein_name"]
    sequence = d["hxrate_output"]["sequence"]
    protein_rt_name = d["hxrate_output"]["protein_rt_name"]
    i_split = len(protein_rt_name.split("_")) // 2
    name_rt_group_pH6 = "_".join(protein_rt_name.split("_")[:i_split])
    name_rt_group_pH9 = "_".join(protein_rt_name.split("_")[i_split:])
    intrinsic_rates, intrinsic_rates_median, reg_intrinsic_rates, netcharge = d["intrinsic_rates"], d[
        "intrinsic_rates_median"], d["reg_intrinsic_rates"], d["netcharge"]
    chain_diag = d["hxrate_output"]["chain_diagnostics"]
    chain_rmse_list = d["hxrate_output"]["chain_diagnostics"]["chain_rmse_list"]
    chain_pass_list = d["hxrate_output"]["chain_diagnostics"]["chain_pass_list"]
    bayes_sample = d["hxrate_output"]["bayes_sample"]
    rates_mean = bayes_sample["rate"]["mean"]
    rates_std = bayes_sample["rate"]["std"]
    rates_ci_5 = bayes_sample["rate"]["ci_5"]
    rates_ci_95 = bayes_sample["rate"]["ci_95"]
    rate_fitting_rmse_total = d["hxrate_output"]["rmse"]["total"]
    rate_fitting_rmse_per_timepoint = d["hxrate_output"]["rmse"]["per_timepoint"]
    merge_factor_mean = bayes_sample["merge_fac"]["mean"][0]
    merge_factor_std = bayes_sample["merge_fac"]["std"][0]
    merge_ci_5 = bayes_sample["merge_fac"]["ci_5"][0]
    merge_ci_95 = bayes_sample["merge_fac"]["ci_95"][0]
    backexchange_per_tp = d["hxrate_output"]["backexchange"]
    backexchange_final_tp = backexchange_per_tp[-1]
    backexchange_n_res_subtract = d["hxrate_output"]["back_exchange_res_subtract"]
    backexchange_ind_label = d["hxrate_output"]["tp_ind_label"]
    timepoints = d["hxrate_output"]["timepoints"]
    timepoints_sort_indices = d["hxrate_output"]["timepoints_sort_indices"]
    free_energy_ = d["free_energy"]
    free_energy = replace_nan_with_zeros(free_energy_)
    dg_max = free_energy[0]
    dg_mean = np.mean(free_energy[:n_highest])
    dg_median = np.median(free_energy[:n_highest])
    file_path = os.path.abspath(f)

    return [protein_name, library, protein_rt_name, name_rt_group_pH6, name_rt_group_pH9, sequence,
            intrinsic_rates, intrinsic_rates_median, reg_intrinsic_rates, netcharge, chain_rmse_list, chain_pass_list,
            rates_mean, rates_std, rates_ci_5, rates_ci_95, rate_fitting_rmse_total, rate_fitting_rmse_per_timepoint,
            merge_factor_mean, merge_factor_std, merge_ci_5, merge_ci_95,
            backexchange_per_tp, backexchange_final_tp, backexchange_n_res_subtract, backexchange_ind_label, timepoints,
            timepoints_sort_indices,
            free_energy_, free_energy, dg_max, dg_mean, dg_median, file_path]

def generate_rates_df(fs, n_highest, library):
    l = []
    for i, f in enumerate(fs):
        if i % 1000 == 0:
            print(f"{i} processed... {len(fs) - i} to go...")
        l.append(extract_info_from_dg_pickle_file(f, n_highest, library))
    return pd.DataFrame(l, columns=["name", "library", "pH", "name_rt-group", "sequence",
                                    "intrinsic_rates", "intrinsic_rates_median", "reg_intrinsic_rates", "netcharge",
                                    "chain_rmse_list", "chain_pass_list",
                                    "rates_mean", "rates_std", "rates_ci_5", "rates_ci_95", "rate_fitting_rmse_total",
                                    "rate_fitting_rmse_per_timepoint",
                                    "backexchange_per_tp", "backexchange_final_tp", "backexchange_n_res_subtract",
                                    "backexchange_ind_label", "timepoints", "timepoints_sort_indices",
                                    "free_energy_", "free_energy", "dg_max", "dg_mean", "dg_median", "file_path"])

def generate_rates_df_with_merge(fs, n_highest, library):
    l = []
    for i, f in enumerate(fs):
        if i % 1000 == 0:
            print(f"{i} processed... {len(fs) - i} to go...")
        l.append(extract_info_from_dg_pickle_file_with_merge(f, n_highest, library))
    return pd.DataFrame(l, columns=["name", "library", "name_rt-group", "name_rt-group_pH6", "name_rt-group_pH9",
                                    "sequence", "intrinsic_rates", "intrinsic_rates_median", "reg_intrinsic_rates",
                                    "netcharge", "chain_rmse_list", "chain_pass_list", "rates_mean", "rates_std",
                                    "rates_ci_5", "rates_ci_95", "rate_fitting_rmse_total",
                                    "rate_fitting_rmse_per_timepoint", "merge_factor_mean", "merge_factor_std",
                                    "merge_factor_ci_5", "merge_factor_ci_95", "backexchange_per_tp",
                                    "backexchange_final_tp", "backexchange_n_res_subtract",
                                    "backexchange_ind_label", "timepoints", "timepoints_sort_indices", "free_energy_",
                                    "free_energy", "dg_max", "dg_mean", "dg_median", "file_path"])

def main():
    parser = argparse.ArgumentParser(description="Process rates files and generate dataframes.")
    parser.add_argument('--n_highest', type=int, default=5, help="Number of highest elements to consider for dg_mean calculation")
    parser.add_argument('--library', type=str, required=True, help="Library identifier")
    parser.add_argument('--merge', action='store_true', help="Include processing with merging information if specified")
    parser.add_argument('--rate_fit_output_folder', type=str, default='rate_fit_output', help="Specify the folder containing rate fitting output")
    parser.add_argument('--po_consolidated_results', type=str, required=True, help="Path to the JSON file with metrics from another pipeline")
    parser.add_argument('--output_nomatches', type=str, required=True, help="Output JSON file name for nomatches", default='consolidated_rates_unmerged.json')
    parser.add_argument('--output_matches', type=str, help="Output JSON file name for matches, if processing with merge", default='consolidated_rates_merged.json')
    parser.add_argument('--output_consolidated_po_and_rates_nomatches', type=str, required=True, help="Output JSON file name for merged nomatches", default='consolidated_po_and_rates_unmerged.json')
    parser.add_argument('--output_consolidated_po_and_rates_matches', type=str, help="Output JSON file name for merged matches, if processing with merge", default='consolidated_po_and_rates_merged.json')

    args = parser.parse_args()

    # Always process nomatches
    fs_nomatches = sorted(glob.glob(f"{args.rate_fit_output_folder}/nomatches/dG/*/*/*pickle"))
    df_nomatches = generate_rates_df(fs_nomatches, args.n_highest, args.library)
    df_nomatches.to_json(args.output_nomatches, orient='records', indent=4)
    print(f"Nomatches output saved to {args.output_nomatches}")

    # Load the other pipeline metrics
    df_po = pd.read_json(args.po_consolidated_results)

    # Merge nomatches with the other pipeline metrics
    df_merged_nomatches = pd.merge(df_nomatches, df_po, how="left", left_on=["name_rt-group", "pH", "library"], right_on=["name_rt-group", "pH", "library"])

    # Select only data from pH6, drop duplicated_timepoints
    df_merged_nomatches.rename({"timepoints_x":"timepoints"}, axis=1, inplace=True)
    df_merged_nomatches.drop(labels="timepoints_y", axis=1, inplace=True)
    df_merged_nomatches = df_merged_nomatches.query("pH == 'pH6'").reset_index(drop=True)

#    pdb.set_trace()

    df_merged_nomatches = process_dataframe_unmerged(df_merged_nomatches)

    df_merged_nomatches.to_json(args.output_consolidated_po_and_rates_nomatches, orient='records', indent=4)
    print(f"Merged nomatches output saved to {args.output_consolidated_po_and_rates_nomatches}")

    # Optionally process matches if --merge is specified
    if args.merge and args.output_matches:
        fs_matches = sorted(glob.glob(f"{args.rate_fit_output_folder}/dG_*/*/*pickle"))
        df_matches = generate_rates_df_with_merge(fs_matches, args.n_highest, args.library)
        df_matches.to_json(args.output_matches, orient='records', indent=4)
        print(f"Matches output saved to {args.output_matches}")

        # Merge matches with the other pipeline metrics
        df_merged_pH6 = pd.merge(df_matches.drop(labels=["name_rt-group"], axis=1),
                                 df_po.query("pH == 'pH6'").drop(labels=["pH", "timepoints"], axis=1),
                                 left_on=["name_rt-group_pH6", "library"],
                                 right_on=["name_rt-group", "library"],
                                 how="left")

        new_cols_pH6 = [i + "_pH6" for i in df_merged_pH6.keys()[len(df_matches.drop(labels=["name_rt-group"], axis=1).columns):]]
        df_merged_pH6.columns = df_matches.drop(labels=["name_rt-group"], axis=1).columns.tolist() + new_cols_pH6

        df_merged_pH6_pH9 = pd.merge(df_merged_pH6,
                                     df_po.query("pH == 'pH9'").drop(labels=["pH", "timepoints"], axis=1),
                                     left_on=["name_rt-group_pH9", "library"],
                                     right_on=["name_rt-group", "library"],
                                     how="left")

        new_cols_pH6_pH9 = [i + "_pH9" for i in df_merged_pH6_pH9.keys()[len(df_merged_pH6.columns):]]
        df_merged_pH6_pH9.columns = df_merged_pH6.columns.tolist() + new_cols_pH6_pH9

        print(f"Duplicated columns: {df_merged_pH6_pH9.columns[df_merged_pH6_pH9.columns.duplicated()]}")
        df_merged_pH6_pH9 = df_merged_pH6_pH9.loc[:, ~df_merged_pH6_pH9.columns.duplicated(keep='first')]

        df_merged_pH6_pH9 = process_dataframe_merged(df_merged_pH6_pH9)

        # df_merged_pH6_pH9[~df_merged_pH6_pH9.isna().any(axis=1)].to_json(args.consolidated_po_and_rates_matches, orient='records', indent=4)
        df_merged_pH6_pH9.to_json(args.output_consolidated_po_and_rates_matches, orient='records', indent=4)
        print(f"Merged matches output saved to {args.output_consolidated_po_and_rates_matches}")

if __name__ == "__main__":
    main()
