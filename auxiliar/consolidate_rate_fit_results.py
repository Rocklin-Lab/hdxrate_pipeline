import argparse
import glob
import os
import numpy as np
import pandas as pd

from populate_stats import process_dataframe_unmerged, process_dataframe_merged

def replace_nan_with_zeros(array):
    return np.where(np.isnan(array) | (array < 0), 0, array)

def extract_info_from_dg_pickle_file(f, n_highest, library):
    d = pd.read_pickle(f)
    pH = "pH6" if 'pH6' in f else ("pH9" if 'pH9' in f else None)
    protein_name = d["hxrate_output"]["protein_name"]
    protein_rt_name = d["hxrate_output"]["protein_rt_name"]
    sequence = d["hxrate_output"]["sequence"]
    intrinsic_rates = d["intrinsic_rates"]
    intrinsic_rates_median = d["intrinsic_rates_median"]
    reg_intrinsic_rates = d["reg_intrinsic_rates"]
    netcharge = d["netcharge"]
    chain_diag = d["hxrate_output"]["chain_diagnostics"]
    chain_rmse_list = chain_diag["chain_rmse_list"]
    chain_pass_list = chain_diag["chain_pass_list"]
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

    return [protein_name, library, pH, protein_rt_name, sequence, intrinsic_rates, intrinsic_rates_median,
            reg_intrinsic_rates, netcharge, chain_rmse_list, chain_pass_list, rates_mean, rates_std,
            rates_ci_5, rates_ci_95, rate_fitting_rmse_total, rate_fitting_rmse_per_timepoint, backexchange_per_tp,
            backexchange_final_tp, backexchange_n_res_subtract, backexchange_ind_label, timepoints,
            timepoints_sort_indices, free_energy_, free_energy, dg_max, dg_mean, dg_median, file_path]

def extract_info_from_dg_pickle_file_with_merge(f, n_highest, library):
    d = pd.read_pickle(f)
    protein_name = d["hxrate_output"]["protein_name"]
    sequence = d["hxrate_output"]["sequence"]
    protein_rt_name = d["hxrate_output"]["protein_rt_name"]
    i_split = len(protein_rt_name.split("_")) // 2
    name_rt_group_pH6 = "_".join(protein_rt_name.split("_")[:i_split])
    name_rt_group_pH9 = "_".join(protein_rt_name.split("_")[i_split:])
    intrinsic_rates = d["intrinsic_rates"]
    intrinsic_rates_median = d["intrinsic_rates_median"]
    reg_intrinsic_rates = d["reg_intrinsic_rates"]
    netcharge = d["netcharge"]
    chain_diag = d["hxrate_output"]["chain_diagnostics"]
    chain_rmse_list = chain_diag["chain_rmse_list"]
    chain_pass_list = chain_diag["chain_pass_list"]
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

    return [protein_name, library, protein_rt_name, name_rt_group_pH6, name_rt_group_pH9, sequence, intrinsic_rates,
            intrinsic_rates_median, reg_intrinsic_rates, netcharge, chain_rmse_list, chain_pass_list, rates_mean,
            rates_std, rates_ci_5, rates_ci_95, rate_fitting_rmse_total, rate_fitting_rmse_per_timepoint,
            merge_factor_mean, merge_factor_std, merge_ci_5, merge_ci_95, backexchange_per_tp, backexchange_final_tp,
            backexchange_n_res_subtract, backexchange_ind_label, timepoints, timepoints_sort_indices, free_energy_,
            free_energy, dg_max, dg_mean, dg_median, file_path]

def generate_rates_df(fs, n_highest, library):
    records = [extract_info_from_dg_pickle_file(f, n_highest, library) for f in fs]
    columns = ["name", "library", "pH", "name_rt-group", "sequence", "intrinsic_rates",
               "intrinsic_rates_median", "reg_intrinsic_rates", "netcharge", "chain_rmse_list",
               "chain_pass_list", "rates_mean", "rates_std", "rates_ci_5", "rates_ci_95",
               "rate_fitting_rmse_total", "rate_fitting_rmse_per_timepoint", "backexchange_per_tp",
               "backexchange_final_tp", "backexchange_n_res_subtract", "backexchange_ind_label",
               "timepoints", "timepoints_sort_indices", "free_energy_", "free_energy", "dg_max",
               "dg_mean", "dg_median", "file_path"]
    return pd.DataFrame(records, columns=columns)

def generate_rates_df_with_merge(fs, n_highest, library):
    records = [extract_info_from_dg_pickle_file_with_merge(f, n_highest, library) for f in fs]
    columns = ["name", "library", "name_rt-group", "name_rt-group_pH6", "name_rt-group_pH9", "sequence",
               "intrinsic_rates", "intrinsic_rates_median", "reg_intrinsic_rates", "netcharge",
               "chain_rmse_list", "chain_pass_list", "rates_mean", "rates_std", "rates_ci_5", "rates_ci_95",
               "rate_fitting_rmse_total", "rate_fitting_rmse_per_timepoint", "merge_factor_mean",
               "merge_factor_std", "merge_factor_ci_5", "merge_factor_ci_95", "backexchange_per_tp",
               "backexchange_final_tp", "backexchange_n_res_subtract", "backexchange_ind_label", "timepoints",
               "timepoints_sort_indices", "free_energy_", "free_energy", "dg_max", "dg_mean", "dg_median",
               "file_path"]
    return pd.DataFrame(records, columns=columns)

def main():
    parser = argparse.ArgumentParser(description="Process rates files and generate dataframes.")
    parser.add_argument('-n', '--n_highest', type=int, default=5, help="Number of highest elements to consider for dg_mean calculation")
    parser.add_argument('-l', '--library', type=str, required=True, help="Library identifier")
    parser.add_argument('-m', '--merge', action='store_true', help="Include processing with merging information if specified")
    parser.add_argument('-po', '--po_results', type=str, required=True, help="Specify PO results path")
    parser.add_argument('--rate_output_dir', type=str, required=True, help="Rate fitting output dir")
    args = parser.parse_args()

    output_dir = os.path.join(args.rate_output_dir, "consolidated_results")
    os.makedirs(output_dir, exist_ok=True)

    fs_nomatches = sorted(glob.glob(os.path.join(args.rate_output_dir, "nomatches", "dG", "*", "*", "*pickle")))
    df_nomatches = generate_rates_df(fs_nomatches, args.n_highest, args.library)
    
    df_po = pd.read_json(args.po_results)
    df_po_ph6 = df_po.query('pH == "pH6"').reset_index(drop=True)
    df_po_ph9 = df_po.query('pH == "pH9"').reset_index(drop=True)

    # Create renaming dictionary
    rename_dict_ph6 = {
        col: f"{col}_pH6" for col in df_po_ph6.columns
        if col not in ['pH', 'library']
    }
    
    # Create renaming dictionary
    rename_dict_ph9 = {
        col: f"{col}_pH9" for col in df_po_ph9.columns
        if col not in ['pH', 'library']
    }

    # Apply renaming
    df_po_ph6 = df_po_ph6.rename(columns=rename_dict_ph6)
    df_po_ph9 = df_po_ph9.rename(columns=rename_dict_ph9)


    df_merged_nomatches = pd.merge(df_nomatches, 
                                   df_po_ph6,
                                   how="left", 
                                   left_on=["name_rt-group", "pH", "library"],
                                   right_on=["name_rt-group_pH6", "pH", "library"]).drop('name_rt-group', axis=1).query('pH == "pH6"').reset_index(drop=True)


    df_merged_nomatches_processed = process_dataframe_unmerged(df_merged_nomatches)
    df_merged_nomatches_processed.to_json(os.path.join(output_dir, "PO_and_RateFittingResults_unmerged.json"), orient='records', indent=4)
    print(f"Merged nomatches output saved to {os.path.join(output_dir, 'PO_and_RateFittingResults_unmerged.json')}")

    if args.merge:
        fs_matches = sorted(glob.glob(os.path.join(args.rate_output_dir, "dG_*", "*", "*pickle")))
        df_matches = generate_rates_df_with_merge(fs_matches, args.n_highest, args.library)

        df_merged_pH6 = pd.merge(
            df_matches.drop(columns=["name_rt-group"]),
            df_po_ph6.drop(columns=["pH"]),
            left_on=["name_rt-group_pH6", "library"],
            right_on=["name_rt-group_pH6", "library"],
            how="left"
        )


        df_merged_pH6_pH9 = pd.merge(
            df_merged_pH6,
            df_po_ph9.drop(columns=["pH"]),
            left_on=["name_rt-group_pH9", "library"],
            right_on=["name_rt-group_pH9", "library"],
            how="left"
        )

        df_merged_pH6_pH9 = process_dataframe_merged(df_merged_pH6_pH9)

        df_merged_pH6_pH9.to_json(os.path.join(output_dir, "PO_and_RateFittingResults_merged.json"), orient='records', indent=4)
        print(f"Merged matches output saved to {os.path.join(output_dir, 'PO_and_RateFittingResults_merged.json')}")

if __name__ == "__main__":
    main()
