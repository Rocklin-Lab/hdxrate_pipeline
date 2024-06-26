import argparse
import pandas as pd


def concatenate_dfs(file_paths):
    dataframes = [pd.read_json(path) for path in file_paths]
    concatenated_df = pd.concat(dataframes, ignore_index=True)
    return concatenated_df


def main():
    parser = argparse.ArgumentParser(description="Concatenate multiple JSON files containing pandas DataFrames.")
    parser.add_argument('json_files', nargs='+', help="List of JSON file paths to concatenate")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output JSON file name")

    args = parser.parse_args()
    result_df = concatenate_dfs(args.json_files)
    result_df.to_json(args.output, orient='records', indent=4)
    print(f"Concatenated DataFrame saved to {args.output}")

    print(f"WARNING: confirm there is a library column and a pH column in the final dataframe and that their values \n\
    correspond to the [low/high]_ph_exp_label and library keys you indicated in the \n\
    config.yml before using it in the next step.")



if __name__ == "__main__":
    main()