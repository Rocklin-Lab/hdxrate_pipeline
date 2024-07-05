# **hxrate**

Set of scripts to obtain set of rates from HX MS data


## **RUN DISTRIBUTED JOBS IN HPC CLUSTERS**

## **RUN VIA SNAKEMAKE**

Create a virtual environment

`conda create -n hdxrate python=3.9`

Activate the environment

`conda activate hdxrate`

Install hxrate package

`python -m pip install hxrate`

Install snakemake

`python -m pip install snakemake`


### **Tips**
You can always do a dry run to make sure the pipeline will work with the way you've set up config file.

for dry run: $ snakemake -s Snakefile -j 1 --dry-run


if there are no errors, you're good to go. If there are errors, make sure everything is set up properly.

If you suspect the error is not your fault and there might be something wrong with the pipeline, contact me and I can help.


## **Single pH data**

configfile: config/config.yml. Config file has all the parameters set to run rate fitting + dg calculation.

Snakefile: workfoler/Snakefile. Make sure to change the path to config file in the begining of the file to the correct one before running


`nohup snakemake -s Snakefile -j 1000 --keep-going --cluster "sbatch -A p30802 -p short -N 1 -n {resources.cpus} --mem=4GB -t 04:00:00" --max-jobs-per-second 3 > nohup_snakefile.out &`


## **Low and high pH data**

configfile: config/config_merge.yml. This config file has all the parameters set to run 1) merge low and high ph data, 2) rate fitting, and 3) dg calculation.

Snakefile: workfoler/Merge_Snakefile & workfoler/Snakefile_nomatches. Run Merge_Snakefile before Snakefile_nomatches.


Merge_Snakefile looks for matches in proteins between two pH runs and continues with merging -> rate fitting -> dg calculation.

Snakefile_nomatches looks for proteins that didn't match and runs rate fitting -> dg calculation

Merge_Snakefile will create backexchange correction files which will be necessary to run Snakefile_nomatches so make sure you run Snakefile_nomatches after backexchange correction files are produced


`nohup snakemake -s Merge_Snakefile -j 1000 --keep-going --cluster "sbatch -A p30802 -p short -N 1 -n {resources.cpus} --mem=4GB -t 04:00:00" --max-jobs-per-second 3 > nohup_merge_snakefile.out &`


`nohup snakemake -s Snakefile_nomatches -j 1000 --keep-going --cluster "sbatch -A p30802 -p short -N 1 -n {resources.cpus} --mem=4GB -t 04:00:00" --max-jobs-per-second 3 > nohup_snakefile_nomatches.out &`



# **Config Files**

## **config.yml**

## filepaths
path_to_repo: path to this repo (path to hxrate)

hx_ms_dist_fpaths_top_dir: path to directory where hx ms files are located

levels_to_fpaths: where the files are present under the directory. if its directly under the top dir, level is 1. if the files are top_dir/PROTEIN_RT_NAME/PROTEIN_RT_NAME_hxms.csv then the level is 2.

library_info_json_fpath: library_info_json_fpath .json from HDX LIMIT Pipeline

protein_rt_column_name: name of column from which protein_rt to be extracted

output_dirpath: output directory path to specify. It doesn't need to actually exist but you need to specify where you would like the output files to go.


## backexchange
backexchange_correction: bool. Whether to calculate backexchange_correction for timepoints

rate_tol: max tolerance for change in mass

min_num_points: at least 5 protein_rt instances need to have similar mass rate to calculate backexchange_correction for that specific timepoint.

change_rate_threshold: threshold for rate of mass change

backexchange_correction_fpath: file path .csv for backexchange_correction if already computed. Keep it empty if not already computed


## exp params
d2o_fraction: fraction of d2o used in experiment. 0.95

d2o_purity: purity of d2o used in experiment. 0.95


## rate_fitting_parameters
adjust_backexchange: Bool. Whether to adjust backexchange is the slowest rate is slower by >= 1.6.

sample_backexchange: Bool. Whether to sample backexchange during rate fitting procedure

num_chains: int. Number of MCMC chains to use during rate fitting optimization

num_warmups: int. Number of initial iteration that doesn't contribute to the samples collected in the MCMC process.

num_samples: int. Number of samples in the MCMC process.


## dg calc parameters
pH: float. ph value

temp: float. temp in Kelvin

nterm: str. Any n terminal addition to the sequence that is not present in the structure file (example: 'HM'). Be careful to use this as it will add nterm to all the proteins in the pipeline. Keep it empty if you don't want additions.

cterm: str. Any c terminal addition to the sequence that is not present in the structure file (example: 'GS'). Be careful to use this as it will add cterm to all the proteins in the pipeline. Keep it empty if you don't want additions.

net_charge_corr: Bool. Whether to use net charge correction while calculating free energy.

anneal_time: float. Time to run the anneal in minutes.

anneal_update_interval: int. Number of steps to save the trajectory information.

## merge parameters
merge_rt_window: window to consider proteins for merging if protein name matches.
