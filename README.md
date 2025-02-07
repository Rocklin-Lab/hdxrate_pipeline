## version 1.25.35

# **HDX Rate Pipeline with Snakemake**

## **Overview**
This pipeline processes HDX MS outputs from `mhdx-pipeline` (https://github.com/Rocklin-Lab/mhdx-pipeline) and computes exchange rates. It supports single-pH (pH6) and dual-pH (pH6 & pH9) experiments.

### **Folder Organization** (for dual-pH experiments)
```
{library_name}/
  ├── {date_ph6_data}_{library_name}_pH6/
  │   ├── {mhdx-pipeline}/
  ├── {date_ph9_data}_{library_name}_pH9/
      ├── {mhdx-pipeline}/
```

---
## **Setup**

### **1. Clone Repositories**
```bash
git clone https://github.com/Rocklin-Lab/hdxrate-tools.git
git clone https://github.com/Rocklin-Lab/hdxrate-pipeline.git
```

### **2. Create and Activate Environment**
```bash
conda create -n hdxrate python=3.9 snakemake=7.26.0
conda activate hdxrate
python -m pip install ./hxrate-tools
```

### **3. Prepare Config File**
Copy the config file to your working directory:
```bash
cp {path-to}/hdxrate-pipeline/config.yml {path-to}/{library_name}
```

### **4. Organize Data**
For **single-pH** experiments, copy the consolidated results:
```bash
cp {library_name}/{date-of-experiment}_{library_name}_pH6/{mhdx-pipeline}/resources/10_ic_time_series/consolidated_results.json .
```
For **dual-pH** experiments, concatenate JSON files:
```bash
python {path-to}/hdxrate-pipeline/auxiliar/concatenate_dataframes.py \
    {library_name}/{date_ph6_data}_{library_name}_pH6/{mhdx-pipeline}/resources/10_ic_time_series/consolidated_results.json \
    {library_name}/{date_ph9_data}_{library_name}_pH9/{mhdx-pipeline}/resources/10_ic_time_series/consolidated_results.json \
    --output {date_ph6_data}_{date_ph9_data}_po_results.json
```
**Important:** Preserve the prefix `{date_ph6_data}_{date_ph9_data}_po_results.json` for post-processing.

---
## **Running the Pipeline**
### **1. Configure `config.yml`**
Modify these key parameters:
```yaml
path_to_repo: "{path-to}/hdxrate-pipeline"
path_to_filtered_data: "{path-to}/{library}/{date_ph6_data}_{date_ph9_data}_po_results.json"
library: "{library-name}"  # Must match exactly in consolidated results
output_dirpath: "{date_ph6_data}_{date_ph9_data}_rate_fit_output"
low_ph_library_info: "{path-to}/{library_name}/{date_ph6_data}_{library_name}_pH6/mhdx-pipeline/resources/7_idotp_filter/checked_library_info.json"
high_ph_library_info: "{path-to}/{library_name}/{date_ph9_data}_{library_name}_pH9/mhdx-pipeline/resources/7_idotp_filter/checked_library_info.json"
```

### **2. Run the Snakemake Workflow**
#### **Step 1: Backexchange Correction**
```bash
snakemake -s {path-to}/hdxrate_pipeline/snakefiles/1_Snakefile_twopHs_bx -j 1000 --keep-going \
    --cluster "sbatch -A p31346 -p short -N 1 -n {resources.cpus} --mem=4GB -t 04:00:00" --max-jobs-per-second 3
```
#### **Step 2: Rate Fitting**
```bash
snakemake -s {path-to}/hxrate_pipeline/snakefiles/2_Snakefile_twopHs_nomatches -j 1000 --keep-going \
    --cluster "sbatch -A p31346 -p short -N 1 -n {resources.cpus} --mem=4GB -t 04:00:00" --max-jobs-per-second 3
```
#### **Step 3: Delta G Calculation**
```bash
snakemake -s {path-to}/hxrate_pipeline/snakefiles/3_Snakefile_twopHs_merge -j 1000 --keep-going \
    --cluster "sbatch -A p30802 -p short -N 1 -n {resources.cpus} --mem=4GB -t 04:00:00" --max-jobs-per-second 3
```
**Note:** If using only pH6 data, replace `twopHs` with `singlepH` in the Snakefile names. Customize `--cluster` specifications for your cluster, or omit in case you are using a local computer

#### **Step 4: Consolidate Results**
```bash
snakemake -s {path-to}/hxrate_pipeline/snakefiles/4_Snakefile_consolidate_results -j 1000 --keep-going \
    --cluster "sbatch -A p30802 -p short -N 1 -n {resources.cpus} --mem=4GB -t 04:00:00" --max-jobs-per-second 3
```

### **Expected final results**
```bash
{output_dirpath}/consolidated_results/unfiltered.json : Unfiltered rate fitting results
{output_dirpath}/consolidated_results/filtered.json : All sequences passing filtering criteria
{output_dirpath}/consolidated_results/deduplicated.json : Single best results per sequence
```

---
## **Advanced options configuration file (`config.yml`)**

### **Backexchange Correction Parameters**
```yaml
backexchange_correction: true  # Enable backexchange correction
rate_tol: 0.06  # Max tolerance for mass change
min_num_points: 5  # Minimum required data points
change_rate_threshold: 0.2  # Threshold for mass rate change
```

### **Experimental Parameters**
```yaml
d2o_fraction: 0.90 # sample dilution into deuterated buffer 1 sample : 9 d2o buffer
d2o_purity: 0.95 # buffer D2O purity
```

### **Rate Fitting Parameters**
```yaml
adjust_backexchange: True  # Adjust if slowest rate is 1.6x slower
sample_backexchange: False  # Sample backexchange during fitting
num_chains: 4  # MCMC chains for rate fitting
num_warmups: 100  # MCMC warmup iterations
num_samples: 250  # Number of samples for MCMC
```

### **Delta G Calculation Parameters**
```yaml
pH: 6.0  # Experiment pH
temp: 293  # Temperature in Kelvin
nterm: ""  # N-terminal addition (optional)
cterm: ""  # C-terminal addition (optional)
net_charge_corr: True  # Use net charge correction
```

### **Merging Parameters**
```yaml
merge_rt_window: 1.0  # Window to merge proteins with matching names in minutes
```

---
## **Final Notes**
- Ensure proper paths are set before running Snakemake commands.
- Adjust computational resources (`--mem=4GB`, `-t 04:00:00`) based on your HPC setup.
- For troubleshooting, check Snakemake logs and intermediate outputs.
- For questions, contact `ajrferrari@gmail.com` or open an issue in the repository.

