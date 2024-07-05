# script to generate job submitting .sh files
import os
import sys
import subprocess
import argparse
import numpy as np
import yaml
import pandas as pd
import random


def gen_bash_script(job_params_dict: dict,
                    hx_rate_params_fpath: str,
                    hx_mass_dist_input_fpath: str,
                    protein_name: str,
                    protein_sequence: str,
                    output_dirpath: str,
                    sbatch_output_fpath: str,
                    sbatch_job_name: str) -> str:
    """
    generate bash script string to be saved for quest jobs
    :param job_params_dict: job params dictionary
    :param hx_rate_params_fpath: hx rate params.yml path
    :param hx_mass_dist_input_fpath: hx mass dist .csv path
    :param protein_name: protein name
    :param protein_sequence: protein sequence
    :param output_dirpath: output dirpath
    :param sbatch_output_fpath:
    :param sbatch_error_fpath:
    :param sbatch_job_name:
    :return: bash script string
    """

    sbatch_string = '#SBATCH '
    new_line = '\n'

    bash_string = ''

    # sbatch strings
    bash_string += '#!/bin/bash' + new_line
    bash_string += sbatch_string + '-A ' + job_params_dict['account'] + new_line
    bash_string += sbatch_string + '-p ' + job_params_dict['partition'] + new_line
    bash_string += sbatch_string + '-t ' + job_params_dict['time'] + new_line
    bash_string += sbatch_string + '-N ' + job_params_dict['num_nodes'] + new_line
    bash_string += sbatch_string + '--mem=' + job_params_dict['memory'] + new_line
    bash_string += sbatch_string + '--ntasks-per-node=' + job_params_dict['ntask_per_node'] + new_line
    if job_params_dict['email'] != 'None':
        if '@' in job_params_dict['email']:
            bash_string += sbatch_string + '--mail-user=' + job_params_dict['email'] + new_line
            bash_string += sbatch_string + '--mail-type=BEGIN,END,FAIL,REQUEUE' + new_line
        else:
            print('Valid email address required (example@email.com)')
            sys.exit(0)

    bash_string += sbatch_string + '--output=' + sbatch_output_fpath + new_line
    bash_string += sbatch_string + '--job-name=' + sbatch_job_name + new_line

    # additional commands
    bash_string += new_line + '#sub commands' + new_line

    for command_line in job_params_dict['sub_commands']:
        bash_string += command_line + new_line + '#' + new_line

    # gen python script
    bash_string += new_line + '#python script command' + new_line
    bash_string += 'python ' + job_params_dict['python_script'] + ' -i ' + hx_mass_dist_input_fpath + ' -s ' + protein_sequence + ' -n ' + protein_name + ' -p ' + hx_rate_params_fpath + ' -o ' + output_dirpath + new_line

    # end of the script
    bash_string += new_line + '#end' + new_line

    return bash_string


def get_jobs_from_sample_list(sample_list_fpath,
                              jobs_params_fpath,
                              sbatch_output_path,
                              run_jobs=False):
    """
    gen bash scripts from the sample script
    :param sample_list_fpath: sample list .csv file path
    :param jobs_params_fpath: sub jobs params .yml file path
    :param run_jobs:bool. if True, will run the command to submit the job
    :return:
    """

    sample_df = pd.read_csv(sample_list_fpath)

    prot_name_list = sample_df.iloc[:, 0].values
    prot_seq_list = sample_df.iloc[:, 1].values
    hx_dist_fpath_list = sample_df.iloc[:, 2].values
    hx_params_fpath_list = sample_df.iloc[:, 3].values
    output_path_list = sample_df.iloc[:, 4].values

    job_params_dict = yaml.load(open(jobs_params_fpath, 'rb'), Loader=yaml.Loader)

    random.seed()
    rand_num = random.randint(1, 9)
    if len(prot_name_list) > 1:
        job_num_range = np.arange(rand_num*1000, rand_num*1000 + len(prot_name_list))
    else:
        job_num_range = [rand_num*1000]

    for ind, (prot_name, prot_seq, hx_dist_fpath, hx_params_fpath, output_path) in enumerate(zip(prot_name_list,
                                                                                                 prot_seq_list,
                                                                                                 hx_dist_fpath_list,
                                                                                                 hx_params_fpath_list,
                                                                                                 output_path_list)):

        job_num = job_num_range[ind]
        job_name = prot_name + '_quest' + str(job_num)

        sbatch_out_file = os.path.join(sbatch_output_path, job_name + '.out')

        bash_script_string = gen_bash_script(job_params_dict=job_params_dict,
                                             hx_rate_params_fpath=hx_params_fpath,
                                             hx_mass_dist_input_fpath=hx_dist_fpath,
                                             protein_name=prot_name,
                                             protein_sequence=prot_seq,
                                             output_dirpath=output_path,
                                             sbatch_output_fpath=sbatch_out_file,
                                             sbatch_job_name=job_name)

        sbatch_sh_file = os.path.join(sbatch_output_path, 'sub_' + job_name + '.sh')

        # write the bash script to sh file
        with open(sbatch_sh_file, 'w') as sh_file:
            sh_file.write(bash_script_string)
            sh_file.close()

        if run_jobs:
            subprocess.run(['sbatch', str(sbatch_sh_file)])
            print('\nran: $sbatch ' + str(sbatch_sh_file) + '\n')


def gen_parser_commands():

    parser = argparse.ArgumentParser(prog='SUB_JOBS', description='Generate sbatch scripts to run hx rate fitting')
    parser.add_argument('-s', '--sample', help='sample list .csv',
                        default='../../workfolder/sample.csv')
    parser.add_argument('-p', '--job_params', help='job params YAML file .yml',
                        default='../../params/sub_jobs_params.yml')
    parser.add_argument('-o', '--output', help='sbatch output path', default='../../workfolder/sub_jobs_')
    parser.add_argument('-r', '--run', action='store_true', help='option to submit the job')

    return parser


def run_sub_job_from_parser(parser):

    options = parser.parse_args()
    sample_list_fpath = options.sample
    job_params_fpath = options.job_params
    sbatch_output_path = options.output
    run_jobs_bool = options.run

    get_jobs_from_sample_list(sample_list_fpath=sample_list_fpath,
                              jobs_params_fpath=job_params_fpath,
                              sbatch_output_path=sbatch_output_path,
                              run_jobs=run_jobs_bool)


if __name__ == '__main__':

    parser_ = gen_parser_commands()
    run_sub_job_from_parser(parser_)
