# generate variable back exchange params and sample file
import os
import copy
import yaml
from methods import normalize_mass_distribution_array, make_new_dir
from hxdata import load_data_from_hdx_ms_dist_
from hx_rate_fit import calc_back_exchange


def gen_var_bkexch_params_files(prot_name, prot_sequence, hx_dist_fpath, params_fpath, output_dirpath, num_samples,
                                sample_size=0.5):
    """

    :param sample_csv_file:
    :param output_dirpath:
    :return:
    """

    params_dict = yaml.load(open(params_fpath, 'rb'), Loader=yaml.Loader)

    timepoints, isotope_dist = load_data_from_hdx_ms_dist_(hx_dist_fpath)

    norm_dist = normalize_mass_distribution_array(isotope_dist)

    backexchange = calc_back_exchange(sequence=prot_sequence,
                                      experimental_isotope_dist=norm_dist[-1],
                                      d2o_fraction=params_dict['d2o_fraction'],
                                      d2o_purity=params_dict['d2o_purity'])

    if num_samples % 2 != 0:
        num_samples += 1
    half_num_samples = int(num_samples/2)

    param_dict_list = []

    for num in range(half_num_samples):
        if num == 0:
            new_param_dict = copy.deepcopy(params_dict)
            new_param_dict['usr_backexchange'] = float(backexchange.backexchange_value)
            param_dict_list.append(new_param_dict)
        else:
            bk_add = backexchange.backexchange_value + (num * sample_size/100)
            add_param_dict = copy.deepcopy(params_dict)
            add_param_dict['usr_backexchange'] = float(bk_add)
            bk_sub = backexchange.backexchange_value - (num * sample_size/100)
            if bk_sub > 0:
                sub_param_dict = copy.deepcopy(params_dict)
                sub_param_dict['usr_backexchange'] = float(bk_sub)
                param_dict_list.append(add_param_dict)
                param_dict_list.append(sub_param_dict)

    # generate new param files from the param dict list

    output_dpath = make_new_dir(output_dirpath)

    sample_csv_header = 'name,sequence,hx_dist_fpath,params_fpath,output_dpath\n'
    sample_data_string = ''

    for ind, new_param_dict in enumerate(param_dict_list):
        param_fpath = os.path.join(output_dpath, 'params_' + str(ind) + '.yml')
        sample_data_string += '{},{},{},{},{}\n'.format(prot_name, prot_sequence, hx_dist_fpath, param_fpath, os.path.join(output_dpath, 'params_' + str(ind)))
        with open(param_fpath, 'w') as paramfile:
            yaml.dump(new_param_dict, paramfile)

    sample_fpath = os.path.join(output_dpath, 'sample_var_bkexch.csv')
    with open(sample_fpath, 'w') as samplefile:
        samplefile.write(sample_csv_header+sample_data_string)
        samplefile.close()


if __name__ == '__main__':

    protein_name = 'HEEH_rd4_0097'
    protein_sequence = 'HMDVEEQIRRLEEVLKKNQPVTWNGTTYTDPNEIKKVIEELRKSM'
    hxdist_path = '../../workfolder/input_hx_dist/HEEH_rd4_0097_hx_mass_dist.csv'
    params_path = '../../params/params.yml'
    output_path = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/var_bkexch'
    number_ = 11
    sample_size_ = 0.5

    gen_var_bkexch_params_files(prot_name=protein_name,
                                prot_sequence=protein_sequence,
                                hx_dist_fpath=hxdist_path,
                                params_fpath=params_path,
                                output_dirpath=output_path,
                                num_samples=number_,
                                sample_size=sample_size_)
