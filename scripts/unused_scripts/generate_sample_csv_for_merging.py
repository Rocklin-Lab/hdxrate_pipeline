import pandas as pd
import numpy as np
import glob
from dataclasses import dataclass
import argparse


@dataclass
class ProteinInfo(object):
    charge_state: float = None
    retention_time: float = None
    im_bin: float = None
    thr_mz: float = None
    exp_mz: float = None
    ppm: float = None


@dataclass
class ProtID(object):
    protein_name: str = None
    protein_name_no_ret: str = None
    protein_sequence: str = None
    retention_time: float = None
    num_charge_states: int = None
    protein_info: list = None


def gen_protein_info(protein_only_df):
    protein_info_list = []

    for ind, (charge, rt, im_bin, thr_mz, exp_mz, ppm) in enumerate(zip(protein_only_df['charge'].values,
                                                                        protein_only_df['RT'].values,
                                                                        protein_only_df['im_mono'].values,
                                                                        protein_only_df['expect_mz'].values,
                                                                        protein_only_df['obs_mz'].values,
                                                                        protein_only_df['ppm'].values)):
        prot_info = ProteinInfo(charge_state=float(charge),
                                retention_time=rt,
                                im_bin=im_bin,
                                thr_mz=thr_mz,
                                exp_mz=exp_mz,
                                ppm=ppm)

        protein_info_list.append(prot_info)

    return protein_info_list


def gen_protein_id_list_from_library_info(library_info):

    prot_names = np.unique(library_info['name'].values)

    list_of_protein_id = []

    for prot_name in prot_names:

        prot_df = library_info[library_info['name'] == prot_name]

        prot_name_no_rt = '_'.join(prot_name.split('_')[:-1])
        protein_object = ProtID(protein_name=prot_name,
                                protein_name_no_ret=prot_name_no_rt,
                                protein_sequence=prot_df['sequence'].values[0],
                                retention_time=np.average(prot_df['rt_group_mean_RT_0_0'].values),
                                num_charge_states=len(prot_df['charge'].values),
                                protein_info=gen_protein_info(protein_only_df=prot_df))

        list_of_protein_id.append(protein_object)

    return list_of_protein_id


def match_prot_ids_(low_ph_library_info, high_ph_library_info, rt_search_window=0.5):
    """

    :param low_ph_library_info:
    :param high_ph_library_info:
    :param rt_search_window:
    :return:
    """

    low_ph_protein_ids = gen_protein_id_list_from_library_info(low_ph_library_info)
    high_ph_protein_ids = gen_protein_id_list_from_library_info(high_ph_library_info)

    low_ph_protein_names = [x.protein_name_no_ret for x in low_ph_protein_ids]

    match_list = []

    for ind, low_ph_prot_name in enumerate(low_ph_protein_names):

        low_ph_prot_id = low_ph_protein_ids[ind]
        high_ph_prot_matches = [x for x in high_ph_protein_ids if x.protein_name_no_ret == low_ph_prot_name]

        if len(high_ph_prot_matches) > 0:

            for high_ph_prot_id in high_ph_prot_matches:

                if abs(high_ph_prot_id.retention_time - low_ph_prot_id.retention_time) <= rt_search_window:

                    store_dict = dict()
                    store_dict['low_prot_id'] = low_ph_prot_id
                    store_dict['high_prot_id'] = high_ph_prot_id
                    match_list.append(store_dict)

    return match_list


def generate_match_sample_list(low_ph_top_dir,
                               high_ph_top_dir,
                               low_ph_level_to_files,
                               high_ph_level_to_files,
                               low_ph_file_id_str,
                               high_ph_file_id_str,
                               low_ph_lib_info_fpath,
                               high_ph_lib_info_fpath,
                               rt_search_window,
                               output_path):

    low_ph_lib_info = pd.read_json(low_ph_lib_info_fpath)
    high_ph_lib_info = pd.read_json(high_ph_lib_info_fpath)

    match_list = match_prot_ids_(low_ph_library_info=low_ph_lib_info,
                                 high_ph_library_info=high_ph_lib_info,
                                 rt_search_window=rt_search_window)

    low_level_fpath_string = '/*'*int(low_ph_level_to_files)
    high_level_fpath_string = '/*'*int(high_ph_level_to_files)

    header = 'protein_name,protein_name_low_ph,protein_name_high_ph,sequence,low_ph_data_fpath,high_ph_data_fpath\n'

    data_string = ''

    for match in match_list:

        low_prot_id = match['low_prot_id']
        high_prot_id = match['high_prot_id']

        low_fpath_search_str = low_ph_top_dir + low_level_fpath_string + low_prot_id.protein_name + low_ph_file_id_str
        high_fpath_search_str = high_ph_top_dir + high_level_fpath_string + high_prot_id.protein_name + high_ph_file_id_str

        low_prot_fpath_list = glob.glob(low_fpath_search_str)
        high_prot_fpath_list = glob.glob(high_fpath_search_str)

        if len(low_prot_fpath_list) != 0:
            if len(high_prot_fpath_list) != 0:
                for ind, low_fpath in enumerate(low_prot_fpath_list):
                    for ind2, high_fpath in enumerate(high_prot_fpath_list):
                        data_string += '{},{},{},{},{},{}\n'.format(low_prot_id.protein_name_no_ret,
                                                                    low_prot_id.protein_name,
                                                                    high_prot_id.protein_name,
                                                                    low_prot_id.protein_sequence,
                                                                    low_fpath,
                                                                    high_fpath)
            else:
                print('high ph file not found for protein id %s. Excluding the protein for merging' % high_prot_id.protein_name)
        else:
            print('low ph file not found for protein id %s. Excluding the protein for merging' % low_prot_id.protein_name)

    with open(output_path, 'w') as outfile:
        outfile.write(header + data_string)
        outfile.close()


def gen_parser_args():

    parser = argparse.ArgumentParser(prog='MERGE SAMPLE', description='Generate Sample list for merging low and high ph data')
    parser.add_argument('-ltd', '--lowtopdir', help='low ph top directory for hxms files')
    parser.add_argument('-llv', '--lowlevel', help='level in which to find low ph hxms files in low top directory')
    parser.add_argument('-htd', '--hightopdir', help='high ph top directory for hxms files')
    parser.add_argument('-hlv', '--highlevel', help='level in which to find high ph hxms files in high top directory')
    parser.add_argument('-lds', '--lowdelimstr', help='low ph hxms files delim string')
    parser.add_argument('-hds', '--highdelimstr', help='high ph hxms files delim string')
    parser.add_argument('-lli', '--lowlibinfo', help='low ph library info .json')
    parser.add_argument('-hli', '--highlibinfo', help='high ph library info .json')
    parser.add_argument('-rtw', '--rtwindow', help='rt window inclusion')
    parser.add_argument('-o', '--outputpath', help='merging list output path .csv')

    return parser


def run_from_parser():

    parser_ = gen_parser_args()
    options = parser_.parse_args()

    generate_match_sample_list(low_ph_top_dir=options.lowtopdir,
                               high_ph_top_dir=options.hightopdir,
                               low_ph_level_to_files=options.lowlevel,
                               high_ph_level_to_files=options.highlevel,
                               low_ph_file_id_str=options.lowdelimstr,
                               high_ph_file_id_str=options.highdelimstr,
                               low_ph_lib_info_fpath=options.lowlibinfo,
                               high_ph_lib_info_fpath=options.highlibinfo,
                               rt_search_window=float(options.rtwindow),
                               output_path=options.outputpath)


if __name__ == '__main__':

    run_from_parser()

    # ph6_libinfo_fpath = '/Users/smd4193/OneDrive - Northwestern University/MS_data/2021_lib15_ph6/config_rt_wide_window_extraction_all_prots_tic_cumsum_warp/resources/library_info/library_info.json'
    # ph7_libinfo_fpath = '/Users/smd4193/OneDrive - Northwestern University/MS_data/2021_lib15_ph7/resources/library_info/library_info.json'
    #
    # low_ph_dir = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/lib15_ph6'
    # high_ph_dir = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/lib15_ph7'
    #
    # low_ph_level = 1
    # high_ph_level = 1
    #
    # rt_search_window = 0.5
    #
    # low_ph_file_id_str = '*winner.cpickle.zlib.csv'
    # high_ph_file_id_str = '*winner.cpickle.zlib.csv'
    #
    # output_path = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/merged_data/merge_sample_list.csv'
    #
    # generate_match_sample_list(low_ph_top_dir=low_ph_dir,
    #                            high_ph_top_dir=high_ph_dir,
    #                            low_ph_level_to_files=low_ph_level,
    #                            high_ph_level_to_files=high_ph_level,
    #                            low_ph_file_id_str=low_ph_file_id_str,
    #                            high_ph_file_id_str=high_ph_file_id_str,
    #                            low_ph_lib_info_fpath=ph6_libinfo_fpath,
    #                            high_ph_lib_info_fpath=ph7_libinfo_fpath,
    #                            rt_search_window=rt_search_window,
    #                            output_path=output_path)
