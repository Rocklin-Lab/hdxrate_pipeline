import os
import pandas as pd
import numpy as np
import argparse


def check_hxdist_fpath(hxms_dir, level_to_fpaths, prot_rt_name, file_end_string):

    file_valid = False

    hx_fpath = hxms_dir + ('/' + prot_rt_name) * level_to_fpaths + file_end_string

    # check if fpath exists
    if os.path.exists(hx_fpath):
        # check if the file size is greater than 0
        if os.stat(hx_fpath).st_size > 0:
            hx_df = pd.read_csv(hx_fpath)
            # check if the file is empty
            if not hx_df.empty:
                file_valid = True

    return hx_fpath, file_valid


def gen_prot_name_from_rt_name(rt_name):
    prot_name_chars = rt_name.split('_')
    del prot_name_chars[-1]
    prot_name = '_'.join(x for x in prot_name_chars)
    if prot_name.endswith('.pdb'):
        prot_name = prot_name.strip('.pdb')
    return prot_name


def find_no_match(prot_name_rt_arr, prot_name_rt_arr_match):

    no_match = []

    for name in prot_name_rt_arr:
        if name not in prot_name_rt_arr_match:
            no_match.append(name)

    return no_match


def gen_match_protein_rt_list(low_ph_libinfo_fpath,
                              high_ph_libinfo_fpath,
                              low_ph_hxms_dir,
                              low_ph_level_to_fpaths,
                              low_ph_hxms_end_string,
                              high_ph_hxms_dir,
                              high_ph_level_to_fpaths,
                              high_ph_hxms_end_string,
                              rt_window=2):

    low_ph_libinfo_df = pd.read_json(low_ph_libinfo_fpath)
    high_ph_libinfo_df = pd.read_json(high_ph_libinfo_fpath)

    # generate matching protein names from the library info first
    low_ph_name_rt_arr = np.unique(low_ph_libinfo_df['name'].values)
    high_ph_name_rt_arr = np.unique(high_ph_libinfo_df['name'].values)

    low_ph_name_list = ['_'.join(x.split('_')[:-1]) for x in low_ph_name_rt_arr]
    high_ph_name_list = ['_'.join(x.split('_')[:-1]) for x in high_ph_name_rt_arr]

    low_ph_match_rtname_list = []
    high_ph_match_rtname_list = []

    for ind, low_ph_name_ in enumerate(low_ph_name_list):

        low_ph_rt_name = low_ph_name_rt_arr[ind]
        low_ph_rt = float(low_ph_rt_name.split('_')[-1])
        high_ph_rt_name_matches = [high_ph_name_rt_arr[ind2] for ind2, x in enumerate(high_ph_name_list) if x == low_ph_name_]

        if len(high_ph_rt_name_matches) > 0:

            for high_ph_rt_name in high_ph_rt_name_matches:

                high_ph_rt = float(high_ph_rt_name.split('_')[-1])

                if abs(high_ph_rt - low_ph_rt) <= rt_window:

                    lowph_fpath, lowph_file_valid = check_hxdist_fpath(hxms_dir=low_ph_hxms_dir,
                                                                       level_to_fpaths=low_ph_level_to_fpaths,
                                                                       prot_rt_name=low_ph_rt_name,
                                                                       file_end_string=low_ph_hxms_end_string)

                    highph_fpath, highph_file_valid = check_hxdist_fpath(hxms_dir=high_ph_hxms_dir,
                                                                         level_to_fpaths=high_ph_level_to_fpaths,
                                                                         prot_rt_name=high_ph_rt_name,
                                                                         file_end_string=high_ph_hxms_end_string)

                    if lowph_file_valid:
                        if highph_file_valid:
                            low_ph_match_rtname_list.append(low_ph_rt_name)
                            high_ph_match_rtname_list.append(high_ph_rt_name)

    return low_ph_match_rtname_list, high_ph_match_rtname_list


def gen_nomatch_protein_list(low_ph_libinfo_fpath,
                             high_ph_libinfo_fpath,
                             low_ph_hxms_dir,
                             low_ph_level_to_fpaths,
                             low_ph_hxms_end_string,
                             high_ph_hxms_dir,
                             high_ph_level_to_fpaths,
                             high_ph_hxms_end_string,
                             rt_window=2,
                             output_path=None,
                             return_flag=False):

    low_ph_libinfo_df = pd.read_json(low_ph_libinfo_fpath)
    high_ph_libinfo_df = pd.read_json(high_ph_libinfo_fpath)

    # generate matching protein names from the library info first
    low_ph_name_rt_arr = np.unique(low_ph_libinfo_df['name'].values)
    high_ph_name_rt_arr = np.unique(high_ph_libinfo_df['name'].values)

    match_lowph_rtname, match_highph_rtname = gen_match_protein_rt_list(low_ph_libinfo_fpath=low_ph_libinfo_fpath,
                                                                        high_ph_libinfo_fpath=high_ph_libinfo_fpath,
                                                                        low_ph_hxms_dir=low_ph_hxms_dir,
                                                                        low_ph_level_to_fpaths=low_ph_level_to_fpaths,
                                                                        low_ph_hxms_end_string=low_ph_hxms_end_string,
                                                                        high_ph_hxms_dir=high_ph_hxms_dir,
                                                                        high_ph_level_to_fpaths=high_ph_level_to_fpaths,
                                                                        high_ph_hxms_end_string=high_ph_hxms_end_string,
                                                                        rt_window=rt_window)

    low_ph_rtname_nomatch_list = find_no_match(prot_name_rt_arr=low_ph_name_rt_arr,
                                               prot_name_rt_arr_match=match_lowph_rtname)
    high_ph_rtname_nomatch_list = find_no_match(prot_name_rt_arr=high_ph_name_rt_arr,
                                                prot_name_rt_arr_match=match_highph_rtname)

    if len(low_ph_rtname_nomatch_list) == 0:
        low_ph_rtname_nomatch_list = [match_lowph_rtname[0]]
    if len(high_ph_rtname_nomatch_list) == 0:
        high_ph_rtname_nomatch_list = [match_highph_rtname[0]]

    hxfpath_list = []
    prot_rt_name_list = []
    prot_name_list = []
    prot_seq_list = []
    ph_str_list = []

    for lowphrtname_nomatch in low_ph_rtname_nomatch_list:
        lowphrt_fpath, lowphrt_file_valid = check_hxdist_fpath(prot_rt_name=lowphrtname_nomatch,
                                                               hxms_dir=low_ph_hxms_dir,
                                                               level_to_fpaths=low_ph_level_to_fpaths,
                                                               file_end_string=low_ph_hxms_end_string)
        if lowphrt_file_valid:
            lowph_prot_seq = low_ph_libinfo_df[low_ph_libinfo_df['name'] == lowphrtname_nomatch]['sequence'].values[0]
            prot_seq_list.append(lowph_prot_seq)
            prot_name_list.append(gen_prot_name_from_rt_name(lowphrtname_nomatch))
            hxfpath_list.append(lowphrt_fpath)
            prot_rt_name_list.append(lowphrtname_nomatch)
            ph_str_list.append('low_ph')

    for highphrtname_nomatch in high_ph_rtname_nomatch_list:
        highphrt_fpath, highphrt_file_valid = check_hxdist_fpath(prot_rt_name=highphrtname_nomatch,
                                                                 hxms_dir=high_ph_hxms_dir,
                                                                 level_to_fpaths=high_ph_level_to_fpaths,
                                                                 file_end_string=high_ph_hxms_end_string)
        if highphrt_file_valid:
            highph_prot_seq = high_ph_libinfo_df[high_ph_libinfo_df['name'] == highphrtname_nomatch]['sequence'].values[0]
            prot_seq_list.append(highph_prot_seq)
            prot_name_list.append(gen_prot_name_from_rt_name(highphrtname_nomatch))
            hxfpath_list.append(highphrt_fpath)
            prot_rt_name_list.append(highphrtname_nomatch)
            ph_str_list.append('high_ph')

    if output_path is not None:

        with open(output_path, 'w') as outfile:

            header = 'ph_str,prot_rt_name,prot_name,prot_seq,hx_ms_fpath\n'
            outfile.write(header)

            for num in range(len(prot_rt_name_list)):

                line = '{},{},{},{},{}\n'.format(ph_str_list[num],
                                                 prot_rt_name_list[num],
                                                 prot_name_list[num],
                                                 prot_seq_list[num],
                                                 hxfpath_list[num])
                outfile.write(line)

            outfile.close()

    if return_flag:
        outdict = dict()
        outdict['prot_rt_name_list'] = prot_rt_name_list
        outdict['prot_name_list'] = prot_name_list
        outdict['prot_seq_list'] = prot_seq_list
        outdict['ph_str_list'] = ph_str_list
        outdict['hx_ms_fpath_list'] = hxfpath_list


def gen_parser():

    parser_ = argparse.ArgumentParser(prog='GEN NO MATCH PROTEIN SAMPLE LIST')

    parser_.add_argument('-llib', '--lowlib', type=str, help='low ph library info filepath .json')
    parser_.add_argument('-hlib', '--highlib', type=str, help='high ph library info filepath .json')
    parser_.add_argument('-lfdir', '--lowfilesdirpath', type=str, help='low ph hxms files dirpath')
    parser_.add_argument('-hfdir', '--highfilesdirpath', type=str, help='high ph hxms files dirpath')
    parser_.add_argument('-lfl', '--lowfilelevels', type=int, help='level under low ph dirpath to find hxms files')
    parser_.add_argument('-hfl', '--highfilelevels', type=int, help='level under high ph dirpath to find hxms files')
    parser_.add_argument('-les', '--lowendstr', type=str, help='low files end string')
    parser_.add_argument('-hes', '--highendstr', type=str, help='high files end string')
    parser_.add_argument('-rw', '--rtwindow', type=float, help='rt window to consider matches')
    parser_.add_argument('-of', '--outputfile', type=str, help='output file path .csv')

    return parser_


def run_from_parser():

    parser = gen_parser()
    options = parser.parse_args()

    gen_nomatch_protein_list(low_ph_libinfo_fpath=options.lowlib,
                             high_ph_libinfo_fpath=options.highlib,
                             low_ph_hxms_dir=options.lowfilesdirpath,
                             low_ph_level_to_fpaths=options.lowfilelevels,
                             low_ph_hxms_end_string=options.lowendstr,
                             high_ph_hxms_dir=options.highfilesdirpath,
                             high_ph_level_to_fpaths=options.highfilelevels,
                             high_ph_hxms_end_string=options.highendstr,
                             rt_window=options.rtwindow,
                             output_path=options.outputfile,
                             return_flag=False)


if __name__ == '__main__':

    run_from_parser()
