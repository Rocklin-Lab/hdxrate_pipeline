import os
import glob
import argparse
import pandas as pd


def gen_sample_list(hxms_dist_fpath_top_dir, file_delimiting_string, library_info_fpath, output_path):
    """
    generate sample list
    :param hxms_dist_fpath_top_dir: top directory where the hx dist files are located
    :param file_delimiting_string: delimiting string anything after prot name
    :param library_info_fpath: library info .json filepath
    :param output_path: output path for sample.csv
    :return: None
    """

    hx_ms_dist_file_list = glob.glob(hxms_dist_fpath_top_dir+'/*/*'+file_delimiting_string)

    library_info_df = pd.read_json(library_info_fpath)

    header = 'name,sequence,hx_dist_fpath\n'

    data_string = ''

    for ind, hx_ms_dist_fpath in enumerate(hx_ms_dist_file_list):

        dpath, hx_ms_dist_fname = os.path.split(hx_ms_dist_fpath)
        prot_name = hx_ms_dist_fname.strip(file_delimiting_string)

        if prot_name in library_info_df['name'].values:

            prot_seq = library_info_df[library_info_df['name'] == prot_name]['sequence'].values[0]

            data_string += '{},{},{}\n'.format(prot_name, prot_seq, hx_ms_dist_fpath)

    with open(output_path, 'w') as outfile:
        outfile.write(header + data_string)
        outfile.close()


def gen_parser_arguments():
    """
    generate commandline arguements to run the hx rate fitting algorithm
    :return:parser
    """
    parser = argparse.ArgumentParser(prog='SAMPLE CSV', description='Generate Sample list for running hx rate fitting')
    parser.add_argument('-t', '--topdir', help='top directory to look for hx ms dist file paths',
                        default='../../workfolder/input_hx_dist')
    parser.add_argument('-d', '--delim', help='delimiting string to look for identifying hx ms dist file path',
                        default='.winner.cpickle.zlib.csv')
    parser.add_argument('-j', '--json', help='library info .json filepath',
                        default='../../library_info.json')
    parser.add_argument('-o', '--outpath', help='output path for sample.csv',
                        default='../../output/sample.csv')
    return parser


def gen_sample_list_from_parser(parser):

    options = parser.parse_args()

    gen_sample_list(hxms_dist_fpath_top_dir=options.topdir,
                    file_delimiting_string=options.delim,
                    library_info_fpath=options.json,
                    output_path=options.outpath)


if __name__ == '__main__':

    parser_ = gen_parser_arguments()
    gen_sample_list_from_parser(parser_)
