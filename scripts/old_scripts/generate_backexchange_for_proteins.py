from hxdata import load_data_from_hdx_ms_dist_, write_backexchange_array, load_tp_dependent_dict
from methods import normalize_mass_distribution_array
from hx_rate_fit import calc_back_exchange
import argparse


def write_tp_backexchange_(sequence: str,
                           hxms_filepath: str,
                           d2o_fraction: float,
                           d2o_purity: float,
                           bkex_tp: float = 1e9,
                           usr_backexchange: float = None,
                           backexch_corr_fpath: str = None,
                           output_path: str = None,
                           return_flag: bool = False):
    # read hxms filepath
    tp, iso_dists = load_data_from_hdx_ms_dist_(fpath=hxms_filepath)
    norm_dists = normalize_mass_distribution_array(mass_dist_array=iso_dists)

    # if correction dictionary is provided, it'll use that to generate backexchange for each timepoint
    bkexch_corr_dict = None
    if backexch_corr_fpath is not None:
        bkexch_corr_dict = load_tp_dependent_dict(backexch_corr_fpath)

    backexchange_object = calc_back_exchange(sequence=sequence,
                                             experimental_isotope_dist=norm_dists[-1],
                                             timepoints_array=tp,
                                             d2o_fraction=d2o_fraction,
                                             d2o_purity=d2o_purity,
                                             bkex_tp=bkex_tp,
                                             usr_backexchange=usr_backexchange,
                                             backexchange_corr_dict=bkexch_corr_dict)
    if output_path is not None:
        write_backexchange_array(timepoints=tp,
                                 backexchange_array=backexchange_object.backexchange_array,
                                 output_path=output_path)

    if return_flag:
        return backexchange_object
    else:
        return None


def gen_parser_args():

    parser = argparse.ArgumentParser(prog='WRITE BACKEXCHANGE ARRAY')
    parser.add_argument('-s', '--sequence', help='protein sequence one letter code')
    parser.add_argument('-i', '--inputfile', help='hxms dist file path')
    parser.add_argument('-df', '--dfrac', help='d2o fraction')
    parser.add_argument('-dp', '--dpur', help='d2o purity')
    parser.add_argument('-ub', '--userbackexchange', help='user input backexchange value', default=None)
    parser.add_argument('-bcf', '--bkexch_corr_fpath', help='backexchange correction filepath', default=None)
    parser.add_argument('-o', '--output_path', help='backexchange array output filepath')

    return parser


def run_from_parser():

    parser_ = gen_parser_args()
    options = parser_.parse_args()

    write_tp_backexchange_(sequence=options.sequence,
                           hxms_filepath=options.inputfile,
                           d2o_fraction=options.dfrac,
                           d2o_purity=options.dpur,
                           usr_backexchange=options.ub,
                           backexch_corr_fpath=options.bcf,
                           output_path=options.output_path,
                           return_flag=False)


if __name__ == '__main__':

    run_from_parser()
