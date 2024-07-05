from hxrate.hxdata import load_pickle_object
from hxrate.openfe import dg_calc


def dg_calc_from_file(hxrate_pickle_fpath: str,
                      temp: float,
                      ph: float,
                      d2o_purity: float,
                      d2o_fraction: float,
                      netcharge_corr: bool = True,
                      min_fe_val: float = -2.0,
                      output_picklepath: str = None,
                      dg_csv_fpath: str = None,
                      dg_plot_fpath: str = None,
                      retun_flag: bool = False):
    """

    :param min_fe_val:
    :param merge_factor_fpath:
    :param hxrate_pickle_fpath:
    :param temp:
    :param ph:
    :param netcharge_corr:
    :param sort_min_val:
    :param output_picklepath:
    :param dg_csv_fpath:
    :param dg_plot_fpath:
    :param retun_flag:
    :return:
    """

    print("AF", temp, ph, d2o_purity, d2o_fraction, netcharge_corr, min_fe_val, output_picklepath)


    hxrate_obj_ = load_pickle_object(hxrate_pickle_fpath)

    dg_output = dg_calc(sequence=hxrate_obj_['sequence'],
                        measured_hx_rates=hxrate_obj_['bayes_sample']['rate']['mean'],
                        temp=temp,
                        ph=ph,
                        netcharge_corr=netcharge_corr,
                        min_fe_val=min_fe_val,
                       d2o_purity=d2o_purity,
                       d2o_fraction=d2o_fraction)

    # save attributes here

    dg_output.hxrate_output = hxrate_obj_

    if output_picklepath is not None:
        dg_output.to_pickle(filepath=output_picklepath)

    if dg_csv_fpath is not None:
        dg_output.to_csv(filepath=dg_csv_fpath)

    if dg_plot_fpath is not None:
        dg_output.plot_dg(filepath=dg_plot_fpath)

    if retun_flag:
        return dg_output


def gen_parser_args():

    import argparse

    parser_ = argparse.ArgumentParser(prog='DG Calculation')
    parser_.add_argument('-i', '--input_', type=str, help='HX rate fit .pickle file path')
    parser_.add_argument('-t', '--temp', type=float, default=295, help='temperature in K')
    parser_.add_argument('-p', '--ph', type=float, default=6.0, help='ph')
    parser_.add_argument('-m', '--minfe', type=float, default=-2.0, help='min fe value')
    parser_.add_argument('-n', '--netcharge', default=True, action=argparse.BooleanOptionalAction)
    parser_.add_argument('-opk', '--output_pickle', type=str, help='dg output .pickle file path')
    parser_.add_argument('-oc', '--output_csv', type=str, help='dg output .csv file path')
    parser_.add_argument('-opd', '--output_pdf', type=str, help='dg output plot .pdf file path')
    parser_.add_argument('-dp', '--d2o_purity', type=float, help='d2o fraction on buffer solution')
    parser_.add_argument('-df', '--d2o_fraction', type=float, help='d2o buffer fraction when mixed with protein')


    return parser_


def run_from_parser():

    parser_ = gen_parser_args()

    options = parser_.parse_args()

    dg_calc_from_file(hxrate_pickle_fpath=options.input_,
                      temp=options.temp,
                      ph=options.ph,
                      netcharge_corr=options.netcharge,
                      min_fe_val=options.minfe,
                      d2o_purity=options.d2o_purity,
                      d2o_fraction=options.d2o_fraction,
                      output_picklepath=options.output_pickle,
                      dg_csv_fpath=options.output_csv,
                      dg_plot_fpath=options.output_pdf,
                      retun_flag=False)


if __name__ == '__main__':

    run_from_parser()
