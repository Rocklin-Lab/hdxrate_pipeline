import os
import glob
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
import Bio
from Bio.PDB import Polypeptide
from scripts.old_scripts.dG_opt import gen_hbond_data


def get_pdb_files(dirpath):

    pdb_files = glob.glob(dirpath + '/*.pdb')

    return pdb_files


def get_matching_fes_file(dirpath, pdb_fpath):

    pdb_fname = os.path.split(pdb_fpath)[-1]
    fes_fname = pdb_fname[0:-3] + 'fes'
    fes_fpath = os.path.join(dirpath, fes_fname)

    if os.path.isfile(fes_fpath):
        fes_fpath = fes_fpath
    else:
        fes_fpath = None

    return fes_fpath


def get_fes(fes_fpath):

    if os.path.isfile(fes_fpath):
        with open(fes_fpath, 'r') as file:
            fes = [float(x.split()[-1]) for x in file.readlines()]

    return fes


def get_fe_dist_pairs(pdb_fpath, fes_fpath):

    fes = get_fes(fes_fpath=fes_fpath)

    hbond_data = gen_hbond_data(pdb_fpath=pdb_fpath)

    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_fpath, pdb_fpath)
    model = structure[0]

    residues = [x for x in model.get_residues()]

    dists = []
    for ind, res1 in enumerate(residues):
        for ind2, res2 in enumerate(residues):
            if res1.get_resname() == 'GLY' or res2.get_resname() == 'GLY':
                continue
            if hbond_data.hbond_bool_num_list[ind] == 0 or hbond_data.hbond_bool_num_list[ind2] == 0 or ind2 <= ind:
                continue
            if fes[ind] < 0 or fes[ind2] < 0:
                continue
            dist = np.linalg.norm(res1['H'].coord.astype('float') - res2['H'].coord.astype('float'))
            dists.append((dist, abs(fes[ind] - fes[ind2]), ind+1, ind2+1))

    print(pdb_fpath, len(dists))

    return dists


def get_list_of_dists(list_of_pdbs_fes_pairs):

    all_dist_data = []

    for pdb_fes_pair in list_of_pdbs_fes_pairs:

        dists = get_fe_dist_pairs(pdb_fpath=pdb_fes_pair[0],
                                  fes_fpath=pdb_fes_pair[1])
        dists = np.array(dists)
        all_dist_data.append(dists)

    all_dist_data = np.concatenate(all_dist_data)

    return all_dist_data


def make_rect(dirpath):

    pdb_files = get_pdb_files(dirpath=dirpath)

    pdb_fes_pairs = []

    for pdb_fpath in pdb_files:

        fes_fpath = get_matching_fes_file(dirpath=dirpath, pdb_fpath=pdb_fpath)
        if fes_fpath is not None:
            pdb_fes_pairs.append((pdb_fpath, fes_fpath))

    all_dist_data = get_list_of_dists(list_of_pdbs_fes_pairs=pdb_fes_pairs)

    dist = all_dist_data[:, 0]
    fe_diff = all_dist_data[:, 1]

    min_dist, max_dist = min(dist), max(dist)
    min_fe_diff, max_fe_diff = min(fe_diff), max(fe_diff)

    print(min_dist, max_dist)
    print(min_fe_diff, max_fe_diff)

    X, Y = np.mgrid[min_dist:max_dist:100j, min_fe_diff:max_fe_diff:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([dist, fe_diff])

    kernel = gaussian_kde(dataset=values)
    Z = np.reshape(kernel(positions).T, X.shape)

    # this ensures Z is in ascending order
    for num1 in range(100):
        for num2 in range(100):
            Z[num1, num2] = max(Z[num1, num2:])

    plt.figure(figsize=(9, 7))
    ax = plt.subplot(111)
    # ax.contourf(np.rot90(Z))
    plt.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
               extent=[min_dist, max_dist, min_fe_diff, max_fe_diff], aspect='auto')

    ax.set_yticklabels([0,1,2,3,4,5,6,7,8,9],fontsize=30)
    ax.set_xticks([5,10,15,20,25])
    ax.set_xticklabels([5,10,15,20,25],fontsize=30)
    plt.savefig(dirpath + '/neighbor_fe_kde.pdf')


if __name__=="__main__":

    nan_arr = np.array([np.nan, np.nan, np.nan])
    max_nan_ = max(nan_arr)

    hx_rate_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/test/rates/PDB3NGP_12.98532_PDB3NGP_13.56409/PDB3NGP_12.98532_PDB3NGP_13.56409_hx_rate.csv'

    hx_fes = "/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/test/fitting_SA_scorefxn/HEEH_rd4_0097_forwardfolding.fes"
    pdb_fpath = "/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/test/fitting_SA_scorefxn/HEEH_rd4_0097_forwardfolding.pdb"

    dirpath_ = "/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/test/fitting_SA_scorefxn"

    make_rect(dirpath=dirpath_)

    # get_fe_dist_pairs(pdb_fpath=pdb_fpath,
    #                   fes_fpath=hx_fes)


