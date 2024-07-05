import os
import pickle
from methods import calculate_intrinsic_exchange_rates_suggie
import Bio
import numpy as np
import pandas as pd
from Bio.PDB import Polypeptide
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from dataclasses import dataclass
import random
from simanneal import Annealer
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import argparse


@dataclass
class HXRate(object):
    mean: np.ndarray
    median: np.ndarray = None
    std: np.ndarray = None
    ci_5: np.ndarray = None
    ci_95: np.ndarray = None


@dataclass
class DgInputMetadata(object):
    nterm: str = None
    cterm: str = None
    pH: float = None
    temp: float = None
    hbond_length: float = None
    hbond_angle: float = None
    min_free_energy: float = None
    net_charge_corr: bool = True
    r_constant: float = None
    min_comp_free_energy: float = None


@dataclass
class EnergyWeights(object):
    pair_energy: float = 50
    full_burial: float = 120
    hbond_burial: float = 14
    hbond_rank: float = 60
    dist_to_sec_struct: float = 60
    dist_to_nonpolar_res: float = 45
    top_energy_std: float = 10
    comp_deltag_rmse: float = 0


@dataclass
class HBond(object):
    hbond_list: list = None
    all_hbond_list: list = None
    hbond_depth_array: np.ndarray = None
    hbond_bool_list: np.ndarray = None
    hbond_bool_num_list: np.ndarray = None
    hbond_burial: np.ndarray = None
    pairlist: list = None
    hbond_dist_list: np.ndarray = None
    hbond_angle_list: np.ndarray = None
    hbond_partner_list: np.ndarray = None


@dataclass
class DGAnnealData(object):
    pair_energy: float = None
    full_burial_corr: float = None
    hbond_burial_corr: float = None
    hbond_rank_factor: float = None
    distance_to_nonpolar_res_corr: float = None
    distance_to_sec_struct_corr: float = None
    top_stdev: float = None
    comp_deltaG_rmse_term: float = None


@dataclass
class TrajData(object):
    step_num: int = None
    anneal_data: object = None


class DgInput(object):

    def __init__(self,
                 hxrate_obj,
                 pH,
                 temp,
                 pdb_fpath,
                 dg_intpol_fpath,
                 comp_dg_fpath=None,
                 nter='',
                 cter='',
                 hbond_length=2.7,
                 hbond_angle=120,
                 min_free_energy=-10,
                 net_charge_corr=True,
                 r_constant=1.9872036e-3,
                 min_comp_free_energy=0.5,
                 sa_energy_weights=None):
        """

        :param hx_rate_fpath:
        :param pH:
        :param temp:
        :param pdb_fpath:
        :param dg_intpol_fpath:
        :param comp_dg_fpath:
        :param nter:
        :param cter:
        :param hbond_length:
        :param hbond_angle:
        :param min_free_energy:
        :param net_charge_corr:
        :param r_constant:
        :param min_comp_free_energy:
        :param sa_energy_weights:
        """

        # save some meta data
        self.metadata = DgInputMetadata(nterm=nter,
                                        cterm=cter,
                                        pH=pH,
                                        temp=temp,
                                        hbond_length=hbond_length,
                                        hbond_angle=hbond_angle,
                                        min_free_energy=min_free_energy,
                                        net_charge_corr=net_charge_corr,
                                        r_constant=r_constant,
                                        min_comp_free_energy=min_comp_free_energy)

        # self.metadata.hxrate_fpath = hx_rate_fpath
        self.metadata.pdb_fpath = pdb_fpath
        self.metadata.deltag_interpol_fpath = dg_intpol_fpath
        self.metadata.comp_dg_fpath = comp_dg_fpath

        # set deltag_interpol
        with open(dg_intpol_fpath, 'rb') as pk_file:
            deltag_interpol = pickle.load(pk_file, encoding='latin1')

        if sa_energy_weights is None:
            self.weights = EnergyWeights()
            self.metadata.weights = 'Generated from EnergyWeights'
        else:
            self.weights = EnergyWeights(pair_energy=sa_energy_weights['pair_energy'],
                                         full_burial=sa_energy_weights['full_burial'],
                                         hbond_burial=sa_energy_weights['hbond_burial'],
                                         hbond_rank=sa_energy_weights['hbond_rank'],
                                         dist_to_sec_struct=sa_energy_weights['dist_to_sec_struct'],
                                         dist_to_nonpolar_res=sa_energy_weights['dist_to_nonpolar_res'],
                                         top_energy_std=sa_energy_weights['top_energy_std'],
                                         comp_deltag_rmse=sa_energy_weights['comp_deltag_rmse'])
            self.metadata.weights = 'Generated from user provided dictionary'

        # get residue info
        structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_fpath, pdb_fpath)
        model = structure[0]
        residues_ob_list = [x for x in model.get_residues()]
        self.raw_sequence = ''.join([Polypeptide.three_to_one(aa.get_resname()) for aa in residues_ob_list])
        self.full_sequence = nter + self.raw_sequence + cter

        self.seq_num = np.zeros(len(self.full_sequence))
        for num in range(len(self.full_sequence)):
            if num > len(nter) - 1:
                if num < len(self.full_sequence) - len(cter):
                    self.seq_num[num] = num - len(nter) + 1

        # get proline res indices
        self.proline_res_ind_raw_seq = get_proline_res_ind(sequence=self.raw_sequence)
        self.proline_res_ind_full_seq = get_proline_res_ind(sequence=self.full_sequence)

        # get intrinsic rate
        self.intrinsic_rates = calculate_intrinsic_exchange_rates_suggie(sequence_str=self.full_sequence,
                                                                         Temperature=temp,
                                                                         pH=pH)
        # first two rates to be zero
        self.intrinsic_rates[1] = 0.0

        self.sec_struct_string_raw_seq = get_dssp_sec_structure_assignment(pdb_filepath=pdb_fpath)
        self.sec_struct_string_full_seq = 'L' * len(nter) + self.sec_struct_string_raw_seq + 'L' * len(cter)

        self.dist_to_sec_struct = distance_to_secondary_structure(
            secondary_structure_string=self.sec_struct_string_full_seq)

        self.dist_to_nonpolar_res = distance_to_nonpolar_residues(sequence=self.full_sequence)

        self.hx_allowed_raw_seq_bool = allowed_res(sequence=self.raw_sequence)
        self.hx_allowed_full_seq_bool = allowed_res(sequence=self.raw_sequence, nterm=nter, cterm=cter)

        self.len_hx_allowed = sum(self.hx_allowed_full_seq_bool)

        self.hx_allowed_seq_index = []
        for ind, hx_allowed_bool in enumerate(self.hx_allowed_full_seq_bool):
            if hx_allowed_bool:
                self.hx_allowed_seq_index.append(ind)
        self.hx_allowed_seq_index = np.array(self.hx_allowed_seq_index)

        # hx rates
        # self.hx_rates = load_hx_rates(fpath=hx_rate_fpath)
        self.hx_rates = hxrate_obj

        # active hx measured and intrinsic rates
        self.active_hx_rates = np.exp(sorted(self.hx_rates.mean[0:self.len_hx_allowed]))
        self.active_hx_intrinsic_rates = self.intrinsic_rates[self.hx_allowed_full_seq_bool]

        # get hbond info
        hbond_data = gen_hbond_data(pdb_fpath=pdb_fpath,
                                    nterm=nter,
                                    cterm=cter,
                                    hbond_length=hbond_length,
                                    hbond_angle=hbond_angle)

        # create free energy grid
        self.net_charge = calc_net_charge_at_ph(protein_sequence=self.full_sequence, pH=pH)

        self.free_energy_grid = gen_free_energy_grid_from_hx_rates(active_hx_rates=self.active_hx_rates,
                                                                   active_intrinsic_hx_rate=
                                                                   self.active_hx_intrinsic_rates,
                                                                   min_free_energy=min_free_energy,
                                                                   temp=temp,
                                                                   r_constant=r_constant,
                                                                   net_charge_corr=net_charge_corr,
                                                                   net_charge=self.net_charge)

        hh_pairlist = []
        for pairlist in hbond_data.pairlist:
            aa1, aa2, hh_dist = pairlist[0], pairlist[1], pairlist[2]
            if hbond_data.hbond_bool_list[aa1] and hbond_data.hbond_bool_list[aa2]:
                hh_pairlist.append([aa1, aa2, hh_dist])
        hh_pairlist = np.array(hh_pairlist)

        if len(hh_pairlist) > 0:

            self.valid = True

            hx_active_index = np.array([sum(self.hx_allowed_full_seq_bool[0:x]) for x in range(len(self.hx_allowed_full_seq_bool))],
                                       dtype=object)
            for num in range(len(hx_active_index)-1):
                if hx_active_index[num+1] == hx_active_index[num]:
                    hx_active_index[num] = None

            self.active_aa1, self.active_aa2 = [], []
            for aa1, aa2 in zip(hh_pairlist[:, 0], hh_pairlist[:, 1]):
                if hx_active_index[int(aa1)] is not None:
                    if hx_active_index[int(aa2)] is not None:
                        self.active_aa1.append(hx_active_index[int(aa1)])
                        self.active_aa2.append(hx_active_index[int(aa2)])
            self.active_aa1, self.active_aa2 = np.array(self.active_aa1), np.array(self.active_aa2)

            self.pair_energies_multigrid = gen_pair_energies_multigrid(len_hx_allowed=self.len_hx_allowed,
                                                                       active_aa1_array=self.active_aa1,
                                                                       active_aa2_array=self.active_aa2,
                                                                       hh_dist_array=hh_pairlist[:, 2],
                                                                       free_energy_grid=self.free_energy_grid,
                                                                       interpol_func=deltag_interpol)

            self.states = list(range(sum(self.hx_allowed_full_seq_bool)))

            self.num_hbonds = sum(hbond_data.hbond_bool_num_list)
            self.num_no_hbonds = len(self.states) - self.num_hbonds
            hbond_rank_best_ranks = np.arange(self.num_no_hbonds, len(self.states))
            self.hbond_rank_best = sum(hbond_rank_best_ranks)
            self.hbond_rank_worst = sum(np.arange(self.num_hbonds))

            self.hbond_bool_num_list = hbond_data.hbond_bool_num_list

            self.hbond_burial = hbond_data.hbond_burial

            if comp_dg_fpath is not None:
                self.deltag_comp_arr = get_comp_dg_array(comp_dg_fpath=comp_dg_fpath,
                                                         sequence=self.raw_sequence,
                                                         min_comp_free_energy=min_comp_free_energy)
            else:
                self.deltag_comp_arr = None

        else:

            self.valid = False


@dataclass
class DGMapOut(object):
    dg_map_valid: bool = False
    protein_name: str = None
    protein_raw_sequence: str = None
    protein_full_sequence: str = None
    nterm_add: str = None
    cterm_add: str = None
    ph: float = None
    temp: float = None
    state: np.ndarray = None
    anneal_data: object = None
    res_num: np.ndarray = None
    dg_array: np.ndarray = None


class DeltaGMapping(Annealer):
    """
    Annealing class to map delta G values for each residue
    """

    def __init__(self, dg_input, traj_fpath=None):

        self.traj_fpath = traj_fpath

        state = dg_input.states

        # randomize state here
        random.shuffle(state)

        self.dg_input = dg_input
        self.weights = self.dg_input.weights

        self.step_num_ = 0

        self.anneal_data = None

        super(DeltaGMapping, self).__init__(state)

    def move(self):
        """
        randomly move and swap two states
        :return:
        """
        a = random.randint(0, len(self.state) - 1)
        b = random.randint(0, len(self.state) - 1)
        self.state[a], self.state[b] = self.state[b], self.state[a]

    def energy(self):
        """
        Calculate the energy ...
        :return: energy
        """

        self.anneal_data = DGAnnealData()
        self.anneal_data.pair_energy = self.weights.pair_energy * self.pair_energy()
        self.anneal_data.full_burial_corr = self.weights.full_burial * self.full_burial_corr()
        self.anneal_data.hbond_burial_corr = self.weights.hbond_burial * self.hbond_burial_corr()
        self.anneal_data.hbond_rank_factor = self.weights.hbond_rank * self.hbond_rank_fact()
        self.anneal_data.distance_to_nonpolar_res_corr = self.weights.dist_to_nonpolar_res * self.distance_to_nonpolar_res_corr()
        self.anneal_data.distance_to_sec_struct_corr = self.weights.dist_to_sec_struct * self.distance_to_sec_struct_corr()
        self.anneal_data.top_stdev = self.weights.top_energy_std * self.top_stdev_energy()

        if self.dg_input.deltag_comp_arr is not None:
            self.anneal_data.comp_deltaG_rmse_term = self.weights.comp_deltag_rmse * self.comp_data_mse()
        else:
            self.anneal_data.comp_deltaG_rmse_term = np.nan

        self.anneal_data.opt_val = np.nansum([x for x in self.anneal_data.__dict__.values()])
        # anneal_data.opt_val = np.nansum([x for x in anneal_data.__dict__.values()])

        # store anneal data after update interval
        if self.traj_fpath is not None:
            self.step_num_ += 1
            if self.step_num_ % self.updates == 0:
                write_traj_file(traj_fpath=self.traj_fpath,
                                init_header=False,
                                step_num=self.step_num_,
                                anneal_data=self.anneal_data)

        return self.anneal_data.opt_val

    def map_energy(self):
        current_free_energy = np.array([self.dg_input.free_energy_grid[ind][eng] for ind, eng in enumerate(self.state)])
        out_array = np.zeros(len(self.dg_input.full_sequence))
        out_array[self.dg_input.hx_allowed_seq_index] = current_free_energy
        # out_array[sa_input['active_hx_rates'] == 0] = np.nan
        if len(self.dg_input.metadata.cterm) > 0:
            out_array[-len(self.dg_input.metadata.cterm):] = -2
        if len(self.dg_input.metadata.nterm) > 0:
            out_array[0:len(self.dg_input.metadata.nterm)] = -2
        return out_array

    def pair_energy(self):
        pair_e = []
        for aa1, aa2 in zip(self.dg_input.active_aa1, self.dg_input.active_aa2):
            # print('indices')
            # print(aa1, aa2, self.state[aa1], self.state[aa2])
            eng = self.dg_input.pair_energies_multigrid[aa1, aa2, self.state[aa1], self.state[aa2]]
            # print('eng')
            # print(eng)
            pair_e.append(eng)
        avg_pair_e = np.average(pair_e)
        return avg_pair_e

    def comp_data_mse(self):
        """
        compare energy dataset from other technique and add to optimization score
        :return: mean square error between
        """
        current_free_energy = np.array([self.dg_input.free_energy_grid[ind][eng] for ind, eng in enumerate(self.state)])

        delta_g_rmse = mean_squared_error(self.dg_input.deltag_comp_arr, current_free_energy)

        return delta_g_rmse

    def full_burial_corr(self):
        """
        calculate the correlation coeffs between an array of free energy grid and burial distribution. Return the Cij
        :return:
        """
        current_free_energy = np.array([self.dg_input.free_energy_grid[ind][eng] for ind, eng in enumerate(self.state)])

        corr = 0
        if min(current_free_energy) == max(current_free_energy):
            corr += 0
        else:
            current_free_energy_fill_min = np.maximum(current_free_energy, 0)
            burial = self.dg_input.hbond_burial[self.dg_input.hx_allowed_seq_index]
            burial_fill_min = np.maximum(burial, 60)
            corr_coef = np.corrcoef(current_free_energy_fill_min, burial_fill_min)
            corr += -corr_coef[0][1]
        return corr

    def hbond_burial_corr(self):
        """
        calculate the correlation coeffs between the burial distance where hbond and hx is allowed and free energy grid.
        Return the Cij
        :return:
        """
        current_free_energy = np.array([self.dg_input.free_energy_grid[ind][eng] for ind, eng in enumerate(self.state)])

        hbond_free_eng = current_free_energy[self.dg_input.hbond_bool_num_list[self.dg_input.hx_allowed_seq_index] == 1]

        corr = 0
        if min(hbond_free_eng) == max(hbond_free_eng):
            corr += 0
        else:
            hbond_free_eng_fill_min = np.maximum(hbond_free_eng, 0)
            hbond_burial = self.dg_input.hbond_burial[
                (self.dg_input.hbond_bool_num_list == 1) & (self.dg_input.hx_allowed_full_seq_bool) == True]
            hbond_burial_fill_min = np.maximum(hbond_burial, 60)
            corr_coef = np.corrcoef(hbond_free_eng_fill_min, hbond_burial_fill_min)
            corr += -corr_coef[0][1]
        return corr

    def distance_to_nonpolar_res_corr(self):

        current_free_energy = np.array([self.dg_input.free_energy_grid[ind][eng] for ind, eng in enumerate(self.state)])

        corr = 0
        if min(current_free_energy) == max(current_free_energy):
            corr += 0
        else:
            current_free_energy_fill_min = np.maximum(current_free_energy, 0)
            dist_non_polar_res_arr = self.dg_input.dist_to_nonpolar_res[self.dg_input.hx_allowed_full_seq_bool]
            corr_coef = np.corrcoef(current_free_energy_fill_min, dist_non_polar_res_arr)
            corr += -corr_coef[0][1]
        return corr

    def distance_to_sec_struct_corr(self):

        current_free_energy = np.array([self.dg_input.free_energy_grid[ind][eng] for ind, eng in enumerate(self.state)])

        corr = 0
        if min(current_free_energy) == max(current_free_energy):
            corr += 0
        else:
            current_free_energy_fill_min = np.maximum(current_free_energy, 0)
            dist_sec_struct_arr = self.dg_input.dist_to_sec_struct[self.dg_input.hx_allowed_full_seq_bool]
            corr_coef = np.corrcoef(current_free_energy_fill_min, dist_sec_struct_arr)
            corr += -corr_coef[0][1]
        return corr

    def hbond_rank_fact(self):

        current_free_energy = np.array([self.dg_input.free_energy_grid[ind][eng] for ind, eng in enumerate(self.state)])

        sort_index = np.argsort(current_free_energy)

        ranks = np.empty(len(current_free_energy), int)
        ranks[sort_index] = np.arange(len(current_free_energy))

        hbond_ranks = ranks[self.dg_input.hbond_bool_num_list[self.dg_input.hx_allowed_seq_index] == 1]

        hbond_rank_factor = (float(sum(hbond_ranks)) - self.dg_input.hbond_rank_worst) / (
                self.dg_input.hbond_rank_best - self.dg_input.hbond_rank_worst)

        return -hbond_rank_factor

    def top_stdev_energy(self):
        """
        top 5% of energy stdev
        :param sa_input_dict:
        :return:
        """

        current_free_energy = np.array([self.dg_input.free_energy_grid[ind][eng] for ind, eng in enumerate(self.state)])

        sorted_curr_free_energy = sorted(current_free_energy)
        num_k_percent = np.ceil(0.05 * len(current_free_energy))
        sorted_curr_free_energy_top = sorted_curr_free_energy[-int(num_k_percent):]

        stdev_full = np.std(sorted_curr_free_energy_top)

        return stdev_full


def get_comp_dg_array(comp_dg_fpath, sequence, min_comp_free_energy=0.5):

    compdf = pd.read_csv(comp_dg_fpath)
    resnum, deltag_val = compdf.iloc[:, 0], compdf.iloc[:, 1]
    res_index = np.subtract(resnum, 1)

    proline_res_ind = get_proline_res_ind(sequence=sequence)

    if proline_res_ind is None:
        comp_res_index_final = res_index
        comp_deltag_final = deltag_val
    else:
        comp_res_index_final = []
        comp_deltag_final = []
        for ind, (res_ind, dg_val) in enumerate(zip(res_index, deltag_val)):
            if res_ind not in proline_res_ind:
                comp_res_index_final.append(res_ind)
                comp_deltag_final.append(dg_val)
        comp_res_index_final = np.array(comp_res_index_final)
        comp_deltag_final = np.array(comp_res_index_final)

    deltag_comp = np.zeros(len(sequence))
    deltag_comp[comp_res_index_final] = comp_deltag_final

    allowed_res_ind = allowed_res(sequence=sequence)

    deltag_comp_allowed = deltag_comp[allowed_res_ind]
    deltag_comp_allowed[deltag_comp_allowed == 0] = min_comp_free_energy

    return deltag_comp_allowed


def gen_pair_energies_multigrid(len_hx_allowed,
                                active_aa1_array,
                                active_aa2_array,
                                hh_dist_array,
                                free_energy_grid,
                                interpol_func):

    pair_energies_multigrid = np.zeros((len_hx_allowed, len_hx_allowed, len_hx_allowed, len_hx_allowed))
    for ind, (aa1, aa2, hhdist) in enumerate(zip(active_aa1_array, active_aa2_array, hh_dist_array)):
        for num1 in range(len_hx_allowed):
            for num2 in range(len_hx_allowed):
                fe1 = free_energy_grid[aa1, num1]
                fe2 = free_energy_grid[aa2, num2]
                fe_dist = abs(fe1 - fe2)
                # if type(fe_dist) == np.ndarray:
                #     sort_fedist = sorted(fe_dist[0])
                #     fe_dist[0] = sort_fedist
                #     # print('heha')
                # else:
                #     pass
                # # print('hhdist: %s\nfedist: %s' % (hhdist, fe_dist))
                interp = interpol_func(hhdist, fe_dist)
                pair_e = -np.log(interp)
                # print('paire: %s\nSuccess' % pair_e)
                pair_energies_multigrid[aa1, aa2, num1, num2] = pair_e

    return pair_energies_multigrid


def calc_net_charge_at_ph(protein_sequence, pH=6.15):
    """
    calculate net charge at certain pH
    :param protein_sequence:
    :param pH:
    :return:
    """
    analyses_seq = ProteinAnalysis(protein_sequence)
    net_charge = analyses_seq.charge_at_pH(pH)
    return net_charge


def corr_fe_with_net_charge(fe, net_charge):
    """

    :param fe:
    :param net_charge:
    :return:
    """
    fe_corr = fe - (net_charge * -0.12)
    return fe_corr


def calc_free_energy_from_hx_rates(hx_instrinsic_rate, hx_meas_rate, temperature=295, r_constant=1.987204e-3):
    """
    calculate free energy of opening (change in protein energy landscape upon opening or unfolding)
    :param hx_instrinsic_rate: hx intrinsic rate
    :param hx_meas_rate: hx measured rate
    :param r_constant: 1.987204e-3 kcal mol-1 K-1
    :param temperature: temperature in Kelvin
    :return: free energy of opening = - R T ln (k_op)
    """
    k_op = hx_meas_rate/(hx_instrinsic_rate - hx_meas_rate)
    free_energy = -r_constant * temperature * np.log(k_op)
    return free_energy


def gen_free_energy_grid_from_hx_rates(active_hx_rates, active_intrinsic_hx_rate, min_free_energy=0,
                                       temp=295, r_constant=1.9872036e-3, net_charge_corr=False, net_charge=None):
    """
    generate free energy grid as N x M were there are N hx rate list and M intrisnsic rate list
    :param active_hx_rates: hx fitted rate (measured) for residues that are allowed to exchange
    :param active_intrinsic_hx_rate: intrinsic hx rate based on sequence for residues that are allowed to exchange
    :return: free energy matrix
    """
    free_energy_arr = np.zeros((len(active_hx_rates), len(active_intrinsic_hx_rate)))

    for ind1, intrinsic_hx_rate in enumerate(active_intrinsic_hx_rate):
        for ind2, fit_hx_rate in enumerate(active_hx_rates):
            if fit_hx_rate >= intrinsic_hx_rate:
                free_energy_arr[ind1, ind2] = min_free_energy
            else:
                fe = calc_free_energy_from_hx_rates(intrinsic_hx_rate, fit_hx_rate, temp,
                                                    r_constant)
                if net_charge_corr:
                    fe_corr = corr_fe_with_net_charge(fe, net_charge)
                    if fe_corr == np.inf:
                        fe_corr = min_free_energy
                    free_energy_arr[ind1, ind2] = fe_corr
                else:
                    if fe == np.inf:
                        fe = min_free_energy
                    free_energy_arr[ind1, ind2] = fe
                # free_energy_arr[ind1, ind2] = calc_free_energy_from_hx_rates(intrinsic_hx_rate, fit_hx_rate, temp,
                #                                                          r_constant)

    return free_energy_arr


def allowed_res(sequence, nterm='', cterm=''):
    """
    gen allowed mask based on sequence
    :param sequence: protein sequence
    :param nterm: nterm addition
    :param cterm: cterm addition
    :return:
    """

    allowed_bool_list = [False if x == 'P' else True for x in sequence]
    if len(nterm) > 0:
        allowed_bool_list = [False for _ in range(len(nterm))] + allowed_bool_list
    if len(cterm) > 0:
        allowed_bool_list = allowed_bool_list + [False for _ in range(len(cterm))]

    allowed_bool_array = np.array(allowed_bool_list)

    # set the first two residues to False
    allowed_bool_array[:2] = False

    return allowed_bool_array


def get_proline_res_ind(sequence):
    """
    get proline res ind
    :param sequence:
    :return:
    """
    aa_list = [str(x) for x in sequence]
    if 'P' in aa_list:
        proline_res_ind = []
        for ind, aa in enumerate(aa_list):
            if aa == 'P':
                proline_res_ind.append(ind)
    else:
        proline_res_ind = None

    return proline_res_ind


def create_aa_dict():
    """
    Create a dictionary containing amino acid three letter, one letter, polar, non polar residue classes
    :return: dictionary
    """
    aa_dict = dict()
    aa_dict['aa_code'] = dict()
    aa_dict['aa_code']['three_letter'] = ['Arg', 'His', 'Lys', 'Asp', 'Glu', 'Ser', 'Thr', 'Asp', 'Gln', 'Cys', 'Gly',
                                          'Pro', 'Ala', 'Ile', 'Leu', 'Met', 'Phe', 'Trp', 'Tyr', 'Val']
    aa_dict['aa_code']['one_letter'] = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'I', 'L', 'M',
                                        'F', 'W', 'Y', 'V']
    # polar amino acids contain polar uncharged and charged side chains. Only in one letter code
    aa_dict['polar_amino_acids'] = ['D', 'E', 'G', 'H', 'K', 'N', 'P', 'Q', 'R', 'S', 'T']
    aa_dict['non_polar_amino_acids'] = ['A', 'C', 'F', 'I', 'L', 'M', 'V', 'W', 'Y']
    return aa_dict


def smooth_dist_cutoff(x):
    """
    exponential function for smooth continuous function
    :param x: variable
    :return: function output
    """
    return 1. / (1. + (np.exp(x - 10)))


def distance_to_nonpolar_residues(sequence):
    """
    Mapping of each residue's nearest distance to a non polar residue. 0 is 1 residue away, -1 is 2 residues away.
    1 indicates that it itself is a non polar residue. positive distances mean distance from nearest polar residues.
    It searches the nearest non polar residues using a window search. Creates a window by iteratively adding an amino
    acid on the either side (n and c terminal ends). When a non polar residue is included in the window, the algorithm
    detects it and finds how many index far it is. That index is then stored in the output array.
    :param sequence: string of input sequence
    :return: array of distances [0,-1,-2,-3,0,1,2,1,0] etc as described above
    """
    aa_dict = create_aa_dict()
    dist_arr = []
    for index in range(len(sequence)):
        current_residue = sequence[index]
        # if the current residue is a polar type, then start generating windows to search for the nearest non polar
        # amino acid and find the distance in terms of index distance
        if current_residue in aa_dict['polar_amino_acids']:
            for num in range(1, len(sequence)):
                window_bounds = [max([0, index - num]), min([len(sequence), index+num+1])]
                window = sequence[window_bounds[0]:window_bounds[1]]
                for non_polar_res in aa_dict['non_polar_amino_acids']:
                    if non_polar_res in window:
                        dist_arr.append(- num + 1)
                        break
                if non_polar_res in window:
                    break
            if num == len(sequence) - 1:
                dist_arr.append(- num + 1)

        else:
            for num in range(1, len(sequence)):
                window_bounds = [max([0, index - num]), min([len(sequence), index+num+1])]
                window = sequence[window_bounds[0]:window_bounds[1]]
                for polar_res in aa_dict['polar_amino_acids']:
                    if polar_res in window:
                        dist_arr.append(num)
                        break
                if polar_res in window:
                    break
            if num == len(sequence) - 1:
                dist_arr.append(num)
    dist_arr = np.array(dist_arr)
    return dist_arr


def get_dssp_sec_structure_assignment(pdb_filepath):
    """
    assigns secondary structure to the given structure object from biopython. Uses the biopython dssp module. It
    requires you to have working dssp executable
    # conda install -c salilab dssp
    Biopython will automatically find dssp and use it. In this case, you don't have to specify the dssp executable path
    # change loop assignments (S, C, T, -) to one letter assignment 'L'
    :param pdb_filepath: filepath for the pdb file
    :param dssp_exec_path: dssp executable path
    :return: dssp assignment. Dssp dict contains the following list for each dict item.
            aa, -> aa one letter code
            ss, -> secondary structure assignment
            acc, -> number of H2O molecules in contact with this residue * 10 or residue water exposed surface in Angstrom **2
            phi, -> backbone torsion angle
            psi, -> backbone torsion angle
            dssp_index, -> dssp index starting with 1
            NH_O_1_relidx, -> Index of O residue in which it forms H bond with NH group
            NH_O_1_energy, -> H bond energy of the previous scenario
            O_NH_1_relidx, -> Index of NH residue in which it forms H bond with OH group
            O_NH_1_energy, -> H bond energy of the previous scenario
            NH_O_2_relidx, -> Index of O residue in which it forms H bond with the NH group
            NH_O_2_energy, -> H bond energy of the previous scenario
            O_NH_2_relidx, -> Index of NH group in which it forms H bond with the OH group
            O_NH_2_energy -> H bond energy of the previous scenario
    """
    dssp_output = dssp_dict_from_pdb_file(pdb_filepath)
    sec_struct_string = ''
    for ind, key_values in enumerate(dssp_output[0].values()):
        sec_struct = key_values[1]
        if sec_struct == 'S' or sec_struct == 'T' or sec_struct == 'C' or sec_struct == '-':
            sec_struct = 'L'
        sec_struct_string += sec_struct
    return sec_struct_string


def distance_to_secondary_structure(secondary_structure_string):
    """
    Mapping of loop to the nearest helix/strand in -ve index distance and helix/strand to the nearest loop in +ve index
    distance
    :param secondary_structure_string: string of secondary structure assignments that contain L, H, and E
    :return: distance_array [0,-1,-2,-3,0,1,2,1,0] etc as described above
    """
    dist_arr = []
    for index in range(len(secondary_structure_string)):
        current_sec_struct = secondary_structure_string[index]
        # if current secondary structure is a loop 'L'
        if current_sec_struct == 'L':
            for num in range(1, len(secondary_structure_string)):
                window_bounds = [max([0, index - num]), min([len(secondary_structure_string), index + num + 1])]
                window = secondary_structure_string[window_bounds[0]:window_bounds[1]]
                # if helix 'H' or beta strand 'E' in the window, then append the distance to the dist_arr
                if 'E' in window or 'H' in window:
                    dist_arr.append(- num + 1)
                    break
            if num == len(secondary_structure_string) - 1:
                dist_arr.append(- num + 1)

        # if current secondary structure is either 'H' or beta strand 'E'
        else:
            for num in range(1, len(secondary_structure_string)):
                window_bounds = [max([0, index - num]), min([len(secondary_structure_string), index + num + 1])]
                window = secondary_structure_string[window_bounds[0]:window_bounds[1]]
                if 'L' in window:
                    dist_arr.append(num)
                    break
            if num == len(secondary_structure_string) - 1:
                dist_arr.append(num)

    dist_arr = np.array(dist_arr)
    return dist_arr


def load_hx_rates(fpath):

    df = pd.read_csv(fpath)
    hxrate = HXRate(mean=df['rate_mean'].values,
                    median=df['rate_median'].values,
                    std=df['rate_std'].values,
                    ci_5=df['rate_5%'].values,
                    ci_95=df['rate_95%'].values)
    return hxrate


def write_traj_file(traj_fpath, init_header=False, step_num=None, anneal_data=None):

    if init_header:
        header = 'step_num,dist_to_sec_struct_corr,dist_to_nonpolar_res_corr,full_burial_corr,hbond_burial_corr,hbond_rank_factor,pair_energy,top_stdev,comp_deltaG_rmse,opt_val\n'
        with open(traj_fpath, 'w') as trajfile:
            trajfile.write(header)
    else:
        line = '{},{},{},{},{},{},{},{},{},{}\n'.format(step_num,
                                                        anneal_data.distance_to_sec_struct_corr,
                                                        anneal_data.distance_to_nonpolar_res_corr,
                                                        anneal_data.full_burial_corr,
                                                        anneal_data.hbond_burial_corr,
                                                        anneal_data.hbond_rank_factor,
                                                        anneal_data.pair_energy,
                                                        anneal_data.top_stdev,
                                                        anneal_data.comp_deltaG_rmse_term,
                                                        anneal_data.opt_val)
        with open(traj_fpath, 'a') as trajfile:
            trajfile.write(line)
            trajfile.close()


def hbond_depth(reisude_list, nterm, cterm, hbond_list, all_hbond_list):

    residues = reisude_list
    hbond_depth_list = [1 if x in all_hbond_list else 0 for x in range(len(residues))]

    for num in range(1, len(residues) - 1):
        if hbond_depth_list[num - 1: num + 2] == [1, 0, 1]:
            hbond_depth_list[num] = 1
            all_hbond_list.append(num)

    num2 = 1
    while hbond_depth_list.count(num2) > 1:
        for index in range(num2, len(residues) - num2):
            if hbond_depth_list[index - 1: index + 2] in [[num2, num2, num2], [num2 + 1, num2, num2]]:
                hbond_depth_list[index] += 1
        num2 += 1

    hbond_depth_arr = np.array(hbond_depth_list) * np.array(
        [1 if x in hbond_list else 0 for x in range(len(residues))])

    # include the nterm and cterm additional residues
    hbond_depth_arr = np.array([0 for x in nterm] + list(hbond_depth_arr) + [0 for x in cterm])

    return hbond_depth_arr, all_hbond_list


def gen_hbond_data(pdb_fpath, nterm='', cterm='', hbond_length=2.7, hbond_angle=120):

    # load the pdb file to a model
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_fpath, pdb_fpath)
    model = structure[0]

    residues = [x for x in model.get_residues()]

    hbond_list = []
    hbond_dists = [None for _ in residues]
    hbond_angles = [None for _ in residues]
    hbond_partners = [None for _ in residues]
    all_hbonds = []
    sc_hbonds = []
    pairlist = []
    atoms_in_radius = [0 for x in nterm]
    for index, res in enumerate(residues):
        atoms_in_radius.append(
            sum([smooth_dist_cutoff(np.linalg.norm(res['N'].coord.astype(float) - atm.coord.astype(float)))
                 for atm in model.get_atoms() if 'H' not in atm.get_name()]))
        if 'H' in res:
            for ind, res2 in enumerate(residues):
                if ind == index:
                    continue
                new_hbond_dist = np.linalg.norm(res['H'].coord.astype(float) - res2['O'].coord.astype(float))
                new_hbond_angle = np.rad2deg(min([Bio.PDB.calc_angle(res['N'].get_vector(),
                                                                     res['H'].get_vector(),
                                                                     res2['O'].get_vector()),
                                                  Bio.PDB.calc_angle(res['H'].get_vector(),
                                                                     res2['O'].get_vector(),
                                                                     res2['C'].get_vector())]))
                if (new_hbond_dist < hbond_length) and (new_hbond_angle > hbond_angle):
                    hbond_list.append(index)
                    all_hbonds.append(index)
                    all_hbonds.append(ind)
                    hbond_partners[index] = ind + len(nterm)
                    hbond_dists[index] = new_hbond_dist
                    hbond_angles[index] = new_hbond_angle

                for o_atom in [x for x in res2 if x.name[0] == 'O' and x.name != 'O']:
                    new_hbond_dist = np.linalg.norm(res['H'].coord.astype(float) - o_atom.coord.astype(float))
                    if new_hbond_dist < hbond_length:
                        sc_hbonds.append(index)

                if index < ind:
                    if 'H' in res and 'H' in res2:
                        h_distance = np.linalg.norm(res['H'].coord.astype(float) - res2['H'].coord.astype(float))
                        if h_distance < 13:
                            pairlist.append((index + len(nterm), ind + len(nterm), h_distance))

    atoms_in_radius += [0 for _ in cterm]

    hbond_list = sorted(set(hbond_list))
    all_hbonds = sorted(set(all_hbonds))

    hbond_dists = np.array([None for _ in nterm] + list(hbond_dists) + [None for _ in cterm])
    hbond_angles = np.array([None for _ in nterm] + list(hbond_angles) + [None for _ in cterm])
    hbond_partners = np.array([None for _ in nterm] + list(hbond_partners) + [None for _ in cterm])

    hbond_depth_array, all_hbond_list = hbond_depth(reisude_list=residues,
                                                    nterm=nterm,
                                                    cterm=cterm,
                                                    hbond_list=hbond_list,
                                                    all_hbond_list=all_hbonds)

    hbond_bool_list = np.array([False if x == 0 else True for x in hbond_depth_array])

    hbond_bool_num_list = []
    for item in hbond_bool_list:
        if item:
            hbond_bool_num_list.append(1)
        else:
            hbond_bool_num_list.append(0)
    hbond_bool_num_list = np.array(hbond_bool_num_list)

    hbond_data = HBond(hbond_list=hbond_list,
                       all_hbond_list=all_hbond_list,
                       hbond_depth_array=hbond_depth_array,
                       hbond_bool_list=hbond_bool_list,
                       hbond_bool_num_list=hbond_bool_num_list,
                       hbond_burial=np.array(atoms_in_radius),
                       pairlist=pairlist,
                       hbond_dist_list=hbond_dists,
                       hbond_angle_list=hbond_angles,
                       hbond_partner_list=hbond_partners)

    return hbond_data


def dg_opt_from_input_obj(dg_input,
                          dg_length_mins,
                          dg_update_interval,
                          traj_fpath):

    # initialize dg anneal
    deltag_anneal = DeltaGMapping(dg_input=dg_input)

    # set auto schedule
    auto_schedule = deltag_anneal.auto(minutes=dg_length_mins)
    auto_schedule['updates'] = dg_update_interval

    # set up writing traj file

    deltag_anneal.traj_fpath = traj_fpath
    if traj_fpath is not None:
        write_traj_file(traj_fpath=traj_fpath, init_header=True)

    # run anneal
    deltag_anneal.set_schedule(auto_schedule)
    deltag_anneal.copy_strategy = "slice"
    deltag_anneal.anneal()

    return deltag_anneal


def get_the_best_anneal_data(dg_anneal):

    # set the state to the best state
    dg_anneal.state = dg_anneal.best_state

    # recalculate the energy with the updated state
    dg_anneal.energy()

    # get the anneal data
    anneal_data = dg_anneal.anneal_data
    return anneal_data


def map_energy_to_res(free_energy_grid, anneal_state, protein_sequence, hx_allowed_index):

    free_energy_arr = np.array([free_energy_grid[ind][eng] for ind, eng in enumerate(anneal_state)])
    out_energy_array = np.empty(len(protein_sequence))
    out_energy_array[:] = np.nan
    out_energy_array[hx_allowed_index] = free_energy_arr

    return out_energy_array


def write_dg_data_to_csv(seq_num, dg_array, output_path):

    header = 'res_num,raw_res_num,dg\n'
    data = ''

    for ind, (res_num, dg) in enumerate(zip(seq_num, dg_array)):
        data += '{},{},{}\n'.format(ind+1, int(res_num), dg)

    with open(output_path, 'w') as outfile:
        outfile.write(header+data)
        outfile.close()


def write_anneal_data_to_csv(anneal_data, output_path):

    header = ','.join([x for x in anneal_data.__dict__.keys()])
    header += '\n'

    data = ','.join([str(x) for x in anneal_data.__dict__.values()])
    data += '\n'

    with open(output_path, 'w') as outfile:
        outfile.write(header+data)
        outfile.close()


def plot_dg_dist_(mapped_energy, res_num, min_free_energy, output_path):

    res_num_re = res_num[res_num > 0]
    mapped_energy_re = mapped_energy[res_num > 0]
    mapped_energy_re[mapped_energy_re == min_free_energy] = np.nan

    num_rows = 2
    num_columns = 1
    fig_size = (7, 3 * num_rows)
    fig = plt.figure(figsize=fig_size)
    gs = fig.add_gridspec(nrows=num_rows, ncols=num_columns)

    ax = fig.add_subplot(gs[0, 0])
    plt.plot(res_num_re, mapped_energy_re, marker='o', ls='-', markerfacecolor='red', markeredgecolor='black', color='red')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xticks(range(int(res_num_re[0]), int(res_num_re[-1])+1, 2))
    plt.xlabel('Residue number')
    plt.ylabel('dG (kcal/mol)')
    plt.grid(alpha=0.25)
    ax.tick_params(length=3, pad=3)

    map_energy_no_nan = mapped_energy_re[~np.isnan(mapped_energy_re)]
    sort_map_energy_no_nan = np.sort(map_energy_no_nan)[::-1]

    ax2 = fig.add_subplot(gs[1, 0])
    plt.plot(sort_map_energy_no_nan, marker='o', ls='-', markerfacecolor='red', markeredgecolor='black', color='red')
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xticks(range(0, len(sort_map_energy_no_nan)+1, 2))
    plt.xlabel('Residue (ranked from most to least stable)')
    plt.ylabel('dG (kcal/mol)')
    plt.grid(alpha=0.25)
    ax2.tick_params(length=3, pad=3)

    plt.subplots_adjust(hspace=0.2, wspace=0.1, top=0.96)

    title_1 = 'dG distribution'

    plt.suptitle(title_1)

    plt.savefig(output_path, bbox_inches="tight")
    plt.close()


def dg_mapping(hx_rate_fpath,
               pdb_fpath,
               dg_intpol_fpath,
               pH,
               temp,
               comp_dg_fpath=None,
               nterm=None,
               cterm=None,
               min_free_energy=-10,
               net_charge_corr=True,
               min_comp_free_energy=0.5,
               sa_energy_weights=None,
               dg_length_mins=3.0,
               dg_update_interval=200,
               traj_fpath=None,
               anneal_data_output=None,
               dg_csv_output=None,
               dg_data_output=None,
               dg_plot_path=None,
               return_flag=False):

    if nterm is None or nterm == 'None':
        nterm = ''
    if cterm is None or cterm == 'None':
        cterm = ''

    hxrate_obj = load_hx_rates(hx_rate_fpath)
    dg_input = DgInput(hxrate_obj=hxrate_obj,
                       pdb_fpath=pdb_fpath,
                       dg_intpol_fpath=dg_intpol_fpath,
                       pH=pH,
                       temp=temp,
                       comp_dg_fpath=comp_dg_fpath,
                       nter=nterm,
                       cter=cterm,
                       min_free_energy=min_free_energy,
                       net_charge_corr=net_charge_corr,
                       min_comp_free_energy=min_comp_free_energy,
                       sa_energy_weights=sa_energy_weights)
    dg_input.metadata.hxrate_fpath = hx_rate_fpath

    # init dg output object
    dg_output = DGMapOut(dg_map_valid=dg_input.valid,
                         protein_name=os.path.split(pdb_fpath)[-1],
                         protein_raw_sequence=dg_input.raw_sequence,
                         protein_full_sequence=dg_input.full_sequence,
                         nterm_add=nterm,
                         cterm_add=cterm,
                         ph=pH,
                         temp=temp,
                         res_num=dg_input.seq_num)

    if dg_input.valid:

        deltag_anneal = dg_opt_from_input_obj(dg_input=dg_input,
                                              dg_length_mins=dg_length_mins,
                                              dg_update_interval=dg_update_interval,
                                              traj_fpath=traj_fpath)

        dg_output.state = deltag_anneal.best_state
        dg_output.anneal_data = get_the_best_anneal_data(dg_anneal=deltag_anneal)
        dg_output.dg_array = map_energy_to_res(free_energy_grid=dg_input.free_energy_grid,
                                               anneal_state=deltag_anneal.best_state,
                                               protein_sequence=dg_input.full_sequence,
                                               hx_allowed_index=dg_input.hx_allowed_seq_index)

    else:

        print("Alert: dG calc not possible due to the structure ... Saving emtpy files")

        dg_array = np.array([np.nan for _ in range(len(dg_input.seq_num))])
        dg_output.dg_array = dg_array

        dg_output.anneal_data = DGAnnealData(pair_energy=np.nan,
                                             full_burial_corr=np.nan,
                                             hbond_burial_corr=np.nan,
                                             hbond_rank_factor=np.nan,
                                             distance_to_nonpolar_res_corr=np.nan,
                                             distance_to_sec_struct_corr=np.nan,
                                             top_stdev=np.nan,
                                             comp_deltaG_rmse_term=np.nan)
        dg_output.anneal_data.opt_val = np.nan

        if traj_fpath is not None:
            write_traj_file(traj_fpath=traj_fpath,
                            init_header=True)

    if dg_csv_output is not None:
        write_dg_data_to_csv(seq_num=dg_output.res_num, dg_array=dg_output.dg_array, output_path=dg_csv_output)

    if anneal_data_output is not None:
        write_anneal_data_to_csv(anneal_data=dg_output.anneal_data, output_path=anneal_data_output)

    if dg_plot_path is not None:
        plot_dg_dist_(mapped_energy=dg_output.dg_array,
                      res_num=dg_output.res_num,
                      min_free_energy=dg_input.metadata.min_free_energy,
                      output_path=dg_plot_path)

    # convert dg output to dict
    dg_output.anneal_data = vars(dg_output.anneal_data)
    dg_output = vars(dg_output)

    if dg_data_output is not None:
        with open(dg_data_output, 'wb') as outfile:
            pickle.dump(dg_output, outfile)

    if return_flag:
        return dg_output


def gen_parser_args():

    parser_ = argparse.ArgumentParser(prog='DGOptimize', description='DGOptimize')
    parser_.add_argument('-hx', '--hxrate_pipeline', default='/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/merge_lib15_ph6_ph7/output/EEHEE_rd4_0642.pdb_9.047_EEHEE_rd4_0642.pdb_9.11541/EEHEE_rd4_0642.pdb_9.047_EEHEE_rd4_0642.pdb_9.11541_hx_rate.csv', help='hx_rate filepath .csv')
    parser_.add_argument('-pdb', '--pdbfpath', default="/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/EEHEE_rd4_0642.pdb", help='pdb filepath .pdb')
    parser_.add_argument('-dip', '--dgintpath', default="/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/refit_new_SA_weights/newrect.pickle", help='dg interpolation filepath .pickle')
    parser_.add_argument('-cf', '--compdgfpath', default=None, help='comp dg file path .csv')
    parser_.add_argument('-p', '--ph', default=6.15, type=float, help='pH')
    parser_.add_argument('-t', '--temp', default=295, type=float, help='temperature in Kelvin')
    parser_.add_argument('-nt', '--nterm', default=None, type=str or None, help='n terminal addition')
    parser_.add_argument('-ct', '--cterm', default=None, type=str or None, help='c terminal addition')
    parser_.add_argument('--netcharge', default=True, action=argparse.BooleanOptionalAction, help='correct fe using net charge from protein sequence')
    parser_.add_argument('-au', '--annealupdate', default=100, type=int, help='Anneal update interval')
    parser_.add_argument('-at', '--annealtime', default=2.0, type=float, help='Anneal time')
    parser_.add_argument('-tf', '--trajfile', default=None, help='trajectory file path .csv')
    parser_.add_argument('-af', '--annealfile', default=None, help='anneal data output file path .csv')
    parser_.add_argument('-df', '--dgfile', default=None, help='dg data output file path .csv')
    parser_.add_argument('-dpf', '--dgpicklefile', default=None, help='dg data output pickle file path .pickle')
    parser_.add_argument('-po', '--plotoutput', default=None, help='dg data plot output path .pdf')

    return parser_


def run_anneal_from_parser():

    parser_ = gen_parser_args()
    options = parser_.parse_args()

    dg_mapping(hx_rate_fpath=options.hxrate,
               pdb_fpath=options.pdbfpath,
               dg_intpol_fpath=options.dgintpath,
               pH=options.ph,
               temp=options.temp,
               comp_dg_fpath=options.compdgfpath,
               nterm=options.nterm,
               cterm=options.cterm,
               net_charge_corr=options.netcharge,
               dg_length_mins=options.annealtime,
               dg_update_interval=options.annealupdate,
               traj_fpath=options.trajfile,
               anneal_data_output=options.annealfile,
               dg_csv_output=options.dgfile,
               dg_data_output=options.dgpicklefile,
               dg_plot_path=options.plotoutput,
               return_flag=False)


if __name__ == '__main__':

    run_anneal_from_parser()

    # pdb_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/test/pdbs/PDB1H8K.pdb'
    #
    # hx_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/test/rates/PDB1H8K_14.70437_PDB1H8K_14.73122/PDB1H8K_14.70437_PDB1H8K_14.73122_hx_rate.csv'
    #
    # dg_interpol_fpath = "../../config/newrect.pickle"
    #
    # dg_mapping(hx_rate_fpath=hx_fpath,
    #            pdb_fpath=pdb_fpath,
    #            dg_intpol_fpath=dg_interpol_fpath,
    #            pH=6.00,
    #            temp=295,
    #            comp_dg_fpath=None,
    #            min_free_energy=-10,
    #            net_charge_corr=True,
    #            min_comp_free_energy=0.5,
    #            sa_energy_weights=None,
    #            dg_length_mins=0.5,
    #            dg_update_interval=100,
    #            traj_fpath=hx_fpath+'_anneal_traj.csv',
    #            anneal_data_output=hx_fpath+'_anneal_data.csv',
    #            dg_csv_output=hx_fpath+'_dg_data.csv',
    #            dg_data_output=hx_fpath+'_dg_data.pickle',
    #            dg_plot_path=hx_fpath+'_dg_data.pdf',
    #            return_flag=False)
    #
    # stop

    # dg_input = DgInput(hx_rate_fpath=hx_fpath,
    #                    pH=6.00,
    #                    temp=295,
    #                    pdb_fpath=pdb_fpath,
    #                    dg_intpol_fpath=dg_interpol_fpath)



    # hbond_obj = gen_hbond_data(pdb_fpath=pdb_fpath)
    #
    # hh_pairlist = []
    # for pairlist in hbond_obj.pairlist:
    #     aa1, aa2, hh_dist = pairlist[0], pairlist[1], pairlist[2]
    #     if hbond_obj.hbond_bool_list[aa1] and hbond_obj.hbond_bool_list[aa2]:
    #         hh_pairlist.append([aa1, aa2, hh_dist])
    # hh_pairlist = np.array(hh_pairlist)

    # print('heho')
    #

    ############################################################

    # merge_names = ['PDB1H8K_14.70437_PDB1H8K_14.73122',
    #                'PDB3NGP_12.98532_PDB3NGP_13.56409',
    #                'PDB1W6X_22.07855_PDB1W6X_22.0786',
    #                'PDB5FWB_16.10001_PDB5FWB_16.10955',
    #                'PDB2J5Y_13.15926_PDB2J5Y_13.15708',
    #                'PDB2KY5_9.238_PDB2KY5_9.23759']
    #
    # # 4th index is the working case
    #
    # num = 4
    #
    # hx_rate_fpath = '/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/test/rates/' + merge_names[num] + '/' + merge_names[num] + '_hx_rate.csv'
    # dg_interpol_fpath = "../../config/newrect.pickle"
    # pdb_fpath = "/Users/smd4193/OneDrive - Northwestern University/hx_ratefit_gabe/hxratefit_new/bayes_opt/test/pdbs/" + merge_names[num].split('_')[0] + ".pdb"
    #
    # dg_mapping(hx_rate_fpath=hx_rate_fpath,
    #            pdb_fpath=pdb_fpath,
    #            dg_intpol_fpath=dg_interpol_fpath,
    #            pH=6.00,
    #            temp=295,
    #            comp_dg_fpath=None,
    #            min_free_energy=-10,
    #            net_charge_corr=True,
    #            min_comp_free_energy=0.5,
    #            sa_energy_weights=None,
    #            dg_length_mins=0.5,
    #            dg_update_interval=100,
    #            traj_fpath=hx_rate_fpath+'_anneal_traj_3.csv',
    #            anneal_data_output=hx_rate_fpath+'_anneal_data_3.csv',
    #            dg_csv_output=hx_rate_fpath+'_dg_data_3.csv',
    #            dg_data_output=hx_rate_fpath+'_dg_data_3.pickle',
    #            dg_plot_path=hx_rate_fpath+'_dg_data_3.pdf',
    #            return_flag=False)
