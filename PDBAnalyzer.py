import logging

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

from HiCExtractor import chr_i_hic
from LogListener import log_linear_regression
from PDBUtils import bp_to_mm, pearson_hic_dist, shrink_ys_to_hist, sort_by_x, log_xy, linear_regression
from Plotter import plot_density, plot_hic_from_dist, plot_ln_hic_from_dist


def dist3d(fst3d: np.ndarray, snd3d: np.ndarray):
    return np.linalg.norm(fst3d - snd3d)


def heatmap(arr, plot_title):
    plt.title(plot_title)
    plt.imshow(arr, cmap='hot', interpolation='nearest')
    plt.show()


class PDBAnalyzer:
    prefix_path = 'pdb_files/'
    prefix_scaled_path = 'scaled_pdb_files/'

    def __init__(self, pdb_path, bins_cnt) -> None:
        self.pdb_path = pdb_path
        self.traj = PDBAnalyzer.extract_traj(self.pdb_path)
        smaller_xyz_why_dont_know = self.traj.xyz[0]
        self.bins_cnt = bins_cnt
        self.xyz = np.array([[0.0, 0.0, 0.0] for _ in range(self.bins_cnt)])
        for i in range(min(len(smaller_xyz_why_dont_know), bins_cnt)):
            self.xyz[i] = smaller_xyz_why_dont_know[i]
        self.dist_matrix = np.array([])

    def count_curve_dist_matrix(self):
        self.dist_matrix = np.array([[0.0 for _ in range(self.bins_cnt)] for _ in range(self.bins_cnt)])
        for i in range(self.bins_cnt):
            for j in range(i + 1, self.bins_cnt):
                fst = self.xyz[i]
                snd = self.xyz[j]
                dist = dist3d(np.array(fst), np.array(snd))
                self.dist_matrix[i][j] = self.dist_matrix[j][i] = dist
        return self.dist_matrix

    @staticmethod
    def extract_traj(path):
        return md.load_pdb(path)

    @staticmethod
    def pdb_curve_length(dists, cnt_bins):
        curve_length = 0
        bin_cur = 0
        bin_next = 1
        for i in range(cnt_bins - 1):
            curve_length += dists[bin_cur][bin_next]
            bin_cur += 1
            bin_next += 1
        return curve_length

    @staticmethod
    def ratio_to_real_size(dists, cnt_bins, real_size):
        curve_length = PDBAnalyzer.pdb_curve_length(dists, cnt_bins)
        return real_size / curve_length

    @staticmethod
    def build_scaled_pdb(traj: md.Trajectory, ratio, pdb_name):
        f_rescale = np.vectorize(lambda t: t * ratio)
        traj.xyz[0] = f_rescale(traj.xyz[0])
        scaled_pdb_name = PDBAnalyzer.prefix_scaled_path + pdb_name
        traj.save_pdb(filename=scaled_pdb_name)


def analyze_chr_i(i, file_name, SOFT_NAME='3DMax', buckets_cnt=1000, known_factor=-1):
    logging.basicConfig(filename='logs/analyzer.log', filemode='w', level=logging.INFO)
    logging.info('\n' + SOFT_NAME + f' chr: {i + 1}\n')

    hic, bins_cnt, bp = chr_i_hic(filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool', i=i)

    len_mm = bp_to_mm(bp)
    ratio, dist = chr_i_scaled_dist(bins_cnt=bins_cnt,
                                    soft_name=SOFT_NAME,
                                    file_name=file_name,
                                    curve_len_mm=len_mm)

    if known_factor != -1:
        pearson = pearson_hic_dist(
            hic=np.array(hic)[:bins_cnt][:bins_cnt],
            dist=(np.array(dist) * ratio)[:bins_cnt][:bins_cnt],
            factor=known_factor
        )
        logging.info(f'Pearson coef: {pearson}')

    dist_xs = []
    hic_ys = []
    for i in range(bins_cnt):
        for j in range(i + 1, bins_cnt):
            dist_xs.append(dist[i][j] * ratio)
            hic_ys.append(hic[i][j])

    shrinked_ys = shrink_ys_to_hist(dist_xs, hic_ys, buckets_cnt)

    plot_density(soft_name=SOFT_NAME, file_name=file_name, xs=dist_xs, buckets_cnt=buckets_cnt)

    xs, ys = sort_by_x(dist_xs, shrinked_ys)
    plot_hic_from_dist(soft_name=SOFT_NAME, file_name=file_name, xs=xs, ys=ys)

    xs, ys = log_xy(xs, ys)
    model, r_sq = linear_regression(xs, ys)
    plot_ln_hic_from_dist(file_name=file_name, model=model, xs=xs, ys=ys, r_sq=r_sq, buckets_cnt=buckets_cnt)

    log_linear_regression(model=model, r_sq=r_sq, soft_name=SOFT_NAME)


def chr_i_scaled_dist(bins_cnt, soft_name, file_name, curve_len_mm):
    analyzer = PDBAnalyzer('pdb_files/' + soft_name + '/' + file_name, bins_cnt)
    dist_mat = analyzer.count_curve_dist_matrix()
    curve_length = PDBAnalyzer.pdb_curve_length(dist_mat, bins_cnt)
    logging.info(f'curve length (in pdb model scales)\n{curve_length}\n')

    ratio = PDBAnalyzer.ratio_to_real_size(dist_mat, bins_cnt, curve_len_mm)
    logging.info(f'ratio to real length [real_len / curve_len] is {ratio}\n')

    PDBAnalyzer.build_scaled_pdb(
        PDBAnalyzer.extract_traj('pdb_files/' + soft_name + '/' + file_name),
        ratio,
        file_name
    )
    heatmap(dist_mat * ratio, plot_title='distance matrix from 3d pdb curve ' + file_name)
    return ratio, dist_mat
