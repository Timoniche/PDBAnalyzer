import math

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt


def dist3d(x1, y1, z1, x2, y2, z2):
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)


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
        self.xyz = np.array([(0.0, 0.0, 0.0) for _ in range(self.bins_cnt)])
        for i in range(min(len(smaller_xyz_why_dont_know), bins_cnt)):
            self.xyz[i] = smaller_xyz_why_dont_know[i]
        self.dist_matrix = np.array([])

    def count_curve_dist_matrix(self):
        self.dist_matrix = np.array([[0.0 for _ in range(self.bins_cnt)] for _ in range(self.bins_cnt)])
        for i in range(self.bins_cnt):
            for j in range(i + 1, self.bins_cnt):
                fst = self.xyz[i]
                snd = self.xyz[j]
                dist = dist3d(fst[0], fst[1], fst[2], snd[0], snd[1], snd[2])
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
