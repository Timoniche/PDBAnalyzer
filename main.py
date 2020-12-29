import math

import numpy as np

from PDBAnalyzer import PDBAnalyzer, heatmap
import logging
import cooler
import matplotlib.pyplot as plt

from PDBUtils import hist, sort_by_x, log_xy, pearson_hic_dist

from sklearn.linear_model import LinearRegression


def main():
    analyze_chr1()


def analyze_chr1(SOFT_NAME='3DMAX'):
    BUCKETS_CNT = 100

    logging.basicConfig(filename='logs/analyzer.log', filemode='w', level=logging.INFO)
    ratio, dist = chr1_scaled_dist_3dmax()
    hic = chr1_hic(filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')

    factor = 0.8
    logging.info(
        f'Pearson coef: {pearson_hic_dist(np.array(hic)[:499][:499], (np.array(dist) * ratio)[:499][:499], factor)}')

    dist_xs = []
    hic_ys = []
    for i in range(499):
        for j in range(i + 1, 499):
            dist_xs.append(dist[i][j] * ratio)
            hic_ys.append(hic[i][j])

    shrinked_ys = hist(dist_xs, hic_ys, BUCKETS_CNT)

    plt.title(SOFT_NAME + ' x density')
    plt.xlabel('3d dist mm')
    plt.ylabel('count of points')
    plt.hist(dist_xs)
    plt.show()

    plt.title(SOFT_NAME + ' y(x) = hic(dist)')
    plt.xlabel('3d dist mm')
    plt.ylabel('hic')
    xs, ys = sort_by_x(dist_xs, shrinked_ys)
    cut_size = -1
    plt.plot(xs[:cut_size], ys[:cut_size], 'o', markersize=5)
    plt.show()

    xs, ys = log_xy(xs, ys)
    plt.title(SOFT_NAME + ' ln( y(x) ) = ln( hic(dist) )')
    plt.xlabel('ln( 3d dist mm )')
    plt.ylabel('ln( hic )')

    xs_linear = np.array(xs).reshape((-1, 1))
    ys_linear = np.array(ys)
    model: LinearRegression
    model = LinearRegression().fit(xs_linear, ys_linear)
    r_sq = model.score(xs_linear, ys_linear)
    logging.info(SOFT_NAME + f' linear regression r_2:\n{r_sq}\n')
    logging.info(SOFT_NAME + f' b0 (intercept):\n{model.intercept_}\n')
    logging.info(SOFT_NAME + f' k (slope):\n{model.coef_}\n')

    xs_ans = np.linspace(min(xs), max(xs), 100)
    ys_ans = model.intercept_ + model.coef_[0] * xs_ans

    plt.plot(xs[:cut_size], ys[:cut_size], 'o', markersize=5)
    plt.plot(xs_ans, ys_ans, color='orange', linewidth=3)
    plt.show()


def chr1_scaled_dist_3dmax():
    file_name = '3DMax_chr1_real.pdb'
    CURVE_LEN_MM = 85
    BINS_CNT = 499

    analyzer = PDBAnalyzer('pdb_files/' + file_name, BINS_CNT)
    dist_mat = analyzer.count_curve_dist_matrix()
    heatmap(dist_mat, plot_title='distance matrix from 3d pdb curve')
    curve_length = PDBAnalyzer.pdb_curve_length(dist_mat, BINS_CNT)
    logging.info(f'curve length (in pdb model scales)\n{curve_length}\n')

    ratio = PDBAnalyzer.ratio_to_real_size(dist_mat, BINS_CNT, CURVE_LEN_MM)
    logging.info(f'ratio to real length [real_len / curve_len] is {ratio}\n')

    PDBAnalyzer.build_scaled_pdb(
        PDBAnalyzer.extract_traj('pdb_files/' + file_name),
        ratio,
        file_name
    )
    return ratio, dist_mat


def chr1_hic(filepath):
    c = cooler.Cooler(filepath)
    bin_size = c.info['bin-size']
    count_bins_chr1 = math.ceil(c.chromsizes[0] / bin_size)
    return c.matrix(balance=False)[:count_bins_chr1, :count_bins_chr1]


if __name__ == '__main__':
    main()
