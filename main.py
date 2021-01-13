import math

import numpy as np

from LogListener import log_linear_regression
from PDBAnalyzer import PDBAnalyzer, heatmap
import logging
import cooler

from PDBUtils import shrink_ys_to_hist, sort_by_x, log_xy, pearson_hic_dist, bp_to_mm, linear_regression

from Plotter import plot_density, plot_hic_from_dist, plot_ln_hic_from_dist


def main():
    analyze_chr_i(0)


def analyze_chr_i(i, SOFT_NAME='3DMax', buckets_cnt=1000, known_factor=-1):
    logging.basicConfig(filename='logs/analyzer.log', filemode='w', level=logging.INFO)
    logging.info('\n' + SOFT_NAME + f' chr: {i + 1}\n')

    hic, bins_cnt, bp = chr_i_hic(filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool', i=i)

    file_name = '3DMax_chr1.pdb'
    len_mm = bp_to_mm(bp)
    ratio, dist = chr1_scaled_dist_3dmax(bins_cnt=bins_cnt,
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


def chr1_scaled_dist_3dmax(bins_cnt, soft_name, file_name, curve_len_mm):
    analyzer = PDBAnalyzer('pdb_files/' + soft_name + '/' + file_name, bins_cnt)
    dist_mat = analyzer.count_curve_dist_matrix()
    heatmap(dist_mat, plot_title='distance matrix from 3d pdb curve ' + file_name)
    curve_length = PDBAnalyzer.pdb_curve_length(dist_mat, bins_cnt)
    logging.info(f'curve length (in pdb model scales)\n{curve_length}\n')

    ratio = PDBAnalyzer.ratio_to_real_size(dist_mat, bins_cnt, curve_len_mm)
    logging.info(f'ratio to real length [real_len / curve_len] is {ratio}\n')

    PDBAnalyzer.build_scaled_pdb(
        PDBAnalyzer.extract_traj('pdb_files/' + soft_name + '/' + file_name),
        ratio,
        file_name
    )
    return ratio, dist_mat


def chr_i_hic(filepath, i):
    c = cooler.Cooler(filepath)
    bin_size = c.info['bin-size']
    count_bins_chr_i = math.ceil(c.chromsizes[i] / bin_size)
    (left_bin_id, right_bin_id) = c.extent('chr' + str(i + 1))
    return (c.matrix(balance=False)[left_bin_id:right_bin_id, left_bin_id:right_bin_id],
            count_bins_chr_i,
            c.chromsizes[i])


if __name__ == '__main__':
    main()
