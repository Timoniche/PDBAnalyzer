import math

from PDBAnalyzer import PDBAnalyzer, heatmap
import logging
import cooler
import matplotlib.pyplot as plt


def shrink(arr, bucket):
    arr_shrinked = []
    for i in range(0, len(arr), bucket):
        v_avg = 0
        v_cnt = 0
        (num_ret, _) = arr[i]
        for b in range(bucket):
            if i + b < len(arr):
                v_cnt += 1
                (num, v) = arr[i + b]
                v_avg += v
        arr_shrinked.append((num_ret, v_avg / v_cnt))
    return arr_shrinked


def log_xy(xs, ys):
    xys = list(zip(xs, ys))
    xys = list(filter(lambda e: e[0] > 1e-9 and e[1] > 1e-9, xys))
    dist_xs = list(map(lambda e: e[0], xys))
    hic_ys = list(map(lambda e: e[1], xys))
    dist_xs = list(map(lambda e: math.log(e), dist_xs))
    hic_ys = list(map(lambda e: math.log(e), hic_ys))
    return dist_xs, hic_ys

def main():
    logging.basicConfig(filename='logs/analyzer.log', filemode='w', level=logging.INFO)
    ratio, dist = chr1_scaled_dist_3dmax()
    hic = chr1_hic(filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')

    dist_xs = []
    hic_ys = []
    for i in range(499):
        for j in range(i + 1, 499):
            dist_xs.append(dist[i][j] * ratio)
            hic_ys.append(hic[i][j])


    dist_xs_num = list(enumerate(dist_xs))
    hic_ys_num = list(enumerate(hic_ys))

    dist_xs_num = sorted(dist_xs_num, key=lambda e: e[1])
    hic_ys_num = sorted(hic_ys_num, key=lambda e: e[1])

    dist_xs_num = shrink(dist_xs_num, 50)
    hic_ys_num = shrink(hic_ys_num, 50)

    dist_xs_num = sorted(dist_xs_num, key=lambda e: e[0])
    hic_ys_num = sorted(hic_ys_num, key=lambda e: e[0])

    dist_xs_num = list(map(lambda e: e[1], dist_xs_num))
    hic_ys_num = list(map(lambda e: e[1], hic_ys_num))


    # xys = list(zip(dist_xs, hic_ys))
    # xys = list(filter(lambda e: e[0] > 1e-9 and e[1] > 35000, xys))
    # dist_xs = list(map(lambda e: e[0], xys))
    # hic_ys = list(map(lambda e: e[1], xys))
    # dist_xs = list(map(lambda e: math.log(e), dist_xs))
    # hic_ys = list(map(lambda e: math.log(e), hic_ys))
    a, b = log_xy(dist_xs_num, hic_ys_num)


    plt.title('y(x) = hic(dist)')
    plt.xlabel('3d dist mm')
    plt.ylabel('hic')
    # plt.plot(dist_xs_num, hic_ys_num, 'o')
    # plt.plot(dist_xs, hic_ys, 'o')
    plt.plot(a, b, 'o')
    plt.show()


def chr1_scaled_dist_3dmax():
    file_name = '3DMax_chr1_real.pdb'
    curve_len_mm = 85
    bins_cnt = 499

    analyzer = PDBAnalyzer('pdb_files/' + file_name, bins_cnt)
    logging.info(f'pdb xyz\'s:\n{analyzer.xyz}\n')
    dist_mat = analyzer.count_curve_dist_matrix()
    logging.info(f'distance matrix from 3d pdb curve:\n{dist_mat}\n')
    heatmap(dist_mat, plot_title='distance matrix from 3d pdb curve')
    curve_length = PDBAnalyzer.pdb_curve_length(dist_mat, bins_cnt)
    logging.info(f'curve length (in pdb model scales)\n{curve_length}\n')

    ratio = PDBAnalyzer.ratio_to_real_size(dist_mat, bins_cnt, curve_len_mm)
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
