import math

import cooler
import numpy as np


def chr_i_hic(filepath, i):
    c = cooler.Cooler(filepath)
    bin_size = c.info['bin-size']
    count_bins_chr_i = math.ceil(c.chromsizes[i] / bin_size)
    (left_bin_id, right_bin_id) = c.extent('chr' + str(i + 1))
    return (c.matrix(balance=False)[left_bin_id:right_bin_id, left_bin_id:right_bin_id],
            count_bins_chr_i,
            c.chromsizes[i])


# simba input
def chr_i_hic_to_npy(filepath, i):
    c = cooler.Cooler(filepath)
    (left_bin_id, right_bin_id) = c.extent('chr' + str(i + 1))
    mat = c.matrix(balance=False)[left_bin_id:right_bin_id, left_bin_id:right_bin_id]
    np.save('npy_hic_inputs/chr' + str(i + 1) + '.npy', mat)


# 3dmax input
def chr_i_hic_to_txt_adj(filepath, i):
    hic_file = open('gen_adj_hics/chr' + str(i + 1) + '.txt', 'w')
    c = cooler.Cooler(filepath)
    bin_size = c.info['bin-size']
    count_bins_chr_i = math.ceil(c.chromsizes[i] / bin_size)
    (left_bin_id, right_bin_id) = c.extent('chr' + str(i + 1))
    mat = c.matrix(balance=False)[left_bin_id:right_bin_id, left_bin_id:right_bin_id]
    for i in range(count_bins_chr_i):
        for j in range(count_bins_chr_i):
            hic_file.write(f'{i + 1} {j + 1} {mat[i][j]}\n')
    hic_file.close()
