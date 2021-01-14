import math

import cooler


def chr_i_hic(filepath, i):
    c = cooler.Cooler(filepath)
    bin_size = c.info['bin-size']
    count_bins_chr_i = math.ceil(c.chromsizes[i] / bin_size)
    (left_bin_id, right_bin_id) = c.extent('chr' + str(i + 1))
    return (c.matrix(balance=False)[left_bin_id:right_bin_id, left_bin_id:right_bin_id],
            count_bins_chr_i,
            c.chromsizes[i])
