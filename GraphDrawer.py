import random

from HiCExtractor import chr_i_hic
from PDBAnalyzer import chr_i_scaled_dist
from PDBUtils import bp_to_mm


SOFT_NAME='Simba3d'
file_name = f'Simba3d_chr1.pdb'
write_to_filepath = 'rand_breakpoints.txt'

def main():
    hic, bins_cnt, bp = chr_i_hic(filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool', i=0)
    len_mm = bp_to_mm(bp)
    ratio, dist = chr_i_scaled_dist(bins_cnt=bins_cnt,
                                    soft_name=SOFT_NAME,
                                    file_name=file_name,
                                    curve_len_mm=len_mm,
                                    dist_plot_on=False)
    dist_real = dist * ratio

    n_breakpnts = 15
    bins = random.sample(range(0, 499), n_breakpnts)

    file = open(write_to_filepath, "w")
    for i in range(n_breakpnts):
        for j in range(i + 1, n_breakpnts):
            near_f = 1 / (dist_real[bins[i]][bins[j]] ** 2)
            file.write(f'{bins[i]} {bins[j]} {near_f}\n')
    file.close()


if __name__ == '__main__':
    main()