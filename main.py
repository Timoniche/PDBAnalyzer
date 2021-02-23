from HiCExtractor import chr_i_hic_to_txt_adj, chr_i_hic_to_npy, chr_i_hic
from PDBAnalyzer import analyze_chr_i, PDBAnalyzer
import numpy as np

def main():
    # analyze_chr_i(1, file_name='3DMax_chr2.pdb')

    gen_all_3dmax_adj()


def gen_all_simba_npy():
    for i in range(21):
        chr_i_hic_to_npy(i=i, filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')


def gen_all_3dmax_adj():
    for i in range(21):
        chr_i_hic_to_txt_adj(filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool', i=i)


def store_distmat(pdbpath, bins_cnt, chr):
    analyzer = PDBAnalyzer(pdbpath, bins_cnt)
    mat = analyzer.count_curve_dist_matrix()
    np.save('dists_npy' + chr + '.npy', mat)


if __name__ == '__main__':
    chr_i = 10
    hic, bins_cnt, bp = chr_i_hic(filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool', i=chr_i - 1)
    store_distmat(pdbpath='pdb_files/Simba3D/Simba3d_chr10.pdb', bins_cnt=bins_cnt, chr=f'chr{chr_i}')
