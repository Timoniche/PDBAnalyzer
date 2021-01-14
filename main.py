from HiCExtractor import chr_i_hic_to_txt_adj, chr_i_hic_to_npy
from PDBAnalyzer import analyze_chr_i


def main():
    # analyze_chr_i(1, file_name='3DMax_chr2.pdb')
    # chr_i_hic_to_txt_adj(filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool', i=5)
    chr_i_hic_to_npy(i=1, filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')


if __name__ == '__main__':
    main()
