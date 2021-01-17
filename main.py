from HiCExtractor import chr_i_hic_to_txt_adj, chr_i_hic_to_npy
from PDBAnalyzer import analyze_chr_i


def main():
    # analyze_chr_i(1, file_name='3DMax_chr2.pdb')

    gen_all_3dmax_adj()


def gen_all_simba_npy():
    for i in range(21):
        chr_i_hic_to_npy(i=i, filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')


def gen_all_3dmax_adj():
    for i in range(21):
        chr_i_hic_to_txt_adj(filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool', i=i)


if __name__ == '__main__':
    main()
