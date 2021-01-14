from HiCExtractor import chr_i_hic_to_txt_adj
from PDBAnalyzer import analyze_chr_i


def main():
    analyze_chr_i(0)
    chr_i_hic_to_txt_adj(filepath='chr/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool', i=1)


if __name__ == '__main__':
    main()
