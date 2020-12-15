import math

from PDBAnalyzer import PDBAnalyzer, heatmap
import logging
import cooler


def main():
    logging.basicConfig(filename='logs/analyzer.log', filemode='w', level=logging.INFO)
    example_3dmax()
    histogram(filepath='data/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')


def example_3dmax():
    file_name = '3DMax_chr1.pdb'
    # file_name = 'Simba3d_chr1.pdb'
    bins_cnt = 1000
    curve_len_mm = 170

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


def histogram(filepath):
    c = cooler.Cooler(filepath)
    bin_size = 500000
    count_bins_chr1 = math.ceil(c.chromsizes[0] / bin_size)
    contact_matrix = c.matrix(balance=False)[:count_bins_chr1, :count_bins_chr1]
    print(contact_matrix)


if __name__ == '__main__':
    main()
