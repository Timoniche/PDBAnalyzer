import matplotlib.pyplot as plt
import numpy as np


def plot_density(soft_name, file_name, xs, buckets_cnt):
    plt.title(soft_name + ' x density ' + file_name)
    plt.xlabel('3d dist mm')
    plt.ylabel('count of points')
    plt.hist(xs, bins=buckets_cnt)
    plt.show()


def plot_hic_from_dist(soft_name, file_name, xs, ys):
    plt.title(soft_name + ' y(x) = hic(dist) ' + file_name)
    plt.xlabel('3d dist mm')
    plt.ylabel('hic')
    cut_size = -1
    plt.plot(xs[:cut_size], ys[:cut_size], 'o', markersize=5)
    plt.show()


def plot_ln_hic_from_dist(file_name, model, xs, ys, r_sq, buckets_cnt):
    plt.title(
        file_name + ' ln(hic(dist)) ' + ' r_sq: {:.2f} '.format(r_sq) + ' k^-1: {:.2f}'.format(1.0 / model.coef_[0]))
    plt.xlabel('ln( 3d dist mm )')
    plt.ylabel('ln( hic )')
    cut_size = -1
    xs_ans = np.linspace(min(xs), max(xs), 100)
    ys_ans = model.intercept_ + model.coef_[0] * xs_ans
    plt.plot(xs[:cut_size], ys[:cut_size], 'o', markersize=5)
    plt.plot(xs_ans, ys_ans, color='orange', linewidth=3)

    counts, bins = np.histogram(xs, bins=buckets_cnt)
    scale = max(ys) / max(counts)
    plt.hist(bins[:-1], bins, weights=scale*counts, alpha=0.2)
    plt.show()
