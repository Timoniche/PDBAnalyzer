import math
import numpy as np
from sklearn.linear_model import LinearRegression


def mm_to_bp(mm):
    return int(2939495.094814116 * float(mm))


def bp_to_mm(bp):
    return float(bp) / 2939495.094814116


def linear_regression(xs, ys):
    xs_linear = np.array(xs).reshape((-1, 1))
    ys_linear = np.array(ys)
    model: LinearRegression
    model = LinearRegression().fit(xs_linear, ys_linear)
    r_sq = model.score(xs_linear, ys_linear)
    return model, r_sq


def pearson_hic_dist(hic: np.ndarray, dist: np.ndarray, factor):
    a, b = hic.shape
    dist_theory = np.zeros(hic.shape)
    for i in range(a):
        for j in range(i + 1, b):
            if hic[i][j] != 0:
                dist_theory[i][j] = dist_theory[j][i] = 1 / (hic[i][j] ** factor)
    return np.corrcoef(dist.flatten(), dist_theory.flatten())[0, 1]


def sort_by_x(xs, ys):
    xys = list(zip(xs, ys))
    xys_sorted = sorted(xys, key=lambda xy: xy[0])
    xs_new, ys_new = zip(*xys_sorted)
    return xs_new, ys_new


def log_xy(xs, ys, ignore_eps=1e-9):
    xys = list(zip(xs, ys))
    xys = list(filter(lambda e: e[0] > ignore_eps and e[1] > ignore_eps, xys))
    xys = list(map(lambda e: (math.log(e[0]), math.log(e[1])), xys))
    xs_log, ys_log = zip(*xys)
    return xs_log, ys_log


def shrink_ys_to_hist(xs, ys, buckets_cnt=10):
    max_x = max(xs)
    buckets = [(0, 0) for _ in range(buckets_cnt)]
    for i in range(len(xs)):
        bucket_idx = int((xs[i] / max_x) * (buckets_cnt - 1))
        (value_pred, cnt_pred) = buckets[bucket_idx]
        buckets[bucket_idx] = (value_pred + ys[i], cnt_pred + 1)

    shrinked_ys = []
    for i in range(len(xs)):
        bucket_idx = int((xs[i] / max_x) * (buckets_cnt - 1))
        (v, cnt) = buckets[bucket_idx]
        shrinked_ys.append(v / cnt)

    return shrinked_ys
