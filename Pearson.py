from math import pow, sqrt
import numpy as np


def calc_average(x):
    return div_by_zero(sum(x), len(x))


def div_by_zero(x, y):
    if y == 0:
        return 0
    else:
        return x / y


def pearson(xs: np.ndarray, ys: np.ndarray):
    x = xs.flatten()
    y = ys.flatten()
    avX = calc_average(x)
    avY = calc_average(y)
    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    for i in range(len(x)):
        sum1 += (x[i] - avX) * (y[i] - avY)
        sum2 += pow((x[i] - avX), 2)
        sum3 += pow((y[i] - avY), 2)
    r = div_by_zero(sum1, (sqrt(sum2) * sqrt(sum3)))
    return r
