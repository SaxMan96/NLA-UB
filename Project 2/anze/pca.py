import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
np.set_printoptions(suppress=True)


def pca_cov(xm, n=None):
    if n is None:
        n = xm.shape[0]
    ym = xm.T / np.sqrt(n - 1)
    um, s, vmt = np.linalg.svd(ym, full_matrices=False)
    s_pow = np.power(s, 2)
    # cov = np.multiply(vmt, s_pow).dot(vmt.T)
    return None, s_pow, vmt, s


def pca_corr(xm, std, n=None, transpose=False):
    if n is None:
        n = xm.shape[0]
    matrix = xm
    if transpose:
        matrix = xm.T
    matrix = matrix / std
    return pca_cov(matrix, n)


def calc_portion_or_variance(s2):
    # accum = np.cumsum(s2)
    return s2 / np.sum(s2)


def get_example1():
    with open('example.dat') as f:
        n = 16
        m = 4
        xm = np.zeros((n, m))
        for idx, line in enumerate(f):
            if idx >= n:
                break
            xm[idx, :] = [int(st) for st in line.strip().split(' ')]
    return xm


def get_RCsGoff():
    data = pd.read_csv("RCsGoff.csv")
    data = data.iloc[:, 1:].to_numpy()
    return data


# X = get_example1()
# C, s2, Vt = pca_corr(X)
# print('Variance accumulated in each of the principal components:', calc_portion_or_variance(s2))
# print('Standard deviation of each of the principal components:',  np.std(Vt, axis=0))
# print('Standard deviation of each of the principal components:',  np.std(Vt.T, axis=0))


X = get_RCsGoff().T

std_arr = np.std(X, axis=0)
mean_arr = np.mean(X, axis=0)
X = X - mean_arr
X = X.T

num = 58581
# C, s2, Vt = pca_corr(X, num, transpose=True)
C, s2, Vt, s = pca_cov(X, num)

print('Variance accumulated in each of the principal components:', calc_portion_or_variance(s2))
pca_space = np.dot(Vt, X)

plt.scatter(pca_space[0], pca_space[1])
plt.show()
