import numpy as np
import scipy as sc
import math
from sklearn.datasets import make_regression
from matplotlib import pyplot as plt


def ls_with_qr(am, bv):
    qm, rm = sc.linalg.qr(am)
    m, n = am.shape
    # print('a:{} b:{} q:{} r:{}'.format(A.shape,b.shape,Q.shape,R.shape))
    yv = qm.T.dot(bv)
    y1 = yv[:n]
    residual = np.linalg.norm(yv[n:], 2)

    print("residual", residual)
    return sc.linalg.solve_triangular(rm[:n], y1)


def ls_with_svd(xm, bv):
    um, sigma, vmt = np.linalg.svd(xm, full_matrices=False)
    n = xm.shape[1]
    r = np.sum(sigma > 1e-10)
    # print(r, len(sigma), n)
    sigma = np.hstack([1 / sigma[:r], np.zeros(n - r)])
    xm_plus = (vmt.T * sigma).dot(um.T)
    return xm_plus.dot(bv)


def get_data_1(polinomial=3):
    xv = []
    bv = []
    with open('../dades.txt') as dades:
        for line in dades:
            x1, x2 = line.split("   ")
            xv.append(float(x1))
            bv.append(float(x2))

    am = np.zeros((len(xv), polinomial))
    for i in range(polinomial):
        am[:, i] = np.power(xv, i)
    return am, bv, xv


def get_data_2():
    n = 15
    m = 11
    am = np.zeros((n, m))
    bv = np.zeros(n)
    with open('./dades_regressio.csv') as f:
        for idx, line in enumerate(f):
            data = line.split(',')
            am[idx, :] = [float(num) for num in data[:-1]]
            bv[idx] = float(data[-1])
    return am, bv


def get_data_generated():
    # generated regression data
    xm, yv, _ = make_regression(
        n_samples=50,
        n_features=1,
        n_informative=1,
        n_targets=1,
        noise=5,
        coef=True,
        random_state=1
    )
    return xm, yv


X, b, y = get_data_1()
w_svd = ls_with_svd(X, b)
w_qr = ls_with_qr(X, b)

error_svd = np.linalg.norm(X.dot(w_svd) - y, ord=2)
error_qr = np.linalg.norm(X.dot(w_qr) - y, ord=2)
print(error_svd, error_qr, abs(error_qr - error_svd))


# plt.scatter(x, b)
# plt.plot(X[0], w * X[0], c='red')
# plt.show()


def f_a(a):
    def f(x):
        acc = 0
        for i in range(len(a)):
            acc += a[i] * np.power(x, i)
        return acc

    return f


line_x = np.array(range(math.ceil(max(y))))
f_approx = f_a(w_svd)

plt.scatter(y, b)
plt.plot(line_x, f_approx(line_x), 'r')
plt.show()
