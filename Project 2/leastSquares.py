from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import scipy.linalg as sc

def load_data_task1():
    dades = pd.read_csv("dades.txt", header=None, sep="   ", engine="python").values
    A = dades[:, 0]
    b = dades[:, 1]
    return A, b


def load_data_task2():
    dades_regressio = pd.read_csv("dades_regressio.csv", header=None).values
    A = dades_regressio[:, :-1]
    b = dades_regressio[:, -1]
    return A, b


def plot_data1(A, solution, b):
    A_sorted, b_sorted = zip(*sorted(zip(A, b)))
    plt.plot(A_sorted, b_sorted)
    A_sorted, solution_sorted = zip(*sorted(zip(A, solution)))
    plt.plot(A_sorted, solution_sorted)
    plt.gcf().set_size_inches(2, 1)
    plt.show()


def generate_poly_matrix(A, dim):
    assert isinstance(dim, int)
    assert dim > 0
    return np.vstack([A ** d for d in range(dim)]).T


def least_squares_qr(A, b, plot=1):
    Q, R = np.linalg.qr(A)
    rank = min(R.shape)
    m, n = A.shape
    assert A.shape[0] == b.shape[0]
    assert A.shape[0] >= A.shape[1]
    y = Q.T @ b
    y1 = y[:n]
    y2 = y[n:]
    R1 = R[:n, :n]
    x = sc.solve_triangular(R1, y1)
    v = np.zeros(n - rank)
    x_qr = np.concatenate((x, v))
    if plot:
        plot_data1(A[:, 1], A @ x, b)
    return x_qr


def least_squares_svd(A, b, plot=1):
    U, S, V = np.linalg.svd(A, full_matrices=False)
    n = A.shape[1]
    S[S < 1e-5] = 0
    r = np.sum(S > 0)
    S = np.hstack([1 / S[:r], np.zeros(n - r)])
    A_plus = (V.T * S).dot(U.T)
    x_svd = A_plus.dot(b)
    if plot:
        plot_data1(A[:, 1], A @ x_svd, b)
    return x_svd
    
# Full Rank Matrix
print("-" * 40)
print("Full Rank Matrix")
A, b = load_data_task1()
for dim in range(2, 6):
    print("Dim =", dim)
    A_poly = generate_poly_matrix(A, dim)
    x_svd = least_squares_svd(A_poly, b, plot=0)
    print("\tError SVD =", np.linalg.norm(A_poly.dot(x_svd) - b))
    x_qr = least_squares_qr(A_poly, b, plot=0)
    print("\tError QR =", np.linalg.norm(A_poly.dot(x_qr) - b))
print("-------- Best Polynomial Dim = 5 ---------")
print("Error SVD =", np.linalg.norm(A_poly.dot(x_svd) - b))
print("Norm(x_svd) =", np.linalg.norm(x_svd))
print("Error QR =", np.linalg.norm(A_poly.dot(x_qr) - b))
print("Norm(x_qr) =", np.linalg.norm(x_qr))
print("Error SVD - Error QR =", np.linalg.norm(A_poly.dot(x_svd) - b) - np.linalg.norm(A_poly.dot(x_qr) - b))
print("Norm(x_svd-x_qr) =", np.linalg.norm(x_svd-x_qr))


# Not Full Rank Matrix
print("-" * 40)
print("Not Full Rank Matrix")
A, b = load_data_task2()
x_svd = least_squares_svd(A, b, plot=0)
x_qr = least_squares_qr(A, b, plot=0)
print("Error SVD =", round(np.linalg.norm(A.dot(x_svd) - b), 4))
print("Norm(x_svd) =", np.linalg.norm(x_svd))
print("Error QR =", round(np.linalg.norm(A.dot(x_qr) - b), 4))
print("Norm(x_qr) =", np.linalg.norm(x_qr))
print("Error SVD - Error QR =", np.linalg.norm(A.dot(x_svd) - b) - np.linalg.norm(A.dot(x_qr) - b))
print("Norm(x_svd-x_qr) =", np.linalg.norm(x_svd-x_qr))