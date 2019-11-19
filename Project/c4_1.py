from utilities import *


def c4_1(n, verbosity=0):
    m = 2 * n
    N = n + 2 * m
    G = np.eye(n)
    C = np.zeros((n, m))
    C[:, :n] = np.eye(n)
    C[:, n:] = -np.eye(n)
    d = np.array([-10] * m)
    g = np.random.rand(n)
    z = np.zeros(N)
    z[-2 * m:] = 1
    small_KKT = create_kkt_c4_1(G, C, z, n, m)
    niter = 100
    loop = True
    rhv = create_rhv_c4_1(z, G, C, g, d, n, m, N)
    iloop = 0
    while loop and iloop < niter:
        loop, z, small_KKT, rhv = algorithm_c4_1(z, small_KKT, rhv, G, C, g, d, n, m, N, verbosity)
        iloop = iloop + 1
    if verbosity>0:
        print('n:{} iters:{} z+g:{}'.format(n, iloop, np.linalg.norm(z[:n] + g)))
    return


if __name__ == '__main__':
    c4_1(3)
