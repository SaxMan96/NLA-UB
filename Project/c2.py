from utilities import *


def c2(n, verbosity=0):
    m = 2 * n
    p = 0
    N = n + p + 2 * m
    A = np.zeros((n, p))
    G = np.eye(n)
    C = np.zeros((n, m))
    C[:, :n] = np.eye(n)
    C[:, n:] = -np.eye(n)
    b = np.zeros(p)
    d = np.array([-10] * m)
    g = np.random.rand(n)
    z = np.zeros(N)
    z[-2 * m:] = 1
    #     def create_KKT(G,A,C,z,n,m,p,N):
    kkt = create_kkt_c2(G, A, C, z, n, m, p, N)
    niter = 100
    loop = True
    rhv = -create_rhv_c2(z, G, A, C, g, b, d, n, m, p, N)
    iloop = 0
    while loop and iloop < niter:
        loop, z, kkt, rhv = algorithm_c2(z, kkt, rhv, G, A, C, g, b, d, n, m, p, N)
        iloop = iloop + 1
    if verbosity:
        
        print('n: {}\niters: {}\nz+g: {}'.format(n, iloop, np.linalg.norm(z[:n] + g)))


if __name__ == '__main__':
    c2(100)
