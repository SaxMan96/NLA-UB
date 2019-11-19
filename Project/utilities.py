import numpy as np


def create_dz_c4_1(z, dz_small, rhv, n, m, N):
    lbd = z[-2 * m:-m]
    s = z[-m:]
    r3 = -rhv[-m:]
    dz = np.zeros(N)
    dz[:n + m] = dz_small
    dz[-m:] = (-r3 - s * dz[-2 * m:-m]) / lbd
    return dz


def create_dz_c4_2(z, dz_small, rhv, C, n, m, N):
    lbd = z[-2 * m:-m]
    s = z[-m:]
    r2 = -rhv[-2 * m:-m]
    r3 = -rhv[-m:]
    dz = np.zeros(N)
    dz[:n] = dz_small
    dz[-2 * m:-m] = (-r3 + lbd * r2) / s - np.matmul(C.T, dz_small) * lbd / s
    dz[-m:] = -r2 + np.matmul(C.T, dz_small)
    return dz


def create_kkt_c2(G, A, C, z, n, m, p, N):
    kkt = np.zeros((N, N))
    kkt[0:n, 0:n] = G
    kkt[n:n + p, 0:n] = -A.T
    kkt[0:n, n:n + p] = -A
    kkt[n + p:n + p + m, 0:n] = -C.T
    kkt[0:n, n + p:n + p + m] = -C
    kkt[-2 * m:-m, -m:] = np.eye(m)
    kkt[-m:, -2 * m:-m] = np.diag(z[-m:])
    kkt[-m:, -m:] = np.diag(z[-2 * m:-m])
    return kkt


# strategy 1
def create_kkt_c4_1(G, C, z, n, m):
    lbd = z[-2 * m:-m]
    s = z[-m:]
    kkt = np.zeros((n + m, n + m))
    kkt[:n, :n] = G
    kkt[:n, n:] = -C
    kkt[n:, :n] = -C.T
    kkt[-m:, -m:] = np.diag(-s / lbd)
    return kkt


def create_kkt_c4_2(G, C, z, m, ):
    lbd = z[-2 * m:-m]
    s = z[-m:]
    kkt = G + np.matmul(C / s * lbd, C.T)
    kkt = np.array(kkt)
    return kkt


def create_rhv_c2(z, G, A, C, g, b, d, n, m, p, N):
    x = z[:n]
    gmm = z[n:n + p]
    lbd = z[-2 * m:-m]
    s = z[-m:]
    rvh = np.zeros(N)
    rvh[:n] = np.matmul(G, x) + g - np.matmul(A, gmm) - np.matmul(C, lbd)
    rvh[n:n + p] = b - np.matmul(A.T, x)
    rvh[-2 * m:-m] = s + d - np.matmul(C.T, x)
    rvh[-m:] = s * lbd
    return rvh


def create_rhv_c4_1(z, G, C, g, d, n, m, N):
    x = z[:n]
    lbd = z[-2 * m:-m]
    s = z[-m:]
    rhv = np.zeros(N)
    rhv[:n] = np.matmul(G, x) + g - np.matmul(C, lbd)
    rhv[-2 * m:-m] = s + d - np.matmul(C.T, x)
    rhv[-m:] = s * lbd
    return -rhv


def create_rhv_small_c4_1(z, rhv, n, m, N):
    lbd = z[-2 * m:-m]
    r1 = -rhv[:n]
    r2 = -rhv[-2 * m:-m]
    r3 = -rhv[-m:]
    small_rhv = np.zeros(n + m)
    small_rhv[:n] = -r1
    small_rhv[-m:] = -r2 + r3 / lbd
    return small_rhv


def create_small_rhv_c4_2(z, rhv, C, n, m, N):
    lbd = z[-2 * m:-m]
    s = z[-m:]
    r1 = -rhv[:n]
    r2 = -rhv[-2 * m:-m]
    r3 = -rhv[-m:]
    small_rhv = -r1 - np.matmul(-C / s, -r3 + r2 * lbd)
    return small_rhv


def lu_np(A):
    # Without pivoting
    n = len(A)
    A = A.astype('float')
    for i in range(n - 1):
        aii = A[i, i]
        for j in range(i + 1, n):
            A[j, i] = A[j, i] / aii
        for j in range(i + 1, n):
            for k in range(i + 1, n):
                A[j, k] = A[j, k] - A[j, i] * A[i, k]
    L = np.tril(A, -1) + np.eye(n)
    U = np.triu(A)
    return L, U


def ldl(A):
    n = len(A)
    L, U = lu_np(A)
    diag_u = np.array([U[i, i] for i in range(n)])
    D = np.diag(diag_u)
    LT = np.array([U[i] / diag_u[i] for i in range(n)])
    return L, D, LT


def cholo(A):
    L, D, LT = ldl(A)
    d = np.sqrt([D[i, i] for i in range(len(D))])
    G = np.array([L[i] * d[i] for i in range(len(d))])
    GT = LT * d
    return G, GT


def strat1_ldl(A, b):
    n = len(A)
    L, D, LT = ldl(A)
    DLTX = frwd_subs(L, b)
    LTX = np.array([DLTX[i] / D[i, i] for i in range(n)])
    dz = bckwd_subs(LT, LTX)
    return dz


def frwd_subs(L, b):
    n = len(L)
    x = np.zeros(n)
    x[0] = b[0] / L[0, 0]
    for i in range(1, n):
        x[i] = b[i]
        for j in range(i):
            x[i] = x[i] - L[i, j] * x[j]
        x[i] = x[i] / L[i, i]
    return x


def bckwd_subs(U, b):
    n = len(U)
    x = np.zeros(n)
    x[-1] = b[-1] / U[-1, -1]
    for i in range(2, n + 1):
        x[-i] = b[-i]
        for j in range(1, i):
            x[-i] = x[-i] - U[-i, -j] * x[-j]
        x[-i] = x[-i] / U[-i, -i]
    return x


def strat2_cholo(A, b):
    n = len(A)
    G, GT = cholo(A)
    GTX = frwd_subs(G, b)
    dz = bckwd_subs(GT, GTX)
    return dz


def newton_step(lamb0, dlamb, s0, ds):
    alp = 1
    idx_lamb0 = np.array(np.where(dlamb < 0))
    if idx_lamb0.size > 0:
        alp = min(alp, np.min(-lamb0[idx_lamb0] / dlamb[idx_lamb0]))
    idx_s0 = np.array(np.where(ds < 0))
    if idx_s0.size > 0:
        alp = min(alp, np.min(-s0[idx_s0] / ds[idx_s0]))
    return alp


def algorithm_c2(z, hessian, rhv, G, A, C, g, b, d, n, m, p, N, verbosity=0):
    epsilon = 1e-16
    lbd = z[-2 * m:-m]
    s = z[-m:]
    # predictor substep
    d_z = np.linalg.solve(hessian, rhv)
    # step-size correction substep
    d_lbd = d_z[-2 * m:-m]
    d_s = d_z[-m:]
    alpha = newton_step(lbd, d_lbd, s, d_s)
    # 3things
    mu = np.matmul(s, lbd) / m
    mu_tilda = np.matmul(s + alpha * d_s, lbd + alpha * d_lbd) / m
    sigma = (mu_tilda / mu) ** 3
    # corrector substep
    rhv[-m:] = rhv[-m:] - d_s * d_lbd + sigma * mu
    d_z = np.linalg.solve(hessian, rhv)
    # step-size correction substep
    d_lbd = d_z[-2 * m:-m]
    d_s = d_z[-m:]
    alpha = newton_step(lbd, d_lbd, s, d_s)
    # update substep
    z = z + 0.95 * alpha * d_z
    hessian = create_kkt_c2(G, A, C, z, n, m, p, N)
    rhv = -create_rhv_c2(z, G, A, C, g, b, d, n, m, p, N)
    norm_rL = np.linalg.norm(rhv[:n])
    norm_rC = np.linalg.norm(rhv[n:n + p])
    if norm_rL < epsilon:
        if verbosity:
            print('norm_rL')
        return False, z, hessian, rhv
    if norm_rC < epsilon and rhv[n:n + p].size > 0:
        if verbosity:
            print('rC:{} shape:{} norm_rC:{}'.format(rhv[n:n + p], rhv[n:n + p].shape, norm_rC))
        return False, z, hessian, rhv
    if abs(mu) < epsilon:
        if verbosity:
            print('mu')
        return False, z, hessian, rhv
    return True, z, hessian, rhv


def algorithm_c4_1(z, hessian, rhv, G, C, g, d, n, m, N, verbosity=0):
    epsilon = 1e-16
    lbd = z[-2 * m:-m]
    s = z[-m:]
    # predictor substep
    small_rhv = create_rhv_small_c4_1(z, rhv, n, m, N)
    d_z_small = strat1_ldl(hessian, small_rhv)
    d_z = create_dz_c4_1(z, d_z_small, rhv, n, m, N)
    # step-size correction substep
    d_lbd = d_z[-2 * m:-m]
    d_s = d_z[-m:]
    alpha = newton_step(lbd, d_lbd, s, d_s)
    # 3things
    mu = np.matmul(s, lbd) / m
    mu_tilda = np.matmul(s + alpha * d_s, lbd + alpha * d_lbd) / m
    sigma = (mu_tilda / mu) ** 3
    # corrector substep
    rhv[-m:] = rhv[-m:] - d_s * d_lbd + sigma * mu
    small_rhv = create_rhv_small_c4_1(z, rhv, n, m, N)
    d_z_small = strat1_ldl(hessian, small_rhv)
    d_z = create_dz_c4_1(z, d_z_small, rhv, n, m, N)
    # step-size correction substep
    d_lbd = d_z[-2 * m:-m]
    d_s = d_z[-m:]
    alpha = newton_step(lbd, d_lbd, s, d_s)
    # update substep
    z = z + 0.95 * alpha * d_z
    hessian = create_kkt_c4_1(G, C, z, n, m)
    rhv = create_rhv_c4_1(z, G, C, g, d, n, m, N)
    norm_rL = np.linalg.norm(rhv[:n])
    # norm_rC = np.linalg.norm(rhv[n:n+p])
    if norm_rL < epsilon:
        if verbosity:
            print('norm_rL')
        return False, z, hessian, rhv
    # if norm_rC<epsilon and rhv[n:n+p].size>0 :
    #    print('rC:{} shape:{} norm_rC:{}'.format(rhv[n:n+p],rhv[n:n+p].shape,norm_rC))
    #   return False, z, hessian, rhv
    if abs(mu) < epsilon:
        if verbosity:
            print('mu')
        return False, z, hessian, rhv
    return True, z, hessian, rhv


def algorithm_c4_2(z, hessian, rhv, G, C, g, d, n, m, N, verbosity=0):
    epsilon = 1e-16
    lbd = z[-2 * m:-m]
    s = z[-m:]
    # predictor substep
    small_rhv = create_small_rhv_c4_2(z, rhv, C, n, m, N)
    d_z_small = strat2_cholo(hessian, small_rhv)
    d_z = create_dz_c4_2(z, d_z_small, rhv, C, n, m, N)
    # step-size correction substep
    d_lbd = d_z[-2 * m:-m]
    d_s = d_z[-m:]
    alpha = newton_step(lbd, d_lbd, s, d_s)
    # 3things
    mu = np.matmul(s, lbd) / m
    mu_tilda = np.matmul(s + alpha * d_s, lbd + alpha * d_lbd) / m
    sigma = (mu_tilda / mu) ** 3
    # corrector substep
    rhv[-m:] = rhv[-m:] - d_s * d_lbd + sigma * mu
    small_rhv = create_small_rhv_c4_2(z, rhv, C, n, m, N)
    d_z_small = strat2_cholo(hessian, small_rhv)
    d_z = create_dz_c4_2(z, d_z_small, rhv, C, n, m, N)
    # step-size correction substep
    d_lbd = d_z[-2 * m:-m]
    d_s = d_z[-m:]
    alpha = newton_step(lbd, d_lbd, s, d_s)
    # update substep
    z = z + 0.95 * alpha * d_z
    hessian = create_kkt_c4_2(G, C, z, m, )
    rhv = create_rhv_c4_1(z, G, C, g, d, n, m, N)
    norm_rL = np.linalg.norm(rhv[:n])
    # norm_rC = np.linalg.norm(rhv[n:n+p])
    if norm_rL < epsilon:
        if verbosity:
            print('norm_rL')
        return False, z, hessian, rhv
    # if norm_rC<epsilon and rhv[n:n+p].size>0 :
    #    print('rC:{} shape:{} norm_rC:{}'.format(rhv[n:n+p],rhv[n:n+p].shape,norm_rC))
    #   return False, z, hessian, rhv
    if abs(mu) < epsilon:
        if verbosity:
            print('mu')
        return False, z, hessian, rhv
    return True, z, hessian, rhv
