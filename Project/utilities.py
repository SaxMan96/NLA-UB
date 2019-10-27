import numpy as np


def create_kkt(G, A, C, z, n, m, p, N):
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


def create_rhv(z, G, A, C, g, b, d, n, m, p, N):
    x = z[:n]
    gmm = z[n:n + p]
    lbd = z[-2 * m:-m]
    s = z[-m:]
    F = np.zeros(N)
    F[:n] = np.matmul(G, x) + g - np.matmul(A, gmm) - np.matmul(C, lbd)
    F[n:n + p] = b - np.matmul(A.T, x)
    F[-2 * m:-m] = s + d - np.matmul(C.T, x)
    F[-m:] = s * lbd
    return F


def newton_step(lamb0, dlamb, s0, ds):
    alp = 1;
    idx_lamb0 = np.array(np.where(dlamb < 0))
    if idx_lamb0.size > 0:
        alp = min(alp, np.min(-lamb0[idx_lamb0] / dlamb[idx_lamb0]))
    idx_s0 = np.array(np.where(ds < 0))
    if idx_s0.size > 0:
        alp = min(alp, np.min(-s0[idx_s0] / ds[idx_s0]))
    return alp


def algorithm(z, hessian, rhv, G, A, C, g, b, d, n, m, p, N):
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
    hessian = create_kkt(G, A, C, z, n, m, p, N)
    rhv = -create_rhv(z, G, A, C, g, b, d, n, m, p, N)
    norm_rL = np.linalg.norm(rhv[:n])
    norm_rC = np.linalg.norm(rhv[n:n + p])
    if norm_rL < epsilon:
        print('norm_rL')
        return False, z, hessian, rhv
    if norm_rC < epsilon and rhv[n:n + p].size > 0:
        print('rC:{} shape:{} norm_rC:{}'.format(rhv[n:n + p], rhv[n:n + p].shape, norm_rC))
        return False, z, hessian, rhv
    if abs(mu) < epsilon:
        print('mu')
        return False, z, hessian, rhv
    return True, z, hessian, rhv
