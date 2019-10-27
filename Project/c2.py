import numpy as np


def create_KKT(G,A,C,z,n,m,p,N):
    KKT = np.zeros((N,N))
    KKT[0:n , 0:n] = G
    KKT[n:n+p , 0:n] = -A.T
    KKT[0:n , n:n+p] = -A
    KKT[n+p:n+p+m , 0:n] = -C.T
    KKT[0:n , n+p:n+p+m] = -C
    KKT[-2*m:-m , -m:] = np.eye(m)
    KKT[-m:, -2*m:-m] = np.diag(z[-m:])
    KKT[-m: , -m:] = np.diag(z[-2*m:-m])
    return KKT

def create_rhv(z,G,A,C,g,b,d,n,m,p,N):
    x = z[:n]
    gmm = z[n:n+p]
    lbd = z[-2*m:-m]
    s = z[-m:]
    F = np.zeros(N)
    F[:n] = np.matmul(G,x) + g - np.matmul(A,gmm) - np.matmul(C,lbd)
    F[n:n+p] = b - np.matmul(A.T,x)
    F[-2*m:-m] = s + d - np.matmul(C.T,x)
    F[-m:] =  s*lbd
    return F

def Newton_step(lamb0,dlamb,s0,ds):
    alp=1;
    idx_lamb0=np.array(np.where(dlamb<0))
    if idx_lamb0.size>0:
        alp = min(alp,np.min(-lamb0[idx_lamb0]/dlamb[idx_lamb0]))
    idx_s0=np.array(np.where(ds<0))
    if idx_s0.size>0:
        alp = min(alp,np.min(-s0[idx_s0]/ds[idx_s0]))
    return alp


def  algorithm(z,hessian,rhv,G,A,C,g,b,d,n,m,p,N):
    epsilon = 1e-16
    lbd = z[-2*m:-m]
    s = z[-m:]
    #predictor substep
    d_z = np.linalg.solve(hessian,rhv)
    #step-size correction substep
    d_lbd = d_z[-2*m:-m]
    d_s = d_z[-m:]
    alpha = Newton_step(lbd,d_lbd,s,d_s)
    #3things
    mu = np.matmul(s,lbd)/m
    mu_tilda = np.matmul(s+alpha*d_s,lbd+alpha*d_lbd)/m
    sigma = (mu_tilda/mu)**3
    #corrector substep
    rhv[-m:] = rhv[-m:] - d_s*d_lbd + sigma*mu
    d_z = np.linalg.solve(hessian,rhv)
    #step-size correction substep
    d_lbd = d_z[-2*m:-m]
    d_s = d_z[-m:]
    alpha = Newton_step(lbd,d_lbd,s,d_s)
    #update substep
    z = z + 0.95*alpha*d_z
    hessian = create_KKT(G,A,C,z,n,m,p,N)
    rhv = -create_rhv(z,G,A,C,g,b,d,n,m,p,N)
    norm_rL = np.linalg.norm(rhv[:n])
    norm_rC = np.linalg.norm(rhv[n:n+p])
    if norm_rL<epsilon:
        print('norm_rL')
        return False, z, hessian, rhv
    if norm_rC<epsilon and rhv[n:n+p].size>0 :
        print('rC:{} shape:{} norm_rC:{}'.format(rhv[n:n+p],rhv[n:n+p].shape,norm_rC))
        return False, z, hessian, rhv
    if  abs(mu)<epsilon:
        print('mu')
        return False, z, hessian, rhv
    return True, z, hessian, rhv

def c2(n):
    m = 2*n
    p = 0
    N = n + p + 2*m
    A = np.zeros((n,p))
    G = np.eye(n)
    C = np.zeros((n,m))
    C[:,:n] = np.eye(n)
    C[:,n:] = -np.eye(n)
    b = np.zeros(p)
    d = np.array([-10]*m)
    g = np.random.rand(n)
    z = np.zeros(N)
    z[-2*m:] = 1
#     def create_KKT(G,A,C,z,n,m,p,N):
    KKT = create_KKT(G,A,C,z,n,m,p,N)
    niter = 100
    loop = True
    rhv = -create_rhv(z,G,A,C,g,b,d,n,m,p,N)
    iloop=0
    while loop and iloop<niter:
        loop, z, KKT, rhv = algorithm(z,KKT,rhv,G,A,C,g,b,d,n,m,p,N)
        iloop = iloop+1
    print('n:{} iters:{} z+g:{}'.format(n,iloop,np.linalg.norm(z[:n]+g)))
    return

#c2(10)
