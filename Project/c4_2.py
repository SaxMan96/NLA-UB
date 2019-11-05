import numpy as np


#strategy 2
def create_matrix(G,C,z,n,m,N):

    lbd = z[-2*m:-m]
    s = z[-m:]
    matrix = G + np.matmul(C/s*lbd,C.T)
    matrix = np.array(matrix)
    #print('matrix:{} {}'.format(type(matrix),matrix.shape))
    return matrix

def create_small_rhv(z,rhv,C,n,m,N):

    lbd = z[-2*m:-m]
    s = z[-m:]
    r1 = -rhv[:n]
    r2 = -rhv[-2*m:-m]
    r3 = -rhv[-m:]
    small_rhv = -r1 - np.matmul(-C/s,-r3+ r2*lbd)
    return small_rhv

def create_rhv(z,G,C,g,d,n,m,N):
    x = z[:n]
    lbd = z[-2*m:-m]
    s = z[-m:]
    rhv = np.zeros(N)
    rhv[:n] = np.matmul(G,x) + g - np.matmul(C,lbd)
    rhv[-2*m:-m] = s + d - np.matmul(C.T,x)
    rhv[-m:] = s*lbd
    return -rhv

def create_dz(z,dz_small,rhv,C,n,m,N):
    lbd = z[-2*m:-m]
    s = z[-m:]
    r2 = -rhv[-2*m:-m]
    r3 = -rhv[-m:]
    #print('dz_small:{} r2:{} r3:{} lbd:{} s:{} C:{}'.format(dz_small.shape,r2.shape,r3.shape,lbd.shape,s.shape,C.shape))
    dz = np.zeros(N)
    dz[:n] = dz_small
    dz[-2*m:-m] =  (-r3 + lbd*r2)/s - np.matmul(C.T,dz_small)*lbd/s
    dz[-m:] = -r2 + np.matmul(C.T,dz_small)
    return dz

def lu_np(A):
#Without pivoting
    n = len(A)
    A =A.astype('float')
    for i in range(n-1):
        aii = A[i,i]
        for j in range(i+1,n):
            A[j,i] = A[j,i]/aii
        for j in range(i+1,n):
            for k in range(i+1,n):
                A[j,k] = A[j,k]-A[j,i]*A[i,k]
    L = np.tril(A,-1)+np.eye(n)
    U = np.triu(A)
    return L,U

def ldl(A):
    n = len(A)
    L,U = lu_np(A)
    diag_U = np.array([U[i,i] for i in range(n)])
    D = np.diag(diag_U)
    LT = np.array([U[i]/diag_U[i] for i in range(n)])
    return L,D,LT

def cholo(A):
    L,D,LT = ldl(A)
    d = np.sqrt([D[i,i] for i in range(len(D))])
    G = np.array([ L[i]*d[i] for i in range(len(d))])
    GT = LT*d
    return G,GT

def frwd_subs(L,b):
    n = len(L)
    x = np.zeros(n)
    #print('b:{} L:{}'.format(type(b),type(L)))
    x[0] = b[0]/L[0,0]
    for i in range(1,n):
        x[i]= b[i]
        for j in range(i):
            x[i] = x[i] - L[i,j]*x[j]
        x[i] = x[i]/L[i,i]
    return x

def bckwd_subs(U,b):
    n = len(U)
    x = np.zeros(n)
    x[-1] =  b[-1]/U[-1,-1]
    for i in range(2,n+1):
        x[-i] = b[-i]
        for j in range(1,i):
            x[-i] = x[-i] - U[-i,-j]*x[-j]
        x[-i] = x[-i]/U[-i,-i]
    return x

def strat2_cholo(A,b):
    n = len(A)
    G,GT = cholo(A)
    GTX= frwd_subs(G,b)
    dz = bckwd_subs(GT,GTX)
    return dz

def Newton_step(lamb0,dlamb,s0,ds):
    alp=1;
    idx_lamb0=np.array(np.where(dlamb<0))
    if idx_lamb0.size>0:
        alp = min(alp,np.min(-lamb0[idx_lamb0]/dlamb[idx_lamb0]))
    idx_s0=np.array(np.where(ds<0))
    if idx_s0.size>0:
        alp = min(alp,np.min(-s0[idx_s0]/ds[idx_s0]))
    return alp


def  algorithm(z,hessian,rhv,G,C,g,d,n,m,N):
    epsilon = 1e-16
    lbd = z[-2*m:-m]
    s = z[-m:]
    #predictor substep
    small_rhv = create_small_rhv(z,rhv,C,n,m,N)
    d_z_small = strat2_cholo(hessian,small_rhv)
    d_z = create_dz(z,d_z_small,rhv,C,n,m,N)
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
    small_rhv = create_small_rhv(z,rhv,C,n,m,N)
    d_z_small = strat2_cholo(hessian,small_rhv)
    d_z = create_dz(z,d_z_small,rhv,C,n,m,N)
    #step-size correction substep
    d_lbd = d_z[-2*m:-m]
    d_s = d_z[-m:]
    alpha = Newton_step(lbd,d_lbd,s,d_s)
    #update substep
    z = z + 0.95*alpha*d_z
    hessian = create_matrix(G,C,z,n,m,N)
    rhv = create_rhv(z,G,C,g,d,n,m,N)
    norm_rL = np.linalg.norm(rhv[:n])
    #norm_rC = np.linalg.norm(rhv[n:n+p])
    if norm_rL<epsilon:
        print('norm_rL')
        return False, z, hessian, rhv
    #if norm_rC<epsilon and rhv[n:n+p].size>0 :
    #    print('rC:{} shape:{} norm_rC:{}'.format(rhv[n:n+p],rhv[n:n+p].shape,norm_rC))
     #   return False, z, hessian, rhv
    if  abs(mu)<epsilon:
        print('mu')
        return False, z, hessian, rhv
    return True, z, hessian, rhv

def c4_2(n):
    m = 2*n
    p = 0
    N = n  + 2*m
    G = np.eye(n)
    C = np.zeros((n,m))
    C[:,:n] = np.eye(n)
    C[:,n:] = -np.eye(n)
    d = np.array([-10]*m)
    g = np.random.rand(n)
    z = np.zeros(N)
    z[-2*m:] = 1
    small_KKT = create_matrix(G,C,z,n,m,N)
    niter = 100
    loop = True
    rhv = create_rhv(z,G,C,g,d,n,m,N)
    iloop=0
    while loop and iloop<niter:
        loop, z, small_KKT, rhv = algorithm(z,small_KKT,rhv,G,C,g,d,n,m,N)
        iloop = iloop+1
    print('n:{} iters:{} z+g:{}'.format(n,iloop,np.linalg.norm(z[:n]+g)))
    return
#c4_2(3)