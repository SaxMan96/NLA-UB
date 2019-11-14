import numpy as np
from scipy import linalg
import scipy

def read_matrix(source,shape,symm = False):
    matrix = np.zeros(shape)
    with open(source,'r') as file:
        a = file.readlines()
    for line in a:
        row, column, value = line.strip().split()
        row = int(row)
        column = int(column)
        value = float(value)
        matrix[row-1,column-1] = value
        if symm == True:
            matrix[column-1, row-1] = value
    return matrix

def read_vector(source,n):
    v = np.zeros(n)
    with open(source,'r') as file:
        a = file.readlines()
    for line in a:
        idx,value = line.strip().split()
        idx = int(idx)
        value = float(value)
        v[idx-1] = value
    return v

def f(z, G, g, n):
    x = z[:n]
    return 0.5*np.matmul(x, G.dot(x)) + g.dot(x)

def constraint(z, A, C, b, d, n):
    x = z[:n]
    eq = np.all( np.matmul(A.T, x) - b < 1e-10 )
    ineq = np.all( np.matmul(C.T,x) - d  > -1e-10)
    return eq, ineq

def create_KKT(G, A, C, z, n, m, p, N):
    KKT = np.zeros((N, N))
    KKT[0:n , 0:n] = G
    KKT[n:n+p , 0:n] = -A.T
    KKT[0:n , n:n+p] = -A
    KKT[n+p:n+p+m , 0:n] = -C.T
    KKT[0:n , n+p:n+p+m] = -C
    KKT[-2*m:-m , -m:] = np.eye(m)
    KKT[-m:, -2*m:-m] = np.diag(z[-m:])
    KKT[-m: , -m:] = np.diag(z[-2*m:-m])
    return KKT

def create_matrix(G, A, C, z, n, m, p):
    lbd = z[-2*m:-m]
    s = z[-m:]
    matrix = np.zeros((n+p+m, n+p+m))
    matrix[:n, :n] = G
    matrix[:n, n:-m] = -A
    matrix[n:-m, :n] = -A.T
    matrix[:n, -m:] = -C
    matrix[-m:, :n] = -C.T
    matrix[-m:, -m:] = np.diag(-s/lbd)
    return matrix

def create_rhv(z, G, A, C, g, b, d, n, m, p, N):
    x = z[:n]
    gmm = z[n:n+p]
    lbd = z[-2*m:-m]
    s = z[-m:]
    F = np.zeros(N)
    F[:n] = np.matmul(G, x) + g - np.matmul(A, gmm) - np.matmul(C, lbd)
    F[n:n+p] = b - np.matmul(A.T, x)
    F[-2*m:-m] = s + d - np.matmul(C.T, x)
    F[-m:] =  s*lbd
    return -F

def create_small_rhv(rhv, z, n, m, p):
    r3 = - rhv[-2*m:-m]
    r4 = - rhv[-m:]
    lbd = z[-2*m:-m]
    small_rhv = np.zeros(n+m+p)
    small_rhv[:n+p] = rhv[:n+p]
    small_rhv[-m:] = - r3 + r4/lbd
    return small_rhv

def create_dz(small_dz, z, rhv, n, m, p, N):
    lbd = z[-2*m:-m]
    s = z[-m:]
    r4 = - rhv[-m:]
    dz = np.zeros(N)
    dz[:-m] = small_dz
    dz[-m:] = - ( r4 + s*dz[-2*m:-m])/lbd
    return dz

def Newton_step(lamb0, dlamb, s0, ds):
    alp=1;
    idx_lamb0 = np.array(np.where(dlamb<0))
    if idx_lamb0.size > 0:
        alp = min(alp, np.min(-lamb0[idx_lamb0]/dlamb[idx_lamb0]))
    idx_s0=np.array(np.where(ds<0))
    if idx_s0.size > 0:
        alp = min(alp, np.min(-s0[idx_s0]/ds[idx_s0]))
    return alp

def frwd_subs(L,b):
    n = len(L)
    x = np.zeros(n)
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

def ldl_solve(A,b): #NOT USED UNFORTUNATELY
    n = len(A)
    L,D,LT = ldl_fact(A)
    DLTX = frwd_subs(L,b)
    LTX = np.array([DLTX[i]/D[i,i] for i in range(n)])
    dz  = bckwd_subs(LT,LTX)
    return dz

def ldl_fact(A):
    A = A.astype('float')
    n = len(A)
    for i in range(n-1):
        Aii  = A[i,i]
        if abs(Aii)==0:
            print('fucking money man')
            return
        for j in range(i+1,n):
            value =  A[j,i]/Aii
            A[j,i] = value
        for j in range(i+1,n):
            for k in range(j,n):
                A[k,j] = A[k,j] - A[j,i]*A[k,i]*Aii
    L = np.tril(A,-1) + np.eye(n)
    D = np.diag(np.diag(A))
    return L,D


def pseudo_ldl_solve(A,b):
    n = len(A)
    L,D,perm = scipy.linalg.ldl(A)
    L = L[perm]
    antiperm = np.zeros(n).astype('int')
    for i in range(n):
        antiperm[perm[i]] = i
    LDLTPTX = b[perm]
    #DLTPTX = frwd_subs(L,LDLTPTX)
    DLTPTX = scipy.linalg.solve_triangular(L,LDLTPTX,lower=True)
    LTPTX = block_diagonal_solver(D,DLTPTX)
    #PTX = bckwd_subs(L.T,LTPTX)
    PTX = scipy.linalg.solve_triangular(L.T,LTPTX,lower=False)
    x = PTX[antiperm]
    return x

def block_diagonal_solver(D,b):
    n = len(D)
    x = np.zeros(n)
    i_to_pass = n
    for i in range(n):
        if i == i_to_pass:
            pass
        else:
            if i==n-1 or abs(D[i,i+1])<1e-16:
                x[i] = b[i]/D[i,i]
            else:
                xi = (b[i]*D[i+1,i+1] - b[i+1]*D[i,i+1])/(D[i,i]*D[i+1,i+1]-D[i+1,i]*D[i,i+1])
                x[i] = xi
                x[i+1] = (b[i]-D[i,i]*xi) / D[i,i+1]
                i_to_pass = i+1
    return x

def  algorithm(z, hessian, rhv, G, A, C, g, b, d, n, m, p, N):
    epsilon = 1e-16
    lbd = z[-2*m:-m]
    s = z[-m:]
    #predictor substep
    small_rhv = create_small_rhv(rhv, z, n, m, p)
    small_dz = pseudo_ldl_solve(hessian, small_rhv)
    dz = create_dz(small_dz, z, rhv, n, m, p, N)
    #step-size correction substep
    d_lbd = dz[-2*m:-m]
    d_s = dz[-m:]
    alpha = Newton_step(lbd, d_lbd, s, d_s)
    #3things
    mu = np.matmul(s, lbd)/m
    mu_tilda = np.matmul(s+alpha*d_s, lbd+alpha*d_lbd)/m
    sigma = (mu_tilda/mu)**3
    #corrector substep
    rhv[-m:] = rhv[-m:] - d_s*d_lbd + sigma*mu
    small_rhv = create_small_rhv(rhv, z, n, m, p)
    small_dz = pseudo_ldl_solve(hessian, small_rhv)
    dz = create_dz(small_dz, z, rhv, n, m, p, N)
    #step-size correction substep
    d_lbd = dz[-2*m:-m]
    d_s = dz[-m:]
    alpha = Newton_step(lbd, d_lbd, s, d_s)
    #update substep
    z = z + 0.95*alpha*dz
    hessian = create_matrix(G, A, C, z, n, m, p)
    rhv = create_rhv(z, G, A, C, g, b, d, n, m, p, N)
    norm_rL = np.linalg.norm(rhv[:n])
    norm_rC = np.linalg.norm(rhv[n:n+p])
    if norm_rL<epsilon:
        print('norm_rL:{}'.format(norm_rL))
        return False, z, hessian, rhv
    if norm_rC<epsilon and rhv[n:n+p].size>0 :
        print('norm_rC:{}'.format(norm_rC))
        return False, z, hessian, rhv
    if  abs(mu)<epsilon:
        print('mu:{}'.format(mu))
        return False, z, hessian, rhv
    return True, z, hessian, rhv

def c6():
    n = 100
    m = 2*n
    p = n//2
    N = n + p + 2*m
    A = read_matrix('optpr1/A.dad', (n, p))
    G = read_matrix('optpr1/G.dad', (n, n), True)
    C = read_matrix('optpr1/C.dad', (n, m))
    b = read_vector('optpr1/b.dad', p)
    d = read_vector('optpr1/d.dad', m)
    g = read_vector('optpr1/g_small.dad', n)
    #g = np.random.rand(n)
    z = np.zeros(N)
    z[n:] = 1
    small_KKT = create_matrix(G, A, C, z, n, m, p)
    niter = 30
    loop = True
    rhv = create_rhv(z, G, A, C, g, b, d, n, m, p, N)
    iloop=0
    while loop and iloop<niter:
        loop, z, small_KKT, rhv = algorithm(z, small_KKT, rhv, G, A, C, g, b, d, n, m, p, N)
        eq, ineq = constraint(z, A, C, b, d, n)
        print('i={} f = {} eq:{} ineq:{} x_norm:{} loop:{}'.format(iloop,f(z, G, g, n), eq, ineq, np.linalg.norm(z[:n]),loop))
        iloop = iloop+1
   # print('n:{} iters:{} z:{} g:{} z+g:{}'.format(n, iloop , -z[:n], g,np.linalg.norm(z[:n]+g)))
    with open('x.dat','w') as file:
        z[:n].tofile(file)
    return

c6()