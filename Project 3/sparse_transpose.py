
"""
Sparse matrix CSR: saxpy transpose A
"""


import numpy as np
import scipy.sparse


def smatvec(a,ja,ia,x,y,n): #computes A*x+y
    for i in range (0,n):
        for j in range (ia[i],ia[i+1]):
            y[i]=y[i]+a[j]*x[ja[j]]
    return y

def transmatvec(a,ja,ia,x,n): #Computes A^t *x
    y=np.zeros(n)
    for i in range (0,n):
        for j in range (ia[i],ia[i+1]):
            y[ja[j]]=y[ja[j]]+a[j]*x[i]
    return y
        

data    = np.array([10,-2,3,9,3,7,8,7,3,8,7,5,8,9,9,13,4,2,-1])   #val
indices = np.array([ 0, 4,0,1,5,1,2,3,0,2,3,4,1,3,4, 5,1,4, 5])   #col_ind
indptr  = np.array([ 0, 2,5,8,12,16,19])                          #row_ptr


#Important:
# Instead of n^2 elements, 
# we save only 2 #nnz + n + 1, where nnz=set of non-zero elements of A

n=6;
A = scipy.sparse.csr_matrix((data,indices,indptr),shape=(n,n))
x=np.ones(n)
y=np.zeros(n) 
yy=np.ones(n) 

print("A:")
print(A.todense())

print("\n")
print("x:")
print(x)
print("y:")
print(y)
print("yy:")
print(yy)

print("\n")
print("A*x:")
print(smatvec(data,indices,indptr,x,y,n))
print("check:")
print(A.dot(x))
print("A*x+yy:")
print(smatvec(data,indices,indptr,x,yy,n))


print("\n")
print("A^t *x:")
print(transmatvec(data,indices,indptr,x,n))
print("check:")
print((A.T).dot(x))
print("\n")
