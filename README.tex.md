# Problems

- Linear System: $Ax = b$
- Least Squares Problem $||Ax - b||_{2}$
- Eigenvalue problem $Ax = \lambda x$, where $x$ - vector $\lambda$ - scalar
- Singular Value problem $A^T Ax = \lambda x$

# Definitions

## Nonsingular matrix
A square matrix that is not singular, i.e. one that has a matrix inverse. Nonsingular matrices are sometimes also called regular matrices. A square matrix is nonsingular if its determinant is nonzero

## Inversing Triangular Matrix

Solve directly:

$$
    \left(\begin{array}{lll}
    {a} & {b} & {c} \\
    {0} & {d} & {e} \\
    {0} & {0} & {f}
    \end{array}\right)
    \left(\begin{array}{lll}
    {x} & {y} & {z} \\
    {0} & {y} & {v} \\
    {0} & {0} & {w}
    \end{array}\right)
    =
    \left(\begin{array}{lll}
    {1} & {0} & {0} \\
    {0} & {1} & {0} \\
    {0} & {0} & {1}
    \end{array}\right)
$$

$$
    \left(\begin{array}{lll}
    {1 / a} & {-b /(a d)} & {(b e-c d) /(a f d)} \\
    {0} & {1 / d} & {-e /(f d)} \\
    {0} & {0} & {1 / f}
    \end{array}\right)
$$
from which we see directly that the matrix is invertible if all $a,d$ and $f$ are different from zero

## Positive Definite matrix

$A$ $n \times n$ symmetric real matrix $M$ is said to be positive definite if $x^{\top} M x>0$ for all non-zero $x$ in $\mathbb{R}^{n}$.

 $M \text { positive definite } \Longleftrightarrow x^{\top} M x>0 \text { for all } x \in \mathbb{R}^{n} \backslash {0}$

## Positive Semidefinite matrix 

$A$ $n \times n$ symmetric real matrix $M$ is said to be positive semidefinite or non-negative definite if $x^{\top} M x \geq 0$ for all $x$ in $\mathbb{R}^{n} .$ Formally,

$M \text { positive definite } \Longleftrightarrow x^{\top} M x>0 \text { for all } x \in \mathbb{R}^{n}$

## Rank of a matrix

The rank of a matrix is defined as 

  - the maximum number of linearly independent column vectors in the matrix or 
  - the maximum number of linearly independent row vectors in the matrix. 

Both definitions are equivalent.

## Eigenvalues and eigenvectors

Eigenvalues (*characteristic roots, proper values, latent roots)* are a special set of scalars associated with a linear system of equations. Each eigenvalue is paired with a corresponding so-called eigenvector. 

  Let $A$ be a linear transformation represented by a matrix A. If there is a vector $\mathbf{X} \in \mathbb{R}^{n} \neq 0$ such that:

$$
  A X=\lambda X
$$

  for some scalar $\lambda,$ then $\lambda$ is called the eigenvalue of A with corresponding (right) eigenvector $\mathbf{X}$.

  Letting A be a $k \times k$ square matrix
$$
  \left[\begin{array}{cccc}
  {a_{11}} & {a_{12}} & {\cdots} & {a_{1 k}} \\
  {a_{21}} & {a_{22}} & {\cdots} & {a_{2 k}} \\
  {\vdots} & {\vdots} & {\ddots} & {\vdots} \\
  {a_{k 1}} & {a_{k 2}} & {\cdots} & {a_{k k}}
  \end{array}\right]
$$
  with eigenvalue $\lambda,$ then the corresponding eigenvectors satisfy
$$
  \left[\begin{array}{cccc}
  {a_{11}} & {a_{12}} & {\cdots} & {a_{1 k}} \\
  {a_{21}} & {a_{22}} & {\cdots} & {a_{2 k}} \\
  {\vdots} & {\vdots} & {\ddots} & {\vdots} \\
  {a_{k 1}} & {a_{k 2}} & {\cdots} & {a_{k k}}
  \end{array}\right]\left[\begin{array}{c}
  {x_{1}} \\
  {x_{2}} \\
  {\vdots} \\
  {x_{k}}
  \end{array}\right]=\lambda\left[\begin{array}{c}
  {x_{1}} \\
  {x_{2}} \\
  {\vdots} \\
  {x_{k}}
  \end{array}\right]
$$
  which is equivalent to the homogeneous system
$$
  \left[\begin{array}{cccc}
  {a_{11}-\lambda} & {a_{12}} & {\cdots} & {a_{1 k}} \\
  {a_{21}} & {a_{22}-\lambda} & {\cdots} & {a_{2 k}} \\
  {\vdots} & {\vdots} & {\ddots} & {\vdots} \\
  {a_{k 1}} & {a_{k 2}} & {\cdots} & {a_{k k}-\lambda}
  \end{array}\right]\left[\begin{array}{c}
  {x_{1}} \\
  {x_{2}} \\
  {\vdots} \\
  {x_{k}}
  \end{array}\right]=\left[\begin{array}{c}
  {0} \\
  {0} \\
  {\vdots} \\
  {0}
  \end{array}\right]
$$
  Equation ( 4) can be written compactly as
$$
  (A-\lambda I) X=0
$$

## Singular Value

## Condition number

## Norms
$$
  \begin{aligned}
  &|\mathbf{x}|_{\infty} \equiv \max \left|x_{i}\right| \\
  &|\mathbf{x}|_{p} \equiv\left(\sum_{i}\left|x_{i}\right|^{p}\right)^{1 / p}
  \end{aligned}
$$

$$
  \begin{array}{llll}
  {\text { name }} & {\text { symbol}} & {\text {  value }} & {\text {    approx. }} \\
  {L^{1}-\text { norm }} & {|\mathbf{x}|_{1}} & {6} & {6.000} \\
  {L^{2}-\text { norm }} & {|\mathbf{x}|_{2}} & {\sqrt{14}} & {3.742} \\
  {L^{3}-\text { norm }} & {|\mathbf{x}|_{3}} & {6^{2 / 3}} & {3.302} \\
  {L^{4}-\text { norm }} & {|\mathbf{x}|_{4}} & {2^{1 / 4} \sqrt{7}} &{3.146} \\
  {L^{\infty}-\text { norm }} & {|\mathbf{x}|_{\infty}} & {3} & {3.000}
  \end{array}
$$

## Rank1/2 approximation

# Matrix Factorization:

A factorization of the matrix $A$ is a representation of $A$ as a product of several "simpler" matrices, which make the problem at hand easier to solve. We give two examples.

$$
\left[\begin{array}{cccc}
  {a_{11}} & {} & {} & {} \\ 
  {a_{21}} &  {a_{22}} & {} & {} \\ 
  {\vdots} & {\vdots} & {\ddots} & {} \\ 
  {a_{n 1}} &   {a_{n 2}} &  {\ldots} & {a_{n n}}
  \end{array}\right]
  \left[\begin{array}{c}
  {x_{1}} \\ 
  {x_{2}} \\ 
  {\vdots} \\
  {x_{n}}
  \end{array}\right]
  =\left[\begin{array}{c}
  {b_{1}} \\ 
  {b_{2}} \\   
  {\vdots} \\ 
  {b_{n}}
  \end{array}\right]
$$

## Forward Substitution


$$
  \text{for } i=1 \text{ to } n \\
  \quad x_{i}=\left(b_{i}-\sum_{k=1}^{i-1} a_{i k} x_{k}\right) / a_{i i} \\
  \text{end for}
$$

## Backward Substitution

  An analogous idea, back substitution, works if $A$ is upper triangular. 
To use this to solve a general system $A x=b$ we need the following matrix factorization, which is just a restatement of Gaussian elimination.

## LU Factorization

LU on example:
- rows operations
  $$
  A = \left[ \begin{array} { r r r } { 1 } & { 4 } & { - 3 } \\ { - 2 } & { 8 } & { 5 } \\ { 3 } & { 4 } & { 7 } \end{array} \right] \begin{matrix}\\\textbf{2}R_{1}+R_{2}\\\textbf{-3}R_{1}+R_{2}\end{matrix} \quad L = \left[ \begin{array} { c c c } { 1 } & { 0 } & { 0 } \\ { ? } & { 1 } & { 0 } \\ { ? } & { ? } & { 1 } \end{array} \right]
  $$
  
- **2** and **-3** goes down to $L$ with changed sign
  $$
U = \left[ \begin{array} { r r r } { 1 } & { 4 } & { - 3 } \\ {0 } & { 16 } & { -1 } \\ { 0 } & { -8 } & { 16 } \end{array} \right] \begin{matrix}\\ \\\textbf{0.5}R_{2}+R_{3}\end{matrix} \quad L = \left[ \begin{array} { c c c } { 1 } & { 0 } & { 0 } \\ { - 2 } & { 1 } & { 0 } \\ { 3 } & { ? } & { 1 } \end{array} \right]
  $$
  
- **0.5** goes to $L$ with minus sign.
  $$
  U = \left[ \begin{array} { r r r } { 1 } & { 4 } & { - 3 } \\ {0 } & { 16 } & { -1 } \\ { 0 } & { 0 } & { 15.5 } \end{array} \right] \quad L = \left[ \begin{array} { c c c } { 1 } & { 0 } & { 0 } \\ { - 2 } & { 1 } & { 0 } \\ { 3 } & { -0.5 } & { 1 } \end{array} \right]
  $$

## Gaussian Elimination

**Partial Pivoting** - Sort rows by its first elements absolute value. Choose $a_{i,k}$ and swap that  $k$-th row with $i$-th row `for i in (1,n) `.

**THEOREM** If the $n-b y-n$ matrix $A$ is *nonsingular*, there exists 

- a permutation matrix $P$ (the identity matrix with its rows permuted), 
- a nonsingular lower triangular matrix $L,$ and 
- a nonsingular upper triangular matrix $U$ 

such that $A=P \cdot L \cdot U .$ To solve $A x=b,$ we solve the equivalent system $P L U x=b$ as follows:

- $\left.L U x=P^{-1} b=P^{T} b \quad \text { (permute entries of } b\right)$
- $U x=L^{-1}\left(P^{T} b\right) \quad \text{(forward substitution),}$
- $x=U^{-1}\left(L^{-1} P^{T} b\right) \quad \text{(back substitution).}$

Important thing $L^{-1} \neq L^{T}$ and  $U^{-1} \neq U^{T}$ - check **inversion of triangular matrix** in section definitions.

## Cholesky Factorization

The Cholesky decomposition is roughly twice as efficient as the LU decomposition for solving systems of linear equations. $\bf{A}$ needs to be symmetric. Every symmetric, positive definite matrix A can be decomposed into a product of a unique lower triangular matrix L and its transpose. $A = LL^{T}$, where $L$ is a lower triangular matrix with real and positive diagonal entries.

### Algorithm:

$$
\begin{aligned}
l_{1,1} &=\sqrt{a_{11}} \\
l_{j, 1} &=\frac{a_{j 1}}{l_{11}}, \quad j \in[2, n] \\
l_{i, i} &=\sqrt{a_{i i}-\sum_{p=1}^{i-1} l_{i p}^{2}}, \quad i \in[2, n] \\
l_{j, i} &=\left(a_{j i}-\sum_{p=1}^{i-1} l_{i p} l_{j p}\right) / l_{i i}, \quad i \in[2, n-1], j \in[i+1, n]
\end{aligned}
$$

### Example

$A = \left[\begin{array}{ccc}
{4} & {12} & {-16} \\
{12} & {37} & {-43} \\
{-16} & {-43} & {98}
\end{array}\right]$ 

- $i = 1:$

  $$
  \begin{aligned}
    & l_{1,1} = \sqrt{a_{1,1}} = 2 \\
    & l_{2,1} = \frac{a_{2,1}}{l_{1,1}} = 6 \\
    & l_{3,1} = \frac{a_{3,1}}{l_{1,1}} = -8
    \end{aligned} \quad L = \left[ \begin{array} { r r r } { 2 } & { 0 } & { 0 } \\ { 6 } & { ? } & { 0 } \\ { - 8 } & { ? } & { ? } \end{array} \right]
  $$

- $i = 2:$

  $$
  \begin{aligned}
    & { l _ { 2,2 } = \sqrt { a _ { 2,2 } -  l _ { 2,1 } ^ { 2 }} }  = \sqrt{37 - 6^{2}} = 1\\ 
    & { l_ { 3,2 } = \left( a _ { 3,2 } - l _ { 2 , 1 } \cdot l _ { 3,1 } \right) / l _ { 2,2 } } = (- 43 - (6 \cdot (- 8))/1 = 5 \end{aligned}
     \quad L = \left[ \begin{array} { r r r } { 2 } & { 0 } & { 0 } \\ { 6 } & { 1 } & { 0 } \\ { - 8 } & { 5 } & { ? } \end{array} \right] 
  $$

- $i = 3$

  $$
  \begin{aligned} 
  l _ { 3,3 } & = \sqrt { a _ { 3,3 } - \left( l _ { 3,1 } ^ { 2 } + l _ { 3,2 } ^ { 2 } \right) } = \\ & = \sqrt { 98 - \left( ( - 8 ) ^ { 2 } + 5 ^ { 2 } \right) } = \\ & = \sqrt { 98 - ( 64 + 25 ) } = \\ & = \sqrt { 98 - 89 } = \sqrt { 9 } = 3 \end{aligned}
  $$

  

  $$L = \left[ \begin{array} { r r r } { 2 } & { 0 } & { 0 } \\ { 6 } & { 1 } & { 0 } \\ { - 8 } & { 5 } & { 3 } \end{array} \right]
\quad L^{T} = \left[ \begin{array} { r r r } { 2 } & { 6 } & { -8 } \\ { 0 } & { 1 } & { 5 } \\ { 0 } & { 0 } & { 3 } \end{array} \right] 
  $$

## QR factorization:

THEOREM $3.1 .$ QR decomposition. Let $A$ be $m-b y-n$ with $m \geq n .$ Suppose that $A$ has full column rank. Then there exists a unique $m-b y-n$ orthogonal matrix $Q\left(Q^{T} Q=I_{n}\right)$ and a unique $n-b y-n$ upper triangular matrix $R$ with positive diagonals $r_{i i}>0$ such that $A=Q R$

### QR Factorization Solves:
- linear equations
- least squares problems
- constrained least squares problems

### Algorithms for QR factorization:

- **Gram-Schmidt** algorithm
  - complexity is $2 m n^{2}$ flops
  - not recommended in practice (sensitive to rounding errors)

- **Householder** algorithm
  - complexity is $2 m n^{2}-(2 / 3) n^{3}$ flops
  - represents $Q$ as a product of elementary orthogonal matrices
  - the most widely used algorithm

### Gram-Schmidt algorithm
Given: $m \times n$ matrix $A$ with linearly independent columns $a_{1}, \ldots, a_{n}$
$$
  \text{Algorithm} \\
  \quad \text{for } k=1 \text{ to } n:
  \quad \begin{aligned}
  \tilde{q}_{1} &= a_{1} \\
  R_{1 k} &=\left\|\tilde{q}_{1}\right\|  \\ 
  q_{1} &=\frac{1}{R_{1 k}}\tilde{q}_{1} \\
  & \vdots \\ 
  R_{k-1, k} &=q_{k-1}^{T} a_{k} \\ 
  \tilde{q}_{k} &=a_{k}-\left(R_{1 k} q_{1}+R_{2 k} q_{2}+\cdots+R_{k-1, k} q_{k-1}\right) \\ 
  R_{k k} &=\left\|\tilde{q}_{k}\right\| \\ 
  q_{k} &=\frac{1}{R_{k k}} \tilde{q}_{k} 
  \end{aligned}
$$

In MATLAB:

```matlab
[m, n] = size(A);
Q = zeros(m,n);
R = zeros(n,n);
for k = 1:n
    R(1:k-1,k) = Q(:,1:k-1)’ * A(:,k);
    v = A(:,k) - Q(:,1:k-1) * R(1:k-1,k);
    R(k,k) = norm(v);
    Q(:,k) = v / R(k,k);
end;
```

![](https://cdn.mathpix.com/snip/images/4FeeEU8Jp1fZcEU4Odpf_u2iL_dlX1kKhjDnlkA9VU4.original.fullsize.png)

### Householder algorithm

TODO - [LINK](http://www.seas.ucla.edu/~vandenbe/133A/lectures/qr.pdf)

## Hessenberg Factorization

The Hessenberg form of $A$ is a factorization $A=Q_{0} \cdot H \cdot Q_{0}^{H}$ where $Q_{0}$ is an orthogonal matrix
and $H$ is upper Hessenberg. How can you compute it in the case when $A$ is a $3 \times 3$ matrix?

A matrix $H=\left(h_{i, j}\right): i, j \in \mathbb{R}^{n \times n}$ is upper *Hessenberg* if all its coefficients below the lower secondary diagonal are zero, that is, if $h_{i, j}=0$ whenever $i \geq j+2$

### Algorithm

- Choose $Q_{1}$ so
  $$
  Q_{1} A=\left[\begin{array}{ccccc}
  {x} & {x} & {x} & {x} & {x} \\
  {x} & {x} & {x} & {x} & {x} \\
  {0} & {x} & {x} & {x} & {x} \\
  {0} & {x} & {x} & {x} & {x} \\
  {0} & {x} & {x} & {x} & {x}
  \end{array}\right] \quad \text { and } \quad A_{1} \equiv Q_{1} A Q_{1}^{T}=\left[\begin{array}{ccccc}
  {x} & {x} & {x} & {x} & {x} \\
  {x} & {x} & {x} & {x} & {x} \\
  {0} & {x} & {x} & {x} & {x} \\
  {0} & {x} & {x} & {x} & {x} \\
  {0} & {x} & {x} & {x} & {x}
  \end{array}\right]
  $$
  $Q_{1}$ leaves the first row of $Q_{1} A$ unchanged, and $Q_{1}^{T}$ leaves the first column of $Q_{1} A Q_{1}^{T}$ unchanged, including the zeros.

- Choose $Q_{2}$ so
  $$
  Q_{2} A_{1}=\left[\begin{array}{ccccc}
  {x} & {x} & {x} & {x} & {x} \\
  {x} & {x} & {x} & {x} & {x} \\
  {0} & {x} & {x} & {x} & {x} \\
  {0} & {0} & {x} & {x} & {x} \\
  {0} & {0} & {x} & {x} & {x}
  \end{array}\right] \text { and } A_{2} \equiv Q_{2} A_{1} Q_{2}^{T}=\left[\begin{array}{ccccc}
  {x} & {x} & {x} & {x} & {x} \\
  {x} & {x} & {x} & {x} & {x} \\
  {0} & {x} & {x} & {x} & {x} \\
  {0} & {0} & {x} & {x} & {x} \\
  {0} & {0} & {x} & {x} & {x}
  \end{array}\right]
  $$
  $Q_{2}$ changes only the last three rows of $A_{1},$ and $Q_{2}^{T}$ leaves the first two columns of $Q_{2} A_{1} Q_{2}^{T}$ unchanged, including the zeros.

- Choose $Q_{3}$ so
  $$
  Q_{3} A_{2}=\left[\begin{array}{ccccc}
  {x} & {x} & {x} & {x} & {x} \\
  {x} & {x} & {x} & {x} & {x} \\
  {o} & {x} & {x} & {x} & {x} \\
  {0} & {o} & {x} & {x} & {x} \\
  {o} & {o} & {o} & {x} & {x}
  \end{array}\right] \quad \text { and } \quad A_{3}=Q_{3} A_{2} Q_{3}^{T}=\left[\begin{array}{ccccc}
  {x} & {x} & {x} & {x} & {x} \\
  {x} & {x} & {x} & {x} & {x} \\
  {0} & {x} & {x} & {x} & {x} \\
  {0} & {0} & {x} & {x} & {x} \\
  {0} & {0} & {0} & {x} & {x}
  \end{array}\right]
  $$
  which is upper Hessenberg. Altogether $A_{3}=\left(Q_{3} Q_{2} Q_{1}\right) \cdot A\left(Q_{3} Q_{2} Q_{1}\right)^{T} \equiv$ $Q A Q^{T}$

ALGORITHM $4.6$ Reduction to upper Hessenberg form:
$$
  \text{if } Q \text{ is desired, set }Q=I\\
$$
$$
  \begin{aligned} 
  \text {for } i &=1: n-2 \\ 
  u_{i} &= \text { House }(A(i+1: n, i)) \\ 
  P_{i} &=I-2 u_{i} u_{i}^{T} & /^{*} Q_{i}=\operatorname{diag}\left(I^{i \times i}, P_{i}\right) * / \\ 
  A(i+1: n, i: n)&=P_{i} \cdot A(i+1: n, i: n) \\ 
  A(1: n, i+1: n)&=A(1: n, i+1: n) \cdot P_{i} \\
  \text{if } Q \text{ is desired} & \\
  Q(i+1: n, i: n) &=P_{i} \cdot Q(i+1: n, i: n) \quad /^{*} Q=Q_{i} \cdot Q^{*} \\
  \quad\text{end if} & \\
  \text{end for} & \\
  \end{aligned}
$$

## Schur Factorization

Schur decomposition allows one to write an arbitrary matrix as unitarily equivalent to an upper triangular matrix whose diagonal elements are the eigenvalues of the original matrix. The Schur decomposition reads as follows: if $A$ is a $n \times n$ square matrix with complex entries, then $A$ can be expressed as:

$$
A=Q U Q^{-1}
$$

- $Q$ is a unitary matrix (so that its inverse $Q^{-1}$ is also the conjugate transpose $Q^{*}$ of $Q$ ), and 

- $U$ is an upper triangular matrix, which is called a **Schur form** of $A$.

- If $A$ is real, its Schur factorization is $A =QUQ^{T}$

- columns of $Q$ are called the **Schur vectors**
- he eigenvalues of $A$ appear on the diagonal of $U$



The first step in a Schur decomposition is a [Hessenberg decomposition](http://mathworld.wolfram.com/HessenbergDecomposition.html).

For example, the Schur decomposition of the matrix
$$
A=\left[\begin{array}{lll}
{3} & {2} & {1} \\
{4} & {2} & {1} \\
{4} & {4} & {0}
\end{array}\right]
$$
is
$$
Q=\left[\begin{array}{ccc}
{0.49857} & {0.76469} & {0.40825} \\ 
{0.57405} & {0.061628} & {-0.81650} \\ 
{0.64953} & {-0.64144} & {0.40825}
\end{array}\right]\\
T=\left[\begin{array}{ccc}
{6.6056} & {4.4907} & {-0.82632} \\
{0.00000} & {-0.60555} & {1.0726} \\
{0.00000} & {0.00000} & {-1.00000}
\end{array}\right]
$$
and the eigenvalues of $T$ are:
$\lambda_{1}=-1, \lambda_{2}=3-\sqrt{13}=-0.60555 \ldots . \lambda_{3}=3+\sqrt{13}=6.6055 \ldots$


# Singular value decomposition (SVD)

# Jacobi method (iterative scheme)

# Gauss-Seidel method (iterative scheme)

# SOR(w) method

**Lecture 1** *(September 13th)***:** **[D]** **1.2** Standard Problems of Numerical Linear Algebra; **1.3** General Techniques.

**Lecture 2** *(September 18th)***: [D]** **2.3** Gaussian Elimination.

**Lecture 3** *(October 2nd)***: [D]** **2.3** Gaussian Elimination (cont.); **1.7** Vector and Matrix Norms; **2.2** Perturbation Theory.

**Lecture 4** *(October 9th)***: [D] 2****.4** Error Analysis; **2.7** Special Linear Systems.

**Lecture 5** *(October 23th)***: [D]** **2.7** Special Linear Systems (cont.).

**Lecture 6** *(October 30th)***: [D]** Linear Least Squares Problems: 3.1 Introduction, 3.2.1 Normal Equations, 3.2.2 QR Decomposition.

**Lecture 7** *(November 6th)***: [D] 3.4** Orthogonal Matrices**, 3.4.1** Householder transformations, **3.4.2** Givens rotations.

**Lecture 8** *(November 13th)***: [D]** Singular Value Decomposition: **3.2.3** SVD and the Least Square Problem, **[S]** Principal Component Analysis.

**Lecture 9** *(November 20th)*: **[BL]** Google's PageRank algorithm.

**Lecture 10** *(November 27th)*: **[D]** Nonsymmetric Eigenvalue Problems: 4.1 Introduction; Algorithms for the Nonsymmetric Eigenproblem: 4.4.1 Power Method; **[E]** Eigenvalues of the Pagerank matrix.

**Lecture 11** *(November 29th)*: **[D]** Nonsymmetric Eigenvalue Problems: 4.2 Canonical Forms, 4.3 Computing Eigenvectors from the Schur Form; Algorithms for the Nonsymmetric Eigenproblem: 4.4.2 Inverse Iteration, 4.4.3 Orthogonal Iteration, 4.4.4 QR Iteration.

**Lecture 12** *(December 4th)*: **[D]** Iterative Methods for Linear Systems: 6.5 Basic Iterative Methods, 6.5.1 Jacobi's Method, 6.5.2 Gauss Seidel Method, 6.5.3 Successive Overrelaxation.

**Lecture 13** *(December 11th)*: [D] Computing the SVD: 5.4 Algorithms for the Singular Value Decomposition, 4.4.7 Tridiagonal and Bidiagonal Reduction, 5.4.1 QR Iteration and Its Variations for the Bidiagonal SVD

**Bibliography**

**[BL]** K. Bryan and T. Leise, *The $25,000,000,000 eigenvector: The linear algebra behind Google*, SIAM Review 48 (2006) 569–581.



**[D]** J. Demmel, *Applied Numerical Linear Algebra*, SIAM 1997.

**[E]** L. Elden, A note on the eigenvalues of the Google matrix, e-print arXiv:math/0401177 (2004), 2 pages.

**[S]** J. Shlens, *A tutorial on principal components analysis*, e-print arXiv:1404:1100 (2014), 12 pages
