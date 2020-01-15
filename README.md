# Problems
<p align="center"><img src="/tex/607d9ff0b7cbe7286243d8abea3971d3.svg?invert_in_darkmode&sanitize=true" align=middle width=0.0pt height=0.0pt/></p>

- Linear System: <img src="/tex/70681e99f542745bf6a0c56bd4600b39.svg?invert_in_darkmode&sanitize=true" align=middle width=50.69621369999999pt height=22.831056599999986pt/>
- Least Squares Problem <img src="/tex/58230adccafdd16765dfcfdcd5a78f0f.svg?invert_in_darkmode&sanitize=true" align=middle width=73.68721634999999pt height=24.65753399999998pt/>
- Eigenvalue problem <img src="/tex/4cc103f3eeb0d48b1c4d816d4ca52803.svg?invert_in_darkmode&sanitize=true" align=middle width=62.62548764999999pt height=22.831056599999986pt/>, where <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/> - vector <img src="/tex/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode&sanitize=true" align=middle width=9.58908224999999pt height=22.831056599999986pt/> - scalar
- Singular Value problem <img src="/tex/0fcfbd523dff9f2e0696b4493cf6da58.svg?invert_in_darkmode&sanitize=true" align=middle width=85.30987904999999pt height=27.6567522pt/>


# Definitions

- **Nonsingular matrix** - A square matrix that is not singular, i.e. one that has a matrix inverse. Nonsingular matrices are sometimes also called regular matrices. A square matrix is nonsingular if its determinant is nonzero

- **Inversing Triangular Matrix**

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
    giving
    $$
    \left(\begin{array}{lll}
    {1 / a} & {-b /(a d)} & {(b e-c d) /(a f d)} \\
    {0} & {1 / d} & {-e /(f d)} \\
    {0} & {0} & {1 / f}
    \end{array}\right)
    $$
  from which we see directly that the matrix is invertible if all <img src="/tex/7704b8b9e1b077e08a7aabdcbd2bde61.svg?invert_in_darkmode&sanitize=true" align=middle width=24.55100009999999pt height=22.831056599999986pt/> and <img src="/tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode&sanitize=true" align=middle width=9.81741584999999pt height=22.831056599999986pt/> are different from zero

- **Positive Definite matrix** - <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> <img src="/tex/3add1221abfa79cb14021bc2dacd5725.svg?invert_in_darkmode&sanitize=true" align=middle width=39.82494449999999pt height=19.1781018pt/> symmetric real matrix <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/> is said to be positive definite if <img src="/tex/79dc7e289f38476fbcb56fae02e83ccd.svg?invert_in_darkmode&sanitize=true" align=middle width=77.76244739999999pt height=27.91243950000002pt/> for all non-zero <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/> in <img src="/tex/3afa189299419b087ab3ae25810cd215.svg?invert_in_darkmode&sanitize=true" align=middle width=19.998202949999992pt height=22.648391699999998pt/>.

  > <img src="/tex/e77c2ab4d1cb17cc5a7788e7ba3e73a6.svg?invert_in_darkmode&sanitize=true" align=middle width=379.7356562999999pt height=27.91243950000002pt/>

- **Positive Semidefinite matrix** - <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> <img src="/tex/3add1221abfa79cb14021bc2dacd5725.svg?invert_in_darkmode&sanitize=true" align=middle width=39.82494449999999pt height=19.1781018pt/> symmetric real matrix <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/> is said to be positive semidefinite or non-negative definite if <img src="/tex/7155bcffbecd42c3e39451a756f0e224.svg?invert_in_darkmode&sanitize=true" align=middle width=77.76244739999999pt height=27.91243950000002pt/> for all <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/> in <img src="/tex/444f1a00d4adf380df349d93dbe734bc.svg?invert_in_darkmode&sanitize=true" align=middle width=25.38634064999999pt height=22.648391699999998pt/> Formally,

  > <img src="/tex/05237438f57a408e67955e466a1ba516.svg?invert_in_darkmode&sanitize=true" align=middle width=362.47532309999997pt height=27.91243950000002pt/>

- **Rank of a matrix** - The rank of a matrix is defined as 

  - the maximum number of linearly independent column vectors in the matrix or 

  - the maximum number of linearly independent row vectors in the matrix. 

    Both definitions are equivalent.

- **Eigenvalues and eigenvectors**- Eigenvalues (*characteristic roots, proper values, latent roots)* are a special set of scalars associated with a linear system of equations. Each eigenvalue is paired with a corresponding so-called eigenvector. 

  Let <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> be a linear transformation represented by a matrix A. If there is a vector <img src="/tex/396f39b6dd4c8d6ac53c15f787027552.svg?invert_in_darkmode&sanitize=true" align=middle width=85.34026049999999pt height=22.831056599999986pt/> such that:

  <img src="/tex/f3bb70f791561eade25e2f232027d264.svg?invert_in_darkmode&sanitize=true" align=middle width=73.65286664999998pt height=22.831056599999986pt/>

  for some scalar <img src="/tex/722d262217f99de10f2a57fe9cb64338.svg?invert_in_darkmode&sanitize=true" align=middle width=14.155307099999991pt height=22.831056599999986pt/> then <img src="/tex/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode&sanitize=true" align=middle width=9.58908224999999pt height=22.831056599999986pt/> is called the eigenvalue of A with corresponding (right) eigenvector <img src="/tex/d05b996d2c08252f77613c25205a0f04.svg?invert_in_darkmode&sanitize=true" align=middle width=14.29216634999999pt height=22.55708729999998pt/>.

  Letting A be a <img src="/tex/1bd1a25ad3542becf0e5739e9fb42d2f.svg?invert_in_darkmode&sanitize=true" align=middle width=38.24192129999999pt height=22.831056599999986pt/> square matrix
  <p align="center"><img src="/tex/b5a0b03d2dc92d81704b5a59ce7ead0f.svg?invert_in_darkmode&sanitize=true" align=middle width=177.55163249999998pt height=88.76800184999999pt/></p>
  with eigenvalue <img src="/tex/722d262217f99de10f2a57fe9cb64338.svg?invert_in_darkmode&sanitize=true" align=middle width=14.155307099999991pt height=22.831056599999986pt/> then the corresponding eigenvectors satisfy
  <p align="center"><img src="/tex/25724d24c06e752352d243020ba84e20.svg?invert_in_darkmode&sanitize=true" align=middle width=326.2159461pt height=88.76800184999999pt/></p>
which is equivalent to the homogeneous system
  <p align="center"><img src="/tex/a756317b45412de29bc75eaef79d36fe.svg?invert_in_darkmode&sanitize=true" align=middle width=392.23733174999995pt height=88.76800184999999pt/></p>
  Equation ( 4) can be written compactly as
  <p align="center"><img src="/tex/b839829210d7ca234a971a5eabd49b05.svg?invert_in_darkmode&sanitize=true" align=middle width=108.35597849999998pt height=16.438356pt/></p>
  
- **Singular Value**

- **Condition number**

- **Norms**
  <p align="center"><img src="/tex/db2a2580836f57ddd04e9e5c401f93ab.svg?invert_in_darkmode&sanitize=true" align=middle width=152.7998373pt height=78.703779pt/></p>

  <p align="center"><img src="/tex/35f4df1d8ea027c65275274ea85fdecf.svg?invert_in_darkmode&sanitize=true" align=middle width=309.42971895pt height=117.1525905pt/></p>
  
- Rank1/2 approximation

# Matrix Factorization:

A factorization of the matrix <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> is a representation of <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> as a product of several "simpler" matrices, which make the problem at hand easier to solve. We give two examples.
<p align="center"><img src="/tex/5309fbf7b6112c3fb2fab4793820496a.svg?invert_in_darkmode&sanitize=true" align=middle width=316.7069796pt height=88.76800184999999pt/></p>
- Forward Substitution

  for <img src="/tex/081b3265e50f611c00daeffa91931873.svg?invert_in_darkmode&sanitize=true" align=middle width=35.80006649999999pt height=21.68300969999999pt/> to <img src="/tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode&sanitize=true" align=middle width=9.86687624999999pt height=14.15524440000002pt/>
  <p align="center"><img src="/tex/1eabfd12ac210a4d1a5e393fb27ad2f9.svg?invert_in_darkmode&sanitize=true" align=middle width=189.94192139999998pt height=49.315569599999996pt/></p>
  

end for

- Backward Substitution

  An analogous idea, back substitution, works if <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> is upper triangular. 
To use this to solve a general system <img src="/tex/594efbce3e38e300428e85e0bdae6d07.svg?invert_in_darkmode&sanitize=true" align=middle width=50.69621369999999pt height=22.831056599999986pt/> we need the following matrix factorization, which is just a restatement of Gaussian elimination.

## LU Factorization
<p align="center"><img src="/tex/ce2b2f16485db094ca7afcbf5d48e8af.svg?invert_in_darkmode&sanitize=true" align=middle width=379.437696pt height=59.1786591pt/></p>
**2** and **-3** goes down to <img src="/tex/ddcb483302ed36a59286424aa5e0be17.svg?invert_in_darkmode&sanitize=true" align=middle width=11.18724254999999pt height=22.465723500000017pt/> with changed sign
<p align="center"><img src="/tex/eab43629b44054198bbc9c76abd17b9f.svg?invert_in_darkmode&sanitize=true" align=middle width=401.3120463pt height=59.1786591pt/></p>
  **0.5** goes to <img src="/tex/ddcb483302ed36a59286424aa5e0be17.svg?invert_in_darkmode&sanitize=true" align=middle width=11.18724254999999pt height=22.465723500000017pt/> with minus sign.
<p align="center"><img src="/tex/5c5dacc34cf2bc4804fcdd2422f1632e.svg?invert_in_darkmode&sanitize=true" align=middle width=346.5777084pt height=59.1786591pt/></p>

<p align="center"><img src="/tex/25f9864740ce9813725ed11ea1dea5cd.svg?invert_in_darkmode&sanitize=true" align=middle width=58.44963299999999pt height=11.232861749999998pt/></p>
## Gaussian Elimination

**Partial Pivoting** - Sort rows by its first elements absolute value. Choose <img src="/tex/e84705df9be792193f6e243da302768c.svg?invert_in_darkmode&sanitize=true" align=middle width=24.51021539999999pt height=14.15524440000002pt/> and swap that  <img src="/tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode&sanitize=true" align=middle width=9.075367949999992pt height=22.831056599999986pt/>-th row with <img src="/tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode&sanitize=true" align=middle width=5.663225699999989pt height=21.68300969999999pt/>-th row `for i in (1,n) `.

**THEOREM** If the <img src="/tex/8c7fe813c58db15b666f150d501557a7.svg?invert_in_darkmode&sanitize=true" align=middle width=75.62013855pt height=22.831056599999986pt/> matrix <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> is *nonsingular*, there exists 

- a permutation matrix <img src="/tex/df5a289587a2f0247a5b97c1e8ac58ca.svg?invert_in_darkmode&sanitize=true" align=middle width=12.83677559999999pt height=22.465723500000017pt/> (the identity matrix with its rows permuted), 
- a nonsingular lower triangular matrix <img src="/tex/d2daf7d13b60e534226525f7dc92de97.svg?invert_in_darkmode&sanitize=true" align=middle width=15.75346574999999pt height=22.465723500000017pt/> and 
- a nonsingular upper triangular matrix <img src="/tex/6bac6ec50c01592407695ef84f457232.svg?invert_in_darkmode&sanitize=true" align=middle width=13.01596064999999pt height=22.465723500000017pt/> 

such that <img src="/tex/53803971358acb1e9e90c324347564ee.svg?invert_in_darkmode&sanitize=true" align=middle width=97.77005534999998pt height=22.465723500000017pt/> To solve <img src="/tex/0e9f8da35c2936ba4468f9f93b549041.svg?invert_in_darkmode&sanitize=true" align=middle width=55.26243689999999pt height=22.831056599999986pt/> we solve the equivalent system <img src="/tex/3fe3b88bb09941f29aac5aca4073d8fa.svg?invert_in_darkmode&sanitize=true" align=middle width=75.40736939999998pt height=22.831056599999986pt/> as follows:

- <img src="/tex/c37437203e3e33eb68c41e294504aa77.svg?invert_in_darkmode&sanitize=true" align=middle width=324.37614329999997pt height=27.94539330000001pt/>
- <img src="/tex/2d2cd1e8c52efbb4779bd385f49a635d.svg?invert_in_darkmode&sanitize=true" align=middle width=304.14248205pt height=27.94539330000001pt/>
- <img src="/tex/8a532061e95bceb81e1aaf3fd728ff9a.svg?invert_in_darkmode&sanitize=true" align=middle width=299.7817251pt height=27.94539330000001pt/>

Important thing <img src="/tex/4c81e1f11c1e21db89203acade86ffb6.svg?invert_in_darkmode&sanitize=true" align=middle width=71.47429245pt height=27.6567522pt/> and  <img src="/tex/51ae9e9184ff89dd7c4b4b0f8cecb0ae.svg?invert_in_darkmode&sanitize=true" align=middle width=75.13168574999999pt height=27.6567522pt/> - check **inversion of triangular matrix** in section definitions.

## Cholesky Factorization

The Cholesky decomposition is roughly twice as efficient as the LU decomposition for solving systems of linear equations. <img src="/tex/71a38687d31c242a33a3ec23e25e3ba0.svg?invert_in_darkmode&sanitize=true" align=middle width=14.29216634999999pt height=22.55708729999998pt/> needs to be symmetric. Every symmetric, positive definite matrix A can be decomposed into a product of a unique lower triangular matrix L and its transpose. <img src="/tex/95a6c12ed4f3720ff49aeb54580956e7.svg?invert_in_darkmode&sanitize=true" align=middle width=66.15463304999999pt height=27.6567522pt/>, where <img src="/tex/ddcb483302ed36a59286424aa5e0be17.svg?invert_in_darkmode&sanitize=true" align=middle width=11.18724254999999pt height=22.465723500000017pt/> is a lower triangular matrix with real and positive diagonal entries.

Cholesky Algorithm:
<p align="center"><img src="/tex/46a9dcf838b81fc0d025f29d9ca50f0b.svg?invert_in_darkmode&sanitize=true" align=middle width=406.9916631pt height=180.72714495pt/></p>
**Example**
<p align="center"><img src="/tex/380d1e5c3be7ad281dada632823f9b5f.svg?invert_in_darkmode&sanitize=true" align=middle width=193.15087605pt height=59.1786591pt/></p>
- <img src="/tex/1246fb43a66896c3cdcdfd81957615a8.svg?invert_in_darkmode&sanitize=true" align=middle width=44.93238914999999pt height=21.68300969999999pt/>
<p align="center"><img src="/tex/679fdcd6f8fcd870707846b6e84ee99a.svg?invert_in_darkmode&sanitize=true" align=middle width=276.2605362pt height=98.8387554pt/></p>
- <img src="/tex/0913491f80945049ed061e8767a50fe3.svg?invert_in_darkmode&sanitize=true" align=middle width=44.93238914999999pt height=21.68300969999999pt/>
<p align="center"><img src="/tex/49b01ce4c91b6fb26dc55b4c2b717b67.svg?invert_in_darkmode&sanitize=true" align=middle width=544.1757419999999pt height=59.1786591pt/></p>
- <img src="/tex/b8fe7a8f30cfdea91cd4bb99033d6e21.svg?invert_in_darkmode&sanitize=true" align=middle width=35.80006649999999pt height=21.68300969999999pt/>
<p align="center"><img src="/tex/b6dbc11f79b366f3a10fd2ab594b1728.svg?invert_in_darkmode&sanitize=true" align=middle width=205.68667514999999pt height=111.27665505pt/></p>

<p align="center"><img src="/tex/fa40a0d636fd317180e831bf13d7648f.svg?invert_in_darkmode&sanitize=true" align=middle width=313.09534245pt height=59.1786591pt/></p>
## QR Factorization

Any real square matrix <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> may be decomposed as
<p align="center"><img src="/tex/b31353482d3c7f99aed080c3d6a7a799.svg?invert_in_darkmode&sanitize=true" align=middle width=59.850326249999995pt height=14.42921205pt/></p>
**Definition** with <img src="/tex/1afcdb0f704394b16fe85fb40c45ca7a.svg?invert_in_darkmode&sanitize=true" align=middle width=12.99542474999999pt height=22.465723500000017pt/> as a vector:

where <img src="/tex/1afcdb0f704394b16fe85fb40c45ca7a.svg?invert_in_darkmode&sanitize=true" align=middle width=12.99542474999999pt height=22.465723500000017pt/> is an orthogonal matrix and <img src="/tex/1e438235ef9ec72fc51ac5025516017c.svg?invert_in_darkmode&sanitize=true" align=middle width=12.60847334999999pt height=22.465723500000017pt/> is an upper triangular matrix.
<p align="center"><img src="/tex/6fbd4573d4f97555024f1429d6d99ca1.svg?invert_in_darkmode&sanitize=true" align=middle width=370.54959974999997pt height=88.76800184999999pt/></p>
vectors <img src="/tex/ea4eed098b4de19e35fb622e5df2383a.svg?invert_in_darkmode&sanitize=true" align=middle width=66.70652339999998pt height=14.15524440000002pt/> are orthonormal <img src="/tex/0e51a2dede42189d77627c4d742822c3.svg?invert_in_darkmode&sanitize=true" align=middle width=14.433101099999991pt height=14.15524440000002pt/> -vectors:
<p align="center"><img src="/tex/f4a1a3dcc122e2298c81eaf36bbc575e.svg?invert_in_darkmode&sanitize=true" align=middle width=207.1345716pt height=19.3534671pt/></p>
**QR factorization solves:**

- linear equations
- least squares problems
- constrained least squares problems

**Algorithms for QR factorization:**

- **Gram-Schmidt** algorithm
  - complexity is <img src="/tex/68decd96d08da9efc67888be6609025c.svg?invert_in_darkmode&sanitize=true" align=middle width=39.071732699999984pt height=26.76175259999998pt/> flops
  - not recommended in practice (sensitive to rounding errors)

- **Householder** algorithm
  - complexity is <img src="/tex/c1368d23c8a63878aa92211eda93579b.svg?invert_in_darkmode&sanitize=true" align=middle width=113.84732204999997pt height=26.76175259999998pt/> flops
  - represents <img src="/tex/1afcdb0f704394b16fe85fb40c45ca7a.svg?invert_in_darkmode&sanitize=true" align=middle width=12.99542474999999pt height=22.465723500000017pt/> as a product of elementary orthogonal matrices
  - the most widely used algorithm

**Gram-Schmidt algorithm**
Given: <img src="/tex/205995f88b807b2f5268f7ef4053f049.svg?invert_in_darkmode&sanitize=true" align=middle width=44.39116769999999pt height=19.1781018pt/> matrix <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> with linearly independent columns <img src="/tex/92abd3df4e519ffaaf3007814390ba78.svg?invert_in_darkmode&sanitize=true" align=middle width=69.40820699999999pt height=14.15524440000002pt/>
<p align="center"><img src="/tex/03dcc267d1b8671c4501cae505f5362f.svg?invert_in_darkmode&sanitize=true" align=middle width=584.02414455pt height=239.31701969999997pt/></p>

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

**Householder algorithm**

TODO - [LINK](http://www.seas.ucla.edu/~vandenbe/133A/lectures/qr.pdf)

## Schur Factorization

Real Schur form

Eigenvals and vects

Algorithm

## Hessenberg Factorization

The Hessenberg form of <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> is a factorization <img src="/tex/65bec5691bc298cf90102fa7f5960db3.svg?invert_in_darkmode&sanitize=true" align=middle width=118.00183229999999pt height=27.6567522pt/> where <img src="/tex/a42991bf8a1397710815cddda20e23ad.svg?invert_in_darkmode&sanitize=true" align=middle width=19.547970749999987pt height=22.465723500000017pt/> is an orthogonal matrix
and <img src="/tex/7b9a0316a2fcd7f01cfd556eedf72e96.svg?invert_in_darkmode&sanitize=true" align=middle width=14.99998994999999pt height=22.465723500000017pt/> is upper Hessenberg. How can you compute it in the case when <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> is a <img src="/tex/9f2b6b0a7f3d99fd3f396a1515926eb3.svg?invert_in_darkmode&sanitize=true" align=middle width=36.52961069999999pt height=21.18721440000001pt/> matrix?

A matrix <img src="/tex/599284dd994b567ae189139c514c7fc0.svg?invert_in_darkmode&sanitize=true" align=middle width=167.52284339999997pt height=26.17730939999998pt/> is upper *Hessenberg* if all its coefficients below the lower secondary diagonal are zero, that is, if <img src="/tex/aafb40d83cf7858c39879e3a77c8bf99.svg?invert_in_darkmode&sanitize=true" align=middle width=55.08934694999999pt height=22.831056599999986pt/> whenever <img src="/tex/dac248d24ec5876d365326897e1e56b9.svg?invert_in_darkmode&sanitize=true" align=middle width=63.601675499999985pt height=21.68300969999999pt/>

Algorithm







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
