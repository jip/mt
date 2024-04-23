require 'math/lapack2'

NB. Description:
NB.   Compute the LU factorization of a general matrix using
NB.   partial pivoting with row interchanges
NB.
NB. Syntax:
NB.   'L1U ipiv'=. xgetrf A
NB. where
NB.   A    - m×n-matrix to be factored as:
NB.            P * L1 * U = A
NB.   L1U  - m×n-matrix, L1 and U combined
NB.   ipiv - k-vector, integer, pivot indices that define P
NB.   L1   - m×k-matrix, unit lower trapezoidal (unit
NB.          diagonal is not stored)
NB.   U    - k×n-matrix, upper trapezoidal
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   m    ≥ 0, the number of rows in A and L1
NB.   n    ≥ 0, the number of columns in A and U
NB.   k    = min(m,n)
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dgetrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  (|:L:0) 3 5 { dgetrf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 00) ; , _1
    NB. (|:) doesn't affect to ipiv
)

zgetrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  (|:L:0) 3 5 { zgetrf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 00) ; , _1
    NB. (|:) doesn't affect to ipiv
)
