require 'math/lapack2'

NB. Description:
NB.   Compute the factorization of a non-complex symmetric
NB.   matrix using the Aasen's algorithm
NB.
NB. Syntax:
NB.   'DT1 ipiv'=. uplo dsytrf_aa A
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.          specifies which triangular part of A is to be
NB.          referenced:
NB.            'L' - lower, the form is:
NB.                    P * L1 * D * L1^T * P^T = A
NB.            'U' - upper, the form is:
NB.                    P * U1 * D * U1^T * P^T = A
NB.   A    - n×n-matrix, real, symmetric or upper or lower
NB.          triangular
NB.   DT1  - n×n-matrix, D and T1 combined
NB.   ipiv - n-vector, integer, pivot indices that define P
NB.   D    - n×n-matrix, real, symmetric tridiagonal
NB.          (opposite diagonal is not stored)
NB.   T1   - n×n-matrix, either L1 or U1
NB.   L1   - n×n-matrix, real, unit lower triangular (unit
NB.          diagonal not stored)
NB.   U1   - n×n-matrix, real, unit upper triangular (unit
NB.          diagonal not stored)
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   n    ≥ 0, the size of A, D, L1, U1 and P
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dsytrf_aa=: 4 : 0
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  NB. lwork=. , 1 >. +: n  NB. minimal
  lwork=. , 1 >. n * >: 64  NB. optimal
  (|: L: 0) 3 5 { dsytrf_aa_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 00) ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to ipiv
)
