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
NB.   A    - n×n-matrix, real, the symmetric or upper or
NB.          lower triangular
NB.   DT1  - n×n-matrix, D and T1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          P
NB.   D    - n×n-matrix, real, the symmetric tridiagonal
NB.          (opposite diagonal not stored)
NB.   T1   - n×n-matrix, either L1 or U1
NB.   L1   - n×n-matrix, real, the unit lower triangular
NB.          (unit diagonal not stored)
NB.   U1   - n×n-matrix, real, the unit upper triangular
NB.          (unit diagonal not stored)
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   n    ≥ 0, the size of A, D, L1, U1 and P
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dsytrf_aa=: 4 : 0
  assert. 'lLuU' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , isreal_jlapack2_) y
  n=. # y
  NB. lwork=. , 1 >. +: n  NB. minimal
  lwork=. , 1 >. n * >: 64  NB. optimal
  cdrc=. dsytrf_aa_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 00) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to ipiv
)
