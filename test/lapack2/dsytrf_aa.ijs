require 'math/lapack2'

NB. Description:
NB.   Compute the factorization of a non-complex symmetric
NB.   matrix using the Aasen's algorithm
NB.
NB. Syntax:
NB.   'DT1 ipiv'=. uplo dsytrf_aa A
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of A only, form is:
NB.                    P * L1 * D * L1^T * P^T = A
NB.            'U' - use upper triangle of A only, form is:
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
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  n=. # y
  NB. lwork=. , 1 >. +: n  NB. minimal
  lwork=. , 1 >. n * >: 64  NB. optimal
  cdrc=. dsytrf_aa_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 00) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to ipiv
)
