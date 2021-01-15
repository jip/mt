require 'math/lapack2'

NB. Description:
NB.   Compute the factorization of a non-complex symmetric
NB.   matrix using the Bunch-Kaufman diagonal pivoting method
NB.
NB. Syntax:
NB.   'DPT1 ipiv'=. uplo dsytrf A
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of A only, form is:
NB.                    PL1 * D * PL1^T = A
NB.            'U' - use upper triangle of A only, form is:
NB.                    PU1 * D * PU1^T = A
NB.   A    - n×n-matrix, real, the symmetric or upper or
NB.          lower triangular
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          permutations
NB.   D    - n×n-matrix, real, the symmetric and block
NB.          diagonal with 1×1 and 2×2 diagonal blocks
NB.          (opposite diagonal not stored)
NB.   PT1  - n×n-matrix, either PL1 or PU1
NB.   PL1  - n×n-matrix, real, the product of permutation and
NB.          unit lower triangular matrices
NB.   PU1  - n×n-matrix, real, the product of permutation and
NB.          unit upper triangular matrices
NB.   n    ≥ 0, the size of A, D, PL1 and PU1
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dsytrf=: 4 : 0
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  n=. # y
  NB. lwork=. , 1  NB. minimal
  lwork=. , 1 >. n * 64  NB. optimal
  cdrc=. dsytrf_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 00) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to ipiv
)
