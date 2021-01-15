require 'math/lapack2'

NB. Description:
NB.   Compute the inverse of a factored symmetric matrix
NB.
NB. Syntax:
NB.   iAA=. uplo dsytri2 DPT1 ; ipiv
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of A only, form is:
NB.                    PL1 * D * PL1^T = A
NB.            'U' - use upper triangle of A only, form is:
NB.                    PU1 * D * PU1^T = A
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          permutations
NB.   D    - n×n-matrix, the symmetric and block diagonal
NB.          with 1×1 and 2×2 diagonal blocks (opposite
NB.          diagonal not stored)
NB.   PT1  - n×n-matrix, either PL1 or PU1
NB.   PL1  - n×n-matrix, the product of permutation and unit
NB.          lower triangular matrices (unit diagonal not
NB.          stored)
NB.   PU1  - n×n-matrix, the product of permutation and unit
NB.          upper triangular matrices (unit diagonal not
NB.          stored)
NB.   iAA  - n×n-matrix, contains the triangular part of iA
NB.          in changed triangle and unchanged elements of
NB.          DPT1 in opposite strict triangle
NB.   iA   - n×n-matrix, the symmetric inversion of A
NB.   A    - n×n-matrix, the symmetric represented in factored
NB.          form by DPT1 and ipiv
NB.   n    ≥ 0, the size of A and iA
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dsytri2=: 4 : 0
  'DPT1 ipiv'=. y
  assert. *./ (1 = # x) , x e. 'lLuU'
  n=. # ipiv
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_ *. n = #) DPT1
  assert. (isvector_jlapack2_ *.                       (-: <.) :: 0:            ) ipiv
  select. 3!:0 DPT1
    case. JCMPX do. DPT1=. 9 o. DPT1
    case. JFL   do.
    case.       do. DPT1=. DPT1 + 0.0
  end.
  if. JINT  ~: 3!:0 ipiv do. ipiv=. <. 9 o. ipiv end.
  lwork=. , (n + 64 + 1) * (64 + 3)  NB. minimal
  cdrc=. dsytri2_jlapack2_ (, x) ; (, n) ; (|: DPT1) ; (, 1 >. n) ; ipiv ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 3 {:: cdrc
)
