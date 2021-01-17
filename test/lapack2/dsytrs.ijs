require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a non-complex
NB.   symmetric matrix using the factorization computed by
NB.   DSYTRF
NB.
NB. Syntax:
NB.   X=. uplo dsytrs DPT1 ; ipiv ; B
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of DPT1 only, form
NB.                  is:
NB.                    PL1 * D * PL1^T = A
NB.            'U' - use upper triangle of DPT1 only, form
NB.                  is:
NB.                    PU1 * D * PU1^T = A
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          permutations
NB.   B    - n×nrhs-matrix, real, RHS
NB.   X    - n×nrhs-matrix, real, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, real, the symmetric, represented in
NB.          the factored form by DPT1 and ipiv
NB.   D    - n×n-matrix, real, the symmetric and block
NB.          diagonal with 1×1 and 2×2 diagonal blocks
NB.          (opposite diagonal not stored)
NB.   PT1  - n×n-matrix, either PL1 or PU1
NB.   PL1  - n×n-matrix, real, the product of permutation and
NB.          unit lower triangular matrices
NB.   PU1  - n×n-matrix, real, the product of permutation and
NB.          unit upper triangular matrices
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dsytrs=: 4 : 0
  'DPT1 ipiv B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_ *. n = #) DPT1
  assert. (isvector_jlapack2_                       *. (-: <.) :: 0:    *. n = #) ipiv
  assert. (ismatrix_jlapack2_ *.                       isreal_jlapack2_         ) B
  select. 3!:0 DPT1
    case. JCMPX do. DPT1=. 9 o. DPT1
    case. JFL   do.
    case.       do. DPT1=. DPT1 + 0.0
  end.
  select. 3!:0 B
    case. JCMPX do. B=. 9 o. B
    case. JFL   do.
    case.       do. B=. B + 0.0
  end.
  if. JINT ~: 3!:0 ipiv do. ipiv=. <. 9 o. ipiv end.
  ld=. , 1 >. n
  cdrc=. dsytrs_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: DPT1) ; ld ; ipiv ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  |: 7 {:: cdrc
)
