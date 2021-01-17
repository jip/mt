require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a non-complex
NB.   symmetric matrix using the factorization computed by
NB.   DSYTRF_AA
NB.
NB. Syntax:
NB.   X=. uplo dsytrs_aa DT1 ; ipiv ; B
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of DT1 only, form is:
NB.                    P * L1 * D * L1^T * P^T = A
NB.            'U' - use upper triangle of DT1 only, form is:
NB.                    P * U1 * D * U1^T * P^T = A
NB.   DT1  - n×n-matrix, D and T1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          P
NB.   B    - n×nrhs-matrix, real, RHS
NB.   X    - n×nrhs-matrix, real, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, real, the symmetric, represented in
NB.          the factored form by DT1 and ipiv
NB.   D    - n×n-matrix, real, the symmetric tridiagonal
NB.          (opposite diagonal not stored)
NB.   T1   - n×n-matrix, either L1 or U1
NB.   L1   - n×n-matrix, real, the unit lower triangular
NB.          (unit diagonal not stored)
NB.   U1   - n×n-matrix, real, the unit upper triangular
NB.          (unit diagonal not stored)
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dsytrs_aa=: 4 : 0
  'DT1 ipiv B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_ *. n = #) DT1
  assert. (isvector_jlapack2_ *.                       (-: <.) :: 0: *.    n = #) ipiv
  assert. (ismatrix_jlapack2_ *.                       isreal_jlapack2_         ) B
  select. 3!:0 DT1
    case. JCMPX do. DT1=. 9 o. DT1
    case. JFL   do.
    case.       do. DT1=. DT1 + 0.0
  end.
  select. 3!:0 B
    case. JCMPX do. B=. 9 o. B
    case. JFL   do.
    case.       do. B=. B + 0.0
  end.
  if. JINT ~: 3!:0 ipiv do. ipiv=. <. 9 o. ipiv end.
  ld=. , 1 >. n
  lwork=. , 1 >. _2 3 p. n  NB. minimal
  cdrc=. dsytrs_aa_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: DT1) ; ld ; ipiv ; (|: B) ; ld ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 7 {:: cdrc
)
