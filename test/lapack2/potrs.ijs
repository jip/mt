require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   (symmetric) positive definite matrix using the
NB.   factorization computed by xPOTRF
NB.
NB. Syntax:
NB.   X=. uplo xpotrs T ; B
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of A only, form is:
NB.                    L * L^H = A
NB.            'U' - use upper triangle of A only, form is:
NB.                    U^H * U = A
NB.   T    - n×n-matrix, L or U factor
NB.   B    - n×nrhs-matrix, RHS
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite represented in the factored form by T
NB.   L    - n×n-matrix, the lower triangular
NB.   U    - n×n-matrix, the upper triangular
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale
NB. - no check for positive definiteness

dpotrs=: 4 : 0
  'T B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_ *. n = #) T
  assert. (ismatrix_jlapack2_ *.                       isreal_jlapack2_         ) B
  select. 3!:0 T
    case. JCMPX do. T=. 9 o. T
    case. JFL   do.
    case.       do. T=. T + 0.0
  end.
  select. 3!:0 B
    case. JCMPX do. B=. 9 o. B
    case. JFL   do.
    case.       do. B=. B + 0.0
  end.
  ld=. , 1 >. n
  cdrc=. dpotrs_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: T) ; ld ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  |: 6 {:: cdrc
)

zpotrs=: 4 : 0
  'T B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. n = #) T
  assert.  ismatrix_jlapack2_                                 B
  if. JCMPX ~: 3!:0 T do. T=. T + 0j0 end.
  if. JCMPX ~: 3!:0 B do. B=. B + 0j0 end.
  ld=. , 1 >. n
  cdrc=. zpotrs_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: T) ; ld ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  |: 6 {:: cdrc
)
