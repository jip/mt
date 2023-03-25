require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   (symmetric) positive definite matrix using the
NB.   factorization computed by xPOTRF
NB.
NB. Syntax:
NB.   X=. uplo xpotrs T ; B
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.           specifies which triangular part of A is to be
NB.           referenced:
NB.            'L' - lower, the form is:
NB.                    L * L^H = A
NB.            'U' - upper, the form is:
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
  assert. 'lLuU' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) T
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  cdrc=. dpotrs_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: T) ; ld ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  |: 6 {:: cdrc
)

zpotrs=: 4 : 0
  'T B'=. y
  'n nrhs'=. $ B
  assert. 'lLuU' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) T
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  cdrc=. zpotrs_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: T) ; ld ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  |: 6 {:: cdrc
)
