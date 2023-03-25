require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   (symmetric) positive definite matrix, and factor the
NB.   last one
NB.
NB. Syntax:
NB.   'T X'=. uplo xposv A ; B
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.           specifies which triangular part of A is to be
NB.           referenced:
NB.            'L' - lower, the form is:
NB.                    L * L^H = A
NB.            'U' - upper, the form is:
NB.                    U^H * U = A
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite to be factored to T
NB.   B    - n×nrhs-matrix, RHS
NB.   T    - n×n-matrix, L or U factor
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   L    - n×n-matrix, the lower triangular
NB.   U    - n×n-matrix, the upper triangular
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale
NB. - no check for positive definiteness

dposv=: 4 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) A
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  cdrc=. dposv_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 4 6 { cdrc
)

zposv=: 4 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) A
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  cdrc=. zposv_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 4 6 { cdrc
)
