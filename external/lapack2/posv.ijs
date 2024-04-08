require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   (symmetric) positive definite matrix, and factor the
NB.   last one
NB.
NB. Syntax:
NB.   'T X'=. uplo xposv AA ; B
NB. where
NB.   uplo - string, case-insensitive, in which the head
NB.           specifies which triangular part of AA is to be
NB.           referenced:
NB.            'L' - lower, the form is:
NB.                    L * L^H = A
NB.            'U' - upper, the form is:
NB.                    U^H * U = A
NB.   AA   - n×n-matrix, contains either LT or UT or both
NB.          part(s) of A
NB.   B    - n×nrhs-matrix, RHS
NB.   T    - n×n-matrix, L or U factor
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, Hermitian (symmetric) positive
NB.          definite, to be factored to T
NB.   L    - n×n-matrix, lower triangular
NB.   U    - n×n-matrix, upper triangular
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale
NB. - no check for positive definiteness

dposv=: 4 : 0
  'AA B'=. y
  'n nrhs'=. $ B
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) AA
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  (|: L: 0) 4 6 { dposv_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: AA) ; ld ; (|: B) ; ld ; , _1
)

zposv=: 4 : 0
  'AA B'=. y
  'n nrhs'=. $ B
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) AA
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  (|: L: 0) 4 6 { zposv_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: AA) ; ld ; (|: B) ; ld ; , _1
)
