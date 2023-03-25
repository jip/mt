require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   matrix using the factorization computed by ZHETRF_AA
NB.
NB. Syntax:
NB.   X=. uplo zhetrs_aa DT1 ; ipiv ; B
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.          specifies which triangular part of DT1 is to be
NB.          referenced:
NB.            'L' - lower, the form is:
NB.                    P * L1 * D * L1^H * P^H = A
NB.            'U' - upper, the form is:
NB.                    P * U1 * D * U1^H * P^H = A
NB.   DT1  - n×n-matrix, D and T1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          P
NB.   B    - n×nrhs-matrix, RHS
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, the Hermitian represented in
NB.          factored form by DT1 and ipiv
NB.   D    - n×n-matrix, the Hermitian tridiagonal (opposite
NB.          diagonal not stored)
NB.   T1   - n×n-matrix, either L1 or U1
NB.   L1   - n×n-matrix, the unit lower triangular (unit
NB.          diagonal not stored)
NB.   U1   - n×n-matrix, the unit upper triangular (unit
NB.          diagonal not stored)
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

zhetrs_aa=: 4 : 0
  'DT1 ipiv B'=. y
  'n nrhs'=. $ B
  assert. 'lLuU' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) DT1
  assert. (isvector_jlapack2_ ,                      n = #) ipiv
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  lwork=. , 1 >. _2 3 p. n  NB. minimal
  cdrc=. zhetrs_aa_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: DT1) ; ld ; ipiv ; (|: B) ; ld ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 7 {:: cdrc
)
