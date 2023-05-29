require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a non-complex
NB.   symmetric matrix using the factorization computed by
NB.   DSYTRF_AA
NB.
NB. Syntax:
NB.   X=. uplo dsytrs_aa DT1 ; ipiv ; B
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.          specifies which triangular part of DT1 is to be
NB.          referenced:
NB.            'L' - lower, the form is:
NB.                    P * L1 * D * L1^T * P^T = A
NB.            'U' - upper, the form is:
NB.                    P * U1 * D * U1^T * P^T = A
NB.   DT1  - n×n-matrix, D and T1 combined
NB.   ipiv - n-vector, integer, pivot indices that define P
NB.   B    - n×nrhs-matrix, real, RHS
NB.   X    - n×nrhs-matrix, real, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, real, symmetric, represented in the
NB.          factored form by DT1 and ipiv
NB.   D    - n×n-matrix, real, symmetric tridiagonal
NB.          (opposite diagonal is not stored)
NB.   T1   - n×n-matrix, either L1 or U1
NB.   L1   - n×n-matrix, real, unit lower triangular (unit
NB.          diagonal is not stored)
NB.   U1   - n×n-matrix, real, unit upper triangular (unit
NB.          diagonal is not stored)
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dsytrs_aa=: 4 : 0
  'DT1 ipiv B'=. y
  'n nrhs'=. $ B
  assert. 'lLuU' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) DT1
  assert. (isvector_jlapack2_ ,                      n = #) ipiv
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  lwork=. , 1 >. _2 3 p. n  NB. minimal
  |: 7 {:: dsytrs_aa_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: DT1) ; ld ; ipiv ; (|: B) ; ld ; (lwork $ 0.0) ; lwork ; , _1
)
