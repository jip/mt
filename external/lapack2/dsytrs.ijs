require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a real
NB.   symmetric matrix using the factorization computed by
NB.   DSYTRF
NB.
NB. Syntax:
NB.   X=. uplo dsytrs DPT1 ; ipiv ; B
NB. where
NB.   uplo - string, case-insensitive, in which the head
NB.          specifies which triangular part of DPT1 is to be
NB.          referenced:
NB.            'L' - lower, the form is:
NB.                    PL1 * D * PL1^T = A
NB.            'U' - upper, the form is:
NB.                    PU1 * D * PU1^T = A
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, pivot indices that define
NB.          permutations
NB.   B    - n×nrhs-matrix, real, RHS
NB.   X    - n×nrhs-matrix, real, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, real, symmetric, represented in the
NB.          factored form by DPT1 and ipiv
NB.   D    - n×n-matrix, real, symmetric and block diagonal
NB.          with 1×1 and 2×2 diagonal blocks (opposite
NB.          diagonal is not stored)
NB.   PT1  - n×n-matrix, either PL1 or PU1
NB.   PL1  - n×n-matrix, real, the product of permutation and
NB.          unit lower triangular matrices
NB.   PU1  - n×n-matrix, real, the product of permutation and
NB.          unit upper triangular matrices
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dsytrs=: 4 : 0
  'DPT1 ipiv B'=. y
  'n nrhs'=. $ B
  assert. 'lLuU' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) DPT1
  assert. (isvector_jlapack2_ ,                      n = #) ipiv
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  |: 7 {:: dsytrs_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: DPT1) ; ld ; ipiv ; (|: B) ; ld ; , _1
)
