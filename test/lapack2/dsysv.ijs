require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a non-complex
NB.   symmetric matrix, and factor the last one
NB.
NB. Syntax:
NB.   'DPT1 ipiv X'=. uplo dsysv AA ; B
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.          specifies which triangular part of AA is to be
NB.          referenced:
NB.            'L' - lower, the form is:
NB.                    PL1 * D * PL1^T = A
NB.            'U' - upper, the form is:
NB.                    PU1 * D * PU1^T = A
NB.   AA   - n×n-matrix, real, contains either lower or upper
NB.          or both part(s) of A
NB.   A    - n×n-matrix, real, symmetric, to be factored to
NB.          DPT1 and ipiv
NB.   B    - n×nrhs-matrix, real, RHS
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, pivot indices that define
NB.          permutations
NB.   X    - n×nrhs-matrix, real, solutions of equation:
NB.            A * X = B
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

dsysv=: 4 : 0
  'AA B'=. y
  'n nrhs'=. $ B
  assert. 'lLuU' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) AA
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. n * 64  NB. optimal
  (|: L: 0) 4 6 7 { dsysv_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: AA) ; ld ; (n $ 00) ; (|: B) ; ld ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to ipiv
)
