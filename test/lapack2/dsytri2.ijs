require 'math/lapack2'

NB. Description:
NB.   Compute the inverse of a factored symmetric matrix
NB.
NB. Syntax:
NB.   iAA=. uplo dsytri2 DPT1 ; ipiv
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.          specifies which triangular part of A is to be
NB.          referenced:
NB.            'L' - lower, the form is:
NB.                    PL1 * D * PL1^T = A
NB.            'U' - upper, the form is:
NB.                    PU1 * D * PU1^T = A
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          permutations
NB.   D    - n×n-matrix, the symmetric and block diagonal
NB.          with 1×1 and 2×2 diagonal blocks (opposite
NB.          diagonal not stored)
NB.   PT1  - n×n-matrix, either PL1 or PU1
NB.   PL1  - n×n-matrix, the product of permutation and unit
NB.          lower triangular matrices (unit diagonal not
NB.          stored)
NB.   PU1  - n×n-matrix, the product of permutation and unit
NB.          upper triangular matrices (unit diagonal not
NB.          stored)
NB.   iAA  - n×n-matrix, contains the triangular part of iA
NB.          in changed triangle and unchanged elements of
NB.          DPT1 in opposite strict triangle
NB.   iA   - n×n-matrix, the symmetric inversion of A
NB.   A    - n×n-matrix, the symmetric represented in factored
NB.          form by DPT1 and ipiv
NB.   n    ≥ 0, the size of A and iA
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dsytri2=: 4 : 0
  'DPT1 ipiv'=. y
  assert. 'lLuU' e.~ {. x
  n=. # ipiv
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) DPT1
  assert.  isvector_jlapack2_                               ipiv
  lwork=. , (n + 64 + 1) * (64 + 3)  NB. minimal
  cdrc=. dsytri2_jlapack2_ (, x) ; (, n) ; (|: DPT1) ; (, 1 >. n) ; ipiv ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 3 {:: cdrc
)
