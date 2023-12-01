require 'math/lapack2'

NB. Description:
NB.   Compute the inverse of a factored Hermitian matrix
NB.
NB. Syntax:
NB.   iAA=. uplo zhetri2 DPT1 ; ipiv
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.          specifies which triangular part of DPT1 is to be
NB.          referenced:
NB.            'L' - lower, the form is:
NB.                    PL1 * D * PL1^H = A
NB.            'U' - upper, the form is:
NB.                    PU1 * D * PU1^H = A
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, pivot indices that define
NB.          permutations
NB.   D    - n×n-matrix, Hermitian and block diagonal with
NB.          1×1 and 2×2 diagonal blocks (opposite diagonal
NB.          is not stored)
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
NB.   iA   - n×n-matrix, the Hermitian inversion of A
NB.   A    - n×n-matrix, Hermitian, to invert, represented in
NB.          factored form by DPT1 and ipiv
NB.   n    ≥ 0, the size of A, iA and iAA
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

zhetri2=: 4 : 0
  'DPT1 ipiv'=. y
  n=. # ipiv
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) DPT1
  assert.  isvector_jlapack2_                               ipiv
  lwork=. , (n + 64 + 1) * (64 + 3)  NB. minimal
  |: 3 {:: zhetri2_jlapack2_ (, x) ; (, n) ; (|: DPT1) ; (, 1 >. n) ; ipiv ; (lwork $ 0j0) ; lwork ; , _1
)
