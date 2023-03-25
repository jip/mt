require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a non-complex
NB.   symmetric matrix, and factor the last one
NB.
NB. Syntax:
NB.   'DPT1 ipiv X'=. uplo dsysv A ; B
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.          specifies which triangular part of A is to be
NB.          referenced:
NB.            'L' - lower, the form is:
NB.                    PL1 * D * PL1^T = A
NB.            'U' - upper, the form is:
NB.                    PU1 * D * PU1^T = A
NB.   A    - n×n-matrix, real, the symmetric to be factored
NB.          to DPT1 and ipiv
NB.   B    - n×nrhs-matrix, real, RHS
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          permutations
NB.   X    - n×nrhs-matrix, real, solutions of equation:
NB.            A * X = B
NB.   D    - n×n-matrix, real, the symmetric and block
NB.          diagonal with 1×1 and 2×2 diagonal blocks
NB.          (opposite diagonal not stored)
NB.   PT1  - n×n-matrix, either PL1 or PU1
NB.   PL1  - n×n-matrix, real, the product of permutation and
NB.          unit lower triangular matrices
NB.   PU1  - n×n-matrix, real, the product of permutation and
NB.          unit upper triangular matrices
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dsysv=: 4 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. 'lLuU' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , isreal_jlapack2_ , n = #) A
  assert. (ismatrix_jlapack2_ ,                      isreal_jlapack2_        ) B
  if. JFL ~: 3!:0 A do. A=. 9 o. A end.
  if. JFL ~: 3!:0 B do. B=. 9 o. B end.
  ld=. , 1 >. n
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. n * 64  NB. optimal
  cdrc=. dsysv_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (n $ 00) ; (|: B) ; ld ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 4 6 7 { cdrc  NB. (|:) doesn't affect to ipiv
)
