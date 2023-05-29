require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a non-complex
NB.   symmetric matrix, and factor the last one
NB.
NB. Syntax:
NB.   'DT1 ipiv X'=. uplo dsysv_aa A ; B
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.          specifies which triangular part of DT1 is to be
NB.          referenced:
NB.            'L' - lower, the form is:
NB.                    P * L1 * D * L1^T * P^T = A
NB.            'U' - upper, the form is:
NB.                    P * U1 * D * U1^T * P^T = A
NB.   A    - n×n-matrix, real, symmetric to be factored to
NB.          DT1 and ipiv
NB.   B    - n×nrhs-matrix, real, RHS
NB.   DT1  - n×n-matrix, D and T1 combined
NB.   ipiv - n-vector, integer, pivot indices that define P
NB.   X    - n×nrhs-matrix, real, solutions of equation:
NB.            A * X = B
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

dsysv_aa=: 4 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. 'lLuU' e.~ {. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) A
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  NB. lwork=. , >./ 1 , (+: n) , (_2 3 p. n)  NB. minimal
  lwork_dsytrf_aa=. 1 >. n * >: 64
  lwork_dsytrs_aa=. 1 >. _2 3 p. n
  lwork=. , lwork_dsytrf_aa >. lwork_dsytrs_aa  NB. optimal
  (|: L: 0) 4 6 7 { dsysv_aa_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (n $ 00) ; (|: B) ; ld ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to ipiv
)
