require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   matrix, and factor the last one
NB.
NB. Syntax:
NB.   'DT1 ipiv X'=. uplo zhesv_aa AA ; B
NB. where
NB.   uplo - string, case-insensitive, in which the head
NB.          specifies which triangular part of AA is to be
NB.          referenced:
NB.            'L' - lower, the form is:
NB.                    P * L1 * D * L1^H * P^H = A
NB.            'U' - upper, the form is:
NB.                    P * U1 * D * U1^H * P^H = A
NB.   AA   - n×n-matrix, contains either LT or UT or both
NB.          part(s) of A
NB.   B    - n×nrhs-matrix, RHS
NB.   DT1  - n×n-matrix, D and T1 combined
NB.   ipiv - n-vector, integer, pivot indices that define P
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, Hermitian, to be factored to DT1 and
NB.          ipiv
NB.   D    - n×n-matrix, Hermitian tridiagonal (opposite
NB.          diagonal is not stored)
NB.   T1   - n×n-matrix, either L1 or U1
NB.   L1   - n×n-matrix, unit lower triangular (unit diagonal
NB.          is not stored)
NB.   U1   - n×n-matrix, unit upper triangular (unit diagonal
NB.          is not stored)
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

zhesv_aa=: 4 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) A
  assert.  ismatrix_jlapack2_                               B
  ld=. , 1 >. n
  NB. lwork=. , >./ 1 , (+: n) , (_2 3 p. n)  NB. minimal
  lwork_zhetrf_aa=. 1 >. n * >: 64
  lwork_zhetrs_aa=. 1 >. _2 3 p. n
  lwork=. , lwork_zhetrf_aa >. lwork_zhetrs_aa  NB. optimal
  (|:L:0) 4 6 7 { zhesv_aa_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (n $ 00) ; (|: B) ; ld ; (lwork $ 0j0) ; lwork ; , _1
    NB. (|:) doesn't affect to ipiv
)
