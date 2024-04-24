require 'math/lapack2'

NB. Description:
NB.   Compute the factorization of a Hermitian matrix using
NB.   the Aasen's algorithm
NB.
NB. Syntax:
NB.   'DT1 ipiv'=. uplo zhetrf_aa AA
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
NB.   DT1  - n×n-matrix, D and T1 combined
NB.   ipiv - n-vector, integer, pivot indices that define P
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
NB.   n    ≥ 0, the size of A, AA, D, L1, U1 and P
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

zhetrf_aa=: 4 : 0
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  NB. lwork=. , 1 >. +: n  NB. minimal
  lwork=. , 1 >. n * >: 64  NB. optimal
  (|:L:0) 3 5 { zhetrf_aa_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 00) ; (lwork $ 0j0) ; lwork ; , _1
    NB. (|:) doesn't affect to ipiv
)
