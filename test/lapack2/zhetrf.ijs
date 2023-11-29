require 'math/lapack2'

NB. Description:
NB.   Compute the factorization of a Hermitian matrix using
NB.   the Bunch-Kaufman diagonal pivoting method
NB.
NB. Syntax:
NB.   'DPT1 ipiv'=. uplo zhetrf AA
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.          specifies which triangular part of AA is to be
NB.          referenced:
NB.            'L' - lower, the form is:
NB.                    PL1 * D * PL1^H = A
NB.            'U' - upper, the form is:
NB.                    PU1 * D * PU1^H = A
NB.   AA   - n×n-matrix, contains either lower or upper or
NB.          both part(s) of A
NB.   A    - n×n-matrix, Hermitian, to be factored to DPT1
NB.          and ipiv
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, pivot indices that define
NB.          permutations
NB.   D    - n×n-matrix, Hermitian and block diagonal with
NB.          1×1 and 2×2 diagonal blocks (opposite diagonal
NB.          is not stored)
NB.   PT1  - n×n-matrix, either PL1 or PU1
NB.   PL1  - n×n-matrix, a product of permutation and unit
NB.          lower triangular matrices (unit diagonal not
NB.          stored)
NB.   PU1  - n×n-matrix, a product of permutation and unit
NB.          upper triangular matrices (unit diagonal not
NB.          stored)
NB.   n    ≥ 0, the size of A, AA, D, PL1 and PU1
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

zhetrf=: 4 : 0
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  NB. lwork=. , 1  NB. minimal
  lwork=. , 1 >. n * 64  NB. optimal
  (|: L: 0) 3 5 { zhetrf_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 00) ; (lwork $ 0j0) ; lwork ; , _1
    NB. (|:) doesn't affect to ipiv
)
