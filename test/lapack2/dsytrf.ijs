require 'math/lapack2'

NB. Description:
NB.   Compute the factorization of a real symmetric matrix
NB.   using the Bunch-Kaufman diagonal pivoting method
NB.
NB. Syntax:
NB.   'DPT1 ipiv'=. uplo dsytrf AA
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
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, pivot indices that define
NB.          permutations
NB.   A    - n×n-matrix, real, symmetric, to factorize
NB.   D    - n×n-matrix, real, symmetric and block diagonal
NB.          with 1×1 and 2×2 diagonal blocks (opposite
NB.          diagonal is not stored)
NB.   PT1  - n×n-matrix, either PL1 or PU1
NB.   PL1  - n×n-matrix, real, the product of permutation and
NB.          unit lower triangular matrices
NB.   PU1  - n×n-matrix, real, the product of permutation and
NB.          unit upper triangular matrices
NB.   n    ≥ 0, the size of A, AA, D, PL1 and PU1
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dsytrf=: 4 : 0
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  NB. lwork=. , 1  NB. minimal
  lwork=. , 1 >. n * 64  NB. optimal
  (|: L: 0) 3 5 { dsytrf_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 00) ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to ipiv
)
