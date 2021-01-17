require 'math/lapack2'

NB. Description:
NB.   Compute the factorization of a Hermitian matrix using
NB.   the Bunch-Kaufman diagonal pivoting method
NB.
NB. Syntax:
NB.   'DPT1 ipiv'=. uplo zhetrf A
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of A only, form is:
NB.                    PL1 * D * PL1^H = A
NB.            'U' - use upper triangle of A only, form is:
NB.                    PU1 * D * PU1^H = A
NB.   A    - n×n-matrix, the Hermitian or upper or lower
NB.          triangular
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          permutations
NB.   D    - n×n-matrix, the Hermitian and block diagonal
NB.          with 1×1 and 2×2 diagonal blocks (opposite
NB.          diagonal not stored)
NB.   PT1  - n×n-matrix, either PL1 or PU1
NB.   PL1  - n×n-matrix, a product of permutation and unit
NB.          lower triangular matrices (unit diagonal not
NB.          stored)
NB.   PU1  - n×n-matrix, a product of permutation and unit
NB.          upper triangular matrices (unit diagonal not
NB.          stored)
NB.   n    ≥ 0, the size of A, D, PL1 and PU1
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

zhetrf=: 4 : 0
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_) y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  n=. # y
  NB. lwork=. , 1  NB. minimal
  lwork=. , 1 >. n * 64  NB. optimal
  cdrc=. zhetrf_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 00) ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to ipiv
)
