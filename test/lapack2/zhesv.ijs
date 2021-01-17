require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   matrix, and factor the last one
NB.
NB. Syntax:
NB.   'DPT1 ipiv X'=. uplo zhesv A ; B
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of A only, form is:
NB.                    PL1 * D * PL1^H = A
NB.            'U' - use upper triangle of A only, form is:
NB.                    PU1 * D * PU1^H = A
NB.   A    - n×n-matrix, the Hermitian to be factored into
NB.          DPT1 and ipiv
NB.   B    - n×nrhs-matrix, RHS
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          permutations
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
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
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

zhesv=: 4 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. n = #) A
  assert.  ismatrix_jlapack2_                                 B
  if. JCMPX ~: 3!:0 A do. A=. A + 0j0 end.
  if. JCMPX ~: 3!:0 B do. B=. B + 0j0 end.
  ld=. , 1 >. n
  nb=. 1:^:((x e. 'uU') *. 64 > n) 64
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. n * nb  NB. optimal
  cdrc=. zhesv_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (n $ 00) ; (|: B) ; ld ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 4 6 7 { cdrc  NB. (|:) doesn't affect to ipiv
)
