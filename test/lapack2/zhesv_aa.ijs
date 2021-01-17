require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   matrix, and factor the last one
NB.
NB. Syntax:
NB.   'DT1 ipiv X'=. uplo zhesv_aa A ; B
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of DT1 only, form is:
NB.                    P * L1 * D * L1^H * P^H = A
NB.            'U' - use upper triangle of DT1 only, form is:
NB.                    P * U1 * D * U1^H * P^H = A
NB.   A    - n×n-matrix, the Hermitian to be factored into
NB.          DT1 and ipiv
NB.   B    - n×nrhs-matrix, RHS
NB.   DT1  - n×n-matrix, D and T1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          P
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   D    - n×n-matrix, the Hermitian tridiagonal (opposite
NB.          diagonal not stored)
NB.   T1   - n×n-matrix, either L1 or U1
NB.   L1   - n×n-matrix, the unit lower triangular (unit
NB.          diagonal not stored)
NB.   U1   - n×n-matrix, the unit upper triangular (unit
NB.          diagonal not stored)
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

zhesv_aa=: 4 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. n = #) A
  assert.  ismatrix_jlapack2_                                 B
  if. JCMPX ~: 3!:0 A do. A=. A + 0j0 end.
  if. JCMPX ~: 3!:0 B do. B=. B + 0j0 end.
  ld=. , 1 >. n
  NB. lwork=. , >./ 1 , (+: n) , (_2 3 p. n)  NB. minimal
  lwork_zhetrf_aa=. 1 >. n * >: 64
  lwork_zhetrs_aa=. 1 >. _2 3 p. n
  lwork=. , lwork_zhetrf_aa >. lwork_zhetrs_aa  NB. optimal
  cdrc=. zhesv_aa_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (n $ 00) ; (|: B) ; ld ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 4 6 7 { cdrc  NB. (|:) doesn't affect to ipiv
)
