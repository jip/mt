require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   matrix using the factorization computed by ZHETRF_AA
NB.
NB. Syntax:
NB.   X=. uplo zhetrs_aa DT1 ; ipiv ; B
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of DT1 only, form is:
NB.                    P * L1 * D * L1^H * P^H = A
NB.            'U' - use upper triangle of DT1 only, form is:
NB.                    P * U1 * D * U1^H * P^H = A
NB.   DT1  - n×n-matrix, D and T1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          P
NB.   B    - n×nrhs-matrix, RHS
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, the Hermitian represented in
NB.          factored form by DT1 and ipiv
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

zhetrs_aa=: 4 : 0
  'DT1 ipiv B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. n = #) DT1
  assert. (isvector_jlapack2_ *. (-: <.) :: 0:      *. n = #) ipiv
  assert.  ismatrix_jlapack2_                                 B
  if. JCMPX ~: 3!:0 DT1  do. DT1=.  DT1 + 0j0    end.
  if. JCMPX ~: 3!:0 B    do. B=.    B   + 0j0    end.
  if. JINT  ~: 3!:0 ipiv do. ipiv=. <. 9 o. ipiv end.
  ld=. , 1 >. n
  lwork=. , 1 >. _2 3 p. n  NB. minimal
  cdrc=. zhetrs_aa_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: DT1) ; ld ; ipiv ; (|: B) ; ld ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 7 {:: cdrc
)
