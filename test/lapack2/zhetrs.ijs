require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   matrix using the factorization computed by ZHETRF
NB.
NB. Syntax:
NB.   X=. uplo zhetrs DPT1 ; ipiv ; B
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of DPT1 only, form
NB.                  is:
NB.                    PL1 * D * PL1^H = A
NB.            'U' - use upper triangle of DPT1 only, form
NB.                  is:
NB.                    PU1 * D * PU1^H = A
NB.   DPT1 - n×n-matrix, D and PT1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          permutations
NB.   B    - n×nrhs-matrix, RHS
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, the Hermitian represented in the
NB.          factored form by DPT1 and ipiv
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

zhetrs=: 4 : 0
  'DPT1 ipiv B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. n = #) DPT1
  assert. (isvector_jlapack2_ *. (-: <.) :: 0:      *. n = #) ipiv
  assert.  ismatrix_jlapack2_                                 B
  if. JCMPX ~: 3!:0 DPT1 do. DPT1=. DPT1 + 0j0   end.
  if. JCMPX ~: 3!:0 B    do. B=.    B    + 0j0   end.
  if. JINT  ~: 3!:0 ipiv do. ipiv=. <. 9 o. ipiv end.
  ld=. , 1 >. n
  cdrc=. zhetrs_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: DPT1) ; ld ; ipiv ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  |: 7 {:: cdrc
)
