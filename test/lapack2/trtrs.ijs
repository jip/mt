require 'math/lapack2'

NB. Description:
NB.   Solve a triangular system
NB.
NB. Syntax:
NB.   X=. (uplo , trans , diag) xtrtrs A ; B
NB. where
NB.   uplo  - scalar, character, case-insensitive:
NB.             'L' - use lower triangle of A only
NB.             'U' - use upper triangle of A only
NB.   trans - scalar, character, case-insensitive, specifies
NB.           the form of the system of equations:
NB.             'N' - A   * X = B  (no transpose)
NB.             'T' - A^T * X = B  (transpose)
NB.             'C' - A^H * X = B  (conjugate transpose)
NB.   diag  - scalar, character, case-insensitive:
NB.             'N' - A is non-unit triangular
NB.             'U' - A is unit triangular, the diagonal
NB.                   elements of A are not
NB.   A     - n×n-matrix, [unit] {lower,upper}-triangular
NB.   B     - n×nrhs-matrix, RHS
NB.   X     - n×nrhs-matrix, solutions
NB.   n     ≥ 0, the order of system
NB.   nrhs  ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale
NB. - trans='T' and trans='C' are identic for dtrtrs

dtrtrs=: 4 : 0
  'uplo trans diag'=. x
  'A B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # uplo ) , uplo  e. 'lLuU'
  assert. *./ (1 = # trans) , trans e. 'nNtTcC'
  assert. *./ (1 = # diag ) , diag  e. 'nNuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_ *. n = #) A
  assert. (ismatrix_jlapack2_ *.                       isreal_jlapack2_         ) B
  select. 3!:0 A
    case. JCMPX do. A=. 9 o. A
    case. JFL   do.
    case.       do. A=. A + 0.0
  end.
  select. 3!:0 B
    case. JCMPX do. B=. 9 o. B
    case. JFL   do.
    case.       do. B=. B + 0.0
  end.
  ld=. , 1 >. n
  cdrc=. dtrtrs_jlapack2_ (, uplo) ; (, trans) ; (, diag) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  |: 8 {:: cdrc
)

ztrtrs=: 4 : 0
  'uplo trans diag'=. x
  'A B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # uplo ) , uplo  e. 'lLuU'
  assert. *./ (1 = # trans) , trans e. 'nNtTcC'
  assert. *./ (1 = # diag ) , diag  e. 'nNuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. n = #) A
  assert.  ismatrix_jlapack2_                                 B
  if. JCMPX ~: 3!:0 A do. A=. A + 0j0 end.
  if. JCMPX ~: 3!:0 B do. B=. B + 0j0 end.
  ld=. , 1 >. n
  cdrc=. ztrtrs_jlapack2_ (, uplo) ; (, trans) ; (, diag) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  |: 8 {:: cdrc
)
