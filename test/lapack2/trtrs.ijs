require 'math/lapack2'

NB. Description:
NB.   Solve a triangular system:
NB.     op(A) * X = B
NB.
NB. Syntax:
NB.   X=. (uplo ; trans ; diag) xtrtrs A ; B
NB. where
NB.   uplo  - literal, case-insensitive, in which the head
NB.            specifies which triangular part of A is to be
NB.            referenced:
NB.             'L' - lower
NB.             'U' - upper
NB.   trans - literal, case-insensitive, in which the head
NB.           specifies the form of op(A):
NB.             'N' - op(A) := A    (no transpose)
NB.             'T' - op(A) := A^T  (transpose)
NB.             'C' - op(A) := A^H  (conjugate transpose)
NB.   diag  - literal, case-insensitive, in which the head
NB.            specifies the form of A:
NB.             'N' - A is non-unit triangular
NB.             'U' - A is unit triangular, the diagonal
NB.                   elements of A are not referenced
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
  assert. 'lLuU'   e.~ {. uplo
  assert. 'nNtTcC' e.~ {. trans
  assert. 'nNuU'   e.~ {. diag
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , isreal_jlapack2_ , n = #) A
  assert. (ismatrix_jlapack2_ ,                      isreal_jlapack2_        ) B
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
  assert. 'lLuU'   e.~ {. uplo
  assert. 'nNtTcC' e.~ {. trans
  assert. 'nNuU'   e.~ {. diag
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) A
  assert.  ismatrix_jlapack2_                               B
  if. JCMPX ~: 3!:0 A do. A=. A + 0j0 end.
  if. JCMPX ~: 3!:0 B do. B=. B + 0j0 end.
  ld=. , 1 >. n
  cdrc=. ztrtrs_jlapack2_ (, uplo) ; (, trans) ; (, diag) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  |: 8 {:: cdrc
)
