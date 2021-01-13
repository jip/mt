require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a general
NB.   square matrix, and factor the last one
NB.
NB. Syntax:
NB.   'L1U ipiv X'=. xgesv A ; B
NB. where
NB.   A    - n×n-matrix to be factored into L1U and ipiv as:
NB.            P * L1 * U = A
NB.   B    - n×nrhs-matrix, RHS
NB.   L1U  - n×n-matrix, L1 and U combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          P
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   L1   - n×n-matrix, the unit lower triangular (unit
NB.          diagonal not stored)
NB.   U    - n×n-matrix, the upper triangular
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dgesv=: 3 : 0
  'A B'=. y
  'n nrhs'=. $ B
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
  cdrc=. dgesv_jlapack2_ (, n) ; (, nrhs) ; (|: A) ; ld ; (n $ 00) ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 6 { cdrc  NB. (|:) doesn't affect to ipiv
)

zgesv=: 3 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. n = #) A
  assert.  ismatrix_jlapack2_                                 B
  if. JCMPX ~: 3!:0 A do. A=. A + 0j0 end.
  if. JCMPX ~: 3!:0 B do. B=. B + 0j0 end.
  ld=. , 1 >. n
  cdrc=. zgesv_jlapack2_ (, n) ; (, nrhs) ; (|: A) ; ld ; (n $ 00) ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 6 { cdrc  NB. (|:) doesn't affect to ipiv
)
