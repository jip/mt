require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   (symmetric) positive definite matrix, and factor the
NB.   last one
NB.
NB. Syntax:
NB.   'T X'=. uplo xposv A ; B
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of A only, form is:
NB.                    L * L^H = A
NB.            'U' - use upper triangle of A only, form is:
NB.                    U^H * U = A
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite to be factored to T
NB.   B    - n×nrhs-matrix, RHS
NB.   T    - n×n-matrix, L or U factor
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   L    - n×n-matrix, the lower triangular
NB.   U    - n×n-matrix, the upper triangular
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale
NB. - no check for positive definiteness

dposv=: 4 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'lLuU'
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
  cdrc=. dposv_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 4 6 { cdrc
)

zposv=: 4 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. n = #) A
  assert.  ismatrix_jlapack2_                                 B
  if. JCMPX ~: 3!:0 A do. A=. A + 0j0 end.
  if. JCMPX ~: 3!:0 B do. B=. B + 0j0 end.
  ld=. , 1 >. n
  cdrc=. zposv_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (|: B) ; ld ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 4 6 { cdrc
)
