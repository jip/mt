require 'math/lapack2'

NB. Description:
NB.   Solves overdetermined or underdetermined linear system
NB.   involving a matrix of full rank, or its
NB.   [conjugate-]transpose
NB.
NB. Syntax:
NB.   X1=. 'N' xgels A1 ; B1  NB. case 1: A1   * X1 = B1
NB.   X2=. 'N' xgels A2 ; B2  NB. case 2: A2   * X2 = B2
NB.   X2=. 'T' dgels A1 ; B2  NB. case 3: A1^T * X2 = B2
NB.   X2=. 'C' zgels A1 ; B2  NB. case 3: A1^H * X2 = B2
NB.   X1=. 'T' dgels A2 ; B1  NB. case 4: A2^T * X1 = B1
NB.   X1=. 'C' zgels A2 ; B1  NB. case 4: A2^H * X1 = B1
NB. where
NB.   trans - scalar, character, case-insensitive, specifies
NB.           the form of the system of equations:
NB.             'N' - A   * X = B  (no transpose)
NB.             'T' - A^T * X = B  (transpose)
NB.             'C' - A^H * X = B  (conjugate transpose)
NB.   A1    - m1×n1-matrix of full rank, m1>=n1, will be
NB.           factored by QR
NB.   B1    - m1×nrhs-matrix, RHS
NB.   X1    - n1×nrhs-matrix, solutions
NB.   A2    - m2×n2-matrix of full rank, m2<n2, will be
NB.           factored by LQ
NB.   B2    - m2×nrhs-matrix, RHS
NB.   X2    - n2×nrhs-matrix, solutions
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dgels=: 4 : 0
  'A B'=. y
  'm n'=. $ A
  'k nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'nNtT'
  assert. (ismatrix_jlapack2_ *. isreal_jlapack2_) A
  assert. (ismatrix_jlapack2_ *. isreal_jlapack2_) B
  if. x e. 'nN' do.
    assert. m = k
    xh=. n
  else.
    assert. n = k
    xh=. m
  end.
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
  lda=. , 1 >. m
  ldb=. , lda >. n
  B=. ldb {. B
  NB. lwork=. , 1 >. nrhs (] + >.) m <. n  NB. minimal
  lwork=. , 1 >. nrhs (] + 32 * >.) m <. n  NB. optimal
  cdrc=. dgels_jlapack2_ (, x) ; (, m) ; (, n) ; (, nrhs) ; (|: A) ; lda ; (|: B) ; ldb ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  xh {. |: 7 {:: cdrc
)

zgels=: 4 : 0
  'A B'=. y
  'm n'=. $ A
  'k nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'nNcC'
  assert. ismatrix_jlapack2_ A
  assert. ismatrix_jlapack2_ B
  if. x e. 'nN' do.
    assert. m = k
    xh=. n
  else.
    assert. n = k
    xh=. m
  end.
  if. JCMPX ~: 3!:0 A do. A=. A + 0j0 end.
  if. JCMPX ~: 3!:0 B do. B=. B + 0j0 end.
  lda=. , 1 >. m
  ldb=. , lda >. n
  B=. ldb {. B
  NB. lwork=. , 1 >. nrhs (] + >.) m <. n  NB. minimal
  lwork=. , 1 >. nrhs (] + 32 * >.) m <. n  NB. optimal
  cdrc=. zgels_jlapack2_ (, x) ; (, m) ; (, n) ; (, nrhs) ; (|: A) ; lda ; (|: B) ; ldb ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  xh {. |: 7 {:: cdrc
)
