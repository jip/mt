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
NB.   trans - string, case-insensitive, in which the head
NB.           specifies the form of the system of equations:
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
NB. - verbs below are loaded into the current locale

dgels=: 4 : 0
  'A B'=. y
  'm n'=. $ A
  'k nrhs'=. $ B
  assert. ismatrix_jlapack2_ A
  assert. ismatrix_jlapack2_ B
  if. 'nN' e.~ {. x do.
    assert. m = k
    xh=. n
  else.
    assert. n = k
    xh=. m
  end.
  lda=. , 1 >. m
  NB. lwork=. , 1 >. nrhs (] + >.) m <. n  NB. minimal
  lwork=. , 1 >. nrhs (] + 32 * >.) m <. n  NB. optimal
  xh {. |: 7 {:: dgels_jlapack2_ (, x) ; (, m) ; (, n) ; (, nrhs) ; (|: A) ; lda ; (|: B) ; (, lda >. n) ; (lwork $ 0.0) ; lwork ; , _1
)

zgels=: 4 : 0
  'A B'=. y
  'm n'=. $ A
  'k nrhs'=. $ B
  assert. ismatrix_jlapack2_ A
  assert. ismatrix_jlapack2_ B
  if. 'nN' e.~ {. x do.
    assert. m = k
    xh=. n
  else.
    assert. n = k
    xh=. m
  end.
  lda=. , 1 >. m
  NB. lwork=. , 1 >. nrhs (] + >.) m <. n  NB. minimal
  lwork=. , 1 >. nrhs (] + 32 * >.) m <. n  NB. optimal
  xh {. |: 7 {:: zgels_jlapack2_ (, x) ; (, m) ; (, n) ; (, nrhs) ; (|: A) ; lda ; (|: B) ; (, lda >. n) ; (lwork $ 0j0) ; lwork ; , _1
)
