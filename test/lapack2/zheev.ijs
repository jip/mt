require 'math/lapack2'

NB. Description:
NB.   Compute the eigenvalues and, optionally, eigenvectors
NB.   for square Hermitian matrix
NB.
NB. Syntax:
NB.   'w V'=. (jobV , uplo) zheev A
NB. where
NB.   jobV - scalar, character, case-insensitive:
NB.            'N' - do not compute V
NB.            'V' - to compute V
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of A only
NB.            'U' - use upper triangle of A only
NB.   A    - n×n-matrix, the Hermitian or lower or upper
NB.          triangular
NB.   w    - n-vector, eigenvalues of A
NB.   V    - n×n-matrix, eigenvectors of A or 0×0-matrix
NB.   n    ≥ 0, the size of A, V and w
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

zheev=: 4 : 0
  'jobV uplo'=. x
  assert. *./ (1 = # jobV) , jobV e. 'nNvV'
  assert. *./ (1 = # uplo) , uplo e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_@((<0 1)&|:)) y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  n=. # y
  NB. lwork=. , 1 >. _1 2 p. n  NB. minimal
  lwork=. , 1 >. n * 33  NB. optimal
  cdrc=. zheev_jlapack2_ (, jobV) ; (, uplo) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 0.0) ; (lwork $ 0j0) ; lwork ; ((1 >. _2 3 p. n) $ 0.0) ; , _1
  assert. 0 = _1 {:: cdrc
  'w V'=. 6 4 { cdrc
  if. jobV e. 'vV' do. V=. |: V else. V=. EMPTY end.
  w ; V
)
