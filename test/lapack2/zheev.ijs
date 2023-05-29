require 'math/lapack2'

NB. Description:
NB.   Compute the eigenvalues and, optionally, eigenvectors
NB.   for square Hermitian matrix
NB.
NB. Syntax:
NB.   'w V'=. (jobV ; uplo) zheev A
NB. where
NB.   jobV - literal, case-insensitive, in which the head
NB.          specifies whether to compute V:
NB.            'N' - to not compute
NB.            'V' - to compute
NB.   uplo - literal, case-insensitive, in which the head
NB.          specifies which triangular part of A is to be
NB.          referenced:
NB.            'L' - lower
NB.            'U' - upper
NB.   A    - n×n-matrix, Hermitian or lower or upper
NB.   w    - n-vector, eigenvalues of A
NB.   V    - n×n-matrix, eigenvectors of A or 0×0-matrix
NB.   n    ≥ 0, the size of A, V and w
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

zheev=: 4 : 0
  'jobV uplo'=. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  NB. lwork=. , 1 >. _1 2 p. n  NB. minimal
  lwork=. , 1 >. n * 33  NB. optimal
  'w V'=. 6 4 { zheev_jlapack2_ (, jobV) ; (, uplo) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 0.0) ; (lwork $ 0j0) ; lwork ; ((1 >. _2 3 p. n) $ 0.0) ; , _1
  if. 'vV' e.~ {. jobV do. V=. |: V else. V=. EMPTY end.
  w ; V
)
