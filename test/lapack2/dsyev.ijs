require 'math/lapack2'

NB. Description:
NB.   Compute the eigenvalues and, optionally, eigenvectors
NB.   for non-complex square symmetric matrix
NB.
NB. Syntax:
NB.   'w V'=. (jobV ; uplo) dsyev A
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
NB.   A    - n×n-matrix, real, symmetric or lower or upper
NB.          triangular
NB.   w    - n-vector, real, eigenvalues of A
NB.   V    - n×n-matrix, eigenvectors of A or 0×0-matrix
NB.   n    ≥ 0, the size of A, V and w
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dsyev=: 4 : 0
  'jobV uplo'=. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  NB. lwork=. , 1 >. _1 3 p. n  NB. minimal
  lwork=. , 1 >. n * 34  NB. optimal
  'w V'=. 6 4 { dsyev_jlapack2_ (, jobV) ; (, uplo) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
  if. 'vV' e.~ {. jobV do. V=. |: V else. V=. EMPTY end.
  w ; V
)
