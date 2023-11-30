require 'math/lapack2'

NB. Description:
NB.   Compute the eigenvalues and, optionally, eigenvectors
NB.   for real square symmetric matrix
NB.
NB. Syntax:
NB.   'w V'=. (jobV ; uplo) dsyev AA
NB. where
NB.   jobV - literal, case-insensitive, in which the head
NB.          specifies whether to compute V:
NB.            'N' - to not compute
NB.            'V' - to compute
NB.   uplo - literal, case-insensitive, in which the head
NB.          specifies which triangular part of AA is to be
NB.          referenced:
NB.            'L' - lower
NB.            'U' - upper
NB.   AA   - n×n-matrix, real, contains either LT or UT or
NB.          both part(s) of A
NB.   w    - n-vector, real, eigenvalues of A
NB.   V    - n×n-matrix, real, eigenvectors of A or
NB.          0×0-matrix
NB.   A    - n×n-matrix, real, symmetric, to be decomposed
NB.   n    ≥ 0, the size of A, AA, V and w
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
