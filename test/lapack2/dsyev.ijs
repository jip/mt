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
NB.   A    - n×n-matrix, real, the symmetric or
NB.          lower or upper triangular
NB.   w    - n-vector, real, eigenvalues of A
NB.   V    - n×n-matrix, eigenvectors of A or 0×0-matrix
NB.   n    ≥ 0, the size of A, V and w
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dsyev=: 4 : 0
  'jobV uplo'=. x
  assert. 'nNvV' e.~ {. jobV
  assert. 'lLuU' e.~ {. uplo
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  n=. # y
  NB. lwork=. , 1 >. _1 3 p. n  NB. minimal
  lwork=. , 1 >. n * 34  NB. optimal
  cdrc=. dsyev_jlapack2_ (, jobV) ; (, uplo) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  'w V'=. 6 4 { cdrc
  if. 'vV' e.~ {. jobV do. V=. |: V else. V=. EMPTY end.
  w ; V
)
