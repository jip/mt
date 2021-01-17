require 'math/lapack2'

NB. Description:
NB.   Compute the LU factorization of a general matrix using
NB.   partial pivoting with row interchanges
NB.
NB. Syntax:
NB.   'L1U ipiv'=. xgetrf A
NB. where
NB.   A    - m×n-matrix to be factored as:
NB.            P * L1 * U = A
NB.   L1U  - m×n-matrix, L1 and U combined
NB.   ipiv - k-vector, integer, the pivot indices that define
NB.          P
NB.   L1   - m×k-matrix, the unit lower trapezoidal (unit
NB.          diagonal not stored)
NB.   U    - k×n-matrix, the upper trapezoidal
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   m    ≥ 0, the number of rows in A and L1
NB.   n    ≥ 0, the number of columns in A and U
NB.   k    = min(m,n)
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dgetrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. (ismatrix_jlapack2_ *. isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  cdrc=. dgetrf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 00) ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to ipiv
)

zgetrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  cdrc=. zgetrf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 00) ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to ipiv
)
