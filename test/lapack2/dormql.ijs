require 'math/lapack2'

NB. Description:
NB.   Multiply a real general matrix by a real orthogonal
NB.   matrix which is defined as a product of elementary
NB.   reflectors
NB.
NB. Syntax:
NB.   B=. (side ; trans) dormql A ; tau ; C
NB. where
NB.   side  - literal, case-insensitive, in which the head
NB.           specifies the side of op(Q):
NB.             'L' - op(Q) * C  (apply op(Q) from the left)
NB.             'R' - C * op(Q)  (apply op(Q) from the right)
NB.   trans - literal, case-insensitive, in which the head
NB.           specifies the form of op(Q):
NB.             'N' - op(Q) := Q    (no transpose)
NB.             'T' - op(Q) := Q^T  (transpose)
NB.   A     - s×k-matrix, contains Qf
NB.   tau   - k-vector, the scalar factors of elementary
NB.           reflectors as returned by DGEQLF
NB.   C     - m×n-matrix, real, the input to be multiplied by
NB.           op(Q)
NB.   B     - m×n-matrix, the result of multiplication:
NB.             Q   * C    if side='L' and trans='N'
NB.             Q^T * C    if side='L' and trans='T'
NB.             C   * Q    if side='R' and trans='N'
NB.             C   * Q^T  if side='R' and trans='T'
NB.   Qf    - s×k-matrix, the unit upper trapezoidal (unit
NB.           diagonal not stored), contains elementary
NB.           reflectors as returned by DGEQLF, with the tau
NB.           it represents the Q in the factored form
NB.   Q     - s×s-matrix, real, orthogonal, which is defined
NB.           as the product of k elementary reflectors H(i)
NB.           of order s:
NB.             Q = Π{H(i),i=k-1:0}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   m     ≥ 0, the number of rows in B and C
NB.   n     ≥ 0, the number of columns in B and C
NB.   s     = m if side='L' or s = n if side='R'
NB.   k     ∈ [0,s], the number of elementary reflectors
NB.           whose product defines the matrix Q
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dormql=: 4 : 0
  'side trans'=. x
  'A tau C'=. y
  'm n'=. sh=. $ C
  's k'=. $ A
  assert.  ismatrix_jlapack2_          A
  assert. (isvector_jlapack2_ , k = #) tau
  assert.  ismatrix_jlapack2_          C
  NB. lwork=. , 1 >. sh {~ 'lL' e.~ {. side  NB. minimal
  nbmax=. 64
  ilaenv=. 32
  ldt=. >: nbmax
  tsize=. nbmax * ldt
  nb=. nbmax <. ilaenv
  lwork=. , tsize + nb * 1 >. sh {~ 'lL' e.~ {. side  NB. optimal
  cdrc=. dormql_jlapack2_ (, side) ; (, trans) ; (, m) ; (, n) ; (, k) ; (|: A) ; (, 1 >. s) ; tau ; (|: C) ; (, 1 >. m) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 9 {:: cdrc
)
