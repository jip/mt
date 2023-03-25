require 'math/lapack2'

NB. Description:
NB.   Multiply a real general matrix by a real orthogonal
NB.   matrix which is defined as a product of elementary
NB.   reflectors
NB.
NB. Syntax:
NB.   B=. (side ; trans) dormlq A ; tau ; C
NB. where
NB.   side  - literal, case-insensitive, in which the head
NB.           specifies the side of op(Q):
NB.             'L' - op(Q) * C  (apply op(Q) from the left)
NB.             'R' - C * op(Q)  (apply op(Q) from the right)
NB.   trans - literal, case-insensitive, in which the head
NB.           specifies the form of op(Q):
NB.             'N' - op(Q) := Q    (no transpose)
NB.             'T' - op(Q) := Q^T  (transpose)
NB.   A     - k×s-matrix, contains Qf
NB.   tau   - k-vector, the scalar factors of elementary
NB.           reflectors as returned by DGELQF
NB.   C     - m×n-matrix, real, the input to be multiplied by
NB.           op(Q)
NB.   B     - m×n-matrix, the result of multiplication:
NB.             Q   * C    if side='L' and trans='N'
NB.             Q^T * C    if side='L' and trans='T'
NB.             C   * Q    if side='R' and trans='N'
NB.             C   * Q^T  if side='R' and trans='T'
NB.   Qf    - k×s-matrix, the unit upper trapezoidal (unit
NB.           diagonal not stored), contains elementary
NB.           reflectors as returned by DGELQF, with the tau
NB.           it represents the Q in the factored form
NB.   Q     - s×s-matrix, real, orthogonal, which is defined
NB.           as the product of k elementary reflectors H(i)
NB.           of order s:
NB.             Q = Π{H(i)',i=k-1:0}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   m     ≥ 0, the number of rows in B and C
NB.   n     ≥ 0, the number of columns in B and C
NB.   s     = m if side='L' or s = n if side='R'
NB.   k     ∈ [0,s], the number of elementary reflectors
NB.           whose product defines the matrix Q
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dormlq=: 4 : 0
  'side trans'=. x
  'A tau C'=. y
  'm n'=. sh=. $ C
  'k s'=. $ A
  assert. 'lLrR' e.~ {. side
  assert. 'nNtT' e.~ {. trans
  assert. s = sh {~ 'rR' e.~ {. side
  assert. (_1 , s) I. k
  assert. (ismatrix_jlapack2_ , isreal_jlapack2_        ) A
  assert. (isvector_jlapack2_ , isreal_jlapack2_ , k = #) tau
  assert. (ismatrix_jlapack2_ , isreal_jlapack2_        ) C
  NB. lwork=. , 1 >. sh {~ 'lL' e.~ {. side  NB. minimal
  nbmax=. 64
  ilaenv=. 32
  ldt=. >: nbmax
  tsize=. nbmax * ldt
  nb=. nbmax <. ilaenv
  lwork=. , tsize + nb * 1 >. sh {~ 'lL' e.~ {. side  NB. optimal
  cdrc=. dormlq_jlapack2_ (, side) ; (, trans) ; (, m) ; (, n) ; (, k) ; (|: A) ; (, 1 >. k) ; tau ; (|: C) ; (, 1 >. m) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 9 {:: cdrc
)
