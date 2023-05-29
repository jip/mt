require 'math/lapack2'

NB. Description:
NB.   Multiply a complex general matrix by a complex unitary
NB.   matrix which is defined as a product of elementary
NB.   reflectors
NB.
NB. Syntax:
NB.   B=. (side ; trans) zunmqr A ; tau ; C
NB. where
NB.   side  - literal, case-insensitive, in which the head
NB.           specifies the side of op(Q):
NB.             'L' - op(Q) * C  (apply op(Q) from the left)
NB.             'R' - C * op(Q)  (apply op(Q) from the right)
NB.   trans - literal, case-insensitive, in which the head
NB.           specifies the form of op(Q):
NB.             'N' - op(Q) := Q    (no transpose)
NB.             'C' - op(Q) := Q^H  (conjugate transpose)
NB.   A     - s×k-matrix, contains Qf
NB.   tau   - k-vector, scalar factors of elementary
NB.           reflectors as returned by ZGEQRF
NB.   C     - m×n-matrix, complex, the input to be multiplied
NB.           by op(Q)
NB.   B     - m×n-matrix, the result of multiplication:
NB.             Q   * C    if side='L' and trans='N'
NB.             Q^H * C    if side='L' and trans='C'
NB.             C   * Q    if side='R' and trans='N'
NB.             C   * Q^H  if side='R' and trans='C'
NB.   Qf    - s×k-matrix, unit lower trapezoidal (unit
NB.           diagonal is not stored), contains elementary
NB.           reflectors as returned by ZGEQRF, with the tau
NB.           it represents the Q in the factored form
NB.   Q     - s×s-matrix, complex, unitary, which is defined
NB.           as the product of k elementary reflectors H(i)
NB.           of order s:
NB.             Q = Π{H(i),i=0:k-1}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   m     ≥ 0, the number of rows in B and C
NB.   n     ≥ 0, the number of columns in B and C
NB.   s     = m if side='L' or s = n if side='R'
NB.   k     ∈ [0,s], the number of elementary reflectors
NB.           whose product defines the matrix Q
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

zunmqr=: 4 : 0
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
  |: 9 {:: zunmqr_jlapack2_ (, side) ; (, trans) ; (, m) ; (, n) ; (, k) ; (|: A) ; (, 1 >. s) ; tau ; (|: C) ; (, 1 >. m) ; (lwork $ 0j0) ; lwork ; , _1
)
