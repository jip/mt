require 'math/lapack2'

NB. Description:
NB.   Multiply a complex general matrix by a complex unitary
NB.   matrix which is defined as a product of elementary
NB.   reflectors
NB.
NB. Syntax:
NB.   B=. (side , trans) zunmlq A ; tau ; C
NB. where
NB.   side  - scalar, character, case-insensitive, specifies
NB.           the side of op(Q):
NB.             'L' - op(Q) * C  (apply op(Q) from the left)
NB.             'R' - C * op(Q)  (apply op(Q) from the right)
NB.   trans - scalar, character, case-insensitive, specifies
NB.           op(Q):
NB.             'N' - op(Q) := Q    (no transpose)
NB.             'C' - op(Q) := Q^H  (conjugate transpose)
NB.   A     - lda×s-matrix, contains Qf
NB.   tau   - k-vector, the scalar factors of elementary
NB.           reflectors as returned by ZGELQF
NB.   C     - m×n-matrix, complex, the input to be multiplied
NB.           by op(Q)
NB.   B     - m×n-matrix, the result of multiplication:
NB.             Q   * C    if side='L' and trans='N'
NB.             Q^H * C    if side='L' and trans='C'
NB.             C   * Q    if side='R' and trans='N'
NB.             C   * Q^H  if side='R' and trans='C'
NB.   Qf    - k×s-matrix, the unit upper trapezoidal (unit
NB.           diagonal not stored), contains elementary
NB.           reflectors as returned by ZGELQF, with the tau
NB.           it represents the Q in the factored form
NB.   Q     - s×s-matrix, complex, unitary, which is defined
NB.           as the product of k elementary reflectors H(i)
NB.           of order s:
NB.             Q = Π{H(i)',i=k-1:0}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   m     ≥ 0, the number of rows in B and C
NB.   n     ≥ 0, the number of columns in B and C
NB.   s     = m if side='L' or s = n if side='R'
NB.   k     ∈ [0,s], the number of elementary reflectors
NB.           whose product defines the matrix Q
NB.   lda   ≥ max(1,k)
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

zunmlq=: 4 : 0
  'side trans'=. x
  'A tau C'=. y
  'm n'=. sh=. $ C
  'lda s'=. $ A
  k=. # tau
  assert. *./ (1 = # side ) , side  e. 'lLrR'
  assert. *./ (1 = # trans) , trans e. 'nNcC'
  assert. s = sh {~ side e. 'rR'
  assert. (_1 , s) I. k
  assert. lda >: 1 >. k
  assert. ismatrix_jlapack2_ A
  assert. isvector_jlapack2_ tau
  assert. ismatrix_jlapack2_ C
  if. JCMPX ~: 3!:0 A   do. A=.   A   + 0j0 end.
  if. JCMPX ~: 3!:0 tau do. tau=. tau + 0j0 end.
  if. JCMPX ~: 3!:0 C   do. C=.   C   + 0j0 end.
  NB. lwork=. , 1 >. sh {~ side e. 'lL'  NB. minimal
  nbmax=. 64
  ilaenv=. 32
  ldt=. >: nbmax
  tsize=. nbmax * ldt
  nb=. nbmax <. ilaenv
  lwork=. , tsize + nb * 1 >. sh {~ side e. 'lL'  NB. optimal
  cdrc=. zunmlq_jlapack2_ (, side) ; (, trans) ; (, m) ; (, n) ; (, k) ; (|: A) ; (, lda) ; tau ; (|: C) ; (, 1 >. m) ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 9 {:: cdrc
)
