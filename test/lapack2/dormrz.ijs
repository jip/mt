require 'math/lapack2'

NB. Description:
NB.   Multiply a real general matrix by a real orthogonal
NB.   matrix which is defined as a product of elementary
NB.   reflectors
NB.
NB. Syntax:
NB.   B=. (side ; trans) dormrz l ; A ; tau ; C
NB. where
NB.   side  - literal, case-insensitive, in which the head
NB.           specifies the side of op(Z):
NB.             'L' - op(Z) * C  (apply op(Z) from the left)
NB.             'R' - C * op(Z)  (apply op(Z) from the right)
NB.   trans - literal, case-insensitive, in which the head
NB.           specifies the form of op(Z):
NB.             'N' - op(Z) := Z    (no transpose)
NB.             'T' - op(Z) := Z^T  (transpose)
NB.   l     ∈ [0,s], the number of columns of the matrix A
NB.           containing the meaningful part of the
NB.           Householder reflectors
NB.   A     - k×s-matrix, contains Zf
NB.   tau   - k-vector, the scalar factors of elementary
NB.           reflectors as returned by DTZRZF
NB.   C     - m×n-matrix, real, the input to be multiplied by
NB.           op(Z)
NB.   B     - m×n-matrix, the result of multiplication:
NB.             Z   * C    if side='L' and trans='N'
NB.             Z^T * C    if side='L' and trans='T'
NB.             C   * Z    if side='R' and trans='N'
NB.             C   * Z^T  if side='R' and trans='T'
NB.   Zf    - k×s-matrix, the unit lower trapezoidal (unit
NB.           diagonal not stored), contains elementary
NB.           reflectors as returned by DTZRZF, with the tau
NB.           it represents the Z in the factored form
NB.   Z     - s×s-matrix, real, orthogonal, which is defined
NB.           as the product of k elementary reflectors H(i)
NB.           of order s:
NB.             Z = Π{H(i)',i=0:k-1}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   m     ≥ 0, the number of rows in B and C
NB.   n     ≥ 0, the number of columns in B and C
NB.   s     = m if side='L' or s = n if side='R'
NB.   k     ∈ [0,s], the number of elementary reflectors
NB.           whose product defines the matrix Z
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dormrz=: 4 : 0
  'side trans'=. x
  'l A tau C'=. y
  'm n'=. sh=. $ C
  'k s'=. $ A
  assert. 'lLrR' e.~ {. side
  assert. 'nNtT' e.~ {. trans
  assert. s = sh {~ 'rR' e.~ {. side
  assert. (_1 , s) I. k , l
  assert. (ismatrix_jlapack2_ , isreal_jlapack2_        ) A
  assert. (isvector_jlapack2_ , isreal_jlapack2_ , k = #) tau
  assert. (ismatrix_jlapack2_ , isreal_jlapack2_        ) C
  if. JFL ~: 3!:0 A   do. A=.   9 o. A   end.
  if. JFL ~: 3!:0 tau do. tau=. 9 o. tau end.
  if. JFL ~: 3!:0 C   do. C=.   9 o. C   end.
  NB. lwork=. , 1 >. sh {~ 'lL' e.~ {. side  NB. minimal
  nbmax=. 64
  ilaenv=. 32
  ldt=. >: nbmax
  tsize=. nbmax * ldt
  nb=. nbmax <. ilaenv
  lwork=. , tsize + nb * 1 >. sh {~ 'lL' e.~ {. side  NB. optimal
  cdrc=. dormrz_jlapack2_ (, side) ; (, trans) ; (, m) ; (, n) ; (, k) ; (, l) ; (|: A) ; (, 1 >. k) ; tau ; (|: C) ; (, 1 >. m) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 10 {:: cdrc
)
