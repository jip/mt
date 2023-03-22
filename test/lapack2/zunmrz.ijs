require 'math/lapack2'

NB. Description:
NB.   Multiply a complex general matrix by a complex unitary
NB.   matrix which is defined as a product of elementary
NB.   reflectors
NB.
NB. Syntax:
NB.   B=. (side , trans) zunmrz l ; A ; tau ; C
NB. where
NB.   side  - scalar, character, case-insensitive, specifies
NB.           the side of op(Z):
NB.             'L' - op(Z) * C  (apply op(Z) from the left)
NB.             'R' - C * op(Z)  (apply op(Z) from the right)
NB.   trans - scalar, character, case-insensitive, specifies
NB.           op(Z):
NB.             'N' - op(Z) := Z    (no transpose)
NB.             'C' - op(Z) := Z^H  (transpose)
NB.   l     ∈ [0,s], the number of columns of the matrix A
NB.           containing the meaningful part of the
NB.           Householder reflectors
NB.   A     - lda×s-matrix, contains Zf
NB.   tau   - k-vector, the scalar factors of elementary
NB.           reflectors as returned by ZTZRZF
NB.   C     - m×n-matrix, complex, the input to be multiplied
NB.           by op(Z)
NB.   B     - m×n-matrix, the result of multiplication:
NB.             Z   * C    if side='L' and trans='N'
NB.             Z^H * C    if side='L' and trans='C'
NB.             C   * Z    if side='R' and trans='N'
NB.             C   * Z^H  if side='R' and trans='C'
NB.   Zf    - k×s-matrix, the unit lower trapezoidal (unit
NB.           diagonal not stored), contains elementary
NB.           reflectors as returned by ZTZRZF, with the tau
NB.           it represents the Z in the factored form
NB.   Z     - s×s-matrix, complex, unitary, which is defined
NB.           as the product of k elementary reflectors H(i)
NB.           of order s:
NB.             Z = Π{H(i)',i=0:k-1}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   m     ≥ 0, the number of rows in B and C
NB.   n     ≥ 0, the number of columns in B and C
NB.   s     = m if side='L' or s = n if side='R'
NB.   k     ∈ [0,s], the number of elementary reflectors
NB.           whose product defines the matrix Z
NB.   lda   ≥ max(1,k)
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

zunmrz=: 4 : 0
  'side trans'=. x
  'l A tau C'=. y
  'm n'=. sh=. $ C
  'lda s'=. $ A
  k=. # tau
  assert. (e.&'lLrR' , #) side
  assert. (e.&'nNcC' , #) trans
  assert. s = sh {~ side e. 'rR'
  assert. (_1 , s) I. k , l
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
  cdrc=. zunmrz_jlapack2_ (, side) ; (, trans) ; (, m) ; (, n) ; (, k) ; (, l) ; (|: A) ; (, lda) ; tau ; (|: C) ; (, 1 >. m) ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 10 {:: cdrc
)
