require 'math/lapack2'

NB. Description:
NB.   Multiply a real general matrix by a real orthogonal
NB.   matrix which is defined as a product of elementary
NB.   reflectors
NB.
NB. Syntax:
NB.   B=. (side ; trans) dormhr ilo ; ihi ; A ; tau ; C
NB. where
NB.   side  - string, case-insensitive, in which the head
NB.           specifies the side of op(Q):
NB.             'L' - op(Q) * C  (apply op(Q) from the left)
NB.             'R' - C * op(Q)  (apply op(Q) from the right)
NB.   trans - string, case-insensitive, in which the head
NB.           specifies the form of op(Q):
NB.             'N' - op(Q) := Q    (no transpose)
NB.             'T' - op(Q) := Q^T  (transpose)
NB.   ilo   ∈ [1,max(1,ihi)], IO starting row and column,
NB.           1-based
NB.   ihi   ∈ [min(ilo,s),s], IO ending row and column,
NB.           1-based
NB.   A     - s×s-matrix, real, contains Qf
NB.   tau   - (s-1)-vector, real, scalar factors of
NB.           elementary reflectors as returned by DGEHRD
NB.   C     - m×n-matrix, real, the input to be multiplied by
NB.           op(Q)
NB.   B     - m×n-matrix, real, the result of multiplication:
NB.             Q   * C    if side='L' and trans='N'
NB.             Q^T * C    if side='L' and trans='T'
NB.             C   * Q    if side='R' and trans='N'
NB.             C   * Q^T  if side='R' and trans='T'
NB.   Qf    - s×s-matrix, columns below the first subdiagonal
NB.           with the tau represent the Q in the factored
NB.           form as returned by DGEHRD
NB.   Q     - s×s-matrix, real, orthogonal, which is defined
NB.           as the product of (ihi-ilo) elementary
NB.           reflectors H(i) of order s:
NB.             Q = Π{H(i),i=ilo:ihi-1}
NB.             H(i) = I - v[i] * τ[i] * v[i]'
NB.           it is equal to the unit matrix except in the
NB.           submatrix Q(ilo+1:ihi,ilo+1:ihi)
NB.   m     ≥ 0, the number of rows in B and C
NB.   n     ≥ 0, the number of columns in B and C
NB.   s     = m if side='L' or s = n if side='R'
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dormhr=: 4 : 0
  'side trans'=. x
  'ilo ihi A tau C'=. y
  'm n'=. sh=. $ C
  s=. # A
  assert. s = sh {~ 'rR' e.~ {. side
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) A
  assert. (isvector_jlapack2_ , (0 >. <: s) = #   ) tau
  assert.  ismatrix_jlapack2_                       C
  NB. lwork=. , 1 >. sh {~ 'lL' e.~ {. side  NB. minimal
  lwork=. , 32 * 1 >. sh {~ 'lL' e.~ {. side  NB. optimal
  |: 10 {:: dormhr_jlapack2_ (, side) ; (, trans) ; (, m) ; (, n) ; (, ilo) ; (, ihi) ; (|: A) ; (, 1 >. s) ; tau ; (|: C) ; (, 1 >. m) ; (lwork $ 0.0) ; lwork ; , _1
)
