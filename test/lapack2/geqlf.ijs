require 'math/lapack2'

NB. Description:
NB.   Compute the QL factorization of a general matrix
NB.
NB. Syntax:
NB.   'QfL tau'=. xgeqlf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   QfL - m×n-matrix, Qf and L combined
NB.   tau - k-vector, the scalar factors of elementary
NB.         reflectors applied to A
NB.   Qf  - m×k-matrix, the unit upper trapezoidal (unit
NB.         diagonal not stored), with the tau it represents
NB.         the Q in the factored form
NB.   L   - k×n-matrix, the lower trapezoidal
NB.   Q   - m×m-matrix, the unitary (orthogonal), which is
NB.         defined as the product of k elementary reflectors
NB.         H(i) of order m:
NB.           Q = Π{H(i),i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   k   = min(m,n)
NB.   m   ≥ 0, the number of rows in A, QfL and Qf
NB.   n   ≥ 0, the number of columns in A, QfL and L, and the
NB.         size of Q
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dgeqlf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. (ismatrix_jlapack2_ *. isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) n  NB. optimal
  cdrc=. dgeqlf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to tau
)

zgeqlf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) n  NB. optimal
  cdrc=. zgeqlf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to tau
)
