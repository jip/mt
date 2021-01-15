require 'math/lapack2'

NB. Description:
NB.   Compute the RQ factorization of a general matrix
NB.
NB. Syntax:
NB.   'RQf tau'=. xgerqf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   RQf - m×n-matrix, R and Qf combined
NB.   tau - k-vector, the scalar factors of elementary
NB.         reflectors applied to A
NB.   R   - m×k-matrix, the upper trapezoidal
NB.   Qf  - k×n-matrix, the unit lower trapezoidal (unit
NB.         diagonal not stored), with the tau it represents
NB.         the Q in the factored form
NB.   Q   - n×n-matrix, the unitary (orthogonal), which is
NB.         defined as the product of m elementary reflectors
NB.         H(i) of order n:
NB.           Q = Π{H(i)',i=0:m-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   k   = min(m,n)
NB.   m   ≥ 0, the number of rows in A, RQf and R
NB.   n   ≥ 0, the number of columns in A, RQf and Qf and Q
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dgerqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. (ismatrix_jlapack2_ *. isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) m  NB. optimal
  cdrc=. dgerqf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to tau
)

zgerqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) m  NB. optimal
  cdrc=. zgerqf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to tau
)
