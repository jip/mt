require 'math/lapack2'

NB. Description:
NB.   Compute the LQ factorization of a general matrix
NB.
NB. Syntax:
NB.   'LQf tau'=. xgelqf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   LQf - m×n-matrix, L and Qf combined
NB.   tau - k-vector, the scalar factors of elementary
NB.         reflectors applied to A
NB.   L   - m×k-matrix, the lower trapezoidal
NB.   Qf  - k×n-matrix, the unit upper trapezoidal (unit
NB.         diagonal not stored), with the tau it represents
NB.         the Q in the factored form
NB.   Q   - n×n-matrix, the unitary (orthogonal), which is
NB.         defined as the product of k elementary reflectors
NB.         H(i) of order n:
NB.           Q = Π{H(i)',i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   k   = min(m,n)
NB.   m   ≥ 0, the number of rows in A, LQf and L
NB.   n   ≥ 0, the number of columns in A, LQf and Qf, and
NB.            the size of Q
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dgelqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. (ismatrix_jlapack2_ *. isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) m  NB. optimal
  cdrc=. dgelqf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to tau
)

zgelqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) m  NB. optimal
  cdrc=. zgelqf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to tau
)
