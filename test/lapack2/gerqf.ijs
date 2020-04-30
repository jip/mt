require 'math/lapack2'

NB. Description:
NB.   Compute the RQ factorization of a general matrix
NB.
NB. Syntax:
NB.   'RQf tau'=. xgerqf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   RQf - m×n-matrix, R and Qf combined
NB.   tau - k-vector, the scalar factors of the elementary
NB.         reflectors applied to A
NB.   R   - m×k-matrix, upper trapezoidal
NB.   Qf  - k×n-matrix, unit lower trapezoidal (unit diagonal
NB.         not stored), with the tau it represents the Q in
NB.         factored form
NB.   Q   - k×n-matrix with orthonormal rows, which is
NB.         defined as the last k rows of a product of m
NB.         elementary reflectors H(i) of order n:
NB.           Q = Π{H(i)',i=0:m-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.           H(k:m-1) ≡ H(u(m-1:k),τ(m-1:k)) = H(0,0) = I
NB.   k   = min(m,n)
NB.   m   ≥ 0
NB.   n   ≥ 0
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
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't effect tau
)

zgerqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) m  NB. optimal
  cdrc=. zgerqf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't effect tau
)
