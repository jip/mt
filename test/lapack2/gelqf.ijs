require 'math/lapack2'

NB. Description:
NB.   Compute the LQ factorization of a general matrix
NB.
NB. Syntax:
NB.   'LQf tau'=. xgelqf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   LQf - m×n-matrix, L and Qf combined
NB.   tau - k-vector, scalar factors of elementary reflectors
NB.         applied to A
NB.   L   - m×k-matrix, lower trapezoidal
NB.   Qf  - k×n-matrix, unit upper trapezoidal (unit diagonal
NB.         is not stored), with the tau it represents the Q
NB.         in the factored form
NB.   Q   - n×n-matrix, unitary (orthogonal), which is
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
NB. - verbs below are loaded into the current locale

dgelqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) m  NB. optimal
  (|: L: 0) 3 5 { dgelqf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)

zgelqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) m  NB. optimal
  (|: L: 0) 3 5 { zgelqf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)
