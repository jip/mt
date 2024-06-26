require 'math/lapack2'

NB. Description:
NB.   Compute the RQ factorization of a general matrix
NB.
NB. Syntax:
NB.   'RQf tau'=. xgerqf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   RQf - m×n-matrix, R and Qf combined
NB.   tau - k-vector, scalar factors of elementary reflectors
NB.         applied to A
NB.   R   - m×k-matrix, upper trapezoidal
NB.   Qf  - k×n-matrix, unit lower trapezoidal (unit diagonal
NB.         is not stored), with the tau it represents the Q
NB.         in the factored form
NB.   Q   - n×n-matrix, unitary (orthogonal), which is
NB.         defined as the product of m elementary reflectors
NB.         H(i) of order n:
NB.           Q = Π{H(i)',i=0:m-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   k   = min(m,n)
NB.   m   ≥ 0, the number of rows in A, RQf and R
NB.   n   ≥ 0, the number of columns in A, RQf and Qf and Q
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dgerqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) m  NB. optimal
  (|:L:0) 3 5 { dgerqf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)

zgerqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) m  NB. optimal
  (|:L:0) 3 5 { zgerqf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)
