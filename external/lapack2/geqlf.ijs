require 'math/lapack2'

NB. Description:
NB.   Compute the QL factorization of a general matrix
NB.
NB. Syntax:
NB.   'QfL tau'=. xgeqlf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   QfL - m×n-matrix, Qf and L combined
NB.   tau - k-vector, scalar factors of elementary
NB.         reflectors applied to A
NB.   Qf  - m×k-matrix, unit upper trapezoidal (unit diagonal
NB.         is not stored), with the tau it represents the Q
NB.         in the factored form
NB.   L   - k×n-matrix, lower trapezoidal
NB.   Q   - m×m-matrix, unitary (orthogonal), which is
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
NB. - verbs below are loaded into the current locale

dgeqlf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) n  NB. optimal
  (|:L:0) 3 5 { dgeqlf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)

zgeqlf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) n  NB. optimal
  (|:L:0) 3 5 { zgeqlf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)
