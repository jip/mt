require 'math/lapack2'

NB. Description:
NB.   Compute the QR factorization of a general matrix
NB.
NB. Syntax:
NB.   'QfR tau'=. xgeqrf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   QfR - m×n-matrix, Qf and R combined
NB.   tau - k-vector, scalar factors of elementary reflectors
NB.         applied to A
NB.   Qf  - m×k-matrix, unit lower trapezoidal (unit diagonal
NB.         is not stored), with the tau it represents the Q
NB.         in the factored form
NB.   R   - k×n-matrix, upper trapezoidal
NB.   Q   - m×m-matrix, unitary (orthogonal), which is
NB.         defined as the product of k elementary reflectors
NB.         H(i) of order m:
NB.           Q = Π{H(i),i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   k   = min(m,n)
NB.   m   ≥ 0, the number of rows in A, QfR and Qf
NB.   n   ≥ 0, the number of columns in A, QfR and R, and the
NB.         size of Q
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dgeqrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) n  NB. optimal
  (|: L: 0) 3 5 { dgeqrf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)

zgeqrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) n  NB. optimal
  (|: L: 0) 3 5 { zgeqrf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)
