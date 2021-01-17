require 'math/lapack2'

NB. Description:
NB.   Compute the QR factorization of a general matrix
NB.
NB. Syntax:
NB.   'QfR tau'=. xgeqrf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   QfR - m×n-matrix, Qf and R combined
NB.   tau - k-vector, the scalar factors of elementary
NB.         reflectors applied to A
NB.   Qf  - m×k-matrix, the unit lower trapezoidal (unit
NB.         diagonal not stored), with the tau it represents
NB.         the Q in the factored form
NB.   R   - k×n-matrix, the upper trapezoidal
NB.   Q   - m×m-matrix, the unitary (orthogonal), which is
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
NB. - the verbs below are loaded into the current locale

dgeqrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. (ismatrix_jlapack2_ *. isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) n  NB. optimal
  cdrc=. dgeqrf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to tau
)

zgeqrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 1 >. 32&*^:(k>128) n  NB. optimal
  cdrc=. zgeqrf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (k $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 { cdrc  NB. (|:) doesn't affect to tau
)
