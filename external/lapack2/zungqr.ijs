require 'math/lapack2'

NB. Description:
NB.   Generate a complex matrix with orthonormal columns from
NB.   output of ZGEQRF
NB.
NB. Syntax:
NB.   Q=. zungqr A ; tau
NB. where
NB.   A   - m×n-matrix, contains Qf
NB.   tau - k-vector, scalar factors of elementary reflectors
NB.         as returned by ZGEQRF
NB.   Q   - m×n-matrix with orthonormal columns, which is
NB.         defined as the first n columns of the product of
NB.         k elementary reflectors H(i) of order m:
NB.           Q = Π{H(i),i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   Qf  - m×n-matrix, unit lower trapezoidal (unit diagonal
NB.         is not stored), contains elementary reflectors as
NB.         returned by ZGEQRF, with the tau it represents
NB.         the Q in the factored form
NB.   m   ≥ n, the number of rows in A, Q and Qf
NB.   n   ≥ 0, the number of columns in A, Q and Qf
NB.   k   ∈ [0,n], the number of elementary reflectors whose
NB.         product defines the matrix Q
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

zungqr=: 3 : 0
  'A tau'=. y
  'm n'=. $ A
  k=. # tau
  assert. ismatrix_jlapack2_ A
  assert. isvector_jlapack2_ tau
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32 *^:(128 < ]) n  NB. optimal
  |: 4 {:: zungqr_jlapack2_ (, m) ; (, n) ; (, k) ; (|: A) ; (, 1 >. m) ; tau ; (lwork $ 0j0) ; lwork ; , _1
)
