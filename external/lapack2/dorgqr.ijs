require 'math/lapack2'

NB. Description:
NB.   Generate a real matrix with orthonormal columns from
NB.   output of DGEQRF
NB.
NB. Syntax:
NB.   Q=. dorgqr A ; tau
NB. where
NB.   A   - m×n-matrix, real, contains Qf
NB.   tau - k-vector, real, scalar factors of elementary
NB.         reflectors as returned by DGEQRF
NB.   Q   - m×n-matrix, real, with orthonormal columns, which
NB.         is defined as the first n columns of the product
NB.         of k elementary reflectors H(i) of order m:
NB.           Q = Π{H(i),i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   Qf  - m×n-matrix, unit lower trapezoidal (unit diagonal
NB.         is not stored), contains elementary reflectors as
NB.         returned by DGEQRF, with the tau it represents
NB.         the Q in the factored form
NB.   m   ≥ n, the number of rows in A, Q and Qf
NB.   n   ≥ 0, the number of columns in A, Q and Qf
NB.   k   ∈ [0,n], the number of elementary reflectors whose
NB.         product defines the matrix Q
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dorgqr=: 3 : 0
  'A tau'=. y
  'm n'=. $ A
  k=. # tau
  assert. ismatrix_jlapack2_ A
  assert. isvector_jlapack2_ tau
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32 *^:(128 < ]) n  NB. optimal
  |: 4 {:: dorgqr_jlapack2_ (, m) ; (, n) ; (, k) ; (|: A) ; (, 1 >. m) ; tau ; (lwork $ 0.0) ; lwork ; , _1
)
