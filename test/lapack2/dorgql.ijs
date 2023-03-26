require 'math/lapack2'

NB. Description:
NB.   Generate a real matrix with orthonormal columns from
NB.   output of DGEQLF
NB.
NB. Syntax:
NB.   Q=. dorgql A ; tau
NB. where
NB.   A   - m×n-matrix, contains Qf
NB.   tau - k-vector, the scalar factors of elementary
NB.         reflectors as returned by DGEQLF
NB.   Q   - m×n-matrix with orthonormal columns, which is
NB.         defined as the last n columns of the product of k
NB.         elementary reflectors H(i) of order m:
NB.           Q = Π{H(i),i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   Qf  - m×n-matrix, the unit upper trapezoidal (unit
NB.         diagonal not stored), contains elementary
NB.         reflectors as returned by DGEQLF, with the tau it
NB.         represents the Q in the factored form
NB.   m   ≥ n, the number of rows in A, Q and Qf
NB.   n   ≥ 0, the number of columns in A, Q and Qf
NB.   k   ∈ [0,n], the number of elementary reflectors whose
NB.         product defines the matrix Q
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dorgql=: 3 : 0
  'A tau'=. y
  'm n'=. $ A
  k=. # tau
  assert. ismatrix_jlapack2_ A
  assert. isvector_jlapack2_ tau
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32 *^:(128 < ]) n  NB. optimal
  |: 4 {:: dorgql_jlapack2_ (, m) ; (, n) ; (, k) ; (|: A) ; (, 1 >. m) ; tau ; (lwork $ 0.0) ; lwork ; , _1
)
