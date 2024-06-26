require 'math/lapack2'

NB. Description:
NB.   Generate a complex matrix with orthonormal rows from
NB.   output of ZGELQF
NB.
NB. Syntax:
NB.   Q=. zunglq A ; tau
NB. where
NB.   A   - m×n-matrix, contains Qf
NB.   tau - k-vector, scalar factors of elementary reflectors
NB.         as returned by ZGELQF
NB.   Q   - m×n-matrix with orthonormal rows, which is
NB.         defined as first m rows of the product of k
NB.         elementary reflectors H(i) of order n:
NB.           Q = Π{H(i)',i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   Qf  - m×n-matrix, unit upper trapezoidal (unit diagonal
NB.         is not stored), contains elementary reflectors as
NB.         returned by ZGELQF, with the tau it represents
NB.         the Q in the factored form
NB.   m   ≥ 0, the number of rows in A, Q and Qf
NB.   n   ≥ m, the number of columns in A, Q and Qf
NB.   k   ∈ [0,m], the number of elementary reflectors whose
NB.         product defines the matrix Q
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

zunglq=: 3 : 0
  'A tau'=. y
  'm n'=. $ A
  k=. # tau
  assert. ismatrix_jlapack2_ A
  assert. isvector_jlapack2_ tau
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32 *^:(128 < ]) m  NB. optimal
  |: 4 {:: zunglq_jlapack2_ (, m) ; (, n) ; (, k) ; (|: A) ; (, 1 >. m) ; tau ; (lwork $ 0j0) ; lwork ; , _1
)
