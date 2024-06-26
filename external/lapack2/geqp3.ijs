require 'math/lapack2'

NB. Description:
NB.   Compute the QR factorization with column pivoting of a
NB.   general matrix
NB.
NB. Syntax:
NB.   'QfR ip tau'=. xgeqp3 A ; pvt
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   pvt - n-vector, boolean, marks leading columns (which
NB.         are permuted to the front of A*P)
NB.   QfR - m×n-matrix, Qf and R combined
NB.   ip  - n-vector, integer, columns inversed permutation
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
NB. Assertions:
NB.   NB. if all columns are leading then the result of
NB.   NB. xGEQP3 matches with xGEQRF
NB.   (dgeqrf -: (< < < 1) { dgeqp3@(; 1 #~ c)) A
NB.   (zgeqrf -: (< < < 1) { zgeqp3@(; 1 #~ c)) A
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dgeqp3=: 3 : 0
  'A pvt'=. y
  k=. <./ 'm n'=. $ A
  assert.  ismatrix_jlapack2_          A
  assert. (isvector_jlapack2_ , n = #) pvt
  NB. lwork=. , >: n  NB. minimal
  lwork=. , 1 >. (32 + 0 2) p. n  NB. optimal
  (|:L:0) 3 5 6 { dgeqp3_jlapack2_ (, m) ; (, n) ; (|: A) ; (, 1 >. m) ; pvt ; (k $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to pvt and tau
)

zgeqp3=: 3 : 0
  'A pvt'=. y
  k=. <./ 'm n'=. $ A
  assert.  ismatrix_jlapack2_          A
  assert. (isvector_jlapack2_ , n = #) pvt
  NB. lwork=. , >: n  NB. minimal
  lwork=. , 1 >. 32 * >: n  NB. optimal
  (|:L:0) 3 5 6 { zgeqp3_jlapack2_ (, m) ; (, n) ; (|: A) ; (, 1 >. m) ; pvt ; (k $ 0j0) ; (lwork $ 0j0) ; lwork ; ((+: n) $ 0.0) ; , _1
    NB. (|:) doesn't affect to pvt and tau
)
