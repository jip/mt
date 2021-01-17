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
NB.   tau - k-vector, the scalar factors of elementary
NB.         reflectors applied to A
NB.   Qf  - m×k-matrix, the unit lower trapezoidal (unit
NB.         diagonal not stored), with the tau it represents
NB.         the Q in the factored form
NB.   R   - k×n-matrix, the upper trapezoidal
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
NB. - the verbs below are loaded into the current locale

dgeqp3=: 3 : 0
  'A pvt'=. y
  k=. <./ 'm n'=. $ A
  assert. (ismatrix_jlapack2_ *. isreal_jlapack2_         ) A
  assert. (isvector_jlapack2_ *. *./@(0 1 e.~ ,)  *. n = #) pvt
  select. 3!:0 A
    case. JCMPX do. A=. 9 o. A
    case. JFL   do.
    case.       do. A=. A + 0.0
  end.
  if. JINT ~: 3!:0 pvt do. pvt=. pvt + 00 end.
  NB. lwork=. , >: n  NB. minimal
  lwork=. , 1 >. (32 + 0 2) p. n  NB. optimal
  cdrc=. dgeqp3_jlapack2_ (, m) ; (, n) ; (|: A) ; (, 1 >. m) ; pvt ; (k $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 6 { cdrc  NB. (|:) doesn't affect to pvt and tau
)

zgeqp3=: 3 : 0
  'A pvt'=. y
  k=. <./ 'm n'=. $ A
  assert.  ismatrix_jlapack2_                               A
  assert. (isvector_jlapack2_ *. *./@(0 1 e.~ ,)  *. n = #) pvt
  if. JCMPX ~: 3!:0 A   do. A=.   A  + 0j0 end.
  if. JINT  ~: 3!:0 pvt do. pvt=. pvt + 00 end.
  NB. lwork=. , >: n  NB. minimal
  lwork=. , 1 >. 32 * >: n  NB. optimal
  cdrc=. zgeqp3_jlapack2_ (, m) ; (, n) ; (|: A) ; (, 1 >. m) ; pvt ; (k $ 0j0) ; (lwork $ 0j0) ; lwork ; ((+: n) $ 0.0) ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 3 5 6 { cdrc  NB. (|:) doesn't affect to pvt and tau
)
