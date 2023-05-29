require 'math/lapack2'

NB. Description:
NB.   Reduce an upper trapezoidal matrix to upper triangular
NB.   form by an unitary (orthogonal) transformations
NB.
NB. Syntax:
NB.   'RZf tau'=. xtzrzf A
NB. where
NB.   A   - m×n-matrix, the leading upper trapezoidal part of
NB.         it must contain the matrix to be factored
NB.   RZf - m×n-matrix, R and Zf combined
NB.   tau - m-vector, scalar factors of elementary
NB.         reflectors applied to A
NB.   R   - m×m-matrix, the upper triangular part contains
NB.         the part factored, the strict lower triangular
NB.         part contains corresp. elements from A unchanged
NB.   Zf  - m×n-matrix, unit upper trapezoidal (unit diagonal
NB.         is not stored), with the tau it represents the Z
NB.         in the factored form
NB.   Z   - n×n-matrix, unitary (orthogonal), which is
NB.         defined as the product of m elementary reflectors
NB.         H(i) of order n:
NB.           Z = Π{H(i)',i=0:m-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   m   ≥ 0, the number of rows in A, RZf and Zf, and the
NB.         size of R
NB.   n   ≥ m, the number of columns in A, RZf and Zf, and
NB.         the size of Z
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dtzrzf=: 3 : 0
  'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32&*`1:@.(m=n) m  NB. optimal
  (|: L: 0) 3 5 { dtzrzf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (m $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)

ztzrzf=: 3 : 0
  'm n'=. $ y
  assert. ismatrix_jlapack2_ y
  NB. lwork=. , 1 >. m  NB. minimal
  lwork=. , 1 >. 32&*`1:@.(m=n) m  NB. optimal
  (|: L: 0) 3 5 { ztzrzf_jlapack2_ (, m) ; (, n) ; (|: y) ; (, 1 >. m) ; (m $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)
