require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a Hermitian
NB.   (symmetric) positive definite tridiagonal matrix, and
NB.   factor the last one as:
NB.     L1 * D * L1^H = A
NB.
NB. Syntax:
NB.   'dD eT1 X'=. xptsv dA ; eA ; B
NB. where
NB.   dA  - n-vector, real, the diagonal of A
NB.   eA  - (n-1)-vector, the subdiagonal of A
NB.   B    - n×nrhs-matrix, RHS
NB.   dD   - n-vector, real, the diagonal of D
NB.   eT1  - (n-1)-vector, the subdiagonal of L1, or the
NB.          superdiagonal of U1
NB.   X    - n×nrhs-matrix, solutions of equation:
NB.            A * X = B
NB.   A    - n×n-matrix, Hermitian (symmetric) positive
NB.          definite tridiagonal represented by dA and eA
NB.   D    - n×n-matrix, diagonal, with dD as diagonal
NB.          elements
NB.   L1   - n×n-matrix, unit lower bidiagonal
NB.   U1   - n×n-matrix, unit upper bidiagonal
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale
NB. - no check for positive definiteness

dptsv=: 3 : 0
  'd e B'=. y
  'n nrhs'=. $ B
  assert. (isvector_jlapack2_ ,          n  = #) d
  assert. (isvector_jlapack2_ , (0 >. <: n) = #) e
  assert.  ismatrix_jlapack2_                    B
  (|: L: 0) 3 4 5 { dptsv_jlapack2_ (, n) ; (, nrhs) ; d ; e ; (|: B) ; (, 1 >. n) ; , _1
    NB. (|:) doesn't affect to d and e
)

zptsv=: 3 : 0
  'd e B'=. y
  'n nrhs'=. $ B
  assert. (isvector_jlapack2_ ,          n  = #) d
  assert. (isvector_jlapack2_ , (0 >. <: n) = #) e
  assert.  ismatrix_jlapack2_                    B
  (|: L: 0) 3 4 5 { zptsv_jlapack2_ (, n) ; (, nrhs) ; d ; e ; (|: B) ; (, 1 >. n) ; , _1
    NB. (|:) doesn't affect to d and e
)
