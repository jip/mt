require 'math/lapack2'

NB. Description:
NB.   Compute the inverse of a factored Hermitian (symmetric)
NB.   positive definite matrix using the factorization
NB.   computed by xPOTRF
NB.
NB. Syntax:
NB.   iAA=. uplo xpotri AA
NB. where
NB.   uplo - literal, case-insensitive, in which the head
NB.           specifies which triangular part of AA is to be
NB.           referenced:
NB.            'L' - lower, the form is:
NB.                    L * L^H = A
NB.            'U' - upper, the form is:
NB.                    U^H * U = A
NB.   AA   - n×n-matrix, contains either lower or upper or
NB.          both part(s) of A
NB.   iAA  - n×n-matrix, contains the triangular part of iA
NB.          in changed triangle and unchanged elements of AA
NB.          in opposite strict triangle
NB.   A    - n×n-matrix, Hermitian (symmetric) positive
NB.          definite, to be factored to either L or U
NB.   iA   - n×n-matrix, Hermitian (symmetric) positive
NB.          definite, the inversion of A
NB.   n    ≥ 0, the size of A, AA, iA and iAA
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dpotri=: 4 : 0
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  |: 3 {:: dpotri_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
)

zpotri=: 4 : 0
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  |: 3 {:: zpotri_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
)
