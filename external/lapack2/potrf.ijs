require 'math/lapack2'

NB. Description:
NB.   Compute the Cholesky factorization of a Hermitian
NB.   (symmetric) positive definite matrix
NB.
NB. Syntax:
NB.   T=. uplo xpotrf AA
NB. where
NB.   uplo - string, case-insensitive, in which the head
NB.           specifies which triangular part of AA is to be
NB.           referenced:
NB.            'L' - lower, the form is:
NB.                    L * L^H = A
NB.            'U' - upper, the form is:
NB.                    U^H * U = A
NB.   AA   - n×n-matrix, contains either LT or UT or both
NB.          part(s) of A
NB.   T    - n×n-matrix, L or U factor combined with
NB.          untouched elements from AA in opposite triangle
NB.   A    - n×n-matrix, Hermitian (symmetric) positive
NB.          definite, to be factored to T
NB.   L    - n×n-matrix, lower triangular
NB.   U    - n×n-matrix, upper triangular
NB.   n    ≥ 0, the size of A, AA, L and U
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale
NB. - no check for positive definiteness

dpotrf=: 4 : 0
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  |: 3 {:: dpotrf_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
)

zpotrf=: 4 : 0
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  |: 3 {:: zpotrf_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; , _1
)
