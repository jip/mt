NB. rcond.ijs
NB. Estimate the reciprocal of the condition number of a
NB. matrix, in either the 1-norm or the ∞-norm
NB.
NB. gecon  adverb to estimate the reciprocal of the
NB.        condition number of a general matrix
NB. hecon  adverb to estimate the reciprocal of the
NB.        condition number of a Hermitian (symmetric)
NB.        matrix
NB. dicon  adverb to estimate the reciprocal of the
NB.        condition number of a diagonalizable matrix
NB. pocon  adverb to estimate the reciprocal of the
NB.        condition number of a Hermitian (symmetric),
NB.        positive definite matrix
NB.
NB. Version: 1.0.0 2009-01-01
NB. Copyright: Igor Zhuravlov, igor at uic.dvgu.ru
NB. License: GNU GPL

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gecon
NB. Adverb to estimate the reciprocal of the condition number
NB. of a general matrix in either the 1-norm or the ∞-norm,
NB. using the LU factorization
NB.   P * L * U = A
NB.
NB. Syntax:
NB.   rcond=. norm gecon A
NB. where
NB.   A     - m×n-matrix
NB.   norm  - monadic verb to calculate matrix norm, either
NB.           norm1 or normi
NB.   rcond - the reciprocal of the condition number of the
NB.           matrix A in norm defined by verb 'norm',
NB.           0≤rcond≤1
NB.   m     ≥ 0
NB.   n     ≥ 0
NB.
NB. Application:
NB.
NB. Notes:
NB. - permutation is ignored since it doesn't changes rcond
NB.
NB. References:
NB.

gecon=: 1 : 0
  LU=. > {: getrfl1u y
  invA=. ((128!:1) utri LU) mp (((128!:1) &. |:) ltri1 LU)
  (% u y) % (u invA)
)

NB. ---------------------------------------------------------
hecon=: [:

NB. ---------------------------------------------------------
dicon=: [:

NB. ---------------------------------------------------------
pocon=: 1 : 0
  L=. potrfl y
  invA=. ((128!:1) ct L) mp (((128!:1) &. |:) L)
  (% u y) % (u invA)
)
