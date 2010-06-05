NB. con.ijs
NB. Condition number
NB.
NB. con  conjunction to estimate the reciprocal of the
NB.      condition number of a matrix in a given norm
NB.
NB. Version: 1.0.0 2009-06-01
NB. Copyright: Igor Zhuravlov, igor at uic.dvgu.ru
NB. License: GNU GPL

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. con
NB. Conjunction to estimate the reciprocal of the condition
NB. number of a matrix in a given norm
NB.
NB. Syntax:
NB.   rcond=. norm con inv A
NB. where
NB.   norm  - monadic verb to calculate matrix norm, usually
NB.           norm1 or normi
NB.   inv   - monadic verb to inverse square matrix, usually
NB.           trtriu, trtriu1, trtril, trtril1, getri, hetri,
NB.           ditri or potri
NB.   A     - n×n-matrix
NB.   rcond - the reciprocal of the condition number of the
NB.           matrix A in norm defined by verb 'norm',
NB.           0≤rcond≤1
NB.   n     ≥ 0
NB.
NB. Applications:
NB.   NB. verb to estimate rcond(y) of a general square
NB.   NB. matrix y in 1-norm
NB.   gecon=: norm1 con getri
NB.   NB. noun, estimated rcond(A) of a Hermitian (symmetric)
NB.   NB. positive definite matrix A in ∞-norm
NB.   rcondA=. normi con potri A
NB.
NB. Notes:
NB. - currently an expensive approach of exact calculation
NB.   is implemented

con=: 2 : '(* & (% @ u)) v'
