NB. Condition number
NB.
NB. con  Conj. to make verb estimating the reciprocal of the
NB.      condition number of a matrix in a given norm
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

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
NB.   norm  - monadic verb to calculate matrix norm
NB.   inv   - monadic verb to inverse square matrix
NB.   A     - n×n-matrix
NB.   rcond - the reciprocal of the condition number of
NB.           matrix A in norm defined by verb 'norm',
NB.           0≤rcond≤1
NB.   n     ≥ 0
NB.
NB. Application:
NB. - verb to estimate rcond(y) of a general square matrix y
NB.   in 1-norm:
NB.     gecon=: norm1 con (getriul1p@getrful1p)
NB. - noun, estimated rcond(A) of a Hermitian (symmetric)
NB.   positive definite matrix A in ∞-norm:
NB.     rcondA=. normi con (potri@potrf) A
NB. - noun, estimated rcond(Q) of a unitary (orthogonal)
NB.   matrix Q in Frobenius-norm
NB.     rcondQ=. norms con ct Q
NB.
NB. Notes:
NB. - currently an expensive approach of exact calculation
NB.   is implemented

con=: 2 : '(* & (% @ u)) v'
