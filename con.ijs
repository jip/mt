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
NB.
NB. Description:
NB.   Conj. to make verb estimating the reciprocal of the
NB.   condition number of a matrix in a given norm
NB.
NB. Syntax:
NB.   vapp=. norm con inv
NB. where
NB.   norm  - monadic verb to calculate norm of matrix, is
NB.           called as:
NB.             normA=. norm A
NB.   inv   - monadic verb to inverse square matrix, is
NB.           called as:
NB.             invA=. inv A
NB.   vapp  - monadic verb to calculate the reciprocal of the
NB.           condition number of a matrix in a given norm,
NB.           is called as:
NB.             rcond=. vapp A
NB.   A     - n×n-matrix
NB.
NB. Application:
NB. - verb to estimate rcond(y) of a general square matrix y
NB.   in 1-norm:
NB.     gecon=: norm1 con (getriul1p@getrful1p)
NB. - noun, an estimated rcond(A) of a Hermitian (symmetric)
NB.   positive definite matrix A in ∞-norm:
NB.     rcondA=. normi con (potril@potrfl) A
NB. - noun, estimated rcond(Q) of a unitary (orthogonal)
NB.   matrix Q in Frobenius-norm
NB.     rcondQ=. norms con ct Q
NB.
NB. TODO:
NB. - implement more practical norm-estimation approach

con=: 2 : '(* & (% @ u)) v'
