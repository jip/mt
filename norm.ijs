NB. norm.ijs
NB. Norms
NB.
NB. norm1   Magnitude-based 1-norm of vector or matrix
NB. normi   Magnitude-based ∞-norm of vector or matrix
NB. norm1t  Taxicab-based 1-norm of vector or matrix
NB. normit  Taxicab-based ∞-norm of vector or matrix
NB. norms   Square-based Euclidean (Frobenius) norm of vector
NB.         (matrix)
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

mocs=: >./ @ (+/   ) @:  NB. vector: sum of, matrix: max of column sums
mors=: >./ @ (+/"_1) @:  NB. vector: max of, matrix: max of row sums

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Magnitude-based norms |y|

norm1=: | mocs           NB. 1-norm of vector or matrix, in case of vector input implements LAPACK's DZSUM1
normi=: | mors           NB. ∞-norm of vector or matrix

NB. ---------------------------------------------------------
NB. Taxicab-based norms |Re(y)| + |Im(y)|

norm1t=: sorim mocs      NB. 1-norm of vector or matrix, in case of vector input implements BLAS's DASUM,DZASUM
normit=: sorim mors      NB. ∞-norm of vector or matrix

NB. ---------------------------------------------------------
NB. Square-based Euclidean (Frobenius) norm of vector
NB. (matrix)
NB.
NB. Note:
NB. - in case of vector input implements BLAS's DZNRM2 and,
NB.   partially, xLASSQ

norms=: (((((+/^:_) &.: *:) @: %) * ]) (>./^:_)) @: | @: +.
