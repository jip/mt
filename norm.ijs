NB. norm.ijs
NB. Norms
NB.
NB. norm1     Magnitude-based 1-norm of vector or matrix
NB. normi     Magnitude-based ∞-norm of vector or matrix
NB. norm1t    Taxicab-based 1-norm of vector or matrix
NB. normit    Taxicab-based ∞-norm of vector or matrix
NB. norms     Square-based (Euclidean/Frobenius) norm of
NB.           vector or matrix
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

norm1=: | mocs           NB. 1-norm of vector or matrix
normi=: | mors           NB. ∞-norm of vector or matrix

NB. ---------------------------------------------------------
NB. Taxicab-based norms |Re(y)| + |Im(y)|

norm1t=: sorim mocs      NB. 1-norm of vector or matrix, for vector input emulates LAPACK's DASUM,DZASUM
normit=: sorim mors      NB. ∞-norm of vector or matrix

NB. ---------------------------------------------------------
NB. Square-based Euclidean (Frobenius) norm of vector
NB. (matrix), emulates LAPACK's DZNRM2 for vector input

norms=: (((((+/^:_) &.: *:) @: %) * ]) (>./^:_)) @: | @: +.
