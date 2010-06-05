NB. Norms
NB.
NB. norm1   Magnitude-based 1-norm of vector (matrix)
NB. normi   Magnitude-based ∞-norm of vector (matrix)
NB. norm1t  Taxicab-based 1-norm of vector (matrix)
NB. normit  Taxicab-based ∞-norm of vector (matrix)
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
NB.
NB. Notes:
NB. - norm1 implements LAPACK's DZSUM1,xLANGE('1')
NB. - to force norm1 act like any of: xLANHE('1') xLANHS('1')
NB.   xLANHT('1') xLANST('1') xLANTR('1'), extraneous values
NB.   in matrix must be zeroed
NB. - normi implements LAPACK's xLANGE('i')
NB. - to force normi act like any of: xLANHE('i') xLANHS('i')
NB.   xLANHT('i') xLANST('i') xLANTR('i'), extraneous values
NB.   in matrix must be zeroed

norm1=: | mocs           NB. 1-norm of vector (matrix)
normi=: | mors           NB. ∞-norm of vector (matrix)

NB. ---------------------------------------------------------
NB. Taxicab-based norms |Re(y)| + |Im(y)|
NB.
NB. Notes:
NB. - norm1t implements BLAS's DASUM,DZASUM

norm1t=: sorim mocs      NB. 1-norm of vector (matrix)
normit=: sorim mors      NB. ∞-norm of vector (matrix)

NB. ---------------------------------------------------------
NB. Square-based Euclidean (Frobenius) norm of vector
NB. (matrix) |y|^2
NB.
NB. Notes:
NB. - implements BLAS's DZNRM2 and partially xLASSQ, LAPACK's
NB.   xLANGE('f')
NB. - to force norms act like any of: xLANHE('f') xLANHS('f')
NB.   xLANHT('f') xLANST('f') xLANTR('f'), extraneous values
NB.   in matrix must be zeroed

norms=: (((((+/^:_) &.: *:) @: %) * ]) (>./^:_)) @: | @: +.
