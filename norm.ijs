NB. Norms
NB.
NB. norm1   Magnitude-based 1-norm of vector (matrix)
NB. normi   Magnitude-based ∞-norm of vector (matrix)
NB. norm1t  Taxicab-based 1-norm of vector (matrix)
NB. normit  Taxicab-based ∞-norm of vector (matrix)
NB. norms   Square-based Euclidean (Frobenius) norm of vector
NB.         (matrix)
NB.
NB. Version: 0.6.0 2010-06-05
NB.
NB. Copyright 2010 Igor Zhuravlov
NB.
NB. This file is part of mt
NB.
NB. mt is free software: you can redistribute it and/or
NB. modify it under the terms of the GNU Lesser General
NB. Public License as published by the Free Software
NB. Foundation, either version 3 of the License, or (at your
NB. option) any later version.
NB.
NB. mt is distributed in the hope that it will be useful, but
NB. WITHOUT ANY WARRANTY; without even the implied warranty
NB. of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
NB. See the GNU Lesser General Public License for more
NB. details.
NB.
NB. You should have received a copy of the GNU Lesser General
NB. Public License along with mt. If not, see
NB. <http://www.gnu.org/licenses/>.

coclass 'mt'

NB. =========================================================
NB. Local definitions

mocs=: >./`0:@.(0=#) @ (+/   ) @:  NB. vector: sum of, matrix: max of column sums
mors=: >./`0:@.(0=#) @ (+/"_1) @:  NB. vector: max of, matrix: max of row sums

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Magnitude-based norms |y|
NB.
NB. Notes:
NB. - norm1 implements LAPACK's DZSUM1, xLANGE('1')
NB. - to force norm1 act like any of: xLANHE('1'),
NB.   xLANHS('1'), xLANHT('1'), xLANST('1'), xLANTR('1');
NB.   extraneous values in matrix must be zeroed
NB. - normi implements LAPACK's xLANGE('i')
NB. - to force normi act like any of: xLANHE('i'),
NB.   xLANHS('i'), xLANHT('i'), xLANST('i'), xLANTR('i');
NB.   extraneous values in matrix must be zeroed

norm1=: | mocs           NB. 1-norm of vector (matrix)
normi=: | mors           NB. ∞-norm of vector (matrix)

NB. ---------------------------------------------------------
NB. Taxicab-based norms |Re(y)| + |Im(y)|
NB.
NB. Notes:
NB. - norm1t implements BLAS's DASUM, DZASUM

norm1t=: sorim mocs      NB. 1-norm of vector (matrix)
normit=: sorim mors      NB. ∞-norm of vector (matrix)

NB. ---------------------------------------------------------
NB. Square-based Euclidean (Frobenius) norm of vector
NB. (matrix) |y|^2
NB.
NB. Notes:
NB. - implements BLAS's DZNRM2 and partially xLASSQ, LAPACK's
NB.   xLANGE('f')
NB. - to force norms act like any of: xLANHE('f'),
NB.   xLANHS('f'), xLANHT('f'), xLANST('f'), xLANTR('f');
NB.   extraneous values in matrix must be zeroed

norms=: (((((+/^:_) &.: *:) @: %) * ]) (>./^:_)) @: | @: +.
