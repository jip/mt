NB. Norms
NB.
NB. norm1   Magnitude-based 1-norm of vector (matrix)
NB. normi   Magnitude-based ∞-norm of vector (matrix)
NB. norm1t  Taxicab-based 1-norm of vector (matrix)
NB. normit  Taxicab-based ∞-norm of vector (matrix)
NB. norms   Square-based Euclidean (Frobenius) norm of vector
NB.         (matrix)
NB.
NB. Version: 0.6.8 2010-11-30
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
NB. - norm1 implements LAPACK's DZSUM1, DLANSY('1'),
NB.   xLANGE('1'), ZLANHE('1')
NB. - to force norm1 act like any of: DLANSB('1'),
NB.   DLANST('1'), xLANGB('1'), xLANGT('1'), xLANHS('1'),
NB.   xLANTB('1'), xLANTR('1'), ZLANHB('1'), ZLANHT('1'),-
NB.   extraneous values in matrix must be zeroed
NB. - normi implements LAPACK's DLANSY('i'), xLANGE('i'),
NB.   ZLANHE('i')
NB. - to force normi act like any of: DLANSB('i'),
NB.   DLANST('i'), xLANGB('i'), xLANGT('i'), xLANHS('i'),
NB.   xLANTB('i'), xLANTR('i'), ZLANHB('i'), ZLANHT('i'),-
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
NB.   DLANSY('f'), xLANGE('f'), ZLANHE('f')
NB. - to force norms act like any of: DLANSB('f'),
NB.   DLANST('f'), xLANGB('f'), xLANGT('f'), xLANHS('f'),
NB.   xLANTB('f'), xLANTR('f'), ZLANHB('f'), ZLANHT('f'),-
NB.   extraneous values in matrix must be zeroed

norms=: (((((+/^:_) &.: *:) @: %) * ]) (>./^:_)) @: | @: +.
