NB. Norms
NB.
NB. norm1    Magnitude-based 1-norm of vector (matrix)
NB. norm1c   Magnitude-based 1-norm of vector (matrix
NB.          columns)
NB. norm1r   Magnitude-based 1-norm of vector (matrix rows)
NB. normi    Magnitude-based ∞-norm of vector (matrix)
NB. normic   Magnitude-based ∞-norm of vector (matrix
NB.          columns)
NB. normir   Magnitude-based ∞-norm of vector (matrix rows)
NB. norm1t   Taxicab-based 1-norm of vector (matrix)
NB. norm1tc  Taxicab-based 1-norm of vector (matrix columns)
NB. norm1tr  Taxicab-based 1-norm of vector (matrix rows)
NB. normit   Taxicab-based ∞-norm of vector (matrix)
NB. normitc  Taxicab-based ∞-norm of vector (matrix columns)
NB. normitr  Taxicab-based ∞-norm of vector (matrix rows)
NB. norms    Square-based Euclidean (Frobenius) norm of
NB.          vector (matrix)
NB. normsc   Square-based Euclidean norm of matrix columns
NB. normsr   Square-based Euclidean norm of matrix rows
NB.
NB. Version: 0.10.0 2017-04-23
NB.
NB. Copyright 2010-2017 Igor Zhuravlov
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

max=: >./`0:@.(0 = #)  NB. max of, 0 for empty list

csum=: +/"2@:          NB. vector: sum of, matrix: column sums
rsum=: +/"1@:          NB. vector: sum of, matrix: row sums

cmax=: max"2@:         NB. vector: max of, matrix: column maximums
rmax=: max"1@:         NB. vector: max of, matrix: row maximums

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. norm1
NB. norm1c
NB. norm1r
NB. normi
NB. normic
NB. normir
NB.
NB. Description:
NB.   Magnitude-based norms |y|
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

norm1=:  | csum      (max@)  NB. 1-norm of vector (matrix)
norm1c=: | csum              NB. 1-norm of vector (matrix columns)
norm1r=: | rsum              NB. 1-norm of vector (matrix rows)

normi=:  | (+/"_1@:) (max@)  NB. ∞-norm of vector (matrix)
normic=: | cmax              NB. ∞-norm of vector (matrix columns)
normir=: | rmax              NB. ∞-norm of vector (matrix rows)

NB. ---------------------------------------------------------
NB. norm1t
NB. norm1tc
NB. norm1tr
NB. normit
NB. normitc
NB. normitr
NB.
NB. Description:
NB.   Taxicab-based norms |Re(y)| + |Im(y)|
NB.
NB. Notes:
NB. - norm1t implements BLAS's DASUM, DZASUM

norm1t=:  sorim csum      (max@)  NB. 1-norm of vector (matrix)
norm1tc=: sorim csum              NB. 1-norm of vector (matrix columns)
norm1tr=: sorim rsum              NB. 1-norm of vector (matrix rows)

normit=:  sorim (+/"_1@:) (max@)  NB. ∞-norm of vector (matrix)
normitc=: sorim cmax              NB. ∞-norm of vector (matrix columns)
normitr=: sorim rmax              NB. ∞-norm of vector (matrix rows)

NB. ---------------------------------------------------------
NB. norms
NB. normsc
NB. normsr
NB.
NB. Description:
NB.   Square-based Euclidean (Frobenius) norm of vector
NB.   (matrix) |y|^2
NB.
NB. Notes:
NB. - norms implements BLAS's DZNRM2 and partially xLASSQ,
NB.   LAPACK's DLANSY('f'), xLANGE('f'), ZLANHE('f')
NB. - to force norms act like any of: DLANSB('f'),
NB.   DLANST('f'), xLANGB('f'), xLANGT('f'), xLANHS('f'),
NB.   xLANTB('f'), xLANTR('f'), ZLANHB('f'), ZLANHT('f'),-
NB.   extraneous values in matrix must be zeroed

norms=:  ((+/          &.:*:@: %    * ]) >./           )@:|@,@:+.  NB. E-norm of vector (F-norm of matrix)
normsc=: ((+/"1@ (+/  )&.:*:@:(%"2) * ]) >./"1@ (>./  ))@:|  @:+.  NB. E-norm of matrix columns
normsr=: ((+/"1@:(+/"2)&.:*:@: %    * ]) >./"1@:(>./"2))@:|  @:+.  NB. E-norm of matrix rows
