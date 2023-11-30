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
NB. normm    Magnitude-based max of modules of elements
NB. norm1t   Taxicab-based 1-norm of vector (matrix)
NB. norm1tc  Taxicab-based 1-norm of vector (matrix columns)
NB. norm1tr  Taxicab-based 1-norm of vector (matrix rows)
NB. normit   Taxicab-based ∞-norm of vector (matrix)
NB. normitc  Taxicab-based ∞-norm of vector (matrix columns)
NB. normitr  Taxicab-based ∞-norm of vector (matrix rows)
NB. normmt   Taxicab-based max of modules of elements
NB. norms    Square-based Euclidean (Frobenius) norm of
NB.          vector (matrix)
NB. normsc   Square-based Euclidean norm of matrix columns
NB. normsr   Square-based Euclidean norm of matrix rows
NB.
NB. Version: 0.13.1 2021-06-06
NB.
NB. Copyright 2010-2021 Igor Zhuravlov
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

csum=: +/"2@:    NB. vector: sum of, matrix: column sums
rsum=: +/"1@:    NB. vector: sum of, matrix: row sums

cmax=: maxc"2@:  NB. vector: max of, matrix: column maximums
rmax=: max "1@:  NB. vector: max of, matrix: row maximums

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. norm1
NB. norm1c
NB. norm1r
NB. normi
NB. normic
NB. normir
NB. normm
NB.
NB. Description:
NB.   Magnitude-based norms |y|
NB.
NB. Assertions:
NB.   0     -: norm1  0 0 $ 0
NB.   0     -: norm1  0 3 $ 0
NB.   0     -: norm1  3 0 $ 0
NB.   ''    -: norm1c 0 0 $ 0
NB.   0 0 0 -: norm1c 0 3 $ 0
NB.   ''    -: norm1c 3 0 $ 0
NB.   ''    -: norm1r 0 0 $ 0
NB.   ''    -: norm1r 0 3 $ 0
NB.   0 0 0 -: norm1r 3 0 $ 0
NB.   0     -: normi  0 0 $ 0
NB.   0     -: normi  0 3 $ 0
NB.   0     -: normi  3 0 $ 0
NB.   ''    -: normic 0 0 $ 0
NB.   0 0 0 -: normic 0 3 $ 0
NB.   ''    -: normic 3 0 $ 0
NB.   ''    -: normir 0 0 $ 0
NB.   ''    -: normir 0 3 $ 0
NB.   0 0 0 -: normir 3 0 $ 0
NB.   (norm1"1@|: -: norm1c)       10 10 ?@$ 0
NB.   (norm1"1@|: -: norm1c) j./ 2 10 10 ?@$ 0
NB.   (norm1"1    -: norm1r)       10 10 ?@$ 0
NB.   (norm1"1    -: norm1r) j./ 2 10 10 ?@$ 0
NB.   (normi"1@|: -: normic)       10 10 ?@$ 0
NB.   (normi"1@|: -: normic) j./ 2 10 10 ?@$ 0
NB.   (normi"1    -: normir)       10 10 ?@$ 0
NB.   (normi"1    -: normir) j./ 2 10 10 ?@$ 0
NB.
NB. Notes:
NB. - norm1 implements LAPACK's DASUM, DZSUM1, DLANSY('1'),
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
NB. - normm implements LAPACK's DLANSY('m'), xLANGE('m'),
NB.   ZLANHE('m')
NB. - to force normm act like any of: DLANSB('m'),
NB.   DLANST('m'), xLANGB('m'), xLANGT('m'), xLANHS('m'),
NB.   xLANTB('m'), xLANTR('m'), ZLANHB('m'), ZLANHT('m'),-
NB.   extraneous values in matrix must be zeroed

norm1=:  norm1c      (max@)    NB. 1-norm of vector (matrix)
norm1c=: | csum                NB. 1-norm of vector (matrix columns)
norm1r=: | rsum                NB. 1-norm of vector (matrix rows)

normi=:  | (+/"_1@:) (max@)    NB. ∞-norm of vector (matrix)
normic=: | cmax                NB. ∞-norm of vector (matrix columns)
normir=: | rmax                NB. ∞-norm of vector (matrix rows)

normm=:  >./@,@:|`0:@.(0 = #)  NB. max of modules of elements of vector (matrix)

NB. ---------------------------------------------------------
NB. norm1t
NB. norm1tc
NB. norm1tr
NB. normit
NB. normitc
NB. normitr
NB. normmt
NB.
NB. Description:
NB.   Taxicab-based norms |Re(y)| + |Im(y)|
NB.
NB. Assertions:
NB.   0     -: norm1t  0 0 $ 0
NB.   0     -: norm1t  0 3 $ 0
NB.   0     -: norm1t  3 0 $ 0
NB.   ''    -: norm1tc 0 0 $ 0
NB.   0 0 0 -: norm1tc 0 3 $ 0
NB.   ''    -: norm1tc 3 0 $ 0
NB.   ''    -: norm1tr 0 0 $ 0
NB.   ''    -: norm1tr 0 3 $ 0
NB.   0 0 0 -: norm1tr 3 0 $ 0
NB.   0     -: normit  0 0 $ 0
NB.   0     -: normit  0 3 $ 0
NB.   0     -: normit  3 0 $ 0
NB.   ''    -: normitc 0 0 $ 0
NB.   0 0 0 -: normitc 0 3 $ 0
NB.   ''    -: normitc 3 0 $ 0
NB.   ''    -: normitr 0 0 $ 0
NB.   ''    -: normitr 0 3 $ 0
NB.   0 0 0 -: normitr 3 0 $ 0
NB.   (norm1t"1@|: -: norm1tc)       10 10 ?@$ 0
NB.   (norm1t"1@|: -: norm1tc) j./ 2 10 10 ?@$ 0
NB.   (norm1t"1    -: norm1tr)       10 10 ?@$ 0
NB.   (norm1t"1    -: norm1tr) j./ 2 10 10 ?@$ 0
NB.   (normit"1@|: -: normitc)       10 10 ?@$ 0
NB.   (normit"1@|: -: normitc) j./ 2 10 10 ?@$ 0
NB.   (normit"1    -: normitr)       10 10 ?@$ 0
NB.   (normit"1    -: normitr) j./ 2 10 10 ?@$ 0
NB.
NB. Notes:
NB. - norm1t implements BLAS' DASUM, DZASUM

norm1t=:  norm1tc         (max@)    NB. 1-norm of vector (matrix)
norm1tc=: sorim csum                NB. 1-norm of vector (matrix columns)
norm1tr=: sorim rsum                NB. 1-norm of vector (matrix rows)

normit=:  sorim (+/"_1@:) (max@)    NB. ∞-norm of vector (matrix)
normitc=: sorim cmax                NB. ∞-norm of vector (matrix columns)
normitr=: sorim rmax                NB. ∞-norm of vector (matrix rows)

normmt=:  >./@,@:sorim`0:@.(0 = #)  NB. max of modules of elements of vector (matrix)

NB. ---------------------------------------------------------
NB. norms
NB. normsc
NB. normsr
NB.
NB. Description:
NB.   Square-based Euclidean (Frobenius) norm of vector
NB.   (matrix) sqrt(|y|^2)
NB.
NB. Assertions:
NB.   0     -: norms  0 0 $ 0
NB.   0     -: norms  0 3 $ 0
NB.   0     -: norms  3 0 $ 0
NB.   ''    -: normsc 0 0 $ 0
NB.   0 0 0 -: normsc 0 3 $ 0
NB.   ''    -: normsc 3 0 $ 0
NB.   ''    -: normsr 0 0 $ 0
NB.   ''    -: normsr 0 3 $ 0
NB.   0 0 0 -: normsr 3 0 $ 0
NB.   (norms"1@|: -: normsc)       10 10 ?@$ 0
NB.   (norms"1@|: -: normsc) j./ 2 10 10 ?@$ 0
NB.   (norms"1    -: normsr)       10 10 ?@$ 0
NB.   (norms"1    -: normsr) j./ 2 10 10 ?@$ 0
NB.
NB. Notes:
NB. - norms implements BLAS' DNRM2, DZNRM2 and models
NB.   xLASSQ, LAPACK's DLANSY('f'), xLANGE('f'), ZLANHE('f')
NB. - to force norms act like any of: DLANSB('f'),
NB.   DLANST('f'), xLANGB('f'), xLANGT('f'), xLANHS('f'),
NB.   xLANTB('f'), xLANTR('f'), ZLANHB('f'), ZLANHT('f'),-
NB.   extraneous values in matrix must be zeroed

norms=:  ((+/  &.:*:@: %    * ]) >./  )@,@:|                                                      @ (+.^:(JCMPX = 3!:0))  NB. E-norm of vector (F-norm of matrix)
normsc=: ((+/  &.:*:@:(%"1) * ]) >./  )  @:|`(((+/@,"2&.:*:@:% * ]) >./@,"2)@(0 2&|: :: ])@:|@:+.)@.     (JCMPX = 3!:0)   NB. E-norm of matrix columns
normsr=: ((+/"1&.:*:@: %    * ]) >./"1)  @:|`(((+/@,"2&.:*:@:% * ]) >./@,"2)              @:|@:+.)@.     (JCMPX = 3!:0)   NB. E-norm of matrix rows
