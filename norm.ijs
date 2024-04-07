NB. Norms
NB.
NB. norm1      Magnitude-based 1-norm of vector (matrix)
NB. norm1c     Magnitude-based 1-norm of vector (matrix
NB.            columns)
NB. norm1r     Magnitude-based 1-norm of vector (matrix rows)
NB. normi      Magnitude-based ∞-norm of vector (matrix)
NB. normic     Magnitude-based ∞-norm of vector (matrix
NB.            columns)
NB. normir     Magnitude-based ∞-norm of vector (matrix rows)
NB. normm      Magnitude-based max of modules of elements
NB. norm1t     Taxicab-based 1-norm of vector (matrix)
NB. norm1tc    Taxicab-based 1-norm of vector (matrix
NB.            columns)
NB. norm1tr    Taxicab-based 1-norm of vector (matrix rows)
NB. normit     Taxicab-based ∞-norm of vector (matrix)
NB. normitc    Taxicab-based ∞-norm of vector (matrix
NB.            columns)
NB. normitr    Taxicab-based ∞-norm of vector (matrix rows)
NB. normmt     Taxicab-based max of modules of elements
NB. norms      Square-based Euclidean (Frobenius) norm of
NB.            vector (matrix)
NB. normsc     Square-based Euclidean norm of matrix columns
NB. normsr     Square-based Euclidean norm of matrix rows
NB.
NB. testnormm  Test magnitude-based norm computing verbs
NB. testnormt  Test taxicab-based norm computing verbs
NB. testnorms  Test square-based norm computing verbs
NB. testnorm   Adv. to make verb to test norm computing
NB.            algorithms by matrix of generator and shape
NB.            given
NB.
NB. Copyright 2010,2011,2013,2017,2018,2020,2021,2023,2024
NB.           Igor Zhuravlov
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

NB. =========================================================
NB. Concepts
NB.
NB. Notation:
NB.   E-norm - Euclidean norm (a.k.a. 2-norm) of vector
NB.   F-norm - Frobenius norm of matrix

NB. =========================================================
NB. Configuration

coclass 'mt'

NB. =========================================================
NB. Local definitions

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
NB.   0  -: norm1  ''
NB.   '' -: norm1c ''
NB.   '' -: norm1r ''
NB.   0  -: normi  ''
NB.   '' -: normic ''
NB.   '' -: normir ''
NB.   0  -: normm  ''
NB.   isnan norm1  _.
NB.   isnan norm1c _.
NB.   isnan norm1r _.
NB.   isnan normi  _.
NB.   isnan normic _.
NB.   isnan normir _.
NB.   isnan normm  _.
NB.   isnan norm1  _.  2
NB.   isnan norm1c _.  2
NB.   isnan norm1r _.  2
NB.   isnan normi  _.  2
NB.   isnan normic _.  2
NB.   isnan normir _.  2
NB.   isnan normm  _.  2
NB.   isnan norm1  _. _.
NB.   isnan norm1c _. _.
NB.   isnan norm1r _. _.
NB.   isnan normi  _. _.
NB.   isnan normic _. _.
NB.   isnan normir _. _.
NB.   isnan normm  _. _.
NB.   1 -: # ~. (norm1 , norm1c , norm1r)       10 ?@$ 0
NB.   1 -: # ~. (norm1 , norm1c , norm1r) j./ 2 10 ?@$ 0
NB.   1 -: # ~. (normi , normic , normir)       10 ?@$ 0
NB.   1 -: # ~. (normi , normic , normir) j./ 2 10 ?@$ 0
NB.   1 -: # ~. norm1 "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. norm1 "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. norm1c"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. norm1c"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. norm1r"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. norm1r"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. normi "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. normi "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. normic"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. normic"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. normir"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. normir"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
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
NB.   0     -: normm  0 0 $ 0
NB.   0     -: normm  0 3 $ 0
NB.   0     -: normm  3 0 $ 0
NB.            isnan norm1  _. (< 1 1)} 3 3 ?@$ 0
NB.   1 1 1 -: isnan norm1c _. (< 1 1)} 3 3 ?@$ 0
NB.   1 1 1 -: isnan norm1r _. (< 1 1)} 3 3 ?@$ 0
NB.            isnan normi  _. (< 1 1)} 3 3 ?@$ 0
NB.   1 1 1 -: isnan normic _. (< 1 1)} 3 3 ?@$ 0
NB.   1 1 1 -: isnan normir _. (< 1 1)} 3 3 ?@$ 0
NB.            isnan normm  _. (< 1 1)} 3 3 ?@$ 0
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

norm1=:   >./  @(+/!.0   )@:|    `     nan @.(isnan@<)` 0:     @.(0 e. $) : [:  NB. 1-norm of vector (matrix)
norm1c=:         +/!.0    @:|    `(c # nan)@.(isnan@<)`(0 #~ c)@.(0 e. $) : [:  NB. 1-norm of vector (matrix columns)
norm1r=:         +/!.0" 1 @:|    `(# # nan)@.(isnan@<)`(0 #~ #)@.(0 e. $) : [:  NB. 1-norm of vector (matrix rows)

normi=:   >./  @(+/!.0"_1)@:|    `     nan @.(isnan@<)` 0:     @.(0 e. $) : [:  NB. ∞-norm of vector (matrix)
normic=:  >./             @:|    `(c # nan)@.(isnan@<)`(0 #~ c)@.(0 e. $) : [:  NB. ∞-norm of vector (matrix columns)
normir=:  >./"1           @:|    `(# # nan)@.(isnan@<)`(0 #~ #)@.(0 e. $) : [:  NB. ∞-norm of vector (matrix rows)

normm=:   >./  @,         @:|    `     nan @.(isnan@<)` 0:     @.(0 e. $) : [:  NB. max of modules of elements of vector (matrix)

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
NB.   0  -: norm1t  ''
NB.   '' -: norm1tc ''
NB.   '' -: norm1tr ''
NB.   0  -: normit  ''
NB.   '' -: normitc ''
NB.   '' -: normitr ''
NB.   0  -: normmt  ''
NB.   isnan norm1t  _.
NB.   isnan norm1tc _.
NB.   isnan norm1tr _.
NB.   isnan normit  _.
NB.   isnan normitc _.
NB.   isnan normitr _.
NB.   isnan normmt  _.
NB.   isnan norm1t  _.  2
NB.   isnan norm1tc _.  2
NB.   isnan norm1tr _.  2
NB.   isnan normit  _.  2
NB.   isnan normitc _.  2
NB.   isnan normitr _.  2
NB.   isnan normmt  _.  2
NB.   isnan norm1t  _. _.
NB.   isnan norm1tc _. _.
NB.   isnan norm1tr _. _.
NB.   isnan normit  _. _.
NB.   isnan normitc _. _.
NB.   isnan normitr _. _.
NB.   isnan normmt  _. _.
NB.   1 -: # ~. (norm1t , norm1tc , norm1tr)       10 ?@$ 0
NB.   1 -: # ~. (norm1t , norm1tc , norm1tr) j./ 2 10 ?@$ 0
NB.   1 -: # ~. (normit , normitc , normitr)       10 ?@$ 0
NB.   1 -: # ~. (normit , normitc , normitr) j./ 2 10 ?@$ 0
NB.   1 -: # ~. norm1t "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. norm1t "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. norm1tc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. norm1tc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. norm1tr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. norm1tr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. normit "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. normit "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. normitc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. normitc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. normitr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. normitr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) j./ 2 10 ?@$ 0
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
NB.   0     -: normmt  0 0 $ 0
NB.   0     -: normmt  0 3 $ 0
NB.   0     -: normmt  3 0 $ 0
NB.            isnan norm1t  _. (< 1 1)} 3 3 ?@$ 0
NB.   1 1 1 -: isnan norm1tc _. (< 1 1)} 3 3 ?@$ 0
NB.   1 1 1 -: isnan norm1tr _. (< 1 1)} 3 3 ?@$ 0
NB.            isnan normit  _. (< 1 1)} 3 3 ?@$ 0
NB.   1 1 1 -: isnan normitc _. (< 1 1)} 3 3 ?@$ 0
NB.   1 1 1 -: isnan normitr _. (< 1 1)} 3 3 ?@$ 0
NB.            isnan normmt  _. (< 1 1)} 3 3 ?@$ 0
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

norm1t=:  >./  @(+/!.0   )@:sorim`     nan @.(isnan@<)` 0:     @.(0 e. $) : [:  NB. 1-norm of vector (matrix)
norm1tc=:        +/!.0    @:sorim`(c # nan)@.(isnan@<)`(0 #~ c)@.(0 e. $) : [:  NB. 1-norm of vector (matrix columns)
norm1tr=:        +/!.0" 1 @:sorim`(# # nan)@.(isnan@<)`(0 #~ #)@.(0 e. $) : [:  NB. 1-norm of vector (matrix rows)

normit=:  >./  @(+/!.0"_1)@:sorim`     nan @.(isnan@<)` 0:     @.(0 e. $) : [:  NB. ∞-norm of vector (matrix)
normitc=: >./             @:sorim`(c # nan)@.(isnan@<)`(0 #~ c)@.(0 e. $) : [:  NB. ∞-norm of vector (matrix columns)
normitr=: >./"1           @:sorim`(# # nan)@.(isnan@<)`(0 #~ #)@.(0 e. $) : [:  NB. ∞-norm of vector (matrix rows)

normmt=:  >./  @,         @:sorim`     nan @.(isnan@<)` 0:     @.(0 e. $) : [:  NB. max of modules of elements of vector (matrix)

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
NB.   0  -: norms  ''
NB.   '' -: normsc ''
NB.   '' -: normsr ''
NB.   isnan norms  _.
NB.   isnan normsc _.
NB.   isnan normsr _.
NB.   isnan norms  _.  2
NB.   isnan normsc _.  2
NB.   isnan normsr _.  2
NB.   isnan norms  _. _.
NB.   isnan normsc _. _.
NB.   isnan normsr _. _.
NB.   1 -: # ~. (norms , normsc , normsr)       10 ?@$ 0
NB.   1 -: # ~. (norms , normsc , normsr) j./ 2 10 ?@$ 0
NB.   1 -: # ~. norms "1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. norms "1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. normsc"1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. normsc"1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #) j./ 2 10 ?@$ 0
NB.   1 -: # ~. normsr"2 (((C."1~ ?~) , (C."1~ ?~) ,: (C."1~ ?~)) #)       10 ?@$ 0
NB.   1 -: # ~. normsr"2 (((C."1~ ?~) , (C."1~ ?~) ,: (C."1~ ?~)) #) j./ 2 10 ?@$ 0
NB.   0     -: norms  0 0 $ 0
NB.   0     -: norms  0 3 $ 0
NB.   0     -: norms  3 0 $ 0
NB.   ''    -: normsc 0 0 $ 0
NB.   0 0 0 -: normsc 0 3 $ 0
NB.   ''    -: normsc 3 0 $ 0
NB.   ''    -: normsr 0 0 $ 0
NB.   ''    -: normsr 0 3 $ 0
NB.   0 0 0 -: normsr 3 0 $ 0
NB.            isnan norms  _. (< 1 1)} 3 3 ?@$ 0
NB.   1 1 1 -: isnan normsc _. (< 1 1)} 3 3 ?@$ 0
NB.   1 1 1 -: isnan normsr _. (< 1 1)} 3 3 ?@$ 0
NB.   (norms"1@|: -: normsc)       10 10 ?@$ 0
NB.   (norms"1@|: -: normsc) j./ 2 10 10 ?@$ 0
NB.   (norms"1    -: normsr)       10 10 ?@$ 0
NB.   (norms"1    -: normsr) j./ 2 10 10 ?@$ 0
NB.
NB. Notes:
NB. - norms models BLAS' DNRM2, DZNRM2, xLASSQ, LAPACK's
NB.   DLANSY('f'), xLANGE('f'), ZLANHE('f')
NB. - to force norms act like any of: DLANSB('f'),
NB.   DLANST('f'), xLANGB('f'), xLANGT('f'), xLANHS('f'),
NB.   xLANTB('f'), xLANTR('f'), ZLANHB('f'), ZLANHT('f'),-
NB.   extraneous values in matrix must be zeroed

norms=:  ((+/!.0&.:*:  @: (% :: 1:  )    * ]) >./  )@,@:|@(+.            ^:(JCMPX = 3!:0))`     nan @.(isnan@<)                    : [:  NB. E-norm of vector (F-norm of matrix)
normsc=: ((+/!.0&.:*:  @:((% :: 1:"0)"1) * ]) >./  )  @:|@((9&o. , 11&o.)^:(JCMPX = 3!:0))`(c # nan)@.(isnan@<)`(0 #~ c)@.(0 e. $) : [:  NB. E-norm of matrix columns
normsr=: ((+/!.0&.:*:"1@: (% :: 1:"0)    * ]) >./"1)  @:|@((,"2@:+.     )^:(JCMPX = 3!:0))`(# # nan)@.(isnan@<)`(0 #~ #)@.(0 e. $) : [:  NB. E-norm of matrix rows

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testnormm
NB.
NB. Description:
NB.   Test magnitude-based norm computing verbs:
NB.   - norm1x (math/mt addon)
NB.   - normix (math/mt addon)
NB.   - normm (math/mt addon)
NB.   by m×n-matrix
NB.
NB. Syntax:
NB.   log=. testnormm A
NB. where
NB.   A   - m×n-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs

testnormm=: 3 : 0
  rcond=. nan`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  log=.          ('norm1'  tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('norm1c' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('norm1r' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('normi'  tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('normic' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('normir' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('normm'  tmonad (]`]`(rcond"_)`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testnormt
NB.
NB. Description:
NB.   Test taxicab-based norm computing verbs:
NB.   - norm1tx (math/mt addon)
NB.   - normitx (math/mt addon)
NB.   - normmt (math/mt addon)
NB.   by m×n-matrix
NB.
NB. Syntax:
NB.   log=. testnormt A
NB. where
NB.   A   - m×n-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs

testnormt=: 3 : 0
  rcond=. nan`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  log=.          ('norm1t'  tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('norm1tc' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('norm1tr' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('normit'  tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('normitc' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('normitr' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('normmt'  tmonad (]`]`(rcond"_)`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testnorms
NB.
NB. Description:
NB.   Test square-based norm computing verbs:
NB.   - normsx (math/mt addon)
NB.   by m×n-matrix
NB.
NB. Syntax:
NB.   log=. testnorms A
NB. where
NB.   A   - m×n-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs

testnorms=: 3 : 0
  rcond=. nan`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  log=.          ('norms'  tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('normsc' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('normsr' tmonad (]`]`(rcond"_)`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testnorm
NB.
NB. Description:
NB.   Adv. to make verb to test norm computing algorithms by
NB.   matrix of generator and shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testnorm) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testnorm_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testnorm_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testnorm_mt_ 150 200

testnorm=: 1 : '(testnorms_mt_ ,&.>~ testnormt_mt_ ,&.>~ testnormm_mt_)@u'
