NB. Norms
NB.
NB. norm1       Magnitude-based 1-norm of vector (matrix)
NB. norm1c      Magnitude-based 1-norm of vector (matrix
NB.             columns)
NB. norm1r      Magnitude-based 1-norm of vector (matrix
NB.             rows)
NB. normi       Magnitude-based ∞-norm of vector (matrix)
NB. normic      Magnitude-based ∞-norm of vector (matrix
NB.             columns)
NB. normir      Magnitude-based ∞-norm of vector (matrix
NB.             rows)
NB. normm       Magnitude-based max of modules of elements
NB. norm1t      Taxicab-based 1-norm of vector (matrix)
NB. norm1tc     Taxicab-based 1-norm of vector (matrix
NB.             columns)
NB. norm1tr     Taxicab-based 1-norm of vector (matrix rows)
NB. normit      Taxicab-based ∞-norm of vector (matrix)
NB. normitc     Taxicab-based ∞-norm of vector (matrix
NB.             columns)
NB. normitr     Taxicab-based ∞-norm of vector (matrix rows)
NB. normmt      Taxicab-based max of modules of elements
NB. norms       Square-based Euclidean (Frobenius) norm of
NB.             vector (matrix)
NB. normsc      Square-based Euclidean norm of matrix columns
NB. normsr      Square-based Euclidean norm of matrix rows
NB.
NB. verifynorm  Verify norm verbs
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
NB. Concepts
NB.
NB. Notation:
NB.   E-norm - Euclidean norm (a.k.a. 2-norm) of vector
NB.   F-norm - Frobenius norm of matrix

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
NB. Verification suite

NB. ---------------------------------------------------------
NB. verifynorm
NB.
NB. Description:
NB.   Nilad to verify norm actors, output result to console
NB.   and return it
NB.
NB. Syntax:
NB.   'probed failed'=. verifynorm ''
NB. where
NB.   probed ≥ 0, assertions probed counter
NB.   failed ≥ 0, assertions failed counter

verifynorm=: 3 : 0

  NB. magnitude-based norms
  NB. - input is an empty vector
  res=.       fassert 0  -: norm1  ''
  res=. res , fassert '' -: norm1c ''
  res=. res , fassert '' -: norm1r ''
  res=. res , fassert 0  -: normi  ''
  res=. res , fassert '' -: normic ''
  res=. res , fassert '' -: normir ''
  res=. res , fassert 0  -: normm  ''
  NB. - input is a vector containing NaN
  res=. res , fassert isnan norm1  _.
  res=. res , fassert isnan norm1c _.
  res=. res , fassert isnan norm1r _.
  res=. res , fassert isnan normi  _.
  res=. res , fassert isnan normic _.
  res=. res , fassert isnan normir _.
  res=. res , fassert isnan normm  _.
  res=. res , fassert isnan norm1  _.  2
  res=. res , fassert isnan norm1c _.  2
  res=. res , fassert isnan norm1r _.  2
  res=. res , fassert isnan normi  _.  2
  res=. res , fassert isnan normic _.  2
  res=. res , fassert isnan normir _.  2
  res=. res , fassert isnan normm  _.  2
  res=. res , fassert isnan norm1  _. _.
  res=. res , fassert isnan norm1c _. _.
  res=. res , fassert isnan norm1r _. _.
  res=. res , fassert isnan normi  _. _.
  res=. res , fassert isnan normic _. _.
  res=. res , fassert isnan normir _. _.
  res=. res , fassert isnan normm  _. _.
  NB. - input is a vector containing ∞
  res=. res , fassert _ -: norm1  __
  res=. res , fassert _ -: norm1   _
  res=. res , fassert _ -: norm1c __
  res=. res , fassert _ -: norm1c  _
  res=. res , fassert _ -: norm1r __
  res=. res , fassert _ -: norm1r  _
  res=. res , fassert _ -: normi  __
  res=. res , fassert _ -: normi   _
  res=. res , fassert _ -: normic __
  res=. res , fassert _ -: normic  _
  res=. res , fassert _ -: normir __
  res=. res , fassert _ -: normir  _
  res=. res , fassert _ -: normm  __
  res=. res , fassert _ -: normm   _
  res=. res , fassert _ -: norm1  __  2
  res=. res , fassert _ -: norm1   _  2
  res=. res , fassert _ -: norm1c __  2
  res=. res , fassert _ -: norm1c  _  2
  res=. res , fassert _ -: norm1r __  2
  res=. res , fassert _ -: norm1r  _  2
  res=. res , fassert _ -: normi  __  2
  res=. res , fassert _ -: normi   _  2
  res=. res , fassert _ -: normic __  2
  res=. res , fassert _ -: normic  _  2
  res=. res , fassert _ -: normir __  2
  res=. res , fassert _ -: normir  _  2
  res=. res , fassert _ -: normm  __  2
  res=. res , fassert _ -: normm  __  2
  res=. res , fassert _ -: norm1  __ __
  res=. res , fassert _ -: norm1  __  _
  res=. res , fassert _ -: norm1   _ __
  res=. res , fassert _ -: norm1   _  _
  res=. res , fassert _ -: norm1c __ __
  res=. res , fassert _ -: norm1c __  _
  res=. res , fassert _ -: norm1c  _ __
  res=. res , fassert _ -: norm1c  _  _
  res=. res , fassert _ -: norm1r __ __
  res=. res , fassert _ -: norm1r __  _
  res=. res , fassert _ -: norm1r  _ __
  res=. res , fassert _ -: norm1r  _  _
  res=. res , fassert _ -: normi  __ __
  res=. res , fassert _ -: normi  __  _
  res=. res , fassert _ -: normi   _ __
  res=. res , fassert _ -: normi   _  _
  res=. res , fassert _ -: normic __ __
  res=. res , fassert _ -: normic __  _
  res=. res , fassert _ -: normic  _ __
  res=. res , fassert _ -: normic  _  _
  res=. res , fassert _ -: normir __ __
  res=. res , fassert _ -: normir __  _
  res=. res , fassert _ -: normir  _ __
  res=. res , fassert _ -: normir  _  _
  res=. res , fassert _ -: normm  __ __
  res=. res , fassert _ -: normm  __  _
  res=. res , fassert _ -: normm   _ __
  res=. res , fassert _ -: normm   _  _
  NB. - all norms are equal for vector input
  res=. res , fassert 1 -: # ~. (norm1 , norm1c , norm1r)  gemat           100
  res=. res , fassert 1 -: # ~. (norm1 , norm1c , norm1r) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. (normi , normic , normir)  gemat           100
  res=. res , fassert 1 -: # ~. (normi , normic , normir) (gemat j. gemat) 100
  NB. - norm should be resistable to permutations
  res=. res , fassert 1 -: # ~. norm1 "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. norm1 "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. norm1c"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. norm1c"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. norm1r"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. norm1r"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. normi "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. normi "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. normic"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. normic"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. normir"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. normir"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  NB. - input is an empty matrix
  res=. res , fassert 0     -: norm1  0 0 $ 0
  res=. res , fassert 0     -: norm1  0 3 $ 0
  res=. res , fassert 0     -: norm1  3 0 $ 0
  res=. res , fassert ''    -: norm1c 0 0 $ 0
  res=. res , fassert 0 0 0 -: norm1c 0 3 $ 0
  res=. res , fassert ''    -: norm1c 3 0 $ 0
  res=. res , fassert ''    -: norm1r 0 0 $ 0
  res=. res , fassert ''    -: norm1r 0 3 $ 0
  res=. res , fassert 0 0 0 -: norm1r 3 0 $ 0
  res=. res , fassert 0     -: normi  0 0 $ 0
  res=. res , fassert 0     -: normi  0 3 $ 0
  res=. res , fassert 0     -: normi  3 0 $ 0
  res=. res , fassert ''    -: normic 0 0 $ 0
  res=. res , fassert 0 0 0 -: normic 0 3 $ 0
  res=. res , fassert ''    -: normic 3 0 $ 0
  res=. res , fassert ''    -: normir 0 0 $ 0
  res=. res , fassert ''    -: normir 0 3 $ 0
  res=. res , fassert 0 0 0 -: normir 3 0 $ 0
  res=. res , fassert 0     -: normm  0 0 $ 0
  res=. res , fassert 0     -: normm  0 3 $ 0
  res=. res , fassert 0     -: normm  3 0 $ 0
  NB. - input is a matrix containing NaN
  res=. res , fassert          isnan norm1  _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 1 1 1 -: isnan norm1c _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 1 1 1 -: isnan norm1r _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert          isnan normi  _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 1 1 1 -: isnan normic _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 1 1 1 -: isnan normir _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert          isnan normm  _. (< 1 1)} 3 3 ?@$ 0
  NB. - input is a matrix containing ∞
  res=. res , fassert         _ = norm1  __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert         _ = norm1   _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = norm1c __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = norm1c  _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = norm1r __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = norm1r  _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert         _ = normi  __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert         _ = normi   _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normic __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normic  _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normir __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normir  _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert         _ = normm  __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert         _ = normm   _ (< 1 1)} 3 3 ?@$ 0
  NB. - row- and and column-wise norms are related
  res=. res , fassert (norm1"1@|: -: norm1c)  gemat           100 100
  res=. res , fassert (norm1"1@|: -: norm1c) (gemat j. gemat) 100 100
  res=. res , fassert (norm1"1    -: norm1r)  gemat           100 100
  res=. res , fassert (norm1"1    -: norm1r) (gemat j. gemat) 100 100
  res=. res , fassert (normi"1@|: -: normic)  gemat           100 100
  res=. res , fassert (normi"1@|: -: normic) (gemat j. gemat) 100 100
  res=. res , fassert (normi"1    -: normir)  gemat           100 100
  res=. res , fassert (normi"1    -: normir) (gemat j. gemat) 100 100

  NB. taxicab-based norms
  NB. - input is an empty vector
  res=. res , fassert 0  -: norm1t  ''
  res=. res , fassert '' -: norm1tc ''
  res=. res , fassert '' -: norm1tr ''
  res=. res , fassert 0  -: normit  ''
  res=. res , fassert '' -: normitc ''
  res=. res , fassert '' -: normitr ''
  res=. res , fassert 0  -: normmt  ''
  NB. - input is a vector containing NaN
  res=. res , fassert isnan norm1t  _.
  res=. res , fassert isnan norm1tc _.
  res=. res , fassert isnan norm1tr _.
  res=. res , fassert isnan normit  _.
  res=. res , fassert isnan normitc _.
  res=. res , fassert isnan normitr _.
  res=. res , fassert isnan normmt  _.
  res=. res , fassert isnan norm1t  _.  2
  res=. res , fassert isnan norm1tc _.  2
  res=. res , fassert isnan norm1tr _.  2
  res=. res , fassert isnan normit  _.  2
  res=. res , fassert isnan normitc _.  2
  res=. res , fassert isnan normitr _.  2
  res=. res , fassert isnan normmt  _.  2
  res=. res , fassert isnan norm1t  _. _.
  res=. res , fassert isnan norm1tc _. _.
  res=. res , fassert isnan norm1tr _. _.
  res=. res , fassert isnan normit  _. _.
  res=. res , fassert isnan normitc _. _.
  res=. res , fassert isnan normitr _. _.
  res=. res , fassert isnan normmt  _. _.
  NB. - input is a vector containing ∞
  res=. res , fassert _ -: norm1t  __
  res=. res , fassert _ -: norm1t   _
  res=. res , fassert _ -: norm1tc __
  res=. res , fassert _ -: norm1tc  _
  res=. res , fassert _ -: norm1tr __
  res=. res , fassert _ -: norm1tr  _
  res=. res , fassert _ -: normit  __
  res=. res , fassert _ -: normit   _
  res=. res , fassert _ -: normitc __
  res=. res , fassert _ -: normitc  _
  res=. res , fassert _ -: normitr __
  res=. res , fassert _ -: normitr  _
  res=. res , fassert _ -: normmt  __
  res=. res , fassert _ -: normmt   _
  res=. res , fassert _ -: norm1t  __  2
  res=. res , fassert _ -: norm1t   _  2
  res=. res , fassert _ -: norm1tc __  2
  res=. res , fassert _ -: norm1tc  _  2
  res=. res , fassert _ -: norm1tr __  2
  res=. res , fassert _ -: norm1tr  _  2
  res=. res , fassert _ -: normit  __  2
  res=. res , fassert _ -: normit   _  2
  res=. res , fassert _ -: normitc __  2
  res=. res , fassert _ -: normitc  _  2
  res=. res , fassert _ -: normitr __  2
  res=. res , fassert _ -: normitr  _  2
  res=. res , fassert _ -: normmt  __  2
  res=. res , fassert _ -: normmt  __  2
  res=. res , fassert _ -: norm1t  __ __
  res=. res , fassert _ -: norm1t  __  _
  res=. res , fassert _ -: norm1t   _ __
  res=. res , fassert _ -: norm1t   _  _
  res=. res , fassert _ -: norm1tc __ __
  res=. res , fassert _ -: norm1tc __  _
  res=. res , fassert _ -: norm1tc  _ __
  res=. res , fassert _ -: norm1tc  _  _
  res=. res , fassert _ -: norm1tr __ __
  res=. res , fassert _ -: norm1tr __  _
  res=. res , fassert _ -: norm1tr  _ __
  res=. res , fassert _ -: norm1tr  _  _
  res=. res , fassert _ -: normit  __ __
  res=. res , fassert _ -: normit  __  _
  res=. res , fassert _ -: normit   _ __
  res=. res , fassert _ -: normit   _  _
  res=. res , fassert _ -: normitc __ __
  res=. res , fassert _ -: normitc __  _
  res=. res , fassert _ -: normitc  _ __
  res=. res , fassert _ -: normitc  _  _
  res=. res , fassert _ -: normitr __ __
  res=. res , fassert _ -: normitr __  _
  res=. res , fassert _ -: normitr  _ __
  res=. res , fassert _ -: normitr  _  _
  res=. res , fassert _ -: normmt  __ __
  res=. res , fassert _ -: normmt  __  _
  res=. res , fassert _ -: normmt   _ __
  res=. res , fassert _ -: normmt   _  _
  NB. - all norms are equal for vector input
  res=. res , fassert 1 -: # ~. (norm1t , norm1tc , norm1tr)  gemat           100
  res=. res , fassert 1 -: # ~. (norm1t , norm1tc , norm1tr) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. (normit , normitc , normitr)  gemat           100
  res=. res , fassert 1 -: # ~. (normit , normitc , normitr) (gemat j. gemat) 100
  NB. - norm should be resistable to permutations
  res=. res , fassert 1 -: # ~. norm1t "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. norm1t "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. norm1tc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. norm1tc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. norm1tr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. norm1tr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. normit "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. normit "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. normitc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. normitc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. normitr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. normitr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
  NB. - input is an empty matrix
  res=. res , fassert 0     -: norm1t  0 0 $ 0
  res=. res , fassert 0     -: norm1t  0 3 $ 0
  res=. res , fassert 0     -: norm1t  3 0 $ 0
  res=. res , fassert ''    -: norm1tc 0 0 $ 0
  res=. res , fassert 0 0 0 -: norm1tc 0 3 $ 0
  res=. res , fassert ''    -: norm1tc 3 0 $ 0
  res=. res , fassert ''    -: norm1tr 0 0 $ 0
  res=. res , fassert ''    -: norm1tr 0 3 $ 0
  res=. res , fassert 0 0 0 -: norm1tr 3 0 $ 0
  res=. res , fassert 0     -: normit  0 0 $ 0
  res=. res , fassert 0     -: normit  0 3 $ 0
  res=. res , fassert 0     -: normit  3 0 $ 0
  res=. res , fassert ''    -: normitc 0 0 $ 0
  res=. res , fassert 0 0 0 -: normitc 0 3 $ 0
  res=. res , fassert ''    -: normitc 3 0 $ 0
  res=. res , fassert ''    -: normitr 0 0 $ 0
  res=. res , fassert ''    -: normitr 0 3 $ 0
  res=. res , fassert 0 0 0 -: normitr 3 0 $ 0
  res=. res , fassert 0     -: normmt  0 0 $ 0
  res=. res , fassert 0     -: normmt  0 3 $ 0
  res=. res , fassert 0     -: normmt  3 0 $ 0
  NB. - input is a matrix containing NaN
  res=. res , fassert          isnan norm1t  _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 1 1 1 -: isnan norm1tc _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 1 1 1 -: isnan norm1tr _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert          isnan normit  _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 1 1 1 -: isnan normitc _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 1 1 1 -: isnan normitr _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert          isnan normmt  _. (< 1 1)} 3 3 ?@$ 0
  NB. - input is a matrix containing ∞
  res=. res , fassert         _ = norm1t  __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert         _ = norm1t   _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = norm1tc __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = norm1tc  _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = norm1tr __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = norm1tr  _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert         _ = normit  __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert         _ = normit   _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normitc __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normitc  _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normitr __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normitr  _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert         _ = normmt  __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert         _ = normmt   _ (< 1 1)} 3 3 ?@$ 0
  NB. - row- and and column-wise norms are related
  res=. res , fassert (norm1t"1@|: -: norm1tc)  gemat           100 100
  res=. res , fassert (norm1t"1@|: -: norm1tc) (gemat j. gemat) 100 100
  res=. res , fassert (norm1t"1    -: norm1tr)  gemat           100 100
  res=. res , fassert (norm1t"1    -: norm1tr) (gemat j. gemat) 100 100
  res=. res , fassert (normit"1@|: -: normitc)  gemat           100 100
  res=. res , fassert (normit"1@|: -: normitc) (gemat j. gemat) 100 100
  res=. res , fassert (normit"1    -: normitr)  gemat           100 100
  res=. res , fassert (normit"1    -: normitr) (gemat j. gemat) 100 100

  NB. square-based norms
  NB. - input is an empty vector
  res=. res , fassert 0  -: norms  ''
  res=. res , fassert '' -: normsc ''
  res=. res , fassert '' -: normsr ''
  NB. - input is a vector containing NaN
  res=. res , fassert isnan norms  _.
  res=. res , fassert isnan normsc _.
  res=. res , fassert isnan normsr _.
  res=. res , fassert isnan norms  _.  2
  res=. res , fassert isnan normsc _.  2
  res=. res , fassert isnan normsr _.  2
  res=. res , fassert isnan norms  _. _.
  res=. res , fassert isnan normsc _. _.
  res=. res , fassert isnan normsr _. _.
  NB. - input is a vector containing ∞
  res=. res , fassert _ -: norms  __
  res=. res , fassert _ -: norms   _
  res=. res , fassert _ -: normsc __
  res=. res , fassert _ -: normsc  _
  res=. res , fassert _ -: normsr __
  res=. res , fassert _ -: normsr  _
  res=. res , fassert _ -: norms  __  2
  res=. res , fassert _ -: norms   _  2
  res=. res , fassert _ -: normsc __  2
  res=. res , fassert _ -: normsc  _  2
  res=. res , fassert _ -: normsr __  2
  res=. res , fassert _ -: normsr  _  2
  res=. res , fassert _ -: norms  __ __
  res=. res , fassert _ -: norms  __  _
  res=. res , fassert _ -: norms   _ __
  res=. res , fassert _ -: norms   _  _
  res=. res , fassert _ -: normsc __ __
  res=. res , fassert _ -: normsc __  _
  res=. res , fassert _ -: normsc  _ __
  res=. res , fassert _ -: normsc  _  _
  res=. res , fassert _ -: normsr __ __
  res=. res , fassert _ -: normsr __  _
  res=. res , fassert _ -: normsr  _ __
  res=. res , fassert _ -: normsr  _  _
  NB. - input is a vector containing big numbers so
  NB.   intermediate calculations may overflow
  res=. res , fassert (2 * %: FP_OVFL) -: norms  _2   * %: FP_OVFL
  res=. res , fassert (2 * %: FP_OVFL) -: norms     2 * %: FP_OVFL
  res=. res , fassert (8 *&%: FP_OVFL) -: norms  _2 2 * %: FP_OVFL
  res=. res , fassert (2 * %: FP_OVFL) -: normsc _2   * %: FP_OVFL
  res=. res , fassert (2 * %: FP_OVFL) -: normsc    2 * %: FP_OVFL
  res=. res , fassert (8 *&%: FP_OVFL) -: normsc _2 2 * %: FP_OVFL
  res=. res , fassert (2 * %: FP_OVFL) -: normsr _2   * %: FP_OVFL
  res=. res , fassert (2 * %: FP_OVFL) -: normsr    2 * %: FP_OVFL
  res=. res , fassert (8 *&%: FP_OVFL) -: normsr _2 2 * %: FP_OVFL
  NB. - input is a vector containing small numbers so
  NB.   intermediate calculations may flush to zero
  res=. res , fassert ( 2        * FP_EPS *&%: FP_UNFL) -: norms  _2   * FP_EPS *&%: FP_UNFL
  res=. res , fassert ( 2        * FP_EPS *&%: FP_UNFL) -: norms     2 * FP_EPS *&%: FP_UNFL
  res=. res , fassert ((2 ^ 3r2) * FP_EPS *&%: FP_UNFL) -: norms  _2 2 * FP_EPS *&%: FP_UNFL
  res=. res , fassert ( 2        * FP_EPS *&%: FP_UNFL) -: normsc _2   * FP_EPS *&%: FP_UNFL
  res=. res , fassert ( 2        * FP_EPS *&%: FP_UNFL) -: normsc    2 * FP_EPS *&%: FP_UNFL
  res=. res , fassert ((2 ^ 3r2) * FP_EPS *&%: FP_UNFL) -: normsc _2 2 * FP_EPS *&%: FP_UNFL
  res=. res , fassert ( 2        * FP_EPS *&%: FP_UNFL) -: normsr _2   * FP_EPS *&%: FP_UNFL
  res=. res , fassert ( 2        * FP_EPS *&%: FP_UNFL) -: normsr    2 * FP_EPS *&%: FP_UNFL
  res=. res , fassert ((2 ^ 3r2) * FP_EPS *&%: FP_UNFL) -: normsr _2 2 * FP_EPS *&%: FP_UNFL
  NB. - all norms are equal for vector input
  res=. res , fassert 1 -: # ~. (norms , normsc , normsr)  gemat           100
  res=. res , fassert 1 -: # ~. (norms , normsc , normsr) (gemat j. gemat) 100
  NB. - norm should be resistable to permutations
  res=. res , fassert 1 -: # ~. norms "1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. norms "1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. normsc"1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. normsc"1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #) (gemat j. gemat) 100
  res=. res , fassert 1 -: # ~. normsr"2 (((C."1~ ?~) , (C."1~ ?~) ,: (C."1~ ?~)) #)  gemat           100
  res=. res , fassert 1 -: # ~. normsr"2 (((C."1~ ?~) , (C."1~ ?~) ,: (C."1~ ?~)) #) (gemat j. gemat) 100
  NB. - input is an empty matrix
  res=. res , fassert 0     -: norms  0 0 $ 0
  res=. res , fassert 0     -: norms  0 3 $ 0
  res=. res , fassert 0     -: norms  3 0 $ 0
  res=. res , fassert ''    -: normsc 0 0 $ 0
  res=. res , fassert 0 0 0 -: normsc 0 3 $ 0
  res=. res , fassert ''    -: normsc 3 0 $ 0
  res=. res , fassert ''    -: normsr 0 0 $ 0
  res=. res , fassert ''    -: normsr 0 3 $ 0
  res=. res , fassert 0 0 0 -: normsr 3 0 $ 0
  NB. - input is a matrix containing NaN
  res=. res , fassert          isnan norms  _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 1 1 1 -: isnan normsc _. (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 1 1 1 -: isnan normsr _. (< 1 1)} 3 3 ?@$ 0
  NB. - input is a matrix containing ∞
  res=. res , fassert         _ = norms  __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert         _ = norms   _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normsc __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normsc  _ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normsr __ (< 1 1)} 3 3 ?@$ 0
  res=. res , fassert 0 1 0 = _ = normsr  _ (< 1 1)} 3 3 ?@$ 0
  NB. - row- and and column-wise norms are related
  res=. res , fassert (norms"1@|: -: normsc)  gemat           100 100
  res=. res , fassert (norms"1@|: -: normsc) (gemat j. gemat) 100 100
  res=. res , fassert (norms"1    -: normsr)  gemat           100 100
  res=. res , fassert (norms"1    -: normsr) (gemat j. gemat) 100 100

  'norm' reportv (# ([ , -) +/) res
)
