NB. Verify norm verbs
NB.
NB. Copyright 2010,2011,2013,2017,2018,2020,2021,2023,2024,
NB.           2025 Igor Zhuravlov
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
NB. Verification suite

NB. magnitude-based norms
NB. - input is an empty vector
0  -: norm1  ''
'' -: norm1c ''
'' -: norm1r ''
0  -: normi  ''
'' -: normic ''
'' -: normir ''
0  -: normm  ''
NB. - input is a vector containing NaN
isnan norm1  _.
isnan norm1c _.
isnan norm1r _.
isnan normi  _.
isnan normic _.
isnan normir _.
isnan normm  _.
isnan norm1  _.  2
isnan norm1c _.  2
isnan norm1r _.  2
isnan normi  _.  2
isnan normic _.  2
isnan normir _.  2
isnan normm  _.  2
isnan norm1  _. _.
isnan norm1c _. _.
isnan norm1r _. _.
isnan normi  _. _.
isnan normic _. _.
isnan normir _. _.
isnan normm  _. _.
NB. - input is a vector containing ∞
_ -: norm1  __
_ -: norm1   _
_ -: norm1c __
_ -: norm1c  _
_ -: norm1r __
_ -: norm1r  _
_ -: normi  __
_ -: normi   _
_ -: normic __
_ -: normic  _
_ -: normir __
_ -: normir  _
_ -: normm  __
_ -: normm   _
_ -: norm1  __  2
_ -: norm1   _  2
_ -: norm1c __  2
_ -: norm1c  _  2
_ -: norm1r __  2
_ -: norm1r  _  2
_ -: normi  __  2
_ -: normi   _  2
_ -: normic __  2
_ -: normic  _  2
_ -: normir __  2
_ -: normir  _  2
_ -: normm  __  2
_ -: normm  __  2
_ -: norm1  __ __
_ -: norm1  __  _
_ -: norm1   _ __
_ -: norm1   _  _
_ -: norm1c __ __
_ -: norm1c __  _
_ -: norm1c  _ __
_ -: norm1c  _  _
_ -: norm1r __ __
_ -: norm1r __  _
_ -: norm1r  _ __
_ -: norm1r  _  _
_ -: normi  __ __
_ -: normi  __  _
_ -: normi   _ __
_ -: normi   _  _
_ -: normic __ __
_ -: normic __  _
_ -: normic  _ __
_ -: normic  _  _
_ -: normir __ __
_ -: normir __  _
_ -: normir  _ __
_ -: normir  _  _
_ -: normm  __ __
_ -: normm  __  _
_ -: normm   _ __
_ -: normm   _  _
NB. - all norms are equal for vector input
1 -: # ~. (norm1 , norm1c , norm1r)  gemat           100
1 -: # ~. (norm1 , norm1c , norm1r) (gemat j. gemat) 100
1 -: # ~. (normi , normic , normir)  gemat           100
1 -: # ~. (normi , normic , normir) (gemat j. gemat) 100
NB. - norm should be resistable to permutations
1 -: # ~. norm1 "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. norm1 "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. norm1c"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. norm1c"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. norm1r"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. norm1r"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. normi "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. normi "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. normic"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. normic"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. normir"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. normir"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
NB. - input is an empty matrix
0     -: norm1  0 0 $ 0
0     -: norm1  0 3 $ 0
0     -: norm1  3 0 $ 0
''    -: norm1c 0 0 $ 0
0 0 0 -: norm1c 0 3 $ 0
''    -: norm1c 3 0 $ 0
''    -: norm1r 0 0 $ 0
''    -: norm1r 0 3 $ 0
0 0 0 -: norm1r 3 0 $ 0
0     -: normi  0 0 $ 0
0     -: normi  0 3 $ 0
0     -: normi  3 0 $ 0
''    -: normic 0 0 $ 0
0 0 0 -: normic 0 3 $ 0
''    -: normic 3 0 $ 0
''    -: normir 0 0 $ 0
''    -: normir 0 3 $ 0
0 0 0 -: normir 3 0 $ 0
0     -: normm  0 0 $ 0
0     -: normm  0 3 $ 0
0     -: normm  3 0 $ 0
NB. - input is a matrix containing NaN
         isnan norm1  _. (< 1 1)} 3 3 ?@$ 0
1 1 1 -: isnan norm1c _. (< 1 1)} 3 3 ?@$ 0
1 1 1 -: isnan norm1r _. (< 1 1)} 3 3 ?@$ 0
         isnan normi  _. (< 1 1)} 3 3 ?@$ 0
1 1 1 -: isnan normic _. (< 1 1)} 3 3 ?@$ 0
1 1 1 -: isnan normir _. (< 1 1)} 3 3 ?@$ 0
         isnan normm  _. (< 1 1)} 3 3 ?@$ 0
NB. - input is a matrix containing ∞
        _ = norm1  __ (< 1 1)} 3 3 ?@$ 0
        _ = norm1   _ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = norm1c __ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = norm1c  _ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = norm1r __ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = norm1r  _ (< 1 1)} 3 3 ?@$ 0
        _ = normi  __ (< 1 1)} 3 3 ?@$ 0
        _ = normi   _ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normic __ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normic  _ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normir __ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normir  _ (< 1 1)} 3 3 ?@$ 0
        _ = normm  __ (< 1 1)} 3 3 ?@$ 0
        _ = normm   _ (< 1 1)} 3 3 ?@$ 0
NB. - row- and and column-wise norms are related
(norm1"1@|: -: norm1c)  gemat           100 100
(norm1"1@|: -: norm1c) (gemat j. gemat) 100 100
(norm1"1    -: norm1r)  gemat           100 100
(norm1"1    -: norm1r) (gemat j. gemat) 100 100
(normi"1@|: -: normic)  gemat           100 100
(normi"1@|: -: normic) (gemat j. gemat) 100 100
(normi"1    -: normir)  gemat           100 100
(normi"1    -: normir) (gemat j. gemat) 100 100

NB. taxicab-based norms
NB. - input is an empty vector
0  -: norm1t  ''
'' -: norm1tc ''
'' -: norm1tr ''
0  -: normit  ''
'' -: normitc ''
'' -: normitr ''
0  -: normmt  ''
NB. - input is a vector containing NaN
isnan norm1t  _.
isnan norm1tc _.
isnan norm1tr _.
isnan normit  _.
isnan normitc _.
isnan normitr _.
isnan normmt  _.
isnan norm1t  _.  2
isnan norm1tc _.  2
isnan norm1tr _.  2
isnan normit  _.  2
isnan normitc _.  2
isnan normitr _.  2
isnan normmt  _.  2
isnan norm1t  _. _.
isnan norm1tc _. _.
isnan norm1tr _. _.
isnan normit  _. _.
isnan normitc _. _.
isnan normitr _. _.
isnan normmt  _. _.
NB. - input is a vector containing ∞
_ -: norm1t  __
_ -: norm1t   _
_ -: norm1tc __
_ -: norm1tc  _
_ -: norm1tr __
_ -: norm1tr  _
_ -: normit  __
_ -: normit   _
_ -: normitc __
_ -: normitc  _
_ -: normitr __
_ -: normitr  _
_ -: normmt  __
_ -: normmt   _
_ -: norm1t  __  2
_ -: norm1t   _  2
_ -: norm1tc __  2
_ -: norm1tc  _  2
_ -: norm1tr __  2
_ -: norm1tr  _  2
_ -: normit  __  2
_ -: normit   _  2
_ -: normitc __  2
_ -: normitc  _  2
_ -: normitr __  2
_ -: normitr  _  2
_ -: normmt  __  2
_ -: normmt  __  2
_ -: norm1t  __ __
_ -: norm1t  __  _
_ -: norm1t   _ __
_ -: norm1t   _  _
_ -: norm1tc __ __
_ -: norm1tc __  _
_ -: norm1tc  _ __
_ -: norm1tc  _  _
_ -: norm1tr __ __
_ -: norm1tr __  _
_ -: norm1tr  _ __
_ -: norm1tr  _  _
_ -: normit  __ __
_ -: normit  __  _
_ -: normit   _ __
_ -: normit   _  _
_ -: normitc __ __
_ -: normitc __  _
_ -: normitc  _ __
_ -: normitc  _  _
_ -: normitr __ __
_ -: normitr __  _
_ -: normitr  _ __
_ -: normitr  _  _
_ -: normmt  __ __
_ -: normmt  __  _
_ -: normmt   _ __
_ -: normmt   _  _
NB. - all norms are equal for vector input
1 -: # ~. (norm1t , norm1tc , norm1tr)  gemat           100
1 -: # ~. (norm1t , norm1tc , norm1tr) (gemat j. gemat) 100
1 -: # ~. (normit , normitc , normitr)  gemat           100
1 -: # ~. (normit , normitc , normitr) (gemat j. gemat) 100
NB. - norm should be resistable to permutations
1 -: # ~. norm1t "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. norm1t "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. norm1tc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. norm1tc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. norm1tr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. norm1tr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. normit "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. normit "1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. normitc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. normitc"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. normitr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #)  gemat           100
1 -: # ~. normitr"1 (((C.~ ?~) , (C.~ ?~) ,: (C.~ ?~)) #) (gemat j. gemat) 100
NB. - input is an empty matrix
0     -: norm1t  0 0 $ 0
0     -: norm1t  0 3 $ 0
0     -: norm1t  3 0 $ 0
''    -: norm1tc 0 0 $ 0
0 0 0 -: norm1tc 0 3 $ 0
''    -: norm1tc 3 0 $ 0
''    -: norm1tr 0 0 $ 0
''    -: norm1tr 0 3 $ 0
0 0 0 -: norm1tr 3 0 $ 0
0     -: normit  0 0 $ 0
0     -: normit  0 3 $ 0
0     -: normit  3 0 $ 0
''    -: normitc 0 0 $ 0
0 0 0 -: normitc 0 3 $ 0
''    -: normitc 3 0 $ 0
''    -: normitr 0 0 $ 0
''    -: normitr 0 3 $ 0
0 0 0 -: normitr 3 0 $ 0
0     -: normmt  0 0 $ 0
0     -: normmt  0 3 $ 0
0     -: normmt  3 0 $ 0
NB. - input is a matrix containing NaN
         isnan norm1t  _. (< 1 1)} 3 3 ?@$ 0
1 1 1 -: isnan norm1tc _. (< 1 1)} 3 3 ?@$ 0
1 1 1 -: isnan norm1tr _. (< 1 1)} 3 3 ?@$ 0
         isnan normit  _. (< 1 1)} 3 3 ?@$ 0
1 1 1 -: isnan normitc _. (< 1 1)} 3 3 ?@$ 0
1 1 1 -: isnan normitr _. (< 1 1)} 3 3 ?@$ 0
         isnan normmt  _. (< 1 1)} 3 3 ?@$ 0
NB. - input is a matrix containing ∞
        _ = norm1t  __ (< 1 1)} 3 3 ?@$ 0
        _ = norm1t   _ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = norm1tc __ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = norm1tc  _ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = norm1tr __ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = norm1tr  _ (< 1 1)} 3 3 ?@$ 0
        _ = normit  __ (< 1 1)} 3 3 ?@$ 0
        _ = normit   _ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normitc __ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normitc  _ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normitr __ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normitr  _ (< 1 1)} 3 3 ?@$ 0
        _ = normmt  __ (< 1 1)} 3 3 ?@$ 0
        _ = normmt   _ (< 1 1)} 3 3 ?@$ 0
NB. - row- and and column-wise norms are related
(norm1t"1@|: -: norm1tc)  gemat           100 100
(norm1t"1@|: -: norm1tc) (gemat j. gemat) 100 100
(norm1t"1    -: norm1tr)  gemat           100 100
(norm1t"1    -: norm1tr) (gemat j. gemat) 100 100
(normit"1@|: -: normitc)  gemat           100 100
(normit"1@|: -: normitc) (gemat j. gemat) 100 100
(normit"1    -: normitr)  gemat           100 100
(normit"1    -: normitr) (gemat j. gemat) 100 100

NB. square-based norms
NB. - input is an empty vector
0  -: norms  ''
'' -: normsc ''
'' -: normsr ''
NB. - input is a vector containing NaN
isnan norms  _.
isnan normsc _.
isnan normsr _.
isnan norms  _.  2
isnan normsc _.  2
isnan normsr _.  2
isnan norms  _. _.
isnan normsc _. _.
isnan normsr _. _.
NB. - input is a vector containing ∞
_ -: norms  __
_ -: norms   _
_ -: normsc __
_ -: normsc  _
_ -: normsr __
_ -: normsr  _
_ -: norms  __  2
_ -: norms   _  2
_ -: normsc __  2
_ -: normsc  _  2
_ -: normsr __  2
_ -: normsr  _  2
_ -: norms  __ __
_ -: norms  __  _
_ -: norms   _ __
_ -: norms   _  _
_ -: normsc __ __
_ -: normsc __  _
_ -: normsc  _ __
_ -: normsc  _  _
_ -: normsr __ __
_ -: normsr __  _
_ -: normsr  _ __
_ -: normsr  _  _
NB. - input is a vector containing big numbers so
NB.   intermediate calculations may overflow
(2 * %: FP_OVFL) -: norms  _2   * %: FP_OVFL
(2 * %: FP_OVFL) -: norms     2 * %: FP_OVFL
(8 *&%: FP_OVFL) -: norms  _2 2 * %: FP_OVFL
(2 * %: FP_OVFL) -: normsc _2   * %: FP_OVFL
(2 * %: FP_OVFL) -: normsc    2 * %: FP_OVFL
(8 *&%: FP_OVFL) -: normsc _2 2 * %: FP_OVFL
(2 * %: FP_OVFL) -: normsr _2   * %: FP_OVFL
(2 * %: FP_OVFL) -: normsr    2 * %: FP_OVFL
(8 *&%: FP_OVFL) -: normsr _2 2 * %: FP_OVFL
NB. - input is a vector containing small numbers so
NB.   intermediate calculations may flush to zero
( 2        * FP_EPS *&%: FP_UNFL) -: norms  _2   * FP_EPS *&%: FP_UNFL
( 2        * FP_EPS *&%: FP_UNFL) -: norms     2 * FP_EPS *&%: FP_UNFL
((2 ^ 3r2) * FP_EPS *&%: FP_UNFL) -: norms  _2 2 * FP_EPS *&%: FP_UNFL
( 2        * FP_EPS *&%: FP_UNFL) -: normsc _2   * FP_EPS *&%: FP_UNFL
( 2        * FP_EPS *&%: FP_UNFL) -: normsc    2 * FP_EPS *&%: FP_UNFL
((2 ^ 3r2) * FP_EPS *&%: FP_UNFL) -: normsc _2 2 * FP_EPS *&%: FP_UNFL
( 2        * FP_EPS *&%: FP_UNFL) -: normsr _2   * FP_EPS *&%: FP_UNFL
( 2        * FP_EPS *&%: FP_UNFL) -: normsr    2 * FP_EPS *&%: FP_UNFL
((2 ^ 3r2) * FP_EPS *&%: FP_UNFL) -: normsr _2 2 * FP_EPS *&%: FP_UNFL
NB. - all norms are equal for vector input
1 -: # ~. (norms , normsc , normsr)  gemat           100
1 -: # ~. (norms , normsc , normsr) (gemat j. gemat) 100
NB. - norm should be resistable to permutations
1 -: # ~. norms "1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #)  gemat           100
1 -: # ~. norms "1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. normsc"1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #)  gemat           100
1 -: # ~. normsc"1 (((C.  ~ ?~) , (C.  ~ ?~) ,: (C.  ~ ?~)) #) (gemat j. gemat) 100
1 -: # ~. normsr"2 (((C."1~ ?~) , (C."1~ ?~) ,: (C."1~ ?~)) #)  gemat           100
1 -: # ~. normsr"2 (((C."1~ ?~) , (C."1~ ?~) ,: (C."1~ ?~)) #) (gemat j. gemat) 100
NB. - input is an empty matrix
0     -: norms  0 0 $ 0
0     -: norms  0 3 $ 0
0     -: norms  3 0 $ 0
''    -: normsc 0 0 $ 0
0 0 0 -: normsc 0 3 $ 0
''    -: normsc 3 0 $ 0
''    -: normsr 0 0 $ 0
''    -: normsr 0 3 $ 0
0 0 0 -: normsr 3 0 $ 0
NB. - input is a matrix containing NaN
         isnan norms  _. (< 1 1)} 3 3 ?@$ 0
1 1 1 -: isnan normsc _. (< 1 1)} 3 3 ?@$ 0
1 1 1 -: isnan normsr _. (< 1 1)} 3 3 ?@$ 0
NB. - input is a matrix containing ∞
        _ = norms  __ (< 1 1)} 3 3 ?@$ 0
        _ = norms   _ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normsc __ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normsc  _ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normsr __ (< 1 1)} 3 3 ?@$ 0
0 1 0 = _ = normsr  _ (< 1 1)} 3 3 ?@$ 0
NB. - row- and and column-wise norms are related
(norms"1@|: -: normsc)  gemat           100 100
(norms"1@|: -: normsc) (gemat j. gemat) 100 100
(norms"1    -: normsr)  gemat           100 100
(norms"1    -: normsr) (gemat j. gemat) 100 100
