NB. Verify struct actors
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
NB. Local definitions

a005=. 0 5 $ 0
a034=. 3 4 $ 0
a045=. 4 5 $ 0
a050=. 5 0 $ 0
a054=. 5 4 $ 0
a055=. 5 5 $ 0
a145=. 4 5 $ 1
a154=. 5 4 $ 1
a155=. 5 5 $ 1
a245=. 4 5 $ 2
ai23=. i. 2 3
ai34=. i. 3 4
ai35=. i. 3 5
ai44=. i. 4 4
ai45=. i. 4 5
ai53=. i. 5 3
ai54=. i. 5 4
ai55=. i. 5 5
ac23=. j./ i. 2 2 3
ac35=. j./ i. 2 3 5
ac44=. j./ i. 2 4 4
ac55=. j./ i. 2 5 5
p=.    2 1 3 0
P=.     P4p p
iP=.   P4ip p
fret=. 1 0 0 1 1 0 1 0 0 1

NB. =========================================================
NB. Verification suite

NB. ft4lisoa
(-: < ft4lisoa) EMPTY
(3 4 $ 0  1  1  1  0 0  1  1  0  0 0  1) -: <  ft4lisoa a034
(3 4 $ 0  0  0  0  1 0  0  0  1  1 0  0) -: >  ft4lisoa a034
(3 4 $ 1  1  1  1  0 1  1  1  0  0 1  1) -: <: ft4lisoa a034
(3 4 $ 1  0  0  0  1 1  0  0  1  1 1  0) -: >: ft4lisoa a034
(3 4 $ 0 _1 _2 _3  1 0 _1 _2  2  1 0 _1) -: -  ft4lisoa a034
(3 4 $ 0  1  2  3 _1 0  1  2 _2 _1 0  1) -: -~ ft4lisoa a034
(3 4 $ 0  1  2  3  1 2  3  4  2  3 4  5) -: +  ft4lisoa a034

NB. mbstencil
(-: 0   &mbstencil) EMPTY
(-: 0 0 &mbstencil) EMPTY
(-: __ _&mbstencil) EMPTY
(3 5 $ 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0) -: 1              mbstencil ai35
(3 5 $ 0 0 1 1 0 0 0 0 1 1 0 0 0 0 1) -: 2 3            mbstencil ai35
(3 5 $ 0 0 1 1 0 1 0 0 1 1 1 1 0 0 1) -: (__ _1 ,: 2 3) mbstencil ai35

NB. mabstencil
(-: 0   &mabstencil) EMPTY
(-: 0 0 &mabstencil) EMPTY
(-: __ _&mabstencil) EMPTY
(3 5 $ 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0) -: 1              mabstencil ai35
(3 5 $ 0 1 1 0 0 1 1 0 0 0 1 0 0 0 0) -: 2 3            mabstencil ai35
(3 5 $ 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1) -: (__ _1 ,: 2 3) mabstencil ai35

NB. diagliso

'' -:         diagliso 0
'' -: 0       diagliso 0
'' -: 0  0    diagliso 0
'' -: 0  0  0 diagliso 0
'' -: 0  0  _ diagliso 0
'' -: 0  0 __ diagliso 0
'' -: 0 _1    diagliso 0
'' -: 0 _1  0 diagliso 0
'' -: 0 _1  _ diagliso 0
'' -: 0 _1 __ diagliso 0

'' -:         diagliso 0 0
'' -: 0       diagliso 0 0
'' -: 0  0    diagliso 0 0
'' -: 0  0  0 diagliso 0 0
'' -: 0  0  _ diagliso 0 0
'' -: 0  0 __ diagliso 0 0
'' -: 0 _1    diagliso 0 0
'' -: 0 _1  0 diagliso 0 0
'' -: 0 _1  _ diagliso 0 0
'' -: 0 _1 __ diagliso 0 0

'' -:         diagliso 0 5
'' -: 0       diagliso 0 5
'' -: 0  0    diagliso 0 5
'' -: 0  0  0 diagliso 0 5
'' -: 0  0  _ diagliso 0 5
'' -: 0  0 __ diagliso 0 5
'' -: 0 _1    diagliso 0 5
'' -: 0 _1  0 diagliso 0 5
'' -: 0 _1  _ diagliso 0 5
'' -: 0 _1 __ diagliso 0 5

'' -:         diagliso 5 0
'' -: 0       diagliso 5 0
'' -: 0  0    diagliso 5 0
'' -: 0  0  0 diagliso 5 0
'' -: 0  0  _ diagliso 5 0
'' -: 0  0 __ diagliso 5 0
'' -: 0 _1    diagliso 5 0
'' -: 0 _1  0 diagliso 5 0
'' -: 0 _1  _ diagliso 5 0
'' -: 0 _1 __ diagliso 5 0

0 6 12 18 -:          diagliso 4 5
0 6 12 18 -:  0       diagliso 4 5
0 6 12 18 -:  0  0    diagliso 4 5
''        -:  0  0  0 diagliso 4 5
0 6       -:  0  0  2 diagliso 4 5
0 6 12 18 -:  0  0  _ diagliso 4 5
6 0       -:  0  0 _2 diagliso 4 5
18 12 6 0 -:  0  0 __ diagliso 4 5
0 6 12 18 -:  0 _1    diagliso 4 5
''        -:  0 _1  0 diagliso 4 5
    12 18 -:  0 _1  2 diagliso 4 5
0 6 12 18 -:  0 _1  _ diagliso 4 5
    18 12 -:  0 _1 _2 diagliso 4 5
18 12 6 0 -:  0 _1 __ diagliso 4 5
  6 12 18 -:  0  1    diagliso 4 5
''        -:  0  1  0 diagliso 4 5
  6 12    -:  0  1  2 diagliso 4 5
  6 12 18 -:  0  1  _ diagliso 4 5
  12 6    -:  0  1 _2 diagliso 4 5
  18 12 6 -:  0  1 __ diagliso 4 5
0 6 12    -:  0 _2    diagliso 4 5
''        -:  0 _2  0 diagliso 4 5
  6 12    -:  0 _2  2 diagliso 4 5
0 6 12    -:  0 _2  _ diagliso 4 5
  12 6    -:  0 _2 _2 diagliso 4 5
12 6 0    -:  0 _2 __ diagliso 4 5
5 11 17   -: _1       diagliso 4 5
5 11 17   -: _1  0    diagliso 4 5
''        -: _1  0  0 diagliso 4 5
5 11      -: _1  0  2 diagliso 4 5
5 11 17   -: _1  0  _ diagliso 4 5
11 5      -: _1  0 _2 diagliso 4 5
17 11 5   -: _1  0 __ diagliso 4 5
5 11 17   -: _1 _1    diagliso 4 5
''        -: _1 _1  0 diagliso 4 5
  11 17   -: _1 _1  2 diagliso 4 5
5 11 17   -: _1 _1  _ diagliso 4 5
  17 11   -: _1 _1 _2 diagliso 4 5
17 11 5   -: _1 _1 __ diagliso 4 5
  11 17   -: _1  1    diagliso 4 5
''        -: _1  1  0 diagliso 4 5
  11 17   -: _1  1  2 diagliso 4 5
  11 17   -: _1  1  _ diagliso 4 5
  17 11   -: _1  1 _2 diagliso 4 5
  17 11   -: _1  1 __ diagliso 4 5
5 11      -: _1 _2    diagliso 4 5
''        -: _1 _2  0 diagliso 4 5
5 11      -: _1 _2  2 diagliso 4 5
5 11      -: _1 _2  _ diagliso 4 5
11 5      -: _1 _2 _2 diagliso 4 5
11 5      -: _1 _2 __ diagliso 4 5
1 7 13 19 -:  1       diagliso 4 5
1 7 13 19 -:  1  0    diagliso 4 5
''        -:  1  0  0 diagliso 4 5
1 7       -:  1  0  2 diagliso 4 5
1 7 13 19 -:  1  0  _ diagliso 4 5
7 1       -:  1  0 _2 diagliso 4 5
19 13 7 1 -:  1  0 __ diagliso 4 5
1 7 13 19 -:  1 _1    diagliso 4 5
''        -:  1 _1  0 diagliso 4 5
    13 19 -:  1 _1  2 diagliso 4 5
1 7 13 19 -:  1 _1  _ diagliso 4 5
    19 13 -:  1 _1 _2 diagliso 4 5
19 13 7 1 -:  1 _1 __ diagliso 4 5
  7 13 19 -:  1  1    diagliso 4 5
''        -:  1  1  0 diagliso 4 5
  7 13    -:  1  1  2 diagliso 4 5
  7 13 19 -:  1  1  _ diagliso 4 5
  13 7    -:  1  1 _2 diagliso 4 5
  19 13 7 -:  1  1 __ diagliso 4 5
1 7 13    -:  1 _2    diagliso 4 5
''        -:  1 _2  0 diagliso 4 5
  7 13    -:  1 _2  2 diagliso 4 5
1 7 13    -:  1 _2  _ diagliso 4 5
  13 7    -:  1 _2 _2 diagliso 4 5
13 7 1    -:  1 _2 __ diagliso 4 5
2 8 14    -:  2       diagliso 4 5
2 8 14    -:  2  0    diagliso 4 5
''        -:  2  0  0 diagliso 4 5
2 8       -:  2  0  2 diagliso 4 5
2 8 14    -:  2  0  _ diagliso 4 5
8 2       -:  2  0 _2 diagliso 4 5
14 8 2    -:  2  0 __ diagliso 4 5
2 8 14    -:  2 _1    diagliso 4 5
''        -:  2 _1  0 diagliso 4 5
  8 14    -:  2 _1  2 diagliso 4 5
2 8 14    -:  2 _1  _ diagliso 4 5
  14 8    -:  2 _1 _2 diagliso 4 5
14 8 2    -:  2 _1 __ diagliso 4 5
  8 14    -:  2  1    diagliso 4 5
''        -:  2  1  0 diagliso 4 5
  8 14    -:  2  1  2 diagliso 4 5
  8 14    -:  2  1  _ diagliso 4 5
  14 8    -:  2  1 _2 diagliso 4 5
  14 8    -:  2  1 __ diagliso 4 5
2 8       -:  2 _2    diagliso 4 5
''        -:  2 _2  0 diagliso 4 5
2 8       -:  2 _2  2 diagliso 4 5
2 8       -:  2 _2  _ diagliso 4 5
8 2       -:  2 _2 _2 diagliso 4 5
8 2       -:  2 _2 __ diagliso 4 5
10 16     -: _2       diagliso 4 5
10 16     -: _2  0    diagliso 4 5
''        -: _2  0  0 diagliso 4 5
10 16     -: _2  0  2 diagliso 4 5
10 16     -: _2  0  _ diagliso 4 5
16 10     -: _2  0 _2 diagliso 4 5
16 10     -: _2  0 __ diagliso 4 5
10 16     -: _2 _1    diagliso 4 5
''        -: _2 _1  0 diagliso 4 5
10 16     -: _2 _1  2 diagliso 4 5
10 16     -: _2 _1  _ diagliso 4 5
16 10     -: _2 _1 _2 diagliso 4 5
16 10     -: _2 _1 __ diagliso 4 5
(, 16)    -: _2  1    diagliso 4 5
''        -: _2  1  0 diagliso 4 5
(, 16)    -: _2  1  2 diagliso 4 5
(, 16)    -: _2  1  _ diagliso 4 5
(, 16)    -: _2  1 _2 diagliso 4 5
(, 16)    -: _2  1 __ diagliso 4 5
(, 10)    -: _2 _2    diagliso 4 5
''        -: _2 _2  0 diagliso 4 5
(, 10)    -: _2 _2  2 diagliso 4 5
(, 10)    -: _2 _2  _ diagliso 4 5
(, 10)    -: _2 _2 _2 diagliso 4 5
(, 10)    -: _2 _2 __ diagliso 4 5

0 5 10 15 -:          diagliso 5 4
0 5 10 15 -:  0       diagliso 5 4
0 5 10 15 -:  0  0    diagliso 5 4
''        -:  0  0  0 diagliso 5 4
0 5       -:  0  0  2 diagliso 5 4
0 5 10 15 -:  0  0  _ diagliso 5 4
5 0       -:  0  0 _2 diagliso 5 4
15 10 5 0 -:  0  0 __ diagliso 5 4
0 5 10 15 -:  0 _1    diagliso 5 4
''        -:  0 _1  0 diagliso 5 4
    10 15 -:  0 _1  2 diagliso 5 4
0 5 10 15 -:  0 _1  _ diagliso 5 4
    15 10 -:  0 _1 _2 diagliso 5 4
15 10 5 0 -:  0 _1 __ diagliso 5 4
  5 10 15 -:  0  1    diagliso 5 4
''        -:  0  1  0 diagliso 5 4
  5 10    -:  0  1  2 diagliso 5 4
  5 10 15 -:  0  1  _ diagliso 5 4
  10 5    -:  0  1 _2 diagliso 5 4
  15 10 5 -:  0  1 __ diagliso 5 4
0 5 10    -:  0 _2    diagliso 5 4
''        -:  0 _2  0 diagliso 5 4
  5 10    -:  0 _2  2 diagliso 5 4
0 5 10    -:  0 _2  _ diagliso 5 4
  10 5    -:  0 _2 _2 diagliso 5 4
10 5 0    -:  0 _2 __ diagliso 5 4
4 9 14 19 -: _1       diagliso 5 4
4 9 14 19 -: _1  0    diagliso 5 4
''        -: _1  0  0 diagliso 5 4
4 9       -: _1  0  2 diagliso 5 4
4 9 14 19 -: _1  0  _ diagliso 5 4
9 4       -: _1  0 _2 diagliso 5 4
19 14 9 4 -: _1  0 __ diagliso 5 4
4 9 14 19 -: _1 _1    diagliso 5 4
''        -: _1 _1  0 diagliso 5 4
    14 19 -: _1 _1  2 diagliso 5 4
4 9 14 19 -: _1 _1  _ diagliso 5 4
    19 14 -: _1 _1 _2 diagliso 5 4
19 14 9 4 -: _1 _1 __ diagliso 5 4
  9 14 19 -: _1  1    diagliso 5 4
''        -: _1  1  0 diagliso 5 4
  9 14    -: _1  1  2 diagliso 5 4
  9 14 19 -: _1  1  _ diagliso 5 4
  14 9    -: _1  1 _2 diagliso 5 4
  19 14 9 -: _1  1 __ diagliso 5 4
4 9 14    -: _1 _2    diagliso 5 4
''        -: _1 _2  0 diagliso 5 4
  9 14    -: _1 _2  2 diagliso 5 4
4 9 14    -: _1 _2  _ diagliso 5 4
  14 9    -: _1 _2 _2 diagliso 5 4
14 9 4    -: _1 _2 __ diagliso 5 4
1 6 11    -:  1       diagliso 5 4
1 6 11    -:  1  0    diagliso 5 4
''        -:  1  0  0 diagliso 5 4
1 6       -:  1  0  2 diagliso 5 4
1 6 11    -:  1  0  _ diagliso 5 4
6 1       -:  1  0 _2 diagliso 5 4
11 6 1    -:  1  0 __ diagliso 5 4
1 6 11    -:  1 _1    diagliso 5 4
''        -:  1 _1  0 diagliso 5 4
  6 11    -:  1 _1  2 diagliso 5 4
1 6 11    -:  1 _1  _ diagliso 5 4
  11 6    -:  1 _1 _2 diagliso 5 4
11 6 1    -:  1 _1 __ diagliso 5 4
  6 11    -:  1  1    diagliso 5 4
''        -:  1  1  0 diagliso 5 4
  6 11    -:  1  1  2 diagliso 5 4
  6 11    -:  1  1  _ diagliso 5 4
  11 6    -:  1  1 _2 diagliso 5 4
  11 6    -:  1  1 __ diagliso 5 4
1 6       -:  1 _2    diagliso 5 4
''        -:  1 _2  0 diagliso 5 4
1 6       -:  1 _2  2 diagliso 5 4
1 6       -:  1 _2  _ diagliso 5 4
6 1       -:  1 _2 _2 diagliso 5 4
6 1       -:  1 _2 __ diagliso 5 4
2 7       -:  2       diagliso 5 4
2 7       -:  2  0    diagliso 5 4
''        -:  2  0  0 diagliso 5 4
2 7       -:  2  0  2 diagliso 5 4
2 7       -:  2  0  _ diagliso 5 4
7 2       -:  2  0 _2 diagliso 5 4
7 2       -:  2  0 __ diagliso 5 4
2 7       -:  2 _1    diagliso 5 4
''        -:  2 _1  0 diagliso 5 4
2 7       -:  2 _1  2 diagliso 5 4
2 7       -:  2 _1  _ diagliso 5 4
7 2       -:  2 _1 _2 diagliso 5 4
7 2       -:  2 _1 __ diagliso 5 4
  (, 7)   -:  2  1    diagliso 5 4
''        -:  2  1  0 diagliso 5 4
  (, 7)   -:  2  1  2 diagliso 5 4
  (, 7)   -:  2  1  _ diagliso 5 4
  (, 7)   -:  2  1 _2 diagliso 5 4
  (, 7)   -:  2  1 __ diagliso 5 4
(, 2)     -:  2 _2    diagliso 5 4
''        -:  2 _2  0 diagliso 5 4
(, 2)     -:  2 _2  2 diagliso 5 4
(, 2)     -:  2 _2  _ diagliso 5 4
(, 2)     -:  2 _2 _2 diagliso 5 4
(, 2)     -:  2 _2 __ diagliso 5 4
8 13 18   -: _2       diagliso 5 4
8 13 18   -: _2  0    diagliso 5 4
''        -: _2  0  0 diagliso 5 4
8 13      -: _2  0  2 diagliso 5 4
8 13 18   -: _2  0  _ diagliso 5 4
13 8      -: _2  0 _2 diagliso 5 4
18 13 8   -: _2  0 __ diagliso 5 4
8 13 18   -: _2 _1    diagliso 5 4
''        -: _2 _1  0 diagliso 5 4
  13 18   -: _2 _1  2 diagliso 5 4
8 13 18   -: _2 _1  _ diagliso 5 4
  18 13   -: _2 _1 _2 diagliso 5 4
18 13 8   -: _2 _1 __ diagliso 5 4
  13 18   -: _2  1    diagliso 5 4
''        -: _2  1  0 diagliso 5 4
  13 18   -: _2  1  2 diagliso 5 4
  13 18   -: _2  1  _ diagliso 5 4
  18 13   -: _2  1 _2 diagliso 5 4
  18 13   -: _2  1 __ diagliso 5 4
8 13      -: _2 _2    diagliso 5 4
''        -: _2 _2  0 diagliso 5 4
8 13      -: _2 _2  2 diagliso 5 4
8 13      -: _2 _2  _ diagliso 5 4
13 8      -: _2 _2 _2 diagliso 5 4
13 8      -: _2 _2 __ diagliso 5 4

0 6 12 18 24 -:          diagliso 5 5
0 6 12 18 24 -:  0       diagliso 5 5
0 6 12 18 24 -:  0  0    diagliso 5 5
''           -:  0  0  0 diagliso 5 5
0 6          -:  0  0  2 diagliso 5 5
0 6 12 18 24 -:  0  0  _ diagliso 5 5
6 0          -:  0  0 _2 diagliso 5 5
24 18 12 6 0 -:  0  0 __ diagliso 5 5
0 6 12 18 24 -:  0 _1    diagliso 5 5
''           -:  0 _1  0 diagliso 5 5
       18 24 -:  0 _1  2 diagliso 5 5
0 6 12 18 24 -:  0 _1  _ diagliso 5 5
       24 18 -:  0 _1 _2 diagliso 5 5
24 18 12 6 0 -:  0 _1 __ diagliso 5 5
  6 12 18 24 -:  0  1    diagliso 5 5
''           -:  0  1  0 diagliso 5 5
  6 12       -:  0  1  2 diagliso 5 5
  6 12 18 24 -:  0  1  _ diagliso 5 5
  12 6       -:  0  1 _2 diagliso 5 5
  24 18 12 6 -:  0  1 __ diagliso 5 5
0 6 12 18    -:  0 _2    diagliso 5 5
''           -:  0 _2  0 diagliso 5 5
    12 18    -:  0 _2  2 diagliso 5 5
0 6 12 18    -:  0 _2  _ diagliso 5 5
    18 12    -:  0 _2 _2 diagliso 5 5
18 12 6 0    -:  0 _2 __ diagliso 5 5
5 11 17 23   -: _1       diagliso 5 5
5 11 17 23   -: _1  0    diagliso 5 5
''           -: _1  0  0 diagliso 5 5
5 11         -: _1  0  2 diagliso 5 5
5 11 17 23   -: _1  0  _ diagliso 5 5
11 5         -: _1  0 _2 diagliso 5 5
23 17 11 5   -: _1  0 __ diagliso 5 5
5 11 17 23   -: _1 _1    diagliso 5 5
''           -: _1 _1  0 diagliso 5 5
     17 23   -: _1 _1  2 diagliso 5 5
5 11 17 23   -: _1 _1  _ diagliso 5 5
     23 17   -: _1 _1 _2 diagliso 5 5
23 17 11 5   -: _1 _1 __ diagliso 5 5
  11 17 23   -: _1  1    diagliso 5 5
''           -: _1  1  0 diagliso 5 5
  11 17      -: _1  1  2 diagliso 5 5
  11 17 23   -: _1  1  _ diagliso 5 5
  17 11      -: _1  1 _2 diagliso 5 5
  23 17 11   -: _1  1 __ diagliso 5 5
5 11 17      -: _1 _2    diagliso 5 5
''           -: _1 _2  0 diagliso 5 5
  11 17      -: _1 _2  2 diagliso 5 5
5 11 17      -: _1 _2  _ diagliso 5 5
  17 11      -: _1 _2 _2 diagliso 5 5
17 11 5      -: _1 _2 __ diagliso 5 5
1 7 13 19    -:  1       diagliso 5 5
1 7 13 19    -:  1  0    diagliso 5 5
''           -:  1  0  0 diagliso 5 5
1 7          -:  1  0  2 diagliso 5 5
1 7 13 19    -:  1  0  _ diagliso 5 5
7 1          -:  1  0 _2 diagliso 5 5
19 13 7 1    -:  1  0 __ diagliso 5 5
1 7 13 19    -:  1 _1    diagliso 5 5
''           -:  1 _1  0 diagliso 5 5
    13 19    -:  1 _1  2 diagliso 5 5
1 7 13 19    -:  1 _1  _ diagliso 5 5
    19 13    -:  1 _1 _2 diagliso 5 5
19 13 7 1    -:  1 _1 __ diagliso 5 5
  7 13 19    -:  1  1    diagliso 5 5
''           -:  1  1  0 diagliso 5 5
  7 13       -:  1  1  2 diagliso 5 5
  7 13 19    -:  1  1  _ diagliso 5 5
  13 7       -:  1  1 _2 diagliso 5 5
  19 13 7    -:  1  1 __ diagliso 5 5
1 7 13       -:  1 _2    diagliso 5 5
''           -:  1 _2  0 diagliso 5 5
  7 13       -:  1 _2  2 diagliso 5 5
1 7 13       -:  1 _2  _ diagliso 5 5
  13 7       -:  1 _2 _2 diagliso 5 5
13 7 1       -:  1 _2 __ diagliso 5 5
2 8 14       -:  2       diagliso 5 5
2 8 14       -:  2  0    diagliso 5 5
''           -:  2  0  0 diagliso 5 5
2 8          -:  2  0  2 diagliso 5 5
2 8 14       -:  2  0  _ diagliso 5 5
8 2          -:  2  0 _2 diagliso 5 5
14 8 2       -:  2  0 __ diagliso 5 5
2 8 14       -:  2 _1    diagliso 5 5
''           -:  2 _1  0 diagliso 5 5
  8 14       -:  2 _1  2 diagliso 5 5
2 8 14       -:  2 _1  _ diagliso 5 5
  14 8       -:  2 _1 _2 diagliso 5 5
14 8 2       -:  2 _1 __ diagliso 5 5
  8 14       -:  2  1    diagliso 5 5
''           -:  2  1  0 diagliso 5 5
  8 14       -:  2  1  2 diagliso 5 5
  8 14       -:  2  1  _ diagliso 5 5
  14 8       -:  2  1 _2 diagliso 5 5
  14 8       -:  2  1 __ diagliso 5 5
2 8          -:  2 _2    diagliso 5 5
''           -:  2 _2  0 diagliso 5 5
2 8          -:  2 _2  2 diagliso 5 5
2 8          -:  2 _2  _ diagliso 5 5
8 2          -:  2 _2 _2 diagliso 5 5
8 2          -:  2 _2 __ diagliso 5 5
10 16 22     -: _2       diagliso 5 5
10 16 22     -: _2  0    diagliso 5 5
''           -: _2  0  0 diagliso 5 5
10 16        -: _2  0  2 diagliso 5 5
10 16 22     -: _2  0  _ diagliso 5 5
16 10        -: _2  0 _2 diagliso 5 5
22 16 10     -: _2  0 __ diagliso 5 5
10 16 22     -: _2 _1    diagliso 5 5
''           -: _2 _1  0 diagliso 5 5
16 22        -: _2 _1  2 diagliso 5 5
10 16 22     -: _2 _1  _ diagliso 5 5
22 16        -: _2 _1 _2 diagliso 5 5
22 16 10     -: _2 _1 __ diagliso 5 5
16 22        -: _2  1    diagliso 5 5
''           -: _2  1  0 diagliso 5 5
16 22        -: _2  1  2 diagliso 5 5
16 22        -: _2  1  _ diagliso 5 5
22 16        -: _2  1 _2 diagliso 5 5
22 16        -: _2  1 __ diagliso 5 5
10 16        -: _2 _2    diagliso 5 5
''           -: _2 _2  0 diagliso 5 5
10 16        -: _2 _2  2 diagliso 5 5
10 16        -: _2 _2  _ diagliso 5 5
16 10        -: _2 _2 _2 diagliso 5 5
16 10        -: _2 _2 __ diagliso 5 5

0 6 12 18 24 -:          diagliso 5
0 6 12 18 24 -:  0       diagliso 5
0 6 12 18 24 -:  0  0    diagliso 5
''           -:  0  0  0 diagliso 5
0 6          -:  0  0  2 diagliso 5
0 6 12 18 24 -:  0  0  _ diagliso 5
6 0          -:  0  0 _2 diagliso 5
24 18 12 6 0 -:  0  0 __ diagliso 5
0 6 12 18 24 -:  0 _1    diagliso 5
''           -:  0 _1  0 diagliso 5
       18 24 -:  0 _1  2 diagliso 5
0 6 12 18 24 -:  0 _1  _ diagliso 5
       24 18 -:  0 _1 _2 diagliso 5
24 18 12 6 0 -:  0 _1 __ diagliso 5
  6 12 18 24 -:  0  1    diagliso 5
''           -:  0  1  0 diagliso 5
  6 12       -:  0  1  2 diagliso 5
  6 12 18 24 -:  0  1  _ diagliso 5
  12 6       -:  0  1 _2 diagliso 5
  24 18 12 6 -:  0  1 __ diagliso 5
0 6 12 18    -:  0 _2    diagliso 5
''           -:  0 _2  0 diagliso 5
    12 18    -:  0 _2  2 diagliso 5
0 6 12 18    -:  0 _2  _ diagliso 5
    18 12    -:  0 _2 _2 diagliso 5
18 12 6 0    -:  0 _2 __ diagliso 5
5 11 17 23   -: _1       diagliso 5
5 11 17 23   -: _1  0    diagliso 5
''           -: _1  0  0 diagliso 5
5 11         -: _1  0  2 diagliso 5
5 11 17 23   -: _1  0  _ diagliso 5
11 5         -: _1  0 _2 diagliso 5
23 17 11 5   -: _1  0 __ diagliso 5
5 11 17 23   -: _1 _1    diagliso 5
''           -: _1 _1  0 diagliso 5
     17 23   -: _1 _1  2 diagliso 5
5 11 17 23   -: _1 _1  _ diagliso 5
     23 17   -: _1 _1 _2 diagliso 5
23 17 11 5   -: _1 _1 __ diagliso 5
  11 17 23   -: _1  1    diagliso 5
''           -: _1  1  0 diagliso 5
  11 17      -: _1  1  2 diagliso 5
  11 17 23   -: _1  1  _ diagliso 5
  17 11      -: _1  1 _2 diagliso 5
  23 17 11   -: _1  1 __ diagliso 5
5 11 17      -: _1 _2    diagliso 5
''           -: _1 _2  0 diagliso 5
  11 17      -: _1 _2  2 diagliso 5
5 11 17      -: _1 _2  _ diagliso 5
  17 11      -: _1 _2 _2 diagliso 5
17 11 5      -: _1 _2 __ diagliso 5
1 7 13 19    -:  1       diagliso 5
1 7 13 19    -:  1  0    diagliso 5
''           -:  1  0  0 diagliso 5
1 7          -:  1  0  2 diagliso 5
1 7 13 19    -:  1  0  _ diagliso 5
7 1          -:  1  0 _2 diagliso 5
19 13 7 1    -:  1  0 __ diagliso 5
1 7 13 19    -:  1 _1    diagliso 5
''           -:  1 _1  0 diagliso 5
    13 19    -:  1 _1  2 diagliso 5
1 7 13 19    -:  1 _1  _ diagliso 5
    19 13    -:  1 _1 _2 diagliso 5
19 13 7 1    -:  1 _1 __ diagliso 5
  7 13 19    -:  1  1    diagliso 5
''           -:  1  1  0 diagliso 5
  7 13       -:  1  1  2 diagliso 5
  7 13 19    -:  1  1  _ diagliso 5
  13 7       -:  1  1 _2 diagliso 5
  19 13 7    -:  1  1 __ diagliso 5
1 7 13       -:  1 _2    diagliso 5
''           -:  1 _2  0 diagliso 5
  7 13       -:  1 _2  2 diagliso 5
1 7 13       -:  1 _2  _ diagliso 5
  13 7       -:  1 _2 _2 diagliso 5
13 7 1       -:  1 _2 __ diagliso 5
2 8 14       -:  2       diagliso 5
2 8 14       -:  2  0    diagliso 5
''           -:  2  0  0 diagliso 5
2 8          -:  2  0  2 diagliso 5
2 8 14       -:  2  0  _ diagliso 5
8 2          -:  2  0 _2 diagliso 5
14 8 2       -:  2  0 __ diagliso 5
2 8 14       -:  2 _1    diagliso 5
''           -:  2 _1  0 diagliso 5
  8 14       -:  2 _1  2 diagliso 5
2 8 14       -:  2 _1  _ diagliso 5
  14 8       -:  2 _1 _2 diagliso 5
14 8 2       -:  2 _1 __ diagliso 5
  8 14       -:  2  1    diagliso 5
''           -:  2  1  0 diagliso 5
  8 14       -:  2  1  2 diagliso 5
  8 14       -:  2  1  _ diagliso 5
  14 8       -:  2  1 _2 diagliso 5
  14 8       -:  2  1 __ diagliso 5
2 8          -:  2 _2    diagliso 5
''           -:  2 _2  0 diagliso 5
2 8          -:  2 _2  2 diagliso 5
2 8          -:  2 _2  _ diagliso 5
8 2          -:  2 _2 _2 diagliso 5
8 2          -:  2 _2 __ diagliso 5
10 16 22     -: _2       diagliso 5
10 16 22     -: _2  0    diagliso 5
''           -: _2  0  0 diagliso 5
10 16        -: _2  0  2 diagliso 5
10 16 22     -: _2  0  _ diagliso 5
16 10        -: _2  0 _2 diagliso 5
22 16 10     -: _2  0 __ diagliso 5
10 16 22     -: _2 _1    diagliso 5
''           -: _2 _1  0 diagliso 5
16 22        -: _2 _1  2 diagliso 5
10 16 22     -: _2 _1  _ diagliso 5
22 16        -: _2 _1 _2 diagliso 5
22 16 10     -: _2 _1 __ diagliso 5
16 22        -: _2  1    diagliso 5
''           -: _2  1  0 diagliso 5
16 22        -: _2  1  2 diagliso 5
16 22        -: _2  1  _ diagliso 5
22 16        -: _2  1 _2 diagliso 5
22 16        -: _2  1 __ diagliso 5
10 16        -: _2 _2    diagliso 5
''           -: _2 _2  0 diagliso 5
10 16        -: _2 _2  2 diagliso 5
10 16        -: _2 _2  _ diagliso 5
16 10        -: _2 _2 _2 diagliso 5
16 10        -: _2 _2 __ diagliso 5

NB. c
0 -: c EMPTY
5 -: c i. 0 5
0 -: c i. 5 0
5 -: c i. 5
0 -: c i. 0
1 -: c i. i. 0
1 -: c 42
4 -: c a034

NB. trace
 0    -: trace EMPTY
18    -: trace ai35
18j63 -: trace ac35

NB. ct
(-: ct) EMPTY
(3 2 $ 0     3   1    4     2    5    ) -: ct ai23
(3 2 $ 0j_6 3j_9 1j_7 4j_10 2j_8 5j_11) -: ct ac23

NB. cp
(-: cp) EMPTY
(3 2 $ 5     2    4     1    3    0   ) -: cp ai23
(3 2 $ 5j_11 2j_8 4j_10 1j_7 3j_9 0j_6) -: cp ac23

NB. fp
(-: ''&fp    ) EMPTY
(-: ''&fp^:_1) EMPTY
(4 4 $ 10     9    11     8    6    5    7    4    14    13    15    12     2    1    3     0   ) -: p fp     ai44
(4 4 $ 10j26  9j25 11j27  8j24 6j22 5j21 7j23 4j20 14j30 13j29 15j31 12j28  2j18 1j17 3j19  0j16) -: p fp     ac44
(4 4 $ 15    13    12    14    7    5    4    6     3     1     0     2    11    9    8    10   ) -: p fp^:_1 ai44
(4 4 $ 15j31 13j29 12j28 14j30 7j23 5j21 4j20 6j22  3j19  1j17  0j16  2j18 11j27 9j25 8j24 10j26) -: p fp^:_1 ac44
ai44 -: p fp^:_1: p fp ai44
ac44 -: p fp^:_1: p fp ac44

NB. P4p, P4ip
EMPTY -:  P4p     ''
EMPTY -: P4ip     ''
''    -:  P4p^:_1 EMPTY
''    -: P4ip^:_1 EMPTY
(4 4 $ 0 0 1 0 0 1 0 0 0 0 0 1 1 0 0 0) -:  P
(4 4 $ 0 0 0 1 0 1 0 0 1 0 0 0 0 0 1 0) -: iP
P -: |: iP
p -:     P4p^:_1  P
p -: /: P4ip^:_1  P
p -: /:  P4p^:_1 iP

NB. icut
(-: ]&.(       fret &(<;.1) :. icut)) i. 10
(-: ]&.((2 # < fret)&(<;.1) :. icut)) i. 10 10
(-: ]&.((3 # < fret)&(<;.1) :. icut)) i. 10 10 10

NB. rt
(-:  _   &rt) ai34
(-: __   &rt) ai34
(-:  _  _&rt) ai34
(-:  _ __&rt) ai34
(-: __  _&rt) ai34
(-: __ __&rt) ai34
(2 4 $ 0 1 2 3 4 5  6  7   ) -:   2     rt ai34
(2 4 $ 4 5 6 7 8 9 10 11   ) -:  _2     rt ai34
(2 4 $ 0 1 2 3 4 5  6  7   ) -:   2  30 rt ai34
(2 4 $ 0 1 2 3 4 5  6  7   ) -:   2   _ rt ai34
(3 3 $ 0 1 2 4 5 6  8  9 10) -:  20   3 rt ai34
(3 3 $ 0 1 2 4 5 6  8  9 10) -:   _   3 rt ai34
(2 4 $ 4 5 6 7 8 9 10 11   ) -:  _2 _30 rt ai34
(2 4 $ 4 5 6 7 8 9 10 11   ) -:  _2  __ rt ai34
(3 3 $ 1 2 3 5 6 7  9 10 11) -: _20  _3 rt ai34
(3 3 $ 1 2 3 5 6 7  9 10 11) -:  __  _3 rt ai34

NB. e0
(-: 0  &e0) EMPTY
(-: 0  &e0) ai34
(-: 0 0&e0) ai34
(1 (< 3 ;&i.       4)} 5 6 $ 0) -:  2    e0 3 4 $ 1
(1 (< 3 ;&i.       4)} 5 7 $ 0) -:  2  3 e0 3 4 $ 1
(1 (< (; <@<) i.   3)} 5 7 $ 0) -:  2 _3 e0 3 4 $ 1
(1 (< ,~ <@<  i.   2)} 5 6 $ 0) -: _2    e0 3 4 $ 1
(1 (< 2 ;&(<@  i.) 4)} 5 7 $ 0) -: _2  3 e0 3 4 $ 1
(1 (< 2 ,&(<@<@i.) 3)} 5 7 $ 0) -: _2 _3 e0 3 4 $ 1

NB. appendx

NB. - one of arguments has rank 1
(,: '') -: ''    appendl EMPTY
(,: '') -: ''    appendr EMPTY
(,: '') -: EMPTY appendl ''
(,: '') -: EMPTY appendr ''
(1 3 $ 3) -: (0 2 $ 0) appendl 3 $ 3
(1 3 $ 3) -: (0 2 $ 0) appendr 3 $ 3
(1 3 $ 3) -: (3 $ 3) appendl 0 2 $ 0
(1 3 $ 3) -: (3 $ 3) appendr 0 2 $ 0
(1 3 $ 2 2 0) -: (0 3 $ 0) appendl 2 $ 2
(1 3 $ 0 2 2) -: (0 3 $ 0) appendr 2 $ 2
(1 3 $ 2 2 0) -: (2 $ 2) appendl 0 3 $ 0
(1 3 $ 0 2 2) -: (2 $ 2) appendr 0 3 $ 0
(3 3 $ 6 3 # 0 3) -: (2 0 $ 0) appendl 3 $ 3
(3 3 $ 6 3 # 0 3) -: (2 0 $ 0) appendr 3 $ 3
(4 2 $ 6 2 # 0 2) -: (3 0 $ 0) appendl 2 $ 2
(4 2 $ 6 2 # 0 2) -: (3 0 $ 0) appendr 2 $ 2
(4 2 $ 2 6 # 2 0) -: (2 $ 2) appendl 3 0 $ 0
(4 2 $ 2 6 # 2 0) -: (2 $ 2) appendr 3 0 $ 0
(3 3 $ 3 6 # 3 0) -: (3 $ 3) appendl 2 0 $ 0
(3 3 $ 3 6 # 3 0) -: (3 $ 3) appendr 2 0 $ 0
(3 3 3 ,  ,:~ 2 2 0) -: (3 $ 3) appendl 2 2 $ 2
(3 3 3 ,  ,:~ 0 2 2) -: (3 $ 3) appendr 2 2 $ 2
(3 3 3 ,~ ,:~ 2 2 0) -: (2 2 $ 2) appendl 3 $ 3
(3 3 3 ,~ ,:~ 0 2 2) -: (2 2 $ 2) appendr 3 $ 3
(4 3 $ 9 2 1 # 3 2 0) -: (3 3 $ 3) appendl 2 $ 2
(4 3 $ 9 1 2 # 3 0 2) -: (3 3 $ 3) appendr 2 $ 2
(4 3 $ 2 1 9 # 2 0 3) -: (2 $ 2) appendl 3 3 $ 3
(4 3 $ 1 2 9 # 0 2 3) -: (2 $ 2) appendr 3 3 $ 3

NB. - both arguments have rank 2
(-: appendl~) EMPTY
(-: appendr~) EMPTY
(3 3 $ 3) -: (0 2 $ 0) appendl 3 3 $ 3
(3 3 $ 3) -: (0 2 $ 0) appendr 3 3 $ 3
(3 3 $ 3) -: (3 3 $ 3) appendl 0 2 $ 0
(3 3 $ 3) -: (3 3 $ 3) appendr 0 2 $ 0
(,:~ 2 2 0) -: (0 3 $ 0) appendl 2 2 $ 2
(,:~ 2 2 0) -: (2 2 $ 2) appendl 0 3 $ 0
(,:~ 0 2 2) -: (0 3 $ 0) appendr 2 2 $ 2
(,:~ 0 2 2) -: (2 2 $ 2) appendr 0 3 $ 0
(5 3 $ 6 9 # 0 3) -: (2 0 $ 0) appendl 3 3 $ 3
(5 2 $ 6 4 # 0 2) -: (3 0 $ 0) appendl 2 2 $ 2
(5 3 $ 6 9 # 0 3) -: (2 0 $ 0) appendr 3 3 $ 3
(5 2 $ 6 4 # 0 2) -: (3 0 $ 0) appendr 2 2 $ 2
(5 2 $ 4 6 # 2 0) -: (2 2 $ 2) appendl 3 0 $ 0
(5 3 $ 9 6 # 3 0) -: (3 3 $ 3) appendl 2 0 $ 0
(5 2 $ 4 6 # 2 0) -: (2 2 $ 2) appendr 3 0 $ 0
(5 3 $ 9 6 # 3 0) -: (3 3 $ 3) appendr 2 0 $ 0
(3 2 # 2 2 0 ,:~ 3 3 3) -: (3 3 $ 3) appendl 2 2 $ 2
(2 3 # 2 2 0 ,:  3 3 3) -: (2 2 $ 2) appendl 3 3 $ 3
(3 2 # 0 2 2 ,:~ 3 3 3) -: (3 3 $ 3) appendr 2 2 $ 2
(2 3 # 0 2 2 ,:  3 3 3) -: (2 2 $ 2) appendr 3 3 $ 3

NB. stitchx

NB. - one of arguments has rank 1
(,. '') -: ''    stitcht EMPTY
(,. '') -: ''    stitchb EMPTY
(,. '') -: EMPTY stitcht ''
(,. '') -: EMPTY stitchb ''
(3 ,.~ 3 2 $ 0) -: (0 2 $ 0) stitcht 3 $ 3
(3 ,.~ 3 2 $ 0) -: (0 2 $ 0) stitchb 3 $ 3
(2 ,.~ 2 3 $ 0) -: (0 3 $ 0) stitcht 2 $ 2
(2 ,.~ 2 3 $ 0) -: (0 3 $ 0) stitchb 2 $ 2
(3 ,.  3 2 $ 0) -: (3 $ 3) stitcht 0 2 $ 0
(3 ,.  3 2 $ 0) -: (3 $ 3) stitchb 0 2 $ 0
(2 ,.  2 3 $ 0) -: (2 $ 2) stitcht 0 3 $ 0
(2 ,.  2 3 $ 0) -: (2 $ 2) stitchb 0 3 $ 0
(,. 3 3 3) -: (2 0 $ 0) stitcht 3 $ 3
(,. 3 3 3) -: (2 0 $ 0) stitchb 3 $ 3
(,. 2 2 0) -: (3 0 $ 0) stitcht 2 $ 2
(,. 0 2 2) -: (3 0 $ 0) stitchb 2 $ 2
(,. 2 2 0) -: (2 $ 2) stitcht 3 0 $ 0
(,. 0 2 2) -: (2 $ 2) stitchb 3 0 $ 0
(,. 3 3 3) -: (3 $ 3) stitcht 2 0 $ 0
(,. 3 3 3) -: (3 $ 3) stitchb 2 0 $ 0
(3 ,.  0 ,~ 2 2 $ 2) -: (3 $ 3) stitcht 2 2 $ 2
(3 ,.  0 ,  2 2 $ 2) -: (3 $ 3) stitchb 2 2 $ 2
(3 ,.~ 0 ,~ 2 2 $ 2) -: (2 2 $ 2) stitcht 3 $ 3
(3 ,.~ 0 ,  2 2 $ 2) -: (2 2 $ 2) stitchb 3 $ 3
(2 2 0 ,.~ 3 3 $ 3) -: (3 3 $ 3) stitcht 2 $ 2
(0 2 2 ,.~ 3 3 $ 3) -: (3 3 $ 3) stitchb 2 $ 2
(2 2 0 ,.  3 3 $ 3) -: (2 $ 2) stitcht 3 3 $ 3
(0 2 2 ,.  3 3 $ 3) -: (2 $ 2) stitchb 3 3 $ 3

NB. - both arguments have rank 2
(-: stitcht~) EMPTY
(-: stitchb~) EMPTY
(3 # ,: 2 3 # 0 3) -: (0 2 $ 0) stitcht 3 3 $ 3
(3 # ,: 2 3 # 0 3) -: (0 2 $ 0) stitchb 3 3 $ 3
(3 # ,: 3 2 # 3 0) -: (3 3 $ 3) stitcht 0 2 $ 0
(3 # ,: 3 2 # 3 0) -: (3 3 $ 3) stitchb 0 2 $ 0
(,:~ 3 2 # 0 2) -: (0 3 $ 0) stitcht 2 2 $ 2
(,:~ 3 2 # 0 2) -: (0 3 $ 0) stitchb 2 2 $ 2
(,:~ 2 3 # 2 0) -: (2 2 $ 2) stitcht 0 3 $ 0
(,:~ 2 3 # 2 0) -: (2 2 $ 2) stitchb 0 3 $ 0
(3 3 $ 3) -: (2 0 $ 0) stitcht 3 3 $ 3
(3 3 $ 3) -: (2 0 $ 0) stitchb 3 3 $ 3
(3 3 $ 3) -: (3 3 $ 3) stitcht 2 0 $ 0
(3 3 $ 3) -: (3 3 $ 3) stitchb 2 0 $ 0
(,.~ 2 2 0) -: (3 0 $ 0) stitcht 2 2 $ 2
(,.~ 2 2 0) -: (2 2 $ 2) stitcht 3 0 $ 0
(,.~ 0 2 2) -: (3 0 $ 0) stitchb 2 2 $ 2
(,.~ 0 2 2) -: (2 2 $ 2) stitchb 3 0 $ 0
((3 3 $ 3) ,.  ,.~ 2 2 0) -: (3 3 $ 3) stitcht 2 2 $ 2
((3 3 $ 3) ,.~ ,.~ 2 2 0) -: (2 2 $ 2) stitcht 3 3 $ 3
((3 3 $ 3) ,.  ,.~ 0 2 2) -: (3 3 $ 3) stitchb 2 2 $ 2
((3 3 $ 3) ,.~ ,.~ 0 2 2) -: (2 2 $ 2) stitchb 3 3 $ 3

NB. ds
(-: ds~) EMPTY
(2 2 $ 2) -: (2 2 $ 2) ds EMPTY
(2 2 $ 2) -: EMPTY     ds 2 2 $ 2
(2 3 # 2 2 0 0 0 ,: 0 0 3 3 3) -: (2 2 $ 2) ds 3 3 $ 3
(3 2 # 3 3 3 0 0 ,: 0 0 0 2 2) -: (3 3 $ 3) ds 2 2 $ 2

NB. diag

'' -:         diag EMPTY
'' -: 0       diag EMPTY
'' -: 0  0    diag EMPTY
'' -: 0  0  0 diag EMPTY
'' -: 0  0  _ diag EMPTY
'' -: 0  0 __ diag EMPTY
'' -: 0 _1    diag EMPTY
'' -: 0 _1  0 diag EMPTY
'' -: 0 _1  _ diag EMPTY
'' -: 0 _1 __ diag EMPTY

'' -:         diag a005
'' -: 0       diag a005
'' -: 0  0    diag a005
'' -: 0  0  0 diag a005
'' -: 0  0  _ diag a005
'' -: 0  0 __ diag a005
'' -: 0 _1    diag a005
'' -: 0 _1  0 diag a005
'' -: 0 _1  _ diag a005
'' -: 0 _1 __ diag a005

'' -:         diag a050
'' -: 0       diag a050
'' -: 0  0    diag a050
'' -: 0  0  0 diag a050
'' -: 0  0  _ diag a050
'' -: 0  0 __ diag a050
'' -: 0 _1    diag a050
'' -: 0 _1  0 diag a050
'' -: 0 _1  _ diag a050
'' -: 0 _1 __ diag a050

0 6 12 18 -:          diag ai45
0 6 12 18 -:  0       diag ai45
0 6 12 18 -:  0  0    diag ai45
''        -:  0  0  0 diag ai45
0 6       -:  0  0  2 diag ai45
0 6 12 18 -:  0  0  _ diag ai45
6 0       -:  0  0 _2 diag ai45
18 12 6 0 -:  0  0 __ diag ai45
0 6 12 18 -:  0 _1    diag ai45
''        -:  0 _1  0 diag ai45
    12 18 -:  0 _1  2 diag ai45
0 6 12 18 -:  0 _1  _ diag ai45
    18 12 -:  0 _1 _2 diag ai45
18 12 6 0 -:  0 _1 __ diag ai45
  6 12 18 -:  0  1    diag ai45
''        -:  0  1  0 diag ai45
  6 12    -:  0  1  2 diag ai45
  6 12 18 -:  0  1  _ diag ai45
  12 6    -:  0  1 _2 diag ai45
  18 12 6 -:  0  1 __ diag ai45
0 6 12    -:  0 _2    diag ai45
''        -:  0 _2  0 diag ai45
  6 12    -:  0 _2  2 diag ai45
0 6 12    -:  0 _2  _ diag ai45
  12 6    -:  0 _2 _2 diag ai45
12 6 0    -:  0 _2 __ diag ai45
5 11 17   -: _1       diag ai45
5 11 17   -: _1  0    diag ai45
''        -: _1  0  0 diag ai45
5 11      -: _1  0  2 diag ai45
5 11 17   -: _1  0  _ diag ai45
11 5      -: _1  0 _2 diag ai45
17 11 5   -: _1  0 __ diag ai45
5 11 17   -: _1 _1    diag ai45
''        -: _1 _1  0 diag ai45
  11 17   -: _1 _1  2 diag ai45
5 11 17   -: _1 _1  _ diag ai45
  17 11   -: _1 _1 _2 diag ai45
17 11 5   -: _1 _1 __ diag ai45
  11 17   -: _1  1    diag ai45
''        -: _1  1  0 diag ai45
  11 17   -: _1  1  2 diag ai45
  11 17   -: _1  1  _ diag ai45
  17 11   -: _1  1 _2 diag ai45
  17 11   -: _1  1 __ diag ai45
5 11      -: _1 _2    diag ai45
''        -: _1 _2  0 diag ai45
5 11      -: _1 _2  2 diag ai45
5 11      -: _1 _2  _ diag ai45
11 5      -: _1 _2 _2 diag ai45
11 5      -: _1 _2 __ diag ai45
1 7 13 19 -:  1       diag ai45
1 7 13 19 -:  1  0    diag ai45
''        -:  1  0  0 diag ai45
1 7       -:  1  0  2 diag ai45
1 7 13 19 -:  1  0  _ diag ai45
7 1       -:  1  0 _2 diag ai45
19 13 7 1 -:  1  0 __ diag ai45
1 7 13 19 -:  1 _1    diag ai45
''        -:  1 _1  0 diag ai45
    13 19 -:  1 _1  2 diag ai45
1 7 13 19 -:  1 _1  _ diag ai45
    19 13 -:  1 _1 _2 diag ai45
19 13 7 1 -:  1 _1 __ diag ai45
  7 13 19 -:  1  1    diag ai45
''        -:  1  1  0 diag ai45
  7 13    -:  1  1  2 diag ai45
  7 13 19 -:  1  1  _ diag ai45
  13 7    -:  1  1 _2 diag ai45
  19 13 7 -:  1  1 __ diag ai45
1 7 13    -:  1 _2    diag ai45
''        -:  1 _2  0 diag ai45
  7 13    -:  1 _2  2 diag ai45
1 7 13    -:  1 _2  _ diag ai45
  13 7    -:  1 _2 _2 diag ai45
13 7 1    -:  1 _2 __ diag ai45
2 8 14    -:  2       diag ai45
2 8 14    -:  2  0    diag ai45
''        -:  2  0  0 diag ai45
2 8       -:  2  0  2 diag ai45
2 8 14    -:  2  0  _ diag ai45
8 2       -:  2  0 _2 diag ai45
14 8 2    -:  2  0 __ diag ai45
2 8 14    -:  2 _1    diag ai45
''        -:  2 _1  0 diag ai45
  8 14    -:  2 _1  2 diag ai45
2 8 14    -:  2 _1  _ diag ai45
  14 8    -:  2 _1 _2 diag ai45
14 8 2    -:  2 _1 __ diag ai45
  8 14    -:  2  1    diag ai45
''        -:  2  1  0 diag ai45
  8 14    -:  2  1  2 diag ai45
  8 14    -:  2  1  _ diag ai45
  14 8    -:  2  1 _2 diag ai45
  14 8    -:  2  1 __ diag ai45
2 8       -:  2 _2    diag ai45
''        -:  2 _2  0 diag ai45
2 8       -:  2 _2  2 diag ai45
2 8       -:  2 _2  _ diag ai45
8 2       -:  2 _2 _2 diag ai45
8 2       -:  2 _2 __ diag ai45
10 16     -: _2       diag ai45
10 16     -: _2  0    diag ai45
''        -: _2  0  0 diag ai45
10 16     -: _2  0  2 diag ai45
10 16     -: _2  0  _ diag ai45
16 10     -: _2  0 _2 diag ai45
16 10     -: _2  0 __ diag ai45
10 16     -: _2 _1    diag ai45
''        -: _2 _1  0 diag ai45
10 16     -: _2 _1  2 diag ai45
10 16     -: _2 _1  _ diag ai45
16 10     -: _2 _1 _2 diag ai45
16 10     -: _2 _1 __ diag ai45
(, 16)    -: _2  1    diag ai45
''        -: _2  1  0 diag ai45
(, 16)    -: _2  1  2 diag ai45
(, 16)    -: _2  1  _ diag ai45
(, 16)    -: _2  1 _2 diag ai45
(, 16)    -: _2  1 __ diag ai45
(, 10)    -: _2 _2    diag ai45
''        -: _2 _2  0 diag ai45
(, 10)    -: _2 _2  2 diag ai45
(, 10)    -: _2 _2  _ diag ai45
(, 10)    -: _2 _2 _2 diag ai45
(, 10)    -: _2 _2 __ diag ai45

0 5 10 15 -:          diag ai54
0 5 10 15 -:  0       diag ai54
0 5 10 15 -:  0  0    diag ai54
''        -:  0  0  0 diag ai54
0 5       -:  0  0  2 diag ai54
0 5 10 15 -:  0  0  _ diag ai54
5 0       -:  0  0 _2 diag ai54
15 10 5 0 -:  0  0 __ diag ai54
0 5 10 15 -:  0 _1    diag ai54
''        -:  0 _1  0 diag ai54
    10 15 -:  0 _1  2 diag ai54
0 5 10 15 -:  0 _1  _ diag ai54
    15 10 -:  0 _1 _2 diag ai54
15 10 5 0 -:  0 _1 __ diag ai54
  5 10 15 -:  0  1    diag ai54
''        -:  0  1  0 diag ai54
  5 10    -:  0  1  2 diag ai54
  5 10 15 -:  0  1  _ diag ai54
  10 5    -:  0  1 _2 diag ai54
  15 10 5 -:  0  1 __ diag ai54
0 5 10    -:  0 _2    diag ai54
''        -:  0 _2  0 diag ai54
  5 10    -:  0 _2  2 diag ai54
0 5 10    -:  0 _2  _ diag ai54
  10 5    -:  0 _2 _2 diag ai54
10 5 0    -:  0 _2 __ diag ai54
4 9 14 19 -: _1       diag ai54
4 9 14 19 -: _1  0    diag ai54
''        -: _1  0  0 diag ai54
4 9       -: _1  0  2 diag ai54
4 9 14 19 -: _1  0  _ diag ai54
9 4       -: _1  0 _2 diag ai54
19 14 9 4 -: _1  0 __ diag ai54
4 9 14 19 -: _1 _1    diag ai54
''        -: _1 _1  0 diag ai54
    14 19 -: _1 _1  2 diag ai54
4 9 14 19 -: _1 _1  _ diag ai54
    19 14 -: _1 _1 _2 diag ai54
19 14 9 4 -: _1 _1 __ diag ai54
  9 14 19 -: _1  1    diag ai54
''        -: _1  1  0 diag ai54
  9 14    -: _1  1  2 diag ai54
  9 14 19 -: _1  1  _ diag ai54
  14 9    -: _1  1 _2 diag ai54
  19 14 9 -: _1  1 __ diag ai54
4 9 14    -: _1 _2    diag ai54
''        -: _1 _2  0 diag ai54
  9 14    -: _1 _2  2 diag ai54
4 9 14    -: _1 _2  _ diag ai54
  14 9    -: _1 _2 _2 diag ai54
14 9 4    -: _1 _2 __ diag ai54
1 6 11    -:  1       diag ai54
1 6 11    -:  1  0    diag ai54
''        -:  1  0  0 diag ai54
1 6       -:  1  0  2 diag ai54
1 6 11    -:  1  0  _ diag ai54
6 1       -:  1  0 _2 diag ai54
11 6 1    -:  1  0 __ diag ai54
1 6 11    -:  1 _1    diag ai54
''        -:  1 _1  0 diag ai54
  6 11    -:  1 _1  2 diag ai54
1 6 11    -:  1 _1  _ diag ai54
  11 6    -:  1 _1 _2 diag ai54
11 6 1    -:  1 _1 __ diag ai54
  6 11    -:  1  1    diag ai54
''        -:  1  1  0 diag ai54
  6 11    -:  1  1  2 diag ai54
  6 11    -:  1  1  _ diag ai54
  11 6    -:  1  1 _2 diag ai54
  11 6    -:  1  1 __ diag ai54
1 6       -:  1 _2    diag ai54
''        -:  1 _2  0 diag ai54
1 6       -:  1 _2  2 diag ai54
1 6       -:  1 _2  _ diag ai54
6 1       -:  1 _2 _2 diag ai54
6 1       -:  1 _2 __ diag ai54
2 7       -:  2       diag ai54
2 7       -:  2  0    diag ai54
''        -:  2  0  0 diag ai54
2 7       -:  2  0  2 diag ai54
2 7       -:  2  0  _ diag ai54
7 2       -:  2  0 _2 diag ai54
7 2       -:  2  0 __ diag ai54
2 7       -:  2 _1    diag ai54
''        -:  2 _1  0 diag ai54
2 7       -:  2 _1  2 diag ai54
2 7       -:  2 _1  _ diag ai54
7 2       -:  2 _1 _2 diag ai54
7 2       -:  2 _1 __ diag ai54
  (, 7)   -:  2  1    diag ai54
''        -:  2  1  0 diag ai54
  (, 7)   -:  2  1  2 diag ai54
  (, 7)   -:  2  1  _ diag ai54
  (, 7)   -:  2  1 _2 diag ai54
  (, 7)   -:  2  1 __ diag ai54
(, 2)     -:  2 _2    diag ai54
''        -:  2 _2  0 diag ai54
(, 2)     -:  2 _2  2 diag ai54
(, 2)     -:  2 _2  _ diag ai54
(, 2)     -:  2 _2 _2 diag ai54
(, 2)     -:  2 _2 __ diag ai54
8 13 18   -: _2       diag ai54
8 13 18   -: _2  0    diag ai54
''        -: _2  0  0 diag ai54
8 13      -: _2  0  2 diag ai54
8 13 18   -: _2  0  _ diag ai54
13 8      -: _2  0 _2 diag ai54
18 13 8   -: _2  0 __ diag ai54
8 13 18   -: _2 _1    diag ai54
''        -: _2 _1  0 diag ai54
  13 18   -: _2 _1  2 diag ai54
8 13 18   -: _2 _1  _ diag ai54
  18 13   -: _2 _1 _2 diag ai54
18 13 8   -: _2 _1 __ diag ai54
  13 18   -: _2  1    diag ai54
''        -: _2  1  0 diag ai54
  13 18   -: _2  1  2 diag ai54
  13 18   -: _2  1  _ diag ai54
  18 13   -: _2  1 _2 diag ai54
  18 13   -: _2  1 __ diag ai54
8 13      -: _2 _2    diag ai54
''        -: _2 _2  0 diag ai54
8 13      -: _2 _2  2 diag ai54
8 13      -: _2 _2  _ diag ai54
13 8      -: _2 _2 _2 diag ai54
13 8      -: _2 _2 __ diag ai54

0 6 12 18 24 -:          diag ai55
0 6 12 18 24 -:  0       diag ai55
0 6 12 18 24 -:  0  0    diag ai55
''           -:  0  0  0 diag ai55
0 6          -:  0  0  2 diag ai55
0 6 12 18 24 -:  0  0  _ diag ai55
6 0          -:  0  0 _2 diag ai55
24 18 12 6 0 -:  0  0 __ diag ai55
0 6 12 18 24 -:  0 _1    diag ai55
''           -:  0 _1  0 diag ai55
       18 24 -:  0 _1  2 diag ai55
0 6 12 18 24 -:  0 _1  _ diag ai55
       24 18 -:  0 _1 _2 diag ai55
24 18 12 6 0 -:  0 _1 __ diag ai55
  6 12 18 24 -:  0  1    diag ai55
''           -:  0  1  0 diag ai55
  6 12       -:  0  1  2 diag ai55
  6 12 18 24 -:  0  1  _ diag ai55
  12 6       -:  0  1 _2 diag ai55
  24 18 12 6 -:  0  1 __ diag ai55
0 6 12 18    -:  0 _2    diag ai55
''           -:  0 _2  0 diag ai55
    12 18    -:  0 _2  2 diag ai55
0 6 12 18    -:  0 _2  _ diag ai55
    18 12    -:  0 _2 _2 diag ai55
18 12 6 0    -:  0 _2 __ diag ai55
5 11 17 23   -: _1       diag ai55
5 11 17 23   -: _1  0    diag ai55
''           -: _1  0  0 diag ai55
5 11         -: _1  0  2 diag ai55
5 11 17 23   -: _1  0  _ diag ai55
11 5         -: _1  0 _2 diag ai55
23 17 11 5   -: _1  0 __ diag ai55
5 11 17 23   -: _1 _1    diag ai55
''           -: _1 _1  0 diag ai55
     17 23   -: _1 _1  2 diag ai55
5 11 17 23   -: _1 _1  _ diag ai55
     23 17   -: _1 _1 _2 diag ai55
23 17 11 5   -: _1 _1 __ diag ai55
  11 17 23   -: _1  1    diag ai55
''           -: _1  1  0 diag ai55
  11 17      -: _1  1  2 diag ai55
  11 17 23   -: _1  1  _ diag ai55
  17 11      -: _1  1 _2 diag ai55
  23 17 11   -: _1  1 __ diag ai55
5 11 17      -: _1 _2    diag ai55
''           -: _1 _2  0 diag ai55
  11 17      -: _1 _2  2 diag ai55
5 11 17      -: _1 _2  _ diag ai55
  17 11      -: _1 _2 _2 diag ai55
17 11 5      -: _1 _2 __ diag ai55
1 7 13 19    -:  1       diag ai55
1 7 13 19    -:  1  0    diag ai55
''           -:  1  0  0 diag ai55
1 7          -:  1  0  2 diag ai55
1 7 13 19    -:  1  0  _ diag ai55
7 1          -:  1  0 _2 diag ai55
19 13 7 1    -:  1  0 __ diag ai55
1 7 13 19    -:  1 _1    diag ai55
''           -:  1 _1  0 diag ai55
    13 19    -:  1 _1  2 diag ai55
1 7 13 19    -:  1 _1  _ diag ai55
    19 13    -:  1 _1 _2 diag ai55
19 13 7 1    -:  1 _1 __ diag ai55
  7 13 19    -:  1  1    diag ai55
''           -:  1  1  0 diag ai55
  7 13       -:  1  1  2 diag ai55
  7 13 19    -:  1  1  _ diag ai55
  13 7       -:  1  1 _2 diag ai55
  19 13 7    -:  1  1 __ diag ai55
1 7 13       -:  1 _2    diag ai55
''           -:  1 _2  0 diag ai55
  7 13       -:  1 _2  2 diag ai55
1 7 13       -:  1 _2  _ diag ai55
  13 7       -:  1 _2 _2 diag ai55
13 7 1       -:  1 _2 __ diag ai55
2 8 14       -:  2       diag ai55
2 8 14       -:  2  0    diag ai55
''           -:  2  0  0 diag ai55
2 8          -:  2  0  2 diag ai55
2 8 14       -:  2  0  _ diag ai55
8 2          -:  2  0 _2 diag ai55
14 8 2       -:  2  0 __ diag ai55
2 8 14       -:  2 _1    diag ai55
''           -:  2 _1  0 diag ai55
  8 14       -:  2 _1  2 diag ai55
2 8 14       -:  2 _1  _ diag ai55
  14 8       -:  2 _1 _2 diag ai55
14 8 2       -:  2 _1 __ diag ai55
  8 14       -:  2  1    diag ai55
''           -:  2  1  0 diag ai55
  8 14       -:  2  1  2 diag ai55
  8 14       -:  2  1  _ diag ai55
  14 8       -:  2  1 _2 diag ai55
  14 8       -:  2  1 __ diag ai55
2 8          -:  2 _2    diag ai55
''           -:  2 _2  0 diag ai55
2 8          -:  2 _2  2 diag ai55
2 8          -:  2 _2  _ diag ai55
8 2          -:  2 _2 _2 diag ai55
8 2          -:  2 _2 __ diag ai55
10 16 22     -: _2       diag ai55
10 16 22     -: _2  0    diag ai55
''           -: _2  0  0 diag ai55
10 16        -: _2  0  2 diag ai55
10 16 22     -: _2  0  _ diag ai55
16 10        -: _2  0 _2 diag ai55
22 16 10     -: _2  0 __ diag ai55
10 16 22     -: _2 _1    diag ai55
''           -: _2 _1  0 diag ai55
16 22        -: _2 _1  2 diag ai55
10 16 22     -: _2 _1  _ diag ai55
22 16        -: _2 _1 _2 diag ai55
22 16 10     -: _2 _1 __ diag ai55
16 22        -: _2  1    diag ai55
''           -: _2  1  0 diag ai55
16 22        -: _2  1  2 diag ai55
16 22        -: _2  1  _ diag ai55
22 16        -: _2  1 _2 diag ai55
22 16        -: _2  1 __ diag ai55
10 16        -: _2 _2    diag ai55
''           -: _2 _2  0 diag ai55
10 16        -: _2 _2  2 diag ai55
10 16        -: _2 _2  _ diag ai55
16 10        -: _2 _2 _2 diag ai55
16 10        -: _2 _2 __ diag ai55

NB. setdiag

NB. empty matrix

(-: (;~ ''       )&setdiag) EMPTY
(-: ('' ; 0      )&setdiag) EMPTY
(-: ('' ; 0  0   )&setdiag) EMPTY
(-: ('' ; 0  0  0)&setdiag) EMPTY
(-: ('' ; 0  0  _)&setdiag) EMPTY
(-: ('' ; 0  0 __)&setdiag) EMPTY
(-: ('' ; 0 _1   )&setdiag) EMPTY
(-: ('' ; 0 _1  0)&setdiag) EMPTY
(-: ('' ; 0 _1  _)&setdiag) EMPTY
(-: ('' ; 0 _1 __)&setdiag) EMPTY

(-: (;~ ''       )&setdiag) a005
(-: ('' ; 0      )&setdiag) a005
(-: ('' ; 0  0   )&setdiag) a005
(-: ('' ; 0  0  0)&setdiag) a005
(-: ('' ; 0  0  _)&setdiag) a005
(-: ('' ; 0  0 __)&setdiag) a005
(-: ('' ; 0 _1   )&setdiag) a005
(-: ('' ; 0 _1  0)&setdiag) a005
(-: ('' ; 0 _1  _)&setdiag) a005
(-: ('' ; 0 _1 __)&setdiag) a005

(-: (;~ ''       )&setdiag) a050
(-: ('' ; 0      )&setdiag) a050
(-: ('' ; 0  0   )&setdiag) a050
(-: ('' ; 0  0  0)&setdiag) a050
(-: ('' ; 0  0  _)&setdiag) a050
(-: ('' ; 0  0 __)&setdiag) a050
(-: ('' ; 0 _1   )&setdiag) a050
(-: ('' ; 0 _1  0)&setdiag) a050
(-: ('' ; 0 _1  _)&setdiag) a050
(-: ('' ; 0 _1 __)&setdiag) a050

NB. scalar e

(_1&(( ,.~      i.         4)}) -: (_1 ; ''      )&setdiag) ai45
(_1&(( ,.~      i.         4)}) -: (_1 ;  0      )&setdiag) ai45
(_1&(( ,.~      i.         4)}) -: (_1 ;  0  0   )&setdiag) ai45
(                               -: (_1 ;  0  0  0)&setdiag) ai45
(_1&(( ,.~         0 1      )}) -: (_1 ;  0  0  2)&setdiag) ai45
(_1&(( ,.~      i.         4)}) -: (_1 ;  0  0  _)&setdiag) ai45
(_1&(( ,.~         0 1      )}) -: (_1 ;  0  0 _2)&setdiag) ai45
(_1&(( ,.~      i.         4)}) -: (_1 ;  0  0 __)&setdiag) ai45
(_1&(( ,.~      i.         4)}) -: (_1 ;  0 _1   )&setdiag) ai45
(                               -: (_1 ;  0 _1  0)&setdiag) ai45
(_1&(( ,.~             2 3  )}) -: (_1 ;  0 _1  2)&setdiag) ai45
(_1&(( ,.~      i.         4)}) -: (_1 ;  0 _1  _)&setdiag) ai45
(_1&(( ,.~             2 3  )}) -: (_1 ;  0 _1 _2)&setdiag) ai45
(_1&(( ,.~      i.         4)}) -: (_1 ;  0 _1 __)&setdiag) ai45
(_1&(( ,.~           1 2 3  )}) -: (_1 ;  0  1   )&setdiag) ai45
(                               -: (_1 ;  0  1  0)&setdiag) ai45
(_1&(( ,.~           1 2    )}) -: (_1 ;  0  1  2)&setdiag) ai45
(_1&(( ,.~           1 2 3  )}) -: (_1 ;  0  1  _)&setdiag) ai45
(_1&(( ,.~           1 2    )}) -: (_1 ;  0  1 _2)&setdiag) ai45
(_1&(( ,.~           1 2 3  )}) -: (_1 ;  0  1 __)&setdiag) ai45
(_1&(( ,.~      i.       3  )}) -: (_1 ;  0 _2   )&setdiag) ai45
(                               -: (_1 ;  0 _2  0)&setdiag) ai45
(_1&(( ,.~           1 2    )}) -: (_1 ;  0 _2  2)&setdiag) ai45
(_1&(( ,.~      i.       3  )}) -: (_1 ;  0 _2  _)&setdiag) ai45
(_1&(( ,.~           1 2    )}) -: (_1 ;  0 _2 _2)&setdiag) ai45
(_1&(( ,.~      i.       3  )}) -: (_1 ;  0 _2 __)&setdiag) ai45
(_1&(((,.~ >: ) i.       3  )}) -: (_1 ; _1      )&setdiag) ai45
(_1&(((,.~ >: ) i.       3  )}) -: (_1 ; _1  0   )&setdiag) ai45
(                               -: (_1 ; _1  0  0)&setdiag) ai45
(_1&(((,.~ >: )    0 1      )}) -: (_1 ; _1  0  2)&setdiag) ai45
(_1&(((,.~ >: ) i.       3  )}) -: (_1 ; _1  0  _)&setdiag) ai45
(_1&(((,.~ >: )    0 1      )}) -: (_1 ; _1  0 _2)&setdiag) ai45
(_1&(((,.~ >: ) i.       3  )}) -: (_1 ; _1  0 __)&setdiag) ai45
(_1&(((,.~ >: ) i.       3  )}) -: (_1 ; _1 _1   )&setdiag) ai45
(                               -: (_1 ; _1 _1  0)&setdiag) ai45
(_1&(((,.~ >: )      1 2    )}) -: (_1 ; _1 _1  2)&setdiag) ai45
(_1&(((,.~ >: ) i.       3  )}) -: (_1 ; _1 _1  _)&setdiag) ai45
(_1&(((,.~ >: )      1 2    )}) -: (_1 ; _1 _1 _2)&setdiag) ai45
(_1&(((,.~ >: ) i.       3  )}) -: (_1 ; _1 _1 __)&setdiag) ai45
(_1&(((,.~ >: )      1 2    )}) -: (_1 ; _1  1   )&setdiag) ai45
(                               -: (_1 ; _1  1  0)&setdiag) ai45
(_1&(((,.~ >: )      1 2    )}) -: (_1 ; _1  1  2)&setdiag) ai45
(_1&(((,.~ >: )      1 2    )}) -: (_1 ; _1  1  _)&setdiag) ai45
(_1&(((,.~ >: )      1 2    )}) -: (_1 ; _1  1 _2)&setdiag) ai45
(_1&(((,.~ >: )      1 2    )}) -: (_1 ; _1  1 __)&setdiag) ai45
(_1&(((,.~ >: )    0 1      )}) -: (_1 ; _1 _2   )&setdiag) ai45
(                               -: (_1 ; _1 _2  0)&setdiag) ai45
(_1&(((,.~ >: )    0 1      )}) -: (_1 ; _1 _2  2)&setdiag) ai45
(_1&(((,.~ >: )    0 1      )}) -: (_1 ; _1 _2  _)&setdiag) ai45
(_1&(((,.~ >: )    0 1      )}) -: (_1 ; _1 _2 _2)&setdiag) ai45
(_1&(((,.~ >: )    0 1      )}) -: (_1 ; _1 _2 __)&setdiag) ai45
(_1&(((,.  >: ) i.         4)}) -: (_1 ;  1      )&setdiag) ai45
(_1&(((,.  >: ) i.         4)}) -: (_1 ;  1  0   )&setdiag) ai45
(                               -: (_1 ;  1  0  0)&setdiag) ai45
(_1&(((,.  >: )    0 1      )}) -: (_1 ;  1  0  2)&setdiag) ai45
(_1&(((,.  >: ) i.         4)}) -: (_1 ;  1  0  _)&setdiag) ai45
(_1&(((,.  >: )    0 1      )}) -: (_1 ;  1  0 _2)&setdiag) ai45
(_1&(((,.  >: ) i.         4)}) -: (_1 ;  1  0 __)&setdiag) ai45
(_1&(((,.  >: ) i.         4)}) -: (_1 ;  1 _1   )&setdiag) ai45
(                               -: (_1 ;  1 _1  0)&setdiag) ai45
(_1&(((,.  >: )        2 3  )}) -: (_1 ;  1 _1  2)&setdiag) ai45
(_1&(((,.  >: ) i.         4)}) -: (_1 ;  1 _1  _)&setdiag) ai45
(_1&(((,.  >: )        2 3  )}) -: (_1 ;  1 _1 _2)&setdiag) ai45
(_1&(((,.  >: ) i.         4)}) -: (_1 ;  1 _1 __)&setdiag) ai45
(_1&(((,.  >: )      1 2 3  )}) -: (_1 ;  1  1   )&setdiag) ai45
(                               -: (_1 ;  1  1  0)&setdiag) ai45
(_1&(((,.  >: )      1 2    )}) -: (_1 ;  1  1  2)&setdiag) ai45
(_1&(((,.  >: )      1 2 3  )}) -: (_1 ;  1  1  _)&setdiag) ai45
(_1&(((,.  >: )      1 2    )}) -: (_1 ;  1  1 _2)&setdiag) ai45
(_1&(((,.  >: )      1 2 3  )}) -: (_1 ;  1  1 __)&setdiag) ai45
(_1&(((,.  >: ) i.       3  )}) -: (_1 ;  1 _2   )&setdiag) ai45
(                               -: (_1 ;  1 _2  0)&setdiag) ai45
(_1&(((,.  >: )      1 2    )}) -: (_1 ;  1 _2  2)&setdiag) ai45
(_1&(((,.  >: ) i.       3  )}) -: (_1 ;  1 _2  _)&setdiag) ai45
(_1&(((,.  >: )      1 2    )}) -: (_1 ;  1 _2 _2)&setdiag) ai45
(_1&(((,.  >: ) i.       3  )}) -: (_1 ;  1 _2 __)&setdiag) ai45
(_1&(((,.  2&+) i.       3  )}) -: (_1 ;  2      )&setdiag) ai45
(_1&(((,.  2&+) i.       3  )}) -: (_1 ;  2  0   )&setdiag) ai45
(                               -: (_1 ;  2  0  0)&setdiag) ai45
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2  0  2)&setdiag) ai45
(_1&(((,.  2&+) i.       3  )}) -: (_1 ;  2  0  _)&setdiag) ai45
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2  0 _2)&setdiag) ai45
(_1&(((,.  2&+) i.       3  )}) -: (_1 ;  2  0 __)&setdiag) ai45
(_1&(((,.  2&+) i.       3  )}) -: (_1 ;  2 _1   )&setdiag) ai45
(                               -: (_1 ;  2 _1  0)&setdiag) ai45
(_1&(((,.  2&+)      1 2    )}) -: (_1 ;  2 _1  2)&setdiag) ai45
(_1&(((,.  2&+) i.       3  )}) -: (_1 ;  2 _1  _)&setdiag) ai45
(_1&(((,.  2&+)      1 2    )}) -: (_1 ;  2 _1 _2)&setdiag) ai45
(_1&(((,.  2&+) i.       3  )}) -: (_1 ;  2 _1 __)&setdiag) ai45
(_1&(((,.  2&+)      1 2    )}) -: (_1 ;  2  1   )&setdiag) ai45
(                               -: (_1 ;  2  1  0)&setdiag) ai45
(_1&(((,.  2&+)      1 2    )}) -: (_1 ;  2  1  2)&setdiag) ai45
(_1&(((,.  2&+)      1 2    )}) -: (_1 ;  2  1  _)&setdiag) ai45
(_1&(((,.  2&+)      1 2    )}) -: (_1 ;  2  1 _2)&setdiag) ai45
(_1&(((,.  2&+)      1 2    )}) -: (_1 ;  2  1 __)&setdiag) ai45
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2 _2   )&setdiag) ai45
(                               -: (_1 ;  2 _2  0)&setdiag) ai45
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2 _2  2)&setdiag) ai45
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2 _2  _)&setdiag) ai45
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2 _2 _2)&setdiag) ai45
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2 _2 __)&setdiag) ai45
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2      )&setdiag) ai45
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2  0   )&setdiag) ai45
(                               -: (_1 ; _2  0  0)&setdiag) ai45
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2  0  2)&setdiag) ai45
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2  0  _)&setdiag) ai45
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2  0 _2)&setdiag) ai45
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2  0 __)&setdiag) ai45
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2 _1   )&setdiag) ai45
(                               -: (_1 ; _2 _1  0)&setdiag) ai45
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2 _1  2)&setdiag) ai45
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2 _1  _)&setdiag) ai45
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2 _1 _2)&setdiag) ai45
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2 _1 __)&setdiag) ai45
(_1&((< 3 1                 )}) -: (_1 ; _2  1   )&setdiag) ai45
(                               -: (_1 ; _2  1  0)&setdiag) ai45
(_1&((< 3 1                 )}) -: (_1 ; _2  1  2)&setdiag) ai45
(_1&((< 3 1                 )}) -: (_1 ; _2  1  _)&setdiag) ai45
(_1&((< 3 1                 )}) -: (_1 ; _2  1 _2)&setdiag) ai45
(_1&((< 3 1                 )}) -: (_1 ; _2  1 __)&setdiag) ai45
(_1&((< 2 0                 )}) -: (_1 ; _2 _2   )&setdiag) ai45
(                               -: (_1 ; _2 _2  0)&setdiag) ai45
(_1&((< 2 0                 )}) -: (_1 ; _2 _2  2)&setdiag) ai45
(_1&((< 2 0                 )}) -: (_1 ; _2 _2  _)&setdiag) ai45
(_1&((< 2 0                 )}) -: (_1 ; _2 _2 _2)&setdiag) ai45
(_1&((< 2 0                 )}) -: (_1 ; _2 _2 __)&setdiag) ai45

(_1&(( ,.~      i.         4)}) -: (_1 ; ''      )&setdiag) ai54
(_1&(( ,.~      i.         4)}) -: (_1 ;  0      )&setdiag) ai54
(_1&(( ,.~      i.         4)}) -: (_1 ;  0  0   )&setdiag) ai54
(                               -: (_1 ;  0  0  0)&setdiag) ai54
(_1&(( ,.~         0 1      )}) -: (_1 ;  0  0  2)&setdiag) ai54
(_1&(( ,.~      i.         4)}) -: (_1 ;  0  0  _)&setdiag) ai54
(_1&(( ,.~         0 1      )}) -: (_1 ;  0  0 _2)&setdiag) ai54
(_1&(( ,.~      i.         4)}) -: (_1 ;  0  0 __)&setdiag) ai54
(_1&(( ,.~      i.         4)}) -: (_1 ;  0 _1   )&setdiag) ai54
(                               -: (_1 ;  0 _1  0)&setdiag) ai54
(_1&(( ,.~             2 3  )}) -: (_1 ;  0 _1  2)&setdiag) ai54
(_1&(( ,.~      i.         4)}) -: (_1 ;  0 _1  _)&setdiag) ai54
(_1&(( ,.~             2 3  )}) -: (_1 ;  0 _1 _2)&setdiag) ai54
(_1&(( ,.~      i.         4)}) -: (_1 ;  0 _1 __)&setdiag) ai54
(_1&(( ,.~           1 2 3  )}) -: (_1 ;  0  1   )&setdiag) ai54
(                               -: (_1 ;  0  1  0)&setdiag) ai54
(_1&(( ,.~           1 2    )}) -: (_1 ;  0  1  2)&setdiag) ai54
(_1&(( ,.~           1 2 3  )}) -: (_1 ;  0  1  _)&setdiag) ai54
(_1&(( ,.~           1 2    )}) -: (_1 ;  0  1 _2)&setdiag) ai54
(_1&(( ,.~           1 2 3  )}) -: (_1 ;  0  1 __)&setdiag) ai54
(_1&(( ,.~      i.       3  )}) -: (_1 ;  0 _2   )&setdiag) ai54
(                               -: (_1 ;  0 _2  0)&setdiag) ai54
(_1&(( ,.~           1 2    )}) -: (_1 ;  0 _2  2)&setdiag) ai54
(_1&(( ,.~      i.       3  )}) -: (_1 ;  0 _2  _)&setdiag) ai54
(_1&(( ,.~           1 2    )}) -: (_1 ;  0 _2 _2)&setdiag) ai54
(_1&(( ,.~      i.       3  )}) -: (_1 ;  0 _2 __)&setdiag) ai54
(_1&(((,.~ >: ) i.         4)}) -: (_1 ; _1      )&setdiag) ai54
(_1&(((,.~ >: ) i.         4)}) -: (_1 ; _1  0   )&setdiag) ai54
(                               -: (_1 ; _1  0  0)&setdiag) ai54
(_1&(((,.~ >: )    0 1      )}) -: (_1 ; _1  0  2)&setdiag) ai54
(_1&(((,.~ >: ) i.         4)}) -: (_1 ; _1  0  _)&setdiag) ai54
(_1&(((,.~ >: )    0 1      )}) -: (_1 ; _1  0 _2)&setdiag) ai54
(_1&(((,.~ >: ) i.         4)}) -: (_1 ; _1  0 __)&setdiag) ai54
(_1&(((,.~ >: ) i.         4)}) -: (_1 ; _1 _1   )&setdiag) ai54
(                               -: (_1 ; _1 _1  0)&setdiag) ai54
(_1&(((,.~ >: )        2 3  )}) -: (_1 ; _1 _1  2)&setdiag) ai54
(_1&(((,.~ >: ) i.         4)}) -: (_1 ; _1 _1  _)&setdiag) ai54
(_1&(((,.~ >: )        2 3  )}) -: (_1 ; _1 _1 _2)&setdiag) ai54
(_1&(((,.~ >: ) i.         4)}) -: (_1 ; _1 _1 __)&setdiag) ai54
(_1&(((,.~ >: )      1 2 3  )}) -: (_1 ; _1  1   )&setdiag) ai54
(                               -: (_1 ; _1  1  0)&setdiag) ai54
(_1&(((,.~ >: )      1 2    )}) -: (_1 ; _1  1  2)&setdiag) ai54
(_1&(((,.~ >: )      1 2 3  )}) -: (_1 ; _1  1  _)&setdiag) ai54
(_1&(((,.~ >: )      1 2    )}) -: (_1 ; _1  1 _2)&setdiag) ai54
(_1&(((,.~ >: )      1 2 3  )}) -: (_1 ; _1  1 __)&setdiag) ai54
(_1&(((,.~ >: ) i.       3  )}) -: (_1 ; _1 _2   )&setdiag) ai54
(                               -: (_1 ; _1 _2  0)&setdiag) ai54
(_1&(((,.~ >: )      1 2    )}) -: (_1 ; _1 _2  2)&setdiag) ai54
(_1&(((,.~ >: ) i.       3  )}) -: (_1 ; _1 _2  _)&setdiag) ai54
(_1&(((,.~ >: )      1 2    )}) -: (_1 ; _1 _2 _2)&setdiag) ai54
(_1&(((,.~ >: ) i.       3  )}) -: (_1 ; _1 _2 __)&setdiag) ai54
(_1&(((,.  >: ) i.       3  )}) -: (_1 ;  1      )&setdiag) ai54
(_1&(((,.  >: ) i.       3  )}) -: (_1 ;  1  0   )&setdiag) ai54
(                               -: (_1 ;  1  0  0)&setdiag) ai54
(_1&(((,.  >: )    0 1      )}) -: (_1 ;  1  0  2)&setdiag) ai54
(_1&(((,.  >: ) i.       3  )}) -: (_1 ;  1  0  _)&setdiag) ai54
(_1&(((,.  >: )    0 1      )}) -: (_1 ;  1  0 _2)&setdiag) ai54
(_1&(((,.  >: ) i.       3  )}) -: (_1 ;  1  0 __)&setdiag) ai54
(_1&(((,.  >: ) i.       3  )}) -: (_1 ;  1 _1   )&setdiag) ai54
(                               -: (_1 ;  1 _1  0)&setdiag) ai54
(_1&(((,.  >: )      1 2    )}) -: (_1 ;  1 _1  2)&setdiag) ai54
(_1&(((,.  >: ) i.       3  )}) -: (_1 ;  1 _1  _)&setdiag) ai54
(_1&(((,.  >: )      1 2    )}) -: (_1 ;  1 _1 _2)&setdiag) ai54
(_1&(((,.  >: ) i.       3  )}) -: (_1 ;  1 _1 __)&setdiag) ai54
(_1&(((,.  >: )      1 2    )}) -: (_1 ;  1  1   )&setdiag) ai54
(                               -: (_1 ;  1  1  0)&setdiag) ai54
(_1&(((,.  >: )      1 2    )}) -: (_1 ;  1  1  2)&setdiag) ai54
(_1&(((,.  >: )      1 2    )}) -: (_1 ;  1  1  _)&setdiag) ai54
(_1&(((,.  >: )      1 2    )}) -: (_1 ;  1  1 _2)&setdiag) ai54
(_1&(((,.  >: )      1 2    )}) -: (_1 ;  1  1 __)&setdiag) ai54
(_1&(((,.  >: )    0 1      )}) -: (_1 ;  1 _2   )&setdiag) ai54
(                               -: (_1 ;  1 _2  0)&setdiag) ai54
(_1&(((,.  >: )    0 1      )}) -: (_1 ;  1 _2  2)&setdiag) ai54
(_1&(((,.  >: )    0 1      )}) -: (_1 ;  1 _2  _)&setdiag) ai54
(_1&(((,.  >: )    0 1      )}) -: (_1 ;  1 _2 _2)&setdiag) ai54
(_1&(((,.  >: )    0 1      )}) -: (_1 ;  1 _2 __)&setdiag) ai54
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2      )&setdiag) ai54
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2  0   )&setdiag) ai54
(                               -: (_1 ;  2  0  0)&setdiag) ai54
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2  0  2)&setdiag) ai54
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2  0  _)&setdiag) ai54
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2  0 _2)&setdiag) ai54
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2  0 __)&setdiag) ai54
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2 _1   )&setdiag) ai54
(                               -: (_1 ;  2 _1  0)&setdiag) ai54
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2 _1  2)&setdiag) ai54
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2 _1  _)&setdiag) ai54
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2 _1 _2)&setdiag) ai54
(_1&(((,.  2&+)    0 1      )}) -: (_1 ;  2 _1 __)&setdiag) ai54
(_1&((< 1 3                 )}) -: (_1 ;  2  1   )&setdiag) ai54
(                               -: (_1 ;  2  1  0)&setdiag) ai54
(_1&((< 1 3                 )}) -: (_1 ;  2  1  2)&setdiag) ai54
(_1&((< 1 3                 )}) -: (_1 ;  2  1  _)&setdiag) ai54
(_1&((< 1 3                 )}) -: (_1 ;  2  1 _2)&setdiag) ai54
(_1&((< 1 3                 )}) -: (_1 ;  2  1 __)&setdiag) ai54
(_1&((< 0 2                 )}) -: (_1 ;  2 _2   )&setdiag) ai54
(                               -: (_1 ;  2 _2  0)&setdiag) ai54
(_1&((< 0 2                 )}) -: (_1 ;  2 _2  2)&setdiag) ai54
(_1&((< 0 2                 )}) -: (_1 ;  2 _2  _)&setdiag) ai54
(_1&((< 0 2                 )}) -: (_1 ;  2 _2 _2)&setdiag) ai54
(_1&((< 0 2                 )}) -: (_1 ;  2 _2 __)&setdiag) ai54
(_1&(((,.~ 2&+) i.       3  )}) -: (_1 ; _2      )&setdiag) ai54
(_1&(((,.~ 2&+) i.       3  )}) -: (_1 ; _2  0   )&setdiag) ai54
(                               -: (_1 ; _2  0  0)&setdiag) ai54
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2  0  2)&setdiag) ai54
(_1&(((,.~ 2&+) i.       3  )}) -: (_1 ; _2  0  _)&setdiag) ai54
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2  0 _2)&setdiag) ai54
(_1&(((,.~ 2&+) i.       3  )}) -: (_1 ; _2  0 __)&setdiag) ai54
(_1&(((,.~ 2&+) i.       3  )}) -: (_1 ; _2 _1   )&setdiag) ai54
(                               -: (_1 ; _2 _1  0)&setdiag) ai54
(_1&(((,.~ 2&+)      1 2    )}) -: (_1 ; _2 _1  2)&setdiag) ai54
(_1&(((,.~ 2&+) i.       3  )}) -: (_1 ; _2 _1  _)&setdiag) ai54
(_1&(((,.~ 2&+)      1 2    )}) -: (_1 ; _2 _1 _2)&setdiag) ai54
(_1&(((,.~ 2&+) i.       3  )}) -: (_1 ; _2 _1 __)&setdiag) ai54
(_1&(((,.~ 2&+)      1 2    )}) -: (_1 ; _2  1   )&setdiag) ai54
(                               -: (_1 ; _2  1  0)&setdiag) ai54
(_1&(((,.~ 2&+)      1 2    )}) -: (_1 ; _2  1  2)&setdiag) ai54
(_1&(((,.~ 2&+)      1 2    )}) -: (_1 ; _2  1  _)&setdiag) ai54
(_1&(((,.~ 2&+)      1 2    )}) -: (_1 ; _2  1 _2)&setdiag) ai54
(_1&(((,.~ 2&+)      1 2    )}) -: (_1 ; _2  1 __)&setdiag) ai54
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2 _2   )&setdiag) ai54
(                               -: (_1 ; _2 _2  0)&setdiag) ai54
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2 _2  2)&setdiag) ai54
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2 _2  _)&setdiag) ai54
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2 _2 _2)&setdiag) ai54
(_1&(((,.~ 2&+)    0 1      )}) -: (_1 ; _2 _2 __)&setdiag) ai54

(_1&(( ,.~      i.           5)}) -: (_1 ; ''      )&setdiag) ai55
(_1&(( ,.~      i.           5)}) -: (_1 ;  0      )&setdiag) ai55
(_1&(( ,.~      i.           5)}) -: (_1 ;  0  0   )&setdiag) ai55
(                                 -: (_1 ;  0  0  0)&setdiag) ai55
(_1&(( ,.~         0 1        )}) -: (_1 ;  0  0  2)&setdiag) ai55
(_1&(( ,.~      i.           5)}) -: (_1 ;  0  0  _)&setdiag) ai55
(_1&(( ,.~         0 1        )}) -: (_1 ;  0  0 _2)&setdiag) ai55
(_1&(( ,.~      i.           5)}) -: (_1 ;  0  0 __)&setdiag) ai55
(_1&(( ,.~      i.           5)}) -: (_1 ;  0 _1   )&setdiag) ai55
(                                 -: (_1 ;  0 _1  0)&setdiag) ai55
(_1&(( ,.~               3 4  )}) -: (_1 ;  0 _1  2)&setdiag) ai55
(_1&(( ,.~      i.           5)}) -: (_1 ;  0 _1  _)&setdiag) ai55
(_1&(( ,.~               3 4  )}) -: (_1 ;  0 _1 _2)&setdiag) ai55
(_1&(( ,.~      i.           5)}) -: (_1 ;  0 _1 __)&setdiag) ai55
(_1&(( ,.~           1 2 3 4  )}) -: (_1 ;  0  1   )&setdiag) ai55
(                                 -: (_1 ;  0  1  0)&setdiag) ai55
(_1&(( ,.~           1 2      )}) -: (_1 ;  0  1  2)&setdiag) ai55
(_1&(( ,.~           1 2 3 4  )}) -: (_1 ;  0  1  _)&setdiag) ai55
(_1&(( ,.~           1 2      )}) -: (_1 ;  0  1 _2)&setdiag) ai55
(_1&(( ,.~           1 2 3 4  )}) -: (_1 ;  0  1 __)&setdiag) ai55
(_1&(( ,.~      i.         4  )}) -: (_1 ;  0 _2   )&setdiag) ai55
(                                 -: (_1 ;  0 _2  0)&setdiag) ai55
(_1&(( ,.~             2 3    )}) -: (_1 ;  0 _2  2)&setdiag) ai55
(_1&(( ,.~      i.         4  )}) -: (_1 ;  0 _2  _)&setdiag) ai55
(_1&(( ,.~             2 3    )}) -: (_1 ;  0 _2 _2)&setdiag) ai55
(_1&(( ,.~      i.         4  )}) -: (_1 ;  0 _2 __)&setdiag) ai55
(_1&(((,.~ >: ) i.         4  )}) -: (_1 ; _1      )&setdiag) ai55
(_1&(((,.~ >: ) i.         4  )}) -: (_1 ; _1  0   )&setdiag) ai55
(                                 -: (_1 ; _1  0  0)&setdiag) ai55
(_1&(((,.~ >: )    0 1        )}) -: (_1 ; _1  0  2)&setdiag) ai55
(_1&(((,.~ >: ) i.         4  )}) -: (_1 ; _1  0  _)&setdiag) ai55
(_1&(((,.~ >: )    0 1        )}) -: (_1 ; _1  0 _2)&setdiag) ai55
(_1&(((,.~ >: ) i.         4  )}) -: (_1 ; _1  0 __)&setdiag) ai55
(_1&(((,.~ >: ) i.         4  )}) -: (_1 ; _1 _1   )&setdiag) ai55
(                                 -: (_1 ; _1 _1  0)&setdiag) ai55
(_1&(((,.~ >: )        2 3    )}) -: (_1 ; _1 _1  2)&setdiag) ai55
(_1&(((,.~ >: ) i.         4  )}) -: (_1 ; _1 _1  _)&setdiag) ai55
(_1&(((,.~ >: )        2 3    )}) -: (_1 ; _1 _1 _2)&setdiag) ai55
(_1&(((,.~ >: ) i.         4  )}) -: (_1 ; _1 _1 __)&setdiag) ai55
(_1&(((,.~ >: )      1 2 3    )}) -: (_1 ; _1  1   )&setdiag) ai55
(                                 -: (_1 ; _1  1  0)&setdiag) ai55
(_1&(((,.~ >: )      1 2      )}) -: (_1 ; _1  1  2)&setdiag) ai55
(_1&(((,.~ >: )      1 2 3    )}) -: (_1 ; _1  1  _)&setdiag) ai55
(_1&(((,.~ >: )      1 2      )}) -: (_1 ; _1  1 _2)&setdiag) ai55
(_1&(((,.~ >: )      1 2 3    )}) -: (_1 ; _1  1 __)&setdiag) ai55
(_1&(((,.~ >: ) i.       3    )}) -: (_1 ; _1 _2   )&setdiag) ai55
(                                 -: (_1 ; _1 _2  0)&setdiag) ai55
(_1&(((,.~ >: )      1 2      )}) -: (_1 ; _1 _2  2)&setdiag) ai55
(_1&(((,.~ >: ) i.       3    )}) -: (_1 ; _1 _2  _)&setdiag) ai55
(_1&(((,.~ >: )      1 2      )}) -: (_1 ; _1 _2 _2)&setdiag) ai55
(_1&(((,.~ >: ) i.       3    )}) -: (_1 ; _1 _2 __)&setdiag) ai55
(_1&(((,.  >: ) i.         4  )}) -: (_1 ;  1      )&setdiag) ai55
(_1&(((,.  >: ) i.         4  )}) -: (_1 ;  1  0   )&setdiag) ai55
(                                 -: (_1 ;  1  0  0)&setdiag) ai55
(_1&(((,.  >: )    0 1        )}) -: (_1 ;  1  0  2)&setdiag) ai55
(_1&(((,.  >: ) i.         4  )}) -: (_1 ;  1  0  _)&setdiag) ai55
(_1&(((,.  >: )    0 1        )}) -: (_1 ;  1  0 _2)&setdiag) ai55
(_1&(((,.  >: ) i.         4  )}) -: (_1 ;  1  0 __)&setdiag) ai55
(_1&(((,.  >: ) i.         4  )}) -: (_1 ;  1 _1   )&setdiag) ai55
(                                 -: (_1 ;  1 _1  0)&setdiag) ai55
(_1&(((,.  >: )        2 3    )}) -: (_1 ;  1 _1  2)&setdiag) ai55
(_1&(((,.  >: ) i.         4  )}) -: (_1 ;  1 _1  _)&setdiag) ai55
(_1&(((,.  >: )        2 3    )}) -: (_1 ;  1 _1 _2)&setdiag) ai55
(_1&(((,.  >: ) i.         4  )}) -: (_1 ;  1 _1 __)&setdiag) ai55
(_1&(((,.  >: )      1 2 3    )}) -: (_1 ;  1  1   )&setdiag) ai55
(                                 -: (_1 ;  1  1  0)&setdiag) ai55
(_1&(((,.  >: )      1 2      )}) -: (_1 ;  1  1  2)&setdiag) ai55
(_1&(((,.  >: )      1 2 3    )}) -: (_1 ;  1  1  _)&setdiag) ai55
(_1&(((,.  >: )      1 2      )}) -: (_1 ;  1  1 _2)&setdiag) ai55
(_1&(((,.  >: )      1 2 3    )}) -: (_1 ;  1  1 __)&setdiag) ai55
(_1&(((,.  >: ) i.       3    )}) -: (_1 ;  1 _2   )&setdiag) ai55
(                                 -: (_1 ;  1 _2  0)&setdiag) ai55
(_1&(((,.  >: )      1 2      )}) -: (_1 ;  1 _2  2)&setdiag) ai55
(_1&(((,.  >: ) i.       3    )}) -: (_1 ;  1 _2  _)&setdiag) ai55
(_1&(((,.  >: )      1 2      )}) -: (_1 ;  1 _2 _2)&setdiag) ai55
(_1&(((,.  >: ) i.       3    )}) -: (_1 ;  1 _2 __)&setdiag) ai55
(_1&(((,.  2&+) i.       3    )}) -: (_1 ;  2      )&setdiag) ai55
(_1&(((,.  2&+) i.       3    )}) -: (_1 ;  2  0   )&setdiag) ai55
(                                 -: (_1 ;  2  0  0)&setdiag) ai55
(_1&(((,.  2&+)    0 1        )}) -: (_1 ;  2  0  2)&setdiag) ai55
(_1&(((,.  2&+) i.       3    )}) -: (_1 ;  2  0  _)&setdiag) ai55
(_1&(((,.  2&+)    0 1        )}) -: (_1 ;  2  0 _2)&setdiag) ai55
(_1&(((,.  2&+) i.       3    )}) -: (_1 ;  2  0 __)&setdiag) ai55
(_1&(((,.  2&+) i.       3    )}) -: (_1 ;  2 _1   )&setdiag) ai55
(                                 -: (_1 ;  2 _1  0)&setdiag) ai55
(_1&(((,.  2&+)      1 2      )}) -: (_1 ;  2 _1  2)&setdiag) ai55
(_1&(((,.  2&+) i.       3    )}) -: (_1 ;  2 _1  _)&setdiag) ai55
(_1&(((,.  2&+)      1 2      )}) -: (_1 ;  2 _1 _2)&setdiag) ai55
(_1&(((,.  2&+) i.       3    )}) -: (_1 ;  2 _1 __)&setdiag) ai55
(_1&(((,.  2&+)      1 2      )}) -: (_1 ;  2  1   )&setdiag) ai55
(                                 -: (_1 ;  2  1  0)&setdiag) ai55
(_1&(((,.  2&+)      1 2      )}) -: (_1 ;  2  1  2)&setdiag) ai55
(_1&(((,.  2&+)      1 2      )}) -: (_1 ;  2  1  _)&setdiag) ai55
(_1&(((,.  2&+)      1 2      )}) -: (_1 ;  2  1 _2)&setdiag) ai55
(_1&(((,.  2&+)      1 2      )}) -: (_1 ;  2  1 __)&setdiag) ai55
(_1&(((,.  2&+)    0 1        )}) -: (_1 ;  2 _2   )&setdiag) ai55
(                                 -: (_1 ;  2 _2  0)&setdiag) ai55
(_1&(((,.  2&+)    0 1        )}) -: (_1 ;  2 _2  2)&setdiag) ai55
(_1&(((,.  2&+)    0 1        )}) -: (_1 ;  2 _2  _)&setdiag) ai55
(_1&(((,.  2&+)    0 1        )}) -: (_1 ;  2 _2 _2)&setdiag) ai55
(_1&(((,.  2&+)    0 1        )}) -: (_1 ;  2 _2 __)&setdiag) ai55
(_1&(((,.~ 2&+) i.       3    )}) -: (_1 ; _2      )&setdiag) ai55
(_1&(((,.~ 2&+) i.       3    )}) -: (_1 ; _2  0   )&setdiag) ai55
(                                 -: (_1 ; _2  0  0)&setdiag) ai55
(_1&(((,.~ 2&+)    0 1        )}) -: (_1 ; _2  0  2)&setdiag) ai55
(_1&(((,.~ 2&+) i.       3    )}) -: (_1 ; _2  0  _)&setdiag) ai55
(_1&(((,.~ 2&+)    0 1        )}) -: (_1 ; _2  0 _2)&setdiag) ai55
(_1&(((,.~ 2&+) i.       3    )}) -: (_1 ; _2  0 __)&setdiag) ai55
(_1&(((,.~ 2&+) i.       3    )}) -: (_1 ; _2 _1   )&setdiag) ai55
(                                 -: (_1 ; _2 _1  0)&setdiag) ai55
(_1&(((,.~ 2&+)      1 2      )}) -: (_1 ; _2 _1  2)&setdiag) ai55
(_1&(((,.~ 2&+) i.       3    )}) -: (_1 ; _2 _1  _)&setdiag) ai55
(_1&(((,.~ 2&+)      1 2      )}) -: (_1 ; _2 _1 _2)&setdiag) ai55
(_1&(((,.~ 2&+) i.       3    )}) -: (_1 ; _2 _1 __)&setdiag) ai55
(_1&(((,.~ 2&+)      1 2      )}) -: (_1 ; _2  1   )&setdiag) ai55
(                                 -: (_1 ; _2  1  0)&setdiag) ai55
(_1&(((,.~ 2&+)      1 2      )}) -: (_1 ; _2  1  2)&setdiag) ai55
(_1&(((,.~ 2&+)      1 2      )}) -: (_1 ; _2  1  _)&setdiag) ai55
(_1&(((,.~ 2&+)      1 2      )}) -: (_1 ; _2  1 _2)&setdiag) ai55
(_1&(((,.~ 2&+)      1 2      )}) -: (_1 ; _2  1 __)&setdiag) ai55
(_1&(((,.~ 2&+)    0 1        )}) -: (_1 ; _2 _2   )&setdiag) ai55
(                                 -: (_1 ; _2 _2  0)&setdiag) ai55
(_1&(((,.~ 2&+)    0 1        )}) -: (_1 ; _2 _2  2)&setdiag) ai55
(_1&(((,.~ 2&+)    0 1        )}) -: (_1 ; _2 _2  _)&setdiag) ai55
(_1&(((,.~ 2&+)    0 1        )}) -: (_1 ; _2 _2 _2)&setdiag) ai55
(_1&(((,.~ 2&+)    0 1        )}) -: (_1 ; _2 _2 __)&setdiag) ai55

NB. vector e

(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ; ''      )&setdiag) ai45
(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0      )&setdiag) ai45
(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0  0   )&setdiag) ai45
(                                        -: (''          ;  0  0  0)&setdiag) ai45
(_1 _2      &(( ,.~         0 1      )}) -: (_1 _2       ;  0  0  2)&setdiag) ai45
(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0  0  _)&setdiag) ai45
(_2 _1      &(( ,.~         0 1      )}) -: (_1 _2       ;  0  0 _2)&setdiag) ai45
(_4 _3 _2 _1&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0  0 __)&setdiag) ai45
(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0 _1   )&setdiag) ai45
(                                        -: (''          ;  0 _1  0)&setdiag) ai45
(_1 _2      &(( ,.~             2 3  )}) -: (_1 _2       ;  0 _1  2)&setdiag) ai45
(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0 _1  _)&setdiag) ai45
(_2 _1      &(( ,.~             2 3  )}) -: (_1 _2       ;  0 _1 _2)&setdiag) ai45
(_4 _3 _2 _1&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0 _1 __)&setdiag) ai45
(_1 _2 _3   &(( ,.~           1 2 3  )}) -: (_1 _2 _3    ;  0  1   )&setdiag) ai45
(                                        -: (''          ;  0  1  0)&setdiag) ai45
(_1 _2      &(( ,.~           1 2    )}) -: (_1 _2       ;  0  1  2)&setdiag) ai45
(_1 _2 _3   &(( ,.~           1 2 3  )}) -: (_1 _2 _3    ;  0  1  _)&setdiag) ai45
(_2 _1      &(( ,.~           1 2    )}) -: (_1 _2       ;  0  1 _2)&setdiag) ai45
(_3 _2 _1   &(( ,.~           1 2 3  )}) -: (_1 _2 _3    ;  0  1 __)&setdiag) ai45
(_1 _2 _3   &(( ,.~      i.       3  )}) -: (_1 _2 _3    ;  0 _2   )&setdiag) ai45
(                                        -: (''          ;  0 _2  0)&setdiag) ai45
(_1 _2      &(( ,.~           1 2    )}) -: (_1 _2       ;  0 _2  2)&setdiag) ai45
(_1 _2 _3   &(( ,.~      i.       3  )}) -: (_1 _2 _3    ;  0 _2  _)&setdiag) ai45
(_2 _1      &(( ,.~           1 2    )}) -: (_1 _2       ;  0 _2 _2)&setdiag) ai45
(_3 _2 _1   &(( ,.~      i.       3  )}) -: (_1 _2 _3    ;  0 _2 __)&setdiag) ai45
(_1 _2 _3   &(((,.~ >: ) i.       3  )}) -: (_1 _2 _3    ; _1      )&setdiag) ai45
(_1 _2 _3   &(((,.~ >: ) i.       3  )}) -: (_1 _2 _3    ; _1  0   )&setdiag) ai45
(                                        -: (''          ; _1  0  0)&setdiag) ai45
(_1 _2      &(((,.~ >: )    0 1      )}) -: (_1 _2       ; _1  0  2)&setdiag) ai45
(_1 _2 _3   &(((,.~ >: ) i.       3  )}) -: (_1 _2 _3    ; _1  0  _)&setdiag) ai45
(_2 _1      &(((,.~ >: )    0 1      )}) -: (_1 _2       ; _1  0 _2)&setdiag) ai45
(_3 _2 _1   &(((,.~ >: ) i.       3  )}) -: (_1 _2 _3    ; _1  0 __)&setdiag) ai45
(_1 _2 _3   &(((,.~ >: ) i.       3  )}) -: (_1 _2 _3    ; _1 _1   )&setdiag) ai45
(                                        -: (''          ; _1 _1  0)&setdiag) ai45
(_1 _2      &(((,.~ >: )      1 2    )}) -: (_1 _2       ; _1 _1  2)&setdiag) ai45
(_1 _2 _3   &(((,.~ >: ) i.       3  )}) -: (_1 _2 _3    ; _1 _1  _)&setdiag) ai45
(_2 _1      &(((,.~ >: )      1 2    )}) -: (_1 _2       ; _1 _1 _2)&setdiag) ai45
(_3 _2 _1   &(((,.~ >: ) i.       3  )}) -: (_1 _2 _3    ; _1 _1 __)&setdiag) ai45
(_1 _2      &(((,.~ >: )      1 2    )}) -: (_1 _2       ; _1  1   )&setdiag) ai45
(                                        -: (''          ; _1  1  0)&setdiag) ai45
(_1 _2      &(((,.~ >: )      1 2    )}) -: (_1 _2       ; _1  1  2)&setdiag) ai45
(_1 _2      &(((,.~ >: )      1 2    )}) -: (_1 _2       ; _1  1  _)&setdiag) ai45
(_2 _1      &(((,.~ >: )      1 2    )}) -: (_1 _2       ; _1  1 _2)&setdiag) ai45
(_2 _1      &(((,.~ >: )      1 2    )}) -: (_1 _2       ; _1  1 __)&setdiag) ai45
(_1 _2      &(((,.~ >: )    0 1      )}) -: (_1 _2       ; _1 _2   )&setdiag) ai45
(                                        -: (''          ; _1 _2  0)&setdiag) ai45
(_1 _2      &(((,.~ >: )    0 1      )}) -: (_1 _2       ; _1 _2  2)&setdiag) ai45
(_1 _2      &(((,.~ >: )    0 1      )}) -: (_1 _2       ; _1 _2  _)&setdiag) ai45
(_2 _1      &(((,.~ >: )    0 1      )}) -: (_1 _2       ; _1 _2 _2)&setdiag) ai45
(_2 _1      &(((,.~ >: )    0 1      )}) -: (_1 _2       ; _1 _2 __)&setdiag) ai45
(_1 _2 _3 _4&(((,.  >: ) i.         4)}) -: (_1 _2 _3 _4 ;  1      )&setdiag) ai45
(_1 _2 _3 _4&(((,.  >: ) i.         4)}) -: (_1 _2 _3 _4 ;  1  0   )&setdiag) ai45
(                                        -: (''          ;  1  0  0)&setdiag) ai45
(_1 _2      &(((,.  >: )    0 1      )}) -: (_1 _2       ;  1  0  2)&setdiag) ai45
(_1 _2 _3 _4&(((,.  >: ) i.         4)}) -: (_1 _2 _3 _4 ;  1  0  _)&setdiag) ai45
(_2 _1      &(((,.  >: )    0 1      )}) -: (_1 _2       ;  1  0 _2)&setdiag) ai45
(_4 _3 _2 _1&(((,.  >: ) i.         4)}) -: (_1 _2 _3 _4 ;  1  0 __)&setdiag) ai45
(_1 _2 _3 _4&(((,.  >: ) i.         4)}) -: (_1 _2 _3 _4 ;  1 _1   )&setdiag) ai45
(                                        -: (''          ;  1 _1  0)&setdiag) ai45
(_1 _2      &(((,.  >: )        2 3  )}) -: (_1 _2       ;  1 _1  2)&setdiag) ai45
(_1 _2 _3 _4&(((,.  >: ) i.         4)}) -: (_1 _2 _3 _4 ;  1 _1  _)&setdiag) ai45
(_2 _1      &(((,.  >: )        2 3  )}) -: (_1 _2       ;  1 _1 _2)&setdiag) ai45
(_4 _3 _2 _1&(((,.  >: ) i.         4)}) -: (_1 _2 _3 _4 ;  1 _1 __)&setdiag) ai45
(_1 _2 _3   &(((,.  >: )      1 2 3  )}) -: (_1 _2 _3    ;  1  1   )&setdiag) ai45
(                                        -: (''          ;  1  1  0)&setdiag) ai45
(_1 _2      &(((,.  >: )      1 2    )}) -: (_1 _2       ;  1  1  2)&setdiag) ai45
(_1 _2 _3   &(((,.  >: )      1 2 3  )}) -: (_1 _2 _3    ;  1  1  _)&setdiag) ai45
(_2 _1      &(((,.  >: )      1 2    )}) -: (_1 _2       ;  1  1 _2)&setdiag) ai45
(_3 _2 _1   &(((,.  >: )      1 2 3  )}) -: (_1 _2 _3    ;  1  1 __)&setdiag) ai45
(_1 _2 _3   &(((,.  >: ) i.       3  )}) -: (_1 _2 _3    ;  1 _2   )&setdiag) ai45
(                                        -: (''          ;  1 _2  0)&setdiag) ai45
(_1 _2      &(((,.  >: )      1 2    )}) -: (_1 _2       ;  1 _2  2)&setdiag) ai45
(_1 _2 _3   &(((,.  >: ) i.       3  )}) -: (_1 _2 _3    ;  1 _2  _)&setdiag) ai45
(_2 _1      &(((,.  >: )      1 2    )}) -: (_1 _2       ;  1 _2 _2)&setdiag) ai45
(_3 _2 _1   &(((,.  >: ) i.       3  )}) -: (_1 _2 _3    ;  1 _2 __)&setdiag) ai45
(_1 _2 _3   &(((,.  2&+) i.       3  )}) -: (_1 _2 _3    ;  2      )&setdiag) ai45
(_1 _2 _3   &(((,.  2&+) i.       3  )}) -: (_1 _2 _3    ;  2  0   )&setdiag) ai45
(                                        -: (''          ;  2  0  0)&setdiag) ai45
(_1 _2      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2  0  2)&setdiag) ai45
(_1 _2 _3   &(((,.  2&+) i.       3  )}) -: (_1 _2 _3    ;  2  0  _)&setdiag) ai45
(_2 _1      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2  0 _2)&setdiag) ai45
(_3 _2 _1   &(((,.  2&+) i.       3  )}) -: (_1 _2 _3    ;  2  0 __)&setdiag) ai45
(_1 _2 _3   &(((,.  2&+) i.       3  )}) -: (_1 _2 _3    ;  2 _1   )&setdiag) ai45
(                                        -: (''          ;  2 _1  0)&setdiag) ai45
(_1 _2      &(((,.  2&+)      1 2    )}) -: (_1 _2       ;  2 _1  2)&setdiag) ai45
(_1 _2 _3   &(((,.  2&+) i.       3  )}) -: (_1 _2 _3    ;  2 _1  _)&setdiag) ai45
(_2 _1      &(((,.  2&+)      1 2    )}) -: (_1 _2       ;  2 _1 _2)&setdiag) ai45
(_3 _2 _1   &(((,.  2&+) i.       3  )}) -: (_1 _2 _3    ;  2 _1 __)&setdiag) ai45
(_1 _2      &(((,.  2&+)      1 2    )}) -: (_1 _2       ;  2  1   )&setdiag) ai45
(                                        -: (''          ;  2  1  0)&setdiag) ai45
(_1 _2      &(((,.  2&+)      1 2    )}) -: (_1 _2       ;  2  1  2)&setdiag) ai45
(_1 _2      &(((,.  2&+)      1 2    )}) -: (_1 _2       ;  2  1  _)&setdiag) ai45
(_2 _1      &(((,.  2&+)      1 2    )}) -: (_1 _2       ;  2  1 _2)&setdiag) ai45
(_2 _1      &(((,.  2&+)      1 2    )}) -: (_1 _2       ;  2  1 __)&setdiag) ai45
(_1 _2      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2 _2   )&setdiag) ai45
(                                        -: (''          ;  2 _2  0)&setdiag) ai45
(_1 _2      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2 _2  2)&setdiag) ai45
(_1 _2      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2 _2  _)&setdiag) ai45
(_2 _1      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2 _2 _2)&setdiag) ai45
(_2 _1      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2 _2 __)&setdiag) ai45
(_1 _2      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2      )&setdiag) ai45
(_1 _2      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2  0   )&setdiag) ai45
(                                        -: (''          ; _2  0  0)&setdiag) ai45
(_1 _2      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2  0  2)&setdiag) ai45
(_1 _2      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2  0  _)&setdiag) ai45
(_2 _1      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2  0 _2)&setdiag) ai45
(_2 _1      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2  0 __)&setdiag) ai45
(_1 _2      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2 _1   )&setdiag) ai45
(                                        -: (''          ; _2 _1  0)&setdiag) ai45
(_1 _2      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2 _1  2)&setdiag) ai45
(_1 _2      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2 _1  _)&setdiag) ai45
(_2 _1      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2 _1 _2)&setdiag) ai45
(_2 _1      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2 _1 __)&setdiag) ai45
(_1         &((< 3 1                 )}) -: ((, _1)      ; _2  1   )&setdiag) ai45
(                                        -: (''          ; _2  1  0)&setdiag) ai45
(_1         &((< 3 1                 )}) -: ((, _1)      ; _2  1  2)&setdiag) ai45
(_1         &((< 3 1                 )}) -: ((, _1)      ; _2  1  _)&setdiag) ai45
(_1         &((< 3 1                 )}) -: ((, _1)      ; _2  1 _2)&setdiag) ai45
(_1         &((< 3 1                 )}) -: ((, _1)      ; _2  1 __)&setdiag) ai45
(_1         &((< 2 0                 )}) -: ((, _1)      ; _2 _2   )&setdiag) ai45
(                                        -: (''          ; _2 _2  0)&setdiag) ai45
(_1         &((< 2 0                 )}) -: ((, _1)      ; _2 _2  2)&setdiag) ai45
(_1         &((< 2 0                 )}) -: ((, _1)      ; _2 _2  _)&setdiag) ai45
(_1         &((< 2 0                 )}) -: ((, _1)      ; _2 _2 _2)&setdiag) ai45
(_1         &((< 2 0                 )}) -: ((, _1)      ; _2 _2 __)&setdiag) ai45

(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ; ''      )&setdiag) ai54
(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0      )&setdiag) ai54
(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0  0   )&setdiag) ai54
(                                        -: (''          ;  0  0  0)&setdiag) ai54
(_1 _2      &(( ,.~         0 1      )}) -: (_1 _2       ;  0  0  2)&setdiag) ai54
(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0  0  _)&setdiag) ai54
(_2 _1      &(( ,.~         0 1      )}) -: (_1 _2       ;  0  0 _2)&setdiag) ai54
(_4 _3 _2 _1&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0  0 __)&setdiag) ai54
(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0 _1   )&setdiag) ai54
(                                        -: (''          ;  0 _1  0)&setdiag) ai54
(_1 _2      &(( ,.~             2 3  )}) -: (_1 _2       ;  0 _1  2)&setdiag) ai54
(_1 _2 _3 _4&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0 _1  _)&setdiag) ai54
(_2 _1      &(( ,.~             2 3  )}) -: (_1 _2       ;  0 _1 _2)&setdiag) ai54
(_4 _3 _2 _1&(( ,.~      i.         4)}) -: (_1 _2 _3 _4 ;  0 _1 __)&setdiag) ai54
(_1 _2 _3   &(( ,.~           1 2 3  )}) -: (_1 _2 _3    ;  0  1   )&setdiag) ai54
(                                        -: (''          ;  0  1  0)&setdiag) ai54
(_1 _2      &(( ,.~           1 2    )}) -: (_1 _2       ;  0  1  2)&setdiag) ai54
(_1 _2 _3   &(( ,.~           1 2 3  )}) -: (_1 _2 _3    ;  0  1  _)&setdiag) ai54
(_2 _1      &(( ,.~           1 2    )}) -: (_1 _2       ;  0  1 _2)&setdiag) ai54
(_3 _2 _1   &(( ,.~           1 2 3  )}) -: (_1 _2 _3    ;  0  1 __)&setdiag) ai54
(_1 _2 _3   &(( ,.~      i.       3  )}) -: (_1 _2 _3    ;  0 _2   )&setdiag) ai54
(                                        -: (''          ;  0 _2  0)&setdiag) ai54
(_1 _2      &(( ,.~           1 2    )}) -: (_1 _2       ;  0 _2  2)&setdiag) ai54
(_1 _2 _3   &(( ,.~      i.       3  )}) -: (_1 _2 _3    ;  0 _2  _)&setdiag) ai54
(_2 _1      &(( ,.~           1 2    )}) -: (_1 _2       ;  0 _2 _2)&setdiag) ai54
(_3 _2 _1   &(( ,.~      i.       3  )}) -: (_1 _2 _3    ;  0 _2 __)&setdiag) ai54
(_1 _2 _3 _4&(((,.~ >: ) i.         4)}) -: (_1 _2 _3 _4 ; _1      )&setdiag) ai54
(_1 _2 _3 _4&(((,.~ >: ) i.         4)}) -: (_1 _2 _3 _4 ; _1  0   )&setdiag) ai54
(                                        -: (''          ; _1  0  0)&setdiag) ai54
(_1 _2      &(((,.~ >: )    0 1      )}) -: (_1 _2       ; _1  0  2)&setdiag) ai54
(_1 _2 _3 _4&(((,.~ >: ) i.         4)}) -: (_1 _2 _3 _4 ; _1  0  _)&setdiag) ai54
(_2 _1      &(((,.~ >: )    0 1      )}) -: (_1 _2       ; _1  0 _2)&setdiag) ai54
(_4 _3 _2 _1&(((,.~ >: ) i.         4)}) -: (_1 _2 _3 _4 ; _1  0 __)&setdiag) ai54
(_1 _2 _3 _4&(((,.~ >: ) i.         4)}) -: (_1 _2 _3 _4 ; _1 _1   )&setdiag) ai54
(                                        -: (''          ; _1 _1  0)&setdiag) ai54
(_1 _2      &(((,.~ >: )        2 3  )}) -: (_1 _2       ; _1 _1  2)&setdiag) ai54
(_1 _2 _3 _4&(((,.~ >: ) i.         4)}) -: (_1 _2 _3 _4 ; _1 _1  _)&setdiag) ai54
(_2 _1      &(((,.~ >: )        2 3  )}) -: (_1 _2       ; _1 _1 _2)&setdiag) ai54
(_4 _3 _2 _1&(((,.~ >: ) i.         4)}) -: (_1 _2 _3 _4 ; _1 _1 __)&setdiag) ai54
(_1 _2 _3   &(((,.~ >: )      1 2 3  )}) -: (_1 _2 _3    ; _1  1   )&setdiag) ai54
(                                        -: (''          ; _1  1  0)&setdiag) ai54
(_1 _2      &(((,.~ >: )      1 2    )}) -: (_1 _2       ; _1  1  2)&setdiag) ai54
(_1 _2 _3   &(((,.~ >: )      1 2 3  )}) -: (_1 _2 _3    ; _1  1  _)&setdiag) ai54
(_2 _1      &(((,.~ >: )      1 2    )}) -: (_1 _2       ; _1  1 _2)&setdiag) ai54
(_3 _2 _1   &(((,.~ >: )      1 2 3  )}) -: (_1 _2 _3    ; _1  1 __)&setdiag) ai54
(_1 _2 _3   &(((,.~ >: ) i.       3  )}) -: (_1 _2 _3    ; _1 _2   )&setdiag) ai54
(                                        -: (''          ; _1 _2  0)&setdiag) ai54
(_1 _2      &(((,.~ >: )      1 2    )}) -: (_1 _2       ; _1 _2  2)&setdiag) ai54
(_1 _2 _3   &(((,.~ >: ) i.       3  )}) -: (_1 _2 _3    ; _1 _2  _)&setdiag) ai54
(_2 _1      &(((,.~ >: )      1 2    )}) -: (_1 _2       ; _1 _2 _2)&setdiag) ai54
(_3 _2 _1   &(((,.~ >: ) i.       3  )}) -: (_1 _2 _3    ; _1 _2 __)&setdiag) ai54
(_1 _2 _3   &(((,.  >: ) i.       3  )}) -: (_1 _2 _3    ;  1      )&setdiag) ai54
(_1 _2 _3   &(((,.  >: ) i.       3  )}) -: (_1 _2 _3    ;  1  0   )&setdiag) ai54
(                                        -: (''          ;  1  0  0)&setdiag) ai54
(_1 _2      &(((,.  >: )    0 1      )}) -: (_1 _2       ;  1  0  2)&setdiag) ai54
(_1 _2 _3   &(((,.  >: ) i.       3  )}) -: (_1 _2 _3    ;  1  0  _)&setdiag) ai54
(_2 _1      &(((,.  >: )    0 1      )}) -: (_1 _2       ;  1  0 _2)&setdiag) ai54
(_3 _2 _1   &(((,.  >: ) i.       3  )}) -: (_1 _2 _3    ;  1  0 __)&setdiag) ai54
(_1 _2 _3   &(((,.  >: ) i.       3  )}) -: (_1 _2 _3    ;  1 _1   )&setdiag) ai54
(                                        -: (''          ;  1 _1  0)&setdiag) ai54
(_1 _2      &(((,.  >: )      1 2    )}) -: (_1 _2       ;  1 _1  2)&setdiag) ai54
(_1 _2 _3   &(((,.  >: ) i.       3  )}) -: (_1 _2 _3    ;  1 _1  _)&setdiag) ai54
(_2 _1      &(((,.  >: )      1 2    )}) -: (_1 _2       ;  1 _1 _2)&setdiag) ai54
(_3 _2 _1   &(((,.  >: ) i.       3  )}) -: (_1 _2 _3    ;  1 _1 __)&setdiag) ai54
(_1 _2      &(((,.  >: )      1 2    )}) -: (_1 _2       ;  1  1   )&setdiag) ai54
(                                        -: (''          ;  1  1  0)&setdiag) ai54
(_1 _2      &(((,.  >: )      1 2    )}) -: (_1 _2       ;  1  1  2)&setdiag) ai54
(_1 _2      &(((,.  >: )      1 2    )}) -: (_1 _2       ;  1  1  _)&setdiag) ai54
(_2 _1      &(((,.  >: )      1 2    )}) -: (_1 _2       ;  1  1 _2)&setdiag) ai54
(_2 _1      &(((,.  >: )      1 2    )}) -: (_1 _2       ;  1  1 __)&setdiag) ai54
(_1 _2      &(((,.  >: )    0 1      )}) -: (_1 _2       ;  1 _2   )&setdiag) ai54
(                                        -: (''          ;  1 _2  0)&setdiag) ai54
(_1 _2      &(((,.  >: )    0 1      )}) -: (_1 _2       ;  1 _2  2)&setdiag) ai54
(_1 _2      &(((,.  >: )    0 1      )}) -: (_1 _2       ;  1 _2  _)&setdiag) ai54
(_2 _1      &(((,.  >: )    0 1      )}) -: (_1 _2       ;  1 _2 _2)&setdiag) ai54
(_2 _1      &(((,.  >: )    0 1      )}) -: (_1 _2       ;  1 _2 __)&setdiag) ai54
(_1 _2      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2      )&setdiag) ai54
(_1 _2      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2  0   )&setdiag) ai54
(                                        -: (''          ;  2  0  0)&setdiag) ai54
(_1 _2      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2  0  2)&setdiag) ai54
(_1 _2      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2  0  _)&setdiag) ai54
(_2 _1      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2  0 _2)&setdiag) ai54
(_2 _1      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2  0 __)&setdiag) ai54
(_1 _2      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2 _1   )&setdiag) ai54
(                                        -: (''          ;  2 _1  0)&setdiag) ai54
(_1 _2      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2 _1  2)&setdiag) ai54
(_1 _2      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2 _1  _)&setdiag) ai54
(_2 _1      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2 _1 _2)&setdiag) ai54
(_2 _1      &(((,.  2&+)    0 1      )}) -: (_1 _2       ;  2 _1 __)&setdiag) ai54
(_1         &((< 1 3                 )}) -: ((, _1)      ;  2  1   )&setdiag) ai54
(                                        -: (''          ;  2  1  0)&setdiag) ai54
(_1         &((< 1 3                 )}) -: ((, _1)      ;  2  1  2)&setdiag) ai54
(_1         &((< 1 3                 )}) -: ((, _1)      ;  2  1  _)&setdiag) ai54
(_1         &((< 1 3                 )}) -: ((, _1)      ;  2  1 _2)&setdiag) ai54
(_1         &((< 1 3                 )}) -: ((, _1)      ;  2  1 __)&setdiag) ai54
(_1         &((< 0 2                 )}) -: ((, _1)      ;  2 _2   )&setdiag) ai54
(                                        -: (''          ;  2 _2  0)&setdiag) ai54
(_1         &((< 0 2                 )}) -: ((, _1)      ;  2 _2  2)&setdiag) ai54
(_1         &((< 0 2                 )}) -: ((, _1)      ;  2 _2  _)&setdiag) ai54
(_1         &((< 0 2                 )}) -: ((, _1)      ;  2 _2 _2)&setdiag) ai54
(_1         &((< 0 2                 )}) -: ((, _1)      ;  2 _2 __)&setdiag) ai54
(_1 _2 _3   &(((,.~ 2&+) i.       3  )}) -: (_1 _2 _3    ; _2      )&setdiag) ai54
(_1 _2 _3   &(((,.~ 2&+) i.       3  )}) -: (_1 _2 _3    ; _2  0   )&setdiag) ai54
(                                        -: (''          ; _2  0  0)&setdiag) ai54
(_1 _2      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2  0  2)&setdiag) ai54
(_1 _2 _3   &(((,.~ 2&+) i.       3  )}) -: (_1 _2 _3    ; _2  0  _)&setdiag) ai54
(_2 _1      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2  0 _2)&setdiag) ai54
(_3 _2 _1   &(((,.~ 2&+) i.       3  )}) -: (_1 _2 _3    ; _2  0 __)&setdiag) ai54
(_1 _2 _3   &(((,.~ 2&+) i.       3  )}) -: (_1 _2 _3    ; _2 _1   )&setdiag) ai54
(                                        -: (''          ; _2 _1  0)&setdiag) ai54
(_1 _2      &(((,.~ 2&+)      1 2    )}) -: (_1 _2       ; _2 _1  2)&setdiag) ai54
(_1 _2 _3   &(((,.~ 2&+) i.       3  )}) -: (_1 _2 _3    ; _2 _1  _)&setdiag) ai54
(_2 _1      &(((,.~ 2&+)      1 2    )}) -: (_1 _2       ; _2 _1 _2)&setdiag) ai54
(_3 _2 _1   &(((,.~ 2&+) i.       3  )}) -: (_1 _2 _3    ; _2 _1 __)&setdiag) ai54
(_1 _2      &(((,.~ 2&+)      1 2    )}) -: (_1 _2       ; _2  1   )&setdiag) ai54
(                                        -: (''          ; _2  1  0)&setdiag) ai54
(_1 _2      &(((,.~ 2&+)      1 2    )}) -: (_1 _2       ; _2  1  2)&setdiag) ai54
(_1 _2      &(((,.~ 2&+)      1 2    )}) -: (_1 _2       ; _2  1  _)&setdiag) ai54
(_2 _1      &(((,.~ 2&+)      1 2    )}) -: (_1 _2       ; _2  1 _2)&setdiag) ai54
(_2 _1      &(((,.~ 2&+)      1 2    )}) -: (_1 _2       ; _2  1 __)&setdiag) ai54
(_1 _2      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2 _2   )&setdiag) ai54
(                                        -: (''          ; _2 _2  0)&setdiag) ai54
(_1 _2      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2 _2  2)&setdiag) ai54
(_1 _2      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2 _2  _)&setdiag) ai54
(_2 _1      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2 _2 _2)&setdiag) ai54
(_2 _1      &(((,.~ 2&+)    0 1      )}) -: (_1 _2       ; _2 _2 __)&setdiag) ai54

(_1 _2 _3 _4 _5&(( ,.~      i.           5)}) -: (_1 _2 _3 _4 _5 ; ''      )&setdiag) ai55
(_1 _2 _3 _4 _5&(( ,.~      i.           5)}) -: (_1 _2 _3 _4 _5 ;  0      )&setdiag) ai55
(_1 _2 _3 _4 _5&(( ,.~      i.           5)}) -: (_1 _2 _3 _4 _5 ;  0  0   )&setdiag) ai55
(                                             -: (''             ;  0  0  0)&setdiag) ai55
(_1 _2         &(( ,.~         0 1        )}) -: (_1 _2          ;  0  0  2)&setdiag) ai55
(_1 _2 _3 _4 _5&(( ,.~      i.           5)}) -: (_1 _2 _3 _4 _5 ;  0  0  _)&setdiag) ai55
(_2 _1         &(( ,.~         0 1        )}) -: (_1 _2          ;  0  0 _2)&setdiag) ai55
(_5 _4 _3 _2 _1&(( ,.~      i.           5)}) -: (_1 _2 _3 _4 _5 ;  0  0 __)&setdiag) ai55
(_1 _2 _3 _4 _5&(( ,.~      i.           5)}) -: (_1 _2 _3 _4 _5 ;  0 _1   )&setdiag) ai55
(                                             -: (''             ;  0 _1  0)&setdiag) ai55
(_1 _2         &(( ,.~               3 4  )}) -: (_1 _2          ;  0 _1  2)&setdiag) ai55
(_1 _2 _3 _4 _5&(( ,.~      i.           5)}) -: (_1 _2 _3 _4 _5 ;  0 _1  _)&setdiag) ai55
(_2 _1         &(( ,.~               3 4  )}) -: (_1 _2          ;  0 _1 _2)&setdiag) ai55
(_5 _4 _3 _2 _1&(( ,.~      i.           5)}) -: (_1 _2 _3 _4 _5 ;  0 _1 __)&setdiag) ai55
(_1 _2 _3 _4   &(( ,.~           1 2 3 4  )}) -: (_1 _2 _3 _4    ;  0  1   )&setdiag) ai55
(                                             -: (''             ;  0  1  0)&setdiag) ai55
(_1 _2         &(( ,.~           1 2      )}) -: (_1 _2          ;  0  1  2)&setdiag) ai55
(_1 _2 _3 _4   &(( ,.~           1 2 3 4  )}) -: (_1 _2 _3 _4    ;  0  1  _)&setdiag) ai55
(_2 _1         &(( ,.~           1 2      )}) -: (_1 _2          ;  0  1 _2)&setdiag) ai55
(_4 _3 _2 _1   &(( ,.~           1 2 3 4  )}) -: (_1 _2 _3 _4    ;  0  1 __)&setdiag) ai55
(_1 _2 _3 _4   &(( ,.~      i.         4  )}) -: (_1 _2 _3 _4    ;  0 _2   )&setdiag) ai55
(                                             -: (''             ;  0 _2  0)&setdiag) ai55
(_1 _2         &(( ,.~             2 3    )}) -: (_1 _2          ;  0 _2  2)&setdiag) ai55
(_1 _2 _3 _4   &(( ,.~      i.         4  )}) -: (_1 _2 _3 _4    ;  0 _2  _)&setdiag) ai55
(_2 _1         &(( ,.~             2 3    )}) -: (_1 _2          ;  0 _2 _2)&setdiag) ai55
(_4 _3 _2 _1   &(( ,.~      i.         4  )}) -: (_1 _2 _3 _4    ;  0 _2 __)&setdiag) ai55
(_1 _2 _3 _4   &(((,.~ >: ) i.         4  )}) -: (_1 _2 _3 _4    ; _1      )&setdiag) ai55
(_1 _2 _3 _4   &(((,.~ >: ) i.         4  )}) -: (_1 _2 _3 _4    ; _1  0   )&setdiag) ai55
(                                             -: (''             ; _1  0  0)&setdiag) ai55
(_1 _2         &(((,.~ >: )    0 1        )}) -: (_1 _2          ; _1  0  2)&setdiag) ai55
(_1 _2 _3 _4   &(((,.~ >: ) i.         4  )}) -: (_1 _2 _3 _4    ; _1  0  _)&setdiag) ai55
(_2 _1         &(((,.~ >: )    0 1        )}) -: (_1 _2          ; _1  0 _2)&setdiag) ai55
(_4 _3 _2 _1   &(((,.~ >: ) i.         4  )}) -: (_1 _2 _3 _4    ; _1  0 __)&setdiag) ai55
(_1 _2 _3 _4   &(((,.~ >: ) i.         4  )}) -: (_1 _2 _3 _4    ; _1 _1   )&setdiag) ai55
(                                             -: (''             ; _1 _1  0)&setdiag) ai55
(_1 _2         &(((,.~ >: )        2 3    )}) -: (_1 _2          ; _1 _1  2)&setdiag) ai55
(_1 _2 _3 _4   &(((,.~ >: ) i.         4  )}) -: (_1 _2 _3 _4    ; _1 _1  _)&setdiag) ai55
(_2 _1         &(((,.~ >: )        2 3    )}) -: (_1 _2          ; _1 _1 _2)&setdiag) ai55
(_4 _3 _2 _1   &(((,.~ >: ) i.         4  )}) -: (_1 _2 _3 _4    ; _1 _1 __)&setdiag) ai55
(_1 _2 _3      &(((,.~ >: )      1 2 3    )}) -: (_1 _2 _3       ; _1  1   )&setdiag) ai55
(                                             -: (''             ; _1  1  0)&setdiag) ai55
(_1 _2         &(((,.~ >: )      1 2      )}) -: (_1 _2          ; _1  1  2)&setdiag) ai55
(_1 _2 _3      &(((,.~ >: )      1 2 3    )}) -: (_1 _2 _3       ; _1  1  _)&setdiag) ai55
(_2 _1         &(((,.~ >: )      1 2      )}) -: (_1 _2          ; _1  1 _2)&setdiag) ai55
(_3 _2 _1      &(((,.~ >: )      1 2 3    )}) -: (_1 _2 _3       ; _1  1 __)&setdiag) ai55
(_1 _2 _3      &(((,.~ >: ) i.       3    )}) -: (_1 _2 _3       ; _1 _2   )&setdiag) ai55
(                                             -: (''             ; _1 _2  0)&setdiag) ai55
(_1 _2         &(((,.~ >: )      1 2      )}) -: (_1 _2          ; _1 _2  2)&setdiag) ai55
(_1 _2 _3      &(((,.~ >: ) i.       3    )}) -: (_1 _2 _3       ; _1 _2  _)&setdiag) ai55
(_2 _1         &(((,.~ >: )      1 2      )}) -: (_1 _2          ; _1 _2 _2)&setdiag) ai55
(_3 _2 _1      &(((,.~ >: ) i.       3    )}) -: (_1 _2 _3       ; _1 _2 __)&setdiag) ai55
(_1 _2 _3 _4   &(((,.  >: ) i.         4  )}) -: (_1 _2 _3 _4    ;  1      )&setdiag) ai55
(_1 _2 _3 _4   &(((,.  >: ) i.         4  )}) -: (_1 _2 _3 _4    ;  1  0   )&setdiag) ai55
(                                             -: (''             ;  1  0  0)&setdiag) ai55
(_1 _2         &(((,.  >: )    0 1        )}) -: (_1 _2          ;  1  0  2)&setdiag) ai55
(_1 _2 _3 _4   &(((,.  >: ) i.         4  )}) -: (_1 _2 _3 _4    ;  1  0  _)&setdiag) ai55
(_2 _1         &(((,.  >: )    0 1        )}) -: (_1 _2          ;  1  0 _2)&setdiag) ai55
(_4 _3 _2 _1   &(((,.  >: ) i.         4  )}) -: (_1 _2 _3 _4    ;  1  0 __)&setdiag) ai55
(_1 _2 _3 _4   &(((,.  >: ) i.         4  )}) -: (_1 _2 _3 _4    ;  1 _1   )&setdiag) ai55
(                                             -: (''             ;  1 _1  0)&setdiag) ai55
(_1 _2         &(((,.  >: )        2 3    )}) -: (_1 _2          ;  1 _1  2)&setdiag) ai55
(_1 _2 _3 _4   &(((,.  >: ) i.         4  )}) -: (_1 _2 _3 _4    ;  1 _1  _)&setdiag) ai55
(_2 _1         &(((,.  >: )        2 3    )}) -: (_1 _2          ;  1 _1 _2)&setdiag) ai55
(_4 _3 _2 _1   &(((,.  >: ) i.         4  )}) -: (_1 _2 _3 _4    ;  1 _1 __)&setdiag) ai55
(_1 _2 _3      &(((,.  >: )      1 2 3    )}) -: (_1 _2 _3       ;  1  1   )&setdiag) ai55
(                                             -: (''             ;  1  1  0)&setdiag) ai55
(_1 _2         &(((,.  >: )      1 2      )}) -: (_1 _2          ;  1  1  2)&setdiag) ai55
(_1 _2 _3      &(((,.  >: )      1 2 3    )}) -: (_1 _2 _3       ;  1  1  _)&setdiag) ai55
(_2 _1         &(((,.  >: )      1 2      )}) -: (_1 _2          ;  1  1 _2)&setdiag) ai55
(_3 _2 _1      &(((,.  >: )      1 2 3    )}) -: (_1 _2 _3       ;  1  1 __)&setdiag) ai55
(_1 _2 _3      &(((,.  >: ) i.       3    )}) -: (_1 _2 _3       ;  1 _2   )&setdiag) ai55
(                                             -: (''             ;  1 _2  0)&setdiag) ai55
(_1 _2         &(((,.  >: )      1 2      )}) -: (_1 _2          ;  1 _2  2)&setdiag) ai55
(_1 _2 _3      &(((,.  >: ) i.       3    )}) -: (_1 _2 _3       ;  1 _2  _)&setdiag) ai55
(_2 _1         &(((,.  >: )      1 2      )}) -: (_1 _2          ;  1 _2 _2)&setdiag) ai55
(_3 _2 _1      &(((,.  >: ) i.       3    )}) -: (_1 _2 _3       ;  1 _2 __)&setdiag) ai55
(_1 _2 _3      &(((,.  2&+) i.       3    )}) -: (_1 _2 _3       ;  2      )&setdiag) ai55
(_1 _2 _3      &(((,.  2&+) i.       3    )}) -: (_1 _2 _3       ;  2  0   )&setdiag) ai55
(                                             -: (''             ;  2  0  0)&setdiag) ai55
(_1 _2         &(((,.  2&+)    0 1        )}) -: (_1 _2          ;  2  0  2)&setdiag) ai55
(_1 _2 _3      &(((,.  2&+) i.       3    )}) -: (_1 _2 _3       ;  2  0  _)&setdiag) ai55
(_2 _1         &(((,.  2&+)    0 1        )}) -: (_1 _2          ;  2  0 _2)&setdiag) ai55
(_3 _2 _1      &(((,.  2&+) i.       3    )}) -: (_1 _2 _3       ;  2  0 __)&setdiag) ai55
(_1 _2 _3      &(((,.  2&+) i.       3    )}) -: (_1 _2 _3       ;  2 _1   )&setdiag) ai55
(                                             -: (''             ;  2 _1  0)&setdiag) ai55
(_1 _2         &(((,.  2&+)      1 2      )}) -: (_1 _2          ;  2 _1  2)&setdiag) ai55
(_1 _2 _3      &(((,.  2&+) i.       3    )}) -: (_1 _2 _3       ;  2 _1  _)&setdiag) ai55
(_2 _1         &(((,.  2&+)      1 2      )}) -: (_1 _2          ;  2 _1 _2)&setdiag) ai55
(_3 _2 _1      &(((,.  2&+) i.       3    )}) -: (_1 _2 _3       ;  2 _1 __)&setdiag) ai55
(_1 _2         &(((,.  2&+)      1 2      )}) -: (_1 _2          ;  2  1   )&setdiag) ai55
(                                             -: (''             ;  2  1  0)&setdiag) ai55
(_1 _2         &(((,.  2&+)      1 2      )}) -: (_1 _2          ;  2  1  2)&setdiag) ai55
(_1 _2         &(((,.  2&+)      1 2      )}) -: (_1 _2          ;  2  1  _)&setdiag) ai55
(_2 _1         &(((,.  2&+)      1 2      )}) -: (_1 _2          ;  2  1 _2)&setdiag) ai55
(_2 _1         &(((,.  2&+)      1 2      )}) -: (_1 _2          ;  2  1 __)&setdiag) ai55
(_1 _2         &(((,.  2&+)    0 1        )}) -: (_1 _2          ;  2 _2   )&setdiag) ai55
(                                             -: (''             ;  2 _2  0)&setdiag) ai55
(_1 _2         &(((,.  2&+)    0 1        )}) -: (_1 _2          ;  2 _2  2)&setdiag) ai55
(_1 _2         &(((,.  2&+)    0 1        )}) -: (_1 _2          ;  2 _2  _)&setdiag) ai55
(_2 _1         &(((,.  2&+)    0 1        )}) -: (_1 _2          ;  2 _2 _2)&setdiag) ai55
(_2 _1         &(((,.  2&+)    0 1        )}) -: (_1 _2          ;  2 _2 __)&setdiag) ai55
(_1 _2 _3      &(((,.~ 2&+) i.       3    )}) -: (_1 _2 _3       ; _2      )&setdiag) ai55
(_1 _2 _3      &(((,.~ 2&+) i.       3    )}) -: (_1 _2 _3       ; _2  0   )&setdiag) ai55
(                                             -: (''             ; _2  0  0)&setdiag) ai55
(_1 _2         &(((,.~ 2&+)    0 1        )}) -: (_1 _2          ; _2  0  2)&setdiag) ai55
(_1 _2 _3      &(((,.~ 2&+) i.       3    )}) -: (_1 _2 _3       ; _2  0  _)&setdiag) ai55
(_2 _1         &(((,.~ 2&+)    0 1        )}) -: (_1 _2          ; _2  0 _2)&setdiag) ai55
(_3 _2 _1      &(((,.~ 2&+) i.       3    )}) -: (_1 _2 _3       ; _2  0 __)&setdiag) ai55
(_1 _2 _3      &(((,.~ 2&+) i.       3    )}) -: (_1 _2 _3       ; _2 _1   )&setdiag) ai55
(                                             -: (''             ; _2 _1  0)&setdiag) ai55
(_1 _2         &(((,.~ 2&+)      1 2      )}) -: (_1 _2          ; _2 _1  2)&setdiag) ai55
(_1 _2 _3      &(((,.~ 2&+) i.       3    )}) -: (_1 _2 _3       ; _2 _1  _)&setdiag) ai55
(_2 _1         &(((,.~ 2&+)      1 2      )}) -: (_1 _2          ; _2 _1 _2)&setdiag) ai55
(_3 _2 _1      &(((,.~ 2&+) i.       3    )}) -: (_1 _2 _3       ; _2 _1 __)&setdiag) ai55
(_1 _2         &(((,.~ 2&+)      1 2      )}) -: (_1 _2          ; _2  1   )&setdiag) ai55
(                                             -: (''             ; _2  1  0)&setdiag) ai55
(_1 _2         &(((,.~ 2&+)      1 2      )}) -: (_1 _2          ; _2  1  2)&setdiag) ai55
(_1 _2         &(((,.~ 2&+)      1 2      )}) -: (_1 _2          ; _2  1  _)&setdiag) ai55
(_2 _1         &(((,.~ 2&+)      1 2      )}) -: (_1 _2          ; _2  1 _2)&setdiag) ai55
(_2 _1         &(((,.~ 2&+)      1 2      )}) -: (_1 _2          ; _2  1 __)&setdiag) ai55
(_1 _2         &(((,.~ 2&+)    0 1        )}) -: (_1 _2          ; _2 _2   )&setdiag) ai55
(                                             -: (''             ; _2 _2  0)&setdiag) ai55
(_1 _2         &(((,.~ 2&+)    0 1        )}) -: (_1 _2          ; _2 _2  2)&setdiag) ai55
(_1 _2         &(((,.~ 2&+)    0 1        )}) -: (_1 _2          ; _2 _2  _)&setdiag) ai55
(_2 _1         &(((,.~ 2&+)    0 1        )}) -: (_1 _2          ; _2 _2 _2)&setdiag) ai55
(_2 _1         &(((,.~ 2&+)    0 1        )}) -: (_1 _2          ; _2 _2 __)&setdiag) ai55

NB. upddiag
(-:     - upddiag ) EMPTY
(-:  1&(- upddiag)) EMPTY
(-: _1&(- upddiag)) EMPTY
(3 4 $ 0  1 2 3  4 _5  6 7 8  9 _10  11) -:        - upddiag ai34
(3 4 $ 0  1 2 3  4 _5  6 7 8  9 _10  11) -:  0     - upddiag ai34
(3 4 $ 0 _1 2 3  4  5 _6 7 8  9  10 _11) -:  1     - upddiag ai34
(3 4 $ 0  1 2 3  4  5 _6 7 8  9  10 _11) -:  1 1   - upddiag ai34
(3 4 $ 0  1 2 3  4  5 _6 7 8  9  10  11) -:  1 1 1 - upddiag ai34
(3 4 $ 0  1 2 3  4  5 _6 7 8  9  10 _11) -:  1 1 2 - upddiag ai34
(3 4 $ 0  1 2 3 _4  5  6 7 8 _9  10  11) -: _1     - upddiag ai34
(3 4 $ 0  1 2 3  4  5  6 7 8 _9  10  11) -: _1 1   - upddiag ai34

NB. bdlpick
(-: bdlpick) EMPTY
(4 5 $ 1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 1 1 0) -: bdlpick a145
(5 4 $ 1 0 0 0 1 1 0 0 0 1 1 0 0 0 1 1 0 0 0 1) -: bdlpick a154
(5 5 $ 1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 1 1) -: bdlpick a155

NB. bdupick
(-: bdupick) EMPTY
(4 5 $ 1 1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 1 1) -: bdupick a145
(5 4 $ 1 1 0 0 0 1 1 0 0 0 1 1 0 0 0 1 0 0 0 0) -: bdupick a154
(5 5 $ 1 1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 1) -: bdupick a155

NB. hslpick
(-: hslpick) EMPTY
(4 5 $ 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1) -: hslpick a145
(5 4 $ 1 1 0 0 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1) -: hslpick a154
(5 5 $ 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1) -: hslpick a155

NB. hsupick
(-: hsupick) EMPTY
(4 5 $ 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1) -: hsupick a145
(5 4 $ 1 1 1 1 1 1 1 1 0 1 1 1 0 0 1 1 0 0 0 1) -: hsupick a154
(5 5 $ 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1) -: hsupick a155

NB. gtpick
(-: gtpick) EMPTY
(4 5 $ 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1) -: gtpick a145
(5 4 $ 1 1 0 0 1 1 1 0 0 1 1 1 0 0 1 1 0 0 0 1) -: gtpick a154
(5 5 $ 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1) -: gtpick a155

NB. trlpick
(-:    trlpick) EMPTY
(-:  0&trlpick) EMPTY
(-:  1&trlpick) EMPTY
(-: _1&trlpick) EMPTY
(4 5 $ 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0) -:    trlpick a145
(4 5 $ 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0) -:  0 trlpick a145
(4 5 $ 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1) -:  1 trlpick a145
(4 5 $ 0 0 0 0 0 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0) -: _1 trlpick a145
(5 4 $ 1 0 0 0 1 1 0 0 1 1 1 0 1 1 1 1 1 1 1 1) -:    trlpick a154
(5 4 $ 1 0 0 0 1 1 0 0 1 1 1 0 1 1 1 1 1 1 1 1) -:  0 trlpick a154
(5 4 $ 1 1 0 0 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1) -:  1 trlpick a154
(5 4 $ 0 0 0 0 1 0 0 0 1 1 0 0 1 1 1 0 1 1 1 1) -: _1 trlpick a154
(5 5 $ 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1) -:    trlpick a155
(5 5 $ 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1) -:  0 trlpick a155
(5 5 $ 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1) -:  1 trlpick a155
(5 5 $ 0 0 0 0 0 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0) -: _1 trlpick a155

NB. trupick
(-:    trupick) EMPTY
(-:  0&trupick) EMPTY
(-:  1&trupick) EMPTY
(-: _1&trupick) EMPTY
(4 5 $ 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1) -:    trupick a145
(4 5 $ 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1) -:  0 trupick a145
(4 5 $ 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1 0 0 0 0 1) -:  1 trupick a145
(4 5 $ 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1) -: _1 trupick a145
(5 4 $ 1 1 1 1 0 1 1 1 0 0 1 1 0 0 0 1 0 0 0 0) -:    trupick a154
(5 4 $ 1 1 1 1 0 1 1 1 0 0 1 1 0 0 0 1 0 0 0 0) -:  0 trupick a154
(5 4 $ 0 1 1 1 0 0 1 1 0 0 0 1 0 0 0 0 0 0 0 0) -:  1 trupick a154
(5 4 $ 1 1 1 1 1 1 1 1 0 1 1 1 0 0 1 1 0 0 0 1) -: _1 trupick a154
(5 5 $ 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1 0 0 0 0 1) -:    trupick a155
(5 5 $ 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1 0 0 0 0 1) -:  0 trupick a155
(5 5 $ 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0) -:  1 trupick a155
(5 5 $ 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1) -: _1 trupick a155

NB. trl1pick
(-:    trl1pick) EMPTY
(-:  0&trl1pick) EMPTY
(-:  1&trl1pick) EMPTY
(-: _1&trl1pick) EMPTY
(4 5 $ 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0) -:    trl1pick a145
(4 5 $ 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0) -:  0 trl1pick a145
(4 5 $ 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1) -:  1 trl1pick a145
(4 5 $ 0 0 0 0 0 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0) -: _1 trl1pick a145
(5 4 $ 1 0 0 0 1 1 0 0 1 1 1 0 1 1 1 1 1 1 1 1) -:    trl1pick a154
(5 4 $ 1 0 0 0 1 1 0 0 1 1 1 0 1 1 1 1 1 1 1 1) -:  0 trl1pick a154
(5 4 $ 1 1 0 0 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1) -:  1 trl1pick a154
(5 4 $ 0 0 0 0 1 0 0 0 1 1 0 0 1 1 1 0 1 1 1 1) -: _1 trl1pick a154
(5 5 $ 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1) -:    trl1pick a155
(5 5 $ 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1) -:  0 trl1pick a155
(5 5 $ 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1) -:  1 trl1pick a155
(5 5 $ 0 0 0 0 0 1 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0) -: _1 trl1pick a155

NB. tru1pick
(-:    tru1pick) EMPTY
(-:  0&tru1pick) EMPTY
(-:  1&tru1pick) EMPTY
(-: _1&tru1pick) EMPTY
(4 5 $ 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1) -:    tru1pick a145
(4 5 $ 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1) -:  0 tru1pick a145
(4 5 $ 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1 0 0 0 0 1) -:  1 tru1pick a145
(4 5 $ 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1) -: _1 tru1pick a145
(5 4 $ 1 1 1 1 0 1 1 1 0 0 1 1 0 0 0 1 0 0 0 0) -:    tru1pick a154
(5 4 $ 1 1 1 1 0 1 1 1 0 0 1 1 0 0 0 1 0 0 0 0) -:  0 tru1pick a154
(5 4 $ 0 1 1 1 0 0 1 1 0 0 0 1 0 0 0 0 0 0 0 0) -:  1 tru1pick a154
(5 4 $ 1 1 1 1 1 1 1 1 0 1 1 1 0 0 1 1 0 0 0 1) -: _1 tru1pick a154
(5 5 $ 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1 0 0 0 0 1) -:    tru1pick a155
(5 5 $ 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1 0 0 0 0 1) -:  0 tru1pick a155
(5 5 $ 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0) -:  1 tru1pick a155
(5 5 $ 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 0 0 0 1 1) -: _1 tru1pick a155

NB. idmat

EMPTY -:         idmat 0
EMPTY -: 0       idmat 0
EMPTY -: 0  0    idmat 0
EMPTY -: 0  0  0 idmat 0
EMPTY -: 0  0  _ idmat 0
EMPTY -: 0  0 __ idmat 0
EMPTY -: 0 _1    idmat 0
EMPTY -: 0 _1  0 idmat 0
EMPTY -: 0 _1  _ idmat 0
EMPTY -: 0 _1 __ idmat 0

EMPTY -:         idmat 0 0
EMPTY -: 0       idmat 0 0
EMPTY -: 0  0    idmat 0 0
EMPTY -: 0  0  0 idmat 0 0
EMPTY -: 0  0  _ idmat 0 0
EMPTY -: 0  0 __ idmat 0 0
EMPTY -: 0 _1    idmat 0 0
EMPTY -: 0 _1  0 idmat 0 0
EMPTY -: 0 _1  _ idmat 0 0
EMPTY -: 0 _1 __ idmat 0 0

a005 -:         idmat 0 5
a005 -: 0       idmat 0 5
a005 -: 0  0    idmat 0 5
a005 -: 0  0  0 idmat 0 5
a005 -: 0  0  _ idmat 0 5
a005 -: 0  0 __ idmat 0 5
a005 -: 0 _1    idmat 0 5
a005 -: 0 _1  0 idmat 0 5
a005 -: 0 _1  _ idmat 0 5
a005 -: 0 _1 __ idmat 0 5

a050 -:         idmat 5 0
a050 -: 0       idmat 5 0
a050 -: 0  0    idmat 5 0
a050 -: 0  0  0 idmat 5 0
a050 -: 0  0  _ idmat 5 0
a050 -: 0  0 __ idmat 5 0
a050 -: 0 _1    idmat 5 0
a050 -: 0 _1  0 idmat 5 0
a050 -: 0 _1  _ idmat 5 0
a050 -: 0 _1 __ idmat 5 0

(1 ( ,.~      i.         4)} a045) -:          idmat 4 5
(1 ( ,.~      i.         4)} a045) -:  0       idmat 4 5
(1 ( ,.~      i.         4)} a045) -:  0  0    idmat 4 5
                             a045  -:  0  0  0 idmat 4 5
(1 ( ,.~      i.     2    )} a045) -:  0  0  2 idmat 4 5
(1 ( ,.~      i.         4)} a045) -:  0  0  _ idmat 4 5
(1 ( ,.~      i.     2    )} a045) -:  0  0 _2 idmat 4 5
(1 ( ,.~      i.         4)} a045) -:  0  0 __ idmat 4 5
(1 ( ,.~      i.         4)} a045) -:  0 _1    idmat 4 5
                             a045  -:  0 _1  0 idmat 4 5
(1 ( ,.~             2 3  )} a045) -:  0 _1  2 idmat 4 5
(1 ( ,.~      i.         4)} a045) -:  0 _1  _ idmat 4 5
(1 ( ,.~             2 3  )} a045) -:  0 _1 _2 idmat 4 5
(1 ( ,.~      i.         4)} a045) -:  0 _1 __ idmat 4 5
(1 ( ,.~           1 2 3  )} a045) -:  0  1    idmat 4 5
                             a045  -:  0  1  0 idmat 4 5
(1 ( ,.~           1 2    )} a045) -:  0  1  2 idmat 4 5
(1 ( ,.~           1 2 3  )} a045) -:  0  1  _ idmat 4 5
(1 ( ,.~           1 2    )} a045) -:  0  1 _2 idmat 4 5
(1 ( ,.~           1 2 3  )} a045) -:  0  1 __ idmat 4 5
(1 ( ,.~      i.       3  )} a045) -:  0 _2    idmat 4 5
                             a045  -:  0 _2  0 idmat 4 5
(1 ( ,.~           1 2    )} a045) -:  0 _2  2 idmat 4 5
(1 ( ,.~      i.       3  )} a045) -:  0 _2  _ idmat 4 5
(1 ( ,.~           1 2    )} a045) -:  0 _2 _2 idmat 4 5
(1 ( ,.~      i.       3  )} a045) -:  0 _2 __ idmat 4 5
(1 ((,.~ >: ) i.       3  )} a045) -: _1       idmat 4 5
(1 ((,.~ >: ) i.       3  )} a045) -: _1  0    idmat 4 5
                             a045  -: _1  0  0 idmat 4 5
(1 ((,.~ >: )    0 1      )} a045) -: _1  0  2 idmat 4 5
(1 ((,.~ >: ) i.       3  )} a045) -: _1  0  _ idmat 4 5
(1 ((,.~ >: )    0 1      )} a045) -: _1  0 _2 idmat 4 5
(1 ((,.~ >: ) i.       3  )} a045) -: _1  0 __ idmat 4 5
(1 ((,.~ >: ) i.       3  )} a045) -: _1 _1    idmat 4 5
                             a045  -: _1 _1  0 idmat 4 5
(1 ((,.~ >: )      1 2    )} a045) -: _1 _1  2 idmat 4 5
(1 ((,.~ >: ) i.       3  )} a045) -: _1 _1  _ idmat 4 5
(1 ((,.~ >: )      1 2    )} a045) -: _1 _1 _2 idmat 4 5
(1 ((,.~ >: ) i.       3  )} a045) -: _1 _1 __ idmat 4 5
(1 ((,.~ >: )      1 2    )} a045) -: _1  1    idmat 4 5
                             a045  -: _1  1  0 idmat 4 5
(1 ((,.~ >: )      1 2    )} a045) -: _1  1  2 idmat 4 5
(1 ((,.~ >: )      1 2    )} a045) -: _1  1  _ idmat 4 5
(1 ((,.~ >: )      1 2    )} a045) -: _1  1 _2 idmat 4 5
(1 ((,.~ >: )      1 2    )} a045) -: _1  1 __ idmat 4 5
(1 ((,.~ >: )    0 1      )} a045) -: _1 _2    idmat 4 5
                             a045  -: _1 _2  0 idmat 4 5
(1 ((,.~ >: )    0 1      )} a045) -: _1 _2  2 idmat 4 5
(1 ((,.~ >: )    0 1      )} a045) -: _1 _2  _ idmat 4 5
(1 ((,.~ >: )    0 1      )} a045) -: _1 _2 _2 idmat 4 5
(1 ((,.~ >: )    0 1      )} a045) -: _1 _2 __ idmat 4 5
(1 ((,.  >: ) i.         4)} a045) -:  1       idmat 4 5
(1 ((,.  >: ) i.         4)} a045) -:  1  0    idmat 4 5
                             a045  -:  1  0  0 idmat 4 5
(1 ((,.  >: )    0 1      )} a045) -:  1  0  2 idmat 4 5
(1 ((,.  >: ) i.         4)} a045) -:  1  0  _ idmat 4 5
(1 ((,.  >: )    0 1      )} a045) -:  1  0 _2 idmat 4 5
(1 ((,.  >: ) i.         4)} a045) -:  1  0 __ idmat 4 5
(1 ((,.  >: ) i.         4)} a045) -:  1 _1    idmat 4 5
                             a045  -:  1 _1  0 idmat 4 5
(1 ((,.  >: )        2 3  )} a045) -:  1 _1  2 idmat 4 5
(1 ((,.  >: ) i.         4)} a045) -:  1 _1  _ idmat 4 5
(1 ((,.  >: )        2 3  )} a045) -:  1 _1 _2 idmat 4 5
(1 ((,.  >: ) i.         4)} a045) -:  1 _1 __ idmat 4 5
(1 ((,.  >: )      1 2 3  )} a045) -:  1  1    idmat 4 5
                             a045  -:  1  1  0 idmat 4 5
(1 ((,.  >: )      1 2    )} a045) -:  1  1  2 idmat 4 5
(1 ((,.  >: )      1 2 3  )} a045) -:  1  1  _ idmat 4 5
(1 ((,.  >: )      1 2    )} a045) -:  1  1 _2 idmat 4 5
(1 ((,.  >: )      1 2 3  )} a045) -:  1  1 __ idmat 4 5
(1 ((,.  >: ) i.       3  )} a045) -:  1 _2    idmat 4 5
                             a045  -:  1 _2  0 idmat 4 5
(1 ((,.  >: )      1 2    )} a045) -:  1 _2  2 idmat 4 5
(1 ((,.  >: ) i.       3  )} a045) -:  1 _2  _ idmat 4 5
(1 ((,.  >: )      1 2    )} a045) -:  1 _2 _2 idmat 4 5
(1 ((,.  >: ) i.       3  )} a045) -:  1 _2 __ idmat 4 5
(1 ((,.  2&+) i.       3  )} a045) -:  2       idmat 4 5
(1 ((,.  2&+) i.       3  )} a045) -:  2  0    idmat 4 5
                             a045  -:  2  0  0 idmat 4 5
(1 ((,.  2&+)    0 1      )} a045) -:  2  0  2 idmat 4 5
(1 ((,.  2&+) i.       3  )} a045) -:  2  0  _ idmat 4 5
(1 ((,.  2&+)    0 1      )} a045) -:  2  0 _2 idmat 4 5
(1 ((,.  2&+) i.       3  )} a045) -:  2  0 __ idmat 4 5
(1 ((,.  2&+) i.       3  )} a045) -:  2 _1    idmat 4 5
                             a045  -:  2 _1  0 idmat 4 5
(1 ((,.  2&+)      1 2    )} a045) -:  2 _1  2 idmat 4 5
(1 ((,.  2&+) i.       3  )} a045) -:  2 _1  _ idmat 4 5
(1 ((,.  2&+)      1 2    )} a045) -:  2 _1 _2 idmat 4 5
(1 ((,.  2&+) i.       3  )} a045) -:  2 _1 __ idmat 4 5
(1 ((,.  2&+)      1 2    )} a045) -:  2  1    idmat 4 5
                             a045  -:  2  1  0 idmat 4 5
(1 ((,.  2&+)      1 2    )} a045) -:  2  1  2 idmat 4 5
(1 ((,.  2&+)      1 2    )} a045) -:  2  1  _ idmat 4 5
(1 ((,.  2&+)      1 2    )} a045) -:  2  1 _2 idmat 4 5
(1 ((,.  2&+)      1 2    )} a045) -:  2  1 __ idmat 4 5
(1 ((,.  2&+)    0 1      )} a045) -:  2 _2    idmat 4 5
                             a045  -:  2 _2  0 idmat 4 5
(1 ((,.  2&+)    0 1      )} a045) -:  2 _2  2 idmat 4 5
(1 ((,.  2&+)    0 1      )} a045) -:  2 _2  _ idmat 4 5
(1 ((,.  2&+)    0 1      )} a045) -:  2 _2 _2 idmat 4 5
(1 ((,.  2&+)    0 1      )} a045) -:  2 _2 __ idmat 4 5
(1 ((,.~ 2&+)    0 1      )} a045) -: _2       idmat 4 5
(1 ((,.~ 2&+)    0 1      )} a045) -: _2  0    idmat 4 5
                             a045  -: _2  0  0 idmat 4 5
(1 ((,.~ 2&+)    0 1      )} a045) -: _2  0  2 idmat 4 5
(1 ((,.~ 2&+)    0 1      )} a045) -: _2  0  _ idmat 4 5
(1 ((,.~ 2&+)    0 1      )} a045) -: _2  0 _2 idmat 4 5
(1 ((,.~ 2&+)    0 1      )} a045) -: _2  0 __ idmat 4 5
(1 ((,.~ 2&+)    0 1      )} a045) -: _2 _1    idmat 4 5
                             a045  -: _2 _1  0 idmat 4 5
(1 ((,.~ 2&+)    0 1      )} a045) -: _2 _1  2 idmat 4 5
(1 ((,.~ 2&+)    0 1      )} a045) -: _2 _1  _ idmat 4 5
(1 ((,.~ 2&+)    0 1      )} a045) -: _2 _1 _2 idmat 4 5
(1 ((,.~ 2&+)    0 1      )} a045) -: _2 _1 __ idmat 4 5
(1 (< 3 1                 )} a045) -: _2  1    idmat 4 5
                             a045  -: _2  1  0 idmat 4 5
(1 (< 3 1                 )} a045) -: _2  1  2 idmat 4 5
(1 (< 3 1                 )} a045) -: _2  1  _ idmat 4 5
(1 (< 3 1                 )} a045) -: _2  1 _2 idmat 4 5
(1 (< 3 1                 )} a045) -: _2  1 __ idmat 4 5
(1 (< 2 0                 )} a045) -: _2 _2    idmat 4 5
                             a045  -: _2 _2  0 idmat 4 5
(1 (< 2 0                 )} a045) -: _2 _2  2 idmat 4 5
(1 (< 2 0                 )} a045) -: _2 _2  _ idmat 4 5
(1 (< 2 0                 )} a045) -: _2 _2 _2 idmat 4 5
(1 (< 2 0                 )} a045) -: _2 _2 __ idmat 4 5

(1 ( ,.~      i.         4)} a054) -:          idmat 5 4
(1 ( ,.~      i.         4)} a054) -:  0       idmat 5 4
(1 ( ,.~      i.         4)} a054) -:  0  0    idmat 5 4
                             a054  -:  0  0  0 idmat 5 4
(1 ( ,.~         0 1      )} a054) -:  0  0  2 idmat 5 4
(1 ( ,.~      i.         4)} a054) -:  0  0  _ idmat 5 4
(1 ( ,.~         0 1      )} a054) -:  0  0 _2 idmat 5 4
(1 ( ,.~      i.         4)} a054) -:  0  0 __ idmat 5 4
(1 ( ,.~      i.         4)} a054) -:  0 _1    idmat 5 4
                             a054  -:  0 _1  0 idmat 5 4
(1 ( ,.~             2 3  )} a054) -:  0 _1  2 idmat 5 4
(1 ( ,.~      i.         4)} a054) -:  0 _1  _ idmat 5 4
(1 ( ,.~             2 3  )} a054) -:  0 _1 _2 idmat 5 4
(1 ( ,.~      i.         4)} a054) -:  0 _1 __ idmat 5 4
(1 ( ,.~           1 2 3  )} a054) -:  0  1    idmat 5 4
                             a054  -:  0  1  0 idmat 5 4
(1 ( ,.~           1 2    )} a054) -:  0  1  2 idmat 5 4
(1 ( ,.~           1 2 3  )} a054) -:  0  1  _ idmat 5 4
(1 ( ,.~           1 2    )} a054) -:  0  1 _2 idmat 5 4
(1 ( ,.~           1 2 3  )} a054) -:  0  1 __ idmat 5 4
(1 ( ,.~     i.        3  )} a054) -:  0 _2    idmat 5 4
                             a054  -:  0 _2  0 idmat 5 4
(1 ( ,.~           1 2    )} a054) -:  0 _2  2 idmat 5 4
(1 ( ,.~     i.        3  )} a054) -:  0 _2  _ idmat 5 4
(1 ( ,.~           1 2    )} a054) -:  0 _2 _2 idmat 5 4
(1 ( ,.~     i.        3  )} a054) -:  0 _2 __ idmat 5 4
(1 ((,.~ >: ) i.         4)} a054) -: _1       idmat 5 4
(1 ((,.~ >: ) i.         4)} a054) -: _1  0    idmat 5 4
                             a054  -: _1  0  0 idmat 5 4
(1 ((,.~ >: )    0 1      )} a054) -: _1  0  2 idmat 5 4
(1 ((,.~ >: ) i.         4)} a054) -: _1  0  _ idmat 5 4
(1 ((,.~ >: )    0 1      )} a054) -: _1  0 _2 idmat 5 4
(1 ((,.~ >: ) i.         4)} a054) -: _1  0 __ idmat 5 4
(1 ((,.~ >: ) i.         4)} a054) -: _1 _1    idmat 5 4
                             a054  -: _1 _1  0 idmat 5 4
(1 ((,.~ >: )        2 3  )} a054) -: _1 _1  2 idmat 5 4
(1 ((,.~ >: ) i.         4)} a054) -: _1 _1  _ idmat 5 4
(1 ((,.~ >: )        2 3  )} a054) -: _1 _1 _2 idmat 5 4
(1 ((,.~ >: ) i.         4)} a054) -: _1 _1 __ idmat 5 4
(1 ((,.~ >: )      1 2 3  )} a054) -: _1  1    idmat 5 4
                             a054  -: _1  1  0 idmat 5 4
(1 ((,.~ >: )      1 2    )} a054) -: _1  1  2 idmat 5 4
(1 ((,.~ >: )      1 2 3  )} a054) -: _1  1  _ idmat 5 4
(1 ((,.~ >: )      1 2    )} a054) -: _1  1 _2 idmat 5 4
(1 ((,.~ >: )      1 2 3  )} a054) -: _1  1 __ idmat 5 4
(1 ((,.~ >: ) i.       3  )} a054) -: _1 _2    idmat 5 4
                             a054  -: _1 _2  0 idmat 5 4
(1 ((,.~ >: )      1 2    )} a054) -: _1 _2  2 idmat 5 4
(1 ((,.~ >: ) i.       3  )} a054) -: _1 _2  _ idmat 5 4
(1 ((,.~ >: )      1 2    )} a054) -: _1 _2 _2 idmat 5 4
(1 ((,.~ >: ) i.       3  )} a054) -: _1 _2 __ idmat 5 4
(1 ((,.  >: ) i.       3  )} a054) -:  1       idmat 5 4
(1 ((,.  >: ) i.       3  )} a054) -:  1  0    idmat 5 4
                             a054  -:  1  0  0 idmat 5 4
(1 ((,.  >: )    0 1      )} a054) -:  1  0  2 idmat 5 4
(1 ((,.  >: ) i.       3  )} a054) -:  1  0  _ idmat 5 4
(1 ((,.  >: )    0 1      )} a054) -:  1  0 _2 idmat 5 4
(1 ((,.  >: ) i.       3  )} a054) -:  1  0 __ idmat 5 4
(1 ((,.  >: ) i.       3  )} a054) -:  1 _1    idmat 5 4
                             a054  -:  1 _1  0 idmat 5 4
(1 ((,.  >: )      1 2    )} a054) -:  1 _1  2 idmat 5 4
(1 ((,.  >: ) i.       3  )} a054) -:  1 _1  _ idmat 5 4
(1 ((,.  >: )      1 2    )} a054) -:  1 _1 _2 idmat 5 4
(1 ((,.  >: ) i.       3  )} a054) -:  1 _1 __ idmat 5 4
(1 ((,.  >: )      1 2    )} a054) -:  1  1    idmat 5 4
                             a054  -:  1  1  0 idmat 5 4
(1 ((,.  >: )      1 2    )} a054) -:  1  1  2 idmat 5 4
(1 ((,.  >: )      1 2    )} a054) -:  1  1  _ idmat 5 4
(1 ((,.  >: )      1 2    )} a054) -:  1  1 _2 idmat 5 4
(1 ((,.  >: )      1 2    )} a054) -:  1  1 __ idmat 5 4
(1 ((,.  >: )    0 1      )} a054) -:  1 _2    idmat 5 4
                             a054  -:  1 _2  0 idmat 5 4
(1 ((,.  >: )    0 1      )} a054) -:  1 _2  2 idmat 5 4
(1 ((,.  >: )    0 1      )} a054) -:  1 _2  _ idmat 5 4
(1 ((,.  >: )    0 1      )} a054) -:  1 _2 _2 idmat 5 4
(1 ((,.  >: )    0 1      )} a054) -:  1 _2 __ idmat 5 4
(1 ((,.  2&+)    0 1      )} a054) -:  2       idmat 5 4
(1 ((,.  2&+)    0 1      )} a054) -:  2  0    idmat 5 4
                             a054  -:  2  0  0 idmat 5 4
(1 ((,.  2&+)    0 1      )} a054) -:  2  0  2 idmat 5 4
(1 ((,.  2&+)    0 1      )} a054) -:  2  0  _ idmat 5 4
(1 ((,.  2&+)    0 1      )} a054) -:  2  0 _2 idmat 5 4
(1 ((,.  2&+)    0 1      )} a054) -:  2  0 __ idmat 5 4
(1 ((,.  2&+)    0 1      )} a054) -:  2 _1    idmat 5 4
                             a054  -:  2 _1  0 idmat 5 4
(1 ((,.  2&+)    0 1      )} a054) -:  2 _1  2 idmat 5 4
(1 ((,.  2&+)    0 1      )} a054) -:  2 _1  _ idmat 5 4
(1 ((,.  2&+)    0 1      )} a054) -:  2 _1 _2 idmat 5 4
(1 ((,.  2&+)    0 1      )} a054) -:  2 _1 __ idmat 5 4
(1 (< 1 3                 )} a054) -:  2  1    idmat 5 4
                             a054  -:  2  1  0 idmat 5 4
(1 (< 1 3                 )} a054) -:  2  1  2 idmat 5 4
(1 (< 1 3                 )} a054) -:  2  1  _ idmat 5 4
(1 (< 1 3                 )} a054) -:  2  1 _2 idmat 5 4
(1 (< 1 3                 )} a054) -:  2  1 __ idmat 5 4
(1 (< 0 2                 )} a054) -:  2 _2    idmat 5 4
                             a054  -:  2 _2  0 idmat 5 4
(1 (< 0 2                 )} a054) -:  2 _2  2 idmat 5 4
(1 (< 0 2                 )} a054) -:  2 _2  _ idmat 5 4
(1 (< 0 2                 )} a054) -:  2 _2 _2 idmat 5 4
(1 (< 0 2                 )} a054) -:  2 _2 __ idmat 5 4
(1 ((,.~ 2&+) i.       3  )} a054) -: _2       idmat 5 4
(1 ((,.~ 2&+) i.       3  )} a054) -: _2  0    idmat 5 4
                             a054  -: _2  0  0 idmat 5 4
(1 ((,.~ 2&+)    0 1      )} a054) -: _2  0  2 idmat 5 4
(1 ((,.~ 2&+) i.       3  )} a054) -: _2  0  _ idmat 5 4
(1 ((,.~ 2&+)    0 1      )} a054) -: _2  0 _2 idmat 5 4
(1 ((,.~ 2&+) i.       3  )} a054) -: _2  0 __ idmat 5 4
(1 ((,.~ 2&+) i.       3  )} a054) -: _2 _1    idmat 5 4
                             a054  -: _2 _1  0 idmat 5 4
(1 ((,.~ 2&+)      1 2    )} a054) -: _2 _1  2 idmat 5 4
(1 ((,.~ 2&+) i.       3  )} a054) -: _2 _1  _ idmat 5 4
(1 ((,.~ 2&+)      1 2    )} a054) -: _2 _1 _2 idmat 5 4
(1 ((,.~ 2&+) i.       3  )} a054) -: _2 _1 __ idmat 5 4
(1 ((,.~ 2&+)      1 2    )} a054) -: _2  1    idmat 5 4
                             a054  -: _2  1  0 idmat 5 4
(1 ((,.~ 2&+)      1 2    )} a054) -: _2  1  2 idmat 5 4
(1 ((,.~ 2&+)      1 2    )} a054) -: _2  1  _ idmat 5 4
(1 ((,.~ 2&+)      1 2    )} a054) -: _2  1 _2 idmat 5 4
(1 ((,.~ 2&+)      1 2    )} a054) -: _2  1 __ idmat 5 4
(1 ((,.~ 2&+)    0 1      )} a054) -: _2 _2    idmat 5 4
                             a054  -: _2 _2  0 idmat 5 4
(1 ((,.~ 2&+)    0 1      )} a054) -: _2 _2  2 idmat 5 4
(1 ((,.~ 2&+)    0 1      )} a054) -: _2 _2  _ idmat 5 4
(1 ((,.~ 2&+)    0 1      )} a054) -: _2 _2 _2 idmat 5 4
(1 ((,.~ 2&+)    0 1      )} a054) -: _2 _2 __ idmat 5 4

(1 ( ,.~      i.           5)} a055) -:          idmat 5 5
(1 ( ,.~      i.           5)} a055) -:  0       idmat 5 5
(1 ( ,.~      i.           5)} a055) -:  0  0    idmat 5 5
                               a055  -:  0  0  0 idmat 5 5
(1 ( ,.~         0 1        )} a055) -:  0  0  2 idmat 5 5
(1 ( ,.~      i.           5)} a055) -:  0  0  _ idmat 5 5
(1 ( ,.~         0 1        )} a055) -:  0  0 _2 idmat 5 5
(1 ( ,.~      i.           5)} a055) -:  0  0 __ idmat 5 5
(1 ( ,.~      i.           5)} a055) -:  0 _1    idmat 5 5
                               a055  -:  0 _1  0 idmat 5 5
(1 ( ,.~               3 4  )} a055) -:  0 _1  2 idmat 5 5
(1 ( ,.~      i.           5)} a055) -:  0 _1  _ idmat 5 5
(1 ( ,.~               3 4  )} a055) -:  0 _1 _2 idmat 5 5
(1 ( ,.~      i.           5)} a055) -:  0 _1 __ idmat 5 5
(1 ( ,.~           1 2 3 4  )} a055) -:  0  1    idmat 5 5
                               a055  -:  0  1  0 idmat 5 5
(1 ( ,.~           1 2      )} a055) -:  0  1  2 idmat 5 5
(1 ( ,.~           1 2 3 4  )} a055) -:  0  1  _ idmat 5 5
(1 ( ,.~           1 2      )} a055) -:  0  1 _2 idmat 5 5
(1 ( ,.~           1 2 3 4  )} a055) -:  0  1 __ idmat 5 5
(1 ( ,.~      i.         4  )} a055) -:  0 _2    idmat 5 5
                               a055  -:  0 _2  0 idmat 5 5
(1 ( ,.~             2 3    )} a055) -:  0 _2  2 idmat 5 5
(1 ( ,.~      i.         4  )} a055) -:  0 _2  _ idmat 5 5
(1 ( ,.~             2 3    )} a055) -:  0 _2 _2 idmat 5 5
(1 ( ,.~      i.         4  )} a055) -:  0 _2 __ idmat 5 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1       idmat 5 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1  0    idmat 5 5
                               a055  -: _1  0  0 idmat 5 5
(1 ((,.~ >: )    0 1        )} a055) -: _1  0  2 idmat 5 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1  0  _ idmat 5 5
(1 ((,.~ >: )    0 1        )} a055) -: _1  0 _2 idmat 5 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1  0 __ idmat 5 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1 _1    idmat 5 5
                               a055  -: _1 _1  0 idmat 5 5
(1 ((,.~ >: )        2 3    )} a055) -: _1 _1  2 idmat 5 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1 _1  _ idmat 5 5
(1 ((,.~ >: )        2 3    )} a055) -: _1 _1 _2 idmat 5 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1 _1 __ idmat 5 5
(1 ((,.~ >: )      1 2 3    )} a055) -: _1  1    idmat 5 5
                               a055  -: _1  1  0 idmat 5 5
(1 ((,.~ >: )      1 2      )} a055) -: _1  1  2 idmat 5 5
(1 ((,.~ >: )      1 2 3    )} a055) -: _1  1  _ idmat 5 5
(1 ((,.~ >: )      1 2      )} a055) -: _1  1 _2 idmat 5 5
(1 ((,.~ >: )      1 2 3    )} a055) -: _1  1 __ idmat 5 5
(1 ((,.~ >: )    0 1 2      )} a055) -: _1 _2    idmat 5 5
                               a055  -: _1 _2  0 idmat 5 5
(1 ((,.~ >: )      1 2      )} a055) -: _1 _2  2 idmat 5 5
(1 ((,.~ >: )    0 1 2      )} a055) -: _1 _2  _ idmat 5 5
(1 ((,.~ >: )      1 2      )} a055) -: _1 _2 _2 idmat 5 5
(1 ((,.~ >: )    0 1 2      )} a055) -: _1 _2 __ idmat 5 5
(1 ((,.  >: ) i.         4  )} a055) -:  1       idmat 5 5
(1 ((,.  >: ) i.         4  )} a055) -:  1  0    idmat 5 5
                               a055  -:  1  0  0 idmat 5 5
(1 ((,.  >: )    0 1        )} a055) -:  1  0  2 idmat 5 5
(1 ((,.  >: ) i.         4  )} a055) -:  1  0  _ idmat 5 5
(1 ((,.  >: )    0 1        )} a055) -:  1  0 _2 idmat 5 5
(1 ((,.  >: ) i.         4  )} a055) -:  1  0 __ idmat 5 5
(1 ((,.  >: ) i.         4  )} a055) -:  1 _1    idmat 5 5
                               a055  -:  1 _1  0 idmat 5 5
(1 ((,.  >: )        2 3    )} a055) -:  1 _1  2 idmat 5 5
(1 ((,.  >: ) i.         4  )} a055) -:  1 _1  _ idmat 5 5
(1 ((,.  >: )        2 3    )} a055) -:  1 _1 _2 idmat 5 5
(1 ((,.  >: ) i.         4  )} a055) -:  1 _1 __ idmat 5 5
(1 ((,.  >: )      1 2 3    )} a055) -:  1  1    idmat 5 5
                               a055  -:  1  1  0 idmat 5 5
(1 ((,.  >: )      1 2      )} a055) -:  1  1  2 idmat 5 5
(1 ((,.  >: )      1 2 3    )} a055) -:  1  1  _ idmat 5 5
(1 ((,.  >: )      1 2      )} a055) -:  1  1 _2 idmat 5 5
(1 ((,.  >: )      1 2 3    )} a055) -:  1  1 __ idmat 5 5
(1 ((,.  >: )    0 1 2      )} a055) -:  1 _2    idmat 5 5
                               a055  -:  1 _2  0 idmat 5 5
(1 ((,.  >: )      1 2      )} a055) -:  1 _2  2 idmat 5 5
(1 ((,.  >: )    0 1 2      )} a055) -:  1 _2  _ idmat 5 5
(1 ((,.  >: )      1 2      )} a055) -:  1 _2 _2 idmat 5 5
(1 ((,.  >: )    0 1 2      )} a055) -:  1 _2 __ idmat 5 5
(1 ((,.  2&+) i.       3    )} a055) -:  2       idmat 5 5
(1 ((,.  2&+) i.       3    )} a055) -:  2  0    idmat 5 5
                               a055  -:  2  0  0 idmat 5 5
(1 ((,.  2&+)    0 1        )} a055) -:  2  0  2 idmat 5 5
(1 ((,.  2&+) i.       3    )} a055) -:  2  0  _ idmat 5 5
(1 ((,.  2&+)    0 1        )} a055) -:  2  0 _2 idmat 5 5
(1 ((,.  2&+) i.       3    )} a055) -:  2  0 __ idmat 5 5
(1 ((,.  2&+) i.       3    )} a055) -:  2 _1    idmat 5 5
                               a055  -:  2 _1  0 idmat 5 5
(1 ((,.  2&+)      1 2      )} a055) -:  2 _1  2 idmat 5 5
(1 ((,.  2&+) i.       3    )} a055) -:  2 _1  _ idmat 5 5
(1 ((,.  2&+)      1 2      )} a055) -:  2 _1 _2 idmat 5 5
(1 ((,.  2&+) i.       3    )} a055) -:  2 _1 __ idmat 5 5
(1 ((,.  2&+)      1 2      )} a055) -:  2  1    idmat 5 5
                               a055  -:  2  1  0 idmat 5 5
(1 ((,.  2&+)      1 2      )} a055) -:  2  1  2 idmat 5 5
(1 ((,.  2&+)      1 2      )} a055) -:  2  1  _ idmat 5 5
(1 ((,.  2&+)      1 2      )} a055) -:  2  1 _2 idmat 5 5
(1 ((,.  2&+)      1 2      )} a055) -:  2  1 __ idmat 5 5
(1 ((,.  2&+)    0 1        )} a055) -:  2 _2    idmat 5 5
                               a055  -:  2 _2  0 idmat 5 5
(1 ((,.  2&+)    0 1        )} a055) -:  2 _2  2 idmat 5 5
(1 ((,.  2&+)    0 1        )} a055) -:  2 _2  _ idmat 5 5
(1 ((,.  2&+)    0 1        )} a055) -:  2 _2 _2 idmat 5 5
(1 ((,.  2&+)    0 1        )} a055) -:  2 _2 __ idmat 5 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2       idmat 5 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2  0    idmat 5 5
                               a055  -: _2  0  0 idmat 5 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2  0  2 idmat 5 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2  0  _ idmat 5 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2  0 _2 idmat 5 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2  0 __ idmat 5 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2 _1    idmat 5 5
                               a055  -: _2 _1  0 idmat 5 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2 _1  2 idmat 5 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2 _1  _ idmat 5 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2 _1 _2 idmat 5 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2 _1 __ idmat 5 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2  1    idmat 5 5
                               a055  -: _2  1  0 idmat 5 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2  1  2 idmat 5 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2  1  _ idmat 5 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2  1 _2 idmat 5 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2  1 __ idmat 5 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2 _2    idmat 5 5
                               a055  -: _2 _2  0 idmat 5 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2 _2  2 idmat 5 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2 _2  _ idmat 5 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2 _2 _2 idmat 5 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2 _2 __ idmat 5 5

(1 ( ,.~      i.           5)} a055) -:          idmat 5
(1 ( ,.~      i.           5)} a055) -:  0       idmat 5
(1 ( ,.~      i.           5)} a055) -:  0  0    idmat 5
                               a055  -:  0  0  0 idmat 5
(1 ( ,.~         0 1        )} a055) -:  0  0  2 idmat 5
(1 ( ,.~      i.           5)} a055) -:  0  0  _ idmat 5
(1 ( ,.~         0 1        )} a055) -:  0  0 _2 idmat 5
(1 ( ,.~      i.           5)} a055) -:  0  0 __ idmat 5
(1 ( ,.~      i.           5)} a055) -:  0 _1    idmat 5
                               a055  -:  0 _1  0 idmat 5
(1 ( ,.~               3 4  )} a055) -:  0 _1  2 idmat 5
(1 ( ,.~      i.           5)} a055) -:  0 _1  _ idmat 5
(1 ( ,.~               3 4  )} a055) -:  0 _1 _2 idmat 5
(1 ( ,.~      i.           5)} a055) -:  0 _1 __ idmat 5
(1 ( ,.~           1 2 3 4  )} a055) -:  0  1    idmat 5
                               a055  -:  0  1  0 idmat 5
(1 ( ,.~           1 2      )} a055) -:  0  1  2 idmat 5
(1 ( ,.~           1 2 3 4  )} a055) -:  0  1  _ idmat 5
(1 ( ,.~           1 2      )} a055) -:  0  1 _2 idmat 5
(1 ( ,.~           1 2 3 4  )} a055) -:  0  1 __ idmat 5
(1 ( ,.~      i.         4  )} a055) -:  0 _2    idmat 5
                               a055  -:  0 _2  0 idmat 5
(1 ( ,.~             2 3    )} a055) -:  0 _2  2 idmat 5
(1 ( ,.~      i.         4  )} a055) -:  0 _2  _ idmat 5
(1 ( ,.~             2 3    )} a055) -:  0 _2 _2 idmat 5
(1 ( ,.~      i.         4  )} a055) -:  0 _2 __ idmat 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1       idmat 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1  0    idmat 5
                               a055  -: _1  0  0 idmat 5
(1 ((,.~ >: )    0 1        )} a055) -: _1  0  2 idmat 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1  0  _ idmat 5
(1 ((,.~ >: )    0 1        )} a055) -: _1  0 _2 idmat 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1  0 __ idmat 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1 _1    idmat 5
                               a055  -: _1 _1  0 idmat 5
(1 ((,.~ >: )        2 3    )} a055) -: _1 _1  2 idmat 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1 _1  _ idmat 5
(1 ((,.~ >: )        2 3    )} a055) -: _1 _1 _2 idmat 5
(1 ((,.~ >: ) i.         4  )} a055) -: _1 _1 __ idmat 5
(1 ((,.~ >: )      1 2 3    )} a055) -: _1  1    idmat 5
                               a055  -: _1  1  0 idmat 5
(1 ((,.~ >: )      1 2      )} a055) -: _1  1  2 idmat 5
(1 ((,.~ >: )      1 2 3    )} a055) -: _1  1  _ idmat 5
(1 ((,.~ >: )      1 2      )} a055) -: _1  1 _2 idmat 5
(1 ((,.~ >: )      1 2 3    )} a055) -: _1  1 __ idmat 5
(1 ((,.~ >: )    0 1 2      )} a055) -: _1 _2    idmat 5
                               a055  -: _1 _2  0 idmat 5
(1 ((,.~ >: )      1 2      )} a055) -: _1 _2  2 idmat 5
(1 ((,.~ >: )    0 1 2      )} a055) -: _1 _2  _ idmat 5
(1 ((,.~ >: )      1 2      )} a055) -: _1 _2 _2 idmat 5
(1 ((,.~ >: )    0 1 2      )} a055) -: _1 _2 __ idmat 5
(1 ((,.  >: ) i.         4  )} a055) -:  1       idmat 5
(1 ((,.  >: ) i.         4  )} a055) -:  1  0    idmat 5
                               a055  -:  1  0  0 idmat 5
(1 ((,.  >: )    0 1        )} a055) -:  1  0  2 idmat 5
(1 ((,.  >: ) i.         4  )} a055) -:  1  0  _ idmat 5
(1 ((,.  >: )    0 1        )} a055) -:  1  0 _2 idmat 5
(1 ((,.  >: ) i.         4  )} a055) -:  1  0 __ idmat 5
(1 ((,.  >: ) i.         4  )} a055) -:  1 _1    idmat 5
                               a055  -:  1 _1  0 idmat 5
(1 ((,.  >: )        2 3    )} a055) -:  1 _1  2 idmat 5
(1 ((,.  >: ) i.         4  )} a055) -:  1 _1  _ idmat 5
(1 ((,.  >: )        2 3    )} a055) -:  1 _1 _2 idmat 5
(1 ((,.  >: ) i.         4  )} a055) -:  1 _1 __ idmat 5
(1 ((,.  >: )      1 2 3    )} a055) -:  1  1    idmat 5
                               a055  -:  1  1  0 idmat 5
(1 ((,.  >: )      1 2      )} a055) -:  1  1  2 idmat 5
(1 ((,.  >: )      1 2 3    )} a055) -:  1  1  _ idmat 5
(1 ((,.  >: )      1 2      )} a055) -:  1  1 _2 idmat 5
(1 ((,.  >: )      1 2 3    )} a055) -:  1  1 __ idmat 5
(1 ((,.  >: )    0 1 2      )} a055) -:  1 _2    idmat 5
                               a055  -:  1 _2  0 idmat 5
(1 ((,.  >: )      1 2      )} a055) -:  1 _2  2 idmat 5
(1 ((,.  >: )    0 1 2      )} a055) -:  1 _2  _ idmat 5
(1 ((,.  >: )      1 2      )} a055) -:  1 _2 _2 idmat 5
(1 ((,.  >: )    0 1 2      )} a055) -:  1 _2 __ idmat 5
(1 ((,.  2&+) i.       3    )} a055) -:  2       idmat 5
(1 ((,.  2&+) i.       3    )} a055) -:  2  0    idmat 5
                               a055  -:  2  0  0 idmat 5
(1 ((,.  2&+)    0 1        )} a055) -:  2  0  2 idmat 5
(1 ((,.  2&+) i.       3    )} a055) -:  2  0  _ idmat 5
(1 ((,.  2&+)    0 1        )} a055) -:  2  0 _2 idmat 5
(1 ((,.  2&+) i.       3    )} a055) -:  2  0 __ idmat 5
(1 ((,.  2&+) i.       3    )} a055) -:  2 _1    idmat 5
                               a055  -:  2 _1  0 idmat 5
(1 ((,.  2&+)      1 2      )} a055) -:  2 _1  2 idmat 5
(1 ((,.  2&+) i.       3    )} a055) -:  2 _1  _ idmat 5
(1 ((,.  2&+)      1 2      )} a055) -:  2 _1 _2 idmat 5
(1 ((,.  2&+) i.       3    )} a055) -:  2 _1 __ idmat 5
(1 ((,.  2&+)      1 2      )} a055) -:  2  1    idmat 5
                               a055  -:  2  1  0 idmat 5
(1 ((,.  2&+)      1 2      )} a055) -:  2  1  2 idmat 5
(1 ((,.  2&+)      1 2      )} a055) -:  2  1  _ idmat 5
(1 ((,.  2&+)      1 2      )} a055) -:  2  1 _2 idmat 5
(1 ((,.  2&+)      1 2      )} a055) -:  2  1 __ idmat 5
(1 ((,.  2&+)    0 1        )} a055) -:  2 _2    idmat 5
                               a055  -:  2 _2  0 idmat 5
(1 ((,.  2&+)    0 1        )} a055) -:  2 _2  2 idmat 5
(1 ((,.  2&+)    0 1        )} a055) -:  2 _2  _ idmat 5
(1 ((,.  2&+)    0 1        )} a055) -:  2 _2 _2 idmat 5
(1 ((,.  2&+)    0 1        )} a055) -:  2 _2 __ idmat 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2       idmat 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2  0    idmat 5
                               a055  -: _2  0  0 idmat 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2  0  2 idmat 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2  0  _ idmat 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2  0 _2 idmat 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2  0 __ idmat 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2 _1    idmat 5
                               a055  -: _2 _1  0 idmat 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2 _1  2 idmat 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2 _1  _ idmat 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2 _1 _2 idmat 5
(1 ((,.~ 2&+) i.       3    )} a055) -: _2 _1 __ idmat 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2  1    idmat 5
                               a055  -: _2  1  0 idmat 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2  1  2 idmat 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2  1  _ idmat 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2  1 _2 idmat 5
(1 ((,.~ 2&+)      1 2      )} a055) -: _2  1 __ idmat 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2 _2    idmat 5
                               a055  -: _2 _2  0 idmat 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2 _2  2 idmat 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2 _2  _ idmat 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2 _2 _2 idmat 5
(1 ((,.~ 2&+)    0 1        )} a055) -: _2 _2 __ idmat 5

NB. diagmat
EMPTY -:     diagmat ''
EMPTY -: 0 0 diagmat ''
(1 1 $ 1) -:     diagmat 1
(1 1 $ 1) -: 0 0 diagmat 1
(3 3 $ 3 0 0 0 5 0 0 0 7) -:     diagmat 3 5 7
(3 3 $ 3 0 0 0 5 0 0 0 7) -: 0 0 diagmat 3 5 7
(3 4 $ 0 3 0 0 0 0 5 0 0 0 0 7) -:  1  0 diagmat 3 5 7
(4 3 $ 0 0 0 3 0 0 0 5 0 0 0 7) -: _1  0 diagmat 3 5 7
(4 3 $ 3 0 0 0 5 0 0 0 7 0 0 0) -:  0  1 diagmat 3 5 7
(3 4 $ 3 0 0 0 0 5 0 0 0 0 7 0) -:  0 _1 diagmat 3 5 7
(4 4 $ 0 3 0 0 0 0 5 0 0 0 0 7 0 0 0 0) -:  1  1 diagmat 3 5 7
(4 4 $ 0 0 0 0 3 0 0 0 0 5 0 0 0 0 7 0) -: _1 _1 diagmat 3 5 7
(3 5 $ 0 3 0 0 0 0 0 5 0 0 0 0 0 7 0) -:  1 _1 diagmat 3 5 7
(5 3 $ 0 0 0 3 0 0 0 5 0 0 0 7 0 0 0) -: _1  1 diagmat 3 5 7

NB. trl

(-:    trl) EMPTY
(-:  0&trl) EMPTY
(-:  1&trl) EMPTY
(-: _1&trl) EMPTY

(-:    trl) 1 1 $ 2
(-:  0&trl) 1 1 $ 2
(-:  1&trl) 1 1 $ 2
(-:  2&trl) 1 1 $ 2
EMPTY -: _1 trl 1 1 $ 2
EMPTY -: _2 trl 1 1 $ 2

(3 3 $ 0 0 0 5 6 0 10 11 12) -:   trl ai35
(3 3 $ 0 0 0 5 6 0 10 11 12) -: 0 trl ai35
(3 4 $ 0 1 0 0 5 6  7  0 10 11 12 13) -: 1 trl ai35
(3 5 $ 0 1 2 0 0 5  6  7  8  0 10 11 12 13 14) -: 2 trl ai35
(3 5 $ 0 1 2 3 0 5  6  7  8  9 10 11 12 13 14) -: 3 trl ai35
(-: 4&trl) ai35
(-: 5&trl) ai35
(-: 6&trl) ai35
(2 2 $ 5 0 10 11) -: _1 trl ai35
(1 1 $ 10) -: _2 trl ai35
EMPTY -: _3 trl ai35
EMPTY -: _4 trl ai35

(5 3 $ 0 0 0 3 4 0 6 7 8 9 10 11 12 13 14) -:   trl ai53
(5 3 $ 0 0 0 3 4 0 6 7 8 9 10 11 12 13 14) -: 0 trl ai53
(5 3 $ 0 1 0 3 4 5 6 7 8 9 10 11 12 13 14) -: 1 trl ai53
(5 3 $ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14) -: 2 trl ai53
(-: 2&trl) ai53
(-: 3&trl) ai53
(-: 4&trl) ai53
(4 3 $  3 0  0  6  7 0  9 10 11 12 13 14) -: _1 trl ai53
(3 3 $  6 0  0  9 10 0 12 13 14) -: _2 trl ai53
(2 2 $  9 0 12 13) -: _3 trl ai53
(1 1 $ 12) -: _4 trl ai53
EMPTY -: _5 trl ai53
EMPTY -: _6 trl ai53

(5 5 $ 0 0 0 0 0 5 6 0 0 0 10 11 12  0  0 15 16 17 18  0 20 21 22 23 24) -:   trl ai55
(5 5 $ 0 0 0 0 0 5 6 0 0 0 10 11 12  0  0 15 16 17 18  0 20 21 22 23 24) -: 0 trl ai55
(5 5 $ 0 1 0 0 0 5 6 7 0 0 10 11 12 13  0 15 16 17 18 19 20 21 22 23 24) -: 1 trl ai55
(5 5 $ 0 1 2 0 0 5 6 7 8 0 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24) -: 2 trl ai55
(5 5 $ 0 1 2 3 0 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24) -: 3 trl ai55
(-: 4&trl) ai55
(-: 5&trl) ai55
(-: 6&trl) ai55
(4 4 $  5 0  0  0 10 11  0  0 15 16 17 0 20 21 22 23) -: _1 trl ai55
(3 3 $ 10 0  0 15 16  0 20 21 22) -: _2 trl ai55
(2 2 $ 15 0 20 21) -: _3 trl ai55
(1 1 $ 20) -: _4 trl ai55
EMPTY -: _5 trl ai55
EMPTY -: _6 trl ai55

NB. tru

(-:    tru) EMPTY
(-:  0&tru) EMPTY
(-:  1&tru) EMPTY
(-: _1&tru) EMPTY

(-:    tru) 1 1 $ 2
(-:  0&tru) 1 1 $ 2
(-: _1&tru) 1 1 $ 2
(-: _2&tru) 1 1 $ 2
EMPTY -: 1 tru 1 1 $ 2
EMPTY -: 2 tru 1 1 $ 2

(3 5 $ 0 1 2 3 4 0 6 7 8 9  0  0 12 13 14) -:   tru ai35
(3 5 $ 0 1 2 3 4 0 6 7 8 9  0  0 12 13 14) -: 0 tru ai35
(3 4 $ 1 2 3 4 0 7 8 9 0 0 13 14) -: 1 tru ai35
(3 3 $ 2 3 4 0 8 9 0 0 14) -: 2 tru ai35
(2 2 $ 3 4 0 9) -: 3 tru ai35
(1 1 $ 4) -: 4 tru ai35
EMPTY -: 5 tru ai35
EMPTY -: 6 tru ai35
(3 5 $ 0 1 2 3 4 5 6 7 8 9 0 11 12 13 14) -: _1 tru ai35
(-: _2&tru) ai35
(-: _3&tru) ai35
(-: _4&tru) ai35

(3 3 $ 0 1 2 0 4 5 0 0 8) -:   tru ai53
(3 3 $ 0 1 2 0 4 5 0 0 8) -: 0 tru ai53
(2 2 $ 1 2 0 5) -: 1 tru ai53
(1 1 $ 2) -: 2 tru ai53
EMPTY -: 3 tru ai53
EMPTY -: 4 tru ai53
(4 3 $ 0 1 2 3 4 5 0 7 8 0  0 11) -: _1 tru ai53
(5 3 $ 0 1 2 3 4 5 6 7 8 0 10 11 0  0 14) -: _2 tru ai53
(5 3 $ 0 1 2 3 4 5 6 7 8 9 10 11 0 13 14) -: _3 tru ai53
(-: _4&tru) ai53
(-: _5&tru) ai53
(-: _6&tru) ai53

(5 5 $ 0 1 2 3 4 0 6 7  8 9  0  0 12 13 14  0 0 0 18 19 0 0 0 0 24) -:   tru ai55
(5 5 $ 0 1 2 3 4 0 6 7  8 9  0  0 12 13 14  0 0 0 18 19 0 0 0 0 24) -: 0 tru ai55
(4 4 $ 1 2 3 4 0 7 8 9  0 0 13 14  0  0  0 19) -: 1 tru ai55
(3 3 $ 2 3 4 0 8 9 0 0 14) -: 2 tru ai55
(2 2 $ 3 4 0 9) -: 3 tru ai55
(1 1 $ 4) -: 4 tru ai55
EMPTY -: 5 tru ai55
EMPTY -: 6 tru ai55
(5 5 $ 0 1 2 3 4 5 6 7 8 9  0 11 12 13 14  0  0 17 18 19 0  0  0 23 24) -: _1 tru ai55
(5 5 $ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14  0 16 17 18 19 0  0 22 23 24) -: _2 tru ai55
(5 5 $ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 0 21 22 23 24) -: _3 tru ai55
(-: _4&tru) ai55
(-: _5&tru) ai55
(-: _6&tru) ai55

NB. trl0

(-:    trl0) EMPTY
(-:  0&trl0) EMPTY
(-:  1&trl0) EMPTY
(-: _1&trl0) EMPTY

(1 1 $ 0) -:   trl0 1 1 $ 2
(1 1 $ 0) -: 0 trl0 1 1 $ 2
(-: 1&trl0) 1 1 $ 2
(-: 2&trl0) 1 1 $ 2
EMPTY -: _1 trl0 1 1 $ 2
EMPTY -: _2 trl0 1 1 $ 2

(3 3 $ 0 0 0 5 0 0 10 11  0) -:   trl0 ai35
(3 3 $ 0 0 0 5 0 0 10 11  0) -: 0 trl0 ai35
(3 4 $ 0 0 0 0 5 6  0  0 10 11 12 0) -: 1 trl0 ai35
(3 5 $ 0 1 0 0 0 5  6  7  0  0 10 11 12 13  0) -: 2 trl0 ai35
(3 5 $ 0 1 2 0 0 5  6  7  8  0 10 11 12 13 14) -: 3 trl0 ai35
(3 5 $ 0 1 2 3 0 5  6  7  8  9 10 11 12 13 14) -: 4 trl0 ai35
(-: 5&trl0) ai35
(-: 6&trl0) ai35
(2 2 $ 0 0 10 0) -: _1 trl0 ai35
(1 1 $ 0) -: _2 trl0 ai35
EMPTY -: _3 trl0 ai35
EMPTY -: _4 trl0 ai35

(5 3 $ 0 0 0 3 0 0 6 7 0 9 10 11 12 13 14) -:   trl0 ai53
(5 3 $ 0 0 0 3 0 0 6 7 0 9 10 11 12 13 14) -: 0 trl0 ai53
(5 3 $ 0 0 0 3 4 0 6 7 8 9 10 11 12 13 14) -: 1 trl0 ai53
(5 3 $ 0 1 0 3 4 5 6 7 8 9 10 11 12 13 14) -: 2 trl0 ai53
(-: 3&trl0) ai53
(-: 4&trl0) ai53
(4 3 $ 0 0  0 6 0 0  9 10 0 12 13 14) -: _1 trl0 ai53
(3 3 $ 0 0  0 9 0 0 12 13 0) -: _2 trl0 ai53
(2 2 $ 0 0 12 0) -: _3 trl0 ai53
(1 1 $ 0) -: _4 trl0 ai53
EMPTY -: _5 trl0 ai53
EMPTY -: _6 trl0 ai53

(5 5 $ 0 0 0 0 0 5 0 0 0 0 10 11  0  0  0 15 16 17  0  0 20 21 22 23  0) -:   trl0 ai55
(5 5 $ 0 0 0 0 0 5 0 0 0 0 10 11  0  0  0 15 16 17  0  0 20 21 22 23  0) -: 0 trl0 ai55
(5 5 $ 0 0 0 0 0 5 6 0 0 0 10 11 12  0  0 15 16 17 18  0 20 21 22 23 24) -: 1 trl0 ai55
(5 5 $ 0 1 0 0 0 5 6 7 0 0 10 11 12 13  0 15 16 17 18 19 20 21 22 23 24) -: 2 trl0 ai55
(5 5 $ 0 1 2 0 0 5 6 7 8 0 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24) -: 3 trl0 ai55
(5 5 $ 0 1 2 3 0 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24) -: 4 trl0 ai55
(-: 5&trl0) ai55
(-: 6&trl0) ai55
(4 4 $ 0 0  0  0 10 0  0  0 15 16 0 0 20 21 22 0) -: _1 trl0 ai55
(3 3 $ 0 0  0 15  0 0 20 21  0) -: _2 trl0 ai55
(2 2 $ 0 0 20  0) -: _3 trl0 ai55
(1 1 $ 0) -: _4 trl0 ai55
EMPTY -: _5 trl0 ai55
EMPTY -: _6 trl0 ai55

NB. tru0

(-:    tru0) EMPTY
(-:  0&tru0) EMPTY
(-:  1&tru0) EMPTY
(-: _1&tru0) EMPTY

(1 1 $ 0) -:   tru0 1 1 $ 2
(1 1 $ 0) -: 0 tru0 1 1 $ 2
EMPTY -: 1 tru0 1 1 $ 2
EMPTY -: 2 tru0 1 1 $ 2
(-: _1&tru0) 1 1 $ 2
(-: _2&tru0) 1 1 $ 2

(3 5 $ 0 1 2 3 4 0 0 7 8 9 0  0 0 13 14) -:   tru0 ai35
(3 5 $ 0 1 2 3 4 0 0 7 8 9 0  0 0 13 14) -: 0 tru0 ai35
(3 4 $ 0 2 3 4 0 0 8 9 0 0 0 14) -: 1 tru0 ai35
(3 3 $ 0 3 4 0 0 9 0 0 0) -: 2 tru0 ai35
(2 2 $ 0 4 0 0) -: 3 tru0 ai35
(1 1 $ 0) -: 4 tru0 ai35
EMPTY -: 5 tru0 ai35
EMPTY -: 6 tru0 ai35
(3 5 $ 0 1 2 3 4 0 6 7 8 9 0  0 12 13 14) -: _1 tru0 ai35
(3 5 $ 0 1 2 3 4 5 6 7 8 9 0 11 12 13 14) -: _2 tru0 ai35
(-: _3&tru0) ai35
(-: _4&tru0) ai35

(3 3 $ 0 1 2 0 0 5 0 0 0) -:   tru0 ai53
(3 3 $ 0 1 2 0 0 5 0 0 0) -: 0 tru0 ai53
(2 2 $ 0 2 0 0) -: 1 tru0 ai53
(1 1 $ 0) -: 2 tru0 ai53
EMPTY -: 3 tru0 ai53
EMPTY -: 4 tru0 ai53
(4 3 $ 0 1 2 0 4 5 0 0 8 0  0  0) -: _1 tru0 ai53
(5 3 $ 0 1 2 3 4 5 0 7 8 0  0 11 0  0  0) -: _2 tru0 ai53
(5 3 $ 0 1 2 3 4 5 6 7 8 0 10 11 0  0 14) -: _3 tru0 ai53
(5 3 $ 0 1 2 3 4 5 6 7 8 9 10 11 0 13 14) -: _4 tru0 ai53
(-: _5&tru0) ai53
(-: _6&tru0) ai53

(5 5 $ 0 1 2 3 4 0 0 7 8 9 0  0 0 13 14 0 0 0 0 19 0 0 0 0 0) -:   tru0 ai55
(5 5 $ 0 1 2 3 4 0 0 7 8 9 0  0 0 13 14 0 0 0 0 19 0 0 0 0 0) -: 0 tru0 ai55
(4 4 $ 0 2 3 4 0 0 8 9 0 0 0 14 0  0  0 0) -: 1 tru0 ai55
(3 3 $ 0 3 4 0 0 9 0 0 0) -: 2 tru0 ai55
(2 2 $ 0 4 0 0) -: 3 tru0 ai55
(1 1 $ 0) -: 4 tru0 ai55
EMPTY -: 5 tru0 ai55
EMPTY -: 6 tru0 ai55
(5 5 $ 0 1 2 3 4 0 6 7 8 9  0  0 12 13 14  0  0  0 18 19 0  0  0  0 24) -: _1 tru0 ai55
(5 5 $ 0 1 2 3 4 5 6 7 8 9  0 11 12 13 14  0  0 17 18 19 0  0  0 23 24) -: _2 tru0 ai55
(5 5 $ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14  0 16 17 18 19 0  0 22 23 24) -: _3 tru0 ai55
(5 5 $ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 0 21 22 23 24) -: _4 tru0 ai55
(-: _5&tru0) ai55
(-: _6&tru0) ai55

NB. trl1

(-:    trl1) EMPTY
(-:  0&trl1) EMPTY
(-:  1&trl1) EMPTY
(-: _1&trl1) EMPTY

(1 1 $ 1) -:   trl1 1 1 $ 2
(1 1 $ 1) -: 0 trl1 1 1 $ 2
(-: 1&trl1) 1 1 $ 2
(-: 2&trl1) 1 1 $ 2
EMPTY -: _1 trl1 1 1 $ 2
EMPTY -: _2 trl1 1 1 $ 2

(3 3 $ 1 0 0 5 1 0 10 11  1) -:   trl1 ai35
(3 3 $ 1 0 0 5 1 0 10 11  1) -: 0 trl1 ai35
(3 4 $ 0 1 0 0 5 6  1  0 10 11 12  1) -: 1 trl1 ai35
(3 5 $ 0 1 1 0 0 5  6  7  1  0 10 11 12 13  1) -: 2 trl1 ai35
(3 5 $ 0 1 2 1 0 5  6  7  8  1 10 11 12 13 14) -: 3 trl1 ai35
(3 5 $ 0 1 2 3 1 5  6  7  8  9 10 11 12 13 14) -: 4 trl1 ai35
(-: 5&trl1) ai35
(-: 6&trl1) ai35
(2 2 $ 1 0 10 1) -: _1 trl1 ai35
(1 1 $ 1) -: _2 trl1 ai35
EMPTY -: _3 trl1 ai35
EMPTY -: _4 trl1 ai35

(5 3 $ 1 0 0 3 1 0 6 7 1 9 10 11 12 13 14) -:   trl1 ai53
(5 3 $ 1 0 0 3 1 0 6 7 1 9 10 11 12 13 14) -: 0 trl1 ai53
(5 3 $ 0 1 0 3 4 1 6 7 8 9 10 11 12 13 14) -: 1 trl1 ai53
(5 3 $ 0 1 1 3 4 5 6 7 8 9 10 11 12 13 14) -: 2 trl1 ai53
(-: 3&trl1) ai53
(-: 4&trl1) ai53
(4 3 $ 1 0  0 6 1 0  9 10 1 12 13 14) -: _1 trl1 ai53
(3 3 $ 1 0  0 9 1 0 12 13 1) -: _2 trl1 ai53
(2 2 $ 1 0 12 1) -: _3 trl1 ai53
(1 1 $ 1) -: _4 trl1 ai53
EMPTY -: _5 trl1 ai53
EMPTY -: _6 trl1 ai53

(5 5 $ 1 0 0 0 0 5 1 0 0 0 10 11  1  0  0 15 16 17  1  0 20 21 22 23  1) -:   trl1 ai55
(5 5 $ 1 0 0 0 0 5 1 0 0 0 10 11  1  0  0 15 16 17  1  0 20 21 22 23  1) -: 0 trl1 ai55
(5 5 $ 0 1 0 0 0 5 6 1 0 0 10 11 12  1  0 15 16 17 18  1 20 21 22 23 24) -: 1 trl1 ai55
(5 5 $ 0 1 1 0 0 5 6 7 1 0 10 11 12 13  1 15 16 17 18 19 20 21 22 23 24) -: 2 trl1 ai55
(5 5 $ 0 1 2 1 0 5 6 7 8 1 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24) -: 3 trl1 ai55
(5 5 $ 0 1 2 3 1 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24) -: 4 trl1 ai55
(-: 5&trl1) ai55
(-: 6&trl1) ai55
(4 4 $ 1 0  0  0 10 1  0  0 15 16 1 0 20 21 22 1) -: _1 trl1 ai55
(3 3 $ 1 0  0 15  1 0 20 21  1) -: _2 trl1 ai55
(2 2 $ 1 0 20  1) -: _3 trl1 ai55
(1 1 $ 1) -: _4 trl1 ai55
EMPTY -: _5 trl1 ai55
EMPTY -: _6 trl1 ai55

NB. tru1

(-:    tru1) EMPTY
(-:  0&tru1) EMPTY
(-:  1&tru1) EMPTY
(-: _1&tru1) EMPTY

(1 1 $ 1) -:   tru1 1 1 $ 2
(1 1 $ 1) -: 0 tru1 1 1 $ 2
EMPTY -: 1 tru1 1 1 $ 2
EMPTY -: 2 tru1 1 1 $ 2
(-: _1&tru1) 1 1 $ 2
(-: _2&tru1) 1 1 $ 2

(3 5 $ 1 1 2 3 4 0 1 7 8 9 0  0 1 13 14) -:   tru1 ai35
(3 5 $ 1 1 2 3 4 0 1 7 8 9 0  0 1 13 14) -: 0 tru1 ai35
(3 4 $ 1 2 3 4 0 1 8 9 0 0 1 14) -: 1 tru1 ai35
(3 3 $ 1 3 4 0 1 9 0 0 1) -: 2 tru1 ai35
(2 2 $ 1 4 0 1) -: 3 tru1 ai35
(1 1 $ 1) -: 4 tru1 ai35
EMPTY -: 5 tru1 ai35
EMPTY -: 6 tru1 ai35
(3 5 $ 0 1 2 3 4 1 6 7 8 9 0  1 12 13 14) -: _1 tru1 ai35
(3 5 $ 0 1 2 3 4 5 6 7 8 9 1 11 12 13 14) -: _2 tru1 ai35
(-: _3&tru1) ai35
(-: _4&tru1) ai35

(3 3 $ 1 1 2 0 1 5 0 0 1) -:   tru1 ai53
(3 3 $ 1 1 2 0 1 5 0 0 1) -: 0 tru1 ai53
(2 2 $ 1 2 0 1) -: 1 tru1 ai53
(1 1 $ 1) -: 2 tru1 ai53
EMPTY -: 3 tru1 ai53
EMPTY -: 4 tru1 ai53
(4 3 $ 0 1 2 1 4 5 0 1 8 0  0  1) -: _1 tru1 ai53
(5 3 $ 0 1 2 3 4 5 1 7 8 0  1 11 0  0  1) -: _2 tru1 ai53
(5 3 $ 0 1 2 3 4 5 6 7 8 1 10 11 0  1 14) -: _3 tru1 ai53
(5 3 $ 0 1 2 3 4 5 6 7 8 9 10 11 1 13 14) -: _4 tru1 ai53
(-: _5&tru1) ai53
(-: _6&tru1) ai53

(5 5 $ 1 1 2 3 4 0 1 7 8 9 0  0 1 13 14 0 0 0 1 19 0 0 0 0 1) -:   tru1 ai55
(5 5 $ 1 1 2 3 4 0 1 7 8 9 0  0 1 13 14 0 0 0 1 19 0 0 0 0 1) -: 0 tru1 ai55
(4 4 $ 1 2 3 4 0 1 8 9 0 0 1 14 0  0  0 1) -: 1 tru1 ai55
(3 3 $ 1 3 4 0 1 9 0 0 1) -: 2 tru1 ai55
(2 2 $ 1 4 0 1) -: 3 tru1 ai55
(1 1 $ 1) -: 4 tru1 ai55
EMPTY -: 5 tru1 ai55
EMPTY -: 6 tru1 ai55
(5 5 $ 0 1 2 3 4 1 6 7 8 9  0  1 12 13 14  0  0  1 18 19 0  0  0  1 24) -: _1 tru1 ai55
(5 5 $ 0 1 2 3 4 5 6 7 8 9  1 11 12 13 14  0  1 17 18 19 0  0  1 23 24) -: _2 tru1 ai55
(5 5 $ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14  1 16 17 18 19 0  1 22 23 24) -: _3 tru1 ai55
(5 5 $ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 1 21 22 23 24) -: _4 tru1 ai55
(-: _5&tru1) ai55
(-: _6&tru1) ai55

NB. sy4gel
(-: sy4gel) EMPTY
(-: sy4gel) 1 1 $ 2
(1 1 $ 2j3) -: sy4gel 1 1 $ 2j3
(4 4 $ 0    4     8    12    4    5    9    13    8    9    10    14    12    13    14    15   ) -: sy4gel ai44
(4 4 $ 0j16 4j20  8j24 12j28 4j20 5j21 9j25 13j29 8j24 9j25 10j26 14j30 12j28 13j29 14j30 15j31) -: sy4gel ac44

NB. sy4geu
(-: sy4geu) EMPTY
(-: sy4geu) 1 1 $ 2
(1 1 $ 2j3) -: sy4geu 1 1 $ 2j3
(4 4 $ 0    1    2    3    1    5    6    7    2    6    10    11    3    7    11    15   ) -: sy4geu ai44
(4 4 $ 0j16 1j17 2j18 3j19 1j17 5j21 6j22 7j23 2j18 6j22 10j26 11j27 3j19 7j23 11j27 15j31) -: sy4geu ac44

NB. he4gel
(-: he4gel) EMPTY
(-: he4gel) 1 1 $ 2
(1 1 $ 2) -: he4gel 1 1 $ 2j3
(4 4 $ 0 4     8     12     4    5 9     13     8    9    10 14     12    13    14    15) -: he4gel ai44
(4 4 $ 0 4j_20 8j_24 12j_28 4j20 5 9j_25 13j_29 8j24 9j25 10 14j_30 12j28 13j29 14j30 15) -: he4gel ac44

NB. he4geu
(-: he4geu) EMPTY
(-: he4geu) 1 1 $ 2
(1 1 $ 2) -: he4geu 1 1 $ 2j3
(4 4 $ 0 1    2    3    1     5 6    7    2     6     10 11    3     7     11     15) -: he4geu ai44
(4 4 $ 0 1j17 2j18 3j19 1j_17 5 6j22 7j23 2j_18 6j_22 10 11j27 3j_19 7j_23 11j_27 15) -: he4geu ac44

NB. ss4gel
(-: ss4gel) EMPTY
(1 1 $ 0) -: ss4gel 1 1 $ 2
(1 1 $ 0) -: ss4gel 1 1 $ 2j3
(4 4 $ 0 _4     _8     _12     4    0 _9     _13     8    9    0 _14     12    13    14    0) -: ss4gel ai44
(4 4 $ 0 _4j_20 _8j_24 _12j_28 4j20 0 _9j_25 _13j_29 8j24 9j25 0 _14j_30 12j28 13j29 14j30 0) -: ss4gel ac44

NB. ss4geu
(-: ss4geu) EMPTY
(1 1 $ 0) -: ss4geu 1 1 $ 2
(1 1 $ 0) -: ss4geu 1 1 $ 2j3
(4 4 $ 0  1    2    3    _1     0 6    7    _2     _6     0 11    _3     _7     _11     0) -: ss4geu ai44
(4 4 $ 0  1j17 2j18 3j19 _1j_17 0 6j22 7j23 _2j_18 _6j_22 0 11j27 _3j_19 _7j_23 _11j_27 0) -: ss4geu ac44

NB. sh4gel
(-: sh4gel) EMPTY
(1 1 $ 0) -: sh4gel 1 1 $ 2
(1 1 $ 0) -: sh4gel 1 1 $ 2j3
(4 4 $ 0 _4    _8    _12    4    0 _9    _13    8    9    0 _14    12    13    14    0) -: sh4gel ai44
(4 4 $ 0 _4j20 _8j24 _12j28 4j20 0 _9j25 _13j29 8j24 9j25 0 _14j30 12j28 13j29 14j30 0) -: sh4gel ac44

NB. sh4geu
(-: sh4geu) EMPTY
(1 1 $ 0) -: sh4geu 1 1 $ 2
(1 1 $ 0) -: sh4geu 1 1 $ 2j3
(4 4 $ 0 1    2    3    _1    0 6    7    _2    _6    0 11    _3    _7    _11    0) -: sh4geu ai44
(4 4 $ 0 1j17 2j18 3j19 _1j17 0 6j22 7j23 _2j18 _6j22 0 11j27 _3j19 _7j23 _11j27 0) -: sh4geu ac44

NB. lxsuy
(-: lxsuy) EMPTY
EMPTY -: 42    lxsuy EMPTY
EMPTY -: EMPTY lxsuy 42
(1 1 $ 1) -: (1 1 $ 1) lxsuy (1 1 $ 2)
(1 1 $ 1) -:        1  lxsuy (1 1 $ 2)
(1 1 $ 1) -: (1 1 $ 1) lxsuy        2
(4 5 $ 1 2 2 2 2 1 1 2 2 2 1 1 1 2 2 1 1 1 1 2) -: a145 lxsuy a245
(4 5 $ 1 2 2 2 2 1 1 2 2 2 1 1 1 2 2 1 1 1 1 2) -: 1    lxsuy a245
(4 5 $ 1 2 2 2 2 1 1 2 2 2 1 1 1 2 2 1 1 1 1 2) -: a145 lxsuy 2

NB. slxuy
(-: slxuy) EMPTY
EMPTY -: 42    slxuy EMPTY
EMPTY -: EMPTY slxuy 42
(1 1 $ 2) -: (1 1 $ 1) slxuy (1 1 $ 2)
(1 1 $ 2) -:        1  slxuy (1 1 $ 2)
(1 1 $ 2) -: (1 1 $ 1) slxuy        2
(4 5 $ 2 2 2 2 2 1 2 2 2 2 1 1 2 2 2 1 1 1 2 2) -: a145 slxuy a245
(4 5 $ 2 2 2 2 2 1 2 2 2 2 1 1 2 2 2 1 1 1 2 2) -: 1    slxuy a245
(4 5 $ 2 2 2 2 2 1 2 2 2 2 1 1 2 2 2 1 1 1 2 2) -: a145 slxuy 2

NB. suxly
(-: suxly) EMPTY
EMPTY -: 42    suxly EMPTY
EMPTY -: EMPTY suxly 42
(1 1 $ 2) -: (1 1 $ 1) suxly (1 1 $ 2)
(1 1 $ 2) -:        1  suxly (1 1 $ 2)
(1 1 $ 2) -: (1 1 $ 1) suxly        2
(4 5 $ 2 1 1 1 1 2 2 1 1 1 2 2 2 1 1 2 2 2 2 1) -: a145 suxly a245
(4 5 $ 2 1 1 1 1 2 2 1 1 1 2 2 2 1 1 2 2 2 2 1) -: 1    suxly a245
(4 5 $ 2 1 1 1 1 2 2 1 1 1 2 2 2 1 1 2 2 2 2 1) -: a145 suxly 2

NB. uxsly
(-: uxsly) EMPTY
EMPTY -: 42    uxsly EMPTY
EMPTY -: EMPTY uxsly 42
(1 1 $ 1) -: (1 1 $ 1) uxsly (1 1 $ 2)
(1 1 $ 1) -:        1  uxsly (1 1 $ 2)
(1 1 $ 1) -: (1 1 $ 1) uxsly        2
(4 5 $ 1 1 1 1 1 2 1 1 1 1 2 2 1 1 1 2 2 2 1 1) -: a145 uxsly a245
(4 5 $ 1 1 1 1 1 2 1 1 1 1 2 2 1 1 1 2 2 2 1 1) -: 1    uxsly a245
(4 5 $ 1 1 1 1 1 2 1 1 1 1 2 2 1 1 1 2 2 2 1 1) -: a145 uxsly 2

NB. slxsuy
(-: 3 slxsuy) EMPTY
EMPTY -: 42    (3 slxsuy) EMPTY
EMPTY -: EMPTY  3 slxsuy  42
(1 1 $ 3) -: (1 1 $ 1)  3 slxsuy  (1 1 $ 2)
(1 1 $ 3) -:        1  (3 slxsuy) (1 1 $ 2)
(1 1 $ 3) -: (1 1 $ 1)  3 slxsuy         2
(4 5 $ 3 2 2 2 2 1 3 2 2 2 1 1 3 2 2 1 1 1 3 2) -: a145  3 slxsuy  a245
(4 5 $ 3 2 2 2 2 1 3 2 2 2 1 1 3 2 2 1 1 1 3 2) -: 1    (3 slxsuy) a245
(4 5 $ 3 2 2 2 2 1 3 2 2 2 1 1 3 2 2 1 1 1 3 2) -: a145  3 slxsuy  2

NB. suxsly
(-: 3 suxsly) EMPTY
EMPTY -: 42    (3 suxsly) EMPTY
EMPTY -: EMPTY  3 suxsly  42
(1 1 $ 3) -: (1 1 $ 1)  3 suxsly  (1 1 $ 2)
(1 1 $ 3) -:        1  (3 suxsly) (1 1 $ 2)
(1 1 $ 3) -: (1 1 $ 1)  3 suxsly         2
(4 5 $ 3 1 1 1 1 2 3 1 1 1 2 2 3 1 1 2 2 2 3 1) -: a145  3 suxsly  a245
(4 5 $ 3 1 1 1 1 2 3 1 1 1 2 2 3 1 1 2 2 2 3 1) -: 1    (3 suxsly) a245
(4 5 $ 3 1 1 1 1 2 3 1 1 1 2 2 3 1 1 2 2 2 3 1) -: a145  3 suxsly  2

NB. po
(-: po) EMPTY
(1 1 $ 1) -: po 1 1 $  1
(1 1 $ 1) -: po 1 1 $ _1
(3 3 $ 122 62 15 62 75 59 15 59 81) -: po (3 3 $ _5 _4 9 _5 _7 1 _7 _4 _4)
(3 3 $ 279 18j4 _121j_98 18j_4 160 7j_100 _121j98 7j100 233) -: po (3 3 $ 5j8 _6j_3 _9j8 _6j_7 0j_1 _5j7 _2j_4 8j7 _6j_8)
