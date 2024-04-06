NB. Verify util verbs
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
NB. Local definitions

delimiters=. LF , ' '
string=. 'foo bar  baz' , LF , 'qux' , LF2 , 'quux' , LF , ' corge ' , LF , 'flob'

NB. =========================================================
NB. Verification suite

NB. verify itself
0 0 1 0 0 -: isnan  _1 _0.0 _. 0.0 1
0 0 0 1 0 -: ispos0 _1 _0.0 _. 0.0 1
0 1 0 0 0 -: isneg0 _1 _0.0 _. 0.0 1

NB. max
 0 -:  max ''
__ -:  max _.   _.
 _ -:  max _.    _ __
 1 -:  max  0    1 _.
 1 -:  max  0    1 _. __
ispos0 max  0.0 _1
isneg0 max _0.0 _1

NB. negneg

       isnan _.   negneg _.
       isnan _.   negneg __
       isnan _.   negneg _1
       isnan _.   negneg _0.0
       isnan _.   negneg  0
       isnan _.   negneg  1
       isnan _.   negneg  _
1 1 -: isnan _.   negneg _. _.
1 1 -: isnan _.   negneg _.  2
1 1 -: isnan _.   negneg __  2
       isnan __   negneg _.
       isnan _1   negneg _.
       isnan _0.0 negneg _.
       isnan  0   negneg _.
       isnan  1   negneg _.
       isnan  _   negneg _.
1 1 -: isnan __   negneg _. _.
1 1 -: isnan _1   negneg _. _.
1 1 -: isnan _0.0 negneg _. _.
1 1 -: isnan  0   negneg _. _.
1 1 -: isnan  1   negneg _. _.
1 1 -: isnan  _   negneg _. _.
1 1 -: isnan __   negneg _.  2
1 1 -: isnan _1   negneg _.  2
1 1 -: isnan _0.0 negneg _.  2
1 1 -: isnan  0   negneg _.  2
1 1 -: isnan  1   negneg _.  2
1 1 -: isnan  _   negneg _.  2

'' -: __   negneg ''
'' -: _1   negneg ''
'' -: _0.0 negneg ''
'' -:  0   negneg ''
'' -:  1   negneg ''
'' -:  _   negneg ''

 _ -: __   negneg __
 _ -: _1   negneg __
 _ -: _0.0 negneg __
__ -:  0   negneg __
__ -:  1   negneg __
__ -:  _   negneg __

 1 -: __   negneg _1
 1 -: _1   negneg _1
 1 -: _0.0 negneg _1
_1 -:  0   negneg _1
_1 -:  1   negneg _1
_1 -:  _   negneg _1

ispos0 __   negneg _0.0
ispos0 _1   negneg _0.0
ispos0 _0.0 negneg _0.0
isneg0  0   negneg _0.0
isneg0  1   negneg _0.0
isneg0  _   negneg _0.0

isneg0 __   negneg 0.0
isneg0 _1   negneg 0.0
isneg0 _0.0 negneg 0.0
ispos0  0   negneg 0.0
ispos0  1   negneg 0.0
ispos0  _   negneg 0.0

_1 -: __   negneg 1
_1 -: _1   negneg 1
_1 -: _0.0 negneg 1
 1 -:  0   negneg 1
 1 -:  1   negneg 1
 1 -:  _   negneg 1

__ -: __   negneg _
__ -: _1   negneg _
__ -: _0.0 negneg _
 _ -:  0   negneg _
 _ -:  1   negneg _
 _ -:  _   negneg _

( _  1 0 0 _1 __&-: *. (0 0 0 1 0 0 -: isneg0) *. 0 0 1 0 0 0 -: ispos0) _1 negneg __ _1 _0.0 0 1 _
(__ _1 0 0  1  _&-: *. (0 0 1 0 0 0 -: isneg0) *. 0 0 0 1 0 0 -: ispos0)  1 negneg __ _1 _0.0 0 1 _

NB. negpos

       isnan _.   negpos _.
       isnan _.   negpos __
       isnan _.   negpos _1
       isnan _.   negpos _0.0
       isnan _.   negpos  0
       isnan _.   negpos  1
       isnan _.   negpos  _
1 1 -: isnan _.   negpos _. _.
1 1 -: isnan _.   negpos _.  2
1 1 -: isnan _.   negpos __  2
       isnan __   negpos _.
       isnan _1   negpos _.
       isnan _0.0 negpos _.
       isnan  0   negpos _.
       isnan  1   negpos _.
       isnan  _   negpos _.
1 1 -: isnan __   negpos _. _.
1 1 -: isnan _1   negpos _. _.
1 1 -: isnan _0.0 negpos _. _.
1 1 -: isnan  0   negpos _. _.
1 1 -: isnan  1   negpos _. _.
1 1 -: isnan  _   negpos _. _.
1 1 -: isnan __   negpos _.  2
1 1 -: isnan _1   negpos _.  2
1 1 -: isnan _0.0 negpos _.  2
1 1 -: isnan  0   negpos _.  2
1 1 -: isnan  1   negpos _.  2
1 1 -: isnan  _   negpos _.  2

'' -: __   negpos ''
'' -: _1   negpos ''
'' -: _0.0 negpos ''
'' -:  0   negpos ''
'' -:  1   negpos ''
'' -:  _   negpos ''

__ -: __   negpos __
__ -: _1   negpos __
__ -: _0.0 negpos __
 _ -:  0   negpos __
 _ -:  1   negpos __
 _ -:  _   negpos __

_1 -: __   negpos _1
_1 -: _1   negpos _1
_1 -: _0.0 negpos _1
 1 -:  0   negpos _1
 1 -:  1   negpos _1
 1 -:  _   negpos _1

isneg0 __   negpos _0.0
isneg0 _1   negpos _0.0
isneg0 _0.0 negpos _0.0
ispos0  0   negpos _0.0
ispos0  1   negpos _0.0
ispos0  _   negpos _0.0

ispos0 __   negpos 0.0
ispos0 _1   negpos 0.0
ispos0 _0.0 negpos 0.0
isneg0  0   negpos 0.0
isneg0  1   negpos 0.0
isneg0  _   negpos 0.0

 1 -: __   negpos 1
 1 -: _1   negpos 1
 1 -: _0.0 negpos 1
_1 -:  0   negpos 1
_1 -:  1   negpos 1
_1 -:  _   negpos 1

 _ -: __   negpos _
 _ -: _1   negpos _
 _ -: _0.0 negpos _
__ -:  0   negpos _
__ -:  1   negpos _
__ -:  _   negpos _

(__ _1 0 0  1  _&-: *. (0 0 1 0 0 0 -: isneg0) *. 0 0 0 1 0 0 -: ispos0) _1 negpos __ _1 _0.0 0 1 _
( _  1 0 0 _1 __&-: *. (0 0 0 1 0 0 -: isneg0) *. 0 0 1 0 0 0 -: ispos0)  1 negpos __ _1 _0.0 0 1 _

NB. copysign

       isnan _.   copysign _.
       isnan _.   copysign __
       isnan _.   copysign _1
       isnan _.   copysign _0.0
       isnan _.   copysign  0
       isnan _.   copysign  1
       isnan _.   copysign  _
1 1 -: isnan _.   copysign _. _.
1 1 -: isnan _.   copysign _.  2
1 1 -: isnan _.   copysign __  2
       isnan __   copysign _.
       isnan _1   copysign _.
       isnan _0.0 copysign _.
       isnan  0   copysign _.
       isnan  1   copysign _.
       isnan  _   copysign _.
1 1 -: isnan __   copysign _. _.
1 1 -: isnan _1   copysign _. _.
1 1 -: isnan _0.0 copysign _. _.
1 1 -: isnan  0   copysign _. _.
1 1 -: isnan  1   copysign _. _.
1 1 -: isnan  _   copysign _. _.
1 1 -: isnan __   copysign _.  2
1 1 -: isnan _1   copysign _.  2
1 1 -: isnan _0.0 copysign _.  2
1 1 -: isnan  0   copysign _.  2
1 1 -: isnan  1   copysign _.  2
1 1 -: isnan  _   copysign _.  2

'' -: __   copysign ''
'' -: _1   copysign ''
'' -: _0.0 copysign ''
'' -:  0   copysign ''
'' -:  1   copysign ''
'' -:  _   copysign ''

__ -: __   copysign __
__ -: _1   copysign __
__ -: _0.0 copysign __
 _ -:  0   copysign __
 _ -:  1   copysign __
 _ -:  _   copysign __

_1 -: __   copysign _1
_1 -: _1   copysign _1
_1 -: _0.0 copysign _1
 1 -:  0   copysign _1
 1 -:  1   copysign _1
 1 -:  _   copysign _1

isneg0 __   copysign _0.0
isneg0 _1   copysign _0.0
isneg0 _0.0 copysign _0.0
ispos0  0   copysign _0.0
ispos0  1   copysign _0.0
ispos0  _   copysign _0.0

isneg0 __   copysign  0.0
isneg0 _1   copysign  0.0
isneg0 _0.0 copysign  0.0
ispos0  0   copysign  0.0
ispos0  1   copysign  0.0
ispos0  _   copysign  0.0

_1 -: __   copysign  1
_1 -: _1   copysign  1
_1 -: _0.0 copysign  1
 1 -:  0   copysign  1
 1 -:  1   copysign  1
 1 -:  _   copysign  1

__ -: __   copysign  _
__ -: _1   copysign  _
__ -: _0.0 copysign  _
 _ -:  0   copysign  _
 _ -:  1   copysign  _
 _ -:  _   copysign  _

(__ _1 0 0 _1 __&-: *. (0 0 1 1 0 0 -: isneg0) *. 0 0 0 0 0 0 -: ispos0) _1 copysign __ _1 _0.0 0 1 _
( _  1 0 0  1  _&-: *. (0 0 0 0 0 0 -: isneg0) *. 0 0 1 1 0 0 -: ispos0)  1 copysign __ _1 _0.0 0 1 _

NB. sorim
'' -: sorim ''
0 2 3 -: sorim 0 2 _3
(_. _ _ 0 2 2 7 7 7 7&-: *. 1 0 0 0 0 0 0 0 0 0 -: isnan) sorim _. __ _ 0 _2 2 _3j_4 _3j4 3j_4 3j4

NB. soris
'' -: soris ''
0 4 9 -: soris 0 2 _3
(_. _ _ 0 4 4 25 25 25 25&-: *. 1 0 0 0 0 0 0 0 0 0 -: isnan) soris _. __ _ 0 _2 2 _3j_4 _3j4 3j_4 3j4

NB. cut3
'' -: cut3 ''
('foo ' ; 'ar  ' ; 'az' , LF , 'qux' , LF2 , 'quux' , LF , ' corge ' , LF , 'flo') -: cut3 string

NB. cut2
(, a:) -: cut2 ''
('foo' ; 'bar' ; '' ; ('baz' , LF , 'qux' , LF2 , 'quux' , LF) ; 'corge' ; LF , 'flob') -: cut2 string
('foo bar  baz' ; 'qux' ; '' ; 'quux' ; ' corge ' ; 'flob') -: LF cut2 string

NB. cut
'' -: cut ''
('foo' ; 'bar' ; ('baz' , LF , 'qux' , LF2 , 'quux' , LF) ; 'corge' ; LF , 'flob') -: cut string
('foo bar  baz' ; 'qux' ; 'quux' ; ' corge ' ; 'flob') -: LF cut string

NB. cutl2
(, a:) -: delimiters cutl2 ''
('foo' ; 'bar' ; '' ; 'baz' ; 'qux' ; '' ; 'quux' ; '' ; 'corge' ; '' ; 'flob') -: delimiters cutl2 string

NB. cutl
'' -: delimiters cutl ''
('foo' ; 'bar' ; 'baz' ; 'qux' ; 'quux' ; 'corge' ; 'flob') -: delimiters cutl string
