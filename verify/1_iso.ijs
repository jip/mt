NB. Verify iso actors
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
NB. Verification suite

NB. liso4th
''    -: 3 liso4th 3
(, 3) -: 4 liso4th 3
3 4   -: 5 liso4th 3

NB. liso4dhs
 4  5  6 -:   liso4dhs  4  3
 4  5  6 -: 1 liso4dhs  4  3
 4  6  8 -: 2 liso4dhs  4  3
_6 _5 _4 -:   liso4dhs _4  3
_6 _5 _4 -: 1 liso4dhs _4  3
_8 _6 _4 -: 2 liso4dhs _4  3
 6  5  4 -:   liso4dhs  4 _3
 6  5  4 -: 1 liso4dhs  4 _3
 8  6  4 -: 2 liso4dhs  4 _3
_4 _5 _6 -:   liso4dhs _4 _3
_4 _5 _6 -: 1 liso4dhs _4 _3
_4 _6 _8 -: 2 liso4dhs _4 _3

NB. iso4riso
riso=. 0 9 3 5 _6 _2 ,: 0 1 2 _3 4 _5
iso=. < '' ; (, 9) ; 3 4 ; 7 6 5 ; _9 _8 _7 _6 ; _2 _3 _4 _5 _6
iso  -: iso4riso     riso
riso -: iso4riso^:_1 iso

NB. liso4riso
sh=. 4 # 10
arr=. i. sh
riso=. 3 5 _6 _2 ,: 2 _3 4 _5
iso=. iso4riso riso
liso=. sh liso4riso riso
(riso ];.0 arr) -: ( iso  {   arr)
(riso ,;.0 arr) -: (liso ({,) arr)

NB. lisoX
arr=. i. 5 6
vec=. _1 _2 _3
vec -: (< 2 ; 3 4 5) { vec ((( 0 lisoE)&c)}) arr
vec -: (< 1 ; 2 3 4) { vec ((( 1 lisoE)&c)}) arr
vec -: (< 1 ; 3 4 5) { vec (((_1 lisoE)&c)}) arr
vec -: (< 2 ; 0 1 2) { vec ((( 0 lisoW)&c)}) arr
vec -: (< 3 ; 0 1 2) { vec ((( 1 lisoW)&c)}) arr
vec -: (< 3 ; 1 2 3) { vec (((_1 lisoW)&c)}) arr
vec -: (< 0 1 2 ; 2) { vec ((( 0 lisoN)&c)}) arr
vec -: (< 0 1 2 ; 3) { vec ((( 1 lisoN)&c)}) arr
vec -: (< 1 2 3 ; 3) { vec (((_1 lisoN)&c)}) arr
vec -: (< 2 3 4 ; 3) { vec ((( 0 lisoS)&c)}) arr
vec -: (< 1 2 3 ; 2) { vec ((( 1 lisoS)&c)}) arr
vec -: (< 2 3 4 ; 2) { vec (((_1 lisoS)&c)}) arr
