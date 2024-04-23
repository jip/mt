NB. Verify quaternion actors
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

'q1 q2 q3 q4'=. 1j2 3j4 , 5j6 7j8 , ,.~ (j.~"0) 1r2 1r8
sqr_big=. 2 *&%: FP_OVFL  NB. its square overflows
sqr_sml=. % sqr_big
q5=. sqr_big * 1j1 1j1
q6=. -: q5
q7=. sqr_sml * 1j1 1j1

NB. =========================================================
NB. Verification suite

NB. qnx

1 -: qn1 q1
2 -: qni q1
3 -: qnj q1
4 -: qnk q1

9j2 3j4 -: 9 qn1 q1
1j9 3j4 -: 9 qni q1
1j2 9j4 -: 9 qnj q1
1j2 3j9 -: 9 qnk q1

NB. qnxx

1j2 -: qn1i q1
1j3 -: qn1j q1
1j4 -: qn1k q1
2j3 -: qnij q1
2j4 -: qnik q1
3j4 -: qnjk q1

0j9 3j4 -: 0j9 qn1i q1
0j2 9j4 -: 0j9 qn1j q1
0j2 3j9 -: 0j9 qn1k q1
1   9j4 -: 0j9 qnij q1
1   3j9 -: 0j9 qnik q1
1j2 0j9 -: 0j9 qnjk q1

NB. qnmarkx
1   0   -: qnmark1 q1
0j2 0   -: qnmarki q1
0   3   -: qnmarkj q1
0   0j4 -: qnmarkk q1

NB. qnmarkxx
1j2 0   -: qnmark1i q1
1   3   -: qnmark1j q1
1   0j4 -: qnmark1k q1
0j2 3   -: qnmarkij q1
0j2 0j4 -: qnmarkik q1
0   3j4 -: qnmarkjk q1

NB. qnmarkxxx
1j2 3   -: qnmark1ij q1
1j2 0j4 -: qnmark1ik q1
1   3j4 -: qnmark1jk q1
0j2 3j4 -: qnmarkijk q1

NB. qnconx
_1j2   3j4  -: qncon1 q1
 1j_2  3j4  -: qnconi q1
 1j2  _3j4  -: qnconj q1
 1j2   3j_4 -: qnconk q1

NB. qnconxx
1j_2 _3j4  -: qnconij q1
1j_2  3j_4 -: qnconik q1
1j2  _3j_4 -: qnconjk q1

NB. qnconv
1j_2 _3j_4 -: qnconv q1

NB. qnlen
NB. - input contains NaN
isnan qnlen  0     0j_.
isnan qnlen  0    _.
isnan qnlen  0j_.  0
isnan qnlen _.     0
isnan qnlen  0    _.j_.
isnan qnlen  0j_.  0j_.
isnan qnlen  0j_. _.
isnan qnlen _.     0j_.
isnan qnlen _.    _.
isnan qnlen _.j_.  0
isnan qnlen _.j_. _.
isnan qnlen _.j_.  0j_.
isnan qnlen _.    _.j_.
isnan qnlen  0j_. _.j_.
isnan qnlen _.j_. _.j_.
1 1 -: isnan qnlen _. 0 ,: 0 _.
1 0 -: isnan qnlen _. 0 ,: q3
0 1 -: isnan qnlen q3   ,: 0 _.
NB. - input contains infinity
_ -: qnlen __j__ __j__
_ -: qnlen __j__ __j_1
_ -: qnlen __j__ __
_ -: qnlen __j__ __j1
_ -: qnlen __j__ __j_
_ -: qnlen __j__ _1j__
_ -: qnlen __j__ _1j_1
_ -: qnlen __j__ _1
_ -: qnlen __j__ _1j1
_ -: qnlen __j__ _1j_
_ -: qnlen __j__  0j__
_ -: qnlen __j__  0j_1
_ -: qnlen __j__  0
_ -: qnlen __j__  0j1
_ -: qnlen __j__  0j_
_ -: qnlen __j__  1j__
_ -: qnlen __j__  1j_1
_ -: qnlen __j__  1
_ -: qnlen __j__  1j1
_ -: qnlen __j__  1j_
_ -: qnlen __j__  _j__
_ -: qnlen __j__  _j_1
_ -: qnlen __j__  _
_ -: qnlen __j__  _j1
_ -: qnlen __j__  _j_
_ -: qnlen __j_1 __j__
_ -: qnlen __j_1 __j_1
_ -: qnlen __j_1 __
_ -: qnlen __j_1 __j1
_ -: qnlen __j_1 __j_
_ -: qnlen __j_1 _1j__
_ -: qnlen __j_1 _1j_1
_ -: qnlen __j_1 _1
_ -: qnlen __j_1 _1j1
_ -: qnlen __j_1 _1j_
_ -: qnlen __j_1  0j__
_ -: qnlen __j_1  0j_1
_ -: qnlen __j_1  0
_ -: qnlen __j_1  0j1
_ -: qnlen __j_1  0j_
_ -: qnlen __j_1  1j__
_ -: qnlen __j_1  1j_1
_ -: qnlen __j_1  1
_ -: qnlen __j_1  1j1
_ -: qnlen __j_1  1j_
_ -: qnlen __j_1  _j__
_ -: qnlen __j_1  _j_1
_ -: qnlen __j_1  _
_ -: qnlen __j_1  _j1
_ -: qnlen __j_1  _j_
_ -: qnlen __    __j__
_ -: qnlen __    __j_1
_ -: qnlen __    __
_ -: qnlen __    __j1
_ -: qnlen __    __j_
_ -: qnlen __    _1j__
_ -: qnlen __    _1j_1
_ -: qnlen __    _1
_ -: qnlen __    _1j1
_ -: qnlen __    _1j_
_ -: qnlen __     0j__
_ -: qnlen __     0j_1
_ -: qnlen __     0
_ -: qnlen __     0j1
_ -: qnlen __     0j_
_ -: qnlen __     1j__
_ -: qnlen __     1j_1
_ -: qnlen __     1
_ -: qnlen __     1j1
_ -: qnlen __     1j_
_ -: qnlen __     _j__
_ -: qnlen __     _j_1
_ -: qnlen __     _
_ -: qnlen __     _j1
_ -: qnlen __     _j_
_ -: qnlen __j1  __j__
_ -: qnlen __j1  __j_1
_ -: qnlen __j1  __
_ -: qnlen __j1  __j1
_ -: qnlen __j1  __j_
_ -: qnlen __j1  _1j__
_ -: qnlen __j1  _1j_1
_ -: qnlen __j1  _1
_ -: qnlen __j1  _1j1
_ -: qnlen __j1  _1j_
_ -: qnlen __j1   0j__
_ -: qnlen __j1   0j_1
_ -: qnlen __j1   0
_ -: qnlen __j1   0j1
_ -: qnlen __j1   0j_
_ -: qnlen __j1   1j__
_ -: qnlen __j1   1j_1
_ -: qnlen __j1   1
_ -: qnlen __j1   1j1
_ -: qnlen __j1   1j_
_ -: qnlen __j1   _j__
_ -: qnlen __j1   _j_1
_ -: qnlen __j1   _
_ -: qnlen __j1   _j1
_ -: qnlen __j1   _j_
_ -: qnlen __j_  __j__
_ -: qnlen __j_  __j_1
_ -: qnlen __j_  __
_ -: qnlen __j_  __j1
_ -: qnlen __j_  __j_
_ -: qnlen __j_  _1j__
_ -: qnlen __j_  _1j_1
_ -: qnlen __j_  _1
_ -: qnlen __j_  _1j1
_ -: qnlen __j_  _1j_
_ -: qnlen __j_   0j__
_ -: qnlen __j_   0j_1
_ -: qnlen __j_   0
_ -: qnlen __j_   0j1
_ -: qnlen __j_   0j_
_ -: qnlen __j_   1j__
_ -: qnlen __j_   1j_1
_ -: qnlen __j_   1
_ -: qnlen __j_   1j1
_ -: qnlen __j_   1j_
_ -: qnlen __j_   _j__
_ -: qnlen __j_   _j_1
_ -: qnlen __j_   _
_ -: qnlen __j_   _j1
_ -: qnlen __j_   _j_
_ -: qnlen _1j__ __j__
_ -: qnlen _1j__ __j_1
_ -: qnlen _1j__ __
_ -: qnlen _1j__ __j1
_ -: qnlen _1j__ __j_
_ -: qnlen _1j__ _1j__
_ -: qnlen _1j__ _1j_1
_ -: qnlen _1j__ _1
_ -: qnlen _1j__ _1j1
_ -: qnlen _1j__ _1j_
_ -: qnlen _1j__  0j__
_ -: qnlen _1j__  0j_1
_ -: qnlen _1j__  0
_ -: qnlen _1j__  0j1
_ -: qnlen _1j__  0j_
_ -: qnlen _1j__  1j__
_ -: qnlen _1j__  1j_1
_ -: qnlen _1j__  1
_ -: qnlen _1j__  1j1
_ -: qnlen _1j__  1j_
_ -: qnlen _1j__  _j__
_ -: qnlen _1j__  _j_1
_ -: qnlen _1j__  _
_ -: qnlen _1j__  _j1
_ -: qnlen _1j__  _j_
_ -: qnlen _1j_1 __j__
_ -: qnlen _1j_1 __j_1
_ -: qnlen _1j_1 __
_ -: qnlen _1j_1 __j1
_ -: qnlen _1j_1 __j_
_ -: qnlen _1j_1 _1j__
_ -: qnlen _1j_1 _1j_
_ -: qnlen _1j_1  0j__
_ -: qnlen _1j_1  0j_
_ -: qnlen _1j_1  1j__
_ -: qnlen _1j_1  1j_
_ -: qnlen _1j_1  _j__
_ -: qnlen _1j_1  _j_1
_ -: qnlen _1j_1  _
_ -: qnlen _1j_1  _j1
_ -: qnlen _1j_1  _j_
_ -: qnlen _1    __j__
_ -: qnlen _1    __j_1
_ -: qnlen _1    __
_ -: qnlen _1    __j1
_ -: qnlen _1    __j_
_ -: qnlen _1    _1j__
_ -: qnlen _1    _1j_
_ -: qnlen _1     0j__
_ -: qnlen _1     0j_
_ -: qnlen _1     1j__
_ -: qnlen _1     1j_
_ -: qnlen _1     _j__
_ -: qnlen _1     _j_1
_ -: qnlen _1     _
_ -: qnlen _1     _j1
_ -: qnlen _1     _j_
_ -: qnlen _1j1  __j__
_ -: qnlen _1j1  __j_1
_ -: qnlen _1j1  __
_ -: qnlen _1j1  __j1
_ -: qnlen _1j1  __j_
_ -: qnlen _1j1  _1j__
_ -: qnlen _1j1  _1j_
_ -: qnlen _1j1   0j__
_ -: qnlen _1j1   0j_
_ -: qnlen _1j1   1j__
_ -: qnlen _1j1   1j_
_ -: qnlen _1j1   _j__
_ -: qnlen _1j1   _j_1
_ -: qnlen _1j1   _
_ -: qnlen _1j1   _j1
_ -: qnlen _1j1   _j_
_ -: qnlen _1j_  __j__
_ -: qnlen _1j_  __j_1
_ -: qnlen _1j_  __
_ -: qnlen _1j_  __j1
_ -: qnlen _1j_  __j_
_ -: qnlen _1j_  _1j__
_ -: qnlen _1j_  _1j_1
_ -: qnlen _1j_  _1
_ -: qnlen _1j_  _1j1
_ -: qnlen _1j_  _1j_
_ -: qnlen _1j_   0j__
_ -: qnlen _1j_   0j_1
_ -: qnlen _1j_   0
_ -: qnlen _1j_   0j1
_ -: qnlen _1j_   0j_
_ -: qnlen _1j_   1j__
_ -: qnlen _1j_   1j_1
_ -: qnlen _1j_   1
_ -: qnlen _1j_   1j1
_ -: qnlen _1j_   1j_
_ -: qnlen _1j_   _j__
_ -: qnlen _1j_   _j_1
_ -: qnlen _1j_   _
_ -: qnlen _1j_   _j1
_ -: qnlen _1j_   _j_
_ -: qnlen  0j__ __j__
_ -: qnlen  0j__ __j_1
_ -: qnlen  0j__ __
_ -: qnlen  0j__ __j1
_ -: qnlen  0j__ __j_
_ -: qnlen  0j__ _1j__
_ -: qnlen  0j__ _1j_1
_ -: qnlen  0j__ _1
_ -: qnlen  0j__ _1j1
_ -: qnlen  0j__ _1j_
_ -: qnlen  0j__  0j__
_ -: qnlen  0j__  0j_1
_ -: qnlen  0j__  0
_ -: qnlen  0j__  0j1
_ -: qnlen  0j__  0j_
_ -: qnlen  0j__  1j__
_ -: qnlen  0j__  1j_1
_ -: qnlen  0j__  1
_ -: qnlen  0j__  1j1
_ -: qnlen  0j__  1j_
_ -: qnlen  0j__  _j__
_ -: qnlen  0j__  _j_1
_ -: qnlen  0j__  _
_ -: qnlen  0j__  _j1
_ -: qnlen  0j__  _j_
_ -: qnlen  0j_1 __j__
_ -: qnlen  0j_1 __j_1
_ -: qnlen  0j_1 __
_ -: qnlen  0j_1 __j1
_ -: qnlen  0j_1 __j_
_ -: qnlen  0j_1 _1j__
_ -: qnlen  0j_1 _1j_
_ -: qnlen  0j_1  0j__
_ -: qnlen  0j_1  0j_
_ -: qnlen  0j_1  1j__
_ -: qnlen  0j_1  1j_
_ -: qnlen  0j_1  _j__
_ -: qnlen  0j_1  _j_1
_ -: qnlen  0j_1  _
_ -: qnlen  0j_1  _j1
_ -: qnlen  0j_1  _j_
_ -: qnlen  0    __j__
_ -: qnlen  0    __j_1
_ -: qnlen  0    __
_ -: qnlen  0    __j1
_ -: qnlen  0    __j_
_ -: qnlen  0    _1j__
_ -: qnlen  0    _1j_
_ -: qnlen  0     0j__
_ -: qnlen  0     0j_
_ -: qnlen  0     1j__
_ -: qnlen  0     1j_
_ -: qnlen  0     _j__
_ -: qnlen  0     _j_1
_ -: qnlen  0     _
_ -: qnlen  0     _j1
_ -: qnlen  0     _j_
_ -: qnlen  0j1  __j__
_ -: qnlen  0j1  __j_1
_ -: qnlen  0j1  __
_ -: qnlen  0j1  __j1
_ -: qnlen  0j1  __j_
_ -: qnlen  0j1  _1j__
_ -: qnlen  0j1  _1j_
_ -: qnlen  0j1   0j__
_ -: qnlen  0j1   0j_
_ -: qnlen  0j1   1j__
_ -: qnlen  0j1   1j_
_ -: qnlen  0j1   _j__
_ -: qnlen  0j1   _j_1
_ -: qnlen  0j1   _
_ -: qnlen  0j1   _j1
_ -: qnlen  0j1   _j_
_ -: qnlen  0j_  __j__
_ -: qnlen  0j_  __j_1
_ -: qnlen  0j_  __
_ -: qnlen  0j_  __j1
_ -: qnlen  0j_  __j_
_ -: qnlen  0j_  _1j__
_ -: qnlen  0j_  _1j_1
_ -: qnlen  0j_  _1
_ -: qnlen  0j_  _1j1
_ -: qnlen  0j_  _1j_
_ -: qnlen  0j_   0j__
_ -: qnlen  0j_   0j_1
_ -: qnlen  0j_   0
_ -: qnlen  0j_   0j1
_ -: qnlen  0j_   0j_
_ -: qnlen  0j_   1j__
_ -: qnlen  0j_   1j_1
_ -: qnlen  0j_   1
_ -: qnlen  0j_   1j1
_ -: qnlen  0j_   1j_
_ -: qnlen  0j_   _j__
_ -: qnlen  0j_   _j_1
_ -: qnlen  0j_   _
_ -: qnlen  0j_   _j1
_ -: qnlen  0j_   _j_
_ -: qnlen  1j__ __j__
_ -: qnlen  1j__ __j_1
_ -: qnlen  1j__ __
_ -: qnlen  1j__ __j1
_ -: qnlen  1j__ __j_
_ -: qnlen  1j__ _1j__
_ -: qnlen  1j__ _1j_1
_ -: qnlen  1j__ _1
_ -: qnlen  1j__ _1j1
_ -: qnlen  1j__ _1j_
_ -: qnlen  1j__  0j__
_ -: qnlen  1j__  0j_1
_ -: qnlen  1j__  0
_ -: qnlen  1j__  0j1
_ -: qnlen  1j__  0j_
_ -: qnlen  1j__  1j__
_ -: qnlen  1j__  1j_1
_ -: qnlen  1j__  1
_ -: qnlen  1j__  1j1
_ -: qnlen  1j__  1j_
_ -: qnlen  1j__  _j__
_ -: qnlen  1j__  _j_1
_ -: qnlen  1j__  _
_ -: qnlen  1j__  _j1
_ -: qnlen  1j__  _j_
_ -: qnlen  1j_1 __j__
_ -: qnlen  1j_1 __j_1
_ -: qnlen  1j_1 __
_ -: qnlen  1j_1 __j1
_ -: qnlen  1j_1 __j_
_ -: qnlen  1j_1 _1j__
_ -: qnlen  1j_1 _1j_
_ -: qnlen  1j_1  0j__
_ -: qnlen  1j_1  0j_
_ -: qnlen  1j_1  1j__
_ -: qnlen  1j_1  1j_
_ -: qnlen  1j_1  _j__
_ -: qnlen  1j_1  _j_1
_ -: qnlen  1j_1  _
_ -: qnlen  1j_1  _j1
_ -: qnlen  1j_1  _j_
_ -: qnlen  1    __j__
_ -: qnlen  1    __j_1
_ -: qnlen  1    __
_ -: qnlen  1    __j1
_ -: qnlen  1    __j_
_ -: qnlen  1    _1j__
_ -: qnlen  1    _1j_
_ -: qnlen  1     0j__
_ -: qnlen  1     0j_
_ -: qnlen  1     1j__
_ -: qnlen  1     1j_
_ -: qnlen  1     _j__
_ -: qnlen  1     _j_1
_ -: qnlen  1     _
_ -: qnlen  1     _j1
_ -: qnlen  1     _j_
_ -: qnlen  1j1  __j__
_ -: qnlen  1j1  __j_1
_ -: qnlen  1j1  __
_ -: qnlen  1j1  __j1
_ -: qnlen  1j1  __j_
_ -: qnlen  1j1  _1j__
_ -: qnlen  1j1  _1j_
_ -: qnlen  1j1   0j__
_ -: qnlen  1j1   0j_
_ -: qnlen  1j1   1j__
_ -: qnlen  1j1   1j_
_ -: qnlen  1j1   _j__
_ -: qnlen  1j1   _j_1
_ -: qnlen  1j1   _
_ -: qnlen  1j1   _j1
_ -: qnlen  1j1   _j_
_ -: qnlen  1j_  __j__
_ -: qnlen  1j_  __j_1
_ -: qnlen  1j_  __
_ -: qnlen  1j_  __j1
_ -: qnlen  1j_  __j_
_ -: qnlen  1j_  _1j__
_ -: qnlen  1j_  _1j_1
_ -: qnlen  1j_  _1
_ -: qnlen  1j_  _1j1
_ -: qnlen  1j_  _1j_
_ -: qnlen  1j_   0j__
_ -: qnlen  1j_   0j_1
_ -: qnlen  1j_   0
_ -: qnlen  1j_   0j1
_ -: qnlen  1j_   0j_
_ -: qnlen  1j_   1j__
_ -: qnlen  1j_   1j_1
_ -: qnlen  1j_   1
_ -: qnlen  1j_   1j1
_ -: qnlen  1j_   1j_
_ -: qnlen  1j_   _j__
_ -: qnlen  1j_   _j_1
_ -: qnlen  1j_   _
_ -: qnlen  1j_   _j1
_ -: qnlen  1j_   _j_
_ -: qnlen  _j__ __j__
_ -: qnlen  _j__ __j_1
_ -: qnlen  _j__ __
_ -: qnlen  _j__ __j1
_ -: qnlen  _j__ __j_
_ -: qnlen  _j__ _1j__
_ -: qnlen  _j__ _1j_1
_ -: qnlen  _j__ _1
_ -: qnlen  _j__ _1j1
_ -: qnlen  _j__ _1j_
_ -: qnlen  _j__  0j__
_ -: qnlen  _j__  0j_1
_ -: qnlen  _j__  0
_ -: qnlen  _j__  0j1
_ -: qnlen  _j__  0j_
_ -: qnlen  _j__  1j__
_ -: qnlen  _j__  1j_1
_ -: qnlen  _j__  1
_ -: qnlen  _j__  1j1
_ -: qnlen  _j__  1j_
_ -: qnlen  _j__  _j__
_ -: qnlen  _j__  _j_1
_ -: qnlen  _j__  _
_ -: qnlen  _j__  _j1
_ -: qnlen  _j__  _j_
_ -: qnlen  _j_1 __j__
_ -: qnlen  _j_1 __j_1
_ -: qnlen  _j_1 __
_ -: qnlen  _j_1 __j1
_ -: qnlen  _j_1 __j_
_ -: qnlen  _j_1 _1j__
_ -: qnlen  _j_1 _1j_1
_ -: qnlen  _j_1 _1
_ -: qnlen  _j_1 _1j1
_ -: qnlen  _j_1 _1j_
_ -: qnlen  _j_1  0j__
_ -: qnlen  _j_1  0j_1
_ -: qnlen  _j_1  0
_ -: qnlen  _j_1  0j1
_ -: qnlen  _j_1  0j_
_ -: qnlen  _j_1  1j__
_ -: qnlen  _j_1  1j_1
_ -: qnlen  _j_1  1
_ -: qnlen  _j_1  1j1
_ -: qnlen  _j_1  1j_
_ -: qnlen  _j_1  _j__
_ -: qnlen  _j_1  _j_1
_ -: qnlen  _j_1  _
_ -: qnlen  _j_1  _j1
_ -: qnlen  _j_1  _j_
_ -: qnlen  _    __j__
_ -: qnlen  _    __j_1
_ -: qnlen  _    __
_ -: qnlen  _    __j1
_ -: qnlen  _    __j_
_ -: qnlen  _    _1j__
_ -: qnlen  _    _1j_1
_ -: qnlen  _    _1
_ -: qnlen  _    _1j1
_ -: qnlen  _    _1j_
_ -: qnlen  _     0j__
_ -: qnlen  _     0j_1
_ -: qnlen  _     0
_ -: qnlen  _     0j1
_ -: qnlen  _     0j_
_ -: qnlen  _     1j__
_ -: qnlen  _     1j_1
_ -: qnlen  _     1
_ -: qnlen  _     1j1
_ -: qnlen  _     1j_
_ -: qnlen  _     _j__
_ -: qnlen  _     _j_1
_ -: qnlen  _     _
_ -: qnlen  _     _j1
_ -: qnlen  _     _j_
_ -: qnlen  _j1  __j__
_ -: qnlen  _j1  __j_1
_ -: qnlen  _j1  __
_ -: qnlen  _j1  __j1
_ -: qnlen  _j1  __j_
_ -: qnlen  _j1  _1j__
_ -: qnlen  _j1  _1j_1
_ -: qnlen  _j1  _1
_ -: qnlen  _j1  _1j1
_ -: qnlen  _j1  _1j_
_ -: qnlen  _j1   0j__
_ -: qnlen  _j1   0j_1
_ -: qnlen  _j1   0
_ -: qnlen  _j1   0j1
_ -: qnlen  _j1   0j_
_ -: qnlen  _j1   1j__
_ -: qnlen  _j1   1j_1
_ -: qnlen  _j1   1
_ -: qnlen  _j1   1j1
_ -: qnlen  _j1   1j_
_ -: qnlen  _j1   _j__
_ -: qnlen  _j1   _j_1
_ -: qnlen  _j1   _
_ -: qnlen  _j1   _j1
_ -: qnlen  _j1   _j_
_ -: qnlen  _j_  __j__
_ -: qnlen  _j_  __j_1
_ -: qnlen  _j_  __
_ -: qnlen  _j_  __j1
_ -: qnlen  _j_  __j_
_ -: qnlen  _j_  _1j__
_ -: qnlen  _j_  _1j_1
_ -: qnlen  _j_  _1
_ -: qnlen  _j_  _1j1
_ -: qnlen  _j_  _1j_
_ -: qnlen  _j_   0j__
_ -: qnlen  _j_   0j_1
_ -: qnlen  _j_   0
_ -: qnlen  _j_   0j1
_ -: qnlen  _j_   0j_
_ -: qnlen  _j_   1j__
_ -: qnlen  _j_   1j_1
_ -: qnlen  _j_   1
_ -: qnlen  _j_   1j1
_ -: qnlen  _j_   1j_
_ -: qnlen  _j_   _j__
_ -: qnlen  _j_   _j_1
_ -: qnlen  _j_   _
_ -: qnlen  _j_   _j1
_ -: qnlen  _j_   _j_
_ _ -: qnlen _   0 ,: 0   __
_ _ -: qnlen _j_ 0 ,: 0j_ __
_ 1 -: qnlen _   0 ,: q3
1 _ -: qnlen q3    ,: 0   __
NB. - edge cases input
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,  -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,  -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen  (-FP_OVFL)              , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen  (-FP_OVFL)              , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen  (-FP_OVFL)              ,  -FP_OVFL
 _               -: qnlen  (-FP_OVFL)              , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen  (-FP_OVFL)              , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,  -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,  -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,  -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,  -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen  (-FP_OVFL)              , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen  (-FP_OVFL)              , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen  (-FP_OVFL)              ,  -FP_UNFL
 FP_OVFL         -: qnlen  (-FP_OVFL)              , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen  (-FP_OVFL)              , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,  -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,  -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   0        j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   0        j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   0        j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   0        j.  FP_OVFL
 _               -: qnlen  (-FP_OVFL)              ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen  (-FP_OVFL)              ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen  (-FP_OVFL)              ,   0        j.  FP_UNFL
 _               -: qnlen  (-FP_OVFL)              ,   0        j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   0        j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   0        j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   0        j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   0        j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen  (-FP_OVFL)              ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen  (-FP_OVFL)              ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen  (-FP_OVFL)              ,   FP_UNFL
 FP_OVFL         -: qnlen  (-FP_OVFL)              ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen  (-FP_OVFL)              ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL
 FP_OVFL         -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen  (-FP_OVFL)              ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen  (-FP_OVFL)              ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen  (-FP_OVFL)              ,   FP_OVFL
 _               -: qnlen  (-FP_OVFL)              ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen  (-FP_OVFL)              ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(sqr_big *    2) -: qnlen ((-sqr_big) j. -sqr_big) , (-sqr_big) j. -sqr_big
(sqr_big *    2) -: qnlen ((-sqr_big) j. -sqr_big) , (-sqr_big) j.  sqr_big
(sqr_big * %: 2) -: qnlen ((-sqr_big) j. -sqr_big) ,   0
(sqr_big *    2) -: qnlen ((-sqr_big) j. -sqr_big) ,   sqr_big  j. -sqr_big
(sqr_big *    2) -: qnlen ((-sqr_big) j. -sqr_big) ,   sqr_big  j.  sqr_big
(sqr_big * %: 2) -: qnlen  (-sqr_big)              ,  -sqr_big
(sqr_big * %: 2) -: qnlen  (-sqr_big)              ,   0        j. -sqr_big
 sqr_big         -: qnlen  (-sqr_big)              ,   0
(sqr_big * %: 2) -: qnlen  (-sqr_big)              ,   0        j.  sqr_big
(sqr_big * %: 2) -: qnlen  (-sqr_big)              ,   sqr_big
(sqr_big *    2) -: qnlen ((-sqr_big) j.  sqr_big) , (-sqr_big) j. -sqr_big
(sqr_big *    2) -: qnlen ((-sqr_big) j.  sqr_big) , (-sqr_big) j.  sqr_big
(sqr_big * %: 2) -: qnlen ((-sqr_big) j.  sqr_big) ,   0
(sqr_big *    2) -: qnlen ((-sqr_big) j.  sqr_big) ,   sqr_big  j. -sqr_big
(sqr_big *    2) -: qnlen ((-sqr_big) j.  sqr_big) ,   sqr_big  j.  sqr_big
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,  -FP_OVFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,  -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen  (-FP_UNFL)              , (-FP_OVFL) j. -FP_OVFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              , (-FP_OVFL) j. -FP_UNFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              ,  -FP_OVFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen  (-FP_UNFL)              , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,  -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,  -FP_OVFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,  -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(FP_UNFL *    2) -: qnlen ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,  -FP_UNFL
(FP_UNFL *    2) -: qnlen ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              , (-FP_UNFL) j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen  (-FP_UNFL)              , (-FP_UNFL) j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen  (-FP_UNFL)              ,  -FP_UNFL
(FP_UNFL * %: 3) -: qnlen  (-FP_UNFL)              , (-FP_UNFL) j.  FP_UNFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              , (-FP_UNFL) j.  FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(FP_UNFL *    2) -: qnlen ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,  -FP_UNFL
(FP_UNFL *    2) -: qnlen ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,  -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   0        j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   0        j.  FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   0        j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   0        j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   0        j.  FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   0        j.  FP_OVFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              ,   0        j. -FP_OVFL
(FP_UNFL * %: 2) -: qnlen  (-FP_UNFL)              ,   0        j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen  (-FP_UNFL)              ,   0        j.  FP_UNFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              ,   0        j.  FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   0        j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   0        j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   0        j.  FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   0        j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   0        j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   0        j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(FP_UNFL *    2) -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL
(FP_UNFL *    2) -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              ,   FP_UNFL  j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen  (-FP_UNFL)              ,   FP_UNFL  j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen  (-FP_UNFL)              ,   FP_UNFL
(FP_UNFL * %: 3) -: qnlen  (-FP_UNFL)              ,   FP_UNFL  j.  FP_UNFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              ,   FP_UNFL  j.  FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(FP_UNFL *    2) -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL
(FP_UNFL *    2) -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen  (-FP_UNFL)              ,   FP_OVFL  j. -FP_OVFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              ,   FP_OVFL  j. -FP_UNFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              ,   FP_OVFL
 FP_OVFL         -: qnlen  (-FP_UNFL)              ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen  (-FP_UNFL)              ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL
 FP_OVFL         -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen (  0        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen (  0        j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen (  0        j. -FP_OVFL) ,  -FP_OVFL
 _               -: qnlen (  0        j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  0        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen (  0        j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen (  0        j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen (  0        j. -FP_OVFL) ,  -FP_UNFL
 FP_OVFL         -: qnlen (  0        j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen (  0        j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen (  0        j. -FP_OVFL) ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen (  0        j. -FP_OVFL) ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen (  0        j. -FP_OVFL) ,   0        j.  FP_UNFL
 _               -: qnlen (  0        j. -FP_OVFL) ,   0        j.  FP_OVFL
 _               -: qnlen (  0        j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen (  0        j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen (  0        j. -FP_OVFL) ,   FP_UNFL
 FP_OVFL         -: qnlen (  0        j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen (  0        j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen (  0        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen (  0        j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen (  0        j. -FP_OVFL) ,   FP_OVFL
 _               -: qnlen (  0        j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  0        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(sqr_big * %: 2) -: qnlen (  0        j. -sqr_big) ,  -sqr_big
(sqr_big * %: 2) -: qnlen (  0        j. -sqr_big) ,   0        j. -sqr_big
 sqr_big         -: qnlen (  0        j. -sqr_big) ,   0
(sqr_big * %: 2) -: qnlen (  0        j. -sqr_big) ,   0        j.  sqr_big
(sqr_big * %: 2) -: qnlen (  0        j. -sqr_big) ,   sqr_big
 _               -: qnlen (  0        j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) ,  -FP_OVFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  0        j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen (  0        j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen (  0        j. -FP_UNFL) ,  -FP_UNFL
(FP_UNFL * %: 3) -: qnlen (  0        j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) ,   0        j. -FP_OVFL
(FP_UNFL * %: 2) -: qnlen (  0        j. -FP_UNFL) ,   0        j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen (  0        j. -FP_UNFL) ,   0        j.  FP_UNFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) ,   0        j.  FP_OVFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen (  0        j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen (  0        j. -FP_UNFL) ,   FP_UNFL
(FP_UNFL * %: 3) -: qnlen (  0        j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen (  0        j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) ,   FP_OVFL
 FP_OVFL         -: qnlen (  0        j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  0        j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(sqr_big * %: 2) -: qnlen    0                     , (-sqr_big) j. -sqr_big
 sqr_big         -: qnlen    0                     ,  -sqr_big
(sqr_big * %: 2) -: qnlen    0                     , (-sqr_big) j.  sqr_big
 sqr_big         -: qnlen    0                     ,   0        j. -sqr_big
 sqr_big         -: qnlen    0                     ,   0        j.  sqr_big
(sqr_big * %: 2) -: qnlen    0                     ,   sqr_big  j. -sqr_big
 sqr_big         -: qnlen    0                     ,   sqr_big
(sqr_big * %: 2) -: qnlen    0                     ,   sqr_big  j.  sqr_big
 _               -: qnlen (  0        j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) ,  -FP_OVFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  0        j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen (  0        j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen (  0        j.  FP_UNFL) ,  -FP_UNFL
(FP_UNFL * %: 3) -: qnlen (  0        j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) ,   0        j. -FP_OVFL
(FP_UNFL * %: 2) -: qnlen (  0        j.  FP_UNFL) ,   0        j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen (  0        j.  FP_UNFL) ,   0        j.  FP_UNFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) ,   0        j.  FP_OVFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen (  0        j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen (  0        j.  FP_UNFL) ,   FP_UNFL
(FP_UNFL * %: 3) -: qnlen (  0        j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen (  0        j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) ,   FP_OVFL
 FP_OVFL         -: qnlen (  0        j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  0        j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(sqr_big * %: 2) -: qnlen (  0        j.  sqr_big) ,  -sqr_big
(sqr_big * %: 2) -: qnlen (  0        j.  sqr_big) ,   0        j. -sqr_big
 sqr_big         -: qnlen (  0        j.  sqr_big) ,   0
(sqr_big * %: 2) -: qnlen (  0        j.  sqr_big) ,   0        j.  sqr_big
(sqr_big * %: 2) -: qnlen (  0        j.  sqr_big) ,   sqr_big
 _               -: qnlen (  0        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen (  0        j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen (  0        j.  FP_OVFL) ,  -FP_OVFL
 _               -: qnlen (  0        j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  0        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen (  0        j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen (  0        j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen (  0        j.  FP_OVFL) ,  -FP_UNFL
 FP_OVFL         -: qnlen (  0        j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen (  0        j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen (  0        j.  FP_OVFL) ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen (  0        j.  FP_OVFL) ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen (  0        j.  FP_OVFL) ,   0        j.  FP_UNFL
 _               -: qnlen (  0        j.  FP_OVFL) ,   0        j.  FP_OVFL
 _               -: qnlen (  0        j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen (  0        j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen (  0        j.  FP_OVFL) ,   FP_UNFL
 FP_OVFL         -: qnlen (  0        j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen (  0        j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen (  0        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen (  0        j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen (  0        j.  FP_OVFL) ,   FP_OVFL
 _               -: qnlen (  0        j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  0        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) ,  -FP_OVFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) ,  -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen    FP_UNFL               , (-FP_OVFL) j. -FP_OVFL
 FP_OVFL         -: qnlen    FP_UNFL               , (-FP_OVFL) j. -FP_UNFL
 FP_OVFL         -: qnlen    FP_UNFL               ,  -FP_OVFL
 FP_OVFL         -: qnlen    FP_UNFL               , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen    FP_UNFL               , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) ,  -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) ,  -FP_OVFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_OVFL) ,  -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(FP_UNFL *    2) -: qnlen (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen (  FP_UNFL  j. -FP_UNFL) ,  -FP_UNFL
(FP_UNFL *    2) -: qnlen (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 FP_OVFL         -: qnlen    FP_UNFL               , (-FP_UNFL) j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen    FP_UNFL               , (-FP_UNFL) j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen    FP_UNFL               ,  -FP_UNFL
(FP_UNFL * %: 3) -: qnlen    FP_UNFL               , (-FP_UNFL) j.  FP_UNFL
 FP_OVFL         -: qnlen    FP_UNFL               , (-FP_UNFL) j.  FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(FP_UNFL *    2) -: qnlen (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen (  FP_UNFL  j.  FP_UNFL) ,  -FP_UNFL
(FP_UNFL *    2) -: qnlen (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_OVFL) ,  -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   0        j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   0        j.  FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   0        j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   0        j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   0        j.  FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   0        j.  FP_OVFL
 FP_OVFL         -: qnlen    FP_UNFL               ,   0        j. -FP_OVFL
(FP_UNFL * %: 2) -: qnlen    FP_UNFL               ,   0        j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen    FP_UNFL               ,   0        j.  FP_UNFL
 FP_OVFL         -: qnlen    FP_UNFL               ,   0        j.  FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   0        j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   0        j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   0        j.  FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   0        j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   0        j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   0        j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(FP_UNFL *    2) -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL
(FP_UNFL *    2) -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 FP_OVFL         -: qnlen    FP_UNFL               ,   FP_UNFL  j. -FP_OVFL
(FP_UNFL * %: 3) -: qnlen    FP_UNFL               ,   FP_UNFL  j. -FP_UNFL
(FP_UNFL * %: 2) -: qnlen    FP_UNFL               ,   FP_UNFL
(FP_UNFL * %: 3) -: qnlen    FP_UNFL               ,   FP_UNFL  j.  FP_UNFL
 FP_OVFL         -: qnlen    FP_UNFL               ,   FP_UNFL  j.  FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(FP_UNFL *    2) -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(FP_UNFL * %: 3) -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL
(FP_UNFL *    2) -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen    FP_UNFL               ,   FP_OVFL  j. -FP_OVFL
 FP_OVFL         -: qnlen    FP_UNFL               ,   FP_OVFL  j. -FP_UNFL
 FP_OVFL         -: qnlen    FP_UNFL               ,   FP_OVFL
 FP_OVFL         -: qnlen    FP_UNFL               ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen    FP_UNFL               ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL
 FP_OVFL         -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,  -FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) ,  -FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen    FP_OVFL               , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen    FP_OVFL               , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen    FP_OVFL               ,  -FP_OVFL
 _               -: qnlen    FP_OVFL               , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen    FP_OVFL               , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) ,  -FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,  -FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,  -FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_OVFL  j. -FP_UNFL) ,  -FP_UNFL
 FP_OVFL         -: qnlen (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen    FP_OVFL               , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen    FP_OVFL               , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen    FP_OVFL               ,  -FP_UNFL
 FP_OVFL         -: qnlen    FP_OVFL               , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen    FP_OVFL               , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_OVFL  j.  FP_UNFL) ,  -FP_UNFL
 FP_OVFL         -: qnlen (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,  -FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   0        j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   0        j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   0        j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   0        j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   0        j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   0        j.  FP_OVFL
 _               -: qnlen    FP_OVFL               ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen    FP_OVFL               ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen    FP_OVFL               ,   0        j.  FP_UNFL
 _               -: qnlen    FP_OVFL               ,   0        j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   0        j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   0        j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   0        j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   0        j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   0        j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   0        j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   0        j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   0        j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL
 FP_OVFL         -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen    FP_OVFL               ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen    FP_OVFL               ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen    FP_OVFL               ,   FP_UNFL
 FP_OVFL         -: qnlen    FP_OVFL               ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen    FP_OVFL               ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
 FP_OVFL         -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
 FP_OVFL         -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL
 FP_OVFL         -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(sqr_big *    2) -: qnlen (  sqr_big  j. -sqr_big) , (-sqr_big) j. -sqr_big
(sqr_big *    2) -: qnlen (  sqr_big  j. -sqr_big) , (-sqr_big) j.  sqr_big
(sqr_big * %: 2) -: qnlen (  sqr_big  j. -sqr_big) ,   0
(sqr_big *    2) -: qnlen (  sqr_big  j. -sqr_big) ,   sqr_big  j. -sqr_big
(sqr_big *    2) -: qnlen (  sqr_big  j. -sqr_big) ,   sqr_big  j.  sqr_big
(sqr_big * %: 2) -: qnlen    sqr_big               ,  -sqr_big
(sqr_big * %: 2) -: qnlen    sqr_big               ,   0        j. -sqr_big
 sqr_big         -: qnlen    sqr_big               ,   0
(sqr_big * %: 2) -: qnlen    sqr_big               ,   0        j.  sqr_big
(sqr_big * %: 2) -: qnlen    sqr_big               ,   sqr_big
(sqr_big *    2) -: qnlen (  sqr_big  j.  sqr_big) , (-sqr_big) j. -sqr_big
(sqr_big *    2) -: qnlen (  sqr_big  j.  sqr_big) , (-sqr_big) j.  sqr_big
(sqr_big * %: 2) -: qnlen (  sqr_big  j.  sqr_big) ,   0
(sqr_big *    2) -: qnlen (  sqr_big  j.  sqr_big) ,   sqr_big  j. -sqr_big
(sqr_big *    2) -: qnlen (  sqr_big  j.  sqr_big) ,   sqr_big  j.  sqr_big
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen    FP_OVFL               ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen    FP_OVFL               ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen    FP_OVFL               ,   FP_OVFL
 _               -: qnlen    FP_OVFL               ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen    FP_OVFL               ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
 _               -: qnlen (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
NB. - input without edge cases
      2  -: qnlen _1j_1 _1j_1
(%:   3) -: qnlen _1j_1 _1
      2  -: qnlen _1j_1 _1j1
(%:   3) -: qnlen _1j_1  0j_1
(%:   2) -: qnlen _1j_1  0
(%:   3) -: qnlen _1j_1  0j1
      2  -: qnlen _1j_1  1j_1
(%:   3) -: qnlen _1j_1  1
      2  -: qnlen _1j_1  1j1
(%:   3) -: qnlen _1    _1j_1
(%:   2) -: qnlen _1    _1
(%:   3) -: qnlen _1    _1j1
(%:   2) -: qnlen _1     0j_1
      1  -: qnlen _1     0
(%:   2) -: qnlen _1     0j1
(%:   3) -: qnlen _1     1j_1
(%:   2) -: qnlen _1     1
(%:   3) -: qnlen _1     1j1
      2  -: qnlen _1j1  _1j_1
(%:   3) -: qnlen _1j1  _1
      2  -: qnlen _1j1  _1j1
(%:   3) -: qnlen _1j1   0j_1
(%:   2) -: qnlen _1j1   0
(%:   3) -: qnlen _1j1   0j1
      2  -: qnlen _1j1   1j_1
(%:   3) -: qnlen _1j1   1
      2  -: qnlen _1j1   1j1
(%:   3) -: qnlen  0j_1 _1j_1
(%:   2) -: qnlen  0j_1 _1
(%:   3) -: qnlen  0j_1 _1j1
(%:   2) -: qnlen  0j_1  0j_1
      1  -: qnlen  0j_1  0
(%:   2) -: qnlen  0j_1  0j1
(%:   3) -: qnlen  0j_1  1j_1
(%:   2) -: qnlen  0j_1  1
(%:   3) -: qnlen  0j_1  1j1
(%:   2) -: qnlen  0    _1j_1
      1  -: qnlen  0    _1
(%:   2) -: qnlen  0    _1j1
      1  -: qnlen  0     0j_1
      0  -: qnlen  0     0
      1  -: qnlen  0     0j1
(%:   2) -: qnlen  0     1j_1
      1  -: qnlen  0     1
(%:   2) -: qnlen  0     1j1
(%:   3) -: qnlen  0j1  _1j_1
(%:   2) -: qnlen  0j1  _1
(%:   3) -: qnlen  0j1  _1j1
(%:   2) -: qnlen  0j1   0j_1
      1  -: qnlen  0j1   0
(%:   2) -: qnlen  0j1   0j1
(%:   3) -: qnlen  0j1   1j_1
(%:   2) -: qnlen  0j1   1
(%:   3) -: qnlen  0j1   1j1
      2  -: qnlen  1j_1 _1j_1
(%:   3) -: qnlen  1j_1 _1
      2  -: qnlen  1j_1 _1j1
(%:   3) -: qnlen  1j_1  0j_1
(%:   2) -: qnlen  1j_1  0
(%:   3) -: qnlen  1j_1  0j1
      2  -: qnlen  1j_1  1j_1
(%:   3) -: qnlen  1j_1  1
      2  -: qnlen  1j_1  1j1
(%:   3) -: qnlen  1    _1j_1
(%:   2) -: qnlen  1    _1
(%:   3) -: qnlen  1    _1j1
(%:   2) -: qnlen  1     0j_1
      1  -: qnlen  1     0
(%:   2) -: qnlen  1     0j1
(%:   3) -: qnlen  1     1j_1
(%:   2) -: qnlen  1     1
(%:   3) -: qnlen  1     1j1
      2  -: qnlen  1j1  _1j_1
(%:   3) -: qnlen  1j1  _1
      2  -: qnlen  1j1  _1j1
(%:   3) -: qnlen  1j1   0j_1
(%:   2) -: qnlen  1j1   0
(%:   3) -: qnlen  1j1   0j1
      2  -: qnlen  1j1   1j_1
(%:   3) -: qnlen  1j1   1
      2  -: qnlen  1j1   1j1
(%:  30) -: qnlen q1
(%: 174) -: qnlen q2
 1       -: qnlen q3
   1r4   -: qnlen       q4
 1 1r4   -: qnlen q3 ,: q4

NB. qnsign
NB. - input contains NaN
1 1 -: isnan qnsign  0     0j_.
1 1 -: isnan qnsign  0    _.
1 1 -: isnan qnsign  0j_.  0
1 1 -: isnan qnsign _.     0
1 1 -: isnan qnsign  0    _.j_.
1 1 -: isnan qnsign  0j_.  0j_.
1 1 -: isnan qnsign  0j_. _.
1 1 -: isnan qnsign _.     0j_.
1 1 -: isnan qnsign _.    _.
1 1 -: isnan qnsign _.j_.  0
1 1 -: isnan qnsign _.j_. _.
1 1 -: isnan qnsign _.j_.  0j_.
1 1 -: isnan qnsign _.    _.j_.
1 1 -: isnan qnsign  0j_. _.j_.
1 1 -: isnan qnsign _.j_. _.j_.
(2 2 $ 1) -: isnan qnsign _. 0 ,: 0 _.
(,.~ 1 0) -: isnan qnsign _. 0 ,: q3
(,.~ 0 1) -: isnan qnsign q3   ,: 0 _.
NB. - result is undefined so NaN error must be throwed
0:@qnsign :: 1: 0 0
NB. - input contains q∞ so NaN error must be throwed
0:@qnsign :: 1: __j__ __j__
0:@qnsign :: 1: __j__ __j_1
0:@qnsign :: 1: __j__ __
0:@qnsign :: 1: __j__ __j1
0:@qnsign :: 1: __j__ __j_
0:@qnsign :: 1: __j__ _1j__
0:@qnsign :: 1: __j__ _1j_1
0:@qnsign :: 1: __j__ _1
0:@qnsign :: 1: __j__ _1j1
0:@qnsign :: 1: __j__ _1j_
0:@qnsign :: 1: __j__  0j__
0:@qnsign :: 1: __j__  0j_1
0:@qnsign :: 1: __j__  0
0:@qnsign :: 1: __j__  0j1
0:@qnsign :: 1: __j__  0j_
0:@qnsign :: 1: __j__  1j__
0:@qnsign :: 1: __j__  1j_1
0:@qnsign :: 1: __j__  1
0:@qnsign :: 1: __j__  1j1
0:@qnsign :: 1: __j__  1j_
0:@qnsign :: 1: __j__  _j__
0:@qnsign :: 1: __j__  _j_1
0:@qnsign :: 1: __j__  _
0:@qnsign :: 1: __j__  _j1
0:@qnsign :: 1: __j__  _j_
0:@qnsign :: 1: __j_1 __j__
0:@qnsign :: 1: __j_1 __j_1
0:@qnsign :: 1: __j_1 __
0:@qnsign :: 1: __j_1 __j1
0:@qnsign :: 1: __j_1 __j_
0:@qnsign :: 1: __j_1 _1j__
0:@qnsign :: 1: __j_1 _1j_
0:@qnsign :: 1: __j_1  0j__
0:@qnsign :: 1: __j_1  0j_
0:@qnsign :: 1: __j_1  1j__
0:@qnsign :: 1: __j_1  1j_
0:@qnsign :: 1: __j_1  _j__
0:@qnsign :: 1: __j_1  _j_1
0:@qnsign :: 1: __j_1  _
0:@qnsign :: 1: __j_1  _j1
0:@qnsign :: 1: __j_1  _j_
0:@qnsign :: 1: __    __j__
0:@qnsign :: 1: __    __j_1
0:@qnsign :: 1: __    __
0:@qnsign :: 1: __    __j1
0:@qnsign :: 1: __    __j_
0:@qnsign :: 1: __    _1j__
0:@qnsign :: 1: __    _1j_
0:@qnsign :: 1: __     0j__
0:@qnsign :: 1: __     0j_
0:@qnsign :: 1: __     1j__
0:@qnsign :: 1: __     1j_
0:@qnsign :: 1: __     _j__
0:@qnsign :: 1: __     _j_1
0:@qnsign :: 1: __     _
0:@qnsign :: 1: __     _j1
0:@qnsign :: 1: __     _j_
0:@qnsign :: 1: __j1  __j__
0:@qnsign :: 1: __j1  __j_1
0:@qnsign :: 1: __j1  __
0:@qnsign :: 1: __j1  __j1
0:@qnsign :: 1: __j1  __j_
0:@qnsign :: 1: __j1  _1j__
0:@qnsign :: 1: __j1  _1j_
0:@qnsign :: 1: __j1   0j__
0:@qnsign :: 1: __j1   0j_
0:@qnsign :: 1: __j1   1j__
0:@qnsign :: 1: __j1   1j_
0:@qnsign :: 1: __j1   _j__
0:@qnsign :: 1: __j1   _j_1
0:@qnsign :: 1: __j1   _
0:@qnsign :: 1: __j1   _j1
0:@qnsign :: 1: __j1   _j_
0:@qnsign :: 1: __j_  __j__
0:@qnsign :: 1: __j_  __j_1
0:@qnsign :: 1: __j_  __
0:@qnsign :: 1: __j_  __j1
0:@qnsign :: 1: __j_  __j_
0:@qnsign :: 1: __j_  _1j__
0:@qnsign :: 1: __j_  _1j_1
0:@qnsign :: 1: __j_  _1
0:@qnsign :: 1: __j_  _1j1
0:@qnsign :: 1: __j_  _1j_
0:@qnsign :: 1: __j_   0j__
0:@qnsign :: 1: __j_   0j_1
0:@qnsign :: 1: __j_   0
0:@qnsign :: 1: __j_   0j1
0:@qnsign :: 1: __j_   0j_
0:@qnsign :: 1: __j_   1j__
0:@qnsign :: 1: __j_   1j_1
0:@qnsign :: 1: __j_   1
0:@qnsign :: 1: __j_   1j1
0:@qnsign :: 1: __j_   1j_
0:@qnsign :: 1: __j_   _j__
0:@qnsign :: 1: __j_   _j_1
0:@qnsign :: 1: __j_   _
0:@qnsign :: 1: __j_   _j1
0:@qnsign :: 1: __j_   _j_
0:@qnsign :: 1: _1j__ __j__
0:@qnsign :: 1: _1j__ __j_1
0:@qnsign :: 1: _1j__ __
0:@qnsign :: 1: _1j__ __j1
0:@qnsign :: 1: _1j__ __j_
0:@qnsign :: 1: _1j__ _1j__
0:@qnsign :: 1: _1j__ _1j_
0:@qnsign :: 1: _1j__  0j__
0:@qnsign :: 1: _1j__  0j_
0:@qnsign :: 1: _1j__  1j__
0:@qnsign :: 1: _1j__  1j_
0:@qnsign :: 1: _1j__  _j__
0:@qnsign :: 1: _1j__  _j_1
0:@qnsign :: 1: _1j__  _
0:@qnsign :: 1: _1j__  _j1
0:@qnsign :: 1: _1j__  _j_
0:@qnsign :: 1: _1j_1 __j__
0:@qnsign :: 1: _1j_1 __j_
0:@qnsign :: 1: _1j_1  _j__
0:@qnsign :: 1: _1j_1  _j_
0:@qnsign :: 1: _1    __j__
0:@qnsign :: 1: _1    __j_
0:@qnsign :: 1: _1     _j__
0:@qnsign :: 1: _1     _j_
0:@qnsign :: 1: _1j1  __j__
0:@qnsign :: 1: _1j1  __j_
0:@qnsign :: 1: _1j1   _j__
0:@qnsign :: 1: _1j1   _j_
0:@qnsign :: 1: _1j_  __j__
0:@qnsign :: 1: _1j_  __j_1
0:@qnsign :: 1: _1j_  __
0:@qnsign :: 1: _1j_  __j1
0:@qnsign :: 1: _1j_  __j_
0:@qnsign :: 1: _1j_  _1j__
0:@qnsign :: 1: _1j_  _1j_
0:@qnsign :: 1: _1j_   0j__
0:@qnsign :: 1: _1j_   0j_
0:@qnsign :: 1: _1j_   1j__
0:@qnsign :: 1: _1j_   1j_
0:@qnsign :: 1: _1j_   _j__
0:@qnsign :: 1: _1j_   _j_1
0:@qnsign :: 1: _1j_   _
0:@qnsign :: 1: _1j_   _j1
0:@qnsign :: 1: _1j_   _j_
0:@qnsign :: 1:  0j__ __j__
0:@qnsign :: 1:  0j__ __j_1
0:@qnsign :: 1:  0j__ __
0:@qnsign :: 1:  0j__ __j1
0:@qnsign :: 1:  0j__ __j_
0:@qnsign :: 1:  0j__ _1j__
0:@qnsign :: 1:  0j__ _1j_
0:@qnsign :: 1:  0j__  0j__
0:@qnsign :: 1:  0j__  0j_
0:@qnsign :: 1:  0j__  1j__
0:@qnsign :: 1:  0j__  1j_
0:@qnsign :: 1:  0j__  _j__
0:@qnsign :: 1:  0j__  _j_1
0:@qnsign :: 1:  0j__  _
0:@qnsign :: 1:  0j__  _j1
0:@qnsign :: 1:  0j__  _j_
0:@qnsign :: 1:  0j_1 __j__
0:@qnsign :: 1:  0j_1 __j_
0:@qnsign :: 1:  0j_1  _j__
0:@qnsign :: 1:  0j_1  _j_
0:@qnsign :: 1:  0    __j__
0:@qnsign :: 1:  0    __j_
0:@qnsign :: 1:  0     _j__
0:@qnsign :: 1:  0     _j_
0:@qnsign :: 1:  0j1  __j__
0:@qnsign :: 1:  0j1  __j_
0:@qnsign :: 1:  0j1   _j__
0:@qnsign :: 1:  0j1   _j_
0:@qnsign :: 1:  0j_  __j__
0:@qnsign :: 1:  0j_  __j_1
0:@qnsign :: 1:  0j_  __
0:@qnsign :: 1:  0j_  __j1
0:@qnsign :: 1:  0j_  __j_
0:@qnsign :: 1:  0j_  _1j__
0:@qnsign :: 1:  0j_  _1j_
0:@qnsign :: 1:  0j_   0j__
0:@qnsign :: 1:  0j_   0j_
0:@qnsign :: 1:  0j_   1j__
0:@qnsign :: 1:  0j_   1j_
0:@qnsign :: 1:  0j_   _j__
0:@qnsign :: 1:  0j_   _j_1
0:@qnsign :: 1:  0j_   _
0:@qnsign :: 1:  0j_   _j1
0:@qnsign :: 1:  0j_   _j_
0:@qnsign :: 1:  1j__ __j__
0:@qnsign :: 1:  1j__ __j_1
0:@qnsign :: 1:  1j__ __
0:@qnsign :: 1:  1j__ __j1
0:@qnsign :: 1:  1j__ __j_
0:@qnsign :: 1:  1j__ _1j__
0:@qnsign :: 1:  1j__ _1j_
0:@qnsign :: 1:  1j__  0j__
0:@qnsign :: 1:  1j__  0j_
0:@qnsign :: 1:  1j__  1j__
0:@qnsign :: 1:  1j__  1j_
0:@qnsign :: 1:  1j__  _j__
0:@qnsign :: 1:  1j__  _j_1
0:@qnsign :: 1:  1j__  _
0:@qnsign :: 1:  1j__  _j1
0:@qnsign :: 1:  1j__  _j_
0:@qnsign :: 1:  1j_1 __j__
0:@qnsign :: 1:  1j_1 __j_
0:@qnsign :: 1:  1j_1  _j__
0:@qnsign :: 1:  1j_1  _j_
0:@qnsign :: 1:  1    __j__
0:@qnsign :: 1:  1    __j_
0:@qnsign :: 1:  1     _j__
0:@qnsign :: 1:  1     _j_
0:@qnsign :: 1:  1j1  __j__
0:@qnsign :: 1:  1j1  __j_
0:@qnsign :: 1:  1j1   _j__
0:@qnsign :: 1:  1j1   _j_
0:@qnsign :: 1:  1j_  __j__
0:@qnsign :: 1:  1j_  __j_1
0:@qnsign :: 1:  1j_  __
0:@qnsign :: 1:  1j_  __j1
0:@qnsign :: 1:  1j_  __j_
0:@qnsign :: 1:  1j_  _1j__
0:@qnsign :: 1:  1j_  _1j_
0:@qnsign :: 1:  1j_   0j__
0:@qnsign :: 1:  1j_   0j_
0:@qnsign :: 1:  1j_   1j__
0:@qnsign :: 1:  1j_   1j_
0:@qnsign :: 1:  1j_   _j__
0:@qnsign :: 1:  1j_   _j_1
0:@qnsign :: 1:  1j_   _
0:@qnsign :: 1:  1j_   _j1
0:@qnsign :: 1:  1j_   _j_
0:@qnsign :: 1:  _j__ __j__
0:@qnsign :: 1:  _j__ __j_1
0:@qnsign :: 1:  _j__ __
0:@qnsign :: 1:  _j__ __j1
0:@qnsign :: 1:  _j__ __j_
0:@qnsign :: 1:  _j__ _1j__
0:@qnsign :: 1:  _j__ _1j_1
0:@qnsign :: 1:  _j__ _1
0:@qnsign :: 1:  _j__ _1j1
0:@qnsign :: 1:  _j__ _1j_
0:@qnsign :: 1:  _j__  0j__
0:@qnsign :: 1:  _j__  0j_1
0:@qnsign :: 1:  _j__  0
0:@qnsign :: 1:  _j__  0j1
0:@qnsign :: 1:  _j__  0j_
0:@qnsign :: 1:  _j__  1j__
0:@qnsign :: 1:  _j__  1j_1
0:@qnsign :: 1:  _j__  1
0:@qnsign :: 1:  _j__  1j1
0:@qnsign :: 1:  _j__  1j_
0:@qnsign :: 1:  _j__  _j__
0:@qnsign :: 1:  _j__  _j_1
0:@qnsign :: 1:  _j__  _
0:@qnsign :: 1:  _j__  _j1
0:@qnsign :: 1:  _j__  _j_
0:@qnsign :: 1:  _j_1 __j__
0:@qnsign :: 1:  _j_1 __j_1
0:@qnsign :: 1:  _j_1 __
0:@qnsign :: 1:  _j_1 __j1
0:@qnsign :: 1:  _j_1 __j_
0:@qnsign :: 1:  _j_1 _1j__
0:@qnsign :: 1:  _j_1 _1j_
0:@qnsign :: 1:  _j_1  0j__
0:@qnsign :: 1:  _j_1  0j_
0:@qnsign :: 1:  _j_1  1j__
0:@qnsign :: 1:  _j_1  1j_
0:@qnsign :: 1:  _j_1  _j__
0:@qnsign :: 1:  _j_1  _j_1
0:@qnsign :: 1:  _j_1  _
0:@qnsign :: 1:  _j_1  _j1
0:@qnsign :: 1:  _j_1  _j_
0:@qnsign :: 1:  _    __j__
0:@qnsign :: 1:  _    __j_1
0:@qnsign :: 1:  _    __
0:@qnsign :: 1:  _    __j1
0:@qnsign :: 1:  _    __j_
0:@qnsign :: 1:  _    _1j__
0:@qnsign :: 1:  _    _1j_
0:@qnsign :: 1:  _     0j__
0:@qnsign :: 1:  _     0j_
0:@qnsign :: 1:  _     1j__
0:@qnsign :: 1:  _     1j_
0:@qnsign :: 1:  _     _j__
0:@qnsign :: 1:  _     _j_1
0:@qnsign :: 1:  _     _
0:@qnsign :: 1:  _     _j1
0:@qnsign :: 1:  _     _j_
0:@qnsign :: 1:  _j1  __j__
0:@qnsign :: 1:  _j1  __j_1
0:@qnsign :: 1:  _j1  __
0:@qnsign :: 1:  _j1  __j1
0:@qnsign :: 1:  _j1  __j_
0:@qnsign :: 1:  _j1  _1j__
0:@qnsign :: 1:  _j1  _1j_
0:@qnsign :: 1:  _j1   0j__
0:@qnsign :: 1:  _j1   0j_
0:@qnsign :: 1:  _j1   1j__
0:@qnsign :: 1:  _j1   1j_
0:@qnsign :: 1:  _j1   _j__
0:@qnsign :: 1:  _j1   _j_1
0:@qnsign :: 1:  _j1   _
0:@qnsign :: 1:  _j1   _j1
0:@qnsign :: 1:  _j1   _j_
0:@qnsign :: 1:  _j_  __j__
0:@qnsign :: 1:  _j_  __j_1
0:@qnsign :: 1:  _j_  __
0:@qnsign :: 1:  _j_  __j1
0:@qnsign :: 1:  _j_  __j_
0:@qnsign :: 1:  _j_  _1j__
0:@qnsign :: 1:  _j_  _1j_1
0:@qnsign :: 1:  _j_  _1
0:@qnsign :: 1:  _j_  _1j1
0:@qnsign :: 1:  _j_  _1j_
0:@qnsign :: 1:  _j_   0j__
0:@qnsign :: 1:  _j_   0j_1
0:@qnsign :: 1:  _j_   0
0:@qnsign :: 1:  _j_   0j1
0:@qnsign :: 1:  _j_   0j_
0:@qnsign :: 1:  _j_   1j__
0:@qnsign :: 1:  _j_   1j_1
0:@qnsign :: 1:  _j_   1
0:@qnsign :: 1:  _j_   1j1
0:@qnsign :: 1:  _j_   1j_
0:@qnsign :: 1:  _j_   _j__
0:@qnsign :: 1:  _j_   _j_1
0:@qnsign :: 1:  _j_   _
0:@qnsign :: 1:  _j_   _j1
0:@qnsign :: 1:  _j_   _j_
NB. - input contains directed infinity and is not trivial
_1     0    -: qnsign ( __        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
_1     0    -: qnsign ( __        j. -FP_OVFL) ,  -FP_OVFL
_1     0    -: qnsign ( __        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
_1     0    -: qnsign ( __        j. -FP_OVFL) ,   0        j. -FP_OVFL
_1     0    -: qnsign ( __        j. -FP_OVFL) ,   0        j.  FP_OVFL
_1     0    -: qnsign ( __        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
_1     0    -: qnsign ( __        j. -FP_OVFL) ,   FP_OVFL
_1     0    -: qnsign ( __        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
_1     0    -: qnsign   __                     , (-FP_OVFL) j. -FP_OVFL
_1     0    -: qnsign   __                     ,  -FP_OVFL
_1     0    -: qnsign   __                     , (-FP_OVFL) j.  FP_OVFL
_1     0    -: qnsign   __                     ,   0        j. -FP_OVFL
_1     0    -: qnsign   __                     ,   0        j.  FP_OVFL
_1     0    -: qnsign   __                     ,   FP_OVFL  j. -FP_OVFL
_1     0    -: qnsign   __                     ,   FP_OVFL
_1     0    -: qnsign   __                     ,   FP_OVFL  j.  FP_OVFL
_1     0    -: qnsign ( __        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
_1     0    -: qnsign ( __        j.  FP_OVFL) ,  -FP_OVFL
_1     0    -: qnsign ( __        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
_1     0    -: qnsign ( __        j.  FP_OVFL) ,   0        j. -FP_OVFL
_1     0    -: qnsign ( __        j.  FP_OVFL) ,   0        j.  FP_OVFL
_1     0    -: qnsign ( __        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
_1     0    -: qnsign ( __        j.  FP_OVFL) ,   FP_OVFL
_1     0    -: qnsign ( __        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
 0j_1  0    -: qnsign ((-FP_OVFL) j. __      ) , (-FP_OVFL) j. -FP_OVFL
 0j_1  0    -: qnsign ((-FP_OVFL) j. __      ) ,  -FP_OVFL
 0j_1  0    -: qnsign ((-FP_OVFL) j. __      ) , (-FP_OVFL) j.  FP_OVFL
 0j_1  0    -: qnsign ((-FP_OVFL) j. __      ) ,   0        j. -FP_OVFL
 0j_1  0    -: qnsign ((-FP_OVFL) j. __      ) ,   0        j.  FP_OVFL
 0j_1  0    -: qnsign ((-FP_OVFL) j. __      ) ,   FP_OVFL  j. -FP_OVFL
 0j_1  0    -: qnsign ((-FP_OVFL) j. __      ) ,   FP_OVFL
 0j_1  0    -: qnsign ((-FP_OVFL) j. __      ) ,   FP_OVFL  j.  FP_OVFL
 0    _1    -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,  __        j. -FP_OVFL
 0    _1    -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,  __
 0    _1    -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,  __        j.  FP_OVFL
 0     0j_1 -: qnsign ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. __
 0     0j1  -: qnsign ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  _
 0     0j_1 -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   0        j. __
 0     0j1  -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  _
 0     0j_1 -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. __
 0     0j1  -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  _
 0     1    -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   _        j. -FP_OVFL
 0     1    -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   _
 0     1    -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   _        j.  FP_OVFL
 0    _1    -: qnsign  (-FP_OVFL)              ,  __        j. -FP_OVFL
 0    _1    -: qnsign  (-FP_OVFL)              ,  __
 0    _1    -: qnsign  (-FP_OVFL)              ,  __        j.  FP_OVFL
 0     0j_1 -: qnsign  (-FP_OVFL)              , (-FP_OVFL) j. __
 0     0j1  -: qnsign  (-FP_OVFL)              , (-FP_OVFL) j.  _
 0     0j_1 -: qnsign  (-FP_OVFL)              ,   0        j. __
 0     0j1  -: qnsign  (-FP_OVFL)              ,   0        j.  _
 0     0j_1 -: qnsign  (-FP_OVFL)              ,   FP_OVFL  j. __
 0     0j1  -: qnsign  (-FP_OVFL)              ,   FP_OVFL  j.  _
 0     1    -: qnsign  (-FP_OVFL)              ,   _        j. -FP_OVFL
 0     1    -: qnsign  (-FP_OVFL)              ,   _
 0     1    -: qnsign  (-FP_OVFL)              ,   _        j.  FP_OVFL
 0    _1    -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,  __        j. -FP_OVFL
 0    _1    -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,  __
 0    _1    -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,  __        j.  FP_OVFL
 0     0j_1 -: qnsign ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. __
 0     0j1  -: qnsign ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  _
 0     0j_1 -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   0        j. __
 0     0j1  -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  _
 0     0j_1 -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. __
 0     0j1  -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  _
 0     1    -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   _        j. -FP_OVFL
 0     1    -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   _
 0     1    -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   _        j.  FP_OVFL
 0j1   0    -: qnsign ((-FP_OVFL) j.  _      ) , (-FP_OVFL) j. -FP_OVFL
 0j1   0    -: qnsign ((-FP_OVFL) j.  _      ) ,  -FP_OVFL
 0j1   0    -: qnsign ((-FP_OVFL) j.  _      ) , (-FP_OVFL) j.  FP_OVFL
 0j1   0    -: qnsign ((-FP_OVFL) j.  _      ) ,   0        j. -FP_OVFL
 0j1   0    -: qnsign ((-FP_OVFL) j.  _      ) ,   0        j.  FP_OVFL
 0j1   0    -: qnsign ((-FP_OVFL) j.  _      ) ,   FP_OVFL  j. -FP_OVFL
 0j1   0    -: qnsign ((-FP_OVFL) j.  _      ) ,   FP_OVFL
 0j1   0    -: qnsign ((-FP_OVFL) j.  _      ) ,   FP_OVFL  j.  FP_OVFL
 0j_1  0    -: qnsign (  0        j. __      ) , (-FP_OVFL) j. -FP_OVFL
 0j_1  0    -: qnsign (  0        j. __      ) ,  -FP_OVFL
 0j_1  0    -: qnsign (  0        j. __      ) , (-FP_OVFL) j.  FP_OVFL
 0j_1  0    -: qnsign (  0        j. __      ) ,   0        j. -FP_OVFL
 0j_1  0    -: qnsign (  0        j. __      ) ,   0        j.  FP_OVFL
 0j_1  0    -: qnsign (  0        j. __      ) ,   FP_OVFL  j. -FP_OVFL
 0j_1  0    -: qnsign (  0        j. __      ) ,   FP_OVFL
 0j_1  0    -: qnsign (  0        j. __      ) ,   FP_OVFL  j.  FP_OVFL
 0    _1    -: qnsign (  0        j. -FP_OVFL) ,  __        j. -FP_OVFL
 0    _1    -: qnsign (  0        j. -FP_OVFL) ,  __
 0    _1    -: qnsign (  0        j. -FP_OVFL) ,  __        j.  FP_OVFL
 0     0j_1 -: qnsign (  0        j. -FP_OVFL) , (-FP_OVFL) j. __
 0     0j1  -: qnsign (  0        j. -FP_OVFL) , (-FP_OVFL) j.  _
 0     0j_1 -: qnsign (  0        j. -FP_OVFL) ,   0        j. __
 0     0j1  -: qnsign (  0        j. -FP_OVFL) ,   0        j.  _
 0     0j_1 -: qnsign (  0        j. -FP_OVFL) ,   FP_OVFL  j. __
 0     0j1  -: qnsign (  0        j. -FP_OVFL) ,   FP_OVFL  j.  _
 0     1    -: qnsign (  0        j. -FP_OVFL) ,   _        j. -FP_OVFL
 0     1    -: qnsign (  0        j. -FP_OVFL) ,   _
 0     1    -: qnsign (  0        j. -FP_OVFL) ,   _        j.  FP_OVFL
 0    _1    -: qnsign (  0        j.  FP_OVFL) ,  __        j. -FP_OVFL
 0    _1    -: qnsign (  0        j.  FP_OVFL) ,  __
 0    _1    -: qnsign (  0        j.  FP_OVFL) ,  __        j.  FP_OVFL
 0     0j_1 -: qnsign (  0        j.  FP_OVFL) , (-FP_OVFL) j. __
 0     0j1  -: qnsign (  0        j.  FP_OVFL) , (-FP_OVFL) j.  _
 0     0j_1 -: qnsign (  0        j.  FP_OVFL) ,   0        j. __
 0     0j1  -: qnsign (  0        j.  FP_OVFL) ,   0        j.  _
 0     0j_1 -: qnsign (  0        j.  FP_OVFL) ,   FP_OVFL  j. __
 0     0j1  -: qnsign (  0        j.  FP_OVFL) ,   FP_OVFL  j.  _
 0     1    -: qnsign (  0        j.  FP_OVFL) ,   _        j. -FP_OVFL
 0     1    -: qnsign (  0        j.  FP_OVFL) ,   _
 0     1    -: qnsign (  0        j.  FP_OVFL) ,   _        j.  FP_OVFL
 0j1   0    -: qnsign (  0        j.  _      ) , (-FP_OVFL) j. -FP_OVFL
 0j1   0    -: qnsign (  0        j.  _      ) ,  -FP_OVFL
 0j1   0    -: qnsign (  0        j.  _      ) , (-FP_OVFL) j.  FP_OVFL
 0j1   0    -: qnsign (  0        j.  _      ) ,   0        j. -FP_OVFL
 0j1   0    -: qnsign (  0        j.  _      ) ,   0        j.  FP_OVFL
 0j1   0    -: qnsign (  0        j.  _      ) ,   FP_OVFL  j. -FP_OVFL
 0j1   0    -: qnsign (  0        j.  _      ) ,   FP_OVFL
 0j1   0    -: qnsign (  0        j.  _      ) ,   FP_OVFL  j.  FP_OVFL
 0j_1  0    -: qnsign (  FP_OVFL  j. __      ) , (-FP_OVFL) j. -FP_OVFL
 0j_1  0    -: qnsign (  FP_OVFL  j. __      ) ,  -FP_OVFL
 0j_1  0    -: qnsign (  FP_OVFL  j. __      ) , (-FP_OVFL) j.  FP_OVFL
 0j_1  0    -: qnsign (  FP_OVFL  j. __      ) ,   0        j. -FP_OVFL
 0j_1  0    -: qnsign (  FP_OVFL  j. __      ) ,   0        j.  FP_OVFL
 0j_1  0    -: qnsign (  FP_OVFL  j. __      ) ,   FP_OVFL  j. -FP_OVFL
 0j_1  0    -: qnsign (  FP_OVFL  j. __      ) ,   FP_OVFL
 0j_1  0    -: qnsign (  FP_OVFL  j. __      ) ,   FP_OVFL  j.  FP_OVFL
 0    _1    -: qnsign (  FP_OVFL  j. -FP_OVFL) ,  __        j. -FP_OVFL
 0    _1    -: qnsign (  FP_OVFL  j. -FP_OVFL) ,  __
 0    _1    -: qnsign (  FP_OVFL  j. -FP_OVFL) ,  __        j.  FP_OVFL
 0     0j_1 -: qnsign (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. __
 0     0j1  -: qnsign (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  _
 0     0j_1 -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   0        j. __
 0     0j1  -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   0        j.  _
 0     0j_1 -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. __
 0     0j1  -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  _
 0     1    -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   _        j. -FP_OVFL
 0     1    -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   _
 0     1    -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   _        j.  FP_OVFL
 0    _1    -: qnsign    FP_OVFL               ,  __        j. -FP_OVFL
 0    _1    -: qnsign    FP_OVFL               ,  __
 0    _1    -: qnsign    FP_OVFL               ,  __        j.  FP_OVFL
 0     0j_1 -: qnsign    FP_OVFL               , (-FP_OVFL) j. __
 0     0j1  -: qnsign    FP_OVFL               , (-FP_OVFL) j.  _
 0     0j_1 -: qnsign    FP_OVFL               ,   0        j. __
 0     0j1  -: qnsign    FP_OVFL               ,   0        j.  _
 0     0j_1 -: qnsign    FP_OVFL               ,   FP_OVFL  j. __
 0     0j1  -: qnsign    FP_OVFL               ,   FP_OVFL  j.  _
 0     1    -: qnsign    FP_OVFL               ,   _        j. -FP_OVFL
 0     1    -: qnsign    FP_OVFL               ,   _
 0     1    -: qnsign    FP_OVFL               ,   _        j.  FP_OVFL
 0    _1    -: qnsign (  FP_OVFL  j.  FP_OVFL) ,  __        j. -FP_OVFL
 0    _1    -: qnsign (  FP_OVFL  j.  FP_OVFL) ,  __
 0    _1    -: qnsign (  FP_OVFL  j.  FP_OVFL) ,  __        j.  FP_OVFL
 0     0j_1 -: qnsign (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. __
 0     0j1  -: qnsign (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  _
 0     0j_1 -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   0        j. __
 0     0j1  -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   0        j.  _
 0     0j_1 -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. __
 0     0j1  -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  _
 0     1    -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   _        j. -FP_OVFL
 0     1    -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   _
 0     1    -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   _        j.  FP_OVFL
 0j1   0    -: qnsign (  FP_OVFL  j.  _      ) , (-FP_OVFL) j. -FP_OVFL
 0j1   0    -: qnsign (  FP_OVFL  j.  _      ) ,  -FP_OVFL
 0j1   0    -: qnsign (  FP_OVFL  j.  _      ) , (-FP_OVFL) j.  FP_OVFL
 0j1   0    -: qnsign (  FP_OVFL  j.  _      ) ,   0        j. -FP_OVFL
 0j1   0    -: qnsign (  FP_OVFL  j.  _      ) ,   0        j.  FP_OVFL
 0j1   0    -: qnsign (  FP_OVFL  j.  _      ) ,   FP_OVFL  j. -FP_OVFL
 0j1   0    -: qnsign (  FP_OVFL  j.  _      ) ,   FP_OVFL
 0j1   0    -: qnsign (  FP_OVFL  j.  _      ) ,   FP_OVFL  j.  FP_OVFL
 1     0    -: qnsign (  _        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 1     0    -: qnsign (  _        j. -FP_OVFL) ,  -FP_OVFL
 1     0    -: qnsign (  _        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 1     0    -: qnsign (  _        j. -FP_OVFL) ,   0        j. -FP_OVFL
 1     0    -: qnsign (  _        j. -FP_OVFL) ,   0        j.  FP_OVFL
 1     0    -: qnsign (  _        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 1     0    -: qnsign (  _        j. -FP_OVFL) ,   FP_OVFL
 1     0    -: qnsign (  _        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
 1     0    -: qnsign    _                     , (-FP_OVFL) j. -FP_OVFL
 1     0    -: qnsign    _                     ,  -FP_OVFL
 1     0    -: qnsign    _                     , (-FP_OVFL) j.  FP_OVFL
 1     0    -: qnsign    _                     ,   0        j. -FP_OVFL
 1     0    -: qnsign    _                     ,   0        j.  FP_OVFL
 1     0    -: qnsign    _                     ,   FP_OVFL  j. -FP_OVFL
 1     0    -: qnsign    _                     ,   FP_OVFL
 1     0    -: qnsign    _                     ,   FP_OVFL  j.  FP_OVFL
 1     0    -: qnsign (  _        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
 1     0    -: qnsign (  _        j.  FP_OVFL) ,  -FP_OVFL
 1     0    -: qnsign (  _        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
 1     0    -: qnsign (  _        j.  FP_OVFL) ,   0        j. -FP_OVFL
 1     0    -: qnsign (  _        j.  FP_OVFL) ,   0        j.  FP_OVFL
 1     0    -: qnsign (  _        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
 1     0    -: qnsign (  _        j.  FP_OVFL) ,   FP_OVFL
 1     0    -: qnsign (  _        j.  FP_OVFL) ,   FP_OVFL j.  FP_OVFL
(1       0       ,: 0       _1     ) -: qnsign (_ , FP_OVFL) ,: FP_OVFL , __
(1r2j1r2 1r2j1r2 ,: 0       _1     ) -: qnsign q3            ,: FP_OVFL , __
(1       0       ,: 1r2j1r2 1r2j1r2) -: qnsign (_ , FP_OVFL) ,: q4
NB. - input is trivial
 _1    0               -: qnsign ( __        j. -FP_OVFL) ,   0
 _1    0               -: qnsign ( __        j. _1      ) ,   0
 _1    0               -: qnsign ( __        j. -FP_UNFL) ,   0
 _1    0               -: qnsign   __                         0
 _1    0               -: qnsign ( __        j.  FP_UNFL) ,   0
 _1    0               -: qnsign ( __        j.  1      ) ,   0
 _1    0               -: qnsign ( __        j.  FP_OVFL) ,   0
  0j_1 0               -: qnsign ((-FP_OVFL) j. __      ) ,   0
(_1j_1 0    % %: 2   ) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   0
( 0,~ _1 j. -%FP_OVFL) -: qnsign ((-FP_OVFL) j. _1      ) ,   0
 _1    0               -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   0
 _1    0               -: qnsign  (-FP_OVFL)              ,   0
 _1    0               -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   0
( 0,~ _1 j.  %FP_OVFL) -: qnsign ((-FP_OVFL) j.  1      ) ,   0
(_1j1  0    % %: 2   ) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   0
  0j1  0               -: qnsign ((-FP_OVFL) j.  _      ) ,   0
  0j_1 0               -: qnsign ((-sqr_big) j. __      ) ,   0
(_1j_1 0    % %: 2   ) -: qnsign ((-sqr_big) j. -sqr_big) ,   0
 _1    0               -: qnsign  (-sqr_big)              ,   0
(_1j1  0    % %: 2   ) -: qnsign ((-sqr_big) j.  sqr_big) ,   0
  0j1  0               -: qnsign ((-sqr_big) j.  _      ) ,   0
  0j_1 0               -: qnsign ( _1        j. __      ) ,   0
( 0,~ _1 j.~-%FP_OVFL) -: qnsign ( _1        j. -FP_OVFL) ,   0
(_1j_1 0    % %: 2   ) -: qnsign ( _1        j. _1      ) ,   0
( 0,~ _1 j. - FP_UNFL) -: qnsign ( _1        j. -FP_UNFL) ,   0
 _1    0               -: qnsign   _1                         0
( 0,~ _1 j.   FP_UNFL) -: qnsign ( _1        j.  FP_UNFL) ,   0
(_1j1  0    % %: 2   ) -: qnsign ( _1        j.  1      ) ,   0
( 0,~  1 j.~-%FP_OVFL) -: qnsign ( _1        j.  FP_OVFL) ,   0
  0j1  0               -: qnsign ( _1        j.  _      ) ,   0
  0j_1 0               -: qnsign ((-FP_UNFL) j. __      ) ,   0
  0j1  0               -: qnsign ((-FP_UNFL) j.  _      ) ,   0
  0j_1 0               -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   0
( 0,~ _1 j.~- FP_UNFL) -: qnsign ((-FP_UNFL) j. _1      ) ,   0
(_1j_1 0    % %: 2   ) -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   0
 _1    0               -: qnsign  (-FP_UNFL)              ,   0
(_1j1  0    % %: 2   ) -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   0
( 0,~  1 j.~- FP_UNFL) -: qnsign ((-FP_UNFL) j.  1      ) ,   0
  0j1  0               -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   0
  0j_1 0               -: qnsign (  0        j. __      ) ,   0
  0j_1 0               -: qnsign (  0        j. -FP_OVFL) ,   0
  0j_1 0               -: qnsign (  0        j. -sqr_big) ,   0
  0j_1 0               -: qnsign (  0        j. _1      ) ,   0
  0j_1 0               -: qnsign (  0        j. -FP_UNFL) ,   0
  0   _1               -: qnsign    0                     ,  __        j. -FP_OVFL
  0   _1               -: qnsign    0                     ,  __        j. _1
  0   _1               -: qnsign    0                     ,  __        j. -FP_UNFL
  0   _1               -: qnsign    0                        __
  0   _1               -: qnsign    0                     ,  __        j.  FP_UNFL
  0   _1               -: qnsign    0                     ,  __        j.  1
  0   _1               -: qnsign    0                     ,  __        j.  FP_OVFL
  0    0j_1            -: qnsign    0                     , (-FP_OVFL) j. __
( 0   _1j_1 % %:2    ) -: qnsign    0                     , (-FP_OVFL) j. -FP_OVFL
( 0,  _1 j. -%FP_OVFL) -: qnsign    0                     , (-FP_OVFL) j. _1
  0   _1               -: qnsign    0                     , (-FP_OVFL) j. -FP_UNFL
  0   _1               -: qnsign    0                     ,  -FP_OVFL
  0   _1               -: qnsign    0                     , (-FP_OVFL) j.  FP_UNFL
( 0,  _1 j.  %FP_OVFL) -: qnsign    0                     , (-FP_OVFL) j.  1
( 0   _1j1  % %:2    ) -: qnsign    0                     , (-FP_OVFL) j.  FP_OVFL
  0    0j1             -: qnsign    0                     , (-FP_OVFL) j.  _
  0    0j_1            -: qnsign    0                     , (-sqr_big) j. __
( 0   _1j_1 % %:2    ) -: qnsign    0                     , (-sqr_big) j. -sqr_big
  0   _1               -: qnsign    0                     ,  -sqr_big
( 0   _1j1  % %:2    ) -: qnsign    0                     , (-sqr_big) j.  sqr_big
  0    0j1             -: qnsign    0                     , (-sqr_big) j.  _
  0    0j_1            -: qnsign    0                     ,  _1        j. __
( 0,  _1 j.~-%FP_OVFL) -: qnsign    0                     ,  _1        j. -FP_OVFL
( 0   _1j_1 % %:2    ) -: qnsign    0                     ,  _1        j. _1
( 0,  _1 j. - FP_UNFL) -: qnsign    0                     ,  _1        j. -FP_UNFL
  0   _1               -: qnsign    0                        _1
( 0,  _1 j.   FP_UNFL) -: qnsign    0                     ,  _1        j.  FP_UNFL
( 0   _1j1  % %:2    ) -: qnsign    0                     ,  _1        j.  1
( 0,   1 j.~-%FP_OVFL) -: qnsign    0                     ,  _1        j.  FP_OVFL
  0    0j1             -: qnsign    0                     ,  _1        j.  _
  0    0j_1            -: qnsign    0                     , (-FP_UNFL) j. __
  0    0j_1            -: qnsign    0                     , (-FP_UNFL) j. -FP_OVFL
( 0,  _1 j.~- FP_UNFL) -: qnsign    0                     , (-FP_UNFL) j. _1
( 0   _1j_1 % %:2    ) -: qnsign    0                     , (-FP_UNFL) j. -FP_UNFL
  0   _1               -: qnsign    0                     ,  -FP_UNFL
( 0   _1j1  % %:2    ) -: qnsign    0                     , (-FP_UNFL) j.  FP_UNFL
( 0,   1 j.~- FP_UNFL) -: qnsign    0                     , (-FP_UNFL) j.  1
  0    0j1             -: qnsign    0                     , (-FP_UNFL) j.  FP_OVFL
  0    0j1             -: qnsign    0                     , (-FP_UNFL) j.  _
  0    0j_1            -: qnsign    0                     ,   0        j. __
  0    0j_1            -: qnsign    0                     ,   0        j. -FP_OVFL
  0    0j_1            -: qnsign    0                     ,   0        j. -sqr_big
  0    0j_1            -: qnsign    0                     ,   0        j. _1
  0    0j_1            -: qnsign    0                     ,   0        j. -FP_UNFL
  0    0j1             -: qnsign    0                     ,   0        j.  FP_UNFL
  0    0j1             -: qnsign    0                     ,   0        j.  1
  0    0j1             -: qnsign    0                     ,   0        j.  sqr_big
  0    0j1             -: qnsign    0                     ,   0        j.  FP_OVFL
  0    0j1             -: qnsign    0                     ,   0        j.  _
  0    0j_1            -: qnsign    0                     ,   FP_UNFL  j. __
  0    0j_1            -: qnsign    0                     ,   FP_UNFL  j. -FP_OVFL
( 0,  _1 j.~  FP_UNFL) -: qnsign    0                     ,   FP_UNFL  j. _1
( 0    1j_1 % %: 2   ) -: qnsign    0                     ,   FP_UNFL  j. -FP_UNFL
  0    1               -: qnsign    0                     ,   FP_UNFL
( 0    1j1  % %: 2   ) -: qnsign    0                     ,   FP_UNFL  j.  FP_UNFL
( 0,   1 j.~  FP_UNFL) -: qnsign    0                     ,   FP_UNFL  j.  1
  0    0j1             -: qnsign    0                     ,   FP_UNFL  j.  FP_OVFL
  0    0j1             -: qnsign    0                     ,   FP_UNFL  j.  _
  0    0j_1            -: qnsign    0                     ,   1        j. __
( 0,  _1 j.~ %FP_OVFL) -: qnsign    0                     ,   1        j. -FP_OVFL
( 0    1j_1 % %: 2   ) -: qnsign    0                     ,   1        j. _1
( 0,   1 j. - FP_UNFL) -: qnsign    0                     ,   1        j. -FP_UNFL
  0    1               -: qnsign    0                         1
( 0,   1 j.   FP_UNFL) -: qnsign    0                     ,   1        j.  FP_UNFL
( 0    1j1  % %: 2   ) -: qnsign    0                     ,   1        j.  1
( 0,   1 j.~ %FP_OVFL) -: qnsign    0                     ,   1        j.  FP_OVFL
  0    0j1             -: qnsign    0                     ,   1        j.  _
  0    0j_1            -: qnsign    0                     ,   sqr_big  j. __
  0    1               -: qnsign    0                     ,   sqr_big
  0    0j1             -: qnsign    0                     ,   sqr_big  j.  _
  0    0j_1            -: qnsign    0                     ,   FP_OVFL  j. __
( 0    1j_1 % %: 2   ) -: qnsign    0                     ,   FP_OVFL  j. -FP_OVFL
( 0,   1 j. -%FP_OVFL) -: qnsign    0                     ,   FP_OVFL  j. _1
  0    1               -: qnsign    0                     ,   FP_OVFL  j. -FP_UNFL
  0    1               -: qnsign    0                     ,   FP_OVFL
  0    1               -: qnsign    0                     ,   FP_OVFL  j.  FP_UNFL
( 0,   1 j.  %FP_OVFL) -: qnsign    0                     ,   FP_OVFL  j.  1
( 0    1j1  % %: 2   ) -: qnsign    0                     ,   FP_OVFL  j.  FP_OVFL
  0    0j1             -: qnsign    0                     ,   FP_OVFL  j.  _
  0    1               -: qnsign    0                     ,   _        j. -FP_OVFL
  0    1               -: qnsign    0                     ,   _        j. _1
  0    1               -: qnsign    0                     ,   _        j. -FP_UNFL
  0    1               -: qnsign    0                         _
  0    1               -: qnsign    0                     ,   _        j.  FP_UNFL
  0    1               -: qnsign    0                     ,   _        j.  1
  0    1               -: qnsign    0                     ,   _        j.  FP_OVFL
  0j1  0               -: qnsign (  0        j.  FP_UNFL) ,   0
  0j1  0               -: qnsign (  0        j.  1      ) ,   0
  0j1  0               -: qnsign (  0        j.  FP_OVFL) ,   0
  0j1  0               -: qnsign (  0        j.  _      ) ,   0
  0j_1 0               -: qnsign (  FP_UNFL  j. __      ) ,   0
  0j_1 0               -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   0
( 0,~ _1 j.~  FP_UNFL) -: qnsign (  FP_UNFL  j. _1      ) ,   0
( 1j_1 0    % %: 2   ) -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   0
  1    0               -: qnsign    FP_UNFL               ,   0
( 1j1  0    % %: 2   ) -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   0
( 0,~  1 j.~  FP_UNFL) -: qnsign (  FP_UNFL  j.  1      ) ,   0
  0j1  0               -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   0
  0j1  0               -: qnsign (  FP_UNFL  j.  _      ) ,   0
  0j_1 0               -: qnsign (  1        j. __      ) ,   0
( 0,~ _1 j.~ %FP_OVFL) -: qnsign (  1        j. -FP_OVFL) ,   0
( 1j_1 0    % %: 2   ) -: qnsign (  1        j. _1      ) ,   0
( 0,~  1 j. - FP_UNFL) -: qnsign (  1        j. -FP_UNFL) ,   0
  1    0               -: qnsign    1                         0
( 0,~  1 j.   FP_UNFL) -: qnsign (  1        j.  FP_UNFL) ,   0
( 1j1  0    % %: 2   ) -: qnsign (  1        j.  1      ) ,   0
( 0,~  1 j.~ %FP_OVFL) -: qnsign (  1        j.  FP_OVFL) ,   0
  0j1  0               -: qnsign (  1        j.  _      ) ,   0
  0j_1 0               -: qnsign (  sqr_big  j. __      ) ,   0
( 1j_1 0    % %: 2   ) -: qnsign (  sqr_big  j. -sqr_big) ,   0
  1    0               -: qnsign    sqr_big               ,   0
( 1j1  0    % %: 2   ) -: qnsign (  sqr_big  j.  sqr_big) ,   0
  0j1  0               -: qnsign (  sqr_big  j.  _      ) ,   0
  0j_1 0               -: qnsign (  FP_OVFL  j. __      ) ,   0
( 1j_1 0    % %: 2   ) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   0
( 0,~  1 j. -%FP_OVFL) -: qnsign (  FP_OVFL  j. _1      ) ,   0
  1    0               -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   0
  1    0               -: qnsign    FP_OVFL               ,   0
  1    0               -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   0
( 0,~  1 j.  %FP_OVFL) -: qnsign (  FP_OVFL  j.  1      ) ,   0
( 1j1  0    % %: 2   ) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   0
  0j1  0               -: qnsign (  FP_OVFL  j.  _      ) ,   0
  1    0               -: qnsign (  _        j. _1      ) ,   0
  1    0               -: qnsign (  _        j. -FP_OVFL) ,   0
  1    0               -: qnsign (  _        j. -FP_UNFL) ,   0
  1    0               -: qnsign    _                         0
  1    0               -: qnsign (  _        j.  FP_UNFL) ,   0
  1    0               -: qnsign (  _        j.  1      ) ,   0
  1    0               -: qnsign (  _        j.  FP_OVFL) ,   0
NB. - edge cases input
( 1j1   1j1  %   _2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(_1j_1 _1    % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(_1j_1 _1    % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,  -FP_OVFL
(_1j_1 _1    % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(_1j_1 _1j1  %    2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(_1    _1j_1 % %: 3) -: qnsign ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(_1    _1    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(_1    _1    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,  -FP_OVFL
(_1    _1    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(_1    _1j1  % %: 3) -: qnsign ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(_1    _1j_1 % %: 3) -: qnsign  (-FP_OVFL)              , (-FP_OVFL) j. -FP_OVFL
(_1    _1    % %: 2) -: qnsign  (-FP_OVFL)              , (-FP_OVFL) j. -FP_UNFL
(_1    _1    % %: 2) -: qnsign  (-FP_OVFL)              ,  -FP_OVFL
(_1    _1    % %: 2) -: qnsign  (-FP_OVFL)              , (-FP_OVFL) j.  FP_UNFL
(_1    _1j1  % %: 3) -: qnsign  (-FP_OVFL)              , (-FP_OVFL) j.  FP_OVFL
(_1    _1j_1 % %: 3) -: qnsign ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(_1    _1    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(_1    _1    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,  -FP_OVFL
(_1    _1    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(_1    _1j1  % %: 3) -: qnsign ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
( 1j_1  1j1  %   _2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(_1j1  _1    % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(_1j1  _1    % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,  -FP_OVFL
(_1j1  _1    % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(_1j1  _1j1  %    2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(_1j_1  0j_1 % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(_1j_1  0    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(_1j_1  0    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,  -FP_UNFL
(_1j_1  0    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(_1j_1  0j1  % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(_1     0j_1 % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
 _1     0            -: qnsign ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
 _1     0            -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,  -FP_UNFL
 _1     0            -: qnsign ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(_1     0j1  % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(_1     0j_1 % %: 2) -: qnsign  (-FP_OVFL)              , (-FP_UNFL) j. -FP_OVFL
 _1     0            -: qnsign  (-FP_OVFL)              , (-FP_UNFL) j. -FP_UNFL
 _1     0            -: qnsign  (-FP_OVFL)              ,  -FP_UNFL
 _1     0            -: qnsign  (-FP_OVFL)              , (-FP_UNFL) j.  FP_UNFL
(_1     0j1  % %: 2) -: qnsign  (-FP_OVFL)              , (-FP_UNFL) j.  FP_OVFL
(_1     0j_1 % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
 _1     0            -: qnsign ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
 _1     0            -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,  -FP_UNFL
 _1     0            -: qnsign ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(_1     0j1  % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(_1j1   0j_1 % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(_1j1   0    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(_1j1   0    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,  -FP_UNFL
(_1j1   0    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(_1j1   0j1  % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(_1j_1  0j_1 % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   0        j. -FP_OVFL
(_1j_1  0    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   0        j. -FP_UNFL
(_1j_1  0    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  FP_UNFL
(_1j_1  0j1  % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  FP_OVFL
(_1     0j_1 % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   0        j. -FP_OVFL
 _1     0            -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   0        j. -FP_UNFL
 _1     0            -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   0        j.  FP_UNFL
(_1     0j1  % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   0        j.  FP_OVFL
(_1     0j_1 % %: 2) -: qnsign  (-FP_OVFL)              ,   0        j. -FP_OVFL
 _1     0            -: qnsign  (-FP_OVFL)              ,   0        j. -FP_UNFL
 _1     0            -: qnsign  (-FP_OVFL)              ,   0        j.  FP_UNFL
(_1     0j1  % %: 2) -: qnsign  (-FP_OVFL)              ,   0        j.  FP_OVFL
(_1     0j_1 % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   0        j. -FP_OVFL
 _1     0            -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   0        j. -FP_UNFL
 _1     0            -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   0        j.  FP_UNFL
(_1     0j1  % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   0        j.  FP_OVFL
(_1j1   0j_1 % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   0        j. -FP_OVFL
(_1j1   0    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   0        j. -FP_UNFL
(_1j1   0    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  FP_UNFL
(_1j1   0j1  % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  FP_OVFL
(_1j_1  0j_1 % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(_1j_1  0    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(_1j_1  0    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL
(_1j_1  0    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(_1j_1  0j1  % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(_1     0j_1 % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
 _1     0            -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
 _1     0            -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL
 _1     0            -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(_1     0j1  % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(_1     0j_1 % %: 2) -: qnsign  (-FP_OVFL)              ,   FP_UNFL  j. -FP_OVFL
 _1     0            -: qnsign  (-FP_OVFL)              ,   FP_UNFL  j. -FP_UNFL
 _1     0            -: qnsign  (-FP_OVFL)              ,   FP_UNFL
 _1     0            -: qnsign  (-FP_OVFL)              ,   FP_UNFL  j.  FP_UNFL
(_1     0j1  % %: 2) -: qnsign  (-FP_OVFL)              ,   FP_UNFL  j.  FP_OVFL
(_1     0j_1 % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
 _1     0            -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
 _1     0            -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL
 _1     0            -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(_1     0j1  % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(_1j1   0j_1 % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(_1j1   0    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(_1j1   0    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL
(_1j1   0    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(_1j1   0j1  % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
( 1j1  _1j1  %   _2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(_1j_1  1    % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(_1j_1  1    % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL
(_1j_1  1    % %: 3) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(_1j_1  1j1  %    2) -: qnsign ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(_1     1j_1 % %: 3) -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(_1     1    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(_1     1    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL
(_1     1    % %: 2) -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(_1     1j1  % %: 3) -: qnsign ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(_1     1j_1 % %: 3) -: qnsign  (-FP_OVFL)              ,   FP_OVFL  j. -FP_OVFL
(_1     1    % %: 2) -: qnsign  (-FP_OVFL)              ,   FP_OVFL  j. -FP_UNFL
(_1     1    % %: 2) -: qnsign  (-FP_OVFL)              ,   FP_OVFL
(_1     1    % %: 2) -: qnsign  (-FP_OVFL)              ,   FP_OVFL  j.  FP_UNFL
(_1     1j1  % %: 3) -: qnsign  (-FP_OVFL)              ,   FP_OVFL  j.  FP_OVFL
(_1     1j_1 % %: 3) -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(_1     1    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(_1     1    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL
(_1     1    % %: 2) -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(_1     1j1  % %: 3) -: qnsign ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(_1j1   1j_1 %    2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(_1j1   1    % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(_1j1   1    % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL
(_1j1   1    % %: 3) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(_1j1   1j1  %    2) -: qnsign ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
( 1j1   1j1  %   _2) -: qnsign ((-sqr_big) j. -sqr_big) , (-sqr_big) j. -sqr_big
( 1j1   1j_1 %   _2) -: qnsign ((-sqr_big) j. -sqr_big) , (-sqr_big) j.  sqr_big
( 1j1  _1j1  %   _2) -: qnsign ((-sqr_big) j. -sqr_big) ,   sqr_big  j. -sqr_big
( 1j1  _1j_1 %   _2) -: qnsign ((-sqr_big) j. -sqr_big) ,   sqr_big  j.  sqr_big
(_1    _1    % %: 2) -: qnsign  (-sqr_big)              ,  -sqr_big
(_1     0j_1 % %: 2) -: qnsign  (-sqr_big)              ,   0        j. -sqr_big
(_1     0j1  % %: 2) -: qnsign  (-sqr_big)              ,   0        j.  sqr_big
(_1     1    % %: 2) -: qnsign  (-sqr_big)              ,   sqr_big
( 1j_1  1j1  %   _2) -: qnsign ((-sqr_big) j.  sqr_big) , (-sqr_big) j. -sqr_big
( 1j_1  1j_1 %   _2) -: qnsign ((-sqr_big) j.  sqr_big) , (-sqr_big) j.  sqr_big
( 1j_1 _1j1  %   _2) -: qnsign ((-sqr_big) j.  sqr_big) ,   sqr_big  j. -sqr_big
( 1j_1 _1j_1 %   _2) -: qnsign ((-sqr_big) j.  sqr_big) ,   sqr_big  j.  sqr_big
( 0j_1 _1j_1 % %: 3) -: qnsign ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
( 0j_1 _1    % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
( 0j_1 _1    % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,  -FP_OVFL
( 0j_1 _1    % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
( 0j_1 _1j1  % %: 3) -: qnsign ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
( 0    _1j_1 % %: 2) -: qnsign ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  0    _1            -: qnsign ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  0    _1            -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,  -FP_OVFL
  0    _1            -: qnsign ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
( 0    _1j1  % %: 2) -: qnsign ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
( 0    _1j_1 % %: 2) -: qnsign  (-FP_UNFL)              , (-FP_OVFL) j. -FP_OVFL
  0    _1            -: qnsign  (-FP_UNFL)              , (-FP_OVFL) j. -FP_UNFL
  0    _1            -: qnsign  (-FP_UNFL)              ,  -FP_OVFL
  0    _1            -: qnsign  (-FP_UNFL)              , (-FP_OVFL) j.  FP_UNFL
( 0    _1j1  % %: 2) -: qnsign  (-FP_UNFL)              , (-FP_OVFL) j.  FP_OVFL
( 0    _1j_1 % %: 2) -: qnsign ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  0    _1            -: qnsign ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  0    _1            -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,  -FP_OVFL
  0    _1            -: qnsign ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
( 0    _1j1  % %: 2) -: qnsign ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
( 0j1  _1j_1 % %: 3) -: qnsign ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
( 0j1  _1    % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
( 0j1  _1    % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,  -FP_OVFL
( 0j1  _1    % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
( 0j1  _1j1  % %: 3) -: qnsign ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
( 0j_1  0j_1 % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  0j_1  0            -: qnsign ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  0j_1  0            -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,  -FP_UNFL
  0j_1  0            -: qnsign ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
( 0j_1  0j1  % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  0     0j_1         -: qnsign ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
( 1j1   1j1  %   _2) -: qnsign ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(_1j_1 _1    % %: 3) -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,  -FP_UNFL
(_1j_1 _1j1  %    2) -: qnsign ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  0     0j1          -: qnsign ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  0     0j_1         -: qnsign  (-FP_UNFL)              , (-FP_UNFL) j. -FP_OVFL
(_1    _1j_1 % %: 3) -: qnsign  (-FP_UNFL)              , (-FP_UNFL) j. -FP_UNFL
(_1    _1    % %: 2) -: qnsign  (-FP_UNFL)              ,  -FP_UNFL
(_1    _1j1  % %: 3) -: qnsign  (-FP_UNFL)              , (-FP_UNFL) j.  FP_UNFL
  0     0j1          -: qnsign  (-FP_UNFL)              , (-FP_UNFL) j.  FP_OVFL
  0     0j_1         -: qnsign ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
( 1j_1  1j1  %   _2) -: qnsign ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(_1j1  _1    % %: 3) -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,  -FP_UNFL
(_1j1  _1j1  %    2) -: qnsign ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  0     0j1          -: qnsign ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
( 0j1   0j_1 % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  0j1   0            -: qnsign ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  0j1   0            -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,  -FP_UNFL
  0j1   0            -: qnsign ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
( 0j1   0j1  % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
( 0j_1  0j_1 % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   0        j. -FP_OVFL
  0j_1  0            -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   0        j. -FP_UNFL
  0j_1  0            -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   0        j.  FP_UNFL
( 0j_1  0j1  % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   0        j.  FP_OVFL
  0     0j_1         -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   0        j. -FP_OVFL
(_1j_1  0j_1 % %: 3) -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   0        j. -FP_UNFL
(_1j_1  0j1  % %: 3) -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   0        j.  FP_UNFL
  0     0j1          -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   0        j.  FP_OVFL
  0     0j_1         -: qnsign  (-FP_UNFL)              ,   0        j. -FP_OVFL
(_1     0j_1 % %: 2) -: qnsign  (-FP_UNFL)              ,   0        j. -FP_UNFL
(_1     0j1  % %: 2) -: qnsign  (-FP_UNFL)              ,   0        j.  FP_UNFL
  0     0j1          -: qnsign  (-FP_UNFL)              ,   0        j.  FP_OVFL
  0     0j_1         -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   0        j. -FP_OVFL
(_1j1   0j_1 % %: 3) -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   0        j. -FP_UNFL
(_1j1   0j1  % %: 3) -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   0        j.  FP_UNFL
  0     0j1          -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   0        j.  FP_OVFL
( 0j1   0j_1 % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   0        j. -FP_OVFL
  0j1   0            -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   0        j. -FP_UNFL
  0j1   0            -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   0        j.  FP_UNFL
( 0j1   0j1  % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   0        j.  FP_OVFL
( 0j_1  0j_1 % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  0j_1  0            -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  0j_1  0            -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL
  0j_1  0            -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
( 0j_1  0j1  % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  0     0j_1         -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
( 1j1  _1j1  %   _2) -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(_1j_1  1    % %: 3) -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL
(_1j_1  1j1  %    2) -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  0     0j1          -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  0     0j_1         -: qnsign  (-FP_UNFL)              ,   FP_UNFL  j. -FP_OVFL
(_1     1j_1 % %: 3) -: qnsign  (-FP_UNFL)              ,   FP_UNFL  j. -FP_UNFL
(_1     1    % %: 2) -: qnsign  (-FP_UNFL)              ,   FP_UNFL
(_1     1j1  % %: 3) -: qnsign  (-FP_UNFL)              ,   FP_UNFL  j.  FP_UNFL
  0     0j1          -: qnsign  (-FP_UNFL)              ,   FP_UNFL  j.  FP_OVFL
  0     0j_1         -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(_1j1   1j_1 %    2) -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(_1j1   1    % %: 3) -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL
(_1j1   1j1  %    2) -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  0     0j1          -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
( 0j1   0j_1 % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  0j1   0            -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  0j1   0            -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL
  0j1   0            -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
( 0j1   0j1  % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
( 0j_1  1j_1 % %: 3) -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
( 0j_1  1    % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
( 0j_1  1    % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL
( 0j_1  1    % %: 2) -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
( 0j_1  1j1  % %: 3) -: qnsign ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
( 0     1j_1 % %: 2) -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  0     1            -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  0     1            -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL
  0     1            -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
( 0     1j1  % %: 2) -: qnsign ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
( 0     1j_1 % %: 2) -: qnsign  (-FP_UNFL)              ,   FP_OVFL  j. -FP_OVFL
  0     1            -: qnsign  (-FP_UNFL)              ,   FP_OVFL  j. -FP_UNFL
  0     1            -: qnsign  (-FP_UNFL)              ,   FP_OVFL
  0     1            -: qnsign  (-FP_UNFL)              ,   FP_OVFL  j.  FP_UNFL
( 0     1j1  % %: 2) -: qnsign  (-FP_UNFL)              ,   FP_OVFL  j.  FP_OVFL
( 0     1j_1 % %: 2) -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  0     1            -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  0     1            -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL
  0     1            -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
( 0     1j1  % %: 2) -: qnsign ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
( 0j1   1j_1 % %: 3) -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
( 0j1   1    % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
( 0j1   1    % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL
( 0j1   1    % %: 2) -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
( 0j1   1j1  % %: 3) -: qnsign ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
( 0j_1 _1j_1 % %: 3) -: qnsign (  0        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
( 0j_1 _1    % %: 2) -: qnsign (  0        j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
( 0j_1 _1    % %: 2) -: qnsign (  0        j. -FP_OVFL) ,  -FP_OVFL
( 0j_1 _1    % %: 2) -: qnsign (  0        j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
( 0j_1 _1j1  % %: 3) -: qnsign (  0        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
( 0j_1 _1    % %: 2) -: qnsign (  0        j. -sqr_big) ,  -sqr_big
( 0j_1  0j_1 % %: 2) -: qnsign (  0        j. -sqr_big) ,   0        j. -sqr_big
( 0j_1  0j1  % %: 2) -: qnsign (  0        j. -sqr_big) ,   0        j.  sqr_big
( 0j_1  1    % %: 2) -: qnsign (  0        j. -sqr_big) ,   sqr_big
( 0    _1j_1 % %: 2) -: qnsign (  0        j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  0    _1            -: qnsign (  0        j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  0    _1            -: qnsign (  0        j. -FP_UNFL) ,  -FP_OVFL
  0    _1            -: qnsign (  0        j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
( 0    _1j1  % %: 2) -: qnsign (  0        j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
( 0    _1j_1 % %: 2) -: qnsign (  0        j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  0    _1            -: qnsign (  0        j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  0    _1            -: qnsign (  0        j.  FP_UNFL) ,  -FP_OVFL
  0    _1            -: qnsign (  0        j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
( 0    _1j1  % %: 2) -: qnsign (  0        j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
( 0j1  _1j_1 % %: 3) -: qnsign (  0        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
( 0j1  _1    % %: 2) -: qnsign (  0        j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
( 0j1  _1    % %: 2) -: qnsign (  0        j.  FP_OVFL) ,  -FP_OVFL
( 0j1  _1    % %: 2) -: qnsign (  0        j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
( 0j1  _1j1  % %: 3) -: qnsign (  0        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
( 0j_1  0j_1 % %: 2) -: qnsign (  0        j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  0j_1  0            -: qnsign (  0        j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  0j_1  0            -: qnsign (  0        j. -FP_OVFL) ,  -FP_UNFL
  0j_1  0            -: qnsign (  0        j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
( 0j_1  0j1  % %: 2) -: qnsign (  0        j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  0     0j_1         -: qnsign (  0        j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
( 0j_1 _1j_1 % %: 3) -: qnsign (  0        j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
( 0j_1 _1    % %: 2) -: qnsign (  0        j. -FP_UNFL) ,  -FP_UNFL
( 0j_1 _1j1  % %: 3) -: qnsign (  0        j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  0     0j1          -: qnsign (  0        j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  0     0j_1         -: qnsign (  0        j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
( 0j1  _1j_1 % %: 3) -: qnsign (  0        j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
( 0j1  _1    % %: 2) -: qnsign (  0        j.  FP_UNFL) ,  -FP_UNFL
( 0j1  _1j1  % %: 3) -: qnsign (  0        j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  0     0j1          -: qnsign (  0        j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
( 0j1  _1    % %: 2) -: qnsign (  0        j.  sqr_big) ,  -sqr_big
( 0j1   0j_1 % %: 2) -: qnsign (  0        j.  sqr_big) ,   0        j. -sqr_big
( 0j1   0j1  % %: 2) -: qnsign (  0        j.  sqr_big) ,   0        j.  sqr_big
( 0j1   1    % %: 2) -: qnsign (  0        j.  sqr_big) ,   sqr_big
( 0j1   0j_1 % %: 2) -: qnsign (  0        j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  0j1   0            -: qnsign (  0        j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  0j1   0            -: qnsign (  0        j.  FP_OVFL) ,  -FP_UNFL
  0j1   0            -: qnsign (  0        j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
( 0j1   0j1  % %: 2) -: qnsign (  0        j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
( 0j_1  0j_1 % %: 2) -: qnsign (  0        j. -FP_OVFL) ,   0        j. -FP_OVFL
  0j_1  0            -: qnsign (  0        j. -FP_OVFL) ,   0        j. -FP_UNFL
  0j_1  0            -: qnsign (  0        j. -FP_OVFL) ,   0        j.  FP_UNFL
( 0j_1  0j1  % %: 2) -: qnsign (  0        j. -FP_OVFL) ,   0        j.  FP_OVFL
  0     0j_1         -: qnsign (  0        j. -FP_UNFL) ,   0        j. -FP_OVFL
( 0j_1  0j_1 % %: 2) -: qnsign (  0        j. -FP_UNFL) ,   0        j. -FP_UNFL
( 0j_1  0j1  % %: 2) -: qnsign (  0        j. -FP_UNFL) ,   0        j.  FP_UNFL
  0     0j1          -: qnsign (  0        j. -FP_UNFL) ,   0        j.  FP_OVFL
  0     0j_1         -: qnsign (  0        j.  FP_UNFL) ,   0        j. -FP_OVFL
( 0j1   0j_1 % %: 2) -: qnsign (  0        j.  FP_UNFL) ,   0        j. -FP_UNFL
( 0j1   0j1  % %: 2) -: qnsign (  0        j.  FP_UNFL) ,   0        j.  FP_UNFL
  0     0j1          -: qnsign (  0        j.  FP_UNFL) ,   0        j.  FP_OVFL
( 0j1   0j_1 % %: 2) -: qnsign (  0        j.  FP_OVFL) ,   0        j. -FP_OVFL
  0j1   0            -: qnsign (  0        j.  FP_OVFL) ,   0        j. -FP_UNFL
  0j1   0            -: qnsign (  0        j.  FP_OVFL) ,   0        j.  FP_UNFL
( 0j1   0j1  % %: 2) -: qnsign (  0        j.  FP_OVFL) ,   0        j.  FP_OVFL
( 0j_1  0j_1 % %: 2) -: qnsign (  0        j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  0j_1  0            -: qnsign (  0        j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  0j_1  0            -: qnsign (  0        j. -FP_OVFL) ,   FP_UNFL
  0j_1  0            -: qnsign (  0        j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
( 0j_1  0j1  % %: 2) -: qnsign (  0        j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  0     0j_1         -: qnsign (  0        j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
( 0j_1  1j_1 % %: 3) -: qnsign (  0        j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
( 0j_1  1    % %: 2) -: qnsign (  0        j. -FP_UNFL) ,   FP_UNFL
( 0j_1  1j1  % %: 3) -: qnsign (  0        j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  0     0j1          -: qnsign (  0        j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  0     0j_1         -: qnsign (  0        j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
( 0j1   1j_1 % %: 3) -: qnsign (  0        j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
( 0j1   1    % %: 2) -: qnsign (  0        j.  FP_UNFL) ,   FP_UNFL
( 0j1   1j1  % %: 3) -: qnsign (  0        j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  0     0j1          -: qnsign (  0        j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
( 0j1   0j_1 % %: 2) -: qnsign (  0        j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  0j1   0            -: qnsign (  0        j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  0j1   0            -: qnsign (  0        j.  FP_OVFL) ,   FP_UNFL
  0j1   0            -: qnsign (  0        j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
( 0j1   0j1  % %: 2) -: qnsign (  0        j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
( 0j_1  1j_1 % %: 3) -: qnsign (  0        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
( 0j_1  1    % %: 2) -: qnsign (  0        j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
( 0j_1  1    % %: 2) -: qnsign (  0        j. -FP_OVFL) ,   FP_OVFL
( 0j_1  1    % %: 2) -: qnsign (  0        j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
( 0j_1  1j1  % %: 3) -: qnsign (  0        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
( 0     1j_1 % %: 2) -: qnsign (  0        j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  0     1            -: qnsign (  0        j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  0     1            -: qnsign (  0        j. -FP_UNFL) ,   FP_OVFL
  0     1            -: qnsign (  0        j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
( 0     1j1  % %: 2) -: qnsign (  0        j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
( 0     1j_1 % %: 2) -: qnsign (  0        j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  0     1            -: qnsign (  0        j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  0     1            -: qnsign (  0        j.  FP_UNFL) ,   FP_OVFL
  0     1            -: qnsign (  0        j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
( 0     1j1  % %: 2) -: qnsign (  0        j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
( 0j1   1j_1 % %: 3) -: qnsign (  0        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
( 0j1   1    % %: 2) -: qnsign (  0        j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
( 0j1   1    % %: 2) -: qnsign (  0        j.  FP_OVFL) ,   FP_OVFL
( 0j1   1    % %: 2) -: qnsign (  0        j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
( 0j1   1j1  % %: 3) -: qnsign (  0        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
( 0j_1 _1j_1 % %: 3) -: qnsign (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
( 0j_1 _1    % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
( 0j_1 _1    % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) ,  -FP_OVFL
( 0j_1 _1    % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
( 0j_1 _1j1  % %: 3) -: qnsign (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
( 0    _1j_1 % %: 2) -: qnsign (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  0    _1            -: qnsign (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  0    _1            -: qnsign (  FP_UNFL  j. -FP_UNFL) ,  -FP_OVFL
  0    _1            -: qnsign (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
( 0    _1j1  % %: 2) -: qnsign (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
( 0    _1j_1 % %: 2) -: qnsign    FP_UNFL               , (-FP_OVFL) j. -FP_OVFL
  0    _1            -: qnsign    FP_UNFL               , (-FP_OVFL) j. -FP_UNFL
  0    _1            -: qnsign    FP_UNFL               ,  -FP_OVFL
  0    _1            -: qnsign    FP_UNFL               , (-FP_OVFL) j.  FP_UNFL
( 0    _1j1  % %: 2) -: qnsign    FP_UNFL               , (-FP_OVFL) j.  FP_OVFL
( 0    _1j_1 % %: 2) -: qnsign (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  0    _1            -: qnsign (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  0    _1            -: qnsign (  FP_UNFL  j.  FP_UNFL) ,  -FP_OVFL
  0    _1            -: qnsign (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
( 0    _1j1  % %: 2) -: qnsign (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
( 0j1  _1j_1 % %: 3) -: qnsign (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
( 0j1  _1    % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
( 0j1  _1    % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) ,  -FP_OVFL
( 0j1  _1    % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
( 0j1  _1j1  % %: 3) -: qnsign (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
( 0j_1  0j_1 % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  0j_1  0            -: qnsign (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  0j_1  0            -: qnsign (  FP_UNFL  j. -FP_OVFL) ,  -FP_UNFL
  0j_1  0            -: qnsign (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
( 0j_1  0j1  % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  0     0j_1         -: qnsign (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(_1j1   1j1  %   _2) -: qnsign (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
( 1j_1 _1    % %: 3) -: qnsign (  FP_UNFL  j. -FP_UNFL) ,  -FP_UNFL
( 1j_1 _1j1  %    2) -: qnsign (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  0     0j1          -: qnsign (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  0     0j_1         -: qnsign    FP_UNFL               , (-FP_UNFL) j. -FP_OVFL
( 1    _1j_1 % %: 3) -: qnsign    FP_UNFL               , (-FP_UNFL) j. -FP_UNFL
( 1    _1    % %: 2) -: qnsign    FP_UNFL               ,  -FP_UNFL
( 1    _1j1  % %: 3) -: qnsign    FP_UNFL               , (-FP_UNFL) j.  FP_UNFL
  0     0j1          -: qnsign    FP_UNFL               , (-FP_UNFL) j.  FP_OVFL
  0     0j_1         -: qnsign (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
( 1j1  _1j_1 %    2) -: qnsign (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
( 1j1  _1    % %: 3) -: qnsign (  FP_UNFL  j.  FP_UNFL) ,  -FP_UNFL
( 1j1  _1j1  %    2) -: qnsign (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  0     0j1          -: qnsign (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
( 0j1   0j_1 % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  0j1   0            -: qnsign (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  0j1   0            -: qnsign (  FP_UNFL  j.  FP_OVFL) ,  -FP_UNFL
  0j1   0            -: qnsign (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
( 0j1   0j1  % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
( 0j_1  0j_1 % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   0        j. -FP_OVFL
  0j_1  0            -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   0        j. -FP_UNFL
  0j_1  0            -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   0        j.  FP_UNFL
( 0j_1  0j1  % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   0        j.  FP_OVFL
  0     0j_1         -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   0        j. -FP_OVFL
( 1j_1  0j_1 % %: 3) -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   0        j. -FP_UNFL
( 1j_1  0j1  % %: 3) -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   0        j.  FP_UNFL
  0     0j1          -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   0        j.  FP_OVFL
  0     0j_1         -: qnsign    FP_UNFL               ,   0        j. -FP_OVFL
( 1     0j_1 % %: 2) -: qnsign    FP_UNFL               ,   0        j. -FP_UNFL
( 1     0j1  % %: 2) -: qnsign    FP_UNFL               ,   0        j.  FP_UNFL
  0     0j1          -: qnsign    FP_UNFL               ,   0        j.  FP_OVFL
  0     0j_1         -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   0        j. -FP_OVFL
( 1j1   0j_1 % %: 3) -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   0        j. -FP_UNFL
( 1j1   0j1  % %: 3) -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   0        j.  FP_UNFL
  0     0j1          -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   0        j.  FP_OVFL
( 0j1   0j_1 % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   0        j. -FP_OVFL
  0j1   0            -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   0        j. -FP_UNFL
  0j1   0            -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   0        j.  FP_UNFL
( 0j1   0j1  % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   0        j.  FP_OVFL
( 0j_1  0j_1 % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  0j_1  0            -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  0j_1  0            -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL
  0j_1  0            -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
( 0j_1  0j1  % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  0     0j_1         -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
( 1j_1  1j_1 %    2) -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
( 1j_1  1    % %: 3) -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL
( 1j_1  1j1  %    2) -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  0     0j1          -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  0     0j_1         -: qnsign    FP_UNFL               ,   FP_UNFL  j. -FP_OVFL
( 1     1j_1 % %: 3) -: qnsign    FP_UNFL               ,   FP_UNFL  j. -FP_UNFL
( 1     1    % %: 2) -: qnsign    FP_UNFL               ,   FP_UNFL
( 1     1j1  % %: 3) -: qnsign    FP_UNFL               ,   FP_UNFL  j.  FP_UNFL
  0     0j1          -: qnsign    FP_UNFL               ,   FP_UNFL  j.  FP_OVFL
  0     0j_1         -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
( 1j1   1j_1 %    2) -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
( 1j1   1    % %: 3) -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL
( 1j1   1j1  %    2) -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  0     0j1          -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
( 0j1   0j_1 % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  0j1   0            -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  0j1   0            -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL
  0j1   0            -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
( 0j1   0j1  % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
( 0j_1  1j_1 % %: 3) -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
( 0j_1  1    % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
( 0j_1  1    % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL
( 0j_1  1    % %: 2) -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
( 0j_1  1j1  % %: 3) -: qnsign (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
( 0     1j_1 % %: 2) -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  0     1            -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  0     1            -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL
  0     1            -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
( 0     1j1  % %: 2) -: qnsign (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
( 0     1j_1 % %: 2) -: qnsign    FP_UNFL               ,   FP_OVFL  j. -FP_OVFL
  0     1            -: qnsign    FP_UNFL               ,   FP_OVFL  j. -FP_UNFL
  0     1            -: qnsign    FP_UNFL               ,   FP_OVFL
  0     1            -: qnsign    FP_UNFL               ,   FP_OVFL  j.  FP_UNFL
( 0     1j1  % %: 2) -: qnsign    FP_UNFL               ,   FP_OVFL  j.  FP_OVFL
( 0     1j_1 % %: 2) -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  0     1            -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  0     1            -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL
  0     1            -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
( 0     1j1  % %: 2) -: qnsign (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
( 0j1   1j_1 % %: 3) -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
( 0j1   1    % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
( 0j1   1    % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL
( 0j1   1    % %: 2) -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
( 0j1   1j1  % %: 3) -: qnsign (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(_1j1   1j1  %   _2) -: qnsign (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
( 1j_1 _1    % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
( 1j_1 _1    % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,  -FP_OVFL
( 1j_1 _1    % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
( 1j_1 _1j1  %    2) -: qnsign (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
( 1    _1j_1 % %: 3) -: qnsign (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
( 1    _1    % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
( 1    _1    % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) ,  -FP_OVFL
( 1    _1    % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
( 1    _1j1  % %: 3) -: qnsign (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
( 1    _1j_1 % %: 3) -: qnsign    FP_OVFL               , (-FP_OVFL) j. -FP_OVFL
( 1    _1    % %: 2) -: qnsign    FP_OVFL               , (-FP_OVFL) j. -FP_UNFL
( 1    _1    % %: 2) -: qnsign    FP_OVFL               ,  -FP_OVFL
( 1    _1    % %: 2) -: qnsign    FP_OVFL               , (-FP_OVFL) j.  FP_UNFL
( 1    _1j1  % %: 3) -: qnsign    FP_OVFL               , (-FP_OVFL) j.  FP_OVFL
( 1    _1j_1 % %: 3) -: qnsign (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
( 1    _1    % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
( 1    _1    % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) ,  -FP_OVFL
( 1    _1    % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
( 1    _1j1  % %: 3) -: qnsign (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
( 1j1  _1j_1 %    2) -: qnsign (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
( 1j1  _1    % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
( 1j1  _1    % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,  -FP_OVFL
( 1j1  _1    % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
( 1j1  _1j1  %    2) -: qnsign (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
( 1j_1  0j_1 % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
( 1j_1  0    % %: 2) -: qnsign (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
( 1j_1  0    % %: 2) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,  -FP_UNFL
( 1j_1  0    % %: 2) -: qnsign (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
( 1j_1  0j1  % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
( 1     0j_1 % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  1     0            -: qnsign (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  1     0            -: qnsign (  FP_OVFL  j. -FP_UNFL) ,  -FP_UNFL
  1     0            -: qnsign (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
( 1     0j1  % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
( 1     0j_1 % %: 2) -: qnsign    FP_OVFL               , (-FP_UNFL) j. -FP_OVFL
  1     0            -: qnsign    FP_OVFL               , (-FP_UNFL) j. -FP_UNFL
  1     0            -: qnsign    FP_OVFL               ,  -FP_UNFL
  1     0            -: qnsign    FP_OVFL               , (-FP_UNFL) j.  FP_UNFL
( 1     0j1  % %: 2) -: qnsign    FP_OVFL               , (-FP_UNFL) j.  FP_OVFL
( 1     0j_1 % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  1     0            -: qnsign (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  1     0            -: qnsign (  FP_OVFL  j.  FP_UNFL) ,  -FP_UNFL
  1     0            -: qnsign (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
( 1     0j1  % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
( 1j1   0j_1 % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
( 1j1   0    % %: 2) -: qnsign (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
( 1j1   0    % %: 2) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,  -FP_UNFL
( 1j1   0    % %: 2) -: qnsign (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
( 1j1   0j1  % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
( 1j_1  0j_1 % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   0        j. -FP_OVFL
( 1j_1  0    % %: 2) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   0        j. -FP_UNFL
( 1j_1  0    % %: 2) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   0        j.  FP_UNFL
( 1j_1  0j1  % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   0        j.  FP_OVFL
( 1     0j_1 % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   0        j. -FP_OVFL
  1     0            -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   0        j. -FP_UNFL
  1     0            -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   0        j.  FP_UNFL
( 1     0j1  % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   0        j.  FP_OVFL
( 1     0j_1 % %: 2) -: qnsign    FP_OVFL               ,   0        j. -FP_OVFL
  1     0            -: qnsign    FP_OVFL               ,   0        j. -FP_UNFL
  1     0            -: qnsign    FP_OVFL               ,   0        j.  FP_UNFL
( 1     0j1  % %: 2) -: qnsign    FP_OVFL               ,   0        j.  FP_OVFL
( 1     0j_1 % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   0        j. -FP_OVFL
  1     0            -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   0        j. -FP_UNFL
  1     0            -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   0        j.  FP_UNFL
( 1     0j1  % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   0        j.  FP_OVFL
( 1j1   0j_1 % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   0        j. -FP_OVFL
( 1j1   0    % %: 2) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   0        j. -FP_UNFL
( 1j1   0    % %: 2) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   0        j.  FP_UNFL
( 1j1   0j1  % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   0        j.  FP_OVFL
( 1j_1  0j_1 % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
( 1j_1  0    % %: 2) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
( 1j_1  0    % %: 2) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL
( 1j_1  0    % %: 2) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
( 1j_1  0j1  % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
( 1     0j_1 % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  1     0            -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  1     0            -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL
  1     0            -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
( 1     0j1  % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
( 1     0j_1 % %: 2) -: qnsign    FP_OVFL               ,   FP_UNFL  j. -FP_OVFL
  1     0            -: qnsign    FP_OVFL               ,   FP_UNFL  j. -FP_UNFL
  1     0            -: qnsign    FP_OVFL               ,   FP_UNFL
  1     0            -: qnsign    FP_OVFL               ,   FP_UNFL  j.  FP_UNFL
( 1     0j1  % %: 2) -: qnsign    FP_OVFL               ,   FP_UNFL  j.  FP_OVFL
( 1     0j_1 % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  1     0            -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  1     0            -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL
  1     0            -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
( 1     0j1  % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
( 1j_1 _1j_1 %    2) -: qnsign (  sqr_big  j. -sqr_big) , (-sqr_big) j. -sqr_big
( 1j_1 _1j1  %    2) -: qnsign (  sqr_big  j. -sqr_big) , (-sqr_big) j.  sqr_big
( 1j_1  1j_1 %    2) -: qnsign (  sqr_big  j. -sqr_big) ,   sqr_big  j. -sqr_big
( 1j_1  1j1  %    2) -: qnsign (  sqr_big  j. -sqr_big) ,   sqr_big  j.  sqr_big
( 1    _1    % %: 2) -: qnsign    sqr_big               ,  -sqr_big
( 1     0j_1 % %: 2) -: qnsign    sqr_big               ,   0        j. -sqr_big
( 1     0j1  % %: 2) -: qnsign    sqr_big               ,   0        j.  sqr_big
( 1     1    % %: 2) -: qnsign    sqr_big               ,   sqr_big
( 1j1  _1j_1 %    2) -: qnsign (  sqr_big  j.  sqr_big) , (-sqr_big) j. -sqr_big
( 1j1  _1j1  %    2) -: qnsign (  sqr_big  j.  sqr_big) , (-sqr_big) j.  sqr_big
( 1j1   1j_1 %    2) -: qnsign (  sqr_big  j.  sqr_big) ,   sqr_big  j. -sqr_big
( 1j1   1j1  %    2) -: qnsign (  sqr_big  j.  sqr_big) ,   sqr_big  j.  sqr_big
( 1j1   0j_1 % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
( 1j1   0    % %: 2) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
( 1j1   0    % %: 2) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL
( 1j1   0    % %: 2) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
( 1j1   0j1  % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
( 1j_1  1j_1 %    2) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
( 1j_1  1    % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
( 1j_1  1    % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL
( 1j_1  1    % %: 3) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
( 1j_1  1j1  %    2) -: qnsign (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
( 1     1j_1 % %: 3) -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
( 1     1    % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
( 1     1    % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL
( 1     1    % %: 2) -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
( 1     1j1  % %: 3) -: qnsign (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
( 1     1j_1 % %: 3) -: qnsign    FP_OVFL               ,   FP_OVFL  j. -FP_OVFL
( 1     1    % %: 2) -: qnsign    FP_OVFL               ,   FP_OVFL  j. -FP_UNFL
( 1     1    % %: 2) -: qnsign    FP_OVFL               ,   FP_OVFL
( 1     1    % %: 2) -: qnsign    FP_OVFL               ,   FP_OVFL  j.  FP_UNFL
( 1     1j1  % %: 3) -: qnsign    FP_OVFL               ,   FP_OVFL  j.  FP_OVFL
( 1     1j_1 % %: 3) -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
( 1     1    % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
( 1     1    % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL
( 1     1    % %: 2) -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
( 1     1j1  % %: 3) -: qnsign (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
( 1j1   1j_1 %    2) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
( 1j1   1    % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
( 1j1   1    % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL
( 1j1   1    % %: 3) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
( 1j1   1j1  %    2) -: qnsign (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
NB. - input without edge cases
( 1j1   1j1  %   _2) -: qnsign _1j_1 _1j_1
(_1j_1 _1    % %: 3) -: qnsign _1j_1 _1
(_1j_1 _1j1  %    2) -: qnsign _1j_1 _1j1
(_1j_1  0j_1 % %: 3) -: qnsign _1j_1  0j_1
(_1j_1  0j1  % %: 3) -: qnsign _1j_1  0j1
(_1j_1  1j_1 %    2) -: qnsign _1j_1  1j_1
(_1j_1  1    % %: 3) -: qnsign _1j_1  1
(_1j_1  1j1  %    2) -: qnsign _1j_1  1j1
(_1    _1j_1 % %: 3) -: qnsign _1    _1j_1
(_1    _1    % %: 2) -: qnsign _1    _1
(_1    _1j1  % %: 3) -: qnsign _1    _1j1
(_1     0j_1 % %: 2) -: qnsign _1     0j_1
(_1     0j1  % %: 2) -: qnsign _1     0j1
(_1     1j_1 % %: 3) -: qnsign _1     1j_1
(_1     1    % %: 2) -: qnsign _1     1
(_1     1j1  % %: 3) -: qnsign _1     1j1
(_1j1  _1j_1 %    2) -: qnsign _1j1  _1j_1
(_1j1  _1    % %: 3) -: qnsign _1j1  _1
(_1j1  _1j1  %    2) -: qnsign _1j1  _1j1
(_1j1   0j_1 % %: 3) -: qnsign _1j1   0j_1
(_1j1   0j1  % %: 3) -: qnsign _1j1   0j1
(_1j1   1j_1 %    2) -: qnsign _1j1   1j_1
(_1j1   1    % %: 3) -: qnsign _1j1   1
(_1j1   1j1  %    2) -: qnsign _1j1   1j1
( 0j_1 _1j_1 % %: 3) -: qnsign  0j_1 _1j_1
( 0j_1 _1    % %: 2) -: qnsign  0j_1 _1
( 0j_1 _1j1  % %: 3) -: qnsign  0j_1 _1j1
( 0j_1  0j_1 % %: 2) -: qnsign  0j_1  0j_1
( 0j_1  0j1  % %: 2) -: qnsign  0j_1  0j1
( 0j_1  1j_1 % %: 3) -: qnsign  0j_1  1j_1
( 0j_1  1    % %: 2) -: qnsign  0j_1  1
( 0j_1  1j1  % %: 3) -: qnsign  0j_1  1j1
( 0j1  _1j_1 % %: 3) -: qnsign  0j1  _1j_1
( 0j1  _1    % %: 2) -: qnsign  0j1  _1
( 0j1  _1j1  % %: 3) -: qnsign  0j1  _1j1
( 0j1   0j_1 % %: 2) -: qnsign  0j1   0j_1
( 0j1   0j1  % %: 2) -: qnsign  0j1   0j1
( 0j1   1j_1 % %: 3) -: qnsign  0j1   1j_1
( 0j1   1    % %: 2) -: qnsign  0j1   1
( 0j1   1j1  % %: 3) -: qnsign  0j1   1j1
( 1j_1 _1j_1 %    2) -: qnsign  1j_1 _1j_1
( 1j_1 _1    % %: 3) -: qnsign  1j_1 _1
( 1j_1 _1j1  %    2) -: qnsign  1j_1 _1j1
( 1j_1  0j_1 % %: 3) -: qnsign  1j_1  0j_1
( 1j_1  0j1  % %: 3) -: qnsign  1j_1  0j1
( 1j_1  1j_1 %    2) -: qnsign  1j_1  1j_1
( 1j_1  1    % %: 3) -: qnsign  1j_1  1
( 1j_1  1j1  %    2) -: qnsign  1j_1  1j1
( 1    _1j_1 % %: 3) -: qnsign  1    _1j_1
( 1    _1    % %: 2) -: qnsign  1    _1
( 1    _1j1  % %: 3) -: qnsign  1    _1j1
( 1     0j_1 % %: 2) -: qnsign  1     0j_1
( 1     0j1  % %: 2) -: qnsign  1     0j1
( 1     1j_1 % %: 3) -: qnsign  1     1j_1
( 1     1    % %: 2) -: qnsign  1     1
( 1     1j1  % %: 3) -: qnsign  1     1j1
( 1j1  _1j_1 %    2) -: qnsign  1j1  _1j_1
( 1j1  _1    % %: 3) -: qnsign  1j1  _1
( 1j1  _1j1  %    2) -: qnsign  1j1  _1j1
( 1j1   0j_1 % %: 3) -: qnsign  1j1   0j_1
( 1j1   0j1  % %: 3) -: qnsign  1j1   0j1
( 1j1   1j_1 %    2) -: qnsign  1j1   1j_1
( 1j1   1    % %: 3) -: qnsign  1j1   1
( 1j1   1j1  %    2) -: qnsign  1j1   1j1
(2 2 $ j.~ 1r2) -: qnsign q3 ,: q4

NB. qninv
NB. - input contains NaN
1 1 -: isnan qninv  0     0j_.
1 1 -: isnan qninv  0    _.
1 1 -: isnan qninv  0j_.  0
1 1 -: isnan qninv _.     0
1 1 -: isnan qninv  0    _.j_.
1 1 -: isnan qninv  0j_.  0j_.
1 1 -: isnan qninv  0j_. _.
1 1 -: isnan qninv _.     0j_.
1 1 -: isnan qninv _.    _.
1 1 -: isnan qninv _.j_.  0
1 1 -: isnan qninv _.j_. _.
1 1 -: isnan qninv _.j_.  0j_.
1 1 -: isnan qninv _.    _.j_.
1 1 -: isnan qninv  0j_. _.j_.
1 1 -: isnan qninv _.j_. _.j_.
(2 2 $ 1) -: isnan qninv _. 0 ,: 0 _.
(,.~ 1 0) -: isnan qninv _. 0 ,: q4
(,.~ 0 1) -: isnan qninv q4   ,: 0 _.
NB. - result is undefined so NaN error must be throwed
0:@qninv :: 1: 0 0
NB. - input contains q∞ so NaN error must be throwed
0:@qninv :: 1: __j__ __j__
0:@qninv :: 1: __j__ __j_1
0:@qninv :: 1: __j__ __
0:@qninv :: 1: __j__ __j1
0:@qninv :: 1: __j__ __j_
0:@qninv :: 1: __j__ _1j__
0:@qninv :: 1: __j__ _1j_1
0:@qninv :: 1: __j__ _1
0:@qninv :: 1: __j__ _1j1
0:@qninv :: 1: __j__ _1j_
0:@qninv :: 1: __j__  0j__
0:@qninv :: 1: __j__  0j_1
0:@qninv :: 1: __j__  0
0:@qninv :: 1: __j__  0j1
0:@qninv :: 1: __j__  0j_
0:@qninv :: 1: __j__  1j__
0:@qninv :: 1: __j__  1j_1
0:@qninv :: 1: __j__  1
0:@qninv :: 1: __j__  1j1
0:@qninv :: 1: __j__  1j_
0:@qninv :: 1: __j__  _j__
0:@qninv :: 1: __j__  _j_1
0:@qninv :: 1: __j__  _
0:@qninv :: 1: __j__  _j1
0:@qninv :: 1: __j__  _j_
0:@qninv :: 1: __j_1 __j__
0:@qninv :: 1: __j_1 __j_1
0:@qninv :: 1: __j_1 __
0:@qninv :: 1: __j_1 __j1
0:@qninv :: 1: __j_1 __j_
0:@qninv :: 1: __j_1 _1j__
0:@qninv :: 1: __j_1 _1j_
0:@qninv :: 1: __j_1  0j__
0:@qninv :: 1: __j_1  0j_
0:@qninv :: 1: __j_1  1j__
0:@qninv :: 1: __j_1  1j_
0:@qninv :: 1: __j_1  _j__
0:@qninv :: 1: __j_1  _j_1
0:@qninv :: 1: __j_1  _
0:@qninv :: 1: __j_1  _j1
0:@qninv :: 1: __j_1  _j_
0:@qninv :: 1: __    __j__
0:@qninv :: 1: __    __j_1
0:@qninv :: 1: __    __
0:@qninv :: 1: __    __j1
0:@qninv :: 1: __    __j_
0:@qninv :: 1: __    _1j__
0:@qninv :: 1: __    _1j_
0:@qninv :: 1: __     0j__
0:@qninv :: 1: __     0j_
0:@qninv :: 1: __     1j__
0:@qninv :: 1: __     1j_
0:@qninv :: 1: __     _j__
0:@qninv :: 1: __     _j_1
0:@qninv :: 1: __     _
0:@qninv :: 1: __     _j1
0:@qninv :: 1: __     _j_
0:@qninv :: 1: __j1  __j__
0:@qninv :: 1: __j1  __j_1
0:@qninv :: 1: __j1  __
0:@qninv :: 1: __j1  __j1
0:@qninv :: 1: __j1  __j_
0:@qninv :: 1: __j1  _1j__
0:@qninv :: 1: __j1  _1j_
0:@qninv :: 1: __j1   0j__
0:@qninv :: 1: __j1   0j_
0:@qninv :: 1: __j1   1j__
0:@qninv :: 1: __j1   1j_
0:@qninv :: 1: __j1   _j__
0:@qninv :: 1: __j1   _j_1
0:@qninv :: 1: __j1   _
0:@qninv :: 1: __j1   _j1
0:@qninv :: 1: __j1   _j_
0:@qninv :: 1: __j_  __j__
0:@qninv :: 1: __j_  __j_1
0:@qninv :: 1: __j_  __
0:@qninv :: 1: __j_  __j1
0:@qninv :: 1: __j_  __j_
0:@qninv :: 1: __j_  _1j__
0:@qninv :: 1: __j_  _1j_1
0:@qninv :: 1: __j_  _1
0:@qninv :: 1: __j_  _1j1
0:@qninv :: 1: __j_  _1j_
0:@qninv :: 1: __j_   0j__
0:@qninv :: 1: __j_   0j_1
0:@qninv :: 1: __j_   0
0:@qninv :: 1: __j_   0j1
0:@qninv :: 1: __j_   0j_
0:@qninv :: 1: __j_   1j__
0:@qninv :: 1: __j_   1j_1
0:@qninv :: 1: __j_   1
0:@qninv :: 1: __j_   1j1
0:@qninv :: 1: __j_   1j_
0:@qninv :: 1: __j_   _j__
0:@qninv :: 1: __j_   _j_1
0:@qninv :: 1: __j_   _
0:@qninv :: 1: __j_   _j1
0:@qninv :: 1: __j_   _j_
0:@qninv :: 1: _1j__ __j__
0:@qninv :: 1: _1j__ __j_1
0:@qninv :: 1: _1j__ __
0:@qninv :: 1: _1j__ __j1
0:@qninv :: 1: _1j__ __j_
0:@qninv :: 1: _1j__ _1j__
0:@qninv :: 1: _1j__ _1j_
0:@qninv :: 1: _1j__  0j__
0:@qninv :: 1: _1j__  0j_
0:@qninv :: 1: _1j__  1j__
0:@qninv :: 1: _1j__  1j_
0:@qninv :: 1: _1j__  _j__
0:@qninv :: 1: _1j__  _j_1
0:@qninv :: 1: _1j__  _
0:@qninv :: 1: _1j__  _j1
0:@qninv :: 1: _1j__  _j_
0:@qninv :: 1: _1j_1 __j__
0:@qninv :: 1: _1j_1 __j_
0:@qninv :: 1: _1j_1  _j__
0:@qninv :: 1: _1j_1  _j_
0:@qninv :: 1: _1    __j__
0:@qninv :: 1: _1    __j_
0:@qninv :: 1: _1     _j__
0:@qninv :: 1: _1     _j_
0:@qninv :: 1: _1j1  __j__
0:@qninv :: 1: _1j1  __j_
0:@qninv :: 1: _1j1   _j__
0:@qninv :: 1: _1j1   _j_
0:@qninv :: 1: _1j_  __j__
0:@qninv :: 1: _1j_  __j_1
0:@qninv :: 1: _1j_  __
0:@qninv :: 1: _1j_  __j1
0:@qninv :: 1: _1j_  __j_
0:@qninv :: 1: _1j_  _1j__
0:@qninv :: 1: _1j_  _1j_
0:@qninv :: 1: _1j_   0j__
0:@qninv :: 1: _1j_   0j_
0:@qninv :: 1: _1j_   1j__
0:@qninv :: 1: _1j_   1j_
0:@qninv :: 1: _1j_   _j__
0:@qninv :: 1: _1j_   _j_1
0:@qninv :: 1: _1j_   _
0:@qninv :: 1: _1j_   _j1
0:@qninv :: 1: _1j_   _j_
0:@qninv :: 1:  0j__ __j__
0:@qninv :: 1:  0j__ __j_1
0:@qninv :: 1:  0j__ __
0:@qninv :: 1:  0j__ __j1
0:@qninv :: 1:  0j__ __j_
0:@qninv :: 1:  0j__ _1j__
0:@qninv :: 1:  0j__ _1j_
0:@qninv :: 1:  0j__  0j__
0:@qninv :: 1:  0j__  0j_
0:@qninv :: 1:  0j__  1j__
0:@qninv :: 1:  0j__  1j_
0:@qninv :: 1:  0j__  _j__
0:@qninv :: 1:  0j__  _j_1
0:@qninv :: 1:  0j__  _
0:@qninv :: 1:  0j__  _j1
0:@qninv :: 1:  0j__  _j_
0:@qninv :: 1:  0j_1 __j__
0:@qninv :: 1:  0j_1 __j_
0:@qninv :: 1:  0j_1  _j__
0:@qninv :: 1:  0j_1  _j_
0:@qninv :: 1:  0    __j__
0:@qninv :: 1:  0    __j_
0:@qninv :: 1:  0     _j__
0:@qninv :: 1:  0     _j_
0:@qninv :: 1:  0j1  __j__
0:@qninv :: 1:  0j1  __j_
0:@qninv :: 1:  0j1   _j__
0:@qninv :: 1:  0j1   _j_
0:@qninv :: 1:  0j_  __j__
0:@qninv :: 1:  0j_  __j_1
0:@qninv :: 1:  0j_  __
0:@qninv :: 1:  0j_  __j1
0:@qninv :: 1:  0j_  __j_
0:@qninv :: 1:  0j_  _1j__
0:@qninv :: 1:  0j_  _1j_
0:@qninv :: 1:  0j_   0j__
0:@qninv :: 1:  0j_   0j_
0:@qninv :: 1:  0j_   1j__
0:@qninv :: 1:  0j_   1j_
0:@qninv :: 1:  0j_   _j__
0:@qninv :: 1:  0j_   _j_1
0:@qninv :: 1:  0j_   _
0:@qninv :: 1:  0j_   _j1
0:@qninv :: 1:  0j_   _j_
0:@qninv :: 1:  1j__ __j__
0:@qninv :: 1:  1j__ __j_1
0:@qninv :: 1:  1j__ __
0:@qninv :: 1:  1j__ __j1
0:@qninv :: 1:  1j__ __j_
0:@qninv :: 1:  1j__ _1j__
0:@qninv :: 1:  1j__ _1j_
0:@qninv :: 1:  1j__  0j__
0:@qninv :: 1:  1j__  0j_
0:@qninv :: 1:  1j__  1j__
0:@qninv :: 1:  1j__  1j_
0:@qninv :: 1:  1j__  _j__
0:@qninv :: 1:  1j__  _j_1
0:@qninv :: 1:  1j__  _
0:@qninv :: 1:  1j__  _j1
0:@qninv :: 1:  1j__  _j_
0:@qninv :: 1:  1j_1 __j__
0:@qninv :: 1:  1j_1 __j_
0:@qninv :: 1:  1j_1  _j__
0:@qninv :: 1:  1j_1  _j_
0:@qninv :: 1:  1    __j__
0:@qninv :: 1:  1    __j_
0:@qninv :: 1:  1     _j__
0:@qninv :: 1:  1     _j_
0:@qninv :: 1:  1j1  __j__
0:@qninv :: 1:  1j1  __j_
0:@qninv :: 1:  1j1   _j__
0:@qninv :: 1:  1j1   _j_
0:@qninv :: 1:  1j_  __j__
0:@qninv :: 1:  1j_  __j_1
0:@qninv :: 1:  1j_  __
0:@qninv :: 1:  1j_  __j1
0:@qninv :: 1:  1j_  __j_
0:@qninv :: 1:  1j_  _1j__
0:@qninv :: 1:  1j_  _1j_
0:@qninv :: 1:  1j_   0j__
0:@qninv :: 1:  1j_   0j_
0:@qninv :: 1:  1j_   1j__
0:@qninv :: 1:  1j_   1j_
0:@qninv :: 1:  1j_   _j__
0:@qninv :: 1:  1j_   _j_1
0:@qninv :: 1:  1j_   _
0:@qninv :: 1:  1j_   _j1
0:@qninv :: 1:  1j_   _j_
0:@qninv :: 1:  _j__ __j__
0:@qninv :: 1:  _j__ __j_1
0:@qninv :: 1:  _j__ __
0:@qninv :: 1:  _j__ __j1
0:@qninv :: 1:  _j__ __j_
0:@qninv :: 1:  _j__ _1j__
0:@qninv :: 1:  _j__ _1j_1
0:@qninv :: 1:  _j__ _1
0:@qninv :: 1:  _j__ _1j1
0:@qninv :: 1:  _j__ _1j_
0:@qninv :: 1:  _j__  0j__
0:@qninv :: 1:  _j__  0j_1
0:@qninv :: 1:  _j__  0
0:@qninv :: 1:  _j__  0j1
0:@qninv :: 1:  _j__  0j_
0:@qninv :: 1:  _j__  1j__
0:@qninv :: 1:  _j__  1j_1
0:@qninv :: 1:  _j__  1
0:@qninv :: 1:  _j__  1j1
0:@qninv :: 1:  _j__  1j_
0:@qninv :: 1:  _j__  _j__
0:@qninv :: 1:  _j__  _j_1
0:@qninv :: 1:  _j__  _
0:@qninv :: 1:  _j__  _j1
0:@qninv :: 1:  _j__  _j_
0:@qninv :: 1:  _j_1 __j__
0:@qninv :: 1:  _j_1 __j_1
0:@qninv :: 1:  _j_1 __
0:@qninv :: 1:  _j_1 __j1
0:@qninv :: 1:  _j_1 __j_
0:@qninv :: 1:  _j_1 _1j__
0:@qninv :: 1:  _j_1 _1j_
0:@qninv :: 1:  _j_1  0j__
0:@qninv :: 1:  _j_1  0j_
0:@qninv :: 1:  _j_1  1j__
0:@qninv :: 1:  _j_1  1j_
0:@qninv :: 1:  _j_1  _j__
0:@qninv :: 1:  _j_1  _j_1
0:@qninv :: 1:  _j_1  _
0:@qninv :: 1:  _j_1  _j1
0:@qninv :: 1:  _j_1  _j_
0:@qninv :: 1:  _    __j__
0:@qninv :: 1:  _    __j_1
0:@qninv :: 1:  _    __
0:@qninv :: 1:  _    __j1
0:@qninv :: 1:  _    __j_
0:@qninv :: 1:  _    _1j__
0:@qninv :: 1:  _    _1j_
0:@qninv :: 1:  _     0j__
0:@qninv :: 1:  _     0j_
0:@qninv :: 1:  _     1j__
0:@qninv :: 1:  _     1j_
0:@qninv :: 1:  _     _j__
0:@qninv :: 1:  _     _j_1
0:@qninv :: 1:  _     _
0:@qninv :: 1:  _     _j1
0:@qninv :: 1:  _     _j_
0:@qninv :: 1:  _j1  __j__
0:@qninv :: 1:  _j1  __j_1
0:@qninv :: 1:  _j1  __
0:@qninv :: 1:  _j1  __j1
0:@qninv :: 1:  _j1  __j_
0:@qninv :: 1:  _j1  _1j__
0:@qninv :: 1:  _j1  _1j_
0:@qninv :: 1:  _j1   0j__
0:@qninv :: 1:  _j1   0j_
0:@qninv :: 1:  _j1   1j__
0:@qninv :: 1:  _j1   1j_
0:@qninv :: 1:  _j1   _j__
0:@qninv :: 1:  _j1   _j_1
0:@qninv :: 1:  _j1   _
0:@qninv :: 1:  _j1   _j1
0:@qninv :: 1:  _j1   _j_
0:@qninv :: 1:  _j_  __j__
0:@qninv :: 1:  _j_  __j_1
0:@qninv :: 1:  _j_  __
0:@qninv :: 1:  _j_  __j1
0:@qninv :: 1:  _j_  __j_
0:@qninv :: 1:  _j_  _1j__
0:@qninv :: 1:  _j_  _1j_1
0:@qninv :: 1:  _j_  _1
0:@qninv :: 1:  _j_  _1j1
0:@qninv :: 1:  _j_  _1j_
0:@qninv :: 1:  _j_   0j__
0:@qninv :: 1:  _j_   0j_1
0:@qninv :: 1:  _j_   0
0:@qninv :: 1:  _j_   0j1
0:@qninv :: 1:  _j_   0j_
0:@qninv :: 1:  _j_   1j__
0:@qninv :: 1:  _j_   1j_1
0:@qninv :: 1:  _j_   1
0:@qninv :: 1:  _j_   1j1
0:@qninv :: 1:  _j_   1j_
0:@qninv :: 1:  _j_   _j__
0:@qninv :: 1:  _j_   _j_1
0:@qninv :: 1:  _j_   _
0:@qninv :: 1:  _j_   _j1
0:@qninv :: 1:  _j_   _j_
NB. - input contains directed infinity and is not trivial
0 0 -: qninv ( __        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv ( __        j. -FP_OVFL) ,  -FP_OVFL
0 0 -: qninv ( __        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv ( __        j. -FP_OVFL) ,   0        j. -FP_OVFL
0 0 -: qninv ( __        j. -FP_OVFL) ,   0        j.  FP_OVFL
0 0 -: qninv ( __        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv ( __        j. -FP_OVFL) ,   FP_OVFL
0 0 -: qninv ( __        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
0 0 -: qninv   __                     , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv   __                     ,  -FP_OVFL
0 0 -: qninv   __                     , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv   __                     ,   0        j. -FP_OVFL
0 0 -: qninv   __                     ,   0        j.  FP_OVFL
0 0 -: qninv   __                     ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv   __                     ,   FP_OVFL
0 0 -: qninv   __                     ,   FP_OVFL  j.  FP_OVFL
0 0 -: qninv ( __        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv ( __        j.  FP_OVFL) ,  -FP_OVFL
0 0 -: qninv ( __        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv ( __        j.  FP_OVFL) ,   0        j. -FP_OVFL
0 0 -: qninv ( __        j.  FP_OVFL) ,   0        j.  FP_OVFL
0 0 -: qninv ( __        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv ( __        j.  FP_OVFL) ,   FP_OVFL
0 0 -: qninv ( __        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. __      ) , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. __      ) ,  -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. __      ) , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. __      ) ,   0        j. -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. __      ) ,   0        j.  FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. __      ) ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. __      ) ,   FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. __      ) ,   FP_OVFL  j.  FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) ,  __        j. -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) ,  __
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) ,  __        j.  FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. __
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  _
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   0        j. __
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  _
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. __
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  _
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   _        j. -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   _
0 0 -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   _        j.  FP_OVFL
0 0 -: qninv  (-FP_OVFL)              ,  __        j. -FP_OVFL
0 0 -: qninv  (-FP_OVFL)              ,  __
0 0 -: qninv  (-FP_OVFL)              ,  __        j.  FP_OVFL
0 0 -: qninv  (-FP_OVFL)              , (-FP_OVFL) j. __
0 0 -: qninv  (-FP_OVFL)              , (-FP_OVFL) j.  _
0 0 -: qninv  (-FP_OVFL)              ,   0        j. __
0 0 -: qninv  (-FP_OVFL)              ,   0        j.  _
0 0 -: qninv  (-FP_OVFL)              ,   FP_OVFL  j. __
0 0 -: qninv  (-FP_OVFL)              ,   FP_OVFL  j.  _
0 0 -: qninv  (-FP_OVFL)              ,   _        j. -FP_OVFL
0 0 -: qninv  (-FP_OVFL)              ,   _
0 0 -: qninv  (-FP_OVFL)              ,   _        j.  FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) ,  __        j. -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) ,  __
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) ,  __        j.  FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. __
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  _
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   0        j. __
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  _
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. __
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  _
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   _        j. -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   _
0 0 -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   _        j.  FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  _      ) , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  _      ) ,  -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  _      ) , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  _      ) ,   0        j. -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  _      ) ,   0        j.  FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  _      ) ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  _      ) ,   FP_OVFL
0 0 -: qninv ((-FP_OVFL) j.  _      ) ,   FP_OVFL  j.  FP_OVFL
0 0 -: qninv (  0        j. __      ) , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv (  0        j. __      ) ,  -FP_OVFL
0 0 -: qninv (  0        j. __      ) , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv (  0        j. __      ) ,   0        j. -FP_OVFL
0 0 -: qninv (  0        j. __      ) ,   0        j.  FP_OVFL
0 0 -: qninv (  0        j. __      ) ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv (  0        j. __      ) ,   FP_OVFL
0 0 -: qninv (  0        j. __      ) ,   FP_OVFL  j.  FP_OVFL
0 0 -: qninv (  0        j. -FP_OVFL) ,  __        j. -FP_OVFL
0 0 -: qninv (  0        j. -FP_OVFL) ,  __
0 0 -: qninv (  0        j. -FP_OVFL) ,  __        j.  FP_OVFL
0 0 -: qninv (  0        j. -FP_OVFL) , (-FP_OVFL) j. __
0 0 -: qninv (  0        j. -FP_OVFL) , (-FP_OVFL) j.  _
0 0 -: qninv (  0        j. -FP_OVFL) ,   0        j. __
0 0 -: qninv (  0        j. -FP_OVFL) ,   0        j.  _
0 0 -: qninv (  0        j. -FP_OVFL) ,   FP_OVFL  j. __
0 0 -: qninv (  0        j. -FP_OVFL) ,   FP_OVFL  j.  _
0 0 -: qninv (  0        j. -FP_OVFL) ,   _        j. -FP_OVFL
0 0 -: qninv (  0        j. -FP_OVFL) ,   _
0 0 -: qninv (  0        j. -FP_OVFL) ,   _        j.  FP_OVFL
0 0 -: qninv (  0        j.  FP_OVFL) ,  __        j. -FP_OVFL
0 0 -: qninv (  0        j.  FP_OVFL) ,  __
0 0 -: qninv (  0        j.  FP_OVFL) ,  __        j.  FP_OVFL
0 0 -: qninv (  0        j.  FP_OVFL) , (-FP_OVFL) j. __
0 0 -: qninv (  0        j.  FP_OVFL) , (-FP_OVFL) j.  _
0 0 -: qninv (  0        j.  FP_OVFL) ,   0        j. __
0 0 -: qninv (  0        j.  FP_OVFL) ,   0        j.  _
0 0 -: qninv (  0        j.  FP_OVFL) ,   FP_OVFL  j. __
0 0 -: qninv (  0        j.  FP_OVFL) ,   FP_OVFL  j.  _
0 0 -: qninv (  0        j.  FP_OVFL) ,   _        j. -FP_OVFL
0 0 -: qninv (  0        j.  FP_OVFL) ,   _
0 0 -: qninv (  0        j.  FP_OVFL) ,   _        j.  FP_OVFL
0 0 -: qninv (  0        j.  _      ) , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv (  0        j.  _      ) ,  -FP_OVFL
0 0 -: qninv (  0        j.  _      ) , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv (  0        j.  _      ) ,   0        j. -FP_OVFL
0 0 -: qninv (  0        j.  _      ) ,   0        j.  FP_OVFL
0 0 -: qninv (  0        j.  _      ) ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv (  0        j.  _      ) ,   FP_OVFL
0 0 -: qninv (  0        j.  _      ) ,   FP_OVFL  j.  FP_OVFL
0 0 -: qninv (  FP_OVFL  j. __      ) , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv (  FP_OVFL  j. __      ) ,  -FP_OVFL
0 0 -: qninv (  FP_OVFL  j. __      ) , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv (  FP_OVFL  j. __      ) ,   0        j. -FP_OVFL
0 0 -: qninv (  FP_OVFL  j. __      ) ,   0        j.  FP_OVFL
0 0 -: qninv (  FP_OVFL  j. __      ) ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv (  FP_OVFL  j. __      ) ,   FP_OVFL
0 0 -: qninv (  FP_OVFL  j. __      ) ,   FP_OVFL  j.  FP_OVFL
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) ,  __        j. -FP_OVFL
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) ,  __
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) ,  __        j.  FP_OVFL
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. __
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  _
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) ,   0        j. __
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) ,   0        j.  _
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. __
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  _
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) ,   _        j. -FP_OVFL
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) ,   _
0 0 -: qninv (  FP_OVFL  j. -FP_OVFL) ,   _        j.  FP_OVFL
0 0 -: qninv    FP_OVFL               ,  __        j. -FP_OVFL
0 0 -: qninv    FP_OVFL               ,  __
0 0 -: qninv    FP_OVFL               ,  __        j.  FP_OVFL
0 0 -: qninv    FP_OVFL               , (-FP_OVFL) j. __
0 0 -: qninv    FP_OVFL               , (-FP_OVFL) j.  _
0 0 -: qninv    FP_OVFL               ,   0        j. __
0 0 -: qninv    FP_OVFL               ,   0        j.  _
0 0 -: qninv    FP_OVFL               ,   FP_OVFL  j. __
0 0 -: qninv    FP_OVFL               ,   FP_OVFL  j.  _
0 0 -: qninv    FP_OVFL               ,   _        j. -FP_OVFL
0 0 -: qninv    FP_OVFL               ,   _
0 0 -: qninv    FP_OVFL               ,   _        j.  FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) ,  __        j. -FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) ,  __
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) ,  __        j.  FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. __
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  _
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) ,   0        j. __
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) ,   0        j.  _
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. __
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  _
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) ,   _        j. -FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) ,   _
0 0 -: qninv (  FP_OVFL  j.  FP_OVFL) ,   _        j.  FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  _      ) , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  _      ) ,  -FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  _      ) , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  _      ) ,   0        j. -FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  _      ) ,   0        j.  FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  _      ) ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  _      ) ,   FP_OVFL
0 0 -: qninv (  FP_OVFL  j.  _      ) ,   FP_OVFL  j.  FP_OVFL
0 0 -: qninv (  _        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv (  _        j. -FP_OVFL) ,  -FP_OVFL
0 0 -: qninv (  _        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv (  _        j. -FP_OVFL) ,   0        j. -FP_OVFL
0 0 -: qninv (  _        j. -FP_OVFL) ,   0        j.  FP_OVFL
0 0 -: qninv (  _        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv (  _        j. -FP_OVFL) ,   FP_OVFL
0 0 -: qninv (  _        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
0 0 -: qninv    _                     , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv    _                     ,  -FP_OVFL
0 0 -: qninv    _                     , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv    _                     ,   0        j. -FP_OVFL
0 0 -: qninv    _                     ,   0        j.  FP_OVFL
0 0 -: qninv    _                     ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv    _                     ,   FP_OVFL
0 0 -: qninv    _                     ,   FP_OVFL  j.  FP_OVFL
0 0 -: qninv (  _        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
0 0 -: qninv (  _        j.  FP_OVFL) ,  -FP_OVFL
0 0 -: qninv (  _        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
0 0 -: qninv (  _        j.  FP_OVFL) ,   0        j. -FP_OVFL
0 0 -: qninv (  _        j.  FP_OVFL) ,   0        j.  FP_OVFL
0 0 -: qninv (  _        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
0 0 -: qninv (  _        j.  FP_OVFL) ,   FP_OVFL
0 0 -: qninv (  _        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(2 2 $ 0         ) -: qninv (_ , FP_OVFL) ,: FP_OVFL , __
(0 ,:  2j_2 _2j_2) -: qninv (_ , FP_OVFL) ,: q4
(0 ,:~ 2j_2 _2j_2) -: qninv  q4           ,: FP_OVFL , __
NB. - input is trivial
 0 0                         -: qninv ( __        j. -FP_OVFL) ,   0
 0 0                         -: qninv ( __        j. _1      ) ,   0
 0 0                         -: qninv ( __        j. -FP_UNFL) ,   0
 0 0                         -: qninv   __                         0
 0 0                         -: qninv ( __        j.  FP_UNFL) ,   0
 0 0                         -: qninv ( __        j.  1      ) ,   0
 0 0                         -: qninv ( __        j.  FP_OVFL) ,   0
 0 0                         -: qninv ((-FP_OVFL) j. __      ) ,   0
(0 ,~ (_1j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   0
(0 ,~  _1         % FP_OVFL) -: qninv ((-FP_OVFL) j. _1      ) ,   0
(0 ,~  _1         % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   0
(0 ,~  _1         % FP_OVFL) -: qninv  (-FP_OVFL)              ,   0
(0 ,~  _1         % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   0
(0 ,~  _1         % FP_OVFL) -: qninv ((-FP_OVFL) j.  1      ) ,   0
(0 ,~ (_1j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   0
 0 0                         -: qninv ((-FP_OVFL) j.  _      ) ,   0
 0 0                         -: qninv ((-sqr_big) j. __      ) ,   0
(0 ,~ (_1j1  % 2) % sqr_big) -: qninv ((-sqr_big) j. -sqr_big) ,   0
(0 ,~  _1         % sqr_big) -: qninv  (-sqr_big)              ,   0
(0 ,~ ( 1j1  %_2) % sqr_big) -: qninv ((-sqr_big) j.  sqr_big) ,   0
 0 0                         -: qninv ((-sqr_big) j.  _      ) ,   0
 0 0                         -: qninv ( _1        j. __      ) ,   0
(0 ,~   0j1       % FP_OVFL) -: qninv ( _1        j. -FP_OVFL) ,   0
(0 ,~  _1j1  % 2           ) -: qninv ( _1        j. _1      ) ,   0
(0 ,~  _1 j.  FP_UNFL      ) -: qninv ( _1        j. -FP_UNFL) ,   0
_1 0                         -: qninv   _1                         0
(0 ,~ - 1 j.  FP_UNFL      ) -: qninv ( _1        j.  FP_UNFL) ,   0
(0 ,~  _1j_1 % 2           ) -: qninv ( _1        j.  1      ) ,   0
(0 ,~   0j_1      % FP_OVFL) -: qninv ( _1        j.  FP_OVFL) ,   0
 0 0                         -: qninv ( _1        j.  _      ) ,   0
 0 0                         -: qninv ((-FP_UNFL) j. __      ) ,   0
 0 0                         -: qninv ((-FP_UNFL) j.  _      ) ,   0
(0 ,~   0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   0
(0 ,~   1 j.~ -FP_UNFL     ) -: qninv ((-FP_UNFL) j. _1      ) ,   0
(0 ,~ (_1j1  % 2) % FP_UNFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   0
(0 ,~  _1         % FP_UNFL) -: qninv  (-FP_UNFL)              ,   0
(0 ,~ (_1j_1 % 2) % FP_UNFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   0
(0 ,~ - 1 j.~ FP_UNFL      ) -: qninv ((-FP_UNFL) j.  1      ) ,   0
(0 ,~   0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   0
 0 0                         -: qninv (  0        j. __      ) ,   0
(0 ,~   0j1       % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   0
(0 ,~   0j1       % sqr_big) -: qninv (  0        j. -sqr_big) ,   0
 0j1 0                       -: qninv (  0        j. _1      ) ,   0
(0 ,~   0j1       % FP_UNFL) -: qninv (  0        j. -FP_UNFL) ,   0
 0 0                         -: qninv    0                     ,  __        j. -FP_OVFL
 0 0                         -: qninv    0                     ,  __        j. _1
 0 0                         -: qninv    0                     ,  __        j. -FP_UNFL
 0 0                         -: qninv    0                        __
 0 0                         -: qninv    0                     ,  __        j.  FP_UNFL
 0 0                         -: qninv    0                     ,  __        j.  1
 0 0                         -: qninv    0                     ,  __        j.  FP_OVFL
 0 0                         -: qninv    0                     , (-FP_OVFL) j. __
(0 , (1j1%sqr_big)% sqr_big) -: qninv    0                     , (-FP_OVFL) j. -FP_OVFL
(0 ,              % FP_OVFL) -: qninv    0                     , (-FP_OVFL) j. _1
(0 ,              % FP_OVFL) -: qninv    0                     , (-FP_OVFL) j. -FP_UNFL
(0 ,              % FP_OVFL) -: qninv    0                     ,  -FP_OVFL
(0 ,              % FP_OVFL) -: qninv    0                     , (-FP_OVFL) j.  FP_UNFL
(0 ,              % FP_OVFL) -: qninv    0                     , (-FP_OVFL) j.  1
(0 ,  ( 1j_1 % 2) % FP_OVFL) -: qninv    0                     , (-FP_OVFL) j.  FP_OVFL
 0 0                         -: qninv    0                     , (-FP_OVFL) j.  _
 0 0                         -: qninv    0                     , (-sqr_big) j. __
(0 ,  ( 1j1  % 2) % sqr_big) -: qninv    0                     , (-sqr_big) j. -sqr_big
(0 ,              % sqr_big) -: qninv    0                     ,  -sqr_big
(0 ,  ( 1j_1 % 2) % sqr_big) -: qninv    0                     , (-sqr_big) j.  sqr_big
 0 0                         -: qninv    0                     , (-sqr_big) j.  _
 0 0                         -: qninv    0                     ,  _1        j. __
(0 ,    0j1       % FP_OVFL) -: qninv    0                     ,  _1        j. -FP_OVFL
(0 ,    1j1  % 2           ) -: qninv    0                     ,  _1        j. _1
(0 ,    1 j.  FP_UNFL      ) -: qninv    0                     ,  _1        j. -FP_UNFL
 0 1                         -: qninv    0                        _1
(0 ,    1 j. -FP_UNFL      ) -: qninv    0                     ,  _1        j.  FP_UNFL
(0 ,    1j_1 % 2           ) -: qninv    0                     ,  _1        j.  1
(0 ,    0j_1      % FP_OVFL) -: qninv    0                     ,  _1        j.  FP_OVFL
 0 0                         -: qninv    0                     ,  _1        j.  _
 0 0                         -: qninv    0                     , (-FP_UNFL) j. __
(0 ,    0j1  % FP_OVFL     ) -: qninv    0                     , (-FP_UNFL) j. -FP_OVFL
(0 ,    1 j.~ FP_UNFL      ) -: qninv    0                     , (-FP_UNFL) j. _1
(0 ,  ( 1j1  % 2) % FP_UNFL) -: qninv    0                     , (-FP_UNFL) j. -FP_UNFL
(0 ,              % FP_UNFL) -: qninv    0                     ,  -FP_UNFL
(0 ,  ( 1j_1 % 2) % FP_UNFL) -: qninv    0                     , (-FP_UNFL) j.  FP_UNFL
(0 ,   _1 j.~ FP_UNFL      ) -: qninv    0                     , (-FP_UNFL) j.  1
(0 ,    0j_1      % FP_OVFL) -: qninv    0                     , (-FP_UNFL) j.  FP_OVFL
 0 0                         -: qninv    0                     , (-FP_UNFL) j.  _
 0 0                         -: qninv    0                     ,   0        j. __
(0 ,    0j1       % FP_OVFL) -: qninv    0                     ,   0        j. -FP_OVFL
(0 ,    0 j. % sqr_big     ) -: qninv    0                     ,   0        j. -sqr_big
 0 0j1                       -: qninv    0                     ,   0        j. _1
(0 ,    0 j.      % FP_UNFL) -: qninv    0                     ,   0        j. -FP_UNFL
(0 ,    0 j. -    % FP_UNFL) -: qninv    0                     ,   0        j.  FP_UNFL
 0 0j_1                      -: qninv    0                     ,   0        j.  1
(0 ,    0j_1      % sqr_big) -: qninv    0                     ,   0        j.  sqr_big
(0 ,    0j_1      % FP_OVFL) -: qninv    0                     ,   0        j.  FP_OVFL
 0 0                         -: qninv    0                     ,   0        j.  _
 0 0                         -: qninv    0                     ,   FP_UNFL  j. __
(0 ,    0j1       % FP_OVFL) -: qninv    0                     ,   FP_UNFL  j. -FP_OVFL
(0 ,    1 j.~ FP_UNFL      ) -: qninv    0                     ,   FP_UNFL  j. _1
(0 ,  (_1j1  % 2) % FP_UNFL) -: qninv    0                     ,   FP_UNFL  j. -FP_UNFL
(0 ,  -           % FP_UNFL) -: qninv    0                     ,   FP_UNFL
(0 ,  ( 1j1  %_2) % FP_UNFL) -: qninv    0                     ,   FP_UNFL  j.  FP_UNFL
(0 ,  - 1 j.~ FP_UNFL      ) -: qninv    0                     ,   FP_UNFL  j.  1
(0 ,    0j_1      % FP_OVFL) -: qninv    0                     ,   FP_UNFL  j.  FP_OVFL
 0 0                         -: qninv    0                     ,   FP_UNFL  j.  _
 0 0                         -: qninv    0                     ,   1        j. __
(0 ,    0j1       % FP_OVFL) -: qninv    0                     ,   1        j. -FP_OVFL
(0 ,   _1j1  % 2           ) -: qninv    0                     ,   1        j. _1
(0 ,   _1 j.  FP_UNFL      ) -: qninv    0                     ,   1        j. -FP_UNFL
 0 _1                        -: qninv    0                         1
(0 ,  - 1 j.  FP_UNFL      ) -: qninv    0                     ,   1        j.  FP_UNFL
(0 ,    1j1  %_2           ) -: qninv    0                     ,   1        j.  1
(0 ,    0j_1      % FP_OVFL) -: qninv    0                     ,   1        j.  FP_OVFL
 0 0                         -: qninv    0                     ,   1        j.  _
 0 0                         -: qninv    0                     ,   sqr_big  j. __
(0 ,  (_1j1  % 2) % sqr_big) -: qninv    0                     ,   sqr_big  j. -sqr_big
(0 ,  -           % sqr_big) -: qninv    0                     ,   sqr_big
(0 ,  (_1j_1 % 2) % sqr_big) -: qninv    0                     ,   sqr_big  j.  sqr_big
 0 0                         -: qninv    0                     ,   sqr_big  j.  _
 0 0                         -: qninv    0                     ,   FP_OVFL  j. __
(0 ,  (_1j1  % 2) % FP_OVFL) -: qninv    0                     ,   FP_OVFL  j. -FP_OVFL
(0 ,  -           % FP_OVFL) -: qninv    0                     ,   FP_OVFL  j. _1
(0 ,  -           % FP_OVFL) -: qninv    0                     ,   FP_OVFL  j. -FP_UNFL
(0 ,  -           % FP_OVFL) -: qninv    0                     ,   FP_OVFL
(0 ,  -           % FP_OVFL) -: qninv    0                     ,   FP_OVFL  j.  FP_UNFL
(0 ,  -           % FP_OVFL) -: qninv    0                     ,   FP_OVFL  j.  1
(0 ,  ( 1j1  %_2) % FP_OVFL) -: qninv    0                     ,   FP_OVFL  j.  FP_OVFL
 0 0                         -: qninv    0                     ,   FP_OVFL  j.  _
 0 0                         -: qninv    0                     ,   _        j. -FP_OVFL
 0 0                         -: qninv    0                     ,   _        j. _1
 0 0                         -: qninv    0                     ,   _        j. -FP_UNFL
 0 0                         -: qninv    0                         _
 0 0                         -: qninv    0                     ,   _        j.  FP_UNFL
 0 0                         -: qninv    0                     ,   _        j.  1
 0 0                         -: qninv    0                     ,   _        j.  FP_OVFL
(0 ,~   0j_1      % FP_UNFL) -: qninv (  0        j.  FP_UNFL) ,   0
 0j_1 0                      -: qninv (  0        j.  1      ) ,   0
(0 ,~   0j_1      % sqr_big) -: qninv (  0        j.  sqr_big) ,   0
(0 ,~   0j_1      % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   0
 0 0                         -: qninv (  0        j.  _      ) ,   0
 0 0                         -: qninv (  FP_UNFL  j. __      ) ,   0
(0 ,~   0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   0
(0 ,~   1 j.~ FP_UNFL      ) -: qninv (  FP_UNFL  j. _1      ) ,   0
(0 ,~ ( 1j1  % 2) % FP_UNFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   0
(0 ,~             % FP_UNFL) -: qninv    FP_UNFL               ,   0
(0 ,~ ( 1j_1 % 2) % FP_UNFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   0
(0 ,~  _1 j.~ FP_UNFL      ) -: qninv (  FP_UNFL  j.  1      ) ,   0
(0 ,~   0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   0
 0 0                         -: qninv (  FP_UNFL  j.  _      ) ,   0
 0 0                         -: qninv (  1        j. __      ) ,   0
(0 ,~   0j1       % FP_OVFL) -: qninv (  1        j. -FP_OVFL) ,   0
(0 ,~   1j1  % 2           ) -: qninv (  1        j. _1      ) ,   0
(0 ,~   1 j.  FP_UNFL      ) -: qninv (  1        j. -FP_UNFL) ,   0
1 0                          -: qninv    1                         0
(0 ,~   1 j. -FP_UNFL      ) -: qninv (  1        j.  FP_UNFL) ,   0
(0 ,~   1j_1 % 2           ) -: qninv (  1        j.  1      ) ,   0
(0 ,~   0j_1      % FP_OVFL) -: qninv (  1        j.  FP_OVFL) ,   0
 0 0                         -: qninv (  1        j.  _      ) ,   0
 0 0                         -: qninv (  sqr_big  j. __      ) ,   0
(0 ,~ ( 1j1  % 2) % sqr_big) -: qninv (  sqr_big  j. -sqr_big) ,   0
(0 ,~             % sqr_big) -: qninv    sqr_big               ,   0
(0 ,~ ( 1j_1 % 2) % sqr_big) -: qninv (  sqr_big  j.  sqr_big) ,   0
 0 0                         -: qninv (  sqr_big  j.  _      ) ,   0
 0 0                         -: qninv (  FP_OVFL  j. __      ) ,   0
(0 ,~ ( 1j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   0
(0 ,~             % FP_OVFL) -: qninv (  FP_OVFL  j. _1      ) ,   0
(0 ,~             % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   0
(0 ,~             % FP_OVFL) -: qninv    FP_OVFL               ,   0
(0 ,~             % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   0
(0 ,~             % FP_OVFL) -: qninv (  FP_OVFL  j.  1      ) ,   0
(0 ,~ ( 1j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   0
 0 0                         -: qninv (  FP_OVFL  j.  _      ) ,   0
 0 0                         -: qninv (  _        j. _1      ) ,   0
 0 0                         -: qninv (  _        j. -FP_OVFL) ,   0
 0 0                         -: qninv (  _        j. -FP_UNFL) ,   0
 0 0                         -: qninv    _                         0
 0 0                         -: qninv (  _        j.  FP_UNFL) ,   0
 0 0                         -: qninv (  _        j.  1      ) ,   0
 0 0                         -: qninv (  _        j.  FP_OVFL) ,   0
NB. - edge cases input
((_1j1   1j1  % 4) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
((_1j1   1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
((_1j1   1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,  -FP_OVFL
((_1j1   1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
((_1j1   1j_1 % 4) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
((_1     1j1  % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
((_1     1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
((_1     1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,  -FP_OVFL
((_1     1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
((_1     1j_1 % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
((_1     1j1  % 3) % FP_OVFL) -: qninv  (-FP_OVFL)              , (-FP_OVFL) j. -FP_OVFL
((_1     1    % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              , (-FP_OVFL) j. -FP_UNFL
((_1     1    % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              ,  -FP_OVFL
((_1     1    % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              , (-FP_OVFL) j.  FP_UNFL
((_1     1j_1 % 3) % FP_OVFL) -: qninv  (-FP_OVFL)              , (-FP_OVFL) j.  FP_OVFL
((_1     1j1  % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
((_1     1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
((_1     1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,  -FP_OVFL
((_1     1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
((_1     1j_1 % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
((_1j_1  1j1  % 4) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
((_1j_1  1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
((_1j_1  1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,  -FP_OVFL
((_1j_1  1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
((_1j_1  1j_1 % 4) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
((_1j1   0j1  % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~  (_1j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~  (_1j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,  -FP_UNFL
(0 ,~  (_1j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
((_1j1   0j_1 % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
((_1     0j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,  -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
((_1     0j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
((_1     0j1  % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              , (-FP_UNFL) j. -FP_OVFL
(0 ,~   _1         % FP_OVFL) -: qninv  (-FP_OVFL)              , (-FP_UNFL) j. -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv  (-FP_OVFL)              ,  -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv  (-FP_OVFL)              , (-FP_UNFL) j.  FP_UNFL
((_1     0j_1 % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              , (-FP_UNFL) j.  FP_OVFL
((_1     0j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,  -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
((_1     0j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
((_1j_1  0j1  % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~  (_1j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~  (_1j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,  -FP_UNFL
(0 ,~  (_1j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
((_1j_1  0j_1 % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
((_1j1   0j1  % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   0        j. -FP_OVFL
(0 ,~  (_1j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   0        j. -FP_UNFL
(0 ,~  (_1j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  FP_UNFL
((_1j1   0j_1 % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  FP_OVFL
((_1     0j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   0        j. -FP_OVFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   0        j. -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   0        j.  FP_UNFL
((_1     0j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   0        j.  FP_OVFL
((_1     0j1  % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              ,   0        j. -FP_OVFL
(0 ,~   _1         % FP_OVFL) -: qninv  (-FP_OVFL)              ,   0        j. -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv  (-FP_OVFL)              ,   0        j.  FP_UNFL
((_1     0j_1 % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              ,   0        j.  FP_OVFL
((_1     0j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   0        j. -FP_OVFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   0        j. -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   0        j.  FP_UNFL
((_1     0j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   0        j.  FP_OVFL
((_1j_1  0j1  % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   0        j. -FP_OVFL
(0 ,~  (_1j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   0        j. -FP_UNFL
(0 ,~  (_1j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  FP_UNFL
((_1j_1  0j_1 % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  FP_OVFL
((_1j1   0j1  % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~  (_1j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~  (_1j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL
(0 ,~  (_1j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
((_1j1   0j_1 % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
((_1     0j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
((_1     0j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
((_1     0j1  % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              ,   FP_UNFL  j. -FP_OVFL
(0 ,~   _1         % FP_OVFL) -: qninv  (-FP_OVFL)              ,   FP_UNFL  j. -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv  (-FP_OVFL)              ,   FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv  (-FP_OVFL)              ,   FP_UNFL  j.  FP_UNFL
((_1     0j_1 % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              ,   FP_UNFL  j.  FP_OVFL
((_1     0j1  % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL
(0 ,~   _1         % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
((_1     0j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
((_1j_1  0j1  % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~  (_1j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~  (_1j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL
(0 ,~  (_1j_1 % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
((_1j_1  0j_1 % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
((_1j1  _1j1  % 4) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
((_1j1  _1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
((_1j1  _1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL
((_1j1  _1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
((_1j1  _1j_1 % 4) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
((_1    _1j1  % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
((_1    _1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
((_1    _1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL
((_1    _1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
((_1    _1j_1 % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
((_1    _1j1  % 3) % FP_OVFL) -: qninv  (-FP_OVFL)              ,   FP_OVFL  j. -FP_OVFL
((_1    _1    % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              ,   FP_OVFL  j. -FP_UNFL
((_1    _1    % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              ,   FP_OVFL
((_1    _1    % 2) % FP_OVFL) -: qninv  (-FP_OVFL)              ,   FP_OVFL  j.  FP_UNFL
((_1    _1j_1 % 3) % FP_OVFL) -: qninv  (-FP_OVFL)              ,   FP_OVFL  j.  FP_OVFL
((_1    _1j1  % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
((_1    _1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
((_1    _1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL
((_1    _1    % 2) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
((_1    _1j_1 % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
((_1j_1 _1j1  % 4) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
((_1j_1 _1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
((_1j_1 _1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL
((_1j_1 _1    % 3) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
((_1j_1 _1j_1 % 4) % FP_OVFL) -: qninv ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
((_1j1   1j1  % 4) % sqr_big) -: qninv ((-sqr_big) j. -sqr_big) , (-sqr_big) j. -sqr_big
((_1j1   1j_1 % 4) % sqr_big) -: qninv ((-sqr_big) j. -sqr_big) , (-sqr_big) j.  sqr_big
((_1j1  _1j1  % 4) % sqr_big) -: qninv ((-sqr_big) j. -sqr_big) ,   sqr_big  j. -sqr_big
((_1j1  _1j_1 % 4) % sqr_big) -: qninv ((-sqr_big) j. -sqr_big) ,   sqr_big  j.  sqr_big
((_1     1    % 2) % sqr_big) -: qninv  (-sqr_big)              ,  -sqr_big
((_1     0j1  % 2) % sqr_big) -: qninv  (-sqr_big)              ,   0        j. -sqr_big
((_1     0j_1 % 2) % sqr_big) -: qninv  (-sqr_big)              ,   0        j.  sqr_big
((_1    _1    % 2) % sqr_big) -: qninv  (-sqr_big)              ,   sqr_big
((_1j_1  1j1  % 4) % sqr_big) -: qninv ((-sqr_big) j.  sqr_big) , (-sqr_big) j. -sqr_big
((_1j_1  1j_1 % 4) % sqr_big) -: qninv ((-sqr_big) j.  sqr_big) , (-sqr_big) j.  sqr_big
((_1j_1 _1j1  % 4) % sqr_big) -: qninv ((-sqr_big) j.  sqr_big) ,   sqr_big  j. -sqr_big
((_1j_1 _1j_1 % 4) % sqr_big) -: qninv ((-sqr_big) j.  sqr_big) ,   sqr_big  j.  sqr_big
(( 0j1   1j1  % 3) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(( 0j1   1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(( 0j1   1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,  -FP_OVFL
(( 0j1   1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(( 0j1   1j_1 % 3) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(0 ,   ( 1j1  % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(0 ,               % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(0 ,               % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,  -FP_OVFL
(0 ,               % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(0 ,   ( 1j_1 % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(0 ,   ( 1j1  % 2) % FP_OVFL) -: qninv  (-FP_UNFL)              , (-FP_OVFL) j. -FP_OVFL
(0 ,               % FP_OVFL) -: qninv  (-FP_UNFL)              , (-FP_OVFL) j. -FP_UNFL
(0 ,               % FP_OVFL) -: qninv  (-FP_UNFL)              ,  -FP_OVFL
(0 ,               % FP_OVFL) -: qninv  (-FP_UNFL)              , (-FP_OVFL) j.  FP_UNFL
(0 ,   ( 1j_1 % 2) % FP_OVFL) -: qninv  (-FP_UNFL)              , (-FP_OVFL) j.  FP_OVFL
(0 ,   ( 1j1  % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(0 ,               % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(0 ,               % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,  -FP_OVFL
(0 ,               % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(0 ,   ( 1j_1 % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(( 0j_1  1j1  % 3) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(( 0j_1  1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(( 0j_1  1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,  -FP_OVFL
(( 0j_1  1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(( 0j_1  1j_1 % 3) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(( 0j1   0j1  % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~    0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,  -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(( 0j1   0j_1 % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
((_1j1   1j1  % 4) % FP_UNFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
((_1j1   1    % 3) % FP_UNFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,  -FP_UNFL
((_1j1   1j_1 % 4) % FP_UNFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv  (-FP_UNFL)              , (-FP_UNFL) j. -FP_OVFL
((_1     1j1  % 3) % FP_UNFL) -: qninv  (-FP_UNFL)              , (-FP_UNFL) j. -FP_UNFL
((_1     1    % 2) % FP_UNFL) -: qninv  (-FP_UNFL)              ,  -FP_UNFL
((_1     1j_1 % 3) % FP_UNFL) -: qninv  (-FP_UNFL)              , (-FP_UNFL) j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv  (-FP_UNFL)              , (-FP_UNFL) j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
((_1j_1  1j1  % 4) % FP_UNFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
((_1j_1  1    % 3) % FP_UNFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,  -FP_UNFL
((_1j_1  1j_1 % 4) % FP_UNFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(( 0j_1  0j1  % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~    0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,  -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(( 0j_1  0j_1 % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(( 0j1   0j1  % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   0        j. -FP_OVFL
(0 ,~    0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   0        j. -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   0        j.  FP_UNFL
(( 0j1   0j_1 % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   0        j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   0        j. -FP_OVFL
((_1j1   0j1  % 3) % FP_UNFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   0        j. -FP_UNFL
((_1j1   0j_1 % 3) % FP_UNFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   0        j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   0        j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv  (-FP_UNFL)              ,   0        j. -FP_OVFL
((_1     0j1  % 2) % FP_UNFL) -: qninv  (-FP_UNFL)              ,   0        j. -FP_UNFL
((_1     0j_1 % 2) % FP_UNFL) -: qninv  (-FP_UNFL)              ,   0        j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv  (-FP_UNFL)              ,   0        j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   0        j. -FP_OVFL
((_1j_1  0j1  % 3) % FP_UNFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   0        j. -FP_UNFL
((_1j_1  0j_1 % 3) % FP_UNFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   0        j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   0        j.  FP_OVFL
(( 0j_1  0j1  % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   0        j. -FP_OVFL
(0 ,~    0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   0        j. -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   0        j.  FP_UNFL
(( 0j_1  0j_1 % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   0        j.  FP_OVFL
(( 0j1   0j1  % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~    0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(( 0j1   0j_1 % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
((_1j1  _1j1  % 4) % FP_UNFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
((_1j1  _1    % 3) % FP_UNFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL
((_1j1  _1j_1 % 4) % FP_UNFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv  (-FP_UNFL)              ,   FP_UNFL  j. -FP_OVFL
((_1    _1j1  % 3) % FP_UNFL) -: qninv  (-FP_UNFL)              ,   FP_UNFL  j. -FP_UNFL
((_1    _1    % 2) % FP_UNFL) -: qninv  (-FP_UNFL)              ,   FP_UNFL
((_1    _1j_1 % 3) % FP_UNFL) -: qninv  (-FP_UNFL)              ,   FP_UNFL  j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv  (-FP_UNFL)              ,   FP_UNFL  j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
((_1j_1 _1j1  % 4) % FP_UNFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
((_1j_1 _1    % 3) % FP_UNFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL
((_1j_1 _1j_1 % 4) % FP_UNFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(( 0j_1  0j1  % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~    0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(( 0j_1  0j_1 % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(( 0j1  _1j1  % 3) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(( 0j1  _1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(( 0j1  _1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL
(( 0j1  _1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(( 0j1  _1j_1 % 3) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(0 ,   (_1j1  % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(0 ,    _1         % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(0 ,   (_1j_1 % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(0 ,   (_1j1  % 2) % FP_OVFL) -: qninv  (-FP_UNFL)              ,   FP_OVFL  j. -FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv  (-FP_UNFL)              ,   FP_OVFL  j. -FP_UNFL
(0 ,    _1         % FP_OVFL) -: qninv  (-FP_UNFL)              ,   FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv  (-FP_UNFL)              ,   FP_OVFL  j.  FP_UNFL
(0 ,   (_1j_1 % 2) % FP_OVFL) -: qninv  (-FP_UNFL)              ,   FP_OVFL  j.  FP_OVFL
(0 ,   (_1j1  % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(0 ,    _1         % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(0 ,   (_1j_1 % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(( 0j_1 _1j1  % 3) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(( 0j_1 _1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(( 0j_1 _1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL
(( 0j_1 _1    % 2) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(( 0j_1 _1j_1 % 3) % FP_OVFL) -: qninv ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(( 0j1   1j1  % 3) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(( 0j1   1    % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(( 0j1   1    % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,  -FP_OVFL
(( 0j1   1    % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(( 0j1   1j_1 % 3) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(( 0j1   1    % 2) % sqr_big) -: qninv (  0        j. -sqr_big) ,  -sqr_big
(( 0j1   0j1  % 2) % sqr_big) -: qninv (  0        j. -sqr_big) ,   0        j. -sqr_big
(( 0j1   0j_1 % 2) % sqr_big) -: qninv (  0        j. -sqr_big) ,   0        j.  sqr_big
(( 0j1  _1    % 2) % sqr_big) -: qninv (  0        j. -sqr_big) ,   sqr_big
(0 ,   ( 1j1  % 2) % FP_OVFL) -: qninv (  0        j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(0 ,               % FP_OVFL) -: qninv (  0        j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(0 ,               % FP_OVFL) -: qninv (  0        j. -FP_UNFL) ,  -FP_OVFL
(0 ,               % FP_OVFL) -: qninv (  0        j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(0 ,   ( 1j_1 % 2) % FP_OVFL) -: qninv (  0        j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(0 ,   ( 1j1  % 2) % FP_OVFL) -: qninv (  0        j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(0 ,               % FP_OVFL) -: qninv (  0        j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(0 ,               % FP_OVFL) -: qninv (  0        j.  FP_UNFL) ,  -FP_OVFL
(0 ,               % FP_OVFL) -: qninv (  0        j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(0 ,   ( 1j_1 % 2) % FP_OVFL) -: qninv (  0        j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(( 0j_1  1j1  % 3) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(( 0j_1  1    % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(( 0j_1  1    % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,  -FP_OVFL
(( 0j_1  1    % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(( 0j_1  1j_1 % 3) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(( 0j1   0j1  % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  0        j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,  -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  0        j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(( 0j1   0j_1 % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  0        j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(( 0j1   1j1  % 3) % FP_UNFL) -: qninv (  0        j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(( 0j1   1    % 2) % FP_UNFL) -: qninv (  0        j. -FP_UNFL) ,  -FP_UNFL
(( 0j1   1j_1 % 3) % FP_UNFL) -: qninv (  0        j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  0        j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  0        j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(( 0j_1  1j1  % 3) % FP_UNFL) -: qninv (  0        j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(( 0j_1  1    % 2) % FP_UNFL) -: qninv (  0        j.  FP_UNFL) ,  -FP_UNFL
(( 0j_1  1j_1 % 3) % FP_UNFL) -: qninv (  0        j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  0        j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(( 0j_1  1    % 2) % sqr_big) -: qninv (  0        j.  sqr_big) ,  -sqr_big
(( 0j_1  0j1  % 2) % sqr_big) -: qninv (  0        j.  sqr_big) ,   0        j. -sqr_big
(( 0j_1  0j_1 % 2) % sqr_big) -: qninv (  0        j.  sqr_big) ,   0        j.  sqr_big
(( 0j_1 _1    % 2) % sqr_big) -: qninv (  0        j.  sqr_big) ,   sqr_big
(( 0j_1  0j1  % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  0        j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,  -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  0        j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(( 0j_1  0j_1 % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(( 0j1   0j1  % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   0        j. -FP_OVFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   0        j. -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   0        j.  FP_UNFL
(( 0j1   0j_1 % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   0        j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  0        j. -FP_UNFL) ,   0        j. -FP_OVFL
(( 0j1   0j1  % 2) % FP_UNFL) -: qninv (  0        j. -FP_UNFL) ,   0        j. -FP_UNFL
(( 0j1   0j_1 % 2) % FP_UNFL) -: qninv (  0        j. -FP_UNFL) ,   0        j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  0        j. -FP_UNFL) ,   0        j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  0        j.  FP_UNFL) ,   0        j. -FP_OVFL
(( 0j_1  0j1  % 2) % FP_UNFL) -: qninv (  0        j.  FP_UNFL) ,   0        j. -FP_UNFL
(( 0j_1  0j_1 % 2) % FP_UNFL) -: qninv (  0        j.  FP_UNFL) ,   0        j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  0        j.  FP_UNFL) ,   0        j.  FP_OVFL
(( 0j_1  0j1  % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   0        j. -FP_OVFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   0        j. -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   0        j.  FP_UNFL
(( 0j_1  0j_1 % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   0        j.  FP_OVFL
(( 0j1   0j1  % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(( 0j1   0j_1 % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  0        j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(( 0j1  _1j1  % 3) % FP_UNFL) -: qninv (  0        j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(( 0j1  _1    % 2) % FP_UNFL) -: qninv (  0        j. -FP_UNFL) ,   FP_UNFL
(( 0j1  _1j_1 % 3) % FP_UNFL) -: qninv (  0        j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  0        j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  0        j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(( 0j_1 _1j1  % 3) % FP_UNFL) -: qninv (  0        j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(( 0j_1 _1    % 2) % FP_UNFL) -: qninv (  0        j.  FP_UNFL) ,   FP_UNFL
(( 0j_1 _1j_1 % 3) % FP_UNFL) -: qninv (  0        j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  0        j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(( 0j_1  0j1  % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(( 0j_1  0j_1 % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(( 0j1  _1j1  % 3) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(( 0j1  _1    % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(( 0j1  _1    % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   FP_OVFL
(( 0j1  _1    % 2) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(( 0j1  _1j_1 % 3) % FP_OVFL) -: qninv (  0        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(0 ,   (_1j1  % 2) % FP_OVFL) -: qninv (  0        j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv (  0        j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(0 ,    _1         % FP_OVFL) -: qninv (  0        j. -FP_UNFL) ,   FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv (  0        j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(0 ,   (_1j_1 % 2) % FP_OVFL) -: qninv (  0        j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(0 ,   (_1j1  % 2) % FP_OVFL) -: qninv (  0        j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv (  0        j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(0 ,    _1         % FP_OVFL) -: qninv (  0        j.  FP_UNFL) ,   FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv (  0        j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(0 ,   (_1j_1 % 2) % FP_OVFL) -: qninv (  0        j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(( 0j_1 _1j1  % 3) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(( 0j_1 _1    % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(( 0j_1 _1    % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   FP_OVFL
(( 0j_1 _1    % 2) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(( 0j_1 _1j_1 % 3) % FP_OVFL) -: qninv (  0        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(( 0j1   1j1  % 3) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(( 0j1   1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(( 0j1   1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,  -FP_OVFL
(( 0j1   1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(( 0j1   1j_1 % 3) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(0 ,   ( 1j1  % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(0 ,               % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(0 ,               % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,  -FP_OVFL
(0 ,               % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(0 ,   ( 1j_1 % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(0 ,   ( 1j1  % 2) % FP_OVFL) -: qninv    FP_UNFL               , (-FP_OVFL) j. -FP_OVFL
(0 ,               % FP_OVFL) -: qninv    FP_UNFL               , (-FP_OVFL) j. -FP_UNFL
(0 ,               % FP_OVFL) -: qninv    FP_UNFL               ,  -FP_OVFL
(0 ,               % FP_OVFL) -: qninv    FP_UNFL               , (-FP_OVFL) j.  FP_UNFL
(0 ,   ( 1j_1 % 2) % FP_OVFL) -: qninv    FP_UNFL               , (-FP_OVFL) j.  FP_OVFL
(0 ,   ( 1j1  % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(0 ,               % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(0 ,               % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,  -FP_OVFL
(0 ,               % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(0 ,   ( 1j_1 % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(( 0j_1  1j1  % 3) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(( 0j_1  1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(( 0j_1  1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,  -FP_OVFL
(( 0j_1  1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(( 0j_1  1j_1 % 3) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(( 0j1   0j1  % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,  -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(( 0j1   0j_1 % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(( 1j1   1j1  % 4) % FP_UNFL) -: qninv (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(( 1j1   1    % 3) % FP_UNFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,  -FP_UNFL
(( 1j1   1j_1 % 4) % FP_UNFL) -: qninv (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv    FP_UNFL               , (-FP_UNFL) j. -FP_OVFL
(( 1     1j1  % 3) % FP_UNFL) -: qninv    FP_UNFL               , (-FP_UNFL) j. -FP_UNFL
(( 1     1    % 2) % FP_UNFL) -: qninv    FP_UNFL               ,  -FP_UNFL
(( 1     1j_1 % 3) % FP_UNFL) -: qninv    FP_UNFL               , (-FP_UNFL) j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv    FP_UNFL               , (-FP_UNFL) j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(( 1j_1  1j1  % 4) % FP_UNFL) -: qninv (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(( 1j_1  1    % 3) % FP_UNFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,  -FP_UNFL
(( 1j_1  1j_1 % 4) % FP_UNFL) -: qninv (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(( 0j_1  0j1  % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,  -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(( 0j_1  0j_1 % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(( 0j1   0j1  % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   0        j. -FP_OVFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   0        j. -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   0        j.  FP_UNFL
(( 0j1   0j_1 % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   0        j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   0        j. -FP_OVFL
(( 1j1   0j1  % 3) % FP_UNFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   0        j. -FP_UNFL
(( 1j1   0j_1 % 3) % FP_UNFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   0        j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   0        j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv    FP_UNFL               ,   0        j. -FP_OVFL
(( 1     0j1  % 2) % FP_UNFL) -: qninv    FP_UNFL               ,   0        j. -FP_UNFL
(( 1     0j_1 % 2) % FP_UNFL) -: qninv    FP_UNFL               ,   0        j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv    FP_UNFL               ,   0        j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   0        j. -FP_OVFL
(( 1j_1  0j1  % 3) % FP_UNFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   0        j. -FP_UNFL
(( 1j_1  0j_1 % 3) % FP_UNFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   0        j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   0        j.  FP_OVFL
(( 0j_1  0j1  % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   0        j. -FP_OVFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   0        j. -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   0        j.  FP_UNFL
(( 0j_1  0j_1 % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   0        j.  FP_OVFL
(( 0j1   0j1  % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL
(0 ,~    0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(( 0j1   0j_1 % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(( 1j1  _1j1  % 4) % FP_UNFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(( 1j1  _1    % 3) % FP_UNFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL
(( 1j1  _1j_1 % 4) % FP_UNFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv    FP_UNFL               ,   FP_UNFL  j. -FP_OVFL
(( 1    _1j1  % 3) % FP_UNFL) -: qninv    FP_UNFL               ,   FP_UNFL  j. -FP_UNFL
(( 1    _1    % 2) % FP_UNFL) -: qninv    FP_UNFL               ,   FP_UNFL
(( 1    _1j_1 % 3) % FP_UNFL) -: qninv    FP_UNFL               ,   FP_UNFL  j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv    FP_UNFL               ,   FP_UNFL  j.  FP_OVFL
(0 ,     0j1       % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(( 1j_1 _1j1  % 4) % FP_UNFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(( 1j_1 _1    % 3) % FP_UNFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL
(( 1j_1 _1j_1 % 4) % FP_UNFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(0 ,     0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(( 0j_1  0j1  % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL
(0 ,~    0j_1      % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(( 0j_1  0j_1 % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(( 0j1  _1j1  % 3) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(( 0j1  _1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(( 0j1  _1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL
(( 0j1  _1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(( 0j1  _1j_1 % 3) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(0 ,   (_1j1  % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(0 ,    _1         % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(0 ,   (_1j_1 % 2) % FP_OVFL) -: qninv (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(0 ,   (_1j1  % 2) % FP_OVFL) -: qninv    FP_UNFL               ,   FP_OVFL  j. -FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv    FP_UNFL               ,   FP_OVFL  j. -FP_UNFL
(0 ,    _1         % FP_OVFL) -: qninv    FP_UNFL               ,   FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv    FP_UNFL               ,   FP_OVFL  j.  FP_UNFL
(0 ,   (_1j_1 % 2) % FP_OVFL) -: qninv    FP_UNFL               ,   FP_OVFL  j.  FP_OVFL
(0 ,   (_1j1  % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(0 ,    _1         % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL
(0 ,    _1         % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(0 ,   (_1j_1 % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(( 0j_1 _1j1  % 3) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(( 0j_1 _1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(( 0j_1 _1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL
(( 0j_1 _1    % 2) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(( 0j_1 _1j_1 % 3) % FP_OVFL) -: qninv (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(( 1j1   1j1  % 4) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(( 1j1   1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(( 1j1   1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,  -FP_OVFL
(( 1j1   1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(( 1j1   1j_1 % 4) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(( 1     1j1  % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(( 1     1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(( 1     1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,  -FP_OVFL
(( 1     1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(( 1     1j_1 % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(( 1     1j1  % 3) % FP_OVFL) -: qninv    FP_OVFL               , (-FP_OVFL) j. -FP_OVFL
(( 1     1    % 2) % FP_OVFL) -: qninv    FP_OVFL               , (-FP_OVFL) j. -FP_UNFL
(( 1     1    % 2) % FP_OVFL) -: qninv    FP_OVFL               ,  -FP_OVFL
(( 1     1    % 2) % FP_OVFL) -: qninv    FP_OVFL               , (-FP_OVFL) j.  FP_UNFL
(( 1     1j_1 % 3) % FP_OVFL) -: qninv    FP_OVFL               , (-FP_OVFL) j.  FP_OVFL
(( 1     1j1  % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(( 1     1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(( 1     1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,  -FP_OVFL
(( 1     1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(( 1     1j_1 % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(( 1j_1  1j1  % 4) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(( 1j_1  1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(( 1j_1  1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,  -FP_OVFL
(( 1j_1  1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(( 1j_1  1j_1 % 4) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(( 1j1   0j1  % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~  ( 1j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~  ( 1j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,  -FP_UNFL
(0 ,~  ( 1j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(( 1j1   0j_1 % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(( 1     0j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,  -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(( 1     0j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(( 1     0j1  % 2) % FP_OVFL) -: qninv    FP_OVFL               , (-FP_UNFL) j. -FP_OVFL
(0 ,~              % FP_OVFL) -: qninv    FP_OVFL               , (-FP_UNFL) j. -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv    FP_OVFL               ,  -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv    FP_OVFL               , (-FP_UNFL) j.  FP_UNFL
(( 1     0j_1 % 2) % FP_OVFL) -: qninv    FP_OVFL               , (-FP_UNFL) j.  FP_OVFL
(( 1     0j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,  -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(( 1     0j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(( 1j_1  0j1  % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
(0 ,~  ( 1j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
(0 ,~  ( 1j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,  -FP_UNFL
(0 ,~  ( 1j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(( 1j_1  0j_1 % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(( 1j1   0j1  % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   0        j. -FP_OVFL
(0 ,~  ( 1j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   0        j. -FP_UNFL
(0 ,~  ( 1j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   0        j.  FP_UNFL
(( 1j1   0j_1 % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   0        j.  FP_OVFL
(( 1     0j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   0        j. -FP_OVFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   0        j. -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   0        j.  FP_UNFL
(( 1     0j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   0        j.  FP_OVFL
(( 1     0j1  % 2) % FP_OVFL) -: qninv    FP_OVFL               ,   0        j. -FP_OVFL
(0 ,~              % FP_OVFL) -: qninv    FP_OVFL               ,   0        j. -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv    FP_OVFL               ,   0        j.  FP_UNFL
(( 1     0j_1 % 2) % FP_OVFL) -: qninv    FP_OVFL               ,   0        j.  FP_OVFL
(( 1     0j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   0        j. -FP_OVFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   0        j. -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   0        j.  FP_UNFL
(( 1     0j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   0        j.  FP_OVFL
(( 1j_1  0j1  % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   0        j. -FP_OVFL
(0 ,~  ( 1j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   0        j. -FP_UNFL
(0 ,~  ( 1j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   0        j.  FP_UNFL
(( 1j_1  0j_1 % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   0        j.  FP_OVFL
(( 1j1   0j1  % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~  ( 1j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~  ( 1j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL
(0 ,~  ( 1j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(( 1j1   0j_1 % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(( 1     0j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(( 1     0j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(( 1     0j1  % 2) % FP_OVFL) -: qninv    FP_OVFL               ,   FP_UNFL  j. -FP_OVFL
(0 ,~              % FP_OVFL) -: qninv    FP_OVFL               ,   FP_UNFL  j. -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv    FP_OVFL               ,   FP_UNFL
(0 ,~              % FP_OVFL) -: qninv    FP_OVFL               ,   FP_UNFL  j.  FP_UNFL
(( 1     0j_1 % 2) % FP_OVFL) -: qninv    FP_OVFL               ,   FP_UNFL  j.  FP_OVFL
(( 1     0j1  % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL
(0 ,~              % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(( 1     0j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(( 1j1   1j1  % 4) % sqr_big) -: qninv (  sqr_big  j. -sqr_big) , (-sqr_big) j. -sqr_big
(( 1j1   1j_1 % 4) % sqr_big) -: qninv (  sqr_big  j. -sqr_big) , (-sqr_big) j.  sqr_big
(( 1j1  _1j1  % 4) % sqr_big) -: qninv (  sqr_big  j. -sqr_big) ,   sqr_big  j. -sqr_big
(( 1j1  _1j_1 % 4) % sqr_big) -: qninv (  sqr_big  j. -sqr_big) ,   sqr_big  j.  sqr_big
(( 1     1    % 2) % sqr_big) -: qninv    sqr_big               ,  -sqr_big
(( 1     0j1  % 2) % sqr_big) -: qninv    sqr_big               ,   0        j. -sqr_big
(( 1     0j_1 % 2) % sqr_big) -: qninv    sqr_big               ,   0        j.  sqr_big
(( 1    _1    % 2) % sqr_big) -: qninv    sqr_big               ,   sqr_big
(( 1j_1  1j1  % 4) % sqr_big) -: qninv (  sqr_big  j.  sqr_big) , (-sqr_big) j. -sqr_big
(( 1j_1  1j_1 % 4) % sqr_big) -: qninv (  sqr_big  j.  sqr_big) , (-sqr_big) j.  sqr_big
(( 1j_1 _1j1  % 4) % sqr_big) -: qninv (  sqr_big  j.  sqr_big) ,   sqr_big  j. -sqr_big
(( 1j_1 _1j_1 % 4) % sqr_big) -: qninv (  sqr_big  j.  sqr_big) ,   sqr_big  j.  sqr_big
(( 1j_1  0j1  % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
(0 ,~  ( 1j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
(0 ,~  ( 1j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL
(0 ,~  ( 1j_1 % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(( 1j_1  0j_1 % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(( 1j1  _1j1  % 4) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(( 1j1  _1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(( 1j1  _1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL
(( 1j1  _1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(( 1j1  _1j_1 % 4) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(( 1    _1j1  % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(( 1    _1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(( 1    _1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL
(( 1    _1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(( 1    _1j_1 % 3) % FP_OVFL) -: qninv (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(( 1    _1j1  % 3) % FP_OVFL) -: qninv    FP_OVFL               ,   FP_OVFL  j. -FP_OVFL
(( 1    _1    % 2) % FP_OVFL) -: qninv    FP_OVFL               ,   FP_OVFL  j. -FP_UNFL
(( 1    _1    % 2) % FP_OVFL) -: qninv    FP_OVFL               ,   FP_OVFL
(( 1    _1    % 2) % FP_OVFL) -: qninv    FP_OVFL               ,   FP_OVFL  j.  FP_UNFL
(( 1    _1j_1 % 3) % FP_OVFL) -: qninv    FP_OVFL               ,   FP_OVFL  j.  FP_OVFL
(( 1    _1j1  % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(( 1    _1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(( 1    _1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL
(( 1    _1    % 2) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(( 1    _1j_1 % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(( 1j_1 _1j1  % 4) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(( 1j_1 _1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(( 1j_1 _1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL
(( 1j_1 _1    % 3) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(( 1j_1 _1j_1 % 4) % FP_OVFL) -: qninv (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
NB. - input without edge cases
(_1j1   1j1  % 4) -: qninv _1j_1 _1j_1
(_1j1   1    % 3) -: qninv _1j_1 _1
(_1j1   1j_1 % 4) -: qninv _1j_1 _1j1
(_1j1   0j1  % 3) -: qninv _1j_1  0j_1
(_1j1   0j_1 % 3) -: qninv _1j_1  0j1
(_1j1  _1j1  % 4) -: qninv _1j_1  1j_1
(_1j1  _1    % 3) -: qninv _1j_1  1
(_1j1  _1j_1 % 4) -: qninv _1j_1  1j1
(_1     1j1  % 3) -: qninv _1    _1j_1
(_1     1    % 2) -: qninv _1    _1
(_1     1j_1 % 3) -: qninv _1    _1j1
(_1     0j1  % 2) -: qninv _1     0j_1
(_1     0j_1 % 2) -: qninv _1     0j1
(_1    _1j1  % 3) -: qninv _1     1j_1
(_1    _1    % 2) -: qninv _1     1
(_1    _1j_1 % 3) -: qninv _1     1j1
(_1j_1  1j1  % 4) -: qninv _1j1  _1j_1
(_1j_1  1    % 3) -: qninv _1j1  _1
(_1j_1  1j_1 % 4) -: qninv _1j1  _1j1
(_1j_1  0j1  % 3) -: qninv _1j1   0j_1
(_1j_1  0j_1 % 3) -: qninv _1j1   0j1
(_1j_1 _1j1  % 4) -: qninv _1j1   1j_1
(_1j_1 _1    % 3) -: qninv _1j1   1
(_1j_1 _1j_1 % 4) -: qninv _1j1   1j1
( 0j1   1j1  % 3) -: qninv  0j_1 _1j_1
( 0j1   1    % 2) -: qninv  0j_1 _1
( 0j1   1j_1 % 3) -: qninv  0j_1 _1j1
( 0j1   0j1  % 2) -: qninv  0j_1  0j_1
( 0j1   0j_1 % 2) -: qninv  0j_1  0j1
( 0j1  _1j1  % 3) -: qninv  0j_1  1j_1
( 0j1  _1    % 2) -: qninv  0j_1  1
( 0j1  _1j_1 % 3) -: qninv  0j_1  1j1
( 0j_1  1j1  % 3) -: qninv  0j1  _1j_1
( 0j_1  1    % 2) -: qninv  0j1  _1
( 0j_1  1j_1 % 3) -: qninv  0j1  _1j1
( 0j_1  0j1  % 2) -: qninv  0j1   0j_1
( 0j_1  0j_1 % 2) -: qninv  0j1   0j1
( 0j_1  1    % 2) -: qninv  0j1  _1
( 0j_1  1j_1 % 3) -: qninv  0j1  _1j1
( 0j_1 _1    % 2) -: qninv  0j1   1
( 1j1   1j1  % 4) -: qninv  1j_1 _1j_1
( 1j1   1    % 3) -: qninv  1j_1 _1
( 1j1   1j_1 % 4) -: qninv  1j_1 _1j1
( 1j1   0j1  % 3) -: qninv  1j_1  0j_1
( 1j1   0j_1 % 3) -: qninv  1j_1  0j1
( 1j1  _1j1  % 4) -: qninv  1j_1  1j_1
( 1j1  _1    % 3) -: qninv  1j_1  1
( 1j1  _1j_1 % 4) -: qninv  1j_1  1j1
( 1     1j1  % 3) -: qninv  1    _1j_1
( 1     1    % 2) -: qninv  1    _1
( 1     1j_1 % 3) -: qninv  1    _1j1
( 1     0j1  % 2) -: qninv  1     0j_1
( 1     0j_1 % 2) -: qninv  1     0j1
( 1    _1j1  % 3) -: qninv  1     1j_1
( 1    _1    % 2) -: qninv  1     1
( 1    _1j_1 % 3) -: qninv  1     1j1
( 1j_1  1j1  % 4) -: qninv  1j1  _1j_1
( 1j_1  1    % 3) -: qninv  1j1  _1
( 1j_1  1j_1 % 4) -: qninv  1j1  _1j1
( 1j_1  0j1  % 3) -: qninv  1j1   0j_1
( 1j_1  0j_1 % 3) -: qninv  1j1   0j1
( 1j_1 _1j1  % 4) -: qninv  1j1   1j_1
( 1j_1 _1    % 3) -: qninv  1j1   1
( 1j_1 _1j_1 % 4) -: qninv  1j1   1j1
(1r2j_1r2 _1r2j_1r2 ,: 2j_2 _2j_2) -: qninv q3 ,: q4
NB. (-: qninv) i. 0 2  NB. this fails until https://github.com/jsoftware/jsource/issues/192 has been fixed

NB. qnmul
NB. - input contains NaN
1 1 -: isnan  0     0j_. qnmul  0     0j_.
1 1 -: isnan  0     0j_. qnmul  0    _.
1 1 -: isnan  0     0j_. qnmul  0j_.  0
1 1 -: isnan  0     0j_. qnmul _.     0
1 1 -: isnan  0    _.    qnmul  0     0j_.
1 1 -: isnan  0    _.    qnmul  0    _.
1 1 -: isnan  0    _.    qnmul  0j_.  0
1 1 -: isnan  0    _.    qnmul _.     0
1 1 -: isnan  0j_.  0    qnmul  0     0j_.
1 1 -: isnan  0j_.  0    qnmul  0    _.
1 1 -: isnan  0j_.  0    qnmul  0j_.  0
1 1 -: isnan  0j_.  0    qnmul _.     0
1 1 -: isnan _.     0    qnmul  0     0j_.
1 1 -: isnan _.     0    qnmul  0    _.
1 1 -: isnan _.     0    qnmul  0j_.  0
1 1 -: isnan _.     0    qnmul _.     0
1 1 -: isnan  0    _.j_. qnmul  0    _.j_.
1 1 -: isnan  0    _.j_. qnmul  0j_.  0j_.
1 1 -: isnan  0    _.j_. qnmul  0j_. _.
1 1 -: isnan  0    _.j_. qnmul _.     0j_.
1 1 -: isnan  0    _.j_. qnmul _.    _.
1 1 -: isnan  0    _.j_. qnmul _.j_.  0
1 1 -: isnan  0j_.  0j_. qnmul  0    _.j_.
1 1 -: isnan  0j_.  0j_. qnmul  0j_.  0j_.
1 1 -: isnan  0j_.  0j_. qnmul  0j_. _.
1 1 -: isnan  0j_.  0j_. qnmul _.     0j_.
1 1 -: isnan  0j_.  0j_. qnmul _.    _.
1 1 -: isnan  0j_.  0j_. qnmul _.j_.  0
1 1 -: isnan  0j_. _.    qnmul  0    _.j_.
1 1 -: isnan  0j_. _.    qnmul  0j_.  0j_.
1 1 -: isnan  0j_. _.    qnmul  0j_. _.
1 1 -: isnan  0j_. _.    qnmul _.     0j_.
1 1 -: isnan  0j_. _.    qnmul _.    _.
1 1 -: isnan  0j_. _.    qnmul _.j_.  0
1 1 -: isnan _.     0j_. qnmul  0    _.j_.
1 1 -: isnan _.     0j_. qnmul  0j_.  0j_.
1 1 -: isnan _.     0j_. qnmul  0j_. _.
1 1 -: isnan _.     0j_. qnmul _.     0j_.
1 1 -: isnan _.     0j_. qnmul _.    _.
1 1 -: isnan _.     0j_. qnmul _.j_.  0
1 1 -: isnan _.    _.    qnmul  0    _.j_.
1 1 -: isnan _.    _.    qnmul  0j_.  0j_.
1 1 -: isnan _.    _.    qnmul  0j_. _.
1 1 -: isnan _.    _.    qnmul _.     0j_.
1 1 -: isnan _.    _.    qnmul _.    _.
1 1 -: isnan _.    _.    qnmul _.j_.  0
1 1 -: isnan _.j_.  0    qnmul  0    _.j_.
1 1 -: isnan _.j_.  0    qnmul  0j_.  0j_.
1 1 -: isnan _.j_.  0    qnmul  0j_. _.
1 1 -: isnan _.j_.  0    qnmul _.     0j_.
1 1 -: isnan _.j_.  0    qnmul _.    _.
1 1 -: isnan _.j_.  0    qnmul _.j_.  0
1 1 -: isnan _.j_. _.    qnmul _.j_. _.
1 1 -: isnan _.j_. _.    qnmul _.j_.  0j_.
1 1 -: isnan _.j_. _.    qnmul _.    _.j_.
1 1 -: isnan _.j_. _.    qnmul  0j_. _.j_.
1 1 -: isnan _.j_.  0j_. qnmul _.j_. _.
1 1 -: isnan _.j_.  0j_. qnmul _.j_.  0j_.
1 1 -: isnan _.j_.  0j_. qnmul _.    _.j_.
1 1 -: isnan _.j_.  0j_. qnmul  0j_. _.j_.
1 1 -: isnan _.    _.j_. qnmul _.j_. _.
1 1 -: isnan _.    _.j_. qnmul _.j_.  0j_.
1 1 -: isnan _.    _.j_. qnmul _.    _.j_.
1 1 -: isnan _.    _.j_. qnmul  0j_. _.j_.
1 1 -: isnan  0j_. _.j_. qnmul _.j_. _.
1 1 -: isnan  0j_. _.j_. qnmul _.j_.  0j_.
1 1 -: isnan  0j_. _.j_. qnmul _.    _.j_.
1 1 -: isnan  0j_. _.j_. qnmul  0j_. _.j_.
1 1 -: isnan _.j_. _.j_. qnmul _.j_. _.j_.
1 1 -: isnan         q1  qnmul _.    _.
1 1 -: isnan         q1  qnmul _.j_. _.j_.
1 1 -: isnan _.    _.    qnmul        q2
1 1 -: isnan _.j_. _.j_. qnmul        q2
1 1 -: isnan (_. qn1 q1) qnmul        q2
1 1 -: isnan (_. qni q1) qnmul        q2
1 1 -: isnan (_. qnj q1) qnmul        q2
1 1 -: isnan (_. qnk q1) qnmul        q2
1 1 -: isnan         q1  qnmul _. qn1 q2
1 1 -: isnan         q1  qnmul _. qni q2
1 1 -: isnan         q1  qnmul _. qnj q2
1 1 -: isnan         q1  qnmul _. qnk q2
(,.~ 0 1) -: isnan  q1               qnmul q3 ,: _. qn1 q4
(,.~ 0 1) -: isnan  q1               qnmul q3 ,: _. qni q4
(,.~ 0 1) -: isnan  q1               qnmul q3 ,: _. qnj q4
(,.~ 0 1) -: isnan  q1               qnmul q3 ,: _. qnk q4
(,.~ 0 1) -: isnan (q1 ,: _. qn1 q2) qnmul q3
(,.~ 0 1) -: isnan (q1 ,: _. qni q2) qnmul q3
(,.~ 0 1) -: isnan (q1 ,: _. qnj q2) qnmul q3
(,.~ 0 1) -: isnan (q1 ,: _. qnk q2) qnmul q3
NB. - edge cases input
(_1j1 1j1 * FP_OVFL) -: q6 qnmul q6
NB. - input without edge cases
0 0 -: q1  qnmul 0 0
0 0 -: 0 0 qnmul q1
q1  -: q1  qnmul 1 0
q1  -: 1 0 qnmul q1
_60j12 30j24 -: q1 qnmul q2
_60j20 14j32 -: q2 qnmul q1
(2 2 $ _4j1 3j2 _1j1r4 3r4j1r2) -:  q1        qnmul q3 ,: q4
(2 2 $ _4j1 3j2 _8j5     7j6  ) -: (q1 ,: q2) qnmul q3
(2 2 $ _4j1 3j2 _2j5r4 7r4j3r2) -: (q1 ,: q2) qnmul q3 ,: q4

NB. qndivl
NB. - input contains NaN
1 1 -: isnan  0     0j_. qndivl  0     0j_.
1 1 -: isnan  0     0j_. qndivl  0    _.
1 1 -: isnan  0     0j_. qndivl  0j_.  0
1 1 -: isnan  0     0j_. qndivl _.     0
1 1 -: isnan  0    _.    qndivl  0     0j_.
1 1 -: isnan  0    _.    qndivl  0    _.
1 1 -: isnan  0    _.    qndivl  0j_.  0
1 1 -: isnan  0    _.    qndivl _.     0
1 1 -: isnan  0j_.  0    qndivl  0     0j_.
1 1 -: isnan  0j_.  0    qndivl  0    _.
1 1 -: isnan  0j_.  0    qndivl  0j_.  0
1 1 -: isnan  0j_.  0    qndivl _.     0
1 1 -: isnan _.     0    qndivl  0     0j_.
1 1 -: isnan _.     0    qndivl  0    _.
1 1 -: isnan _.     0    qndivl  0j_.  0
1 1 -: isnan _.     0    qndivl _.     0
1 1 -: isnan  0    _.j_. qndivl  0    _.j_.
1 1 -: isnan  0    _.j_. qndivl  0j_.  0j_.
1 1 -: isnan  0    _.j_. qndivl  0j_. _.
1 1 -: isnan  0    _.j_. qndivl _.     0j_.
1 1 -: isnan  0    _.j_. qndivl _.    _.
1 1 -: isnan  0    _.j_. qndivl _.j_.  0
1 1 -: isnan  0j_.  0j_. qndivl  0    _.j_.
1 1 -: isnan  0j_.  0j_. qndivl  0j_.  0j_.
1 1 -: isnan  0j_.  0j_. qndivl  0j_. _.
1 1 -: isnan  0j_.  0j_. qndivl _.     0j_.
1 1 -: isnan  0j_.  0j_. qndivl _.    _.
1 1 -: isnan  0j_.  0j_. qndivl _.j_.  0
1 1 -: isnan  0j_. _.    qndivl  0    _.j_.
1 1 -: isnan  0j_. _.    qndivl  0j_.  0j_.
1 1 -: isnan  0j_. _.    qndivl  0j_. _.
1 1 -: isnan  0j_. _.    qndivl _.     0j_.
1 1 -: isnan  0j_. _.    qndivl _.    _.
1 1 -: isnan  0j_. _.    qndivl _.j_.  0
1 1 -: isnan _.     0j_. qndivl  0    _.j_.
1 1 -: isnan _.     0j_. qndivl  0j_.  0j_.
1 1 -: isnan _.     0j_. qndivl  0j_. _.
1 1 -: isnan _.     0j_. qndivl _.     0j_.
1 1 -: isnan _.     0j_. qndivl _.    _.
1 1 -: isnan _.     0j_. qndivl _.j_.  0
1 1 -: isnan _.    _.    qndivl  0    _.j_.
1 1 -: isnan _.    _.    qndivl  0j_.  0j_.
1 1 -: isnan _.    _.    qndivl  0j_. _.
1 1 -: isnan _.    _.    qndivl _.     0j_.
1 1 -: isnan _.    _.    qndivl _.    _.
1 1 -: isnan _.    _.    qndivl _.j_.  0
1 1 -: isnan _.j_.  0    qndivl  0    _.j_.
1 1 -: isnan _.j_.  0    qndivl  0j_.  0j_.
1 1 -: isnan _.j_.  0    qndivl  0j_. _.
1 1 -: isnan _.j_.  0    qndivl _.     0j_.
1 1 -: isnan _.j_.  0    qndivl _.    _.
1 1 -: isnan _.j_.  0    qndivl _.j_.  0
1 1 -: isnan _.j_. _.    qndivl _.j_. _.
1 1 -: isnan _.j_. _.    qndivl _.j_.  0j_.
1 1 -: isnan _.j_. _.    qndivl _.    _.j_.
1 1 -: isnan _.j_. _.    qndivl  0j_. _.j_.
1 1 -: isnan _.j_.  0j_. qndivl _.j_. _.
1 1 -: isnan _.j_.  0j_. qndivl _.j_.  0j_.
1 1 -: isnan _.j_.  0j_. qndivl _.    _.j_.
1 1 -: isnan _.j_.  0j_. qndivl  0j_. _.j_.
1 1 -: isnan _.    _.j_. qndivl _.j_. _.
1 1 -: isnan _.    _.j_. qndivl _.j_.  0j_.
1 1 -: isnan _.    _.j_. qndivl _.    _.j_.
1 1 -: isnan _.    _.j_. qndivl  0j_. _.j_.
1 1 -: isnan  0j_. _.j_. qndivl _.j_. _.
1 1 -: isnan  0j_. _.j_. qndivl _.j_.  0j_.
1 1 -: isnan  0j_. _.j_. qndivl _.    _.j_.
1 1 -: isnan  0j_. _.j_. qndivl  0j_. _.j_.
1 1 -: isnan _.j_. _.j_. qndivl _.j_. _.j_.
1 1 -: isnan         q1  qndivl _.    _.
1 1 -: isnan         q1  qndivl _.j_. _.j_.
1 1 -: isnan _.    _.    qndivl        q2
1 1 -: isnan _.j_. _.j_. qndivl        q2
1 1 -: isnan (_. qn1 q1) qndivl        q2
1 1 -: isnan (_. qni q1) qndivl        q2
1 1 -: isnan (_. qnj q1) qndivl        q2
1 1 -: isnan (_. qnk q1) qndivl        q2
1 1 -: isnan         q1  qndivl _. qn1 q2
1 1 -: isnan         q1  qndivl _. qni q2
1 1 -: isnan         q1  qndivl _. qnj q2
1 1 -: isnan         q1  qndivl _. qnk q2
(,.~ 0 1) -: isnan  q1               qndivl q3 ,: _. qn1 q4
(,.~ 0 1) -: isnan  q1               qndivl q3 ,: _. qni q4
(,.~ 0 1) -: isnan  q1               qndivl q3 ,: _. qnj q4
(,.~ 0 1) -: isnan  q1               qndivl q3 ,: _. qnk q4
(,.~ 0 1) -: isnan (q1 ,: _. qn1 q2) qndivl q3
(,.~ 0 1) -: isnan (q1 ,: _. qni q2) qndivl q3
(,.~ 0 1) -: isnan (q1 ,: _. qnj q2) qndivl q3
(,.~ 0 1) -: isnan (q1 ,: _. qnk q2) qndivl q3
NB. - denominator is zero so NaN error must be throwed
q1 0:@qndivl :: 1: 0 0
NB. - input contains q∞ so NaN error must be throwed
__j__ __j__ 0:@qndivl :: 1 q1
__j__ __j_1 0:@qndivl :: 1 q1
__j__ __    0:@qndivl :: 1 q1
__j__ __j1  0:@qndivl :: 1 q1
__j__ __j_  0:@qndivl :: 1 q1
__j__ _1j__ 0:@qndivl :: 1 q1
__j__ _1j_1 0:@qndivl :: 1 q1
__j__ _1    0:@qndivl :: 1 q1
__j__ _1j1  0:@qndivl :: 1 q1
__j__ _1j_  0:@qndivl :: 1 q1
__j__  0j__ 0:@qndivl :: 1 q1
__j__  0j_1 0:@qndivl :: 1 q1
__j__  0    0:@qndivl :: 1 q1
__j__  0j1  0:@qndivl :: 1 q1
__j__  0j_  0:@qndivl :: 1 q1
__j__  1j__ 0:@qndivl :: 1 q1
__j__  1j_1 0:@qndivl :: 1 q1
__j__  1    0:@qndivl :: 1 q1
__j__  1j1  0:@qndivl :: 1 q1
__j__  1j_  0:@qndivl :: 1 q1
__j__  _j__ 0:@qndivl :: 1 q1
__j__  _j_1 0:@qndivl :: 1 q1
__j__  _    0:@qndivl :: 1 q1
__j__  _j1  0:@qndivl :: 1 q1
__j__  _j_  0:@qndivl :: 1 q1
__j_1 __j__ 0:@qndivl :: 1 q1
__j_1 __j_1 0:@qndivl :: 1 q1
__j_1 __    0:@qndivl :: 1 q1
__j_1 __j1  0:@qndivl :: 1 q1
__j_1 __j_  0:@qndivl :: 1 q1
__j_1 _1j__ 0:@qndivl :: 1 q1
__j_1 _1j_  0:@qndivl :: 1 q1
__j_1  0j__ 0:@qndivl :: 1 q1
__j_1  0j_  0:@qndivl :: 1 q1
__j_1  1j__ 0:@qndivl :: 1 q1
__j_1  1j_  0:@qndivl :: 1 q1
__j_1  _j__ 0:@qndivl :: 1 q1
__j_1  _j_1 0:@qndivl :: 1 q1
__j_1  _    0:@qndivl :: 1 q1
__j_1  _j1  0:@qndivl :: 1 q1
__j_1  _j_  0:@qndivl :: 1 q1
__    __j__ 0:@qndivl :: 1 q1
__    __j_1 0:@qndivl :: 1 q1
__    __    0:@qndivl :: 1 q1
__    __j1  0:@qndivl :: 1 q1
__    __j_  0:@qndivl :: 1 q1
__    _1j__ 0:@qndivl :: 1 q1
__    _1j_  0:@qndivl :: 1 q1
__     0j__ 0:@qndivl :: 1 q1
__     0j_  0:@qndivl :: 1 q1
__     1j__ 0:@qndivl :: 1 q1
__     1j_  0:@qndivl :: 1 q1
__     _j__ 0:@qndivl :: 1 q1
__     _j_1 0:@qndivl :: 1 q1
__     _    0:@qndivl :: 1 q1
__     _j1  0:@qndivl :: 1 q1
__     _j_  0:@qndivl :: 1 q1
__j1  __j__ 0:@qndivl :: 1 q1
__j1  __j_1 0:@qndivl :: 1 q1
__j1  __    0:@qndivl :: 1 q1
__j1  __j1  0:@qndivl :: 1 q1
__j1  __j_  0:@qndivl :: 1 q1
__j1  _1j__ 0:@qndivl :: 1 q1
__j1  _1j_  0:@qndivl :: 1 q1
__j1   0j__ 0:@qndivl :: 1 q1
__j1   0j_  0:@qndivl :: 1 q1
__j1   1j__ 0:@qndivl :: 1 q1
__j1   1j_  0:@qndivl :: 1 q1
__j1   _j__ 0:@qndivl :: 1 q1
__j1   _j_1 0:@qndivl :: 1 q1
__j1   _    0:@qndivl :: 1 q1
__j1   _j1  0:@qndivl :: 1 q1
__j1   _j_  0:@qndivl :: 1 q1
__j_  __j__ 0:@qndivl :: 1 q1
__j_  __j_1 0:@qndivl :: 1 q1
__j_  __    0:@qndivl :: 1 q1
__j_  __j1  0:@qndivl :: 1 q1
__j_  __j_  0:@qndivl :: 1 q1
__j_  _1j__ 0:@qndivl :: 1 q1
__j_  _1j_1 0:@qndivl :: 1 q1
__j_  _1    0:@qndivl :: 1 q1
__j_  _1j1  0:@qndivl :: 1 q1
__j_  _1j_  0:@qndivl :: 1 q1
__j_   0j__ 0:@qndivl :: 1 q1
__j_   0j_1 0:@qndivl :: 1 q1
__j_   0    0:@qndivl :: 1 q1
__j_   0j1  0:@qndivl :: 1 q1
__j_   0j_  0:@qndivl :: 1 q1
__j_   1j__ 0:@qndivl :: 1 q1
__j_   1j_1 0:@qndivl :: 1 q1
__j_   1    0:@qndivl :: 1 q1
__j_   1j1  0:@qndivl :: 1 q1
__j_   1j_  0:@qndivl :: 1 q1
__j_   _j__ 0:@qndivl :: 1 q1
__j_   _j_1 0:@qndivl :: 1 q1
__j_   _    0:@qndivl :: 1 q1
__j_   _j1  0:@qndivl :: 1 q1
__j_   _j_  0:@qndivl :: 1 q1
_1j__ __j__ 0:@qndivl :: 1 q1
_1j__ __j_1 0:@qndivl :: 1 q1
_1j__ __    0:@qndivl :: 1 q1
_1j__ __j1  0:@qndivl :: 1 q1
_1j__ __j_  0:@qndivl :: 1 q1
_1j__ _1j__ 0:@qndivl :: 1 q1
_1j__ _1j_  0:@qndivl :: 1 q1
_1j__  0j__ 0:@qndivl :: 1 q1
_1j__  0j_  0:@qndivl :: 1 q1
_1j__  1j__ 0:@qndivl :: 1 q1
_1j__  1j_  0:@qndivl :: 1 q1
_1j__  _j__ 0:@qndivl :: 1 q1
_1j__  _j_1 0:@qndivl :: 1 q1
_1j__  _    0:@qndivl :: 1 q1
_1j__  _j1  0:@qndivl :: 1 q1
_1j__  _j_  0:@qndivl :: 1 q1
_1j_1 __j__ 0:@qndivl :: 1 q1
_1j_1 __j_  0:@qndivl :: 1 q1
_1j_1  _j__ 0:@qndivl :: 1 q1
_1j_1  _j_  0:@qndivl :: 1 q1
_1    __j__ 0:@qndivl :: 1 q1
_1    __j_  0:@qndivl :: 1 q1
_1     _j__ 0:@qndivl :: 1 q1
_1     _j_  0:@qndivl :: 1 q1
_1j1  __j__ 0:@qndivl :: 1 q1
_1j1  __j_  0:@qndivl :: 1 q1
_1j1   _j__ 0:@qndivl :: 1 q1
_1j1   _j_  0:@qndivl :: 1 q1
_1j_  __j__ 0:@qndivl :: 1 q1
_1j_  __j_1 0:@qndivl :: 1 q1
_1j_  __    0:@qndivl :: 1 q1
_1j_  __j1  0:@qndivl :: 1 q1
_1j_  __j_  0:@qndivl :: 1 q1
_1j_  _1j__ 0:@qndivl :: 1 q1
_1j_  _1j_  0:@qndivl :: 1 q1
_1j_   0j__ 0:@qndivl :: 1 q1
_1j_   0j_  0:@qndivl :: 1 q1
_1j_   1j__ 0:@qndivl :: 1 q1
_1j_   1j_  0:@qndivl :: 1 q1
_1j_   _j__ 0:@qndivl :: 1 q1
_1j_   _j_1 0:@qndivl :: 1 q1
_1j_   _    0:@qndivl :: 1 q1
_1j_   _j1  0:@qndivl :: 1 q1
_1j_   _j_  0:@qndivl :: 1 q1
 0j__ __j__ 0:@qndivl :: 1 q1
 0j__ __j_1 0:@qndivl :: 1 q1
 0j__ __    0:@qndivl :: 1 q1
 0j__ __j1  0:@qndivl :: 1 q1
 0j__ __j_  0:@qndivl :: 1 q1
 0j__ _1j__ 0:@qndivl :: 1 q1
 0j__ _1j_  0:@qndivl :: 1 q1
 0j__  0j__ 0:@qndivl :: 1 q1
 0j__  0j_  0:@qndivl :: 1 q1
 0j__  1j__ 0:@qndivl :: 1 q1
 0j__  1j_  0:@qndivl :: 1 q1
 0j__  _j__ 0:@qndivl :: 1 q1
 0j__  _j_1 0:@qndivl :: 1 q1
 0j__  _    0:@qndivl :: 1 q1
 0j__  _j1  0:@qndivl :: 1 q1
 0j__  _j_  0:@qndivl :: 1 q1
 0j_1 __j__ 0:@qndivl :: 1 q1
 0j_1 __j_  0:@qndivl :: 1 q1
 0j_1  _j__ 0:@qndivl :: 1 q1
 0j_1  _j_  0:@qndivl :: 1 q1
 0    __j__ 0:@qndivl :: 1 q1
 0    __j_  0:@qndivl :: 1 q1
 0     _j__ 0:@qndivl :: 1 q1
 0     _j_  0:@qndivl :: 1 q1
 0j1  __j__ 0:@qndivl :: 1 q1
 0j1  __j_  0:@qndivl :: 1 q1
 0j1   _j__ 0:@qndivl :: 1 q1
 0j1   _j_  0:@qndivl :: 1 q1
 0j_  __j__ 0:@qndivl :: 1 q1
 0j_  __j_1 0:@qndivl :: 1 q1
 0j_  __    0:@qndivl :: 1 q1
 0j_  __j1  0:@qndivl :: 1 q1
 0j_  __j_  0:@qndivl :: 1 q1
 0j_  _1j__ 0:@qndivl :: 1 q1
 0j_  _1j_  0:@qndivl :: 1 q1
 0j_   0j__ 0:@qndivl :: 1 q1
 0j_   0j_  0:@qndivl :: 1 q1
 0j_   1j__ 0:@qndivl :: 1 q1
 0j_   1j_  0:@qndivl :: 1 q1
 0j_   _j__ 0:@qndivl :: 1 q1
 0j_   _j_1 0:@qndivl :: 1 q1
 0j_   _    0:@qndivl :: 1 q1
 0j_   _j1  0:@qndivl :: 1 q1
 0j_   _j_  0:@qndivl :: 1 q1
 1j__ __j__ 0:@qndivl :: 1 q1
 1j__ __j_1 0:@qndivl :: 1 q1
 1j__ __    0:@qndivl :: 1 q1
 1j__ __j1  0:@qndivl :: 1 q1
 1j__ __j_  0:@qndivl :: 1 q1
 1j__ _1j__ 0:@qndivl :: 1 q1
 1j__ _1j_  0:@qndivl :: 1 q1
 1j__  0j__ 0:@qndivl :: 1 q1
 1j__  0j_  0:@qndivl :: 1 q1
 1j__  1j__ 0:@qndivl :: 1 q1
 1j__  1j_  0:@qndivl :: 1 q1
 1j__  _j__ 0:@qndivl :: 1 q1
 1j__  _j_1 0:@qndivl :: 1 q1
 1j__  _    0:@qndivl :: 1 q1
 1j__  _j1  0:@qndivl :: 1 q1
 1j__  _j_  0:@qndivl :: 1 q1
 1j_1 __j__ 0:@qndivl :: 1 q1
 1j_1 __j_  0:@qndivl :: 1 q1
 1j_1  _j__ 0:@qndivl :: 1 q1
 1j_1  _j_  0:@qndivl :: 1 q1
 1    __j__ 0:@qndivl :: 1 q1
 1    __j_  0:@qndivl :: 1 q1
 1     _j__ 0:@qndivl :: 1 q1
 1     _j_  0:@qndivl :: 1 q1
 1j1  __j__ 0:@qndivl :: 1 q1
 1j1  __j_  0:@qndivl :: 1 q1
 1j1   _j__ 0:@qndivl :: 1 q1
 1j1   _j_  0:@qndivl :: 1 q1
 1j_  __j__ 0:@qndivl :: 1 q1
 1j_  __j_1 0:@qndivl :: 1 q1
 1j_  __    0:@qndivl :: 1 q1
 1j_  __j1  0:@qndivl :: 1 q1
 1j_  __j_  0:@qndivl :: 1 q1
 1j_  _1j__ 0:@qndivl :: 1 q1
 1j_  _1j_  0:@qndivl :: 1 q1
 1j_   0j__ 0:@qndivl :: 1 q1
 1j_   0j_  0:@qndivl :: 1 q1
 1j_   1j__ 0:@qndivl :: 1 q1
 1j_   1j_  0:@qndivl :: 1 q1
 1j_   _j__ 0:@qndivl :: 1 q1
 1j_   _j_1 0:@qndivl :: 1 q1
 1j_   _    0:@qndivl :: 1 q1
 1j_   _j1  0:@qndivl :: 1 q1
 1j_   _j_  0:@qndivl :: 1 q1
 _j__ __j__ 0:@qndivl :: 1 q1
 _j__ __j_1 0:@qndivl :: 1 q1
 _j__ __    0:@qndivl :: 1 q1
 _j__ __j1  0:@qndivl :: 1 q1
 _j__ __j_  0:@qndivl :: 1 q1
 _j__ _1j__ 0:@qndivl :: 1 q1
 _j__ _1j_1 0:@qndivl :: 1 q1
 _j__ _1    0:@qndivl :: 1 q1
 _j__ _1j1  0:@qndivl :: 1 q1
 _j__ _1j_  0:@qndivl :: 1 q1
 _j__  0j__ 0:@qndivl :: 1 q1
 _j__  0j_1 0:@qndivl :: 1 q1
 _j__  0    0:@qndivl :: 1 q1
 _j__  0j1  0:@qndivl :: 1 q1
 _j__  0j_  0:@qndivl :: 1 q1
 _j__  1j__ 0:@qndivl :: 1 q1
 _j__  1j_1 0:@qndivl :: 1 q1
 _j__  1    0:@qndivl :: 1 q1
 _j__  1j1  0:@qndivl :: 1 q1
 _j__  1j_  0:@qndivl :: 1 q1
 _j__  _j__ 0:@qndivl :: 1 q1
 _j__  _j_1 0:@qndivl :: 1 q1
 _j__  _    0:@qndivl :: 1 q1
 _j__  _j1  0:@qndivl :: 1 q1
 _j__  _j_  0:@qndivl :: 1 q1
 _j_1 __j__ 0:@qndivl :: 1 q1
 _j_1 __j_1 0:@qndivl :: 1 q1
 _j_1 __    0:@qndivl :: 1 q1
 _j_1 __j1  0:@qndivl :: 1 q1
 _j_1 __j_  0:@qndivl :: 1 q1
 _j_1 _1j__ 0:@qndivl :: 1 q1
 _j_1 _1j_  0:@qndivl :: 1 q1
 _j_1  0j__ 0:@qndivl :: 1 q1
 _j_1  0j_  0:@qndivl :: 1 q1
 _j_1  1j__ 0:@qndivl :: 1 q1
 _j_1  1j_  0:@qndivl :: 1 q1
 _j_1  _j__ 0:@qndivl :: 1 q1
 _j_1  _j_1 0:@qndivl :: 1 q1
 _j_1  _    0:@qndivl :: 1 q1
 _j_1  _j1  0:@qndivl :: 1 q1
 _j_1  _j_  0:@qndivl :: 1 q1
 _    __j__ 0:@qndivl :: 1 q1
 _    __j_1 0:@qndivl :: 1 q1
 _    __    0:@qndivl :: 1 q1
 _    __j1  0:@qndivl :: 1 q1
 _    __j_  0:@qndivl :: 1 q1
 _    _1j__ 0:@qndivl :: 1 q1
 _    _1j_  0:@qndivl :: 1 q1
 _     0j__ 0:@qndivl :: 1 q1
 _     0j_  0:@qndivl :: 1 q1
 _     1j__ 0:@qndivl :: 1 q1
 _     1j_  0:@qndivl :: 1 q1
 _     _j__ 0:@qndivl :: 1 q1
 _     _j_1 0:@qndivl :: 1 q1
 _     _    0:@qndivl :: 1 q1
 _     _j1  0:@qndivl :: 1 q1
 _     _j_  0:@qndivl :: 1 q1
 _j1  __j__ 0:@qndivl :: 1 q1
 _j1  __j_1 0:@qndivl :: 1 q1
 _j1  __    0:@qndivl :: 1 q1
 _j1  __j1  0:@qndivl :: 1 q1
 _j1  __j_  0:@qndivl :: 1 q1
 _j1  _1j__ 0:@qndivl :: 1 q1
 _j1  _1j_  0:@qndivl :: 1 q1
 _j1   0j__ 0:@qndivl :: 1 q1
 _j1   0j_  0:@qndivl :: 1 q1
 _j1   1j__ 0:@qndivl :: 1 q1
 _j1   1j_  0:@qndivl :: 1 q1
 _j1   _j__ 0:@qndivl :: 1 q1
 _j1   _j_1 0:@qndivl :: 1 q1
 _j1   _    0:@qndivl :: 1 q1
 _j1   _j1  0:@qndivl :: 1 q1
 _j1   _j_  0:@qndivl :: 1 q1
 _j_  __j__ 0:@qndivl :: 1 q1
 _j_  __j_1 0:@qndivl :: 1 q1
 _j_  __    0:@qndivl :: 1 q1
 _j_  __j1  0:@qndivl :: 1 q1
 _j_  __j_  0:@qndivl :: 1 q1
 _j_  _1j__ 0:@qndivl :: 1 q1
 _j_  _1j_1 0:@qndivl :: 1 q1
 _j_  _1    0:@qndivl :: 1 q1
 _j_  _1j1  0:@qndivl :: 1 q1
 _j_  _1j_  0:@qndivl :: 1 q1
 _j_   0j__ 0:@qndivl :: 1 q1
 _j_   0j_1 0:@qndivl :: 1 q1
 _j_   0    0:@qndivl :: 1 q1
 _j_   0j1  0:@qndivl :: 1 q1
 _j_   0j_  0:@qndivl :: 1 q1
 _j_   1j__ 0:@qndivl :: 1 q1
 _j_   1j_1 0:@qndivl :: 1 q1
 _j_   1    0:@qndivl :: 1 q1
 _j_   1j1  0:@qndivl :: 1 q1
 _j_   1j_  0:@qndivl :: 1 q1
 _j_   _j__ 0:@qndivl :: 1 q1
 _j_   _j_1 0:@qndivl :: 1 q1
 _j_   _    0:@qndivl :: 1 q1
 _j_   _j1  0:@qndivl :: 1 q1
 _j_   _j_  0:@qndivl :: 1 q1
q1 0:@qndivl :: 1: __j__ __j__
q1 0:@qndivl :: 1: __j__ __j_1
q1 0:@qndivl :: 1: __j__ __
q1 0:@qndivl :: 1: __j__ __j1
q1 0:@qndivl :: 1: __j__ __j_
q1 0:@qndivl :: 1: __j__ _1j__
q1 0:@qndivl :: 1: __j__ _1j_1
q1 0:@qndivl :: 1: __j__ _1
q1 0:@qndivl :: 1: __j__ _1j1
q1 0:@qndivl :: 1: __j__ _1j_
q1 0:@qndivl :: 1: __j__  0j__
q1 0:@qndivl :: 1: __j__  0j_1
q1 0:@qndivl :: 1: __j__  0
q1 0:@qndivl :: 1: __j__  0j1
q1 0:@qndivl :: 1: __j__  0j_
q1 0:@qndivl :: 1: __j__  1j__
q1 0:@qndivl :: 1: __j__  1j_1
q1 0:@qndivl :: 1: __j__  1
q1 0:@qndivl :: 1: __j__  1j1
q1 0:@qndivl :: 1: __j__  1j_
q1 0:@qndivl :: 1: __j__  _j__
q1 0:@qndivl :: 1: __j__  _j_1
q1 0:@qndivl :: 1: __j__  _
q1 0:@qndivl :: 1: __j__  _j1
q1 0:@qndivl :: 1: __j__  _j_
q1 0:@qndivl :: 1: __j_1 __j__
q1 0:@qndivl :: 1: __j_1 __j_1
q1 0:@qndivl :: 1: __j_1 __
q1 0:@qndivl :: 1: __j_1 __j1
q1 0:@qndivl :: 1: __j_1 __j_
q1 0:@qndivl :: 1: __j_1 _1j__
q1 0:@qndivl :: 1: __j_1 _1j_
q1 0:@qndivl :: 1: __j_1  0j__
q1 0:@qndivl :: 1: __j_1  0j_
q1 0:@qndivl :: 1: __j_1  1j__
q1 0:@qndivl :: 1: __j_1  1j_
q1 0:@qndivl :: 1: __j_1  _j__
q1 0:@qndivl :: 1: __j_1  _j_1
q1 0:@qndivl :: 1: __j_1  _
q1 0:@qndivl :: 1: __j_1  _j1
q1 0:@qndivl :: 1: __j_1  _j_
q1 0:@qndivl :: 1: __    __j__
q1 0:@qndivl :: 1: __    __j_1
q1 0:@qndivl :: 1: __    __
q1 0:@qndivl :: 1: __    __j1
q1 0:@qndivl :: 1: __    __j_
q1 0:@qndivl :: 1: __    _1j__
q1 0:@qndivl :: 1: __    _1j_
q1 0:@qndivl :: 1: __     0j__
q1 0:@qndivl :: 1: __     0j_
q1 0:@qndivl :: 1: __     1j__
q1 0:@qndivl :: 1: __     1j_
q1 0:@qndivl :: 1: __     _j__
q1 0:@qndivl :: 1: __     _j_1
q1 0:@qndivl :: 1: __     _
q1 0:@qndivl :: 1: __     _j1
q1 0:@qndivl :: 1: __     _j_
q1 0:@qndivl :: 1: __j1  __j__
q1 0:@qndivl :: 1: __j1  __j_1
q1 0:@qndivl :: 1: __j1  __
q1 0:@qndivl :: 1: __j1  __j1
q1 0:@qndivl :: 1: __j1  __j_
q1 0:@qndivl :: 1: __j1  _1j__
q1 0:@qndivl :: 1: __j1  _1j_
q1 0:@qndivl :: 1: __j1   0j__
q1 0:@qndivl :: 1: __j1   0j_
q1 0:@qndivl :: 1: __j1   1j__
q1 0:@qndivl :: 1: __j1   1j_
q1 0:@qndivl :: 1: __j1   _j__
q1 0:@qndivl :: 1: __j1   _j_1
q1 0:@qndivl :: 1: __j1   _
q1 0:@qndivl :: 1: __j1   _j1
q1 0:@qndivl :: 1: __j1   _j_
q1 0:@qndivl :: 1: __j_  __j__
q1 0:@qndivl :: 1: __j_  __j_1
q1 0:@qndivl :: 1: __j_  __
q1 0:@qndivl :: 1: __j_  __j1
q1 0:@qndivl :: 1: __j_  __j_
q1 0:@qndivl :: 1: __j_  _1j__
q1 0:@qndivl :: 1: __j_  _1j_1
q1 0:@qndivl :: 1: __j_  _1
q1 0:@qndivl :: 1: __j_  _1j1
q1 0:@qndivl :: 1: __j_  _1j_
q1 0:@qndivl :: 1: __j_   0j__
q1 0:@qndivl :: 1: __j_   0j_1
q1 0:@qndivl :: 1: __j_   0
q1 0:@qndivl :: 1: __j_   0j1
q1 0:@qndivl :: 1: __j_   0j_
q1 0:@qndivl :: 1: __j_   1j__
q1 0:@qndivl :: 1: __j_   1j_1
q1 0:@qndivl :: 1: __j_   1
q1 0:@qndivl :: 1: __j_   1j1
q1 0:@qndivl :: 1: __j_   1j_
q1 0:@qndivl :: 1: __j_   _j__
q1 0:@qndivl :: 1: __j_   _j_1
q1 0:@qndivl :: 1: __j_   _
q1 0:@qndivl :: 1: __j_   _j1
q1 0:@qndivl :: 1: __j_   _j_
q1 0:@qndivl :: 1: _1j__ __j__
q1 0:@qndivl :: 1: _1j__ __j_1
q1 0:@qndivl :: 1: _1j__ __
q1 0:@qndivl :: 1: _1j__ __j1
q1 0:@qndivl :: 1: _1j__ __j_
q1 0:@qndivl :: 1: _1j__ _1j__
q1 0:@qndivl :: 1: _1j__ _1j_
q1 0:@qndivl :: 1: _1j__  0j__
q1 0:@qndivl :: 1: _1j__  0j_
q1 0:@qndivl :: 1: _1j__  1j__
q1 0:@qndivl :: 1: _1j__  1j_
q1 0:@qndivl :: 1: _1j__  _j__
q1 0:@qndivl :: 1: _1j__  _j_1
q1 0:@qndivl :: 1: _1j__  _
q1 0:@qndivl :: 1: _1j__  _j1
q1 0:@qndivl :: 1: _1j__  _j_
q1 0:@qndivl :: 1: _1j_1 __j__
q1 0:@qndivl :: 1: _1j_1 __j_
q1 0:@qndivl :: 1: _1j_1  _j__
q1 0:@qndivl :: 1: _1j_1  _j_
q1 0:@qndivl :: 1: _1    __j__
q1 0:@qndivl :: 1: _1    __j_
q1 0:@qndivl :: 1: _1     _j__
q1 0:@qndivl :: 1: _1     _j_
q1 0:@qndivl :: 1: _1j1  __j__
q1 0:@qndivl :: 1: _1j1  __j_
q1 0:@qndivl :: 1: _1j1   _j__
q1 0:@qndivl :: 1: _1j1   _j_
q1 0:@qndivl :: 1: _1j_  __j__
q1 0:@qndivl :: 1: _1j_  __j_1
q1 0:@qndivl :: 1: _1j_  __
q1 0:@qndivl :: 1: _1j_  __j1
q1 0:@qndivl :: 1: _1j_  __j_
q1 0:@qndivl :: 1: _1j_  _1j__
q1 0:@qndivl :: 1: _1j_  _1j_
q1 0:@qndivl :: 1: _1j_   0j__
q1 0:@qndivl :: 1: _1j_   0j_
q1 0:@qndivl :: 1: _1j_   1j__
q1 0:@qndivl :: 1: _1j_   1j_
q1 0:@qndivl :: 1: _1j_   _j__
q1 0:@qndivl :: 1: _1j_   _j_1
q1 0:@qndivl :: 1: _1j_   _
q1 0:@qndivl :: 1: _1j_   _j1
q1 0:@qndivl :: 1: _1j_   _j_
q1 0:@qndivl :: 1:  0j__ __j__
q1 0:@qndivl :: 1:  0j__ __j_1
q1 0:@qndivl :: 1:  0j__ __
q1 0:@qndivl :: 1:  0j__ __j1
q1 0:@qndivl :: 1:  0j__ __j_
q1 0:@qndivl :: 1:  0j__ _1j__
q1 0:@qndivl :: 1:  0j__ _1j_
q1 0:@qndivl :: 1:  0j__  0j__
q1 0:@qndivl :: 1:  0j__  0j_
q1 0:@qndivl :: 1:  0j__  1j__
q1 0:@qndivl :: 1:  0j__  1j_
q1 0:@qndivl :: 1:  0j__  _j__
q1 0:@qndivl :: 1:  0j__  _j_1
q1 0:@qndivl :: 1:  0j__  _
q1 0:@qndivl :: 1:  0j__  _j1
q1 0:@qndivl :: 1:  0j__  _j_
q1 0:@qndivl :: 1:  0j_1 __j__
q1 0:@qndivl :: 1:  0j_1 __j_
q1 0:@qndivl :: 1:  0j_1  _j__
q1 0:@qndivl :: 1:  0j_1  _j_
q1 0:@qndivl :: 1:  0    __j__
q1 0:@qndivl :: 1:  0    __j_
q1 0:@qndivl :: 1:  0     _j__
q1 0:@qndivl :: 1:  0     _j_
q1 0:@qndivl :: 1:  0j1  __j__
q1 0:@qndivl :: 1:  0j1  __j_
q1 0:@qndivl :: 1:  0j1   _j__
q1 0:@qndivl :: 1:  0j1   _j_
q1 0:@qndivl :: 1:  0j_  __j__
q1 0:@qndivl :: 1:  0j_  __j_1
q1 0:@qndivl :: 1:  0j_  __
q1 0:@qndivl :: 1:  0j_  __j1
q1 0:@qndivl :: 1:  0j_  __j_
q1 0:@qndivl :: 1:  0j_  _1j__
q1 0:@qndivl :: 1:  0j_  _1j_
q1 0:@qndivl :: 1:  0j_   0j__
q1 0:@qndivl :: 1:  0j_   0j_
q1 0:@qndivl :: 1:  0j_   1j__
q1 0:@qndivl :: 1:  0j_   1j_
q1 0:@qndivl :: 1:  0j_   _j__
q1 0:@qndivl :: 1:  0j_   _j_1
q1 0:@qndivl :: 1:  0j_   _
q1 0:@qndivl :: 1:  0j_   _j1
q1 0:@qndivl :: 1:  0j_   _j_
q1 0:@qndivl :: 1:  1j__ __j__
q1 0:@qndivl :: 1:  1j__ __j_1
q1 0:@qndivl :: 1:  1j__ __
q1 0:@qndivl :: 1:  1j__ __j1
q1 0:@qndivl :: 1:  1j__ __j_
q1 0:@qndivl :: 1:  1j__ _1j__
q1 0:@qndivl :: 1:  1j__ _1j_
q1 0:@qndivl :: 1:  1j__  0j__
q1 0:@qndivl :: 1:  1j__  0j_
q1 0:@qndivl :: 1:  1j__  1j__
q1 0:@qndivl :: 1:  1j__  1j_
q1 0:@qndivl :: 1:  1j__  _j__
q1 0:@qndivl :: 1:  1j__  _j_1
q1 0:@qndivl :: 1:  1j__  _
q1 0:@qndivl :: 1:  1j__  _j1
q1 0:@qndivl :: 1:  1j__  _j_
q1 0:@qndivl :: 1:  1j_1 __j__
q1 0:@qndivl :: 1:  1j_1 __j_
q1 0:@qndivl :: 1:  1j_1  _j__
q1 0:@qndivl :: 1:  1j_1  _j_
q1 0:@qndivl :: 1:  1    __j__
q1 0:@qndivl :: 1:  1    __j_
q1 0:@qndivl :: 1:  1     _j__
q1 0:@qndivl :: 1:  1     _j_
q1 0:@qndivl :: 1:  1j1  __j__
q1 0:@qndivl :: 1:  1j1  __j_
q1 0:@qndivl :: 1:  1j1   _j__
q1 0:@qndivl :: 1:  1j1   _j_
q1 0:@qndivl :: 1:  1j_  __j__
q1 0:@qndivl :: 1:  1j_  __j_1
q1 0:@qndivl :: 1:  1j_  __
q1 0:@qndivl :: 1:  1j_  __j1
q1 0:@qndivl :: 1:  1j_  __j_
q1 0:@qndivl :: 1:  1j_  _1j__
q1 0:@qndivl :: 1:  1j_  _1j_
q1 0:@qndivl :: 1:  1j_   0j__
q1 0:@qndivl :: 1:  1j_   0j_
q1 0:@qndivl :: 1:  1j_   1j__
q1 0:@qndivl :: 1:  1j_   1j_
q1 0:@qndivl :: 1:  1j_   _j__
q1 0:@qndivl :: 1:  1j_   _j_1
q1 0:@qndivl :: 1:  1j_   _
q1 0:@qndivl :: 1:  1j_   _j1
q1 0:@qndivl :: 1:  1j_   _j_
q1 0:@qndivl :: 1:  _j__ __j__
q1 0:@qndivl :: 1:  _j__ __j_1
q1 0:@qndivl :: 1:  _j__ __
q1 0:@qndivl :: 1:  _j__ __j1
q1 0:@qndivl :: 1:  _j__ __j_
q1 0:@qndivl :: 1:  _j__ _1j__
q1 0:@qndivl :: 1:  _j__ _1j_1
q1 0:@qndivl :: 1:  _j__ _1
q1 0:@qndivl :: 1:  _j__ _1j1
q1 0:@qndivl :: 1:  _j__ _1j_
q1 0:@qndivl :: 1:  _j__  0j__
q1 0:@qndivl :: 1:  _j__  0j_1
q1 0:@qndivl :: 1:  _j__  0
q1 0:@qndivl :: 1:  _j__  0j1
q1 0:@qndivl :: 1:  _j__  0j_
q1 0:@qndivl :: 1:  _j__  1j__
q1 0:@qndivl :: 1:  _j__  1j_1
q1 0:@qndivl :: 1:  _j__  1
q1 0:@qndivl :: 1:  _j__  1j1
q1 0:@qndivl :: 1:  _j__  1j_
q1 0:@qndivl :: 1:  _j__  _j__
q1 0:@qndivl :: 1:  _j__  _j_1
q1 0:@qndivl :: 1:  _j__  _
q1 0:@qndivl :: 1:  _j__  _j1
q1 0:@qndivl :: 1:  _j__  _j_
q1 0:@qndivl :: 1:  _j_1 __j__
q1 0:@qndivl :: 1:  _j_1 __j_1
q1 0:@qndivl :: 1:  _j_1 __
q1 0:@qndivl :: 1:  _j_1 __j1
q1 0:@qndivl :: 1:  _j_1 __j_
q1 0:@qndivl :: 1:  _j_1 _1j__
q1 0:@qndivl :: 1:  _j_1 _1j_
q1 0:@qndivl :: 1:  _j_1  0j__
q1 0:@qndivl :: 1:  _j_1  0j_
q1 0:@qndivl :: 1:  _j_1  1j__
q1 0:@qndivl :: 1:  _j_1  1j_
q1 0:@qndivl :: 1:  _j_1  _j__
q1 0:@qndivl :: 1:  _j_1  _j_1
q1 0:@qndivl :: 1:  _j_1  _
q1 0:@qndivl :: 1:  _j_1  _j1
q1 0:@qndivl :: 1:  _j_1  _j_
q1 0:@qndivl :: 1:  _    __j__
q1 0:@qndivl :: 1:  _    __j_1
q1 0:@qndivl :: 1:  _    __
q1 0:@qndivl :: 1:  _    __j1
q1 0:@qndivl :: 1:  _    __j_
q1 0:@qndivl :: 1:  _    _1j__
q1 0:@qndivl :: 1:  _    _1j_
q1 0:@qndivl :: 1:  _     0j__
q1 0:@qndivl :: 1:  _     0j_
q1 0:@qndivl :: 1:  _     1j__
q1 0:@qndivl :: 1:  _     1j_
q1 0:@qndivl :: 1:  _     _j__
q1 0:@qndivl :: 1:  _     _j_1
q1 0:@qndivl :: 1:  _     _
q1 0:@qndivl :: 1:  _     _j1
q1 0:@qndivl :: 1:  _     _j_
q1 0:@qndivl :: 1:  _j1  __j__
q1 0:@qndivl :: 1:  _j1  __j_1
q1 0:@qndivl :: 1:  _j1  __
q1 0:@qndivl :: 1:  _j1  __j1
q1 0:@qndivl :: 1:  _j1  __j_
q1 0:@qndivl :: 1:  _j1  _1j__
q1 0:@qndivl :: 1:  _j1  _1j_
q1 0:@qndivl :: 1:  _j1   0j__
q1 0:@qndivl :: 1:  _j1   0j_
q1 0:@qndivl :: 1:  _j1   1j__
q1 0:@qndivl :: 1:  _j1   1j_
q1 0:@qndivl :: 1:  _j1   _j__
q1 0:@qndivl :: 1:  _j1   _j_1
q1 0:@qndivl :: 1:  _j1   _
q1 0:@qndivl :: 1:  _j1   _j1
q1 0:@qndivl :: 1:  _j1   _j_
q1 0:@qndivl :: 1:  _j_  __j__
q1 0:@qndivl :: 1:  _j_  __j_1
q1 0:@qndivl :: 1:  _j_  __
q1 0:@qndivl :: 1:  _j_  __j1
q1 0:@qndivl :: 1:  _j_  __j_
q1 0:@qndivl :: 1:  _j_  _1j__
q1 0:@qndivl :: 1:  _j_  _1j_1
q1 0:@qndivl :: 1:  _j_  _1
q1 0:@qndivl :: 1:  _j_  _1j1
q1 0:@qndivl :: 1:  _j_  _1j_
q1 0:@qndivl :: 1:  _j_   0j__
q1 0:@qndivl :: 1:  _j_   0j_1
q1 0:@qndivl :: 1:  _j_   0
q1 0:@qndivl :: 1:  _j_   0j1
q1 0:@qndivl :: 1:  _j_   0j_
q1 0:@qndivl :: 1:  _j_   1j__
q1 0:@qndivl :: 1:  _j_   1j_1
q1 0:@qndivl :: 1:  _j_   1
q1 0:@qndivl :: 1:  _j_   1j1
q1 0:@qndivl :: 1:  _j_   1j_
q1 0:@qndivl :: 1:  _j_   _j__
q1 0:@qndivl :: 1:  _j_   _j_1
q1 0:@qndivl :: 1:  _j_   _
q1 0:@qndivl :: 1:  _j_   _j1
q1 0:@qndivl :: 1:  _j_   _j_
NB. - input contains directed infinity
__j_   _j_  -: (__ qn1 q1) qndivl        q2
 _j__ __j__ -: ( _ qn1 q1) qndivl        q2
__j__ __j_  -: (__ qni q1) qndivl        q2
 _j_   _j__ -: ( _ qni q1) qndivl        q2
__j_  __j__ -: (__ qnj q1) qndivl        q2
 _j__  _j_  -: ( _ qnj q1) qndivl        q2
__j__  _j__ -: (__ qnk q1) qndivl        q2
 _j_  __j_  -: ( _ qnk q1) qndivl        q2
 0     0    -:         q1  qndivl __ qn1 q2
 0     0    -:         q1  qndivl __ qni q2
 0     0    -:         q1  qndivl __ qnj q2
 0     0    -:         q1  qndivl __ qnk q2
 0     0    -:         q1  qndivl  _ qn1 q2
 0     0    -:         q1  qndivl  _ qni q2
 0     0    -:         q1  qndivl  _ qnj q2
 0     0    -:         q1  qndivl  _ qnk q2
 0     0    -: (__ qn1 q1) qndivl __ qn1 q2
 0     0    -: (__ qn1 q1) qndivl __ qni q2
 0     0    -: (__ qn1 q1) qndivl __ qnj q2
 0     0    -: (__ qn1 q1) qndivl __ qnk q2
 0     0    -: (__ qn1 q1) qndivl  _ qn1 q2
 0     0    -: (__ qn1 q1) qndivl  _ qni q2
 0     0    -: (__ qn1 q1) qndivl  _ qnj q2
 0     0    -: (__ qn1 q1) qndivl  _ qnk q2
 0     0    -: ( _ qn1 q1) qndivl __ qn1 q2
 0     0    -: ( _ qn1 q1) qndivl __ qni q2
 0     0    -: ( _ qn1 q1) qndivl __ qnj q2
 0     0    -: ( _ qn1 q1) qndivl __ qnk q2
 0     0    -: ( _ qn1 q1) qndivl  _ qn1 q2
 0     0    -: ( _ qn1 q1) qndivl  _ qni q2
 0     0    -: ( _ qn1 q1) qndivl  _ qnj q2
 0     0    -: ( _ qn1 q1) qndivl  _ qnk q2
(5j1 0j2 ,: 0) -: q1 qndivl q3 ,: _ qn1 q4
(5j1 0j2 ,: 0) -: q1 qndivl q3 ,: _ qni q4
(5j1 0j2 ,: 0) -: q1 qndivl q3 ,: _ qnj q4
(5j1 0j2 ,: 0) -: q1 qndivl q3 ,: _ qnk q4
(5j1 0j2 ,: _j__ __j__) -: (q1 ,: _ qn1 q2) qndivl q3
(5j1 0j2 ,: _j_   _j__) -: (q1 ,: _ qni q2) qndivl q3
(5j1 0j2 ,: _j__  _j_ ) -: (q1 ,: _ qnj q2) qndivl q3
(5j1 0j2 ,: _j_  __j_ ) -: (q1 ,: _ qnk q2) qndivl q3
NB. - edge cases input
1     0j1    -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivl (FP_OVFL  j.  FP_OVFL) , 0
1     0j_1   -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivl  0                     , FP_OVFL j. FP_OVFL
0.5   0j_0.5 -: ((FP_OVFL j. FP_OVFL) , 0                 ) qndivl (FP_OVFL  j.  FP_OVFL) , FP_OVFL j. FP_OVFL
0.5   0j0.5  -: ( 0                   , FP_OVFL j. FP_OVFL) qndivl (FP_OVFL  j.  FP_OVFL) , FP_OVFL j. FP_OVFL
_j_   _j_    -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivl  FP_UNFL               , 0
_j__ __j_    -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivl (0        j. FP_UNFL)  , 0
_j_  __j__   -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivl  0                     , FP_UNFL j.   0
_j__  _j__   -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivl  0                     , 0       j. FP_UNFL
0     0      -: ((FP_UNFL j. FP_UNFL) , FP_UNFL j. FP_UNFL) qndivl (FP_OVFL  j.  FP_OVFL) , FP_OVFL j. FP_OVFL
0     0      -: ((FP_UNFL j. FP_UNFL) , FP_UNFL j. FP_UNFL) qndivl (FP_OVFL  j.  FP_OVFL) , 0
0     0      -: ((FP_UNFL j. FP_UNFL) , FP_UNFL j. FP_UNFL) qndivl  0                     , FP_OVFL j. FP_OVFL
0     0      -: ((FP_UNFL j. FP_UNFL) , 0                 ) qndivl (FP_OVFL  j.  FP_OVFL) , FP_OVFL j. FP_OVFL
0     0      -: ( 0                   , FP_UNFL j. FP_UNFL) qndivl (FP_OVFL  j.  FP_OVFL) , FP_OVFL j. FP_OVFL
1     0j1    -: ((FP_UNFL j. FP_UNFL) , FP_UNFL j. FP_UNFL) qndivl (FP_UNFL  j.  FP_UNFL) , 0
1     0j_1   -: ((FP_UNFL j. FP_UNFL) , FP_UNFL j. FP_UNFL) qndivl  0                     , FP_UNFL j. FP_UNFL
0.5   0j_0.5 -: ((FP_UNFL j. FP_UNFL) , 0                 ) qndivl (FP_UNFL  j.  FP_UNFL) , FP_UNFL j. FP_UNFL
0.5   0j0.5  -: ( 0                   , FP_UNFL j. FP_UNFL) qndivl (FP_UNFL  j.  FP_UNFL) , FP_UNFL j. FP_UNFL
_     0            -: q5 qndivl q7
(0 ,~ FP_UNFL % 8) -: q7 qndivl q5
NB. - input without edge cases
 0 0 -:   0 0 qndivl  q1
 1 0 -:   q1  qndivl  q1
 1 0 -:   q2  qndivl  q2
_1 0 -:   q1  qndivl -q1
_1 0 -: (-q1) qndivl  q1
_1 0 -:   q2  qndivl -q2
_1 0 -: (-q2) qndivl  q2
q1 -: q1 qndivl 1 0
q2 -: _60j20 14j32 qndivl q1
q1 -: _60j12 30j24 qndivl q2
(1 0 ,: q1) -:  q1        qndivl q1  ,: 1 0
(q1  ,: q2) -: (q1 ,: q2) qndivl 1 0
(1 0 ,: q2) -: (q1 ,: q2) qndivl q1  ,: 1 0
(-: qndivl~) i. 0 2

NB. qndivr
NB. - input contains NaN
1 1 -: isnan  0     0j_. qndivr  0     0j_.
1 1 -: isnan  0     0j_. qndivr  0    _.
1 1 -: isnan  0     0j_. qndivr  0j_.  0
1 1 -: isnan  0     0j_. qndivr _.     0
1 1 -: isnan  0    _.    qndivr  0     0j_.
1 1 -: isnan  0    _.    qndivr  0    _.
1 1 -: isnan  0    _.    qndivr  0j_.  0
1 1 -: isnan  0    _.    qndivr _.     0
1 1 -: isnan  0j_.  0    qndivr  0     0j_.
1 1 -: isnan  0j_.  0    qndivr  0    _.
1 1 -: isnan  0j_.  0    qndivr  0j_.  0
1 1 -: isnan  0j_.  0    qndivr _.     0
1 1 -: isnan _.     0    qndivr  0     0j_.
1 1 -: isnan _.     0    qndivr  0    _.
1 1 -: isnan _.     0    qndivr  0j_.  0
1 1 -: isnan _.     0    qndivr _.     0
1 1 -: isnan  0    _.j_. qndivr  0    _.j_.
1 1 -: isnan  0    _.j_. qndivr  0j_.  0j_.
1 1 -: isnan  0    _.j_. qndivr  0j_. _.
1 1 -: isnan  0    _.j_. qndivr _.     0j_.
1 1 -: isnan  0    _.j_. qndivr _.    _.
1 1 -: isnan  0    _.j_. qndivr _.j_.  0
1 1 -: isnan  0j_.  0j_. qndivr  0    _.j_.
1 1 -: isnan  0j_.  0j_. qndivr  0j_.  0j_.
1 1 -: isnan  0j_.  0j_. qndivr  0j_. _.
1 1 -: isnan  0j_.  0j_. qndivr _.     0j_.
1 1 -: isnan  0j_.  0j_. qndivr _.    _.
1 1 -: isnan  0j_.  0j_. qndivr _.j_.  0
1 1 -: isnan  0j_. _.    qndivr  0    _.j_.
1 1 -: isnan  0j_. _.    qndivr  0j_.  0j_.
1 1 -: isnan  0j_. _.    qndivr  0j_. _.
1 1 -: isnan  0j_. _.    qndivr _.     0j_.
1 1 -: isnan  0j_. _.    qndivr _.    _.
1 1 -: isnan  0j_. _.    qndivr _.j_.  0
1 1 -: isnan _.     0j_. qndivr  0    _.j_.
1 1 -: isnan _.     0j_. qndivr  0j_.  0j_.
1 1 -: isnan _.     0j_. qndivr  0j_. _.
1 1 -: isnan _.     0j_. qndivr _.     0j_.
1 1 -: isnan _.     0j_. qndivr _.    _.
1 1 -: isnan _.     0j_. qndivr _.j_.  0
1 1 -: isnan _.    _.    qndivr  0    _.j_.
1 1 -: isnan _.    _.    qndivr  0j_.  0j_.
1 1 -: isnan _.    _.    qndivr  0j_. _.
1 1 -: isnan _.    _.    qndivr _.     0j_.
1 1 -: isnan _.    _.    qndivr _.    _.
1 1 -: isnan _.    _.    qndivr _.j_.  0
1 1 -: isnan _.j_.  0    qndivr  0    _.j_.
1 1 -: isnan _.j_.  0    qndivr  0j_.  0j_.
1 1 -: isnan _.j_.  0    qndivr  0j_. _.
1 1 -: isnan _.j_.  0    qndivr _.     0j_.
1 1 -: isnan _.j_.  0    qndivr _.    _.
1 1 -: isnan _.j_.  0    qndivr _.j_.  0
1 1 -: isnan _.j_. _.    qndivr _.j_. _.
1 1 -: isnan _.j_. _.    qndivr _.j_.  0j_.
1 1 -: isnan _.j_. _.    qndivr _.    _.j_.
1 1 -: isnan _.j_. _.    qndivr  0j_. _.j_.
1 1 -: isnan _.j_.  0j_. qndivr _.j_. _.
1 1 -: isnan _.j_.  0j_. qndivr _.j_.  0j_.
1 1 -: isnan _.j_.  0j_. qndivr _.    _.j_.
1 1 -: isnan _.j_.  0j_. qndivr  0j_. _.j_.
1 1 -: isnan _.    _.j_. qndivr _.j_. _.
1 1 -: isnan _.    _.j_. qndivr _.j_.  0j_.
1 1 -: isnan _.    _.j_. qndivr _.    _.j_.
1 1 -: isnan _.    _.j_. qndivr  0j_. _.j_.
1 1 -: isnan  0j_. _.j_. qndivr _.j_. _.
1 1 -: isnan  0j_. _.j_. qndivr _.j_.  0j_.
1 1 -: isnan  0j_. _.j_. qndivr _.    _.j_.
1 1 -: isnan  0j_. _.j_. qndivr  0j_. _.j_.
1 1 -: isnan _.j_. _.j_. qndivr _.j_. _.j_.
1 1 -: isnan         q1  qndivr _.    _.
1 1 -: isnan         q1  qndivr _.j_. _.j_.
1 1 -: isnan _.    _.    qndivr        q2
1 1 -: isnan _.j_. _.j_. qndivr        q2
1 1 -: isnan (_. qn1 q1) qndivr        q2
1 1 -: isnan (_. qni q1) qndivr        q2
1 1 -: isnan (_. qnj q1) qndivr        q2
1 1 -: isnan (_. qnk q1) qndivr        q2
1 1 -: isnan         q1  qndivr _. qn1 q2
1 1 -: isnan         q1  qndivr _. qni q2
1 1 -: isnan         q1  qndivr _. qnj q2
1 1 -: isnan         q1  qndivr _. qnk q2
(,.~ 0 1) -: isnan  q1               qndivr q3 ,: _. qn1 q4
(,.~ 0 1) -: isnan  q1               qndivr q3 ,: _. qni q4
(,.~ 0 1) -: isnan  q1               qndivr q3 ,: _. qnj q4
(,.~ 0 1) -: isnan  q1               qndivr q3 ,: _. qnk q4
(,.~ 0 1) -: isnan (q1 ,: _. qn1 q2) qndivr q3
(,.~ 0 1) -: isnan (q1 ,: _. qni q2) qndivr q3
(,.~ 0 1) -: isnan (q1 ,: _. qnj q2) qndivr q3
(,.~ 0 1) -: isnan (q1 ,: _. qnk q2) qndivr q3
NB. - denominator is zero so NaN error must be throwed
q1 0:@qndivr :: 1: 0 0
NB. - input contains q∞ so NaN error must be throwed
__j__ __j__ 0:@qndivr :: 1 q1
__j__ __j_1 0:@qndivr :: 1 q1
__j__ __    0:@qndivr :: 1 q1
__j__ __j1  0:@qndivr :: 1 q1
__j__ __j_  0:@qndivr :: 1 q1
__j__ _1j__ 0:@qndivr :: 1 q1
__j__ _1j_1 0:@qndivr :: 1 q1
__j__ _1    0:@qndivr :: 1 q1
__j__ _1j1  0:@qndivr :: 1 q1
__j__ _1j_  0:@qndivr :: 1 q1
__j__  0j__ 0:@qndivr :: 1 q1
__j__  0j_1 0:@qndivr :: 1 q1
__j__  0    0:@qndivr :: 1 q1
__j__  0j1  0:@qndivr :: 1 q1
__j__  0j_  0:@qndivr :: 1 q1
__j__  1j__ 0:@qndivr :: 1 q1
__j__  1j_1 0:@qndivr :: 1 q1
__j__  1    0:@qndivr :: 1 q1
__j__  1j1  0:@qndivr :: 1 q1
__j__  1j_  0:@qndivr :: 1 q1
__j__  _j__ 0:@qndivr :: 1 q1
__j__  _j_1 0:@qndivr :: 1 q1
__j__  _    0:@qndivr :: 1 q1
__j__  _j1  0:@qndivr :: 1 q1
__j__  _j_  0:@qndivr :: 1 q1
__j_1 __j__ 0:@qndivr :: 1 q1
__j_1 __j_1 0:@qndivr :: 1 q1
__j_1 __    0:@qndivr :: 1 q1
__j_1 __j1  0:@qndivr :: 1 q1
__j_1 __j_  0:@qndivr :: 1 q1
__j_1 _1j__ 0:@qndivr :: 1 q1
__j_1 _1j_  0:@qndivr :: 1 q1
__j_1  0j__ 0:@qndivr :: 1 q1
__j_1  0j_  0:@qndivr :: 1 q1
__j_1  1j__ 0:@qndivr :: 1 q1
__j_1  1j_  0:@qndivr :: 1 q1
__j_1  _j__ 0:@qndivr :: 1 q1
__j_1  _j_1 0:@qndivr :: 1 q1
__j_1  _    0:@qndivr :: 1 q1
__j_1  _j1  0:@qndivr :: 1 q1
__j_1  _j_  0:@qndivr :: 1 q1
__    __j__ 0:@qndivr :: 1 q1
__    __j_1 0:@qndivr :: 1 q1
__    __    0:@qndivr :: 1 q1
__    __j1  0:@qndivr :: 1 q1
__    __j_  0:@qndivr :: 1 q1
__    _1j__ 0:@qndivr :: 1 q1
__    _1j_  0:@qndivr :: 1 q1
__     0j__ 0:@qndivr :: 1 q1
__     0j_  0:@qndivr :: 1 q1
__     1j__ 0:@qndivr :: 1 q1
__     1j_  0:@qndivr :: 1 q1
__     _j__ 0:@qndivr :: 1 q1
__     _j_1 0:@qndivr :: 1 q1
__     _    0:@qndivr :: 1 q1
__     _j1  0:@qndivr :: 1 q1
__     _j_  0:@qndivr :: 1 q1
__j1  __j__ 0:@qndivr :: 1 q1
__j1  __j_1 0:@qndivr :: 1 q1
__j1  __    0:@qndivr :: 1 q1
__j1  __j1  0:@qndivr :: 1 q1
__j1  __j_  0:@qndivr :: 1 q1
__j1  _1j__ 0:@qndivr :: 1 q1
__j1  _1j_  0:@qndivr :: 1 q1
__j1   0j__ 0:@qndivr :: 1 q1
__j1   0j_  0:@qndivr :: 1 q1
__j1   1j__ 0:@qndivr :: 1 q1
__j1   1j_  0:@qndivr :: 1 q1
__j1   _j__ 0:@qndivr :: 1 q1
__j1   _j_1 0:@qndivr :: 1 q1
__j1   _    0:@qndivr :: 1 q1
__j1   _j1  0:@qndivr :: 1 q1
__j1   _j_  0:@qndivr :: 1 q1
__j_  __j__ 0:@qndivr :: 1 q1
__j_  __j_1 0:@qndivr :: 1 q1
__j_  __    0:@qndivr :: 1 q1
__j_  __j1  0:@qndivr :: 1 q1
__j_  __j_  0:@qndivr :: 1 q1
__j_  _1j__ 0:@qndivr :: 1 q1
__j_  _1j_1 0:@qndivr :: 1 q1
__j_  _1    0:@qndivr :: 1 q1
__j_  _1j1  0:@qndivr :: 1 q1
__j_  _1j_  0:@qndivr :: 1 q1
__j_   0j__ 0:@qndivr :: 1 q1
__j_   0j_1 0:@qndivr :: 1 q1
__j_   0    0:@qndivr :: 1 q1
__j_   0j1  0:@qndivr :: 1 q1
__j_   0j_  0:@qndivr :: 1 q1
__j_   1j__ 0:@qndivr :: 1 q1
__j_   1j_1 0:@qndivr :: 1 q1
__j_   1    0:@qndivr :: 1 q1
__j_   1j1  0:@qndivr :: 1 q1
__j_   1j_  0:@qndivr :: 1 q1
__j_   _j__ 0:@qndivr :: 1 q1
__j_   _j_1 0:@qndivr :: 1 q1
__j_   _    0:@qndivr :: 1 q1
__j_   _j1  0:@qndivr :: 1 q1
__j_   _j_  0:@qndivr :: 1 q1
_1j__ __j__ 0:@qndivr :: 1 q1
_1j__ __j_1 0:@qndivr :: 1 q1
_1j__ __    0:@qndivr :: 1 q1
_1j__ __j1  0:@qndivr :: 1 q1
_1j__ __j_  0:@qndivr :: 1 q1
_1j__ _1j__ 0:@qndivr :: 1 q1
_1j__ _1j_  0:@qndivr :: 1 q1
_1j__  0j__ 0:@qndivr :: 1 q1
_1j__  0j_  0:@qndivr :: 1 q1
_1j__  1j__ 0:@qndivr :: 1 q1
_1j__  1j_  0:@qndivr :: 1 q1
_1j__  _j__ 0:@qndivr :: 1 q1
_1j__  _j_1 0:@qndivr :: 1 q1
_1j__  _    0:@qndivr :: 1 q1
_1j__  _j1  0:@qndivr :: 1 q1
_1j__  _j_  0:@qndivr :: 1 q1
_1j_1 __j__ 0:@qndivr :: 1 q1
_1j_1 __j_  0:@qndivr :: 1 q1
_1j_1  _j__ 0:@qndivr :: 1 q1
_1j_1  _j_  0:@qndivr :: 1 q1
_1    __j__ 0:@qndivr :: 1 q1
_1    __j_  0:@qndivr :: 1 q1
_1     _j__ 0:@qndivr :: 1 q1
_1     _j_  0:@qndivr :: 1 q1
_1j1  __j__ 0:@qndivr :: 1 q1
_1j1  __j_  0:@qndivr :: 1 q1
_1j1   _j__ 0:@qndivr :: 1 q1
_1j1   _j_  0:@qndivr :: 1 q1
_1j_  __j__ 0:@qndivr :: 1 q1
_1j_  __j_1 0:@qndivr :: 1 q1
_1j_  __    0:@qndivr :: 1 q1
_1j_  __j1  0:@qndivr :: 1 q1
_1j_  __j_  0:@qndivr :: 1 q1
_1j_  _1j__ 0:@qndivr :: 1 q1
_1j_  _1j_  0:@qndivr :: 1 q1
_1j_   0j__ 0:@qndivr :: 1 q1
_1j_   0j_  0:@qndivr :: 1 q1
_1j_   1j__ 0:@qndivr :: 1 q1
_1j_   1j_  0:@qndivr :: 1 q1
_1j_   _j__ 0:@qndivr :: 1 q1
_1j_   _j_1 0:@qndivr :: 1 q1
_1j_   _    0:@qndivr :: 1 q1
_1j_   _j1  0:@qndivr :: 1 q1
_1j_   _j_  0:@qndivr :: 1 q1
 0j__ __j__ 0:@qndivr :: 1 q1
 0j__ __j_1 0:@qndivr :: 1 q1
 0j__ __    0:@qndivr :: 1 q1
 0j__ __j1  0:@qndivr :: 1 q1
 0j__ __j_  0:@qndivr :: 1 q1
 0j__ _1j__ 0:@qndivr :: 1 q1
 0j__ _1j_  0:@qndivr :: 1 q1
 0j__  0j__ 0:@qndivr :: 1 q1
 0j__  0j_  0:@qndivr :: 1 q1
 0j__  1j__ 0:@qndivr :: 1 q1
 0j__  1j_  0:@qndivr :: 1 q1
 0j__  _j__ 0:@qndivr :: 1 q1
 0j__  _j_1 0:@qndivr :: 1 q1
 0j__  _    0:@qndivr :: 1 q1
 0j__  _j1  0:@qndivr :: 1 q1
 0j__  _j_  0:@qndivr :: 1 q1
 0j_1 __j__ 0:@qndivr :: 1 q1
 0j_1 __j_  0:@qndivr :: 1 q1
 0j_1  _j__ 0:@qndivr :: 1 q1
 0j_1  _j_  0:@qndivr :: 1 q1
 0    __j__ 0:@qndivr :: 1 q1
 0    __j_  0:@qndivr :: 1 q1
 0     _j__ 0:@qndivr :: 1 q1
 0     _j_  0:@qndivr :: 1 q1
 0j1  __j__ 0:@qndivr :: 1 q1
 0j1  __j_  0:@qndivr :: 1 q1
 0j1   _j__ 0:@qndivr :: 1 q1
 0j1   _j_  0:@qndivr :: 1 q1
 0j_  __j__ 0:@qndivr :: 1 q1
 0j_  __j_1 0:@qndivr :: 1 q1
 0j_  __    0:@qndivr :: 1 q1
 0j_  __j1  0:@qndivr :: 1 q1
 0j_  __j_  0:@qndivr :: 1 q1
 0j_  _1j__ 0:@qndivr :: 1 q1
 0j_  _1j_  0:@qndivr :: 1 q1
 0j_   0j__ 0:@qndivr :: 1 q1
 0j_   0j_  0:@qndivr :: 1 q1
 0j_   1j__ 0:@qndivr :: 1 q1
 0j_   1j_  0:@qndivr :: 1 q1
 0j_   _j__ 0:@qndivr :: 1 q1
 0j_   _j_1 0:@qndivr :: 1 q1
 0j_   _    0:@qndivr :: 1 q1
 0j_   _j1  0:@qndivr :: 1 q1
 0j_   _j_  0:@qndivr :: 1 q1
 1j__ __j__ 0:@qndivr :: 1 q1
 1j__ __j_1 0:@qndivr :: 1 q1
 1j__ __    0:@qndivr :: 1 q1
 1j__ __j1  0:@qndivr :: 1 q1
 1j__ __j_  0:@qndivr :: 1 q1
 1j__ _1j__ 0:@qndivr :: 1 q1
 1j__ _1j_  0:@qndivr :: 1 q1
 1j__  0j__ 0:@qndivr :: 1 q1
 1j__  0j_  0:@qndivr :: 1 q1
 1j__  1j__ 0:@qndivr :: 1 q1
 1j__  1j_  0:@qndivr :: 1 q1
 1j__  _j__ 0:@qndivr :: 1 q1
 1j__  _j_1 0:@qndivr :: 1 q1
 1j__  _    0:@qndivr :: 1 q1
 1j__  _j1  0:@qndivr :: 1 q1
 1j__  _j_  0:@qndivr :: 1 q1
 1j_1 __j__ 0:@qndivr :: 1 q1
 1j_1 __j_  0:@qndivr :: 1 q1
 1j_1  _j__ 0:@qndivr :: 1 q1
 1j_1  _j_  0:@qndivr :: 1 q1
 1    __j__ 0:@qndivr :: 1 q1
 1    __j_  0:@qndivr :: 1 q1
 1     _j__ 0:@qndivr :: 1 q1
 1     _j_  0:@qndivr :: 1 q1
 1j1  __j__ 0:@qndivr :: 1 q1
 1j1  __j_  0:@qndivr :: 1 q1
 1j1   _j__ 0:@qndivr :: 1 q1
 1j1   _j_  0:@qndivr :: 1 q1
 1j_  __j__ 0:@qndivr :: 1 q1
 1j_  __j_1 0:@qndivr :: 1 q1
 1j_  __    0:@qndivr :: 1 q1
 1j_  __j1  0:@qndivr :: 1 q1
 1j_  __j_  0:@qndivr :: 1 q1
 1j_  _1j__ 0:@qndivr :: 1 q1
 1j_  _1j_  0:@qndivr :: 1 q1
 1j_   0j__ 0:@qndivr :: 1 q1
 1j_   0j_  0:@qndivr :: 1 q1
 1j_   1j__ 0:@qndivr :: 1 q1
 1j_   1j_  0:@qndivr :: 1 q1
 1j_   _j__ 0:@qndivr :: 1 q1
 1j_   _j_1 0:@qndivr :: 1 q1
 1j_   _    0:@qndivr :: 1 q1
 1j_   _j1  0:@qndivr :: 1 q1
 1j_   _j_  0:@qndivr :: 1 q1
 _j__ __j__ 0:@qndivr :: 1 q1
 _j__ __j_1 0:@qndivr :: 1 q1
 _j__ __    0:@qndivr :: 1 q1
 _j__ __j1  0:@qndivr :: 1 q1
 _j__ __j_  0:@qndivr :: 1 q1
 _j__ _1j__ 0:@qndivr :: 1 q1
 _j__ _1j_1 0:@qndivr :: 1 q1
 _j__ _1    0:@qndivr :: 1 q1
 _j__ _1j1  0:@qndivr :: 1 q1
 _j__ _1j_  0:@qndivr :: 1 q1
 _j__  0j__ 0:@qndivr :: 1 q1
 _j__  0j_1 0:@qndivr :: 1 q1
 _j__  0    0:@qndivr :: 1 q1
 _j__  0j1  0:@qndivr :: 1 q1
 _j__  0j_  0:@qndivr :: 1 q1
 _j__  1j__ 0:@qndivr :: 1 q1
 _j__  1j_1 0:@qndivr :: 1 q1
 _j__  1    0:@qndivr :: 1 q1
 _j__  1j1  0:@qndivr :: 1 q1
 _j__  1j_  0:@qndivr :: 1 q1
 _j__  _j__ 0:@qndivr :: 1 q1
 _j__  _j_1 0:@qndivr :: 1 q1
 _j__  _    0:@qndivr :: 1 q1
 _j__  _j1  0:@qndivr :: 1 q1
 _j__  _j_  0:@qndivr :: 1 q1
 _j_1 __j__ 0:@qndivr :: 1 q1
 _j_1 __j_1 0:@qndivr :: 1 q1
 _j_1 __    0:@qndivr :: 1 q1
 _j_1 __j1  0:@qndivr :: 1 q1
 _j_1 __j_  0:@qndivr :: 1 q1
 _j_1 _1j__ 0:@qndivr :: 1 q1
 _j_1 _1j_  0:@qndivr :: 1 q1
 _j_1  0j__ 0:@qndivr :: 1 q1
 _j_1  0j_  0:@qndivr :: 1 q1
 _j_1  1j__ 0:@qndivr :: 1 q1
 _j_1  1j_  0:@qndivr :: 1 q1
 _j_1  _j__ 0:@qndivr :: 1 q1
 _j_1  _j_1 0:@qndivr :: 1 q1
 _j_1  _    0:@qndivr :: 1 q1
 _j_1  _j1  0:@qndivr :: 1 q1
 _j_1  _j_  0:@qndivr :: 1 q1
 _    __j__ 0:@qndivr :: 1 q1
 _    __j_1 0:@qndivr :: 1 q1
 _    __    0:@qndivr :: 1 q1
 _    __j1  0:@qndivr :: 1 q1
 _    __j_  0:@qndivr :: 1 q1
 _    _1j__ 0:@qndivr :: 1 q1
 _    _1j_  0:@qndivr :: 1 q1
 _     0j__ 0:@qndivr :: 1 q1
 _     0j_  0:@qndivr :: 1 q1
 _     1j__ 0:@qndivr :: 1 q1
 _     1j_  0:@qndivr :: 1 q1
 _     _j__ 0:@qndivr :: 1 q1
 _     _j_1 0:@qndivr :: 1 q1
 _     _    0:@qndivr :: 1 q1
 _     _j1  0:@qndivr :: 1 q1
 _     _j_  0:@qndivr :: 1 q1
 _j1  __j__ 0:@qndivr :: 1 q1
 _j1  __j_1 0:@qndivr :: 1 q1
 _j1  __    0:@qndivr :: 1 q1
 _j1  __j1  0:@qndivr :: 1 q1
 _j1  __j_  0:@qndivr :: 1 q1
 _j1  _1j__ 0:@qndivr :: 1 q1
 _j1  _1j_  0:@qndivr :: 1 q1
 _j1   0j__ 0:@qndivr :: 1 q1
 _j1   0j_  0:@qndivr :: 1 q1
 _j1   1j__ 0:@qndivr :: 1 q1
 _j1   1j_  0:@qndivr :: 1 q1
 _j1   _j__ 0:@qndivr :: 1 q1
 _j1   _j_1 0:@qndivr :: 1 q1
 _j1   _    0:@qndivr :: 1 q1
 _j1   _j1  0:@qndivr :: 1 q1
 _j1   _j_  0:@qndivr :: 1 q1
 _j_  __j__ 0:@qndivr :: 1 q1
 _j_  __j_1 0:@qndivr :: 1 q1
 _j_  __    0:@qndivr :: 1 q1
 _j_  __j1  0:@qndivr :: 1 q1
 _j_  __j_  0:@qndivr :: 1 q1
 _j_  _1j__ 0:@qndivr :: 1 q1
 _j_  _1j_1 0:@qndivr :: 1 q1
 _j_  _1    0:@qndivr :: 1 q1
 _j_  _1j1  0:@qndivr :: 1 q1
 _j_  _1j_  0:@qndivr :: 1 q1
 _j_   0j__ 0:@qndivr :: 1 q1
 _j_   0j_1 0:@qndivr :: 1 q1
 _j_   0    0:@qndivr :: 1 q1
 _j_   0j1  0:@qndivr :: 1 q1
 _j_   0j_  0:@qndivr :: 1 q1
 _j_   1j__ 0:@qndivr :: 1 q1
 _j_   1j_1 0:@qndivr :: 1 q1
 _j_   1    0:@qndivr :: 1 q1
 _j_   1j1  0:@qndivr :: 1 q1
 _j_   1j_  0:@qndivr :: 1 q1
 _j_   _j__ 0:@qndivr :: 1 q1
 _j_   _j_1 0:@qndivr :: 1 q1
 _j_   _    0:@qndivr :: 1 q1
 _j_   _j1  0:@qndivr :: 1 q1
 _j_   _j_  0:@qndivr :: 1 q1
q1 0:@qndivr :: 1: __j__ __j__
q1 0:@qndivr :: 1: __j__ __j_1
q1 0:@qndivr :: 1: __j__ __
q1 0:@qndivr :: 1: __j__ __j1
q1 0:@qndivr :: 1: __j__ __j_
q1 0:@qndivr :: 1: __j__ _1j__
q1 0:@qndivr :: 1: __j__ _1j_1
q1 0:@qndivr :: 1: __j__ _1
q1 0:@qndivr :: 1: __j__ _1j1
q1 0:@qndivr :: 1: __j__ _1j_
q1 0:@qndivr :: 1: __j__  0j__
q1 0:@qndivr :: 1: __j__  0j_1
q1 0:@qndivr :: 1: __j__  0
q1 0:@qndivr :: 1: __j__  0j1
q1 0:@qndivr :: 1: __j__  0j_
q1 0:@qndivr :: 1: __j__  1j__
q1 0:@qndivr :: 1: __j__  1j_1
q1 0:@qndivr :: 1: __j__  1
q1 0:@qndivr :: 1: __j__  1j1
q1 0:@qndivr :: 1: __j__  1j_
q1 0:@qndivr :: 1: __j__  _j__
q1 0:@qndivr :: 1: __j__  _j_1
q1 0:@qndivr :: 1: __j__  _
q1 0:@qndivr :: 1: __j__  _j1
q1 0:@qndivr :: 1: __j__  _j_
q1 0:@qndivr :: 1: __j_1 __j__
q1 0:@qndivr :: 1: __j_1 __j_1
q1 0:@qndivr :: 1: __j_1 __
q1 0:@qndivr :: 1: __j_1 __j1
q1 0:@qndivr :: 1: __j_1 __j_
q1 0:@qndivr :: 1: __j_1 _1j__
q1 0:@qndivr :: 1: __j_1 _1j_
q1 0:@qndivr :: 1: __j_1  0j__
q1 0:@qndivr :: 1: __j_1  0j_
q1 0:@qndivr :: 1: __j_1  1j__
q1 0:@qndivr :: 1: __j_1  1j_
q1 0:@qndivr :: 1: __j_1  _j__
q1 0:@qndivr :: 1: __j_1  _j_1
q1 0:@qndivr :: 1: __j_1  _
q1 0:@qndivr :: 1: __j_1  _j1
q1 0:@qndivr :: 1: __j_1  _j_
q1 0:@qndivr :: 1: __    __j__
q1 0:@qndivr :: 1: __    __j_1
q1 0:@qndivr :: 1: __    __
q1 0:@qndivr :: 1: __    __j1
q1 0:@qndivr :: 1: __    __j_
q1 0:@qndivr :: 1: __    _1j__
q1 0:@qndivr :: 1: __    _1j_
q1 0:@qndivr :: 1: __     0j__
q1 0:@qndivr :: 1: __     0j_
q1 0:@qndivr :: 1: __     1j__
q1 0:@qndivr :: 1: __     1j_
q1 0:@qndivr :: 1: __     _j__
q1 0:@qndivr :: 1: __     _j_1
q1 0:@qndivr :: 1: __     _
q1 0:@qndivr :: 1: __     _j1
q1 0:@qndivr :: 1: __     _j_
q1 0:@qndivr :: 1: __j1  __j__
q1 0:@qndivr :: 1: __j1  __j_1
q1 0:@qndivr :: 1: __j1  __
q1 0:@qndivr :: 1: __j1  __j1
q1 0:@qndivr :: 1: __j1  __j_
q1 0:@qndivr :: 1: __j1  _1j__
q1 0:@qndivr :: 1: __j1  _1j_
q1 0:@qndivr :: 1: __j1   0j__
q1 0:@qndivr :: 1: __j1   0j_
q1 0:@qndivr :: 1: __j1   1j__
q1 0:@qndivr :: 1: __j1   1j_
q1 0:@qndivr :: 1: __j1   _j__
q1 0:@qndivr :: 1: __j1   _j_1
q1 0:@qndivr :: 1: __j1   _
q1 0:@qndivr :: 1: __j1   _j1
q1 0:@qndivr :: 1: __j1   _j_
q1 0:@qndivr :: 1: __j_  __j__
q1 0:@qndivr :: 1: __j_  __j_1
q1 0:@qndivr :: 1: __j_  __
q1 0:@qndivr :: 1: __j_  __j1
q1 0:@qndivr :: 1: __j_  __j_
q1 0:@qndivr :: 1: __j_  _1j__
q1 0:@qndivr :: 1: __j_  _1j_1
q1 0:@qndivr :: 1: __j_  _1
q1 0:@qndivr :: 1: __j_  _1j1
q1 0:@qndivr :: 1: __j_  _1j_
q1 0:@qndivr :: 1: __j_   0j__
q1 0:@qndivr :: 1: __j_   0j_1
q1 0:@qndivr :: 1: __j_   0
q1 0:@qndivr :: 1: __j_   0j1
q1 0:@qndivr :: 1: __j_   0j_
q1 0:@qndivr :: 1: __j_   1j__
q1 0:@qndivr :: 1: __j_   1j_1
q1 0:@qndivr :: 1: __j_   1
q1 0:@qndivr :: 1: __j_   1j1
q1 0:@qndivr :: 1: __j_   1j_
q1 0:@qndivr :: 1: __j_   _j__
q1 0:@qndivr :: 1: __j_   _j_1
q1 0:@qndivr :: 1: __j_   _
q1 0:@qndivr :: 1: __j_   _j1
q1 0:@qndivr :: 1: __j_   _j_
q1 0:@qndivr :: 1: _1j__ __j__
q1 0:@qndivr :: 1: _1j__ __j_1
q1 0:@qndivr :: 1: _1j__ __
q1 0:@qndivr :: 1: _1j__ __j1
q1 0:@qndivr :: 1: _1j__ __j_
q1 0:@qndivr :: 1: _1j__ _1j__
q1 0:@qndivr :: 1: _1j__ _1j_
q1 0:@qndivr :: 1: _1j__  0j__
q1 0:@qndivr :: 1: _1j__  0j_
q1 0:@qndivr :: 1: _1j__  1j__
q1 0:@qndivr :: 1: _1j__  1j_
q1 0:@qndivr :: 1: _1j__  _j__
q1 0:@qndivr :: 1: _1j__  _j_1
q1 0:@qndivr :: 1: _1j__  _
q1 0:@qndivr :: 1: _1j__  _j1
q1 0:@qndivr :: 1: _1j__  _j_
q1 0:@qndivr :: 1: _1j_1 __j__
q1 0:@qndivr :: 1: _1j_1 __j_
q1 0:@qndivr :: 1: _1j_1  _j__
q1 0:@qndivr :: 1: _1j_1  _j_
q1 0:@qndivr :: 1: _1    __j__
q1 0:@qndivr :: 1: _1    __j_
q1 0:@qndivr :: 1: _1     _j__
q1 0:@qndivr :: 1: _1     _j_
q1 0:@qndivr :: 1: _1j1  __j__
q1 0:@qndivr :: 1: _1j1  __j_
q1 0:@qndivr :: 1: _1j1   _j__
q1 0:@qndivr :: 1: _1j1   _j_
q1 0:@qndivr :: 1: _1j_  __j__
q1 0:@qndivr :: 1: _1j_  __j_1
q1 0:@qndivr :: 1: _1j_  __
q1 0:@qndivr :: 1: _1j_  __j1
q1 0:@qndivr :: 1: _1j_  __j_
q1 0:@qndivr :: 1: _1j_  _1j__
q1 0:@qndivr :: 1: _1j_  _1j_
q1 0:@qndivr :: 1: _1j_   0j__
q1 0:@qndivr :: 1: _1j_   0j_
q1 0:@qndivr :: 1: _1j_   1j__
q1 0:@qndivr :: 1: _1j_   1j_
q1 0:@qndivr :: 1: _1j_   _j__
q1 0:@qndivr :: 1: _1j_   _j_1
q1 0:@qndivr :: 1: _1j_   _
q1 0:@qndivr :: 1: _1j_   _j1
q1 0:@qndivr :: 1: _1j_   _j_
q1 0:@qndivr :: 1:  0j__ __j__
q1 0:@qndivr :: 1:  0j__ __j_1
q1 0:@qndivr :: 1:  0j__ __
q1 0:@qndivr :: 1:  0j__ __j1
q1 0:@qndivr :: 1:  0j__ __j_
q1 0:@qndivr :: 1:  0j__ _1j__
q1 0:@qndivr :: 1:  0j__ _1j_
q1 0:@qndivr :: 1:  0j__  0j__
q1 0:@qndivr :: 1:  0j__  0j_
q1 0:@qndivr :: 1:  0j__  1j__
q1 0:@qndivr :: 1:  0j__  1j_
q1 0:@qndivr :: 1:  0j__  _j__
q1 0:@qndivr :: 1:  0j__  _j_1
q1 0:@qndivr :: 1:  0j__  _
q1 0:@qndivr :: 1:  0j__  _j1
q1 0:@qndivr :: 1:  0j__  _j_
q1 0:@qndivr :: 1:  0j_1 __j__
q1 0:@qndivr :: 1:  0j_1 __j_
q1 0:@qndivr :: 1:  0j_1  _j__
q1 0:@qndivr :: 1:  0j_1  _j_
q1 0:@qndivr :: 1:  0    __j__
q1 0:@qndivr :: 1:  0    __j_
q1 0:@qndivr :: 1:  0     _j__
q1 0:@qndivr :: 1:  0     _j_
q1 0:@qndivr :: 1:  0j1  __j__
q1 0:@qndivr :: 1:  0j1  __j_
q1 0:@qndivr :: 1:  0j1   _j__
q1 0:@qndivr :: 1:  0j1   _j_
q1 0:@qndivr :: 1:  0j_  __j__
q1 0:@qndivr :: 1:  0j_  __j_1
q1 0:@qndivr :: 1:  0j_  __
q1 0:@qndivr :: 1:  0j_  __j1
q1 0:@qndivr :: 1:  0j_  __j_
q1 0:@qndivr :: 1:  0j_  _1j__
q1 0:@qndivr :: 1:  0j_  _1j_
q1 0:@qndivr :: 1:  0j_   0j__
q1 0:@qndivr :: 1:  0j_   0j_
q1 0:@qndivr :: 1:  0j_   1j__
q1 0:@qndivr :: 1:  0j_   1j_
q1 0:@qndivr :: 1:  0j_   _j__
q1 0:@qndivr :: 1:  0j_   _j_1
q1 0:@qndivr :: 1:  0j_   _
q1 0:@qndivr :: 1:  0j_   _j1
q1 0:@qndivr :: 1:  0j_   _j_
q1 0:@qndivr :: 1:  1j__ __j__
q1 0:@qndivr :: 1:  1j__ __j_1
q1 0:@qndivr :: 1:  1j__ __
q1 0:@qndivr :: 1:  1j__ __j1
q1 0:@qndivr :: 1:  1j__ __j_
q1 0:@qndivr :: 1:  1j__ _1j__
q1 0:@qndivr :: 1:  1j__ _1j_
q1 0:@qndivr :: 1:  1j__  0j__
q1 0:@qndivr :: 1:  1j__  0j_
q1 0:@qndivr :: 1:  1j__  1j__
q1 0:@qndivr :: 1:  1j__  1j_
q1 0:@qndivr :: 1:  1j__  _j__
q1 0:@qndivr :: 1:  1j__  _j_1
q1 0:@qndivr :: 1:  1j__  _
q1 0:@qndivr :: 1:  1j__  _j1
q1 0:@qndivr :: 1:  1j__  _j_
q1 0:@qndivr :: 1:  1j_1 __j__
q1 0:@qndivr :: 1:  1j_1 __j_
q1 0:@qndivr :: 1:  1j_1  _j__
q1 0:@qndivr :: 1:  1j_1  _j_
q1 0:@qndivr :: 1:  1    __j__
q1 0:@qndivr :: 1:  1    __j_
q1 0:@qndivr :: 1:  1     _j__
q1 0:@qndivr :: 1:  1     _j_
q1 0:@qndivr :: 1:  1j1  __j__
q1 0:@qndivr :: 1:  1j1  __j_
q1 0:@qndivr :: 1:  1j1   _j__
q1 0:@qndivr :: 1:  1j1   _j_
q1 0:@qndivr :: 1:  1j_  __j__
q1 0:@qndivr :: 1:  1j_  __j_1
q1 0:@qndivr :: 1:  1j_  __
q1 0:@qndivr :: 1:  1j_  __j1
q1 0:@qndivr :: 1:  1j_  __j_
q1 0:@qndivr :: 1:  1j_  _1j__
q1 0:@qndivr :: 1:  1j_  _1j_
q1 0:@qndivr :: 1:  1j_   0j__
q1 0:@qndivr :: 1:  1j_   0j_
q1 0:@qndivr :: 1:  1j_   1j__
q1 0:@qndivr :: 1:  1j_   1j_
q1 0:@qndivr :: 1:  1j_   _j__
q1 0:@qndivr :: 1:  1j_   _j_1
q1 0:@qndivr :: 1:  1j_   _
q1 0:@qndivr :: 1:  1j_   _j1
q1 0:@qndivr :: 1:  1j_   _j_
q1 0:@qndivr :: 1:  _j__ __j__
q1 0:@qndivr :: 1:  _j__ __j_1
q1 0:@qndivr :: 1:  _j__ __
q1 0:@qndivr :: 1:  _j__ __j1
q1 0:@qndivr :: 1:  _j__ __j_
q1 0:@qndivr :: 1:  _j__ _1j__
q1 0:@qndivr :: 1:  _j__ _1j_1
q1 0:@qndivr :: 1:  _j__ _1
q1 0:@qndivr :: 1:  _j__ _1j1
q1 0:@qndivr :: 1:  _j__ _1j_
q1 0:@qndivr :: 1:  _j__  0j__
q1 0:@qndivr :: 1:  _j__  0j_1
q1 0:@qndivr :: 1:  _j__  0
q1 0:@qndivr :: 1:  _j__  0j1
q1 0:@qndivr :: 1:  _j__  0j_
q1 0:@qndivr :: 1:  _j__  1j__
q1 0:@qndivr :: 1:  _j__  1j_1
q1 0:@qndivr :: 1:  _j__  1
q1 0:@qndivr :: 1:  _j__  1j1
q1 0:@qndivr :: 1:  _j__  1j_
q1 0:@qndivr :: 1:  _j__  _j__
q1 0:@qndivr :: 1:  _j__  _j_1
q1 0:@qndivr :: 1:  _j__  _
q1 0:@qndivr :: 1:  _j__  _j1
q1 0:@qndivr :: 1:  _j__  _j_
q1 0:@qndivr :: 1:  _j_1 __j__
q1 0:@qndivr :: 1:  _j_1 __j_1
q1 0:@qndivr :: 1:  _j_1 __
q1 0:@qndivr :: 1:  _j_1 __j1
q1 0:@qndivr :: 1:  _j_1 __j_
q1 0:@qndivr :: 1:  _j_1 _1j__
q1 0:@qndivr :: 1:  _j_1 _1j_
q1 0:@qndivr :: 1:  _j_1  0j__
q1 0:@qndivr :: 1:  _j_1  0j_
q1 0:@qndivr :: 1:  _j_1  1j__
q1 0:@qndivr :: 1:  _j_1  1j_
q1 0:@qndivr :: 1:  _j_1  _j__
q1 0:@qndivr :: 1:  _j_1  _j_1
q1 0:@qndivr :: 1:  _j_1  _
q1 0:@qndivr :: 1:  _j_1  _j1
q1 0:@qndivr :: 1:  _j_1  _j_
q1 0:@qndivr :: 1:  _    __j__
q1 0:@qndivr :: 1:  _    __j_1
q1 0:@qndivr :: 1:  _    __
q1 0:@qndivr :: 1:  _    __j1
q1 0:@qndivr :: 1:  _    __j_
q1 0:@qndivr :: 1:  _    _1j__
q1 0:@qndivr :: 1:  _    _1j_
q1 0:@qndivr :: 1:  _     0j__
q1 0:@qndivr :: 1:  _     0j_
q1 0:@qndivr :: 1:  _     1j__
q1 0:@qndivr :: 1:  _     1j_
q1 0:@qndivr :: 1:  _     _j__
q1 0:@qndivr :: 1:  _     _j_1
q1 0:@qndivr :: 1:  _     _
q1 0:@qndivr :: 1:  _     _j1
q1 0:@qndivr :: 1:  _     _j_
q1 0:@qndivr :: 1:  _j1  __j__
q1 0:@qndivr :: 1:  _j1  __j_1
q1 0:@qndivr :: 1:  _j1  __
q1 0:@qndivr :: 1:  _j1  __j1
q1 0:@qndivr :: 1:  _j1  __j_
q1 0:@qndivr :: 1:  _j1  _1j__
q1 0:@qndivr :: 1:  _j1  _1j_
q1 0:@qndivr :: 1:  _j1   0j__
q1 0:@qndivr :: 1:  _j1   0j_
q1 0:@qndivr :: 1:  _j1   1j__
q1 0:@qndivr :: 1:  _j1   1j_
q1 0:@qndivr :: 1:  _j1   _j__
q1 0:@qndivr :: 1:  _j1   _j_1
q1 0:@qndivr :: 1:  _j1   _
q1 0:@qndivr :: 1:  _j1   _j1
q1 0:@qndivr :: 1:  _j1   _j_
q1 0:@qndivr :: 1:  _j_  __j__
q1 0:@qndivr :: 1:  _j_  __j_1
q1 0:@qndivr :: 1:  _j_  __
q1 0:@qndivr :: 1:  _j_  __j1
q1 0:@qndivr :: 1:  _j_  __j_
q1 0:@qndivr :: 1:  _j_  _1j__
q1 0:@qndivr :: 1:  _j_  _1j_1
q1 0:@qndivr :: 1:  _j_  _1
q1 0:@qndivr :: 1:  _j_  _1j1
q1 0:@qndivr :: 1:  _j_  _1j_
q1 0:@qndivr :: 1:  _j_   0j__
q1 0:@qndivr :: 1:  _j_   0j_1
q1 0:@qndivr :: 1:  _j_   0
q1 0:@qndivr :: 1:  _j_   0j1
q1 0:@qndivr :: 1:  _j_   0j_
q1 0:@qndivr :: 1:  _j_   1j__
q1 0:@qndivr :: 1:  _j_   1j_1
q1 0:@qndivr :: 1:  _j_   1
q1 0:@qndivr :: 1:  _j_   1j1
q1 0:@qndivr :: 1:  _j_   1j_
q1 0:@qndivr :: 1:  _j_   _j__
q1 0:@qndivr :: 1:  _j_   _j_1
q1 0:@qndivr :: 1:  _j_   _
q1 0:@qndivr :: 1:  _j_   _j1
q1 0:@qndivr :: 1:  _j_   _j_
NB. - input contains directed infinity
__j_   _j_  -: (__ qn1 q1) qndivr        q2
 _j__ __j__ -: ( _ qn1 q1) qndivr        q2
__j__  _j__ -: (__ qni q1) qndivr        q2
 _j_  __j_  -: ( _ qni q1) qndivr        q2
__j__ __j_  -: (__ qnj q1) qndivr        q2
 _j_   _j__ -: ( _ qnj q1) qndivr        q2
__j_  __j__ -: (__ qnk q1) qndivr        q2
 _j__  _j_  -: ( _ qnk q1) qndivr        q2
 0     0    -:         q1  qndivr __ qn1 q2
 0     0    -:         q1  qndivr __ qni q2
 0     0    -:         q1  qndivr __ qnj q2
 0     0    -:         q1  qndivr __ qnk q2
 0     0    -:         q1  qndivr  _ qn1 q2
 0     0    -:         q1  qndivr  _ qni q2
 0     0    -:         q1  qndivr  _ qnj q2
 0     0    -:         q1  qndivr  _ qnk q2
 0     0    -: (__ qn1 q1) qndivr __ qn1 q2
 0     0    -: (__ qn1 q1) qndivr __ qni q2
 0     0    -: (__ qn1 q1) qndivr __ qnj q2
 0     0    -: (__ qn1 q1) qndivr __ qnk q2
 0     0    -: (__ qn1 q1) qndivr  _ qn1 q2
 0     0    -: (__ qn1 q1) qndivr  _ qni q2
 0     0    -: (__ qn1 q1) qndivr  _ qnj q2
 0     0    -: (__ qn1 q1) qndivr  _ qnk q2
 0     0    -: ( _ qn1 q1) qndivr __ qn1 q2
 0     0    -: ( _ qn1 q1) qndivr __ qni q2
 0     0    -: ( _ qn1 q1) qndivr __ qnj q2
 0     0    -: ( _ qn1 q1) qndivr __ qnk q2
 0     0    -: ( _ qn1 q1) qndivr  _ qn1 q2
 0     0    -: ( _ qn1 q1) qndivr  _ qni q2
 0     0    -: ( _ qn1 q1) qndivr  _ qnj q2
 0     0    -: ( _ qn1 q1) qndivr  _ qnk q2
(5 2j1 ,: 0) -: q1 qndivr q3 ,: _ qn1 q4
(5 2j1 ,: 0) -: q1 qndivr q3 ,: _ qni q4
(5 2j1 ,: 0) -: q1 qndivr q3 ,: _ qnj q4
(5 2j1 ,: 0) -: q1 qndivr q3 ,: _ qnk q4
(5 2j1 ,: _j__ __j__) -: (q1 ,: _ qn1 q2) qndivr q3
(5 2j1 ,: _j_  __j_ ) -: (q1 ,: _ qni q2) qndivr q3
(5 2j1 ,: _j_   _j__) -: (q1 ,: _ qnj q2) qndivr q3
(5 2j1 ,: _j__  _j_ ) -: (q1 ,: _ qnk q2) qndivr q3
NB. - edge cases input
1     1    -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivr (FP_OVFL  j.  FP_OVFL) , 0
1    _1    -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivr  0                     , FP_OVFL j. FP_OVFL
0.5  _0.5  -: ((FP_OVFL j. FP_OVFL) , 0                 ) qndivr (FP_OVFL  j.  FP_OVFL) , FP_OVFL j. FP_OVFL
0.5   0.5  -: ( 0                   , FP_OVFL j. FP_OVFL) qndivr (FP_OVFL  j.  FP_OVFL) , FP_OVFL j. FP_OVFL
_j_   _j_  -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivr  FP_UNFL               , 0
_j__  _j__ -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivr (0        j. FP_UNFL)  , 0
_j__ __j_  -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivr  0                     , FP_UNFL j.   0
_j_  __j__ -: ((FP_OVFL j. FP_OVFL) , FP_OVFL j. FP_OVFL) qndivr  0                     , 0       j. FP_UNFL
0     0    -: ((FP_UNFL j. FP_UNFL) , FP_UNFL j. FP_UNFL) qndivr (FP_OVFL  j.  FP_OVFL) , FP_OVFL j. FP_OVFL
0     0    -: ((FP_UNFL j. FP_UNFL) , FP_UNFL j. FP_UNFL) qndivr (FP_OVFL  j.  FP_OVFL) , 0
0     0    -: ((FP_UNFL j. FP_UNFL) , FP_UNFL j. FP_UNFL) qndivr  0                     , FP_OVFL j. FP_OVFL
0     0    -: ((FP_UNFL j. FP_UNFL) , 0                 ) qndivr (FP_OVFL  j.  FP_OVFL) , FP_OVFL j. FP_OVFL
0     0    -: ( 0                   , FP_UNFL j. FP_UNFL) qndivr (FP_OVFL  j.  FP_OVFL) , FP_OVFL j. FP_OVFL
1     1    -: ((FP_UNFL j. FP_UNFL) , FP_UNFL j. FP_UNFL) qndivr (FP_UNFL  j.  FP_UNFL) , 0
1    _1    -: ((FP_UNFL j. FP_UNFL) , FP_UNFL j. FP_UNFL) qndivr  0                     , FP_UNFL j. FP_UNFL
0.5  _0.5  -: ((FP_UNFL j. FP_UNFL) , 0                 ) qndivr (FP_UNFL  j.  FP_UNFL) , FP_UNFL j. FP_UNFL
0.5   0.5  -: ( 0                   , FP_UNFL j. FP_UNFL) qndivr (FP_UNFL  j.  FP_UNFL) , FP_UNFL j. FP_UNFL
_     0            -: q5 qndivl q7
(0 ,~ FP_UNFL % 8) -: q7 qndivl q5
NB. - input without edge cases
 0 0 -:   0 0 qndivr  q1
 1 0 -:   q1  qndivr  q1
 1 0 -:   q2  qndivr  q2
_1 0 -:   q1  qndivr -q1
_1 0 -: (-q1) qndivr  q1
_1 0 -:   q2  qndivr -q2
_1 0 -: (-q2) qndivr  q2
q1 -: q1 qndivr 1 0
q2 -: _60j12 30j24 qndivr q1
q1 -: _60j20 14j32 qndivr q2
(1 0 ,: q1) -:  q1        qndivr q1  ,: 1 0
(q1  ,: q2) -: (q1 ,: q2) qndivr 1 0
(1 0 ,: q2) -: (q1 ,: q2) qndivr q1  ,: 1 0
NB. (-: qndivr~) i. 0 2  NB. this fails until https://github.com/jsoftware/jsource/issues/192 has been fixed

NB. qnf
 1x1                                       0                                        -: ^ qnf qnmark1  q1
_1.13120438375681354j2.47172667200481877   0                                        -: ^ qnf qnmark1i q1
 1.69392272368329944j_0.78955962454155892 _1.18433943681233833j_1.57911924908311785 -: ^ qnf          q1

0                                       0                                       -: ^. qnf qnmark1  q1
0.80471895621705025j1.10714871779409041 0                                       -: ^. qnf qnmark1i q1
1.70059869083107773j0.51519029266408511 0.77278543899612773j1.03038058532817023 -: ^. qnf          q1
