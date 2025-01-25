NB. Verify rot verbs
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

'q3 q4'=. ,.~ (j.~"0) 1r2 1r8
sqr_big=. 2 *&%: FP_OVFL               NB. its square overflows
cs=. qnsign , 0 1 1 ]`(j./)/. gemat 3  NB. random rotator

NB. =========================================================
NB. Verification suite

NB. lartg
NB. - input contains NaN
1 1 -: isnan lartg  0     0j_.
1 1 -: isnan lartg  0    _.     NB. LAPACK's DLARTG gives (c=0,s=1,r=NaN), ZLARTG gives (c=0,s=NaN+i*NaN,r=NaN+i*0)
1 1 -: isnan lartg  0j_.  0
1 1 -: isnan lartg _.     0     NB. LAPACK's xLARTG gives (c=1,s=0,r=NaN)
1 1 -: isnan lartg  0    _.j_.
1 1 -: isnan lartg  0j_.  0j_.
1 1 -: isnan lartg  0j_. _.
1 1 -: isnan lartg _.     0j_.
1 1 -: isnan lartg _.    _.
1 1 -: isnan lartg _.j_.  0
1 1 -: isnan lartg _.j_. _.
1 1 -: isnan lartg _.j_.  0j_.
1 1 -: isnan lartg _.    _.j_.
1 1 -: isnan lartg  0j_. _.j_.
1 1 -: isnan lartg _.j_. _.j_.
(2 2 $ 1) -: isnan lartg _. 0 ,: 0 _.
(,.~ 1 0) -: isnan lartg _. 0 ,: q3
(,.~ 0 1) -: isnan lartg q3   ,: 0 _.
NB. - input contains q∞ so NaN error must be throwed
0:@lartg :: 1: __j__ __j__
0:@lartg :: 1: __j__ __j_1
0:@lartg :: 1: __j__ __
0:@lartg :: 1: __j__ __j1
0:@lartg :: 1: __j__ __j_
0:@lartg :: 1: __j__ _1j__
0:@lartg :: 1: __j__ _1j_1
0:@lartg :: 1: __j__ _1
0:@lartg :: 1: __j__ _1j1
0:@lartg :: 1: __j__ _1j_
0:@lartg :: 1: __j__  0j__
0:@lartg :: 1: __j__  0j_1
0:@lartg :: 1: __j__  0
0:@lartg :: 1: __j__  0j1
0:@lartg :: 1: __j__  0j_
0:@lartg :: 1: __j__  1j__
0:@lartg :: 1: __j__  1j_1
0:@lartg :: 1: __j__  1
0:@lartg :: 1: __j__  1j1
0:@lartg :: 1: __j__  1j_
0:@lartg :: 1: __j__  _j__
0:@lartg :: 1: __j__  _j_1
0:@lartg :: 1: __j__  _
0:@lartg :: 1: __j__  _j1
0:@lartg :: 1: __j__  _j_
0:@lartg :: 1: __j_1 __j__
0:@lartg :: 1: __j_1 __j_1
0:@lartg :: 1: __j_1 __
0:@lartg :: 1: __j_1 __j1
0:@lartg :: 1: __j_1 __j_
0:@lartg :: 1: __j_1 _1j__
0:@lartg :: 1: __j_1 _1j_
0:@lartg :: 1: __j_1  0j__
0:@lartg :: 1: __j_1  0j_
0:@lartg :: 1: __j_1  1j__
0:@lartg :: 1: __j_1  1j_
0:@lartg :: 1: __j_1  _j__
0:@lartg :: 1: __j_1  _j_1
0:@lartg :: 1: __j_1  _
0:@lartg :: 1: __j_1  _j1
0:@lartg :: 1: __j_1  _j_
0:@lartg :: 1: __    __j__
0:@lartg :: 1: __    __j_1
0:@lartg :: 1: __    __
0:@lartg :: 1: __    __j1
0:@lartg :: 1: __    __j_
0:@lartg :: 1: __    _1j__
0:@lartg :: 1: __    _1j_
0:@lartg :: 1: __     0j__
0:@lartg :: 1: __     0j_
0:@lartg :: 1: __     1j__
0:@lartg :: 1: __     1j_
0:@lartg :: 1: __     _j__
0:@lartg :: 1: __     _j_1
0:@lartg :: 1: __     _
0:@lartg :: 1: __     _j1
0:@lartg :: 1: __     _j_
0:@lartg :: 1: __j1  __j__
0:@lartg :: 1: __j1  __j_1
0:@lartg :: 1: __j1  __
0:@lartg :: 1: __j1  __j1
0:@lartg :: 1: __j1  __j_
0:@lartg :: 1: __j1  _1j__
0:@lartg :: 1: __j1  _1j_
0:@lartg :: 1: __j1   0j__
0:@lartg :: 1: __j1   0j_
0:@lartg :: 1: __j1   1j__
0:@lartg :: 1: __j1   1j_
0:@lartg :: 1: __j1   _j__
0:@lartg :: 1: __j1   _j_1
0:@lartg :: 1: __j1   _
0:@lartg :: 1: __j1   _j1
0:@lartg :: 1: __j1   _j_
0:@lartg :: 1: __j_  __j__
0:@lartg :: 1: __j_  __j_1
0:@lartg :: 1: __j_  __
0:@lartg :: 1: __j_  __j1
0:@lartg :: 1: __j_  __j_
0:@lartg :: 1: __j_  _1j__
0:@lartg :: 1: __j_  _1j_1
0:@lartg :: 1: __j_  _1
0:@lartg :: 1: __j_  _1j1
0:@lartg :: 1: __j_  _1j_
0:@lartg :: 1: __j_   0j__
0:@lartg :: 1: __j_   0j_1
0:@lartg :: 1: __j_   0
0:@lartg :: 1: __j_   0j1
0:@lartg :: 1: __j_   0j_
0:@lartg :: 1: __j_   1j__
0:@lartg :: 1: __j_   1j_1
0:@lartg :: 1: __j_   1
0:@lartg :: 1: __j_   1j1
0:@lartg :: 1: __j_   1j_
0:@lartg :: 1: __j_   _j__
0:@lartg :: 1: __j_   _j_1
0:@lartg :: 1: __j_   _
0:@lartg :: 1: __j_   _j1
0:@lartg :: 1: __j_   _j_
0:@lartg :: 1: _1j__ __j__
0:@lartg :: 1: _1j__ __j_1
0:@lartg :: 1: _1j__ __
0:@lartg :: 1: _1j__ __j1
0:@lartg :: 1: _1j__ __j_
0:@lartg :: 1: _1j__ _1j__
0:@lartg :: 1: _1j__ _1j_
0:@lartg :: 1: _1j__  0j__
0:@lartg :: 1: _1j__  0j_
0:@lartg :: 1: _1j__  1j__
0:@lartg :: 1: _1j__  1j_
0:@lartg :: 1: _1j__  _j__
0:@lartg :: 1: _1j__  _j_1
0:@lartg :: 1: _1j__  _
0:@lartg :: 1: _1j__  _j1
0:@lartg :: 1: _1j__  _j_
0:@lartg :: 1: _1j_1 __j__
0:@lartg :: 1: _1j_1 __j_
0:@lartg :: 1: _1j_1  _j__
0:@lartg :: 1: _1j_1  _j_
0:@lartg :: 1: _1    __j__
0:@lartg :: 1: _1    __j_
0:@lartg :: 1: _1     _j__
0:@lartg :: 1: _1     _j_
0:@lartg :: 1: _1j1  __j__
0:@lartg :: 1: _1j1  __j_
0:@lartg :: 1: _1j1   _j__
0:@lartg :: 1: _1j1   _j_
0:@lartg :: 1: _1j_  __j__
0:@lartg :: 1: _1j_  __j_1
0:@lartg :: 1: _1j_  __
0:@lartg :: 1: _1j_  __j1
0:@lartg :: 1: _1j_  __j_
0:@lartg :: 1: _1j_  _1j__
0:@lartg :: 1: _1j_  _1j_
0:@lartg :: 1: _1j_   0j__
0:@lartg :: 1: _1j_   0j_
0:@lartg :: 1: _1j_   1j__
0:@lartg :: 1: _1j_   1j_
0:@lartg :: 1: _1j_   _j__
0:@lartg :: 1: _1j_   _j_1
0:@lartg :: 1: _1j_   _
0:@lartg :: 1: _1j_   _j1
0:@lartg :: 1: _1j_   _j_
0:@lartg :: 1:  0j__ __j__
0:@lartg :: 1:  0j__ __j_1
0:@lartg :: 1:  0j__ __
0:@lartg :: 1:  0j__ __j1
0:@lartg :: 1:  0j__ __j_
0:@lartg :: 1:  0j__ _1j__
0:@lartg :: 1:  0j__ _1j_
0:@lartg :: 1:  0j__  0j__
0:@lartg :: 1:  0j__  0j_
0:@lartg :: 1:  0j__  1j__
0:@lartg :: 1:  0j__  1j_
0:@lartg :: 1:  0j__  _j__
0:@lartg :: 1:  0j__  _j_1
0:@lartg :: 1:  0j__  _
0:@lartg :: 1:  0j__  _j1
0:@lartg :: 1:  0j__  _j_
0:@lartg :: 1:  0j_1 __j__
0:@lartg :: 1:  0j_1 __j_
0:@lartg :: 1:  0j_1  _j__
0:@lartg :: 1:  0j_1  _j_
0:@lartg :: 1:  0    __j__
0:@lartg :: 1:  0    __j_
0:@lartg :: 1:  0     _j__
0:@lartg :: 1:  0     _j_
0:@lartg :: 1:  0j1  __j__
0:@lartg :: 1:  0j1  __j_
0:@lartg :: 1:  0j1   _j__
0:@lartg :: 1:  0j1   _j_
0:@lartg :: 1:  0j_  __j__
0:@lartg :: 1:  0j_  __j_1
0:@lartg :: 1:  0j_  __
0:@lartg :: 1:  0j_  __j1
0:@lartg :: 1:  0j_  __j_
0:@lartg :: 1:  0j_  _1j__
0:@lartg :: 1:  0j_  _1j_
0:@lartg :: 1:  0j_   0j__
0:@lartg :: 1:  0j_   0j_
0:@lartg :: 1:  0j_   1j__
0:@lartg :: 1:  0j_   1j_
0:@lartg :: 1:  0j_   _j__
0:@lartg :: 1:  0j_   _j_1
0:@lartg :: 1:  0j_   _
0:@lartg :: 1:  0j_   _j1
0:@lartg :: 1:  0j_   _j_
0:@lartg :: 1:  1j__ __j__
0:@lartg :: 1:  1j__ __j_1
0:@lartg :: 1:  1j__ __
0:@lartg :: 1:  1j__ __j1
0:@lartg :: 1:  1j__ __j_
0:@lartg :: 1:  1j__ _1j__
0:@lartg :: 1:  1j__ _1j_
0:@lartg :: 1:  1j__  0j__
0:@lartg :: 1:  1j__  0j_
0:@lartg :: 1:  1j__  1j__
0:@lartg :: 1:  1j__  1j_
0:@lartg :: 1:  1j__  _j__
0:@lartg :: 1:  1j__  _j_1
0:@lartg :: 1:  1j__  _
0:@lartg :: 1:  1j__  _j1
0:@lartg :: 1:  1j__  _j_
0:@lartg :: 1:  1j_1 __j__
0:@lartg :: 1:  1j_1 __j_
0:@lartg :: 1:  1j_1  _j__
0:@lartg :: 1:  1j_1  _j_
0:@lartg :: 1:  1    __j__
0:@lartg :: 1:  1    __j_
0:@lartg :: 1:  1     _j__
0:@lartg :: 1:  1     _j_
0:@lartg :: 1:  1j1  __j__
0:@lartg :: 1:  1j1  __j_
0:@lartg :: 1:  1j1   _j__
0:@lartg :: 1:  1j1   _j_
0:@lartg :: 1:  1j_  __j__
0:@lartg :: 1:  1j_  __j_1
0:@lartg :: 1:  1j_  __
0:@lartg :: 1:  1j_  __j1
0:@lartg :: 1:  1j_  __j_
0:@lartg :: 1:  1j_  _1j__
0:@lartg :: 1:  1j_  _1j_
0:@lartg :: 1:  1j_   0j__
0:@lartg :: 1:  1j_   0j_
0:@lartg :: 1:  1j_   1j__
0:@lartg :: 1:  1j_   1j_
0:@lartg :: 1:  1j_   _j__
0:@lartg :: 1:  1j_   _j_1
0:@lartg :: 1:  1j_   _
0:@lartg :: 1:  1j_   _j1
0:@lartg :: 1:  1j_   _j_
0:@lartg :: 1:  _j__ __j__
0:@lartg :: 1:  _j__ __j_1
0:@lartg :: 1:  _j__ __
0:@lartg :: 1:  _j__ __j1
0:@lartg :: 1:  _j__ __j_
0:@lartg :: 1:  _j__ _1j__
0:@lartg :: 1:  _j__ _1j_1
0:@lartg :: 1:  _j__ _1
0:@lartg :: 1:  _j__ _1j1
0:@lartg :: 1:  _j__ _1j_
0:@lartg :: 1:  _j__  0j__
0:@lartg :: 1:  _j__  0j_1
0:@lartg :: 1:  _j__  0
0:@lartg :: 1:  _j__  0j1
0:@lartg :: 1:  _j__  0j_
0:@lartg :: 1:  _j__  1j__
0:@lartg :: 1:  _j__  1j_1
0:@lartg :: 1:  _j__  1
0:@lartg :: 1:  _j__  1j1
0:@lartg :: 1:  _j__  1j_
0:@lartg :: 1:  _j__  _j__
0:@lartg :: 1:  _j__  _j_1
0:@lartg :: 1:  _j__  _
0:@lartg :: 1:  _j__  _j1
0:@lartg :: 1:  _j__  _j_
0:@lartg :: 1:  _j_1 __j__
0:@lartg :: 1:  _j_1 __j_1
0:@lartg :: 1:  _j_1 __
0:@lartg :: 1:  _j_1 __j1
0:@lartg :: 1:  _j_1 __j_
0:@lartg :: 1:  _j_1 _1j__
0:@lartg :: 1:  _j_1 _1j_
0:@lartg :: 1:  _j_1  0j__
0:@lartg :: 1:  _j_1  0j_
0:@lartg :: 1:  _j_1  1j__
0:@lartg :: 1:  _j_1  1j_
0:@lartg :: 1:  _j_1  _j__
0:@lartg :: 1:  _j_1  _j_1
0:@lartg :: 1:  _j_1  _
0:@lartg :: 1:  _j_1  _j1
0:@lartg :: 1:  _j_1  _j_
0:@lartg :: 1:  _    __j__
0:@lartg :: 1:  _    __j_1
0:@lartg :: 1:  _    __
0:@lartg :: 1:  _    __j1
0:@lartg :: 1:  _    __j_
0:@lartg :: 1:  _    _1j__
0:@lartg :: 1:  _    _1j_
0:@lartg :: 1:  _     0j__
0:@lartg :: 1:  _     0j_
0:@lartg :: 1:  _     1j__
0:@lartg :: 1:  _     1j_
0:@lartg :: 1:  _     _j__
0:@lartg :: 1:  _     _j_1
0:@lartg :: 1:  _     _
0:@lartg :: 1:  _     _j1
0:@lartg :: 1:  _     _j_
0:@lartg :: 1:  _j1  __j__
0:@lartg :: 1:  _j1  __j_1
0:@lartg :: 1:  _j1  __
0:@lartg :: 1:  _j1  __j1
0:@lartg :: 1:  _j1  __j_
0:@lartg :: 1:  _j1  _1j__
0:@lartg :: 1:  _j1  _1j_
0:@lartg :: 1:  _j1   0j__
0:@lartg :: 1:  _j1   0j_
0:@lartg :: 1:  _j1   1j__
0:@lartg :: 1:  _j1   1j_
0:@lartg :: 1:  _j1   _j__
0:@lartg :: 1:  _j1   _j_1
0:@lartg :: 1:  _j1   _
0:@lartg :: 1:  _j1   _j1
0:@lartg :: 1:  _j1   _j_
0:@lartg :: 1:  _j_  __j__
0:@lartg :: 1:  _j_  __j_1
0:@lartg :: 1:  _j_  __
0:@lartg :: 1:  _j_  __j1
0:@lartg :: 1:  _j_  __j_
0:@lartg :: 1:  _j_  _1j__
0:@lartg :: 1:  _j_  _1j_1
0:@lartg :: 1:  _j_  _1
0:@lartg :: 1:  _j_  _1j1
0:@lartg :: 1:  _j_  _1j_
0:@lartg :: 1:  _j_   0j__
0:@lartg :: 1:  _j_   0j_1
0:@lartg :: 1:  _j_   0
0:@lartg :: 1:  _j_   0j1
0:@lartg :: 1:  _j_   0j_
0:@lartg :: 1:  _j_   1j__
0:@lartg :: 1:  _j_   1j_1
0:@lartg :: 1:  _j_   1
0:@lartg :: 1:  _j_   1j1
0:@lartg :: 1:  _j_   1j_
0:@lartg :: 1:  _j_   _j__
0:@lartg :: 1:  _j_   _j_1
0:@lartg :: 1:  _j_   _
0:@lartg :: 1:  _j_   _j1
0:@lartg :: 1:  _j_   _j_
NB. - input contains directed infinity and is not trivial
NB.   so NaN error must be throwed
0:@lartg :: 1  ( __        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1  ( __        j. -FP_OVFL) ,  -FP_OVFL
0:@lartg :: 1  ( __        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1  ( __        j. -FP_OVFL) ,   0        j. -FP_OVFL
0:@lartg :: 1  ( __        j. -FP_OVFL) ,   0        j.  FP_OVFL
0:@lartg :: 1  ( __        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1  ( __        j. -FP_OVFL) ,   FP_OVFL
0:@lartg :: 1  ( __        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
0:@lartg :: 1:   __                     , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1:   __                     ,  -FP_OVFL
0:@lartg :: 1:   __                     , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1:   __                     ,   0        j. -FP_OVFL
0:@lartg :: 1:   __                     ,   0        j.  FP_OVFL
0:@lartg :: 1:   __                     ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1:   __                     ,   FP_OVFL
0:@lartg :: 1:   __                     ,   FP_OVFL  j.  FP_OVFL
0:@lartg :: 1  ( __        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1  ( __        j.  FP_OVFL) ,  -FP_OVFL
0:@lartg :: 1  ( __        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1  ( __        j.  FP_OVFL) ,   0        j. -FP_OVFL
0:@lartg :: 1  ( __        j.  FP_OVFL) ,   0        j.  FP_OVFL
0:@lartg :: 1  ( __        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1  ( __        j.  FP_OVFL) ,   FP_OVFL
0:@lartg :: 1  ( __        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. __      ) , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,  -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. __      ) , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,   0        j. -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,   0        j.  FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,   FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,   FP_OVFL  j.  FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,  __        j. -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,  __
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,  __        j.  FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. __
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  _
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   0        j. __
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  _
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. __
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  _
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   _        j. -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   _
0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   _        j.  FP_OVFL
0:@lartg :: 1   (-FP_OVFL)              ,  __        j. -FP_OVFL
0:@lartg :: 1   (-FP_OVFL)              ,  __
0:@lartg :: 1   (-FP_OVFL)              ,  __        j.  FP_OVFL
0:@lartg :: 1   (-FP_OVFL)              , (-FP_OVFL) j. __
0:@lartg :: 1   (-FP_OVFL)              , (-FP_OVFL) j.  _
0:@lartg :: 1   (-FP_OVFL)              ,   0        j. __
0:@lartg :: 1   (-FP_OVFL)              ,   0        j.  _
0:@lartg :: 1   (-FP_OVFL)              ,   FP_OVFL  j. __
0:@lartg :: 1   (-FP_OVFL)              ,   FP_OVFL  j.  _
0:@lartg :: 1   (-FP_OVFL)              ,   _        j. -FP_OVFL
0:@lartg :: 1   (-FP_OVFL)              ,   _
0:@lartg :: 1   (-FP_OVFL)              ,   _        j.  FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,  __        j. -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,  __
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,  __        j.  FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. __
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  _
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   0        j. __
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  _
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. __
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  _
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   _        j. -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   _
0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   _        j.  FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  _      ) , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,  -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  _      ) , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,   0        j. -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,   0        j.  FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,   FP_OVFL
0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,   FP_OVFL  j.  FP_OVFL
0:@lartg :: 1  (  0        j. __      ) , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1  (  0        j. __      ) ,  -FP_OVFL
0:@lartg :: 1  (  0        j. __      ) , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1  (  0        j. __      ) ,   0        j. -FP_OVFL
0:@lartg :: 1  (  0        j. __      ) ,   0        j.  FP_OVFL
0:@lartg :: 1  (  0        j. __      ) ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1  (  0        j. __      ) ,   FP_OVFL
0:@lartg :: 1  (  0        j. __      ) ,   FP_OVFL  j.  FP_OVFL
0:@lartg :: 1  (  0        j. -FP_OVFL) ,  __        j. -FP_OVFL
0:@lartg :: 1  (  0        j. -FP_OVFL) ,  __
0:@lartg :: 1  (  0        j. -FP_OVFL) ,  __        j.  FP_OVFL
0:@lartg :: 1  (  0        j. -FP_OVFL) , (-FP_OVFL) j. __
0:@lartg :: 1  (  0        j. -FP_OVFL) , (-FP_OVFL) j.  _
0:@lartg :: 1  (  0        j. -FP_OVFL) ,   0        j. __
0:@lartg :: 1  (  0        j. -FP_OVFL) ,   0        j.  _
0:@lartg :: 1  (  0        j. -FP_OVFL) ,   FP_OVFL  j. __
0:@lartg :: 1  (  0        j. -FP_OVFL) ,   FP_OVFL  j.  _
0:@lartg :: 1  (  0        j. -FP_OVFL) ,   _        j. -FP_OVFL
0:@lartg :: 1  (  0        j. -FP_OVFL) ,   _
0:@lartg :: 1  (  0        j. -FP_OVFL) ,   _        j.  FP_OVFL
0:@lartg :: 1  (  0        j.  FP_OVFL) ,  __        j. -FP_OVFL
0:@lartg :: 1  (  0        j.  FP_OVFL) ,  __
0:@lartg :: 1  (  0        j.  FP_OVFL) ,  __        j.  FP_OVFL
0:@lartg :: 1  (  0        j.  FP_OVFL) , (-FP_OVFL) j. __
0:@lartg :: 1  (  0        j.  FP_OVFL) , (-FP_OVFL) j.  _
0:@lartg :: 1  (  0        j.  FP_OVFL) ,   0        j. __
0:@lartg :: 1  (  0        j.  FP_OVFL) ,   0        j.  _
0:@lartg :: 1  (  0        j.  FP_OVFL) ,   FP_OVFL  j. __
0:@lartg :: 1  (  0        j.  FP_OVFL) ,   FP_OVFL  j.  _
0:@lartg :: 1  (  0        j.  FP_OVFL) ,   _        j. -FP_OVFL
0:@lartg :: 1  (  0        j.  FP_OVFL) ,   _
0:@lartg :: 1  (  0        j.  FP_OVFL) ,   _        j.  FP_OVFL
0:@lartg :: 1  (  0        j.  _      ) , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1  (  0        j.  _      ) ,  -FP_OVFL
0:@lartg :: 1  (  0        j.  _      ) , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1  (  0        j.  _      ) ,   0        j. -FP_OVFL
0:@lartg :: 1  (  0        j.  _      ) ,   0        j.  FP_OVFL
0:@lartg :: 1  (  0        j.  _      ) ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1  (  0        j.  _      ) ,   FP_OVFL
0:@lartg :: 1  (  0        j.  _      ) ,   FP_OVFL  j.  FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. __      ) , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. __      ) ,  -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. __      ) , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. __      ) ,   0        j. -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. __      ) ,   0        j.  FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. __      ) ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. __      ) ,   FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. __      ) ,   FP_OVFL  j.  FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,  __        j. -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,  __
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,  __        j.  FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. __
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  _
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   0        j. __
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   0        j.  _
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. __
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  _
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   _        j. -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   _
0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   _        j.  FP_OVFL
0:@lartg :: 1:    FP_OVFL               ,  __        j. -FP_OVFL
0:@lartg :: 1:    FP_OVFL               ,  __
0:@lartg :: 1:    FP_OVFL               ,  __        j.  FP_OVFL
0:@lartg :: 1:    FP_OVFL               , (-FP_OVFL) j. __
0:@lartg :: 1:    FP_OVFL               , (-FP_OVFL) j.  _
0:@lartg :: 1:    FP_OVFL               ,   0        j. __
0:@lartg :: 1:    FP_OVFL               ,   0        j.  _
0:@lartg :: 1:    FP_OVFL               ,   FP_OVFL  j. __
0:@lartg :: 1:    FP_OVFL               ,   FP_OVFL  j.  _
0:@lartg :: 1:    FP_OVFL               ,   _        j. -FP_OVFL
0:@lartg :: 1:    FP_OVFL               ,   _
0:@lartg :: 1:    FP_OVFL               ,   _        j.  FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,  __        j. -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,  __
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,  __        j.  FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. __
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  _
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   0        j. __
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   0        j.  _
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. __
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  _
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   _        j. -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   _
0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   _        j.  FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  _      ) , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,  -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  _      ) , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,   0        j. -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,   0        j.  FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,   FP_OVFL
0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,   FP_OVFL  j.  FP_OVFL
0:@lartg :: 1  (  _        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1  (  _        j. -FP_OVFL) ,  -FP_OVFL
0:@lartg :: 1  (  _        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1  (  _        j. -FP_OVFL) ,   0        j. -FP_OVFL
0:@lartg :: 1  (  _        j. -FP_OVFL) ,   0        j.  FP_OVFL
0:@lartg :: 1  (  _        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1  (  _        j. -FP_OVFL) ,   FP_OVFL
0:@lartg :: 1  (  _        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
0:@lartg :: 1:    _                     , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1:    _                     ,  -FP_OVFL
0:@lartg :: 1:    _                     , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1:    _                     ,   0        j. -FP_OVFL
0:@lartg :: 1:    _                     ,   0        j.  FP_OVFL
0:@lartg :: 1:    _                     ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1:    _                     ,   FP_OVFL
0:@lartg :: 1:    _                     ,   FP_OVFL  j.  FP_OVFL
0:@lartg :: 1  (  _        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
0:@lartg :: 1  (  _        j.  FP_OVFL) ,  -FP_OVFL
0:@lartg :: 1  (  _        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
0:@lartg :: 1  (  _        j.  FP_OVFL) ,   0        j. -FP_OVFL
0:@lartg :: 1  (  _        j.  FP_OVFL) ,   0        j.  FP_OVFL
0:@lartg :: 1  (  _        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
0:@lartg :: 1  (  _        j.  FP_OVFL) ,   FP_OVFL
0:@lartg :: 1  (  _        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
NB. - input contains directed infinity and is trivial
1  0    -: lartg ( __        j. -FP_OVFL) ,    0
1  0    -: lartg ( __        j. _1      ) ,    0
1  0    -: lartg ( __        j. -FP_UNFL) ,    0
1  0    -: lartg ( __        j.  FP_UNFL) ,    0
1  0    -: lartg ( __        j.  1      ) ,    0
1  0    -: lartg ( __        j.  FP_OVFL) ,    0
1  0    -: lartg ((-FP_OVFL) j. __      ) ,    0
1  0    -: lartg ((-FP_OVFL) j.  _      ) ,    0
1  0    -: lartg ( _1        j. __      ) ,    0
1  0    -: lartg ( _1        j.  _      ) ,    0
1  0    -: lartg ((-FP_UNFL) j. __      ) ,    0
1  0    -: lartg ((-FP_UNFL) j.  _      ) ,    0
0 _1    -: lartg    0                     , ( __        j. -FP_OVFL)
0 _1    -: lartg    0                     , ( __        j. _1      )
0 _1    -: lartg    0                     , ( __        j. -FP_UNFL)
0 _1    -: lartg    0                     , ( __        j.  FP_UNFL)
0 _1    -: lartg    0                     , ( __        j.  1      )
0 _1    -: lartg    0                     , ( __        j.  FP_OVFL)
0  0j1  -: lartg    0                     , ((-FP_OVFL) j. __      )
0  0j_1 -: lartg    0                     , ((-FP_OVFL) j.  _      )
0  0j1  -: lartg    0                     , ( _1        j. __      )
0  0j_1 -: lartg    0                     , ( _1        j.  _      )
0  0j1  -: lartg    0                     , ((-FP_UNFL) j. __      )
0  0j_1 -: lartg    0                     , ((-FP_UNFL) j.  _      )
0  0j1  -: lartg    0                     , (  FP_UNFL  j. __      )
0  0j_1 -: lartg    0                     , (  FP_UNFL  j.  _      )
0  0j1  -: lartg    0                     , (  1        j. __      )
0  0j_1 -: lartg    0                     , (  1        j.  _      )
0  0j1  -: lartg    0                     , (  FP_OVFL  j. __      )
0  0j_1 -: lartg    0                     , (  FP_OVFL  j.  _      )
0  1    -: lartg    0                     , (  _        j. -FP_OVFL)
0  1    -: lartg    0                     , (  _        j. _1      )
0  1    -: lartg    0                     , (  _        j. -FP_UNFL)
0  1    -: lartg    0                     , (  _        j.  FP_UNFL)
0  1    -: lartg    0                     , (  _        j.  1      )
0  1    -: lartg    0                     , (  _        j.  FP_OVFL)
1  0    -: lartg (  FP_UNFL  j. __      ) ,    0
1  0    -: lartg (  FP_UNFL  j.  _      ) ,    0
1  0    -: lartg (  1        j. __      ) ,    0
1  0    -: lartg (  1        j.  _      ) ,    0
1  0    -: lartg (  FP_OVFL  j. __      ) ,    0
1  0    -: lartg (  FP_OVFL  j.  _      ) ,    0
1  0    -: lartg (  _        j. -FP_OVFL) ,    0
1  0    -: lartg (  _        j. _1      ) ,    0
1  0    -: lartg (  _        j. -FP_UNFL) ,    0
1  0    -: lartg (  _        j.  FP_UNFL) ,    0
1  0    -: lartg (  _        j.  1      ) ,    0
1  0    -: lartg (  _        j.  FP_OVFL) ,    0
(0       1  ,: 1     0 ) -: lartg (0 , _ j. FP_OVFL) ,: (FP_OVFL j. __) , 0
((% %: 2 2) ,: 1     0 ) -: lartg q3                 ,: (FP_OVFL j. __) , 0
(0       1  ,: % %: 2 2) -: lartg (0 , _ j. FP_OVFL) ,: q4
NB. - input is imaginary infinity and consequently is
NB.   trivial
1  0    -: lartg __     0
1  0    -: lartg  0j__  0
0 _1    -: lartg  0    __
0  0j1  -: lartg  0     0j__
0  0j_1 -: lartg  0     0j_
0  1    -: lartg  0     _
1  0    -: lartg  0j_   0
1  0    -: lartg  _     0
NB. - input is trivial
 1    0                      -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   0
 1    0                      -: lartg ((-FP_OVFL) j. _1      ) ,   0
 1    0                      -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   0
 1    0                      -: lartg  (-FP_OVFL)              ,   0
 1    0                      -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   0
 1    0                      -: lartg ((-FP_OVFL) j.  1      ) ,   0
 1    0                      -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   0
 1    0                      -: lartg ((-sqr_big) j. -sqr_big) ,   0
 1    0                      -: lartg  (-sqr_big)              ,   0
 1    0                      -: lartg ((-sqr_big) j.  sqr_big) ,   0
 1    0                      -: lartg ( _1        j. -FP_OVFL) ,   0
 1    0                      -: lartg ( _1        j. _1      ) ,   0
 1    0                      -: lartg ( _1        j. -FP_UNFL) ,   0
 1    0                      -: lartg   _1                         0
 1    0                      -: lartg ( _1        j.  FP_UNFL) ,   0
 1    0                      -: lartg ( _1        j.  1      ) ,   0
 1    0                      -: lartg ( _1        j.  FP_OVFL) ,   0
 1    0                      -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   0
 1    0                      -: lartg ((-FP_UNFL) j. _1      ) ,   0
 1    0                      -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   0
 1    0                      -: lartg  (-FP_UNFL)              ,   0
 1    0                      -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   0
 1    0                      -: lartg ((-FP_UNFL) j.  1      ) ,   0
 1    0                      -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   0
 1    0                      -: lartg (  0        j. -FP_OVFL) ,   0
 1    0                      -: lartg (  0        j. -sqr_big) ,   0
 1    0                      -: lartg (  0        j. _1      ) ,   0
 1    0                      -: lartg (  0        j. -FP_UNFL) ,   0
(0 , _1j1  % %: 2          ) -: lartg    0                     , (-FP_OVFL) j. -FP_OVFL
(0 , _1    j.   FP_UNFL % 4) -: lartg    0                     , (-FP_OVFL) j. _1        NB. LAPACK's DLARTG: (c=0,s=-1,r=FP_OVFL), ZLARTG: (c=0,s=-1+i*FP_UNFL/4,r=FP_OVFL)
 0   _1                      -: lartg    0                     , (-FP_OVFL) j. -FP_UNFL
 0   _1                      -: lartg    0                     ,  -FP_OVFL
 0   _1                      -: lartg    0                     , (-FP_OVFL) j.  FP_UNFL
(0 , _1    j.   FP_UNFL %_4) -: lartg    0                     , (-FP_OVFL) j.  1        NB. LAPACK's DLARTG: (c=0,s=-1,r=FP_OVFL), ZLARTG: (c=0,s=-1-i*FP_UNFL/4,r=FP_OVFL)
(0 , _1j_1 % %: 2          ) -: lartg    0                     , (-FP_OVFL) j.  FP_OVFL
(0 , _1j1  % %: 2          ) -: lartg    0                     , (-sqr_big) j. -sqr_big
 0   _1                      -: lartg    0                     ,  -sqr_big
(0 , _1j_1 % %: 2          ) -: lartg    0                     , (-sqr_big) j.  sqr_big
(0 ,  1    j.~  FP_UNFL %_4) -: lartg    0                     ,  _1        j. -FP_OVFL
(0 , _1j1  % %: 2          ) -: lartg    0                     ,  _1        j. _1
(0 , _1    j.   FP_UNFL    ) -: lartg    0                     ,  _1        j. -FP_UNFL  NB. LAPACK's DLARTG: (c=0,s=-1,r=1), ZLARTG: (c=0,s=-1+i*FP_UNFL,r=1)
 0   _1                      -: lartg    0                        _1
(0 , _1    j.  -FP_UNFL    ) -: lartg    0                     ,  _1        j.  FP_UNFL  NB. LAPACK's DLARTG: (c=0,s=-1,r=1), ZLARTG: (c=0,s=-1-i*FP_UNFL,r=1)
(0 , _1j_1 % %: 2          ) -: lartg    0                     ,  _1        j.  1
(0 , _1    j.~  FP_UNFL %_4) -: lartg    0                     ,  _1        j.  FP_OVFL
 0    0j1                    -: lartg    0                     , (-FP_UNFL) j. -FP_OVFL
(0 ,  1    j.~ -FP_UNFL    ) -: lartg    0                     , (-FP_UNFL) j. _1
(0 , _1j1  % %: 2          ) -: lartg    0                     , (-FP_UNFL) j. -FP_UNFL
 0   _1                      -: lartg    0                     ,  -FP_UNFL
(0 , _1j_1 % %: 2          ) -: lartg    0                     , (-FP_UNFL) j.  FP_UNFL
(0 , _1    j.~ -FP_UNFL    ) -: lartg    0                     , (-FP_UNFL) j.  1
 0    0j_1                   -: lartg    0                     , (-FP_UNFL) j.  FP_OVFL
 0    0j1                    -: lartg    0                     ,   0        j. -FP_OVFL
 0    0j1                    -: lartg    0                     ,   0        j. -sqr_big
 0    0j1                    -: lartg    0                     ,   0        j. _1
 0    0j1                    -: lartg    0                     ,   0        j. -FP_UNFL
 1    0                      -: lartg    0                         0                     NB. test branching for zero input
 0    0j_1                   -: lartg    0                     ,   0        j.  FP_UNFL
 0    0j_1                   -: lartg    0                     ,   0        j.  1
 0    0j_1                   -: lartg    0                     ,   0        j.  sqr_big
 0    0j_1                   -: lartg    0                     ,   0        j.  FP_OVFL
 0    0j1                    -: lartg    0                     ,   FP_UNFL  j. -FP_OVFL
(0 ,  1    j.~  FP_UNFL    ) -: lartg    0                     ,   FP_UNFL  j. _1
(0 ,  1j1   % %: 2         ) -: lartg    0                     ,   FP_UNFL  j. -FP_UNFL
 0    1                      -: lartg    0                     ,   FP_UNFL
(0 ,  1j_1  % %: 2         ) -: lartg    0                     ,   FP_UNFL  j.  FP_UNFL
(0 , _1    j.~  FP_UNFL    ) -: lartg    0                     ,   FP_UNFL  j.  1
 0    0j_1                   -: lartg    0                     ,   FP_UNFL  j.  FP_OVFL
(0 ,  1    j.~  FP_UNFL % 4) -: lartg    0                     ,   1        j. -FP_OVFL
(0 ,  1j1   % %: 2         ) -: lartg    0                     ,   1        j. _1
(0 ,  1    j.   FP_UNFL    ) -: lartg    0                     ,   1        j. -FP_UNFL  NB. LAPACK's DLARTG: (c=0,s= 1,r=1), ZLARTG: (c=0,s= 1+i*FP_UNFL,r=1)
 0    1                      -: lartg    0                         1
(0 ,  1    j.  -FP_UNFL    ) -: lartg    0                     ,   1        j.  FP_UNFL  NB. LAPACK's DLARTG: (c=0,s= 1,r=1), ZLARTG: (c=0,s= 1-i*FP_UNFL,r=1)
(0 ,  1j_1  % %: 2         ) -: lartg    0                     ,   1        j.  1
(0 , _1    j.~  FP_UNFL % 4) -: lartg    0                     ,   1        j.  FP_OVFL
 0    1                      -: lartg    0                     ,   sqr_big
(0 ,  1j1   % %: 2         ) -: lartg    0                     ,   FP_OVFL  j. -FP_OVFL
(0 ,  1    j.   FP_UNFL % 4) -: lartg    0                     ,   FP_OVFL  j. _1        NB. LAPACK's DLARTG: (c=0,s= 1,r=FP_OVFL), ZLARTG: (c=0,s= 1+i*FP_UNFL/4,r=FP_OVFL)
 0    1                      -: lartg    0                     ,   FP_OVFL  j. -FP_UNFL
 0    1                      -: lartg    0                     ,   FP_OVFL
 0    1                      -: lartg    0                     ,   FP_OVFL  j.  FP_UNFL
(0 ,  1    j.   FP_UNFL %_4) -: lartg    0                     ,   FP_OVFL  j.  1        NB. LAPACK's DLARTG: (c=0,s= 1,r=FP_OVFL), ZLARTG: (c=0,s= 1-i*FP_UNFL/4,r=FP_OVFL)
(0 ,  1j_1  % %: 2         ) -: lartg    0                     ,   FP_OVFL  j.  FP_OVFL
 1    0                      -: lartg (  0        j.  FP_UNFL) ,   0
 1    0                      -: lartg (  0        j.  1      ) ,   0
 1    0                      -: lartg (  0        j.  FP_OVFL) ,   0
 1    0                      -: lartg (  FP_UNFL  j. -FP_OVFL) ,   0
 1    0                      -: lartg (  FP_UNFL  j. _1      ) ,   0
 1    0                      -: lartg (  FP_UNFL  j. -FP_UNFL) ,   0
 1    0                      -: lartg    FP_UNFL               ,   0
 1    0                      -: lartg (  FP_UNFL  j.  FP_UNFL) ,   0
 1    0                      -: lartg (  FP_UNFL  j.  1      ) ,   0
 1    0                      -: lartg (  FP_UNFL  j.  FP_OVFL) ,   0
 1    0                      -: lartg (  1        j. -FP_OVFL) ,   0
 1    0                      -: lartg (  1        j. _1      ) ,   0
 1    0                      -: lartg (  1        j. -FP_UNFL) ,   0
 1    0                      -: lartg    1                         0
 1    0                      -: lartg (  1        j.  FP_UNFL) ,   0
 1    0                      -: lartg (  1        j.  1      ) ,   0
 1    0                      -: lartg (  1        j.  FP_OVFL) ,   0
 1    0                      -: lartg (  sqr_big  j. -sqr_big) ,   0
 1    0                      -: lartg    sqr_big               ,   0
 1    0                      -: lartg (  sqr_big  j.  sqr_big) ,   0
 1    0                      -: lartg (  FP_OVFL  j. -FP_OVFL) ,   0
 1    0                      -: lartg (  FP_OVFL  j. _1      ) ,   0
 1    0                      -: lartg (  FP_OVFL  j. -FP_UNFL) ,   0
 1    0                      -: lartg    FP_OVFL               ,   0
 1    0                      -: lartg (  FP_OVFL  j.  FP_UNFL) ,   0
 1    0                      -: lartg (  FP_OVFL  j.  1      ) ,   0
 1    0                      -: lartg (  FP_OVFL  j.  FP_OVFL) ,   0
NB. - edge cases input
(1  1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(2  1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(2  1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,  -FP_OVFL
(2  1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,  -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   0        j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   0        j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  FP_UNFL
(2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  FP_OVFL
(2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL
(2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(1 _1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(1  1j_1 % %: 3) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(1  1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(1  1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,  -FP_OVFL
(1  1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(1  1j1  % %: 3) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,  -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   0        j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   0        j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   0        j.  FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   0        j.  FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(1 _1j_1 % %: 3) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(1 _1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(1 _1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL
(1 _1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(1 _1j1  % %: 3) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(1  1j_1 % %: 3) -: lartg  (-FP_OVFL)              , (-FP_OVFL) j. -FP_OVFL
(1  1    % %: 2) -: lartg  (-FP_OVFL)              , (-FP_OVFL) j. -FP_UNFL
(1  1    % %: 2) -: lartg  (-FP_OVFL)              ,  -FP_OVFL
(1  1    % %: 2) -: lartg  (-FP_OVFL)              , (-FP_OVFL) j.  FP_UNFL
(1  1j1  % %: 3) -: lartg  (-FP_OVFL)              , (-FP_OVFL) j.  FP_OVFL
(1 ,   %FP_OVFL) -: lartg  (-FP_OVFL)              ,  _1
(1  0j_1 % %: 2) -: lartg  (-FP_OVFL)              , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg  (-FP_OVFL)              , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg  (-FP_OVFL)              ,  -FP_UNFL
 1  0            -: lartg  (-FP_OVFL)              , (-FP_UNFL) j.  FP_UNFL
(1  0j1  % %: 2) -: lartg  (-FP_OVFL)              , (-FP_UNFL) j.  FP_OVFL
(1  0j_1 % %: 2) -: lartg  (-FP_OVFL)              ,   0        j. -FP_OVFL
 1  0            -: lartg  (-FP_OVFL)              ,   0        j. -FP_UNFL
 1  0            -: lartg  (-FP_OVFL)              ,   0        j.  FP_UNFL
(1  0j1  % %: 2) -: lartg  (-FP_OVFL)              ,   0        j.  FP_OVFL
(1  0j_1 % %: 2) -: lartg  (-FP_OVFL)              ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg  (-FP_OVFL)              ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg  (-FP_OVFL)              ,   FP_UNFL
 1  0            -: lartg  (-FP_OVFL)              ,   FP_UNFL  j.  FP_UNFL
(1  0j1  % %: 2) -: lartg  (-FP_OVFL)              ,   FP_UNFL  j.  FP_OVFL
(1 , - %FP_OVFL) -: lartg  (-FP_OVFL)              ,   1
(1 _1j_1 % %: 3) -: lartg  (-FP_OVFL)              ,   FP_OVFL  j. -FP_OVFL
(1 _1    % %: 2) -: lartg  (-FP_OVFL)              ,   FP_OVFL  j. -FP_UNFL
(1 _1    % %: 2) -: lartg  (-FP_OVFL)              ,   FP_OVFL
(1 _1    % %: 2) -: lartg  (-FP_OVFL)              ,   FP_OVFL  j.  FP_UNFL
(1 _1j1  % %: 3) -: lartg  (-FP_OVFL)              ,   FP_OVFL  j.  FP_OVFL
(1  1j_1 % %: 3) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(1  1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(1  1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,  -FP_OVFL
(1  1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(1  1j1  % %: 3) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,  -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   0        j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   0        j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   0        j.  FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   0        j.  FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(1 _1j_1 % %: 3) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(1 _1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(1 _1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL
(1 _1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(1 _1j1  % %: 3) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,  -FP_OVFL
(2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(1  1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,  -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(2  1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   0        j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   0        j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  FP_UNFL
(2  1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  FP_OVFL
(2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL
 1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(2  1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(1 _1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL
(2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(1  1    % %: 2) -: lartg ((-sqr_big) j. -sqr_big) , (-sqr_big) j. -sqr_big
(1  0j1  % %: 2) -: lartg ((-sqr_big) j. -sqr_big) , (-sqr_big) j.  sqr_big
(1  0j_1 % %: 2) -: lartg ((-sqr_big) j. -sqr_big) ,   sqr_big  j. -sqr_big
(1 _1    % %: 2) -: lartg ((-sqr_big) j. -sqr_big) ,   sqr_big  j.  sqr_big
(1  1    % %: 2) -: lartg  (-sqr_big)              ,  -sqr_big
(1  0j_1 % %: 2) -: lartg  (-sqr_big)              ,   0        j. -sqr_big
(1  0j1  % %: 2) -: lartg  (-sqr_big)              ,   0        j.  sqr_big
(1 _1    % %: 2) -: lartg  (-sqr_big)              ,   sqr_big
(1  0j_1 % %: 2) -: lartg ((-sqr_big) j.  sqr_big) , (-sqr_big) j. -sqr_big
(1  1    % %: 2) -: lartg ((-sqr_big) j.  sqr_big) , (-sqr_big) j.  sqr_big
(1 _1    % %: 2) -: lartg ((-sqr_big) j.  sqr_big) ,   sqr_big  j. -sqr_big
(1  0j1  % %: 2) -: lartg ((-sqr_big) j.  sqr_big) ,   sqr_big  j.  sqr_big
(1 ,~  %FP_OVFL) -: lartg   _1                     ,  -FP_OVFL
(1 ,    FP_UNFL) -: lartg   _1                     ,  -FP_UNFL
(1 , -  FP_UNFL) -: lartg   _1                     ,   FP_UNFL
(_1,~  %FP_OVFL) -: lartg   _1                     ,   FP_OVFL
(1  1j1  % %: 3) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(1  0j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,  -FP_OVFL
(1  0j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(1 _1j1  % %: 3) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(1  1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,  -FP_UNFL
 1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(1 _1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(1  1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   0        j. -FP_OVFL
 1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   0        j. -FP_UNFL
 1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   0        j.  FP_UNFL
(1 _1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   0        j.  FP_OVFL
(1  1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL
 1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(1 _1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(1  1j_1 % %: 3) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(1 _1j_1 % %: 3) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
 0  1            -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(0  1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(0  1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,  -FP_OVFL
(0  1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 0  0j1          -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(1  1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(2  1j1  % %: 6) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,  -FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   0        j. -FP_OVFL
(2  1j_1 % %: 6) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   0        j. -FP_UNFL
(2 _1j1  % %: 6) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   0        j.  FP_UNFL
(0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   0        j.  FP_OVFL
(0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(2 _1j_1 % %: 6) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL
(1 _1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 0  0j_1         -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL
(0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 0 _1            -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(0  1j_1 % %: 2) -: lartg  (-FP_UNFL)              , (-FP_OVFL) j. -FP_OVFL
 0  1            -: lartg  (-FP_UNFL)              , (-FP_OVFL) j. -FP_UNFL
 0  1            -: lartg  (-FP_UNFL)              ,  -FP_OVFL
 0  1            -: lartg  (-FP_UNFL)              , (-FP_OVFL) j.  FP_UNFL
(0  1j1  % %: 2) -: lartg  (-FP_UNFL)              , (-FP_OVFL) j.  FP_OVFL
(1 ,~   FP_UNFL) -: lartg  (-FP_UNFL)              ,  _1
 0  0j_1         -: lartg  (-FP_UNFL)              , (-FP_UNFL) j. -FP_OVFL
(1  1j_1 % %: 3) -: lartg  (-FP_UNFL)              , (-FP_UNFL) j. -FP_UNFL
(1  1    % %: 2) -: lartg  (-FP_UNFL)              ,  -FP_UNFL
(1  1j1  % %: 3) -: lartg  (-FP_UNFL)              , (-FP_UNFL) j.  FP_UNFL
 0  0j1          -: lartg  (-FP_UNFL)              , (-FP_UNFL) j.  FP_OVFL
 0  0j_1         -: lartg  (-FP_UNFL)              ,   0        j. -FP_OVFL
(1  0j_1 % %: 2) -: lartg  (-FP_UNFL)              ,   0        j. -FP_UNFL
(1  0j1  % %: 2) -: lartg  (-FP_UNFL)              ,   0        j.  FP_UNFL
 0  0j1          -: lartg  (-FP_UNFL)              ,   0        j.  FP_OVFL
 0  0j_1         -: lartg  (-FP_UNFL)              ,   FP_UNFL  j. -FP_OVFL
(1 _1j_1 % %: 3) -: lartg  (-FP_UNFL)              ,   FP_UNFL  j. -FP_UNFL
(1 _1    % %: 2) -: lartg  (-FP_UNFL)              ,   FP_UNFL
(1 _1j1  % %: 3) -: lartg  (-FP_UNFL)              ,   FP_UNFL  j.  FP_UNFL
 0  0j1          -: lartg  (-FP_UNFL)              ,   FP_UNFL  j.  FP_OVFL
(_1,~   FP_UNFL) -: lartg  (-FP_UNFL)              ,   1
(0 _1j_1 % %: 2) -: lartg  (-FP_UNFL)              ,   FP_OVFL  j. -FP_OVFL
 0 _1            -: lartg  (-FP_UNFL)              ,   FP_OVFL  j. -FP_UNFL
 0 _1            -: lartg  (-FP_UNFL)              ,   FP_OVFL
 0 _1            -: lartg  (-FP_UNFL)              ,   FP_OVFL  j.  FP_UNFL
(0 _1j1  % %: 2) -: lartg  (-FP_UNFL)              ,   FP_OVFL  j.  FP_OVFL
 0  0j_1         -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,  -FP_OVFL
(0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 0  1            -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(2  1j_1 % %: 6) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,  -FP_UNFL
(1  1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(0  1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   0        j. -FP_OVFL
(2 _1j_1 % %: 6) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   0        j. -FP_UNFL
(2  1j1  % %: 6) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   0        j.  FP_UNFL
(0  1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   0        j.  FP_OVFL
(0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(1 _1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(2 _1j1  % %: 6) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(0  1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 0 _1            -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL
(0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 0  0j1          -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(1 _1j_1 % %: 3) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,  -FP_OVFL
(1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(1  1j_1 % %: 3) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(1 _1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,  -FP_UNFL
 1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(1  1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(1 _1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   0        j. -FP_OVFL
 1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   0        j. -FP_UNFL
 1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   0        j.  FP_UNFL
(1  1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   0        j.  FP_OVFL
(1 _1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL
 1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(1  1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(1 _1j1  % %: 3) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(1  0j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(1  0j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL
(1  0j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(1  1j1  % %: 3) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(1  1j1  % %: 3) -: lartg (  0        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(1  0j1  % %: 2) -: lartg (  0        j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(1  0j1  % %: 2) -: lartg (  0        j. -FP_OVFL) ,  -FP_OVFL
(1  0j1  % %: 2) -: lartg (  0        j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(1 _1j1  % %: 3) -: lartg (  0        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(1  1    % %: 2) -: lartg (  0        j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg (  0        j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg (  0        j. -FP_OVFL) ,  -FP_UNFL
 1  0            -: lartg (  0        j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(1 _1    % %: 2) -: lartg (  0        j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(1  1    % %: 2) -: lartg (  0        j. -FP_OVFL) ,   0        j. -FP_OVFL
 1  0            -: lartg (  0        j. -FP_OVFL) ,   0        j. -FP_UNFL
 1  0            -: lartg (  0        j. -FP_OVFL) ,   0        j.  FP_UNFL
(1 _1    % %: 2) -: lartg (  0        j. -FP_OVFL) ,   0        j.  FP_OVFL
(1  1    % %: 2) -: lartg (  0        j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg (  0        j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg (  0        j. -FP_OVFL) ,   FP_UNFL
 1  0            -: lartg (  0        j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(1 _1    % %: 2) -: lartg (  0        j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(1  1j_1 % %: 3) -: lartg (  0        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(1  0j_1 % %: 2) -: lartg (  0        j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(1  0j_1 % %: 2) -: lartg (  0        j. -FP_OVFL) ,   FP_OVFL
(1  0j_1 % %: 2) -: lartg (  0        j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(1 _1j_1 % %: 3) -: lartg (  0        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(1  0j1  % %: 2) -: lartg (  0        j. -sqr_big) ,  -sqr_big
(1  1    % %: 2) -: lartg (  0        j. -sqr_big) ,   0        j. -sqr_big
(1 _1    % %: 2) -: lartg (  0        j. -sqr_big) ,   0        j.  sqr_big
(1  0j_1 % %: 2) -: lartg (  0        j. -sqr_big) ,   sqr_big
(0  1j1  % %: 2) -: lartg (  0        j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 0  0j1          -: lartg (  0        j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 0  0j1          -: lartg (  0        j. -FP_UNFL) ,  -FP_OVFL
 0  0j1          -: lartg (  0        j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(0 _1j1  % %: 2) -: lartg (  0        j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 0  1            -: lartg (  0        j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(1  1j1  % %: 3) -: lartg (  0        j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(1  0j1  % %: 2) -: lartg (  0        j. -FP_UNFL) ,  -FP_UNFL
(1 _1j1  % %: 3) -: lartg (  0        j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 0 _1            -: lartg (  0        j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 0  1            -: lartg (  0        j. -FP_UNFL) ,   0        j. -FP_OVFL
(1  1    % %: 2) -: lartg (  0        j. -FP_UNFL) ,   0        j. -FP_UNFL
(1 _1    % %: 2) -: lartg (  0        j. -FP_UNFL) ,   0        j.  FP_UNFL
 0 _1            -: lartg (  0        j. -FP_UNFL) ,   0        j.  FP_OVFL
 0  1            -: lartg (  0        j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(1  1j_1 % %: 3) -: lartg (  0        j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(1  0j_1 % %: 2) -: lartg (  0        j. -FP_UNFL) ,   FP_UNFL
(1 _1j_1 % %: 3) -: lartg (  0        j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 0 _1            -: lartg (  0        j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(0  1j_1 % %: 2) -: lartg (  0        j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 0  0j_1         -: lartg (  0        j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 0  0j_1         -: lartg (  0        j. -FP_UNFL) ,   FP_OVFL
 0  0j_1         -: lartg (  0        j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(0 _1j_1 % %: 2) -: lartg (  0        j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(0 _1j_1 % %: 2) -: lartg (  0        j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
 0  0j_1         -: lartg (  0        j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
 0  0j_1         -: lartg (  0        j.  FP_UNFL) ,  -FP_OVFL
 0  0j_1         -: lartg (  0        j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(0  1j_1 % %: 2) -: lartg (  0        j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
 0 _1            -: lartg (  0        j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(1 _1j_1 % %: 3) -: lartg (  0        j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(1  0j_1 % %: 2) -: lartg (  0        j.  FP_UNFL) ,  -FP_UNFL
(1  1j_1 % %: 3) -: lartg (  0        j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
 0  1            -: lartg (  0        j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
 0 _1            -: lartg (  0        j.  FP_UNFL) ,   0        j. -FP_OVFL
(1 _1    % %: 2) -: lartg (  0        j.  FP_UNFL) ,   0        j. -FP_UNFL
(1  1    % %: 2) -: lartg (  0        j.  FP_UNFL) ,   0        j.  FP_UNFL
 0  1            -: lartg (  0        j.  FP_UNFL) ,   0        j.  FP_OVFL
 0 _1            -: lartg (  0        j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(1 _1j1  % %: 3) -: lartg (  0        j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(1  0j1  % %: 2) -: lartg (  0        j.  FP_UNFL) ,   FP_UNFL
(1  1j1  % %: 3) -: lartg (  0        j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
 0  1            -: lartg (  0        j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(0 _1j1  % %: 2) -: lartg (  0        j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
 0  0j1          -: lartg (  0        j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
 0  0j1          -: lartg (  0        j.  FP_UNFL) ,   FP_OVFL
 0  0j1          -: lartg (  0        j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(0  1j1  % %: 2) -: lartg (  0        j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(1  0j_1 % %: 2) -: lartg (  0        j.  sqr_big) ,  -sqr_big
(1 _1    % %: 2) -: lartg (  0        j.  sqr_big) ,   0        j. -sqr_big
(1  1    % %: 2) -: lartg (  0        j.  sqr_big) ,   0        j.  sqr_big
(1  0j1  % %: 2) -: lartg (  0        j.  sqr_big) ,   sqr_big
(1 _1j_1 % %: 3) -: lartg (  0        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(1  0j_1 % %: 2) -: lartg (  0        j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(1  0j_1 % %: 2) -: lartg (  0        j.  FP_OVFL) ,  -FP_OVFL
(1  0j_1 % %: 2) -: lartg (  0        j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(1  1j_1 % %: 3) -: lartg (  0        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(1 _1    % %: 2) -: lartg (  0        j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg (  0        j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg (  0        j.  FP_OVFL) ,  -FP_UNFL
 1  0            -: lartg (  0        j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(1  1    % %: 2) -: lartg (  0        j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(1 _1    % %: 2) -: lartg (  0        j.  FP_OVFL) ,   0        j. -FP_OVFL
 1  0            -: lartg (  0        j.  FP_OVFL) ,   0        j. -FP_UNFL
 1  0            -: lartg (  0        j.  FP_OVFL) ,   0        j.  FP_UNFL
(1  1    % %: 2) -: lartg (  0        j.  FP_OVFL) ,   0        j.  FP_OVFL
(1 _1    % %: 2) -: lartg (  0        j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg (  0        j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg (  0        j.  FP_OVFL) ,   FP_UNFL
 1  0            -: lartg (  0        j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(1  1    % %: 2) -: lartg (  0        j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(1 _1j1  % %: 3) -: lartg (  0        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(1  0j1  % %: 2) -: lartg (  0        j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(1  0j1  % %: 2) -: lartg (  0        j.  FP_OVFL) ,   FP_OVFL
(1  0j1  % %: 2) -: lartg (  0        j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(1  1j1  % %: 3) -: lartg (  0        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(1  1j1  % %: 3) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(1  0j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,  -FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(1 _1j1  % %: 3) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(1  1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,  -FP_UNFL
 1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(1 _1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(1  1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   0        j. -FP_OVFL
 1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,   0        j. -FP_UNFL
 1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,   0        j.  FP_UNFL
(1 _1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   0        j.  FP_OVFL
(1  1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL
 1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(1 _1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(1  1j_1 % %: 3) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(1  0j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL
(1  0j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(1 _1j_1 % %: 3) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
 0  0j1          -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(0 _1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(0 _1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,  -FP_OVFL
(0 _1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 0 _1            -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(0  1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(2 _1j1  % %: 6) -: lartg (  FP_UNFL  j. -FP_UNFL) ,  -FP_UNFL
(1 _1    % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(0  1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   0        j. -FP_OVFL
(2  1j1  % %: 6) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   0        j. -FP_UNFL
(2 _1j_1 % %: 6) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   0        j.  FP_UNFL
(0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   0        j.  FP_OVFL
(0  1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(1  1    % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(2  1j_1 % %: 6) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 0  1            -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(0  1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(0  1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL
(0  1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 0  0j_1         -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(0 _1j1  % %: 2) -: lartg    FP_UNFL               , (-FP_OVFL) j. -FP_OVFL
 0 _1            -: lartg    FP_UNFL               , (-FP_OVFL) j. -FP_UNFL
 0 _1            -: lartg    FP_UNFL               ,  -FP_OVFL
 0 _1            -: lartg    FP_UNFL               , (-FP_OVFL) j.  FP_UNFL
(0 _1j_1 % %: 2) -: lartg    FP_UNFL               , (-FP_OVFL) j.  FP_OVFL
(_1,~   FP_UNFL) -: lartg    FP_UNFL               ,  _1
 0  0j1          -: lartg    FP_UNFL               , (-FP_UNFL) j. -FP_OVFL
(1 _1j1  % %: 3) -: lartg    FP_UNFL               , (-FP_UNFL) j. -FP_UNFL
(1 _1    % %: 2) -: lartg    FP_UNFL               ,  -FP_UNFL
(1 _1j_1 % %: 3) -: lartg    FP_UNFL               , (-FP_UNFL) j.  FP_UNFL
 0  0j_1         -: lartg    FP_UNFL               , (-FP_UNFL) j.  FP_OVFL
 0  0j1          -: lartg    FP_UNFL               ,   0        j. -FP_OVFL
(1  0j1  % %: 2) -: lartg    FP_UNFL               ,   0        j. -FP_UNFL
(1  0j_1 % %: 2) -: lartg    FP_UNFL               ,   0        j.  FP_UNFL
 0  0j_1         -: lartg    FP_UNFL               ,   0        j.  FP_OVFL
 0  0j1          -: lartg    FP_UNFL               ,   FP_UNFL  j. -FP_OVFL
(1  1j1  % %: 3) -: lartg    FP_UNFL               ,   FP_UNFL  j. -FP_UNFL
(1  1    % %: 2) -: lartg    FP_UNFL               ,   FP_UNFL
(1  1j_1 % %: 3) -: lartg    FP_UNFL               ,   FP_UNFL  j.  FP_UNFL
 0  0j_1         -: lartg    FP_UNFL               ,   FP_UNFL  j.  FP_OVFL
(1 ,~   FP_UNFL) -: lartg    FP_UNFL               ,   1
(0  1j1  % %: 2) -: lartg    FP_UNFL               ,   FP_OVFL  j. -FP_OVFL
 0  1            -: lartg    FP_UNFL               ,   FP_OVFL  j. -FP_UNFL
 0  1            -: lartg    FP_UNFL               ,   FP_OVFL
 0  1            -: lartg    FP_UNFL               ,   FP_OVFL  j.  FP_UNFL
(0  1j_1 % %: 2) -: lartg    FP_UNFL               ,   FP_OVFL  j.  FP_OVFL
 0 _1            -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,  -FP_OVFL
(0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
 0  0j_1         -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(0 _1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
(1 _1    % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
(2 _1j_1 % %: 6) -: lartg (  FP_UNFL  j.  FP_UNFL) ,  -FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(0  1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(0 _1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   0        j. -FP_OVFL
(2 _1j1  % %: 6) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   0        j. -FP_UNFL
(2  1j_1 % %: 6) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   0        j.  FP_UNFL
(0  1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   0        j.  FP_OVFL
(0 _1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
(2  1j1  % %: 6) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL
(1  1    % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(0  1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
 0  0j1          -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(0  1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(0  1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL
(0  1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
 0  1            -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(1 _1j_1 % %: 3) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(1  0j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,  -FP_OVFL
(1  0j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(1  1j_1 % %: 3) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(1 _1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,  -FP_UNFL
 1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(1  1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(1 _1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   0        j. -FP_OVFL
 1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,   0        j. -FP_UNFL
 1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,   0        j.  FP_UNFL
(1  1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   0        j.  FP_OVFL
(1 _1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL
 1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(1  1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(1 _1j1  % %: 3) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(1  0j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(1  1j1  % %: 3) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(_1,~  %FP_OVFL) -: lartg    1                     ,  -FP_OVFL
( 1, -  FP_UNFL) -: lartg    1                     ,  -FP_UNFL
( 1,    FP_UNFL) -: lartg    1                     ,   FP_UNFL
( 1,~  %FP_OVFL) -: lartg    1                     ,   FP_OVFL
(1  0j1  % %: 2) -: lartg (  sqr_big  j. -sqr_big) , (-sqr_big) j. -sqr_big
(1 _1    % %: 2) -: lartg (  sqr_big  j. -sqr_big) , (-sqr_big) j.  sqr_big
(1  1    % %: 2) -: lartg (  sqr_big  j. -sqr_big) ,   sqr_big  j. -sqr_big
(1  0j_1 % %: 2) -: lartg (  sqr_big  j. -sqr_big) ,   sqr_big  j.  sqr_big
(1 _1    % %: 2) -: lartg    sqr_big               ,  -sqr_big
(1  0j1  % %: 2) -: lartg    sqr_big               ,   0        j. -sqr_big
(1  0j_1 % %: 2) -: lartg    sqr_big               ,   0        j.  sqr_big
(1  1    % %: 2) -: lartg    sqr_big               ,   sqr_big
(1 _1    % %: 2) -: lartg (  sqr_big  j.  sqr_big) , (-sqr_big) j. -sqr_big
(1  0j_1 % %: 2) -: lartg (  sqr_big  j.  sqr_big) , (-sqr_big) j.  sqr_big
(1  0j1  % %: 2) -: lartg (  sqr_big  j.  sqr_big) ,   sqr_big  j. -sqr_big
(1  1    % %: 2) -: lartg (  sqr_big  j.  sqr_big) ,   sqr_big  j.  sqr_big
(1  0j1  % %: 2) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(2 _1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(2 _1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,  -FP_OVFL
(2 _1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(1 _1    % %: 2) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(2  1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,  -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(2  1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   0        j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,   0        j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,   0        j.  FP_UNFL
(2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   0        j.  FP_OVFL
(2  1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL
 1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(1  1    % %: 2) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(2  1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(2  1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL
(2  1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
(1 _1j1  % %: 3) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(1 _1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(1 _1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,  -FP_OVFL
(1 _1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(1 _1j_1 % %: 3) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,  -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   0        j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,   0        j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,   0        j.  FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   0        j.  FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL
 1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(1  1j1  % %: 3) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(1  1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(1  1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL
(1  1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(1  1j_1 % %: 3) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(1 _1j1  % %: 3) -: lartg    FP_OVFL               , (-FP_OVFL) j. -FP_OVFL
(1 _1    % %: 2) -: lartg    FP_OVFL               , (-FP_OVFL) j. -FP_UNFL
(1 _1    % %: 2) -: lartg    FP_OVFL               ,  -FP_OVFL
(1 _1    % %: 2) -: lartg    FP_OVFL               , (-FP_OVFL) j.  FP_UNFL
(1 _1j_1 % %: 3) -: lartg    FP_OVFL               , (-FP_OVFL) j.  FP_OVFL
(1 , - %FP_OVFL) -: lartg    FP_OVFL               ,  _1
(1  0j1  % %: 2) -: lartg    FP_OVFL               , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg    FP_OVFL               , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg    FP_OVFL               ,  -FP_UNFL
 1  0            -: lartg    FP_OVFL               , (-FP_UNFL) j.  FP_UNFL
(1  0j_1 % %: 2) -: lartg    FP_OVFL               , (-FP_UNFL) j.  FP_OVFL
(1  0j1  % %: 2) -: lartg    FP_OVFL               ,   0        j. -FP_OVFL
 1  0            -: lartg    FP_OVFL               ,   0        j. -FP_UNFL
 1  0            -: lartg    FP_OVFL               ,   0        j.  FP_UNFL
(1  0j_1 % %: 2) -: lartg    FP_OVFL               ,   0        j.  FP_OVFL
(1  0j1  % %: 2) -: lartg    FP_OVFL               ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg    FP_OVFL               ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg    FP_OVFL               ,   FP_UNFL
 1  0            -: lartg    FP_OVFL               ,   FP_UNFL  j.  FP_UNFL
(1  0j_1 % %: 2) -: lartg    FP_OVFL               ,   FP_UNFL  j.  FP_OVFL
(1 ,   %FP_OVFL) -: lartg    FP_OVFL               ,   1
(1  1j1  % %: 3) -: lartg    FP_OVFL               ,   FP_OVFL  j. -FP_OVFL
(1  1    % %: 2) -: lartg    FP_OVFL               ,   FP_OVFL  j. -FP_UNFL
(1  1    % %: 2) -: lartg    FP_OVFL               ,   FP_OVFL
(1  1    % %: 2) -: lartg    FP_OVFL               ,   FP_OVFL  j.  FP_UNFL
(1  1j_1 % %: 3) -: lartg    FP_OVFL               ,   FP_OVFL  j.  FP_OVFL
(1 _1j1  % %: 3) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
(1 _1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
(1 _1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,  -FP_OVFL
(1 _1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
(1 _1j_1 % %: 3) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,  -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   0        j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,   0        j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,   0        j.  FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   0        j.  FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL
 1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
(1  1j1  % %: 3) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
(1  1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
(1  1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL
(1  1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
(1  1j_1 % %: 3) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
(1 _1    % %: 2) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
(2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
(2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,  -FP_OVFL
(2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
(1  0j_1 % %: 2) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
(2 _1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,  -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
(2  1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
(2 _1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   0        j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,   0        j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,   0        j.  FP_UNFL
(2  1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   0        j.  FP_OVFL
(2 _1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
 1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
 1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL
 1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
(2  1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
(1  0j1  % %: 2) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
(2  1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
(2  1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL
(2  1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
(1  1    % %: 2) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
NB. - input without edge cases
(1  1    % %: 2) -: lartg _1j_1 _1j_1
(2  1j1  % %: 6) -: lartg _1j_1 _1
(1  0j1  % %: 2) -: lartg _1j_1 _1j1
(2  1j_1 % %: 6) -: lartg _1j_1  0j_1
(2 _1j1  % %: 6) -: lartg _1j_1  0j1
(1  0j_1 % %: 2) -: lartg _1j_1  1j_1
(2 _1j_1 % %: 6) -: lartg _1j_1  1
(1 _1    % %: 2) -: lartg _1j_1  1j1
(1  1j_1 % %: 3) -: lartg _1    _1j_1
(1  1    % %: 2) -: lartg _1    _1
(1  1j1  % %: 3) -: lartg _1    _1j1
(1  0j_1 % %: 2) -: lartg _1     0j_1
(1  0j1  % %: 2) -: lartg _1     0j1
(1 _1j_1 % %: 3) -: lartg _1     1j_1
(1 _1    % %: 2) -: lartg _1     1
(1 _1j1  % %: 3) -: lartg _1     1j1
(1  0j_1 % %: 2) -: lartg _1j1  _1j_1
(2  1j_1 % %: 6) -: lartg _1j1  _1
(1  1    % %: 2) -: lartg _1j1  _1j1
(2 _1j_1 % %: 6) -: lartg _1j1   0j_1
(2  1j1  % %: 6) -: lartg _1j1   0j1
(1 _1    % %: 2) -: lartg _1j1   1j_1
(2 _1j1  % %: 6) -: lartg _1j1   1
(1  0j1  % %: 2) -: lartg _1j1   1j1
(1  1j1  % %: 3) -: lartg  0j_1 _1j_1
(1  0j1  % %: 2) -: lartg  0j_1 _1
(1 _1j1  % %: 3) -: lartg  0j_1 _1j1
(1  1    % %: 2) -: lartg  0j_1  0j_1
(1 _1    % %: 2) -: lartg  0j_1  0j1
(1  1j_1 % %: 3) -: lartg  0j_1  1j_1
(1  0j_1 % %: 2) -: lartg  0j_1  1
(1 _1j_1 % %: 3) -: lartg  0j_1  1j1
(1 _1j_1 % %: 3) -: lartg  0j1  _1j_1
(1  0j_1 % %: 2) -: lartg  0j1  _1
(1  1j_1 % %: 3) -: lartg  0j1  _1j1
(1 _1    % %: 2) -: lartg  0j1   0j_1
(1  1    % %: 2) -: lartg  0j1   0j1
(1 _1j1  % %: 3) -: lartg  0j1   1j_1
(1  0j1  % %: 2) -: lartg  0j1   1
(1  1j1  % %: 3) -: lartg  0j1   1j1
(1  0j1  % %: 2) -: lartg  1j_1 _1j_1
(2 _1j1  % %: 6) -: lartg  1j_1 _1
(1 _1    % %: 2) -: lartg  1j_1 _1j1
(2  1j1  % %: 6) -: lartg  1j_1  0j_1
(2 _1j_1 % %: 6) -: lartg  1j_1  0j1
(1  1    % %: 2) -: lartg  1j_1  1j_1
(2  1j_1 % %: 6) -: lartg  1j_1  1
(1  0j_1 % %: 2) -: lartg  1j_1  1j1
(1 _1j1  % %: 3) -: lartg  1    _1j_1
(1 _1    % %: 2) -: lartg  1    _1
(1 _1j_1 % %: 3) -: lartg  1    _1j1
(1  0j1  % %: 2) -: lartg  1     0j_1
(1  0j_1 % %: 2) -: lartg  1     0j1
(1  1j1  % %: 3) -: lartg  1     1j_1
(1  1    % %: 2) -: lartg  1     1
(1  1j_1 % %: 3) -: lartg  1     1j1
(1 _1    % %: 2) -: lartg  1j1  _1j_1
(2 _1j_1 % %: 6) -: lartg  1j1  _1
(1  0j_1 % %: 2) -: lartg  1j1  _1j1
(2 _1j1  % %: 6) -: lartg  1j1   0j_1
(2  1j_1 % %: 6) -: lartg  1j1   0j1
(1  0j1  % %: 2) -: lartg  1j1   1j_1
(2  1j1  % %: 6) -: lartg  1j1   1
(1  1    % %: 2) -: lartg  1j1   1j1
(2 2 $   % %: 2) -: lartg q3 ,: q4
NB. c is real
0 -:!.0 qni lartg _1j_1 _1j_1
0 -:!.0 qni lartg _1j_1 _1j1
0 -:!.0 qni lartg _1j_1  1j_1
0 -:!.0 qni lartg _1j_1  1j1
0 -:!.0 qni lartg _1j1  _1j_1
0 -:!.0 qni lartg _1j1  _1j1
0 -:!.0 qni lartg _1j1   1j_1
0 -:!.0 qni lartg _1j1   1j1
0 -:!.0 qni lartg  1j_1 _1j_1
0 -:!.0 qni lartg  1j_1 _1j1
0 -:!.0 qni lartg  1j_1  1j_1
0 -:!.0 qni lartg  1j_1  1j1
0 -:!.0 qni lartg  1j1  _1j_1
0 -:!.0 qni lartg  1j1  _1j1
0 -:!.0 qni lartg  1j1   1j_1
0 -:!.0 qni lartg  1j1   1j1

NB. rot
NB. - input contains NaN
1 1 -: isnan  0     0j_. rot  0     0j_.
1 1 -: isnan  0     0j_. rot  0    _.
1 1 -: isnan  0     0j_. rot  0j_.  0
1 1 -: isnan  0     0j_. rot _.     0
1 1 -: isnan  0    _.    rot  0     0j_.
1 1 -: isnan  0    _.    rot  0    _.
1 1 -: isnan  0    _.    rot  0j_.  0
1 1 -: isnan  0    _.    rot _.     0
1 1 -: isnan  0j_.  0    rot  0     0j_.
1 1 -: isnan  0j_.  0    rot  0    _.
1 1 -: isnan  0j_.  0    rot  0j_.  0
1 1 -: isnan  0j_.  0    rot _.     0
1 1 -: isnan _.     0    rot  0     0j_.
1 1 -: isnan _.     0    rot  0    _.
1 1 -: isnan _.     0    rot  0j_.  0
1 1 -: isnan _.     0    rot _.     0
1 1 -: isnan  0    _.j_. rot  0    _.j_.
1 1 -: isnan  0    _.j_. rot  0j_.  0j_.
1 1 -: isnan  0    _.j_. rot  0j_. _.
1 1 -: isnan  0    _.j_. rot _.     0j_.
1 1 -: isnan  0    _.j_. rot _.    _.
1 1 -: isnan  0    _.j_. rot _.j_.  0
1 1 -: isnan  0j_.  0j_. rot  0    _.j_.
1 1 -: isnan  0j_.  0j_. rot  0j_.  0j_.
1 1 -: isnan  0j_.  0j_. rot  0j_. _.
1 1 -: isnan  0j_.  0j_. rot _.     0j_.
1 1 -: isnan  0j_.  0j_. rot _.    _.
1 1 -: isnan  0j_.  0j_. rot _.j_.  0
1 1 -: isnan  0j_. _.    rot  0    _.j_.
1 1 -: isnan  0j_. _.    rot  0j_.  0j_.
1 1 -: isnan  0j_. _.    rot  0j_. _.
1 1 -: isnan  0j_. _.    rot _.     0j_.
1 1 -: isnan  0j_. _.    rot _.    _.
1 1 -: isnan  0j_. _.    rot _.j_.  0
1 1 -: isnan _.     0j_. rot  0    _.j_.
1 1 -: isnan _.     0j_. rot  0j_.  0j_.
1 1 -: isnan _.     0j_. rot  0j_. _.
1 1 -: isnan _.     0j_. rot _.     0j_.
1 1 -: isnan _.     0j_. rot _.    _.
1 1 -: isnan _.     0j_. rot _.j_.  0
1 1 -: isnan _.    _.    rot  0    _.j_.
1 1 -: isnan _.    _.    rot  0j_.  0j_.
1 1 -: isnan _.    _.    rot  0j_. _.
1 1 -: isnan _.    _.    rot _.     0j_.
1 1 -: isnan _.    _.    rot _.    _.
1 1 -: isnan _.    _.    rot _.j_.  0
1 1 -: isnan _.j_.  0    rot  0    _.j_.
1 1 -: isnan _.j_.  0    rot  0j_.  0j_.
1 1 -: isnan _.j_.  0    rot  0j_. _.
1 1 -: isnan _.j_.  0    rot _.     0j_.
1 1 -: isnan _.j_.  0    rot _.    _.
1 1 -: isnan _.j_.  0    rot _.j_.  0
1 1 -: isnan _.j_. _.    rot _.j_. _.
1 1 -: isnan _.j_. _.    rot _.j_.  0j_.
1 1 -: isnan _.j_. _.    rot _.    _.j_.
1 1 -: isnan _.j_. _.    rot  0j_. _.j_.
1 1 -: isnan _.j_.  0j_. rot _.j_. _.
1 1 -: isnan _.j_.  0j_. rot _.j_.  0j_.
1 1 -: isnan _.j_.  0j_. rot _.    _.j_.
1 1 -: isnan _.j_.  0j_. rot  0j_. _.j_.
1 1 -: isnan _.    _.j_. rot _.j_. _.
1 1 -: isnan _.    _.j_. rot _.j_.  0j_.
1 1 -: isnan _.    _.j_. rot _.    _.j_.
1 1 -: isnan _.    _.j_. rot  0j_. _.j_.
1 1 -: isnan  0j_. _.j_. rot _.j_. _.
1 1 -: isnan  0j_. _.j_. rot _.j_.  0j_.
1 1 -: isnan  0j_. _.j_. rot _.    _.j_.
1 1 -: isnan  0j_. _.j_. rot  0j_. _.j_.
1 1 -: isnan _.j_. _.j_. rot _.j_. _.j_.
1 1 -: isnan         cs  rot  0     0j_.
1 1 -: isnan         cs  rot  0    _.
1 1 -: isnan         cs  rot  0j_.  0
1 1 -: isnan         cs  rot _.     0
1 1 -: isnan         cs  rot  0    _.j_.
1 1 -: isnan         cs  rot  0j_.  0j_.
1 1 -: isnan         cs  rot  0j_. _.
1 1 -: isnan         cs  rot _.     0j_.
1 1 -: isnan         cs  rot _.    _.
1 1 -: isnan         cs  rot _.j_.  0
1 1 -: isnan         cs  rot _.j_. _.
1 1 -: isnan         cs  rot _.j_.  0j_.
1 1 -: isnan         cs  rot _.    _.j_.
1 1 -: isnan         cs  rot  0j_. _.j_.
1 1 -: isnan         cs  rot _.j_. _.j_.
1 1 -: isnan         cs  rot _. qn1 q3
1 1 -: isnan         cs  rot _. qni q3
1 1 -: isnan         cs  rot _. qnj q3
1 1 -: isnan         cs  rot _. qnk q3
1 1 -: isnan (_. qn1 cs) rot        q3
1 1 -: isnan (_. qni cs) rot        q3
1 1 -: isnan (_. qnj cs) rot        q3
1 1 -: isnan (_. qnk cs) rot        q3
1 1 -: isnan  0     0j_. rot        q3
1 1 -: isnan  0    _.    rot        q3
1 1 -: isnan  0j_.  0    rot        q3
1 1 -: isnan _.     0    rot        q3
1 1 -: isnan  0    _.j_. rot        q3
1 1 -: isnan  0j_.  0j_. rot        q3
1 1 -: isnan  0j_. _.    rot        q3
1 1 -: isnan _.     0j_. rot        q3
1 1 -: isnan _.    _.    rot        q3
1 1 -: isnan _.j_.  0    rot        q3
1 1 -: isnan _.j_. _.    rot        q3
1 1 -: isnan _.j_.  0j_. rot        q3
1 1 -: isnan _.    _.j_. rot        q3
1 1 -: isnan  0j_. _.j_. rot        q3
1 1 -: isnan _.j_. _.j_. rot        q3
(,.~ 0 1) -: isnan  cs               rot q3 ,: _. qn1 q4
(,.~ 0 1) -: isnan  cs               rot q3 ,: _. qni q4
(,.~ 0 1) -: isnan  cs               rot q3 ,: _. qnj q4
(,.~ 0 1) -: isnan  cs               rot q3 ,: _. qnk q4
(,.~ 0 1) -: isnan (cs ,: _. qn1 cs) rot q3
(,.~ 0 1) -: isnan (cs ,: _. qni cs) rot q3
(,.~ 0 1) -: isnan (cs ,: _. qnj cs) rot q3
(,.~ 0 1) -: isnan (cs ,: _. qnk cs) rot q3
NB. - input without edge cases
0 0 -: cs rot 0 0
           q3  -:  1     0    rot q3
(-         q3) -: _1     0    rot q3
( +        q3) -:  0j_1  0    rot q3
(- &.(1&{) q3) -:  0     1    rot q3
(- &.(0&{) q3) -:  0    _1    rot q3
(-+        q3) -:  0     0j1  rot q3
( +        q3) -:  0     0j_1 rot q3
((,: -)  q3) -: (1 0 ,: _1 0) rot q3
(q3 ,:   q4) -:  1 0          rot q3 ,: q4
(q3 ,: - q4) -: (1 0 ,: _1 0) rot q3 ,: q4
