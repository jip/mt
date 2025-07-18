NB. Verify iso actors
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

NB. 1) tests for liofmax are inspired by
NB.    https://doi.org/10.48550/arXiv.2207.09281, see
NB.    https://github.com/Reference-LAPACK/lapack/pull/1116
NB. 2) tests for liolmax are derived from liofmax's ones by
NB.    row-wise reversing
n=. 6  NB. make testset for N=6
iso=. (0 , 1 , <.@-: , <:) n
NB. IDAMAX tests
A=. (* _1&^) i. n
A1a=. iso _."_`[}"0 _ A
A1b=. (_2 ]\ 0 1 0 2 0 3 1 2 1 3 2 3 { iso) _. _."_`[}"1 _ A
A1c=. (_3 ]\ 0 1 2 0 1 3 0 2 3 1 2 3 { iso) _. _. _."_`[}"1 _ A
A1d=. _."0 A
A2ia=. 1 0 0 0 _:`[}"0 1 A1a
A2iia=. 1 0 0 0 __:`[}"0 1 A1a
A2iiia=. _1 _1 _1 _2 _:`[}"0 1 A1a
A2_iiia=. _1 _1 _1 _2 __:`[}"0 1 A1a
A2iva1=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (_ _"_)`[}"1 1 A1a
A2iva2=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (_ __"_)`[}"1 1 A1a
A2iva3=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (__ _"_)`[}"1 1 A1a
A2iva4=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (__ __"_)`[}"1 1 A1a
A2va=. (-. isnan A1a) (4 : 'x} y')"_1 A1a ,:"1 (]`-"0) _"0 A
A2ib=. 2 1 1 0 0 0 _:`[}"0 1 A1b
A2iib=. 2 1 1 0 0 0 __:`[}"0 1 A1b
A2iiib=. _1 _1 _2 _1 _2 _2 _:`[}"0 1 A1b
A2_iiib=. _1 _1 _2 _1 _2 _2 __:`[}"0 1 A1b
A2ivb1=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (_ _"_)`[}"1 1 A1b
A2ivb2=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (_ __"_)`[}"1 1 A1b
A2ivb3=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (__ _"_)`[}"1 1 A1b
A2ivb4=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (__ __"_)`[}"1 1 A1b
A2vb=. (-. isnan A1b) (4 : 'x} y')"_1 A1b ,:"1 (]`-"0) _"0 A
A2ic=. 2 2 1 0 _:`[}"0 1 A1c
A2iic=. 2 2 1 0 __:`[}"0 1 A1c
A2iiic=. _1 _2 _2 _2 _:`[}"0 1 A1c
A2_iiic=. _1 _2 _2 _2 __:`[}"0 1 A1c
A2ivc1=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (_ _"_)`[}"1 1 A1c
A2ivc2=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (_ __"_)`[}"1 1 A1c
A2ivc3=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (__ _"_)`[}"1 1 A1c
A2ivc4=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (__ __"_)`[}"1 1 A1c
A2vc=. (-. isnan A1c) (4 : 'x} y')"_1 A1c ,:"1 (]`-"0) _"0 A
A3a=. (n {."1 I.^:_1:"0 iso) (4 : 'x} y')"1 2 A ,:"1 (]`-"0) _"0 A
A3b=. (n {."1 (I.^:_1:"1) _2 ]\ 0 1 0 2 0 3 1 2 1 3 2 3 { iso) (4 : 'x} y')"1 2 A ,:"1 (]`-"0) _"0 A
A3c=. (n {."1 (I.^:_1:"1)_3 ]\ 0 1 2 0 1 3 0 2 3 1 2 3 { iso) (4 : 'x} y')"1 2 A ,:"1 (]`-"0) _"0 A
NB. IZAMAX tests
in=. i. n
A4A=. (2 | in)} ((FP_OVFL * j.~@(+&2 % +&3)) ,: (j.~ -)) in
A4B=. (0 1 ; 2 3 ; 4 5) C. A4A
A4C=. (2 | in)} ((FP_OVFL * j.~@(8&- % 9&-)) ,: (j.~ -)) in
A4D=. (0 1 ; 2 3 ; 4 5) C. A4C
A5A1a=. iso _."_`[}"0 _ A4A
A5A1b=. (_2 ]\ 0 1 0 2 0 3 1 2 1 3 2 3 { iso) _. _."_`[}"1 _ A4A
A5A1c=. (_3 ]\ 0 1 2 0 1 3 0 2 3 1 2 3 { iso) _. _. _."_`[}"1 _ A4A
A5A2ia=. 1 0 0 0 _:`[}"0 1 A5A1a
A5A2iia=. 1 0 0 0 __:`[}"0 1 A5A1a
A5A2iiia=. _1 _1 _1 _2 _:`[}"0 1 A5A1a
A5A2_iiia=. _1 _1 _1 _2 __:`[}"0 1 A5A1a
A5A2iva1=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (_ _"_)`[}"1 1 A5A1a
A5A2iva2=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (_ __"_)`[}"1 1 A5A1a
A5A2iva3=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (__ _"_)`[}"1 1 A5A1a
A5A2iva4=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (__ __"_)`[}"1 1 A5A1a
A5A2va=. (-. isnan A5A1a) (4 : 'x} y')"_1 A5A1a ,:"1 (]`-"0) _"0 A4A
A5A2ib=. 2 1 1 0 0 0 _:`[}"0 1 A5A1b
A5A2iib=. 2 1 1 0 0 0 __:`[}"0 1 A5A1b
A5A2iiib=. _1 _1 _2 _1 _2 _2 _:`[}"0 1 A5A1b
A5A2_iiib=. _1 _1 _2 _1 _2 _2 __:`[}"0 1 A5A1b
A5A2ivb1=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (_ _"_)`[}"1 1 A5A1b
A5A2ivb2=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (_ __"_)`[}"1 1 A5A1b
A5A2ivb3=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (__ _"_)`[}"1 1 A5A1b
A5A2ivb4=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (__ __"_)`[}"1 1 A5A1b
A5A2vb=. (-. isnan A5A1b) (4 : 'x} y')"_1 A5A1b ,:"1 (]`-"0) _"0 A4A
A5A2ic=. 2 2 1 0 _:`[}"0 1 A5A1c
A5A2iic=. 2 2 1 0 __:`[}"0 1 A5A1c
A5A2iiic=. _1 _2 _2 _2 _:`[}"0 1 A5A1c
A5A2_iiic=. _1 _2 _2 _2 __:`[}"0 1 A5A1c
A5A2ivc1=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (_ _"_)`[}"1 1 A5A1c
A5A2ivc2=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (_ __"_)`[}"1 1 A5A1c
A5A2ivc3=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (__ _"_)`[}"1 1 A5A1c
A5A2ivc4=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (__ __"_)`[}"1 1 A5A1c
A5A2vc=. (-. isnan A5A1c) (4 : 'x} y')"_1 A5A1c ,:"1 (]`-"0) _"0 A4A
A5A3a=. (n {."1 I.^:_1:"0 iso) (4 : 'x} y')"1 2 A4A ,:"1 (]`-"0) _"0 A4A
A5A3b=. (n {."1 (I.^:_1:"1) _2 ]\ 0 1 0 2 0 3 1 2 1 3 2 3 { iso) (4 : 'x} y')"1 2 A4A ,:"1 (]`-"0) _"0 A4A
A5A3c=. (n {."1 (I.^:_1:"1)_3 ]\ 0 1 2 0 1 3 0 2 3 1 2 3 { iso) (4 : 'x} y')"1 2 A4A ,:"1 (]`-"0) _"0 A4A
A5B1a=. iso _."_`[}"0 _ A4B
A5B1b=. (_2 ]\ 0 1 0 2 0 3 1 2 1 3 2 3 { iso) _. _."_`[}"1 _ A4B
A5B1c=. (_3 ]\ 0 1 2 0 1 3 0 2 3 1 2 3 { iso) _. _. _."_`[}"1 _ A4B
A5B2ia=. 1 0 0 0 _:`[}"0 1 A5B1a
A5B2iia=. 1 0 0 0 __:`[}"0 1 A5B1a
A5B2iiia=. _1 _1 _1 _2 _:`[}"0 1 A5B1a
A5B2_iiia=. _1 _1 _1 _2 __:`[}"0 1 A5B1a
A5B2iva1=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (_ _"_)`[}"1 1 A5B1a
A5B2iva2=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (_ __"_)`[}"1 1 A5B1a
A5B2iva3=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (__ _"_)`[}"1 1 A5B1a
A5B2iva4=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (__ __"_)`[}"1 1 A5B1a
A5B2va=. (-. isnan A5B1a) (4 : 'x} y')"_1 A5B1a ,:"1 (]`-"0) _"0 A4B
A5B2ib=. 2 1 1 0 0 0 _:`[}"0 1 A5B1b
A5B2iib=. 2 1 1 0 0 0 __:`[}"0 1 A5B1b
A5B2iiib=. _1 _1 _2 _1 _2 _2 _:`[}"0 1 A5B1b
A5B2_iiib=. _1 _1 _2 _1 _2 _2 __:`[}"0 1 A5B1b
A5B2ivb1=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (_ _"_)`[}"1 1 A5B1b
A5B2ivb2=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (_ __"_)`[}"1 1 A5B1b
A5B2ivb3=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (__ _"_)`[}"1 1 A5B1b
A5B2ivb4=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (__ __"_)`[}"1 1 A5B1b
A5B2vb=. (-. isnan A5B1b) (4 : 'x} y')"_1 A5B1b ,:"1 (]`-"0) _"0 A4B
A5B2ic=. 2 2 1 0 _:`[}"0 1 A5B1c
A5B2iic=. 2 2 1 0 __:`[}"0 1 A5B1c
A5B2iiic=. _1 _2 _2 _2 _:`[}"0 1 A5B1c
A5B2_iiic=. _1 _2 _2 _2 __:`[}"0 1 A5B1c
A5B2ivc1=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (_ _"_)`[}"1 1 A5B1c
A5B2ivc2=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (_ __"_)`[}"1 1 A5B1c
A5B2ivc3=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (__ _"_)`[}"1 1 A5B1c
A5B2ivc4=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (__ __"_)`[}"1 1 A5B1c
A5B2vc=. (-. isnan A5B1c) (4 : 'x} y')"_1 A5B1c ,:"1 (]`-"0) _"0 A4B
A5B3a=. (n {."1 I.^:_1:"0 iso) (4 : 'x} y')"1 2 A4B ,:"1 (]`-"0) _"0 A4B
A5B3b=. (n {."1 (I.^:_1:"1) _2 ]\ 0 1 0 2 0 3 1 2 1 3 2 3 { iso) (4 : 'x} y')"1 2 A4B ,:"1 (]`-"0) _"0 A4B
A5B3c=. (n {."1 (I.^:_1:"1)_3 ]\ 0 1 2 0 1 3 0 2 3 1 2 3 { iso) (4 : 'x} y')"1 2 A4B ,:"1 (]`-"0) _"0 A4B
A5C1a=. iso _."_`[}"0 _ A4C
A5C1b=. (_2 ]\ 0 1 0 2 0 3 1 2 1 3 2 3 { iso) _. _."_`[}"1 _ A4C
A5C1c=. (_3 ]\ 0 1 2 0 1 3 0 2 3 1 2 3 { iso) _. _. _."_`[}"1 _ A4C
A5C2ia=. 1 0 0 0 _:`[}"0 1 A5C1a
A5C2iia=. 1 0 0 0 __:`[}"0 1 A5C1a
A5C2iiia=. _1 _1 _1 _2 _:`[}"0 1 A5C1a
A5C2_iiia=. _1 _1 _1 _2 __:`[}"0 1 A5C1a
A5C2iva1=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (_ _"_)`[}"1 1 A5C1a
A5C2iva2=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (_ __"_)`[}"1 1 A5C1a
A5C2iva3=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (__ _"_)`[}"1 1 A5C1a
A5C2iva4=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (__ __"_)`[}"1 1 A5C1a
A5C2va=. (-. isnan A5C1a) (4 : 'x} y')"_1 A5C1a ,:"1 (]`-"0) _"0 A4C
A5C2ib=. 2 1 1 0 0 0 _:`[}"0 1 A5C1b
A5C2iib=. 2 1 1 0 0 0 __:`[}"0 1 A5C1b
A5C2iiib=. _1 _1 _2 _1 _2 _2 _:`[}"0 1 A5C1b
A5C2_iiib=. _1 _1 _2 _1 _2 _2 __:`[}"0 1 A5C1b
A5C2ivb1=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (_ _"_)`[}"1 1 A5C1b
A5C2ivb2=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (_ __"_)`[}"1 1 A5C1b
A5C2ivb3=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (__ _"_)`[}"1 1 A5C1b
A5C2ivb4=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (__ __"_)`[}"1 1 A5C1b
A5C2vb=. (-. isnan A5C1b) (4 : 'x} y')"_1 A5C1b ,:"1 (]`-"0) _"0 A4C
A5C2ic=. 2 2 1 0 _:`[}"0 1 A5C1c
A5C2iic=. 2 2 1 0 __:`[}"0 1 A5C1c
A5C2iiic=. _1 _2 _2 _2 _:`[}"0 1 A5C1c
A5C2_iiic=. _1 _2 _2 _2 __:`[}"0 1 A5C1c
A5C2ivc1=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (_ _"_)`[}"1 1 A5C1c
A5C2ivc2=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (_ __"_)`[}"1 1 A5C1c
A5C2ivc3=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (__ _"_)`[}"1 1 A5C1c
A5C2ivc4=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (__ __"_)`[}"1 1 A5C1c
A5C2vc=. (-. isnan A5C1c) (4 : 'x} y')"_1 A5C1c ,:"1 (]`-"0) _"0 A4C
A5C3a=. (n {."1 I.^:_1:"0 iso) (4 : 'x} y')"1 2 A4C ,:"1 (]`-"0) _"0 A4C
A5C3b=. (n {."1 (I.^:_1:"1) _2 ]\ 0 1 0 2 0 3 1 2 1 3 2 3 { iso) (4 : 'x} y')"1 2 A4C ,:"1 (]`-"0) _"0 A4C
A5C3c=. (n {."1 (I.^:_1:"1)_3 ]\ 0 1 2 0 1 3 0 2 3 1 2 3 { iso) (4 : 'x} y')"1 2 A4C ,:"1 (]`-"0) _"0 A4C
A5D1a=. iso _."_`[}"0 _ A4D
A5D1b=. (_2 ]\ 0 1 0 2 0 3 1 2 1 3 2 3 { iso) _. _."_`[}"1 _ A4D
A5D1c=. (_3 ]\ 0 1 2 0 1 3 0 2 3 1 2 3 { iso) _. _. _."_`[}"1 _ A4D
A5D2ia=. 1 0 0 0 _:`[}"0 1 A5D1a
A5D2iia=. 1 0 0 0 __:`[}"0 1 A5D1a
A5D2iiia=. _1 _1 _1 _2 _:`[}"0 1 A5D1a
A5D2_iiia=. _1 _1 _1 _2 __:`[}"0 1 A5D1a
A5D2iva1=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (_ _"_)`[}"1 1 A5D1a
A5D2iva2=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (_ __"_)`[}"1 1 A5D1a
A5D2iva3=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (__ _"_)`[}"1 1 A5D1a
A5D2iva4=. (_2 ]\ 1 _1 0 _1 0 _1 0 _2) (__ __"_)`[}"1 1 A5D1a
A5D2va=. (-. isnan A5D1a) (4 : 'x} y')"_1 A5D1a ,:"1 (]`-"0) _"0 A4D
A5D2ib=. 2 1 1 0 0 0 _:`[}"0 1 A5D1b
A5D2iib=. 2 1 1 0 0 0 __:`[}"0 1 A5D1b
A5D2iiib=. _1 _1 _2 _1 _2 _2 _:`[}"0 1 A5D1b
A5D2_iiib=. _1 _1 _2 _1 _2 _2 __:`[}"0 1 A5D1b
A5D2ivb1=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (_ _"_)`[}"1 1 A5D1b
A5D2ivb2=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (_ __"_)`[}"1 1 A5D1b
A5D2ivb3=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (__ _"_)`[}"1 1 A5D1b
A5D2ivb4=. (_2 ]\ 2 _1 1 _1 1 _2 0 _1 0 _2 0 _2) (__ __"_)`[}"1 1 A5D1b
A5D2vb=. (-. isnan A5D1b) (4 : 'x} y')"_1 A5D1b ,:"1 (]`-"0) _"0 A4D
A5D2ic=. 2 2 1 0 _:`[}"0 1 A5D1c
A5D2iic=. 2 2 1 0 __:`[}"0 1 A5D1c
A5D2iiic=. _1 _2 _2 _2 _:`[}"0 1 A5D1c
A5D2_iiic=. _1 _2 _2 _2 __:`[}"0 1 A5D1c
A5D2ivc1=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (_ _"_)`[}"1 1 A5D1c
A5D2ivc2=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (_ __"_)`[}"1 1 A5D1c
A5D2ivc3=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (__ _"_)`[}"1 1 A5D1c
A5D2ivc4=. (_2 ]\ 2 _1 2 _2 1 _2 0 _2) (__ __"_)`[}"1 1 A5D1c
A5D2vc=. (-. isnan A5D1c) (4 : 'x} y')"_1 A5D1c ,:"1 (]`-"0) _"0 A4D
A5D3a=. (n {."1 I.^:_1:"0 iso) (4 : 'x} y')"1 2 A4D ,:"1 (]`-"0) _"0 A4D
A5D3b=. (n {."1 (I.^:_1:"1) _2 ]\ 0 1 0 2 0 3 1 2 1 3 2 3 { iso) (4 : 'x} y')"1 2 A4D ,:"1 (]`-"0) _"0 A4D
A5D3c=. (n {."1 (I.^:_1:"1)_3 ]\ 0 1 2 0 1 3 0 2 3 1 2 3 { iso) (4 : 'x} y')"1 2 A4D ,:"1 (]`-"0) _"0 A4D

NB. =========================================================
NB. Verification suite

NB. liofmax
0 -: liofmax ''
0 -: liofmax _.
0 -: liofmax __
0 -: liofmax 0
0 -: liofmax _
0 -: liofmax _. _.
0 -: liofmax _. __
0 -: liofmax _.  0
0 -: liofmax _.  _
1 -: liofmax __ _.
0 -: liofmax __ __
0 -: liofmax __  0
0 -: liofmax __  _
1 -: liofmax  0 _.
1 -: liofmax  0 __
0 -: liofmax  0  0
1 -: liofmax  0  _
1 -: liofmax  _ _.
0 -: liofmax  _ __
0 -: liofmax  _  0
0 -: liofmax  _  _
0 -: liofmax i: 3
NB. IDAMAX tests
5 -: liofmax A
0 1 3 5 -: liofmax"1 A1a
0 0 0 1 1 3 -: liofmax"1 A1b
0 0 0 1 -: liofmax"1 A1c
0 -: liofmax A1d
0 1 3 5 -: liofmax"1 A2ia
0 1 3 5 -: liofmax"1 A2iia
0 1 3 5 -: liofmax"1 A2iiia
0 1 3 5 -: liofmax"1 A2_iiia
0 1 3 5 -: liofmax"1 A2iva1
0 1 3 5 -: liofmax"1 A2iva2
0 1 3 5 -: liofmax"1 A2iva3
0 1 3 5 -: liofmax"1 A2iva4
0 1 3 5 -: liofmax"1 A2va
0 0 0 1 1 3 -: liofmax"1 A2ib
0 0 0 1 1 3 -: liofmax"1 A2iib
0 0 0 1 1 3 -: liofmax"1 A2iiib
0 0 0 1 1 3 -: liofmax"1 A2_iiib
0 0 0 1 1 3 -: liofmax"1 A2ivb1
0 0 0 1 1 3 -: liofmax"1 A2ivb2
0 0 0 1 1 3 -: liofmax"1 A2ivb3
0 0 0 1 1 3 -: liofmax"1 A2ivb4
0 0 0 1 1 3 -: liofmax"1 A2vb
0 0 0 1 -: liofmax"1 A2ic
0 0 0 1 -: liofmax"1 A2iic
0 0 0 1 -: liofmax"1 A2iiic
0 0 0 1 -: liofmax"1 A2_iiic
0 0 0 1 -: liofmax"1 A2ivc1
0 0 0 1 -: liofmax"1 A2ivc2
0 0 0 1 -: liofmax"1 A2ivc3
0 0 0 1 -: liofmax"1 A2ivc4
0 0 0 1 -: liofmax"1 A2vc
0 1 3 5 -: liofmax"1 A3a
0 0 0 1 1 3 -: liofmax"1 A3b
0 0 0 1 -: liofmax"1 A3c
NB. IZAMAX tests
4 -: liofmax A4A
5 -: liofmax A4B
0 -: liofmax A4C
1 -: liofmax A4D
0 1 3 5 -: liofmax"1 A5A1a
0 0 0 1 1 3 -: liofmax"1 A5A1b
0 0 0 1 -: liofmax"1 A5A1c
0 1 3 5 -: liofmax"1 A5A2ia
0 1 3 5 -: liofmax"1 A5A2iia
0 1 3 5 -: liofmax"1 A5A2iiia
0 1 3 5 -: liofmax"1 A5A2_iiia
0 1 3 5 -: liofmax"1 A5A2iva1
0 1 3 5 -: liofmax"1 A5A2iva2
0 1 3 5 -: liofmax"1 A5A2iva3
0 1 3 5 -: liofmax"1 A5A2iva4
0 1 3 5 -: liofmax"1 A5A2va
0 0 0 1 1 3 -: liofmax"1 A5A2ib
0 0 0 1 1 3 -: liofmax"1 A5A2iib
0 0 0 1 1 3 -: liofmax"1 A5A2iiib
0 0 0 1 1 3 -: liofmax"1 A5A2_iiib
0 0 0 1 1 3 -: liofmax"1 A5A2ivb1
0 0 0 1 1 3 -: liofmax"1 A5A2ivb2
0 0 0 1 1 3 -: liofmax"1 A5A2ivb3
0 0 0 1 1 3 -: liofmax"1 A5A2ivb4
0 0 0 1 1 3 -: liofmax"1 A5A2vb
0 0 0 1 -: liofmax"1 A5A2ic
0 0 0 1 -: liofmax"1 A5A2iic
0 0 0 1 -: liofmax"1 A5A2iiic
0 0 0 1 -: liofmax"1 A5A2_iiic
0 0 0 1 -: liofmax"1 A5A2ivc1
0 0 0 1 -: liofmax"1 A5A2ivc2
0 0 0 1 -: liofmax"1 A5A2ivc3
0 0 0 1 -: liofmax"1 A5A2ivc4
0 0 0 1 -: liofmax"1 A5A2vc
0 1 3 5 -: liofmax"1 A5A3a
0 0 0 1 1 3 -: liofmax"1 A5A3b
0 0 0 1 -: liofmax"1 A5A3c
0 1 3 5 -: liofmax"1 A5B1a
0 0 0 1 1 3 -: liofmax"1 A5B1b
0 0 0 1 -: liofmax"1 A5B1c
0 1 3 5 -: liofmax"1 A5B2ia
0 1 3 5 -: liofmax"1 A5B2iia
0 1 3 5 -: liofmax"1 A5B2iiia
0 1 3 5 -: liofmax"1 A5B2_iiia
0 1 3 5 -: liofmax"1 A5B2iva1
0 1 3 5 -: liofmax"1 A5B2iva2
0 1 3 5 -: liofmax"1 A5B2iva3
0 1 3 5 -: liofmax"1 A5B2iva4
0 1 3 5 -: liofmax"1 A5B2va
0 0 0 1 1 3 -: liofmax"1 A5B2ib
0 0 0 1 1 3 -: liofmax"1 A5B2iib
0 0 0 1 1 3 -: liofmax"1 A5B2iiib
0 0 0 1 1 3 -: liofmax"1 A5B2_iiib
0 0 0 1 1 3 -: liofmax"1 A5B2ivb1
0 0 0 1 1 3 -: liofmax"1 A5B2ivb2
0 0 0 1 1 3 -: liofmax"1 A5B2ivb3
0 0 0 1 1 3 -: liofmax"1 A5B2ivb4
0 0 0 1 1 3 -: liofmax"1 A5B2vb
0 0 0 1 -: liofmax"1 A5B2ic
0 0 0 1 -: liofmax"1 A5B2iic
0 0 0 1 -: liofmax"1 A5B2iiic
0 0 0 1 -: liofmax"1 A5B2_iiic
0 0 0 1 -: liofmax"1 A5B2ivc1
0 0 0 1 -: liofmax"1 A5B2ivc2
0 0 0 1 -: liofmax"1 A5B2ivc3
0 0 0 1 -: liofmax"1 A5B2ivc4
0 0 0 1 -: liofmax"1 A5B2vc
0 1 3 5 -: liofmax"1 A5B3a
0 0 0 1 1 3 -: liofmax"1 A5B3b
0 0 0 1 -: liofmax"1 A5B3c
0 1 3 5 -: liofmax"1 A5C1a
0 0 0 1 1 3 -: liofmax"1 A5C1b
0 0 0 1 -: liofmax"1 A5C1c
0 1 3 5 -: liofmax"1 A5C2ia
0 1 3 5 -: liofmax"1 A5C2iia
0 1 3 5 -: liofmax"1 A5C2iiia
0 1 3 5 -: liofmax"1 A5C2_iiia
0 1 3 5 -: liofmax"1 A5C2iva1
0 1 3 5 -: liofmax"1 A5C2iva2
0 1 3 5 -: liofmax"1 A5C2iva3
0 1 3 5 -: liofmax"1 A5C2iva4
0 1 3 5 -: liofmax"1 A5C2va
0 0 0 1 1 3 -: liofmax"1 A5C2ib
0 0 0 1 1 3 -: liofmax"1 A5C2iib
0 0 0 1 1 3 -: liofmax"1 A5C2iiib
0 0 0 1 1 3 -: liofmax"1 A5C2_iiib
0 0 0 1 1 3 -: liofmax"1 A5C2ivb1
0 0 0 1 1 3 -: liofmax"1 A5C2ivb2
0 0 0 1 1 3 -: liofmax"1 A5C2ivb3
0 0 0 1 1 3 -: liofmax"1 A5C2ivb4
0 0 0 1 1 3 -: liofmax"1 A5C2vb
0 0 0 1 -: liofmax"1 A5C2ic
0 0 0 1 -: liofmax"1 A5C2iic
0 0 0 1 -: liofmax"1 A5C2iiic
0 0 0 1 -: liofmax"1 A5C2_iiic
0 0 0 1 -: liofmax"1 A5C2ivc1
0 0 0 1 -: liofmax"1 A5C2ivc2
0 0 0 1 -: liofmax"1 A5C2ivc3
0 0 0 1 -: liofmax"1 A5C2ivc4
0 0 0 1 -: liofmax"1 A5C2vc
0 1 3 5 -: liofmax"1 A5C3a
0 0 0 1 1 3 -: liofmax"1 A5C3b
0 0 0 1 -: liofmax"1 A5C3c
0 1 3 5 -: liofmax"1 A5D1a
0 0 0 1 1 3 -: liofmax"1 A5D1b
0 0 0 1 -: liofmax"1 A5D1c
0 1 3 5 -: liofmax"1 A5D2ia
0 1 3 5 -: liofmax"1 A5D2iia
0 1 3 5 -: liofmax"1 A5D2iiia
0 1 3 5 -: liofmax"1 A5D2_iiia
0 1 3 5 -: liofmax"1 A5D2iva1
0 1 3 5 -: liofmax"1 A5D2iva2
0 1 3 5 -: liofmax"1 A5D2iva3
0 1 3 5 -: liofmax"1 A5D2iva4
0 1 3 5 -: liofmax"1 A5D2va
0 0 0 1 1 3 -: liofmax"1 A5D2ib
0 0 0 1 1 3 -: liofmax"1 A5D2iib
0 0 0 1 1 3 -: liofmax"1 A5D2iiib
0 0 0 1 1 3 -: liofmax"1 A5D2_iiib
0 0 0 1 1 3 -: liofmax"1 A5D2ivb1
0 0 0 1 1 3 -: liofmax"1 A5D2ivb2
0 0 0 1 1 3 -: liofmax"1 A5D2ivb3
0 0 0 1 1 3 -: liofmax"1 A5D2ivb4
0 0 0 1 1 3 -: liofmax"1 A5D2vb
0 0 0 1 -: liofmax"1 A5D2ic
0 0 0 1 -: liofmax"1 A5D2iic
0 0 0 1 -: liofmax"1 A5D2iiic
0 0 0 1 -: liofmax"1 A5D2_iiic
0 0 0 1 -: liofmax"1 A5D2ivc1
0 0 0 1 -: liofmax"1 A5D2ivc2
0 0 0 1 -: liofmax"1 A5D2ivc3
0 0 0 1 -: liofmax"1 A5D2ivc4
0 0 0 1 -: liofmax"1 A5D2vc
0 1 3 5 -: liofmax"1 A5D3a
0 0 0 1 1 3 -: liofmax"1 A5D3b
0 0 0 1 -: liofmax"1 A5D3c

NB. liolmax
0 -: liolmax ''
0 -: liolmax _.
0 -: liolmax __
0 -: liolmax 0
0 -: liolmax _
1 -: liolmax _. _.
0 -: liolmax _. __
0 -: liolmax _.  0
0 -: liolmax _.  _
1 -: liolmax __ _.
1 -: liolmax __ __
0 -: liolmax __  0
1 -: liolmax __  _
1 -: liolmax  0 _.
1 -: liolmax  0 __
1 -: liolmax  0  0
1 -: liolmax  0  _
1 -: liolmax  _ _.
1 -: liolmax  _ __
0 -: liolmax  _  0
1 -: liolmax  _  _
6 -: liolmax i: 3
NB. tests for real datatype
0 -: liolmax |. A
5 4 2 0 -: liolmax@|."1 A1a
5 5 5 4 4 2 -: liolmax@|."1 A1b
5 5 5 4 -: liolmax@|."1 A1c
5 -: liolmax |. A1d
5 4 2 0 -: liolmax@|."1 A2ia
5 4 2 0 -: liolmax@|."1 A2iia
5 4 2 0 -: liolmax@|."1 A2iiia
5 4 2 0 -: liolmax@|."1 A2_iiia
5 4 2 0 -: liolmax@|."1 A2iva1
5 4 2 0 -: liolmax@|."1 A2iva2
5 4 2 0 -: liolmax@|."1 A2iva3
5 4 2 0 -: liolmax@|."1 A2iva4
5 4 2 0 -: liolmax@|."1 A2va
5 5 5 4 4 2 -: liolmax@|."1 A2ib
5 5 5 4 4 2 -: liolmax@|."1 A2iib
5 5 5 4 4 2 -: liolmax@|."1 A2iiib
5 5 5 4 4 2 -: liolmax@|."1 A2_iiib
5 5 5 4 4 2 -: liolmax@|."1 A2ivb1
5 5 5 4 4 2 -: liolmax@|."1 A2ivb2
5 5 5 4 4 2 -: liolmax@|."1 A2ivb3
5 5 5 4 4 2 -: liolmax@|."1 A2ivb4
5 5 5 4 4 2 -: liolmax@|."1 A2vb
5 5 5 4 -: liolmax@|."1 A2ic
5 5 5 4 -: liolmax@|."1 A2iic
5 5 5 4 -: liolmax@|."1 A2iiic
5 5 5 4 -: liolmax@|."1 A2_iiic
5 5 5 4 -: liolmax@|."1 A2ivc1
5 5 5 4 -: liolmax@|."1 A2ivc2
5 5 5 4 -: liolmax@|."1 A2ivc3
5 5 5 4 -: liolmax@|."1 A2ivc4
5 5 5 4 -: liolmax@|."1 A2vc
5 4 2 0 -: liolmax@|."1 A3a
5 5 5 4 4 2 -: liolmax@|."1 A3b
5 5 5 4 -: liolmax@|."1 A3c
NB. tests for complex datatype
1 -: liolmax |. A4A
0 -: liolmax |. A4B
5 -: liolmax |. A4C
4 -: liolmax |. A4D
5 4 2 0 -: liolmax@|."1 A5A1a
5 5 5 4 4 2 -: liolmax@|."1 A5A1b
5 5 5 4 -: liolmax@|."1 A5A1c
5 4 2 0 -: liolmax@|."1 A5A2ia
5 4 2 0 -: liolmax@|."1 A5A2iia
5 4 2 0 -: liolmax@|."1 A5A2iiia
5 4 2 0 -: liolmax@|."1 A5A2_iiia
5 4 2 0 -: liolmax@|."1 A5A2iva1
5 4 2 0 -: liolmax@|."1 A5A2iva2
5 4 2 0 -: liolmax@|."1 A5A2iva3
5 4 2 0 -: liolmax@|."1 A5A2iva4
5 4 2 0 -: liolmax@|."1 A5A2va
5 5 5 4 4 2 -: liolmax@|."1 A5A2ib
5 5 5 4 4 2 -: liolmax@|."1 A5A2iib
5 5 5 4 4 2 -: liolmax@|."1 A5A2iiib
5 5 5 4 4 2 -: liolmax@|."1 A5A2_iiib
5 5 5 4 4 2 -: liolmax@|."1 A5A2ivb1
5 5 5 4 4 2 -: liolmax@|."1 A5A2ivb2
5 5 5 4 4 2 -: liolmax@|."1 A5A2ivb3
5 5 5 4 4 2 -: liolmax@|."1 A5A2ivb4
5 5 5 4 4 2 -: liolmax@|."1 A5A2vb
5 5 5 4 -: liolmax@|."1 A5A2ic
5 5 5 4 -: liolmax@|."1 A5A2iic
5 5 5 4 -: liolmax@|."1 A5A2iiic
5 5 5 4 -: liolmax@|."1 A5A2_iiic
5 5 5 4 -: liolmax@|."1 A5A2ivc1
5 5 5 4 -: liolmax@|."1 A5A2ivc2
5 5 5 4 -: liolmax@|."1 A5A2ivc3
5 5 5 4 -: liolmax@|."1 A5A2ivc4
5 5 5 4 -: liolmax@|."1 A5A2vc
5 4 2 0 -: liolmax@|."1 A5A3a
5 5 5 4 4 2 -: liolmax@|."1 A5A3b
5 5 5 4 -: liolmax@|."1 A5A3c
5 4 2 0 -: liolmax@|."1 A5B1a
5 5 5 4 4 2 -: liolmax@|."1 A5B1b
5 5 5 4 -: liolmax@|."1 A5B1c
5 4 2 0 -: liolmax@|."1 A5B2ia
5 4 2 0 -: liolmax@|."1 A5B2iia
5 4 2 0 -: liolmax@|."1 A5B2iiia
5 4 2 0 -: liolmax@|."1 A5B2_iiia
5 4 2 0 -: liolmax@|."1 A5B2iva1
5 4 2 0 -: liolmax@|."1 A5B2iva2
5 4 2 0 -: liolmax@|."1 A5B2iva3
5 4 2 0 -: liolmax@|."1 A5B2iva4
5 4 2 0 -: liolmax@|."1 A5B2va
5 5 5 4 4 2 -: liolmax@|."1 A5B2ib
5 5 5 4 4 2 -: liolmax@|."1 A5B2iib
5 5 5 4 4 2 -: liolmax@|."1 A5B2iiib
5 5 5 4 4 2 -: liolmax@|."1 A5B2_iiib
5 5 5 4 4 2 -: liolmax@|."1 A5B2ivb1
5 5 5 4 4 2 -: liolmax@|."1 A5B2ivb2
5 5 5 4 4 2 -: liolmax@|."1 A5B2ivb3
5 5 5 4 4 2 -: liolmax@|."1 A5B2ivb4
5 5 5 4 4 2 -: liolmax@|."1 A5B2vb
5 5 5 4 -: liolmax@|."1 A5B2ic
5 5 5 4 -: liolmax@|."1 A5B2iic
5 5 5 4 -: liolmax@|."1 A5B2iiic
5 5 5 4 -: liolmax@|."1 A5B2_iiic
5 5 5 4 -: liolmax@|."1 A5B2ivc1
5 5 5 4 -: liolmax@|."1 A5B2ivc2
5 5 5 4 -: liolmax@|."1 A5B2ivc3
5 5 5 4 -: liolmax@|."1 A5B2ivc4
5 5 5 4 -: liolmax@|."1 A5B2vc
5 4 2 0 -: liolmax@|."1 A5B3a
5 5 5 4 4 2 -: liolmax@|."1 A5B3b
5 5 5 4 -: liolmax@|."1 A5B3c
5 4 2 0 -: liolmax@|."1 A5C1a
5 5 5 4 4 2 -: liolmax@|."1 A5C1b
5 5 5 4 -: liolmax@|."1 A5C1c
5 4 2 0 -: liolmax@|."1 A5C2ia
5 4 2 0 -: liolmax@|."1 A5C2iia
5 4 2 0 -: liolmax@|."1 A5C2iiia
5 4 2 0 -: liolmax@|."1 A5C2_iiia
5 4 2 0 -: liolmax@|."1 A5C2iva1
5 4 2 0 -: liolmax@|."1 A5C2iva2
5 4 2 0 -: liolmax@|."1 A5C2iva3
5 4 2 0 -: liolmax@|."1 A5C2iva4
5 4 2 0 -: liolmax@|."1 A5C2va
5 5 5 4 4 2 -: liolmax@|."1 A5C2ib
5 5 5 4 4 2 -: liolmax@|."1 A5C2iib
5 5 5 4 4 2 -: liolmax@|."1 A5C2iiib
5 5 5 4 4 2 -: liolmax@|."1 A5C2_iiib
5 5 5 4 4 2 -: liolmax@|."1 A5C2ivb1
5 5 5 4 4 2 -: liolmax@|."1 A5C2ivb2
5 5 5 4 4 2 -: liolmax@|."1 A5C2ivb3
5 5 5 4 4 2 -: liolmax@|."1 A5C2ivb4
5 5 5 4 4 2 -: liolmax@|."1 A5C2vb
5 5 5 4 -: liolmax@|."1 A5C2ic
5 5 5 4 -: liolmax@|."1 A5C2iic
5 5 5 4 -: liolmax@|."1 A5C2iiic
5 5 5 4 -: liolmax@|."1 A5C2_iiic
5 5 5 4 -: liolmax@|."1 A5C2ivc1
5 5 5 4 -: liolmax@|."1 A5C2ivc2
5 5 5 4 -: liolmax@|."1 A5C2ivc3
5 5 5 4 -: liolmax@|."1 A5C2ivc4
5 5 5 4 -: liolmax@|."1 A5C2vc
5 4 2 0 -: liolmax@|."1 A5C3a
5 5 5 4 4 2 -: liolmax@|."1 A5C3b
5 5 5 4 -: liolmax@|."1 A5C3c
5 4 2 0 -: liolmax@|."1 A5D1a
5 5 5 4 4 2 -: liolmax@|."1 A5D1b
5 5 5 4 -: liolmax@|."1 A5D1c
5 4 2 0 -: liolmax@|."1 A5D2ia
5 4 2 0 -: liolmax@|."1 A5D2iia
5 4 2 0 -: liolmax@|."1 A5D2iiia
5 4 2 0 -: liolmax@|."1 A5D2_iiia
5 4 2 0 -: liolmax@|."1 A5D2iva1
5 4 2 0 -: liolmax@|."1 A5D2iva2
5 4 2 0 -: liolmax@|."1 A5D2iva3
5 4 2 0 -: liolmax@|."1 A5D2iva4
5 4 2 0 -: liolmax@|."1 A5D2va
5 5 5 4 4 2 -: liolmax@|."1 A5D2ib
5 5 5 4 4 2 -: liolmax@|."1 A5D2iib
5 5 5 4 4 2 -: liolmax@|."1 A5D2iiib
5 5 5 4 4 2 -: liolmax@|."1 A5D2_iiib
5 5 5 4 4 2 -: liolmax@|."1 A5D2ivb1
5 5 5 4 4 2 -: liolmax@|."1 A5D2ivb2
5 5 5 4 4 2 -: liolmax@|."1 A5D2ivb3
5 5 5 4 4 2 -: liolmax@|."1 A5D2ivb4
5 5 5 4 4 2 -: liolmax@|."1 A5D2vb
5 5 5 4 -: liolmax@|."1 A5D2ic
5 5 5 4 -: liolmax@|."1 A5D2iic
5 5 5 4 -: liolmax@|."1 A5D2iiic
5 5 5 4 -: liolmax@|."1 A5D2_iiic
5 5 5 4 -: liolmax@|."1 A5D2ivc1
5 5 5 4 -: liolmax@|."1 A5D2ivc2
5 5 5 4 -: liolmax@|."1 A5D2ivc3
5 5 5 4 -: liolmax@|."1 A5D2ivc4
5 5 5 4 -: liolmax@|."1 A5D2vc
5 4 2 0 -: liolmax@|."1 A5D3a
5 5 5 4 4 2 -: liolmax@|."1 A5D3b
5 5 5 4 -: liolmax@|."1 A5D3c

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
