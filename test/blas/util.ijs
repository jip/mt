NB. Utilities
NB.
NB. basicxxx  Utilities to either check or modify argument
NB.
NB. Version: 0.14.0 2023-03-21
NB.
NB. Copyright 2010-2023 Igor Zhuravlov
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
NB. Configuration

coclass 'mtbla'
coinsert 'mt'

NB. =========================================================
NB. Interface

issquare=: =/@$  NB. same as in the (math/lapack2) addon

NB. check
NB. - ranks
basiccr0=: 0 2 2         -: #@$S:0
basiccr1=: 2 1 0         -: #@$S:0
basiccr2=: 0 1 0 2       -: #@$S:0
basiccr3=: 0 2 0 2       -: #@$S:0
basiccr4=: 0 2 2 0 2     -: #@$S:0
basiccr5=: 0 1 0 1 0 2   -: #@$S:0
basiccr6=: 0 2 1 0 0 1 0 -: #@$S:0
NB. - shape
basiccs0=: issquare@(0&{::)
basiccs1=: issquare@(1&{::)
basiccs3=: issquare@(3&{::)
basiccs4=: issquare@(4&{::)
basiccs5=: issquare@(5&{::)
NB. - compare shapes
basiccmp=: -:/@($L:0)@(1 2&{)

NB. modify
NB. - conjugate under ISO specified
basiccj0=: 1      &(+&.> upd)
basiccj1=: 0 1 3  &(+&.> upd)
basiccj2=: 0 2 4 5&(+&.> upd)
NB. - swap elements
basicswp=: (< 1 2)&C.
