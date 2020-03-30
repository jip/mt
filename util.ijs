NB. Utilities
NB.
NB. negneg    Conditional negate
NB. negpos    Conditional negate
NB. copysign  Copy sign
NB. sorim     Sum of real and imaginary parts' modules
NB. soris     Sum of real and imaginary parts' squares
NB. fmtlog    Format log string
NB.
NB. Version: 0.10.2 2017-10-19
NB.
NB. Copyright 2010-2017 Igor Zhuravlov
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
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Miscellaneous

negneg=: -@]^:(0>[)                                             NB. if x<0 then -y else y endif
negpos=: -@]^:(0<:[)                                            NB. if x≥0 then -y else y endif

copysign=: -@]^:((=-)&*)                                        NB. if x<0 then -|y| else |y| endif

sorim=: +/"1@:| @:+.                                            NB. sum of real and imaginary parts' modules, |Re(y)| + |Im(y)|
soris=: +/"1@:*:@:+.                                            NB. sum of real and imaginary parts' squares, Re(y)^2 + Im(y)^2

fmtlog=: ;@:(40 17 17 17 17 _16&(({.{.@('d<n/a>'&(8!:2)))&.>))  NB. log string format
