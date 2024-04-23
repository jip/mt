NB. Apply a plane rotation
NB.
NB. drot  Apply a plane rotation to elements of vectors
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
NB. Configuration

coclass 'mtbla'

NB. =========================================================
NB. Includes

require 'math/mt/external/blas/blas'

NB. ---------------------------------------------------------
NB. drot
NB.
NB. Description:
NB.   Applies a plane rotation to elements of vectors
NB.
NB. Syntax:
NB.   'ox oy'=. drot ix ; iy ; c ; s
NB. where
NB.   ix - n-vector, 1st components of vectors to be rotated
NB.   iy - n-vector, 2nd components of vectors to be rotated
NB.   c  - scalar, the cosine of the plane rotation
NB.   s  - scalar, the sine of the plane rotation
NB.   ox - n-vector, 1st components of vectors rotated
NB.   oy - n-vector, 2nd components of vectors rotated
NB.
NB. Notes:
NB. - verb below is loaded into the current locale

drot=: 3 : 0
  'x y c s'=. y
  assert x =/&# y
  assert c , &# s
  2 4 { drot_cd_mtbla_ (, # x) ; (, x) ; (, 1) ; (, y) ; (, 1) ; (, c) ; (, s)
)
