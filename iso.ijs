NB. ISO
NB.
NB. liofmax    lIO 1st element with maximum sum of real and
NB.            imagine parts' modules
NB. liolmax    lIO last element with maximum sum of real and
NB.            imagine parts' modules
NB. th2liso    Generate lISO from tail and head
NB. dhs2liso   Generate lISO from head, size and optional
NB.            delta
NB. riso2iso   Convert rISO to ISO
NB. riso2liso  Convert rISO to lISO
NB. lisoX      lISO vector laying between diagonal and
NB.            matrix edge
NB.
NB. Version: 0.11.0 2021-01-17
NB.
NB. Copyright 2010-2021 Igor Zhuravlov
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
NB. Concepts
NB.
NB. IO   - index of
NB. lIO  - linear IO, is an integer
NB. ISO  - indices of
NB. lISO - linear ISO, is a vector of integers
NB. rISO - rectangular ISO, for r-rank array is a 2×r-array
NB.        of integers ((head0,head1,...),:(size0,size1,...))
NB.
NB. Following are equivalents:
NB.   (3 5 _7 ,: 2 _3 4) ];.0 brick
NB.   (< 3 4 ; 7 6 5 ; _10 _9 _8 _7) { brick
NB.   (riso2iso 3 5 _7 ,: 2 _3 4) { brick
NB.   (iso2riso < 3 4 ; 7 6 5 ; _10 _9 _8 _7) ];.0 brick
NB.
NB. Following are equivalents:
NB.   (0 1 ; 1 2 ; 2 3) { i. 3 4
NB.   (4 liso2iso 1 6 11) { i. 3 4
NB.   (4 iso2liso 0 1 ; 1 2 ; 2 3) ({,) i. 3 4
NB.   1 6 11 ({,) i. 3 4

NB. =========================================================
NB. Local definitions

NB. convert rISO to opened (non-boxed) ISO
riso2oiso=: <@dhs2liso/"1@:|:

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Miscellaneous

NB. lIO 1st element e with max(|Re(e)|+|Im(e)|) from list y
NB. implements BLAS's IxAMAX
liofmax=: (i.>./)@sorim

NB. lIO last element e with max(|Re(e)|+|Im(e)|) from list y
liolmax=: (i:>./)@sorim

NB. ---------------------------------------------------------
NB. th2liso
NB.
NB. Description:
NB.   Generate lISO from tail and head
NB.
NB. Syntax:
NB.   liso=. t th2liso h
NB. where
NB.   h    - integer, head of liso
NB.   t    - integer, tail of liso
NB.   liso - (x-y)-vector of integers, lISO from head h to
NB.          tail (t-1) with delta=1:
NB.            h (h+1) ... (t-1)
NB.
NB. Notes:
NB. - monadic case is possible, though awkward:
NB.     _3 _2 _1 -: th2liso _3
NB.     5 4 3    -: th2liso  3

th2liso=: ] + i.@-

NB. ---------------------------------------------------------
NB. dhs2liso
NB.
NB. Description:
NB.   Generate lISO from head, size and delta
NB.
NB. Syntax:
NB.   liso=. [d] dhs2liso (h,s)
NB. where
NB.   h    - integer, head of liso, if h<0 then liso is
NB.          pointed to h, otherwise away from h
NB.   s    - integer, size of liso, if s<0 then liso's order
NB.          is reversed
NB.   d    ≥ 0 integer, optional delta of liso, default is 1
NB.   liso - |s|-vector of integers
NB.
NB. Examples:
NB.    2 dhs2liso 4 3                 2 dhs2liso _4 3
NB. 4 6 8                          _8 _6 _4
NB.    2 dhs2liso 4 _3                2 dhs2liso _4 _3
NB. 8 6 4                          _4 _6 _8
NB.
NB. Notes:
NB. - monadic case models rISO in (u;.0) with following
NB.   difference: s cannot be ±∞

dhs2liso=: 1&$: :({.@] + (negneg~ {.) * i.@(negneg/)@])

NB. ---------------------------------------------------------
NB. riso2iso
NB.
NB. Description:
NB.   Convert rISO to ISO
NB.
NB. Syntax:
NB.   iso=. riso2iso riso
NB.
NB. Notes:
NB. - riso with columns count less than array's rank is
NB.   indexing the slice

riso2iso=: <"1@riso2oiso

NB. ---------------------------------------------------------
NB. riso2liso
NB.
NB. Description:
NB.   Convert rISO to lISO
NB.
NB. Syntax:
NB.   liso=. sh riso2liso riso
NB. where
NB.   riso  - 2×r-array of integers, rISO subarray:
NB.             (2,r) $ from[0:r-1],size[0:r-1]
NB.   sh    - r-array of integers, shape of array to explore:
NB.             Size[0:r-1]
NB.   liso  - |Π{size[i],i=0:r-1}|-array of integers, rowwise
NB.           lISO subarray elements
NB.   r     ≥ 0, integer, rank of array to explore
NB.
NB. Formula:
NB.   liso[k] := Σ{Π{Size[j],j=i+1:r-1}*(n[k][i]-(n[k][i+1]<0 ? 1 : 0)),i=0:r-2} + n[k][r-1]
NB. where
NB.   k       = 0:|Π{size[i],i=0:r-1}|-1, IO liso' item
NB.   n[k][i] - i-th axis' IO for k-th liso' item
NB.
NB. Assertions:
NB.   (liso ({,) array) -: (riso ,;.0 array)
NB. where
NB.   riso=. 2 4 $ 7 _3 7 _3 2 2 _2 _2
NB.   sh=. 10 11 12 13
NB.   array=. i. sh
NB.   liso=. sh riso2liso riso

riso2liso=: */\.@(1&(|.!.1))@[ +/@:* (+ 1&(|.!.0)@(0&>))@|:@:>@,@{@riso2oiso@]

NB. ---------------------------------------------------------
NB. lisoE
NB. lisoN
NB. lisoS
NB. lisoW
NB.
NB. Description:
NB.   lISO vector laying between diagonal and matrix edge
NB.   in any of one cardinal direction: east, north, south or
NB.   west; and having optional gap between diagonal, at head
NB.   or tail
NB.
NB. Syntax:
NB.   vapp=. gap lisoX
NB. where
NB.   lisoX - adv., any of: lisoE lisoN lisoS lisoW
NB.   gap   - integer, negative value means "from
NB.           head", otherwise "from tail"
NB.   vapp  - dyad to return liso; is called as:
NB.             liso=. l vapp n
NB.   liso  - l-vector of integers, lISO v in ravelled A
NB.   v     - l-vector from A:
NB.             v -: liso ({,) A
NB.   A     - m×n-matrix
NB.
NB. Examples:
NB.    '***' ((((0 lisoE)&c)}),.' ',.(((1 lisoE)&c)}),.' ',.(((_1 lisoE)&c)})) 5 6$'-'
NB. ------ ------ ------
NB. ------ --***- ---***
NB. ---*** ------ ------
NB. ------ ------ ------
NB. ------ ------ ------
NB.    '***' ((((0 lisoN)&c)}),.' ',.(((1 lisoN)&c)}),.' ',.(((_1 lisoN)&c)})) 5 6$'-'
NB. --*--- ---*-- ------
NB. --*--- ---*-- ---*--
NB. --*--- ---*-- ---*--
NB. ------ ------ ---*--
NB. ------ ------ ------
NB.    '***' ((((0 lisoS)&c)}),.' ',.(((1 lisoS)&c)}),.' ',.(((_1 lisoS)&c)})) 5 6$'-'
NB. ------ ------ ------
NB. ------ --*--- ------
NB. ---*-- --*--- --*---
NB. ---*-- --*--- --*---
NB. ---*-- ------ --*---
NB.    '***' ((((0 lisoW)&c)}),.' ',.(((1 lisoW)&c)}),.' ',.(((_1 lisoW)&c)})) 5 6$'-'
NB. ------ ------ ------
NB. ------ ------ ------
NB. ***--- ------ ------
NB. ------ ***--- -***--
NB. ------ ------ ------

lisoE=: 1 : 'dhs2liso_mt_@(((_1 - 0 >. m) -  (* ((<: | m)&+))~) , [)'
lisoW=: 1 : 'dhs2liso_mt_@(((     0 <. m) -~ (* ((<: | m)&+))~) , [)'

lisoN=: 1 : '] dhs2liso_mt_ (((-~ ((<: | m)&+))~ (*&(0 <. m))) , [)'
lisoS=: 1 : '] dhs2liso_mt_ (((-~ ((-  | m)&-))~ (*&(0 >. m))) , [)'
