NB. ISO
NB.
NB. liofmax    lIO 1st element with maximum sum of real and
NB.            imagine parts' modules
NB. liolmax    lIO last element with maximum sum of real and
NB.            imagine parts' modules
NB. liso4th    Generate lISO from tail and head
NB. liso4dhs   Generate lISO from head, size and optional
NB.            delta
NB. iso4riso   Convert rISO to/from ISO
NB. liso4riso  Convert rISO to lISO
NB. lisoX      Adv. to make dyad to compute lISO vector
NB.            laying between diagonal and matrix edge
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
NB. Concepts
NB.
NB. Notation:
NB.   IO   - index of
NB.   lIO  - linear IO, is an integer
NB.   ISO  - indices of
NB.   lISO - linear ISO, is a vector of integers
NB.   rISO - rectangular ISO, for r-rank array is a 2×r-array
NB.          of integers ((head0,head1,...),:(size0,size1,...))
NB.
NB. Assertions:
NB.   riso -: iso4riso^:_1 iso               NB. ISO converted back to the same rISO
NB.   (-: iso4riso^:_1@ iso4riso     ) riso  NB. rISO -> ISO -> back to the same rISO
NB.   (-: iso4riso    @(iso4riso^:_1)) iso   NB. ISO -> rISO -> back to the same ISO
NB.   (riso ];.0 arr) -: (iso { arr)         NB. riso and iso both are pointing to the same subarray
NB.   (riso ,;.0 arr) -: (liso ({,) arr)     NB. riso and liso both are pointing to the same subarray
NB. where
NB.   sh=. 4 # 10                            NB. some shape
NB.   arr=. i. sh                            NB. some array of shape sh is storing lISO values
NB.   riso=. 3 5 _6 _2 ,: 2 _3 4 _5          NB. some rISO within arr
NB.   iso=. iso4riso riso                    NB. ISO representation of riso
NB.   liso=. sh liso4riso riso               NB. lISO representation of riso

NB. =========================================================
NB. Configuration

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. convert rISO to opened (non-boxed) ISO
oiso4riso=: <@liso4dhs"1@:|:

NB. ---------------------------------------------------------
NB. lioxmax
NB.
NB. Description:
NB.   Conj. to make monad to locate extreme element e with
NB.   max(|Re(e)|+|Im(e)|) in numeric vector
NB.
NB. Syntax:
NB.   io=. (iof lioxmax xtr) vec
NB. where
NB.   iof - a monad to locate y in x, is either (i.) to
NB.         locate 1st element or (i:) to locate last one
NB.   xtr - a dyad to choose extreme element from x and y, is
NB.         either (<.) or (>.), is determined by iof
NB.   vec - a numeric vector
NB.   io  - a lIO of extreme element

lioxmax=: 2 : 0
  i=. (isnan y) u 1
  if. i < # y do.
    i return.
  end.
  if. JCMPX = 3!:0 y do.
    c=. |: | +. y
    i=. v/ c (u"1) _
    if. i < # y do.
      i return.
    end.
    s=. +/ c
    i=. s u _
    if. i < # y do.
      s=. +/ c % 4
    end.
  else.
    s=. | y
    i=. s u _
    if. i < # y do.
      i return.
    end.
  end.
  i=. (u >./) s
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad      Locates lIO
NB. liofmax    1st element
NB. liolmax    last element
NB.
NB. Description:
NB.   Locate extreme element e with
NB.   max(|Re(e)|+|Im(e)|) in numeric vector
NB.
NB. Syntax:
NB.   io=. lioxmax vec
NB.
NB. Assertions:
NB.   (liofmax -: <:@# - liolmax@|.) vec
NB.
NB. Notes:
NB. - liofmax implements BLAS' IxAMAX

liofmax=: i. lioxmax <.
liolmax=: i: lioxmax >.

NB. ---------------------------------------------------------
NB. liso4th
NB.
NB. Description:
NB.   Generate lISO from tail and head
NB.
NB. Syntax:
NB.   liso=. t liso4th h
NB. where
NB.   h    - integer, head of liso
NB.   t    ≥ h, integer, tail of liso
NB.   liso - (x-y)-vector of integers, lISO from head h to
NB.          tail (t-1) with delta=1:
NB.            h (h+1) ... (t-1)
NB.
NB. Examples:
NB.      3 liso4th 3           $ 3 liso4th 3
NB.                         0
NB.      4 liso4th 3           5 liso4th 3
NB.   3                     3 4

liso4th=: [: : (] + i.@-)

NB. ---------------------------------------------------------
NB. liso4dhs
NB.
NB. Description:
NB.   Generate lISO from head, size and delta
NB.
NB. Syntax:
NB.   liso=. [d] liso4dhs (h,s)
NB. where
NB.   h    - integer, head of liso, if h<0 then liso is
NB.          pointed to h, otherwise away from h
NB.   s    - integer, size of liso, if s<0 then liso's order
NB.          is reversed
NB.   d    ≥ 0 integer, optional delta of liso, default is 1
NB.   liso - |s|-vector of integers
NB.
NB. Examples:
NB.      liso4dhs 4 3              1 liso4dhs 4 3
NB.   4 5 6                     4 5 6
NB.      2 liso4dhs 4 3            2 liso4dhs _4 3
NB.   4 6 8                     _8 _6 _4
NB.      2 liso4dhs 4 _3           2 liso4dhs _4 _3
NB.   8 6 4                     _4 _6 _8
NB.
NB. Notes:
NB. - monadic case models rISO in (u;.0) with following
NB.   difference: s cannot be ±∞

liso4dhs=: 1&$: :({.@] + (negneg~ {.) * i.@(negneg/)@])

NB. ---------------------------------------------------------
NB. iso4riso
NB.
NB. Description:
NB.   Convert rISO to/from ISO
NB.
NB. Syntax:
NB.   iso=.  iso4riso     riso
NB.   riso=. iso4riso^:_1 iso
NB.
NB. Assertions:
NB.   iso  -: iso4riso     riso
NB.   riso -: iso4riso^:_1 iso
NB. where
NB.      ] riso=. 3 5 _6 _2 ,: 2 _3 4 _5
NB.   3  5 _6 _2
NB.   2 _3  4 _5
NB.      ] iso=. < 3 4 ; 7 6 5 ; _9 _8 _7 _6 ; _2 _3 _4 _5 _6
NB.   +--------------------------------------+
NB.   |+---+-----+-----------+--------------+|
NB.   ||3 4|7 6 5|_9 _8 _7 _6|_2 _3 _4 _5 _6||
NB.   |+---+-----+-----------+--------------+|
NB.   +--------------------------------------+
NB.
NB. Notes:
NB. - riso with columns count less than array's rank is
NB.   indexing the slice

iso4riso=: <"1@oiso4riso :. (3 : 0)
  'not representable as rISO' assert (('' -: -.&_1 1) *. 2 > #@~.)S:0 (2&(-/\)^:(0 < #))L:0 y
    NB. each open element is either (1) empty list or (2) scalar or (3) range which is (3a) continuous and either (3b1) increasing or (3b2) decreasing
  dm=. (-/@:|@(2&{.)`_1:@.(2 > #)       L:0) y  NB. deltas of modules
  d=.  (-/   @(2&{.)`_1:@.(2 > #)       L:0) y  NB. deltas
  h=. dm ([:`(_1 { ::0 ])`(0 { ::0 ])@.[S:0) y  NB. heads
  s=. d  ([:`(-@#@]     )`(#@]      )@.[S:0) y  NB. sizes
  h ,: s
)

NB. ---------------------------------------------------------
NB. liso4riso
NB.
NB. Description:
NB.   Convert rISO to lISO
NB.
NB. Syntax:
NB.   liso=. sh liso4riso riso
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
NB.   (riso ,;.0 arr) -: (liso ({,) arr)
NB. where
NB.   sh=. 4 # 10                    NB. some shape
NB.   arr=. sh ?@$ 0                 NB. some array of shape sh
NB.   riso=. 3 5 _6 _2 ,: 2 _3 4 _5  NB. some rISO within arr
NB.   liso=. sh liso4riso riso       NB. lISO representation of riso

liso4riso=: */\.@(1&(|.!.1))@[ +/@:* (+ 1&(|.!.0)@(0&>))@|:@:>@,@{@oiso4riso@]

NB. ---------------------------------------------------------
NB. lisoE
NB. lisoN
NB. lisoS
NB. lisoW
NB.
NB. Description:
NB.   Adv. to make dyad to compute lISO vector laying between
NB.   diagonal and matrix edge in any of one cardinal
NB.   direction: east, north, south or west; and having
NB.   optional gap to diagonal, at head or tail
NB.
NB. Syntax:
NB.   liso=. l (gap lisoE) n
NB.   liso=. l (gap lisoN) n
NB.   liso=. l (gap lisoS) n
NB.   liso=. l (gap lisoW) n
NB. where
NB.   gap  - integer, negative value means "from
NB.          head", otherwise "from tail"
NB.   liso - l-vector of integers, lISO v in ravelled A
NB.   v    - l-vector from A:
NB.            v -: liso ({,) A
NB.   A    - m×n-matrix
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

lisoE=: 1 : 'liso4dhs_mt_@(((_1 - 0 >. m) -  (* ((<: | m)&+))~) , [)'
lisoW=: 1 : 'liso4dhs_mt_@(((     0 <. m) -~ (* ((<: | m)&+))~) , [)'

lisoN=: 1 : '] liso4dhs_mt_ (((-~ ((<: | m)&+))~ (*&(0 <. m))) , [)'
lisoS=: 1 : '] liso4dhs_mt_ (((-~ ((-  | m)&-))~ (*&(0 >. m))) , [)'
