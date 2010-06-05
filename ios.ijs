NB. IOS
NB.
NB. liofmax    lIO 1st element with maximum sum of real and
NB.            imagine parts' modules
NB. liolmax    lIO last element with maximum sum of real and
NB.            imagine parts' modules
NB. lios2cp    Convert lIOS to cycle permutation
NB. th2lios    Generate lIOS from tail and head
NB. dhs2lios   Generate lIOS from head, size and optional
NB.            delta
NB. rios2ios   Convert rIOS to IOS
NB. rios2lios  Convert rIOS to lIOS
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Concepts
NB.
NB. IO   - index of
NB. lIO  - linear IO, is an integer
NB. IOS  - indices of
NB. lIOS - linear IOS, is a vector of integers
NB. rIOS - rectangular IOS, for r-rank array is a 2×r-array
NB.        of integers ((head0,head1,...),:(size0,size1,...))
NB.
NB. Following are equivalents:
NB.   (3 5 _7,:2 _3 4) (] ;. 0) report
NB.   (< 3 4;7 6 5;_10 _9 _8 _7) { report
NB.   (rios2ios (3 5 _7,:2 _3 4)) { report
NB.   (ios2rios (< 3 4;7 6 5;_10 _9 _8 _7)) (] ;. 0) report
NB.
NB. Following are equivalents:
NB.   (0 1 ; 1 2 ; 2 3) { i. 3 4
NB.   (4 lios2ios 1 6 11) { i. 3 4
NB.   (4 ios2lios (0 1 ; 1 2 ; 2 3)) ({,) i. 3 4
NB.   1 6 11 ({,) i. 3 4

NB. =========================================================
NB. Local definitions

NB. convert rIOS to opened (non-boxed) IOS
rios2oios=: < @ dhs2lios/ " 1 @: |:

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Miscellaneous

NB. lIO 1st element e with max(|Re(e)|+|Im(e)|) from list y
NB. implements BLAS's IxAMAX
liofmax=: (i.>./) @ sorim

NB. lIO last element e with max(|Re(e)|+|Im(e)|) from list y
liolmax=: (i:>./) @ sorim

NB. ---------------------------------------------------------
NB. lios2cp
NB.
NB. Description:
NB.   Convert lIOS to cycle permutation
NB.
NB. Syntax:
NB.   cp=. io0 lios2cp io1
NB. where
NB.   io0,io1 - lIOS
NB.   cp      - cycle permutation
NB.
NB. References:
NB. [1] [Jprogramming] Swapping array elements
NB.     Roger Hui, Mon May 11 06:23:07 HKT 2009
NB.     http://www.jsoftware.com/pipermail/programming/2009-May/014682.html

lios2cp=: < @ ~. @ ,

NB. ---------------------------------------------------------
NB. th2lios
NB.
NB. Description:
NB.   Generate lIOS from tail and head
NB.
NB. Syntax:
NB.   lios=. t th2lios h
NB. where
NB.   h    - integer, head of lios
NB.   t    - integer, tail of lios
NB.   lios - (x-y)-vector of integers, lIOS from head h to
NB.          tail (t-1) with delta=1:
NB.            h (h+1) ... (t-1)
NB.
NB. Notes:
NB. - monadic case is possible, though awkward:
NB.     _3 _2 _1 -: th2lios _3
NB.     5 4 3    -: th2lios  3

th2lios=: ] + (i. @ -)

NB. ---------------------------------------------------------
NB. dhs2lios
NB.
NB. Description:
NB.   Generate lIOS from head, size and optional delta
NB.
NB. Syntax:
NB.   lios=. [d] dhs2lios (h,s)
NB. where
NB.   h    - integer, head of lios, if h<0 then lios is
NB.          pointed to h, otherwise away from h
NB.   s    - integer, size of lios, if s<0 then lios's order
NB.          is reversed
NB.   d    - optional non-negative integer, delta of lios,
NB.          default is 1
NB.   lios - |s|-vector of integers, lIOS of solid part of
NB.          diagonal
NB.
NB. Notes:
NB. - monadic case models rIOS in (u;.0) with following
NB.   difference: s cannot be ±∞

dhs2lios=: 1&$: :({.@] + (condneg~ {.) * i.@(condneg/)@])

NB. ---------------------------------------------------------
NB. rios2ios
NB.
NB. Description:
NB.   Convert rIOS to IOS
NB.
NB. Syntax:
NB.   ios=. rios2ios rios
NB.
NB. Notes:
NB. - rios with columns count less than array's rank is
NB.   indexing the slice

rios2ios=: < " 1 @ rios2oios

NB. ---------------------------------------------------------
NB. rios2lios
NB.
NB. Description:
NB.   Convert rIOS to lIOS
NB.
NB. Syntax:
NB.   lios=. sh rios2lios rios
NB. where
NB.   rios  - 2×r-array of integers, rIOS of subarray:
NB.             (2,r) $ from[0:r-1],size[0:r-1]
NB.   sh    - r-array of integers, shape of array to explore:
NB.             Size[0:r-1]
NB.   lios  - |Π{size[i],i=0:r-1}|-array of integers, rowwise
NB.           lIOS of subarray elements
NB.   r     ≥ 0, integer, rank of array to explore
NB.
NB. Formula:
NB.   lios[k] := Σ{Π{Size[j],j=i+1:r-1}*(n[k][i]-(n[k][i+1]<0 ? 1 : 0)),i=0:r-2} + n[k][r-1]
NB. where
NB.   k       = 0:|Π{size[i],i=0:r-1}|-1, IO lios' item
NB.   n[k][i] - i-th axis' IO for k-th lios' item
NB.
NB. Assertion:
NB.   (lios ({,) array) -: (rios (, ;. 0) array)
NB. where
NB.   rios=. 2 4 $ 7 _3 7 _3 2 2 _2 _2
NB.   sh=. 10 11 12 13
NB.   array=. i. sh
NB.   lios=. sh rios2lios rios

rios2lios=: ((*/\.) @ (1 & (|.!.1)) @ [) (+/ @: *) ((+ (1 & (|.!.0) @ (0 & >))) @ |: @: > @ , @ { @ rios2oios @ ])
