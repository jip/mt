NB. util.ijs
NB. Utilities
NB.
NB. sgn        Simplified signum
NB. condneg    Conditional negate
NB. copysign   Copy sign
NB. fmtlog     Format log string
NB. gi         Conj. to evoke n-th verb from gerund m
NB. iofmax     IO 1st element with maximum sum of real and
NB.            imagine parts' modules
NB. iolmax     IO last element with maximum sum of real and
NB.            imagine parts' modules
NB. ms         Minimum in sum of vectors
NB. ios2cp     Make cycle permutation from indices
NB.
NB. norm1      Magnitude-based 1-norm of vector or matrix
NB. normi      Magnitude-based ∞-norm of vector or matrix
NB. norm1t     Taxicab-based 1-norm of vector or matrix
NB. normit     Taxicab-based ∞-norm of vector or matrix
NB. norms      Square-based (Euclidean/Frobenius) norm of
NB.            vector or matrix
NB.
NB. hds2ios    Form IOS from head, delta and size
NB. ht2ios     Form IOS from head and tail
NB. hs2ios     Form IOS from head and size
NB. rios2ios   Convert rIOS to IOS
NB. rios2lios  Convert rIOS to lIOS
NB.
NB. tmonad     Template conj. to make verbs to test
NB.            computational monad
NB. tdyad      Template conj. to make verbs to test
NB.            computational dyad
NB. dbg        Conj. to show verb's input and output
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

rios2oios=: < @ hs2ios/ " 1 @: |:        NB. convert rIOS to opened (non-boxed) IOS
gshapes=: $`($(;<)($ L: 0)) @. (0 < L.)  NB. get shapes, boxes are accepted, too
sorim=: +/"1 @: | @: +.                  NB. sum of real and imaginary parts' modules
mocs=: >./ @ (+/   ) @:                  NB. vector: sum of, matrix: max of column sums
mors=: >./ @ (+/"_1) @:                  NB. vector: max of, matrix: max of row sums

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Miscellaneous

sgn=: 0&(<: - >)                                        NB. if y<0 then -1 else 1 endif
condneg=: -@]^:(0>[)                                    NB. if x<0 then -y else y endif
copysign=: -@]^:((=-)&*)                                NB. if x<0 then -|y| else |y| endif
fmtlog=: '%-25S %-12g %-12g %-12g %-12g %12d' vsprintf  NB. Format log string
gi=: 2 : '(n{m)`:6'                                     NB. Conj. to evoke n-th verb from gerund m: m[n]

NB. IO 1st element e with max(|Re(e)|+|Im(e)|) from list y
NB. emulates LAPACK's IxAMAX
iofmax=: (i.>./) @ sorim

NB. IO last element e with max(|Re(e)|+|Im(e)|) from list y
iolmax=: (i:>./) @ sorim

NB. ---------------------------------------------------------
NB. ms
NB.
NB. Description: Minimum in [sum of] vector[s]
NB. Syntax:      k=. [(delta0,delta1,...)] ms (value0,value1,...)
NB. where        default deltai is 0
NB. Formula:     k = min(delta0+value0,delta1+value1,...)
NB. Notes:       is memo, since repetitive calls are expected

ms=: <./@:(] :+)M.

NB. ---------------------------------------------------------
NB. ios2cp
NB.
NB. Description:
NB.   Make cycle permutation from indices x and y
NB.
NB. Syntax:
NB.   cp=. io0 ios2cp io1
NB. where
NB.   io0 io1 - IOS
NB.   cp      - cycle permutation
NB.
NB. References:
NB. [1] [Jprogramming] Swapping array elements
NB.     Roger Hui, Mon May 11 06:23:07 HKT 2009
NB.     http://www.jsoftware.com/pipermail/programming/2009-May/014682.html
NB.
NB. TODO:
NB. - if y is empty then output must be a:
NB.
NB. Consider:
NB.   ios2cp=: < @ (, ` (, @ ]) @. =)
NB.   ios2cp=: < @ , @ (, ^: ~:)
NB.   ios2cp=: (, ^: (-. @ -:)) & ,

ios2cp=: < @ ~. @ ,

NB. ---------------------------------------------------------
NB. Norms

NB. Magnitude-based norms |y|
norm1=: | mocs       NB. 1-norm of vector or matrix
normi=: | mors       NB. ∞-norm of vector or matrix

NB. Taxicab-based norms |Re(y)| + |Im(y)|
norm1t=: sorim mocs  NB. 1-norm of vector or matrix
normit=: sorim mors  NB. ∞-norm of vector or matrix

NB. Square-based (Euclidean/Frobenius) norm of vector or matrix
NB. for vector input emulates LAPACK's DZNRM2
norms=: (((((+/^:_) &.: *:) @: %) * ]) (>./^:_)) @: | @: +.

NB. ---------------------------------------------------------
NB. IOS generators

NB. y[2]-vector of integers from head y[0] by delta y[1]
hds2ios=: + ` (* i.)/

NB. (x-y)-vector of integers from head y to tail (x-1):
NB. y (y+1) ... (x-1)
NB. monadic case is possible:
NB.   _3 _2 _1 -: ht2ios _3
NB.   5 4 3    -: ht2ios  3
ht2ios=: ] + (i. @ -)

NB. y-vector of integers from head x; models rIOS in (u;.0)
NB. TODO:
NB. - current_hs2ios ===      hs2ios (head,size)
NB.                  === 1    hs2ios (head,size)
NB. - riosed_hds2ios === step hs2ios (head,size)

hs2ios=: [ + ((] * i. @ *) sgn)~

NB. lios=. [delta] dhs2ios (head,size)
dhs2ios1=: ({. + ({: (] * i. @ *) (sgn @ {.))) : ({.@] + (* ({: (]*i.@*) (sgn@{.))))
dhs2ios2=: (1&$:) : ({.@] + (* ({: (] * i. @ *) (sgn @ {.))))
dhs2ios3=: ({. + ({: (] * i. @ *) (sgn @ {.))) : ({.@] + (*^:(1~:[) ({: (] * i. @ *) (sgn @ {.))))
dhs2ios4=: (1&$:) : ({.@] + (*^:(1~:[) ({: (] * i. @ *) (sgn @ {.))))
dhs2ios5=: ({. + ({. condneg ((i. @ condneg)/))) : ({.@] + ({.@] condneg (* ((i. @ condneg)/))))
dhs2ios6=: (1&$:) : ({.@] + ({.@] condneg (* ((i. @ condneg)/))))
dhs2ios7=: ({. + ({. condneg ((i. @ condneg)/))) : ({.@] + ({.@] condneg (*^:(1~:[) ((i. @ condneg)/))))
dhs2ios8=: (1&$:) : ({.@] + ({.@] condneg (*^:(1~:[) ((i. @ condneg)/))))
dhs2ios9=: ({. + ({. condneg ((i. @ condneg)/))) : ({.@] + (condneg~ {.) * i.@(condneg/)@])
dhs2ios10=: (1&$:) : ({.@] + (condneg~ {.) * i.@(condneg/)@])

testdhs2ios=: 3 : 0
  'hpsp hpsn hnsp hnsn'=. (1 1 , 1 _1 , _1 1 ,: _1 _1) (*"1) 1 {:: 'd hs'=. ({. ; }.) y
  e=. i. 0

  smoutput 'dhs2ios2' ; e ; hpsp ; (  (dhs2ios1 -: dhs2ios2) hpsp)
  smoutput 'dhs2ios2' ; e ; hpsn ; (  (dhs2ios1 -: dhs2ios2) hpsn)
  smoutput 'dhs2ios2' ; e ; hnsp ; (  (dhs2ios1 -: dhs2ios2) hnsp)
  smoutput 'dhs2ios2' ; e ; hnsn ; (  (dhs2ios1 -: dhs2ios2) hnsn)

  smoutput 'dhs2ios2' ; 1 ; hpsp ; (1 (dhs2ios1 -: dhs2ios2) hpsp)
  smoutput 'dhs2ios2' ; 1 ; hpsn ; (1 (dhs2ios1 -: dhs2ios2) hpsn)
  smoutput 'dhs2ios2' ; 1 ; hnsp ; (1 (dhs2ios1 -: dhs2ios2) hnsp)
  smoutput 'dhs2ios2' ; 1 ; hnsn ; (1 (dhs2ios1 -: dhs2ios2) hnsn)

  smoutput 'dhs2ios2' ; d ; hpsp ; (1 (dhs2ios1 -: dhs2ios2) hpsp)
  smoutput 'dhs2ios2' ; d ; hpsn ; (1 (dhs2ios1 -: dhs2ios2) hpsn)
  smoutput 'dhs2ios2' ; d ; hnsp ; (1 (dhs2ios1 -: dhs2ios2) hnsp)
  smoutput 'dhs2ios2' ; d ; hnsn ; (1 (dhs2ios1 -: dhs2ios2) hnsn)


  smoutput 'dhs2ios3' ; e ; hpsp ; (  (dhs2ios1 -: dhs2ios3) hpsp)
  smoutput 'dhs2ios3' ; e ; hpsn ; (  (dhs2ios1 -: dhs2ios3) hpsn)
  smoutput 'dhs2ios3' ; e ; hnsp ; (  (dhs2ios1 -: dhs2ios3) hnsp)
  smoutput 'dhs2ios3' ; e ; hnsn ; (  (dhs2ios1 -: dhs2ios3) hnsn)

  smoutput 'dhs2ios3' ; 1 ; hpsp ; (1 (dhs2ios1 -: dhs2ios3) hpsp)
  smoutput 'dhs2ios3' ; 1 ; hpsn ; (1 (dhs2ios1 -: dhs2ios3) hpsn)
  smoutput 'dhs2ios3' ; 1 ; hnsp ; (1 (dhs2ios1 -: dhs2ios3) hnsp)
  smoutput 'dhs2ios3' ; 1 ; hnsn ; (1 (dhs2ios1 -: dhs2ios3) hnsn)

  smoutput 'dhs2ios3' ; d ; hpsp ; (1 (dhs2ios1 -: dhs2ios3) hpsp)
  smoutput 'dhs2ios3' ; d ; hpsn ; (1 (dhs2ios1 -: dhs2ios3) hpsn)
  smoutput 'dhs2ios3' ; d ; hnsp ; (1 (dhs2ios1 -: dhs2ios3) hnsp)
  smoutput 'dhs2ios3' ; d ; hnsn ; (1 (dhs2ios1 -: dhs2ios3) hnsn)


  smoutput 'dhs2ios4' ; e ; hpsp ; (  (dhs2ios1 -: dhs2ios4) hpsp)
  smoutput 'dhs2ios4' ; e ; hpsn ; (  (dhs2ios1 -: dhs2ios4) hpsn)
  smoutput 'dhs2ios4' ; e ; hnsp ; (  (dhs2ios1 -: dhs2ios4) hnsp)
  smoutput 'dhs2ios4' ; e ; hnsn ; (  (dhs2ios1 -: dhs2ios4) hnsn)

  smoutput 'dhs2ios4' ; 1 ; hpsp ; (1 (dhs2ios1 -: dhs2ios4) hpsp)
  smoutput 'dhs2ios4' ; 1 ; hpsn ; (1 (dhs2ios1 -: dhs2ios4) hpsn)
  smoutput 'dhs2ios4' ; 1 ; hnsp ; (1 (dhs2ios1 -: dhs2ios4) hnsp)
  smoutput 'dhs2ios4' ; 1 ; hnsn ; (1 (dhs2ios1 -: dhs2ios4) hnsn)

  smoutput 'dhs2ios4' ; d ; hpsp ; (1 (dhs2ios1 -: dhs2ios4) hpsp)
  smoutput 'dhs2ios4' ; d ; hpsn ; (1 (dhs2ios1 -: dhs2ios4) hpsn)
  smoutput 'dhs2ios4' ; d ; hnsp ; (1 (dhs2ios1 -: dhs2ios4) hnsp)
  smoutput 'dhs2ios4' ; d ; hnsn ; (1 (dhs2ios1 -: dhs2ios4) hnsn)


  smoutput 'dhs2ios5' ; e ; hpsp ; (  (dhs2ios1 -: dhs2ios5) hpsp)
  smoutput 'dhs2ios5' ; e ; hpsn ; (  (dhs2ios1 -: dhs2ios5) hpsn)
  smoutput 'dhs2ios5' ; e ; hnsp ; (  (dhs2ios1 -: dhs2ios5) hnsp)
  smoutput 'dhs2ios5' ; e ; hnsn ; (  (dhs2ios1 -: dhs2ios5) hnsn)

  smoutput 'dhs2ios5' ; 1 ; hpsp ; (1 (dhs2ios1 -: dhs2ios5) hpsp)
  smoutput 'dhs2ios5' ; 1 ; hpsn ; (1 (dhs2ios1 -: dhs2ios5) hpsn)
  smoutput 'dhs2ios5' ; 1 ; hnsp ; (1 (dhs2ios1 -: dhs2ios5) hnsp)
  smoutput 'dhs2ios5' ; 1 ; hnsn ; (1 (dhs2ios1 -: dhs2ios5) hnsn)

  smoutput 'dhs2ios5' ; d ; hpsp ; (1 (dhs2ios1 -: dhs2ios5) hpsp)
  smoutput 'dhs2ios5' ; d ; hpsn ; (1 (dhs2ios1 -: dhs2ios5) hpsn)
  smoutput 'dhs2ios5' ; d ; hnsp ; (1 (dhs2ios1 -: dhs2ios5) hnsp)
  smoutput 'dhs2ios5' ; d ; hnsn ; (1 (dhs2ios1 -: dhs2ios5) hnsn)


  smoutput 'dhs2ios6' ; e ; hpsp ; (  (dhs2ios1 -: dhs2ios6) hpsp)
  smoutput 'dhs2ios6' ; e ; hpsn ; (  (dhs2ios1 -: dhs2ios6) hpsn)
  smoutput 'dhs2ios6' ; e ; hnsp ; (  (dhs2ios1 -: dhs2ios6) hnsp)
  smoutput 'dhs2ios6' ; e ; hnsn ; (  (dhs2ios1 -: dhs2ios6) hnsn)

  smoutput 'dhs2ios6' ; 1 ; hpsp ; (1 (dhs2ios1 -: dhs2ios6) hpsp)
  smoutput 'dhs2ios6' ; 1 ; hpsn ; (1 (dhs2ios1 -: dhs2ios6) hpsn)
  smoutput 'dhs2ios6' ; 1 ; hnsp ; (1 (dhs2ios1 -: dhs2ios6) hnsp)
  smoutput 'dhs2ios6' ; 1 ; hnsn ; (1 (dhs2ios1 -: dhs2ios6) hnsn)

  smoutput 'dhs2ios6' ; d ; hpsp ; (1 (dhs2ios1 -: dhs2ios6) hpsp)
  smoutput 'dhs2ios6' ; d ; hpsn ; (1 (dhs2ios1 -: dhs2ios6) hpsn)
  smoutput 'dhs2ios6' ; d ; hnsp ; (1 (dhs2ios1 -: dhs2ios6) hnsp)
  smoutput 'dhs2ios6' ; d ; hnsn ; (1 (dhs2ios1 -: dhs2ios6) hnsn)


  smoutput 'dhs2ios7' ; e ; hpsp ; (  (dhs2ios1 -: dhs2ios7) hpsp)
  smoutput 'dhs2ios7' ; e ; hpsn ; (  (dhs2ios1 -: dhs2ios7) hpsn)
  smoutput 'dhs2ios7' ; e ; hnsp ; (  (dhs2ios1 -: dhs2ios7) hnsp)
  smoutput 'dhs2ios7' ; e ; hnsn ; (  (dhs2ios1 -: dhs2ios7) hnsn)

  smoutput 'dhs2ios7' ; 1 ; hpsp ; (1 (dhs2ios1 -: dhs2ios7) hpsp)
  smoutput 'dhs2ios7' ; 1 ; hpsn ; (1 (dhs2ios1 -: dhs2ios7) hpsn)
  smoutput 'dhs2ios7' ; 1 ; hnsp ; (1 (dhs2ios1 -: dhs2ios7) hnsp)
  smoutput 'dhs2ios7' ; 1 ; hnsn ; (1 (dhs2ios1 -: dhs2ios7) hnsn)

  smoutput 'dhs2ios7' ; d ; hpsp ; (1 (dhs2ios1 -: dhs2ios7) hpsp)
  smoutput 'dhs2ios7' ; d ; hpsn ; (1 (dhs2ios1 -: dhs2ios7) hpsn)
  smoutput 'dhs2ios7' ; d ; hnsp ; (1 (dhs2ios1 -: dhs2ios7) hnsp)
  smoutput 'dhs2ios7' ; d ; hnsn ; (1 (dhs2ios1 -: dhs2ios7) hnsn)


  smoutput 'dhs2ios8' ; e ; hpsp ; (  (dhs2ios1 -: dhs2ios8) hpsp)
  smoutput 'dhs2ios8' ; e ; hpsn ; (  (dhs2ios1 -: dhs2ios8) hpsn)
  smoutput 'dhs2ios8' ; e ; hnsp ; (  (dhs2ios1 -: dhs2ios8) hnsp)
  smoutput 'dhs2ios8' ; e ; hnsn ; (  (dhs2ios1 -: dhs2ios8) hnsn)

  smoutput 'dhs2ios8' ; 1 ; hpsp ; (1 (dhs2ios1 -: dhs2ios8) hpsp)
  smoutput 'dhs2ios8' ; 1 ; hpsn ; (1 (dhs2ios1 -: dhs2ios8) hpsn)
  smoutput 'dhs2ios8' ; 1 ; hnsp ; (1 (dhs2ios1 -: dhs2ios8) hnsp)
  smoutput 'dhs2ios8' ; 1 ; hnsn ; (1 (dhs2ios1 -: dhs2ios8) hnsn)

  smoutput 'dhs2ios8' ; d ; hpsp ; (1 (dhs2ios1 -: dhs2ios8) hpsp)
  smoutput 'dhs2ios8' ; d ; hpsn ; (1 (dhs2ios1 -: dhs2ios8) hpsn)
  smoutput 'dhs2ios8' ; d ; hnsp ; (1 (dhs2ios1 -: dhs2ios8) hnsp)
  smoutput 'dhs2ios8' ; d ; hnsn ; (1 (dhs2ios1 -: dhs2ios8) hnsn)


  smoutput 'dhs2ios9' ; e ; hpsp ; (  (dhs2ios1 -: dhs2ios9) hpsp)
  smoutput 'dhs2ios9' ; e ; hpsn ; (  (dhs2ios1 -: dhs2ios9) hpsn)
  smoutput 'dhs2ios9' ; e ; hnsp ; (  (dhs2ios1 -: dhs2ios9) hnsp)
  smoutput 'dhs2ios9' ; e ; hnsn ; (  (dhs2ios1 -: dhs2ios9) hnsn)

  smoutput 'dhs2ios9' ; 1 ; hpsp ; (1 (dhs2ios1 -: dhs2ios9) hpsp)
  smoutput 'dhs2ios9' ; 1 ; hpsn ; (1 (dhs2ios1 -: dhs2ios9) hpsn)
  smoutput 'dhs2ios9' ; 1 ; hnsp ; (1 (dhs2ios1 -: dhs2ios9) hnsp)
  smoutput 'dhs2ios9' ; 1 ; hnsn ; (1 (dhs2ios1 -: dhs2ios9) hnsn)

  smoutput 'dhs2ios9' ; d ; hpsp ; (1 (dhs2ios1 -: dhs2ios9) hpsp)
  smoutput 'dhs2ios9' ; d ; hpsn ; (1 (dhs2ios1 -: dhs2ios9) hpsn)
  smoutput 'dhs2ios9' ; d ; hnsp ; (1 (dhs2ios1 -: dhs2ios9) hnsp)
  smoutput 'dhs2ios9' ; d ; hnsn ; (1 (dhs2ios1 -: dhs2ios9) hnsn)


  smoutput 'dhs2ios10' ; e ; hpsp ; (  (dhs2ios1 -: dhs2ios10) hpsp)
  smoutput 'dhs2ios10' ; e ; hpsn ; (  (dhs2ios1 -: dhs2ios10) hpsn)
  smoutput 'dhs2ios10' ; e ; hnsp ; (  (dhs2ios1 -: dhs2ios10) hnsp)
  smoutput 'dhs2ios10' ; e ; hnsn ; (  (dhs2ios1 -: dhs2ios10) hnsn)

  smoutput 'dhs2ios10' ; 1 ; hpsp ; (1 (dhs2ios1 -: dhs2ios10) hpsp)
  smoutput 'dhs2ios10' ; 1 ; hpsn ; (1 (dhs2ios1 -: dhs2ios10) hpsn)
  smoutput 'dhs2ios10' ; 1 ; hnsp ; (1 (dhs2ios1 -: dhs2ios10) hnsp)
  smoutput 'dhs2ios10' ; 1 ; hnsn ; (1 (dhs2ios1 -: dhs2ios10) hnsn)

  smoutput 'dhs2ios10' ; d ; hpsp ; (1 (dhs2ios1 -: dhs2ios10) hpsp)
  smoutput 'dhs2ios10' ; d ; hpsn ; (1 (dhs2ios1 -: dhs2ios10) hpsn)
  smoutput 'dhs2ios10' ; d ; hnsp ; (1 (dhs2ios1 -: dhs2ios10) hnsp)
  smoutput 'dhs2ios10' ; d ; hnsn ; (1 (dhs2ios1 -: dhs2ios10) hnsn)

  EMPTY
)

timingdhs2ios=: 3 : 0

  NB. output: verb;(time,space);IO_time_graded_down;IO_space_graded_down

  NB. timing: head,size

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios00a=. 5 hs2ios 100000';'ios1a=. dhs2ios1 5 100000';'ios2a=. dhs2ios2 5 100000';'ios3a=. dhs2ios3 5 100000';'ios4a=. dhs2ios4 5 100000';'ios5a=. dhs2ios5 5 100000';'ios6a=. dhs2ios6 5 100000';'ios7a=. dhs2ios7 5 100000';'ios8a=. dhs2ios8 5 100000';'ios9a=. dhs2ios9 5 100000';'ios10a=. dhs2ios10 5 100000'

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios00b=. 5 hs2ios _100000';'ios1b=. dhs2ios1 5 _100000';'ios2b=. dhs2ios2 5 _100000';'ios3b=. dhs2ios3 5 _100000';'ios4b=. dhs2ios4 5 _100000';'ios5b=. dhs2ios5 5 _100000';'ios6b=. dhs2ios6 5 _100000';'ios7b=. dhs2ios7 5 _100000';'ios8b=. dhs2ios8 5 _100000';'ios9b=. dhs2ios9 5 _100000';'ios10b=. dhs2ios10 5 _100000'

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios00c=. _5 hs2ios 100000';'ios1c=. dhs2ios1 _5 100000';'ios2c=. dhs2ios2 _5 100000';'ios3c=. dhs2ios3 _5 100000';'ios4c=. dhs2ios4 _5 100000';'ios5c=. dhs2ios5 _5 100000';'ios6c=. dhs2ios6 _5 100000';'ios7c=. dhs2ios7 _5 100000';'ios8c=. dhs2ios8 _5 100000';'ios9c=. dhs2ios9 _5 100000';'ios10c=. dhs2ios10 _5 100000'

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios00d=. _5 hs2ios _100000';'ios1d=. dhs2ios1 _5 _100000';'ios2d=. dhs2ios2 _5 _100000';'ios3d=. dhs2ios3 _5 _100000';'ios4d=. dhs2ios4 _5 _100000';'ios5d=. dhs2ios5 _5 _100000';'ios6d=. dhs2ios6 _5 _100000';'ios7d=. dhs2ios7 _5 _100000';'ios8d=. dhs2ios8 _5 _100000';'ios9d=. dhs2ios9 _5 _100000';'ios10d=. dhs2ios10 _5 _100000'

  NB. timing: head,1,size

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios0e=. hds2ios 5 1 100000';'ios1e=. 1 dhs2ios1 5 100000';'ios2e=. 1 dhs2ios2 5 100000';'ios3e=. 1 dhs2ios3 5 100000';'ios4e=. 1 dhs2ios4 5 100000';'ios5e=. 1 dhs2ios5 5 100000';'ios6e=. 1 dhs2ios6 5 100000';'ios7e=. 1 dhs2ios7 5 100000';'ios8e=. 1 dhs2ios8 5 100000';'ios9e=. 1 dhs2ios9 5 100000';'ios10e=. 1 dhs2ios10 5 100000'

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios1f=. 1 dhs2ios1 5 _100000';'ios2f=. 1 dhs2ios2 5 _100000';'ios3f=. 1 dhs2ios3 5 _100000';'ios4f=. 1 dhs2ios4 5 _100000';'ios5f=. 1 dhs2ios5 5 _100000';'ios6f=. 1 dhs2ios6 5 _100000';'ios7f=. 1 dhs2ios7 5 _100000';'ios8f=. 1 dhs2ios8 5 _100000';'ios9f=. 1 dhs2ios9 5 _100000';'ios10f=. 1 dhs2ios10 5 _100000'

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios0g=. hds2ios _5 1 100000';'ios1g=. 1 dhs2ios1 _5 100000';'ios2g=. 1 dhs2ios2 _5 100000';'ios3g=. 1 dhs2ios3 _5 100000';'ios4g=. 1 dhs2ios4 _5 100000';'ios5g=. 1 dhs2ios5 _5 100000';'ios6g=. 1 dhs2ios6 _5 100000';'ios7g=. 1 dhs2ios7 _5 100000';'ios8g=. 1 dhs2ios8 _5 100000';'ios9g=. 1 dhs2ios9 _5 100000';'ios10g=. 1 dhs2ios10 _5 100000'

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios1h=. 1 dhs2ios1 _5 _100000';'ios2h=. 1 dhs2ios2 _5 _100000';'ios3h=. 1 dhs2ios3 _5 _100000';'ios4h=. 1 dhs2ios4 _5 _100000';'ios5h=. 1 dhs2ios5 _5 _100000';'ios6h=. 1 dhs2ios6 _5 _100000';'ios7h=. 1 dhs2ios7 _5 _100000';'ios8h=. 1 dhs2ios8 _5 _100000';'ios9h=. 1 dhs2ios9 _5 _100000';'ios10h=. 1 dhs2ios10 _5 _100000'

  NB. timing: head,10,size

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios0i=. hds2ios 5 10 100000';'ios1i=. 10 dhs2ios1 5 100000';'ios2i=. 10 dhs2ios2 5 100000';'ios3i=. 10 dhs2ios3 5 100000';'ios4i=. 10 dhs2ios4 5 100000';'ios5i=. 10 dhs2ios5 5 100000';'ios6i=. 10 dhs2ios6 5 100000';'ios7i=. 10 dhs2ios7 5 100000';'ios8i=. 10 dhs2ios8 5 100000';'ios9i=. 10 dhs2ios9 5 100000';'ios10i=. 10 dhs2ios10 5 100000'

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios1j=. 10 dhs2ios1 5 _100000';'ios2j=. 10 dhs2ios2 5 _100000';'ios3j=. 10 dhs2ios3 5 _100000';'ios4j=. 10 dhs2ios4 5 _100000';'ios5j=. 10 dhs2ios5 5 _100000';'ios6j=. 10 dhs2ios6 5 _100000';'ios7j=. 10 dhs2ios7 5 _100000';'ios8j=. 10 dhs2ios8 5 _100000';'ios9j=. 10 dhs2ios9 5 _100000';'ios10j=. 10 dhs2ios10 5 _100000'

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios0k=. hds2ios _5 10 100000';'ios1k=. 10 dhs2ios1 _5 100000';'ios2k=. 10 dhs2ios2 _5 100000';'ios3k=. 10 dhs2ios3 _5 100000';'ios4k=. 10 dhs2ios4 _5 100000';'ios5k=. 10 dhs2ios5 _5 100000';'ios6k=. 10 dhs2ios6 _5 100000';'ios7k=. 10 dhs2ios7 _5 100000';'ios8k=. 10 dhs2ios8 _5 100000';'ios9k=. 10 dhs2ios9 _5 100000';'ios10k=. 10 dhs2ios10 _5 100000'

  smoutput (,. (<"0@((/:"1)&.|:)@:>@:({:"1)))@:((; (40&timespacex))&>) 'ios1l=. dhs2ios1 _5 _100000';'ios2l=. dhs2ios2 _5 _100000';'ios3l=. dhs2ios3 _5 _100000';'ios4l=. dhs2ios4 _5 _100000';'ios5l=. dhs2ios5 _5 _100000';'ios6l=. dhs2ios6 _5 _100000';'ios7l=. dhs2ios7 _5 _100000';'ios8l=. dhs2ios8 _5 _100000';'ios9l=. dhs2ios9 _5 _100000';'ios10l=. dhs2ios10 _5 _100000'

  EMPTY
)

NB. ---------------------------------------------------------
NB. IOS converters

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Convert rIOS to IOS
NB. Note: IOS with length less than array's rank indexes
NB.       the slice

rios2ios=: < " 1 @ rios2oios

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. rios2lios
NB. Convert rIOS to lIOS
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
NB. Formulae:
NB.   lios[k] := Σ{Π{Size[j],j=i+1:r-1}*(n[k][i]-(n[k][i+1]<0 ? 1 : 0)),i=0:r-2} + n[k][r-1]
NB. where
NB.   k       = 0:|Π{size[i],i=0:r-1}|-1, IO lios' item
NB.   n[k][i] - i-th axis' IO for k-th lios' item
NB.
NB. If:
NB.   rios=. 2 4 $ 7 _3 7 _3 2 2 _2 _2
NB.   sh=. 10 11 12 13
NB.   array=. i. sh
NB.   lios=. sh rios2lios rios
NB. then
NB.   (lios ({,) array) -: (rios (, ;. 0) array)

rios2lios=: ((*/\.) @ (1 & (|.!.1)) @ [) (+/ @: *) ((+ (1 & (|.!.0) @ (0 & >))) @ |: @: > @ , @ { @ rios2oios @ ])

NB. ---------------------------------------------------------
NB. tmonad
NB. tdyad
NB. Template conj. to make monad to test computational verb
NB.
NB. Syntax:
NB.   vtestm=. mname tmonad vgety`vgeto`vrcond`vferr`vberr
NB.   vtestd=. dname tdyad vgetx`vgety`vgeto`vrcond`vferr`vberr
NB. where
NB.   vgetx  - monad to extract left argument for vd; is
NB.            called as:
NB.              argx=. vgetx y
NB.   vgety  - monad to extract right argument for vm or vd;
NB.            is called as:
NB.              argy=. vgety y
NB.   vgeto  - monad to extract output from ret;
NB.            is called as:
NB.              out=. vgeto ret
NB.   vrcond - dyad to find rcond; is called as:
NB.              rcond=. y vrcond out
NB.   vferr  - dyad to find ferr; is called as:
NB.              ferr=. y vferr out
NB.   vberr  - dyad to find berr; is called as:
NB.              berr=. y vberr out
NB.   mname  - literal, the name of monad vm to test
NB.   dname  - literal, the name of dyad vd to test
NB.   vtestm - monad to test monad vm and to log result:
NB.              mname rcond ferr berr time space
NB.            on the screen, in the global var TESTLOG and,
NB.            optionally, in the log file; is called as:
NB.              vtestm y
NB.   vtestd - monad to test dyad vd and to log result:
NB.              dname rcond ferr berr time space
NB.            on the screen, in the global var TESTLOG and,
NB.            optionally, in the log file; is called as:
NB.              vtestd y
NB.   vm     - monad to test; is called as:
NB.              ret=. vm argy
NB.   vd     - dyad to test; is called as:
NB.              ret=. argx vd argy
NB.   y      - some input for vtestm or vtestd
NB.   argx   - some left argument for vd
NB.   argy   - some right argument for vm or vd
NB.   ret    - some output from vm or vd
NB.   out    - rectified ret, i.e. filtered output
NB.   ferr   ≥ 0 or +∞ or indeterminate, the relative forward
NB.            error
NB.   berr   ≥ 0 or +∞ or indeterminate, the relative
NB.            backward error
NB.   rcond  ≥ 0, the estimated reciprocal of the condition
NB.            number of the input matrix; +∞ if matrix is
NB.            singular; indeterminate if matrix is
NB.            non-square
NB.
NB. Application 1:
NB.   NB. to estimate rcond in 1-norm
NB.   vrcond=. ((_."_)`(norm1 con getri) @. (=/@$))@[
NB.   NB. to calc. berr, assuming:
NB.   NB.   berr := ||A - realA||_1 / (m * ε * ||A||_1)
NB.   vberr=. ((- (% & norm1) [) % (FP_EPS * (norm1 * #) @ [)) unmqr
NB.   NB. let's test geqrf
NB.   ('geqrf' tmonad ]`]`vrcond`(_."_)`vberr) A
NB.
NB. Application 2:
NB.   NB. to estimate rcond in ∞-norm
NB.   vrcond=. ((_."_)`(normi con getri) @. (=/@$)) @ (0 {:: [)
NB.   NB. to calc. ferr, assuming:
NB.   NB.   ferr := ||x - realx||_inf / ||realx||_inf
NB.   vferr=. ((- (% & normi) [) (1&{::))~
NB.   NB. to calc. componentwise berr [LUG 75], assuming:
NB.   NB.   berr := max_i(|b - A * realx|_i / (|A| * |realx| + |b|)_i)
NB.   vberr=. ((mp & >/@[) (|@-) (0 {:: [) mp ]) (>./ @ %) (((0 {:: [) (mp & |) ]) + (|@mp & >/@[))
NB.   NB. let's test getrs
NB.   ('getrs' tdyad (0&{::)`(mp & >/)`]`vrcond`vferr`vberr) (A;x)

tmonad=: 2 : 0
  '`vgety vgeto vrcond vferr vberr'=. n
  argy=. vgety y
  try. 't s'=. timespacex 'ret=. ' , m , ' argy'      catch. t=. s=. ret=. _. end.
  try. out=. vgeto ret                                catch. out=. _.         end.
  try. rcond=. y vrcond out                           catch. rcond=. _        end.
  try. ferr=. y vferr out                             catch. ferr=. _.        end.
  try. berr=. y vberr out                             catch. berr=. _.        end.
  logline=. fmtlog m ; rcond ; ferr ; berr ; t ; s
  (logline , LF) ((1!:3) ^: (0 < (#@]))) TESTLOGFILE
  TESTLOG=: TESTLOG , logline
  logline (1!:2) 2
  EMPTY
)

tdyad=: 2 : 0
  '`vgetx vgety vgeto vrcond vferr vberr'=. n
  argx=. vgetx y
  argy=. vgety y
  try. 't s'=. timespacex 'ret=. argx ' , m , ' argy' catch. t=. s=. ret=. _. end.
  try. out=. vgeto ret                                catch. out=. _.         end.
  try. rcond=. y vrcond out                           catch. rcond=. _        end.
  try. ferr=. y vferr out                             catch. ferr=. _.        end.
  try. berr=. y vberr out                             catch. berr=. _.        end.
  logline=. fmtlog m ; rcond ; ferr ; berr ; t ; s
  (logline , LF) ((1!:3) ^: (0 < (#@]))) TESTLOGFILE
  TESTLOG=: TESTLOG , logline
  logline (1!:2) 2
  EMPTY
)

NB. ---------------------------------------------------------
NB. dbg
NB. Conj. to show verb's input and output
NB.
NB. Syntax:
NB.   vapp=. v dbg title
NB. where
NB.   title - any literal to name v
NB.   v     - verb to switch to debug mode
NB.   vapp  - being verb v augmented with output of incoming
NB.           parameters and outcoming result
NB.
NB. Application:
NB. - to debug verb '*' in verb (+/ .*) try:
NB.   C=. A (+/ .(* dbg '*')) B

dbg=: 2 : 0
  smoutput 'dbg' ; (n , ' [MONAD] ' , (": u b. 0)) ; 'y' ; (gshapes y) ; < y
  o=. u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (gshapes o) ; < o
  o
:
  smoutput 'dbg' ; 'x' ; (gshapes x) ; x ; (n , ' [DYAD] ' , (": u b. 0)) ; 'y' ; (gshapes y) ; < y
  o=. x u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (gshapes o) ; < o
  o
)

dbg2=: 2 : 0
  smoutput 'dbg' ; (n , ' [MONAD] ' , (": u b. 0)) ; 'y' ; (gshapes y)
  o=. u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (gshapes o)
  o
:
  smoutput 'dbg' ; 'x' ; (gshapes x) ; (n , ' [DYAD] ' , (": u b. 0)) ; 'y' ; (gshapes y)
  o=. x u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (gshapes o)
  o
)
