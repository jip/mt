NB. util.ijs
NB. Utilities
NB.
NB. ii2cp      Make cycle permutation from indices x and y
NB.
NB. sgn        Simplified signum
NB. condneg    Conditional negate
NB. copysign   Copy sign
NB.
NB. hds2ios    IOS from head, delta and size
NB. ht2ios     IOS from head and tail
NB. hs2ios     IOS from head and size
NB. rios2ios   Convert rIOS to IOS
NB. rios2lios  Convert rIOS to lIOS
NB.
NB. norm1      Magnitude-based 1-norm of vector or matrix
NB. normi      Magnitude-based ∞-norm of vector or matrix
NB. norm1t     Taxicab-based 1-norm of vector or matrix
NB. normit     Taxicab-based ∞-norm of vector or matrix
NB. norms      Square-based (Euclidean/Frobenius) norm of
NB.            vector or matrix
NB.
NB. fmtlog     Format log string
NB.
NB. tmonad     Template conj. to make verbs to test
NB.            computational monad
NB. tdyad      Template conj. to make verbs to test
NB.            computational dyad
NB. dbg        Conj. to sjow verb's input and output
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

rios2oios=: < @ hs2ios/ " 1 @: |:  NB. convert rIOS to opened (non-boxed) IOS
gshapes=: $`($(;<)($ L: 0)) @. (0 < L.)   NB. get shapes, boxes are accepted, too

NB. template adverbs to form norm verbs
mocs=: >./ @ (+/     ) @:          NB. vector: sum of, matrix: max of column sums
mors=: >./ @ (+/ " _1) @:          NB. vector: max of, matrix: max of row sums

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Miscellaneous

ii2cp=: < @ (, ` (, @ ]) @. =)     NB. make cycle permutation from indices x and y (NB!: if y is empty then output must be a:)
                                   NB. CHECKME: ii2cp=: < @ , @ (, ^: ~:)
                                   NB. CHECKME: ii2cp=: (, ^: (-. @ -:)) & ,
                                   NB. sw=: <@~.@[ C. ]
                                   NB.   [Jprogramming] Swapping array elements
                                   NB.   Roger Hui, Mon May 11 06:23:07 HKT 2009
                                   NB.   http://www.jsoftware.com/pipermail/programming/2009-May/014682.html

sgn=: 0 & (<: - >)                 NB. if y<0 then -1 else 1 endif
condneg=: (* sgn)~                 NB. if x<0 then -y else y endif
copysign=: condneg |               NB. if x<0 then -|y| else |y| endif

fmtlog=: '%-25S %-12g %-12g %-12g %-12g %12d' vsprintf  NB. Format log string

NB. ---------------------------------------------------------
NB. Norms

NB. Magnitude-based norms |y|
norm1=: | mocs                     NB. 1-norm of vector or matrix
normi=: | mors                     NB. ∞-norm of vector or matrix

NB. Taxicab-based norms |Re(y)| + |Im(y)|
norm1t=: (+/ " 1 @: | @: +.) mocs  NB. 1-norm of vector or matrix
normit=: (+/ " 1 @: | @: +.) mors  NB. ∞-norm of vector or matrix

NB. Square-based (Euclidean/Frobenius) norm of vector or matrix
NB. for vector input emulates LAPACK's DZNRM2 and SCNRM2
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
hs2ios=: [ + ((] * i. @ *) sgn)~

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
NB.   vberr=. ((- (% & norm1) [) % (FP_EPS * # @ [)) unmqr
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
