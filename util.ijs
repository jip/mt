NB. util.ijs
NB. Utilities
NB.
NB. sgn        Simplified signum
NB. condneg    Conditional negate
NB. copysign   Copy sign
NB. fmtlog     Format log string
NB. gi         Conj. to evoke n-th verb from gerund m
NB. ms         Minimum in sum of vectors
NB.
NB. norm1      Magnitude-based 1-norm of vector or matrix
NB. normi      Magnitude-based ∞-norm of vector or matrix
NB. norm1t     Taxicab-based 1-norm of vector or matrix
NB. normit     Taxicab-based ∞-norm of vector or matrix
NB. norms      Square-based (Euclidean/Frobenius) norm of
NB.            vector or matrix
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

gshapes=: $`($(;<)($ L: 0)) @. (0 < L.)  NB. get shapes, boxes are accepted, too
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
