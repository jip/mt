NB. util.ijs
NB. Utilities
NB.
NB. ii2cp     Make cycle permutation from indices x and y
NB.
NB. sgn       Simplified signum
NB. condneg   Conditional negate
NB. copysign  Copy sign
NB.
NB. hds2ios   IOS from head, delta and size
NB. ht2ios    IOS from head and tail
NB. hs2ios    IOS from head and size
NB. rios2ios  Convert rIOS to IOS
NB.
NB. step      Template adv. to make verbs of single iteration
NB.
NB. prn       Formatted console output
NB. dbg       Conj. to show verb's input and output
NB.
NB. norm1     Magnitude-based 1-norm of vector or matrix
NB. normi     Magnitude-based ∞-norm of vector or matrix
NB. norm1t    Taxicab-based 1-norm of vector or matrix
NB. normit    Taxicab-based ∞-norm of vector or matrix
NB. norms     Square-based (Euclidean/Frobenius) norm of vector or matrix
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. template adverbs to form norm verbs
mocs=: >./ @ (+/     ) @:          NB. vector: sum of, matrix: max of column sums
mors=: >./ @ (+/ " _1) @:          NB. vector: max of, matrix: max of row sums

NB. =========================================================
NB. Interface

ii2cp=: < @ (, ` (, @ ]) @. =)     NB. make cycle permutation from indices x and y (NB!: if y is empty then output must be a:)
                                   NB. CHECKME: ii2cp=: < @ , @ (, ^: ~:)
                                   NB. CHECKME: ii2cp=: (, ^: (-. @ -:)) & ,

sgn=: 0 & (<: - >)                 NB. if y<0 then -1 else 1 endif
condneg=: (* sgn)~                 NB. if x<0 then -y else y endif
copysign=: condneg |               NB. if x<0 then -|y| else |y| endif

NB. ---------------------------------------------------------
NB. Norms

NB. Magnitude-based norms |y|
norm1=: | mocs                     NB. 1-norm of vector or matrix
normi=: | mors                     NB. ∞-norm of vector or matrix

NB. Taxicab-based norms |Re(y)| + |Im(y)|
norm1t=: (+/ " 1 @: | @: +.) mocs  NB. 1-norm of vector or matrix
normit=: (+/ " 1 @: | @: +.) mors  NB. ∞-norm of vector or matrix

NB. Square-based (Euclidean/Frobenius) norm of vector or matrix
NB. for vector input emulates LAPACK's DZNRM2
norms=: (((((+/^:_) &.: *:) @: %) * ]) (>./^:_)) @: | @: +.

NB. ---------------------------------------------------------
NB. IOS generators

NB. y[2]-vector of integers from head y[0] by delta y[1]
hds2ios=: + ` (* i.)/

NB. (x-y)-vector of integers from head y to tail (x-1):
NB. y (y+1) ... (x-1)
ht2ios=: ] + (i. @ -)

NB. y-vector of integers from head x; models rIOS in (u;.0)
hs2ios=: [ + ((] * i. @ *) sgn)~

NB. ---------------------------------------------------------
NB. IOS converters

NB. Convert rIOS to IOS
NB. Note: IOS with length less than array's rank indexes
NB.       the slice
rios2ios=: < " 1 @ (< @ hs2ios/ " 1 @: |:)

NB. ---------------------------------------------------------
NB. step
NB. Template adv. to make verbs of single iteration
NB.
NB. Syntax:
NB.   vstep=: vchange step
NB. where
NB.   vchange - verb to change matrix, is called as:
NB.               Ai1=. riosi vchange Ai
NB.   vstep   - verb to do single iteration, is called as:
NB.               'Ai1 riosi1'=. drios vstep (Ai ; riosi)
NB.   Ai      - matrix A(i) to update before i-th
NB.             iteration
NB.   riosi   - matrix rios(i) of rIOS (rIOSs) for i-th
NB.             iteration
NB.   drios   - difference between rIOS (rIOSs) at
NB.             consequent iterations: rios(i+1)-rios(i)
NB.   Ai1     - matrix A(i+1) after i-th iteration
NB.   riosi1  - matrix rios(i+1) of rIOSs for (i+1)-th
NB.             iteration
NB.   i       ≥ 0
NB.
NB. Algorithm:
NB.   1) Ai1=. riosi vchange Ai
NB.   2) Tmp=. ((< Ai1) 0} Input
NB.   3) riosi1=. drios + riosi
NB.   4) Output=. (< riosi1) 1} Tmp
NB.   where
NB.     Input -: (Ai ; riosi)
NB.     Tmp -: (Ai1 ; riosi)
NB.     Output -: (Ai1 ; riosi1)

step=: 1 : '(< @ (+ (1&{::))) 1} (((1 {:: ]) (< @ u) (0 {:: ])) 0} ])'

NB. ---------------------------------------------------------
NB. prn
NB. Formatted console output

fmtlog=: '%-25S %-12g %-12g %-12g %-12g %12d' & vsprintf

NB. ---------------------------------------------------------
NB. tmonad1pass
NB. tmonad2pass
NB. tdyad1pass
NB. tdyad2pass
NB.
NB. Test single computational verb
NB.
NB. Syntax:
NB.   vapp=. v2test`vi2test`vferr`vberr`vnorm test1
NB. where
NB.   v2test  - monad to test; is called as:
NB.               out=. v2test in
NB.   vi2test - v2test's inversion; is called as:
NB.               in=. vi2test out
NB.   vferr   - dyad to find ferr of some realvalue
NB.             relatively to some modelvalue; is called as:
NB.               ferr=. modelvalue vferr realvalue
NB.   vberr   - dyad to find berr of some realvalue
NB.             relatively to some modelvalue; is called as:
NB.               berr=. modelvalue vberr realvalue
NB.   vnorm   - monad to find norm; is called as:
NB.               norm=. vnorm array
NB.   vapp    - dyad to try to estimate reciprocal of
NB.             condition number of the input matrix (only
NB.             for square matrices), ferr, berr, execution
NB.             time and space of verb v2test, and optionally
NB.             save result into log file and/or log array
NB.             and/or console; is called as:
NB.               vname vapp in
NB.   in      - some input for v2test and output for vi2test
NB.   out     - some output for v2test and input for vi2test
NB.   ferr    ≥ 0, relative forward error
NB.   berr    ≥ 0, relative backward error
NB.   norm    ≥ 0, matrix or vector norm
NB.   vname   - literal array, the name of v2test
NB.
NB. Application:
NB.   'geqrf' (geqrf`unmqr`vferr`vberr`vnorm test1) A
NB.
NB. TODO:
NB.   try. catch.

tmonad2pass=: 1 : 0
:
  '`v2test vi2test vferr vberr vnorm'=. m
  try.
    modelout=. v2test y
    modeliny=. vi2test modelout
    rcond=. (_."_)`(vnorm con getri) @. (=/@$) modeliny
    't s'=. timespacex 'realout=. ' , x , ' modeliny'
    realiny=. vi2test realout
    try. ferr=. modeliny vferr realiny catch. ferr=. _. end.
    try. berr=. modelout vberr realout catch. berr=. _. end.
  catch.
    rcond=. berr=. ferr=. t=. s=. _.
  end.
  logline=. fmtlog x ; rcond ; berr ; ferr ; t ; s
  logline ((1!:3) ^: (0 < (#@]))) TESTLOGFILE
  TESTLOG=: TESTLOG , logline
  logline (1!:2) 2
  EMPTY
)

tmonad1pass=: 1 : 0
:
  '`v2test vi2test vferr vberr vnorm'=. m
  try.
    modeliny=. y
    rcond=. (_."_)`(vnorm con getri) @. (=/@$) modeliny
    't s'=. timespacex 'realout=. ' , x , ' modeliny'
    realiny=. vi2test realout
    ferr=. _.
    try. berr=. modelout vberr realout catch. berr=. _. end.
  catch.
    rcond=. berr=. ferr=. t=. s=. _.
  end.
  logline=. fmtlog x ; rcond ; berr ; ferr ; t ; s
  logline ((1!:3) ^: (0 < (#@]))) TESTLOGFILE
  TESTLOG=: TESTLOG , logline
  logline (1!:2) 2
  EMPTY
)

tdyad2pass=: [:

tdyad1pass=: [:

NB. ---------------------------------------------------------
NB. dbg
NB. Conj. to show verb's input and output
NB.
NB. Application:
NB. - to debug verb '*' in verb (+/ .*) try:
NB.   C=. A (+/ .(* dbg '*')) B

dbg=: 2 : 0
  smoutput 'dbg' ; (n , ' [MONAD] ' , (": u b. 0)) ; 'y' ; (($`($ L: 0) @. (0 < L.)) y) ; < y
  o=. u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (($`($ L: 0) @. (0 < L.)) o) ; < o
  o
:
  smoutput 'dbg' ; 'x' ; (($`($ L: 0) @. (0 < L.)) x) ; x ; (n , ' [DYAD] ' , (": u b. 0)) ; 'y' ; (($`($ L: 0) @. (0 < L.)) y) ; < y
  o=. x u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (($`($ L: 0) @. (0 < L.)) o) ; < o
  o
)
