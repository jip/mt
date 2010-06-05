NB. gq.ijs
NB. Generate Q from LQ QL QR RQ output
NB.
NB. unglq  Generate a matrix with orthonormal rows from
NB.        output of gelqf
NB. ungql  Generate a matrix with orthonormal columns from
NB.        output of geqlf
NB. ungqr  Generate a matrix with orthonormal columns from
NB.        output of geqrf
NB. ungrq  Generate a matrix with orthonormal rows from
NB.        output of gerqf
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. Block size limits for blocked code
UNGQBSMIN=: 2
UNGQBSMAX=: 32

NB. ---------------------------------------------------------
NB. Differences between rIOSs at consequent iterations:
NB.   rios(i+1)-rios(i)

NB. for non-blocked code
UNGL2DCIOS=: 3 2 2 $ 0 0 0 0 _1 _1 1 1 _1 0 0 _1  NB. I,eC,rto0
UNG2LDCIOS=: 3 2 2 $ 0 0 0 0 0 0 1 1 1 1 _1 0     NB. I,eC,cto0
UNG2RDCIOS=: 3 2 2 $ 0 0 0 0 _1 _1 1 1 0 _1 _1 0  NB. I,eC,cto0
UNGR2DCIOS=: 3 2 $ 1 0 0 0 0 0 0 1 1 1 1 0 _1     NB. I,eC,rto0

NB. ---------------------------------------------------------
NB. mkrios0ungl2
NB. mkrios0geql2
NB. mkrios0geqr2
NB. mkrios0gerq2
NB.
NB. Create cIOS at 0-th iteration for corresp. method
NB.
NB. Syntax:
NB.   rios0=. mkrios0ungl2 (p,m,n1)
NB.   rios0=. mkrios0geql2 (m,n,k,p)
NB.   rios0=. mkrios0geqr2 (m,n,k,p)
NB.   rios0=. mkrios0gerq2 (m,n,k,p)
NB. where
NB.   kn    - 2-vector of integers (k,n), shape of matrix to
NB.           generate
NB.   rios0 - 6×2-table rios(0), cIOSs corresponding to
NB.           iteration 0, see ung*step verbs

mkrios0ungl2=: 3 : 0
  'p1 k'=. 0 1 1 (<:`(<./))/. 'p m n'=. 0 0 _1 + y
  3 2 2 $ p,0,(k-p),n,p1,p1,(k-p1),(n+1-p1),p1,0 1,p1
)

mkrios0ung2l=: 3 : 0
  'p1 p2'=. (- & 1 2 @ {:) 'm n k p'=. y
  3 2 2 $ 1,(n-k),m,(k-p),0,(n-k),(m-p2),(k-p1),(m-p2),(n-p),p1,1
)

mkrios0ung2r=: 3 : 0
  p1=. <: {: 'm n k p'=. y
  3 2 2 $ 0,p,m,(k-p),p1,p1,(m+1-p1),(k-p1),0,p1,p1,1
)

mkrios0ungr2=: 3 : 0
  'p1 p2'=. (- & 1 2 @ {:) 'm n k p'=. y
  3 2 2 $ (m-k),1,(k-p),n,(m-k),0,(k-p1),(n-p2),(m-p),(n-p2),1,p1
)

NB. ---------------------------------------------------------
NB. Template adv. to form verbs ung((2[lr])|([lr]2))step
NB. algo:
NB. LQ2 eC == (β,v,τ) , (C ,. trash)
NB.     r2z before z

ungq2step=: 2 : '([ (0: setir 2) (((((- @ (m gi 1)) ((>:@[) (m gi 2) *) ]) (n}) ([ - (m gi 1) * (m gi 0))) (n&{)) updir 1)) step'

ungl2step=: (larfr (+ @ (1 0 & (0 _1}))))`(+ @ (_1 { ]))`( 0}) ungq2step IOSFR
ung2lstep=: (larfl (    (0 1 & (0 _1}))))`(    ( 0 { ]))`(_1}) ungq2step IOSLC
ung2rstep=: (larfl (    (1 0 & (0 _1}))))`(    (_1 { ]))`( 0}) ungq2step IOSFC
ungr2step=: (larfr (+ @ (0 1 & (0 _1}))))`(+ @ ( 0 { ]))`(_1}) ungq2step IOSLR

NB. ---------------------------------------------------------
NB. ungl2
NB. Generate a matrix with orthonormal rows from output of
NB. gelq2 or gelqf (non-blocked version)
NB.
NB. Syntax:
NB.   eQ=. [p] ungl2 LQf
NB. where
NB.   LQf - m×(n+1)-matrix, output of gelqf
NB.   p   - integer in range 0:k, default is k, the number of
NB.         elementary reflectors, whose product defines the
NB.         matrix Q
NB.   eQ  - m×(n+1)-matrix containing Q:
NB.           Q -: (k , n) {. eQ
NB.   Q   - k×n-matrix with orthonormal rows, which is
NB.         defined as the first k rows of a product of p
NB.         elementary reflectors of order n:
NB.           Q = H(p)' ... H(2)' H(1)'
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'm n'=. $ A
NB.   k=. m <. n
NB.   LQf=. gelq2 A
NB.   L=. trl 0 _1 }. LQf
NB.   Q=. (k,n) {. ungl2 LQf
NB.   Q2=. ungl2 tru1 LQf
NB. then
NB.   Q -: Q2
NB.   I -: (mp ct) Q
NB.   A -: L mp Q
NB.   (-: (((trl @ (0 _1 & }.)) mp (ungl2 @ tru1)) @ gelq2)) A
NB.
NB. Algo:
NB.   p-th diagonal := 1

ungl2=: ($:~ (<./ @ (0 _1 & +) @ $)) : '[ ((0&{::) @ (UNGL2DCIOS & ungl2step)) ((mkrios0ung2l @ (, $)) ((((0 ({,) [) idmat ((<0 1) { [)) setir 0) ; [) ])'

NB. ---------------------------------------------------------
NB. ung2l
NB. Generate a matrix with orthonormal columns from output of
NB. geql2 or geqlf (non-blocked version)
NB.
NB. Syntax:
NB.   Q=. [p] ung2l QfL
NB. where
NB.   QfL - (m+1)×n-matrix, output of geql2 or geqlf
NB.   p   - integer in range 0:k, default is k, the number of
NB.         elementary reflectors, whose product defines the
NB.         matrix Q
NB.   Q   - m×k-matrix Q with orthonormal columns, which is
NB.         defined as the last k columns of a product of p
NB.         elementary reflectors of order m:
NB.           Q = H(p) ... H(2) H(1)
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'm n'=. $ A
NB.   QfL=. geql2 A
NB.   Q=. ung2l QfL
NB.   L=. (n - m) trl }. QfL
NB. then
NB.   I -: (mp~ ct) Q
NB.   A -: Q mp L
NB.   (-: ((ung2l mp (((n - m) & trl) @ }.)) @ geql2)) A

ung2l=: ($:~ (<./ @ (_1 0 & +) @ $)) :(4 : 0)
  k=. <./ 'm n'=. _1 0 + $ y
  rios0=. mkrios0ung2l (m , n , k , x)
  y=. rios0 (((7 6 (>:@-/@({,)) [) idmat ((<0 1) { [)) setir 0) y  NB. 1+(k+1-p)-(m+2-p)=(k-m)-th diagonal := 1
  (- (m , k)) {. 0 {:: x (UNG2LDCIOS & ung2lstep) (y ; rios0)
)

NB. ---------------------------------------------------------
NB. ung2r
NB. Generate a matrix with orthonormal columns from output of
NB. geqr2 or geqrf (non-blocked version)
NB.
NB. Syntax:
NB.   Q=. [p] ung2r QfR
NB. where
NB.   QfR - (m+1)×n-matrix, output of geqr2 or geqrf
NB.   p   - integer in range 0:k, default is k, the number of
NB.         elementary reflectors, whose product defines the
NB.         matrix Q
NB.   Q   - m×k-matrix Q with orthonormal columns, which is
NB.         defined as the first k columns of a product of p
NB.         elementary reflectors of order m:
NB.           Q = H(1) H(2) ... H(p)
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'm n'=. $ A
NB.   QfR=. geqr2 A
NB.   Q=. ung2r QfR
NB.   R=. tru }: QfR
NB. then
NB.   I -: (mp~ ct) Q
NB.   A -: Q mp R
NB.   (-: ((ung2r mp (tru @ }:)) @ geqr2)) A

ung2r=: ($:~ (<./ @ (_1 0 & +) @ $)) :(4 : 0)
  k=. <./ 'm n'=. _1 0 + $ y
  rios0=. mkrios0ung2к (m , n , k , x)
  y=. rios0 (((1 ({,) [) idmat ((<0 1) { [)) setir 0) y  NB. p-th diagonal := 1
  (m , k) {. 0 {:: x (UNG2RDCIOS & ung2rstep) (y ; rios0)
)

NB. ---------------------------------------------------------
NB. ungr2
NB. Generate a matrix with orthonormal rows from output of
NB. gerq2 or gerqf (non-blocked version)
NB.
NB. Syntax:
NB.   Q=. [p] ungr2 RQf
NB. where
NB.   RQf - m×(n+1)-matrix, output of gerq2 or gerqf
NB.   p   - integer in range 0:k, default is k, the number of
NB.         elementary reflectors, whose product defines the
NB.         matrix Q
NB.   Q   - k×n-matrix Q with orthonormal rows, which is
NB.         defined as the last k rows of a product of p
NB.         elementary reflectors of order n:
NB.           Q = H(1)' H(2)' ... H(p)'
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'm n'=. $ A
NB.   RQf=. gerq2 A
NB.   R=. (n - m) trl 0 1 }. RQf
NB.   Q=. ungr2 RQf
NB. then
NB.   I -: (mp ct) Q
NB.   A -: R mp Q
NB.   (-: (((((n - m) & tru) @ (0 1 & }.)) mp ungr2) @ gerq2)) A

ungr2=: ($:~ (<./ @ (0 _1 & +) @ $)) :(4 : 0)
  k=. <./ 'm n'=. 0 _1 + $ y
  rios0=. mkrios0ungr2 (m , n , k , x)
  y=. rios0 (((7 2 (-&2@-/@({,)) [) idmat ((<0 1) { [)) setir 0) y  NB. ((n+2-p)-(k-p))-2=(n-k)-th diagonal := 1
  (- (k , n)) {. 0 {:: x (UNGR2DCIOS & ungr2step) (y ; rios0)
)

NB. =========================================================
NB. Interface

NB. for blocked code (see ung{lq,ql,qr,rq}step)
NB. in 'block size' units
UNGLQDCIOS=: 4 2 2 $ 0 0 0 0 0 0 0 0 _1 _1 0 1 _1 _1 1 1  NB. LQ0,T,LQi,C

NB. rios0=. mkrios0unglq (t,l,m,n,bs,iters)
mkrios0unglq=: 3 : 0
  't l m n bs iters'=. y
  ibs=. bs*iter
  4 2 2 $ ibs,ibs,(m-ibs),(n+1-ibs),0,(n+1),bs,bs,(ibs-bs),(ibs-bs),bs,(n+1+bs-ibs),ibs,(ibs-bs),(m-ibs),(n+1+bs-ibs)
)

unglqstep=: [:

NB. 'bs iters'=. sha ungenv riosQf
NB. where shc - 2-vector to convert Qf shape to Q shape
NB.       riosQf -: ((topQf,leftQf),:(heightQf,widthQf))
NB. Note: ($ Q) -: (shc + (heightQf,widthQf))
NB.
NB. Formula:
NB.   k := min(sha + (heightQf,widthQf))     NB. min(heightQ,widthQ)
NB.   bs := max(UNGQBSMIN,min(UNGQBSMAX,k))  NB. adjusted block size
NB.   iters := ⌊k % bs⌋                      NB. # of iterations

ungenv=: ((] , (<:@%)) (UNGQBSMIN >. UNGQBSMAX <. ])) @ (<./ @ (+ {:))

NB. Q=.           unglq Qf
NB. Aupd=. riosQf unglq A

unglq=: ((0 _1 + $) {. ((0 0 ,: $) $: ])) : (4 : 0)
  'bs iters'=. 0 _1 ungenv x
  drios=. bs * UNGLQDCIOS           NB. rIOS increment between iterations
  iters (0 {:: (drios & unglqstep)) (x ((mkrios0unglq $) (([ (0: setir 1) ((0 { [) ungl2 ])) ; [) ]) y)
)



iters (((0,(-bs)) & }.) @ ((1 {:: ]) (gelq2 updir 5) (0 {:: ])) @ (drios & unglqstep)) (y ; rios0)
  
  k=. <./ 'm n'=. 0 _1 + $ y
  bs=. UNGQBSMIN >. UNGQBSMAX <. k  NB. adjusted block size
  y=. (m , (n+1+bs)) {. y           NB. allocate space for T, filled by zeros
  (0 0 ,: (m,n)) unglq y
:
  k=. <./ _2 {. 't l m n'=. , x     NB. top,left,height,width of Qf
  bs=. UNGQBSMIN >. UNGQBSMAX <. k  NB. adjusted block size
  iters=. <. k % bs                 NB. # of iterations
  rios0=. mkrios0unglq (,x) , bs , iters
  drios=. bs * UNGLQDCIOS
   (y ; rios0)
  iters (((0,(-bs)) & }.) @ ((1 {:: ]) (gelq2 updir 5) (0 {:: ])) @ (drios & unglqstep)) (y ; rios0)
)

--------
gelqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  bs=. GEQFBSMIN >. GEQFBSMAX <. k  NB. adjusted block size
  y=. (m , (n+1+bs)) {. y           NB. allocate space for τ and T, filled by zeros
  (m,n,bs,iters) gelqf y
:
  'bs iters'=. _2 {. x
  rios0=. mkrios0gelqf x
  drios=. bs * GELQFDRIOS
  iters (((0,(-bs)) & }.) @ ((1 {:: ]) (gelq2 updir 5) (0 {:: ])) @ (drios & gelqfstep)) (y ; rios0)
)
--------

NB. stubs
ungql=: ung2l
ungqr=: ung2r
ungrq=: ungr2

NB. TODO:
NB. - template adv. ung2

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB.*tgq v test ung*

tgq=: 3 : 0
)

NB. ---------------------------------------------------------
NB.*testgq v test Q genberators

testgq=: 3 : 0
  Ac6x6=. 6 6 $ 8j_7 _8j_8 1j_6 _4j1 4j_2 _7j_2 _4j9 _8j_2 _1j_4 1j_7 5j_8 6j_9 5j_9 1j_5 _9j_7 _6j7 _1j_5 1j8 5j5 5j6 _3j_7 2j2 _8j_7 9j8 _2j6 8j_1 _6j7 _5j_3 _5j_6 _8 1j7 _7j8 6j2 3 8j_5 _1
  Ac6x4=. 6 4 $ _3j1 _5j_2 1 _9j6 8j_4 _2j2 _8j_1 _6j6 _2j5 _8j5 _5j_1 3j3 _7j_6 _7j_8 _3j_4 4j_2 5j4 8j8 6 3j6 _2j_6 9j5 _9 _1j1
  Ac4x6=. 4 6 $ _3j_5 2j2 _1j_4 _5j_8 _4j8 8j_5 8j_2 2j4 _5j_3 _6j_5 8j4 9j_1 7j2 _6j_5 7j_1 _5j6 _1 _1j9 1j_5 _4j_7 _6 _1j_9 _6j1 _7j1

  smoutput 'LQ 6x6' ; (; (((trl @ (0 _1 & }.)) mp ungl2) @ gelq2)) Ac6x6
  smoutput 'LQ 6x4' ; (; (((trl @ (0 _1 & }.)) mp ungl2) @ gelq2)) Ac6x4
  smoutput 'LQ 4x6' ; (; (((trl @ (0 _1 & }.)) mp ungl2) @ gelq2)) Ac4x6

  smoutput 'QL 6x6' ; (; ((ung2lx mp (((6 - 6) & trl) @ }.)) @ geql2)) Ac6x6
  smoutput 'QL 6x4' ; (; ((ung2lx mp (((4 - 6) & trl) @ }.)) @ geql2)) Ac6x4
  smoutput 'QL 4x6' ; (; ((ung2lx mp (((6 - 4) & trl) @ }.)) @ geql2)) Ac4x6

  smoutput 'QR 6x6' ; (; ((ung2rx mp (tru @ }:)) @ geqr2)) Ac6x6
  smoutput 'QR 6x4' ; (; ((ung2rx mp (tru @ }:)) @ geqr2)) Ac6x4
  smoutput 'QR 4x6' ; (; ((ung2rx mp (tru @ }:)) @ geqr2)) Ac4x6

  smoutput 'RQ 6x6' ; (; (((((6 - 6) & tru) @ (0 1 & }.)) mp ungr2) @ gerq2)) Ac6x6
  smoutput 'RQ 6x4' ; (; (((((4 - 6) & tru) @ (0 1 & }.)) mp ungr2) @ gerq2)) Ac6x4
  smoutput 'RQ 4x6' ; (; (((((6 - 4) & tru) @ (0 1 & }.)) mp ungr2) @ gerq2)) Ac4x6
)
