NB. qf.ijs
NB. Orthogonal factorizations LQ QL QR RQ
NB.
NB. gelqf  LQ factorization of a general matrix
NB. geqlf  QL factorization of a general matrix
NB. geqrf  QR factorization of a general matrix
NB. gerqf  RQ factorization of a general matrix
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. Block size limits for blocked code
GEQFBSMIN=: 2
GEQFBSMAX=: 32

NB. ---------------------------------------------------------
NB. Differences between rIOSs at consequent iterations:
NB.   rios(i+1)-rios(i)

NB. for non-blocked code, reduces height and width by 1
GEQ2DRIOS=: 2 2 $ 0 0 _1 _1

NB. for blocked code (see ge{lq,ql,qr,rq}fstep)
NB. in 'block size' units
GELQFDRIOS=: 6 2 2 $ 1 1 0 _1 1 1 0 _1 1 0 0 0 0 0 0 0 1 1 _1 _1 0 0 0 0
GEQLFDRIOS=: 6 2 2 $ 0 _1 _1 0 0 _1 _1 0 0 _1 0 0 0 0 0 0 0 0 _1 _1 0 0 0 0
GEQRFDRIOS=: 6 2 2 $ 1 1 _1 0 1 1 _1 0 0 1 0 0 0 0 0 0 1 1 _1 _1 0 0 0 0
GERQFDRIOS=: 6 2 2 $ _1 0 0 _1 _1 0 0 _1 _1 0 0 0 0 0 0 0 0 0 _1 _1 0 0 0 0

NB. ---------------------------------------------------------
NB. geq2step
NB. Template conj. to make dyads ge*2step
NB.
NB. Syntax:
NB.   vstep=. vref`vmul`vtau geq2step iosyz
NB. where
NB.   vref  - monad to generate an elementary reflector, see
NB.           larfg* larfp*; is called as:
NB.             z=. vref y
NB.   vmul  - dyad to extract v from z and to product
NB.           (v*v'*Ci) or (Ci*v*v'); see larfl, larfr; is
NB.           called as:
NB.             aCi1=. aCi vmul z
NB.   vtau  - monad to extract τ from z and pre-process it;
NB.           is called as:
NB.             tau=. vtau z
NB.   iosyz - IOS vectors y and z in matrix aCi
NB.   vstep - verb to perform single step of Q-factorization;
NB.           see ge*2step; is called as:
NB.             'aAi1 riosi1'=. drios vstep (aAi ; riosi)
NB.
NB. Algorithm:
NB.   1) y=. iosyz { aCi
NB.   2) z=. vref y
NB.   3) tmp=. aCi - (vtau z) * (aCi vmul z)
NB.   4) aCi1=. z iosyz } tmp
NB. where
NB.   y    - k-vector or k×1-matrix to reflect:
NB.            (α,x,0) or (0,x,α)
NB.   z    - k-vector or k×1-matrix (the same shape as of
NB.          y), the result of y reflection:
NB.            (β,v,τ) or (τ,v,β)
NB.   aCi  - matrix, being matrix Ci, augmented by vector y
NB.          and zero vector, before i-th iteration:
NB.            for LQ: (m-i-1)×(n-i)-matrix (y,(Ci,.0)) == A[i:i+m-1,i:n]
NB.            for QL: (m-i)×(n-i-1)-matrix ((0,Ci),.y) == A[0:m-i,0:n-i-1]
NB.            for QR: (m-i)×(n-i-1)-matrix (y,.(Ci,0)) == A[i:i+m,i:n-1]
NB.            for RQ: (m-i-1)×(n-i)-matrix ((0,.Ci),y) == A[0:m-i-1,0:n-i]
NB.   Ci   - matrix C(i) to update
NB.   aCi1 - matrix, being matrix Ci1, augmented by vector
NB.          z and zero vector, after i-th iteration:
NB.            for LQ: (m-i-1)×(n-i)-matrix (z,(Ci1,.0)) == A[i:i+m-1,i:n]
NB.            for QL: (m-i)×(n-i-1)-matrix ((0,Ci1),.z) == A[0:m-i,0:n-i-1]
NB.            for QR: (m-i)×(n-i-1)-matrix (z,.(Ci1,0)) == A[i:i+m,i:n-1]
NB.            for RQ: (m-i-1)×(n-i)-matrix ((0,.Ci1),z) == A[0:m-i-1,0:n-i]
NB.   Ci1  - matrix C(i+1), being updated C(i) by matrix H:
NB.            for LQ and RQ: Ci1 := Ci*H
NB.            for QL and QR: Ci1 := H'*Ci
NB.   H    - matrix, an elementary reflector formed as:
NB.            H := I - τ * (1,v) * (1,v)'

geq2step=: 2 : '(((] (n}) ([ - (m gi 2) * (m gi 1))) ((m gi 0) @ (n&{))) updr) step'

NB. ---------------------------------------------------------
NB. gelq2step
NB. geql2step
NB. geqr2step
NB. gerq2step
NB. Single step of ge*2
NB.
NB. Syntax:
NB.   'aAi1 riosi1'=. drios gelq2step (aAi ; riosi)
NB.   'aAi1 riosi1'=. drios geql2step (aAi ; riosi)
NB.   'aAi1 riosi1'=. drios geqr2step (aAi ; riosi)
NB.   'aAi1 riosi1'=. drios gerq2step (aAi ; riosi)
NB. where
NB.   aAi    - matrix, being matrix Ai augmented by vector
NB.            τi, before i-th iteration:
NB.              for LQ: m×(n+1)-matrix (Ai ,. τi)
NB.              for QL: (m+1)×n-matrix (τi , Ai)
NB.              for QR: (m+1)×n-matrix (Ai , τi)
NB.              for RQ: m×(n+1)-matrix (τi ,. Ai)
NB.   Ai     - m×n-matrix A(i) containing matrix aCi
NB.   τi     - vector, being vector τ[0:k-1] filled by zeros
NB.            upto proper length, before i-th iteration:
NB.              for LQ: m-vector (τ[0:i-2],0[i:m-1])
NB.              for QL: n-vector (0[0:n-i],τ[k-i+1:k-1])
NB.              for QR: n-vector (τ[0:i-2],0[i:n-1])
NB.              for RQ: m-vector (0[0:m-i],τ[k-i+1:k-1])
NB.   riosi  - 2×2-matrix, Ci's rIOS rios(i) for i-th
NB.            iteration
NB.   aAi1   - matrix, being matrix Ai1 augmented by vector
NB.            τi1, after i-th iteration:
NB.              for LQ: m×(n+1)-matrix (Ai1 ,. τi1)
NB.              for QL: (m+1)×n-matrix (τi1 , Ai1)
NB.              for QR: (m+1)×n-matrix (Ai1 , τi1)
NB.              for RQ: m×(n+1)-matrix (τi1 ,. Ai1)
NB.   Ai1    - m×n-matrix A(i+1) containing matrix aCi1
NB.   τi1    - vector, being vector τi after i-th iteration:
NB.              for LQ: m-vector (τ[0:i-1],0[i+1:m-1])
NB.              for QL: n-vector (0[0:n-i-1],τ[k-i:k-1])
NB.              for QR: n-vector (τ[0:i-1],0[i+1:n-1])
NB.              for RQ: m-vector (0[0:m-i-1],τ[k-i:k-1])
NB.   drios  - difference between rIOSs at consequent
NB.            iterations: rios(i+1)-rios(i)
NB.   riosi1 - 2×2-matrix, Ci1's rIOS rios(i+1) for (i+1)-th
NB.            iteration
NB.   i      - integer in range 0:k-1, iteration number
NB.   k      = min(m,n)

gelq2step=: larfgfc`(larfr (+ @ (1 0 & (0 _1}))))`     (_1 { ])  geq2step IOSFR
geql2step=: larfgb `(larfl      (0 1 & (0 _1})) )`(+ @ ( 0 { ])) geq2step IOSLC
geqr2step=: larfgf `(larfl      (1 0 & (0 _1})) )`(+ @ (_1 { ])) geq2step IOSFC
gerq2step=: larfgbc`(larfr (+ @ (0 1 & (0 _1}))))`     ( 0 { ])  geq2step IOSLR

NB. ---------------------------------------------------------
NB. geq2
NB. Template adv. to make monads ge*2
NB.
NB. Syntax:
NB.   vapp=. vrsh`vrios0`vstep geq2
NB. where
NB.   vrsh   - monad to restore A's shape; is called as:
NB.              sh=. vrsh ($ aA)
NB.   vrios0 - monad to form rios(0); is called as:
NB.              rios0=. vrios0 ($ aA)
NB.   vstep  - dyad to perform single step of
NB.            Q-factorization; see ge*2step; is called as:
NB.              'aAi1 riosi1'=. drios vstep (aAi ; riosi)
NB.   vapp   - monad to do non-blocked Q-factorization; is
NB.            called as:
NB.              B=. vapp aA
NB.
NB. Algorithm:
NB.   1) k=. <./ vrsh ($ aA)
NB.   2) rios0=. vrios0 ($ aA)
NB.   3) B=. 0 {:: (k vstep (aA ; rios0))
NB. where
NB.   aA########################
NB.   B

geq2=: 1 : '((<./ @ (m gi 0) @ ]) (0 {:: (m gi 2)) (; (m gi 1))) $'

NB. ---------------------------------------------------------
NB. gelq2
NB. geql2
NB. geqr2
NB. gerq2
NB. emulate xGE{LQ,QL,QR,RQ}2
NB.
NB. TODO:
NB. - consider delay tau conjugation to final batch processing
NB. - ensure geq*2 accepts eA and returns eAupd

gelq2=: (0 _1 & +)`(_1 _1 & ,:)`(GEQ2DRIOS & gelq2step) geq2
geql2=: (_1 0 & +)`( 0  0 & ,:)`(GEQ2DRIOS & geql2step) geq2
geqr2=: (_1 0 & +)`(_1 _1 & ,:)`(GEQ2DRIOS & geqr2step) geq2
gerq2=: (0 _1 & +)`( 0  0 & ,:)`(GEQ2DRIOS & gerq2step) geq2

NB.   Vτ     - V and τ aggregated
NB.   V      - current block within input matrix
NB.   τ      - vector of τ[i] corresp. to V
NB.   T      - triangular factor of complex block reflector H
NB.   C      - not yet processed part of input matrix to
NB.            update
NB.   Vτlast - last piece of input matrix to apply
NB.            non-blocked code {lq2 ql2 qr2 rq2}

NB. rios0=. mkrios0gexxf m,n,bs,iter
NB. Vτ,V,τ,T,C,Vτlast

mkrios0gelqf=: 3 : 0
  'm n bs iter'=. y
  ibs=. bs*iter
  6 2 2 $ 0 0,bs,(n+1),0 0,bs,n,0,n,bs,1,0,(n+1),bs,bs,bs,0,(m-bs),n,ibs,ibs,(m-ibs),(n+1-ibs)
)

mkrios0geqlf=: 3 : 0
  'm n bs iter'=. y
  ibs=. bs*iter
  6 2 2 $ bs,_1,(m+1),bs,(bs+1),_1,m,bs,bs,_1 1,bs,0 _1,bs,bs,(bs+1),0,m,(n-bs),bs,0,(m+1-ibs),(n-ibs)
)

mkrios0geqrf=: 3 : 0
  'm n bs iter'=. y
  ibs=. bs*iter
  6 2 2 $ 0 0,(m+1),bs,0 0,m,bs,m,0 1,bs,(m+1),0,bs,bs,0,bs,m,(n-bs),ibs,ibs,(m+1-ibs),(n-ibs)
)

mkrios0gerqf=: 3 : 0
  'm n bs iter'=. y
  ibs=. bs*iter
  6 2 2 $ _1,bs,bs,(n+1),_1,(bs+1),bs,n,_1,bs,bs,1 _1 0,bs,bs,0,(bs+1),(m-bs),n,0,bs,(m-ibs),(n+1-ibs)
)

gelqfstep=: ([ ([ (3 1 4 larfbrnfr) (1 2 larftfr)) (gelq2 updir 0)) step
geqlfstep=: ([ ([ (3 1 4 larfblcbc) (1 2 larftbc)) (geql2 updir 0)) step
geqrfstep=: ([ ([ (3 1 4 larfblcfc) (1 2 larftfc)) (geqr2 updir 0)) step
gerqfstep=: ([ ([ (3 1 4 larfbrnbr) (1 2 larftbr)) (gerq2 updir 0)) step

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gelqf
NB. geqlf
NB. geqrf
NB. gerqf
NB. Ambivalent verbs to perform LQ QL QR RQ factorization
NB.
NB. Syntax:
NB.   LQf=. gelqf A
NB.   LQf=. (m,n,bs,iters) gelqf ((m,(n+1+bs)) {. A)
NB.   QfL=. geqlf A
NB.   QfL=. (m,n,bs,iters) geqlf (((-(m+1+bs)),n) {. A)
NB.   QfR=. geqrf A
NB.   QfR=. (m,n,bs,iters) geqrf (((m+1+bs),n) {. A)
NB.   RQf=. gerqf A
NB.   RQf=. (m,n,bs,iters) gerqf ((m,(n+1+bs)) {. A)
NB.
NB. Notes:
NB. - emulate xGE{LQ,QL,QR,RQ}F from LAPACK with following
NB.   differences:
NB.   1) ...
NB.
NB. TODO:
NB. - geqrfr & Co.: QR of submatrix defined by rIOS, when
NB.   space for τ and T is pre-allocated
NB. - geqrfri & Co.: QR of submatrix defined by rIOS
NB.   indirectly, when space for τ and T is pre-allocated
NB. - geqrfp & Co.: version using larfp to produce
NB.   non-negative diagonal
NB. - >>>>>>>>>>> ensure geq*f monad: A->LQf, dyad: (A,riosB)->Aupd

gelqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  bs=. GEQFBSMIN >. GEQFBSMAX <. k  NB. adjusted block size
  iters=. <. k % bs                 NB. # of iterations
  y=. (m , (n+1+bs)) {. y           NB. allocate space for τ and T, filled by zeros
  (m,n,bs,iters) gelqf y
:
  'bs iters'=. _2 {. x
  rios0=. mkrios0gelqf x
  drios=. bs * GELQFDRIOS
  iters (((0,(-bs)) & }.) @ ((1 {:: ]) (gelq2 updir 5) (0 {:: ])) @ (drios & gelqfstep)) (y ; rios0)
)

geqlf=: 3 : 0
  k=. <./ 'm n'=. $ y
  bs=. GEQFBSMIN >. GEQFBSMAX <. k  NB. adjusted block size
  iters=. <. k % bs                 NB. # of iterations
  y=. ((-(m+1+bs)) , n) {. y        NB. allocate space for τ and T, filled by zeros
  (m,n,bs,iters) geqlf y
:
  'bs iters'=. _2 {. x
  rios0=. mkrios0geqlf x
  drios=. bs * GEQLFDRIOS
  iters (((bs,0) & }.) @ ((1 {:: ]) (geql2 updir 5) (0 {:: ])) @ (drios & geqlfstep)) (y ; rios0)
)

geqrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  bs=. GEQFBSMIN >. GEQFBSMAX <. k  NB. adjusted block size
  iters=. <. k % bs                 NB. # of iterations
  y=. ((m+1+bs) , n) {. y           NB. allocate space for τ and T, filled by zeros
  (m,n,bs,iters) geqrf y
:
  'bs iters'=. _2 {. x
  rios0=. mkrios0geqrf x
  drios=. bs * GEQRFDRIOS
  iters ((((-bs),0) & }.) @ ((1 {:: ]) (geqr2 updir 5) (0 {:: ])) @ (drios & geqrfstep)) (y ; rios0)
)

gerqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  bs=. GEQFBSMIN >. GEQFBSMAX <. k  NB. adjusted block size
  iters=. <. k % bs                 NB. # of iterations
  y=. (m , (-(n+1+bs))) {. y        NB. allocate space for τ and T, filled by zeros
  (m,n,bs,iters) gerqf y
:
  'bs iters'=. _2 {. x
  rios0=. mkrios0gerqf x
  drios=. bs * GERQFDRIOS
  iters (((0,bs) & }.) @ ((1 {:: ]) (gerq2 updir 5) (0 {:: ])) @ (drios & gerqfstep)) (y ; rios0)
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. tgeqf
NB. Test orthogonal factorization algorithms:
NB. - non-blocked: gelq2 geql2 geqr2 gerq2
NB. - built-in: 128!:0
NB. - blocked from LAPACK: gelqf geqlf geqrf gerqf
NB. - blocked: gelqf geqlf geqrf gerqf
NB. by general matrix
NB.
NB. Syntax: tgeqf A
NB. where A - general m×n-matrix
NB.
NB. TODO:
NB. - split on tgeqf_{nonblocked,builtin,lapack,blocked}

tgeqf=: 3 : 0
  require '~addons/math/lapack/lapack.ijs'
  NB. need_jlapack_ 'gelqf geqlf geqrf gerqf'
  load '/home/jip/j602-user/projects/lapack/gelqf.ijs'
  load '/home/jip/j602-user/projects/lapack/geqlf.ijs'
  load '/home/jip/j602-user/projects/lapack/geqrf.ijs'
  load '/home/jip/j602-user/projects/lapack/gerqf.ijs'

  rcond=. ((_."_)`(norm1 con getir) @. (=/@$)) y

  ('gelq2' tmonad (,.&0 )`]`(rcond"_)`(_."_)`((norm1@(- ((         trl   @( 0 _1&}.)) mp  unglq)))%(FP_EPS*(#*norm1)@[))) y
  ('geql2' tmonad (0 &, )`]`(rcond"_)`(_."_)`((norm1@(- ((((-~/@$) trl ])@( 1  0&}.)) mp~ ungql)))%(FP_EPS*(#*norm1)@[))) y
  ('geqr2' tmonad (, &0 )`]`(rcond"_)`(_."_)`((norm1@(- ((         tru   @(_1  0&}.)) mp~ ungqr)))%(FP_EPS*(#*norm1)@[))) y
  ('gerq2' tmonad (0 &,.)`]`(rcond"_)`(_."_)`((norm1@(- ((((-~/@$) tru ])@( 0  1&}.)) mp  ungrq)))%(FP_EPS*(#*norm1)@[))) y

  ('128!:0' tmonad ]`]`(rcond"_)`(_."_)`((norm1@(- (mp & >/)))%(FP_EPS*(#*norm1)@[))) y

  ('2b1110&gelqf_jlapack_' tmonad ]`({. , <@((2&{::) ,.  (1&{::)))`(rcond"_)`(_."_)`((norm1@(- ((mp  unglq) & > /)))%((FP_EPS*#*norm1)@[))) y
  ('2b0111&geqlf_jlapack_' tmonad ]`({: , <@((0&{::) ,~  (1&{::)))`(rcond"_)`(_."_)`((norm1@(- ((mp~ ungql) & > /)))%((FP_EPS*#*norm1)@]))) y
  ('2b0111&geqrf_jlapack_' tmonad ]`({: , <@((0&{::) ,   (1&{::)))`(rcond"_)`(_."_)`((norm1@(- ((mp~ ungqr) & > /)))%((FP_EPS*#*norm1)@]))) y
  ('2b1110&gerqf_jlapack_' tmonad ]`({. , <@((2&{::) ,.~ (1&{::)))`(rcond"_)`(_."_)`((norm1@(- ((mp  ungrq) & > /)))%((FP_EPS*#*norm1)@]))) y

  ('gelqf' tmonad ]`]`(rcond"_)`(_."_)`((norm1@(- ((         trl   @( 0 _1&}.)) mp  unglq)))%((FP_EPS*#*norm1)@[))) y
  ('geqlf' tmonad ]`]`(rcond"_)`(_."_)`((norm1@(- ((((-~/@$) trl ])@( 1  0&}.)) mp~ ungql)))%((FP_EPS*#*norm1)@[))) y
  ('geqrf' tmonad ]`]`(rcond"_)`(_."_)`((norm1@(- ((         tru   @(_1  0&}.)) mp~ ungqr)))%((FP_EPS*#*norm1)@[))) y
  ('gerqf' tmonad ]`]`(rcond"_)`(_."_)`((norm1@(- ((((-~/@$) tru ])@( 0  1&}.)) mp  ungrq)))%((FP_EPS*#*norm1)@[))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testqf
NB. Test orthogonal factorizations
NB. Syntax: mkge testqf (m,n)

testqf=: 1 : 'EMPTY [ tgeqf @ u'
