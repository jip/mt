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

NB. Nouns, differences between rIOSs at consequent
NB. iterations: rios(i+1)-rios(i)

GELQ2DRIOS=: 2 2 $ 0 0 _1 _1  NB. reduce height and width by 1
GEQL2DRIOS=: 2 2 $ 0 0 _1 _1  NB. reduce height and width by 1
GEQR2DRIOS=: 2 2 $ 0 0 _1 _1  NB. reduce height and width by 1
GERQ2DRIOS=: 2 2 $ 0 0 _1 _1  NB. reduce height and width by 1

NB. ---------------------------------------------------------
NB. gerf0    Template conj. to make verbs to generate and
NB.          conditionally apply an elementary reflector to a
NB.          matrix
NB. gerf02   Template conj. to make verbs to generate and
NB.          conditionally apply an elementary reflector to a
NB.          matrix
NB.
NB. ---------------------------------------------------------
NB. gerf0
NB. Template conj. to make verbs to generate and
NB. conditionally apply an elementary reflector to a matrix
NB.
NB. Syntax:
NB.   vcondapp=. (vref`vapp) gerf0 ios
NB. where
NB.   ios      - 2-vector of integers (iosYZ,iosT),
NB.              indices in rIOS bundle (see layout below)
NB.   vref     - verb to reflect vector y (see larfg or
NB.              larfp)
NB.   vapp     - verb to apply an elementary reflector to
NB.              submatrix (rIOSs are taken from x, matrix is
NB.              in y) from either the left or the right; is
NB.              called when τ≠0
NB.   vcondapp - verb to do steps:
NB.                1) reflect vector ((0{x){y) by verb vref
NB.                   and store result into vector ((1{x){y)
NB.                2) if ((2{x){y)≠0 then call verb (x vapp y)
NB.
NB. Storage layout for rios:
NB.   iosYZ{rios - rIOS of vector y or z (see larfr or larfl)
NB.   iosT{rios  - rIOS of scalar τ (see larfr or larfl)
NB.
NB. Application:
NB. - let rios=. (riosYZ , riosT , riosR ,: riosL)
NB.   is rIOS bundle for vectors y (or z), scalar τ, and
NB.   submatrices R and L, respectively; vector y is stored
NB.   vertically, α is in head, τ is in tail
NB. - to reflect in the matrix A the vector y by larfg, then
NB.   to apply vector v from the left (emulate main loop of
NB.   xGEQR2 in LAPACK version 3.1.1):
NB.     vref=: 0 _1 & larfg                    NB. verb to reflect Y to Z
NB.     vapp=: 0 3 larfL                       NB. verb to apply v from the left to L
NB.     Aupd=. rios ((vref`vapp) gerf0 0 1) A  NB. apply gerund and ios to do the job
NB. - to do the same as above, but apply v from the right,
NB.   then from the left (emulate main loop of xGEHD2 in
NB.   LAPACK version 3.1.1):
NB.     vref=: 0 _1 & larfg                    NB. verb to reflect Y to Z
NB.     vapp=: 0 2 3 larfRL                    NB. verb to apply v from the right to R, then from the left to L
NB.     Aupd=. rios ((vref`vapp) gerf0 0 1) A  NB. apply gerund and ios to do the job
NB. - do the same as above, but with matrix shrinking:
NB.     vref=: 0 _1 & larfg                    NB. verb to reflect Y to Z
NB.     vapp=: 0 1 2 3 0 4 5 6 7 larfRLs       NB. verb to apply shrinked v from the right to shrinked R, then from the left to shrinked L
NB.     Aupd=. rios ((vref`vapp) gerf0 0 1) A  NB. apply gerund and ios to do the job

NB. ---------------------------------------------------------
NB. gerf02
NB. Template conj. to make verbs to generate and
NB. conditionally apply an elementary reflector to a matrix
NB.
NB. Syntax:
NB.   vcondapp=. (vref`vneg`vapp) gerf02 ios
NB. where
NB.   ios      - 2-vector of integers (iosYZ,iosT),
NB.              indices in rIOS bundle (see layout below)
NB.   vref     - verb to reflect vector y (see larfg or
NB.              larfp)
NB.   vneg     - verb to negate either a 1st row (left side
NB.              case), or a 1st column (right side case); is
NB.              called when τ=2; rIOSs are taken from x,
NB.              matrix is in y
NB.   vapp     - verb to apply an elementary reflector to a
NB.              submatrix from either the left or the right;
NB.              is called when τ≠0 and τ≠2; rIOSs are taken
NB.              from x, matrix is in y
NB.   vcondapp - verb to do steps:
NB.                1) reflect vector ((0{x){y) by verb vref
NB.                   and store result into vector ((1{x){y)
NB.                2a) if ((2{x){y)=0 then do nothing
NB.                2b) if ((2{x){y)=2 then call verb (x vneg y)
NB.                2c) othewise call verb (x vapp y)
NB.
NB. Storage layout for rios:
NB.   iosYZ{rios - rIOS of vector y or z (see larfr or larfl)
NB.   iosT{rios  - rIOS of scalar τ (see larfr or larfl)
NB.
NB. Application:
NB. - let rios=. (riosYZ , riosT , riosR ,: riosL)
NB.   is rIOS bundle for vectors y, z, scalar τ, and
NB.   submatrices R and L, respectively; vector y is stored
NB.   vertically, α is in head, τ is in tail
NB. - to reflect in the matrix A the vector y by larfp, then
NB.   either to negate 1st row of submatrix L (if τ=2), or to
NB.   apply shrinked vector v from the left (if τ≠0 and τ≠2)
NB.   to the shrinked submatrix L (emulate main loop of
NB.   xGEQR2 in LAPACK version 3.2):
NB.     vref=: 0 _1 & larfp                          NB. verb to reflect Y to Z
NB.     vneg=: 3 nfv 0                               NB. verb to negate L's 1st row
NB.     vapp=: 0 1 3 0 6 7 larfLs                    NB. verb to apply v from the left to L
NB.     Aupd=. rios ((vref`vneg`vapp) gerf02 0 1) A  NB. apply gerund and ios to do the job
NB. - to reflect in the matrix A the vector y by larfp, then
NB.   either to negate 1st column of submatrix R and then to
NB.   negate 1st row of submatrix L (if τ=2), or to apply
NB.   shrinked vector v to the shrinked submatrix R from the
NB.   right, then to the shrinked submatrix L from the left
NB.   (if τ≠0 and τ≠2) (emulate main loop of xGEHD2 in
NB.   LAPACK version 3.2):
NB.     vref=: 0 _1 & larfp                          NB. verb to reflect Y to Z
NB.     vneg=: [ (3 nfv 0) (2 nfv 1)                 NB. verb to negate R's 1st col, then L's 1st row
NB.     vapp=: 0 1 2 3 0 4 5 6 7 larfRLs             NB. verb to apply shrinked v from the right to shrinked R, then from the left to shrinked L
NB.     Aupd=. rios ((vref`vneg`vapp) gerf02 0 1) A  NB. apply gerund and ios to do the job
NB. - to do the same as previous, but without shrinking
NB.   (optimize speed for non-band matrices):
NB.     vref=: 0 _1 & larfp                          NB. verb to reflect Y to Z
NB.     vneg=: [ (3 nfv 0) (2 nfv 1)                 NB. verb to negate R's 1st col, then L's 1st row
NB.     vapp=: 0 2 3 larfRL                          NB. verb to apply v from the right to R, then from the left to L
NB.     Aupd=. rios ((vref`vneg`vapp) gerf02 0 1) A  NB. apply gerund and ios to do the job

NB. ---------------------------------------------------------
NB. geq2step
NB.
NB. Template adv. to form verbs ge*2step
NB.
NB. Syntax:
NB.   vstep=. vref`vmul`vtau geq2step iosyz
NB. where
NB.   vref  - monad to generate an elementary reflector, see
NB.           larfg* larfp*; is called as:
NB.             z=. vref y
NB.   vmul  - dyad to extract v from z and to product
NB.           (v*v'*C) or (C*v*v'); see larfl, larfr; is
NB.           called as:
NB.             eCupd=. eC vmul z
NB.   vstep - verb to perform single step of Q-factorization;
NB.           see ge*2step; is called as:
NB.             'Ai1 riosi1'=. drios vstep (Ai ; riosi)
NB.   iosyz - IOS vectors y and z in matrix eC

NB. geq2step=: 1 : '([ ((m gi 1) ^: (0 ~: (0 ({,) (1 fromci)))) ((m gi 0) updci 0)) step'     NB. gerf0
NB. geq2step=: 1 : '([ ((]`(}.u)) @. (0 2 i. (0 ({,) (1 fromci)))) ((m gi 0) updci 0)) step'  NB. gerf02
NB. geq2step=: 1 : '([ (m gi 1) ((m gi 0) updri 0)) step'                                     NB. gerf

geq2step=: 2 : '(((] (n}) ([ - (m gi 2) * (m gi 1))) ((m gi 0) @ (n&{))) updr) step'          NB. aggregated gerf

NB. ---------------------------------------------------------
NB. gelq2step
NB. geql2step
NB. geqr2step
NB. gerq2step
NB. Single step of corresp. method
NB.
NB. Syntax:
NB.   'Ai1 riosi1'=. drios gelq2step (Ai ; riosi)
NB.   'Ai1 riosi1'=. drios geql2step (Ai ; riosi)
NB.   'Ai1 riosi1'=. drios geqr2step (Ai ; riosi)
NB.   'Ai1 riosi1'=. drios gerq2step (Ai ; riosi)
NB. where
NB.   Ai    - (m+1)×n-matrix A(i) to update before i-th
NB.           iteration (i=0:min(m,n)-1)
NB.   riosi - 3×2×2-matrix rios(i) of rIOSs (see struct.ijs)
NB.           for i-th iteration; rows (0:2) contains:
NB.             0 - 2×2-matrix riosYZ for (m-i+1)-vector
NB.                   Y=(α[i],x[i][1:m-(i+1)],0)
NB.                 to reflect, and for vector
NB.                   Z=(β[i],v[i][1:m-(i+1)],τ[i])
NB.                 of reflection result
NB.             1 - 2×2-matrix riosT for scalar τ[i]
NB.             2 - 2×2-matrix riosL for (m-i)×(n-(i+1))-
NB.                 matrix L to apply an elementary reflector
NB.                 from the left
NB.   drios - difference between rIOSs at consequent
NB.           iterations: rios(i+1)-rios(i)
NB.   Ai1   - (m+1)×n-matrix A(i+1) after i-th iteration
NB.   riosi1 - 3×2-matrix rios(i+1) of rIOSs for (i+1)-th
NB.           iteration

NB.           (), performs actions:
NB.             1) extract A(i) and rios(i) and supply its to
NB.                gerf0; the last is configured to use
NB.                gerund (vref`vapp), to get vectors y(i)
NB.                and z(i)'s rIOS from rios(i)[0], scalar
NB.                τ(i)'s rIOS from rios(i)[1];
NB.             2) box output A(i+1) and write it into 0-th
NB.                item of input;
NB.             3) adjust rios(i) by Δrios;
NB.             4) box output rios(i+1) and write it into
NB.                1-th item of input.

NB.   z     - vector, reflection result of y
NB.   eC    - matrix C augmented by vector y and zero vector:
NB.             ((α,x,0),(C,.0)) - for LQ
NB.             (0,(C,.y)) - for QL
NB.             ((y,.C),0) - for QR
NB.             (0,.(C,y)) - for RQ
NB.   eCupd - updated matrix eC:
NB.             ((β,v,τ),((C*H),.0)) - for LQ
NB.             (0,(C,.y)) - for QL
NB.             ((y,.C),0) - for QR
NB.             (0,.(C,y)) - for RQ
NB.           τ

NB. (C,.y) before update, ((H'*C),.z) after update
NB. (y,.C) before update, (z,.(H'*C)) after update
NB. (C,y) before update, ((C*H),z) after update

NB. geqr2step=: (vref`vmul`vtau geq2step iosY) updr

gelq2step=: larfgfc`(larfr (+ @ (1 0 & (0 _1}))))`     (_1 { ])  geq2step IOSFR
geql2step=: larfgb `(larfl      (0 1 & (0 _1})) )`(+ @ ( 0 { ])) geq2step IOSLC
geqr2step=: larfgf `(larfl      (1 0 & (0 _1})) )`(+ @ (_1 { ])) geq2step IOSFC
gerq2step=: larfgbc`(larfr (+ @ (0 1 & (0 _1}))))`     ( 0 { ])  geq2step IOSLR

NB. ---------------------------------------------------------
NB. geq2
NB. Template adv. to form verbs ge*2

NB. Algo: ((min/ shape) (0 {:: ge*2step) (A ; cios0)) shape
geq2=: 1 : '((<./ @ (m gi 0) @ ]) (0 {:: (m gi 2)) (; (m gi 1))) $'

NB. ---------------------------------------------------------
NB. gelq2
NB. geql2
NB. geqr2
NB. gerq2
NB. emulate xGELQ2
NB. LQf=. gelq2 A

NB. TODO: delay tau conjugation to final batch processing

gelq2=: (0 _1 & +)`(_1 _1 & ,:)`(GELQ2DRIOS & gelq2step) geq2
geql2=: (_1 0 & +)`( 0  0 & ,:)`(GEQL2DRIOS & geql2step) geq2
geqr2=: (_1 0 & +)`(_1 _1 & ,:)`(GEQR2DRIOS & geqr2step) geq2
gerq2=: (0 _1 & +)`( 0  0 & ,:)`(GERQ2DRIOS & gerq2step) geq2

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gelqf
NB. geqlf
NB. geqrf
NB. gerqf
NB. emulate xGELQF

geqlf=: geql2  NB. stub for a while

QFBSMIN=: 2  NB. min block size for qf algorithms
QFBSMAX=: 3  NB. max block size for qf algorithms

NB. Vτ,V,τ,T,C,Vτlast; in 'block size' units

GELQFDRIOS=: 6 2 2 $ 1 1 0 _1 1 1 0 _1 1 0 0 0 0 0 0 0 1 1 _1 _1 0 0 0 0
GEQLFDRIOS=: 6 2 2 $ _1 _1 _1 0 _1 _1 _1 0 0 _1 0 0 0 0 0 0 0 0 _1 _1 0 0 0 0
GEQRFDRIOS=: 6 2 2 $ 1 1 _1 0 1 1 _1 0 0 1 0 0 0 0 0 0 1 1 _1 _1 0 0 0 0

NB. rios0=. mkrios0gexxf m,n,bs,iter

mkrios0gelqf=: 3 : 0
  'm n bs iter'=. y
  ibs=. bs*iter
  6 2 2 $ 0 0,bs,(n+1),0 0,bs,n,0,n,bs,1,0,(n+1),bs,bs,bs,0,(m-bs),n,ibs,ibs,(m-ibs),(n+1-ibs)
)

mkrios0geqlf=: 3 : 0
  'm n bs iter'=. y
  ibs=. bs*iter
  6 2 2 $ _1 _1,(m+1),bs,_1 _1,m,bs,bs,_1 1,bs,0,_1,bs,bs,(bs+1),0,m,(n-bs),bs,0,(m+1-ibs),(n-ibs)
)

mkrios0geqrf=: 3 : 0
  'm n bs iter'=. y
  ibs=. bs*iter
  6 2 2 $ 0 0,(m+1),bs,0 0,m,bs,m,0 1,bs,(m+1),0,bs,bs,0,bs,m,(n-bs),ibs,ibs,(m+1-ibs),(n-ibs)
)

gelqfstep=: ([ ([ (3 1 4 larfbrnfr) (1 2 larftfr)) (gelq2 updri 0)) step
geqlfstep=: (([ ([ (3 1 4 larfblcbc) (1 2 larftbc)) (geql2 updri 0)) dbg 'step') step
geqrfstep=: ([ ([ (3 1 4 larfblcfc) (1 2 larftfc)) (geqr2 updri 0)) step

NB. ambivalent:
NB. LQf=.                gelqf A
NB. LQf=. (m,n,bs,iters) gelqf ((m , (n+1+bs)) {. A)

gelqf=: 3 : 0
  k=. <./ 'm n'=. $ y
  bs=. QFBSMIN >. QFBSMAX <. k  NB. adjusted block size
  iters=. <. k % bs             NB. # of iterations
  y=. (m , (n+1+bs)) {. y       NB. allocate space for τ and T, filled by zeros
  (m,n,bs,iters) gelqf y
:
  'bs iters'=. _2 {. x
  rios0=. mkrios0gelqf x
  drios=. bs * GELQFDRIOS
  iters (((0,(-bs)) & }.) @ ((1 {:: ]) (gelq2 updri 5) (0 {:: ])) @ (drios & gelqfstep)) (y ; rios0)
)

NB. ambivalent:
NB. QfL=.                geqlf A
NB. QfL=. (m,n,bs,iters) geqlf (((-(m+1+bs)) , n) {. A)

geqlf=: 3 : 0
  k=. <./ 'm n'=. $ y
  bs=. QFBSMIN >. QFBSMAX <. k  NB. adjusted block size
  iters=. <. k % bs             NB. # of iterations
  y=. ((-(m+1+bs)) , n) {. y    NB. allocate space for τ and T, filled by zeros
  (m,n,bs,iters) geqlf y
:
  'bs iters'=. _2 {. x
  rios0=. mkrios0geqlf x
  drios=. bs * GEQLFDRIOS
  iters (((bs,0) & }.) @ ((1 {:: ]) (((geql2 dbg 'geql2') updri 5) dbg 'last QL2') (0 {:: ])) @ (drios & geqlfstep)) (y ; rios0)
)

NB. ambivalent:
NB. QfR=.                geqrf A
NB. QfR=. (m,n,bs,iters) geqrf (((m+1+bs) , n) {. A)

geqrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  bs=. QFBSMIN >. QFBSMAX <. k  NB. adjusted block size
  iters=. <. k % bs             NB. # of iterations
  y=. ((m+1+bs) , n) {. y       NB. allocate space for τ and T, filled by zeros
  (m,n,bs,iters) geqrf y
:
  'bs iters'=. _2 {. x
  rios0=. mkrios0geqrf x
  drios=. bs * GEQRFDRIOS
  iters ((((-bs),0) & }.) @ ((1 {:: ]) (geqr2 updri 5) (0 {:: ])) @ (drios & geqrfstep)) (y ; rios0)
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. tgeqf2
NB. Test orthogonal factorizations gelq2 geql2 geqr2 gerq2
NB. Syntax: tgeq2 A

tgeqf2=: 3 : 0
  'gelq2' ??? tqf (y ,. 0)
  'geql2' ??? tqf (0 ,  y)
  'geqr2' ??? tqf (y ,  0)
  'gerq2' ??? tqf (0 ,. y)
  EMPTY
)

NB. ---------------------------------------------------------
NB. tgeqf
NB. Test orthogonal factorizations gelqf geqlf geqrf gerqf
NB. Syntax: tgeqf A

tgeqf=: 3 : 0
  'gelqf'          ??? tqf y
  'gelqf_jlapack_' ??? tqf y
  'geqlf'          ??? tqf y
  'geqlf_jlapack_' ??? tqf y
  'geqrf'          ??? tqf y
  if. >:/ $ y do. 
    '128!:0' ??? tmonad y
  end.
  'geqrf_jlapack_' ??? tqf y
  'gerqf'          ??? tqf y
  'gerqf_jlapack_' ??? tqf y
  EMPTY
)

NB. ---------------------------------------------------------
NB. testqf
NB. Test orthogonal factorizations
NB. Syntax: testqf (m,n)

testqf=: 1 : 0
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'gelqf geqlf geqrf gerqf'
  T=. _.   NB. make triangular (trapezoidal) matrix with non-neg diag
  Qf=. _.  NB. make matrix w orth rows/cols in factorized form
  (tgeqf2 @ u) T;Qf
  (tgeqf @ u)  T;Qf
  EMPTY
)
