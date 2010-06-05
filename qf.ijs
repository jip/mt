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

GELQ2DRIOS=: 3 2 2 $ 1 0 0 _1 1 0 0 0 0 0 _1 _1    NB. z,τ,C*H
GEQL2DRIOS=: 3 2 2 $ 0 _1 _1 0 0 _1 0 0 0 0 _1 _1  NB. z,τ,H'*C
NB. ### GEQR2DRIOS=: 3 2 2 $ 0 1 _1 0 0 1 0 0 0 0 _1 _1    NB. z,τ,H'*C
GEQR2DRIOS=: 2 2 $ 0 0 _1 _1  NB. extended C to apply reflector
GERQ2DRIOS=: 3 2 2 $ _1 0 _1 0 0 0 0 _1 0 0 _1 _1  NB. z,τ,C*H

NB. ---------------------------------------------------------
NB. mkrios0gelq2
NB. mkrios0geql2
NB. mkrios0geqr2
NB. mkrios0gerq2
NB.
NB. Create rIOS at 0-th iteration for corresp. method
NB.
NB. Syntax:
NB.   rios0=. mkrios0gelq2 mn
NB.   rios0=. mkrios0geql2 mn
NB.   rios0=. mkrios0geqr2 mn
NB.   rios0=. mkrios0gerq2 mn
NB. where
NB.   mn    - 2-vector of integers (m,n), shape of matrix to
NB.           factorize
NB.   rios0 - 3×2×2-table rios(0), rIOSs corresponding to
NB.           iteration 0, see ge*2step verbs

mkrios0gelq2=: 3 : 0
  'm1 n1'=. _1 1 + 'm n'=. y
  3 2 2 $ 0 _1 1,n1,0 _1 1 1 _1 _2,m1,n
)

mkrios0geql2=: 3 : 0
  'm1 n1'=. 1 _1 + 'm n'=. y
  3 2 2 $ 0 _1,m1,1 0 _1 1 1 1 0,m,n1
)

mkrios0geqr2=: 3 : 0
  'm1 n1'=. 1 _1 + 'm n'=. y
  3 2 2 $ _1 0,m1,1 _1 0 1 1 _2 _1,m,n1
)

mkrios0gerq2=: 3 : 0
  'm1 n1'=. _1 1 + 'm n'=. y
  3 2 2 $ _1 0 1,n1,_1 0 1 1 0 1,m1,n
)

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
NB.   vstep=. vref`vapp geq2step
NB. where
NB.   vref  - verb to generate an elementary reflector, see
NB.           larfg* larfp*; is called as:
NB.             z=. vref y
NB.   vapp  - verb to apply an elementary reflector, see
NB.           larfL* larfR*; is called as:
NB.             Aupd=. rios vapp A
NB.   vstep - verb to perform single step of Q-factorization;
NB.           see ge*2step; is called as:
NB.             'Ai1 riosi1'=. drios gelq2step (Ai ; riosi)

NB. geq2step=: 1 : '([ ((m gi 1) ^: (0 ~: (0 ({,) (1 fromci)))) ((m gi 0) updci 0)) step'     NB. gerf0
NB. geq2step=: 1 : '([ ((]`(}.u)) @. (0 2 i. (0 ({,) (1 fromci)))) ((m gi 0) updci 0)) step'  NB. gerf02
NB. geq2step=: 1 : '([ (m gi 1) ((m gi 0) updri 0)) step'   NB. gerf
NB. geq2step=: 1 : '([ (m gi 1) ((m gi 0) updri 0)) step'   NB. gerfx

NB. ---------------------------------------------------------
NB. gelq2step
NB. geql2step
NB. geqr2step
NB. gerq2step
NB.
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

NB. ### gelq2step=: (larfgfc`(0 2 larfRfcs)) geq2step
NB. ### geql2step=: (larfgb`(0 2 larfLbsc)) geq2step
NB. ### geqr2step=: (larfgf`(0 2 larfLfsc)) geq2step
NB. ### gerq2step=: (larfgbc`(0 2 larfRbcs)) geq2step

NB. ---------------------------------------------------------
NB. Template adv. to form verbs ge*2

NB. Algo: (min/ shape) (0 {:: ge*2step) ((0 mix A) ; (mkcios0 shape))
NB. geq2=: 1 : '(<./ @ $) (0 {:: (m gi 2)) ((0 (m gi 0) ]) ; ((m gi 1) @ $))'

NB. Algo: ((min/ shape) (0 {:: ge*2step) (A ; cios0)) shape
geq2=: 1 : '((<./ @ (m gi 0) @ ]) (0 {:: (m gi 2)) (; (m gi 1))) $'

IOSFC=: < a: ; 0     NB. IOS 1st column
IOSLC=: < a: ; _1    NB. IOS last column
IOSFR=: < 0 ; < a:   NB. IOS 1st row
IOSLR=: < _1 ; < a:  NB. IOS last row

NB. - - - poligon - - - 

NB. QR2

NB. 'Ai1 riosi1'=. drios geqr2step (Ai ; riosi)
geq2step=: 2 : '((((((] ((n}) dbg ''wr Z'') (([ - ((m gi 3) dbg ''vtau'') * (((m gi 1) dbg ''vmul'') ((m gi 2) dbg ''vv''))) dbg ''H*C'')) dbg ''f(C,Z)'') ((m gi 0) @ ((n&{) dbg ''get Y for ref''))) dbg ''vjob'') updr) dbg ''updr'') step'

NB. geqr2step=: (vref`vmul`vv`vtau geq2step iosY) updr
geqr2step=: larfgf`larfl`(0 & (_1}))`(_1 { ]) geq2step IOSFC

NB. - - - /poligon - - - 

NB. ---------------------------------------------------------
NB. gelq2
NB. geql2
NB. geqr2
NB. gerq2
NB. emulate xGELQ2
NB. LQf=. gelq2 A

NB. ### gelq2=: (0 _1 & +)`mkrios0gelq2`(GELQ2DRIOS & gelq2step) geq2
NB. ### geql2=: (_1 0 & +)`mkrios0geql2`(GEQL2DRIOS & geql2step) geq2
geqr2=: (_1 0 & +)`(_1 _1 & ,:)`(GEQR2DRIOS & geqr2step) geq2
NB. ### gerq2=: (0 _1 & +)`mkrios0gerq2`(GERQ2DRIOS & gerq2step) geq2

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gelqf
NB. geqlf
NB. geqrf
NB. gerqf
NB. emulate xGELQF

gelqf=: gelq2  NB. stub for a while
geqlf=: geql2  NB. stub for a while

QFBSMIN=: 2  NB. min block size for qf algorithms
QFBSMAX=: 3  NB. max block size for qf algorithms

GEQRFDRIOS=: 3 2 2 $ 1 1 _1 0 0 0 0 0 1 1 _1 _1  NB. Z,T,CbyH; divided by block size

NB. rios0=. mkrios0geqrf m,n,bs
mkrios0geqrf=: 3 : 0
  'm n bs'=. y
  3 2 2 $ 0 0,m,bs,(m+1),0,bs,bs,0,bs,m,(n-bs)
)

NB. ### geqrfstep=: ([ ([ ((((0 1 2 larfblcfc) dbg 'larfB') upd3ri 0 1 2) dbg 'H''*C') ((larftfc mapri 0 1) dbg 'larfT')) ((geqr2 updri 0) dbg 'QR2')) dbg 'step'

geqrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  bs=. QFBSMIN >. QFBSMAX <. k  NB. adjusted block size
  i=. <. k % bs                 NB. # of iterations
  y=. (m+1+bs) {. y             NB. allocate space for τ and T, filled by zeros
  rios0=. mkrios0geqrf m,n,bs
  drios=. bs * GEQRFDRIOS
  ((-(bs+1)),0) }. 0 {:: i (drios & geqrfstep) (y ; rios0)
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
