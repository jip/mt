NB. Triangular factorization
NB.
NB. getrfxxxx  Triangular factorization with partial pivoting
NB.            of a general matrix
NB. hetrfpx    Triangular factorization with full pivoting of
NB.            a Hermitian (symmetric) matrix
NB. potrfx     Cholesky factorization of a Hermitian
NB.            (symmetric) positive definite matrix
NB. pttrfx     Triangular factorization of a Hermitian
NB.            (symmetric) positive definite tridiagonal
NB.            matrix
NB.
NB. testgetrf  Test getrfxxxx by general matrix given
NB. testhetrf  Test hetrfpx by Hermitian (symmetric) matrix
NB.            given
NB. testpotrf  Test potrfx by Hermitian (symmetric) positive
NB.            definite matrix given
NB. testpttrf  Test pttrfx by Hermitian (symmetric) positive
NB.            definite tridiagonal matrix given
NB. testtrf    Adv. to make verb to test xxtrfxxxx by matrix
NB.            of generator and shape given
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

TRFNB=: 3  NB. block size limit, >0

NB. ---------------------------------------------------------
NB. hetf2pl
NB.
NB. Description:
NB.   Partial triangular factorization with full pivoting of
NB.   a Hermitian (symmetric) matrix:
NB.     P * L1 * T * L1^H * P^_1 = A
NB.   by non-blocked version of Aasen's algorithm
NB.
NB. Syntax:
NB.     'ip L1 T'=. (l1 ; nb) hetf2pl A
NB. where
NB.   A   - n×n-matrix to factorize, Hermitian (symmetric)
NB.   l1  - n-vector, 1st column of L1
NB.   nb  > 0, columns quantity to factorize
NB.   ip  - n-vector, full inversed permutation of A
NB.   L1  - n×n-matrix, unit lower triangular
NB.   T   - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   P   - n×n-matrix, full permutation of A
NB.
NB. Assertion:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A -: clean p sp T (] mp (mp ct)) L1
NB.   A -: clean P (mp mp |:@[) T (] mp (mp ct)) L1
NB. where
NB.   'ip L1 T'=. hetf2pl A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.
NB. References:
NB. [1] Miroslav Rozloznik, Gil Shklarski, Sivan Toledo.
NB.     Partitioned triangular tridiagonalization.
NB.     26 September 2007.
NB.     http://www.cs.cas.cz/miro/rst08.pdf
NB.
NB. TODO:
NB. - T would be sparse

appendl=: , `([, (({."1~    ({:@$) )~))`(({."1~    ({:@$) ), ]) @. (*@-&({:@$))
appendr=: , `([, (({."1~ (-@({:@$)))~))`(({."1~ (-@({:@$))), ]) @. (*@-&({:@$))

stitcht=: ,.`([,.(({.  ~        #  )~))`(({.~          #  ),.]) @. (*@-&    # )
stitchb=: ,.`([,.(({.  ~ (-@    # ))~))`(({.~ (-@      # )),.]) @. (*@-&    # )

NB. ---------------------------------------------------------
NB. lahefpl
NB.
NB. Description:
NB.   Partial factorization of a Hermitian (symmetric) matrix
NB.   using the combination of Parlett and Reid, and Aasen
NB.   methods [1].
NB.
NB. Syntax:
NB.   'ipo L1o l0o l1o t0o t1o Ao Ho h0o h1o'=. lahefpl (ipi;L1i;l0i;l1i;t0i;t1i;Ai;Hi;h0i;h1i)
NB. where
NB.   ipi  - (n-k)-vector (i. (n-k)), pre-allocated space for
NB.          inversed permutation after i-th step of
NB.          partitioned algorithm and before (i+1)-th one
NB.   L1i  - (n-k)×1-matrix L1 after i-th step of partitioned
NB.          algorithm  and before (i+1)-th one, unit lower
NB.          triangular
NB.   l0i  - ##########(n-k)-vector (i. (n-k)), pre-allocated space for
NB.          inversed permutation after i-th step of
NB.          partitioned algorithm and before (i+1)-th one
NB.   Ai   - (n-k)×(n-k)-matrix A after i-th step and before
NB.          (i+1)-th one, Hermitian (symmetric)
NB.   li   - (n-j)-vector, last column under diagonal in
NB.          L1(i)
NB.
NB.  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:
NB.        ( 0  U22 ) (  0   D  ) ( U12' U22' )
NB.
NB.  A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  if UPLO = 'L'
NB.        ( L21  I ) (  0  A22 ) (  0    I   )
NB.
NB. Storage layout:
NB.   
NB.
NB.
NB.
NB.
NB.
NB. Algorithm:
NB.
NB.
NB.
NB.
NB.   'ipi1 L1i1 Ti1 Ai1 Hi1 hi1 li1'=. lahefplstep (ipi;L1i;Ti;Ai;Hi;hi;li)
NB.
NB.
NB. References:
NB. [1] Miroslav Rozloznik, Gil Shklarski, Sivan Toledo.
NB.     Partitioned triangular tridiagonalization.
NB.     26 September 2007.
NB.     http://www.cs.cas.cz/miro/rst08.pdf
NB. 

NB.    clean lahefplc_mt_ ((i.10);(0 {."1 HE10);(10 ($!.0) 1);(,0);(0 ({,) HE10);(i.0);HE10;(0 {."1 HE10);({."1 HE10);((< (<0);0) { HE10))
NB.    clean lahefpur_mt_ ((i.10);(0 {. HE10);(10 ($!.0) 1);(,0);(0 ({,) HE10);(i.0);HE10;(0 {. HE10);({. HE10);((< 0;(<<0)) { HE10))
NB. 
NB.    'ip U1 T'=. clean hetrfpur_mt_ HE10
NB.    HE10 -: clean (/: ip) sp_mt_ U1 (mp~ mp~ (ct_mt_ @ [)) T
NB. 1
NB.    'ip L1 T'=. clean hetrfplc_mt_ HE10
NB.    HE10 -: clean (/: ip) sp_mt_ L1 (mp mp (ct_mt_ @ [)) T
NB. 1
NB.    HE1000=. (_1 1 0 4 _6 4 & (gemat_mt_ j. gemat_mt_)) hemat_mt_ 1000 1000

NB. column-wise traversing

lahefplc=: (3 : 0) ^: (TRFNB<.(#@(0 & {::)))
  'ip L1 l0 l1 t0 t1 A H h0 h1'=. y
  L1=. L1 stitchb l0
  l0=. }. l0
  'n j'=. $ L1
  l1=. h1 - l0 * ({: t0)
  H=. H stitchb h0                NB. h0 instead of (({: t1) , h0)
  dip0=. < 0 , ((i.>./)@soris) l1  NB. non-standard cycle permutation!
  dip=. j (+ &. >) dip0            NB. non-standard cycle permutation!
  ip=. dip (C. :: ]) ip
  L1=. dip (C. :: ]) L1
  l0=. dip0 (C. :: ]) l0
  l1=. dip0 (C. :: ]) l1
  t01=. + t10=. (0&{ :: ]) l1
  l1=. l1 % t10                    NB. 1 at head is not guaranteed!
  A=. dip (sp :: ]) A
  H=. dip (C. :: ]) H
  h0=. ((j (] dhs2lios (((* >:)),-~)) n) ({,) A) - (j }. H) mp (+ j ({ :: 0:) L1)  NB. if j is n then use 0: instead of (j{L1)
  h1=. h0 - l0 * t01
  t11=. 9 o. 0 ({ :: ]) h1         NB. CHECKME: is Re() necessary?; if h1 is (i.0) then use h1 instead of ({.h1)
  t0=. t0 , t11
  t1=. t1 , t01
  ip ; L1 ; l1 ; l0 ; t0 ; t1 ; A ; H ; h0 ; (}. h1)
)

NB. 'ip L1 T'=. hetrfpl A
NB. 'ip L1 T permA H'=. step (ip;L1;T;A;H)

NB. 'ipi L1i Ti Ai Hi hi'=. (9 # 0) hetf2pl_mt_ HEci10
NB. ((((] dhs2lios (_1,-))/@$) ({,) ]) L1i) hetf2pl_mt_ (((_1 ({,) Ti) , hi) (< a: ; 0)} (2 2 }. HEci10))

hetrfplc=: 3 : 0
  n=. # y
  ip=. i. n
  t1=. i. 0
  L1=. 0 {."1 y
  l0i=. n ($!.0) 1
  t0=. 0 ({,) y
  for_k. n (] dhs2lios (0,(>.@%))) TRFNB do.
    'ipi L1i l0i l1i t0i t1i y Hi h0i h1i'=. lahefplc ((i. # y);(0 {."1 y);l0i;(,0);({:t0);(i.0);y;(0 {."1 y);({."1 y);((< (<0);0) { y))
    dip=. k (i.@[ , +) ipi         NB. force permutation to act on tail part only
    ip=. dip C. ip
    L1=. (dip C. L1) stitchb L1i
NB.    l0i=. 1 (0}) l0i            NB. to guarantee 1 at head
    t0=. t0 , (}. t0i)
    t1=. t1 , t1i
    y=. ((2 # TRFNB) }. y) - ((Hi ((mp ct) &(TRFNB&}.)) L1i) + ((l1i * {: t1i) */ + l0i))
  end.
  T=. ((+t1);_1) setdiag (t1;1) setdiag (t0;0) setdiag (2 # n) $ 0
  ip ; L1 ; T
)

NB. row-wise traversing

NB. clean lahefpur_mt_ ((i.5);(0 {. HE5);(5 ($!.0) 1);(,0);(0 ({,) HE5);(i.0);HE5;(0 {. HE5);({. HE5);((< 0;(<<0)) { HE5))

lahefpur=: ((3 : 0) dbg 'step') ^: (TRFNB<.(#@(0 & {::)))
  'ip U1 u0 u1 t0 t1 A H h0 h1'=. y
  U1=. U1 appendr u0
  u0=. }. u0
  'j n'=. $ U1
  u1=. h1 - u0 * ({: t0)
  H=. H appendr h0                 NB. h0 instead of (({: t1) , h0)
  dip0=. < 0 , ((i.>./)@soris) u1  NB. non-standard cycle permutation!
  dip=. j (+ &. >) dip0            NB. non-standard cycle permutation!
  ip=. dip (C. :: ]) ip
  U1=. dip (C."1 :: ]) U1
  u0=. dip0 (C. :: ]) u0
  u1=. dip0 (C. :: ]) u1
  t10=. + t01=. (0&{ :: ]) u1
  u1=. u1 % t01                    NB. 1 at head is not guaranteed!
  A=. dip (sp :: ]) A
  H=. dip (C."1 :: ]) H
  h0=. ((j (dhs2lios@(((* >:)),-~)) n) (({,) dbg 'arg0') A) (- dbg 'h0=-') (+ j ({"1 :: 0:) U1) (mp dbg 'h0=mp') (j }."1  H)  NB. if j is n then use 0: instead of (j {"1 U1)
  h1=. h0 - t10 * u0
  t11=. 9 o. 0 ({ :: ]) h1         NB. CHECKME: is Re() necessary?; if h1 is (i.0) then use h1 instead of ({.h1)
  t0=. t0 , t11
  t1=. t1 , t10
  ip ; U1 ; u1 ; u0 ; t0 ; t1 ; A ; H ; h0 ; (}. h1)
)

NB. 'ip U1 T'=. hetrfpu A
NB. 'ip U1 T permA H'=. step (ip;U1;T;A;H)

NB. 'ipi U1i Ti Ai Hi hi'=. (9 # 0) hetf2pl_mt_ HEci10
NB. ((((] dhs2lios (_1,-))/@$) ({,) ]) U1i) hetf2pl_mt_ (((_1 ({,) Ti) , hi) (< a: ; 0)} (2 2 }. HEci10))

hetrfpur=: 3 : 0
  n=. # y
  ip=. i. n
  t1=. i. 0
  U1=. 0 {. y
  u0i=. n ($!.0) 1
  t0=. 0 ({,) y
  for_k. n (] dhs2lios (0,(>.@%))) TRFNB do.
    'ipi U1i u0i u1i t0i t1i y Hi h0i h1i'=. (lahefpur dbg 'lahefpur') ((i. # y);(0 {. y);u0i;(,0);({:t0);(i.0);y;(0 {. y);({. y);((< 0;(<<0)) { y))
    dip=. k (i.@[ , +) ipi         NB. force permutation to act on tail part only
    ip=. dip C. ip
    U1=. (dip C."1 U1) appendr U1i
NB.    u0i=. 1 (0}) u0i            NB. to guarantee 1 at head
    t0=. t0 , (}. t0i)
    t1=. t1 , t1i
    y=. ((2 # TRFNB) }. y) (- dbg '-') ((Hi (((mp dbg 'mp')~ ct) &(TRFNB&(}."1))) U1i) (+ dbg '+r1u') (((+ u0i) (*/ dbg '*/r1u') ({: t1i) (* dbg '*r1u') u1i)))
  end.
  T=. ((+t1);1) setdiag (t1;_1) setdiag (t0;0) setdiag (2 # n) $ 0
  ip ; U1 ; T
)

NB. ---------------------------------------------------------
NB. row-wise traversing with fused A,L1,H

NB. ---------------------------------------------------------

NB. clean lahefplro_mt_ ((i.5);HE5;(4 # 0);(,0);(0 ({,) HE5);(i.0);({. HE5);(}. {. HE5))

lahefplro=: ((3 : 0) dbg 'step') ^: (TRFNB<.(#@(2 & {::)))
  'ip A l0 l1 t0 t1 h0 h1'=. y     NB. n n*n n-(j+1) n-(j+1) k+(j+1) k+j n-j
  A=. l0 ((1 liosS)&c) } A             NB. n dhs2lios j*(n+1)+n,n-(j+1)
  l1=. (+ h1) - l0 * {: t0
  A=. h0 ((0 liosE)&c) } A             NB.   dhs2lios j*(n+1),n-j
  dip0=. < 0 , ((i.>./)@soris) l1  NB. non-standard cycle permutation!
  dip=. (A -&# l0) (+ &. >) dip0   NB. non-standard cycle permutation!
  ip=. dip (C. :: ]) ip
  l0=. dip0 (C. :: ]) l0
  l1=. dip0 (C. :: ]) l1
  t01=. + t10=. (0&{ :: ]) l1
  l1=. l1 % t10                    NB. 1 at head is not guaranteed!
  A=. dip (sp :: ]) A
  h0=. l0 (((((0 liosE) dbg '0liosE')&c) (({,) dbg 'arg0') ]) (- dbg 'h0=-') ((((-~ (1 liosW) ]) dbg '1liosW')&c) ({,) ]) (mp dbg 'h0=mp') (((-~ , -@[)&#) {. ])) A
  h1=. h0 - t10 * + l0
  t11=. 9 o. 0 ({ :: ]) h1         NB. CHECKME: is Re() necessary?; if h1 is (i.0) then use h1 instead of ({.h1)
  t0=. t0 , t11
  t1=. t1 , t10
  ip ; A ; (}. l1) ; l0 ; t0 ; t1 ; h0 ; (}. h1)
)

NB. ---------------------------------------------------------

NB. 'ip A l0'=. lahefplr_mt_ ip;A;l0
NB. clean lahefplr_mt_ (i.n);HE;(n ($!.0) 1)

lahefplr=: (3 : 0) ^: (TRFNB<.(#@(2 & {::)))
  'ip A l0'=. y                    NB. n n*n n-j, j=0:n-1
  j=. A -&# l0
  a=. j { A
  h0=. (j }. a) - (j {. a) mp ((j,(-#l0)) {. A)
  h1=. h0 - ({. l0) * (+ l0 (((_1 liosS)&c ({,) ]) :: (0 (0}) [)) A)
  l1=. (+ }. h1) - (}. l0) * ({. h1)
  dip0=. < 0 , ((i.>./)@soris) l1  NB. non-standard cycle permutation!
  dip=. (>:j) +&.> dip0            NB. non-standard cycle permutation!
  ip=. dip C. :: ] ip
  l1=. dip0 C. :: ] l1
  A=. ((({. h1) 0 } h0) , l0) ((((0 liosE),(] (((-~ {.)`0:`])}) (0 liosS)))~ -:)~&c } :: ((}.~ (>:@-:@#))@[ ((_1 liosS)&c }) ])) A
  A=. dip sp :: ] A
  l1=. 0 ({`[`(] % {))} :: ] l1
  ip ; A ; l1
)

lahefplrd=: ((3 : 0) dbg 'step') ^: (TRFNB<.(#@(2 & {::)))
  'ip A l0'=. y                    NB. n n*n n-(j+?), j=1:?
  j=. A ((-&#) dbg 'j') l0
  a=. j ({ dbg 'a') A
  h0=. (j }. a) (- dbg 'h0=asfx-') (j {. a) (mp dbg 'apfx mp subH') ((j,(-#l0)) {. A)
  h1=. h0 (- dbg 'h1=}.h0-') ({. l0) (* dbg '({.*+}.)l0') (+ l0 (((((_1 liosS)&c) dbg 'lios(l0)') (({,) dbg '{,') ]) :: ((0 (0}) [) dbg 'FAILED0')) A)
  l1=. (+ }. h1) (- dbg 'l1=+h1-') (}. l0) (* dbg '}.l0 * {.h1') ({. h1)
  dip0=. < 0 , (((i.>./)@soris) dbg 'p') l1  NB. non-standard cycle permutation!
  dip=. (>:j) ((+ &. >) dbg 'j,dip0 -> dip') dip0   NB. non-standard cycle permutation!
  ip=. dip ((C. :: ]) dbg 'C.ip') ip
  l1=. dip0 ((C. :: ]) dbg 'C.l1') l1
  A=. ((({. h1) 0 } h0) , l0) ((((((((0 liosE) dbg '0liosE'),(] ((((-~ {.)`0:`])}) dbg '}liosS') (0 liosS)))~ -:)~&c) }) dbg 'l0,h0->A') :: (((}.~ (>:@-:@#))@[ (((_1 liosS)&c) }) ]) dbg 'FAILED1')) A NB. 
  A=. dip ((sp :: ]) dbg 'sp(A)') A
  l1=. 0 ((({`[`(] % {))} :: ]) dbg 'l1[1:] /= l1[0]') l1
  ip ; A ; l1
)

NB. 'ip L1 T'=. hetrfplr A

hetrfplr=: 3 : 0
  n=. # y
  ip=. i. n
  l0=. n ($!.0) 1
  for. i. n (>.@%) TRFNB do.
    ios=. TRFNB (rt (; dbg 'ios') }.) y ([ th2lios -)&# l0
    'ip y l0'=. (lahefplrd dbg 'lahefplr') (ip;y;l0)
    subA=. ios ((<@(1 1{[)){]) y
    subL=. ios ((<@(1 0{[)){]) y
    subH=. ios ((<@(0 1{[)){]) y
    l1=. l0 ((_1 liosS)&c ({,) ]) y
    rank1upd=. ((0 { :: ] l0) (+@* dbg 't10*l0') (1 (0}) :: ] l0)) (*/ dbg 'conj(t10l0)*/l1') l1
    subAupd=. subA (- dbg '-') ((subL (mp dbg 'mp') subH) (+ dbg 'LH+r1u') rank1upd)
    y=. subAupd (((< 1 1 { ios) }) dbg 'subA}A') y
    NB. y=. ios (((<@(0 0{[)){]) (- dbg '-') (((<@(0 1{[)){]) (mp dbg 'mp') ((<@(1 0{[)){])))`(<@(0 0{[))`] } y
  end.
  L1=. trl1 y
  t0=. diag y
  t1=. 1 diag y
  T=. ((+t1);1) setdiag (t1;_1) setdiag (t0;0) setdiag (2 # n) $ 0
  ip ; L1 ; T
)


NB. ---------------------------------------------------------
NB. hetf2pu
NB.
NB. Description:
NB.   Triangular factorization with full pivoting of a
NB.   Hermitian (symmetric) matrix:
NB.     P * U1 * T * U1^H * P^_1 = A
NB.   by non-blocked version of Aasen's algorithm
NB.
NB. Syntax:
NB.     'ip U1 T'=. hetf2pu A
NB. where
NB.   A   - n×n-matrix to factorize, Hermitian (symmetric)
NB.   ip  - n-vector, full inversed permutation of A
NB.   P   - n×n-matrix, full permutation of A
NB.   U1  - n×n-matrix, unit upper triangular
NB.   T   - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.
NB. Assertion:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A -: clean p sp T (] mp (mp ct)) U1
NB.   A -: clean P (mp mp |:@[) T (] mp (mp ct)) U1
NB. where
NB.   'ip U1 T'=. hetf2pu A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.
NB. TODO:
NB. - T would be sparse

hetf2pu=: 3 : 0
  n=. # y
  U1=. (- $ y) {. 1                NB. last U1's column in Aasen method is e(n)-vector
  T=. ((- @ $) {. (_1 _1&{.)) y    NB. T[n-1,n-1]=A[n-1,n-1]
  ip=. i. n
  h=. i. 1                         NB. h[0:n-j-1]=H[j+1:n-1,j+1], may be defined arbitrary before the 1st iteration only
  ios=. i. <: n                    NB. 0:j
  for_j. i. 1 - n do.              NB. n-2:0
    a=. (< ios ; (>: j)) { y       NB. A[0:j,j+1]
    lum=. ((>:j) ([ , -) n) {. U1  NB. U1[0:j,j+1:n-1]
    v=. a - lum mp h               NB. v[0:j]=A[0:j,j+1]-U1[0:j,j+1:n-1]*H[j+1:n-1,j+1]=U1[0:j,j]*T[j,j+1] non-pivoted yet
    q=. liolmax v                  NB. IO pivot from tail
    dip=. < q , j                  NB. any[j]↔any[q], non-standard cycle permutation!
    v=. dip (C. :: ]) v            NB. v[_1]↔v[q]
    y=. dip (sp :: ]) y            NB. A[j,0:j]↔A[q,0:j], A[0:j,j]↔A[0:j,q]
    ip=. dip (C. :: ]) ip          NB. ip[j]↔ip[q]
    to=. {: v                      NB. T[j,j+1]
    lu=. q { lum                   NB. U1[j,j+1:n-1] after pivoting
    luecto=. to * + {. lu          NB. conj(U1[j,j+1])*T[j,j+1]
    h=. to (((+ +)~ {.)`0:`]) } ((2 # >: j - n) {. T) mp + lu  NB. h[0:n-j-1]=H[j+1:n-1,j]=T[j+1:n-1,j+1:n-1]*conj(U1[j,j+1:n-1])
    U1=. (1 (_1}) v % to) (< ios ; j) } dip C. U1              NB. U1[j,0:n-1]↔U1[q,0:n-1], U1[0:j,j]=v[0:j]/v[_1], v[_1] may be 0
    td=. 9 o. ((< 2 # j) { y) - (luecto + lu mp h)             NB. T[j,j]=Re(A[j,j]-U1[j,j+1:n-1]*H[j+1:n-1,j]-conj(U1[j,j+1])*T[j,j+1])
    T=. (td (,~ (, +)) to) (_2 <\ 0 1 1 0 0 0 + j) } T         NB. TODO: amend by lIOS NB. batch write diagonal and off-diagonal elements T[j,j+1] T[j+1,j] T[j,j]
    h=. h ,~ luecto + td           NB. h[0:j]=H[0:j,j]=T[0:J,0:J]*conj(U1[j,0:j])
    ios=. }: ios                   NB. 0:j-1
  end.
  ip ; U1 ; T
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. getrflu1p
NB.
NB. Description:
NB.   Triangular factorization with partial pivoting of a
NB.   general matrix:
NB.     L * U1 * P = A
NB.
NB. Syntax:
NB.     'ip LU1'=. getrflu1p A
NB. where
NB.   A   - m×n-matrix to factorize
NB.   ip  - n-vector, columns inversed permutation of A
NB.   LU1 - m×n-matrix, lower triangle contains L, and strict
NB.         upper triangle contains U1 without unit diagonal
NB.   P   - n×n-matrix, columns permutation of A
NB.   L   - m×k-matrix, lower triangular
NB.   U1  - k×n-matrix, unit upper triangular
NB.   k   = min(m,n)
NB.
NB. Assertion:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A -: clean p {"1 L mp U1
NB.   A -: clean p C."1 L mp U1
NB.   A -: clean ip C.^:_1"1 L mp U1
NB.   A -: clean L mp U1 mp P
NB. where
NB.   'ip LU1'=. getrflu1p A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.   L=. trl LU1
NB.   U1=. tru1 LU1

getrflu1p=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. n) ; y
  elseif. 1=m do.
    dip=. < 0 , liofmax {. y                                       NB. non-standard cycle permutation!
    ip=. dip (C. :: ]) i. n
    y=. ((] 0:} %) (0&({,))) dip ((C."1) :: ]) y                   NB. permute single row, scale by head, keep head unscaled
    ip ; y
  elseif. do.
    k=. n (<. >.@-:) m
    'pia Afa'=. getrflu1p k {. y                                   NB. factorize 1st block recursively
    y=. pia (C."1) k }. y                                          NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trsmxu1 & (k & ({."1))) y                          NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. getrflu1p y ((- (Afba & mp)) & (k & (}."1))) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. (i. k) , (k + pib)                                      NB. apply 2nd block's permutation to 1st block
    (dpib (C."1) pia) ; ((dpib (C."1) Afa) , (Afba ,. Afbb))       NB. assemble solution
  end.
)

NB. ---------------------------------------------------------
NB. getrfpl1u
NB.
NB. Description:
NB.   Triangular factorization with partial pivoting of a
NB.   general matrix:
NB.     P * L1 * U = A
NB.
NB. Syntax:
NB.     'ip L1U'=. getrfpl1u A
NB. where
NB.   A   - m×n-matrix to factorize
NB.   ip  - m-vector, rows inversed permutation of A
NB.   L1U - m×n-matrix, upper triangle contains U, and strict
NB.         lower triangle contains L1 without unit diagonal
NB.   P   - n×n-matrix, rows permutation of A
NB.   L1  - m×k-matrix, unit lower triangular
NB.   U   - k×n-matrix, upper triangular
NB.   k   = min(m,n)
NB.
NB. Assertion:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A -: clean p { L1 mp U
NB.   A -: clean p C. L1 mp U
NB.   A -: clean ip C.^:_1 L1 mp U
NB.   A -: clean P mp L1 mp U
NB. where
NB.   'ip L1U'=. getrfpl1u A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.   L1=. trl1 L1U
NB.   U=. tru L1U
NB.
NB. Notes:
NB. - models LAPACK's xGETRF

getrfpl1u=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. m) ; y
  elseif. 1 = n do.
    dip=. < 0 , liofmax y                                      NB. non-standard cycle permutation!
    ip=. dip (C. :: ]) i. m
    y=. ((] 0:} %) (0&({,))) dip (C. :: ]) y                   NB. permute single column, scale by head, keep head unscaled
    ip ; y
  elseif. do.
    k=. m (<. >.@-:) n
    'pia Afa'=. getrfpl1u k {."1 y                             NB. factorize 1st block recursively
    y=. pia C. k }."1 y                                        NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trsml1x & (k & {.)) y                          NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. getrfpl1u y ((- (mp & Afba)) & (k & }.)) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. (i. k) , (k + pib)                                  NB. apply 2nd block's permutation to 1st block
    (dpib C. pia) ; ((dpib C. Afa) ,. (Afba , Afbb))           NB. assemble solution
  end.
)

NB. ---------------------------------------------------------
NB. getrfpu1l
NB.
NB. Description:
NB.   Triangular factorization with partial pivoting of a
NB.   general matrix:
NB.     P * U1 * L = A
NB.
NB. Syntax:
NB.     'ip U1L'=. getrfpu1l A
NB. where
NB.   ip  - m-vector, rows inversed permutation of A
NB.   A   - m×n-matrix to factorize
NB.   U1L - m×n-matrix, lower triangle contains L, and strict
NB.         upper triangle contains U1 without unit diagonal
NB.   P   - n×n-matrix, rows permutation of A
NB.   L   - m×k-matrix, lower triangular
NB.   U1  - k×n-matrix, unit upper triangular
NB.   k   = min(m,n)
NB.
NB. Assertion:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A -: clean p { U1 mp L
NB.   A -: clean p C. U1 mp L
NB.   A -: clean ip C.^:_1 U1 mp L
NB.   A -: clean P mp U1 mp L
NB. where
NB.   'ip U1L'=. getrfpu1l A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.   U1=. (tru1~ (-~/@$)) U1L
NB.   L=. (trl~ (-~/@$)) U1L

getrfpu1l=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. m) ; y
  elseif. 1 = n do.
    dip=. < _1 , liolmax y                                        NB. non-standard cycle permutation!
    ip=. dip (C. :: ]) i. m
    y=. ((] _1:} %) (_1&({,))) dip (C. :: ]) y                    NB. permute single column, scale by head, keep head unscaled
    ip ; y
  elseif. do.
    k=. m (<. >.@-:) n
    'pia Afa'=. getrfpu1l (-k) {."1 y                             NB. factorize 1st block recursively
    y=. pia C. (-k) }."1 y                                        NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trsmu1x & ((-k) & {.)) y                          NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. getrfpu1l y ((- (mp & Afba)) & ((-k) & }.)) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. pib , ((m-k) + i. k)                                   NB. apply 2nd block permutation to 1st block
    (dpib C. pia) ; ((Afbb , Afba) ,. (dpib C. Afa))              NB. assemble solution
  end.
)
  
NB. ---------------------------------------------------------
NB. getrful1p
NB.
NB. Description:
NB.   Triangular factorization with partial pivoting of a
NB.   general matrix:
NB.     U * L1 * P = A
NB.
NB. Syntax:
NB.     'ip UL1'=. getrful1p A
NB. where
NB.   A   - m×n-matrix to factorize
NB.   ip  - n-vector, columns inversed permutation of A
NB.   UL1 - m×n-matrix, upper triangle contains U, and strict
NB.         lower triangle contains L1 without unit diagonal
NB.   P   - n×n-matrix, columns permutation of A
NB.   L   - m×k-matrix, unit lower triangular
NB.   U1  - k×n-matrix, upper triangular
NB.   k   = min(m,n)
NB.
NB. Assertion:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A -: clean p {"1 U mp L1
NB.   A -: clean p C."1 U mp L1
NB.   A -: clean ip C.^:_1"1 U mp L1
NB.   A -: clean U mp L1 mp iP           NB. why iP ???
NB. where
NB.   'ip UL1'=. getrful1p A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.   U=. (tru~ (-~/@$)) UL1
NB.   L1=. (trl1~ (-~/@$)) UL1

getrful1p=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. n) ; y
  elseif. 1 = m do.
    dip=. < _1 , liolmax {. y                                         NB. non-standard cycle permutation!
    ip=. dip (C. :: ]) i. n
    y=. ((] _1:} %) (_1&({,))) dip (C."1 :: ]) y                      NB. permute single row, scale by head, keep head unscaled
    ip ; y
  elseif. do.
    k=. n (<. >.@-:) m
    'pia Afa'=. getrful1p (-k) {. y                                   NB. factorize 1st block recursively
    y=. pia (C."1) (-k) }. y                                          NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trsmxl1 & ((-k) & ({."1))) y                          NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. getrful1p y ((- (Afba & mp)) & ((-k) & (}."1))) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. pib , ((n-k) + i. k)                                       NB. apply 2nd block permutation to 1st block
    (dpib C."1 pia) ; ((Afbb ,. Afba) , (dpib C."1 Afa))              NB. assemble solution
  end.
)

NB. ---------------------------------------------------------
NB. hetrfpl
NB.
NB. Description:
NB.   Triangular factorization with full pivoting of a
NB.   Hermitian (symmetric) matrix:
NB.     P * L1 * T * L1^H * P^_1 = A
NB.   by blocked version of Aasen's algorithm
NB.
NB. Syntax:
NB.     'ip L1 T'=. hetrfpl A
NB. where
NB.   A   - n×n-matrix to factorize, Hermitian (symmetric)
NB.   ip  - n-vector, full inversed permutation of A
NB.   P   - n×n-matrix, full permutation of A
NB.   L1  - n×n-matrix, unit lower triangular
NB.   T   - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.
NB. Assertion:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A -: clean p sp T (] mp (mp ct)) L1
NB.   A -: clean P (mp mp |:@[) T (] mp (mp ct)) L1
NB. where
NB.   'ip L1 T'=. hetrfpl A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.
NB. References:
NB. [1] Miroslav Rozloznik, Gil Shklarski, Sivan Toledo.
NB.     Partitioned triangular tridiagonalization.
NB.     26 September 2007.
NB.     http://www.cs.cas.cz/miro/rst08.pdf
NB.
NB. TODO:
NB. - implement partitioned algorithm
NB. - T would be sparse

NB. ---------------------------------------------------------
NB. hetrfpu
NB.
NB. Description:
NB.   Triangular factorization with full pivoting of a
NB.   Hermitian (symmetric) matrix:
NB.     P * U1 * T * U1^H * P^_1 = A
NB.   by blocked version of Aasen's algorithm
NB.
NB. Syntax:
NB.     'ip U1 T'=. hetrfpu A
NB. where
NB.   A   - n×n-matrix to factorize, Hermitian (symmetric)
NB.   ip  - n-vector, full inversed permutation of A
NB.   P   - n×n-matrix, full permutation of A
NB.   U1  - n×n-matrix, unit upper triangular
NB.   T   - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.
NB. Assertion:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A -: clean p sp T (] mp (mp ct)) U1
NB.   A -: clean P (mp mp |:@[) T (] mp (mp ct)) U1
NB. where
NB.   'ip U1 T'=. hetrfpu A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.
NB. TODO:
NB. - implement partitioned algorithm
NB. - T would be sparse

hetrfpu=: hetf2pu

NB. ---------------------------------------------------------
NB. potrfl
NB.
NB. Description:
NB.   Cholesky factorization of a Hermitian (symmetric)
NB.   positive definite matrix:
NB.     L * L^H = A
NB.
NB. Syntax:
NB.   L=. potrfl A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.   L - n×n-matrix, lower triangular with positive diagonal
NB.       entries, Cholesky triangle
NB.
NB. Assertion:
NB.   A -: clean (mp ct) L
NB. where
NB.   L=. potrfl A

potrfl=: 3 : 0
  n=. # y
  if. n > 1 do.
    k=. >. -: n
    Ta=. potrfl (2 # k) {. y     NB. recursively factorize square matrix from top left or bottom right corner
    Ac=. (k ([ , -) n) {. y      NB. off-diagonal part of input matrix
    Tb=. ct Tbh=. Ta trsmlx Ac   NB. off-diagonal part of output matrix
    Aa=. (2 # k) }. y            NB. square matrix from opposite corner on diagonal of input matrix
    Tc=. potrfl Aa - Tb mp Tbh   NB. recursively factorize square matrix from opposite corner on diagonal
    Ta 0 append (Tb ,. Tc)       NB. assemble output as triangular matrix
  else.
    %: y
  end.
)

potrflf=: %:`((0:`0:`0:`]`]`(potrflf@(0 0 & {::))`,`(ct@])`trsmlx`(0 1 & {::)`[`mp`[`]`(1 1 & {::)`(_1 stitch)`(potrflf@:-~)`]`[`0:`0: fork6)@((<;.1)~ (;~@((0) 1:`(, >.@-:)`(#~)} #))))@.(1<#)

NB. ---------------------------------------------------------
NB. potrfu
NB.
NB. Description:
NB.   Cholesky factorization of a Hermitian (symmetric)
NB.   positive definite matrix:
NB.     U * U^H = A
NB.
NB. Syntax:
NB.   U=. potrfu A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.   U - n×n-matrix, upper triangular with positive diagonal
NB.       entries, Cholesky triangle
NB.
NB. Assertion:
NB.   A -: clean (mp ct) U
NB. where
NB.   U=. potrfu A

potrfu=: 3 : 0
  n=. # y
  if. n > 1 do.
    k=. >. -: n
    Ta=. potrfu (2 # k) }. y     NB. recursively factorize square matrix from top left or bottom right corner
    Ac=. (k ([ , -) n) }. y      NB. off-diagonal part of input matrix
    Tb=. ct Tbh=. Ta trsmux Ac   NB. off-diagonal part of output matrix
    Aa=. (2 # k) {. y            NB. square matrix from opposite corner on diagonal of input matrix
    Tc=. potrfu Aa - Tb mp Tbh   NB. recursively factorize square matrix from opposite corner on diagonal
    (Tc ,. Tb) _1 append Ta      NB. assemble output as triangular matrix
  else.
    %: y
  end.
)

NB. ---------------------------------------------------------
NB. pttrfl
NB.
NB. Description:
NB.   Triangular factorization of a Hermitian (symmetric)
NB.   positive definite tridiagonal matrix:
NB.     L1 * D * L1^H = A
NB.
NB. Syntax:
NB.   'L1 D'=. pttrfl A
NB. where
NB.   A  - n×n-matrix, Hermitian (symmetric) positive
NB.        definite tridiagonal
NB.   L1 - n×n-matrix, unit lower bidiangonal
NB.   D  - n×n-matrix, diagonal with positive diagonal
NB.        entries
NB.
NB. Algorithm:
NB.   In:  A
NB.   Out: L1 D
NB.   0) extract main diagonal d and suberdiagonal e from A
NB.   1) prepare input:
NB.        dee2 := d ,. (0,e) ,. (0,|e|^2)
NB.   2) do iterations k=1:n-1 by reversed suffix scan:
NB.        de=. u~/\.&.|. dee2
NB.      to find :
NB.        d[k] := d[k] - |e[k-1]|^2 / d[k-1]
NB.        e[k-1] := e[k-1] / d[k-1]
NB.   3) extract d - D's main diagonal, and e - L1's
NB.      subdiagonal from de:
NB.        d=. {."1 de
NB.        e=. }. 1 {"1 de
NB.   4) form output matrices:
NB.        L1=. (e;_1) setdiag idmat $ A
NB.        D=. diagmat d
NB.   5) link matrices L1 and D to form output:
NB.        L1 ; D
NB.
NB. Assertion:
NB.   A -: clean L1 (mp mp (ct@[)) D
NB. where
NB.   'L1 D'=. pttrfl A
NB.
NB. Notes:
NB. - implements LAPACK's xPTTRF
NB. - if A is indefinite then factors may have unacceptably
NB.   large elements
NB.
NB. References:
NB. [1] G. H. Golub and C. F. Van Loan, Matrix Computations,
NB.     Johns Hopkins University Press, Baltimore, Md, USA,
NB.     3rd edition, 1996, p. 157
NB.
NB. TODO:
NB. - L1 and D would be sparse

pttrfl=: ({."1 (((setdiag idmat@#)~ ;&_1) ; diagmat@[) }.@(1&({"1)))@(({.@] ((- {:) , {.@]) ((% {.)~ }.))~/\.&.|.)@(diag (,. (,. soris)@(0&,)) _1&diag)

NB. ---------------------------------------------------------
NB. pttrfu
NB.
NB. Description:
NB.   Triangular factorization of a Hermitian (symmetric)
NB.   positive definite tridiagonal matrix:
NB.     U1 * D * U1^H = A
NB.
NB. Syntax:
NB.   'U1 D'=. pttrfu A
NB. where
NB.   A  - n×n-matrix, Hermitian (symmetric) positive
NB.        definite tridiagonal
NB.   U1 - n×n-matrix, unit upper bidiangonal
NB.   D  - n×n-matrix, diagonal with positive diagonal
NB.        entries
NB.
NB. Algorithm:
NB.   In:  A
NB.   Out: U1 D
NB.   0) extract main diagonal d and superdiagonal e from A
NB.   1) prepare input:
NB.        dee2 := d ,. (e,0) ,. ((|e|^2),0)
NB.   2) do iterations k=n-2:0 by suffix scan:
NB.        de=. u~/\. dee2
NB.      to find :
NB.        d[k] := d[k] - |e[k]|^2 / d[k+1]
NB.        e[k] := e[k] / d[k+1]
NB.   3) extract d - D's main diagonal, and e - L1's
NB.      subdiagonal from de:
NB.        d=. {."1 de
NB.        e=. }: 1 {"1 de
NB.   4) form output matrices:
NB.        U1=. (e;1) setdiag idmat $ A
NB.        D=. diagmat d
NB.   5) link matrices U1 and D to form output:
NB.        U1 ; D
NB.
NB. Assertion:
NB.   A -: clean U1 (mp mp (ct@[)) D
NB. where
NB.   'U1 D'=. pttrfu A
NB.
NB. Notes:
NB. - if A is indefinite then factors may have unacceptably
NB.   large elements
NB.
NB. TODO:
NB. - U1 and D would be sparse

pttrfu=: ({."1 (((setdiag idmat@#)~ ;&1) ; diagmat@[) }:@(1&({"1)))@(({.@[ ((- {:) , {.@]) ((%~ }.)~ {.))/\.)@(diag (,. (,. soris)@(,&0)) 1&diag)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgetrf
NB.
NB. Description:
NB.   Test:
NB.   - lud (math/misc)
NB.   - getrf (math/lapack)
NB.   - getrfxxxx (math/mt)
NB.   by general matrix given
NB.
NB. Syntax:
NB.   testgetrf A
NB. where
NB.   A - m×n-matrix
NB.
NB. Formula:
NB. - for L * U1 * P = A :
NB.     berr := ||L * U1 * P - A|| / (ε * ||A|| * m)
NB. - for P * L1 * U = A :
NB.     berr := ||P * L1 * U - A|| / (ε * ||A|| * n)
NB. - for P * U1 * L = A :
NB.     berr := ||P * U1 * L - A|| / (ε * ||A|| * n)
NB. - for U * L1 * P = A :
NB.     berr := ||U * L1 * P - A|| / (ε * ||A|| * m)

testgetrf=: 3 : 0
EMPTY
return.
  require '~addons/math/misc/matfacto.ijs'
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'getrf'

  rcond=. ((_."_)`(norm1 con (getriul1p@getrful1p)) @. (=/@$)) y  NB. meaninigful for square matrices only

  ('lud' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (((mp & >)/ @ }:) %. (2 & {::) ))) % (FP_EPS*((norm1*c)@[)))))) y

  ('getrf_jlapack_' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (((mp & >)/ @ }:) invperm_jlapack_ (2 & {::) ))) % (FP_EPS*((norm1*c)@[)))))) y

  ('getrflu1p'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((0 & {::) C.^:_1"1 (( trl            mp  tru1          )@(1 & {::))))) % (FP_EPS*((norm1*#)@[)))))) y
  ('getrfpl1u'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((0 & {::) C.^:_1   (( trl1           mp  tru           )@(1 & {::))))) % (FP_EPS*((norm1*c)@[)))))) y
  ('getrfpu1l'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((0 & {::) C.^:_1   (((tru1~ (-~/@$)) mp (trl ~ (-~/@$)))@(1 & {::))))) % (FP_EPS*((norm1*c)@[)))))) y
  ('getrful1p'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((0 & {::) C.^:_1"1 (((tru ~ (-~/@$)) mp (trl1~ (-~/@$)))@(1 & {::))))) % (FP_EPS*((norm1*#)@[)))))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhetrf
NB.
NB. Description:
NB.   Test hetrfpx by Hermitian (symmetric) matrix given
NB.
NB. Syntax:
NB.   testhetrf A
NB. where
NB.   A - n×n-matrix, Hermitian
NB.
NB. Formula:
NB. - for P * L * T * L^H * P^_1 = A :
NB.     berr := ||P * L * T * L^H * P^_1 - A|| / (ε * ||A|| * n)
NB. - for P * U * T * U^H * P^_1 = A :
NB.     berr := ||P * U * T * U^H * P^_1 - A|| / (ε * ||A|| * n)

testhetrf=: 3 : 0
EMPTY
return.
  rcond=. (norm1 con (hetripl@hetrfpl)) y

  ('hetrfpl' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((mp mp ct@[)&>/@}. (sp~ /:) (0&({::)))))) % (FP_EPS*((norm1*c)@[))))) y
  ('hetrfpu' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((mp mp ct@[)&>/@}. (sp~ /:) (0&({::)))))) % (FP_EPS*((norm1*#)@[))))) y

EMPTY
)

NB. ---------------------------------------------------------
NB. testpotrf
NB.
NB. Description:
NB.   Test:
NB.   - choleski (math/misc addon)
NB.   - potrf (math/lapack addon)
NB.   - potrfx (math/mt addon)
NB.   by Hermitian (symmetric) positive definite matrix given
NB.
NB. Syntax:
NB.   testpotrf A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive
NB.       definite
NB.
NB. Formula:
NB. - for L * L^H = A :
NB.     berr := ||L * L^H - A|| / (ε * ||A|| * n)
NB. - for U * U^H = A :
NB.     berr := ||U * U^H - A|| / (ε * ||A|| * n)

testpotrf=: 3 : 0
  require '~addons/math/misc/matfacto.ijs'
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'potrf'

  rcond=. (norm1 con (potril@potrfl)) y

  ('choleski_mt_'       tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (mp ct)))) % (FP_EPS*((norm1*c)@[))))) y

  ('potrf_jlapack_' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (mp ct)))) % (FP_EPS*((norm1*c)@[))))) y

  ('potrfl'         tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (mp ct)))) % (FP_EPS*((norm1*c)@[))))) y
  ('potrflf'        tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (mp ct)))) % (FP_EPS*((norm1*c)@[))))) y
  ('potrfu'         tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (mp ct)))) % (FP_EPS*((norm1*#)@[))))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpttrf
NB.
NB. Description:
NB.   Test pttrfpx by Hermitian (symmetric) positive definite
NB.   tridiagonal matrix given
NB.
NB. Syntax:
NB.   testpttrf A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive
NB.       definite tridiagonal
NB.
NB. Formula:
NB. - for L1 * D * L1^H = A :
NB.     berr := ||L1 * D * L1^H - A|| / (ε * ||A|| * n)
NB. - for U1 * D * U1^H = A :
NB.     berr := ||U1 * D * U1^H - A|| / (ε * ||A|| * n)

testpttrf=: 3 : 0
  rcond=. (norm1 con (pttril@pttrfl)) y

  ('pttrfl' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((mp mp (ct@[))&>/)))) % (FP_EPS*((norm1*c)@[))))) y
  ('pttrfu' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((mp mp (ct@[))&>/)))) % (FP_EPS*((norm1*#)@[))))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testtrf
NB.
NB. Description:
NB.   Adv. to make verb to test xxtrfxxxx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testtrf
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     (? @ $ 0:) testtrf_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testtrf_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testtrf_mt_ 150 200

testtrf=: 1 : 'EMPTY_mt_ [ (((testpttrf_mt_ @ (u ptmat_mt_)) [ (testpotrf_mt_ @ (u pomat_mt_)) [ (testhetrf_mt_ @ (u hemat_mt_))) ^: (=/)) [ testgetrf_mt_ @ u'


NB. ip=. 0 4 3 1 2
NB. L1=. 5 5 $ 1 0 0 0 0 0 1 0 0 0 0 _0.7 1 0 0 0 0.1 0.8 1 0 0 _0.1 0.2 0.6 1
NB. T=. 5 5 $ 1 0.3 0 0 0 0.3 1.6 0.7 0 0 0 0.7 0.6 0.9 0 0 0 0.9 1.4 0.3 0 0 0 0.3 1.2
NB. HE5=. 5 5 $ 1 0.03 _0.03 _0.21 0.3 0.03 3.352 1.79 0.946 0.72 _0.03 1.79 2.292 0.604 _0.02 _0.21 0.946 0.604 0.404 _0.42 0.3 0.72 _0.02 _0.42 1.6

NB. L1=. trl1_mt_ 0 (< a: ; 0) } 0.1 * _8 + ? 5 5 $ 18
NB. T=. (+ ct_mt_) bdlpick_mt_ trl_mt_ 0.1 * >: ? 5 5 $ 9
NB. ip=. (?@!@<: A. i.) 5
NB. HE5=. clean (/: ip) sp_mt_ L1 (mp mp (ct_mt_ @ [)) T

NB. L1=. trl1_mt_ 0 (< a: ; 0) } 0.1 * _8 + ? 6 6 $ 18
NB. T=. (+ ct_mt_) bdlpick_mt_ trl_mt_ 0.1 * >: ? 6 6 $ 9
NB. ip=. (?@!@<: A. i.) 6                                  NB. <: to force ip[0]==0
NB. HE6=. clean (/: ip) sp_mt_ L1 (mp mp (ct_mt_ @ [)) T

NB. L1=. trl1_mt_ 0 (< a: ; 0) } j./ 0.1 * _7 + ? 2 10 10 $ 14
NB. T=. (+ ct_mt_) bdlpick_mt_ trl_mt_ j./ 0.1 * >: ? 2 10 10 $ 9
NB. ip=. (?@!@<: A. i.) 10
NB. HE10=. clean (/: ip) sp_mt_ L1 (mp mp (ct_mt_ @ [)) T
NB. ip ; L1 ; (i.0) ; (i.0) ; (diag_mt_ T) ; (1 diag_mt_ T) ; (ip sp_mt_ HE10) ; (L1 mp T) ; (i.0) ; (i.0)
NB. clean lahefpl_mt_ ((i.10);(0 {."1 HE10);(10 ($!.0) 1);(,0);(0 ({,) HE10);(i.0);HE10;(0 {."1 HE10);({."1 HE10);((< (<0);0) { HE10))
NB. clean hetrfpl_mt_ HE10
