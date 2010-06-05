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

NB. NB. gi        Conj. to evoke m-th verb from gerund n
NB. gi=: 2 : '(m{n)`:6'                    NB. Conj. to evoke m-th verb from gerund n: n[m]
NB. 
NB. tile=: 1 : '[`((0 gi m)~)`((1 gi m)&#)`(2 gi m)`(0 gi m)`] fork3'
NB. 
NB. stitcht=: {.  `   >. `,. tile
NB. stitchb=: {.  `(-@>.)`,. tile
NB. appendl=: {."1`   >. `,  tile
NB. appendr=: {."1`(-@>.)`,  tile

NB. stitcht=: ,.`([,.(({.~    # )~))`(({.~    # ),.]) @. (*@-&#)
NB. stitchb=: ,.`([,.(({.~ (-@#))~))`(({.~ (-@#)),.]) @. (*@-&#)

stitchrb=: [ ,. ({.~ (-@#))~            NB. stitch right aligned down to the left

NB. ---------------------------------------------------------
NB. xsh
NB.
NB. Description:
NB.   Extend [dimensions within] shape
NB.
NB. Syntax:
NB.   eA=. dsh xsh A
NB. where
NB.   A   - sh-array to extend by zeros
NB.   dsh - r-vector or scalar, integer, delta of shape,
NB.         negative element forces corresponding dimension
NB.         extending from the head
NB.   eA  - (|dsh|+sh)-array, being A extended by zeros
NB.   sh  - r-vector, the shape of A
NB.   r   ≥ 0, the rank of A and eA
NB.
NB. Examples:
NB.    0 xsh 2 2 $ 2          _1 xsh 2 2 $ 2
NB. 2 2                    0 0 0
NB. 2 2                    0 2 2
NB.                        0 2 2
NB.    1 3 xsh 2 2 $ 2        1 _3 xsh 2 2 $ 2
NB. 2 2 0 0 0              0 0 0 2 2
NB. 2 2 0 0 0              0 0 0 2 2
NB. 0 0 0 0 0              0 0 0 0 0
NB.    _1 3 xsh 2 2 $ 2       _1 _3 xsh 2 2 $ 2
NB. 0 0 0 0 0              0 0 0 0 0
NB. 2 2 0 0 0              0 0 0 2 2
NB. 2 2 0 0 0              0 0 0 2 2

xsh=. [: : (((condneg"0 + [) $) {. ])

jomat=: (((, 0:)~ #) xsh [) stitchrb ]  NB. make jordain matrix from blocks x and y

NB. 'ip L1 T permHE5 H'=. 0 0 0 0 hetf2pl HE5
NB. 'ip L1 T permA H'=. l0 hetf2pl A
NB. 'ip L1 T permA H h'=. step (ip;L1;T;A;H;h)

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

NB. lahefpl_mt_ ((i.5);(0 {."1 HE5);(5 ($!.0) 1);(4 # 0);(0 ({,) HE5);(i.0);HE5;(1 {."1 HE5);(4 # 0);((< (<0);0) { HE5))
NB. lahefpl_mt_ ((i.10);(0 {."1 HEci10);(10 ($!.0) 1);(9 # 0);(0 ({,) HEci10);(i.0);HEci10;(1 {."1 HEci10);(9 # 0);((< (<0);0) { HEci10))

lahefpl=: (3 : 0) ^: (TRFNB<.(#@(3 & {::)))
  'ip L1 l0 l1 t0 t1 A H h0 h1'=. y
  L1=. L1 stitchrb l0
  l0=. }. l0
  'n j'=. $ L1
  l1=. h1 - l0 * ({: t0)
  q=. liofmax l1
  dip0=. 0 lios2cp q
  dip=. j (+ &. >) dip0
  l0=. dip0 C. l0
  l1=. dip0 C. l1
  ip=. dip C. ip
  L1=. dip C. L1
  H=. dip C. H
  A=. dip sp A
  t01=. + t10=. {. l1
  l1=. l1 % t10
  h0=. ((j (] dhs2lios (((* >:)),-~)) n) ({,) A) - (j }. H) mp (+ j { L1)
  H=. H stitchrb (t01 , h0)
  h1=. h0 - l0 * t01
  t11=. {. h1
  t0=. t0 , t11
  t1=. t1 , t01
  ip ; L1 ; l1 ; l0 ; t0 ; t1 ; A ; H ; (}. h0) ; (}. h1)
)

NB. 'ip L1 T'=. hetrfpl A
NB. 'ip L1 T permA H'=. step (ip;L1;T;A;H)

NB. 'ipi L1i Ti Ai Hi hi'=. (9 # 0) hetf2pl_mt_ HEci10
NB. ((((] dhs2lios (_1,-))/@$) ({,) ]) L1i) hetf2pl_mt_ (((_1 ({,) Ti) , hi) (< a: ; 0)} (2 2 }. HEci10))

hetrfpl=: 3 : 0
  n=. # y
  ip=. i. n
  L1=. 0 {."1 y
  t0=. 0 ({,) y
  t1=. i. 0
  h1i=. (< (<0);0) { y
  l0i=. n ($!.0) 1
  l1i=. h0i=. }. l0i                             NB. (<: n) # 0
  for_k. n (] dhs2lios (0,(>.@%))) TRFNB do.
    'ipi L1i l0i l1i t0i t1i y Hi h0i h1i'=. (lahefpl dbg 'lahefpl') ((i.#y);(0 {."1 y);l0i;l1i;({:t0);(i.0);y;(,.({:t0),h1i);h0i;h1i)
    dip=. (i.k),(k+ipi)
    ip=. dip (C. dbg 'C.ip') ip
    L1=. (dip (C. dbg 'C.L1') L1) (stitchrb dbg 'L1 stitchb L1i') L1i
    t0=. t0 (, dbg 't0') (}. t0i)
    t1=. t1 (, dbg 't1') t1i
    y=. ((2 # TRFNB) }. y) (- dbg 'Aupd') ((Hi ((((mp ct) dbg 'mpct')&(TRFNB&}.)) dbg 'Hi*L1i') L1i) (+ dbg '+') ((l0i (* dbg '*') {: t1i) (*/ dbg '*/') l1i))
  end.
  T=. ((+t1);_1) setdiag (t1;1) setdiag (t0;0) (2 # n) $ 0
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
    v=. (j lios2cp q) C. v         NB. v[_1]↔v[q]
    dip=. q lios2cp j              NB. any[j]↔any[q]
    y=. dip sp y                   NB. A[j,0:j]↔A[q,0:j], A[0:j,j]↔A[0:j,q]
    ip=. dip C. ip                 NB. ip[j]↔ip[q]
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
    dip=. 0 lios2cp liofmax {. y
    ip=. dip C. i. n
    y=. ((] 0:} %) (0&({,))) dip (C."1) y                          NB. permute single row, scale by head, keep head unscaled
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
    dip=. 0 lios2cp liofmax y
    ip=. dip C. i. m
    y=. ((] 0:} %) (0&({,))) dip C. y                          NB. permute single column, scale by head, keep head unscaled
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
    dip=. (<:m) lios2cp liolmax y
    ip=. dip C. i. m
    y=. ((] _1:} %) (_1&({,))) dip C. y                           NB. permute single column, scale by head, keep head unscaled
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
    dip=. (<:n) lios2cp liolmax {. y
    ip=. dip C. i. n
    y=. ((] _1:} %) (_1&({,))) dip C."1 y                             NB. permute single row, scale by head, keep head unscaled
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



NB.    ] L1=. 1 0 0 0 0 (< a: ; 0) } trl1_mt_ 0.1 * j./ _8 + ? 2 5 5 $ 18
NB. 1       0        0        0 0
NB. 0       1        0        0 0
NB. 0  0j_0.2        1        0 0
NB. 0     0.9 _0.5j0.7        1 0
NB. 0 0.9j0.1 _0.4j0.8 _0.4j0.2 1
NB.    | L1
NB. 1        0        0        0 0
NB. 0        1        0        0 0
NB. 0      0.2        1        0 0
NB. 0      0.9 0.860233        1 0
NB. 0 0.905539 0.894427 0.447214 1
NB.    ] T=. (+ ct_mt_) bdlpick_mt_ trl_mt_ 0.1 * j./ >: ? 2 5 5 $ 9
NB.     1.6 0.5j_0.9        0        0        0
NB. 0.5j0.9      0.2 0.6j_0.5        0        0
NB.       0  0.6j0.5      1.6 0.7j_0.6        0
NB.       0        0  0.7j0.6      0.6 0.3j_0.8
NB.       0        0        0  0.3j0.8      0.8
NB.    ip=. 0 1 3 4 2
NB.    P=. ip2P_mt_ ip
NB.    ] HE5=. clean P mp L1 mp T mp (ct_mt_ L1) mp |: P
NB.       1.6   0.5j_0.9   0.36j_0.86    0.18j0.1   0.45j_0.81
NB.   0.5j0.9        0.2   _0.46j_0.3   0.6j_0.46  _0.47j_0.17
NB. 0.36j0.86  _0.46j0.3        1.508 _0.51j0.698   0.624j1.91
NB. 0.18j_0.1   0.6j0.46 _0.51j_0.698       1.408 0.406j_1.176
NB. 0.45j0.81 _0.47j0.17  0.624j_1.91 0.406j1.176        0.916
NB.    hetf2plo_mt_ HE5
NB. ┌─────────┬─────────────────────────────┬───────────────────────────────────────────┐
NB. │0 1 3 4 2│1       0        0        0 0│    1.6 0.5j_0.9        0        0        0│
NB. │         │0       1        0        0 0│0.5j0.9      0.2 0.6j_0.5        0        0│
NB. │         │0  0j_0.2        1        0 0│      0  0.6j0.5      1.6 0.7j_0.6        0│
NB. │         │0     0.9 _0.5j0.7        1 0│      0        0  0.7j0.6      0.6 0.3j_0.8│
NB. │         │0 0.9j0.1 _0.4j0.8 _0.4j0.2 1│      0        0        0  0.3j0.8      0.8│
NB. └─────────┴─────────────────────────────┴───────────────────────────────────────────┘
NB.    1 0 0 0 0 hetf2plb_mt_ HE5
