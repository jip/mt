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

TRFNB=: 64  NB. block size limit

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


stitchrb=: [ ,. ({.~ (-@#))~

NB. 'ip L1 T permA H'=. l0 hetf2pl A
NB. 'ip L1 T permA H h'=. step (ip;L1;T;A;H;h)

hetf2pla=: 5 {. ((3 : 0) ^: (0>.(<: TRFNB)<.(<:@#@(3 & {::)))) @ ((i.@#@]);(,.@[);(1 1 {. ]);];(1 {."1 ]);(}.@({."1 @ ])))
  'ip L1 T A H h'=. y
  'n j'=. $ L1
  l0=. ((j dhs2lios (_1&,)) (n-j)) ({,) L1
  l1=. h - l0 * (_1 ({,) T)
  q=. liofmax l1
  dip0=. 0 lios2cp q
  dip=. j (+ &. >) dip0              NB. j ([ lios2cp +) q
  l1=. dip0 C. l1
  A=. dip sp A
  ip=. dip C. ip
  l0=. dip0 C. l0
  H=. dip C. H
  L1=. dip C. L1
  t01=. + t10=. {. l1
  l1=. l1 % t10
  h=. ((< (n th2lios j) ; j) { A) - (j }. H) mp (+ j { L1)
  L1=. L1 stitchrb l1
  H=. H stitchrb (t01 , h)
  h=. h - l0 * t01
  t11=. {. h
  T=. (t11,t10,t01) (_1 _2,(_1 - #@])) } 1 xsh T
  ip ; L1 ; T ; A ; H ; (}. h)
)

xsh=: 1 : '{.~ (m+$)'

NB. 'rawip L1 T permsubA H l0'=. l0 hetf2pl A
NB. 'rawip L100 T permsubA H00 l0 h H10 L110'=. step (rawip;L100;T;A;H00;l0;h;H10;L110)

hetf2plb=: 6 {. ((3 : 0) ^: (0>.(<: TRFNB)<.(#@(3 & {::)))) @ (0:;(1 1 $ [);(1 1 {. ]);(1 1 }. ]);(1 1 {. ]);(}.@[);(}.@:({."1)@]);(}.@(1 {."1 ]));(,.@}.@[))
  'rawip L100 T A H00 l0 h H10 L110'=. y
  l1=. h - l0 * (_1 ({,) T)
  q=. liofmax l1
  dip=. 0 lios2cp q
  l1=. dip C. l1
  A=. dip sp A
  rawip=. rawip , q
  l0=. dip C. l0
  H10=. dip C. H10
  L110=. dip C. L110
  t01=. + t10=. {. l1
  l1=. l1 % t10
  h=. ({."1 A) - H10 mp (+ {. L110)
  H10=. H10 ,. h
  H00=. (t01 _1: } 0 1 xsh H00) , ({. H10)
  H10=. }. H10
  L110=. L110 ,. l1
  L100=. ({. L110) _1 } 1 xsh L100
  L110=. }. L110
  h=. h - l0 * t01
  t11=. {. h
  T=. (t11,t10,t01) (_1 _2,(_1 - #@])) } 1 xsh T
  rawip ; L100 ; T ; (1 1 }. A) ; H00 ; (}. l1) ; (}. h) ; H10 ; L110
)

hetf2plo=: 3 : 0
  n=. # y
  L1=. ($ y) {. 1              NB. 1st L1's column in Aasen method is e(1)-vector
  T=. ($ {. (1 1 {. ])) y      NB. T[0,0]=A[0,0]
  ip=. i. n
  h=. i. 1                     NB. h[0:j-1]=H[0:j-1,j-1], may be defined arbitrary before the 1st iteration only
  ios=. n th2lios 1            NB. j:n-1
  for_j. (n <. TRFNB) th2lios 1 do.       NB. 1:min(n,NB)-1
    a=. (< ios ; (<: j)) { y   NB. A[j:n-1,j-1]
    lum=. (j (- , [) n) {. L1  NB. L1[j:n-1,0:j-1]
    v=. a - lum mp h           NB. v[0:n-j-1]=A[j:n-1,j-1]-L1[j:n-1,0:j-1]*H[0:j-1,j-1]=L1[j:n-1,j]*T[j,j-1] non-pivoted yet
    q=. liofmax v              NB. IO pivot from head
    v=. (0 lios2cp q) C. v     NB. v[0]↔v[q]
    dip=. q (+ lios2cp ]) j    NB. any[j]↔any[j+q]
    y=. dip sp y               NB. A[j,j:n-1]↔A[j+q,j:n-1], A[j:n-1,j]↔A[j:n-1,j+q]
    ip=. dip C. ip             NB. ip[j]↔ip[j+q]
    to=. {. v                  NB. T[j,j-1]
    lu=. q { lum               NB. L1[j,0:j-1] after pivoting
    luecto=. to * + {: lu      NB. conj(L1[j,j-1])*T[j,j-1]
    h=. to (((+ +)~ {:)`_1:`]) } ((2 # j) {. T) mp + lu   NB. h[0:j-1]=H[0:j-1,j]=T[0:j-1,0:j-1]*conj(L1[j,0:j-1])
    L1=. (1 (0}) v % to) (< ios ; j) } dip C. L1          NB. L1[j,0:n-1]↔L1[j+q,0:n-1], L1[j:n-1,j]=v[0:n-j-1]/v[0], v[0] may be 0
    td=. 9 o. ((< 2 # j) { y) - (luecto + lu mp h)        NB. T[j,j]=Re(A[j,j]-L1[j,0:j-1]*H[0:j-1,j]-conj(L1[j,j-1])*T[j,j-1])
    T=. (td (,~ (, +)) to) (_2 <\ 0 _1 _1 0 0 0 + j) } T  NB. TODO: amend by lIOS NB. batch write diagonal and off-diagonal elements T[j,j-1] T[j-1,j] T[j,j]
    h=. h , luecto + td        NB. h[0:j]=H[0:j,j]=T[0:J,0:J]*conj(L1[j,0:j])
    ios=. }. ios               NB. j+1:n-1
  end.
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

hetrfpl=: 3 : 0
  n=. # y
  L1=. idmat n
  T=. ($ y) $ 0
  ip=. i. n
  l1=. n {. 1
  for_i. TRFNB dhs2lios (0, (<. n % TRFNB)) do.
    n=. # y
    nb=. TRFNB <. # y
    'ipi L1i Ti'=. (l1 ; nb) hetf2pl y
    T=. ((0 0,nb) (diag ; (i 1} [)) Ti) setdiag T
    k=. nb <. (n-i)
    d=. (1 0,k) diag Ti
    T=. (d;(1,i,k))) setdiag T
    T=. ((+d);(_1,i,k))) setdiag T           NB. TODO: amend by lIOS
    L1=. L1i (< (n th2lios i) ; (dhs2lios (i,nb))) } L1
    l1=. ((n-i) dhs2lios (_1,(n-(i+nb)))) ({,) L1i
    ipi=. (i. i) , (i+ipi)
    y=. ipi sp y
    L1=. ipi C. L1
    ip=. ipi C. ip
    A11=. A11 - H*L1^H - L1*T*L1^H       NB. FIXME!
    recombine A = (  A00  A01  )
                  (  A10  A11  )
  end.
  ip ; L1 ; T
)

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
