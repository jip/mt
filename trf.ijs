NB. trf.ijs
NB. Triangular factorization of matrix
NB.
NB. getrfpl1u  Factorization with partial pivoting of a
NB.            general matrix: inv(P) * L1 * U = A, where P
NB.            is permutation matrix, L1 is unit lower
NB.            triangular and U is upper triangular
NB. getrflu1p  Factorization with partial pivoting of a
NB.            general matrix: L * U1 * inv(P) = A, where P
NB.            is permutation matrix, L is lower triangular
NB.            and U1 is unit upper triangular
NB. getrfpu1l  Factorization with partial pivoting of a
NB.            general matrix: inv(P) * U1 * L = A, where P
NB.            is permutation matrix, U1 is unit upper
NB.            triangular and L is lower triangular
NB. getrful1p  Factorization with partial pivoting of a
NB.            general matrix: U * L1 * inv(P) = A, where P
NB.            is permutation matrix, U is upper triangular
NB.            and L1 is unit lower triangular
NB.
NB. getrfl1u   Factorization without pivoting of a general
NB.            matrix: L1 * U = A, where L1 is unit lower
NB.            triangular and U is upper triangular
NB. getrflu1   Factorization without pivoting of a general
NB.            matrix: L * U1 = A, where L is lower
NB.            triangular and U1 is unit upper triangular
NB. getrfu1l   Factorization without pivoting of a general
NB.            matrix: U1 * L = A, where U1 is unit upper
NB.            triangular and L is lower triangular
NB. getrful1   Factorization without pivoting of a general
NB.            matrix: U * L1 = A, where U is upper
NB.            triangular and L1 is unit lower triangular
NB.
NB. hetrfpu    Factorization with full pivoting of a
NB.            Hermitian (symmetric) matrix:
NB.            inv(P) * U1 * T * U1' * P = A, where P is
NB.            permutation matrix, U1 is unit upper
NB.            triangular and T is tridiagonal
NB. hetrfpl    Factorization with full pivoting of a
NB.            Hermitian (symmetric) matrix:
NB.            inv(P) * L1 * T * L1' * P = A, where P is
NB.            permutation matrix, L1 is unit lower
NB.            triangular and T is tridiagonal
NB.
NB. hetrfu     Factorization without pivoting of a Hermitian
NB.            (symmetric) matrix: U1 * T * U1' = A, where P
NB.            is permutation matrix, U1 is unit upper
NB.            triangular and T is tridiagonal
NB. hetrfl     Factorization without pivoting of a Hermitian
NB.            (symmetric) matrix: L1 * T * L1' = A, where P
NB.            is permutation matrix, L1 is unit lower
NB.            triangular and T is tridiagonal
NB.
NB. potrfu     Cholesky factorization of a Hermitian
NB.            (symmetric) positive definite matrix:
NB.            U * U' = A, where U is upper triangular
NB. potrfl     Cholesky factorization of a Hermitian
NB.            (symmetric) positive definite matrix:
NB.            L * L' = A, where L is lower triangular
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-10

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

TRFNB=: 32   NB. block size limit

NB. ---------------------------------------------------------

iofmaxm=: (i.>./) @: |  NB. IO 1st element with max magnitude from list y
iolmaxm=: (i:>./) @: |  NB. IO last element with max magnitude from list y

NB. ---------------------------------------------------------
NB. Description:
NB.   Single step of non-blocked version of algorithms
NB.
NB. Syntax
NB.   'pi1 pfxi1 sfxLi1 sfxRi1=. getf2xxxxstep (pi ; pfxi ; sfxLi ; sfxRi)
NB. where
NB.   pi     - p(i), m-vector after i-th and before (i+1)-th
NB.            step, permutation
NB.   pfxi   - pfx(i), i*n-matrix, first i rows of matrix
NB.            L1U(i)
NB.   sfxLi  - sfxL(i), the left part of 
NB.   sfxRi  - sfxR(i), the right part of 
NB.   pi1    - p(i+1), m-vector after (i+1)-th step,
NB.            permutation
NB.   pfxi1  - pfx(i+1), (i+1)*n-matrix, first i rows of
NB.            matrix L1U(i+1)
NB.   sfxLi1 - sfxL(i+1), the left part of 
NB.   sfxRi1 - sfxR(i+1), the right part of 

getf2pl1ustepa=: 3 : 0
  'p pfx sfxL sfxR'=. y
  c=. {."1 sfxR
  jp=. iofmaxm c
  r=. jp { sfxR
  c=. (}. ({. c) jp } c) % ({. r)
  ((jp (+ ii2cp ]) # pfx) C. p) ; (pfx , ((jp { sfxL) , r)) ; ((}. ({. sfxL) jp } sfxL) ,. c) ; ((1 1 }. ({. sfxR) jp } sfxR) - (c */ }. r))
)

NB. ---------------------------------------------------------
NB. Description:
NB.   Non-blocked version of algorithms
NB.
NB. Syntax
NB.   'p L1U'=. getf2pl1u A
NB. where
NB.   A   - m*n-matrix

getf2pl1ua=: ((0&{) , ((,&.>)/@(1 2&{))) @ (getf2pl1ustepa ^: ((<./@$)`((i.@#);(0&{.);(_ 0&{.);])))

NB. ---------------------------------------------------------
NB. Description:
NB.   Single step of algorithms
NB.
NB. Syntax
NB.   'pi1 pfxi1 sfxLi1 sfxRi1'=. getf2xxxxstep (pi ; pfxi ; sfxLi ; sfxRi)
NB. where
NB.   pfxi   - pfx(i), ##########
NB.
NB. a00 a00 a00 a00      a00 a00 a00 a00
NB. a10 a11 a12 a12      a00 a00 a00 a00
NB. a10 a21 a22 a22      a10 a10 a11 a12
NB. a10 a21 a22 a22      a10 a10 a21 a22
NB.

getrfpl1ustepa=: 3 : 0
  'p pfx sfxL sfxR'=. y
  nj=. -~/ 'j n'=. $ pfx
  bs=. TRFNB <. (# sfxL) <. nj
  'pi L1U'=. getf2pl1ua bs {."1 sfxR
  sfxL=. pi C. sfxL
  A1222=. pi C. bs }."1 sfxR
  A11=. bs {. L1U
  A21=. bs }. L1U
  A12=. A11 trtrsl1x (bs {. A1222)
  p=. ((i.j),(j+pi)) C. p
  pfx=. pfx , ((bs {. sfxL) ,. A11 ,. A12)
  sfxL=. (bs }. sfxL) ,. A21
  sfxR=. (bs }. A1222) - A21 mp A12
  p ; pfx ; sfxL ; sfxR
)

NB. (-: (clean@((/: @ (0&({::))) { ((trl1 mp tru)@(1&({::))))@getrfpl1ua)) A
NB. 'ip L1U'=. getrfpl1ua A

getrfpl1ua=: getf2pl1ua`(((0&{),((,&.>)/@(1 2&{))) @ (getrfpl1ustepa ^: ((>.@(TRFNB %~ (<./@$)))`((i.@#);(0&{.);(_ 0 {. ]);]))))@.(TRFNB<(<./@$))

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

getf2pl1ustepb=: (((1;0)&({::)) (([ C.~ (<@~.@(],]+({~(<@(_1&,)))))))~ ((<@}.) ((([ (0}"0 1) (] - ((0:`(_1,({:@[))`[}) */ (({~ {:)~))))~ ((({`[`([ _1} (]%{))})~ iofmaxm)@:({."1))) upd1) (0&({::)))) ; (     }.&.>@}.)

getf2pl1ub=: 0 {:: (getf2pl1ustepb ^: ((_1 0&(ms $))`(];((;&i.)/@$))))


NB. Aperm=. invpraw laswp A
laswp=: (C.^:_1~ <@~."1@((+ ,. ]) i.@#))~

NB. A=. invpraw swp2fw Aperm
swp2fw=: (C.~ <@~."1@((+ ,. ]) i.@#))~

NB.        0     1        2     3   4          5          6
NB. input: 0:j-1 j:j+bs-1 j:m-1 j:m j+bs-1:m-1 j+bs-1:n-1 _1
NB.
NB. ioseY: 3;1
NB. iospi: 6;1
NB. iosL: 2;0
NB. iosC: 2;5
NB. iosA11: 1;1
NB. iosA12: 1;5
NB. iosA21: 4;1
NB. iosA22: 4;5
NB.
NB.         0     1     2    3    4      5      6      7
NB. output: ioseY iospi iosL iosC iosA11 iosA12 iosA21 iosA22
NB.
NB. 'ioseY iospi iosL iosC iosA11 iosA12 iosA21 iosA22'=. ios getrfpl1umkiosb A

getrfpl1umkiosb=: <"1@((4 2,7 2,3 1,3 6,2 2,2 6,5 2,:5 6)&{)

NB.         0         1               2         3       4              5              6
NB. input:  0:j-1     j:j+NB-1        j:m-1     j:m     j+NB-1:m-1     j+NB-1:n-1     _1
NB.
NB.         0         1               2         3       4              5              6
NB. output: 0:j+NB-1  j+NB:j+NB+bs-1  j+NB:m-1  j+NB:m  j+NB+nb-1:m-1  j+NB+nb-1:n-1  _1
NB.
NB. 'i0 i1 i2 i3 i4 i5 i6'=. getrfpl1uchgiosb (i0;i1;i2;i3;i4;i5;i6)

getrfpl1uchgiosb=: 3 : 0
  'i0 i1 i2 i3 i4 i5 i6'=. }. y
  nb=. TRFNB <. (# i4) <. (# i5)
  (i0 , i1) ; (nb rt i4) ; (TRFNB }. i2) ; (TRFNB }. i3) ; (TRFNB }. i4) ; (TRFNB }. i5) ; i6
)

getrfpl1ustepb=: (getrfpl1umkiosb ([ ([ ([ ([ (-`mp map3i 7 6 5 7) (trtrsl1x map2i 4 5 5)) (laswp map2i 1 3 3)) (laswp map2i 1 2 2)) (getf2pl1ub upd1i 0)) (0&({::))) ; getrfpl1uchgiosb

NB. ios=. getrfpl1uios0b ((m+1),n)
getrfpl1uios0b=: 3 : 0
  m=. <: {. 'm1 n'=. y
  nb=. TRFNB <. m <. n
  (i.0);(i.nb);(i.m);(i.m1);(m ht2ios nb);(n ht2ios nb);_1
)

NB. 'invpraw L1 U'=. ({: ; ((trl1 ; tru)@}:))@getrfpl1ub A
NB. (-: (clean@(((_1 0&ms@$) {. {:) swp2fw ((trl1 mp tru)@}:))@getrfpl1ub)) A

getrfpl1ub=: (getf2pl1ub`(0 {:: (getrfpl1ustepb ^: ((>.@(TRFNB %~ (_1 0&(ms $))))`(];(getrfpl1uios0b@$)))))@.(TRFNB<(_1 0&(ms $)))) @ (,&0)

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

NB. 'rawpi1 pfxi1 sfxLi1 sfxRi1'=. getf2pl1ustepc (rawpi ; pfxi ; sfxLi ; sfxRi)

getf2pl1ustepc=: 3 : 0
  'rawp pfx sfxL sfxR'=. y
  c=. {."1 sfxR
  jp=. iofmaxm c
  r=. jp { sfxR
  c=. c % ({. r)
  rawp=. rawp , jp
  pfx=. pfx , ((jp { sfxL) , r)
  sfxL=. }. jp ({.@])`[`]} sfxL ,. c
  sfxR=. 1 1 }. jp ({.@])`[`]} sfxR - (c */ r)
  rawp ; pfx ; sfxL ; sfxR
)

NB. 'p L1U'=. getf2pl1u A

getf2pl1uc=: ((0&{) , ((,&.>)/@(1 2&{))) @ (getf2pl1ustepc ^: ((<./@$)`(((i.0)"_);(0&{.);(_ 0&{.);])))

NB. 'pi1 pfxi1 sfxLi1 sfxRi1'=. getf2xxxxstep (pi ; pfxi ; sfxLi ; sfxRi)

getrfpl1ustepc=: 3 : 0
  'p pfx sfxL sfxR'=. y
  nj=. -~/ 'j n'=. $ pfx
  bs=. TRFNB <. (# sfxL) <. nj
  'pi L1U'=. getf2pl1uc bs {."1 sfxR
  sfxL=. pi laswp sfxL
  A1222=. pi laswp bs }."1 sfxR
  A11=. bs {. L1U
  A21=. bs }. L1U
  A12=. A11 trtrsl1x (bs {. A1222)
  p=. p , pi
  pfx=. pfx , ((bs {. sfxL) ,. A11 ,. A12)
  sfxL=. (bs }. sfxL) ,. A21
  sfxR=. (bs }. A1222) - A21 mp A12
  p ; pfx ; sfxL ; sfxR
)

NB. (-: (clean@((/: @ (0&({::))) { ((trl1 mp tru)@(1&({::))))@getrfpl1ua)) A
NB. 'ip L1U'=. getrfpl1ua A

getrfpl1uc=: getf2pl1uc`(((0&{),((,&.>)/@(1 2&{))) @ (getrfpl1ustepc ^: ((>.@(TRFNB %~ (<./@$)))`((i.@#);(0&{.);(_ 0 {. ]);]))))@.(TRFNB<(<./@$))

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

getrfpl1ud=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. m) ; y
  elseif. 1 = n do.
    dpi=. 0 ii2cp iofmaxm y
    pi=. dpi C. i. m
    y=. ((] 0:} %) (0: ({,) ])) dpi C. y                        NB. permute single column (row), scale by head, keep head unscaled
    pi ; y
  elseif. do.
    k=. m (<. >.@-:) n
    'pia Afa'=. getrfpl1ud k {."1 y                             NB. factorize 1st block recursively
    y=. pia C. k }."1 y                                         NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trtrsl1x & (k & {.)) y                          NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. getrfpl1ud y ((- (mp & Afba)) & (k & }.)) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. (i. k) , (k + pib)                                   NB. apply 2nd block's permutation to 1st block
    (dpib C. pia) ; ((dpib C. Afa) ,. (Afba , Afbb))            NB. assemble solution
  end.
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. getrfpl1u
NB.
NB. Description:
NB.   Triangular factorization of a general matrix:
NB.     P * L1 * U = A
NB.
NB. Syntax
NB.     'p L1U'=. getrfpl1u A
NB. where
NB.   A  - m*n matrix
NB.
NB. Notes:
NB. - this is a right-looking version of algorithm


NB. =========================================================
NB. Test suite

NB. #########################################################
NB. #########################################################
NB. #########################################################
NB. #########################################################
NB. #########################################################
NB. #########################################################

NB. ---------------------------------------------------------
NB. getrfp
NB. Template adverb to make verbs to triangular factorization
NB. with partial pivoting of a general matrix
NB.
NB. Syntax:
NB.   u0=. u0`u1`u2`u3`u4`u5`u6`u7`u8`u9`u10 getrfp
NB. where
NB.   u0  - factorization verb to do recursive call
NB.         warning: avoid to change u0's name class!
NB.   u1  - select number of rows (columns) from shape
NB.   u2  - find cycle permutation to swap pivot and head
NB.         (tail)
NB.   u3  - IO head (tail) in vector to pivot
NB.   u4  - permute columns (rows)
NB.   u5  - to form misc. rectangular IOS and deltas
NB.   u6  - solve unit triangular system for 2nd block's 1st
NB.         sub-block
NB.   u7  - multiply 1st block's 2nd sub-block by 2nd block's
NB.         1st sub-block
NB.   u8  - unite 2nd block's permutation parts
NB.   u9  - unite 1st and 2nd blocks
NB.   u10 - unite 2nd block's 1st and 2nd sub-blocks

getrfp=: 1 : 0
  '`u0 u1 u2 u3 u4 u5 u6 u7 u8 u9 u10'=. u
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. (m u1 n)) ; y
  elseif. 1 = (m u1~ n) do.
    dpi=. sh u2 y
    pi=. dpi C. i. (m u1 n)
    y=. ((] u3 } %) (u3 ({,) ])) dpi u4 y                NB. permute single column (row), scale by head, keep head unscaled
    pi ; y
  elseif. do.
    k=. m (u1 (<. >.@-:) u1~) n
    'i0 i1 i2 i3 d0 d1'=. u5 (0 _ 1 _1 * k) , sh
    'pia Afa'=. u0 i0 {. y                               NB. factorize 1st block recursively
    y=. pia u4 i1 }. y                                   NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (u6 & (i2 & {.)) y                        NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. u0 y ((- (u7 & Afba)) & (i3 & }.)) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. ((-/ d0) + i. k) u8 ((+/ d1) + pib)           NB. apply 2nd block's permutation to 1st block
    (dpib u4 pia) ; ((dpib u4 Afa) u9 (Afba u10 Afbb))   NB. assemble solution
  end.
)

NB. ---------------------------------------------------------
NB. getrf
NB. Template adverb to make verbs to triangular factorization
NB. without pivoting of a general matrix
NB.
NB. Syntax:
NB.   u0=. u0`u1`u2`u3`u4`u5`u6`u7`u8`u9`u10 getrf
NB. where
NB.   u0 - factorization verb to do recursive call
NB.        warning: avoid to change u0's name class!
NB.   u1 - select number of rows (columns) from shape
NB.   u2 - IO head (tail) in vector to pivot
NB.   u3 - permute columns (rows)
NB.   u4 - to form misc. rectangular IOS and deltas
NB.   u5 - solve unit triangular system for 2nd block's 1st
NB.        sub-block
NB.   u6 - multiply 1st block's 2nd sub-block by 2nd block's
NB.        1st sub-block
NB.   u7 - unite 1st and 2nd blocks
NB.   u8 - unite 2nd block's 1st and 2nd sub-blocks

getrf=: 1 : 0
  '`u0 u1 u2 u3 u4 u5 u6 u7 u8'=. u
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    y
  elseif. 1 = (m u1~ n) do.
    ((] u2 } %) (u2 ({,) ])) y                     NB. scale by head, keep head unscaled
  elseif. do.
    k=. m (u1 (<. >.@-:) u1~) n
    'i0 i1 i2 i3 d0 d1'=. u4 (0 _ 1 _1 * k) , sh
    Afa=. u0 i0 {. y                               NB. factorize 1st block recursively
    y=. pia u3 i1 }. y                             NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (u5 & (i2 & {.)) y                  NB. calculate 2nd block's 1st sub-block
    Afbb=. u0 y ((- (u6 & Afba)) & (i3 & }.)) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    Afa u7 (Afba u8 Afbb)                          NB. assemble solution
  end.
)

NB. ---------------------------------------------------------
NB. hetrfp
NB. Template adverb to make verbs to triangular factorization
NB. with full pivoting of a Hermitian (symmetric) matrix
NB.
NB. Syntax:
NB.   u0=. u1`u2`u3`u4`u5`u6`u7`u8`u9 hetrfp
NB. where
NB.   u0 - resulting factorization verb
NB.   u1 - to form iteration values
NB.   u2 - A's column IO to take from
NB.   u3 - select L (U) submatrix
NB.   u4 - IO pivoting element
NB.   u5 - IO element to pivot
NB.   u6 - make boxed cycle permutation to pivot
NB.   u7 - IO extreme element
NB.   u8 - IO opposite extreme element
NB.   u9 - append (prepend) to form next h
NB.
NB. Notes:
NB. - only A's lower (upper) triangle is referenced
NB. - T should be sparse

hetrfp=: 1 : 0
  '`u1 u2 u3 u4 u5 u6 u7 u8 u9'=. u
  n=. # y
  UL=. (u1 $ y) {. 1                   NB. 1st L's (last U's) column in Aasen method is e(1)-vector (e(n)-vector)
  T=. ((u1 @ $) {. ((u1 1 1) {. ])) y  NB. T[0,0]=A[0,0] (T[n-1,n-1]=A[n-1,n-1])
  pi=. i. n
  h=. i. 1                             NB. h[0:j-1]=H[0:j-1,j-1] (h[0:n-j-1]=H[j+1:n-1,j+1]), may be defined arbitrary before the 1st iteration only
  ios=. (u1 1) }. i. n                 NB. j:n-1 (0:j)
  for_j. }. i. u1 n do.                NB. 1:n-1 (n-2:0)
    a=. (< ios ; (u2 j)) { y           NB. A[j:n-1,j-1] (A[0:j,j+1])
    lum=. (n ((- u9 [)~ u3) j) {. UL   NB. L[j:n-1,0:j-1] (U[0:j,j+1:n-1])
    v=. a - lum mp h                   NB. v[0:n-j-1]=A[j:n-1,j-1]-UL[j:n-1,0:j-1]*H[0:j-1,j-1]=UL[j:n-1,j]*T[j,j-1] (v[0:j]=A[0:j,j+1]-U[0:j,j+1:n-1]*H[j+1:n-1,j+1]=U[0:j,j]*T[j,j+1]) non-pivoted yet
    q=. u4 v                           NB. IO pivot from head (tail)
    v=. ((u5 j) ii2cp q) C. v          NB. v[0]↔v[q] (v[_1]↔v[q])
    dpi=. q (u6 ii2cp ]) j             NB. any[j]↔any[j+q] (any[j]↔any[q])
    y=. dpi pt y                       NB. A[j,j:n-1]↔A[j+q,j:n-1], A[j:n-1,j]↔A[j:n-1,j+q] (A[j,0:j]↔A[q,0:j], A[0:j,j]↔A[0:j,q])
    pi=. dpi C. pi                     NB. pi[j]↔pi[j+q] (pi[j]↔pi[q])
    to=. ({~ u7) v                     NB. T[j,j-1] (T[j,j+1])
    lu=. q { lum                       NB. L[j,0:j-1] (U[j,j+1:n-1]) after pivoting
    luecto=. (+ ({~ u8) lu) * to       NB. conj(L[j,j-1])*T[j,j-1] (conj(U[j,j+1])*T[j,j+1])
    h=. to (((+ +)~ ({~ u8))`u8`]) } ((n (2 $ (u1~ u3)) j) {. T) mp + lu  NB. h[0:j-1]=H[0:j-1,j]=T[0:j-1,0:j-1]*conj(UL[j,0:j-1]) (h[0:n-j-1]=H[j+1:n-1,j]=T[j+1:n-1,j+1:n-1]*conj(U[j,j+1:n-1]))
    UL=. (1 u7 } v % to) (< ios ; j) } dpi C. UL                          NB. L[j,0:n-1]↔L[j+q,0:n-1], L[j:n-1,j]=v[0:n-j-1]/v[0], v[0] may be 0 (U[j,0:n-1]↔U[q,0:n-1], U[0:j,j]=v[0:j]/v[_1], v[_1] may be 0)
    td=. 9 o. ((< 2 $ j) { y) - ((lu mp h) + luecto)                      NB. T[j,j]=Re(A[j,j]-L[j,0:j-1]*H[0:j-1,j]-conj(L[j,j-1])*T[j,j-1]) (T[j,j]=Re(A[j,j]-U[j,j+1:n-1]*H[j+1:n-1,j]-conj(U[j,j+1])*T[j,j+1])
    T=. (td (,~ (, +)) to) (_2 <\ u1 0 _1 _1 0 0 0 + j) } T               NB. batch write diagonal and off-diagonal elements T[j,j-1] T[j-1,j] T[j,j] (T[j,j+1] T[j+1,j] T[j,j])
    h=. h u9 luecto + td               NB. h[0:j]=H[0:j,j]=T[0:J,0:J]*conj(L[j,0:j]) (h[0:j]=H[0:j,j]=T[0:J,0:J]*conj(U[j,0:j]))
    ios=. (u1 1) }. ios                NB. j+1:n-1 (0:j-1)
  end.
  UL ; T ; pi
)

NB. ---------------------------------------------------------
NB. hetrf
NB. Template adverb to make verbs to triangular factorization
NB. of a Hermitian (symmetric) matrix
NB.
NB. Syntax:
NB.   u0=. u1`u2`u3`u4`u5`u6`u7 hetrf
NB. where
NB.   u0 - resulting factorization verb
NB.   u1 - to form iteration values
NB.   u2 - A's column IO to take from
NB.   u3 - select L (U) submatrix
NB.   u4 - IO extreme element
NB.   u5 - IO opposite extreme element
NB.   u6 - append (prepend) to form next h
NB.
NB. Notes:
NB. - only A's lower (upper) triangle is referenced
NB. - T should be sparse

hetrf=: 1 : 0
  '`u1 u2 u3 u4 u5 u6 u7'=. u
  n=. # y
  UL=. (u1 $ y) {. 1                                  NB. 1st L's (last U's) column in Aasen method is e1-vector (e_(n-1)-vector)
  T=. ((u1 @ $) {. ((u1 1 1) {. ])) y                 NB. T[0,0]=A[0,0] (T[n-1,n-1]=A[n-1,n-1])
  h=. i. 1                                            NB. h[0:j-1]=H[0:j-1,j-1] (h[0:n-j-1]=H[j+1:n-1,j+1]), may be defined arbitrary before the 1st iteration only
  ios=. (u1 1) }. i. n                                NB. j:n-1 (0:j)
  for_j. }. i. u1 n do.                               NB. 1:n-1 (n-2:0)
    a=. (< ios ; (u2 j)) { y                          NB. A[j:n-1,j-1] (A[0:j,j+1])
    lum=. (n ((- u6 [)~ u3) j) {. UL                  NB. L[j:n-1,0:j-1] (U[0:j,j+1:n-1])
    v=. a - lum mp h                                  NB. v[0:n-j-1]=A[j:n-1,j-1]-UL[j:n-1,0:j-1]*H[0:j-1,j-1]=UL[j:n-1,j]*T[j,j-1] (v[0:j]=A[0:j,j+1]-U[0:j,j+1:n-1]*H[j+1:n-1,j+1]=U[0:j,j]*T[j,j+1])
    to=. ({~ u4) v                                    NB. T[j,j-1] (T[j,j+1])
    lu=. ({~ u4) lum                                  NB. L[j,0:j-1] (U[j,j+1:n-1])
    luecto=. (+ ({~ u5) lu) * to                      NB. conj(L[j,j-1])*T[j,j-1] (conj(U[j,j+1])*T[j,j+1])
    h=. to (((+ +)~ ({~ u5))`u5`]) } ((n (2 $ (u1~ u3)) j) {. T) mp + lu  NB. h[0:j-1]=H[0:j-1,j]=T[0:j-1,0:j-1]*conj(UL[j,0:j-1]) (h[0:n-j-1]=H[j+1:n-1,j]=T[j+1:n-1,j+1:n-1]*conj(U[j,j+1:n-1]))
    UL=. (1 u4 } v % to) (< ios ; j) } UL                                 NB. L[j:n-1,j]=v[0:n-j-1]/v[0], v[0] may be 0 (U[0:j,j]=v[0:j]/v[_1], v[_1] may be 0)
    td=. 9 o. ((< 2 $ j) { y) - ((lu mp h) + luecto)                      NB. T[j,j]=Re(A[j,j]-L[j,0:j-1]*H[0:j-1,j]-conj(L[j,j-1])*T[j,j-1]) (T[j,j]=Re(A[j,j]-U[j,j+1:n-1]*H[j+1:n-1,j]-conj(U[j,j+1])*T[j,j+1])
    T=. (td (,~ (, +)) to) ((_2 <\ u1 0 _1 _1 0 0 0 + j) j) } T           NB. batch write diagonal and off-diagonal elements T[j,j-1] T[j-1,j] T[j,j] (T[j,j+1] T[j+1,j] T[j,j])
    h=. h u6 luecto + td                              NB. h[0:j]=H[0:j,j]=T[0:J,0:J]*conj(L[j,0:j]) (h[0:j]=H[0:j,j]=T[0:J,0:J]*conj(U[j,0:j]))
    ios=. (u1 1) }. ios                               NB. j+1:n-1 (0:j-1)
  end.
  UL ; T
)

NB. ---------------------------------------------------------
NB. potrf
NB. Template adverb to make Cholesky factorization verbs of
NB. a Hermitian (symmetric) positive definite matrix
NB.
NB. Syntax:
NB.   v0=. v0`u1`u2`u3`u4`u5 potrf
NB. where
NB.   v0 - factorization verb for recursive call
NB.   u1 - either }. for upper triangular factor
NB.        or {. for lower triangular factor
NB.   u2 - either trtrsux for upper triangular factor, or
NB.        trtrslx for lower triangular factor
NB.   u3 - either {. for upper triangular factor
NB.        or }. for lower triangular factor
NB.   u4 - either ,. for upper triangular factor
NB.        or ,.~ for [unit] lower triangular factor
NB.   u5 - either (_1 append) for upper triangular factor
NB.        or (0 append~) for lower triangular factor

potrf=: 1 : 0
  '`u0 u1 u2 u3 u4 u5'=. u
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    Ta=. u0 (2 $ p) u1 y     NB. recursively factorize square matrix from top left or bottom right corner
    Ac=. (p , (p - n)) u1 y  NB. off-diagonal part of input matrix
    Tb=. ct Tbh=. Ta u2 Ac   NB. off-diagonal part of output matrix
    Aa=. (2 $ p) u3 y        NB. square matrix from opposite corner on diagonal of input matrix
    Tc=. u0 Aa - Tb mp Tbh   NB. recursively factorize square matrix from opposite corner on diagonal
    (Tc u4 Tb) u5 Ta         NB. assemble output as triangular matrix
  else.
    %: y
  end.
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:          Solves:                  Syntax:
NB. getrfpl1u      inv(P) * L1 * U = A      'pi Af'=. getrfpl1u A
NB. getrflu1p      L * U1 * inv(P) = A      'pi Af'=. getrflu1p A
NB. getrfpu1l      inv(P) * U1 * L = A      'pi Af'=. getrfpu1l A
NB. getrful1p      U * L1 * inv(P) = A      'pi Af'=. getrful1p A
NB.
NB. Description:
NB.   triangular factorization with partial pivoting of a
NB.   general matrix
NB. where:
NB.   A  - m×n-matrix, containing either U, U1, L or L1
NB.   Af - m×n-matrix, contains either U and L1, or U1 and L
NB.   pi - m-vector or n-vector, inversed rows (columns)
NB.        permutation of A
NB.   U  - ?×?-matrix, upper triangular factor
NB.   U1 - ?×?-matrix, unit upper triangular matrix (diagonal is not saved)
NB.   L  - ?×?-matrix, lower triangular matrix
NB.   L1 - ?×?-matrix, unit lower triangular matrix (diagonal is not saved)
NB.   m  ≥ 0
NB.   n  ≥ 0
NB.
NB. If:
NB.   'pi Af'=. getrfpl1u A
NB.   p=. /: pi
NB.   Pi=. ((i.#pi)=/pi)
NB.   P=. ((i.#p)=/p)
NB.   L1=. trl1 Af
NB.   U=. tru Af
NB. then (with appropriate comparison tolerance)
NB.   Pi -: %. P
NB.   Pi -: |: P
NB.   A -: p { L1 mp U
NB.   A -: p C. L1 mp U
NB.   A -: P mp L1 mp U
NB.
NB. Notes:
NB.   - FLOPs: 
NB.   - recursive calls: 3*k*(2^k)-2^(k+1)+3,
NB.     where k = ⌈log_2(n)⌉

getrfpl1u=: getrfpl1u`[`((ii2cp~ 0:   )~ iofmaxm)` 0:` C.   `((4 2,0 2,2 1,2 0,0 0,:0 2)&{)`trtrsl1x` mp  ` ,  ` ,.  ` ,    getrfp
getrflu1p=: getrflu1p`]`((ii2cp~ 0:   )~ iofmaxm)` 0:`(C."1)`((2 5,2 0,1 2,0 2,0 0,:0 2)&{)`trtrsxu1`(mp~)` ,  ` ,   ` ,.   getrfp
getrfpu1l=: getrfpu1l`[`((ii2cp~ <:@[/)~ iolmaxm)`_1:` C.   `((4 3,0 3,3 1,3 0,4 2,:0 0)&{)`trtrsu1x` mp  `(,~)`(,.~)`(,~ ) getrfp
getrful1p=: getrful1p`]`((ii2cp~ <:@]/)~ iolmaxm)`_1:`(C."1)`((3 5,3 0,1 3,0 3,5 2,:0 0)&{)`trtrsxl1`(mp~)`(,~)`(,~ )`(,.~) getrfp

NB. ---------------------------------------------------------
NB. Verb:          Solves:                  Syntax:
NB. getrfl1u       L1 * U = A               Af=. getrfl1u A
NB. getrflu1       L * U1 = A               Af=. getrflu1 A
NB. getrfu1l       U1 * L = A               Af=. getrfu1l A
NB. getrful1       U * L1 = A               Af=. getrful1 A
NB.
NB. Description:
NB.   triangular factorization without pivoting of a
NB.   general matrix
NB. where
NB.   A  - m×n-matrix, containing either U, U1, L or L1
NB.   Af - m×n-matrix, contains either U and L1, or U1 and L
NB.   U  - ?×?-matrix, upper triangular factor
NB.   U1 - ?×?-matrix, unit upper triangular matrix (diagonal is not saved)
NB.   L  - ?×?-matrix, lower triangular matrix
NB.   L1 - ?×?-matrix, unit lower triangular matrix (diagonal is not saved)
NB.   m  ≥ 0
NB.   n  ≥ 0
NB.
NB. If: <<<<<<<<<<<<<<<<<<
NB.   'pi Af'=. getrfpl1u A
NB.   p=. /: pi
NB.   Pi=. ((i.#pi)=/pi)
NB.   P=. ((i.#p)=/p)
NB.   L1=. trl1 Af
NB.   U=. tru Af
NB. then (with appropriate comparison tolerance)
NB.   Pi -: %. P
NB.   Pi -: |: P
NB.   A -: p { L1 mp U
NB.   A -: p C. L1 mp U
NB.   A -: P mp L1 mp U
NB.
NB. Notes:
NB.   - FLOPs: 
NB.   - recursive calls: 3*k*(2^k)-2^(k+1)+3,
NB.     where k = ⌈log_2(n)⌉

getrfl1u=: getrfpl1u`[` 0:` C.   `((4 2,0 2,2 1,2 0,0 0,:0 2)&{)`trtrsl1x` mp  ` ,.  ` ,    getrf
getrflu1=: getrflu1p`]` 0:`(C."1)`((2 5,2 0,1 2,0 2,0 0,:0 2)&{)`trtrsxu1`(mp~)` ,   ` ,.   getrf
getrfu1l=: getrfpu1l`[`_1:` C.   `((4 3,0 3,3 1,3 0,4 2,:0 0)&{)`trtrsu1x` mp  `(,.~)`(,~ ) getrf
getrful1=: getrful1p`]`_1:`(C."1)`((3 5,3 0,1 3,0 3,5 2,:0 0)&{)`trtrsxl1`(mp~)`(,~ )`(,.~) getrf

NB. ---------------------------------------------------------
NB. Verb:      Solves:                          Syntax:
NB. hetrfpl    inv(P) * L1 * T * L1' * P = A    'L T pi'=. hetrfpl A
NB. hetrfpu    inv(P) * U1 * T * U1' * P = A    'U T pi'=. hetrfpu A
NB.
NB. Description:
NB.   Factorization of a Hermitian (symmetric) matrix with
NB.   full pivoting by Aasen's method
NB. where
NB.   A  - n×n-matrix, Hermitian (symmetric)
NB.   pi - n-vector, inversed rows and columns permutation
NB.        of A
NB.   U1 - n×n-matrix, unit upper triangular
NB.   L1 - n×n-matrix, unit lower triangular
NB.   T  - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   n  ≥ 0
NB.
NB. If:
NB.   'U T pi'=. hetrfpu A
NB.   p=. /: pi
NB.   Pi=. ((i.#pi)=/pi)
NB.   P=. ((i.#p)=/p)
NB.   U1=. tru1 U
NB. then (with appropriate comparison tolerance)
NB.   Pi -: %. P
NB.   Pi -: |: P
NB.   U1 -: U
NB.   A -: Pi mp U mp T mp ct U mp P
NB.   A -: (/: pi) pt U mp T mp ct U
NB.
NB. Notes:
NB.   - FLOPs:

hetrfpl=: [`<:`] `iofmaxm`0:`+`0: `_1:`,    hetrfp
hetrfpu=: -`>:`>:`iolmaxm`] `[`_1:`0: `(,~) hetrfp

NB. ---------------------------------------------------------
NB. Verb:      Solves:                          Syntax:
NB. hetrfl     L1 * T * L1' = A                 'L T'=. hetrfl A
NB. hetrfu     U1 * T * U1' = A                 'U T'=. hetrfu A
NB.
NB. Description:
NB.   Factorization of a Hermitian (symmetric) matrix by
NB.   Aasen's method
NB. where
NB.   A  - n×n-matrix, Hermitian (symmetric)
NB.   U1 - n×n-matrix, unit upper triangular
NB.   L1 - n×n-matrix, unit lower triangular
NB.   T  - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   n  ≥ 0
NB.
NB. If:
NB.   'U T'=. hetrfu A
NB.   U1=. tru1 U
NB. then (with appropriate comparison tolerance)
NB.   U1 -: U
NB.   A -: U mp T mp ct U
NB.
NB. Notes:
NB.   - FLOPs:

hetrfl=: [`<:`] `0: `_1:`,    hetrf
hetrfu=: -`>:`>:`_1:`0: `(,~) hetrf

NB. ---------------------------------------------------------
NB. potrfu
NB. potrfl
NB. Hermitian (symmetric) positive definite matrix Cholesky
NB. recursive factorization by a triangular factor
NB.   U * U' = A
NB.   L * L' = A
NB.
NB. Syntax:
NB.   U=. potrfu A
NB.   L=. potrfl A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric), positive
NB.       definite
NB.   U  - n×n-matrix, upper triangular Cholesky factor
NB.   L  - n×n-matrix, lower triangular Cholesky factor
NB.   n      ≥ 0
NB.
NB. If:
NB.   U=. potrfl A
NB.   L=. potrfl A
NB. then (with appropriate comparison tolerance)
NB.   A -: (mp ct) U
NB.   A -: (mp ct) L
NB. where
NB.   A - positive definite matrix
NB.   U - upper triangular matrix
NB.   L - lower triangular matrix
NB.
NB. Notes:
NB.   - FLOPs: 
NB.   - recursive calls: 3*k*(2^k)-2^(k+1)+3,
NB.     where k = ⌈log_2(n)⌉

potrfu=: potrfu`}.`trtrsux`{.` ,.  `(_1 append ) potrf
potrfl=: potrfl`{.`trtrslx`}.`(,.~)`( 0 append~) potrf

NB. =========================================================
NB. Test suite

NB. name vextract ttrf A
NB. TODO: forward error (getrf has two: for L and U?)

ttrf=: 1 : 0
:
  't s'=. timespacex 'out=. ' , x , ' y'
  be=. (((norm1 (y - u out)) % ({: $ y)) % (norm1 y)) % FP_EPS  NB. backward error
  prn x ; ((_."_)`(norm1 con getri) @. (=/@$) y) ; be ; (i.0) ; t ; s
)

NB. tgetrf A

tgetrf=: 3 : 0
  y=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y

  'getrfpl1u'       (((C.  ~ /:)~ (trl1 mp tru )) & > /)               ttrf y
  'getrflu1p'       (((C."1~ /:)~ (trl  mp tru1)) & > /)               ttrf y
  'getrfpu1l'       (((C.  ~ /:)~ (tru1 mp trl )) & > /)               ttrf y
  'getrful1p'       (((C."1~ /:)~ (tru  mp trl1)) & > /)               ttrf y
  'getrf_jlapack_' (((mp & >)/ @ (2 & {.)) invperm_jlapack_ (2 & {::)) ttrf y
  EMPTY
)

NB. thetrf A===============

thetrfp=: 3 : 0
  if. (i. 0) -: $ y do.
    y=. (_1 1 0 16 _6 4 & (gemat j. gemat)) hemat y
  end.
  'L T pi'=. hetrfpl y
  Asol=: (/: pi) pt (L mp T mp ct L)
  smoutput 'A' ; ($ y) ; y ; 'L' ; ($ L) ; L ; 'T' ; ($ T) ; T ; 'Asol' ; ($ Asol) ; Asol
  smoutput 'A -: Asol' ; (y -: Asol) ; '||A - Asol||/||A|| in 1-norm' ; ((y - Asol) ((% (FP_SFMIN & >.)) & norm1) y) ; '||A - Asol||/||A|| in i-norm' ; ((y - Asol) ((% (FP_SFMIN & >.)) & normi) y)
  'U T pi'=. hetrfpu y
  Asol=: (/: pi) pt (U mp T mp ct U)
  smoutput 'A' ; ($ y) ; y ; 'U' ; ($ U) ; U ; 'T' ; ($ T) ; T ; 'Asol' ; ($ Asol) ; Asol
  smoutput 'A -: Asol' ; (y -: Asol) ; '||A - Asol||/||A|| in 1-norm' ; ((y - Asol) ((% (FP_SFMIN & >.)) & norm1) y) ; '||A - Asol||/||A|| in i-norm' ; ((y - Asol) ((% (FP_SFMIN & >.)) & normi) y)
)



NB. tpotrf A

tpotrf=: 3 : 0
  y=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y

  'potrfu'         ((mp ct) @ tru) ttrf y
  'potrfl'         ((mp ct) @ trl) ttrf y
  'potrf_jlapack_' ((mp ct) @ trl) ttrf y
  EMPTY
)

NB. ---------------------------------------------------------
NB. testtrf
NB. Adverb to test triangular factorization algorithms
NB.
NB. Syntax:
NB.   mkge testtrf m,n
NB. where
NB.   m,n  - 2-vector of integers, shape of random matrices
NB.          to test algorithms
NB.   mkge - monadic verb to generate random non-singular
NB.          general y-matrix (shape is taken from y)
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   cocurrent 'mt'
NB.   (_1 1 0 16 _6 4 & gemat) testtrf 500 500
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testtrf 500 500
NB.
NB. Notes:
NB. - if m≠n then both thetrf and tpotrf are skipped

testtrf=: 1 : 0
  (tgetrf @ u)                      y
  (thetrf @ (u hemat) @ {. ^: (=/)) y
  (tpotrf @ (u pomat) @ {. ^: (=/)) y
  EMPTY
)
