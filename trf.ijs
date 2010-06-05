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
NB. hetrfu     UDU' factorization of a Hermitian (symmetric)
NB.            matrix, where U is unit upper triangular and
NB.            D is diagonal
NB. hetrfl     LDL' factorization of a Hermitian (symmetric)
NB.            matrix, where L is unit lower triangular and
NB.            D is diagonal
NB.
NB. hetrfu     UDU' factorization of a Hermitian (symmetric)
NB.            matrix, where U is unit upper triangular and
NB.            D is diagonal
NB. hetrfl     LDL' factorization of a Hermitian (symmetric)
NB.            matrix, where L is unit lower triangular and
NB.            D is diagonal
NB.
NB. potrfl     LL' (Cholesky) factorization of a Hermitian
NB.            (symmetric) positive definite matrix for lower
NB.            triangular factor
NB. potrfu     UU' (Cholesky) factorization of a Hermitian
NB.            (symmetric) positive definite matrix for upper
NB.            triangular factor
NB.
NB. Copyright (C) 2009  Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

iofmaxm=: (i.>./) @: |  NB. IO 1st element with max magnitude from list y
iolmaxm=: (i:>./) @: |  NB. IO last element with max magnitude from list y

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
NB.   u0=. u1`u2`u3`u4`u5`u6`u7`u8`u9`u10`u11`u12`u13`u14`u15`u16 hetrfp
NB. where
NB.   u0  - factorization verb to find inv(P)*L*D*L'*P=A
NB.         (inv(P)*U*D*U'*P=A)
NB.   u1  - allocate L (U) and initialize 1st (last) column
NB.   u2  - allocate T and initialize T[0,0] (T[n-1,n-1])
NB.   u3  - cut off 1st (last) element
NB.   u4  - to form iteration values
NB.   u5  - A's column IO to take from
NB.   u6  - select L (U) submatrix
NB.   u7  - IO pivoting element
NB.   u8  - IO element to pivot
NB.   u9  - make boxed cycle permutation to pivot
NB.   u10 - element off the T's diagonal
NB.   u11 - extreme last (1st) element
NB.   u12 - IO extreme last (1st) element
NB.   u13 - select T submatrix
NB.   u14 - IO extreme element to replace by 1
NB.   u15 - IOS to write into T
NB.   u16 - append (prepend) to form next h
NB.
NB. Notes:
NB. - only A's lower (upper) triangle is referenced
NB. - T should be sparse

hetrfp=: 1 : 0
  '`u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16'=. u
  n=. # y
  UL=. (u1 y) {. 1                                          NB. 1st (last) column in Aasen method is e1-vector (e_(n-1)-vector)
  T=. (u1 {. u2) y                                          NB. T[0,0]=A[0,0] (T[n-1,n-1]=A[n-1,n-1])
  pi=. i. n
  h=. i. 1                                                  NB. h[0:j-1]=H[0:j-1,j-1] (h[0:n-j-1]=H[j+1:n-1,j+1]), may be defined arbitrary before the 1st iteration only
  ios=. u3 i. n                                             NB. j..(n-1) (0..j)
  for_j. }. i. u4 n do.                                     NB. 1..(n-1) ((n-2)..0)
    a=. (< ios ; (u5 j)) { y                                NB. A[j:n-1,j-1] (A[0:j,j+1])
    lum=. (j u6 n) {. UL                                    NB. L[j:n-1,0:j-1] (U[0:j,j+1:n-1])
    v=. a - lum mp h                                        NB. v[0:n-j-1]=A[j:n-1,j-1]-UL[j:n-1,0:j-1]*H[0:j-1,j-1]=UL[j:n-1,j]*T[j,j-1] (v[0:j]=A[0:j,j+1]-U[0:j,j+1:n-1]*H[j+1:n-1,j+1]=U[0:j,j]*T[j,j+1]) non-pivoted yet
    q=. u7 v                                                NB. IO pivot from head (tail)
    v=. ((u8 j) ii2cp q) C. v                               NB. v[0]↔v[q] (v[_1]↔v[q])
    dpi=. q u9 j                                            NB. any[j]↔any[j+q] (any[j]↔any[q])
    y=. dpi pt y                                            NB. A[j,j:n-1]↔A[j+q,j:n-1], A[j:n-1,j]↔A[j:n-1,j+q] (A[j,0:j]↔A[q,0:j], A[0:j,j]↔A[0:j,q])
    pi=. dpi C. pi                                          NB. pi[j]↔pi[j+q] (pi[j]↔pi[q])
    to=. u10 v                                              NB. T[j,j-1] (T[j,j+1])
    lu=. q { lum                                            NB. L[j,0:j-1] (U[j,j+1:n-1]) after pivoting
    luceto=. (+ u11 lu) * to                                NB. conj(L[j,j-1])*T[j,j-1] (conj(U[j,j+1])*T[j,j+1])
    h=. to (((+ +)~ u11)`u12`]) } ((j u13 n) {. T) mp + lu  NB. h[0:j-1]=H[0:j-1,j]=T[0:j-1,0:j-1]*conj(UL[j,0:j-1]) (h[0:n-j-1]=H[j+1:n-1,j]=T[j+1:n-1,j+1:n-1]*conj(U[j,j+1:n-1]))
    UL=. (1 u14 } v % to) (< ios ; j) } dpi C. UL           NB. L[j,0:n-1]↔L[j+q,0:n-1], L[j:n-1,j]=v[0:n-j-1]/v[0], v[0] may be 0 (U[j,0:n-1]↔U[q,0:n-1], U[0:j,j]=v[0:j]/v[_1], v[_1] may be 0)
    td=. 9 o. ((< 2 $ j) { y) - ((lu mp h) + luceto)        NB. T[j,j]=Re(A[j,j]-L[j,0:j-1]*H[0:j-1,j]-conj(L[j,j-1])*T[j,j-1]) (T[j,j]=Re(A[j,j]-U[j,j+1:n-1]*H[j+1:n-1,j]-conj(U[j,j+1])*T[j,j+1])
    T=. (td (,~ (, +)) to) (u15 j) } T                      NB. batch write T[j,j-1] T[j-1,j] T[j,j] (T[j,j+1] T[j+1,j] T[j,j])
    h=. h u16 luceto + td                                   NB. h[0:j]=H[0:j,j]=T[0:J,0:J]*conj(L[j,0:j]) (h[0:j]=H[0:j,j]=T[0:J,0:J]*conj(U[j,0:j]))
    ios=. u3 ios                                            NB. (j+1)..(n-1) (0..(j-1))
  end.
  UL ; T ; pi
)

NB. ---------------------------------------------------------
NB. potrf
NB. Template adverb to make Cholesky factorization verbs of
NB. a Hermitian (symmetric) positive definite matrix
NB.
NB. Syntax:
NB.   vpotrf=. vpotrf`u1`u2`u3`u4`u5 potrf
NB. where
NB.   vpotrf - factorization verb for recursive call
NB.   u1     - either }. for upper triangular factor
NB.            or {. for lower triangular factor
NB.   u2     - either trtrsux for upper triangular factor, or
NB.            trtrslx for lower triangular factor
NB.   u3     - either {. for upper triangular factor
NB.            or }. for lower triangular factor
NB.   u4     - either ,. for upper triangular factor
NB.            or ,.~ for [unit] lower triangular factor
NB.   u5     - either (_1 append) for upper triangular factor
NB.            or (0 append~) for lower triangular factor

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
NB. then
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
NB. then
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
NB.   Factorization of a  Hermitian (symmetric) matrix full
NB.   pivoting by Aasen's method
NB. where
NB.   A  - n×n-matrix, Hermitian (symmetric)
NB.   pi - n-vector, inversed rows and columns permutation
NB.        of A
NB.   U1 - n×n-matrix, unit upper triangular matrix
NB.   L1 - n×n-matrix, unit lower triangular matrix
NB.   D  - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   n  ≥ 0

hetrfpl=: $    `(1 1&{.)  `}.`]`<:`(-,[)        `iofmaxm`0:`(] ii2cp +)`{.`{:`_1:`(2&$@[)     `0: `(_2 <\ 0 _1 _1 0 0 0&+)`,    hetrfp
hetrfpu=: (-@$)`(_1 _1&{.)`}:`-`>:`(((],-~)>:)~)`iolmaxm`] `ii2cp      `{:`{.`0: `(2&$@(- <:))`_1:`(_2 <\ 0 1 1 0 0 0&+)  `(,~) hetrfp

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
NB. then
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

  'getrfpl1u'       (((C.~ /:)~ (trl1 mp tru)) & > /)                   ttrf y
  'getrf_jlapack_' (((mp & >)/ @ (2 & {.)) invperm_jlapack_ (2 & {::)) ttrf y
  EMPTY
)

NB. thetrf A===============

thetrfp=: 3 : 0
  if. (i. 0) -: $ y do.
    y=. (_1 1 0 16 _6 4 & (gemat j. gemat)) hemat y
  end.
  'L T pi'=: hetrfpl y
  Asol=: (/: pi) pt (L mp T mp ct L)
  smoutput 'A' ; ($ y) ; y ; 'L' ; ($ L) ; L ; 'T' ; ($ T) ; T ; 'Asol' ; ($ Asol) ; Asol
  smoutput 'A -: Asol' ; (y -: Asol) ; '||A - Asol||/||A|| in 1-norm' ; ((y - Asol) ((% (FP_SFMIN & >.)) & norm1) y) ; '||A - Asol||/||A|| in i-norm' ; ((y - Asol) ((% (FP_SFMIN & >.)) & normi) y)
  'U T pi'=: hetrfpu y
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
