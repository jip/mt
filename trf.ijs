NB. trf.ijs
NB. Triangular factorization of matrix
NB.
NB. getrflu1p  Triangular factorization with partial pivoting
NB.            of a general matrix: L * U1 * P = A, where P
NB.            is permutation matrix, L is lower triangular
NB.            and U1 is unit upper triangular
NB. getrfpl1u  Triangular factorization with partial pivoting
NB.            of a general matrix: P * L1 * U = A, where P
NB.            is permutation matrix, L1 is unit lower
NB.            triangular and U is upper triangular
NB. getrfpu1l  Triangular factorization with partial pivoting
NB.            of a general matrix: P * U1 * L = A, where P
NB.            is permutation matrix, U1 is unit upper
NB.            triangular and L is lower triangular
NB. getrful1p  Triangular factorization with partial pivoting
NB.            of a general matrix: U * L1 * P = A, where P
NB.            is is permutation matrix, U is upper
NB.            triangular and L1 is unit lower triangular
NB.
NB. hetrfpl    Triangular factorization with full pivoting of
NB.            a Hermitian (symmetric) matrix:
NB.            P * L1 * T * L1' * P' = A, where P is
NB.            permutation matrix, L1 is unit lower
NB.            triangular and T is tridiagonal
NB. hetrfpu    Triangular factorization with full pivoting of
NB.            a Hermitian (symmetric) matrix:
NB.            P * U1 * T * U1' * P' = A, where P is
NB.            permutation matrix, U1 is unit upper
NB.            triangular and T is tridiagonal
NB.
NB. potrfl     Cholesky factorization of a Hermitian
NB.            (symmetric) positive definite matrix:
NB.            L * L' = A, where L is lower triangular
NB. potrfu     Cholesky factorization of a Hermitian
NB.            (symmetric) positive definite matrix:
NB.            U * U' = A, where U is upper triangular
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-10

coclass 'mt'

NB. =========================================================
NB. Local definitions

iofmaxm=: (i.>./) @: |  NB. IO 1st element with max magnitude from list y
iolmaxm=: (i:>./) @: |  NB. IO last element with max magnitude from list y

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
NB. Syntax
NB.     'ip LU1'=. getrflu1p A
NB. where
NB.   A   - m×n-matrix to factorize
NB.   ip  - n-vector, inversed columns permutation of A
NB.   LU1 - m×n matrix, lower triangle contains L, and strict
NB.         upper triangle contains U1, diagonal isn't stored
NB.   L   - m×k-matrix, lower triangular
NB.   U1  - k×n-matrix, unit upper triangular
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'ip LU1'=. getrflu1p A
NB.   p=. /: pi
NB.   iP=. p2P ip
NB.   P=. p2P p
NB. then
NB.   P -: %. iP
NB.   P -: |: iP
NB.   ((-:!.(2^_34)) ((/: @ (0&({::))) {"1 ((trl mp tru1)@(1&({::))))@getrflu1p) A
NB.   ((-:!.(2^_34)) ((/: @ (0&({::))) C."1 ((trl mp tru1)@(1&({::))))@getrflu1p) A
NB.   ((-:!.(2^_34)) (((trl mp tru1)@(1&({::))) mp (ip2P @ (0&({::))))@getrflu1p) A

getrflu1p=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. n) ; y
  elseif. 1=m do.
    dpi=. 0 ios2cp iofmaxm y
    pi=. dpi C. i. n
    y=. ((] 0:} %) (0&({,))) dpi (C."1) y                          NB. permute single column (row), scale by head, keep head unscaled
    pi ; y
  elseif. do.
    k=. n (<. >.@-:) m
    'pia Afa'=. getrflu1p k {. y                                   NB. factorize 1st block recursively
    y=. pia (C."1) k }. y                                          NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trtrsxu1 & (k & ({."1))) y                         NB. calculate 2nd block's 1st sub-block
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
NB. Syntax
NB.     'ip L1U'=. getrfpl1u A
NB. where
NB.   A   - m×n-matrix to factorize
NB.   ip  - m-vector, inversed rows permutation of A
NB.   L1U - m×n matrix, upper triangle contains U, and strict
NB.         lower triangle contains L1, diagonal isn't stored
NB.   L1  - m×k-matrix, unit lower triangular
NB.   U   - k×n-matrix, upper triangular
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'ip L1U'=. getrfpl1u A
NB.   p=. /: pi
NB.   iP=. p2P ip
NB.   P=. p2P p
NB. then
NB.   P -: %. iP
NB.   P -: |: iP
NB.   ((-:!.(2^_34)) ((/: @ (0&({::))) { ((trl1 mp tru)@(1&({::))))@getrfpl1u) A
NB.   ((-:!.(2^_34)) ((/: @ (0&({::))) C. ((trl1 mp tru)@(1&({::))))@getrfpl1u) A
NB.   ((-:!.(2^_34)) ((ip2P @ (0&({::))) mp ((trl1 mp tru)@(1&({::))))@getrfpl1u) A

getrfpl1u=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. m) ; y
  elseif. 1 = n do.
    dpi=. 0 ios2cp iofmaxm y
    pi=. dpi C. i. m
    y=. ((] 0:} %) (0&({,))) dpi C. y                          NB. permute single column (row), scale by head, keep head unscaled
    pi ; y
  elseif. do.
    k=. m (<. >.@-:) n
    'pia Afa'=. getrfpl1u k {."1 y                             NB. factorize 1st block recursively
    y=. pia C. k }."1 y                                        NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trtrsl1x & (k & {.)) y                         NB. calculate 2nd block's 1st sub-block
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
NB. Syntax
NB.     'ip U1L'=. getrfpu1l A
NB. where
NB.   ip  - m-vector, inversed rows permutation of A
NB.   A   - m×n-matrix to factorize
NB.   U1L - m×n matrix, lower triangle contains L, and strict
NB.         upper triangle contains U1, diagonal isn't stored
NB.   L   - m×k-matrix, lower triangular
NB.   U1  - k×n-matrix, unit upper triangular
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'ip U1L'=. getrfpu1l A
NB.   p=. /: pi
NB.   iP=. p2P ip
NB.   P=. p2P p
NB. then
NB.   P -: %. iP
NB.   P -: |: iP
NB.   ((-:!.(2^_34)) ((/: @ (0&({::))) { (((tru1~ (-~/@$)) mp (trl~ (-~/@$)))@(1&({::))))@getrfpu1l) A
NB.   ((-:!.(2^_34)) ((/: @ (0&({::))) C. (((tru1~ (-~/@$)) mp (trl~ (-~/@$)))@(1&({::))))@getrfpu1l) A
NB.   ((-:!.(2^_34)) ((ip2P @ (0&({::))) mp (((tru1~ (-~/@$)) mp (trl~ (-~/@$)))@(1&({::))))@getrfpu1l) A

getrfpu1l=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. m) ; y
  elseif. 1 = n do.
    dpi=. (<:m) ios2cp iolmaxm y
    pi=. dpi C. i. m
    y=. ((] _1:} %) (_1&({,))) dpi C. y                           NB. permute single column (row), scale by head, keep head unscaled
    pi ; y
  elseif. do.
    k=. m (<. >.@-:) n
    'pia Afa'=. getrfpu1l (-k) {."1 y                             NB. factorize 1st block recursively
    y=. pia C. (-k) }."1 y                                        NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trtrsu1x & ((-k) & {.)) y                         NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. getrfpu1l y ((- (mp & Afba)) & ((-k) & }.)) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. pib , ((m-k) + i. k)                                   NB. FIXME! remove 2nd; apply 2nd block's permutation to 1st block
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
NB. Syntax
NB.     'ip UL1'=. getrful1p A
NB. where
NB.   A   - m×n-matrix to factorize
NB.   ip  - n-vector, inversed columns permutation of A
NB.   UL1 - m×n matrix, upper triangle contains U, and strict
NB.         lower triangle contains L1, diagonal isn't stored
NB.   L   - m×k-matrix, unit lower triangular
NB.   U1  - k×n-matrix, upper triangular
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'p UL1'=. getrful1p A
NB.   P=. p2P p
NB.   ip=. /: p
NB.   iP=. p2P ip
NB. then
NB.   P -: %. iP
NB.   P -: |: iP
NB.   ((-:!.(2^_34)) ((/: @ (0&({::))) {"1 (((tru~ (-~/@$)) mp (trl1~ (-~/@$)))@(1&({::))))@getrful1p) A
NB.   ((-:!.(2^_34)) ((/: @ (0&({::))) C."1 (((tru~ (-~/@$)) mp (trl1~ (-~/@$)))@(1&({::))))@getrful1p) A
NB.   ((-:!.(2^_34)) ((((tru~ (-~/@$)) mp (trl1~ (-~/@$)))@(1&({::))) mp (p2P @ (0&({::))))@getrful1p) A

getrful1p=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. n) ; y
  elseif. 1 = m do.
    dpi=. (<:n) ios2cp iolmaxm y
    pi=. dpi C. i. n
    y=. ((] _1:} %) (_1&({,))) dpi (C."1) y                           NB. permute single column (row), scale by head, keep head unscaled
    pi ; y
  elseif. do.
    k=. n (<. >.@-:) m
    'pia Afa'=. getrful1p (-k) {. y                                 NB. factorize 1st block recursively
    y=. pia (C."1) (-k) }. y                                        NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trtrsxl1 & ((-k) & ({."1))) y                         NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. getrful1p y ((- (Afba & mp)) & ((-k) & (}."1))) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. pib , ((n-k) + i. k)                                       NB. FIXME! remove 2nd; apply 2nd block's permutation to 1st block
    (dpib (C."1) pia) ; ((Afbb ,. Afba) , (dpib (C."1) Afa))          NB. assemble solution
  end.
)

NB. ---------------------------------------------------------
NB. hetrfpl
NB.
NB. Description:
NB.   Triangular factorization with full pivoting of a
NB.   Hermitian (symmetric) matrix:
NB.     P * L1 * T * L1' * P' = A
NB.
NB. Syntax:
NB.     'ip L1 T'=. hetrfpl A
NB. where
NB.   A   - n×n-matrix to factorize, Hermitian (symmetric)
NB.   ip  - n-vector, inversed full permutation of A
NB.   L1  - n×n-matrix, unit lower triangular
NB.   T   - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'ip L1 T'=. hetrfpl A
NB.   P=. p2P p
NB.   ip=. /: p
NB.   iP=. p2P ip
NB. then
NB.   P -: %. iP
NB.   P -: |: iP
NB.   ((-:!.(2^_34)) ((/: @ (0&({::))) pp ((2&({::)) (] mp (mp ct)) (1&({::))))@hetrfpl) A
NB.   ((-:!.(2^_34)) ((ip2P @ (0&({::))) (mp mp (ct@[)) ((2&({::)) (] mp (mp ct)) (1&({::))))@hetrfpl) A
NB.
NB. Notes:
NB. - non-blocked Aasen's method is used
NB.
NB. TODO:
NB. - T should be sparse

hetrfpl=: 3 : 0
  n=. # y
  L1=. ($ y) {. 1              NB. 1st L1's column in Aasen method is e(1)-vector
  T=. ($ {. (1 1 {. ])) y     NB. T[0,0]=A[0,0]
  ip=. i. n
  h=. i. 1                    NB. h[0:j-1]=H[0:j-1,j-1], may be defined arbitrary before the 1st iteration only
  ios=. 1 }. i. n             NB. j:n-1
  for_j. }. i. n do.          NB. 1:n-1
    a=. (< ios ; (<: j)) { y  NB. A[j:n-1,j-1]
    lum=. (j (- , [) n) {. L1  NB. L1[j:n-1,0:j-1]
    v=. a - lum mp h          NB. v[0:n-j-1]=A[j:n-1,j-1]-L1[j:n-1,0:j-1]*H[0:j-1,j-1]=L1[j:n-1,j]*T[j,j-1] non-pivoted yet
    q=. iofmaxm v             NB. IO pivot from head
    v=. (0 ios2cp q) C. v     NB. v[0]↔v[q]
    dip=. q (+ ios2cp ]) j    NB. any[j]↔any[j+q]
    y=. dip pp y              NB. A[j,j:n-1]↔A[j+q,j:n-1], A[j:n-1,j]↔A[j:n-1,j+q]
    ip=. dip C. ip            NB. ip[j]↔ip[j+q]
    to=. {. v                 NB. T[j,j-1]
    lu=. q { lum              NB. L1[j,0:j-1] after pivoting
    luecto=. to * + {: lu     NB. conj(L1[j,j-1])*T[j,j-1]
    h=. to (((+ +)~ ({~ _1:))`_1:`]) } ((2 # j) {. T) mp + lu  NB. h[0:j-1]=H[0:j-1,j]=T[0:j-1,0:j-1]*conj(L1[j,0:j-1])
    L1=. (1 0:} v % to) (< ios ; j) } dip C. L1                  NB. L1[j,0:n-1]↔L1[j+q,0:n-1], L1[j:n-1,j]=v[0:n-j-1]/v[0], v[0] may be 0
    td=. 9 o. ((< 2 # j) { y) - (luecto + lu mp h)             NB. T[j,j]=Re(A[j,j]-L1[j,0:j-1]*H[0:j-1,j]-conj(L1[j,j-1])*T[j,j-1])
    T=. (td (,~ (, +)) to) (_2 <\ 0 _1 _1 0 0 0 + j) } T     NB. batch write diagonal and off-diagonal elements T[j,j-1] T[j-1,j] T[j,j]
    h=. h , luecto + td       NB. h[0:j]=H[0:j,j]=T[0:J,0:J]*conj(L1[j,0:j])
    ios=. }. ios              NB. j+1:n-1
  end.
  ip ; L1 ; T
)

NB. ---------------------------------------------------------
NB. hetrfpu
NB.
NB. Description:
NB.   Triangular factorization with full pivoting of a
NB.   Hermitian (symmetric) matrix:
NB.     P * U1 * T * U1' * P' = A
NB.
NB. Syntax:
NB.     'ip U1 T'=. hetrfpu A
NB. where
NB.   A   - n×n-matrix to factorize, Hermitian (symmetric)
NB.   ip  - n-vector, inversed full permutation of A
NB.   U1  - n×n-matrix, unit upper triangular
NB.   T   - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'ip U1 T'=. hetrfpu A
NB.   P=. p2P p
NB.   ip=. /: p
NB.   iP=. p2P ip
NB. then
NB.   P -: %. iP
NB.   P -: |: iP
NB.   ((-:!.(2^_34)) ((/: @ (0&({::))) pp ((2&({::)) (] mp (mp ct)) (1&({::))))@hetrfpu) A
NB.   ((-:!.(2^_34)) ((ip2P @ (0&({::))) (mp mp (ct@[)) ((2&({::)) (] mp (mp ct)) (1&({::))))@hetrfpu) A
NB.
NB. Notes:
NB. - non-blocked Aasen's method is used
NB.
NB. TODO:
NB. - T should be sparse

hetrfpu=: 3 : 0
  n=. # y
  U1=. (- $ y) {. 1                NB. last U1's column in Aasen method is e(n)-vector
  T=. ((- @ $) {. (_1 _1&{.)) y   NB. T[n-1,n-1]=A[n-1,n-1]
  ip=. i. n
  h=. i. 1                        NB. h[0:n-j-1]=H[j+1:n-1,j+1], may be defined arbitrary before the 1st iteration only
  ios=. }: i. n                   NB. 0:j
  for_j. }. i. - n do.            NB. n-2:0
    a=. (< ios ; (>: j)) { y      NB. A[0:j,j+1]
    lum=. ((>:j) ([ , -) n) {. U1  NB. U1[0:j,j+1:n-1]
    v=. a - lum mp h              NB. v[0:j]=A[0:j,j+1]-U1[0:j,j+1:n-1]*H[j+1:n-1,j+1]=U1[0:j,j]*T[j,j+1] non-pivoted yet
    q=. iolmaxm v                 NB. IO pivot from tail
    v=. (j ios2cp q) C. v         NB. v[_1]↔v[q]
    dip=. q ios2cp j              NB. any[j]↔any[q]
    y=. dip pp y                  NB. A[j,0:j]↔A[q,0:j], A[0:j,j]↔A[0:j,q]
    ip=. dip C. ip                NB. ip[j]↔ip[q]
    to=. {: v                     NB. T[j,j+1]
    lu=. q { lum                  NB. U1[j,j+1:n-1] after pivoting
    luecto=. to * + {. lu         NB. conj(U1[j,j+1])*T[j,j+1]
    h=. to (((+ +)~ ({~ 0:))`0:`]) } ((2 # >: j - n) {. T) mp + lu  NB. h[0:n-j-1]=H[j+1:n-1,j]=T[j+1:n-1,j+1:n-1]*conj(U1[j,j+1:n-1])
    U1=. (1 _1:} v % to) (< ios ; j) } dip C. U1                      NB. U1[j,0:n-1]↔U1[q,0:n-1], U1[0:j,j]=v[0:j]/v[_1], v[_1] may be 0
    td=. 9 o. ((< 2 $ j) { y) - (luecto + lu mp h)                  NB. T[j,j]=Re(A[j,j]-U1[j,j+1:n-1]*H[j+1:n-1,j]-conj(U1[j,j+1])*T[j,j+1])
    T=. (td (,~ (, +)) to) (_2 <\ 0 1 1 0 0 0 + j) } T          NB. batch write diagonal and off-diagonal elements T[j,j+1] T[j+1,j] T[j,j]
    h=. h ,~ luecto + td          NB. h[0:j]=H[0:j,j]=T[0:J,0:J]*conj(U1[j,0:j])
    ios=. }: ios                  NB. 0:j-1
  end.
  ip ; U1 ; T
)

NB. ---------------------------------------------------------
NB. potrfl
NB.
NB. Description:
NB.   Cholesky factorization of a Hermitian (symmetric)
NB.   positive definite matrix:
NB.     L * L' = A
NB.
NB. Syntax:
NB.   L=. potrfl A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.   L - n×n-matrix, lower triangular Cholesky factor
NB.
NB. If:
NB.   L=. potrfl A
NB. then
NB.   ((-:!.(2^_34)) (mp ct)@potrfl) A

potrfl=: 3 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    Ta=. potrfl (2 # p) {. y     NB. recursively factorize square matrix from top left or bottom right corner
    Ac=. (p ([ , -) n) {. y      NB. off-diagonal part of input matrix
    Tb=. ct Tbh=. Ta trtrslx Ac  NB. off-diagonal part of output matrix
    Aa=. (2 # p) }. y            NB. square matrix from opposite corner on diagonal of input matrix
    Tc=. potrfl Aa - Tb mp Tbh   NB. recursively factorize square matrix from opposite corner on diagonal
    Ta 0 append (Tb ,. Tc)       NB. assemble output as triangular matrix
  else.
    %: y
  end.
)

NB. ---------------------------------------------------------
NB. potrfu
NB.
NB. Description:
NB.   Cholesky factorization of a Hermitian (symmetric)
NB.   positive definite matrix:
NB.     U * U' = A
NB.
NB. Syntax:
NB.   U=. potrfu A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.   U - n×n-matrix, upper triangular Cholesky factor
NB.
NB. If:
NB.   U=. potrfu A
NB. then
NB.   ((-:!.(2^_34)) (mp ct)@potrfu) A

potrfu=: 3 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    Ta=. potrfu (2 # p) }. y     NB. recursively factorize square matrix from top left or bottom right corner
    Ac=. (p ([ , -) n) }. y      NB. off-diagonal part of input matrix
    Tb=. ct Tbh=. Ta trtrsux Ac  NB. off-diagonal part of output matrix
    Aa=. (2 # p) {. y            NB. square matrix from opposite corner on diagonal of input matrix
    Tc=. potrfu Aa - Tb mp Tbh   NB. recursively factorize square matrix from opposite corner on diagonal
    (Tc ,. Tb) _1 append Ta      NB. assemble output as triangular matrix
  else.
    %: y
  end.
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. tgetrf
NB.
NB. Description:
NB.   Test triangular factorization algorithms
NB.   - lud (math/misc)
NB.   - getrf (math/lapack)
NB.   - getrflu1p getrfpl1u getrfpu1l getrful1p (math/mt)
NB.   by general matrix given
NB.
NB. Syntax:
NB.   tgetrf A
NB. where
NB.   A - m×n-matrix
NB.
NB. Formula:
NB. - L*U1*P=A : berr := ||L*U1*P - A ||/(ε*||A||*m)
NB. - P*L1*U=A : berr := ||P*L1*U - A ||/(ε*||A||*n)
NB. - P*U1*L=A : berr := ||P*U1*L - A ||/(ε*||A||*n)
NB. - U*L1*P=A : berr := ||U*L1*P - A ||/(ε*||A||*m)

tgetrf=: 3 : 0
  require '~addons/math/misc/matfacto.ijs'
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'getrf'

  rcond=. ((_."_)`(norm1 con getri) @. (=/@$)) y  NB. meaninigful for square matrices only

  ('getrf_jlapack_' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (((mp & >)/ @ }:) invperm_jlapack_                           (2 & {::) ))) % (FP_EPS*((norm1*c)@[)))))) y
  ('lud'            tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (((mp & >)/ @ }:) %.                                         (2 & {::) ))) % (FP_EPS*((norm1*c)@[)))))) y

  ('getrflu1p'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((/: @ (0 & {::)) C."1 (( trl            mp  tru1          )@(1 & {::))))) % (FP_EPS*((norm1*#)@[)))))) y
  ('getrfpl1u'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((/: @ (0 & {::)) C.   (( trl1           mp  tru           )@(1 & {::))))) % (FP_EPS*((norm1*c)@[)))))) y
  ('getrfpu1l'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((/: @ (0 & {::)) C.   (((tru1~ (-~/@$)) mp (trl ~ (-~/@$)))@(1 & {::))))) % (FP_EPS*((norm1*c)@[)))))) y
  ('getrful1p'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((/: @ (0 & {::)) C."1 (((tru ~ (-~/@$)) mp (trl1~ (-~/@$)))@(1 & {::))))) % (FP_EPS*((norm1*#)@[)))))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. thetrf
NB.
NB. Description:
NB.   Test triangular factorization algorithms
NB.   - hetrfpl hetrfpu (math/mt)
NB.   by Hermitian (symmetric) matrix given
NB.
NB. Syntax:
NB.   thetrf A
NB. where
NB.   A - n×n-matrix, Hermitian
NB.
NB. Formula:
NB. - P'*L*T*L'*P=A : berr := ||P'*L*T*L'*P - A ||/(ε*||A||*n)
NB. - P'*U*T*U'*P=A : berr := ||P'*U*T*U'*P - A ||/(ε*||A||*n)

thetrf=: 3 : 0
  rcond=. _.  NB. rcond=. (norm1 con hetri) y

  ('hetrfpl' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((mp mp ct@[)&>/@}. (pp~ /:) (0&({::)))))) % (FP_EPS*((norm1*c)@[))))) y
  ('hetrfpu' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((mp mp ct@[)&>/@}. (pp~ /:) (0&({::)))))) % (FP_EPS*((norm1*#)@[))))) y

EMPTY
)

NB. ---------------------------------------------------------
NB. tpotrf
NB.
NB. Description:
NB.   Test triangular factorization algorithms
NB.   - choleski (math/misc addon)
NB.   - potrf (math/lapack addon)
NB.   - potrfpl potrfpu (math/mt addon)
NB.   by Hermitian (symmetric) positive definite matrix given
NB.
NB. Syntax:
NB.   tpotrf A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive
NB.       definite
NB.
NB. Formula:
NB. - L*L'=A : berr := ||L*L' - A ||/(ε*||A||*n)
NB. - U*U'=A : berr := ||U*U' - A ||/(ε*||A||*n)

tpotrf=: 3 : 0
  require '~addons/math/misc/matfacto.ijs'
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'potrf'

  rcond=. _.  NB. rcond=. (norm1 con potri) y

  ('potrf_jlapack_' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (mp ct)))) % (FP_EPS*((norm1*c)@[))))) y
  ('choleski'       tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (mp ct)))) % (FP_EPS*((norm1*c)@[))))) y
  ('potrfl'         tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (mp ct)))) % (FP_EPS*((norm1*c)@[))))) y
  ('potrfu'         tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (mp ct)))) % (FP_EPS*((norm1*#)@[))))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testtrf
NB.
NB. Description:
NB.   Adv. to make verb to test TRF algorithms by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testtrf
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.            mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testtrf 150 100

testtrf=: 1 : 'EMPTY [ (((tpotrf @ (u pomat)) [ (thetrf @ (u hemat))) ^: (=/)) [ tgetrf @ u'
