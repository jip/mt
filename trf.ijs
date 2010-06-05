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
NB.            is permutation matrix, U is upper
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

NB. ---------------------------------------------------------
NB. hetf2pl
NB.
NB. Description:
NB.   Triangular factorization with full pivoting of a
NB.   Hermitian (symmetric) matrix:
NB.     P * L1 * T * L1' * P' = A
NB.    by non-blocked version of Aasen's algorithm
NB.
NB. Syntax:
NB.     'ip L1 T'=. hetf2pl A
NB. where
NB.   A   - n×n-matrix to factorize, Hermitian (symmetric)
NB.   ip  - n-vector, full inversed permutation of A
NB.   L1  - n×n-matrix, unit lower triangular
NB.   T   - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   k   = min(m,n)
NB.
NB. Assertions:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A (-:!.(2^_34)) p sp T (] mp (mp ct)) L1
NB.   A (-:!.(2^_34)) P (mp mp |:@[) T (] mp (mp ct)) L1
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
NB. - T should be sparse

hetf2pl=: 3 : 0
  n=. # y
  L1=. ($ y) {. 1              NB. 1st L1's column in Aasen method is e(1)-vector
  T=. ($ {. (1 1 {. ])) y      NB. T[0,0]=A[0,0]
  ip=. i. n
  h=. i. 1                     NB. h[0:j-1]=H[0:j-1,j-1], may be defined arbitrary before the 1st iteration only
  ios=. n ht2lios 1            NB. j:n-1
  for_j. n ht2lios 1 do.       NB. 1:n-1
    a=. (< ios ; (<: j)) { y   NB. A[j:n-1,j-1]
    lum=. (j (- , [) n) {. L1  NB. L1[j:n-1,0:j-1]
    v=. a - lum mp h           NB. v[0:n-j-1]=A[j:n-1,j-1]-L1[j:n-1,0:j-1]*H[0:j-1,j-1]=L1[j:n-1,j]*T[j,j-1] non-pivoted yet
    q=. iofmax v               NB. IO pivot from head
    v=. (0 ios2cp q) C. v      NB. v[0]↔v[q]
    dip=. q (+ ios2cp ]) j     NB. any[j]↔any[j+q]
    y=. dip sp y               NB. A[j,j:n-1]↔A[j+q,j:n-1], A[j:n-1,j]↔A[j:n-1,j+q]
    ip=. dip C. ip             NB. ip[j]↔ip[j+q]
    to=. {. v                  NB. T[j,j-1]
    lu=. q { lum               NB. L1[j,0:j-1] after pivoting
    luecto=. to * + {: lu      NB. conj(L1[j,j-1])*T[j,j-1]
    h=. to (((+ +)~ {:)`_1:`]) } ((2 # j) {. T) mp + lu   NB. h[0:j-1]=H[0:j-1,j]=T[0:j-1,0:j-1]*conj(L1[j,0:j-1])
    L1=. (1 (0}) v % to) (< ios ; j) } dip C. L1          NB. L1[j,0:n-1]↔L1[j+q,0:n-1], L1[j:n-1,j]=v[0:n-j-1]/v[0], v[0] may be 0
    td=. 9 o. ((< 2 # j) { y) - (luecto + lu mp h)        NB. T[j,j]=Re(A[j,j]-L1[j,0:j-1]*H[0:j-1,j]-conj(L1[j,j-1])*T[j,j-1])
    T=. (td (,~ (, +)) to) (_2 <\ 0 _1 _1 0 0 0 + j) } T  NB. batch write diagonal and off-diagonal elements T[j,j-1] T[j-1,j] T[j,j]
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
NB.     P * U1 * T * U1' * P' = A
NB.    by non-blocked version of Aasen's algorithm
NB.
NB. Syntax:
NB.     'ip U1 T'=. hetf2pu A
NB. where
NB.   A   - n×n-matrix to factorize, Hermitian (symmetric)
NB.   ip  - n-vector, full inversed permutation of A
NB.   U1  - n×n-matrix, unit upper triangular
NB.   T   - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   k   = min(m,n)
NB.
NB. Assertions:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A (-:!.(2^_34)) p sp T (] mp (mp ct)) U1
NB.   A (-:!.(2^_34)) P (mp mp |:@[) T (] mp (mp ct)) U1
NB. where
NB.   'ip U1 T'=. hetf2pu A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.
NB. TODO:
NB. - T should be sparse

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
    q=. iolmax v                   NB. IO pivot from tail
    v=. (j ios2cp q) C. v          NB. v[_1]↔v[q]
    dip=. q ios2cp j               NB. any[j]↔any[q]
    y=. dip sp y                   NB. A[j,0:j]↔A[q,0:j], A[0:j,j]↔A[0:j,q]
    ip=. dip C. ip                 NB. ip[j]↔ip[q]
    to=. {: v                      NB. T[j,j+1]
    lu=. q { lum                   NB. U1[j,j+1:n-1] after pivoting
    luecto=. to * + {. lu          NB. conj(U1[j,j+1])*T[j,j+1]
    h=. to (((+ +)~ {.)`0:`]) } ((2 # >: j - n) {. T) mp + lu  NB. h[0:n-j-1]=H[j+1:n-1,j]=T[j+1:n-1,j+1:n-1]*conj(U1[j,j+1:n-1])
    U1=. (1 (_1}) v % to) (< ios ; j) } dip C. U1              NB. U1[j,0:n-1]↔U1[q,0:n-1], U1[0:j,j]=v[0:j]/v[_1], v[_1] may be 0
    td=. 9 o. ((< 2 # j) { y) - (luecto + lu mp h)             NB. T[j,j]=Re(A[j,j]-U1[j,j+1:n-1]*H[j+1:n-1,j]-conj(U1[j,j+1])*T[j,j+1])
    T=. (td (,~ (, +)) to) (_2 <\ 0 1 1 0 0 0 + j) } T         NB. batch write diagonal and off-diagonal elements T[j,j+1] T[j+1,j] T[j,j]
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
NB.   L   - m×k-matrix, lower triangular
NB.   U1  - k×n-matrix, unit upper triangular
NB.   k   = min(m,n)
NB.
NB. Assertions:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A (-:!.(2^_34)) p {"1 L mp U1
NB.   A (-:!.(2^_34)) p C."1 L mp U1
NB.   A (-:!.(2^_34)) ip C.^:_1"1 L mp U1
NB.   A (-:!.(2^_34)) L mp U1 mp P
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
    dip=. 0 ios2cp iofmax {. y
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
NB.   L1  - m×k-matrix, unit lower triangular
NB.   U   - k×n-matrix, upper triangular
NB.   k   = min(m,n)
NB.
NB. Assertions:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A (-:!.(2^_34)) p { L1 mp U
NB.   A (-:!.(2^_34)) p C. L1 mp U
NB.   A (-:!.(2^_34)) ip C.^:_1 L1 mp U
NB.   A (-:!.(2^_34)) P mp L1 mp U
NB. where
NB.   'ip L1U'=. getrfpl1u A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.   L1=. trl1 L1U
NB.   U=. tru L1U

getrfpl1u=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. m) ; y
  elseif. 1 = n do.
    dip=. 0 ios2cp iofmax y
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
NB.   L   - m×k-matrix, lower triangular
NB.   U1  - k×n-matrix, unit upper triangular
NB.   k   = min(m,n)
NB.
NB. Assertions:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A (-:!.(2^_34)) p { U1 mp L
NB.   A (-:!.(2^_34)) p C. U1 mp L
NB.   A (-:!.(2^_34)) ip C.^:_1 U1 mp L
NB.   A (-:!.(2^_34)) P mp U1 mp L
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
    dip=. (<:m) ios2cp iolmax y
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
NB.   L   - m×k-matrix, unit lower triangular
NB.   U1  - k×n-matrix, upper triangular
NB.   k   = min(m,n)
NB.
NB. Assertions:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A (-:!.(2^_34)) p {"1 U mp L1
NB.   A (-:!.(2^_34)) p C."1 U mp L1
NB.   A (-:!.(2^_34)) ip C.^:_1"1 U mp L1
NB.   A (-:!.(2^_34)) U mp L1 mp iP
NB. where
NB.   'p UL1'=. getrful1p A
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
    dip=. (<:n) ios2cp iolmax {. y
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
NB.     P * L1 * T * L1' * P' = A
NB.    by blocked version of Aasen's algorithm
NB.
NB. Syntax:
NB.     'ip L1 T'=. hetrfpl A
NB. where
NB.   A   - n×n-matrix to factorize, Hermitian (symmetric)
NB.   ip  - n-vector, full inversed permutation of A
NB.   L1  - n×n-matrix, unit lower triangular
NB.   T   - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   k   = min(m,n)
NB.
NB. Assertions:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A (-:!.(2^_34)) p sp T (] mp (mp ct)) L1
NB.   A (-:!.(2^_34)) P (mp mp |:@[) T (] mp (mp ct)) L1
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
NB. - T should be sparse

hetrfpl=: hetf2pl

NB. ---------------------------------------------------------
NB. hetrfpu
NB.
NB. Description:
NB.   Triangular factorization with full pivoting of a
NB.   Hermitian (symmetric) matrix:
NB.     P * U1 * T * U1' * P' = A
NB.    by blocked version of Aasen's algorithm
NB.
NB. Syntax:
NB.     'ip U1 T'=. hetrfpu A
NB. where
NB.   A   - n×n-matrix to factorize, Hermitian (symmetric)
NB.   ip  - n-vector, full inversed permutation of A
NB.   U1  - n×n-matrix, unit upper triangular
NB.   T   - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   k   = min(m,n)
NB.
NB. Assertions:
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A (-:!.(2^_34)) p sp T (] mp (mp ct)) U1
NB.   A (-:!.(2^_34)) P (mp mp |:@[) T (] mp (mp ct)) U1
NB. where
NB.   'ip U1 T'=. hetrfpu A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.
NB. TODO:
NB. - implement partitioned algorithm
NB. - T should be sparse

hetrfpu=: hetf2pu

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
NB.   L - n×n-matrix, lower triangular with positive diagonal
NB.       entries, Cholesky triangle
NB.
NB. Assertions:
NB.   A (-:!.(2^_34)) (mp ct) L
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
NB.   U - n×n-matrix, upper triangular with positive diagonal
NB.       entries, Cholesky triangle
NB.
NB. Assertions:
NB.   A (-:!.(2^_34)) (mp ct) U
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

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgetrf
NB.
NB. Description:
NB.   Test triangular factorization algorithms:
NB.   - lud (math/misc)
NB.   - getrf (math/lapack)
NB.   - getrflu1p getrfpl1u getrfpu1l getrful1p (math/mt)
NB.   by general matrix given
NB.
NB. Syntax:
NB.   testgetrf A
NB. where
NB.   A - m×n-matrix
NB.
NB. Formula:
NB. - L*U1*P=A : berr := ||L*U1*P - A ||/(ε*||A||*m)
NB. - P*L1*U=A : berr := ||P*L1*U - A ||/(ε*||A||*n)
NB. - P*U1*L=A : berr := ||P*U1*L - A ||/(ε*||A||*n)
NB. - U*L1*P=A : berr := ||U*L1*P - A ||/(ε*||A||*m)

testgetrf=: 3 : 0
  require '~addons/math/misc/matfacto.ijs'
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'getrf'

  rcond=. ((_."_)`(norm1 con getri) @. (=/@$)) y  NB. meaninigful for square matrices only

  ('getrf_jlapack_' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (((mp & >)/ @ }:) invperm_jlapack_                           (2 & {::) ))) % (FP_EPS*((norm1*c)@[)))))) y
  ('lud'            tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (((mp & >)/ @ }:) %.                                         (2 & {::) ))) % (FP_EPS*((norm1*c)@[)))))) y

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
NB.   Test triangular factorization algorithms:
NB.   - hetrfpl hetrfpu (math/mt)
NB.   by Hermitian (symmetric) matrix given
NB.
NB. Syntax:
NB.   testhetrf A
NB. where
NB.   A - n×n-matrix, Hermitian
NB.
NB. Formula:
NB. - P*L*T*L'*P'=A : berr := ||P*L*T*L'*P' - A ||/(ε*||A||*n)
NB. - P*U*T*U'*P'=A : berr := ||P*U*T*U'*P' - A ||/(ε*||A||*n)

testhetrf=: 3 : 0
  rcond=. _.  NB. rcond=. (norm1 con hetri) y

  ('hetrfpl' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((mp mp ct@[)&>/@}. (sp~ /:) (0&({::)))))) % (FP_EPS*((norm1*c)@[))))) y
  ('hetrfpu' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((mp mp ct@[)&>/@}. (sp~ /:) (0&({::)))))) % (FP_EPS*((norm1*#)@[))))) y

EMPTY
)

NB. ---------------------------------------------------------
NB. testpotrf
NB.
NB. Description:
NB.   Test triangular factorization algorithms:
NB.   - choleski (math/misc addon)
NB.   - potrf (math/lapack addon)
NB.   - potrfpl potrfpu (math/mt addon)
NB.   by Hermitian (symmetric) positive definite matrix given
NB.
NB. Syntax:
NB.   testpotrf A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive
NB.       definite
NB.
NB. Formula:
NB. - L*L'=A : berr := ||L*L' - A ||/(ε*||A||*n)
NB. - U*U'=A : berr := ||U*U' - A ||/(ε*||A||*n)

testpotrf=: 3 : 0
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

testtrf=: 1 : 'EMPTY_mt_ [ (((testpotrf_mt_ @ (u pomat_mt_)) [ (testhetrf_mt_ @ (u hemat_mt_))) ^: (=/)) [ testgetrf_mt_ @ u'
