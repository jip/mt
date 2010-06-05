NB. trf.ijs
NB. Triangular factorization of matrix
NB.
NB. getrfl1u   LU factorization of a general matrix, where
NB.            L is unit lower triangular and U is upper
NB.            triangular
NB. getrflu1   LU factorization of a general matrix, where
NB.            L is lower triangular and U is unit upper
NB.            triangular
NB. getrfu1l   UL factorization of a general matrix, where
NB.            U is unit upper triangular and L is lower
NB.            triangular
NB. getrful1   UL factorization of a general matrix, where
NB.            U is upper triangular and L is unit lower
NB.            triangular
NB. potrfl     Cholesky factorization of a Hermitian
NB.            (symmetric) positive definite matrix for lower
NB.            triangular factor
NB. potrfu     Cholesky factorization of a Hermitian
NB.            (symmetric) positive definite matrix for upper
NB.            triangular factor
NB.
NB. Version: 1.0.0 2008-11-01
NB. Copyright: Igor Zhuravlov, igor at uic.dvgu.ru
NB. License: GNU GPL

coclass 'mt'

NB. =========================================================
NB. Local definitions

iomaxm=: (i.>./) @: |  NB. IO element with max magnitude from list y

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. getrfl1u
NB. LU factorization of a general matrix with rows pivoting,
NB. where L is unit lower triangular and U is upper
NB. triangular
NB.   P * L * U = A
NB.
NB. Syntax:
NB.   'p L1U'=. getrfl1u A
NB. where
NB.   A   - m×n-matrix
NB.   p  - m-vector, inversed rows permutation of A
NB.   L1U - m×n-matrix, strict lower triangle keeps L (unit
NB.         diagonal elements are not saved), upper triangle
NB.         keeps U
NB.   L   - m×min(m,n)-matrix, unit lower triangular Cholesky
NB.         factor
NB.   U   - min(m,n)×n-matrix, upper triangular Cholesky
NB.         factor
NB.   m   ≥ 0
NB.   n   ≥ 0
NB.
NB. If:
NB.   'pi L1U'=. getrfl1u A
NB.   p=. /: pi
NB.   Pi=. ((i.#pi)=/pi)
NB.   P=. ((i.#p)=/p)
NB.   L1=. ltri1 L1U
NB.   U=. utri L1U
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

getrfl1u=: 3 : 0
  'm n'=. mn=. $ y
  if. n > 1 do.
    k=. m <. n2=. >. -: n
    'p1 LU1'=. getrfl1u (m , k) {. y                        NB. factorize A's 1st block column recursively
    y=. p1 C. (0 , k) }. y                                  NB. apply p1 to A's 2nd block column
    U12=. y (getrsl1x & (k & {.)) LU1                       NB. solve L11*U12=A12 for U12
    'p2 LU22'=. getrfl1u y ((- (mp & U12)) & (k & }.)) LU1  NB. factorize updated A22 recursively
    dp=. (i. k) , (k + p2)
    p=. dp C. p1
    p ; (dp C. LU1) ,. (U12 , LU22)
  else.
    dp=. 0 ii2cp iomaxm y                       NB. find pivot and form cycle permutation
    p=. dp C. i. m                              NB. adjust p by p2: ((p2 { (iosr { p)) iosr } p)
    y=. ({. 0 } (%"1 {.)) dp C. y               NB. apply rows permutation starting from j-th row, then scale by head
    p ; y
  end.
)

NB. ---------------------------------------------------------
NB. potrfl
NB. Hermitian (symmetric) positive definite matrix Cholesky
NB. recursive factorization by a lower triangular factor
NB.   L * L' = A
NB.
NB. Syntax:
NB.   L=. potrfl A
NB. where
NB.   A - N×N-matrix, Hermitian (symmetric), positive
NB.       definite
NB.   L  - N×N-matrix, lower triangular Cholesky factor
NB.
NB. If:
NB.   L=. potrfl A
NB. then
NB.   A -: (mp ct) L
NB. where
NB.   A - positive definite matrix
NB.   L - lower triangular matrix
NB.
NB. Notes:
NB.   - FLOPs: 
NB.   - recursive calls: 3*k*(2^k)-2^(k+1)+3,
NB.     where k = ⌈log_2(n)⌉

potrfl=: 3 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    q=. n - p
    L0=. potrfl (2 $ p) {. y
    Y=. (p , -q) {. y
    L1=. ct L1h=. Y getrslx L0
    Z=. (2 $ p) }. y
    L2=. potrfl Z - L1 mp L1h
    L0 , L1 ,. L2
  else.
    %: y
  end.
)

NB. ---------------------------------------------------------
NB. potrfu
NB. Hermitian (symmetric) positive definite matrix Cholesky
NB. recursive factorization by a upper triangular factor
NB.   U * U' = A
NB.
NB. Syntax:
NB.   U=. potrfu A
NB. where
NB.   A - N×N-matrix, Hermitian (symmetric), positive
NB.       definite
NB.   U  - N×N-matrix, upper triangular Cholesky factor
NB.
NB. If:
NB.   U=. potrfl A
NB. then
NB.   A -: (mp ct) U
NB. where
NB.   A - positive definite matrix
NB.   U - upper triangular matrix
NB.
NB. Notes:
NB.   - FLOPs: 
NB.   - recursive calls: 3*k*(2^k)-2^(k+1)+3,
NB.     where k  = ⌈log_2(⌈n/nb⌉)⌉
NB.           nb = CPU_CACHE, block size

potrfu=: 3 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    q=. n - p
    U2=. potrfu (2 $ p) }. y
    Yh=. (p , -q) }. y
    U1=. ct U1h=. Yh getrsux U2
    X=. (2 $ p) {. y
    U0=. potrfu X - U1 mp U1h
    (n {. U0) ,. (U1 , U2)
  else.
    %: y
  end.
)

NB. =========================================================
NB. Test suite

NB. 'name error time space'=. name vextract ttrf A
ttrf=: 1 : 0
:
  't s'=. timespacex 'out=. ' , x , ' y'
  be=. (((norm1 (y - u out)) % ({: $ y)) % (norm1 y)) % FP_EPS  NB. backward error
  dbprn x ; be ; (i.0) ; t ; s
)

NB. r=. tgetrf A
NB. where r is boxed table with columns: name error time space
tgetrf=: 3 : 0
  r=.  ,: 'getrfl1u'       (((C.~ /:)~ (ltri1 mp utri)) & > /)                 ttrf y
  r=. r , 'getrf_jlapack_' (((mp & >)/ @ (2 & {.)) invperm_jlapack_ (2 & {::)) ttrf y
)

NB. r=. tpotrf A
NB. where r is boxed table with columns: name error time space
tpotrf=: 3 : 0
  r=.  ,: 'potrfu'         ((mp ct) @ utri) ttrf y
  r=. r , 'potrfl'         ((mp ct) @ ltri) ttrf y
  r=. r , 'potrf_jlapack_' ((mp ct) @ ltri) ttrf y
)

NB. ---------------------------------------------------------
NB. testtrf
NB. Adverb to test triangular factorization algorithms
NB.
NB. Syntax:
NB.   r=. mkge testtrf m,n
NB. where
NB.   m,n  - 2-vector of integers, shape of random matrices
NB.          to test algorithms; if m≠n then algorithms that
NB.          accept square matrices only are skipped
NB.   mkge - monadic verb to generate random non-singular
NB.          general y-matrix (shape is taken from y)
NB.   r    - boxed table with 4 columns: 'algorithm name'
NB.          'error' 'time, sec.' 'space, bytes'
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   r=. (_1 1 0 26 _16 8 & gemat) testtrf 7 7

testtrf=: 1 : 0
  r=.     (tpotrf @ (u pomat) @ {. ^: (=/)) y     NB. FIXME!
  r=. r , (tgetrf @ u)                      y
)
