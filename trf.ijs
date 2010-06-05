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
NB. hetrfu     UDU' factorization of a Hermitian (symmetric)
NB.            matrix, where U is unit upper triangular and
NB.            D is diagonal
NB. hetrfl     LDL' factorization of a Hermitian (symmetric)
NB.            matrix, where L is unit lower triangular and
NB.            D is diagonal
NB. ditrfu     ? factorization of a diagonalizable
NB.            matrix, where ?
NB. ditrfl     ? factorization of a diagonalizable
NB.            matrix, where ?
NB. potrfl     LL' (Cholesky) factorization of a Hermitian
NB.            (symmetric) positive definite matrix for lower
NB.            triangular factor
NB. potrfu     UU' (Cholesky) factorization of a Hermitian
NB.            (symmetric) positive definite matrix for upper
NB.            triangular factor
NB.
NB. Version: 1.0.0 2009-06-01
NB. Copyright: Igor Zhuravlov, igor at uic.dvgu.ru
NB. License: GNU GPL

coclass 'mt'

NB. =========================================================
NB. Local definitions

iomaxm=: (i.>./) @: |  NB. IO element with max magnitude from list y

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
NB. getrfl1u
NB. LU factorization of a general matrix with rows pivoting,
NB. where L is unit lower triangular and U is upper
NB. triangular:
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
NB.   L1=. trl1 L1U
NB.   U=. tru L1U
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
  'm n'=. sh=. $ y
  if. 0 e. sh do. (i. m) ; y return. end.
  if. n > 1 do.
    k=. m <. >. -: n
    'p LBC'=. getrfl1u (m , k) {. y                           NB. factorize left block column recursively
    y=. p C. (0 , k) }. y                                     NB. apply p1 to A's 2nd block column
    U12=. LBC (trtrsl1x & (k & {.)) y                         NB. solve L11*U12=A12 for U12
    if. k < m do.
      'p2 LU22'=. getrfl1u y ((- (mp & U12)) & (k & }.)) LBC  NB. factorize updated A22 recursively
      dp=. (i. k) , (k + p2)
      p=. dp C. p
      LBC=. dp C. LBC
      U12=. U12 , LU22
    end.
    p ; (LBC ,. U12)
  else.
    dp=. 0 ii2cp iomaxm y          NB. find pivot and form cycle permutation
    p=. dp C. i. m                 NB. adjust p by p2: ((p2 { (iosr { p)) iosr } p)
    y=. ({. 0 } (%"1 {.)) dp C. y  NB. apply rows permutation starting from j-th row, then scale by head
    p ; y
  end.
)

NB. ---------------------------------------------------------
NB. hetrfl
NB. Factorization of a  Hermitian (symmetric) matrix:
NB.    L * D * L' = A
NB. where L is unit lower triangular and D is diagonal

hetrfl=: 3 : 0
)

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

potrfu=: potrfu`}.`trtrsux`{.` ,.  `(_1 append ) potrf  NB. 
potrfl=: potrfl`{.`trtrslx`}.`(,.~)`( 0 append~) potrf  NB. 

NB. =========================================================
NB. Test suite

NB. 'name rcond bwerror fwerror time space'=. name vextract ttrf A
NB. TODO: forward error (getrf has two: for L and U?)
ttrf=: 1 : 0
:
  't s'=. timespacex 'out=. ' , x , ' y'
  be=. (((norm1 (y - u out)) % ({: $ y)) % (norm1 y)) % FP_EPS  NB. backward error
  dbprn x ; ((_."_)`(norm1 con getri) @. (=/@$) y) ; be ; (i.0) ; t ; s
)

NB. r=. tgetrf A
NB. where r is boxed table with columns: name rcond bwerror fwerror time space
tgetrf=: 3 : 0
  y=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y

  r=.  ,: 'getrfl1u'       (((C.~ /:)~ (trl1 mp tru)) & > /)                   ttrf y
  r=. r , 'getrf_jlapack_' (((mp & >)/ @ (2 & {.)) invperm_jlapack_ (2 & {::)) ttrf y
)

NB. r=. tpotrf A
NB. where r is boxed table with columns: name rcond bwerror fwerror time space
tpotrf=: 3 : 0
  y=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y

  r=.  ,: 'potrfu'         ((mp ct) @ tru) ttrf y
  r=. r , 'potrfl'         ((mp ct) @ trl) ttrf y
  r=. r , 'potrf_jlapack_' ((mp ct) @ trl) ttrf y
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
NB.   r    - boxed table with 6 columns:
NB.          - algorithm name
NB.          - the estimated reciprocal of the condition
NB.            number of a matrix in 1-norm
NB.          - relative backward error
NB.          - relative forward error
NB.          - execution time, sec.
NB.          - execution space, bytes
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   cocurrent 'mt'
NB.   r=. (_1 1 0 16 _6 4 & gemat) testtrf 500 500
NB.   r=. (_1 1 0 16 _6 4 & (gemat j. gemat)) testtrf 500 500

testtrf=: 1 : 0
  r=.     (tgetrf @ u)                      y
  r=. r , (tpotrf @ (u pomat) @ {. ^: (=/)) y
)
