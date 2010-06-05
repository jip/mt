NB. sv.ijs
NB. Solve linear system
NB.
NB. getrsu1x Solve system U*X=B, where U is an unit upper
NB.          triangular matrix
NB. getrsux  Solve system U*X=B, where U is an upper
NB.          triangular matrix
NB. getrsxu1 Solve system X*U=B, where U is an unit upper
NB.          triangular matrix
NB. getrsxu  Solve system X*U=B, where U is an upper
NB.          triangular matrix
NB. getrsl1x Solve system L*X=B, where L is an unit lower
NB.          triangular matrix
NB. getrslx  Solve system L*X=B, where L is a lower
NB.          triangular matrix
NB. getrsxl1 Solve system X*L=B, where L is an unit lower
NB.          triangular matrix
NB. getrsxl  Solve system X*L=B, where L is a lower
NB.          triangular matrix
NB. gesv     Solve system A*X=B via LU factorization with
NB.          partial pivoting
NB. disv     Solve system A*x=b, where A is a diagonalizable
NB.          matrix
NB. hesv     Solve system A*x=b, where A is a Hermitian
NB.          (symmetric) matrix
NB. posv     Solve system A*x=b, where A is a Hermitian
NB.          (symmetric) positive definite matrix

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. getrsu1x
NB. getrsux
NB. getrsxu1
NB. getrsxu
NB. getrsl1x
NB. getrslx
NB. getrsxl1
NB. getrsxl
NB. Solve systems: U1*X=B, U*X=B, X*U1=B, X*U=B, L1*X=B,
NB. L*X=B, X*L1=B, X*L=B; where U is an upper triangular
NB. matrix, L is an lower triangular matrix, U1 is an
NB. upper triangular matrix with units on diagonal, L1 is an
NB. lower triangular matrix with units on diagonal.
NB.
NB. Syntax:
NB.   X=. B getrsu1x LU1
NB.   X=. B getrsux L1U
NB.   X=. B getrsxu1 LU1
NB.   X=. B getrsxu L1U
NB.   X=. B getrsl1x L1U
NB.   X=. B getrslx LU1
NB.   X=. B getrsxl1 L1U
NB.   X=. B getrsxl LU1
NB. where <<<<<<<<<<<<<<<<<<<<<<<<<<
NB.   LU - N×N-matrix, the factors L and U from
NB.        factorization
NB.   B  - N-vector or N×NRHS-matrix, RHS
NB.   X  - same shape as B, solution
NB.   N    >= 0
NB.   NRHS >= 0
NB.
NB. If:
NB.   
NB. then
NB.   
NB.
NB. Notes:
NB. - strict lower triangle is not referenced
NB. - reuse x to save result
NB. - result is identical to LAPACK's dgetrs/zgetrs
NB. - mp involves excessive pair A[i][i]*z[i] (z[i]==0) to avoid increment
NB.
NB. TODO:
NB. - 

getrsu1x=: 4 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    U22=. (2 $ p) }. y
    B2=. p }. x
    X2=. B2 getrsu1x U22
    U12=. (p , (p - n)) {. y
    U11=. (2 $ p) {. y
    B1=. p {. x
    ((B1 - U12 mp X2) getrsu1x U11) , X2
  else.
    x
  end.
)

getrsux=: 4 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    U22=. (2 $ p) }. y
    B2=. p }. x
    X2=. B2 getrsux U22
    U12=. (p , (p - n)) {. y
    U11=. (2 $ p) {. y
    B1=. p {. x
    ((B1 - U12 mp X2) getrsux U11) , X2
  else.
    x % 0 ({,) y
  end.
)

getrsxu1=: 4 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    U11=. (2 $ p) {. y
    B1=. p {. x
    X1=. B1 getrsxu1 U11
    U12=. (p , (p - n)) {. y
    U22=. (2 $ p) }. y
    B2=. p }. x
    X1 , (B2 - X1 mp U12) getrsxu1 U22
  else.
    x
  end.
)

getrsxu=: 4 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    U11=. (2 $ p) {. y
    B1=. p {. x
    X1=. B1 getrsxu U11
    U12=. (p , (p - n)) {. y
    U22=. (2 $ p) }. y
    B2=. p }. x
    X1 , (B2 - X1 mp U12) getrsxu U22
  else.
    x % 0 ({,) y
  end.
)

getrsl1x=: 4 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    L11=. (2 $ p) {. y
    B1=. p {. x
    X1=. B1 getrsl1x L11
    L21=. ((p - n) , p) {. y
    L22=. (2 $ p) }. y
    B2=. p }. x
    X1 , (B2 - L21 mp X1) getrsl1x L22
  else.
    x
  end.
)

getrslx=: 4 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    L11=. (2 $ p) {. y
    B1=. p {. x
    X1=. B1 getrslx L11
    L21=. ((p - n) , p) {. y
    L22=. (2 $ p) }. y
    B2=. p }. x
    X1 , (B2 - L21 mp X1) getrslx L22
  else.
    x % 0 ({,) y
  end.
)

getrsxl1=: 4 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    L22=. (2 $ p) }. y
    B2=. p }. x
    X2=. B2 getrsxl1 L22
    L21=. ((p - n) , p) {. y
    L11=. (2 $ p) {. y
    B1=. p {. x
    ((B1 - X2 mp L21) getrsxl1 L11) , X2
  else.
    x
  end.
)

getrsxl=: 4 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    L22=. (2 $ p) }. y
    B2=. p }. x
    X2=. B2 getrsxl L22
    L21=. ((p - n) , p) {. y
    L11=. (2 $ p) {. y
    B1=. p {. x
    ((B1 - X2 mp L21) getrsxl L11) , X2
  else.
    x % 0 ({,) y
  end.
)

NB. ---------------------------------------------------------
NB. gesv                                                  1 2
NB. Solve system A*X=B, where A is a general matrix
NB.
NB. Syntax:
NB.   X=. B gesv A
NB. where
NB.   A - N×N-matrix
NB.   B  - N-vector or N×NRHS-matrix, RHS
NB.   X  - same shape as B, solution
NB.   N    >= 0
NB.   NRHS >= 0
NB.
NB. If:
NB.   X=. B gesv A
NB. then
NB.   B -: A mp X
NB.
NB. Notes:
NB. - result is identical to LAPACK's dgesv/zgesv

gesv=: (((0 {:: ]) C. [) (getrsl1x getrsux ]) (1 {:: ])) getrfl1u

NB. ---------------------------------------------------------
NB. disv                                                  1 1
NB. Solve system A*x=b, where A is a diagonalizable matrix
NB.
NB. Syntax:
NB.   x=. b disv (rv ; ev ; rvi)
NB. where
NB.   rv  - N×N-matrix, right eigenvectors of A
NB.   ev  - N-vector, eigenvalues of A
NB.   rvi - N×N-matrix, inversion of rv
NB.   x   - N-vector, solution
NB.   N  >= 0

disv=: [: : ((0 {:: ]) % (1 {:: ])) mp (mp~ (2 & {::)) " 1 1

NB. ---------------------------------------------------------
NB. hesv                                                  1 1
NB. Solve system A*x=b, where A is a Hermitian matrix
NB.
NB. Syntax:
NB.   x=. b hesv (rv ; ev)
NB. where
NB.   rv - N×N-matrix, right eigenvectors of A
NB.   ev - N-vector, eigenvalues of A
NB.   x  - N-vector, solution
NB.   N >= 0

hesv=: [: : ((0 {:: ]) % (1 {:: ])) mp (mp~ (+ @ |: @ (0 & {::))) " 1 1

NB. ---------------------------------------------------------
NB. posv
NB. Solve system A*x=b, where A is a Hermitian (symmetric)
NB. positive definite matrix
NB.
NB. Syntax:
NB.   X=. B posv A
NB. where
NB.   A    - n×n-matrix, right eigenvectors of A
NB.   B    - n-vector or n×nrhs-matrix, RHS
NB.   X    - same shape as B, solution
NB.   n    ≥ 0
NB.   nrhs ≥ 0

posv=: [:

NB. =========================================================
NB. Test suite

NB. 'name errors time space'=. name vrcond ttrs A;X;B
ttrs=: 1 : 0
:
  'A X B'=. y
  't s'=. timespacex 'Xsol=. B ' , x , ' A'
  be=. (((norm1t (B - A mp Xsol)) % (norm1 A)) % (norm1t Xsol)) % FP_EPS  NB. backward error
  fe=. (((normit (X - Xsol)) * (u A)) % (normit X)) % FP_EPS              NB. forward error
  dbprn x ; be ; fe ; t ; s
)

NB. r=. tgetrs A;X
NB. where r is boxed table with columns: name error time space
tgetrs=: 3 : 0
  'A X'=. y
  A=. A + diagmat (# A) $ 1
  'BAX  BXA'=.  X (mp~ ; mp) A
  'BU1X BXU1'=. X (mp~ ; mp) U1=. utri1 A
  'BUX  BXU'=.  X (mp~ ; mp) U=.  utri  A
  'BL1X BXL1'=. X (mp~ ; mp) L1=. ltri1 A
  'BLX  BXL'=.  X (mp~ ; mp) L=.  ltri  A

  r=.  ,: 'getrsu1x'               norm1 gecon ttrs (U1;X;BU1X)
  r=. r , 'getrsux'                norm1 gecon ttrs (U;X;BUX)
  r=. r , 'getrsxu1'               norm1 gecon ttrs (U1;X;BXU1)
  r=. r , 'getrsxu'                norm1 gecon ttrs (U;X;BXU)
  r=. r , 'getrsl1x'               norm1 gecon ttrs (L1;X;BL1X)
  r=. r , 'getrslx'                norm1 gecon ttrs (L;X;BLX)
  r=. r , 'getrsxl1'               norm1 gecon ttrs (L1;X;BXL1)
  r=. r , 'getrsxl'                norm1 gecon ttrs (L;X;BXL)
  r=. r , '%.'                     norm1 gecon ttrs (U;X;BUX)
  r=. r , '(mp~ (128!:1))'         norm1 gecon ttrs (U;X;BUX)
  r=. r , 'gesv'                   norm1 gecon ttrs (A;X;BAX)
  r=. r , '(gesv_jlapack_ @ (;~))' norm1 gecon ttrs (A;X;BAX)
)

NB. r=. tpotrs A
NB. where r is boxed table with columns: name error time space
tpotrs=: 3 : 0
  0 5 $ a:
)

NB. ---------------------------------------------------------
NB. testtrs
NB. Adverb to test triangular solver algorithms
NB.
NB. Syntax:
NB.   r=. mkge testtrs m,n
NB. where
NB.   m,n  - 2-vector of integers, shape of random matrices
NB.          to test algorithms; if m≠n then algorithms that
NB.          accept square matrices only are skipped
NB.   mkge - monadic verb to generate random non-singular
NB.          general y-matrix (shape is taken from y)
NB.   r    - boxed table with 4 columns: 'algorithm name'
NB.          'errors' 'time, sec.' 'space, bytes'
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   r=. (_1 1 0 26 _16 8 & gemat) testtrs 7 7

testtrs=: 1 : 0
  r=.     (tpotrs @ ((u pomat) ; u) @ {. ^: (=/)) y     NB. FIXME!
  r=. r , (tgetrs @ (u ; (u @ {.)))               y
)
