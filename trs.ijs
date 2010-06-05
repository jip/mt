NB. trs.ijs
NB. Solve monomial equation via triangular factorization
NB.
NB. trtrsux   Solve equation U*X=B, where U is an upper
NB.           triangular matrix
NB. trtrsu1x  Solve equation U1*X=B, where U1 is an unit
NB.           upper triangular matrix
NB. trtrsxu   Solve equation X*U=B, where U is an upper
NB.           triangular matrix
NB. trtrsxu1  Solve equation X*U1=B, where U1 is an unit
NB.           upper triangular matrix
NB. trtrslx   Solve equation L*X=B, where L is a lower
NB.           triangular matrix
NB. trtrsl1x  Solve equation L1*X=B, where L1 is an unit
NB.           lower triangular matrix
NB. trtrsxl   Solve equation X*L=B, where L is a lower
NB.           triangular matrix
NB. trtrsxl1  Solve equation X*L1=B, where L1 is an unit
NB.           lower triangular matrix
NB. getrs     Solve equation A*X=B, where A is a general
NB.           matrix
NB. hetrs      Solve equation A*X=B, where A is a Hermitian
NB.           (symmetric) matrix
NB. ditrs      Solve equation A*X=B, where A is a
NB.           diagonalizable matrix
NB. potrs      Solve equation A*X=B, where A is a Hermitian
NB.           (symmetric) positive definite matrix
NB.
NB. Version: 1.0.0 2009-06-01
NB. Copyright: Igor Zhuravlov, igor at uic.dvgu.ru
NB. License: GNU GPL

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. trtrs
NB. Template adverb to make triangular solver verbs
NB.
NB. Syntax:
NB.   vtrtrs=. vtrtrs`u1`u2`u3`u4`u5`u6 trtrs
NB. where
NB.   vtrtrs - solver verb for recursive call
NB.   u1     - either }. for systems: U1*X=B, U*X=B, X*L1=B, X*L=B
NB.            or {. for systems: X*U1=B, X*U=B, L1*X=B, L*X=B
NB.   u2     - either , for systems: U1*X=B, U*X=B, X*U1=B, X*U=B
NB.            or ,~ for systems: L1*X=B, L*X=B, X*L1=B, X*L=B
NB.   u3     - either {. for systems: U1*X=B, U*X=B, X*L1=B, X*L=B
NB.            or }. for systems: X*U1=B, X*U=B, L1*X=B, L*X=B
NB.   u4     - either mp for systems: U1*X=B, U*X=B, L1*X=B, L*X=B
NB.            or mp~ for systems: X*U1=B, X*U=B, X*L1=B, X*L=B
NB.   u5     - either , for systems: X*U1=B, X*U=B, L1*X=B, L*X=B
NB.            or ,~ for systems: U1*X=B, U*X=B, X*L1=B, X*L=B
NB.   u6     - either [ for systems: U1*X=B, X*U1=B, L1*X=B, X*L1=B
NB.            or (% 0&({,)) for systems: U*X=B, X*U=B, L*X=B, X*L=B

trtrs=: 1 : 0
:
  '`u0 u1 u2 u3 u4 u5 u6'=. u
  n=. # x
  if. n > 1 do.
    p=. >. -: n
    Aa=. (2 $ p) u1 x              NB. square matrix from top left or bottom right corner
    Ba=. p u1 y                    NB. RHS part corresponding to Aa
    Xa=. Aa u0 Ba                  NB. find corresponding solution recursively
    Ab=. (2 $ p) u3 x              NB. square matrix opposite to Aa on diagonal
    Bb=. p u3 y                    NB. RHS part corresponding to Ab
    Ac=. (p u2 (p - n)) {. x       NB. either top right part (tru case), or bottom left part (trl case)
    Xa u5 (Ab u0 (Bb - Ac u4 Xa))  NB. find corresponding solution recursively and form complete answer
  else.
    y u6 x
  end.
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. trtrsux
NB. trtrsu1x 
NB. trtrsxu
NB. trtrsxu1
NB. trtrslx
NB. trtrsl1x
NB. trtrsxl
NB. trtrsxl1
NB. Various triangular solvers to solve systems:
NB.   T*X=B or X*T=B
NB. where T is either unit triangular (diagonal is not saved)
NB. or triangular, B is RHS and X is solution
NB.
NB. Syntax:
NB.   X=. T solver B
NB. where
NB.   solver - any verb of: trtrsu1x trtrsux trtrsxu1 trtrsxu
NB.            trtrsl1x trtrslx trtrsxl1 trtrsxl
NB.   B      - n-vector or n×nrhs-matrix, RHS
NB.   T      - m×n-matrix, containing triangular n×n-matrix
NB.            of the system
NB.   X      - same shape as B, the solution
NB.   m      ≥ 0
NB.   n      ≥ 0
NB.   nrhs   ≥ 0
NB.
NB. If:
NB.   BUX=. U mp X
NB.   X=. <<<<<<<<<<<<<<
NB. then
NB.   

erase 'trtrsux';'trtrsu1x';'trtrsux';'trtrsxu1';'trtrsxu';'trtrsl1x';'trtrslx';'trtrsxl1';'trtrsxl'

trtrsux=:  trtrsux `}.` ,  `{.` mp  `(,~)`(% 0&({,)) trtrs  NB. U *X=B, where U  is      upper triangular
trtrsu1x=: trtrsu1x`}.` ,  `{.` mp  `(,~)`[          trtrs  NB. U1*X=B, where U1 is unit upper triangular
trtrsxu=:  trtrsxu `{.` ,  `}.`(mp~)` ,  `(% 0&({,)) trtrs  NB. X*U =B, where U  is      upper triangular
trtrsxu1=: trtrsxu1`{.` ,  `}.`(mp~)` ,  `[          trtrs  NB. X*U1=B, where U1 is unit upper triangular
trtrslx=:  trtrslx `{.`(,~)`}.` mp  ` ,  `(% 0&({,)) trtrs  NB. L *X=B, where L  is      lower triangular
trtrsl1x=: trtrsl1x`{.`(,~)`}.` mp  ` ,  `[          trtrs  NB. L1*X=B, where L1 is unit lower triangular
trtrsxl=:  trtrsxl `}.`(,~)`{.`(mp~)`(,~)`(% 0&({,)) trtrs  NB. X*L =B, where L  is      lower triangular
trtrsxl1=: trtrsxl1`}.`(,~)`{.`(mp~)`(,~)`[          trtrs  NB. X*L1=B, where L1 is unit lower triangular

NB. ---------------------------------------------------------
NB. getrs
NB. Solve system A*X=B, where A is a general matrix via LU
NB. factorization
NB.
NB. Syntax:
NB.   X=. A getrs B
NB. where
NB.   B    - n-vector or n×nrhs-matrix, RHS
NB.   A    - n×n-matrix
NB.   X    - same shape as B, the solution
NB.   n    ≥ 0
NB.   nrhs ≥ 0
NB.
NB. If:
NB.   X=. A getrs B
NB. then
NB.   B -: A mp X
NB.
NB. TODO:
NB. - getrsax  A * X = B (done)
NB. - getrsxa  X * A = B
NB. - getrsatx A^T * X = B
NB. - getrsxat X * A^T = B
NB. - getrsahx A^H * X = B
NB. - getrsxah X * A^H = B

getrs=: (((1 {:: [) ([ trtrsux trtrsl1x) ((0 {:: [) C. ]))~ getrfl1u)~

NB. ---------------------------------------------------------
NB. hetrs
NB. Solve system A*X=B, where A is a Hermitian (symmetric)
NB. matrix
NB.
NB. Syntax:
NB.   X=. A hetrs B
NB. where
NB.   B    - n-vector or n×nrhs-matrix, RHS
NB.   A    - n×n-matrix, Hermitian (symmetric)
NB.   X    - same shape as B, the solution
NB.   n    ≥ 0
NB.   nrhs ≥ 0
NB.
NB. TODO:
NB. - hetrsax  A * X = B
NB. - hetrsxa  X * A = B
NB. - hetrsatx A^T * X = B
NB. - hetrsxat X * A^T = B
NB. - hetrsahx A^H * X = B
NB. - hetrsxah X * A^H = B

hetrs=: [:

NB. ---------------------------------------------------------
NB. ditrs
NB. Solve system A*X=B, where A is a diagonalizable matrix
NB.
NB. Syntax:
NB.   x=. A ditrs B
NB. where
NB.   B    - n-vector or n×nrhs-matrix, RHS
NB.   A    - n×n-matrix, diagonalizable
NB.   X    - same shape as B, the solution
NB.   n    ≥ 0
NB.   nrhs ≥ 0
NB.
NB. TODO:
NB. - ditrsax  A * X = B
NB. - ditrsxa  X * A = B
NB. - ditrsatx A^T * X = B
NB. - ditrsxat X * A^T = B
NB. - ditrsahx A^H * X = B
NB. - ditrsxah X * A^H = B

ditrs=: [:

NB. ---------------------------------------------------------
NB. potrs
NB. Solve system A*X=B, where A is a Hermitian (symmetric)
NB. positive definite matrix
NB.
NB. Syntax:
NB.   X=. A potrs B
NB. where
NB.   B    - n-vector or n×nrhs-matrix, RHS
NB.   A    - n×n-matrix, a Hermitian (symmetric) positive
NB.          definite
NB.   X    - same shape as B, the solution
NB.   n    ≥ 0
NB.   nrhs ≥ 0
NB.
NB. TODO:
NB. - potrsax  A * X = B
NB. - potrsxa  X * A = B
NB. - potrsatx A^T * X = B
NB. - potrsxat X * A^T = B
NB. - potrsahx A^H * X = B
NB. - potrsxah X * A^H = B

potrs=: [:

NB. =========================================================
NB. Test suite

NB. 'name rcond bwerror fwerror time space'=. name vmp ttrs A;X;B;rcond
ttrs=: 1 : 0
:
  'A X B rcond'=. y
  't s'=. timespacex 'Xsol=. A ' , x , ' B'
  be=. (((norm1t (B - A u Xsol)) % (norm1 A)) % (norm1t Xsol)) % FP_EPS  NB. backward error
  fe=. (((normit (X - Xsol)) * rcond) % (normit X)) % FP_EPS             NB. forward error
  dbprn x ; rcond ; be ; fe ; t ; s
)

NB. r=. ttrtrs A;X
NB. where r is boxed table with columns: name rcond bwerror fwerror time space
ttrtrs=: 3 : 0
  'A X'=. y
  'm n'=. $ A
  se=. 0 j. (m - n)
  mn=. m <. n
  A=. (+ (se diagmat ((mn&$)@:(10&*)@:((*@diag) * (>./@:|@,))))) A

  'UX  XU'=.  X (mp~ ,: mp) U=.  tru  A
  'U1X XU1'=. X (mp~ ,: mp) U1=. tru1 A
  'LX  XL'=.  X (mp~ ,: mp) L=.  trl  A
  'L1X XL1'=. X (mp~ ,: mp) L1=. trl1 A

  rcondU=.  norm1 con trtriu  U
  rcondU1=. norm1 con trtriu1 U1
  rcondL=.  norm1 con trtril  L
  rcondL1=. norm1 con trtril1 L1

  r=.  ,: 'trtrsux'         mp  ttrs (U ;X;UX ;rcondU )
  r=. r , '(mp~ (128!:1))~' mp  ttrs (U ;X;UX ;rcondU )
  r=. r , 'trtrsu1x'        mp  ttrs (U1;X;U1X;rcondU1)
  r=. r , 'trtrsxu'         mp~ ttrs (U ;X;XU ;rcondU )
  r=. r , '(mp (128!:1))~'  mp~ ttrs (U ;X;XU ;rcondU )
  r=. r , 'trtrsxu1'        mp~ ttrs (U1;X;XU1;rcondU1)
  r=. r , 'trtrslx'         mp  ttrs (L ;X;LX ;rcondL )
  r=. r , 'trtrsl1x'        mp  ttrs (L1;X;L1X;rcondL1)
  r=. r , 'trtrsxl'         mp~ ttrs (L ;X;XL ;rcondL )
  r=. r , 'trtrsxl1'        mp~ ttrs (L1;X;XL1;rcondL1)
)

NB. r=. tgetrs A;X
NB. where r is boxed table with columns: name rcond bwerror fwerror time space
tgetrs=: 3 : 0
  'A X'=. y
  A=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) A
  'AX XA'=.  X (mp~ ; mp) A
  rcondA=. norm1 con getri A

  r=.  ,: '%.~'                 mp ttrs (A;X;AX;rcondA)
  r=. r , 'getrs'               mp ttrs (A;X;AX;rcondA)
  r=. r , '(gesv_jlapack_ @ ;)' mp ttrs (A;X;AX;rcondA)
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
NB.   r=. (_1 1 0 16 _6 4 & gemat) testtrs 500 500
NB.   r=. (_1 1 0 16 _6 4 & (gemat j. gemat)) testtrs 500 500

testtrs=: 1 : 0
  r=.     ttrtrs @ (u ; (u @ {.))       y
  r=. r , tgetrs @ (u ; (u @ {.))       y
NB.  r=. r , thetrs @ ((u hemat) ; u) @ {. y
NB.  r=. r , tditrs @ ((u dimat) ; u) @ {. y
NB.  r=. r , tpotrs @ ((u pomat) ; u) @ {. y
)
