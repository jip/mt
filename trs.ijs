NB. trs.ijs
NB. Solve linear monomial equation via triangular
NB. factorization
NB.
NB. trtrsux    Solve equation U * X = B, where U is an upper
NB.            triangular matrix
NB. trtrsutx   Solve equation U^T * X = B, where U is an
NB.            upper triangular matrix
NB. trtrsuhx   Solve equation U^H * X = B, where U is an
NB.            upper triangular matrix
NB. trtrsu1x   Solve equation U1 * X = B, where U1 is an unit
NB.            upper triangular matrix
NB. trtrsu1tx  Solve equation U1^T * X = B, where U1 is an
NB.            unit upper triangular matrix
NB. trtrsu1hx  Solve equation U1^H * X = B, where U1 is an
NB.            unit upper triangular matrix
NB. trtrsxu    Solve equation X * U = B, where U is an upper
NB.            triangular matrix
NB. trtrsxut   Solve equation X * U^T = B, where U is an
NB.            upper triangular matrix
NB. trtrsxuh   Solve equation X * U^H = B, where U is an
NB.            upper triangular matrix
NB. trtrsxu1   Solve equation X * U1 = B, where U1 is an unit
NB.            upper triangular matrix
NB. trtrsxu1t  Solve equation X * U1^T = B, where U1 is an
NB.            unit upper triangular matrix
NB. trtrsxu1h  Solve equation X * U1^H = B, where U1 is an
NB.            unit upper triangular matrix
NB. trtrslx    Solve equation L * X = B, where L is a lower
NB.            triangular matrix
NB. trtrsltx   Solve equation L^T * X = B, where L is a lower
NB.            triangular matrix
NB. trtrslhx   Solve equation L^H * X = B, where L is a lower
NB.            triangular matrix
NB. trtrsl1x   Solve equation L1 * X = B, where L1 is an unit
NB.            lower triangular matrix
NB. trtrsl1tx  Solve equation L1^T * X = B, where L1 is an
NB.            unit lower triangular matrix
NB. trtrsl1hx  Solve equation L1^H * X = B, where L1 is an
NB.            unit lower triangular matrix
NB. trtrsxl    Solve equation X * L = B, where L is a lower
NB.            triangular matrix
NB. trtrsxlt   Solve equation X * L^T = B, where L is a lower
NB.            triangular matrix
NB. trtrsxlh   Solve equation X * L^H = B, where L is a lower
NB.            triangular matrix
NB. trtrsxl1   Solve equation X * L1 = B, where L1 is an unit
NB.            lower triangular matrix
NB. trtrsxl1t  Solve equation X * L1^T = B, where L1 is an
NB.            unit lower triangular matrix
NB. trtrsxl1h  Solve equation X * L1^H = B, where L1 is an
NB.            unit lower triangular matrix
NB.
NB. getrsax    Solve equation A * X = B, where A is a general
NB.            matrix
NB. getrsatx   Solve equation A^T * X = B, where A is a
NB.            general matrix
NB. getrsahx   Solve equation A^H * X = B, where A is a
NB.            general matrix
NB. getrsxa    Solve equation X * A = B, where A is a general
NB.            matrix
NB. getrsxat   Solve equation X * A^T = B, where A is a
NB.            general matrix
NB. getrsxah   Solve equation X * A^H = B, where A is a
NB.            general matrix
NB.
NB. hetrsax    Solve equation A * X = B, where A is a
NB.            Hermitian (symmetric) matrix
NB. hetrsatx   Solve equation A^T * X = B, where A is a
NB.            Hermitian (symmetric) matrix
NB. hetrsxa    Solve equation X * A = B, where A is a
NB.            Hermitian (symmetric) matrix
NB. hetrsxat   Solve equation X * A^T = B, where A is a
NB.            Hermitian (symmetric) matrix
NB.
NB. potrs      Solve equation A*X=B, where A is a Hermitian
NB.            (symmetric) positive definite matrix
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. trtrs
NB. Template adverb to make triangular solver verbs
NB.
NB. Syntax:
NB.   u0=. u0`u1`u2`u3`u4`u5`u6 trtrs
NB. where
NB.   u0 - solver verb to do recursive call
NB.        warning: avoid to change u0's name class!
NB.   u1 - extract either the square matrix Aa from top left
NB.        or bottom right corner, or the corresponding RHS
NB.        part Ba
NB.   u2 - extract either the square matrix Ab opposite to
NB.        Aa on diagonal, or the corresponding RHS part Bb
NB.   u3 - combine intervals to extract either top right part
NB.        (tru case), or bottom left part (trl case) from
NB.        input matrix A
NB.   u4 - [commuted] matrix product, optionally pre-process
NB.        Xa and post-process result
NB.   u5 - assemble solution X from halves Xa and Xb
NB.   u6 - optionally scale output, when diagonal is not unit

trtrs=: 1 : 0
:
  '`u0 u1 u2 u3 u4 u5 u6'=. u
  n=. # x
  if. n > 1 do.
    p=. >. -: n
    Aa=. (2 $ p) u1 x
    Ba=. p u1 y
    Xa=. Aa u0 Ba
    Ab=. (2 $ p) u2 x
    Bb=. p u2 y
    Ac=. (p u3 (p - n)) {. x
    Xa u5 (Ab u0 (Bb - Ac u4 Xa))
  else.
    y u6 x
  end.
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:          Solves:           Syntax:
NB. trtrsux        U    * X = B      X=. A trtrsux B
NB. trtrsutx       U^T  * X = B      X=. A trtrsutx B
NB. trtrsuhx       U^H  * X = B      X=. A trtrsuhx B
NB. trtrsu1x       U1   * X = B      X=. A trtrsu1x B
NB. trtrsu1tx      U1^T * X = B      X=. A trtrsu1tx B
NB. trtrsu1hx      U1^H * X = B      X=. A trtrsu1hx B
NB. trtrslx        L    * X = B      X=. A trtrslx B
NB. trtrsltx       L^T  * X = B      X=. A trtrsltx B
NB. trtrslhx       L^H  * X = B      X=. A trtrslhx B
NB. trtrsl1x       L1   * X = B      X=. A trtrsl1x B
NB. trtrsltx       L1^T * X = B      X=. A trtrsl1tx B
NB. trtrslhx       L1^H * X = B      X=. A trtrsl1hx B
NB.
NB. Description:
NB.   solve triangular system
NB. where:
NB.   A    - n×n-matrix, containing either U, U1, L or L1
NB.   U    - n×n upper triangular matrix
NB.   U1   - n×n unit upper triangular matrix (diagonal is
NB.          not saved)
NB.   L    - n×n lower triangular matrix
NB.   L1   - n×n unit lower triangular matrix (diagonal is
NB.          not saved)
NB.   B    - n-vector or n×nrhs-matrix, the RHS
NB.   X    - same shape as B, the solution
NB.   n    ≥ 0
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - opposite triangle is not referenced
NB. - unit diagonal is not referenced

trtrsux=:   trtrsux  `}.`{.` ,  `     mp      `(,~)`(% 0&({,)) trtrs
trtrsutx=:  trtrsutx `{.`}.` ,  `(|:@(mp~ |:))` ,  `(% 0&({,)) trtrs
trtrsuhx=:  trtrsuhx `{.`}.` ,  `(ct@(mp~ ct))` ,  `(% 0&({,)) trtrs
trtrsu1x=:  trtrsu1x `}.`{.` ,  `     mp      `(,~)`[          trtrs
trtrsu1tx=: trtrsu1tx`{.`}.` ,  `(|:@(mp~ |:))` ,  `[          trtrs
trtrsu1hx=: trtrsu1hx`{.`}.` ,  `(ct@(mp~ ct))` ,  `[          trtrs
trtrslx=:   trtrslx  `{.`}.`(,~)`     mp      ` ,  `(% 0&({,)) trtrs
trtrsltx=:  trtrsltx `}.`{.`(,~)`(|:@(mp~ |:))`(,~)`(% 0&({,)) trtrs
trtrslhx=:  trtrslhx `}.`{.`(,~)`(ct@(mp~ ct))`(,~)`(% 0&({,)) trtrs
trtrsl1x=:  trtrsl1x `{.`}.`(,~)`     mp      ` ,  `[          trtrs
trtrsl1tx=: trtrsl1tx`}.`{.`(,~)`(|:@(mp~ |:))`(,~)`[          trtrs
trtrsl1hx=: trtrsl1hx`}.`{.`(,~)`(ct@(mp~ ct))`(,~)`[          trtrs

NB. ---------------------------------------------------------
NB. Verb:          Solves:           Syntax:
NB. trtrsxu        X * U    = B      X=. A trtrsxu B
NB. trtrsxut       X * U^T  = B      X=. A trtrsxut B
NB. trtrsxuh       X * U^H  = B      X=. A trtrsxuh B
NB. trtrsxu1       X * U1   = B      X=. A trtrsxu1 B
NB. trtrsxu1t      X * U1^T = B      X=. A trtrsxu1t B
NB. trtrsxu1h      X * U1^H = B      X=. A trtrsxu1h B
NB. trtrsxl        X * L    = B      X=. A trtrsxl B
NB. trtrsxlt       X * L^T  = B      X=. A trtrsxlt B
NB. trtrsxlh       X * L^H  = B      X=. A trtrsxlh B
NB. trtrsxl1       X * L1   = B      X=. A trtrsxl1 B
NB. trtrsxl1t      X * L1^T = B      X=. A trtrsxl1t B
NB. trtrsxl1h      X * L1^H = B      X=. A trtrsxl1h B
NB.
NB. Description:
NB.   solve triangular system
NB. where:
NB.   A    - n×n-matrix, containing either U, U1, L or L1
NB.   U    - n×n upper triangular matrix
NB.   U1   - n×n unit upper triangular matrix (diagonal is
NB.          not saved)
NB.   L    - n×n lower triangular matrix
NB.   L1   - n×n unit lower triangular matrix (diagonal is
NB.          not saved)
NB.   B    - n-vector or nrhs×n-matrix, the RHS
NB.   X    - same shape as B, the solution
NB.   n    ≥ 0
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - opposite triangle is not referenced
NB. - unit diagonal is not referenced

trtrsxu=:   |:@trtrsutx |:
trtrsxut=:  |:@trtrsux |:
trtrsxuh=:  ct@trtrsux ct
trtrsxu1=:  |:@trtrsu1tx |:
trtrsxu1t=: |:@trtrsu1x |:
trtrsxu1h=: ct@trtrsu1x ct
trtrsxl=:   |:@trtrsltx |:
trtrsxlt=:  |:@trtrslx |:
trtrsxlt=:  ct@trtrslx ct
trtrsxl1=:  |:@trtrsl1tx |:
trtrsxl1t=: |:@trtrsl1x |:
trtrsxl1t=: ct@trtrsl1x ct

NB. ---------------------------------------------------------
NB. Verb:          Solves:              Syntax:
NB. getrsax        A   * X = B          X=. A getrsax B
NB. getrsatx       A^T * X = B          X=. A getrsatx B
NB. getrsahx       A^H * X = B          X=. A getrsahx B
NB. getrsxa        X * A   = B          X=. A getrsxa B
NB. getrsxat       X * A^T = B          X=. A getrsxat B
NB. getrsxah       X * A^H = B          X=. A getrsxah B
NB.
NB. Description:
NB.   solve system via LU factorization
NB. where:
NB.   A    - n×n-matrix
NB.   B    - n-vector or n×nrhs-matrix, the RHS
NB.   X    - same shape as B, the solution
NB.   n    ≥ 0
NB.   nrhs ≥ 0

getrsax=:  (((1 {:: ]) ([ trtrsux trtrsl1x) (C.~ (0 & {::))) getrfpl1u)~
getrsatx=: (((0 {:: ]) C. (([ trtrsltx trtrsu1tx)~ (1 & {::))) getrflu1)~
getrsahx=: (((0 {:: ]) C. (([ trtrslhx trtrsu1hx)~ (1 & {::))) getrflu1)~
getrsxa=:  (((0 {:: ]) C. (([ trtrsxl trtrsxu1)~ (1 & {::))) getrflu1)~
getrsxat=: (((1 {:: ]) ([ trtrsxut trtrsxl1t) (C.~ (0 & {::))) getrfpl1u)~
getrsxah=: (((1 {:: ]) ([ trtrsxuh trtrsxl1h) (C.~ (0 & {::))) getrfpl1u)~

NB. ---------------------------------------------------------
NB. Verb:          Solves:              Syntax:
NB. hetrsax        A   * X   = B        X=. A hetrsax B
NB. hetrsatx       A^T * X   = B        X=. A hetrsatx B
NB. hetrsxa        X   * A   = B        X=. A hetrsxa B
NB. hetrsxat       X   * A^T = B        X=. A hetrsxat B
NB.
NB. Description:
NB.   solve system via ? factorization
NB. where:
NB.   A    - n×n Hermitian (symmetric) matrix
NB.   B    - n-vector or n×nrhs-matrix, the RHS
NB.   X    - same shape as B, the solution
NB.   n    ≥ 0
NB.   nrhs ≥ 0

hetrsax=:  [:
hetrsatx=: [:
hetrsahx=: [:
hetrsxa=:  [:
hetrsxat=: [:
hetrsxah=: [:

NB. ---------------------------------------------------------
NB. Verb:          Solves:              Syntax:
NB. potrsax        A   * X = B          X=. A potrsax B
NB. potrsatx       A^T * X = B          X=. A potrsatx B
NB. potrsxa        X * A   = B          X=. A potrsxa B
NB. potrsxat       X * A^T = B          X=. A potrsxat B
NB.
NB. Description:
NB.   solve system via Cholesky factorization
NB. where:
NB.   A    - n×n Hermitian (symmetric) positive definite
NB.          matrix
NB.   B    - n-vector or n×nrhs-matrix, the RHS
NB.   X    - same shape as B, the solution
NB.   n    ≥ 0
NB.   nrhs ≥ 0

potrsax=:  (([ trtrslhx trtrslx)~ potrfl)~
potrsatx=: (([ +@trtrslhx trtrslx)~ potrfl)~ +
potrsxa=:  (([ trtrsxl trtrsxlh)~ potrfl)~
potrsxat=: (([ trtrsxl trtrslhx)~ potrfl)~ +          NB. CHECKME!

NB. =========================================================
NB. Test suite

NB. name vmp ttrs A;X;B;rcond

ttrs=: 1 : 0
:
  'A X B rcond'=. y
  't s'=. timespacex 'Xsol=. A ' , x , ' B'
  be=. (((norm1t (B - A u Xsol)) % (norm1 A)) % (norm1t Xsol)) % FP_EPS  NB. backward error
  fe=. (((normit (X - Xsol)) * rcond) % (normit X)) % FP_EPS             NB. forward error
  prn x ; rcond ; be ; fe ; t ; s
)

NB. ttrtrs A;X
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

  'trtrsux'         mp  ttrs (U ;X;UX ;rcondU )  NB. <<<<<< TODO: send A, not U
  '(mp~ (128!:1))~' mp  ttrs (U ;X;UX ;rcondU )
  'trtrsu1x'        mp  ttrs (U1;X;U1X;rcondU1)
  'trtrsxu'         mp~ ttrs (U ;X;XU ;rcondU )
  '(mp (128!:1))~'  mp~ ttrs (U ;X;XU ;rcondU )
  'trtrsxu1'        mp~ ttrs (U1;X;XU1;rcondU1)
  'trtrslx'         mp  ttrs (L ;X;LX ;rcondL )
  'trtrsl1x'        mp  ttrs (L1;X;L1X;rcondL1)
  'trtrsxl'         mp~ ttrs (L ;X;XL ;rcondL )
  'trtrsxl1'        mp~ ttrs (L1;X;XL1;rcondL1)
  EMPTY
)

NB. tgetrs A;X

tgetrs=: 3 : 0
  'A X'=. y
  A=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) A
  'AX XA'=.  X (mp~ ; mp) A
  rcondA=. norm1 con getri A
  '%.~'                 mp ttrs (A;X;AX;rcondA)
  'getrsax'             mp ttrs (A;X;AX;rcondA)
  '(gesv_jlapack_ @ ;)' mp ttrs (A;X;AX;rcondA)
  EMPTY
)

NB. ---------------------------------------------------------
NB. testtrs
NB. Adverb to test triangular solver algorithms
NB.
NB. Syntax:
NB.   mkge testtrs m,n
NB. where
NB.   m,n  - 2-vector of integers, shape of random matrices
NB.          to test algorithms; tests are run only when m=n
NB.   mkge - monadic verb to generate random non-singular
NB.          general y-matrix (shape is taken from y)
NB.
NB. Application:
NB.   NB. with limited random matrix values' amplitudes
NB.   cocurrent 'mt'
NB.   (_1 1 0 16 _6 4 & gemat) testtrs 500 500
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testtrs 500 500
NB.
NB. Notes:
NB. - diagonalizable matrices are processed the same way as
NB.   general matrices

testtrs=: 1 : 0
  ttrtrs @ (u ; (u @ {.))       y  NB. FIXME! add: (^: (=/))
  tgetrs @ (u ; (u @ {.))       y
  thetrs @ ((u hemat) ; u) @ {. y
  tpotrs @ ((u pomat) ; u) @ {. y
)
