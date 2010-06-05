NB. trs.ijs
NB. Solve matrix linear monomial equation with general matrix
NB. via triangular factorization
NB.
NB. getrsax   Solve equation A * X = B, where A is a general
NB.           matrix
NB. getrsahx  Solve equation A^H * X = B, where A is a
NB.           general matrix
NB. getrsatx  Solve equation A^T * X = B, where A is a
NB.           general matrix
NB. getrsxa   Solve equation X * A = B, where A is a general
NB.           matrix
NB. getrsxah  Solve equation X * A^H = B, where A is a
NB.           general matrix
NB. getrsxat  Solve equation X * A^T = B, where A is a
NB.           general matrix
NB.
NB. hetrsax   Solve equation A * X = B, where A is a
NB.           Hermitian (symmetric) matrix
NB. hetrsatx  Solve equation A^T * X = B, where A is a
NB.           Hermitian (symmetric) matrix
NB. hetrsxa   Solve equation X * A = B, where A is a
NB.           Hermitian (symmetric) matrix
NB. hetrsxat  Solve equation X * A^T = B, where A is a
NB.           Hermitian (symmetric) matrix
NB.
NB. potrsax   Solve equation A * X = B, where A is a
NB.           Hermitian (symmetric) positive definite matrix
NB. potrsatx  Solve equation A^T * X = B, where A is a
NB.           Hermitian (symmetric) positive definite matrix
NB. potrsxa   Solve equation X * A = B, where A is a
NB.           Hermitian (symmetric) positive definite matrix
NB. potrsxat  Solve equation X * A^T = B, where A is a
NB.           Hermitian (symmetric) positive definite matrix
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:          Solves:           Syntax:
NB. getrsax        A   * X = B       Xv=. A getrsax  Bv
NB. getrsahx       A^H * X = B       Xv=. A getrsahx Bv
NB. getrsatx       A^T * X = B       Xv=. A getrsatx Bv
NB. getrsxa        X * A   = B       Xh=. A getrsxa  Bh
NB. getrsxah       X * A^H = B       Xh=. A getrsxah Bh
NB. getrsxat       X * A^T = B       Xh=. A getrsxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general matrix via
NB.   triangular factorization:
NB.     P * L1 * U = A
NB. where:
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - based on (P*L1*U) variant of factorization as the
NB.   fastest among getrfxxxx
NB.
NB. TODO:
NB. - emulate LAPACK's xGETRS

getrsax=:  (((1 {:: ]) ([ trsmux trsml1x) (C.~ (0 & {::))) getrfpl1u)~
getrsahx=: (((0 {:: ]) C. (([ trsmlhx trsmu1hx)~ (1 & {::))) getrflu1)~
getrsatx=: (((0 {:: ]) C. (([ trsmltx trsmu1tx)~ (1 & {::))) getrflu1)~
getrsxa=:  (((0 {:: ]) C. (([ trsmxl trsmxu1)~ (1 & {::))) getrflu1)~
getrsxah=: (((1 {:: ]) ([ trsmxuh trsmxl1h) (C.~ (0 & {::))) getrfpl1u)~
getrsxat=: (((1 {:: ]) ([ trsmxut trsmxl1t) (C.~ (0 & {::))) getrfpl1u)~

NB. ---------------------------------------------------------
NB. Verb:          Solves:           Syntax:
NB. hetrsax        A   * X   = B     X=. A hetrsax  B
NB. hetrsatx       A^T * X   = B     X=. A hetrsatx B
NB. hetrsxa        X   * A   = B     X=. A hetrsxa  B
NB. hetrsxat       X   * A^T = B     X=. A hetrsxat B
NB.
NB. Description:
NB.   Solve Hermitian (symmetric) system via factorization:
NB.     P * L1 * T * L1' * P' = A
NB. where:
NB.   A    - n×n-matrix, Hermitian (symmetric)
NB.   B    - n-vector or n×nrhs-matrix, the RHS
NB.   X    - same shape as B, the solution
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - based on (P*L1*T*L1'*P') variant of factorization as
NB.   the fastest among hetrfxxxx
NB.
NB. TODO:
NB. - emulate LAPACK's xHETRS

hetrsax=:  (((0 {:: ]) C.^:_1 ((1 {:: ]) trsml1hx ((2 {:: ]) httrstx  ((1 {:: ]) trsml1x  (C.    ~ (0 & {::)))))) hetrfpl)~
hetrsatx=: (((0 {:: ]) C.^:_1 ((1 {:: ]) trsml1tx ((2 {:: ]) httrsttx ((1 {:: ]) trsml1cx (C.    ~ (0 & {::)))))) hetrfpl)~
hetrsxa=:  (((0 {:: ]) C.     ((1 {:: ]) trsmxl1  ((2 {:: ]) httrsxt  ((1 {:: ]) trsmxl1h (C.^:_1~ (0 & {::)))))) hetrfpl)~
hetrsxat=: (((0 {:: ]) C.     ((1 {:: ]) trsmxl1c ((2 {:: ]) httrsxtt ((1 {:: ]) trsmxl1t (C.^:_1~ (0 & {::)))))) hetrfpl)~

NB. ---------------------------------------------------------
NB. Verb:          Solves:           Syntax:
NB. potrsax        A   * X = B       X=. A potrsax  B
NB. potrsatx       A^T * X = B       X=. A potrsatx B
NB. potrsxa        X * A   = B       X=. A potrsxa  B
NB. potrsxat       X * A^T = B       X=. A potrsxat B
NB.
NB. Description:
NB.   Solve Hermitian (symmetric) positive definite system
NB.   via Cholesky factorization:
NB.     L * L' = A
NB. where:
NB.   A    - n×n Hermitian (symmetric) positive definite
NB.          matrix
NB.   B    - n-vector or n×nrhs-matrix, the RHS
NB.   X    - same shape as B, the solution
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - based on (L*L) variant of factorization as the
NB.   fastest among potrfxxxx
NB.
NB. TODO:
NB. - emulate LAPACK's xPOTRS

potrsax=:  (([   trsmlhx trsmlx)~ potrfl)~
potrsatx=: (([ +@trsmlhx trsmlx)~ potrfl)~ +
potrsxa=:  (([ trsmxl trsmxlh)~ potrfl)~
potrsxat=: (([ trsmxl trsmlhx)~ potrfl)~ +          NB. CHECKME!

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgetrs
NB.
NB. Description:
NB.   Test linear monomial equation solving algorithms:
NB.   - %. (built-in)
NB.   - getrsmxxxx (math/mt addon)
NB.   by general matrix given
NB.
NB. Syntax:
NB.   testgetrs (A;X)
NB. where
NB.   A - n×n-matrix
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (||op(A)|| * ||X||*eps))

testgetrs=: 3 : 0
  'A X'=. y
  conA=. (norm1 con getri) A

  ('%.'       tdyad (( mp       & >/)`(0&{::)`]`(conA"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((( mp       & >/)@[) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X)
  ('getrsax'  tdyad ((0&{::)`( mp       & >/)`]`(conA"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((( mp       & >/)@[) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X)
  ('getrsahx' tdyad ((0&{::)`((mp~ ct)~ & >/)`]`(conA"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ ct)~ & >/)@[) - ((mp~ ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((ct A);X)
  ('getrsatx' tdyad ((0&{::)`((mp~ |:)~ & >/)`]`(conA"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ |:)~ & >/)@[) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((|: A);X)
  ('getrsxa'  tdyad ((0&{::)`( mp~      & >/)`]`(conA"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((( mp     ~ & >/)@[) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (A;X)
  ('getrsxah' tdyad ((0&{::)`((mp  ct)~ & >/)`]`(conA"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  ct)~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((ct A);X)
  ('getrsxat' tdyad ((0&{::)`((mp  |:)~ & >/)`]`(conA"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  |:)~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((|: A);X)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhetrs
NB.
NB. Description:
NB.   Test linear monomial equation solving algorithms
NB.   hetrsxxx by Hermitian (symmetric) matrix given
NB.
NB. Syntax:
NB.   testhetrs (A;X)
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric)
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:########################################
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (||op(A)|| * ||X||*eps))

testgetrs=: 3 : 0
  'A X'=. y
  conA=. (norm1 con getri) A

  ('%.'       tdyad (( mp       & >/)`(0&{::)`]`(conA"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((( mp       & >/)@[) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X)
  ('getrsax'  tdyad ((0&{::)`( mp       & >/)`]`(conA"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((( mp       & >/)@[) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X)
  ('getrsahx' tdyad ((0&{::)`((mp~ ct)~ & >/)`]`(conA"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ ct)~ & >/)@[) - ((mp~ ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((ct A);X)
  ('getrsatx' tdyad ((0&{::)`((mp~ |:)~ & >/)`]`(conA"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ |:)~ & >/)@[) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((|: A);X)
  ('getrsxa'  tdyad ((0&{::)`( mp~      & >/)`]`(conA"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((( mp     ~ & >/)@[) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (A;X)
  ('getrsxah' tdyad ((0&{::)`((mp  ct)~ & >/)`]`(conA"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  ct)~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((ct A);X)
  ('getrsxat' tdyad ((0&{::)`((mp  |:)~ & >/)`]`(conA"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  |:)~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((|: A);X)

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

testtrs=: 1 : 'EMPTY_mt_ [ ((testpotrs_mt_ @ ((u pomat_mt_) ; u)) [ ((testhetrs_mt_ @ ((u hemat_mt_) ; u))) [ (testgetrs_mt_ @ (u ; u))) ^: (=/)'
