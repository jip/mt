NB. Solve linear monomial equation
NB.
NB. Interface:
NB.   gesvxxx   Solve equation (op(A) * X = B) or
NB.             (X * op(A) = B), where A is a general matrix;
NB.             op(A) is either A itself, or A^T (the
NB.             transposition of A), or A^H (the conjugate
NB.             transposition of A); B is known right-hand
NB.             side (RHS), X is unknown solution
NB.   hesvxxx   Solve equation (op(A) * X = B) or
NB.             (X * op(A) = B), where A is a Hermitian
NB.             (symmetric) matrix; op(A) is either A itself,
NB.             or A^T (the transposition of A); B is known
NB.             right-hand side (RHS), X is unknown solution
NB.   posvxxx   Solve equation (op(A) * X = B) or
NB.             (X * op(A) = B), where A is a Hermitian
NB.             (symmetric) positive definite matrix; op(A)
NB.             is either A itself, or A^T (the transposition
NB.             of A); B is known right-hand side (RHS), X is
NB.             unknown solution
NB.   ptsvx     Solve equation (op(A) * X = B) or
NB.             (X * op(A) = B), where A is a Hermitian
NB.             (symmetric) positive definite tridiagonal
NB.             matrix; op(A) is either A itself, or A^T (the
NB.             transposition of A); B is known right-hand
NB.             side (RHS), X is unknown solution
NB.
NB. Test suite:
NB.   testgesv  Test linear monomial equation solving
NB.             algorithms by general matrix given
NB.   testhesv  Test linear monomial equation solving
NB.             algorithms by Hermitian (symmetric) matrix
NB.             given
NB.   testposv  Test linear monomial equation solving
NB.             algorithms by Hermitian (symmetric) positive
NB.             definite matrix given
NB.   testsv    Adv. to make verb to test triangular solver
NB.             algorithms by matrix of generator and shape
NB.             given
NB.
NB. Requisites:
NB.   Copyright (C) 2010 Igor Zhuravlov
NB.   For license terms, see the file COPYING in this distribution
NB.   Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:          Solves:         Syntax:
NB. gesvax         A   * X = B     Xv=. A gesvax  Bv
NB. gesvahx        A^H * X = B     Xv=. A gesvahx Bv
NB. gesvatx        A^T * X = B     Xv=. A gesvatx Bv
NB. gesvxa         X * A   = B     Xh=. A gesvxa  Bh
NB. gesvxah        X * A^H = B     Xh=. A gesvxah Bh
NB. gesvxat        X * A^T = B     Xh=. A gesvxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general matrix A
NB.   by triangular factorization:
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
NB. - implements LAPACK's xGESV
NB. - based on (P*L1*U) variant of factorization as the
NB.   fastest among getrfxxxx
NB.
NB. TODO:

gesvax=:  (getrsax ~ getrfpl1u)~
gesvahx=: (getrsahx~ getrfpl1u)~
gesvatx=: (getrsatx~ getrfpl1u)~
gesvxa=:  (getrsxa ~ getrfpl1u)~
gesvxah=: (getrsxah~ getrfpl1u)~
gesvxat=: (getrsxat~ getrfpl1u)~

NB. ---------------------------------------------------------
NB. Verb:          Solves:         Syntax:
NB. hesvax         A   * X = B     Xv=. A hesvax  Bv
NB. hesvatx        A^T * X = B     Xv=. A hesvatx Bv
NB. hesvxa         X * A   = B     Xh=. A hesvxa  Bh
NB. hesvxat        X * A^T = B     Xh=. A hesvxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) matrix, represented in factored form:
NB.     P * L1 * T * L1' * P' = A
NB. where:
NB.   A    - n×n-matrix, Hermitian (symmetric)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - implements LAPACK's xHESV
NB. - based on (P*L1*T*L1'*P') variant of factorization as
NB.   the fastest among hetrfxxxx

hesvax=:  (hetrsax ~ hetrfpl)~
hesvatx=: (hetrsatx~ hetrfpl)~
hesvxa=:  (hetrsxa ~ hetrfpl)~
hesvxat=: (hetrsxat~ hetrfpl)~

NB. ---------------------------------------------------------
NB. Verb:          Solves:         Syntax:
NB. posvax         A   * X = B     Xv=. A posvax  Bv
NB. posvatx        A^T * X = B     Xv=. A posvatx Bv
NB. posvxa         X * A   = B     Xh=. A posvxa  Bh
NB. posvxat        X * A^T = B     Xh=. A posvxat Bh
NB.
NB. Description:
NB.   Solve Hermitian (symmetric) positive definite system
NB.   by Cholesky factorization:
NB.     L * L' = A
NB. where:
NB.   A    - n×n Hermitian (symmetric) positive definite
NB.          matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - implements LAPACK's xPOSV
NB. - based on (L*L') variant of factorization as the
NB.   fastest among potrfx

posvax=:  (potrsax ~ potrfl)~
posvatx=: (potrsatx~ potrfl)~
posvxa=:  (potrsxa ~ potrfl)~
posvxat=: (potrsxat~ potrfl)~

NB. ---------------------------------------------------------
NB. Verb:          Solves:         Syntax:
NB. ptsvax         A   * X = B     Xv=. A ptsvax  Bv
NB. ptsvatx        A^T * X = B     Xv=. A ptsvatx Bv
NB. ptsvxa         X * A   = B     Xh=. A ptsvxa  Bh
NB. ptsvxat        X * A^T = B     Xh=. A ptsvxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) positive definite tridiagonal matrix A via
NB.   factorization:
NB.     L1 * D * L^H = A
NB. where:
NB.   A    - n×n-matrix, Hermitian (symmetric) positive
NB.          definite tridiagonal
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   L1   - n×n-matrix, unit lower bidiangonal
NB.   D    - n×n-matrix, diagonal with positive diagonal
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - implements LAPACK's xPTSV
NB. - based on (L*D*L^H) variant of factorization as the
NB.   fastest among pttrfx

ptsvax=:  pttrfl @ [ pttrsax ]  NB. ##################
ptsvatx=: (potrsatx~ potrfl)~
ptsvxa=:  (potrsxa ~ potrfl)~
ptsvxat=: (potrsxat~ potrfl)~

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgesv
NB.
NB. Description:
NB.   Test linear monomial equation solving algorithms:
NB.   - %. (built-in)
NB.   - gesv (math/lapack addon)
NB.   - gesvxxx (math/mt addon)
NB.   by general matrix given
NB.
NB. Syntax:
NB.   testgesv (A;X)
NB. where
NB.   A - n×n-matrix
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (||op(A)|| * ||X||*eps))

testgesv=: 3 : 0
  'A X'=. y
  'conA conAh conAt'=. (norm1 con getri)"2 (] , ct ,: |:) A

  ('%.' tdyad ((mp & >/)`(0 & {::)`]`(conA "_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp & >/)@[) - (mp~ (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (A;X)

  ('gesv_jlapack_' tmonad (((0 & {::);(mp & >/))`]`(conA"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp & >/)@[) - (mp~ (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (A;X)

  ('gesvax'        tdyad  ((0 & {::)`(mp  & >/)`]`(conA "_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@[) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X)
  ('gesvahx'       tdyad  ((0 & {::)`(mp  & >/)`]`(conAh"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@[) - ((mp~ ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((ct A);X)
  ('gesvatx'       tdyad  ((0 & {::)`(mp  & >/)`]`(conAt"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@[) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((|: A);X)
  ('gesvxa'        tdyad  ((0 & {::)`(mp~ & >/)`]`(conA "_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@[) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (    A ;X)
  ('gesvxah'       tdyad  ((0 & {::)`(mp~ & >/)`]`(conAh"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((ct A);X)
  ('gesvxat'       tdyad  ((0 & {::)`(mp~ & >/)`]`(conAt"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((|: A);X)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhesv
NB.
NB. Description:
NB.   Test linear monomial equation solving algorithms
NB.   hesvxxx by Hermitian (symmetric) matrix given
NB.
NB. Syntax:
NB.   testhesv (A;X)
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric)
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (||op(A)|| * ||X||*eps))

testhesv=: 3 : 0
  'A X'=. y
  'conA conAt'=. (norm1 con hetri)"2 (] ,: |:) A

  ('hesvax'  tdyad ((0 & {::)`(mp  & >/)`]`(conA "_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@[) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X)
  ('hesvatx' tdyad ((0 & {::)`(mp  & >/)`]`(conAt"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@[) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((|: A);X)
  ('hesvxa'  tdyad ((0 & {::)`(mp~ & >/)`]`(conA "_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@[) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (    A ;X)
  ('hesvxat' tdyad ((0 & {::)`(mp~ & >/)`]`(conAt"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((|: A);X)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testposv
NB.
NB. Description:
NB.   Test linear monomial equation solving algorithms
NB.   posvxxx by Hermitian (symmetric) positive definite
NB.   matrix given
NB.
NB. Syntax:
NB.   testposv (A;X)
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (||op(A)|| * ||X||*eps))

testposv=: 3 : 0
  'A X'=. y
  'conA conAt'=. (norm1 con potri)"2 (] ,: |:) A

  ('posvax'  tdyad ((0 & {::)`(mp  & >/)`]`(conA "_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@[) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X)
  ('posvatx' tdyad ((0 & {::)`(mp  & >/)`]`(conAt"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@[) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((|: A);X)
  ('posvxa'  tdyad ((0 & {::)`(mp~ & >/)`]`(conA "_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@[) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (    A ;X)
  ('posvxat' tdyad ((0 & {::)`(mp~ & >/)`]`(conAt"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((|: A);X)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testsv
NB.
NB. Description:
NB.   Adv. to make verb to test triangular solver algorithms
NB.   by matrix of generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testsv
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random square real matrix with limited values'
NB.   amplitudes:
NB.     (_1 1 0 16 _6 4 & gemat_mt_) testsv_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testsv_mt_ 150 200

testsv=: 1 : 'EMPTY_mt_ [ ((testptsv_mt_ @ ((u ptmat_mt_) ; u)) [ (testposv_mt_ @ ((u pomat_mt_) ; u)) [ ((testhesv_mt_ @ ((u hemat_mt_) ; u))) [ (testgesv_mt_ @ (u ; u))) ^: (=/)'
