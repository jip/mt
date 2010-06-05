NB. trs.ijs
NB. Solve linear monomial equation from triangular
NB. factorization
NB.
NB. getrsxxx  Solve equation (op(A) * X = B) or
NB.           (X * op(A) = B), where A is a general matrix,
NB.           represented in factored form; op(A) is either
NB.           A itself, or A^T (the transposition of A), or
NB.           A^H (the conjugate transposition of A); B is
NB.           known right-hand side (RHS), X is unknown
NB.           solution
NB. hetrsxxx  Solve equation (op(A) * X = B) or
NB.           (X * op(A) = B), where A is a Hermitian
NB.           (symmetric) matrix, represented in factored
NB.           form; op(A) is either A itself, or A^T (the
NB.           transposition of A); B is known right-hand
NB.           side (RHS), X is unknown solution
NB. potrsxxx  Solve equation (op(A) * X = B) or
NB.           (X * op(A) = B), where A is a Hermitian
NB.           (symmetric) positive definite matrix,
NB.           represented as Cholesky triangle; op(A) is
NB.           either A itself, or A^T (the transposition of
NB.           A); B is known right-hand side (RHS), X is
NB.           unknown solution
NB. pttrsxxx  Solve equation (op(A) * X = B) or
NB.           (X * op(A) = B), where A is a Hermitian
NB.           (symmetric) positive definite tridiagonal
NB.           matrix, represented as superdiagonal linked to
NB.           diagonal; op(A) is either A itself, or A^T (the
NB.           transposition of A); B is known right-hand side
NB.           (RHS), X is unknown solution
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
NB. Verb:          Solves:         Syntax:
NB. getrsax        A   * X = B     Xv=. (ip;L1U) getrsax  Bv
NB. getrsahx       A^H * X = B     Xv=. (ip;L1U) getrsahx Bv
NB. getrsatx       A^T * X = B     Xv=. (ip;L1U) getrsatx Bv
NB. getrsxa        X * A   = B     Xh=. (ip;L1U) getrsxa  Bh
NB. getrsxah       X * A^H = B     Xh=. (ip;L1U) getrsxah Bh
NB. getrsxat       X * A^T = B     Xh=. (ip;L1U) getrsxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general matrix A,
NB.   represented in factored form:
NB.     P * L1 * U = A
NB. where:
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   ip   - n-vector, rows inversed permutation of A
NB.   L1U  - n×n-matrix, upper triangle contains U, and
NB.          strict lower triangle contains L1 without unit
NB.          diagonal
NB.   L1   - n×n-matrix, unit lower triangular
NB.   U    - n×n-matrix, upper triangular
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - models LAPACK's xGETRS
NB. - based on (P*L1*U) variant of factorization as the
NB.   fastest among getrfxxxx
NB.
NB. TODO:
NB. - implement LAPACK's xGETRS and choose the best

getrsax=:  (1 {:: [) ([ trsmux  trsml1x ) ((0 {:: [) C.       ])
getrsxah=: (1 {:: [) ([ trsmxuh trsmxl1h) ((0 {:: [) C.^:_1"1 ])
getrsxah=: (1 {:: [) ([ trsmxut trsmxl1t) ((0 {:: [) C.^:_1"1 ])

getrsahx=: (0 {:: [) C. ((] trsml1hx trsmuhx~) (1 & {::))~
getrsatx=: (0 {:: [) C. ((] trsml1tx trsmutx~) (1 & {::))~
getrsxa=:  (0 {:: [) C. ((] trsmxl1  trsmxu ~) (1 & {::))~

NB. ---------------------------------------------------------
NB. Verb:          Solves:         Syntax:
NB. hetrsax        A   * X = B     Xv=. (ip;L1;T) hetrsax  Bv
NB. hetrsatx       A^T * X = B     Xv=. (ip;L1;T) hetrsatx Bv
NB. hetrsxa        X * A   = B     Xh=. (ip;L1;T) hetrsxa  Bh
NB. hetrsxat       X * A^T = B     Xh=. (ip;L1;T) hetrsxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) matrix, represented in factored form:
NB.     P * L1 * T * L1^H * P^_1 = A
NB. where:
NB.   A    - n×n-matrix, Hermitian (symmetric)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   ip   - n-vector, full inversed permutation of A
NB.   L1   - n×n-matrix, unit lower triangular
NB.   T    - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - based on (P * L1 * T * L1^H * P^_1 = A) variant of
NB.   factorization as the fastest among hetrfxxxx
NB.
NB. TODO:
NB. - implement LAPACK's xHETRS and choose the best

hetrsax=:   (0 {:: ])    C.^:_1  ((1 {:: ]) trsml1hx ((2 {:: [) pttrsax  ((1 {:: [) trsml1x  ((0 {:: [) C.       ]))))
hetrsatx=: ((0 {:: ]) +@(C.^:_1) ((1 {:: ]) trsml1hx ((2 {:: [) pttrsax  ((1 {:: [) trsml1x  ((0 {:: [) C.       ]))))) +
hetrsxa=:   (0 {:: ])    C."1    ((1 {:: ]) trsmxl1  ((2 {:: ]) pttrsxa  ((1 {:: ]) trsmxl1h ((0 {:: [) C.^:_1"1 ]))))
hetrsxat=: ((0 {:: ])    C."1    ((1 {:: ]) trsmxl1  ((2 {:: ]) pttrsxa  ((1 {:: ]) trsmxl1h ((0 {:: [) C.^:_1"1 ]))))) +

NB. ---------------------------------------------------------
NB. Verb:          Solves:         Syntax:
NB. potrsax        A   * X = B     Xv=. L potrsax  Bv
NB. potrsatx       A^T * X = B     Xv=. L potrsatx Bv
NB. potrsxa        X * A   = B     Xh=. L potrsxa  Bh
NB. potrsxat       X * A^T = B     Xh=. L potrsxat Bh
NB.
NB. Description:
NB.   Solve Hermitian (symmetric) positive definite system
NB.   via Cholesky factorization:
NB.     L * L^H = A
NB. where:
NB.   A    - n×n Hermitian (symmetric) positive definite
NB.          matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   L    - n×n-matrix, lower triangular with positive
NB.          diagonal entries, Cholesky triangle
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - based on (L * L^H = A) variant of factorization as the
NB.   fastest among potrfx
NB.
NB. TODO:
NB. - implement LAPACK's xPOTRS and choose the best

potrsax=:   [   trsmlhx trsmlx
potrsatx=: ([ +@trsmlhx trsmlx ) +
potrsxa=:   [   trsmxl  trsmxlh
potrsxat=: ([ +@trsmxl  trsmxlh) +

pttrfsaxstep1=: 3 : 0
  'ein Bin Bout'=. y
  (}. ein) ; (}. Bin) ; (Bout , (({. Bin) - ({. ein) * ({: Bout)))
)

pttrfsaxstep2=: 3 : 0
  'ein din Bin Bout'=. y
  (}: ein) ; (}: din) ; (}: Bin) ; (((Bin (% & {:) din) - (({: ein) * ({. Bout))) , Bout)
)

NB. X=. (L1 ,: D) pttrfsax B
pttrfsax=: 4 : 0
  'L1 D'=. x
  e=. _1 diag L1
  d=. diag D
  y=. _1 {:: pttrfsaxstep1 ^: (<: # L1) (e ; (}. y) ; (1 {. y))
  y=. _1 {:: pttrfsaxstep2 ^: (<: # L1) (e ; (}: d) ; (}: y) ; (y (% & (_1&{.)) d))
)


NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgetrs
NB.
NB. Description:
NB.   Test linear monomial equation solving algorithms
NB.   getrsxxx by general matrix given
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
  'conA conAh conAt'=. (norm1 con getri)"2 (] , ct ,: |:) A
  Af=. getrfpl1u A

  ('getrsax'  tdyad ((_2&{.)`((mp  & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X;Af)
  ('getrsahx' tdyad ((_2&{.)`((mp  & >/)@(2&{.))`]`(conAh"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ((mp~ ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((ct A);X;Af)
  ('getrsatx' tdyad ((_2&{.)`((mp  & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((|: A);X;Af)
  ('getrsxa'  tdyad ((_2&{.)`((mp~ & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (    A ;X;Af)
  ('getrsxah' tdyad ((_2&{.)`((mp~ & >/)@(2&{.))`]`(conAh"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((ct A);X;Af)
  ('getrsxat' tdyad ((_2&{.)`((mp~ & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((|: A);X;Af)

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
NB. Formula:
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (||op(A)|| * ||X||*eps))

testhetrs=: 3 : 0
  'A X'=. y
  'conA conAt'=. (norm1 con hetri)"2 (] ,: |:) A
  Af=. hetrfpl A

  ('hetrsax'  tdyad ((_3&{.)`((mp  & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X;Af)
  ('hetrsatx' tdyad ((_3&{.)`((mp  & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((|: A);X;Af)
  ('hetrsxa'  tdyad ((_3&{.)`((mp~ & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (    A ;X;Af)
  ('hetrsxat' tdyad ((_3&{.)`((mp~ & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((|: A);X;Af)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpotrs
NB.
NB. Description:
NB.   Test linear monomial equation solving algorithms
NB.   potrsxxx by Hermitian (symmetric) positive definite
NB.   matrix given
NB.
NB. Syntax:
NB.   testpotrs (A;X)
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (||op(A)|| * ||X||*eps))

testpotrs=: 3 : 0
  'A X'=. y
  'conA conAt'=. (norm1 con potri)"2 (] ,: |:) A
  L=. potrfpl A

  ('potrsax'  tdyad ((2 & {::)`((mp  & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X;L)
  ('potrsatx' tdyad ((2 & {::)`((mp  & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((|: A);X;L)
  ('potrsxa'  tdyad ((2 & {::)`((mp~ & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (    A ;X;L)
  ('potrsxat' tdyad ((2 & {::)`((mp~ & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((|: A);X;L)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testtrs
NB.
NB. Description:
NB.   Adv. to make verb to test triangular solver algorithms
NB.   by matrix of generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testtrs
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
NB.     (_1 1 0 16 _6 4 & gemat_mt_) testtrs_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testtrs_mt_ 150 200

testtrs=: 1 : 'EMPTY_mt_ [ ((testpotrs_mt_ @ ((u pomat_mt_) ; u)) [ ((testhetrs_mt_ @ ((u hemat_mt_) ; u))) [ (testgetrs_mt_ @ (u ; u))) ^: (=/)'
