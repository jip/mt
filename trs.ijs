NB. Solve linear monomial equation from triangular
NB. factorization
NB.
NB. getrsxxx   Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a general matrix,
NB.            represented in factored form; op(A) is either
NB.            A itself, or A^T (the transposition of A), or
NB.            A^H (the conjugate transposition of A); B is
NB.            known right-hand side (RHS), X is unknown
NB.            solution
NB. hetrsxxx   Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a Hermitian
NB.            (symmetric) matrix, represented in factored
NB.            form; op(A) is either A itself, or A^T (the
NB.            transposition of A); B is known right-hand
NB.            side (RHS), X is unknown solution
NB. potrsxxx   Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a Hermitian
NB.            (symmetric) positive definite matrix,
NB.            represented in factored form; op(A) is either
NB.            A itself, or A^T (the transposition of A); B
NB.            is known right-hand side (RHS), X is unknown
NB.            solution
NB. pttrsxxx   Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a Hermitian
NB.            (symmetric) positive definite tridiagonal
NB.            matrix, represented in factored form; op(A) is
NB.            either A itself, or A^T (the transposition of
NB.            A); B is known right-hand side (RHS), X is
NB.            unknown solution
NB.
NB. testgetrs  Test getrsxxx by general matrix given
NB. testhetrs  Test hetrsxxx by Hermitian (symmetric) matrix
NB.            given
NB. testpotrs  Test potrsxxx by Hermitian (symmetric)
NB.            positive definite matrix given
NB. testpttrs  Test pttrsxxx by Hermitian (symmetric)
NB.            positive definite tridiagonal matrix given
NB. testtrs    Adv. to make verb to test xxtrsxxx by matrix
NB.            of generator and shape given
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
NB.   P    - n×n-matrix, rows permutation of A
NB.   L1   - n×n-matrix, unit lower triangular
NB.   U    - n×n-matrix, upper triangular
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - models LAPACK's xGETRS
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
NB.   P    - n×n-matrix, full permutation of A
NB.   L1   - n×n-matrix, unit lower triangular
NB.   T    - n×n-matrix, Hermitian (symmetric) 3-diagonal
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - models LAPACK's xHETRS
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
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) positive definite matrix A, represented in
NB.   factored form:
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
NB. - models LAPACK's xPOTRS
NB.
NB. TODO:
NB. - implement LAPACK's xPOTRS and choose the best

potrsax=:   [   trsmlhx trsmlx
potrsatx=: ([ +@trsmlhx trsmlx ) +
potrsxa=:   [   trsmxl  trsmxlh
potrsxat=: ([ +@trsmxl  trsmxlh) +

NB. ---------------------------------------------------------
NB. Verb:          Solves:         Syntax:
NB. pttrsax        A   * X = B     Xv=. (L1;D) pttrsax  Bv
NB. pttrsatx       A^T * X = B     Xv=. (L1;D) pttrsatx Bv
NB. pttrsxa        X * A   = B     Xh=. (L1;D) pttrsxa  Bh
NB. pttrsxat       X * A^T = B     Xh=. (L1;D) pttrsxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) positive definite tridiagonal matrix A,
NB.   represented in factored form:
NB.     L1 * D * L1^H = A
NB. where:
NB.   A    - n×n Hermitian (symmetric) positive definite
NB.          tridiagonal matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   L1   - n×n-matrix, unit lower bidiangonal
NB.   D    - n×n-matrix, diagonal with positive diagonal
NB.          entries
NB.   nrhs ≥ 0
NB.
NB. Algorithm for pttrsax:
NB.   In:  L1 D Bv
NB.   Out: Xv
NB.   0) extract main diagonal d from D and subdiagonal e
NB.      from L1
NB.   1) prepare input for iterations:
NB.        ebin=. (e ,. (}. Bv))
NB.        bout=. {. Bv
NB.   2) start iterations k=1:n-1 by Power (^:)
NB.      on (ebin;bout) :
NB.      2.0) extract input for current k-th iteration:
NB.             (e[k] b[k]) := ebin[k]
NB.      2.1) extract b[k-1] produced during previous
NB.           (k-1)-th iteration:
NB.             b[k-1] := bout[-1]
NB.      2.2) find new b[k]:
NB.             b[k] := b[k] - b[k-1]*e[k-1]
NB.      2.3) recombine (shift splitting edge) for next
NB.           iteration:
NB.             ebin=. }. ebin
NB.             bout := bout , b[k]
NB.      then rewrite b by bout
NB.   3) prepare input for iterations:
NB.        bdein=. ((}: b) ,. (}: d) ,. conj(e))
NB.        bout := b[n-1] / d[n-1]
NB.   4) start iterations k=n-2:0 by Power (^:)
NB.      on (bdein;bout) :
NB.      4.0) extract input for current k-th iteration:
NB.             (b[k] d[k] e[k]) := bdein[k]
NB.      4.1) extract b[k-1] produced during previous
NB.           (k-1)-th iteration:
NB.             b[k+1] := bout[0]
NB.      4.2) find new b[k]:
NB.             b[k] := b[k]/d[k] - b[k+1]*e[k]
NB.      4.3) recombine (shift splitting edge) for next
NB.           iteration:
NB.             bdein=. }: bdein
NB.             bout := b[k] , bout
NB.      then rewrite b by bout
NB.   5) return b
NB.
NB. Assertions:
NB.   A -: clean L1 (mp mp (ct@[)) D
NB. where
NB.   'L1 D'=. pttrfl A
NB.
NB. References:
NB. [1] G. H. Golub and C. F. Van Loan, Matrix Computations,
NB.     Johns Hopkins University Press, Baltimore, Md, USA,
NB.     3rd edition, 1996, p. 157
NB.
NB. Notes:
NB. - 'continued fractions' approach is useless here since
NB.   infix scan is non-consequtive
NB. - implements LAPACK's xPTTRS

pttrsax=: 4 : 0
  'L1 D'=. x
  e=. _1 diag L1
  d=. diag D
  y=. _1 {:: ((}.@[ ; ] , (({.@[) ((}.@[) - ]) ((0 ({,) [) * {:@]))) & >/) ^: (# e) (e ((,. }.) ; 1 {. ]) y)
  y=. _1 {:: (((}:@[) ; (({:@[) (((_2}.[)%(_2{[))-({:@[*])) ({.@])) , ]) & >/) ^: (# e) (((+e) ,.~ y (,. & }:) d) ; (y (% & (_1&{.)) d))
)

pttrsatx=: + @ (pttrsax +)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgetrs
NB.
NB. Description:
NB.   Test getrsxxx by general matrix given
NB.
NB. Syntax:
NB.   testgetrs (A;X)
NB. where
NB.   A - n×n-matrix
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (||ε * op(A)|| * ||X||))

testgetrs=: 3 : 0
  'A X'=. y
  'conA conAh conAt'=. 3 # _. NB. (norm1 con (getriul1p@getrful1p))"2 (] , ct ,: |:) A
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
NB.   Test hetrsxxx by Hermitian (symmetric) matrix given
NB.
NB. Syntax:
NB.   testhetrs (A;X)
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric)
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (ε * ||op(A)|| * ||X||))

testhetrs=: 3 : 0
  'A X'=. y
  'conA conAt'=. (norm1 con (hetripl@hetrfpl))"2 (] ,: |:) A
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
NB.   Test potrsxxx by Hermitian (symmetric) positive
NB.   definite matrix given
NB.
NB. Syntax:
NB.   testpotrs (A;X)
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (ε * ||op(A)|| * ||X||))

testpotrs=: 3 : 0
  'A X'=. y
  'conA conAt'=. 2 # _. NB. (norm1 con (potri@potrf))"2 (] ,: |:) A
  L=. potrfl A

  ('potrsax'  tdyad ((2 & {::)`((mp  & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X;L)
  ('potrsatx' tdyad ((2 & {::)`((mp  & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((|: A);X;L)
  ('potrsxa'  tdyad ((2 & {::)`((mp~ & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (    A ;X;L)
  ('potrsxat' tdyad ((2 & {::)`((mp~ & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((|: A);X;L)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpttrs
NB.
NB. Description:
NB.   Test pttrsxxx by Hermitian (symmetric) positive
NB.   definite tridiagonal matrix given
NB.
NB. Syntax:
NB.   testpttrs (A;X)
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.       tridiagonal
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (ε * ||op(A)|| * ||X||))

testpttrs=: 3 : 0
  'A X'=. y
  'conA conAt'=. 2 # _. NB. (norm1 con (pttri@pttrf))"2 (] ,: |:) A NB. ##################
  'L1 D'=. pttrfpl A

  ('pttrsax'  tdyad ((2 & {::)`((mp  & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X;L)
  ('pttrsatx' tdyad ((2 & {::)`((mp  & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((|: A);X;L)
  ('pttrsxa'  tdyad ((2 & {::)`((mp~ & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (    A ;X;L)
  ('pttrsxat' tdyad ((2 & {::)`((mp~ & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((|: A);X;L)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testtrs
NB.
NB. Description:
NB.   Adv. to make verb to test xxtrsxxx by matrix of
NB.   generator and shape given
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
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     (? @ $ 0:) testtrs_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testtrs_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testtrs_mt_ 150 200

testtrs=: 1 : 'EMPTY_mt_ [ ((testpttrs_mt_ @ ((u ptmat_mt_) ; u)) [ (testpotrs_mt_ @ ((u pomat_mt_) ; u)) [ ((testhetrs_mt_ @ ((u hemat_mt_) ; u))) [ (testgetrs_mt_ @ (u ; u))) ^: (=/)'
