NB. Solve linear monomial equation by triangular
NB. factorization
NB.
NB. getrsxxxxxx  Solve equation (op(A) * X = B) or
NB.              (X * op(A) = B), where A is a general
NB.              matrix, represented in factored form; op(A)
NB.              is either A itself, or A^T (the
NB.              transposition of A), or A^H (the conjugate
NB.              transposition of A); B is known right-hand
NB.              side (RHS), X is unknown solution
NB. hetrsxxxx    Solve equation (op(A) * X = B) or
NB.              (X * op(A) = B), where A is a Hermitian
NB.              (symmetric) matrix, represented in factored
NB.              form; op(A) is either A itself, or A^T (the
NB.              transposition of A); B is known right-hand
NB.              side (RHS), X is unknown solution
NB. potrsxxx     Solve equation (op(A) * X = B) or
NB.              (X * op(A) = B), where A is a Hermitian
NB.              (symmetric) positive definite matrix,
NB.              represented in factored form; op(A) is
NB.              either A itself, or A^T (the transposition
NB.              of A); B is known right-hand side (RHS), X
NB.              is unknown solution
NB. pttrsxxx     Solve equation (op(A) * X = B) or
NB.              (X * op(A) = B), where A is a Hermitian
NB.              (symmetric) positive definite tridiagonal
NB.              matrix, represented in factored form; op(A)
NB.              is either A itself, or A^T (the
NB.              transposition of A); B is known right-hand
NB.              side (RHS), X is unknown solution
NB.
NB. testgetrs    Test getrsxxxxxx by general matrix given
NB. testhetrs    Test hetrsxxxx by Hermitian (symmetric)
NB.              matrix given
NB. testpotrs    Test potrsxxx by Hermitian (symmetric)
NB.              positive definite matrix given
NB. testpttrs    Test pttrsxxx by Hermitian (symmetric)
NB.              positive definite tridiagonal matrix given
NB. testtrs      Adv. to make verb to test xxtrsxxx by matrix
NB.              of generator and shape given
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
NB. Verb:         Solves:       Syntax:
NB. getrslu1px    A   * X = B   Xv=. (ip;LU1) getrslu1px  Bv
NB. getrslu1phx   A^H * X = B   Xv=. (ip;LU1) getrslu1phx Bv
NB. getrslu1ptx   A^T * X = B   Xv=. (ip;LU1) getrslu1ptx Bv
NB. getrsxlu1p    X * A   = B   Xh=. (ip;LU1) getrsxlu1p  Bh
NB. getrsxlu1ph   X * A^H = B   Xh=. (ip;LU1) getrsxlu1ph Bh
NB. getrsxlu1pt   X * A^T = B   Xh=. (ip;LU1) getrsxlu1pt Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general square
NB.   matrix A, represented in factored form:
NB.     L * U1 * P = A
NB. where:
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   ip   - n-vector, columns inversed permutation of A
NB.   LU1  - n×n-matrix, lower triangle contains L, and
NB.          strict upper triangle contains U1 without unit
NB.          diagonal
NB.   P    - n×n-matrix, columns permutation of A
NB.   L    - n×n-matrix, lower triangular
NB.   U1   - n×n-matrix, unit upper triangular
NB.   nrhs ≥ 0, number of RHSs

getrslu1px=:  (0 {:: [) C.       ((] trsmu1x  trsmlx ~) (1 & {::))~
getrsxlu1ph=: (0 {:: [) C.^:_1"1 ((] trsmxu1h trsmxlh~) (1 & {::))~
getrsxlu1pt=: (0 {:: [) C.^:_1"1 ((] trsmxu1t trsmxlt~) (1 & {::))~

getrslu1phx=: (1 {:: [) ([ trsmlhx trsmu1hx) ((0 {:: [) C.^:_1   ])
getrslu1ptx=: (1 {:: [) ([ trsmltx trsmu1tx) ((0 {:: [) C.^:_1   ])
getrsxlu1p=:  (1 {:: [) ([ trsmxl  trsmxu1 ) ((0 {:: [) C.       ])

NB. ---------------------------------------------------------
NB. Verb:         Solves:       Syntax:
NB. getrspl1ux    A   * X = B   Xv=. (ip;L1U) getrspl1ux  Bv
NB. getrspl1uhx   A^H * X = B   Xv=. (ip;L1U) getrspl1uhx Bv
NB. getrspl1utx   A^T * X = B   Xv=. (ip;L1U) getrspl1utx Bv
NB. getrsxpl1u    X * A   = B   Xh=. (ip;L1U) getrsxpl1u  Bh
NB. getrsxpl1uh   X * A^H = B   Xh=. (ip;L1U) getrsxpl1uh Bh
NB. getrsxpl1ut   X * A^T = B   Xh=. (ip;L1U) getrsxpl1ut Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general square
NB.   matrix A, represented in factored form:
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
NB.   nrhs ≥ 0, number of RHSs
NB.
NB. Notes:
NB. - implements LAPACK's xGETRS

getrspl1ux=:  (1 {:: [) ([ trsmux  trsml1x ) ((0 {:: [) C.       ])
getrsxpl1uh=: (1 {:: [) ([ trsmxuh trsmxl1h) ((0 {:: [) C.^:_1"1 ])
getrsxpl1ut=: (1 {:: [) ([ trsmxut trsmxl1t) ((0 {:: [) C.^:_1"1 ])

getrspl1uhx=: (0 {:: [) C.   ((] trsml1hx trsmuhx~) (1 & {::))~
getrspl1utx=: (0 {:: [) C.   ((] trsml1tx trsmutx~) (1 & {::))~
getrsxpl1u=:  (0 {:: [) C."1 ((] trsmxl1  trsmxu ~) (1 & {::))~

NB. ---------------------------------------------------------
NB. Verb:         Solves:       Syntax:
NB. getrspu1lx    A   * X = B   Xv=. (ip;U1L) getrspu1lx  Bv
NB. getrspu1lhx   A^H * X = B   Xv=. (ip;U1L) getrspu1lhx Bv
NB. getrspu1ltx   A^T * X = B   Xv=. (ip;U1L) getrspu1ltx Bv
NB. getrsxpu1l    X * A   = B   Xh=. (ip;U1L) getrsxpu1l  Bh
NB. getrsxpu1lh   X * A^H = B   Xh=. (ip;U1L) getrsxpu1lh Bh
NB. getrsxpu1lt   X * A^T = B   Xh=. (ip;U1L) getrsxpu1lt Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general square
NB.   matrix A, represented in factored form:
NB.     P * U1 * L = A
NB. where:
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   ip   - n-vector, rows inversed permutation of A
NB.   U1L  - n×n-matrix, lower triangle contains L, and
NB.          strict upper triangle contains U1 without unit
NB.          diagonal
NB.   P    - n×n-matrix, rows permutation of A
NB.   L    - n×n-matrix, lower triangular
NB.   U1   - n×n-matrix, unit upper triangular
NB.   nrhs ≥ 0, number of RHSs

getrspu1lx=:  (1 {:: [) ([ trsmlx  trsmu1x ) ((0 {:: [) C.       ])
getrsxpu1lh=: (1 {:: [) ([ trsmxlh trsmxu1h) ((0 {:: [) C.^:_1"1 ])
getrsxpu1lt=: (1 {:: [) ([ trsmxlt trsmxu1t) ((0 {:: [) C.^:_1"1 ])

getrspu1lhx=: (0 {:: [) C.   ((] trsmu1hx trsmlhx~) (1 & {::))~
getrspu1ltx=: (0 {:: [) C.   ((] trsmu1tx trsmltx~) (1 & {::))~
getrsxpu1l=:  (0 {:: [) C."1 ((] trsmxu1  trsmxl ~) (1 & {::))~

NB. ---------------------------------------------------------
NB. Verb:         Solves:       Syntax:
NB. getrsul1px    A   * X = B   Xv=. (ip;UL1) getrsul1px  Bv
NB. getrsul1phx   A^H * X = B   Xv=. (ip;UL1) getrsul1phx Bv
NB. getrsul1ptx   A^T * X = B   Xv=. (ip;UL1) getrsul1ptx Bv
NB. getrsxul1p    X * A   = B   Xh=. (ip;UL1) getrsxul1p  Bh
NB. getrsxul1ph   X * A^H = B   Xh=. (ip;UL1) getrsxul1ph Bh
NB. getrsxul1pt   X * A^T = B   Xh=. (ip;UL1) getrsxul1pt Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general square
NB.   matrix A, represented in factored form:
NB.     U * L1 * P = A
NB. where:
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   ip   - n-vector, columns inversed permutation of A
NB.   UL1  - n×n-matrix, upper triangle contains U, and
NB.          strict lower triangle contains L1 without unit
NB.          diagonal
NB.   P    - n×n-matrix, columns permutation of A
NB.   L1   - n×n-matrix, unit lower triangular
NB.   U    - n×n-matrix, upper triangular
NB.   nrhs ≥ 0, number of RHSs

getrsul1px=:  (0 {:: [) C.       ((] trsml1x  trsmux ~) (1 & {::))~
getrsxul1ph=: (0 {:: [) C.^:_1"1 ((] trsmxl1h trsmxuh~) (1 & {::))~
getrsxul1pt=: (0 {:: [) C.^:_1"1 ((] trsmxl1t trsmxut~) (1 & {::))~

getrsul1phx=: (1 {:: [) ([ trsmuhx trsml1hx) ((0 {:: [) C.^:_1   ])
getrsul1ptx=: (1 {:: [) ([ trsmutx trsml1tx) ((0 {:: [) C.^:_1   ])
getrsxul1p=:  (1 {:: [) ([ trsmxu  trsmxl1 ) ((0 {:: [) C.       ])

NB. ---------------------------------------------------------
NB. Verb:         Solves:       Syntax:
NB. hetrsplx      A   * X = B   Xv=. (ip;L1;T) hetrsplx  Bv
NB. hetrspltx     A^T * X = B   Xv=. (ip;L1;T) hetrspltx Bv
NB. hetrsxpl      X * A   = B   Xh=. (ip;L1;T) hetrsxpl  Bh
NB. hetrsxplt     X * A^T = B   Xh=. (ip;L1;T) hetrsxplt Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) matrix A, represented in factored form:
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
NB.   T    - n×n-matrix, Hermitian (symmetric) tridiagonal
NB.   nrhs ≥ 0, number of RHSs
NB.
NB. Notes:
NB. - is similar to LAPACK's xHETRS, but uses another
NB.   factorization, see hetrfx

hetrsplx=:   (0 {:: ])    C.^:_1  ((1 {:: ]) trsml1hx ((2 {:: [) pttrsax  ((1 {:: [) trsml1x  ((0 {:: [) C.       ]))))
hetrspltx=: ((0 {:: ]) +@(C.^:_1) ((1 {:: ]) trsml1hx ((2 {:: [) pttrsax  ((1 {:: [) trsml1x  ((0 {:: [) C.       ]))))) +
hetrsxpl=:   (0 {:: ])    C."1    ((1 {:: ]) trsmxl1  ((2 {:: ]) pttrsxa  ((1 {:: ]) trsmxl1h ((0 {:: [) C.^:_1"1 ]))))
hetrsxplt=: ((0 {:: ])    C."1    ((1 {:: ]) trsmxl1  ((2 {:: ]) pttrsxa  ((1 {:: ]) trsmxl1h ((0 {:: [) C.^:_1"1 ]))))) +

NB. ---------------------------------------------------------
NB. Verb:         Solves:       Syntax:
NB. hetrspux      A   * X = B   Xv=. (ip;U1;T) hetrspux  Bv
NB. hetrsputx     A^T * X = B   Xv=. (ip;U1;T) hetrsputx Bv
NB. hetrsxpu      X * A   = B   Xh=. (ip;U1;T) hetrsxpu  Bh
NB. hetrsxput     X * A^T = B   Xh=. (ip;U1;T) hetrsxput Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) matrix A, represented in factored form:
NB.     P * U1 * T * U1^H * P^_1 = A
NB. where:
NB.   A    - n×n-matrix, Hermitian (symmetric)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   ip   - n-vector, full inversed permutation of A
NB.   P    - n×n-matrix, full permutation of A
NB.   U1   - n×n-matrix, unit upper triangular
NB.   T    - n×n-matrix, Hermitian (symmetric) tridiagonal
NB.   nrhs ≥ 0

hetrspux=:   (0 {:: ])    C.^:_1  ((1 {:: ]) trsmu1hx ((2 {:: [) pttrsax  ((1 {:: [) trsmu1x  ((0 {:: [) C.       ]))))
hetrsputx=: ((0 {:: ]) +@(C.^:_1) ((1 {:: ]) trsmu1hx ((2 {:: [) pttrsax  ((1 {:: [) trsmu1x  ((0 {:: [) C.       ]))))) +
hetrsxpu=:   (0 {:: ])    C."1    ((1 {:: ]) trsmxu1  ((2 {:: ]) pttrsxa  ((1 {:: ]) trsmxu1h ((0 {:: [) C.^:_1"1 ]))))
hetrsxput=: ((0 {:: ])    C."1    ((1 {:: ]) trsmxu1  ((2 {:: ]) pttrsxa  ((1 {:: ]) trsmxu1h ((0 {:: [) C.^:_1"1 ]))))) +

NB. ---------------------------------------------------------
NB. Verb:         Solves:       Syntax:
NB. potrslx       A   * X = B   Xv=. L potrslx  Bv
NB. potrsltx      A^T * X = B   Xv=. L potrsltx Bv
NB. potrsxl       X * A   = B   Xh=. L potrsxl  Bh
NB. potrsxlt      X * A^T = B   Xh=. L potrsxlt Bh
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
NB.   nrhs ≥ 0, number of RHSs
NB.
NB. Notes:
NB. - implements LAPACK's xPOTRS

potrslx=:   [   trsmlhx trsmlx
potrsltx=: ([ +@trsmlhx trsmlx ) +
potrsxl=:   [   trsmxl  trsmxlh
potrsxlt=: ([ +@trsmxl  trsmxlh) +

NB. ---------------------------------------------------------
NB. Verb:         Solves:       Syntax:
NB. potrsux       A   * X = B   Xv=. U potrsux  Bv
NB. potrsutx      A^T * X = B   Xv=. U potrsutx Bv
NB. potrsxu       X * A   = B   Xh=. U potrsxu  Bh
NB. potrsxut      X * A^T = B   Xh=. U potrsxut Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) positive definite matrix A, represented in
NB.   factored form:
NB.     U * U^H = A
NB. where:
NB.   A    - n×n Hermitian (symmetric) positive definite
NB.          matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   U    - n×n-matrix, upper triangular with positive
NB.          diagonal entries, Cholesky triangle
NB.   nrhs ≥ 0, number of RHSs

potrsux=:   [   trsmuhx trsmux
potrsutx=: ([ +@trsmuhx trsmux ) +
potrsxu=:   [   trsmxu  trsmxuh
potrsxut=: ([ +@trsmxu  trsmxuh) +

NB. ---------------------------------------------------------
NB. Verb:         Solves:       Syntax:
NB. pttrslx       A   * X = B   Xv=. (L1;D) pttrslx  Bv
NB. pttrsltx      A^T * X = B   Xv=. (L1;D) pttrsltx Bv
NB. pttrsxl       X * A   = B   Xh=. (L1;D) pttrsxl  Bh
NB. pttrsxlt      X * A^T = B   Xh=. (L1;D) pttrsxlt Bh
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
NB. Algorithm for pttrslx:
NB.   In:  L1 D Bv
NB.   Out: Xv
NB.   0) extract main diagonal d from D and subdiagonal e
NB.      from L1
NB.   1) prepare input:
NB.        be=. Bv ,. (0,e)
NB.   2) do iterations k=1:n-1 by reversed suffix scan:
NB.        btrash=. u/\.&.|. be
NB.      to find :
NB.        b[k] := b[k] - b[k-1]*e[k-1]
NB.   3) cut off trash column to extract updated Bv:
NB.        b=. (}:"1) btrash
NB.   4) prepare intermediate input:
NB.        bde=. ((}: b) , (({: b) % ({: d))) ,. d ,. ((conj(e),0)
NB.   5) do iterations k=n-2:0 by suffix scan:
NB.        btrash=. u/\. bde
NB.      to find :
NB.        b[k] := b[k]/d[k] - b[k+1]*e[k]
NB.   6) cut off two last columns of trash to extract raw Xv
NB.      and re-shape to Bv's shape:
NB.        Xv=. ($ Bv) ($,) _2 }."1 btrash
NB.
NB. Assertions:
NB.   Xv -: clean L1D pttrslx Bv
NB. where
NB.   L1D=. pttrfl A
NB.   Bv=. A mp Xv
NB.
NB. Notes:
NB. - implements LAPACK's xPTTS2(0)
NB. - if A is singular then solution Xx will be wrong
NB. - if A is indefinite then solution Xx may be wrong
NB.
NB. References:
NB. [1] G. H. Golub and C. F. Van Loan, Matrix Computations,
NB.     Johns Hopkins University Press, Baltimore, Md, USA,
NB.     3rd edition, 1996, p. 157

pttrslx=:  $@] ($,) (_2 }."1 (   ]`((}:"1)@((}:@[(-,0:)((* {:)~ }:))/\.&.|.)@(,. (0&,)))`(_1   diag (0 {:: [))`(((_2 (}.%{)[)(-,0 0"_)((* {:)~ (_2&}.)))/\. @ (_2 ({:@] % ({,))`_1:`]} ,.))`((,. (+,0:))~)`(diag@(1 {:: [)) fork3))
pttrsltx=: $@] ($,) (_2 }."1 (   ]`((}:"1)@((}:@[(-,0:)((* {:)~ }:))/\.&.|.)@(,. (0&,)))`(_1 +@diag (0 {:: [))`(((_2 (}.%{)[)(-,0 0"_)((* {:)~ (_2&}.)))/\. @ (_2 ({:@] % ({,))`_1:`]} ,.))`((,. (+,0:))~)`(diag@(1 {:: [)) fork3))  NB. pttrsltx=: (<@:+@(0 {:: [))`0:`[} pttrslx ]
pttrsxl=:  $@] ($,) (_2 }."1 (|:@]`((}:"1)@((}:@[(-,0:)((* {:)~ }:))/\.&.|.)@(,. (0&,)))`(_1 +@diag (0 {:: [))`(((_2 (}.%{)[)(-,0 0"_)((* {:)~ (_2&}.)))/\. @ (_2 ({:@] % ({,))`_1:`]} ,.))`((,. (+,0:))~)`(diag@(1 {:: [)) fork3))  NB. pttrsxl=: pttrsltx |:
pttrsxlt=: $@] ($,) (_2 }."1 (|:@]`((}:"1)@((}:@[(-,0:)((* {:)~ }:))/\.&.|.)@(,. (0&,)))`(_1   diag (0 {:: [))`(((_2 (}.%{)[)(-,0 0"_)((* {:)~ (_2&}.)))/\. @ (_2 ({:@] % ({,))`_1:`]} ,.))`((,. (+,0:))~)`(diag@(1 {:: [)) fork3))  NB. pttrsxlt=: pttrslx |:

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
  'LU1ip ipL1U ipU1L UL1ip'=. (getrflu1p ; getrfpl1u ; getrfpu1l ; < @getrful1p) A
##############
  ('getrspl1ux'  tdyad ((_2&{.)`((mp  & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (    A ;X;ipL1U)
  ('getrspl1uhx' tdyad ((_2&{.)`((mp  & >/)@(2&{.))`]`(conAh"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ((mp~ ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((ct A);X;ipL1U)
  ('getrspl1utx' tdyad ((_2&{.)`((mp  & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1@|:)) [) (1 & {::))~))`(normi@((norm1t"1@|:@(((mp  & >/)@(2 {. [)) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) ((|: A);X;ipL1U)
  ('getrsxpl1u'  tdyad ((_2&{.)`((mp~ & >/)@(2&{.))`]`(conA "_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (    A ;X;ipL1U)
  ('getrsxpl1uh' tdyad ((_2&{.)`((mp~ & >/)@(2&{.))`]`(conAh"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((ct A);X;ipL1U)
  ('getrsxpl1ut' tdyad ((_2&{.)`((mp~ & >/)@(2&{.))`]`(conAt"_)`(normi@(((- (% & (normi"1   )) [) (1 & {::))~))`(normi@((norm1t"1   @(((mp~ & >/)@(2 {. [)) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) ((|: A);X;ipL1U)

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
  'conA conAt'=. 2 # _. NB. (norm1 con (potril@potrfl))"2 (] ,: |:) A
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
  'conA conAt'=. 2 # _. NB. (norm1 con (pttril@pttrfl))"2 (] ,: |:) A NB. ##################
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
