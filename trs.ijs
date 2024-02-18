NB. Solve linear monomial equation by triangular
NB. factorization
NB.
NB. getrsxxxxxx  Solve equation (op(A) * X = B) or
NB.              (X * op(A) = B), where A is a general square
NB.              matrix, represented in factored form; op(A)
NB.              is either A itself, or A^T (the
NB.              transposition of A), or A^H (the conjugate
NB.              transposition of A); B is known right-hand
NB.              sides (RHS), X is unknown solutions
NB. hetrsxxxx    Solve equation (op(A) * X = B) or
NB.              (X * op(A) = B), where A is a Hermitian
NB.              (symmetric) matrix, represented in factored
NB.              form; op(A) is either A itself, or A^T (the
NB.              transposition of A); B is known right-hand
NB.              sides (RHS), X is unknown solutions
NB. potrsxxx     Solve equation (op(A) * X = B) or
NB.              (X * op(A) = B), where A is a Hermitian
NB.              (symmetric) positive definite matrix,
NB.              represented in factored form; op(A) is
NB.              either A itself, or A^T (the transposition
NB.              of A); B is known right-hand sides (RHS), X
NB.              is unknown solutions
NB. pttrsxxx     Solve equation (op(A) * X = B) or
NB.              (X * op(A) = B), where A is a Hermitian
NB.              (symmetric) positive definite tridiagonal
NB.              matrix, represented in factored form; op(A)
NB.              is either A itself, or A^T (the
NB.              transposition of A); B is known right-hand
NB.              sides (RHS), X is unknown solutions
NB. trtrsxxx     Solve equation (op(A) * X = B) or
NB.              (X * op(A) = B), where A is either unit or
NB.              non-unit, either lower or upper, triangular
NB.              matrix; op(A) is either A itself, or A^T,
NB.              the transposition of A, or A^H, the conjugate
NB.              transposition of A; B is known right-hand
NB.              sides (RHS), X is unknown solutions
NB.
NB. testgetrs    Test getrsxxxxxx by general square matrix
NB. testhetrs    Test hetrsxxxx by Hermitian (symmetric)
NB.              matrix
NB. testpotrs    Test potrsxxx by Hermitian (symmetric)
NB.              positive definite matrix
NB. testpttrs    Test pttrsxxx by Hermitian (symmetric)
NB.              positive definite tridiagonal matrix
NB. testtrtrs    Test trtrsxxx by triangular matrix
NB. testtrs      Adv. to make verb to test xxtrsxxxxxx by
NB.              matrix of generator and shape given
NB.
NB. Version: 0.13.2 2021-06-24
NB.
NB. Copyright 2010-2021 Igor Zhuravlov
NB.
NB. This file is part of mt
NB.
NB. mt is free software: you can redistribute it and/or
NB. modify it under the terms of the GNU Lesser General
NB. Public License as published by the Free Software
NB. Foundation, either version 3 of the License, or (at your
NB. option) any later version.
NB.
NB. mt is distributed in the hope that it will be useful, but
NB. WITHOUT ANY WARRANTY; without even the implied warranty
NB. of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
NB. See the GNU Lesser General Public License for more
NB. details.
NB.
NB. You should have received a copy of the GNU Lesser General
NB. Public License along with mt. If not, see
NB. <http://www.gnu.org/licenses/>.

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb           Solves         Syntax
NB. getrslu1px     A   * X = B    Xv=. LU1p getrslu1px  Bv
NB. getrslu1pcx    A^H * X = B    Xv=. LU1p getrslu1pcx Bv
NB. getrslu1ptx    A^T * X = B    Xv=. LU1p getrslu1ptx Bv
NB. getrsxlu1p     X * A   = B    Xh=. LU1p getrsxlu1p  Bh
NB. getrsxlu1pc    X * A^H = B    Xh=. LU1p getrsxlu1pc Bh
NB. getrsxlu1pt    X * A^T = B    Xh=. LU1p getrsxlu1pt Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general square
NB.   matrix A, represented in factored form:
NB.     L * U1 * P = A
NB. where
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   LU1p - 2-vector of boxes, the output of getrflu1p, the
NB.          matrix A represented in factored form
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   Xv -: clean LU1p getrslu1px      A  mp  Xv
NB.   Xv -: clean LU1p getrslu1pcx (ct A) mp  Xv
NB.   Xv -: clean LU1p getrslu1ptx (|: A) mp  Xv
NB.   Xh -: clean LU1p getrsxlu1p      A  mp~ Xh
NB.   Xh -: clean LU1p getrsxlu1pc (ct A) mp~ Xh
NB.   Xh -: clean LU1p getrsxlu1pt (|: A) mp~ Xh
NB. where
NB.   LU1p=. getrflu1p A

getrslu1px=:  (0 {:: [) C.^:_1   ((] trsmlunu trsmllnn~) 1&{::)~
getrsxlu1pc=: (0 {:: [) C.^:_1"1 ((] trsmrucu trsmrlcn~) 1&{::)~
getrsxlu1pt=: (0 {:: [) C.^:_1"1 ((] trsmrutu trsmrltn~) 1&{::)~

getrslu1pcx=: (1 {:: [) ([ trsmllcn trsmlucu) ((0 {:: [) C.   ])
getrslu1ptx=: (1 {:: [) ([ trsmlltn trsmlutu) ((0 {:: [) C.   ])
getrsxlu1p=:  (1 {:: [) ([ trsmrlnn trsmrunu) ((0 {:: [) C."1 ])

NB. ---------------------------------------------------------
NB. Verb           Solves         Syntax
NB. getrspl1ux     A   * X = B    Xv=. pL1U getrspl1ux  Bv
NB. getrspl1ucx    A^H * X = B    Xv=. pL1U getrspl1ucx Bv
NB. getrspl1utx    A^T * X = B    Xv=. pL1U getrspl1utx Bv
NB. getrsxpl1u     X * A   = B    Xh=. pL1U getrsxpl1u  Bh
NB. getrsxpl1uc    X * A^H = B    Xh=. pL1U getrsxpl1uc Bh
NB. getrsxpl1ut    X * A^T = B    Xh=. pL1U getrsxpl1ut Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general square
NB.   matrix A, represented in factored form:
NB.     P * L1 * U = A
NB. where
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   pL1U - 2-vector of boxes, the output of getrfpl1u, the
NB.          matrix A represented in factored form
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   Xv -: clean pL1U getrspl1ux      A  mp  Xv
NB.   Xv -: clean pL1U getrspl1ucx (ct A) mp  Xv
NB.   Xv -: clean pL1U getrspl1utx (|: A) mp  Xv
NB.   Xh -: clean pL1U getrsxpl1u      A  mp~ Xh
NB.   Xh -: clean pL1U getrsxpl1uc (ct A) mp~ Xh
NB.   Xh -: clean pL1U getrsxpl1ut (|: A) mp~ Xh
NB. where
NB.   pL1U=. getrfpl1u A
NB.
NB. Notes:
NB. - getrspxxxX implement LAPACK's xGETRS

getrspl1ux=:  (1 {:: [) ([ trsmlunn trsmllnu) ((0 {:: [) C.   ])
getrsxpl1uc=: (1 {:: [) ([ trsmrucn trsmrlcu) ((0 {:: [) C."1 ])
getrsxpl1ut=: (1 {:: [) ([ trsmrutn trsmrltu) ((0 {:: [) C."1 ])

getrspl1ucx=: (0 {:: [) C.^:_1   ((] trsmllcu trsmlucn~) 1&{::)~
getrspl1utx=: (0 {:: [) C.^:_1   ((] trsmlltu trsmlutn~) 1&{::)~
getrsxpl1u=:  (0 {:: [) C.^:_1"1 ((] trsmrlnu trsmrunn~) 1&{::)~

NB. ---------------------------------------------------------
NB. Verb           Solves         Syntax
NB. getrspu1lx     A   * X = B    Xv=. pU1L getrspu1lx  Bv
NB. getrspu1lcx    A^H * X = B    Xv=. pU1L getrspu1lcx Bv
NB. getrspu1ltx    A^T * X = B    Xv=. pU1L getrspu1ltx Bv
NB. getrsxpu1l     X * A   = B    Xh=. pU1L getrsxpu1l  Bh
NB. getrsxpu1lc    X * A^H = B    Xh=. pU1L getrsxpu1lc Bh
NB. getrsxpu1lt    X * A^T = B    Xh=. pU1L getrsxpu1lt Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general square
NB.   matrix A, represented in factored form:
NB.     P * U1 * L = A
NB. where
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   pU1L - 2-vector of boxes, the output of getrfpu1l, the
NB.          matrix A represented in factored form
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   Xv -: clean pU1L getrspu1lx      A  mp  Xv
NB.   Xv -: clean pU1L getrspu1lcx (ct A) mp  Xv
NB.   Xv -: clean pU1L getrspu1ltx (|: A) mp  Xv
NB.   Xh -: clean pU1L getrsxpu1l      A  mp~ Xh
NB.   Xh -: clean pU1L getrsxpu1lc (ct A) mp~ Xh
NB.   Xh -: clean pU1L getrsxpu1lt (|: A) mp~ Xh
NB. where
NB.   pU1L=. getrfpu1l A

getrspu1lx=:  (1 {:: [) ([ trsmllnn trsmlunu) ((0 {:: [) C.   ])
getrsxpu1lc=: (1 {:: [) ([ trsmrlcn trsmrucu) ((0 {:: [) C."1 ])
getrsxpu1lt=: (1 {:: [) ([ trsmrltn trsmrutu) ((0 {:: [) C."1 ])

getrspu1lcx=: (0 {:: [) C.^:_1   ((] trsmlucu trsmllcn~) 1&{::)~
getrspu1ltx=: (0 {:: [) C.^:_1   ((] trsmlutu trsmlltn~) 1&{::)~
getrsxpu1l=:  (0 {:: [) C.^:_1"1 ((] trsmrunu trsmrlnn~) 1&{::)~

NB. ---------------------------------------------------------
NB. Verb           Solves         Syntax
NB. getrsul1px     A   * X = B    Xv=. UL1p getrsul1px  Bv
NB. getrsul1pcx    A^H * X = B    Xv=. UL1p getrsul1pcx Bv
NB. getrsul1ptx    A^T * X = B    Xv=. UL1p getrsul1ptx Bv
NB. getrsxul1p     X * A   = B    Xh=. UL1p getrsxul1p  Bh
NB. getrsxul1pc    X * A^H = B    Xh=. UL1p getrsxul1pc Bh
NB. getrsxul1pt    X * A^T = B    Xh=. UL1p getrsxul1pt Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general square
NB.   matrix A, represented in factored form:
NB.     U * L1 * P = A
NB. where
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   UL1p - 2-vector of boxes, the output of getrful1p, the
NB.          matrix A represented in factored form
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   Xv -: clean UL1p getrsul1px      A  mp  Xv
NB.   Xv -: clean UL1p getrsul1pcx (ct A) mp  Xv
NB.   Xv -: clean UL1p getrsul1ptx (|: A) mp  Xv
NB.   Xh -: clean UL1p getrsxul1p      A  mp~ Xh
NB.   Xh -: clean UL1p getrsxul1pc (ct A) mp~ Xh
NB.   Xh -: clean UL1p getrsxul1pt (|: A) mp~ Xh
NB. where
NB.   UL1p=. getrful1p A

getrsul1px=:  (0 {:: [) C.^:_1   ((] trsmllnu trsmlunn~) 1&{::)~
getrsxul1pc=: (0 {:: [) C.^:_1"1 ((] trsmrlcu trsmrucn~) 1&{::)~
getrsxul1pt=: (0 {:: [) C.^:_1"1 ((] trsmrltu trsmrutn~) 1&{::)~

getrsul1pcx=: (1 {:: [) ([ trsmlucn trsmllcu) ((0 {:: [) C.   ])
getrsul1ptx=: (1 {:: [) ([ trsmlutn trsmlltu) ((0 {:: [) C.   ])
getrsxul1p=:  (1 {:: [) ([ trsmrunn trsmrlnu) ((0 {:: [) C."1 ])

NB. ---------------------------------------------------------
NB. Verb         Solves         Syntax
NB. hetrsplx     A   * X = B    Xv=. pL1T hetrsplx  Bv
NB. hetrspltx    A^T * X = B    Xv=. pL1T hetrspltx Bv
NB. hetrsxpl     X * A   = B    Xh=. pL1T hetrsxpl  Bh
NB. hetrsxplt    X * A^T = B    Xh=. pL1T hetrsxplt Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) matrix A, represented in factored form:
NB.     P * L1 * T * L1^H * P^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   pL1T - 3-vector of boxes, the output of hetrfpl, the
NB.          matrix A represented in factored form
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   Xv -: clean pL1T hetrsplx      A  mp  Xv
NB.   Xv -: clean pL1T hetrspltx (|: A) mp  Xv
NB.   Xh -: clean pL1T hetrsxpl      A  mp~ Xh
NB.   Xh -: clean pL1T hetrsxplt (|: A) mp~ Xh
NB. where
NB.   pL1T=. hetrfpl A
NB.
NB. Notes:
NB. - models LAPACK's DSYTRS_AA('L') and ZHETRS_AA('L')

hetrsplx=:   (0 {:: [)    C.^:_1    (1 {:: [) trsmllcu (2 {:: [) gtsvax (1 {:: [) trsmllnu (0 {:: [) C.   ]
hetrspltx=: ((0 {:: [) +@(C.^:_1  ) (1 {:: [) trsmllcu (2 {:: [) gtsvax (1 {:: [) trsmllnu (0 {:: [) C.   ]) +
hetrsxpl=:   (0 {:: [)    C.^:_1"1  (1 {:: [) trsmrlnu (2 {:: [) gtsvxa (1 {:: [) trsmrlcu (0 {:: [) C."1 ]
hetrsxplt=: ((0 {:: [) +@(C.^:_1"1) (1 {:: [) trsmrlnu (2 {:: [) gtsvxa (1 {:: [) trsmrlcu (0 {:: [) C."1 ]) +

NB. ---------------------------------------------------------
NB. Verb         Solves         Syntax
NB. hetrspux     A   * X = B    Xv=. pU1T hetrspux  Bv
NB. hetrsputx    A^T * X = B    Xv=. pU1T hetrsputx Bv
NB. hetrsxpu     X * A   = B    Xh=. pU1T hetrsxpu  Bh
NB. hetrsxput    X * A^T = B    Xh=. pU1T hetrsxput Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) matrix A, represented in factored form:
NB.     P * U1 * T * U1^H * P^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   pU1T - 3-vector of boxes, the output of hetrfpu, the
NB.          matrix A represented in factored form
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   Xv -: clean pU1T hetrspux      A  mp  Xv
NB.   Xv -: clean pU1T hetrsputx (|: A) mp  Xv
NB.   Xh -: clean pU1T hetrsxpu      A  mp~ Xh
NB.   Xh -: clean pU1T hetrsxput (|: A) mp~ Xh
NB. where
NB.   pU1T=. hetrfpu A

hetrspux=:   (0 {:: [)    C.^:_1    (1 {:: [) trsmlucu (2 {:: [) gtsvax (1 {:: [) trsmlunu (0 {:: [) C.   ]
hetrsputx=: ((0 {:: [) +@(C.^:_1  ) (1 {:: [) trsmlucu (2 {:: [) gtsvax (1 {:: [) trsmlunu (0 {:: [) C.   ]) +
hetrsxpu=:   (0 {:: [)    C.^:_1"1  (1 {:: [) trsmrunu (2 {:: [) gtsvxa (1 {:: [) trsmrucu (0 {:: [) C."1 ]
hetrsxput=: ((0 {:: [) +@(C.^:_1"1) (1 {:: [) trsmrunu (2 {:: [) gtsvxa (1 {:: [) trsmrucu (0 {:: [) C."1 ]) +

NB. ---------------------------------------------------------
NB. Verb        Solves         Syntax
NB. potrslx     A   * X = B    Xv=. L potrslx  Bv
NB. potrsltx    A^T * X = B    Xv=. L potrsltx Bv
NB. potrsxl     X * A   = B    Xh=. L potrsxl  Bh
NB. potrsxlt    X * A^T = B    Xh=. L potrsxlt Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) positive definite matrix A, represented in
NB.   factored form:
NB.     L * L^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   L    - n×n-matrix, the lower triangular with positive
NB.          diagonal entries, the Cholesky triangle, the
NB.          output of potrfl, the matrix A represented in
NB.          factored form
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   Xv -: clean L potrslx      A  mp  Xv
NB.   Xv -: clean L potrsltx (|: A) mp  Xv
NB.   Xh -: clean L potrsxl      A  mp~ Xh
NB.   Xh -: clean L potrsxlt (|: A) mp~ Xh
NB. where
NB.   L=. potrfl A
NB.
NB. Notes:
NB. - implements LAPACK's xPOTRS('L')

potrslx=:   [   trsmllcn trsmllnn
potrsltx=: ([ +@trsmllcn trsmllnn) +
potrsxl=:   [   trsmrlnn trsmrlcn
potrsxlt=: ([ +@trsmrlnn trsmrlcn) +

NB. ---------------------------------------------------------
NB. Verb        Solves         Syntax
NB. potrsux     A   * X = B    Xv=. U potrsux  Bv
NB. potrsutx    A^T * X = B    Xv=. U potrsutx Bv
NB. potrsxu     X * A   = B    Xh=. U potrsxu  Bh
NB. potrsxut    X * A^T = B    Xh=. U potrsxut Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) positive definite matrix A, represented in
NB.   factored form:
NB.     U * U^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   U    - n×n-matrix, the upper triangular with positive
NB.          diagonal entries, the Cholesky triangle, the
NB.          output of potrfu, the matrix A represented in
NB.          factored form
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   Xv -: clean U potrsux      A  mp  Xv
NB.   Xv -: clean U potrsutx (|: A) mp  Xv
NB.   Xh -: clean U potrsxu      A  mp~ Xh
NB.   Xh -: clean U potrsxut (|: A) mp~ Xh
NB. where
NB.   U=. potrfu A

potrsux=:   [   trsmlucn trsmlunn
potrsutx=: ([ +@trsmlucn trsmlunn) +
potrsxu=:   [   trsmrunn trsmrucn
potrsxut=: ([ +@trsmrunn trsmrucn) +

NB. ---------------------------------------------------------
NB. Verb        Solves         Syntax
NB. pttrslx     A   * X = B    Xv=. L1D pttrslx  Bv
NB. pttrsltx    A^T * X = B    Xv=. L1D pttrsltx Bv
NB. pttrsxl     X * A   = B    Xh=. L1D pttrsxl  Bh
NB. pttrsxlt    X * A^T = B    Xh=. L1D pttrsxlt Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) positive definite tridiagonal matrix A,
NB.   represented in factored form [1]:
NB.     L1 * D * L1^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite tridiagonal
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   L1D  - 2-vector of boxes, the output of pttrfl, the
NB.          matrix A represented in factored form
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Algorithm for pttrslx:
NB.   In:  L1 D Bv
NB.   Out: Xv
NB.   1) extract main diagonal d from D and subdiagonal e
NB.      from L1
NB.   2) prepare input:
NB.        be=. Bv ,. (0,e)
NB.   3) do iterations k=1:n-1 by reversed suffix scan:
NB.        btrash=. u/\.&.|. be
NB.      to find :
NB.        b[k] := b[k] - b[k-1]*e[k-1]
NB.      for non-empty be only
NB.   4) cut off trash column to extract updated Bv:
NB.        b=. (}:"1) btrash
NB.   5) prepare intermediate input:
NB.        bde=. ((}: b) , (({: b) % ({: d))) ,. d ,. (conj(e),0)
NB.   6) do iterations k=n-2:0 by suffix scan:
NB.        btrash=. u/\. bde
NB.      to find :
NB.        b[k] := b[k]/d[k] - b[k+1]*conj(e[k])
NB.   7) cut off two last columns of trash to extract raw Xv
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
NB. [1] G. H. Golub, C. F. Van Loan. Matrix Computations.
NB.     Johns Hopkins University Press, Baltimore, Md, USA,
NB.     3rd edition, 1996, p. 157.
NB.
NB. TODO:
NB. - L1 and D would be sparse

pttrslx=:  $@] ($,) _2 }."1 ((diag@(1 {:: [))`(stitcht +)`(_1   diag (0 {:: [))`((((_2 (}. % {) [) (- , 0 0"_) ((* {:)~ _2&}.))/\.    @(_2  ({:@] % ({,))`_1:`]} ,.~))^:(0 < #@]))`(}:"1@((}:@[ (- , 0:) ((* {:)~ }:))/\.&.|.)@stitchb~)`] fork3)
pttrsltx=: $@] ($,) _2 }."1 ((diag@(1 {:: [))`(stitcht +)`(_1 +@diag (0 {:: [))`((((_2 (}. % {) [) (- , 0 0"_) ((* {:)~ _2&}.))/\.    @(_2  ({:@] % ({,))`_1:`]} ,.~))^:(0 < #@]))`(}:"1@((}:@[ (- , 0:) ((* {:)~ }:))/\.&.|.)@stitchb~)`] fork3)
pttrsxl=:  pttrsltx&.(a:`|:)
pttrsxlt=: pttrslx &.(a:`|:)

NB. ---------------------------------------------------------
NB. Verb        Solves         Syntax
NB. pttrsux     A   * X = B    Xv=. U1D pttrsux  Bv
NB. pttrsutx    A^T * X = B    Xv=. U1D pttrsutx Bv
NB. pttrsxu     X * A   = B    Xh=. U1D pttrsxu  Bh
NB. pttrsxut    X * A^T = B    Xh=. U1D pttrsxut Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) positive definite tridiagonal matrix A,
NB.   represented in factored form, based on [1]:
NB.     U1 * D * U1^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite tridiagonal
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   U1D  - 2-vector of boxes, the output of pttrfu, the
NB.          matrix A represented in factored form
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Algorithm for pttrsux:
NB.   In:  U1 D Bv
NB.   Out: Xv
NB.   1) extract main diagonal d from D and superdiagonal e
NB.      from U1
NB.   2) prepare input:
NB.        be=. Bv ,. (e,0)
NB.   3) do iterations k=n-2:0 by suffix scan:
NB.        btrash=. u/\. be
NB.      to find :
NB.        b[k] := b[k] - b[k+1]*e[k]
NB.      for non-empty be only
NB.   4) cut off trash column to extract updated Bv:
NB.        b=. (}:"1) btrash
NB.   5) prepare intermediate input:
NB.        bde=. ((({. b) % ({. d)) , (}. b)) ,. d ,. (0,conj(e))
NB.   6) do iterations k=1:n-1 by reversed suffix scan:
NB.        btrash=. u/\.&.|. bde
NB.      to find :
NB.        b[k] := b[k]/d[k] - b[k-1]*conj(e[k-1])
NB.   7) cut off two last columns of trash to extract raw Xv
NB.      and re-shape to Bv's shape:
NB.        Xv=. ($ Bv) ($,) _2 }."1 btrash
NB.
NB. Assertions:
NB.   Xv -: clean U1D pttrsux Bv
NB. where
NB.   U1D=. pttrfu A
NB.   Bv=. A mp Xv
NB.
NB. Notes:
NB. - if A is singular then solution Xx will be wrong
NB. - if A is indefinite then solution Xx may be wrong
NB.
NB. References:
NB. [1] G. H. Golub, C. F. Van Loan. Matrix Computations.
NB.     Johns Hopkins University Press, Baltimore, Md, USA,
NB.     3rd edition, 1996, p. 157.
NB.
NB. TODO:
NB. - U1 and D would be sparse

pttrsux=:  $@] ($,) _2 }."1 ((diag@(1 {:: [))`(stitchb +)`( 1   diag (0 {:: [))`((((_2 (}. % {) [) (- , 0 0"_) ((* {:)~ _2&}.))/\.&.|.@(c@] ({.@] % ({,))` 0:`]} ,.~))^:(0 < #@]))`(}:"1@((}:@[ (- , 0:) ((* {:)~ }:))/\.    )@stitcht~)`] fork3)
pttrsutx=: $@] ($,) _2 }."1 ((diag@(1 {:: [))`(stitchb +)`( 1 +@diag (0 {:: [))`((((_2 (}. % {) [) (- , 0 0"_) ((* {:)~ _2&}.))/\.&.|.@(c@] ({.@] % ({,))` 0:`]} ,.~))^:(0 < #@]))`(}:"1@((}:@[ (- , 0:) ((* {:)~ }:))/\.    )@stitcht~)`] fork3)
pttrsxu=:  pttrsutx&.(a:`|:)
pttrsxut=: pttrsux &.(a:`|:)

NB. ---------------------------------------------------------
NB. Verb        Reads in A    Solves          Syntax
NB. trtrslnn     LT           L    * X = B    X=. A trtrslnn B
NB. trtrslnu    SLT           L1   * X = B    X=. A trtrslnu B
NB. trtrsltn     LT           L^T  * X = B    X=. A trtrsltn B
NB. trtrsltu    SLT           L1^T * X = B    X=. A trtrsltu B
NB. trtrslcn     LT           L^H  * X = B    X=. A trtrslcn B
NB. trtrslcu    SLT           L1^H * X = B    X=. A trtrslcu B
NB. trtrsunn     UT           U    * X = B    X=. A trtrsunn B
NB. trtrsunu    SUT           U1   * X = B    X=. A trtrsunu B
NB. trtrsutn     UT           U^T  * X = B    X=. A trtrsutn B
NB. trtrsutu    SUT           U1^T * X = B    X=. A trtrsutu B
NB. trtrsucn     UT           U^H  * X = B    X=. A trtrsucn B
NB. trtrsucu    SUT           U1^H * X = B    X=. A trtrsucu B
NB.
NB. Description:
NB.   Solve the linear monomial matrix equation:
NB.     op(A) * X = B
NB. where
NB.   A - triangular and non-singular
NB.
NB. Syntax:
NB.   X=. B trtrsxxx A
NB. where
NB.   A  - m×m-matrix, contains either L, L1, U or U1
NB.        (unit diagonal is not stored)
NB.   B  - m×n-matrix or m-vector, the RHS
NB.   X  - m×n-matrix or m-vector, the solution(s)
NB.   m  ≥ 0, the size of A and the number of rows in B and X
NB.   n  ≥ 0, the number of columns in B and X
NB.
NB. Notes:
NB. - models LAPACK's xTRTRS when B and X are 2-rank

trtrslnn=: [: : ((trsmllnn~ ([ 'matrix is singular' assert 0 *./@:~: diag))~)
trtrslnu=: [: :   trsmllnu
trtrsltn=: [: : ((trsmlltn~ ([ 'matrix is singular' assert 0 *./@:~: diag))~)
trtrsltu=: [: :   trsmlltu
trtrslcn=: [: : ((trsmllcn~ ([ 'matrix is singular' assert 0 *./@:~: diag))~)
trtrslcu=: [: :   trsmllcu
trtrsunn=: [: : ((trsmlunn~ ([ 'matrix is singular' assert 0 *./@:~: diag))~)
trtrsunu=: [: :   trsmlunu
trtrsutn=: [: : ((trsmlutn~ ([ 'matrix is singular' assert 0 *./@:~: diag))~)
trtrsutu=: [: :   trsmlutu
trtrsucn=: [: : ((trsmlucn~ ([ 'matrix is singular' assert 0 *./@:~: diag))~)
trtrsucu=: [: :   trsmlucu

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgetrs
NB.
NB. Description:
NB.   Test:
NB.   - xGETRS (math/lapack2 addon)
NB.   - getrsxxxxxx (math/mt addon)
NB.   by general square matrix
NB.
NB. Syntax:
NB.   log=. testgetrs (X ; A)
NB. where
NB.   X   - m×n-matrix, exact solutions
NB.   A   - k×k-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - models LAPACK's xCHKGE

testgetrs=: 3 : 0
  load_mttmp_ 'math/mt/test/lapack2/getrf'
  load_mttmp_ 'math/mt/test/lapack2/getrs'

  'X A'=. y
  'm n'=. $ X
  Am=. (2 # m) {. A
  An=. (2 # n) {. A

  'rcondAm rcondAcm rcondAtm'=. (gecon1 , gecon1@ct , gecon1@|:) Am
  'rcondAn rcondAcn rcondAtn'=. (geconi , geconi@ct , geconi@|:) An

  'normiAtm norm1Atm'=. 'normiAcm norm1Acm'=. 'norm1Am normiAm'=. (norm1 , normi) Am
  'normiAtn norm1Atn'=. 'normiAcn norm1Acn'=. 'norm1An normiAn'=. (norm1 , normi) An

  'LU1pm pL1Um pU1Lm UL1pm'=. (getrflu1p ; getrfpl1u ; getrfpu1l ; <@getrful1p) Am
  'LU1pn pL1Un pU1Ln UL1pn'=. (getrflu1p ; getrfpl1u ; getrfpu1l ; <@getrful1p) An

  NB. matrix X

  Bax=.      Am  mp X
  Bacx=. (ct Am) mp X
  Batx=. (|: Am) mp X

  Bxa=.  X mp    An
  Bxac=. X mp ct An
  Bxat=. X mp |: An

  NB. vferrx
  vferrr=: normitc t04m  NB. for X at right side
  vferrl=: normitr t04m  NB. for X at left side

  NB. vberrxx for X at right side
  vberrax=:   mp~     t02m norm1tc
  vberracx=: (mp~ ct) t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc

  NB. vberrxx for Xh at left side
  vberrxa=:   mp      t02m norm1tr
  vberrxac=: (mp  ct) t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  log=.          ('''n''&dgetrs_mttmp_' tmonad (        (1&{ ,~ dgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondAm  ; norm1Am
  log=. log lcat ('''c''&dgetrs_mttmp_' tmonad (        (1&{ ,~ dgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrr`vberracx)) Am ; Bacx ; X  ; rcondAcm ; norm1Acm
  log=. log lcat ('''t''&dgetrs_mttmp_' tmonad (        (1&{ ,~ dgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondAtm ; norm1Atm

  log=. log lcat ('''n''&zgetrs_mttmp_' tmonad (        (1&{ ,~ zgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondAm  ; norm1Am
  log=. log lcat ('''c''&zgetrs_mttmp_' tmonad (        (1&{ ,~ zgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrr`vberracx)) Am ; Bacx ; X  ; rcondAcm ; norm1Acm
  log=. log lcat ('''t''&zgetrs_mttmp_' tmonad (        (1&{ ,~ zgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondAtm ; norm1Atm

  log=. log lcat ('getrslu1px'          tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondAm  ; norm1Am  ; LU1pm
  log=. log lcat ('getrslu1pcx'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberracx)) Am ; Bacx ; X  ; rcondAcm ; norm1Acm ; LU1pm
  log=. log lcat ('getrslu1ptx'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondAtm ; norm1Atm ; LU1pm
  log=. log lcat ('getrsxlu1p'          tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxa )) An ; Bxa  ; X  ; rcondAn  ; normiAn  ; LU1pn
  log=. log lcat ('getrsxlu1pc'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxac)) An ; Bxac ; X  ; rcondAcn ; normiAcn ; LU1pn
  log=. log lcat ('getrsxlu1pt'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxat)) An ; Bxat ; X  ; rcondAtn ; normiAtn ; LU1pn

  log=. log lcat ('getrspl1ux'          tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondAm  ; norm1Am  ; pL1Um
  log=. log lcat ('getrspl1ucx'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberracx)) Am ; Bacx ; X  ; rcondAcm ; norm1Acm ; pL1Um
  log=. log lcat ('getrspl1utx'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondAtm ; norm1Atm ; pL1Um
  log=. log lcat ('getrsxpl1u'          tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxa )) An ; Bxa  ; X  ; rcondAn  ; normiAn  ; pL1Un
  log=. log lcat ('getrsxpl1uc'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxac)) An ; Bxac ; X  ; rcondAcn ; normiAcn ; pL1Un
  log=. log lcat ('getrsxpl1ut'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxat)) An ; Bxat ; X  ; rcondAtn ; normiAtn ; pL1Un

  log=. log lcat ('getrspu1lx'          tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondAm  ; norm1Am  ; pU1Lm
  log=. log lcat ('getrspu1lcx'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberracx)) Am ; Bacx ; X  ; rcondAcm ; norm1Acm ; pU1Lm
  log=. log lcat ('getrspu1ltx'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondAtm ; norm1Atm ; pU1Lm
  log=. log lcat ('getrsxpu1l'          tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxa )) An ; Bxa  ; X  ; rcondAn  ; normiAn  ; pU1Ln
  log=. log lcat ('getrsxpu1lc'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxac)) An ; Bxac ; X  ; rcondAcn ; normiAcn ; pU1Ln
  log=. log lcat ('getrsxpu1lt'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxat)) An ; Bxat ; X  ; rcondAtn ; normiAtn ; pU1Ln

  log=. log lcat ('getrsul1px'          tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondAm  ; norm1Am  ; UL1pm
  log=. log lcat ('getrsul1pcx'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberracx)) Am ; Bacx ; X  ; rcondAcm ; norm1Acm ; UL1pm
  log=. log lcat ('getrsul1ptx'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondAtm ; norm1Atm ; UL1pm
  log=. log lcat ('getrsxul1p'          tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxa )) An ; Bxa  ; X  ; rcondAn  ; normiAn  ; UL1pn
  log=. log lcat ('getrsxul1pc'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxac)) An ; Bxac ; X  ; rcondAcn ; normiAcn ; UL1pn
  log=. log lcat ('getrsxul1pt'         tdyad  ((_2&{.)`(1&{::                       )`]`(3&{::)`vferrl`vberrxat)) An ; Bxat ; X  ; rcondAtn ; normiAtn ; UL1pn

  NB. vector x

  'xm xn'=. ({."1 ; {.) X

  bax=.      Am  mp xm
  bacx=. (ct Am) mp xm
  batx=. (|: Am) mp xm

  bxa=.  xn mp    An
  bxac=. xn mp ct An
  bxat=. xn mp |: An

  NB. vberrxx for x at right side
  vberrax=:   mp~     t02v
  vberracx=: (mp~ ct) t02v
  vberratx=: (mp~ |:) t02v

  NB. vberrxx for x at left side
  vberrxa=:   mp      t02v
  vberrxac=: (mp  ct) t02v
  vberrxat=: (mp  |:) t02v

  log=. log lcat ('getrslu1px'          tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrax )) Am ; bax  ; xm ; rcondAm  ; norm1Am  ; LU1pm
  log=. log lcat ('getrslu1pcx'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberracx)) Am ; bacx ; xm ; rcondAcm ; norm1Acm ; LU1pm
  log=. log lcat ('getrslu1ptx'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberratx)) Am ; batx ; xm ; rcondAtm ; norm1Atm ; LU1pm
  log=. log lcat ('getrsxlu1p'          tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxa )) An ; bxa  ; xn ; rcondAn  ; normiAn  ; LU1pn
  log=. log lcat ('getrsxlu1pc'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxac)) An ; bxac ; xn ; rcondAcn ; normiAcn ; LU1pn
  log=. log lcat ('getrsxlu1pt'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxat)) An ; bxat ; xn ; rcondAtn ; normiAtn ; LU1pn

  log=. log lcat ('getrspl1ux'          tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrax )) Am ; bax  ; xm ; rcondAm  ; norm1Am  ; pL1Um
  log=. log lcat ('getrspl1ucx'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberracx)) Am ; bacx ; xm ; rcondAcm ; norm1Acm ; pL1Um
  log=. log lcat ('getrspl1utx'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberratx)) Am ; batx ; xm ; rcondAtm ; norm1Atm ; pL1Um
  log=. log lcat ('getrsxpl1u'          tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxa )) An ; bxa  ; xn ; rcondAn  ; normiAn  ; pL1Un
  log=. log lcat ('getrsxpl1uc'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxac)) An ; bxac ; xn ; rcondAcn ; normiAcn ; pL1Un
  log=. log lcat ('getrsxpl1ut'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxat)) An ; bxat ; xn ; rcondAtn ; normiAtn ; pL1Un

  log=. log lcat ('getrspu1lx'          tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrax )) Am ; bax  ; xm ; rcondAm  ; norm1Am  ; pU1Lm
  log=. log lcat ('getrspu1lcx'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberracx)) Am ; bacx ; xm ; rcondAcm ; norm1Acm ; pU1Lm
  log=. log lcat ('getrspu1ltx'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberratx)) Am ; batx ; xm ; rcondAtm ; norm1Atm ; pU1Lm
  log=. log lcat ('getrsxpu1l'          tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxa )) An ; bxa  ; xn ; rcondAn  ; normiAn  ; pU1Ln
  log=. log lcat ('getrsxpu1lc'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxac)) An ; bxac ; xn ; rcondAcn ; normiAcn ; pU1Ln
  log=. log lcat ('getrsxpu1lt'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxat)) An ; bxat ; xn ; rcondAtn ; normiAtn ; pU1Ln

  log=. log lcat ('getrsul1px'          tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrax )) Am ; bax  ; xm ; rcondAm  ; norm1Am  ; UL1pm
  log=. log lcat ('getrsul1pcx'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberracx)) Am ; bacx ; xm ; rcondAcm ; norm1Acm ; UL1pm
  log=. log lcat ('getrsul1ptx'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberratx)) Am ; batx ; xm ; rcondAtm ; norm1Atm ; UL1pm
  log=. log lcat ('getrsxul1p'          tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxa )) An ; bxa  ; xn ; rcondAn  ; normiAn  ; UL1pn
  log=. log lcat ('getrsxul1pc'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxac)) An ; bxac ; xn ; rcondAcn ; normiAcn ; UL1pn
  log=. log lcat ('getrsxul1pt'         tdyad ((_2&{.)`(1&{::                        )`]`(3&{::)`t04v  `vberrxat)) An ; bxat ; xn ; rcondAtn ; normiAtn ; UL1pn

  coerase < 'mttmp'
  erase 'vferrr vferrl vberrax vberracx vberratx vberrxa vberrxac vberrxat'

  log
)

NB. ---------------------------------------------------------
NB. testhetrs
NB.
NB. Description:
NB.   Test:
NB.   - DSYTRS DSYTRS_AA ZHETRS ZHETRS_AA (math/lapack2
NB.     addon)
NB.   - hetrsxxxx (math/mt addon)
NB.   by Hermitian (symmetric) matrix
NB.
NB. Syntax:
NB.   log=. testhetrs (X ; A)
NB. where
NB.   X   - m×n-matrix, exact solutions
NB.   A   - k×k-matrix, the Hermitian (symmetric)
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - models LAPACK's DCHKSY_AA, ZCHKHE_AA with the following
NB.   difference:
NB.   - ferr is computed, too by call to t04x

testhetrs=: 3 : 0
  load_mttmp_ 'math/mt/test/lapack2/dsytrf'
  load_mttmp_ 'math/mt/test/lapack2/dsytrf_aa'
  load_mttmp_ 'math/mt/test/lapack2/zhetrf'
  load_mttmp_ 'math/mt/test/lapack2/zhetrf_aa'
  load_mttmp_ 'math/mt/test/lapack2/dsytrs'
  load_mttmp_ 'math/mt/test/lapack2/dsytrs_aa'
  load_mttmp_ 'math/mt/test/lapack2/zhetrs'
  load_mttmp_ 'math/mt/test/lapack2/zhetrs_aa'

  'X A'=. y
  'm n'=. $ X
  Am=. (2 # m) {. A
  An=. (2 # n) {. A

  rcondm=. hecon1 Am
  rcondn=. heconi An

  'normiAtm norm1Atm'=. 'norm1Am normiAm'=. (norm1 , normi) Am
  'normiAtn norm1Atn'=. 'norm1An normiAn'=. (norm1 , normi) An

  'pL1Tm pU1Tm'=. (hetrfpl ; <@hetrfpu) Am
  'pL1Tn pU1Tn'=. (hetrfpl ; <@hetrfpu) An

  NB. matrix X

  Bax=.      Am  mp X
  Batx=. (|: Am) mp X

  Bxa=.  X mp    An
  Bxat=. X mp |: An

  NB. vferrx
  vferrr=: normitc t04m  NB. for X at right side
  vferrl=: normitr t04m  NB. for X at left side

  NB. vberrxx for X at right side
  vberrax=:   mp~     t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc

  NB. vberrxx for X at left side
  vberrxa=:   mp      t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  log=.          ('''l''&dsytrs_mttmp_'    tmonad (        (1&{ ,~ 'l' dsytrf_mttmp_    0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am
  log=. log lcat ('''u''&dsytrs_mttmp_'    tmonad (        (1&{ ,~ 'u' dsytrf_mttmp_    0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am
  log=. log lcat ('''l''&zhetrs_mttmp_'    tmonad (        (1&{ ,~ 'l' zhetrf_mttmp_    0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am
  log=. log lcat ('''u''&zhetrs_mttmp_'    tmonad (        (1&{ ,~ 'u' zhetrf_mttmp_    0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am

  log=. log lcat ('''l''&dsytrs_aa_mttmp_' tmonad (        (1&{ ,~ 'l' dsytrf_aa_mttmp_ 0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am
  log=. log lcat ('''u''&dsytrs_aa_mttmp_' tmonad (        (1&{ ,~ 'u' dsytrf_aa_mttmp_ 0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am
  log=. log lcat ('''l''&zhetrs_aa_mttmp_' tmonad (        (1&{ ,~ 'l' zhetrf_aa_mttmp_ 0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am
  log=. log lcat ('''u''&zhetrs_aa_mttmp_' tmonad (        (1&{ ,~ 'u' zhetrf_aa_mttmp_ 0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am

  log=. log lcat ('hetrsplx'               tdyad  ((_3&{.)`(1&{::                            )`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am  ; pL1Tm
  log=. log lcat ('hetrspltx'              tdyad  ((_3&{.)`(1&{::                            )`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondm ; norm1Atm ; pL1Tm
  log=. log lcat ('hetrsxpl'               tdyad  ((_3&{.)`(1&{::                            )`]`(3&{::)`vferrl`vberrxa )) An ; Bxa  ; X  ; rcondn ; normiAn  ; pL1Tn
  log=. log lcat ('hetrsxplt'              tdyad  ((_3&{.)`(1&{::                            )`]`(3&{::)`vferrl`vberrxat)) An ; Bxat ; X  ; rcondn ; normiAtn ; pL1Tn

  log=. log lcat ('hetrspux'               tdyad  ((_3&{.)`(1&{::                            )`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am  ; pU1Tm
  log=. log lcat ('hetrsputx'              tdyad  ((_3&{.)`(1&{::                            )`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondm ; norm1Atm ; pU1Tm
  log=. log lcat ('hetrsxpu'               tdyad  ((_3&{.)`(1&{::                            )`]`(3&{::)`vferrl`vberrxa )) An ; Bxa  ; X  ; rcondn ; normiAn  ; pU1Tn
  log=. log lcat ('hetrsxput'              tdyad  ((_3&{.)`(1&{::                            )`]`(3&{::)`vferrl`vberrxat)) An ; Bxat ; X  ; rcondn ; normiAtn ; pU1Tn

  NB. vector x

  'xm xn'=. ({."1 ; {.) X

  bax=.      Am  mp xm
  batx=. (|: Am) mp xm

  bxa=.  xn mp    An
  bxat=. xn mp |: An

  NB. vberrxx for x at right side
  vberrax=:   mp~     t02v
  vberratx=: (mp~ |:) t02v

  NB. vberrxx for x at left side
  vberrxa=:   mp      t02v
  vberrxat=: (mp  |:) t02v

  log=. log lcat ('hetrsplx'               tdyad ((_3&{.)`(1&{::                             )`]`(3&{::)`t04v  `vberrax )) Am ; bax  ; xm ; rcondm ; norm1Am  ; pL1Tm
  log=. log lcat ('hetrspltx'              tdyad ((_3&{.)`(1&{::                             )`]`(3&{::)`t04v  `vberratx)) Am ; batx ; xm ; rcondm ; norm1Atm ; pL1Tm
  log=. log lcat ('hetrsxpl'               tdyad ((_3&{.)`(1&{::                             )`]`(3&{::)`t04v  `vberrxa )) An ; bxa  ; xn ; rcondn ; normiAn  ; pL1Tn
  log=. log lcat ('hetrsxplt'              tdyad ((_3&{.)`(1&{::                             )`]`(3&{::)`t04v  `vberrxat)) An ; bxat ; xn ; rcondn ; normiAtn ; pL1Tn

  log=. log lcat ('hetrspux'               tdyad ((_3&{.)`(1&{::                             )`]`(3&{::)`t04v  `vberrax )) Am ; bax  ; xm ; rcondm ; norm1Am  ; pU1Tm
  log=. log lcat ('hetrsputx'              tdyad ((_3&{.)`(1&{::                             )`]`(3&{::)`t04v  `vberratx)) Am ; batx ; xm ; rcondm ; norm1Atm ; pU1Tm
  log=. log lcat ('hetrsxpu'               tdyad ((_3&{.)`(1&{::                             )`]`(3&{::)`t04v  `vberrxa )) An ; bxa  ; xn ; rcondn ; normiAn  ; pU1Tn
  log=. log lcat ('hetrsxput'              tdyad ((_3&{.)`(1&{::                             )`]`(3&{::)`t04v  `vberrxat)) An ; bxat ; xn ; rcondn ; normiAtn ; pU1Tn

  coerase < 'mttmp'
  erase 'vferrr vferrl vberrax vberratx vberrxa vberrxat'

  log
)

NB. ---------------------------------------------------------
NB. testpotrs
NB.
NB. Description:
NB.   Test:
NB.   - xPOTRS (math/lapack2 addon)
NB.   - potrsxxx (math/mt addon)
NB.   by Hermitian (symmetric) positive definite matrix
NB.
NB. Syntax:
NB.   log=. testpotrs (X ; A)
NB. where
NB.   X   - m×n-matrix, exact solutions
NB.   A   - k×k-matrix, the Hermitian (symmetric) positive
NB.         definite
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - models LAPACK's xCHKPO

testpotrs=: 3 : 0
  load_mttmp_ 'math/mt/test/lapack2/potrf'
  load_mttmp_ 'math/mt/test/lapack2/potrs'

  'X A'=. y
  'm n'=. $ X
  Am=. (2 # m) {. A
  An=. (2 # n) {. A

  rcondm=. pocon1 Am
  rcondn=. poconi An

  'normiAtm norm1Atm'=. 'norm1Am normiAm'=. (norm1 , normi) Am
  'normiAtn norm1Atn'=. 'norm1An normiAn'=. (norm1 , normi) An

  'Lm Um'=. (potrfl ,: potrfu) Am
  'Ln Un'=. (potrfl ,: potrfu) An

  NB. matrix X

  Bax=.      Am  mp X
  Batx=. (|: Am) mp X

  Bxa=.  X mp    An
  Bxat=. X mp |: An

  NB. vferrx
  vferrr=: normitc t04m  NB. for X at right side
  vferrl=: normitr t04m  NB. for X at left side

  NB. vberrxx for X at right side
  vberrax=:   mp~     t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc

  NB. vberrxx for X at left side
  vberrxa=:   mp      t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  log=.          ('''l''&dpotrs_mttmp_' tmonad (        (1&{ ;~ 'l' dpotrf_mttmp_ 0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am
  log=. log lcat ('''u''&dpotrs_mttmp_' tmonad (        (1&{ ;~ 'u' dpotrf_mttmp_ 0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am
  log=. log lcat ('''l''&zpotrs_mttmp_' tmonad (        (1&{ ;~ 'l' zpotrf_mttmp_ 0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am
  log=. log lcat ('''u''&zpotrs_mttmp_' tmonad (        (1&{ ;~ 'u' zpotrf_mttmp_ 0&{::)`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am

  log=. log lcat ('potrslx'             tdyad  ((5&{::)`(1&{::                         )`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am  ; Lm
  log=. log lcat ('potrsltx'            tdyad  ((5&{::)`(1&{::                         )`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondm ; norm1Atm ; Lm
  log=. log lcat ('potrsxl'             tdyad  ((5&{::)`(1&{::                         )`]`(3&{::)`vferrl`vberrxa )) An ; Bxa  ; X  ; rcondn ; normiAn  ; Ln
  log=. log lcat ('potrsxlt'            tdyad  ((5&{::)`(1&{::                         )`]`(3&{::)`vferrl`vberrxat)) An ; Bxat ; X  ; rcondn ; normiAtn ; Ln

  log=. log lcat ('potrsux'             tdyad  ((5&{::)`(1&{::                         )`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am  ; Um
  log=. log lcat ('potrsutx'            tdyad  ((5&{::)`(1&{::                         )`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondm ; norm1Atm ; Um
  log=. log lcat ('potrsxu'             tdyad  ((5&{::)`(1&{::                         )`]`(3&{::)`vferrl`vberrxa )) An ; Bxa  ; X  ; rcondn ; normiAn  ; Un
  log=. log lcat ('potrsxut'            tdyad  ((5&{::)`(1&{::                         )`]`(3&{::)`vferrl`vberrxat)) An ; Bxat ; X  ; rcondn ; normiAtn ; Un

  NB. vector x

  'xm xn'=. ({."1 ; {.) X

  bax=.      Am  mp xm
  batx=. (|: Am) mp xm

  bxa=.  xn mp    An
  bxat=. xn mp |: An

  NB. vberrxx for x at right side
  vberrax=:   mp~     t02v
  vberratx=: (mp~ |:) t02v

  NB. vberrxx for x at left side
  vberrxa=:   mp      t02v
  vberrxat=: (mp  |:) t02v

  log=. log lcat ('potrslx'             tdyad ((5&{::)`(1&{::                          )`]`(3&{::)`t04v  `vberrax )) Am ; bax  ; xm ; rcondm ; norm1Am  ; Lm
  log=. log lcat ('potrsltx'            tdyad ((5&{::)`(1&{::                          )`]`(3&{::)`t04v  `vberratx)) Am ; batx ; xm ; rcondm ; norm1Atm ; Lm
  log=. log lcat ('potrsxl'             tdyad ((5&{::)`(1&{::                          )`]`(3&{::)`t04v  `vberrxa )) An ; bxa  ; xn ; rcondn ; normiAn  ; Ln
  log=. log lcat ('potrsxlt'            tdyad ((5&{::)`(1&{::                          )`]`(3&{::)`t04v  `vberrxat)) An ; bxat ; xn ; rcondn ; normiAtn ; Ln

  log=. log lcat ('potrsux'             tdyad ((5&{::)`(1&{::                          )`]`(3&{::)`t04v  `vberrax )) Am ; bax  ; xm ; rcondm ; norm1Am  ; Um
  log=. log lcat ('potrsutx'            tdyad ((5&{::)`(1&{::                          )`]`(3&{::)`t04v  `vberratx)) Am ; batx ; xm ; rcondm ; norm1Atm ; Um
  log=. log lcat ('potrsxu'             tdyad ((5&{::)`(1&{::                          )`]`(3&{::)`t04v  `vberrxa )) An ; bxa  ; xn ; rcondn ; normiAn  ; Un
  log=. log lcat ('potrsxut'            tdyad ((5&{::)`(1&{::                          )`]`(3&{::)`t04v  `vberrxat)) An ; bxat ; xn ; rcondn ; normiAtn ; Un

  coerase < 'mttmp'
  erase 'vferrr vferrl vberrax vberratx vberrxa vberrxat'

  log
)

NB. ---------------------------------------------------------
NB. testpttrs
NB.
NB. Description:
NB.   Test:
NB.   - xPTTRS (math/lapack2 addon)
NB.   - pttrsxxx (math/mt addon)
NB.   by Hermitian (symmetric) positive definite tridiagonal
NB.   matrix
NB.
NB. Syntax:
NB.   log=. testpttrs (X ; A)
NB. where
NB.   X   - m×n-matrix, exact solutions
NB.   A   - k×k-matrix, the Hermitian (symmetric) positive
NB.         definite tridiagonal
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - models LAPACK's xCHKPT
NB.
NB. TODO:
NB. - A would be sparse

testpttrs=: 3 : 0
  load_mttmp_ 'math/mt/test/lapack2/pttrf'
  load_mttmp_ 'math/mt/test/lapack2/pttrs'

  'X A'=. y
  'm n'=. $ X
  Am=. (2 # m) {. A
  An=. (2 # n) {. A

  rcondm=. ptcon1 Am
  rcondn=. ptconi An

  'normiAtm norm1Atm'=. 'norm1Am normiAm'=. (norm1 , normi) Am
  'normiAtn norm1Atn'=. 'norm1An normiAn'=. (norm1 , normi) An

  'L1Dm U1Dm'=. (pttrfl ; <@pttrfu) Am
  'L1Dn U1Dn'=. (pttrfl ; <@pttrfu) An

  NB. matrix X

  Bax=.      Am  mp X
  Batx=. (|: Am) mp X

  Bxa=.  X mp    An
  Bxat=. X mp |: An

  NB. vferrx
  vferrr=: normitc t04m  NB. for X at right side
  vferrl=: normitr t04m  NB. for X at left side

  NB. vberrxx for X at right side
  vberrax=:   mp~     t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc

  NB. vberrxx for X at left side
  vberrxa=:   mp      t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  log=.          ('dpttrs_mttmp_'       tmonad (        (1&{ ,~ dpttrf_mttmp_@(diag ; _1&diag)@(0&{::))`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am
  log=. log lcat ('''l''&zpttrs_mttmp_' tmonad (        (1&{ ,~ zpttrf_mttmp_@(diag ; _1&diag)@(0&{::))`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am

  log=. log lcat ('pttrslx'             tdyad  ((_2&{.)`(1&{::                                        )`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am  ; L1Dm
  log=. log lcat ('pttrsltx'            tdyad  ((_2&{.)`(1&{::                                        )`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondm ; norm1Atm ; L1Dm
  log=. log lcat ('pttrsxl'             tdyad  ((_2&{.)`(1&{::                                        )`]`(3&{::)`vferrl`vberrxa )) An ; Bxa  ; X  ; rcondn ; normiAn  ; L1Dn
  log=. log lcat ('pttrsxlt'            tdyad  ((_2&{.)`(1&{::                                        )`]`(3&{::)`vferrl`vberrxat)) An ; Bxat ; X  ; rcondn ; normiAtn ; L1Dn

  log=. log lcat ('pttrsux'             tdyad  ((_2&{.)`(1&{::                                        )`]`(3&{::)`vferrr`vberrax )) Am ; Bax  ; X  ; rcondm ; norm1Am  ; U1Dm
  log=. log lcat ('pttrsutx'            tdyad  ((_2&{.)`(1&{::                                        )`]`(3&{::)`vferrr`vberratx)) Am ; Batx ; X  ; rcondm ; norm1Atm ; U1Dm
  log=. log lcat ('pttrsxu'             tdyad  ((_2&{.)`(1&{::                                        )`]`(3&{::)`vferrl`vberrxa )) An ; Bxa  ; X  ; rcondn ; normiAn  ; U1Dn
  log=. log lcat ('pttrsxut'            tdyad  ((_2&{.)`(1&{::                                        )`]`(3&{::)`vferrl`vberrxat)) An ; Bxat ; X  ; rcondn ; normiAtn ; U1Dn

  NB. vector x

  'xm xn'=. ({."1 ; {.) X

  bax=.      Am  mp xm
  batx=. (|: Am) mp xm

  bxa=.  xn mp    An
  bxat=. xn mp |: An

  NB. vberrxx for x at right side
  vberrax=:   mp~     t02v
  vberratx=: (mp~ |:) t02v

  NB. vberrxx for x at left side
  vberrxa=:   mp      t02v
  vberrxat=: (mp  |:) t02v

  log=. log lcat ('pttrslx'             tdyad ((_2&{.)`(1&{::                                         )`]`(3&{::)`t04v  `vberrax )) Am ; bax  ; xm ; rcondm ; norm1Am  ; L1Dm
  log=. log lcat ('pttrsltx'            tdyad ((_2&{.)`(1&{::                                         )`]`(3&{::)`t04v  `vberratx)) Am ; batx ; xm ; rcondm ; norm1Atm ; L1Dm
  log=. log lcat ('pttrsxl'             tdyad ((_2&{.)`(1&{::                                         )`]`(3&{::)`t04v  `vberrxa )) An ; bxa  ; xn ; rcondn ; normiAn  ; L1Dn
  log=. log lcat ('pttrsxlt'            tdyad ((_2&{.)`(1&{::                                         )`]`(3&{::)`t04v  `vberrxat)) An ; bxat ; xn ; rcondn ; normiAtn ; L1Dn

  log=. log lcat ('pttrsux'             tdyad ((_2&{.)`(1&{::                                         )`]`(3&{::)`t04v  `vberrax )) Am ; bax  ; xm ; rcondm ; norm1Am  ; U1Dm
  log=. log lcat ('pttrsutx'            tdyad ((_2&{.)`(1&{::                                         )`]`(3&{::)`t04v  `vberratx)) Am ; batx ; xm ; rcondm ; norm1Atm ; U1Dm
  log=. log lcat ('pttrsxu'             tdyad ((_2&{.)`(1&{::                                         )`]`(3&{::)`t04v  `vberrxa )) An ; bxa  ; xn ; rcondn ; normiAn  ; U1Dn
  log=. log lcat ('pttrsxut'            tdyad ((_2&{.)`(1&{::                                         )`]`(3&{::)`t04v  `vberrxat)) An ; bxat ; xn ; rcondn ; normiAtn ; U1Dn

  coerase < 'mttmp'
  erase 'vferrr vferrl vberrax vberratx vberrxa vberrxat'

  log
)

NB. ---------------------------------------------------------
NB. testtrtrs
NB.
NB. Description:
NB.   Test:
NB.   - xTRTRS (math/lapack2 addon)
NB.   - trtrsxxx (math/mt addon)
NB.   by triangular matrix
NB.
NB. Syntax:
NB.   log=. testtrtrs (X ; A)
NB. where
NB.   X   - m×n-matrix, exact solutions
NB.   A   - m×m-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.
NB. Notes:
NB. - models LAPACK's xCHKTR

testtrtrs=: 3 : 0
  load_mttmp_ 'math/mt/test/lapack2/trtrs'

  'X A'=. y

  rcondL=.  trlcon1  L=.           trlpick A
  rcondU=.  trucon1  U=.           trupick A
  rcondL1=. trl1con1 L1=. (1 ; '') setdiag L
  rcondU1=. tru1con1 U1=. (1 ; '') setdiag U

  'norm1L  normiL '=. (norm1 , normi) L
  'norm1L1 normiL1'=. (norm1 , normi) L1
  'norm1U  normiU '=. (norm1 , normi) U
  'norm1U1 normiU1'=. (norm1 , normi) U1

  NB. matrix X

  Blnn=.     L   mp X
  Blnu=.     L1  mp X
  Blcn=. (ct L ) mp X
  Blcu=. (ct L1) mp X
  Bltn=. (|: L ) mp X
  Bltu=. (|: L1) mp X

  Bunn=.     U   mp X
  Bunu=.     U1  mp X
  Bucn=. (ct U ) mp X
  Bucu=. (ct U1) mp X
  Butn=. (|: U ) mp X
  Butu=. (|: U1) mp X

  vferrr=: normitc t04m

  vberrlnn=: (mp~    trlpick ) t02m norm1tc
  vberrlnu=: (mp~    trl1pick) t02m norm1tc
  vberrlcn=: (mp~ ct@trlpick ) t02m norm1tc
  vberrlcu=: (mp~ ct@trl1pick) t02m norm1tc
  vberrltn=: (mp~ |:@trlpick ) t02m norm1tc
  vberrltu=: (mp~ |:@trl1pick) t02m norm1tc
  vberrunn=: (mp~    trupick ) t02m norm1tc
  vberrunu=: (mp~    tru1pick) t02m norm1tc
  vberrucn=: (mp~ ct@trupick ) t02m norm1tc
  vberrucu=: (mp~ ct@tru1pick) t02m norm1tc
  vberrutn=: (mp~ |:@trupick ) t02m norm1tc
  vberrutu=: (mp~ |:@tru1pick) t02m norm1tc

  log=.          ('''lnn''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrlnn                )) A ; Blnn               ; X  ; rcondL  ; norm1L
  log=. log lcat ('''lnu''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrlnu                )) A ; Blnu               ; X  ; rcondL1 ; norm1L1
  log=. log lcat ('''ltn''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrltn                )) A ; Bltn               ; X  ; rcondL  ; normiL
  log=. log lcat ('''ltu''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrltu                )) A ; Bltu               ; X  ; rcondL1 ; normiL1
  log=. log lcat ('''lcn''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrlcn                )) A ; Blcn               ; X  ; rcondL  ; normiL
  log=. log lcat ('''lcu''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrlcu                )) A ; Blcu               ; X  ; rcondL1 ; normiL1

  log=. log lcat ('''unn''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrunn                )) A ; Bunn               ; X  ; rcondL  ; norm1U
  log=. log lcat ('''unu''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrunu                )) A ; Bunu               ; X  ; rcondL1 ; norm1U1
  log=. log lcat ('''utn''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrutn                )) A ; Butn               ; X  ; rcondL  ; normiU
  log=. log lcat ('''utu''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrutu                )) A ; Butu               ; X  ; rcondL1 ; normiU1
  log=. log lcat ('''ucn''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrucn                )) A ; Bucn               ; X  ; rcondL  ; normiU
  log=. log lcat ('''ucu''&dtrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrucu                )) A ; Bucu               ; X  ; rcondL1 ; normiU1

  log=. log lcat ('''lnn''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrlnn                )) A ; Blnn               ; X  ; rcondL  ; norm1L
  log=. log lcat ('''lnu''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrlnu                )) A ; Blnu               ; X  ; rcondL1 ; norm1L1
  log=. log lcat ('''ltn''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrltn                )) A ; Bltn               ; X  ; rcondL  ; normiL
  log=. log lcat ('''ltu''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrltu                )) A ; Bltu               ; X  ; rcondL1 ; normiL1
  log=. log lcat ('''lcn''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrlcn                )) A ; Blcn               ; X  ; rcondL  ; normiL
  log=. log lcat ('''lcu''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrlcu                )) A ; Blcu               ; X  ; rcondL1 ; normiL1

  log=. log lcat ('''unn''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrunn                )) A ; Bunn               ; X  ; rcondL  ; norm1U
  log=. log lcat ('''unu''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrunu                )) A ; Bunu               ; X  ; rcondL1 ; norm1U1
  log=. log lcat ('''utn''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrutn                )) A ; Butn               ; X  ; rcondL  ; normiU
  log=. log lcat ('''utu''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrutu                )) A ; Butu               ; X  ; rcondL1 ; normiU1
  log=. log lcat ('''ucn''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrucn                )) A ; Bucn               ; X  ; rcondL  ; normiU
  log=. log lcat ('''ucu''&ztrtrs_mttmp_' tmonad (        (2&{. )`]`(3&{::)`vferrr`vberrucu                )) A ; Bucu               ; X  ; rcondL1 ; normiU1

  log=. log lcat ('trtrslnn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrlnn                )) A ; Blnn               ; X  ; rcondL  ; norm1L
  log=. log lcat ('trtrslnu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrlnu                )) A ; Blnu               ; X  ; rcondL1 ; norm1L1
  log=. log lcat ('trtrsltn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrltn                )) A ; Bltn               ; X  ; rcondL  ; normiL
  log=. log lcat ('trtrsltu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrltu                )) A ; Bltu               ; X  ; rcondL1 ; normiL1
  log=. log lcat ('trtrslcn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrlcn                )) A ; Blcn               ; X  ; rcondL  ; normiL
  log=. log lcat ('trtrslcu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrlcu                )) A ; Blcu               ; X  ; rcondL1 ; normiL1

  log=. log lcat ('trtrsunn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrunn                )) A ; Bunn               ; X  ; rcondU  ; norm1U
  log=. log lcat ('trtrsunu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrunu                )) A ; Bunu               ; X  ; rcondU1 ; norm1U1
  log=. log lcat ('trtrsutn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrutn                )) A ; Butn               ; X  ; rcondU  ; normiU
  log=. log lcat ('trtrsutu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrutu                )) A ; Butu               ; X  ; rcondU1 ; normiU1
  log=. log lcat ('trtrsucn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrucn                )) A ; Bucn               ; X  ; rcondU  ; normiU
  log=. log lcat ('trtrsucu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrr`vberrucu                )) A ; Bucu               ; X  ; rcondU1 ; normiU1

  NB. vector x

  xm=. {."1 X

  log=. log lcat ('trtrslnn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~    trlpick ) t02v))) A ; (L   mp       xm ) ; xm ; rcondL  ; norm1L
  log=. log lcat ('trtrslnu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~    trl1pick) t02v))) A ; (L1  mp       xm ) ; xm ; rcondL1 ; norm1L1
  log=. log lcat ('trtrslcn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~ ct@trlpick ) t02v))) A ; (L  (mp~ ct)~ xm ) ; xm ; rcondL  ; normiL
  log=. log lcat ('trtrslcu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~ ct@trl1pick) t02v))) A ; (L1 (mp~ ct)~ xm ) ; xm ; rcondL1 ; normiL1
  log=. log lcat ('trtrsltn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~ |:@trlpick ) t02v))) A ; (L  (mp~ |:)~ xm ) ; xm ; rcondL  ; normiL
  log=. log lcat ('trtrsltu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~ |:@trl1pick) t02v))) A ; (L1 (mp~ |:)~ xm ) ; xm ; rcondL1 ; normiL1

  log=. log lcat ('trtrsunn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~    trupick ) t02v))) A ; (U   mp       xm ) ; xm ; rcondU  ; norm1U
  log=. log lcat ('trtrsunu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~    tru1pick) t02v))) A ; (U1  mp       xm ) ; xm ; rcondU1 ; norm1U1
  log=. log lcat ('trtrsucn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~ ct@trupick ) t02v))) A ; (U  (mp~ ct)~ xm ) ; xm ; rcondU  ; normiU
  log=. log lcat ('trtrsucu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~ ct@tru1pick) t02v))) A ; (U1 (mp~ ct)~ xm ) ; xm ; rcondU1 ; normiU1
  log=. log lcat ('trtrsutn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~ |:@trupick ) t02v))) A ; (U  (mp~ |:)~ xm ) ; xm ; rcondU  ; normiU
  log=. log lcat ('trtrsutu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`t04v  `((mp~ |:@tru1pick) t02v))) A ; (U1 (mp~ |:)~ xm ) ; xm ; rcondU1 ; normiU1

  coerase < 'mttmp'
  erase 'vferrr vberrlnn vberrlnu vberrltn vberrltu vberrlcn vberrlcu vberrunn vberrunu vberrutn vberrutu vberrucn vberrucu vberrxnx vberrxcx vberrxtx'

  log
)

NB. ---------------------------------------------------------
NB. testtrs
NB.
NB. Description:
NB.   Adv. to make verb to test xxtrsxxxxxx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testtrs) (m,n)
NB. where
NB.   mkmat - monad to generate a material for matrix; is
NB.           called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random real matrices with elements distributed
NB.   uniformly with support (0,1), system is either
NB.   100×100-matrix or 150×150-matrix and RHS is either
NB.   100-vector or 150-vector or 100×150-matrix:
NB.     log=. ?@$&0 testtrs_mt_ 100 150
NB. - test by random real matrices with elements with
NB.   limited value's amplitude, system is 150×150-matrix and
NB.   RHS is either 150-vector or 150×150-matrix:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testtrs_mt_ 150 150
NB. - test by random complex matrices, system is either
NB.   200×200-matrix or 150×150-matrix and RHS is either
NB.   200-vector or 150-vector or 200×150-matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testtrs_mt_ 200 150

testtrs=: 1 : 'testtrtrs_mt_@(u ; u@(2 # {.)) ,&.>~ (testpttrs_mt_@(u@{. ; (u ptmat2_mt_)@{:) ,&.>~ testpotrs_mt_@(u@{. ; (u pomat_mt_)@{:) ,&.>~ testhetrs_mt_@(u@{. ; (u hemat_mt_)@{:) ,&.>~ testgetrs_mt_@:(<@u"1))@(,: >./)'
