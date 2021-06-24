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
NB.
NB. testgetrs1   Test getrsxxxxxx by general square matrix
NB.              and single RHS
NB. testgetrs3   Test getrsxxxxxx by general square matrix
NB.              and multiple RHS
NB. testhetrs1   Test hetrsxxxx by Hermitian (symmetric)
NB.              matrix and single RHS
NB. testhetrs3   Test hetrsxxxx by Hermitian (symmetric)
NB.              matrix and multiple RHS
NB. testpotrs1   Test potrsxxx by Hermitian (symmetric)
NB.              positive definite matrix and single RHS
NB. testpotrs3   Test potrsxxx by Hermitian (symmetric)
NB.              positive definite matrix and multiple RHS
NB. testpttrs1   Test pttrsxxx by Hermitian (symmetric)
NB.              positive definite tridiagonal matrix and
NB.              single RHS
NB. testpttrs3   Test pttrsxxx by Hermitian (symmetric)
NB.              positive definite tridiagonal matrix and
NB.              multiple RHS
NB. testtrs      Adv. to make verb to test xxtrsxxx by matrix
NB.              of generator and shape given
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
NB. Verb:          Solves:        Syntax:
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
NB.   Xv   - same shape as Bv, the solution[s]
NB.   Xh   - same shape as Bh, the solution[s]
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
NB. Verb:          Solves:        Syntax:
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
NB.   Xv   - same shape as Bv, the solution[s]
NB.   Xh   - same shape as Bh, the solution[s]
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
NB. - implements LAPACK's xGETRS

getrspl1ux=:  (1 {:: [) ([ trsmlunn trsmllnu) ((0 {:: [) C.   ])
getrsxpl1uc=: (1 {:: [) ([ trsmrucn trsmrlcu) ((0 {:: [) C."1 ])
getrsxpl1ut=: (1 {:: [) ([ trsmrutn trsmrltu) ((0 {:: [) C."1 ])

getrspl1ucx=: (0 {:: [) C.^:_1   ((] trsmllcu trsmlucn~) 1&{::)~
getrspl1utx=: (0 {:: [) C.^:_1   ((] trsmlltu trsmlutn~) 1&{::)~
getrsxpl1u=:  (0 {:: [) C.^:_1"1 ((] trsmrlnu trsmrunn~) 1&{::)~

NB. ---------------------------------------------------------
NB. Verb:          Solves:        Syntax:
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
NB.   Xv   - same shape as Bv, the solution[s]
NB.   Xh   - same shape as Bh, the solution[s]
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
NB. Verb:          Solves:        Syntax:
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
NB.   Xv   - same shape as Bv, the solution[s]
NB.   Xh   - same shape as Bh, the solution[s]
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
NB. Verb:        Solves:        Syntax:
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
NB.   Xv   - same shape as Bv, the solution[s]
NB.   Xh   - same shape as Bh, the solution[s]
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
NB. Verb:        Solves:        Syntax:
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
NB.   Xv   - same shape as Bv, the solution[s]
NB.   Xh   - same shape as Bh, the solution[s]
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
NB. Verb:       Solves:        Syntax:
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
NB.   Xv   - same shape as Bv, the solution[s]
NB.   Xh   - same shape as Bh, the solution[s]
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
NB. Verb:       Solves:        Syntax:
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
NB.   Xv   - same shape as Bv, the solution[s]
NB.   Xh   - same shape as Bh, the solution[s]
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
NB. Verb:       Solves:        Syntax:
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
NB.   Xv   - same shape as Bv, the solution[s]
NB.   Xh   - same shape as Bh, the solution[s]
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
NB. Verb:       Solves:        Syntax:
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
NB.   Xv   - same shape as Bv, the solution[s]
NB.   Xh   - same shape as Bh, the solution[s]
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

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgetrs1
NB.
NB. Description:
NB.   Test getrsxxx by general square matrix and single RHS
NB.
NB. Syntax:
NB.   testgetrs1 (A ; x)
NB. where
NB.   A - n×n-matrix
NB.   x - n-vector, the exact solution

testgetrs1=: 3 : 0
  'A x'=. y

  'rcondA rcondAc rcondAt'=. (gecon1 , gecon1@ct , gecon1@|:) A

  'normiAt norm1At'=. 'normiAc norm1Ac'=. 'norm1A normiA'=. (norm1 , normi) A

  'LU1p pL1U pU1L UL1p'=. (getrflu1p ; getrfpl1u ; getrfpu1l ; <@getrful1p) A

  Bax=.      A  mp x
  Bacx=. (ct A) mp x
  Batx=. (|: A) mp x

  Bxa=.  x mp    A
  Bxac=. x mp ct A
  Bxat=. x mp |: A

  NB. vberrX for x at right side
  vberrax=:   mp~     t02v
  vberracx=: (mp~ ct) t02v
  vberratx=: (mp~ |:) t02v
  NB. vberrX for x at left side
  vberrxa=:   mp      t02v
  vberrxac=: (mp  ct) t02v
  vberrxat=: (mp  |:) t02v

  ('getrslu1px'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrax )) A ; Bax  ; x ; rcondA  ; norm1A  ; LU1p
  ('getrslu1pcx' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberracx)) A ; Bacx ; x ; rcondAc ; norm1Ac ; LU1p
  ('getrslu1ptx' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberratx)) A ; Batx ; x ; rcondAt ; norm1At ; LU1p
  ('getrsxlu1p'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxa )) A ; Bxa  ; x ; rcondA  ; normiA  ; LU1p
  ('getrsxlu1pc' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxac)) A ; Bxac ; x ; rcondAc ; normiAc ; LU1p
  ('getrsxlu1pt' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxat)) A ; Bxat ; x ; rcondAt ; normiAt ; LU1p

  ('getrspl1ux'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrax )) A ; Bax  ; x ; rcondA  ; norm1A  ; pL1U
  ('getrspl1ucx' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberracx)) A ; Bacx ; x ; rcondAc ; norm1Ac ; pL1U
  ('getrspl1utx' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberratx)) A ; Batx ; x ; rcondAt ; norm1At ; pL1U
  ('getrsxpl1u'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxa )) A ; Bxa  ; x ; rcondA  ; normiA  ; pL1U
  ('getrsxpl1uc' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxac)) A ; Bxac ; x ; rcondAc ; normiAc ; pL1U
  ('getrsxpl1ut' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxat)) A ; Bxat ; x ; rcondAt ; normiAt ; pL1U

  ('getrspu1lx'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrax )) A ; Bax  ; x ; rcondA  ; norm1A  ; pU1L
  ('getrspu1lcx' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberracx)) A ; Bacx ; x ; rcondAc ; norm1Ac ; pU1L
  ('getrspu1ltx' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberratx)) A ; Batx ; x ; rcondAt ; norm1At ; pU1L
  ('getrsxpu1l'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxa )) A ; Bxa  ; x ; rcondA  ; normiA  ; pU1L
  ('getrsxpu1lc' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxac)) A ; Bxac ; x ; rcondAc ; normiAc ; pU1L
  ('getrsxpu1lt' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxat)) A ; Bxat ; x ; rcondAt ; normiAt ; pU1L

  ('getrsul1px'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrax )) A ; Bax  ; x ; rcondA  ; norm1A  ; UL1p
  ('getrsul1pcx' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberracx)) A ; Bacx ; x ; rcondAc ; norm1Ac ; UL1p
  ('getrsul1ptx' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberratx)) A ; Batx ; x ; rcondAt ; norm1At ; UL1p
  ('getrsxul1p'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxa )) A ; Bxa  ; x ; rcondA  ; normiA  ; UL1p
  ('getrsxul1pc' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxac)) A ; Bxac ; x ; rcondAc ; normiAc ; UL1p
  ('getrsxul1pt' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxat)) A ; Bxat ; x ; rcondAt ; normiAt ; UL1p

  coerase < 'mttmp'
  erase 'vberrax vberracx vberratx vberrxa vberrxac vberrxat'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testgetrs3
NB.
NB. Description:
NB.   Test:
NB.   - xGETRS (math/lapack2 addon)
NB.   - getrsxxx (math/mt addon)
NB.   by general square matrix and multiple RHS
NB.
NB. Syntax:
NB.   testgetrs3 (A ; X)
NB. where
NB.   A - n×n-matrix
NB.   X - n×3-matrix, exact solutions

testgetrs3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/getrf'
  load_mttmp_ :: ] 'math/mt/test/lapack2/getrs'

  'A Xv'=. y
  Xh=. |: Xv

  'rcondA rcondAc rcondAt'=. (gecon1 , gecon1@ct , gecon1@|:) A

  'normiAt norm1At'=. 'normiAc norm1Ac'=. 'norm1A normiA'=. (norm1 , normi) A

  'LU1p pL1U pU1L UL1p'=. (getrflu1p ; getrfpl1u ; getrfpu1l ; <@getrful1p) A

  Bax=.      A  mp Xv
  Bacx=. (ct A) mp Xv
  Batx=. (|: A) mp Xv

  Bxa=.  Xh mp    A
  Bxac=. Xh mp ct A
  Bxat=. Xh mp |: A

  vferrv=: normitc t04m  NB. for Xv
  vferrh=: normitr t04m  NB. for Xh

  NB. vberrX for Xv at right side
  vberrax=:   mp~     t02m norm1tc
  vberracx=: (mp~ ct) t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc
  NB. vberrX for Xh at left side
  vberrxa=:   mp      t02m norm1tr
  vberrxac=: (mp  ct) t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  ('''n''&dgetrs_mttmp_' tmonad ((1&{ ,~ dgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcondA  ; norm1A
  ('''c''&dgetrs_mttmp_' tmonad ((1&{ ,~ dgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrv`vberracx)) A ; Bacx ; Xv ; rcondAc ; norm1Ac
  ('''t''&dgetrs_mttmp_' tmonad ((1&{ ,~ dgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcondAt ; norm1At

  ('''n''&zgetrs_mttmp_' tmonad ((1&{ ,~ zgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcondA  ; norm1A
  ('''c''&zgetrs_mttmp_' tmonad ((1&{ ,~ zgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrv`vberracx)) A ; Bacx ; Xv ; rcondAc ; norm1Ac
  ('''t''&zgetrs_mttmp_' tmonad ((1&{ ,~ zgetrf_mttmp_@(0&{::))`]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcondAt ; norm1At

  ('getrslu1px'          tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcondA  ; norm1A  ; LU1p
  ('getrslu1pcx'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberracx)) A ; Bacx ; Xv ; rcondAc ; norm1Ac ; LU1p
  ('getrslu1ptx'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcondAt ; norm1At ; LU1p
  ('getrsxlu1p'          tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxa )) A ; Bxa  ; Xh ; rcondA  ; normiA  ; LU1p
  ('getrsxlu1pc'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxac)) A ; Bxac ; Xh ; rcondAc ; normiAc ; LU1p
  ('getrsxlu1pt'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxat)) A ; Bxat ; Xh ; rcondAt ; normiAt ; LU1p

  ('getrspl1ux'          tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcondA  ; norm1A  ; pL1U
  ('getrspl1ucx'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberracx)) A ; Bacx ; Xv ; rcondAc ; norm1Ac ; pL1U
  ('getrspl1utx'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcondAt ; norm1At ; pL1U
  ('getrsxpl1u'          tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxa )) A ; Bxa  ; Xh ; rcondA  ; normiA  ; pL1U
  ('getrsxpl1uc'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxac)) A ; Bxac ; Xh ; rcondAc ; normiAc ; pL1U
  ('getrsxpl1ut'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxat)) A ; Bxat ; Xh ; rcondAt ; normiAt ; pL1U

  ('getrspu1lx'          tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcondA  ; norm1A  ; pU1L
  ('getrspu1lcx'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberracx)) A ; Bacx ; Xv ; rcondAc ; norm1Ac ; pU1L
  ('getrspu1ltx'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcondAt ; norm1At ; pU1L
  ('getrsxpu1l'          tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxa )) A ; Bxa  ; Xh ; rcondA  ; normiA  ; pU1L
  ('getrsxpu1lc'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxac)) A ; Bxac ; Xh ; rcondAc ; normiAc ; pU1L
  ('getrsxpu1lt'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxat)) A ; Bxat ; Xh ; rcondAt ; normiAt ; pU1L

  ('getrsul1px'          tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcondA  ; norm1A  ; UL1p
  ('getrsul1pcx'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberracx)) A ; Bacx ; Xv ; rcondAc ; norm1Ac ; UL1p
  ('getrsul1ptx'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcondAt ; norm1At ; UL1p
  ('getrsxul1p'          tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxa )) A ; Bxa  ; Xh ; rcondA  ; normiA  ; UL1p
  ('getrsxul1pc'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxac)) A ; Bxac ; Xh ; rcondAc ; normiAc ; UL1p
  ('getrsxul1pt'         tdyad  ((_2&{.)`(1&{::)               `]`(3&{::)`vferrh`vberrxat)) A ; Bxat ; Xh ; rcondAt ; normiAt ; UL1p

  coerase < 'mttmp'
  erase 'vferrv vferrh vberrax vberracx vberratx vberrxa vberrxac vberrxat'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhetrs1
NB.
NB. Description:
NB.   Test hetrsxxxx by Hermitian (symmetric) matrix and
NB.   single RHS
NB.
NB. Syntax:
NB.   testhetrs1 (A ; x)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric)
NB.   x - n-vector, the exact solution
NB.
NB. Notes:
NB. - models LAPACK's DCHKSY_AA, ZCHKHE_AA with the following
NB.   difference:
NB.   - ferr is computed, too by call to t04v

testhetrs1=: 3 : 0
  'A x'=. y

  rcond=. hecon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  'pL1T pU1T'=. (hetrfpl ; <@hetrfpu) A

  Bax=.      A  mp x
  Batx=. (|: A) mp x
  Bxa=.  x mp    A
  Bxat=. x mp |: A

  NB. vberrX for x at right side
  vberrax=:   mp~     t02v
  vberratx=: (mp~ |:) t02v
  NB. vberrX for x at left side
  vberrxa=:   mp      t02v
  vberrxat=: (mp  |:) t02v

  ('hetrsplx'  tdyad ((_3&{.)`(1&{::)`]`(3&{::)`t04v`vberrax )) A ; Bax  ; x ; rcond ; norm1A  ; pL1T
  ('hetrspltx' tdyad ((_3&{.)`(1&{::)`]`(3&{::)`t04v`vberratx)) A ; Batx ; x ; rcond ; norm1At ; pL1T
  ('hetrsxpl'  tdyad ((_3&{.)`(1&{::)`]`(3&{::)`t04v`vberrxa )) A ; Bxa  ; x ; rcond ; normiA  ; pL1T
  ('hetrsxplt' tdyad ((_3&{.)`(1&{::)`]`(3&{::)`t04v`vberrxat)) A ; Bxat ; x ; rcond ; normiAt ; pL1T

  ('hetrspux'  tdyad ((_3&{.)`(1&{::)`]`(3&{::)`t04v`vberrax )) A ; Bax  ; x ; rcond ; norm1A  ; pU1T
  ('hetrsputx' tdyad ((_3&{.)`(1&{::)`]`(3&{::)`t04v`vberratx)) A ; Batx ; x ; rcond ; norm1At ; pU1T
  ('hetrsxpu'  tdyad ((_3&{.)`(1&{::)`]`(3&{::)`t04v`vberrxa )) A ; Bxa  ; x ; rcond ; normiA  ; pU1T
  ('hetrsxput' tdyad ((_3&{.)`(1&{::)`]`(3&{::)`t04v`vberrxat)) A ; Bxat ; x ; rcond ; normiAt ; pU1T

  coerase < 'mttmp'
  erase 'vberrax vberratx vberrxa vberrxat'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhetrs3
NB.
NB. Description:
NB.   Test:
NB.   - DSYTRS DSYTRS_AA ZHETRS ZHETRS_AA (math/lapack2
NB.     addon)
NB.   - hetrsxxxx (math/mt addon)
NB.   by Hermitian (symmetric) matrix and multiple RHS
NB.
NB. Syntax:
NB.   testhetrs3 (A ; X)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric)
NB.   X - n×3-matrix, exact solutions
NB.
NB. Notes:
NB. - models LAPACK's DCHKSY_AA, ZCHKHE_AA with the following
NB.   difference:
NB.   - ferr is computed, too by call to t04m

testhetrs3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/dsytrf'
  load_mttmp_ :: ] 'math/mt/test/lapack2/dsytrf_aa'
  load_mttmp_ :: ] 'math/mt/test/lapack2/zhetrf'
  load_mttmp_ :: ] 'math/mt/test/lapack2/zhetrf_aa'
  load_mttmp_ :: ] 'math/mt/test/lapack2/dsytrs'
  load_mttmp_ :: ] 'math/mt/test/lapack2/dsytrs_aa'
  load_mttmp_ :: ] 'math/mt/test/lapack2/zhetrs'
  load_mttmp_ :: ] 'math/mt/test/lapack2/zhetrs_aa'

  'A Xv'=. y
  Xh=. |: Xv

  rcond=. hecon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  'pL1T pU1T'=. (hetrfpl ; <@hetrfpu) A

  Bax=.      A  mp Xv
  Batx=. (|: A) mp Xv
  Bxa=.  Xh mp    A
  Bxat=. Xh mp |: A

  vferrv=: normitc t04m  NB. for Xv
  vferrh=: normitr t04m  NB. for Xh

  NB. vberrX for Xv at right side
  vberrax=:   mp~     t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc
  NB. vberrX for Xh at left side
  vberrxa=:   mp      t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  ('''l''&dsytrs_mttmp_'    tmonad ((1&{ ,~ 'l' dsytrf_mttmp_    0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A
  ('''u''&dsytrs_mttmp_'    tmonad ((1&{ ,~ 'u' dsytrf_mttmp_    0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A
  ('''l''&zhetrs_mttmp_'    tmonad ((1&{ ,~ 'l' zhetrf_mttmp_    0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A
  ('''u''&zhetrs_mttmp_'    tmonad ((1&{ ,~ 'u' zhetrf_mttmp_    0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A

  ('''l''&dsytrs_aa_mttmp_' tmonad ((1&{ ,~ 'l' dsytrf_aa_mttmp_ 0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A
  ('''u''&dsytrs_aa_mttmp_' tmonad ((1&{ ,~ 'u' dsytrf_aa_mttmp_ 0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A
  ('''l''&zhetrs_aa_mttmp_' tmonad ((1&{ ,~ 'l' zhetrf_aa_mttmp_ 0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A
  ('''u''&zhetrs_aa_mttmp_' tmonad ((1&{ ,~ 'u' zhetrf_aa_mttmp_ 0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A

  ('hetrsplx'               tdyad  ((_3&{.)`(1&{::)                    `]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A  ; pL1T
  ('hetrspltx'              tdyad  ((_3&{.)`(1&{::)                    `]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcond ; norm1At ; pL1T
  ('hetrsxpl'               tdyad  ((_3&{.)`(1&{::)                    `]`(3&{::)`vferrh`vberrxa )) A ; Bxa  ; Xh ; rcond ; normiA  ; pL1T
  ('hetrsxplt'              tdyad  ((_3&{.)`(1&{::)                    `]`(3&{::)`vferrh`vberrxat)) A ; Bxat ; Xh ; rcond ; normiAt ; pL1T

  ('hetrspux'               tdyad  ((_3&{.)`(1&{::)                    `]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A  ; pU1T
  ('hetrsputx'              tdyad  ((_3&{.)`(1&{::)                    `]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcond ; norm1At ; pU1T
  ('hetrsxpu'               tdyad  ((_3&{.)`(1&{::)                    `]`(3&{::)`vferrh`vberrxa )) A ; Bxa  ; Xh ; rcond ; normiA  ; pU1T
  ('hetrsxput'              tdyad  ((_3&{.)`(1&{::)                    `]`(3&{::)`vferrh`vberrxat)) A ; Bxat ; Xh ; rcond ; normiAt ; pU1T

  coerase < 'mttmp'
  erase 'vferrv vferrh vberrax vberratx vberrxa vberrxat'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpotrs1
NB.
NB. Description:
NB.   Test potrsxxx by Hermitian (symmetric) positive
NB.   definite matrix and single RHS
NB.
NB. Syntax:
NB.   testpotrs1 (A ; x)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric) positive
NB.       definite
NB.   x - n-vector, the exact solution

testpotrs1=: 3 : 0
  'A x'=. y

  rcond=. pocon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  'L U'=. (potrfl ,: potrfu) A

  Bax=.      A  mp x
  Batx=. (|: A) mp x
  Bxa=.  x mp    A
  Bxat=. x mp |: A

  NB. vberrX for x at right side
  vberrax=:   mp~     t02v
  vberratx=: (mp~ |:) t02v
  NB. vberrX for x at left side
  vberrxa=:   mp      t02v
  vberrxat=: (mp  |:) t02v

  ('potrslx'  tdyad ((5&{::)`(1&{::)`]`(3&{::)`t04v`vberrax )) A ; Bax  ; x ; rcond ; norm1A  ; L
  ('potrsltx' tdyad ((5&{::)`(1&{::)`]`(3&{::)`t04v`vberratx)) A ; Batx ; x ; rcond ; norm1At ; L
  ('potrsxl'  tdyad ((5&{::)`(1&{::)`]`(3&{::)`t04v`vberrxa )) A ; Bxa  ; x ; rcond ; normiA  ; L
  ('potrsxlt' tdyad ((5&{::)`(1&{::)`]`(3&{::)`t04v`vberrxat)) A ; Bxat ; x ; rcond ; normiAt ; L

  ('potrsux'  tdyad ((5&{::)`(1&{::)`]`(3&{::)`t04v`vberrax )) A ; Bax  ; x ; rcond ; norm1A  ; U
  ('potrsutx' tdyad ((5&{::)`(1&{::)`]`(3&{::)`t04v`vberratx)) A ; Batx ; x ; rcond ; norm1At ; U
  ('potrsxu'  tdyad ((5&{::)`(1&{::)`]`(3&{::)`t04v`vberrxa )) A ; Bxa  ; x ; rcond ; normiA  ; U
  ('potrsxut' tdyad ((5&{::)`(1&{::)`]`(3&{::)`t04v`vberrxat)) A ; Bxat ; x ; rcond ; normiAt ; U

  coerase < 'mttmp'
  erase 'vberrax vberratx vberrxa vberrxat'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpotrs3
NB.
NB. Description:
NB.   Test:
NB.   - xPOTRS (math/lapack2 addon)
NB.   - potrsxxx (math/mt addon)
NB.   by Hermitian (symmetric) positive definite matrix and
NB.   multiple RHS
NB.
NB. Syntax:
NB.   testpotrs3 (A ; X)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric) positive
NB.       definite
NB.   X - n×3-matrix, exact solutions

testpotrs3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/potrf'
  load_mttmp_ :: ] 'math/mt/test/lapack2/potrs'

  'A Xv'=. y
  Xh=. |: Xv

  rcond=. pocon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  'L U'=. (potrfl ,: potrfu) A

  Bax=.      A  mp Xv
  Batx=. (|: A) mp Xv
  Bxa=.  Xh mp    A
  Bxat=. Xh mp |: A

  vferrv=: normitc t04m  NB. for Xv
  vferrh=: normitr t04m  NB. for Xh

  NB. vberrX for Xv at right side
  vberrax=:   mp~     t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc
  NB. vberrX for Xh at left side
  vberrxa=:   mp      t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  ('''l''&dpotrs_mttmp_' tmonad ((1&{ ;~ 'l' dpotrf_mttmp_ 0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A
  ('''u''&dpotrs_mttmp_' tmonad ((1&{ ;~ 'u' dpotrf_mttmp_ 0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A
  ('''l''&zpotrs_mttmp_' tmonad ((1&{ ;~ 'l' zpotrf_mttmp_ 0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A
  ('''u''&zpotrs_mttmp_' tmonad ((1&{ ;~ 'u' zpotrf_mttmp_ 0&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A

  ('potrslx'             tdyad  ((5&{::)`(1&{::)                 `]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A  ; L
  ('potrsltx'            tdyad  ((5&{::)`(1&{::)                 `]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcond ; norm1At ; L
  ('potrsxl'             tdyad  ((5&{::)`(1&{::)                 `]`(3&{::)`vferrh`vberrxa )) A ; Bxa  ; Xh ; rcond ; normiA  ; L
  ('potrsxlt'            tdyad  ((5&{::)`(1&{::)                 `]`(3&{::)`vferrh`vberrxat)) A ; Bxat ; Xh ; rcond ; normiAt ; L

  ('potrsux'             tdyad  ((5&{::)`(1&{::)                 `]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A  ; U
  ('potrsutx'            tdyad  ((5&{::)`(1&{::)                 `]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcond ; norm1At ; U
  ('potrsxu'             tdyad  ((5&{::)`(1&{::)                 `]`(3&{::)`vferrh`vberrxa )) A ; Bxa  ; Xh ; rcond ; normiA  ; U
  ('potrsxut'            tdyad  ((5&{::)`(1&{::)                 `]`(3&{::)`vferrh`vberrxat)) A ; Bxat ; Xh ; rcond ; normiAt ; U

  coerase < 'mttmp'
  erase 'vferrv vferrh vberrax vberratx vberrxa vberrxat'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpttrs1
NB.
NB. Description:
NB.   Test pttrsxxx by Hermitian (symmetric) positive
NB.   definite tridiagonal matrix and single RHS
NB.
NB. Syntax:
NB.   testpttrs1 (A ; x)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric) positive
NB.       definite tridiagonal
NB.   x - n-vector, the exact solution
NB.
NB. TODO:
NB. - A would be sparse

testpttrs1=: 3 : 0
  'A x'=. y

  rcond=. ptcon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  'L1D U1D'=. (pttrfl ; <@pttrfu) A

  Bax=.      A  mp x
  Batx=. (|: A) mp x
  Bxa=.  x mp    A
  Bxat=. x mp |: A

  NB. vberrX for x at right side
  vberrax=:   mp~     t02v
  vberratx=: (mp~ |:) t02v
  NB. vberrX for x at left side
  vberrxa=:   mp      t02v
  vberrxat=: (mp  |:) t02v

  ('pttrslx'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrax )) A ; Bax  ; x ; rcond ; norm1A  ; L1D
  ('pttrsltx' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberratx)) A ; Batx ; x ; rcond ; norm1At ; L1D
  ('pttrsxl'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxa )) A ; Bxa  ; x ; rcond ; normiA  ; L1D
  ('pttrsxlt' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxat)) A ; Bxat ; x ; rcond ; normiAt ; L1D

  ('pttrsux'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrax )) A ; Bax  ; x ; rcond ; norm1A  ; U1D
  ('pttrsutx' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberratx)) A ; Batx ; x ; rcond ; norm1At ; U1D
  ('pttrsxu'  tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxa )) A ; Bxa  ; x ; rcond ; normiA  ; U1D
  ('pttrsxut' tdyad ((_2&{.)`(1&{::)`]`(3&{::)`t04v`vberrxat)) A ; Bxat ; x ; rcond ; normiAt ; U1D

  coerase < 'mttmp'
  erase 'vberrax vberratx vberrxa vberrxat'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpttrs3
NB.
NB. Description:
NB.   Test:
NB.   - xPTTRS (math/lapack2 addon)
NB.   - pttrsxxx (math/mt addon)
NB.   by Hermitian (symmetric) positive definite tridiagonal
NB.   matrix and multiple RHS
NB.
NB. Syntax:
NB.   testpttrs3 (A ; X)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric) positive
NB.       definite tridiagonal
NB.   X - n×3-matrix, the exact solutions
NB.
NB. TODO:
NB. - A would be sparse

testpttrs3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/pttrf'
  load_mttmp_ :: ] 'math/mt/test/lapack2/pttrs'

  'A Xv'=. y
  Xh=. |: Xv

  rcond=. ptcon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  'L1D U1D'=. (pttrfl ; <@pttrfu) A

  Bax=.      A  mp Xv
  Batx=. (|: A) mp Xv
  Bxa=.  Xh mp    A
  Bxat=. Xh mp |: A

  vferrv=: normitc t04m  NB. for Xv
  vferrh=: normitr t04m  NB. for Xh

  NB. vberrX for Xv at right side
  vberrax=:   mp~     t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc
  NB. vberrX for Xh at left side
  vberrxa=:   mp      t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  ('dpttrs_mttmp_'       tmonad ((1&{ ,~ dpttrf_mttmp_@(diag ; _1&diag)@(0&{::))`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A
  ('''l''&zpttrs_mttmp_' tmonad ((1&{ ,~ zpttrf_mttmp_@(diag ; _1&diag)@(0&{::))`]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A

  ('pttrslx'             tdyad  ((_2&{.)`(1&{::)                                `]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A  ; L1D
  ('pttrsltx'            tdyad  ((_2&{.)`(1&{::)                                `]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcond ; norm1At ; L1D
  ('pttrsxl'             tdyad  ((_2&{.)`(1&{::)                                `]`(3&{::)`vferrh`vberrxa )) A ; Bxa  ; Xh ; rcond ; normiA  ; L1D
  ('pttrsxlt'            tdyad  ((_2&{.)`(1&{::)                                `]`(3&{::)`vferrh`vberrxat)) A ; Bxat ; Xh ; rcond ; normiAt ; L1D

  ('pttrsux'             tdyad  ((_2&{.)`(1&{::)                                `]`(3&{::)`vferrv`vberrax )) A ; Bax  ; Xv ; rcond ; norm1A  ; U1D
  ('pttrsutx'            tdyad  ((_2&{.)`(1&{::)                                `]`(3&{::)`vferrv`vberratx)) A ; Batx ; Xv ; rcond ; norm1At ; U1D
  ('pttrsxu'             tdyad  ((_2&{.)`(1&{::)                                `]`(3&{::)`vferrh`vberrxa )) A ; Bxa  ; Xh ; rcond ; normiA  ; U1D
  ('pttrsxut'            tdyad  ((_2&{.)`(1&{::)                                `]`(3&{::)`vferrh`vberrxat)) A ; Bxat ; Xh ; rcond ; normiAt ; U1D

  coerase < 'mttmp'
  erase 'vferrv vferrh vberrax vberratx vberrxa vberrxat'

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
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testtrs_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testtrs_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testtrs_mt_ 150 150
NB.
NB. Notes:
NB. - nrhs=3 is assumed

testtrs=: 1 : 'EMPTY [ ((((u ptmat2_mt_)@# (testpttrs3_mt_@; [ testpttrs1_mt_@(; {:"1)) ]) [ ((u pomat_mt_)@# (testpotrs3_mt_@; [ testpotrs1_mt_@(; {:"1)) ]) [ ((u hemat_mt_)@# (testhetrs3_mt_@; [ testhetrs1_mt_@(; {:"1)) ]))@u@({. , 3:) [ ((}."1 (testgetrs3_mt_@; [ testgetrs1_mt_@(; {:"1)) {."1)~ _3:)@u@(+&0 3))^:(=/)'
