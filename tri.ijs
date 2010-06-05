NB. Inverse by triangular factorization
NB.
NB. trtrixx    Inverse triangular matrix
NB. getrixxxx  Inverse general matrix
NB. hetripx    Inverse Hermitian (symmetric) matrix
NB. potrix     Inverse Hermitian (symmetric) positive
NB.            definite matrix
NB. pttrix     Inverse Hermitian (symmetric) positive
NB.            definite tridiagonal matrix
NB.
NB. testtrtri  Test trtrixx by triangular matrix given
NB. testgetri  Test getrixxxx by general matrix given
NB. testhetri  Test hetripx by Hermitian (symmetric) matrix
NB.            given
NB. testpotri  Test potrix by Hermitian (symmetric) positive
NB.            definite matrix given
NB. testpttri  Test pttrix by Hermitian (symmetric) positive
NB.            definite tridiagonal matrix given
NB. testtri    Adv. to make verb to test xxtrixx by matrix of
NB.            generator and shape given
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

TRINB=: 64  NB. block size limit

NB. ---------------------------------------------------------
NB. getripl1ustep
NB.
NB. Description:
NB.   Single step of non-blocked version of getripl1u
NB.
NB. Syntax:
NB.   'pfxi1 sfxi1'=. getripl1ustep (pfxi ; sfxi)
NB. where
NB.   pfxi  - pfx(i), n×(j+TRINB)-matrix after i-th step and
NB.           before (i+1)-th one, contains not yet processed
NB.           part
NB.   sfxi  - sfx(i), n×(n-j-TRINB)-matrix after i-th step
NB.           and before (i+1)-th one, contains already
NB.           processed part
NB.   pfxi1 - pfx(i+1), n×j-matrix, being pfx(i) after
NB.           (i+1)-th step
NB.   sfxi1 - sfx(i+1), n×(n-j)-matrix, being sfx(i) after
NB.           (i+1)-th step
NB.
NB. References:
NB. [1] J. DuCroz, N. Higham. Stability of Methods for Matrix
NB.     Inversion, UT-CS-90-119, October, 1990.
NB.     LAPACK Working Note 27.
NB.     http://www.netlib.org/lapack/lawns/downloads/

getripl1ustep=: 3 : 0
  'pfx sfx'=. y
  'n j'=. $ pfx
  j=. j - TRINB
  L0=. (,.~ j , TRINB) (];.0) pfx
  L1=. (- (c sfx),TRINB) {. pfx
  C=. (-TRINB) {."1 pfx
  C=. (((i. n) - j) >/ i. TRINB) } C ,: 0  NB. spec code
  C=. C - sfx mp L1
  C=. L0 trsmxl1 C
  ((-TRINB) }."1 pfx) ; (C ,. sfx)
)

NB. ---------------------------------------------------------
NB. getriul1pstep
NB.
NB. Description:
NB.   Single step of non-blocked version of getriul1p
NB.
NB. Syntax:
NB.   'pfxi1 sfxi1'=. getristep (pfxi ; sfxi)
NB. where
NB.   pfxi  - pfx(i), j×n-matrix after i-th step and before
NB.           (i+1)-th one, contains already processed part
NB.   sfxi  - sfx(i), (n-j)×n-matrix after i-th step and
NB.           before (i+1)-th one, contains not yet processed
NB.           part
NB.   pfxi1 - pfx(i+1), (j+TRINB)×n-matrix, being pfx(i)
NB.           after (i+1)-th step
NB.   sfxi1 - sfx(i+1), (n-j-TRINB)×n-matrix, being sfx(i)
NB.           after (i+1)-th step
NB.
NB. References:
NB. [1] J. DuCroz, N. Higham. Stability of Methods for Matrix
NB.     Inversion, UT-CS-90-119, October, 1990.
NB.     LAPACK Working Note 27.
NB.     http://www.netlib.org/lapack/lawns/downloads/

getriul1pstep=: 3 : 0
  'pfx sfx'=. y
  'j n'=. $ pfx
  L0=. ((0 , j) ,: (2 # TRINB)) (];.0) sfx
  L1=. (TRINB , j) {. sfx
  R=. TRINB {. sfx
  R=. ((i. TRINB) >/ ((i. n) - j)) } R ,: 0  NB. spec code
  R=. R - L1 mp pfx
  R=. L0 trsml1x R
  (pfx , R) ; (TRINB }. sfx)
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:          Solves:              Syntax:
NB. trtril         L  * iL  = I         iL=.  trtril  A
NB. trtril1        L1 * iL1 = I         iL1=. trtril1 A
NB. trtriu         U  * iU  = I         iU=.  trtriu  A
NB. trtriu1        U1 * iU1 = I         iU1=. trtriu1 A
NB.
NB. Description:
NB.   Inverse triangular matrix
NB. where:
NB.   A   - n×n-matrix, containing triangular matrix
NB.   U   - n×n upper triangular matrix
NB.   U1  - n×n unit upper triangular matrix (diagonal is not
NB.         saved)
NB.   L   - n×n lower triangular matrix
NB.   L1  - n×n unit lower triangular matrix (diagonal is not
NB.         saved)
NB.   iU  - n×n upper triangular matrix, inversion of U
NB.   iU1 - n×n unit upper triangular matrix (diagonal is not
NB.         saved), inversion of U1
NB.   iL  - n×n lower triangular matrix, inversion of L
NB.   iL1 - n×n unit lower triangular matrix (diagonal is not
NB.         saved), inversion of L1
NB.   n   ≥ 0
NB.
NB. Algorithm for trtriu:
NB.   In: A
NB.   Out: iU
NB.   0) if 1 < # A
NB.      0.0) then
NB.           0.0.0) form (#A)-vector of zeros:
NB.                    fret=. (#A) $ 0
NB.           0.0.1) find splitting edge:
NB.                    k=. >. -: # A
NB.           0.0.2) mark intervals:
NB.                    fret=. 1 (0,k)} fret
NB.           0.0.3) cut A by fret to block matrix bA:
NB.                    bA = (  A00  A01  )  k
NB.                         (  A10  A11  )  n-k
NB.                            k    n-k
NB.           0.0.4) apply trtriu itself to A00 and A11:
NB.                    iA00=. $: A00
NB.                    iA11=. $: A11
NB.           0.0.5) replace A00 and A11 by iA00 and iA11,
NB.                  respectively, in bA
NB.           0.0.6) inverse A01:
NB.                    iA01=. - iA00 mp A01 mp iA11
NB.           0.0.7) replace A01 by iA01 in bA
NB.           0.0.8) assemble iA from block matrix:
NB.                    iA=. icut bA
NB.     0.1) else return reciprocal:
NB.            iA=. % A
NB.
NB. Notes:
NB. - models LAPACK's xTRTRI, but uses blocked not
NB.   partitioned algorithm
NB. - opposite triangle is not referenced
NB. - unit diagonal is not referenced
NB.
NB. References:
NB. [1] JWiki/Essays/Triangular Matrix Inverse
NB.     Roger Hui
NB.     http://www.jsoftware.com/jwiki/Essays/Triangular%20Matrix%20Inverse
NB.
NB. TODO:
NB. - fret would be sparse

trtril=:  %     `(icut@((2:}~ <@:-@((1 1 {:: ]) mp (1 0 {:: ]) mp 0 0 {:: ]))@(<@$:`<`<;.1~ ;~@((0) 1:`(, >.@-:)`(#~)} #))))@.(1 < #)
trtril1=: (1:"0)`(icut@((2:}~ <@:-@((1 1 {:: ]) mp (1 0 {:: ]) mp 0 0 {:: ]))@(<@$:`<`<;.1~ ;~@((0) 1:`(, >.@-:)`(#~)} #))))@.(1 < #)
trtriu=:  %     `(icut@((1:}~ <@:-@((0 0 {:: ]) mp (0 1 {:: ]) mp 1 1 {:: ]))@(<@$:`<`<;.1~ ;~@((0) 1:`(, >.@-:)`(#~)} #))))@.(1 < #)
trtriu1=: (1:"0)`(icut@((1:}~ <@:-@((0 0 {:: ]) mp (0 1 {:: ]) mp 1 1 {:: ]))@(<@$:`<`<;.1~ ;~@((0) 1:`(, >.@-:)`(#~)} #))))@.(1 < #)

NB. ---------------------------------------------------------
NB. Verb:      Factorization used:  Syntax:
NB. getripl1u  P * L1 * U = A       iA=. getripl1u (ip ; L1U)
NB. getriul1p  U * L1 * P = A       iA=. getriul1p (ip ; UL1)
NB.
NB. Description:
NB.   Inverse a general matrix using the triangular
NB.   factorization with partial pivoting
NB. where
NB.   ip  - n-vector, partial inversed permutation of A, as
NB.         returned by getrfxxxx
NB.   L1U - n×n-matrix, contains U and L1, as returned by
NB.         getrfpl1u
NB.   UL1 - n×n-matrix, contains U and L1, as returned by
NB.         getrful1p
NB.   iA  - n×n-matrix, inversion of A
NB.   P   - n×n-matrix, partial permutation of A
NB.   L1  - n×n-matrix, unit lower triangular
NB.   U   - n×n-matrix, upper triangular
NB.   A   - n×n-matrix to inverse
NB.
NB. Notes:
NB. - getripl1u implements LAPACK's xGETRI
NB. - getriul1p models LAPACK's xGETRI, but uses row blocks
NB.   and another triangular factorization
NB.
NB. References:
NB. [1] J. DuCroz, N. Higham. Stability of Methods for Matrix
NB.     Inversion, UT-CS-90-119, October, 1990.
NB.     LAPACK Working Note 27.
NB.     http://www.netlib.org/lapack/lawns/downloads/

getripl1u=: 3 : 0
  'ip L1U'=. y
  n=. # L1U
  y=. trtriu L1U
  y=. (>/~ i. n) } y ,: L1U  NB. spec code
  I=. 0 >. <. (n-1) % TRINB
  ip (C.^:_1"1) 1 {:: getripl1ustep ^: I ((TRINB * I) (({."1) ; (-@[ (trl1 trsmxl1 tru) (}."1))) y)
)

getriul1p=: 3 : 0
  'ip UL1'=. y
  n=. # UL1
  y=. trtriu UL1
  y=. (>/~ i. n) } y ,: UL1          NB. spec code
  rn=. (TRINB+1) (| &. <:) (1 >. n)  NB. reminder of n, (n-nn): if. n>0 then. (NB|(n-1))+1 else. 0 end.
  I=. <. TRINB %~ <: (1 >. n)        NB. if. n>0 then. max(0,floor((n-1)%NB)) else. 0 end.
  ip (C.^:_1) 0 {:: getriul1pstep ^: I (rn ((((2 # [) {. ]) trsml1x (tru@{.)) ; }.) y)
)

NB. ---------------------------------------------------------
NB. Verb:    Factorization used:           Syntax:
NB. hetripl  P * L1 * T * L1^H * P^_1 = A  iA=. hetripl (ip ; L1 ; T)
NB. hetripu  P * U1 * T * U1^H * P^_1 = A  iA=. hetripu (ip ; U1 ; T)
NB.
NB. Description:
NB.   Inverse Hermitian (symmetric) matrix using the
NB.   triangular factorization with full pivoting
NB. where
NB.   ip - n-vector, full inversed permutation of A, as
NB.        returned by hetrfpx
NB.   L1 - n×n-matrix, unit lower triangular, as returned by
NB.        hetrfpl
NB.   U1 - n×n-matrix, unit upper triangular, as returned by
NB.        hetrfpu
NB.   T  - n×n-matrix, Hermitian (symmetric) 3-diagonal, as
NB.        returned by hetrfpx
NB.   iA - n×n-matrix, inversion of A
NB.   P  - n×n-matrix, full permutation of A
NB.   A  - n×n-matrix to inverse, Hermitian (symmetric)
NB.
NB. Notes:
NB. - models LAPACK's xHETRI, but uses Aasen factorization
NB.   instead of Bunch-Kaufman one

hetripl=: (0 & {::) sp (httril @ httrfl @ (2 & {::)) ((ct @ ]) mp mp) (trtril1 @ (1 & {::))  NB. IMPLEMENTME
hetripu=: (0 & {::) sp (httril @ httrfl @ (2 & {::)) ((ct @ ]) mp mp) (trtriu1 @ (1 & {::))  NB. IMPLEMENTME

NB. ---------------------------------------------------------
NB. Verb:     Factorization used:          Syntax:
NB. potril    L * L^H = A                  iA=. potril L
NB. potriu    U * U^H = A                  iA=. potriu U
NB.
NB. Description:
NB.   Inverse Hermitian (symmetric) positive definite matrix
NB.   using the triangular factorization
NB. where
NB.   L  - n×n-matrix, lower triangular, Cholesky triangle as
NB.        returned by potrfl
NB.   U  - n×n-matrix, upper triangular, Cholesky triangle as
NB.        returned by potrfu
NB.   iA - n×n-matrix, inversion of A
NB.   A  - n×n-matrix to inverse, Hermitian (symmetric)
NB.        positive definite
NB.
NB. Notes:
NB. - models LAPACK's xPOTRI, but uses direct, not iterative
NB.   matrix product

potril=: (mp~ ct) @ trtril
potriu=: (mp~ ct) @ trtriu

NB. ---------------------------------------------------------
NB. Verb:     Factorization used:        Syntax:
NB. pttril    L1 * D * L1^H = A          iA=. pttril (L1 ; D)
NB. pttriu    U1 * D * U1^H = A          iA=. pttriu (U1 ; D)
NB.
NB. Description:
NB.   Inverse Hermitian (symmetric) positive definite
NB.   tridiagonal matrix using the triangular factorization
NB. where
NB.   L1 - n×n-matrix, unit lower triangular, as returned by
NB.        pttrfl
NB.   U1 - n×n-matrix, unit upper triangular, as returned by
NB.        pttrfu
NB.   D  - n×n-matrix, diagonal with positive diagonal
NB.        entries, as returned by pttrfx
NB.   iA - n×n-matrix, inversion of A
NB.   A  - n×n-matrix to inverse, Hermitian (symmetric)
NB.        positive definite tridiagonal

pttril=: pttrsax idmat @ # @ (0 & {::)
pttriu=: [:                             NB. pttrs version accepting U1*D*U1^H is not implemented yet

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testtrtri
NB.
NB. Description:
NB.   Test:
NB.   - 128!:1 (built-in)
NB.   - trtrixx (math/mt addon)
NB.   by triangular matrix given
NB.
NB. Syntax:
NB.   testtrtri A
NB. where
NB.   A - n×n-matrix, lower triangular
NB.
NB. Formula:
NB. - berr := ||I - A * A^_1|| / (ε * ||A|| * ||A^_1|| * n)

testtrtri=: 3 : 0
  L1=. |: U1=. tru1 U=. |: y
  rcondU=.  (norm1 con trtriu ) U
  rcondU1=. (norm1 con trtriu1) U1
  rcondL=.  (norm1 con trtril ) y
  rcondL1=. (norm1 con trtril1) L1

  ('(128!:1)'  tmonad (]`]`(rcondU "_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) U

  ('trtril'    tmonad (]`]`(rcondL "_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) y
  ('trtril1'   tmonad (]`]`(rcondL1"_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) L1
  ('trtriu'    tmonad (]`]`(rcondU "_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) U
  ('trtriu1'   tmonad (]`]`(rcondU1"_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) U1

  EMPTY
)

NB. ---------------------------------------------------------
NB. testgetri
NB.
NB. Description:
NB.   Test getrixxxx by general matrix given
NB.
NB. Syntax:
NB.   testgetri A
NB. where
NB.   A - n×n-matrix
NB.
NB. Formula:
NB. - berr := ||I - A * A^_1|| / (ε * ||A|| * ||A^_1|| * n)

testgetri=: 3 : 0
  rcond=. (norm1 con (getriul1p@getrful1p)) y

  ('%.'        tmonad (] `]`(rcond"_)`(_."_)`(           (norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@])))   ))               y

  ('getripl1u' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; getrfpl1u) y
  ('getriul1p' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; getrful1p) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhetri
NB.
NB. Description:
NB.   Test hetripx by Hermitian (symmetric) matrix given
NB.
NB. Syntax:
NB.   testhetri A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric)
NB.
NB. Formula:
NB. - berr := ||I - A * A^_1|| / (ε * ||A|| * ||A^_1|| * n)

testhetri=: 3 : 0
  rcond=. _. NB. (norm1 con (hetripl@hetrfpl)) y

  ('hetripl' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; hetrfpl) y
  ('hetripu' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; hetrfpu) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpotri
NB.
NB. Description:
NB.   Test hetripx by Hermitian (symmetric) positive definite
NB.   matrix given
NB.
NB. Syntax:
NB.   testpotri A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.
NB. Formula:
NB. - berr := ||I - A * A^_1|| / (ε * ||A|| * ||A^_1|| * n)

testpotri=: 3 : 0
  rcond=. (norm1 con (potril@potrfl)) y

  ('potril' tmonad ((1 & {::)`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; potrfl) y
  ('potriu' tmonad ((1 & {::)`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; potrfu) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpttri
NB.
NB. Description:
NB.   Test pttrix by Hermitian (symmetric) positive definite
NB.   tridiagonal matrix given
NB.
NB. Syntax:
NB.   testpttri A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive
NB.       definite tridiagonal
NB.
NB. Formula:
NB. - berr := ||I - A * A^_1|| / (ε * ||A|| * ||A^_1|| * n)

testpttri=: 3 : 0
  rcond=. (norm1 con (pttril@pttrfl)) y

  ('pttril' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; pttrfl) y
  ('pttriu' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; pttrfu) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testtri
NB.
NB. Description:
NB.   Adv. to make verb to test triangular inversion
NB.   algorithms by matrix of generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testtri
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
NB.     (? @ $ 0:) testtri_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testtri_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testtri_mt_ 150 200

testtri=: 1 : 'EMPTY_mt_ [ ((testpttri_mt_ @ (u ptmat_mt_)) [ (testpotri_mt_ @ (u pomat_mt_)) [ (testhetri_mt_ @ (u hemat_mt_)) [ (testgetri_mt_ @ u) [ (testtrtri_mt_ @ (u trlmat_mt_))) ^: (=/)'
