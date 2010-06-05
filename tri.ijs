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
NB. trtril         L  * invL  = I       invL=.  trtril  A
NB. trtril1        L1 * invL1 = I       invL1=. trtril1 A
NB. trtriu         U  * invU  = I       invU=.  trtriu  A
NB. trtriu1        U1 * invU1 = I       invU1=. trtriu1 A
NB.
NB. Description:
NB.   Inverse triangular matrix
NB. where:
NB.   A     - n×n-matrix, containing triangular matrix
NB.   U     - n×n upper triangular matrix
NB.   U1    - n×n unit upper triangular matrix (diagonal is
NB.           not saved)
NB.   L     - n×n lower triangular matrix
NB.   L1    - n×n unit lower triangular matrix (diagonal is
NB.           not saved)
NB.   invU  - n×n upper triangular matrix, inversion of U
NB.   invU1 - n×n unit upper triangular matrix (diagonal is
NB.           not saved), inversion of U1
NB.   invL  - n×n lower triangular matrix, inversion of L
NB.   invL1 - n×n unit lower triangular matrix (diagonal is
NB.           not saved), inversion of L1
NB.   n     ≥ 0
NB.
NB. Notes:
NB. - models LAPACK's xTRTRI, but uses recursive not
NB.   partitioned algorithm
NB. - opposite triangle is not referenced
NB. - unit diagonal is not referenced

trtril=: 3 : 0
  n=. # y
  if. n > 1 do.
    k=. >. -: n
    Ta=. (2 # k) }. y
    invTa=. trtril Ta
    invTb=. trtril (2 # k) {. y
    Ac=. (k ,~ (k - n)) {. y
    invAc=. - invTa mp Ac mp invTb
    (invTa ,.~ invAc) ( 0 append~) invTb
  else.
    % y
  end.
)

trtril1=: 3 : 0
  n=. # y
  if. n > 1 do.
    k=. >. -: n
    Ta=. (2 # k) }. y
    invTa=. trtril1 Ta
    invTb=. trtril1 (2 # k) {. y
    Ac=. (k ,~ (k - n)) {. y
    invAc=. - invTa mp Ac mp invTb
    (invTa ,.~ invAc) ( 0 append~) invTb
  else.
    (1:"0) y
  end.
)

trtriu=: 3 : 0
  n=. # y
  if. n > 1 do.
    k=. >. -: n
    Ta=. (2 # k) {. y
    invTa=. trtriu Ta
    invTb=. trtriu (2 # k) }. y
    Ac=. (k , (k - n)) {. y
    invAc=. - invTa mp Ac mp invTb
    (invTa ,. invAc) (_1 append) invTb
  else.
    % y
  end.
)

trtriu1=: 3 : 0
  n=. # y
  if. n > 1 do.
    k=. >. -: n
    Ta=. (2 # k) {. y
    invTa=. trtriu1 Ta
    invTb=. trtriu1 (2 # k) }. y
    Ac=. (k , (k - n)) {. y
    invAc=. - invTa mp Ac mp invTb
    (invTa ,. invAc) (_1 append) invTb
  else.
    (1:"0) y
  end.
)

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
NB.   ip  - n-vector, full inversed permutation of A, as
NB.         returned by hetrfpx
NB.   L1  - n×n-matrix, unit lower triangular, as returned by
NB.         hetrfpl
NB.   U1  - n×n-matrix, unit upper triangular, as returned by
NB.         hetrfpu
NB.   T   - n×n-matrix, Hermitian (symmetric) 3-diagonal, as
NB.         returned by hetrfpx
NB.   iA  - n×n-matrix, inversion of A
NB.   P   - n×n-matrix, full permutation of A
NB.   A   - n×n-matrix to inverse
NB.
NB. Notes:
NB. - models LAPACK's xHETRI, but uses Aasen factorization
NB.   instead of Bunch-Kaufman one

hetripl=: (0 & {::) sp (pttril @ pttrfl @ (2 & {::)) ((ct @ ]) mp mp) (trtril1 @ (1 & {::))
hetripu=: (0 & {::) sp (pttril @ pttrfl @ (2 & {::)) ((ct @ ]) mp mp) (trtriu1 @ (1 & {::))

NB. quick-and-dirty iterative variant - doesn't work

NB. 'pfxipi1 pfxL1i1 pfxTi1 sfxipi1 sfxL1i1 sfxTi1'=. hetripl (pfxipi;pfxL1i;pfxTi;sfxipi;sfxL1i;sfxTi)
zzzhetriplstep=: 3 : 0
  'pfxipi pfxL1i pfxTi sfxipi sfxL1i sfxTi'=. y
  'n k'=. $ pfxL1i
  akk=. % {. 0 _1 1 diag pfxTi        NB. last from main diagonal
  ak1n=. (< _1 ; n ht2lios k) { pfxL1i
  ak1nupd=. sfxL1i mp - ak1n
  akkupd=. akk - 9 o. (+ ak1n) mp ak1nupd
  ipk=. {: pfxipi
  dp0=. 0 (lios2cp`(a:"_) @. (*. & (0&=))) ipk - (k-1)
  dpk=. (k-1) lios2cp ipk
  pfxipi=. }: pfxipi
  pfxL1i=. }:"1 pfxL1i
  pfxTi=. }:"1 pfxTi
  sfxipi=. ipk , sfxipi
  sfxL1i=. (1 , + ak1nupd) , ak1nupd ,. sfxL1i
  tdiag=. akk , diag sfxTi
  sfxTi=. ((dp0 (C. dbg 'C.') tdiag) ; a:) (setdiag dbg 'setdiag') (2 # 1+k-n) {. sfxTi
  pfxipi ; pfxL1i ; pfxTi ; sfxipi ; sfxL1i ; sfxTi
)

NB. iA=. hetripl (ip ; L1 ; T)
zzzhetripl=: 3 : 0
  'ip L1 T'=. y
  (hetriplstep dbg 'STEP') ^: (# ip) y , (i. 0) ; EMPTY ; EMPTY
)

NB. ---------------------------------------------------------
NB. potri
NB. Inverse a Hermitian (symmetric) positive definite matrix
NB. using the Cholesky factorization:
NB.   L * L' = A
NB.   inv(A) = inv(L)' * inv(L)

potri=: (mp~ ct) @ trtril

NB. ---------------------------------------------------------
NB. pttril
NB. iA=. pttril (L1;D)

pttril=: pttrsax idmat @ (0 & {::)

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
  rcondU=. (norm1 con trtriu) U
  rcondU1=. (norm1 con trtriu1) U1
  rcondL=. (norm1 con trtril) y
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
  ipL1U=. getrfpl1u y
  ipUL1=. getrful1p y

  ('%.'        tmonad (]`]`(rcond"_)`(_."_)`(  (norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@])))                                     )) y

  ('getripl1u' tmonad (]`]`(rcond"_)`(_."_)`((((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@])))~ ((C.^:_1   (trl1 mp tru )) & >/))~))) ipL1U
  ('getriul1p' tmonad (]`]`(rcond"_)`(_."_)`((((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@])))~ ((C.^:_1"1 (tru  mp trl1)) & >/))~))) ipUL1

  EMPTY
)


NB. ---------------------------------------------------------
NB. testhetri
NB. Test inverse algorithms with random Hermitian (symmetric)
NB. matrix y
NB.
NB. thetri A

testhetri=: 3 : 0
  n=. # y
  GE=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y
  rcondHE=. norm1 con hetri HE=. ge2he GE
  'hetri' ttri (HE;rcondHE)
  EMPTY
)

NB. ---------------------------------------------------------
NB. testpotri
NB. Test inverse algorithms with random Hermitian (symmetric)
NB. positive definite matrix y
NB.
NB. tpotri A

testpotri=: 3 : 0
  n=. # y
  GE=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y
  rcondPO=. norm1 con potri   PO=. (mp ct) GE  NB. >>> ###### PO=. ge2po GE
  'potri'    ttri (PO;rcondPO)
  EMPTY
)

NB. ---------------------------------------------------------
NB. testpttri
NB.
NB. Description:
NB.   Test triangular inversion algorithms:
NB.   - pttri (math/mt addon)
NB.   by Hermitian (symmetric) positive definite tridiagonal
NB.   matrix given
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
  rcond=. (norm1 con (pttri@pttrfl)) y

  ('pttri' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((mp mp (ct@[))&>/)))) % (FP_EPS*((norm1*c)@[))))) y

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
