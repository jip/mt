NB. Inverse by triangular factorization
NB.
NB. trtrixx    Inverse triangular matrix
NB. getri      Inverse general matrix
NB. hetri      Inverse Hermitian (symmetric) matrix
NB. potri      Inverse Hermitian (symmetric) positive
NB.            definite matrix
NB. pttri      Inverse Hermitian (symmetric) positive
NB.            definite tridiagonal matrix
NB.
NB. testtrtri  Test trtrixx by triangular matrix given
NB. testgetri  Test getri by general matrix given
NB. testhetri  Test hetri by Hermitian (symmetric) matrix
NB.            given
NB. testpotri  Test potri by Hermitian (symmetric) positive
NB.            definite matrix given
NB. testpttri  Test pttri by Hermitian (symmetric) positive
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

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:          Solves:              Syntax:
NB. trtril         L  * invL  = I       invL=.  trtril A
NB. trtril1        L1 * invL1 = I       invL1=. trtril1 A
NB. trtriu         U  * invU  = I       invU=.  trtriu A
NB. trtriu1        U1 * invU1 = I       invU1=. trtriu1 A
NB.
NB. Description:
NB.   inverse triangular matrix
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
NB. - opposite triangle is not referenced
NB. - unit diagonal is not referenced
NB.
NB. References:
NB. [1] J. DuCroz, N. Higham. Stability of Methods for Matrix
NB.     Inversion, UT-CS-90-119, October, 1990.
NB.     LAPACK Working Note 27.
NB.     http://www.netlib.org/lapack/lawns/downloads/

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
NB. getri
NB. Inverse a general matrix using the LU factorization
NB.   inv(P) * L1 * U = A
NB.   inv(A) = inv(U) * inv(L1) * P
NB.
NB. References:
NB. [1] J. DuCroz, N. Higham. Stability of Methods for Matrix
NB.     Inversion, UT-CS-90-119, October, 1990.
NB.     LAPACK Working Note 27.
NB.     http://www.netlib.org/lapack/lawns/downloads/

getri=: ((C."1~ /:)~ (trtriu mp trtril1)) & >/

TRINB=: 64

NB. implements LAPACK's xGETRI with column blocks, P*L1*U=A
getri2step=: 3 : 0
  'pfx sfx'=. y
  n=. # pfx
  j=. (c pfx) - TRINB
  L0=. (,.~ j , TRINB) (];.0) pfx
  L1=. (- (c sfx),TRINB) {. pfx
  bcol=. (-TRINB) {."1 pfx
  bcol=. (((i. n) - j) >/ i. TRINB) } bcol ,: 0       NB. spec code
  bcol=. bcol - sfx mp L1
  bcol=. L0 trsmxl1 bcol
  ((-TRINB) }."1 pfx) ; (bcol ,. sfx)
)

getri2=: 3 : 0
  'ip L1U'=. y
  n=. # L1U
  y=. trtriu L1U
  y=. (>/~ i. n) } y ,: L1U       NB. spec code
  I=. 0 >. <. (n-1) % TRINB
  ip (C.^:_1"1) 1 {:: getri2step ^: I ((TRINB * I) (({."1) ; (-@[ (trl1 trsmxl1 tru) (}."1))) y)
)

NB. models LAPACK's xGETRI with row blocks, U*L1*P=A
NB. TODO: x m&v y ↔ m&v^:x y
getri3step=: 3 : 0
  'pfx sfx'=. y
  'j n'=. $ pfx
  L0=. ((0 , j) ,: (2 # TRINB)) (];.0) sfx
  L1=. (TRINB , j) {. sfx
  R=. TRINB {. sfx
  R=. ((i. TRINB) >/ ((i. n) - j)) } R ,: 0       NB. spec code
  R=. R - L1 mp pfx
  R=. L0 trsml1x R
  (pfx , R) ; (TRINB }. sfx)
)

getri3=: 3 : 0
  'ip UL1'=. y
  n=. # UL1
  y=. trtriu UL1
  y=. (>/~ i. n) } y ,: UL1              NB. spec code
  rn=. (TRINB+1) (| &. <:) (1 >. n)      NB. reminder of n, (n-nn): if. n>0 then. (NB|(n-1))+1 else. 0 end.
  I=. <. TRINB %~ <: (1 >. n)            NB. if. n>0 then. max(0,floor((n-1)%NB)) else. 0 end.
  ip (C.^:_1) 0 {:: getri3step ^: I (rn ((((2 # [) {. ]) trsml1x (tru@{.)) ; }.) y)
)

NB. x m&v y ↔ m&v^:x y
getri4step=: 4 : 0
  'pfx sfx'=. y
  'j n'=. $ pfx
  L0=. ((0 , j) ,: (2 # TRINB)) (];.0) sfx
  L1=. (TRINB , j) {. sfx
  R=. TRINB {. sfx
  R=. ((i. TRINB) >/ ((i. n) - j)) } R ,: 0       NB. spec code
  R=. R - L1 mp pfx
  R=. L0 trsml1x R
  (pfx , R) ; (TRINB }. sfx)
)

getri4=: 3 : 0
  'ip UL1'=. y
  n=. # UL1
  y=. trtriu UL1
  y=. (>/~ i. n) } y ,: UL1              NB. spec code
  rn=. (TRINB+1) (| &. <:) (1 >. n)      NB. reminder of n, (n-nn): if. n>0 then. (NB|(n-1))+1 else. 0 end.
  I=. <. TRINB %~ <: (1 >. n)            NB. if. n>0 then. max(0,floor((n-1)%NB)) else. 0 end.
  ip (C.^:_1) 0 {:: I (0&getri4step) (rn ((((2 # [) {. ]) trsml1x (tru@{.)) ; }.) y)
)

NB. ---------------------------------------------------------
hetri=: getri  NB. stub for a while

NB. ---------------------------------------------------------
NB. potri
NB. Inverse a Hermitian (symmetric) positive definite matrix
NB. using the Cholesky factorization:
NB.   L * L' = A
NB.   inv(A) = inv(L)' * inv(L)

potri=: (mp~ ct) @ trtril

NB. ---------------------------------------------------------
NB. pttri
NB. invA=. pttri (L1;D)

pttri=: pttrs idmat @ (0 & {::)

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
NB.   Test getri by general matrix given
NB.
NB. Syntax:
NB.   testgetri A
NB. where
NB.   A - n×n-matrix
NB.
NB. Formula:
NB. - berr := ||I - A * A^_1|| / (ε * ||A|| * ||A^_1|| * n)

testgetri=: 3 : 0
  rcond=. (norm1 con getri) y
  ipL1U=. getrfpl1u y
  ipUL1=. getrful1p y             NB. TODO: save back to y

  ('%.'     tmonad (]`]`(rcond"_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) y

  ('getri'  tmonad (]`]`(rcond"_)`(_."_)`((((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@])))~ ((C.^:_1   (trl1 mp tru )) & >/))~))) ipL1U
  ('getri2' tmonad (]`]`(rcond"_)`(_."_)`((((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@])))~ ((C.^:_1   (trl1 mp tru )) & >/))~))) ipL1U
  ('getri3' tmonad (]`]`(rcond"_)`(_."_)`((((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@])))~ ((C.^:_1"1 (tru  mp trl1)) & >/))~))) ipUL1
  ('getri4' tmonad (]`]`(rcond"_)`(_."_)`((((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@])))~ ((C.^:_1"1 (tru  mp trl1)) & >/))~))) ipUL1

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
