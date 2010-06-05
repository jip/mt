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

NB. ---------------------------------------------------------
NB. Blocked code constants

TRINB=: 64  NB. block size limit

NB. ---------------------------------------------------------
NB. getristep
NB.
NB. Description:
NB.   Single step of non-blocked version of getri
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

getristep=: 3 : 0
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
NB.
NB. Description:
NB.   Inverse a general matrix using the UL factorization:
NB.     U * L1 * P = A
NB.
NB. Syntax:
NB.   iA=. getri (ip ; UL1)
NB. where
NB.   ip  - n-vector, columns inversed permutation of A, as
NB.         returned by getrful1p
NB.   UL1 - n×n-matrix, contains U and L1, as returned by
NB.         getrful1p
NB.   iA  - n×n-matrix, inversion of A
NB.   P   - n×n-matrix, columns permutation of A
NB.   L   - n×n-matrix, unit lower triangular
NB.   U1  - n×n-matrix, upper triangular
NB.   A   - n×n-matrix to inverse
NB.
NB. Notes:
NB. - models LAPACK's xGETRI, but uses row blocks and another
NB.   triangular factorization
NB.
NB. References:
NB. [1] J. DuCroz, N. Higham. Stability of Methods for Matrix
NB.     Inversion, UT-CS-90-119, October, 1990.
NB.     LAPACK Working Note 27.
NB.     http://www.netlib.org/lapack/lawns/downloads/

getri=: 3 : 0
  'ip UL1'=. y
  n=. # UL1
  y=. trtriu UL1
  y=. (>/~ i. n) } y ,: UL1          NB. spec code
  rn=. (TRINB+1) (| &. <:) (1 >. n)  NB. reminder of n, (n-nn): if. n>0 then. (NB|(n-1))+1 else. 0 end.
  I=. <. TRINB %~ <: (1 >. n)        NB. if. n>0 then. max(0,floor((n-1)%NB)) else. 0 end.
  ip (C.^:_1) 0 {:: getristep ^: I (rn ((((2 # [) {. ]) trsml1x (tru@{.)) ; }.) y)
)

NB. non-blocked
getrinb=: ((C.^:_1) (0 {:: (((3 : 0) ^: (# @ (1 & {::))) @ (3 : 0)))) & >/
  n=. # y
  iUL1=. trtriu y
  y=. (>/~ i. n) } iUL1 ,: y          NB. spec code
  ((0 & {.) ; ]) y
)
  'pfx sfx'=. y
  'j n'=. $ pfx
  r=. {. sfx
  l=. j {. r
  r=. (j <: (i. n)) } 0 ,: r       NB. spec code
  (pfx , (r - l mp pfx)) ; (}. sfx)
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
