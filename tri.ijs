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
NB. testtrtri  Test triangular inversion algorithms by
NB.            triangular matrix given
NB. testgetri  Test triangular inversion algorithms by
NB.            general matrix given
NB. testhetri  Test triangular inversion algorithms by
NB.            Hermitian (symmetric) matrix given
NB. testpotri  Test triangular inversion algorithms by
NB.            Hermitian (symmetric) positive definite
NB.            matrix given
NB. testpttri  Test triangular inversion algorithms by
NB.            Hermitian (symmetric) positive definite
NB.            tridiagonal matrix given
NB. testtri    Adv. to make verb to test triangular
NB.            inversion algorithms by matrix of generator
NB.            and shape given
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. trtri
NB. Template adverb to make triangular matrix inversion verbs
NB.
NB. Syntax:
NB.   u0=. u0`u1`u2`u3`u4`u5`u6 trtri
NB. where
NB.   u0 - inversion verb for recursive call
NB.   u1 - extract the square triangular matrix Ta from
NB.        either top left or bottom right corner
NB.   u2 - extract the square triangular matrix Tb opposite
NB.        to Ta on diagonal
NB.   u3 - combine intervals to extract non-zero part either
NB.        from top right part (tru case), or bottom left
NB.        part (trl case) of input matrix A
NB.   u4 - combine inverses of both invTa and invAc
NB.   u5 - assemble output as triangular matrix
NB.   u6 - inverse 1×1 sized system
NB.
NB. References:
NB. [1] E. E. Tyrtyshnikov "Matrix analysis and linear
NB.     algebra", Lecture notes, 2004-2005, pp. 284-285
NB.     (Тыртышников Е.Е. "Матричный анализ и линейная
NB.     алгебра", лекции, 2004-2005, стр. 284-285)
NB.     http://www.inm.ras.ru/vtm/lection/all.pdf

trtri=: 1 : 0
  '`u0 u1 u2 u3 u4 u5 u6'=. u
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    Ta=. (2 $ p) u1 y
    invTa=. u0 Ta
    invTb=. u0 (2 $ p) u2 y
    Ac=. (p u3 (p - n)) {. y
    invAc=. - invTa mp Ac mp invTb
    (invTa u4 invAc) u5 invTb
  else.
    u6 y
  end.
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:          Solves:              Syntax:
NB. trtriu         U  * invU  = I       invU=. trtriu A
NB. trtriu1        U1 * invU1 = I       invU1=. trtriu1 A
NB. trtril         L  * invL  = I       invL=. trtril A
NB. trtril1        L1 * invL1 = I       invL1=. trtril1 A
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

trtriu=:  trtriu `{.`}.` ,  ` ,.  `(_1 append) `%      trtri
trtriu1=: trtriu1`{.`}.` ,  ` ,.  `(_1 append) `(1:"0) trtri
trtril=:  trtril `}.`{.`(,~)`(,.~)`( 0 append~)`%      trtri
trtril1=: trtril1`}.`{.`(,~)`(,.~)`( 0 append~)`(1:"0) trtri

NB. ---------------------------------------------------------
NB. getri
NB. Inverse a general matrix using the LU factorization
NB.   inv(P) * L1 * U = A
NB.   inv(A) = inv(U) * inv(L1) * P

getri=: ((C."1~ /:)~ (trtriu mp trtril1)) & >/ @ getrfpl1u

NB. ##############
NB.    ('';1 1 0 1 1) <;.1 i.3 5
NB. ┌──┬─────┬──┬──┐
NB. │ 0│ 1  2│ 3│ 4│
NB. │ 5│ 6  7│ 8│ 9│
NB. │10│11 12│13│14│
NB. └──┴─────┴──┴──┘

TRNB=: 3 NB. 64

getri2step=: 3 : 0
  'pfx bcol sfx'=. y
  j=. c pfx
  jb=. c bcol
  L0=. ((j+i.jb);(i.jb)) { bcol
  L1=. ((j+jb-n),jb) {. bcol
  bcol=. bcol - sfx mp L1
  bcol=. L0 trsmxl1 bcol
  pfx=. NB }."1 pfx
  bcol=. (-TRNB) {."1 pfx
  sfx=. bcol ,. sfx
  pfx ; bcol ; sfx
)

getri2=: 3 : 0
  'ip L1U'=. y
  n=. # L1U
  nn=. (>: dbg 'nn') TRNB * <. (n-1) % TRNB
  I=. (<. dbg 'I') (n+(TRNB-1)) % TRNB
  jb=. TRNB (<. dbg 'jb') (>: n-nn)
  ios=. (jb+nn-1) (ht2lios dbg 'ht2lios') nn-1
  ip C.^:_1"1 (,. &: >)/ getri2step ^: I (nn {."1 L1U) ; (ios ({"1 dbg '{"1') L1U) ; ((nn+jb-n) {."1 L1U)
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

NB. name ttri A;rcondA
testtri=: 4 : 0
  'A rcondA'=. y
  n=. # A
  I=. idmat n
  't s'=. timespacex 'invA=. ' , x , ' A'
  be=. (((norm1 (I - invA mp A)) * rcondA) % n) % FP_EPS  NB. backward error
  fe=. _.                                                 NB. forward error
  prn x ; rcondA ; be ; fe ; t ; s
)

NB. ---------------------------------------------------------
NB. testtrtri
NB. Test inverse algorithms with random triangular matrix
NB.
NB. ttrtri A

testtrtri=: 3 : 0
  n=. # y
  GE=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y

  rcondU=.  norm1 con trtriu  U=.  tru  GE
  rcondU1=. norm1 con trtriu1 U1=. tru1 GE
  rcondL=.  norm1 con trtril  L=.  trl  GE
  rcondL1=. norm1 con trtril1 L1=. trl1 GE

  'trtriu'   t1tri (U ;rcondU )
  '(128!:1)' t1tri (U ;rcondU )
  'trtriu1'  t1tri (U1;rcondU1)
  'trtril'   t1tri (L ;rcondL )
  'trtril1'  t1tri (L1;rcondL1)
  EMPTY
)

NB. ---------------------------------------------------------
NB. testgetri
NB. Test inverse algorithms with random general matrix y
NB.
NB. tgetri A

testgetri=: 3 : 0
  n=. # y
  GE=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y
  rcondGE=. norm1 con getri   GE
  'getri' ttri (GE;rcondGE)
  'gefri' ttri (GE;rcondGE)
  '%.'    ttri (GE;rcondGE)
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
NB. - test by random square real matrix with limited values'
NB.   amplitudes:
NB.     (_1 1 0 16 _6 4 & gemat_mt_) testtri_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testtri_mt_ 150 200

testtri=: 1 : 'EMPTY_mt_ [ ((testpttri_mt_ @ (u ptmat_mt_)) [ (testpotri_mt_ @ (u pomat_mt_)) [ (testhetri_mt_ @ (u hemat_mt_)) [ (testgetri_mt_ @ u) [ (testtrtri_mt_ @ (u trumat_mt_))) ^: (=/)'
