NB. fri.ijs
NB. Compute the inverse using the Frobenius formula
NB.
NB. gefri  Inverse a general matrix
NB. hefri  Inverse a Hermitian (symmetric) matrix
NB. pofri  Inverse a Hermitian (symmetric) positive definite
NB.        matrix
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
NB. gefri
NB. Inverse a general matrix using the Frobenius formula
NB.
NB. Syntax:
NB.   invA=. gefri A
NB. where
NB.   A    - n×n-matrix
NB.   invA - n×n-matrix, inversion of A
NB.   n    ≥ 0
NB.
NB. Notes:
NB. - is less expensive than getri (inversion by PLU) when
NB.   7*n^(log2(7)) < 2*n^3 (i.e. when n>667)
NB.
NB. References:
NB. [1] E. E. Tyrtyshnikov "Matrix analysis and linear
NB.     algebra", Lecture notes, 2004-2005, pp. 284-285
NB.     (Тыртышников Е.Е. "Матричный анализ и линейная
NB.     алгебра", лекции, 2004-2005, стр. 284-285)
NB.     http://www.inm.ras.ru/vtm/lection/all.pdf

gefri=: 3 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    invA11=. gefri (2 $ p) {. y
    A12=. (p , (p - n)) {. y
    A21=. ((p - n) , p) {. y
    A22=. (2 $ p) }. y
    invW=. gefri (A21 mp invA11 mp A12) - A22
    M1=. invA11 mp A12
    M2=. A21 mp invA11
    M1invW=. M1 mp invW
    ((invA11 + (M1invW mp M2)) ,. M1invW) , ((invW mp M2) ,. (- invW))
  else.
    % y
  end.
)

NB. ---------------------------------------------------------
NB. hefri
NB. Inverse a Hermitian (symmetric) matrix using the
NB. Frobenius formula
NB.
NB. Syntax:
NB.   invA=. hefri A
NB. where
NB.   A    - n×n Hermitian (symmetric) matrix
NB.   invA - n×n Hermitian (symmetric) matrix, inversion of A
NB.   n    ≥ 0
NB.
NB. Notes:
NB. - strict lower triangle is not used in calculations
NB. - is less expensive than getri (inversion by PLU) when
NB.   7*n^(log2(7)) < 2*n^3 (i.e. when n>667)
NB.
NB. References:
NB. [1] E. E. Tyrtyshnikov "Matrix analysis and linear
NB.     algebra", Lecture notes, 2004-2005, pp. 284-285
NB.     (Тыртышников Е.Е. "Матричный анализ и линейная
NB.     алгебра", лекции, 2004-2005, стр. 284-285)
NB.     http://www.inm.ras.ru/vtm/lection/all.pdf

hefri=: 3 : 0
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    invA11=. hefri (2 $ p) {. y
    A12=. (p , (p - n)) {. y
    A22=. (2 $ p) }. y
    invW=. hefri ((ct A12) mp invA11 mp A12) - A22
    M1=. invA11 mp A12
    M1invW=. M1 mp invW
    ((invA11 - (M1invW mp ct M1)) ,. M1invW) , ((ct M1invW) ,. (- invW))
  else.
    % y
  end.
)

NB. =========================================================
NB. Test suite

NB. name tfri A;rcondA
tfri=: 4 : 0
  'A rcondA'=. y
  n=. # A
  I=. idmat n
  't s'=. timespacex 'invA=. ' , x , ' A'
  be=. (((norm1 (I - invA mp A)) * rcondA) % n) % FP_EPS  NB. backward error
  fe=. _.                                                 NB. forward error
  prn x ; rcondA ; be ; fe ; t ; s
)

NB. ---------------------------------------------------------
NB. ttrfri
NB. Test inverse algorithms with random triangular matrix
NB.
NB. ttrfri A

testtrfri=: 3 : 0
  n=. # y
  GE=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y

  rcondU=.  norm1 con trfriu  U=.  tru  GE
  rcondU1=. norm1 con trfriu1 U1=. tru1 GE
  rcondL=.  norm1 con trfril  L=.  trl  GE
  rcondL1=. norm1 con trfril1 L1=. trl1 GE

  'trfriu'   t1fri (U ;rcondU )
  '(128!:1)' t1fri (U ;rcondU )
  'trfriu1'  t1fri (U1;rcondU1)
  'trfril'   t1fri (L ;rcondL )
  'trfril1'  t1fri (L1;rcondL1)
  EMPTY
)

NB. ---------------------------------------------------------
NB. tgefri
NB. Test inverse algorithms with random general matrix y
NB.
NB. tgefri A

testgefri=: 3 : 0
  n=. # y
  GE=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y
  rcondGE=. norm1 con gefri   GE
  'gefri' tfri (GE;rcondGE)
  'gefri' tfri (GE;rcondGE)
  '%.'    tfri (GE;rcondGE)
  EMPTY
)

NB. ---------------------------------------------------------
NB. thefri
NB. Test inverse algorithms with random Hermitian (symmetric)
NB. matrix y
NB.
NB. thefri A

testhefri=: 3 : 0
  n=. # y
  GE=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y
  rcondHE=. norm1 con hefri HE=. ge2he GE
  'hefri' tfri (HE;rcondHE)
  EMPTY
)

NB. ---------------------------------------------------------
NB. testfri
NB. Adverb to test inverse algorithms with random matrices of
NB. generator and shape given
NB.
NB. Syntax:
NB.   mkge testfri m,n
NB. where
NB.   m,n  - 2-vector of integers, shape of random matrices
NB.          to test algorithms; if m≠n then algorithms that
NB.          accept square matrices only are skipped
NB.   mkge - monadic verb to generate random non-singular
NB.          general y-matrix (shape is taken from y)
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     (? @ $ 0:) testfri_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 16 _6 4 & gemat_mt_) testfri_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testfri_mt_ 150 200
NB.
NB. Notes:
NB. - diagonalizable matrices are processed the same way as
NB.   general matrices

testfri=: 1 : 0
  (testgefri @  u             ^: (=/)) y
  (testhefri @ (u hemat) @ {. ^: (=/)) y
  EMPTY
)
