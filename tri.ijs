NB. tri.ijs
NB. Compute the inverse of a matrix using the triangular
NB. factorization
NB.
NB. trtriu   Inverse an upper triangular matrix
NB. trtriu1  Inverse an unit upper triangular matrix
NB. trtril   Inverse a lower triangular matrix
NB. trtril1  Inverse an unit lower triangular matrix
NB. getri    Inverse a general matrix
NB. hetri    Inverse a Hermitian (symmetric) matrix
NB. ditri    Inverse a diagonalizable matrix
NB. potri    Inverse a Hermitian (symmetric) positive
NB.          definite matrix
NB.
NB. Version: 1.0.0 2009-06-01
NB. Copyright: Igor Zhuravlov, igor at uic.dvgu.ru
NB. License: GNU GPL

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. trtri
NB. Template adverb to make triangular matrix inversion verbs
NB.
NB. Syntax:
NB.   vtrtri=. vtrtri`u1`u2`u3`u4`u5`u6 trtri
NB. where
NB.   vtrtri - inversion verb for recursive call
NB.   u1     - either {. for [unit] upper triangular matrix
NB.            or }. for [unit] lower triangular matrix
NB.   u2     - either }. for [unit] upper triangular matrix
NB.            or {. for [unit] lower triangular matrix
NB.   u3     - either , for [unit] upper triangular matrix
NB.            or ,~ for [unit] lower triangular matrix
NB.   u4     - either trtrsux for upper triangular matrix, or
NB.            trtrsu1x for unit upper triangular matrix, or
NB.            trtrslx for lower triangular matrix, or
NB.            trtrsu1x for unit lower triangular matrix
NB.   u5     - either ,. for [unit] upper triangular matrix
NB.            or ,.~ for [unit] lower triangular matrix
NB.   u6     - either (_1 append) for [unit] upper triangular
NB.            matrix or (0 append~) for [unit] lower
NB.            triangular matrix

trtri=: 1 : 0
  '`u0 u1 u2 u3 u4 u5 u6'=. u
  n=. # y
  if. n > 1 do.
    p=. >. -: n
    Ta=. (2 $ p) u1 y              NB. square matrix from top left or bottom right corner
    invTa=. u0 Ta                  NB. inverse it recursively
    invTb=. u0 (2 $ p) u2 y        NB. inverse the square matrix from opposite corner recursively
    Ac=. (p u3 (p - n)) {. y       NB. non-zero off-diagonal part of input matrix
    invAc=. Ta u4 (- Ac mp invTb)  NB. invert it
    (invTa u5 invAc) u6 invTb      NB. assemble output as triangular matrix
  else.
    % y
  end.
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. trtriu
NB. trtriu1
NB. trtril
NB. trtril1
NB. Various invertors for triangular matrix
NB.
NB. Syntax:
NB.   invT=. invertor A
NB. where
NB.   invertor - any verb of: trtriu trtriu1 trtril trtril1
NB.   A        - n×n-matrix, containing triangular matrix T
NB.   invT     - n×n-matrix, inversion of T
NB.   n      ≥ 0
NB.
NB. If:
NB.   invT=. <<<<<<<<<<<<
NB.   X=. 
NB. then
NB.   ?

erase 'trtriu';'trtriu1';'trtril';'trtril1'

trtriu=:  trtriu `{.`}.` ,  `trtrsux ` ,.  `(_1 append)  trtri  NB. inverse      upper triangular matrix
trtriu1=: trtriu1`{.`}.` ,  `trtrsu1x` ,.  `(_1 append)  trtri  NB. inverse unit upper triangular matrix
trtril=:  trtril `}.`{.`(,~)`trtrslx `(,.~)`( 0 append~) trtri  NB. inverse      lower triangular matrix
trtril1=: trtril1`}.`{.`(,~)`trtrsl1x`(,.~)`( 0 append~) trtri  NB. inverse unit lower triangular matrix

NB. ---------------------------------------------------------
NB. getri
NB. Inverse a general matrix using the LU factorization
NB.   inv(P) * L1 * U = A
NB.   inv(A) = inv(U) * inv(L1) * P

getri=: ((C."1~ /:)~ ((trtriu@tru) mp (trtril1@trl1))) & >/ @ getrfl1u

NB. ---------------------------------------------------------
NB. potri
NB. Inverse a Hermitian (symmetric) positive definite matrix
NB. using the Cholesky factorization:
NB.   L * L' = A
NB.   inv(A) = inv(L') * inv(L)

potri=: (mp~ ct) @ trtril @ potrfl

NB. =========================================================
NB. Test suite

NB. 'name rcond bwerror fwerror time space'=. name t1tri A;rcondA
t1tri=: 4 : 0
  'A rcondA'=. y
  n=. # A
  I=. idmat n
  't s'=. timespacex 'invA=. ' , x , ' A'
  be=. (((norm1 (I - invA mp A)) * rcondA) % n) % FP_EPS  NB. backward error
  fe=. i.0                                                NB. forward error
  dbprn x ; rcondA ; be ; fe ; t ; s
)

NB. ---------------------------------------------------------
NB. ttri
NB. Test inverse algorithms with random matrix y and its
NB. derivatives
NB.
NB. r=. ttri A
NB. where A is a general square matrix
NB.       r is boxed table with columns: name rcond bwerror fwerror time space
ttri=: 3 : 0
  n=. # y
  GE=. (sdiag~ (# $ ((10&*)@:((*@diag) * (>./@:|@,))))) y

  rcondU=.  norm1 con trtriu  U=.  tru  GE
  rcondU1=. norm1 con trtriu1 U1=. tru1 GE
  rcondL=.  norm1 con trtril  L=.  trl  GE
  rcondL1=. norm1 con trtril1 L1=. trl1 GE
  rcondGE=. norm1 con getri   GE
NB.  rcondHE=. norm1 con hetri   HE=. ge2he GE
NB.  rcondDI=. norm1 con ditri   DI=. ge2di GE
  rcondPO=. norm1 con potri   PO=. (mp ct) GE  NB. >>> ###### PO=. ge2po GE

  r=.  ,: 'trtriu'   t1tri (U ;rcondU )
  r=. r , '(128!:1)' t1tri (U ;rcondU )
  r=. r , 'trtriu1'  t1tri (U1;rcondU1)
  r=. r , 'trtril'   t1tri (L ;rcondL )
  r=. r , 'trtril1'  t1tri (L1;rcondL1)
  r=. r , 'getri'    t1tri (GE;rcondGE)
  r=. r , '%.'       t1tri (GE;rcondGE)
NB.  r=. r , 'hotri'    t1tri (HE;rcondHE)
NB.  r=. r , 'ditri'    t1tri (DI;rcondDI)
  r=. r , 'potri'    t1tri (PO;rcondPO)
)

NB. ---------------------------------------------------------
NB. testtri
NB. Adverb to test inverse algorithms with random matrices of
NB. shape y
NB.
NB. Syntax:
NB.   r=. mkge testtri m,n
NB. where
NB.   m,n  - 2-vector of integers, shape of random matrices
NB.          to test algorithms; if m≠n then algorithms that
NB.          accept square matrices only are skipped
NB.   mkge - monadic verb to generate random non-singular
NB.          general y-matrix (shape is taken from y)
NB.   r    - boxed table with 6 columns:
NB.          - algorithm name
NB.          - the estimated reciprocal of the condition
NB.            number of a matrix in 1-norm
NB.          - relative backward error
NB.          - relative forward error
NB.          - execution time, sec.
NB.          - execution space, bytes
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   cocurrent 'mt'
NB.   r=. (_1 1 0 16 _6 4 & gemat) testtri 500 500
NB.   r=. (_1 1 0 16 _6 4 & (gemat j. gemat)) testtri 500 500

testtri=: 1 : 'ttri @ u ^: (=/)'
