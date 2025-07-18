NB. Solve linear monomial equation
NB.
NB. gesvxxx   Solve the equation (op(A) * X = B) or
NB.           (X * op(A) = B), where A is a square matrix;
NB.           op(A) is either A itself, or A^T (the
NB.           transposition of A), or A^H (the conjugate
NB.           transposition of A); B is known right-hand
NB.           sides (RHS), X is unknown solutions
NB. gtsvxxx   Solve the equation (op(A) * X = B) or
NB.           (X * op(A) = B), where A is a tridiagonal
NB.           matrix; op(A) is either A itself, or A^T (the
NB.           transposition of A), or A^H (the conjugate
NB.           transposition of A); B is known right-hand
NB.           sides (RHS), X is unknown solutions
NB. hesvxxx   Solve the equation (op(A) * X = B) or
NB.           (X * op(A) = B), where A is a Hermitian
NB.           (symmetric) matrix; op(A) is either A itself,
NB.           or A^T (the transposition of A); B is known
NB.           right-hand sides (RHS), X is unknown solutions
NB. posvxxx   Solve the equation (op(A) * X = B) or
NB.           (X * op(A) = B), where A is a Hermitian
NB.           (symmetric) positive definite matrix; op(A) is
NB.           either A itself, or A^T (the transposition of
NB.           A); B is known right-hand sides (RHS), X is
NB.           unknown solutions
NB. ptsvxxx   Solve the equation (op(A) * X = B) or
NB.           (X * op(A) = B), where A is a Hermitian
NB.           (symmetric) positive definite tridiagonal
NB.           matrix; op(A) is either A itself, or A^T (the
NB.           transposition of A); B is known right-hand
NB.           sides (RHS), X is unknown solutions
NB.
NB. testgesv  Test gesvxxx by square matrix
NB. testgtsv  Test gtsvxxx by tridiagonal matrix
NB. testhesv  Test hesvxxx by Hermitian (symmetric) matrix
NB. testposv  Test posvxxx by Hermitian (symmetric) positive
NB.           definite matrix
NB. testptsv  Test ptsvxxx by Hermitian (symmetric) positive
NB.           definite tridiagonal matrix
NB. testsv    Adv. to make verb to test xxsvxxx by matrix of
NB.           generator and shape given
NB.
NB. Copyright 2010,2011,2013,2017,2018,2020,2021,2023,2024,
NB.           2025 Igor Zhuravlov
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

NB. =========================================================
NB. Configuration

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. gtsv
NB.
NB. Description:
NB.   Solves the linear monomial equation:
NB.     A * X = B
NB.   with tridiagonal matrix A via triangular factorization:
NB.     P * L1 * U = A
NB.
NB. Syntax:
NB.   Xv=. ds gtsv Bv
NB. where
NB.   ds   - (n+1)×3-matrix defined as:
NB.            ds -: (dA , trash0) ,. (duA , trash1 , trash2) ,. (dlA , trash3 , trash4)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   dA   - n-vector, diagonal of A
NB.   duA  - (n-1)-vector, superdiagonal of A
NB.   dlA  - (n-1)-vector, subdiagonal of A
NB.   P    - n×n-matrix, rows permutation of A
NB.   L1   - n×n-matrix, the unit lower triangular
NB.   U    - n×n-matrix, the upper triangular
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - gtsv implements LAPACK's xGTSV

gtsv=: 4 : 0
  v=. ,`,:@.(1 < # $ y)
  'n n1 n2 n3'=. (# x) - 1 2 3 4
  dk1=. k=. {. kk1=. 0 1  NB. (k,k1)
  while. k < n1 do.
    'dk duk dlk dk1 duk1 dlk1'=. , kk1 { x
    if. 0 = dlk do.
      NB. subdiagonal is 0, no elimination is required
      if. 0 = dk do.
        NB. diagonal is 0, a unique solution can't be found
        ($ y) $ _. return.
      end.
    elseif. dk >:&sorim dlk do.
      NB. no rows swapping required
      mul=. dlk % dk
      dk1=. dk1 - mul * duk
      y=. (- mul&*)~/\&.(kk1&{) y
      if. k < n2 do.
        dlk=. 0
      end.
    else.
      NB. swap rows k and k1
      mul=. dk % dlk
      dk=. dlk
      'duk dk1'=. duk (] , (- mul&*)) dk1
      if. k < n2 do.
        'dlk duk1'=. duk1 ([ , *) - mul
      end.
      y=. (] v (- mul&*))/&.(kk1&{) y
    end.
    x=. ((dk , duk , dlk) ,: dk1 , duk1 , dlk1) kk1} x
    k=. {. kk1=. >: kk1
  end.
  if. 0 = dk1 do.
    ($ y) $ _. return.
  end.
  NB. back solve with the matrix U from the factorization
  y=. (({: y) % dk1) n1} y
  if. 1 < n do.
    y=. (((n2 { y) - duk * {: y) % dk) n2} y
  end.
  'k k1 k2'=. n3 , n2 , n1
  while. 0 <: k do.
    'dk duk dlk'=. k { x
    y=. ((-`((+ duk&*)~ dlk&*)/ (k , k1 , k2) { y) % dk) k} y
    'k k1 k2'=. k1 (,~ (,~ <:)) k
  end.
  y
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Dyad       Solves         Syntax
NB. gesvax     A   * X = B    Xv=. A gesvax  Bv
NB. gesvacx    A^H * X = B    Xv=. A gesvacx Bv
NB. gesvatx    A^T * X = B    Xv=. A gesvatx Bv
NB. gesvxa     X * A   = B    Xh=. A gesvxa  Bh
NB. gesvxac    X * A^H = B    Xh=. A gesvxac Bh
NB. gesvxat    X * A^T = B    Xh=. A gesvxat Bh
NB.
NB. Description:
NB.   Solve the linear monomial equation with square matrix A
NB.   via triangular factorization:
NB.     P * L1 * U = A
NB. where
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   P    - n×n-matrix, rows permutation of A
NB.   L1   - n×n-matrix, the unit lower triangular
NB.   U    - n×n-matrix, the upper triangular
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   A (] -: clean@([ gesvax   mp      )) Xv
NB.   A (] -: clean@([ gesvacx (mp~ ct)~)) Xv
NB.   A (] -: clean@([ gesvatx (mp~ |:)~)) Xv
NB.   A (] -: clean@([ gesvxa   mp~     )) Xh
NB.   A (] -: clean@([ gesvxac (mp  ct)~)) Xh
NB.   A (] -: clean@([ gesvxat (mp  |:)~)) Xh
NB.
NB. Notes:
NB. - gesvax implements LAPACK's xGESV

gesvax=:  (getrspl1ux ~ getrfpl1u)~
gesvacx=: (getrspl1ucx~ getrfpl1u)~
gesvatx=: (getrspl1utx~ getrfpl1u)~
gesvxa=:  (getrsxpl1u ~ getrfpl1u)~
gesvxac=: (getrsxpl1uc~ getrfpl1u)~
gesvxat=: (getrsxpl1ut~ getrfpl1u)~

NB. ---------------------------------------------------------
NB. Dyad       Solves         Syntax
NB. gtsvax     A   * X = B    Xv=. A gtsvax  Bv
NB. gtsvacx    A^H * X = B    Xv=. A gtsvacx Bv
NB. gtsvatx    A^T * X = B    Xv=. A gtsvatx Bv
NB. gtsvxa     X * A   = B    Xh=. A gtsvxa  Bh
NB. gtsvxac    X * A^H = B    Xh=. A gtsvxac Bh
NB. gtsvxat    X * A^T = B    Xh=. A gtsvxat Bh
NB.
NB. Description:
NB.   Solve the linear monomial equation with tridiagonal
NB.   matrix A via triangular factorization:
NB.     P * L1 * U = A
NB. where
NB.   A    - n×n-matrix, the tridiagonal
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   P    - n×n-matrix, rows permutation of A
NB.   L1   - n×n-matrix, the unit lower triangular
NB.   U    - n×n-matrix, the upper triangular
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   A (] -: clean@([ gtsvax   mp      )) Xv
NB.   A (] -: clean@([ gtsvacx (mp~ ct)~)) Xv
NB.   A (] -: clean@([ gtsvatx (mp~ |:)~)) Xv
NB.   A (] -: clean@([ gtsvxa   mp~     )) Xh
NB.   A (] -: clean@([ gtsvxac (mp  ct)~)) Xh
NB.   A (] -: clean@([ gtsvxat (mp  |:)~)) Xh
NB.
NB. Notes:
NB. - gtsvax models LAPACK's xGTSV

gtsvax=:  (gtsv~     0  1 _1&({"1) @(($,)~ >:@$))~^:(*@#@])
gtsvacx=: (gtsv~ +@:(0 _1  1&({"1))@(($,)~ >:@$))~^:(*@#@])
gtsvatx=: (gtsv~     0 _1  1&({"1) @(($,)~ >:@$))~^:(*@#@])

gtsvxa=:  gtsvatx                                 &.(a:`|:)
gtsvxac=: (gtsv~ +@:(0  1 _1&({"1))@(($,)~ >:@$))~&.(a:`|:)^:(*@#@])
gtsvxat=: gtsvax                                  &.(a:`|:)

NB. ---------------------------------------------------------
NB. Dyad       Solves         Syntax
NB. hesvax     A   * X = B    Xv=. A hesvax  Bv
NB. hesvatx    A^T * X = B    Xv=. A hesvatx Bv
NB. hesvxa     X * A   = B    Xh=. A hesvxa  Bh
NB. hesvxat    X * A^T = B    Xh=. A hesvxat Bh
NB.
NB. Description:
NB.   Solve the linear monomial equation with Hermitian
NB.   (symmetric) matrix A via triangular factorization:
NB.     P * L1 * T * L1^H * P^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   P    - n×n-matrix, the full permutation of A
NB.   L1   - n×n-matrix, the unit lower triangular
NB.   T    - n×n-matrix, the Hermitian (symmetric)
NB.          tridiagonal
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   A (] -: clean@([ hesvax   mp      )) Xv
NB.   A (] -: clean@([ hesvatx (mp~ |:)~)) Xv
NB.   A (] -: clean@([ hesvxa   mp~     )) Xh
NB.   A (] -: clean@([ hesvxat (mp  |:)~)) Xh
NB.
NB. Notes:
NB. - implements LAPACK's DSYSV_AA('L'), ZHESV_AA('L')

hesvax=:  (hetrsplx ~ hetrfpl)~
hesvatx=: (hetrspltx~ hetrfpl)~
hesvxa=:  (hetrsxpl ~ hetrfpl)~
hesvxat=: (hetrsxplt~ hetrfpl)~

NB. ---------------------------------------------------------
NB. Dyad       Solves         Syntax
NB. posvax     A   * X = B    Xv=. A posvax  Bv
NB. posvatx    A^T * X = B    Xv=. A posvatx Bv
NB. posvxa     X * A   = B    Xh=. A posvxa  Bh
NB. posvxat    X * A^T = B    Xh=. A posvxat Bh
NB.
NB. Description:
NB.   Solve the linear monomial equation with Hermitian
NB.   (symmetric) positive definite matrix A via Cholesky
NB.   factorization:
NB.     L * L^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   L    - n×n-matrix, the lower triangular with positive
NB.          diagonal entries, the Cholesky triangle
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   A (] -: clean@([ posvax   mp      )) Xv
NB.   A (] -: clean@([ posvatx (mp~ |:)~)) Xv
NB.   A (] -: clean@([ posvxa   mp~     )) Xh
NB.   A (] -: clean@([ posvxat (mp  |:)~)) Xh
NB.
NB. Notes:
NB. - implements LAPACK's xPOSV('L')

posvax=:  (potrslx ~ potrfl)~
posvatx=: (potrsltx~ potrfl)~
posvxa=:  (potrsxl ~ potrfl)~
posvxat=: (potrsxlt~ potrfl)~

NB. ---------------------------------------------------------
NB. Dyad       Solves         Syntax
NB. ptsvax     A   * X = B    Xv=. A ptsvax  Bv
NB. ptsvatx    A^T * X = B    Xv=. A ptsvatx Bv
NB. ptsvxa     X * A   = B    Xh=. A ptsvxa  Bh
NB. ptsvxat    X * A^T = B    Xh=. A ptsvxat Bh
NB.
NB. Description:
NB.   Solve the linear monomial equation with Hermitian
NB.   (symmetric) positive definite tridiagonal matrix A via
NB.   factorization:
NB.     L1 * D * L1^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite tridiagonal
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - the same shape as Bv, the solution(s)
NB.   Xh   - the same shape as Bh, the solution(s)
NB.   L1   - n×n-matrix, the unit lower bidiangonal
NB.   D    - n×n-matrix, diagonal with positive diagonal
NB.          entries
NB.   n    ≥ 0, the size of A
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Assertions:
NB.   A (] -: clean@([ ptsvax   mp      )) Xv
NB.   A (] -: clean@([ ptsvatx (mp~ |:)~)) Xv
NB.   A (] -: clean@([ ptsvxa   mp~     )) Xh
NB.   A (] -: clean@([ ptsvxat (mp  |:)~)) Xh
NB.
NB. Notes:
NB. - implements LAPACK's xPTSV
NB. - if A is singular then solution Xx will be wrong
NB. - if A is indefinite then solution Xx may be wrong

ptsvax=:  (pttrslx ~ pttrfl)~
ptsvatx=: (pttrsltx~ pttrfl)~
ptsvxa=:  (pttrsxl ~ pttrfl)~
ptsvxat=: (pttrsxlt~ pttrfl)~

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgesv
NB.
NB. Description:
NB.   Test:
NB.   - %. (built-in)
NB.   - xGESV (math/lapack2 addon)
NB.   - gesvxxx (math/mt addon)
NB.   by square matrix
NB.
NB. Syntax:
NB.   log=. testgesv (X ; A)
NB. where
NB.   X   - m×n-matrix, exact solutions
NB.   A   - k×k-matrix
NB.   log - 6-vector of boxes, test log
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - models LAPACK's xDRVGE

testgesv=: 3 : 0
  _1 cocreate < 'mttmp'
  load_mttmp_ 'math/mt/external/lapack2/gesv'

  'X A'=. y
  'm n'=. $ X
  Am=. (2 # m) {. A
  An=. (2 # n) {. A

  'rcondAm rcondAcm rcondAtm'=. (gecon1 , gecon1@ct , gecon1@|:) Am
  'rcondAn rcondAcn rcondAtn'=. (gecon1 , gecon1@ct , gecon1@|:) An

  'normiAtm norm1Atm'=. 'normiAcm norm1Acm'=. 'norm1Am normiAm'=. (norm1 , normi) Am
  'normiAtn norm1Atn'=. 'normiAcn norm1Acn'=. 'norm1An normiAn'=. (norm1 , normi) An

  NB. matrix X

  Bax=. A mp X

  NB. vferrx
  vferrr=: normitc t04m  NB. for X at right side
  vferrl=: normitr t04m  NB. for X at left side

  NB. vberrX for X at right side
  vberrax=:   mp~     t02m norm1tc
  vberracx=: (mp~ ct) t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc

  NB. vberrX for X at left side
  vberrxa=:   mp      t02m norm1tr
  vberrxac=: (mp  ct) t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  NB. Aapprox=. calcA (L1U ; ipiv ; trash)
  calcA=: ((C.~ makeper_jlapack2_)~ trl1pick mp trupick)~&>/@}:

  log=.          ('%.'           tdyad  ((1&{::)`(0&{::)`]`(3&{::)` vferrr       `  vberrax                                             )) Am ; Bax           ; X  ; rcondAm  ; norm1Am

  log=. log lcat ('dgesv_mttmp_' tmonad (        (2&{. )`]`(3&{::)`(vferrr 2&{::)`((vberrax 2&{::) >. (((norm1 get01 #)~ 0 4&{)~ calcA)))) Am ; Bax           ; X  ; rcondAm  ; norm1Am
  log=. log lcat ('zgesv_mttmp_' tmonad (        (2&{. )`]`(3&{::)`(vferrr 2&{::)`((vberrax 2&{::) >. (((norm1 get01 #)~ 0 4&{)~ calcA)))) Am ; Bax           ; X  ; rcondAm  ; norm1Am

  log=. log lcat ('gesvax'       tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrr       `  vberrax                                             )) Am ; Bax           ; X  ; rcondAm  ; norm1Am
  log=. log lcat ('gesvacx'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrr       `  vberracx                                            )) Am ; (X  mp~ ct A) ; X  ; rcondAcm ; norm1Acm
  log=. log lcat ('gesvatx'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrr       `  vberratx                                            )) Am ; (X  mp~ |: A) ; X  ; rcondAtm ; norm1Atm
  log=. log lcat ('gesvxa'       tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrl       `  vberrxa                                             )) An ; (X  mp     A) ; X  ; rcondAn  ; normiAn
  log=. log lcat ('gesvxac'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrl       `  vberrxac                                            )) An ; (X  mp  ct A) ; X  ; rcondAcn ; normiAcn
  log=. log lcat ('gesvxat'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrl       `  vberrxat                                            )) An ; (X  mp  |: A) ; X  ; rcondAtn ; normiAtn

  NB. vector x

  'xm xn'=. ({."1 ; {.) X

  log=. log lcat ('%.'           tdyad  ((1&{::)`(0&{::)`]`(3&{::)` t04v         `( mp~     t02v)                                       )) Am ; (xm mp~    A) ; xm ; rcondAm  ; norm1Am

  log=. log lcat ('gesvax'       tdyad  ((0&{::)`(1&{::)`]`(3&{::)` t04v         `( mp~     t02v)                                       )) Am ; (xm mp~    A) ; xm ; rcondAm  ; norm1Am
  log=. log lcat ('gesvacx'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` t04v         `((mp~ ct) t02v)                                       )) Am ; (xm mp~ ct A) ; xm ; rcondAcm ; norm1Acm
  log=. log lcat ('gesvatx'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` t04v         `((mp~ |:) t02v)                                       )) Am ; (xm mp~ |: A) ; xm ; rcondAtm ; norm1Atm
  log=. log lcat ('gesvxa'       tdyad  ((0&{::)`(1&{::)`]`(3&{::)` t04v         `( mp      t02v)                                       )) An ; (xn mp     A) ; xn ; rcondAn  ; normiAn
  log=. log lcat ('gesvxac'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` t04v         `((mp  ct) t02v)                                       )) An ; (xn mp  ct A) ; xn ; rcondAcn ; normiAcn
  log=. log lcat ('gesvxat'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` t04v         `((mp  |:) t02v)                                       )) An ; (xn mp  |: A) ; xn ; rcondAtn ; normiAtn

  coerase < 'mttmp'
  erase 'vferrr vferrl vberrax vberracx vberratx vberrxa vberrxac vberrxat calcA'

  log
)

NB. ---------------------------------------------------------
NB. testgtsv
NB.
NB. Description:
NB.   Test:
NB.   - xGTSV (math/lapack2 addon)
NB.   - gtsvxxx (math/mt addon)
NB.   by tridiagonal matrix
NB.
NB. Syntax:
NB.   log=. testgtsv (X ; A)
NB. where
NB.   X   - m×n-matrix, exact solutions
NB.   A   - k×k-matrix, the tridiagonal
NB.   log - 6-vector of boxes, test log
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - models LAPACK's xDRVGT

testgtsv=: 3 : 0
  _1 cocreate < 'mttmp'
  load_mttmp_ 'math/mt/external/lapack2/gtsv'

  'X A'=. y
  'm n'=. $ X
  Am=. (2 # m) {. A
  An=. (2 # n) {. A

  'rcondAm rcondAcm rcondAtm'=. (gecon1 , gecon1@ct , gecon1@|:) Am
  'rcondAn rcondAcn rcondAtn'=. (gecon1 , gecon1@ct , gecon1@|:) An

  'normiAtm norm1Atm'=. 'normiAcm norm1Acm'=. 'norm1Am normiAm'=. (norm1 , normi) Am
  'normiAtn norm1Atn'=. 'normiAcn norm1Acn'=. 'norm1An normiAn'=. (norm1 , normi) An

  NB. matrix X

  Bax=. A mp X

  NB. vferrx
  vferrr=: normitc t04m  NB. for X at right side
  vferrl=: normitr t04m  NB. for X at left side

  NB. vberrX for X at right side
  vberrax=:   mp~     t02m norm1tc
  vberracx=: (mp~ ct) t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc

  NB. vberrX for X at left side
  vberrxa=:   mp      t02m norm1tr
  vberrxac=: (mp  ct) t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  log=.          ('dgtsv_mttmp_' tmonad (        ((_1&diag ; diag ; 1&diag)@(0&{::) , 1&{)`]`(3&{::)`vferrr`vberrax        )) Am ; Bax           ; X  ; rcondAm  ; norm1Am
  log=. log lcat ('zgtsv_mttmp_' tmonad (        ((_1&diag ; diag ; 1&diag)@(0&{::) , 1&{)`]`(3&{::)`vferrr`vberrax        )) Am ; Bax           ; X  ; rcondAm  ; norm1Am

  log=. log lcat ('gtsvax'       tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`vferrr`vberrax        )) Am ; Bax           ; X  ; rcondAm  ; norm1Am
  log=. log lcat ('gtsvacx'      tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`vferrr`vberracx       )) Am ; (X  mp~ ct A) ; X  ; rcondAcm ; norm1Acm
  log=. log lcat ('gtsvatx'      tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`vferrr`vberratx       )) Am ; (X  mp~ |: A) ; X  ; rcondAtm ; norm1Atm
  log=. log lcat ('gtsvxa'       tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`vferrl`vberrxa        )) An ; (X  mp     A) ; X  ; rcondAn  ; normiAn
  log=. log lcat ('gtsvxac'      tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`vferrl`vberrxac       )) An ; (X  mp  ct A) ; X  ; rcondAcn ; normiAcn
  log=. log lcat ('gtsvxat'      tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`vferrl`vberrxat       )) An ; (X  mp  |: A) ; X  ; rcondAtn ; normiAtn

  NB. vector x

  'xm xn'=. ({."1 ; {.) X

  log=. log lcat ('gtsvax'       tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`t04v  `( mp~     t02v))) Am ; (xm mp~    A) ; xm ; rcondAm  ; norm1Am
  log=. log lcat ('gtsvacx'      tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`t04v  `((mp~ ct) t02v))) Am ; (xm mp~ ct A) ; xm ; rcondAcm ; norm1Acm
  log=. log lcat ('gtsvatx'      tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`t04v  `((mp~ |:) t02v))) Am ; (xm mp~ |: A) ; xm ; rcondAtm ; norm1Atm
  log=. log lcat ('gtsvxa'       tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`t04v  `( mp      t02v))) An ; (xn mp     A) ; xn ; rcondAn  ; normiAn
  log=. log lcat ('gtsvxac'      tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`t04v  `((mp  ct) t02v))) An ; (xn mp  ct A) ; xn ; rcondAcn ; normiAcn
  log=. log lcat ('gtsvxat'      tdyad  ((0&{::)`(1&{::                                  )`]`(3&{::)`t04v  `((mp  |:) t02v))) An ; (xn mp  |: A) ; xn ; rcondAtn ; normiAtn

  coerase < 'mttmp'
  erase 'vferrr vferrl vberrax vberracx vberratx vberrxa vberrxac vberrxat'

  log
)

NB. ---------------------------------------------------------
NB. testhesv
NB.
NB. Description:
NB.   Test:
NB.   - DSYSV DSYSV_AA ZHESV ZHESV_AA (math/lapack2 addon)
NB.   - hesvxxx (math/mt addon)
NB.   by Hermitian (symmetric) matrix
NB.
NB. Syntax:
NB.   log=. testhesv (X ; A)
NB. where
NB.   X   - m×n-matrix, exact solutions
NB.   A   - k×k-matrix, the Hermitian (symmetric)
NB.   log - 6-vector of boxes, test log
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - models LAPACK's DDRVSY and ZDRVHE with the following
NB.   difference:
NB.   - no berrA calc for LAPACK's DSYSV and ZHESV yet since
NB.     its output is intricate

testhesv=: 3 : 0
  _1 cocreate < 'mttmp'
  load_mttmp_ 'math/mt/external/lapack2/dsysv'
  load_mttmp_ 'math/mt/external/lapack2/dsysv_aa'
  load_mttmp_ 'math/mt/external/lapack2/zhesv'
  load_mttmp_ 'math/mt/external/lapack2/zhesv_aa'

  'X A'=. y
  'm n'=. $ X
  Am=. (2 # m) {. A
  An=. (2 # n) {. A

  rcondAm=. hecon1 Am
  rcondAn=. hecon1 An

  'normiAtm norm1Atm'=. 'norm1Am normiAm'=. (norm1 , normi) Am
  'normiAtn norm1Atn'=. 'norm1An normiAn'=. (norm1 , normi) An

  NB. matrix X

  Bax=. A mp X

  NB. vferrx
  vferrr=: normitc t04m  NB. for X at right side
  vferrl=: normitr t04m  NB. for X at left side

  NB. vberrX for X at right side
  vberrax=:  (mp~   ) t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc

  NB. vberrX for X at left side
  vberrxa=:  (mp    ) t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  NB. Aapprox=. calcAxx (DT1 ; ipiv ; trash)
  calcAdl=: makeper_jlapack2_@(1&{::) fp ((setdiag~  1 ;~    _1&diag )@bdlpick (mp~ mp  |:@]) trl1pick@:(|.!.0"1))@(0&{::)
  calcAdu=: makeper_jlapack2_@(1&{::) fp ((setdiag~ _1 ;~     1&diag )@bdupick (mp  mp~ |:@]) tru1pick@:(|.!.0  ))@(0&{::)
  calcAzl=: makeper_jlapack2_@(1&{::) fp ((setdiag~  1 ;~ +@(_1&diag))@bdlpick (mp~ mp  ct@]) trl1pick@:(|.!.0"1))@(0&{::)
  calcAzu=: makeper_jlapack2_@(1&{::) fp ((setdiag~ _1 ;~ +@( 1&diag))@bdupick (mp  mp~ ct@]) tru1pick@:(|.!.0  ))@(0&{::)

  log=.          ('''l''&dsysv_mttmp_'    tmonad (        (2&{. )`(2&{::)`(3&{::)` vferrr       `  vberrax                                      )) Am ; Bax           ; X  ; rcondAm ; norm1Am
  log=. log lcat ('''u''&dsysv_mttmp_'    tmonad (        (2&{. )`(2&{::)`(3&{::)` vferrr       `  vberrax                                      )) Am ; Bax           ; X  ; rcondAm ; norm1Am
  log=. log lcat ('''l''&zhesv_mttmp_'    tmonad (        (2&{. )`(2&{::)`(3&{::)` vferrr       `  vberrax                                      )) Am ; Bax           ; X  ; rcondAm ; norm1Am
  log=. log lcat ('''u''&zhesv_mttmp_'    tmonad (        (2&{. )`(2&{::)`(3&{::)` vferrr       `  vberrax                                      )) Am ; Bax           ; X  ; rcondAm ; norm1Am

  log=. log lcat ('''l''&dsysv_aa_mttmp_' tmonad (        (2&{. )`]      `(3&{::)`(vferrr 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAdl)))) Am ; Bax           ; X  ; rcondAm ; norm1Am
  log=. log lcat ('''u''&dsysv_aa_mttmp_' tmonad (        (2&{. )`]      `(3&{::)`(vferrr 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAdu)))) Am ; Bax           ; X  ; rcondAm ; norm1Am
  log=. log lcat ('''l''&zhesv_aa_mttmp_' tmonad (        (2&{. )`]      `(3&{::)`(vferrr 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAzl)))) Am ; Bax           ; X  ; rcondAm ; norm1Am
  log=. log lcat ('''u''&zhesv_aa_mttmp_' tmonad (        (2&{. )`]      `(3&{::)`(vferrr 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAzu)))) Am ; Bax           ; X  ; rcondAm ; norm1Am

  log=. log lcat ('hesvax'                tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` vferrr       `  vberrax                                      )) Am ; Bax           ; X  ; rcondAm ; norm1Am
  log=. log lcat ('hesvatx'               tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` vferrr       `  vberratx                                     )) Am ; (X  mp~ |: A) ; X  ; rcondAm ; norm1Atm
  log=. log lcat ('hesvxa'                tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` vferrl       `  vberrxa                                      )) An ; (X  mp     A) ; X  ; rcondAn ; normiAn
  log=. log lcat ('hesvxat'               tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` vferrl       `  vberrxat                                     )) An ; (X  mp  |: A) ; X  ; rcondAn ; normiAtn

  NB. vector x

  'xm xn'=. ({."1 ; {.) X

  log=. log lcat ('hesvax'                tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` t04v         `( mp~     t02v                                ))) Am ; (xm mp~    A) ; xm ; rcondAm ; norm1Am
  log=. log lcat ('hesvatx'               tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` t04v         `((mp~ |:) t02v                                ))) Am ; (xm mp~ |: A) ; xm ; rcondAm ; norm1Atm
  log=. log lcat ('hesvxa'                tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` t04v         `( mp      t02v                                ))) An ; (xn mp     A) ; xn ; rcondAn ; normiAn
  log=. log lcat ('hesvxat'               tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` t04v         `((mp  |:) t02v                                ))) An ; (xn mp  |: A) ; xn ; rcondAn ; normiAtn

  coerase < 'mttmp'
  erase 'vferrr vferrl vberrax vberratx vberrxa vberrxat calcAdl calcAdu calcAzl calcAzu'

  log
)

NB. ---------------------------------------------------------
NB. testposv
NB.
NB. Description:
NB.   Test:
NB.   - xPOSV (math/lapack2 addon)
NB.   - posvxxx (math/mt addon)
NB.   by Hermitian (symmetric) positive definite matrix
NB.
NB. Syntax:
NB.   log=. testposv (X ; A)
NB. where
NB.   X   - m×n-matrix, exact solutions
NB.   A   - k×k-matrix, the Hermitian (symmetric) positive
NB.         definite
NB.   log - 6-vector of boxes, test log
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - models LAPACK's xDRVPO

testposv=: 3 : 0
  _1 cocreate < 'mttmp'
  load_mttmp_ 'math/mt/external/lapack2/posv'

  'X A'=. y
  'm n'=. $ X
  Am=. (2 # m) {. A
  An=. (2 # n) {. A

  rcondAm=. pocon1 Am
  rcondAn=. pocon1 An

  'normiAtm norm1Atm'=. 'norm1Am normiAm'=. (norm1 , normi) Am
  'normiAtn norm1Atn'=. 'norm1An normiAn'=. (norm1 , normi) An

  NB. matrix X

  Bax=. A mp X

  NB. vferrx
  vferrr=: normitc t04m  NB. for X at right side
  vferrl=: normitr t04m  NB. for X at left side

  NB. vberrX for X at right side
  vberrax=:  (mp~   ) t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc

  NB. vberrX for X at left side
  vberrxa=:  (mp    ) t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  NB. Aapprox=. calcAxx (T ; trash)  NB. where T is either L or U
  calcAdl=: (mp  |:)@trlpick@(0&{::)
  calcAdu=: (mp~ |:)@trupick@(0&{::)
  calcAzl=: (mp  ct)@trlpick@(0&{::)
  calcAzu=: (mp~ ct)@trupick@(0&{::)

  log=.          ('''l''&dposv_mttmp_' tmonad (        (2&{. )`]`(3&{::)`(vferrr 1&{::)`((vberrax  1&{::) >. ((het01~ 0 4&{)~ calcAdl)))) Am ; Bax          ; X  ; rcondAm ; norm1Am
  log=. log lcat ('''u''&dposv_mttmp_' tmonad (        (2&{. )`]`(3&{::)`(vferrr 1&{::)`((vberrax  1&{::) >. ((het01~ 0 4&{)~ calcAdu)))) Am ; Bax          ; X  ; rcondAm ; norm1Am
  log=. log lcat ('''l''&zposv_mttmp_' tmonad (        (2&{. )`]`(3&{::)`(vferrr 1&{::)`((vberrax  1&{::) >. ((het01~ 0 4&{)~ calcAzl)))) Am ; Bax          ; X  ; rcondAm ; norm1Am
  log=. log lcat ('''u''&zposv_mttmp_' tmonad (        (2&{. )`]`(3&{::)`(vferrr 1&{::)`((vberrax  1&{::) >. ((het01~ 0 4&{)~ calcAzu)))) Am ; Bax          ; X  ; rcondAm ; norm1Am

  log=. log lcat ('posvax'             tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrr       `  vberrax                                      )) Am ; Bax          ; X  ; rcondAm ; norm1Am
  log=. log lcat ('posvatx'            tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrr       `  vberratx                                     )) Am ; (X mp~ |: A) ; X  ; rcondAm ; norm1Atm
  log=. log lcat ('posvxa'             tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrl       `  vberrxa                                      )) An ; (X mp     A) ; X  ; rcondAn ; normiAn
  log=. log lcat ('posvxat'            tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrl       `  vberrxat                                     )) An ; (X mp  |: A) ; X  ; rcondAn ; normiAtn

  NB. vector x

  'xm xn'=. ({."1 ; {.) X

  log=. log lcat ('posvax'             tdyad  ((0&{::)`(1&{::)`]`(3&{::)` t04v         `( mp~     t02v                                ))) Am ; (xm mp~    A) ; xm ; rcondAm ; norm1Am
  log=. log lcat ('posvatx'            tdyad  ((0&{::)`(1&{::)`]`(3&{::)` t04v         `((mp~ |:) t02v                                ))) Am ; (xm mp~ |: A) ; xm ; rcondAm ; norm1Atm
  log=. log lcat ('posvxa'             tdyad  ((0&{::)`(1&{::)`]`(3&{::)` t04v         `( mp      t02v                                ))) An ; (xn mp     A) ; xn ; rcondAn ; normiAn
  log=. log lcat ('posvxat'            tdyad  ((0&{::)`(1&{::)`]`(3&{::)` t04v         `((mp  |:) t02v                                ))) An ; (xn mp  |: A) ; xn ; rcondAn ; normiAtn

  coerase < 'mttmp'
  erase 'vferrr vferrl vberrax vberratx vberrxa vberrxat calcAdl calcAdu calcAzl calcAzu'

  log
)

NB. ---------------------------------------------------------
NB. testptsv
NB.
NB. Description:
NB.   Test:
NB.   - xPTSV (math/lapack2 addon)
NB.   - ptsvxxx (math/mt addon)
NB.   by Hermitian (symmetric) positive definite tridiagonal
NB.   matrix
NB.
NB. Syntax:
NB.   log=. testptsv (X ; A)
NB. where
NB.   X   - m×n-matrix, exact solutions
NB.   A   - k×k-matrix, the Hermitian (symmetric) positive
NB.         definite tridiagonal
NB.   log - 6-vector of boxes, test log
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - models LAPACK's xDRVPT

testptsv=: 3 : 0
  _1 cocreate < 'mttmp'
  load_mttmp_ 'math/mt/external/lapack2/ptsv'

  'X A'=. y
  'm n'=. $ X
  Am=. (2 # m) {. A
  An=. (2 # n) {. A

  rcondAm=. ptcon1 Am
  rcondAn=. ptcon1 An

  'normiAtm norm1Atm'=. 'norm1Am normiAm'=. (norm1 , normi) Am
  'normiAtn norm1Atn'=. 'norm1An normiAn'=. (norm1 , normi) An

  NB. matrix X

  Bax=. A mp X

  NB. vferrx
  vferrr=: normitc t04m  NB. for X at right side
  vferrl=: normitr t04m  NB. for X at left side

  NB. vberrX for X at right side
  vberrax=:  (mp~   ) t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc

  NB. vberrX for X at left side
  vberrxa=:  (mp    ) t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  NB. Aapprox=. calcAxx (dD ; eT1 ; trash)
  calcAdl=: (((setdiag~ ;&_1)~ idmat@#)~ (mp mp |:@[) diagmat@[)&>/@}:
  calcAzl=: (((setdiag~ ;&_1)~ idmat@#)~ (mp mp ct@[) diagmat@[)&>/@}:

  log=.          ('dptsv_mttmp_' tmonad (        ((diag ; _1&diag)@(0&{::) , 1&{)`]`(3&{::)`(vferrr 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAdl)))) Am ; Bax          ; X  ; rcondAm ; norm1Am
  log=. log lcat ('zptsv_mttmp_' tmonad (        ((diag ; _1&diag)@(0&{::) , 1&{)`]`(3&{::)`(vferrr 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAzl)))) Am ; Bax          ; X  ; rcondAm ; norm1Am

  log=. log lcat ('ptsvax'       tdyad  ((0&{::)`(1&{::                         )`]`(3&{::)` vferrr       `  vberrax                                      )) Am ; Bax          ; X  ; rcondAm ; norm1Am
  log=. log lcat ('ptsvatx'      tdyad  ((0&{::)`(1&{::                         )`]`(3&{::)` vferrr       `  vberratx                                     )) Am ; (X mp~ |: A) ; X  ; rcondAm ; norm1Atm
  log=. log lcat ('ptsvxa'       tdyad  ((0&{::)`(1&{::                         )`]`(3&{::)` vferrr       `  vberrxa                                      )) An ; (X mp     A) ; X  ; rcondAn ; normiAn
  log=. log lcat ('ptsvxat'      tdyad  ((0&{::)`(1&{::                         )`]`(3&{::)` vferrr       `  vberrxat                                     )) An ; (X mp  |: A) ; X  ; rcondAn ; normiAtn

  NB. vector x

  'xm xn'=. ({."1 ; {.) X

  log=. log lcat ('ptsvax'       tdyad  ((0&{::)`(1&{::                         )`]`(3&{::)` t04v         `( mp~     t02v                                ))) Am ; (xm mp~    A) ; xm ; rcondAm ; norm1Am
  log=. log lcat ('ptsvatx'      tdyad  ((0&{::)`(1&{::                         )`]`(3&{::)` t04v         `((mp~ |:) t02v                                ))) Am ; (xm mp~ |: A) ; xm ; rcondAm ; norm1Atm
  log=. log lcat ('ptsvxa'       tdyad  ((0&{::)`(1&{::                         )`]`(3&{::)` t04v         `( mp      t02v                                ))) An ; (xn mp     A) ; xn ; rcondAn ; normiAn
  log=. log lcat ('ptsvxat'      tdyad  ((0&{::)`(1&{::                         )`]`(3&{::)` t04v         `((mp  |:) t02v                                ))) An ; (xn mp  |: A) ; xn ; rcondAn ; normiAtn

  coerase < 'mttmp'
  erase 'vferrr vferrl vberrax vberratx vberrxa vberrxat calcAdl calcAzl'

  log
)

NB. ---------------------------------------------------------
NB. testsv
NB.
NB. Description:
NB.   Adv. to make verb to test xxsvxxx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testsv) (m,n)
NB. where
NB.   mkmat - monad to generate a material for matrix; is
NB.           called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of mat
NB.   log   - 6-vector of boxes, test log
NB.
NB. Application:
NB. - test by random real matrices with elements distributed
NB.   uniformly with support (0,1), system is either
NB.   100×100-matrix or 150×150-matrix and RHS is either
NB.   100-vector or 150-vector or 100×150-matrix:
NB.     log=. ?@$&0 testsv_mt_ 100 150
NB. - test by random real matrices with elements with
NB.   limited value's amplitude, system is 150×150-matrix and
NB.   RHS is either 150-vector or 150×150-matrix:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testsv_mt_ 150 150
NB. - test by random complex matrices, system is either
NB.   200×200-matrix or 150×150-matrix and RHS is either
NB.   200-vector or 150-vector or 200×150-matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testsv_mt_ 200 150

testsv=: 1 : 'lcat_mt_@(testgesv_mt_@:(<@u"1)`(testgtsv_mt_@(u@{. ; gtpick_mt_@u@{:))`(testhesv_mt_@(u@{. ; (u hemat_mt_)@{:))`(testposv_mt_@(u@{. ; (u pomat_mt_)@{:))`(testptsv_mt_@(u@{. ; (u ptmat2_mt_)@{:))`:0)@(,: >./)'
