NB. Solve linear monomial equation
NB.
NB. gesvxxx    Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a general square
NB.            matrix; op(A) is either A itself, or A^T (the
NB.            transposition of A), or A^H (the conjugate
NB.            transposition of A); B is known right-hand
NB.            sides (RHS), X is unknown solutions
NB. gtsvxxx    Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a tridiagonal
NB.            matrix; op(A) is either A itself, or A^T (the
NB.            transposition of A), or A^H (the conjugate
NB.            transposition of A); B is known right-hand
NB.            sides (RHS), X is unknown solutions
NB. hesvxxx    Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a Hermitian
NB.            (symmetric) matrix; op(A) is either A itself,
NB.            or A^T (the transposition of A); B is known
NB.            right-hand sides (RHS), X is unknown solutions
NB. posvxxx    Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a Hermitian
NB.            (symmetric) positive definite matrix; op(A) is
NB.            either A itself, or A^T (the transposition of
NB.            A); B is known right-hand sides (RHS), X is
NB.            unknown solutions
NB. ptsvxxx    Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a Hermitian
NB.            (symmetric) positive definite tridiagonal
NB.            matrix; op(A) is either A itself, or A^T (the
NB.            transposition of A); B is known right-hand
NB.            sides (RHS), X is unknown solutions
NB.
NB. testgesv1  Test gesvxxx by general square matrix and
NB.            single RHS
NB. testgesv3  Test gesvxxx by general square matrix and
NB.            multiple RHS
NB. testgtsv1  Test gtsvxxx by tridiagonal matrix and single
NB.            RHS
NB. testgtsv3  Test gtsvxxx by tridiagonal matrix and multiple
NB.            RHS
NB. testhesv1  Test hesvxxx by Hermitian (symmetric) matrix
NB.            and single RHS
NB. testhesv3  Test hesvxxx by Hermitian (symmetric) matrix
NB.            and multiple RHS
NB. testposv1  Test posvxxx by Hermitian (symmetric) positive
NB.            definite matrix and single RHS
NB. testposv3  Test posvxxx by Hermitian (symmetric) positive
NB.            definite matrix and multiple RHS
NB. testptsv1  Test ptsvxxx by Hermitian (symmetric) positive
NB.            definite tridiagonal matrix and single RHS
NB. testptsv3  Test ptsvxxx by Hermitian (symmetric) positive
NB.            definite tridiagonal matrix and multiple RHS
NB. testsv     Adv. to make verb to test xxsvxxx by matrix of
NB.            generator and shape given
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

NB. ---------------------------------------------------------
NB. gtsv
NB.
NB. Description:
NB.   Solve linear monomial equation:
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
NB.   Xv   - same shape as Bv, solutions
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
  'k k1'=. 0 1
  while. k < n1 do.
    'dk duk dlk dk1 duk1 dlk1'=. , (k , k1) { x
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
      y=. k1 -&(mul * k { y) upd y
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
      y=. (k , k1) (] v (- mul&*))/ upd y
    end.
    x=. ((dk , duk , dlk) ,: dk1 , duk1 , dlk1) (k , k1)} x
    'k k1'=. (, >:) k1
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
NB. Verb:      Solves:        Syntax:
NB. gesvax     A   * X = B    Xv=. A gesvax  Bv
NB. gesvacx    A^H * X = B    Xv=. A gesvacx Bv
NB. gesvatx    A^T * X = B    Xv=. A gesvatx Bv
NB. gesvxa     X * A   = B    Xh=. A gesvxa  Bh
NB. gesvxac    X * A^H = B    Xh=. A gesvxac Bh
NB. gesvxat    X * A^T = B    Xh=. A gesvxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general square
NB.   matrix A via triangular factorization:
NB.     P * L1 * U = A
NB. where
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, solutions
NB.   Xh   - same shape as Bh, solutions
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
NB. Verb:      Solves:        Syntax:
NB. gtsvax     A   * X = B    Xv=. A gtsvax  Bv
NB. gtsvacx    A^H * X = B    Xv=. A gtsvacx Bv
NB. gtsvatx    A^T * X = B    Xv=. A gtsvatx Bv
NB. gtsvxa     X * A   = B    Xh=. A gtsvxa  Bh
NB. gtsvxac    X * A^H = B    Xh=. A gtsvxac Bh
NB. gtsvxat    X * A^T = B    Xh=. A gtsvxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with tridiagonal
NB.   matrix A via triangular factorization:
NB.     P * L1 * U = A
NB. where
NB.   A    - n×n-matrix, the tridiagonal
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, solutions
NB.   Xh   - same shape as Bh, solutions
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

gtsvax=:  (gtsv~     0  1 _1&({"1) @(($,)~ >:@$))~
gtsvacx=: (gtsv~ +@:(0 _1  1&({"1))@(($,)~ >:@$))~
gtsvatx=: (gtsv~     0 _1  1&({"1) @(($,)~ >:@$))~

gtsvxa=:  gtsvatx                                 &.(a:`|:)
gtsvxac=: (gtsv~ +@:(0  1 _1&({"1))@(($,)~ >:@$))~&.(a:`|:)
gtsvxat=: gtsvax                                  &.(a:`|:)

NB. ---------------------------------------------------------
NB. Verb:      Solves:        Syntax:
NB. hesvax     A   * X = B    Xv=. A hesvax  Bv
NB. hesvatx    A^T * X = B    Xv=. A hesvatx Bv
NB. hesvxa     X * A   = B    Xh=. A hesvxa  Bh
NB. hesvxat    X * A^T = B    Xh=. A hesvxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) matrix A via triangular factorization:
NB.     P * L1 * T * L1^H * P^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, solutions
NB.   Xh   - same shape as Bh, solutions
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
NB. Verb:      Solves:        Syntax:
NB. posvax     A   * X = B    Xv=. A posvax  Bv
NB. posvatx    A^T * X = B    Xv=. A posvatx Bv
NB. posvxa     X * A   = B    Xh=. A posvxa  Bh
NB. posvxat    X * A^T = B    Xh=. A posvxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) positive definite matrix A via Cholesky
NB.   factorization:
NB.     L * L^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, solutions
NB.   Xh   - same shape as Bh, solutions
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
NB. Verb:      Solves:        Syntax:
NB. ptsvax     A   * X = B    Xv=. A ptsvax  Bv
NB. ptsvatx    A^T * X = B    Xv=. A ptsvatx Bv
NB. ptsvxa     X * A   = B    Xh=. A ptsvxa  Bh
NB. ptsvxat    X * A^T = B    Xh=. A ptsvxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) positive definite tridiagonal matrix A via
NB.   factorization:
NB.     L1 * D * L1^H = A
NB. where
NB.   A    - n×n-matrix, the Hermitian (symmetric) positive
NB.          definite tridiagonal
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, solutions
NB.   Xh   - same shape as Bh, solutions
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
NB. testgesv1
NB.
NB. Description:
NB.   Test:
NB.   - %. (built-in)
NB.   - gesvxxx (math/mt addon)
NB.   by general square matrix and single RHS
NB.
NB. Syntax:
NB.   testgesv1 (A ; x)
NB. where
NB.   A - n×n-matrix
NB.   x - n-vector, the exact solution

testgesv1=: 3 : 0
  'A x'=. y

  'rcondA rcondAc rcondAt'=. (gecon1 , gecon1@ct , gecon1@|:) A

  'normiAt norm1At'=. 'normiAc norm1Ac'=. 'norm1A normiA'=. (norm1 , normi) A

  ('%.'      tdyad ((1&{::)`(0&{::)`]`(3&{::)`t04v`( mp~     t02v))) A ; (    A  mp x) ; x ; rcondA  ; norm1A

  ('gesvax'  tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`( mp~     t02v))) A ; (    A  mp x) ; x ; rcondA  ; norm1A
  ('gesvacx' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp~ ct) t02v))) A ; ((ct A) mp x) ; x ; rcondAc ; norm1Ac
  ('gesvatx' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp~ |:) t02v))) A ; ((|: A) mp x) ; x ; rcondAt ; norm1At
  ('gesvxa'  tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`( mp      t02v))) A ; (x mp    A  ) ; x ; rcondA  ; normiA
  ('gesvxac' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp  ct) t02v))) A ; (x mp ct A  ) ; x ; rcondAc ; normiAc
  ('gesvxat' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp  |:) t02v))) A ; (x mp |: A  ) ; x ; rcondAt ; normiAt

  EMPTY
)

NB. ---------------------------------------------------------
NB. testgesv3
NB.
NB. Description:
NB.   Test:
NB.   - %. (built-in)
NB.   - xGESV (math/lapack2 addon)
NB.   - gesvxxx (math/mt addon)
NB.   by general square matrix and multiple RHS
NB.
NB. Syntax:
NB.   testgesv3 (A ; X)
NB. where
NB.   A - n×n-matrix
NB.   X - n×3-matrix, exact solutions

testgesv3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/gesv'

  'A Xv'=. y
  Xh=. |: Xv

  'rcondA rcondAc rcondAt'=. (gecon1 , gecon1@ct , gecon1@|:) A

  'normiAt norm1At'=. 'normiAc norm1Ac'=. 'norm1A normiA'=. (norm1 , normi) A

  Bax=. A mp Xv

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

  NB. Aapprox=. calcA (L1U ; ipiv ; trash)
  calcA=: ((C.~ makeper_jlapack2_)~ trl1pick mp trupick)~&>/@}:

  ('%.'           tdyad  ((1&{::)`(0&{::)`]`(3&{::)` vferrv       `  vberrax                                             )) A ; Bax            ; Xv ; rcondA  ; norm1A

  ('dgesv_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`(vferrv 2&{::)`((vberrax 2&{::) >. (((norm1 get01 #)~ 0 4&{)~ calcA)))) A ; Bax            ; Xv ; rcondA  ; norm1A
  ('zgesv_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`(vferrv 2&{::)`((vberrax 2&{::) >. (((norm1 get01 #)~ 0 4&{)~ calcA)))) A ; Bax            ; Xv ; rcondA  ; norm1A

  ('gesvax'       tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrv       `  vberrax                                             )) A ; Bax            ; Xv ; rcondA  ; norm1A
  ('gesvacx'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrv       `  vberracx                                            )) A ; ((ct A) mp Xv) ; Xv ; rcondAc ; norm1Ac
  ('gesvatx'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrv       `  vberratx                                            )) A ; ((|: A) mp Xv) ; Xv ; rcondAt ; norm1At
  ('gesvxa'       tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrh       `  vberrxa                                             )) A ; (Xh mp    A  ) ; Xh ; rcondA  ; normiA
  ('gesvxac'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrh       `  vberrxac                                            )) A ; (Xh mp ct A  ) ; Xh ; rcondAc ; normiAc
  ('gesvxat'      tdyad  ((0&{::)`(1&{::)`]`(3&{::)` vferrh       `  vberrxat                                            )) A ; (Xh mp |: A  ) ; Xh ; rcondAt ; normiAt

  coerase < 'mttmp'
  erase 'vferrv vferrh vberrax vberracx vberratx vberrxa vberrxac vberrxat calcA'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testgtsv1
NB.
NB. Description:
NB.   Test gtsvxxx by tridiagonal matrix and single RHS
NB.
NB. Syntax:
NB.   testgtsv1 (A ; x)
NB. where
NB.   A - n×n-matrix, the tridiagonal
NB.   x - n-vector, the exact solution

testgtsv1=: 3 : 0
  'A x'=. y

  'rcondA rcondAc rcondAt'=. (gecon1 , gecon1@ct , gecon1@|:) A

  'normiAt norm1At'=. 'normiAc norm1Ac'=. 'norm1A normiA'=. (norm1 , normi) A

  ('gtsvax'  tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`( mp~     t02v))) A ; (    A  mp x) ; x ; rcondA  ; norm1A
  ('gtsvacx' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp~ ct) t02v))) A ; ((ct A) mp x) ; x ; rcondAc ; norm1Ac
  ('gtsvatx' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp~ |:) t02v))) A ; ((|: A) mp x) ; x ; rcondAt ; norm1At
  ('gtsvxa'  tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`( mp      t02v))) A ; (x mp    A  ) ; x ; rcondA  ; normiA
  ('gtsvxac' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp  ct) t02v))) A ; (x mp ct A  ) ; x ; rcondAc ; normiAc
  ('gtsvxat' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp  |:) t02v))) A ; (x mp |: A  ) ; x ; rcondAt ; normiAt

  EMPTY
)

NB. ---------------------------------------------------------
NB. testgtsv3
NB.
NB. Description:
NB.   Test:
NB.   - xGTSV (math/lapack2 addon)
NB.   - gtsvxxx (math/mt addon)
NB.   by tridiagonal matrix and multiple RHS
NB.
NB. Syntax:
NB.   testgtsv3 (A ; X)
NB. where
NB.   A - n×n-matrix, the tridiagonal
NB.   X - n×3-matrix, exact solutions

testgtsv3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/gtsv'

  'A Xv'=. y
  Xh=. |: Xv

  'rcondA rcondAc rcondAt'=. (gecon1 , gecon1@ct , gecon1@|:) A

  'normiAt norm1At'=. 'normiAc norm1Ac'=. 'norm1A normiA'=. (norm1 , normi) A

  Bax=. A mp Xv

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

  ('dgtsv_mttmp_' tmonad (((_1&diag ; diag ; 1&diag)@(0&{::) , 1&{)`        ]`(3&{::)`vferrv`vberrax )) A ; Bax            ; Xv ; rcondA  ; norm1A
  ('zgtsv_mttmp_' tmonad (((_1&diag ; diag ; 1&diag)@(0&{::) , 1&{)`        ]`(3&{::)`vferrv`vberrax )) A ; Bax            ; Xv ; rcondA  ; norm1A

  ('gtsvax'       tdyad  ((                           0&{::       )`(1&{::)`]`(3&{::)`vferrv`vberrax )) A ; Bax            ; Xv ; rcondA  ; norm1A
  ('gtsvacx'      tdyad  ((                           0&{::       )`(1&{::)`]`(3&{::)`vferrv`vberracx)) A ; ((ct A) mp Xv) ; Xv ; rcondAc ; norm1Ac
  ('gtsvatx'      tdyad  ((                           0&{::       )`(1&{::)`]`(3&{::)`vferrv`vberratx)) A ; ((|: A) mp Xv) ; Xv ; rcondAt ; norm1At
  ('gtsvxa'       tdyad  ((                           0&{::       )`(1&{::)`]`(3&{::)`vferrh`vberrxa )) A ; (Xh mp    A  ) ; Xh ; rcondA  ; normiA
  ('gtsvxac'      tdyad  ((                           0&{::       )`(1&{::)`]`(3&{::)`vferrh`vberrxac)) A ; (Xh mp ct A  ) ; Xh ; rcondAc ; normiAc
  ('gtsvxat'      tdyad  ((                           0&{::       )`(1&{::)`]`(3&{::)`vferrh`vberrxat)) A ; (Xh mp |: A  ) ; Xh ; rcondAt ; normiAt

  coerase < 'mttmp'
  erase 'vferrv vferrh vberrax vberracx vberratx vberrxa vberrxac vberrxat'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhesv1
NB.
NB. Description:
NB.   Test hesvxxx by Hermitian (symmetric) matrix and single
NB.   RHS
NB.
NB. Syntax:
NB.   testhesv1 (A ; x)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric)
NB.   x - n-vector, the exact solution

testhesv1=: 3 : 0
  'A x'=. y

  rcondA=. hecon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  ('hesvax'  tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`( mp~     t02v))) A ; (    A  mp x) ; x ; rcondA ; norm1A
  ('hesvatx' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp~ |:) t02v))) A ; ((|: A) mp x) ; x ; rcondA ; norm1At
  ('hesvxa'  tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`( mp      t02v))) A ; (x mp    A  ) ; x ; rcondA ; normiA
  ('hesvxat' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp  |:) t02v))) A ; (x mp |: A  ) ; x ; rcondA ; normiAt

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhesv3
NB.
NB. Description:
NB.   Test:
NB.   - DSYSV DSYSV_AA ZHESV ZHESV_AA (math/lapack2 addon)
NB.   - hesvxxx (math/mt addon)
NB.   by Hermitian (symmetric) matrix and multiple RHS
NB.
NB. Syntax:
NB.   testhesv3 (A ; X)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric)
NB.   X - n×3-matrix, exact solutions
NB.
NB. Notes:
NB. - no berrA calc for LAPACK's DSYSV and ZHESV yet since
NB.   its output is intricate

testhesv3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/dsysv'
  load_mttmp_ :: ] 'math/mt/test/lapack2/dsysv_aa'
  load_mttmp_ :: ] 'math/mt/test/lapack2/zhesv'
  load_mttmp_ :: ] 'math/mt/test/lapack2/zhesv_aa'

  'A Xv'=. y
  Xh=. |: Xv

  rcondA=. hecon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  Bax=. A mp Xv

  vferrv=: normitc t04m  NB. for Xv
  vferrh=: normitr t04m  NB. for Xh

  NB. vberrX for Xv at right side
  vberrax=:  (mp~   ) t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc
  NB. vberrX for Xh at left side
  vberrxa=:  (mp    ) t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  NB. Aapprox=. calcAxx (DT1 ; ipiv ; trash)
  calcAdl=: makeper_jlapack2_@(1&{::) fp ((setdiag~  1 ;~    _1&diag )@bdlpick (mp~ mp  |:@]) trl1pick@:(|.!.0"1))@(0&{::)
  calcAdu=: makeper_jlapack2_@(1&{::) fp ((setdiag~ _1 ;~     1&diag )@bdupick (mp  mp~ |:@]) tru1pick@:(|.!.0  ))@(0&{::)
  calcAzl=: makeper_jlapack2_@(1&{::) fp ((setdiag~  1 ;~ +@(_1&diag))@bdlpick (mp~ mp  ct@]) trl1pick@:(|.!.0"1))@(0&{::)
  calcAzu=: makeper_jlapack2_@(1&{::) fp ((setdiag~ _1 ;~ +@( 1&diag))@bdupick (mp  mp~ ct@]) tru1pick@:(|.!.0  ))@(0&{::)

  ('''l''&dsysv_mttmp_'    tmonad ((2&{. )        `(2&{::)`(3&{::)` vferrv       `  vberrax                                      )) A ; Bax            ; Xv ; rcondA ; norm1A
  ('''u''&dsysv_mttmp_'    tmonad ((2&{. )        `(2&{::)`(3&{::)` vferrv       `  vberrax                                      )) A ; Bax            ; Xv ; rcondA ; norm1A
  ('''l''&zhesv_mttmp_'    tmonad ((2&{. )        `(2&{::)`(3&{::)` vferrv       `  vberrax                                      )) A ; Bax            ; Xv ; rcondA ; norm1A
  ('''u''&zhesv_mttmp_'    tmonad ((2&{. )        `(2&{::)`(3&{::)` vferrv       `  vberrax                                      )) A ; Bax            ; Xv ; rcondA ; norm1A

  ('''l''&dsysv_aa_mttmp_' tmonad ((2&{. )        `]      `(3&{::)`(vferrv 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAdl)))) A ; Bax            ; Xv ; rcondA ; norm1A
  ('''u''&dsysv_aa_mttmp_' tmonad ((2&{. )        `]      `(3&{::)`(vferrv 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAdu)))) A ; Bax            ; Xv ; rcondA ; norm1A
  ('''l''&zhesv_aa_mttmp_' tmonad ((2&{. )        `]      `(3&{::)`(vferrv 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAzl)))) A ; Bax            ; Xv ; rcondA ; norm1A
  ('''u''&zhesv_aa_mttmp_' tmonad ((2&{. )        `]      `(3&{::)`(vferrv 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAzu)))) A ; Bax            ; Xv ; rcondA ; norm1A

  ('hesvax'                tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` vferrv       `  vberrax                                      )) A ; Bax            ; Xv ; rcondA ; norm1A
  ('hesvatx'               tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` vferrv       `  vberratx                                     )) A ; ((|: A) mp Xv) ; Xv ; rcondA ; norm1At
  ('hesvxa'                tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` vferrh       `  vberrxa                                      )) A ; (Xh mp    A  ) ; Xh ; rcondA ; normiA
  ('hesvxat'               tdyad  ((0&{::)`(1&{::)`]      `(3&{::)` vferrh       `  vberrxat                                     )) A ; (Xh mp |: A  ) ; Xh ; rcondA ; normiAt

  coerase < 'mttmp'
  erase 'vferrv vferrh vberrax vberratx vberrxa vberrxat'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testposv1
NB.
NB. Description:
NB.   Test posvxxx by Hermitian (symmetric) positive definite
NB.   matrix and single RHS
NB.
NB. Syntax:
NB.   testposv1 (A ; x)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric) positive
NB.       definite
NB.   x - n-vector, the exact solution

testposv1=: 3 : 0
  'A x'=. y

  rcondA=. pocon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  ('posvax'  tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`( mp~     t02v))) A ; (    A  mp x) ; x ; rcondA ; norm1A
  ('posvatx' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp~ |:) t02v))) A ; ((|: A) mp x) ; x ; rcondA ; norm1At
  ('posvxa'  tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`( mp      t02v))) A ; (x mp    A  ) ; x ; rcondA ; normiA
  ('posvxat' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp  |:) t02v))) A ; (x mp |: A  ) ; x ; rcondA ; normiAt

  EMPTY
)

NB. ---------------------------------------------------------
NB. testposv3
NB.
NB. Description:
NB.   Test:
NB.   - xPOSV (math/lapack2 addon)
NB.   - posvxxx (math/mt addon)
NB.   by Hermitian (symmetric) positive definite matrix and
NB.   multiple RHS
NB.
NB. Syntax:
NB.   testposv3 (A ; X)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric) positive
NB.       definite
NB.   X - n×3-matrix, exact solutions

testposv3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/posv'

  'A Xv'=. y
  Xh=. |: Xv

  rcondA=. pocon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  Bax=. A mp Xv

  vferrv=: normitc t04m  NB. for Xv
  vferrh=: normitr t04m  NB. for Xh

  NB. vberrX for Xv at right side
  vberrax=:  (mp~   ) t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc
  NB. vberrX for Xh at left side
  vberrxa=:  (mp    ) t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  NB. Aapprox=. calcAxx (T ; trash)  NB. where T is either L or U
  calcAdl=: (mp  |:)@trlpick@(0&{::)
  calcAdu=: (mp~ |:)@trupick@(0&{::)
  calcAzl=: (mp  ct)@trlpick@(0&{::)
  calcAzu=: (mp~ ct)@trupick@(0&{::)

  ('''l''&dposv_mttmp_' tmonad ((2&{. )       `]`(3&{::)`(vferrv 1&{::)`((vberrax  1&{::) >. ((het01~ 0 4&{)~ calcAdl)))) A ; Bax            ; Xv ; rcondA ; norm1A
  ('''u''&dposv_mttmp_' tmonad ((2&{. )       `]`(3&{::)`(vferrv 1&{::)`((vberrax  1&{::) >. ((het01~ 0 4&{)~ calcAdu)))) A ; Bax            ; Xv ; rcondA ; norm1A
  ('''l''&zposv_mttmp_' tmonad ((2&{. )       `]`(3&{::)`(vferrv 1&{::)`((vberrax  1&{::) >. ((het01~ 0 4&{)~ calcAzl)))) A ; Bax            ; Xv ; rcondA ; norm1A
  ('''u''&zposv_mttmp_' tmonad ((2&{. )       `]`(3&{::)`(vferrv 1&{::)`((vberrax  1&{::) >. ((het01~ 0 4&{)~ calcAzu)))) A ; Bax            ; Xv ; rcondA ; norm1A

  ('posvax'             tdyad ((0&{::)`(1&{::)`]`(3&{::)` vferrv       `  vberrax                                      )) A ; Bax            ; Xv ; rcondA ; norm1A
  ('posvatx'            tdyad ((0&{::)`(1&{::)`]`(3&{::)` vferrv       `  vberratx                                     )) A ; ((|: A) mp Xv) ; Xv ; rcondA ; norm1At
  ('posvxa'             tdyad ((0&{::)`(1&{::)`]`(3&{::)` vferrh       `  vberrxa                                      )) A ; (Xh mp    A  ) ; Xh ; rcondA ; normiA
  ('posvxat'            tdyad ((0&{::)`(1&{::)`]`(3&{::)` vferrh       `  vberrxat                                     )) A ; (Xh mp |: A  ) ; Xh ; rcondA ; normiAt

  coerase < 'mttmp'
  erase 'vferrv vferrh vberrax vberratx vberrxa vberrxat calcAdl calcAdu calcAzl calcAzu'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testptsv1
NB.
NB. Description:
NB.   Test ptsvxxx by Hermitian (symmetric) positive definite
NB.   tridiagonal matrix and single RHS
NB.
NB. Syntax:
NB.   testptsv1 (A ; x)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric) positive
NB.       definite tridiagonal
NB.   x - n-vector, the exact solution

testptsv1=: 3 : 0
  'A x'=. y

  rcondA=. ptcon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  ('ptsvax'  tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`( mp~     t02v))) A ; (    A  mp x) ; x ; rcondA ; norm1A
  ('ptsvatx' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp~ |:) t02v))) A ; ((|: A) mp x) ; x ; rcondA ; norm1At
  ('ptsvxa'  tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`( mp      t02v))) A ; (x mp    A  ) ; x ; rcondA ; normiA
  ('ptsvxat' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`((mp  |:) t02v))) A ; (x mp |: A  ) ; x ; rcondA ; normiAt

  EMPTY
)

NB. ---------------------------------------------------------
NB. testptsv3
NB.
NB. Description:
NB.   Test:
NB.   - xPTSV (math/lapack2 addon)
NB.   - ptsvxxx (math/mt addon)
NB.   by Hermitian (symmetric) positive definite tridiagonal
NB.   matrix and multiple RHS
NB.
NB. Syntax:
NB.   testptsv3 (A ; X)
NB. where
NB.   A - n×n-matrix, the Hermitian (symmetric) positive
NB.       definite tridiagonal
NB.   X - n×3-matrix, exact solutions

testptsv3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/ptsv'

  'A Xv'=. y
  Xh=. |: Xv

  rcondA=. ptcon1 A

  'normiAt norm1At'=. 'norm1A normiA'=. (norm1 , normi) A

  Bax=. A mp Xv

  vferrv=: normitc t04m  NB. for Xv
  vferrh=: normitr t04m  NB. for Xh

  NB. vberrX for Xv at right side
  vberrax=:  (mp~   ) t02m norm1tc
  vberratx=: (mp~ |:) t02m norm1tc
  NB. vberrX for Xh at left side
  vberrxa=:  (mp    ) t02m norm1tr
  vberrxat=: (mp  |:) t02m norm1tr

  NB. Aapprox=. calcAxx (dD ; eT1 ; trash)
  calcAdl=: (((setdiag~ ;&_1)~ idmat@#)~ (mp mp |:@[) diagmat@[)&>/@}:
  calcAzl=: (((setdiag~ ;&_1)~ idmat@#)~ (mp mp ct@[) diagmat@[)&>/@}:

  ('dptsv_mttmp_' tmonad (((diag ; _1&diag)@(0&{::) , 1&{)`]`(3&{::)`(vferrv 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAdl)))) A ; Bax            ; Xv ; rcondA ; norm1A
  ('zptsv_mttmp_' tmonad (((diag ; _1&diag)@(0&{::) , 1&{)`]`(3&{::)`(vferrv 2&{::)`((vberrax  2&{::) >. ((het01~ 0 4&{)~ calcAzl)))) A ; Bax            ; Xv ; rcondA ; norm1A

  ('ptsvax'       tdyad  ((0&{::)`(1&{::)                 `]`(3&{::)` vferrv       `  vberrax                                      )) A ; Bax            ; Xv ; rcondA ; norm1A
  ('ptsvatx'      tdyad  ((0&{::)`(1&{::)                 `]`(3&{::)` vferrv       `  vberratx                                     )) A ; ((|: A) mp Xv) ; Xv ; rcondA ; norm1At
  ('ptsvxa'       tdyad  ((0&{::)`(1&{::)                 `]`(3&{::)` vferrv       `  vberrxa                                      )) A ; (Xh mp    A  ) ; Xh ; rcondA ; normiA
  ('ptsvxat'      tdyad  ((0&{::)`(1&{::)                 `]`(3&{::)` vferrv       `  vberrxat                                     )) A ; (Xh mp |: A  ) ; Xh ; rcondA ; normiAt

  coerase < 'mttmp'
  erase 'vferrv vferrh vberrax vberratx vberrxa vberrxat calcAdl calcAzl'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testsv
NB.
NB. Description:
NB.   Adv. to make verb to test xxsvxxx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testsv
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
NB.     ?@$&0 testsv_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testsv_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testsv_mt_ 150 150
NB.
NB. Notes:
NB. - nrhs=3 is assumed

testsv=: 1 : 'EMPTY [ ((((u ptmat2_mt_)@# (testptsv3_mt_@; [ testptsv1_mt_@(; {:"1)) ]) [ ((u pomat_mt_)@# (testposv3_mt_@; [ testposv1_mt_@(; {:"1)) ]) [ ((u hemat_mt_)@# (testhesv3_mt_@; [ testhesv1_mt_@(; {:"1)) ]))@u@({. , 3:) [ ((}."1 ((gtpick_mt_@[ (testgtsv3_mt_@; [ testgtsv1_mt_@(; {:"1)) ]) [ testgesv3_mt_@; [ testgesv1_mt_@(; {:"1)) {."1)~ _3:)@u@(+&0 3))^:(=/)'
