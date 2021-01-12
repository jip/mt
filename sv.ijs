NB. Solve linear monomial equation
NB.
NB. gesvxxx    Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a general matrix;
NB.            op(A) is either A itself, or A^T (the
NB.            transposition of A), or A^H (the conjugate
NB.            transposition of A); B is known right-hand
NB.            side (RHS), X is unknown solution
NB. hesvxxx    Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a Hermitian
NB.            (symmetric) matrix; op(A) is either A itself,
NB.            or A^T (the transposition of A); B is known
NB.            right-hand side (RHS), X is unknown solution
NB. posvxxx    Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a Hermitian
NB.            (symmetric) positive definite matrix; op(A) is
NB.            either A itself, or A^T (the transposition of
NB.            A); B is known right-hand side (RHS), X is
NB.            unknown solution
NB. ptsvxxx    Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is a Hermitian
NB.            (symmetric) positive definite tridiagonal
NB.            matrix; op(A) is either A itself, or A^T (the
NB.            transposition of A); B is known right-hand
NB.            side (RHS), X is unknown solution
NB.
NB. testgesv1  Test gesvxxx by square matrix and single RHS
NB. testgesv3  Test gesvxxx by square matrix and multiple RHS
NB. testhesv   Test hesvxxx by Hermitian (symmetric) matrix
NB. testposv   Test posvxxx by Hermitian (symmetric) positive
NB.            definite matrix
NB. testptsv   Test ptsvxxx by Hermitian (symmetric) positive
NB.            definite tridiagonal matrix
NB. testsv     Adv. to make verb to test xxsvxxx by matrix of
NB.            generator and shape given
NB.
NB. Version: 0.10.0 2017-04-23
NB.
NB. Copyright 2010-2017 Igor Zhuravlov
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
NB. Verb:      Solves:        Syntax:
NB. gesvax     A   * X = B    Xv=. A gesvax  Bv
NB. gesvacx    A^H * X = B    Xv=. A gesvacx Bv
NB. gesvatx    A^T * X = B    Xv=. A gesvatx Bv
NB. gesvxa     X * A   = B    Xh=. A gesvxa  Bh
NB. gesvxac    X * A^H = B    Xh=. A gesvxac Bh
NB. gesvxat    X * A^T = B    Xh=. A gesvxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with general matrix A
NB.   via triangular factorization:
NB.     P * L1 * U = A
NB. where:
NB.   A    - n×n-matrix
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   P    - n×n-matrix, rows permutation of A
NB.   L1   - n×n-matrix, unit lower triangular
NB.   U    - n×n-matrix, upper triangular
NB.   nrhs ≥ 0
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
NB. hesvax     A   * X = B    Xv=. A hesvax  Bv
NB. hesvatx    A^T * X = B    Xv=. A hesvatx Bv
NB. hesvxa     X * A   = B    Xh=. A hesvxa  Bh
NB. hesvxat    X * A^T = B    Xh=. A hesvxat Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with Hermitian
NB.   (symmetric) matrix A via triangular factorization:
NB.     P * L1 * T * L1^H * P^H = A
NB. where:
NB.   A    - n×n-matrix, Hermitian (symmetric)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   P    - n×n-matrix, full permutation of A
NB.   L1   - n×n-matrix, unit lower triangular
NB.   T    - n×n-matrix, Hermitian (symmetric) tridiagonal
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - implements LAPACK's DSYSV('L'), ZHESV('L')

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
NB. where:
NB.   A    - n×n-matrix, Hermitian (symmetric) positive
NB.          definite
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   L    - n×n-matrix, lower triangular with positive
NB.          diagonal entries, Cholesky triangle
NB.   nrhs ≥ 0
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
NB. where:
NB.   A    - n×n-matrix, Hermitian (symmetric) positive
NB.          definite tridiagonal
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   L1   - n×n-matrix, unit lower bidiangonal
NB.   D    - n×n-matrix, diagonal with positive diagonal
NB.          entries
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - implements LAPACK's xPTSV

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
NB.   by square matrix and single RHS
NB.
NB. Syntax:
NB.   testgesv3 (A ; x)
NB. where
NB.   A - n×n-matrix
NB.   xX - n-vector, exact solution
NB.
NB. Formula:
NB. - see testgesv3

testgesv1=: 3 : 0
  'A x'=. y

  'rcondA rcondAc rcondAt'=. (gecon1 , gecon1@ct , gecon1@|:) A
  'norm1A normiA'=. (norm1 , normi) A

  NB. ferr=. (A ; B ; X ; rcondA ; normA) vferr Xapprox
  vferr_mttmp_=. %&FP_EPS^:(1 > FP_EPS&*)@((3 {:: [) *`%/@,`(%@FP_EPS)@.(>/@:*@]) ] (- ,&normitr_mt_ ]) 2 {:: [)`(%@FP_EPS)@.(0 = 3 {:: [)`0:@.(0 = #@])

  NB. berrX=. (A ; B ; X ; rcondA ; normA) (calcB aberrX) Xapprox
  aberrX=. 1 : '((FP_EPS , 4 {:: [) %~/@,@,.`(%@FP_EPS)@.(0 ([ = {) ]) ] ,&norm1tc_mt_ (u 0&{::)~ - 1 {:: [)`(%@FP_EPS)@.(0 = 4 {:: [)`0:@.(0 = #@])'
  NB. vberr for x at right side
  vberrax_mttmp_=.   mp_mt_~         aberrX
  vberracx_mttmp_=. (mp_mt_~ ct_mt_) aberrX
  vberratx_mttmp_=. (mp_mt_~ |:    ) aberrX
  NB. vberr for x at left side
  vberrxa_mttmp_=.   mp_mt_          aberrX
  vberrxac_mttmp_=. (mp_mt_  ct_mt_) aberrX
  vberrxat_mttmp_=. (mp_mt_  |:    ) aberrX

  ('%.'      tdyad  (1&{::`(0&{::)`]`(rcondA "_)`vferr_mttmp_`vberrax_mttmp_ )) A ; (    A  mp x) ; x ; rcondA  ; norm1A

  ('gesvax'  tdyad  (0&{::`(1&{::)`]`(rcondA "_)`vferr_mttmp_`vberrax_mttmp_ )) A ; (    A  mp x) ; x ; rcondA  ; norm1A
  ('gesvacx' tdyad  (0&{::`(1&{::)`]`(rcondAc"_)`vferr_mttmp_`vberracx_mttmp_)) A ; ((ct A) mp x) ; x ; rcondAc ; norm1A
  ('gesvatx' tdyad  (0&{::`(1&{::)`]`(rcondAt"_)`vferr_mttmp_`vberratx_mttmp_)) A ; ((|: A) mp x) ; x ; rcondAt ; norm1A
  ('gesvxa'  tdyad  (0&{::`(1&{::)`]`(rcondA "_)`vferr_mttmp_`vberrxa_mttmp_ )) A ; (x mp    A  ) ; x ; rcondA  ; norm1A
  ('gesvxac' tdyad  (0&{::`(1&{::)`]`(rcondAc"_)`vferr_mttmp_`vberrxac_mttmp_)) A ; (x mp ct A  ) ; x ; rcondAc ; norm1A
  ('gesvxat' tdyad  (0&{::`(1&{::)`]`(rcondAt"_)`vferr_mttmp_`vberrxat_mttmp_)) A ; (x mp |: A  ) ; x ; rcondAt ; norm1A

  coerase < 'mttmp'

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
NB.   by square matrix and multiple RHS
NB.
NB. Syntax:
NB.   testgesv3 (A ; X)
NB. where
NB.   A - n×n-matrix
NB.   X - n×3-matrix, exact solutions
NB.
NB. Formula:
NB. - ferr for ax case:
NB.   foreach i-th pair (X,Xapprox) from nrhs solutions do
NB.     if 0=n or 0=nrhs then
NB.       ferr[i] := 0
NB.     elseif 0=rcond(op(A)) or 0=||X|| and 0<||X-Xapprox|| then
NB.       ferr[i] := 1 / FP_EPS
NB.     else
NB.       ferr[i] := (||X - Xapprox|| / ||X||) * rcond(op(A))
NB.     endif
NB.   endfor
NB.   ferr := max(ferr[i])
NB.   if ferr * FP_EPS < 1 then
NB.     ferr := ferr / FP_EPS
NB.   endif
NB.   ||vector|| := normit(vector)
NB. - berr for ax case:
NB.   berr := max(berrA,berrX)
NB.   berrA := ((||P * L1 * U - A||_1 / n) / ||A||_1) / FP_EPS
NB.   foreach i-th computed solution X from nrhs solutions do
NB.     if 0=m or 0=n or 0=nrhs then
NB.       berrX[i] := 0
NB.     elseif 0=||op(A)||_1 or 0=min(||X||) then
NB.       berrX[i] := 1 / FP_EPS
NB.     else
NB.       berrX[i] := ((||B - op(A) * X|| / ||op(A)||_1) / ||X||) / FP_EPS
NB.     endif
NB.   endfor
NB.   berrX := max(berrX[i])
NB.   ||vector|| := norm1t(vector)
NB.
NB. Notes:
NB. - models LAPACK's xGET01, xGET02 and xGET04

testgesv3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/gesv'

  'A Xv'=. y
  Xh=. |: Xv
  'rcondA rcondAc rcondAt'=. (gecon1 , gecon1@ct , gecon1@|:) A
  'norm1A normiA'=. (norm1 , normi) A
  Bax=. A mp Xv

  NB. ferr=. (A ; B ; X ; rcondA ; normA) (normitx aferr) Xapprox
  aferr=. 1 : '%&FP_EPS^:(1 > FP_EPS&*)@max_mt_@((3 {:: [) >/@:*@}.`(*`%/ ,: %@FP_EPS)}@, ] (- ,:&u ]) 2 {:: [)`(%@FP_EPS)@.(0 = 3 {:: [)`0:@.(0 e. $@])'
  vferrv_mttmp_=. normitc_mt_ aferr  NB. for Xv
  vferrh_mttmp_=. normitr_mt_ aferr  NB. for Xh

  NB. berrX=. (A ; B ; X ; rcondA ; normA) (calcB cberrX norm1tx) Xapprox
  cberrX=. 2 : 'max_mt_@((FP_EPS , 4 {:: [) %~/@,@,.`(%@FP_EPS)@.(0 ([ = {) ])"1 ] ,.&v (u 0&{::)~ - 1 {:: [)`(%@FP_EPS)@.(0 = 4 {:: [)`0:@.(0 e. $@])'
  NB. vberr for Xv at right side
  vberrax_mttmp_=.   mp_mt_~         cberrX norm1tc_mt_
  vberracx_mttmp_=. (mp_mt_~ ct_mt_) cberrX norm1tc_mt_
  vberratx_mttmp_=. (mp_mt_~ |:    ) cberrX norm1tc_mt_
  NB. vberr for Xh at left side
  vberrxa_mttmp_=.   mp_mt_          cberrX norm1tr_mt_
  vberrxac_mttmp_=. (mp_mt_  ct_mt_) cberrX norm1tr_mt_
  vberrxat_mttmp_=. (mp_mt_  |:    ) cberrX norm1tr_mt_

  NB. Aapprox := P * L1 * U
  NB. Aapprox=. calcA_mttmp_ (L1U ; ipiv ; trash)
  calcA_mttmp_=. ((C.~ makeper_jlapack2_)~ trl1pick mp trupick)~&>/@}:

  NB. berrA=. (A ; B ; X ; rcondA ; normA) vberrA_mttmp_ Aapprox
  vberrA_mttmp_=. ((FP_EPS , #@]) %~/@,@,.`(%@FP_EPS)@.(</@:*@]) ((] ,&norm1_mt_ -) 0&{::)~)`0:@.(0 = #@])

  ('%.'           tdyad  (1&{::`(0&{::)`]`(rcondA "_)` vferrv_mttmp_       `  vberrax_mttmp_                                        )) A ; Bax            ; Xv ; rcondA  ; norm1A

  ('dgesv_mttmp_' tmonad (2&{. `        ]`(rcondA "_)`(vferrv_mttmp_ 2&{::)`((vberrax_mttmp_ 2&{::) >. (vberrA_mttmp_ calcA_mttmp_)))) A ; Bax            ; Xv ; rcondA  ; norm1A
  ('zgesv_mttmp_' tmonad (2&{. `        ]`(rcondA "_)`(vferrv_mttmp_ 2&{::)`((vberrax_mttmp_ 2&{::) >. (vberrA_mttmp_ calcA_mttmp_)))) A ; Bax            ; Xv ; rcondA  ; norm1A

  ('gesvax'       tdyad  (0&{::`(1&{::)`]`(rcondA "_)` vferrv_mttmp_       `  vberrax_mttmp_                                        )) A ; Bax            ; Xv ; rcondA  ; norm1A
  ('gesvacx'      tdyad  (0&{::`(1&{::)`]`(rcondAc"_)` vferrv_mttmp_       `  vberracx_mttmp_                                       )) A ; ((ct A) mp Xv) ; Xv ; rcondAc ; norm1A
  ('gesvatx'      tdyad  (0&{::`(1&{::)`]`(rcondAt"_)` vferrv_mttmp_       `  vberratx_mttmp_                                       )) A ; ((|: A) mp Xv) ; Xv ; rcondAt ; norm1A
  ('gesvxa'       tdyad  (0&{::`(1&{::)`]`(rcondA "_)` vferrh_mttmp_       `  vberrxa_mttmp_                                        )) A ; (Xh mp    A  ) ; Xh ; rcondA  ; norm1A
  ('gesvxac'      tdyad  (0&{::`(1&{::)`]`(rcondAc"_)` vferrh_mttmp_       `  vberrxac_mttmp_                                       )) A ; (Xh mp ct A  ) ; Xh ; rcondAc ; norm1A
  ('gesvxat'      tdyad  (0&{::`(1&{::)`]`(rcondAt"_)` vferrh_mttmp_       `  vberrxat_mttmp_                                       )) A ; (Xh mp |: A  ) ; Xh ; rcondAt ; norm1A

  coerase < 'mttmp'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhesv
NB.
NB. Description:
NB.   Test hesvxxx by Hermitian (symmetric) matrix
NB.
NB. Syntax:
NB.   testhesv (A;X)
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric)
NB.   X - n×n-matrix, exact solutions
NB.
NB. Formula:
NB.   ferr := max(||X - Xapprox|| / ||Xapprox||)
NB.   berr := max(||B - op(A) * X|| / (FP_EPS * ||op(A)|| * ||X||))

testhesv=: 3 : 0
  'A X'=. y
  'rcondA rcondAt'=. (hecon1 , hecon1@|:) A

  ('hesvax'  tdyad ((    0&{:: )`(mp &>/)`]`(rcondA "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(mp &>/@[ - (mp~ 0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tc@])))) y
  ('hesvatx' tdyad ((|:@(0&{::))`(mp &>/)`]`(rcondAt"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(mp &>/@[ - (mp~ 0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tc@])))) (|: A) ; X
  ('hesvxa'  tdyad ((    0&{:: )`(mp~&>/)`]`(rcondA "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(mp~&>/@[ - (mp  0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tr@])))) y
  ('hesvxat' tdyad ((|:@(0&{::))`(mp~&>/)`]`(rcondAt"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(mp~&>/@[ - (mp  0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tr@])))) (|: A) ; X

  EMPTY
)

NB. ---------------------------------------------------------
NB. testposv
NB.
NB. Description:
NB.   Test posvxxx by Hermitian (symmetric) positive definite
NB.   matrix
NB.
NB. Syntax:
NB.   testposv (A;X)
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.   X - n×n-matrix, exact solutions
NB.
NB. Formula:
NB.   ferr := max(||X - Xapprox|| / ||Xapprox||)
NB.   berr := max(||B - op(A) * X|| / (FP_EPS * ||op(A)|| * ||X||))

testposv=: 3 : 0
  'A X'=. y
  'rcondA rcondAt'=. (pocon1 , pocon1@|:) A

  ('posvax'  tdyad ((    0&{:: )`(mp &>/)`]`(rcondA "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(mp &>/@[ - (mp~ 0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tc@])))) y
  ('posvatx' tdyad ((|:@(0&{::))`(mp &>/)`]`(rcondAt"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(mp &>/@[ - (mp~ 0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tc@])))) (|: A) ; X
  ('posvxa'  tdyad ((    0&{:: )`(mp~&>/)`]`(rcondA "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(mp~&>/@[ - (mp  0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tr@])))) y
  ('posvxat' tdyad ((|:@(0&{::))`(mp~&>/)`]`(rcondAt"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(mp~&>/@[ - (mp  0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tr@])))) (|: A) ; X

  EMPTY
)

NB. ---------------------------------------------------------
NB. testptsv
NB.
NB. Description:
NB.   Test ptsvxxx by Hermitian (symmetric) positive definite
NB.   tridiagonal matrix
NB.
NB. Syntax:
NB.   testptsv (A;X)
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.       tridiagonal
NB.   X - n×n-matrix, exact solutions
NB.
NB. Formula:
NB.   ferr := max(||X - Xapprox|| / ||Xapprox||)
NB.   berr := max(||B - op(A) * X|| / (FP_EPS * ||op(A)|| * ||X||))

testptsv=: 3 : 0
  'A X'=. y
  'rcondA rcondAt'=. (ptcon1 , ptcon1@|:) A

  ('ptsvax'  tdyad ((    0&{:: )`(mp &>/)`]`(rcondA "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(mp &>/@[ - (mp~ 0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tc@])))) y
  ('ptsvatx' tdyad ((|:@(0&{::))`(mp &>/)`]`(rcondAt"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(mp &>/@[ - (mp~ 0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tc@])))) (|: A) ; X
  ('ptsvxa'  tdyad ((    0&{:: )`(mp~&>/)`]`(rcondA "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(mp~&>/@[ - (mp  0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tr@])))) y
  ('ptsvxat' tdyad ((|:@(0&{::))`(mp~&>/)`]`(rcondAt"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(mp~&>/@[ - (mp  0&{::)~) % (FP_EPS * 1:^:(0&=)@norm1@(0 {:: [)) * norm1tr@])))) (|: A) ; X

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

testsv=: 1 : 'EMPTY [ (testptsv_mt_@((u ptmat2_mt_) ; u) [ testposv_mt_@((u pomat_mt_) ; u) [ testhesv_mt_@((u hemat_mt_) ; u) [ (testgesv3_mt_@((}."1 ; {."1)~ _3:) [ testgesv1_mt_@(_3&(}."1) ; {:"1))@u@(+&0 3))^:(=/)'
