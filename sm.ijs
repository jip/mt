NB. Solve linear monomial equation with triangular matrix
NB.
NB. trsmxxxx   Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is either unit or
NB.            non-unit, either lower or upper, triangular
NB.            matrix; op(A) is either A itself, or A^T, the
NB.            transposition of A, or A^H, the conjugate
NB.            transposition of A; B is known right-hand side
NB.            (RHS), X is unknown solution
NB.
NB. testtrsm1  Test trsmxxxx by triangular matrix and single
NB.            RHS
NB. testtrsm3  Test trsmxxxx by triangular matrix and
NB.            multiple RHS
NB. testsm     Adv. to make verb to test trsmxxxx by matrix of
NB.            generator and shape given
NB.
NB. Version: 0.10.5 2020-03-30
NB.
NB. Copyright 2010-2020 Igor Zhuravlov
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
NB. Verb:       Solves:         Syntax:
NB. trsmllcn    L^H  * X = B    Xv=. A trsmllcn Bv
NB. trsmllcu    L1^H * X = B    Xv=. A trsmllcu Bv
NB. trsmllnn    L    * X = B    Xv=. A trsmllnn Bv
NB. trsmllnu    L1   * X = B    Xv=. A trsmllnu Bv
NB. trsmlltn    L^T  * X = B    Xv=. A trsmlltn Bv
NB. trsmlltu    L1^T * X = B    Xv=. A trsmlltu Bv
NB. trsmlucn    U^H  * X = B    Xv=. A trsmlucn Bv
NB. trsmlucu    U1^H * X = B    Xv=. A trsmlucu Bv
NB. trsmlunn    U    * X = B    Xv=. A trsmlunn Bv
NB. trsmlunu    U1   * X = B    Xv=. A trsmlunu Bv
NB. trsmlutn    U^T  * X = B    Xv=. A trsmlutn Bv
NB. trsmlutu    U1^T * X = B    Xv=. A trsmlutu Bv
NB. trsmrlcn    X * L^H  = B    Xh=. A trsmrlcn Bh
NB. trsmrlcu    X * L1^H = B    Xh=. A trsmrlcu Bh
NB. trsmrlnn    X * L    = B    Xh=. A trsmrlnn Bh
NB. trsmrlnu    X * L1   = B    Xh=. A trsmrlnu Bh
NB. trsmrltn    X * L^T  = B    Xh=. A trsmrltn Bh
NB. trsmrltu    X * L1^T = B    Xh=. A trsmrltu Bh
NB. trsmrucn    X * U^H  = B    Xh=. A trsmrucn Bh
NB. trsmrucu    X * U1^H = B    Xh=. A trsmrucu Bh
NB. trsmrunn    X * U    = B    Xh=. A trsmrunn Bh
NB. trsmrunu    X * U1   = B    Xh=. A trsmrunu Bh
NB. trsmrutn    X * U^T  = B    Xh=. A trsmrutn Bh
NB. trsmrutu    X * U1^T = B    Xh=. A trsmrutu Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with triangular matrix
NB. where:
NB.   A    - n×n-matrix, containing either U, U1, L or L1
NB.   U    - n×n-matrix, upper triangular
NB.   U1   - n×n-matrix, unit upper triangular (diagonal is
NB.          not saved)
NB.   L    - n×n-matrix, lower triangular
NB.   L1   - n×n-matrix, unit lower triangular (diagonal is
NB.          not saved)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - opposite triangle is not referenced
NB. - unit diagonal is not referenced
NB. - models BLAS's xTRSV for vector RHS
NB. - models BLAS's xTRSM for matrix RHS with the following
NB.   difference: alpha parameter is assumed to be always
NB.   equal to 1; to workaround this limitation use the
NB.   following pattern:
NB.     X=. A trsmxxxx alpha*B
NB.
NB. TODO:
NB. - replace column-wise algos by row-wise

trsmllcn=: (((((_1 - #@]) {    1 {:: [ ) - ] mp~ ((_1 lisoS)&# 0&{::)~ +@({,) 0 {:: [) % ((_1 - (*>:))&# 0&{::)~ +@({,) 0 {:: [) ,      ])^:(;`(#@])`(0 {.   ]))  NB. liso(li)=(n dhs2liso (-1-i  ,i)), lio(lii)=-1-i*(n+1)
trsmllcu=: ( (((_1 - #@]) {    1 {:: [ ) - ] mp~ ((_1 lisoS)&# 0&{::)~ +@({,) 0 {:: [                                          ) ,      ])^:(;`(#@])`(0 {.   ]))  NB. liso(li)=(n dhs2liso (-1-i  ,i))
trsmllnn=: ((((      #@]  {    1 {:: [ ) - ] mp~ (( 1 lisoW)&# 0&{::)~   ({,) 0 {:: [) % (      (*>:) &# 0&{::)~   ({,) 0 {:: [) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. liso(li)=(1 dhs2liso (   i*n,i)), lio(lii)=   i*(n+1)
trsmllnu=: ( ((      #@]  {    1 {:: [ ) - ] mp~ (( 1 lisoW)&# 0&{::)~   ({,) 0 {:: [                                          ) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. liso(li)=(1 dhs2liso (   i*n,i))
trsmlltn=: (((((_1 - #@]) {    1 {:: [ ) - ] mp~ ((_1 lisoS)&# 0&{::)~   ({,) 0 {:: [) % ((_1 - (*>:))&# 0&{::)~   ({,) 0 {:: [) ,      ])^:(;`(#@])`(0 {.   ]))  NB. liso(li)=(n dhs2liso (-1-i  ,i)), lio(lii)=-1-i*(n+1)
trsmlltu=: ( (((_1 - #@]) {    1 {:: [ ) - ] mp~ ((_1 lisoS)&# 0&{::)~   ({,) 0 {:: [                                          ) ,      ])^:(;`(#@])`(0 {.   ]))  NB. liso(li)=(n dhs2liso (-1-i  ,i))

trsmlucn=: ((((      #@]  {    1 {:: [ ) - ] mp~ (( 1 lisoN)&# 0&{::)~ +@({,) 0 {:: [) % (      (*>:) &# 0&{::)~ +@({,) 0 {:: [) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. liso(ui)=(n dhs2liso (   i  ,i)), lio(lii)=   i*(n+1)
trsmlucu=: ( ((      #@]  {    1 {:: [ ) - ] mp~ (( 1 lisoN)&# 0&{::)~ +@({,) 0 {:: [                                          ) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. liso(ui)=(n dhs2liso (   i  ,i))
trsmlunn=: (((((_1 - #@]) {    1 {:: [ ) - ] mp~ ((_1 lisoE)&# 0&{::)~   ({,) 0 {:: [) % ((_1 - (*>:))&# 0&{::)~   ({,) 0 {:: [) ,      ])^:(;`(#@])`(0 {.   ]))  NB. liso(ui)=(1 dhs2liso (-1-i*n,i)), lio(lii)=-1-i*(n+1)
trsmlunu=: ( (((_1 - #@]) {    1 {:: [ ) - ] mp~ ((_1 lisoE)&# 0&{::)~   ({,) 0 {:: [                                          ) ,      ])^:(;`(#@])`(0 {.   ]))  NB. liso(ui)=(1 dhs2liso (-1-i*n,i))
trsmlutn=: ((((      #@]  {    1 {:: [ ) - ] mp~ (( 1 lisoN)&# 0&{::)~   ({,) 0 {:: [) % (      (*>:) &# 0&{::)~   ({,) 0 {:: [) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. liso(ui)=(n dhs2liso (   i  ,i)), lio(lii)=   i*(n+1)
trsmlutu=: ( ((      #@]  {    1 {:: [ ) - ] mp~ (( 1 lisoN)&# 0&{::)~   ({,) 0 {:: [                                          ) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. liso(ui)=(n dhs2liso (   i  ,i))

trsmrlcn=: ((((      c@]  {"1 (1 {:: [)) - ] mp  (( 1 lisoW)&c 0&{::)~ +@({,) 0 {:: [) % (      (*>:) &c 0&{::)~ +@({,) 0 {:: [) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(1 dhs2liso (   i*n,i)), lio(lii)=   i*(n+1)
trsmrlcu=: ( ((      c@]  {"1 (1 {:: [)) - ] mp  (( 1 lisoW)&c 0&{::)~ +@({,) 0 {:: [                                          ) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(1 dhs2liso (   i*n,i))
trsmrlnn=: (((((_1 - c@]) {"1 (1 {:: [)) - ] mp  ((_1 lisoS)&c 0&{::)~   ({,) 0 {:: [) % ((_1 - (*>:))&c 0&{::)~   ({,) 0 {:: [) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(n dhs2liso (-1-i  ,i)), lio(lii)=-1-i*(n+1)
trsmrlnu=: ( (((_1 - c@]) {"1 (1 {:: [)) - ] mp  ((_1 lisoS)&c 0&{::)~   ({,) 0 {:: [                                          ) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(n dhs2liso (-1-i  ,i))
trsmrltn=: ((((      c@]  {"1 (1 {:: [)) - ] mp  (( 1 lisoW)&c 0&{::)~   ({,) 0 {:: [) % (      (*>:) &c 0&{::)~   ({,) 0 {:: [) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(1 dhs2liso (   i*n,i)), lio(lii)=   i*(n+1)
trsmrltu=: ( ((      c@]  {"1 (1 {:: [)) - ] mp  (( 1 lisoW)&c 0&{::)~   ({,) 0 {:: [                                          ) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(1 dhs2liso (   i*n,i))

trsmrucn=: (((((_1 - c@]) {"1 (1 {:: [)) - ] mp  ((_1 lisoE)&c 0&{::)~ +@({,) 0 {:: [) % ((_1 - (*>:))&c 0&{::)~ +@({,) 0 {:: [) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(1 dhs2liso (-1-i*n,i)), lio(lii)=-1-i*(n+1)
trsmrucu=: ( (((_1 - c@]) {"1 (1 {:: [)) - ] mp  ((_1 lisoE)&c 0&{::)~ +@({,) 0 {:: [                                          ) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(1 dhs2liso (-1-i*n,i))
trsmrunn=: ((((      c@]  {"1 (1 {:: [)) - ] mp  (( 1 lisoN)&c 0&{::)~   ({,) 0 {:: [) % (      (*>:) &c 0&{::)~   ({,) 0 {:: [) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(n dhs2liso (   i  ,i)), lio(lii)=   i*(n+1)
trsmrunu=: ( ((      c@]  {"1 (1 {:: [)) - ] mp  (( 1 lisoN)&c 0&{::)~   ({,) 0 {:: [                                          ) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(n dhs2liso (   i  ,i))
trsmrutn=: (((((_1 - c@]) {"1 (1 {:: [)) - ] mp  ((_1 lisoE)&c 0&{::)~   ({,) 0 {:: [) % ((_1 - (*>:))&c 0&{::)~   ({,) 0 {:: [) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(1 dhs2liso (-1-i*n,i)), lio(lii)=-1-i*(n+1)
trsmrutu=: ( (((_1 - c@]) {"1 (1 {:: [)) - ] mp  ((_1 lisoE)&c 0&{::)~   ({,) 0 {:: [                                          ) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. liso(li)=(1 dhs2liso (-1-i*n,i))

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testtrsm1
NB.
NB. Description:
NB.   Test trsmxxxx (math/mt addon) by triangular matrix and
NB.   single RHS
NB.
NB. Syntax:
NB.   testtrsm1 (A ; x)
NB. where
NB.   A - n×n-matrix, triangular
NB.   x - n-vector, exact solution
NB.
NB. Formula:
NB. - see testtrsm3

testtrsm1=: 3 : 0
  'A x'=. y

  rcondL=.  trlcon1  L=.  trlpick  A
  rcondU=.  trucon1  U=.  trupick  A
  rcondL1=. trl1con1 L1=. (1 ; '') setdiag L
  rcondU1=. tru1con1 U1=. (1 ; '') setdiag U

  'normAxlnn normAxlxn'=. (norm1 , normi) L
  'normAxlnu normAxlxu'=. (norm1 , normi) L1
  'normAxunn normAxuxn'=. (norm1 , normi) U
  'normAxunu normAxuxu'=. (norm1 , normi) U1

  NB. ferr=. (A ; B ; X ; rcondA ; normA) vferr Xapprox
  vferr_mttmp_=. %&FP_EPS^:(1 > FP_EPS&*)@((3 {:: [) *`%/@,`(%@FP_EPS)@.(>/@:*@]) ] (- ,&normitr_mt_ ]) 2 {:: [)`(%@FP_EPS)@.(0 = 3 {:: [)`0:@.(0 = #@])

  NB. berr=. (A ; B ; X ; rcondA ; normA) (calcB aberr) Xapprox
  aberr=. 1 : '((FP_EPS , 4 {:: [) %~/@,@,.`(%@FP_EPS)@.(0 ([ = {) ]) ] ,&norm1tc_mt_ (u 0&{::)~ - 1 {:: [)`(%@FP_EPS)@.(0 = 4 {:: [)`0:@.(0 = #@])'
  NB. vberr for x at right side
  vberrlxcx_mttmp_=. (mp_mt_~ ct_mt_) aberr
  vberrlxnx_mttmp_=.  mp_mt_~         aberr
  vberrlxtx_mttmp_=. (mp_mt_~ |:    ) aberr
  NB. vberr for x at left side
  vberrrxcx_mttmp_=. (mp_mt_  ct_mt_) aberr
  vberrrxnx_mttmp_=.  mp_mt_          aberr
  vberrrxtx_mttmp_=. (mp_mt_  |:    ) aberr

  ('trsmllcn' tdyad (0&{::`(1&{::)`]`(rcondL "_)`vferr_mttmp_`vberrlxcx_mttmp_)) L  ; (L  (mp~ ct)~ x ) ; x  ; rcondL  ; normAxlxn
  ('trsmllcu' tdyad (0&{::`(1&{::)`]`(rcondL1"_)`vferr_mttmp_`vberrlxcx_mttmp_)) L1 ; (L1 (mp~ ct)~ x ) ; x  ; rcondL1 ; normAxlxu
  ('trsmllnn' tdyad (0&{::`(1&{::)`]`(rcondL "_)`vferr_mttmp_`vberrlxnx_mttmp_)) L  ; (L   mp       x ) ; x  ; rcondL  ; normAxlnn
  ('trsmllnu' tdyad (0&{::`(1&{::)`]`(rcondL1"_)`vferr_mttmp_`vberrlxnx_mttmp_)) L1 ; (L1  mp       x ) ; x  ; rcondL1 ; normAxlnu
  ('trsmlltn' tdyad (0&{::`(1&{::)`]`(rcondL "_)`vferr_mttmp_`vberrlxtx_mttmp_)) L  ; (L  (mp~ |:)~ x ) ; x  ; rcondL  ; normAxlxn
  ('trsmlltu' tdyad (0&{::`(1&{::)`]`(rcondL1"_)`vferr_mttmp_`vberrlxtx_mttmp_)) L1 ; (L1 (mp~ |:)~ x ) ; x  ; rcondL1 ; normAxlxu

  ('trsmlucn' tdyad (0&{::`(1&{::)`]`(rcondU "_)`vferr_mttmp_`vberrlxcx_mttmp_)) U  ; (U  (mp~ ct)~ x ) ; x  ; rcondU  ; normAxuxn
  ('trsmlucu' tdyad (0&{::`(1&{::)`]`(rcondU1"_)`vferr_mttmp_`vberrlxcx_mttmp_)) U1 ; (U1 (mp~ ct)~ x ) ; x  ; rcondU1 ; normAxuxu
  ('trsmlunn' tdyad (0&{::`(1&{::)`]`(rcondU "_)`vferr_mttmp_`vberrlxnx_mttmp_)) U  ; (U   mp       x ) ; x  ; rcondU  ; normAxunn
  ('trsmlunu' tdyad (0&{::`(1&{::)`]`(rcondU1"_)`vferr_mttmp_`vberrlxnx_mttmp_)) U1 ; (U1  mp       x ) ; x  ; rcondU1 ; normAxunu
  ('trsmlutn' tdyad (0&{::`(1&{::)`]`(rcondU "_)`vferr_mttmp_`vberrlxtx_mttmp_)) U  ; (U  (mp~ |:)~ x ) ; x  ; rcondU  ; normAxuxn
  ('trsmlutu' tdyad (0&{::`(1&{::)`]`(rcondU1"_)`vferr_mttmp_`vberrlxtx_mttmp_)) U1 ; (U1 (mp~ |:)~ x ) ; x  ; rcondU1 ; normAxuxu

  ('trsmrlcn' tdyad (0&{::`(1&{::)`]`(rcondL "_)`vferr_mttmp_`vberrrxcx_mttmp_)) L  ; (x  (mp  ct)  L ) ; x  ; rcondL  ; normAxlxn
  ('trsmrlcu' tdyad (0&{::`(1&{::)`]`(rcondL1"_)`vferr_mttmp_`vberrrxcx_mttmp_)) L1 ; (x  (mp  ct)  L1) ; x  ; rcondL1 ; normAxlxu
  ('trsmrlnn' tdyad (0&{::`(1&{::)`]`(rcondL "_)`vferr_mttmp_`vberrrxnx_mttmp_)) L  ; (x   mp       L ) ; x  ; rcondL  ; normAxlnn
  ('trsmrlnu' tdyad (0&{::`(1&{::)`]`(rcondL1"_)`vferr_mttmp_`vberrrxnx_mttmp_)) L1 ; (x   mp       L1) ; x  ; rcondL1 ; normAxlnu
  ('trsmrltn' tdyad (0&{::`(1&{::)`]`(rcondL "_)`vferr_mttmp_`vberrrxtx_mttmp_)) L  ; (x  (mp  |:)  L ) ; x  ; rcondL  ; normAxlxn
  ('trsmrltu' tdyad (0&{::`(1&{::)`]`(rcondL1"_)`vferr_mttmp_`vberrrxtx_mttmp_)) L1 ; (x  (mp  |:)  L1) ; x  ; rcondL1 ; normAxlxu

  ('trsmrucn' tdyad (0&{::`(1&{::)`]`(rcondU "_)`vferr_mttmp_`vberrrxcx_mttmp_)) U  ; (x  (mp  ct)  U ) ; x  ; rcondU  ; normAxuxn
  ('trsmrucu' tdyad (0&{::`(1&{::)`]`(rcondU1"_)`vferr_mttmp_`vberrrxcx_mttmp_)) U1 ; (x  (mp  ct)  U1) ; x  ; rcondU1 ; normAxuxu
  ('trsmrunn' tdyad (0&{::`(1&{::)`]`(rcondU "_)`vferr_mttmp_`vberrrxnx_mttmp_)) U  ; (x   mp       U ) ; x  ; rcondU  ; normAxunn
  ('trsmrunu' tdyad (0&{::`(1&{::)`]`(rcondU1"_)`vferr_mttmp_`vberrrxnx_mttmp_)) U1 ; (x   mp       U1) ; x  ; rcondU1 ; normAxunu
  ('trsmrutn' tdyad (0&{::`(1&{::)`]`(rcondU "_)`vferr_mttmp_`vberrrxtx_mttmp_)) U  ; (x  (mp  |:)  U ) ; x  ; rcondU  ; normAxuxn
  ('trsmrutu' tdyad (0&{::`(1&{::)`]`(rcondU1"_)`vferr_mttmp_`vberrrxtx_mttmp_)) U1 ; (x  (mp  |:)  U1) ; x  ; rcondU1 ; normAxuxu

  coerase < 'mttmp'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testtrsm3
NB.
NB. Description:
NB.   Test:
NB.   - xTRTRS (math/lapack2 addon)
NB.   - trsmxxxx (math/mt addon)
NB.   by triangular matrix and multiple RHS
NB.
NB. Syntax:
NB.   testtrsm3 (A ; X)
NB. where
NB.   A - n×n-matrix, triangular
NB.   X - n×3-matrix, exact solutions
NB.
NB. Formula:
NB. - ferr for lnn case:
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
NB. - berr for lnn case:
NB.   foreach i-th computed solution X from nrhs solutions do
NB.     if 0=n or 0=nrhs then
NB.       berr[i] := 0
NB.     elseif 0=||op(A)||_1 or 0=||X|| then
NB.       berr[i] := 1 / FP_EPS
NB.     else
NB.       berr[i] := ((||B - op(A) * X|| / ||op(A)||_1) / ||X||) / FP_EPS
NB.     endif
NB.   endfor
NB.   berr := max(berr[i])
NB.   ||vector|| := norm1t(vector)
NB.
NB. Notes:
NB. - models LAPACK's xTRT02 and xGET04

testtrsm3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/trtrs'

  'A Xv'=. y
  Xh=. |: Xv

  rcondL=.  trlcon1  L=.  trlpick  A
  rcondU=.  trucon1  U=.  trupick  A
  rcondL1=. trl1con1 L1=. (1 ; '') setdiag L
  rcondU1=. tru1con1 U1=. (1 ; '') setdiag U

  'normAxlnn normAxlxn'=. (norm1 , normi) L
  'normAxlnu normAxlxu'=. (norm1 , normi) L1
  'normAxunn normAxuxn'=. (norm1 , normi) U
  'normAxunu normAxuxu'=. (norm1 , normi) U1

  Bllcn=. (ct L ) mp Xv
  Bllcu=. (ct L1) mp Xv
  Bllnn=.     L   mp Xv
  Bllnu=.     L1  mp Xv
  Blltn=. (|: L ) mp Xv
  Blltu=. (|: L1) mp Xv

  Blucn=. (ct U ) mp Xv
  Blucu=. (ct U1) mp Xv
  Blunn=.     U   mp Xv
  Blunu=.     U1  mp Xv
  Blutn=. (|: U ) mp Xv
  Blutu=. (|: U1) mp Xv

  Brlcn=. Xh mp ct L
  Brlcu=. Xh mp ct L1
  Brlnn=. Xh mp    L
  Brlnu=. Xh mp    L1
  Brltn=. Xh mp |: L
  Brltu=. Xh mp |: L1

  Brucn=. Xh mp ct U
  Brucu=. Xh mp ct U1
  Brunn=. Xh mp    U
  Brunu=. Xh mp    U1
  Brutn=. Xh mp |: U
  Brutu=. Xh mp |: U1

  NB. ferr=. (A ; B ; X ; rcondA ; normA) (normitx aferr) Xapprox
  aferr=. 1 : '%&FP_EPS^:(1 > FP_EPS&*)@max_mt_@((3 {:: [) >/@:*@}.`(*`%/ ,: %@FP_EPS)}@, ] (- ,:&u ]) 2 {:: [)`(%@FP_EPS)@.(0 = 3 {:: [)`0:@.(0 e. $@])'
  vferrv_mttmp_=. normitc_mt_ aferr  NB. for Xv
  vferrh_mttmp_=. normitr_mt_ aferr  NB. for Xh

  NB. berr=. (A ; B ; X ; rcondA ; normA) (calcB cberr norm1tx) Xapprox
  cberr=. 2 : 'max_mt_@((FP_EPS , 4 {:: [) %~/@,@,.`(%@FP_EPS)@.(0 ([ = {) ])"1 ] ,.&v (u 0&{::)~ - 1 {:: [)`(%@FP_EPS)@.(0 = 4 {:: [)`0:@.(0 e. $@])'
  NB. vberr for Xv at right side and a triangular matrix embedded in A
  vberrlcn_mttmp_=.  (mp_mt_~ ct_mt_@trlpick_mt_ ) cberr norm1tc_mt_
  vberrlcu_mttmp_=.  (mp_mt_~ ct_mt_@trl1pick_mt_) cberr norm1tc_mt_
  vberrlnn_mttmp_=.  (mp_mt_~        trlpick_mt_ ) cberr norm1tc_mt_
  vberrlnu_mttmp_=.  (mp_mt_~        trl1pick_mt_) cberr norm1tc_mt_
  vberrltn_mttmp_=.  (mp_mt_~ |:    @trlpick_mt_ ) cberr norm1tc_mt_
  vberrltu_mttmp_=.  (mp_mt_~ |:    @trl1pick_mt_) cberr norm1tc_mt_
  vberrucn_mttmp_=.  (mp_mt_~ ct_mt_@trupick_mt_ ) cberr norm1tc_mt_
  vberrucu_mttmp_=.  (mp_mt_~ ct_mt_@tru1pick_mt_) cberr norm1tc_mt_
  vberrunn_mttmp_=.  (mp_mt_~        trupick_mt_ ) cberr norm1tc_mt_
  vberrunu_mttmp_=.  (mp_mt_~        tru1pick_mt_) cberr norm1tc_mt_
  vberrutn_mttmp_=.  (mp_mt_~ |:    @trupick_mt_ ) cberr norm1tc_mt_
  vberrutu_mttmp_=.  (mp_mt_~ |:    @tru1pick_mt_) cberr norm1tc_mt_
  NB. vberr for Xv at right side and a triangular matrix A
  vberrlxcx_mttmp_=. (mp_mt_~ ct_mt_             ) cberr norm1tc_mt_
  vberrlxnx_mttmp_=. (mp_mt_~                    ) cberr norm1tc_mt_
  vberrlxtx_mttmp_=. (mp_mt_~ |:                 ) cberr norm1tc_mt_
  NB. vberr for Xh at left side and a triangular matrix A
  vberrrxcx_mttmp_=. (mp_mt_  ct_mt_             ) cberr norm1tr_mt_
  vberrrxnx_mttmp_=. (mp_mt_                     ) cberr norm1tr_mt_
  vberrrxtx_mttmp_=. (mp_mt_  |:                 ) cberr norm1tr_mt_

  ('''lcn''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrlcn_mttmp_ )) A  ; Bllcn ; Xv ; rcondL  ; normAxlxn
  ('''lcu''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrlcu_mttmp_ )) A  ; Bllcu ; Xv ; rcondL1 ; normAxlxu
  ('''lnn''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrlnn_mttmp_ )) A  ; Bllnn ; Xv ; rcondL  ; normAxlnn
  ('''lnu''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrlnu_mttmp_ )) A  ; Bllnu ; Xv ; rcondL1 ; normAxlnu
  ('''ltn''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrltn_mttmp_ )) A  ; Blltn ; Xv ; rcondL  ; normAxlxn
  ('''ltu''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrltu_mttmp_ )) A  ; Blltu ; Xv ; rcondL1 ; normAxlxu

  ('''ucn''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrucn_mttmp_ )) A  ; Blucn ; Xv ; rcondL  ; normAxuxn
  ('''ucu''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrucu_mttmp_ )) A  ; Blucu ; Xv ; rcondL1 ; normAxuxu
  ('''unn''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrunn_mttmp_ )) A  ; Blunn ; Xv ; rcondL  ; normAxunn
  ('''unu''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrunu_mttmp_ )) A  ; Blunu ; Xv ; rcondL1 ; normAxunu
  ('''utn''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrutn_mttmp_ )) A  ; Blutn ; Xv ; rcondL  ; normAxuxn
  ('''utu''&dtrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrutu_mttmp_ )) A  ; Blutu ; Xv ; rcondL1 ; normAxuxu

  ('''lcn''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrlcn_mttmp_ )) A  ; Bllcn ; Xv ; rcondL  ; normAxlxn
  ('''lcu''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrlcu_mttmp_ )) A  ; Bllcu ; Xv ; rcondL1 ; normAxlxu
  ('''lnn''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrlnn_mttmp_ )) A  ; Bllnn ; Xv ; rcondL  ; normAxlnn
  ('''lnu''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrlnu_mttmp_ )) A  ; Bllnu ; Xv ; rcondL1 ; normAxlnu
  ('''ltn''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrltn_mttmp_ )) A  ; Blltn ; Xv ; rcondL  ; normAxlxn
  ('''ltu''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrltu_mttmp_ )) A  ; Blltu ; Xv ; rcondL1 ; normAxlxu

  ('''ucn''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrucn_mttmp_ )) A  ; Blucn ; Xv ; rcondL  ; normAxuxn
  ('''ucu''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrucu_mttmp_ )) A  ; Blucu ; Xv ; rcondL1 ; normAxuxu
  ('''unn''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrunn_mttmp_ )) A  ; Blunn ; Xv ; rcondL  ; normAxunn
  ('''unu''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrunu_mttmp_ )) A  ; Blunu ; Xv ; rcondL1 ; normAxunu
  ('''utn''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrutn_mttmp_ )) A  ; Blutn ; Xv ; rcondL  ; normAxuxn
  ('''utu''&ztrtrs_mttmp_' tmonad (2&{. `        ]`(rcondL "_)`vferrv_mttmp_`vberrutu_mttmp_ )) A  ; Blutu ; Xv ; rcondL1 ; normAxuxu

  ('trsmllcn'              tdyad  (0&{::`(1&{::)`]`(rcondL "_)`vferrv_mttmp_`vberrlxcx_mttmp_)) L  ; Bllcn ; Xv ; rcondL  ; normAxlxn
  ('trsmllcu'              tdyad  (0&{::`(1&{::)`]`(rcondL1"_)`vferrv_mttmp_`vberrlxcx_mttmp_)) L1 ; Bllcu ; Xv ; rcondL1 ; normAxlxu
  ('trsmllnn'              tdyad  (0&{::`(1&{::)`]`(rcondL "_)`vferrv_mttmp_`vberrlxnx_mttmp_)) L  ; Bllnn ; Xv ; rcondL  ; normAxlnn
  ('trsmllnu'              tdyad  (0&{::`(1&{::)`]`(rcondL1"_)`vferrv_mttmp_`vberrlxnx_mttmp_)) L1 ; Bllnu ; Xv ; rcondL1 ; normAxlnu
  ('trsmlltn'              tdyad  (0&{::`(1&{::)`]`(rcondL "_)`vferrv_mttmp_`vberrlxtx_mttmp_)) L  ; Blltn ; Xv ; rcondL  ; normAxlxn
  ('trsmlltu'              tdyad  (0&{::`(1&{::)`]`(rcondL1"_)`vferrv_mttmp_`vberrlxtx_mttmp_)) L1 ; Blltu ; Xv ; rcondL1 ; normAxlxu

  ('trsmlucn'              tdyad  (0&{::`(1&{::)`]`(rcondU "_)`vferrv_mttmp_`vberrlxcx_mttmp_)) U  ; Blucn ; Xv ; rcondU  ; normAxuxn
  ('trsmlucu'              tdyad  (0&{::`(1&{::)`]`(rcondU1"_)`vferrv_mttmp_`vberrlxcx_mttmp_)) U1 ; Blucu ; Xv ; rcondU1 ; normAxuxu
  ('trsmlunn'              tdyad  (0&{::`(1&{::)`]`(rcondU "_)`vferrv_mttmp_`vberrlxnx_mttmp_)) U  ; Blunn ; Xv ; rcondU  ; normAxunn
  ('trsmlunu'              tdyad  (0&{::`(1&{::)`]`(rcondU1"_)`vferrv_mttmp_`vberrlxnx_mttmp_)) U1 ; Blunu ; Xv ; rcondU1 ; normAxunu
  ('trsmlutn'              tdyad  (0&{::`(1&{::)`]`(rcondU "_)`vferrv_mttmp_`vberrlxtx_mttmp_)) U  ; Blutn ; Xv ; rcondU  ; normAxuxn
  ('trsmlutu'              tdyad  (0&{::`(1&{::)`]`(rcondU1"_)`vferrv_mttmp_`vberrlxtx_mttmp_)) U1 ; Blutu ; Xv ; rcondU1 ; normAxuxu

  ('trsmrlcn'              tdyad  (0&{::`(1&{::)`]`(rcondL "_)`vferrh_mttmp_`vberrrxcx_mttmp_)) L  ; Brlcn ; Xh ; rcondL  ; normAxlxn
  ('trsmrlcu'              tdyad  (0&{::`(1&{::)`]`(rcondL1"_)`vferrh_mttmp_`vberrrxcx_mttmp_)) L1 ; Brlcu ; Xh ; rcondL1 ; normAxlxu
  ('trsmrlnn'              tdyad  (0&{::`(1&{::)`]`(rcondL "_)`vferrh_mttmp_`vberrrxnx_mttmp_)) L  ; Brlnn ; Xh ; rcondL  ; normAxlnn
  ('trsmrlnu'              tdyad  (0&{::`(1&{::)`]`(rcondL1"_)`vferrh_mttmp_`vberrrxnx_mttmp_)) L1 ; Brlnu ; Xh ; rcondL1 ; normAxlnu
  ('trsmrltn'              tdyad  (0&{::`(1&{::)`]`(rcondL "_)`vferrh_mttmp_`vberrrxtx_mttmp_)) L  ; Brltn ; Xh ; rcondL  ; normAxlxn
  ('trsmrltu'              tdyad  (0&{::`(1&{::)`]`(rcondL1"_)`vferrh_mttmp_`vberrrxtx_mttmp_)) L1 ; Brltu ; Xh ; rcondL1 ; normAxlxu

  ('trsmrucn'              tdyad  (0&{::`(1&{::)`]`(rcondU "_)`vferrh_mttmp_`vberrrxcx_mttmp_)) U  ; Brucn ; Xh ; rcondU  ; normAxuxn
  ('trsmrucu'              tdyad  (0&{::`(1&{::)`]`(rcondU1"_)`vferrh_mttmp_`vberrrxcx_mttmp_)) U1 ; Brucu ; Xh ; rcondU1 ; normAxuxu
  ('trsmrunn'              tdyad  (0&{::`(1&{::)`]`(rcondU "_)`vferrh_mttmp_`vberrrxnx_mttmp_)) U  ; Brunn ; Xh ; rcondU  ; normAxunn
  ('trsmrunu'              tdyad  (0&{::`(1&{::)`]`(rcondU1"_)`vferrh_mttmp_`vberrrxnx_mttmp_)) U1 ; Brunu ; Xh ; rcondU1 ; normAxunu
  ('trsmrutn'              tdyad  (0&{::`(1&{::)`]`(rcondU "_)`vferrh_mttmp_`vberrrxtx_mttmp_)) U  ; Brutn ; Xh ; rcondU  ; normAxuxn
  ('trsmrutu'              tdyad  (0&{::`(1&{::)`]`(rcondU1"_)`vferrh_mttmp_`vberrrxtx_mttmp_)) U1 ; Brutu ; Xh ; rcondU1 ; normAxuxu

  coerase < 'mttmp'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testsm
NB.
NB. Description:
NB.   Adv. to make verb to test trsmxxxx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testsm
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
NB.     ?@$&0 testsm_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testsm_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testsm_mt_ 150 150

testsm=: 1 : 'EMPTY [ (testtrsm3_mt_@((}."1 ; {."1)~ _3:) [ testtrsm1_mt_@(_3&(}."1) ; {:"1))@u@(+&0 3)^:(=/)'
