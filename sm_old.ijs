NB. Solve linear monomial equation with triangular matrix
NB.
NB. trsmxxxx  Solve equation (op(A) * X = B) or
NB.           (X * op(A) = B), where A is either unit or
NB.           non-unit, either lower or upper, triangular
NB.           matrix; op(A) is either A itself, or A^T, the
NB.           transposition of A, or A^H, the conjugate
NB.           transposition of A; B is known right-hand side
NB.           (RHS), X is unknown solution
NB.
NB. testtrsm_old  Test trsmxxxx by triangular matrix given
NB. testsm_old    Adv. to make verb to test trsmxxxx by matrix of
NB.           generator and shape given
NB.
NB. Version: 0.7.0 2011-08-06
NB.
NB. Copyright 2010-2011 Igor Zhuravlov
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
NB. Verb:          Solves:           Syntax:
NB. trsmlx_old         L    * X = B      Xv=. A trsmlx_old   Bv
NB. trsml1x_old        L1   * X = B      Xv=. A trsml1x_old  Bv
NB. trsml1hx_old       L1^H * X = B      Xv=. A trsml1hx_old Bv
NB. trsml1tx_old       L1^T * X = B      Xv=. A trsml1tx_old Bv
NB. trsmlhx_old        L^H  * X = B      Xv=. A trsmlhx_old  Bv
NB. trsmltx_old        L^T  * X = B      Xv=. A trsmltx_old  Bv
NB. trsmux_old         U    * X = B      Xv=. A trsmux_old   Bv
NB. trsmu1x_old        U1   * X = B      Xv=. A trsmu1x_old  Bv
NB. trsmu1hx_old       U1^H * X = B      Xv=. A trsmu1hx_old Bv
NB. trsmu1tx_old       U1^T * X = B      Xv=. A trsmu1tx_old Bv
NB. trsmuhx_old        U^H  * X = B      Xv=. A trsmuhx_old  Bv
NB. trsmutx_old        U^T  * X = B      Xv=. A trsmutx_old  Bv
NB. trsmxl_old         X * L    = B      Xh=. A trsmxl_old   Bh
NB. trsmxl1_old        X * L1   = B      Xh=. A trsmxl1_old  Bh
NB. trsmxl1h_old       X * L1^H = B      Xh=. A trsmxl1h_old Bh
NB. trsmxl1t_old       X * L1^T = B      Xh=. A trsmxl1t_old Bh
NB. trsmxlh_old        X * L^H  = B      Xh=. A trsmxlh_old  Bh
NB. trsmxlt_old        X * L^T  = B      Xh=. A trsmxlt_old  Bh
NB. trsmxu_old         X * U    = B      Xh=. A trsmxu_old   Bh
NB. trsmxu1_old        X * U1   = B      Xh=. A trsmxu1_old  Bh
NB. trsmxu1h_old       X * U1^H = B      Xh=. A trsmxu1h_old Bh
NB. trsmxu1t_old       X * U1^T = B      Xh=. A trsmxu1t_old Bh
NB. trsmxuh_old        X * U^H  = B      Xh=. A trsmxuh_old  Bh
NB. trsmxut_old        X * U^T  = B      Xh=. A trsmxut_old  Bh
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
NB. - model BLAS's xTRSV for vector RHS
NB. - model BLAS's xTRSM for matrix RHS with following
NB.   difference: alpha parameter is assumed to be always
NB.   equal to 1; to workaround this limitation use the
NB.   following pattern:
NB.     X=. A trsmxxxx alpha*B
NB.
NB. TODO:
NB. - replace column-wise algos by row-wise

trsmlx_old=:   ((((    #@]  {   (1 {:: [)) - ] mp~ (( 1 liosW)&# 0&{::)~   ({,) 0 {:: [) % (    (*>:) &# 0&{::)~   ({,) 0 {:: [) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. lios(li)=(1 dhs2lios (   i*n,i)), lio(lii)=   i*(n+1)
trsml1x_old=:  ( ((    #@]  {   (1 {:: [)) - ] mp~ (( 1 liosW)&# 0&{::)~   ({,) 0 {:: [                                        ) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. lios(li)=(1 dhs2lios (   i*n,i))
trsml1hx_old=: ( (((_1-#@]) {   (1 {:: [)) - ] mp~ ((_1 liosS)&# 0&{::)~ +@({,) 0 {:: [                                        ) ,      ])^:(;`(#@])`(0 {.   ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i))
trsml1tx_old=: ( (((_1-#@]) {   (1 {:: [)) - ] mp~ ((_1 liosS)&# 0&{::)~   ({,) 0 {:: [                                        ) ,      ])^:(;`(#@])`(0 {.   ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i))
trsmlhx_old=:  (((((_1-#@]) {   (1 {:: [)) - ] mp~ ((_1 liosS)&# 0&{::)~ +@({,) 0 {:: [) % ((_1-(*>:))&# 0&{::)~ +@({,) 0 {:: [) ,      ])^:(;`(#@])`(0 {.   ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i)), lio(lii)=-1-i*(n+1)
trsmltx_old=:  (((((_1-#@]) {   (1 {:: [)) - ] mp~ ((_1 liosS)&# 0&{::)~   ({,) 0 {:: [) % ((_1-(*>:))&# 0&{::)~   ({,) 0 {:: [) ,      ])^:(;`(#@])`(0 {.   ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i)), lio(lii)=-1-i*(n+1)
trsmux_old=:   (((((_1-#@]) {   (1 {:: [)) - ] mp~ ((_1 liosE)&# 0&{::)~   ({,) 0 {:: [) % ((_1-(*>:))&# 0&{::)~   ({,) 0 {:: [) ,      ])^:(;`(#@])`(0 {.   ]))  NB. lios(ui)=(1 dhs2lios (-1-i*n,i)), lio(lii)=-1-i*(n+1)
trsmu1x_old=:  ( (((_1-#@]) {   (1 {:: [)) - ] mp~ ((_1 liosE)&# 0&{::)~   ({,) 0 {:: [                                        ) ,      ])^:(;`(#@])`(0 {.   ]))  NB. lios(ui)=(1 dhs2lios (-1-i*n,i))
trsmu1hx_old=: ( ((    #@]  {   (1 {:: [)) - ] mp~ (( 1 liosN)&# 0&{::)~ +@({,) 0 {:: [                                        ) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. lios(ui)=(n dhs2lios (   i  ,i))
trsmu1tx_old=: ( ((    #@]  {   (1 {:: [)) - ] mp~ (( 1 liosN)&# 0&{::)~   ({,) 0 {:: [                                        ) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. lios(ui)=(n dhs2lios (   i  ,i))
trsmuhx_old=:  ((((    #@]  {   (1 {:: [)) - ] mp~ (( 1 liosN)&# 0&{::)~ +@({,) 0 {:: [) % (    (*>:) &# 0&{::)~ +@({,) 0 {:: [) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. lios(ui)=(n dhs2lios (   i  ,i)), lio(lii)=   i*(n+1)
trsmutx_old=:  ((((    #@]  {   (1 {:: [)) - ] mp~ (( 1 liosN)&# 0&{::)~   ({,) 0 {:: [) % (    (*>:) &# 0&{::)~   ({,) 0 {:: [) ,    ~ ])^:(;`(#@])`(0 {.   ]))  NB. lios(ui)=(n dhs2lios (   i  ,i)), lio(lii)=   i*(n+1)
trsmxl_old=:   (((((_1-c@]) {"1 (1 {:: [)) - ] mp  ((_1 liosS)&c 0&{::)~   ({,) 0 {:: [) % ((_1-(*>:))&c 0&{::)~   ({,) 0 {:: [) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i)), lio(lii)=-1-i*(n+1)
trsmxl1_old=:  ( (((_1-c@]) {"1 (1 {:: [)) - ] mp  ((_1 liosS)&c 0&{::)~   ({,) 0 {:: [                                        ) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i))
trsmxl1h_old=: ( ((    c@]  {"1 (1 {:: [)) - ] mp  (( 1 liosW)&c 0&{::)~ +@({,) 0 {:: [                                        ) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (   i*n,i))
trsmxl1t_old=: ( ((    c@]  {"1 (1 {:: [)) - ] mp  (( 1 liosW)&c 0&{::)~   ({,) 0 {:: [                                        ) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (   i*n,i))
trsmxlh_old=:  ((((    c@]  {"1 (1 {:: [)) - ] mp  (( 1 liosW)&c 0&{::)~ +@({,) 0 {:: [) % (    (*>:) &c 0&{::)~ +@({,) 0 {:: [) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (   i*n,i)), lio(lii)=   i*(n+1)
trsmxlt_old=:  ((((    c@]  {"1 (1 {:: [)) - ] mp  (( 1 liosW)&c 0&{::)~   ({,) 0 {:: [) % (    (*>:) &c 0&{::)~   ({,) 0 {:: [) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (   i*n,i)), lio(lii)=   i*(n+1)
trsmxu_old=:   ((((    c@]  {"1 (1 {:: [)) - ] mp  (( 1 liosN)&c 0&{::)~   ({,) 0 {:: [) % (    (*>:) &c 0&{::)~   ({,) 0 {:: [) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(n dhs2lios (   i  ,i)), lio(lii)=   i*(n+1)
trsmxu1_old=:  ( ((    c@]  {"1 (1 {:: [)) - ] mp  (( 1 liosN)&c 0&{::)~   ({,) 0 {:: [                                        ) ,"1 0~ ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(n dhs2lios (   i  ,i))
trsmxu1h_old=: ( (((_1-c@]) {"1 (1 {:: [)) - ] mp  ((_1 liosE)&c 0&{::)~ +@({,) 0 {:: [                                        ) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (-1-i*n,i))
trsmxu1t_old=: ( (((_1-c@]) {"1 (1 {:: [)) - ] mp  ((_1 liosE)&c 0&{::)~   ({,) 0 {:: [                                        ) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (-1-i*n,i))
trsmxuh_old=:  (((((_1-c@]) {"1 (1 {:: [)) - ] mp  ((_1 liosE)&c 0&{::)~ +@({,) 0 {:: [) % ((_1-(*>:))&c 0&{::)~ +@({,) 0 {:: [) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (-1-i*n,i)), lio(lii)=-1-i*(n+1)
trsmxut_old=:  (((((_1-c@]) {"1 (1 {:: [)) - ] mp  ((_1 liosE)&c 0&{::)~   ({,) 0 {:: [) % ((_1-(*>:))&c 0&{::)~   ({,) 0 {:: [) ,"0 1  ])^:(;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (-1-i*n,i)), lio(lii)=-1-i*(n+1)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testtrsm_old
NB.
NB. Description:
NB.   Test:
NB.   - trtrs (math/lapack addon)
NB.   - trsmxxxx (math/mt addon)
NB.   by triangular matrix given
NB.
NB. Syntax:
NB.   testtrsm_old (A;X)
NB. where
NB.   A - n×n-matrix, triangular
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB.   ferr := max(||X - exactX|| / ||X||)
NB.   berr := max(||B - op(A) * X|| / (FP_EPS * ||op(A)|| * ||X||))

testtrsm_old=: 3 : 0
  require :: ] '~addons/math/lapack/lapack.ijs'
  need_jlapack_ :: ] 'trtrs'

  'A X'=. y
  'L L1 U U1'=. bT=. (trl ; trl1 ; tru ; tru1) A
  'conL conL1 conU conU1'=. ((trlcon1&.>)`(trl1con1&.>)`(trucon1&.>)`(tru1con1&.>)) ag bT

  ('trtrs_jlapack_' tmonad (({.,(mp&.>/))`]`(conU "_)`(normi@(((- %&normic [) 1&{::)~))`(normi@(norm1tc@((mp&>/)@[ - (mp~ 0&{::)~) % (FP_EPS*norm1@(0 {:: [))*norm1tc@])))) (U;X)

  ('trsmlx_old'   tdyad ((0&{::)`( mp      &>/)`]`(conL "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(( mp      &>/)@[ - ( mp~     0&{::)~) % (FP_EPS*norm1@(0 {:: [))*norm1tc@])))) (L ;X)
  ('trsml1x_old'  tdyad ((0&{::)`( mp      &>/)`]`(conL1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(( mp      &>/)@[ - ( mp~     0&{::)~) % (FP_EPS*norm1@(0 {:: [))*norm1tc@])))) (L1;X)
  ('trsml1hx_old' tdyad ((0&{::)`((mp~ ct)~&>/)`]`(conL1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ ct)~&>/)@[ - ((mp~ ct) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tc@])))) (L1;X)
  ('trsml1tx_old' tdyad ((0&{::)`((mp~ |:)~&>/)`]`(conL1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ |:)~&>/)@[ - ((mp~ |:) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tc@])))) (L1;X)
  ('trsmlhx_old'  tdyad ((0&{::)`((mp~ ct)~&>/)`]`(conL "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ ct)~&>/)@[ - ((mp~ ct) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tc@])))) (L ;X)
  ('trsmltx_old'  tdyad ((0&{::)`((mp~ |:)~&>/)`]`(conL "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ |:)~&>/)@[ - ((mp~ |:) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tc@])))) (L ;X)
  ('trsmux_old'   tdyad ((0&{::)`( mp      &>/)`]`(conU "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(( mp      &>/)@[ - ( mp~     0&{::)~) % (FP_EPS*norm1@(0 {:: [))*norm1tc@])))) (U ;X)
  ('trsmu1x_old'  tdyad ((0&{::)`( mp      &>/)`]`(conU1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(( mp      &>/)@[ - ( mp~     0&{::)~) % (FP_EPS*norm1@(0 {:: [))*norm1tc@])))) (U1;X)
  ('trsmu1hx_old' tdyad ((0&{::)`((mp~ ct)~&>/)`]`(conU1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ ct)~&>/)@[ - ((mp~ ct) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tc@])))) (U1;X)
  ('trsmu1tx_old' tdyad ((0&{::)`((mp~ |:)~&>/)`]`(conU1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ |:)~&>/)@[ - ((mp~ |:) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tc@])))) (U1;X)
  ('trsmuhx_old'  tdyad ((0&{::)`((mp~ ct)~&>/)`]`(conU "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ ct)~&>/)@[ - ((mp~ ct) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tc@])))) (U ;X)
  ('trsmutx_old'  tdyad ((0&{::)`((mp~ |:)~&>/)`]`(conU "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ |:)~&>/)@[ - ((mp~ |:) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tc@])))) (U ;X)
  ('trsmxl_old'   tdyad ((0&{::)`( mp~     &>/)`]`(conL "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(( mp     ~&>/)@[ - ( mp      0&{::)~) % (FP_EPS*norm1@(0 {:: [))*norm1tr@])))) (L ;X)
  ('trsmxl1_old'  tdyad ((0&{::)`( mp~     &>/)`]`(conL1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(( mp     ~&>/)@[ - ( mp      0&{::)~) % (FP_EPS*norm1@(0 {:: [))*norm1tr@])))) (L1;X)
  ('trsmxl1h_old' tdyad ((0&{::)`((mp  ct)~&>/)`]`(conL1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  ct)~&>/)@[ - ((mp  ct) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tr@])))) (L1;X)
  ('trsmxl1t_old' tdyad ((0&{::)`((mp  |:)~&>/)`]`(conL1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  |:)~&>/)@[ - ((mp  |:) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tr@])))) (L1;X)
  ('trsmxlh_old'  tdyad ((0&{::)`((mp  ct)~&>/)`]`(conL "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  ct)~&>/)@[ - ((mp  ct) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tr@])))) (L ;X)
  ('trsmxlt_old'  tdyad ((0&{::)`((mp  |:)~&>/)`]`(conL "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  |:)~&>/)@[ - ((mp  |:) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tr@])))) (L ;X)
  ('trsmxu_old'   tdyad ((0&{::)`( mp~     &>/)`]`(conU "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(( mp     ~&>/)@[ - ( mp      0&{::)~) % (FP_EPS*norm1@(0 {:: [))*norm1tr@])))) (U ;X)
  ('trsmxu1_old'  tdyad ((0&{::)`( mp~     &>/)`]`(conU1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(( mp     ~&>/)@[ - ( mp      0&{::)~) % (FP_EPS*norm1@(0 {:: [))*norm1tr@])))) (U1;X)
  ('trsmxu1h_old' tdyad ((0&{::)`((mp  ct)~&>/)`]`(conU1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  ct)~&>/)@[ - ((mp  ct) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tr@])))) (U1;X)
  ('trsmxu1t_old' tdyad ((0&{::)`((mp  |:)~&>/)`]`(conU1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  |:)~&>/)@[ - ((mp  |:) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tr@])))) (U1;X)
  ('trsmxuh_old'  tdyad ((0&{::)`((mp  ct)~&>/)`]`(conU "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  ct)~&>/)@[ - ((mp  ct) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tr@])))) (U ;X)
  ('trsmxut_old'  tdyad ((0&{::)`((mp  |:)~&>/)`]`(conU "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  |:)~&>/)@[ - ((mp  |:) 0&{::)~) % (FP_EPS*normi@(0 {:: [))*norm1tr@])))) (U ;X)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testsm_old
NB.
NB. Description:
NB.   Adv. to make verb to test trsmxxxx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testsm_old
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
NB.     ?@$&0 testsm_old_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testsm_old_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testsm_old_mt_ 150 150
NB.
NB. Notes:
NB. - trsmxxxx are impractical for large matrices

testsm_old=: 1 : 'EMPTY_mt_ [ testtrsm_old_mt_@(u ; u)^:(=/ *. 500 >: {.)'
