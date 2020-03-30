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
NB. testtrsm  Test trsmxxxx by triangular matrix
NB. testsm    Adv. to make verb to test trsmxxxx by matrix of
NB.           generator and shape given
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
NB. Verb:          Solves:           Syntax:
NB. trsmllcn       L^H  * X = B      Xv=. A trsmllcn Bv
NB. trsmllcu       L1^H * X = B      Xv=. A trsmllcu Bv
NB. trsmllnn       L    * X = B      Xv=. A trsmllnn Bv
NB. trsmllnu       L1   * X = B      Xv=. A trsmllnu Bv
NB. trsmlltn       L^T  * X = B      Xv=. A trsmlltn Bv
NB. trsmlltu       L1^T * X = B      Xv=. A trsmlltu Bv
NB. trsmlucn       U^H  * X = B      Xv=. A trsmlucn Bv
NB. trsmlucu       U1^H * X = B      Xv=. A trsmlucu Bv
NB. trsmlunn       U    * X = B      Xv=. A trsmlunn Bv
NB. trsmlunu       U1   * X = B      Xv=. A trsmlunu Bv
NB. trsmlutn       U^T  * X = B      Xv=. A trsmlutn Bv
NB. trsmlutu       U1^T * X = B      Xv=. A trsmlutu Bv
NB. trsmrlcn       X * L^H  = B      Xh=. A trsmrlcn Bh
NB. trsmrlcu       X * L1^H = B      Xh=. A trsmrlcu Bh
NB. trsmrlnn       X * L    = B      Xh=. A trsmrlnn Bh
NB. trsmrlnu       X * L1   = B      Xh=. A trsmrlnu Bh
NB. trsmrltn       X * L^T  = B      Xh=. A trsmrltn Bh
NB. trsmrltu       X * L1^T = B      Xh=. A trsmrltu Bh
NB. trsmrucn       X * U^H  = B      Xh=. A trsmrucn Bh
NB. trsmrucu       X * U1^H = B      Xh=. A trsmrucu Bh
NB. trsmrunn       X * U    = B      Xh=. A trsmrunn Bh
NB. trsmrunu       X * U1   = B      Xh=. A trsmrunu Bh
NB. trsmrutn       X * U^T  = B      Xh=. A trsmrutn Bh
NB. trsmrutu       X * U1^T = B      Xh=. A trsmrutu Bh
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
NB. testtrsm
NB.
NB. Description:
NB.   Test:
NB.   - trtrs (math/lapack addon)
NB.   - trsmxxxx (math/mt addon)
NB.   by triangular matrix
NB.
NB. Syntax:
NB.   testtrsm (A;X)
NB. where
NB.   A - n×n-matrix, triangular
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB.   ferr := max(||X - exactX|| / ||X||)
NB.   berr := max(||B - op(A) * X|| / (FP_EPS * ||op(A)|| * ||X||))

testtrsm=: 3 : 0
  require :: ] '~addons/math/lapack/lapack.ijs'
  need_jlapack_ :: ] 'trtrs'

  'A X'=. y
  conL=.  trlcon1  L=.  trl  A
  conL1=. trl1con1 L1=. trl1 A
  conU=.  trucon1  U=.  tru  A
  conU1=. tru1con1 U1=. tru1 A

  ('trtrs_jlapack_' tmonad (({.,(mp&.>/))`]`(conU "_)`(normi@(((- %&normic [) 1&{::)~))`(normi@(norm1tc@(( mp      &>/)@[ - ( mp~     0&{::)~) % (FP_EPS * (1:`]@.*)@norm1@(0 {:: [)) * norm1tc@])))) U  ; X

  ('trsmllcn' tdyad ((0&{::)`((mp~ ct)~&>/)`]`(conL "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ ct)~&>/)@[ - ((mp~ ct) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tc@])))) L  ; X
  ('trsmllcu' tdyad ((0&{::)`((mp~ ct)~&>/)`]`(conL1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ ct)~&>/)@[ - ((mp~ ct) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tc@])))) L1 ; X
  ('trsmllnn' tdyad ((0&{::)`( mp      &>/)`]`(conL "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(( mp      &>/)@[ - ( mp~     0&{::)~) % (FP_EPS * (1:`]@.*)@norm1@(0 {:: [)) * norm1tc@])))) L  ; X
  ('trsmllnu' tdyad ((0&{::)`( mp      &>/)`]`(conL1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(( mp      &>/)@[ - ( mp~     0&{::)~) % (FP_EPS * (1:`]@.*)@norm1@(0 {:: [)) * norm1tc@])))) L1 ; X
  ('trsmlltn' tdyad ((0&{::)`((mp~ |:)~&>/)`]`(conL "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ |:)~&>/)@[ - ((mp~ |:) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tc@])))) L  ; X
  ('trsmlltu' tdyad ((0&{::)`((mp~ |:)~&>/)`]`(conL1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ |:)~&>/)@[ - ((mp~ |:) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tc@])))) L1 ; X

  ('trsmlucn' tdyad ((0&{::)`((mp~ ct)~&>/)`]`(conU "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ ct)~&>/)@[ - ((mp~ ct) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tc@])))) U  ; X
  ('trsmlucu' tdyad ((0&{::)`((mp~ ct)~&>/)`]`(conU1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ ct)~&>/)@[ - ((mp~ ct) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tc@])))) U1 ; X
  ('trsmlunn' tdyad ((0&{::)`( mp      &>/)`]`(conU "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(( mp      &>/)@[ - ( mp~     0&{::)~) % (FP_EPS * (1:`]@.*)@norm1@(0 {:: [)) * norm1tc@])))) U  ; X
  ('trsmlunu' tdyad ((0&{::)`( mp      &>/)`]`(conU1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(( mp      &>/)@[ - ( mp~     0&{::)~) % (FP_EPS * (1:`]@.*)@norm1@(0 {:: [)) * norm1tc@])))) U1 ; X
  ('trsmlutn' tdyad ((0&{::)`((mp~ |:)~&>/)`]`(conU "_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ |:)~&>/)@[ - ((mp~ |:) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tc@])))) U  ; X
  ('trsmlutu' tdyad ((0&{::)`((mp~ |:)~&>/)`]`(conU1"_)`(normi@((- %&normic [) 1&{::)~)`(normi@(norm1tc@(((mp~ |:)~&>/)@[ - ((mp~ |:) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tc@])))) U1 ; X

  ('trsmrlcn' tdyad ((0&{::)`((mp  ct)~&>/)`]`(conL "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  ct)~&>/)@[ - ((mp  ct) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tr@])))) L  ; X
  ('trsmrlcu' tdyad ((0&{::)`((mp  ct)~&>/)`]`(conL1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  ct)~&>/)@[ - ((mp  ct) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tr@])))) L1 ; X
  ('trsmrlnn' tdyad ((0&{::)`( mp~     &>/)`]`(conL "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(( mp     ~&>/)@[ - ( mp      0&{::)~) % (FP_EPS * (1:`]@.*)@norm1@(0 {:: [)) * norm1tr@])))) L  ; X
  ('trsmrlnu' tdyad ((0&{::)`( mp~     &>/)`]`(conL1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(( mp     ~&>/)@[ - ( mp      0&{::)~) % (FP_EPS * (1:`]@.*)@norm1@(0 {:: [)) * norm1tr@])))) L1 ; X
  ('trsmrltn' tdyad ((0&{::)`((mp  |:)~&>/)`]`(conL "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  |:)~&>/)@[ - ((mp  |:) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tr@])))) L  ; X
  ('trsmrltu' tdyad ((0&{::)`((mp  |:)~&>/)`]`(conL1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  |:)~&>/)@[ - ((mp  |:) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tr@])))) L1 ; X

  ('trsmrucn' tdyad ((0&{::)`((mp  ct)~&>/)`]`(conU "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  ct)~&>/)@[ - ((mp  ct) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tr@])))) U  ; X
  ('trsmrucu' tdyad ((0&{::)`((mp  ct)~&>/)`]`(conU1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  ct)~&>/)@[ - ((mp  ct) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tr@])))) U1 ; X
  ('trsmrunn' tdyad ((0&{::)`( mp~     &>/)`]`(conU "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(( mp     ~&>/)@[ - ( mp      0&{::)~) % (FP_EPS * (1:`]@.*)@norm1@(0 {:: [)) * norm1tr@])))) U  ; X
  ('trsmrunu' tdyad ((0&{::)`( mp~     &>/)`]`(conU1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(( mp     ~&>/)@[ - ( mp      0&{::)~) % (FP_EPS * (1:`]@.*)@norm1@(0 {:: [)) * norm1tr@])))) U1 ; X
  ('trsmrutn' tdyad ((0&{::)`((mp  |:)~&>/)`]`(conU "_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  |:)~&>/)@[ - ((mp  |:) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tr@])))) U  ; X
  ('trsmrutu' tdyad ((0&{::)`((mp  |:)~&>/)`]`(conU1"_)`(normi@((- %&normir [) 1&{::)~)`(normi@(norm1tr@(((mp  |:)~&>/)@[ - ((mp  |:) 0&{::)~) % (FP_EPS * (1:`]@.*)@normi@(0 {:: [)) * norm1tr@])))) U1 ; X

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
NB.
NB. Notes:
NB. - trsmxxxx are impractical for large matrices

testsm=: 1 : 'EMPTY [ testtrsm_mt_@(u ; u)^:(=/ *. 500 >: {.)'
