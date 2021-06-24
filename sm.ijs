NB. Solve linear monomial equation with triangular matrix
NB.
NB. trsmxxxx   Solve equation (op(A) * X = B) or
NB.            (X * op(A) = B), where A is either unit or
NB.            non-unit, either lower or upper, triangular
NB.            matrix; op(A) is either A itself, or A^T, the
NB.            transposition of A, or A^H, the conjugate
NB.            transposition of A; B is known right-hand
NB.            sides (RHS), X is unknown solutions
NB.
NB. testtrsm1  Test trsmxxxx by triangular matrix and single
NB.            RHS
NB. testtrsm3  Test trsmxxxx by triangular matrix and
NB.            multiple RHS
NB. testsm     Adv. to make verb to test trsmxxxx by matrix of
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
NB. where
NB.   A    - n×n-matrix, containing either U, U1, L or L1
NB.   U    - n×n-matrix, the upper triangular
NB.   U1   - n×n-matrix, the unit upper triangular (diagonal
NB.          is not stored)
NB.   L    - n×n-matrix, the lower triangular
NB.   L1   - n×n-matrix, the unit lower triangular (diagonal
NB.          is not stored)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, solutions
NB.   Xh   - same shape as Bh, solutions
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
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
NB. - models LAPACK's xTRTRS for matrix RHS with the
NB.   following difference: no check for A singularity
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
NB.   Test trsmxxxx by triangular matrix and single RHS
NB.
NB. Syntax:
NB.   testtrsm1 (A ; x)
NB. where
NB.   A - n×n-matrix, triangular
NB.   x - n-vector, the exact solution

testtrsm1=: 3 : 0
  'A x'=. y

  rcondL=.  trlcon1  L=.  trlpick A
  rcondU=.  trucon1  U=.  trupick A
  rcondL1=. trl1con1 L1=. (1 ; '') setdiag L
  rcondU1=. tru1con1 U1=. (1 ; '') setdiag U

  'norm1L  normiL '=. (norm1 , normi) L
  'norm1L1 normiL1'=. (norm1 , normi) L1
  'norm1U  normiU '=. (norm1 , normi) U
  'norm1U1 normiU1'=. (norm1 , normi) U1

  NB. vberrX for x at right side
  vberrlxnx=:  mp~     t02v
  vberrlxcx=: (mp~ ct) t02v
  vberrlxtx=: (mp~ |:) t02v
  NB. vberrX for x at left side
  vberrrxnx=:  mp      t02v
  vberrrxcx=: (mp  ct) t02v
  vberrrxtx=: (mp  |:) t02v

  ('trsmllnn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxnx)) L  ; (L   mp       x ) ; x  ; rcondL  ; norm1L
  ('trsmllnu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxnx)) L1 ; (L1  mp       x ) ; x  ; rcondL1 ; norm1L1
  ('trsmllcn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxcx)) L  ; (L  (mp~ ct)~ x ) ; x  ; rcondL  ; normiL
  ('trsmllcu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxcx)) L1 ; (L1 (mp~ ct)~ x ) ; x  ; rcondL1 ; normiL1
  ('trsmlltn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxtx)) L  ; (L  (mp~ |:)~ x ) ; x  ; rcondL  ; normiL
  ('trsmlltu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxtx)) L1 ; (L1 (mp~ |:)~ x ) ; x  ; rcondL1 ; normiL1

  ('trsmlunn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxnx)) U  ; (U   mp       x ) ; x  ; rcondU  ; norm1U
  ('trsmlunu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxnx)) U1 ; (U1  mp       x ) ; x  ; rcondU1 ; norm1U1
  ('trsmlucn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxcx)) U  ; (U  (mp~ ct)~ x ) ; x  ; rcondU  ; normiU
  ('trsmlucu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxcx)) U1 ; (U1 (mp~ ct)~ x ) ; x  ; rcondU1 ; normiU1
  ('trsmlutn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxtx)) U  ; (U  (mp~ |:)~ x ) ; x  ; rcondU  ; normiU
  ('trsmlutu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrlxtx)) U1 ; (U1 (mp~ |:)~ x ) ; x  ; rcondU1 ; normiU1

  ('trsmrlnn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxnx)) L  ; (x   mp       L ) ; x  ; rcondL  ; normiL
  ('trsmrlnu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxnx)) L1 ; (x   mp       L1) ; x  ; rcondL1 ; normiL1
  ('trsmrlcn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxcx)) L  ; (x  (mp  ct)  L ) ; x  ; rcondL  ; norm1L
  ('trsmrlcu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxcx)) L1 ; (x  (mp  ct)  L1) ; x  ; rcondL1 ; norm1L1
  ('trsmrltn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxtx)) L  ; (x  (mp  |:)  L ) ; x  ; rcondL  ; norm1L
  ('trsmrltu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxtx)) L1 ; (x  (mp  |:)  L1) ; x  ; rcondL1 ; norm1L1

  ('trsmrunn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxnx)) U  ; (x   mp       U ) ; x  ; rcondU  ; normiU
  ('trsmrunu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxnx)) U1 ; (x   mp       U1) ; x  ; rcondU1 ; normiU1
  ('trsmrucn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxcx)) U  ; (x  (mp  ct)  U ) ; x  ; rcondU  ; norm1U
  ('trsmrucu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxcx)) U1 ; (x  (mp  ct)  U1) ; x  ; rcondU1 ; norm1U1
  ('trsmrutn' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxtx)) U  ; (x  (mp  |:)  U ) ; x  ; rcondU  ; norm1U
  ('trsmrutu' tdyad ((0&{::)`(1&{::)`]`(3&{::)`t04v`vberrrxtx)) U1 ; (x  (mp  |:)  U1) ; x  ; rcondU1 ; norm1U1

  coerase < 'mttmp'
  erase 'vberrlxnx vberrlxcx vberrlxtx vberrrxnx vberrrxcx vberrrxtx'

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

testtrsm3=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/trtrs'

  'A Xv'=. y
  Xh=. |: Xv

  rcondL=.  trlcon1  L=.  trlpick  A
  rcondU=.  trucon1  U=.  trupick  A
  rcondL1=. trl1con1 L1=. (1 ; '') setdiag L
  rcondU1=. tru1con1 U1=. (1 ; '') setdiag U

  'norm1L  normiL '=. (norm1 , normi) L
  'norm1L1 normiL1'=. (norm1 , normi) L1
  'norm1U  normiU '=. (norm1 , normi) U
  'norm1U1 normiU1'=. (norm1 , normi) U1

  Bllnn=.     L   mp Xv
  Bllnu=.     L1  mp Xv
  Bllcn=. (ct L ) mp Xv
  Bllcu=. (ct L1) mp Xv
  Blltn=. (|: L ) mp Xv
  Blltu=. (|: L1) mp Xv

  Blunn=.     U   mp Xv
  Blunu=.     U1  mp Xv
  Blucn=. (ct U ) mp Xv
  Blucu=. (ct U1) mp Xv
  Blutn=. (|: U ) mp Xv
  Blutu=. (|: U1) mp Xv

  Brlnn=. Xh mp    L
  Brlnu=. Xh mp    L1
  Brlcn=. Xh mp ct L
  Brlcu=. Xh mp ct L1
  Brltn=. Xh mp |: L
  Brltu=. Xh mp |: L1

  Brunn=. Xh mp    U
  Brunu=. Xh mp    U1
  Brucn=. Xh mp ct U
  Brucu=. Xh mp ct U1
  Brutn=. Xh mp |: U
  Brutu=. Xh mp |: U1

  vferrv=: normitc t04m  NB. for Xv
  vferrh=: normitr t04m  NB. for Xh

  NB. vberrX for Xv at right side and a triangular matrix embedded in A
  vberrlnn=:  (mp~    trlpick ) t02m norm1tc
  vberrlnu=:  (mp~    trl1pick) t02m norm1tc
  vberrlcn=:  (mp~ ct@trlpick ) t02m norm1tc
  vberrlcu=:  (mp~ ct@trl1pick) t02m norm1tc
  vberrltn=:  (mp~ |:@trlpick ) t02m norm1tc
  vberrltu=:  (mp~ |:@trl1pick) t02m norm1tc
  vberrunn=:  (mp~    trupick ) t02m norm1tc
  vberrunu=:  (mp~    tru1pick) t02m norm1tc
  vberrucn=:  (mp~ ct@trupick ) t02m norm1tc
  vberrucu=:  (mp~ ct@tru1pick) t02m norm1tc
  vberrutn=:  (mp~ |:@trupick ) t02m norm1tc
  vberrutu=:  (mp~ |:@tru1pick) t02m norm1tc
  NB. vberrX for Xv at right side and a triangular matrix A
  vberrlxnx=:  mp~              t02m norm1tc
  vberrlxcx=: (mp~ ct         ) t02m norm1tc
  vberrlxtx=: (mp~ |:         ) t02m norm1tc
  NB. vberrX for Xh at left side and a triangular matrix A
  vberrrxnx=:  mp               t02m norm1tr
  vberrrxcx=: (mp  ct         ) t02m norm1tr
  vberrrxtx=: (mp  |:         ) t02m norm1tr

  ('''lnn''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrlnn )) A  ; Bllnn ; Xv ; rcondL  ; norm1L
  ('''lnu''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrlnu )) A  ; Bllnu ; Xv ; rcondL1 ; norm1L1
  ('''lcn''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrlcn )) A  ; Bllcn ; Xv ; rcondL  ; normiL
  ('''lcu''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrlcu )) A  ; Bllcu ; Xv ; rcondL1 ; normiL1
  ('''ltn''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrltn )) A  ; Blltn ; Xv ; rcondL  ; normiL
  ('''ltu''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrltu )) A  ; Blltu ; Xv ; rcondL1 ; normiL1

  ('''unn''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrunn )) A  ; Blunn ; Xv ; rcondL  ; norm1U
  ('''unu''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrunu )) A  ; Blunu ; Xv ; rcondL1 ; norm1U1
  ('''ucn''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrucn )) A  ; Blucn ; Xv ; rcondL  ; normiU
  ('''ucu''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrucu )) A  ; Blucu ; Xv ; rcondL1 ; normiU1
  ('''utn''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrutn )) A  ; Blutn ; Xv ; rcondL  ; normiU
  ('''utu''&dtrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrutu )) A  ; Blutu ; Xv ; rcondL1 ; normiU1

  ('''lnn''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrlnn )) A  ; Bllnn ; Xv ; rcondL  ; norm1L
  ('''lnu''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrlnu )) A  ; Bllnu ; Xv ; rcondL1 ; norm1L1
  ('''lcn''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrlcn )) A  ; Bllcn ; Xv ; rcondL  ; normiL
  ('''lcu''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrlcu )) A  ; Bllcu ; Xv ; rcondL1 ; normiL1
  ('''ltn''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrltn )) A  ; Blltn ; Xv ; rcondL  ; normiL
  ('''ltu''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrltu )) A  ; Blltu ; Xv ; rcondL1 ; normiL1

  ('''unn''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrunn )) A  ; Blunn ; Xv ; rcondL  ; norm1U
  ('''unu''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrunu )) A  ; Blunu ; Xv ; rcondL1 ; norm1U1
  ('''ucn''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrucn )) A  ; Blucn ; Xv ; rcondL  ; normiU
  ('''ucu''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrucu )) A  ; Blucu ; Xv ; rcondL1 ; normiU1
  ('''utn''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrutn )) A  ; Blutn ; Xv ; rcondL  ; normiU
  ('''utu''&ztrtrs_mttmp_' tmonad ((2&{. )`        ]`(3&{::)`vferrv`vberrutu )) A  ; Blutu ; Xv ; rcondL1 ; normiU1

  ('trsmllnn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxnx)) L  ; Bllnn ; Xv ; rcondL  ; norm1L
  ('trsmllnu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxnx)) L1 ; Bllnu ; Xv ; rcondL1 ; norm1L1
  ('trsmllcn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxcx)) L  ; Bllcn ; Xv ; rcondL  ; normiL
  ('trsmllcu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxcx)) L1 ; Bllcu ; Xv ; rcondL1 ; normiL1
  ('trsmlltn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxtx)) L  ; Blltn ; Xv ; rcondL  ; normiL
  ('trsmlltu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxtx)) L1 ; Blltu ; Xv ; rcondL1 ; normiL1

  ('trsmlunn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxnx)) U  ; Blunn ; Xv ; rcondU  ; norm1U
  ('trsmlunu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxnx)) U1 ; Blunu ; Xv ; rcondU1 ; norm1U1
  ('trsmlucn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxcx)) U  ; Blucn ; Xv ; rcondU  ; normiU
  ('trsmlucu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxcx)) U1 ; Blucu ; Xv ; rcondU1 ; normiU1
  ('trsmlutn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxtx)) U  ; Blutn ; Xv ; rcondU  ; normiU
  ('trsmlutu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrv`vberrlxtx)) U1 ; Blutu ; Xv ; rcondU1 ; normiU1

  ('trsmrlnn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxnx)) L  ; Brlnn ; Xh ; rcondL  ; normiL
  ('trsmrlnu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxnx)) L1 ; Brlnu ; Xh ; rcondL1 ; normiL1
  ('trsmrlcn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxcx)) L  ; Brlcn ; Xh ; rcondL  ; norm1L
  ('trsmrlcu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxcx)) L1 ; Brlcu ; Xh ; rcondL1 ; norm1L1
  ('trsmrltn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxtx)) L  ; Brltn ; Xh ; rcondL  ; norm1L
  ('trsmrltu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxtx)) L1 ; Brltu ; Xh ; rcondL1 ; norm1L1

  ('trsmrunn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxnx)) U  ; Brunn ; Xh ; rcondU  ; normiU
  ('trsmrunu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxnx)) U1 ; Brunu ; Xh ; rcondU1 ; normiU1
  ('trsmrucn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxcx)) U  ; Brucn ; Xh ; rcondU  ; norm1U
  ('trsmrucu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxcx)) U1 ; Brucu ; Xh ; rcondU1 ; norm1U1
  ('trsmrutn'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxtx)) U  ; Brutn ; Xh ; rcondU  ; norm1U
  ('trsmrutu'              tdyad  ((0&{::)`(1&{::)`]`(3&{::)`vferrh`vberrrxtx)) U1 ; Brutu ; Xh ; rcondU1 ; norm1U1

  coerase < 'mttmp'
  erase 'vferrv vferrh vberrlnn vberrlnu vberrlcn vberrlcu vberrltn vberrltu vberrunn vberrunu vberrucn vberrucu vberrutn vberrutu vberrlxnx vberrlxcx vberrlxtx vberrrxnx vberrrxcx vberrrxtx'

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

testsm=: 1 : 'EMPTY [ ((}."1 (testtrsm3_mt_@; [ testtrsm1_mt_@(; {:"1)) {."1)~ _3:)@u@(+&0 3)^:(=/)'
