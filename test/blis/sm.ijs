NB. Solve matrix equation
NB.
NB. xtrsmxxxx  Solve matrix equation with triangular matrix
NB.
NB. Version: 0.14.0 2023-07-04
NB.
NB. Copyright 2023 Igor Zhuravlov
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

coclass 'mtbli'

NB. =========================================================
NB. Includes

require 'math/mt/test/blis/blis'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Dyad         API       Domain
NB. trsmcore     object    mixed
NB. dtrsmcore    typed     real
NB. ztrsmcore    typed     complex
NB.
NB. Description:
NB.   Solves the matrix equation:
NB.     op(A) * X = alpha * B  (1)
NB.   or
NB.     X * op(A) = alpha * B  (2)
NB.   where A is triangular, and op(A) is either A, A^T,
NB.   conj(A) or A^H
NB.
NB. Syntax:
NB.   X=. (side ; uplo ; trans ; diag) xtrsmcore alpha ; AA ; B
NB. where
NB.   side  - 4-vector of chars, side_t as bitstring,
NB.           defines the side of A:
NB.             LEFT               NB. to solve (1) (op(A) is on the left)
NB.             RIGHT              NB. to solve (2) (op(A) is on the right)
NB.   uplo  - 4-vector of chars, uplo_t as bitstring,
NB.           defines the part of A to be read:
NB.             LOWER              NB. LT
NB.             UPPER              NB. UT
NB.   trans - 4-vector of chars, trans_t as bitstring,
NB.           defines the form of op(A):
NB.             NO_TRANSPOSE       NB. op(A) := A        (no transpose)
NB.             TRANSPOSE          NB. op(A) := A^T      (transpose)
NB.             CONJ_NO_TRANSPOSE  NB. op(A) := conj(A)  (conjugate)
NB.             CONJ_TRANSPOSE     NB. op(A) := A^H      (conjugate and transpose)
NB.   diag  - 4-vector of chars, diag_t as bitstring,
NB.           defines the form of A:
NB.             NONUNIT_DIAG       NB. A is either L or U
NB.             UNIT_DIAG          NB. A is either L1 or U1, diagonal
NB.                                NB.   elements of A are not referenced
NB.   alpha - scalar
NB.   AA    - k×k-matrix, contains either non-zero or both
NB.           part(s) of A
NB.   B     - m×n-matrix, non-scaled RHS
NB.   X     - m×n-matrix, solutions
NB.   A     - k×k-matrix, triangular
NB.   m     ≥ 0, the number of rows in B and X
NB.   n     ≥ 0, the number of columns in B and X
NB.   k     = m for (1) or k = n for (2)

trsmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  y=. memu&.>&.(_1&{) y                   NB. unalias B
  objs=. obja L: 0 y                      NB. allocate BLIS objects bonded to J nouns
  (TRIANGULAR 0} x) obj_set_struc`obj_set_uplo`obj_set_conjtrans`obj_set_diag"1 0 (1;1) {:: objs
    NB. set A structure, select A part, transpose A optionally, set A diag type
    NB. note: changing the object is being the op with side-effect so
    NB.       it is possible to do it in parallel as here
  trsm_cd ({. _2 ic {. x) ; <@{:L:0 objs  NB. call bli_trsm() with all but head params marked as object address
  objf L: 0 objs                          NB. free BLIS objects
  _1 {:: y                                NB. return changed copy of B
)

dtrsmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  'side uplo trans diag'=. _2 ic , x
  'alpha AA B'=. y
  mn=. (RIGHT -: {. x) { 'm n'=. $ B
  11 {:: dtrsm_cd side ; uplo ; trans ; diag ; m ; n ; (, alpha) ; AA ; mn ; 1 ; B ; n ; 1
)

ztrsmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  'side uplo trans diag'=. _2 ic , x
  'alpha AA B'=. y
  mn=. (RIGHT -: {. x) { 'm n'=. $ B
  11 {:: ztrsm_cd side ; uplo ; trans ; diag ; m ; n ; (, alpha) ; AA ; mn ; 1 ; B ; n ; 1
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad        API       Domain     Side    A     Reads in A    op1(A)
NB. trsmllnn     object    mixed      (1)     L      LT                A
NB. trsmllnu     object    mixed      (1)     L1    SLT                A
NB. trsmlltn     object    mixed      (1)     L      LT                A^T
NB. trsmlltu     object    mixed      (1)     L1    SLT                A^T
NB. trsmlljn     object    mixed      (1)     L      LT           conj(A)
NB. trsmllju     object    mixed      (1)     L1    SLT           conj(A)
NB. trsmllcn     object    mixed      (1)     L      LT                A^H
NB. trsmllcu     object    mixed      (1)     L1    SLT                A^H
NB. trsmlunn     object    mixed      (1)     U      UT                A
NB. trsmlunu     object    mixed      (1)     U1    SUT                A
NB. trsmlutn     object    mixed      (1)     U      UT                A^T
NB. trsmlutu     object    mixed      (1)     U1    SUT                A^T
NB. trsmlujn     object    mixed      (1)     U      UT           conj(A)
NB. trsmluju     object    mixed      (1)     U1    SUT           conj(A)
NB. trsmlucn     object    mixed      (1)     U      UT                A^H
NB. trsmlucu     object    mixed      (1)     U1    SUT                A^H
NB. trsmrlnn     object    mixed      (2)     L     LT                 A
NB. trsmrlnu     object    mixed      (2)     L1    SLT                A
NB. trsmrltn     object    mixed      (2)     L      LT                A^T
NB. trsmrltu     object    mixed      (2)     L1    SLT                A^T
NB. trsmrljn     object    mixed      (2)     L      LT           conj(A)
NB. trsmrlju     object    mixed      (2)     L1    SLT           conj(A)
NB. trsmrlcn     object    mixed      (2)     L      LT                A^H
NB. trsmrlcu     object    mixed      (2)     L1    SLT                A^H
NB. trsmrunn     object    mixed      (2)     U     UT                 A
NB. trsmrunu     object    mixed      (2)     U1    SUT                A
NB. trsmrutn     object    mixed      (2)     U      UT                A^T
NB. trsmrutu     object    mixed      (2)     U1    SUT                A^T
NB. trsmrujn     object    mixed      (2)     U      UT           conj(A)
NB. trsmruju     object    mixed      (2)     U1    SUT           conj(A)
NB. trsmrucn     object    mixed      (2)     U      UT                A^H
NB. trsmrucu     object    mixed      (2)     U1    SUT                A^H
NB. dtrsmllnn    typed     real       (1)     L      LT                A
NB. dtrsmllnu    typed     real       (1)     L1    SLT                A
NB. dtrsmlltn    typed     real       (1)     L      LT                A^T
NB. dtrsmlltu    typed     real       (1)     L1    SLT                A^T
NB. dtrsmlunn    typed     real       (1)     U      UT                A
NB. dtrsmlunu    typed     real       (1)     U1    SUT                A
NB. dtrsmlutn    typed     real       (1)     U      UT                A^T
NB. dtrsmlutu    typed     real       (1)     U1    SUT                A^T
NB. dtrsmrlnn    typed     real       (2)     L     LT                 A
NB. dtrsmrlnu    typed     real       (2)     L1    SLT                A
NB. dtrsmrltn    typed     real       (2)     L      LT                A^T
NB. dtrsmrltu    typed     real       (2)     L1    SLT                A^T
NB. dtrsmrunn    typed     real       (2)     U     UT                 A
NB. dtrsmrunu    typed     real       (2)     U1    SUT                A
NB. dtrsmrutn    typed     real       (2)     U      UT                A^T
NB. dtrsmrutu    typed     real       (2)     U1    SUT                A^T
NB. ztrsmllnn    typed     complex    (1)     L      LT                A
NB. ztrsmllnu    typed     complex    (1)     L1    SLT                A
NB. ztrsmlltn    typed     complex    (1)     L      LT                A^T
NB. ztrsmlltu    typed     complex    (1)     L1    SLT                A^T
NB. ztrsmlljn    typed     complex    (1)     L      LT           conj(A)
NB. ztrsmllju    typed     complex    (1)     L1    SLT           conj(A)
NB. ztrsmllcn    typed     complex    (1)     L      LT                A^H
NB. ztrsmllcu    typed     complex    (1)     L1    SLT                A^H
NB. ztrsmlunn    typed     complex    (1)     U      UT                A
NB. ztrsmlunu    typed     complex    (1)     U1    SUT                A
NB. ztrsmlutn    typed     complex    (1)     U      UT                A^T
NB. ztrsmlutu    typed     complex    (1)     U1    SUT                A^T
NB. ztrsmlujn    typed     complex    (1)     U      UT           conj(A)
NB. ztrsmluju    typed     complex    (1)     U1    SUT           conj(A)
NB. ztrsmlucn    typed     complex    (1)     U      UT                A^H
NB. ztrsmlucu    typed     complex    (1)     U1    SUT                A^H
NB. ztrsmrlnn    typed     complex    (2)     L     LT                 A
NB. ztrsmrlnu    typed     complex    (2)     L1    SLT                A
NB. ztrsmrltn    typed     complex    (2)     L      LT                A^T
NB. ztrsmrltu    typed     complex    (2)     L1    SLT                A^T
NB. ztrsmrljn    typed     complex    (2)     L      LT           conj(A)
NB. ztrsmrlju    typed     complex    (2)     L1    SLT           conj(A)
NB. ztrsmrlcn    typed     complex    (2)     L      LT                A^H
NB. ztrsmrlcu    typed     complex    (2)     L1    SLT                A^H
NB. ztrsmrunn    typed     complex    (2)     U     UT                 A
NB. ztrsmrunu    typed     complex    (2)     U1    SUT                A
NB. ztrsmrutn    typed     complex    (2)     U      UT                A^T
NB. ztrsmrutu    typed     complex    (2)     U1    SUT                A^T
NB. ztrsmrujn    typed     complex    (2)     U      UT           conj(A)
NB. ztrsmruju    typed     complex    (2)     U1    SUT           conj(A)
NB. ztrsmrucn    typed     complex    (2)     U      UT                A^H
NB. ztrsmrucu    typed     complex    (2)     U1    SUT                A^H
NB.
NB. Description:
NB.   Solves the matrix equation:
NB.     op(A) * X = alpha * B  (1)
NB.   or
NB.     X * op(A) = alpha * B  (2)
NB.   where A is triangular, and op(A) is either A, A^T,
NB.   conj(A) or A^H
NB.
NB. Syntax:
NB.   X=. xtrsmxxxx alpha ; AA ; B
NB. where
NB.   alpha - scalar
NB.   AA    - k×k-matrix, contains either non-zero or both
NB.           part(s) of A
NB.   B     - m×n-matrix, non-scaled RHS
NB.   X     - m×n-matrix, solutions
NB.   A     - k×k-matrix, triangular
NB.   m     ≥ 0, the number of rows in B and X
NB.   n     ≥ 0, the number of columns in B and X
NB.   k     = m for xtrsmlxxx or k = n for xtrsmrxxx
NB.
NB. Notes:
NB. - monad        provides BLIS'
NB.   trsmxxxx     bli_trsm (...)
NB.   dtrsmxxxx    bli_dtrsm(...)
NB.   ztrsmxxxx    bli_ztrsm(...)

trsmllnn=:  (LEFT  , LOWER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmllnu=:  (LEFT  , LOWER ,      NO_TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmlltn=:  (LEFT  , LOWER ,         TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmlltu=:  (LEFT  , LOWER ,         TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmlljn=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmllju=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmllcn=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmllcu=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmlunn=:  (LEFT  , UPPER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmlunu=:  (LEFT  , UPPER ,      NO_TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmlutn=:  (LEFT  , UPPER ,         TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmlutu=:  (LEFT  , UPPER ,         TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmlujn=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmluju=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmlucn=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmlucu=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmrlnn=:  (RIGHT , LOWER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmrlnu=:  (RIGHT , LOWER ,      NO_TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmrltn=:  (RIGHT , LOWER ,         TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmrltu=:  (RIGHT , LOWER ,         TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmrljn=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmrlju=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmrlcn=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmrlcu=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmrunn=:  (RIGHT , UPPER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmrunu=:  (RIGHT , UPPER ,      NO_TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmrutn=:  (RIGHT , UPPER ,         TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmrutu=:  (RIGHT , UPPER ,         TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmrujn=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmruju=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,:    UNIT_DIAG)& trsmcore
trsmrucn=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE ,: NONUNIT_DIAG)& trsmcore
trsmrucu=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE ,:    UNIT_DIAG)& trsmcore

dtrsmllnn=: (LEFT  , LOWER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)&dtrsmcore
dtrsmllnu=: (LEFT  , LOWER ,      NO_TRANSPOSE ,:    UNIT_DIAG)&dtrsmcore
dtrsmlltn=: (LEFT  , LOWER ,         TRANSPOSE ,: NONUNIT_DIAG)&dtrsmcore
dtrsmlltu=: (LEFT  , LOWER ,         TRANSPOSE ,:    UNIT_DIAG)&dtrsmcore
dtrsmlunn=: (LEFT  , UPPER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)&dtrsmcore
dtrsmlunu=: (LEFT  , UPPER ,      NO_TRANSPOSE ,:    UNIT_DIAG)&dtrsmcore
dtrsmlutn=: (LEFT  , UPPER ,         TRANSPOSE ,: NONUNIT_DIAG)&dtrsmcore
dtrsmlutu=: (LEFT  , UPPER ,         TRANSPOSE ,:    UNIT_DIAG)&dtrsmcore
dtrsmrlnn=: (RIGHT , LOWER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)&dtrsmcore
dtrsmrlnu=: (RIGHT , LOWER ,      NO_TRANSPOSE ,:    UNIT_DIAG)&dtrsmcore
dtrsmrltn=: (RIGHT , LOWER ,         TRANSPOSE ,: NONUNIT_DIAG)&dtrsmcore
dtrsmrltu=: (RIGHT , LOWER ,         TRANSPOSE ,:    UNIT_DIAG)&dtrsmcore
dtrsmrunn=: (RIGHT , UPPER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)&dtrsmcore
dtrsmrunu=: (RIGHT , UPPER ,      NO_TRANSPOSE ,:    UNIT_DIAG)&dtrsmcore
dtrsmrutn=: (RIGHT , UPPER ,         TRANSPOSE ,: NONUNIT_DIAG)&dtrsmcore
dtrsmrutu=: (RIGHT , UPPER ,         TRANSPOSE ,:    UNIT_DIAG)&dtrsmcore

ztrsmllnn=: (LEFT  , LOWER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmllnu=: (LEFT  , LOWER ,      NO_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmlltn=: (LEFT  , LOWER ,         TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmlltu=: (LEFT  , LOWER ,         TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmlljn=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmllju=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmllcn=: (LEFT  , LOWER ,    CONJ_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmllcu=: (LEFT  , LOWER ,    CONJ_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmlunn=: (LEFT  , UPPER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmlunu=: (LEFT  , UPPER ,      NO_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmlutn=: (LEFT  , UPPER ,         TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmlutu=: (LEFT  , UPPER ,         TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmlujn=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmluju=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmlucn=: (LEFT  , UPPER ,    CONJ_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmlucu=: (LEFT  , UPPER ,    CONJ_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmrlnn=: (RIGHT , LOWER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmrlnu=: (RIGHT , LOWER ,      NO_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmrltn=: (RIGHT , LOWER ,         TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmrltu=: (RIGHT , LOWER ,         TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmrljn=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmrlju=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmrlcn=: (RIGHT , LOWER ,    CONJ_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmrlcu=: (RIGHT , LOWER ,    CONJ_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmrunn=: (RIGHT , UPPER ,      NO_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmrunu=: (RIGHT , UPPER ,      NO_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmrutn=: (RIGHT , UPPER ,         TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmrutu=: (RIGHT , UPPER ,         TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmrujn=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmruju=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
ztrsmrucn=: (RIGHT , UPPER ,    CONJ_TRANSPOSE ,: NONUNIT_DIAG)&ztrsmcore
ztrsmrucu=: (RIGHT , UPPER ,    CONJ_TRANSPOSE ,:    UNIT_DIAG)&ztrsmcore
