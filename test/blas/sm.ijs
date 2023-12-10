NB. Solve matrix equation
NB.
NB. xtrsmxxxx  Solve matrix equation with triangular matrix
NB.
NB. Version: 0.14.0 2023-03-21
NB.
NB. Copyright 2010-2023 Igor Zhuravlov
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

coclass 'mtbla'

NB. =========================================================
NB. Includes

require 'math/mt/test/blas/blas'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Dyad         Domain
NB. dtrsmcore    real
NB. ztrsmcore    complex
NB.
NB. Description:
NB.   Solves the matrix equation:
NB.     op(A) * X = alpha * B  (1)
NB.   or
NB.     X * op(A) = alpha * B  (2)
NB.   with transposed matrices, where A is triangular, and
NB.   op(A) is either A, A^T or A^H
NB.
NB. Syntax:
NB.   Xt=. (side ; uplo ; trans ; diag) xtrsmcore alpha ; AAt ; Bt
NB. where
NB.   side  - string, case-insensitive, in which the head
NB.           specifies the side of op(A):
NB.             'L'  NB. (1): op(A) is on the left
NB.             'R'  NB. (2): op(A) is on the right
NB.   uplo  - string, case-insensitive, in which the head
NB.           specifies whether the matrix A is upper or
NB.           lower triangular:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   trans - string, case-insensitive, in which the head
NB.           specifies the form of op(A):
NB.             'N'  NB. op(A) := A    (no transpose)
NB.             'T'  NB. op(A) := A^T  (transpose)
NB.             'C'  NB. op(A) := A^T  (transpose)           for dtrsmcore
NB.                  NB.       := A^H  (conjugate transpose) for ztrsmcore
NB.   diag  - string, case-insensitive, in which the head
NB.            specifies the form of A:
NB.             'N'  NB. A is either L or U
NB.             'U'  NB. A is either L1 or U1, diagonal
NB.                  NB.   elements of A are not referenced
NB.   alpha - scalar
NB.   AAt   - k×k-matrix, contains either non-zero or both
NB.           part(s) of A^T
NB.   Bt    - n×m-matrix, B^T, non-scaled RHS
NB.   Xt    - n×m-matrix, X^T, solutions
NB.           transposed
NB.   A     - k×k-matrix, triangular
NB.   m     ≥ 0, the number of columns in Bt and Xt
NB.   n     ≥ 0, the number of rows in Bt and Xt
NB.   k     = m if side='L' or k = n otherwise
NB.
NB. Notes:
NB. - operate on transposed matrix to avoid transposition

dtrsmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  'side uplo trans diag'=. x
  'alpha At Bt'=. y
  'n m'=. $ Bt
  10 {:: dtrsmcd (, side) ; (, uplo) ; (, trans) ; (, diag) ; (, m) ; (, n) ; (, alpha) ; At ; (, 1 >. c At) ; Bt ; , 1 >. m
)

ztrsmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  'side uplo trans diag'=. x
  'alpha At Bt'=. y
  'n m'=. $ Bt
  10 {:: ztrsmcd (, side) ; (, uplo) ; (, trans) ; (, diag) ; (, m) ; (, n) ; (, alpha) ; At ; (, 1 >. c At) ; Bt ; , 1 >. m
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad        Domain     Side    A     Reads in A    op(A)
NB. dtrsmllnn    real       (1)     L      LT           A
NB. dtrsmllnu    real       (1)     L1    SLT           A
NB. dtrsmlltn    real       (1)     L      LT           A^T
NB. dtrsmlltu    real       (1)     L1    SLT           A^T
NB. dtrsmlunn    real       (1)     U      UT           A
NB. dtrsmlunu    real       (1)     U1    SUT           A
NB. dtrsmlutn    real       (1)     U      UT           A^T
NB. dtrsmlutu    real       (1)     U1    SUT           A^T
NB. dtrsmrlnn    real       (2)     L      LT           A
NB. dtrsmrlnu    real       (2)     L1    SLT           A
NB. dtrsmrltn    real       (2)     L      LT           A^T
NB. dtrsmrltu    real       (2)     L1    SLT           A^T
NB. dtrsmrunn    real       (2)     U      UT           A
NB. dtrsmrunu    real       (2)     U1    SUT           A
NB. dtrsmrutn    real       (2)     U      UT           A^T
NB. dtrsmrutu    real       (2)     U1    SUT           A^T
NB. ztrsmllnn    complex    (1)     L      LT           A
NB. ztrsmllnu    complex    (1)     L1    SLT           A
NB. ztrsmlltn    complex    (1)     L      LT           A^T
NB. ztrsmlltu    complex    (1)     L1    SLT           A^T
NB. ztrsmllcn    complex    (1)     L      LT           A^H
NB. ztrsmllcu    complex    (1)     L1    SLT           A^H
NB. ztrsmlunn    complex    (1)     U      UT           A
NB. ztrsmlunu    complex    (1)     U1    SUT           A
NB. ztrsmlutn    complex    (1)     U      UT           A^T
NB. ztrsmlutu    complex    (1)     U1    SUT           A^T
NB. ztrsmlucn    complex    (1)     U      UT           A^H
NB. ztrsmlucu    complex    (1)     U1    SUT           A^H
NB. ztrsmrlnn    complex    (2)     L      LT           A
NB. ztrsmrlnu    complex    (2)     L1    SLT           A
NB. ztrsmrltn    complex    (2)     L      LT           A^T
NB. ztrsmrltu    complex    (2)     L1    SLT           A^T
NB. ztrsmrlcn    complex    (2)     L      LT           A^H
NB. ztrsmrlcu    complex    (2)     L1    SLT           A^H
NB. ztrsmrunn    complex    (2)     U      UT           A
NB. ztrsmrunu    complex    (2)     U1    SUT           A
NB. ztrsmrutn    complex    (2)     U      UT           A^T
NB. ztrsmrutu    complex    (2)     U1    SUT           A^T
NB. ztrsmrucn    complex    (2)     U      UT           A^H
NB. ztrsmrucu    complex    (2)     U1    SUT           A^H
NB.
NB. Description:
NB.   Solves the matrix equation:
NB.     op(A) * X = alpha * B  (1)
NB.   or
NB.     X * op(A) = alpha * B  (2)
NB.   where A is triangular, and op(A) is either A, A^T or
NB.   A^H
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
NB. - monad        provides BLAS'
NB.   dtrsmllnn    DTRSM('L','L','N','N',...)
NB.   dtrsmllnu    DTRSM('L','L','N','U',...)
NB.   dtrsmlltn    DTRSM('L','L','T','N',...), DTRSM('L','L','C','N',...)
NB.   dtrsmlltu    DTRSM('L','L','T','U',...), DTRSM('L','L','C','U',...)
NB.   dtrsmlunn    DTRSM('L','U','N','N',...)
NB.   dtrsmlunu    DTRSM('L','U','N','U',...)
NB.   dtrsmlutn    DTRSM('L','U','T','N',...), DTRSM('L','U','C','N',...)
NB.   dtrsmlutu    DTRSM('L','U','T','U',...), DTRSM('L','U','C','U',...)
NB.   dtrsmrlnn    DTRSM('R','L','N','N',...)
NB.   dtrsmrlnu    DTRSM('R','L','N','U',...)
NB.   dtrsmrltn    DTRSM('R','L','T','N',...), DTRSM('R','L','C','N',...)
NB.   dtrsmrltu    DTRSM('R','L','T','U',...), DTRSM('R','L','C','U',...)
NB.   dtrsmrunn    DTRSM('R','U','N','N',...)
NB.   dtrsmrunu    DTRSM('R','U','N','U',...)
NB.   dtrsmrutn    DTRSM('R','U','T','N',...), DTRSM('R','U','C','N',...)
NB.   dtrsmrutu    DTRSM('R','U','T','U',...), DTRSM('R','U','C','U',...)
NB.   ztrsmllnn    ZTRSM('L','L','N','N',...)
NB.   ztrsmllnu    ZTRSM('L','L','N','U',...)
NB.   ztrsmlltn    ZTRSM('L','L','T','N',...)
NB.   ztrsmlltu    ZTRSM('L','L','T','U',...)
NB.   ztrsmllcn    ZTRSM('L','L','C','N',...)
NB.   ztrsmllcu    ZTRSM('L','L','C','U',...)
NB.   ztrsmlunn    ZTRSM('L','U','N','N',...)
NB.   ztrsmlunu    ZTRSM('L','U','N','U',...)
NB.   ztrsmlutn    ZTRSM('L','U','T','N',...)
NB.   ztrsmlutu    ZTRSM('L','U','T','U',...)
NB.   ztrsmlucn    ZTRSM('L','U','C','N',...)
NB.   ztrsmlucu    ZTRSM('L','U','C','U',...)
NB.   ztrsmrlnn    ZTRSM('R','L','N','N',...)
NB.   ztrsmrlnu    ZTRSM('R','L','N','U',...)
NB.   ztrsmrltn    ZTRSM('R','L','T','N',...)
NB.   ztrsmrltu    ZTRSM('R','L','T','U',...)
NB.   ztrsmrlcn    ZTRSM('R','L','C','N',...)
NB.   ztrsmrlcu    ZTRSM('R','L','C','U',...)
NB.   ztrsmrunn    ZTRSM('R','U','N','N',...)
NB.   ztrsmrunu    ZTRSM('R','U','N','U',...)
NB.   ztrsmrutn    ZTRSM('R','U','T','N',...)
NB.   ztrsmrutu    ZTRSM('R','U','T','U',...)
NB.   ztrsmrucn    ZTRSM('R','U','C','N',...)
NB.   ztrsmrucu    ZTRSM('R','U','C','U',...)

dtrsmllnn=: 'runn'&dtrsmcore
dtrsmllnu=: 'runu'&dtrsmcore
dtrsmlltn=: 'rutn'&dtrsmcore
dtrsmlltu=: 'rutu'&dtrsmcore
dtrsmlunn=: 'rlnn'&dtrsmcore
dtrsmlunu=: 'rlnu'&dtrsmcore
dtrsmlutn=: 'rltn'&dtrsmcore
dtrsmlutu=: 'rltu'&dtrsmcore
dtrsmrlnn=: 'lunn'&dtrsmcore
dtrsmrlnu=: 'lunu'&dtrsmcore
dtrsmrltn=: 'lutn'&dtrsmcore
dtrsmrltu=: 'lutu'&dtrsmcore
dtrsmrunn=: 'llnn'&dtrsmcore
dtrsmrunu=: 'llnu'&dtrsmcore
dtrsmrutn=: 'lltn'&dtrsmcore
dtrsmrutu=: 'lltu'&dtrsmcore

ztrsmllnn=: 'runn'&ztrsmcore
ztrsmllnu=: 'runu'&ztrsmcore
ztrsmlltn=: 'rutn'&ztrsmcore
ztrsmlltu=: 'rutu'&ztrsmcore
ztrsmllcn=: 'rucn'&ztrsmcore
ztrsmllcu=: 'rucu'&ztrsmcore
ztrsmlunn=: 'rlnn'&ztrsmcore
ztrsmlunu=: 'rlnu'&ztrsmcore
ztrsmlutn=: 'rltn'&ztrsmcore
ztrsmlutu=: 'rltu'&ztrsmcore
ztrsmlucn=: 'rlcn'&ztrsmcore
ztrsmlucu=: 'rlcu'&ztrsmcore
ztrsmrlnn=: 'lunn'&ztrsmcore
ztrsmrlnu=: 'lunu'&ztrsmcore
ztrsmrltn=: 'lutn'&ztrsmcore
ztrsmrltu=: 'lutu'&ztrsmcore
ztrsmrlcn=: 'lucn'&ztrsmcore
ztrsmrlcu=: 'lucu'&ztrsmcore
ztrsmrunn=: 'llnn'&ztrsmcore
ztrsmrunu=: 'llnu'&ztrsmcore
ztrsmrutn=: 'lltn'&ztrsmcore
ztrsmrutu=: 'lltu'&ztrsmcore
ztrsmrucn=: 'llcn'&ztrsmcore
ztrsmrucu=: 'llcu'&ztrsmcore
