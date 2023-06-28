NB. Solve equation
NB.
NB. xtrsvxxx  Solve equation with triangular matrix
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

require 'math/mt/test/lapack2/blas/blas'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. dtrsvcore
NB. ztrsvcore
NB.
NB. Description:
NB.   Solves the equation:
NB.     op(A) * x = b
NB.   with transposed matrix, where A is triangular
NB.
NB. Syntax:
NB.   x=. (uplo ; trans ; diag) xtrsvcore At ; b ; incb
NB. where
NB.   uplo  - literal, case-insensitive, in which the head
NB.           specifies whether the matrix A is upper or
NB.           lower triangular:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   trans - literal, case-insensitive, in which the head
NB.           specifies the form of op(A):
NB.             'N'  NB. op(A) := A    (no transpose)
NB.             'T'  NB. op(A) := A^T  (transpose)
NB.             'C'  NB. op(A) := A^T  (transpose)           for dtrsvcore
NB.                  NB.       := A^H  (conjugate transpose) for ztrsvcore
NB.   diag  - literal, case-insensitive, in which the head
NB.            specifies the form of A:
NB.             'N'  NB. A is either L or U
NB.             'U'  NB. A is either L1 or U1, diagonal
NB.                  NB.   elements of A are not referenced
NB.   At    - n×n-matrix, A^T
NB.   b     - (1+(n-1)*|incx|)-vector, the RHS
NB.   incb  ≠ 0, the increment for the elements of b and x
NB.   x     - the same shape as b, the solution
NB.   n     ≥ 0, the size of matrix A and vectors b and x
NB.
NB. Notes:
NB. - operate on transposed matrix to avoid transposition

dtrsvcore=: (4 : 0) ([ assert@(basiccs0 , basiccr1))
  'uplo trans diag'=. x
  'At b incb'=. y
  n=. # At
  7 {:: dtrsvcd (, uplo) ; (, trans) ; (, diag) ; (, n) ; At ; (, 1 >. n) ; b ; , incb
)

ztrsvcore=: (4 : 0) ([ assert@(basiccs0 , basiccr1))
  'uplo trans diag'=. x
  'At b incb'=. y
  n=. # At
  7 {:: ztrsvcd (, uplo) ; (, trans) ; (, diag) ; (, n) ; At ; (, 1 >. n) ; b ; , incb
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad       Reads in A    Solves
NB. dtrsvlnn     LT           L    * x = b
NB. dtrsvlnu    SLT           L1   * x = b
NB. dtrsvltn     LT           L ^T * x = b
NB. dtrsvltu    SLT           L1^T * x = b
NB. dtrsvunn     UT           U    * x = b
NB. dtrsvunu    SUT           U1   * x = b
NB. dtrsvutn     UT           U ^T * x = b
NB. dtrsvutu    SUT           U1^T * x = b
NB. ztrsvlnn     LT           L    * x = b
NB. ztrsvlnu    SLT           L1   * x = b
NB. ztrsvltn     LT           L ^T * x = b
NB. ztrsvltu    SLT           L1^T * x = b
NB. ztrsvlcn     LT           L ^H * x = b
NB. ztrsvlcu    SLT           L1^H * x = b
NB. ztrsvunn     UT           U    * x = b
NB. ztrsvunu    SUT           U1   * x = b
NB. ztrsvutn     UT           U ^T * x = b
NB. ztrsvutu    SUT           U1^T * x = b
NB. ztrsvucn     UT           U ^H * x = b
NB. ztrsvucu    SUT           U1^H * x = b
NB.
NB. Description:
NB.   Solves the equation:
NB.     op(A) * x = b
NB.   where A is triangular
NB.
NB. Syntax:
NB.   x=. xtrsvxxx A ; b ; incb
NB. where
NB.   A    - n×n-matrix, contains either L, L1, U or U1
NB.          (unit diagonal is not stored)
NB.   b    - (1+(n-1)*|incb|)-vector, the RHS
NB.   incb ≠ 0, the increment for the elements of b and x
NB.   x    - the same shape as b, the solution
NB.   n    ≥ 0, the size of matrix A, vectors b and x
NB.
NB. Notes:
NB. - monad       provides BLAS'
NB.   dtrsvlnn    DTRSV('L','N','N',...)
NB.   dtrsvlnu    DTRSV('L','N','U',...)
NB.   dtrsvltn    DTRSV('L','T','N',...), DTRSV('L','C','N',...)
NB.   dtrsvltu    DTRSV('L','T','U',...), DTRSV('L','C','U',...)
NB.   dtrsvunn    DTRSV('U','N','N',...)
NB.   dtrsvunu    DTRSV('U','N','U',...)
NB.   dtrsvutn    DTRSV('U','T','N',...), DTRSV('U','C','N',...)
NB.   dtrsvutu    DTRSV('U','T','U',...), DTRSV('U','C','U',...)
NB.   ztrsvlnn    ZTRSV('L','N','N',...)
NB.   ztrsvlnu    ZTRSV('L','N','U',...)
NB.   ztrsvltn    ZTRSV('L','T','N',...)
NB.   ztrsvltu    ZTRSV('L','T','U',...)
NB.   ztrsvlcn    ZTRSV('L','C','N',...)
NB.   ztrsvlcu    ZTRSV('L','C','U',...)
NB.   ztrsvunn    ZTRSV('U','N','N',...)
NB.   ztrsvunu    ZTRSV('U','N','U',...)
NB.   ztrsvutn    ZTRSV('U','T','N',...)
NB.   ztrsvutu    ZTRSV('U','T','U',...)
NB.   ztrsvucn    ZTRSV('U','C','N',...)
NB.   ztrsvucu    ZTRSV('U','C','U',...)

dtrsvlnn=: 'utn'&   dtrsvcore
dtrsvlnu=: 'utu'&   dtrsvcore
dtrsvltn=: 'unn'&   dtrsvcore
dtrsvltu=: 'unu'&   dtrsvcore
dtrsvunn=: 'ltn'&   dtrsvcore
dtrsvunu=: 'ltu'&   dtrsvcore
dtrsvutn=: 'lnn'&   dtrsvcore
dtrsvutu=: 'lnu'&   dtrsvcore

ztrsvlnn=: 'utn'&   ztrsvcore
ztrsvlnu=: 'utu'&   ztrsvcore
ztrsvltn=: 'unn'&   ztrsvcore
ztrsvltu=: 'unu'&   ztrsvcore
ztrsvlcn=: 'unn'&(+@ztrsvcore basiccj0)
ztrsvlcu=: 'unu'&(+@ztrsvcore basiccj0)
ztrsvunn=: 'ltn'&   ztrsvcore
ztrsvunu=: 'ltu'&   ztrsvcore
ztrsvutn=: 'lnn'&   ztrsvcore
ztrsvutu=: 'lnu'&   ztrsvcore
ztrsvucn=: 'lnn'&(+@ztrsvcore basiccj0)
ztrsvucu=: 'lnu'&(+@ztrsvcore basiccj0)
