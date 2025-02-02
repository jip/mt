NB. Solve equation
NB.
NB. xtrsvxxx  Solve equation with triangular matrix
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

coclass 'mtbla'

NB. =========================================================
NB. Includes

require 'math/mt/external/blas/blas'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Dyad         Domain
NB. dtrsvcore    real
NB. ztrsvcore    complex
NB.
NB. Description:
NB.   Solves the equation:
NB.     op(A) * x = b
NB.   with transposed matrix, where A is triangular, and
NB.   op(A) is either A, A^T or A^H
NB.
NB. Syntax:
NB.   x=. (uplo ; trans ; diag) xtrsvcore AAt ; b ; incb
NB. where
NB.   uplo  - string, case-insensitive, in which the head
NB.           specifies whether the matrix A is upper or
NB.           lower triangular:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   trans - string, case-insensitive, in which the head
NB.           specifies the form of op(A):
NB.             'N'  NB. op(A) := A    (no transpose)
NB.             'T'  NB. op(A) := A^T  (transpose)
NB.             'C'  NB. op(A) := A^T  (transpose)           for dtrsvcore
NB.                  NB.       := A^H  (conjugate transpose) for ztrsvcore
NB.   diag  - string, case-insensitive, in which the head
NB.            specifies the form of A:
NB.             'N'  NB. A is either L or U
NB.             'U'  NB. A is either L1 or U1, diagonal
NB.                  NB.   elements of A are not referenced
NB.   AAt   - n×n-matrix, contains either non-zero or both
NB.           part(s) of A^T
NB.   b     - (1+(n-1)*|incx|)-vector, the RHS
NB.   incb  ≠ 0, the increment for the elements of b and x
NB.   x     - the same shape as b, the solution
NB.   A     - n×n-matrix, triangular
NB.   n     ≥ 0, the size of A and AAt
NB.
NB. Notes:
NB. - operate on transposed matrix to avoid transposition

dtrsvcore=: (4 : 0) ([ assert@(basiccs0 , basiccr1))
  'uplo trans diag'=. x
  'At b incb'=. y
  n=. # At
  7 {:: dtrsv_cd (, uplo) ; (, trans) ; (, diag) ; (, n) ; At ; (, 1 >. n) ; b ; , incb
)

ztrsvcore=: (4 : 0) ([ assert@(basiccs0 , basiccr1))
  'uplo trans diag'=. x
  'At b incb'=. y
  n=. # At
  7 {:: ztrsv_cd (, uplo) ; (, trans) ; (, diag) ; (, n) ; At ; (, 1 >. n) ; b ; , incb
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad       Domain     A     Reads in A    op(A)
NB. dtrsvlnn    real       L      LT           A
NB. dtrsvlnu    real       L1    SLT           A
NB. dtrsvltn    real       L      LT           A^T
NB. dtrsvltu    real       L1    SLT           A^T
NB. dtrsvunn    real       U      UT           A
NB. dtrsvunu    real       U1    SUT           A
NB. dtrsvutn    real       U      UT           A^T
NB. dtrsvutu    real       U1    SUT           A^T
NB. ztrsvlnn    complex    L      LT           A
NB. ztrsvlnu    complex    L1    SLT           A
NB. ztrsvltn    complex    L      LT           A^T
NB. ztrsvltu    complex    L1    SLT           A^T
NB. ztrsvlcn    complex    L      LT           A^H
NB. ztrsvlcu    complex    L1    SLT           A^H
NB. ztrsvunn    complex    U      UT           A
NB. ztrsvunu    complex    U1    SUT           A
NB. ztrsvutn    complex    U      UT           A^T
NB. ztrsvutu    complex    U1    SUT           A^T
NB. ztrsvucn    complex    U      UT           A^H
NB. ztrsvucu    complex    U1    SUT           A^H
NB.
NB. Description:
NB.   Solves the equation:
NB.     op(A) * x = b
NB.   where A is triangular, and op(A) is either A, A^T or
NB.   A^H
NB.
NB. Syntax:
NB.   x=. xtrsvxxx AA ; b ; incb
NB. where
NB.   AA   - n×n-matrix, contains either non-zero or both
NB.          part(s) of A
NB.   b    - (1+(n-1)*|incb|)-vector, the RHS
NB.   incb ≠ 0, the increment for the elements of b and x
NB.   x    - the same shape as b, the solution
NB.   A    - n×n-matrix, triangular
NB.   n    ≥ 0, the size of A and AA
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
