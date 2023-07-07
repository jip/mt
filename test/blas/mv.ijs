NB. Matrix-vector operations
NB.
NB. xgemvx    Matrix-vector operation with general matrix
NB. dsymvx    Matrix-vector operation with symmetric matrix
NB. zhemvx    Matrix-vector operation with Hermitian matrix
NB. xtrmvxxx  Matrix-vector operation with triangular matrix
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
NB. dgemvcore
NB. zgemvcore
NB.
NB. Description:
NB.   Performs the matrix-vector operation:
NB.     y := alpha * op(A) * x + beta * y
NB.   with transposed matrix
NB.
NB. Syntax:
NB.   yupd=. trans xgemvcore alpha ; At ; x ; incx ; beta ; y ; incy
NB. where
NB.   trans - literal, case-insensitive, in which the head
NB.           specifies the form of op(A):
NB.             'N'  NB. op(A) := A    (no transpose)
NB.             'T'  NB. op(A) := A^T  (transpose)
NB.             'C'  NB. op(A) := A^T  (transpose)           for dgemvcore
NB.                  NB.       := A^H  (conjugate transpose) for zgemvcore
NB.   alpha - scalar
NB.   At    - n×m-matrix, A^T
NB.   x     - (1+(kx-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   beta  - scalar
NB.   y     - (1+(ky-1)*|incy|)-vector
NB.   incy  ≠ 0, the increment for the elements of y
NB.   yupd  - an updated y
NB.   kx    = n if trans='N' or kx = m otherwise
NB.   ky    = m if trans='N' or ky = n otherwise
NB.   m     ≥ 0, the number of rows in A
NB.   n     ≥ 0, the number of columns in A
NB.
NB. Notes:
NB. - operate on transposed matrix to avoid transposition

dgemvcore=: (4 : 0) ([ assert@basiccr6)
  'alpha At xx incx beta y incy'=. y
  'n m'=. $ At
  10 {:: dgemvcd (, x) ; (, m) ; (, n) ; (, alpha) ; At ; (, 1 >. m) ; xx ; (, incx) ; (, beta) ; y ; , incy
)

zgemvcore=: (4 : 0) ([ assert@basiccr6)
  'alpha At xx incx beta y incy'=. y
  'n m'=. $ At
  10 {:: zgemvcd (, x) ; (, m) ; (, n) ; (, alpha) ; At ; (, 1 >. m) ; xx ; (, incx) ; (, beta) ; y ; , incy
)

NB. ---------------------------------------------------------
NB. dsymvcore
NB. zhemvcore
NB.
NB. Description:
NB.   Performs the matrix-vector operation:
NB.     y := alpha * A * x + beta * y
NB.   with transposed matrix, where A is Hermitian
NB.   (symmetric)
NB.
NB. Syntax:
NB.   yupd=. uplo xxxmvcore alpha ; At ; x ; incx ; beta ; y ; incy
NB. where
NB.   uplo  - literal, case-insensitive, in which the head
NB.           specifies which triangular part of A is to be
NB.           referenced:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   alpha - scalar
NB.   At    - n×n-matrix with real diagonal, A^T
NB.   x     - (1+(n-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   beta  - scalar
NB.   y     - (1+(n-1)*|incy|)-vector
NB.   incy  ≠ 0, the increment for the elements of y
NB.   yupd  - an updated y
NB.   n     ≥ 0, the size of A
NB.
NB. Notes:
NB. - operates on transposed matrix to avoid transposition

dsymvcore=: (4 : 0) ([ assert@(basiccs1 , basiccr6))
  'alpha At xx incx beta y incy'=. y
  n=. # At
  9 {:: dsymvcd (, x) ; (, n) ; (, alpha) ; At ; (, 1 >. n) ; xx ; (, incx) ; (, beta) ; y ; , incy
)

zhemvcore=: (4 : 0) ([ assert@(basiccs1 , basiccr6))
  'alpha At xx incx beta y incy'=. y
  n=. # At
  9 {:: zhemvcd (, x) ; (, n) ; (, alpha) ; At ; (, 1 >. n) ; xx ; (, incx) ; (, beta) ; y ; , incy
)

NB. ---------------------------------------------------------
NB. dtrmvcore
NB. ztrmvcore
NB.
NB. Description:
NB.   Performs the matrix-vector operation:
NB.     x := op(A) * x
NB.   with transposed matrix, where A is triangular
NB.
NB. Syntax:
NB.   xupd=. (uplo ; trans ; diag) xtrmvcore At ; x ; incx
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
NB.             'C'  NB. op(A) := A^T  (transpose)           for dtrmvcore
NB.                  NB.       := A^H  (conjugate transpose) for ztrmvcore
NB.   diag  - literal, case-insensitive, in which the head
NB.           specifies the form of A:
NB.             'N'  NB. A is either L or U
NB.             'U'  NB. A is either L1 or U1, diagonal
NB.                  NB.   elements of A are not referenced
NB.   At    - n×n-matrix, A^T
NB.   x     - (1+(n-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   xupd  - an updated x
NB.   n     ≥ 0, the size of A
NB.
NB. Notes:
NB. - operate on transposed matrix to avoid transposition

dtrmvcore=: (4 : 0) ([ assert@(basiccs0 , basiccr1))
  'uplo trans diag'=. x
  'At y incy'=. y
  n=. # At
  7 {:: dtrmvcd (, uplo) ; (, trans) ; (, diag) ; (, n) ; At ; (, 1 >. n) ; y ; , incy
)

ztrmvcore=: (4 : 0) ([ assert@(basiccs0 , basiccr1))
  'uplo trans diag'=. x
  'At y incy'=. y
  n=. # At
  7 {:: ztrmvcd (, uplo) ; (, trans) ; (, diag) ; (, n) ; At ; (, 1 >. n) ; y ; , incy
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad     Operation
NB. dgemvn    C := alpha * A   * x + beta * y
NB. dgemvt    C := alpha * A^T * x + beta * y
NB. zgemvn    C := alpha * A   * x + beta * y
NB. zgemvt    C := alpha * A^T * x + beta * y
NB. zgemvc    C := alpha * A^H * x + beta * y
NB.
NB. Description:
NB.   Performs the matrix-vector operation:
NB.     y := alpha * op(A) * x + beta * y
NB.
NB. Syntax:
NB.   yupd=. xgemvx alpha ; A ; x ; incx ; beta ; y ; incy
NB. where
NB.   alpha - scalar
NB.   A     - m×n-matrix
NB.   x     - (1+(kx-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   beta  - scalar
NB.   y     - (1+(ky-1)*|incy|)-vector
NB.   incy  ≠ 0, the increment for the elements of y
NB.   yupd  - an updated y
NB.   kx    = n for xgemvn or kx = m otherwise
NB.   ky    = m for xgemvn or ky = n otherwise
NB.   m     ≥ 0, the number of rows in A
NB.   n     ≥ 0, the number of columns in A
NB.
NB. Notes:
NB. - monad     provides BLAS'
NB.   dgemvn    DGEMV('N',...)
NB.   dgemvt    DGEMV('T',...), DGEMV('C',...)
NB.   zgemvn    ZGEMV('N',...)
NB.   zgemvt    ZGEMV('T',...)
NB.   zgemvc    ZGEMV('C',...)

dgemvn=: 't'&   dgemvcore
dgemvt=: 'n'&   dgemvcore

zgemvn=: 't'&   zgemvcore
zgemvt=: 'n'&   zgemvcore
zgemvc=: 'n'&(+@zgemvcore basiccj2)

NB. ---------------------------------------------------------
NB. Monad     A            Reads in A
NB. dsymvl    symmetric    LT
NB. dsymvu    symmetric    UT
NB. zhemvl    Hermitian    LT
NB. zhemvu    Hermitian    UT
NB.
NB. Description:
NB.   Performs the matrix-vector operation:
NB.     y := alpha * A * x + beta * y
NB.   where A is Hermitian (symmetric)
NB.
NB. Syntax:
NB.   yupd=. xxxmvx alpha ; A ; x ; incx ; beta ; y ; incy
NB. where
NB.   alpha - scalar
NB.   A     - n×n-matrix with real diagonal
NB.   x     - (1+(n-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   beta  - scalar
NB.   y     - (1+(n-1)*|incy|)-vector
NB.   incy  ≠ 0, the increment for the elements of y
NB.   yupd  - an updated y
NB.   n     ≥ 0, the size of A
NB.
NB. Notes:
NB. - monad     provides BLAS'
NB.   dsymvl    DSYMV('L',...)
NB.   dsymvu    DSYMV('U',...)
NB.   zhemvl    ZHEMV('L',...)
NB.   zhemvu    ZHEMV('U',...)

dsymvl=: 'u'&dsymvcore
dsymvu=: 'l'&dsymvcore

zhemvl=: 'u'&zhemvcore
zhemvu=: 'l'&zhemvcore

NB. ---------------------------------------------------------
NB. Monad       A     Reads in A    Operation
NB. dtrmvlnn    L      LT           x := A   * x
NB. dtrmvlnu    L1    SLT           x := A   * x
NB. dtrmvltn    L      LT           x := A^T * x
NB. dtrmvltu    L1    SLT           x := A^T * x
NB. dtrmvunn    U      UT           x := A   * x
NB. dtrmvunu    U1    SUT           x := A   * x
NB. dtrmvutn    U      UT           x := A^T * x
NB. dtrmvutu    U1    SUT           x := A^T * x
NB. ztrmvlnn    L      LT           x := A   * x
NB. ztrmvlnu    L1    SLT           x := A   * x
NB. ztrmvltn    L      LT           x := A^T * x
NB. ztrmvltu    L1    SLT           x := A^T * x
NB. ztrmvlcn    L      LT           x := A^H * x
NB. ztrmvlcu    L1    SLT           x := A^H * x
NB. ztrmvunn    U      UT           x := A   * x
NB. ztrmvunu    U1    SUT           x := A   * x
NB. ztrmvutn    U      UT           x := A^T * x
NB. ztrmvutu    U1    SUT           x := A^T * x
NB. ztrmvucn    U      UT           x := A^H * x
NB. ztrmvucu    U1    SUT           x := A^H * x
NB.
NB. Description:
NB.   Performs the matrix-vector operation:
NB.     x := op(A) * x
NB.   where A is triangular
NB.
NB. Syntax:
NB.   xupd=. xtrmvxxx A ; x ; incx
NB. where
NB.   A    - n×n-matrix, contains either L, L1, U or U1
NB.          (unit diagonal is not stored)
NB.   x    - (1+(n-1)*|incx|)-vector
NB.   incx ≠ 0, the increment for the elements of x
NB.   xupd - an updated x
NB.   n    ≥ 0, the size of A
NB.
NB. Notes:
NB. - monad       provides BLAS'
NB.   dtrmvlnn    DTRMV('L','N','N',...)
NB.   dtrmvlnu    DTRMV('L','N','U',...)
NB.   dtrmvltn    DTRMV('L','T','N',...), DTRMV('L','C','N',...)
NB.   dtrmvltu    DTRMV('L','T','U',...), DTRMV('L','C','U',...)
NB.   dtrmvunn    DTRMV('U','N','N',...)
NB.   dtrmvunu    DTRMV('U','N','U',...)
NB.   dtrmvutn    DTRMV('U','T','N',...), DTRMV('U','C','N',...)
NB.   dtrmvutu    DTRMV('U','T','U',...), DTRMV('U','C','U',...)
NB.   ztrmvlnn    ZTRMV('L','N','N',...)
NB.   ztrmvlnu    ZTRMV('L','N','U',...)
NB.   ztrmvltn    ZTRMV('L','T','N',...)
NB.   ztrmvltu    ZTRMV('L','T','U',...)
NB.   ztrmvlcn    ZTRMV('L','C','N',...)
NB.   ztrmvlcu    ZTRMV('L','C','U',...)
NB.   ztrmvunn    ZTRMV('U','N','N',...)
NB.   ztrmvunu    ZTRMV('U','N','U',...)
NB.   ztrmvutn    ZTRMV('U','T','N',...)
NB.   ztrmvutu    ZTRMV('U','T','U',...)
NB.   ztrmvucn    ZTRMV('U','C','N',...)
NB.   ztrmvucu    ZTRMV('U','C','U',...)

dtrmvlnn=: 'utn'&   dtrmvcore
dtrmvlnu=: 'utu'&   dtrmvcore
dtrmvltn=: 'unn'&   dtrmvcore
dtrmvltu=: 'unu'&   dtrmvcore
dtrmvunn=: 'ltn'&   dtrmvcore
dtrmvunu=: 'ltu'&   dtrmvcore
dtrmvutn=: 'lnn'&   dtrmvcore
dtrmvutu=: 'lnu'&   dtrmvcore

ztrmvlnn=: 'utn'&   ztrmvcore
ztrmvlnu=: 'utu'&   ztrmvcore
ztrmvltn=: 'unn'&   ztrmvcore
ztrmvltu=: 'unu'&   ztrmvcore
ztrmvlcn=: 'unn'&(+@ztrmvcore basiccj0)
ztrmvlcu=: 'unu'&(+@ztrmvcore basiccj0)
ztrmvunn=: 'ltn'&   ztrmvcore
ztrmvunu=: 'ltu'&   ztrmvcore
ztrmvutn=: 'lnn'&   ztrmvcore
ztrmvutu=: 'lnu'&   ztrmvcore
ztrmvucn=: 'lnn'&(+@ztrmvcore basiccj0)
ztrmvucu=: 'lnu'&(+@ztrmvcore basiccj0)
