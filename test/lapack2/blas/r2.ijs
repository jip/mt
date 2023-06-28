NB. Rank 2 operations
NB.
NB. dsyr2x  Rank 2 operation with symmetric matrix
NB. zher2x  Rank 2 operation with Hermitian matrix
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
NB. Dyad         A            Operation
NB. dsyr2core    symmetric    A := alpha * x * y^T +      alpha  * y * x^T + A
NB. zher2core    Hermitian    A := alpha * x * y^H + conj(alpha) * y * x^H + A
NB.
NB. Description:
NB.   Performs the hermitian (symmetric) rank 2 operation
NB.   with transposed matrix, where A is Hermitian
NB.   (symmetric)
NB.
NB. Syntax:
NB.   Aupdt=. uplo xxxr2core alpha ; x ; incx ; y ; incy ; At
NB. where
NB.   uplo  - literal, case-insensitive, in which the head
NB.           specifies which triangular part of A is to be
NB.           referenced:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   alpha - scalar
NB.   x     - (1+(n-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   y     - (1+(n-1)*|incy|)-vector
NB.   incy  ≠ 0, the increment for the elements of y
NB.   At    - n×n-matrix with real diagonal, A^T
NB.   Aupdt - At with either LT (if uplo='U') or UT (if
NB.           uplo='L') updated
NB.   n     ≥ 0, the size of A and Aupdt
NB.
NB. Notes:
NB. - operate on transposed matrices to avoid transposition

dsyr2core=: (4 : 0) ([ assert@(basiccs5 , basiccr5))
  'alpha xx incx y incy At'=. y
  n=. # At
  8 {:: dsyr2cd (, x) ; (, n) ; (, alpha) ; xx ; (, incx) ; y ; (, incy) ; At ; , 1 >. n
)

zher2core=: (4 : 0) ([ assert@(basiccs5 , basiccr5))
  'alpha xx incx y incy At'=. y
  n=. # At
  8 {:: zher2cd (, x) ; (, n) ; (, alpha) ; xx ; (, incx) ; y ; (, incy) ; At ; , 1 >. n
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad     A            R/W in A    Operation
NB. dsyr2l    symmetric    LT          A := alpha * x * y^T +      alpha  * y * x^T + A
NB. dsyr2u    symmetric    UT          A := alpha * x * y^T +      alpha  * y * x^T + A
NB. zher2l    Hermitian    LT          A := alpha * x * y^H + conj(alpha) * y * x^H + A
NB. zher2u    Hermitian    UT          A := alpha * x * y^H + conj(alpha) * y * x^H + A
NB.
NB. Description:
NB.   Performs the hermitian (symmetric) rank 2 operation
NB.   where A is Hermitian (symmetric)
NB.
NB. Syntax:
NB.   Aupd=. xxxr2x alpha ; x ; incx ; y ; incy ; A
NB. where
NB.   alpha - scalar
NB.   x     - (1+(n-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   y     - (1+(n-1)*|incy|)-vector
NB.   incy  ≠ 0, the increment for the elements of y
NB.   A     - n×n-matrix with real diagonal
NB.   Aupd  - A with either LT (for xxxr2l) or UT (for
NB.           xxxr2u) updated
NB.   n     ≥ 0, the size of A and Aupd
NB.
NB. Notes:
NB. - dsyr2l provides BLAS' DSYR2('L',...)
NB. - dsyr2u provides BLAS' DSYR2('U',...)
NB. - zher2l provides BLAS' ZHER2('L',...)
NB. - zher2u provides BLAS' ZHER2('U',...)

dsyr2l=: 'u'& dsyr2core
dsyr2u=: 'l'& dsyr2core

zher2l=: 'u'&(zher2core basiccj1)
zher2u=: 'l'&(zher2core basiccj1)
