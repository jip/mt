NB. Rank 1 operations
NB.
NB. xgerx  Rank 1 operation with general matrix
NB. dsyrx  Rank 1 operation with symmetric matrix
NB. zherx  Rank 1 operation with Hermitian matrix
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
NB. Dyad        Domain     A     op(x)
NB. dsyrcore    real       SY    x^T
NB. zhercore    complex    HE    x^H
NB.
NB. Description:
NB.   Performs the hermitian (symmetric) rank 1 operation:
NB.     A := alpha * x * op(x) + A
NB.   with transposed matrix, where A is Hermitian
NB.   (symmetric)
NB.
NB. Syntax:
NB.   AAupdt=. uplo xxxrcore alpha ; x ; incx ; AAt
NB. where
NB.   uplo  - string, case-insensitive, in which the head
NB.           specifies which triangular part of A is to be
NB.           referenced:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   alpha  - scalar, real
NB.   x      - (1+(n-1)*|incx|)-vector
NB.   incx   ≠ 0, the increment for the elements of x
NB.   AAt    - n×n-matrix, contains either LT or UT or both
NB.            part(s) of A^T
NB.   AAupdt - AAt with either LT (if uplo='U') or UT (if
NB.            uplo='L') updated
NB.   A      - n×n-matrix, Hermitian (symmetric)
NB.   n      ≥ 0, the size of A, AAt and AAupdt
NB.
NB. Notes:
NB. - operate on transposed matrices to avoid transposition

dsyrcore=: (4 : 0) ([ assert@(basiccs3 , basiccr2))
  'alpha y incy AAt'=. y
  n=. # AAt
  6 {:: dsyrcd (, x) ; (, n) ; (, alpha) ; y ; (, incy) ; AAt ; , 1 >. n
)

zhercore=: (4 : 0) ([ assert@(basiccs3 , basiccr2))
  'alpha y incy AAt'=. y
  n=. # AAt
  6 {:: zhercd (, x) ; (, n) ; (, alpha) ; y ; (, incy) ; AAt ; , 1 >. n
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad    Domain     op(y)
NB. dger     real       y^T
NB. zgerc    complex    y^H
NB. zgeru    complex    y^T
NB.
NB. Description:
NB.   Performs the rank 1 operation:
NB.     A := alpha * x * op(y) + A
NB.
NB. Syntax:
NB.   Aupd=. xgerx alpha ; x ; incx ; y ; incy ; A
NB. where
NB.   alpha - scalar
NB.   x     - (1+(m-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   y     - (1+(n-1)*|incy|)-vector
NB.   incy  ≠ 0, the increment for the elements of y
NB.   A     - m×n-matrix
NB.   Aupd  - an updated A
NB.   m     ≥ 0, the number of rows in A and Aupd
NB.   n     ≥ 0, the number of columns in A in Aupd
NB.
NB. Notes:
NB. - dger  provides BLAS' DGER(...)
NB. - zgerc provides BLAS' ZGERC(...)
NB. - zgeru provides BLAS' ZGERU(...)

dger=: (3 : 0)@([ assert@basiccr5)
  'alpha x incx y incy A'=. y
  'm n'=. $ A
  assert (0 >. >: (<: m) * | incx) = # x
  assert (0 >. >: (<: n) * | incy) = # y
  8 {:: dgercd (, n) ; (, m) ; (, alpha) ;    y  ; (, incy) ;    x  ; (, incx) ; A ; , 1 >. n
)

zgerc=: (3 : 0)@([ assert@basiccr5)
  'alpha x incx y incy A'=. y
  'm n'=. $ A
  assert (0 >. >: (<: m) * | incx) = # x
  assert (0 >. >: (<: n) * | incy) = # y
  8 {:: zgerccd (, n) ; (, m) ; (, alpha) ; (+ y) ; (, incy) ; (+ x) ; (, incx) ; A ; , 1 >. n
)

zgeru=: (3 : 0)@([ assert@basiccr5)
  'alpha x incx y incy A'=. y
  'm n'=. $ A
  assert (0 >. >: (<: m) * | incx) = # x
  assert (0 >. >: (<: n) * | incy) = # y
  8 {:: zgerucd (, n) ; (, m) ; (, alpha) ;    y  ; (, incy) ;    x  ; (, incx) ; A ; , 1 >. n
)

NB. ---------------------------------------------------------
NB. Monad    Domain     A     R/W in A    op(x)
NB. dsyrl    real       SY    LT          x^T
NB. dsyru    real       SY    UT          x^T
NB. zherl    complex    HE    LT          x^H
NB. zheru    complex    HE    UT          x^H
NB.
NB. Description:
NB.   Performs the hermitian (symmetric) rank 1 operation:
NB.     A := alpha * x * op(x) + A
NB.   where A is Hermitian (symmetric)
NB.
NB. Syntax:
NB.   AAupd=. xxxrx alpha ; x ; incx ; AA
NB. where
NB.   alpha - scalar, real
NB.   x     - (1+(n-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   AA    - n×n-matrix, contains either LT or UT or both
NB.           part(s) of A
NB.   AAupd - AA with either LT (for xxxrl) or UT (for xxxru)
NB.           updated
NB.   A     - n×n-matrix, Hermitian (symmetric)
NB.   n     ≥ 0, the size of A, AA and AAupd
NB.
NB. Notes:
NB. - monad    provides BLAS'
NB.   dsyrl    DSYR('L',...)
NB.   dsyru    DSYR('U',...)
NB.   zherl    ZHER('L',...)
NB.   zheru    ZHER('U',...)

dsyrl=: 'u'& dsyrcore
dsyru=: 'l'& dsyrcore

zherl=: 'u'&(zhercore basiccj0)
zheru=: 'l'&(zhercore basiccj0)
