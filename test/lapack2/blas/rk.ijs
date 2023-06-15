NB. Rank k operations
NB.
NB. xsyrkxx  Rank k operation with symmetric matrix
NB. zherkxx  Rank k operation with Hermitian matrix
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

coclass 'mttst'

NB. =========================================================
NB. Includes

require 'math/mt/test/lapack2/blas/blas'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. dsyrkcore
NB. zsyrkcore
NB.
NB. Description:
NB.   Performs symmetric rank k operations:
NB.     C := alpha * A * A^T + beta * C  (1)
NB.   or
NB.     C := alpha * A^T * A + beta * C  (2)
NB.   with transposed matrices, where C is symmetric
NB.
NB. Syntax:
NB.   Cupdt=. (uplo ; trans) xsyrkcore alpha ; At ; beta ; Ct
NB. where
NB.   uplo  - literal, case-insensitive, in which the head
NB.           specifies which triangular part of C is to be
NB.           referenced:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   trans - literal, case-insensitive, in which the head
NB.           specifies the operation:
NB.             'N'  NB. (1)
NB.             'T'  NB. (2)
NB.             'C'  NB. (2) for dsyrkcore
NB.   alpha - scalar
NB.   At    - ka×na-matrix, A^T
NB.   beta  - scalar
NB.   Ct    - n×n-matrix, C^T
NB.   Cupdt - Ct with either LT (if uplo='U') or UT (if
NB.           uplo='L') updated
NB.   n     ≥ 0, the size of Ct and Cupdt and the number of
NB.           columns or rows in At
NB.   k     ≥ 0, the number of rows or columns in At
NB.   ka    = k if trans='N' or ka = n otherwise
NB.   na    = n if trans='N' or na = k otherwise
NB.
NB. Notes:
NB. - operate on transposed matrices to avoid transposition

dsyrkcore=: (4 : 0) ([ assert@(basiccs3 , basiccr3))
  'uplo trans'=. x
  'alpha At beta Ct'=. y
  n=. # Ct
  k=. (-. 'nN' e.~ {. trans) { $ At
  9 {:: dsyrkcd (, uplo) ; (, trans) ; (, n) ; (, k) ; (, alpha) ; At ; (, 1 >. c At) ; (, beta) ; Ct ; , 1 >. n
)

zsyrkcore=: (4 : 0) ([ assert@(basiccs3 , basiccr3))
  'uplo trans'=. x
  'alpha At beta Ct'=. y
  n=. # Ct
  k=. (-. 'nN' e.~ {. trans) { $ At
  9 {:: zsyrkcd (, uplo) ; (, trans) ; (, n) ; (, k) ; (, alpha) ; At ; (, 1 >. c At) ; (, beta) ; Ct ; , 1 >. n
)

NB. ---------------------------------------------------------
NB. zherkcore
NB.
NB. Description:
NB.   Performs hermitian rank k operations:
NB.     C := alpha * A * A^H + beta * C  (1)
NB.   or
NB.     C := alpha * A^H * A + beta * C  (2)
NB.   with transposed matrices, where C is Hermitian
NB.
NB. Syntax:
NB.   Cupdt=. (uplo ; trans) zherkcore alpha ; At ; beta ; Ct
NB. where
NB.   uplo  - literal, case-insensitive, in which the head
NB.           specifies which triangular part of C is to be
NB.           referenced:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   trans - literal, case-insensitive, in which the head
NB.           specifies the operation:
NB.             'N'  NB. (1)
NB.             'C'  NB. (2)
NB.   alpha - scalar, real
NB.   At    - ka×na-matrix, A^T
NB.   beta  - scalar, real
NB.   Ct    - n×n-matrix with real diagonal, C^T
NB.   Cupdt - Ct with either LT (if uplo='U') or UT (if
NB.           uplo='L') updated
NB.   n     ≥ 0, the size of Ct and Cupdt and the number of
NB.           columns or rows in At
NB.   k     ≥ 0, the number of rows or columns in At
NB.   ka    = k if trans='N' or ka = n otherwise
NB.   na    = n if trans='N' or na = k otherwise

zherkcore=: (4 : 0) ([ assert@(basiccs3 , basiccr3))
  'uplo trans'=. x
  'alpha At beta Ct'=. y
  n=. # Ct
  k=. (-. 'nN' e.~ {. trans) { $ At
  9 {:: zherkcd (, uplo) ; (, trans) ; (, n) ; (, k) ; (, alpha) ; At ; (, 1 >. c At) ; (, beta) ; Ct ; , 1 >. n
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad      R/W in C    Operation
NB. dsyrkln    LT          C := alpha * A   * A^T + beta * C
NB. dsyrklt    LT          C := alpha * A^T * A   + beta * C
NB. dsyrkun    UT          C := alpha * A   * A^T + beta * C
NB. dsyrkut    UT          C := alpha * A^T * A   + beta * C
NB. zsyrkln    LT          C := alpha * A   * A^T + beta * C
NB. zsyrklt    LT          C := alpha * A^T * A   + beta * C
NB. zsyrkun    UT          C := alpha * A   * A^T + beta * C
NB. zsyrkut    UT          C := alpha * A^T * A   + beta * C
NB.
NB. Description:
NB.   Performs symmetric rank k operations:
NB.     C := alpha * A * A^T + beta * C
NB.   or
NB.     C := alpha * A^T * A + beta * C
NB.   where C is symmetric
NB.
NB. Syntax:
NB.   Cupd=. xsyrkxx alpha ; A ; beta ; C
NB. where
NB.   alpha - scalar
NB.   A     - na×ka-matrix
NB.   beta  - scalar
NB.   C     - n×n-matrix
NB.   Cupd  - C with either LT (for xsyrklx) or UT (for
NB.           xsyrkux) updated
NB.   n     ≥ 0, the size of C and Cupd and the number of
NB.           rows or columns in A
NB.   k     ≥ 0, the number of columns or rows in A
NB.   ka    = k for xsyrkxn or ka = n otherwise
NB.   na    = n for xsyrkxn or na = k otherwise
NB.
NB. Notes:
NB. - monad      provides BLAS'
NB.   dsyrkln    DSYRK('L','N',...)
NB.   dsyrklt    DSYRK('L','T',...), DSYRK('L','C',...)
NB.   dsyrkun    DSYRK('U','N',...)
NB.   dsyrkut    DSYRK('U','T',...), DSYRK('U','C',...)
NB.   zsyrkln    ZSYRK('L','N',...)
NB.   zsyrklt    ZSYRK('L','T',...)
NB.   zsyrkun    ZSYRK('U','N',...)
NB.   zsyrkut    ZSYRK('U','T',...)

dsyrkln=: 'ut'&dsyrkcore
dsyrklt=: 'un'&dsyrkcore
dsyrkun=: 'lt'&dsyrkcore
dsyrkut=: 'ln'&dsyrkcore

zsyrkln=: 'ut'&zsyrkcore
zsyrklt=: 'un'&zsyrkcore
zsyrkun=: 'lt'&zsyrkcore
zsyrkut=: 'ln'&zsyrkcore

NB. ---------------------------------------------------------
NB. Monad      R/W in C    Operation
NB. zherkln    LT          C := alpha * A   * A^H + beta * C
NB. zherklt    LT          C := alpha * A^H * A   + beta * C
NB. zherkun    UT          C := alpha * A   * A^H + beta * C
NB. zherkut    UT          C := alpha * A^H * A   + beta * C
NB.
NB. Description:
NB.   Performs hermitian rank k operations:
NB.     C := alpha * A * A^H + beta * C
NB.   or
NB.     C := alpha * A^H * A + beta * C
NB.   where C is Hermitian
NB.
NB. Syntax:
NB.   Cupd=. zherkxx alpha ; A ; beta ; C
NB. where
NB.   alpha - scalar, real
NB.   A     - na×ka-matrix
NB.   beta  - scalar, real
NB.   C     - n×n-matrix with real diagonal
NB.   Cupd  - C with either LT (for zherklx) or UT (for
NB.           zherkux) updated
NB.   n     ≥ 0, the size of C and Cupd and the number of
NB.           rows or columns in A
NB.   k     ≥ 0, the number of columns or rows in A
NB.   ka    = k for zherkxn or ka = n otherwise
NB.   na    = n for zherkxn or na = k otherwise
NB.
NB. Notes:
NB. - monad      provides BLAS'
NB.   zherkln    ZHERK('L','N',...)
NB.   zherklc    ZHERK('L','C',...)
NB.   zherkun    ZHERK('U','N',...)
NB.   zherkuc    ZHERK('U','C',...)

zherkln=: 'uc'&zherkcore
zherklc=: 'un'&zherkcore
zherkun=: 'lc'&zherkcore
zherkuc=: 'ln'&zherkcore
