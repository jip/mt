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

coclass 'mtbla'

NB. =========================================================
NB. Includes

require 'math/mt/test/blas/blas'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Dyad         Domain
NB. dsyrkcore    real
NB. zsyrkcore    complex
NB.
NB. Description:
NB.   Performs the symmetric rank k operation:
NB.     C := alpha * A * A^T + beta * C  (1)
NB.   or
NB.     C := alpha * A^T * A + beta * C  (2)
NB.   with transposed matrices, where C is symmetric
NB.
NB. Syntax:
NB.   CCupdt=. (uplo ; trans) xsyrkcore alpha ; At ; beta ; CCt
NB. where
NB.   uplo   - literal, case-insensitive, in which the head
NB.            specifies which triangular part of C is to be
NB.            referenced:
NB.              'L'  NB. LT
NB.              'U'  NB. UT
NB.   trans  - literal, case-insensitive, in which the head
NB.            specifies the operation:
NB.              'N'  NB. (1)
NB.              'T'  NB. (2)
NB.              'C'  NB. (2) for dsyrkcore
NB.   alpha  - scalar
NB.   At     - ka×na-matrix, A^T
NB.   beta   - scalar
NB.   CCt    - n×n-matrix, contains either LT or UT or both
NB.            part(s) of C^T
NB.   CCupdt - CCt with either LT (if uplo='U') or UT (if
NB.            uplo='L') updated
NB.   C      - n×n-matrix, symmetric
NB.   n      ≥ 0, the size of C, CCt and CCupdt and the
NB.            number of columns or rows in At
NB.   k      ≥ 0, the number of rows or columns in At
NB.   ka     = k if trans='N' or ka = n otherwise
NB.   na     = n if trans='N' or na = k otherwise
NB.
NB. Notes:
NB. - operate on transposed matrices to avoid transposition

dsyrkcore=: (4 : 0) ([ assert@(basiccs3 , basiccr3))
  'uplo trans'=. x
  'alpha At beta CCt'=. y
  n=. # CCt
  k=. (-. 'nN' e.~ {. trans) { $ At
  9 {:: dsyrkcd (, uplo) ; (, trans) ; (, n) ; (, k) ; (, alpha) ; At ; (, 1 >. c At) ; (, beta) ; CCt ; , 1 >. n
)

zsyrkcore=: (4 : 0) ([ assert@(basiccs3 , basiccr3))
  'uplo trans'=. x
  'alpha At beta CCt'=. y
  n=. # CCt
  k=. (-. 'nN' e.~ {. trans) { $ At
  9 {:: zsyrkcd (, uplo) ; (, trans) ; (, n) ; (, k) ; (, alpha) ; At ; (, 1 >. c At) ; (, beta) ; CCt ; , 1 >. n
)

NB. ---------------------------------------------------------
NB. zherkcore
NB.
NB. Description:
NB.   Performs the hermitian rank k operation:
NB.     C := alpha * A * A^H + beta * C  (1)
NB.   or
NB.     C := alpha * A^H * A + beta * C  (2)
NB.   with transposed matrices, where C is Hermitian
NB.
NB. Syntax:
NB.   CCupdt=. (uplo ; trans) zherkcore alpha ; At ; beta ; CCt
NB. where
NB.   uplo   - literal, case-insensitive, in which the head
NB.            specifies which triangular part of C is to be
NB.            referenced:
NB.              'L'  NB. LT
NB.              'U'  NB. UT
NB.   trans  - literal, case-insensitive, in which the head
NB.            specifies the operation:
NB.              'N'  NB. (1)
NB.              'C'  NB. (2)
NB.   alpha  - scalar, real
NB.   At     - ka×na-matrix, A^T
NB.   beta   - scalar, real
NB.   CCt    - n×n-matrix, contains either LT or UT or both
NB.            part(s) of C^T
NB.   CCupdt - CCt with either LT (if uplo='U') or UT (if
NB.            uplo='L') updated
NB.   C      - n×n-matrix, Hermitian
NB.   n      ≥ 0, the size of C, CCt and CCupdt and the
NB.            number of columns or rows in At
NB.   k      ≥ 0, the number of rows or columns in At
NB.   ka     = k if trans='N' or ka = n otherwise
NB.   na     = n if trans='N' or na = k otherwise

zherkcore=: (4 : 0) ([ assert@(basiccs3 , basiccr3))
  'uplo trans'=. x
  'alpha At beta CCt'=. y
  n=. # CCt
  k=. (-. 'nN' e.~ {. trans) { $ At
  9 {:: zherkcd (, uplo) ; (, trans) ; (, n) ; (, k) ; (, alpha) ; At ; (, 1 >. c At) ; (, beta) ; CCt ; , 1 >. n
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad      Domain     R/W in C    op1(A)    op2(A)
NB. dsyrkln    real       LT          A         A^T
NB. dsyrklt    real       LT          A^T       A
NB. dsyrkun    real       UT          A         A^T
NB. dsyrkut    real       UT          A^T       A
NB. zsyrkln    complex    LT          A         A^T
NB. zsyrklt    complex    LT          A^T       A
NB. zsyrkun    complex    UT          A         A^T
NB. zsyrkut    complex    UT          A^T       A
NB.
NB. Description:
NB.   Performs the symmetric rank k operation:
NB.     C := alpha * op1(A) * op2(A) + beta * C
NB.   where C is symmetric and opX(A) is either A or A^T
NB.
NB. Syntax:
NB.   CCupd=. xsyrkxx alpha ; A ; beta ; CC
NB. where
NB.   alpha - scalar
NB.   A     - na×ka-matrix
NB.   beta  - scalar
NB.   CC    - n×n-matrix, contains either LT or UT or both
NB.           part(s) of C
NB.   CCupd - CC with either LT (for xsyrklx) or UT (for
NB.           xsyrkux) updated
NB.   C     - n×n-matrix, symmetric
NB.   n     ≥ 0, the size of C, CC and CCupd and the number
NB.           of rows or columns in A
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
NB. Monad      Domain     R/W in C    op1(A)    op2(A)
NB. zherkln    complex    LT          A         A^H
NB. zherklt    complex    LT          A^H       A
NB. zherkun    complex    UT          A         A^H
NB. zherkut    complex    UT          A^H       A
NB.
NB. Description:
NB.   Performs the hermitian rank k operation:
NB.     C := alpha * op1(A) * op2(A) + beta * C
NB.   where C is Hermitian
NB.
NB. Syntax:
NB.   CCupd=. zherkxx alpha ; A ; beta ; CC
NB. where
NB.   alpha - scalar, real
NB.   A     - na×ka-matrix
NB.   beta  - scalar, real
NB.   CC    - n×n-matrix, contains either LT or UT or both
NB.           part(s) of C
NB.   CCupd - CC with either LT (for zherklx) or UT (for
NB.           zherkux) updated
NB.   C     - n×n-matrix, Hermitian
NB.   n     ≥ 0, the size of C, CC and CCupd and the number
NB.           of rows or columns in A
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
