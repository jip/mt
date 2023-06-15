NB. Rank 2k operations
NB.
NB. xsyr2kxx  Rank 2k operation with symmetric matrix
NB. zher2kxx  Rank 2k operation with Hermitian matrix
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
NB. dsyr2kcore
NB. zsyr2kcore
NB.
NB. Description:
NB.   Performs symmetric rank 2k operations:
NB.     C := alpha * A * B^T + alpha * B * A^T + beta * C  (1)
NB.   or
NB.     C := alpha * A^T * B + alpha * B^T * A + beta * C  (2)
NB.   with transposed matrices, where C is symmetric
NB.
NB. Syntax:
NB.   Cupdt=. (uplo ; trans) xsyr2kcore alpha ; At ; Bt ; beta ; Ct
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
NB.             'C'  NB. (2) for dsyr2kcore
NB.   alpha - scalar
NB.   At    - kab×nab-matrix, A^T
NB.   Bt    - kab×nab-matrix, B^T
NB.   beta  - scalar
NB.   Ct    - n×n-matrix, C^T
NB.   Cupdt - Ct with either LT (if uplo='U') or UT (if
NB.           uplo='L') updated
NB.   n     ≥ 0, the size of Ct and Cupdt and the number of
NB.           columns or rows in At and Bt
NB.   k     ≥ 0, the number of rows or columns in At and Bt
NB.   kab   = k if trans='N' or kab = n otherwise
NB.   nab   = n if trans='N' or nab = k otherwise
NB.
NB. Notes:
NB. - operate on transposed matrices to avoid transposition

dsyr2kcore=: (4 : 0) basicswp@([ assert@(basiccmp , basiccs4 , basiccr4))
  'uplo trans'=. x
  'alpha At Bt beta Ct'=. y
  n=. # Ct
  k=. (-. 'nN' e.~ {. trans) { $ At
  ld=. , 1 >. c At
  11 {:: dsyr2kcd (, uplo) ; (, trans) ; (, n) ; (, k) ; (, alpha) ; At ; ld ; Bt ; ld ; (, beta) ; Ct ; , 1 >. n
)

zsyr2kcore=: (4 : 0) basicswp@([ assert@(basiccmp , basiccs4 , basiccr4))
  'uplo trans'=. x
  'alpha At Bt beta Ct'=. y
  n=. # Ct
  k=. (-. 'nN' e.~ {. trans) { $ At
  ld=. , 1 >. c At
  11 {:: zsyr2kcd (, uplo) ; (, trans) ; (, n) ; (, k) ; (, alpha) ; At ; ld ; Bt ; ld ; (, beta) ; Ct ; , 1 >. n
)

NB. ---------------------------------------------------------
NB. zher2kcore
NB.
NB. Description:
NB.   Performs hermitian rank 2k operations:
NB.     C := alpha * A * B^H + conj(alpha) * B * A^H + beta * C  (1)
NB.   or
NB.     C := alpha * A^H * B + conj(alpha) * B^H * A + beta * C  (2)
NB.   with transposed matrices, where C is Hermitian
NB.
NB. Syntax:
NB.   Cupdt=. (uplo ; trans) zher2kcore alpha ; At ; Bt ; beta ; Ct
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
NB.   alpha - scalar
NB.   At    - kab×nab-matrix, A^T
NB.   Bt    - kab×nab-matrix, B^T
NB.   beta  - scalar, real
NB.   Ct    - n×n-matrix with real diagonal, C^T
NB.   Cupdt - Ct with either LT (if uplo='U') or UT (if
NB.           uplo='L') updated
NB.   n     ≥ 0, the size of Ct and Cupdt and the number of
NB.           columns or rows in At and Bt
NB.   k     ≥ 0, the number of rows or columns in At and Bt
NB.   kab   = k if trans='N' or kab = n otherwise
NB.   nab   = n if trans='N' or nab = k otherwise
NB.
NB. Notes:
NB. - operates on transposed matrices to avoid transposition

zher2kcore=: (4 : 0) basicswp@([ assert@(basiccmp , basiccs4 , basiccr4))
  'uplo trans'=. x
  'alpha At Bt beta Ct'=. y
  n=. # Ct
  k=. (-. 'nN' e.~ {. trans) { $ At
  ld=. , 1 >. c At
  11 {:: zher2kcd (, uplo) ; (, trans) ; (, n) ; (, k) ; (, alpha) ; At ; ld ; Bt ; ld ; (, beta) ; Ct ; , 1 >. n
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad       R/W in C    Operation
NB. dsyr2kln    LT          C := alpha * A   * B^T + alpha * B   * A^T + beta * C
NB. dsyr2klt    LT          C := alpha * A^T * B   + alpha * B^T * A   + beta * C
NB. dsyr2kun    UT          C := alpha * A   * B^T + alpha * B   * A^T + beta * C
NB. dsyr2kut    UT          C := alpha * A^T * B   + alpha * B^T * A   + beta * C
NB. zsyr2kln    LT          C := alpha * A   * B^T + alpha * B   * A^T + beta * C
NB. zsyr2klt    LT          C := alpha * A^T * B   + alpha * B^T * A   + beta * C
NB. zsyr2kun    UT          C := alpha * A   * B^T + alpha * B   * A^T + beta * C
NB. zsyr2kut    UT          C := alpha * A^T * B   + alpha * B^T * A   + beta * C
NB.
NB. Description:
NB.   Performs symmetric rank 2k operations:
NB.     C := alpha * A * B^T + alpha * B * A^T + beta * C
NB.   or
NB.     C := alpha * A^T * B + alpha * B^T * A + beta * C
NB.   where C is symmetric
NB.
NB. Syntax:
NB.   Cupd=. xsyr2kxx alpha ; A ; B ; beta ; C
NB. where
NB.   alpha - scalar
NB.   A     - nab×kab-matrix
NB.   B     - nab×kab-matrix
NB.   beta  - scalar
NB.   C     - n×n-matrix
NB.   Cupd  - C with either LT (for xsyr2klx) or UT (for
NB.           xsyr2kux) updated
NB.   n     ≥ 0, the size of C and Cupd and the number of
NB.           rows or columns in A and B
NB.   k     ≥ 0, the number of columns or rows in A and B
NB.   kab   = k for xsyr2kxn or kab = n otherwise
NB.   nab   = n for xsyr2kxn or nab = k otherwise
NB.
NB. Notes:
NB. - monad       provides BLAS'
NB.   dsyr2kln    DSYR2K('L','N',...)
NB.   dsyr2klt    DSYR2K('L','T',...), DSYR2K('L','C',...)
NB.   dsyr2kun    DSYR2K('U','N',...)
NB.   dsyr2kut    DSYR2K('U','T',...), DSYR2K('U','C',...)
NB.   zsyr2kln    ZSYR2K('L','N',...)
NB.   zsyr2klt    ZSYR2K('L','T',...)
NB.   zsyr2kun    ZSYR2K('U','N',...)
NB.   zsyr2kut    ZSYR2K('U','T',...)

dsyr2kln=: 'ut'&dsyr2kcore
dsyr2klt=: 'un'&dsyr2kcore
dsyr2kun=: 'lt'&dsyr2kcore
dsyr2kut=: 'ln'&dsyr2kcore

zsyr2kln=: 'ut'&zsyr2kcore
zsyr2klt=: 'un'&zsyr2kcore
zsyr2kun=: 'lt'&zsyr2kcore
zsyr2kut=: 'ln'&zsyr2kcore

NB. ---------------------------------------------------------
NB. Monad       R/W in C    Operation
NB. zher2kln    LT          C := alpha * A   * B^H + alpha * B   * A^H + beta * C
NB. zher2klc    LT          C := alpha * A^H * B   + alpha * B^H * A   + beta * C
NB. zher2kun    UT          C := alpha * A   * B^H + alpha * B   * A^H + beta * C
NB. zher2kuc    UT          C := alpha * A^H * B   + alpha * B^H * A   + beta * C
NB.
NB. Description:
NB.   Performs hermitian rank 2k operations:
NB.     C := alpha * A * B^H + conj(alpha) * B * A^H + beta * C
NB.   or
NB.     C := alpha * A^H * B + conj(alpha) * B^H * A + beta * C
NB.   where C is Hermitian
NB.
NB. Syntax:
NB.   Cupd=. zher2kxx alpha ; A ; B ; beta ; C
NB. where
NB.   alpha - scalar
NB.   A     - nab×kab-matrix
NB.   B     - nab×kab-matrix
NB.   beta  - scalar, real
NB.   C     - n×n-matrix with real diagonal
NB.   Cupd  - C with either LT (for zher2klx) or UT (for
NB.           zher2kux) updated
NB.   n     ≥ 0, the size of C and Cupd and the number of
NB.           rows or columns in A and B
NB.   k     ≥ 0, the number of columns or rows in A and B
NB.   kab   = k for zher2kxn or kab = n otherwise
NB.   nab   = n for zher2kxn or nab = k otherwise
NB.
NB. Notes:
NB. - monad       provides BLAS'
NB.   zher2kln    ZHER2K('L','N',...)
NB.   zher2klc    ZHER2K('L','C',...)
NB.   zher2kun    ZHER2K('U','N',...)
NB.   zher2kuc    ZHER2K('U','C',...)

zher2kln=: 'uc'&zher2kcore
zher2klc=: 'un'&zher2kcore
zher2kun=: 'lc'&zher2kcore
zher2kuc=: 'ln'&zher2kcore
