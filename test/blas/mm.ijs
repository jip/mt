NB. Matrix-matrix operations
NB.
NB. xgemmxx    Matrix-matrix operation with general matrix
NB. xsymmxx    Matrix-matrix operation with symmetric matrix
NB. zhemmxx    Matrix-matrix operation with Hermitian matrix
NB. xtrmmxxxx  Matrix-matrix operation with triangular matrix
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
NB. dgemmcore
NB. zgemmcore
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C
NB.   with transposed matrices, where opX(M) is either M, M^T
NB.   or M^H
NB.
NB. Syntax:
NB.   Cupdt=. (transA ; transB) xgemmcore alpha ; At ; Bt ; beta ; Ct
NB. where
NB.   transA - literal, case-insensitive, in which the head
NB.            specifies the form of op1(A):
NB.              'N'  NB. op1(A) := A    (no transpose)
NB.              'T'  NB. op1(A) := A^T  (transpose)
NB.              'C'  NB. op1(A) := A^T  (transpose)           for dgemmcore
NB.                   NB.        := A^H  (conjugate transpose) for zgemmcore
NB.   transB - literal, case-insensitive, in which the head
NB.            specifies the form of op2(B):
NB.              'N'  NB. op2(B) := B    (no transpose)
NB.              'T'  NB. op2(B) := B^T  (transpose)
NB.              'C'  NB. op2(B) := B^T  (transpose)           for dgemmcore
NB.                   NB.        := B^H  (conjugate transpose) for zgemmcore
NB.   alpha  - scalar
NB.   At     - ka×ma-matrix, A^T
NB.   Bt     - nb×kb-matrix, B^T
NB.   beta   - scalar
NB.   Ct     - n×m-matrix, C^T
NB.   Cupdt  - an updated Ct
NB.   m      ≥ 0, the number of rows in C and op1(A)
NB.   n      ≥ 0, the number of columns in C and op2(B)
NB.   k      ≥ 0, the number of columns in op1(A) and the
NB.            number of rows in op2(B)
NB.   ma     = m if transA='N' or ma = k otherwise
NB.   ka     = k if transA='N' or ka = m otherwise
NB.   kb     = k if transB='N' or kb = n otherwise
NB.   nb     = n if transB='N' or nb = k otherwise
NB.
NB. Notes:
NB. - operate on transposed matrices to avoid transposition

dgemmcore=: (4 : 0) basicswp@([ assert@basiccr4)
  'transA transB'=. x
  'alpha At Bt beta Ct'=. y
  'n m'=. $ Ct
  k=. (-. 'nN' e.~ {. transA) { $ At
  12 {:: dgemmcd (, transA) ; (, transB) ; (, m) ; (, n) ; (, k) ; (, alpha) ; At ; (, 1 >. c At) ; Bt ; (, 1 >. c Bt) ; (, beta) ; Ct ; , 1 >. m
)

zgemmcore=: (4 : 0) basicswp@([ assert@basiccr4)
  'transA transB'=. x
  'alpha At Bt beta Ct'=. y
  'n m'=. $ Ct
  k=. (-. 'nN' e.~ {. transA) { $ At
  12 {:: zgemmcd (, transA) ; (, transB) ; (, m) ; (, n) ; (, k) ; (, alpha) ; At ; (, 1 >. c At) ; Bt ; (, 1 >. c Bt) ; (, beta) ; Ct ; , 1 >. m
)

NB. ---------------------------------------------------------
NB. Dyad         A            diag(A)
NB. dsymmcore    symmetric    real
NB. zsymmcore    symmetric    any
NB. zhemmcore    Hermitian    real
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * A * B + beta * C  (1)
NB.   or
NB.     C := alpha * B * A + beta * C  (2)
NB.   with transposed matrices, where A is Hermitian
NB.   (symmetric)
NB.
NB. Syntax:
NB.   Cupdt=. (side ; uplo) xxxmmcore alpha ; AAt ; Bt ; beta ; Ct
NB. where
NB.   side  - literal, case-insensitive, in which the head
NB.           specifies the side of A:
NB.             'L'  NB. to perform (1) (A on the left)
NB.             'R'  NB. to perform (2) (A on the right)
NB.   uplo  - literal, case-insensitive, in which the head
NB.           specifies which triangular part of A is to be
NB.           referenced:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   alpha - scalar
NB.   AAt   - ma×ma-matrix, contains either lower or upper or
NB.           both part(s) of A^T
NB.   Bt    - n×m-matrix, B^T
NB.   beta  - scalar
NB.   Ct    - n×m-matrix, C^T
NB.   Cupdt - an updated Ct
NB.   A     - ma×ma-matrix, Hermitian (symmetric)
NB.   m     ≥ 0, the number of rows in C and B
NB.   n     ≥ 0, the number of columns in C and B
NB.   ma    = m if side='L' or ma = n otherwise
NB.
NB. Notes:
NB. - operate on transposed matrices to avoid transposition

dsymmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  'side uplo'=. x
  'alpha AAt Bt beta Ct'=. y
  'n m'=. $ Ct
  ld=. , 1 >. m
  11 {:: dsymmcd (, side) ; (, uplo) ; (, m) ; (, n) ; (, alpha) ; AAt ; (, 1 >. c AAt) ; Bt ; ld ; (, beta) ; Ct ; ld
)

zsymmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  'side uplo'=. x
  'alpha AAt Bt beta Ct'=. y
  'n m'=. $ Ct
  ld=. , 1 >. m
  11 {:: zsymmcd (, side) ; (, uplo) ; (, m) ; (, n) ; (, alpha) ; AAt ; (, 1 >. c AAt) ; Bt ; ld ; (, beta) ; Ct ; ld
)

zhemmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  'side uplo'=. x
  'alpha AAt Bt beta Ct'=. y
  'n m'=. $ Ct
  ld=. , 1 >. m
  11 {:: zhemmcd (, side) ; (, uplo) ; (, m) ; (, n) ; (, alpha) ; AAt ; (, 1 >. c AAt) ; Bt ; ld ; (, beta) ; Ct ; ld
)

NB. ---------------------------------------------------------
NB. dtrmmcore
NB. ztrmmcore
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     B := alpha * op(A) * B  (1)
NB.   or
NB.     B := alpha * B * op(A)  (2)
NB.   with transposed matrices, where A is triangular, and
NB.   op(A) is either A, A^T or A^H
NB.
NB. Syntax:
NB.   Bupdt=. (side ; uplo ; trans ; diag) xtrmmcore alpha ; AAt ; Bt
NB. where
NB.   side  - literal, case-insensitive, in which the head
NB.           specifies the side of op(A):
NB.             'L'  NB. to perform (1) (op(A) on the left)
NB.             'R'  NB. to perform (2) (op(A) on the right)
NB.   uplo  - literal, case-insensitive, in which the head
NB.           specifies whether the matrix A is upper or
NB.           lower triangular:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   trans - literal, case-insensitive, in which the head
NB.           specifies the form of op(A):
NB.             'N'  NB. op(A) := A    (no transpose)
NB.             'T'  NB. op(A) := A^T  (transpose)
NB.             'C'  NB. op(A) := A^T  (transpose)           for dtrmmcore
NB.                  NB. op(A) := A^H  (conjugate transpose) for ztrmmcore
NB.   diag  - literal, case-insensitive, in which the head
NB.           specifies the form of A:
NB.             'N'  NB. A is either L or U
NB.             'U'  NB. A is either L1 or U1, diagonal
NB.                  NB.   elements of A are not referenced
NB.   alpha - scalar
NB.   AAt   - k×k-matrix, contains either non-zero or both
NB.           triangular parts of A^T
NB.   Bt    - n×m-matrix, B^T
NB.   Bupdt - an updated Bt
NB.   A     - k×k-matrix, triangular
NB.   m     ≥ 0, the number of rows in B
NB.   n     ≥ 0, the number of columns in B
NB.   k     = m if side='L' or k = n otherwise
NB.
NB. Notes:
NB. - operate on transposed matrices to avoid transposition

dtrmmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  'side uplo trans diag'=. x
  'alpha AAt Bt'=. y
  'n m'=. $ Bt
  10 {:: dtrmmcd (, side) ; (, uplo) ; (, trans) ; (, diag) ; (, m) ; (, n) ; (, alpha) ; AAt ; (, 1 >. c AAt) ; Bt ; , 1 >. m
)

ztrmmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  'side uplo trans diag'=. x
  'alpha AAt Bt'=. y
  'n m'=. $ Bt
  10 {:: ztrmmcd (, side) ; (, uplo) ; (, trans) ; (, diag) ; (, m) ; (, n) ; (, alpha) ; AAt ; (, 1 >. c AAt) ; Bt ; , 1 >. m
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad      Operation
NB. dgemmnn    C := alpha * A   * B   + beta * C
NB. dgemmnt    C := alpha * A   * B^T + beta * C
NB. dgemmtn    C := alpha * A^T * B   + beta * C
NB. dgemmtt    C := alpha * A^T * B^T + beta * C
NB. zgemmnn    C := alpha * A   * B   + beta * C
NB. zgemmnt    C := alpha * A   * B^T + beta * C
NB. zgemmnc    C := alpha * A   * B^H + beta * C
NB. zgemmtn    C := alpha * A^T * B   + beta * C
NB. zgemmtt    C := alpha * A^T * B^T + beta * C
NB. zgemmtc    C := alpha * A^T * B^H + beta * C
NB. zgemmcn    C := alpha * A^H * B   + beta * C
NB. zgemmct    C := alpha * A^H * B^T + beta * C
NB. zgemmcc    C := alpha * A^H * B^H + beta * C
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C
NB.   where opX(M) is either M, M^T or M^H
NB.
NB. Syntax:
NB.   Cupd=. xgemmxx alpha ; A ; B ; beta ; C
NB. where
NB.   alpha - scalar
NB.   A     - ma×ka-matrix
NB.   B     - kb×nb-matrix
NB.   beta  - scalar
NB.   C     - m×n-matrix
NB.   Cupd  - an updated C
NB.   m     ≥ 0, the number of rows in C, Cupd and op1(A)
NB.   n     ≥ 0, the number of columns in C, Cupd and op2(B)
NB.   k     ≥ 0, the number of columns in op1(A) and the
NB.           number of rows in op2(B)
NB.   ma    = m for xgemmnx or ma = k otherwise
NB.   ka    = k for xgemmnx or ka = m otherwise
NB.   kb    = k for xgemmxn or kb = n otherwise
NB.   nb    = n for xgemmxn or nb = k otherwise
NB.
NB. Notes:
NB. - monad      provides BLAS'
NB.   dgemmnn    DGEMM('N','N',...)
NB.   dgemmnt    DGEMM('N','T',...), DGEMM('N','C',...)
NB.   dgemmtn    DGEMM('T','N',...), DGEMM('C','N',...)
NB.   dgemmtt    DGEMM('T','T',...), DGEMM('C','T',...), DGEMM('T','C',...), DGEMM('C','C',...)
NB.   zgemmnn    ZGEMM('N','N',...)
NB.   zgemmnt    ZGEMM('N','T',...)
NB.   zgemmnc    ZGEMM('N','C',...)
NB.   zgemmtn    ZGEMM('T','N',...)
NB.   zgemmtt    ZGEMM('T','T',...)
NB.   zgemmtc    ZGEMM('T','C',...)
NB.   zgemmcn    ZGEMM('C','N',...)
NB.   zgemmct    ZGEMM('C','T',...)
NB.   zgemmcc    ZGEMM('C','C',...)

dgemmnn=: 'nn'&dgemmcore
dgemmnt=: 'tn'&dgemmcore
dgemmtn=: 'nt'&dgemmcore
dgemmtt=: 'tt'&dgemmcore

zgemmnn=: 'nn'&zgemmcore
zgemmnt=: 'tn'&zgemmcore
zgemmnc=: 'cn'&zgemmcore
zgemmtn=: 'nt'&zgemmcore
zgemmtt=: 'tt'&zgemmcore
zgemmtc=: 'ct'&zgemmcore
zgemmcn=: 'nc'&zgemmcore
zgemmct=: 'tc'&zgemmcore
zgemmcc=: 'cc'&zgemmcore

NB. ---------------------------------------------------------
NB. Monad      A            diag(A)   Reads in A    Operation
NB. dsymmll    symmetric    real      LT            C := alpha * A * B + beta * C
NB. dsymmlu    symmetric    real      UT            C := alpha * A * B + beta * C
NB. dsymmrl    symmetric    real      LT            C := alpha * B * A + beta * C
NB. dsymmru    symmetric    real      UT            C := alpha * B * A + beta * C
NB. zsymmll    symmetric    any       LT            C := alpha * A * B + beta * C
NB. zsymmlu    symmetric    any       UT            C := alpha * A * B + beta * C
NB. zsymmrl    symmetric    any       LT            C := alpha * B * A + beta * C
NB. zsymmru    symmetric    any       UT            C := alpha * B * A + beta * C
NB. zhemmll    Hermitian    real      LT            C := alpha * A * B + beta * C
NB. zhemmlu    Hermitian    real      UT            C := alpha * A * B + beta * C
NB. zhemmrl    Hermitian    real      LT            C := alpha * B * A + beta * C
NB. zhemmru    Hermitian    real      UT            C := alpha * B * A + beta * C
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * A * B + beta * C
NB.   or
NB.     C := alpha * B * A + beta * C
NB.   where A is Hermitian (symmetric)
NB.
NB. Syntax:
NB.   Cupd=. xxxmmxx alpha ; AA ; B ; beta ; C
NB. where
NB.   alpha - scalar
NB.   AA    - ma×ma-matrix, contains either lower or upper or
NB.           both part(s) of A
NB.   A     - ma×ma-matrix, Hermitian (symmetric)
NB.   B     - m×n-matrix
NB.   beta  - scalar
NB.   C     - m×n-matrix
NB.   Cupd  - an updated C
NB.   m     ≥ 0, the number of rows in C, Cupd and B
NB.   n     ≥ 0, the number of columns in C, Cupd and B
NB.   ma    = m for xxxmmlx or ma = n otherwise
NB.
NB. Notes:
NB. - monad      provides BLAS'
NB.   dsymmll    DSYMM('L','L',...)
NB.   dsymmlu    DSYMM('L','U',...)
NB.   dsymmrl    DSYMM('R','L',...)
NB.   dsymmru    DSYMM('R','U',...)
NB.   zsymmll    ZSYMM('L','L',...)
NB.   zsymmlu    ZSYMM('L','U',...)
NB.   zsymmrl    ZSYMM('R','L',...)
NB.   zsymmru    ZSYMM('R','U',...)
NB.   zhemmll    ZHEMM('L','L',...)
NB.   zhemmlu    ZHEMM('L','U',...)
NB.   zhemmrl    ZHEMM('R','L',...)
NB.   zhemmru    ZHEMM('R','U',...)

dsymmll=: 'ru'&dsymmcore
dsymmlu=: 'rl'&dsymmcore
dsymmrl=: 'lu'&dsymmcore
dsymmru=: 'll'&dsymmcore

zsymmll=: 'ru'&zsymmcore
zsymmlu=: 'rl'&zsymmcore
zsymmrl=: 'lu'&zsymmcore
zsymmru=: 'll'&zsymmcore

zhemmll=: 'ru'&zhemmcore
zhemmlu=: 'rl'&zhemmcore
zhemmrl=: 'lu'&zhemmcore
zhemmru=: 'll'&zhemmcore

NB. ---------------------------------------------------------
NB. Monad        Reads in A    Operation
NB. dtrmmllnn     LT           B := alpha * L    * B
NB. dtrmmllnu    SLT           B := alpha * L1   * B
NB. dtrmmlltn     LT           B := alpha * L ^T * B
NB. dtrmmlltu    SLT           B := alpha * L1^T * B
NB. dtrmmlunn     UT           B := alpha * U    * B
NB. dtrmmlunu    SUT           B := alpha * U1   * B
NB. dtrmmlutn     UT           B := alpha * U ^T * B
NB. dtrmmlutu    SUT           B := alpha * U1^T * B
NB. dtrmmrlnn     LT           B := alpha * B * L
NB. dtrmmrlnu    SLT           B := alpha * B * L1
NB. dtrmmrltn     LT           B := alpha * B * L ^T
NB. dtrmmrltu    SLT           B := alpha * B * L1^T
NB. dtrmmrunn     UT           B := alpha * B * U
NB. dtrmmrunu    SUT           B := alpha * B * U1
NB. dtrmmrutn     UT           B := alpha * B * U ^T
NB. dtrmmrutu    SUT           B := alpha * B * U1^T
NB. ztrmmllnn     LT           B := alpha * L    * B
NB. ztrmmllnu    SLT           B := alpha * L1   * B
NB. ztrmmlltn     LT           B := alpha * L ^T * B
NB. ztrmmlltu    SLT           B := alpha * L1^T * B
NB. ztrmmllcn     LT           B := alpha * L ^H * B
NB. ztrmmllcu    SLT           B := alpha * L1^H * B
NB. ztrmmlunn     UT           B := alpha * U    * B
NB. ztrmmlunu    SUT           B := alpha * U1   * B
NB. ztrmmlutn     UT           B := alpha * U ^T * B
NB. ztrmmlutu    SUT           B := alpha * U1^T * B
NB. ztrmmlucn     UT           B := alpha * U ^H * B
NB. ztrmmlucu    SUT           B := alpha * U1^H * B
NB. ztrmmrlnn     LT           B := alpha * B * L
NB. ztrmmrlnu    SLT           B := alpha * B * L1
NB. ztrmmrltn     LT           B := alpha * B * L ^T
NB. ztrmmrltu    SLT           B := alpha * B * L1^T
NB. ztrmmrlcn     LT           B := alpha * B * L ^H
NB. ztrmmrlcu    SLT           B := alpha * B * L1^H
NB. ztrmmrunn     UT           B := alpha * B * U
NB. ztrmmrunu    SUT           B := alpha * B * U1
NB. ztrmmrutn     UT           B := alpha * B * U ^T
NB. ztrmmrutu    SUT           B := alpha * B * U1^T
NB. ztrmmrucn     UT           B := alpha * B * U ^H
NB. ztrmmrucu    SUT           B := alpha * B * U1^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     B := alpha * op(A) * B
NB.   or
NB.     B := alpha * B * op(A)
NB.   where A is triangular, and op(A) is either A, A^T or
NB.   A^H
NB.
NB. Syntax:
NB.   Bupd=. xtrmmxxxx alpha ; AA ; B
NB. where
NB.   alpha - scalar
NB.   AA    - k×k-matrix, contains either non-zero or both
NB.           part(s) of A
NB.   A     - k×k-matrix, triangular
NB.   B     - m×n-matrix
NB.   Bupd  - an updated B
NB.   m     ≥ 0, the number of rows in B and Bupd
NB.   n     ≥ 0, the number of columns in B and Bupd
NB.   k     = m for xtrmmlxxx or k = n xtrmmrxxx
NB.
NB. Notes:
NB. - monad        provides BLAS'
NB.   dtrmmllnn    DTRMM('L','L','N','N',...)
NB.   dtrmmllnu    DTRMM('L','L','N','U',...)
NB.   dtrmmllnu    DTRMM('L','L','N','U',...)
NB.   dtrmmlltn    DTRMM('L','L','T','N',...), DTRMM('L','L','C','N',...)
NB.   dtrmmlltu    DTRMM('L','L','T','U',...), DTRMM('L','L','C','U',...)
NB.   dtrmmlunn    DTRMM('L','U','N','N',...)
NB.   dtrmmlunu    DTRMM('L','U','N','U',...)
NB.   dtrmmlutn    DTRMM('L','U','T','N',...), DTRMM('L','U','C','N',...)
NB.   dtrmmlutu    DTRMM('L','U','T','U',...), DTRMM('L','U','C','U',...)
NB.   dtrmmrlnn    DTRMM('R','L','N','N',...)
NB.   dtrmmrlnu    DTRMM('R','L','N','U',...)
NB.   dtrmmrltn    DTRMM('R','L','T','N',...), DTRMM('R','L','C','N',...)
NB.   dtrmmrltu    DTRMM('R','L','T','U',...), DTRMM('R','L','C','U',...)
NB.   dtrmmrunn    DTRMM('R','U','N','N',...)
NB.   dtrmmrunu    DTRMM('R','U','N','U',...)
NB.   dtrmmrutn    DTRMM('R','U','T','N',...), DTRMM('R','U','C','N',...)
NB.   dtrmmrutu    DTRMM('R','U','T','U',...), DTRMM('R','U','C','U',...)
NB.   ztrmmllnn    ZTRMM('L','L','N','N',...)
NB.   ztrmmllnu    ZTRMM('L','L','N','U',...)
NB.   ztrmmlltn    ZTRMM('L','L','T','N',...)
NB.   ztrmmlltu    ZTRMM('L','L','T','U',...)
NB.   ztrmmllcn    ZTRMM('L','L','C','N',...)
NB.   ztrmmllcu    ZTRMM('L','L','C','U',...)
NB.   ztrmmlunn    ZTRMM('L','U','N','N',...)
NB.   ztrmmlunu    ZTRMM('L','U','N','U',...)
NB.   ztrmmlutn    ZTRMM('L','U','T','N',...)
NB.   ztrmmlutu    ZTRMM('L','U','T','U',...)
NB.   ztrmmlucn    ZTRMM('L','U','C','N',...)
NB.   ztrmmlucu    ZTRMM('L','U','C','U',...)
NB.   ztrmmrlnn    ZTRMM('R','L','N','N',...)
NB.   ztrmmrlnu    ZTRMM('R','L','N','U',...)
NB.   ztrmmrltn    ZTRMM('R','L','T','N',...)
NB.   ztrmmrltu    ZTRMM('R','L','T','U',...)
NB.   ztrmmrlcn    ZTRMM('R','L','C','N',...)
NB.   ztrmmrlcu    ZTRMM('R','L','C','U',...)
NB.   ztrmmrunn    ZTRMM('R','U','N','N',...)
NB.   ztrmmrunu    ZTRMM('R','U','N','U',...)
NB.   ztrmmrutn    ZTRMM('R','U','T','N',...)
NB.   ztrmmrutu    ZTRMM('R','U','T','U',...)
NB.   ztrmmrucn    ZTRMM('R','U','C','N',...)
NB.   ztrmmrucu    ZTRMM('R','U','C','U',...)

dtrmmllnn=: 'runn'&dtrmmcore
dtrmmllnu=: 'runu'&dtrmmcore
dtrmmlltn=: 'rutn'&dtrmmcore
dtrmmlltu=: 'rutu'&dtrmmcore
dtrmmlunn=: 'rlnn'&dtrmmcore
dtrmmlunu=: 'rlnu'&dtrmmcore
dtrmmlutn=: 'rltn'&dtrmmcore
dtrmmlutu=: 'rltu'&dtrmmcore
dtrmmrlnn=: 'lunn'&dtrmmcore
dtrmmrlnu=: 'lunu'&dtrmmcore
dtrmmrltn=: 'lutn'&dtrmmcore
dtrmmrltu=: 'lutu'&dtrmmcore
dtrmmrunn=: 'llnn'&dtrmmcore
dtrmmrunu=: 'llnu'&dtrmmcore
dtrmmrutn=: 'lltn'&dtrmmcore
dtrmmrutu=: 'lltu'&dtrmmcore

ztrmmllnn=: 'runn'&ztrmmcore
ztrmmllnu=: 'runu'&ztrmmcore
ztrmmlltn=: 'rutn'&ztrmmcore
ztrmmlltu=: 'rutu'&ztrmmcore
ztrmmllcn=: 'rucn'&ztrmmcore
ztrmmllcu=: 'rucu'&ztrmmcore
ztrmmlunn=: 'rlnn'&ztrmmcore
ztrmmlunu=: 'rlnu'&ztrmmcore
ztrmmlutn=: 'rltn'&ztrmmcore
ztrmmlutu=: 'rltu'&ztrmmcore
ztrmmlucn=: 'rlcn'&ztrmmcore
ztrmmlucu=: 'rlcu'&ztrmmcore
ztrmmrlnn=: 'lunn'&ztrmmcore
ztrmmrlnu=: 'lunu'&ztrmmcore
ztrmmrltn=: 'lutn'&ztrmmcore
ztrmmrltu=: 'lutu'&ztrmmcore
ztrmmrlcn=: 'lucn'&ztrmmcore
ztrmmrlcu=: 'lucu'&ztrmmcore
ztrmmrunn=: 'llnn'&ztrmmcore
ztrmmrunu=: 'llnu'&ztrmmcore
ztrmmrutn=: 'lltn'&ztrmmcore
ztrmmrutu=: 'lltu'&ztrmmcore
ztrmmrucn=: 'llcn'&ztrmmcore
ztrmmrucu=: 'llcu'&ztrmmcore
