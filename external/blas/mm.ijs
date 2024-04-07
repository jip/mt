NB. Matrix-matrix operations
NB.
NB. xgemmxx    Matrix-matrix operation with general matrix
NB. xsymmxx    Matrix-matrix operation with symmetric matrix
NB. zhemmxx    Matrix-matrix operation with Hermitian matrix
NB. xtrmmxxxx  Matrix-matrix operation with triangular matrix
NB.
NB. Copyright 2010,2011,2013,2017,2018,2020,2021,2023,2024
NB.           Igor Zhuravlov
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
NB. dgemmcore    real
NB. zgemmcore    complex
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
NB.   transA - string, case-insensitive, in which the head
NB.            specifies the form of op1(A):
NB.              'N'  NB. op1(A) := A    (no transpose)
NB.              'T'  NB. op1(A) := A^T  (transpose)
NB.              'C'  NB. op1(A) := A^T  (transpose)           for dgemmcore
NB.                   NB.        := A^H  (conjugate transpose) for zgemmcore
NB.   transB - string, case-insensitive, in which the head
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
  12 {:: dgemm_cd (, transA) ; (, transB) ; (, m) ; (, n) ; (, k) ; (, alpha) ; At ; (, 1 >. c At) ; Bt ; (, 1 >. c Bt) ; (, beta) ; Ct ; , 1 >. m
)

zgemmcore=: (4 : 0) basicswp@([ assert@basiccr4)
  'transA transB'=. x
  'alpha At Bt beta Ct'=. y
  'n m'=. $ Ct
  k=. (-. 'nN' e.~ {. transA) { $ At
  12 {:: zgemm_cd (, transA) ; (, transB) ; (, m) ; (, n) ; (, k) ; (, alpha) ; At ; (, 1 >. c At) ; Bt ; (, 1 >. c Bt) ; (, beta) ; Ct ; , 1 >. m
)

NB. ---------------------------------------------------------
NB. Dyad         Domain     A
NB. dsymmcore    real       SY
NB. zsymmcore    complex    SY
NB. zhemmcore    complex    HE
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
NB.   side  - string, case-insensitive, in which the head
NB.           specifies the side of A:
NB.             'L'  NB. to perform (1) (A on the left)
NB.             'R'  NB. to perform (2) (A on the right)
NB.   uplo  - string, case-insensitive, in which the head
NB.           specifies which triangular part of A is to be
NB.           referenced:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   alpha - scalar
NB.   AAt   - mn×mn-matrix, contains either LT or UT or both
NB.           part(s) of A^T
NB.   Bt    - n×m-matrix, B^T
NB.   beta  - scalar
NB.   Ct    - n×m-matrix, C^T
NB.   Cupdt - an updated Ct
NB.   A     - mn×mn-matrix, Hermitian (symmetric)
NB.   m     ≥ 0, the number of rows in C and B
NB.   n     ≥ 0, the number of columns in C and B
NB.   mn    = m if side='L' or mn = n otherwise
NB.
NB. Notes:
NB. - operate on transposed matrices to avoid transposition

dsymmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  'side uplo'=. x
  'alpha AAt Bt beta Ct'=. y
  'n m'=. $ Ct
  ld=. , 1 >. m
  11 {:: dsymm_cd (, side) ; (, uplo) ; (, m) ; (, n) ; (, alpha) ; AAt ; (, 1 >. c AAt) ; Bt ; ld ; (, beta) ; Ct ; ld
)

zsymmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  'side uplo'=. x
  'alpha AAt Bt beta Ct'=. y
  'n m'=. $ Ct
  ld=. , 1 >. m
  11 {:: zsymm_cd (, side) ; (, uplo) ; (, m) ; (, n) ; (, alpha) ; AAt ; (, 1 >. c AAt) ; Bt ; ld ; (, beta) ; Ct ; ld
)

zhemmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  'side uplo'=. x
  'alpha AAt Bt beta Ct'=. y
  'n m'=. $ Ct
  ld=. , 1 >. m
  11 {:: zhemm_cd (, side) ; (, uplo) ; (, m) ; (, n) ; (, alpha) ; AAt ; (, 1 >. c AAt) ; Bt ; ld ; (, beta) ; Ct ; ld
)

NB. ---------------------------------------------------------
NB. Dyad         Domain
NB. dtrmmcore    real
NB. ztrmmcore    complex
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
NB.   side  - string, case-insensitive, in which the head
NB.           specifies the side of op(A):
NB.             'L'  NB. to perform (1) (op(A) on the left)
NB.             'R'  NB. to perform (2) (op(A) on the right)
NB.   uplo  - string, case-insensitive, in which the head
NB.           specifies whether the matrix A is upper or
NB.           lower triangular:
NB.             'L'  NB. LT
NB.             'U'  NB. UT
NB.   trans - string, case-insensitive, in which the head
NB.           specifies the form of op(A):
NB.             'N'  NB. op(A) := A    (no transpose)
NB.             'T'  NB. op(A) := A^T  (transpose)
NB.             'C'  NB. op(A) := A^T  (transpose)           for dtrmmcore
NB.                  NB. op(A) := A^H  (conjugate transpose) for ztrmmcore
NB.   diag  - string, case-insensitive, in which the head
NB.           specifies the form of A:
NB.             'N'  NB. A is either L or U
NB.             'U'  NB. A is either L1 or U1, diagonal
NB.                  NB.   elements of A are not referenced
NB.   alpha - scalar
NB.   AAt   - mn×mn-matrix, contains either non-zero or both
NB.           triangular parts of A^T
NB.   Bt    - n×m-matrix, B^T
NB.   Bupdt - an updated Bt
NB.   A     - mn×mn-matrix, triangular
NB.   m     ≥ 0, the number of rows in B
NB.   n     ≥ 0, the number of columns in B
NB.   mn    = m if side='L' or mn = n otherwise
NB.
NB. Notes:
NB. - operate on transposed matrices to avoid transposition

dtrmmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  'side uplo trans diag'=. x
  'alpha AAt Bt'=. y
  'n m'=. $ Bt
  10 {:: dtrmm_cd (, side) ; (, uplo) ; (, trans) ; (, diag) ; (, m) ; (, n) ; (, alpha) ; AAt ; (, 1 >. c AAt) ; Bt ; , 1 >. m
)

ztrmmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  'side uplo trans diag'=. x
  'alpha AAt Bt'=. y
  'n m'=. $ Bt
  10 {:: ztrmm_cd (, side) ; (, uplo) ; (, trans) ; (, diag) ; (, m) ; (, n) ; (, alpha) ; AAt ; (, 1 >. c AAt) ; Bt ; , 1 >. m
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad      Domain     op1(A)    op2(B)
NB. dgemmnn    real       A         B
NB. dgemmnt    real       A         B^T
NB. dgemmtn    real       A^T       B
NB. dgemmtt    real       A^T       B^T
NB. zgemmnn    complex    A         B
NB. zgemmnt    complex    A         B^T
NB. zgemmnc    complex    A         B^H
NB. zgemmtn    complex    A^T       B
NB. zgemmtt    complex    A^T       B^T
NB. zgemmtc    complex    A^T       B^H
NB. zgemmcn    complex    A^H       B
NB. zgemmct    complex    A^H       B^T
NB. zgemmcc    complex    A^H       B^H
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
NB. Monad      Domain     A     Reads in A    Side
NB. dsymmll    real       SY    LT            (1)
NB. dsymmlu    real       SY    UT            (1)
NB. dsymmrl    real       SY    LT            (2)
NB. dsymmru    real       SY    UT            (2)
NB. zsymmll    complex    SY    LT            (1)
NB. zsymmlu    complex    SY    UT            (1)
NB. zsymmrl    complex    SY    LT            (2)
NB. zsymmru    complex    SY    UT            (2)
NB. zhemmll    complex    HE    LT            (1)
NB. zhemmlu    complex    HE    UT            (1)
NB. zhemmrl    complex    HE    LT            (2)
NB. zhemmru    complex    HE    UT            (2)
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * A * B + beta * C  (1)
NB.   or
NB.     C := alpha * B * A + beta * C  (2)
NB.   where A is Hermitian (symmetric)
NB.
NB. Syntax:
NB.   Cupd=. xxxmmxx alpha ; AA ; B ; beta ; C
NB. where
NB.   alpha - scalar
NB.   AA    - mn×mn-matrix, contains either LT or UT or both
NB.           part(s) of A
NB.   B     - m×n-matrix
NB.   beta  - scalar
NB.   C     - m×n-matrix
NB.   Cupd  - an updated C
NB.   A     - mn×mn-matrix, Hermitian (symmetric)
NB.   m     ≥ 0, the number of rows in C, Cupd and B
NB.   n     ≥ 0, the number of columns in C, Cupd and B
NB.   mn    = m for xxxmmlx or mn = n otherwise
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
NB. Monad        Domain     Side    A     Reads in A    op(A)
NB. dtrmmllnn    real       (1)     L      LT           A
NB. dtrmmllnu    real       (1)     L1    SLT           A
NB. dtrmmlltn    real       (1)     L      LT           A^T
NB. dtrmmlltu    real       (1)     L1    SLT           A^T
NB. dtrmmlunn    real       (1)     U      UT           A
NB. dtrmmlunu    real       (1)     U1    SUT           A
NB. dtrmmlutn    real       (1)     U      UT           A^T
NB. dtrmmlutu    real       (1)     U1    SUT           A^T
NB. dtrmmrlnn    real       (2)     L      LT           A
NB. dtrmmrlnu    real       (2)     L1    SLT           A
NB. dtrmmrltn    real       (2)     L      LT           A^T
NB. dtrmmrltu    real       (2)     L1    SLT           A^T
NB. dtrmmrunn    real       (2)     U      UT           A
NB. dtrmmrunu    real       (2)     U1    SUT           A
NB. dtrmmrutn    real       (2)     U      UT           A^T
NB. dtrmmrutu    real       (2)     U1    SUT           A^T
NB. ztrmmllnn    complex    (1)     L      LT           A
NB. ztrmmllnu    complex    (1)     L1    SLT           A
NB. ztrmmlltn    complex    (1)     L      LT           A^T
NB. ztrmmlltu    complex    (1)     L1    SLT           A^T
NB. ztrmmllcn    complex    (1)     L      LT           A^H
NB. ztrmmllcu    complex    (1)     L1    SLT           A^H
NB. ztrmmlunn    complex    (1)     U      UT           A
NB. ztrmmlunu    complex    (1)     U1    SUT           A
NB. ztrmmlutn    complex    (1)     U      UT           A^T
NB. ztrmmlutu    complex    (1)     U1    SUT           A^T
NB. ztrmmlucn    complex    (1)     U      UT           A^H
NB. ztrmmlucu    complex    (1)     U1    SUT           A^H
NB. ztrmmrlnn    complex    (2)     L      LT           A
NB. ztrmmrlnu    complex    (2)     L1    SLT           A
NB. ztrmmrltn    complex    (2)     L      LT           A^T
NB. ztrmmrltu    complex    (2)     L1    SLT           A^T
NB. ztrmmrlcn    complex    (2)     L      LT           A^H
NB. ztrmmrlcu    complex    (2)     L1    SLT           A^H
NB. ztrmmrunn    complex    (2)     U      UT           A
NB. ztrmmrunu    complex    (2)     U1    SUT           A
NB. ztrmmrutn    complex    (2)     U      UT           A^T
NB. ztrmmrutu    complex    (2)     U1    SUT           A^T
NB. ztrmmrucn    complex    (2)     U      UT           A^H
NB. ztrmmrucu    complex    (2)     U1    SUT           A^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     B := alpha * op(A) * B  (1)
NB.   or
NB.     B := alpha * B * op(A)  (2)
NB.   where A is triangular, and op(A) is either A, A^T or
NB.   A^H
NB.
NB. Syntax:
NB.   Bupd=. xtrmmxxxx alpha ; AA ; B
NB. where
NB.   alpha - scalar
NB.   AA    - k×k-matrix, contains either non-zero or both
NB.           part(s) of A
NB.   B     - m×n-matrix
NB.   Bupd  - an updated B
NB.   A     - k×k-matrix, triangular
NB.   m     ≥ 0, the number of rows in B and Bupd
NB.   n     ≥ 0, the number of columns in B and Bupd
NB.   k     = m for xtrmmlxxx or k = n for xtrmmrxxx
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
