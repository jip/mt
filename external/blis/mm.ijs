NB. Matrix-matrix operations
NB.
NB. xgemmxxx    Matrix-matrix operation with general matrix
NB. xsymmxxxx   Matrix-matrix operation with symmetric matrix
NB. xhemmxxxx   Matrix-matrix operation with Hermitian matrix
NB. xtrmmxxxxx  Matrix-matrix operation with triangular matrix
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

coclass 'mtbli'

NB. =========================================================
NB. Includes

require 'math/mt/external/blis/blis'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Dyad          API       Domain
NB. gemmcore      object    mixed
NB. dgemmcore     typed     real
NB. zgemmcore     typed     complex
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C
NB.   where opX(M) is either M, M^T, conj(M) or M^H
NB.
NB. Syntax:
NB.   Cupd=. (transA , transB) xgemmcore alpha ; A ; B ; beta ; C
NB. where
NB.   transA - trans_t scalar, defines the form of op1(A):
NB.              NO_TRANSPOSE       NB. op1(A) := A        (no transpose)
NB.              TRANSPOSE          NB. op1(A) := A^T      (transpose)
NB.              CONJ_NO_TRANSPOSE  NB. op1(A) := conj(A)  (conjugate)
NB.              CONJ_TRANSPOSE     NB. op1(A) := A^H      (conjugate and transpose)
NB.   transB - trans_t scalar, defines the form of op2(B):
NB.              NO_TRANSPOSE       NB. op1(B) := B        (no transpose)
NB.              TRANSPOSE          NB. op1(B) := B^T      (transpose)
NB.              CONJ_NO_TRANSPOSE  NB. op1(B) := conj(B)  (conjugate)
NB.              CONJ_TRANSPOSE     NB. op1(B) := B^H      (conjugate and transpose)
NB.   alpha  - scalar
NB.   A      - ma×ka-matrix
NB.   B      - kb×nb-matrix
NB.   beta   - scalar
NB.   C      - m×n-matrix
NB.   Cupd   - an updated C
NB.   m      ≥ 0, the number of rows in C, Cupd and op1(A)
NB.   n      ≥ 0, the number of columns in C, Cupd and op2(B)
NB.   k      ≥ 0, the number of columns in op1(A) and the
NB.            number of rows in op2(B)
NB.   ma     = m if (op1(A) is either A or conj(A)) or ma = k
NB.            otherwise
NB.   ka     = k if (op1(A) is either A or conj(A)) or ka = m
NB.            otherwise
NB.   kb     = k if (op1(B) is either B or conj(B)) or kb = n
NB.            otherwise
NB.   nb     = n if (op1(B) is either B or conj(B)) or nb = k
NB.            otherwise

gemmcore=: 4 : 0
  y=. 1&memu&.>&.(_1&{) y  NB. unalias C
  objs=. obja L: 0 y       NB. allocate BLIS objects bonded to J nouns
  x obj_set_conjtrans"0 (((;"0&1) 1 2)) {:: objs
    NB. transpose A and/or B optionally
  gemm_cd <@{:L:0 objs     NB. call bli_gemm() with each param marked as object address
  objf L: 0 objs           NB. free BLIS objects
  _1 {:: y                 NB. return changed copy of C
)

dgemmcore=: 4 : 0
  'transA transB'=. x
  'ka nb n'=. 1 2 4 { c S: 0 'alpha A B beta C'=. y
  k=. (transA e. NO_TRANSPOSE , CONJ_NO_TRANSPOSE) { $ A
  14 {:: dgemm_cd transA ; transB ; (# C) ; n ; k ; (, alpha) ; A ; ka ; 1 ; B ; nb ; 1 ; (, beta) ; C ; n ; 1
)

zgemmcore=: 4 : 0
  'transA transB'=. x
  'ka nb n'=. 1 2 4 { c S: 0 'alpha A B beta C'=. y
  k=. (transA e. NO_TRANSPOSE , CONJ_NO_TRANSPOSE) { $ A
  14 {:: zgemm_cd transA ; transB ; (# C) ; n ; k ; (, alpha) ; A ; ka ; 1 ; B ; nb ; 1 ; (, beta) ; C ; n ; 1
)

NB. ---------------------------------------------------------
NB. Dyad          API       Domain
NB. gemmtcore     object    mixed
NB. dgemmtcore    typed     real
NB. zgemmtcore    typed     complex
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C
NB.   updating C in triangular part, where C is square,
NB.   opX(M) is either M, M^T, conj(M) or M^H
NB.
NB. Syntax:
NB.   Cupd=. (transA , transB , uploC) xgemmtcore alpha ; A ; B ; beta ; C
NB. where
NB.   transA - trans_t scalar, defines the form of op1(A):
NB.              NO_TRANSPOSE       NB. op1(A) := A        (no transpose)
NB.              TRANSPOSE          NB. op1(A) := A^T      (transpose)
NB.              CONJ_NO_TRANSPOSE  NB. op1(A) := conj(A)  (conjugate)
NB.              CONJ_TRANSPOSE     NB. op1(A) := A^H      (conjugate and transpose)
NB.   transB - trans_t scalar, defines the form of op2(B):
NB.              NO_TRANSPOSE       NB. op1(B) := B        (no transpose)
NB.              TRANSPOSE          NB. op1(B) := B^T      (transpose)
NB.              CONJ_NO_TRANSPOSE  NB. op1(B) := conj(B)  (conjugate)
NB.              CONJ_TRANSPOSE     NB. op1(B) := B^H      (conjugate and transpose)
NB.   uploC  - uplo_t scalar, defines the part of C to be
NB.            updated:
NB.              LOWER              NB. LT
NB.              UPPER              NB. UT
NB.   alpha  - scalar
NB.   A      - ma×ka-matrix
NB.   B      - kb×mb-matrix
NB.   beta   - scalar
NB.   C      - m×m-matrix
NB.   Cupd   - C updated in part specified by uploC
NB.   m      ≥ 0, the size of C and Cupd, the number of rows
NB.            in op1(A) and the number of columns in op2(B)
NB.   k      ≥ 0, the number of columns in op1(A) and the
NB.            number of rows in op2(B)
NB.   ma     = m if (op1(A) is either A or conj(A)) or ma = k
NB.            otherwise
NB.   ka     = k if (op1(A) is either A or conj(A)) or ka = m
NB.            otherwise
NB.   kb     = k if (op1(B) is either B or conj(B)) or kb = m
NB.            otherwise
NB.   mb     = m if (op1(B) is either B or conj(B)) or mb = k
NB.            otherwise

gemmtcore=: (4 : 0) ([ assert@basiccr4)
  y=. 1&memu&.>&.(_1&{) y  NB. unalias C
  objs=. obja L: 0 y       NB. allocate BLIS objects bonded to J nouns
  x (2 1 # obj_set_conjtrans`obj_set_uplo)"0 (((;"0&1) 1 2 4)) {:: objs
    NB. transpose A and/or B optionally, select C part
  gemmt_cd <@{:L:0 objs    NB. call bli_gemmt() with each param marked as object address
  objf L: 0 objs           NB. free BLIS objects
  _1 {:: y                 NB. return changed copy of C
)

dgemmtcore=: (4 : 0) ([ assert@basiccr4)
  'transA transB uploC'=. x
  'ka mb m'=. 1 2 4 { c S: 0 'alpha A B beta C'=. y
  k=. (transA e. NO_TRANSPOSE , CONJ_NO_TRANSPOSE) { $ A
  14 {:: dgemmt_cd uploC ; transA ; transB ; m ; k ; (, alpha) ; A ; ka ; 1 ; B ; mb ; 1 ; (, beta) ; C ; m ; 1
)

zgemmtcore=: (4 : 0) ([ assert@basiccr4)
  'transA transB uploC'=. x
  'ka mb m'=. 1 2 4 { c S: 0 'alpha A B beta C'=. y
  k=. (transA e. NO_TRANSPOSE , CONJ_NO_TRANSPOSE) { $ A
  14 {:: zgemmt_cd uploC ; transA ; transB ; m ; k ; (, alpha) ; A ; ka ; 1 ; B ; mb ; 1 ; (, beta) ; C ; m ; 1
)

NB. ---------------------------------------------------------
NB. Dyad          API       Domain     A
NB. symmcore      object    mixed      SY
NB. hemmcore      object    mixed      HE
NB. dsymmcore     typed     real       SY
NB. zsymmcore     typed     complex    SY
NB. zhemmcore     typed     complex    HE
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C  (1)
NB.   or
NB.     C := alpha * op2(B) * op1(A) + beta * C  (2)
NB.   where A is Hermitian (symmetric), op1(A) is either A or
NB.   conj(A), and op2(B) is either B, B^T, conj(B) or B^H
NB.
NB. Syntax:
NB.   Cupd=. (sideA , uploA , conjA , transB) xxxmmcore alpha ; AA ; B ; beta ; C
NB. where
NB.   sideA  - side_t scalar, defines the side of A:
NB.              LEFT               NB. to perform (1) (A on the left)
NB.              RIGHT              NB. to perform (2) (A on the right)
NB.   uploA  - uplo_t scalar, defines the part of A to be
NB.            read:
NB.              LOWER              NB. LT
NB.              UPPER              NB. UT
NB.   conjA  - conj_t scalar, defines the form of op1(A):
NB.              NO_CONJUGATE       NB. op1(A) := A        (no conjugate)
NB.              CONJUGATE          NB. op1(A) := conj(A)  (conjugate)
NB.   transB - trans_t scalar, defines the form of op2(B):
NB.              NO_TRANSPOSE       NB. op1(B) := B        (no transpose)
NB.              TRANSPOSE          NB. op1(B) := B^T      (transpose)
NB.              CONJ_NO_TRANSPOSE  NB. op1(B) := conj(B)  (conjugate)
NB.              CONJ_TRANSPOSE     NB. op1(B) := B^H      (conjugate and transpose)
NB.   alpha  - scalar
NB.   AA     - mn×mn-matrix, contains either LT or UT or both
NB.            part(s) of A
NB.   B      - m×n-matrix
NB.   beta   - scalar
NB.   C      - m×n-matrix
NB.   Cupd   - an updated C
NB.   A      - mn×mn-matrix, Hermitian (symmetric)
NB.   m      ≥ 0, the number of rows in B, C and Cupd, the
NB.            size of A and AA when (sideA == LEFT)
NB.   n      ≥ 0, the number of columns in B, C and Cupd, the
NB.            size of A and AA when (sideA != LEFT)
NB.   mn     = m if (sideA == LEFT) or mn = n otherwise

symmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  y=. 1&memu&.>&.(_1&{) y        NB. unalias C
  objs=. obja L: 0 y             NB. allocate BLIS objects bonded to J nouns
  (SYMMETRIC 0} x) obj_set_struc`obj_set_uplo`obj_set_conj`obj_set_conjtrans"0 (((;"0&1) 1 1 1 2)) {:: objs
    NB. set A structure, select A part, conjugate A optionally, transpose B optionally
    NB. note: changing the object is being the op with side-effect so
    NB.       it is possible to do it in parallel as here
  symm_cd ({. x) ; <@{:L:0 objs  NB. call bli_symm() with all but head params marked as object address
  objf L: 0 objs                 NB. free BLIS objects
  _1 {:: y                       NB. return changed copy of C
)

hemmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  y=. 1&memu&.>&.(_1&{) y        NB. unalias C
  objs=. obja L: 0 y             NB. allocate BLIS objects bonded to J nouns
  (HERMITIAN 0} x) obj_set_struc`obj_set_uplo`obj_set_conj`obj_set_conjtrans"0 (((;"0&1) 1 1 1 2)) {:: objs
    NB. set A structure, select A part, conjugate A optionally, transpose B optionally
    NB. note: changing the object is being the op with side-effect so
    NB.       it is possible to do it in parallel as here
  hemm_cd ({. x) ; <@{:L:0 objs  NB. call bli_hemm() with all but head params marked as object address
  objf L: 0 objs                 NB. free BLIS objects
  _1 {:: y                       NB. return changed copy of C
)

dsymmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  'sideA uploA conjA transB'=. x
  'alpha AA B beta C'=. y
  mn=. (RIGHT = sideA) { 'm n'=. $ C
  15 {:: dsymm_cd sideA ; uploA ; conjA ; transB ; m ; n ; (, alpha) ; AA ; mn ; 1 ; B ; n ; 1 ; (, beta) ; C ; n ; 1
)

zsymmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  'sideA uploA conjA transB'=. x
  'alpha AA B beta C'=. y
  mn=. (RIGHT = sideA) { 'm n'=. $ C
  15 {:: zsymm_cd sideA ; uploA ; conjA ; transB ; m ; n ; (, alpha) ; AA ; mn ; 1 ; B ; n ; 1 ; (, beta) ; C ; n ; 1
)

zhemmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  'sideA uploA conjA transB'=. x
  'alpha AA B beta C'=. y
  mn=. (RIGHT = sideA) { 'm n'=. $ C
  15 {:: zhemm_cd sideA ; uploA ; conjA ; transB ; m ; n ; (, alpha) ; AA ; mn ; 1 ; B ; n ; 1 ; (, beta) ; C ; n ; 1
)

NB. ---------------------------------------------------------
NB. Dyad         API        Domain
NB. trmmcore     object     mixed
NB. dtrmmcore    typed      real
NB. ztrmmcore    typed      complex
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     B := alpha * op(A) * B  (1)
NB.   or
NB.     B := alpha * B * op(A)  (2)
NB.   where A is triangular, and op(A) is either A, A^T,
NB.   conj(A) or A^H
NB.
NB. Syntax:
NB.   Bupd=. (sideA , uploA , transA , diagA) xtrmmcore alpha ; AA ; B
NB. where
NB.   sideA  - side_t scalar, defines the side of A:
NB.              LEFT               NB. to perform (1) (A on the left)
NB.              RIGHT              NB. to perform (2) (A on the right)
NB.   uploA  - uplo_t scalar, defines the part of A to be
NB.            read:
NB.              LOWER              NB. LT
NB.              UPPER              NB. UT
NB.   transA - trans_t scalar, defines the form of op(A):
NB.              NO_TRANSPOSE       NB. op(A) := A        (no transpose)
NB.              TRANSPOSE          NB. op(A) := A^T      (transpose)
NB.              CONJ_NO_TRANSPOSE  NB. op(A) := conj(A)  (conjugate)
NB.              CONJ_TRANSPOSE     NB. op(A) := A^H      (conjugate and transpose)
NB.   diagA  - diag_t scalar, defines the form of A:
NB.              NONUNIT_DIAG       NB. A is either L or U
NB.              UNIT_DIAG          NB. A is either L1 or U1, diagonal
NB.                                 NB.   elements of A are not referenced
NB.   alpha  - scalar
NB.   AA     - mn×mn-matrix, contains either non-zero or both
NB.            part(s) of A
NB.   B      - m×n-matrix
NB.   Bupd   - an updated B
NB.   A      - mn×mn-matrix, triangular
NB.   m      ≥ 0, the number of rows in B and Bupd, the size
NB.            of A and AA when (sideA == LEFT)
NB.   n      ≥ 0, the number of columns in B and Bupd, the
NB.            size of A and AA when (sideA != LEFT)
NB.   mn     = m if (sideA == LEFT) or mn = n otherwise

trmmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  y=. 1&memu&.>&.(_1&{) y        NB. unalias B
  objs=. obja L: 0 y             NB. allocate BLIS objects bonded to J nouns
  (TRIANGULAR 0} x) obj_set_struc`obj_set_uplo`obj_set_conjtrans`obj_set_diag"0 (1;1) {:: objs
    NB. set A structure, select A part, transpose A optionally, set A diag type
    NB. note: changing the object is being the op with side-effect so
    NB.       it is possible to do it in parallel as here
  trmm_cd ({. x) ; <@{:L:0 objs  NB. call bli_trmm() with all but head params marked as object address
  objf L: 0 objs                 NB. free BLIS objects
  _1 {:: y                       NB. return changed copy of B
)

dtrmmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  'sideA uploA transA diagA'=. x
  'alpha AA B'=. y
  mn=. (RIGHT = sideA) { 'm n'=. $ B
  11 {:: dtrmm_cd sideA ; uploA ; transA ; diagA ; m ; n ; (, alpha) ; AA ; mn ; 1 ; B ; n ; 1
)

ztrmmcore=: (4 : 0) ([ assert@(basiccs1 , basiccr0))
  'sideA uploA transA diagA'=. x
  'alpha AA B'=. y
  mn=. (RIGHT = sideA) { 'm n'=. $ B
  11 {:: ztrmm_cd sideA ; uploA ; transA ; diagA ; m ; n ; (, alpha) ; AA ; mn ; 1 ; B ; n ; 1
)

NB. ---------------------------------------------------------
NB. Dyad          API        Domain
NB. trmm3core     object     mixed
NB. dtrmm3core    typed      real
NB. ztrmm3core    typed      complex
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C  (1)
NB.   or
NB.     C := alpha * op2(B) * op1(A) + beta * C  (2)
NB.   where A is triangular, and opX(M) is either M, M^T,
NB.   conj(M) or M^H
NB.
NB. Syntax:
NB.   Cupd=. (sideA , uploA , transA , diagA , transB) xtrmm3core alpha ; AA ; B ; beta ; C
NB. where
NB.   sideA  - side_t scalar, defines the side of A:
NB.              LEFT               NB. to perform (1) (A on the left)
NB.              RIGHT              NB. to perform (2) (A on the right)
NB.   uploA  - uplo_t scalar, defines the part of A to be read:
NB.              LOWER              NB. LT
NB.              UPPER              NB. UT
NB.   transA - trans_t scalar, defines the form of op1(A):
NB.              NO_TRANSPOSE       NB. op1(A) := A        (no transpose)
NB.              TRANSPOSE          NB. op1(A) := A^T      (transpose)
NB.              CONJ_NO_TRANSPOSE  NB. op1(A) := conj(A)  (conjugate)
NB.              CONJ_TRANSPOSE     NB. op1(A) := A^H      (conjugate and transpose)
NB.   diagA  - diag_t scalar, defines the form of A:
NB.              NONUNIT_DIAG       NB. A is either L or U
NB.              UNIT_DIAG          NB. A is either L1 or U1, diagonal
NB.                                 NB.   elements of A are not referenced
NB.   transB - trans_t scalar, defines the form of op2(B):
NB.              NO_TRANSPOSE       NB. op2(B) := B        (no transpose)
NB.              TRANSPOSE          NB. op2(B) := B^T      (transpose)
NB.              CONJ_NO_TRANSPOSE  NB. op2(B) := conj(B)  (conjugate)
NB.              CONJ_TRANSPOSE     NB. op2(B) := B^H      (conjugate and transpose)
NB.   alpha  - scalar
NB.   AA     - mn×mn-matrix, contains either non-zero or both
NB.            part(s) of A
NB.   B      - mb×nb-matrix
NB.   beta   - scalar
NB.   C      - m×n-matrix
NB.   Cupd   - an updated C
NB.   A      - mn×mn-matrix, triangular
NB.   m      ≥ 0, the size of A and AA for trmmlxxxx, and the
NB.            number of rows in op2(B), C and Cupd
NB.   n      ≥ 0, the size of A and AA for trmmrxxxx, and the
NB.            number of columns in op2(B), C and Cupd
NB.   mn     = m if (sideA == LEFT) or mn = n otherwise
NB.   mb     = m if (op1(B) is either B or conj(B)) or mb = n
NB.            otherwise
NB.   nb     = n if (op1(B) is either B or conj(B)) or nb = m
NB.            otherwise

trmm3core=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  y=. 1&memu&.>&.(_1&{) y         NB. unalias C
  objs=. obja L: 0 y              NB. allocate BLIS objects bonded to J nouns
  (TRIANGULAR 0} x) obj_set_struc`obj_set_uplo`obj_set_conjtrans`obj_set_diag`obj_set_conjtrans"0 (((;"0&1) 4 1 # 1 2)) {:: objs
    NB. set A structure, select A part, transpose A optionally, set A diag type, transpose B optionally
    NB. note: changing the object is being the op with side-effect so
    NB.       it is possible to do it in parallel as here
  trmm3_cd ({. x) ; <@{:L:0 objs  NB. call bli_trmm() with all but head params marked as object address
  objf L: 0 objs                  NB. free BLIS objects
  _1 {:: y                        NB. return changed copy of C
)

dtrmm3core=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  'sideA uploA transA diagA transB'=. x
  'alpha AA B beta C'=. y
  mn=. (RIGHT = sideA) { 'm n'=. $ C
  nb=. (transB e. NO_TRANSPOSE , CONJ_NO_TRANSPOSE) { $ B
  16 {:: dtrmm3_cd sideA ; uploA ; transA ; diagA ; transB ; m ; n ; (, alpha) ; AA ; mn ; 1 ; B ; nb ; 1 ; (, beta) ; C ; n ; 1
)

ztrmm3core=: (4 : 0) ([ assert@(basiccs1 , basiccr4))
  'sideA uploA transA diagA transB'=. x
  'alpha AA B beta C'=. y
  mn=. (RIGHT = sideA) { 'm n'=. $ C
  nb=. (transB e. NO_TRANSPOSE , CONJ_NO_TRANSPOSE) { $ B
  16 {:: ztrmm3_cd sideA ; uploA ; transA ; diagA ; transB ; m ; n ; (, alpha) ; AA ; mn ; 1 ; B ; nb ; 1 ; (, beta) ; C ; n ; 1
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Monad      API       Domain     op1(A)     op2(B)
NB. gemmnn     object    mixed      A          B
NB. gemmnt     object    mixed      A          B^T
NB. gemmnj     object    mixed      A          conj(B)
NB. gemmnc     object    mixed      A          B^H
NB. gemmtn     object    mixed      A^T        B
NB. gemmtt     object    mixed      A^T        B^T
NB. gemmtj     object    mixed      A^T        conj(B)
NB. gemmtc     object    mixed      A^T        B^H
NB. gemmjn     object    mixed      conj(A)    B
NB. gemmjt     object    mixed      conj(A)    B^T
NB. gemmjj     object    mixed      conj(A)    conj(B)
NB. gemmjc     object    mixed      conj(A)    B^H
NB. gemmcn     object    mixed      A^H        B
NB. gemmct     object    mixed      A^H        B^T
NB. gemmcj     object    mixed      A^H        conj(B)
NB. gemmcc     object    mixed      A^H        B^H
NB. dgemmnn    typed     real       A          B
NB. dgemmnt    typed     real       A          B^T
NB. dgemmtn    typed     real       A^T        B
NB. dgemmtt    typed     real       A^T        B^T
NB. zgemmnn    typed     complex    A          B
NB. zgemmnt    typed     complex    A          B^T
NB. zgemmnj    typed     complex    A          conj(B)
NB. zgemmnc    typed     complex    A          B^H
NB. zgemmtn    typed     complex    A^T        B
NB. zgemmtt    typed     complex    A^T        B^T
NB. zgemmtj    typed     complex    A^T        conj(B)
NB. zgemmtc    typed     complex    A^T        B^H
NB. zgemmjn    typed     complex    conj(A)    B
NB. zgemmjt    typed     complex    conj(A)    B^T
NB. zgemmjj    typed     complex    conj(A)    conj(B)
NB. zgemmjc    typed     complex    conj(A)    B^H
NB. zgemmcn    typed     complex    A^H        B
NB. zgemmct    typed     complex    A^H        B^T
NB. zgemmcj    typed     complex    A^H        conj(B)
NB. zgemmcc    typed     complex    A^H        B^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C
NB.   where opX(M) is either M, M^T, conj(M) or M^H
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
NB.   m     ≥ 0, the number of rows in C and Cupd, the number
NB.           of rows (for xgemmnx and xgemmjx) or columns
NB.           (for xgemmtx and xgemmcx) in A
NB.   n     ≥ 0, the number of columns in C and Cupd, the
NB.           number of rows (for xgemmxt and xgemmxc) or
NB.           columns (for xgemmxn and xgemmxj) in B
NB.   k     ≥ 0, the number of rows (for xgemmtx and xgemmcx)
NB.           or columns (for xgemmnx and xgemmjx) in A, the
NB.           number of rows (for xgemmxn and xgemmxj) or
NB.           columns (for xgemmxt and xgemmxc) in B
NB.   ma    = m for xgemmnx and xgemmjx or ma = k otherwise
NB.   ka    = k for xgemmnx and xgemmjx or ka = m otherwise
NB.   kb    = k for xgemmxn and xgemmxj or kb = n otherwise
NB.   nb    = n for xgemmxn and xgemmxj or nb = k otherwise
NB.
NB. Notes:
NB. - monad      provides BLIS'
NB.   gemmxx     bli_gemm (...)
NB.   dgemmxx    bli_dgemm(...)
NB.   zgemmnt    bli_zgemm(...)

gemmnn=:  (     NO_TRANSPOSE ,      NO_TRANSPOSE)& gemmcore
gemmnt=:  (     NO_TRANSPOSE ,         TRANSPOSE)& gemmcore
gemmnj=:  (     NO_TRANSPOSE , CONJ_NO_TRANSPOSE)& gemmcore
gemmnc=:  (     NO_TRANSPOSE ,    CONJ_TRANSPOSE)& gemmcore
gemmtn=:  (        TRANSPOSE ,      NO_TRANSPOSE)& gemmcore
gemmtt=:  (        TRANSPOSE ,         TRANSPOSE)& gemmcore
gemmtj=:  (        TRANSPOSE , CONJ_NO_TRANSPOSE)& gemmcore
gemmtc=:  (        TRANSPOSE ,    CONJ_TRANSPOSE)& gemmcore
gemmjn=:  (CONJ_NO_TRANSPOSE ,      NO_TRANSPOSE)& gemmcore
gemmjt=:  (CONJ_NO_TRANSPOSE ,         TRANSPOSE)& gemmcore
gemmjj=:  (CONJ_NO_TRANSPOSE , CONJ_NO_TRANSPOSE)& gemmcore
gemmjc=:  (CONJ_NO_TRANSPOSE ,    CONJ_TRANSPOSE)& gemmcore
gemmcn=:  (   CONJ_TRANSPOSE ,      NO_TRANSPOSE)& gemmcore
gemmct=:  (   CONJ_TRANSPOSE ,         TRANSPOSE)& gemmcore
gemmcj=:  (   CONJ_TRANSPOSE , CONJ_NO_TRANSPOSE)& gemmcore
gemmcc=:  (   CONJ_TRANSPOSE ,    CONJ_TRANSPOSE)& gemmcore

dgemmnn=: (     NO_TRANSPOSE ,      NO_TRANSPOSE)&dgemmcore
dgemmnt=: (     NO_TRANSPOSE ,         TRANSPOSE)&dgemmcore
dgemmtn=: (        TRANSPOSE ,      NO_TRANSPOSE)&dgemmcore
dgemmtt=: (        TRANSPOSE ,         TRANSPOSE)&dgemmcore

zgemmnn=: (     NO_TRANSPOSE ,      NO_TRANSPOSE)&zgemmcore
zgemmnt=: (     NO_TRANSPOSE ,         TRANSPOSE)&zgemmcore
zgemmnj=: (     NO_TRANSPOSE , CONJ_NO_TRANSPOSE)&zgemmcore
zgemmnc=: (     NO_TRANSPOSE ,    CONJ_TRANSPOSE)&zgemmcore
zgemmtn=: (        TRANSPOSE ,      NO_TRANSPOSE)&zgemmcore
zgemmtt=: (        TRANSPOSE ,         TRANSPOSE)&zgemmcore
zgemmtj=: (        TRANSPOSE , CONJ_NO_TRANSPOSE)&zgemmcore
zgemmtc=: (        TRANSPOSE ,    CONJ_TRANSPOSE)&zgemmcore
zgemmjn=: (CONJ_NO_TRANSPOSE ,      NO_TRANSPOSE)&zgemmcore
zgemmjt=: (CONJ_NO_TRANSPOSE ,         TRANSPOSE)&zgemmcore
zgemmjj=: (CONJ_NO_TRANSPOSE , CONJ_NO_TRANSPOSE)&zgemmcore
zgemmjc=: (CONJ_NO_TRANSPOSE ,    CONJ_TRANSPOSE)&zgemmcore
zgemmcn=: (   CONJ_TRANSPOSE ,      NO_TRANSPOSE)&zgemmcore
zgemmct=: (   CONJ_TRANSPOSE ,         TRANSPOSE)&zgemmcore
zgemmcj=: (   CONJ_TRANSPOSE , CONJ_NO_TRANSPOSE)&zgemmcore
zgemmcc=: (   CONJ_TRANSPOSE ,    CONJ_TRANSPOSE)&zgemmcore

NB. ---------------------------------------------------------
NB. Monad       API       Domain     R/W in C    op1(A)     op2(B)
NB. gemmlnn     object    mixed      LT          A          B
NB. gemmlnt     object    mixed      LT          A          B^T
NB. gemmlnj     object    mixed      LT          A          conj(B)
NB. gemmlnc     object    mixed      LT          A          B^H
NB. gemmltn     object    mixed      LT          A^T        B
NB. gemmltt     object    mixed      LT          A^T        B^T
NB. gemmltj     object    mixed      LT          A^T        conj(B)
NB. gemmltc     object    mixed      LT          A^T        B^H
NB. gemmljn     object    mixed      LT          conj(A)    B
NB. gemmljt     object    mixed      LT          conj(A)    B^T
NB. gemmljj     object    mixed      LT          conj(A)    conj(B)
NB. gemmljc     object    mixed      LT          conj(A)    B^H
NB. gemmlcn     object    mixed      LT          A^H        B
NB. gemmlct     object    mixed      LT          A^H        B^T
NB. gemmlcj     object    mixed      LT          A^H        conj(B)
NB. gemmlcc     object    mixed      LT          A^H        B^H
NB. gemmunn     object    mixed      UT          A          B
NB. gemmunt     object    mixed      UT          A          B^T
NB. gemmunj     object    mixed      UT          A          conj(B)
NB. gemmunc     object    mixed      UT          A          B^H
NB. gemmutn     object    mixed      UT          A^T        B
NB. gemmutt     object    mixed      UT          A^T        B^T
NB. gemmutj     object    mixed      UT          A^T        conj(B)
NB. gemmutc     object    mixed      UT          A^T        B^H
NB. gemmujn     object    mixed      UT          conj(A)    B
NB. gemmujt     object    mixed      UT          conj(A)    B^T
NB. gemmujj     object    mixed      UT          conj(A)    conj(B)
NB. gemmujc     object    mixed      UT          conj(A)    B^H
NB. gemmucn     object    mixed      UT          A^H        B
NB. gemmuct     object    mixed      UT          A^H        B^T
NB. gemmucj     object    mixed      UT          A^H        conj(B)
NB. gemmucc     object    mixed      UT          A^H        B^H
NB. dgemmlnn    typed     real       LT          A          B
NB. dgemmlnt    typed     real       LT          A          B^T
NB. dgemmltn    typed     real       LT          A^T        B
NB. dgemmltt    typed     real       LT          A^T        B^T
NB. dgemmunn    typed     real       UT          A          B
NB. dgemmunt    typed     real       UT          A          B^T
NB. dgemmutn    typed     real       UT          A^T        B
NB. dgemmutt    typed     real       UT          A^T        B^T
NB. zgemmlnn    typed     complex    LT          A          B
NB. zgemmlnt    typed     complex    LT          A          B^T
NB. zgemmlnj    typed     complex    LT          A          conj(B)
NB. zgemmlnc    typed     complex    LT          A          B^H
NB. zgemmltn    typed     complex    LT          A^T        B
NB. zgemmltt    typed     complex    LT          A^T        B^T
NB. zgemmltj    typed     complex    LT          A^T        conj(B)
NB. zgemmltc    typed     complex    LT          A^T        B^H
NB. zgemmljn    typed     complex    LT          conj(A)    B
NB. zgemmljt    typed     complex    LT          conj(A)    B^T
NB. zgemmljj    typed     complex    LT          conj(A)    conj(B)
NB. zgemmljc    typed     complex    LT          conj(A)    B^H
NB. zgemmlcn    typed     complex    LT          A^H        B
NB. zgemmlct    typed     complex    LT          A^H        B^T
NB. zgemmlcj    typed     complex    LT          A^H        conj(B)
NB. zgemmlcc    typed     complex    LT          A^H        B^H
NB. zgemmunn    typed     complex    UT          A          B
NB. zgemmunt    typed     complex    UT          A          B^T
NB. zgemmunj    typed     complex    UT          A          conj(B)
NB. zgemmunc    typed     complex    UT          A          B^H
NB. zgemmutn    typed     complex    UT          A^T        B
NB. zgemmutt    typed     complex    UT          A^T        B^T
NB. zgemmutj    typed     complex    UT          A^T        conj(B)
NB. zgemmutc    typed     complex    UT          A^T        B^H
NB. zgemmujn    typed     complex    UT          conj(A)    B
NB. zgemmujt    typed     complex    UT          conj(A)    B^T
NB. zgemmujj    typed     complex    UT          conj(A)    conj(B)
NB. zgemmujc    typed     complex    UT          conj(A)    B^H
NB. zgemmucn    typed     complex    UT          A^H        B
NB. zgemmuct    typed     complex    UT          A^H        B^T
NB. zgemmucj    typed     complex    UT          A^H        conj(B)
NB. zgemmucc    typed     complex    UT          A^H        B^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C
NB.   updating C in triangular part, where opX(M) is either
NB.   M, M^T, conj(M) or M^H
NB.
NB. Syntax:
NB.   Cupd=. xgemmxxx alpha ; A ; B ; beta ; C
NB. where
NB.   alpha - scalar
NB.   A     - ma×ka-matrix
NB.   B     - kb×mb-matrix
NB.   beta  - scalar
NB.   C     - m×m-matrix
NB.   Cupd  - C updated either in LT (for xgemmlxx) or in UT
NB.           (for xgemmuxx)
NB.   m     ≥ 0, the size of C and Cupd, the number of rows
NB.           (for xgemmxnx and xgemmxjx) or columns (for
NB.           xgemmxtx and xgemmxcx) in A, the number of rows
NB.           (for xgemmxxt and xgemmxxc) or columns (for
NB.           xgemmxxn and xgemmxxj) in B
NB.   k     ≥ 0, the number of rows (for xgemmxtx and
NB.           xgemmxcx) or columns (for xgemmxnx and
NB.           xgemmxjx) in A, the number of rows (for
NB.           xgemmxxn and xgemmxxj) or columns (for xgemmxxt
NB.           and xgemmxxc) in B
NB.   ma    = m for xgemmxnx and xgemmxjx or ma = k otherwise
NB.   ka    = k for xgemmxnx and xgemmxjx or ka = m otherwise
NB.   kb    = k for xgemmxxn and xgemmxxj or kb = m otherwise
NB.   mb    = m for xgemmxxn and xgemmxxj or mb = k otherwise
NB.
NB. Notes:
NB. - monad       provides BLIS'
NB.   gemmxxx     bli_gemmt (...)
NB.   dgemmxxx    bli_dgemmt(...)
NB.   zgemmxxx    bli_zgemmt(...)

gemmlnn=:  (     NO_TRANSPOSE ,      NO_TRANSPOSE , LOWER)& gemmtcore
gemmlnt=:  (     NO_TRANSPOSE ,         TRANSPOSE , LOWER)& gemmtcore
gemmlnj=:  (     NO_TRANSPOSE , CONJ_NO_TRANSPOSE , LOWER)& gemmtcore
gemmlnc=:  (     NO_TRANSPOSE ,    CONJ_TRANSPOSE , LOWER)& gemmtcore
gemmltn=:  (        TRANSPOSE ,      NO_TRANSPOSE , LOWER)& gemmtcore
gemmltt=:  (        TRANSPOSE ,         TRANSPOSE , LOWER)& gemmtcore
gemmltj=:  (        TRANSPOSE , CONJ_NO_TRANSPOSE , LOWER)& gemmtcore
gemmltc=:  (        TRANSPOSE ,    CONJ_TRANSPOSE , LOWER)& gemmtcore
gemmljn=:  (CONJ_NO_TRANSPOSE ,      NO_TRANSPOSE , LOWER)& gemmtcore
gemmljt=:  (CONJ_NO_TRANSPOSE ,         TRANSPOSE , LOWER)& gemmtcore
gemmljj=:  (CONJ_NO_TRANSPOSE , CONJ_NO_TRANSPOSE , LOWER)& gemmtcore
gemmljc=:  (CONJ_NO_TRANSPOSE ,    CONJ_TRANSPOSE , LOWER)& gemmtcore
gemmlcn=:  (   CONJ_TRANSPOSE ,      NO_TRANSPOSE , LOWER)& gemmtcore
gemmlct=:  (   CONJ_TRANSPOSE ,         TRANSPOSE , LOWER)& gemmtcore
gemmlcj=:  (   CONJ_TRANSPOSE , CONJ_NO_TRANSPOSE , LOWER)& gemmtcore
gemmlcc=:  (   CONJ_TRANSPOSE ,    CONJ_TRANSPOSE , LOWER)& gemmtcore
gemmunn=:  (     NO_TRANSPOSE ,      NO_TRANSPOSE , UPPER)& gemmtcore
gemmunt=:  (     NO_TRANSPOSE ,         TRANSPOSE , UPPER)& gemmtcore
gemmunj=:  (     NO_TRANSPOSE , CONJ_NO_TRANSPOSE , UPPER)& gemmtcore
gemmunc=:  (     NO_TRANSPOSE ,    CONJ_TRANSPOSE , UPPER)& gemmtcore
gemmutn=:  (        TRANSPOSE ,      NO_TRANSPOSE , UPPER)& gemmtcore
gemmutt=:  (        TRANSPOSE ,         TRANSPOSE , UPPER)& gemmtcore
gemmutj=:  (        TRANSPOSE , CONJ_NO_TRANSPOSE , UPPER)& gemmtcore
gemmutc=:  (        TRANSPOSE ,    CONJ_TRANSPOSE , UPPER)& gemmtcore
gemmujn=:  (CONJ_NO_TRANSPOSE ,      NO_TRANSPOSE , UPPER)& gemmtcore
gemmujt=:  (CONJ_NO_TRANSPOSE ,         TRANSPOSE , UPPER)& gemmtcore
gemmujj=:  (CONJ_NO_TRANSPOSE , CONJ_NO_TRANSPOSE , UPPER)& gemmtcore
gemmujc=:  (CONJ_NO_TRANSPOSE ,    CONJ_TRANSPOSE , UPPER)& gemmtcore
gemmucn=:  (   CONJ_TRANSPOSE ,      NO_TRANSPOSE , UPPER)& gemmtcore
gemmuct=:  (   CONJ_TRANSPOSE ,         TRANSPOSE , UPPER)& gemmtcore
gemmucj=:  (   CONJ_TRANSPOSE , CONJ_NO_TRANSPOSE , UPPER)& gemmtcore
gemmucc=:  (   CONJ_TRANSPOSE ,    CONJ_TRANSPOSE , UPPER)& gemmtcore

dgemmlnn=: (     NO_TRANSPOSE ,      NO_TRANSPOSE , LOWER)&dgemmtcore
dgemmlnt=: (     NO_TRANSPOSE ,         TRANSPOSE , LOWER)&dgemmtcore
dgemmltn=: (        TRANSPOSE ,      NO_TRANSPOSE , LOWER)&dgemmtcore
dgemmltt=: (        TRANSPOSE ,         TRANSPOSE , LOWER)&dgemmtcore
dgemmunn=: (     NO_TRANSPOSE ,      NO_TRANSPOSE , UPPER)&dgemmtcore
dgemmunt=: (     NO_TRANSPOSE ,         TRANSPOSE , UPPER)&dgemmtcore
dgemmutn=: (        TRANSPOSE ,      NO_TRANSPOSE , UPPER)&dgemmtcore
dgemmutt=: (        TRANSPOSE ,         TRANSPOSE , UPPER)&dgemmtcore

zgemmlnn=: (     NO_TRANSPOSE ,      NO_TRANSPOSE , LOWER)&zgemmtcore
zgemmlnt=: (     NO_TRANSPOSE ,         TRANSPOSE , LOWER)&zgemmtcore
zgemmlnj=: (     NO_TRANSPOSE , CONJ_NO_TRANSPOSE , LOWER)&zgemmtcore
zgemmlnc=: (     NO_TRANSPOSE ,    CONJ_TRANSPOSE , LOWER)&zgemmtcore
zgemmltn=: (        TRANSPOSE ,      NO_TRANSPOSE , LOWER)&zgemmtcore
zgemmltt=: (        TRANSPOSE ,         TRANSPOSE , LOWER)&zgemmtcore
zgemmltj=: (        TRANSPOSE , CONJ_NO_TRANSPOSE , LOWER)&zgemmtcore
zgemmltc=: (        TRANSPOSE ,    CONJ_TRANSPOSE , LOWER)&zgemmtcore
zgemmljn=: (CONJ_NO_TRANSPOSE ,      NO_TRANSPOSE , LOWER)&zgemmtcore
zgemmljt=: (CONJ_NO_TRANSPOSE ,         TRANSPOSE , LOWER)&zgemmtcore
zgemmljj=: (CONJ_NO_TRANSPOSE , CONJ_NO_TRANSPOSE , LOWER)&zgemmtcore
zgemmljc=: (CONJ_NO_TRANSPOSE ,    CONJ_TRANSPOSE , LOWER)&zgemmtcore
zgemmlcn=: (   CONJ_TRANSPOSE ,      NO_TRANSPOSE , LOWER)&zgemmtcore
zgemmlct=: (   CONJ_TRANSPOSE ,         TRANSPOSE , LOWER)&zgemmtcore
zgemmlcj=: (   CONJ_TRANSPOSE , CONJ_NO_TRANSPOSE , LOWER)&zgemmtcore
zgemmlcc=: (   CONJ_TRANSPOSE ,    CONJ_TRANSPOSE , LOWER)&zgemmtcore
zgemmunn=: (     NO_TRANSPOSE ,      NO_TRANSPOSE , UPPER)&zgemmtcore
zgemmunt=: (     NO_TRANSPOSE ,         TRANSPOSE , UPPER)&zgemmtcore
zgemmunj=: (     NO_TRANSPOSE , CONJ_NO_TRANSPOSE , UPPER)&zgemmtcore
zgemmunc=: (     NO_TRANSPOSE ,    CONJ_TRANSPOSE , UPPER)&zgemmtcore
zgemmutn=: (        TRANSPOSE ,      NO_TRANSPOSE , UPPER)&zgemmtcore
zgemmutt=: (        TRANSPOSE ,         TRANSPOSE , UPPER)&zgemmtcore
zgemmutj=: (        TRANSPOSE , CONJ_NO_TRANSPOSE , UPPER)&zgemmtcore
zgemmutc=: (        TRANSPOSE ,    CONJ_TRANSPOSE , UPPER)&zgemmtcore
zgemmujn=: (CONJ_NO_TRANSPOSE ,      NO_TRANSPOSE , UPPER)&zgemmtcore
zgemmujt=: (CONJ_NO_TRANSPOSE ,         TRANSPOSE , UPPER)&zgemmtcore
zgemmujj=: (CONJ_NO_TRANSPOSE , CONJ_NO_TRANSPOSE , UPPER)&zgemmtcore
zgemmujc=: (CONJ_NO_TRANSPOSE ,    CONJ_TRANSPOSE , UPPER)&zgemmtcore
zgemmucn=: (   CONJ_TRANSPOSE ,      NO_TRANSPOSE , UPPER)&zgemmtcore
zgemmuct=: (   CONJ_TRANSPOSE ,         TRANSPOSE , UPPER)&zgemmtcore
zgemmucj=: (   CONJ_TRANSPOSE , CONJ_NO_TRANSPOSE , UPPER)&zgemmtcore
zgemmucc=: (   CONJ_TRANSPOSE ,    CONJ_TRANSPOSE , UPPER)&zgemmtcore

NB. ---------------------------------------------------------
NB. Monad        API       Domain     A     Side    Reads in A    op1(A)     op2(B)
NB. symmllnn     object    mixed      SY    (1)     LT                 A          B
NB. symmllnt     object    mixed      SY    (1)     LT                 A          B^T
NB. symmllnj     object    mixed      SY    (1)     LT                 A     conj(B)
NB. symmllnc     object    mixed      SY    (1)     LT                 A          B^H
NB. symmlljn     object    mixed      SY    (1)     LT            conj(A)         B
NB. symmlljt     object    mixed      SY    (1)     LT            conj(A)         B^T
NB. symmlljj     object    mixed      SY    (1)     LT            conj(A)    conj(B)
NB. symmlljc     object    mixed      SY    (1)     LT            conj(A)         B^H
NB. symmlunn     object    mixed      SY    (1)     UT                 A          B
NB. symmlunt     object    mixed      SY    (1)     UT                 A          B^T
NB. symmlunj     object    mixed      SY    (1)     UT                 A     conj(B)
NB. symmlunc     object    mixed      SY    (1)     UT                 A          B^H
NB. symmlujn     object    mixed      SY    (1)     UT            conj(A)         B
NB. symmlujt     object    mixed      SY    (1)     UT            conj(A)         B^T
NB. symmlujj     object    mixed      SY    (1)     UT            conj(A)    conj(B)
NB. symmlujc     object    mixed      SY    (1)     UT            conj(A)         B^H
NB. symmrlnn     object    mixed      SY    (2)     LT                 A          B
NB. symmrlnt     object    mixed      SY    (2)     LT                 A          B^T
NB. symmrlnj     object    mixed      SY    (2)     LT                 A     conj(B)
NB. symmrlnc     object    mixed      SY    (2)     LT                 A          B^H
NB. symmrljn     object    mixed      SY    (2)     LT            conj(A)         B
NB. symmrljt     object    mixed      SY    (2)     LT            conj(A)         B^T
NB. symmrljj     object    mixed      SY    (2)     LT            conj(A)    conj(B)
NB. symmrljc     object    mixed      SY    (2)     LT            conj(A)         B^H
NB. symmrunn     object    mixed      SY    (2)     UT                 A          B
NB. symmrunt     object    mixed      SY    (2)     UT                 A          B^T
NB. symmrunj     object    mixed      SY    (2)     UT                 A     conj(B)
NB. symmrunc     object    mixed      SY    (2)     UT                 A          B^H
NB. symmrujn     object    mixed      SY    (2)     UT            conj(A)         B
NB. symmrujt     object    mixed      SY    (2)     UT            conj(A)         B^T
NB. symmrujj     object    mixed      SY    (2)     UT            conj(A)    conj(B)
NB. symmrujc     object    mixed      SY    (2)     UT            conj(A)         B^H
NB. dsymmllnn    typed     real       SY    (1)     LT                 A          B
NB. dsymmllnt    typed     real       SY    (1)     LT                 A          B^T
NB. dsymmlunn    typed     real       SY    (1)     UT                 A          B
NB. dsymmlunt    typed     real       SY    (1)     UT                 A          B^T
NB. dsymmrlnn    typed     real       SY    (1)     LT                 A          B
NB. dsymmrlnt    typed     real       SY    (1)     LT                 A          B^T
NB. dsymmrunn    typed     real       SY    (1)     UT                 A          B
NB. dsymmrunt    typed     real       SY    (1)     UT                 A          B^T
NB. zsymmllnn    typed     complex    SY    (1)     LT                 A          B
NB. zsymmllnt    typed     complex    SY    (1)     LT                 A          B^T
NB. zsymmllnj    typed     complex    SY    (1)     LT                 A     conj(B)
NB. zsymmllnc    typed     complex    SY    (1)     LT                 A          B^H
NB. zsymmlljn    typed     complex    SY    (1)     LT            conj(A)         B
NB. zsymmlljt    typed     complex    SY    (1)     LT            conj(A)         B^T
NB. zsymmlljj    typed     complex    SY    (1)     LT            conj(A)    conj(B)
NB. zsymmlljc    typed     complex    SY    (1)     LT            conj(A)         B^H
NB. zsymmlunn    typed     complex    SY    (1)     UT                 A          B
NB. zsymmlunt    typed     complex    SY    (1)     UT                 A          B^T
NB. zsymmlunj    typed     complex    SY    (1)     UT                 A     conj(B)
NB. zsymmlunc    typed     complex    SY    (1)     UT                 A          B^H
NB. zsymmlujn    typed     complex    SY    (1)     UT            conj(A)         B
NB. zsymmlujt    typed     complex    SY    (1)     UT            conj(A)         B^T
NB. zsymmlujj    typed     complex    SY    (1)     UT            conj(A)    conj(B)
NB. zsymmlujc    typed     complex    SY    (1)     UT            conj(A)         B^H
NB. zsymmrlnn    typed     complex    SY    (2)     LT                 A          B
NB. zsymmrlnt    typed     complex    SY    (2)     LT                 A          B^T
NB. zsymmrlnj    typed     complex    SY    (2)     LT                 A     conj(B)
NB. zsymmrlnc    typed     complex    SY    (2)     LT                 A          B^H
NB. zsymmrljn    typed     complex    SY    (2)     LT            conj(A)         B
NB. zsymmrljt    typed     complex    SY    (2)     LT            conj(A)         B^T
NB. zsymmrljj    typed     complex    SY    (2)     LT            conj(A)    conj(B)
NB. zsymmrljc    typed     complex    SY    (2)     LT            conj(A)         B^H
NB. zsymmrunn    typed     complex    SY    (2)     UT                 A          B
NB. zsymmrunt    typed     complex    SY    (2)     UT                 A          B^T
NB. zsymmrunj    typed     complex    SY    (2)     UT                 A     conj(B)
NB. zsymmrunc    typed     complex    SY    (2)     UT                 A          B^H
NB. zsymmrujn    typed     complex    SY    (2)     UT            conj(A)         B
NB. zsymmrujt    typed     complex    SY    (2)     UT            conj(A)         B^T
NB. zsymmrujj    typed     complex    SY    (2)     UT            conj(A)    conj(B)
NB. zsymmrujc    typed     complex    SY    (2)     UT            conj(A)         B^H
NB. hemmllnn     object    mixed      HE    (1)     LT                 A          B
NB. hemmllnt     object    mixed      HE    (1)     LT                 A          B^T
NB. hemmllnj     object    mixed      HE    (1)     LT                 A     conj(B)
NB. hemmllnc     object    mixed      HE    (1)     LT                 A          B^H
NB. hemmlljn     object    mixed      HE    (1)     LT            conj(A)         B
NB. hemmlljt     object    mixed      HE    (1)     LT            conj(A)         B^T
NB. hemmlljj     object    mixed      HE    (1)     LT            conj(A)    conj(B)
NB. hemmlljc     object    mixed      HE    (1)     LT            conj(A)         B^H
NB. hemmlunn     object    mixed      HE    (1)     UT                 A          B
NB. hemmlunt     object    mixed      HE    (1)     UT                 A          B^T
NB. hemmlunj     object    mixed      HE    (1)     UT                 A     conj(B)
NB. hemmlunc     object    mixed      HE    (1)     UT                 A          B^H
NB. hemmlujn     object    mixed      HE    (1)     UT            conj(A)         B
NB. hemmlujt     object    mixed      HE    (1)     UT            conj(A)         B^T
NB. hemmlujj     object    mixed      HE    (1)     UT            conj(A)    conj(B)
NB. hemmlujc     object    mixed      HE    (1)     UT            conj(A)         B^H
NB. hemmrlnn     object    mixed      HE    (2)     LT                 A          B
NB. hemmrlnt     object    mixed      HE    (2)     LT                 A          B^T
NB. hemmrlnj     object    mixed      HE    (2)     LT                 A     conj(B)
NB. hemmrlnc     object    mixed      HE    (2)     LT                 A          B^H
NB. hemmrljn     object    mixed      HE    (2)     LT            conj(A)         B
NB. hemmrljt     object    mixed      HE    (2)     LT            conj(A)         B^T
NB. hemmrljj     object    mixed      HE    (2)     LT            conj(A)    conj(B)
NB. hemmrljc     object    mixed      HE    (2)     LT            conj(A)         B^H
NB. hemmrunn     object    mixed      HE    (2)     UT                 A          B
NB. hemmrunt     object    mixed      HE    (2)     UT                 A          B^T
NB. hemmrunj     object    mixed      HE    (2)     UT                 A     conj(B)
NB. hemmrunc     object    mixed      HE    (2)     UT                 A          B^H
NB. hemmrujn     object    mixed      HE    (2)     UT            conj(A)         B
NB. hemmrujt     object    mixed      HE    (2)     UT            conj(A)         B^T
NB. hemmrujj     object    mixed      HE    (2)     UT            conj(A)    conj(B)
NB. hemmrujc     object    mixed      HE    (2)     UT            conj(A)         B^H
NB. zhemmllnn    typed     complex    HE    (1)     LT                 A          B
NB. zhemmllnt    typed     complex    HE    (1)     LT                 A          B^T
NB. zhemmllnj    typed     complex    HE    (1)     LT                 A     conj(B)
NB. zhemmllnc    typed     complex    HE    (1)     LT                 A          B^H
NB. zhemmlljn    typed     complex    HE    (1)     LT            conj(A)         B
NB. zhemmlljt    typed     complex    HE    (1)     LT            conj(A)         B^T
NB. zhemmlljj    typed     complex    HE    (1)     LT            conj(A)    conj(B)
NB. zhemmlljc    typed     complex    HE    (1)     LT            conj(A)         B^H
NB. zhemmlunn    typed     complex    HE    (1)     UT                 A          B
NB. zhemmlunt    typed     complex    HE    (1)     UT                 A          B^T
NB. zhemmlunj    typed     complex    HE    (1)     UT                 A     conj(B)
NB. zhemmlunc    typed     complex    HE    (1)     UT                 A          B^H
NB. zhemmlujn    typed     complex    HE    (1)     UT            conj(A)         B
NB. zhemmlujt    typed     complex    HE    (1)     UT            conj(A)         B^T
NB. zhemmlujj    typed     complex    HE    (1)     UT            conj(A)    conj(B)
NB. zhemmlujc    typed     complex    HE    (1)     UT            conj(A)         B^H
NB. zhemmrlnn    typed     complex    HE    (2)     LT                 A          B
NB. zhemmrlnt    typed     complex    HE    (2)     LT                 A          B^T
NB. zhemmrlnj    typed     complex    HE    (2)     LT                 A     conj(B)
NB. zhemmrlnc    typed     complex    HE    (2)     LT                 A          B^H
NB. zhemmrljn    typed     complex    HE    (2)     LT            conj(A)         B
NB. zhemmrljt    typed     complex    HE    (2)     LT            conj(A)         B^T
NB. zhemmrljj    typed     complex    HE    (2)     LT            conj(A)    conj(B)
NB. zhemmrljc    typed     complex    HE    (2)     LT            conj(A)         B^H
NB. zhemmrunn    typed     complex    HE    (2)     UT                 A          B
NB. zhemmrunt    typed     complex    HE    (2)     UT                 A          B^T
NB. zhemmrunj    typed     complex    HE    (2)     UT                 A     conj(B)
NB. zhemmrunc    typed     complex    HE    (2)     UT                 A          B^H
NB. zhemmrujn    typed     complex    HE    (2)     UT            conj(A)         B
NB. zhemmrujt    typed     complex    HE    (2)     UT            conj(A)         B^T
NB. zhemmrujj    typed     complex    HE    (2)     UT            conj(A)    conj(B)
NB. zhemmrujc    typed     complex    HE    (2)     UT            conj(A)         B^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C  (1)
NB.   or
NB.     C := alpha * op2(B) * op1(A) + beta * C  (2)
NB.   where A is Hermitian (symmetric), op1(A) is either A or
NB.   conj(A), and op2(B) is either B, B^T, conj(B) or B^H
NB.
NB. Syntax:
NB.   Cupd=. xxxmmxxxx alpha ; AA ; B ; beta ; C
NB. where
NB.   alpha - scalar
NB.   AA    - mn×mn-matrix, contains either LT or UT or both
NB.           part(s) of A
NB.   B     - m×n-matrix
NB.   beta  - scalar
NB.   C     - m×n-matrix
NB.   Cupd  - an updated C
NB.   A     - mn×mn-matrix, Hermitian (symmetric)
NB.   m     ≥ 0, the number of rows in B, C and Cupd, the
NB.           size of A and AA for xxxmmlxxx
NB.   n     ≥ 0, the number of columns in B, C and Cupd, the
NB.           size of A and AA for xxxmmrxxx
NB.   mn    = m for xxxmmlxxx or mn = n for xxxmmrxxx
NB.
NB. Notes:
NB. - monad       provides BLIS'
NB.   symmxxxx    bli_symm (...)
NB.   dsymmxxx    bli_dsymm(...)
NB.   zsymmxxx    bli_zsymm(...)
NB.   hemmxxxx    bli_hemm (...)
NB.   zhemmxxx    bli_zhemm(...)

symmllnn=:  (LEFT  , LOWER ,      NO_CONJUGATE ,      NO_TRANSPOSE)& symmcore
symmllnt=:  (LEFT  , LOWER ,      NO_CONJUGATE ,         TRANSPOSE)& symmcore
symmllnj=:  (LEFT  , LOWER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)& symmcore
symmllnc=:  (LEFT  , LOWER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)& symmcore
symmlljn=:  (LEFT  , LOWER ,         CONJUGATE ,      NO_TRANSPOSE)& symmcore
symmlljt=:  (LEFT  , LOWER ,         CONJUGATE ,         TRANSPOSE)& symmcore
symmlljj=:  (LEFT  , LOWER ,         CONJUGATE , CONJ_NO_TRANSPOSE)& symmcore
symmlljc=:  (LEFT  , LOWER ,         CONJUGATE ,    CONJ_TRANSPOSE)& symmcore
symmlunn=:  (LEFT  , UPPER ,      NO_CONJUGATE ,      NO_TRANSPOSE)& symmcore
symmlunt=:  (LEFT  , UPPER ,      NO_CONJUGATE ,         TRANSPOSE)& symmcore
symmlunj=:  (LEFT  , UPPER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)& symmcore
symmlunc=:  (LEFT  , UPPER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)& symmcore
symmlujn=:  (LEFT  , UPPER ,         CONJUGATE ,      NO_TRANSPOSE)& symmcore
symmlujt=:  (LEFT  , UPPER ,         CONJUGATE ,         TRANSPOSE)& symmcore
symmlujj=:  (LEFT  , UPPER ,         CONJUGATE , CONJ_NO_TRANSPOSE)& symmcore
symmlujc=:  (LEFT  , UPPER ,         CONJUGATE ,    CONJ_TRANSPOSE)& symmcore
symmrlnn=:  (RIGHT , LOWER ,      NO_CONJUGATE ,      NO_TRANSPOSE)& symmcore
symmrlnt=:  (RIGHT , LOWER ,      NO_CONJUGATE ,         TRANSPOSE)& symmcore
symmrlnj=:  (RIGHT , LOWER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)& symmcore
symmrlnc=:  (RIGHT , LOWER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)& symmcore
symmrljn=:  (RIGHT , LOWER ,         CONJUGATE ,      NO_TRANSPOSE)& symmcore
symmrljt=:  (RIGHT , LOWER ,         CONJUGATE ,         TRANSPOSE)& symmcore
symmrljj=:  (RIGHT , LOWER ,         CONJUGATE , CONJ_NO_TRANSPOSE)& symmcore
symmrljc=:  (RIGHT , LOWER ,         CONJUGATE ,    CONJ_TRANSPOSE)& symmcore
symmrunn=:  (RIGHT , UPPER ,      NO_CONJUGATE ,      NO_TRANSPOSE)& symmcore
symmrunt=:  (RIGHT , UPPER ,      NO_CONJUGATE ,         TRANSPOSE)& symmcore
symmrunj=:  (RIGHT , UPPER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)& symmcore
symmrunc=:  (RIGHT , UPPER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)& symmcore
symmrujn=:  (RIGHT , UPPER ,         CONJUGATE ,      NO_TRANSPOSE)& symmcore
symmrujt=:  (RIGHT , UPPER ,         CONJUGATE ,         TRANSPOSE)& symmcore
symmrujj=:  (RIGHT , UPPER ,         CONJUGATE , CONJ_NO_TRANSPOSE)& symmcore
symmrujc=:  (RIGHT , UPPER ,         CONJUGATE ,    CONJ_TRANSPOSE)& symmcore

dsymmllnn=: (LEFT  , LOWER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&dsymmcore
dsymmllnt=: (LEFT  , LOWER ,      NO_CONJUGATE ,         TRANSPOSE)&dsymmcore
dsymmlunn=: (LEFT  , UPPER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&dsymmcore
dsymmlunt=: (LEFT  , UPPER ,      NO_CONJUGATE ,         TRANSPOSE)&dsymmcore
dsymmrlnn=: (RIGHT , LOWER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&dsymmcore
dsymmrlnt=: (RIGHT , LOWER ,      NO_CONJUGATE ,         TRANSPOSE)&dsymmcore
dsymmrunn=: (RIGHT , UPPER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&dsymmcore
dsymmrunt=: (RIGHT , UPPER ,      NO_CONJUGATE ,         TRANSPOSE)&dsymmcore

zsymmllnn=: (LEFT  , LOWER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&zsymmcore
zsymmllnt=: (LEFT  , LOWER ,      NO_CONJUGATE ,         TRANSPOSE)&zsymmcore
zsymmllnj=: (LEFT  , LOWER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)&zsymmcore
zsymmllnc=: (LEFT  , LOWER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)&zsymmcore
zsymmlljn=: (LEFT  , LOWER ,         CONJUGATE ,      NO_TRANSPOSE)&zsymmcore
zsymmlljt=: (LEFT  , LOWER ,         CONJUGATE ,         TRANSPOSE)&zsymmcore
zsymmlljj=: (LEFT  , LOWER ,         CONJUGATE , CONJ_NO_TRANSPOSE)&zsymmcore
zsymmlljc=: (LEFT  , LOWER ,         CONJUGATE ,    CONJ_TRANSPOSE)&zsymmcore
zsymmlunn=: (LEFT  , UPPER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&zsymmcore
zsymmlunt=: (LEFT  , UPPER ,      NO_CONJUGATE ,         TRANSPOSE)&zsymmcore
zsymmlunj=: (LEFT  , UPPER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)&zsymmcore
zsymmlunc=: (LEFT  , UPPER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)&zsymmcore
zsymmlujn=: (LEFT  , UPPER ,         CONJUGATE ,      NO_TRANSPOSE)&zsymmcore
zsymmlujt=: (LEFT  , UPPER ,         CONJUGATE ,         TRANSPOSE)&zsymmcore
zsymmlujj=: (LEFT  , UPPER ,         CONJUGATE , CONJ_NO_TRANSPOSE)&zsymmcore
zsymmlujc=: (LEFT  , UPPER ,         CONJUGATE ,    CONJ_TRANSPOSE)&zsymmcore
zsymmrlnn=: (RIGHT , LOWER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&zsymmcore
zsymmrlnt=: (RIGHT , LOWER ,      NO_CONJUGATE ,         TRANSPOSE)&zsymmcore
zsymmrlnj=: (RIGHT , LOWER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)&zsymmcore
zsymmrlnc=: (RIGHT , LOWER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)&zsymmcore
zsymmrljn=: (RIGHT , LOWER ,         CONJUGATE ,      NO_TRANSPOSE)&zsymmcore
zsymmrljt=: (RIGHT , LOWER ,         CONJUGATE ,         TRANSPOSE)&zsymmcore
zsymmrljj=: (RIGHT , LOWER ,         CONJUGATE , CONJ_NO_TRANSPOSE)&zsymmcore
zsymmrljc=: (RIGHT , LOWER ,         CONJUGATE ,    CONJ_TRANSPOSE)&zsymmcore
zsymmrunn=: (RIGHT , UPPER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&zsymmcore
zsymmrunt=: (RIGHT , UPPER ,      NO_CONJUGATE ,         TRANSPOSE)&zsymmcore
zsymmrunj=: (RIGHT , UPPER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)&zsymmcore
zsymmrunc=: (RIGHT , UPPER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)&zsymmcore
zsymmrujn=: (RIGHT , UPPER ,         CONJUGATE ,      NO_TRANSPOSE)&zsymmcore
zsymmrujt=: (RIGHT , UPPER ,         CONJUGATE ,         TRANSPOSE)&zsymmcore
zsymmrujj=: (RIGHT , UPPER ,         CONJUGATE , CONJ_NO_TRANSPOSE)&zsymmcore
zsymmrujc=: (RIGHT , UPPER ,         CONJUGATE ,    CONJ_TRANSPOSE)&zsymmcore

hemmllnn=:  (LEFT  , LOWER ,      NO_CONJUGATE ,      NO_TRANSPOSE)& hemmcore
hemmllnt=:  (LEFT  , LOWER ,      NO_CONJUGATE ,         TRANSPOSE)& hemmcore
hemmllnj=:  (LEFT  , LOWER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)& hemmcore
hemmllnc=:  (LEFT  , LOWER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)& hemmcore
hemmlljn=:  (LEFT  , LOWER ,         CONJUGATE ,      NO_TRANSPOSE)& hemmcore
hemmlljt=:  (LEFT  , LOWER ,         CONJUGATE ,         TRANSPOSE)& hemmcore
hemmlljj=:  (LEFT  , LOWER ,         CONJUGATE , CONJ_NO_TRANSPOSE)& hemmcore
hemmlljc=:  (LEFT  , LOWER ,         CONJUGATE ,    CONJ_TRANSPOSE)& hemmcore
hemmlunn=:  (LEFT  , UPPER ,      NO_CONJUGATE ,      NO_TRANSPOSE)& hemmcore
hemmlunt=:  (LEFT  , UPPER ,      NO_CONJUGATE ,         TRANSPOSE)& hemmcore
hemmlunj=:  (LEFT  , UPPER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)& hemmcore
hemmlunc=:  (LEFT  , UPPER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)& hemmcore
hemmlujn=:  (LEFT  , UPPER ,         CONJUGATE ,      NO_TRANSPOSE)& hemmcore
hemmlujt=:  (LEFT  , UPPER ,         CONJUGATE ,         TRANSPOSE)& hemmcore
hemmlujj=:  (LEFT  , UPPER ,         CONJUGATE , CONJ_NO_TRANSPOSE)& hemmcore
hemmlujc=:  (LEFT  , UPPER ,         CONJUGATE ,    CONJ_TRANSPOSE)& hemmcore
hemmrlnn=:  (RIGHT , LOWER ,      NO_CONJUGATE ,      NO_TRANSPOSE)& hemmcore
hemmrlnt=:  (RIGHT , LOWER ,      NO_CONJUGATE ,         TRANSPOSE)& hemmcore
hemmrlnj=:  (RIGHT , LOWER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)& hemmcore
hemmrlnc=:  (RIGHT , LOWER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)& hemmcore
hemmrljn=:  (RIGHT , LOWER ,         CONJUGATE ,      NO_TRANSPOSE)& hemmcore
hemmrljt=:  (RIGHT , LOWER ,         CONJUGATE ,         TRANSPOSE)& hemmcore
hemmrljj=:  (RIGHT , LOWER ,         CONJUGATE , CONJ_NO_TRANSPOSE)& hemmcore
hemmrljc=:  (RIGHT , LOWER ,         CONJUGATE ,    CONJ_TRANSPOSE)& hemmcore
hemmrunn=:  (RIGHT , UPPER ,      NO_CONJUGATE ,      NO_TRANSPOSE)& hemmcore
hemmrunt=:  (RIGHT , UPPER ,      NO_CONJUGATE ,         TRANSPOSE)& hemmcore
hemmrunj=:  (RIGHT , UPPER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)& hemmcore
hemmrunc=:  (RIGHT , UPPER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)& hemmcore
hemmrujn=:  (RIGHT , UPPER ,         CONJUGATE ,      NO_TRANSPOSE)& hemmcore
hemmrujt=:  (RIGHT , UPPER ,         CONJUGATE ,         TRANSPOSE)& hemmcore
hemmrujj=:  (RIGHT , UPPER ,         CONJUGATE , CONJ_NO_TRANSPOSE)& hemmcore
hemmrujc=:  (RIGHT , UPPER ,         CONJUGATE ,    CONJ_TRANSPOSE)& hemmcore

zhemmllnn=: (LEFT  , LOWER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&zhemmcore
zhemmllnt=: (LEFT  , LOWER ,      NO_CONJUGATE ,         TRANSPOSE)&zhemmcore
zhemmllnj=: (LEFT  , LOWER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)&zhemmcore
zhemmllnc=: (LEFT  , LOWER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)&zhemmcore
zhemmlljn=: (LEFT  , LOWER ,         CONJUGATE ,      NO_TRANSPOSE)&zhemmcore
zhemmlljt=: (LEFT  , LOWER ,         CONJUGATE ,         TRANSPOSE)&zhemmcore
zhemmlljj=: (LEFT  , LOWER ,         CONJUGATE , CONJ_NO_TRANSPOSE)&zhemmcore
zhemmlljc=: (LEFT  , LOWER ,         CONJUGATE ,    CONJ_TRANSPOSE)&zhemmcore
zhemmlunn=: (LEFT  , UPPER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&zhemmcore
zhemmlunt=: (LEFT  , UPPER ,      NO_CONJUGATE ,         TRANSPOSE)&zhemmcore
zhemmlunj=: (LEFT  , UPPER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)&zhemmcore
zhemmlunc=: (LEFT  , UPPER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)&zhemmcore
zhemmlujn=: (LEFT  , UPPER ,         CONJUGATE ,      NO_TRANSPOSE)&zhemmcore
zhemmlujt=: (LEFT  , UPPER ,         CONJUGATE ,         TRANSPOSE)&zhemmcore
zhemmlujj=: (LEFT  , UPPER ,         CONJUGATE , CONJ_NO_TRANSPOSE)&zhemmcore
zhemmlujc=: (LEFT  , UPPER ,         CONJUGATE ,    CONJ_TRANSPOSE)&zhemmcore
zhemmrlnn=: (RIGHT , LOWER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&zhemmcore
zhemmrlnt=: (RIGHT , LOWER ,      NO_CONJUGATE ,         TRANSPOSE)&zhemmcore
zhemmrlnj=: (RIGHT , LOWER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)&zhemmcore
zhemmrlnc=: (RIGHT , LOWER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)&zhemmcore
zhemmrljn=: (RIGHT , LOWER ,         CONJUGATE ,      NO_TRANSPOSE)&zhemmcore
zhemmrljt=: (RIGHT , LOWER ,         CONJUGATE ,         TRANSPOSE)&zhemmcore
zhemmrljj=: (RIGHT , LOWER ,         CONJUGATE , CONJ_NO_TRANSPOSE)&zhemmcore
zhemmrljc=: (RIGHT , LOWER ,         CONJUGATE ,    CONJ_TRANSPOSE)&zhemmcore
zhemmrunn=: (RIGHT , UPPER ,      NO_CONJUGATE ,      NO_TRANSPOSE)&zhemmcore
zhemmrunt=: (RIGHT , UPPER ,      NO_CONJUGATE ,         TRANSPOSE)&zhemmcore
zhemmrunj=: (RIGHT , UPPER ,      NO_CONJUGATE , CONJ_NO_TRANSPOSE)&zhemmcore
zhemmrunc=: (RIGHT , UPPER ,      NO_CONJUGATE ,    CONJ_TRANSPOSE)&zhemmcore
zhemmrujn=: (RIGHT , UPPER ,         CONJUGATE ,      NO_TRANSPOSE)&zhemmcore
zhemmrujt=: (RIGHT , UPPER ,         CONJUGATE ,         TRANSPOSE)&zhemmcore
zhemmrujj=: (RIGHT , UPPER ,         CONJUGATE , CONJ_NO_TRANSPOSE)&zhemmcore
zhemmrujc=: (RIGHT , UPPER ,         CONJUGATE ,    CONJ_TRANSPOSE)&zhemmcore

NB. ---------------------------------------------------------
NB. Monad        API       Domain     Side    A     Reads in A    op(A)
NB. trmmllnn     object    mixed      (1)     L      LT                A
NB. trmmllnu     object    mixed      (1)     L1    SLT                A
NB. trmmlltn     object    mixed      (1)     L      LT                A^T
NB. trmmlltu     object    mixed      (1)     L1    SLT                A^T
NB. trmmlljn     object    mixed      (1)     L      LT           conj(A)
NB. trmmllju     object    mixed      (1)     L1    SLT           conj(A)
NB. trmmllcn     object    mixed      (1)     L      LT                A^H
NB. trmmllcu     object    mixed      (1)     L1    SLT                A^H
NB. trmmlunn     object    mixed      (1)     U      UT                A
NB. trmmlunu     object    mixed      (1)     U1    SUT                A
NB. trmmlutn     object    mixed      (1)     U      UT                A^T
NB. trmmlutu     object    mixed      (1)     U1    SUT                A^T
NB. trmmlujn     object    mixed      (1)     U      UT           conj(A)
NB. trmmluju     object    mixed      (1)     U1    SUT           conj(A)
NB. trmmlucn     object    mixed      (1)     U      UT                A^H
NB. trmmlucu     object    mixed      (1)     U1    SUT                A^H
NB. trmmrlnn     object    mixed      (2)     L     LT                 A
NB. trmmrlnu     object    mixed      (2)     L1    SLT                A
NB. trmmrltn     object    mixed      (2)     L      LT                A^T
NB. trmmrltu     object    mixed      (2)     L1    SLT                A^T
NB. trmmrljn     object    mixed      (2)     L      LT           conj(A)
NB. trmmrlju     object    mixed      (2)     L1    SLT           conj(A)
NB. trmmrlcn     object    mixed      (2)     L      LT                A^H
NB. trmmrlcu     object    mixed      (2)     L1    SLT                A^H
NB. trmmrunn     object    mixed      (2)     U     UT                 A
NB. trmmrunu     object    mixed      (2)     U1    SUT                A
NB. trmmrutn     object    mixed      (2)     U      UT                A^T
NB. trmmrutu     object    mixed      (2)     U1    SUT                A^T
NB. trmmrujn     object    mixed      (2)     U      UT           conj(A)
NB. trmmruju     object    mixed      (2)     U1    SUT           conj(A)
NB. trmmrucn     object    mixed      (2)     U      UT                A^H
NB. trmmrucu     object    mixed      (2)     U1    SUT                A^H
NB. dtrmmllnn    typed     real       (1)     L      LT                A
NB. dtrmmllnu    typed     real       (1)     L1    SLT                A
NB. dtrmmlltn    typed     real       (1)     L      LT                A^T
NB. dtrmmlltu    typed     real       (1)     L1    SLT                A^T
NB. dtrmmlunn    typed     real       (1)     U      UT                A
NB. dtrmmlunu    typed     real       (1)     U1    SUT                A
NB. dtrmmlutn    typed     real       (1)     U      UT                A^T
NB. dtrmmlutu    typed     real       (1)     U1    SUT                A^T
NB. dtrmmrlnn    typed     real       (2)     L     LT                 A
NB. dtrmmrlnu    typed     real       (2)     L1    SLT                A
NB. dtrmmrltn    typed     real       (2)     L      LT                A^T
NB. dtrmmrltu    typed     real       (2)     L1    SLT                A^T
NB. dtrmmrunn    typed     real       (2)     U     UT                 A
NB. dtrmmrunu    typed     real       (2)     U1    SUT                A
NB. dtrmmrutn    typed     real       (2)     U      UT                A^T
NB. dtrmmrutu    typed     real       (2)     U1    SUT                A^T
NB. ztrmmllnn    typed     complex    (1)     L      LT                A
NB. ztrmmllnu    typed     complex    (1)     L1    SLT                A
NB. ztrmmlltn    typed     complex    (1)     L      LT                A^T
NB. ztrmmlltu    typed     complex    (1)     L1    SLT                A^T
NB. ztrmmlljn    typed     complex    (1)     L      LT           conj(A)
NB. ztrmmllju    typed     complex    (1)     L1    SLT           conj(A)
NB. ztrmmllcn    typed     complex    (1)     L      LT                A^H
NB. ztrmmllcu    typed     complex    (1)     L1    SLT                A^H
NB. ztrmmlunn    typed     complex    (1)     U      UT                A
NB. ztrmmlunu    typed     complex    (1)     U1    SUT                A
NB. ztrmmlutn    typed     complex    (1)     U      UT                A^T
NB. ztrmmlutu    typed     complex    (1)     U1    SUT                A^T
NB. ztrmmlujn    typed     complex    (1)     U      UT           conj(A)
NB. ztrmmluju    typed     complex    (1)     U1    SUT           conj(A)
NB. ztrmmlucn    typed     complex    (1)     U      UT                A^H
NB. ztrmmlucu    typed     complex    (1)     U1    SUT                A^H
NB. ztrmmrlnn    typed     complex    (2)     L     LT                 A
NB. ztrmmrlnu    typed     complex    (2)     L1    SLT                A
NB. ztrmmrltn    typed     complex    (2)     L      LT                A^T
NB. ztrmmrltu    typed     complex    (2)     L1    SLT                A^T
NB. ztrmmrljn    typed     complex    (2)     L      LT           conj(A)
NB. ztrmmrlju    typed     complex    (2)     L1    SLT           conj(A)
NB. ztrmmrlcn    typed     complex    (2)     L      LT                A^H
NB. ztrmmrlcu    typed     complex    (2)     L1    SLT                A^H
NB. ztrmmrunn    typed     complex    (2)     U     UT                 A
NB. ztrmmrunu    typed     complex    (2)     U1    SUT                A
NB. ztrmmrutn    typed     complex    (2)     U      UT                A^T
NB. ztrmmrutu    typed     complex    (2)     U1    SUT                A^T
NB. ztrmmrujn    typed     complex    (2)     U      UT           conj(A)
NB. ztrmmruju    typed     complex    (2)     U1    SUT           conj(A)
NB. ztrmmrucn    typed     complex    (2)     U      UT                A^H
NB. ztrmmrucu    typed     complex    (2)     U1    SUT                A^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     B := alpha * op(A) * B  (1)
NB.   or
NB.     B := alpha * B * op(A)  (2)
NB.   where A is triangular, and op(A) is either A, A^T,
NB.   conj(A) or A^H
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
NB. - monad        provides BLIS'
NB.   trmmxxxx     bli_trmm (...)
NB.   dtrmmxxxx    bli_dtrmm(...)
NB.   ztrmmxxxx    bli_ztrmm(...)

trmmllnn=:  (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmllnu=:  (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmlltn=:  (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmlltu=:  (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmlljn=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmllju=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmllcn=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmllcu=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmlunn=:  (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmlunu=:  (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmlutn=:  (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmlutu=:  (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmlujn=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmluju=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmlucn=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmlucu=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmrlnn=:  (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmrlnu=:  (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmrltn=:  (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmrltu=:  (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmrljn=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmrlju=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmrlcn=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmrlcu=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmrunn=:  (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmrunu=:  (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmrutn=:  (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmrutu=:  (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmrujn=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmruju=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG)& trmmcore
trmmrucn=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG)& trmmcore
trmmrucu=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG)& trmmcore

dtrmmllnn=: (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG)&dtrmmcore
dtrmmllnu=: (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG)&dtrmmcore
dtrmmlltn=: (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG)&dtrmmcore
dtrmmlltu=: (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG)&dtrmmcore
dtrmmlunn=: (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG)&dtrmmcore
dtrmmlunu=: (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG)&dtrmmcore
dtrmmlutn=: (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG)&dtrmmcore
dtrmmlutu=: (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG)&dtrmmcore
dtrmmrlnn=: (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG)&dtrmmcore
dtrmmrlnu=: (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG)&dtrmmcore
dtrmmrltn=: (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG)&dtrmmcore
dtrmmrltu=: (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG)&dtrmmcore
dtrmmrunn=: (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG)&dtrmmcore
dtrmmrunu=: (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG)&dtrmmcore
dtrmmrutn=: (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG)&dtrmmcore
dtrmmrutu=: (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG)&dtrmmcore

ztrmmllnn=: (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmllnu=: (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmlltn=: (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmlltu=: (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmlljn=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmllju=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmllcn=: (LEFT  , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmllcu=: (LEFT  , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmlunn=: (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmlunu=: (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmlutn=: (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmlutu=: (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmlujn=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmluju=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmlucn=: (LEFT  , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmlucu=: (LEFT  , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmrlnn=: (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmrlnu=: (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmrltn=: (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmrltu=: (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmrljn=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmrlju=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmrlcn=: (RIGHT , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmrlcu=: (RIGHT , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmrunn=: (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmrunu=: (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmrutn=: (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmrutu=: (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmrujn=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmruju=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore
ztrmmrucn=: (RIGHT , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG)&ztrmmcore
ztrmmrucu=: (RIGHT , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG)&ztrmmcore

NB. ---------------------------------------------------------
NB. Monad         API       Domain     Side    A     Reads in A    op1(A)      op2(B)
NB. trmmllnnn     object    mixed      (1)     L      LT                A           B
NB. trmmllnnt     object    mixed      (1)     L      LT                A           B^T
NB. trmmllnnj     object    mixed      (1)     L      LT                A      conj(B)
NB. trmmllnnc     object    mixed      (1)     L      LT                A           B^H
NB. trmmllnun     object    mixed      (1)     L1    SLT                A           B
NB. trmmllnut     object    mixed      (1)     L1    SLT                A           B^T
NB. trmmllnuj     object    mixed      (1)     L1    SLT                A      conj(B)
NB. trmmllnuc     object    mixed      (1)     L1    SLT                A           B^H
NB. trmmlltnn     object    mixed      (1)     L      LT                A^T         B
NB. trmmlltnt     object    mixed      (1)     L      LT                A^T         B^T
NB. trmmlltnj     object    mixed      (1)     L      LT                A^T    conj(B)
NB. trmmlltnc     object    mixed      (1)     L      LT                A^T         B^H
NB. trmmlltun     object    mixed      (1)     L1    SLT                A^T         B
NB. trmmlltut     object    mixed      (1)     L1    SLT                A^T         B^T
NB. trmmlltuj     object    mixed      (1)     L1    SLT                A^T    conj(B)
NB. trmmlltuc     object    mixed      (1)     L1    SLT                A^T         B^H
NB. trmmlljnn     object    mixed      (1)     L      LT           conj(A)          B
NB. trmmlljnt     object    mixed      (1)     L      LT           conj(A)          B^T
NB. trmmlljnj     object    mixed      (1)     L      LT           conj(A)     conj(B)
NB. trmmlljnc     object    mixed      (1)     L      LT           conj(A)          B^H
NB. trmmlljun     object    mixed      (1)     L1    SLT           conj(A)          B
NB. trmmlljut     object    mixed      (1)     L1    SLT           conj(A)          B^T
NB. trmmlljuj     object    mixed      (1)     L1    SLT           conj(A)     conj(B)
NB. trmmlljuc     object    mixed      (1)     L1    SLT           conj(A)          B^H
NB. trmmllcnn     object    mixed      (1)     L      LT                A^H         B
NB. trmmllcnt     object    mixed      (1)     L      LT                A^H         B^T
NB. trmmllcnj     object    mixed      (1)     L      LT                A^H    conj(B)
NB. trmmllcnc     object    mixed      (1)     L      LT                A^H         B^H
NB. trmmllcun     object    mixed      (1)     L1    SLT                A^H         B
NB. trmmllcut     object    mixed      (1)     L1    SLT                A^H         B^T
NB. trmmllcuj     object    mixed      (1)     L1    SLT                A^H    conj(B)
NB. trmmllcuc     object    mixed      (1)     L1    SLT                A^H         B^H
NB. trmmlunnn     object    mixed      (1)     U      UT                A           B
NB. trmmlunnt     object    mixed      (1)     U      UT                A           B^T
NB. trmmlunnj     object    mixed      (1)     U      UT                A      conj(B)
NB. trmmlunnc     object    mixed      (1)     U      UT                A           B^H
NB. trmmlunun     object    mixed      (1)     U1    SUT                A           B
NB. trmmlunut     object    mixed      (1)     U1    SUT                A           B^T
NB. trmmlunuj     object    mixed      (1)     U1    SUT                A      conj(B)
NB. trmmlunuc     object    mixed      (1)     U1    SUT                A           B^H
NB. trmmlutnn     object    mixed      (1)     U      UT                A^T         B
NB. trmmlutnt     object    mixed      (1)     U      UT                A^T         B^T
NB. trmmlutnj     object    mixed      (1)     U      UT                A^T    conj(B)
NB. trmmlutnc     object    mixed      (1)     U      UT                A^T         B^H
NB. trmmlutun     object    mixed      (1)     U1    SUT                A^T         B
NB. trmmlutut     object    mixed      (1)     U1    SUT                A^T         B^T
NB. trmmlutuj     object    mixed      (1)     U1    SUT                A^T    conj(B)
NB. trmmlutuc     object    mixed      (1)     U1    SUT                A^T         B^H
NB. trmmlujnn     object    mixed      (1)     U      UT           conj(A)          B
NB. trmmlujnt     object    mixed      (1)     U      UT           conj(A)          B^T
NB. trmmlujnj     object    mixed      (1)     U      UT           conj(A)     conj(B)
NB. trmmlujnc     object    mixed      (1)     U      UT           conj(A)          B^H
NB. trmmlujun     object    mixed      (1)     U1    SUT           conj(A)          B
NB. trmmlujut     object    mixed      (1)     U1    SUT           conj(A)          B^T
NB. trmmlujuj     object    mixed      (1)     U1    SUT           conj(A)     conj(B)
NB. trmmlujuc     object    mixed      (1)     U1    SUT           conj(A)          B^H
NB. trmmlucnn     object    mixed      (1)     U      UT                A^H         B
NB. trmmlucnt     object    mixed      (1)     U      UT                A^H         B^T
NB. trmmlucnj     object    mixed      (1)     U      UT                A^H    conj(B)
NB. trmmlucnc     object    mixed      (1)     U      UT                A^H         B^H
NB. trmmlucun     object    mixed      (1)     U1    SUT                A^H         B
NB. trmmlucut     object    mixed      (1)     U1    SUT                A^H         B^T
NB. trmmlucuj     object    mixed      (1)     U1    SUT                A^H    conj(B)
NB. trmmlucuc     object    mixed      (1)     U1    SUT                A^H         B^H
NB. trmmrlnnn     object    mixed      (2)     L      LT                A           B
NB. trmmrlnnt     object    mixed      (2)     L      LT                A           B^T
NB. trmmrlnnj     object    mixed      (2)     L      LT                A      conj(B)
NB. trmmrlnnc     object    mixed      (2)     L      LT                A           B^H
NB. trmmrlnun     object    mixed      (2)     L1    SLT                A           B
NB. trmmrlnut     object    mixed      (2)     L1    SLT                A           B^T
NB. trmmrlnuj     object    mixed      (2)     L1    SLT                A      conj(B)
NB. trmmrlnuc     object    mixed      (2)     L1    SLT                A           B^H
NB. trmmrltnn     object    mixed      (2)     L      LT                A^T         B
NB. trmmrltnt     object    mixed      (2)     L      LT                A^T         B^T
NB. trmmrltnj     object    mixed      (2)     L      LT                A^T    conj(B)
NB. trmmrltnc     object    mixed      (2)     L      LT                A^T         B^H
NB. trmmrltun     object    mixed      (2)     L1    SLT                A^T         B
NB. trmmrltut     object    mixed      (2)     L1    SLT                A^T         B^T
NB. trmmrltuj     object    mixed      (2)     L1    SLT                A^T    conj(B)
NB. trmmrltuc     object    mixed      (2)     L1    SLT                A^T         B^H
NB. trmmrljnn     object    mixed      (2)     L      LT           conj(A)          B
NB. trmmrljnt     object    mixed      (2)     L      LT           conj(A)          B^T
NB. trmmrljnj     object    mixed      (2)     L      LT           conj(A)     conj(B)
NB. trmmrljnc     object    mixed      (2)     L      LT           conj(A)          B^H
NB. trmmrljun     object    mixed      (2)     L1    SLT           conj(A)          B
NB. trmmrljut     object    mixed      (2)     L1    SLT           conj(A)          B^T
NB. trmmrljuj     object    mixed      (2)     L1    SLT           conj(A)     conj(B)
NB. trmmrljuc     object    mixed      (2)     L1    SLT           conj(A)          B^H
NB. trmmrlcnn     object    mixed      (2)     L      LT                A^H         B
NB. trmmrlcnt     object    mixed      (2)     L      LT                A^H         B^T
NB. trmmrlcnj     object    mixed      (2)     L      LT                A^H    conj(B)
NB. trmmrlcnc     object    mixed      (2)     L      LT                A^H         B^H
NB. trmmrlcun     object    mixed      (2)     L1    SLT                A^H         B
NB. trmmrlcut     object    mixed      (2)     L1    SLT                A^H         B^T
NB. trmmrlcuj     object    mixed      (2)     L1    SLT                A^H    conj(B)
NB. trmmrlcuc     object    mixed      (2)     L1    SLT                A^H         B^H
NB. trmmrunnn     object    mixed      (2)     U     UT                 A           B
NB. trmmrunnt     object    mixed      (2)     U     UT                 A           B^T
NB. trmmrunnj     object    mixed      (2)     U     UT                 A      conj(B)
NB. trmmrunnc     object    mixed      (2)     U     UT                 A           B^H
NB. trmmrunun     object    mixed      (2)     U1    SUT                A           B
NB. trmmrunut     object    mixed      (2)     U1    SUT                A           B^T
NB. trmmrunuj     object    mixed      (2)     U1    SUT                A      conj(B)
NB. trmmrunuc     object    mixed      (2)     U1    SUT                A           B^H
NB. trmmrutnn     object    mixed      (2)     U      UT                A^T         B
NB. trmmrutnt     object    mixed      (2)     U      UT                A^T         B^T
NB. trmmrutnj     object    mixed      (2)     U      UT                A^T    conj(B)
NB. trmmrutnc     object    mixed      (2)     U      UT                A^T         B^H
NB. trmmrutun     object    mixed      (2)     U1    SUT                A^T         B
NB. trmmrutut     object    mixed      (2)     U1    SUT                A^T         B^T
NB. trmmrutuj     object    mixed      (2)     U1    SUT                A^T    conj(B)
NB. trmmrutuc     object    mixed      (2)     U1    SUT                A^T         B^H
NB. trmmrujnn     object    mixed      (2)     U      UT           conj(A)          B
NB. trmmrujnt     object    mixed      (2)     U      UT           conj(A)          B^T
NB. trmmrujnj     object    mixed      (2)     U      UT           conj(A)     conj(B)
NB. trmmrujnc     object    mixed      (2)     U      UT           conj(A)          B^H
NB. trmmrujun     object    mixed      (2)     U1    SUT           conj(A)          B
NB. trmmrujut     object    mixed      (2)     U1    SUT           conj(A)          B^T
NB. trmmrujuj     object    mixed      (2)     U1    SUT           conj(A)     conj(B)
NB. trmmrujuc     object    mixed      (2)     U1    SUT           conj(A)          B^H
NB. trmmrucnn     object    mixed      (2)     U      UT                A^H         B
NB. trmmrucnt     object    mixed      (2)     U      UT                A^H         B^T
NB. trmmrucnj     object    mixed      (2)     U      UT                A^H    conj(B)
NB. trmmrucnc     object    mixed      (2)     U      UT                A^H         B^H
NB. trmmrucun     object    mixed      (2)     U1    SUT                A^H         B
NB. trmmrucut     object    mixed      (2)     U1    SUT                A^H         B^T
NB. trmmrucuj     object    mixed      (2)     U1    SUT                A^H    conj(B)
NB. trmmrucuc     object    mixed      (2)     U1    SUT                A^H         B^H
NB. dtrmmllnnn    typed     real       (1)     L      LT                A           B
NB. dtrmmllnnt    typed     real       (1)     L      LT                A           B^T
NB. dtrmmllnun    typed     real       (1)     L1    SLT                A           B
NB. dtrmmllnut    typed     real       (1)     L1    SLT                A           B^T
NB. dtrmmlltnn    typed     real       (1)     L      LT                A^T         B
NB. dtrmmlltnt    typed     real       (1)     L      LT                A^T         B^T
NB. dtrmmlltun    typed     real       (1)     L1    SLT                A^T         B
NB. dtrmmlltut    typed     real       (1)     L1    SLT                A^T         B^T
NB. dtrmmlunnn    typed     real       (1)     U      UT                A           B
NB. dtrmmlunnt    typed     real       (1)     U      UT                A           B^T
NB. dtrmmlunun    typed     real       (1)     U1    SUT                A           B
NB. dtrmmlunut    typed     real       (1)     U1    SUT                A           B^T
NB. dtrmmlutnn    typed     real       (1)     U      UT                A^T         B
NB. dtrmmlutnt    typed     real       (1)     U      UT                A^T         B^T
NB. dtrmmlutun    typed     real       (1)     U1    SUT                A^T         B
NB. dtrmmlutut    typed     real       (1)     U1    SUT                A^T         B^T
NB. dtrmmrlnnn    typed     real       (2)     L     LT                 A           B
NB. dtrmmrlnnt    typed     real       (2)     L     LT                 A           B^T
NB. dtrmmrlnun    typed     real       (2)     L1    SLT                A           B
NB. dtrmmrlnut    typed     real       (2)     L1    SLT                A           B^T
NB. dtrmmrltnn    typed     real       (2)     L      LT                A^T         B
NB. dtrmmrltnt    typed     real       (2)     L      LT                A^T         B^T
NB. dtrmmrltun    typed     real       (2)     L1    SLT                A^T         B
NB. dtrmmrltut    typed     real       (2)     L1    SLT                A^T         B^T
NB. dtrmmrunnn    typed     real       (2)     U     UT                 A           B
NB. dtrmmrunnt    typed     real       (2)     U     UT                 A           B^T
NB. dtrmmrunun    typed     real       (2)     U1    SUT                A           B
NB. dtrmmrunut    typed     real       (2)     U1    SUT                A           B^T
NB. dtrmmrutnn    typed     real       (2)     U      UT                A^T         B
NB. dtrmmrutnt    typed     real       (2)     U      UT                A^T         B^T
NB. dtrmmrutun    typed     real       (2)     U1    SUT                A^T         B
NB. dtrmmrutut    typed     real       (2)     U1    SUT                A^T         B^T
NB. ztrmmllnnn    typed     complex    (1)     L      LT                A           B
NB. ztrmmllnnt    typed     complex    (1)     L      LT                A           B^T
NB. ztrmmllnnj    typed     complex    (1)     L      LT                A      conj(B)
NB. ztrmmllnnc    typed     complex    (1)     L      LT                A           B^H
NB. ztrmmllnun    typed     complex    (1)     L1    SLT                A           B
NB. ztrmmllnut    typed     complex    (1)     L1    SLT                A           B^T
NB. ztrmmllnuj    typed     complex    (1)     L1    SLT                A      conj(B)
NB. ztrmmllnuc    typed     complex    (1)     L1    SLT                A           B^H
NB. ztrmmlltnn    typed     complex    (1)     L      LT                A^T         B
NB. ztrmmlltnt    typed     complex    (1)     L      LT                A^T         B^T
NB. ztrmmlltnj    typed     complex    (1)     L      LT                A^T    conj(B)
NB. ztrmmlltnc    typed     complex    (1)     L      LT                A^T         B^H
NB. ztrmmlltun    typed     complex    (1)     L1    SLT                A^T         B
NB. ztrmmlltut    typed     complex    (1)     L1    SLT                A^T         B^T
NB. ztrmmlltuj    typed     complex    (1)     L1    SLT                A^T    conj(B)
NB. ztrmmlltuc    typed     complex    (1)     L1    SLT                A^T         B^H
NB. ztrmmlljnn    typed     complex    (1)     L      LT           conj(A)          B
NB. ztrmmlljnt    typed     complex    (1)     L      LT           conj(A)          B^T
NB. ztrmmlljnj    typed     complex    (1)     L      LT           conj(A)     conj(B)
NB. ztrmmlljnc    typed     complex    (1)     L      LT           conj(A)          B^H
NB. ztrmmlljun    typed     complex    (1)     L1    SLT           conj(A)          B
NB. ztrmmlljut    typed     complex    (1)     L1    SLT           conj(A)          B^T
NB. ztrmmlljuj    typed     complex    (1)     L1    SLT           conj(A)     conj(B)
NB. ztrmmlljuc    typed     complex    (1)     L1    SLT           conj(A)          B^H
NB. ztrmmllcnn    typed     complex    (1)     L      LT                A^H         B
NB. ztrmmllcnt    typed     complex    (1)     L      LT                A^H         B^T
NB. ztrmmllcnj    typed     complex    (1)     L      LT                A^H    conj(B)
NB. ztrmmllcnc    typed     complex    (1)     L      LT                A^H         B^H
NB. ztrmmllcun    typed     complex    (1)     L1    SLT                A^H         B
NB. ztrmmllcut    typed     complex    (1)     L1    SLT                A^H         B^T
NB. ztrmmllcuj    typed     complex    (1)     L1    SLT                A^H    conj(B)
NB. ztrmmllcuc    typed     complex    (1)     L1    SLT                A^H         B^H
NB. ztrmmlunnn    typed     complex    (1)     U      UT                A           B
NB. ztrmmlunnt    typed     complex    (1)     U      UT                A           B^T
NB. ztrmmlunnj    typed     complex    (1)     U      UT                A      conj(B)
NB. ztrmmlunnc    typed     complex    (1)     U      UT                A           B^H
NB. ztrmmlunun    typed     complex    (1)     U1    SUT                A           B
NB. ztrmmlunut    typed     complex    (1)     U1    SUT                A           B^T
NB. ztrmmlunuj    typed     complex    (1)     U1    SUT                A      conj(B)
NB. ztrmmlunuc    typed     complex    (1)     U1    SUT                A           B^H
NB. ztrmmlutnn    typed     complex    (1)     U      UT                A^T         B
NB. ztrmmlutnt    typed     complex    (1)     U      UT                A^T         B^T
NB. ztrmmlutnj    typed     complex    (1)     U      UT                A^T    conj(B)
NB. ztrmmlutnc    typed     complex    (1)     U      UT                A^T         B^H
NB. ztrmmlutun    typed     complex    (1)     U1    SUT                A^T         B
NB. ztrmmlutut    typed     complex    (1)     U1    SUT                A^T         B^T
NB. ztrmmlutuj    typed     complex    (1)     U1    SUT                A^T    conj(B)
NB. ztrmmlutuc    typed     complex    (1)     U1    SUT                A^T         B^H
NB. ztrmmlujnn    typed     complex    (1)     U      UT           conj(A)          B
NB. ztrmmlujnt    typed     complex    (1)     U      UT           conj(A)          B^T
NB. ztrmmlujnj    typed     complex    (1)     U      UT           conj(A)     conj(B)
NB. ztrmmlujnc    typed     complex    (1)     U      UT           conj(A)          B^H
NB. ztrmmlujun    typed     complex    (1)     U1    SUT           conj(A)          B
NB. ztrmmlujut    typed     complex    (1)     U1    SUT           conj(A)          B^T
NB. ztrmmlujuj    typed     complex    (1)     U1    SUT           conj(A)     conj(B)
NB. ztrmmlujuc    typed     complex    (1)     U1    SUT           conj(A)          B^H
NB. ztrmmlucnn    typed     complex    (1)     U      UT                A^H         B
NB. ztrmmlucnt    typed     complex    (1)     U      UT                A^H         B^T
NB. ztrmmlucnj    typed     complex    (1)     U      UT                A^H    conj(B)
NB. ztrmmlucnc    typed     complex    (1)     U      UT                A^H         B^H
NB. ztrmmlucun    typed     complex    (1)     U1    SUT                A^H         B
NB. ztrmmlucut    typed     complex    (1)     U1    SUT                A^H         B^T
NB. ztrmmlucuj    typed     complex    (1)     U1    SUT                A^H    conj(B)
NB. ztrmmlucuc    typed     complex    (1)     U1    SUT                A^H         B^H
NB. ztrmmrlnnn    typed     complex    (2)     L      LT                A           B
NB. ztrmmrlnnt    typed     complex    (2)     L      LT                A           B^T
NB. ztrmmrlnnj    typed     complex    (2)     L      LT                A      conj(B)
NB. ztrmmrlnnc    typed     complex    (2)     L      LT                A           B^H
NB. ztrmmrlnun    typed     complex    (2)     L1    SLT                A           B
NB. ztrmmrlnut    typed     complex    (2)     L1    SLT                A           B^T
NB. ztrmmrlnuj    typed     complex    (2)     L1    SLT                A      conj(B)
NB. ztrmmrlnuc    typed     complex    (2)     L1    SLT                A           B^H
NB. ztrmmrltnn    typed     complex    (2)     L      LT                A^T         B
NB. ztrmmrltnt    typed     complex    (2)     L      LT                A^T         B^T
NB. ztrmmrltnj    typed     complex    (2)     L      LT                A^T    conj(B)
NB. ztrmmrltnc    typed     complex    (2)     L      LT                A^T         B^H
NB. ztrmmrltun    typed     complex    (2)     L1    SLT                A^T         B
NB. ztrmmrltut    typed     complex    (2)     L1    SLT                A^T         B^T
NB. ztrmmrltuj    typed     complex    (2)     L1    SLT                A^T    conj(B)
NB. ztrmmrltuc    typed     complex    (2)     L1    SLT                A^T         B^H
NB. ztrmmrljnn    typed     complex    (2)     L      LT           conj(A)          B
NB. ztrmmrljnt    typed     complex    (2)     L      LT           conj(A)          B^T
NB. ztrmmrljnj    typed     complex    (2)     L      LT           conj(A)     conj(B)
NB. ztrmmrljnc    typed     complex    (2)     L      LT           conj(A)          B^H
NB. ztrmmrljun    typed     complex    (2)     L1    SLT           conj(A)          B
NB. ztrmmrljut    typed     complex    (2)     L1    SLT           conj(A)          B^T
NB. ztrmmrljuj    typed     complex    (2)     L1    SLT           conj(A)     conj(B)
NB. ztrmmrljuc    typed     complex    (2)     L1    SLT           conj(A)          B^H
NB. ztrmmrlcnn    typed     complex    (2)     L      LT                A^H         B
NB. ztrmmrlcnt    typed     complex    (2)     L      LT                A^H         B^T
NB. ztrmmrlcnj    typed     complex    (2)     L      LT                A^H    conj(B)
NB. ztrmmrlcnc    typed     complex    (2)     L      LT                A^H         B^H
NB. ztrmmrlcun    typed     complex    (2)     L1    SLT                A^H         B
NB. ztrmmrlcut    typed     complex    (2)     L1    SLT                A^H         B^T
NB. ztrmmrlcuj    typed     complex    (2)     L1    SLT                A^H    conj(B)
NB. ztrmmrlcuc    typed     complex    (2)     L1    SLT                A^H         B^H
NB. ztrmmrunnn    typed     complex    (2)     U     UT                 A           B
NB. ztrmmrunnt    typed     complex    (2)     U     UT                 A           B^T
NB. ztrmmrunnj    typed     complex    (2)     U     UT                 A      conj(B)
NB. ztrmmrunnc    typed     complex    (2)     U     UT                 A           B^H
NB. ztrmmrunun    typed     complex    (2)     U1    SUT                A           B
NB. ztrmmrunut    typed     complex    (2)     U1    SUT                A           B^T
NB. ztrmmrunuj    typed     complex    (2)     U1    SUT                A      conj(B)
NB. ztrmmrunuc    typed     complex    (2)     U1    SUT                A           B^H
NB. ztrmmrutnn    typed     complex    (2)     U      UT                A^T         B
NB. ztrmmrutnt    typed     complex    (2)     U      UT                A^T         B^T
NB. ztrmmrutnj    typed     complex    (2)     U      UT                A^T    conj(B)
NB. ztrmmrutnc    typed     complex    (2)     U      UT                A^T         B^H
NB. ztrmmrutun    typed     complex    (2)     U1    SUT                A^T         B
NB. ztrmmrutut    typed     complex    (2)     U1    SUT                A^T         B^T
NB. ztrmmrutuj    typed     complex    (2)     U1    SUT                A^T    conj(B)
NB. ztrmmrutuc    typed     complex    (2)     U1    SUT                A^T         B^H
NB. ztrmmrujnn    typed     complex    (2)     U      UT           conj(A)          B
NB. ztrmmrujnt    typed     complex    (2)     U      UT           conj(A)          B^T
NB. ztrmmrujnj    typed     complex    (2)     U      UT           conj(A)     conj(B)
NB. ztrmmrujnc    typed     complex    (2)     U      UT           conj(A)          B^H
NB. ztrmmrujun    typed     complex    (2)     U1    SUT           conj(A)          B
NB. ztrmmrujut    typed     complex    (2)     U1    SUT           conj(A)          B^T
NB. ztrmmrujuj    typed     complex    (2)     U1    SUT           conj(A)     conj(B)
NB. ztrmmrujuc    typed     complex    (2)     U1    SUT           conj(A)          B^H
NB. ztrmmrucnn    typed     complex    (2)     U      UT                A^H         B
NB. ztrmmrucnt    typed     complex    (2)     U      UT                A^H         B^T
NB. ztrmmrucnj    typed     complex    (2)     U      UT                A^H    conj(B)
NB. ztrmmrucnc    typed     complex    (2)     U      UT                A^H         B^H
NB. ztrmmrucun    typed     complex    (2)     U1    SUT                A^H         B
NB. ztrmmrucut    typed     complex    (2)     U1    SUT                A^H         B^T
NB. ztrmmrucuj    typed     complex    (2)     U1    SUT                A^H    conj(B)
NB. ztrmmrucuc    typed     complex    (2)     U1    SUT                A^H         B^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C  (1)
NB.   or
NB.     C := alpha * op2(B) * op1(A) + beta * C  (2)
NB.   where A is triangular, and opX(M) is either M, M^T,
NB.   conj(M) or M^H
NB.
NB. Syntax:
NB.   Сupd=. xtrmmxxxxx alpha ; AA ; B ; beta ; C
NB. where
NB.   alpha - scalar
NB.   AA    - mn×mn-matrix, contains either non-zero or both
NB.           part(s) of A
NB.   B     - mb×nb-matrix
NB.   beta  - scalar
NB.   C     - m×n-matrix
NB.   Cupd  - an updated C
NB.   A     - mn×mn-matrix, triangular
NB.   m     ≥ 0, the size of A and AA for trmmlxxxx, and the
NB.           number of rows in op2(B), C and Cupd
NB.   n     ≥ 0, the size of A and AA for trmmrxxxx, and the
NB.           number of columns in op2(B), C and Cupd
NB.   mn    = m for xtrmmlxxxx or mn = n for xtrmmrxxxx
NB.   mb    = m for xtrmmxxxxn and xtrmmxxxxj, or mb = n
NB.           otherwise
NB.   nb    = n for xtrmmxxxxn and xtrmmxxxxj, or nb = m
NB.           otherwise
NB.
NB. Notes:
NB. - monad         provides BLIS'
NB.   trmmxxxxx     bli_trmm3 (...)
NB.   dtrmmxxxxx    bli_dtrmm3(...)
NB.   ztrmmxxxxx    bli_ztrmm3(...)

trmmllnnn=:  (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmllnnt=:  (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmllnnj=:  (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmllnnc=:  (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmllnun=:  (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmllnut=:  (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmllnuj=:  (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmllnuc=:  (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlltnn=:  (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlltnt=:  (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlltnj=:  (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlltnc=:  (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlltun=:  (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlltut=:  (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlltuj=:  (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlltuc=:  (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlljnn=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlljnt=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlljnj=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlljnc=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlljun=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlljut=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlljuj=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlljuc=:  (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmllcnn=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmllcnt=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmllcnj=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmllcnc=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmllcun=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmllcut=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmllcuj=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmllcuc=:  (LEFT  , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlunnn=:  (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlunnt=:  (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlunnj=:  (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlunnc=:  (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlunun=:  (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlunut=:  (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlunuj=:  (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlunuc=:  (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlutnn=:  (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlutnt=:  (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlutnj=:  (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlutnc=:  (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlutun=:  (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlutut=:  (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlutuj=:  (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlutuc=:  (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlujnn=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlujnt=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlujnj=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlujnc=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlujun=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlujut=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlujuj=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlujuc=:  (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlucnn=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlucnt=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlucnj=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlucnc=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmlucun=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmlucut=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmlucuj=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmlucuc=:  (LEFT  , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrlnnn=:  (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrlnnt=:  (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrlnnj=:  (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrlnnc=:  (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrlnun=:  (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrlnut=:  (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrlnuj=:  (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrlnuc=:  (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrltnn=:  (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrltnt=:  (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrltnj=:  (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrltnc=:  (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrltun=:  (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrltut=:  (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrltuj=:  (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrltuc=:  (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrljnn=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrljnt=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrljnj=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrljnc=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrljun=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrljut=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrljuj=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrljuc=:  (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrlcnn=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrlcnt=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrlcnj=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrlcnc=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrlcun=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrlcut=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrlcuj=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrlcuc=:  (RIGHT , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrunnn=:  (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrunnt=:  (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrunnj=:  (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrunnc=:  (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrunun=:  (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrunut=:  (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrunuj=:  (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrunuc=:  (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrutnn=:  (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrutnt=:  (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrutnj=:  (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrutnc=:  (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrutun=:  (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrutut=:  (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrutuj=:  (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrutuc=:  (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrujnn=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrujnt=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrujnj=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrujnc=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrujun=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrujut=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrujuj=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrujuc=:  (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrucnn=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrucnt=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrucnj=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrucnc=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core
trmmrucun=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)& trmm3core
trmmrucut=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)& trmm3core
trmmrucuj=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)& trmm3core
trmmrucuc=:  (RIGHT , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)& trmm3core

dtrmmllnnn=: (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmllnnt=: (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmllnun=: (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmllnut=: (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmlltnn=: (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmlltnt=: (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmlltun=: (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmlltut=: (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmlunnn=: (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmlunnt=: (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmlunun=: (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmlunut=: (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmlutnn=: (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmlutnt=: (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmlutun=: (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmlutut=: (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmrlnnn=: (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmrlnnt=: (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmrlnun=: (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmrlnut=: (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmrltnn=: (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmrltnt=: (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmrltun=: (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmrltut=: (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmrunnn=: (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmrunnt=: (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmrunun=: (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmrunut=: (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmrutnn=: (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmrutnt=: (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&dtrmm3core
dtrmmrutun=: (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&dtrmm3core
dtrmmrutut=: (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&dtrmm3core

ztrmmllnnn=: (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmllnnt=: (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmllnnj=: (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmllnnc=: (LEFT  , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmllnun=: (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmllnut=: (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmllnuj=: (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmllnuc=: (LEFT  , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlltnn=: (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlltnt=: (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlltnj=: (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlltnc=: (LEFT  , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlltun=: (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlltut=: (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlltuj=: (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlltuc=: (LEFT  , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlljnn=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlljnt=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlljnj=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlljnc=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlljun=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlljut=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlljuj=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlljuc=: (LEFT  , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmllcnn=: (LEFT  , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmllcnt=: (LEFT  , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmllcnj=: (LEFT  , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmllcnc=: (LEFT  , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmllcun=: (LEFT  , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmllcut=: (LEFT  , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmllcuj=: (LEFT  , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmllcuc=: (LEFT  , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlunnn=: (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlunnt=: (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlunnj=: (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlunnc=: (LEFT  , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlunun=: (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlunut=: (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlunuj=: (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlunuc=: (LEFT  , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlutnn=: (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlutnt=: (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlutnj=: (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlutnc=: (LEFT  , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlutun=: (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlutut=: (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlutuj=: (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlutuc=: (LEFT  , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlujnn=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlujnt=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlujnj=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlujnc=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlujun=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlujut=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlujuj=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlujuc=: (LEFT  , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlucnn=: (LEFT  , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlucnt=: (LEFT  , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlucnj=: (LEFT  , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlucnc=: (LEFT  , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmlucun=: (LEFT  , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmlucut=: (LEFT  , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmlucuj=: (LEFT  , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmlucuc=: (LEFT  , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrlnnn=: (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrlnnt=: (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrlnnj=: (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrlnnc=: (RIGHT , LOWER ,      NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrlnun=: (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrlnut=: (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrlnuj=: (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrlnuc=: (RIGHT , LOWER ,      NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrltnn=: (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrltnt=: (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrltnj=: (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrltnc=: (RIGHT , LOWER ,         TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrltun=: (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrltut=: (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrltuj=: (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrltuc=: (RIGHT , LOWER ,         TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrljnn=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrljnt=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrljnj=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrljnc=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrljun=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrljut=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrljuj=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrljuc=: (RIGHT , LOWER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrlcnn=: (RIGHT , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrlcnt=: (RIGHT , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrlcnj=: (RIGHT , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrlcnc=: (RIGHT , LOWER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrlcun=: (RIGHT , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrlcut=: (RIGHT , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrlcuj=: (RIGHT , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrlcuc=: (RIGHT , LOWER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrunnn=: (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrunnt=: (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrunnj=: (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrunnc=: (RIGHT , UPPER ,      NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrunun=: (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrunut=: (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrunuj=: (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrunuc=: (RIGHT , UPPER ,      NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrutnn=: (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrutnt=: (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrutnj=: (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrutnc=: (RIGHT , UPPER ,         TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrutun=: (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrutut=: (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrutuj=: (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrutuc=: (RIGHT , UPPER ,         TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrujnn=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrujnt=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrujnj=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrujnc=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrujun=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrujut=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrujuj=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrujuc=: (RIGHT , UPPER , CONJ_NO_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrucnn=: (RIGHT , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrucnt=: (RIGHT , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrucnj=: (RIGHT , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrucnc=: (RIGHT , UPPER ,    CONJ_TRANSPOSE , NONUNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
ztrmmrucun=: (RIGHT , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,      NO_TRANSPOSE)&ztrmm3core
ztrmmrucut=: (RIGHT , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,         TRANSPOSE)&ztrmm3core
ztrmmrucuj=: (RIGHT , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG , CONJ_NO_TRANSPOSE)&ztrmm3core
ztrmmrucuc=: (RIGHT , UPPER ,    CONJ_TRANSPOSE ,    UNIT_DIAG ,    CONJ_TRANSPOSE)&ztrmm3core
