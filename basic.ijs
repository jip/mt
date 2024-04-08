NB. Basic linear algebra operations
NB.
NB. trsmxxxx        Solve linear monomial matrix equation
NB.                 with triangular matrix
NB.
NB. testbasicger    Test rank 1 operations with general
NB.                 matrix
NB. testbasicher    Test rank 1 operations with Hermitian
NB.                 (symmetric) matrix
NB. testbasicr      Test rank 1 operations
NB.
NB. testbasicher2   Test rank 2 operations with Hermitian
NB.                 (symmetric) matrix
NB. testbasicr2     Test rank 2 operations
NB.
NB. testbasicsyrk   Test rank k operations with symmetric
NB.                 matrix
NB. testbasicherk   Test rank k operations with Hermitian
NB.                 matrix
NB. testbasicrk     Test rank k operations
NB.
NB. testbasicsyr2k  Test rank 2k operations with symmetric
NB.                 matrix
NB. testbasicher2k  Test rank 2k operations with Hermitian
NB.                 matrix
NB. testbasicr2k    Test rank 2k operations
NB.
NB. testbasicgemv   Test matrix-vector operations with
NB.                 general matrix
NB. testbasichemv   Test matrix-vector operations with
NB.                 Hermitian (symmetric) matrix
NB. testbasictrmv   Test matrix-vector operations with
NB.                 triangular matrix
NB. testbasicmv     Test matrix-vector operations
NB.
NB. testbasicgemm   Test matrix-matrix operations with
NB.                 general matrix
NB. testbasicsymm   Test matrix-matrix operations with
NB.                 symmetric matrix
NB. testbasichemm   Test matrix-matrix operations with
NB.                 Hermitian matrix
NB. testbasictrmm   Test matrix-matrix operations with
NB.                 triangular matrix
NB. testbasicmm     Test matrix-matrix operations
NB.
NB. testbasictrsv   Test equation solver with triangular
NB.                 matrix
NB. testbasicsv     Test equation solver
NB.
NB. testbasictrsm   Test matrix equation solver with
NB.                 triangular matrix
NB. testbasicsm     Test matrix equation solver
NB.
NB. testbasic       Adv. to make verb to test basic
NB.                 operations all levels
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
NB. Concepts
NB.
NB. Conventions:
NB. 1) test suite here is aimed to benchmark BLAS and BLIS
NB.    subroutines
NB. 2) ZCHKx is not tested with real alpha and beta, just as
NB.    BLAS counterparts do

NB. =========================================================
NB. Configuration

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. extract
NB.
NB. Description:
NB.   Extract evenly spaced elements from vector
NB.
NB. Syntax:
NB.   vec=. inc extract svec
NB. where
NB.   inc  ≠ 0, the increment for the elements of vec in svec
NB.   svec - vector with evenly spaced elements of vec
NB.   vec  - vector to extract
NB.
NB. Examples:
NB.      _2 extract_mt_ 11 22 33 44 55
NB.   55 33 11
NB.
NB. Notes:
NB. - the following sentences are equivalent:
NB.     {."1 (- | inc) ]\ svec
NB.     inc (] #~ 1 ,~ 1 #~ ((% |)~ <:@#) # 1 j. <:@|@[) svec

extract=: 4 : '|.^:(0 > x) {."1 (- | x) ]\ y'

NB. ---------------------------------------------------------
NB. expand
NB.
NB. Description:
NB.   Conj. to make monad to expand vector in boxed list
NB.
NB. Syntax:
NB.   out=. (iovec expand ioinc) inp
NB. where
NB.   iovec - IO vec in inp and xvec in out
NB.   ioinc - IO inc in inp and out
NB.   vec   - n-vector, the following holds:
NB.             vec -: iovec {:: inp
NB.   inc   ≠ 0, the increment for the elements of vec, the
NB.           following holds:
NB.             inc -: ioinc {:: inp
NB.             inc -: ioinc {:: out
NB.   xvec  - (1+(n-1)*|inc|)-vector, expanded by random
NB.           number, the following holds:
NB.             xvec -: iovec {:: out
NB.   inp   - m-vector of boxes
NB.   out   - inp with vec replaced by xvec
NB.
NB. Examples:
NB.      (1 expand_mt_ 3) (_. ; 11 22 33 ; _. ; _2)
NB.   +--+----------------+--+--+
NB.   |_.|11 0.1 22 0.1 33|_.|_2|
NB.   +--+----------------+--+--+
NB. where
NB.   0.1      - random number
NB.   33 22 11 - effective vector (reversed since inc < 0)
NB.
NB. Notes:
NB. - random filler is useful to detect errors as early as
NB.   possible
NB. - filler must be non-zero to avoid problems in chkxxx
NB.   which replace 0 by 1 before comparison and division
NB. - negative inc is used by caller to reverse vec

expand=: 2 : '<@((1 j. <:@|@(n&{::)) ((#!.(?0) }:) , {:@]) m&{::) m} ]'

NB. ---------------------------------------------------------
NB. shrink
NB.
NB. Description:
NB.   Conj. to make conj. to make monad to shrink matrix in
NB.   boxed list
NB.
NB. Syntax:
NB.   out=. (getsz (iosmp shrink iomat) setsz) inp
NB. where
NB.   getsz - monad to get size of sample matrix
NB.   setsz - monad to set size of shrinked matrix
NB.   iosmp - IO sample matrix in inp
NB.   iomat - IO matrix to shrink in inp and shrinked matrix
NB.           in out
NB.   inp   - m-vector of boxes
NB.   out   - inp with matrix shrinked instead of original
NB.           one
NB.
NB. Examples:
NB.      (# (1 shrink_mt_ 2) ({."1)) 0.7 ; (i. 2 5) ; (i. 4 3) ; 1.1 ; (i. 2 3)
NB.   +---+---------+----+---+-----+
NB.   |0.7|0 1 2 3 4|0  1|1.1|0 1 2|
NB.   |   |5 6 7 8 9|3  4|   |3 4 5|
NB.   |   |         |6  7|   |     |
NB.   |   |         |9 10|   |     |
NB.   +---+---------+----+---+-----+
NB. where
NB.   #       - monad to take size (rows here) of sample
NB.   {."1    - monad to set size (columns here) of matrix
NB.   1       - IO sample matrix
NB.   2       - IO matrix to shrink
NB.   i. 2 5  - sample matrix
NB.   i. 4 3  - matrix to shrink
NB. to shrink (i. 4 3) to (2 {."1 i. 4 3)
NB.
NB. Notes:
NB. - shrink is a double conjunction

shrink=: {{2 : ('(((u@(' , (": m) , '&{::)) <@:v (' , (": n) , '&{::)) ' , (": n) , '} ])')}}

NB. ---------------------------------------------------------
NB. reverse
NB.
NB. Description:
NB.   Conj. to make monad to reverse vector in boxed list
NB.
NB. Syntax:
NB.   out=. (iovec reverse ioinc) inp
NB. where
NB.   iovec - IO vec in inp and rvec in out
NB.   ioinc - IO inc in inp and out
NB.   vec   - n-vector, the following holds:
NB.             vec -: iovec {:: inp
NB.   inc   ≠ 0, the increment for the elements of vec, the
NB.           following holds:
NB.             inc -: ioinc {:: inp
NB.             inc -: ioinc {:: out
NB.   rvec  - vec reversed iif (inc < 0), the following
NB.           holds:
NB.             rvec -: iovec {:: out
NB.   inp   - m-vector of boxes
NB.   out   - inp with vec replaced by rvec
NB.
NB. Examples:
NB.      (1 reverse_mt_ 3) (_. ; 11 22 33 44 55 ; _. ; _2)
NB.   +--+----------------+--+--+
NB.   |_.|55 44 33 22 11|_.|_2|
NB.   +--+----------------+--+--+
NB. where
NB.   55 33 11 - effective vector (reversed since inc < 0)

reverse=: 2 : '[:`]`(|.&.>&.(m&{))@.(*@(n&{::))'

NB. ---------------------------------------------------------
NB. ger
NB.
NB. Description:
NB.   Adv. to make monad to perform the rank 1 operation:
NB.     A := alpha * x * op(y) + A
NB. where
NB.   op(y) - either y^T or y^H
NB.
NB. Syntax:
NB.   Aupd=. (mul ger) alpha ; x ; incx ; y ; incy ; A
NB. where
NB.   mul - dyad to compute the product (alpha * op(y)), is
NB.         called as:
NB.           product=. alpha mul y
NB.         and is one of:
NB.           *      NB. op(y) = y^T
NB.           (* +)  NB. op(y) = y^H

ger=: 1 : 0
  'alpha xx incx y incy A'=. y
  xx=. incx extract_mt_ xx
  y=.  incy extract_mt_ y
  A=. A + xx */ alpha u y
)

NB. ---------------------------------------------------------
NB. Monad    op(y)
NB. gerc     y^H
NB. geru     y^T
NB.
NB. Description:
NB.   Performs the rank 1 operation:
NB.     A := alpha * x * op(y) + A
NB. where
NB.   op(y) - either y^T or y^H
NB.
NB. Syntax:
NB.   Aupd=. gerx alpha ; x ; incx ; y ; incy ; A
NB. where
NB.   alpha - scalar
NB.   x     - (1+(m-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   y     - (1+(n-1)*|incy|)-vector
NB.   incy  ≠ 0, the increment for the elements of y
NB.   A     - m×n-matrix
NB.   Aupd  - an updated A
NB.   m     ≥ 0, the number of rows in A and Aupd
NB.   n     ≥ 0, the number of columns in A and Aupd
NB.
NB. Notes:
NB. - gerc implements BLAS' ZGERC
NB. - geru implements BLAS' DGER ZGERU
NB. - reference implementation

gerc=: (* +) ger
geru=:  *    ger

NB. ---------------------------------------------------------
NB. her
NB.
NB. Description:
NB.   Conj. to make monad to perform the hermitian
NB.   (symmetric) rank 1 operation:
NB.     A := alpha * x * op(x) + A
NB. where
NB.   A     - Hermitian (symmetric)
NB.   op(x) - either x^T or x^H
NB.
NB. Syntax:
NB.   AAupd=. (ref her mul) alpha ; x ; incx ; AA
NB. where
NB.   ref - monad to pick a triangular part, is one of:
NB.           trlpick  NB. LT
NB.           trupick  NB. UT
NB.   mul - dyad to define the form of op(A), is one of:
NB.           *        NB. the symmetric operation: op(A) := A^T
NB.           (* +)    NB. the hermitian operation: op(A) := A^H

her=: 2 : 0
  'alpha y incy AA'=. y
  y=. incy extract_mt_ y
  AA=. AA + u alpha (] */ v) y
)

NB. ---------------------------------------------------------
NB. Monad    A     R/W in A    op(x)
NB. syrl     SY    LT          x^T
NB. syru     SY    UT          x^T
NB. herl     HE    LT          x^H
NB. heru     HE    UT          x^H
NB.
NB. Description:
NB.   Performs the hermitian (symmetric) rank 1 operation:
NB.     A := alpha * x * op(x) + A
NB. where
NB.   A - Hermitian (symmetric)
NB.
NB. Syntax:
NB.   AAupd=. xxrx alpha ; x ; incx ; AA
NB. where
NB.   alpha - scalar, real
NB.   x     - (1+(n-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   AA    - n×n-matrix, contains either LT or UT or both
NB.           part(s) of A
NB.   AAupd - AA with either LT (for xxrl) or UT (for xxru)
NB.           updated
NB.   A     - n×n-matrix, Hermitian (symmetric)
NB.   n     ≥ 0, the size of A, AA and AAupd
NB.
NB. Notes:
NB. - syrl models BLAS' DSYR('L',...) with the following
NB.   extension: A can have complex datatype, too
NB. - syru models BLAS' DSYR('U',...) with the following
NB.   extension: A can have complex datatype, too
NB. - herl models BLAS' ZHER('L',...)
NB. - heru models BLAS' ZHER('U',...)
NB. - reference implementation

syrl=: trlpick_mt_ her  *
syru=: trupick_mt_ her  *

herl=: trlpick_mt_ her (* +)
heru=: trupick_mt_ her (* +)

NB. ---------------------------------------------------------
NB. her2
NB.
NB. Description:
NB.   Conj. to make monad to perform the hermitian
NB.   (symmetric) rank 2 operation:
NB.     A := alpha * x * op1(y) + op2(alpha) * y * op1(x) + A
NB. where
NB.   A          - Hermitian (symmetric)
NB.   op1(x)     - either x^T or x^H
NB.   op2(alpha) - either alpha or conj(alpha)
NB.
NB. Syntax:
NB.   AAupd=. (ref her2 trans) alpha ; x ; incx ; y ; incy ; AA
NB. where
NB.   ref   - monad to pick a triangular part, is one of:
NB.             trlpick  NB. LT
NB.             trupick  NB. UT
NB.   trans - monad to define the form of op1(v) and op2(s),
NB.           is one of:
NB.             |:       NB. the symmetric operation: op1(v) = v^T, op2(s) = s
NB.             ct       NB. the hermitian operation: op1(v) = v^H, op2(s) = conj(s)

her2=: 2 : 0
  'alpha xx incx y incy AA'=. y
  xx=. incx extract_mt_ xx
  y=.  incy extract_mt_ y
  AA=. AA + u (+ v)~ xx */ alpha * v y
)

NB. ---------------------------------------------------------
NB. Monad    A     R/W in A    op1(v)    op2(alpha)
NB. syr2l    SY    LT          v^T       alpha
NB. syr2u    SY    UT          v^T       alpha
NB. her2l    HE    LT          v^H       conj(alpha)
NB. her2u    HE    UT          v^H       conj(alpha)
NB.
NB. Description:
NB.   Performs the hermitian (symmetric) rank 2 operation:
NB.     A := alpha * x * op1(y) + op2(alpha) * y * op1(x) + A
NB. where
NB.   A - Hermitian (symmetric)
NB.
NB. Syntax:
NB.   AAupd=. xxr2x alpha ; x ; incx ; y ; incy ; AA
NB. where
NB.   alpha - scalar
NB.   x     - (1+(n-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   y     - (1+(n-1)*|incy|)-vector
NB.   incy  ≠ 0, the increment for the elements of y
NB.   AA    - n×n-matrix, contains either LT or UT or both
NB.           part(s) of A
NB.   AAupd - AA with either LT (for xxr2l) or UT (for
NB.           xxr2u) updated
NB.   A     - n×n-matrix, Hermitian (symmetric)
NB.   n     ≥ 0, the size of A, AA and AAupd
NB.
NB. Notes:
NB. - syr2l models BLAS' DSYR2('L',...) with the following
NB.   extension: A can have complex datatype, too
NB. - syr2u models BLAS' DSYR2('U',...) with the following
NB.   extension: A can have complex datatype, too
NB. - her2l models BLAS' ZHER2('L',...)
NB. - her2u models BLAS' ZHER2('U',...)
NB. - reference implementation

syr2l=: trlpick_mt_ her2 |:
syr2u=: trupick_mt_ her2 |:

her2l=: trlpick_mt_ her2 ct_mt_
her2u=: trupick_mt_ her2 ct_mt_

NB. ---------------------------------------------------------
NB. herk
NB.
NB. Description:
NB.   Conj. to make monad to perform the hermitian
NB.   (symmetric) rank k operation:
NB.     C := alpha * A * op(A) + beta * C
NB.   or
NB.     C := alpha * op(A) * A + beta * C
NB. where
NB.   C     - Hermitian (symmetric)
NB.   op(A) - either A^T or A^H
NB.
NB. Syntax:
NB.   CCupd=. (ctp herk mul) alpha ; A ; beta ; CC
NB. where
NB.   ctp - dyad to compose a matrix from triangular parts,
NB.         is one of:
NB.           slxuy  NB. take SLT from CC, and UT from the matrix computed
NB.           suxly  NB. take SUT from CC, and LT from the matrix computed
NB.   mul - dyad to compute the product either (A * op(A)) or
NB.         (op(A) * A), is called as:
NB.           product=. mul A

herk=: 2 : '3&{:: u (0&{:: (* v~) 1&{::) + 2&{:: * 3&{::'

NB. ---------------------------------------------------------
NB. Monad     C     alpha,beta    R/W in C    op1(A)    op2(A)
NB. syrkln    SY    any           LT          A         A^T
NB. syrklt    SY    any           LT          A^T       A
NB. syrkun    SY    any           UT          A         A^T
NB. syrkut    SY    any           UT          A^T       A
NB. herkln    HE    real          LT          A         A^H
NB. herklc    HE    real          LT          A^H       A
NB. herkun    HE    real          UT          A         A^H
NB. herkuc    HE    real          UT          A^H       A
NB.
NB. Description:
NB.   Performs the hermitian (symmetric) rank k operation:
NB.     C := alpha * op1(A) * op2(A) + beta * C
NB. where
NB.   C - Hermitian (symmetric)
NB.
NB. Syntax:
NB.   CCupd=. xxrkxx alpha ; A ; beta ; CC
NB. where
NB.   alpha - scalar, must be real for herkxx
NB.   A     - na×ka-matrix
NB.   beta  - scalar, must be real for herkxx
NB.   CC    - n×n-matrix, contains either LT or UT or both
NB.           part(s) of C
NB.   CCupd - CC with either LT (for xxrklx) or UT (for
NB.           xxrkux) updated
NB.   C     - n×n-matrix, Hermitian (symmetric)
NB.   n     ≥ 0, the size of C, CC and CCupd and the number
NB.           of rows or columns in A
NB.   k     ≥ 0, the number of columns or rows in A
NB.   ka    = k for xxrkxn or ka = n otherwise
NB.   na    = n for xxrkxn or na = k otherwise
NB.
NB. Notes:
NB. - monad     models BLAS'
NB.   syrkln    xSYRK('L','N',...)
NB.   syrklt    xSYRK('L','T',...) and DSYRK('L','C',...)
NB.   syrkun    xSYRK('U','N',...)
NB.   syrkut    xSYRK('U','T',...) and DSYRK('U','C',...)
NB.   herkln    ZHERK('L','N',...)
NB.   herklc    ZHERK('L','C',...)
NB.   herkun    ZHERK('U','N',...)
NB.   herkuc    ZHERK('U','C',...)
NB. - reference implementation

syrkln=: suxly_mt_ herk (mp_mt_  |:    )
syrklt=: suxly_mt_ herk (mp_mt_~ |:    )
syrkun=: slxuy_mt_ herk (mp_mt_  |:    )
syrkut=: slxuy_mt_ herk (mp_mt_~ |:    )

herkln=: suxly_mt_ herk (mp_mt_  ct_mt_)
herklc=: suxly_mt_ herk (mp_mt_~ ct_mt_)
herkun=: slxuy_mt_ herk (mp_mt_  ct_mt_)
herkuc=: slxuy_mt_ herk (mp_mt_~ ct_mt_)

NB. ---------------------------------------------------------
NB. her2k
NB.
NB. Description:
NB.   Conj. to make monad to perform the hermitian
NB.   (symmetric) rank 2k operation:
NB.     C := alpha * A * op1(B) + op2(alpha) * B * op1(A) + beta * C  (1)
NB.   or
NB.     C := alpha * op1(A) * B + op2(alpha) * op1(B) * A + beta * C  (2)
NB. where
NB.   C          - Hermitian (symmetric)
NB.   op1(M)     - either M^T or M^H
NB.   op2(alpha) - either alpha or conj(alpha)
NB.
NB. Syntax:
NB.   CCupd=. (ctp`trans her2k kind) alpha ; A ; B ; beta ; CC
NB. where
NB.   ctp   - dyad to compose a matrix from triangular parts,
NB.           is one of:
NB.             slxuy  NB. take SLT from CC, and UT from the matrix computed
NB.             suxly  NB. take SUT from CC, and LT from the matrix computed
NB.   trans - monad to define the form of op1(M) and op2(s),
NB.           is one of:
NB.             |:     NB. the symmetric operation: op1(M) = M^T, op2(s) = s
NB.             ct     NB. the hermitian operation: op1(M) = M^H, op2(s) = conj(s)
NB.   kind  - boolean scalar to define operation:
NB.             0      NB. (2)
NB.             1      NB. (1)
NB.
NB. Notes:
NB. - her2k's design solves a problem: how to allow C1 to see
NB.   V2 in the train (V0 C1 V2 A3), the solution is: send V2
NB.   not (V2 A3) into C1, implement A3 functionality inside
NB.   C1 inline, use switch N4 to control A3 behavior, the
NB.   resulting train becomes (V0`V2 C1 N4)

her2k=: 2 : '4&{:: m@.0 (0&{:: (+ m@.1)@:* 1&{:: (mp_mt_~ m@.1)~`(mp_mt_ m@.1)@.n 2&{::) + 3&{:: * 4&{::'

NB. ---------------------------------------------------------
NB. Monad      C     beta    R/W in C    op1(M)    op2(M)    op3(alpha)
NB. syr2kln    SY    any     LT          M         M^T       alpha
NB. syr2klt    SY    any     LT          M^T       M         alpha
NB. syr2kun    SY    any     UT          M         M^T       alpha
NB. syr2kut    SY    any     UT          M^T       M         alpha
NB. her2kln    HE    real    LT          M         M^H       conj(alpha)
NB. her2klc    HE    real    LT          M^H       M         conj(alpha)
NB. her2kun    HE    real    UT          M         M^H       conj(alpha)
NB. her2kuc    HE    real    UT          M^H       M         conj(alpha)
NB.
NB. Description:
NB.   Performs the hermitian (symmetric) rank 2k operation:
NB.     C := alpha * op1(A) * op2(B) + op3(alpha) * op1(B) * op2(A) + beta * C
NB. where
NB.   C - Hermitian (symmetric)
NB.
NB. Syntax:
NB.   CCupd=. xxr2kxx alpha ; A ; B ; beta ; CC
NB. where
NB.   alpha - scalar
NB.   A     - nab×kab-matrix
NB.   B     - nab×kab-matrix
NB.   beta  - scalar, must be real for her2kxx
NB.   CC    - n×n-matrix, contains either LT or UT or both
NB.           part(s) of C
NB.   CCupd - CC with either LT (for xxr2klx) or UT (for
NB.           xxr2kux) updated
NB.   C     - n×n-matrix, Hermitian (symmetric)
NB.   n     ≥ 0, the size of C, CC and CCupd and the number
NB.           of rows or columns in A and B
NB.   k     ≥ 0, the number of columns or rows in A and B
NB.   kab   = k for xxr2kxn or kab = n otherwise
NB.   nab   = n for xxr2kxn or nab = k otherwise
NB.
NB. Notes:
NB. - monad      models BLAS'
NB.   syr2kln    xSYR2K('L','N',...)
NB.   syr2klt    xSYR2K('L','T',...) and DSYR2K('L','C',...)
NB.   syr2kun    xSYR2K('U','N',...)
NB.   syr2kut    xSYR2K('U','T',...) and DSYR2K('U','C',...)
NB.   her2kln    ZHER2K('L','N',...)
NB.   her2klc    ZHER2K('L','C',...)
NB.   her2kun    ZHER2K('U','N',...)
NB.   her2kuc    ZHER2K('U','C',...)
NB. - reference implementation

syr2kln=: suxly_mt_`|:     her2k 1
syr2klt=: suxly_mt_`|:     her2k 0
syr2kun=: slxuy_mt_`|:     her2k 1
syr2kut=: slxuy_mt_`|:     her2k 0

her2kln=: suxly_mt_`ct_mt_ her2k 1
her2klc=: suxly_mt_`ct_mt_ her2k 0
her2kun=: slxuy_mt_`ct_mt_ her2k 1
her2kuc=: slxuy_mt_`ct_mt_ her2k 0

NB. ---------------------------------------------------------
NB. gemv
NB.
NB. Description:
NB.   Adv. to make monad to perform the matrix-vector
NB.   operation:
NB.     y := alpha * op(A) * x + beta * y
NB.
NB. Syntax:
NB.   yupd=. (trans gemv) alpha ; A ; x ; incx ; beta ; y ; incy
NB. where
NB.   trans - monad to define the form of op(A), is one of:
NB.             ]   NB. op(A) := A
NB.             |:  NB. op(A) := A^T
NB.             ct  NB. op(A) := A^H

gemv=: 1 : 0
  'alpha A xx incx beta ybak incy'=. y
  xx=. incx extract_mt_ xx
  y=.  incy extract_mt_ ybak
  y=.  ((u A) mp_mt_ alpha * xx) + beta * y
  y=.  y (incy ([ ((* |)~ i.) negneg_mt_) #@[)} ybak
)

NB. ---------------------------------------------------------
NB. Monad    op(A)
NB. gemvn    A
NB. gemvt    A^T
NB. gemvc    A^H
NB.
NB. Description:
NB.   Performs the matrix-vector operation:
NB.     y := alpha * op(A) * x + beta * y
NB.
NB. Syntax:
NB.   yupd=. gemvx alpha ; A ; x ; incx ; beta ; y ; incy
NB. where
NB.   alpha - scalar
NB.   A     - m×n-matrix
NB.   x     - (1+(kx-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   beta  - scalar
NB.   y     - (1+(ky-1)*|incy|)-vector
NB.   incy  ≠ 0, the increment for the elements of y
NB.   yupd  - an updated y
NB.   kx    = n for gemvn or kx = m otherwise
NB.   ky    = m for gemvn or ky = n otherwise
NB.   m     ≥ 0, the number of rows in A
NB.   n     ≥ 0, the number of columns in A
NB.
NB. Notes:
NB. - monad    models BLAS'
NB.   gemvn    xGEMV('N',...)
NB.   gemvt    xGEMV('T',...)
NB.   gemvc    DGEMV('T',...) and ZGEMV('C',...)
NB. - reference implementation

gemvn=: ]      gemv
gemvt=: |:     gemv
gemvc=: ct_mt_ gemv

NB. ---------------------------------------------------------
NB. hemv
NB.
NB. Description:
NB.   Adv. to make monad to perform the matrix-vector
NB.   operation:
NB.     y := alpha * A * x + beta * y
NB. where
NB.   A - Hermitian (symmetric)
NB.
NB. Syntax:
NB.   yupd=. (ref hemv) alpha ; AA ; x ; incx ; beta ; y ; incy
NB. where
NB.   ref - monad to restore A from triangular part, is one
NB.         of:
NB.           sy4gel  NB. LT, A is symmetric
NB.           sy4geu  NB. UT, A is symmetric
NB.           he4gel  NB. LT, A is Hermitian
NB.           he4geu  NB. UT, A is Hermitian

hemv=: gemv

NB. ---------------------------------------------------------
NB. Monad    A     Reads in A
NB. symvl    SY    LT
NB. symvu    SY    UT
NB. hemvl    HE    LT
NB. hemvu    HE    UT
NB.
NB. Description:
NB.   Performs the matrix-vector operation:
NB.     y := alpha * A * x + beta * y
NB. where
NB.   A - Hermitian (symmetric)
NB.
NB. Syntax:
NB.   yupd=. xxmvx alpha ; AA ; x ; incx ; beta ; y ; incy
NB. where
NB.   alpha - scalar
NB.   AA    - n×n-matrix, contains either LT or UT or both
NB.           part(s) of A
NB.   x     - (1+(n-1)*|incx|)-vector
NB.   incx  ≠ 0, the increment for the elements of x
NB.   beta  - scalar
NB.   y     - (1+(n-1)*|incy|)-vector
NB.   incy  ≠ 0, the increment for the elements of y
NB.   yupd  - an updated y
NB.   A     - n×n-matrix, Hermitian (symmetric)
NB.   n     ≥ 0, the size of A and AA
NB.
NB. Notes:
NB. - symvl models BLAS' DSYMV('L',...) with the following
NB.   extension: A can have complex datatype, too
NB. - symvu models BLAS' DSYMV('U',...) with the following
NB.   extension: A can have complex datatype, too
NB. - hemvl models BLAS' ZHEMV('L',...)
NB. - hemvu models BLAS' ZHEMV('U',...)
NB. - reference implementation

symvl=: sy4gel_mt_ hemv
symvu=: sy4geu_mt_ hemv

hemvl=: he4gel_mt_ hemv
hemvu=: he4geu_mt_ hemv

NB. ---------------------------------------------------------
NB. trmv
NB.
NB. Description:
NB.   Adv. to make monad to perform the matrix-vector
NB.   operation:
NB.     x := op(A) * x
NB. where
NB.   A - triangular
NB.
NB. Syntax:
NB.   xupd=. ((mp_mt_~ trans@ref) trmv) AA ; x ; incx
NB. where
NB.   ref   - monad to restore A from triangular part, is one
NB.           of:
NB.             trlpick   NB.  LT, A is L
NB.             trl1pick  NB. SLT, A is L1
NB.             trupick   NB.  UT, A is U
NB.             tru1pick  NB. SUT, A is U1
NB.   trans - monad to define the form of op(A), is one of:
NB.             ]         NB. op(A) := A
NB.             |:        NB. op(A) := A^T
NB.             ct        NB. op(A) := A^H

trmv=: 1 : 0
  'AA ybak incy'=. y
  y=. incy extract_mt_ ybak
  y=. y u AA
  y=. y (incy ([ ((* |)~ i.) negneg_mt_) #@[)} ybak
)

NB. ---------------------------------------------------------
NB. Monad      A     Reads in A    op(A)
NB. trmvlnn    L      LT           A
NB. trmvlnu    L1    SLT           A
NB. trmvltn    L      LT           A^T
NB. trmvltu    L1    SLT           A^T
NB. trmvlcn    L      LT           A^H
NB. trmvlcu    L1    SLT           A^H
NB. trmvunn    U      UT           A
NB. trmvunu    U1    SUT           A
NB. trmvutn    U      UT           A^T
NB. trmvutu    U1    SUT           A^T
NB. trmvucn    U      UT           A^H
NB. trmvucu    U1    SUT           A^H
NB.
NB. Description:
NB.   Performs the matrix-vector operation:
NB.     x := op(A) * x
NB. where
NB.   A - triangular
NB.
NB. Syntax:
NB.   xupd=. trmvxxx AA ; x ; incx
NB. where
NB.   AA   - n×n-matrix, contains either non-zero or both
NB.          part(s) of A
NB.   x    - (1+(n-1)*|incx|)-vector
NB.   incx ≠ 0, the increment for the elements of x
NB.   xupd - an updated x
NB.   A    - n×n-matrix, triangular
NB.   n    ≥ 0, the size of A and AA
NB.
NB. Notes:
NB. - monad      models BLAS'
NB.   trmvlnn    xTRMV('L','N','N',...)
NB.   trmvlnu    xTRMV('L','N','U',...)
NB.   trmvltn    xTRMV('L','T','N',...)
NB.   trmvltu    xTRMV('L','T','U',...)
NB.   trmvlcn    xTRMV('L','C','N',...)
NB.   trmvlcu    xTRMV('L','C','U',...)
NB.   trmvunn    xTRMV('U','N','N',...)
NB.   trmvunu    xTRMV('U','N','U',...)
NB.   trmvutn    xTRMV('U','T','N',...)
NB.   trmvutu    xTRMV('U','T','U',...)
NB.   trmvucn    xTRMV('U','C','N',...)
NB.   trmvucu    xTRMV('U','C','U',...)
NB. - reference implementation

trmvlnn=: (mp_mt_~        trlpick_mt_ ) trmv
trmvlnu=: (mp_mt_~        trl1pick_mt_) trmv
trmvltn=: (mp_mt_~ |:    @trlpick_mt_ ) trmv
trmvltu=: (mp_mt_~ |:    @trl1pick_mt_) trmv
trmvlcn=: (mp_mt_~ ct_mt_@trlpick_mt_ ) trmv
trmvlcu=: (mp_mt_~ ct_mt_@trl1pick_mt_) trmv
trmvunn=: (mp_mt_~        trupick_mt_ ) trmv
trmvunu=: (mp_mt_~        tru1pick_mt_) trmv
trmvutn=: (mp_mt_~ |:    @trupick_mt_ ) trmv
trmvutu=: (mp_mt_~ |:    @tru1pick_mt_) trmv
trmvucn=: (mp_mt_~ ct_mt_@trupick_mt_ ) trmv
trmvucu=: (mp_mt_~ ct_mt_@tru1pick_mt_) trmv

NB. ---------------------------------------------------------
NB. gemm
NB.
NB. Description:
NB.   Adv. to make monad to perform the matrix-matrix
NB.   operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C
NB. where
NB.   opX(M) - either M, M^T, conj(M) or M^H
NB.
NB. Syntax:
NB.   Cupd=. (mul gemm) alpha ; A ; B ; beta ; C
NB. where
NB.   mul - dyad to compute the product (op1(A) * op2(B)), is
NB.         called as:
NB.           product=. A mul B

gemm=: 1 : '(0&{:: * 1&{:: u 2&{::) + 3&{:: * 4&{::'

NB. ---------------------------------------------------------
NB. Monad     op1(A)      op2(B)
NB. gemmnn         A           B
NB. gemmnt         A           B^T
NB. gemmnj         A      conj(B)
NB. gemmnc         A           B^H
NB. gemmtn         A^T         B
NB. gemmtt         A^T         B^T
NB. gemmtj         A^T    conj(B)
NB. gemmtc         A^T         B^H
NB. gemmjn    conj(A)          B
NB. gemmjt    conj(A)          B^T
NB. gemmjj    conj(A)     conj(B)
NB. gemmjc    conj(A)          B^H
NB. gemmcn         A^H         B
NB. gemmct         A^H         B^T
NB. gemmcj         A^H    conj(B)
NB. gemmcc         A^H         B^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C
NB.
NB. Syntax:
NB.   Cupd=. gemmxx alpha ; A ; B ; beta ; C
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
NB.   ma    = m for gemmnx and gemmjx, or ma = k otherwise
NB.   ka    = k for gemmnx and gemmjx, or ka = m otherwise
NB.   kb    = k for gemmxn and gemmxj, or kb = n otherwise
NB.   nb    = n for gemmxn and gemmxj, or nb = k otherwise
NB.
NB. Notes:
NB. - monad     models BLAS'
NB.   gemmnn    xGEMM('N','N',...)
NB.   gemmnt    xGEMM('N','T',...)
NB.   gemmnc    xGEMM('N','C',...)
NB.   gemmtn    xGEMM('T','N',...)
NB.   gemmtt    xGEMM('T','T',...)
NB.   gemmtc    xGEMM('T','C',...)
NB.   gemmcn    xGEMM('C','N',...)
NB.   gemmct    xGEMM('C','T',...)
NB.   gemmcc    xGEMM('C','C',...)
NB. - reference implementation

gemmnn=:   mp_mt_                   gemm
gemmnt=:  (mp_mt_  |:    )          gemm
gemmnj=:  (mp_mt_  +     )          gemm
gemmnc=:  (mp_mt_  ct_mt_)          gemm

gemmtn=:  (mp_mt_~ |:    )~         gemm
gemmtt=:   mp_mt_& |:               gemm
gemmtj=: ((mp_mt_~ |:    )~ +     ) gemm
gemmtc=: ((mp_mt_~ |:    )~ ct_mt_) gemm

gemmjn=:  (mp_mt_~ +     )~         gemm
gemmjt=: ((mp_mt_~ +     )~ |:    ) gemm
gemmjj=:   mp_mt_&:+                gemm
gemmjc=: ((mp_mt_~ +     )~ ct_mt_) gemm

gemmcn=:  (mp_mt_~ ct_mt_)~         gemm
gemmct=: ((mp_mt_~ ct_mt_)~ |:    ) gemm
gemmcj=: ((mp_mt_~ ct_mt_)~ +     ) gemm
gemmcc=:   mp_mt_& ct_mt_           gemm

NB. ---------------------------------------------------------
NB. hemm
NB.
NB. Description:
NB.   Adv. to make monad to perform the matrix-matrix
NB.   operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C  (1)
NB.   or
NB.     C := alpha * op2(B) * op1(A) + beta * C  (2)
NB. where
NB.   A      - Hermitian (symmetric)
NB.   op1(A) - either A or conj(A)
NB.   op2(B) - either B, B^T, conj(B) or B^H
NB.
NB. Syntax:
NB.   Cupd=. (mul hemm) alpha ; AA ; B ; beta ; C
NB. where
NB.   mul - dyad to restore A from triangular part, and to
NB.         compute the product either (op1(A) * op2(B)) or
NB.         (op2(B) * op1(A)), is called as:
NB.           product=. AA mul B

hemm=: gemm

NB. ---------------------------------------------------------
NB. Monad       A     Side    Reads in A    op1(A)     op2(B)
NB. symmllnn    SY    (1)     LT                 A          B
NB. symmllnt    SY    (1)     LT                 A          B^T
NB. symmllnj    SY    (1)     LT                 A     conj(B)
NB. symmllnc    SY    (1)     LT                 A          B^H
NB. symmlljn    SY    (1)     LT            conj(A)         B
NB. symmlljt    SY    (1)     LT            conj(A)         B^T
NB. symmlljj    SY    (1)     LT            conj(A)    conj(B)
NB. symmlljc    SY    (1)     LT            conj(A)         B^H
NB. symmlunn    SY    (1)     UT                 A          B
NB. symmlunt    SY    (1)     UT                 A          B^T
NB. symmlunj    SY    (1)     UT                 A     conj(B)
NB. symmlunc    SY    (1)     UT                 A          B^H
NB. symmlujn    SY    (1)     UT            conj(A)         B
NB. symmlujt    SY    (1)     UT            conj(A)         B^T
NB. symmlujj    SY    (1)     UT            conj(A)    conj(B)
NB. symmlujc    SY    (1)     UT            conj(A)         B^H
NB. symmrlnn    SY    (2)     LT                 A          B
NB. symmrlnt    SY    (2)     LT                 A          B^T
NB. symmrlnj    SY    (2)     LT                 A     conj(B)
NB. symmrlnc    SY    (2)     LT                 A          B^H
NB. symmrljn    SY    (2)     LT            conj(A)         B
NB. symmrljt    SY    (2)     LT            conj(A)         B^T
NB. symmrljj    SY    (2)     LT            conj(A)    conj(B)
NB. symmrljc    SY    (2)     LT            conj(A)         B^H
NB. symmrunn    SY    (2)     UT                 A          B
NB. symmrunt    SY    (2)     UT                 A          B^T
NB. symmrunj    SY    (2)     UT                 A     conj(B)
NB. symmrunc    SY    (2)     UT                 A          B^H
NB. symmrujn    SY    (2)     UT            conj(A)         B
NB. symmrujt    SY    (2)     UT            conj(A)         B^T
NB. symmrujj    SY    (2)     UT            conj(A)    conj(B)
NB. symmrujc    SY    (2)     UT            conj(A)         B^H
NB. hemmllnn    HE    (1)     LT                 A          B
NB. hemmllnt    HE    (1)     LT                 A          B^T
NB. hemmllnj    HE    (1)     LT                 A     conj(B)
NB. hemmllnc    HE    (1)     LT                 A          B^H
NB. hemmlljn    HE    (1)     LT            conj(A)         B
NB. hemmlljt    HE    (1)     LT            conj(A)         B^T
NB. hemmlljj    HE    (1)     LT            conj(A)    conj(B)
NB. hemmlljc    HE    (1)     LT            conj(A)         B^H
NB. hemmlunn    HE    (1)     UT                 A          B
NB. hemmlunt    HE    (1)     UT                 A          B^T
NB. hemmlunj    HE    (1)     UT                 A     conj(B)
NB. hemmlunc    HE    (1)     UT                 A          B^H
NB. hemmlujn    HE    (1)     UT            conj(A)         B
NB. hemmlujt    HE    (1)     UT            conj(A)         B^T
NB. hemmlujj    HE    (1)     UT            conj(A)    conj(B)
NB. hemmlujc    HE    (1)     UT            conj(A)         B^H
NB. hemmrlnn    HE    (2)     LT                 A          B
NB. hemmrlnt    HE    (2)     LT                 A          B^T
NB. hemmrlnj    HE    (2)     LT                 A     conj(B)
NB. hemmrlnc    HE    (2)     LT                 A          B^H
NB. hemmrljn    HE    (2)     LT            conj(A)         B
NB. hemmrljt    HE    (2)     LT            conj(A)         B^T
NB. hemmrljj    HE    (2)     LT            conj(A)    conj(B)
NB. hemmrljc    HE    (2)     LT            conj(A)         B^H
NB. hemmrunn    HE    (2)     UT                 A          B
NB. hemmrunt    HE    (2)     UT                 A          B^T
NB. hemmrunj    HE    (2)     UT                 A     conj(B)
NB. hemmrunc    HE    (2)     UT                 A          B^H
NB. hemmrujn    HE    (2)     UT            conj(A)         B
NB. hemmrujt    HE    (2)     UT            conj(A)         B^T
NB. hemmrujj    HE    (2)     UT            conj(A)    conj(B)
NB. hemmrujc    HE    (2)     UT            conj(A)         B^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C  (1)
NB.   or
NB.     C := alpha * op2(B) * op1(A) + beta * C  (2)
NB. where
NB.   A      - Hermitian (symmetric)
NB.   op1(A) - either A or conj(A)
NB.   op2(B) - either B, B^T, conj(B) or B^H
NB.
NB. Syntax:
NB.   Cupd=. xxmmxxxx alpha ; AA ; B ; beta ; C
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
NB.   mn    = m for xxmmlxxx or mn = n for xxmmrxxx
NB.
NB. Notes:
NB. - monad       models BLAS'
NB.   symmllnn    xSYMM('L','L',...)
NB.   symmlunn    xSYMM('L','U',...)
NB.   symmrlnn    xSYMM('R','L',...)
NB.   symmrunn    xSYMM('R','U',...)
NB.   hemmllnn    ZHEMM('L','L',...)
NB.   hemmlunn    ZHEMM('L','U',...)
NB.   hemmrlnn    ZHEMM('R','L',...)
NB.   hemmrunn    ZHEMM('R','U',...)
NB. - reference implementation

symmllnn=:  (mp_mt_~    sy4gel_mt_)~         hemm
symmllnt=: ((mp_mt_~    sy4gel_mt_)~ |:    ) hemm
symmllnj=: ((mp_mt_~    sy4gel_mt_)~ +     ) hemm
symmllnc=: ((mp_mt_~    sy4gel_mt_)~ ct_mt_) hemm

symmlljn=:  (mp_mt_~  +@sy4gel_mt_)~         hemm
symmlljt=: ((mp_mt_~  +@sy4gel_mt_)~ |:    ) hemm
symmlljj=:  (mp_mt_~&:+ sy4gel_mt_)~         hemm
symmlljc=: ((mp_mt_~  +@sy4gel_mt_)~ ct_mt_) hemm

symmlunn=:  (mp_mt_~    sy4geu_mt_)~         hemm
symmlunt=: ((mp_mt_~    sy4geu_mt_)~ |:    ) hemm
symmlunj=: ((mp_mt_~    sy4geu_mt_)~ +     ) hemm
symmlunc=: ((mp_mt_~    sy4geu_mt_)~ ct_mt_) hemm

symmlujn=:  (mp_mt_~  +@sy4geu_mt_)~         hemm
symmlujt=: ((mp_mt_~  +@sy4geu_mt_)~ |:    ) hemm
symmlujj=:  (mp_mt_~&:+ sy4geu_mt_)~         hemm
symmlujc=: ((mp_mt_~  +@sy4geu_mt_)~ ct_mt_) hemm

symmrlnn=:  (mp_mt_     sy4gel_mt_)~         hemm
symmrlnt=: ((mp_mt_     sy4gel_mt_)~ |:    ) hemm
symmrlnj=: ((mp_mt_     sy4gel_mt_)~ +     ) hemm
symmrlnc=: ((mp_mt_     sy4gel_mt_)~ ct_mt_) hemm

symmrljn=:  (mp_mt_   +@sy4gel_mt_)~         hemm
symmrljt=: ((mp_mt_   +@sy4gel_mt_)~ |:    ) hemm
symmrljj=:  (mp_mt_ &:+ sy4gel_mt_)~         hemm
symmrljc=: ((mp_mt_   +@sy4gel_mt_)~ ct_mt_) hemm

symmrunn=:  (mp_mt_     sy4geu_mt_)~         hemm
symmrunt=: ((mp_mt_     sy4geu_mt_)~ |:    ) hemm
symmrunj=: ((mp_mt_     sy4geu_mt_)~ +     ) hemm
symmrunc=: ((mp_mt_     sy4geu_mt_)~ ct_mt_) hemm

symmrujn=:  (mp_mt_   +@sy4geu_mt_)~         hemm
symmrujt=: ((mp_mt_   +@sy4geu_mt_)~ |:    ) hemm
symmrujj=:  (mp_mt_ &:+ sy4geu_mt_)~         hemm
symmrujc=: ((mp_mt_   +@sy4geu_mt_)~ ct_mt_) hemm

hemmllnn=:  (mp_mt_~    he4gel_mt_)~         hemm
hemmllnt=: ((mp_mt_~    he4gel_mt_)~ |:    ) hemm
hemmllnj=: ((mp_mt_~    he4gel_mt_)~ +     ) hemm
hemmllnc=: ((mp_mt_~    he4gel_mt_)~ ct_mt_) hemm

hemmlljn=:  (mp_mt_~  +@he4gel_mt_)~         hemm
hemmlljt=: ((mp_mt_~  +@he4gel_mt_)~ |:    ) hemm
hemmlljj=:  (mp_mt_~&:+ he4gel_mt_)~         hemm
hemmlljc=: ((mp_mt_~  +@he4gel_mt_)~ ct_mt_) hemm

hemmlunn=:  (mp_mt_~    he4geu_mt_)~         hemm
hemmlunt=: ((mp_mt_~    he4geu_mt_)~ |:    ) hemm
hemmlunj=: ((mp_mt_~    he4geu_mt_)~ +     ) hemm
hemmlunc=: ((mp_mt_~    he4geu_mt_)~ ct_mt_) hemm

hemmlujn=:  (mp_mt_~  +@he4geu_mt_)~         hemm
hemmlujt=: ((mp_mt_~  +@he4geu_mt_)~ |:    ) hemm
hemmlujj=:  (mp_mt_~&:+ he4geu_mt_)~         hemm
hemmlujc=: ((mp_mt_~  +@he4geu_mt_)~ ct_mt_) hemm

hemmrlnn=:  (mp_mt_     he4gel_mt_)~         hemm
hemmrlnt=: ((mp_mt_     he4gel_mt_)~ |:    ) hemm
hemmrlnj=: ((mp_mt_     he4gel_mt_)~ +     ) hemm
hemmrlnc=: ((mp_mt_     he4gel_mt_)~ ct_mt_) hemm

hemmrljn=:  (mp_mt_   +@he4gel_mt_)~         hemm
hemmrljt=: ((mp_mt_   +@he4gel_mt_)~ |:    ) hemm
hemmrljj=:  (mp_mt_ &:+ he4gel_mt_)~         hemm
hemmrljc=: ((mp_mt_   +@he4gel_mt_)~ ct_mt_) hemm

hemmrunn=:  (mp_mt_     he4geu_mt_)~         hemm
hemmrunt=: ((mp_mt_     he4geu_mt_)~ |:    ) hemm
hemmrunj=: ((mp_mt_     he4geu_mt_)~ +     ) hemm
hemmrunc=: ((mp_mt_     he4geu_mt_)~ ct_mt_) hemm

hemmrujn=:  (mp_mt_   +@he4geu_mt_)~         hemm
hemmrujt=: ((mp_mt_   +@he4geu_mt_)~ |:    ) hemm
hemmrujj=:  (mp_mt_ &:+ he4geu_mt_)~         hemm
hemmrujc=: ((mp_mt_   +@he4geu_mt_)~ ct_mt_) hemm

NB. ---------------------------------------------------------
NB. trmm
NB.
NB. Description:
NB.   Adv. to make monad to perform the matrix-matrix
NB.   operation:
NB.     B := alpha * op(A) * B  (1)
NB.   or
NB.     B := alpha * B * op(A)  (2)
NB. where
NB.   A     - triangular
NB.   op(A) - either A, A^T, conj(A) or A^H
NB.
NB. Syntax:
NB.   Bupd=. ((mul trans@ref) trmm) alpha ; AA ; B
NB. where
NB.   ref   - monad to restore A from triangular part, is one
NB.           of:
NB.             trlpick   NB.  LT, A is L
NB.             trl1pick  NB. SLT, A is L1
NB.             trupick   NB.  UT, A is U
NB.             tru1pick  NB. SUT, A is U1
NB.   trans - monad to define the form of op(A), is one of:
NB.             ]         NB. op(A) :=      A
NB.             |:        NB. op(A) :=      A^T
NB.             +         NB. op(A) := conj(A)
NB.             ct        NB. op(A) :=      A^H
NB.   mul   - dyad to compute the product either (op(A) * B) or
NB.           (B * op(A)), is called as:
NB.             product=. B mul opA
NB.           and is one of:
NB.             mp~       NB. to compute (1)
NB.             mp        NB. to compute (2)

trmm=: 1 : '(0&{:: * 2&{::) u 1&{::'

NB. ---------------------------------------------------------
NB. Monad       Side    A     Reads in A    op(A)
NB. trmmllnn    (1)     L      LT                A
NB. trmmllnu    (1)     L1    SLT                A
NB. trmmlltn    (1)     L      LT                A^T
NB. trmmlltu    (1)     L1    SLT                A^T
NB. trmmlljn    (1)     L      LT           conj(A)
NB. trmmllju    (1)     L1    SLT           conj(A)
NB. trmmllcn    (1)     L      LT                A^H
NB. trmmllcu    (1)     L1    SLT                A^H
NB. trmmlunn    (1)     U      UT                A
NB. trmmlunu    (1)     U1    SUT                A
NB. trmmlutn    (1)     U      UT                A^T
NB. trmmlutu    (1)     U1    SUT                A^T
NB. trmmlujn    (1)     U      UT           conj(A)
NB. trmmluju    (1)     U1    SUT           conj(A)
NB. trmmlucn    (1)     U      UT                A^H
NB. trmmlucu    (1)     U1    SUT                A^H
NB. trmmrlnn    (2)     L      LT                A
NB. trmmrlnu    (2)     L1    SLT                A
NB. trmmrltn    (2)     L      LT                A^T
NB. trmmrltu    (2)     L1    SLT                A^T
NB. trmmrljn    (2)     L      LT           conj(A)
NB. trmmrlju    (2)     L1    SLT           conj(A)
NB. trmmrlcn    (2)     L      LT                A^H
NB. trmmrlcu    (2)     L1    SLT                A^H
NB. trmmrunn    (2)     U      UT                A
NB. trmmrunu    (2)     U1    SUT                A
NB. trmmrutn    (2)     U      UT                A^T
NB. trmmrutu    (2)     U1    SUT                A^T
NB. trmmrujn    (2)     U      UT           conj(A)
NB. trmmruju    (2)     U1    SUT           conj(A)
NB. trmmrucn    (2)     U      UT                A^H
NB. trmmrucu    (2)     U1    SUT                A^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     B := alpha * op(A) * B  (1)
NB.   or
NB.     B := alpha * B * op(A)  (2)
NB. where
NB.   A     - triangular
NB.   op(A) - either A, A^T, conj(A) or A^H
NB.
NB. Syntax:
NB.   Bupd=. trmmxxxx alpha ; AA ; B
NB. where
NB.   alpha - scalar
NB.   AA    - k×k-matrix, contains either non-zero or both
NB.           part(s) of A
NB.   B     - m×n-matrix
NB.   Bupd  - an updated B
NB.   A     - k×k-matrix, triangular
NB.   m     ≥ 0, the number of rows in B and Bupd
NB.   n     ≥ 0, the number of columns in B and Bupd
NB.   k     = m for trmmlxxx or k = n for trmmrxxx
NB.
NB. Notes:
NB. - monad       models BLAS'
NB.   trmmllnn    xTRMM('L','L','N','N',...)
NB.   trmmllnu    xTRMM('L','L','N','U',...)
NB.   trmmlltn    xTRMM('L','L','T','N',...)
NB.   trmmlltu    xTRMM('L','L','T','U',...)
NB.   trmmllcn    xTRMM('L','L','C','N',...)
NB.   trmmllcu    xTRMM('L','L','C','U',...)
NB.   trmmlunn    xTRMM('L','U','N','N',...)
NB.   trmmlunu    xTRMM('L','U','N','U',...)
NB.   trmmlutn    xTRMM('L','U','T','N',...)
NB.   trmmlutu    xTRMM('L','U','T','U',...)
NB.   trmmlucn    xTRMM('L','U','C','N',...)
NB.   trmmlucu    xTRMM('L','U','C','U',...)
NB.   trmmrlnn    xTRMM('R','L','N','N',...)
NB.   trmmrlnu    xTRMM('R','L','N','U',...)
NB.   trmmrltn    xTRMM('R','L','T','N',...)
NB.   trmmrltu    xTRMM('R','L','T','U',...)
NB.   trmmrlcn    xTRMM('R','L','C','N',...)
NB.   trmmrlcu    xTRMM('R','L','C','U',...)
NB.   trmmrunn    xTRMM('R','U','N','N',...)
NB.   trmmrunu    xTRMM('R','U','N','U',...)
NB.   trmmrutn    xTRMM('R','U','T','N',...)
NB.   trmmrutu    xTRMM('R','U','T','U',...)
NB.   trmmrucn    xTRMM('R','U','C','N',...)
NB.   trmmrucu    xTRMM('R','U','C','U',...)
NB. - reference implementation

trmmllnn=: (mp_mt_~        trlpick_mt_ ) trmm
trmmllnu=: (mp_mt_~        trl1pick_mt_) trmm
trmmlltn=: (mp_mt_~ |:    @trlpick_mt_ ) trmm
trmmlltu=: (mp_mt_~ |:    @trl1pick_mt_) trmm
trmmlljn=: (mp_mt_~ +     @trlpick_mt_ ) trmm
trmmllju=: (mp_mt_~ +     @trl1pick_mt_) trmm
trmmllcn=: (mp_mt_~ ct_mt_@trlpick_mt_ ) trmm
trmmllcu=: (mp_mt_~ ct_mt_@trl1pick_mt_) trmm

trmmlunn=: (mp_mt_~        trupick_mt_ ) trmm
trmmlunu=: (mp_mt_~        tru1pick_mt_) trmm
trmmlutn=: (mp_mt_~ |:    @trupick_mt_ ) trmm
trmmlutu=: (mp_mt_~ |:    @tru1pick_mt_) trmm
trmmlujn=: (mp_mt_~ +     @trupick_mt_ ) trmm
trmmluju=: (mp_mt_~ +     @tru1pick_mt_) trmm
trmmlucn=: (mp_mt_~ ct_mt_@trupick_mt_ ) trmm
trmmlucu=: (mp_mt_~ ct_mt_@tru1pick_mt_) trmm

trmmrlnn=: (mp_mt_         trlpick_mt_ ) trmm
trmmrlnu=: (mp_mt_         trl1pick_mt_) trmm
trmmrltn=: (mp_mt_  |:    @trlpick_mt_ ) trmm
trmmrltu=: (mp_mt_  |:    @trl1pick_mt_) trmm
trmmrljn=: (mp_mt_  +     @trlpick_mt_ ) trmm
trmmrlju=: (mp_mt_  +     @trl1pick_mt_) trmm
trmmrlcn=: (mp_mt_  ct_mt_@trlpick_mt_ ) trmm
trmmrlcu=: (mp_mt_  ct_mt_@trl1pick_mt_) trmm

trmmrunn=: (mp_mt_         trupick_mt_ ) trmm
trmmrunu=: (mp_mt_         tru1pick_mt_) trmm
trmmrutn=: (mp_mt_  |:    @trupick_mt_ ) trmm
trmmrutu=: (mp_mt_  |:    @tru1pick_mt_) trmm
trmmrujn=: (mp_mt_  +     @trupick_mt_ ) trmm
trmmruju=: (mp_mt_  +     @tru1pick_mt_) trmm
trmmrucn=: (mp_mt_  ct_mt_@trupick_mt_ ) trmm
trmmrucu=: (mp_mt_  ct_mt_@tru1pick_mt_) trmm

NB. ---------------------------------------------------------
NB. trmm3
NB.
NB. Description:
NB.   Adv. to make monad to perform the matrix-matrix
NB.   operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C  (1)
NB.   or
NB.     C := alpha * op2(B) * op1(A) + beta * C  (2)
NB. where
NB.   A      - triangular
NB.   opX(M) - either M, M^T, conj(M) or M^H
NB.
NB. Syntax:
NB.   Cupd=. (mul trmm3) alpha ; AA ; B ; beta ; C
NB. where
NB.   mul - dyad to restore A from AA triangular part, and to
NB.         compute the product either (op1(A) * op2(B)) or
NB.         (op2(B) * op1(A)), is called as:
NB.           product=. AA mul B

trmm3=: gemm

NB. ---------------------------------------------------------
NB. Monad        Side    A     Reads in A    op1(A)      op2(B)
NB. trmmllnnn    (1)     L      LT                A           B
NB. trmmllnnt    (1)     L      LT                A           B^T
NB. trmmllnnj    (1)     L      LT                A      conj(B)
NB. trmmllnnc    (1)     L      LT                A           B^H
NB. trmmllnun    (1)     L1    SLT                A           B
NB. trmmllnut    (1)     L1    SLT                A           B^T
NB. trmmllnuj    (1)     L1    SLT                A      conj(B)
NB. trmmllnuc    (1)     L1    SLT                A           B^H
NB. trmmlltnn    (1)     L      LT                A^T         B
NB. trmmlltnt    (1)     L      LT                A^T         B^T
NB. trmmlltnj    (1)     L      LT                A^T    conj(B)
NB. trmmlltnc    (1)     L      LT                A^T         B^H
NB. trmmlltun    (1)     L1    SLT                A^T         B
NB. trmmlltut    (1)     L1    SLT                A^T         B^T
NB. trmmlltuj    (1)     L1    SLT                A^T    conj(B)
NB. trmmlltuc    (1)     L1    SLT                A^T         B^H
NB. trmmlljnn    (1)     L      LT           conj(A)          B
NB. trmmlljnt    (1)     L      LT           conj(A)          B^T
NB. trmmlljnj    (1)     L      LT           conj(A)     conj(B)
NB. trmmlljnc    (1)     L      LT           conj(A)          B^H
NB. trmmlljun    (1)     L1    SLT           conj(A)          B
NB. trmmlljut    (1)     L1    SLT           conj(A)          B^T
NB. trmmlljuj    (1)     L1    SLT           conj(A)     conj(B)
NB. trmmlljuc    (1)     L1    SLT           conj(A)          B^H
NB. trmmllcnn    (1)     L      LT                A^H         B
NB. trmmllcnt    (1)     L      LT                A^H         B^T
NB. trmmllcnj    (1)     L      LT                A^H    conj(B)
NB. trmmllcnc    (1)     L      LT                A^H         B^H
NB. trmmllcun    (1)     L1    SLT                A^H         B
NB. trmmllcut    (1)     L1    SLT                A^H         B^T
NB. trmmllcuj    (1)     L1    SLT                A^H    conj(B)
NB. trmmllcuc    (1)     L1    SLT                A^H         B^H
NB. trmmlunnn    (1)     U      UT                A           B
NB. trmmlunnt    (1)     U      UT                A           B^T
NB. trmmlunnj    (1)     U      UT                A      conj(B)
NB. trmmlunnc    (1)     U      UT                A           B^H
NB. trmmlunun    (1)     U1    SUT                A           B
NB. trmmlunut    (1)     U1    SUT                A           B^T
NB. trmmlunuj    (1)     U1    SUT                A      conj(B)
NB. trmmlunuc    (1)     U1    SUT                A           B^H
NB. trmmlutnn    (1)     U      UT                A^T         B
NB. trmmlutnt    (1)     U      UT                A^T         B^T
NB. trmmlutnj    (1)     U      UT                A^T    conj(B)
NB. trmmlutnc    (1)     U      UT                A^T         B^H
NB. trmmlutun    (1)     U1    SUT                A^T         B
NB. trmmlutut    (1)     U1    SUT                A^T         B^T
NB. trmmlutuj    (1)     U1    SUT                A^T    conj(B)
NB. trmmlutuc    (1)     U1    SUT                A^T         B^H
NB. trmmlujnn    (1)     U      UT           conj(A)          B
NB. trmmlujnt    (1)     U      UT           conj(A)          B^T
NB. trmmlujnj    (1)     U      UT           conj(A)     conj(B)
NB. trmmlujnc    (1)     U      UT           conj(A)          B^H
NB. trmmlujun    (1)     U1    SUT           conj(A)          B
NB. trmmlujut    (1)     U1    SUT           conj(A)          B^T
NB. trmmlujuj    (1)     U1    SUT           conj(A)     conj(B)
NB. trmmlujuc    (1)     U1    SUT           conj(A)          B^H
NB. trmmlucnn    (1)     U      UT                A^H         B
NB. trmmlucnt    (1)     U      UT                A^H         B^T
NB. trmmlucnj    (1)     U      UT                A^H    conj(B)
NB. trmmlucnc    (1)     U      UT                A^H         B^H
NB. trmmlucun    (1)     U1    SUT                A^H         B
NB. trmmlucut    (1)     U1    SUT                A^H         B^T
NB. trmmlucuj    (1)     U1    SUT                A^H    conj(B)
NB. trmmlucuc    (1)     U1    SUT                A^H         B^H
NB. trmmrlnnn    (2)     L      LT                A           B
NB. trmmrlnnt    (2)     L      LT                A           B^T
NB. trmmrlnnj    (2)     L      LT                A      conj(B)
NB. trmmrlnnc    (2)     L      LT                A           B^H
NB. trmmrlnun    (2)     L1    SLT                A           B
NB. trmmrlnut    (2)     L1    SLT                A           B^T
NB. trmmrlnuj    (2)     L1    SLT                A      conj(B)
NB. trmmrlnuc    (2)     L1    SLT                A           B^H
NB. trmmrltnn    (2)     L      LT                A^T         B
NB. trmmrltnt    (2)     L      LT                A^T         B^T
NB. trmmrltnj    (2)     L      LT                A^T    conj(B)
NB. trmmrltnc    (2)     L      LT                A^T         B^H
NB. trmmrltun    (2)     L1    SLT                A^T         B
NB. trmmrltut    (2)     L1    SLT                A^T         B^T
NB. trmmrltuj    (2)     L1    SLT                A^T    conj(B)
NB. trmmrltuc    (2)     L1    SLT                A^T         B^H
NB. trmmrljnn    (2)     L      LT           conj(A)          B
NB. trmmrljnt    (2)     L      LT           conj(A)          B^T
NB. trmmrljnj    (2)     L      LT           conj(A)     conj(B)
NB. trmmrljnc    (2)     L      LT           conj(A)          B^H
NB. trmmrljun    (2)     L1    SLT           conj(A)          B
NB. trmmrljut    (2)     L1    SLT           conj(A)          B^T
NB. trmmrljuj    (2)     L1    SLT           conj(A)     conj(B)
NB. trmmrljuc    (2)     L1    SLT           conj(A)          B^H
NB. trmmrlcnn    (2)     L      LT                A^H         B
NB. trmmrlcnt    (2)     L      LT                A^H         B^T
NB. trmmrlcnj    (2)     L      LT                A^H    conj(B)
NB. trmmrlcnc    (2)     L      LT                A^H         B^H
NB. trmmrlcun    (2)     L1    SLT                A^H         B
NB. trmmrlcut    (2)     L1    SLT                A^H         B^T
NB. trmmrlcuj    (2)     L1    SLT                A^H    conj(B)
NB. trmmrlcuc    (2)     L1    SLT                A^H         B^H
NB. trmmrunnn    (2)     U     UT                 A           B
NB. trmmrunnt    (2)     U     UT                 A           B^T
NB. trmmrunnj    (2)     U     UT                 A      conj(B)
NB. trmmrunnc    (2)     U     UT                 A           B^H
NB. trmmrunun    (2)     U1    SUT                A           B
NB. trmmrunut    (2)     U1    SUT                A           B^T
NB. trmmrunuj    (2)     U1    SUT                A      conj(B)
NB. trmmrunuc    (2)     U1    SUT                A           B^H
NB. trmmrutnn    (2)     U      UT                A^T         B
NB. trmmrutnt    (2)     U      UT                A^T         B^T
NB. trmmrutnj    (2)     U      UT                A^T    conj(B)
NB. trmmrutnc    (2)     U      UT                A^T         B^H
NB. trmmrutun    (2)     U1    SUT                A^T         B
NB. trmmrutut    (2)     U1    SUT                A^T         B^T
NB. trmmrutuj    (2)     U1    SUT                A^T    conj(B)
NB. trmmrutuc    (2)     U1    SUT                A^T         B^H
NB. trmmrujnn    (2)     U      UT           conj(A)          B
NB. trmmrujnt    (2)     U      UT           conj(A)          B^T
NB. trmmrujnj    (2)     U      UT           conj(A)     conj(B)
NB. trmmrujnc    (2)     U      UT           conj(A)          B^H
NB. trmmrujun    (2)     U1    SUT           conj(A)          B
NB. trmmrujut    (2)     U1    SUT           conj(A)          B^T
NB. trmmrujuj    (2)     U1    SUT           conj(A)     conj(B)
NB. trmmrujuc    (2)     U1    SUT           conj(A)          B^H
NB. trmmrucnn    (2)     U      UT                A^H         B
NB. trmmrucnt    (2)     U      UT                A^H         B^T
NB. trmmrucnj    (2)     U      UT                A^H    conj(B)
NB. trmmrucnc    (2)     U      UT                A^H         B^H
NB. trmmrucun    (2)     U1    SUT                A^H         B
NB. trmmrucut    (2)     U1    SUT                A^H         B^T
NB. trmmrucuj    (2)     U1    SUT                A^H    conj(B)
NB. trmmrucuc    (2)     U1    SUT                A^H         B^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * op1(A) * op2(B) + beta * C  (1)
NB.   or
NB.     C := alpha * op2(B) * op1(A) + beta * C  (2)
NB. where
NB.   A      - triangular
NB.   opX(M) - either M, M^T, conj(M) or M^H
NB.
NB. Syntax:
NB.   Cupd=. trmmxxxxx alpha ; AA ; B ; beta ; C
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
NB. - monad        models BLIS'
NB.   trmmxxxxx    bli_xtrmm3 (...)

trmmllnnn=:  (mp_mt_~         trlpick_mt_ )~         trmm3
trmmllnnt=: ((mp_mt_~         trlpick_mt_ )~ |:    ) trmm3
trmmllnnj=: ((mp_mt_~         trlpick_mt_ )~ +     ) trmm3
trmmllnnc=: ((mp_mt_~         trlpick_mt_ )~ ct_mt_) trmm3

trmmllnun=:  (mp_mt_~         trl1pick_mt_)~         trmm3
trmmllnut=: ((mp_mt_~         trl1pick_mt_)~ |:    ) trmm3
trmmllnuj=: ((mp_mt_~         trl1pick_mt_)~ +     ) trmm3
trmmllnuc=: ((mp_mt_~         trl1pick_mt_)~ ct_mt_) trmm3

trmmlltnn=:  (mp_mt_~  |:    @trlpick_mt_ )~         trmm3
trmmlltnt=:  (mp_mt_~& |:     trlpick_mt_ )~         trmm3
trmmlltnj=: ((mp_mt_~  |:    @trlpick_mt_ )~ +     ) trmm3
trmmlltnc=: ((mp_mt_~  |:    @trlpick_mt_ )~ ct_mt_) trmm3

trmmlltun=:  (mp_mt_~  |:    @trl1pick_mt_)~         trmm3
trmmlltut=:  (mp_mt_~& |:     trl1pick_mt_)~         trmm3
trmmlltuj=: ((mp_mt_~  |:    @trl1pick_mt_)~ +     ) trmm3
trmmlltuc=: ((mp_mt_~  |:    @trl1pick_mt_)~ ct_mt_) trmm3

trmmlljnn=:  (mp_mt_~  +     @trlpick_mt_ )~         trmm3
trmmlljnt=: ((mp_mt_~  +     @trlpick_mt_ )~ |:    ) trmm3
trmmlljnj=:  (mp_mt_~&:+      trlpick_mt_ )~         trmm3
trmmlljnc=: ((mp_mt_~  +     @trlpick_mt_ )~ ct_mt_) trmm3

trmmlljun=:  (mp_mt_~  +     @trl1pick_mt_)~         trmm3
trmmlljut=: ((mp_mt_~  +     @trl1pick_mt_)~ |:    ) trmm3
trmmlljuj=:  (mp_mt_~&:+      trl1pick_mt_)~         trmm3
trmmlljuc=: ((mp_mt_~  +     @trl1pick_mt_)~ ct_mt_) trmm3

trmmllcnn=:  (mp_mt_~  ct_mt_@trlpick_mt_ )~         trmm3
trmmllcnt=: ((mp_mt_~  ct_mt_@trlpick_mt_ )~ |:    ) trmm3
trmmllcnj=: ((mp_mt_~  ct_mt_@trlpick_mt_ )~ +     ) trmm3
trmmllcnc=:  (mp_mt_~& ct_mt_ trlpick_mt_ )~         trmm3

trmmllcun=:  (mp_mt_~  ct_mt_@trl1pick_mt_)~         trmm3
trmmllcut=: ((mp_mt_~  ct_mt_@trl1pick_mt_)~ |:    ) trmm3
trmmllcuj=: ((mp_mt_~  ct_mt_@trl1pick_mt_)~ +     ) trmm3
trmmllcuc=:  (mp_mt_~& ct_mt_ trl1pick_mt_)~         trmm3

trmmlunnn=:  (mp_mt_~         trupick_mt_ )~         trmm3
trmmlunnt=: ((mp_mt_~         trupick_mt_ )~ |:    ) trmm3
trmmlunnj=: ((mp_mt_~         trupick_mt_ )~ +     ) trmm3
trmmlunnc=: ((mp_mt_~         trupick_mt_ )~ ct_mt_) trmm3

trmmlunun=:  (mp_mt_~         tru1pick_mt_)~         trmm3
trmmlunut=: ((mp_mt_~         tru1pick_mt_)~ |:    ) trmm3
trmmlunuj=: ((mp_mt_~         tru1pick_mt_)~ +     ) trmm3
trmmlunuc=: ((mp_mt_~         tru1pick_mt_)~ ct_mt_) trmm3

trmmlutnn=:  (mp_mt_~  |:    @trupick_mt_ )~         trmm3
trmmlutnt=:  (mp_mt_~& |:     trupick_mt_ )~         trmm3
trmmlutnj=: ((mp_mt_~  |:    @trupick_mt_ )~ +     ) trmm3
trmmlutnc=: ((mp_mt_~  |:    @trupick_mt_ )~ ct_mt_) trmm3

trmmlutun=:  (mp_mt_~  |:    @tru1pick_mt_)~         trmm3
trmmlutut=:  (mp_mt_~& |:     tru1pick_mt_)~         trmm3
trmmlutuj=: ((mp_mt_~  |:    @tru1pick_mt_)~ +     ) trmm3
trmmlutuc=: ((mp_mt_~  |:    @tru1pick_mt_)~ ct_mt_) trmm3

trmmlujnn=:  (mp_mt_~  +     @trupick_mt_ )~         trmm3
trmmlujnt=: ((mp_mt_~  +     @trupick_mt_ )~ |:    ) trmm3
trmmlujnj=:  (mp_mt_~&:+      trupick_mt_ )~         trmm3
trmmlujnc=: ((mp_mt_~  +     @trupick_mt_ )~ ct_mt_) trmm3

trmmlujun=:  (mp_mt_~  +     @tru1pick_mt_)~         trmm3
trmmlujut=: ((mp_mt_~  +     @tru1pick_mt_)~ |:    ) trmm3
trmmlujuj=:  (mp_mt_~&:+      tru1pick_mt_)~         trmm3
trmmlujuc=: ((mp_mt_~  +     @tru1pick_mt_)~ ct_mt_) trmm3

trmmlucnn=:  (mp_mt_~  ct_mt_@trupick_mt_ )~         trmm3
trmmlucnt=: ((mp_mt_~  ct_mt_@trupick_mt_ )~ |:    ) trmm3
trmmlucnj=: ((mp_mt_~  ct_mt_@trupick_mt_ )~ +     ) trmm3
trmmlucnc=:  (mp_mt_~& ct_mt_ trupick_mt_ )~         trmm3

trmmlucun=:  (mp_mt_~  ct_mt_@tru1pick_mt_)~         trmm3
trmmlucut=: ((mp_mt_~  ct_mt_@tru1pick_mt_)~ |:    ) trmm3
trmmlucuj=: ((mp_mt_~  ct_mt_@tru1pick_mt_)~ +     ) trmm3
trmmlucuc=:  (mp_mt_~& ct_mt_ tru1pick_mt_)~         trmm3

trmmrlnnn=:  (mp_mt_          trlpick_mt_ )~         trmm3
trmmrlnnt=: ((mp_mt_          trlpick_mt_ )~ |:    ) trmm3
trmmrlnnj=: ((mp_mt_          trlpick_mt_ )~ +     ) trmm3
trmmrlnnc=: ((mp_mt_          trlpick_mt_ )~ ct_mt_) trmm3

trmmrlnun=:  (mp_mt_          trl1pick_mt_)~         trmm3
trmmrlnut=: ((mp_mt_          trl1pick_mt_)~ |:    ) trmm3
trmmrlnuj=: ((mp_mt_          trl1pick_mt_)~ +     ) trmm3
trmmrlnuc=: ((mp_mt_          trl1pick_mt_)~ ct_mt_) trmm3

trmmrltnn=:  (mp_mt_   |:    @trlpick_mt_ )~         trmm3
trmmrltnt=:  (mp_mt_ & |:     trlpick_mt_ )~         trmm3
trmmrltnj=: ((mp_mt_   |:    @trlpick_mt_ )~ +     ) trmm3
trmmrltnc=: ((mp_mt_   |:    @trlpick_mt_ )~ ct_mt_) trmm3

trmmrltun=:  (mp_mt_   |:    @trl1pick_mt_)~         trmm3
trmmrltut=:  (mp_mt_ & |:     trl1pick_mt_)~         trmm3
trmmrltuj=: ((mp_mt_   |:    @trl1pick_mt_)~ +     ) trmm3
trmmrltuc=: ((mp_mt_   |:    @trl1pick_mt_)~ ct_mt_) trmm3

trmmrljnn=:  (mp_mt_   +     @trlpick_mt_ )~         trmm3
trmmrljnt=: ((mp_mt_   +     @trlpick_mt_ )~ |:    ) trmm3
trmmrljnj=:  (mp_mt_ &:+      trlpick_mt_ )~         trmm3
trmmrljnc=: ((mp_mt_   +     @trlpick_mt_ )~ ct_mt_) trmm3

trmmrljun=:  (mp_mt_   +     @trl1pick_mt_)~         trmm3
trmmrljut=: ((mp_mt_   +     @trl1pick_mt_)~ |:    ) trmm3
trmmrljuj=:  (mp_mt_ &:+      trl1pick_mt_)~         trmm3
trmmrljuc=: ((mp_mt_   +     @trl1pick_mt_)~ ct_mt_) trmm3

trmmrlcnn=:  (mp_mt_   ct_mt_@trlpick_mt_ )~         trmm3
trmmrlcnt=: ((mp_mt_   ct_mt_@trlpick_mt_ )~ |:    ) trmm3
trmmrlcnj=: ((mp_mt_   ct_mt_@trlpick_mt_ )~ +     ) trmm3
trmmrlcnc=:  (mp_mt_ & ct_mt_ trlpick_mt_ )~         trmm3

trmmrlcun=:  (mp_mt_   ct_mt_@trl1pick_mt_)~         trmm3
trmmrlcut=: ((mp_mt_   ct_mt_@trl1pick_mt_)~ |:    ) trmm3
trmmrlcuj=: ((mp_mt_   ct_mt_@trl1pick_mt_)~ +     ) trmm3
trmmrlcuc=:  (mp_mt_ & ct_mt_ trl1pick_mt_)~         trmm3

trmmrunnn=:  (mp_mt_          trupick_mt_ )~         trmm3
trmmrunnt=: ((mp_mt_          trupick_mt_ )~ |:    ) trmm3
trmmrunnj=: ((mp_mt_          trupick_mt_ )~ +     ) trmm3
trmmrunnc=: ((mp_mt_          trupick_mt_ )~ ct_mt_) trmm3

trmmrunun=:  (mp_mt_          tru1pick_mt_)~         trmm3
trmmrunut=: ((mp_mt_          tru1pick_mt_)~ |:    ) trmm3
trmmrunuj=: ((mp_mt_          tru1pick_mt_)~ +     ) trmm3
trmmrunuc=: ((mp_mt_          tru1pick_mt_)~ ct_mt_) trmm3

trmmrutnn=:  (mp_mt_   |:    @trupick_mt_ )~         trmm3
trmmrutnt=:  (mp_mt_ & |:     trupick_mt_ )~         trmm3
trmmrutnj=: ((mp_mt_   |:    @trupick_mt_ )~ +     ) trmm3
trmmrutnc=: ((mp_mt_   |:    @trupick_mt_ )~ ct_mt_) trmm3

trmmrutun=:  (mp_mt_   |:    @tru1pick_mt_)~         trmm3
trmmrutut=:  (mp_mt_ & |:     tru1pick_mt_)~         trmm3
trmmrutuj=: ((mp_mt_   |:    @tru1pick_mt_)~ +     ) trmm3
trmmrutuc=: ((mp_mt_   |:    @tru1pick_mt_)~ ct_mt_) trmm3

trmmrujnn=:  (mp_mt_   +     @trupick_mt_ )~         trmm3
trmmrujnt=: ((mp_mt_   +     @trupick_mt_ )~ |:    ) trmm3
trmmrujnj=:  (mp_mt_ &:+      trupick_mt_ )~         trmm3
trmmrujnc=: ((mp_mt_   +     @trupick_mt_ )~ ct_mt_) trmm3

trmmrujun=:  (mp_mt_   +     @tru1pick_mt_)~         trmm3
trmmrujut=: ((mp_mt_   +     @tru1pick_mt_)~ |:    ) trmm3
trmmrujuj=:  (mp_mt_ &:+      tru1pick_mt_)~         trmm3
trmmrujuc=: ((mp_mt_   +     @tru1pick_mt_)~ ct_mt_) trmm3

trmmrucnn=:  (mp_mt_   ct_mt_@trupick_mt_ )~         trmm3
trmmrucnt=: ((mp_mt_   ct_mt_@trupick_mt_ )~ |:    ) trmm3
trmmrucnj=: ((mp_mt_   ct_mt_@trupick_mt_ )~ +     ) trmm3
trmmrucnc=:  (mp_mt_ & ct_mt_ trupick_mt_ )~         trmm3

trmmrucun=:  (mp_mt_   ct_mt_@tru1pick_mt_)~         trmm3
trmmrucut=: ((mp_mt_   ct_mt_@tru1pick_mt_)~ |:    ) trmm3
trmmrucuj=: ((mp_mt_   ct_mt_@tru1pick_mt_)~ +     ) trmm3
trmmrucuc=:  (mp_mt_ & ct_mt_ tru1pick_mt_)~         trmm3

NB. ---------------------------------------------------------
NB. trsv
NB.
NB. Description:
NB.   Conj. to make monad to solve the equation:
NB.     op(A) * x = b
NB. where
NB.   A - triangular
NB.
NB. Syntax:
NB.   x=. (sol trsv) AA ; b ; incb
NB. where
NB.   sol   - dyad to solve equation with triangular matrix,
NB.           is called as:
NB.             x=. b sol A
NB.           e.g.
NB.             (%. trans@ref)  NB. to exploit built-in Matrix Divide (%.)
NB.             trsmllnn~       NB. to exploit mt's (trsmllnn)
NB.   ref   - monad to restore A from triangular part, is one
NB.           of:
NB.             trlpick         NB.  LT, A is L
NB.             trl1pick        NB. SLT, A is L1
NB.             trupick         NB.  UT, A is U
NB.             tru1pick        NB. SUT, A is U1
NB.   trans - monad to define the form of op(A), is one of:
NB.             ]               NB. op(A) := A
NB.             |:              NB. op(A) := A^T
NB.             ct              NB. op(A) := A^H

trsv=: trmv

NB. ---------------------------------------------------------
NB. Monad      A     Reads in A    op(A)
NB. trsvlnn    L      LT           A
NB. trsvlnu    L1    SLT           A
NB. trsvltn    L      LT           A^T
NB. trsvltu    L1    SLT           A^T
NB. trsvlcn    L      LT           A^H
NB. trsvlcu    L1    SLT           A^H
NB. trsvunn    U      UT           A
NB. trsvunu    U1    SUT           A
NB. trsvutn    U      UT           A^T
NB. trsvutu    U1    SUT           A^T
NB. trsvucn    U      UT           A^H
NB. trsvucu    U1    SUT           A^H
NB.
NB. Description:
NB.   Solves the equation:
NB.     op(A) * x = b
NB. where
NB.   A - triangular
NB.
NB. Syntax:
NB.   x=. trsvxxx AA ; b ; incb
NB. where
NB.   AA   - n×n-matrix, contains either non-zero or both
NB.          part(s) of A
NB.   b    - (1+(n-1)*|incb|)-vector, the RHS
NB.   incb ≠ 0, the increment for the elements of b and x
NB.   x    - the same shape as b, the solution
NB.   A    - n×n-matrix, triangular
NB.   n    ≥ 0, the size of A
NB.
NB. Notes:
NB. - monad      models BLAS'
NB.   trsvlnn    xTRSV('L','N','N',...)
NB.   trsvlnu    xTRSV('L','N','U',...)
NB.   trsvltn    xTRSV('L','T','N',...)
NB.   trsvltu    xTRSV('L','T','U',...)
NB.   trsvlcn    xTRSV('L','C','N',...)
NB.   trsvlcu    xTRSV('L','C','U',...)
NB.   trsvunn    xTRSV('U','N','N',...)
NB.   trsvunu    xTRSV('U','N','U',...)
NB.   trsvutn    xTRSV('U','T','N',...)
NB.   trsvutu    xTRSV('U','T','U',...)
NB.   trsvucn    xTRSV('U','C','N',...)
NB.   trsvucu    xTRSV('U','C','U',...)
NB. - reference implementation

trsvlnn=: trsmllnn_mt_~ trsv
trsvlnu=: trsmllnu_mt_~ trsv
trsvltn=: trsmlltn_mt_~ trsv
trsvltu=: trsmlltu_mt_~ trsv
trsvlcn=: trsmllcn_mt_~ trsv
trsvlcu=: trsmllcu_mt_~ trsv
trsvunn=: trsmlunn_mt_~ trsv
trsvunu=: trsmlunu_mt_~ trsv
trsvutn=: trsmlutn_mt_~ trsv
trsvutu=: trsmlutu_mt_~ trsv
trsvucn=: trsmlucn_mt_~ trsv
trsvucu=: trsmlucu_mt_~ trsv

NB. ---------------------------------------------------------
NB. trsm
NB.
NB. Description:
NB.   Adv. to make ambivalent verb to solve the linear
NB.   monomial matrix equation:
NB.     op(A) * X = alpha * B
NB.   or
NB.     X * op(A) = alpha * B
NB. where
NB.   A - triangular
NB.
NB. Syntax:
NB.   X=.    (mdiv trsm) alpha ; AA ; B  (1)
NB.   X=. AA (mdiv trsm)              B  (2)
NB. where
NB.   mdiv - dyad, defines what equation is to be solved, is
NB.          called as:
NB.            X=. (alpha * B) mdiv A  NB. for (1)
NB.            X=.          B  mdiv A  NB. for (2)

trsm=: 1 : '(u trmm_mt_) : (u~)'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb        Side    A     Reads in A    op(A)
NB. trsmllnn    (1)     L      LT                A
NB. trsmllnu    (1)     L1    SLT                A
NB. trsmlltn    (1)     L      LT                A^T
NB. trsmlltu    (1)     L1    SLT                A^T
NB. trsmlljn    (1)     L      LT           conj(A)
NB. trsmllju    (1)     L1    SLT           conj(A)
NB. trsmllcn    (1)     L      LT                A^H
NB. trsmllcu    (1)     L1    SLT                A^H
NB. trsmlunn    (1)     U      UT                A
NB. trsmlunu    (1)     U1    SUT                A
NB. trsmlutn    (1)     U      UT                A^T
NB. trsmlutu    (1)     U1    SUT                A^T
NB. trsmlujn    (1)     U      UT           conj(A)
NB. trsmluju    (1)     U1    SUT           conj(A)
NB. trsmlucn    (1)     U      UT                A^H
NB. trsmlucu    (1)     U1    SUT                A^H
NB. trsmrlnn    (2)     L      LT                A
NB. trsmrlnu    (2)     L1    SLT                A
NB. trsmrltn    (2)     L      LT                A^T
NB. trsmrltu    (2)     L1    SLT                A^T
NB. trsmrljn    (2)     L      LT           conj(A)
NB. trsmrlju    (2)     L1    SLT           conj(A)
NB. trsmrlcn    (2)     L      LT                A^H
NB. trsmrlcu    (2)     L1    SLT                A^H
NB. trsmrunn    (2)     U      UT                A
NB. trsmrunu    (2)     U1    SUT                A
NB. trsmrutn    (2)     U      UT                A^T
NB. trsmrutu    (2)     U1    SUT                A^T
NB. trsmrujn    (2)     U      UT           conj(A)
NB. trsmruju    (2)     U1    SUT           conj(A)
NB. trsmrucn    (2)     U      UT                A^H
NB. trsmrucu    (2)     U1    SUT                A^H
NB.
NB. Description:
NB.   Ambivalent verb to solve the linear monomial matrix
NB.   equation:
NB.     op(A) * X = alpha * B  (1)
NB.   or
NB.     X * op(A) = alpha * B  (2)
NB. where
NB.   A     - triangular
NB.   op(A) - either A, A^T, conj(A) or A^H
NB.
NB. Syntax:
NB.   X=.    trsmxxxx alpha ; AA ; B
NB.   X=. AA trsmxxxx              B
NB. where
NB.   alpha - scalar, is supposed to be 1 in dyadic case
NB.   AA    - k×k-matrix, contains either non-zero or both
NB.           part(s) of A
NB.   B     - l-vector or m×n-matrix, the RHS
NB.   X     - the same shape as B, the solution(s)
NB.   A     - k×k-matrix, triangular
NB.   m     ≥ 0, the number of rows in B and X
NB.   n     ≥ 0, the number of columns in B and X
NB.   k     = m for trsmlxxx or k = n for trsmrxxx
NB.   l     = n for trsmlxxx or l = m for trsmrxxx
NB.
NB. Notes:
NB. - for 2-rank B and X:
NB.   - monadic case:
NB.     - verb       models BLAS'
NB.       trsmllnn   xTRSM('L','L','N','N',...)
NB.       trsmllnu   xTRSM('L','L','N','U',...)
NB.       trsmlltn   xTRSM('L','L','T','N',...)
NB.       trsmlltu   xTRSM('L','L','T','U',...)
NB.       trsmllcn   xTRSM('L','L','C','N',...)
NB.       trsmllcu   xTRSM('L','L','C','U',...)
NB.       trsmlunn   xTRSM('L','U','N','N',...)
NB.       trsmlunu   xTRSM('L','U','N','U',...)
NB.       trsmlutn   xTRSM('L','U','T','N',...)
NB.       trsmlutu   xTRSM('L','U','T','U',...)
NB.       trsmlucn   xTRSM('L','U','C','N',...)
NB.       trsmlucu   xTRSM('L','U','C','U',...)
NB.       trsmrlnn   xTRSM('R','L','N','N',...)
NB.       trsmrlnu   xTRSM('R','L','N','U',...)
NB.       trsmrltn   xTRSM('R','L','T','N',...)
NB.       trsmrltu   xTRSM('R','L','T','U',...)
NB.       trsmrlcn   xTRSM('R','L','C','N',...)
NB.       trsmrlcu   xTRSM('R','L','C','U',...)
NB.       trsmrunn   xTRSM('R','U','N','N',...)
NB.       trsmrunu   xTRSM('R','U','N','U',...)
NB.       trsmrutn   xTRSM('R','U','T','N',...)
NB.       trsmrutu   xTRSM('R','U','T','U',...)
NB.       trsmrucn   xTRSM('R','U','C','N',...)
NB.       trsmrucu   xTRSM('R','U','C','U',...)
NB.     - reference implementation
NB.   - dyadic case models BLAS' xTRSM with (alpha=1)
NB. - for 1-rank B and X:
NB.   - monadic case with (alpha=1) simulates BLAS' xTRSV
NB.     with (incb=1)
NB.   - dyadic case simulates BLAS' xTRSV with (incb=1)
NB. - (alpha=0) means (X -: 0"0 B)

trsmllnn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i.   # y do. z=. z ,      ((i {   y) - (    i       {. ai    ) mp  z) % i { ai=. i {   x end.   z}}
trsmlltn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i. - # y do. z=. z ,~     ((i {   y) - ((>: i)      }. ai    ) mp  z) % i { ai=. i {"1 x end.   z}}
trsmlljn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i.   # y do. z=. z ,      ((i {   y) - (    i       {. ai    ) mp  z) % i { ai=. i {   x end. + z}}
trsmllcn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i. - # y do. z=. z ,~     ((i {   y) - ((>: i)      }. ai    ) mp  z) % i { ai=. i {"1 x end. + z}}

trsmlunn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i. - # y do. z=. z ,~     ((i {   y) - ((>: i)      }. ai    ) mp  z) % i { ai=. i {   x end.   z}}
trsmlutn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i.   # y do. z=. z ,      ((i {   y) - (    i       {. ai    ) mp  z) % i { ai=. i {"1 x end.   z}}
trsmlujn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i. - # y do. z=. z ,~     ((i {   y) - ((>: i)      }. ai    ) mp  z) % i { ai=. i {   x end. + z}}
trsmlucn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i.   # y do. z=. z ,      ((i {   y) - (    i       {. ai    ) mp  z) % i { ai=. i {"1 x end. + z}}

trsmrlnn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i. - c y do. z=. z ,~"1 0 ((i {"1 y) - ((>: i)      }. ai    ) mp~ z) % i { ai=. i {"1 x end.   z}}
trsmrltn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i.   c y do. z=. z , "1 0 ((i {"1 y) - (    i       {. ai    ) mp~ z) % i { ai=. i {   x end.   z}}
trsmrljn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1 y=. + y for_i. i. - c y do. z=. z ,~"1 0 ((i {"1 y) - ((>: i)      }. ai    ) mp~ z) % i { ai=. i {"1 x end. + z}}
trsmrlcn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1 y=. + y for_i. i.   c y do. z=. z , "1 0 ((i {"1 y) - (    i       {. ai    ) mp~ z) % i { ai=. i {   x end. + z}}

trsmrunn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i.   c y do. z=. z , "1 0 ((i {"1 y) - (    i       {. ai    ) mp~ z) % i { ai=. i {"1 x end.   z}}
trsmrutn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i. - c y do. z=. z ,~"1 0 ((i {"1 y) - ((>: i)      }. ai    ) mp~ z) % i { ai=. i {   x end.   z}}
trsmrujn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1 y=. + y for_i. i.   c y do. z=. z , "1 0 ((i {"1 y) - (    i       {. ai    ) mp~ z) % i { ai=. i {"1 x end. + z}}
trsmrucn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1 y=. + y for_i. i. - c y do. z=. z ,~"1 0 ((i {"1 y) - ((>: i)      }. ai    ) mp~ z) % i { ai=. i {   x end. + z}}

trsmllnu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i.   # y do. z=. z ,       (i {   y) - (    i (   [ {. {  ) x) mp  z                     end.   z}}
trsmlltu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i. - # y do. z=. z ,~      (i {   y) - (    i (>:@[ }. {"1) x) mp  z                     end.   z}}
trsmllju=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i.   # y do. z=. z ,       (i {   y) - (    i (   [ {. {  ) x) mp  z                     end. + z}}
trsmllcu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i. - # y do. z=. z ,~      (i {   y) - (    i (>:@[ }. {"1) x) mp  z                     end. + z}}

trsmlunu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i. - # y do. z=. z ,~      (i {   y) - (    i (>:@[ }. {  ) x) mp  z                     end.   z}}
trsmlutu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i.   # y do. z=. z ,       (i {   y) - (    i (   [ {. {"1) x) mp  z                     end.   z}}
trsmluju=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i. - # y do. z=. z ,~      (i {   y) - (    i (>:@[ }. {  ) x) mp  z                     end. + z}}
trsmlucu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i.   # y do. z=. z ,       (i {   y) - (    i (   [ {. {"1) x) mp  z                     end. + z}}

trsmrlnu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i. - c y do. z=. z ,~"1 0  (i {"1 y) - (    i (>:@[ }. {"1) x) mp~ z                     end.   z}}
trsmrltu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i.   c y do. z=. z , "1 0  (i {"1 y) - (    i (   [ {. {  ) x) mp~ z                     end.   z}}
trsmrlju=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1 y=. + y for_i. i. - c y do. z=. z ,~"1 0  (i {"1 y) - (    i (>:@[ }. {"1) x) mp~ z                     end. + z}}
trsmrlcu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1 y=. + y for_i. i.   c y do. z=. z , "1 0  (i {"1 y) - (    i (   [ {. {  ) x) mp~ z                     end. + z}}

trsmrunu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i.   c y do. z=. z , "1 0  (i {"1 y) - (    i (   [ {. {"1) x) mp~ z                     end.   z}}
trsmrutu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i. - c y do. z=. z ,~"1 0  (i {"1 y) - (    i (>:@[ }. {  ) x) mp~ z                     end.   z}}
trsmruju=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1 y=. + y for_i. i.   c y do. z=. z , "1 0  (i {"1 y) - (    i (   [ {. {"1) x) mp~ z                     end. + z}}
trsmrucu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1 y=. + y for_i. i. - c y do. z=. z ,~"1 0  (i {"1 y) - (    i (>:@[ }. {  ) x) mp~ z                     end. + z}}

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testbasicger
NB.
NB. Description:
NB.   Test rank 1 operations:
NB.   - DGER ZGERC ZGERU (BLAS)
NB.   by vectors and general matrix
NB.
NB. Syntax:
NB.   log=. testbasicger x ; y ; A
NB. where
NB.   x   - m-vector
NB.   y   - n-vector
NB.   A   - m×n-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs

testbasicger=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  inc=. 1 2 _1 _2
  'x y A'=. y
  argsd=. { (<"0 dcoeff) ; (< x) ; (<"0 inc) ; (< y) ; (<"0 inc) ; < < A
  argsz=. { (<"0 zcoeff) ; (< x) ; (<"0 inc) ; (< y) ; (<"0 inc) ; < < A

  NB. for every i feed the tuple (alpha_i ; expanded_x_i ; incx_i ; expanded_y_i ; incy_i ; A) to tmonad
  log=.          ('dger_mtbla_'  tmonad (]`]`nan`nan`(geru chk4r)))@(3 expand 4)@(1 expand 2)@>"0 argsd
  log=. log lcat ('zgeru_mtbla_' tmonad (]`]`nan`nan`(geru chk4r)))@(3 expand 4)@(1 expand 2)@>"0 argsz
  log=. log lcat ('zgerc_mtbla_' tmonad (]`]`nan`nan`(gerc chk4r)))@(3 expand 4)@(1 expand 2)@>"0 argsz
)

NB. ---------------------------------------------------------
NB. testbasicher
NB.
NB. Description:
NB.   Test hermitian (symmetric) rank 1 operations:
NB.   - DSYR ZHER (BLAS)
NB.   by vector and Hermitian (symmetric) matrix
NB.
NB. Syntax:
NB.   log=. testbasicher x ; AA
NB. where
NB.   x   - n-vector
NB.   AA  - n×n-matrix, A material with real diagonal
NB.   log - 6-vector of boxes, test log, see test.ijs

testbasicher=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  inc=. 1 2 _1 _2
  'x AA'=. y
  args=. { (<"0 dcoeff) ; (< x) ; (<"0 inc) ; < < AA

  NB. for every i feed the tuple (alpha_i ; expanded_x_i ; incx_i ; AA) to tmonad
  log=.          ('dsyrl_mtbla_' tmonad (]`]`nan`nan`(syrl chk5r trlpick)))@(1 expand 2)@>"0 args
  log=. log lcat ('dsyru_mtbla_' tmonad (]`]`nan`nan`(syru chk5r trupick)))@(1 expand 2)@>"0 args
  log=. log lcat ('zherl_mtbla_' tmonad (]`]`nan`nan`(herl chk5r trlpick)))@(1 expand 2)@>"0 args
  log=. log lcat ('zheru_mtbla_' tmonad (]`]`nan`nan`(heru chk5r trupick)))@(1 expand 2)@>"0 args
)

NB. ---------------------------------------------------------
NB. testbasicr
NB.
NB. Description:
NB.   Adv. to make verb to test basic rank 1 operations
NB.
NB. Syntax:
NB.   log=. (mkmat testbasicr) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testbasicr_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testbasicr_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testbasicr_mt_ 150 200

testbasicr=: 1 : 'nolog_mt_`(testbasicher_mt_@(u@{. ; (9&o. upddiag_mt_)@u))@.(=/) ,&.>~ testbasicger_mt_@(u@{. ; u@{: ; u) [ require@''math/mt/external/blas/r'''

NB. ---------------------------------------------------------
NB. testbasicher2
NB.
NB. Description:
NB.   Test hermitian (symmetric) rank 2 operations:
NB.   - DSYR2 ZHER2 (BLAS)
NB.   by vectors and Hermitian (symmetric) matrix
NB.
NB. Syntax:
NB.   log=. testbasicher2 x ; y ; AA
NB. where
NB.   x   - n-vector
NB.   y   - n-vector
NB.   AA  - n×n-matrix, A material with real diagonal
NB.   log - 6-vector of boxes, test log, see test.ijs

testbasicher2=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  inc=. 1 2 _1 _2
  'x y AA'=. y
  argsd=. { (<"0 dcoeff) ; (< x) ; (<"0 inc) ; (< y) ; (<"0 inc) ; < < AA
  argsz=. { (<"0 zcoeff) ; (< x) ; (<"0 inc) ; (< y) ; (<"0 inc) ; < < AA

  NB. for every i feed the tuple (alpha_i ; expanded_x_i ; incx_i ; expanded_y_i ; incy_i ; AA) to tmonad
  log=.          ('dsyr2l_mtbla_' tmonad (]`]`nan`nan`(syr2l chk6r2 trlpick)))@(3 expand 4)@(1 expand 2)@>"0 argsd
  log=. log lcat ('dsyr2u_mtbla_' tmonad (]`]`nan`nan`(syr2u chk6r2 trupick)))@(3 expand 4)@(1 expand 2)@>"0 argsd
  log=. log lcat ('zher2l_mtbla_' tmonad (]`]`nan`nan`(her2l chk6r2 trlpick)))@(3 expand 4)@(1 expand 2)@>"0 argsz
  log=. log lcat ('zher2u_mtbla_' tmonad (]`]`nan`nan`(her2u chk6r2 trupick)))@(3 expand 4)@(1 expand 2)@>"0 argsz
)

NB. ---------------------------------------------------------
NB. testbasicr2
NB.
NB. Description:
NB.   Adv. to make verb to test basic rank 2 operations
NB.
NB. Syntax:
NB.   log=. (mkmat testbasicr2) (n,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (n,n)
NB.   (n,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testbasicr2_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testbasicr2_mt_ 200 200
NB. - test by random square complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testbasicr2_mt_ 250 250

testbasicr2=: 1 : 'testbasicher2_mt_@(u@{. ; u@{: ; (9&o. upddiag_mt_)@u) [ require@''math/mt/external/blas/r2'''

NB. ---------------------------------------------------------
NB. testbasicsyrk
NB.
NB. Description:
NB.   Test symmetric rank k operations:
NB.   - xSYRK (BLAS)
NB.   by general and symmetric matrices
NB.
NB. Syntax:
NB.   log=. testbasicsyrk A ; CC
NB. where
NB.   A   - m×n-matrix
NB.   CC  - k×k-matrix, C material
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)

testbasicsyrk=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  'A CC'=. y
  mn=. <./ 'm n'=. $ A
  CCmn=. (2 # mn) {. CC
  argsdl=. { (<"0 dcoeff) ; (< A) ; (<"0 dcoeff) ; < < CCmn [^:(m < n) CC
  argsdg=. { (<"0 dcoeff) ; (< A) ; (<"0 dcoeff) ; < < CCmn [^:(m > n) CC
  argszl=. { (<"0 zcoeff) ; (< A) ; (<"0 zcoeff) ; < < CCmn [^:(m < n) CC
  argszg=. { (<"0 zcoeff) ; (< A) ; (<"0 zcoeff) ; < < CCmn [^:(m > n) CC

  NB. for every i feed the tuple (alpha_i ; A ; beta_i ; CC) to tmonad
  log=.          ('dsyrkln_mtbla_' tmonad (]`]`nan`nan`(syrkln chk4rk trlpick)))@>"0 argsdl
  log=. log lcat ('dsyrklt_mtbla_' tmonad (]`]`nan`nan`(syrklt chk4rk trlpick)))@>"0 argsdg
  log=. log lcat ('dsyrkun_mtbla_' tmonad (]`]`nan`nan`(syrkun chk4rk trupick)))@>"0 argsdl
  log=. log lcat ('dsyrkut_mtbla_' tmonad (]`]`nan`nan`(syrkut chk4rk trupick)))@>"0 argsdg
  log=. log lcat ('zsyrkln_mtbla_' tmonad (]`]`nan`nan`(syrkln chk4rk trlpick)))@>"0 argszl
  log=. log lcat ('zsyrklt_mtbla_' tmonad (]`]`nan`nan`(syrklt chk4rk trlpick)))@>"0 argszg
  log=. log lcat ('zsyrkun_mtbla_' tmonad (]`]`nan`nan`(syrkun chk4rk trupick)))@>"0 argszl
  log=. log lcat ('zsyrkut_mtbla_' tmonad (]`]`nan`nan`(syrkut chk4rk trupick)))@>"0 argszg
)

NB. ---------------------------------------------------------
NB. testbasicherk
NB.
NB. Description:
NB.   Test the hermitian rank k operation:
NB.   - ZHERK (BLAS)
NB.   by general and Hermitian matrices
NB.
NB. Syntax:
NB.   log=. testbasicherk A ; CC
NB. where
NB.   A   - m×n-matrix
NB.   CC  - k×k-matrix, C material with real diagonal
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)

testbasicherk=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  'A CC'=. y
  mn=. <./ 'm n'=. $ A
  CCmn=. (2 # mn) {. CC
  argsl=. { (<"0 dcoeff) ; (< A) ; (<"0 dcoeff) ; < < CCmn [^:(m < n) CC
  argsg=. { (<"0 dcoeff) ; (< A) ; (<"0 dcoeff) ; < < CCmn [^:(m > n) CC

  NB. for every i feed the tuple (alpha_i ; A ; beta_i ; CC) to tmonad
  log=.          ('zherkln_mtbla_' tmonad (]`]`nan`nan`(herkln chk4rk trlpick)))@>"0 argsl
  log=. log lcat ('zherklc_mtbla_' tmonad (]`]`nan`nan`(herklc chk4rk trlpick)))@>"0 argsg
  log=. log lcat ('zherkun_mtbla_' tmonad (]`]`nan`nan`(herkun chk4rk trupick)))@>"0 argsl
  log=. log lcat ('zherkuc_mtbla_' tmonad (]`]`nan`nan`(herkuc chk4rk trupick)))@>"0 argsg
)

NB. ---------------------------------------------------------
NB. testbasicrk
NB.
NB. Description:
NB.   Adv. to make verb to test basic rank k operations
NB.
NB. Syntax:
NB.   log=. (mkmat testbasicrk) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testbasicrk_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testbasicrk_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testbasicrk_mt_ 150 200

testbasicrk=: 1 : 'testbasicherk_mt_@(u ; (9&o. upddiag_mt_)@u@(2 # >./)) ,&.>~ testbasicsyrk_mt_@(u ; u@(2 # >./)) [ require@''math/mt/external/blas/rk'''

NB. ---------------------------------------------------------
NB. testbasicsyr2k
NB.
NB. Description:
NB.   Test symmetric rank 2k operations:
NB.   - xSYR2K (BLAS)
NB.   by general and symmetric matrices
NB.
NB. Syntax:
NB.   log=. testbasicsyr2k A ; B ; CC
NB. where
NB.   A   - m×n-matrix
NB.   B   - m×n-matrix
NB.   CC  - k×k-matrix, C material
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)

testbasicsyr2k=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  'A B CC'=. y
  mn=. <./ 'm n'=. $ A
  CCmn=. (2 # mn) {. CC
  argsdl=. { (<"0 dcoeff) ; (< A) ; (< B) ; (<"0 dcoeff) ; < < CCmn [^:(m < n) CC
  argsdg=. { (<"0 dcoeff) ; (< A) ; (< B) ; (<"0 dcoeff) ; < < CCmn [^:(m > n) CC
  argszl=. { (<"0 zcoeff) ; (< A) ; (< B) ; (<"0 zcoeff) ; < < CCmn [^:(m < n) CC
  argszg=. { (<"0 zcoeff) ; (< A) ; (< B) ; (<"0 zcoeff) ; < < CCmn [^:(m > n) CC

  NB. for every i feed the tuple (alpha_i ; A ; B ; beta_i ; CC) to tmonad
  log=.          ('dsyr2kln_mtbla_' tmonad (]`]`nan`nan`(syr2kln chk5r2k trlpick)))@>"0 argsdl
  log=. log lcat ('dsyr2klt_mtbla_' tmonad (]`]`nan`nan`(syr2klt chk5r2k trlpick)))@>"0 argsdg
  log=. log lcat ('dsyr2kun_mtbla_' tmonad (]`]`nan`nan`(syr2kun chk5r2k trupick)))@>"0 argsdl
  log=. log lcat ('dsyr2kut_mtbla_' tmonad (]`]`nan`nan`(syr2kut chk5r2k trupick)))@>"0 argsdg
  log=. log lcat ('zsyr2kln_mtbla_' tmonad (]`]`nan`nan`(syr2kln chk5r2k trlpick)))@>"0 argszl
  log=. log lcat ('zsyr2klt_mtbla_' tmonad (]`]`nan`nan`(syr2klt chk5r2k trlpick)))@>"0 argszg
  log=. log lcat ('zsyr2kun_mtbla_' tmonad (]`]`nan`nan`(syr2kun chk5r2k trupick)))@>"0 argszl
  log=. log lcat ('zsyr2kut_mtbla_' tmonad (]`]`nan`nan`(syr2kut chk5r2k trupick)))@>"0 argszg
)

NB. ---------------------------------------------------------
NB. testbasicher2k
NB.
NB. Description:
NB.   Test the hermitian rank 2k operation:
NB.   - ZHER2K (BLAS)
NB.   by general and Hermitian matrices
NB.
NB. Syntax:
NB.   log=. testbasicher2k A ; B ; CC
NB. where
NB.   A   - m×n-matrix
NB.   B   - m×n-matrix
NB.   CC  - k×k-matrix, C material with real diagonal
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)

testbasicher2k=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  'A B CC'=. y
  mn=. <./ 'm n'=. $ A
  CCmn=. (2 # mn) {. CC
  argsl=. { (<"0 zcoeff) ; (< A) ; (< B) ; (<"0 dcoeff) ; < < CCmn [^:(m < n) CC
  argsg=. { (<"0 zcoeff) ; (< A) ; (< B) ; (<"0 dcoeff) ; < < CCmn [^:(m > n) CC

  NB. for every i feed the tuple (alpha_i ; A ; B ; beta_i ; CC) to tmonad
  log=.          ('zher2kln_mtbla_' tmonad (]`]`nan`nan`(her2kln chk5r2k trlpick)))@>"0 argsl
  log=. log lcat ('zher2klc_mtbla_' tmonad (]`]`nan`nan`(her2klc chk5r2k trlpick)))@>"0 argsg
  log=. log lcat ('zher2kun_mtbla_' tmonad (]`]`nan`nan`(her2kun chk5r2k trupick)))@>"0 argsl
  log=. log lcat ('zher2kuc_mtbla_' tmonad (]`]`nan`nan`(her2kuc chk5r2k trupick)))@>"0 argsg
)

NB. ---------------------------------------------------------
NB. testbasicr2k
NB.
NB. Description:
NB.   Adv. to make verb to test basic rank 2k operations
NB.
NB. Syntax:
NB.   log=. (mkmat testbasicr2k) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testbasicrk_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testbasicrk_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testbasicrk_mt_ 150 200

testbasicr2k=: 1 : 'testbasicher2k_mt_@(u ; u ; (9&o. upddiag_mt_)@u@(2 # >./)) ,&.>~ testbasicsyr2k_mt_@(u ; u ; u@(2 # >./)) [ require@''math/mt/external/blas/r2k'''

NB. ---------------------------------------------------------
NB. testbasicgemv
NB.
NB. Description:
NB.   Test matrix-vector operations:
NB.   - (+/ .*) (built-in)
NB.   - mp (math/misc/mathutil addon)
NB.   - xGEMV (BLAS)
NB.   by general matrix and vectors
NB.
NB. Syntax:
NB.   log=. testbasicgemv A ; x ; y
NB. where
NB.   A   - m×n-matrix
NB.   x   - k-vector
NB.   y   - k-vector
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)

testbasicgemv=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  inc=. 1 2 _1 _2
  'A x y'=. y
  'm n'=. $ A
  xm=. m {. x
  xn=. n {. x
  ym=. m {. y
  yn=. n {. y

  NB. test for the case: ('alpha beta inc'=. 1 0 1) and (op(A) = A)
  log=.          ('(+/ .*)'        tdyad  ((0&{::)`(1&{::)`0:`nan`nan`0:            ))                                                   A  ;    xn
  log=. log lcat ('mp'             tdyad  ((0&{::)`(1&{::)`0:`nan`nan`0:            ))                                                   A  ;    xn

  NB. for every i feed the tuple (alpha_i ; A ; expanded_x_i ; incx_i ; beta_i ; expanded_y_i ; incy_i) to tmonad
  log=. log lcat ('dgemvn_mtbla_' tmonad (         ]      `] `nan`nan`(gemvn chk1mv)))@(5 expand 6)@(2 expand 3)@>"0 { (<"0 dcoeff) ; (< A) ; (< xn) ; (<"0 inc) ; (<"0 dcoeff) ; (< ym) ; < <"0 inc
  log=. log lcat ('dgemvt_mtbla_' tmonad (         ]      `] `nan`nan`(gemvt chk1mv)))@(5 expand 6)@(2 expand 3)@>"0 { (<"0 dcoeff) ; (< A) ; (< xm) ; (<"0 inc) ; (<"0 dcoeff) ; (< yn) ; < <"0 inc
  log=. log lcat ('zgemvn_mtbla_' tmonad (         ]      `] `nan`nan`(gemvn chk1mv)))@(5 expand 6)@(2 expand 3)@>"0 { (<"0 zcoeff) ; (< A) ; (< xn) ; (<"0 inc) ; (<"0 zcoeff) ; (< ym) ; < <"0 inc
  log=. log lcat ('zgemvt_mtbla_' tmonad (         ]      `] `nan`nan`(gemvt chk1mv)))@(5 expand 6)@(2 expand 3)@>"0 { (<"0 zcoeff) ; (< A) ; (< xm) ; (<"0 inc) ; (<"0 zcoeff) ; (< ym) ; < <"0 inc
  log=. log lcat ('zgemvc_mtbla_' tmonad (         ]      `] `nan`nan`(gemvc chk1mv)))@(5 expand 6)@(2 expand 3)@>"0 { (<"0 zcoeff) ; (< A) ; (< xm) ; (<"0 inc) ; (<"0 zcoeff) ; (< ym) ; < <"0 inc
)

NB. ---------------------------------------------------------
NB. testbasichemv
NB.
NB. Description:
NB.   Test hermitian (symmetric) matrix-vector operations:
NB.   - DSYMV ZHEMV (BLAS)
NB.   by Hermitian (symmetric) matrix and vectors
NB.
NB. Syntax:
NB.   log=. testbasichemv AA ; x ; y
NB. where
NB.   AA  - n×n-matrix, A material with real diagonal
NB.   x   - n-vector
NB.   y   - n-vector
NB.   log - 6-vector of boxes, test log, see test.ijs

testbasichemv=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  inc=. 1 2 _1 _2
  'AA x y'=. y
  argsd=. { (<"0 dcoeff) ; (< AA) ; (< x) ; (<"0 inc) ; (<"0 dcoeff) ; (< y) ; < <"0 inc
  argsz=. { (<"0 zcoeff) ; (< AA) ; (< x) ; (<"0 inc) ; (<"0 zcoeff) ; (< y) ; < <"0 inc

  NB. for every i feed the tuple (alpha_i ; AA ; expanded_x_i ; incx_i ; beta_i ; expanded_y_i ; incy_i) to tmonad
  log=.          ('dsymvl_mtbla_' tmonad (]`]`nan`nan`(symvl chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 argsd
  log=. log lcat ('dsymvu_mtbla_' tmonad (]`]`nan`nan`(symvu chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 argsd
  log=. log lcat ('zhemvl_mtbla_' tmonad (]`]`nan`nan`(hemvl chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 argsz
  log=. log lcat ('zhemvu_mtbla_' tmonad (]`]`nan`nan`(hemvu chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 argsz
)

NB. ---------------------------------------------------------
NB. testbasictrmv
NB.
NB. Description:
NB.   Test matrix-vector operations:
NB.   - xTRMV (BLAS)
NB.   by triangular matrix and vector
NB.
NB. Syntax:
NB.   log=. testbasictrmv AA ; x
NB. where
NB.   AA  - n×n-matrix, A material
NB.   x   - n-vector
NB.   log - 6-vector of boxes, test log, see test.ijs

testbasictrmv=: 3 : 0
  inc=. 1 2 _1 _2
  'AA y'=. y
  args=. { (< AA) ; (< y) ; < <"0 inc

  NB. for every i feed the tuple (AA ; expanded_x_i ; incx_i) to tmonad
  log=.          ('dtrmvlnn_mtbla_' tmonad (]`]`nan`nan`(trmvlnn chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('dtrmvlnu_mtbla_' tmonad (]`]`nan`nan`(trmvlnu chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('dtrmvltn_mtbla_' tmonad (]`]`nan`nan`(trmvltn chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('dtrmvltu_mtbla_' tmonad (]`]`nan`nan`(trmvltu chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('dtrmvunn_mtbla_' tmonad (]`]`nan`nan`(trmvunn chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('dtrmvunu_mtbla_' tmonad (]`]`nan`nan`(trmvunu chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('dtrmvutn_mtbla_' tmonad (]`]`nan`nan`(trmvutn chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('dtrmvutu_mtbla_' tmonad (]`]`nan`nan`(trmvutu chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvlnn_mtbla_' tmonad (]`]`nan`nan`(trmvlnn chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvlnu_mtbla_' tmonad (]`]`nan`nan`(trmvlnu chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvltn_mtbla_' tmonad (]`]`nan`nan`(trmvltn chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvltu_mtbla_' tmonad (]`]`nan`nan`(trmvltu chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvlcn_mtbla_' tmonad (]`]`nan`nan`(trmvlcn chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvlcu_mtbla_' tmonad (]`]`nan`nan`(trmvlcu chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvunn_mtbla_' tmonad (]`]`nan`nan`(trmvunn chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvunu_mtbla_' tmonad (]`]`nan`nan`(trmvunu chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvutn_mtbla_' tmonad (]`]`nan`nan`(trmvutn chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvutu_mtbla_' tmonad (]`]`nan`nan`(trmvutu chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvucn_mtbla_' tmonad (]`]`nan`nan`(trmvucn chk3mv)))@(1 expand 2)@>"0 args
  log=. log lcat ('ztrmvucu_mtbla_' tmonad (]`]`nan`nan`(trmvucu chk3mv)))@(1 expand 2)@>"0 args
)

NB. ---------------------------------------------------------
NB. testbasicmv
NB.
NB. Description:
NB.   Adv. to make verb to test basic matrix-vector operations
NB.
NB. Syntax:
NB.   log=. (mkmat testbasicmv) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testbasicmv_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testbasicmv_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testbasicmv_mt_ 150 200

testbasicmv=: 1 : 'nolog_mt_`(testbasictrmv_mt_@(u ; u@{.) ,&.>~ testbasichemv_mt_@((9&o. upddiag_mt_)@u ; (u ; u)@{.))@.(=/) ,&.>~ testbasicgemv_mt_@(u ; (u ; u)@(>./)) [ require@''math/mt/external/blas/mv'''

NB. ---------------------------------------------------------
NB. testbasicgemm
NB.
NB. Description:
NB.   Test matrix-matrix operations:
NB.   - (+/ .*) (built-in)
NB.   - mp (math/misc/mathutil addon)
NB.   - xGEMM (BLAS)
NB.   - bli_xgemm (BLIS)
NB.   by general matrices
NB.
NB. Syntax:
NB.   log=. testbasicgemm As ; Bs ; C
NB. where
NB.   As  - m×(m+n)-matrix, A material
NB.   Bs  - (m+n)×n-matrix, B material
NB.   C   - m×n-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.
NB. Notes:
NB. - For real matrices and complex coefficients bli_xgemm
NB.   saves the real part of result in C and imagine part
NB.   will be lost. That is why it's sufficient to test
NB.   bli_xgemm by arguments of the same datatype only.

testbasicgemm=: 3 : 0
  'As Bs C'=. y
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  acoeff=. /:~ ~. zcoeff ,^:(JCMPX = 3!:0 C) dcoeff
  'm n'=. $ C
  ks=. /:~ ~. m (0 1 , (, >.@-:)@(, , +)) n  NB. 0,1,⌈m/2⌉,⌈n/2⌉,⌈(m+n)/2⌉,m,n,m+n
  As=. ks <@:({."0 1)"0 _ As                 NB. As[i] is m×k[i]-matrix
  argsdnn=. { (<"0 dcoeff) ;         As  ; (<    Bs) ; (<"0 dcoeff) ; < < C
  argsdnt=. { (<"0 dcoeff) ;         As  ; (< |: Bs) ; (<"0 dcoeff) ; < < C
  argsdtn=. { (<"0 dcoeff) ; (|: L:0 As) ; (<    Bs) ; (<"0 dcoeff) ; < < C
  argsdtt=. { (<"0 dcoeff) ; (|: L:0 As) ; (< |: Bs) ; (<"0 dcoeff) ; < < C
  NB. BLIS doesn't support mixed-datatype functionality with
  NB. non-zero imaginary part of alpha parameter yet, so
  NB. force alpha to be always real
  argsann=. { (<"0 dcoeff) ;         As  ; (<    Bs) ; (<"0 acoeff) ; < < C
  argsant=. { (<"0 dcoeff) ;         As  ; (< |: Bs) ; (<"0 acoeff) ; < < C
  argsanj=. { (<"0 dcoeff) ;         As  ; (< +  Bs) ; (<"0 acoeff) ; < < C
  argsanc=. { (<"0 dcoeff) ;         As  ; (< ct Bs) ; (<"0 acoeff) ; < < C
  argsatn=. { (<"0 dcoeff) ; (|: L:0 As) ; (<    Bs) ; (<"0 acoeff) ; < < C
  argsatt=. { (<"0 dcoeff) ; (|: L:0 As) ; (< |: Bs) ; (<"0 acoeff) ; < < C
  argsatj=. { (<"0 dcoeff) ; (|: L:0 As) ; (< +  Bs) ; (<"0 acoeff) ; < < C
  argsatc=. { (<"0 dcoeff) ; (|: L:0 As) ; (< ct Bs) ; (<"0 acoeff) ; < < C
  argsajn=. { (<"0 dcoeff) ; (+  L:0 As) ; (<    Bs) ; (<"0 acoeff) ; < < C
  argsajt=. { (<"0 dcoeff) ; (+  L:0 As) ; (< |: Bs) ; (<"0 acoeff) ; < < C
  argsajj=. { (<"0 dcoeff) ; (+  L:0 As) ; (< +  Bs) ; (<"0 acoeff) ; < < C
  argsajc=. { (<"0 dcoeff) ; (+  L:0 As) ; (< ct Bs) ; (<"0 acoeff) ; < < C
  argsacn=. { (<"0 dcoeff) ; (ct L:0 As) ; (<    Bs) ; (<"0 acoeff) ; < < C
  argsact=. { (<"0 dcoeff) ; (ct L:0 As) ; (< |: Bs) ; (<"0 acoeff) ; < < C
  argsacj=. { (<"0 dcoeff) ; (ct L:0 As) ; (< +  Bs) ; (<"0 acoeff) ; < < C
  argsacc=. { (<"0 dcoeff) ; (ct L:0 As) ; (< ct Bs) ; (<"0 acoeff) ; < < C

  NB. note: A_i and B_i shapes are related; to get them to
  NB.       match, a full fixed Bs was feeded to Catalogue
  NB.       ({) and now will be shrinked to the shape
  NB.       suitable before call to tmonad

  NB. test for the case: ('alpha beta'=. 1.0 0.0) and (op(A) = A)
  log=.          ('(+/ .*)'         tdyad  ((0&{::)`(1&{::)`0:`nan`nan`0:                           ))@(c (0 shrink 1)  {.   )@>"0 {                        As  ;  < <  Bs
  log=. log lcat ('mp'              tdyad  ((0&{::)`(1&{::)`0:`nan`nan`0:                           ))@(c (0 shrink 1)  {.   )@>"0 {                        As  ;  < <  Bs

  NB. for every i feed the tuple (alpha_i ; A_i ; B_i ; beta_i ; C) to tmonad

  log=. log lcat ('dgemmnn_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmnn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsdnn
  log=. log lcat ('dgemmnt_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmnt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsdnt
  log=. log lcat ('dgemmtn_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmtn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsdtn
  log=. log lcat ('dgemmtt_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmtt  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsdtt

  log=. log lcat ('zgemmnn_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmnn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 { (<"0 zcoeff) ;         As  ; (<    Bs) ; (<"0 zcoeff) ; < < C
  log=. log lcat ('zgemmnt_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmnt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ;         As  ; (< |: Bs) ; (<"0 zcoeff) ; < < C
  log=. log lcat ('zgemmnc_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmnc  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ;         As  ; (< ct Bs) ; (<"0 zcoeff) ; < < C
  log=. log lcat ('zgemmtn_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmtn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 { (<"0 zcoeff) ; (|: L:0 As) ; (<    Bs) ; (<"0 zcoeff) ; < < C
  log=. log lcat ('zgemmtt_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmtt  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (|: L:0 As) ; (< |: Bs) ; (<"0 zcoeff) ; < < C
  log=. log lcat ('zgemmtc_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmtc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (|: L:0 As) ; (< ct Bs) ; (<"0 zcoeff) ; < < C
  log=. log lcat ('zgemmcn_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmcn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 { (<"0 zcoeff) ; (ct L:0 As) ; (<    Bs) ; (<"0 zcoeff) ; < < C
  log=. log lcat ('zgemmct_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmct  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (ct L:0 As) ; (< |: Bs) ; (<"0 zcoeff) ; < < C
  log=. log lcat ('zgemmcc_mtbla_'  tmonad (        ]      `] `nan`nan`(             gemmcc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (ct L:0 As) ; (< ct Bs) ; (<"0 zcoeff) ; < < C

  log=. log lcat ('gemmnn_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmnn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  log=. log lcat ('gemmnt_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmnt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  log=. log lcat ('gemmnj_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmnj  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  log=. log lcat ('gemmnc_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmnc  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  log=. log lcat ('gemmtn_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmtn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  log=. log lcat ('gemmtt_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmtt  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  log=. log lcat ('gemmtj_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmtj  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  log=. log lcat ('gemmtc_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmtc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  log=. log lcat ('gemmjn_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmjn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  log=. log lcat ('gemmjt_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmjt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  log=. log lcat ('gemmjj_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmjj  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  log=. log lcat ('gemmjc_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmjc  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  log=. log lcat ('gemmcn_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmcn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  log=. log lcat ('gemmct_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmct  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  log=. log lcat ('gemmcj_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmcj  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  log=. log lcat ('gemmcc_mtbli_'   tmonad (        ]      `] `nan`nan`(             gemmcc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc

  log=. log lcat ('dgemmnn_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmnn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsdnn
  log=. log lcat ('dgemmnt_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmnt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsdnt
  log=. log lcat ('dgemmtn_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmtn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsdtn
  log=. log lcat ('dgemmtt_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmtt  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsdtt

  log=. log lcat ('zgemmnn_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmnn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  log=. log lcat ('zgemmnt_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmnt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  log=. log lcat ('zgemmnj_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmnj  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  log=. log lcat ('zgemmnc_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmnc  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  log=. log lcat ('zgemmtn_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmtn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  log=. log lcat ('zgemmtt_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmtt  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  log=. log lcat ('zgemmtj_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmtj  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  log=. log lcat ('zgemmtc_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmtc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  log=. log lcat ('zgemmjn_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmjn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  log=. log lcat ('zgemmjt_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmjt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  log=. log lcat ('zgemmjj_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmjj  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  log=. log lcat ('zgemmjc_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmjc  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  log=. log lcat ('zgemmcn_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmcn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  log=. log lcat ('zgemmct_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmct  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  log=. log lcat ('zgemmcj_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmcj  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  log=. log lcat ('zgemmcc_mtbli_'  tmonad (        ]      `] `nan`nan`(             gemmcc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc
)

NB. ---------------------------------------------------------
NB. testbasicgemmt
NB.
NB. Description:
NB.   Test matrix-matrix operations:
NB.   - bli_xgemmt (BLIS)
NB.   by general matrices
NB.
NB. Syntax:
NB.   testbasicgemmt As ; Bs ; C
NB. where
NB.   As  - m×(m+n)-matrix, A material
NB.   Bs  - (m+n)×m-matrix, B material
NB.   C   - m×m-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.
NB. Notes:
NB. - bli_xgemmt requres A,B,C to be of the same datatype.
NB. - For real matrices and complex coefficients bi_xgemmt
NB.   saves the real part of result in C and imagine part
NB.   will be lost. That is why it's sufficient to test
NB.   bli_xgemmt by arguments of the same datatype only.

testbasicgemmt=: 3 : 0
  NB. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NB. bli_gemmt fails with message:
  NB.   libblis: frame/base/bli_prune.c (line 130):
  NB.   libblis: Requested functionality not yet implemented.
  NB.   libblis: Aborting.
  NB. when is executed under JE
  nolog '' return.
  NB. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  'As Bs C'=. y
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  acoeff=. /:~ ~. zcoeff ,^:(JCMPX = 3!:0 C) dcoeff
  'n m'=. (-/ , ]) $ Bs
  ks=. /:~ ~. m (0 1 , (, >.@-:)@(, , +)) n  NB. 0,1,⌈m/2⌉,⌈n/2⌉,⌈(m+n)/2⌉,m,n,m+n
  As=. ks <@:({."0 1)"0 _ As                 NB. As[i] is m×k[i]-matrix
  argsdnn=. { (<"0 dcoeff) ;         As  ; (<    Bs) ; (<"0 dcoeff) ; < < C
  argsdnt=. { (<"0 dcoeff) ;         As  ; (< |: Bs) ; (<"0 dcoeff) ; < < C
  argsdtn=. { (<"0 dcoeff) ; (|: L:0 As) ; (<    Bs) ; (<"0 dcoeff) ; < < C
  argsdtt=. { (<"0 dcoeff) ; (|: L:0 As) ; (< |: Bs) ; (<"0 dcoeff) ; < < C
  argsann=. { (<"0 acoeff) ;         As  ; (<    Bs) ; (<"0 acoeff) ; < < C
  argsant=. { (<"0 acoeff) ;         As  ; (< |: Bs) ; (<"0 acoeff) ; < < C
  argsanj=. { (<"0 acoeff) ;         As  ; (< +  Bs) ; (<"0 acoeff) ; < < C
  argsanc=. { (<"0 acoeff) ;         As  ; (< ct Bs) ; (<"0 acoeff) ; < < C
  argsatn=. { (<"0 acoeff) ; (|: L:0 As) ; (<    Bs) ; (<"0 acoeff) ; < < C
  argsatt=. { (<"0 acoeff) ; (|: L:0 As) ; (< |: Bs) ; (<"0 acoeff) ; < < C
  argsatj=. { (<"0 acoeff) ; (|: L:0 As) ; (< +  Bs) ; (<"0 acoeff) ; < < C
  argsatc=. { (<"0 acoeff) ; (|: L:0 As) ; (< ct Bs) ; (<"0 acoeff) ; < < C
  argsajn=. { (<"0 acoeff) ; (+  L:0 As) ; (<    Bs) ; (<"0 acoeff) ; < < C
  argsajt=. { (<"0 acoeff) ; (+  L:0 As) ; (< |: Bs) ; (<"0 acoeff) ; < < C
  argsajj=. { (<"0 acoeff) ; (+  L:0 As) ; (< +  Bs) ; (<"0 acoeff) ; < < C
  argsajc=. { (<"0 acoeff) ; (+  L:0 As) ; (< ct Bs) ; (<"0 acoeff) ; < < C
  argsacn=. { (<"0 acoeff) ; (ct L:0 As) ; (<    Bs) ; (<"0 acoeff) ; < < C
  argsact=. { (<"0 acoeff) ; (ct L:0 As) ; (< |: Bs) ; (<"0 acoeff) ; < < C
  argsacj=. { (<"0 acoeff) ; (ct L:0 As) ; (< +  Bs) ; (<"0 acoeff) ; < < C
  argsacc=. { (<"0 acoeff) ; (ct L:0 As) ; (< ct Bs) ; (<"0 acoeff) ; < < C

  NB. note: A_i and B_i shapes are related; to get them to
  NB.       match, a full fixed Bs was feeded to Catalogue
  NB.       ({) and now will be shrinked to the shape
  NB.       suitable before call to tmonad

  NB. for every i feed the tuple (alpha_i ; A_i ; B_i ; beta_i ; C) to tmonad

  log=.          ('gemmlnn_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  log=. log lcat ('gemmlnt_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  log=. log lcat ('gemmlnj_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmnj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  log=. log lcat ('gemmlnc_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmnc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  log=. log lcat ('gemmltn_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  log=. log lcat ('gemmltt_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  log=. log lcat ('gemmltj_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmtj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  log=. log lcat ('gemmltc_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmtc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  log=. log lcat ('gemmljn_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmjn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  log=. log lcat ('gemmljt_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmjt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  log=. log lcat ('gemmljj_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmjj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  log=. log lcat ('gemmljc_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmjc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  log=. log lcat ('gemmlcn_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmcn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  log=. log lcat ('gemmlct_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmct) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  log=. log lcat ('gemmlcj_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmcj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  log=. log lcat ('gemmlcc_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmcc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc
  log=. log lcat ('gemmunn_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  log=. log lcat ('gemmunt_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  log=. log lcat ('gemmunj_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmnj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  log=. log lcat ('gemmunc_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmnc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  log=. log lcat ('gemmutn_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  log=. log lcat ('gemmutt_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  log=. log lcat ('gemmutj_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmtj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  log=. log lcat ('gemmutc_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmtc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  log=. log lcat ('gemmujn_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmjn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  log=. log lcat ('gemmujt_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmjt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  log=. log lcat ('gemmujj_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmjj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  log=. log lcat ('gemmujc_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmjc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  log=. log lcat ('gemmucn_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmcn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  log=. log lcat ('gemmuct_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmct) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  log=. log lcat ('gemmucj_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmcj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  log=. log lcat ('gemmucc_mtbli_'  tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmcc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc

  log=. log lcat ('dgemmlnn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsdnn
  log=. log lcat ('dgemmlnt_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsdnt
  log=. log lcat ('dgemmltn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsdtn
  log=. log lcat ('dgemmltt_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsdtt
  log=. log lcat ('dgemmunn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsdnn
  log=. log lcat ('dgemmunt_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsdnt
  log=. log lcat ('dgemmutn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsdtn
  log=. log lcat ('dgemmutt_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsdtt

  log=. log lcat ('zgemmlnn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  log=. log lcat ('zgemmlnt_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  log=. log lcat ('zgemmlnj_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmnj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  log=. log lcat ('zgemmlnc_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmnc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  log=. log lcat ('zgemmltn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  log=. log lcat ('zgemmltt_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  log=. log lcat ('zgemmltj_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmtj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  log=. log lcat ('zgemmltc_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmtc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  log=. log lcat ('zgemmljn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmjn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  log=. log lcat ('zgemmljt_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmjt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  log=. log lcat ('zgemmljj_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmjj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  log=. log lcat ('zgemmljc_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmjc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  log=. log lcat ('zgemmlcn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmcn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  log=. log lcat ('zgemmlct_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmct) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  log=. log lcat ('zgemmlcj_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmcj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  log=. log lcat ('zgemmlcc_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: suxly gemmcc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc
  log=. log lcat ('zgemmunn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  log=. log lcat ('zgemmunt_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  log=. log lcat ('zgemmunj_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmnj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  log=. log lcat ('zgemmunc_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmnc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  log=. log lcat ('zgemmutn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  log=. log lcat ('zgemmutt_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  log=. log lcat ('zgemmutj_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmtj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  log=. log lcat ('zgemmutc_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmtc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  log=. log lcat ('zgemmujn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmjn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  log=. log lcat ('zgemmujt_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmjt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  log=. log lcat ('zgemmujj_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmjj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  log=. log lcat ('zgemmujc_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmjc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  log=. log lcat ('zgemmucn_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmcn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  log=. log lcat ('zgemmuct_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmct) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  log=. log lcat ('zgemmucj_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmcj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  log=. log lcat ('zgemmucc_mtbli_' tmonad (        ]      `] `nan`nan`((4&{:: slxuy gemmcc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc
)

NB. ---------------------------------------------------------
NB. testbasicsymm
NB.
NB. Description:
NB.   Test symmetric matrix-vector operations:
NB.   - xSYMM (BLAS)
NB.   - bli_xsymm (BLIS)
NB.   by symmetric matrix
NB.
NB. Syntax:
NB.   log=. testbasicsymm AA ; B ; C
NB. where
NB.   AA  - k×k-matrix, A material
NB.   B   - m×n-matrix
NB.   C   - m×n-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - For real matrices and complex coefficients bli_xsymm
NB.   saves the real part of result in C and imagine part
NB.   will be lost. That is why it's sufficient to test
NB.   bli_xsymm by arguments of the same datatype only.

testbasicsymm=: 3 : 0
  'AA B C'=. y
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  acoeff=. /:~ ~. zcoeff ,^:(JCMPX = 3!:0 C) dcoeff
  'm n'=. $ C
  Am=. (2 # m) {. AA
  An=. (2 # n) {. AA
  argsdm=. { (<"0 dcoeff) ; (< Am) ; (< B) ; (<"0 dcoeff) ; < < C
  argsdn=. { (<"0 dcoeff) ; (< An) ; (< B) ; (<"0 dcoeff) ; < < C
  argszm=. { (<"0 zcoeff) ; (< Am) ; (< B) ; (<"0 zcoeff) ; < < C
  argszn=. { (<"0 zcoeff) ; (< An) ; (< B) ; (<"0 zcoeff) ; < < C
  argsam=. { (<"0 acoeff) ; (< Am) ; (< B) ; (<"0 acoeff) ; < < C
  argsan=. { (<"0 acoeff) ; (< An) ; (< B) ; (<"0 acoeff) ; < < C

  NB. for every i feed the tuple (alpha_i ; AA ; B ; beta_i ; C) to tmonad

  log=.          ('dsymmll_mtbla_'   tmonad (]`]`nan`nan`(symmllnn chk2mm)))@>"0 argsdm
  log=. log lcat ('dsymmlu_mtbla_'   tmonad (]`]`nan`nan`(symmlunn chk2mm)))@>"0 argsdm
  log=. log lcat ('dsymmrl_mtbla_'   tmonad (]`]`nan`nan`(symmrlnn chk2mm)))@>"0 argsdn
  log=. log lcat ('dsymmru_mtbla_'   tmonad (]`]`nan`nan`(symmrunn chk2mm)))@>"0 argsdn

  log=. log lcat ('zsymmll_mtbla_'   tmonad (]`]`nan`nan`(symmllnn chk2mm)))@>"0 argszm
  log=. log lcat ('zsymmlu_mtbla_'   tmonad (]`]`nan`nan`(symmlunn chk2mm)))@>"0 argszm
  log=. log lcat ('zsymmrl_mtbla_'   tmonad (]`]`nan`nan`(symmrlnn chk2mm)))@>"0 argszn
  log=. log lcat ('zsymmru_mtbla_'   tmonad (]`]`nan`nan`(symmrunn chk2mm)))@>"0 argszn

  log=. log lcat ('symmllnn_mtbli_'  tmonad (]`]`nan`nan`(symmllnn chk2mm)))@>"0 argsam
  log=. log lcat ('symmllnt_mtbli_'  tmonad (]`]`nan`nan`(symmllnt chk2mm)))@>"0 argsam
  log=. log lcat ('symmllnj_mtbli_'  tmonad (]`]`nan`nan`(symmllnj chk2mm)))@>"0 argsam
  log=. log lcat ('symmllnc_mtbli_'  tmonad (]`]`nan`nan`(symmllnc chk2mm)))@>"0 argsam
  log=. log lcat ('symmlljn_mtbli_'  tmonad (]`]`nan`nan`(symmlljn chk2mm)))@>"0 argsam
  log=. log lcat ('symmlljt_mtbli_'  tmonad (]`]`nan`nan`(symmlljt chk2mm)))@>"0 argsam
  log=. log lcat ('symmlljj_mtbli_'  tmonad (]`]`nan`nan`(symmlljj chk2mm)))@>"0 argsam
  log=. log lcat ('symmlljc_mtbli_'  tmonad (]`]`nan`nan`(symmlljc chk2mm)))@>"0 argsam
  log=. log lcat ('symmlunn_mtbli_'  tmonad (]`]`nan`nan`(symmlunn chk2mm)))@>"0 argsam
  log=. log lcat ('symmlunt_mtbli_'  tmonad (]`]`nan`nan`(symmlunt chk2mm)))@>"0 argsam
  log=. log lcat ('symmlunj_mtbli_'  tmonad (]`]`nan`nan`(symmlunj chk2mm)))@>"0 argsam
  log=. log lcat ('symmlunc_mtbli_'  tmonad (]`]`nan`nan`(symmlunc chk2mm)))@>"0 argsam
  log=. log lcat ('symmlujn_mtbli_'  tmonad (]`]`nan`nan`(symmlujn chk2mm)))@>"0 argsam
  log=. log lcat ('symmlujt_mtbli_'  tmonad (]`]`nan`nan`(symmlujt chk2mm)))@>"0 argsam
  log=. log lcat ('symmlujj_mtbli_'  tmonad (]`]`nan`nan`(symmlujj chk2mm)))@>"0 argsam
  log=. log lcat ('symmlujc_mtbli_'  tmonad (]`]`nan`nan`(symmlujc chk2mm)))@>"0 argsam
  log=. log lcat ('symmrlnn_mtbli_'  tmonad (]`]`nan`nan`(symmrlnn chk2mm)))@>"0 argsan
  log=. log lcat ('symmrlnt_mtbli_'  tmonad (]`]`nan`nan`(symmrlnt chk2mm)))@>"0 argsan
  log=. log lcat ('symmrlnj_mtbli_'  tmonad (]`]`nan`nan`(symmrlnj chk2mm)))@>"0 argsan
  log=. log lcat ('symmrlnc_mtbli_'  tmonad (]`]`nan`nan`(symmrlnc chk2mm)))@>"0 argsan
  log=. log lcat ('symmrljn_mtbli_'  tmonad (]`]`nan`nan`(symmrljn chk2mm)))@>"0 argsan
  log=. log lcat ('symmrljt_mtbli_'  tmonad (]`]`nan`nan`(symmrljt chk2mm)))@>"0 argsan
  log=. log lcat ('symmrljj_mtbli_'  tmonad (]`]`nan`nan`(symmrljj chk2mm)))@>"0 argsan
  log=. log lcat ('symmrljc_mtbli_'  tmonad (]`]`nan`nan`(symmrljc chk2mm)))@>"0 argsan
  log=. log lcat ('symmrunn_mtbli_'  tmonad (]`]`nan`nan`(symmrunn chk2mm)))@>"0 argsan
  log=. log lcat ('symmrunt_mtbli_'  tmonad (]`]`nan`nan`(symmrunt chk2mm)))@>"0 argsan
  log=. log lcat ('symmrunj_mtbli_'  tmonad (]`]`nan`nan`(symmrunj chk2mm)))@>"0 argsan
  log=. log lcat ('symmrunc_mtbli_'  tmonad (]`]`nan`nan`(symmrunc chk2mm)))@>"0 argsan
  log=. log lcat ('symmrujn_mtbli_'  tmonad (]`]`nan`nan`(symmrujn chk2mm)))@>"0 argsan
  log=. log lcat ('symmrujt_mtbli_'  tmonad (]`]`nan`nan`(symmrujt chk2mm)))@>"0 argsan
  log=. log lcat ('symmrujj_mtbli_'  tmonad (]`]`nan`nan`(symmrujj chk2mm)))@>"0 argsan
  log=. log lcat ('symmrujc_mtbli_'  tmonad (]`]`nan`nan`(symmrujc chk2mm)))@>"0 argsan

  log=. log lcat ('dsymmllnn_mtbli_' tmonad (]`]`nan`nan`(symmllnn chk2mm)))@>"0 argsdm
  log=. log lcat ('dsymmllnt_mtbli_' tmonad (]`]`nan`nan`(symmllnt chk2mm)))@>"0 argsdm
  log=. log lcat ('dsymmlunn_mtbli_' tmonad (]`]`nan`nan`(symmlunn chk2mm)))@>"0 argsdm
  log=. log lcat ('dsymmlunt_mtbli_' tmonad (]`]`nan`nan`(symmlunt chk2mm)))@>"0 argsdm
  log=. log lcat ('dsymmrlnn_mtbli_' tmonad (]`]`nan`nan`(symmrlnn chk2mm)))@>"0 argsdn
  log=. log lcat ('dsymmrlnt_mtbli_' tmonad (]`]`nan`nan`(symmrlnt chk2mm)))@>"0 argsdn
  log=. log lcat ('dsymmrunn_mtbli_' tmonad (]`]`nan`nan`(symmrunn chk2mm)))@>"0 argsdn
  log=. log lcat ('dsymmrunt_mtbli_' tmonad (]`]`nan`nan`(symmrunt chk2mm)))@>"0 argsdn

  log=. log lcat ('zsymmllnn_mtbli_' tmonad (]`]`nan`nan`(symmllnn chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmllnt_mtbli_' tmonad (]`]`nan`nan`(symmllnt chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmllnj_mtbli_' tmonad (]`]`nan`nan`(symmllnj chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmllnc_mtbli_' tmonad (]`]`nan`nan`(symmllnc chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlljn_mtbli_' tmonad (]`]`nan`nan`(symmlljn chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlljt_mtbli_' tmonad (]`]`nan`nan`(symmlljt chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlljj_mtbli_' tmonad (]`]`nan`nan`(symmlljj chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlljc_mtbli_' tmonad (]`]`nan`nan`(symmlljc chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlunn_mtbli_' tmonad (]`]`nan`nan`(symmlunn chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlunt_mtbli_' tmonad (]`]`nan`nan`(symmlunt chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlunj_mtbli_' tmonad (]`]`nan`nan`(symmlunj chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlunc_mtbli_' tmonad (]`]`nan`nan`(symmlunc chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlujn_mtbli_' tmonad (]`]`nan`nan`(symmlujn chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlujt_mtbli_' tmonad (]`]`nan`nan`(symmlujt chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlujj_mtbli_' tmonad (]`]`nan`nan`(symmlujj chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmlujc_mtbli_' tmonad (]`]`nan`nan`(symmlujc chk2mm)))@>"0 argsam
  log=. log lcat ('zsymmrlnn_mtbli_' tmonad (]`]`nan`nan`(symmrlnn chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrlnt_mtbli_' tmonad (]`]`nan`nan`(symmrlnt chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrlnj_mtbli_' tmonad (]`]`nan`nan`(symmrlnj chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrlnc_mtbli_' tmonad (]`]`nan`nan`(symmrlnc chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrljn_mtbli_' tmonad (]`]`nan`nan`(symmrljn chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrljt_mtbli_' tmonad (]`]`nan`nan`(symmrljt chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrljj_mtbli_' tmonad (]`]`nan`nan`(symmrljj chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrljc_mtbli_' tmonad (]`]`nan`nan`(symmrljc chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrunn_mtbli_' tmonad (]`]`nan`nan`(symmrunn chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrunt_mtbli_' tmonad (]`]`nan`nan`(symmrunt chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrunj_mtbli_' tmonad (]`]`nan`nan`(symmrunj chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrunc_mtbli_' tmonad (]`]`nan`nan`(symmrunc chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrujn_mtbli_' tmonad (]`]`nan`nan`(symmrujn chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrujt_mtbli_' tmonad (]`]`nan`nan`(symmrujt chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrujj_mtbli_' tmonad (]`]`nan`nan`(symmrujj chk2mm)))@>"0 argsan
  log=. log lcat ('zsymmrujc_mtbli_' tmonad (]`]`nan`nan`(symmrujc chk2mm)))@>"0 argsan
)

NB. ---------------------------------------------------------
NB. testbasichemm
NB.
NB. Description:
NB.   Test hermitian matrix-vector operation:
NB.   - ZHEMM (BLAS)
NB.   - bli_xhemm (BLIS)
NB.   by Hermitian matrix
NB.
NB. Syntax:
NB.   log=. testbasichemm AA ; B ; C
NB. where
NB.   AA  - k×k-matrix, A material with real diagonal
NB.   B   - m×n-matrix
NB.   C   - m×n-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - For real matrices and complex coefficients bli_xhemm
NB.   saves the real part of result in C and imagine part
NB.   will be lost. That is why it's sufficient to test
NB.   bli_xhemm by arguments of the same datatype only.

testbasichemm=: 3 : 0
  'AA B C'=. y
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  acoeff=. /:~ ~. zcoeff ,^:(JCMPX = 3!:0 C) dcoeff
  'm n'=. $ C
  Am=. (2 # m) {. AA
  An=. (2 # n) {. AA
  argszm=. { (<"0 zcoeff) ; (< Am) ; (< B) ; (<"0 zcoeff) ; < < C
  argszn=. { (<"0 zcoeff) ; (< An) ; (< B) ; (<"0 zcoeff) ; < < C
  argsam=. { (<"0 acoeff) ; (< Am) ; (< B) ; (<"0 acoeff) ; < < C
  argsan=. { (<"0 acoeff) ; (< An) ; (< B) ; (<"0 acoeff) ; < < C

  NB. for every i feed the tuple (alpha_i ; Ax ; B ; beta_i ; C) to tmonad

  log=.          ('zhemmll_mtbla_'   tmonad (]`]`nan`nan`(hemmllnn chk2mm)))@>"0 argszm
  log=. log lcat ('zhemmlu_mtbla_'   tmonad (]`]`nan`nan`(hemmlunn chk2mm)))@>"0 argszm
  log=. log lcat ('zhemmrl_mtbla_'   tmonad (]`]`nan`nan`(hemmrlnn chk2mm)))@>"0 argszn
  log=. log lcat ('zhemmru_mtbla_'   tmonad (]`]`nan`nan`(hemmrunn chk2mm)))@>"0 argszn

  log=. log lcat ('hemmllnn_mtbli_'  tmonad (]`]`nan`nan`(hemmllnn chk2mm)))@>"0 argsam
  log=. log lcat ('hemmllnt_mtbli_'  tmonad (]`]`nan`nan`(hemmllnt chk2mm)))@>"0 argsam
  log=. log lcat ('hemmllnj_mtbli_'  tmonad (]`]`nan`nan`(hemmllnj chk2mm)))@>"0 argsam
  log=. log lcat ('hemmllnc_mtbli_'  tmonad (]`]`nan`nan`(hemmllnc chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlljn_mtbli_'  tmonad (]`]`nan`nan`(hemmlljn chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlljt_mtbli_'  tmonad (]`]`nan`nan`(hemmlljt chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlljj_mtbli_'  tmonad (]`]`nan`nan`(hemmlljj chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlljc_mtbli_'  tmonad (]`]`nan`nan`(hemmlljc chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlunn_mtbli_'  tmonad (]`]`nan`nan`(hemmlunn chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlunt_mtbli_'  tmonad (]`]`nan`nan`(hemmlunt chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlunj_mtbli_'  tmonad (]`]`nan`nan`(hemmlunj chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlunc_mtbli_'  tmonad (]`]`nan`nan`(hemmlunc chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlujn_mtbli_'  tmonad (]`]`nan`nan`(hemmlujn chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlujt_mtbli_'  tmonad (]`]`nan`nan`(hemmlujt chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlujj_mtbli_'  tmonad (]`]`nan`nan`(hemmlujj chk2mm)))@>"0 argsam
  log=. log lcat ('hemmlujc_mtbli_'  tmonad (]`]`nan`nan`(hemmlujc chk2mm)))@>"0 argsam
  log=. log lcat ('hemmrlnn_mtbli_'  tmonad (]`]`nan`nan`(hemmrlnn chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrlnt_mtbli_'  tmonad (]`]`nan`nan`(hemmrlnt chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrlnj_mtbli_'  tmonad (]`]`nan`nan`(hemmrlnj chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrlnc_mtbli_'  tmonad (]`]`nan`nan`(hemmrlnc chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrljn_mtbli_'  tmonad (]`]`nan`nan`(hemmrljn chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrljt_mtbli_'  tmonad (]`]`nan`nan`(hemmrljt chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrljj_mtbli_'  tmonad (]`]`nan`nan`(hemmrljj chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrljc_mtbli_'  tmonad (]`]`nan`nan`(hemmrljc chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrunn_mtbli_'  tmonad (]`]`nan`nan`(hemmrunn chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrunt_mtbli_'  tmonad (]`]`nan`nan`(hemmrunt chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrunj_mtbli_'  tmonad (]`]`nan`nan`(hemmrunj chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrunc_mtbli_'  tmonad (]`]`nan`nan`(hemmrunc chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrujn_mtbli_'  tmonad (]`]`nan`nan`(hemmrujn chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrujt_mtbli_'  tmonad (]`]`nan`nan`(hemmrujt chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrujj_mtbli_'  tmonad (]`]`nan`nan`(hemmrujj chk2mm)))@>"0 argsan
  log=. log lcat ('hemmrujc_mtbli_'  tmonad (]`]`nan`nan`(hemmrujc chk2mm)))@>"0 argsan

  log=. log lcat ('zhemmllnn_mtbli_' tmonad (]`]`nan`nan`(hemmllnn chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmllnt_mtbli_' tmonad (]`]`nan`nan`(hemmllnt chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmllnj_mtbli_' tmonad (]`]`nan`nan`(hemmllnj chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmllnc_mtbli_' tmonad (]`]`nan`nan`(hemmllnc chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlljn_mtbli_' tmonad (]`]`nan`nan`(hemmlljn chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlljt_mtbli_' tmonad (]`]`nan`nan`(hemmlljt chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlljj_mtbli_' tmonad (]`]`nan`nan`(hemmlljj chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlljc_mtbli_' tmonad (]`]`nan`nan`(hemmlljc chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlunn_mtbli_' tmonad (]`]`nan`nan`(hemmlunn chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlunt_mtbli_' tmonad (]`]`nan`nan`(hemmlunt chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlunj_mtbli_' tmonad (]`]`nan`nan`(hemmlunj chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlunc_mtbli_' tmonad (]`]`nan`nan`(hemmlunc chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlujn_mtbli_' tmonad (]`]`nan`nan`(hemmlujn chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlujt_mtbli_' tmonad (]`]`nan`nan`(hemmlujt chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlujj_mtbli_' tmonad (]`]`nan`nan`(hemmlujj chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmlujc_mtbli_' tmonad (]`]`nan`nan`(hemmlujc chk2mm)))@>"0 argsam
  log=. log lcat ('zhemmrlnn_mtbli_' tmonad (]`]`nan`nan`(hemmrlnn chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrlnt_mtbli_' tmonad (]`]`nan`nan`(hemmrlnt chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrlnj_mtbli_' tmonad (]`]`nan`nan`(hemmrlnj chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrlnc_mtbli_' tmonad (]`]`nan`nan`(hemmrlnc chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrljn_mtbli_' tmonad (]`]`nan`nan`(hemmrljn chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrljt_mtbli_' tmonad (]`]`nan`nan`(hemmrljt chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrljj_mtbli_' tmonad (]`]`nan`nan`(hemmrljj chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrljc_mtbli_' tmonad (]`]`nan`nan`(hemmrljc chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrunn_mtbli_' tmonad (]`]`nan`nan`(hemmrunn chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrunt_mtbli_' tmonad (]`]`nan`nan`(hemmrunt chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrunj_mtbli_' tmonad (]`]`nan`nan`(hemmrunj chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrunc_mtbli_' tmonad (]`]`nan`nan`(hemmrunc chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrujn_mtbli_' tmonad (]`]`nan`nan`(hemmrujn chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrujt_mtbli_' tmonad (]`]`nan`nan`(hemmrujt chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrujj_mtbli_' tmonad (]`]`nan`nan`(hemmrujj chk2mm)))@>"0 argsan
  log=. log lcat ('zhemmrujc_mtbli_' tmonad (]`]`nan`nan`(hemmrujc chk2mm)))@>"0 argsan
)

NB. ---------------------------------------------------------
NB. testbasictrmm
NB.
NB. Description:
NB.   Test matrix-matrix operations:
NB.   - xTRMM (BLAS)
NB.   - bli_xtrmm bli_xtrmm3 (BLIS)
NB.   by triangular matrix
NB.
NB. Syntax:
NB.   log=. testbasictrmm AA ; B ; С
NB. where
NB.   AA  - k×k-matrix, A material
NB.   B   - m×n-matrix
NB.   C   - m×n-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - C is used to test bli_xtrmm3 only
NB. - For real matrices and complex coefficients bli_xtrmm
NB.   saves the real part of result in C and imagine part
NB.   will be lost. That is why it's sufficient to test
NB.   bli_xtrmm by arguments of the same datatype only.

testbasictrmm=: 3 : 0
  'AA B C'=. y
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  acoeff=. /:~ ~. zcoeff ,^:(JCMPX = 3!:0 C) dcoeff
  'm n'=. $ C
  Am=. (2 # m) {. AA
  An=. (2 # n) {. AA
  Bt=. |: B
  argsdm=.  { (<"0 dcoeff) ; (< Am) ; < < B
  argsdn=.  { (<"0 dcoeff) ; (< An) ; < < B
  argszm=.  { (<"0 zcoeff) ; (< Am) ; < < B
  argszn=.  { (<"0 zcoeff) ; (< An) ; < < B
  argsam=.  { (<"0 acoeff) ; (< Am) ; < < B
  argsan=.  { (<"0 acoeff) ; (< An) ; < < B
  argsdmm=. { (<"0 dcoeff) ; (< Am) ; (<  B ) ; (<"0 dcoeff) ; < < C
  argsdmn=. { (<"0 dcoeff) ; (< Am) ; (<  Bt) ; (<"0 dcoeff) ; < < C
  argsdnm=. { (<"0 dcoeff) ; (< An) ; (<  B ) ; (<"0 dcoeff) ; < < C
  argsdnn=. { (<"0 dcoeff) ; (< An) ; (<  Bt) ; (<"0 dcoeff) ; < < C
  argszmm=. { (<"0 zcoeff) ; (< Am) ; (<  B ) ; (<"0 zcoeff) ; < < C
  argszmn=. { (<"0 zcoeff) ; (< Am) ; (<  Bt) ; (<"0 zcoeff) ; < < C
  argsznm=. { (<"0 zcoeff) ; (< An) ; (<  B ) ; (<"0 zcoeff) ; < < C
  argsznn=. { (<"0 zcoeff) ; (< An) ; (<  Bt) ; (<"0 zcoeff) ; < < C
  argsamm=. { (<"0 acoeff) ; (< Am) ; (<  B ) ; (<"0 acoeff) ; < < C
  argsamn=. { (<"0 acoeff) ; (< Am) ; (<  Bt) ; (<"0 acoeff) ; < < C
  argsanm=. { (<"0 acoeff) ; (< An) ; (<  B ) ; (<"0 acoeff) ; < < C
  argsann=. { (<"0 acoeff) ; (< An) ; (<  Bt) ; (<"0 acoeff) ; < < C

  NB. for every i feed the tuple (alpha_i ; Ax ; B) to tmonad

  log=.          ('dtrmmllnn_mtbla_'  tmonad (]`]`nan`nan`(trmmllnn  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmllnu_mtbla_'  tmonad (]`]`nan`nan`(trmmllnu  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlltn_mtbla_'  tmonad (]`]`nan`nan`(trmmlltn  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlltu_mtbla_'  tmonad (]`]`nan`nan`(trmmlltu  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlunn_mtbla_'  tmonad (]`]`nan`nan`(trmmlunn  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlunu_mtbla_'  tmonad (]`]`nan`nan`(trmmlunu  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlutn_mtbla_'  tmonad (]`]`nan`nan`(trmmlutn  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlutu_mtbla_'  tmonad (]`]`nan`nan`(trmmlutu  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmrlnn_mtbla_'  tmonad (]`]`nan`nan`(trmmrlnn  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrlnu_mtbla_'  tmonad (]`]`nan`nan`(trmmrlnu  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrltn_mtbla_'  tmonad (]`]`nan`nan`(trmmrltn  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrltu_mtbla_'  tmonad (]`]`nan`nan`(trmmrltu  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrunn_mtbla_'  tmonad (]`]`nan`nan`(trmmrunn  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrunu_mtbla_'  tmonad (]`]`nan`nan`(trmmrunu  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrutn_mtbla_'  tmonad (]`]`nan`nan`(trmmrutn  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrutu_mtbla_'  tmonad (]`]`nan`nan`(trmmrutu  chk3mm)))@>"0 argsdn

  log=. log lcat ('ztrmmllnn_mtbla_'  tmonad (]`]`nan`nan`(trmmllnn  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmllnu_mtbla_'  tmonad (]`]`nan`nan`(trmmllnu  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmlltn_mtbla_'  tmonad (]`]`nan`nan`(trmmlltn  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmlltu_mtbla_'  tmonad (]`]`nan`nan`(trmmlltu  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmllcn_mtbla_'  tmonad (]`]`nan`nan`(trmmllcn  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmllcu_mtbla_'  tmonad (]`]`nan`nan`(trmmllcu  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmlunn_mtbla_'  tmonad (]`]`nan`nan`(trmmlunn  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmlunu_mtbla_'  tmonad (]`]`nan`nan`(trmmlunu  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmlutn_mtbla_'  tmonad (]`]`nan`nan`(trmmlutn  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmlutu_mtbla_'  tmonad (]`]`nan`nan`(trmmlutu  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmlucn_mtbla_'  tmonad (]`]`nan`nan`(trmmlucn  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmlucu_mtbla_'  tmonad (]`]`nan`nan`(trmmlucu  chk3mm)))@>"0 argszm
  log=. log lcat ('ztrmmrlnn_mtbla_'  tmonad (]`]`nan`nan`(trmmrlnn  chk3mm)))@>"0 argszn
  log=. log lcat ('ztrmmrlnu_mtbla_'  tmonad (]`]`nan`nan`(trmmrlnu  chk3mm)))@>"0 argszn
  log=. log lcat ('ztrmmrltn_mtbla_'  tmonad (]`]`nan`nan`(trmmrltn  chk3mm)))@>"0 argszn
  log=. log lcat ('ztrmmrltu_mtbla_'  tmonad (]`]`nan`nan`(trmmrltu  chk3mm)))@>"0 argszn
  log=. log lcat ('ztrmmrlcn_mtbla_'  tmonad (]`]`nan`nan`(trmmrlcn  chk3mm)))@>"0 argszn
  log=. log lcat ('ztrmmrlcu_mtbla_'  tmonad (]`]`nan`nan`(trmmrlcu  chk3mm)))@>"0 argszn
  log=. log lcat ('ztrmmrunn_mtbla_'  tmonad (]`]`nan`nan`(trmmrunn  chk3mm)))@>"0 argszn
  log=. log lcat ('ztrmmrunu_mtbla_'  tmonad (]`]`nan`nan`(trmmrunu  chk3mm)))@>"0 argszn
  log=. log lcat ('ztrmmrutn_mtbla_'  tmonad (]`]`nan`nan`(trmmrutn  chk3mm)))@>"0 argszn
  log=. log lcat ('ztrmmrutu_mtbla_'  tmonad (]`]`nan`nan`(trmmrutu  chk3mm)))@>"0 argszn
  log=. log lcat ('ztrmmrucn_mtbla_'  tmonad (]`]`nan`nan`(trmmrucn  chk3mm)))@>"0 argszn
  log=. log lcat ('ztrmmrucu_mtbla_'  tmonad (]`]`nan`nan`(trmmrucu  chk3mm)))@>"0 argszn

  log=. log lcat ('trmmllnn_mtbli_'   tmonad (]`]`nan`nan`(trmmllnn  chk3mm)))@>"0 argsam
  log=. log lcat ('trmmllnu_mtbli_'   tmonad (]`]`nan`nan`(trmmllnu  chk3mm)))@>"0 argsam
  log=. log lcat ('trmmlltn_mtbli_'   tmonad (]`]`nan`nan`(trmmlltn  chk3mm)))@>"0 argsam
  log=. log lcat ('trmmlltu_mtbli_'   tmonad (]`]`nan`nan`(trmmlltu  chk3mm)))@>"0 argsam
  log=. log lcat ('trmmlunn_mtbli_'   tmonad (]`]`nan`nan`(trmmlunn  chk3mm)))@>"0 argsam
  log=. log lcat ('trmmlunu_mtbli_'   tmonad (]`]`nan`nan`(trmmlunu  chk3mm)))@>"0 argsam
  log=. log lcat ('trmmlutn_mtbli_'   tmonad (]`]`nan`nan`(trmmlutn  chk3mm)))@>"0 argsam
  log=. log lcat ('trmmlutu_mtbli_'   tmonad (]`]`nan`nan`(trmmlutu  chk3mm)))@>"0 argsam
  log=. log lcat ('trmmrlnn_mtbli_'   tmonad (]`]`nan`nan`(trmmrlnn  chk3mm)))@>"0 argsan
  log=. log lcat ('trmmrlnu_mtbli_'   tmonad (]`]`nan`nan`(trmmrlnu  chk3mm)))@>"0 argsan
  log=. log lcat ('trmmrltn_mtbli_'   tmonad (]`]`nan`nan`(trmmrltn  chk3mm)))@>"0 argsan
  log=. log lcat ('trmmrltu_mtbli_'   tmonad (]`]`nan`nan`(trmmrltu  chk3mm)))@>"0 argsan
  log=. log lcat ('trmmrunn_mtbli_'   tmonad (]`]`nan`nan`(trmmrunn  chk3mm)))@>"0 argsan
  log=. log lcat ('trmmrunu_mtbli_'   tmonad (]`]`nan`nan`(trmmrunu  chk3mm)))@>"0 argsan
  log=. log lcat ('trmmrutn_mtbli_'   tmonad (]`]`nan`nan`(trmmrutn  chk3mm)))@>"0 argsan
  log=. log lcat ('trmmrutu_mtbli_'   tmonad (]`]`nan`nan`(trmmrutu  chk3mm)))@>"0 argsan

  log=. log lcat ('dtrmmllnn_mtbli_'  tmonad (]`]`nan`nan`(trmmllnn  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmllnu_mtbli_'  tmonad (]`]`nan`nan`(trmmllnu  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlltn_mtbli_'  tmonad (]`]`nan`nan`(trmmlltn  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlltu_mtbli_'  tmonad (]`]`nan`nan`(trmmlltu  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlunn_mtbli_'  tmonad (]`]`nan`nan`(trmmlunn  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlunu_mtbli_'  tmonad (]`]`nan`nan`(trmmlunu  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlutn_mtbli_'  tmonad (]`]`nan`nan`(trmmlutn  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmlutu_mtbli_'  tmonad (]`]`nan`nan`(trmmlutu  chk3mm)))@>"0 argsdm
  log=. log lcat ('dtrmmrlnn_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnn  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrlnu_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnu  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrltn_mtbli_'  tmonad (]`]`nan`nan`(trmmrltn  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrltu_mtbli_'  tmonad (]`]`nan`nan`(trmmrltu  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrunn_mtbli_'  tmonad (]`]`nan`nan`(trmmrunn  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrunu_mtbli_'  tmonad (]`]`nan`nan`(trmmrunu  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrutn_mtbli_'  tmonad (]`]`nan`nan`(trmmrutn  chk3mm)))@>"0 argsdn
  log=. log lcat ('dtrmmrutu_mtbli_'  tmonad (]`]`nan`nan`(trmmrutu  chk3mm)))@>"0 argsdn

  log=. log lcat ('ztrmmllnn_mtbli_'  tmonad (]`]`nan`nan`(trmmllnn  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmllnu_mtbli_'  tmonad (]`]`nan`nan`(trmmllnu  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmlltn_mtbli_'  tmonad (]`]`nan`nan`(trmmlltn  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmlltu_mtbli_'  tmonad (]`]`nan`nan`(trmmlltu  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmllcn_mtbli_'  tmonad (]`]`nan`nan`(trmmllcn  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmllcu_mtbli_'  tmonad (]`]`nan`nan`(trmmllcu  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmlunn_mtbli_'  tmonad (]`]`nan`nan`(trmmlunn  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmlunu_mtbli_'  tmonad (]`]`nan`nan`(trmmlunu  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmlutn_mtbli_'  tmonad (]`]`nan`nan`(trmmlutn  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmlutu_mtbli_'  tmonad (]`]`nan`nan`(trmmlutu  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmlucn_mtbli_'  tmonad (]`]`nan`nan`(trmmlucn  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmlucu_mtbli_'  tmonad (]`]`nan`nan`(trmmlucu  chk3mm)))@>"0 argsam
  log=. log lcat ('ztrmmrlnn_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnn  chk3mm)))@>"0 argsan
  log=. log lcat ('ztrmmrlnu_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnu  chk3mm)))@>"0 argsan
  log=. log lcat ('ztrmmrltn_mtbli_'  tmonad (]`]`nan`nan`(trmmrltn  chk3mm)))@>"0 argsan
  log=. log lcat ('ztrmmrltu_mtbli_'  tmonad (]`]`nan`nan`(trmmrltu  chk3mm)))@>"0 argsan
  log=. log lcat ('ztrmmrlcn_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcn  chk3mm)))@>"0 argsan
  log=. log lcat ('ztrmmrlcu_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcu  chk3mm)))@>"0 argsan
  log=. log lcat ('ztrmmrunn_mtbli_'  tmonad (]`]`nan`nan`(trmmrunn  chk3mm)))@>"0 argsan
  log=. log lcat ('ztrmmrunu_mtbli_'  tmonad (]`]`nan`nan`(trmmrunu  chk3mm)))@>"0 argsan
  log=. log lcat ('ztrmmrutn_mtbli_'  tmonad (]`]`nan`nan`(trmmrutn  chk3mm)))@>"0 argsan
  log=. log lcat ('ztrmmrutu_mtbli_'  tmonad (]`]`nan`nan`(trmmrutu  chk3mm)))@>"0 argsan
  log=. log lcat ('ztrmmrucn_mtbli_'  tmonad (]`]`nan`nan`(trmmrucn  chk3mm)))@>"0 argsan
  log=. log lcat ('ztrmmrucu_mtbli_'  tmonad (]`]`nan`nan`(trmmrucu  chk3mm)))@>"0 argsan

  NB. note: chk2mm accepts input in format compatible with
  NB.       BLIS' trmm3, so let's employ it

  log=. log lcat ('trmmllnnn_mtbli_'  tmonad (]`]`nan`nan`(trmmllnnn chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmllnnt_mtbli_'  tmonad (]`]`nan`nan`(trmmllnnt chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmllnnj_mtbli_'  tmonad (]`]`nan`nan`(trmmllnnj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmllnnc_mtbli_'  tmonad (]`]`nan`nan`(trmmllnnc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmllnun_mtbli_'  tmonad (]`]`nan`nan`(trmmllnun chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmllnut_mtbli_'  tmonad (]`]`nan`nan`(trmmllnut chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmllnuj_mtbli_'  tmonad (]`]`nan`nan`(trmmllnuj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmllnuc_mtbli_'  tmonad (]`]`nan`nan`(trmmllnuc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlltnn_mtbli_'  tmonad (]`]`nan`nan`(trmmlltnn chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlltnt_mtbli_'  tmonad (]`]`nan`nan`(trmmlltnt chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlltnj_mtbli_'  tmonad (]`]`nan`nan`(trmmlltnj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlltnc_mtbli_'  tmonad (]`]`nan`nan`(trmmlltnc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlltun_mtbli_'  tmonad (]`]`nan`nan`(trmmlltun chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlltut_mtbli_'  tmonad (]`]`nan`nan`(trmmlltut chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlltuj_mtbli_'  tmonad (]`]`nan`nan`(trmmlltuj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlltuc_mtbli_'  tmonad (]`]`nan`nan`(trmmlltuc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlljnn_mtbli_'  tmonad (]`]`nan`nan`(trmmlljnn chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlljnt_mtbli_'  tmonad (]`]`nan`nan`(trmmlljnt chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlljnj_mtbli_'  tmonad (]`]`nan`nan`(trmmlljnj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlljnc_mtbli_'  tmonad (]`]`nan`nan`(trmmlljnc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlljun_mtbli_'  tmonad (]`]`nan`nan`(trmmlljun chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlljut_mtbli_'  tmonad (]`]`nan`nan`(trmmlljut chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlljuj_mtbli_'  tmonad (]`]`nan`nan`(trmmlljuj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlljuc_mtbli_'  tmonad (]`]`nan`nan`(trmmlljuc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmllcnn_mtbli_'  tmonad (]`]`nan`nan`(trmmllcnn chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmllcnt_mtbli_'  tmonad (]`]`nan`nan`(trmmllcnt chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmllcnj_mtbli_'  tmonad (]`]`nan`nan`(trmmllcnj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmllcnc_mtbli_'  tmonad (]`]`nan`nan`(trmmllcnc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmllcun_mtbli_'  tmonad (]`]`nan`nan`(trmmllcun chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmllcut_mtbli_'  tmonad (]`]`nan`nan`(trmmllcut chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmllcuj_mtbli_'  tmonad (]`]`nan`nan`(trmmllcuj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmllcuc_mtbli_'  tmonad (]`]`nan`nan`(trmmllcuc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlunnn_mtbli_'  tmonad (]`]`nan`nan`(trmmlunnn chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlunnt_mtbli_'  tmonad (]`]`nan`nan`(trmmlunnt chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlunnj_mtbli_'  tmonad (]`]`nan`nan`(trmmlunnj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlunnc_mtbli_'  tmonad (]`]`nan`nan`(trmmlunnc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlunun_mtbli_'  tmonad (]`]`nan`nan`(trmmlunun chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlunut_mtbli_'  tmonad (]`]`nan`nan`(trmmlunut chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlunuj_mtbli_'  tmonad (]`]`nan`nan`(trmmlunuj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlunuc_mtbli_'  tmonad (]`]`nan`nan`(trmmlunuc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlutnn_mtbli_'  tmonad (]`]`nan`nan`(trmmlutnn chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlutnt_mtbli_'  tmonad (]`]`nan`nan`(trmmlutnt chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlutnj_mtbli_'  tmonad (]`]`nan`nan`(trmmlutnj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlutnc_mtbli_'  tmonad (]`]`nan`nan`(trmmlutnc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlutun_mtbli_'  tmonad (]`]`nan`nan`(trmmlutun chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlutut_mtbli_'  tmonad (]`]`nan`nan`(trmmlutut chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlutuj_mtbli_'  tmonad (]`]`nan`nan`(trmmlutuj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlutuc_mtbli_'  tmonad (]`]`nan`nan`(trmmlutuc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlujnn_mtbli_'  tmonad (]`]`nan`nan`(trmmlujnn chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlujnt_mtbli_'  tmonad (]`]`nan`nan`(trmmlujnt chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlujnj_mtbli_'  tmonad (]`]`nan`nan`(trmmlujnj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlujnc_mtbli_'  tmonad (]`]`nan`nan`(trmmlujnc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlujun_mtbli_'  tmonad (]`]`nan`nan`(trmmlujun chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlujut_mtbli_'  tmonad (]`]`nan`nan`(trmmlujut chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlujuj_mtbli_'  tmonad (]`]`nan`nan`(trmmlujuj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlujuc_mtbli_'  tmonad (]`]`nan`nan`(trmmlujuc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlucnn_mtbli_'  tmonad (]`]`nan`nan`(trmmlucnn chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlucnt_mtbli_'  tmonad (]`]`nan`nan`(trmmlucnt chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlucnj_mtbli_'  tmonad (]`]`nan`nan`(trmmlucnj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlucnc_mtbli_'  tmonad (]`]`nan`nan`(trmmlucnc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlucun_mtbli_'  tmonad (]`]`nan`nan`(trmmlucun chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlucut_mtbli_'  tmonad (]`]`nan`nan`(trmmlucut chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmlucuj_mtbli_'  tmonad (]`]`nan`nan`(trmmlucuj chk2mm)))@>"0 argsamm
  log=. log lcat ('trmmlucuc_mtbli_'  tmonad (]`]`nan`nan`(trmmlucuc chk2mm)))@>"0 argsamn
  log=. log lcat ('trmmrlnnn_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnnn chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrlnnt_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnnt chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrlnnj_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnnj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrlnnc_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnnc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrlnun_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnun chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrlnut_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnut chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrlnuj_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnuj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrlnuc_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnuc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrltnn_mtbli_'  tmonad (]`]`nan`nan`(trmmrltnn chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrltnt_mtbli_'  tmonad (]`]`nan`nan`(trmmrltnt chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrltnj_mtbli_'  tmonad (]`]`nan`nan`(trmmrltnj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrltnc_mtbli_'  tmonad (]`]`nan`nan`(trmmrltnc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrltun_mtbli_'  tmonad (]`]`nan`nan`(trmmrltun chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrltut_mtbli_'  tmonad (]`]`nan`nan`(trmmrltut chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrltuj_mtbli_'  tmonad (]`]`nan`nan`(trmmrltuj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrltuc_mtbli_'  tmonad (]`]`nan`nan`(trmmrltuc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrljnn_mtbli_'  tmonad (]`]`nan`nan`(trmmrljnn chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrljnt_mtbli_'  tmonad (]`]`nan`nan`(trmmrljnt chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrljnj_mtbli_'  tmonad (]`]`nan`nan`(trmmrljnj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrljnc_mtbli_'  tmonad (]`]`nan`nan`(trmmrljnc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrljun_mtbli_'  tmonad (]`]`nan`nan`(trmmrljun chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrljut_mtbli_'  tmonad (]`]`nan`nan`(trmmrljut chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrljuj_mtbli_'  tmonad (]`]`nan`nan`(trmmrljuj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrljuc_mtbli_'  tmonad (]`]`nan`nan`(trmmrljuc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrlcnn_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcnn chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrlcnt_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcnt chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrlcnj_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcnj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrlcnc_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcnc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrlcun_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcun chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrlcut_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcut chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrlcuj_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcuj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrlcuc_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcuc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrunnn_mtbli_'  tmonad (]`]`nan`nan`(trmmrunnn chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrunnt_mtbli_'  tmonad (]`]`nan`nan`(trmmrunnt chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrunnj_mtbli_'  tmonad (]`]`nan`nan`(trmmrunnj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrunnc_mtbli_'  tmonad (]`]`nan`nan`(trmmrunnc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrunun_mtbli_'  tmonad (]`]`nan`nan`(trmmrunun chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrunut_mtbli_'  tmonad (]`]`nan`nan`(trmmrunut chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrunuj_mtbli_'  tmonad (]`]`nan`nan`(trmmrunuj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrunuc_mtbli_'  tmonad (]`]`nan`nan`(trmmrunuc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrutnn_mtbli_'  tmonad (]`]`nan`nan`(trmmrutnn chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrutnt_mtbli_'  tmonad (]`]`nan`nan`(trmmrutnt chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrutnj_mtbli_'  tmonad (]`]`nan`nan`(trmmrutnj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrutnc_mtbli_'  tmonad (]`]`nan`nan`(trmmrutnc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrutun_mtbli_'  tmonad (]`]`nan`nan`(trmmrutun chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrutut_mtbli_'  tmonad (]`]`nan`nan`(trmmrutut chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrutuj_mtbli_'  tmonad (]`]`nan`nan`(trmmrutuj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrutuc_mtbli_'  tmonad (]`]`nan`nan`(trmmrutuc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrujnn_mtbli_'  tmonad (]`]`nan`nan`(trmmrujnn chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrujnt_mtbli_'  tmonad (]`]`nan`nan`(trmmrujnt chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrujnj_mtbli_'  tmonad (]`]`nan`nan`(trmmrujnj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrujnc_mtbli_'  tmonad (]`]`nan`nan`(trmmrujnc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrujun_mtbli_'  tmonad (]`]`nan`nan`(trmmrujun chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrujut_mtbli_'  tmonad (]`]`nan`nan`(trmmrujut chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrujuj_mtbli_'  tmonad (]`]`nan`nan`(trmmrujuj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrujuc_mtbli_'  tmonad (]`]`nan`nan`(trmmrujuc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrucnn_mtbli_'  tmonad (]`]`nan`nan`(trmmrucnn chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrucnt_mtbli_'  tmonad (]`]`nan`nan`(trmmrucnt chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrucnj_mtbli_'  tmonad (]`]`nan`nan`(trmmrucnj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrucnc_mtbli_'  tmonad (]`]`nan`nan`(trmmrucnc chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrucun_mtbli_'  tmonad (]`]`nan`nan`(trmmrucun chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrucut_mtbli_'  tmonad (]`]`nan`nan`(trmmrucut chk2mm)))@>"0 argsann
  log=. log lcat ('trmmrucuj_mtbli_'  tmonad (]`]`nan`nan`(trmmrucuj chk2mm)))@>"0 argsanm
  log=. log lcat ('trmmrucuc_mtbli_'  tmonad (]`]`nan`nan`(trmmrucuc chk2mm)))@>"0 argsann

  log=. log lcat ('dtrmmllnnn_mtbli_' tmonad (]`]`nan`nan`(trmmllnnn chk2mm)))@>"0 argsdmm
  log=. log lcat ('dtrmmllnnt_mtbli_' tmonad (]`]`nan`nan`(trmmllnnt chk2mm)))@>"0 argsdmn
  log=. log lcat ('dtrmmllnun_mtbli_' tmonad (]`]`nan`nan`(trmmllnun chk2mm)))@>"0 argsdmm
  log=. log lcat ('dtrmmllnut_mtbli_' tmonad (]`]`nan`nan`(trmmllnut chk2mm)))@>"0 argsdmn
  log=. log lcat ('dtrmmlltnn_mtbli_' tmonad (]`]`nan`nan`(trmmlltnn chk2mm)))@>"0 argsdmm
  log=. log lcat ('dtrmmlltnt_mtbli_' tmonad (]`]`nan`nan`(trmmlltnt chk2mm)))@>"0 argsdmn
  log=. log lcat ('dtrmmlltun_mtbli_' tmonad (]`]`nan`nan`(trmmlltun chk2mm)))@>"0 argsdmm
  log=. log lcat ('dtrmmlltut_mtbli_' tmonad (]`]`nan`nan`(trmmlltut chk2mm)))@>"0 argsdmn
  log=. log lcat ('dtrmmlunnn_mtbli_' tmonad (]`]`nan`nan`(trmmlunnn chk2mm)))@>"0 argsdmm
  log=. log lcat ('dtrmmlunnt_mtbli_' tmonad (]`]`nan`nan`(trmmlunnt chk2mm)))@>"0 argsdmn
  log=. log lcat ('dtrmmlunun_mtbli_' tmonad (]`]`nan`nan`(trmmlunun chk2mm)))@>"0 argsdmm
  log=. log lcat ('dtrmmlunut_mtbli_' tmonad (]`]`nan`nan`(trmmlunut chk2mm)))@>"0 argsdmn
  log=. log lcat ('dtrmmlutnn_mtbli_' tmonad (]`]`nan`nan`(trmmlutnn chk2mm)))@>"0 argsdmm
  log=. log lcat ('dtrmmlutnt_mtbli_' tmonad (]`]`nan`nan`(trmmlutnt chk2mm)))@>"0 argsdmn
  log=. log lcat ('dtrmmlutun_mtbli_' tmonad (]`]`nan`nan`(trmmlutun chk2mm)))@>"0 argsdmm
  log=. log lcat ('dtrmmlutut_mtbli_' tmonad (]`]`nan`nan`(trmmlutut chk2mm)))@>"0 argsdmn
  log=. log lcat ('dtrmmrlnnn_mtbli_' tmonad (]`]`nan`nan`(trmmrlnnn chk2mm)))@>"0 argsdnm
  log=. log lcat ('dtrmmrlnnt_mtbli_' tmonad (]`]`nan`nan`(trmmrlnnt chk2mm)))@>"0 argsdnn
  log=. log lcat ('dtrmmrlnun_mtbli_' tmonad (]`]`nan`nan`(trmmrlnun chk2mm)))@>"0 argsdnm
  log=. log lcat ('dtrmmrlnut_mtbli_' tmonad (]`]`nan`nan`(trmmrlnut chk2mm)))@>"0 argsdnn
  log=. log lcat ('dtrmmrltnn_mtbli_' tmonad (]`]`nan`nan`(trmmrltnn chk2mm)))@>"0 argsdnm
  log=. log lcat ('dtrmmrltnt_mtbli_' tmonad (]`]`nan`nan`(trmmrltnt chk2mm)))@>"0 argsdnn
  log=. log lcat ('dtrmmrltun_mtbli_' tmonad (]`]`nan`nan`(trmmrltun chk2mm)))@>"0 argsdnm
  log=. log lcat ('dtrmmrltut_mtbli_' tmonad (]`]`nan`nan`(trmmrltut chk2mm)))@>"0 argsdnn
  log=. log lcat ('dtrmmrunnn_mtbli_' tmonad (]`]`nan`nan`(trmmrunnn chk2mm)))@>"0 argsdnm
  log=. log lcat ('dtrmmrunnt_mtbli_' tmonad (]`]`nan`nan`(trmmrunnt chk2mm)))@>"0 argsdnn
  log=. log lcat ('dtrmmrunun_mtbli_' tmonad (]`]`nan`nan`(trmmrunun chk2mm)))@>"0 argsdnm
  log=. log lcat ('dtrmmrunut_mtbli_' tmonad (]`]`nan`nan`(trmmrunut chk2mm)))@>"0 argsdnn
  log=. log lcat ('dtrmmrutnn_mtbli_' tmonad (]`]`nan`nan`(trmmrutnn chk2mm)))@>"0 argsdnm
  log=. log lcat ('dtrmmrutnt_mtbli_' tmonad (]`]`nan`nan`(trmmrutnt chk2mm)))@>"0 argsdnn
  log=. log lcat ('dtrmmrutun_mtbli_' tmonad (]`]`nan`nan`(trmmrutun chk2mm)))@>"0 argsdnm
  log=. log lcat ('dtrmmrutut_mtbli_' tmonad (]`]`nan`nan`(trmmrutut chk2mm)))@>"0 argsdnn

  log=. log lcat ('ztrmmllnnn_mtbli_' tmonad (]`]`nan`nan`(trmmllnnn chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmllnnt_mtbli_' tmonad (]`]`nan`nan`(trmmllnnt chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmllnnj_mtbli_' tmonad (]`]`nan`nan`(trmmllnnj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmllnnc_mtbli_' tmonad (]`]`nan`nan`(trmmllnnc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmllnun_mtbli_' tmonad (]`]`nan`nan`(trmmllnun chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmllnut_mtbli_' tmonad (]`]`nan`nan`(trmmllnut chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmllnuj_mtbli_' tmonad (]`]`nan`nan`(trmmllnuj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmllnuc_mtbli_' tmonad (]`]`nan`nan`(trmmllnuc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlltnn_mtbli_' tmonad (]`]`nan`nan`(trmmlltnn chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlltnt_mtbli_' tmonad (]`]`nan`nan`(trmmlltnt chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlltnj_mtbli_' tmonad (]`]`nan`nan`(trmmlltnj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlltnc_mtbli_' tmonad (]`]`nan`nan`(trmmlltnc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlltun_mtbli_' tmonad (]`]`nan`nan`(trmmlltun chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlltut_mtbli_' tmonad (]`]`nan`nan`(trmmlltut chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlltuj_mtbli_' tmonad (]`]`nan`nan`(trmmlltuj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlltuc_mtbli_' tmonad (]`]`nan`nan`(trmmlltuc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlljnn_mtbli_' tmonad (]`]`nan`nan`(trmmlljnn chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlljnt_mtbli_' tmonad (]`]`nan`nan`(trmmlljnt chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlljnj_mtbli_' tmonad (]`]`nan`nan`(trmmlljnj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlljnc_mtbli_' tmonad (]`]`nan`nan`(trmmlljnc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlljun_mtbli_' tmonad (]`]`nan`nan`(trmmlljun chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlljut_mtbli_' tmonad (]`]`nan`nan`(trmmlljut chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlljuj_mtbli_' tmonad (]`]`nan`nan`(trmmlljuj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlljuc_mtbli_' tmonad (]`]`nan`nan`(trmmlljuc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmllcnn_mtbli_' tmonad (]`]`nan`nan`(trmmllcnn chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmllcnt_mtbli_' tmonad (]`]`nan`nan`(trmmllcnt chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmllcnj_mtbli_' tmonad (]`]`nan`nan`(trmmllcnj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmllcnc_mtbli_' tmonad (]`]`nan`nan`(trmmllcnc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmllcun_mtbli_' tmonad (]`]`nan`nan`(trmmllcun chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmllcut_mtbli_' tmonad (]`]`nan`nan`(trmmllcut chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmllcuj_mtbli_' tmonad (]`]`nan`nan`(trmmllcuj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmllcuc_mtbli_' tmonad (]`]`nan`nan`(trmmllcuc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlunnn_mtbli_' tmonad (]`]`nan`nan`(trmmlunnn chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlunnt_mtbli_' tmonad (]`]`nan`nan`(trmmlunnt chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlunnj_mtbli_' tmonad (]`]`nan`nan`(trmmlunnj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlunnc_mtbli_' tmonad (]`]`nan`nan`(trmmlunnc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlunun_mtbli_' tmonad (]`]`nan`nan`(trmmlunun chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlunut_mtbli_' tmonad (]`]`nan`nan`(trmmlunut chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlunuj_mtbli_' tmonad (]`]`nan`nan`(trmmlunuj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlunuc_mtbli_' tmonad (]`]`nan`nan`(trmmlunuc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlutnn_mtbli_' tmonad (]`]`nan`nan`(trmmlutnn chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlutnt_mtbli_' tmonad (]`]`nan`nan`(trmmlutnt chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlutnj_mtbli_' tmonad (]`]`nan`nan`(trmmlutnj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlutnc_mtbli_' tmonad (]`]`nan`nan`(trmmlutnc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlutun_mtbli_' tmonad (]`]`nan`nan`(trmmlutun chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlutut_mtbli_' tmonad (]`]`nan`nan`(trmmlutut chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlutuj_mtbli_' tmonad (]`]`nan`nan`(trmmlutuj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlutuc_mtbli_' tmonad (]`]`nan`nan`(trmmlutuc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlujnn_mtbli_' tmonad (]`]`nan`nan`(trmmlujnn chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlujnt_mtbli_' tmonad (]`]`nan`nan`(trmmlujnt chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlujnj_mtbli_' tmonad (]`]`nan`nan`(trmmlujnj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlujnc_mtbli_' tmonad (]`]`nan`nan`(trmmlujnc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlujun_mtbli_' tmonad (]`]`nan`nan`(trmmlujun chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlujut_mtbli_' tmonad (]`]`nan`nan`(trmmlujut chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlujuj_mtbli_' tmonad (]`]`nan`nan`(trmmlujuj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlujuc_mtbli_' tmonad (]`]`nan`nan`(trmmlujuc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlucnn_mtbli_' tmonad (]`]`nan`nan`(trmmlucnn chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlucnt_mtbli_' tmonad (]`]`nan`nan`(trmmlucnt chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlucnj_mtbli_' tmonad (]`]`nan`nan`(trmmlucnj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlucnc_mtbli_' tmonad (]`]`nan`nan`(trmmlucnc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlucun_mtbli_' tmonad (]`]`nan`nan`(trmmlucun chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlucut_mtbli_' tmonad (]`]`nan`nan`(trmmlucut chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmlucuj_mtbli_' tmonad (]`]`nan`nan`(trmmlucuj chk2mm)))@>"0 argszmm
  log=. log lcat ('ztrmmlucuc_mtbli_' tmonad (]`]`nan`nan`(trmmlucuc chk2mm)))@>"0 argszmn
  log=. log lcat ('ztrmmrlnnn_mtbli_' tmonad (]`]`nan`nan`(trmmrlnnn chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrlnnt_mtbli_' tmonad (]`]`nan`nan`(trmmrlnnt chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrlnnj_mtbli_' tmonad (]`]`nan`nan`(trmmrlnnj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrlnnc_mtbli_' tmonad (]`]`nan`nan`(trmmrlnnc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrlnun_mtbli_' tmonad (]`]`nan`nan`(trmmrlnun chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrlnut_mtbli_' tmonad (]`]`nan`nan`(trmmrlnut chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrlnuj_mtbli_' tmonad (]`]`nan`nan`(trmmrlnuj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrlnuc_mtbli_' tmonad (]`]`nan`nan`(trmmrlnuc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrltnn_mtbli_' tmonad (]`]`nan`nan`(trmmrltnn chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrltnt_mtbli_' tmonad (]`]`nan`nan`(trmmrltnt chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrltnj_mtbli_' tmonad (]`]`nan`nan`(trmmrltnj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrltnc_mtbli_' tmonad (]`]`nan`nan`(trmmrltnc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrltun_mtbli_' tmonad (]`]`nan`nan`(trmmrltun chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrltut_mtbli_' tmonad (]`]`nan`nan`(trmmrltut chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrltuj_mtbli_' tmonad (]`]`nan`nan`(trmmrltuj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrltuc_mtbli_' tmonad (]`]`nan`nan`(trmmrltuc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrljnn_mtbli_' tmonad (]`]`nan`nan`(trmmrljnn chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrljnt_mtbli_' tmonad (]`]`nan`nan`(trmmrljnt chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrljnj_mtbli_' tmonad (]`]`nan`nan`(trmmrljnj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrljnc_mtbli_' tmonad (]`]`nan`nan`(trmmrljnc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrljun_mtbli_' tmonad (]`]`nan`nan`(trmmrljun chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrljut_mtbli_' tmonad (]`]`nan`nan`(trmmrljut chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrljuj_mtbli_' tmonad (]`]`nan`nan`(trmmrljuj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrljuc_mtbli_' tmonad (]`]`nan`nan`(trmmrljuc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrlcnn_mtbli_' tmonad (]`]`nan`nan`(trmmrlcnn chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrlcnt_mtbli_' tmonad (]`]`nan`nan`(trmmrlcnt chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrlcnj_mtbli_' tmonad (]`]`nan`nan`(trmmrlcnj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrlcnc_mtbli_' tmonad (]`]`nan`nan`(trmmrlcnc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrlcun_mtbli_' tmonad (]`]`nan`nan`(trmmrlcun chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrlcut_mtbli_' tmonad (]`]`nan`nan`(trmmrlcut chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrlcuj_mtbli_' tmonad (]`]`nan`nan`(trmmrlcuj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrlcuc_mtbli_' tmonad (]`]`nan`nan`(trmmrlcuc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrunnn_mtbli_' tmonad (]`]`nan`nan`(trmmrunnn chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrunnt_mtbli_' tmonad (]`]`nan`nan`(trmmrunnt chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrunnj_mtbli_' tmonad (]`]`nan`nan`(trmmrunnj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrunnc_mtbli_' tmonad (]`]`nan`nan`(trmmrunnc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrunun_mtbli_' tmonad (]`]`nan`nan`(trmmrunun chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrunut_mtbli_' tmonad (]`]`nan`nan`(trmmrunut chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrunuj_mtbli_' tmonad (]`]`nan`nan`(trmmrunuj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrunuc_mtbli_' tmonad (]`]`nan`nan`(trmmrunuc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrutnn_mtbli_' tmonad (]`]`nan`nan`(trmmrutnn chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrutnt_mtbli_' tmonad (]`]`nan`nan`(trmmrutnt chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrutnj_mtbli_' tmonad (]`]`nan`nan`(trmmrutnj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrutnc_mtbli_' tmonad (]`]`nan`nan`(trmmrutnc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrutun_mtbli_' tmonad (]`]`nan`nan`(trmmrutun chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrutut_mtbli_' tmonad (]`]`nan`nan`(trmmrutut chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrutuj_mtbli_' tmonad (]`]`nan`nan`(trmmrutuj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrutuc_mtbli_' tmonad (]`]`nan`nan`(trmmrutuc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrujnn_mtbli_' tmonad (]`]`nan`nan`(trmmrujnn chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrujnt_mtbli_' tmonad (]`]`nan`nan`(trmmrujnt chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrujnj_mtbli_' tmonad (]`]`nan`nan`(trmmrujnj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrujnc_mtbli_' tmonad (]`]`nan`nan`(trmmrujnc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrujun_mtbli_' tmonad (]`]`nan`nan`(trmmrujun chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrujut_mtbli_' tmonad (]`]`nan`nan`(trmmrujut chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrujuj_mtbli_' tmonad (]`]`nan`nan`(trmmrujuj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrujuc_mtbli_' tmonad (]`]`nan`nan`(trmmrujuc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrucnn_mtbli_' tmonad (]`]`nan`nan`(trmmrucnn chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrucnt_mtbli_' tmonad (]`]`nan`nan`(trmmrucnt chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrucnj_mtbli_' tmonad (]`]`nan`nan`(trmmrucnj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrucnc_mtbli_' tmonad (]`]`nan`nan`(trmmrucnc chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrucun_mtbli_' tmonad (]`]`nan`nan`(trmmrucun chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrucut_mtbli_' tmonad (]`]`nan`nan`(trmmrucut chk2mm)))@>"0 argsznn
  log=. log lcat ('ztrmmrucuj_mtbli_' tmonad (]`]`nan`nan`(trmmrucuj chk2mm)))@>"0 argsznm
  log=. log lcat ('ztrmmrucuc_mtbli_' tmonad (]`]`nan`nan`(trmmrucuc chk2mm)))@>"0 argsznn
)

NB. ---------------------------------------------------------
NB. testbasicmm
NB.
NB. Description:
NB.   Adv. to make verb to test basic matrix-matrix operations
NB.
NB. Syntax:
NB.   log=. (mkmat testbasicmm) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testbasicmm_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testbasicmm_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testbasicmm_mt_ 150 200

testbasicmm=: 1 : 'testbasictrmm_mt_@(u@(2 # >./) ; u ; u) ,&.>~ testbasichemm_mt_@((9&o. upddiag_mt_)@u@(2 # >./) ; u ; u) ,&.>~ testbasicsymm_mt_@(u@(2 # >./) ; u ; u) ,&.>~ testbasicgemmt_mt_@(u@(+/\) ; u@|.@(+/\) ; u@(2 # {.)) ,&.>~ testbasicgemm_mt_@(u@(+/\) ; u@(+/\.) ; u) [ require@''math/mt/external/blis/mm'' [ require@''math/mt/external/blas/mm'''

NB. ---------------------------------------------------------
NB. testbasictrsv
NB.
NB. Description:
NB.   Test equation solvers:
NB.   - xTRSV (BLAS)
NB.   by triangular matrix and vector
NB.
NB. Syntax:
NB.   log=. testbasictrsv AA ; b
NB. where
NB.   AA  - n×n-matrix, A material
NB.   x   - n-vector, the RHS
NB.   log - 6-vector of boxes, test log, see test.ijs

testbasictrsv=: 3 : 0
  inc=. 1 2 _1 _2
  'AA y'=. y
  args=. (1 reverse 2)@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc  NB. 4×3-matrix of boxes, each row is argument to tmonad

  NB. for every i feed the tuple (AA ; expanded_b_i ; incb_i) to tmonad
  log=.          ('dtrsvlnn_mtbla_' tmonad (]`]`nan`nan`(trmvlnn chk3sv)))"1 args
  log=. log lcat ('dtrsvlnu_mtbla_' tmonad (]`]`nan`nan`(trmvlnu chk3sv)))"1 args
  log=. log lcat ('dtrsvltn_mtbla_' tmonad (]`]`nan`nan`(trmvltn chk3sv)))"1 args
  log=. log lcat ('dtrsvltu_mtbla_' tmonad (]`]`nan`nan`(trmvltu chk3sv)))"1 args
  log=. log lcat ('dtrsvunn_mtbla_' tmonad (]`]`nan`nan`(trmvunn chk3sv)))"1 args
  log=. log lcat ('dtrsvunu_mtbla_' tmonad (]`]`nan`nan`(trmvunu chk3sv)))"1 args
  log=. log lcat ('dtrsvutn_mtbla_' tmonad (]`]`nan`nan`(trmvutn chk3sv)))"1 args
  log=. log lcat ('dtrsvutu_mtbla_' tmonad (]`]`nan`nan`(trmvutu chk3sv)))"1 args
  log=. log lcat ('ztrsvlnn_mtbla_' tmonad (]`]`nan`nan`(trmvlnn chk3sv)))"1 args
  log=. log lcat ('ztrsvlnu_mtbla_' tmonad (]`]`nan`nan`(trmvlnu chk3sv)))"1 args
  log=. log lcat ('ztrsvltn_mtbla_' tmonad (]`]`nan`nan`(trmvltn chk3sv)))"1 args
  log=. log lcat ('ztrsvltu_mtbla_' tmonad (]`]`nan`nan`(trmvltu chk3sv)))"1 args
  log=. log lcat ('ztrsvlcn_mtbla_' tmonad (]`]`nan`nan`(trmvlcn chk3sv)))"1 args
  log=. log lcat ('ztrsvlcu_mtbla_' tmonad (]`]`nan`nan`(trmvlcu chk3sv)))"1 args
  log=. log lcat ('ztrsvunn_mtbla_' tmonad (]`]`nan`nan`(trmvunn chk3sv)))"1 args
  log=. log lcat ('ztrsvunu_mtbla_' tmonad (]`]`nan`nan`(trmvunu chk3sv)))"1 args
  log=. log lcat ('ztrsvutn_mtbla_' tmonad (]`]`nan`nan`(trmvutn chk3sv)))"1 args
  log=. log lcat ('ztrsvutu_mtbla_' tmonad (]`]`nan`nan`(trmvutu chk3sv)))"1 args
  log=. log lcat ('ztrsvucn_mtbla_' tmonad (]`]`nan`nan`(trmvucn chk3sv)))"1 args
  log=. log lcat ('ztrsvucu_mtbla_' tmonad (]`]`nan`nan`(trmvucu chk3sv)))"1 args
)

NB. ---------------------------------------------------------
NB. testbasicsv
NB.
NB. Description:
NB.   Adv. to make verb to test equation solvers
NB.
NB. Syntax:
NB.   log=. (mkmat testbasicsv) (n,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (n,n)
NB.   (n,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testbasicsv_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testbasicsv_mt_ 200 200
NB. - test by random square complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testbasicsv_mt_ 250 250

testbasicsv=: 1 : 'testbasictrsv_mt_@(u ; u@{.) [ require@''math/mt/external/blas/sv'''

NB. ---------------------------------------------------------
NB. testbasictrsm
NB.
NB. Description:
NB.   Test matrix equation solvers:
NB.   - trsmxxxx (math/mt addon)
NB.   - xTRSM (BLAS)
NB.   - bli_xtrsm (BLIS)
NB.   by triangular matrix
NB.
NB. Syntax:
NB.   log=. testbasictrsm AA ; B
NB. where
NB.   AA  - k×k-matrix, A material
NB.   B   - m×n-matrix, RHS
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   k   = max(m,n)
NB.
NB. Notes:
NB. - For real matrices and complex coefficients bli_xtrsm
NB.   saves the real part of result in B and imagine part
NB.   will be lost. That is why it's sufficient to test
NB.   bli_xtrsm by arguments of the same datatype only.

testbasictrsm=: 3 : 0
  'AA B'=. y
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  acoeff=. /:~ ~. zcoeff ,^:(JCMPX = 3!:0 B) dcoeff
  'm n'=. $ B
  Am=. (2 # m) {. AA
  An=. (2 # n) {. AA
  'bm bn'=. ({."1 ; {.) B
  argsdm=. { (<"0 dcoeff) ; (< Am) ; < < B
  argsdn=. { (<"0 dcoeff) ; (< An) ; < < B
  argszm=. { (<"0 zcoeff) ; (< Am) ; < < B
  argszn=. { (<"0 zcoeff) ; (< An) ; < < B
  argsam=. { (<"0 acoeff) ; (< Am) ; < < B
  argsan=. { (<"0 acoeff) ; (< An) ; < < B
  argsbm=. { (<"0 acoeff) ; (< Am) ; < < bm
  argsbn=. { (<"0 acoeff) ; (< An) ; < < bn

  NB. for every i feed the tuple (alpha_i ; Ax ; B) to tmonad

  NB. BLAS' staff

  log=.          ('dtrsmllnn_mtbla_' tmonad (]`]`nan`nan`(trmmllnn chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmllnu_mtbla_' tmonad (]`]`nan`nan`(trmmllnu chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlltn_mtbla_' tmonad (]`]`nan`nan`(trmmlltn chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlltu_mtbla_' tmonad (]`]`nan`nan`(trmmlltu chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlunn_mtbla_' tmonad (]`]`nan`nan`(trmmlunn chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlunu_mtbla_' tmonad (]`]`nan`nan`(trmmlunu chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlutn_mtbla_' tmonad (]`]`nan`nan`(trmmlutn chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlutu_mtbla_' tmonad (]`]`nan`nan`(trmmlutu chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmrlnn_mtbla_' tmonad (]`]`nan`nan`(trmmrlnn chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrlnu_mtbla_' tmonad (]`]`nan`nan`(trmmrlnu chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrltn_mtbla_' tmonad (]`]`nan`nan`(trmmrltn chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrltu_mtbla_' tmonad (]`]`nan`nan`(trmmrltu chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrunn_mtbla_' tmonad (]`]`nan`nan`(trmmrunn chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrunu_mtbla_' tmonad (]`]`nan`nan`(trmmrunu chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrutn_mtbla_' tmonad (]`]`nan`nan`(trmmrutn chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrutu_mtbla_' tmonad (]`]`nan`nan`(trmmrutu chk3sm)))@>"0 argsdn

  log=. log lcat ('ztrsmllnn_mtbla_' tmonad (]`]`nan`nan`(trmmllnn chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmllnu_mtbla_' tmonad (]`]`nan`nan`(trmmllnu chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmlltn_mtbla_' tmonad (]`]`nan`nan`(trmmlltn chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmlltu_mtbla_' tmonad (]`]`nan`nan`(trmmlltu chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmllcn_mtbla_' tmonad (]`]`nan`nan`(trmmllcn chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmllcu_mtbla_' tmonad (]`]`nan`nan`(trmmllcu chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmlunn_mtbla_' tmonad (]`]`nan`nan`(trmmlunn chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmlunu_mtbla_' tmonad (]`]`nan`nan`(trmmlunu chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmlutn_mtbla_' tmonad (]`]`nan`nan`(trmmlutn chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmlutu_mtbla_' tmonad (]`]`nan`nan`(trmmlutu chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmlucn_mtbla_' tmonad (]`]`nan`nan`(trmmlucn chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmlucu_mtbla_' tmonad (]`]`nan`nan`(trmmlucu chk3sm)))@>"0 argszm
  log=. log lcat ('ztrsmrlnn_mtbla_' tmonad (]`]`nan`nan`(trmmrlnn chk3sm)))@>"0 argszn
  log=. log lcat ('ztrsmrlnu_mtbla_' tmonad (]`]`nan`nan`(trmmrlnu chk3sm)))@>"0 argszn
  log=. log lcat ('ztrsmrltn_mtbla_' tmonad (]`]`nan`nan`(trmmrltn chk3sm)))@>"0 argszn
  log=. log lcat ('ztrsmrltu_mtbla_' tmonad (]`]`nan`nan`(trmmrltu chk3sm)))@>"0 argszn
  log=. log lcat ('ztrsmrlcn_mtbla_' tmonad (]`]`nan`nan`(trmmrlcn chk3sm)))@>"0 argszn
  log=. log lcat ('ztrsmrlcu_mtbla_' tmonad (]`]`nan`nan`(trmmrlcu chk3sm)))@>"0 argszn
  log=. log lcat ('ztrsmrunn_mtbla_' tmonad (]`]`nan`nan`(trmmrunn chk3sm)))@>"0 argszn
  log=. log lcat ('ztrsmrunu_mtbla_' tmonad (]`]`nan`nan`(trmmrunu chk3sm)))@>"0 argszn
  log=. log lcat ('ztrsmrutn_mtbla_' tmonad (]`]`nan`nan`(trmmrutn chk3sm)))@>"0 argszn
  log=. log lcat ('ztrsmrutu_mtbla_' tmonad (]`]`nan`nan`(trmmrutu chk3sm)))@>"0 argszn
  log=. log lcat ('ztrsmrucn_mtbla_' tmonad (]`]`nan`nan`(trmmrucn chk3sm)))@>"0 argszn
  log=. log lcat ('ztrsmrucu_mtbla_' tmonad (]`]`nan`nan`(trmmrucu chk3sm)))@>"0 argszn

  NB. BLIS' staff

  log=. log lcat ('trsmllnn_mtbli_'  tmonad (]`]`nan`nan`(trmmllnn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmllnu_mtbli_'  tmonad (]`]`nan`nan`(trmmllnu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlltn_mtbli_'  tmonad (]`]`nan`nan`(trmmlltn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlltu_mtbli_'  tmonad (]`]`nan`nan`(trmmlltu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlljn_mtbli_'  tmonad (]`]`nan`nan`(trmmlljn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmllju_mtbli_'  tmonad (]`]`nan`nan`(trmmllju chk3sm)))@>"0 argsam
  log=. log lcat ('trsmllcn_mtbli_'  tmonad (]`]`nan`nan`(trmmllcn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmllcu_mtbli_'  tmonad (]`]`nan`nan`(trmmllcu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlunn_mtbli_'  tmonad (]`]`nan`nan`(trmmlunn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlunu_mtbli_'  tmonad (]`]`nan`nan`(trmmlunu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlutn_mtbli_'  tmonad (]`]`nan`nan`(trmmlutn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlutu_mtbli_'  tmonad (]`]`nan`nan`(trmmlutu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlujn_mtbli_'  tmonad (]`]`nan`nan`(trmmlujn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmluju_mtbli_'  tmonad (]`]`nan`nan`(trmmluju chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlucn_mtbli_'  tmonad (]`]`nan`nan`(trmmlucn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlucu_mtbli_'  tmonad (]`]`nan`nan`(trmmlucu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmrlnn_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrlnu_mtbli_'  tmonad (]`]`nan`nan`(trmmrlnu chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrltn_mtbli_'  tmonad (]`]`nan`nan`(trmmrltn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrltu_mtbli_'  tmonad (]`]`nan`nan`(trmmrltu chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrljn_mtbli_'  tmonad (]`]`nan`nan`(trmmrljn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrlju_mtbli_'  tmonad (]`]`nan`nan`(trmmrlju chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrlcn_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrlcu_mtbli_'  tmonad (]`]`nan`nan`(trmmrlcu chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrunn_mtbli_'  tmonad (]`]`nan`nan`(trmmrunn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrunu_mtbli_'  tmonad (]`]`nan`nan`(trmmrunu chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrutn_mtbli_'  tmonad (]`]`nan`nan`(trmmrutn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrutu_mtbli_'  tmonad (]`]`nan`nan`(trmmrutu chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrujn_mtbli_'  tmonad (]`]`nan`nan`(trmmrujn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmruju_mtbli_'  tmonad (]`]`nan`nan`(trmmruju chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrucn_mtbli_'  tmonad (]`]`nan`nan`(trmmrucn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrucu_mtbli_'  tmonad (]`]`nan`nan`(trmmrucu chk3sm)))@>"0 argsan

  log=. log lcat ('dtrsmllnn_mtbli_' tmonad (]`]`nan`nan`(trmmllnn chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmllnu_mtbli_' tmonad (]`]`nan`nan`(trmmllnu chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlltn_mtbli_' tmonad (]`]`nan`nan`(trmmlltn chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlltu_mtbli_' tmonad (]`]`nan`nan`(trmmlltu chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlunn_mtbli_' tmonad (]`]`nan`nan`(trmmlunn chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlunu_mtbli_' tmonad (]`]`nan`nan`(trmmlunu chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlutn_mtbli_' tmonad (]`]`nan`nan`(trmmlutn chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmlutu_mtbli_' tmonad (]`]`nan`nan`(trmmlutu chk3sm)))@>"0 argsdm
  log=. log lcat ('dtrsmrlnn_mtbli_' tmonad (]`]`nan`nan`(trmmrlnn chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrlnu_mtbli_' tmonad (]`]`nan`nan`(trmmrlnu chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrltn_mtbli_' tmonad (]`]`nan`nan`(trmmrltn chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrltu_mtbli_' tmonad (]`]`nan`nan`(trmmrltu chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrunn_mtbli_' tmonad (]`]`nan`nan`(trmmrunn chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrunu_mtbli_' tmonad (]`]`nan`nan`(trmmrunu chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrutn_mtbli_' tmonad (]`]`nan`nan`(trmmrutn chk3sm)))@>"0 argsdn
  log=. log lcat ('dtrsmrutu_mtbli_' tmonad (]`]`nan`nan`(trmmrutu chk3sm)))@>"0 argsdn

  log=. log lcat ('ztrsmllnn_mtbli_' tmonad (]`]`nan`nan`(trmmllnn chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmllnu_mtbli_' tmonad (]`]`nan`nan`(trmmllnu chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmlltn_mtbli_' tmonad (]`]`nan`nan`(trmmlltn chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmlltu_mtbli_' tmonad (]`]`nan`nan`(trmmlltu chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmlljn_mtbli_' tmonad (]`]`nan`nan`(trmmlljn chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmllju_mtbli_' tmonad (]`]`nan`nan`(trmmllju chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmllcn_mtbli_' tmonad (]`]`nan`nan`(trmmllcn chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmllcu_mtbli_' tmonad (]`]`nan`nan`(trmmllcu chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmlunn_mtbli_' tmonad (]`]`nan`nan`(trmmlunn chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmlunu_mtbli_' tmonad (]`]`nan`nan`(trmmlunu chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmlutn_mtbli_' tmonad (]`]`nan`nan`(trmmlutn chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmlutu_mtbli_' tmonad (]`]`nan`nan`(trmmlutu chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmlujn_mtbli_' tmonad (]`]`nan`nan`(trmmlujn chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmluju_mtbli_' tmonad (]`]`nan`nan`(trmmluju chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmlucn_mtbli_' tmonad (]`]`nan`nan`(trmmlucn chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmlucu_mtbli_' tmonad (]`]`nan`nan`(trmmlucu chk3sm)))@>"0 argsam
  log=. log lcat ('ztrsmrlnn_mtbli_' tmonad (]`]`nan`nan`(trmmrlnn chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrlnu_mtbli_' tmonad (]`]`nan`nan`(trmmrlnu chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrltn_mtbli_' tmonad (]`]`nan`nan`(trmmrltn chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrltu_mtbli_' tmonad (]`]`nan`nan`(trmmrltu chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrljn_mtbli_' tmonad (]`]`nan`nan`(trmmrljn chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrlju_mtbli_' tmonad (]`]`nan`nan`(trmmrlju chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrlcn_mtbli_' tmonad (]`]`nan`nan`(trmmrlcn chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrlcu_mtbli_' tmonad (]`]`nan`nan`(trmmrlcu chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrunn_mtbli_' tmonad (]`]`nan`nan`(trmmrunn chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrunu_mtbli_' tmonad (]`]`nan`nan`(trmmrunu chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrutn_mtbli_' tmonad (]`]`nan`nan`(trmmrutn chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrutu_mtbli_' tmonad (]`]`nan`nan`(trmmrutu chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrujn_mtbli_' tmonad (]`]`nan`nan`(trmmrujn chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmruju_mtbli_' tmonad (]`]`nan`nan`(trmmruju chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrucn_mtbli_' tmonad (]`]`nan`nan`(trmmrucn chk3sm)))@>"0 argsan
  log=. log lcat ('ztrsmrucu_mtbli_' tmonad (]`]`nan`nan`(trmmrucu chk3sm)))@>"0 argsan

  NB. mt staff

  NB. monadic trsmxxxx, 1-rank b and x
  log=. log lcat ('trsmllnn'         tmonad (]`]`nan`nan`(trmmllnn chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmllnu'         tmonad (]`]`nan`nan`(trmmllnu chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmlltn'         tmonad (]`]`nan`nan`(trmmlltn chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmlltu'         tmonad (]`]`nan`nan`(trmmlltu chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmlljn'         tmonad (]`]`nan`nan`(trmmlljn chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmllju'         tmonad (]`]`nan`nan`(trmmllju chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmllcn'         tmonad (]`]`nan`nan`(trmmllcn chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmllcu'         tmonad (]`]`nan`nan`(trmmllcu chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmlunn'         tmonad (]`]`nan`nan`(trmmlunn chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmlunu'         tmonad (]`]`nan`nan`(trmmlunu chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmlutn'         tmonad (]`]`nan`nan`(trmmlutn chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmlutu'         tmonad (]`]`nan`nan`(trmmlutu chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmlujn'         tmonad (]`]`nan`nan`(trmmlujn chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmluju'         tmonad (]`]`nan`nan`(trmmluju chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmlucn'         tmonad (]`]`nan`nan`(trmmlucn chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmlucu'         tmonad (]`]`nan`nan`(trmmlucu chk3sm)))@>"0 argsbm
  log=. log lcat ('trsmrlnn'         tmonad (]`]`nan`nan`(trmmrlnn chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrlnu'         tmonad (]`]`nan`nan`(trmmrlnu chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrltn'         tmonad (]`]`nan`nan`(trmmrltn chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrltu'         tmonad (]`]`nan`nan`(trmmrltu chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrljn'         tmonad (]`]`nan`nan`(trmmrljn chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrlju'         tmonad (]`]`nan`nan`(trmmrlju chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrlcn'         tmonad (]`]`nan`nan`(trmmrlcn chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrlcu'         tmonad (]`]`nan`nan`(trmmrlcu chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrunn'         tmonad (]`]`nan`nan`(trmmrunn chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrunu'         tmonad (]`]`nan`nan`(trmmrunu chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrutn'         tmonad (]`]`nan`nan`(trmmrutn chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrutu'         tmonad (]`]`nan`nan`(trmmrutu chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrujn'         tmonad (]`]`nan`nan`(trmmrujn chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmruju'         tmonad (]`]`nan`nan`(trmmruju chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrucn'         tmonad (]`]`nan`nan`(trmmrucn chk3sm)))@>"0 argsbn
  log=. log lcat ('trsmrucu'         tmonad (]`]`nan`nan`(trmmrucu chk3sm)))@>"0 argsbn

  NB. monadic trsmxxxx, 2-rank B and X
  log=. log lcat ('trsmllnn'         tmonad (]`]`nan`nan`(trmmllnn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmllnu'         tmonad (]`]`nan`nan`(trmmllnu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlltn'         tmonad (]`]`nan`nan`(trmmlltn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlltu'         tmonad (]`]`nan`nan`(trmmlltu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlljn'         tmonad (]`]`nan`nan`(trmmlljn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmllju'         tmonad (]`]`nan`nan`(trmmllju chk3sm)))@>"0 argsam
  log=. log lcat ('trsmllcn'         tmonad (]`]`nan`nan`(trmmllcn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmllcu'         tmonad (]`]`nan`nan`(trmmllcu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlunn'         tmonad (]`]`nan`nan`(trmmlunn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlunu'         tmonad (]`]`nan`nan`(trmmlunu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlutn'         tmonad (]`]`nan`nan`(trmmlutn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlutu'         tmonad (]`]`nan`nan`(trmmlutu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlujn'         tmonad (]`]`nan`nan`(trmmlujn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmluju'         tmonad (]`]`nan`nan`(trmmluju chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlucn'         tmonad (]`]`nan`nan`(trmmlucn chk3sm)))@>"0 argsam
  log=. log lcat ('trsmlucu'         tmonad (]`]`nan`nan`(trmmlucu chk3sm)))@>"0 argsam
  log=. log lcat ('trsmrlnn'         tmonad (]`]`nan`nan`(trmmrlnn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrlnu'         tmonad (]`]`nan`nan`(trmmrlnu chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrltn'         tmonad (]`]`nan`nan`(trmmrltn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrltu'         tmonad (]`]`nan`nan`(trmmrltu chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrljn'         tmonad (]`]`nan`nan`(trmmrljn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrlju'         tmonad (]`]`nan`nan`(trmmrlju chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrlcn'         tmonad (]`]`nan`nan`(trmmrlcn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrlcu'         tmonad (]`]`nan`nan`(trmmrlcu chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrunn'         tmonad (]`]`nan`nan`(trmmrunn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrunu'         tmonad (]`]`nan`nan`(trmmrunu chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrutn'         tmonad (]`]`nan`nan`(trmmrutn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrutu'         tmonad (]`]`nan`nan`(trmmrutu chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrujn'         tmonad (]`]`nan`nan`(trmmrujn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmruju'         tmonad (]`]`nan`nan`(trmmruju chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrucn'         tmonad (]`]`nan`nan`(trmmrucn chk3sm)))@>"0 argsan
  log=. log lcat ('trsmrucu'         tmonad (]`]`nan`nan`(trmmrucu chk3sm)))@>"0 argsan

  NB. dyadic trsmxxxx
  NB. note:
  NB. - t02x expect in (0{::x) a matrix with opposite
  NB.   triangle zeroed and ignore element in (3{::x), withal
  NB.   trsmxxxx should be tested by a matrix with a noise in
  NB.   opposite (not referenced) triangle, so we send A in
  NB.   (3{::x) and extract it by vgetx as the left argument
  NB.   for trsmxxxx

  'norm1Lm  normiLm'=.  (norm1 , normi) Lm=.           trlpick Am
  'norm1Um  normiUm'=.  (norm1 , normi) Um=.           trupick Am
  'norm1Ln  normiLn'=.  (norm1 , normi) Ln=.           trlpick An
  'norm1Un  normiUn'=.  (norm1 , normi) Un=.           trupick An

  'norm1L1m normiL1m'=. (norm1 , normi) L1m=. (1 ; '') setdiag Lm
  'norm1U1m normiU1m'=. (norm1 , normi) U1m=. (1 ; '') setdiag Um
  'norm1L1n normiL1n'=. (norm1 , normi) L1n=. (1 ; '') setdiag Ln
  'norm1U1n normiU1n'=. (norm1 , normi) U1n=. (1 ; '') setdiag Un

  NB. 1-rank b and x
  NB. note:
  NB. - we use RHS vectors bm and bn here as solution vectors
  NB.   xm and xn for (2{::x), RHS is computed explicitely
  NB.   and is supplied in (1{::x)
  log=. log lcat ('trsmllnn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~    trlpick ) t02v))) Lm  ; (Lm   mp       bm ) ; bm ; Am ; norm1Lm
  log=. log lcat ('trsmllnu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~    trl1pick) t02v))) L1m ; (L1m  mp       bm ) ; bm ; Am ; norm1L1m
  log=. log lcat ('trsmlltn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ |:@trlpick ) t02v))) Lm  ; (Lm  (mp~ ct)~ bm ) ; bm ; Am ; normiLm
  log=. log lcat ('trsmlltu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ |:@trl1pick) t02v))) L1m ; (L1m (mp~ ct)~ bm ) ; bm ; Am ; normiL1m
  log=. log lcat ('trsmlljn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ + @trlpick ) t02v))) Lm  ; (Lm  (mp~ + )~ bm ) ; bm ; Am ; normiLm
  log=. log lcat ('trsmllju'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ + @trl1pick) t02v))) L1m ; (L1m (mp~ + )~ bm ) ; bm ; Am ; normiL1m
  log=. log lcat ('trsmllcn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ ct@trlpick ) t02v))) Lm  ; (Lm  (mp~ |:)~ bm ) ; bm ; Am ; normiLm
  log=. log lcat ('trsmllcu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ ct@trl1pick) t02v))) L1m ; (L1m (mp~ |:)~ bm ) ; bm ; Am ; normiL1m
  log=. log lcat ('trsmlunn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~    trupick ) t02v))) Um  ; (Um   mp       bm ) ; bm ; Am ; norm1Um
  log=. log lcat ('trsmlunu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~    tru1pick) t02v))) U1m ; (U1m  mp       bm ) ; bm ; Am ; norm1U1m
  log=. log lcat ('trsmlutn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ |:@trupick ) t02v))) Um  ; (Um  (mp~ ct)~ bm ) ; bm ; Am ; normiUm
  log=. log lcat ('trsmlutu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ |:@tru1pick) t02v))) U1m ; (U1m (mp~ ct)~ bm ) ; bm ; Am ; normiU1m
  log=. log lcat ('trsmlujn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ + @trupick ) t02v))) Um  ; (Um  (mp~ + )~ bm ) ; bm ; Am ; normiUm
  log=. log lcat ('trsmluju'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ + @tru1pick) t02v))) U1m ; (U1m (mp~ + )~ bm ) ; bm ; Am ; normiU1m
  log=. log lcat ('trsmlucn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ ct@trupick ) t02v))) Um  ; (Um  (mp~ |:)~ bm ) ; bm ; Am ; normiUm
  log=. log lcat ('trsmlucu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ ct@tru1pick) t02v))) U1m ; (U1m (mp~ |:)~ bm ) ; bm ; Am ; normiU1m
  log=. log lcat ('trsmrlnn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp     trlpick ) t02v))) Ln  ; (bn   mp       Ln ) ; bn ; An ; normiLn
  log=. log lcat ('trsmrlnu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp     trl1pick) t02v))) L1n ; (bn   mp       L1n) ; bn ; An ; normiL1n
  log=. log lcat ('trsmrltn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  |:@trlpick ) t02v))) Ln  ; (bn  (mp  ct)  Ln ) ; bn ; An ; norm1Ln
  log=. log lcat ('trsmrltu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  |:@trl1pick) t02v))) L1n ; (bn  (mp  ct)  L1n) ; bn ; An ; norm1L1n
  log=. log lcat ('trsmrljn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  + @trlpick ) t02v))) Ln  ; (bn  (mp  + )  Ln ) ; bn ; An ; norm1Ln
  log=. log lcat ('trsmrlju'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  + @trl1pick) t02v))) L1n ; (bn  (mp  + )  L1n) ; bn ; An ; norm1L1n
  log=. log lcat ('trsmrlcn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  ct@trlpick ) t02v))) Ln  ; (bn  (mp  |:)  Ln ) ; bn ; An ; norm1Ln
  log=. log lcat ('trsmrlcu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  ct@trl1pick) t02v))) L1n ; (bn  (mp  |:)  L1n) ; bn ; An ; norm1L1n
  log=. log lcat ('trsmrunn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp     trupick ) t02v))) Un  ; (bn   mp       Un ) ; bn ; An ; normiUn
  log=. log lcat ('trsmrunu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp     tru1pick) t02v))) U1n ; (bn   mp       U1n) ; bn ; An ; normiU1n
  log=. log lcat ('trsmrutn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  |:@trupick ) t02v))) Un  ; (bn  (mp  ct)  Un ) ; bn ; An ; norm1Un
  log=. log lcat ('trsmrutu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  |:@tru1pick) t02v))) U1n ; (bn  (mp  ct)  U1n) ; bn ; An ; norm1U1n
  log=. log lcat ('trsmrujn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  + @trupick ) t02v))) Un  ; (bn  (mp  + )  Un ) ; bn ; An ; norm1Un
  log=. log lcat ('trsmruju'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  + @tru1pick) t02v))) U1n ; (bn  (mp  + )  U1n) ; bn ; An ; norm1U1n
  log=. log lcat ('trsmrucn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  ct@trupick ) t02v))) Un  ; (bn  (mp  |:)  Un ) ; bn ; An ; norm1Un
  log=. log lcat ('trsmrucu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  ct@tru1pick) t02v))) U1n ; (bn  (mp  |:)  U1n) ; bn ; An ; norm1U1n

  NB. dyadic trsmxxxx, 2-rank B and X
  NB. notes:
  NB. - we use RHS matrix B here as solution vector X for
  NB.   (2{::x), RHS is computed explicitely and is supplied
  NB.   in (1{::x)
  log=. log lcat ('trsmllnn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~    trlpick ) t02m norm1tc))) Lm  ; (Lm   mp       B ) ; B  ; Am ; norm1Lm
  log=. log lcat ('trsmllnu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~    trl1pick) t02m norm1tc))) L1m ; (L1m  mp       B ) ; B  ; Am ; norm1L1m
  log=. log lcat ('trsmlltn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ |:@trlpick ) t02m norm1tc))) Lm  ; (Lm  (mp~ ct)~ B ) ; B  ; Am ; normiLm
  log=. log lcat ('trsmlltu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ |:@trl1pick) t02m norm1tc))) L1m ; (L1m (mp~ ct)~ B ) ; B  ; Am ; normiL1m
  log=. log lcat ('trsmlljn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ + @trlpick ) t02m norm1tc))) Lm  ; (Lm  (mp~ + )~ B ) ; B  ; Am ; normiLm
  log=. log lcat ('trsmllju'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ + @trl1pick) t02m norm1tc))) L1m ; (L1m (mp~ + )~ B ) ; B  ; Am ; normiL1m
  log=. log lcat ('trsmllcn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ ct@trlpick ) t02m norm1tc))) Lm  ; (Lm  (mp~ |:)~ B ) ; B  ; Am ; normiLm
  log=. log lcat ('trsmllcu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ ct@trl1pick) t02m norm1tc))) L1m ; (L1m (mp~ |:)~ B ) ; B  ; Am ; normiL1m
  log=. log lcat ('trsmlunn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~    trupick ) t02m norm1tc))) Um  ; (Um   mp       B ) ; B  ; Am ; norm1Um
  log=. log lcat ('trsmlunu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~    tru1pick) t02m norm1tc))) U1m ; (U1m  mp       B ) ; B  ; Am ; norm1U1m
  log=. log lcat ('trsmlutn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ |:@trupick ) t02m norm1tc))) Um  ; (Um  (mp~ ct)~ B ) ; B  ; Am ; normiUm
  log=. log lcat ('trsmlutu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ |:@tru1pick) t02m norm1tc))) U1m ; (U1m (mp~ ct)~ B ) ; B  ; Am ; normiU1m
  log=. log lcat ('trsmlujn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ + @trupick ) t02m norm1tc))) Um  ; (Um  (mp~ + )~ B ) ; B  ; Am ; normiUm
  log=. log lcat ('trsmluju'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ + @tru1pick) t02m norm1tc))) U1m ; (U1m (mp~ + )~ B ) ; B  ; Am ; normiU1m
  log=. log lcat ('trsmlucn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ ct@trupick ) t02m norm1tc))) Um  ; (Um  (mp~ |:)~ B ) ; B  ; Am ; normiUm
  log=. log lcat ('trsmlucu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp~ ct@tru1pick) t02m norm1tc))) U1m ; (U1m (mp~ |:)~ B ) ; B  ; Am ; normiU1m
  log=. log lcat ('trsmrlnn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp     trlpick ) t02m normitc))) Ln  ; (B   mp       Ln ) ; B  ; An ; normiLn
  log=. log lcat ('trsmrlnu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp     trl1pick) t02m normitc))) L1n ; (B   mp       L1n) ; B  ; An ; normiL1n
  log=. log lcat ('trsmrltn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  |:@trlpick ) t02m normitc))) Ln  ; (B  (mp  ct)  Ln ) ; B  ; An ; norm1Ln
  log=. log lcat ('trsmrltu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  |:@trl1pick) t02m normitc))) L1n ; (B  (mp  ct)  L1n) ; B  ; An ; norm1L1n
  log=. log lcat ('trsmrljn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  + @trlpick ) t02m normitc))) Ln  ; (B  (mp  + )  Ln ) ; B  ; An ; norm1Ln
  log=. log lcat ('trsmrlju'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  + @trl1pick) t02m normitc))) L1n ; (B  (mp  + )  L1n) ; B  ; An ; norm1L1n
  log=. log lcat ('trsmrlcn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  ct@trlpick ) t02m normitc))) Ln  ; (B  (mp  |:)  Ln ) ; B  ; An ; norm1Ln
  log=. log lcat ('trsmrlcu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  ct@trl1pick) t02m normitc))) L1n ; (B  (mp  |:)  L1n) ; B  ; An ; norm1L1n
  log=. log lcat ('trsmrunn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp     trupick ) t02m normitc))) Un  ; (B   mp       Un ) ; B  ; An ; normiUn
  log=. log lcat ('trsmrunu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp     tru1pick) t02m normitc))) U1n ; (B   mp       U1n) ; B  ; An ; normiU1n
  log=. log lcat ('trsmrutn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  |:@trupick ) t02m normitc))) Un  ; (B  (mp  ct)  Un ) ; B  ; An ; norm1Un
  log=. log lcat ('trsmrutu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  |:@tru1pick) t02m normitc))) U1n ; (B  (mp  ct)  U1n) ; B  ; An ; norm1U1n
  log=. log lcat ('trsmrujn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  + @trupick ) t02m normitc))) Un  ; (B  (mp  + )  Un ) ; B  ; An ; norm1Un
  log=. log lcat ('trsmruju'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  + @tru1pick) t02m normitc))) U1n ; (B  (mp  + )  U1n) ; B  ; An ; norm1U1n
  log=. log lcat ('trsmrucn'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  ct@trupick ) t02m normitc))) Un  ; (B  (mp  |:)  Un ) ; B  ; An ; norm1Un
  log=. log lcat ('trsmrucu'         tdyad  ((3&{::)`(1&{::)`]`nan`nan`((mp  ct@tru1pick) t02m normitc))) U1n ; (B  (mp  |:)  U1n) ; B  ; An ; norm1U1n
)

NB. ---------------------------------------------------------
NB. testbasicsm
NB.
NB. Description:
NB.   Adv. to make verb to test basic matrix equation solvers
NB.
NB. Syntax:
NB.   log=. (mkmat testbasicsm) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testbasicsm_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testbasicsm_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testbasicsm_mt_ 150 200

testbasicsm=: 1 : 'testbasictrsm_mt_@(u@(2 # >./) ; u) [ require@''math/mt/external/blis/sm'' [ require@''math/mt/external/blas/sm'''

NB. ---------------------------------------------------------
NB. testbasic
NB.
NB. Description:
NB.   Adv. to make verb to test basic operations all levels
NB.
NB. Syntax:
NB.   log=. (mkmat testbasic) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testbasic_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testbasic_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testbasic_mt_ 150 200

testbasic=: 1 : '(u testbasicsm_mt_) ,&.>~ nolog_mt_`(u testbasicsv_mt_)@.(=/) ,&.>~ (u testbasicmm_mt_) ,&.>~ (u testbasicmv_mt_) ,&.>~ (u testbasicr2k_mt_) ,&.>~ (u testbasicrk_mt_) ,&.>~ nolog_mt_`(u testbasicr2_mt_)@.(=/) ,&.>~ (u testbasicr_mt_)'
