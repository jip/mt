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

reverse=: 2 : '[:`]`(m&(|.&.> upd_mt_))@.(*@(n&{::))'

NB. ---------------------------------------------------------
NB. ger
NB.
NB. Description:
NB.   Adv. to make monad to perform the rank 1 operation:
NB.     A := alpha * x * op(y) + A
NB.   where op(y) is either y^T or y^H
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
NB.   where op(y) is either y^T or y^H
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
NB.   where A is Hermitian (symmetric), and op(x) is either
NB.   x^T or x^H
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
NB.   where A is Hermitian (symmetric)
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
NB.   where A is Hermitian (symmetric), op1(x) is either x^T
NB.   or x^H, and op2(alpha) is either alpha or conj(alpha)
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
NB.   where A is Hermitian (symmetric)
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
NB.   where C is Hermitian (symmetric), and op(A) is either
NB.   A^T or A^H
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
NB.   where C is Hermitian (symmetric)
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
NB.   where C is Hermitian (symmetric), op1(M) is either M^T
NB.   or M^H and op2(alpha) is either alpha or conj(alpha)
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
NB.   where C is Hermitian (symmetric)
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
NB.   where A is Hermitian (symmetric)
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
NB.   where A is Hermitian (symmetric)
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
NB.   where A is triangular
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
NB.   where A is triangular
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
NB.   where opX(M) is either M, M^T, conj(M) or M^H
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
NB.   where A is Hermitian (symmetric), op1(A) is either A or
NB.   conj(A), and op2(B) is either B, B^T, conj(B) or B^H
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
NB.   where A is Hermitian (symmetric), op1(A) is either A or
NB.   conj(A), and op2(B) is either B, B^T, conj(B) or B^H
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
NB.   where A is triangular, and op(A) is either A, A^T,
NB.   conj(A) or A^H
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
NB.   where A is triangular, and op(A) is either A, A^T,
NB.   conj(A) or A^H
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
NB.   where A is triangular, and opX(M) is either M, M^T,
NB.   conj(M) or M^H
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
NB.   where A is triangular, and opX(M) is either M, M^T,
NB.   conj(M) or M^H
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
NB.   where A is triangular
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
NB.   where A is triangular
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
NB.   where A is triangular
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
NB.   where A is triangular, and op(A) is either A, A^T,
NB.   conj(A) or A^H
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
NB.   testbasicger x ; y ; A
NB. where
NB.   x - m-vector
NB.   y - n-vector
NB.   A - m×n-matrix

testbasicger=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  inc=. 1 2 _1 _2
  'x y A'=. y
  argsd=. { (<"0 dcoeff) ; (< x) ; (<"0 inc) ; (< y) ; (<"0 inc) ; < < A
  argsz=. { (<"0 zcoeff) ; (< x) ; (<"0 inc) ; (< y) ; (<"0 inc) ; < < A

  NB. for every i feed the tuple (alpha_i ; expanded_x_i ; incx_i ; expanded_y_i ; incy_i ; A) to tmonad
  ('dger_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(geru chk4r)))@(3 expand 4)@(1 expand 2)@>"0 argsd
  ('zgeru_mtbla_' tmonad (]`]`(_."_)`(_."_)`(geru chk4r)))@(3 expand 4)@(1 expand 2)@>"0 argsz
  ('zgerc_mtbla_' tmonad (]`]`(_."_)`(_."_)`(gerc chk4r)))@(3 expand 4)@(1 expand 2)@>"0 argsz

  EMPTY
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
NB.   testbasicher x ; AA
NB. where
NB.   x  - n-vector
NB.   AA - n×n-matrix, A material with real diagonal

testbasicher=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  inc=. 1 2 _1 _2
  'x AA'=. y
  args=. { (<"0 dcoeff) ; (< x) ; (<"0 inc) ; < < AA

  NB. for every i feed the tuple (alpha_i ; expanded_x_i ; incx_i ; AA) to tmonad
  ('dsyrl_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrl chk5r trlpick)))@(1 expand 2)@>"0 args
  ('dsyru_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syru chk5r trupick)))@(1 expand 2)@>"0 args
  ('zherl_mtbla_' tmonad (]`]`(_."_)`(_."_)`(herl chk5r trlpick)))@(1 expand 2)@>"0 args
  ('zheru_mtbla_' tmonad (]`]`(_."_)`(_."_)`(heru chk5r trupick)))@(1 expand 2)@>"0 args

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbasicr
NB.
NB. Description:
NB.   Adv. to make verb to test basic rank 1 operations
NB.
NB. Syntax:
NB.   vtest=. mkmat testbasicr
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms; is called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testbasicr_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testbasicr_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbasicr_mt_ 150 200

testbasicr=: 1 : 'EMPTY [ testbasicher_mt_@(u@{. ; (9&o. upddiag_mt_)@u)^:(=/) [ testbasicger_mt_@(u@{. ; u@{: ; u) [ load@''math/mt/test/blas/r'''

NB. ---------------------------------------------------------
NB. testbasicher2
NB.
NB. Description:
NB.   Test hermitian (symmetric) rank 2 operations:
NB.   - DSYR2 ZHER2 (BLAS)
NB.   by vectors and Hermitian (symmetric) matrix
NB.
NB. Syntax:
NB.   testbasicher2 x ; y ; AA
NB. where
NB.   x  - n-vector
NB.   y  - n-vector
NB.   AA - n×n-matrix, A material with real diagonal

testbasicher2=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  inc=. 1 2 _1 _2
  'x y AA'=. y
  argsd=. { (<"0 dcoeff) ; (< x) ; (<"0 inc) ; (< y) ; (<"0 inc) ; < < AA
  argsz=. { (<"0 zcoeff) ; (< x) ; (<"0 inc) ; (< y) ; (<"0 inc) ; < < AA

  NB. for every i feed the tuple (alpha_i ; expanded_x_i ; incx_i ; expanded_y_i ; incy_i ; AA) to tmonad
  ('dsyr2l_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2l chk6r2 trlpick)))@(3 expand 4)@(1 expand 2)@>"0 argsd
  ('dsyr2u_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2u chk6r2 trupick)))@(3 expand 4)@(1 expand 2)@>"0 argsd
  ('zher2l_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2l chk6r2 trlpick)))@(3 expand 4)@(1 expand 2)@>"0 argsz
  ('zher2u_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2u chk6r2 trupick)))@(3 expand 4)@(1 expand 2)@>"0 argsz

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbasicr2
NB.
NB. Description:
NB.   Adv. to make verb to test basic rank 2 operations
NB.
NB. Syntax:
NB.   vtest=. mkmat testbasicr2
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (n,n)
NB.   vtest - monad to test algorithms; is called as:
NB.             vtest (n,n)
NB.   (n,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testbasicr2_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testbasicr2_mt_ 200 200
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbasicr2_mt_ 250 250

testbasicr2=: 1 : 'EMPTY [ testbasicher2_mt_@(u@{. ; u@{: ; (9&o. upddiag_mt_)@u) [ load@''math/mt/test/blas/r2'''

NB. ---------------------------------------------------------
NB. testbasicsyrk
NB.
NB. Description:
NB.   Test symmetric rank k operations:
NB.   - xSYRK (BLAS)
NB.   by general and symmetric matrices
NB.
NB. Syntax:
NB.   testbasicsyrk A ; CC
NB. where
NB.   A  - m×n-matrix
NB.   CC - k×k-matrix, C material
NB.   k  = max(m,n)

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
  ('dsyrkln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkln chk4rk trlpick)))@>"0 argsdl
  ('dsyrklt_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrklt chk4rk trlpick)))@>"0 argsdg
  ('dsyrkun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkun chk4rk trupick)))@>"0 argsdl
  ('dsyrkut_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkut chk4rk trupick)))@>"0 argsdg
  ('zsyrkln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkln chk4rk trlpick)))@>"0 argszl
  ('zsyrklt_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrklt chk4rk trlpick)))@>"0 argszg
  ('zsyrkun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkun chk4rk trupick)))@>"0 argszl
  ('zsyrkut_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkut chk4rk trupick)))@>"0 argszg

  EMPTY
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
NB.   testbasicherk A ; CC
NB. where
NB.   A  - m×n-matrix
NB.   CC - k×k-matrix, C material with real diagonal
NB.   k  = max(m,n)

testbasicherk=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  'A CC'=. y
  mn=. <./ 'm n'=. $ A
  CCmn=. (2 # mn) {. CC
  argsl=. { (<"0 dcoeff) ; (< A) ; (<"0 dcoeff) ; < < CCmn [^:(m < n) CC
  argsg=. { (<"0 dcoeff) ; (< A) ; (<"0 dcoeff) ; < < CCmn [^:(m > n) CC

  NB. for every i feed the tuple (alpha_i ; A ; beta_i ; CC) to tmonad
  ('zherkln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(herkln chk4rk trlpick)))@>"0 argsl
  ('zherklc_mtbla_' tmonad (]`]`(_."_)`(_."_)`(herklc chk4rk trlpick)))@>"0 argsg
  ('zherkun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(herkun chk4rk trupick)))@>"0 argsl
  ('zherkuc_mtbla_' tmonad (]`]`(_."_)`(_."_)`(herkuc chk4rk trupick)))@>"0 argsg

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbasicrk
NB.
NB. Description:
NB.   Adv. to make verb to test basic rank k operations
NB.
NB. Syntax:
NB.   vtest=. mkmat testbasicrk
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms; is called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testbasicrk_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testbasicrk_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbasicrk_mt_ 150 200

testbasicrk=: 1 : 'EMPTY [ testbasicherk_mt_@(u ; (9&o. upddiag_mt_)@u@(2 # >./)) [ testbasicsyrk_mt_@(u ; u@(2 # >./)) [ load@''math/mt/test/blas/rk'''

NB. ---------------------------------------------------------
NB. testbasicsyr2k
NB.
NB. Description:
NB.   Test symmetric rank 2k operations:
NB.   - xSYR2K (BLAS)
NB.   by general and symmetric matrices
NB.
NB. Syntax:
NB.   testbasicsyr2k A ; B ; CC
NB. where
NB.   A  - m×n-matrix
NB.   B  - m×n-matrix
NB.   CC - k×k-matrix, C material
NB.   k  = max(m,n)

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
  ('dsyr2kln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kln chk5r2k trlpick)))@>"0 argsdl
  ('dsyr2klt_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2klt chk5r2k trlpick)))@>"0 argsdg
  ('dsyr2kun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kun chk5r2k trupick)))@>"0 argsdl
  ('dsyr2kut_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kut chk5r2k trupick)))@>"0 argsdg
  ('zsyr2kln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kln chk5r2k trlpick)))@>"0 argszl
  ('zsyr2klt_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2klt chk5r2k trlpick)))@>"0 argszg
  ('zsyr2kun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kun chk5r2k trupick)))@>"0 argszl
  ('zsyr2kut_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kut chk5r2k trupick)))@>"0 argszg

  EMPTY
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
NB.   testbasicher2k A ; B ; CC
NB. where
NB.   A  - m×n-matrix
NB.   B  - m×n-matrix
NB.   CC - k×k-matrix, C material with real diagonal
NB.   k  = max(m,n)

testbasicher2k=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  'A B CC'=. y
  mn=. <./ 'm n'=. $ A
  CCmn=. (2 # mn) {. CC
  argsl=. { (<"0 zcoeff) ; (< A) ; (< B) ; (<"0 dcoeff) ; < < CCmn [^:(m < n) CC
  argsg=. { (<"0 zcoeff) ; (< A) ; (< B) ; (<"0 dcoeff) ; < < CCmn [^:(m > n) CC

  NB. for every i feed the tuple (alpha_i ; A ; B ; beta_i ; CC) to tmonad
  ('zher2kln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2kln chk5r2k trlpick)))@>"0 argsl
  ('zher2klc_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2klc chk5r2k trlpick)))@>"0 argsg
  ('zher2kun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2kun chk5r2k trupick)))@>"0 argsl
  ('zher2kuc_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2kuc chk5r2k trupick)))@>"0 argsg

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbasicr2k
NB.
NB. Description:
NB.   Adv. to make verb to test basic rank 2k operations
NB.
NB. Syntax:
NB.   vtest=. mkmat testbasicr2k
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms; is called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testbasicrk_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testbasicrk_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbasicrk_mt_ 150 200

testbasicr2k=: 1 : 'EMPTY [ testbasicher2k_mt_@(u ; u ; (9&o. upddiag_mt_)@u@(2 # >./)) [ testbasicsyr2k_mt_@(u ; u ; u@(2 # >./)) [ load@''math/mt/test/blas/r2k'''

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
NB.   testbasicgemv A ; x ; y
NB. where
NB.   A - m×n-matrix
NB.   x - k-vector
NB.   y - k-vector
NB.   k = max(m,n)

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
  ('(+/ .*)'        tdyad  ((0&{::)`(1&{::)`0:`(_."_)`(_."_)`0:            ))                                                   A  ;    xn
  ('mp'             tdyad  ((0&{::)`(1&{::)`0:`(_."_)`(_."_)`0:            ))                                                   A  ;    xn

  NB. for every i feed the tuple (alpha_i ; A ; expanded_x_i ; incx_i ; beta_i ; expanded_y_i ; incy_i) to tmonad
  ('dgemvn_mtbla_' tmonad (         ]      `] `(_."_)`(_."_)`(gemvn chk1mv)))@(5 expand 6)@(2 expand 3)@>"0 { (<"0 dcoeff) ; (< A) ; (< xn) ; (<"0 inc) ; (<"0 dcoeff) ; (< ym) ; < <"0 inc
  ('dgemvt_mtbla_' tmonad (         ]      `] `(_."_)`(_."_)`(gemvt chk1mv)))@(5 expand 6)@(2 expand 3)@>"0 { (<"0 dcoeff) ; (< A) ; (< xm) ; (<"0 inc) ; (<"0 dcoeff) ; (< yn) ; < <"0 inc
  ('zgemvn_mtbla_' tmonad (         ]      `] `(_."_)`(_."_)`(gemvn chk1mv)))@(5 expand 6)@(2 expand 3)@>"0 { (<"0 zcoeff) ; (< A) ; (< xn) ; (<"0 inc) ; (<"0 zcoeff) ; (< ym) ; < <"0 inc
  ('zgemvt_mtbla_' tmonad (         ]      `] `(_."_)`(_."_)`(gemvt chk1mv)))@(5 expand 6)@(2 expand 3)@>"0 { (<"0 zcoeff) ; (< A) ; (< xm) ; (<"0 inc) ; (<"0 zcoeff) ; (< ym) ; < <"0 inc
  ('zgemvc_mtbla_' tmonad (         ]      `] `(_."_)`(_."_)`(gemvc chk1mv)))@(5 expand 6)@(2 expand 3)@>"0 { (<"0 zcoeff) ; (< A) ; (< xm) ; (<"0 inc) ; (<"0 zcoeff) ; (< ym) ; < <"0 inc

  EMPTY
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
NB.   testbasichemv AA ; x ; y
NB. where
NB.   AA - n×n-matrix, A material with real diagonal
NB.   x  - n-vector
NB.   y  - n-vector

testbasichemv=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  inc=. 1 2 _1 _2
  'AA x y'=. y
  argsd=. { (<"0 dcoeff) ; (< AA) ; (< x) ; (<"0 inc) ; (<"0 dcoeff) ; (< y) ; < <"0 inc
  argsz=. { (<"0 zcoeff) ; (< AA) ; (< x) ; (<"0 inc) ; (<"0 zcoeff) ; (< y) ; < <"0 inc

  NB. for every i feed the tuple (alpha_i ; AA ; expanded_x_i ; incx_i ; beta_i ; expanded_y_i ; incy_i) to tmonad
  ('dsymvl_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symvl chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 argsd
  ('dsymvu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symvu chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 argsd
  ('zhemvl_mtbla_' tmonad (]`]`(_."_)`(_."_)`(hemvl chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 argsz
  ('zhemvu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(hemvu chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 argsz

  EMPTY
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
NB.   testbasictrmv AA ; x
NB. where
NB.   AA - n×n-matrix, A material
NB.   x  - n-vector

testbasictrmv=: 3 : 0
  inc=. 1 2 _1 _2
  'AA y'=. y
  args=. { (< AA) ; (< y) ; < <"0 inc

  NB. for every i feed the tuple (AA ; expanded_x_i ; incx_i) to tmonad
  ('dtrmvlnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnn chk3mv)))@(1 expand 2)@>"0 args
  ('dtrmvlnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnu chk3mv)))@(1 expand 2)@>"0 args
  ('dtrmvltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltn chk3mv)))@(1 expand 2)@>"0 args
  ('dtrmvltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltu chk3mv)))@(1 expand 2)@>"0 args
  ('dtrmvunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunn chk3mv)))@(1 expand 2)@>"0 args
  ('dtrmvunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunu chk3mv)))@(1 expand 2)@>"0 args
  ('dtrmvutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutn chk3mv)))@(1 expand 2)@>"0 args
  ('dtrmvutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutu chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvlnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnn chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvlnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnu chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltn chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltu chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvlcn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlcn chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvlcu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlcu chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunn chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunu chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutn chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutu chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvucn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvucn chk3mv)))@(1 expand 2)@>"0 args
  ('ztrmvucu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvucu chk3mv)))@(1 expand 2)@>"0 args

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbasicmv
NB.
NB. Description:
NB.   Adv. to make verb to test basic matrix-vector operations
NB.
NB. Syntax:
NB.   vtest=. mkmat testbasicmv
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms; is called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testbasicmv_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testbasicmv_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbasicmv_mt_ 150 200

testbasicmv=: 1 : 'EMPTY [ (testbasictrmv_mt_@(u ; u@{.) [ testbasichemv_mt_@((9&o. upddiag_mt_)@u ; (u ; u)@{.))^:(=/) [ testbasicgemv_mt_@(u ; (u ; u)@(>./)) [ load@''math/mt/test/blas/mv'''

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
NB.   testbasicgemm As ; Bs ; C
NB. where
NB.   As - m×(m+n)-matrix, A material
NB.   Bs - (m+n)×n-matrix, B material
NB.   C  - m×n-matrix
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
  ('(+/ .*)'         tdyad  ((0&{::)`(1&{::)`0:`(_."_)`(_."_)`0:                           ))@(c (0 shrink 1)  {.   )@>"0 {                        As  ;  < <  Bs
  ('mp'              tdyad  ((0&{::)`(1&{::)`0:`(_."_)`(_."_)`0:                           ))@(c (0 shrink 1)  {.   )@>"0 {                        As  ;  < <  Bs

  NB. for every i feed the tuple (alpha_i ; A_i ; B_i ; beta_i ; C) to tmonad

  ('dgemmnn_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsdnn
  ('dgemmnt_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsdnt
  ('dgemmtn_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsdtn
  ('dgemmtt_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtt  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsdtt

  ('zgemmnn_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 { (<"0 zcoeff) ;         As  ; (<    Bs) ; (<"0 zcoeff) ; < < C
  ('zgemmnt_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ;         As  ; (< |: Bs) ; (<"0 zcoeff) ; < < C
  ('zgemmnc_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnc  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ;         As  ; (< ct Bs) ; (<"0 zcoeff) ; < < C
  ('zgemmtn_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 { (<"0 zcoeff) ; (|: L:0 As) ; (<    Bs) ; (<"0 zcoeff) ; < < C
  ('zgemmtt_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtt  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (|: L:0 As) ; (< |: Bs) ; (<"0 zcoeff) ; < < C
  ('zgemmtc_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (|: L:0 As) ; (< ct Bs) ; (<"0 zcoeff) ; < < C
  ('zgemmcn_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmcn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 { (<"0 zcoeff) ; (ct L:0 As) ; (<    Bs) ; (<"0 zcoeff) ; < < C
  ('zgemmct_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmct  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (ct L:0 As) ; (< |: Bs) ; (<"0 zcoeff) ; < < C
  ('zgemmcc_mtbla_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmcc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (ct L:0 As) ; (< ct Bs) ; (<"0 zcoeff) ; < < C

  ('gemmnn_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  ('gemmnt_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  ('gemmnj_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnj  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  ('gemmnc_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnc  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  ('gemmtn_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  ('gemmtt_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtt  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  ('gemmtj_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtj  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  ('gemmtc_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  ('gemmjn_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmjn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  ('gemmjt_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmjt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  ('gemmjj_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmjj  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  ('gemmjc_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmjc  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  ('gemmcn_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmcn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  ('gemmct_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmct  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  ('gemmcj_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmcj  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  ('gemmcc_mtbli_'   tmonad (        ]      `] `(_."_)`(_."_)`(             gemmcc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc

  ('dgemmnn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsdnn
  ('dgemmnt_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsdnt
  ('dgemmtn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsdtn
  ('dgemmtt_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtt  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsdtt

  ('zgemmnn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  ('zgemmnt_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  ('zgemmnj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnj  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  ('zgemmnc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmnc  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  ('zgemmtn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  ('zgemmtt_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtt  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  ('zgemmtj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtj  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  ('zgemmtc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmtc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  ('zgemmjn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmjn  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  ('zgemmjt_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmjt  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  ('zgemmjj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmjj  chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  ('zgemmjc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmjc  chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  ('zgemmcn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmcn  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  ('zgemmct_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmct  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  ('zgemmcj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmcj  chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  ('zgemmcc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`(             gemmcc  chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc

  EMPTY
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
NB.   As - m×(m+n)-matrix, A material
NB.   Bs - (m+n)×m-matrix, B material
NB.   C  - m×m-matrix
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
  EMPTY return.
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

  ('gemmlnn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  ('gemmlnt_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  ('gemmlnj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmnj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  ('gemmlnc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmnc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  ('gemmltn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  ('gemmltt_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  ('gemmltj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmtj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  ('gemmltc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmtc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  ('gemmljn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmjn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  ('gemmljt_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmjt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  ('gemmljj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmjj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  ('gemmljc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmjc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  ('gemmlcn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmcn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  ('gemmlct_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmct) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  ('gemmlcj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmcj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  ('gemmlcc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmcc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc
  ('gemmunn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  ('gemmunt_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  ('gemmunj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmnj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  ('gemmunc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmnc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  ('gemmutn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  ('gemmutt_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  ('gemmutj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmtj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  ('gemmutc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmtc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  ('gemmujn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmjn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  ('gemmujt_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmjt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  ('gemmujj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmjj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  ('gemmujc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmjc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  ('gemmucn_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmcn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  ('gemmuct_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmct) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  ('gemmucj_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmcj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  ('gemmucc_mtbli_'  tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmcc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc

  ('dgemmlnn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsdnn
  ('dgemmlnt_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsdnt
  ('dgemmltn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsdtn
  ('dgemmltt_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsdtt
  ('dgemmunn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsdnn
  ('dgemmunt_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsdnt
  ('dgemmutn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsdtn
  ('dgemmutt_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsdtt

  ('zgemmlnn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  ('zgemmlnt_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  ('zgemmlnj_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmnj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  ('zgemmlnc_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmnc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  ('zgemmltn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  ('zgemmltt_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  ('zgemmltj_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmtj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  ('zgemmltc_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmtc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  ('zgemmljn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmjn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  ('zgemmljt_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmjt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  ('zgemmljj_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmjj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  ('zgemmljc_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmjc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  ('zgemmlcn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmcn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  ('zgemmlct_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmct) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  ('zgemmlcj_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmcj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  ('zgemmlcc_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: suxly gemmcc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc
  ('zgemmunn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmnn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsann
  ('zgemmunt_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmnt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsant
  ('zgemmunj_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmnj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsanj
  ('zgemmunc_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmnc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsanc
  ('zgemmutn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmtn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatn
  ('zgemmutt_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmtt) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatt
  ('zgemmutj_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmtj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsatj
  ('zgemmutc_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmtc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsatc
  ('zgemmujn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmjn) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajn
  ('zgemmujt_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmjt) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajt
  ('zgemmujj_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmjj) chk1mm)))@(c (1 shrink 2)  {.   )@>"0 argsajj
  ('zgemmujc_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmjc) chk1mm)))@(c (1 shrink 2) ({."1))@>"0 argsajc
  ('zgemmucn_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmcn) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacn
  ('zgemmuct_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmct) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsact
  ('zgemmucj_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmcj) chk1mm)))@(# (1 shrink 2)  {.   )@>"0 argsacj
  ('zgemmucc_mtbli_' tmonad (        ]      `] `(_."_)`(_."_)`((4&{:: slxuy gemmcc) chk1mm)))@(# (1 shrink 2) ({."1))@>"0 argsacc

  EMPTY
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
NB.   testbasicsymm AA ; B ; C
NB. where
NB.   AA - k×k-matrix, A material
NB.   B  - m×n-matrix
NB.   C  - m×n-matrix
NB.   k  = max(m,n)
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

  ('dsymmll_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(symmllnn chk2mm trlpick)))@>"0 argsdm
  ('dsymmlu_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(symmlunn chk2mm trupick)))@>"0 argsdm
  ('dsymmrl_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(symmrlnn chk2mm trlpick)))@>"0 argsdn
  ('dsymmru_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(symmrunn chk2mm trupick)))@>"0 argsdn

  ('zsymmll_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(symmllnn chk2mm trlpick)))@>"0 argszm
  ('zsymmlu_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(symmlunn chk2mm trupick)))@>"0 argszm
  ('zsymmrl_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(symmrlnn chk2mm trlpick)))@>"0 argszn
  ('zsymmru_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(symmrunn chk2mm trupick)))@>"0 argszn

  ('symmllnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmllnn chk2mm trlpick)))@>"0 argsam
  ('symmllnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmllnt chk2mm trlpick)))@>"0 argsam
  ('symmllnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmllnj chk2mm trlpick)))@>"0 argsam
  ('symmllnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmllnc chk2mm trlpick)))@>"0 argsam
  ('symmlljn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlljn chk2mm trlpick)))@>"0 argsam
  ('symmlljt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlljt chk2mm trlpick)))@>"0 argsam
  ('symmlljj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlljj chk2mm trlpick)))@>"0 argsam
  ('symmlljc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlljc chk2mm trlpick)))@>"0 argsam
  ('symmlunn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlunn chk2mm trupick)))@>"0 argsam
  ('symmlunt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlunt chk2mm trupick)))@>"0 argsam
  ('symmlunj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlunj chk2mm trupick)))@>"0 argsam
  ('symmlunc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlunc chk2mm trupick)))@>"0 argsam
  ('symmlujn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlujn chk2mm trupick)))@>"0 argsam
  ('symmlujt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlujt chk2mm trupick)))@>"0 argsam
  ('symmlujj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlujj chk2mm trupick)))@>"0 argsam
  ('symmlujc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmlujc chk2mm trupick)))@>"0 argsam
  ('symmrlnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrlnn chk2mm trlpick)))@>"0 argsan
  ('symmrlnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrlnt chk2mm trlpick)))@>"0 argsan
  ('symmrlnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrlnj chk2mm trlpick)))@>"0 argsan
  ('symmrlnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrlnc chk2mm trlpick)))@>"0 argsan
  ('symmrljn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrljn chk2mm trlpick)))@>"0 argsan
  ('symmrljt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrljt chk2mm trlpick)))@>"0 argsan
  ('symmrljj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrljj chk2mm trlpick)))@>"0 argsan
  ('symmrljc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrljc chk2mm trlpick)))@>"0 argsan
  ('symmrunn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrunn chk2mm trupick)))@>"0 argsan
  ('symmrunt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrunt chk2mm trupick)))@>"0 argsan
  ('symmrunj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrunj chk2mm trupick)))@>"0 argsan
  ('symmrunc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrunc chk2mm trupick)))@>"0 argsan
  ('symmrujn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrujn chk2mm trupick)))@>"0 argsan
  ('symmrujt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrujt chk2mm trupick)))@>"0 argsan
  ('symmrujj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrujj chk2mm trupick)))@>"0 argsan
  ('symmrujc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(symmrujc chk2mm trupick)))@>"0 argsan

  ('dsymmllnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmllnn chk2mm trlpick)))@>"0 argsdm
  ('dsymmllnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmllnt chk2mm trlpick)))@>"0 argsdm
  ('dsymmlunn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlunn chk2mm trupick)))@>"0 argsdm
  ('dsymmlunt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlunt chk2mm trupick)))@>"0 argsdm
  ('dsymmrlnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrlnn chk2mm trlpick)))@>"0 argsdn
  ('dsymmrlnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrlnt chk2mm trlpick)))@>"0 argsdn
  ('dsymmrunn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrunn chk2mm trupick)))@>"0 argsdn
  ('dsymmrunt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrunt chk2mm trupick)))@>"0 argsdn

  ('zsymmllnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmllnn chk2mm trlpick)))@>"0 argsam
  ('zsymmllnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmllnt chk2mm trlpick)))@>"0 argsam
  ('zsymmllnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmllnj chk2mm trlpick)))@>"0 argsam
  ('zsymmllnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmllnc chk2mm trlpick)))@>"0 argsam
  ('zsymmlljn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlljn chk2mm trlpick)))@>"0 argsam
  ('zsymmlljt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlljt chk2mm trlpick)))@>"0 argsam
  ('zsymmlljj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlljj chk2mm trlpick)))@>"0 argsam
  ('zsymmlljc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlljc chk2mm trlpick)))@>"0 argsam
  ('zsymmlunn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlunn chk2mm trupick)))@>"0 argsam
  ('zsymmlunt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlunt chk2mm trupick)))@>"0 argsam
  ('zsymmlunj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlunj chk2mm trupick)))@>"0 argsam
  ('zsymmlunc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlunc chk2mm trupick)))@>"0 argsam
  ('zsymmlujn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlujn chk2mm trupick)))@>"0 argsam
  ('zsymmlujt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlujt chk2mm trupick)))@>"0 argsam
  ('zsymmlujj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlujj chk2mm trupick)))@>"0 argsam
  ('zsymmlujc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmlujc chk2mm trupick)))@>"0 argsam
  ('zsymmrlnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrlnn chk2mm trlpick)))@>"0 argsan
  ('zsymmrlnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrlnt chk2mm trlpick)))@>"0 argsan
  ('zsymmrlnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrlnj chk2mm trlpick)))@>"0 argsan
  ('zsymmrlnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrlnc chk2mm trlpick)))@>"0 argsan
  ('zsymmrljn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrljn chk2mm trlpick)))@>"0 argsan
  ('zsymmrljt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrljt chk2mm trlpick)))@>"0 argsan
  ('zsymmrljj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrljj chk2mm trlpick)))@>"0 argsan
  ('zsymmrljc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrljc chk2mm trlpick)))@>"0 argsan
  ('zsymmrunn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrunn chk2mm trupick)))@>"0 argsan
  ('zsymmrunt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrunt chk2mm trupick)))@>"0 argsan
  ('zsymmrunj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrunj chk2mm trupick)))@>"0 argsan
  ('zsymmrunc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrunc chk2mm trupick)))@>"0 argsan
  ('zsymmrujn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrujn chk2mm trupick)))@>"0 argsan
  ('zsymmrujt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrujt chk2mm trupick)))@>"0 argsan
  ('zsymmrujj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrujj chk2mm trupick)))@>"0 argsan
  ('zsymmrujc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(symmrujc chk2mm trupick)))@>"0 argsan

  EMPTY
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
NB.   testbasichemm AA ; B ; C
NB. where
NB.   AA - k×k-matrix, A material with real diagonal
NB.   B  - m×n-matrix
NB.   C  - m×n-matrix
NB.   k  = max(m,n)
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

  ('zhemmll_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(hemmllnn chk2mm trlpick)))@>"0 argszm
  ('zhemmlu_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(hemmlunn chk2mm trupick)))@>"0 argszm
  ('zhemmrl_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(hemmrlnn chk2mm trlpick)))@>"0 argszn
  ('zhemmru_mtbla_'   tmonad (]`]`(_."_)`(_."_)`(hemmrunn chk2mm trupick)))@>"0 argszn

  ('hemmllnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmllnn chk2mm trlpick)))@>"0 argsam
  ('hemmllnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmllnt chk2mm trlpick)))@>"0 argsam
  ('hemmllnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmllnj chk2mm trlpick)))@>"0 argsam
  ('hemmllnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmllnc chk2mm trlpick)))@>"0 argsam
  ('hemmlljn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlljn chk2mm trlpick)))@>"0 argsam
  ('hemmlljt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlljt chk2mm trlpick)))@>"0 argsam
  ('hemmlljj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlljj chk2mm trlpick)))@>"0 argsam
  ('hemmlljc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlljc chk2mm trlpick)))@>"0 argsam
  ('hemmlunn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlunn chk2mm trupick)))@>"0 argsam
  ('hemmlunt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlunt chk2mm trupick)))@>"0 argsam
  ('hemmlunj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlunj chk2mm trupick)))@>"0 argsam
  ('hemmlunc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlunc chk2mm trupick)))@>"0 argsam
  ('hemmlujn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlujn chk2mm trupick)))@>"0 argsam
  ('hemmlujt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlujt chk2mm trupick)))@>"0 argsam
  ('hemmlujj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlujj chk2mm trupick)))@>"0 argsam
  ('hemmlujc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmlujc chk2mm trupick)))@>"0 argsam
  ('hemmrlnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrlnn chk2mm trlpick)))@>"0 argsan
  ('hemmrlnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrlnt chk2mm trlpick)))@>"0 argsan
  ('hemmrlnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrlnj chk2mm trlpick)))@>"0 argsan
  ('hemmrlnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrlnc chk2mm trlpick)))@>"0 argsan
  ('hemmrljn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrljn chk2mm trlpick)))@>"0 argsan
  ('hemmrljt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrljt chk2mm trlpick)))@>"0 argsan
  ('hemmrljj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrljj chk2mm trlpick)))@>"0 argsan
  ('hemmrljc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrljc chk2mm trlpick)))@>"0 argsan
  ('hemmrunn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrunn chk2mm trupick)))@>"0 argsan
  ('hemmrunt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrunt chk2mm trupick)))@>"0 argsan
  ('hemmrunj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrunj chk2mm trupick)))@>"0 argsan
  ('hemmrunc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrunc chk2mm trupick)))@>"0 argsan
  ('hemmrujn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrujn chk2mm trupick)))@>"0 argsan
  ('hemmrujt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrujt chk2mm trupick)))@>"0 argsan
  ('hemmrujj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrujj chk2mm trupick)))@>"0 argsan
  ('hemmrujc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(hemmrujc chk2mm trupick)))@>"0 argsan

  ('zhemmllnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmllnn chk2mm trlpick)))@>"0 argsam
  ('zhemmllnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmllnt chk2mm trlpick)))@>"0 argsam
  ('zhemmllnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmllnj chk2mm trlpick)))@>"0 argsam
  ('zhemmllnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmllnc chk2mm trlpick)))@>"0 argsam
  ('zhemmlljn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlljn chk2mm trlpick)))@>"0 argsam
  ('zhemmlljt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlljt chk2mm trlpick)))@>"0 argsam
  ('zhemmlljj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlljj chk2mm trlpick)))@>"0 argsam
  ('zhemmlljc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlljc chk2mm trlpick)))@>"0 argsam
  ('zhemmlunn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlunn chk2mm trupick)))@>"0 argsam
  ('zhemmlunt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlunt chk2mm trupick)))@>"0 argsam
  ('zhemmlunj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlunj chk2mm trupick)))@>"0 argsam
  ('zhemmlunc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlunc chk2mm trupick)))@>"0 argsam
  ('zhemmlujn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlujn chk2mm trupick)))@>"0 argsam
  ('zhemmlujt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlujt chk2mm trupick)))@>"0 argsam
  ('zhemmlujj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlujj chk2mm trupick)))@>"0 argsam
  ('zhemmlujc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmlujc chk2mm trupick)))@>"0 argsam
  ('zhemmrlnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrlnn chk2mm trlpick)))@>"0 argsan
  ('zhemmrlnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrlnt chk2mm trlpick)))@>"0 argsan
  ('zhemmrlnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrlnj chk2mm trlpick)))@>"0 argsan
  ('zhemmrlnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrlnc chk2mm trlpick)))@>"0 argsan
  ('zhemmrljn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrljn chk2mm trlpick)))@>"0 argsan
  ('zhemmrljt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrljt chk2mm trlpick)))@>"0 argsan
  ('zhemmrljj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrljj chk2mm trlpick)))@>"0 argsan
  ('zhemmrljc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrljc chk2mm trlpick)))@>"0 argsan
  ('zhemmrunn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrunn chk2mm trupick)))@>"0 argsan
  ('zhemmrunt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrunt chk2mm trupick)))@>"0 argsan
  ('zhemmrunj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrunj chk2mm trupick)))@>"0 argsan
  ('zhemmrunc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrunc chk2mm trupick)))@>"0 argsan
  ('zhemmrujn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrujn chk2mm trupick)))@>"0 argsan
  ('zhemmrujt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrujt chk2mm trupick)))@>"0 argsan
  ('zhemmrujj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrujj chk2mm trupick)))@>"0 argsan
  ('zhemmrujc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(hemmrujc chk2mm trupick)))@>"0 argsan

  EMPTY
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
NB.   testbasictrmm AA ; B ; С
NB. where
NB.   AA - k×k-matrix, A material
NB.   B  - m×n-matrix
NB.   C  - m×n-matrix
NB.   k  = max(m,n)
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

  ('dtrmmllnn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnn  chk3mm trlpick )))@>"0 argsdm
  ('dtrmmllnu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnu  chk3mm trlpick )))@>"0 argsdm
  ('dtrmmlltn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltn  chk3mm trlpick )))@>"0 argsdm
  ('dtrmmlltu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltu  chk3mm trlpick )))@>"0 argsdm
  ('dtrmmlunn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunn  chk3mm trupick )))@>"0 argsdm
  ('dtrmmlunu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunu  chk3mm trupick )))@>"0 argsdm
  ('dtrmmlutn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutn  chk3mm trupick )))@>"0 argsdm
  ('dtrmmlutu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutu  chk3mm trupick )))@>"0 argsdm
  ('dtrmmrlnn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnn  chk3mm trlpick )))@>"0 argsdn
  ('dtrmmrlnu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnu  chk3mm trlpick )))@>"0 argsdn
  ('dtrmmrltn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltn  chk3mm trlpick )))@>"0 argsdn
  ('dtrmmrltu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltu  chk3mm trlpick )))@>"0 argsdn
  ('dtrmmrunn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunn  chk3mm trupick )))@>"0 argsdn
  ('dtrmmrunu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunu  chk3mm trupick )))@>"0 argsdn
  ('dtrmmrutn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutn  chk3mm trupick )))@>"0 argsdn
  ('dtrmmrutu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutu  chk3mm trupick )))@>"0 argsdn

  ('ztrmmllnn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnn  chk3mm trlpick )))@>"0 argszm
  ('ztrmmllnu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnu  chk3mm trlpick )))@>"0 argszm
  ('ztrmmlltn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltn  chk3mm trlpick )))@>"0 argszm
  ('ztrmmlltu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltu  chk3mm trlpick )))@>"0 argszm
  ('ztrmmllcn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcn  chk3mm trlpick )))@>"0 argszm
  ('ztrmmllcu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcu  chk3mm trlpick )))@>"0 argszm
  ('ztrmmlunn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunn  chk3mm trupick )))@>"0 argszm
  ('ztrmmlunu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunu  chk3mm trupick )))@>"0 argszm
  ('ztrmmlutn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutn  chk3mm trupick )))@>"0 argszm
  ('ztrmmlutu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutu  chk3mm trupick )))@>"0 argszm
  ('ztrmmlucn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucn  chk3mm trupick )))@>"0 argszm
  ('ztrmmlucu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucu  chk3mm trupick )))@>"0 argszm
  ('ztrmmrlnn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnn  chk3mm trlpick )))@>"0 argszn
  ('ztrmmrlnu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnu  chk3mm trlpick )))@>"0 argszn
  ('ztrmmrltn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltn  chk3mm trlpick )))@>"0 argszn
  ('ztrmmrltu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltu  chk3mm trlpick )))@>"0 argszn
  ('ztrmmrlcn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcn  chk3mm trlpick )))@>"0 argszn
  ('ztrmmrlcu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcu  chk3mm trlpick )))@>"0 argszn
  ('ztrmmrunn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunn  chk3mm trupick )))@>"0 argszn
  ('ztrmmrunu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunu  chk3mm trupick )))@>"0 argszn
  ('ztrmmrutn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutn  chk3mm trupick )))@>"0 argszn
  ('ztrmmrutu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutu  chk3mm trupick )))@>"0 argszn
  ('ztrmmrucn_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucn  chk3mm trupick )))@>"0 argszn
  ('ztrmmrucu_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucu  chk3mm trupick )))@>"0 argszn

  ('trmmllnn_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmllnn  chk3mm trlpick )))@>"0 argsam
  ('trmmllnu_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmllnu  chk3mm trlpick )))@>"0 argsam
  ('trmmlltn_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmlltn  chk3mm trlpick )))@>"0 argsam
  ('trmmlltu_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmlltu  chk3mm trlpick )))@>"0 argsam
  ('trmmlunn_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmlunn  chk3mm trupick )))@>"0 argsam
  ('trmmlunu_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmlunu  chk3mm trupick )))@>"0 argsam
  ('trmmlutn_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmlutn  chk3mm trupick )))@>"0 argsam
  ('trmmlutu_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmlutu  chk3mm trupick )))@>"0 argsam
  ('trmmrlnn_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmrlnn  chk3mm trlpick )))@>"0 argsan
  ('trmmrlnu_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmrlnu  chk3mm trlpick )))@>"0 argsan
  ('trmmrltn_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmrltn  chk3mm trlpick )))@>"0 argsan
  ('trmmrltu_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmrltu  chk3mm trlpick )))@>"0 argsan
  ('trmmrunn_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmrunn  chk3mm trupick )))@>"0 argsan
  ('trmmrunu_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmrunu  chk3mm trupick )))@>"0 argsan
  ('trmmrutn_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmrutn  chk3mm trupick )))@>"0 argsan
  ('trmmrutu_mtbli_'   tmonad (]`]`(_."_)`(_."_)`(trmmrutu  chk3mm trupick )))@>"0 argsan

  ('dtrmmllnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnn  chk3mm trlpick )))@>"0 argsdm
  ('dtrmmllnu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnu  chk3mm trlpick )))@>"0 argsdm
  ('dtrmmlltn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltn  chk3mm trlpick )))@>"0 argsdm
  ('dtrmmlltu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltu  chk3mm trlpick )))@>"0 argsdm
  ('dtrmmlunn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunn  chk3mm trupick )))@>"0 argsdm
  ('dtrmmlunu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunu  chk3mm trupick )))@>"0 argsdm
  ('dtrmmlutn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutn  chk3mm trupick )))@>"0 argsdm
  ('dtrmmlutu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutu  chk3mm trupick )))@>"0 argsdm
  ('dtrmmrlnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnn  chk3mm trlpick )))@>"0 argsdn
  ('dtrmmrlnu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnu  chk3mm trlpick )))@>"0 argsdn
  ('dtrmmrltn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltn  chk3mm trlpick )))@>"0 argsdn
  ('dtrmmrltu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltu  chk3mm trlpick )))@>"0 argsdn
  ('dtrmmrunn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunn  chk3mm trupick )))@>"0 argsdn
  ('dtrmmrunu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunu  chk3mm trupick )))@>"0 argsdn
  ('dtrmmrutn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutn  chk3mm trupick )))@>"0 argsdn
  ('dtrmmrutu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutu  chk3mm trupick )))@>"0 argsdn

  ('ztrmmllnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnn  chk3mm trlpick )))@>"0 argsam
  ('ztrmmllnu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnu  chk3mm trlpick )))@>"0 argsam
  ('ztrmmlltn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltn  chk3mm trlpick )))@>"0 argsam
  ('ztrmmlltu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltu  chk3mm trlpick )))@>"0 argsam
  ('ztrmmllcn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcn  chk3mm trlpick )))@>"0 argsam
  ('ztrmmllcu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcu  chk3mm trlpick )))@>"0 argsam
  ('ztrmmlunn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunn  chk3mm trupick )))@>"0 argsam
  ('ztrmmlunu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunu  chk3mm trupick )))@>"0 argsam
  ('ztrmmlutn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutn  chk3mm trupick )))@>"0 argsam
  ('ztrmmlutu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutu  chk3mm trupick )))@>"0 argsam
  ('ztrmmlucn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucn  chk3mm trupick )))@>"0 argsam
  ('ztrmmlucu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucu  chk3mm trupick )))@>"0 argsam
  ('ztrmmrlnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnn  chk3mm trlpick )))@>"0 argsan
  ('ztrmmrlnu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnu  chk3mm trlpick )))@>"0 argsan
  ('ztrmmrltn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltn  chk3mm trlpick )))@>"0 argsan
  ('ztrmmrltu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltu  chk3mm trlpick )))@>"0 argsan
  ('ztrmmrlcn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcn  chk3mm trlpick )))@>"0 argsan
  ('ztrmmrlcu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcu  chk3mm trlpick )))@>"0 argsan
  ('ztrmmrunn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunn  chk3mm trupick )))@>"0 argsan
  ('ztrmmrunu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunu  chk3mm trupick )))@>"0 argsan
  ('ztrmmrutn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutn  chk3mm trupick )))@>"0 argsan
  ('ztrmmrutu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutu  chk3mm trupick )))@>"0 argsan
  ('ztrmmrucn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucn  chk3mm trupick )))@>"0 argsan
  ('ztrmmrucu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucu  chk3mm trupick )))@>"0 argsan

  NB. note: chk2mm accepts input in format compatible with
  NB.       BLIS' trmm3, so let's employ it

  ('trmmllnnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnnn chk2mm trlpick )))@>"0 argsamm
  ('trmmllnnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnnt chk2mm trlpick )))@>"0 argsamn
  ('trmmllnnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnnj chk2mm trlpick )))@>"0 argsamm
  ('trmmllnnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnnc chk2mm trlpick )))@>"0 argsamn
  ('trmmllnun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnun chk2mm trl1pick)))@>"0 argsamm
  ('trmmllnut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnut chk2mm trl1pick)))@>"0 argsamn
  ('trmmllnuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnuj chk2mm trl1pick)))@>"0 argsamm
  ('trmmllnuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnuc chk2mm trl1pick)))@>"0 argsamn
  ('trmmlltnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltnn chk2mm trlpick )))@>"0 argsamm
  ('trmmlltnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltnt chk2mm trlpick )))@>"0 argsamn
  ('trmmlltnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltnj chk2mm trlpick )))@>"0 argsamm
  ('trmmlltnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltnc chk2mm trlpick )))@>"0 argsamn
  ('trmmlltun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltun chk2mm trl1pick)))@>"0 argsamm
  ('trmmlltut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltut chk2mm trl1pick)))@>"0 argsamn
  ('trmmlltuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltuj chk2mm trl1pick)))@>"0 argsamm
  ('trmmlltuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltuc chk2mm trl1pick)))@>"0 argsamn
  ('trmmlljnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlljnn chk2mm trlpick )))@>"0 argsamm
  ('trmmlljnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlljnt chk2mm trlpick )))@>"0 argsamn
  ('trmmlljnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlljnj chk2mm trlpick )))@>"0 argsamm
  ('trmmlljnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlljnc chk2mm trlpick )))@>"0 argsamn
  ('trmmlljun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlljun chk2mm trl1pick)))@>"0 argsamm
  ('trmmlljut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlljut chk2mm trl1pick)))@>"0 argsamn
  ('trmmlljuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlljuj chk2mm trl1pick)))@>"0 argsamm
  ('trmmlljuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlljuc chk2mm trl1pick)))@>"0 argsamn
  ('trmmllcnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcnn chk2mm trlpick )))@>"0 argsamm
  ('trmmllcnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcnt chk2mm trlpick )))@>"0 argsamn
  ('trmmllcnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcnj chk2mm trlpick )))@>"0 argsamm
  ('trmmllcnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcnc chk2mm trlpick )))@>"0 argsamn
  ('trmmllcun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcun chk2mm trl1pick)))@>"0 argsamm
  ('trmmllcut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcut chk2mm trl1pick)))@>"0 argsamn
  ('trmmllcuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcuj chk2mm trl1pick)))@>"0 argsamm
  ('trmmllcuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcuc chk2mm trl1pick)))@>"0 argsamn
  ('trmmlunnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunnn chk2mm trupick )))@>"0 argsamm
  ('trmmlunnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunnt chk2mm trupick )))@>"0 argsamn
  ('trmmlunnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunnj chk2mm trupick )))@>"0 argsamm
  ('trmmlunnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunnc chk2mm trupick )))@>"0 argsamn
  ('trmmlunun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunun chk2mm tru1pick)))@>"0 argsamm
  ('trmmlunut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunut chk2mm tru1pick)))@>"0 argsamn
  ('trmmlunuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunuj chk2mm tru1pick)))@>"0 argsamm
  ('trmmlunuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunuc chk2mm tru1pick)))@>"0 argsamn
  ('trmmlutnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutnn chk2mm trupick )))@>"0 argsamm
  ('trmmlutnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutnt chk2mm trupick )))@>"0 argsamn
  ('trmmlutnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutnj chk2mm trupick )))@>"0 argsamm
  ('trmmlutnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutnc chk2mm trupick )))@>"0 argsamn
  ('trmmlutun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutun chk2mm tru1pick)))@>"0 argsamm
  ('trmmlutut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutut chk2mm tru1pick)))@>"0 argsamn
  ('trmmlutuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutuj chk2mm tru1pick)))@>"0 argsamm
  ('trmmlutuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutuc chk2mm tru1pick)))@>"0 argsamn
  ('trmmlujnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlujnn chk2mm trupick )))@>"0 argsamm
  ('trmmlujnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlujnt chk2mm trupick )))@>"0 argsamn
  ('trmmlujnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlujnj chk2mm trupick )))@>"0 argsamm
  ('trmmlujnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlujnc chk2mm trupick )))@>"0 argsamn
  ('trmmlujun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlujun chk2mm tru1pick)))@>"0 argsamm
  ('trmmlujut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlujut chk2mm tru1pick)))@>"0 argsamn
  ('trmmlujuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlujuj chk2mm tru1pick)))@>"0 argsamm
  ('trmmlujuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlujuc chk2mm tru1pick)))@>"0 argsamn
  ('trmmlucnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucnn chk2mm trupick )))@>"0 argsamm
  ('trmmlucnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucnt chk2mm trupick )))@>"0 argsamn
  ('trmmlucnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucnj chk2mm trupick )))@>"0 argsamm
  ('trmmlucnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucnc chk2mm trupick )))@>"0 argsamn
  ('trmmlucun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucun chk2mm tru1pick)))@>"0 argsamm
  ('trmmlucut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucut chk2mm tru1pick)))@>"0 argsamn
  ('trmmlucuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucuj chk2mm tru1pick)))@>"0 argsamm
  ('trmmlucuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucuc chk2mm tru1pick)))@>"0 argsamn
  ('trmmrlnnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnnn chk2mm trlpick )))@>"0 argsanm
  ('trmmrlnnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnnt chk2mm trlpick )))@>"0 argsann
  ('trmmrlnnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnnj chk2mm trlpick )))@>"0 argsanm
  ('trmmrlnnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnnc chk2mm trlpick )))@>"0 argsann
  ('trmmrlnun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnun chk2mm trl1pick)))@>"0 argsanm
  ('trmmrlnut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnut chk2mm trl1pick)))@>"0 argsann
  ('trmmrlnuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnuj chk2mm trl1pick)))@>"0 argsanm
  ('trmmrlnuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnuc chk2mm trl1pick)))@>"0 argsann
  ('trmmrltnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltnn chk2mm trlpick )))@>"0 argsanm
  ('trmmrltnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltnt chk2mm trlpick )))@>"0 argsann
  ('trmmrltnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltnj chk2mm trlpick )))@>"0 argsanm
  ('trmmrltnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltnc chk2mm trlpick )))@>"0 argsann
  ('trmmrltun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltun chk2mm trl1pick)))@>"0 argsanm
  ('trmmrltut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltut chk2mm trl1pick)))@>"0 argsann
  ('trmmrltuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltuj chk2mm trl1pick)))@>"0 argsanm
  ('trmmrltuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltuc chk2mm trl1pick)))@>"0 argsann
  ('trmmrljnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrljnn chk2mm trlpick )))@>"0 argsanm
  ('trmmrljnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrljnt chk2mm trlpick )))@>"0 argsann
  ('trmmrljnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrljnj chk2mm trlpick )))@>"0 argsanm
  ('trmmrljnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrljnc chk2mm trlpick )))@>"0 argsann
  ('trmmrljun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrljun chk2mm trl1pick)))@>"0 argsanm
  ('trmmrljut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrljut chk2mm trl1pick)))@>"0 argsann
  ('trmmrljuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrljuj chk2mm trl1pick)))@>"0 argsanm
  ('trmmrljuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrljuc chk2mm trl1pick)))@>"0 argsann
  ('trmmrlcnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcnn chk2mm trlpick )))@>"0 argsanm
  ('trmmrlcnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcnt chk2mm trlpick )))@>"0 argsann
  ('trmmrlcnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcnj chk2mm trlpick )))@>"0 argsanm
  ('trmmrlcnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcnc chk2mm trlpick )))@>"0 argsann
  ('trmmrlcun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcun chk2mm trl1pick)))@>"0 argsanm
  ('trmmrlcut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcut chk2mm trl1pick)))@>"0 argsann
  ('trmmrlcuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcuj chk2mm trl1pick)))@>"0 argsanm
  ('trmmrlcuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcuc chk2mm trl1pick)))@>"0 argsann
  ('trmmrunnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunnn chk2mm trupick )))@>"0 argsanm
  ('trmmrunnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunnt chk2mm trupick )))@>"0 argsann
  ('trmmrunnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunnj chk2mm trupick )))@>"0 argsanm
  ('trmmrunnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunnc chk2mm trupick )))@>"0 argsann
  ('trmmrunun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunun chk2mm tru1pick)))@>"0 argsanm
  ('trmmrunut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunut chk2mm tru1pick)))@>"0 argsann
  ('trmmrunuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunuj chk2mm tru1pick)))@>"0 argsanm
  ('trmmrunuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunuc chk2mm tru1pick)))@>"0 argsann
  ('trmmrutnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutnn chk2mm trupick )))@>"0 argsanm
  ('trmmrutnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutnt chk2mm trupick )))@>"0 argsann
  ('trmmrutnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutnj chk2mm trupick )))@>"0 argsanm
  ('trmmrutnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutnc chk2mm trupick )))@>"0 argsann
  ('trmmrutun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutun chk2mm tru1pick)))@>"0 argsanm
  ('trmmrutut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutut chk2mm tru1pick)))@>"0 argsann
  ('trmmrutuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutuj chk2mm tru1pick)))@>"0 argsanm
  ('trmmrutuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutuc chk2mm tru1pick)))@>"0 argsann
  ('trmmrujnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrujnn chk2mm trupick )))@>"0 argsanm
  ('trmmrujnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrujnt chk2mm trupick )))@>"0 argsann
  ('trmmrujnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrujnj chk2mm trupick )))@>"0 argsanm
  ('trmmrujnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrujnc chk2mm trupick )))@>"0 argsann
  ('trmmrujun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrujun chk2mm tru1pick)))@>"0 argsanm
  ('trmmrujut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrujut chk2mm tru1pick)))@>"0 argsann
  ('trmmrujuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrujuj chk2mm tru1pick)))@>"0 argsanm
  ('trmmrujuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrujuc chk2mm tru1pick)))@>"0 argsann
  ('trmmrucnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucnn chk2mm trupick )))@>"0 argsanm
  ('trmmrucnt_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucnt chk2mm trupick )))@>"0 argsann
  ('trmmrucnj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucnj chk2mm trupick )))@>"0 argsanm
  ('trmmrucnc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucnc chk2mm trupick )))@>"0 argsann
  ('trmmrucun_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucun chk2mm tru1pick)))@>"0 argsanm
  ('trmmrucut_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucut chk2mm tru1pick)))@>"0 argsann
  ('trmmrucuj_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucuj chk2mm tru1pick)))@>"0 argsanm
  ('trmmrucuc_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucuc chk2mm tru1pick)))@>"0 argsann

  ('dtrmmllnnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnnn chk2mm trlpick )))@>"0 argsdmm
  ('dtrmmllnnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnnt chk2mm trlpick )))@>"0 argsdmn
  ('dtrmmllnun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnun chk2mm trl1pick)))@>"0 argsdmm
  ('dtrmmllnut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnut chk2mm trl1pick)))@>"0 argsdmn
  ('dtrmmlltnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltnn chk2mm trlpick )))@>"0 argsdmm
  ('dtrmmlltnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltnt chk2mm trlpick )))@>"0 argsdmn
  ('dtrmmlltun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltun chk2mm trl1pick)))@>"0 argsdmm
  ('dtrmmlltut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltut chk2mm trl1pick)))@>"0 argsdmn
  ('dtrmmlunnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunnn chk2mm trupick )))@>"0 argsdmm
  ('dtrmmlunnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunnt chk2mm trupick )))@>"0 argsdmn
  ('dtrmmlunun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunun chk2mm tru1pick)))@>"0 argsdmm
  ('dtrmmlunut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunut chk2mm tru1pick)))@>"0 argsdmn
  ('dtrmmlutnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutnn chk2mm trupick )))@>"0 argsdmm
  ('dtrmmlutnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutnt chk2mm trupick )))@>"0 argsdmn
  ('dtrmmlutun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutun chk2mm tru1pick)))@>"0 argsdmm
  ('dtrmmlutut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutut chk2mm tru1pick)))@>"0 argsdmn
  ('dtrmmrlnnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnnn chk2mm trlpick )))@>"0 argsdnm
  ('dtrmmrlnnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnnt chk2mm trlpick )))@>"0 argsdnn
  ('dtrmmrlnun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnun chk2mm trl1pick)))@>"0 argsdnm
  ('dtrmmrlnut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnut chk2mm trl1pick)))@>"0 argsdnn
  ('dtrmmrltnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltnn chk2mm trlpick )))@>"0 argsdnm
  ('dtrmmrltnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltnt chk2mm trlpick )))@>"0 argsdnn
  ('dtrmmrltun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltun chk2mm trl1pick)))@>"0 argsdnm
  ('dtrmmrltut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltut chk2mm trl1pick)))@>"0 argsdnn
  ('dtrmmrunnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunnn chk2mm trupick )))@>"0 argsdnm
  ('dtrmmrunnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunnt chk2mm trupick )))@>"0 argsdnn
  ('dtrmmrunun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunun chk2mm tru1pick)))@>"0 argsdnm
  ('dtrmmrunut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunut chk2mm tru1pick)))@>"0 argsdnn
  ('dtrmmrutnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutnn chk2mm trupick )))@>"0 argsdnm
  ('dtrmmrutnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutnt chk2mm trupick )))@>"0 argsdnn
  ('dtrmmrutun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutun chk2mm tru1pick)))@>"0 argsdnm
  ('dtrmmrutut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutut chk2mm tru1pick)))@>"0 argsdnn

  ('ztrmmllnnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnnn chk2mm trlpick )))@>"0 argszmm
  ('ztrmmllnnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnnt chk2mm trlpick )))@>"0 argszmn
  ('ztrmmllnnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnnj chk2mm trlpick )))@>"0 argszmm
  ('ztrmmllnnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnnc chk2mm trlpick )))@>"0 argszmn
  ('ztrmmllnun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnun chk2mm trl1pick)))@>"0 argszmm
  ('ztrmmllnut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnut chk2mm trl1pick)))@>"0 argszmn
  ('ztrmmllnuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnuj chk2mm trl1pick)))@>"0 argszmm
  ('ztrmmllnuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnuc chk2mm trl1pick)))@>"0 argszmn
  ('ztrmmlltnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltnn chk2mm trlpick )))@>"0 argszmm
  ('ztrmmlltnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltnt chk2mm trlpick )))@>"0 argszmn
  ('ztrmmlltnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltnj chk2mm trlpick )))@>"0 argszmm
  ('ztrmmlltnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltnc chk2mm trlpick )))@>"0 argszmn
  ('ztrmmlltun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltun chk2mm trl1pick)))@>"0 argszmm
  ('ztrmmlltut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltut chk2mm trl1pick)))@>"0 argszmn
  ('ztrmmlltuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltuj chk2mm trl1pick)))@>"0 argszmm
  ('ztrmmlltuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltuc chk2mm trl1pick)))@>"0 argszmn
  ('ztrmmlljnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlljnn chk2mm trlpick )))@>"0 argszmm
  ('ztrmmlljnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlljnt chk2mm trlpick )))@>"0 argszmn
  ('ztrmmlljnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlljnj chk2mm trlpick )))@>"0 argszmm
  ('ztrmmlljnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlljnc chk2mm trlpick )))@>"0 argszmn
  ('ztrmmlljun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlljun chk2mm trl1pick)))@>"0 argszmm
  ('ztrmmlljut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlljut chk2mm trl1pick)))@>"0 argszmn
  ('ztrmmlljuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlljuj chk2mm trl1pick)))@>"0 argszmm
  ('ztrmmlljuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlljuc chk2mm trl1pick)))@>"0 argszmn
  ('ztrmmllcnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllcnn chk2mm trlpick )))@>"0 argszmm
  ('ztrmmllcnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllcnt chk2mm trlpick )))@>"0 argszmn
  ('ztrmmllcnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllcnj chk2mm trlpick )))@>"0 argszmm
  ('ztrmmllcnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllcnc chk2mm trlpick )))@>"0 argszmn
  ('ztrmmllcun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllcun chk2mm trl1pick)))@>"0 argszmm
  ('ztrmmllcut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllcut chk2mm trl1pick)))@>"0 argszmn
  ('ztrmmllcuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllcuj chk2mm trl1pick)))@>"0 argszmm
  ('ztrmmllcuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllcuc chk2mm trl1pick)))@>"0 argszmn
  ('ztrmmlunnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunnn chk2mm trupick )))@>"0 argszmm
  ('ztrmmlunnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunnt chk2mm trupick )))@>"0 argszmn
  ('ztrmmlunnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunnj chk2mm trupick )))@>"0 argszmm
  ('ztrmmlunnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunnc chk2mm trupick )))@>"0 argszmn
  ('ztrmmlunun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunun chk2mm tru1pick)))@>"0 argszmm
  ('ztrmmlunut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunut chk2mm tru1pick)))@>"0 argszmn
  ('ztrmmlunuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunuj chk2mm tru1pick)))@>"0 argszmm
  ('ztrmmlunuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunuc chk2mm tru1pick)))@>"0 argszmn
  ('ztrmmlutnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutnn chk2mm trupick )))@>"0 argszmm
  ('ztrmmlutnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutnt chk2mm trupick )))@>"0 argszmn
  ('ztrmmlutnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutnj chk2mm trupick )))@>"0 argszmm
  ('ztrmmlutnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutnc chk2mm trupick )))@>"0 argszmn
  ('ztrmmlutun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutun chk2mm tru1pick)))@>"0 argszmm
  ('ztrmmlutut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutut chk2mm tru1pick)))@>"0 argszmn
  ('ztrmmlutuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutuj chk2mm tru1pick)))@>"0 argszmm
  ('ztrmmlutuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutuc chk2mm tru1pick)))@>"0 argszmn
  ('ztrmmlujnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlujnn chk2mm trupick )))@>"0 argszmm
  ('ztrmmlujnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlujnt chk2mm trupick )))@>"0 argszmn
  ('ztrmmlujnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlujnj chk2mm trupick )))@>"0 argszmm
  ('ztrmmlujnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlujnc chk2mm trupick )))@>"0 argszmn
  ('ztrmmlujun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlujun chk2mm tru1pick)))@>"0 argszmm
  ('ztrmmlujut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlujut chk2mm tru1pick)))@>"0 argszmn
  ('ztrmmlujuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlujuj chk2mm tru1pick)))@>"0 argszmm
  ('ztrmmlujuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlujuc chk2mm tru1pick)))@>"0 argszmn
  ('ztrmmlucnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlucnn chk2mm trupick )))@>"0 argszmm
  ('ztrmmlucnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlucnt chk2mm trupick )))@>"0 argszmn
  ('ztrmmlucnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlucnj chk2mm trupick )))@>"0 argszmm
  ('ztrmmlucnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlucnc chk2mm trupick )))@>"0 argszmn
  ('ztrmmlucun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlucun chk2mm tru1pick)))@>"0 argszmm
  ('ztrmmlucut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlucut chk2mm tru1pick)))@>"0 argszmn
  ('ztrmmlucuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlucuj chk2mm tru1pick)))@>"0 argszmm
  ('ztrmmlucuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlucuc chk2mm tru1pick)))@>"0 argszmn
  ('ztrmmrlnnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnnn chk2mm trlpick )))@>"0 argsznm
  ('ztrmmrlnnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnnt chk2mm trlpick )))@>"0 argsznn
  ('ztrmmrlnnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnnj chk2mm trlpick )))@>"0 argsznm
  ('ztrmmrlnnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnnc chk2mm trlpick )))@>"0 argsznn
  ('ztrmmrlnun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnun chk2mm trl1pick)))@>"0 argsznm
  ('ztrmmrlnut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnut chk2mm trl1pick)))@>"0 argsznn
  ('ztrmmrlnuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnuj chk2mm trl1pick)))@>"0 argsznm
  ('ztrmmrlnuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnuc chk2mm trl1pick)))@>"0 argsznn
  ('ztrmmrltnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltnn chk2mm trlpick )))@>"0 argsznm
  ('ztrmmrltnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltnt chk2mm trlpick )))@>"0 argsznn
  ('ztrmmrltnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltnj chk2mm trlpick )))@>"0 argsznm
  ('ztrmmrltnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltnc chk2mm trlpick )))@>"0 argsznn
  ('ztrmmrltun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltun chk2mm trl1pick)))@>"0 argsznm
  ('ztrmmrltut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltut chk2mm trl1pick)))@>"0 argsznn
  ('ztrmmrltuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltuj chk2mm trl1pick)))@>"0 argsznm
  ('ztrmmrltuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltuc chk2mm trl1pick)))@>"0 argsznn
  ('ztrmmrljnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrljnn chk2mm trlpick )))@>"0 argsznm
  ('ztrmmrljnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrljnt chk2mm trlpick )))@>"0 argsznn
  ('ztrmmrljnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrljnj chk2mm trlpick )))@>"0 argsznm
  ('ztrmmrljnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrljnc chk2mm trlpick )))@>"0 argsznn
  ('ztrmmrljun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrljun chk2mm trl1pick)))@>"0 argsznm
  ('ztrmmrljut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrljut chk2mm trl1pick)))@>"0 argsznn
  ('ztrmmrljuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrljuj chk2mm trl1pick)))@>"0 argsznm
  ('ztrmmrljuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrljuc chk2mm trl1pick)))@>"0 argsznn
  ('ztrmmrlcnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcnn chk2mm trlpick )))@>"0 argsznm
  ('ztrmmrlcnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcnt chk2mm trlpick )))@>"0 argsznn
  ('ztrmmrlcnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcnj chk2mm trlpick )))@>"0 argsznm
  ('ztrmmrlcnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcnc chk2mm trlpick )))@>"0 argsznn
  ('ztrmmrlcun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcun chk2mm trl1pick)))@>"0 argsznm
  ('ztrmmrlcut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcut chk2mm trl1pick)))@>"0 argsznn
  ('ztrmmrlcuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcuj chk2mm trl1pick)))@>"0 argsznm
  ('ztrmmrlcuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcuc chk2mm trl1pick)))@>"0 argsznn
  ('ztrmmrunnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunnn chk2mm trupick )))@>"0 argsznm
  ('ztrmmrunnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunnt chk2mm trupick )))@>"0 argsznn
  ('ztrmmrunnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunnj chk2mm trupick )))@>"0 argsznm
  ('ztrmmrunnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunnc chk2mm trupick )))@>"0 argsznn
  ('ztrmmrunun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunun chk2mm tru1pick)))@>"0 argsznm
  ('ztrmmrunut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunut chk2mm tru1pick)))@>"0 argsznn
  ('ztrmmrunuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunuj chk2mm tru1pick)))@>"0 argsznm
  ('ztrmmrunuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunuc chk2mm tru1pick)))@>"0 argsznn
  ('ztrmmrutnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutnn chk2mm trupick )))@>"0 argsznm
  ('ztrmmrutnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutnt chk2mm trupick )))@>"0 argsznn
  ('ztrmmrutnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutnj chk2mm trupick )))@>"0 argsznm
  ('ztrmmrutnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutnc chk2mm trupick )))@>"0 argsznn
  ('ztrmmrutun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutun chk2mm tru1pick)))@>"0 argsznm
  ('ztrmmrutut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutut chk2mm tru1pick)))@>"0 argsznn
  ('ztrmmrutuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutuj chk2mm tru1pick)))@>"0 argsznm
  ('ztrmmrutuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutuc chk2mm tru1pick)))@>"0 argsznn
  ('ztrmmrujnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrujnn chk2mm trupick )))@>"0 argsznm
  ('ztrmmrujnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrujnt chk2mm trupick )))@>"0 argsznn
  ('ztrmmrujnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrujnj chk2mm trupick )))@>"0 argsznm
  ('ztrmmrujnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrujnc chk2mm trupick )))@>"0 argsznn
  ('ztrmmrujun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrujun chk2mm tru1pick)))@>"0 argsznm
  ('ztrmmrujut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrujut chk2mm tru1pick)))@>"0 argsznn
  ('ztrmmrujuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrujuj chk2mm tru1pick)))@>"0 argsznm
  ('ztrmmrujuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrujuc chk2mm tru1pick)))@>"0 argsznn
  ('ztrmmrucnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrucnn chk2mm trupick )))@>"0 argsznm
  ('ztrmmrucnt_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrucnt chk2mm trupick )))@>"0 argsznn
  ('ztrmmrucnj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrucnj chk2mm trupick )))@>"0 argsznm
  ('ztrmmrucnc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrucnc chk2mm trupick )))@>"0 argsznn
  ('ztrmmrucun_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrucun chk2mm tru1pick)))@>"0 argsznm
  ('ztrmmrucut_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrucut chk2mm tru1pick)))@>"0 argsznn
  ('ztrmmrucuj_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrucuj chk2mm tru1pick)))@>"0 argsznm
  ('ztrmmrucuc_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrucuc chk2mm tru1pick)))@>"0 argsznn

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbasicmm
NB.
NB. Description:
NB.   Adv. to make verb to test basic matrix-matrix operations
NB.
NB. Syntax:
NB.   vtest=. mkmat testbasicmm
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms; is called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testbasicmm_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testbasicmm_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbasicmm_mt_ 150 200

testbasicmm=: 1 : 'EMPTY [ testbasictrmm_mt_@(u@(2 # >./) ; u ; u) [ testbasichemm_mt_@((9&o. upddiag_mt_)@u@(2 # >./) ; u ; u) [ testbasicsymm_mt_@(u@(2 # >./) ; u ; u) [ testbasicgemmt_mt_@(u@(+/\) ; u@|.@(+/\) ; u@(2 # {.)) [ testbasicgemm_mt_@(u@(+/\) ; u@(+/\.) ; u) [ load@''math/mt/test/blis/mm'' [ load@''math/mt/test/blas/mm'''

NB. ---------------------------------------------------------
NB. testbasictrsv
NB.
NB. Description:
NB.   Test equation solvers:
NB.   - xTRSV (BLAS)
NB.   by triangular matrix and vector
NB.
NB. Syntax:
NB.   testbasictrsv AA ; b
NB. where
NB.   AA - n×n-matrix, A material
NB.   x  - n-vector, the RHS

testbasictrsv=: 3 : 0
  inc=. 1 2 _1 _2
  'AA y'=. y
  args=. (1 reverse 2)@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc  NB. 4×3-matrix of boxes, each row is argument to tmonad

  NB. for every i feed the tuple (AA ; expanded_b_i ; incb_i) to tmonad
  ('dtrsvlnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnn chk3sv)))"1 args
  ('dtrsvlnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnu chk3sv)))"1 args
  ('dtrsvltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltn chk3sv)))"1 args
  ('dtrsvltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltu chk3sv)))"1 args
  ('dtrsvunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunn chk3sv)))"1 args
  ('dtrsvunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunu chk3sv)))"1 args
  ('dtrsvutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutn chk3sv)))"1 args
  ('dtrsvutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutu chk3sv)))"1 args
  ('ztrsvlnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnn chk3sv)))"1 args
  ('ztrsvlnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnu chk3sv)))"1 args
  ('ztrsvltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltn chk3sv)))"1 args
  ('ztrsvltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltu chk3sv)))"1 args
  ('ztrsvlcn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlcn chk3sv)))"1 args
  ('ztrsvlcu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlcu chk3sv)))"1 args
  ('ztrsvunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunn chk3sv)))"1 args
  ('ztrsvunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunu chk3sv)))"1 args
  ('ztrsvutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutn chk3sv)))"1 args
  ('ztrsvutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutu chk3sv)))"1 args
  ('ztrsvucn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvucn chk3sv)))"1 args
  ('ztrsvucu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvucu chk3sv)))"1 args

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbasicsv
NB.
NB. Description:
NB.   Adv. to make verb to test equation solvers
NB.
NB. Syntax:
NB.   vtest=. mkmat testbasicsv
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (n,n)
NB.   vtest - monad to test algorithms; is called as:
NB.             vtest (n,n)
NB.   (n,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testbasicsv_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testbasicsv_mt_ 200 200
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbasicsv_mt_ 250 250

testbasicsv=: 1 : 'EMPTY [ testbasictrsv_mt_@(u ; u@{.) [ load@''math/mt/test/blas/sv'''

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
NB.   testbasictrsm AA ; B
NB. where
NB.   AA - k×k-matrix, A material
NB.   B  - m×n-matrix, RHS
NB.   k  = max(m,n)
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

  ('dtrsmllnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllnn chk3sm)))@>"0 argsdm
  ('dtrsmllnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllnu chk3sm)))@>"0 argsdm
  ('dtrsmlltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlltn chk3sm)))@>"0 argsdm
  ('dtrsmlltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlltu chk3sm)))@>"0 argsdm
  ('dtrsmlunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlunn chk3sm)))@>"0 argsdm
  ('dtrsmlunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlunu chk3sm)))@>"0 argsdm
  ('dtrsmlutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlutn chk3sm)))@>"0 argsdm
  ('dtrsmlutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlutu chk3sm)))@>"0 argsdm
  ('dtrsmrlnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnn chk3sm)))@>"0 argsdn
  ('dtrsmrlnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnu chk3sm)))@>"0 argsdn
  ('dtrsmrltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrltn chk3sm)))@>"0 argsdn
  ('dtrsmrltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrltu chk3sm)))@>"0 argsdn
  ('dtrsmrunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrunn chk3sm)))@>"0 argsdn
  ('dtrsmrunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrunu chk3sm)))@>"0 argsdn
  ('dtrsmrutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrutn chk3sm)))@>"0 argsdn
  ('dtrsmrutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrutu chk3sm)))@>"0 argsdn

  ('ztrsmllnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllnn chk3sm)))@>"0 argszm
  ('ztrsmllnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllnu chk3sm)))@>"0 argszm
  ('ztrsmlltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlltn chk3sm)))@>"0 argszm
  ('ztrsmlltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlltu chk3sm)))@>"0 argszm
  ('ztrsmllcn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllcn chk3sm)))@>"0 argszm
  ('ztrsmllcu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllcu chk3sm)))@>"0 argszm
  ('ztrsmlunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlunn chk3sm)))@>"0 argszm
  ('ztrsmlunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlunu chk3sm)))@>"0 argszm
  ('ztrsmlutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlutn chk3sm)))@>"0 argszm
  ('ztrsmlutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlutu chk3sm)))@>"0 argszm
  ('ztrsmlucn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlucn chk3sm)))@>"0 argszm
  ('ztrsmlucu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlucu chk3sm)))@>"0 argszm
  ('ztrsmrlnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnn chk3sm)))@>"0 argszn
  ('ztrsmrlnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnu chk3sm)))@>"0 argszn
  ('ztrsmrltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrltn chk3sm)))@>"0 argszn
  ('ztrsmrltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrltu chk3sm)))@>"0 argszn
  ('ztrsmrlcn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcn chk3sm)))@>"0 argszn
  ('ztrsmrlcu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcu chk3sm)))@>"0 argszn
  ('ztrsmrunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrunn chk3sm)))@>"0 argszn
  ('ztrsmrunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrunu chk3sm)))@>"0 argszn
  ('ztrsmrutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrutn chk3sm)))@>"0 argszn
  ('ztrsmrutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrutu chk3sm)))@>"0 argszn
  ('ztrsmrucn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrucn chk3sm)))@>"0 argszn
  ('ztrsmrucu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrucu chk3sm)))@>"0 argszn

  NB. BLIS' staff

  ('trsmllnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnn chk3sm)))@>"0 argsam
  ('trsmllnu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllnu chk3sm)))@>"0 argsam
  ('trsmlltn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltn chk3sm)))@>"0 argsam
  ('trsmlltu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlltu chk3sm)))@>"0 argsam
  ('trsmlljn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlljn chk3sm)))@>"0 argsam
  ('trsmllju_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllju chk3sm)))@>"0 argsam
  ('trsmllcn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcn chk3sm)))@>"0 argsam
  ('trsmllcu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmllcu chk3sm)))@>"0 argsam
  ('trsmlunn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunn chk3sm)))@>"0 argsam
  ('trsmlunu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlunu chk3sm)))@>"0 argsam
  ('trsmlutn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutn chk3sm)))@>"0 argsam
  ('trsmlutu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlutu chk3sm)))@>"0 argsam
  ('trsmlujn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlujn chk3sm)))@>"0 argsam
  ('trsmluju_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmluju chk3sm)))@>"0 argsam
  ('trsmlucn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucn chk3sm)))@>"0 argsam
  ('trsmlucu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmlucu chk3sm)))@>"0 argsam
  ('trsmrlnn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnn chk3sm)))@>"0 argsan
  ('trsmrlnu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlnu chk3sm)))@>"0 argsan
  ('trsmrltn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltn chk3sm)))@>"0 argsan
  ('trsmrltu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrltu chk3sm)))@>"0 argsan
  ('trsmrljn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrljn chk3sm)))@>"0 argsan
  ('trsmrlju_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlju chk3sm)))@>"0 argsan
  ('trsmrlcn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcn chk3sm)))@>"0 argsan
  ('trsmrlcu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrlcu chk3sm)))@>"0 argsan
  ('trsmrunn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunn chk3sm)))@>"0 argsan
  ('trsmrunu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrunu chk3sm)))@>"0 argsan
  ('trsmrutn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutn chk3sm)))@>"0 argsan
  ('trsmrutu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrutu chk3sm)))@>"0 argsan
  ('trsmrujn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrujn chk3sm)))@>"0 argsan
  ('trsmruju_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmruju chk3sm)))@>"0 argsan
  ('trsmrucn_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucn chk3sm)))@>"0 argsan
  ('trsmrucu_mtbli_'  tmonad (]`]`(_."_)`(_."_)`(trmmrucu chk3sm)))@>"0 argsan

  ('dtrsmllnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnn chk3sm)))@>"0 argsdm
  ('dtrsmllnu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnu chk3sm)))@>"0 argsdm
  ('dtrsmlltn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltn chk3sm)))@>"0 argsdm
  ('dtrsmlltu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltu chk3sm)))@>"0 argsdm
  ('dtrsmlunn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunn chk3sm)))@>"0 argsdm
  ('dtrsmlunu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunu chk3sm)))@>"0 argsdm
  ('dtrsmlutn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutn chk3sm)))@>"0 argsdm
  ('dtrsmlutu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutu chk3sm)))@>"0 argsdm
  ('dtrsmrlnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnn chk3sm)))@>"0 argsdn
  ('dtrsmrlnu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnu chk3sm)))@>"0 argsdn
  ('dtrsmrltn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltn chk3sm)))@>"0 argsdn
  ('dtrsmrltu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltu chk3sm)))@>"0 argsdn
  ('dtrsmrunn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunn chk3sm)))@>"0 argsdn
  ('dtrsmrunu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunu chk3sm)))@>"0 argsdn
  ('dtrsmrutn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutn chk3sm)))@>"0 argsdn
  ('dtrsmrutu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutu chk3sm)))@>"0 argsdn

  ('ztrsmllnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnn chk3sm)))@>"0 argsam
  ('ztrsmllnu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllnu chk3sm)))@>"0 argsam
  ('ztrsmlltn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltn chk3sm)))@>"0 argsam
  ('ztrsmlltu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlltu chk3sm)))@>"0 argsam
  ('ztrsmlljn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlljn chk3sm)))@>"0 argsam
  ('ztrsmllju_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllju chk3sm)))@>"0 argsam
  ('ztrsmllcn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllcn chk3sm)))@>"0 argsam
  ('ztrsmllcu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmllcu chk3sm)))@>"0 argsam
  ('ztrsmlunn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunn chk3sm)))@>"0 argsam
  ('ztrsmlunu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlunu chk3sm)))@>"0 argsam
  ('ztrsmlutn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutn chk3sm)))@>"0 argsam
  ('ztrsmlutu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlutu chk3sm)))@>"0 argsam
  ('ztrsmlujn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlujn chk3sm)))@>"0 argsam
  ('ztrsmluju_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmluju chk3sm)))@>"0 argsam
  ('ztrsmlucn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlucn chk3sm)))@>"0 argsam
  ('ztrsmlucu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmlucu chk3sm)))@>"0 argsam
  ('ztrsmrlnn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnn chk3sm)))@>"0 argsan
  ('ztrsmrlnu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnu chk3sm)))@>"0 argsan
  ('ztrsmrltn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltn chk3sm)))@>"0 argsan
  ('ztrsmrltu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrltu chk3sm)))@>"0 argsan
  ('ztrsmrljn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrljn chk3sm)))@>"0 argsan
  ('ztrsmrlju_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlju chk3sm)))@>"0 argsan
  ('ztrsmrlcn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcn chk3sm)))@>"0 argsan
  ('ztrsmrlcu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcu chk3sm)))@>"0 argsan
  ('ztrsmrunn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunn chk3sm)))@>"0 argsan
  ('ztrsmrunu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrunu chk3sm)))@>"0 argsan
  ('ztrsmrutn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutn chk3sm)))@>"0 argsan
  ('ztrsmrutu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrutu chk3sm)))@>"0 argsan
  ('ztrsmrujn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrujn chk3sm)))@>"0 argsan
  ('ztrsmruju_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmruju chk3sm)))@>"0 argsan
  ('ztrsmrucn_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrucn chk3sm)))@>"0 argsan
  ('ztrsmrucu_mtbli_' tmonad (]`]`(_."_)`(_."_)`(trmmrucu chk3sm)))@>"0 argsan

  NB. mt staff

  NB. monadic trsmxxxx, 1-rank b and x
  ('trsmllnn'         tmonad (]`]`(_."_)`(_."_)`(trmmllnn chk3sm)))@>"0 argsbm
  ('trsmllnu'         tmonad (]`]`(_."_)`(_."_)`(trmmllnu chk3sm)))@>"0 argsbm
  ('trsmlltn'         tmonad (]`]`(_."_)`(_."_)`(trmmlltn chk3sm)))@>"0 argsbm
  ('trsmlltu'         tmonad (]`]`(_."_)`(_."_)`(trmmlltu chk3sm)))@>"0 argsbm
  ('trsmlljn'         tmonad (]`]`(_."_)`(_."_)`(trmmlljn chk3sm)))@>"0 argsbm
  ('trsmllju'         tmonad (]`]`(_."_)`(_."_)`(trmmllju chk3sm)))@>"0 argsbm
  ('trsmllcn'         tmonad (]`]`(_."_)`(_."_)`(trmmllcn chk3sm)))@>"0 argsbm
  ('trsmllcu'         tmonad (]`]`(_."_)`(_."_)`(trmmllcu chk3sm)))@>"0 argsbm
  ('trsmlunn'         tmonad (]`]`(_."_)`(_."_)`(trmmlunn chk3sm)))@>"0 argsbm
  ('trsmlunu'         tmonad (]`]`(_."_)`(_."_)`(trmmlunu chk3sm)))@>"0 argsbm
  ('trsmlutn'         tmonad (]`]`(_."_)`(_."_)`(trmmlutn chk3sm)))@>"0 argsbm
  ('trsmlutu'         tmonad (]`]`(_."_)`(_."_)`(trmmlutu chk3sm)))@>"0 argsbm
  ('trsmlujn'         tmonad (]`]`(_."_)`(_."_)`(trmmlujn chk3sm)))@>"0 argsbm
  ('trsmluju'         tmonad (]`]`(_."_)`(_."_)`(trmmluju chk3sm)))@>"0 argsbm
  ('trsmlucn'         tmonad (]`]`(_."_)`(_."_)`(trmmlucn chk3sm)))@>"0 argsbm
  ('trsmlucu'         tmonad (]`]`(_."_)`(_."_)`(trmmlucu chk3sm)))@>"0 argsbm
  ('trsmrlnn'         tmonad (]`]`(_."_)`(_."_)`(trmmrlnn chk3sm)))@>"0 argsbn
  ('trsmrlnu'         tmonad (]`]`(_."_)`(_."_)`(trmmrlnu chk3sm)))@>"0 argsbn
  ('trsmrltn'         tmonad (]`]`(_."_)`(_."_)`(trmmrltn chk3sm)))@>"0 argsbn
  ('trsmrltu'         tmonad (]`]`(_."_)`(_."_)`(trmmrltu chk3sm)))@>"0 argsbn
  ('trsmrljn'         tmonad (]`]`(_."_)`(_."_)`(trmmrljn chk3sm)))@>"0 argsbn
  ('trsmrlju'         tmonad (]`]`(_."_)`(_."_)`(trmmrlju chk3sm)))@>"0 argsbn
  ('trsmrlcn'         tmonad (]`]`(_."_)`(_."_)`(trmmrlcn chk3sm)))@>"0 argsbn
  ('trsmrlcu'         tmonad (]`]`(_."_)`(_."_)`(trmmrlcu chk3sm)))@>"0 argsbn
  ('trsmrunn'         tmonad (]`]`(_."_)`(_."_)`(trmmrunn chk3sm)))@>"0 argsbn
  ('trsmrunu'         tmonad (]`]`(_."_)`(_."_)`(trmmrunu chk3sm)))@>"0 argsbn
  ('trsmrutn'         tmonad (]`]`(_."_)`(_."_)`(trmmrutn chk3sm)))@>"0 argsbn
  ('trsmrutu'         tmonad (]`]`(_."_)`(_."_)`(trmmrutu chk3sm)))@>"0 argsbn
  ('trsmrujn'         tmonad (]`]`(_."_)`(_."_)`(trmmrujn chk3sm)))@>"0 argsbn
  ('trsmruju'         tmonad (]`]`(_."_)`(_."_)`(trmmruju chk3sm)))@>"0 argsbn
  ('trsmrucn'         tmonad (]`]`(_."_)`(_."_)`(trmmrucn chk3sm)))@>"0 argsbn
  ('trsmrucu'         tmonad (]`]`(_."_)`(_."_)`(trmmrucu chk3sm)))@>"0 argsbn

  NB. monadic trsmxxxx, 2-rank B and X
  ('trsmllnn'         tmonad (]`]`(_."_)`(_."_)`(trmmllnn chk3sm)))@>"0 argsam
  ('trsmllnu'         tmonad (]`]`(_."_)`(_."_)`(trmmllnu chk3sm)))@>"0 argsam
  ('trsmlltn'         tmonad (]`]`(_."_)`(_."_)`(trmmlltn chk3sm)))@>"0 argsam
  ('trsmlltu'         tmonad (]`]`(_."_)`(_."_)`(trmmlltu chk3sm)))@>"0 argsam
  ('trsmlljn'         tmonad (]`]`(_."_)`(_."_)`(trmmlljn chk3sm)))@>"0 argsam
  ('trsmllju'         tmonad (]`]`(_."_)`(_."_)`(trmmllju chk3sm)))@>"0 argsam
  ('trsmllcn'         tmonad (]`]`(_."_)`(_."_)`(trmmllcn chk3sm)))@>"0 argsam
  ('trsmllcu'         tmonad (]`]`(_."_)`(_."_)`(trmmllcu chk3sm)))@>"0 argsam
  ('trsmlunn'         tmonad (]`]`(_."_)`(_."_)`(trmmlunn chk3sm)))@>"0 argsam
  ('trsmlunu'         tmonad (]`]`(_."_)`(_."_)`(trmmlunu chk3sm)))@>"0 argsam
  ('trsmlutn'         tmonad (]`]`(_."_)`(_."_)`(trmmlutn chk3sm)))@>"0 argsam
  ('trsmlutu'         tmonad (]`]`(_."_)`(_."_)`(trmmlutu chk3sm)))@>"0 argsam
  ('trsmlujn'         tmonad (]`]`(_."_)`(_."_)`(trmmlujn chk3sm)))@>"0 argsam
  ('trsmluju'         tmonad (]`]`(_."_)`(_."_)`(trmmluju chk3sm)))@>"0 argsam
  ('trsmlucn'         tmonad (]`]`(_."_)`(_."_)`(trmmlucn chk3sm)))@>"0 argsam
  ('trsmlucu'         tmonad (]`]`(_."_)`(_."_)`(trmmlucu chk3sm)))@>"0 argsam
  ('trsmrlnn'         tmonad (]`]`(_."_)`(_."_)`(trmmrlnn chk3sm)))@>"0 argsan
  ('trsmrlnu'         tmonad (]`]`(_."_)`(_."_)`(trmmrlnu chk3sm)))@>"0 argsan
  ('trsmrltn'         tmonad (]`]`(_."_)`(_."_)`(trmmrltn chk3sm)))@>"0 argsan
  ('trsmrltu'         tmonad (]`]`(_."_)`(_."_)`(trmmrltu chk3sm)))@>"0 argsan
  ('trsmrljn'         tmonad (]`]`(_."_)`(_."_)`(trmmrljn chk3sm)))@>"0 argsan
  ('trsmrlju'         tmonad (]`]`(_."_)`(_."_)`(trmmrlju chk3sm)))@>"0 argsan
  ('trsmrlcn'         tmonad (]`]`(_."_)`(_."_)`(trmmrlcn chk3sm)))@>"0 argsan
  ('trsmrlcu'         tmonad (]`]`(_."_)`(_."_)`(trmmrlcu chk3sm)))@>"0 argsan
  ('trsmrunn'         tmonad (]`]`(_."_)`(_."_)`(trmmrunn chk3sm)))@>"0 argsan
  ('trsmrunu'         tmonad (]`]`(_."_)`(_."_)`(trmmrunu chk3sm)))@>"0 argsan
  ('trsmrutn'         tmonad (]`]`(_."_)`(_."_)`(trmmrutn chk3sm)))@>"0 argsan
  ('trsmrutu'         tmonad (]`]`(_."_)`(_."_)`(trmmrutu chk3sm)))@>"0 argsan
  ('trsmrujn'         tmonad (]`]`(_."_)`(_."_)`(trmmrujn chk3sm)))@>"0 argsan
  ('trsmruju'         tmonad (]`]`(_."_)`(_."_)`(trmmruju chk3sm)))@>"0 argsan
  ('trsmrucn'         tmonad (]`]`(_."_)`(_."_)`(trmmrucn chk3sm)))@>"0 argsan
  ('trsmrucu'         tmonad (]`]`(_."_)`(_."_)`(trmmrucu chk3sm)))@>"0 argsan

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
  ('trsmllnn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    trlpick ) t02v))) Lm  ; (Lm   mp       bm ) ; bm ; Am ; norm1Lm
  ('trsmllnu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    trl1pick) t02v))) L1m ; (L1m  mp       bm ) ; bm ; Am ; norm1L1m
  ('trsmlltn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@trlpick ) t02v))) Lm  ; (Lm  (mp~ ct)~ bm ) ; bm ; Am ; normiLm
  ('trsmlltu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@trl1pick) t02v))) L1m ; (L1m (mp~ ct)~ bm ) ; bm ; Am ; normiL1m
  ('trsmlljn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ + @trlpick ) t02v))) Lm  ; (Lm  (mp~ + )~ bm ) ; bm ; Am ; normiLm
  ('trsmllju'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ + @trl1pick) t02v))) L1m ; (L1m (mp~ + )~ bm ) ; bm ; Am ; normiL1m
  ('trsmllcn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trlpick ) t02v))) Lm  ; (Lm  (mp~ |:)~ bm ) ; bm ; Am ; normiLm
  ('trsmllcu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trl1pick) t02v))) L1m ; (L1m (mp~ |:)~ bm ) ; bm ; Am ; normiL1m
  ('trsmlunn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    trupick ) t02v))) Um  ; (Um   mp       bm ) ; bm ; Am ; norm1Um
  ('trsmlunu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    tru1pick) t02v))) U1m ; (U1m  mp       bm ) ; bm ; Am ; norm1U1m
  ('trsmlutn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@trupick ) t02v))) Um  ; (Um  (mp~ ct)~ bm ) ; bm ; Am ; normiUm
  ('trsmlutu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@tru1pick) t02v))) U1m ; (U1m (mp~ ct)~ bm ) ; bm ; Am ; normiU1m
  ('trsmlujn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ + @trupick ) t02v))) Um  ; (Um  (mp~ + )~ bm ) ; bm ; Am ; normiUm
  ('trsmluju'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ + @tru1pick) t02v))) U1m ; (U1m (mp~ + )~ bm ) ; bm ; Am ; normiU1m
  ('trsmlucn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trupick ) t02v))) Um  ; (Um  (mp~ |:)~ bm ) ; bm ; Am ; normiUm
  ('trsmlucu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@tru1pick) t02v))) U1m ; (U1m (mp~ |:)~ bm ) ; bm ; Am ; normiU1m
  ('trsmrlnn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trlpick ) t02v))) Ln  ; (bn   mp       Ln ) ; bn ; An ; normiLn
  ('trsmrlnu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trl1pick) t02v))) L1n ; (bn   mp       L1n) ; bn ; An ; normiL1n
  ('trsmrltn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trlpick ) t02v))) Ln  ; (bn  (mp  ct)  Ln ) ; bn ; An ; norm1Ln
  ('trsmrltu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trl1pick) t02v))) L1n ; (bn  (mp  ct)  L1n) ; bn ; An ; norm1L1n
  ('trsmrljn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  + @trlpick ) t02v))) Ln  ; (bn  (mp  + )  Ln ) ; bn ; An ; norm1Ln
  ('trsmrlju'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  + @trl1pick) t02v))) L1n ; (bn  (mp  + )  L1n) ; bn ; An ; norm1L1n
  ('trsmrlcn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trlpick ) t02v))) Ln  ; (bn  (mp  |:)  Ln ) ; bn ; An ; norm1Ln
  ('trsmrlcu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trl1pick) t02v))) L1n ; (bn  (mp  |:)  L1n) ; bn ; An ; norm1L1n
  ('trsmrunn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trupick ) t02v))) Un  ; (bn   mp       Un ) ; bn ; An ; normiUn
  ('trsmrunu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     tru1pick) t02v))) U1n ; (bn   mp       U1n) ; bn ; An ; normiU1n
  ('trsmrutn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trupick ) t02v))) Un  ; (bn  (mp  ct)  Un ) ; bn ; An ; norm1Un
  ('trsmrutu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@tru1pick) t02v))) U1n ; (bn  (mp  ct)  U1n) ; bn ; An ; norm1U1n
  ('trsmrujn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  + @trupick ) t02v))) Un  ; (bn  (mp  + )  Un ) ; bn ; An ; norm1Un
  ('trsmruju'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  + @tru1pick) t02v))) U1n ; (bn  (mp  + )  U1n) ; bn ; An ; norm1U1n
  ('trsmrucn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trupick ) t02v))) Un  ; (bn  (mp  |:)  Un ) ; bn ; An ; norm1Un
  ('trsmrucu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@tru1pick) t02v))) U1n ; (bn  (mp  |:)  U1n) ; bn ; An ; norm1U1n

  NB. dyadic trsmxxxx, 2-rank B and X
  NB. notes:
  NB. - we use RHS matrix B here as solution vector X for
  NB.   (2{::x), RHS is computed explicitely and is supplied
  NB.   in (1{::x)
  ('trsmllnn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    trlpick ) t02m norm1tc))) Lm  ; (Lm   mp       B ) ; B  ; Am ; norm1Lm
  ('trsmllnu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    trl1pick) t02m norm1tc))) L1m ; (L1m  mp       B ) ; B  ; Am ; norm1L1m
  ('trsmlltn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@trlpick ) t02m norm1tc))) Lm  ; (Lm  (mp~ ct)~ B ) ; B  ; Am ; normiLm
  ('trsmlltu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@trl1pick) t02m norm1tc))) L1m ; (L1m (mp~ ct)~ B ) ; B  ; Am ; normiL1m
  ('trsmllcn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trlpick ) t02m norm1tc))) Lm  ; (Lm  (mp~ |:)~ B ) ; B  ; Am ; normiLm
  ('trsmllcu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trl1pick) t02m norm1tc))) L1m ; (L1m (mp~ |:)~ B ) ; B  ; Am ; normiL1m
  ('trsmlunn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    trupick ) t02m norm1tc))) Um  ; (Um   mp       B ) ; B  ; Am ; norm1Um
  ('trsmlunu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    tru1pick) t02m norm1tc))) U1m ; (U1m  mp       B ) ; B  ; Am ; norm1U1m
  ('trsmlutn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@trupick ) t02m norm1tc))) Um  ; (Um  (mp~ ct)~ B ) ; B  ; Am ; normiUm
  ('trsmlutu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@tru1pick) t02m norm1tc))) U1m ; (U1m (mp~ ct)~ B ) ; B  ; Am ; normiU1m
  ('trsmlucn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trupick ) t02m norm1tc))) Um  ; (Um  (mp~ |:)~ B ) ; B  ; Am ; normiUm
  ('trsmlucu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@tru1pick) t02m norm1tc))) U1m ; (U1m (mp~ |:)~ B ) ; B  ; Am ; normiU1m
  ('trsmrlnn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trlpick ) t02m normitc))) Ln  ; (B   mp       Ln ) ; B  ; An ; normiLn
  ('trsmrlnu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trl1pick) t02m normitc))) L1n ; (B   mp       L1n) ; B  ; An ; normiL1n
  ('trsmrltn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trlpick ) t02m normitc))) Ln  ; (B  (mp  ct)  Ln ) ; B  ; An ; norm1Ln
  ('trsmrltu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trl1pick) t02m normitc))) L1n ; (B  (mp  ct)  L1n) ; B  ; An ; norm1L1n
  ('trsmrlcn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trlpick ) t02m normitc))) Ln  ; (B  (mp  |:)  Ln ) ; B  ; An ; norm1Ln
  ('trsmrlcu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trl1pick) t02m normitc))) L1n ; (B  (mp  |:)  L1n) ; B  ; An ; norm1L1n
  ('trsmrunn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trupick ) t02m normitc))) Un  ; (B   mp       Un ) ; B  ; An ; normiUn
  ('trsmrunu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     tru1pick) t02m normitc))) U1n ; (B   mp       U1n) ; B  ; An ; normiU1n
  ('trsmrutn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trupick ) t02m normitc))) Un  ; (B  (mp  ct)  Un ) ; B  ; An ; norm1Un
  ('trsmrutu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@tru1pick) t02m normitc))) U1n ; (B  (mp  ct)  U1n) ; B  ; An ; norm1U1n
  ('trsmrucn0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trupick ) t02m normitc))) Un  ; (B  (mp  |:)  Un ) ; B  ; An ; norm1Un
  ('trsmrucu0'        tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@tru1pick) t02m normitc))) U1n ; (B  (mp  |:)  U1n) ; B  ; An ; norm1U1n

  ('trsmllnn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    trlpick ) t02m norm1tc))) Lm  ; (Lm   mp       B ) ; B  ; Am ; norm1Lm
  ('trsmllnu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    trl1pick) t02m norm1tc))) L1m ; (L1m  mp       B ) ; B  ; Am ; norm1L1m
  ('trsmlltn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@trlpick ) t02m norm1tc))) Lm  ; (Lm  (mp~ ct)~ B ) ; B  ; Am ; normiLm
  ('trsmlltu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@trl1pick) t02m norm1tc))) L1m ; (L1m (mp~ ct)~ B ) ; B  ; Am ; normiL1m
  ('trsmlljn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ + @trlpick ) t02m norm1tc))) Lm  ; (Lm  (mp~ + )~ B ) ; B  ; Am ; normiLm
  ('trsmllju'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ + @trl1pick) t02m norm1tc))) L1m ; (L1m (mp~ + )~ B ) ; B  ; Am ; normiL1m
  ('trsmllcn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trlpick ) t02m norm1tc))) Lm  ; (Lm  (mp~ |:)~ B ) ; B  ; Am ; normiLm
  ('trsmllcu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trl1pick) t02m norm1tc))) L1m ; (L1m (mp~ |:)~ B ) ; B  ; Am ; normiL1m
  ('trsmlunn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    trupick ) t02m norm1tc))) Um  ; (Um   mp       B ) ; B  ; Am ; norm1Um
  ('trsmlunu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    tru1pick) t02m norm1tc))) U1m ; (U1m  mp       B ) ; B  ; Am ; norm1U1m
  ('trsmlutn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@trupick ) t02m norm1tc))) Um  ; (Um  (mp~ ct)~ B ) ; B  ; Am ; normiUm
  ('trsmlutu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@tru1pick) t02m norm1tc))) U1m ; (U1m (mp~ ct)~ B ) ; B  ; Am ; normiU1m
  ('trsmlujn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ + @trupick ) t02m norm1tc))) Um  ; (Um  (mp~ + )~ B ) ; B  ; Am ; normiUm
  ('trsmluju'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ + @tru1pick) t02m norm1tc))) U1m ; (U1m (mp~ + )~ B ) ; B  ; Am ; normiU1m
  ('trsmlucn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trupick ) t02m norm1tc))) Um  ; (Um  (mp~ |:)~ B ) ; B  ; Am ; normiUm
  ('trsmlucu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@tru1pick) t02m norm1tc))) U1m ; (U1m (mp~ |:)~ B ) ; B  ; Am ; normiU1m
  ('trsmrlnn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trlpick ) t02m normitc))) Ln  ; (B   mp       Ln ) ; B  ; An ; normiLn
  ('trsmrlnu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trl1pick) t02m normitc))) L1n ; (B   mp       L1n) ; B  ; An ; normiL1n
  ('trsmrltn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trlpick ) t02m normitc))) Ln  ; (B  (mp  ct)  Ln ) ; B  ; An ; norm1Ln
  ('trsmrltu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trl1pick) t02m normitc))) L1n ; (B  (mp  ct)  L1n) ; B  ; An ; norm1L1n
  ('trsmrljn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  + @trlpick ) t02m normitc))) Ln  ; (B  (mp  + )  Ln ) ; B  ; An ; norm1Ln
  ('trsmrlju'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  + @trl1pick) t02m normitc))) L1n ; (B  (mp  + )  L1n) ; B  ; An ; norm1L1n
  ('trsmrlcn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trlpick ) t02m normitc))) Ln  ; (B  (mp  |:)  Ln ) ; B  ; An ; norm1Ln
  ('trsmrlcu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trl1pick) t02m normitc))) L1n ; (B  (mp  |:)  L1n) ; B  ; An ; norm1L1n
  ('trsmrunn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trupick ) t02m normitc))) Un  ; (B   mp       Un ) ; B  ; An ; normiUn
  ('trsmrunu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     tru1pick) t02m normitc))) U1n ; (B   mp       U1n) ; B  ; An ; normiU1n
  ('trsmrutn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trupick ) t02m normitc))) Un  ; (B  (mp  ct)  Un ) ; B  ; An ; norm1Un
  ('trsmrutu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@tru1pick) t02m normitc))) U1n ; (B  (mp  ct)  U1n) ; B  ; An ; norm1U1n
  ('trsmrujn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  + @trupick ) t02m normitc))) Un  ; (B  (mp  + )  Un ) ; B  ; An ; norm1Un
  ('trsmruju'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  + @tru1pick) t02m normitc))) U1n ; (B  (mp  + )  U1n) ; B  ; An ; norm1U1n
  ('trsmrucn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trupick ) t02m normitc))) Un  ; (B  (mp  |:)  Un ) ; B  ; An ; norm1Un
  ('trsmrucu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@tru1pick) t02m normitc))) U1n ; (B  (mp  |:)  U1n) ; B  ; An ; norm1U1n

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbasicsm
NB.
NB. Description:
NB.   Adv. to make verb to test basic matrix equation solvers
NB.
NB. Syntax:
NB.   vtest=. mkmat testbasicsm
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms; is called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testbasicsm_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testbasicsm_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbasicsm_mt_ 150 200

testbasicsm=: 1 : 'EMPTY [ testbasictrsm_mt_@(u@(2 # >./) ; u) [ load@''math/mt/test/blis/sm'' [ load@''math/mt/test/blas/sm'''

NB. ---------------------------------------------------------
NB. testbasic
NB.
NB. Description:
NB.   Adv. to make verb to test basic operations all levels
NB.
NB. Syntax:
NB.   vtest=. mkmat testbasic
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms; is called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testbasic_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testbasic_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbasic_mt_ 150 200

testbasic=: 1 : 'EMPTY [ u testbasicsm_mt_ [ (u testbasicsv_mt_)^:(=/) [ u testbasicmm_mt_ [ u testbasicmv_mt_ [ u testbasicr2k_mt_ [ u testbasicrk_mt_ [ (u testbasicr2_mt_)^:(=/) [ u testbasicr_mt_'
