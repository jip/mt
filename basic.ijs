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
NB. 1) test suite here is aimed to benchmark BLAS subroutines
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
NB.   AAupd=. (cmp her mul) alpha ; x ; incx ; AA
NB. where
NB.   cmp - dyad to define which triangular part of A is to
NB.         be referenced, is one of:
NB.           >:     NB. LT
NB.           <:     NB. UT
NB.   mul - dyad to define the form of op(A), is one of:
NB.           *      NB. the symmetric operation: op(A) := A^T
NB.           (* +)  NB. the hermitian operation: op(A) := A^H

her=: 2 : 0
  'alpha y incy AA'=. y
  y=. incy extract_mt_ y
  AA=. AA + u/~&i.@#`(0&,:)} alpha (] */ v) y
)

NB. ---------------------------------------------------------
NB. Monad    A            R/W in A    op(x)
NB. syrl     symmetric    LT          x^T
NB. syru     symmetric    UT          x^T
NB. herl     Hermitian    LT          x^H
NB. heru     Hermitian    UT          x^H
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

syrl=: >: her  *
syru=: <: her  *

herl=: >: her (* +)
heru=: <: her (* +)

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
NB.   AAupd=. (cmp her2 trans) alpha ; x ; incx ; y ; incy ; AA
NB. where
NB.   cmp   - dyad to define which triangular part of A is to
NB.           be referenced, is one of:
NB.             >:      NB. LT
NB.             <:      NB. UT
NB.   trans - monad to define the form of op1(v) and op2(s),
NB.           is one of:
NB.             |:      NB. the symmetric operation: op1(v) = v^T, op2(s) = s
NB.             ct      NB. the hermitian operation: op1(v) = v^H, op2(s) = conj(s)

her2=: 2 : 0
  'alpha xx incx y incy AA'=. y
  xx=. incx extract_mt_ xx
  y=.  incy extract_mt_ y
  AA=. AA + u/~&i.@#`(0&,:)} (+ v)~ xx */ alpha * v y
)

NB. ---------------------------------------------------------
NB. Monad    A            R/W in A    op1(v)    op2(alpha)
NB. syr2l    symmetric    LT          v^T       alpha
NB. syr2u    symmetric    UT          v^T       alpha
NB. her2l    Hermitian    LT          v^H       conj(alpha)
NB. her2u    Hermitian    UT          v^H       conj(alpha)
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

syr2l=: >: her2 |:
syr2u=: <: her2 |:

her2l=: >: her2 ct_mt_
her2u=: <: her2 ct_mt_

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
NB.   CCupd=. (cmp herk mul) alpha ; A ; beta ; CC
NB. where
NB.   cmp - dyad to define which triangular part of C is to
NB.         be referenced, is one of:
NB.           >:  NB. LT
NB.           <:  NB. UT
NB.   mul - dyad to compute the product either (A * op(A)) or
NB.         (op(A) * A), is called as:
NB.           product=. mul A

herk=: 2 : 0
  'alpha A beta CC'=. y
  CC=. u/~&i.@c_mt_`]} CC ,: (alpha * v~ A) + beta * CC
)

NB. ---------------------------------------------------------
NB. Monad     C            alpha,beta    R/W in C    op1(A)    op2(A)
NB. syrkln    symmetric    any           LT          A         A^T
NB. syrklt    symmetric    any           LT          A^T       A
NB. syrkun    symmetric    any           UT          A         A^T
NB. syrkut    symmetric    any           UT          A^T       A
NB. herkln    Hermitian    real          LT          A         A^H
NB. herklc    Hermitian    real          LT          A^H       A
NB. herkun    Hermitian    real          UT          A         A^H
NB. herkuc    Hermitian    real          UT          A^H       A
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

syrkln=: >: herk (mp_mt_  |:    )
syrklt=: >: herk (mp_mt_~ |:    )
syrkun=: <: herk (mp_mt_  |:    )
syrkut=: <: herk (mp_mt_~ |:    )

herkln=: >: herk (mp_mt_  ct_mt_)
herklc=: >: herk (mp_mt_~ ct_mt_)
herkun=: <: herk (mp_mt_  ct_mt_)
herkuc=: <: herk (mp_mt_~ ct_mt_)

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
NB.   CCupd=. (cmp`trans her2k kind) alpha ; A ; B ; beta ; CC
NB. where
NB.   cmp   - dyad to define which triangular part of C is to
NB.           be referenced, is one of:
NB.             >:      NB. LT
NB.             <:      NB. UT
NB.   trans - monad to define the form of op1(M) and op2(s),
NB.           is one of:
NB.             |:      NB. the symmetric operation: op1(M) = M^T, op2(s) = s
NB.             ct      NB. the hermitian operation: op1(M) = M^H, op2(s) = conj(s)
NB.   kind  - boolean scalar to define operation:
NB.             0       NB. (2)
NB.             1       NB. (1)
NB.
NB. Notes:
NB. - her2k's design solves a problem: how to allow C1 to see
NB.   V2 in the train (V0 C1 V2 A3), the solution is: send V2
NB.   not (V2 A3) into C1, implement A3 functionality inside
NB.   C1 inline, use switch N4 to control A3 behavior, the
NB.   resulting train becomes (V0`V2 C1 N4)

her2k=: 2 : 0
  'alpha A B beta CC'=. y
  CC=. m@.0/~&i.@c_mt_`]} CC ,: ((+ m@.1) alpha * A (mp_mt_~ m@.1)~`(mp_mt_ m@.1)@.n B) + beta * CC
)

NB. ---------------------------------------------------------
NB. Monad      C            beta    R/W in C    op1(M)    op2(M)    op3(alpha)
NB. syr2kln    symmetric    any     LT          M         M^T       alpha
NB. syr2klt    symmetric    any     LT          M^T       M         alpha
NB. syr2kun    symmetric    any     UT          M         M^T       alpha
NB. syr2kut    symmetric    any     UT          M^T       M         alpha
NB. her2kln    Hermitian    real    LT          M         M^H       conj(alpha)
NB. her2klc    Hermitian    real    LT          M^H       M         conj(alpha)
NB. her2kun    Hermitian    real    UT          M         M^H       conj(alpha)
NB. her2kuc    Hermitian    real    UT          M^H       M         conj(alpha)
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

syr2kln=: >:`|:     her2k 1
syr2klt=: >:`|:     her2k 0
syr2kun=: <:`|:     her2k 1
syr2kut=: <:`|:     her2k 0

her2kln=: >:`ct_mt_ her2k 1
her2klc=: >:`ct_mt_ her2k 0
her2kun=: <:`ct_mt_ her2k 1
her2kuc=: <:`ct_mt_ her2k 0

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
NB. Monad    A            Reads in A
NB. symvl    symmetric    LT
NB. symvu    symmetric    UT
NB. hemvl    Hermitian    LT
NB. hemvu    Hermitian    UT
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
NB. Monad     op1(A)    op2(B)
NB. gemmnn    A         B
NB. gemmnt    A         B^T
NB. gemmnc    A         B^H
NB. gemmtn    A^T       B
NB. gemmtt    A^T       B^T
NB. gemmtc    A^T       B^H
NB. gemmcn    A^H       B
NB. gemmct    A^H       B^T
NB. gemmcc    A^H       B^H
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
NB.   ma    = m for gemmnx or ma = k otherwise
NB.   ka    = k for gemmnx or ka = m otherwise
NB.   kb    = k for gemmxn or kb = n otherwise
NB.   nb    = n for gemmxn or nb = k otherwise
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
gemmnc=:  (mp_mt_  ct_mt_)          gemm
gemmtn=:  (mp_mt_~ |:    )~         gemm
gemmtt=:   mp_mt_& |:               gemm
gemmtc=: ((mp_mt_~ |:    )~ ct_mt_) gemm
gemmcn=:  (mp_mt_~ ct_mt_)~         gemm
gemmct=: ((mp_mt_~ ct_mt_)~ |:    ) gemm
gemmcc=:   mp_mt_& ct_mt_           gemm

NB. ---------------------------------------------------------
NB. hemm
NB.
NB. Description:
NB.   Adv. to make monad to perform the matrix-matrix
NB.   operation:
NB.     C := alpha * A * B + beta * C
NB.   or
NB.     C := alpha * B * A + beta * C
NB.   where A is Hermitian (symmetric)
NB.
NB. Syntax:
NB.   Cupd=. (mul hemm) alpha ; AA ; B ; beta ; C
NB. where
NB.   mul - dyad to read triangular part of A and to compute
NB.         the product either (A * B) or (B * A), is called
NB.         as:
NB.           product=. AA mul B

hemm=: gemm

NB. ---------------------------------------------------------
NB. Monad     A            Side    Reads in A
NB. symmll    symmetric    (1)     LT
NB. symmlu    symmetric    (1)     UT
NB. symmrl    symmetric    (2)     LT
NB. symmru    symmetric    (2)     UT
NB. hemmll    Hermitian    (1)     LT
NB. hemmlu    Hermitian    (1)     UT
NB. hemmrl    Hermitian    (2)     LT
NB. hemmru    Hermitian    (2)     UT
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     C := alpha * A * B + beta * C  (1)
NB.   or
NB.     C := alpha * B * A + beta * C  (2)
NB.   where A is Hermitian (symmetric)
NB.
NB. Syntax:
NB.   Cupd=. xxmmxx alpha ; AA ; B ; beta ; C
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
NB.   mn    = m for xxmmlx or mn = n for xxmmrx
NB.
NB. Notes:
NB. - monad     models BLAS'
NB.   symmll    xSYMM('L','L',...)
NB.   symmlu    xSYMM('L','U',...)
NB.   symmrl    xSYMM('R','L',...)
NB.   symmru    xSYMM('R','U',...)
NB.   hemmll    ZHEMM('L','L',...)
NB.   hemmlu    ZHEMM('L','U',...)
NB.   hemmrl    ZHEMM('R','L',...)
NB.   hemmru    ZHEMM('R','U',...)
NB. - reference implementation

symmll=: (mp_mt_~ sy4gel_mt_)~ hemm
symmlu=: (mp_mt_~ sy4geu_mt_)~ hemm
symmrl=: (mp_mt_  sy4gel_mt_)~ hemm
symmru=: (mp_mt_  sy4geu_mt_)~ hemm

hemmll=: (mp_mt_~ he4gel_mt_)~ hemm
hemmlu=: (mp_mt_~ he4geu_mt_)~ hemm
hemmrl=: (mp_mt_  he4gel_mt_)~ hemm
hemmru=: (mp_mt_  he4geu_mt_)~ hemm

NB. ---------------------------------------------------------
NB. trmm
NB.
NB. Description:
NB.   Adv. to make monad to perform the matrix-matrix
NB.   operation:
NB.     B := alpha * op(A) * B  (1)
NB.   or
NB.     B := alpha * B * op(A)  (2)
NB.   where A is triangular, and op(A) is either A, A^T or
NB.   A^H
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
NB.             ]         NB. op(A) := A
NB.             |:        NB. op(A) := A^T
NB.             ct        NB. op(A) := A^H
NB.   mul   - dyad to compute the product either (op(A) * B) or
NB.           (B * op(A)), is called as:
NB.             product=. B mul opA
NB.           and is one of:
NB.             mp~       NB. to compute (1)
NB.             mp        NB. to compute (2)

trmm=: 1 : '(0&{:: * 2&{::) u 1&{::'

NB. ---------------------------------------------------------
NB. Monad       Side    A     Reads in A    op(A)
NB. trmmllnn    (1)     L      LT           A
NB. trmmllnu    (1)     L1    SLT           A
NB. trmmlltn    (1)     L      LT           A^T
NB. trmmlltu    (1)     L1    SLT           A^T
NB. trmmllcn    (1)     L      LT           A^H
NB. trmmllcu    (1)     L1    SLT           A^H
NB. trmmlunn    (1)     U      UT           A
NB. trmmlunu    (1)     U1    SUT           A
NB. trmmlutn    (1)     U      UT           A^T
NB. trmmlutu    (1)     U1    SUT           A^T
NB. trmmlucn    (1)     U      UT           A^H
NB. trmmlucu    (1)     U1    SUT           A^H
NB. trmmrlnn    (2)     L      LT           A
NB. trmmrlnu    (2)     L1    SLT           A
NB. trmmrltn    (2)     L      LT           A^T
NB. trmmrltu    (2)     L1    SLT           A^T
NB. trmmrlcn    (2)     L      LT           A^H
NB. trmmrlcu    (2)     L1    SLT           A^H
NB. trmmrunn    (2)     U      UT           A
NB. trmmrunu    (2)     U1    SUT           A
NB. trmmrutn    (2)     U      UT           A^T
NB. trmmrutu    (2)     U1    SUT           A^T
NB. trmmrucn    (2)     U      UT           A^H
NB. trmmrucu    (2)     U1    SUT           A^H
NB.
NB. Description:
NB.   Performs the matrix-matrix operation:
NB.     B := alpha * op(A) * B  (1)
NB.   or
NB.     B := alpha * B * op(A)  (2)
NB.   where A is triangular
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
trmmllcn=: (mp_mt_~ ct_mt_@trlpick_mt_ ) trmm
trmmllcu=: (mp_mt_~ ct_mt_@trl1pick_mt_) trmm
trmmlunn=: (mp_mt_~        trupick_mt_ ) trmm
trmmlunu=: (mp_mt_~        tru1pick_mt_) trmm
trmmlutn=: (mp_mt_~ |:    @trupick_mt_ ) trmm
trmmlutu=: (mp_mt_~ |:    @tru1pick_mt_) trmm
trmmlucn=: (mp_mt_~ ct_mt_@trupick_mt_ ) trmm
trmmlucu=: (mp_mt_~ ct_mt_@tru1pick_mt_) trmm
trmmrlnn=: (mp_mt_         trlpick_mt_ ) trmm
trmmrlnu=: (mp_mt_         trl1pick_mt_) trmm
trmmrltn=: (mp_mt_  |:    @trlpick_mt_ ) trmm
trmmrltu=: (mp_mt_  |:    @trl1pick_mt_) trmm
trmmrlcn=: (mp_mt_  ct_mt_@trlpick_mt_ ) trmm
trmmrlcu=: (mp_mt_  ct_mt_@trl1pick_mt_) trmm
trmmrunn=: (mp_mt_         trupick_mt_ ) trmm
trmmrunu=: (mp_mt_         tru1pick_mt_) trmm
trmmrutn=: (mp_mt_  |:    @trupick_mt_ ) trmm
trmmrutu=: (mp_mt_  |:    @tru1pick_mt_) trmm
trmmrucn=: (mp_mt_  ct_mt_@trupick_mt_ ) trmm
trmmrucu=: (mp_mt_  ct_mt_@tru1pick_mt_) trmm

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
NB. trsmllnn    (1)     L      LT           A
NB. trsmllnu    (1)     L1    SLT           A
NB. trsmlltn    (1)     L      LT           A^T
NB. trsmlltu    (1)     L1    SLT           A^T
NB. trsmllcn    (1)     L      LT           A^H
NB. trsmllcu    (1)     L1    SLT           A^H
NB. trsmlunn    (1)     U      UT           A
NB. trsmlunu    (1)     U1    SUT           A
NB. trsmlutn    (1)     U      UT           A^T
NB. trsmlutu    (1)     U1    SUT           A^T
NB. trsmlucn    (1)     U      UT           A^H
NB. trsmlucu    (1)     U1    SUT           A^H
NB. trsmrlnn    (2)     L      LT           A
NB. trsmrlnu    (2)     L1    SLT           A
NB. trsmrltn    (2)     L      LT           A^T
NB. trsmrltu    (2)     L1    SLT           A^T
NB. trsmrlcn    (2)     L      LT           A^H
NB. trsmrlcu    (2)     L1    SLT           A^H
NB. trsmrunn    (2)     U      UT           A
NB. trsmrunu    (2)     U1    SUT           A
NB. trsmrutn    (2)     U      UT           A^T
NB. trsmrutu    (2)     U1    SUT           A^T
NB. trsmrucn    (2)     U      UT           A^H
NB. trsmrucu    (2)     U1    SUT           A^H
NB.
NB. Description:
NB.   Ambivalent verb to solve the linear monomial matrix
NB.   equation:
NB.     op(A) * X = alpha * B  (1)
NB.   or
NB.     X * op(A) = alpha * B  (2)
NB.   where A is triangular
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
trsmllcn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i. - # y do. z=. z ,~     ((i {   y) - ((>: i)      }. ai    ) mp  z) % i { ai=. i {"1 x end. + z}}

trsmlunn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i. - # y do. z=. z ,~     ((i {   y) - ((>: i)      }. ai    ) mp  z) % i { ai=. i {   x end.   z}}
trsmlutn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i.   # y do. z=. z ,      ((i {   y) - (    i       {. ai    ) mp  z) % i { ai=. i {"1 x end.   z}}
trsmlucn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i.   # y do. z=. z ,      ((i {   y) - (    i       {. ai    ) mp  z) % i { ai=. i {"1 x end. + z}}

trsmrlnn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i. - c y do. z=. z ,~"1 0 ((i {"1 y) - ((>: i)      }. ai    ) mp~ z) % i { ai=. i {"1 x end.   z}}
trsmrltn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i.   c y do. z=. z , "1 0 ((i {"1 y) - (    i       {. ai    ) mp~ z) % i { ai=. i {   x end.   z}}
trsmrlcn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1 y=. + y for_i. i.   c y do. z=. z , "1 0 ((i {"1 y) - (    i       {. ai    ) mp~ z) % i { ai=. i {   x end. + z}}

trsmrunn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i.   c y do. z=. z , "1 0 ((i {"1 y) - (    i       {. ai    ) mp~ z) % i { ai=. i {"1 x end.   z}}
trsmrutn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i. - c y do. z=. z ,~"1 0 ((i {"1 y) - ((>: i)      }. ai    ) mp~ z) % i { ai=. i {   x end.   z}}
trsmrucn=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1 y=. + y for_i. i. - c y do. z=. z ,~"1 0 ((i {"1 y) - ((>: i)      }. ai    ) mp~ z) % i { ai=. i {   x end. + z}}

trsmllnu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i.   # y do. z=. z ,       (i {   y) - (    i (   [ {. {  ) x) mp  z                     end.   z}}
trsmlltu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i. - # y do. z=. z ,~      (i {   y) - (    i (>:@[ }. {"1) x) mp  z                     end.   z}}
trsmllcu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i. - # y do. z=. z ,~      (i {   y) - (    i (>:@[ }. {"1) x) mp  z                     end. + z}}

trsmlunu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i. - # y do. z=. z ,~      (i {   y) - (    i (>:@[ }. {  ) x) mp  z                     end.   z}}
trsmlutu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.         y for_i. i.   # y do. z=. z ,       (i {   y) - (    i (   [ {. {"1) x) mp  z                     end.   z}}
trsmlucu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {.   y=. + y for_i. i.   # y do. z=. z ,       (i {   y) - (    i (   [ {. {"1) x) mp  z                     end. + z}}

trsmrlnu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i. - c y do. z=. z ,~"1 0  (i {"1 y) - (    i (>:@[ }. {"1) x) mp~ z                     end.   z}}
trsmrltu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i.   c y do. z=. z , "1 0  (i {"1 y) - (    i (   [ {. {  ) x) mp~ z                     end.   z}}
trsmrlcu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1 y=. + y for_i. i.   c y do. z=. z , "1 0  (i {"1 y) - (    i (   [ {. {  ) x) mp~ z                     end. + z}}

trsmrunu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i.   c y do. z=. z , "1 0  (i {"1 y) - (    i (   [ {. {"1) x) mp~ z                     end.   z}}
trsmrutu=: (1&{:: $: (0&{:: * 2&{::)) : {{z=. 0 {."1       y for_i. i. - c y do. z=. z ,~"1 0  (i {"1 y) - (    i (>:@[ }. {  ) x) mp~ z                     end.   z}}
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

  NB. for every i feed the tuple (alpha_i ; expanded_x_i ; incx_i ; expanded_y_i ; incy_i ; A) to tmonad
  ('dger_mtbla_'  tmonad (]`]`(_."_)`(_."_)`(geru chk4r)))@(3 expand 4)@(1 expand 2)@>"0 { (<"0 dcoeff) ; (< x) ; (<"0 inc) ; (< y) ; (<"0 inc) ; < < A
  ('zgerc_mtbla_' tmonad (]`]`(_."_)`(_."_)`(gerc chk4r)))@(3 expand 4)@(1 expand 2)@>"0 { (<"0 zcoeff) ; (< x) ; (<"0 inc) ; (< y) ; (<"0 inc) ; < < A
  ('zgeru_mtbla_' tmonad (]`]`(_."_)`(_."_)`(geru chk4r)))@(3 expand 4)@(1 expand 2)@>"0 { (<"0 zcoeff) ; (< x) ; (<"0 inc) ; (< y) ; (<"0 inc) ; < < A

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

  NB. for every i feed the tuple (alpha_i ; expanded_x_i ; incx_i ; AA) to tmonad
  ('dsyrl_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrl chk5r trlpick)))@(1 expand 2)@>"0 { ((<"0) dcoeff) ; (< x) ; ((<"0) inc) ; < < AA
  ('dsyru_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syru chk5r trupick)))@(1 expand 2)@>"0 { ((<"0) dcoeff) ; (< x) ; ((<"0) inc) ; < < AA
  ('zherl_mtbla_' tmonad (]`]`(_."_)`(_."_)`(herl chk5r trlpick)))@(1 expand 2)@>"0 { ((<"0) dcoeff) ; (< x) ; ((<"0) inc) ; < < AA
  ('zheru_mtbla_' tmonad (]`]`(_."_)`(_."_)`(heru chk5r trupick)))@(1 expand 2)@>"0 { ((<"0) dcoeff) ; (< x) ; ((<"0) inc) ; < < AA

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

  NB. for every i feed the tuple (alpha_i ; expanded_x_i ; incx_i ; expanded_y_i ; incy_i ; AA) to tmonad
  ('dsyr2l_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2l chk6r2 trlpick)))@(3 expand 4)@(1 expand 2)@>"0 { ((<"0) dcoeff) ; (< x) ; ((<"0) inc) ; (< y) ; ((<"0) inc) ; < < AA
  ('dsyr2u_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2u chk6r2 trupick)))@(3 expand 4)@(1 expand 2)@>"0 { ((<"0) dcoeff) ; (< x) ; ((<"0) inc) ; (< y) ; ((<"0) inc) ; < < AA
  ('zher2l_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2l chk6r2 trlpick)))@(3 expand 4)@(1 expand 2)@>"0 { ((<"0) zcoeff) ; (< x) ; ((<"0) inc) ; (< y) ; ((<"0) inc) ; < < AA
  ('zher2u_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2u chk6r2 trupick)))@(3 expand 4)@(1 expand 2)@>"0 { ((<"0) zcoeff) ; (< x) ; ((<"0) inc) ; (< y) ; ((<"0) inc) ; < < AA

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

  NB. for every i feed the tuple (alpha_i ; A ; beta_i ; CC) to tmonad
  ('dsyrkln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkln chk4rk trlpick)))@>"0 { ((<"0) dcoeff) ; (< A) ; ((<"0) dcoeff) ; < < CCmn [^:(m < n) CC
  ('dsyrklt_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrklt chk4rk trlpick)))@>"0 { ((<"0) dcoeff) ; (< A) ; ((<"0) dcoeff) ; < < CCmn [^:(m > n) CC
  ('dsyrkun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkun chk4rk trupick)))@>"0 { ((<"0) dcoeff) ; (< A) ; ((<"0) dcoeff) ; < < CCmn [^:(m < n) CC
  ('dsyrkut_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkut chk4rk trupick)))@>"0 { ((<"0) dcoeff) ; (< A) ; ((<"0) dcoeff) ; < < CCmn [^:(m > n) CC
  ('zsyrkln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkln chk4rk trlpick)))@>"0 { ((<"0) zcoeff) ; (< A) ; ((<"0) zcoeff) ; < < CCmn [^:(m < n) CC
  ('zsyrklt_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrklt chk4rk trlpick)))@>"0 { ((<"0) zcoeff) ; (< A) ; ((<"0) zcoeff) ; < < CCmn [^:(m > n) CC
  ('zsyrkun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkun chk4rk trupick)))@>"0 { ((<"0) zcoeff) ; (< A) ; ((<"0) zcoeff) ; < < CCmn [^:(m < n) CC
  ('zsyrkut_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syrkut chk4rk trupick)))@>"0 { ((<"0) zcoeff) ; (< A) ; ((<"0) zcoeff) ; < < CCmn [^:(m > n) CC

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

  NB. for every i feed the tuple (alpha_i ; A ; beta_i ; CC) to tmonad
  ('zherkln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(herkln chk4rk trlpick)))@>"0 { ((<"0) dcoeff) ; (< A) ; ((<"0) dcoeff) ; < < CCmn [^:(m < n) CC
  ('zherklc_mtbla_' tmonad (]`]`(_."_)`(_."_)`(herklc chk4rk trlpick)))@>"0 { ((<"0) dcoeff) ; (< A) ; ((<"0) dcoeff) ; < < CCmn [^:(m > n) CC
  ('zherkun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(herkun chk4rk trupick)))@>"0 { ((<"0) dcoeff) ; (< A) ; ((<"0) dcoeff) ; < < CCmn [^:(m < n) CC
  ('zherkuc_mtbla_' tmonad (]`]`(_."_)`(_."_)`(herkuc chk4rk trupick)))@>"0 { ((<"0) dcoeff) ; (< A) ; ((<"0) dcoeff) ; < < CCmn [^:(m > n) CC

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

  NB. for every i feed the tuple (alpha_i ; A ; B ; beta_i ; CC) to tmonad
  ('dsyr2kln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kln chk5r2k trlpick)))@>"0 { ((<"0) dcoeff) ; (< A) ; (< B) ; ((<"0) dcoeff) ; < < CCmn [^:(m < n) CC
  ('dsyr2klt_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2klt chk5r2k trlpick)))@>"0 { ((<"0) dcoeff) ; (< A) ; (< B) ; ((<"0) dcoeff) ; < < CCmn [^:(m > n) CC
  ('dsyr2kun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kun chk5r2k trupick)))@>"0 { ((<"0) dcoeff) ; (< A) ; (< B) ; ((<"0) dcoeff) ; < < CCmn [^:(m < n) CC
  ('dsyr2kut_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kut chk5r2k trupick)))@>"0 { ((<"0) dcoeff) ; (< A) ; (< B) ; ((<"0) dcoeff) ; < < CCmn [^:(m > n) CC
  ('zsyr2kln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kln chk5r2k trlpick)))@>"0 { ((<"0) zcoeff) ; (< A) ; (< B) ; ((<"0) zcoeff) ; < < CCmn [^:(m < n) CC
  ('zsyr2klt_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2klt chk5r2k trlpick)))@>"0 { ((<"0) zcoeff) ; (< A) ; (< B) ; ((<"0) zcoeff) ; < < CCmn [^:(m > n) CC
  ('zsyr2kun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kun chk5r2k trupick)))@>"0 { ((<"0) zcoeff) ; (< A) ; (< B) ; ((<"0) zcoeff) ; < < CCmn [^:(m < n) CC
  ('zsyr2kut_mtbla_' tmonad (]`]`(_."_)`(_."_)`(syr2kut chk5r2k trupick)))@>"0 { ((<"0) zcoeff) ; (< A) ; (< B) ; ((<"0) zcoeff) ; < < CCmn [^:(m > n) CC

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

  NB. for every i feed the tuple (alpha_i ; A ; B ; beta_i ; CC) to tmonad
  ('zher2kln_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2kln chk5r2k trlpick)))@>"0 { ((<"0) zcoeff) ; (< A) ; (< B) ; ((<"0) dcoeff) ; < < CCmn [^:(m < n) CC
  ('zher2klc_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2klc chk5r2k trlpick)))@>"0 { ((<"0) zcoeff) ; (< A) ; (< B) ; ((<"0) dcoeff) ; < < CCmn [^:(m > n) CC
  ('zher2kun_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2kun chk5r2k trupick)))@>"0 { ((<"0) zcoeff) ; (< A) ; (< B) ; ((<"0) dcoeff) ; < < CCmn [^:(m < n) CC
  ('zher2kuc_mtbla_' tmonad (]`]`(_."_)`(_."_)`(her2kuc chk5r2k trupick)))@>"0 { ((<"0) zcoeff) ; (< A) ; (< B) ; ((<"0) dcoeff) ; < < CCmn [^:(m > n) CC

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

  NB. for every i feed the tuple (alpha_i ; AA ; expanded_x_i ; incx_i ; beta_i ; expanded_y_i ; incy_i) to tmonad
  ('dsymvl_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symvl chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 { ((<"0) dcoeff) ; (< AA) ; (< x) ; ((<"0) inc) ; ((<"0) dcoeff) ; (< y) ; < <"0 inc
  ('dsymvu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symvu chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 { ((<"0) dcoeff) ; (< AA) ; (< x) ; ((<"0) inc) ; ((<"0) dcoeff) ; (< y) ; < <"0 inc
  ('zhemvl_mtbla_' tmonad (]`]`(_."_)`(_."_)`(hemvl chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 { ((<"0) zcoeff) ; (< AA) ; (< x) ; ((<"0) inc) ; ((<"0) zcoeff) ; (< y) ; < <"0 inc
  ('zhemvu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(hemvu chk2mv)))@(5 expand 6)@(2 expand 3)@>"0 { ((<"0) zcoeff) ; (< AA) ; (< x) ; ((<"0) inc) ; ((<"0) zcoeff) ; (< y) ; < <"0 inc

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

  NB. for every i feed the tuple (AA ; expanded_x_i ; incx_i) to tmonad
  ('dtrmvlnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnn chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('dtrmvlnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnu chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('dtrmvltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltn chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('dtrmvltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltu chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('dtrmvunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunn chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('dtrmvunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunu chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('dtrmvutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutn chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('dtrmvutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutu chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvlnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnn chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvlnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlnu chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltn chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvltu chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvlcn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlcn chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvlcu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvlcu chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunn chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvunu chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutn chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvutu chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvucn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvucn chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc
  ('ztrmvucu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmvucu chk3mv)))@(1 expand 2)@>"0 { (< AA) ; (< y) ; < <"0 inc

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
NB.   by general matrices
NB.
NB. Syntax:
NB.   testbasicgemv C ; As ; Bs
NB. where
NB.   C  - m×n-matrix
NB.   As - m×(m+n)-matrix, A material
NB.   Bs - (m+n)×n-matrix, B material

testbasicgemm=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  'C As B'=. y
  'm n'=. $ C
  ks=. /:~ ~. m (0 1 , (, >.@-:)@(, , +)) n  NB. 0,1,⌈m/2⌉,⌈n/2⌉,⌈(m+n)/2⌉,m,n,m+n
  As=. ks <@:({."0 1)"0 _ As                     NB. As[i] is m×k[i]-matrix

  NB. test for the case: ('alpha beta'=. 1.0 0.0) and (op(A) = A)
  ('(+/ .*)'         tdyad  ((0&{::)`(1&{::)`0:`(_."_)`(_."_)`0:                           ))@(c (0 shrink 1)  {.   )@>"0 {                        As  ;  < <  Bs
  ('mp'              tdyad  ((0&{::)`(1&{::)`0:`(_."_)`(_."_)`0:                           ))@(c (0 shrink 1)  {.   )@>"0 {                        As  ;  < <  Bs

  NB. for every i feed the tuple (alpha_i ; A_i ; B_i ; beta_i ; C) to tmonad
  NB. note: A_i and B_i shapes are related; to emulate this,
  NB.       a full fixed B is feeded to Catalogue ({) and
  NB.       then is shrinked to the shape suitable before
  NB.       call to tmonad
  ('dgemmnn_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmnn chk1mm)))@(c (1 shrink 2)  {.   )@>"0 { (<"0 dcoeff) ;         As  ; (<    B) ; (<"0 dcoeff) ; < < C
  ('dgemmnt_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmnt chk1mm)))@(c (1 shrink 2) ({."1))@>"0 { (<"0 dcoeff) ;         As  ; (< |: B) ; (<"0 dcoeff) ; < < C
  ('dgemmtn_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmtn chk1mm)))@(# (1 shrink 2)  {.   )@>"0 { (<"0 dcoeff) ; (|: L:0 As) ; (<    B) ; (<"0 dcoeff) ; < < C
  ('dgemmtt_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmtt chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 dcoeff) ; (|: L:0 As) ; (< |: B) ; (<"0 dcoeff) ; < < C
  ('zgemmnn_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmnn chk1mm)))@(c (1 shrink 2)  {.   )@>"0 { (<"0 zcoeff) ;         As  ; (<    B) ; (<"0 zcoeff) ; < < C
  ('zgemmnt_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmnt chk1mm)))@(c (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ;         As  ; (< |: B) ; (<"0 zcoeff) ; < < C
  ('zgemmnc_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmnc chk1mm)))@(c (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ;         As  ; (< ct B) ; (<"0 zcoeff) ; < < C
  ('zgemmtn_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmtn chk1mm)))@(# (1 shrink 2)  {.   )@>"0 { (<"0 zcoeff) ; (|: L:0 As) ; (<    B) ; (<"0 zcoeff) ; < < C
  ('zgemmtt_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmtt chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (|: L:0 As) ; (< |: B) ; (<"0 zcoeff) ; < < C
  ('zgemmtc_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmtc chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (|: L:0 As) ; (< ct B) ; (<"0 zcoeff) ; < < C
  ('zgemmcn_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmcn chk1mm)))@(# (1 shrink 2)  {.   )@>"0 { (<"0 zcoeff) ; (ct L:0 As) ; (<    B) ; (<"0 zcoeff) ; < < C
  ('zgemmct_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmct chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (ct L:0 As) ; (< |: B) ; (<"0 zcoeff) ; < < C
  ('zgemmcc_mtbla_' tmonad (        ]      `] `(_."_)`(_."_)`(gemmcc chk1mm)))@(# (1 shrink 2) ({."1))@>"0 { (<"0 zcoeff) ; (ct L:0 As) ; (< ct B) ; (<"0 zcoeff) ; < < C

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbasicsymm
NB.
NB. Description:
NB.   Test symmetric matrix-vector operations:
NB.   - xSYMM (BLAS)
NB.   by symmetric matrix
NB.
NB. Syntax:
NB.   testbasicsymm AA ; B ; C
NB. where
NB.   AA - k×k-matrix, A material
NB.   B  - m×n-matrix
NB.   C  - m×n-matrix
NB.   k  = max(m,n)

testbasicsymm=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  'AA B C'=. y
  'm n'=. $ C
  Am=. (2 # m) {. AA
  An=. (2 # n) {. AA

  NB. for every i feed the tuple (alpha_i ; AA ; B ; beta_i ; C) to tmonad
  ('dsymmll_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symmll chk2mm trlpick)))@>"0 { (<"0 dcoeff) ; (< Am) ; (< B) ; (<"0 dcoeff) ; < < C
  ('dsymmlu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symmlu chk2mm trupick)))@>"0 { (<"0 dcoeff) ; (< Am) ; (< B) ; (<"0 dcoeff) ; < < C
  ('dsymmrl_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symmrl chk2mm trlpick)))@>"0 { (<"0 dcoeff) ; (< An) ; (< B) ; (<"0 dcoeff) ; < < C
  ('dsymmru_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symmru chk2mm trupick)))@>"0 { (<"0 dcoeff) ; (< An) ; (< B) ; (<"0 dcoeff) ; < < C
  ('zsymmll_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symmll chk2mm trlpick)))@>"0 { (<"0 zcoeff) ; (< Am) ; (< B) ; (<"0 zcoeff) ; < < C
  ('zsymmlu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symmlu chk2mm trupick)))@>"0 { (<"0 zcoeff) ; (< Am) ; (< B) ; (<"0 zcoeff) ; < < C
  ('zsymmrl_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symmrl chk2mm trlpick)))@>"0 { (<"0 zcoeff) ; (< An) ; (< B) ; (<"0 zcoeff) ; < < C
  ('zsymmru_mtbla_' tmonad (]`]`(_."_)`(_."_)`(symmru chk2mm trupick)))@>"0 { (<"0 zcoeff) ; (< An) ; (< B) ; (<"0 zcoeff) ; < < C

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbasichemm
NB.
NB. Description:
NB.   Test hermitian matrix-vector operation:
NB.   - ZHEMM (BLAS)
NB.   by Hermitian matrix
NB.
NB. Syntax:
NB.   testbasichemm AA ; B ; C
NB. where
NB.   AA - k×k-matrix, A material with real diagonal
NB.   B  - m×n-matrix
NB.   C  - m×n-matrix
NB.   k  = max(m,n)

testbasichemm=: 3 : 0
  zcoeff=. 0j0 1j0 0.7j_0.9
  'AA B C'=. y
  'm n'=. $ C
  Am=. (2 # m) {. AA
  An=. (2 # n) {. AA

  NB. for every i feed the tuple (alpha_i ; AA ; B ; beta_i ; C) to tmonad
  ('zhemmll_mtbla_' tmonad (]`]`(_."_)`(_."_)`(hemmll chk2mm trlpick)))@>"0 { (<"0 zcoeff) ; (< Am) ; (< B) ; (<"0 zcoeff) ; < < C
  ('zhemmlu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(hemmlu chk2mm trupick)))@>"0 { (<"0 zcoeff) ; (< Am) ; (< B) ; (<"0 zcoeff) ; < < C
  ('zhemmrl_mtbla_' tmonad (]`]`(_."_)`(_."_)`(hemmrl chk2mm trlpick)))@>"0 { (<"0 zcoeff) ; (< An) ; (< B) ; (<"0 zcoeff) ; < < C
  ('zhemmru_mtbla_' tmonad (]`]`(_."_)`(_."_)`(hemmru chk2mm trupick)))@>"0 { (<"0 zcoeff) ; (< An) ; (< B) ; (<"0 zcoeff) ; < < C

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbasictrmm
NB.
NB. Description:
NB.   Test matrix-matrix operations:
NB.   - xTRMM (BLAS)
NB.   by triangular matrix
NB.
NB. Syntax:
NB.   testbasictrmm AA ; B
NB. where
NB.   AA - k×k-matrix, A material
NB.   B  - m×n-matrix
NB.   k  = max(m,n)

testbasictrmm=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  'AA B'=. y
  'm n'=. $ B
  Am=. (2 # m) {. AA
  An=. (2 # n) {. AA

  argsdm=. { (<"0 dcoeff) ; (< Am) ; < < B
  argsdn=. { (<"0 dcoeff) ; (< An) ; < < B
  argszm=. { (<"0 zcoeff) ; (< Am) ; < < B
  argszn=. { (<"0 zcoeff) ; (< An) ; < < B

  NB. for every i feed the tuple (alpha_i ; A ; B) to tmonad
  ('dtrmmllnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllnn chk3mm trlpick)))@>"0 argsdm
  ('dtrmmllnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllnu chk3mm trlpick)))@>"0 argsdm
  ('dtrmmlltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlltn chk3mm trlpick)))@>"0 argsdm
  ('dtrmmlltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlltu chk3mm trlpick)))@>"0 argsdm
  ('dtrmmlunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlunn chk3mm trupick)))@>"0 argsdm
  ('dtrmmlunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlunu chk3mm trupick)))@>"0 argsdm
  ('dtrmmlutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlutn chk3mm trupick)))@>"0 argsdm
  ('dtrmmlutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlutu chk3mm trupick)))@>"0 argsdm
  ('dtrmmrlnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnn chk3mm trlpick)))@>"0 argsdn
  ('dtrmmrlnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnu chk3mm trlpick)))@>"0 argsdn
  ('dtrmmrltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrltn chk3mm trlpick)))@>"0 argsdn
  ('dtrmmrltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrltu chk3mm trlpick)))@>"0 argsdn
  ('dtrmmrunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrunn chk3mm trupick)))@>"0 argsdn
  ('dtrmmrunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrunu chk3mm trupick)))@>"0 argsdn
  ('dtrmmrutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrutn chk3mm trupick)))@>"0 argsdn
  ('dtrmmrutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrutu chk3mm trupick)))@>"0 argsdn
  ('ztrmmllnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllnn chk3mm trlpick)))@>"0 argszm
  ('ztrmmllnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllnu chk3mm trlpick)))@>"0 argszm
  ('ztrmmlltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlltn chk3mm trlpick)))@>"0 argszm
  ('ztrmmlltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlltu chk3mm trlpick)))@>"0 argszm
  ('ztrmmllcn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllcn chk3mm trlpick)))@>"0 argszm
  ('ztrmmllcu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmllcu chk3mm trlpick)))@>"0 argszm
  ('ztrmmlunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlunn chk3mm trupick)))@>"0 argszm
  ('ztrmmlunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlunu chk3mm trupick)))@>"0 argszm
  ('ztrmmlutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlutn chk3mm trupick)))@>"0 argszm
  ('ztrmmlutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlutu chk3mm trupick)))@>"0 argszm
  ('ztrmmlucn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlucn chk3mm trupick)))@>"0 argszm
  ('ztrmmlucu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmlucu chk3mm trupick)))@>"0 argszm
  ('ztrmmrlnn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnn chk3mm trlpick)))@>"0 argszn
  ('ztrmmrlnu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlnu chk3mm trlpick)))@>"0 argszn
  ('ztrmmrltn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrltn chk3mm trlpick)))@>"0 argszn
  ('ztrmmrltu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrltu chk3mm trlpick)))@>"0 argszn
  ('ztrmmrlcn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcn chk3mm trlpick)))@>"0 argszn
  ('ztrmmrlcu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrlcu chk3mm trlpick)))@>"0 argszn
  ('ztrmmrunn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrunn chk3mm trupick)))@>"0 argszn
  ('ztrmmrunu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrunu chk3mm trupick)))@>"0 argszn
  ('ztrmmrutn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrutn chk3mm trupick)))@>"0 argszn
  ('ztrmmrutu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrutu chk3mm trupick)))@>"0 argszn
  ('ztrmmrucn_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrucn chk3mm trupick)))@>"0 argszn
  ('ztrmmrucu_mtbla_' tmonad (]`]`(_."_)`(_."_)`(trmmrucu chk3mm trupick)))@>"0 argszn

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

testbasicmm=: 1 : 'EMPTY [ testbasictrmm_mt_@(u@(2 # >./) ; u) [ testbasichemm_mt_@((9&o. upddiag_mt_)@u@(2 # >./) ; u ; u) [ testbasicsymm_mt_@(u@(2 # >./) ; u ; u) [ testbasicgemm_mt_@(u ; u@(+/\) ; u@(+/\.)) [ load@''math/mt/test/blas/mm'''

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
NB.   by triangular matrix
NB.
NB. Syntax:
NB.   testbasictrsm AA ; B
NB. where
NB.   AA - k×k-matrix, A material
NB.   B  - m×n-matrix, RHS
NB.   k  = max(m,n)

testbasictrsm=: 3 : 0
  dcoeff=. 0.0 1.0 0.7
  zcoeff=. 0j0 1j0 0.7j_0.9
  'AA B'=. y
  'm n'=. $ B
  Am=. (2 # m) {. AA
  An=. (2 # n) {. AA
  argsdm=. { ((<"0) dcoeff) ; (< Am) ; < < B
  argsdn=. { ((<"0) dcoeff) ; (< An) ; < < B
  argszm=. { ((<"0) zcoeff) ; (< Am) ; < < B
  argszn=. { ((<"0) zcoeff) ; (< An) ; < < B

  NB. BLAS' staff
  NB. for every i feed the tuple (alpha_i ; AA ; B) to tmonad
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

  acoeff=. /:~ ~. dcoeff , zcoeff
  argsBm=. { (<"0 acoeff) ; (< Am) ; < < B
  argsBn=. { (<"0 acoeff) ; (< An) ; < < B
  'bm bn'=. ({."1 ; {.) B
  argsbm=. { (<"0 acoeff) ; (< Am) ; < < bm
  argsbn=. { (<"0 acoeff) ; (< An) ; < < bn

  NB. monadic trsmxxxx, 1-rank b and x
  ('trsmllnn'         tmonad (]`]`(_."_)`(_."_)`(trmmllnn chk3sm)))@>"0 argsbm
  ('trsmllnu'         tmonad (]`]`(_."_)`(_."_)`(trmmllnu chk3sm)))@>"0 argsbm
  ('trsmlltn'         tmonad (]`]`(_."_)`(_."_)`(trmmlltn chk3sm)))@>"0 argsbm
  ('trsmlltu'         tmonad (]`]`(_."_)`(_."_)`(trmmlltu chk3sm)))@>"0 argsbm
  ('trsmllcn'         tmonad (]`]`(_."_)`(_."_)`(trmmllcn chk3sm)))@>"0 argsbm
  ('trsmllcu'         tmonad (]`]`(_."_)`(_."_)`(trmmllcu chk3sm)))@>"0 argsbm
  ('trsmlunn'         tmonad (]`]`(_."_)`(_."_)`(trmmlunn chk3sm)))@>"0 argsbm
  ('trsmlunu'         tmonad (]`]`(_."_)`(_."_)`(trmmlunu chk3sm)))@>"0 argsbm
  ('trsmlutn'         tmonad (]`]`(_."_)`(_."_)`(trmmlutn chk3sm)))@>"0 argsbm
  ('trsmlutu'         tmonad (]`]`(_."_)`(_."_)`(trmmlutu chk3sm)))@>"0 argsbm
  ('trsmlucn'         tmonad (]`]`(_."_)`(_."_)`(trmmlucn chk3sm)))@>"0 argsbm
  ('trsmlucu'         tmonad (]`]`(_."_)`(_."_)`(trmmlucu chk3sm)))@>"0 argsbm
  ('trsmrlnn'         tmonad (]`]`(_."_)`(_."_)`(trmmrlnn chk3sm)))@>"0 argsbn
  ('trsmrlnu'         tmonad (]`]`(_."_)`(_."_)`(trmmrlnu chk3sm)))@>"0 argsbn
  ('trsmrltn'         tmonad (]`]`(_."_)`(_."_)`(trmmrltn chk3sm)))@>"0 argsbn
  ('trsmrltu'         tmonad (]`]`(_."_)`(_."_)`(trmmrltu chk3sm)))@>"0 argsbn
  ('trsmrlcn'         tmonad (]`]`(_."_)`(_."_)`(trmmrlcn chk3sm)))@>"0 argsbn
  ('trsmrlcu'         tmonad (]`]`(_."_)`(_."_)`(trmmrlcu chk3sm)))@>"0 argsbn
  ('trsmrunn'         tmonad (]`]`(_."_)`(_."_)`(trmmrunn chk3sm)))@>"0 argsbn
  ('trsmrunu'         tmonad (]`]`(_."_)`(_."_)`(trmmrunu chk3sm)))@>"0 argsbn
  ('trsmrutn'         tmonad (]`]`(_."_)`(_."_)`(trmmrutn chk3sm)))@>"0 argsbn
  ('trsmrutu'         tmonad (]`]`(_."_)`(_."_)`(trmmrutu chk3sm)))@>"0 argsbn
  ('trsmrucn'         tmonad (]`]`(_."_)`(_."_)`(trmmrucn chk3sm)))@>"0 argsbn
  ('trsmrucu'         tmonad (]`]`(_."_)`(_."_)`(trmmrucu chk3sm)))@>"0 argsbn

  NB. monadic trsmxxxx, 2-rank B and X
  ('trsmllnn'         tmonad (]`]`(_."_)`(_."_)`(trmmllnn chk3sm)))@>"0 argsBm
  ('trsmllnu'         tmonad (]`]`(_."_)`(_."_)`(trmmllnu chk3sm)))@>"0 argsBm
  ('trsmlltn'         tmonad (]`]`(_."_)`(_."_)`(trmmlltn chk3sm)))@>"0 argsBm
  ('trsmlltu'         tmonad (]`]`(_."_)`(_."_)`(trmmlltu chk3sm)))@>"0 argsBm
  ('trsmllcn'         tmonad (]`]`(_."_)`(_."_)`(trmmllcn chk3sm)))@>"0 argsBm
  ('trsmllcu'         tmonad (]`]`(_."_)`(_."_)`(trmmllcu chk3sm)))@>"0 argsBm
  ('trsmlunn'         tmonad (]`]`(_."_)`(_."_)`(trmmlunn chk3sm)))@>"0 argsBm
  ('trsmlunu'         tmonad (]`]`(_."_)`(_."_)`(trmmlunu chk3sm)))@>"0 argsBm
  ('trsmlutn'         tmonad (]`]`(_."_)`(_."_)`(trmmlutn chk3sm)))@>"0 argsBm
  ('trsmlutu'         tmonad (]`]`(_."_)`(_."_)`(trmmlutu chk3sm)))@>"0 argsBm
  ('trsmlucn'         tmonad (]`]`(_."_)`(_."_)`(trmmlucn chk3sm)))@>"0 argsBm
  ('trsmlucu'         tmonad (]`]`(_."_)`(_."_)`(trmmlucu chk3sm)))@>"0 argsBm
  ('trsmrlnn'         tmonad (]`]`(_."_)`(_."_)`(trmmrlnn chk3sm)))@>"0 argsBn
  ('trsmrlnu'         tmonad (]`]`(_."_)`(_."_)`(trmmrlnu chk3sm)))@>"0 argsBn
  ('trsmrltn'         tmonad (]`]`(_."_)`(_."_)`(trmmrltn chk3sm)))@>"0 argsBn
  ('trsmrltu'         tmonad (]`]`(_."_)`(_."_)`(trmmrltu chk3sm)))@>"0 argsBn
  ('trsmrlcn'         tmonad (]`]`(_."_)`(_."_)`(trmmrlcn chk3sm)))@>"0 argsBn
  ('trsmrlcu'         tmonad (]`]`(_."_)`(_."_)`(trmmrlcu chk3sm)))@>"0 argsBn
  ('trsmrunn'         tmonad (]`]`(_."_)`(_."_)`(trmmrunn chk3sm)))@>"0 argsBn
  ('trsmrunu'         tmonad (]`]`(_."_)`(_."_)`(trmmrunu chk3sm)))@>"0 argsBn
  ('trsmrutn'         tmonad (]`]`(_."_)`(_."_)`(trmmrutn chk3sm)))@>"0 argsBn
  ('trsmrutu'         tmonad (]`]`(_."_)`(_."_)`(trmmrutu chk3sm)))@>"0 argsBn
  ('trsmrucn'         tmonad (]`]`(_."_)`(_."_)`(trmmrucn chk3sm)))@>"0 argsBn
  ('trsmrucu'         tmonad (]`]`(_."_)`(_."_)`(trmmrucu chk3sm)))@>"0 argsBn

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
  ('trsmllcn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trlpick ) t02v))) Lm  ; (Lm  (mp~ |:)~ bm ) ; bm ; Am ; normiLm
  ('trsmllcu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trl1pick) t02v))) L1m ; (L1m (mp~ |:)~ bm ) ; bm ; Am ; normiL1m
  ('trsmlunn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    trupick ) t02v))) Um  ; (Um   mp       bm ) ; bm ; Am ; norm1Um
  ('trsmlunu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    tru1pick) t02v))) U1m ; (U1m  mp       bm ) ; bm ; Am ; norm1U1m
  ('trsmlutn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@trupick ) t02v))) Um  ; (Um  (mp~ ct)~ bm ) ; bm ; Am ; normiUm
  ('trsmlutu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@tru1pick) t02v))) U1m ; (U1m (mp~ ct)~ bm ) ; bm ; Am ; normiU1m
  ('trsmlucn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trupick ) t02v))) Um  ; (Um  (mp~ |:)~ bm ) ; bm ; Am ; normiUm
  ('trsmlucu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@tru1pick) t02v))) U1m ; (U1m (mp~ |:)~ bm ) ; bm ; Am ; normiU1m
  ('trsmrlnn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trlpick ) t02v))) Ln  ; (bn   mp       Ln ) ; bn ; An ; normiLn
  ('trsmrlnu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trl1pick) t02v))) L1n ; (bn   mp       L1n) ; bn ; An ; normiL1n
  ('trsmrltn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trlpick ) t02v))) Ln  ; (bn  (mp  ct)  Ln ) ; bn ; An ; norm1Ln
  ('trsmrltu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trl1pick) t02v))) L1n ; (bn  (mp  ct)  L1n) ; bn ; An ; norm1L1n
  ('trsmrlcn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trlpick ) t02v))) Ln  ; (bn  (mp  |:)  Ln ) ; bn ; An ; norm1Ln
  ('trsmrlcu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trl1pick) t02v))) L1n ; (bn  (mp  |:)  L1n) ; bn ; An ; norm1L1n
  ('trsmrunn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trupick ) t02v))) Un  ; (bn   mp       Un ) ; bn ; An ; normiUn
  ('trsmrunu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     tru1pick) t02v))) U1n ; (bn   mp       U1n) ; bn ; An ; normiU1n
  ('trsmrutn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trupick ) t02v))) Un  ; (bn  (mp  ct)  Un ) ; bn ; An ; norm1Un
  ('trsmrutu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@tru1pick) t02v))) U1n ; (bn  (mp  ct)  U1n) ; bn ; An ; norm1U1n
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
  ('trsmllcn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trlpick ) t02m norm1tc))) Lm  ; (Lm  (mp~ |:)~ B ) ; B  ; Am ; normiLm
  ('trsmllcu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trl1pick) t02m norm1tc))) L1m ; (L1m (mp~ |:)~ B ) ; B  ; Am ; normiL1m
  ('trsmlunn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    trupick ) t02m norm1tc))) Um  ; (Um   mp       B ) ; B  ; Am ; norm1Um
  ('trsmlunu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~    tru1pick) t02m norm1tc))) U1m ; (U1m  mp       B ) ; B  ; Am ; norm1U1m
  ('trsmlutn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@trupick ) t02m norm1tc))) Um  ; (Um  (mp~ ct)~ B ) ; B  ; Am ; normiUm
  ('trsmlutu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ |:@tru1pick) t02m norm1tc))) U1m ; (U1m (mp~ ct)~ B ) ; B  ; Am ; normiU1m
  ('trsmlucn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@trupick ) t02m norm1tc))) Um  ; (Um  (mp~ |:)~ B ) ; B  ; Am ; normiUm
  ('trsmlucu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp~ ct@tru1pick) t02m norm1tc))) U1m ; (U1m (mp~ |:)~ B ) ; B  ; Am ; normiU1m
  ('trsmrlnn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trlpick ) t02m normitc))) Ln  ; (B   mp       Ln ) ; B  ; An ; normiLn
  ('trsmrlnu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trl1pick) t02m normitc))) L1n ; (B   mp       L1n) ; B  ; An ; normiL1n
  ('trsmrltn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trlpick ) t02m normitc))) Ln  ; (B  (mp  ct)  Ln ) ; B  ; An ; norm1Ln
  ('trsmrltu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trl1pick) t02m normitc))) L1n ; (B  (mp  ct)  L1n) ; B  ; An ; norm1L1n
  ('trsmrlcn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trlpick ) t02m normitc))) Ln  ; (B  (mp  |:)  Ln ) ; B  ; An ; norm1Ln
  ('trsmrlcu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  ct@trl1pick) t02m normitc))) L1n ; (B  (mp  |:)  L1n) ; B  ; An ; norm1L1n
  ('trsmrunn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     trupick ) t02m normitc))) Un  ; (B   mp       Un ) ; B  ; An ; normiUn
  ('trsmrunu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp     tru1pick) t02m normitc))) U1n ; (B   mp       U1n) ; B  ; An ; normiU1n
  ('trsmrutn'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@trupick ) t02m normitc))) Un  ; (B  (mp  ct)  Un ) ; B  ; An ; norm1Un
  ('trsmrutu'         tdyad  ((3&{::)`(1&{::)`]`(_."_)`(_."_)`((mp  |:@tru1pick) t02m normitc))) U1n ; (B  (mp  ct)  U1n) ; B  ; An ; norm1U1n
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

testbasicsm=: 1 : 'EMPTY [ testbasictrsm_mt_@(u@(2 # >./) ; u) [ load@''math/mt/test/blas/sm'''

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
