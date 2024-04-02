NB. Structure handlers
NB.
NB. c         Columns in noun
NB. trace     Matrix trace
NB. ct        Conjugate transpose
NB. cp        Conjugate pertranspose
NB. fp        Full (symmetric) permutation
NB. P4p       Transform permutation vector to/from
NB.           permutation matrix
NB. P4ip      Transform inversed permutation vector to/from
NB.           permutation matrix
NB. rt        Restrained Take
NB. icut      Inversed cut
NB.
NB. upd       Adv. to update subarray by a monad
NB. e0        Extend matrix by zeros
NB. appendx   Enhance built-in Append verb (,)
NB. stitchx   Enhance built-in Stitch verb (,.)
NB. ds        Direct sum of matrices A⊕B
NB.
NB. diag      Return a solid part of diagonal
NB. setdiag   Assign value(s) to a solid part of diagonal
NB. upddiag   Adv. to make verbs to update a solid part of
NB.           diagonal
NB.
NB. bdlpick   Zeroize elements outside lower bidiagonal part
NB.           of the matrix
NB. bdupick   Zeroize elements outside upper bidiagonal part
NB.           of the matrix
NB. hslpick   Zeroize elements outside lower Hessenberg part
NB.           of the matrix
NB. hsupick   Zeroize elements outside upper Hessenberg part
NB.           of the matrix
NB. gtpick    Zeroize elements outside tridiagonal part of
NB.           the matrix
NB. trlpick   Zeroize elements outside lower trapezoidal part
NB.           of the matrix
NB. trupick   Zeroize elements outside upper trapezoidal part
NB.           of the matrix
NB. trl1pick  Zeroize elements outside lower trapezoidal part
NB.           of the matrix and set diagonal to 1
NB. tru1pick  Zeroize elements outside upper trapezoidal part
NB.           of the matrix and set diagonal to 1
NB.
NB. idmat     Make identity matrix with units on solid part
NB.           of diagonal
NB. diagmat   Make diagonal matrix
NB. trl       Extract lower trapezoidal matrix
NB. tru       Extract upper trapezoidal matrix
NB. trl0      Extract strictly lower trapezoidal matrix
NB. tru0      Extract strictly upper trapezoidal matrix
NB. trl1      Extract unit lower trapezoidal matrix
NB. tru1      Extract unit upper trapezoidal matrix
NB. xx4gex    Compose structured matrix from SLT (SUT) part
NB.           and diagonal of general square matrix
NB. xxxxxy    Compose matrix from triangular parts of general
NB.           matrices
NB. sxxsxy    Adv. to make dyad to compose matrix from strict
NB.           triangular parts of general matrices
NB. po        Make Hermitian (symmetric) positive definite
NB.           matrix from general square invertible one
NB.
NB. Version: 0.12.0 2021-02-01
NB.
NB. Copyright 2005-2021 Oleg Kobchenko, Roger Hui, Igor Zhuravlov
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

coclass 'mt'

NB. =========================================================
NB. Concepts
NB.
NB. Notation:
NB.   array ranks:
NB.     l-vector   - a vector of length l  e.g. "2-vector v"
NB.     s-matrix   - a matrix of shape  s  e.g. "2×3-matrix M"
NB.     s-brick    - a brick  of shape  s  e.g. "2×3×4-brick B"
NB.     r-rank     - an array of rank   r  e.g. "3-rank array A"
NB.   matrix types:
NB.     BD,BDL,BDU - [{lower,upper}] bidiagonal
NB.     DI         - diagonalizable
NB.     GE         - general
NB.     GG         - general-general pair, generalized form
NB.     GT         - general tridiagonal
NB.     HE         - Hermitian for complex data type
NB.                  (symmetric for float data type)
NB.     HG         - Hessenberg-triangular pair, generalized
NB.                  Hessenberg form
NB.     HS,HSL,HSU - [{lower,upper}] Hessenberg
NB.     HT         - Hermitian (symmetric) tridiagonal
NB.     OR         - orthogonal
NB.     PO         - Hermitian (symmetric) positive definite
NB.     PT         - Hermitian (symmetric) positive definite
NB.                  tridiagonal
NB.     SY         - symmetric
NB.     TG         - triangular-triangular pair, generalized
NB.                  Schur form
NB.     TR         - triangular
NB.     TZ         - trapezoidal
NB.     UN         - unitary for complex data type
NB.                  (orthogonal for float data type)
NB.   triangular matrix [parts]:
NB.     L          -        lower triangular         matrix
NB.     L1         - unit   lower triangular         matrix
NB.     U          -        upper triangular         matrix
NB.     U1         - unit   upper triangular         matrix
NB.     LT         -        lower triangular part of matrix
NB.     SLT        - strict lower triangular part of matrix
NB.     UT         -        upper triangular part of matrix
NB.     SUT        - strict upper triangular part of matrix
NB.
NB. Notes:
NB. - unit diagonal in L1 and U1 is usually neither stored
NB.   nor referenced

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. ft4lisoa
NB.
NB. Description:
NB.   Functional table for linear ISO, all axes
NB.
NB. Syntax:
NB.   ft=. (op ft4lisoa) a
NB. where
NB.   op - dyad to compute a value from axes IO
NB.   a  - m×n-matrix
NB.   ft - the same shape as a, functional table:
NB.          ft -: m op/&i. n
NB.
NB. Examples:
NB.      < ft4lisoa 3 4 $ 0        > ft4lisoa 3 4 $ 0
NB.   0 1 1 1                   0 0 0 0
NB.   0 0 1 1                   1 0 0 0
NB.   0 0 0 1                   1 1 0 0
NB.
NB.      <: ft4lisoa 3 4 $ 0       >: ft4lisoa 3 4 $ 0
NB.   1 1 1 1                   1 0 0 0
NB.   0 1 1 1                   1 1 0 0
NB.   0 0 1 1                   1 1 1 0
NB.
NB.      - ft4lisoa 3 4 $ 0        -~ ft4lisoa 3 4 $ 0
NB.   0 _1 _2 _3                 0  1 2 3
NB.   1  0 _1 _2                _1  0 1 2
NB.   2  1  0 _1                _2 _1 0 1
NB.
NB.      + ft4lisoa 3 4 $ 0
NB.   0 1 2 3
NB.   1 2 3 4
NB.   2 3 4 5

ft4lisoa=: /(&i.)/(@$)

NB. ---------------------------------------------------------
NB. Miscellaneous

NB. conj. to extract matrix circumscribing the trapezoidal
NB. matrix starting from diagonal number x in the matrix y
trcut=: 2 : '(u@(0 >. <./"1)@v $) {. ]'

NB. extract upper trapezoidal matrix
trucut=: ]`-"0 trcut (] ,. (-~ {:))

NB. extract lower trapezoidal matrix
trlcut=: -`]"0 trcut ((+ {.) ,. ])

NB. conj. to extract trapezoidal matrix starting from
NB. diagonal number x in the circumscribing matrix y
tr=: 2 : '0&$: :([ ((u~ (-~ ft4lisoa)) {.`(0 ,: {:)}@,: ]) v)'

NB. ---------------------------------------------------------
NB. mxbstencil
NB.
NB. Description:
NB.   Adv. to make dyad returning [multi-][anti-]band
NB.   stencil for matrix
NB.
NB. Syntax:
NB.   S=. bs (vmix mxbstencil) A
NB. where
NB.   vmix - dyad to mix lISO x and y, is either (-~) for
NB.          band, or (+) for anti-band stencils, is called
NB.          as:
NB.            mix=. lIOrow vmix lIOcolumn
NB.   bs   - k×2-matrix of (b)s, or single b, or d, defines
NB.          [anti-]bands to stencil
NB.   b    - 2-vector (h,t), defines one [anti-]band to
NB.          stencil
NB.   h    - integer in range [-∞,t], lIO head of
NB.          [anti-]diagonal
NB.   t    - integer in range [h,+∞], lIO tail of
NB.          [anti-]diagonal
NB.   d    - integer in range [-∞,+∞], lIO single
NB.          [anti-]diagonal to stencil
NB.   A    - m×n-matrix
NB.   S    - m×n-matrix, boolean, having 1s on [anti-]band(s)
NB.
NB. Examples:
NB. - see mbstencil, mabstencil

mxbstencil=: 1 : '(+./^:(_2 + #@$)@:((1=I.)"1 2)~ -&1 0"1)~ (u ft4lisoa)'

NB. ---------------------------------------------------------
NB. mbstencil
NB. mabstencil
NB.
NB. Description:
NB.   [Multi-]band and [multi-]anti-band stencils for matrix
NB.
NB. Syntax:
NB.   S=. bs mbstencil  A
NB.   S=. bs mabstencil A
NB. where
NB.   bs - k×2-matrix of (b)s, or single b, or d, defines
NB.        [anti-]bands to stencil
NB.   A  - m×n-matrix
NB.   S  - m×n-matrix, boolean, having 1s on [anti-]band(s)
NB.   b  - 2-vector (h,t), defines one [anti-]band to stencil
NB.   h  - integer in range [-∞,t], defines lIO head of
NB.        [anti-]diagonal
NB.   t  - integer in range [h,+∞], defines lIO tail of
NB.        [anti-]diagonal
NB.   d  - integer in range [-∞,+∞], defines one
NB.        [anti-]diagonal to stencil
NB.
NB. Examples:
NB.    1 mbstencil i. 3 5                    1 mabstencil i. 3 5
NB. 0 1 0 0 0                             0 0 0 1 0
NB. 0 0 1 0 0                             0 0 1 0 0
NB. 0 0 0 1 0                             0 1 0 0 0
NB.
NB.    2 3 mbstencil i. 3 5                  2 3 mabstencil i. 3 5
NB. 0 0 1 1 0                             0 1 1 0 0
NB. 0 0 0 1 1                             1 1 0 0 0
NB. 0 0 0 0 1                             1 0 0 0 0
NB.
NB.    (__ _1 ,: 2 3) mbstencil i. 3 5       (__ _1 ,: 2 3) mabstencil i. 3 5
NB. 0 0 1 1 0                             0 1 1 0 0
NB. 1 0 0 1 1                             1 1 0 0 1
NB. 1 1 0 0 1                             1 0 0 1 1

mbstencil=:                   -~ mxbstencil
mabstencil=: (|."1@:-~ <:@c) (+  mxbstencil) ]

NB. ---------------------------------------------------------
NB. diagliso
NB.
NB. Description:
NB.   Return lISO solid part of diagonal of matrix
NB.
NB. Syntax:
NB.   liso=. [(d[,h[,s]])] diagliso [m,]n
NB. where
NB.   m    ≥ 0, integer, optional rows in matrix, default is
NB.          n
NB.   n    ≥ 0, integer, columns in matrix
NB.   d    - integer in range [1-m,n-1], optional lIO
NB.          diagonal, default is 0 (main diagonal)
NB.   h    - integer in range [-S,S-1], optional lIO extreme
NB.          element of solid part of diagonal, default is 0
NB.          (take from head)
NB.   s    - integer in range [-S,S] or ±∞, optional size of
NB.          solid part of diagonal, default is +∞ (all
NB.          elements in forward direction)
NB.   liso - min(S,|s|)-vector of integers, lISO solid
NB.          part of diagonal
NB.   S    ≥ 0, the length of diagonal
NB.
NB. Formulae:
NB. - the whole diagonal's lIO extreme element:
NB.     H := (d ≥ 0) ? d : (-n*d)
NB. - the whole diagonal's size:
NB.     S := max(0,min(m,n,⌊(n+m-|n-m-2*d|)/2⌋))
NB.
NB. Notes:
NB. - (h,s) pair defines raveled rISO solid part of
NB.   diagonal

diagliso=: 0 0 _&$: :(4 : 0)
  'd h s'=. x=. ((i. 3) < (# x))} 0 0 _ ,: x  NB. in-place op
  'm n'=. y=. 2 $ y
  H=. n (-@*^:(0 > ])) d
  S=. 0 >. <./ y , <. -: (n + m - | n - m + +: d)
  (h ,: s <. S) ];.0 (>: n) liso4dhs H , S
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Miscellaneous

c=: {:!.1@$      NB. Columns in noun

trace=: +/@diag  NB. Matrix trace

ct=: +@:|:       NB. Conjugate transpose
cp=: ct&.|.      NB. Conjugate pertranspose

NB. Do/undo full (symmetric) permutation
NB. Syntax:
NB.   Aperm=. p fp     A
NB.   A=.     p fp^:_1 Aperm
fp=: ([ C."1 C.) :. ([ C.^:_1"1 C.^:_1)

NB. Transform permutation vector to/from permutation matrix,
NB. to permute rows by y or columns by (/: y)
P4p=:  ({    =) :. (     i.&1"1 )

NB. Transform inversed permutation vector to/from permutation
NB. matrix, or permutation vector to/from inversed
NB. permutation matrix, to permute rows by (/: y) or columns
NB. by y
P4ip=: ({^:_1=) :. (/:@:(i.&1"1))

NB. ---------------------------------------------------------
NB. icut
NB.
NB. Description:
NB.   Inversed cut to model <;.1^:_1
NB.
NB. Syntax:
NB.   A=. icut bA
NB. where
NB.   bA - block array
NB.   A  - sh-array
NB.
NB. Assertions:
NB.   A -: icut fret <;.1 A
NB. where
NB.   A    - some array
NB.   fret - some fret
NB.
NB. References:
NB. [1] Roger Hui. JWiki/Essays/Block Matrix Inverse.
NB.     2005-11-24 03:53:19.
NB.     http://code.jsoftware.com/wiki/Essays/Block%20Matrix%20Inverse
NB.
NB. TODO:
NB. - fret would be sparse

icut=: [: > 3 : ',"(#$y)&.>/y'^:(#@$)

NB. ---------------------------------------------------------
NB. rt
NB.
NB. Description:
NB.   Restrained Take. Just like built-in Take verb ({.), but
NB.   without overtake feature. Overtaking value means "all
NB.   elements along this axis".
NB.
NB. Examples:
NB.    2 rt i. 3 4                  _2 rt i. 3 4
NB. 0 1 2 3                      4 5  6  7
NB. 4 5 6 7                      8 9 10 11
NB.
NB.    2 30 rt i. 3 4               2 _ rt i. 3 4
NB. 0 1 2 3                      0 1 2 3
NB. 4 5 6 7                      4 5 6 7
NB.
NB.    20 3 rt i. 3 4               _ 3 rt i. 3 4
NB. 0 1  2                       0 1  2
NB. 4 5  6                       4 5  6
NB. 8 9 10                       8 9 10
NB.
NB.    _2 _30 rt i. 3 4             _2 __ rt i. 3 4
NB. 4 5  6  7                    4 5  6  7
NB. 8 9 10 11                    8 9 10 11
NB.
NB.    _20 _3 rt i. 3 4             __ _3 rt i. 3 4
NB. 1  2  3                      1  2  3
NB. 5  6  7                      5  6  7
NB. 9 10 11                      9 10 11

rt=: (*@[ * |@[ <. (({.~ #)~ $)) {. ]

NB. ---------------------------------------------------------
NB. upd
NB.
NB. Description:
NB.   Adv. to update subarray by a monad
NB.
NB. Syntax:
NB.   vapp=. u upd
NB. where
NB.   u    - monad to update subA; is called as:
NB.            subAupd=. u subA
NB.   vapp - verb to update A; is called as:
NB.            Aupd=. iso vapp A
NB.   iso  - ISO subA in the A
NB.   subA - subarray in the A
NB.   A    - array
NB.   Aupd - A with subA replaced by subAupd
NB.
NB. Assertions:
NB.   Aupd -: iso vapp A
NB. where
NB.   vapp=. u upd
NB.   subA=. iso { A
NB.   subAupd=. u subA
NB.   Aupd=. subAupd iso} A
NB.
NB. Examples:
NB.    1 ('*'"_ upd) 4 5$'-'             1 2 ('*'"_ upd) 4 5$'-'
NB. -----                             -----
NB. *****                             *****
NB. -----                             *****
NB. -----                             -----
NB.
NB.    (<1 2) ('*'"_ upd) 4 5$'-'        (1 2;2 3) ('*'"_ upd) 4 5$'-'
NB. -----                             -----
NB. --*--                             --*--
NB. -----                             ---*-
NB. -----                             -----
NB.
NB. References:
NB. [1] Oleg Kobchenko. [Jprogramming] Transform to Amend.
NB.     2007-03-02 22:15:54 HKT.
NB.     http://www.jsoftware.com/pipermail/programming/2007-March/005415.html

upd=: (@:{) (`[) (`])}

NB. ---------------------------------------------------------
NB. e0
NB.
NB. Description:
NB.   Extend matrix by blanks
NB.
NB. Syntax:
NB.   eA=. sh e0 A
NB. where
NB.   A  - matrix to extend
NB.   sh - scalar or 2-vector, integer, extended size
NB.   eA - extend A
NB.
NB. Examples:
NB.    2 e0 3 4 $ 1        2 3 e0 3 4 $ 1        2 _3 e0 3 4 $ 1
NB. 1 1 1 1 0 0         1 1 1 1 0 0 0         0 0 0 1 1 1 1
NB. 1 1 1 1 0 0         1 1 1 1 0 0 0         0 0 0 1 1 1 1
NB. 1 1 1 1 0 0         1 1 1 1 0 0 0         0 0 0 1 1 1 1
NB. 0 0 0 0 0 0         0 0 0 0 0 0 0         0 0 0 0 0 0 0
NB. 0 0 0 0 0 0         0 0 0 0 0 0 0         0 0 0 0 0 0 0
NB.
NB.    _2 e0 3 4 $ 1       _2 3 e0 3 4 $ 1       _2 _3 e0 3 4 $ 1
NB. 0 0 0 0 0 0         0 0 0 0 0 0 0         0 0 0 0 0 0 0
NB. 0 0 0 0 0 0         0 0 0 0 0 0 0         0 0 0 0 0 0 0
NB. 0 0 1 1 1 1         1 1 1 1 0 0 0         0 0 0 1 1 1 1
NB. 0 0 1 1 1 1         1 1 1 1 0 0 0         0 0 0 1 1 1 1
NB. 0 0 1 1 1 1         1 1 1 1 0 0 0         0 0 0 1 1 1 1
NB.       _1 _2 e0 3 4 $ 'abcdefghijkl'
NB.
NB.   abcd
NB.   efgh
NB.   ijkl

e0=: ([ + (negneg"0 $)) {. ]

NB. ---------------------------------------------------------
NB. appendl
NB. appendr
NB.
NB. Description:
NB.   Enhance built-in Append verb (,)
NB.
NB. Syntax:
NB.   B=. A0 appendx A1
NB.
NB. Examples:
NB.    (3 3$3) appendl (2 2$2)      (3 3$3) appendr (2 2$2)
NB. 3 3 3                        3 3 3
NB. 3 3 3                        3 3 3
NB. 3 3 3                        3 3 3
NB. 2 2 0                        0 2 2
NB. 2 2 0                        0 2 2
NB.
NB.    (2 2$2) appendl (3 3$3)      (2 2$2) appendr (3 3$3)
NB. 2 2 0                        0 2 2
NB. 2 2 0                        0 2 2
NB. 3 3 3                        3 3 3
NB. 3 3 3                        3 3 3
NB. 3 3 3                        3 3 3
NB.
NB. Notes:
NB. - at most one of A0, A1 can be 1-rank array (i.e. vector)

appendl=: , `([, ({."1~   c)~)`(({."1~   c), ])@.(*@-&c)
appendr=: , `([, ({."1~ -@c)~)`(({."1~ -@c), ])@.(*@-&c)

NB. ---------------------------------------------------------
NB. stitcht
NB. stitchb
NB.
NB. Description:
NB.   Enhance built-in Stitch verb (,.)
NB.
NB. Syntax:
NB.   B=. A0 stitchx A1
NB.
NB. Examples:
NB.    (3 3$3) stitcht (2 2$2)      (3 3$3) stitchb (2 2$2)
NB. 3 3 3 2 2                    3 3 3 0 0
NB. 3 3 3 2 2                    3 3 3 2 2
NB. 3 3 3 0 0                    3 3 3 2 2
NB.
NB.    (2 2$2) stitcht (3 3$3)      (2 2$2) stitchb (3 3$3)
NB. 2 2 3 3 3                    0 0 3 3 3
NB. 2 2 3 3 3                    2 2 3 3 3
NB. 0 0 3 3 3                    2 2 3 3 3
NB.
NB. Notes:
NB. - 1-rank arrays (i.e. vectors) are also acceptable

stitcht=: ,.`([,.({.  ~   #)~)`(({.  ~   #),.])@.(*@-&#)
stitchb=: ,.`([,.({.  ~ -@#)~)`(({.  ~ -@#),.])@.(*@-&#)

NB. ---------------------------------------------------------
NB. ds
NB.
NB. Description:
NB.   Direct sum of matrices A⊕B
NB.
NB. Syntax:
NB.   C=. A ds B
NB. where
NB.   A - ma×na-matrix
NB.   B - mb×nb-matrix
NB.   C - (ma+mb)×(na+nb)-matrix:
NB.         C = (  A 0  )
NB.             (  0 B  )

ds=: (+&c {."1 [) appendr ]

NB. ---------------------------------------------------------
NB. diag
NB.
NB. Description:
NB.   Return a solid part of diagonal of matrix
NB.
NB. Syntax:
NB.   e=. [(d[,h[,s]])] diag A
NB. where
NB.   A - m×n-matrix
NB.   d - integer in range [1-m,n-1], optional lIO diagonal,
NB.       default is 0 (main diagonal)
NB.   h - integer in range [-S,S-1], optional lIO extreme
NB.       element of solid part of diagonal, default is 0
NB.       (take from head)
NB.   s - integer in range [-S,S] or ±∞, optional size of
NB.       solid part of diagonal, default is +∞ (all elements
NB.       in forward direction)
NB.   e - min(S,|s|)-vector, elements from the solid part of
NB.       diagonal
NB.   S ≥ 0, the length of diagonal

diag=: ((<0 1)&|:) :((diagliso $) ({,) ])

NB. ---------------------------------------------------------
NB. setdiag
NB.
NB. Description:
NB.   Assign value(s) to a solid part of diagonal
NB.
NB. Syntax:
NB.   Aupd=. (e;[d[,h[,s]]]) setdiag A
NB. where
NB.   A    - m×n-matrix to change
NB.   e    - {0,1}-rank array, value(s) to assign
NB.   d    - integer in range [1-m,n-1], optional lIO
NB.          diagonal, default is 0 (main diagonal)
NB.   h    - integer in range [-S,S-1], optional lIO extreme
NB.          element of solid part of diagonal, default is 0
NB.          (take from head)
NB.   s    - integer in range [-S,S] or ±∞ when e is scalar,
NB.          or any from set {±k,±∞} when e is vector;
NB.          optional size of solid part of diagonal, default
NB.          is +∞ (all elements in forward direction)
NB.   Aupd - m×n-matrix A with value(s) e assigned to solid
NB.          part of d-th diagonal
NB.   S    ≥ 0, the length of d-th diagonal
NB.   k    ≤ S, the length of vector e
NB.
NB. Examples:
NB.    (2;a:) setdiag 4 4 $ 0         (2;_1 1 1) setdiag 4 4 $ 0
NB. 2 0 0 0                        0 0 0 0
NB. 0 2 0 0                        0 0 0 0
NB. 0 0 2 0                        0 2 0 0
NB. 0 0 0 2                        0 0 0 0
NB.
NB.    (2;_1) setdiag 4 4 $ 0         (1 2 3;_1) setdiag 4 4 $ 0
NB. 0 0 0 0                        0 0 0 0
NB. 2 0 0 0                        1 0 0 0
NB. 0 2 0 0                        0 2 0 0
NB. 0 0 2 0                        0 0 3 0
NB.
NB.    (2;_1 1) setdiag 4 4 $ 0       (1 2 3;_1 _1 _3) setdiag 4 4 $ 0
NB. 0 0 0 0                        0 0 0 0
NB. 0 0 0 0                        3 0 0 0
NB. 0 2 0 0                        0 2 0 0
NB. 0 0 2 0                        0 0 1 0

setdiag=: 4 : 0
  'e dhs'=. x
  dhs=. ((i. 3) < (# dhs))} 0 0 _ ,: dhs  NB. assign defaults, in-place op
  liso=. dhs diagliso $ y
  e (liso"_)} y
)

NB. ---------------------------------------------------------
NB. upddiag
NB.
NB. Description:
NB.   Adv. to make verbs to update a solid part of diagonal
NB.
NB. Syntax:
NB.   Aupd=. [(d,[h[,s]])] (u upddiag) A
NB. where
NB.   u    - monad to change elements; is called as:
NB.            eupd=. u e
NB.   d    - integer in range [1-m,n-1], optional lIO
NB.          diagonal, default is 0 (main diagonal)
NB.   h    - integer in range [-S,S-1], optional lIO extreme
NB.          element of solid part of diagonal, default is 0
NB.          (take from head)
NB.   s    - integer in range [-S,S] or ±∞, optional size of
NB.          solid part of diagonal, default is +∞ (all
NB.          elements in forward direction)
NB.   A    - m×n-matrix to update
NB.   Aupd - A with solid part of d-th diagonal updated by
NB.          monad u
NB.   S    ≥ 0, the length of d-th diagonal

upddiag=: 1 : 'diagliso_mt_^:(1:`(] $)) (u {{(u x ({,) y) (x"_)} y}}) ]'

NB. ---------------------------------------------------------
NB. bdlpick
NB.
NB. Description:
NB.   Zeroize elements outside lower bidiagonal part of the
NB.   matrix
NB.
NB. Syntax:
NB.   B=. bdlpick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   B - m×n-matrix, the lower bidiagonal
NB.
NB. TODO:
NB. - B would be sparse

bdlpick=: _1 0&mbstencil`(0&,:)}

NB. ---------------------------------------------------------
NB. bdupick
NB.
NB. Description:
NB.   Zeroize elements outside upper bidiagonal part of the
NB.   matrix
NB.
NB. Syntax:
NB.   B=. bdupick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   B - m×n-matrix, the upper bidiagonal
NB.
NB. TODO:
NB. - B would be sparse

bdupick=: 0 1&mbstencil`(0&,:)}

NB. ---------------------------------------------------------
NB. hslpick
NB.
NB. Description:
NB.   Zeroize elements outside lower Hessenberg part of the
NB.   matrix
NB.
NB. Syntax:
NB.   B=. hslpick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   B - m×n-matrix, the lower Hessenberg

hslpick=: __ 1&mbstencil`(0&,:)}

NB. ---------------------------------------------------------
NB. hsupick
NB.
NB. Description:
NB.   Zeroize elements outside upper Hessenberg part of the
NB.   matrix
NB.
NB. Syntax:
NB.   B=. hsupick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   B - m×n-matrix, the upper Hessenberg

hsupick=: _1 _&mbstencil`(0&,:)}

NB. ---------------------------------------------------------
NB. gtpick
NB.
NB. Description:
NB.   Zeroize elements outside tridiagonal part of the matrix
NB.
NB. Syntax:
NB.   B=. gtpick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   B - m×n-matrix, tridiagonal
NB.
NB. TODO:
NB. - B would be sparse

gtpick=: _1 1&mbstencil`(0&,:)}

NB. ---------------------------------------------------------
NB. trlpick
NB.
NB. Description:
NB.   Zeroize elements outside lower trapezoidal part of the
NB.   matrix
NB.
NB. Syntax:
NB.   B=. [d] trlpick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   d - integer in range [-∞,+∞], optional lIO last
NB.       non-zero diagonal, default is 0
NB.   B - m×n-matrix, the lower trapezoidal

trlpick=: (>: ft4lisoa)`(0&,:)} :(4 : '((__ , x) mbstencil ])`(0 ,: ])} y')

NB. ---------------------------------------------------------
NB. trupick
NB.
NB. Description:
NB.   Zeroize elements outside upper trapezoidal part of the
NB.   matrix
NB.
NB. Syntax:
NB.   B=. [d] trupick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   d - integer in range [-∞,+∞], lIO first non-zero
NB.       diagonal, default is 0
NB.   B - m×n-matrix, the upper trapezoidal

trupick=: (<: ft4lisoa)`(0&,:)} :(4 : '((x , _) mbstencil ])`(0 ,: ])} y')

NB. ---------------------------------------------------------
NB. trl1pick
NB.
NB. Description:
NB.   Zeroize elements outside lower trapezoidal part of the
NB.   matrix and set diagonal to 1
NB.
NB. Syntax:
NB.   B=. [d] trl1pick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   d - integer in range [-∞,+∞], optional lIO last
NB.       non-zero diagonal, default is 0
NB.   B - m×n-matrix, the lower trapezoidal with unit on
NB.       diagonal d

trl1pick=: 0&$: :((1 ; [) setdiag trlpick)

NB. ---------------------------------------------------------
NB. tru1pick
NB.
NB. Description:
NB.   Zeroize elements outside upper trapezoidal part of the
NB.   matrix and set diagonal to 1
NB.
NB. Syntax:
NB.   B=. [d] trupick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   d - integer in range [-∞,+∞], optional lIO first
NB.       non-zero diagonal, default is 0
NB.   B - m×n-matrix, the upper trapezoidal with unit on
NB.       diagonal d

tru1pick=: 0&$: :((1 ; [) setdiag trupick)

NB. ---------------------------------------------------------
NB. idmat
NB.
NB. Description:
NB.   Make identity matrix with units on solid part of
NB.   diagonal
NB.
NB. Syntax:
NB.   I=. [(d[,h[,s]])] idmat [m,]n
NB. where
NB.   m ≥ 0, integer, optional rows in matrix I, default is n
NB.   n ≥ 0, integer, columns in matrix I
NB.   d - integer in range [1-m,n-1], optional lIO diagonal,
NB.       default is 0 (main diagonal)
NB.   h - integer in range [-S,S-1], optional lIO extreme
NB.       element of solid part of diagonal, default is 0
NB.       (take from head)
NB.   s - integer in range [-S,S] or ±∞, optional size of
NB.       solid part of diagonal, default is +∞ (all elements
NB.       in forward direction)
NB.   I - m×n-matrix of zeros with unit assigned to solid
NB.       part of d-th diagonal
NB.   S ≥ 0, the length of d-th diagonal
NB.
NB. Examples:
NB.    idmat 3                      idmat 3 4
NB. 1 0 0                        1 0 0 0
NB. 0 1 0                        0 1 0 0
NB. 0 0 1                        0 0 1 0
NB.
NB.    1 idmat 3 4                  _1 idmat 3 4
NB. 0 1 0 0                      0 0 0 0
NB. 0 0 1 0                      1 0 0 0
NB. 0 0 0 1                      0 1 0 0
NB.
NB. TODO:
NB. - I would be sparse

idmat=: a:&$: :((1 ; [) setdiag (0 $~ 2 $ ]))

NB. ---------------------------------------------------------
NB. diagmat
NB.
NB. Description:
NB.   Make diagonal matrix
NB.
NB. Syntax:
NB.   D=. [(h,t)] diagmat e
NB. where
NB.   e - S-vector, new values for diagonal
NB.   h - integer in range [1-m,n-1], optional lIO diagonal
NB.       of v's head, relatively to top left corner, default
NB.       is 0
NB.   t - integer in range [1-m,n-1], optional lIO diagonal
NB.       of v's tail, relatively to bottom right corner,
NB.       default is 0
NB.   D - m×n-matrix of zeros with vector e assigned to h-th
NB.       diagonal
NB.   S ≥ 0, the length of h-th diagonal
NB.
NB. Algorithm:
NB.   1) find D shape
NB.   2) generate lIO h-th diagonal
NB.   3) write e into matrix of zeros
NB.
NB. Examples:
NB.    diagmat 3 5 7            0 0 diagmat 3 5 7
NB. 3 0 0                    3 0 0
NB. 0 5 0                    0 5 0
NB. 0 0 7                    0 0 7
NB.
NB.    1 0 diagmat 3 5 7        _1 0 diagmat 3 5 7
NB. 0 3 0 0                  0 0 0
NB. 0 0 5 0                  3 0 0
NB. 0 0 0 7                  0 5 0
NB.                          0 0 7
NB.
NB.    0 1 diagmat 3 5 7        0 _1 diagmat 3 5 7
NB. 3 0 0                    3 0 0 0
NB. 0 5 0                    0 5 0 0
NB. 0 0 7                    0 0 7 0
NB. 0 0 0
NB.
NB.    1 1 diagmat 3 5 7        _1 _1 diagmat 3 5 7
NB. 0 3 0 0                  0 0 0 0
NB. 0 0 5 0                  3 0 0 0
NB. 0 0 0 7                  0 5 0 0
NB. 0 0 0 0                  0 0 7 0
NB.
NB.    1 _1 diagmat 3 5 7       _1 1 diagmat 3 5 7
NB. 0 3 0 0 0                0 0 0
NB. 0 0 5 0 0                3 0 0
NB. 0 0 0 7 0                0 5 0
NB.                          0 0 7
NB.                          0 0 0
NB.
NB. TODO:
NB. - D would be sparse

diagmat=: 0 0&$: :((; {.)~ setdiag ((0 $~ (+ (2&(|.@}. - {.)@(0&(<. , >.)))))~ #))

NB. ---------------------------------------------------------
NB. trl
NB.
NB. Description:
NB.   Extract lower trapezoidal matrix with optional
NB.   shrinking
NB.
NB. Examples:
NB.    trl >: i. 3 4         0 trl >: i. 3 4
NB. 1  0  0               1  0  0
NB. 5  6  0               5  6  0
NB. 9 10 11               9 10 11
NB.
NB.    1 trl >: i. 3 4       _1 trl >: i. 3 4
NB. 1  2  0  0            5  0
NB. 5  6  7  0            9 10
NB. 9 10 11 12
NB.
NB.    1 trl >: i. 4 3       _1 trl >: i. 4 3
NB.  1  2  0               4  0  0
NB.  4  5  6               7  8  0
NB.  7  8  9              10 11 12
NB. 10 11 12

trl=: (>:~ 0&>.) tr trlcut

NB. ---------------------------------------------------------
NB. tru
NB.
NB. Description:
NB.   Extract upper trapezoidal matrix with optional
NB.   shrinking
NB.
NB. Examples:
NB.    tru >: i. 3 4         0 tru >: i. 3 4
NB. 1 2  3  4             1 2  3  4
NB. 0 6  7  8             0 6  7  8
NB. 0 0 11 12             0 0 11 12
NB.
NB.    1 tru >: i. 3 4       _1 tru >: i. 3 4
NB. 2 3  4                1  2  3  4
NB. 0 7  8                5  6  7  8
NB. 0 0 12                0 10 11 12
NB.
NB.    1 tru >: i. 4 3       _1 tru >: i. 4 3
NB. 2 3                   1 2  3
NB. 0 6                   4 5  6
NB.                       0 8  9
NB.                       0 0 12

tru=: (<:~ 0&<.) tr trucut

NB. ---------------------------------------------------------
NB. trl0
NB.
NB. Description:
NB.   Extract strictly lower trapezoidal matrix with optional
NB.   shrinking
NB.
NB. Examples:
NB.    trl0 >: i. 4 3         0 trl0 >: i. 4 3
NB.  0  0  0                0  0  0
NB.  4  0  0                4  0  0
NB.  7  8  0                7  8  0
NB. 10 11 12               10 11 12
NB.
NB.    1 trl0 >: i. 4 3       _1 trl0 >: i. 4 3
NB.  1  0  0                0  0 0
NB.  4  5  0                7  0 0
NB.  7  8  9               10 11 0
NB. 10 11 12
NB.
NB.    1 trl0 >: i. 3 4       _1 trl0 >: i. 3 4
NB. 1  0  0 0              0 0
NB. 5  6  0 0              9 0
NB. 9 10 11 0

trl0=: (>~ 0&>.) tr trlcut

NB. ---------------------------------------------------------
NB. tru0
NB.
NB. Description:
NB.   Extract strictly upper trapezoidal matrix with optional
NB.   shrinking
NB.
NB. Examples:
NB.    tru0 >: i. 3 4         0 tru0 >: i. 3 4
NB. 0 2 3  4               0 2 3  4
NB. 0 0 7  8               0 0 7  8
NB. 0 0 0 12               0 0 0 12
NB.
NB.    1 tru0 >: i. 3 4       _1 tru0 >: i. 3 4
NB. 0 3 4                  1 2  3  4
NB. 0 0 8                  0 6  7  8
NB. 0 0 0                  0 0 11 12
NB.
NB.    1 tru0 >: i. 4 3       _1 tru0 >: i. 4 3
NB. 0 3                    1 2 3
NB. 0 0                    0 5 6
NB.                        0 0 9
NB.                        0 0 0

tru0=: (<~ 0&<.) tr trucut

NB. ---------------------------------------------------------
NB. trl1
NB.
NB. Description:
NB.   Extract unit lower trapezoidal matrix with optional
NB.   shrinking
NB.
NB. Examples:
NB.    trl1 >: i. 4 3         0 trl1 >: i. 4 3
NB.  1  0  0                1  0  0
NB.  4  1  0                4  1  0
NB.  7  8  1                7  8  1
NB. 10 11 12               10 11 12
NB.
NB.    1 trl1 >: i. 4 3       _1 trl1 >: i. 4 3
NB.  1  1  0                1  0 0
NB.  4  5  1                7  1 0
NB.  7  8  9               10 11 1
NB. 10 11 12
NB.
NB.    1 trl1 >: i. 3 4       _1 trl1 >: i. 3 4
NB. 1  1  0 0              1 0
NB. 5  6  1 0              9 1
NB. 9 10 11 1

trl1=: 0&$: :([ trl (1 ; [) setdiag ])

NB. ---------------------------------------------------------
NB. tru1
NB.
NB. Description:
NB.   Extract unit upper trapezoidal matrix with optional
NB.   shrinking
NB.
NB. Examples:
NB.    tru1 >: i. 3 4         0 tru1 >: i. 3 4
NB. 1 2 3  4               1 2 3  4
NB. 0 1 7  8               0 1 7  8
NB. 0 0 1 12               0 0 1 12
NB.
NB.    1 tru1 >: i. 3 4       _1 tru1 >: i. 3 4
NB. 1 3 4                  1 2  3  4
NB. 0 1 8                  1 6  7  8
NB. 0 0 1                  0 1 11 12
NB.
NB.    1 tru1 >: i. 4 3       _1 tru1 >: i. 4 3
NB. 1 3                    1 2 3
NB. 0 1                    1 5 6
NB.                        0 1 9
NB.                        0 0 1

tru1=: 0&$: :([ tru (1 ; [) setdiag ])

NB. ---------------------------------------------------------
NB. Verb      Reads    Overwrites    Diagonal        Composes matrix
NB. sy4gel    SLT      SUT           leave as is     SY
NB. sy4geu    SUT      SLT           leave as is     SY
NB. he4gel    SLT       UT           reificate       HE
NB. he4geu    SUT       LT           reificate       HE
NB. ss4gel    SLT       UT           zeroize         skew-symmetric
NB. ss4geu    SUT       LT           zeroize         skew-symmetric
NB. sh4gel    SLT       UT           zeroize         skew-Hermitian
NB. sh4geu    SUT       LT           zeroize         skew-Hermitian
NB.
NB. Description:
NB.   Compose structured matrix from SLT (SUT) part and
NB.   diagonal of general square matrix
NB.
NB. Syntax:
NB.   S=. xx4gex G
NB. where
NB.   G - n×n-matrix
NB.   S - the same shape as G, structured

sy4gel=:                  (</~@i.@#)`(,:   |:)}
sy4geu=:                  (>/~@i.@#)`(,:   |:)}
he4gel=: 0 (9&o. upddiag) (</~@i.@#)`(,:   ct)}
he4geu=: 0 (9&o. upddiag) (>/~@i.@#)`(,:   ct)}
ss4gel=: (0 ; a:) setdiag (</~@i.@#)`(,: -@|:)}
ss4geu=: (0 ; a:) setdiag (>/~@i.@#)`(,: -@|:)}
sh4gel=: (0 ; a:) setdiag (</~@i.@#)`(,: -@ct)}
sh4geu=: (0 ; a:) setdiag (>/~@i.@#)`(,: -@ct)}

NB. ---------------------------------------------------------
NB. Actor     P.o.S.     x arg goes to    y arg goes to    Diagonal comes from
NB. lxsuy     verb        LT              SUT              x
NB. slxuy     verb       SLT               UT              y
NB. suxly     verb       SUT               LT              y
NB. uxsly     verb        UT              SLT              x
NB. slxsuy    adverb     SLT              SUT              m
NB. suxsly    adverb     SUT              SLT              m
NB.
NB. Description:
NB.   Compose matrix from triangular parts of general
NB.   matrices
NB.
NB. Syntax:
NB.   C=. A    xxxxxy  B
NB.   D=. A (d sxxsxy) B
NB. where
NB.   A - m×n-matrix or scalar
NB.   B - m×n-matrix or scalar
NB.   d - scalar to place in diagonal
NB.   C - m×n-matrix composed from A and B triangular parts
NB.   D - m×n-matrix composed from A and B strict triangular
NB.       parts and d as diagonal elements
NB.
NB. Notes:
NB. - at most one of A, B can be scalar
NB.
NB. Examples:
NB.      ] 'X Y'=. 4 6 ;/@:($"1 0) 'xy'
NB.   +------+------+
NB.   |xxxxxx|yyyyyy|
NB.   |xxxxxx|yyyyyy|
NB.   |xxxxxx|yyyyyy|
NB.   |xxxxxx|yyyyyy|
NB.   +------+------+
NB.      ] A=. (1 1 ,: 4 6) ];.0 i. 10 10
NB.   11 12 13 14 15 16
NB.   21 22 23 24 25 26
NB.   31 32 33 34 35 36
NB.   41 42 43 44 45 46
NB.
NB.      X lxsuy Y       X slxuy Y       X suxly Y       X uxsly Y
NB.   xyyyyy          yyyyyy          yxxxxx          xxxxxx
NB.   xxyyyy          xyyyyy          yyxxxx          yxxxxx
NB.   xxxyyy          xxyyyy          yyyxxx          yyxxxx
NB.   xxxxyy          xxxyyy          yyyyxx          yyyxxx
NB.
NB.      X '\' slxsuy Y                  X '\' suxsly Y
NB.   \yyyyy                          \xxxxx
NB.   x\yyyy                          y\xxxx
NB.   xx\yyy                          yy\xxx
NB.   xxx\yy                          yyy\xx
NB.
NB.      NB. simulate trlpick
NB.      trlpick A          A lxsuy 0          0 suxly A
NB.   11  0  0  0 0 0    11  0  0  0 0 0    11  0  0  0 0 0
NB.   21 22  0  0 0 0    21 22  0  0 0 0    21 22  0  0 0 0
NB.   31 32 33  0 0 0    31 32 33  0 0 0    31 32 33  0 0 0
NB.   41 42 43 44 0 0    41 42 43 44 0 0    41 42 43 44 0 0
NB.
NB.      NB. simulate trl1pick
NB.      trl1pick A         A 1 slxsuy 0       0 (1 suxsly) A
NB.    1  0  0 0 0 0      1  0  0 0 0 0      1  0  0 0 0 0
NB.   21  1  0 0 0 0     21  1  0 0 0 0     21  1  0 0 0 0
NB.   31 32  1 0 0 0     31 32  1 0 0 0     31 32  1 0 0 0
NB.   41 42 43 1 0 0     41 42 43 1 0 0     41 42 43 1 0 0

lxsuy=: < /&i./@}.@$}@,:
slxuy=: <:/&i./@}.@$}@,:
suxly=: >:/&i./@}.@$}@,:
uxsly=: > /&i./@}.@$}@,:

slxsuy=: 1 : '*@:  (-/&i./)@}.@$}@(m , ,:)'
suxsly=: 1 : '*@:-@(-/&i./)@}.@$}@(m , ,:)'

NB. ---------------------------------------------------------
NB. po
NB.
NB. Description:
NB.   Make Hermitian (symmetric) positive definite matrix
NB.   from general square invertible one
NB.
NB. Syntax:
NB.   P=. po G
NB. where
NB.   G - n×n-matrix, invertible
NB.   H - n×n-matrix, the Hermitian (symmetric)

po=: mp ct
