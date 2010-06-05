NB. struct.ijs
NB. Structure handlers
NB.
NB. c         Columns in matrix
NB. trace     Matrix trace
NB. ct        Conjugate transpose
NB. pt        Pertranspose
NB. cpt       Conjugate pertranspose
NB. pp        Apply permutation to both rows and columns
NB. p2P       Transform permutation vector to permutation
NB.           matrix
NB. ip2P      Transform inversed permutation vector to
NB.           permutation matrix
NB. rt        Restrained Take
NB.
NB. upd1      Adv. to update subarray by a monad
NB. append    Template adv. to make verbs to enhance append
NB.           built-in verb (,)
NB. stitch    Template adv. to make verbs to enhance stitch
NB.           built-in verb (,.)
NB.
NB. diag      Return a solid part of diagonal
NB. setdiag   Assign scalar value to a solid part of diagonal
NB. upddiag   Template adv. to make verbs to update a solid
NB.           part of diagonal
NB.
NB. idmat     Make identity matrix with units on solid part
NB.           of diagonal
NB. diagmat   Make diagonal matrix
NB.
NB. tru       Extract upper triangular (trapezoidal) matrix
NB. trl       Extract lower triangular (trapezoidal) matrix
NB. tru0      Extract strictly upper triangular (trapezoidal)
NB.           matrix
NB. trl0      Extract strictly lower triangular (trapezoidal)
NB.           matrix
NB. tru1      Extract unit upper triangular (trapezoidal)
NB.           matrix
NB. trl1      Extract unit lower triangular (trapezoidal)
NB.           matrix
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Concepts
NB.
NB. IO   - index of
NB. lIO  - linear IO, is an integer
NB. IOS  - indices of
NB. lIOS - linear IOS, is a vector of integers
NB. rIOS - rectangular IOS, for r-rank array is a 2×r-array
NB.        of integers ((from0,from1,...),:(size0,size1,...))
NB.
NB. Following are equivalents:
NB.   (3 5 _7,:2 _3 4) (] ;. 0) report
NB.   (< 3 4;7 6 5;_10 _9 _8 _7) { report
NB.   (rios2ios (3 5 _7,:2 _3 4)) { report
NB.   (ios2rios (< 3 4;7 6 5;_10 _9 _8 _7)) (] ;. 0) report
NB.
NB. Following are equivalents:
NB.   (0 1 ; 1 2 ; 2 3) { i. 3 4
NB.   (4 lios2ios 1 6 11) { i. 3 4
NB.   (4 ios2lios (0 1 ; 1 2 ; 2 3)) ({,) i. 3 4
NB.   1 6 11 ({,) i. 3 4

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Misc.

NB. convert shape y to IOS differences table
sh2id=: {. -~/&i. {:

NB. template conj. to extract rectangular matrix
NB. circumscribing the triangular (trapezoidal) matrix
NB. starting from diagonal number x in the matrix y
trcut=: 2 : '((m & *) @: (<./ " 1) @ v $) {. ]'

NB. extract upper triangular (trapezoidal) matrix
trucut=: 1 _1 trcut (] ,. (-~ {:))

NB. extract lower triangular (trapezoidal) matrix
trlcut=: _1 1 trcut ((+ {.) ,. ])

NB. template conj. to extract triangular (trapezoidal)
NB. matrix starting from diagonal number x in the rectangular
NB. circumscribing matrix y
tr=: 2 : '0&$: : ([ (] * (u~ sh2id@$)) v)'

NB. ---------------------------------------------------------
NB. diaglios
NB.
NB. Description:
NB.   Return lIOS of solid part of diagonal of rectangular
NB.   matrix
NB.
NB. Syntax:
NB.   lios=. [(d[,f[,s]])] diaglios [m,]n
NB. where
NB.   m    ≥ 0, integer, optional rows in matrix, default is
NB.          n
NB.   n    ≥ 0, integer, columns in matrix
NB.   d    - integer in range [1-m,n-1], optional IO
NB.          diagonal, default is 0 (main diagonal)
NB.   f    - integer in range [-S,S-1], optional IO extreme
NB.          element of solid part within diagonal,
NB.          default is 0 (take from head)
NB.   s    - integer in range [-S,S] or ±∞, optional size of
NB.          solid part within diagonal, default is +∞ (all
NB.          elements in forward direction)
NB.   lios - min(S,|s|)-vector of integers, lIOS of solid
NB.          part of diagonal
NB.   S    ≥ 0, the length of diagonal
NB.
NB. Formulae:
NB. - the whole diagonal's IO extreme element:
NB.     F := (d ≥ 0) ? d : (-n*d)
NB. - the whole diagonal's size:
NB.     S := max(0,min(m,n,⌊(n+m-|n-m-2*d|)/2⌋))
NB.
NB. Notes:
NB. - (f,s) pair defines raveled rIOS of solid part within diagonal

diaglios=: (0 0 _&$:) :(4 : 0)
  'd f s'=. x=. ((i. 3) < (# x)) } 0 0 _ ,: x  NB. in-place op
  'm n'=. y=. 2 $ y
  F=. n (- @ * ^: (0 > ])) d
  S=. 0 >. <./ y , <. -: (n + m - | n - m + +: d)
  (f ,: (s <. S)) (] ;. 0) hds2ios F , (>: n) , S
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Misc.

c=: 1{$              NB. Columns in matrix
trace=: +/ @ diag    NB. matrix trace
ct=: + @ |:          NB. conjugate transpose
pt=: |. @ |: @ |.    NB. pertranspose
cpt=: + @ pt         NB. conjugate pertranspose
pp=: [ C."1 C.       NB. apply permutation x to both rows and columns of table y
p2P=: =/ (i. @ #)    NB. transform permutation vector to permutation matrix
ip2P=: =/~ (i. @ #)  NB. transform inversed permutation vector to permutation
                     NB.   matrix, or permutation vector to inversed permutation matrix

NB. ---------------------------------------------------------
NB. rt
NB.
NB. Description:
NB.   Restrained Take. Just like built-in Take verb ({.), but
NB.   without overtake feature. Overtaking sizes mean "all
NB.   elements along this axis".
NB.
NB. Examples:
NB.    2 rt i. 3 4                  _2 rt i. 3 4
NB. 0 1 2 3                      4 5  6  7
NB. 4 5 6 7                      8 9 10 11
NB.    2 30 rt i. 3 4               2 _ rt i. 3 4
NB. 0 1 2 3                      0 1 2 3
NB. 4 5 6 7                      4 5 6 7
NB.    20 3 rt i. 3 4               _ 3 rt i. 3 4
NB. 0 1  2                       0 1  2
NB. 4 5  6                       4 5  6
NB. 8 9 10                       8 9 10
NB.    _2 _30 rt i. 3 4             _2 __ rt i. 3 4
NB. 4 5  6  7                    4 5  6  7
NB. 8 9 10 11                    8 9 10 11
NB.    _20 _3 rt i. 3 4             __ _3 rt i. 3 4
NB. 1  2  3                      1  2  3
NB. 5  6  7                      5  6  7
NB. 9 10 11                      9 10 11

rt=: ((({.~#)~$),:[)({.~(>|)/`(_,:{:)})~^:(((>+./@:*. _~:])|)/@[)]

NB. ---------------------------------------------------------
NB. upd1
NB.
NB. Description:
NB.   Adv. to update subarray by a monad
NB.
NB. Syntax:
NB.   vapp=. u upd1
NB. where
NB.   u       - monad to update subA; is called as:
NB.               subAupd=. u subA
NB.   vapp    - verb to update A; is called as:
NB.               Aupd=. ios vapp A
NB.   ios     - IOS of subA in the A
NB.   subA    - subarray in the A
NB.   A       - array
NB.   Aupd    - A with subA being replaced by subAupd
NB.
NB. If:
NB.   vapp=. u upd1
NB.   subA=. ios { A
NB.   subAupd=. u subA
NB.   Aupd=. subAupd ios } A
NB. then
NB.   Aupd -: ios vapp A
NB.
NB. Example:######################
NB. - the following replaces 87 by _8525 in array (i. 10 10):
NB.     (5 2 2 $ 4 9 1 1 3 5 1 1 2 4 1 1 1 _1 1 1 _2 _3 1 1) (+:`+`-:`-`*:`*`%: map4ri 0 1 2 3 4) i. 10 10
NB.   since
NB.     _8525 -: (+: 19) + (-: 24) - (*: 35) * (%: 49)
NB.
NB. References:
NB. [1] [Jprogramming] Transform to Amend
NB.     Dan Bron, Sat Mar 3 03:26:44 HKT 2007
NB.     http://www.jsoftware.com/pipermail/programming/2007-March/005425.html

upd1=: (@:{) (`[) (`]) }

NB. ---------------------------------------------------------
NB. append
NB.
NB. Description:
NB.   Template adv. to make verbs to enhance append built-in
NB.   verb (,)
NB.
NB. Examples:
NB.    (3 3$3) 0 append (2 2$2)     (3 3$3) _1 append (2 2$2)
NB. 3 3 3                        3 3 3
NB. 3 3 3                        3 3 3
NB. 3 3 3                        3 3 3
NB. 2 2 0                        0 2 2
NB. 2 2 0                        0 2 2
NB.    (2 2$2) 0 append (3 3$3)     (2 2$2) _1 append (3 3$3)
NB. 2 2 0                        0 2 2
NB. 2 2 0                        0 2 2
NB. 3 3 3                        3 3 3
NB. 3 3 3                        3 3 3
NB. 3 3 3                        3 3 3
NB.
NB. Notes:
NB. - 1-rank arrays (i.e. vectors) are also acceptable

append=: 1 : '(,`([,({."1~ ((-^:m)@c))~)`(({."1~ ((-^:m)@c)),]) @. ((*@-)&c)) & (,: ^: (2 > (#@$)))'

NB. ---------------------------------------------------------
NB. stitch
NB.
NB. Description:
NB.   Template adv. to make verbs to enhance stitch built-in
NB.   verb (,.)
NB.
NB. Examples:
NB.    (3 3$3) 0 stitch (2 2$2)     (3 3$3) _1 stitch (2 2$2)
NB. 3 3 3 2 2                    3 3 3 0 0
NB. 3 3 3 2 2                    3 3 3 2 2
NB. 3 3 3 0 0                    3 3 3 2 2
NB.    (2 2$2) 0 stitch (3 3$3)     (2 2$2) _1 stitch (3 3$3)
NB. 2 2 3 3 3                    0 0 3 3 3
NB. 2 2 3 3 3                    2 2 3 3 3
NB. 0 0 3 3 3                    2 2 3 3 3
NB.
NB. Notes:
NB. - 1-rank arrays (i.e. vectors) are also acceptable

stitch=: 1 : '((({.~ #),.])`([,.({.~ #)~)@.(>&#))`((({.~ (-@#)),.])`([,.({.~ (-@#))~)@.(>&#))@.(m"_)'

NB. ---------------------------------------------------------
NB. diag
NB.
NB. Description:
NB.   Return a solid part of diagonal of rectangular matrix
NB.
NB. Syntax:
NB.   e=. [(d[,f[,s]])] diag A
NB. where
NB.   A - m×n-matrix
NB.   d - integer in range [1-m,n-1], optional IO diagonal,
NB.       default is 0 (main diagonal)
NB.   f - integer in range [-S,S-1], optional IO extreme
NB.       element of solid part within diagonal, default is 0
NB.       (take from head)
NB.   s - integer in range [-S,S] or ±∞, optional size of 
NB.       solid part within diagonal, default is +∞ (all
NB.       elements in forward direction)
NB.   e - min(S,|s|)-vector, elements from the solid part of
NB.       diagonal
NB.   S ≥ 0, the length of diagonal

diag=: ((<0 1)&|:) :((diaglios $) ({,) ])

NB. ---------------------------------------------------------
NB. setdiag
NB.
NB. Description:
NB.   Assign scalar value to a solid part of diagonal
NB.
NB. Syntax:
NB.   Aupd=. (e[,d[,f[,s]]]) setdiag A
NB. where
NB.   A    - m×n-matrix to change
NB.   e    - scalar, the value to assign
NB.   d    - integer in range [1-m,n-1], optional IO
NB.          diagonal, default is 0 (main diagonal)
NB.   f    - integer in range [-S,S-1], optional IO extreme
NB.          element of solid part within diagonal, default
NB.          is 0 (take from head)
NB.   s    - integer in range [-S,S] or ±∞, optional size of
NB.          solid part within diagonal, default is +∞ (all
NB.          elements in forward direction)
NB.   Aupd - m×n-matrix A with value e assigned to solid part
NB.          within d-th diagonal
NB.   S    ≥ 0, the length of d-th diagonal

setdiag=: 4 : 0
  'e dfs'=. ({. ; }.) x=. ((i. 4) < (# x)) } 0 0 0 _ ,: x  NB. in-place op
  lios=. dfs (diaglios $) y
  e (lios"_) } y
)

NB. ---------------------------------------------------------
NB. upddiag
NB.
NB. Description:
NB.   Template adv. to make verbs to update a solid part of
NB.   diagonal
NB.
NB. Syntax:
NB.   vapp=. u upddiag
NB. where
NB.   u    - monad to change elements; is called as:
NB.            eupd=. u e
NB.   vapp - ambivalent verb to update a solid part within
NB.          diagonal of matrix A by monad u; is called
NB.          as:
NB.             Aupd=. [(d,[f[,s]])] vapp A
NB.   d    - integer in range [1-m,n-1], optional IO
NB.          diagonal, default is 0 (main diagonal)
NB.   f    - integer in range [-S,S-1], optional IO extreme
NB.          element of solid part within diagonal, default
NB.          is 0 (take from head)
NB.   s    - integer in range [-S,S] or ±∞, optional size of
NB.          solid part within diagonal, default is +∞ (all
NB.          elements in forward direction)
NB.   A    - m×n-matrix to update
NB.   Aupd - m×n-matrix A with solid part within d-th
NB.          diagonal being updated by monad u
NB.   S    ≥ 0, the length of d-th diagonal
NB.
NB. TODO:
NB. - [Jgeneral] duce/fold in J as an adverb or conjuction
NB.   Henry Rich, Sun Nov 25 14:07:46 HKT 2007
NB.   http://www.jsoftware.com/pipermail/general/2007-November/031233.html

upddiag=: 1 : 0
  lios=. diaglios $ y
  e=. lios ({,) y
  (u e) (lios"_) } y
:
  lios=. x diaglios $ y
  e=. lios ({,) y
  (u e) (lios"_) } y
)

NB. ---------------------------------------------------------
NB. idmat
NB.
NB. Description:
NB.   Make identity matrix with units on solid part of
NB.   diagonal
NB.
NB. Syntax:
NB.   I=. [(d[,f[,s]])] idmat [m,]n
NB. where
NB.   m    ≥ 0, integer, optional rows in matrix I, default
NB.          is n
NB.   n    ≥ 0, integer, columns in matrix I
NB.   d    - integer in range [1-m,n-1], optional IO
NB.          diagonal, default is 0 (main diagonal)
NB.   f    - integer in range [-S,S-1], optional IO extreme
NB.          element of solid part within diagonal, default
NB.          is 0 (take from head)
NB.   s    - integer in range [-S,S] or ±∞, optional size of
NB.          solid part within diagonal, default is +∞ (all
NB.          elements in forward direction)
NB.   I    - m×n-matrix of zeros with unit assigned to solid
NB.          part within d-th diagonal
NB.   S    ≥ 0, the length of d-th diagonal
NB.
NB. Examples:
NB.    idmat 3                      idmat 3 4
NB. 1 0 0                        1 0 0 0
NB. 0 1 0                        0 1 0 0
NB. 0 0 1                        0 0 1 0
NB.    1 idmat 3 4                  _1 idmat 3 4
NB. 0 1 0 0                      0 0 0 0
NB. 0 0 1 0                      1 0 0 0
NB. 0 0 0 1                      0 1 0 0

idmat=: (0 0 _&$:) :((1 , [) setdiag (0 $~ 2 $ ]))

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
NB.   h - integer in range [1-m,n-1], IO diagonal of v's
NB.       head, relatively to top left corner, default is 0
NB.   t - integer in range [1-m,n-1], IO diagonal of v's
NB.       tail, relatively to bottom right corner, default is
NB.       0
NB.   D - m×n-matrix of zeros with vector e assigned to h-th
NB.       diagonal
NB.   S ≥ 0, the length of h-th diagonal
NB.
NB. Examples:
NB.    diagmat 3 5 7                0 0 diagmat 3 5 7
NB. 3 0 0                        3 0 0
NB. 0 5 0                        0 5 0
NB. 0 0 7                        0 0 7
NB.    1 0 diagmat 3 5 7            _1 0 diagmat 3 5 7
NB. 0 3 0 0                      0 0 0
NB. 0 0 5 0                      3 0 0
NB. 0 0 0 7                      0 5 0
NB.                              0 0 7
NB.    0 1 diagmat 3 5 7            0 _1 diagmat 3 5 7
NB. 3 0 0                        3 0 0 0
NB. 0 5 0                        0 5 0 0
NB. 0 0 7                        0 0 7 0
NB. 0 0 0

diagmat=: (0 & $:) :(4 : 0)
  sh=. (#y) + (2&(|.@}. - {.)@(0&(<. , >.))) x  NB. find D shape
  lios=. ({. x) diaglios sh                     NB. lIOS for h-th diagonal
  y (lios"_) } sh $ 0                           NB. write e into matrix of zeros
)

NB. ---------------------------------------------------------
NB. trl
NB.
NB. Description:
NB.   Extract lower triangular (trapezoidal) matrix with
NB.   optional shrinking
NB.
NB. Examples:
NB.    trl >: i. 3 4                0 trl >: i. 3 4
NB. 1  0  0                      1  0  0
NB. 5  6  0                      5  6  0
NB. 9 10 11                      9 10 11
NB.    1 trl >: i. 3 4              _1 trl >: i. 3 4
NB. 1  2  0  0                   5  0
NB. 5  6  7  0                   9 10
NB. 9 10 11 12
NB.    1 trl >: i. 4 3              _1 trl >: i. 4 3
NB.  1  2  0                      4  0  0
NB.  4  5  6                      7  8  0
NB.  7  8  9                     10 11 12
NB. 10 11 12

trl=: (>:~ 0&>.) tr trlcut

NB. ---------------------------------------------------------
NB. tru
NB.
NB. Description:
NB.   Extract upper triangular (trapezoidal) matrix with
NB.   optional shrinking
NB.
NB. Examples:
NB.    tru >: i. 3 4                0 tru >: i. 3 4
NB. 1 2  3  4                    1 2  3  4
NB. 0 6  7  8                    0 6  7  8
NB. 0 0 11 12                    0 0 11 12
NB.    1 tru >: i. 3 4              _1 tru >: i. 3 4
NB. 2 3  4                       1  2  3  4
NB. 0 7  8                       5  6  7  8
NB. 0 0 12                       0 10 11 12
NB.    1 tru >: i. 4 3              _1 tru >: i. 4 3
NB. 2 3                          1 2  3
NB. 0 6                          4 5  6
NB.                              0 8  9
NB.                              0 0 12

tru=: (<:~ 0&<.) tr trucut

NB. ---------------------------------------------------------
NB. trl0
NB.
NB. Description:
NB.   Extract strictly lower triangular (trapezoidal) matrix
NB.   with optional shrinking
NB.
NB. Examples:
NB.    trl0 >: i. 4 3               0 trl0 >: i. 4 3
NB.  0  0  0                      0  0  0
NB.  4  0  0                      4  0  0
NB.  7  8  0                      7  8  0
NB. 10 11 12                     10 11 12
NB.    1 trl0 >: i. 4 3             _1 trl0 >: i. 4 3
NB.  1  0  0                      0  0 0
NB.  4  5  0                      7  0 0
NB.  7  8  9                     10 11 0
NB. 10 11 12
NB.    1 trl0 >: i. 3 4             _1 trl0 >: i. 3 4
NB. 1  0  0 0                    0 0
NB. 5  6  0 0                    9 0
NB. 9 10 11 0

trl0=: (>~ 0&>.) tr trlcut

NB. ---------------------------------------------------------
NB. tru0
NB.
NB. Description:
NB.   Extract strictly upper triangular (trapezoidal) matrix
NB.   with optional shrinking
NB.
NB. Examples:
NB.    tru0 >: i. 3 4               0 tru0 >: i. 3 4
NB. 0 2 3  4                     0 2 3  4
NB. 0 0 7  8                     0 0 7  8
NB. 0 0 0 12                     0 0 0 12
NB.    1 tru0 >: i. 3 4             _1 tru0 >: i. 3 4
NB. 0 3 4                        1 2  3  4
NB. 0 0 8                        0 6  7  8
NB. 0 0 0                        0 0 11 12
NB.    1 tru0 >: i. 4 3             _1 tru0 >: i. 4 3
NB. 0 3                          1 2 3
NB. 0 0                          0 5 6
NB.                              0 0 9
NB.                              0 0 0

tru0=: (<~ 0&<.) tr trucut

NB. ---------------------------------------------------------
NB. trl1
NB.
NB. Description:
NB.   Extract unit lower triangular (trapezoidal) matrix with
NB.   optional shrinking
NB.
NB. Examples:
NB.    trl1 >: i. 4 3               0 trl1 >: i. 4 3
NB.  1  0  0                      1  0  0
NB.  4  1  0                      4  1  0
NB.  7  8  1                      7  8  1
NB. 10 11 12                     10 11 12
NB.    1 trl1 >: i. 4 3             _1 trl1 >: i. 4 3
NB.  1  1  0                      1  0 0
NB.  4  5  1                      7  1 0
NB.  7  8  9                     10 11 1
NB. 10 11 12
NB.    1 trl1 >: i. 3 4             _1 trl1 >: i. 3 4
NB. 1  1  0 0                    1 0
NB. 5  6  1 0                    9 1
NB. 9 10 11 1

trl1=: (0&$:) :([ trl ((1 , [) setdiag ]))

NB. ---------------------------------------------------------
NB. tru1
NB.
NB. Description:
NB.   Extract unit upper triangular (trapezoidal) matrix with
NB.   optional shrinking
NB.
NB. Examples:
NB.    tru1 >: i. 3 4               0 tru1 >: i. 3 4
NB. 1 2 3  4                     1 2 3  4
NB. 0 1 7  8                     0 1 7  8
NB. 0 0 1 12                     0 0 1 12
NB.    1 tru1 >: i. 3 4             _1 tru1 >: i. 3 4
NB. 1 3 4                        1 2  3  4
NB. 0 1 8                        1 6  7  8
NB. 0 0 1                        0 1 11 12
NB.    1 tru1 >: i. 4 3             _1 tru1 >: i. 4 3
NB. 1 3                          1 2 3
NB. 0 1                          1 5 6
NB.                              0 1 9
NB.                              0 0 1

tru1=: (0&$:) :([ tru ((1 , [) setdiag ]))
