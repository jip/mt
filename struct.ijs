NB. struct.ijs
NB. Structure handlers
NB.
NB. IOSFC     Noun, IOS 1st column
NB. IOSLC     Noun, IOS last column
NB. IOSFR     Noun, IOS 1st row
NB. IOSLR     Noun, IOS last row
NB.
NB. trace     Matrix trace
NB. ct        Conjugate transpose
NB. pt        Pertranspose
NB. cpt       Conjugate pertranspose
NB. updl      Conj. to update array accepting linear IOS
NB. gi        Conj. to evoke n-th verb from gerund m
NB.
NB. uncut     To frame matrix by border of zeros
NB. append    Template adv. to make verbs to enhance append
NB.           verb
NB. stitch    Template adv. to make verbs to enhance stitch
NB.           verb
NB.
NB. diag      Return solid part of diagonal
NB. setdiag0  Assign scalar or vector value to all elements
NB.           of 0-th diagonal
NB. setdiag   Assign scalar value to solid part of diagonal
NB. upddiag0  Template adv. to make verbs to update all
NB.           elements of 0-th diagonal
NB. upddiag   Template adv. to make verbs to update solid
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
NB. diaglios
NB. Return lIOS of solid part of diagonal of rectangular
NB. matrix
NB.
NB. Syntax:
NB.   lios=. [d[,f,s]] diaglios [m,]n
NB. where
NB.   m    ≥ 0, integer, optional rows in matrix, default is
NB.          n
NB.   n    ≥ 0, integer, columns in matrix
NB.   d    - integer, IO diagonal, default is 0 (main
NB.          diagonal)
NB.   f    - integer in range [-min(m,n):min(m,n)-1], IO 1st
NB.          element of diagonal's fragment to be returned,
NB.          default is 0 (take from head)
NB.   s    - integer in range [-min(m,n):min(m,n)] or ±∞,
NB.          size of diagonal's fragment to be returned,
NB.          default is _ (all elements in forward direction)
NB.   lios - (0:min(m,n))-vector of integers, lIOS of
NB.          d-th diagonal elements in matrix of shape (m,n)
NB.
NB. Formulae:
NB. - the whole diagonal's IO 1st element:
NB.     F := (d ≥ 0) ? d : (-n*d)
NB. - the whole diagonal's size:
NB.     S := max(0,min(m,n,⌊(n+m-|n-m-2*d|)/2⌋))
NB.
NB. Notes:
NB. - (f,s) pair defines raveled rIOS of diagonal's fragment
NB.   to be returned

diaglios=: (0 0 _&$:) :(4 : 0)
  'd f s'=. 3 {. x , 0 _
  'm n'=. ({. , {:) y
  F=. n (-@* ^: (0 > ])) d
  S=. 0 >. <./ y , <. -: (n + m - | n - m + +: d)
  (f ,: (s <. S)) (] ;. 0) hds2ios F , (>: n) , S
)

NB. ---------------------------------------------------------
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

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Nouns, IOS a priori known submatrices

IOSFC=: < a: ; 0   NB. IOS 1st column
IOSLC=: < a: ; _1  NB. IOS last column
IOSFR=: 0          NB. IOS 1st row, or (< 0 ; < a:)
IOSLR=: _1         NB. IOS last row, or (< _1 ; < a:)

NB. ---------------------------------------------------------
NB. Misc.

trace=: +/ @ diag                        NB. matrix trace
ct=: + @ |:                              NB. conjugate transpose
pt=: |. @ |: @ |.                        NB. pertranspose
cpt=: + @ pt                             NB. conjugate pertranspose
updl=: 2 : '((n"_)})~ (u @ (n & ({,)))'  NB. Conj. to update array accepting linear IOS
gi=: 2 : '(n{m)`:6'                      NB. Conj. to evoke n-th verb from gerund m: m[n]

NB. ---------------------------------------------------------
NB. uncut
NB. To frame matrix by border of zeros
NB.
NB. Syntax:
NB.   A=. ((t,l),:(b,r)) uncut subA
NB. where
NB.   t l b r ≥ 0, integers, border widths at top, left,
NB.             bottom and right, respectively
NB.   subA    - h×w-matrix to border
NB.   A       - (t+h+b)×(l+w+r)-matrix such that
NB.               subA -: ((t,h),:(l,w)) (] ;. 0) A
NB.             and other elements are zeros
NB.
NB. Note:
NB. - this verb is some kind of the cut inversion

uncut=: ((((_1 1 * (+"1)) (+/\))~ $)) (({:@[) {. (({.~ {.)~)) ]

NB. ---------------------------------------------------------
NB. append
NB. Template adv. to make verbs to enhance append verb
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

append=: 1 : ',`((({."_1~ (-@{:@$)),])`([,({."_1~ (-@{:@$))~)@.(>&({:@$)))@.(m"_)'

NB. ---------------------------------------------------------
NB. stitch
NB. Template adv. to make verbs to enhance stitch verb
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

stitch=: 1 : '((({.~ #),.])`([,.({.~ #)~)@.(>&#))`((({.~ (-@#)),.])`([,.({.~ (-@#))~)@.(>&#))@.(m"_)'

NB. ---------------------------------------------------------
NB. diag
NB. Ambivalent verb to return a solid part of diagonal
NB.
NB. Syntax:
NB.   e=. [d[,f,s]] diag A
NB. where
NB.   A - m×n-matrix
NB.   d - integer, optional IO diagonal, default is 0 (main
NB.       diagonal)
NB.   f - integer in range [-min(m,n):min(m,n)-1], optional
NB.       IO 1st element of diagonal's fragment to be
NB.       returned, default is 0 (take from head)
NB.   s - integer in range [-min(m,n):min(m,n)] or ±∞,
NB.       optional size of diagonal's fragment to be
NB.       returned, default is +∞ (all elements in forward
NB.       direction)
NB.   e - (0:min(m,n))-vector, d-th diagonal elements from
NB.       matrix A

diag=: ((<0 1)&|:) :((diaglios $) ({,) ])

NB. ---------------------------------------------------------
NB. setdiag0
NB. Assign scalar or vector value to all elements of 0-th
NB. diagonal
NB.
NB. Syntax:
NB.   Aupd=. e setdiag0 A
NB. where
NB.   A    - m×n-matrix to change
NB.   e    - scalar or min(m,n)-vector, new diagonal
NB.          element(s)
NB.   Aupd - m×n-matrix A with 0-th diagonal replaced by e
NB.
NB. TODO:
NB. - change solid part possibility, too

setdiag0=: (diaglios @ $ @ ]) }

NB. ---------------------------------------------------------
NB. setdiag
NB. Assign scalar value to a solid part of any diagonal
NB.
NB. Syntax:
NB.   Aupd=. (e[,d[,f[,s]]]) setdiag A
NB. where
NB.   A    - m×n-matrix to change
NB.   e    - scalar, diagonal elements new value
NB.   d    - integer in range [1-m,n-1], optional IO
NB.          diagonal to change, default is 0
NB.   f    - integer in range [-k,k-1], optional IO extreme
NB.          element of solid part within d-th diagonal to
NB.          change, default is 0
NB.   s    - integer in range [-k,k-1] or ±∞, optional size
NB.          of solid part within diagonal to change,
NB.          default is 0
NB.   Aupd - m×n-matrix A with d-th diagonal elements defined
NB.          by rIOS (f,:s), being replaced by value v
NB.   k    - integer, the length of d-th diagonal

setdiag=: 4 : 0
  'v d f s'=. x=. ((i. 4) < (# x)) } 0 0 0 _ ,: x
  lios=. (f ,: s) (] ;. 0) d (diaglios $) y
  v (lios"_) } y
)

NB. ---------------------------------------------------------
NB. upddiag0
NB. Template adv. to make verbs to update a solid part of
NB. 0-th diagonal
NB.
NB. Syntax:
NB.   vapp=. u upddiag0
NB. where
NB.   u    - monad to change elements; is called as:
NB.            eupd=. u e
NB.   vapp - ambivalent verb to update a solid part of
NB.          0-th diagonal of matrix A by monad u; is called
NB.          as:
NB.            Aupd=. [d[,f,s]] vapp A
NB.   d    - integer, optional IO diagonal, default is 0
NB.          (main diagonal)
NB.   f    - integer in range [-min(m,n):min(m,n)-1],
NB.          optional IO 1st element of diagonal's fragment
NB.          to be returned, default is 0 (take from head)
NB.   s    - integer in range [-min(m,n):min(m,n)] or ±∞,
NB.          optional size of diagonal's fragment to be
NB.          returned, default is +∞ (all elements in forward
NB.          direction)
NB.   A    - m×n-matrix to update
NB.   Aupd - m×n-matrix A with 0-th diagonal elements being
NB.          replaced by value(s) eupd
NB.
NB. Note:
NB. - may be used to map any diagonal to 0-th one

upddiag0=: 1 : '(u @ diag) setdiag0 ]'

NB. ---------------------------------------------------------
NB. upddiag
NB. Template adv. to make verbs to update a solid part of
NB. diagonal
NB.
NB. Syntax:
NB.   vapp=. u upddiag
NB. where
NB.   u    - monad to change elements; is called as:
NB.            eupd=. u e
NB.   vapp - ambivalent verb to update a solid part of
NB.          d-th diagonal of matrix A by monad u; is called
NB.          as:
NB.             Aupd=. [d,[f,s]] vapp A
NB.   d    - integer, optional IO diagonal, default is 0
NB.          (main diagonal)
NB.   f    - integer in range [-min(m,n):min(m,n)-1],
NB.          optional IO 1st element of diagonal's fragment
NB.          to be returned, default is 0 (take from head)
NB.   s    - integer in range [-min(m,n):min(m,n)] or ±∞,
NB.          optional size of diagonal's fragment to be
NB.          returned, default is +∞ (all elements in forward
NB.          direction)
NB.   A    - m×n-matrix to update
NB.   Aupd - m×n-matrix A with d-th diagonal elements defined
NB.          by rIOS (f,:s), being replaced by value(s) eupd

upddiag=: 1 : 0
  (u upddiag0) y
:
  lios=. x (diaglios $) y
  (u (lios ({,) y)) (lios"_) } y
)

NB. ---------------------------------------------------------
NB. idmat
NB. Make rectangular identity matrix with shifted diagonal,
NB. partially filled by units
NB.
NB. Syntax:
NB.   I=. [d[,f,s]] idmat [m,]n
NB.
NB. Examples:
NB.    idmat 3                        idmat 3 4
NB. 1 0 0                          1 0 0 0
NB. 0 1 0                          0 1 0 0
NB. 0 0 1                          0 0 1 0
NB.    1 idmat 3 4                    _1 idmat 3 4
NB. 0 1 0 0                        0 0 0 0
NB. 0 0 1 0                        1 0 0 0
NB. 0 0 0 1                        0 1 0 0

idmat=: (0 & $:) :((setdiag~ (1 & ,))~ (({. , {:) $ 0:))

NB. ---------------------------------------------------------
NB. diagmat
NB. Make rectangular diagonal matrix with y on diagonal
NB.
NB. Syntax:
NB.   D=. [se] diagmat v
NB. where
NB.   v  - min(m,n)-vector, new values for diagonal
NB.   se - complex number (s j. e):
NB.        s - integer, diagonal number for v's head,
NB.            relatively to top left corner
NB.        e - integer, diagonal number for v's tail,
NB.            relatively to bottom right corner
NB.        default se=0j0, i.e. square matrix 0-th diagonal
NB.   D  - m×n-matrix, diagonal rectangular with v on diagonal
NB.
NB. Examples:
NB.    diagmat 3 5 7                  0j0 diagmat 3 5 7
NB. 3 0 0                          3 0 0
NB. 0 5 0                          0 5 0
NB. 0 0 7                          0 0 7
NB.    1j0 diagmat 3 5 7              _1j0 diagmat 3 5 7
NB. 0 3 0 0                        0 0 0
NB. 0 0 5 0                        3 0 0
NB. 0 0 0 7                        0 5 0
NB.                                0 0 7
NB.    0j1 diagmat 3 5 7              0j_1 diagmat 3 5 7
NB. 3 0 0                          3 0 0 0
NB. 0 5 0                          0 5 0 0
NB. 0 0 7                          0 0 7 0
NB. 0 0 0

diagmat=: (0 & $:) :(4 : 0)
  'm r'=. sh=. (#y) + (2&(|.@}. - {.)@(0&(<. , >.)@+.)) x  NB. find D shape
  d=. 9 o. x                                               NB. IO diagonal
  lios=. d diaglios sh                                     NB. lIOS for d
  y (lios " _) } sh $ 0                                    NB. write v into matrix of zeros
)

NB. ---------------------------------------------------------
NB. trl
NB. Extract lower triangular (trapezoidal) matrix with
NB. optional shrinking
NB.
NB. Examples:
NB.    trl i. 3 4                     0 trl i. 3 4
NB. 0 0  0                         0 0  0
NB. 4 5  0                         4 5  0
NB. 8 9 10                         8 9 10
NB.    1 trl i. 3 4                   _1 trl i. 3 4
NB. 0 1  0  0                      4 0
NB. 4 5  6  0                      8 9
NB. 8 9 10 11
NB.    1 trl i. 4 3                   _1 trl i. 4 3
NB. 0  1  0                        3  0  0
NB. 3  4  5                        6  7  0
NB. 6  7  8                        9 10 11
NB. 9 10 11

trl=: (>:~ 0&>.) tr trlcut

NB. ---------------------------------------------------------
NB. tru
NB. Extract upper triangular (trapezoidal) matrix with
NB. optional shrinking
NB.
NB. Examples:
NB.    tru i. 3 4                     0 tru i. 3 4
NB. 0 1  2  3                      0 1  2  3
NB. 0 5  6  7                      0 5  6  7
NB. 0 0 10 11                      0 0 10 11
NB.    1 tru i. 3 4                   _1 tru i. 3 4
NB. 1 2  3                         0 1  2  3
NB. 0 6  7                         4 5  6  7
NB. 0 0 11                         0 9 10 11
NB.    1 tru i. 4 3                   _1 tru i. 4 3
NB. 1 2                            0 1  2
NB. 0 5                            3 4  5
NB.                                0 7  8
NB.                                0 0 11

tru=: (<:~ 0&<.) tr trucut

NB. ---------------------------------------------------------
NB. trl0
NB. Extract strictly lower triangular (trapezoidal) matrix
NB. with optional shrinking
NB.
NB. Examples:
NB.    trl0 i. 4 3                0 trl0 i. 4 3
NB. 0  0  0                        0  0  0
NB. 3  0  0                        3  0  0
NB. 6  7  0                        6  7  0
NB. 9 10 11                        9 10 11
NB.    1 trl0 i. 4 3                  _1 trl0 i. 4 3
NB. 0  0  0                        0  0 0
NB. 3  4  0                        6  0 0
NB. 6  7  8                        9 10 0
NB. 9 10 11
NB.    1 trl0 i. 3 4                  _1 trl0 i. 3 4
NB. 0 0  0 0                       0 0
NB. 4 5  0 0                       8 0
NB. 8 9 10 0

trl0=: (>~ 0&>.) tr trlcut

NB. ---------------------------------------------------------
NB. tru0
NB. Extract strictly upper triangular (trapezoidal) matrix
NB. with optional shrinking
NB.
NB. Examples:
NB.    tru0 i. 3 4                    0 tru0 i. 3 4
NB. 0 1 2  3                       0 1 2  3
NB. 0 0 6  7                       0 0 6  7
NB. 0 0 0 11                       0 0 0 11
NB.    1 tru0 i. 3 4                  _1 tru0 i. 3 4
NB. 0 2 3                          0 1  2  3
NB. 0 0 7                          0 5  6  7
NB. 0 0 0                          0 0 10 11
NB.    1 tru0 i. 4 3                  _1 tru0 i. 4 3
NB. 0 2                            0 1 2
NB. 0 0                            0 4 5
NB.                                0 0 8
NB.                                0 0 0

tru0=: (<~ 0&<.) tr trucut

NB. ---------------------------------------------------------
NB. trl1
NB. Extract unit lower triangular (trapezoidal) matrix with
NB. optional shrinking

trl1=: (trl @ (1 & setdiag0)) : ([ trl ((1 , [) setdiag ]))

NB. ---------------------------------------------------------
NB. tru1
NB. Extract unit upper triangular (trapezoidal) matrix with
NB. optional shrinking

tru1=: (tru @ (1 & setdiag0)) : ([ tru ((1 , [) setdiag ]))
