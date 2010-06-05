NB. Structure handlers
NB.
NB. c         Columns in scalar or array
NB. trace     Matrix trace
NB. ct        Conjugate transpose
NB. cp        Conjugate pertranspose
NB. sp        Symmetric permutation
NB. p2P       Transform permutation vector to permutation
NB.           matrix
NB. ip2P      Transform inversed permutation vector to
NB.           permutation matrix
NB. rt        Restrained Take
NB. icut      Inversed cut
NB.
NB. upd1      Adv. to update subarray by a monad
NB. append    Template adv. to make verbs to enhance append
NB.           built-in verb (,)
NB. stitch    Template adv. to make verbs to enhance stitch
NB.           built-in verb (,.)
NB.
NB. diag      Return a solid part of diagonal
NB. setdiag   Assign value[s] to a solid part of diagonal
NB. upddiag   Template adv. to make verbs to update a solid
NB.           part of diagonal
NB.
NB. bdlpick   Zeroize elements located outside lower
NB.           bidiagonal part of the matrix
NB. bdupick   Zeroize elements located outside upper
NB.           bidiagonal part of the matrix
NB. hslpick   Zeroize elements located outside lower
NB.           Hessenberg part of the matrix
NB. hsupick   Zeroize elements located outside lower
NB.           Hessenberg part of the matrix
NB. tdpick    Zeroize elements located outside tridiagonal
NB.           part of the matrix
NB. trlpick   Zeroize elements located outside lower
NB.           triangular part of the matrix
NB. trupick   Zeroize elements located outside upper
NB.           triangular part of the matrix
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
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Misc.

NB. convert table y to table of diagonals
t2td=: 1 : '({. u/&i. {:)@$'

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
tr=: 2 : '0&$: : ([ (] * (u~ (-~ t2td))) v)'

NB. ---------------------------------------------------------
NB. mxbstencil
NB.
NB. Description:
NB.   Template adv. to make verbs returning
NB.   [multi-][anti-]band stencil for rectangular matrix
NB.
NB. Syntax:
NB.   vapp=. vmix mxbstencil
NB. where
NB.   vmix - dyad to mix lIOS x and y, is either (-~) for
NB.          band, or (+) for anti-band stencils, is evoked
NB.          as:
NB.            mix=. lIOrow vmix lIOcolumn
NB.   vapp - dyad to make multi-[anti-]band stencil, is
NB.          evoked as:
NB.            s=. bs vapp A
NB.   bs   - k×2-matrix of (b)s, or single b, or d, defines
NB.          [anti-]bands to stencil
NB.   b    - 2-vector (h,t), defines one [anti-]band to
NB.          stencil
NB.   h    - integer in range [-∞,t], lIO head
NB.          [anti-]diagonal
NB.   t    - integer in range [h,+∞], lIO tail
NB.          [anti-]diagonal
NB.   d    - integer in range [-∞,+∞], lIO one
NB.          [anti-]diagonal to stencil
NB.   A    - m×n-matrix
NB.   s    - m×n-matrix, boolean, having 1s on [anti-]band[s]
NB.
NB. Examples:
NB. - see mbstencil, mabstencil

mxbstencil=: 1 : '(+./^:(_2+(#@$))@:((1=I.)"1 2)~ (-&1 0"1))~ (u t2td)'

NB. ---------------------------------------------------------
NB. mbstencil
NB. mabstencil
NB.
NB. Description:
NB.   [Multi-]band and [multi-]anti-band stencils for
NB.   rectangular matrix
NB.
NB. Syntax:
NB.   s=. bs mbstencil  A
NB.   s=. bs mabstencil A
NB. where
NB.   bs   - k×2-matrix of (b)s, or single b, or d, defines
NB.          [anti-]bands to stencil
NB.   A    - m×n-matrix
NB.   s    - m×n-matrix, boolean, having 1s on [anti-]band[s]
NB.   b    - 2-vector (h,t), defines one [anti-]band to
NB.          stencil
NB.   h    - integer in range [-∞,t], defines lIO head
NB.          [anti-]diagonal
NB.   t    - integer in range [h,+∞], defines lIO tail
NB.          [anti-]diagonal
NB.   d    - integer in range [-∞,+∞], defines one
NB.          [anti-]diagonal to stencil
NB.
NB. Examples:
NB.    2 mbstencil i. 3 5                    2 mabstencil i. 3 5
NB. 0 0 1 0 0                             0 0 1 0 0
NB. 0 0 0 1 0                             0 1 0 0 0
NB. 0 0 0 0 1                             1 0 0 0 0
NB.    2 3 mbstencil i. 3 5                  2 3 mabstencil i. 3 5
NB. 0 0 1 1 0                             0 1 1 0 0
NB. 0 0 0 1 1                             1 1 0 0 0
NB. 0 0 0 0 1                             1 0 0 0 0
NB.    (__ _1 ,: 2 3) mbstencil i. 3 5       (__ _1 ,: 2 3) mabstencil i. 3 5
NB. 0 0 1 1 0                             0 1 1 0 0
NB. 1 0 0 1 1                             1 1 0 0 1
NB. 1 1 0 0 1                             1 0 0 1 1

mbstencil=:                       -~ mxbstencil
mabstencil=: ((|."1)@:-~ (<:@c)) (+  mxbstencil) ]

NB. ---------------------------------------------------------
NB. diaglios
NB.
NB. Description:
NB.   Return lIOS of solid part of diagonal of rectangular
NB.   matrix
NB.
NB. Syntax:
NB.   lios=. [(d[,h[,s]])] diaglios [m,]n
NB. where
NB.   m    ≥ 0, integer, optional rows in matrix, default is
NB.          n
NB.   n    ≥ 0, integer, columns in matrix
NB.   d    - integer in range [1-m,n-1], optional IO
NB.          diagonal, default is 0 (main diagonal)
NB.   h    - integer in range [-S,S-1], optional IO extreme
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
NB.     H := (d ≥ 0) ? d : (-n*d)
NB. - the whole diagonal's size:
NB.     S := max(0,min(m,n,⌊(n+m-|n-m-2*d|)/2⌋))
NB.
NB. Notes:
NB. - (h,s) pair defines raveled rIOS of solid part within
NB.   diagonal

diaglios=: (0 0 _&$:) :(4 : 0)
  'd h s'=. x=. ((i. 3) < (# x)) } 0 0 _ ,: x  NB. in-place op
  'm n'=. y=. 2 $ y
  H=. n (- @ * ^: (0 > ])) d
  S=. 0 >. <./ y , <. -: (n + m - | n - m + +: d)
  (h ,: (s <. S)) (] ;. 0) (>: n) dhs2lios H , S
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Misc.

c=: {:!.1 @ $      NB. Columns in scalar or array

trace=: +/ @ diag  NB. Matrix trace

ct=: + @ |:        NB. Conjugate transpose
cp=: |. @ ct @ |.  NB. Conjugate pertranspose

sp=: [ C."1 C.     NB. Symmetric permutation

p2P=:  {    =      NB. Transform permutation vector to
                   NB.   permutation matrix
ip2P=: {^:_1=      NB. Transform inversed permutation vector
                   NB.   to permutation matrix, or
                   NB.   permutation vector to inversed
                   NB.   permutation matrix

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
NB. Assertion:
NB.   A -: icut fret <;.1 A
NB. where
NB.   A    - some array
NB.   fret - some fret
NB.
NB. References:
NB. [1] JWiki/Essays/ Block Matrix Inverse
NB.     Roger Hui
NB.     http://www.jsoftware.com/jwiki/Essays/Block%20Matrix%20Inverse

icut=: [: > 3 : ',"(#$y)&.>/y'^:(#@$)

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

rt=: (*@[ * |@[ <. (({.~ #)~ $)) {. ]

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

NB. USEME:
appendl=: , `([, (({."1~    c )~))`(({."1~    c ), ]) @. (*@-&c)
appendr=: , `([, (({."1~ (-@c))~))`(({."1~ (-@c)), ]) @. (*@-&c)

stitcht=: ,.`([,.(({.  ~    # )~))`(({.  ~    # ),.]) @. (*@-&#)
stitchb=: ,.`([,.(({.  ~ (-@#))~))`(({.  ~ (-@#)),.]) @. (*@-&#)


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
NB.   e=. [(d[,h[,s]])] diag A
NB. where
NB.   A - m×n-matrix
NB.   d - integer in range [1-m,n-1], optional IO diagonal,
NB.       default is 0 (main diagonal)
NB.   h - integer in range [-S,S-1], optional IO extreme
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
NB.   Assign value[s] to a solid part of diagonal
NB.
NB. Syntax:
NB.   Aupd=. (e;[d[,h[,s]]]) setdiag A
NB. where
NB.   A    - m×n-matrix to change
NB.   e    - scalar or k-vector, value[s] to assign
NB.   d    - integer in range [1-m,n-1], optional IO
NB.          diagonal, default is 0 (main diagonal)
NB.   h    - integer in range [-S,S-1], optional IO extreme
NB.          element of solid part within diagonal, default
NB.          is 0 (take from head)
NB.   s    - integer in range [-S,S] or ±∞ when e is scalar,
NB.          or any from set {±k,±∞} when e is vector;
NB.          optional size of solid part within diagonal,
NB.          default is +∞ (all elements in forward
NB.          direction)
NB.   Aupd - m×n-matrix A with value[s] e assigned to solid
NB.          part within d-th diagonal
NB.   S    ≥ 0, the length of d-th diagonal
NB.   k    ≤ S, the length of vector e
NB.
NB. Examples:
NB. 
NB.    (2; a:) setdiag 4 4 $ 0         (2; _1 1 1) setdiag 4 4 $ 0
NB. 2 0 0 0                         0 0 0 0
NB. 0 2 0 0                         0 0 0 0
NB. 0 0 2 0                         0 2 0 0
NB. 0 0 0 2                         0 0 0 0
NB.    (2; _1) setdiag 4 4 $ 0         (1 2 3; _1) setdiag 4 4 $ 0
NB. 0 0 0 0                         0 0 0 0
NB. 2 0 0 0                         1 0 0 0
NB. 0 2 0 0                         0 2 0 0
NB. 0 0 2 0                         0 0 3 0
NB.    (2; _1 1) setdiag 4 4 $ 0       (1 2 3; _1 _1 _3) setdiag 4 4 $ 0
NB. 0 0 0 0                         0 0 0 0
NB. 0 0 0 0                         3 0 0 0
NB. 0 2 0 0                         0 2 0 0
NB. 0 0 2 0                         0 0 1 0

setdiag=: 4 : 0
  'e dhs'=. x
  dhs=. ((i. 3) < (# dhs)) } 0 0 _ ,: dhs  NB. assign defaults, in-place op
  lios=. dhs diaglios $ y
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
NB.             Aupd=. [(d,[h[,s]])] vapp A
NB.   d    - integer in range [1-m,n-1], optional IO
NB.          diagonal, default is 0 (main diagonal)
NB.   h    - integer in range [-S,S-1], optional IO extreme
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
  lios=. diaglios_mt_ $ y
  e=. lios ({,) y
  (u e) (lios"_) } y
:
  lios=. x diaglios_mt_ $ y
  e=. lios ({,) y
  (u e) (lios"_) } y
)

NB. ---------------------------------------------------------
NB. bdlpick
NB.
NB. Description:
NB.   Zeroize elements located outside lower bidiagonal part
NB.   of the matrix
NB.
NB. Syntax:
NB.   B=. bdlpick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   B - m×n-matrix, lower bidiagonal
NB.
NB. TODO:
NB. - B would be sparse

bdlpick=: * _1 0 & mbstencil

NB. ---------------------------------------------------------
NB. bdupick
NB.
NB. Description:
NB.   Zeroize elements located outside upper bidiagonal part
NB.   of the matrix
NB.
NB. Syntax:
NB.   B=. bdupick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   B - m×n-matrix, upper bidiagonal
NB.
NB. TODO:
NB. - B would be sparse

bdupick=: * 0 1 & mbstencil

NB. ---------------------------------------------------------
NB. hslpick
NB.
NB. Description:
NB.   Zeroize elements located outside lower Hessenberg part
NB.   of the matrix
NB.
NB. Syntax:
NB.   B=. hslpick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   B - m×n-matrix, lower Hessenberg

hslpick=: * __ 1 & mbstencil

NB. ---------------------------------------------------------
NB. hsupick
NB.
NB. Description:
NB.   Zeroize elements located outside upper Hessenberg part
NB.   of the matrix
NB.
NB. Syntax:
NB.   B=. hsupick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   B - m×n-matrix, upper Hessenberg

hsupick=: * _1 _ & mbstencil

NB. ---------------------------------------------------------
NB. tdpick
NB.
NB. Description:
NB.   Zeroize elements located outside tridiagonal part of
NB.   the matrix
NB.
NB. Syntax:
NB.   B=. tdpick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   B - m×n-matrix, tridiagonal
NB.
NB. TODO:
NB. - B would be sparse

tdpick=: * _1 1 & mbstencil

NB. ---------------------------------------------------------
NB. trlpick
NB.
NB. Description:
NB.   Zeroize elements located outside lower triangular part
NB.   of the matrix
NB.
NB. Syntax:
NB.   B=. [d] trlpick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   d - integer in range [-∞,+∞], lIO last non-zero
NB.       diagonal
NB.   B - m×n-matrix, lower triangular

trlpick=: (0&$:) : (((__ , [) mbstencil ]) * ])

NB. ---------------------------------------------------------
NB. trupick
NB.
NB. Description:
NB.   Zeroize elements located outside upper triangular part
NB.   of the matrix
NB.
NB. Syntax:
NB.   B=. [d] trupick A
NB. where
NB.   A - m×n-matrix, contains B
NB.   d - integer in range [-∞,+∞], lIO first non-zero
NB.       diagonal
NB.   B - m×n-matrix, upper triangular

trupick=: (0&$:) : (((_ ,~ [) mbstencil ]) * ])

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
NB.   m    ≥ 0, integer, optional rows in matrix I, default
NB.          is n
NB.   n    ≥ 0, integer, columns in matrix I
NB.   d    - integer in range [1-m,n-1], optional IO
NB.          diagonal, default is 0 (main diagonal)
NB.   h    - integer in range [-S,S-1], optional IO extreme
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
NB.
NB. TODO:
NB. - I would be sparse

idmat=: (a:&$:) :((1;[) setdiag (0 $~ 2 $ ]))

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
NB. Algorithm:
NB.   1) find D shape
NB.   2) generate lIOS for h-th diagonal
NB.   3) write e into matrix of zeros
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
NB.
NB. TODO:
NB. - D would be sparse

diagmat=: (0 0&$:) :((; {.)~ setdiag ((0 $~ ((+ (2&(|.@}. - {.)@(0&(<. , >.))))))~ #))

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

trl1=: (0&$:) :([ trl ((1;[) setdiag ]))

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

tru1=: (0&$:) :([ tru ((1;[) setdiag ]))
