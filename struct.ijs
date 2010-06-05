NB. struct.ijs
NB. Structure handlers
NB.
NB. ioscv     Locate constant rows or columns in matrix
NB. ct1       Count trailing 1s in Boolean vector
NB. ct0       Count trailing 0s in vector
NB. ct0c      Count trailing zero columns in matrix
NB. ct0r      Count trailing zero rows in matrix
NB. trace     Matrix trace
NB. ct        Conjugate transpose
NB.
NB. fromc     Model direct 'from' accepting cIOS
NB. fromi     Adv. to model indirect 'from'
NB. fromci    Adv. to model 'fromi' accepting cIOS
NB. amendci   Adv. to model indirect 'amend' accepting cIOS
NB. upd       Adv. to update subarray by monad
NB. updl      Conj. to model 'upd' accepting lIOS
NB. updi      Conj. to model indirect 'upd'
NB. updci     Conj. to model 'updi' accepting cIOS
NB. upd2i     Conj. to indirect update subarray by dyad
NB. upd2ci    Conj. to model 'upd2i' accepting cIOS
NB. upd3i     Conj. to indirect update subarray by triad
NB. upd3ci    Conj. to model 'upd3i' accepting cIOS
NB. mapi      Conj. to indirect map subarray to another one
NB.           by monad
NB. mapci     Conj. to model 'mapi' accepting cIOS
NB. mapli     Conj. to model 'mapi' accepting lIOS
NB. map2i     Conj. to indirect map two subarrays to another
NB.           one by dyad
NB. map2ci    Conj. to model 'map2i' accepting cIOS
NB.
NB. append    Template adv. to make verbs to enhance append
NB. stitch    Template adv. to make verbs to enhance stitch
NB.
NB. rws       Rows weighted sum
NB. sdiag     Add element[s from] x to diagonal of matrix y
NB.
NB. idmat     Make rectangular identity matrix with shifted
NB.           diagonal
NB. diagmat   Make rectangular diagonal matrix with shifted
NB.           diagonal
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
NB. Indirect extractors

NB. Adv. to extract m-th from x: x[m]
mx=: 1 : 'm { ['

NB. Conj. to extract n-th from m-th from x: x[m[n]]
nmx=: 2 : '(n { m) { ['

NB. Conj. to evoke n-th verb from gerund m: m[n]
gi=: 2 : '(n{m)`:6'

NB. Conj. to apply u to submatrix from y, defined by n-th
NB. rIOS from x: u(y[x[n]])
uci=: 2 : '(n mx) (u ;. 0) ]'

NB. Conj. to supply n-th verb from gerund to 'cut':
NB. u[n](y[x])
gic=: 2 : '(m gi n) ;. 0'

NB. Conj. to implement indirect version of gic:
NB. u[n[1]](y[x[n[0]]])
gici=: 2 : '((0{n) mx) (m gic (1{n)) ]'

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

IOSFC=: < a: ; 0     NB. IOS 1st column
IOSLC=: < a: ; _1    NB. IOS last column
IOSFR=: < 0 ; < a:   NB. IOS 1st row
IOSLR=: < _1 ; < a:  NB. IOS last row

NB. ---------------------------------------------------------
NB. Misc.

ioscv=: *./ .=            NB. locate constant vector in matrix (column: 'm & ioscv', row: 'ioscv & n')
ct1=: +/ @: (*./\.)       NB. count trailing 1s in Boolean vector y
ct0=: ct1 @ (0&=)         NB. count trailing 0s in vector y
ct0c=: ct1 @ (0 & ioscv)  NB. count trailing zero columns in matrix y
ct0r=: ct1 @ (ioscv & 0)  NB. count trailing zero rows in matrix y

trace=: +/ @ diag         NB. matrix trace
ct=: + @ |:               NB. conjugate transpose

NB. ---------------------------------------------------------
NB. IOS explorers
NB. monad     dyad    triad   tetrad     pentad     hexad     heptad    octad    ennead decade
NB.           couplet triplet quadruplet quintuplet sextuplet septuplet octuplet
NB. singleton duet    trio    quartet    quintet    sextet    septet    octet

NB. Adv. to update submatrix
NB. Syntax: Aupd=. ios (vapp upd) A
upd=: 1 : '(u@{)`[`]}'

NB. Conj. to model 'upd' accepting lIOS
updl=: 2 : '((n"_)})~ (u @ (n & ({,)))'

NB. Conj. to model 'upd' accepting rIOS
NB. y[x[n]] := u(y[x[n]])
NB. Aupd=. rios (u updr) A
updr=: 1 : '(u ;. 0)`(rios2ios@[)`] }'

NB. Conj. to model 'updi' accepting rIOS
NB. y[x[n]] := u(y[x[n]])
NB. Aupd=. rios (u updri io) A
updri=: 2 : '(n mx) (u updr) ]'

NB. Conj. to model 'upd2i' accepting rIOS
NB. y[x[n[1]]] := u[1](u[2](y[x[n[1]]]),u[0](y[x[n[0]]])))
upd2ri=: 2 : '((m gici ((1{n),2)) (m gi 1) (m gici ((0{n),0)))`(rios2ios@(n nmx 1))`] }'

NB. ---------------------------------------------------------
NB. upd3ri
NB. Template conj. to make verbs to update subarray by pentad
NB. indirectly, addressing mode is rIOS
NB.
NB. Syntax:
NB.   vapp=. u0`u1`u2`u3`u4 upd3ri io0 io1 io2
NB. where
NB.   rios     - k×2×2-array of integers, vector of rIOSs
NB.   io0..io2 - integers, IOs in rios
NB.   A        - input matrix
NB.   A0..A2   - 2-rank arrays, submatrices in A:
NB.                A0 -: (io0{rios) (] ;. 0) A
NB.                A1 -: (io1{rios) (] ;. 0) A
NB.                A2 -: (io2{rios) (] ;. 0) A
NB.   u0 u2 u4 - monads to adjust A2 A1 A0, respectively
NB.   u1 u3    - dyads to glue staff
NB.   vapp     - verb to update A2, is called as:
NB.                Aupd=. rios vapp A
NB.   Aupd     - A with submatrix A2 replaced by value:
NB.               (u0 A2) u1 (u2 A1) u3 (u4 A0)
NB.   k        > 0
NB.
NB. Example:
NB. - the following replaces 24 by _1220 in array (i. 10 10):
NB.     (3 2 2 $ 4 9 1 1 3 5 1 1 2 4 1 1) (-:`-`*:`+`%: upd3ri 0 1 2) i. 10 10
NB.   since
NB.     _1220 -: (-: 24) - (*: 35) + (%: 49)

upd3ri=: 2 : '((m gici ((2{n),0)) (m gi 1) (m gici ((1{n),2)) (m gi 3) (m gici ((0{n),4)))`(rios2ios@(n nmx 2))`] }'

NB. Conj. to 'map' set of elements to another one, accepting IOS
NB. y[x[1]] := u(y[x[0]])
NB. Aupd=. IOS (u map (ioFrom,ioTo)) A

map=: 2 : '(u @ (n nmx 0))`(n nmx 1)`] }'

NB. Adv. to model 'map' accepting lIOS
NB. y[u[1](u[0](y[x]),y)] := u[0](y[x])
NB. Aupd=. lios (u0`u1 mapl) A
mapl=: 1 : '((m gi 0) @ ({,)) ((m gi 1) }) ]'

NB. FIXME! NB. Adv. to model 'mapl' with indirect access to lIOS
NB. FIXME! NB. y[u[1](x[1])] := u[0](y[x[0]])
NB. FIXME! NB. vapp=. u0`u1 mapli
NB. FIXME! NB. Aupd=. liosa vapp A
NB. FIXME! mapli=: 1 : '((((m gi 0) @ ({,))~ (0&{))~) (((m gi 1) @ (1 mx)) }) ]'

NB. Adv. to model 'map' accepting rIOS to read and
NB. using lIOS to write
NB. y[u[1](u[0](y[x]),y)] := u[0](y[x])
NB. Aupd=. rIOS (u0`u1 maprl) A
maprl=: 1 : '(m gic 0) ((m gi 1) }) ]'

NB. Conj. to model 'maprl' with indirect access to rIOS
NB. y[u[1](u[0](y[x[n]]),y)] := u[0](y[x[n]])
NB. Aupd=. rIOSs (u0`u1 maprli ioFrom) A
maprli=: 2 : '(m gici (n,0)) ((m gi 1) }) ]'

NB. Conj. to model 'mapi' accepting rIOS
NB. y[u[1](y[x[n]],y)] := u[0](y[x[n]])
NB. Aupd=. rios (u mapri (ioFrom,ioTo)) A
mapri=: 2 : '(u uci (0{n))`(rios2ios@(n nmx 1))`] }'

NB. Conj. to model 'upd3i' accepting rIOS
NB. y[x[n[2]]] := u[3](u[4](y[x[n[2]]]),u[1](u[2](y[x[n[1]]]),u[0](y[x[n[0]]])))
map3ri=: 2 : '((m gici ((2{n),4)) (m gi 3) (m gici ((1{n),2)) (m gi 1) (m gici ((0{n),0)))`(rios2ios@(n nmx 3))`] }'

NB. ---------------------------------------------------------
NB. map4ri
NB. Template conj. to make verbs to map 4 subarrays by heptad
NB. to another subarray indirectly, addressing mode is rIOS
NB.
NB. Syntax:
NB.   vapp=. u0`u1`u2`u3`u4`u5`u6 map4ri io0 io1 io2 io3 io4 io5
NB. where
NB.   rios        - k×2×2-array of integers, vector of rIOSs
NB.   io0..io5    - integers, IOs in rios
NB.   A           - input matrix
NB.   A0..A4      - 2-rank arrays, submatrices in A:
NB.                   A0 -: (io0{rios) (] ;. 0) A
NB.                   A1 -: (io1{rios) (] ;. 0) A
NB.                   A2 -: (io2{rios) (] ;. 0) A
NB.                   A3 -: (io3{rios) (] ;. 0) A
NB.                   A4 -: (io4{rios) (] ;. 0) A
NB.   u0 u2 u4 u6 - monads to adjust A3 A2 A1 A0, respectively
NB.   u1 u3 u5    - dyads to glue staff
NB.   vapp        - verb to update A4, is called as:
NB.                   Aupd=. rios vapp A
NB.   Aupd        - A with submatrix A4 replaced by value:
NB.                  (u6 A3) u5 (u4 A2) u3 (u2 A1) u1 (u0 A0)
NB.   k           > 0
NB.
NB. Example:
NB. - the following replaces 87 by _8525 in array (i. 10 10):
NB.     (5 2 2 $ 4 9 1 1 3 5 1 1 2 4 1 1 1 _1 1 1 _2 _3 1 1) (+:`+`-:`-`*:`*`%: map4ri 0 1 2 3 4) i. 10 10
NB.   since
NB.     _8525 -: (+: 19) + (-: 24) - (*: 35) * (%: 49)

map4ri=: 2 : '((m gici ((3{n),0)) (m gi 1) ((m gici ((2{n),2)) (m gi 3) (m gici ((1{n),4)) (m gi 5) (m gici ((0{n),6))))`(rios2ios@(n nmx 4))`] }'

NB. ---------------------------------------------------------
NB. append
NB. Template adv. to make verbs to enhance append
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
NB. Template adv. to make verbs to enhance stitch
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
NB. rws
NB. Rows weighted sum
NB.
NB. Syntax:
NB.   sum=. weight rows rws array
NB.
NB. Application:
NB. - to subtract from row 3 twiced row 4:
NB.   ((1 _2) 3 4 rws (i. 5 5)

rws=: 1 : '+/ @ (* (m & {))'

NB. ---------------------------------------------------------
NB. diaglios
NB. Specified diagonal element's lIOS of rectangular matrix
NB.
NB. Syntax:
NB.   lios=. [d[,f,s]] diaglios (m,n)
NB. where
NB.   m    ≥ 0, integer, rows in matrix
NB.   n    ≥ 0, integer, columns in matrix
NB.   d    - integer in range [-(m-1):n-1], IO diagonal,
NB.          default is 0 (main diagonal)
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
NB.     S := min(m,n,⌊(n+m-|n-m-2*d|)/2⌋)
NB.
NB. Notes:
NB. - (f,s) pair defines rIOS of diagonal's fragment to be
NB.   returned

diaglios=: (0 0 _&$:) :(4 : 0)
  'd f s'=. 3 {. x , 0 _
  'm n'=. y
  F=. n (-@* ^: (0 > ])) d
  S=. <./ y , <. -: (n + m - | n - m + +: d)
  (f,:s) (] ;. 0) hds2ios F , (>: n) , S
)

NB. ---------------------------------------------------------
NB. diag
NB. Return possibly partial diagonal of rectangular matrix
NB.
NB. Syntax:
NB.   e=. [d[,f,s]] diag A
NB. where
NB.   A - m×n-matrix
NB.   d - integer in range [-(m-1):n-1], IO diagonal,
NB.       default is 0 (main diagonal)
NB.   f - integer in range [-min(m,n):min(m,n)-1], IO 1st
NB.       element of diagonal's fragment to be returned,
NB.       default is 0 (take from head)
NB.   s - integer in range [-min(m,n):min(m,n)] or ±∞, size
NB.       of diagonal's fragment to be returned, default is
NB.       _ (all elements in forward direction)
NB.   e - (0:min(m,n))-vector, d-th diagonal elements from
NB.       matrix A
NB.   m   ≥ 0
NB.   n   ≥ 0

diag=: ((<0 1)&|:) :((diaglios $) ({,) ])

NB. ---------------------------------------------------------
NB. sdiag
NB. Shift diagonal of y by values from x, i.e. make x*I+y
NB. from scalar or vector x and matrix y
NB.
NB. Syntax:
NB.   s=. x sdiag y
NB. where
NB.   y - n×n-matrix
NB.   x - numeric scalar of n-vector, shift for y's diagonal
NB.   s - n×n-matrix, equals to (y+x*idmat(#y))
NB.   n >= 0
NB.
NB. TODO:
NB. - implement for non-square matrices via diag[lios]:
NB.   sdiag=: blah-blah-blah updl

sdiag=: (+ diag) ((>: * i.) @ # @ ]) } ]  NB. FIXME! via diaglios

NB. ---------------------------------------------------------
NB. idmat
NB. Make rectangular identity matrix with shifted diagonal
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

idmat=: (0 & $:) :(= sh2id)  NB. FIXME! via diaglios

NB. ---------------------------------------------------------
NB. diagmat
NB. Make rectangular diagonal matrix with y on diagonal
NB.
NB. Syntax:
NB.   D=. [se] diagmat d
NB. where
NB.   d  - k-vector, new values for diagonal
NB.   se - complex number (s j. e):
NB.        s - integer, diagonal number for d's head,
NB.            relatively to top left corner
NB.        e - integer, diagonal number for d's tail,
NB.            relatively to bottom right corner
NB.        default se=0j0
NB.   D  - diagonal rectangular m×n-matrix with d on diagonal
NB.   k  = min(m,n), along with se defines indirectly the
NB.        shape (m,n) of D
NB.   m  ≥ 0, rows in matrix D
NB.   n  ≥ 0, columns in matrix D
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

NB. FIXME! via diaglios
diagmat=: (0 & $:) :(4 : 0)
  'r c'=. shape=. (#y) + (2&(|.@}. - {.)@(0&(<. , >.)@+.)) x  NB. find D shape
  shift=. c ((* |) ^: (0 > ])) (9 o. x)                       NB. find d head linear IO
  ios=. shift + (>: c) * i. # y                               NB. linear IOS for d
  y (ios " _) } shape $ 0                                     NB. write d into matrix of zeros
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
NB.    trl0_mt_ i. 4 3                0 trl0_mt_ i. 4 3
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

trl1=: (0 & $:) :((0 >. [) (] + (idmat $)) trl0)

NB. ---------------------------------------------------------
NB. tru1
NB. Extract unit upper triangular (trapezoidal) matrix with
NB. optional shrinking

tru1=: (0 & $:) :((0 <. [) (] + (idmat $)) tru0)
