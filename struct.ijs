NB. struct.ijs
NB. Structure handlers
NB.
NB. ioscv     Locate constant rows or columns in matrix
NB. ct1       Count trailing 1s in Boolean vector
NB. ct0       Count trailing 0s in vector
NB. ct0c      Count trailing zero columns in matrix
NB. ct0r      Count trailing zero rows in matrix
NB. nfv       Template conj. to make verbs negating 1st
NB.           vector (either row or column)
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

NB. Conj. to extract m-th from n-th from x: x[n[m]]
mnx=: 2 : '(m { n) { ['

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

ioscv=: *./ .=            NB. locate constant vector in matrix (column: 'm & ioscv', row: 'ioscv & n')
ct1=: +/ @: (*./\.)       NB. count trailing 1s in Boolean vector y
ct0=: ct1 @ (0&=)         NB. count trailing 0s in vector y
ct0c=: ct1 @ (0 & ioscv)  NB. count trailing zero columns in matrix y
ct0r=: ct1 @ (ioscv & 0)  NB. count trailing zero rows in matrix y

trace=: +/ @ diag         NB. matrix trace
ct=: + @ |:               NB. conjugate transpose

NB. ---------------------------------------------------------
NB. nfv
NB. Template conj. to make verbs negating 1st vector (either
NB. row or column)
NB.
NB. Syntax:
NB.   vneg=. iocios nfv iovh
NB. where
NB.   iovh   - integer in range (-r:r-1), IO in cIOS
NB.            (cios[iocios]) to select axis to negate
NB.   iocios - integer, IO in cIOS bundle (cios) to select
NB.   vneg   - verb to negate, is called as: (cios vneg A)
NB.
NB. Algorighm for verb vneg:
NB.   1) take cIOS from x[m]
NB.   2) set imaginary part to 1 in its n-th item
NB.   3) convert transformed cIOS to IOS
NB.   4) use this IOS to update by verb (-) submatrix in y
NB.
NB. Application:
NB. - let cios=. (ciosYZ , ciosT , ciosR ,: ciosL) is cIOS
NB.   bundle for vectors y and z, scalar τ, and submatrices
NB.   R and L, respectively; vector y is stored vertically
NB. - to make verb to negate 1st column of submatrix R:
NB.     vneg=: 2 nfv 1
NB. - to make verb to negate 1st row of submatrix L:
NB.     vneg=: 3 nfv 0
NB. - to make verb to negate 1st column of submatrix R and
NB.   then to negate 1st row of submatrix L
NB.     vneg=: [ (3 nfv 0) (2 nfv 1)

nfv=: 2 : '((- upd)~ (cios2ios @ ((9&o. j. 1:) updl n) @ (m&{)))~'

NB. ---------------------------------------------------------
NB. IOS explorers
NB. monad     dyad    triad   tetrad     pentad     hexad     heptad    octad    ennead decade
NB.           couplet triplet quadruplet quintuplet sextuplet septuplet octuplet
NB. singleton duet    trio    quartet    quintet    sextet    septet    octet

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
upd2ri=: 2 : '((u gici ((1{n),2)) (u gi 1) (u gici ((0{n),0)))`(rios2ios@(1 mnx n))`] }'

NB. ---------------------------------------------------------
NB. upd3ri
NB. Template conj. to make verbs to update subarray by pentad
NB. indirectly, addressing mode is rIOS
NB.
NB. Syntax:
NB.   vapp=. u0`u1`u2`u3`u4 upd3ri io0 io1 io2
NB. where
NB.   rios        - k×2×2-array of integers, vector of rIOSs
NB.   io0 io1 io2 - integers, IOs in rios
NB.   A           - input matrix
NB.   A0 A1 A2    - 2-rank arrays, submatrices in A:
NB.                   A0 -: (io0{rios) (] ;. 0) A
NB.                   A1 -: (io1{rios) (] ;. 0) A
NB.                   A2 -: (io2{rios) (] ;. 0) A
NB.   u0 u2 u4    - monads to adjust A0 A1 A2, respectively
NB.   u1 u3       - dyads to glue staff
NB.   vapp        - verb to update A2, is called as:
NB.                   Aupd=. rios vapp A
NB.   Aupd        - A with submatrix A2 replaced by value:
NB.                  (u4 A2) u3 (u2 A1) u1 (u0 A0)
NB.   k           > 0
NB.
NB. Example:
NB. - following replaces 24 by _1220 in array (i. 10 10):
NB.     (3 2 2 $ 4 9 1 1 3 5 1 1 2 4 1 1) (%:`+`*:`-`-: upd3ri 0 1 2) i. 10 10
NB.   since
NB.     _1220 -: (-: 24) - (*: 35) + (%: 49)

upd3ri=: 2 : '((u gici ((2{n),4)) (u gi 3) (u gici ((1{n),2)) (u gi 1) (u gici ((0{n),0)))`(rios2ios@(2 mnx n))`] }'

NB. Adv. to model 'map' accepting lIOS
NB. y[u[1](u[0](y[x]),y)] := u[0](y[x])
NB. Aupd=. lios (u0`u1 mapl) A
mapl=: 1 : '((u gi 0) @ ({,)) ((u gi 1) }) ]'

NB. Adv. to model 'map' accepting rIOS to read and
NB. using lIOS to write
NB. y[u[1](u[0](y[x]),y)] := u[0](y[x])
NB. Aupd=. rIOS (u0`u1 maprl) A
maprl=: 1 : '(u gic 0) ((u gi 1) }) ]'

NB. Conj. to model 'mapi' accepting rIOS
NB. y[u[1](y[x[n]],y)] := u[0](y[x[n]])
NB. Aupd=. rios (u mapri (ioFrom,ioTo)) A
mapri=: 2 : '(u uci (0{n))`(rios2ios@(1 mnx n))`] }'

NB. Conj. to model 'upd3i' accepting rIOS
NB. y[x[n[2]]] := u[3](u[4](y[x[n[2]]]),u[1](u[2](y[x[n[1]]]),u[0](y[x[n[0]]])))
map3ri=: 2 : '((u gici ((2{n),4)) (u gi 3) (u gici ((1{n),2)) (u gi 1) (u gici ((0{n),0)))`(rios2ios@(3 mnx n))`] }'

NB. Conj. to model 'upd4i' accepting rIOS
NB. y[x[n[4]]] := u[5](u[6](y[x[n[3]]]),u[3](u[4](y[x[n[2]]]),u[1](u[2](y[x[n[1]]]),u[0](y[x[n[0]]]))))
map4ri=: 2 : '((u gici ((3{n),6)) (u gi 5) ((u gici ((2{n),4)) (u gi 3) (u gici ((1{n),2)) (u gi 1) (u gici ((0{n),0))))`(rios2ios@(4 mnx n))`] }'

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
NB. - implement for non-square matrices

sdiag=: (+ diag) ((>: * i.) @ # @ ]) } ]

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

idmat=: (0 & $:) :(= sh2id)

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
