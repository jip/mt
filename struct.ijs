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
NB. IO   = index of
NB. IOS  = indices of
NB. lIOS = linear IOS
NB. rIOS = rectangular IOS
NB. cIOS = complex IOS
NB.
NB. Following are equivalents:
NB.   (3 5 _7 ,: 2 _3 4) ] ;. 0 report
NB.   (< 3 4 ; 7 6 5 ; _10 _9 _8 _7) { report
NB.   (cios2ios 3j2 5j_3 _7j4) { report
NB.   (cios2rios 3j2 5j_3 _7j4) ] ;. 0 report
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

NB. m-th from y from x: x[y[m]]
myx=: 1 : '{~ (m & {)'

NB. m-th from n from y from x: x[y[n[m]]]
mnyx=: 2 : '{~ ((m { n) & {)'

NB. m-th from x: x[m]
mx=: 1 : 'm { ['

NB. m-th from n from x: x[n[m]]
mnx=: 2 : '(m { n) { ['

NB. aggregated extractor: u(x[y[n[0]]],x[y[n[1]]])
mnyx01=: 2 : '(0 mnyx n) u (1 mnyx n)'

NB. TODO: mnyx012

NB. ---------------------------------------------------------
NB. Wrappers

NB. cIOS wrapper: u(cios2ios(x),y)
cioswrap=: 1 : '(u~ cios2ios)~'

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
NB. TODO: sink cios2ios to low-level indirect extractors

NB. Model direct 'from' accepting cIOS: y[x]
fromc=: ({~ cios2ios)~

NB. Adv. to model indirect 'from': y[x[m]]
fromi=: 1 : '(m myx)~'

NB. Adv. to model 'fromi' accepting cIOS
fromci=: 1 : '(m fromi) cioswrap'

NB. Adv. to model indirect 'amend' accepting cIOS
NB. y[x[n]] := u(x,y)
amendci=: 2 : 'u`(cios2ios @ (n mx))`] }'

NB. Adv. to update subarray by monad: y[x] := u(y[x])
upd=: 1 : '(u @ {)`[`]}'

NB. Conj. to model 'upd' accepting lIOS
updl=: 2 : '((n"_)})~ (u @ (n & ({,)))'

NB. Conj. to model indirect 'upd':
NB. y[x[n]] := u(y[x[n]])
updi=: 2 : '((u @ (n myx))~)`(n mx)`] }'

NB. Conj. to model 'updi' accepting cIOS
updci=: 2 : '(u updi n) cioswrap'

NB. Conj. to indirect update subarray by dyad:
NB. y[x[n[1]]] := u(y[x[n[0]]],y[x[n[1]]])
upd2i=: 2 : '((u mnyx01 n)~)`(1 mnx n)`] }'

NB. Conj. to model 'upd2i' accepting cIOS
upd2ci=: 2 : '(u upd2i n) cioswrap'

NB. Conj. to indirect update subarray by triad:
NB. y[x[n[2]]] := u[1](u[0](y[x[n[0]]],y[x[n[1]]]),y[x[n[2]]])
upd3i=: 2 : '(((((0{u)`:6) mnyx01 n) ((1{u)`:6) (2 mnyx n))~)`(2 mnx n)`] }'

NB. Conj. to model 'upd3i' accepting cIOS
upd3ci=: 2 : '(u upd3i n) cioswrap'

NB. Conj. to indirect map subarray to another one by monad:
NB. y[x[n[1]]] := u(y[x[n[0]]])
mapi=: 2 : '((u @ (0 mnyx n))~)`(1 mnx n)`] }'

NB. Conj. to model 'mapi' accepting cIOS
mapci=: 2 : '(u mapi n) cioswrap'

NB. Conj. to model 'mapi' accepting lIOS
NB. y[x[n[1]]] := u(y[x[n[0]]])
NB. FIXME!
NB. mapli=: 2 : '((0 mnx n) (u @ ({,)) ]) ((1 mnx n) }) ]'
NB. mapli=: 2 : '((0 {:: [) (u ;. 0) ]) ((1 {:: [) }) ]'    NB. (rIOSget ; lIOSwrite) (u mapli 0 1) A

NB. Conj. to indirect map two subarrays to another one by
NB. dyad: y[x[n[2]]] := u(y[x[n[0]]],y[x[n[1]]])
map2i=: 2 : '((u mnyx01 n)~)`(2 mnx n)`] }'

NB. Conj. to model 'map2i' accepting cIOS
map2ci=: 2 : '(u map2i n) cioswrap'

NB. Conj. to indirect map three subarrays to another one by
NB. triad:
NB. y[x[n[3]]] := u[1](u[0](y[x[n[0]]],y[x[n[1]]]),y[x[n[2]]])
map3i=: 2 : '(((((0{u)`:6) mnyx01 n) ((1{u)`:6) (2 mnyx n))~)`(3 mnx n)`] }'

NB. Conj. to model 'map3i' accepting cIOS
map3ci=: 2 : '(u map3i n) cioswrap'

NB. Conj. to indirect map four subarrays to another one by
NB. tetrad:
NB. y[x[n[4]]] := u[2](u[1](u[0](y[x[n[0]]],y[x[n[1]]]),y[x[n[2]]]),y[x[n[3]]])
map4i=: 2 : '((((((0{u)`:6) mnyx01 n) ((1{u)`:6) (2 mnyx n)) ((2{u)`:6) (3 mnyx n))~)`(4 mnx n)`] }'

NB. Conj. to model 'map3i' accepting cIOS
map4ci=: 2 : '(u map4i n) cioswrap'

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
