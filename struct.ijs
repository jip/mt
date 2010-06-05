NB. struct.ijs
NB. Structure handlers
NB.
NB. ioscv     Locate constant vector in matrix
NB. ct1       Count trailing 1s in Boolean vector
NB. ct0       Count trailing 0s in vector
NB. ct0c      Count trailing zero columns in matrix
NB. ct0r      Count trailing zero rows in matrix
NB.
NB. hds2ios   IOS from head, delta and size
NB. ht2ios    IOS from head and tail
NB. hs2ios    IOS from head and size
NB. cios2ios  Convert complex IOS (cIOS) to IOS
NB. from      Adv. to model indirect 'from'
NB. cfrom     Adv. to model 'from' accepting cIOS
NB. upd       Conj. to update subarray by monad
NB. cupd      Conj. to model 'upd' accepting cIOS
NB. upd2      Conj. to update subarray by dyad
NB. cupd2     Conj. to model 'upd2' accepting cIOS
NB. map       Conj. to map subarray to another one by monad
NB. cmap      Conj. to model 'map' accepting cIOS
NB. map2      Conj. to map two subarrays to another one by
NB.           dyad
NB. cmap2     Conj. to model 'map2' accepting cIOS
NB.
NB. append   Enhance append (A)
NB. stitch   Enhance stitch (A)
NB.
NB. rws      Rows weighted sum
NB.
NB. idmat    Make rectangular identity matrix with shifted
NB.          diagonal
NB. diagmat  Make rectangular diagonal matrix with shifted
NB.          diagonal
NB.
NB. tru      Extract upper triangular (trapezoidal) matrix
NB. trl      Extract lower triangular (trapezoidal) matrix
NB. tru0     Extract strictly upper triangular (trapezoidal)
NB.          matrix
NB. trl0     Extract strictly lower triangular (trapezoidal)
NB.          matrix
NB. tru1     Extract unit upper triangular (trapezoidal)
NB.          matrix
NB. trl1     Extract unit lower triangular (trapezoidal)
NB.          matrix
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

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

NB. indirect extractors
myx=: 1 : '{~ (m & {)'                  NB. indirect m from y from x: (m{y){x
mnyx=: 2 : '{~ ((m { n) & {)'           NB. indirect m from n from y from x: ((m{n){y){x
mx=: 1 : 'm { ['                        NB. indirect m from x: m{x
mnx=: 2 : '(m { n) { ['                 NB. indirect m from n from x: (m{n){x
mnyx01=: 2 : '(0 mnyx n) u (1 mnyx n)'  NB. aggregated extractor

NB. cIOS wrapper
cioswrap=: 1 : '(u~ cios2ios)~'

NB. if (v is noun) then: update items defined by (v) in (y) by verb (u)
NB. else: update items defined by linear IOS (x v y) in (y) by verb (u)
uvfyu=: 2 : '(u @ (v { ])) v} ]'

NB. =========================================================
NB. Interface

ioscv=: *./ .=            NB. locate constant vector in matrix (column: 'm & ioscv', row: 'ioscv & n')
ct1=: +/ @: (*./\.)       NB. count trailing 1s in Boolean vector y
ct0=: ct1 @ (0&=)         NB. count trailing 0s in vector y
ct0c=: ct1 @ (0 & ioscv)  NB. count trailing zero columns in matrix y
ct0r=: ct1 @ (ioscv & 0)  NB. count trailing zero rows in matrix y

NB. ---------------------------------------------------------
NB. Linear IOS / Rectangular IOS (rIOS) / Complex IOS (cIOS)
NB.
NB. Following are equivalents:
NB.    (3 5 _7 ,: 2 _3 4) ] ;. 0 report
NB.    (< 3 4 ; 7 6 5 ; _10 _9 _8 _7) { report
NB.    (cios2ios 3j2 5j_3 _7j4) { report
NB.    (cios2rios 3j2 5j_3 _7j4) ] ;. 0 report

NB. generators
hds2ios=: + ` (* i.)/             NB. (2{y)-vector of integers from head (0{y) by delta (1{y)
ht2ios=: ] + (i. @ -)             NB. (x-y)-vector of integers from head y to tail (x-1): y (y+1) ... (x-1)
hs2ios=: [ + ((] * i. @ *) sgn)~  NB. y-vector of integers from head x of size y, models verb's (u;.0) rIOS

NB. converters
cios2ios=: < " 1 @ (< @ hs2ios/ " 1 @: +.)  NB. convert cIOS to IOS; side effect: output is incorrect for IOS with length less than array's rank

NB. explorers
NB. FIXME: sink cios2ios to low-level individual selectors
from=: 1 : '(m myx)~'                          NB. model indirect 'from': (m{x){y
cfrom=: 1 : '(m from) cioswrap'                NB. model 'from' accepting cIOS
upd=: 2 : '((u @ (n myx))~)`(n mx)`] }'        NB. update subarray by monad: (u ((n{x){y)) (n{x) } y
cupd=: 2 : '(u upd n) cioswrap'                NB. model 'upd' accepting cIOS
updu=: 1 : '(u @ {)`[`]}'                      NB. update in y items defined by x by verb u
upd2=: 2 : '((u mnyx01 n)~)`(1 mnx n)`] }'     NB. update subarray by dyad: (((({.n){x){y) u ((({:n){x){y)) (({:n){x) } y
cupd2=: 2 : '(u upd2 n) cioswrap'              NB. model 'upd2' accepting cIOS
upd3=: 2 : '(((((0{u)`:6) mnyx01 n) ((1{u)`:6) (2 mnyx n))~)`(2 mnx n)`] }'  NB. update subarray by gerund
cupd3=: 2 : '(u upd3 n) cioswrap'              NB. model 'upd3' accepting cIOS
map=: 2 : '((u @ (0 mnyx n))~)`(1 mnx n)`] }'  NB. map subarray to another one by monad: (u ((({.n){x){y)) (({:n){x) } y
cmap=: 2 : '(u map n) cioswrap'                NB. model 'map' accepting cIOS
map2=: 2 : '((u mnyx01 n)~)`(2 mnx n)`] }'     NB. map two subarrays to another one by dyad: (((({.n){x){y) u (((1{n){x){y)) (({:n){x) } y
cmap2=: 2 : '(u map2 n) cioswrap'              NB. model 'map2' accepting cIOS

NB. ---------------------------------------------------------
NB. nfv
NB. Template conj to make verbs negating 1st vector (either
NB. row or column)
NB.
NB. Syntax:
NB.   vneg=. iocios nfv iovh
NB. where
NB.   iovh   - integer ub rabge (-r:r-1), IO in cIOS
NB.            (iocios{cios) to select axis to negate
NB.   iocios - integer, IO in cIOS bundle (cios) to select
NB.   vneg   - verb to negate, is called as: (cios vneg A),
NB.            implements steps:
NB.              1) take 2-vector (cIOS) (m{x)
NB.              2) set imaginary part to 1 in its n-th item
NB.              3) convert transformed cIOS to IOS
NB.              4) use this IOS to update by verb (-)
NB.                 submatrix in y
NB.
NB. Application:
NB. - let cios=. (ciosY , ciosZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y, z, scalar τ, and
NB.   submatrices R and L, respectively; vector y is stored
NB.   vertically
NB. - to make verb to negate 1st column of submatrix R:
NB.     vneg=: 3 nfv 1
NB. - to make verb to negate 1st row of submatrix L:
NB.     vneg=: 4 nfv 0
NB. - to make verb to negate 1st column of submatrix R and
NB.   then to negate 1st row of submatrix L
NB.     vneg=: [ (4 nfv 0) (3 nfv 1)

nfv=: 2 : '((- updu)~ (cios2ios @ ((9&o. j. 1:) uvfyu n) @ (m&{)))~'

NB. ---------------------------------------------------------
NB. append
NB. Adverb to enhance append
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
NB. Adverb to enhance stitch
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
