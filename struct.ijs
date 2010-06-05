NB. struct.ijs
NB. Matrix structure handlers
NB.
NB. append   adverb to enhance append
NB. stitch   adverb to enhance stitch
NB.
NB. idmat    make rectangular identity matrix with shifted diagonal
NB. diagmat  make rectangular diagonal matrix with shifted diagonal
NB.
NB. tru     extract upper triangular (trapezoidal) matrix
NB. trl     extract lower triangular (trapezoidal) matrix
NB. tru0    extract strictly upper triangular (trapezoidal) matrix
NB. trl0    extract strictly lower triangular (trapezoidal) matrix
NB. tru1    extract unit upper triangular (trapezoidal) matrix
NB. trl1    extract unit lower triangular (trapezoidal) matrix
NB.
NB. Version: 1.0.0 2009-06-01
NB. Copyright: Igor Zhuravlov, igor at uic.dvgu.ru
NB. License: GNU GPL

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

NB. =========================================================
NB. Interface

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
NB.   k  ≥ 0, along with se defines indirectly the shape
NB.        (m,n) of D
NB.   m  ≥ 0, rows in D
NB.   n  ≥ 0, columns in D
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
NB. tru1                                               _ _ _
NB. Extract unit upper triangular (trapezoidal) matrix with
NB. optional shrinking

tru1=: (0 & $:) :((0 <. [) (] + (idmat $)) tru0)
