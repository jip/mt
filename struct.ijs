NB. struct.ijs
NB. Matrix structure handlers
NB.
NB. idmat    make rectangular identity matrix with shifted diagonal
NB. diagmat  make rectangular diagonal matrix with shifted diagonal
NB. ltri     extract lower triangular (trapezoidal) matrix
NB. utri     extract upper triangular (trapezoidal) matrix
NB. ltri0    extract strictly lower triangular (trapezoidal) matrix
NB. utri0    extract strictly upper triangular (trapezoidal) matrix
NB. ltri1    extract unit lower triangular (trapezoidal) matrix
NB. utri1    extract unit upper triangular (trapezoidal) matrix

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. convert shape y to IOS differences table
sh2id=: {. -~/&i. {:

NB. template conj. to extract rectangular matrix
NB. circumscribing the triangular (trapezoidal) matrix
NB. starting from diagonal number x in the matrix y
tricut=: 2 : '((m & *) @: (<./ " 1) @ v $) {. ]'

NB. extract lower triangular (trapezoidal) matrix
ltricut=: _1 1 tricut ((+ {.) ,. ])

NB. extract upper triangular (trapezoidal) matrix
utricut=: 1 _1 tricut (] ,. (-~ {:))

NB. template conj. to extract triangular (trapezoidal)
NB. matrix starting from diagonal number x in the rectangular
NB. circumscribing matrix y
tri=: 2 : '0&$: : ([ (] * (u~ sh2id@$)) v)'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. idmat                                               _ _ _
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
NB. diagmat                                             _ _ _
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
NB. ltri                                                _ _ _
NB. Extract lower triangular (trapezoidal) matrix with
NB. optional shrinking
NB.
NB. Examples:
NB.    ltri i. 3 4                    0 ltri i. 3 4
NB. 0 0  0                         0 0  0
NB. 4 5  0                         4 5  0
NB. 8 9 10                         8 9 10
NB.    1 ltri i. 3 4                  _1 ltri i. 3 4
NB. 0 1  0  0                      4 0
NB. 4 5  6  0                      8 9
NB. 8 9 10 11
NB.    1 ltri i. 4 3                  _1 ltri i. 4 3
NB. 0  1  0                        3  0  0
NB. 3  4  5                        6  7  0
NB. 6  7  8                        9 10 11
NB. 9 10 11

ltri=: (>:~ 0&>.) tri ltricut

NB. ---------------------------------------------------------
NB. utri                                                _ _ _
NB. Extract upper triangular (trapezoidal) matrix with
NB. optional shrinking
NB.
NB. Examples:
NB.    utri i. 3 4                    0 utri i. 3 4
NB. 0 1  2  3                      0 1  2  3
NB. 0 5  6  7                      0 5  6  7
NB. 0 0 10 11                      0 0 10 11
NB.    1 utri i. 3 4                  _1 utri i. 3 4
NB. 1 2  3                         0 1  2  3
NB. 0 6  7                         4 5  6  7
NB. 0 0 11                         0 9 10 11
NB.    1 utri i. 4 3                  _1 utri i. 4 3
NB. 1 2                            0 1  2
NB. 0 5                            3 4  5
NB.                                0 7  8
NB.                                0 0 11

utri=: (<:~ 0&<.) tri utricut

NB. ---------------------------------------------------------
NB. ltri0                                               _ _ _
NB. Extract strictly lower triangular (trapezoidal) matrix
NB. with optional shrinking
NB.
NB. Examples:
NB.    ltri0_mt_ i. 4 3               0 ltri0_mt_ i. 4 3
NB. 0  0  0                        0  0  0
NB. 3  0  0                        3  0  0
NB. 6  7  0                        6  7  0
NB. 9 10 11                        9 10 11
NB.    1 ltri0 i. 4 3                 _1 ltri0 i. 4 3
NB. 0  0  0                        0  0 0
NB. 3  4  0                        6  0 0
NB. 6  7  8                        9 10 0
NB. 9 10 11
NB.    1 ltri0 i. 3 4                 _1 ltri0 i. 3 4
NB. 0 0  0 0                       0 0
NB. 4 5  0 0                       8 0
NB. 8 9 10 0

ltri0=: (>~ 0&>.) tri ltricut

NB. ---------------------------------------------------------
NB. utri0                                               _ _ _
NB. Extract strictly upper triangular (trapezoidal) matrix
NB. with optional shrinking
NB.
NB. Examples:
NB.    utri0 i. 3 4                   0 utri0 i. 3 4
NB. 0 1 2  3                       0 1 2  3
NB. 0 0 6  7                       0 0 6  7
NB. 0 0 0 11                       0 0 0 11
NB.    1 utri0 i. 3 4                 _1 utri0 i. 3 4
NB. 0 2 3                          0 1  2  3
NB. 0 0 7                          0 5  6  7
NB. 0 0 0                          0 0 10 11
NB.    1 utri0 i. 4 3                 _1 utri0 i. 4 3
NB. 0 2                            0 1 2
NB. 0 0                            0 4 5
NB.                                0 0 8
NB.                                0 0 0

utri0=: (<~ 0&<.) tri utricut

NB. ---------------------------------------------------------
NB. ltri1                                               _ _ _
NB. Extract unit lower triangular (trapezoidal) matrix with
NB. optional shrinking

ltri1=: (0 & $:) :((0 >. [) (] + (idmat $)) ltri0)

NB. ---------------------------------------------------------
NB. utri1                                               _ _ _
NB. Extract unit upper triangular (trapezoidal) matrix with
NB. optional shrinking

utri1=: (0 & $:) :((0 <. [) (] + (idmat $)) utri0)
