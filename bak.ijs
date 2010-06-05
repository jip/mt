NB. Restore original eigenvectors by backward transformation
NB. from a balanced matrix or matrix pair
NB.
NB. gebakxp    Undo permutations maden by gebalxp
NB. gebakxsx   Undo scaling maden by gebals
NB. gebakxx    Form eigenvectors by backward transformation
NB.            of the matrix balanced by gebalx
NB.
NB. ggbakxp    Undo permutations maden by ggbalxp
NB. ggbakxsx   Undo scaling maden by ggbalxs
NB. ggbakxx    Form eigenvectors by backward transformation
NB.            of the pair of matrices balanced by ggbalx
NB.
NB. testgebak  Test gebakxx by general matrix given
NB. testggbak  Test ggbakxx by pair of general matrices given
NB. testbak    Adv. to make verb to test gxbakxx by matrix or
NB.            matrix pair of generator and shape given
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gebaklsl
NB. gebaklsr
NB. gebakusl
NB. gebakusr
NB.
NB. Description:
NB.   Do backward scaling for left or right eigenvectors
NB.
NB. Syntax:
NB.   'B p'=. gebakxsx A ; p ; d
NB. where
NB.   A    - n×n-matrix, eigenvectors to be transformed
NB.   p    - some parameter, transmitted immutably
NB.   d    - n-vector, diagonal of scaling matrix D
NB.   B    - n×n-matrix, transformed eigenvectors
NB.
NB. Notes:
NB. - models LAPACK's xGEBAK with 'S' option

gebaklsl=: ((0 & {::) %"1 (2 & {::)) ; (1 & {::)
gebaklsr=: ((0 & {::) *"1 (2 & {::)) ; (1 & {::)
gebakusl=: ((0 & {::) %   (2 & {::)) ; (1 & {::)
gebakusr=: ((0 & {::) *   (2 & {::)) ; (1 & {::)

NB. ---------------------------------------------------------
NB. gebaklp
NB. gebakup
NB.
NB. Description:
NB.   Do backward permutation for left or right eigenvectors
NB.
NB. Syntax:
NB.   C=. gebakxp B ; p
NB. where
NB.   B  - n×n-matrix, eigenvectors to be transformed
NB.   p  - n-vector, permutation of B
NB.   C  - n×n-matrix, transformed eigenvectors
NB.
NB. Notes:
NB. - models LAPACK's xGEBAK with 'P' option

gebaklp=: C."1~ & >
gebakup=: C.  ~ & >

NB. ---------------------------------------------------------
NB. Verb:      Balancer used:     Eigenvectors to form:
NB. gebakll    geball (lower)     left
NB. gebaklr    geball (lower)     right
NB. gebakul    gebalu (upper)     left
NB. gebakur    gebalu (upper)     right
NB.
NB. Description:
NB.   Form eigenvectors of a general matrix by backward
NB.   transformation on the computed eigenvectors, as
NB.   returned by hsein or trevc. This involves, first,
NB.   backward balance, and second, backward permutation
NB.
NB. Syntax:
NB.   C=. gebakxx A ; p ; d
NB. where
NB.   A - n×n-matrix, eigenvectors to be transformed
NB.   p - n-vector, permutation of A
NB.   d - n-vector, diagonal of scaling matrix D
NB.   C - n×n-matrix, transformed eigenvectors
NB.
NB. Notes:
NB. - models LAPACK's xGEBAK with 'B' option

gebakll=: gebaklp @ gebaklsl
gebaklr=: gebaklp @ gebaklsr
gebakul=: gebakup @ gebakusl
gebakur=: gebakup @ gebakusr

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgebak
NB.
NB. Description:
NB.   Test gebakxx by general matrix given
NB.
NB. Syntax:
NB.   testgebak A
NB. where
NB.   A - n×n-matrix

testgebak=: 3 : 0
  rcond=. (norm1 con (getrilu1p@geballu1p)) y

  ('gebakll' tmonad ((];(i.;(0&,);($&1))@#)`]`(rcond"_)`(_."_)`(_."_))) y
  ('gebaklr' tmonad ((];(i.;(0&,);($&1))@#)`]`(rcond"_)`(_."_)`(_."_))) y
  ('gebakul' tmonad ((];(i.;(0&,);($&1))@#)`]`(rcond"_)`(_."_)`(_."_))) y
  ('gebakur' tmonad ((];(i.;(0&,);($&1))@#)`]`(rcond"_)`(_."_)`(_."_))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbak
NB.
NB. Description:
NB.   Adv. to make verb to test gebakxx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testbak
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     (? @ $ 0:) testbak_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testbak_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbak_mt_ 150 200

testbal=: 1 : 'EMPTY_mt_ [ (testgebak_mt_ @ u ^: (=/))'
