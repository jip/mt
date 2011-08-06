NB. Restore original eigenvectors by backward transformation
NB. from a balanced matrix or matrix pair
NB.
NB. gxbakxp    Undo permutations maden by gxbalxp
NB. gxbakxsx   Undo scaling maden by gxbals
NB. gxbakxx    Form eigenvectors by backward transformation
NB.            of the matrix balanced by gxbalx
NB.
NB. testgxbak  Test gxbakxx by general matrix given
NB. testbak    Adv. to make verb to test gxbakxx by matrix of
NB.            generator and shape given
NB.
NB. Version: 0.6.8 2010-10-30
NB.
NB. Copyright 2010 Igor Zhuravlov
NB.
NB. This file is part of mt
NB.
NB. mt is free software: you can redistribute it and/or
NB. modify it under the terms of the GNU Lesser General
NB. Public License as published by the Free Software
NB. Foundation, either version 3 of the License, or (at your
NB. option) any later version.
NB.
NB. mt is distributed in the hope that it will be useful, but
NB. WITHOUT ANY WARRANTY; without even the implied warranty
NB. of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
NB. See the GNU Lesser General Public License for more
NB. details.
NB.
NB. You should have received a copy of the GNU Lesser General
NB. Public License along with mt. If not, see
NB. <http://www.gnu.org/licenses/>.

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gxbaklsl
NB. gxbaklsr
NB. gxbakusl
NB. gxbakusr
NB.
NB. Description:
NB.   Do backward scaling for left or right eigenvectors
NB.
NB. Syntax:
NB.   'B p'=. gxbakxsx A ; p ; d
NB. where
NB.   A    - n×n-matrix, eigenvectors to be transformed
NB.   p    - some parameter, transmitted immutably
NB.   d    - n-vector, diagonal of scaling matrix D
NB.   B    - n×n-matrix, transformed eigenvectors
NB.
NB. Notes:
NB. - gxbakusl models LAPACK's xGxBAK('S','L')
NB. - gxbakusr models LAPACK's xGxBAK('S','R')

gxbaklsl=: ((0 & {::) %"1 (2 & {::)) ; (1 & {::)
gxbaklsr=: ((0 & {::) *"1 (2 & {::)) ; (1 & {::)
gxbakusl=: ((0 & {::) %   (2 & {::)) ; (1 & {::)
gxbakusr=: ((0 & {::) *   (2 & {::)) ; (1 & {::)

NB. ---------------------------------------------------------
NB. gxbaklp
NB. gxbakup
NB.
NB. Description:
NB.   Do backward permutation for left or right eigenvectors
NB.
NB. Syntax:
NB.   C=. gxbakxp B ; p
NB. where
NB.   B  - n×n-matrix, eigenvectors to be transformed
NB.   p  - n-vector, permutation of B
NB.   C  - n×n-matrix, transformed eigenvectors
NB.
NB. Notes:
NB. - gxbakup models LAPACK's xGxBAK('P')

gxbaklp=: C."1~ & >/
gxbakup=: C.  ~ & >/

NB. ---------------------------------------------------------
NB. Verb:      Balancer used:     Eigenvectors to form:
NB. gxbakll    gxball (lower)     left
NB. gxbaklr    gxball (lower)     right
NB. gxbakul    gxbalu (upper)     left
NB. gxbakur    gxbalu (upper)     right
NB.
NB. Description:
NB.   Form eigenvectors of a general matrix by backward
NB.   transformation on the computed eigenvectors, as
NB.   returned by hsein or trevc. This involves, first,
NB.   backward balance, and second, backward permutation
NB.
NB. Syntax:
NB.   C=. gxbakxx A ; p ; d
NB. where
NB.   A - n×n-matrix, eigenvectors to be transformed
NB.   p - n-vector, permutation of A
NB.   d - n-vector, diagonal of scaling matrix D
NB.   C - n×n-matrix, transformed eigenvectors
NB.
NB. Notes:
NB. - gxbakul models LAPACK's xGxBAK('B','L')
NB. - gxbakur models LAPACK's xGxBAK('B','R')

gxbakll=: gxbaklp @ gxbaklsl
gxbaklr=: gxbaklp @ gxbaklsr
gxbakul=: gxbakup @ gxbakusl
gxbakur=: gxbakup @ gxbakusr

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgxbak
NB.
NB. Description:
NB.   Test gxbakxx by general matrix given
NB.
NB. Syntax:
NB.   testgxbak A
NB. where
NB.   A - n×n-matrix

testgxbak=: 3 : 0
  rcond=. gecon1 y

  ('gxbakll' tmonad ((];(i.;($&1))@#)`]`(rcond"_)`(_."_)`(_."_))) y
  ('gxbaklr' tmonad ((];(i.;($&1))@#)`]`(rcond"_)`(_."_)`(_."_))) y
  ('gxbakul' tmonad ((];(i.;($&1))@#)`]`(rcond"_)`(_."_)`(_."_))) y
  ('gxbakur' tmonad ((];(i.;($&1))@#)`]`(rcond"_)`(_."_)`(_."_))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbak
NB.
NB. Description:
NB.   Adv. to make verb to test gxbakxx by matrix of
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
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testbak_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testbak_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbak_mt_ 150 150

testbak=: 1 : 'EMPTY_mt_ [ (testgxbak_mt_ @ u ^: (=/))'
