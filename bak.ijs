NB. Restore original eigenvectors by backward transformation
NB.
NB. gebakxp    Undo permutations maden by gebalxp
NB. gebakxsx   Undo scaling maden by gebals
NB. gebakxx    Form eigenvectors by backward transformation
NB.            of the matrix balanced by gebalx
NB.
NB. testgebak  Test gebakxx by square matrix
NB. testbak    Adv. to make verb to test gebakxx by matrix of
NB.            generator and shape given
NB.
NB. Copyright 2010,2011,2013,2017,2018,2020,2021,2023,2024,
NB.           2025 Igor Zhuravlov
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

NB. =========================================================
NB. Configuration

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
NB.   p    - some argument, transmitted immutably
NB.   d    - n-vector, diagonal of scaling matrix D
NB.   B    - n×n-matrix, transformed eigenvectors
NB.
NB. Application:
NB. - models xGGBAK('S','L')
NB.     'CD pl'=. gebakxsx EF ; ({. plr) ; ({. dlr)
NB. - models xGGBAK('S','R')
NB.     'CD pr'=. gebakxsx EF ; ({: plr) ; ({: dlr)
NB.
NB. Notes:
NB. - gebakusl models LAPACK's xGxBAK('S','L')
NB. - gebakusr models LAPACK's xGxBAK('S','R')

gebaklsl=: (0&{:: (%"1) 2&{::) ; 1&{::
gebaklsr=: (0&{:: (*"1) 2&{::) ; 1&{::
gebakusl=: (0&{:: (%"2) 2&{::) ; 1&{::
gebakusr=: (0&{:: (*"2) 2&{::) ; 1&{::

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
NB. Application:
NB. - models xGGBAK('P','L')
NB.     AB=. gebakup CD ; ({. plr)
NB. - models xGGBAK('P','R')
NB.     AB=. gebakup CD ; ({: plr)
NB.
NB. Notes:
NB. - gebakup models LAPACK's xGxBAK('P')

gebaklp=: C."1~&>/
gebakup=: C."2~&>/

NB. ---------------------------------------------------------
NB. Monad      Balancer used     Eigenvectors to form
NB. gebakll    geball (lower)    left
NB. gebaklr    geball (lower)    right
NB. gebakul    gebalu (upper)    left
NB. gebakur    gebalu (upper)    right
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
NB. Application:
NB. - models xGGBAK('B','L')
NB.     AB=. gebakul EF ; ({. plr) ; ({. dlr)
NB. - models xGGBAK('B','R')
NB.     AB=. gebakur EF ; ({: plr) ; ({: dlr)
NB.
NB. Notes:
NB. - gebakul models LAPACK's xGxBAK('B','L')
NB. - gebakur models LAPACK's xGxBAK('B','R')

gebakll=: gebaklp@gebaklsl
gebaklr=: gebaklp@gebaklsr
gebakul=: gebakup@gebakusl
gebakur=: gebakup@gebakusr

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgebak
NB.
NB. Description:
NB.   Test gebakxx by square matrix
NB.
NB. Syntax:
NB.   log=. testgebak A
NB. where
NB.   A   - n×n-matrix
NB.   log - 6-vector of boxes, test log

testgebak=: 3 : 0
  'rcondl rcondu'=. (geconi , gecon1) y

  log=.          ('gebakll' tmonad ((] ; (i. ; $&1)@#)`]`(rcondl"_)`nan`nan)) y
  log=. log lcat ('gebaklr' tmonad ((] ; (i. ; $&1)@#)`]`(rcondl"_)`nan`nan)) y
  log=. log lcat ('gebakul' tmonad ((] ; (i. ; $&1)@#)`]`(rcondu"_)`nan`nan)) y
  log=. log lcat ('gebakur' tmonad ((] ; (i. ; $&1)@#)`]`(rcondu"_)`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testbak
NB.
NB. Description:
NB.   Adv. to make verb to test gebakxx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testbak) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log
NB.
NB. Application:
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testbak_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testbak_mt_ 150 150
NB. - test by random square complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testbak_mt_ 150 150

testbak=: 1 : 'nolog_mt_`(testgebak_mt_@u)@.(=/)'
