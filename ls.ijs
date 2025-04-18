NB. Solve overdetermined or underdetermined linear monomial
NB. equation
NB.
NB. gelsxxx    Solve overdetermined or underdetermined linear
NB.            system involving a matrix of full rank, or its
NB.            [conjugate-]transpose
NB.
NB. testgels1  Test gelsxxx by general matrix and single RHS
NB. testgels3  Test gelsxxx by general matrix and multiple
NB.            RHS
NB. testgels   Adv. to make verb to test gelsxxx by a matrix
NB.            generated by generator given
NB. testls     Adv. to make verb to test xxlsxxx by matrix of
NB.            datatype and shape given
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

NB. ---------------------------------------------------------
NB. dlarnv2
NB. zlarnv2
NB.
NB. Description:
NB.   Generate a random array of shape given
NB.
NB. Syntax:
NB.   R=. dlarnv2 sh
NB.   C=. zlarnv2 sh
NB. where
NB.   R  - sh-array of real values r ~ U(-1,1)
NB.   C  - sh-array of complex values c having:
NB.          Re(c) ~ U(-1,1)
NB.          Im(c) ~ U(-1,1)
NB.   sh - vector of non-negative integers, the shape of
NB.        matrix R or C
NB.
NB. Notes:
NB. - for scalar sh:
NB.   - dlarnv2 models LAPACK's DLARNV(2)
NB.   - zlarnv2 models LAPACK's ZLARNV(2)

dlarnv2=: _1 1& randu
zlarnv2=: _1 1&(randu j. randu)

NB. ---------------------------------------------------------
NB. dqrt131
NB. zqrt131
NB.
NB. Description:
NB.   Generate a full-rank matrix scaled normally to have
NB.   norm in range [safe_min/precision , precision/safe_min]
NB.
NB. Syntax:
NB.   A=. dqrt131 (m,n)
NB.   B=. zqrt131 (m,n)
NB. where
NB.   A - m×n-matrix of full rank, real
NB.   B - m×n-matrix of full rank, complex
NB.
NB. Notes:
NB. - dqrt131 models LAPACK's DQRT13(1)
NB. - zqrt131 models LAPACK's ZQRT13(1)

dqrt131=: (setdiag~ '' ;~ diag ([ +  negneg        ) norm1tc@(({."1~ #)^:(</@$)))@dlarnv2
zqrt131=: (setdiag~ '' ;~ diag ([ + (negneg~ 9&o.)~) norm1tc@(({."1~ #)^:(</@$)))@zlarnv2

NB. ---------------------------------------------------------
NB. dqrt132
NB. zqrt132
NB.
NB. Description:
NB.   Generate a full-rank matrix scaled up to have large
NB.   norm
NB.
NB. Syntax:
NB.   A=. dqrt132 (m,n)
NB.   B=. zqrt132 (m,n)
NB. where
NB.   A - m×n-matrix of full rank and large norm, real
NB.   B - m×n-matrix of full rank and large norm, complex
NB.
NB. Notes:
NB. - dqrt132 models LAPACK's DQRT13(2)
NB. - zqrt132 models LAPACK's ZQRT13(2)

dqrt132=: (scl~ (FP_EPS % FP_SFMIN) ,~ normm)@dqrt131
zqrt132=: (scl~ (FP_EPS % FP_SFMIN) ,~ normm)@zqrt131

NB. ---------------------------------------------------------
NB. dqrt133
NB. zqrt133
NB.
NB. Description:
NB.   Generate a full-rank matrix scaled down to have small
NB.   norm
NB.
NB. Syntax:
NB.   A=. dqrt133 (m,n)
NB.   B=. zqrt133 (m,n)
NB. where
NB.   A - m×n-matrix of full rank and small norm, real
NB.   B - m×n-matrix of full rank and small norm, complex
NB.
NB. Notes:
NB. - dqrt133 models LAPACK's DQRT13(3)
NB. - zqrt133 models LAPACK's ZQRT13(3)

dqrt133=: (scl~ (FP_SFMIN % FP_EPS) ,~ normm)@dqrt131
zqrt133=: (scl~ (FP_SFMIN % FP_EPS) ,~ normm)@zqrt131

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gelsax
NB. gelsacx
NB.
NB. Description:
NB.   For m×n-matrix A of full rank:
NB.   1) if m>=n then gelsax finds the least squares solution
NB.      of an overdetermined system, i.e., solve the least
NB.      squares problem:
NB.        min || B - A*X ||
NB.   2) if m<n then gelsax finds the minimum norm solution
NB.      of an undetermined system
NB.        A * X = B
NB.   3) if m>=n then gelsacx finds the minimum norm solution
NB.      of an undetermined system
NB.        A^H * X = B
NB.   4) if m<n then gelsacx finds the least squares solution
NB.      of an overdetermined system, i.e., solve the least
NB.      squares problem:
NB.        min || B - A^H * X ||
NB.
NB. Syntax:
NB.   X1=. A1 gelsax  B1  NB. case 1: A1   * X1 = B1
NB.   X2=. A2 gelsax  B2  NB. case 2: A2   * X2 = B2
NB.   X2=. A1 gelsacx B2  NB. case 3: A1^H * X2 = B2
NB.   X1=. A2 gelsacx B1  NB. case 4: A2^H * X1 = B1
NB. where
NB.   A1 - m1×n1-matrix of full rank, m1>=n1, will be
NB.        factored by QR
NB.   B1 - m1-vector or m1×nrhs-matrix, the RHS
NB.   X1 - n1-vector or n1×nrhs-matrix, same rank as B1,
NB.        the solution(s)
NB.   A2 - m2×n2-matrix of full rank, m2<n2, will be
NB.        factored by LQ
NB.   B2 - m2-vector or m2×nrhs-matrix, the RHS
NB.   X2 - n2-vector or n2×nrhs-matrix, same rank as B2,
NB.        the solution(s)
NB.
NB. Algorithm:
NB.   In: A, B
NB.   Out: X
NB.   1) scale A and/or B
NB.   2) if (A == 0) then (X := 0) and return
NB.   3.1) if case 1 then:
NB.        3.1.1) compute QR factorization of A (Q * R = A)
NB.               without forming Q explicitly
NB.        3.1.2) multiply (Q^H * B) without forming Q
NB.               explicitly
NB.        3.1.3) extract m first rows from product (Q^H * B)
NB.               to get RHS
NB.        3.1.4) solve equation (R * X = RHS) for X
NB.   3.2) if case 2 then:
NB.        3.2.1) compute LQ factorization of A (L * Q = A)
NB.               without forming Q explicitly
NB.        3.2.2) solve equation (L * (Q*X) = B) for (Q*X)
NB.        3.2.3) expand the product (Q*X) by zeros from m to
NB.               n rows to get eQX
NB.        3.2.4) multiply (Q^H * eQX) without forming Q
NB.               explicitly to get X
NB.   3.3) if case 3 then:
NB.        3.3.1) compute QR factorization of A (Q * R = A)
NB.               without forming Q explicitly
NB.        3.3.2) solve the equation (R^H * (Q^H * X) = B)
NB.               for (Q^H * X)
NB.        3.3.3) expand the solution (Q^H * X) by zeros from
NB.               n to m rows to get eQhX
NB.        3.3.4) multiply (Q * eQhX) without forming Q
NB.               explicitly to get X
NB.   3.4) if case 4 then:
NB.        3.4.1) compute LQ factorization of A (L * Q = A)
NB.               without forming Q explicitly
NB.        3.4.2) multiply (Q * B) without forming Q
NB.               explicitly
NB.        3.4.3) shrink the product (Q * B) from n to m rows
NB.               to get RHS
NB.        3.4.4) solve the equation (L^H * X = RHS) for X
NB.   4) undo scaling
NB.
NB. Notes:
NB. - gelsax models LAPACK's xGELS('N')
NB. - gelsacx models LAPACK's DGELS('T') and ZGELS('C')

gelsax=: 4 : 0
  'smlnum bignum'=. FP_SFMIN (% , %~) FP_PREC
  segs=. 0 , (smlnum * 1 - FP_EPS) , bignum
  'm n'=. $ x
  xsh=. n 0} $ y
  anrm=. normm x
  bnrm=. normm y
  NB. verbs to [un]scale array norm if its max element is outside range [smlnum,bignum]
  scla=. [:`((anrm , smlnum)&scl)`]`((anrm , bignum)&scl)@.(segs I. anrm)
  sclb=. ] `((bnrm , smlnum)&scl)`]`((bnrm , bignum)&scl)@.(segs I. bnrm)
  if. m < n do.
    NB. case 2
    y=. x sclb^:_1@scla@((([ unmlqlc n {. (trtrslnn~ m&({."1))~) sclb)~ gelqf@scla)~ :: (xsh $ 0:) y
  else.
    NB. case 1
    y=. x sclb^:_1@scla@((([ trtrsunn&(n&{.) unmqrlc) sclb)~ geqrf@scla)~ :: (xsh $ 0:) y
  end.
)

gelsacx=: 4 : 0
  'smlnum bignum'=. FP_SFMIN (% , %~) FP_PREC
  segs=. 0 , (smlnum * 1 - FP_EPS) , bignum
  'm n'=. $ x
  xsh=. m 0} $ y
  anrm=. normm x
  bnrm=. normm y
  NB. verbs to [un]scale array norm if its max element is outside range [smlnum,bignum]
  scla=. [:`((anrm , smlnum)&scl)`]`((anrm , bignum)&scl)@.(segs I. anrm)
  sclb=. ] `((bnrm , smlnum)&scl)`]`((bnrm , bignum)&scl)@.(segs I. bnrm)
  if. m < n do.
    NB. case 4
    y=. x sclb^:_1@scla@((((m ({."1) [) trtrslcn m {. unmlqln) sclb)~ gelqf@scla)~ :: (xsh $ 0:) y
  else.
    NB. case 3
    y=. x sclb^:_1@scla@((([ unmqrln m {. (trtrsucn~ n& {.   )~) sclb)~ geqrf@scla)~ :: (xsh $ 0:) y
  end.
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgels1
NB.
NB. Description:
NB.   Test gelsxxx by general matrix and single RHS
NB.
NB. Syntax:
NB.   log=. testgels1 (A ; x)
NB. where
NB.   A   - m×n-matrix
NB.   x   - max(m,n)-vector, a material for the exact
NB.         solution values
NB.   log - 6-vector of boxes, test log

testgels1=: 3 : 0
  'A x'=. y

  rcondA=. nan`gecon1@.(=/@$) A  NB. meaninigful for square matrices only
  'norm1A normiA'=. (norm1 , normi) A
  'm n'=. $ A
  if. m = # x do.
    NB. A is tall
    xn=. n {. xc=. x
  else.
    NB. A is wide
    xc=. m {. xn=. x
  end.

  log=.          ('gelsax'  tdyad ((0&{::)`(1&{::)`]`(2&{::)`nan`(( mp~     qrt16v) >. ( mp~     qrt171)`qrt14@.(</@$@(0 {:: [))))) A ; (    A  mp xn) ; rcondA ; norm1A
  log=. log lcat ('gelsacx' tdyad ((0&{::)`(1&{::)`]`(2&{::)`nan`(((mp~ ct) qrt16v) >. qrt14`((mp~ ct) qrt171)@.(</@$@(0 {:: [))))) A ; ((ct A) mp xc) ; rcondA ; normiA
)

NB. ---------------------------------------------------------
NB. testgels3
NB.
NB. Description:
NB.   Test:
NB.   - xGELS (math/lapack2 addon)
NB.   - gelsxxx (math/mt addon)
NB.   by general matrix and multiple RHS
NB.
NB. Syntax:
NB.   log=. testgels3 (A ; X)
NB. where
NB.   A   - m×n-matrix
NB.   X   - max(m,n)×3-matrix, a material for exact solutions
NB.         values
NB.   log - 6-vector of boxes, test log

testgels3=: 3 : 0
  _1 cocreate < 'mttmp'
  load_mttmp_ 'math/mt/external/lapack2/gels'

  'A X'=. y
  rcondA=. nan`gecon1@.(=/@$) A  NB. meaninigful for square matrices only
  'norm1A normiA'=. (norm1 , normi) A
  'm n'=. $ A
  if. m = # X do.
    NB. A is tall
    Xn=. n {. Xc=. X
  else.
    NB. A is wide
    Xc=. m {. Xn=. X
  end.
  Bax=.      A  mp Xn
  Bacx=. (ct A) mp Xc

  log=.          ('dgels_mttmp_' tdyad (('n'"_)`(2&{. )`]`(2&{::)`nan`(( mp~     qrt16m) >. ( mp~     qrt171)`qrt14@.(</@$@(0 {:: [))))) A ; Bax  ; rcondA ; norm1A
  log=. log lcat ('dgels_mttmp_' tdyad (('t'"_)`(2&{. )`]`(2&{::)`nan`(((mp~ ct) qrt16m) >. qrt14`((mp~ ct) qrt171)@.(</@$@(0 {:: [))))) A ; Bacx ; rcondA ; normiA
  log=. log lcat ('zgels_mttmp_' tdyad (('n'"_)`(2&{. )`]`(2&{::)`nan`(( mp~     qrt16m) >. ( mp~     qrt171)`qrt14@.(</@$@(0 {:: [))))) A ; Bax  ; rcondA ; norm1A
  log=. log lcat ('zgels_mttmp_' tdyad (('c'"_)`(2&{. )`]`(2&{::)`nan`(((mp~ ct) qrt16m) >. qrt14`((mp~ ct) qrt171)@.(</@$@(0 {:: [))))) A ; Bacx ; rcondA ; normiA

  log=. log lcat ('gelsax'       tdyad ((0&{::)`(1&{::)`]`(2&{::)`nan`(( mp~     qrt16m) >. ( mp~     qrt171)`qrt14@.(</@$@(0 {:: [))))) A ; Bax  ; rcondA ; norm1A
  log=. log lcat ('gelsacx'      tdyad ((0&{::)`(1&{::)`]`(2&{::)`nan`(((mp~ ct) qrt16m) >. qrt14`((mp~ ct) qrt171)@.(</@$@(0 {:: [))))) A ; Bacx ; rcondA ; normiA

  coerase < 'mttmp'

  log
)

NB. ---------------------------------------------------------
NB. testgels
NB.
NB. Description:
NB.   Adv. to make verb to test gelsxxx by a matrix generated
NB.   by generator given
NB.
NB. Syntax:
NB.   log=. (mkmat testgels) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of non-negative integers, the shape of
NB.           matrix mat
NB.   log   - 6-vector of boxes, test log

testgels=: 1 : '(testgels3_mt_@; lcat_mt_~ testgels1_mt_@(; {:"1))~ u'

NB. ---------------------------------------------------------
NB. testls
NB.
NB. Description:
NB.   Adv. to make verb to test xxlsxxx by matrix of
NB.   datatype and shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testls) (m,n)
NB. where
NB.   mkmat - monad to generate a datatype; is called as:
NB.             dat=. mkmat 1
NB.           and is used to detect dat datatype only
NB.   (m,n) - 2-vector of integers, the shape of matrix to
NB.           generate
NB.   log   - 6-vector of boxes, test log
NB.
NB. Application:
NB. - test by random real wide matrices:
NB.     log=. ?@$&0 testls_mt_ 100 150
NB. - test by random complex tall matrices:
NB.     log=. (gemat_mt_ j. gemat_mt_) testls_mt_ 150 100
NB.
NB. Notes:
NB. - nrhs=3 is assumed
NB. - models part of LAPACK's xDRVLS which tests xGELS

testls=: 1 : 'lcat_mt_@((((dqrt131_mt_ testgels_mt_)`(dqrt132_mt_ testgels_mt_)`(dqrt133_mt_ testgels_mt_)`:0)~ (dlarnv2_mt_ % 3:)@(>./ , 3:))`(((zqrt131_mt_ testgels_mt_)`(zqrt132_mt_ testgels_mt_)`(zqrt133_mt_ testgels_mt_)`:0)~ (zlarnv2_mt_ % 3:)@(>./ , 3:))@.(JCMPX = 3!:0@u@1))'
