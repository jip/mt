NB. Condition number
NB.
NB. con        Conj. to make monad estimating the reciprocal
NB.            of the condition number of a matrix in a given
NB.            norm
NB. xxconx     Calculate reciprocal of the condition number
NB.            of a matrix in a given norm
NB. laic1x     Apply one step of incremental condition
NB.            estimation
NB.
NB. testcontr  Test reciprocal of the condition number
NB.            computing verbs by triangular matrix
NB. testconge  Test reciprocal of the condition number
NB.            computing verbs by square matrix
NB. testconhe  Test reciprocal of the condition number
NB.            computing verbs by Hermitian (symmetric)
NB.            matrix
NB. testconpo  Test reciprocal of the condition number
NB.            computing verbs by Hermitian (symmetric)
NB.            positive definite matrix
NB. testconpt  Test reciprocal of the condition number
NB.            computing verbs by Hermitian (symmetric)
NB.            positive definite tridiagonal matrix
NB. testconun  Test reciprocal of the condition number
NB.            computing verbs by unitary (orthogonal) matrix
NB. testcon    Adv. to make verb to test reciprocal of the
NB.            condition number computing algorithms by
NB.            matrix of generator and shape given
NB.
NB. Copyright 2010,2011,2013,2017,2018,2020,2021,2023,2024
NB.           Igor Zhuravlov
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
NB. stardot
NB. Extend monad *. for quaternions

stardot=: ;/@*.`(qnlen ; qnsign)@.(2 = #)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. con
NB.
NB. Description:
NB.   Conj. to make monad estimating the reciprocal of the
NB.   condition number of a matrix in a given norm
NB.
NB. Syntax:
NB.   rcondA=. (norm con inv) A
NB. where
NB.   norm - monad to calculate norm of matrix, is called as:
NB.            normA=. norm A
NB.   inv  - monad to inverse square matrix, is called as:
NB.            invA=. inv A
NB.   A    - n×n-matrix
NB.
NB. TODO:
NB. - implement more practical norm-estimation approach

con=: 2 : '*&(%@u) v'

NB. ---------------------------------------------------------
NB. trl1con1
NB. trl1coni
NB. trlcon1
NB. trlconi
NB. tru1con1
NB. tru1coni
NB. trucon1
NB. truconi
NB. gecon1
NB. geconi
NB. hecon1
NB. heconi
NB. pocon1
NB. poconi
NB. ptcon1
NB. ptconi
NB. uncon1
NB. unconi
NB.
NB. Description:
NB.   Calculate reciprocal of the condition number of a
NB.   matrix in a given norm
NB.
NB. Syntax:
NB.   rcondR=. trxxconx R
NB.   rcondG=. geconx G
NB.   rcondH=. heconx H
NB.   rcondP=. poconx P
NB.   rcondT=. ptconx T
NB.   rcondQ=. unconx Q
NB. where
NB.   R      - n×n-matrix of type: triangular
NB.   rcondR ≥ 0, reciprocal of the condition number of R
NB.   G      - n×n-matrix of type: general, band,
NB.            tridiagonal, triangular or triangular band
NB.   rcondG ≥ 0, reciprocal of the condition number of G
NB.   H      - n×n-matrix of type: Hermitian (symmetric)
NB.   rcondH ≥ 0, reciprocal of the condition number of H
NB.   P      - n×n-matrix of type: Hermitian (symmetric)
NB.            positive definite
NB.   rcondP ≥ 0, reciprocal of the condition number of P
NB.   T      - n×n-matrix of type: Hermitian (symmetric)
NB.            positive definite tridiagonal
NB.   rcondT ≥ 0, reciprocal of the condition number of T
NB.   Q      - n×n-matrix, the unitary (orthogonal)
NB.   rcondQ ≥ 0, reciprocal of the condition number of Q in
NB.            1-norm
NB.
NB. Notes:
NB. - extraneous values in band, tridiagonal, triangular and
NB.   triangular band matrices must be zeroed
NB. - trl1con1 simulates LAPACK's xTRCON('1','L','U')
NB. - trl1coni simulates LAPACK's xTRCON('i','L','U')
NB. - trlcon1 simulates LAPACK's xTRCON('1','L','N')
NB. - trlconi simulates LAPACK's xTRCON('i','L','N')
NB. - tru1con1 simulates LAPACK's xTRCON('1','U','U')
NB. - tru1coni simulates LAPACK's xTRCON('i','U','U')
NB. - trucon1 simulates LAPACK's xTRCON('1','U','N')
NB. - truconi simulates LAPACK's xTRCON('i','U','N')
NB. - gecon1 simulates LAPACK's xGECON('1'), xGBCON('1'),
NB.   xGTCON('1'), xTBCON('1')
NB. - geconi simulates LAPACK's xGECON('i'), xGBCON('i'),
NB.   xGTCON('i'), xTBCON('i')
NB. - hecon1 simulates LAPACK's DSYCON('1'), ZHECON('1')
NB. - heconi simulates LAPACK's DSYCON('i'), ZHECON('i')
NB. - pocon1 simulates LAPACK's xPBCON('1'), xPOCON('1')
NB. - poconi simulates LAPACK's xPBCON('i'), xPOCON('i')
NB. - ptcon1 simulates LAPACK's xPTCON('1')
NB. - ptconi simulates LAPACK's xPTCON('i')

trl1con1=: 1:`(norm1 con trtril1 :: 0)@.(*@#)
trl1coni=: 1:`(normi con trtril1 :: 0)@.(*@#)
trlcon1=:  1:`(norm1 con trtril  :: 0)@.(*@#)
trlconi=:  1:`(normi con trtril  :: 0)@.(*@#)
tru1con1=: 1:`(norm1 con trtriu1 :: 0)@.(*@#)
tru1coni=: 1:`(normi con trtriu1 :: 0)@.(*@#)
trucon1=:  1:`(norm1 con trtriu  :: 0)@.(*@#)
truconi=:  1:`(normi con trtriu  :: 0)@.(*@#)

gecon1=: 1:`(norm1 con (getrilu1p@getrflu1p) :: 0)@.(*@#)
geconi=: 1:`(normi con (getrilu1p@getrflu1p) :: 0)@.(*@#)

hecon1=: 1:`(norm1 con (hetripl@hetrfpl) :: 0)@.(*@#)
heconi=: 1:`(normi con (hetripl@hetrfpl) :: 0)@.(*@#)

pocon1=: 1:`(norm1 con (potril@potrfl) :: 0)@.(*@#)
poconi=: 1:`(normi con (potril@potrfl) :: 0)@.(*@#)

ptcon1=: 1:`(norm1 con pttril :: 0)@.(*@#)
ptconi=: 1:`(normi con pttril :: 0)@.(*@#)

uncon1=: 1:`(norm1 con ct :: 0)@.(*@#)
unconi=: 1:`(normi con ct :: 0)@.(*@#)

NB. ---------------------------------------------------------
NB. laic11
NB. laic12
NB.
NB. Description:
NB.   Apply one step of incremental condition estimation in
NB.   its simplest version. Let ix, twonorm(ix) = 1, be an
NB.   approximate singular vector of a lower triangular
NB.   j×j-matrix iL, such that
NB.     twonorm(iL*ix) = isest
NB.   Then laic1x computes osest, s, c such that the vector
NB.          [ s*ix ]
NB.     ox = [  c   ]
NB.   is an approximate singular vector of
NB.          [ iL     0  ]
NB.     oL = [ w'  gamma ]
NB.   in the sense that
NB.     twonorm(oL*ox) = osest.
NB.   Note that [s c]' and osest^2 is an eigenpair of the
NB.   system
NB.     diag(isest^2, 0) + [alpha gamma] * [ conjg(alpha) ]
NB.                                        [ conjg(gamma) ]
NB.   where alpha = conjg(ix)'*w.
NB.
NB. Syntax:
NB.   'osest cs'=. laic1x isest;ga
NB. where
NB.   isest -: norms iL mp ix
NB.   ga    -: gamma , alpha
NB.   osest -: norms oL mp ox
NB.   cs    -: c , s
NB.
NB. Notes:
NB. - laic11 models LAPACK's xLAIC1(1) and computes largest
NB.   singular value
NB. - laic12 models LAPACK's xLAIC1(2) and computes smallest
NB.   singular value

laic11=: 3 : 0
  'absest absga'=. (|L:0) 'sest ga'=. y
  'absg absa'=. absga
  NB. special cases
  if. 0 = sest do.
    if. (>./) absga do.
      stardot ga
    else.
      0 ; 1 0
    end.
  elseif. absg <: FP_EPS * absest do.
    (| absest j. absa);0 1
  elseif. absa <: FP_EPS * absest do.
    if. absg <: absest do.
      absest ; 0 1
    else.
      absg ; 1 0
    end.
  elseif. (+./) absest <: FP_EPS * absga do.
    stardot ga
  else.
    NB. normal case
    b=. -: 1 - +/ 't c'=. *: absga % absest
    if. 0 < b do.
      t=. c % b + %: c + *: b
    else.
      t=. (%: c + *: b) - b
    end.
    (absest * %: >: t) ; qnsign (ga % absest) % (,~ <:) -t
  end.
)

laic12=: 3 : 0
  'sest ga'=. y
  absest=. | sest
  'absa absg'=. absag=. | ag=. |. ga
  NB. special cases
  if. 0 = sest do.
    if. (>./) absag do.
      0 ; qnsign qnconij ag
    else.
      0 ; 0 1
    end.
  elseif. absg <: FP_EPS * absest do.
    absg ; 1 0
  elseif. absa <: FP_EPS * absest do.
    if. absg <: absest do.
      absg ; 1 0
    else.
      absest ; 0 1
    end.
  elseif. (+./) absest <: FP_EPS * absag do.
    scl=. (>:&.*:) tmp=. %/ 'm M'=. absa (<. , >.) absg
    cs=. scl %~ M %~ qnconij ag
    if. absg > absa do.
      (absest % scl) ; cs
    else.
      (absest * tmp % scl) ; cs
    end.
  else.
    NB. normal case
    norma=. (>:@(+/)@(* {.) >. +/@(* {:)) zeta=. absag % absest
    if. 0 <: >: +: ((- * +)/) zeta do.
      b=. -: >: +/ 't c'=. *: zeta
      t=. c % b + %: | (*: b) - c
      (absest * %: t + 4 * norma * *: FP_EPS) ; qnsign (ga % absest) % (, >:) -t
    else.
      b=. -: <: +/ 'c t'=. *: zeta
      if. 0 <: b do.
        t=. - c % b + %: c + *: b
      else.
        t=. b - %: c + *: b
      end.
      (absest * %: >: t + 4 * norma * *: FP_EPS) ; qnsign (ga % absest) % (,~ <:) -t
    end.
  end.
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testcontr
NB.
NB. Description:
NB.   Test reciprocal of the condition number computing
NB.   verbs:
NB.   - trl1con1 (math/mt addon)
NB.   - trl1coni (math/mt addon)
NB.   - trlcon1 (math/mt addon)
NB.   - trlconi (math/mt addon)
NB.   - tru1con1 (math/mt addon)
NB.   - tru1coni (math/mt addon)
NB.   - trucon1 (math/mt addon)
NB.   - truconi (math/mt addon)
NB.   by triangular matrix
NB.
NB. Syntax:
NB.   log=. testcontr G
NB. where
NB.   G   - n×n-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs

testcontr=: 3 : 0
  log=.          ('trl1con1' tmonad (trl1pick`]`(trl1con1@trl1pick)`nan`nan)) y
  log=. log lcat ('trl1coni' tmonad (trl1pick`]`(trl1coni@trl1pick)`nan`nan)) y
  log=. log lcat ('trlcon1'  tmonad (trlpick `]`(trlcon1 @trlpick )`nan`nan)) y
  log=. log lcat ('trlconi'  tmonad (trlpick `]`(trlconi @trlpick )`nan`nan)) y
  log=. log lcat ('tru1con1' tmonad (tru1pick`]`(tru1con1@tru1pick)`nan`nan)) y
  log=. log lcat ('tru1coni' tmonad (tru1pick`]`(tru1coni@tru1pick)`nan`nan)) y
  log=. log lcat ('trucon1'  tmonad (trupick `]`(trucon1 @trupick )`nan`nan)) y
  log=. log lcat ('truconi'  tmonad (trupick `]`(truconi @trupick )`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testconge
NB.
NB. Description:
NB.   Test reciprocal of the condition number computing
NB.   verbs:
NB.   - gecon1 (math/mt addon)
NB.   - geconi (math/mt addon)
NB.   by square matrix
NB.
NB. Syntax:
NB.   log=. testconge G
NB. where
NB.   G   - n×n-matrix
NB.   log - 6-vector of boxes, test log, see test.ijs

testconge=: 3 : 0
  rcond=. geconi y

  log=.          ('gecon1' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('geconi' tmonad (]`]`(rcond"_)`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testconhe
NB.
NB. Description:
NB.   Test reciprocal of the condition number computing
NB.   verbs:
NB.   - hecon1 (math/mt addon)
NB.   - heconi (math/mt addon)
NB.   by Hermitian (symmetric) matrix
NB.
NB. Syntax:
NB.   log=. testconhe H
NB. where
NB.   H   - n×n-matrix, Hermitian (symmetric)
NB.   log - 6-vector of boxes, test log, see test.ijs

testconhe=: 3 : 0
  rcond=. heconi y

  log=.          ('hecon1' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('heconi' tmonad (]`]`(rcond"_)`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testconpo
NB.
NB. Description:
NB.   Test reciprocal of the condition number computing
NB.   verbs:
NB.   - pocon1 (math/mt addon)
NB.   - poconi (math/mt addon)
NB.   by Hermitian (symmetric) positive definite matrix
NB.
NB. Syntax:
NB.   log=. testconpo P
NB. where
NB.   P   - n×n-matrix, Hermitian (symmetric) positive
NB.         definite
NB.   log - 6-vector of boxes, test log, see test.ijs

testconpo=: 3 : 0
  rcond=. poconi y

  log=.          ('pocon1' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('poconi' tmonad (]`]`(rcond"_)`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testconpt
NB.
NB. Description:
NB.   Test reciprocal of the condition number computing
NB.   verbs:
NB.   - ptcon1 (math/mt addon)
NB.   - ptconi (math/mt addon)
NB.   by Hermitian (symmetric) positive definite tridiagonal
NB.   matrix
NB.
NB. Syntax:
NB.   log=. testconpt T
NB. where
NB.   T   - n×n-matrix, Hermitian (symmetric) positive
NB.         definite tridiagonal
NB.   log - 6-vector of boxes, test log, see test.ijs

testconpt=: 3 : 0
  rcond=. ptconi y

  log=.          ('ptcon1' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('ptconi' tmonad (]`]`(rcond"_)`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testconun
NB.
NB. Description:
NB.   Test reciprocal of the condition number computing
NB.   verbs:
NB.   - ptcon1 (math/mt addon)
NB.   - ptconi (math/mt addon)
NB.   by unitary (orthogonal) matrix
NB.
NB. Syntax:
NB.   log=. testconun Q
NB. where
NB.   Q   - n×n-matrix, unitary (orthogonal)
NB.   log - 6-vector of boxes, test log, see test.ijs

testconun=: 3 : 0
  rcond=. unconi y

  log=.          ('uncon1' tmonad (]`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('unconi' tmonad (]`]`(rcond"_)`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testcon
NB.
NB. Description:
NB.   Adv. to make verb to test reciprocal of the condition
NB.   number computing algorithms by matrix of generator and
NB.   shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testcon) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testcon_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testcon_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testcon_mt_ 150 200

testcon=: 1 : 'nolog_mt_`(testconun_mt_@((randnr_mt_ unmat_mt_)`(randnc_mt_ unmat_mt_)@.(JCMPX = (3!:0)@u@1:)) ,&.>~ testconpt_mt_@(u ptmat2_mt_) ,&.>~ testconpo_mt_@(u pomat_mt_) ,&.>~ testconhe_mt_@(u hemat_mt_) ,&.>~ (testconge_mt_ ,&.>~ testcontr_mt_)@u)@.(=/)'
