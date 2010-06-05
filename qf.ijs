NB. qf.ijs
NB. Orthogonal factorizations LQ QL QR RQ
NB.
NB. gelqf  LQ factorization of a general matrix
NB. geqlf  QL factorization of a general matrix
NB. geqrf  QR factorization of a general matrix
NB. gerqf  RQ factorization of a general matrix
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

GEQFBS=: 32  NB. block size limit

NB. ---------------------------------------------------------
NB. Algorithm    Non-blocked version    Blocked version
NB. LQ           gelq2                  gelq3
NB. QL           geql2                  geql3
NB. QR           geqr2                  geqr3
NB. RQ           gerq2                  gerq3
NB.
NB. Description:
NB.   Q-factorization of the augmented input matrix
NB.
NB. Syntax:
NB.   LQf=. gelqx (A ,. 0)
NB.   QfL=. geqlx (0 ,  A)
NB.   QfR=. geqrx (A ,  0)
NB.   RQf=. gerqx (0 ,. A)
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   LQf - m×(n+1)-matrix, combined lower triangular
NB.         m×k-matrix L and unit upper triangular
NB.         k×(n+1)-matrix Qf (unit diagonal not stored)
NB.   QfL - (m+1)×n-matrix, combined unit upper triangular
NB.         (m+1)×k-matrix Qf (unit diagonal not stored) and
NB.         lower triangular k×n-matrix L
NB.   QfR - (m+1)×n-matrix, combined unit lower triangular
NB.         (m+1)×k-matrix Qf (unit diagonal not stored) and
NB.         upper triangular k×n-matrix R
NB.   RQf - m×(n+1)-matrix, combined upper triangular
NB.         m×k-matrix R and unit lower triangular
NB.         k×(n+1)-matrix Qf (unit diagonal not stored)
NB.   Qf  - the matrix Q represented in factored form
NB.   Q   - matrix with orthonormal rows or columns which is
NB.         defined as the product of k elementary reflectors
NB.   k   = min(m,n)
NB.
NB. Algorithm:
NB.   Input: h×w-matrix eA
NB.   Output: h×w-matrix eAf
NB.   Note: gexx2 operates on vectors, gexx3 operates on
NB.         GEQFBS-wide blocks
NB.   1) if min(adjustedh,adjustedw)=0 then return eA
NB.   2) extract first/last row/column vector/block from eA
NB.      and reflect it (gexx2):
NB.        z=. (larfxxx @ (ios & {)) eA
NB.      or (gexx3):
NB.        z=. (gexx2 @ (ios & {)) eA
NB.   3) replace β in z by 1 (gexx2):
NB.        z=. (1 & (io })) z
NB.      or diagonal of βs by 1 (gexx3):
NB.        z=. [diag] trx z
NB.   4) apply an elementary reflector (packed in z) to eA
NB.      (gexx2):
NB.        eAupd=. z larfxxxx eA
NB.      (there is overhead here to re-calculate z), or
NB.      or generate a block reflector from z and apply it to
NB.      eA (gexx3):
NB.        eAupd=. z larfbxxxx eA
NB.   5) extract two submatrices S0 and S1 from eAupd
NB.   6) supply S1 to recursive call to itself:
NB.        S1upd=. $: S1
NB.   8) combine output from z, S0 and S1upd
NB.
NB. Notes:
NB. - ge{lq,ql,qr,rq}2 emulates LAPACK's xGE{LQ,QL,QR,RQ}2
NB.   respectively
NB. - input's and output's shapes are the same
NB. - ge{lq,ql,qr,rq}2 and ge{lq,ql,qr,rq}3 respectively are
NB.   topologic equivalents

gelq2=: ((] ,   (((( 1 ,~ (1 -  #)       ) {. ]) ,.  ($: @ ( 1  1 & }.))) @ (larfrnfr~ (1 & ( 0 }))))) (larfgfc @ (IOSFR & {))) ^: (*./ @ (0 < (0 _1 + $)))
geql2=: ((] ,.~ ((((_1 ,  (1 -~ ({: @ $))) {. ]) ,~  ($: @ (_1 _1 & }.))) @ (larflcbc~ (1 & (_1 }))))) (larfgb  @ (IOSLC & {))) ^: (*./ @ (0 < (_1 0 + $)))
geqr2=: ((] ,.  (((( 1 ,  (1 -  ({: @ $))) {. ]) ,   ($: @ ( 1  1 & }.))) @ (larflcfc~ (1 & ( 0 }))))) (larfgf  @ (IOSFC & {))) ^: (*./ @ (0 < (_1 0 + $)))
gerq2=: ((] ,~  ((((_1 ,~ (1 -~ #)       ) {. ]) ,.~ ($: @ (_1 _1 & }.))) @ (larfrnbr~ (1 & (_1 }))))) (larfgbc @ (IOSLR & {))) ^: (*./ @ (0 < (0 _1 + $)))

gelq3=: (((     GEQFBS    <.   #     ) }. ]) (] ,   ((((_ , (GEQFBS    <.  ({:@$))) {. ]) ,.  ($: @ ((0 ,   GEQFBS ) }. ]))) @ (larfbrnfr~          tru1   ))) (gelq2 @ ((     GEQFBS    <.  #      ) {. ]))) ^: (*./ @ (0 < (0 _1 + $)))
geql3=: (((0 , (GEQFBS (-@<.) ({:@$))) }. ]) (] ,.~ ((((     GEQFBS (-@<.)     #  ) {. ]) ,~  ($: @ (     (-GEQFBS)  }. ]))) @ (larfblcbc~ ((-~/@$) tru1 ])))) (geql2 @ ((_ , (GEQFBS (-@<.) ({:@$))) {. ]))) ^: (*./ @ (0 < (_1 0 + $)))
geqr3=: (((0 , (GEQFBS    <.  ({:@$))) }. ]) (] ,.  ((((     GEQFBS    <.      #  ) {. ]) ,   ($: @ (       GEQFBS   }. ]))) @ (larfblcfc~          trl1   ))) (geqr2 @ ((_ , (GEQFBS    <.  ({:@$))) {. ]))) ^: (*./ @ (0 < (_1 0 + $)))
gerq3=: (((     GEQFBS (-@<.)  #     ) }. ]) (] ,~  ((((_ , (GEQFBS (-@<.) ({:@$))) {. ]) ,.~ ($: @ ((0 , (-GEQFBS)) }. ]))) @ (larfbrnbr~ ((-~/@$) trl1 ])))) (gerq2 @ ((     GEQFBS (-@<.) #      ) {. ]))) ^: (*./ @ (0 < (0 _1 + $)))

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gelqf
NB.
NB. Description:
NB.   LQ factorization of a general matrix
NB.
NB. Syntax:
NB.   LQf=. gelqf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   LQf - m×(n+1)-matrix, combined lower triangular
NB.         m×k-matrix L and unit upper triangular
NB.         k×(n+1)-matrix Qf (unit diagonal not stored)
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'm n'=. $ A
NB.   k=. m <. n
NB.   LQf=. gelqf A
NB.   L=. trl (0 _1 }. LQf)
NB.   Q=. unglq LQf
NB. then
NB.   Q -: unglq (k {. LQf)
NB.   I -: (mp ct) Q
NB.   A -: L mp Q
NB.   (-: (((trl @ (0 _1 & }.)) mp unglq) @ gelqf)) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGELQF
NB.
NB. TODO:
NB. - non-negative diagonal

gelqf=: gelq3 @ (,. & 0)

NB. ---------------------------------------------------------
NB. geqlf
NB.
NB. Description:
NB.   QL factorization of a general matrix
NB.
NB. Syntax:
NB.   QfL=. geqlf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   QfL - (m+1)×n-matrix, combined unit upper triangular
NB.         (m+1)×k-matrix Qf (unit diagonal not stored) and
NB.         lower triangular k×n-matrix L
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'm n'=. $ A
NB.   k=. m <. n
NB.   QfL=. geqlf A
NB.   Q=. ungql QfL
NB.   L=. (n - m) trl (}. QfL)
NB. then
NB.   Q -: ungql (((m+1),(-n)) {. QfL)
NB.   I -: (mp~ ct) Q
NB.   A -: Q mp L
NB.   (-: ((ungql mp (((n - m) & trl) @ }.)) @ geqlf)) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGEQLF
NB.
NB. TODO:
NB. - non-negative diagonal

geqlf=: geql3 @ (0 & ,)

NB. ---------------------------------------------------------
NB. geqrf
NB.
NB. Description:
NB.   QR factorization of a general matrix
NB.
NB. Syntax:
NB.   QfR=. geqrf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   QfR - (m+1)×n-matrix, combined unit lower triangular
NB.         (m+1)×k-matrix Qf (unit diagonal not stored) and
NB.         upper triangular k×n-matrix R
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'm n'=. $ A
NB.   k=. m <. n
NB.   QfR=. geqrf A
NB.   Q=. ungqr QfR
NB.   R=. tru (}: QfR)
NB. then
NB.   Q -: ungql (((m+1),n) {. QfR)
NB.   I -: (mp~ ct) Q
NB.   A -: Q mp R
NB.   (-: ((ungqr mp (tru @ }:)) @ geqrf)) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGEQRF
NB.
NB. TODO:
NB. - non-negative diagonal

geqrf=: geqr3 @ (, & 0)

NB. ---------------------------------------------------------
NB. gerqf
NB.
NB. Description:
NB.   RQ factorization of a general matrix
NB.
NB. Syntax:
NB.   RQf=. gerqf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   RQf - m×(n+1)-matrix, combined upper triangular
NB.         m×k-matrix R and unit lower triangular
NB.         k×(n+1)-matrix Qf (unit diagonal not stored)
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'm n'=. $ A
NB.   k=. m <. n
NB.   RQf=. gerqf A
NB.   R=. (n - m) trl (0 1 }. RQf)
NB.   Q=. ungrq RQf
NB. then
NB.   Q -: ungrq (((-k),(n+1)) {. RQf)
NB.   I -: (mp ct) Q
NB.   A -: R mp Q
NB.   (-: (((((n - m) & tru) @ (0 1 & }.)) mp ungrq) @ gerqf)) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGERQF
NB.
NB. TODO:
NB. - non-negative diagonal

gerqf=: gerq3 @ (0 & ,.)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. tgeqf
NB. Test orthogonal factorization algorithms:
NB. - built-in: 128!:0
NB. - LAPACK: gelqf geqlf geqrf gerqf
NB. - mt package: gelqf geqlf geqrf gerqf
NB. by general matrix
NB.
NB. Syntax: tgeqf A
NB. where A - general m×n-matrix
NB.
NB. TODO:
NB. - split on tgeqf_{nonblocked,builtin,lapack,blocked}

tgeqf=: 3 : 0
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'gelqf geqlf geqrf gerqf'

  rcond=. ((_."_)`(norm1 con getri) @. (=/@$)) y

  ('128!:0' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(- (mp & >/)))%(FP_EPS*(#*norm1)@[)))) y

  ('2b1110 & gelqf_jlapack_' tmonad (]`({. , (,.  &. > / @ }.))`(rcond"_)`(_."_)`((norm1@(- ((mp  unglq) & > /)))%((FP_EPS*#*norm1)@[)))) y
  ('2b0111 & geqlf_jlapack_' tmonad (]`({: , (,~  &. > / @ }:))`(rcond"_)`(_."_)`((norm1@(- ((mp~ ungql) & > /)))%((FP_EPS*#*norm1)@[)))) y
  ('2b0111 & geqrf_jlapack_' tmonad (]`({: , (,   &. > / @ }:))`(rcond"_)`(_."_)`((norm1@(- ((mp~ ungqr) & > /)))%((FP_EPS*#*norm1)@[)))) y
  ('2b1110 & gerqf_jlapack_' tmonad (]`({. , (,.~ &. > / @ }.))`(rcond"_)`(_."_)`((norm1@(- ((mp  ungrq) & > /)))%((FP_EPS*#*norm1)@[)))) y

  ('gelqf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(- ((         trl   @( 0 _1&}.)) mp  unglq)))%((FP_EPS*#*norm1)@[)))) y  NB. berr := ||A-L*Q||/ε*m*||A||
  ('geqlf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(- ((((-~/@$) trl ])@( 1  0&}.)) mp~ ungql)))%((FP_EPS*#*norm1)@[)))) y  NB. berr := ||A-Q*L||/ε*m*||A||
  ('geqrf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(- ((         tru   @(_1  0&}.)) mp~ ungqr)))%((FP_EPS*#*norm1)@[)))) y  NB. berr := ||A-Q*R||/ε*m*||A||
  ('gerqf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(- ((((-~/@$) tru ])@( 0  1&}.)) mp  ungrq)))%((FP_EPS*#*norm1)@[)))) y  NB. berr := ||A-R*Q||/ε*m*||A||

  EMPTY
)

NB. ---------------------------------------------------------
NB. testqf
NB. Test orthogonal factorizations
NB. Syntax: mkge testqf (m,n)
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testqf 150 100

testqf=: 1 : 'EMPTY [ tgeqf @ u'
