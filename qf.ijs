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

QFNB=: 32   NB. block size limit
QFNX=: 128  NB. crossover point, QFNX ≥ QFNB

NB. ---------------------------------------------------------
NB. gelq2
NB. geql2
NB. geqr2
NB. gerq2
NB.
NB. Description:
NB.   Q-factorization of the augmented input matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   LQf=. gelq2 (A ,. 0)
NB.   QfL=. geql2 (0 ,  A)
NB.   QfR=. geqr2 (A ,  0)
NB.   RQf=. gerq2 (0 ,. A)
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
NB.   1) if min(adjustedh,adjustedw)=0 then return eA
NB.   2) extract first/last row/column vector from eA and
NB.      reflect it:
NB.        z=. (larfxxx @ (ios & {)) eA
NB.   3) replace β in z by 1:
NB.        z=. (1 & (io })) z
NB.   4) apply an elementary reflector (packed in z) to eA:
NB.        eAupd=. z larfxxxx eA
NB.      (there is overhead here to re-calculate z) - CHECKME
NB.   5) extract two submatrices S0 and S1 from eAupd
NB.   6) supply S1 to recursive call to itself:
NB.        S1upd=. $: S1
NB.   8) combine output from z, S0 and S1upd
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - ge{lq,ql,qr,rq}2 emulates LAPACK's xGE{LQ,QL,QR,RQ}2
NB.   respectively
NB. - ge{lq,ql,qr,rq}2 and ge{lq,ql,qr,rq}3 respectively are
NB.   topologic equivalents
NB. - if triangular matrix diagonal's non-negativity is not
NB.   required, then larfp* may be replaced by faster larfg*

NB. ---------------------------------------------------------
NB. gelq3
NB. geql3
NB. geqr3
NB. gerq3
NB.
NB. Description:
NB.   Q-factorization of the augmented input matrix by
NB.   blocked version of algorithm
NB.
NB. Syntax:
NB.   LQf=. gelq3 (A ,. 0)
NB.   QfL=. geql3 (0 ,  A)
NB.   QfR=. geqr3 (A ,  0)
NB.   RQf=. gerq3 (0 ,. A)
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
NB.   1) if min(adjustedh,adjustedw)=0 then return eA
NB.   2) extract first/last row/column block from eA and
NB.      reflect it:
NB.        Vtau=. (gexx2 @ (ios & {)) eA
NB.   3) make Vtau unit triangular:
NB.        Vtau=. [diag] trx vTau
NB.   4) generate a block reflector from Vtau and apply it to
NB.      eA without block extracted:
NB.        eAupd=. Vtau larfbxxxx eA
NB.   5) extract two submatrices S0 and S1 from eAupd
NB.   6) supply S1 to recursive call to itself:
NB.        S1upd=. $: S1
NB.   8) combine output from Vtau, S0 and S1upd
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - ge{lq,ql,qr,rq}2 and ge{lq,ql,qr,rq}3 respectively are
NB.   topologic equivalents

NB. =========================================================
NB. Interface

NB. ####### TODO ###############
NB. If:
NB.   2 -: # $ A
NB. then (with appropriate comparison tolerance)
NB.   (] -: clean @ ((         trl   @( 0 _1&}.)) mp  unglq)@gelqf) A
NB.   (] -: clean @ ((((-~/@$) trl ])@( 1  0&}.)) mp~ ungql)@geqlf) A
NB.   (] -: clean @ ((         tru   @(_1  0&}.)) mp~ ungqr)@geqrf) A
NB.   (] -: clean @ ((((-~/@$) tru ])@( 0  1&}.)) mp  ungrq)@gerqf) A

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
NB. then (with appropriate comparison tolerance)
NB.   Q -: unglq (k {. LQf)
NB.   I -: (mp ct) Q
NB.   A -: L mp Q
NB.   (-: (((trl @ (0 _1 & }.)) mp unglq) @ gelqf)) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGELQF
NB.
NB. TODO:
NB. - consider >. instead of <.

qfi=: (<.@(QFNB %~ (_1+QFNB-QFNX)&+))M.

NB. usage: (iv4C`ih4C`iv4z`ih4z geqfios)
geqfios=: 1 : '_2 <@,\ (2{.m)/. , (2}.m)/. , ,.'

gelq2ios=: ((        }.&.>)`]`(        {.&.>)`]) geqfios
gelqfios=: ((  QFNB &}.&.>)`]`(  QFNB &{.&.>)`]) geqfios

geql2ios=: (]`(        }:&.>)`]`(        {:&.>)) geqfios
geqlfios=: (]`((-QFNB)&}.&.>)`]`((-QFNB)&{.&.>)) geqfios

geqr2ios=: (]`(        }.&.>)`]`(        {.&.>)) geqfios
geqrfios=: (]`(  QFNB &}.&.>)`]`(  QFNB &{.&.>)) geqfios

gerq2ios=: ((        }:&.>)`]`(        {:&.>)`]) geqfios
gerqfios=: (((-QFNB)&}.&.>)`]`((-QFNB)&{.&.>)`]) geqfios


NB. (] -: (clean@(trl mp unglq)@gelq2@(,.&0))) A
NB. (] -: (clean@(trl mp unglq)@gelqf)) A

gelq2step=: ((gelq2ios@}.) (((],(larfrnfr ~(1&(0})))) larfpfc) map2i 0 1 2) (0&({::))) ; (     }.&.>@}.)
gelq2=: 0 {:: (gelq2step ^: ((0 _1&(ms $))`(];((;&i.)/@$))))
gelqfstep=: ((gelqfios@}.) (((],(larfbrnfr~tru1    )) gelq2  ) map2i 0 1 2) (0&({::))) ; (QFNB&}.&.>@}.)
gelqf=: (gelq2`(((<@}.) (gelq2 upd1) (0&({::)))@(gelqfstep ^: ((qfi@(0 _1&(ms $)))`(];((;&i.)/@$)))))@.(*./@(QFNX<(0 _1+$))))@(,. & 0)


NB. (] -: (clean@((trl~ (-~/ @ $)) mp~ ungql)@geql2@(,~&0))) A
NB. (] -: (clean@((trl~ (-~/ @ $)) mp~ ungql)@geqlf)) A

geql2step=: ((geql2ios@}.) (((],.~(larflcbc ~(1&(_1})))) larfpb ) map2i 0 1 2) (0&({::))) ; (     }:&.>@}.)
geql2=: 0 {:: (geql2step ^: ((_1 0&(ms $))`(];((;&i.)/@$))))
geqlfstep=: ((geqlfios@}.) (((],.~(larfblcbc~(tru1~(-~/@$)))) geql2  ) map2i 0 1 2) (0&({::))) ; ((-QFNB)&}.&.>@}.)
geqlf=: (geql2`(((<@}.) (geql2 upd1) (0&({::)))@(geqlfstep ^: ((qfi@(_1 0&(ms $)))`(];((;&i.)/@$)))))@.(*./@(QFNX<(_1 0+$))))@(,~ & 0)


NB. (] -: (clean@(tru mp~ ungqr)@geqr2@(,&0))) A
NB. (] -: (clean@(tru mp~ ungqr)@geqrf)) A

geqr2step=: ((geqr2ios@}.) (((],.(larflcfc ~(1&(0})))) larfpf ) map2i 0 1 2) (0&({::))) ; (     }.&.>@}.)
geqr2=: 0 {:: (geqr2step ^: ((_1 0&(ms $))`(];((;&i.)/@$))))
geqrfstep=: ((geqrfios@}.) (((],.(larfblcfc~trl1    )) geqr2  ) map2i 0 1 2) (0&({::))) ; (QFNB&}.&.>@}.)
geqrf=: (geqr2`(((<@}.) (geqr2 upd1) (0&({::)))@(geqrfstep ^: ((qfi@(_1 0&(ms $)))`(];((;&i.)/@$)))))@.(*./@(QFNX<(_1 0+$))))@(, & 0)


NB. (] -: (clean@((tru~ (-~/ @ $)) mp ungrq)@gerq2@(,.~&0))) A
NB. (] -: (clean@((tru~ (-~/ @ $)) mp ungrq)@gerqf)) A

gerq2step=: ((gerq2ios@}.) (((],~(larfrnbr ~(1&(_1})))) larfpbc) map2i 0 1 2) (0&({::))) ; (     }:&.>@}.)
gerq2=: 0 {:: (gerq2step ^: ((0 _1&(ms $))`(];((;&i.)/@$))))
gerqfstep=: ((gerqfios@}.) (((],~(larfbrnbr~(trl1~(-~/@$)))) gerq2  ) map2i 0 1 2) (0&({::))) ; ((-QFNB)&}.&.>@}.)
gerqf=: (gerq2`(((<@}.) (gerq2 upd1) (0&({::)))@(gerqfstep ^: ((qfi@(0 _1&(ms $)))`(];((;&i.)/@$)))))@.(*./@(QFNX<(0 _1+$))))@(,.~ & 0)

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
NB. then (with appropriate comparison tolerance)
NB.   Q -: ungql (((m+1),(-n)) {. QfL)
NB.   I -: (mp~ ct) Q
NB.   A -: Q mp L
NB.   (-: ((ungql mp (((n - m) & trl) @ }.)) @ geqlf)) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGEQLF


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
NB. then (with appropriate comparison tolerance)
NB.   Q -: ungql (((m+1),n) {. QfR)
NB.   I -: (mp~ ct) Q
NB.   A -: Q mp R
NB.   (-: ((ungqr mp (tru @ }:)) @ geqrf)) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGEQRF


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
NB. then (with appropriate comparison tolerance)
NB.   Q -: ungrq (((-k),(n+1)) {. RQf)
NB.   I -: (mp ct) Q
NB.   A -: R mp Q
NB.   (-: (((((n - m) & tru) @ (0 1 & }.)) mp ungrq) @ gerqf)) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGERQF


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

tgeqf=: 3 : 0
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'gelqf geqlf geqrf gerqf'

  rcond=. ((_."_)`(norm1 con getri) @. (=/@$)) y

  ('128!:0' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(- (mp & >/)))%(FP_EPS*(#*norm1)@[)))) y

  ('2b1110 & gelqf_jlapack_' tmonad (]`({. , (,.  &. > / @ }.))`(rcond"_)`(_."_)`((norm1@(- ((mp  unglq) & > /)))%((FP_EPS*c*norm1)@[)))) y
  ('2b0111 & geqlf_jlapack_' tmonad (]`({: , (,~  &. > / @ }:))`(rcond"_)`(_."_)`((norm1@(- ((mp~ ungql) & > /)))%((FP_EPS*#*norm1)@[)))) y
  ('2b0111 & geqrf_jlapack_' tmonad (]`({: , (,   &. > / @ }:))`(rcond"_)`(_."_)`((norm1@(- ((mp~ ungqr) & > /)))%((FP_EPS*#*norm1)@[)))) y
  ('2b1110 & gerqf_jlapack_' tmonad (]`({. , (,.~ &. > / @ }.))`(rcond"_)`(_."_)`((norm1@(- ((mp  ungrq) & > /)))%((FP_EPS*c*norm1)@[)))) y

  ('gelqf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(( trl         @:(}:"1)@]) - unmlqrc~))%((FP_EPS*c*norm1)@[)))) y  NB. berr := ||A-L*Q||/(ε*n*||A||)
  ('geqlf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(((trl~(-~/@$))@  }.   @]) - unmqllc~))%((FP_EPS*#*norm1)@[)))) y  NB. berr := ||A-Q*L||/(ε*m*||A||)
  ('geqrf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(( tru         @  }:   @]) - unmqrlc~))%((FP_EPS*#*norm1)@[)))) y  NB. berr := ||A-Q*R||/(ε*m*||A||)
  ('gerqf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(((tru~(-~/@$))@:(}."1)@]) - unmrqrc~))%((FP_EPS*c*norm1)@[)))) y  NB. berr := ||A-R*Q||/(ε*n*||A||)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testqf
NB.
NB. Description:
NB.   Test orthogonal factorization algorithms by general
NB.   matrix of given size
NB.
NB. Syntax:
NB.   mkge testqf (m,n)
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testqf 150 100

testqf=: 1 : 'EMPTY [ tgeqf @ u'
