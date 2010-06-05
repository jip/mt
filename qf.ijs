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
NB. qfi
NB.
NB. Description: Number of iterations
NB. Syntax:      iters=. ungi k
NB. where        k = min(rows,columns)
NB. Formula:     iters = max(0,⌊(k+NB-NX-1)/NB⌋)
NB. Notes:       is memo, since repetitive calls are expected

qfi=: (0 >. <.@(QFNB %~ (_1+QFNB-QFNX)&+))M.

NB. ---------------------------------------------------------
NB. gelq2
NB.
NB. Description:
NB.   LQ-factorization of the augmented input matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   LQf=. gelq2 eA
NB. where
NB.   eA  - m×(n+1)-matrix, being A augmented by trash vector
NB.   A   - m×n-matrix, the input to factorize
NB.   LQf - m×(n+1)-matrix, combined lower triangular
NB.         m×k-matrix L and unit upper triangular
NB.         k×(n+1)-matrix Qf (unit diagonal not stored)
NB.   Qf  - the matrix Q represented in factored form
NB.   Q   - matrix with orthonormal rows which is defined as
NB.         the product of k elementary reflectors
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   input:                output:
NB.   (  pfx sfxT  ) k      (  pfx   sfxT    ) k+nb
NB.   (  pfx sfxB  ) m-k    (  pfx   sfxB    ) m-k-nb
NB.      k   n+1-k             k+nb  n+1-k-nb
NB.
NB. Notes:
NB. - gelq2 emulates LAPACK's xGELQ2
NB. - gelq2 and gelqf are topologic equivalents
NB. - if triangular matrix diagonal's non-negativity is not
NB.   required, then larfp* may be replaced by faster larfg*

gelq2=: ((0&({::)) 0 stitch (1&({::))) @ ((3 : 0) ^: ((0 _1&(ms $))`((0&({."1));(0&{.);])))
  'pfx sfxT sfxB'=. y
  z=. larfpfc {. sfxB
  sfxT=. sfxT , z
  sfxB=. (1 (0}) z) larfrnfr }. sfxB
  (pfx ,. (sfxT (, &: ({."1)) sfxB)) ; (sfxT (; &: (}."1)) sfxB)
)

NB. ---------------------------------------------------------
NB. geql2
NB.
NB. Description:
NB.   QL-factorization of the augmented input matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   QfL=. geql2 eA
NB. where
NB.   eA  - (m+1)×n-matrix, being A augmented by trash vector
NB.   A   - m×n-matrix, the input to factorize
NB.   QfL - (m+1)×n-matrix, combined unit upper triangular
NB.         (m+1)×k-matrix Qf (unit diagonal not stored) and
NB.         lower triangular k×n-matrix L
NB.   Qf  - the matrix Q represented in factored form
NB.   Q   - matrix with orthonormal columns which is defined
NB.         as the product of k elementary reflectors
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   input:                   output:
NB.   (  pfxL pfxR  ) m+1-k    (  pfxL    pfxR  ) m+1-k-nb
NB.   (  sfx  sfx   ) k        (  sfx     sfx   ) k+nb
NB.      n-k  k                   n-k-nb  k+nb
NB.
NB. Notes:
NB. - geql2 emulates LAPACK's xGEQL2
NB. - geql2 and geqlf are topologic equivalents
NB. - if triangular matrix diagonal's non-negativity is not
NB.   required, then larfp* may be replaced by faster larfg*

geql2=: ((1&({::)) _1 append (2&({::))) @ ((3 : 0) ^: ((_1 0&(ms $))`(];(0&({."1));(0&{.))))
  'pfxL pfxR sfx'=. y
  z=. larfpb {:"1 pfxL
  pfxR=. z ,. pfxR
  pfxL=. (1 (_1}) z) larflcbc }:"1 pfxL
  (pfxL (; & }:) pfxR) , < ((pfxL (, & {:) pfxR) , sfx)
)

NB. ---------------------------------------------------------
NB. geqr2
NB.
NB. Description:
NB.   QR-factorization of the augmented input matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   QfR=. geqr2 eA
NB. where
NB.   eA  - (m+1)×n-matrix, being A augmented by trash vector
NB.   A   - m×n-matrix, the input to factorize
NB.   QfR - (m+1)×n-matrix, combined unit lower triangular
NB.         (m+1)×k-matrix Qf (unit diagonal not stored) and
NB.         upper triangular k×n-matrix R
NB.   Qf  - the matrix Q represented in factored form
NB.   Q   - matrix with orthonormal columns which is defined
NB.         as the product of k elementary reflectors
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   input:                    output:
NB.   (  pfx   pfx   ) k        (  pfx   pfx   ) k+nb
NB.   (  sfxL  sfxR  ) m+1-k    (  sfxL  sfxR  ) m+1-k-nb
NB.      k     n-k                 k+nb  n-k-nb
NB.
NB. Notes:
NB. - geqr2 emulates LAPACK's xGEQR2
NB. - gerq2 and geqrf are topologic equivalents
NB. - if triangular matrix diagonal's non-negativity is not
NB.   required, then larfp* may be replaced by faster larfg*

geqr2=: ((0&({::)) , (1&({::))) @ ((3 : 0) ^: ((_1 0&(ms $))`((0&{.);(0&({."1));])))
  'pfx sfxL sfxR'=. y
  z=. larfpf {."1 sfxR
  sfxL=. sfxL ,. z
  sfxR=. (1 (0}) z) larflcfc }."1 sfxR
  (pfx , (sfxL (, & {.) sfxR)) ; (sfxL (; & }.) sfxR)
)

NB. ---------------------------------------------------------
NB. gerq2
NB.
NB. Description:
NB.   RQ-factorization of the augmented input matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   RQf=. gerq2 eA
NB. where
NB.   eA  - m×(n+1)-matrix, being A augmented by trash vector
NB.   A   - m×n-matrix, the input to factorize
NB.   RQf - m×(n+1)-matrix, combined upper triangular
NB.         m×k-matrix R and unit lower triangular
NB.         k×(n+1)-matrix Qf (unit diagonal not stored)
NB.   Qf  - the matrix Q represented in factored form
NB.   Q   - matrix with orthonormal rows which is defined as
NB.         the product of k elementary reflectors
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   input:                  output:
NB.   (  pfxT   sfx  ) m-k    (  pfxT      sfx  ) m-k-nb
NB.   (  pfxB   sfx  ) k      (  pfxB      sfx  ) k+nb
NB.      n+1-k  k                n+1-k-nb  k+nb
NB.
NB. Notes:
NB. - gerq2 emulates LAPACK's xGERQ2
NB. - gerq2 and gerqf are topologic equivalents
NB. - if triangular matrix diagonal's non-negativity is not
NB.   required, then larfp* may be replaced by faster larfg*

gerq2=: ((1&({::)) _1 stitch (2&({::))) @ ((3 : 0) ^: ((0 _1&(ms $))`(];(0&{.);(0&({."1)))))
  'pfxT pfxB sfx'=. y
  z=. larfpbc {: pfxT
  pfxB=. z , pfxB
  pfxT=. (1 (_1}) z) larfrnbr }: pfxT
  (pfxT (; &: (}:"1)) pfxB) , < ((pfxT (, &: ({:"1)) pfxB) ,. sfx)
)

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
NB. then (with appropriate comparison tolerance)
NB.   Q -: unglq (k {. LQf)
NB.   I -: clean (mp ct) Q
NB.   A -: L mp Q
NB.   (] -: clean @ ((         trl   @( 0 _1&}.)) mp  unglq)@gelqf) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGELQF

gelqf=: ((0&({::)) ,. (1&({::)) , (gelq2 @ (2&({::)))) @ ((3 : 0) ^: ((qfi@(0 _1&(ms $)))`((0&({."1));(0&{.);]))) @ (,. & 0)
  'pfx sfxT sfxB'=. y
  nb=. QFNB <. k=. <./ 0 _1 ms $ sfxB
  Z=. gelq2 nb {. sfxB
  sfxT=. sfxT , Z
  sfxB=. (tru1 Z) larfbrnfr nb }. sfxB
  (pfx ,. (sfxT (, &: (nb & ({."1))) sfxB)) ; (sfxT (; &: (nb&(}."1))) sfxB)
)

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
NB.   I -: clean (mp~ ct) Q
NB.   A -: Q mp L
NB.   (] -: clean @ ((((-~/@$) trl ])@( 1  0&}.)) mp~ ungql)@geqlf) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGEQLF

geqlf=: (((geql2 @ (0&({::))) ,. (1&({::))) , (2&({::))) @ ((3 : 0) ^: ((qfi@(_1 0&(ms $)))`(];(0&({."1));(0&{.)))) @ (0 & ,)
  'pfxL pfxR sfx'=. y
  nnb=. - QFNB <. k=. <./ _1 0 ms $ pfxL
  Z=. geql2 nnb {."1 pfxL
  pfxR=. Z ,. pfxR
  pfxL=. ((tru1~ (-~/@$)) Z) larfblcbc nnb }."1 pfxL
  (pfxL (; & (nnb&}.)) pfxR) , < ((pfxL (,. & (nnb & {.)) pfxR) , sfx)
)

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
NB.   I -: clean (mp~ ct) Q
NB.   A -: Q mp R
NB.   (] -: clean @ ((         tru   @(_1  0&}.)) mp~ ungqr)@geqrf) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGEQRF

geqrf=: ((0&({::)) , (1&({::)) ,. (geqr2 @ (2&({::)))) @ ((3 : 0) ^: ((qfi@(_1 0&(ms $)))`((0&{.);(0&({."1));]))) @ (, & 0)
  'pfx sfxL sfxR'=. y
  nb=. QFNB <. k=. <./ _1 0 ms $ sfxR
  Z=. geqr2 nb {."1 sfxR
  sfxL=. sfxL ,. Z
  sfxR=. (trl1 Z) larfblcfc nb }."1 sfxR
  (pfx , (sfxL (,. & (nb & {.)) sfxR)) ; (sfxL (; & (nb&}.)) sfxR)
)

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
NB.   I -: clean (mp ct) Q
NB.   A -: R mp Q
NB.   (] -: clean @ ((((-~/@$) tru ])@( 0  1&}.)) mp  ungrq)@gerqf) A
NB.
NB. Notes:
NB. - emulates LAPACK's xGERQF

NB. RQf=. gerqf A
gerqf=: (((gerq2 @ (0&({::))) , (1&({::))) ,. (2&({::))) @ ((3 : 0) ^: ((qfi@(0 _1&(ms $)))`(];(0&{.);(0&({."1))))) @ (0 & ,.)
  'pfxT pfxB sfx'=. y
  nnb=. - QFNB <. k=. <./ 0 _1 ms $ pfxT
  Z=. gerq2 nnb {. pfxT
  pfxB=. Z , pfxB
  pfxT=. ((trl1~ (-~/@$)) Z) larfbrnbr nnb }. pfxT
  (pfxT (; &: (nnb&(}."1))) pfxB) , < ((pfxT (, &: (nnb & ({."1))) pfxB) ,. sfx)
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. tgeqf
NB. Test orthogonal factorization algorithms:
NB. - built-in: 128!:0
NB. - LAPACK addon: gelqf geqlf geqrf gerqf
NB. - mt addon: gelqf geqlf geqrf gerqf
NB. by matrix given
NB.
NB. Syntax: tgeqf A
NB. where A - general m×n-matrix
NB.
NB. Formula:
NB. - berr for LQ: berr := ||A-L*Q||/(ε*n*||A||)
NB. - berr for QL: berr := ||A-Q*L||/(ε*m*||A||)
NB. - berr for QR: berr := ||A-Q*R||/(ε*m*||A||)
NB. - berr for RQ: berr := ||A-R*Q||/(ε*n*||A||)

tgeqf=: 3 : 0
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'gelqf geqlf geqrf gerqf'

  rcond=. ((_."_)`(norm1 con getri) @. (=/@$)) y  NB. meaninigful for square matrices only

  ('128!:0' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(- (mp & >/)))%(FP_EPS*(#*norm1)@[)))) y

  ('2b1110 & gelqf_jlapack_' tmonad (]`({. , (,.  &. > / @ }.))`(rcond"_)`(_."_)`((norm1@(- ((mp  unglq) & > /)))%((FP_EPS*c*norm1)@[)))) y
  ('2b0111 & geqlf_jlapack_' tmonad (]`({: , (,~  &. > / @ }:))`(rcond"_)`(_."_)`((norm1@(- ((mp~ ungql) & > /)))%((FP_EPS*#*norm1)@[)))) y
  ('2b0111 & geqrf_jlapack_' tmonad (]`({: , (,   &. > / @ }:))`(rcond"_)`(_."_)`((norm1@(- ((mp~ ungqr) & > /)))%((FP_EPS*#*norm1)@[)))) y
  ('2b1110 & gerqf_jlapack_' tmonad (]`({. , (,.~ &. > / @ }.))`(rcond"_)`(_."_)`((norm1@(- ((mp  ungrq) & > /)))%((FP_EPS*c*norm1)@[)))) y

  ('gelqf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(( trl         @:(}:"1)@]) (- dbg2 'lq-') unmlqrc~))%((FP_EPS*c*norm1)@[)))) y
  ('geqlf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(((trl~(-~/@$))@  }.   @]) (- dbg2 'ql-') unmqllc~))%((FP_EPS*#*norm1)@[)))) y
  ('geqrf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(( tru         @  }:   @]) (- dbg2 'qr-') unmqrlc~))%((FP_EPS*#*norm1)@[)))) y
  ('gerqf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(((tru~(-~/@$))@:(}."1)@]) (- dbg2 'rq-') unmrqrc~))%((FP_EPS*c*norm1)@[)))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testqf
NB.
NB. Description:
NB.   Test orthogonal factorization algorithms by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   mkge testqf (m,n)
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testqf 150 100

testqf=: 1 : 'EMPTY [ tgeqf @ u'
