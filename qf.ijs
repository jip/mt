NB. Orthogonal factorization
NB.
NB. gelqf     LQ factorization of a general matrix
NB. geqlf     QL factorization of a general matrix
NB. geqrf     QR factorization of a general matrix
NB. gerqf     RQ factorization of a general matrix
NB.
NB. tzlzf     LQ factorization of a trapezoidal matrix
NB. tzzlf     QL factorization of a trapezoidal matrix
NB. tzzrf     QR factorization of a trapezoidal matrix
NB. tzrzf     RQ factorization of a trapezoidal matrix
NB.
NB. testgeqf  Test gexxf by general matrix given
NB. testtzqf  Test tzxxf by trapezoidal matrix
NB. testqf    Adv. to make verb to test gexxf and tzxxf by
NB.           matrix of generator and shape given
NB.
NB. Version: 0.9.0 2012-11-15
NB.
NB. Copyright 2010-2012 Igor Zhuravlov
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

NB. ---------------------------------------------------------
NB. Blocked code constants

QFNB=: 32   NB. block size limit
QFNX=: 128  NB. crossover point, QFNX ≥ QFNB

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
NB. - models LAPACK's xGELQ2
NB. - gelq2 and gelqf are topologic equivalents
NB. - if triangular matrix diagonal's non-negativity is
NB.   required, then larfg* should be replaced by larfp*

gelq2=: 3 : 0
  pfx=. 0 {."1 y
  sfxT=. 0 {. y
  while. -. 0 e. 0 _1 + $ y do.
    z=. larfgfc {. y
    sfxT=. sfxT , z
    y=. ((1) 0} z) larfrnfr }. y
    pfx=. pfx ,. sfxT ,&:({."1) y
    sfxT=. }."1 sfxT
    y=. }."1 y
  end.
  pfx stitcht sfxT
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
NB. - models LAPACK's xGEQL2
NB. - geql2 and geqlf are topologic equivalents
NB. - if triangular matrix diagonal's non-negativity is
NB.   required, then larfg* should be replaced by larfp*

geql2=: 3 : 0
  pfxR=. 0 {."1 y
  sfx=. 0 {. y
  while. -. 0 e. _1 0 + $ y do.
    z=. larfgb {:"1 y
    pfxR=. z ,. pfxR
    y=. ((1) _1} z) larflcbc }:"1 y
    sfx=. sfx ,~ y ,&{: pfxR
    y=. }: y
    pfxR=. }: pfxR
  end.
  pfxR appendr sfx
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
NB. - models LAPACK's xGEQR2
NB. - gerq2 and geqrf are topologic equivalents
NB. - if triangular matrix diagonal's non-negativity is
NB.   required, then larfg* should be replaced by larfp*

geqr2=: 3 : 0
  pfx=. 0 {. y
  sfxL=. 0 {."1 y
  while. -. 0 e. _1 0 + $ y do.
    z=. larfgf {."1 y
    sfxL=. sfxL ,. z
    y=. ((1) 0} z) larflcfc }."1 y
    pfx=. pfx , sfxL ,&{. y
    sfxL=. }. sfxL
    y=. }. y
  end.
  pfx , sfxL
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
NB. - models LAPACK's xGERQ2
NB. - gerq2 and gerqf are topologic equivalents
NB. - if triangular matrix diagonal's non-negativity is
NB.   required, then larfg* should be replaced by larfp*

gerq2=: 3 : 0
  pfxB=. 0 {. y
  sfx=. 0 {."1 y
  while. -. 0 e. 0 _1 + $ y do.
    z=. larfgbc {: y
    pfxB=. z , pfxB
    y=. ((1) _1} z) larfrnbr }: y
    sfx=. sfx ,.~ y ,&:({:"1) pfxB
    y=. }:"1 y
    pfxB=. }:"1 pfxB
  end.
  pfxB stitchb sfx
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
NB.   LQf - m×(n+1)-matrix, combined L and Qf (unit
NB.         diagonal not stored)
NB.   L   - m×k-matrix, lower triangular
NB.   Qf  - k×(n+1)-matrix, unit upper triangular, the Q
NB.         represented in factored form
NB.   Q   - k×n-matrix with orthonormal rows, which is
NB.         defined as the first k rows of a product of m
NB.         elementary reflectors of order n:
NB.           Q = Π{H(i)',i=m-1:0}
NB.         where
NB.           H(m-1:k)≡H(v(m-1:k),τ(m-1:k))=H(0,0)=I
NB.   k   = min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   I -: po Q
NB.   A -: L mp Q
NB.   (] -: ((         trl   @:(}:"1)) mp  unglq)@gelqf) A
NB. where
NB.   LQf=. gelqf A
NB.   L=. (trl@:(}:"1)) LQf
NB.   Q=. unglq LQf
NB.
NB. Notes:
NB. - models LAPACK's xGELQF

gelqf=: 3 : 0
  y=. y ,. 0
  pfx=. 0 {."1 y
  sfxT=. 0 {. y
  while. -. 0 e. QFNX < 0 _1 + $ y do.
    nb=. <./ QFNB , 0 _1 + $ y
    Z=. gelq2 nb {. y
    sfxT=. sxT , Z
    y=. (tru1 Z) larfbrnfr nb }. y
    pfx=. pfx ,. sfxT ,&(nb&({."1)) y
    sfxT=. nb }."1 sfxT
    y=. nb }."1 y
  end.
  pfx ,. sfxT , gelq2 y
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
NB.   QfL - (m+1)×n-matrix, combined Qf (unit diagonal not
NB.         stored) and L
NB.   Qf  - (m+1)×k-matrix, unit upper triangular, the Q
NB.         represented in factored form
NB.   L   - k×n-matrix, lower triangular
NB.   Q   - m×k-matrix with orthonormal columns, which is
NB.         defined as the last k columns of a product of n
NB.         elementary reflectors of order m:
NB.           Q = Π{H(i),i=n-1:0}
NB.         where
NB.           H(n-1:k)≡H(v(n-1:k),τ(n-1:k))=H(0,0)=I
NB.   k   = min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   I -: (mp~ ct) Q
NB.   A -: Q mp L
NB.   (] -: (((trl~ -~/@$)@  }.   ) mp~ ungql)@geqlf) A
NB. where
NB.   QfL=. geqlf A
NB.   Q=. ungql QfL
NB.   L=. ((trl~ -~/@$)@}.) QfL
NB.
NB. Notes:
NB. - models LAPACK's xGEQLF

geqlf=: 3 : 0
  y=. 0 , y
  pfxR=. 0 {."1 y
  sfx=. 0 {. y
  while. -. 0 e. QFNX < _1 0 + $ y do.
    nb=. - <./ QFNB , _1 0 + $ y
    Z=. geql2 nb {."1 y
    pfxR=. Z ,. pfxR
    y=. ((tru1~ -~/@$) Z) larfblcbc nb }."1 y
    sfx=. sfx ,~ y ,.&(nb&{.) pfxR
    y=. nb }. y
    pfxR=. nb }. pfxR
  end.
  sfx ,~ pfxR ,.~ geql2 y
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
NB.   QfR - (m+1)×n-matrix, combined Qf (unit diagonal not
NB.         stored) and R
NB.   Qf  - (m+1)×k-matrix, unit lower triangular, the Q
NB.         represented in factored form
NB.   R   - k×n-matrix, upper triangular
NB.   Q   - m×k-matrix with orthonormal columns, which is
NB.         defined as the first k columns of a product of n
NB.         elementary reflectors of order m:
NB.           Q = Π{H(i),i=0:n-1}
NB.         where
NB.           H(k:n-1)≡H(v(k:n-1),τ(k:n-1))=H(0,0)=I
NB.   k   = min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   I -: (mp~ ct) Q
NB.   A -: Q mp R
NB.   (] -: ((         tru   @  }:   ) mp~ ungqr)@geqrf) A
NB. where
NB.   QfR=. geqrf A
NB.   Q=. ungqr QfR
NB.   R=. tru }: QfR
NB.
NB. Notes:
NB. - models LAPACK's xGEQRF

geqrf=: 3 : 0
  y=. y , 0
  pfx=. 0 {. y
  sfxL=. 0 {."1 y
  while. -. 0 e. QFNX < _1 0 + $ y do.
    nb=. <./ QFNB , _1 0 + $ y
    Z=. geqr2 nb {."1 y
    sfxL=. sfxL ,. Z
    y=. (trl1 Z) larfblcfc nb }."1 y
    pfx=. pfx , sfxL ,.&(nb&{.) y
    sfxL=. nb }. sfxL
    y=. nb }. y
  end.
  pfx , sfxL ,. geqr2 y
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
NB.   RQf - m×(n+1)-matrix, combined R and Qf (unit diagonal
NB.         not stored)
NB.   R   - m×k-matrix, upper triangular
NB.   Qf  - k×(n+1)-matrix, unit lower triangular, the Q
NB.         represented in factored form
NB.   Q   - k×n-matrix with orthonormal rows, which is
NB.         defined as the last k rows of a product of m
NB.         elementary reflectors of order n:
NB.           Q = Π{H(i)',i=0:m-1}
NB.         where
NB.           H(k:m-1)≡H(v(k:m-1),τ(k:m-1))=H(0,0)=I
NB.   k   = min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   I -: po Q
NB.   A -: R mp Q
NB.   (] -: (((tru~ -~/@$)@:(}."1)) mp  ungrq)@gerqf) A
NB. where
NB.   RQf=. gerqf A
NB.   R=. ((tru~ -~/@$)@:(}."1)) RQf
NB.   Q=. ungrq RQf
NB.
NB. Notes:
NB. - models LAPACK's xGERQF

gerqf=: 3 : 0
  y=. 0 ,. y
  pfxB=. 0 {. y
  sfx=. 0 {."1 y
  while. -. 0 e. QFNX < 0 _1 + $ y do.
    nb=. - <./ QFNB , 0 _1 + $ y
    Z=. gerq2 nb {. y
    pfxB=. Z , pfxB
    y=. ((trl1~ -~/@$) Z) larfbrnbr nb }. y
    sfx=. sfx ,.~ y ,&(nb&({."1)) pfxB
    y=. nb }."1 y
    pfxB=. nb }."1 pfxB
  end.
  sfx ,.~ pfxB ,~ gerq2 y
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgeqf
NB.
NB. Description:
NB.   Test orthogonal factorization algorithms:
NB.   - 128!:0 (built-in)
NB.   - gelqf geqlf geqrf gerqf (math/lapack addon)
NB.   - gelqf geqlf geqrf gerqf (math/mt addon)
NB.   by general matrix given
NB.
NB. Syntax:
NB.   testgeqf A
NB. where
NB.   A - m×n-matrix
NB.
NB. Formula:
NB. - LQ berr := max( ||L - A * Q^H|| / (FP_EPS * ||A|| * n), ||Q * Q^H - I|| / (FP_EPS * n) )
NB. - QL berr := max( ||L - Q^H * A|| / (FP_EPS * ||A|| * m), ||Q^H * Q - I|| / (FP_EPS * m) )
NB. - QR berr := max( ||R - Q^H * A|| / (FP_EPS * ||A|| * m), ||Q^H * Q - I|| / (FP_EPS * m) )
NB. - RQ berr := max( ||R - A * Q^H|| / (FP_EPS * ||A|| * n), ||Q * Q^H - I|| / (FP_EPS * n) )

testgeqf=: 3 : 0
  require :: ] '~addons/math/lapack/lapack.ijs'
  need_jlapack_ :: ] 'gelqf geqlf geqrf gerqf'

  rcond=. (_."_)`gecon1@.(=/@$) y  NB. meaninigful for square matrices only

  ('128!:0'                tmonad (]`]                 `(rcond"_)`(_."_)`((norm1@((1 {:: ])              - (( mp~ ct                   )  0&{::)) % (FP_EPS * norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@(<: upddiag)@(mp~ ct       )@(0&{::))))) y

  ('2b1110&gelqf_jlapack_' tmonad (]`({. ,  ,. &.>/@}.)`(rcond"_)`(_."_)`((norm1@((0 {:: ])              - (((   <./ @$@]) {."1 unmlqrc)~ 1&{::)) % (FP_EPS * norm1 * c)@[) >. ((% FP_EPS * c)~ norm1@(<: upddiag)@(unmlqrc unglq)@(1&{::))))) y
  ('2b0111&geqlf_jlapack_' tmonad (]`({: ,~ , ~&.>/@}:)`(rcond"_)`(_."_)`((norm1@((1 {:: ])              - (((-@(<./)@$@]) {.   unmqllc)~ 0&{::)) % (FP_EPS * norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@(<: upddiag)@(unmqllc ungql)@(0&{::))))) y
  ('2b0111&geqrf_jlapack_' tmonad (]`({: ,~ ,  &.>/@}:)`(rcond"_)`(_."_)`((norm1@((1 {:: ])              - (((   <./ @$@]) {.   unmqrlc)~ 0&{::)) % (FP_EPS * norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@(<: upddiag)@(unmqrlc ungqr)@(0&{::))))) y
  ('2b1110&gerqf_jlapack_' tmonad (]`({. ,  ,.~&.>/@}.)`(rcond"_)`(_."_)`((norm1@((0 {:: ])              - (((-@(<./)@$@]) {."1 unmrqrc)~ 1&{::)) % (FP_EPS * norm1 * c)@[) >. ((% FP_EPS * c)~ norm1@(<: upddiag)@(unmrqrc ungrq)@(1&{::))))) y

  ('gelqf'                 tmonad (]`]                 `(rcond"_)`(_."_)`((norm1@( trl        @:(}:"1)@] -  ((   <./ @$@]) {."1 unmlqrc)~       ) % (FP_EPS * norm1 * c)@[) >. ((% FP_EPS * c)~ norm1@(<: upddiag)@(unmlqrc unglq)        )))) y
  ('geqlf'                 tmonad (]`]                 `(rcond"_)`(_."_)`((norm1@((trl~ -~/@$)@  }.   @] -  ((-@(<./)@$@]) {.   unmqllc)~       ) % (FP_EPS * norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@(<: upddiag)@(unmqllc ungql)        )))) y
  ('geqrf'                 tmonad (]`]                 `(rcond"_)`(_."_)`((norm1@( tru        @  }:   @] -  ((   <./ @$@]) {.   unmqrlc)~       ) % (FP_EPS * norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@(<: upddiag)@(unmqrlc ungqr)        )))) y
  ('gerqf'                 tmonad (]`]                 `(rcond"_)`(_."_)`((norm1@((tru~ -~/@$)@:(}."1)@] -  ((-@(<./)@$@]) {."1 unmrqrc)~       ) % (FP_EPS * norm1 * c)@[) >. ((% FP_EPS * c)~ norm1@(<: upddiag)@(unmrqrc ungrq)        )))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testqf
NB.
NB. Description:
NB.   Adv. to make verb to test gexxx by matrix of generator
NB.   and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testqf
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
NB.     ?@$&0 testqf_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testqf_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testqf_mt_ 150 200

testqf=: 1 : 'EMPTY_mt_ [ testgeqf_mt_@u'
