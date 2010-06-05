NB. ref.ijs
NB. Reflections
NB.
NB. larfg      Dyad to generate an elementary reflector
NB. larfp      Dyad to generate an elementary reflector with
NB.            β≥0
NB. larfxxx    Monads to generate an elementary reflector
NB. larftxx    Monads to form the triangular factor of a
NB.            block reflector
NB. larfxxxx   Dyads to apply a reflector or its transpose to
NB.            a matrix from either the left or the right
NB. larfbxxxx  Dyads to build a block reflector by larftxx
NB.            and to apply it or its transpose to a matrix
NB.            from either the left or the right
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. larft
NB.
NB. Description:
NB.   Template adv. to make monads to form the triangular
NB.   factor T of a block reflector H, which is defined as a
NB.   product of k elementary reflectors.
NB.
NB. Syntax:
NB.   vapp=. vV`vP`vTi`va0`vs1 larft iost
NB. where
NB.   vV   - verb to prepare V0r; input must be a triangular
NB.          (trapezoidal) matrix with zeroed tau elements;
NB.          vectors within output must get horizontal
NB.          orientation; is called as:
NB.            V0r=. vV V0
NB.   vP   - verb to calculate coefficients matrix P; is
NB.          called as:
NB.            P=. tau vP V0r
NB.   vTi  - verb to stitch new column while building matrix
NB.          T on current iteration; is called as:
NB.            Tiupd=. P vTi Ti
NB.   va0  - verb to append or prepend Tiupd with row of
NB.          values from scalar x; is called as:
NB.            Tiupd0=. (0 & va0) Tiupd
NB.   vs1  - verb to place scalar x in the diagonal element
NB.          of column just added by va0; is called as:
NB.            Tiupd1=. (1 & vs1) Tiupd0
NB.   iost - IOS tau in Vtau
NB.   vapp - dyad to form the triangular factor T; is called
NB.          as:
NB.            T=. vapp Vtau
NB.   Vtau - matrix V with appended or stitched vector tau
NB.   V0   - unit triangular (trapezoidal) matrix with zeroed
NB.          tau elements
NB.   V0r  - k×(n+1)-matrix, unit triangular (trapezoidal)
NB.          with rows - elementary reflectors v[i], i=0:k-1,
NB.          and zeroed tau elements
NB.   P    - coefficients, it is either the k×k-matrix
NB.          (-tau)*Pfwd for forward direction, or
NB.          k×(k-1)-matrix (0 1 }. ((-tau)*Pbwd)) for
NB.          backward direction (item-by-row product in both
NB.          cases)
NB.   Pfwd - k×k-matrix ('*' means trash):
NB.            [ *        0        0        ... 0          0   ]
NB.            [ p[1,0]   *        0        ... 0          0   ]
NB.            [ p[2,0]   p[2,1]   *        ... 0          0   ]
NB.            [ ...      ...      ...      ... ...        ... ]
NB.            [ p[k-1,0] p[k-1,1] p[k-1,2] ... p[k-1,k-2] *   ]
NB.          where for i=0:k-1
NB.            p[i] = (p[i,0:i-1] , *) = V0r[0:i,0:n-1] * V0r'[i,0:n-1]
NB.   Pbwd - k×k-matrix ('*' means trash):
NB.            [ *   0        0        0        ... 0          ]
NB.            [ *   p[1,0]   0        0        ... 0          ]
NB.            [ *   p[2,0]   p[2,1]   0        ... 0          ]
NB.            [ ... ...      ...      ...      ... ...        ]
NB.            [ *   p[k-1,0] p[k-1,1] p[k-1,2] ... p[k-1,k-2] ]
NB.          where for i=0:k-1, j=k-1-i
NB.            p[i] = (* , p[i,0:i-1]) = V0r[j:k-1,0:n-1] * V0r'[j,0:n-1]
NB.   tau  - k-vector of τ[0:k-1], corresp. to V,
NB.            tau -: iost { Vtau
NB.   T    - k×k-matrix, upper or lower triangular
NB.
NB. Algorithm:
NB.   1) prepare V0: V0=. iot 0} Vtau
NB.   2) prepare V0r: V0r=. vV V0
NB.   3) extract tau: tau=. iot { Vtau
NB.   4) prepare coefficients matrix P: P=. tau vP V0r
NB.   5) let T[i]=0×0-matrix, start loop i=0:k-1:
NB.      5.1) extract vector p[i]: pi=P[i,0:i-1]
NB.      5.2) calc new column: c[i]: ci=. pi mp Ti
NB.      5.3) stitch c[i] to T[i]: Ti=. Ti vTi ci
NB.      5.4) append zero row to T[i]: Ti=. (0 & va0) Ti
NB.      5.5) write 1 in the diagonal element of appended
NB.           zero row: Ti=. (1 & vs1) Ti
NB.   6) multiply T[i] by tau (item-by-row): T=. tau * Ti
NB.
NB. Notes:
NB. - emulates LAPACK's xLARFT, but without zeros scanning in
NB.   v[i]

larft=: 2 : '(n&{) ([ * (EMPTY (((1 & (m gi 4)) @ (0 & (m gi 3)) @ (] (m gi 2) (] mp (([ {. { )~ #)))) ^: (#@[))~ (m gi 1))) ((m gi 0) @ (0 & (n })))'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. larfg
NB. larfp
NB.
NB. Description:
NB.   Generate an elementary reflector H of order n such that
NB.   H'*y = β*e1, where H=I-v*τ*v', H'*H=I. H is represented
NB.   in factored form by n-vector v and scalar τ. Vector v
NB.   can have either forward (α in head) or backward (α in
NB.   tail) direction. Input and output vector directions are
NB.   the same.
NB.
NB. Syntax:
NB.   z=. ios larfg y
NB.   z=. ios larfp y
NB. where
NB.   ios - 2-vector of integers (ioa,iot)
NB.   ioa - lIO α in y
NB.   iot - lIO pre-allocated scalar for τ in y
NB.   y   - (n+1)-vector or (n+1)×1-matrix or 1×(n+1)-matrix
NB.         having scalar α at index ioa, any scalar at index
NB.         iot, and vector x[1:n-1] in the rest elements,
NB.         vector to reflect is (α,x[1:n-1]), α ∊ ℂ,
NB.         x ∊ ℂⁿ⁻¹
NB.   z   - (n+1)-vector or (n+1)×1-matrix or 1×(n+1)-matrix
NB.         (y and z shapes are match) having scalar β at
NB.         index ioa, scalar τ at index iot, and vector
NB.         v[1:n-1] in the rest elements, reflection result
NB.         is vector (1,v[1:n-1]), 1 is not stored, β ∊ ℝ,
NB.         (larfp provides β≥0), v ∊ ℂⁿ⁻¹, τ ∊ ℂ
NB.
NB. Applications:
NB. - reflect vector (α,x) by larfg and store τ at tail:
NB.     z=. 0 _1 larfg (α,x,0)
NB.     v=. 1 (0: }) }: z
NB.     'beta tau'=. 0 _1 ({,) z
NB. - reflect vector (x,α) by larfp and store τ at head:
NB.     z=. _1 0 larfp (0,x,α)
NB.     v=. 1 (_1: }) }. z
NB.     'beta tau'=. _1 0 ({,) z
NB.
NB. References:
NB. [1] James W. Demmel, Mark Hoemmen, Yozo Hida, and E.
NB.     Jason Riedy. (2008) Non-Negative Diagonals and High
NB.     Performance on Low-Profile Matrices from Householder
NB.     QR. UCB/EECS-2008-76, May 30, 2008
NB.     LAPACK Working Note 203.
NB.     http://www.netlib.org/lapack/lawns/downloads/
NB.
NB. Notes:
NB. - IEEE floating point configuration is encoded implicitly
NB. - larfg emulates LAPACK's xLARFG
NB. - larfp emulates LAPACK's xLARFP
NB. - larfp provides τ=2 ↔ ||x||_2 ≈ 0
NB. - not zeroing v in larfp in case τ=2 relies on fact:
NB.   'comparison tolerance tol>0'; otherwise (tol=0) x
NB.   should be filled by zeros

larfg=: 4 : 0
  'ioalpha iotau'=. x
  alpha=. ioalpha ({,) y
  xnorm=. norms (0 (x"_) } y)                               NB. ||x||_2
  if. xnorm (+. & (0 & ~:)) (11 o. alpha) do.
    beta=. - (9 o. alpha) condneg norms alpha , xnorm       NB. β=-copysign(||y||_2,Re(α))
    y=. beta (iotau " _) } y                                NB. write in-place β
    if. FP_SFMIN > | beta do.
      y=. y % FP_SFMIN                                      NB. scale (α,x[1],...,x[n-1],β)
    end.
    dzeta=. -/ x ({,) y                                     NB. ζ=α_scaled-β_scaled
    tau=. (- dzeta) % (iotau ({,) y)                        NB. τ=-ζ/β_scaled
    y=. y % dzeta                                           NB. y=y/ζ
    y=. (beta , tau) (x " _) } y                            NB. replace α_scaled by β_unscaled and β_scaled by τ
  else.
    y=. 0 (iotau " _) } y
  end.
  y
)

larfp=: 4 : 0
  'ioalpha iotau'=. x
  alpha=. ioalpha ({,) y
  xnorm=. norms (0 (x"_) } y)                               NB. ||x||_2
  if. (0 = xnorm) do.
    y=. ((| , (1 - *)) alpha) (x " _) } y                   NB. replace in-place α by |β| and τ by (1-α/|α|)
  else.
    beta=. - (9 o. alpha) condneg norms alpha , xnorm       NB. β=-copysign(||y||_2,Re(α))
    y=. (| beta) (iotau " _) } y                            NB. write in-place |β|
    if. FP_SFMIN > | beta do.
      y=. y % FP_SFMIN                                      NB. scale (α,x[1],...,x[n-1],|β|)
      xnorm=. xnorm % FP_SFMIN
    end.
    if. 0 <: beta do.
      dzeta=. -/ x ({,) y                                   NB. ζ=α_scaled-|β_scaled|
      tau=. (- dzeta) % (iotau ({,) y)                      NB. τ=-ζ/|β_scaled|
    else.
      beta=. - beta                                         NB. |β_unscaled|
      'realpha imalpha'=. +. ioalpha ({,) y                 NB. Re(α_scaled) , Im(α_scaled)
      gamma=. realpha + iotau ({,) y                        NB. γ=Re(α_scaled)+|β_scaled|
      delta=. (imalpha , xnorm) (- @ (+/) @ ([ * %)) gamma  NB. δ=-(Im(α_scaled)*(Im(α_scaled)/γ)+||x||_2*(||x||_2/γ))
      dzeta=. delta j. imalpha                              NB. ζ=δ+i*Im(α_scaled)
      tau=. - dzeta % iotau ({,) y                          NB. τ=-ζ/|β_scaled|
    end.
    y=. y % dzeta
    y=. (beta , tau) (x " _) } y                            NB. replace α_scaled by |β_unscaled| and |β_scaled| by τ
  end.
)

NB. ---------------------------------------------------------
NB. Verb       Input             Output                   β
NB. larfgf     (α x[1:n-1] 0)    (β v[1:n-1] τ)          ∊ℝ
NB. larfgfc    (α x[1:n-1] 0)    (β v[1:n-1] conj(τ))    ∊ℝ
NB. larfgb     (0 x[1:n-1] α)    (τ v[1:n-1] β)          ∊ℝ
NB. larfgbc    (0 x[1:n-1] α)    (conj(τ) v[1:n-1] β)    ∊ℝ
NB. larfpf     (α x[1:n-1] 0)    (β v[1:n-1] τ)          ≥0
NB. larfpfc    (α x[1:n-1] 0)    (β v[1:n-1] conj(τ))    ≥0
NB. larfpb     (0 x[1:n-1] α)    (τ v[1:n-1] β)          ≥0
NB. larfpbc    (0 x[1:n-1] α)    (conj(τ) v[1:n-1] β)    ≥0
NB.
NB. Description:
NB.   Monads to generate an elementary reflector, see larfg,
NB.   larfp for details.

larfgf=: 0 _1 & larfg
larfgfc=: (+ updl _1) @ larfgf
larfgb=: _1 0 & larfg
larfgbc=: (+ updl 0) @ larfgb

larfpf=: 0 _1 & larfp
larfpfc=: (+ updl _1) @ larfpf
larfpb=: _1 0 & larfp
larfpbc=: (+ updl 0) @ larfpb

NB. ---------------------------------------------------------
NB. Verb       Direction    Layout
NB. larftbc    backward     columnwise
NB. larftbr    backward     rowwise
NB. larftfc    forward      columnwise
NB. larftfr    forward      rowwise
NB.
NB. Description:
NB.   Monads to form the triangular factor T of a block
NB.   reflector H. If direction is forward, then:
NB.     H = H(0) * H(1) * ... * H(k-1) and T is upper
NB.   triangular, otherwise direction is backward, and:
NB.     H = H(k-1) * ... * H(1) * H(0) and T is lower
NB.   triangular. If layout is columnwise, then:
NB.     H = I - V * T * V' ,
NB.   otherwise layout is rowwise, and:
NB.     H = I - V' * T * V .
NB.
NB. Syntax:
NB.   Tl=. larftbc Vbct
NB.   Tl=. larftbr Vbrt
NB.   Tu=. larftfc Vfct
NB.   Tu=. larftfr Vfrt
NB. where
NB.   Vbct - (m+1)×k-matrix (tau , Vbc)
NB.   Vbc  - m×k-matrix, unit upper triangular (trapezoidal)
NB.          with 1s on (k-m)-th diagonal and 0s below
NB.   Vbrt - k×(n+1)-matrix (tau ,. Vbr)
NB.   Vbr  - k×n-matrix, unit lower triangular (trapezoidal)
NB.          with 1s on (n-k)-th diagonal and 0s above
NB.   Vfct - (m+1)×k-matrix (Vfc , tau)
NB.   Vfc  - m×k-matrix, unit lower triangular (trapezoidal)
NB.          with 1s on 0-th diagonal and 0s above
NB.   Vfrt - k×(n+1)-matrix (Vfr ,. tau)
NB.   Vfr  - k×n-matrix, unit upper triangular (trapezoidal)
NB.          with 1s on 0-th diagonal and 0s below
NB.   tau  - k-vector τ[0:k-1] corresp. to Vbc, Vbr, Vfc or
NB.          Vfr
NB.   Tl   - k×k-matrix, lower triangular
NB.   Tu   - k×k-matrix, upper triangular

larftbc=: ct`(|.@(0 1&}.)@((* -)~ ((mp ct) {.)\.))`(,.~)` ,  `( 0:}) larft IOSFR
larftbr=: ] `(|.@(0 1&}.)@((* -)~ ((mp ct) {.)\.))`(,.~)` ,  `( 0:}) larft IOSFC
larftfc=: ct`(             (* -)~ ((mp ct) {:)\  )` ,.  `(,~)`(_1:}) larft IOSLR
larftfr=: ] `(             (* -)~ ((mp ct) {:)\  )` ,.  `(,~)`(_1:}) larft IOSLC

NB. ---------------------------------------------------------
NB. Verb       Action   Side   Tran  Dir  Layout    eC    Used in
NB. larflcbc   H' * C   left   ct    bwd  col-wise  0, C  geql2
NB. larflcbr   H' * C   left   ct    bwd  rowwise   0, C
NB. larflcfc   H' * C   left   ct    fwd  col-wise  C, 0  geqr2,gehd2u
NB. larflcfr   H' * C   left   ct    fwd  rowwise   C, 0
NB. larflnbc   H  * C   left   none  bwd  col-wise  0, C  ung2l
NB. larflnbr   H  * C   left   none  bwd  rowwise   0, C
NB. larflnfc   H  * C   left   none  fwd  col-wise  C, 0  ung2r
NB. larflnfr   H  * C   left   none  fwd  rowwise   C, 0  gehd2l
NB. larfrcbc   C  * H'  right  ct    bwd  col-wise  0,.C
NB. larfrcbr   C  * H'  right  ct    bwd  rowwise   0,.C  ungr2
NB. larfrcfc   C  * H'  right  ct    fwd  col-wise  C,.0
NB. larfrcfr   C  * H'  right  ct    fwd  rowwise   C,.0  ungl2
NB. larfrnbc   C  * H   right  none  bwd  col-wise  0,.C
NB. larfrnbr   C  * H   right  none  bwd  rowwise   0,.C  gerq2
NB. larfrnfc   C  * H   right  none  fwd  col-wise  C,.0  gehd2u,gehd2l
NB. larfrnfr   C  * H   right  none  fwd  rowwise   C,.0  gelq2,unml2
NB.
NB. Description:
NB.   Dyads to apply a reflector or its transpose to a
NB.   matrix, from either the left or the right
NB.
NB. Syntax:
NB.   eCupd=. vtau larfxxxx eC
NB. where
NB.   eC    - matrix C to update with appended or stitched
NB.           trash vector
NB.   vtau  - v with appended tau
NB.   eCupd - being updated matrix C with appended or
NB.           stitched modified trash vector
NB.   v     - vector with 1 at head (forward direction) or
NB.           tail (backward direction)
NB.   tau   - scalar τ
NB.
NB. Notes:
NB. - emulates LAPACK's xLARF
NB. - larfxxxx and larfbxxxx are topological equivalents

larflcbc=: ] - [ */ (mp~ (+ @ ((0 & ( 0 })) * {.)))~    NB. C - v * ((v * τ)' * C)
larflcbr=: ] - (+ @ (* {.) @ [) */ ((0 ( 0) } [) mp ])  NB. C - (τ * v)' * (v * C)
larflcfc=: ] - [ */ (mp~ (+ @ ((0 & (_1 })) * {:)))~    NB. C - v * ((v * τ)' * C)
larflcfr=: ] - (+ @ (* {:) @ [) */ ((0 (_1) } [) mp ])  NB. C - (τ * v)' * (v * C)

larflnbc=: ] - [ */ (mp~ ({. * (+ @ (0 & ( 0 })))))~    NB. C - v * ((τ * v') * C)
larflnbr=: ] - ((+ * {.) @ [) */ (mp~ (0 & ( 0 })))~    NB. C - (v' * τ) * (v * C)
larflnfc=: ] - [ */ (mp~ ({: * (+ @ (0 & (_1 })))))~    NB. C - v * ((τ * v') * C)
larflnfr=: ] - ((+ * {:) @ [) */ (mp~ (0 & (_1 })))~    NB. C - (v' * τ) * (v * C)

larfrcbc=: ] - (mp (0 & ( 0 })))~ */ (+ @ (* {.) @ [)   NB. C - (C * v) * (v * τ)'
larfrcbr=: ] - (mp (+ @ ({. * (0 & ( 0 })))))~ */ [     NB. C - (C * (τ * v)') * v
larfrcfc=: ] - (mp (0 & (_1 })))~ */ (+ @ (* {:) @ [)   NB. C - (C * v) * (v * τ)'
larfrcfr=: ] - (mp (+ @ ({: * (0 & (_1 })))))~ */ [     NB. C - (C * (τ * v)') * v

larfrnbc=: ] - (mp ((0 & ( 0 })) * {.))~ */ (+ @ [)     NB. C - (C * (v * τ)) * v'
larfrnbr=: ] - (mp ((+ @ (0 & ( 0 }))) * {.))~ */ [     NB. C - (C * (v' * τ)) * v
larfrnfc=: ] - (mp ((0 & (_1 })) * {:))~ */ (+ @ [)     NB. C - (C * (v * τ)) * v'
larfrnfr=: ] - (mp ((+ @ (0 & (_1 }))) * {:))~ */ [     NB. C - (C * (v' * τ)) * v

NB. ---------------------------------------------------------
NB. Verb       Action   Side   Tran  Dir  Layout    eC    Used in
NB. larfblcbc  H' * C   left   ct    bwd  col-wise  0, C  geqlf
NB. larfblcbr  H' * C   left   ct    bwd  rowwise   0, C
NB. larfblcfc  H' * C   left   ct    fwd  col-wise  C, 0  geqrf,gehrdu
NB. larfblcfr  H' * C   left   ct    fwd  rowwise   C, 0  unmlq
NB. larfblnbc  H  * C   left   none  bwd  col-wise  0, C  ungql
NB. larfblnbr  H  * C   left   none  bwd  rowwise   0, C
NB. larfblnfc  H  * C   left   none  fwd  col-wise  C, 0  ungqr
NB. larfblnfr  H  * C   left   none  fwd  rowwise   C, 0
NB. larfbrcbc  C  * H'  right  ct    bwd  col-wise  0,.C
NB. larfbrcbr  C  * H'  right  ct    bwd  rowwise   0,.C  ungrq
NB. larfbrcfc  C  * H'  right  ct    fwd  col-wise  C,.0
NB. larfbrcfr  C  * H'  right  ct    fwd  rowwise   C,.0  unglq
NB. larfbrnbc  C  * H   right  none  bwd  col-wise  0,.C
NB. larfbrnbr  C  * H   right  none  bwd  rowwise   0,.C  gerqf
NB. larfbrnfc  C  * H   right  none  fwd  col-wise  C,.0
NB. larfbrnfr  C  * H   right  none  fwd  rowwise   C,.0  gelqf,unmlq
NB.
NB. Description:
NB.   Dyads to build and apply a block reflector or its
NB.   transpose to a matrix, from either the left or the
NB.   right
NB.
NB. Syntax:
NB.   eCupd=. Vtau larfbxxxx eC
NB. where
NB.   eC    - matrix C to update with appended or stitched
NB.           trash vector
NB.   Vtau  - matrix V with appended or stitched vector tau
NB.   eCupd - being updated matrix C with appended or
NB.           stitched modified trash vector
NB.   V     - unit triangular (trapezoidal) matrix
NB.   tau   - k-vector τ[0:k-1] corresp. to V
NB.
NB. Notes:
NB. - emulates LAPACK's sequence of calls to xLARFT and then
NB.   to xLARFB
NB. - larfxxxx and larfbxxxx are topological equivalents

larfblcbc=: ] - [ mp (mp~ (ct @ ((0 & (IOSFR })) mp larftbc)))~   NB. C - V * ((V * T)' * C)
larfblcbr=: ] - (ct @ (mp~ larftbr) @ [) mp ((0 IOSFC } [) mp ])  NB. C - (T * V)' * (V * C)
larfblcfc=: ] - [ mp (mp~ (ct @ ((0 & (IOSLR })) mp larftfc)))~   NB. C - V * ((V * T)' * C)
larfblcfr=: ] - (ct @ (mp~ larftfr) @ [) mp ((0 IOSLC } [) mp ])  NB. C - (T * V)' * (V * C)

larfblnbc=: ] - [ mp (mp~ (larftbc mp (ct @ (0 & (IOSFR })))))~   NB. C - V * ((T * V') * C)
larfblnbr=: ] - ((ct mp larftbr) @ [) mp (mp~ (0 & (IOSFC })))~   NB. C - (V' * T) * (V * C)
larfblnfc=: ] - [ mp (mp~ (larftfc mp (ct @ (0 & (IOSLR })))))~   NB. C - V * ((T * V') * C)
larfblnfr=: ] - ((ct mp larftfr) @ [) mp (mp~ (0 & (IOSLC })))~   NB. C - (V' * T) * (V * C)

larfbrcbc=: ] - (mp (0 & (IOSFR })))~ mp (ct @ (mp larftbc) @ [)  NB. C - (C * V) * (V * T)'
larfbrcbr=: ] - (mp (ct @ (larftbr mp (0 & (IOSFC })))))~ mp [    NB. C - (C * (T * V)') * V
larfbrcfc=: ] - (mp (0 & (IOSLR })))~ mp (ct @ (mp larftfc) @ [)  NB. C - (C * V) * (V * T)'
larfbrcfr=: ] - (mp (ct @ (larftfr mp (0 & (IOSLC })))))~ mp [    NB. C - (C * (T * V)') * V

larfbrnbc=: ] - (mp ((0 & (IOSFR })) mp larftbc))~ mp (ct @ [)    NB. C - (C * (V * T)) * V'
larfbrnbr=: ] - (mp ((ct @ (0 & (IOSFC }))) mp larftbr))~ mp [    NB. C - (C * (V' * T)) * V
larfbrnfc=: ] - (mp ((0 & (IOSLR })) mp larftfc))~ mp (ct @ [)    NB. C - (C * (V * T)) * V'
larfbrnfr=: ] - (mp ((ct @ (0 & (IOSLC }))) mp larftfr))~ mp [    NB. C - (C * (V' * T)) * V

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB.*tlarfg v test larfg

tlarfg=: 3 : 0
  EMPTY
)

NB. ---------------------------------------------------------
NB.*tlarfp v test larfp

tlarfp=: 3 : 0
  EMPTY
)

NB. ---------------------------------------------------------
NB.*tref v test reflectors with data suppplied

tref=: 3 : 0
  EMPTY
)

NB. ---------------------------------------------------------
NB.*testref v test reflectors

testref=: 3 : 0
  EMPTY
)
