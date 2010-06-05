NB. ref.ijs
NB. Reflections
NB.
NB. larfg      Dyad to generate an elementary reflector
NB. larfp      Dyad to generate an elementary reflector with
NB.            β≥0
NB. larfxxx    Monads to generate an elementary reflector
NB. larftxx    Monads to form the triangular factor of a
NB.            block reflector
NB. larfxxxx   Dyads to apply an elementary reflector or its
NB.            transpose to a matrix from either the left or
NB.            the right
NB. larfbxxxx  Dyads to build a block reflector by larftxx
NB.            and to apply it or its transpose to a matrix
NB.            from either the left or the right
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Concepts

NB. ---------------------------------------------------------
NB. larftxx
NB.
NB. Description:
NB.   Monads to form the triangular factor T of a block
NB.   reflector H, which is defined as a product of k
NB.   elementary reflectors.
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
NB. Notes:
NB. - emulates LAPACK's xLARFT, but without zeros scanning in
NB.   v[i]

NB. =========================================================
NB. Local definitions

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
larfgfc=: _1 + upd1 larfgf
larfgb=: _1 0 & larfg
larfgbc=: 0 + upd1 larfgb

larfpf=: 0 _1 & larfp
larfpfc=: _1 + upd1 larfpf
larfpb=: _1 0 & larfp
larfpbc=: 0 + upd1 larfpb

NB. ---------------------------------------------------------
NB. larftbc
NB.
NB. Description:
NB.   Monad to form the triangular factor T of a block
NB.   reflector H:
NB.     H = H(k-1) * ... * H(1) * H(0) = I - V * T * V' ,
NB.   where T is lower triangular.
NB.
NB. Syntax:
NB.   T=. larftbc eV
NB. where
NB.   eV  - (m+1)×k-matrix (tau,V)
NB.   V   - m×k-matrix, unit upper triangular (trapezoidal)
NB.         with 1s on (k-m)-th diagonal and 0s below
NB.   tau - k-vector τ[0:k-1] corresp. to V
NB.   T   - k×k-matrix, lower triangular

larftbc=: (         0 &{) ([ * (EMPTY (((1 & ( 0:})) @ (0&,) @ (] (,.~) (] mp (([ {. {)~ #)))) ^: (#@[))~ (|.@(1&(}."1))@((* -)~ ((mp ct) {.)\.)))) (ct@(0&(         0 })))
larftbr=: ((< a: ;  0)&{) ([ * (EMPTY (((1 & ( 0:})) @ (0&,) @ (] (,.~) (] mp (([ {. {)~ #)))) ^: (#@[))~ (|.@(1&(}."1))@((* -)~ ((mp ct) {.)\.)))) (    0&((< a: ;  0)}) )
larftfc=: (        _1 &{) ([ * (EMPTY (((1 & (_1:})) @ (,&0) @ (]  ,.   (] mp (([ {. {)~ #)))) ^: (#@[))~ (               (* -)~ ((mp ct) {:)\  ))) (ct@(0&(        _1 })))
larftfr=: ((< a: ; _1)&{) ([ * (EMPTY (((1 & (_1:})) @ (,&0) @ (]  ,.   (] mp (([ {. {)~ #)))) ^: (#@[))~ (               (* -)~ ((mp ct) {:)\  ))) (    0&((< a: ; _1)}) )

NB. ---------------------------------------------------------
NB. larftbr
NB.
NB. Description:
NB.   Monad to form the triangular factor T of a block
NB.   reflector H:
NB.     H = H(k-1) * ... * H(1) * H(0) = I - V' * T * V ,
NB.   where T is lower triangular.
NB.
NB. Syntax:
NB.   T=. larftbr eV
NB. where
NB.   eV  - k×(n+1)-matrix (tau,.V)
NB.   V   - k×n-matrix, unit lower triangular (trapezoidal)
NB.         with 1s on (n-k)-th diagonal and 0s above
NB.   tau - k-vector τ[0:k-1] corresp. to V
NB.   T   - k×k-matrix, lower triangular


NB. ---------------------------------------------------------
NB. larftfc
NB.
NB. Description:
NB.   Monad to form the triangular factor T of a block
NB.   reflector H:
NB.     H = H(0) * H(1) * ... * H(k-1) = I - V * T * V' ,
NB.   where T is upper triangular.
NB.
NB. Syntax:
NB.   T=. larftfc eV
NB. where
NB.   eV  - (m+1)×k-matrix (V,tau)
NB.   V   - m×k-matrix, unit lower triangular (trapezoidal)
NB.         with 1s on 0-th diagonal and 0s above
NB.   tau - k-vector τ[0:k-1] corresp. to V
NB.   T   - k×k-matrix, upper triangular


NB. ---------------------------------------------------------
NB. larftfr
NB.
NB. Description:
NB.   Monad to form the triangular factor T of a block
NB.   reflector H:
NB.     H = H(0) * H(1) * ... * H(k-1) = I - V' * T * V ,
NB.   where T is upper triangular.
NB.
NB. Syntax:
NB.   T=. larftfr eV
NB. where
NB.   eV  - k×(n+1)-matrix (V,.tau)
NB.   V   - k×n-matrix, unit upper triangular (trapezoidal)
NB.         with 1s on 0-th diagonal and 0s below
NB.   tau - k-vector τ[0:k-1] corresp. to V
NB.   T   - k×k-matrix, upper triangular


NB. ---------------------------------------------------------
NB. Verb       Action    Side    Tran   Dir   Layout     eC
NB. larflcbc   H' * C    left    ct     bwd   col-wise   0, C
NB. larflcbr   H' * C    left    ct     bwd   rowwise    0, C
NB. larflcfc   H' * C    left    ct     fwd   col-wise   C, 0
NB. larflcfr   H' * C    left    ct     fwd   rowwise    C, 0
NB. larflnbc   H  * C    left    none   bwd   col-wise   0, C
NB. larflnbr   H  * C    left    none   bwd   rowwise    0, C
NB. larflnfc   H  * C    left    none   fwd   col-wise   C, 0
NB. larflnfr   H  * C    left    none   fwd   rowwise    C, 0
NB. larfrcbc   C  * H'   right   ct     bwd   col-wise   0,.C
NB. larfrcbr   C  * H'   right   ct     bwd   rowwise    0,.C
NB. larfrcfc   C  * H'   right   ct     fwd   col-wise   C,.0
NB. larfrcfr   C  * H'   right   ct     fwd   rowwise    C,.0
NB. larfrnbc   C  * H    right   none   bwd   col-wise   0,.C
NB. larfrnbr   C  * H    right   none   bwd   rowwise    0,.C
NB. larfrnfc   C  * H    right   none   fwd   col-wise   C,.0
NB. larfrnfr   C  * H    right   none   fwd   rowwise    C,.0
NB.
NB. Description:
NB.   Dyads to apply an elementary reflector H or its
NB.   transpose H' to a matrix, from either the left or the
NB.   right. H is defined by pair (v,τ) .
NB.
NB. Syntax:
NB.   eCupd=. vtau larfxxxx eC
NB. where
NB.   eC    - matrix C to update, augmented by trash vector
NB.   vtau  - vector v augmented by scalar τ
NB.   eCupd - being updated matrix C , augmented by modified
NB.           trash vector
NB.   v     - vector with 1 at head (forward direction) or
NB.           tail (backward direction)
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
NB. Verb       Action    Side    Tran   Dir   Layout     eC
NB. larfblcbc  H' * C    left    ct     bwd   col-wise   0, C
NB. larfblcbr  H' * C    left    ct     bwd   rowwise    0, C
NB. larfblcfc  H' * C    left    ct     fwd   col-wise   C, 0
NB. larfblcfr  H' * C    left    ct     fwd   rowwise    C, 0
NB. larfblnbc  H  * C    left    none   bwd   col-wise   0, C
NB. larfblnbr  H  * C    left    none   bwd   rowwise    0, C
NB. larfblnfc  H  * C    left    none   fwd   col-wise   C, 0
NB. larfblnfr  H  * C    left    none   fwd   rowwise    C, 0
NB. larfbrcbc  C  * H'   right   ct     bwd   col-wise   0,.C
NB. larfbrcbr  C  * H'   right   ct     bwd   rowwise    0,.C
NB. larfbrcfc  C  * H'   right   ct     fwd   col-wise   C,.0
NB. larfbrcfr  C  * H'   right   ct     fwd   rowwise    C,.0
NB. larfbrnbc  C  * H    right   none   bwd   col-wise   0,.C
NB. larfbrnbr  C  * H    right   none   bwd   rowwise    0,.C
NB. larfbrnfc  C  * H    right   none   fwd   col-wise   C,.0
NB. larfbrnfr  C  * H    right   none   fwd   rowwise    C,.0
NB.
NB. Description:
NB.   Dyads to build and apply a block reflector H or its
NB.   transpose H' to a matrix, from either the left or the
NB.   right. H is defined by pair (V,Τ) , where Τ is the
NB.   triangular factor produced from pair (V,tau) by
NB.   larftxx .
NB.
NB. Syntax:
NB.   eCupd=. Vtau larfbxxxx eC
NB. where
NB.   eC    - matrix C to update, augmented by trash vector
NB.   Vtau  - matrix V augmented by vector tau
NB.   eCupd - being updated matrix C , augmented by modified
NB.           trash vector
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
NB.
NB. Notes:
NB. - foregoing verbs are tested implicitly in gq, mq, qf
NB.   tests


NB. ---------------------------------------------------------
NB. trefrft
NB.
NB. Description:
NB.   Test algorithms forming the triangular factor of a
NB.   block reflector, by a general matrix
NB.
NB. Syntax:
NB.   trefrft (A;trash)
NB. where
NB.   A - m×n-matrix, is used to get Qf

trefrft=: 3 : 0
  AC=: y=. 0 {:: y
  rcond=. ((_."_)`(norm1 con getri) @. (=/@$)) y

  ('larftbc' tmonad (geqlf`]`(rcond"_)`(_."_)`(_."_))) y
  ('larftbr' tmonad (gerqf`]`(rcond"_)`(_."_)`(_."_))) y
  ('larftfc' tmonad (geqrf`]`(rcond"_)`(_."_)`(_."_))) y
  ('larftfr' tmonad (gelqf`]`(rcond"_)`(_."_)`(_."_))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. trefrfb
NB.
NB. Description:
NB.   Test algorithms applying a block reflector, by a
NB.   general matrix
NB.
NB. Syntax:
NB.   trefrfb (A;C)
NB. where
NB.   A - m×n-matrix, is used to get Qf
NB.   C - m×n-matrix, is used as multiplier

trefrfb=: 3 : 0
  'A C'=. y
  rcond=. (norm1 con getri) C
  'LQf QfL QfR RQf'=. (gelqf ; geqlf ; geqrf ; gerqf) A

  ('larfblcbc' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (QfL;(    C , ~0))
  ('larfblcbr' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (RQf;((ct C), ~0))
  ('larfblcfc' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (QfR;(    C ,  0))
  ('larfblcfr' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (LQf;((ct C),  0))

  ('larfblnbc' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (QfL;(    C , ~0))
  ('larfblnbr' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (RQf;((ct C), ~0))
  ('larfblnfc' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (QfR;(    C ,  0))
  ('larfblnfr' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (LQf;((ct C),  0))

  ('larfbrcbc' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (QfL;((ct C),.~0))
  ('larfbrcbr' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (RQf;(    C ,.~0))
  ('larfbrcfc' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (QfR;((ct C),. 0))
  ('larfbrcfr' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (LQf;(    C ,. 0))

  ('larfbrnbc' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (QfL;((ct C),.~0))
  ('larfbrnbr' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (RQf;(    C ,.~0))
  ('larfbrnfc' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (QfR;((ct C),. 0))
  ('larfbrnfr' tdyad ((0&({::))`(1&({::))`]`(rcond"_)`(_."_)`(_."_))) (LQf;(    C ,. 0))

  EMPTY
)

NB. ---------------------------------------------------------
NB. testref
NB.
NB. Description:
NB.   Adv. to make verb to test block reflection algorithms
NB.   by matrices of generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testref
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.            mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testref 150 100
NB.
NB. Notes:
NB. - non-blocked larfxxxx algos are tested implicitly in gq,
NB.   mq, qf tests

testref=: 1 : 'EMPTY [ (trefrfb [ trefrft) @ (u ; u)'
