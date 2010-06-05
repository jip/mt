NB. ref.ijs
NB. Reflections
NB.
NB. larfg      Generate an elementary reflector
NB. larfp      Generate an elementary reflector with β≥0
NB.
NB. larfl      Multiply matrix by vectors from the left
NB. larfr      Multiply matrix by vectors from the right
NB.
NB. larftbc    Monad to form the triangular factor of a block
NB.            reflector from input vectors which have
NB.            backward direction and columnwise layout
NB. larftbr    Monad to form the triangular factor of a block
NB.            reflector from input vectors which have
NB.            backward direction and rowwise layout
NB. larftfc    Monad to form the triangular factor of a block
NB.            reflector from input vectors which have
NB.            forward direction and columnwise layout
NB. larftfr    Monad to form the triangular factor of a block
NB.            reflector from input vectors which have
NB.            forward direction and rowwise layout
NB.
NB. larfblcbc  Dyad to build a block reflector by larftbc and
NB.            to apply its transpose to a submatrix from the
NB.            left
NB. larfblcfc  Dyad to build a block reflector by larftfc and
NB.            to apply its transpose to a submatrix from the
NB.            left
NB. larfbrcfr  Dyad to build a block reflector by larftfr and
NB.            to apply it to a submatrix from the right
NB. larfbrnbr  Dyad to build a block reflector by larftbr and
NB.            to apply it to a submatrix from the right
NB. larfbrnfr  Dyad to build a block reflector by larftfr and
NB.            to apply it to a submatrix from the right
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. return y without elements with linear IOS from x
ywolx=: (({,)~ (<^:3))~

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
NB.      a) extract vector p[i]: pi=P[i,0:i-1]
NB.      b) calc new column: c[i]: ci=. pi mp Ti
NB.      c) stitch c[i] to T[i]: Ti=. Ti vTi ci
NB.      d) append zero row to T[i]: Ti=. (0 & va0) Ti
NB.      e) write 1 in the diagonal element of appended zero
NB.         row: Ti=. (1 & vs1) Ti
NB.   6) multiply T[i] by tau (item-by-row): T=. tau * Ti
NB.
NB. Notes:
NB. - emulates LAPACK's LARFT, but without zeros scanning in
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
NB.   H'*y = β*e1, where H=I-τ*v*v', H'*H=I. H is represented
NB.   in factored form by n-vector v and scalar τ. Input and
NB.   output arrays have the same shape and may be either a
NB.   (n+1)-vector or (n+1)×1-matrix or 1×(n+1)-matrix.
NB.   Vector v can have either forward (α in head) or
NB.   backward (α in tail) direction. Input and output vector
NB.   directions are the same.
NB.
NB. Syntax:
NB.   z=. ios larfg y
NB.   z=. ios larfp y
NB. where
NB.   ios - 2-vector of integers (ioa,iot)
NB.   ioa - lIO α in y
NB.   iot - lIO pre-allocated zero scalar for τ in y
NB.   y   - (n+1)-vector or (n+1)×1-matrix or 1×(n+1)-matrix
NB.         having scalar α at index ioa, scalar 0 at index
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
NB. - define monad to get input: (α x[1:n-1] 0), and to
NB.   output: (β v[1:n-1] τ):
NB.     larfgf=: 0 _1 & larfg
NB. - define monad to get input: (0 x[1:n-1] α), and to
NB.   output: (τ v[1:n-1] β):
NB.     larfgb=: _1 0 & larfg
NB. - define monad to get input: (α x[1:n-1] 0), and to
NB.   output: (β v[1:n-1] conj(τ)):
NB.     larfgfc=: (+ updl _1) @ larfgf
NB. - define monad to get input: (0 x[1:n-1] α), and to
NB.   output: (conj(τ) v[1:n-1] β):
NB.     larfgbc=: (+ updl 0) @ larfgb
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
  xnorm=. norms (x ywolx y)                                 NB. ||x||_2
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
  end.
  y
)

larfp=: 4 : 0
  'ioalpha iotau'=. x
  alpha=. ioalpha ({,) y
  xnorm=. norms (x ywolx y)                                 NB. ||x||_2
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
NB. larfl
NB. Multiply matrix A by vectors from the left:
NB.   B=. v * (v' * A)
NB.
NB. Syntax:
NB.   B=. A larfl v
NB. where
NB.   A - m×n-matrix
NB.   v - m-vector
NB.   B - m×n-matrix

larfl=: ] */ (mp~ +)

NB. ---------------------------------------------------------
NB. larfr
NB. Multiply matrix A by vectors from the right:
NB.   B=. (A * v) * v'
NB.
NB. Syntax:
NB.   B=. A larfr v
NB. where
NB.   A - m×n-matrix
NB.   v - n-vector
NB.   B - m×n-matrix

larfr=: mp */ (+@])

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
NB.     H = H(1) H(2) ... H(k) and T is upper triangular,
NB.   otherwise direction is backward, and:
NB.     H = H(k) ... H(2) H(1) and T is lower triangular.
NB.   If layout is columnwise, then:
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
NB. Verb       Action   Side   Transp  Direction  Layout      Used in
NB. larfblcbc  H' * C   left   ct      backward   columnwise  geqlf
NB. larfblcbr  H' * C   left   ct      backward   rowwise
NB. larfblcfc  H' * C   left   ct      forward    columnwise  geqrf
NB. larfblcfr  H' * C   left   ct      forward    rowwise
NB. larfblnbc  H  * C   left   none    backward   columnwise
NB. larfblnbr  H  * C   left   none    backward   rowwise
NB. larfblnfc  H  * C   left   none    forward    columnwise
NB. larfblnfr  H  * C   left   none    forward    rowwise
NB. larfbrcbc  C  * H'  right  ct      backward   columnwise
NB. larfbrcbr  C  * H'  right  ct      backward   rowwise
NB. larfbrcfc  C  * H'  right  ct      forward    columnwise
NB. larfbrcfr  C  * H'  right  ct      forward    rowwise     unglq
NB. larfbrnbc  C  * H   right  none    backward   columnwise
NB. larfbrnbr  C  * H   right  none    backward   rowwise     gerqf
NB. larfbrnfc  C  * H   right  none    forward    columnwise
NB. larfbrnfr  C  * H   right  none    forward    rowwise     gelqf
NB.
NB. Description:
NB.   Dyads to build and apply a block reflector H or its
NB.   transpose H' to a submatrix C, from either the left or
NB.   the right.
NB.
NB. Syntax:
NB.   eCupd=. eC larfbxxxx Vtau
NB. where
NB.   Vtau  - matrix V with appended or stitched vector tau
NB.   V     - unit triangular (trapezoidal) matrix
NB.   tau   - k-vector τ[0:k-1] corresp. to V
NB.   eC    - matrix C to update with appended or stitched
NB.           trash vector
NB.   eCupd - being updated matrix C with appended or
NB.           stitched trash vector
NB.
NB. Algorithm:
NB.   1) compute negated T: Tneg=. - larftxx Vtau
NB.   2) compute negated part of H without I added:
NB.        Hnegpart = op1(V) * op2(Tneg) * op1(V')
NB.      where op(Z) is either Z or Z'
NB.   3) compute op2(H) = (I+Hnegpart) by incrementing the
NB.      0-th diagonal:
NB.        op2(H) =: (>: upddiag0) Hnegpart
NB.
NB. Notes:
NB. - emulates LAPACK's sequence of calls to LARFT and then
NB.   to LARFB

larfblcbc=: [ (mp~) (((0&(IOSFR }))@]) (    [  ((>: upddiag0)@ mp  ) (ct@ mp  )) (-@larftbc@]))  NB. (I + V * (V * (-T))') * C
larfblcfc=: [ (mp~) (((0&(IOSLR }))@]) (    [  ((>: upddiag0)@ mp  ) (ct@ mp  )) (-@larftfc@]))  NB. (I + V * (V * (-T))') * C
larfbrcfr=: [  mp   (((0&(IOSLC }))@]) (    [  ((>: upddiag0)@(mp~)) (ct@(mp~))) (-@larftfr@]))  NB. C  * (I + ((-T) * V)' * V)
larfbrnbr=: [  mp   (((0&(IOSFC }))@]) ((ct@[) ((>: upddiag0)@ mp  ) (    mp~ )) (-@larftbr@]))  NB. C * (I + V' * (-T) * V)
larfbrnfr=: [  mp   (((0&(IOSLC }))@]) ((ct@[) ((>: upddiag0)@ mp  ) (    mp~ )) (-@larftfr@]))  NB. C * (I + V' * (-T) * V)

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
