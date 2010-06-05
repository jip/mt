NB. ---------------------------------------------------------
NB. larfg
NB. larfp
NB.
NB. Generate an elementary reflector H of order n such that
NB. H'*y = β*e1, where H=I-τ*v*v', H'*H=I
NB.
NB. Syntax:
NB.   z=. larfg y
NB.   z=. larfp y
NB. where
NB.   y - n-vector or n×1-matrix (α,x[1:n-1]) to reflect,
NB.       α ∊ ℂ, x ∊ ℂⁿ⁻¹
NB.   z - (n+1)-vector or (n+1)×1-matrix (y and z ranks are
NB.       match) (β,v[1:n-1],τ) of reflection result,
NB.       β ∊ ℝ, (β≥0 for larfp only), v ∊ ℂⁿ⁻¹, τ ∊ ℂ
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
NB. - IEEE fp cfg features are implicitly encoded here
NB. - larfg emulates LAPACK's xLARFG
NB. - larfp emulates LAPACK's xLARFP
NB. - both larfg (<-CHECKME) and larfp provide: τ=2 ↔ ||x||_2 = 0
NB. - not zeroing v in larfp in case τ=2 (x≈0) relies on
NB.   fact: 'comparison tolerance tol>0'; otherwise
NB.   (tol=0) x should be filled by zeros (see larf* gerf*)

larfg=: 3 : 0
  alpha=. {. y
  xnorm=. norms }. y                                        NB. ||x||_2
  if. (0 = xnorm) *. (0 = 11 o. alpha) do.
    y=. y , 0                                               NB. append in-place τ=0
  else.
    beta=. - (9 o. alpha) condneg norms alpha , xnorm       NB. β=-copysign(||y||_2,Re(α))
    y=. y , beta                                            NB. append in-place β
    if. FP_SFMIN > | beta do.
      y=. y % FP_SFMIN                                      NB. scale (α,x[1],...,x[n-1],β)
    end.
    tau=. ({: (- % [) {.) y                                 NB. τ=(β_scaled-α_scaled)/β_scaled
    y=. (% ({. ({.@:-) {:)) y                               NB. y=y/(α_scaled-β_scaled)
    y=. (beta , tau) (0 _1 " _) } y                         NB. replace α_scaled by β_unscaled and β_scaled by τ
  end.
)

larfp=: 3 : 0
  alpha=. {. y
  xnorm=. norms }. y                                        NB. ||x||_2
  if. (0 = xnorm) do.
    y=. (| alpha) 0 } y                                     NB. replace in-place α by β
    y=. y , (1 - * alpha)                                   NB. append in-place τ=1-α/|α|
  else.
    beta=. - (9 o. alpha) condneg norms alpha , xnorm       NB. β=-copysign(||y||_2,Re(α))
    y=. y , | beta                                          NB. append in-place |β|
    if. FP_SFMIN > | beta do.
      y=. y % FP_SFMIN                                      NB. scale (α,x[1],...,x[n-1],|β|)
      xnorm=. xnorm % FP_SFMIN
    end.
    if. 0 <: beta do.
      dzeta=. ({. - {:) y                                   NB. ζ=α_scaled-|β_scaled|
      tau=. ({: (- % [) {.) y                               NB. τ=(|β_scaled|-α_scaled)/|β_scaled|
    else.
      beta=. - beta                                         NB. |β_unscaled|
      'realpha imalpha'=. +. 0 ({,) y                       NB. Re(α_scaled) , Im(α_scaled)
      gamma=. realpha + _1 ({,) y                           NB. γ=Re(α_scaled)+|β_scaled|
      delta=. (imalpha , xnorm) (- @ (+/) @ ([ * %)) gamma  NB. δ=-(Im(α_scaled)*(Im(α_scaled)/γ)+||x||_2*(||x||_2/γ))
      dzeta=. delta j. imalpha                              NB. ζ=δ+i*Im(α_scaled)
      tau=. - dzeta % {: y                                  NB. τ=-ζ/|β_scaled|
    end.
    y=. y % dzeta
    y=. (beta , tau) (0 _1 " _) } y                         NB. replace α_scaled by |β_unscaled| and |β_scaled| by τ
  end.
)
