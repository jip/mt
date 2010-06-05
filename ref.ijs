NB. ref.ijs
NB. Reflections
NB.
NB. larfg   Generate an elementary reflector
NB. larfl   Apply an elementary reflector to a matrix from
NB.         the left
NB. larfr   Apply an elementary reflector to a matrix from
NB.         the right
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. larfg
NB. Generates an elementary reflector H such that
NB. H'*x = β*e1, where x=(α,x[1],...,x[n-1]), H=I-τ*v*v', H'*H=I,
NB. τ=(β-α)/β, v=(x-β*e1)/(α-β)=(1,v[1],...,v[n-1]), x,v ∊ ℂⁿ,
NB. α,τ ∊ ℂ, β ∊ ℝ
NB.
NB. Syntax:
NB.   y=. larfg x
NB. where
NB.   y=(β,v[1],...,v[n-1],τ)
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

NB. emulate LAPACK's xLARFG

larfg=: 3 : 0
  alpha=. {. y
  xnorm=. norms }. y                                        NB. ||x-α*e1||_2
  if. (0 = xnorm) *. (0 = 11 o. alpha) do.
    y=. y , 0                                               NB. append in-place τ=0
  else.
    beta=. - (9 o. alpha) condneg norms alpha , xnorm       NB. β=-copysign(||x||_2,Re(α))
    y=. y , beta                                            NB. append in-place β
    if. FP_SFMIN > | beta do.
      y=. y % FP_SFMIN                                      NB. scale (α,x[1],...,x[n-1],β)
    end.
    tau=. ({: (- % [) {.) y                                 NB. τ=(β_scaled-α_scaled)/β_scaled
    y=. (% ({. - {:)) y                                     NB. x=x/(α_scaled-β_scaled)
    y=. (beta , tau) 0 _1 } y                               NB. replace in-place α_scaled by β_unscaled and β_scaled by τ
  end.
)

NB. emulate LAPACK's xLARFP

larfp=: 3 : 0
yold=. y
  alpha=. {. y
  xnorm=. norms }. y                                        NB. ||x-α*e1||_2
NB.  if. (0 = xnorm) do.
NB.    y=. (| alpha) 0 } y                                     NB. replace in-place α by β
NB.    y=. y , (1 - * alpha)                                   NB. append in-place τ=1-α/|α|
NB.  else.
    beta=. - (9 o. alpha) condneg norms alpha , xnorm       NB. β=-copysign(||x||_2,Re(α))
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
      'realpha imalpha'=. +. {. y                           NB. Re(α_scaled) , Im(α_scaled)
      gamma=. realpha + {: y                                NB. γ=Re(α_scaled)+|β_scaled|
      delta=. (imalpha , xnorm) (- @ (+/) @ ([ * %)) gamma  NB. δ=-(Im(α_scaled)*(Im(α_scaled)/γ)+||x-α*e1||_2*(||x-α*e1||_2/γ))
      dzeta=. delta j. imalpha                              NB. ζ=δ+i*Im(α_scaled)
      tau=. - dzeta % {: y                                  NB. τ=-ζ/|β_scaled|
    end.
    y=. y % dzeta
    y=. (beta , tau) 0 _1 } y                               NB. replace in-place α_scaled by |β_unscaled| and |β_scaled| by τ
smoutput 'LARFP' ; '[a;x]' ; yold ; '[b;v;t]' ; y
NB.  end.
)

NB. verb_to_apply_H=: larflr larf1 cios

larf1=: 2 : '((u ((cios2ios n) { ]))`(cios2ios @: (n"_))`]) }'

NB. H=. tauv2H tau ; v

tauv2H=: 3 : 0
  'tau v'=. y
  n=. # v
  (idmat n) - tau * v */ + v
)

testlarfg=: 3 : 0
  try.
    'tau_enh y_enh'=. ({: ; }:) larfp y
NB. smoutput 'tau_enh' ; (0j313 ": ,. +. tau_enh)
    beta_enh=. {. y_enh
    v_enh=. 1 (0) } y_enh
NB. smoutput 'v_enh' ; (0j313 ": ,. , +. v_enh)
    H_enh=. tauv2H tau_enh ; v_enh
    Y_enh=. (ct H_enh) mp y
NB. smoutput 'Y_enh' ; (0j313 ": ,. , +. Y_enh)
    error_enh=. (}. Y_enh) (% & norms) y
  catch.
    beta_enh=. _.
    error_enh=. _.
  end.

  try.
    'tau_enh2 y_enh2'=. ({: ; }:) larfp2 y
NB. smoutput 'tau_enh2' ; (0j313 ": ,. +. tau_enh2)
    beta_enh2=. {. y_enh2
    v_enh2=. 1 (0) } y_enh2
NB. smoutput 'v_enh2' ; (0j313 ": ,. , +. v_enh2)
    H_enh2=. tauv2H tau_enh2 ; v_enh2
    Y_enh2=. (ct H_enh2) mp y
NB. smoutput 'Y_enh2' ; (0j313 ": ,. , +. Y_enh2)
    error_enh2=. (}. Y_enh2) (% & norms) y
  catch.
    beta_enh2=. _.
    error_enh2=. _.
  end.

  error_enh ([,],(|@-)) error_enh2
)

NB. Template adverb to make verbs applying an elementary
NB. reflector H = I - tau*v*v' to a matrix A
NB. larfx=. v1 larf
NB. Note: consider remove condition since is called only when τ≠0

NB. larf=: 1 : '(] - ({: @ [) * ((u }:)~)) ^: (0 ~: {: @ [)'  NB. check τ≠0, no assign v[0]=1
larf=: 1 : '] - ({: @ [) * ((u (1 & (0 }) @ }:))~)'  NB. no check τ≠0, assign v[0]=1

NB. apply an elementary reflector H from left:
NB.   H*A = A - tau * v * (v' * A)
NB. emulate LAPACK's xLARF('L')
NB. HbyA=. (v,τ) larfl A

larfl=: (] */ (mp~ +)) larf

NB. apply an elementary reflector H from right:
NB.   A*H = A - tau * (A * v) * v'
NB. emulate LAPACK's xLARF('R')
NB. AbyH=. (v,τ) larfr A

larfr=: (mp */ (+@])) larf

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB.*tlarfg v test larfg

tlarfg=: 3 : 0
)

NB. ---------------------------------------------------------
NB.*tlarf v test larf

tlarf=: 3 : 0
)

NB. ---------------------------------------------------------
NB.*testref v test reflectors

testref=: 3 : 0
)

NB. ((# #: i.@<.@(^~)) #) { ]
NB. generate all allocations of length x from vector y
