NB. ref.ijs
NB. Reflections
NB.
NB. larfgf   Generate an elementary reflector with forward
NB.          direction
NB. larfgb   Generate an elementary reflector with backward
NB.          direction
NB. larfgfc  Generate an elementary reflector with forward
NB.          direction and conjugated τ
NB. larfgbc  Generate an elementary reflector with backward
NB.          direction and conjugated τ
NB. larfpf   Generate an elementary reflector with forward
NB.          direction and β≥0
NB. larfpb   Generate an elementary reflector with backward
NB.          direction and β≥0
NB. larfpfc  Generate an elementary reflector with forward
NB.          direction, β≥0 and conjugated τ
NB. larfpbc  Generate an elementary reflector with backward
NB.          direction, β≥0 and conjugated τ
NB.
NB. larflsc  Template adv. to make verbs to apply an
NB.          elementary reflector to a matrix from the left,
NB.          τ is pre-conjugated
NB. larfrcs  Template adv. to make verbs to apply an
NB.          elementary reflector to a matrix from the right,
NB.          v is pre-conjugated
NB.########
NB. larfL    Template adv. to make verbs to update submatrix
NB.          by larfl
NB. larfR    Template adv. to make verbs to update submatrix
NB.          by larfr
NB. larfRL   Template adv. to make verbs to update one
NB.          submatrix by larfr, then another one by larfl
NB. larfLs   Template adv. to make verbs to update shrinked
NB.          submatrix by larfl
NB. larfRs   Template adv. to make verbs to update shrinked
NB.          submatrix by larfr
NB. larfRLs  Template adv. to make verbs to update one
NB.          shrinked submatrix by larfr, then another one by
NB.          larfl
NB.
NB. gerf0    Template conj. to make verbs to generate and
NB.          conditionally apply an elementary reflector to a
NB.          matrix
NB. gerf02   Template conj. to make verbs to generate and
NB.          conditionally apply an elementary reflector to a
NB.          matrix
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. update n-th from x (warning!: linear IOS) by (x u y)
unxu=: 2 : '((n({,)[) - u)(n"_)}['

NB. return y without elements with linear IOS from x
ywolx=: (({,)~ (<^:3))~

NB. ---------------------------------------------------------
NB. - extract n-th cIOS cios from x, then extract by fromci
NB.   from y submatrix defined by cios
NB. - apply u
NB. - make imaginary
cfuj=: 2 : 'j. @ u @ (n fromci)'

NB. ---------------------------------------------------------
NB. larfg
NB. larfp
NB.
NB. Generate an elementary reflector H of order n such that
NB. H'*y = β*e1, where H=I-τ*v*v', H'*H=I
NB.
NB. Syntax:
NB.   z=. ios larfg y
NB.   z=. ios larfp y
NB. where
NB.   ios - 2-vector of integers (ioa,iot)
NB.   ioa - IO α in y
NB.   iot - IO pre-allocated zero scalar for τ in y
NB.   y   - (n+1)-vector or (n+1)×1-matrix or 1×(n+1)-matrix
NB.         having scalar α at index ioa, scalar 0 at index
NB.         iot, and vector x[1:n-1] in the rest elements,
NB.         vector to reflect is (α,x[1],...,x[n-1]), α ∊ ℂ,
NB.         x ∊ ℂⁿ⁻¹
NB.   z   - (n+1)-vector or (n+1)×1-matrix or 1×(n+1)-matrix
NB.         (y and z ranks are match) having scalar β at
NB.         index ioa, scalar τ at index iot, and vector
NB.         v[1:n-1] in the rest elements, reflection result
NB.         is vector (1,v[1],...,v[n-1]), 1 is not stored,
NB.         β ∊ ℝ, (β≥0 for larfp only), v ∊ ℂⁿ⁻¹, τ ∊ ℂ
NB.
NB. Applications:
NB. - reflect vector (α,x) by larfg and store τ at tail:
NB.   z=. 0 _1 larfg (α,x,0)
NB.   v=. 1 , (0 _1 ywox z)
NB.   'beta tau'=. 0 _1 { z
NB. - reflect vector (x,α) by larfp and store τ at head:
NB.   z=. _1 0 larfp (0,x,α)
NB.   v=. 1 ,~ (0 _1 ywox z)
NB.   'beta tau'=. _1 0 { z
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
NB. - larfp provides τ=2 ↔ ||x||_2 = 0
NB. - not zeroing v in larfp in case τ=2 (x≈0) relies on
NB.   fact: 'comparison tolerance tol>0'; otherwise
NB.   (tol=0) x should be filled by zeros (see larf* gerf*)

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
NB. larf
NB. Template conj. to make verbs to apply an elementary
NB. reflector H = I - τ*v*v' to a matrix A from either the
NB. left or the right
NB.
NB. Syntax:
NB.   vapp=. (vmul`vv`vtau) larf (iob,iot)
NB. where
NB.   vmul - either larfl or larfr; is called as: (A vmul v)
NB.   vv   - verb to pre-process v; is called as: (vv v)
NB.   vtau - verb to pre-process τ; is called as: (vtau τ)
NB.   vapp - verb to apply an elementary reflector; is called
NB.          as: (z vapp A)
NB.   A    - m×n-matrix
NB.   z    - either (k+1)-vector, or 1×(k+1)-matrix, or
NB.          (k+1)×1-matrix, contains scalar β at index iob,
NB.          scalar τ at index iot, and all but one elements
NB.          of vector v in the rest elements
NB.   v    - k-vector, is produced from z by removing τ and
NB.          replacing β by value 1
NB.   k    = either m (if vmul is larfl) or n (if vmul is
NB.          larfr)
NB.
NB. Notes:
NB. - pre-ravel z
NB. - no check τ≠0
NB. - assign (1 iob } z)

larf=: 2 : '(([ - (((2{u)`:6) @ (({:n) { ])) * ((((0{u)`:6) (1 & ((({.n)"_) }) @: ((1{u)`:6) @: (({:n) & ywolx))))) ,)~'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. larfgf
NB. larfgb
NB. larfgfc
NB. larfgbc
NB. larfpf
NB. larfpb
NB. larfpfc
NB. larfpbc
NB.
NB. Shortcuts to generate an elementary reflector H of order
NB. n such that H'*y = β*e1, where H=I-τ*v*v', H'*H=I. H is
NB. represented in factored form by n-vector v and scalar τ.
NB. Vector v can have either forward (α in head) or backward
NB. (α in tail) direction. Input and output vector directions
NB. are the same.

larfgf=: 0 _1 & larfg           NB. (β,v[1:n-1],τ)=. larfgf (α,x[1:n-1],0)
larfgb=: _1 0 & larfg           NB. (τ,v[1:n-1],β)=. larfgb (0,x[1:n-1],α)
larfgfc=: (+ updl _1) @ larfgf  NB. (β,v[1:n-1],conj(τ))=. larfgfc (α,x[1:n-1],0)
larfgbc=: (+ updl 0) @ larfgb   NB. (conj(τ),v[1:n-1],β)=. larfgbc (0,x[1:n-1],α)

larfpf=: 0 _1 & larfp           NB. (β,v[1:n-1],τ)=. larfpf (α,x[1:n-1],0)
larfpb=: _1 0 & larfp           NB. (τ,v[1:n-1],β)=. larfpb (0,x[1:n-1],α)
larfpfc=: (+ updl _1) @ larfpf  NB. (β,v[1:n-1],conj(τ))=. larfpfc (α,x[1:n-1],0)
larfpbc=: (+ updl 0) @ larfpb   NB. (conj(τ),v[1:n-1],β)=. larfpbc (0,x[1:n-1],α)

NB. ---------------------------------------------------------
NB. larflsc
NB. larflss
NB. Template adv. to make verbs to apply an elementary
NB. reflector H to a matrix A from the left: ########
NB.   H*A = A - conj(τ) * v * (v' * A)
NB.
NB. Syntax:
NB.   vapp=. (iob,iot) larflsc
NB. where
NB.   iob  - integer, IO β in z
NB.   iot  - integer, IO τ in z
NB.   vapp - verb to apply an elementary reflector; is called
NB.          as: (HA=. z vapp A)
NB.   z    - either (m+1)-vector, or 1×(m+1)-matrix, or
NB.          (m+1)×1-matrix, contains scalar β at index iob,
NB.          scalar τ at index iot, and all but one elements
NB.          of vector v in the rest elements
NB.   A    - m×n-matrix
NB.   v    - m-vector, is produced from z by removing τ and
NB.          replacing β by value 1
NB.   HA   - m×n-matrix H'*A
NB.   H    - m×m-matrix H
NB.   m    ≥ 0
NB.   n    ≥ 0
NB.
NB. Notes:
NB. - differs from LAPACK version 3.2's xLARF('L') in not
NB.   scanning trailing zeros in v

larflsc=: (larfl`]`+) larf  NB. pre-processing: the same v, conjugated τ
larflss=: (larfl`]`]) larf  NB. pre-processing: the same v, the same τ

NB. ---------------------------------------------------------
NB. larfrss
NB. larfrcs
NB. Template adv. to make verbs to apply an elementary
NB. reflector H to a matrix A from the right:
NB.   A*H = A - τ * (A * v) * v'
NB.
NB. Syntax:
NB.   vapp=. (iob,iot) larfrss
NB.   vapp=. (iob,iot) larfrcs
NB. where
NB.   iob  - integer, IO β in z
NB.   iot  - integer, IO τ in z
NB.   vapp - verb to apply an elementary reflector; is called
NB.          as: (AH=. z vapp A)
NB.   z    - either (n+1)-vector, or 1×(n+1)-matrix, or
NB.          (n+1)×1-matrix, contains scalar β at index iob,
NB.          scalar τ at index iot, and all but one elements
NB.          of vector v in the rest elements
NB.   A    - m×n-matrix
NB.   v    - n-vector, is produced from z by removing τ and
NB.          replacing β by value 1
NB.   HA   - m×n-matrix A*H
NB.   H    - n×n-matrix H
NB.   m    ≥ 0
NB.   n    ≥ 0
NB.
NB. Notes:
NB. - differs from LAPACK version 3.2's xLARF('R') in not
NB.   scanning trailing zeros in v

larfrcc=: (larfr`+`+) larf  NB. pre-processing: conjugate v, conjugate τ
larfrcs=: (larfr`+`]) larf  NB. pre-processing: conjugate v, the same τ
NB. larfrss=: (larfr`]`]) larf  NB. pre-processing: the same v, same τ

NB. ---------------------------------------------------------
NB. Adverb:     Syntax:                Apply:
NB. larfLbcc    vapp=. ios larfLbcc    conjugate(H(τ,v,β)) * Asub
NB. larfLbcs    vapp=. ios larfLbcs    transpose(H(τ,v,β)) * Asub
NB. larfLbsc    vapp=. ios larfLbsc    H(τ,v,β)' * Asub
NB. larfLbss    vapp=. ios larfLbss    H(τ,v,β) * Asub
NB. larfLfcc    vapp=. ios larfLfcc    conjugate(H(β,v,τ)) * Asub
NB. larfLfcs    vapp=. ios larfLfcs    transpose(H(β,v,τ)) * Asub
NB. larfLfsc    vapp=. ios larfLfsc    H(β,v,τ)' * Asub
NB. larfLfss    vapp=. ios larfLfss    H(β,v,τ) * Asub
NB. larfRfcc    vapp=. ios larfRfcs    Asub * conjugate(H(β,v,τ))
NB. larfRfcs    vapp=. ios larfRfcs    Asub * transpose(H(β,v,τ))
NB. larfRfsc    vapp=. ios larfRfsc    Asub * H(β,v,τ)'
NB. larfRfss    vapp=. ios larfRfss    Asub * H(β,v,τ)
NB. larfRbcc    vapp=. ios larfRbcs    Asub * conjugate(H(τ,v,β))  Q: cungl2 says A*H'
NB. larfRbcs    vapp=. ios larfRbcs    Asub * transpose(H(τ,v,β))
NB. larfRbsc    vapp=. ios larfRbsc    Asub * H(τ,v,β)'
NB. larfRbss    vapp=. ios larfRbss    Asub * H(τ,v,β)
NB.
NB. Template adv. to make verbs to update submatrix by an
NB. elementary reflector from either the left or the right.
NB. Reflector's layout is either a row or a column.
NB. Reflector's direction is either forward (β in the head)
NB. or backward (β in the tail). Vector v or scalar τ may be
NB. pre-conjugated or not independently of each other.
NB.
NB. Syntax:
NB.   Aupd=. cios vapp A
NB. ##########vapp=. vref larfS (ioz,ios,iob,iot)
NB. where
NB.   vref    - verb to apply an elementary reflector; is
NB.             called as: (subAupd=. z vref subA)
NB.   ioz     - integer, IO vector z's cIOS in cIOS bundle
NB.             cios
NB.   ios     - integer, IO submatrix subA's cIOS in cIOS
NB.             bundle cios
NB.   iob     - integer, IO scalar β in vector z
NB.   iot     - integer, IO scalar τ in vector z
NB.   vapp    - verb to apply an elementary reflector to
NB.             submatrix from either the left or the right;
NB.             is called as: ()
NB.   z       - either (l+1)-vector, or 1×(l+1)-matrix, or
NB.             (k+1)×1-matrix, contains scalar β at index iob,
NB.             scalar τ at index iot, and all but one elements
NB.             of vector v in the rest elements
NB.   cios    - k×2 matrix of cIOSs
NB.   A       - m×n-matrix
NB.   subA    - m1×n1-matrix to update, solid part of A
NB.   Aupd    - m×n-matrix, A with updated subA
NB.   subAupd - m1×n1-matrix, updated subA, solid part of Aupd
NB.   k       ≥ max(iosz,ioss)+1
NB.   m       ≥ 0
NB.   n       ≥ 0

larfLfsc=: (0 _1 larflsc) upd2ci  NB. from the left; forward; the same v; conjugated τ
larfLfss=: (0 _1 larflss) upd2ci  NB. from the left; forward; the same v; the same τ

larfLbsc=: (_1 0 larflsc) upd2ci  NB. from the left; backward; the same v; conjugated τ
larfLbss=: (_1 0 larflss) upd2ci  NB. from the left; backward; the same v; conjugated τ

larfRfcc=: (0 _1 larfrcc) upd2ci  NB. from the right; forward; conjugated v; conjugated τ
larfRfcs=: (0 _1 larfrcs) upd2ci  NB. from the right; forward; conjugated v; the same τ

larfRbcc=: (_1 0 larfrcc) upd2ci  NB. from the right; backward; conjugated v; conjugated τ
larfRbcs=: (_1 0 larfrcs) upd2ci  NB. from the right; backward; conjugated v; the same τ

NB. ######---------------------------------------------------------
NB. larfL
NB. Template adv. to make verbs to update submatrix by larfl
NB.
NB. Syntax:
NB.   vapp=. iosv larfL ioss
NB. where
NB.   iosv - 2-vector of integers (iosZ,iosL)
NB.   ioss - 2-vector of integers (ioa,iot)
NB.   vapp - verb to update submatrix L by larfl; usage:
NB.            Aupd=. cios vapp A
NB.          where
NB.            A    - m×n-matrix
NB.            cios - k×2 matrix of cIOSs (see layout below)
NB.            Aupd - m×n-matrix, A with updated submatrix L
NB.            k    ≥ max(iosZ,iosL)+1
NB.            m    ≥ 0
NB.            n    ≥ 0
NB.
NB. Storage layout for cios:
NB.   iosZ{cios - cIOS of vector z (see larfl)
NB.   iosL{cios - cIOS of submatrix L in matrix A to update
NB.               by larfl
NB.
NB. Application:
NB. - let cios=. (ciosYZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y and z (see larfl), scalar
NB.   τ, and submatrices R, L, respectively; scalars α and τ
NB.   are in head and tail, respectively, of vector z; then
NB.   to apply larfl to submatrix L:
NB.     Aupd=. cios (1 4 larfL 0 _1) A


NB. ---------------------------------------------------------
NB. larfR
NB. Template adv. to make verbs to update submatrix by larfr
NB.
NB. Syntax:
NB.   vapp=. ios larfR
NB. where
NB.   ios  - 2-vector of integers (iosZ,iosR)
NB.   vapp - verb to update submatrix R by larfr; usage:
NB.            Aupd=. cios vapp A
NB.          where
NB.            A    - m×n-matrix
NB.            cios - k×2 matrix of cIOSs (see layout below)
NB.            Aupd - m×n-matrix, A with updated submatrix R
NB.            k    ≥ max(iosZ,iosR)+1
NB.            m    ≥ 0
NB.            n    ≥ 0
NB.
NB. Storage layout for cios:
NB.   iosZ{cios - cIOS of vector z (see larfr)
NB.   iosR{cios - cIOS of submatrix R in matrix A to update
NB.               by larfr
NB.
NB. Application:
NB. - let cios=. (ciosY , ciosZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y, z (see larfr), scalar τ,
NB.   and submatrices R, L, respectively; then to apply larfr
NB.   to submatrix R:
NB.     Aupd=. cios (1 3 larfR) A

NB. ---------------------------------------------------------
NB. larfRL
NB. Template adv. to make verbs to update one submatrix by
NB. larfr, then another one by larfl
NB.
NB. Syntax:
NB.   vapp=. ios larfRL
NB. where
NB.   ios  - 3-vector of integers (iosZ,iosR,iosL)
NB.   vapp - verb to update submatrix R by larfr, then
NB.          submatrix L by larfl; usage:
NB.            Aupd=. cios vapp A
NB.          where
NB.            A    - m×n-matrix
NB.            cios - k×2 matrix of cIOSs (see layout below)
NB.            Aupd - m×n-matrix, A with updated submatrices
NB.                   R and L
NB.            k    ≥ max(iosZ,iosR,iosL)+1
NB.            m    ≥ 0
NB.            n    ≥ 0
NB.
NB. Storage layout for cios:
NB.   iosZ{cios - cIOS of vector z (see larfr or larfl)
NB.   iosR{cios - cIOS of submatrix R in matrix A to update
NB.               by larfr
NB.   iosL{cios - cIOS of submatrix L in matrix A to update
NB.               by larfl
NB.
NB. Application:
NB. - let cios=. (ciosY , ciosZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y, z (see larfr), scalar τ,
NB.   and submatrices R, L, respectively; then to apply larfr
NB.   to submatrix R and then larfl to submatrix L:
NB.     Aupd=. cios (1 3 4 larfRL) A

NB. larfRL=: 1 : '[ ((0 2{m) larfL 0 _1) ((0 1{m) larfR 0 _1)'

NB. ---------------------------------------------------------
NB. larfLs
NB. Template adv. to make verbs to update shrinked submatrix
NB. by larfl. If m1-vector v has zr trailing zeros, and
NB. (m1-zr)×n1-submatrix L[0:m1-zr-1,:] has zc trailing zero
NB. columns, then v[0:m1-zr-1] and L[0:m1-zr-1,0:n1-zc-1]
NB. are supplied to larfl instead of v[:] and L[:,:].
NB.
NB. Syntax:
NB.   vapp=. ios larfLs
NB. where
NB.   ios  - 6-vector of integers (iosY,iosT,iosL,ioYs,ioLv,ioLh)
NB.   vapp - verb to update shrinked submatrix L by larfl;
NB.          usage:
NB.            Aupd=. cios vapp A
NB.          where
NB.            A    - m×n-matrix
NB.            cios - k×2 matrix of cIOSs (see layout below)
NB.            Aupd - m×n-matrix, A with updated shrinked
NB.                   submatrix L
NB.            k    ≥ max(iosY,iosT,iosL,⌊max(ioYs,ioLv,ioLh)/2⌋)+1
NB.            m    ≥ 0
NB.            n    ≥ 0
NB.
NB. Storage layout for cios:
NB.   iosY{cios    - cIOS of vector y (see larfl)
NB.   iosT{cios    - cIOS of scalar τ (see larfl)
NB.   iosL{cios    - cIOS of submatrix L in matrix A to
NB.                  update by larfl
NB.   ioYs({,)cios - cIOS half where size of vector y is
NB.                  stored, being (0{cIOS) if y is stored
NB.                  vertically, or (1{cIOS)) if horizontally
NB.   ioLv({,)cios - cIOS half (0{cIOS) of L where vertical
NB.                  axis is stored
NB.   ioLh({,)cios - cIOS half (1{cIOS) of L where horizontal
NB.                  axis is stored
NB.
NB. Application:
NB. - let cios=. (ciosY , ciosZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y, z (see larfl), scalar τ,
NB.   and submatrices R, L, respectively; vector y is stored
NB.   vertically; then to apply larfl to shrinked submatrix
NB.   L:
NB.     Aupd=. cios (0 2 4 0 8 9 larfLs) A
NB.
NB. References:
NB. [1] James W. Demmel, Mark Hoemmen, Yozo Hida, and E.
NB.     Jason Riedy. (2008) Non-Negative Diagonals and High
NB.     Performance on Low-Profile Matrices from Householder
NB.     QR. UCB/EECS-2008-76, May 30, 2008
NB.     LAPACK Working Note 203.
NB.     http://www.netlib.org/lapack/lawns/downloads/
NB.
NB. TODO:
NB. - CHECKME! after larf[gp] monad -> dyad

larfLs=: 1 : '(((ct0 cfuj (0{m)) unxu (3 4{m)) ((ct0c cfuj (2{m)) unxu (5{m)) ]) ((,`larfl) upd3ci (0 1 2{m)) ]'

NB. ---------------------------------------------------------
NB. larfRs
NB. Template adv. to make verbs to update shrinked submatrix
NB. by larfr. If n1-vector v has zc trailing zeros, and
NB. m1×(n1-zc)-submatrix R[:,0:n1-zc-1] has zr trailing zero
NB. rows, then v[0:n1-zc-1] and R[0:m1-zr-1,0:n1-zc-1]
NB. are supplied to larfr instead of v[:] and R[:,:].
NB.
NB. Syntax:
NB.   vapp=. ios larfRs
NB. where
NB.   ios  - 6-vector of integers (iosY,iosT,iosR,ioYs,ioRv,ioRh)
NB.   vapp - verb to update shrinked submatrix R by larfr;
NB.          usage:
NB.            Aupd=. cios vapp A
NB.          where
NB.            A    - m×n-matrix
NB.            cios - k×2 matrix of cIOSs (see layout below)
NB.            Aupd - m×n-matrix, A with updated shrinked
NB.                   submatrix R
NB.            k    ≥ max(iosY,iosT,iosR,⌊max(ioYs,ioRv,ioRh)/2⌋)+1
NB.            m    ≥ 0
NB.            n    ≥ 0
NB.
NB. Storage layout for cios:
NB.   iosY{cios    - cIOS of vector y (see larfr)
NB.   iosT{cios    - cIOS of scalar τ (see larfr)
NB.   iosR{cios    - cIOS of submatrix R in matrix A to
NB.                  update by larfr
NB.   ioYs({,)cios - cIOS half where size of vector y is
NB.                  stored, being (0{cIOS) if y is stored
NB.                  vertically, or (1{cIOS)) if horizontally
NB.   ioRv({,)cios - cIOS half (0{cIOS) of R where vertical
NB.                  axis is stored
NB.   ioRh({,)cios - cIOS half (1{cIOS) of R where horizontal
NB.                  axis is stored
NB.
NB. Application:
NB. - let cios=. (ciosY , ciosZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y, z (see larfr), scalar τ,
NB.   and submatrices R, L, respectively; vector y is stored
NB.   vertically; then to apply larfr to shrinked submatrix
NB.   R:
NB.     Aupd=. cios (0 2 3 0 7 6 larfRs) A
NB.
NB. References:
NB. [1] James W. Demmel, Mark Hoemmen, Yozo Hida, and E.
NB.     Jason Riedy. (2008) Non-Negative Diagonals and High
NB.     Performance on Low-Profile Matrices from Householder
NB.     QR. UCB/EECS-2008-76, May 30, 2008
NB.     LAPACK Working Note 203.
NB.     http://www.netlib.org/lapack/lawns/downloads/
NB.
NB. TODO:
NB. - CHECKME! after larf[gp] monad -> dyad

larfRs=: 1 : '(((ct0 cfuj (0{m)) unxu (3 4{m)) ((ct0r cfuj (2{m)) unxu (5{m)) ]) ((,`larfr) upd3ci (0 1 2{m)) ]'

NB. ---------------------------------------------------------
NB. larfRLs
NB. Template adv. to make verbs to update one shrinked
NB. submatrix by larfr, then another one by larfl. If
NB. n1-vector v has z trailing zeros, and m1×(n1-z)-submatrix
NB. R[:,0:n1-z-1] has zr trailing zero rows, then v[0:n1-z-1]
NB. and R[0:m1-zr-1,0:n1-z-1] are supplied to larfr instead
NB. of v[:] and R[:,:]. Then, if (n1-z)×k1-submatrix
NB. L[0:n1-z-1,:] has zc trailing zero columns, then
NB. v[0:n1-z-1] and L[0:n1-z-1,0:k1-zc-1] are supplied to
NB. larfl instead of v[:] and L[:,:].
NB.
NB. Syntax:
NB.   vapp=. ios larfRLs
NB. where
NB.   ios  - 9-vector of integers (iosY,iosT,iosR,iosL,ioYs,ioRv,ioRh,ioLv,ioLh)
NB.   vapp - verb to update shrinked submatrix R by larfr,
NB.          then shrinked submatrix L by larfl; usage:
NB.            Aupd=. cios vapp A
NB.          where
NB.            A    - m×n-matrix
NB.            cios - k×2 matrix of cIOSs (see layout below)
NB.            Aupd - m×n-matrix, A with updated shrinked
NB.                   submatrices R and L
NB.            k    ≥ max(iosY,iosT,iosR,iosL,⌊max(ioYs,ioRv,ioRh,ioLv,ioLh)/2⌋)+1
NB.            m    ≥ 0
NB.            n    ≥ 0
NB.
NB. Storage layout for cios:
NB.   iosY{cios    - cIOS of vector y (see larfr or larfl)
NB.   iosT{cios    - cIOS of scalar τ (see larfr or larfl)
NB.   iosR{cios    - cIOS of submatrix R in matrix A to
NB.                  update by larfr
NB.   iosL{cios    - cIOS of submatrix L in matrix A to
NB.                  update by larfl
NB.   ioYs({,)cios - cIOS half where size of vector y is
NB.                  stored, being (0{cIOS) if y is stored
NB.                  vertically, or (1{cIOS)) if horizontally
NB.   ioRv({,)cios - cIOS half (0{cIOS) of R where vertical
NB.                  axis is stored
NB.   ioRh({,)cios - cIOS half (1{cIOS) of R where horizontal
NB.                  axis is stored
NB.   ioLv({,)cios - cIOS half (0{cIOS) of L where vertical
NB.                  axis is stored
NB.   ioLh({,)cios - cIOS half (1{cIOS) of L where horizontal
NB.                  axis is stored
NB.
NB. Application:
NB. - let cios=. (ciosY , ciosZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y, z (see larfr), scalar τ,
NB.   and submatrices R, L, respectively; vector y is stored
NB.   vertically; then to apply larfr to shrinked submatrix
NB.   R, then larfl to shrinked submatrix L:
NB.     Aupd=. cios (0 2 3 4 0 6 7 8 9 larfRLs) A
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
NB. - aviods scanning V twice
NB.
NB. TODO:
NB. - CHECKME! after larf[gp] monad -> dyad

larfRLs=: 1 : '(((ct0 cfuj (0{m)) unxu (4 6{m)) ((ct0r cfuj (2{m)) unxu (5{m)) ]) (((ct0c cfuj (3{m)) unxu (8{m)) ((,`larfl) upd3ci (0 1 3{m)) ((,`larfr) upd3ci (0 1 2{m))) ]'

NB. ---------------------------------------------------------
NB. gerf0
NB. Template conj. to make verbs to generate and
NB. conditionally apply an elementary reflector to a matrix
NB.
NB. Syntax:
NB.   vcondapp=. (vref`vapp) gerf0 ios
NB. where
NB.   ios      - 2-vector of integers (iosYZ,iosT),
NB.              indices in cIOS bundle (see layout below)
NB.   vref     - verb to reflect vector y (see larfg or
NB.              larfp)
NB.   vapp     - verb to apply an elementary reflector to
NB.              submatrix (cIOSs are taken from x, matrix is
NB.              in y) from either the left or the right; is
NB.              called when τ≠0
NB.   vcondapp - verb to do steps:
NB.                1) reflect vector ((0{x){y) by verb vref
NB.                   and store result into vector ((1{x){y)
NB.                2) if ((2{x){y)≠0 then call verb (x vapp y)
NB.
NB. Storage layout for cios:
NB.   iosYZ{cios - cIOS of vector y or z (see larfr or larfl)
NB.   iosT{cios  - cIOS of scalar τ (see larfr or larfl)
NB.
NB. Application:
NB. - let cios=. (ciosYZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y (or z), scalar τ, and
NB.   submatrices R and L, respectively; vector y is stored
NB.   vertically, α is in head, τ is in tail
NB. - to reflect in the matrix A the vector y by larfg, then
NB.   to apply vector v from the left (emulate main loop of
NB.   xGEQR2 in LAPACK version 3.1.1):
NB.     vref=: 0 _1 & larfg                    NB. verb to reflect Y to Z
NB.     vapp=: 0 3 larfL                       NB. verb to apply v from the left to L
NB.     Aupd=. cios ((vref`vapp) gerf0 0 1) A  NB. apply gerund and ios to do the job
NB. - to do the same as above, but apply v from the right,
NB.   then from the left (emulate main loop of xGEHD2 in
NB.   LAPACK version 3.1.1):
NB.     vref=: 0 _1 & larfg                    NB. verb to reflect Y to Z
NB.     vapp=: 0 2 3 larfRL                    NB. verb to apply v from the right to R, then from the left to L
NB.     Aupd=. cios ((vref`vapp) gerf0 0 1) A  NB. apply gerund and ios to do the job
NB. - do the same as above, but with matrix shrinking:
NB.     vref=: 0 _1 & larfg                    NB. verb to reflect Y to Z
NB.     vapp=: 0 1 2 3 0 4 5 6 7 larfRLs       NB. verb to apply shrinked v from the right to shrinked R, then from the left to shrinked L
NB.     Aupd=. cios ((vref`vapp) gerf0 0 1) A  NB. apply gerund and ios to do the job

gerf0=: 2 : '[ (((1{u)`:6) ^: (0 ~: (0 ({,) ((1{n) fromci)))) (((0{u)`:6) updci (0{n))'

NB. ---------------------------------------------------------
NB. gerf02
NB. Template conj. to make verbs to generate and
NB. conditionally apply an elementary reflector to a matrix
NB.
NB. Syntax:
NB.   vcondapp=. (vref`vneg`vapp) gerf02 ios
NB. where
NB.   ios      - 2-vector of integers (iosYZ,iosT),
NB.              indices in cIOS bundle (see layout below)
NB.   vref     - verb to reflect vector y (see larfg or
NB.              larfp)
NB.   vneg     - verb to negate either a 1st row (left side
NB.              case), or a 1st column (right side case); is
NB.              called when τ=2; cIOSs are taken from x,
NB.              matrix is in y
NB.   vapp     - verb to apply an elementary reflector to a
NB.              submatrix from either the left or the right;
NB.              is called when τ≠0 and τ≠2; cIOSs are taken
NB.              from x, matrix is in y
NB.   vcondapp - verb to do steps:
NB.                1) reflect vector ((0{x){y) by verb vref
NB.                   and store result into vector ((1{x){y)
NB.                2a) if ((2{x){y)=0 then do nothing
NB.                2b) if ((2{x){y)=2 then call verb (x vneg y)
NB.                2c) othewise call verb (x vapp y)
NB.
NB. Storage layout for cios:
NB.   iosYZ{cios - cIOS of vector y or z (see larfr or larfl)
NB.   iosT{cios  - cIOS of scalar τ (see larfr or larfl)
NB.
NB. Application:
NB. - let cios=. (ciosYZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y, z, scalar τ, and
NB.   submatrices R and L, respectively; vector y is stored
NB.   vertically, α is in head, τ is in tail
NB. - to reflect in the matrix A the vector y by larfp, then
NB.   either to negate 1st row of submatrix L (if τ=2), or to
NB.   apply shrinked vector v from the left (if τ≠0 and τ≠2)
NB.   to the shrinked submatrix L (emulate main loop of
NB.   xGEQR2 in LAPACK version 3.2):
NB.     vref=: 0 _1 & larfp                          NB. verb to reflect Y to Z
NB.     vneg=: 3 nfv 0                               NB. verb to negate L's 1st row
NB.     vapp=: 0 1 3 0 6 7 larfLs                    NB. verb to apply v from the left to L
NB.     Aupd=. cios ((vref`vneg`vapp) gerf02 0 1) A  NB. apply gerund and ios to do the job
NB. - to reflect in the matrix A the vector y by larfp, then
NB.   either to negate 1st column of submatrix R and then to
NB.   negate 1st row of submatrix L (if τ=2), or to apply
NB.   shrinked vector v to the shrinked submatrix R from the
NB.   right, then to the shrinked submatrix L from the left
NB.   (if τ≠0 and τ≠2) (emulate main loop of xGEHD2 in
NB.   LAPACK version 3.2):
NB.     vref=: 0 _1 & larfp                          NB. verb to reflect Y to Z
NB.     vneg=: [ (3 nfv 0) (2 nfv 1)                 NB. verb to negate R's 1st col, then L's 1st row
NB.     vapp=: 0 1 2 3 0 4 5 6 7 larfRLs             NB. verb to apply shrinked v from the right to shrinked R, then from the left to shrinked L
NB.     Aupd=. cios ((vref`vneg`vapp) gerf02 0 1) A  NB. apply gerund and ios to do the job
NB. - to do the same as previous, but without shrinking
NB.   (optimize speed for non-band matrices):
NB.     vref=: 0 _1 & larfp                          NB. verb to reflect Y to Z
NB.     vneg=: [ (3 nfv 0) (2 nfv 1)                 NB. verb to negate R's 1st col, then L's 1st row
NB.     vapp=: 0 2 3 larfRL                          NB. verb to apply v from the right to R, then from the left to L
NB.     Aupd=. cios ((vref`vneg`vapp) gerf02 0 1) A  NB. apply gerund and ios to do the job

gerf02=: 2 : '[ ((]`(}.u)) @. (0 2 i. (0 ({,) ((1{n) fromci)))) (((0{u)`:6) updci (0{n))'

NB. ---------------------------------------------------------
NB. larft
NB. Form the triangular factor T of a complex block reflector
NB. H of order n, which is defined as a product of k
NB. elementary reflectors. If direction is forward, then:
NB.   H = H(1) H(2) ... H(k) and T is upper triangular
NB. If direction is backward, then:
NB.   H = H(k) ... H(2) H(1) and T is lower triangular
NB. If layout is columnwise, then:
NB.   H = I - V * T * V'
NB. If layout is rowwise, then:
NB.   H = I - V' * T * V
NB.
NB. Syntax:
NB.   T=. (iociosV)larft V
NB. where
NB.
NB. TODO: scan trailing zeros

LARFTFCDCIOS=: 5 2 $ 1j_1 0j1 1j_1 1 0 1 0j1 0j1 0j1 1  NB. A,v,τ,Ttru0,Tcol

NB. cios1=. mkcios1larftfc (t,h,l,w)
mkcios1larftfc=: 3 : 0
  th1=. >: (+/) }: 't1 h1 l1'=. 1 _1 1 + }: 't h l w'=. y
  5 2 $ (t1,l,t1,l1,(t+h),l1,th1,0,th1,1) j. (h1,1,h1,1 1 1 1 1 1 1)
)

NB. Tr * ((v' * A) * (+ τ))'
NB. larftfcstep=: (((((mp~ ((1 & (0:}))@ct)) dbg 'y''*x')`((ct@(* (+@(0&({,))))) dbg '(*+t)''')`((mp~) dbg 'mp~') map4ci 0 1 2 3 4) dbg 'map4')) step
larftfcstep=: ((mp~ ((1 & (0:}))@ct))`(ct@(* (+@(0&({,)))))`(mp~) map4ci 0 1 2 3 4) step

NB. VT=. ciosV larftfc A
larftfc=: 4 : 0
  'm n'=. $ y
  'th th1'=. 0 1 + (+/) 2 {. 't h l w'=. thlw=. , +. x
  liostau=. (l+n*th) hs2ios w         NB. current block's τ-s vector's lIOS
  liosdiag=. hds2ios (n*th1),(>:n),w  NB. diagonal's lIOS of T=. (_1 1 * k) {. A
  tau=.liostau ({,) y                 NB. FIXME! 
  y=. tau ((liosdiag"_) }) y          NB. FIXME! copy τ-s vector to T's diagonal
  cios1=. mkcios1larftfc thlw
  0 {:: (<:w) (LARFTFCDCIOS & larftfcstep) (y ; cios1)
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB.*tlarfg v test larfg

tlarfg=: 3 : 0
)

NB. ---------------------------------------------------------
NB.*tlarfp v test larfp

tlarfp=: 3 : 0
)

NB. ---------------------------------------------------------
NB.*tref v test reflectors with data suppplied

tref=: 3 : 0
)

NB. ---------------------------------------------------------
NB.*testref v test reflectors

testref=: 3 : 0
)
