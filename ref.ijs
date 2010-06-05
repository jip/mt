NB. ref.ijs
NB. Reflections
NB.
NB. larfg    Generate an elementary reflector
NB. larfp    Generate an elementary reflector
NB.
NB. larfl    Apply an elementary reflector from the left
NB. larfr    Apply an elementary reflector from the right
NB.
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
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. - extract n-th cIOS cios from x, then extract by cfrom
NB.   from y submatrix defined by cios
NB. - apply u
NB. - make imaginary
cfuj=: 2 : 'j. @ u @ (n cfrom)'

NB. update n-th from x (warning!: linear IOS) by (x u y)
unxu=: 2 : '((n({,)[) - u)(n"_)}['

NB. ---------------------------------------------------------
NB. larf
NB. Template conj. to make verbs to apply an elementary
NB. reflector H = I - τ*v*v' to a matrix A from either the
NB. left or the right
NB.
NB. Syntax:
NB.   vapp=. vmul larf vtau
NB. where
NB.   vmul - verb to calculate (v*v'*A) or (A*v*v'); is
NB.          called as: (A vmul (β,v,τ))
NB.   vtau - verb to pre-process τ; is called as: (vtau τ)
NB.   vapp - verb to apply an elementary reflector; is called
NB.          as: ((β,v,τ) vapp A)
NB.   β    - any scalar, is not used
NB.   v    - (n-1)-vector of v[1:n-1]; v[0]=1 is not stored
NB.   τ    - scalar
NB.
NB. Notes:
NB. - pre-ravel (β,v,τ)
NB. - no check τ≠0
NB. - assign v[0]=1
NB. - pre-process τ by verb v

larf=: 2 : '(([ - (v @ (_1 { ])) * ((u (1 & (0 }) @ }:)))) ,)~'

NB. =========================================================
NB. Interface

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

NB. ---------------------------------------------------------
NB. larfl
NB. Apply an elementary reflector H' to a matrix A from left:
NB.   H'*A = A - conj(tau) * v * (v' * A)
NB.
NB. Syntax:
NB.   HA=. z larfl A
NB. where
NB.   A  - n×n-matrix
NB.   z  - (n+1)-vector (β,v[1:n-1],τ), where β is any
NB.        scalar
NB.   HA - n×n-matrix H'*A
NB.
NB. Notes:
NB. - differs from LAPACK's xLARF('L') in:
NB.   - τ conjugation
NB.   - not scans trailing zeros in v

larfl=: (] */ (mp~ +)) larf +

NB. ---------------------------------------------------------
NB. larfr
NB. Apply an elementary reflector H to a matrix A from right:
NB.   A*H = A - tau * (A * v) * v'
NB.
NB. Syntax:
NB.   AH=. z larfr A
NB. where
NB.   A  - n×n-matrix
NB.   z  - (n+1)-vector (β,v[1:n-1],τ), where β is any
NB.        scalar
NB.   AH - n×n-matrix A*H
NB.
NB. Notes:
NB. - differs from LAPACK's xLARF('R') in:
NB.   - not scans trailing zeros in v

larfr=: (mp */ (+@])) larf ]

NB. ---------------------------------------------------------
NB. larfL
NB. Template adv. to make verbs to update submatrix by larfl
NB.
NB. Syntax:
NB.   vapp=. ios larfL
NB. where
NB.   ios  - 2-vector of integers (iosZ,iosL)
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
NB. - let cios=. (ciosY , ciosZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y, z (see larfl), scalar τ,
NB.   and submatrices R, L, respectively; then to apply larfl
NB.   to submatrix L:
NB.     Aupd=. cios (1 4 larfL) A

larfL=: 1 : 'larfl cupd2 m'

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

larfR=: 1 : 'larfr cupd2 m'

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

larfRL=: 1 : '[ ((0 2{m) larfL) ((0 1{m) larfR)'

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

larfLs=: 1 : '(((ct0 cfuj (0{m)) unxu (3 4{m)) ((ct0c cfuj (2{m)) unxu (5{m)) ]) ((,`larfl) cupd3 (0 1 2{m)) ]'

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

larfRs=: 1 : '(((ct0 cfuj (0{m)) unxu (3 4{m)) ((ct0r cfuj (2{m)) unxu (5{m)) ]) ((,`larfr) cupd3 (0 1 2{m)) ]'

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

larfRLs=: 1 : '(((ct0 cfuj (0{m)) unxu (4 6{m)) ((ct0r cfuj (2{m)) unxu (5{m)) ]) (((ct0c cfuj (3{m)) unxu (8{m)) ((,`larfl) cupd3 (0 1 3{m)) ((,`larfr) cupd3 (0 1 2{m))) ]'

NB. ---------------------------------------------------------
NB. gerf0
NB. Template conj. to make verbs to generate and
NB. conditionally apply an elementary reflector to a matrix
NB.
NB. Syntax:
NB.   vcondapp=. (vref`vapp) gerf0 ios
NB. where
NB.   ios      - 3-vector of integers (iosY,iosZ,iosT),
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
NB.   iosY{cios - cIOS of vector y (see larfr or larfl)
NB.   iosZ{cios - cIOS of vector z (see larfr or larfl)
NB.   iosT{cios - cIOS of scalar τ (see larfr or larfl)
NB.
NB. Application:
NB. - let cios=. (ciosY , ciosZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y, z, scalar τ, and
NB.   submatrices R and L, respectively; vector y is stored
NB.   vertically
NB. - to reflect in the matrix A the vector y by larfg, then
NB.   to apply vector v from the left (emulate main loop of
NB.   xGEQR2 in LAPACK version 3.1.1):
NB.     vapp=: 1 4 larfL                          NB. verb to apply v from the left to L
NB.     Aupd=. cios ((larfg`vapp) gerf0 0 1 2) A  NB. apply gerund and ios to do the job
NB. - to do the same as above, but apply v from the right,
NB.   then from the left (emulate main loop of xGEHD2 in
NB.   LAPACK version 3.1.1):
NB.     vapp=: 1 3 4 larfRL                       NB. verb to apply v from the right to R, then from the left to L
NB.     Aupd=. cios ((larfg`vapp) gerf0 0 1 2) A  NB. apply gerund and ios to do the job
NB. - do the same as above, but with matrix shrinking:
NB.     vapp=: 0 2 3 4 0 6 7 8 9 larfRLs          NB. verb to apply shrinked v from the right to shrinked R, then from the left to shrinked L
NB.     Aupd=. cios ((larfg`vapp) gerf0 0 1 2) A  NB. apply gerund and ios to do the job

gerf0=: 2 : '[ (((1{u)`:6) ^: (0 ~: (0 ({,) ((2{n) cfrom)))) (((0{u)`:6) cmap (0 1{n))'

NB. ---------------------------------------------------------
NB. gerf02
NB. Template conj. to make verbs to generate and
NB. conditionally apply an elementary reflector to a matrix
NB.
NB. Syntax:
NB.   vcondapp=. (vref`vneg`vapp) gerf02 ios
NB. where
NB.   ios      - 3-vector of integers (iosY,iosZ,iosT),
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
NB.   iosY{cios - cIOS of vector y (see larfr or larfl)
NB.   iosZ{cios - cIOS of vector z (see larfr or larfl)
NB.   iosT{cios - cIOS of scalar τ (see larfr or larfl)
NB.
NB. Application:
NB. - let cios=. (ciosY , ciosZ , ciosT , ciosR ,: ciosL)
NB.   is cIOS bundle for vectors y, z, scalar τ, and
NB.   submatrices R and L, respectively; vector y is stored
NB.   vertically
NB. - to reflect in the matrix A the vector y by larfp, then
NB.   either to negate 1st row of submatrix L (if τ=2), or to
NB.   apply shrinked vector v from the left (if τ≠0 and τ≠2)
NB.   to the shrinked submatrix L (emulate main loop of
NB.   xGEQR2 in LAPACK version 3.2):
NB.     vneg=: 4 nfv 0                                  NB. verb to negate L's 1st row
NB.     vapp=: 0 2 4 0 8 9 larfLs                       NB. verb to apply v from the left to L
NB.     Aupd=. cios ((larfp`vneg`vapp) gerf02 0 1 2) A  NB. apply gerund and ios to do the job
NB. - to reflect in the matrix A the vector y by larfp, then
NB.   either to negate 1st column of submatrix R and then to
NB.   negate 1st row of submatrix L (if τ=2), or to apply
NB.   shrinked vector v to the shrinked submatrix R from the
NB.   right, then to the shrinked submatrix L from the left
NB.   (if τ≠0 and τ≠2) (emulate main loop of xGEHD2 in
NB.   LAPACK version 3.2):
NB.     vneg=: [ (4 nfv 0) (3 nfv 1)                    NB. verb to negate R's 1st col, then L's 1st row
NB.     vapp=: 0 2 3 4 0 6 7 8 9 larfRLs                NB. verb to apply shrinked v from the right to shrinked R, then from the left to shrinked L
NB.     Aupd=. cios ((larfp`vneg`vapp) gerf02 0 1 2) A  NB. apply gerund and ios to do the job
NB. - to do the same as previous, but without shrinking
NB.   (optimize speed for non-band matrices):
NB.     vneg=: [ (4 nfv 0) (3 nfv 1)                    NB. verb to negate R's 1st col, then L's 1st row
NB.     vapp=: 1 3 4 larfRL                             NB. verb to apply v from the right to R, then from the left to L
NB.     Aupd=. cios ((larfp`vneg`vapp) gerf02 0 1 2) A  NB. apply gerund and ios to do the job

gerf02=: 2 : '[ ((]`(}.u)) @. (0 2 i. (0 ({,) ((2{n) cfrom)))) (((0{u)`:6) cmap (0 1{n))'

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
