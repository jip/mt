NB. Reflection
NB.
NB. larfx      Generate an elementary reflector
NB. larfxxx    Generate an elementary reflector
NB. larxtxx    Form the triangular factor of a block
NB.            reflector
NB. larxxxxx   Apply an elementary reflector or its transpose
NB.            to a matrix from either the left or the right
NB. larxbxxxx  Build and apply a block reflector or its
NB.            transpose to a matrix from either the left or
NB.            the right
NB. refga      Conj. to make monad to generate and apply
NB.            reflector
NB.
NB. testlarfg  Test larfx by general vector
NB. testlarf   Test larfxxxx by general matrix
NB. testlarz   Test larzxxxx by general matrix
NB. testlarft  Test larftxx by general matrix
NB. testlarzt  Test larztxx by general matrix
NB. testlarfb  Test larfbxxxx by general matrix
NB. testlarzb  Test larzbxxxx by general matrix
NB. testref    Adv. to make verb to test larxxxxxx by matrix
NB.            of generator and shape given
NB.
NB. Copyright 2010,2011,2013,2017,2018,2020,2021,2023,2024,
NB.           2025 Igor Zhuravlov
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

NB. =========================================================
NB. Configuration

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Constants

NB. Safe minimum
REFSAFMIN=: FP_SFMIN % FP_EPS

NB. ---------------------------------------------------------
NB. larxtx
NB.
NB. Description:
NB.   Conj. to make monad to form the triangular factor Τ of
NB.   a block reflector
NB.
NB. Syntax:
NB.   T=. (ioTau larxtx mprod) VTau
NB. where
NB.   ioTau - IOs Tau in VTau:
NB.             Tau -: ioTau { VTau
NB.   mprod - dyad to multiply matrices, is called as:
NB.             M3=. M1 mprod M2
NB.   VTau  - (m+1)×k- or k×(n+1)-matrix, V and Tau combined
NB.   T     - k×k-matrix, triangular
NB.   V     - m×k- or k×n-matrix, the unit trapezoidal
NB.   Tau   - k-vector, τ[0:k-1] corresp. to V

larxtf=: 2 : '(] (] ((0 ,~ [) ,.  (mp }:) , _1 { ]) ({~ (<@;~ i.@>:)@#       ))^:(-&#)  1  1&rt)@(-@(m&{) *"1 (_1;a:) tru@setdiag (v ct)@(0&(m})))'
larxtb=: 2 : '(] (] ((0 ,  [) ,.~ (mp }.) ,~ 0 { ]) ((<@;~ <@i.)@<:@(-&#) { [))^:(-&#) _1 _1&rt)@(-@(m&{) *"1 (_1;a:) trl@setdiag (v ct)@(0&(m})))'

NB. ---------------------------------------------------------
NB. larxxxxx
NB.
NB. Description:
NB.   Adv. to make dyad to apply an elementary reflector H
NB.   or its transpose H' to a matrix, from either the left
NB.   or the right. H is defined by pair (v,τ).
NB.
NB. Syntax:
NB.   eCupd=. Vtau (iotau larxxxxx) eC
NB. where
NB.   iotau - IO tau in Vtau:
NB.             tau -: iotau { Vtau
NB.   Vtau  - vector V augmented by scalar τ
NB.   eC    - matrix C augmented by trash vector
NB.   eCupd - an updated eC

larxlcxc=: 1 : '] - [ */ (mp~ +@(0&(m}) * m&{))~'  NB. C - v * ((v * τ)' * C)
larxlcxr=: 1 : '] - +@(* m&{)@[ */ (0 m} [) mp ]'  NB. C - (τ * v)' * (v * C)

larxlnxc=: 1 : '] - [ */ (mp~ m&{ * +@(0&(m})))~'  NB. C - v * ((τ * v') * C)
larxlnxr=: 1 : '] - (+ * m&{)@[ */ (mp~ 0&(m}))~'  NB. C - (v' * τ) * (v * C)

larxrcxc=: 1 : '] - (mp 0&(m}))~ */ +@(* m&{)@['   NB. C - (C * v) * (v * τ)'
larxrcxr=: 1 : '] - (mp +@(m&{ * 0&(m})))~ */ ['   NB. C - (C * (τ * v)') * v

larxrnxc=: 1 : '] - (mp 0&(m}) * m&{)~ */ +@['     NB. C - (C * (v * τ)) * v'
larxrnxr=: 1 : '] - (mp +@(0&(m})) * m&{)~ */ ['   NB. C - (C * (v' * τ)) * v

NB. ---------------------------------------------------------
NB. larxbxxxx
NB.
NB. Description:
NB.   Conj. to make dyad to apply an elementary reflector H
NB.   or its transpose H' to a matrix, from either the left
NB.   or the right. H is defined by pair (V,Τ).
NB.
NB. Syntax:
NB.   eCupd=. VTau (ioTau larxbxxxx makeT) eC
NB. where
NB.   ioTau - IO tau in VTau:
NB.             Tau -: ioTau { VTau
NB.   makeT - monad to make Τ, usually one of larxtxx, is
NB.           called as:
NB.             T=. makeT VTau
NB.   VTau  - matrix V augmented by vector τ
NB.   eC    - matrix C augmented by trash vector
NB.   eCupd - an updated eC

larxblcxc=: 2 : '] - [ mp (mp~ ct@(0&(m}) mp v))~'   NB. C - V * ((V * Τ)' * C)
larxblcxr=: 2 : '] - ct@(mp~ v)@[ mp (0 m} [) mp ]'  NB. C - (Τ * V)' * (V * C)

larxblnxc=: 2 : '] - [ mp (mp~ v mp ct@(0&(m})))~'   NB. C - V * ((Τ * V') * C)
larxblnxr=: 2 : '] - (ct mp v)@[ mp (mp~ 0&(m}))~'   NB. C - (V' * Τ) * (V * C)

larxbrcxc=: 2 : '] - (mp 0&(m}))~ mp ct@(mp v)@['    NB. C - (C * V) * (V * Τ)'
larxbrcxr=: 2 : '] - (mp ct@(v mp 0&(m})))~ mp ['    NB. C - (C * (Τ * V)') * V

larxbrnxc=: 2 : '] - (mp 0&(m}) mp v)~ mp ct@['      NB. C - (C * (V * Τ)) * V'
larxbrnxr=: 2 : '] - (mp ct@(0&(m})) mp v)~ mp ['    NB. C - (C * (V' * Τ)) * V

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. larfg
NB. larfp
NB.
NB. Description:
NB.   Generate an elementary reflector H of order n such that
NB.     H^H * (α,x) = (β,0)
NB. where
NB.   H - represented in factored form by n-vector (1,v) and
NB.       scalar τ
NB.
NB. Syntax:
NB.   z=. iso larfg y
NB.   z=. iso larfp y
NB. where
NB.   iso - 2-vector of integers (ioa,iot)
NB.   ioa - lIO α in y
NB.   iot - lIO pre-allocated scalar in y
NB.   y   - (n+1)-vector having scalar α ∈ ℂ at index ioa,
NB.         any scalar at index iot, and vector x ∈ ℂ^(n-1)
NB.         in the rest elements, vector to reflect is:
NB.           (<<< iot) { y
NB.   z   - (n+1)-vector having scalar β ∈ ℝ (larfp [1]
NB.         provides β≥0) at index ioa, scalar τ ∈ ℂ at index
NB.         iot, and vector v ∈ ℂ^(n-1) in the rest elements,
NB.         reflected vector is:
NB.           beta ioa} n $ 0
NB.   n   > 0, the length of α and x combined
NB.
NB. Formula:
NB.   H = I - (1,v) * τ * (1,v)^H
NB.
NB. Application:
NB. - reflect vector (α,x) by larfg and store τ at tail:
NB.     z=. 0 _1 larfg (alpha , x , _.)
NB.     v=. (<<<0 _1) { z
NB.     'beta tau'=. 0 _1 { z
NB. - reflect vector (x,α) by larfp and store τ at head:
NB.     z=. _1 0 larfp (_. , x , alpha)
NB.     v=. (<<<_1 0) { z
NB.     'beta tau'=. _1 0 { z
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.    (n {. beta) -: H (mp~ ct)~ }: y  NB. H^H * (α,x) = (β,0)
NB.    I           -:   (mp~ ct)~ H     NB. H^H * H = I
NB. where
NB.   x     - (n-1)-vector
NB.   alpha - scalar
NB.   y=. alpha , x , _.
NB.   z=. 0 _1 larfg y
NB.   v=. (<<<0 _1) { z
NB.   'beta tau'=. 0 _1 { z
NB.   I=. idmat n
NB.   H=. I - tau (] */ (* +)) (1,v)
NB.
NB. Notes:
NB. - IEEE floating point configuration is encoded implicitly
NB. - larfg models LAPACK's xLARFG
NB. - larfp models LAPACK's xLARFP
NB. - larfp provides τ=2 ↔ ||x|| ≈ 0
NB. - not zeroing v in larfp in case τ=2 relies on fact:
NB.   'comparison tolerance tol>0'; otherwise (tol=0) x
NB.   would be filled by zeros
NB.
NB. References:
NB. [1] James W. Demmel, Mark Hoemmen, Yozo Hida, E. Jason
NB.     Riedy. Non-Negative Diagonals and High Performance on
NB.     Low-Profile Matrices from Householder QR.
NB.     UCB/EECS-2008-76, May 30, 2008.
NB.     LAPACK Working Note 203
NB.     http://www.netlib.org/lapack/lawns/downloads/

larfg=: 4 : 0
  'ioa iot'=. x
  alpha=. ioa { y
  y=. 0 iot} y                            NB. τ := 0
  ynorm=. norms y
  if. ynorm (=!.0) | 9 o. alpha do.       NB. ||y|| == ||(α,x,0)|| == ||α|| and α ∈ ℝ ?
    y                                     NB. (α,0,0) i.e. H==I, τ==0, β==α, v==0
  else.
    if. REFSAFMIN > ynorm do.             NB. xnorm, β may be inaccurate; scale x and recompute them
      y=. y % REFSAFMIN                   NB. (α_scaled,x_scaled,0)
      beta=. (9 o. alpha) negpos norms y  NB. use Re(α) instead Re(α_ascaled) since sign(Re(α)) == sign(Re(α_scaled)); |β_scaled| ∈ [REFSAFMIN,1)
      dzeta=. beta - ioa { y              NB. ζ := β_scaled-α_scaled
      tau=. dzeta % beta                  NB. τ := ζ/β_scaled
      beta=. REFSAFMIN * beta             NB. unscale β; if α is subnormal, it may lose relative accuracy
    else.
      beta=. (9 o. alpha) negpos ynorm    NB. β := -copysign(||y||,Re(α)), since ||y|| ≥ 0
      dzeta=. beta - alpha
      tau=. dzeta % beta
    end.
    y=. y % - dzeta                       NB. z := (trash,v,0)
    y=. (beta , tau) x} y                 NB. z := (β_scaled,v,τ)
  end.
)

larfp=: 4 : 0
  'ioa iot'=. x
  alpha=. ioa { y
  xnorm=. norms 0 x} y                                    NB. ||x||
  if. 0 = xnorm do.
    y=. ((| , 1 - *) alpha) x} y                          NB. replace α by |β| and τ by (1-α/|α|)
  else.
    beta=. (9 o. alpha) negpos norms alpha , xnorm        NB. β := -copysign(||y||,Re(α))
    y=. (| beta) iot} y
    if. FP_SFMIN > | beta do.
      y=. y % FP_SFMIN                                    NB. scale (α,x[1],...,x[n-1],|β|)
      xnorm=. xnorm % FP_SFMIN
    end.
    if. 0 <: beta do.
      dzeta=. -/ x { y                                    NB. ζ := α_scaled-|β_scaled|
      tau=. - dzeta % iot { y                             NB. τ := -ζ/|β_scaled|
    else.
      beta=. - beta                                       NB. |β_unscaled|
      'realpha imalpha'=. +. ioa { y                      NB. Re(α_scaled) , Im(α_scaled)
      gamma=. realpha + iot { y                           NB. γ := Re(α_scaled)+|β_scaled|
      delta=. - (imalpha , xnorm) ([ +/@:*"1!.0 %) gamma  NB. δ := -(Im(α_scaled)*(Im(α_scaled)/γ)+||x||*(||x||/γ))
      dzeta=. delta j. imalpha                            NB. ζ := δ+i*Im(α_scaled)
      tau=. - dzeta % iot { y                             NB. τ := -ζ/|β_scaled|
    end.
    y=. y % dzeta
    y=. (beta , tau) x} y                                 NB. replace α_scaled by |β_unscaled| and |β_scaled| by τ
  end.
)

NB. ---------------------------------------------------------
NB. Monad      Input             Output                  β
NB. larfgf     (α x[1:n-1] 0)    (β v[1:n-1] τ)          ∈ ℝ
NB. larfgfc    (α x[1:n-1] 0)    (β v[1:n-1] conj(τ))    ∈ ℝ
NB. larfgb     (0 x[1:n-1] α)    (τ v[1:n-1] β)          ∈ ℝ
NB. larfgbc    (0 x[1:n-1] α)    (conj(τ) v[1:n-1] β)    ∈ ℝ
NB. larfpf     (α x[1:n-1] 0)    (β v[1:n-1] τ)          ≥ 0
NB. larfpfc    (α x[1:n-1] 0)    (β v[1:n-1] conj(τ))    ≥ 0
NB. larfpb     (0 x[1:n-1] α)    (τ v[1:n-1] β)          ≥ 0
NB. larfpbc    (0 x[1:n-1] α)    (conj(τ) v[1:n-1] β)    ≥ 0
NB.
NB. Description:
NB.   Generate an elementary reflector, see larfg, larfp for
NB.   details.

larfgf=:  0 _1&larfg
larfgb=: _1  0&larfg
larfpf=:  0 _1&larfp
larfpb=: _1  0&larfp

larfgfc=: +&.(_1&{)@larfgf
larfgbc=: +&.( 0&{)@larfgb
larfpfc=: +&.(_1&{)@larfpf
larfpbc=: +&.( 0&{)@larfpb

NB. ---------------------------------------------------------
NB. larftbc
NB. larztbc
NB.
NB. Description:
NB.   Form the triangular factor Τ of a block reflector H:
NB.     H = H(k-1) * ... * H(1) * H(0) = I - V * Τ * V'
NB. where
NB.   Τ - lower triangular
NB.
NB. Syntax:
NB.   T=. larxtbc VTau
NB. where
NB.   VTau - (m+1)×k-matrix (tau,V)
NB.   V    - m×k-matrix, the unit upper trapezoidal with 1s
NB.          on (k-m)-th diagonal and 0s below
NB.   tau  - k-vector τ[0:k-1] corresp. to V
NB.   T    - k×k-matrix, the lower triangular
NB.
NB. Notes:
NB. - larftbc models LAPACK's xLARFT('B','C')
NB. - larztbc models LAPACK's xLARZT('B','C')

larftbc=:          0  larxtb (mp~)
larztbc=:         _1  larxtb (mp~)

NB. ---------------------------------------------------------
NB. larftbr
NB. larztbr
NB.
NB. Description:
NB.   Form the triangular factor Τ of a block reflector H:
NB.     H = H(k-1) * ... * H(1) * H(0) = I - V' * Τ * V
NB. where
NB.   Τ - lower triangular
NB.
NB. Syntax:
NB.   T=. larxtbr VTau
NB. where
NB.   VTau - k×(n+1)-matrix (tau,.V)
NB.   V    - k×n-matrix, the unit lower trapezoidal with 1s
NB.          on (n-k)-th diagonal and 0s above
NB.   tau  - k-vector τ[0:k-1] corresp. to V
NB.   T    - k×k-matrix, the lower triangular
NB.
NB. Notes:
NB. - larftbr models LAPACK's xLARFT('B','R')
NB. - larztbr models LAPACK's xLARZT('B','R')

larftbr=: (< a: ;  0) larxtb mp
larztbr=: (< a: ; _1) larxtb mp

NB. ---------------------------------------------------------
NB. larftfc
NB. larztfc
NB.
NB. Description:
NB.   Form the triangular factor Τ of a block reflector H:
NB.     H = H(0) * H(1) * ... * H(k-1) = I - V * Τ * V'
NB. where
NB.   Τ - upper triangular
NB.
NB. Syntax:
NB.   T=. larxtfc VTau
NB. where
NB.   VTau - (m+1)×n-matrix (V,tau)
NB.   V    - m×n-matrix, the unit lower trapezoidal with 1s
NB.          on 0-th diagonal and 0s above
NB.   tau  - n-vector τ[0:n-1] corresp. to V
NB.   T    - n×n-matrix, the upper triangular
NB.
NB. Notes:
NB. - larftfc models LAPACK's xLARFT('F','C')
NB. - larztfc models LAPACK's xLARZT('F','C')

larftfc=:         _1  larxtf (mp~)
larztfc=:          0  larxtf (mp~)

NB. ---------------------------------------------------------
NB. larftfr
NB. larztfr
NB.
NB. Description:
NB.   Form the triangular factor Τ of a block reflector H:
NB.     H = H(0) * H(1) * ... * H(k-1) = I - V' * Τ * V
NB. where
NB.   Τ - upper triangular
NB.
NB. Syntax:
NB.   T=. larxtfr VTau
NB. where
NB.   VTau - k×(n+1)-matrix (V,.tau)
NB.   V    - k×n-matrix, the unit upper trapezoidal with 1s
NB.          on 0-th diagonal and 0s below
NB.   tau  - k-vector τ[0:k-1] corresp. to V
NB.   T    - k×k-matrix, the upper triangular
NB.
NB. Notes:
NB. - larftfr models LAPACK's xLARFT('F','R')
NB. - larztfr models LAPACK's xLARZT('F','R')

larftfr=: (< a: ; _1) larxtf mp
larztfr=: (< a: ;  0) larxtf mp

NB. ---------------------------------------------------------
NB. Dyad        Action     Side     Tran    Dir    Layout        eC
NB. larflcbc    H' * C     left     ct      bwd    columnwise    0, C
NB. larflcbr    H' * C     left     ct      bwd    rowwise       0, C
NB. larflcfc    H' * C     left     ct      fwd    columnwise    C, 0
NB. larflcfr    H' * C     left     ct      fwd    rowwise       C, 0
NB. larflnbc    H  * C     left     none    bwd    columnwise    0, C
NB. larflnbr    H  * C     left     none    bwd    rowwise       0, C
NB. larflnfc    H  * C     left     none    fwd    columnwise    C, 0
NB. larflnfr    H  * C     left     none    fwd    rowwise       C, 0
NB. larfrcbc    C  * H'    right    ct      bwd    columnwise    0,.C
NB. larfrcbr    C  * H'    right    ct      bwd    rowwise       0,.C
NB. larfrcfc    C  * H'    right    ct      fwd    columnwise    C,.0
NB. larfrcfr    C  * H'    right    ct      fwd    rowwise       C,.0
NB. larfrnbc    C  * H     right    none    bwd    columnwise    0,.C
NB. larfrnbr    C  * H     right    none    bwd    rowwise       0,.C
NB. larfrnfc    C  * H     right    none    fwd    columnwise    C,.0
NB. larfrnfr    C  * H     right    none    fwd    rowwise       C,.0
NB.
NB. Description:
NB.   Apply an elementary reflector H or its transpose H' to
NB.   a matrix, from either the left or the right. H is
NB.   defined by pair (v,τ).
NB.
NB. Syntax:
NB.   eCupd=. vtau larfxxxx eC
NB. where
NB.   eC    - matrix C augmented by trash vector
NB.   vtau  - vector v augmented by scalar τ
NB.   eCupd - an updated eC
NB.   v     - vector with 1 at head (forward direction) or
NB.           tail (backward direction)
NB.
NB. Notes:
NB. - models LAPACK's xLARF
NB. - larfxxxx and larfbxxxx are topological equivalents
NB. - larfxxxx and larzxxxx are equivalent up to τ position
NB. - if τ=0 then v can have any element values

larflcbc=:  0 larxlcxc
larflcbr=:  0 larxlcxr
larflcfc=: _1 larxlcxc
larflcfr=: _1 larxlcxr

larflnbc=:  0 larxlnxc
larflnbr=:  0 larxlnxr
larflnfc=: _1 larxlnxc
larflnfr=: _1 larxlnxr

larfrcbc=:  0 larxrcxc
larfrcbr=:  0 larxrcxr
larfrcfc=: _1 larxrcxc
larfrcfr=: _1 larxrcxr

larfrnbc=:  0 larxrnxc
larfrnbr=:  0 larxrnxr
larfrnfc=: _1 larxrnxc
larfrnfr=: _1 larxrnxr

NB. ---------------------------------------------------------
NB. Dyad        Action     Side     Tran    Dir    Layout        eC
NB. larzlcbc    H' * C     left     ct      bwd    columnwise    C, 0
NB. larzlcbr    H' * C     left     ct      bwd    rowwise       C, 0
NB. larzlcfc    H' * C     left     ct      fwd    columnwise    0, C
NB. larzlcfr    H' * C     left     ct      fwd    rowwise       0, C
NB. larzlnbc    H  * C     left     none    bwd    columnwise    C, 0
NB. larzlnbr    H  * C     left     none    bwd    rowwise       C, 0
NB. larzlnfc    H  * C     left     none    fwd    columnwise    0, C
NB. larzlnfr    H  * C     left     none    fwd    rowwise       0, C
NB. larzrcbc    C  * H'    right    ct      bwd    columnwise    C,.0
NB. larzrcbr    C  * H'    right    ct      bwd    rowwise       C,.0
NB. larzrcfc    C  * H'    right    ct      fwd    columnwise    0,.C
NB. larzrcfr    C  * H'    right    ct      fwd    rowwise       0,.C
NB. larzrnbc    C  * H     right    none    bwd    columnwise    C,.0
NB. larzrnbr    C  * H     right    none    bwd    rowwise       C,.0
NB. larzrnfc    C  * H     right    none    fwd    columnwise    0,.C
NB. larzrnfr    C  * H     right    none    fwd    rowwise       0,.C
NB.
NB. Description:
NB.   Apply an elementary reflector H or its transpose H' to
NB.   a matrix, from either the left or the right. H is
NB.   defined by pair (v,τ).
NB.
NB. Syntax:
NB.   eCupd=. vtau larzxxxx eC
NB. where
NB.   eC    - matrix C augmented by trash vector
NB.   vtau  - vector v augmented by scalar τ
NB.   eCupd - an updated eC
NB.   v     - vector with 1 at head (backward direction) or
NB.           tail (forward direction), and 0s in atoms to be
NB.           ignored
NB.
NB. Notes:
NB. - models LAPACK's xLARZ
NB. - larzxxxx and larzbxxxx are topological equivalents
NB. - larzxxxx and larfxxxx are equivalent up to τ position
NB. - if τ=0 then v can have any element values

larzlcbc=: _1 larxlcxc
larzlcbr=: _1 larxlcxr
larzlcfc=:  0 larxlcxc
larzlcfr=:  0 larxlcxr

larzlnbc=: _1 larxlnxc
larzlnbr=: _1 larxlnxr
larzlnfc=:  0 larxlnxc
larzlnfr=:  0 larxlnxr

larzrcbc=: _1 larxrcxc
larzrcbr=: _1 larxrcxr
larzrcfc=:  0 larxrcxc
larzrcfr=:  0 larxrcxr

larzrnbc=: _1 larxrnxc
larzrnbr=: _1 larxrnxr
larzrnfc=:  0 larxrnxc
larzrnfr=:  0 larxrnxr

NB. ---------------------------------------------------------
NB. Dyad         Action     Side     Tran    Dir    Layout        eC
NB. larfblcbc    H' * C     left     ct      bwd    columnwise    0, C
NB. larfblcbr    H' * C     left     ct      bwd    rowwise       0, C
NB. larfblcfc    H' * C     left     ct      fwd    columnwise    C, 0
NB. larfblcfr    H' * C     left     ct      fwd    rowwise       C, 0
NB. larfblnbc    H  * C     left     none    bwd    columnwise    0, C
NB. larfblnbr    H  * C     left     none    bwd    rowwise       0, C
NB. larfblnfc    H  * C     left     none    fwd    columnwise    C, 0
NB. larfblnfr    H  * C     left     none    fwd    rowwise       C, 0
NB. larfbrcbc    C  * H'    right    ct      bwd    columnwise    0,.C
NB. larfbrcbr    C  * H'    right    ct      bwd    rowwise       0,.C
NB. larfbrcfc    C  * H'    right    ct      fwd    columnwise    C,.0
NB. larfbrcfr    C  * H'    right    ct      fwd    rowwise       C,.0
NB. larfbrnbc    C  * H     right    none    bwd    columnwise    0,.C
NB. larfbrnbr    C  * H     right    none    bwd    rowwise       0,.C
NB. larfbrnfc    C  * H     right    none    fwd    columnwise    C,.0
NB. larfbrnfr    C  * H     right    none    fwd    rowwise       C,.0
NB.
NB. Description:
NB.   Build and apply a block reflector H or its transpose H'
NB.   to a matrix, from either the left or the right. H is
NB.   defined by pair (V,Τ), where Τ is the triangular factor
NB.   produced from pair (V,Τ) by larftxx.
NB.
NB. Syntax:
NB.   eCupd=. VTau larfbxxxx eC
NB. where
NB.   eC    - matrix C augmented by trash vector
NB.   VTau  - matrix V augmented by vector Tau
NB.   eCupd - an updated eC
NB.   V     - unit trapezoidal matrix
NB.   Tau   - k-vector τ[0:k-1] corresp. to V
NB.
NB. Notes:
NB. - models sequence of calls to LAPACK's xLARFT and then
NB.   to xLARFB
NB. - larfxxxx and larfbxxxx are topological equivalents
NB. - larfbxxxx and larzbxxxx are equivalent up to Tau
NB.   position

larfblcbc=:          0  larxblcxc larftbc
larfblcbr=: (< a: ;  0) larxblcxr larftbr
larfblcfc=:         _1  larxblcxc larftfc
larfblcfr=: (< a: ; _1) larxblcxr larftfr

larfblnbc=:          0  larxblnxc larftbc
larfblnbr=: (< a: ;  0) larxblnxr larftbr
larfblnfc=:         _1  larxblnxc larftfc
larfblnfr=: (< a: ; _1) larxblnxr larftfr

larfbrcbc=:          0  larxbrcxc larftbc
larfbrcbr=: (< a: ;  0) larxbrcxr larftbr
larfbrcfc=:         _1  larxbrcxc larftfc
larfbrcfr=: (< a: ; _1) larxbrcxr larftfr

larfbrnbc=:          0  larxbrnxc larftbc
larfbrnbr=: (< a: ;  0) larxbrnxr larftbr
larfbrnfc=:         _1  larxbrnxc larftfc
larfbrnfr=: (< a: ; _1) larxbrnxr larftfr

NB. ---------------------------------------------------------
NB. Dyad         Action     Side     Tran    Dir    Layout        eC
NB. larzblcbc    H' * C     left     ct      bwd    columnwise    C, 0
NB. larzblcbr    H' * C     left     ct      bwd    rowwise       C, 0
NB. larzblcfc    H' * C     left     ct      fwd    columnwise    0, C
NB. larzblcfr    H' * C     left     ct      fwd    rowwise       0, C
NB. larzblnbc    H  * C     left     none    bwd    columnwise    C, 0
NB. larzblnbr    H  * C     left     none    bwd    rowwise       C, 0
NB. larzblnfc    H  * C     left     none    fwd    columnwise    0, C
NB. larzblnfr    H  * C     left     none    fwd    rowwise       0, C
NB. larzbrcbc    C  * H'    right    ct      bwd    columnwise    C,.0
NB. larzbrcbr    C  * H'    right    ct      bwd    rowwise       C,.0
NB. larzbrcfc    C  * H'    right    ct      fwd    columnwise    0,.C
NB. larzbrcfr    C  * H'    right    ct      fwd    rowwise       0,.C
NB. larzbrnbc    C  * H     right    none    bwd    columnwise    C,.0
NB. larzbrnbr    C  * H     right    none    bwd    rowwise       C,.0
NB. larzbrnfc    C  * H     right    none    fwd    columnwise    0,.C
NB. larzbrnfr    C  * H     right    none    fwd    rowwise       0,.C
NB.
NB. Description:
NB.   Build and apply a block reflector H or its transpose H'
NB.   to a matrix, from either the left or the right. H is
NB.   defined by pair (V,Τ), where Τ is the triangular factor
NB.   produced from pair (V,Τ) by larztxx.
NB.
NB. Syntax:
NB.   eCupd=. VTau larzbxxxx eC
NB. where
NB.   eC    - matrix C augmented by trash vector
NB.   VTau  - matrix V augmented by vector Tau
NB.   eCupd - an updated eC
NB.   V     - unit trapezoidal matrix, with 0s in atoms to be
NB.           ignored
NB.   Tau   - k-vector τ[0:k-1] corresp. to V
NB.
NB. Notes:
NB. - models sequence of calls to LAPACK's xLARZT and then
NB.   to xLARZB
NB. - larzxxxx and larzbxxxx are topological equivalents
NB. - larzbxxxx and larfbxxxx are equivalent up to Tau
NB.   position

larzblcbc=:         _1  larxblcxc larztbc
larzblcbr=: (< a: ; _1) larxblcxr larztbr
larzblcfc=:          0  larxblcxc larztfc
larzblcfr=: (< a: ;  0) larxblcxr larztfr

larzblnbc=:         _1  larxblnxc larztbc
larzblnbr=: (< a: ; _1) larxblnxr larztbr
larzblnfc=:          0  larxblnxc larztfc
larzblnfr=: (< a: ;  0) larxblnxr larztfr

larzbrcbc=:         _1  larxbrcxc larztbc
larzbrcbr=: (< a: ; _1) larxbrcxr larztbr
larzbrcfc=:          0  larxbrcxc larztfc
larzbrcfr=: (< a: ;  0) larxbrcxr larztfr

larzbrnbc=:         _1  larxbrnxc larztbc
larzbrnbr=: (< a: ; _1) larxbrnxr larztbr
larzbrnfc=:          0  larxbrnxc larztfc
larzbrnfr=: (< a: ;  0) larxbrnxr larztfr

NB. ---------------------------------------------------------
NB. refga
NB.
NB. Description:
NB.   Conj. to make monad to generate and apply reflector
NB.
NB. Syntax:
NB.   'Aupd vtau'=. (larfxxx refga larfxxxx) A ; isosubA ; isoy ; isoa
NB. where
NB.   larfxxx  - monad to generate a reflector; is called as:
NB.                z=. iso larfxxx y
NB.   larfxxxx - dyad to apply a reflector; is called as:
NB.                subAupd=. vtau larfxxxx subA
NB.   z        - vector, source to produce vector vtau
NB.   A        - m×n-matrix to update, is augmented by trash
NB.              vector according to larfxxxx
NB.   Aupd     - A with subA replaced by subAupd
NB.   subA     - submatrix of A to apply reflection
NB.   subAupd  - subA reflected
NB.   isosubA  - ISO subA (subAupd) within A (Aupd)
NB.   isoy     - ISO within subA of vector y which  defines
NB.              reflector
NB.   ioa      - IO within z of scalar α, usually 0 or _1
NB.   vtau     - vector, produced from z, defines
NB.              reflection matrix

refga=: 2 : 0
  'A isosubA isoy ioa'=. y
  subA=. isosubA { A
  vtau=. 1 ioa} u isoy { subA
  ((vtau v subA) isosubA} A) ; vtau
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testlarfg
NB.
NB. Description:
NB.   Test larfx by general vector
NB.
NB. Syntax:
NB.   log=. testlarfg ey
NB. where
NB.   ey  - (n+1)-vector
NB.   log - 6-vector of boxes, test log

testlarfg=: 3 : 0
  iso=. (2 ?@$ <:@#) y

  log=.          ('larfg' tdyad ((0&{::)`(1&{::)`]`nan`nan`nan)) iso ; y
  log=. log lcat ('larfp' tdyad ((0&{::)`(1&{::)`]`nan`nan`nan)) iso ; y
)

NB. ---------------------------------------------------------
NB. testlarf
NB.
NB. Description:
NB.   Test larfxxxx by general matrix
NB.
NB. Syntax:
NB.   log=. testlarf (trash ; C)
NB. where
NB.   C   - m×n-matrix, is used as multiplier, the 1st row or
NB.         column is used to form reflector
NB.   log - 6-vector of boxes, test log

testlarf=: 3 : 0
  y=. 1 {:: y
  rcond=. nan`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  log=.          ('larflcbc' tdyad ((1&(_1})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan))     y  , ~ 0
  log=. log lcat ('larflcbr' tdyad ((1&(_1})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan)) (ct y) , ~ 0
  log=. log lcat ('larflcfc' tdyad ((1&( 0})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan))     y  ,   0
  log=. log lcat ('larflcfr' tdyad ((1&( 0})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan)) (ct y) ,   0

  log=. log lcat ('larflnbc' tdyad ((1&(_1})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan))     y  , ~ 0
  log=. log lcat ('larflnbr' tdyad ((1&(_1})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan)) (ct y) , ~ 0
  log=. log lcat ('larflnfc' tdyad ((1&( 0})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan))     y  ,   0
  log=. log lcat ('larflnfr' tdyad ((1&( 0})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan)) (ct y) ,   0

  log=. log lcat ('larfrcbc' tdyad ((1&(_1})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan)) (ct y) ,.~ 0
  log=. log lcat ('larfrcbr' tdyad ((1&(_1})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan))     y  ,.~ 0
  log=. log lcat ('larfrcfc' tdyad ((1&( 0})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan)) (ct y) ,.  0
  log=. log lcat ('larfrcfr' tdyad ((1&( 0})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan))     y  ,.  0

  log=. log lcat ('larfrnbc' tdyad ((1&(_1})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan)) (ct y) ,.~ 0
  log=. log lcat ('larfrnbr' tdyad ((1&(_1})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan))     y  ,.~ 0
  log=. log lcat ('larfrnfc' tdyad ((1&( 0})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan)) (ct y) ,.  0
  log=. log lcat ('larfrnfr' tdyad ((1&( 0})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan))     y  ,.  0
)

NB. ---------------------------------------------------------
NB. testlarz
NB.
NB. Description:
NB.   Test larzxxxx by general matrix
NB.
NB. Syntax:
NB.   log=. testlarz (trash ; C)
NB. where
NB.   C   - m×n-matrix, is used as multiplier, the 1st row or
NB.         column is used to form reflector
NB.   log - 6-vector of boxes, test log

testlarz=: 3 : 0
  y=. 1 {:: y
  rcond=. nan`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  log=.          ('larzlcbc' tdyad ((1&( 0})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan))     y  ,   0
  log=. log lcat ('larzlcbr' tdyad ((1&( 0})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan)) (ct y) ,   0
  log=. log lcat ('larzlcfc' tdyad ((1&(_1})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan))     y  , ~ 0
  log=. log lcat ('larzlcfr' tdyad ((1&(_1})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan)) (ct y) , ~ 0

  log=. log lcat ('larzlnbc' tdyad ((1&( 0})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan))     y  ,   0
  log=. log lcat ('larzlnbr' tdyad ((1&( 0})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan)) (ct y) ,   0
  log=. log lcat ('larzlnfc' tdyad ((1&(_1})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan))     y  , ~ 0
  log=. log lcat ('larzlnfr' tdyad ((1&(_1})^:(1 < #)@:({."1))`]`]`(rcond"_)`nan`nan)) (ct y) , ~ 0

  log=. log lcat ('larzrcbc' tdyad ((1&( 0})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan)) (ct y) ,.  0
  log=. log lcat ('larzrcbr' tdyad ((1&( 0})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan))     y  ,.  0
  log=. log lcat ('larzrcfc' tdyad ((1&(_1})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan)) (ct y) ,.~ 0
  log=. log lcat ('larzrcfr' tdyad ((1&(_1})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan))     y  ,.~ 0

  log=. log lcat ('larzrnbc' tdyad ((1&( 0})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan)) (ct y) ,.  0
  log=. log lcat ('larzrnbr' tdyad ((1&( 0})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan))     y  ,.  0
  log=. log lcat ('larzrnfc' tdyad ((1&(_1})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan)) (ct y) ,.~ 0
  log=. log lcat ('larzrnfr' tdyad ((1&(_1})^:(1 < #)@  {.   )`]`]`(rcond"_)`nan`nan))     y  ,.~ 0
)

NB. ---------------------------------------------------------
NB. testlarft
NB.
NB. Description:
NB.   Test larftxx by general matrix
NB.
NB. Syntax:
NB.   log=. testlarft (A ; trash)
NB. where
NB.   A   - m×n-matrix, is used to form Qf
NB.   log - 6-vector of boxes, test log

testlarft=: 3 : 0
  y=. 0 {:: y
  rcond=. nan`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  log=.          ('larftbc' tmonad (((tru1~ -~/@$)@geqlf)`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('larftbr' tmonad (((trl1~ -~/@$)@gerqf)`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('larftfc' tmonad (( trl1        @geqrf)`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('larftfr' tmonad (( tru1        @gelqf)`]`(rcond"_)`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testlarzt
NB.
NB. Description:
NB.   Test larztxx by general matrix
NB.
NB. Syntax:
NB.   log=. testlarzt (A ; trash)
NB. where
NB.   A   - m×n-matrix, is used to form Qf
NB.   log - 6-vector of boxes, test log

testlarzt=: 3 : 0
  y=. 0 {:: y
  rcond=. nan`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  log=.          ('larztbc' tmonad ((((idmat@]`(           i. @])`[)} c)@tzzlf@ trl        )`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('larztbr' tmonad ((((idmat@]`(a: <@;     i. @])`[)} #)@tzrzf@ tru        )`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('larztfc' tmonad ((((idmat@]`(       (-~ i.)@])`[)} c)@tzzrf@(tru~ -~/@$))`]`(rcond"_)`nan`nan)) y
  log=. log lcat ('larztfr' tmonad ((((idmat@]`(a: <@; (-~ i.)@])`[)} #)@tzlzf@(trl~ -~/@$))`]`(rcond"_)`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testlarfb
NB.
NB. Description:
NB.   Test larfbxxxx by general matrix
NB.
NB. Syntax:
NB.   log=. testlarfb (A ; C)
NB. where
NB.   A   - m×n-matrix, is used to form Qf
NB.   C   - m×n-matrix, is used as multiplier
NB.   log - 6-vector of boxes, test log

testlarfb=: 3 : 0
  'A C'=. y
  rcond=. nan`geconi@.(=/@$) C  NB. meaninigful for square matrices only

  Qfbc=. (tru1~ -~/@$) geqlf A
  Qfbr=. (trl1~ -~/@$) gerqf A
  Qffc=.  trl1         geqrf A
  Qffr=.  tru1         gelqf A

  log=.          ('larfblcbc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbc ;     C  , ~ 0
  log=. log lcat ('larfblcbr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbr ; (ct C) , ~ 0
  log=. log lcat ('larfblcfc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffc ;     C  ,   0
  log=. log lcat ('larfblcfr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffr ; (ct C) ,   0

  log=. log lcat ('larfblnbc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbc ;     C  , ~ 0
  log=. log lcat ('larfblnbr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbr ; (ct C) , ~ 0
  log=. log lcat ('larfblnfc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffc ;     C  ,   0
  log=. log lcat ('larfblnfr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffr ; (ct C) ,   0

  log=. log lcat ('larfbrcbc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbc ; (ct C) ,.~ 0
  log=. log lcat ('larfbrcbr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbr ;     C  ,.~ 0
  log=. log lcat ('larfbrcfc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffc ; (ct C) ,.  0
  log=. log lcat ('larfbrcfr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffr ;     C  ,.  0

  log=. log lcat ('larfbrnbc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbc ; (ct C) ,.~ 0
  log=. log lcat ('larfbrnbr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbr ;     C  ,.~ 0
  log=. log lcat ('larfbrnfc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffc ; (ct C) ,.  0
  log=. log lcat ('larfbrnfr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffr ;     C  ,.  0
)

NB. ---------------------------------------------------------
NB. testlarzb
NB.
NB. Description:
NB.   Test larzbxxxx by general matrix
NB.
NB. Syntax:
NB.   log=. testlarzb (A ; C)
NB. where
NB.   A   - m×n-matrix, is used to form Qf
NB.   C   - m×n-matrix, is used as multiplier
NB.   log - 6-vector of boxes, test log

testlarzb=: 3 : 0
  'A C'=. y
  rcond=. nan`geconi@.(=/@$) C  NB. meaninigful for square matrices only

  I=. idmat k=. <./ $ A
  Qfbc=. I (           i.  k)} tzzlf  trl         A
  Qfbr=. I (< a: ;     i.  k)} tzrzf  tru         A
  Qffc=. I (       (-~ i.) k)} tzzrf (tru~ -~/@$) A
  Qffr=. I (< a: ; (-~ i.) k)} tzlzf (trl~ -~/@$) A

  log=.          ('larzblcbc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbc ;     C  ,   0
  log=. log lcat ('larzblcbr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbr ; (ct C) ,   0
  log=. log lcat ('larzblcfc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffc ;     C  , ~ 0
  log=. log lcat ('larzblcfr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffr ; (ct C) , ~ 0

  log=. log lcat ('larzblnbc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbc ;     C  ,   0
  log=. log lcat ('larzblnbr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbr ; (ct C) ,   0
  log=. log lcat ('larzblnfc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffc ;     C  , ~ 0
  log=. log lcat ('larzblnfr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffr ; (ct C) , ~ 0

  log=. log lcat ('larzbrcbc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbc ; (ct C) ,.  0
  log=. log lcat ('larzbrcbr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbr ;     C  ,.  0
  log=. log lcat ('larzbrcfc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffc ; (ct C) ,.~ 0
  log=. log lcat ('larzbrcfr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffr ;     C  ,.~ 0

  log=. log lcat ('larzbrnbc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbc ; (ct C) ,.  0
  log=. log lcat ('larzbrnbr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qfbr ;     C  ,.  0
  log=. log lcat ('larzbrnfc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffc ; (ct C) ,.~ 0
  log=. log lcat ('larzbrnfr' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`nan`nan)) Qffr ;     C  ,.~ 0
)

NB. ---------------------------------------------------------
NB. testref
NB.
NB. Description:
NB.   Adv. to make verb to test larfxxxxx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testref) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testref_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testref_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testref_mt_ 150 200
NB.
NB. Notes:
NB. - non-blocked larfxxxx algos are tested implicitly in
NB.   testgq, testmq, testqf
NB. - larxtxx and larxbxxxx are impractical for large
NB.   matrices

testref=: 1 : 'lcat_mt_@(testlarfg_mt_@u@(2 + {.)`(nolog_mt_`((lcat_mt_@(testlarf_mt_`testlarz_mt_`testlarft_mt_`testlarzt_mt_`testlarfb_mt_`testlarzb_mt_`:0))@(u ; u))@.(200 >: <./))`:0)'
