NB. Eigenvectors
NB.
NB. tgevcxx     Some or all of the left and/or right
NB.             eigenvectors of generalized Schur form
NB. tgevcxxb    Backtransformed left and/or right
NB.             eigenvectors of generalized Schur form
NB.
NB. testtgevc   Test tgevcxxx by general matrices given
NB. testevc     Adv. to make verb to test tgevcxxx by
NB.             matrices of generator and shape given
NB.
NB. Version: 0.6.8 2010-11-30
NB.
NB. Copyright 2010 Igor Zhuravlov
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
NB. tgevcli
NB. tgevcui
NB.
NB. Description:
NB.   Calculate initial parameters for tgevcxxx
NB.
NB. Syntax:
NB.   'small big bignum d0 d1 d2 abnorm abrwork abscale cond1 cond2 abcoeff abcoeffa d bigd'=. tgevcxi SP
NB. where
NB.   SP       - 2×n×n-matrix (S,:P), generalized Schur form,
NB.              produced by hgexxsxx
NB.   small    ≥ 0
NB.   big      > 0
NB.   bignum   > 0
NB.   d0       - 2×n-matrix, diagonals of S,P
NB.   d1       - 2×n-matrix
NB.   d2       - 2×n-matrix
NB.   abnorm   - 2-vector, being (norm1t"2 SP)
NB.   abrwork  - 2×n-matrix, norm1t either of rows of strict
NB.              lower triangular part of S,P (tgevcli), or
NB.              of columns of strict upper triangular part
NB.              of S,P (tgevcui)
NB.   abscale  - 2-vector, some parameters of S,P
NB.   cond1    - n-vector, pre-calculated part of some
NB.              condition
NB.   cond2    - n-vector, some pre-calculated condition
NB.   abcoeff  - 2×n-matrix
NB.   abcoeffa - 2×n-matrix
NB.   d        - n-vector
NB.   bigd     - n-vector, the scaled d

tgevcui=: 3 : 0
  bignum=. % FP_SFMIN * c y
  small=. % FP_PREC * bignum
  big=. % small
  d0=. diag"2 y
  d1=. 0 1 ]`(9&o.) ag d0
  d2=. 0 1 sorim`| ag d1
  temp=. norm1tc"2 y                                                           NB. 2×n-matrix
  abnorm=. >./"1 temp
  abrwork=. temp - sorim d0
  abscale=. % FP_SFMIN >. abnorm
  temp=. % (>./) FP_SFMIN ,: abscale * d2                                      NB. n-vector
  sba=. |. abscale * temp *"1 d1                                               NB. 2×n-matrix
  abcoeff=. abscale * sba
  NB. scale to avoid underflow
  lsab=. *./ 0 1 (>:&FP_SFMIN)`(<&small) ag 0 1 (|`sorim ag)"2 sba ,: abcoeff  NB. 2×n-matrix
  scale=. >./ lsab } 1 ,: (big <. abnorm) * small % 0 1 |`sorim ag sba         NB. n-vector
  scale=. (+./ lsab) } scale ,: scale <. % FP_SFMIN * (>./) 1 , 0 1 |`sorim ag abcoeff
  abcoeff=. lsab } (abcoeff *"1 scale) ,: (abscale * scale *"1 sba)
  abcoeffa=. 0 1 |`sorim ag abcoeff                                            NB. coeffs for triangular solvers
  cond1=. +/ abcoeffa * abrwork
  dmin=. >./ FP_SFMIN , FP_PREC * abnorm * abcoeffa
  d=. + (-/) abcoeff * d0
  d=. (dmin >: sorim d) } d ,: dmin
  abs1d=. sorim d
  cond2=. 1 > abs1d
  bigd=. bignum * abs1d
  small ; big ; bignum ; d0 ; d1 ; d2 ; abnorm ; abrwork ; abscale ; cond1 ; cond2 ; abcoeff ; abcoeffa ; d ; bigd
)

NB. ---------------------------------------------------------
NB. tgevclx
NB. tgevcux
NB.
NB. Description:
NB.   Compute some or all of right eigenvectors
NB.
NB. Syntax:
NB.   vapp=. vbt tgevc?x  #################
NB.   X=.    (ios ; init) tgevc?x A ,: B
NB. where
NB.   ios   - k-vector, lIOS eigenvectors to compute
NB.   init  - 15-vector of boxes, the output of tgevc?i
NB.   AB    - one of following:
NB.           - 2×n×n-matrix SP, generalized Schur form,
NB.             output of hgexxsxx
NB.           - 3×n×n-matrix SP , Q
NB.           - 3×n×n-matrix SP , Z
NB.           - 4×n×n-matrix SP , Q ,: Z
NB.   X     - n×k-matrix, some or all of right eigenvectors
NB.   k     - integer in range [0,n]
NB.
NB. Notes:
NB. - wildcard char is '?' here instead of 'x' used in names

tgevcux=: 1 : 0
:
  'ios small big bignum d0 d1 d2 abnorm abrwork abscale cond1 cond2 abcoeff abcoeffa d bigd'=. x
  n=. c y
  k=. # ios
  X=. (n , 0) $ 0
  je=. <: k
  while. je >: 0 do.
    if. *./ FP_SFMIN >: (je { ios) {"1 d2 do.
      NB. singular matrix pencil - return unit eigenvector
      X=. (1 je } n $ 0) ,. X
    else.
      NB. non-singular eigenvalue: triangular solve of:
      NB.   (a*A - b*B) * x = 0 ,
      NB. columnwise
      NB. work[0:j-1] contains sums w
      NB. work[j+1:je] contains x
      work=. 1 ,~ (-/) ((je { ios) {"1 abcoeff) * (((0 0,]),:(2,],1:)) je { ios) ({.@(2 0 1&|:)) ;. 0 y
      j=. <: je { ios
      while. j >: 0 do.
        NB. form:
        NB.   x[j] = - w[j] / d
        NB. with scaling and perturbation of the denominator
        abs1wj=. sorim j { work
        if. (abs1wj >: (j { bigd)) *. (j { cond2) do.
          work=. work % abs1wj
        end.
        work=. j (-@(%&(j{d))) upd work
        abs1wj=. sorim j { work
        if. j > 0 do.
          NB. w = w + x[j] * (a*S[:,j] - b*P[:,j]) with scaling
          if. ((abcoeffa mp & (j&({"1)) abrwork) >: (bignum % abs1wj)) *. (1 < abs1wj) do.
            work=. work % abs1wj
          end.
          workadd=. (((je { ios) {"1 abcoeff) * j { work) * (((0 0,]),:(2,],1:)) j) ({.@(2 0 1&|:)) ;. 0 y
          work=. (i. j) (+`-/@(,&workadd)) upd work
        end.
        j=. <: j
      end.
      NB. optional back transforming
      work=. ((>: je) ({."1) _1 { y) u work  NB. when backtransform, (ios -: i. n) <=> (je -: je { ios)
      NB. scale eigenvector and stitch X to it
      xmax=. normit work
      if. FP_SFMIN < xmax do.
        X=. (work % xmax) stitcht X
      else.
        X=. 0 ,. X
      end.
    end.
    je=. <: je
  end.
  X
)

NB. ---------------------------------------------------------
NB. tgevcly
NB. tgevcuy
NB.
NB. Description:
NB.   Compute some or all of left eigenvectors
NB.
NB. Syntax:
NB.   vapp=. vbt tgevcuy  #################
NB.   Y=.    (ios ; init) tgevcuy A ,: B
NB. where
NB.   ios   - k-vector, lIOS eigenvectors to compute
NB.   init  - 15-vector of boxes, the output of tgevcxi
NB.   AB    - one of following:
NB.           - 2×n×n-matrix SP, generalized Schur form,
NB.             output of hgexxsxx
NB.           - 3×n×n-matrix SP , Q
NB.           - 3×n×n-matrix SP , Z
NB.           - 4×n×n-matrix SP , Q ,: Z
NB.   Y     - n×k-matrix, some or all of left eigenvectors
NB.   k     - integer in range [0,n]

tgevcuy=: 1 : 0
:
  'ios small big bignum d0 d1 d2 abnorm abrwork abscale cond1 cond2 abcoeff abcoeffa d bigd'=. x
  n=. c y
  k=. # ios
  Y=. (n , 0) $ 0
  je=. 0
  while. je < k do.
    if. *./ FP_SFMIN >: (je { ios) {"1 d2 do.
      NB. singular matrix pencil - return unit eigenvector
      Y=. Y ,. (1 je } n $ 0)
    else.
      NB. non-singular eigenvalue: triangular solve of:
      NB.   (a*A - b*B)^H * y = 0 ,
      NB. rowwise in (a*A - b*B)^H , or columnwise in (a*A - b*B)
      work=. 1
      xmax=. 1
      j=. >: je { ios
      while. j < n do.
        NB. compute:
        NB.         j-1
        NB.   sum = Σ  conjg(a*S[k,j] - b*P[k,j]) * x[k] ,
        NB.         k=je
        NB. scale if necessary
        if. (j { cond1) > (bignum % xmax) do.
          work=. work % xmax
          xmax=. 1
        end.
        sum=. -/ ((je { ios) {"1 abcoeff) * work mp"1 (j ((0,,~),:(2,-,1:)) je { ios) (+@{.@(2 0 1&|:)) ;. 0 y
        NB. form:
        NB.   x[j] = - sum / conjg(a*S[j,j] - b*P[j,j])
        NB. with scaling and perturbation of the denominator
        abs1sum=. sorim sum
        if. (abs1sum >: j { bigd) *. (j { cond2) do.
          work=. work % abs1sum
          xmax=. xmax % abs1sum
          sum=. sum % abs1sum
        end.
        workj=. - sum % j { d
        work=. work , workj
        xmax=. xmax >. sorim workj
        j=. >: j
      end.
      NB. optional back transforming
      work=. (je (}."1) 2 ({ :: 0:) y) u work  NB. when backtransform, (ios -: i. n) <=> (je -: je { ios)
      NB. scale eigenvector and stitch to Y
      xmax=. normit work
      if. FP_SFMIN < xmax do.
        Y=. Y stitchb work % xmax
      else.
        Y=. Y ,. 0
      end.
    end.
    je=. >: je
  end.
  Y
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB.   Y=.     [ios] tgevcxl  S ,: P
NB.   X=.     [ios] tgevcxr  S ,: P
NB.   'Y X'=. [ios] tgevcxb  S ,: P

tgevcll=: ($:~ (i.@c)) : ((; tgevcli) (] tgevcly) ])
tgevclr=: ($:~ (i.@c)) : ((; tgevcli) (] tgevclx) ])
tgevclb=: ($:~ (i.@c)) : ((; tgevcli) ((] tgevcly) ,: (] tgevclx)) ])

tgevcul=: ($:~ (i.@c)) : ((; tgevcui) (] tgevcuy) ])
tgevcur=: ($:~ (i.@c)) : ((; tgevcui) (] tgevcux) ])
tgevcub=: ($:~ (i.@c)) : ((; tgevcui) ((] tgevcuy) ,: (] tgevcux)) ])

NB. ---------------------------------------------------------
NB.   Y=.     tgevcxlb S , P ,: Q
NB.   X=.     tgevcxrb S , P ,: Z
NB.   'Y X'=. tgevcxbb S , P , Q ,: Z

tgevcllb=: (((i.@c) (; tgevcli) (2&{.)) (mp tgevcly) ]) : [:
tgevclrb=: (((i.@c) (; tgevcli) (2&{.)) (mp tgevclx) ]) : [:
tgevclbb=: (((i.@c) (; tgevcli) (2&{.)) ((mp tgevcly) ,: (mp tgevclx)) ]) : [:

tgevculb=: (((i.@c) (; tgevcui) (2&{.)) (mp tgevcuy) ]) : [:
tgevcurb=: (((i.@c) (; tgevcui) (2&{.)) (mp tgevcux) ]) : [:
tgevcubb=: (((i.@c) (; tgevcui) (2&{.)) ((mp tgevcuy) ,: (mp tgevcux)) ]) : [:

NB. =========================================================
NB. Test suite
