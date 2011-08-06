NB. Eigenvectors
NB.
NB. tgevcxx     Some or all of the right and/or left
NB.             eigenvectors of generalized Schur form
NB. tgevcxxb    Backtransformed right and/or left
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
NB.   'ios small big bignum d0 d1 d2 abnorm abrwork abscale cond1 cond2 abcoeff abcoeffa d bigd'=. ios tgevcxi SP
NB. where
NB.   SP       - 2×n×n-matrix (S,:P), generalized Schur form,
NB.              produced by hgexxsxx
NB.   ios      - k-vector, lIOS eigenvectors to process
NB.   small    ≥ 0
NB.   big      > 0
NB.   bignum   > 0
NB.   d0       - 2×k-matrix, diagonals of S,P
NB.   d1       - 2×k-matrix
NB.   d2       - 2×k-matrix
NB.   abnorm   - 2-vector, some parameters of S,P
NB.   abrwork  - 2×k-matrix, norms either of rows of strict
NB.              lower triangular part of S,P (tgevcli), or
NB.              of columns of strict upper triangular part
NB.              of S,P (tgevcui)
NB.   abscale  - 2-vector, some parameters of S,P
NB.   cond1    - k-vector, pre-calculated part of some
NB.              condition
NB.   cond2    - k-vector, some pre-calculated condition
NB.   abcoeff  - 2×k-matrix
NB.   abcoeffa - 2×k-matrix
NB.   d        - k-vector
NB.   bigd     - k-vector, the scaled d
NB.   k        - integer in range [0,n]

tgevcui=: 4 : 0
  bignum=. % FP_SFMIN * c y
  small=. % FP_PREC * bignum
  big=. % small
  d0=. x {"1 diag"2 y
  d1=. (]`(9&o.) ag) d0
  d2=. (sorim`| ag) d1
  temp=. x {"1 norm1tc"2 y                                                                 NB. 2×k-matrix
  abnorm=. >./"1 temp
  abrwork=. temp - sorim d0
  abscale=. % FP_SFMIN >. abnorm
  temp=. % (>./) FP_SFMIN ,: abscale * d2                                                  NB. k-vector
  sba=. |. abscale * temp *"1 d1                                                           NB. 2×k-matrix
  abcoeff=. abscale * sba
  NB. scale
  lsab=. *./ ((>:&FP_SFMIN)`(<&small) ag) (|`sorim ag)"2 sba ,: abcoeff                    NB. 2×k-matrix
  scale=. (*./ lsab) } (] ,: (>./)) lsab } 1 , (big <. abnorm) * small % (|`sorim ag) sba  NB. k-vector
  scale=. (+./ lsab) } scale ,: scale <. % FP_SFMIN * (>./) 1 , (|`sorim ag) abcoeff
  abcoeff=. lsab } (abcoeff * scale) ,: (abscale * scale *"1 sba)
  abcoeffa=. |`sorim ag abcoeff
  cond1=. +/ abcoeffa * abrwork
  dmin=. >./ FP_SFMIN , FP_PREC * abnorm * abcoeffa
  d=. + (+/) (1 _1 * abcoeff) * d0
  d=. (dmin >: sorim d) } d ,: dmin
  abs1d=. sorim d
  cond2=. 1 > abs1d
  bigd=. bignum * abs1d
  x ; small ; big ; bignum ; d0 ; d1 ; d2 ; abnorm ; abrwork ; abscale ; cond1 ; cond2 ; abcoeff ; abcoeffa ; d ; bigd
)

NB. ---------------------------------------------------------
NB. tgevclx
NB. tgevcux
NB.
NB. Description:
NB.   Compute some or all of right eigenvectors
NB.
NB. Syntax:
NB.   vapp=. vbt tgevcux  #################
NB.   subX=. (ios ; init) tgevclx ios {   SP
NB.   subX=. (ios ; init) tgevcux ios {"1 SP
NB.   X=.    (ios ; init) tgevcux AB
NB. where
NB.   ios   - k-vector, 
NB.   init  - 15-vector of boxes, the output of tgevcxi
NB.   AB    - one of following:
NB.           - 2×n×n-matrix SP, generalized Schur form
NB.             (output of hgexxsxx)
NB.           - 3×n×n-matrix SP , Q
NB.           - 3×n×k-matrix SP , Z
NB.           - 4×n×n-matrix SP , Q ,: Z
NB.   subX  - n×k-matrix, right eigenvectors
NB.   X     - n×n-matrix, right eigenvectors
NB.   k     - integer in range [0,n]

tgevcux=: 1 : 0
:
  'ios small big bignum d0 d1 d2 abnorm abrwork abscale cond1 cond2 abcoeff abcoeffa d bigd'=. x
  'n k'=. }. $ y
  X=. (n , 0) $ 0
  je=. <: k
  while. je >: 0 do.
    if. *./ FP_SFMIN >: je {"1 d2 do.
      X=. (1 (c X) } n $ 0) ,. X
    else.
      NB. triang solver
      work=. 1 ,~ +/ (1 _1 * abcoeff) * (((0 0,]),:(2,],1:)) je { x) ({.@(2 0 1&|:)) ,. 0 y
      j=. <: je { x
      while. j >: 0 do.
        NB. form x[j]
        if. ((sorim j { work) >: (j { bigd)) *. (j { cond2) do.  NB. FIXME!!! (j { work) may outbound
          work=. work % sorim j { work
        end.
        work=. j (-@(%&(j{d))) upd work
        if. j > 0 do.
          if. ((abcoeffa mp j {"1 abrwork) >: (bignum % sorim j { work)) *. (1 < sorim j { work) do.
            work=. work % sorim j { work
          end.
          workadd=. (1 _1 * abcoeff * j { work) * (((0 0,]),:(2,],1:)) j) ] ;. 0 y
          work=. (i. j) (+/@(,&workadd)) upd work
        end.
        j=. <: j
      end.
      NB. back transform
      work=. ((je { x) ({."1) _1 { y) u work  NB. when backtransform, (ios -: i. n) <=> (je -: je { x)
      NB. copy and scale eigenvector into X
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
NB.   subY=. (ios ; init) tgevcly ios {   SP
NB.   subY=. (ios ; init) tgevcuy ios {"1 SP
NB.   Y=.    (ios ; init) tgevcuy AB
NB. where
NB.   ios   - k-vector, 
NB.   init  - 15-vector of boxes, the output of tgevcxi
NB.   AB    - one of following:
NB.           - 2×n×n-matrix SP, generalized Schur form
NB.             (output of hgexxsxx)
NB.           - 3×n×n-matrix SP , Q
NB.           - 3×n×k-matrix SP , Z
NB.           - 4×n×n-matrix SP , Q ,: Z
NB.   subY  - n×k-matrix, left eigenvectors
NB.   Y     - n×n-matrix, left eigenvectors
NB.   k     - integer in range [0,n]

tgevcuy=: 1 : 0
:
  'ios small big bignum d0 d1 d2 abnorm abrwork abscale cond1 cond2 abcoeff abcoeffa d bigd'=. x
  'n k'=. }. $ y
  Y=. (n , 0) $ 0
  je=. 0
  while. je < k do.
    if. *./ FP_SFMIN >: je {"1 d2 do.
      Y=. Y ,. (1 (c Y) } n $ 0)
    else.
      NB. triang solver
      work=. 1
      xmax=. 1
      j=. >: je { x
      while. j < n do.
        NB. compute sum, scale if necessary
        if. (j { cond1) > (bignum % xmax) do.  NB. FIXME!!! (j { cond1) may outbound
          work=. work % xmax
          xmax=. 1
        end.
        sum=. +/ (1 _1 * abcoeff) * work mp"1 (j ((0,,~),:(2,-,1:)) je { x) (+@{.@(2 0 1&|:)) ,. 0 y
        NB. form x[j]
        if. ((sorim sum) >: (j { bigd)) *. (j { cond2) do.
          abs1sum=. sorim sum
          work=. work % abs1sum
          xmax=. xmax % abs1sum
          sum=. sum % abs1sum
        end.
        workj=. - sum % j { d
        work=. work , workj
        xmax=. xmax >. sorim workj
        j=. >: j
      end.
      NB. back transform
      work=. ((je { x) (}."1) 2 ({ :: 0:) y) u work  NB. when backtransform, (ios -: i. n) <=> (je -: je { x)
      NB. copy and scale eigenvector into Y
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
NB.   Y=.  [ios] tgevcxl  S ,: P
NB.   X=.  [ios] tgevcxr  S ,: P
NB.   YX=. [ios] tgevcxb  S ,: P

tgevcll=: ($:~ (i.@c)) : (tgevcli (] tgevcly) {)
tgevclr=: ($:~ (i.@c)) : (tgevcli (] tgevclx) {)
tgevclb=: ($:~ (i.@c)) : (tgevcli ((] tgevcly) ,: (] tgevclx)) {)

tgevcul=: ($:~ (i.@c)) : (tgevcui (] tgevcuy) ({"1))
tgevcur=: ($:~ (i.@c)) : (tgevcui (] tgevcux) ({"1))
tgevcub=: ($:~ (i.@c)) : (tgevcui ((] tgevcuy) ,: (] tgevcux)) ({"1))

NB. ---------------------------------------------------------
NB.   Y=.  tgevclxb S , P ,: Q
NB.   X=.  tgevclxb S , P ,: Z
NB.   YX=. tgevclxb S , P , Q ,: Z

tgevcllb=: (((tgevcli~ (i.@c)) @ (2&{.)) (mp tgevcly) ]) : [:
tgevclrb=: (((tgevcli~ (i.@c)) @ (2&{.)) (mp tgevclx) ]) : [:
tgevclbb=: (((tgevcli~ (i.@c)) @ (2&{.)) ((mp tgevcly) ,: (mp tgevclx)) ]) : [:

tgevculb=: (((tgevcui~ (i.@c)) @ (2&{.)) (mp tgevcuy) ]) : [:
tgevcurb=: (((tgevcui~ (i.@c)) @ (2&{.)) (mp tgevcux) ]) : [:
tgevcubb=: (((tgevcui~ (i.@c)) @ (2&{.)) ((mp tgevcuy) ,: (mp tgevcux)) ]) : [:

NB. =========================================================
NB. Test suite
