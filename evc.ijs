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
NB. tgevci
NB.
NB. Description:
NB.   Calculate initial parameters for tgevclxx
NB.
NB. Syntax:
NB.   'bignum d0 d1 d2 abnorm abrwork abscale cond1 cond2 abcoeff abcoeffa d'=. ios tgevci SP
NB. where
NB.   SP       - 2×n×n-matrix (S,:P), generalized Schur form,
NB.              produced by hgezqsxx
NB.   ios      - k-vector, lIOS eigenvectors to compute
NB.   bignum   > 0
NB.   d0       - 2×n-matrix, diagonals of S and P
NB.   d1       - 2×n-matrix
NB.   d2       - 2×n-matrix
NB.   abnorm   - 2-vector, being (norm1t"2 SP)
NB.   abrwork  - 2×n-matrix, norm1t of rows of strict lower
NB.              triangular part of S and P
NB.   abscale  - 2-vector
NB.   cond1    - n-vector, pre-calculated part of some
NB.              condition
NB.   cond2    - k×n-matrix, pre-calculated part of some
NB.              condition
NB.   abcoeff  - 2×n-matrix, (a,b) coeffs for pencil a*S-b*P
NB.   abcoeffa - 2×n-matrix, coeffs for triangular solvers
NB.   d        - k×n-matrix

tgevci=: 4 : 0
  bignum=. % FP_SFMIN * c y
  small=. % FP_PREC * bignum
  d0=. diag"2 y
  d1=. 0 1 ]`(9&o.) ag d0
  d2=. 0 1 sorim`| ag d1
  temp=. norm1tr"2 y                                                           NB. 2×n-matrix
  abnorm=. >./"1 temp
  abrwork=. temp - sorim d0
  abscale=. % FP_SFMIN >. abnorm
  temp=. % (>./) FP_SFMIN , abscale * d2                                       NB. n-vector
  sba=. |. abscale * temp *"1 d1                                               NB. 2×n-matrix
  abcoeff=. abscale * sba
  NB. scale to avoid underflow
  lsab=. *./ 0 1 (>:&FP_SFMIN)`(<&small) ag 0 1 (|`sorim ag)"2 sba ,: abcoeff  NB. 2×n-matrix
  scale=. >./ lsab } 1 ,: ((% small) <. abnorm) * small % 0 1 |`sorim ag sba   NB. n-vector
  scale=. (+./ lsab) } scale ,: scale <. % FP_SFMIN * (>./) 1 , 0 1 |`sorim ag abcoeff
  abcoeff=. lsab } (abcoeff *"1 scale) ,: (abscale * scale *"1 sba)
  abcoeffa=. 0 1 |`sorim ag abcoeff
  cond1=. +/ abcoeffa * abrwork
  dmin=. (# x) # ,: >./ FP_SFMIN , FP_PREC * abnorm * abcoeffa                 NB. k×n-matrix
  d=. (x { |: abcoeff) mp 0 1 ]`- ag d0
  d=. (dmin < sorim d) } dmin ,: d
  cond2=. bignum (((1>]),.*) sorim) d

  bignum ; d0 ; d1 ; d2 ; abnorm ; abrwork ; abscale ; cond1 ; cond2 ; abcoeff ; abcoeffa ; d
)

NB. ---------------------------------------------------------
NB. tgevcly
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
NB.
NB. Assertions:
NB.   (tgevcurb SP , Z) -: ct (tgevcllb ct"2 SP , Q)
NB.   (tgevculb SP , Z) -: ct (tgevclrb ct"2 SP , Q)

tgevcly=: 4 : 0
  'ios bignum d0 d1 d2 abnorm abrwork abscale cond1 cond2 abcoeff abcoeffa d'=. x
  n=. c y
  k=. # ios
  W=. (0,n) $ 0
  je=. <: k
  while. je >: 0 do.
    if. *./ FP_SFMIN >: (je { ios) {"1 d2 do.
      NB. singular matrix pencil - return unit eigenvector
      work=. 1 je } n $ 0
    else.
      NB. non-singular eigenvalue: triangular solve of:
      NB.   x * (a*A - b*B)^H = 0 ,
      NB. columnwise in (a*A - b*B)^H, rowwise in (a*A - b*B)
      NB. work[0:j-1] contains sums w
      NB. work[j+1:je] contains x
      work=. 1 ,~ (-/) ((je { ios) {"1 abcoeff) * (((0,],0:),:(2 1,])) je { ios) ({.@(1 0 2&|:)) ;. 0 y
      di=. je { d
      j=. <: je { ios
      while. j >: 0 do.
        NB. form:
        NB.   x[j] = - w[j] / di
        NB. with scaling and perturbation of the denominator
        abs1wj=. sorim j { work
        if. *.`<:/ (j { cond2) , abs1wj do.
          work=. work % abs1wj
        end.
        work=. j (-@(%&(j{di))) upd work
        abs1wj=. sorim j { work
        if. j > 0 do.
          NB. w = w + x[j] * (a*S[:,j] - b*P[:,j]) with scaling
          if. ((abcoeffa mp & (j&({"1)) abrwork) >: (bignum % abs1wj)) *. (1 < abs1wj) do.
            work=. work % abs1wj
          end.
          workadd=. (((je { ios) {"1 abcoeff) * j { work) * (((0,],0:),:(2 1,])) j) ({.@(1 0 2&|:)) ;. 0 y
          work=. (i. j) (+`-/@(,&workadd)) upd work
        end.
        j=. <: j
      end.
    end.
    je=. <: je
    W=. work , W
  end.
)

NB. scale
NB. V=. tgevcs W

tgevcs=: 3 : 0
  norm=. normitr y
  ios=. (#y) #"0 FP_SFMIN < norm
  y=. y % norm
  y=. ios } 0 ,: y
)

NB. NB. optional back transforming
NB. work=. ((>: je) {. _1 { y) mp~ work  NB. when backtransform, (ios -: i. n) <=> (je -: je { ios)
NB. NB. scale eigenvector and append W to it
NB. W=. ((0:`%@.(FP_SFMIN<]) normit) work) appendl W

NB. ---------------------------------------------------------
NB. tgevclx
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

tgevclx=: 4 : 0
  'ios bignum d0 d1 d2 abnorm abrwork abscale cond1 cond2 abcoeff abcoeffa d'=. x
  d=. + d
  n=. c y
  k=. # ios
  W=. (0,n) $ 0
  je=. 0
  while. je < k do.
    if. *./ FP_SFMIN >: (je { ios) {"1 d2 do.
      NB. singular matrix pencil - return unit eigenvector
      work=. 1 je } n $ 0
    else.
      NB. non-singular eigenvalue: triangular solve of:
      NB.   y * (a*A - b*B) = 0 ,
      NB. colwise in (a*A - b*B) , or rowwise in (a*A - b*B)^H
      work=. 1
      xmax=. 1
      di=. je { d
      j=. >: je { ios
      while. j < n do.
        NB. compute:
        NB.         j-1
        NB.   sum = Σ  conjg(a*S[j,k] - b*P[j,k]) * x[k] ,
        NB.         k=je
        NB. scale if necessary
        if. (j { cond1) > (bignum % xmax) do.
          work=. work % xmax
          xmax=. 1
        end.
        sum=. -/ (0 1 ]`+ ag (je { ios) {"1 abcoeff) * work mp"1 (j ((0,, ),:(2 1,- )) je { ios) (+@{.@(1 0 2&|:)) ;. 0 y
        NB. form:
        NB.   x[j] = - sum / conjg(a*S[j,j] - b*P[j,j])
        NB. with scaling and perturbation of the denominator
        abs1sum=. sorim sum
        if. *.`<:/ (j { cond2) , abs1sum do.
          work=. work % abs1sum
          xmax=. xmax % abs1sum
          sum=. sum % abs1sum
        end.
        workj=. - sum % j { di
        work=. work , workj
        xmax=. xmax >. sorim workj
        j=. >: j
      end.
      work=. (-n) {. work
    end.
    je=. >: je
    W=. W , work
  end.
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB.   Y=.     [ios] tgevcxl  S ,: P
NB.   X=.     [ios] tgevcxr  S ,: P
NB.   'Y X'=. [ios] tgevcxb  S ,: P

tgevcll=:  ($:~ i.@c) : (([;tgevci) tgevcs  @ tgevcly           ])
tgevclr=:  ($:~ i.@c) : (([;tgevci) tgevcs  @          tgevclx  ])
tgevclb=:  ($:~ i.@c) : (([;tgevci) tgevcs"2@(tgevcly,:tgevclx) ])

tgevcul=: (($:~ i.@c) : (([;tgevci) tgevcs  @ tgevclx           ])) &.: (|:"2)
tgevcur=: (($:~ i.@c) : (([;tgevci) tgevcs  @(         tgevcly) ])) &.: (|:"2)
tgevcub=: (($:~ i.@c) : (([;tgevci) tgevcs"2@(tgevclx,:tgevcly) ])) &.: (|:"2)

NB. ---------------------------------------------------------
NB.   Y=.     tgevcxlb S , P ,: Q
NB.   X=.     tgevcxrb S , P ,: Z
NB.   'Y X'=. tgevcxbb S , P , Q ,: Z
NB.
NB. Assertions:

tgevcllb=:  (i.@c ([;tgevci) 2&{.) tgevcs  @( tgevcly mp 2{]                    ) ]
tgevclrb=:  (i.@c ([;tgevci) 2&{.) tgevcs  @(                   tgevclx mp _1{] ) ]
tgevclbb=:  (i.@c ([;tgevci) 2&{.) tgevcs"2@((tgevcly mp 2{]),:(tgevclx mp _1{])) ]

tgevculb=: ((i.@c ([;tgevci) 2&{.) tgevcs  @( tgevclx mp 2{]                    ) ]) &.: (|:"2)
tgevcurb=: ((i.@c ([;tgevci) 2&{.) tgevcs  @(                   tgevcly mp _1{] ) ]) &.: (|:"2)
tgevcubb=: ((i.@c ([;tgevci) 2&{.) tgevcs"2@((tgevclx mp 2{]),:(tgevcly mp _1{])) ]) &.: (|:"2)

NB. =========================================================
NB. Test suite
