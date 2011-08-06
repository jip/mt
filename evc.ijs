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
NB.   'bignum d2 abrwork cond1 cond2 abcoeff abcoeffa d'=. ios tgevci SP
NB. where
NB.   SP       - 2×n×n-matrix (S,:P), generalized Schur form,
NB.              produced by hgezqsxx
NB.   ios      - k-vector, lIOS eigenvectors to compute
NB.   bignum   > 0
NB.   d2       - n×2-matrix
NB.   abrwork  - n×2-matrix, stitched norm1t of rows of
NB.              strict lower triangular part of S and P
NB.   cond1    - n-vector, pre-calculated part of some
NB.              condition
NB.   cond2    - k×n-matrix, pre-calculated part of some
NB.              condition
NB.   abcoeff  - n×2-matrix, (a,b) coeffs for pencil a*S-b*P
NB.   abcoeffa - n×2-matrix, coeffs for triangular solvers
NB.   d        - k×n-matrix

tgevci=: 4 : 0
  bignum=. % FP_SFMIN * c y
  small=. % FP_PREC * bignum
  d0=. diag"2 y                                                                NB. 2×n-matrix
  d1=. 0 1 ]`(9&o.) ag d0                                                      NB. 2×n-matrix
  d2=. 0 1 sorim`| ag d1
  temp=. norm1tr"2 y                                                           NB. 2×n-matrix
  abnorm=. >./"1 temp                                                          NB. 2-vector
  abrwork=. temp - sorim d0
  abscale=. % FP_SFMIN >. abnorm                                               NB. 2-vector
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

  bignum ; (|:d2) ; (|:abrwork) ; cond1 ; cond2 ; (|:abcoeff) ; (|:abcoeffa) ; d
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
  'ios bignum d2 abrwork cond1 cond2 abcoeff abcoeffa d'=. x
  n=. c y
  k=. # ios
  W=. (0,n) $ 0
  je=. <: k
  while. je >: 0 do.
    if. *./ FP_SFMIN >: (je { ios) { d2 do.
      NB. singular matrix pencil - return unit eigenvector
      work=. 1 je } n $ 0
    else.
      NB. non-singular eigenvalue: triangular solve of:
      NB.   x * (a*A - b*B)^H = 0 ,
      NB. columnwise in (a*A - b*B)^H, rowwise in (a*A - b*B)
      NB. work[0:j-1] contains sums w
      NB. work[j+1:je] contains x
      work=. 1 ,~ (-/) ((je { ios) { abcoeff) * (((0,],0:),:(2 1,])) je { ios) ({.@(1 0 2&|:)) ;. 0 y
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
          if. ((abcoeffa mp & (j&{) abrwork) >: (bignum % abs1wj)) *. (1 < abs1wj) do.
            work=. work % abs1wj
          end.
          workadd=. (((je { ios) { abcoeff) * j { work) * (((0,],0:),:(2 1,])) j) ({.@(1 0 2&|:)) ;. 0 y
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
  'ios bignum d2 abrwork cond1 cond2 abcoeff abcoeffa d'=. x
  d=. + d
  n=. c y
  k=. # ios
  W=. (0,n) $ 0
  je=. 0
  while. je < k do.
    if. *./ FP_SFMIN >: (je { ios) { d2 do.
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
        sum=. -/ (0 1 ]`+ ag (je { ios) { abcoeff) * work mp (j ((0,, ),:(2 1,- )) je { ios) (+@{.@(0&|:)) ;. 0 y
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

tgevcul=: ($:~ i.@c) : ((([;tgevci) tgevcs  @ tgevclx           ]) &.: (|:"2))
tgevcur=: ($:~ i.@c) : ((([;tgevci) tgevcs  @          tgevcly  ]) &.: (|:"2))
tgevcub=: ($:~ i.@c) : ((([;tgevci) tgevcs"2@(tgevclx,:tgevcly) ]) &.: (|:"2))

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

NB. ---------------------------------------------------------
NB. testtgevc
NB.
NB. Description:
NB.   Test tgevcxxx by general matrices given
NB.
NB. Syntax:
NB.   testtgevc AB
NB. where
NB.   AB - 2×n×n-report
NB.
NB. Formula:
NB.   berr := max(berr0,berr1)
NB. where
NB.   ||M|| := max(||M||_1 , FP_SFMIN)
NB.   ||v|| := max(|Re(v(i))|+Im(v(i))|)
NB.   β     - machine precision
NB.   α(i)  - i-th eigenvalue, also i-th element on S
NB.           diagonal
NB.   β(i)  - i-th eigenvalue, also i-th element on P
NB.           diagonal
NB.   l(i)  - i-th left eigenvector, also i-th column of L
NB.   lb(i) - i-th back transformed left eigenvector, also
NB.           i-th column of LB
NB.   r(i)  - i-th right eigenvector, also i-th column of R
NB.   rb(i) - i-th back transformed right eigenvector, also
NB.           i-th column of RB
NB.   S P dQ1 dZ1'=. (0,n) hgexxsvv H , T , ,:~ I
NB.   - tgevcll:
NB.       berr0 := 
NB.       berr1 := (max(||l(i)||) - 1) / (ulp * n)
NB.   - tgevclr:
NB.       berr0 := 
NB.       berr1 := (max(||r(i)||) - 1) / (ulp * n)
NB.   - tgevclb:
NB.       berr0 := berr(tgevcll)
NB.       berr1 := berr(tgevclr)
NB.   - tgevcllb:
NB.       berr0 := 
NB.       berr1 := (max(||lb(i)||) - 1) / (ulp * n)
NB.   - tgevclrb:
NB.       berr0 := 
NB.       berr1 := (max(||rb(i)||) - 1) / (ulp * n)
NB.   - tgevclbb:
NB.       berr0 := berr(tgevcllb)
NB.       berr1 := berr(tgevclrb)
NB.   - tgevcul:
NB.       berr0 := max(||(β(i)*S - α(i)*P)^H * l(i)|| / (ulp * max(||β(i)*S||,||α(i)*P||)))
NB.       berr1 := (max(||l(i)||) - 1) / (ulp * n)
NB.   - tgevcur:
NB.       berr0 := max(||(β(i)*S - α(i)*P)   * r(i)|| / (ulp * max(||β(i)*S||,||α(i)*P||)))
NB.       berr1 := (max(||r(i)||) - 1) / (ulp * n)
NB.   - tgevcub:
NB.       berr0 := berr(tgevcul)
NB.       berr1 := berr(tgevcur)
NB.   - tgevculb:
NB.       berr0 := max(||(β(i)*H - α(i)*T)^H * l(i)|| / (ulp * max(||β(i)*H||,||α(i)*T||)))
NB.       berr1 := (max(||lb(i)||) - 1) / (ulp * n)
NB.   - tgevcurb:
NB.       berr0 := max(||(β(i)*H - α(i)*T)   * r(i)|| / (ulp * max(||β(i)*H||,||α(i)*T||)))
NB.       berr1 := (max(||rb(i)||) - 1) / (ulp * n)
NB.   - tgevcubb:
NB.       berr0 := berr(tgevculb)
NB.       berr1 := berr(tgevcurb)

############
   I=. idmat_mt_ 7
   'H T'=. 0 7 gghrdunn_mt_ ((unmqrlc_mt_~ ,: tru_mt_@]) geqrf_mt_)/ AB
   'S P Q2 Z2'=. 0 7 hgeqzsvv_mt_ H , T , ,:~ I

   L=. tgevcul_mt_ S ,: P
   (ct_mt_ ((i { diag_mt_ P) * S) - ((i { diag_mt_ S) * P)) mp i {"1 L
   LB=. tgevculb_mt_ S , P ,: Q2
   (ct_mt_ ((i { diag_mt_ P) * H) - ((i { diag_mt_ S) * T)) mp i {"1 LB
   R=. tgevcur_mt_ S ,: P
   (((i { diag_mt_ P) * S) - ((i { diag_mt_ S) * P)) mp i {"1 R
   RB=. tgevcurb_mt_ S , P ,: Z2
   (((i { diag_mt_ P) * H) - ((i { diag_mt_ S) * T)) mp i {"1 RB
############
   NB. HT=. mkHT AB
   mkHT=. geqrf@(1{]) ((0,#@]) gghrdunn unmqrlc ,: tru@[) (0{])
   NB. SPQ2Z2=. mkSPQ2Z2 HT
   mkSPQ2Z2=. ((0,]) hgeqzsvv (, (,:~)@idmat)) c
############
   berr=. (vberr vgeto @ tgevcul @ vgety) AB                     berr <- AB
   vgety=. mkSPQ2Z2 @ mkHT                                       SPQ2Z2 <- AB
   tgevcul                                                       L <- SPQ2Z2
  BAD
############
y=. SPQ2Z2
   SPQ2Z2=. mkSPQ2Z2 @ mkHT AB
   vgety=. 2&{.                                                  SP <- SPQ2Z2
   tgevcul                                                       Y <- SP
   vgeto=. ]                                                     Y <- Y
   vberr=.                                                       berr <- SPQ2Z2 vberr Y
   vEnorm=. (normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@]
   
    ba=. |:@|.@:(diag"2) SP
    SP *"_ 1 ba
############
NB.   - tgevcul:
NB.       berr0 := max(||(β(i)*S - α(i)*P)^H * l(i)|| / (ulp * max(||β(i)*S||,||α(i)*P||)))
NB.       berr1 := (max(||l(i)||) - 1) / (ulp * n)

   vberr=. ((vberr0>.vberr1)~(2&{.))~                                  berr=.  SPxxxx vberr  Y
   vberr1=. (normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@]                 berr1=. xxxx vberr1 Y
   vberr0=. normir@:(((vnom%vden)~vbSaP)~)                             berr0=. SP vberr0 Y
   vbSaP=. *"_ 1|:@|.@:(diag"2)                                        bSaP=. vbSaP SP
   vden=. (FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)            den=. bSaP vden Y
   vnom=. norm1r@:((mp"1 2(-/"3))~ct)                                  nom=. bSaP vnom Y

(((normir@:((((norm1r@:((mp"1 2(-/"3))~ct))%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~(*"_ 1|:@|.@:(diag"2)))~))>.((normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@]))~(2&{.))~
############
NB.   - tgevcur:
NB.       berr0 := max(||(β(i)*S - α(i)*P)   * r(i)|| / (ulp * max(||β(i)*S||,||α(i)*P||)))
NB.       berr1 := (max(||r(i)||) - 1) / (ulp * n)

   vberr=. ((vberr0>.vberr1)~(2&{.))~                                  berr=.  SPxxxx vberr  Y
   vberr1=. (normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@]                 berr1=. xxxx vberr1 X
   vberr0=. normir@:(((vnom%vden)~vbSaP)~)                             berr0=. SP vberr0 X
   vbSaP=. *"_ 1|:@|.@:(diag"2)                                        bSaP=. vbSaP SP
   vden=. (FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)            den=. bSaP vden X
   vnom=. norm1r@:((mp"2 1~(-/"3))~|:)                                 nom=. bSaP vnom X

(((normir@:((((norm1r@:((mp"2 1~(-/"3))~|:))%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~(*"_ 1|:@|.@:(diag"2)))~))>.((normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@]))~(2&{.))~
############
NB.   - tgevcub:

   vberr=. ((vberr0>.vberr1)~(2&{.))~                                  berr=.  SPxxxx vberr  Y
   vberr1=. (>./@:normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@]            berr1=. xxxx vberr1 YX
   vberr0=. >./@:normir@:(((vnom%vden)~vbSaP)~)                        berr0=. SP vberr0 YX
   vbSaP=. *"_ 1|:@|.@:(diag"2)                                        bSaP=. vbSaP SP
   vden=. (FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)            den=. bSaP vden YX
   vnom=. ((norm1r@:(mp"1 2~ct@{.)>.norm1r@:((mp"2 1)|:@{:))~(-/"3))~  nom=. bSaP vnom YX

(((>./@:normir@:((((((norm1r@:(mp"1 2~ct@{.)>.norm1r@:((mp"2 1)|:@{:))~(-/"3))~)%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~(*"_ 1|:@|.@:(diag"2)))~))>.((>./@:normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@]))~(2&{.))~
############
NB.   - tgevculb:
NB.       berr0 := max(||(β(i)*H - α(i)*T)^H * l(i)|| / (ulp * max(||β(i)*H||,||α(i)*T||)))
NB.       berr1 := (max(||lb(i)||) - 1) / (ulp * n)

the same!
############
testtgevc=: 3 : 0
  berrul=:  (((     normir@:((((  norm1r@:((mp"1 2 (-/"3))~ct   )                                           )%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~(*"_ 1|:@|.@:(diag"2)))~))>.((     normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@]))~(2&{.))~ f.
  berrur=:  (((     normir@:((((                                   norm1r@:((mp"2 1~(-/"3))~|:   )          )%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~(*"_ 1|:@|.@:(diag"2)))~))>.((     normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@]))~(2&{.))~ f.
  berrub=:  (((>./@:normir@:((((((norm1r@:( mp"1 2        ~ct@{.)>.norm1r@:((mp"2 1       ) |:@{:))~(-/"3))~)%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~(*"_ 1|:@|.@:(diag"2)))~))>.((>./@:normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@]))~(2&{.))~ f.

  n=. c y
  hs=. 0,n
  HTl=. hs gghrdlnn ((unmqllc~ ,: trl@]) geqlf)/ y           NB. FIXME
  HTu=. hs gghrdunn ((unmqrlc~ ,: tru@]) geqrf)/ y
  HTSPQZl=. HTl , hs hgezqsvv HTl , ,:~ idmat n
  HTSPQZu=. HTu , hs hgeqzsvv HTu , ,:~ idmat n
  SPl=. 2 3 { HTSPQZl
  SPu=. 2 3 { HTSPQZu
  rcondl=. <./ trlcon1"2 SPl
  rcondu=. <./ trucon1"2 SPu

  ('tgevcll'  tmonad (]`]          `(rcondl"_)`(_."_)`berrul)) SPl
  ('tgevclr'  tmonad (]`]          `(rcondl"_)`(_."_)`berrur)) SPl
  ('tgevclb'  tmonad (]`]          `(rcondl"_)`(_."_)`berrub)) SPl
  ('tgevcllb' tmonad (]`(2 3 4  &{)`(rcondl"_)`(_."_)`berrul)) HTSPQZl
  ('tgevclrb' tmonad (]`(2 3   5&{)`(rcondl"_)`(_."_)`berrur)) HTSPQZl
  ('tgevclbb' tmonad (]`(2 3 4 5&{)`(rcondl"_)`(_."_)`berrub)) HTSPQZl
  ('tgevcul'  tmonad (]`]          `(rcondu"_)`(_."_)`berrul)) SPu
  ('tgevcur'  tmonad (]`]          `(rcondu"_)`(_."_)`berrur)) SPu
  ('tgevcub'  tmonad (]`]          `(rcondu"_)`(_."_)`berrub)) SPu
  ('tgevculb' tmonad (]`(2 3 4  &{)`(rcondu"_)`(_."_)`berrul)) HTSPQZu
  ('tgevcurb' tmonad (]`(2 3   5&{)`(rcondu"_)`(_."_)`berrur)) HTSPQZu
  ('tgevcubb' tmonad (]`(2 3 4 5&{)`(rcondu"_)`(_."_)`berrub)) HTSPQZu

  erase 'berrll berrlr berrlb berrul berrur berrub'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testevc
NB.
NB. Description:
NB.   Adv. to make verb to test tgevcxxx by matrices of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testevc
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testevc_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testevc_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testevc_mt_ 150 150

testevc=: 1 : 'EMPTY_mt_ [ (testtgevc_mt_ @ u @ (2&,)) ^: (=/)'
