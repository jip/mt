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
NB.   Calculate initial parameters for tgevcly and tgevclx
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
NB.   Compute some or all of left eigenvectors
NB.
NB. Syntax:
NB.   Yns=. (ios ; init) tgevcly SP
NB. where
NB.   ios   - k-vector, lIOS eigenvectors to compute
NB.   init  - boxed 8-vector, the output of tgevci
NB.   SP    - 2×n×n-matrix (S,:P), generalized Schur form,
NB.           produced by hgezqsxx
NB.   Yns   - k×n-matrix, some or all of left eigenvectors,
NB.           non-scaled
NB.   k     - integer in range [0,n]

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
      NB.   y * (a*A - b*B) = 0  (rowwise)
      NB. work[0:j-1] contains sums w
      NB. work[j+1:je] contains y
      work=. 1 ,~ (-/) ((je { ios) { abcoeff) * (((0,],0:),:(2 1,])) je { ios) ({.@(1 0 2&|:)) ;. 0 y
      di=. je { d
      j=. <: je { ios
      while. j >: 0 do.
        NB. form:
        NB.   y[j] = - w[j] / di
        NB. with scaling and perturbation of the denominator
        abs1wj=. sorim j { work
        if. *.`<:/ (j { cond2) , abs1wj do.
          work=. work % abs1wj
        end.
        work=. j (-@(%&(j{di))) upd work
        abs1wj=. sorim j { work
        if. j > 0 do.
          NB. w = w + y[j] * (a*S[:,j] - b*P[:,j]) with scaling
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

NB. ---------------------------------------------------------
NB. tgevclx
NB.
NB. Description:
NB.   Compute some or all of right eigenvectors
NB.
NB. Syntax:
NB.   Xns=. (ios ; init) tgevclx SP
NB. where
NB.   ios   - k-vector, lIOS eigenvectors to compute
NB.   init  - boxed 8-vector, the output of tgevci
NB.   SP    - 2×n×n-matrix (S,:P), generalized Schur form,
NB.           produced by hgezqsxx
NB.   Xns   - k×n-matrix, some or all of right eigenvectors,
NB.           non-scaled
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
      NB.   x * (a*A - b*B)^H = 0 ,
      NB. colwise in (a*A - b*B)^H , or rowwise in (a*A - b*B)
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

NB. ---------------------------------------------------------
NB. tgevcs
NB.
NB. Description:
NB.   Scale eigenvectors
NB.
NB. Syntax:
NB.   V=. tgevcs W
NB. where
NB.   W - k×n-matrix, some or all of left or right
NB.       eigenvectors, non-scaled
NB.   V - k×n-matrix, scaled W

tgevcs=: 3 : 0
  norm=. normitr y
  ios=. (#y) #"0 FP_SFMIN < norm
  y=. y % norm
  y=. ios } 0 ,: y
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. tgevcll
NB. tgevclr
NB. tgevclb
NB.
NB. Description:
NB.   Compute some or all of left eigenvectors Y and/or right
NB.   eigenvectors X:
NB.     E2 * Y * S = E1 * Y * P
NB.     S * X^H * E2 = P * X^H * E1
NB.   corresponding to eigenvalues represented as a pair of
NB.   values (α[i],β[i]):
NB.     E1 = diagmat(α[0:n-1]) = diagmat(diag(S))
NB.     E2 = diagmat(β[0:n-1]) = diagmat(diag(P))
NB.   for matrix pair (S,P) of lower triangular matrices
NB.   produced by the generalized Schur factorization
NB.   hgezqsxx:
NB.     Q^H * S * Z = A
NB.     Q^H * P * Z = B
NB.
NB. Syntax:
NB.   Y=.     [ios] tgevcll SP
NB.   X=.     [ios] tgevclr SP
NB.   'Y X'=. [ios] tgevclb SP
NB. where
NB.   ios - k-vector, lIOS eigenvectors to compute, default
NB.         is "all eigenvectors"
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgezqsxx
NB.   Y   - k×n-matrix, some or all of left eigenvectors
NB.   X   - k×n-matrix, some or all of right eigenvectors
NB.
NB. Assertions:
NB.   (tgevcll -: tgevcur &.: (ct"2)) SP
NB.   (tgevclr -: tgevcul &.: (ct"2)) SP
NB.   (tgevclb -: tgevcub &.: (ct"2)) SP

tgevcll=:  ($:~ i.@c) : (([;tgevci) tgevcs  @ tgevcly           ])
tgevclr=:  ($:~ i.@c) : (([;tgevci) tgevcs  @          tgevclx  ])
tgevclb=:  ($:~ i.@c) : (([;tgevci) tgevcs"2@(tgevcly,:tgevclx) ])

NB. ---------------------------------------------------------
NB. tgevcul
NB. tgevcur
NB. tgevcub
NB.
NB. Description:
NB.   Compute some or all of left eigenvectors Y and/or right
NB.   eigenvectors X:
NB.     E2 * Y^H * S = E1 * Y^H * P
NB.     S * X * E2 = P * X * E1
NB.   corresponding to eigenvalues represented as a pair of
NB.   values (α[i],β[i]):
NB.     E1 = diagmat(α[0:n-1]) = diagmat(diag(S))
NB.     E2 = diagmat(β[0:n-1]) = diagmat(diag(P))
NB.   for matrix pair (S,P) of upper triangular matrices
NB.   produced by the generalized Schur factorization
NB.   hgeqzsxx
NB.     Q * S * Z^H = A
NB.     Q * P * Z^H = B
NB.
NB. Syntax:
NB.   Y=.     [ios] tgevcul SP
NB.   X=.     [ios] tgevcur SP
NB.   'Y X'=. [ios] tgevcub SP
NB. where
NB.   ios - k-vector, lIOS eigenvectors to compute, default
NB.         is "all eigenvectors"
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgeqzsxx
NB.   Y   - n×k-matrix, some or all of left eigenvectors
NB.   X   - n×k-matrix, some or all of right eigenvectors
NB.
NB. Assertions:
NB.   (tgevcul -: tgevclr &.: (ct"2)) SP
NB.   (tgevcur -: tgevcll &.: (ct"2)) SP
NB.   (tgevcub -: tgevclb &.: (ct"2)) SP
NB.
NB. Notes:
NB. - tgevcul implements LAPACK's xTGEVC('L','S')
NB. - tgevcur implements LAPACK's xTGEVC('R','S')
NB. - tgevcub implements LAPACK's xTGEVC('B','S')

tgevcul=: ($:~ i.@c) : ((([;tgevci) tgevcs  @ tgevclx           ]) &.: (|:"2))
tgevcur=: ($:~ i.@c) : ((([;tgevci) tgevcs  @          tgevcly  ]) &.: (|:"2))
tgevcub=: ($:~ i.@c) : ((([;tgevci) tgevcs"2@(tgevclx,:tgevcly) ]) &.: (|:"2))

NB. ---------------------------------------------------------
NB. tgevcllb
NB. tgevclrb
NB. tgevclbb
NB.
NB. Description:
NB.   Compute left eigenvectors Y*Q and/or right eigenvectors
NB.   X*Z:
NB.     E2 * (Y * Q) * S = E1 * (Y * Q) * P
NB.     S * (X * Z)^H * E2 = P * (X * Z)^H * E1
NB.   corresponding to eigenvalues represented as a pair of
NB.   values (α[i],β[i]):
NB.     E1 = diagmat(α[0:n-1]) = diagmat(diag(S))
NB.     E2 = diagmat(β[0:n-1]) = diagmat(diag(P))
NB.   for matrix pair (S,P) of lower triangular matrices
NB.   produced by the generalized Schur factorization
NB.   hgezqsxx:
NB.     Q^H * S * Z = A
NB.     Q^H * P * Z = B
NB.
NB. Syntax:
NB.   YQ=.      tgevcllb SP , Q
NB.   XZ=.      tgevclrb SP , Z
NB.   'YQ XZ'=. tgevclbb SP , Q ,: Z
NB. where
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgezqsxx
NB.   YQ  - n×n-matrix, left eigenvectors Y*Q
NB.   XZ  - n×n-matrix, right eigenvectors X*Z
NB.
NB. Assertions:
NB.   (tgevcllb -: tgevcurb &.: (ct"2)) SP , Q
NB.   (tgevclrb -: tgevculb &.: (ct"2)) SP , Z
NB.   (tgevclbb -: tgevcubb &.: (ct"2)) SP , Q ,: Z

tgevcllb=:  (i.@c ([;tgevci) 2&{.) tgevcs  @( tgevcly mp 2{]                    ) ]
tgevclrb=:  (i.@c ([;tgevci) 2&{.) tgevcs  @(                   tgevclx mp _1{] ) ]
tgevclbb=:  (i.@c ([;tgevci) 2&{.) tgevcs"2@((tgevcly mp 2{]),:(tgevclx mp _1{])) ]

NB. ---------------------------------------------------------
NB. tgevculb
NB. tgevcurb
NB. tgevcubb
NB.
NB. Description:
NB.   Compute left eigenvectors Q*Y and/or right eigenvectors
NB.   Z*X:
NB.     E2 * (Q * Y)^H * S = E1 * (Q * Y)^H * P
NB.     S * (Z * X) * E2 = P * (Z * X) * E1
NB.   corresponding to eigenvalues represented as a pair of
NB.   values (α[i],β[i]):
NB.     E1 = diagmat(α[0:n-1]) = diagmat(diag(S))
NB.     E2 = diagmat(β[0:n-1]) = diagmat(diag(P))
NB.   for matrix pair (S,P) of upper triangular matrices
NB.   produced by the generalized Schur factorization
NB.   hgeqzsxx
NB.     Q * S * Z^H = A
NB.     Q * P * Z^H = B
NB.
NB. Syntax:
NB.   QY=.      tgevculb SP , Q
NB.   ZX=.      tgevcurb SP , Z
NB.   'QY ZX'=. tgevcubb SP , Q ,: Z
NB. where
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgeqzsxx
NB.   QY  - n×n-matrix, left eigenvectors Q*Y
NB.   ZX  - n×n-matrix, right eigenvectors Z*X
NB.
NB. Assertions:
NB.   (tgevculb -: tgevclrb &.: (ct"2)) SP , Q
NB.   (tgevcurb -: tgevcllb &.: (ct"2)) SP , Z
NB.   (tgevcubb -: tgevclbb &.: (ct"2)) SP , Q ,: Z
NB.
NB. Notes:
NB. - tgevculb implements LAPACK's xTGEVC('L','B')
NB. - tgevcurb implements LAPACK's xTGEVC('R','B')
NB. - tgevcubb implements LAPACK's xTGEVC('B','B')

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

testtgevc=: 3 : 0
  berrll=:  (     normir@:((((                                   norm1r@:((mp"2 1~(-/"3))~     )          )%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~((_2 _1&{)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.((     normir@:<:@:normitr%(FP_EPS*FP_BASE)*c)@])
  berrlr=:  (     normir@:((((  norm1r@:((mp"1 2 (-/"3))~+    )                                           )%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~((_2 _1&{)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.((     normir@:<:@:normitr%(FP_EPS*FP_BASE)*c)@])
  berrlb=:  (>./@:normir@:((((((norm1r@:( mp"1 2        ~+ @{.)>.norm1r@:((mp"2 1       )    {:))~(-/"3))~)%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~((_2 _1&{)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.((>./@:normir@:<:@:normitr%(FP_EPS*FP_BASE)*c)@])
  berrul=:  (     normir@:((((  norm1r@:((mp"1 2 (-/"3))~ct   )                                           )%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~((_2 _1&{)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.((     normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@])
  berrur=:  (     normir@:((((                                   norm1r@:((mp"2 1~(-/"3))~|:   )          )%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~((_2 _1&{)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.((     normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@])
  berrub=:  (>./@:normir@:((((((norm1r@:( mp"1 2        ~ct@{.)>.norm1r@:((mp"2 1       ) |:@{:))~(-/"3))~)%((FP_EPS*FP_BASE)*(FP_SFMIN>.(>./"1)@:(norm1"2)@[)))~((_2 _1&{)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.((>./@:normir@:<:@:normitc%(FP_EPS*FP_BASE)*c)@])

  rcondl=. <./ trlcon1"2 SPl=. 2 {. SPQZHTl=. (([ (((0,[) hgezqsvv (, ,:~@idmat)~) , ]) ((gghrdlnn~ (0&,))~((unmlqrc~ ,: trl@]) gelqf)/))~ c) y
  rcondu=. <./ trucon1"2 SPu=. 2 {. SPQZHTu=. (([ (((0,[) hgeqzsvv (, ,:~@idmat)~) , ]) ((gghrdunn~ (0&,))~((unmqrlc~ ,: tru@]) geqrf)/))~ c) y

  ('tgevcll'  tmonad (]          `]`(rcondl"_)`(_."_)`berrll)) SPl
  ('tgevclr'  tmonad (]          `]`(rcondl"_)`(_."_)`berrlr)) SPl
  ('tgevclb'  tmonad (]          `]`(rcondl"_)`(_."_)`berrlb)) SPl
  ('tgevcllb' tmonad ((0 1 2  &{)`]`(rcondl"_)`(_."_)`berrll)) SPQZHTl
  ('tgevclrb' tmonad ((0 1   3&{)`]`(rcondl"_)`(_."_)`berrlr)) SPQZHTl
  ('tgevclbb' tmonad ((0 1 2 3&{)`]`(rcondl"_)`(_."_)`berrlb)) SPQZHTl
  ('tgevcul'  tmonad (]          `]`(rcondu"_)`(_."_)`berrul)) SPu
  ('tgevcur'  tmonad (]          `]`(rcondu"_)`(_."_)`berrur)) SPu
  ('tgevcub'  tmonad (]          `]`(rcondu"_)`(_."_)`berrub)) SPu
  ('tgevculb' tmonad ((0 1 2  &{)`]`(rcondu"_)`(_."_)`berrul)) SPQZHTu
  ('tgevcurb' tmonad ((0 1   3&{)`]`(rcondu"_)`(_."_)`berrur)) SPQZHTu
  ('tgevcubb' tmonad ((0 1 2 3&{)`]`(rcondu"_)`(_."_)`berrub)) SPQZHTu

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
