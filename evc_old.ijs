NB. Eigenvectors
NB.
NB. tgevcxx     Some or all of the left and/or right
NB.             eigenvectors of generalized Schur form
NB. tgevcxxb    Backtransformed left and/or right
NB.             eigenvectors of generalized Schur form
NB.
NB. testtgevc_old   Test tgevcxxx by general matrices given
NB. testevc_old     Adv. to make verb to test tgevcxxx by
NB.             matrices of generator and shape given
NB.
NB. Version: 0.7.0 2011-08-06
NB.
NB. Copyright 2011 Igor Zhuravlov
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
NB. tgevci_old
NB.
NB. Description:
NB.   Calculate initial parameters for tgevcly_old and tgevclx_old
NB.
NB. Syntax:
NB.   'bignum d2 abrwork cond1 cond2 abcoeff abcoeffa d'=. ios tgevci_old SP
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

tgevci_old=: 4 : 0
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
NB. tgevcly_old
NB.
NB. Description:
NB.   Compute some or all of non-scaled left eigenvectors
NB.
NB. Syntax:
NB.   W=. (ios ; init) tgevcly_old SP
NB. where
NB.   ios   - k-vector, lIOS eigenvectors to compute
NB.   init  - boxed 8-vector, the output of tgevci_old
NB.   SP    - 2×n×n-matrix (S,:P), generalized Schur form,
NB.           produced by hgezqsxx
NB.   W     - k×n-matrix, some or all of left eigenvectors,
NB.           non-scaled
NB.   k     - integer in range [0,n]

tgevcly_old=: 4 : 0
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
          if. ((abcoeffa mp&(j&{) abrwork) >: (bignum % abs1wj)) *. (1 < abs1wj) do.
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
  W
)

NB. ---------------------------------------------------------
NB. tgevclx_old
NB.
NB. Description:
NB.   Compute some or all of non-scaled right eigenvectors
NB.
NB. Syntax:
NB.   W=. (ios ; init) tgevclx_old SP
NB. where
NB.   ios   - k-vector, lIOS eigenvectors to compute
NB.   init  - boxed 8-vector, the output of tgevci_old
NB.   SP    - 2×n×n-matrix (S,:P), generalized Schur form,
NB.           produced by hgezqsxx
NB.   W     - k×n-matrix, some or all of right eigenvectors,
NB.           non-scaled
NB.   k     - integer in range [0,n]

tgevclx_old=: 4 : 0
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
      NB. columnwise in (a*A - b*B)^H , or rowwise in
      NB. (a*A - b*B)
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
  W
)

NB. ---------------------------------------------------------
NB. tgevcs_old
NB.
NB. Description:
NB.   Scale left or right eigenvectors
NB.
NB. Syntax:
NB.   V=. tgevcs_old W
NB. where
NB.   W - k×n-matrix, some or all of left or right
NB.       eigenvectors, non-scaled
NB.   V - k×n-matrix, scaled W

tgevcs_old=: 3 : 0
  norm=. normitr y
  ios=. (#y) #"0 FP_SFMIN < norm
  y=. y % norm
  y=. ios } 0 ,: y
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. tgevcll_old
NB. tgevclr_old
NB. tgevclb_old
NB.
NB. Description:
NB.   Compute some or all of left eigenvectors Y:
NB.     E2 * Y * S = E1 * Y * P
NB.   and/or right eigenvectors X:
NB.     S * X^H * E2 = P * X^H * E1
NB.   for matrix pair (S,P) of lower triangular matrices
NB.   produced by the generalized Schur factorization
NB.   hgezqsxx:
NB.     Q^H * S * Z = A
NB.     Q^H * P * Z = B
NB.   Each i-th eigenvector (row) from Y and X has a
NB.   corresponding eigenvalue represented as a pair of i-th
NB.   diagonal elements in matrices E1, E2:
NB.     E1=. diagmat(diag(S))
NB.     E2=. diagmat(diag(P))
NB.
NB. Syntax:
NB.   Y=.     [ios] tgevcll_old SP
NB.   X=.     [ios] tgevclr_old SP
NB.   'Y X'=. [ios] tgevclb_old SP
NB. where
NB.   ios - k-vector, optional lIOS eigenvectors to compute,
NB.         default is "all eigenvectors"
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgezqsxx
NB.   Y   - k×n-matrix, some or all of left eigenvectors
NB.   X   - k×n-matrix, some or all of right eigenvectors
NB.   k   - integer in range [0,n]
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (tgevcll_old -:      tgevcur_old&.:(|:"2)) SP
NB.   (tgevclr_old -:      tgevcul_old&.:(|:"2)) SP
NB.   (tgevclb_old -: 1 A. tgevcub_old&.:(|:"2)) SP
NB.   (E2 mp Y mp S) -: (E1 mp Y mp P)
NB.   (S mp (ct X) mp E2) -: (P mp (ct X) mp E1)
NB. where
NB.   'Y X'=. tgevclb_old SP
NB.   'E1 E2'=. diagmat@diag"2 SP

tgevcll_old=:  ($:~ i.@c) : (([;tgevci_old) tgevcs_old  @ tgevcly_old           ])
tgevclr_old=:  ($:~ i.@c) : (([;tgevci_old) tgevcs_old  @          tgevclx_old  ])
tgevclb_old=:  ($:~ i.@c) : (([;tgevci_old) tgevcs_old"2@(tgevcly_old,:tgevclx_old) ])

NB. ---------------------------------------------------------
NB. tgevcul_old
NB. tgevcur_old
NB. tgevcub_old
NB.
NB. Description:
NB.   Compute some or all of left eigenvectors Y:
NB.     E2 * Y^H * S = E1 * Y^H * P
NB.   and/or right eigenvectors X:
NB.     S * X * E2 = P * X * E1
NB.   for matrix pair (S,P) of upper triangular matrices
NB.   produced by the generalized Schur factorization
NB.   hgeqzsxx:
NB.     Q * S * Z^H = A
NB.     Q * P * Z^H = B
NB.   Each i-th eigenvector (column) from Y and X has a
NB.   corresponding eigenvalue represented as a pair of i-th
NB.   diagonal elements in matrices E1, E2:
NB.     E1=. diagmat(diag(S))
NB.     E2=. diagmat(diag(P))
NB.
NB. Syntax:
NB.   Y=.     [ios] tgevcul_old SP
NB.   X=.     [ios] tgevcur_old SP
NB.   'Y X'=. [ios] tgevcub_old SP
NB. where
NB.   ios - k-vector, optional lIOS eigenvectors to compute,
NB.         default is "all eigenvectors"
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgeqzsxx
NB.   Y   - n×k-matrix, some or all of left eigenvectors
NB.   X   - n×k-matrix, some or all of right eigenvectors
NB.   k   - integer in range [0,n]
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (tgevcul_old -:      tgevclr_old&.:(|:"2)) SP
NB.   (tgevcur_old -:      tgevcll_old&.:(|:"2)) SP
NB.   (tgevcub_old -: 1 A. tgevclb_old&.:(|:"2)) SP
NB.   (E2 mp (ct Y) mp S) -: (E1 mp (ct Y) mp P)
NB.   (S mp X mp E2) -: (P mp X mp E1)
NB. where
NB.   'Y X'=. tgevcub_old SP
NB.   'E1 E2'=. diagmat@diag"2 SP
NB.
NB. Notes:
NB. - tgevcul_old models LAPACK's xTGEVC('L','S')
NB. - tgevcur_old models LAPACK's xTGEVC('R','S')
NB. - tgevcub_old models LAPACK's xTGEVC('B','S')

tgevcul_old=: ($:~ i.@c) : ((([;tgevci_old) tgevcs_old  @ tgevclx_old           ])&.:(|:"2))
tgevcur_old=: ($:~ i.@c) : ((([;tgevci_old) tgevcs_old  @          tgevcly_old  ])&.:(|:"2))
tgevcub_old=: ($:~ i.@c) : ((([;tgevci_old) tgevcs_old"2@(tgevclx_old,:tgevcly_old) ])&.:(|:"2))

NB. ---------------------------------------------------------
NB. tgevcllb_old
NB. tgevclrb_old
NB. tgevclbb_old
NB.
NB. Description:
NB.   Compute left eigenvectors Y*Q:
NB.     E2 * (Y * Q) * A = E1 * (Y * Q) * B
NB.   and/or right eigenvectors X*Z:
NB.     A * (X * Z)^H * E2 = B * (X * Z)^H * E1
NB.   for matrix pair (A,B) and matrices Q, Z produced by the
NB.   generalized Schur factorization hgezqsxx:
NB.     Q^H * S * Z = A
NB.     Q^H * P * Z = B
NB.   Each i-th eigenvector (row) from Y*Q and X*Z has a
NB.   corresponding eigenvalue represented as a pair of i-th
NB.   diagonal elements in matrices E1, E2:
NB.     E1=. diagmat(diag(S))
NB.     E2=. diagmat(diag(P))
NB.
NB. Syntax:
NB.   YQ=.      tgevcllb_old SP , Q
NB.   XZ=.      tgevclrb_old SP , Z
NB.   'YQ XZ'=. tgevclbb_old SP , Q ,: Z
NB. where
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgezqsxx
NB.   YQ  - n×n-matrix, left eigenvectors Y*Q
NB.   XZ  - n×n-matrix, right eigenvectors X*Z
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (tgevcllb_old -: tgevcurb_old&.:(|:"2        )) SP , Q2
NB.   (tgevclrb_old -: tgevculb_old&.:(|:"2        )) SP , Z2
NB.   (tgevclbb_old -: tgevcubb_old&.:(|:"2@:(1&A.))) SP , Q2 ,: Z2
NB.   D                         -: B mp Z0
NB.   D                         -: BZ0f unmlqrn B
NB.   A                         -: BZ0f unmlqrc C
NB.   Q1                        -: dQ0 mp Q0
NB.   Q1                        -: dQ0
NB.   Z1                        -: dZ0 mp Z0
NB.   Q2                        -: dQ1 mp Q1
NB.   Z2                        -: dZ1 mp Z1
NB.   dQ1dQ0                    -: dQ1 mp dQ0
NB.   dZ1dZ0                    -: dZ1 mp dZ0
NB.   (E2 mp Y mp S)            -: (E1 mp Y mp P)
NB.   (S mp (ct X) mp E2)       -: (P mp (ct X) mp E1)
NB.   (E2 mp YdQ1 mp H)         -: (E1 mp YdQ1 mp T)
NB.   (H mp (ct XdZ1) mp E2)    -: (T mp (ct XdZ1) mp E1)
NB.   (E2 mp YdQ1dQ0 mp A)      -: (E1 mp YdQ1dQ0 mp B)
NB.   (A mp (ct XdZ1dZ0) mp E2) -: (B mp (ct XdZ1dZ0) mp E1)
NB.   (E2 mp YQ2 mp C)          -: (E1 mp YQ2 mp D)
NB.   (C mp (ct XZ2) mp E2)     -: (D mp (ct XZ2) mp E1)
NB. where
NB.   C - n×n-matrix, general
NB.   D - n×n-matrix, general
NB.   n=. # C
NB.   hs=. 0 , n
NB.   I=. idmat n
NB.   BZ0f=. gelqf D
NB.   B=. trl }:"1 BZ0f
NB.   Q0=. I
NB.   Z0=. unglq BZ0f
NB.   A=. C mp ct Z0
NB.   'H T Q1 Z1'=. hs gghrdlvv_old A , B , Q0 ,: Z0
NB.   'H T dQ0 dZ0'=. hs gghrdlvv_old A , B , ,:~ I
NB.   'S P Q2 Z2'=. hs hgezqsvv_old H , T , Q1 ,: Z1
NB.   'S P dQ1 dZ1'=. hs hgezqsvv_old H , T , ,:~ I
NB.   'S P dQ1dQ0 dZ1dZ0'=. hs hgezqsvv_old H , T , dQ0 ,: dZ0
NB.   'Y X'=. tgevclb_old S ,: P
NB.   'YdQ1 XdZ1'=. tgevclbb_old S , P , dQ1 ,: dZ1
NB.   'YdQ1dQ0 XdZ1dZ0'=. tgevclbb_old S , P , dQ1dQ0 ,: dZ1dZ0
NB.   'YQ2 XZ2'=. tgevclbb_old S , P , Q2 ,: Z2
NB.   'E1 E2'=. diagmat@diag"2 S ,: P

tgevcllb_old=:  (i.@c ([;tgevci_old) 2&{.) tgevcs_old  @( tgevcly_old mp 2{]                    ) ]
tgevclrb_old=:  (i.@c ([;tgevci_old) 2&{.) tgevcs_old  @(                   tgevclx_old mp _1{] ) ]
tgevclbb_old=:  (i.@c ([;tgevci_old) 2&{.) tgevcs_old"2@((tgevcly_old mp 2{]),:(tgevclx_old mp _1{])) ]

NB. ---------------------------------------------------------
NB. tgevculb_old
NB. tgevcurb_old
NB. tgevcubb_old
NB.
NB. Description:
NB.   Compute left eigenvectors Q*Y:
NB.     E2 * (Q * Y)^H * A = E1 * (Q * Y)^H * B
NB.   and/or right eigenvectors Z*X:
NB.     A * (Z * X) * E2 = B * (Z * X) * E1
NB.   for matrix pair (A,B) and matrices Q, Z produced by the
NB.   generalized Schur factorization hgeqzsxx:
NB.     Q * S * Z^H = A
NB.     Q * P * Z^H = B
NB.   Each i-th eigenvector (column) from Q*Y and Z*X has a
NB.   corresponding eigenvalue represented as a pair of i-th
NB.   diagonal elements in matrices E1, E2:
NB.     E1=. diagmat(diag(S))
NB.     E2=. diagmat(diag(P))
NB.
NB. Syntax:
NB.   QY=.      tgevculb_old SP , Q
NB.   ZX=.      tgevcurb_old SP , Z
NB.   'QY ZX'=. tgevcubb_old SP , Q ,: Z
NB. where
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgeqzsxx
NB.   QY  - n×n-matrix, left eigenvectors Q*Y
NB.   ZX  - n×n-matrix, right eigenvectors Z*X
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (tgevculb_old -: tgevclrb_old&.:(|:"2        )) SP , Q2
NB.   (tgevcurb_old -: tgevcllb_old&.:(|:"2        )) SP , Z2
NB.   (tgevcubb_old -: tgevclbb_old&.:(|:"2@:(1&A.))) SP , Q2 ,: Z2
NB.   D                         -: Q0 mp B
NB.   D                         -: Q0fB unmqrln B
NB.   A                         -: Q0fB unmqrlc C
NB.   Q1                        -: Q0 mp dQ0
NB.   Z1                        -: Z0 mp dZ0
NB.   Z1                        -: dZ0
NB.   Q2                        -: Q1 mp dQ1
NB.   Z2                        -: Z1 mp dZ1
NB.   dQ0dQ1                    -: dQ0 mp dQ1
NB.   dZ0dZ1                    -: dZ0 mp dZ1
NB.   (E2 mp (ct Y) mp S)       -: (E1 mp (ct Y) mp P)
NB.   (S mp X mp E2)            -: (P mp X mp E1)
NB.   (E2 mp (ct dQ1Y) mp H)    -: (E1 mp (ct dQ1Y) mp T)
NB.   (H mp dZ1X mp E2)         -: (T mp dZ1X mp E1)
NB.   (E2 mp (ct dQ0dQ1Y) mp A) -: (E1 mp (ct dQ0dQ1Y) mp B)
NB.   (A mp dZ0dZ1X mp E2)      -: (B mp dZ0dZ1X mp E1)
NB.   (E2 mp (ct Q2Y) mp C)     -: (E1 mp (ct Q2Y) mp D)
NB.   (C mp Z2X mp E2)          -: (D mp Z2X mp E1)
NB. where
NB.   C - n×n-matrix, general
NB.   D - n×n-matrix, general
NB.   n=. # C
NB.   hs=. 0 , n
NB.   I=. idmat n
NB.   Q0fB=. geqrf D
NB.   B=. tru }: Q0fB
NB.   Q0=. ungqr Q0fB
NB.   Z0=. I
NB.   A=. Q0 (mp~ ct)~ C
NB.   'H T Q1 Z1'=. hs gghrduvv_old A , B , Q0 ,: Z0
NB.   'H T dQ0 dZ0'=. hs gghrduvv_old A , B , ,:~ I
NB.   'S P Q2 Z2'=. hs hgeqzsvv_old H , T , Q1 ,: Z1
NB.   'S P dQ1 dZ1'=. hs hgeqzsvv_old H , T , ,:~ I
NB.   'S P dQ0dQ1 dZ0dZ1'=. hs hgeqzsvv_old H , T , dQ0 ,: dZ0
NB.   'Y X'=. tgevcub_old S ,: P
NB.   'dQ1Y dZ1X'=. tgevcubb_old S , P , dQ1 ,: dZ1
NB.   'dQ0dQ1Y dZ0dZ1X'=. tgevcubb_old S , P , dQ0dQ1 ,: dZ0dZ1
NB.   'Q2Y Z2X'=. tgevcubb_old S , P , Q2 ,: Z2
NB.   'E1 E2'=. diagmat@diag"2 S ,: P
NB.
NB. Notes:
NB. - tgevculb_old models LAPACK's xTGEVC('L','B')
NB. - tgevcurb_old models LAPACK's xTGEVC('R','B')
NB. - tgevcubb_old models LAPACK's xTGEVC('B','B')

tgevculb_old=: ((i.@c ([;tgevci_old) 2&{.) tgevcs_old  @( tgevclx_old mp 2{]                    ) ])&.:(|:"2)
tgevcurb_old=: ((i.@c ([;tgevci_old) 2&{.) tgevcs_old  @(                   tgevcly_old mp _1{] ) ])&.:(|:"2)
tgevcubb_old=: ((i.@c ([;tgevci_old) 2&{.) tgevcs_old"2@((tgevclx_old mp 2{]),:(tgevcly_old mp _1{])) ])&.:(|:"2)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testtgevc_old
NB.
NB. Description:
NB.   Test tgevcxxx by general matrices given
NB.
NB. Syntax:
NB.   testtgevc_old AB
NB. where
NB.   AB - 2×n×n-report (A,:B)
NB.
NB. Formula:
NB.   berr := max(berr0,berr1)
NB. where
NB.   ||M|| := max(||M||_1 , FP_SFMIN)
NB.   ||v|| := max(|Re(v(i))|+|Im(v(i))|)
NB.   e1(i) - i-th eigenvalue, also i-th element on S
NB.           diagonal
NB.   e2(i) - i-th eigenvalue, also i-th element on P
NB.           diagonal
NB.   l(i)  - i-th left eigenvector
NB.   lb(i) - i-th back transformed left eigenvector
NB.   r(i)  - i-th right eigenvector
NB.   rb(i) - i-th back transformed right eigenvector
NB.   - tgevcll_old:
NB.       berr0 := max(||l(i) * (e2(i)*S - e1(i)*P)  || / (FP_PREC * max(|| e2(i)*S   ||,|| e1(i)*P   ||)))
NB.       berr1 := max(| ||l(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevclr_old:
NB.       berr0 := max(||r(i) * (e2(i)*S - e1(i)*P)^H|| / (FP_PREC * max(||(e2(i)*S)^H||,||(e1(i)*P)^H||)))
NB.       berr1 := max(| ||r(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevclb_old:
NB.       berr0 := berr(tgevcll_old)
NB.       berr1 := berr(tgevclr_old)
NB.   - tgevcllb_old:
NB.       berr0 := max(||l(i) * (e2(i)*H - e1(i)*T)  || / (FP_PREC * max(|| e2(i)*H   ||,|| e1(i)*T   ||)))
NB.       berr1 := max(| ||lb(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevclrb_old:
NB.       berr0 := max(||r(i) * (e2(i)*H - e1(i)*T)^H|| / (FP_PREC * max(||(e2(i)*H)^H||,||(e1(i)*T)^H||)))
NB.       berr1 := max(| ||rb(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevclbb_old:
NB.       berr0 := berr(tgevcllb_old)
NB.       berr1 := berr(tgevclrb_old)
NB.   - tgevcul_old:
NB.       berr0 := max(||(e2(i)*S - e1(i)*P)^H * l(i)|| / (FP_PREC * max(||(e2(i)*S)^H||,||(e1(i)*P)^H||)))
NB.       berr1 := max(| ||l(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevcur_old:
NB.       berr0 := max(||(e2(i)*S - e1(i)*P)   * r(i)|| / (FP_PREC * max(|| e2(i)*S   ||,|| e1(i)*P   ||)))
NB.       berr1 := max(| ||r(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevcub_old:
NB.       berr0 := berr(tgevcul_old)
NB.       berr1 := berr(tgevcur_old)
NB.   - tgevculb_old:
NB.       berr0 := max(||(e2(i)*H - e1(i)*T)^H * l(i)|| / (FP_PREC * max(||(e2(i)*H)^H||,||(e1(i)*T)^H||)))
NB.       berr1 := max(| ||lb(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevcurb_old:
NB.       berr0 := max(||(e2(i)*H - e1(i)*T)   * r(i)|| / (FP_PREC * max(|| e2(i)*H   ||,|| e1(i)*T   ||)))
NB.       berr1 := max(| ||rb(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevcubb_old:
NB.       berr0 := berr(tgevculb_old)
NB.       berr1 := berr(tgevcurb_old)
NB.
NB. Notes:
NB. - berrxx are non-iterative and are require O(N^3) RAM

testtgevc_old=: 3 : 0
  vberrll=:  (     normir@:((((  norm1r@:((mp"1 2 (-/"3))~     )                                           )      % (FP_PREC*(FP_SFMIN>.    (>./"1)@:( norm1       "2)@[)))~((_2&{.)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.(     normir@:<:@:normitr%FP_PREC*c)@]
  vberrlr=:  (     normir@:((((                                   norm1r@:((mp"2 1~(-/"3))~+    )          )      % (FP_PREC*(FP_SFMIN>.    (>./"1)@:(       normi "2)@[)))~((_2&{.)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.(     normir@:<:@:normitr%FP_PREC*c)@]
  vberrlb=:  (>./@:normir@:((((((norm1r@:( mp"1 2        ~   {.),:norm1r@:((mp"2 1       ) + @{:))~(-/"3))~)(>./@:%)(FP_PREC*(FP_SFMIN>.|:@:(>./"2)@:((norm1,normi)"2)@[)))~((_2&{.)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.(>./@:normir@:<:@:normitr%FP_PREC*c)@]

  vberrul=:  (     normir@:((((  norm1r@:((mp"1 2 (-/"3))~ct   )                                           )      % (FP_PREC*(FP_SFMIN>.    (>./"1)@:( normi       "2)@[)))~((_2&{.)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.(     normir@:<:@:normitc%FP_PREC*c)@]
  vberrur=:  (     normir@:((((                                   norm1r@:((mp"2 1~(-/"3))~|:   )          )      % (FP_PREC*(FP_SFMIN>.    (>./"1)@:(       norm1 "2)@[)))~((_2&{.)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.(     normir@:<:@:normitc%FP_PREC*c)@]
  vberrub=:  (>./@:normir@:((((((norm1r@:( mp"1 2        ~ct@{.),:norm1r@:((mp"2 1       ) |:@{:))~(-/"3))~)(>./@:%)(FP_PREC*(FP_SFMIN>.|:@:(>./"2)@:((normi,norm1)"2)@[)))~((_2&{.)*"_ 1|:@|.@:(diag"2)@(2&{.)))~))>.(>./@:normir@:<:@:normitc%FP_PREC*c)@]

  rcondl=. <./ trlcon1"2 SPl=. 2 {. SPQZHTl=. (([ (((0,[) hgezqsvv_old (, ,:~@idmat)~) , ]) ((gghrdlnn_old~ 0&,)~((unmlqrc~ ,: trl@:(}:"1)@]) gelqf)/))~ c) y
  rcondu=. <./ trucon1"2 SPu=. 2 {. SPQZHTu=. (([ (((0,[) hgeqzsvv_old (, ,:~@idmat)~) , ]) ((gghrdunn_old~ 0&,)~((unmqrlc~ ,: tru@  }:   @]) geqrf)/))~ c) y

  ('tgevcll_old'  tmonad (]          `]`(rcondl"_)`(_."_)`vberrll)) SPl
  ('tgevclr_old'  tmonad (]          `]`(rcondl"_)`(_."_)`vberrlr)) SPl
  ('tgevclb_old'  tmonad (]          `]`(rcondl"_)`(_."_)`vberrlb)) SPl
  ('tgevcllb_old' tmonad ((0 1 2  &{)`]`(rcondl"_)`(_."_)`vberrll)) SPQZHTl
  ('tgevclrb_old' tmonad ((0 1   3&{)`]`(rcondl"_)`(_."_)`vberrlr)) SPQZHTl
  ('tgevclbb_old' tmonad ((0 1 2 3&{)`]`(rcondl"_)`(_."_)`vberrlb)) SPQZHTl
  ('tgevcul_old'  tmonad (]          `]`(rcondu"_)`(_."_)`vberrul)) SPu
  ('tgevcur_old'  tmonad (]          `]`(rcondu"_)`(_."_)`vberrur)) SPu
  ('tgevcub_old'  tmonad (]          `]`(rcondu"_)`(_."_)`vberrub)) SPu
  ('tgevculb_old' tmonad ((0 1 2  &{)`]`(rcondu"_)`(_."_)`vberrul)) SPQZHTu
  ('tgevcurb_old' tmonad ((0 1   3&{)`]`(rcondu"_)`(_."_)`vberrur)) SPQZHTu
  ('tgevcubb_old' tmonad ((0 1 2 3&{)`]`(rcondu"_)`(_."_)`vberrub)) SPQZHTu

  erase 'vberrll vberrlr vberrlb vberrul vberrur vberrub'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testevc_old
NB.
NB. Description:
NB.   Adv. to make verb to test tgevcxxx by matrices of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testevc_old
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
NB.     _1 1 0 4 _6 4&gemat_mt_ testevc_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testevc_mt_ 150 150

testevc_old=: 1 : 'EMPTY_mt_ [ testtgevc_old_mt_@u@(2&,)^:(=/)'
