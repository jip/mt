NB. Eigenvectors
NB.
NB. tgevcxx     Some or all of the left and/or right
NB.             eigenvectors of generalized Schur form
NB. tgevcxxb    Backtransformed left and/or right
NB.             eigenvectors of generalized Schur form
NB.
NB. testtgevc   Test tgevcxxx by square matrices
NB. testevc     Adv. to make verb to test tgevcxxx by
NB.             matrices of generator and shape given
NB.
NB. Version: 0.11.0 2021-01-17
NB.
NB. Copyright 2011-2021 Igor Zhuravlov
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
NB.   'bignum d2 abrwork cond1 cond2 abcoeff abcoeffa d'=. iso tgevci SP
NB. where
NB.   SP       - 2×n×n-matrix (S,:P), generalized Schur form,
NB.              produced by hgezqsxx
NB.   iso      - k-vector, lISO eigenvectors to compute
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
  d0=. diag"2 y                                                         NB. 2×n-matrix
  d1=. ]`(9&o.)"1 d0                                                    NB. 2×n-matrix
  d2=. sorim`|"1 d1
  temp=. norm1tr"2 y                                                    NB. 2×n-matrix
  abnorm=. >./"1 temp                                                   NB. 2-vector
  abrwork=. temp - sorim d0
  abscale=. % FP_SFMIN >. abnorm                                        NB. 2-vector
  temp=. % (>./) FP_SFMIN , abscale * d2                                NB. n-vector
  sba=. |. abscale * temp *"1 d1                                        NB. 2×n-matrix
  abcoeff=. abscale * sba                                               NB. 2×n-matrix
  NB. scale to avoid underflow
  lsab=. *./ >:&FP_SFMIN`(<&small)"2 |`sorim"1"2 sba ,: abcoeff         NB. 2×n-matrix
  scale=. >./ lsab} 1 ,: ((% small) <. abnorm) * small % |`sorim"1 sba  NB. n-vector
  scale=. (+./ lsab)} scale ,: scale <. % FP_SFMIN * >./ 1 , |`sorim"1 abcoeff
  abcoeff=. lsab} (abcoeff *"1 scale) ,: (abscale * scale *"1 sba)
  abcoeffa=. |`sorim"1 abcoeff
  cond1=. +/ abcoeffa * abrwork
  dmin=. (# x) # ,: >./ FP_SFMIN , FP_PREC * abnorm * abcoeffa          NB. k×n-matrix
  d=. (x { |: abcoeff) mp ]`-"1 d0
  d=. (dmin < sorim d)} dmin ,: d
  cond2=. bignum (((1 > ]) ,. *) sorim) d

  bignum ; (|: d2) ; (|: abrwork) ; cond1 ; cond2 ; (|: abcoeff) ; (|: abcoeffa) ; d
)

NB. ---------------------------------------------------------
NB. tgevcly
NB.
NB. Description:
NB.   Compute some or all of non-scaled left eigenvectors
NB.
NB. Syntax:
NB.   W=. (iso ; init) tgevcly SP
NB. where
NB.   iso   - k-vector, lISO eigenvectors to compute
NB.   init  - boxed 8-vector, the output of tgevci
NB.   SP    - 2×n×n-matrix (S,:P), generalized Schur form,
NB.           produced by hgezqsxx
NB.   W     - k×n-matrix, some or all of left eigenvectors,
NB.           non-scaled
NB.   k     - integer in range [0,n]

tgevcly=: 4 : 0
  'iso bignum d2 abrwork cond1 cond2 abcoeff abcoeffa d'=. x
  n=. c y
  k=. # iso
  W=. (0,n) $ 0
  je=. <: k
  while. je >: 0 do.
    if. *./ FP_SFMIN >: (je { iso) { d2 do.
      NB. singular matrix pencil - return unit eigenvector
      work=. 1 je} n $ 0
    else.
      NB. non-singular eigenvalue: triangular solve of:
      NB.   y * (a*A - b*B) = 0  (rowwise)
      NB. work[0:j-1] contains sums w
      NB. work[j+1:je] contains y
      work=. 1 ,~ -/ ((je { iso) { abcoeff) * (((0 , ] , 0:) ,: 2 1 , ]) je { iso) {.@(1 0 2&|:);.0 y
      di=. je { d
      j=. <: je { iso
      while. j >: 0 do.
        NB. form:
        NB.   y[j] = - w[j] / di
        NB. with scaling and perturbation of the denominator
        abs1wj=. sorim j { work
        if. *.`<:/ (j { cond2) , abs1wj do.
          work=. work % abs1wj
        end.
        work=. j -@(%&(j{di)) upd work
        abs1wj=. sorim j { work
        if. j > 0 do.
          NB. w = w + y[j] * (a*S[:,j] - b*P[:,j]) with scaling
          if. ((abcoeffa mp&(j&{) abrwork) >: (bignum % abs1wj)) *. (1 < abs1wj) do.
            work=. work % abs1wj
          end.
          workadd=. (((je { iso) { abcoeff) * j { work) * (((0 , ] , 0:) ,: 2 1 , ]) j) {.@(1 0 2&|:);.0 y
          work=. (i. j) +`-/@(,&workadd) upd work
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
NB. tgevclx
NB.
NB. Description:
NB.   Compute some or all of non-scaled right eigenvectors
NB.
NB. Syntax:
NB.   W=. (iso ; init) tgevclx SP
NB. where
NB.   iso   - k-vector, lISO eigenvectors to compute
NB.   init  - boxed 8-vector, the output of tgevci
NB.   SP    - 2×n×n-matrix (S,:P), generalized Schur form,
NB.           produced by hgezqsxx
NB.   W     - k×n-matrix, some or all of right eigenvectors,
NB.           non-scaled
NB.   k     - integer in range [0,n]

tgevclx=: 4 : 0
  'iso bignum d2 abrwork cond1 cond2 abcoeff abcoeffa d'=. x
  d=. + d
  n=. c y
  k=. # iso
  W=. (0,n) $ 0
  je=. 0
  while. je < k do.
    if. *./ FP_SFMIN >: (je { iso) { d2 do.
      NB. singular matrix pencil - return unit eigenvector
      work=. 1 je} n $ 0
    else.
      NB. non-singular eigenvalue: triangular solve of:
      NB.   x * (a*A - b*B)^H = 0 ,
      NB. columnwise in (a*A - b*B)^H , or rowwise in
      NB. (a*A - b*B)
      work=. 1
      xmax=. 1
      di=. je { d
      j=. >: je { iso
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
        sum=. -/ (]`+"0 (je { iso) { abcoeff) * work mp (j ((0 , ,) ,: (2 1 , -)) je { iso) +@{.@(0&|:);.0 y
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
      work=. (- n) {. work
    end.
    je=. >: je
    W=. W , work
  end.
  W
)

NB. ---------------------------------------------------------
NB. tgevcs
NB.
NB. Description:
NB.   Scale left or right eigenvectors
NB.
NB. Syntax:
NB.   V=. tgevcs W
NB. where
NB.   W - k×n-matrix, some or all of left or right
NB.       eigenvectors, non-scaled
NB.   V - k×n-matrix, scaled W

tgevcs=: 3 : 0
  norm=. normitr y
  iso=. (#y) #"0 FP_SFMIN < norm
  y=. y % norm
  y=. iso} 0 ,: y
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. tgevcll
NB. tgevclr
NB. tgevclb
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
NB.   Y=.     [iso] tgevcll SP
NB.   X=.     [iso] tgevclr SP
NB.   'Y X'=. [iso] tgevclb SP
NB. where
NB.   iso - k-vector, optional lISO eigenvectors to compute,
NB.         default is "all eigenvectors"
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgezqsxx
NB.   Y   - k×n-matrix, some or all of left eigenvectors
NB.   X   - k×n-matrix, some or all of right eigenvectors
NB.   k   - integer in range [0,n]
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (tgevcll -:      tgevcur&.:(|:"2)) SP
NB.   (tgevclr -:      tgevcul&.:(|:"2)) SP
NB.   (tgevclb -: 1 A. tgevcub&.:(|:"2)) SP
NB.   (E2 mp Y mp S) -: (E1 mp Y mp P)
NB.   (S mp (ct X) mp E2) -: (P mp (ct X) mp E1)
NB. where
NB.   'Y X'=. tgevclb SP
NB.   'E1 E2'=. diagmat@diag"2 SP

tgevcll=:  ($:~ i.@c) : (([ ; tgevci) tgevcs  @ tgevcly             ])
tgevclr=:  ($:~ i.@c) : (([ ; tgevci) tgevcs  @            tgevclx  ])
tgevclb=:  ($:~ i.@c) : (([ ; tgevci) tgevcs"2@(tgevcly ,: tgevclx) ])

NB. ---------------------------------------------------------
NB. tgevcul
NB. tgevcur
NB. tgevcub
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
NB.   Y=.     [iso] tgevcul SP
NB.   X=.     [iso] tgevcur SP
NB.   'Y X'=. [iso] tgevcub SP
NB. where
NB.   iso - k-vector, optional lISO eigenvectors to compute,
NB.         default is "all eigenvectors"
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgeqzsxx
NB.   Y   - n×k-matrix, some or all of left eigenvectors
NB.   X   - n×k-matrix, some or all of right eigenvectors
NB.   k   - integer in range [0,n]
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (tgevcul -:      tgevclr&.:(|:"2)) SP
NB.   (tgevcur -:      tgevcll&.:(|:"2)) SP
NB.   (tgevcub -: 1 A. tgevclb&.:(|:"2)) SP
NB.   (E2 mp (ct Y) mp S) -: (E1 mp (ct Y) mp P)
NB.   (S mp X mp E2) -: (P mp X mp E1)
NB. where
NB.   'Y X'=. tgevcub SP
NB.   'E1 E2'=. diagmat@diag"2 SP
NB.
NB. Notes:
NB. - tgevcul models LAPACK's xTGEVC('L','S')
NB. - tgevcur models LAPACK's xTGEVC('R','S')
NB. - tgevcub models LAPACK's xTGEVC('B','S')

tgevcul=: ($:~ i.@c) : ((([ ; tgevci) tgevcs  @ tgevclx             ])&.:(|:"2))
tgevcur=: ($:~ i.@c) : ((([ ; tgevci) tgevcs  @            tgevcly  ])&.:(|:"2))
tgevcub=: ($:~ i.@c) : ((([ ; tgevci) tgevcs"2@(tgevclx ,: tgevcly) ])&.:(|:"2))

NB. ---------------------------------------------------------
NB. tgevcllb
NB. tgevclrb
NB. tgevclbb
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
NB.   YQ=.      tgevcllb SP , Q
NB.   XZ=.      tgevclrb SP , Z
NB.   'YQ XZ'=. tgevclbb SP , Q ,: Z
NB. where
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgezqsxx
NB.   YQ  - n×n-matrix, left eigenvectors Y*Q
NB.   XZ  - n×n-matrix, right eigenvectors X*Z
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (tgevcllb -: tgevcurb&.:(|:"2        )) SP , Q2
NB.   (tgevclrb -: tgevculb&.:(|:"2        )) SP , Z2
NB.   (tgevclbb -: tgevcubb&.:(|:"2@:(1&A.))) SP , Q2 ,: Z2
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
NB.   'H T Q1 Z1'=. hs gghrdlvv A , B , Q0 ,: Z0
NB.   'H T dQ0 dZ0'=. hs gghrdlvv A , B , ,:~ I
NB.   'S P Q2 Z2'=. hs hgezqsvv H , T , Q1 ,: Z1
NB.   'S P dQ1 dZ1'=. hs hgezqsvv H , T , ,:~ I
NB.   'S P dQ1dQ0 dZ1dZ0'=. hs hgezqsvv H , T , dQ0 ,: dZ0
NB.   'Y X'=. tgevclb S ,: P
NB.   'YdQ1 XdZ1'=. tgevclbb S , P , dQ1 ,: dZ1
NB.   'YdQ1dQ0 XdZ1dZ0'=. tgevclbb S , P , dQ1dQ0 ,: dZ1dZ0
NB.   'YQ2 XZ2'=. tgevclbb S , P , Q2 ,: Z2
NB.   'E1 E2'=. diagmat@diag"2 S ,: P

tgevcllb=:  (i.@c ([ ; tgevci) 2&{.) tgevcs  @( tgevcly mp 2 { ]                        ) ]
tgevclrb=:  (i.@c ([ ; tgevci) 2&{.) tgevcs  @(                       tgevclx mp _1 { ] ) ]
tgevclbb=:  (i.@c ([ ; tgevci) 2&{.) tgevcs"2@((tgevcly mp 2 { ]) ,: (tgevclx mp _1 { ])) ]

NB. ---------------------------------------------------------
NB. tgevculb
NB. tgevcurb
NB. tgevcubb
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
NB.   QY=.      tgevculb SP , Q
NB.   ZX=.      tgevcurb SP , Z
NB.   'QY ZX'=. tgevcubb SP , Q ,: Z
NB. where
NB.   SP  - 2×n×n-matrix (S,:P), generalized Schur form,
NB.         produced by hgeqzsxx
NB.   QY  - n×n-matrix, left eigenvectors Q*Y
NB.   ZX  - n×n-matrix, right eigenvectors Z*X
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (tgevculb -: tgevclrb&.:(|:"2        )) SP , Q2
NB.   (tgevcurb -: tgevcllb&.:(|:"2        )) SP , Z2
NB.   (tgevcubb -: tgevclbb&.:(|:"2@:(1&A.))) SP , Q2 ,: Z2
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
NB.   'H T Q1 Z1'=. hs gghrduvv A , B , Q0 ,: Z0
NB.   'H T dQ0 dZ0'=. hs gghrduvv A , B , ,:~ I
NB.   'S P Q2 Z2'=. hs hgeqzsvv H , T , Q1 ,: Z1
NB.   'S P dQ1 dZ1'=. hs hgeqzsvv H , T , ,:~ I
NB.   'S P dQ0dQ1 dZ0dZ1'=. hs hgeqzsvv H , T , dQ0 ,: dZ0
NB.   'Y X'=. tgevcub S ,: P
NB.   'dQ1Y dZ1X'=. tgevcubb S , P , dQ1 ,: dZ1
NB.   'dQ0dQ1Y dZ0dZ1X'=. tgevcubb S , P , dQ0dQ1 ,: dZ0dZ1
NB.   'Q2Y Z2X'=. tgevcubb S , P , Q2 ,: Z2
NB.   'E1 E2'=. diagmat@diag"2 S ,: P
NB.
NB. Notes:
NB. - tgevculb models LAPACK's xTGEVC('L','B')
NB. - tgevcurb models LAPACK's xTGEVC('R','B')
NB. - tgevcubb models LAPACK's xTGEVC('B','B')

tgevculb=: ((i.@c ([ ; tgevci) 2&{.) tgevcs  @( tgevclx mp 2{]                        ) ])&.:(|:"2)
tgevcurb=: ((i.@c ([ ; tgevci) 2&{.) tgevcs  @(                     tgevcly mp _1 { ] ) ])&.:(|:"2)
tgevcubb=: ((i.@c ([ ; tgevci) 2&{.) tgevcs"2@((tgevclx mp 2{]) ,: (tgevcly mp _1 { ])) ])&.:(|:"2)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testtgevc
NB.
NB. Description:
NB.   Test tgevcxxx by pair of square matrices
NB.
NB. Syntax:
NB.   testtgevc AB
NB. where
NB.   AB - 2×n×n-brick

testtgevc=: 3 : 0
  rcondl=. <./ trlconi"2 SPl=. 2 {. SPQZHTl=. (([ (((0,[) hgezqsvv (, ,:~@idmat)~) , ]) ((gghrdlnn~ 0&,)~ ((unmlqrc~ ,: trl@:(}:"1)@]) gelqf)/))~ c) y
  rcondu=. <./ trucon1"2 SPu=. 2 {. SPQZHTu=. (([ (((0,[) hgeqzsvv (, ,:~@idmat)~) , ]) ((gghrdunn~ 0&,)~ ((unmqrlc~ ,: tru@  }:   @]) geqrf)/))~ c) y

  ('tgevcll'  tmonad (]          `]`(rcondl"_)`(_."_)`t52ll)) SPl
  ('tgevclr'  tmonad (]          `]`(rcondl"_)`(_."_)`t52lr)) SPl
  ('tgevclb'  tmonad (]          `]`(rcondl"_)`(_."_)`t52lb)) SPl
  ('tgevcllb' tmonad ((0 1 2  &{)`]`(rcondl"_)`(_."_)`t52ll)) SPQZHTl
  ('tgevclrb' tmonad ((0 1   3&{)`]`(rcondl"_)`(_."_)`t52lr)) SPQZHTl
  ('tgevclbb' tmonad ((0 1 2 3&{)`]`(rcondl"_)`(_."_)`t52lb)) SPQZHTl
  ('tgevcul'  tmonad (]          `]`(rcondu"_)`(_."_)`t52ul)) SPu
  ('tgevcur'  tmonad (]          `]`(rcondu"_)`(_."_)`t52ur)) SPu
  ('tgevcub'  tmonad (]          `]`(rcondu"_)`(_."_)`t52ub)) SPu
  ('tgevculb' tmonad ((0 1 2  &{)`]`(rcondu"_)`(_."_)`t52ul)) SPQZHTu
  ('tgevcurb' tmonad ((0 1   3&{)`]`(rcondu"_)`(_."_)`t52ur)) SPQZHTu
  ('tgevcubb' tmonad ((0 1 2 3&{)`]`(rcondu"_)`(_."_)`t52ub)) SPQZHTu

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
NB.     _1 1 0 4 _6 4&gemat_mt_ testevc_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testevc_mt_ 150 150

testevc=: 1 : 'EMPTY [ testtgevc_mt_@u@(2&,)^:(=/)'
