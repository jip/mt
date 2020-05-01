NB. Eigenvalues and eigenvectors
NB.
NB. ggevxxx   Eigenvalues and, optionally, eigenvectors of
NB.           pair of matrices
NB.
NB. testgeev  Test geevxxx by square matrix
NB. testheev  Test heevxx by Hermitian (symmetric) matrix
NB. testggev  Test ggevxxx by pair of square matrices
NB. testev    Adv. to make verb to test xxevxxx by matrices
NB.           of generator and shape given
NB.
NB. Version: 0.10.5 2020-03-30
NB.
NB. Copyright 2010-2020 Igor Zhuravlov
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
NB. Concepts
NB.
NB. I have a dream (2010-11-02):
NB.   ggev=: gizmo @ tgevc @ hgeqz @ gghrd @ qr &. bal &. scl
NB. Is it possible?

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Constants

NB. Scaling limits
EVBIGNUM=: % EVSMLNUM=: (%: FP_SFMIN) % FP_PREC

NB. Scaling factors
EVSCL=: 1 , EVSMLNUM , 1 , EVBIGNUM

NB. ---------------------------------------------------------
NB. ggevi
NB.
NB. Description:
NB.   Adv. to make verb to calculate initial parameters for
NB.   ggevxxx
NB.
NB. Syntax:
NB.   vapp=. ggbalp ggevi
NB. where
NB.   ggbalp  - monad to permute matrix pair (A,B) to isolate
NB.             eigenvalues, is either ggballp or ggbalup, is
NB.             called as:
NB.               'CD plr hs'=. ggbalp AB
NB.   vapp    - monad to calculate initial parameters for
NB.             ggevxxx, is called as:
NB.               'abnrmio ABupd plr hs'=. vapp AB
NB.   AB      - 2×n×n-matrix, matrix pair (A,B)
NB.   abnrmio -:abnrm ,. abio
NB.   abnrm   - 2-vector, norms of A and B
NB.   abio    - 2-vector of integers, defines both necessity
NB.             and value of scaling for A and B
NB.   ABupd   - 2×n×n-matrix, scaled and permuted A and B
NB.   plr     - 2×n-matrix of integers, permutations of A and
NB.             B, produced by ggbalp
NB.   hs      - 2-vector of integers, defines submatrices
NB.             position, produced by ggbalp

ggevi=: 1 : '(,.(0,(EVSMLNUM_mt_*1-FP_EPS_mt_),EVBIGNUM_mt_)&I.)@:(normm_mt_"2) ([ ; u@:(scl_mt_^:((,{&EVSCL_mt_)/@[`({&0 1 0 1@{:@[)`])"1 2)) ]'

NB. ---------------------------------------------------------
NB. drgev
NB.
NB. Description:
NB.   Adv. to make verb to compute backward error of
NB.   eigenvectors produced by generalized eigenvalue problem
NB.   solvers
NB.
NB. Syntax:
NB.   vberr=. mmul`vmul`norma`normb drgev
NB. where
NB.   mmul     - dyad to multiply matrices; is called as:
NB.                M3=. M1 mmul M2
NB.   vmul     - dyad to multiply matrix by diagonal matrix
NB.              represented as vector; is called as:
NB.                M2=. v1 vmul M1
NB.   norma    - monad to compute norm of matrix; is called
NB.              as:
NB.                normM=. norma M1
NB.   normb    - monad to compute norms of matrix rows
NB.              (columns); is called as:
NB.                v2=. normb M1
NB.   vberr    - dyad to compute berr; is called as:
NB.                berr=. AB vberr (e1e2 ; V)
NB.   AB       - 2×n×n-brick, matrix pair (A,B) to
NB.              eigen-decompose
NB.   e1e2     - 2×n-matrix, laminated vectors of
NB.              generalized eigenvalues α and β
NB.   V        - either L or R
NB.   M1,M2,M3 - n×n-matrix
NB.   v1,v2    - n-vector
NB.   L        - n×n-matrix, columns with left eigenvectors
NB.   R        - n×n-matrix, columns with right eigenvectors
NB.   berr     ≥ 0, float, scalar, backward error
NB.   normM    ≥ 0, float, scalar, norm of matrix
NB.
NB. Application:
NB. - to compute berr for L from ggevlvx:
NB.     NB. berr=. AB vberrlL (e1e2 ; L)
NB.     vberrlL=:  mp_mt_         "2` *   `normi_mt_`normitr_mt_ drgev
NB. - to compute berr for R from ggevlxv:
NB.     NB. berr=. AB vberrlR (e1e2 ; R)
NB.     vberrlR=: (mp_mt_  ct_mt_)"2`(*"1)`normi_mt_`normitr_mt_ drgev
NB. - to compute berr for L from ggevuvx:
NB.     NB. berr=. AB vberruL (e1e2 ; L)
NB.     vberruL=: (mp_mt_~ ct_mt_)"2` *   `norm1_mt_`normitc_mt_ drgev
NB. - to compute berr for R from ggevuxv:
NB.     NB. berr=. AB vberruR (e1e2 ; R)
NB.     vberruR=: mp_mt_          "2`(*"1)`norm1_mt_`normitc_mt_ drgev
NB.
NB. Formula:
NB.   berr := max(berr0,berr1)
NB.   - for ggevlxx:
NB.     ||L|| := max(||L||_inf , FP_PREC)
NB.     ||R|| := max(||R||_inf , FP_PREC)
NB.     - for L:
NB.         berr0 := (||(C2 * (C1 * E2)) * L   * A - (C2 * (C1 * E2)) * L   * B||_1 / ||R||) / FP_PREC
NB.         berr1 := normit(normitc(L) - 1) / (FP_PREC * n)
NB.     - for R:
NB.         berr0 := (||A * R^H * ((E2 * С1) * С2) - B * R^H * ((E1 * С1) * С2)||_1 / ||R||) / FP_PREC
NB.         berr1 := normit(normitc(R) - 1) / (FP_PREC * n)
NB.     - for L and R:
NB.         berr0 := max(berr0(ggevlvx),berr0(ggevlxv))
NB.         berr1 := max(berr1(ggevlvx),berr1(ggevlxv))
NB.   - for ggevuxx:
NB.     ||L|| := max(||L||_1 , FP_PREC)
NB.     ||R|| := max(||R||_1 , FP_PREC)
NB.     - for L:
NB.         berr0 := (||(C2 * (C1 * E2)) * L^H * A - (C2 * (C1 * E2)) * L^H * B||_1 / ||R||) / FP_PREC
NB.         berr1 := normit(normitc(L) - 1) / (FP_PREC * n)
NB.     - for R:
NB.         berr0 := (||A * R   * ((E2 * С1) * С2) - B * R   * ((E1 * С1) * С2)||_1 / ||R||) / FP_PREC
NB.         berr1 := normit(normitc(R) - 1) / (FP_PREC * n)
NB.     - for L and R:
NB.         berr0 := max(berr0(ggevuvx),berr0(ggevuxv))
NB.         berr1 := max(berr1(ggevuvx),berr1(ggevuxv))
NB.   C1 := (diag(coeff1))^_1
NB.   C2 := (diag(coeff2))^_1
NB.   coeff1(i) := if(sorim(α(i)) > (1/FP_SFMIN) / ||B||
NB.                or sorim(β(i)) > (1/FP_SFMIN) / ||A||
NB.                or max(sorim(α(i)),sorim(β(i))) < 1)
NB.                then max(max(sorim(α(i)),sorim(β(i))),FP_SFMIN)
NB.                else 1
NB.   coeff2(i) := max(sorim(α(i))*||B||,sorim(β(i))*||A||, FP_SFMIN)
NB.   ||A|| := max(||A||_1 , FP_SFMIN)
NB.   ||B|| := max(||B||_1 , FP_SFMIN)
NB. where
NB.   v(i) - i-th element of vector v
NB.
NB. Notes:
NB. - models LAPACK's xDRGEV and xGET52

drgev=: 1 : 0
  '`mmul vmul norma normb'=. m
  'e1e2 V'=. y
  n=. # V
  if. 0 = n do. 0 return. end.
  banorm=. |. FP_SFMIN_mt_ >. norma"2 x                       NB. 2-vector, float
  alfbetmax=. (% FP_SFMIN_mt_) % 1 >. banorm                  NB. 2-vector, float
  abs1ab=. sorim_mt_"1 e1e2                                   NB. 2×n-matrix, float
  abmax=. >./ abs1ab                                          NB. n-vector, float
  cond=. +./ (abs1ab > alfbetmax) , 1 > abmax                 NB. n-vector, boolean
  e1e2=. e1e2 %"1 cond} 1 ,: FP_SFMIN_mt_ >. abmax            NB. 2×n-matrix
  abcoeff=. (|. e1e2) %"1 >./ FP_SFMIN_mt_ , abs1ab * banorm  NB. 2×n-matrix
  Err=. -/ abcoeff vmul x mmul V                              NB. n×n-matrix
  errnrm=. (norma Err) % FP_PREC_mt_ >. norma V               NB. scalar, float
  result1=. errnrm % FP_PREC_mt_                              NB. scalar, float
  enrmer=. normir_mt_ <: normb V                              NB. scalar, float
  result2=. enrmer % FP_PREC_mt_ * n                          NB. scalar, float
  result1 >. result2                                          NB. scalar, float
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. ggevlnn
NB. ggevlnv
NB. ggevlvn
NB. ggevlvv
NB.
NB. Description:
NB.   Generalized nonsymmetric eigenvalue problem (GNEP):
NB.   find eigenvalue vectors e1, e2 and, optionally, left
NB.   eigenvectors L:
NB.     E2 * L * A = E1 * L * B                           (1)
NB.   and/or right eigenvectors R:
NB.     A * R^H * E2 = B * R^H * E1                       (2)
NB.   of pair of matrices (A,B). To avoid overflow,
NB.   eigenvalues of the matrix pair (A,B) are computed as a
NB.   pair of values. Each i-th eigenvector (row) from L and
NB.   R has a corresponding eigenvalue represented as a pair
NB.   of i-th elements from e1 and e2:
NB.     E1=. diagmat e1
NB.     E2=. diagmat e2
NB.   If E2 is nonsingular then:
NB.     E=. diagmat e1%e2
NB.   is a diagonal matrix of eigenvalues, and GNEP (1), (2)
NB.   can be expressed as:
NB.     L * A = E * L * B                                 (3)
NB.     A * R^H = B * R^H * E                             (4)
NB.   and if E1 is nonsingular then:
NB.     E=. diagmat e2%e1
NB.   is a diagonal matrix of eigenvalues, and GNEP (1), (2)
NB.   can be expressed as:
NB.     E * L * A = L * B                                 (5)
NB.     A * R^H * E = B * R^H * E                         (6)
NB.   Eigenvectors are normalized to have taxicab-based
NB.   ∞-norm equal to 1
NB.
NB. Syntax:
NB.   e1e2=.      ggevlnn AB
NB.   'e1e2 R'=.  ggevlnv AB
NB.   'e1e2 L'=.  ggevlvn AB
NB.   'e1e2 LR'=. ggevlvv AB
NB. where
NB.   AB   - 2×n×n-matrix, matrix pair (A,B):
NB.            AB -: A ,: B
NB.   e1e2 - 2×n-matrix of eigenvalues e1 and e2:
NB.            e1e2 -: e1 ,: e2
NB.   L    - n×n-matrix, left eigenvectors (rows)
NB.   R    - n×n-matrix, right eigenvectors (rows)
NB.   LR   - 2×n×n-matrix, left and right eigenvectors:
NB.            LR -: L ,: R
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (ggevlnn -:                    +@ggevunn@:(ct"2)) A ,: B
NB.   (ggevlnn -:&(/:~@(%/))           ggevunn        ) A ,: B
NB.   (ggevlnv -: +&.>`(ct      &.>)"0@ggevuvn@:(ct"2)) A ,: B
NB.   (ggevlvn -: +&.>`(ct      &.>)"0@ggevunv@:(ct"2)) A ,: B
NB.   (ggevlvv -: +&.>`(ct"2@:|.&.>)"0@ggevuvv@:(ct"2)) A ,: B
NB.   (E2 mp L mp A) -: (E1 mp L mp B)
NB.   (A mp (ct R) mp E2) -: (B mp (ct R) mp E1)
NB. where
NB.   A - n×n-matrix, general
NB.   B - n×n-matrix, general
NB.   'e1e2 LR'=. ggevlvv A ,: B
NB.   'E1 E2'=. diagmat"1 e1e2
NB.   'L R'=. LR
NB.
NB. Application:
NB. - simulate LAPACK's xGEEV('N','N'):
NB.     NB. e=. geevlnn A
NB.     geevlnn=: {.@ggevlnn@(,: idmat@c)
NB. - simulate LAPACK's xGEEV('N','V') (see notes):
NB.     NB. 'e R'=. geevlnv A
NB.     geevlnv=: {.&.>`(((* *@+@((i. >./)"1@sorim{"0 1]))%normsr)      &.>)"0@ggevlnv@(,: idmat@c)
NB. - simulate LAPACK's xGEEV('V','N') (see notes):
NB.     NB. 'e L'=. geevlvn A
NB.     geevlvn=: {.&.>`(((* *@+@((i. >./)"1@sorim{"0 1]))%normsr)      &.>)"0@ggevlvn@(,: idmat@c)
NB. - simulate LAPACK's xGEEV('V','V') (see notes):
NB.     NB. 'e LR'=. geevlvv A
NB.     geevlvv=: {.&.>`(((* *@+@((i. >./)"1@sorim{"0 1]))%normsr)    "2&.>)"0@ggevlvv@(,: idmat@c)
NB. - simulate LAPACK's xHEEV('N'):
NB.     NB. e=. heevln A
NB.     heevln=: 9 o. {.@ggevlnn@(,: idmat@c)
NB. - simulate LAPACK's xHEEV('V') (see notes):
NB.     NB. 'e V'=. heevlv A
NB.     heevlv=: (9 o. {.)&.>`((%  %:@diag@(mp ct))&.>)"0@ggevlvn@(,: idmat@c)
NB.
NB. Notes:
NB. - eigenvectors from LAPACK's xGEEV are normalized to have
NB.   Euclidean norm equal to 1 and largest component real
NB. - eigenvectors from LAPACK's xHEEV are orthonormal

ggevlnn=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi y
  y=. (<0 1;;~dhs2liso hs) ([ ((gghrdlnn~0,c) upd) ((unmlqrc~,:trl@:(}:"1)@])gelqf)/@{`[`]}) y
  e1e2=. hs hgezqenn y
  e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
)

ggevlnv=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi y
  y=. (0 1;(<i.{.hs);dhs2liso hs) ((unmlqrc~,(trl@:(}:"1),:unglq)@])gelqf)/@({~<)~`((<0 1 2)<@(0})[)`((, idmat@c)@])} y
  y=. (gghrdlnv~0,c) y
  y=. hs hgezqsnv y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    y=. _. $~ 2 #c y
  else.
    y=. tgevclrb y
    y=. gebaklp y ; {: plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr) y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
  end.
  e1e2 ; y
)

ggevlvn=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi y
  y=. (<0 1;(<i.{.hs);dhs2liso hs) ((unmlqrc~,:trl@:(}:"1)@])gelqf)/@{`[`]} y
  y=. (((0,]) gghrdlvn (,idmat)) c) y
  y=. hs hgezqsvn y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    y=. _. $~ 2 #c y
  else.
    y=. tgevcllb y
    y=. gebaklp y ; {. plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr) y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
  end.
  e1e2 ; y
)

ggevlvv=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi y
  y=. (0 1;(<i.{.hs);dhs2liso hs) ((unmlqrc~,(trl@:(}:"1),:unglq)@])gelqf)/@({~<)~`((<0 1 3)<@(0})[)`((, ,:~@idmat@c)@])} y
  y=. (gghrdlvv~0,c) y
  y=. hs hgezqsvv y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    y=. _. $~ $ y
  else.
    y=. tgevclbb y
    y=. y gebaklp@;"2 1 plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr)"2 y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
  end.
  e1e2 ; y
)

NB. ---------------------------------------------------------
NB. ggevunn
NB. ggevunv
NB. ggevuvn
NB. ggevuvv
NB.
NB. Description:
NB.   Generalized nonsymmetric eigenvalue problem (GNEP):
NB.   find eigenvalue vectors e1, e2 and, optionally, left
NB.   eigenvectors L:
NB.     E2 * L^H * A = E1 * L^H * B                       (7)
NB.   and/or right eigenvectors R:
NB.     A * R * E2 = B * R * E1                           (8)
NB.   of pair of matrices (A,B). To avoid overflow,
NB.   eigenvalues of the matrix pair (A,B) are computed as a
NB.   pair of values. Each i-th eigenvector (column) from L
NB.   and R has a corresponding eigenvalue represented as a
NB.   pair of i-th elements from e1 and e2:
NB.     E1=. diagmat e1
NB.     E2=. diagmat e2
NB.   If E2 is nonsingular then:
NB.     E=. diagmat e1%e2
NB.   is a diagonal matrix of eigenvalues, and GNEP (7), (8)
NB.   can be expressed as:
NB.     L^H * A = E * L^H * B                             (9)
NB.     A * R = B * R * E                                (10)
NB.   and if E1 is nonsingular then:
NB.     E=. diagmat e2%e1
NB.   is a diagonal matrix of eigenvalues, and GNEP (7), (8)
NB.   can be expressed as:
NB.     E * L^H * A = L^H * B                            (11)
NB.     A * R * E = B * R * E                            (12)
NB.   Eigenvectors are normalized to have taxicab-based
NB.   ∞-norm equal to 1
NB.
NB. Syntax:
NB.   e1e2=.      ggevunn AB
NB.   'e1e2 R'=.  ggevunv AB
NB.   'e1e2 L'=.  ggevuvn AB
NB.   'e1e2 LR'=. ggevuvv AB
NB. where
NB.   AB   - 2×n×n-matrix, matrix pair (A,B):
NB.            AB -: A ,: B
NB.   e1e2 - 2×n-matrix of eigenvalues e1 and e2:
NB.            e1e2 -: e1 ,: e2
NB.   L    - n×n-matrix, left eigenvectors (columns)
NB.   R    - n×n-matrix, right eigenvectors (columns)
NB.   LR   - 2×n×n-matrix, left and right eigenvectors:
NB.            LR -: L ,: R
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (ggevunn -:                    +@ggevlnn@:(ct"2)) A ,: B
NB.   (ggevunn -:&(/:~@(%/))           ggevlnn        ) A ,: B
NB.   (ggevunv -: +&.>`(ct      &.>)"0@ggevlvn@:(ct"2)) A ,: B
NB.   (ggevuvn -: +&.>`(ct      &.>)"0@ggevlnv@:(ct"2)) A ,: B
NB.   (ggevuvv -: +&.>`(ct"2@:|.&.>)"0@ggevlvv@:(ct"2)) A ,: B
NB.   (E2 mp (ct L) mp A) -: (E1 mp (ct L) mp B)
NB.   (A mp R mp E2) -: (B mp R mp E1)
NB. where
NB.   A - n×n-matrix, general
NB.   B - n×n-matrix, general
NB.   'e1e2 LR'=. ggevuvv A ,: B
NB.   'E1 E2'=. diagmat"1 e1e2
NB.   'L R'=. LR
NB.
NB. Application:
NB. - simulate LAPACK's xGEEV('N','N'):
NB.     NB. e=. geevunn A
NB.     geevunn=: {.@ggevunn@(,: idmat@c)
NB. - simulate LAPACK's xGEEV('N','V') (see notes):
NB.     NB. 'e R'=. geevunv A
NB.     geevunv=: {.&.>`(((* *@+@((i. >./)"1@sorim{"0 1]))%normsr)&.|:  &.>)"0@ggevunv@(,: idmat@c)
NB. - simulate LAPACK's xGEEV('V','N') (see notes):
NB.     NB. 'e L'=. geevuvn A
NB.     geevuvn=: {.&.>`(((* *@+@((i. >./)"1@sorim{"0 1]))%normsr)&.|:  &.>)"0@ggevuvn@(,: idmat@c)
NB. - simulate LAPACK's xGEEV('V','V') (see notes):
NB.     NB. 'e LR'=. geevuvv A
NB.     geevuvv=: {.&.>`(((* *@+@((i. >./)"1@sorim{"0 1]))%normsr)&.|:"2&.>)"0@ggevuvv@(,: idmat@c)
NB. - simulate LAPACK's xHEEV('N'):
NB.     NB. e=. heevun A
NB.     heevun=: 9 o. {.@ggevunn@(,: idmat@c)
NB. - simulate LAPACK's xHEEV('V') (see notes):
NB.     NB. 'e V'=. heevuv A
NB.     heevuv=: (9 o. {.)&.>`((%"1%:@diag@(mp~ct))&.>)"0@ggevunv@(,: idmat@c)
NB.
NB. Notes:
NB. - ggevunn models LAPACK's xGGEV('N','N')
NB. - ggevunv models LAPACK's xGGEV('N','V')
NB. - ggevuvn models LAPACK's xGGEV('V','N')
NB. - ggevuvv models LAPACK's xGGEV('V','V')
NB. - eigenvectors from LAPACK's xGEEV are normalized to have
NB.   Euclidean norm equal to 1 and largest component real
NB. - eigenvectors from LAPACK's xHEEV are orthonormal

ggevunn=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi y
  y=. (<0 1;;~dhs2liso hs) ([ ((gghrdunn~0,c) upd) ((unmqrlc~,:tru@}:@])geqrf)/@{`[`]}) y
  e1e2=. hs hgeqzenn y
  e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
)

ggevuvn=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi y
  y=. (0 1;(dhs2liso hs);<<i.{.hs) ((unmqrlc~,(tru@}:,:ungqr)@])geqrf)/@({~<)~`((<0 1 2)<@(0})[)`((, idmat@c)@])} y
  y=. (gghrduvn~0,c) y
  y=. hs hgeqzsvn y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    y=. _. $~ 2 #c y
  else.
    y=. tgevculb y
    y=. gebakup y ; {. plr
    y=. (%"1 (EVSMLNUM&>`(,:&1))}@:normitc) y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
  end.
  e1e2 ; y
)

ggevunv=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi y
  y=. (<0 1;(dhs2liso hs);<<i.{.hs) ((unmqrlc~,:tru@}:@])geqrf)/@{`[`]} y
  y=. (((0,]) gghrdunv (,idmat)) c) y
  y=. hs hgeqzsnv y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    y=. _. $~ 2 #c y
  else.
    y=. tgevcurb y
    y=. gebakup y ; {: plr
    y=. (%"1 (EVSMLNUM&>`(,:&1))}@:normitc) y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
  end.
  e1e2 ; y
)

ggevuvv=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi y
  y=. (0 1;(dhs2liso hs);<<i.{.hs) ((unmqrlc~,(tru@}:,:ungqr)@])geqrf)/@({~<)~`((<0 1 2)<@(0})[)`((, ,:~@idmat@c)@])} y
  y=. (gghrduvv~0,c) y
  y=. hs hgeqzsvv y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    y=. _. $~ $ y
  else.
    y=. tgevcubb y
    y=. y gebakup@;"2 1 plr
    y=. (%"1 (EVSMLNUM&>`(,:&1))}@:normitc)"2 y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
  end.
  e1e2 ; y
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgeev
NB.
NB. Description:
NB.   Test xGEEV (math/lapack2) by square matrix
NB.
NB. Syntax:
NB.   testgeev A
NB. where
NB.   A - n×n-matrix
NB.
NB. Formula:
NB.   berru    := max(berruL,berruR,berrul,berrur)
NB.   berruL   := min((||L^H * A - W * L^H||_1 / max(||L||_1 , FP_PREC)) / ||A|| , 1) / FP_PREC
NB.   berruR   := min((||A   * R - R * W  ||_1 / max(||R||_1 , FP_PREC)) / ||A|| , 1) / FP_PREC
NB.   ||A||    := max(||A||_1 , FP_SFMIN)
NB.
NB.   vrmax[j] := max(|Re(X[i,j])|)
NB.                i | Im(X[i,j]) == 0
NB.
NB.   vmax[j]  := max(|X[:,j]|)
NB.                i
NB.
NB.   if ∃ j | vrmax[j] / vmax[j] < 1 - 2 * FP_PREC then
NB.     berrux := 1 / FP_PREC
NB.   else
NB.     berrux := max(min(| ||X[:,j]||_E - 1 | , 1) / FP_PREC)
NB.               j
NB.   endif
NB. where
NB.   A      - matrix to eigen-decompose
NB.   W      - diagonal matrix of eigenvalues
NB.   L      - matrix of left eigenvectors
NB.   R      - matrix of right eigenvectors
NB.   berrux - either berrul or berrur
NB.   X      - either L or R
NB.
NB. Notes:
NB. - vberrux models LAPACK's xDRVEV and xGET22

testgeev=: 3 : 0
  load_mttmp_ :: ] '~addons/math/mt/test/lapack2/geev.ijs'

  rcond=. gecon1 y

  vberruX_mttmp_=. 2 : 'FP_PREC_mt_ %~ 1 <. (u % FP_PREC_mt_ >. norm1_mt_@v) % FP_SFMIN_mt_ >. norm1_mt_@['
  vberruL_mttmp_=. normi_mt_@(((mp_mt_~ ct_mt_)~ 1&{::) - (1 {:: ]) *"1 +@(0 {:: ])) vberruX_mttmp_ (1 {:: ])
  vberruR_mttmp_=. norm1_mt_@(( mp_mt_           2&{::) - (2 {:: ]) *"1   (0 {:: ])) vberruX_mttmp_ (2 {:: ])
  vberrux_mttmp_=. (% FP_PREC)"_`(>./@(FP_PREC %~ 1 ([ <. |@:-) normsc_mt_))@.((1 _2 p. FP_PREC) *./@:<: (%~&(>./@:|) ({."1 (* 0 = *) {:"1)@:+.))
  vberrul_mttmp_=. vberrux_mttmp_@(1 {:: ])
  vberrur_mttmp_=. vberrux_mttmp_@(2 {:: ])

  ('''nn''&dgeev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(_."_)                                                                )) y
  ('''nv''&dgeev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(                                    vberruR_mttmp_ >. vberrur_mttmp_))) y
  ('''vn''&dgeev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(vberruL_mttmp_ >. vberrul_mttmp_                                    ))) y
  ('''vv''&dgeev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(vberruL_mttmp_ >. vberrul_mttmp_ >. vberruR_mttmp_ >. vberrur_mttmp_))) y

  ('''nn''&zgeev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(_."_)                                                                )) y
  ('''nv''&zgeev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(                                    vberruR_mttmp_ >. vberrur_mttmp_))) y
  ('''vn''&zgeev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(vberruL_mttmp_ >. vberrul_mttmp_                                    ))) y
  ('''vv''&zgeev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(vberruL_mttmp_ >. vberrul_mttmp_ >. vberruR_mttmp_ >. vberrur_mttmp_))) y

  coerase < 'mttmp'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testheev
NB.
NB. Description:
NB.   Test DSYEV and ZHEEV (math/lapack2) by Hermitian
NB.   (symmetric) matrix
NB.
NB. Syntax:
NB.   testheev A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric)
NB.
NB. Formula:
NB.   berru := max(berru0,berru1)
NB.   ||A|| := max(||A||_1 , FP_SFMIN)
NB.   ||W|| := ||A - V * E * V^H||_1
NB.   if ||A|| > ||W|| then
NB.     berru0 := (||W|| / ||A||) / (FP_PREC * n)
NB.   elseif ||A|| < 1 then
NB.     berru0 := (min(||W|| , n * ||A||) / ||A||) / (FP_PREC * n)
NB.   else
NB.     berru0 := min(||W|| / ||A|| , n) / (FP_PREC * n)
NB.   endif
NB.   berru1 := min(||V * V^H - I||_1 , n) / (FP_PREC * n)
NB. where
NB.   A - matrix to eigen-decompose
NB.   W - diagonal matrix of eigenvalues
NB.   V - matrix of eigenvectors
NB.
NB. Notes:
NB. - models LAPACK's DSYT21 and ZHET21

testheev=: 3 : 0
  load_mttmp_ :: ] '~addons/math/mt/test/lapack2/dsyev.ijs'
  load_mttmp_ :: ] '~addons/math/mt/test/lapack2/zheev.ijs'

  rcond=. hecon1 y
  n=. # y
  nulp=. FP_PREC * n

  vberru0_mttmp_=. (- (0&{:: (] mp_mt_ (* ct_mt_)) 1&{::)) ((nulp %~ n <. %)`(nulp %~ (<. n&*) % ])@.(1 > ])`(nulp %~ %)@.< FP_SFMIN&>.)&norm1_mt_ [
  vberru1_mttmp_=. nulp %~ n <. norm1_mt_@(<: upddiag_mt_)@(mp_mt_ ct_mt_)@(1 {:: ])

  ('''nl''&dsyev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(_."_                            ))) y
  ('''nu''&dsyev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(_."_                            ))) y
  ('''vl''&dsyev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(vberru0_mttmp_ >. vberru1_mttmp_))) y
  ('''vu''&dsyev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(vberru0_mttmp_ >. vberru1_mttmp_))) y

  ('''nl''&zheev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(_."_                            ))) y
  ('''nu''&zheev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(_."_                            ))) y
  ('''vl''&zheev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(vberru0_mttmp_ >. vberru1_mttmp_))) y
  ('''vu''&zheev_mttmp_' tmonad (]`]`(rcond"_)`(_."_)`(vberru0_mttmp_ >. vberru1_mttmp_))) y

  coerase < 'mttmp'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testggev
NB.
NB. Description:
NB.   Test
NB.   - xGGEV (math/lapack2)
NB.   - ggevxxx (math/mt)
NB.   by pair of square matrices
NB.
NB. Syntax:
NB.   testggev AB
NB. where
NB.   AB - 2×n×n-brick

testggev=: 3 : 0
  load_mttmp_ :: ] '~addons/math/mt/test/lapack2/ggev.ijs'

  rcond=. <./ gecon1"2 y
  n=. c y

  vberrlL_mttmp_=.  mp_mt_~        "2` *   `normi_mt_`normitr_mt_ drgev_mt_
  vberrlR_mttmp_=. (mp_mt_  ct_mt_)"2`(*"1)`normi_mt_`normitr_mt_ drgev_mt_
  vberruL_mttmp_=. (mp_mt_~ ct_mt_)"2` *   `norm1_mt_`normitc_mt_ drgev_mt_
  vberruR_mttmp_=.  mp_mt_         "2`(*"1)`norm1_mt_`normitc_mt_ drgev_mt_

  ('''nn''&dggev_mttmp_' tmonad (]`]                                  `(rcond"_)`(_."_)`(_."_                                         ))) y
  ('''nv''&dggev_mttmp_' tmonad (]`((0&{:: ,: 1&{::) ;         3&{:: )`(rcond"_)`(_."_)`                         vberruR_mttmp_        )) y
  ('''vn''&dggev_mttmp_' tmonad (]`((0&{:: ,: 1&{::) ; 2&{::         )`(rcond"_)`(_."_)`  vberruL_mttmp_                               )) y
  ('''vv''&dggev_mttmp_' tmonad (]`((0&{:: ,: 1&{::) ; 2&{:: ; 3&{:: )`(rcond"_)`(_."_)`((vberruL_mttmp_ }:) >. (vberruR_mttmp_ 0 2&{)))) y

  ('''nn''&zggev_mttmp_' tmonad (]`]                                  `(rcond"_)`(_."_)`(_."_                                         ))) y
  ('''nv''&zggev_mttmp_' tmonad (]`((0&{:: ,: 1&{::) ;         3&{:: )`(rcond"_)`(_."_)`                         vberruR_mttmp_        )) y
  ('''vn''&zggev_mttmp_' tmonad (]`((0&{:: ,: 1&{::) ; 2&{::         )`(rcond"_)`(_."_)`  vberruL_mttmp_                               )) y
  ('''vv''&zggev_mttmp_' tmonad (]`((0&{:: ,: 1&{::) ; 2&{:: ; 3&{:: )`(rcond"_)`(_."_)`((vberruL_mttmp_ }:) >. (vberruR_mttmp_ 0 2&{)))) y

  ('ggevlnn'             tmonad (]`]                                  `(rcond"_)`(_."_)`(_."_                                         ))) y
  ('ggevlnv'             tmonad (]`]                                  `(rcond"_)`(_."_)`                         vberrlR_mttmp_        )) y
  ('ggevlvn'             tmonad (]`]                                  `(rcond"_)`(_."_)`  vberrlL_mttmp_                               )) y
  ('ggevlvv'             tmonad (]`(0&{:: ; (1 ; 0)&{:: ; (1 ; 1)&{::)`(rcond"_)`(_."_)`((vberrlL_mttmp_ }:) >. (vberrlR_mttmp_ 0 2&{)))) y

  ('ggevunn'             tmonad (]`]                                  `(rcond"_)`(_."_)`(_."_                                         ))) y
  ('ggevunv'             tmonad (]`]                                  `(rcond"_)`(_."_)`                         vberruR_mttmp_        )) y
  ('ggevuvn'             tmonad (]`]                                  `(rcond"_)`(_."_)`  vberruL_mttmp_                               )) y
  ('ggevuvv'             tmonad (]`(0&{:: ; (1 ; 0)&{:: ; (1 ; 1)&{::)`(rcond"_)`(_."_)`((vberruL_mttmp_ }:) >. (vberruR_mttmp_ 0 2&{)))) y

  coerase < 'mttmp'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testev
NB.
NB. Description:
NB.   Adv. to make verb to test ggevxxx by matrices of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testev
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
NB.     ?@$&0 testev_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testev_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testev_mt_ 150 150

testev=: 1 : 'EMPTY [ (testggev_mt_@u@(2&,) [ testheev_mt_@(u hemat_mt_) [ testgeev_mt_@u)^:(=/)'
