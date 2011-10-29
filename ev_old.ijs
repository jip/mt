NB. Eigenvalues and eigenvectors
NB.
NB. ggevxxx   Eigenvalues and, optionally, eigenvectors of
NB.           pair of matrices
NB.
NB. testgeev_old  Test geevxxx by general matrix given
NB. testheev_old  Test heevxx by Hermitian (symmetric) matrix
NB.           given
NB. testggev_old  Test ggevxxx by general matrices given
NB. testev_old    Adv. to make verb to test xxevxxx by matrices
NB.           of generator and shape given
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
NB. Constants

NB. Scaling limits
EVBIGNUM=: % EVSMLNUM=: (%: FP_SFMIN) % FP_PREC

NB. Scaling factors
EVSCL=: 1 , EVSMLNUM , 1 , EVBIGNUM

NB. ---------------------------------------------------------
NB. ggevi_old
NB.
NB. Description:
NB.   Adv. to make verb to calculate initial parameters for
NB.   ggevxxx
NB.
NB. Syntax:
NB.   vapp=. ggbalp ggevi_old
NB. where
NB.   ggbalp  - monadic verb to permute matrix pair (A,B) to
NB.             isolate eigenvalues, is either ggballp or
NB.             ggbalup, is called as:
NB.               'CD plr hs'=. ggbalp AB
NB.   vapp    - monadic verb to calculate initial parameters
NB.             for ggevxxx, is called as:
NB.               'abnrmio ABupd plr hs'=. vapp AB
NB.   AB      - 2×n×n-matrix, matrix pair (A,B)
NB.   abnrmio -: abnrm ,. abio
NB.   abnrm   - 2-vector, norms of A and B
NB.   abio    - 2-vector of integers, defines both necessity
NB.             and value of scaling for A and B
NB.   ABupd   - 2×n×n-matrix, scaled and permuted A and B
NB.   plr     - 2×n-matrix of integers, permutations of A and
NB.             B, produced by ggbalp
NB.   hs      - 2-vector of integers, defines submatrices
NB.             position, produced by ggbalp

ggevi_old=: 1 : '(,.(0,(EVSMLNUM_mt_*1-FP_EPS_mt_),EVBIGNUM_mt_)&I.)@:(>./@,"2)@:| ([ ; u@:(scl_mt_^:((,{&EVSCL_mt_)/@[`({&0 1 0 1@{:@[)`])"1 2)) ]'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. ggevlnn_old
NB. ggevlnv_old
NB. ggevlvn_old
NB. ggevlvv_old
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
NB.   e1e2=.      ggevlnn_old AB
NB.   'e1e2 R'=.  ggevlnv_old AB
NB.   'e1e2 L'=.  ggevlvn_old AB
NB.   'e1e2 LR'=. ggevlvv_old AB
NB. where
NB.   AB   - 2×n×n-matrix, matrix pair (A,B):
NB.            AB -: A ,: B
NB.   e1e2 - 2×n-matrix of eigenvalues e1 and e2:
NB.            e1e2 -: e1 ,: e2
NB.   L    - n×n-matrix. left eigenvectors (rows)
NB.   R    - n×n-matrix. right eigenvectors (rows)
NB.   LR   - 2×n×n-matrix. left and right eigenvectors:
NB.            LR -: L ,: R
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (ggevlnn_old -:                           +@:ggevunn_old@:(ct"2) ) A ,: B
NB.   (ggevlnn_old -:&(/:~@(%/))                   ggevunn_old         ) A ,: B
NB.   (ggevlnv_old -: (0 1 (+&.>)`(ct&.>)       ag ggevuvn_old@:(ct"2))) A ,: B
NB.   (ggevlvn_old -: (0 1 (+&.>)`(ct&.>)       ag ggevunv_old@:(ct"2))) A ,: B
NB.   (ggevlvv_old -: (0 1 (+&.>)`(ct"2@:|.&.>) ag ggevuvv_old@:(ct"2))) A ,: B
NB.   (E2 mp L mp A) -: (E1 mp L mp B)
NB.   (A mp (ct R) mp E2) -: (B mp (ct R) mp E1)
NB. where
NB.   A - n×n-matrix, general
NB.   B - n×n-matrix, general
NB.   'e1e2 LR'=. ggevlvv_old A ,: B
NB.   'E1 E2'=. diagmat"1 e1e2
NB.   'L R'=. LR
NB.
NB. Application:
NB. - simulate LAPACK's xGEEV('N','N'):
NB.     NB. e=. geevlnn A
NB.     geevlnn=: {.@ggevlnn_old@(,:(idmat@c))
NB. - simulate LAPACK's xGEEV('N','V') (see notes):
NB.     NB. 'e R'=. geevlnv A
NB.     geevlnv=: 0 1 ({.&.>)`(((**@+@((i.>./)"1@sorim{"0 1]))%norms"1)      &.>)ag ggevlnv_old@(,:(idmat@c))
NB. - simulate LAPACK's xGEEV('V','N') (see notes):
NB.     NB. 'e L'=. geevlvn A
NB.     geevlvn=: 0 1 ({.&.>)`(((**@+@((i.>./)"1@sorim{"0 1]))%norms"1)      &.>)ag ggevlvn_old@(,:(idmat@c))
NB. - simulate LAPACK's xGEEV('V','V') (see notes):
NB.     NB. 'e LR'=. geevlvv A
NB.     geevlvv=: 0 1 ({.&.>)`(((**@+@((i.>./)"1@sorim{"0 1]))%norms"1)    "2&.>)ag ggevlvv_old@(,:(idmat@c))
NB. - simulate LAPACK's xHEEV('N'):
NB.     NB. e=. heevln A
NB.     heevln=: 9 o.{.@ggevlnn_old@(,:(idmat@c))
NB. - simulate LAPACK's xHEEV('V') (see notes):
NB.     NB. 'e V'=. heevlv A
NB.     heevlv=: 0 1 ((9 o.{.)&.>)`((%  %:@diag@(mp ct))&.>)ag ggevlvn_old@(,:(idmat@c))
NB.
NB. Notes:
NB. - eigenvectors from LAPACK's xGEEV are normalized to have
NB.   Euclidean norm equal to 1 and largest component real
NB. - eigenvectors from LAPACK's xHEEV are orthonormal

ggevlnn_old=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi_old y
  y=. (<0 1;;~dhs2lios hs) ([ ((gghrdlnn_old~0,c) upd) ((unmlqrc~,:trl@:(}:"1)@])gelqf)/@{`[`] }) y
  e1e2=. hs hgezqenn_old y
  e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
)

ggevlnv_old=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi_old y
  y=. (0 1;(<i.{.hs);dhs2lios hs) ((unmlqrc~,(trl@:(}:"1),:unglq)@])gelqf)/@({~<)~`((<0 1 2)<@(0})[)`((, idmat@c)@]) } y
  y=. (gghrdlnv_old~0,c) y
  y=. hs hgezqsnv_old y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2;_.$~2#c y
  else.
    y=. tgevclrb_old y
    y=. gebaklp y ; {: plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr) y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevlvn_old=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi_old y
  y=. (<0 1;(<i.{.hs);dhs2lios hs) ((unmlqrc~,:trl@:(}:"1)@])gelqf)/@{`[`] } y
  y=. (((0,]) gghrdlvn_old (,idmat)) c) y
  y=. hs hgezqsvn_old y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2;_.$~2#c y
  else.
    y=. tgevcllb_old y
    y=. gebaklp y ; {. plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr) y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevlvv_old=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi_old y
  y=. (0 1;(<i.{.hs);dhs2lios hs) ((unmlqrc~,(trl@:(}:"1),:unglq)@])gelqf)/@({~<)~`((<0 1 3)<@(0})[)`((, ,:~@idmat@c)@]) } y
  y=. (gghrdlvv_old~0,c) y
  y=. hs hgezqsvv_old y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2;_.$~2,2#c y
  else.
    y=. tgevclbb_old y
    y=. y gebaklp@;"2 1 plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr)"2 y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

NB. ---------------------------------------------------------
NB. ggevunn_old
NB. ggevunv_old
NB. ggevuvn_old
NB. ggevuvv_old
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
NB.   e1e2=.      ggevunn_old AB
NB.   'e1e2 R'=.  ggevunv_old AB
NB.   'e1e2 L'=.  ggevuvn_old AB
NB.   'e1e2 LR'=. ggevuvv_old AB
NB. where
NB.   AB   - 2×n×n-matrix, matrix pair (A,B):
NB.            AB -: A ,: B
NB.   e1e2 - 2×n-matrix of eigenvalues e1 and e2:
NB.            e1e2 -: e1 ,: e2
NB.   L    - n×n-matrix. left eigenvectors (columns)
NB.   R    - n×n-matrix. right eigenvectors (columns)
NB.   LR   - 2×n×n-matrix. left and right eigenvectors:
NB.            LR -: L ,: R
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (ggevunn_old -:                           +@:ggevlnn_old@:(ct"2) ) A ,: B
NB.   (ggevunn_old -:&(/:~@(%/))                   ggevlnn_old         ) A ,: B
NB.   (ggevunv_old -: (0 1 (+&.>)`(ct&.>)       ag ggevlvn_old@:(ct"2))) A ,: B
NB.   (ggevuvn_old -: (0 1 (+&.>)`(ct&.>)       ag ggevlnv_old@:(ct"2))) A ,: B
NB.   (ggevuvv_old -: (0 1 (+&.>)`(ct"2@:|.&.>) ag ggevlvv_old@:(ct"2))) A ,: B
NB.   (E2 mp (ct L) mp A) -: (E1 mp (ct L) mp B)
NB.   (A mp R mp E2) -: (B mp R mp E1)
NB. where
NB.   A - n×n-matrix, general
NB.   B - n×n-matrix, general
NB.   'e1e2 LR'=. ggevuvv_old A ,: B
NB.   'E1 E2'=. diagmat"1 e1e2
NB.   'L R'=. LR
NB.
NB. Application:
NB. - simulate LAPACK's xGEEV('N','N'):
NB.     NB. e=. geevunn A
NB.     geevunn=: {.@ggevunn_old@(,:(idmat@c))
NB. - simulate LAPACK's xGEEV('N','V') (see notes):
NB.     NB. 'e R'=. geevunv A
NB.     geevunv=: 0 1 ({.&.>)`(((**@+@((i.>./)"1@sorim{"0 1]))%norms"1)&.|:  &.>)ag ggevunv_old@(,:(idmat@c))
NB. - simulate LAPACK's xGEEV('V','N') (see notes):
NB.     NB. 'e L'=. geevuvn A
NB.     geevuvn=: 0 1 ({.&.>)`(((**@+@((i.>./)"1@sorim{"0 1]))%norms"1)&.|:  &.>)ag ggevuvn_old@(,:(idmat@c))
NB. - simulate LAPACK's xGEEV('V','V') (see notes):
NB.     NB. 'e LR'=. geevuvv A
NB.     geevuvv=: 0 1 ({.&.>)`(((**@+@((i.>./)"1@sorim{"0 1]))%norms"1)&.|:"2&.>)ag ggevuvv_old@(,:(idmat@c))
NB. - simulate LAPACK's xHEEV('N'):
NB.     NB. e=. heevun A
NB.     heevun=: 9 o.{.@ggevunn_old@(,:(idmat@c))
NB. - simulate LAPACK's xHEEV('V') (see notes):
NB.     NB. 'e V'=. heevuv A
NB.     heevuv=: 0 1 ((9 o.{.)&.>)`((%"1%:@diag@(mp~ct))&.>)ag ggevunv_old@(,:(idmat@c))
NB.
NB. Notes:
NB. - ggevunn_old models LAPACK's xGGEV('N','N')
NB. - ggevunv_old models LAPACK's xGGEV('N','V')
NB. - ggevuvn_old models LAPACK's xGGEV('V','N')
NB. - ggevuvv_old models LAPACK's xGGEV('V','V')
NB. - eigenvectors from LAPACK's xGEEV are normalized to have
NB.   Euclidean norm equal to 1 and largest component real
NB. - eigenvectors from LAPACK's xHEEV are orthonormal

ggevunn_old=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi_old y
  y=. (<0 1;;~dhs2lios hs) ([ ((gghrdunn_old~0,c) upd) ((unmqrlc~,:tru@}:@])geqrf)/@{`[`] }) y
  e1e2=. hs hgeqzenn_old y
  e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
)

ggevuvn_old=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi_old y
  y=. (0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,(tru@}:,:ungqr)@])geqrf)/@({~<)~`((<0 1 2)<@(0})[)`((, idmat@c)@]) } y
  y=. (gghrduvn_old~0,c) y
  y=. hs hgeqzsvn_old y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2;_.$~2#c y
  else.
    y=. tgevculb_old y
    y=. gebakup y ; {. plr
    y=. (%"1 (EVSMLNUM&>`(,:&1))}@:normitc) y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevunv_old=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi_old y
  y=. (<0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,:tru@}:@])geqrf)/@{`[`] } y
  y=. (((0,]) gghrdunv_old (,idmat)) c) y
  y=. hs hgeqzsnv_old y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2;_.$~2#c y
  else.
    y=. tgevcurb_old y
    y=. gebakup y ; {: plr
    y=. (%"1 (EVSMLNUM&>`(,:&1))}@:normitc) y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevuvv_old=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi_old y
  y=. (0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,(tru@}:,:ungqr)@])geqrf)/@({~<)~`((<0 1 2)<@(0})[)`((, ,:~@idmat@c)@]) } y
  y=. (gghrduvv_old~0,c) y
  y=. hs hgeqzsvv_old y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2;_.$~2,2#c y
  else.
    y=. tgevcubb_old y
    y=. y gebakup@;"2 1 plr
    y=. (%"1 (EVSMLNUM&>`(,:&1))}@:normitc)"2 y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgeev_old
NB.
NB. Description:
NB.   Test geev (math/lapack) by general matrix given
NB.
NB. Syntax:
NB.   testgeev_old A
NB. where
NB.   A - n×n-matrix
NB.
NB. Formula:
NB.   berr := max(berrL,berrR,berrl,berrr)
NB. where
NB.   ||A|| := max(||A||_1 , FP_SFMIN)
NB.   L - matrix of left eigenvectors
NB.   R - matrix of right eigenvectors
NB.   E - diagonal matrix of eigenvalues
NB.   - geev_jlapack_:
NB.       berrL := min((||L^H * A - E * L^H||_1 / max(||L||_1 , FP_PREC)) / ||A|| , 1) / FP_PREC
NB.       berrR := min((||A   * R - R * E  ||_1 / max(||R||_1 , FP_PREC)) / ||A|| , 1) / FP_PREC
NB.       if max(|Re(X[i,j])|) / max(|X[i,j]|) < 1 - 2*FP_PREC then
NB.         berrx := 1/FP_PREC
NB.       else
NB.         berrx := max(min(1/FP_PREC , | ||X[:,j]||_E - 1 |) / FP_PREC
NB.       endif
NB.       where berrx is either berrl or berrr, and corresponding X is either L or R

testgeev_old=: 3 : 0
  require :: ] '~addons/math/lapack/lapack.ijs'
  need_jlapack_ :: ] 'geev'

  rcond=. gecon1 y

  vberruL=: FP_PREC%~1<.norm1@[%~(((mp~ ct)~(1;0)&{::)normi@:-((0+@{::])*"1(1;0){::])) % FP_PREC>.norm1@((1;0){::])
  vberruR=: FP_PREC%~1<.norm1@[%~(( mp      (1;1)&{::)norm1@:-((0  {::])*"1(1;1){::])) % FP_PREC>.norm1@((1;1){::])
  vberrux=: 1(((%FP_PREC)>./"1@:<.FP_PREC%~|@:<:@:(norms"1@|:"2))([`((%FP_PREC)"_)@.])"0((1-+:FP_PREC)>(|@(9&o.))%&:(>./"1)|)@:(,"2))@{::]
  vberruvv=: >./@(vberruL,vberruR,vberrux)  NB. STUDYME: verb is invisible from base locale when erase below is removed (???)

  ('geev_jlapack_' tmonad (]`(1&({::) ; (0,:2)&({::))`(rcond"_)`(_."_)`vberruvv)) y

  erase 'vberruL vberruR vberrux vberruvv'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testheev_old
NB.
NB. Description:
NB.   Test heev (math/lapack) by Hermitian (symmetric) matrix
NB.   given
NB.
NB. Syntax:
NB.   testheev_old A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric)
NB.
NB. Formula:
NB.   berr := max(berr0,berr1)
NB. where
NB.   ||A|| := max(||A||_1 , FP_SFMIN)
NB.   V - matrix of eigenvectors
NB.   E - diagonal matrix of eigenvalues
NB.   - heev_jlapack_:
NB.       wnorm := ||A - V * E * V^H||_1
NB.       if ||A|| > wnorm then
NB.         berr0 := (wnorm / ||A||) / (FP_PREC * n)
NB.       elseif ||A|| < 1 then
NB.         berr0 := (min(wnorm , n * ||A||) / ||A||) / (FP_PREC * n)
NB.       else
NB.         berr0 := min(wnorm / ||A|| , n) / (FP_PREC * n)
NB.       endif
NB.       berr1 := min(||V * V^H - I||_1 , n) / (FP_PREC * n)

testheev_old=: 3 : 0
  require :: ] '~addons/math/lapack/lapack.ijs'
  need_jlapack_ :: ] 'heev'

  rcond=. hecon1 y

  vberru0=: [((((%~{.)<.1{])`(((0{])<.(*{:))%[)@.(1>[)`(%~{.)@.(>{.)%FP_PREC*1{])~(FP_SFMIN>.{.))~&(norm1,#)(-(]mp(*ct))&>/)
  vberru1=: 1((#<.norm1@(<:upddiag)@(mp ct))%(FP_PREC*#))@{::]
  vberruv=: vberru0>.vberru1

  ('heev_jlapack_' tmonad (]`]`(rcond"_)`(_."_)`vberruv)) y

  erase 'vberru0 vberru1 vberruv'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testggev_old
NB.
NB. Description:
NB.   Test ggevxxx by general matrices given
NB.
NB. Syntax:
NB.   testggev_old AB
NB. where
NB.   AB - 2×n×n-report (A,:B)
NB.
NB. Formula:
NB.   berr := max(berr0,berr1)
NB. where
NB.   ||M|| := max(||M||_1 , FP_SFMIN)
NB.   ||v|| := max(|Re(v(i))|+|Im(v(i))|)
NB.   α(i)  - i-th eigenvalue, also i-th element on S
NB.           diagonal
NB.   β(i)  - i-th eigenvalue, also i-th element on P
NB.           diagonal
NB.   l(i)  - i-th left eigenvector
NB.   lb(i) - i-th back transformed left eigenvector
NB.   r(i)  - i-th right eigenvector
NB.   rb(i) - i-th back transformed right eigenvector
NB.   - ggevlvn_old:
NB.       berr0 := max(||l(i) * (β(i)*A - α(i)*B)  || / (FP_PREC * max(|| β(i)*A   ||,|| α(i)*B   ||)))
NB.       berr1 := max(| ||l(i)|| - 1 |) / (FP_PREC * n)
NB.   - ggevlnv_old:
NB.       berr0 := max(||r(i) * (β(i)*A - α(i)*B)^H|| / (FP_PREC * max(||(β(i)*A)^H||,||(α(i)*B)^H||)))
NB.       berr1 := max(| ||r(i)|| - 1 |) / (FP_PREC * n)
NB.   - ggevlvv_old:
NB.       berr0 := berr(ggevlvn_old)
NB.       berr1 := berr(ggevlnv_old)
NB.   - ggevuvn_old:
NB.       berr0 := max(||(β(i)*A - α(i)*B)^H * l(i)|| / (FP_PREC * max(||(β(i)*A)^H||,||(α(i)*B)^H||)))
NB.       berr1 := max(| ||l(i)|| - 1 |) / (FP_PREC * n)
NB.   - ggevunv_old:
NB.       berr0 := max(||(β(i)*A - α(i)*B)   * r(i)|| / (FP_PREC * max(|| β(i)*A   ||,|| α(i)*B   ||)))
NB.       berr1 := max(| ||r(i)|| - 1 |) / (FP_PREC * n)
NB.   - ggevuvv_old:
NB.       berr0 := berr(ggevuvn_old)
NB.       berr1 := berr(ggevunv_old)
NB.
NB. Notes:
NB. - berrxxx are non-iterative and are require O(N^3) RAM

testggev_old=: 3 : 0
  vberrlvn=:  (     normir@:((*"_ 1|:@|.@(0&{::))((  norm1r@:((mp"1 2 (-/"3))~     )                                           )      % FP_PREC*(FP_SFMIN>.    (>./"1)@:( norm1       "2)@[))(1 {:: ])))>.(     normir@:<:@:normitr%FP_PREC*c)@(1 {:: ])
  vberrlnv=:  (     normir@:((*"_ 1|:@|.@(0&{::))((                                   norm1r@:((mp"2 1~(-/"3))~+    )          )      % FP_PREC*(FP_SFMIN>.    (>./"1)@:(       normi "2)@[))(1 {:: ])))>.(     normir@:<:@:normitr%FP_PREC*c)@(1 {:: ])
  vberrlvv=:  (>./@:normir@:((*"_ 1|:@|.@(0&{::))((((norm1r@:( mp"1 2        ~   {.),:norm1r@:((mp"2 1       ) + @{:))~(-/"3))~)(>./@:%)FP_PREC*(FP_SFMIN>.|:@:(>./"1)@:((norm1,normi)"2)@[))(1 {:: ])))>.(>./@:normir@:<:@:normitr%FP_PREC*c)@(1 {:: ])

  vberruvn=:  (     normir@:((*"_ 1|:@|.@(0&{::))((  norm1r@:((mp"1 2 (-/"3))~ct   )                                           )      % FP_PREC*(FP_SFMIN>.    (>./"1)@:( normi       "2)@[))(1 {:: ])))>.(     normir@:<:@:normitc%FP_PREC*c)@(1 {:: ])
  vberrunv=:  (     normir@:((*"_ 1|:@|.@(0&{::))((                                   norm1r@:((mp"2 1~(-/"3))~|:   )          )      % FP_PREC*(FP_SFMIN>.    (>./"1)@:(       norm1 "2)@[))(1 {:: ])))>.(     normir@:<:@:normitc%FP_PREC*c)@(1 {:: ])
  vberruvv=:  (>./@:normir@:((*"_ 1|:@|.@(0&{::))((((norm1r@:( mp"1 2        ~ct@{.),:norm1r@:((mp"2 1       ) |:@{:))~(-/"3))~)(>./@:%)FP_PREC*(FP_SFMIN>.|:@:(>./"1)@:((normi,norm1)"2)@[))(1 {:: ])))>.(>./@:normir@:<:@:normitc%FP_PREC*c)@(1 {:: ])

  rcond=. <./ gecon1"2 y

  ('ggevlnn_old' tmonad (]`]`(rcond"_)`(_."_)`(_."_)  )) y
  ('ggevlnv_old' tmonad (]`]`(rcond"_)`(_."_)`vberrlnv)) y
  ('ggevlvn_old' tmonad (]`]`(rcond"_)`(_."_)`vberrlvn)) y
  ('ggevlvv_old' tmonad (]`]`(rcond"_)`(_."_)`vberrlvv)) y
  ('ggevunn_old' tmonad (]`]`(rcond"_)`(_."_)`(_."_)  )) y
  ('ggevunv_old' tmonad (]`]`(rcond"_)`(_."_)`vberrunv)) y
  ('ggevuvn_old' tmonad (]`]`(rcond"_)`(_."_)`vberruvn)) y
  ('ggevuvv_old' tmonad (]`]`(rcond"_)`(_."_)`vberruvv)) y

  erase 'vberrlnv vberrlvn vberrlvv vberrunv vberruvn vberruvv'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testev_old
NB.
NB. Description:
NB.   Adv. to make verb to test ggevxxx by matrices of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testev_old
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

testev_old=: 1 : 'EMPTY_mt_ [ (testggev_old_mt_@u@(2&,) [ testheev_old_mt_@(u hemat_mt_) [ testgeev_old_mt_@u)^:(=/)'
