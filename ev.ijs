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
NB. Concepts
NB.
NB. I have a dream (2010-11-02):
NB.   ggev=: gizmo @ tgevc @ hgeqz @ gghrd @ qr &. bal &. scl
NB. Is it possible?

NB. =========================================================
NB. Configuration

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
NB. ggevi
NB.
NB. Description:
NB.   Adv. to make monad to calculate initial arguments for
NB.   ggevxxx
NB.
NB. Syntax:
NB.   'abnrmio ABupd plr hs'=. (ggbalp ggevi) AB
NB. where
NB.   ggbalp  - monad to permute matrix pair (A,B) to isolate
NB.             eigenvalues, is either ggballp or ggbalup, is
NB.             called as:
NB.               'CD plr hs'=. ggbalp AB
NB.   AB      - 2×n×n-matrix, matrix pair (A,B)
NB.   abnrmio -:abnrm ,. abio
NB.   abnrm   - 2-vector, norms of A and B
NB.   abio    - 2-vector of integers, defines both necessity
NB.             and value of scaling for A and B
NB.   ABupd   - an updated AB, contains scaled and permuted A
NB.             and B
NB.   plr     - 2×n-matrix of integers, permutations of A and
NB.             B, produced by ggbalp
NB.   hs      - 2-vector of integers, defines submatrices
NB.             position, produced by ggbalp

ggevi=: 1 : '(,.(0,(EVSMLNUM_mt_*1-FP_EPS_mt_),EVBIGNUM_mt_)&I.)@:(normm_mt_"2) ([ ; u@:(scl_mt_^:((,{&EVSCL_mt_)/@[`({&0 1 0 1@{:@[)`])"1 2)) ]'

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
  iso=. < 0 1 ; ;~ liso4dhs hs
  y=. (gghrdlnn~ 0 , c)&.(iso&{) ((unmlqrc~ ,: trl@:(}:"1)@]) gelqf)/&.(iso&{) y
  e1e2=. hs hgezqenn y
  e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
)

ggevlnv=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi y
  y=. (0 1;(<i.{.hs);liso4dhs hs) ((unmlqrc~,(trl@:(}:"1),:unglq)@])gelqf)/@({~<)~`(((<0 1 2))<@(0})[)`((, idmat@c)@])} y
  y=. (gghrdlnv~0,c) y
  y=. hs hgezqsnv y
  e1e2=. 2 {. diag"2 y
  if. isnan < e1e2 do.
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
  iso=. < 0 1 ; (< i. {. hs) ; liso4dhs hs
  y=. ((unmlqrc~ ,: trl@:(}:"1)@]) gelqf)/&.(iso&{) y
  y=. (((0,]) gghrdlvn (,idmat)) c) y
  y=. hs hgezqsvn y
  e1e2=. 2 {. diag"2 y
  if. isnan < e1e2 do.
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
  y=. (0 1;(<i.{.hs);liso4dhs hs) ((unmlqrc~,(trl@:(}:"1),:unglq)@])gelqf)/@({~<)~`(((<0 1 3))<@(0})[)`((, ,:~@idmat@c)@])} y
  y=. (gghrdlvv~0,c) y
  y=. hs hgezqsvv y
  e1e2=. 2 {. diag"2 y
  if. isnan < e1e2 do.
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
  iso=. < 0 1 ; ;~ liso4dhs hs
  y=. (gghrdunn~ 0 , c)&.(iso&{) ((unmqrlc~ ,: tru@}:@]) geqrf)/&.(iso&{) y
  e1e2=. hs hgeqzenn y
  e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
)

ggevuvn=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi y
  y=. (0 1;(liso4dhs hs);<<i.{.hs) ((unmqrlc~,(tru@}:,:ungqr)@])geqrf)/@({~<)~`(((<0 1 2))<@(0})[)`((, idmat@c)@])} y
  y=. (gghrduvn~0,c) y
  y=. hs hgeqzsvn y
  e1e2=. 2 {. diag"2 y
  if. isnan < e1e2 do.
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
  iso=. < 0 1 ; (liso4dhs hs) ; < < i. {. hs
  y=. ((unmqrlc~ ,: tru@}:@]) geqrf)/&.(iso&{) y
  y=. (((0,]) gghrdunv (,idmat)) c) y
  y=. hs hgeqzsnv y
  e1e2=. 2 {. diag"2 y
  if. isnan < e1e2 do.
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
  y=. (0 1;(liso4dhs hs);<<i.{.hs) ((unmqrlc~,(tru@}:,:ungqr)@])geqrf)/@({~<)~`(((<0 1 2))<@(0})[)`((, ,:~@idmat@c)@])} y
  y=. (gghrduvv~0,c) y
  y=. hs hgeqzsvv y
  e1e2=. 2 {. diag"2 y
  if. isnan < e1e2 do.
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
NB.   Test xGEEV (math/lapack2 addon) by square matrix
NB.
NB. Syntax:
NB.   log=. testgeev A
NB. where
NB.   A   - n×n-matrix
NB.   log - 6-vector of boxes, test log

testgeev=: 3 : 0
  _1 cocreate < 'mttmp'
  load_mttmp_ 'math/mt/external/lapack2/geev'

  rcondu=. gecon1 y

  'norml normr'=. (normi , norm1) y

  log=.          ('''nn''&dgeev_mttmp_' tmonad (]      `]`(rcondu"_)`nan`nan                               )) y
  log=. log lcat ('''nv''&dgeev_mttmp_' tmonad ((0&{::)`]`(rcondu"_)`nan`(                  t22r >. drvevr))) y ; _.    ; normr
  log=. log lcat ('''vn''&dgeev_mttmp_' tmonad ((0&{::)`]`(rcondu"_)`nan`(t22l >. drvevl                  ))) y ; norml
  log=. log lcat ('''vv''&dgeev_mttmp_' tmonad ((0&{::)`]`(rcondu"_)`nan`(t22l >. drvevl >. t22r >. drvevr))) y ; norml ; normr

  log=. log lcat ('''nn''&zgeev_mttmp_' tmonad (]      `]`(rcondu"_)`nan`nan                               )) y
  log=. log lcat ('''nv''&zgeev_mttmp_' tmonad ((0&{::)`]`(rcondu"_)`nan`(                  t22r >. drvevr))) y ; _.    ; normr
  log=. log lcat ('''vn''&zgeev_mttmp_' tmonad ((0&{::)`]`(rcondu"_)`nan`(t22l >. drvevl                  ))) y ; norml
  log=. log lcat ('''vv''&zgeev_mttmp_' tmonad ((0&{::)`]`(rcondu"_)`nan`(t22l >. drvevl >. t22r >. drvevr))) y ; norml ; normr

  coerase < 'mttmp'

  log
)

NB. ---------------------------------------------------------
NB. testheev
NB.
NB. Description:
NB.   Test DSYEV and ZHEEV (math/lapack2 addon) by Hermitian
NB.   (symmetric) matrix
NB.
NB. Syntax:
NB.   log=. testheev A
NB. where
NB.   A   - n×n-matrix, the Hermitian (symmetric)
NB.   log - 6-vector of boxes, test log

testheev=: 3 : 0
  _1 cocreate < 'mttmp'
  load_mttmp_ 'math/mt/external/lapack2/dsyev'
  load_mttmp_ 'math/mt/external/lapack2/zheev'

  rcondl=. heconi y

  norml=. normi y

  log=.          ('''nl''&dsyev_mttmp_' tmonad ((0&{::)`]`(rcondl"_)`nan`nan )) y ; norml
  log=. log lcat ('''nu''&dsyev_mttmp_' tmonad ((0&{::)`]`(rcondl"_)`nan`nan )) y ; norml
  log=. log lcat ('''vl''&dsyev_mttmp_' tmonad ((0&{::)`]`(rcondl"_)`nan`t211)) y ; norml
  log=. log lcat ('''vu''&dsyev_mttmp_' tmonad ((0&{::)`]`(rcondl"_)`nan`t211)) y ; norml

  log=. log lcat ('''nl''&zheev_mttmp_' tmonad ((0&{::)`]`(rcondl"_)`nan`nan )) y ; norml
  log=. log lcat ('''nu''&zheev_mttmp_' tmonad ((0&{::)`]`(rcondl"_)`nan`nan )) y ; norml
  log=. log lcat ('''vl''&zheev_mttmp_' tmonad ((0&{::)`]`(rcondl"_)`nan`t211)) y ; norml
  log=. log lcat ('''vu''&zheev_mttmp_' tmonad ((0&{::)`]`(rcondl"_)`nan`t211)) y ; norml

  coerase < 'mttmp'

  log
)

NB. ---------------------------------------------------------
NB. testggev
NB.
NB. Description:
NB.   Test:
NB.   - xGGEV (math/lapack2 addon)
NB.   - ggevxxx (math/mt addon)
NB.   by pair of square matrices
NB.
NB. Syntax:
NB.   log=. testggev AB
NB. where
NB.   AB  - 2×n×n-brick
NB.   log - 6-vector of boxes, test log

testggev=: 3 : 0
  _1 cocreate < 'mttmp'
  load_mttmp_ 'math/mt/external/lapack2/ggev'

  'rcondl rcondu'=. <./ (geconi , gecon1)"2 y

  'norml normu'=. |: (normi , norm1)"2 y

  vberrlL=:  mp~    "2` *   `normi`normitr drgev
  vberrlR=: (mp  ct)"2`(*"1)`normi`normitr drgev
  vberruL=: (mp~ ct)"2` *   `norm1`normitc drgev
  vberruR=:  mp     "2`(*"1)`norm1`normitc drgev

  log=.          ('''nn''&dggev_mttmp_' tmonad ((0&{::)`]                                  `(rcondu"_)`nan`nan                              )) y ; normu
  log=. log lcat ('''nv''&dggev_mttmp_' tmonad ((0&{::)`((0&{:: ,: 1&{::) ;         3&{:: )`(rcondu"_)`nan`                  vberruR        )) y ; normu
  log=. log lcat ('''vn''&dggev_mttmp_' tmonad ((0&{::)`((0&{:: ,: 1&{::) ; 2&{::         )`(rcondu"_)`nan`  vberruL                        )) y ; normu
  log=. log lcat ('''vv''&dggev_mttmp_' tmonad ((0&{::)`((0&{:: ,: 1&{::) ; 2&{:: ; 3&{:: )`(rcondu"_)`nan`((vberruL }:) >. (vberruR 0 2&{)))) y ; normu

  log=. log lcat ('''nn''&zggev_mttmp_' tmonad ((0&{::)`]                                  `(rcondu"_)`nan`nan                              )) y ; normu
  log=. log lcat ('''nv''&zggev_mttmp_' tmonad ((0&{::)`((0&{:: ,: 1&{::) ;         3&{:: )`(rcondu"_)`nan`                  vberruR        )) y ; normu
  log=. log lcat ('''vn''&zggev_mttmp_' tmonad ((0&{::)`((0&{:: ,: 1&{::) ; 2&{::         )`(rcondu"_)`nan`  vberruL                        )) y ; normu
  log=. log lcat ('''vv''&zggev_mttmp_' tmonad ((0&{::)`((0&{:: ,: 1&{::) ; 2&{:: ; 3&{:: )`(rcondu"_)`nan`((vberruL }:) >. (vberruR 0 2&{)))) y ; normu

  log=. log lcat ('ggevlnn'             tmonad ((0&{::)`]                                  `(rcondl"_)`nan`nan                              )) y ; norml
  log=. log lcat ('ggevlnv'             tmonad ((0&{::)`]                                  `(rcondl"_)`nan`                  vberrlR        )) y ; norml
  log=. log lcat ('ggevlvn'             tmonad ((0&{::)`]                                  `(rcondl"_)`nan`  vberrlL                        )) y ; norml
  log=. log lcat ('ggevlvv'             tmonad ((0&{::)`(0&{:: ; (1 ; 0)&{:: ; (1 ; 1)&{::)`(rcondl"_)`nan`((vberrlL }:) >. (vberrlR 0 2&{)))) y ; norml

  log=. log lcat ('ggevunn'             tmonad ((0&{::)`]                                  `(rcondu"_)`nan`nan                              )) y ; normu
  log=. log lcat ('ggevunv'             tmonad ((0&{::)`]                                  `(rcondu"_)`nan`                  vberruR        )) y ; normu
  log=. log lcat ('ggevuvn'             tmonad ((0&{::)`]                                  `(rcondu"_)`nan`  vberruL                        )) y ; normu
  log=. log lcat ('ggevuvv'             tmonad ((0&{::)`(0&{:: ; (1 ; 0)&{:: ; (1 ; 1)&{::)`(rcondu"_)`nan`((vberruL }:) >. (vberruR 0 2&{)))) y ; normu

  coerase < 'mttmp'
  erase 'vberrlL vberrlR vberruL vberruR'

  log
)

NB. ---------------------------------------------------------
NB. testev
NB.
NB. Description:
NB.   Adv. to make verb to test ggevxxx by matrices of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testev) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log
NB.
NB. Application:
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testev_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testev_mt_ 150 150
NB. - test by random square complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testev_mt_ 150 150

testev=: 1 : 'nolog_mt_`(lcat_mt_@(testgeev_mt_@u`(testheev_mt_@(u hemat_mt_))`(testggev_mt_@u@(2&,))`:0))@.(=/)'
