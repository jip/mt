NB. Eigenvalues and, optionally, eigenvectors
NB.
NB. ggevxxx   Eigenvalues and, optionally, eigenvectors of
NB.           pair of matrices
NB.
NB. testggev  Test ggevxxx by general matrices given
NB. testev    Adv. to make verb to test gxevxxx by matrices
NB.           of generator and shape given
NB.
NB. Version: 0.6.8 2011-07-14
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

EVBIGNUM=: % EVSMLNUM=: (%: FP_SFMIN) % FP_PREC

NB. ---------------------------------------------------------
NB. ggevli
NB. ggevui
NB.
NB. Description:
NB.   Calculate initial parameters for ggevxxx
NB.
NB. Syntax:
NB.   'abnrm abio ABupd plr hs'=. ggevxi AB
NB. where
NB.   AB    - 2×n×n-matrix, matrix pair (A,B)
NB.   abnrm - 2-vector, norms of A and B
NB.   abio  - 2-vector of integers, defines necessity and
NB.           value of scaling for A and B
NB.   ABupd - 2×n×n-matrix, scaled and balanced version on AB
NB.   plr   - 2×n-matrix of integers, permutations of A and
NB.           B, produced by ggbalxp
NB.   hs    - 2-vector of integers, defines submatrices
NB.           position, produced by ggbalxp
NB.
NB. TODO:
NB. - re-factorize to tacit adverb ggevi

ggevli=: 3 : 0
  abnrm=. (>./@,)"2 | y
  abio=. (0,(EVSMLNUM*1-FP_EPS),EVBIGNUM) I. abnrm
  abnrm ; abio ; ggballp (abnrm,.abio) scl^:((,{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 2 y
)

ggevui=: 3 : 0
  abnrm=. (>./@,)"2 | y
  abio=. (0,(EVSMLNUM*1-FP_EPS),EVBIGNUM) I. abnrm
  abnrm ; abio ; ggbalup (abnrm,.abio) scl^:((,{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 2 y
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
NB.   find eigenvalues e1, e2 and, optionally, left
NB.   eigenvectors L:
NB.     E2 * L * A = E1 * L * B                           (1)
NB.   and/or right eigenvectors R:
NB.     A * R^H * E2 = B * R^H * E1                       (2)
NB.   of pair of matrices (A,B). To avoid overflow,
NB.   eigenvalues of the matrix pair (A,B) are computed as a
NB.   pair of values. Each i-th eigenvector (row) from L and
NB.   R has a corresponding eigenvalue represented as a pair
NB.   of i-th elements from vectors e1 and e2:
NB.     E1=. diagmat(e1)
NB.     E2=. diagmat(e2)
NB.   If E2 is nonsingular then:
NB.     E=. diagmat(e1%e2)
NB.   is a diagonal matrix of eigenvalues, and GNEP (1), (2)
NB.   can be expressed as:
NB.     L * A = E * L * B                                 (3)
NB.     A * R^H = B * R^H * E                             (4)
NB.   and if E1 is nonsingular then:
NB.     E=. diagmat(e2%e1)
NB.   is a diagonal matrix of eigenvalues, and GNEP (1), (2)
NB.   can be expressed as:
NB.     E * L * A = L * B                                 (5)
NB.     A * R^H * E = B * R^H * E                         (6)
NB.
NB. Syntax:
NB.   e1e2=.      ggevlnn AB
NB.   'e1e2 R'=.  ggevlnv AB
NB.   'e1e2 L'=.  ggevlvn AB
NB.   'e1e2 LR'=. ggevlvv AB
NB. where
NB.   AB    - 2×n×n-matrix, matrix pair (A,B):
NB.             AB -: A ,: B
NB.   e1e2  - 2×n-matrix of eigenvalues e1 and e2:
NB.             e1e2 -: e1 ,: e2
NB.   L     - n×n-matrix of left eigenvectors
NB.   R     - n×n-matrix of right eigenvectors
NB.   LR    - 2×n×n-matrix of left and right eigenvectors:
NB.             LR -: L ,: R
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (ggevlnn -: ggevunn &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevlnv -: ggevunv &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevlvn -: ggevuvn &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevlvv -: ggevuvv &.: (ct"2)) A ,: B                   CHECKME!
NB.   (E2 mp L mp C) -: (E1 mp L mp D)
NB.   (C mp (ct R) mp E2) -: (D mp (ct R) mp E1)
NB. where
NB.   A - n×n-matrix, general
NB.   B - n×n-matrix, general
NB.   'e1e2 LR'=. ggevlvv A ,: B
NB.   'L R'=. LR

ggevlnn=: 3 : 0
  'abnrm abio y plr hs'=. ggevli y
  y=. (<0 1;;~dhs2lios hs) ([ (((0(0})hs)&gghrdlnn) upd) ((unmlqrc~,:trl@])gelqf)/@{`[`] }) y
  e1e2=. hs hgezqenn y
  e1e2=. (abnrm,.abio) scl^:((,~{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
)

ggevlnv=: 3 : 0
  'abnrm abio y plr hs'=. ggevli y
  y=. (0 1;(i.{.hs);<<dhs2lios hs) ((unmlqrc~,trl@],:unglq@])gelqf)/@({~<)~`((<0 1 2)<@(0})[)`((, idmat@c)@]) } y
  y=. (0,c y) gghrdlnv y
  y=. hs hgezqsnv y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2 ; _. $~ 2 # c y
  else.
    y=. tgevclrb y
    y=. gebaklp y ; {: plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr) y
    e1e2=. (abnrm,.abio) scl^:((,~{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevlvn=: 3 : 0
  'abnrm abio y plr hs'=. ggevli y
  y=. (<0 1;(i.{.hs);<<dhs2lios hs) ((unmlqrc~,:trl@])gelqf)/@{`[`] } y
  y=. (((0,]) gghrdlvn (,idmat)) c) y
  y=. hs hgezqsvn y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2 ; _. $~ 2 # c y
  else.
    y=. tgevculb y
    y=. gebaklp y ; {. plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr) y
    e1e2=. (abnrm,.abio) scl^:((,~{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevlvv=: 3 : 0
  'abnrm abio y plr hs'=. ggevli y
  y=. (0 1;(i.{.hs);<<dhs2lios hs) ((unmlqrc~,trl@],:unglq@])gelqf)/@({~<)~`((<0 1 3)<@(0})[)`((, ,:~@idmat@c)@]) } y
  y=. (0,c y) gghrdlvv y
  y=. hs hgezqsvv y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2 ; _. $~ 2 # c y
  else.
    y=. tgevcubb y
    y=. y gebaklp@;"2 1 plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr)"2 y
    e1e2=. (abnrm,.abio) scl^:((,~{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

NB. ---------------------------------------------------------
NB. ggevunn
NB. ggevunv
NB. ggevuvn
NB. ggevuvv
NB.
NB. Description:
NB.   Generalized nonsymmetric eigenvalue problem (GNEP):
NB.   find eigenvalues e1, e2 and, optionally, left
NB.   eigenvectors L:
NB.     E2 * L^H * A = E1 * L^H * B                       (7)
NB.   and/or right eigenvectors R:
NB.     A * R * E2 = B * R * E1                           (8)
NB.   of pair of matrices (A,B). To avoid overflow,
NB.   eigenvalues of the matrix pair (A,B) are computed as a
NB.   pair of values. Each i-th eigenvector (column) from L
NB.   and R has a corresponding eigenvalue represented as a
NB.   pair of i-th elements from vectors e1 and e2:
NB.     E1=. diagmat(e1)
NB.     E2=. diagmat(e2)
NB.   If E2 is nonsingular then:
NB.     E=. diagmat(e1%e2)
NB.   is a diagonal matrix of eigenvalues, and GNEP (7), (8)
NB.   can be expressed as:
NB.     L^H * A = E * L^H * B                             (9)
NB.     A * R = B * R * E                                (10)
NB.   and if E1 is nonsingular then:
NB.     E=. diagmat(e2%e1)
NB.   is a diagonal matrix of eigenvalues, and GNEP (7), (8)
NB.   can be expressed as:
NB.     E * L^H * A = L^H * B                            (11)
NB.     A * R * E = B * R * E                            (12)
NB.
NB. Syntax:
NB.   e1e2=.      ggevunn AB
NB.   'e1e2 R'=.  ggevunv AB
NB.   'e1e2 L'=.  ggevuvn AB
NB.   'e1e2 LR'=. ggevuvv AB
NB. where
NB.   AB    - 2×n×n-matrix, matrix pair (A,B):
NB.             AB -: A ,: B
NB.   e1e2  - 2×n-matrix of eigenvalues e1 and e2:
NB.             e1e2 -: e1 ,: e2
NB.   L     - n×n-matrix of left eigenvectors
NB.   R     - n×n-matrix of right eigenvectors
NB.   LR    - 2×n×n-matrix of left and right eigenvectors:
NB.             LR -: L ,: R
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (ggevunn -: ggevlnn &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevunv -: ggevlnv &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevuvn -: ggevlvn &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevuvv -: ggevlvv &.: (ct"2)) A ,: B                   CHECKME!
NB.   (E2 mp (ct L) mp C) -: (E1 mp (ct L) mp D)
NB.   (C mp R mp E2) -: (D mp R mp E1)
NB. where
NB.   A - n×n-matrix, general
NB.   B - n×n-matrix, general
NB.   'e1e2 LR'=. ggevuvv A ,: B
NB.   'L R'=. LR
NB.
NB. Notes:
NB. - ggevunn models LAPACK's xGGEV('N','N')
NB. - ggevunv models LAPACK's xGGEV('N','V')
NB. - ggevuvn models LAPACK's xGGEV('V','N')
NB. - ggevuvv models LAPACK's xGGEV('V','V')

ggevunn=: 3 : 0
  'abnrm abio y plr hs'=. ggevui y
  y=. (<0 1;;~dhs2lios hs) ([ (((0(0})hs)&gghrdunn) upd) ((unmqrlc~,:tru@])geqrf)/@{`[`] }) y
  e1e2=. hs hgeqzenn y
  e1e2=. (abnrm,.abio) scl^:((,~{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
)

ggevuvn=: 3 : 0
  'abnrm abio y plr hs'=. ggevui y
  y=. (0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,tru@],:ungqr@])geqrf)/@({~<)~`((<0 1 2)<@(0})[)`((, idmat@c)@]) } y
  y=. (0,c y) gghrduvn y
  y=. hs hgeqzsvn y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2 ; _. $~ 2 # c y
  else.
    y=. tgevculb y
    y=. gebakup y ; {. plr
    y=. (%"1 (EVSMLNUM&>`(,:&1))}@:normitc) y
    e1e2=. (abnrm,.abio) scl^:((,~{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevunv=: 3 : 0
  'abnrm abio y plr hs'=. ggevui y
  y=. (<0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,:tru@])geqrf)/@{`[`] } y
  y=. (((0,]) gghrdunv (,idmat)) c) y
  y=. hs hgeqzsnv y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2 ; _. $~ 2 # c y
  else.
    y=. tgevcurb y
    y=. gebakup y ; {: plr
    y=. (%"1 (EVSMLNUM&>`(,:&1))}@:normitc) y
    e1e2=. (abnrm,.abio) scl^:((,~{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevuvv=: 3 : 0
  'abnrm abio y plr hs'=. ggevui y
  y=. (0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,tru@],:ungqr@])geqrf)/@({~<)~`((<0 1 2)<@(0})[)`((, ,:~@idmat@c)@]) } y
  y=. (0,c y) gghrduvv y
  y=. hs hgeqzsvv y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2 ; _. $~ 2 # c y
  else.
    y=. tgevcubb y
    y=. y gebakup@;"2 1 plr
    y=. (%"1 (EVSMLNUM&>`(,:&1))}@:normitc)"2 y
    e1e2=. (abnrm,.abio) scl^:((,~{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

NB. =========================================================
NB. Test suite
