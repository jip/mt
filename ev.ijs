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

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. ggevlnn
NB. ggevlnv
NB. ggevlvn
NB. ggevlvv
NB.
NB. Description:
NB.   Find eigenvalues e1, e2 and, optionally, left
NB.   eigenvectors L:
NB.     E2 * L * A = E1 * L * B
NB.   and/or right eigenvectors R:
NB.     A * R^H * E2 = B * R^H * E1
NB.   of pair of matrices (A,B). Each i-th eigenvector
NB.   (column) from L and R has a corresponding eigenvalue
NB.   represented as a pair of i-th elements from e1 and e2:
NB.     E1 = diagmat(e1)
NB.     E2 = diagmat(e2)
NB.
NB. Syntax:
NB.   e1e2=.      ggevlnn AB
NB.   'e1e2 R'=.  ggevlnv AB
NB.   'e1e2 L'=.  ggevlvn AB
NB.   'e1e2 LR'=. ggevlvv AB
NB. where
NB.   AB    -: A ,: B
NB.   e1e2  -: e1 ,: e2
NB.   LR    -: L ,: R
NB.   L
NB.   R
NB.   A
NB.   B
NB.   e1
NB.   e2

ggevunn=: 3 : 0
  abnrm=. (>./@,)"2 | y
  abio=. (0,(EVSMLNUM*1-FP_EPS),EVBIGNUM) I. abnrm
  y=. (abnrm,.abio) scl^:((,{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 2 y
  'y plr hs'=. ggbalup y
  y=. (<0 1;;~dhs2lios hs) ([ (((0(0})hs)&gghrdunn) upd) ((unmqrlc~,:tru@])geqrf)/@{`[`] }) y
  e1e2=. hs hgeqzenn y
  e1e2=. (abnrm,.abio) scl^:((,~{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
)

NB. (E2 mp (ct L) mp C) -: (E1 mp (ct L) mp D)
ggevuvn=: 3 : 0
  abnrm=. (>./@,)"2 | y
  abio=. (0,(EVSMLNUM*1-FP_EPS),EVBIGNUM) I. abnrm
  y=. (abnrm,.abio) scl^:((,{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 2 y
  'y plr hs'=. ggbalup y
  y=. (0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,tru@],:ungqr@])geqrf)/@({~<)~`((<0 1 2)<@(0})[)`((, idmat@c)@]) } y
  y=. (0,c y) gghrduvn y  NB. y -: H , T ,: Q1
                          NB. Q1 = Q0 * dQ0
  y=. hs hgeqzsvn y       NB. y -: S , P ,: Q2
                          NB. Q2 = Q0 * dQ0 * dQ1
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

NB. (C mp R mp E2) -:&{. (D mp R mp E1)
ggevunv=: 3 : 0
  abnrm=. (>./@,)"2 | y
  abio=. (0,(EVSMLNUM*1-FP_EPS),EVBIGNUM) I. abnrm
  y=. (abnrm,.abio) scl^:((,{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 2 y
  'y plr hs'=. ggbalup y
  y=. (<0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,:tru@])geqrf)/@{`[`] } y
  y=. (((0,]) gghrdunv (,idmat)) c) y  NB. y -: H , T ,: dZ0
                                       NB. Z1 = Z0 * dZ0
                                       NB. Z0 = I
  y=. hs hgeqzsnv y                    NB. y -: S , P ,: Z2
                                       NB. Z2 = Z0 * dZ0 * dZ1
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
  abnrm=. (>./@,)"2 | y
  abio=. (0,(EVSMLNUM*1-FP_EPS),EVBIGNUM) I. abnrm
  y=. (abnrm,.abio) scl^:((,{&(1,EVSMLNUM,1,EVBIGNUM))/@[`({&0 1 0 1@{:@[)`])"1 2 y
  'y plr hs'=. ggbalup y
  y=. (0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,tru@],:ungqr@])geqrf)/@({~<)~`((<0 1 2)<@(0})[)`((, ,:~@idmat@c)@]) } y
  y=. (0,c y) gghrduvv y  NB. y -: H , T , Q1 ,: Z1
                          NB. Q1 = Q0 * dQ0
                          NB. Z1 = Z0 * dZ0
                          NB. Z0 = I
  y=. hs hgeqzsvv y       NB. y -: S , P , Q2 ,: Z2
                          NB. Q2 = Q0 * dQ0 * dQ1
                          NB. Z2 = Z0 * dZ0 * dZ1
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
