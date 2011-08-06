NB. Eigenvalues and, optionally, eigenvectors
NB.
NB. ggevxxx   Eigenvalues and, optionally, eigenvectors of
NB.           pair of matrices
NB.
NB. testggev  Test ggevxxx by general matrices given
NB. testev    Adv. to make verb to test gxevxxx by matrices
NB.           of generator and shape given
NB.
NB. Version: 0.6.8 2010-10-30
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

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:          Side:        Syntax:
NB. ggevlnn        lower        ab=.       ggevlnn AB
NB. ggevlnv        lower        'ab VR'=.  ggevlnv AB
NB. ggevlvn        lower        'ab VL'=.  ggevlvn AB
NB. ggevlvv        lower        'ab VLR'=. ggevlnv AB
NB. ggevunn        upper        ab=.       ggevunn AB
NB. ggevunv        upper        'ab VR'=.  ggevunv AB
NB. ggevuvn        upper        'ab VL'=.  ggevuvn AB
NB. ggevuvv        upper        'ab VLR'=. ggevunv AB
NB.
NB. Description:
NB.   Find eigenvalues and, optionally, eigenvectors of
NB.   pair of matrices
NB. where
NB.   AB    -: A ,: B
NB.   ab    -: alpha ,: beta
NB.   VRL   -: VL ,: VR
NB.   VL
NB.   VR
NB.   A
NB.   B
NB.   alpha
NB.   beta

EV_BIGNUM=. % EV_SMLNUM=. (%: FP_SFMIN) % FP_PREC

ggevunn=: 3 : 0
  'anrm bnrm'=. abnrm=. (>./ @ ,)"2 | AB
  abio=. (0 , (EV_SMLNUM - ?????????) , EV_BIGNUM) I. abnrm  NB. FIXME!
  abnrmto=. abio { 1 , EV_SMLNUM , 1 , EV_BIGNUM
  ilabscl=. abio { 0 1 0 1
  AB=. (abnrm ,. abnrmto) (scl ^: ilabscl)"1 2 AB            NB. FIXME!
  'AB plr hs'=. bbgalup AB
  QfR=. (1 1 ,. ,.~ hs) geqrf@{: ;. 0 AB
  QhA=. QfR unmqrlc (0 1 ,. ,.~ ns) {. ;. 0 AB
  iosHT=. < 0 2 ; ;~ dhs2lios hs
  AB=. iosHT (gghrdu @ ((0,s)&;) upd AB
  ab=. hgeqzue hs ; AB
  ab=. (abnrmto ,. abnrm) (scl ^: ilabscl)"1 1 ab            NB. FIXME! TODO: ^:_1
)

NB. =========================================================
NB. Test suite
