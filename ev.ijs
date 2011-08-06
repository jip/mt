NB. Eigenvalues and eigenvectors
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

NB. Scaling limits
EVBIGNUM=: % EVSMLNUM=: (%: FP_SFMIN) % FP_PREC

NB. Scaling factors
EVSCL=: 1,EVSMLNUM,1,EVBIGNUM

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

ggevi=: 1 : '(,.(0,(EVSMLNUM_mt_*1-FP_EPS_mt_),EVBIGNUM_mt_)&I.)@:(>./@,"2)@:| ([ ; u@:(scl_mt_^:((,{&EVSCL_mt_)/@[`({&0 1 0 1@{:@[)`])"1 2)) ]'

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
NB.   L     - n×n-matrix. left eigenvectors (rows)
NB.   R     - n×n-matrix. right eigenvectors (rows)
NB.   LR    - 2×n×n-matrix. left and right eigenvectors:
NB.             LR -: L ,: R
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (ggevlnn -: ggevunn &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevlnv -: ggevunv &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevlvn -: ggevuvn &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevlvv -: ggevuvv &.: (ct"2)) A ,: B                   CHECKME!
NB.   (E2 mp L mp A) -: (E1 mp L mp B)
NB.   (A mp (ct R) mp E2) -: (B mp (ct R) mp E1)
NB. where
NB.   A - n×n-matrix, general
NB.   B - n×n-matrix, general
NB.   'e1e2 LR'=. ggevlvv A ,: B
NB.   'E1 E2'=. diagmat"1 e1e2
NB.   'L R'=. LR

ggevlnn=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi y
  y=. (<0 1;;~dhs2lios hs) ([ ((gghrdlnn~0,c) upd) ((unmlqrc~,:trl@:(}:"1)@])gelqf)/@{`[`] }) y
  e1e2=. hs hgezqenn y
  e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
)

ggevlnv=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi y
  y=. (0 1;(<i.{.hs);dhs2lios hs) ((unmlqrc~,(trl@:(}:"1),:unglq)@])gelqf)/@({~<)~`((<0 1 2)<@(0})[)`((, idmat@c)@]) } y
  y=. (gghrdlnv~0,c) y
  y=. hs hgezqsnv y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2 ; _. $~ 2 # c y
  else.
    y=. tgevclrb y
    y=. gebaklp y ; {: plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr) y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevlvn=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi y
  y=. (<0 1;(<i.{.hs);dhs2lios hs) ((unmlqrc~,:trl@:(}:"1)@])gelqf)/@{`[`] } y
  y=. (((0,]) gghrdlvn (,idmat)) c) y
  y=. hs hgezqsvn y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2 ; _. $~ 2 # c y
  else.
    y=. tgevcllb y
    y=. gebaklp y ; {. plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr) y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevlvv=: 3 : 0
  'abnrmio y plr hs'=. ggballp ggevi y
  y=. (0 1;(<i.{.hs);dhs2lios hs) ((unmlqrc~,(trl@:(}:"1),:unglq)@])gelqf)/@({~<)~`((<0 1 3)<@(0})[)`((, ,:~@idmat@c)@]) } y
  y=. (gghrdlvv~0,c) y
  y=. hs hgezqsvv y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2 ; _. $~ 2 # c y
  else.
    y=. tgevclbb y
    y=. y gebaklp@;"2 1 plr
    y=. (% (EVSMLNUM&>`(,:&1))}@:normitr)"2 y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
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
NB.   L     - n×n-matrix. left eigenvectors (columns)
NB.   R     - n×n-matrix. right eigenvectors (columns)
NB.   LR    - 2×n×n-matrix. left and right eigenvectors:
NB.             LR -: L ,: R
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (ggevunn -: ggevlnn &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevunv -: ggevlnv &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevuvn -: ggevlvn &.: (ct"2)) A ,: B                   CHECKME!
NB.   (ggevuvv -: ggevlvv &.: (ct"2)) A ,: B                   CHECKME!
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
NB. - model LAPACK's xGEEV('N','N'):
NB.     geevlnn=: ggevlnn @ (, (idmat @ c))
NB.     geevunn=: ggevunn @ (, (idmat @ c))
NB. - model LAPACK's xGEEV('N','V'):
NB.     geevlnv=: ggevlnv @ (, (idmat @ c))
NB.     geevunv=: ggevunv @ (, (idmat @ c))
NB. - model LAPACK's xGEEV('V','N'):
NB.     geevlvn=: ggevlvn @ (, (idmat @ c))
NB.     geevuvn=: ggevuvn @ (, (idmat @ c))
NB. - model LAPACK's xGEEV('V','V'):
NB.     geevlvv=: ggevlvv @ (, (idmat @ c))
NB.     geevuvv=: ggevuvv @ (, (idmat @ c))
NB.
NB. Notes:
NB. - ggevunn models LAPACK's xGGEV('N','N')
NB. - ggevunv models LAPACK's xGGEV('N','V')
NB. - ggevuvn models LAPACK's xGGEV('V','N')
NB. - ggevuvv models LAPACK's xGGEV('V','V')

ggevunn=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi y
  y=. (<0 1;;~dhs2lios hs) ([ ((gghrdunn~0,c) upd) ((unmqrlc~,:tru@}:@])geqrf)/@{`[`] }) y
  e1e2=. hs hgeqzenn y
  e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
)

ggevuvn=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi y
  y=. (0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,(tru@}:,:ungqr)@])geqrf)/@({~<)~`((<0 1 2)<@(0})[)`((, idmat@c)@]) } y
  y=. (gghrduvn~0,c) y
  y=. hs hgeqzsvn y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2 ; _. $~ 2 # c y
  else.
    y=. tgevculb y
    y=. gebakup y ; {. plr
    y=. (%"1 (EVSMLNUM&>`(,:&1))}@:normitc) y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevunv=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi y
  y=. (<0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,:tru@}:@])geqrf)/@{`[`] } y
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
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

ggevuvv=: 3 : 0
  'abnrmio y plr hs'=. ggbalup ggevi y
  y=. (0 1;(dhs2lios hs);<<i.{.hs) ((unmqrlc~,(tru@}:,:ungqr)@])geqrf)/@({~<)~`((<0 1 2)<@(0})[)`((, ,:~@idmat@c)@]) } y
  y=. (gghrduvv~0,c) y
  y=. hs hgeqzsvv y
  e1e2=. 2 {. diag"2 y
  if. 128!:5 < e1e2 do.
    NB. non-converged
    e1e2 ; _. $~ 2 # c y
  else.
    y=. tgevcubb y
    y=. y gebakup@;"2 1 plr
    y=. (%"1 (EVSMLNUM&>`(,:&1))}@:normitc)"2 y
    e1e2=. abnrmio scl^:((,~{&EVSCL)/@[`({&0 1 0 1@{:@[)`])"1 1 e1e2
    e1e2 ; y
  end.
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testggev
NB.
NB. Description:
NB.   Test ggevxxx by general matrices given
NB.
NB. Syntax:
NB.   testggev AB
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
NB.   - ggevlvn:
NB.       berr0 := max(||l(i) * (β(i)*A - α(i)*B)  || / (FP_PREC * max(|| β(i)*A   ||,|| α(i)*B   ||)))
NB.       berr1 := max(| ||l(i)|| - 1 |) / (FP_PREC * n)
NB.   - ggevlnv:
NB.       berr0 := max(||r(i) * (β(i)*A - α(i)*B)^H|| / (FP_PREC * max(||(β(i)*A)^H||,||(α(i)*B)^H||)))
NB.       berr1 := max(| ||r(i)|| - 1 |) / (FP_PREC * n)
NB.   - ggevlvv:
NB.       berr0 := berr(ggevlvn)
NB.       berr1 := berr(ggevlnv)
NB.   - ggevuvn:
NB.       berr0 := max(||(β(i)*A - α(i)*B)^H * l(i)|| / (FP_PREC * max(||(β(i)*A)^H||,||(α(i)*B)^H||)))
NB.       berr1 := max(| ||l(i)|| - 1 |) / (FP_PREC * n)
NB.   - ggevunv:
NB.       berr0 := max(||(β(i)*A - α(i)*B)   * r(i)|| / (FP_PREC * max(|| β(i)*A   ||,|| α(i)*B   ||)))
NB.       berr1 := max(| ||r(i)|| - 1 |) / (FP_PREC * n)
NB.   - ggevuvv:
NB.       berr0 := berr(ggevuvn)
NB.       berr1 := berr(ggevunv)
NB.
NB. Notes:
NB. - berrxxx are non-iterative and are require O(N^3) RAM

testggev=: 3 : 0
  vberrlvn=:  (     normir@:((*"_ 1|:@|.@(0&{::))((  norm1r@:((mp"1 2 (-/"3))~     )                                           )      % (FP_PREC*(FP_SFMIN>.    (>./"1)@:( norm1       "2)@[)))(1 {:: ])))>.((     normir@:<:@:normitr%FP_PREC*c)@(1 {:: ]))
  vberrlnv=:  (     normir@:((*"_ 1|:@|.@(0&{::))((                                   norm1r@:((mp"2 1~(-/"3))~+    )          )      % (FP_PREC*(FP_SFMIN>.    (>./"1)@:(       normi "2)@[)))(1 {:: ])))>.((     normir@:<:@:normitr%FP_PREC*c)@(1 {:: ]))
  vberrlvv=:  (>./@:normir@:((*"_ 1|:@|.@(0&{::))((((norm1r@:( mp"1 2        ~   {.),:norm1r@:((mp"2 1       ) + @{:))~(-/"3))~)(>./@:%)(FP_PREC*(FP_SFMIN>.|:@:(>./"1)@:((norm1,normi)"2)@[)))(1 {:: ])))>.((>./@:normir@:<:@:normitr%FP_PREC*c)@(1 {:: ]))

  vberruvn=:  (     normir@:((*"_ 1|:@|.@(0&{::))((  norm1r@:((mp"1 2 (-/"3))~ct   )                                           )      % (FP_PREC*(FP_SFMIN>.    (>./"1)@:( normi       "2)@[)))(1 {:: ])))>.((     normir@:<:@:normitc%FP_PREC*c)@(1 {:: ]))
  vberrunv=:  (     normir@:((*"_ 1|:@|.@(0&{::))((                                   norm1r@:((mp"2 1~(-/"3))~|:   )          )      % (FP_PREC*(FP_SFMIN>.    (>./"1)@:(       norm1 "2)@[)))(1 {:: ])))>.((     normir@:<:@:normitc%FP_PREC*c)@(1 {:: ]))
  vberruvv=:  (>./@:normir@:((*"_ 1|:@|.@(0&{::))((((norm1r@:( mp"1 2        ~ct@{.),:norm1r@:((mp"2 1       ) |:@{:))~(-/"3))~)(>./@:%)(FP_PREC*(FP_SFMIN>.|:@:(>./"1)@:((normi,norm1)"2)@[)))(1 {:: ])))>.((>./@:normir@:<:@:normitc%FP_PREC*c)@(1 {:: ]))
gy=: y
  rcond=. <./ gecon1"2 y

  ('ggevlnn' tmonad (]`]`(rcond"_)`(_."_)`(_."_)  )) y
  ('ggevlnv' tmonad (]`]`(rcond"_)`(_."_)`vberrlnv)) y
  ('ggevlvn' tmonad (]`]`(rcond"_)`(_."_)`vberrlvn)) y
  ('ggevlvv' tmonad (]`]`(rcond"_)`(_."_)`vberrlvv)) y
  ('ggevunn' tmonad (]`]`(rcond"_)`(_."_)`(_."_)  )) y
  ('ggevunv' tmonad (]`]`(rcond"_)`(_."_)`vberrunv)) y
  ('ggevuvn' tmonad (]`]`(rcond"_)`(_."_)`vberruvn)) y
  ('ggevuvv' tmonad (]`]`(rcond"_)`(_."_)`vberruvv)) y

  erase 'vberrlnv vberrlvn vberrlvv vberrunv vberruvn vberruvv'

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
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testev_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testev_mt_ 150 150

testev=: 1 : 'EMPTY_mt_ [ (testggev_mt_ @ u @ (2&,)) ^: (=/)'