NB. Rotation
NB.
NB. lartg      Generates a plane rotation of a 2-vector
NB. rot        Applies a plane rotation(s) to a 2-vector(s)
NB. rotga      Adv. to make monad to generate and apply
NB.            rotation
NB. rotsclx    Update array by rotations and scalings
NB.            accumulated
NB.
NB. testlartg  Test lartg by vectors
NB. testlartv  Test rot by vectors
NB. testrot    Adv. to make verb to test rotation algorithms by
NB.            vectors of generator given
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
NB. References:
NB. [1] David S. Bindel, James W. Demmel, W. Kahan, Osni A.
NB.     Marques. On Computing Givens rotations reliably and
NB.     efficiently. UT-CS-00-449, October 2000. LAPACK
NB.     Working Note 148
NB.     http://www.netlib.org/lapack/lawns/downloads/
NB. [2] Edward Anderson. Discontinuous Plane Rotations and
NB.     the Symmetric Eigenvalue Problem. University of
NB.     Tennessee, UT-CS-00-454, December 4, 2000.
NB.     LAPACK Working Note 150
NB.     http://www.netlib.org/lapack/lawns/downloads/
NB. [3] Edward Anderson. 2017. Algorithm 978: Safe scaling in
NB.     the level 1 BLAS. ACM Trans. Math. Softw. 44, 1
NB.     (July 2017).
NB.     https://doi.org/10.1145/3061665
NB. [4] Weslley da Silva Pereira, Ali Lotfi, Julien Langou.
NB.     2022. Numerical analysis of Givens rotation. arXiv
NB.     https://doi.org/10.48550/arXiv.2211.04010

NB. =========================================================
NB. Configuration

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Signum, monad, sgn(0)=1

NB. y is real or complex
sgnc=: *!.0`1:@.(0&(=!.0))

NB. y is quaternion
NB. throws NaN error for y containing either quaternion
NB. infinity or directed infinity with both atoms non-zero
sgnq=: qnsign`(1 0"_)@.(0 0&(-:!.0))@(dbsig@33^:(-.@(0&(e.!.0)) *. (_&e. > (, _) -: -.!.0&0)@:|@(,@:+.^:(JCMPX = 3!:0))))

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. lartg
NB.
NB. Description:
NB.   Generate left multiplier cs to make plane rotation
NB.   from fg to r0:
NB.     r0 := cs * fg
NB.   such that cs's 1st Cayley-Dickson half is real; r0's
NB.   2nd Cayley-Dickson half is 0. fg, cs and r0 are of the
NB.   same type: either complex or quaternion.
NB.
NB. Formula:
NB.   cs := sgn(cd1st(fg)) * sgn(qnconik(fg))
NB. where
NB.   cd1st()   - extract 1st half from Cayley-Dickson pair:
NB.                 cd1st(a + b*i) = a
NB.                 cd1st(a + b*i + c*j + d*k) = a + b*i
NB.   qnconik() - conjugate i and k components:
NB.                 qnconik(a + b*i + c*j + d*k) = a - b*i + c*j - d*k
NB.   sgn(x)    = x/|x|, if x ≠ 0
NB.             = 1,     if x = 0
NB.
NB. Syntax:
NB.   cs=. lartg fg
NB. where
NB.   fg - 2-vector (f,g) to rotate
NB.   cs - 2-vector (c,s), the rotator, c is real
NB.   r0 - 2-vector (r,0), the rotated fg
NB.          r0 := cs * fg
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (r,0) -: (c,s) rot (f,g)
NB.   (r,0) -: (f,g) (qnmul qnconj) (c,s)  NB. equiv. to (qnmul qnconij) since c is real
NB.   (0,r) -: (c,s) (rot~ qnconj)~ (g,f)  NB. equiv. to (rot~ qnconij)~ since c is real
NB.   (0,r) -: (g,f) qnmul (c,s)
NB.   (-: +) c
NB. where
NB.   'c s'=. lartg (f,g)
NB.   r=. (c,s) mp (f,g)
NB.
NB. Notes:
NB. - simulates LAPACK's xLARTG
NB. - [G]SEP requires plane rotation to be continuous [2]:
NB.     cs := sgn(cd1st(cd1st(fg))) * sgn(cd1st(fg)) * sgn(qnconik(fg))
NB.   To achieve this, use modified definition:
NB.     lartg=: ((*&sgnc 9&o.)@{. * sgnq@:+)"1

lartg=: 9&o.&.(0&{)@(sgnc@{. * sgnq@:+)`(2 # nan)@.(isnan@<)"1 : [:

NB. ---------------------------------------------------------
NB. rot
NB.
NB. Description:
NB.   Applies plane rotation(s) cs to pair(s) ixy, anyone of:
NB.     oxy    := cs    * ixy
NB.   or
NB.     oxy[i] := cs[i] * ixy
NB.   or
NB.     oxy[i] := cs    * ixy[i]
NB.   or
NB.     oxy[i] := cs[i] * ixy[i]
NB.   for i=0:n-1
NB.
NB. Syntax:
NB.     oxy=. cs rot ixy
NB. where
NB.   ixy - 2-vector of Cayley-Dickson halves (ix,iy) or
NB.         n×2-matrix of laminated Cayley-Dickson halves
NB.         (ix[i],iy[i]), defines pair(s) to rotate
NB.         (complex or quaternion)
NB.   cs  - 2-vector of Cayley-Dickson halves (c,s) or
NB.         n×2-matrix of laminated Cayley-Dickson halves
NB.         (c[i],s[i]), defines rotator(s) (complex or
NB.         quaternion)
NB.   oxy - has the shape of either ixy or cs which has
NB.         greater rank, the rotated pair(s)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   cs (rot -: (qnmul qnconj)~) ixy
NB.
NB. Application:
NB. - with 2-rank ixy:
NB.     sentence    rank(cs)    implements BLAS's/LAPACK's
NB.     --------    --------    --------------------------
NB.     rot         1           xROT   with INCX=INCY     =1  (goes along columns)
NB.     rot&.|:     1           xROT   with INCX=INCY     =LD (goes along rows   )
NB.     rot         2           xLARTV with INCX=INCY=INCC=1  (goes along columns)
NB.     rot&.|:     2           xLARTV with INCX=INCY=INCC=LD (goes along rows   )
NB.
NB. Notes:
NB. - implements BLAS's DROT and LAPACK's ZROT and xLARTV
NB. - resembles qnmul

rot=: [: : (,.~@:(+./"1)@isnan`(,: nan)}@(mp"2 1~ (,:"1 (+@:-@:({:"1) ,. {."1)))~)

NB. ---------------------------------------------------------
NB. rotga
NB.
NB. Description:
NB.   Adv. to make monad to generate and apply rotation
NB.
NB. Syntax:
NB.   'Aupd cs'=. (vrota rotga) A ; isosubA ; isofg
NB. where
NB.   vrota   - dyad to apply rotation; is called as:
NB.               subAupd=. cs vrota subA
NB.             and is any of:
NB.               rot      NB. apply rotation to rows
NB.               rot&.|:  NB. apply rotation to columns
NB.   cs      - 2-vector (c,s), output of lartg, defines
NB.             rotation matrix
NB.   A       - m×n-matrix to update
NB.   Aupd    - A with subA replaced by subAupd
NB.   subA    - 2×any-matrix or any×2-matrix, array of
NB.             2-vectors to apply rotation
NB.   subAupd - subA rotated
NB.   isosubA - ISO subA (subAupd) within A (Aupd)
NB.   isofg   - ISO within subA of 2-vector (f,g) which
NB.             defines rotation

rotga=: 1 : 0
  'A isosubA isofg'=. y
  subA=. isosubA { A
  cs=. lartg_mt_ isofg { subA
  ((cs u subA) isosubA} A) ; cs
)

NB. ---------------------------------------------------------
NB. rotscll
NB. rotsclu
NB.
NB. Description:
NB.   Update A by rotations and scalings accumulated in dA
NB.
NB. Syntax:
NB.   Aupd=. A rotsclx dA
NB. where
NB.   dA      - any×4-matrix, where each row is 4-vector of
NB.             values, either:
NB.               0 0 0 0         NB. no action
NB.             or:
NB.               m , io , 0 0    NB. defines scaling
NB.             or:
NB.               cs , iof , iog  NB. defines rotation
NB.             accumulates scalings and rotations
NB.   A       - n×n-matrix or (i.0)
NB.   Aupd    - either (i.0) when A -: (i.0) , or n×n-matrix
NB.             (A*dA) otherwise
NB.   m       - multiplier to scale either row (rotscll) or
NB.             column (rotsclu)
NB.   io      - IO either row (rotscll) or column (rotsclu)
NB.             to scale
NB.   cs      - 2-vector (c,s), output of lartg, defines
NB.             rotation matrix
NB.   iof,iog - ISO either rows (rotscll) or columns
NB.             (rotsclu) which contain 2-vectors (f,g) to
NB.             rotate, iof≠iog
NB.
NB. TODO:
NB. - aggregate non-intersecting groups of vectors to change
NB.   them simultaneously

rotscll=: 4 : 0
  i=. 0
  while. i < # y do.                NB. traverse dA rows down
    'cs iofg'=. _2 ]\ i { y
    if. 0 0 -: iofg do.
      if. -. 0 0 -: cs do.
        'm io'=. cs
        x=. m&*&.(io&{) x           NB. do scale
      end.
    else.
      x=. cs&(rot&.|:)&.(iofg&{) x  NB. do rotation
    end.
    i=. >: i
  end.
  x
)

rotsclu=: 4 : 0
  i=. 0
  while. i < # y do.                   NB. traverse dA rows down
    'cs iofg'=. _2 ]\ i { y
    if. 0 0 -: iofg do.
      if. -. 0 0 -: cs do.
        'm io'=. cs
        x=. m&*&.((< a: ; io)&{) x     NB. do scale
      end.
    else.
      x=. cs&rot&.((< a: ; iofg)&{) x  NB. do rotation
    end.
    i=. >: i
  end.
  x
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testlartg
NB.
NB. Description:
NB.   Test:
NB.   - xLARTG (math/lapack2 addon)
NB.   - lartg (math/mt addon)
NB.   by vectors
NB.
NB. Syntax:
NB.   log=. testlartg FG
NB. where
NB.   FG  - m×2-matrix, m (f,g) pairs to test
NB.   log - 6-vector of boxes, test log
NB.
NB. Formula:
NB.   relative error in the singular value:
NB.     (sqrt(c^2 + |s|^2) - 1) % ε                       (1)
NB.   relative backward error:
NB.     ||Q^H(r,0)^T-(f,g)^T||_2 / (||(f,g)||_2 * ε)      (2)
NB.   where
NB.     Q = ( c       s )
NB.         (-conj(s) c )
NB.     c ∈ ℝ
NB.     Q * ( f g )^T = ( r 0 )^T
NB.   relative error of multiple rotations:
NB.     ((Π σ_i, i=1:M) - 1) % ε                          (3)
NB.   where
NB.     σ_i = sqrt(c_i^2 + |s_i|^2)
NB.
NB. Algorithm to calculate errors:
NB.   In:  FG, CS
NB.   where
NB.     CS=. lartg FG
NB.     CS -: C ,. S
NB.   Out: fwderr, bwderr
NB.   1) find R:
NB.        R=. CS mp"1 FG
NB.   2) exclude rows containing NaN or ∞ in (CS ,. R)
NB.   3) calculate (1) for each (c,s) pair, exclude NaN and
NB.      ∞, find maximum, put into the "fwd error" column
NB.   4) calculate (2) for each (f,g,c,s,r) quintet, exclude
NB.      NaN and ∞, find maximum, put into the "bwd error"
NB.      column
NB.
NB. Application:
NB. - test by testset from [1]
NB.     zpow=. <. 1r2 + 1r4 * 2 ^. FP_EPS % FP_SFMIN  NB. S,C: 26; D,Z: 242; z == 2 ^ zpow
NB.     ptsa=. 2 ^ zpow * , (i: 4) +/ (16 %~ i: 1)    NB. 27 test points absolute values which power of 2 are: (1) below by 1r16, (2) equal, (3) above by 1r16
NB.     pts=. (0 , (, -)) ptsa                        NB. 55 test points
NB.     FG=.           (2 permrep # pts) {"1 pts      NB. for S,D only: 55^2 real    test pairs r with values                              from pts vector
NB.     FG=. _2 j./\"1 (4 permrep # pts) {"1 pts      NB. for C,Z only: 55^4 complex test pairs   with numbers c of components Re(c),Im(c) from pts vector
NB.     'C S'=. |: CS=. lartg FG                      NB. (f,g) -> (c,s)
NB.     R=. CS mp"1 FG                                NB. Q*(f,g)^T = (r,0)^T
NB. - test by testset from [3]
NB.     ptsa=. (_ , (, %)) FP_SFMIN , FP_PREC , 3     NB. we use 3 instead of 1 here as Medium value
NB.     pts=. (0 _. , (, -)) ptsa
NB.     FG=.           (2 permrep # pts) {"1 pts      NB. for S,D only: 16^2 real    test pairs r with values                              from pts vector
NB.     FG=. _2 j./\"1 (4 permrep # pts) {"1 pts      NB. for C,Z only: 16^4 complex test pairs   with numbers c of components Re(c),Im(c) from pts vector
NB.     'C S'=. |: CS=. lartg FG                      NB. (f,g) -> (c,s)
NB.     R=. CS mp"1 FG                                NB. Q*(f,g)^T = (r,0)^T
NB.     assert (F +.&isnan   G) *.  isnan           R
NB.     assert (F +.&(_ = |) G) *. (isnan +. _ = |) R
NB.     assert (C (,. ,:"1 (,.~ +@:-)) S) mp"2 1 FG) =!.0"1 R ,. 0
NB.     assert 1 =!.0 (9 o. C) +&.*: | S
NB.     assert 0 <: 9 o. C
NB. - test by testset from [4]
NB.     rhoMax=. IF64 { 51 484.5
NB.     mklen=. 2 ^ (_1 1 * rhoMax)&randu
NB.     'M N'=. MN=. 1e5 1e3                          NB. as in [4] for multiple rotations
NB.     FG=.  mklen                 MN , 2            NB. for S,D only: M×N real    random pairs (f,g)
NB.     FG=. (mklen r. 0 2p1&randu) MN , 2            NB. for C,Z only: M×N complex random pairs (f,g)
NB.     CS=. lartg FG                                 NB. M×N pairs (c,s)
NB.     R=. CS mp"1 FG                                NB. M×N scalars r, Q*(f,g)^T = (r,0)^T
NB.     NB. 1. accuracy of a single rotation
NB.     N1=. 1e6
NB.     NB. 1a. only N1 pairs are used in [4] for a single rotation test
NB.     FG1=. N1 ({. ,/) FG                           NB. select N1 pairs
NB.     CS1=. N1 ({. ,/) CS                           NB. select N1 pairs
NB.     R1=.  N1 ({. , ) R                            NB. select N1 scalars
NB.     NB. 1b. relative error in the singular value (1):
NB.     resv=. FP_EPS %~ <: ((+&.*: 9&o.)~ |)/ |: CS1
NB.     plot resv
NB.     echo (mean , stddev      )   resv             NB. avg( err ),std( err )
NB.     echo (mean , stddev , max) | resv             NB. avg(|err|),std(|err|),max(|err|)
NB.     NB. 1c. relative backward error (2):
NB.     rbe=. 2 (}."1 ((- (FP_EPS %~ %)&:norms"1 ]) (({."1 ((,. -) ,:"1 (,.~ +)) {:"1)@:(}:"1) mp"2 1 {:"1 ,. 0:))~ {."1) FG1 ,. CS1 ,. R1
NB.     plot rbe
NB.     echo (mean , stddev      )   rbe              NB. avg( err ),std( err )
NB.     echo (mean , stddev , max) | rbe              NB. avg(|err|),std(|err|),max(|err|)
NB.     NB. 2. accuracy of multiple rotations on 2×2-matrices
NB.     svs=. ((+&.*: 9&o.)~ |)/@|:"2 CS              NB. singular values, M×N-matrix
NB.     mux=. mean svs                                NB. E[svs], N-vector
NB.     sigmax=. stddev svs                           NB. σ_svs, N-vector
NB.     NB. 2a. estimate (singular value)^M
NB.     muy=. mux
NB.     sigmay=. muy * %: <: ^ M * *: sigmax % mux
NB.     NB. 2b. relative error of multiple rotations (3):
NB.     rep=. FP_EPS %~ <: */ svs
NB.     repx=. steps (min,max,50"_) rey
NB.     plot repx ([ ; }:@histogram) rey              NB. https://code.jsoftware.com/wiki/Essays/Histogram

testlartg=: 3 : 0
  _1 cocreate < 'mttmp'
  load_mttmp_ 'math/mt/external/lapack2/lartg'

  NB. exclude items containing NaN or ∞ from vector or matrix
  xicni=: #~ -.@(+./)"1@(__&= ,. _&= ,. isnan)

  NB. compute relative error in the singular value (1)
  NB. err1=. vresv C ,. S
  vresv=: FP_EPS >./@xicni@:%~ <:@(((+&.*: 9&o.)~ |)/)@|:

  NB. compute relative backward error (2)
  NB.   err2=. 2 vrbe (F ,. G ,. C ,. S ,. R)
  vrbe=: >./@xicni@(}."1 ((- (FP_EPS %~ %)&:norms"1 ]) (({."1 ((,. -) ,:"1 (,.~ +)) {:"1)@:(}:"1) mp"2 1 {:"1 ,. 0:))~ {."1)

  log=.          ('dlartg_mttmp_"1' tmonad (]`]`nan`(vresv@:>@(2 {."1 ]))`(2 vrbe xicni@(,. >       )))) y
  log=. log lcat ('zlartg_mttmp_"1' tmonad (]`]`nan`(vresv@:>@(2 {."1 ]))`(2 vrbe xicni@(,. >       )))) y
  log=. log lcat ('lartg'           tmonad (]`]`nan`(vresv@           ] )`(2 vrbe xicni@(,. ,. mp"1~)))) y

  coerase < 'mttmp'
  erase 'xicni vresv vrbe'

  log
)

NB. ---------------------------------------------------------
NB. testlartv
NB.
NB. Description:
NB.   Test:
NB.   - DROT (test BLAS in math/mt addon)
NB.   - ZROT (math/lapack2 addon)
NB.   - xLARTV (math/lapack2 addon)
NB.   - rot (math/mt addon)
NB.   by scalars and vectors
NB.
NB. Syntax:
NB.   log=. testlartv (X ,. Y ,. F ,. G)
NB. where
NB.   (X ,. Y) - m×2-matrix, m (x,y) pairs to be rotated
NB.   (F ,. G) - m×2-matrix, m (f,g) pairs to make rotators
NB.   log      - 6-vector of boxes, test log

testlartv=: 3 : 0
  _1 cocreate < 'mttmp'
  load        'math/mt/external/blas/drot'
  load_mttmp_ 'math/mt/external/lapack2/lartv'
  load_mttmp_ 'math/mt/external/lapack2/zrot'

  xycs=. ;/ |: y=. 2 ({."1 ,. lartg@:(}."1)) y

  log=.          ('drot_mtbla_"1' tmonad (]`]`nan`nan`nan)) y
  log=. log lcat ('zrot_mttmp_"1' tmonad (]`]`nan`nan`nan)) y

  log=. log lcat ('dlartv_mttmp_' tmonad (]`]`nan`nan`nan)) xycs
  log=. log lcat ('zlartv_mttmp_' tmonad (]`]`nan`nan`nan)) xycs

  log=. log lcat ('rot' tdyad (2&(}."1)`(2&({."1))`]`nan`nan`nan)) y
  log=. log lcat ('rot' tdyad (2&(}."1)`(2&({."1))`]`nan`nan`nan)) y
  log=. log lcat ('rot' tdyad (2&(}."1)`(2&({."1))`]`nan`nan`nan)) y
  log=. log lcat ('rot' tdyad (2&(}."1)`(2&({."1))`]`nan`nan`nan)) y

  coerase < 'mttmp'

  log
)

NB. ---------------------------------------------------------
NB. testrot
NB.
NB. Description:
NB.   Adv. to make verb to test rotation algorithms by matrix
NB.   of generator and shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testrot) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testrot_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testrot_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testrot_mt_ 150 200

testrot=: 1 : 'testlartv_mt_@u@(4&(1})) lcat_mt_~ testlartg_mt_@u@(2&(1}))'
