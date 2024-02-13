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
NB. testrot    Test rotation algorithms by predefined matrix
NB.
NB. verifyrot  Verify rot verbs
NB.
NB. Version: 0.11.0 2021-01-17
NB.
NB. Copyright 2010-2021 Igor Zhuravlov
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

lartg=: 9&o.&.(0&{)@(sgnc@{. * sgnq@:+)`(_.j_. _.j_."_)@.(isnan@<)"1 : [:

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

rot=: [: : (,.~@:(+./"1)@isnan`(,:&_.j_.)}@(mp"2 1~ (,:"1 (+@:-@:({:"1) ,. {."1)))~)

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
NB.   Test lartg by vectors
NB.
NB. Syntax:
NB.   log=. testlartg FG
NB. where
NB.   FG  - n×2-matrix of laminated (f,g) pairs to test
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.
NB. Algorithm for calculating backward error:
NB.   In:  FG, CS
NB.   where
NB.        CS=. lartg FG
NB.        CS -: C ,. S
NB.   Out: maxberr
NB.   1) find R:
NB.        R=. CS mp"1 FG
NB.   2) find exact solution for each row by Algorithm 1 [1]:
NB.        CSexact=. algo1"1 FG
NB.      where
NB.        CSexact -: Cexact ,. Sexact
NB.   3) find Rexact:
NB.        Rexact=. CSexact mp"1 FG
NB.   4) combine R and Rexact:
NB.        Rboth=. R ,. Rexact
NB.   5) exclude rows containing NaNs:
NB.        Rboth=. xrNaN Rboth
NB.   6) exclude rows containing ±∞:
NB.        Rboth=. xrInf Rboth
NB.   7) calculate backward error for each pair
NB.      (R[i],Rexact[i]) [1]:
NB.        BErr[i] := |R[i] - Rexact[i]| / max(FP_EPS * |Rexact[i]|, FP_SFMIN * FP_PREC)
NB.      where
NB.        R[i]      - approximation computed by lartg,
NB.                      |R[i] - Rexact[i]| ≤ FP_OVFL
NB.        Rexact[i] - exact value computed by Algorithm 1
NB.                    [1],
NB.                      |Rexact[i]| ≤ FP_OVFL
NB.   8) exclude +∞ from vector BErr:
NB.        BErr=. xeInf BErr
NB.   9) find backward error:
NB.        maxberr = max(BErr[:])

testlartg=: 3 : 0
  NB. implement Algorithm 1 [1]
  algo1=: 3 : 'if. 0 = {: y do. (sgn , 0:) {. y elseif. 0 = {. y do. (0 , sgn@+) {: y else. try. ((| f),((sgn f) * (+ g))) % %: +/ soris ''f g''=. y catch. 2 # _. end. end.'

  NB. exclude rows containing NaN from the table y
  xrNaN=: #~ +:/"1@isnan

  NB. exclude rows containing ±∞ from the table y
  xrInf=: #~ -.@(+./)@|:@:(_ = |)

  NB. exclude elements ±∞ from the table y
  xeInf=: #~ _ ~: |

  NB. backward error calculator:
  vberrlartg=: (mp algo1)"1@[ (|@:- >./@:xeInf@:% (FP_SFMIN * FP_PREC) >. FP_EPS * |@[)/@|:@xrInf@xrNaN@,. mp"1

  log=. ('lartg' tmonad (]`]`(_."_)`(_."_)`vberrlartg)) y

  erase 'algo1 xrNaN xrInf xeInf vberrlartg'

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
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by 200 random real 2-vectors with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testrot_mt_ 200 150
NB. - test by 200 random real 2-vectors with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testrot_mt_ 200 200
NB. - test by 150 random complex 2-vectors:
NB.     log=. (gemat_mt_ j. gemat_mt_) testrot_mt_ 150 200

testrot=: 1 : 'testlartg_mt_@u@({. , 2:)'

NB. =========================================================
NB. Verification suite

NB. ---------------------------------------------------------
NB. verifyrot
NB.
NB. Description:
NB.   Nilad to verify rot actors, output result to console
NB.   and return it
NB.
NB. Syntax:
NB.   'probed failed'=. verifyrot ''
NB. where
NB.   probed ≥ 0, assertions probed counter
NB.   failed ≥ 0, assertions failed counter

verifyrot=: 3 : 0
  'q3 q4'=. ,.~ (j.~"0) 1r2 1r8
  sqr_big=. 2 *&%: FP_OVFL               NB. its square overflows
  cs=. qnsign , 0 1 1 ]`(j./)/. gemat 3  NB. random rotator

  NB. lartg
  NB. - input contains NaN
  res=.       fassert 1 1 -: isnan lartg  0     0j_.
  res=. res , fassert 1 1 -: isnan lartg  0    _.     NB. LAPACK's DLARTG gives (c=0,s=1,r=NaN), ZLARTG gives (c=0,s=NaN+i*NaN,r=NaN+i*0)
  res=. res , fassert 1 1 -: isnan lartg  0j_.  0
  res=. res , fassert 1 1 -: isnan lartg _.     0     NB. LAPACK's xLARTG gives (c=1,s=0,r=NaN)
  res=. res , fassert 1 1 -: isnan lartg  0    _.j_.
  res=. res , fassert 1 1 -: isnan lartg  0j_.  0j_.
  res=. res , fassert 1 1 -: isnan lartg  0j_. _.
  res=. res , fassert 1 1 -: isnan lartg _.     0j_.
  res=. res , fassert 1 1 -: isnan lartg _.    _.
  res=. res , fassert 1 1 -: isnan lartg _.j_.  0
  res=. res , fassert 1 1 -: isnan lartg _.j_. _.
  res=. res , fassert 1 1 -: isnan lartg _.j_.  0j_.
  res=. res , fassert 1 1 -: isnan lartg _.    _.j_.
  res=. res , fassert 1 1 -: isnan lartg  0j_. _.j_.
  res=. res , fassert 1 1 -: isnan lartg _.j_. _.j_.
  res=. res , fassert (2 2 $ 1) -: isnan lartg _. 0 ,: 0 _.
  res=. res , fassert (,.~ 1 0) -: isnan lartg _. 0 ,: q3
  res=. res , fassert (,.~ 0 1) -: isnan lartg q3   ,: 0 _.
  NB. - input contains q∞ so NaN error must be throwed
  res=. res , fassert 0:@lartg :: 1: __j__ __j__
  res=. res , fassert 0:@lartg :: 1: __j__ __j_1
  res=. res , fassert 0:@lartg :: 1: __j__ __
  res=. res , fassert 0:@lartg :: 1: __j__ __j1
  res=. res , fassert 0:@lartg :: 1: __j__ __j_
  res=. res , fassert 0:@lartg :: 1: __j__ _1j__
  res=. res , fassert 0:@lartg :: 1: __j__ _1j_1
  res=. res , fassert 0:@lartg :: 1: __j__ _1
  res=. res , fassert 0:@lartg :: 1: __j__ _1j1
  res=. res , fassert 0:@lartg :: 1: __j__ _1j_
  res=. res , fassert 0:@lartg :: 1: __j__  0j__
  res=. res , fassert 0:@lartg :: 1: __j__  0j_1
  res=. res , fassert 0:@lartg :: 1: __j__  0
  res=. res , fassert 0:@lartg :: 1: __j__  0j1
  res=. res , fassert 0:@lartg :: 1: __j__  0j_
  res=. res , fassert 0:@lartg :: 1: __j__  1j__
  res=. res , fassert 0:@lartg :: 1: __j__  1j_1
  res=. res , fassert 0:@lartg :: 1: __j__  1
  res=. res , fassert 0:@lartg :: 1: __j__  1j1
  res=. res , fassert 0:@lartg :: 1: __j__  1j_
  res=. res , fassert 0:@lartg :: 1: __j__  _j__
  res=. res , fassert 0:@lartg :: 1: __j__  _j_1
  res=. res , fassert 0:@lartg :: 1: __j__  _
  res=. res , fassert 0:@lartg :: 1: __j__  _j1
  res=. res , fassert 0:@lartg :: 1: __j__  _j_
  res=. res , fassert 0:@lartg :: 1: __j_1 __j__
  res=. res , fassert 0:@lartg :: 1: __j_1 __j_1
  res=. res , fassert 0:@lartg :: 1: __j_1 __
  res=. res , fassert 0:@lartg :: 1: __j_1 __j1
  res=. res , fassert 0:@lartg :: 1: __j_1 __j_
  res=. res , fassert 0:@lartg :: 1: __j_1 _1j__
  res=. res , fassert 0:@lartg :: 1: __j_1 _1j_
  res=. res , fassert 0:@lartg :: 1: __j_1  0j__
  res=. res , fassert 0:@lartg :: 1: __j_1  0j_
  res=. res , fassert 0:@lartg :: 1: __j_1  1j__
  res=. res , fassert 0:@lartg :: 1: __j_1  1j_
  res=. res , fassert 0:@lartg :: 1: __j_1  _j__
  res=. res , fassert 0:@lartg :: 1: __j_1  _j_1
  res=. res , fassert 0:@lartg :: 1: __j_1  _
  res=. res , fassert 0:@lartg :: 1: __j_1  _j1
  res=. res , fassert 0:@lartg :: 1: __j_1  _j_
  res=. res , fassert 0:@lartg :: 1: __    __j__
  res=. res , fassert 0:@lartg :: 1: __    __j_1
  res=. res , fassert 0:@lartg :: 1: __    __
  res=. res , fassert 0:@lartg :: 1: __    __j1
  res=. res , fassert 0:@lartg :: 1: __    __j_
  res=. res , fassert 0:@lartg :: 1: __    _1j__
  res=. res , fassert 0:@lartg :: 1: __    _1j_
  res=. res , fassert 0:@lartg :: 1: __     0j__
  res=. res , fassert 0:@lartg :: 1: __     0j_
  res=. res , fassert 0:@lartg :: 1: __     1j__
  res=. res , fassert 0:@lartg :: 1: __     1j_
  res=. res , fassert 0:@lartg :: 1: __     _j__
  res=. res , fassert 0:@lartg :: 1: __     _j_1
  res=. res , fassert 0:@lartg :: 1: __     _
  res=. res , fassert 0:@lartg :: 1: __     _j1
  res=. res , fassert 0:@lartg :: 1: __     _j_
  res=. res , fassert 0:@lartg :: 1: __j1  __j__
  res=. res , fassert 0:@lartg :: 1: __j1  __j_1
  res=. res , fassert 0:@lartg :: 1: __j1  __
  res=. res , fassert 0:@lartg :: 1: __j1  __j1
  res=. res , fassert 0:@lartg :: 1: __j1  __j_
  res=. res , fassert 0:@lartg :: 1: __j1  _1j__
  res=. res , fassert 0:@lartg :: 1: __j1  _1j_
  res=. res , fassert 0:@lartg :: 1: __j1   0j__
  res=. res , fassert 0:@lartg :: 1: __j1   0j_
  res=. res , fassert 0:@lartg :: 1: __j1   1j__
  res=. res , fassert 0:@lartg :: 1: __j1   1j_
  res=. res , fassert 0:@lartg :: 1: __j1   _j__
  res=. res , fassert 0:@lartg :: 1: __j1   _j_1
  res=. res , fassert 0:@lartg :: 1: __j1   _
  res=. res , fassert 0:@lartg :: 1: __j1   _j1
  res=. res , fassert 0:@lartg :: 1: __j1   _j_
  res=. res , fassert 0:@lartg :: 1: __j_  __j__
  res=. res , fassert 0:@lartg :: 1: __j_  __j_1
  res=. res , fassert 0:@lartg :: 1: __j_  __
  res=. res , fassert 0:@lartg :: 1: __j_  __j1
  res=. res , fassert 0:@lartg :: 1: __j_  __j_
  res=. res , fassert 0:@lartg :: 1: __j_  _1j__
  res=. res , fassert 0:@lartg :: 1: __j_  _1j_1
  res=. res , fassert 0:@lartg :: 1: __j_  _1
  res=. res , fassert 0:@lartg :: 1: __j_  _1j1
  res=. res , fassert 0:@lartg :: 1: __j_  _1j_
  res=. res , fassert 0:@lartg :: 1: __j_   0j__
  res=. res , fassert 0:@lartg :: 1: __j_   0j_1
  res=. res , fassert 0:@lartg :: 1: __j_   0
  res=. res , fassert 0:@lartg :: 1: __j_   0j1
  res=. res , fassert 0:@lartg :: 1: __j_   0j_
  res=. res , fassert 0:@lartg :: 1: __j_   1j__
  res=. res , fassert 0:@lartg :: 1: __j_   1j_1
  res=. res , fassert 0:@lartg :: 1: __j_   1
  res=. res , fassert 0:@lartg :: 1: __j_   1j1
  res=. res , fassert 0:@lartg :: 1: __j_   1j_
  res=. res , fassert 0:@lartg :: 1: __j_   _j__
  res=. res , fassert 0:@lartg :: 1: __j_   _j_1
  res=. res , fassert 0:@lartg :: 1: __j_   _
  res=. res , fassert 0:@lartg :: 1: __j_   _j1
  res=. res , fassert 0:@lartg :: 1: __j_   _j_
  res=. res , fassert 0:@lartg :: 1: _1j__ __j__
  res=. res , fassert 0:@lartg :: 1: _1j__ __j_1
  res=. res , fassert 0:@lartg :: 1: _1j__ __
  res=. res , fassert 0:@lartg :: 1: _1j__ __j1
  res=. res , fassert 0:@lartg :: 1: _1j__ __j_
  res=. res , fassert 0:@lartg :: 1: _1j__ _1j__
  res=. res , fassert 0:@lartg :: 1: _1j__ _1j_
  res=. res , fassert 0:@lartg :: 1: _1j__  0j__
  res=. res , fassert 0:@lartg :: 1: _1j__  0j_
  res=. res , fassert 0:@lartg :: 1: _1j__  1j__
  res=. res , fassert 0:@lartg :: 1: _1j__  1j_
  res=. res , fassert 0:@lartg :: 1: _1j__  _j__
  res=. res , fassert 0:@lartg :: 1: _1j__  _j_1
  res=. res , fassert 0:@lartg :: 1: _1j__  _
  res=. res , fassert 0:@lartg :: 1: _1j__  _j1
  res=. res , fassert 0:@lartg :: 1: _1j__  _j_
  res=. res , fassert 0:@lartg :: 1: _1j_1 __j__
  res=. res , fassert 0:@lartg :: 1: _1j_1 __j_
  res=. res , fassert 0:@lartg :: 1: _1j_1  _j__
  res=. res , fassert 0:@lartg :: 1: _1j_1  _j_
  res=. res , fassert 0:@lartg :: 1: _1    __j__
  res=. res , fassert 0:@lartg :: 1: _1    __j_
  res=. res , fassert 0:@lartg :: 1: _1     _j__
  res=. res , fassert 0:@lartg :: 1: _1     _j_
  res=. res , fassert 0:@lartg :: 1: _1j1  __j__
  res=. res , fassert 0:@lartg :: 1: _1j1  __j_
  res=. res , fassert 0:@lartg :: 1: _1j1   _j__
  res=. res , fassert 0:@lartg :: 1: _1j1   _j_
  res=. res , fassert 0:@lartg :: 1: _1j_  __j__
  res=. res , fassert 0:@lartg :: 1: _1j_  __j_1
  res=. res , fassert 0:@lartg :: 1: _1j_  __
  res=. res , fassert 0:@lartg :: 1: _1j_  __j1
  res=. res , fassert 0:@lartg :: 1: _1j_  __j_
  res=. res , fassert 0:@lartg :: 1: _1j_  _1j__
  res=. res , fassert 0:@lartg :: 1: _1j_  _1j_
  res=. res , fassert 0:@lartg :: 1: _1j_   0j__
  res=. res , fassert 0:@lartg :: 1: _1j_   0j_
  res=. res , fassert 0:@lartg :: 1: _1j_   1j__
  res=. res , fassert 0:@lartg :: 1: _1j_   1j_
  res=. res , fassert 0:@lartg :: 1: _1j_   _j__
  res=. res , fassert 0:@lartg :: 1: _1j_   _j_1
  res=. res , fassert 0:@lartg :: 1: _1j_   _
  res=. res , fassert 0:@lartg :: 1: _1j_   _j1
  res=. res , fassert 0:@lartg :: 1: _1j_   _j_
  res=. res , fassert 0:@lartg :: 1:  0j__ __j__
  res=. res , fassert 0:@lartg :: 1:  0j__ __j_1
  res=. res , fassert 0:@lartg :: 1:  0j__ __
  res=. res , fassert 0:@lartg :: 1:  0j__ __j1
  res=. res , fassert 0:@lartg :: 1:  0j__ __j_
  res=. res , fassert 0:@lartg :: 1:  0j__ _1j__
  res=. res , fassert 0:@lartg :: 1:  0j__ _1j_
  res=. res , fassert 0:@lartg :: 1:  0j__  0j__
  res=. res , fassert 0:@lartg :: 1:  0j__  0j_
  res=. res , fassert 0:@lartg :: 1:  0j__  1j__
  res=. res , fassert 0:@lartg :: 1:  0j__  1j_
  res=. res , fassert 0:@lartg :: 1:  0j__  _j__
  res=. res , fassert 0:@lartg :: 1:  0j__  _j_1
  res=. res , fassert 0:@lartg :: 1:  0j__  _
  res=. res , fassert 0:@lartg :: 1:  0j__  _j1
  res=. res , fassert 0:@lartg :: 1:  0j__  _j_
  res=. res , fassert 0:@lartg :: 1:  0j_1 __j__
  res=. res , fassert 0:@lartg :: 1:  0j_1 __j_
  res=. res , fassert 0:@lartg :: 1:  0j_1  _j__
  res=. res , fassert 0:@lartg :: 1:  0j_1  _j_
  res=. res , fassert 0:@lartg :: 1:  0    __j__
  res=. res , fassert 0:@lartg :: 1:  0    __j_
  res=. res , fassert 0:@lartg :: 1:  0     _j__
  res=. res , fassert 0:@lartg :: 1:  0     _j_
  res=. res , fassert 0:@lartg :: 1:  0j1  __j__
  res=. res , fassert 0:@lartg :: 1:  0j1  __j_
  res=. res , fassert 0:@lartg :: 1:  0j1   _j__
  res=. res , fassert 0:@lartg :: 1:  0j1   _j_
  res=. res , fassert 0:@lartg :: 1:  0j_  __j__
  res=. res , fassert 0:@lartg :: 1:  0j_  __j_1
  res=. res , fassert 0:@lartg :: 1:  0j_  __
  res=. res , fassert 0:@lartg :: 1:  0j_  __j1
  res=. res , fassert 0:@lartg :: 1:  0j_  __j_
  res=. res , fassert 0:@lartg :: 1:  0j_  _1j__
  res=. res , fassert 0:@lartg :: 1:  0j_  _1j_
  res=. res , fassert 0:@lartg :: 1:  0j_   0j__
  res=. res , fassert 0:@lartg :: 1:  0j_   0j_
  res=. res , fassert 0:@lartg :: 1:  0j_   1j__
  res=. res , fassert 0:@lartg :: 1:  0j_   1j_
  res=. res , fassert 0:@lartg :: 1:  0j_   _j__
  res=. res , fassert 0:@lartg :: 1:  0j_   _j_1
  res=. res , fassert 0:@lartg :: 1:  0j_   _
  res=. res , fassert 0:@lartg :: 1:  0j_   _j1
  res=. res , fassert 0:@lartg :: 1:  0j_   _j_
  res=. res , fassert 0:@lartg :: 1:  1j__ __j__
  res=. res , fassert 0:@lartg :: 1:  1j__ __j_1
  res=. res , fassert 0:@lartg :: 1:  1j__ __
  res=. res , fassert 0:@lartg :: 1:  1j__ __j1
  res=. res , fassert 0:@lartg :: 1:  1j__ __j_
  res=. res , fassert 0:@lartg :: 1:  1j__ _1j__
  res=. res , fassert 0:@lartg :: 1:  1j__ _1j_
  res=. res , fassert 0:@lartg :: 1:  1j__  0j__
  res=. res , fassert 0:@lartg :: 1:  1j__  0j_
  res=. res , fassert 0:@lartg :: 1:  1j__  1j__
  res=. res , fassert 0:@lartg :: 1:  1j__  1j_
  res=. res , fassert 0:@lartg :: 1:  1j__  _j__
  res=. res , fassert 0:@lartg :: 1:  1j__  _j_1
  res=. res , fassert 0:@lartg :: 1:  1j__  _
  res=. res , fassert 0:@lartg :: 1:  1j__  _j1
  res=. res , fassert 0:@lartg :: 1:  1j__  _j_
  res=. res , fassert 0:@lartg :: 1:  1j_1 __j__
  res=. res , fassert 0:@lartg :: 1:  1j_1 __j_
  res=. res , fassert 0:@lartg :: 1:  1j_1  _j__
  res=. res , fassert 0:@lartg :: 1:  1j_1  _j_
  res=. res , fassert 0:@lartg :: 1:  1    __j__
  res=. res , fassert 0:@lartg :: 1:  1    __j_
  res=. res , fassert 0:@lartg :: 1:  1     _j__
  res=. res , fassert 0:@lartg :: 1:  1     _j_
  res=. res , fassert 0:@lartg :: 1:  1j1  __j__
  res=. res , fassert 0:@lartg :: 1:  1j1  __j_
  res=. res , fassert 0:@lartg :: 1:  1j1   _j__
  res=. res , fassert 0:@lartg :: 1:  1j1   _j_
  res=. res , fassert 0:@lartg :: 1:  1j_  __j__
  res=. res , fassert 0:@lartg :: 1:  1j_  __j_1
  res=. res , fassert 0:@lartg :: 1:  1j_  __
  res=. res , fassert 0:@lartg :: 1:  1j_  __j1
  res=. res , fassert 0:@lartg :: 1:  1j_  __j_
  res=. res , fassert 0:@lartg :: 1:  1j_  _1j__
  res=. res , fassert 0:@lartg :: 1:  1j_  _1j_
  res=. res , fassert 0:@lartg :: 1:  1j_   0j__
  res=. res , fassert 0:@lartg :: 1:  1j_   0j_
  res=. res , fassert 0:@lartg :: 1:  1j_   1j__
  res=. res , fassert 0:@lartg :: 1:  1j_   1j_
  res=. res , fassert 0:@lartg :: 1:  1j_   _j__
  res=. res , fassert 0:@lartg :: 1:  1j_   _j_1
  res=. res , fassert 0:@lartg :: 1:  1j_   _
  res=. res , fassert 0:@lartg :: 1:  1j_   _j1
  res=. res , fassert 0:@lartg :: 1:  1j_   _j_
  res=. res , fassert 0:@lartg :: 1:  _j__ __j__
  res=. res , fassert 0:@lartg :: 1:  _j__ __j_1
  res=. res , fassert 0:@lartg :: 1:  _j__ __
  res=. res , fassert 0:@lartg :: 1:  _j__ __j1
  res=. res , fassert 0:@lartg :: 1:  _j__ __j_
  res=. res , fassert 0:@lartg :: 1:  _j__ _1j__
  res=. res , fassert 0:@lartg :: 1:  _j__ _1j_1
  res=. res , fassert 0:@lartg :: 1:  _j__ _1
  res=. res , fassert 0:@lartg :: 1:  _j__ _1j1
  res=. res , fassert 0:@lartg :: 1:  _j__ _1j_
  res=. res , fassert 0:@lartg :: 1:  _j__  0j__
  res=. res , fassert 0:@lartg :: 1:  _j__  0j_1
  res=. res , fassert 0:@lartg :: 1:  _j__  0
  res=. res , fassert 0:@lartg :: 1:  _j__  0j1
  res=. res , fassert 0:@lartg :: 1:  _j__  0j_
  res=. res , fassert 0:@lartg :: 1:  _j__  1j__
  res=. res , fassert 0:@lartg :: 1:  _j__  1j_1
  res=. res , fassert 0:@lartg :: 1:  _j__  1
  res=. res , fassert 0:@lartg :: 1:  _j__  1j1
  res=. res , fassert 0:@lartg :: 1:  _j__  1j_
  res=. res , fassert 0:@lartg :: 1:  _j__  _j__
  res=. res , fassert 0:@lartg :: 1:  _j__  _j_1
  res=. res , fassert 0:@lartg :: 1:  _j__  _
  res=. res , fassert 0:@lartg :: 1:  _j__  _j1
  res=. res , fassert 0:@lartg :: 1:  _j__  _j_
  res=. res , fassert 0:@lartg :: 1:  _j_1 __j__
  res=. res , fassert 0:@lartg :: 1:  _j_1 __j_1
  res=. res , fassert 0:@lartg :: 1:  _j_1 __
  res=. res , fassert 0:@lartg :: 1:  _j_1 __j1
  res=. res , fassert 0:@lartg :: 1:  _j_1 __j_
  res=. res , fassert 0:@lartg :: 1:  _j_1 _1j__
  res=. res , fassert 0:@lartg :: 1:  _j_1 _1j_
  res=. res , fassert 0:@lartg :: 1:  _j_1  0j__
  res=. res , fassert 0:@lartg :: 1:  _j_1  0j_
  res=. res , fassert 0:@lartg :: 1:  _j_1  1j__
  res=. res , fassert 0:@lartg :: 1:  _j_1  1j_
  res=. res , fassert 0:@lartg :: 1:  _j_1  _j__
  res=. res , fassert 0:@lartg :: 1:  _j_1  _j_1
  res=. res , fassert 0:@lartg :: 1:  _j_1  _
  res=. res , fassert 0:@lartg :: 1:  _j_1  _j1
  res=. res , fassert 0:@lartg :: 1:  _j_1  _j_
  res=. res , fassert 0:@lartg :: 1:  _    __j__
  res=. res , fassert 0:@lartg :: 1:  _    __j_1
  res=. res , fassert 0:@lartg :: 1:  _    __
  res=. res , fassert 0:@lartg :: 1:  _    __j1
  res=. res , fassert 0:@lartg :: 1:  _    __j_
  res=. res , fassert 0:@lartg :: 1:  _    _1j__
  res=. res , fassert 0:@lartg :: 1:  _    _1j_
  res=. res , fassert 0:@lartg :: 1:  _     0j__
  res=. res , fassert 0:@lartg :: 1:  _     0j_
  res=. res , fassert 0:@lartg :: 1:  _     1j__
  res=. res , fassert 0:@lartg :: 1:  _     1j_
  res=. res , fassert 0:@lartg :: 1:  _     _j__
  res=. res , fassert 0:@lartg :: 1:  _     _j_1
  res=. res , fassert 0:@lartg :: 1:  _     _
  res=. res , fassert 0:@lartg :: 1:  _     _j1
  res=. res , fassert 0:@lartg :: 1:  _     _j_
  res=. res , fassert 0:@lartg :: 1:  _j1  __j__
  res=. res , fassert 0:@lartg :: 1:  _j1  __j_1
  res=. res , fassert 0:@lartg :: 1:  _j1  __
  res=. res , fassert 0:@lartg :: 1:  _j1  __j1
  res=. res , fassert 0:@lartg :: 1:  _j1  __j_
  res=. res , fassert 0:@lartg :: 1:  _j1  _1j__
  res=. res , fassert 0:@lartg :: 1:  _j1  _1j_
  res=. res , fassert 0:@lartg :: 1:  _j1   0j__
  res=. res , fassert 0:@lartg :: 1:  _j1   0j_
  res=. res , fassert 0:@lartg :: 1:  _j1   1j__
  res=. res , fassert 0:@lartg :: 1:  _j1   1j_
  res=. res , fassert 0:@lartg :: 1:  _j1   _j__
  res=. res , fassert 0:@lartg :: 1:  _j1   _j_1
  res=. res , fassert 0:@lartg :: 1:  _j1   _
  res=. res , fassert 0:@lartg :: 1:  _j1   _j1
  res=. res , fassert 0:@lartg :: 1:  _j1   _j_
  res=. res , fassert 0:@lartg :: 1:  _j_  __j__
  res=. res , fassert 0:@lartg :: 1:  _j_  __j_1
  res=. res , fassert 0:@lartg :: 1:  _j_  __
  res=. res , fassert 0:@lartg :: 1:  _j_  __j1
  res=. res , fassert 0:@lartg :: 1:  _j_  __j_
  res=. res , fassert 0:@lartg :: 1:  _j_  _1j__
  res=. res , fassert 0:@lartg :: 1:  _j_  _1j_1
  res=. res , fassert 0:@lartg :: 1:  _j_  _1
  res=. res , fassert 0:@lartg :: 1:  _j_  _1j1
  res=. res , fassert 0:@lartg :: 1:  _j_  _1j_
  res=. res , fassert 0:@lartg :: 1:  _j_   0j__
  res=. res , fassert 0:@lartg :: 1:  _j_   0j_1
  res=. res , fassert 0:@lartg :: 1:  _j_   0
  res=. res , fassert 0:@lartg :: 1:  _j_   0j1
  res=. res , fassert 0:@lartg :: 1:  _j_   0j_
  res=. res , fassert 0:@lartg :: 1:  _j_   1j__
  res=. res , fassert 0:@lartg :: 1:  _j_   1j_1
  res=. res , fassert 0:@lartg :: 1:  _j_   1
  res=. res , fassert 0:@lartg :: 1:  _j_   1j1
  res=. res , fassert 0:@lartg :: 1:  _j_   1j_
  res=. res , fassert 0:@lartg :: 1:  _j_   _j__
  res=. res , fassert 0:@lartg :: 1:  _j_   _j_1
  res=. res , fassert 0:@lartg :: 1:  _j_   _
  res=. res , fassert 0:@lartg :: 1:  _j_   _j1
  res=. res , fassert 0:@lartg :: 1:  _j_   _j_
  NB. - input contains directed infinity and is not trivial
  NB.   so NaN error must be throwed
  res=. res , fassert 0:@lartg :: 1  ( __        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j. -FP_OVFL) ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j. -FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j. -FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j. -FP_OVFL) ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1:   __                     , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1:   __                     ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1:   __                     , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1:   __                     ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1:   __                     ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1:   __                     ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1:   __                     ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1:   __                     ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j.  FP_OVFL) ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j.  FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j.  FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j.  FP_OVFL) ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ( __        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. __      ) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. __      ) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. __      ) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,  __        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,  __
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,  __        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. __
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  _
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   0        j. __
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  _
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. __
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  _
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   _        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   _
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j. -FP_OVFL) ,   _        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              ,  __        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              ,  __
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              ,  __        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              , (-FP_OVFL) j. __
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              , (-FP_OVFL) j.  _
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              ,   0        j. __
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              ,   0        j.  _
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              ,   FP_OVFL  j. __
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              ,   FP_OVFL  j.  _
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              ,   _        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              ,   _
  res=. res , fassert 0:@lartg :: 1   (-FP_OVFL)              ,   _        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,  __        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,  __
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,  __        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. __
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  _
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   0        j. __
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  _
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. __
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  _
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   _        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   _
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  FP_OVFL) ,   _        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  _      ) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  _      ) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1  ((-FP_OVFL) j.  _      ) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. __      ) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. __      ) ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. __      ) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. __      ) ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. __      ) ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. __      ) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. __      ) ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. __      ) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) ,  __        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) ,  __
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) ,  __        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) , (-FP_OVFL) j. __
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) , (-FP_OVFL) j.  _
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) ,   0        j. __
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) ,   0        j.  _
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) ,   FP_OVFL  j. __
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) ,   FP_OVFL  j.  _
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) ,   _        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) ,   _
  res=. res , fassert 0:@lartg :: 1  (  0        j. -FP_OVFL) ,   _        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) ,  __        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) ,  __
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) ,  __        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) , (-FP_OVFL) j. __
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) , (-FP_OVFL) j.  _
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) ,   0        j. __
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) ,   0        j.  _
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) ,   FP_OVFL  j. __
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) ,   FP_OVFL  j.  _
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) ,   _        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) ,   _
  res=. res , fassert 0:@lartg :: 1  (  0        j.  FP_OVFL) ,   _        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  _      ) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  _      ) ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  _      ) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  _      ) ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  _      ) ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  _      ) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  _      ) ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  0        j.  _      ) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. __      ) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. __      ) ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. __      ) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. __      ) ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. __      ) ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. __      ) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. __      ) ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. __      ) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,  __        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,  __
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,  __        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. __
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  _
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   0        j. __
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   0        j.  _
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. __
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  _
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   _        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   _
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j. -FP_OVFL) ,   _        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               ,  __        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               ,  __
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               ,  __        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               , (-FP_OVFL) j. __
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               , (-FP_OVFL) j.  _
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               ,   0        j. __
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               ,   0        j.  _
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               ,   FP_OVFL  j. __
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               ,   FP_OVFL  j.  _
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               ,   _        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               ,   _
  res=. res , fassert 0:@lartg :: 1:    FP_OVFL               ,   _        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,  __        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,  __
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,  __        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. __
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  _
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   0        j. __
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   0        j.  _
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. __
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  _
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   _        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   _
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  FP_OVFL) ,   _        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  _      ) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  _      ) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  FP_OVFL  j.  _      ) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j. -FP_OVFL) ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j. -FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j. -FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j. -FP_OVFL) ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    _                     , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    _                     ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    _                     , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    _                     ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    _                     ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    _                     ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    _                     ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1:    _                     ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j.  FP_OVFL) ,  -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j.  FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j.  FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j.  FP_OVFL) ,   FP_OVFL
  res=. res , fassert 0:@lartg :: 1  (  _        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  NB. - input contains directed infinity and is trivial
  res=. res , fassert 1  0    -: lartg ( __        j. -FP_OVFL) ,    0
  res=. res , fassert 1  0    -: lartg ( __        j. _1      ) ,    0
  res=. res , fassert 1  0    -: lartg ( __        j. -FP_UNFL) ,    0
  res=. res , fassert 1  0    -: lartg ( __        j.  FP_UNFL) ,    0
  res=. res , fassert 1  0    -: lartg ( __        j.  1      ) ,    0
  res=. res , fassert 1  0    -: lartg ( __        j.  FP_OVFL) ,    0
  res=. res , fassert 1  0    -: lartg ((-FP_OVFL) j. __      ) ,    0
  res=. res , fassert 1  0    -: lartg ((-FP_OVFL) j.  _      ) ,    0
  res=. res , fassert 1  0    -: lartg ( _1        j. __      ) ,    0
  res=. res , fassert 1  0    -: lartg ( _1        j.  _      ) ,    0
  res=. res , fassert 1  0    -: lartg ((-FP_UNFL) j. __      ) ,    0
  res=. res , fassert 1  0    -: lartg ((-FP_UNFL) j.  _      ) ,    0
  res=. res , fassert 0 _1    -: lartg    0                     , ( __        j. -FP_OVFL)
  res=. res , fassert 0 _1    -: lartg    0                     , ( __        j. _1      )
  res=. res , fassert 0 _1    -: lartg    0                     , ( __        j. -FP_UNFL)
  res=. res , fassert 0 _1    -: lartg    0                     , ( __        j.  FP_UNFL)
  res=. res , fassert 0 _1    -: lartg    0                     , ( __        j.  1      )
  res=. res , fassert 0 _1    -: lartg    0                     , ( __        j.  FP_OVFL)
  res=. res , fassert 0  0j1  -: lartg    0                     , ((-FP_OVFL) j. __      )
  res=. res , fassert 0  0j_1 -: lartg    0                     , ((-FP_OVFL) j.  _      )
  res=. res , fassert 0  0j1  -: lartg    0                     , ( _1        j. __      )
  res=. res , fassert 0  0j_1 -: lartg    0                     , ( _1        j.  _      )
  res=. res , fassert 0  0j1  -: lartg    0                     , ((-FP_UNFL) j. __      )
  res=. res , fassert 0  0j_1 -: lartg    0                     , ((-FP_UNFL) j.  _      )
  res=. res , fassert 0  0j1  -: lartg    0                     , (  FP_UNFL  j. __      )
  res=. res , fassert 0  0j_1 -: lartg    0                     , (  FP_UNFL  j.  _      )
  res=. res , fassert 0  0j1  -: lartg    0                     , (  1        j. __      )
  res=. res , fassert 0  0j_1 -: lartg    0                     , (  1        j.  _      )
  res=. res , fassert 0  0j1  -: lartg    0                     , (  FP_OVFL  j. __      )
  res=. res , fassert 0  0j_1 -: lartg    0                     , (  FP_OVFL  j.  _      )
  res=. res , fassert 0  1    -: lartg    0                     , (  _        j. -FP_OVFL)
  res=. res , fassert 0  1    -: lartg    0                     , (  _        j. _1      )
  res=. res , fassert 0  1    -: lartg    0                     , (  _        j. -FP_UNFL)
  res=. res , fassert 0  1    -: lartg    0                     , (  _        j.  FP_UNFL)
  res=. res , fassert 0  1    -: lartg    0                     , (  _        j.  1      )
  res=. res , fassert 0  1    -: lartg    0                     , (  _        j.  FP_OVFL)
  res=. res , fassert 1  0    -: lartg (  FP_UNFL  j. __      ) ,    0
  res=. res , fassert 1  0    -: lartg (  FP_UNFL  j.  _      ) ,    0
  res=. res , fassert 1  0    -: lartg (  1        j. __      ) ,    0
  res=. res , fassert 1  0    -: lartg (  1        j.  _      ) ,    0
  res=. res , fassert 1  0    -: lartg (  FP_OVFL  j. __      ) ,    0
  res=. res , fassert 1  0    -: lartg (  FP_OVFL  j.  _      ) ,    0
  res=. res , fassert 1  0    -: lartg (  _        j. -FP_OVFL) ,    0
  res=. res , fassert 1  0    -: lartg (  _        j. _1      ) ,    0
  res=. res , fassert 1  0    -: lartg (  _        j. -FP_UNFL) ,    0
  res=. res , fassert 1  0    -: lartg (  _        j.  FP_UNFL) ,    0
  res=. res , fassert 1  0    -: lartg (  _        j.  1      ) ,    0
  res=. res , fassert 1  0    -: lartg (  _        j.  FP_OVFL) ,    0
  res=. res , fassert (0       1  ,: 1     0 ) -: lartg (0 , _ j. FP_OVFL) ,: (FP_OVFL j. __) , 0
  res=. res , fassert ((% %: 2 2) ,: 1     0 ) -: lartg q3                 ,: (FP_OVFL j. __) , 0
  res=. res , fassert (0       1  ,: % %: 2 2) -: lartg (0 , _ j. FP_OVFL) ,: q4
  NB. - input is imaginary infinity and consequently is
  NB.   trivial
  res=. res , fassert 1  0    -: lartg __     0
  res=. res , fassert 1  0    -: lartg  0j__  0
  res=. res , fassert 0 _1    -: lartg  0    __
  res=. res , fassert 0  0j1  -: lartg  0     0j__
  res=. res , fassert 0  0j_1 -: lartg  0     0j_
  res=. res , fassert 0  1    -: lartg  0     _
  res=. res , fassert 1  0    -: lartg  0j_   0
  res=. res , fassert 1  0    -: lartg  _     0
  NB. - input is trivial
  res=. res , fassert  1    0                      -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg ((-FP_OVFL) j. _1      ) ,   0
  res=. res , fassert  1    0                      -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg  (-FP_OVFL)              ,   0
  res=. res , fassert  1    0                      -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg ((-FP_OVFL) j.  1      ) ,   0
  res=. res , fassert  1    0                      -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg ((-sqr_big) j. -sqr_big) ,   0
  res=. res , fassert  1    0                      -: lartg  (-sqr_big)              ,   0
  res=. res , fassert  1    0                      -: lartg ((-sqr_big) j.  sqr_big) ,   0
  res=. res , fassert  1    0                      -: lartg ( _1        j. -FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg ( _1        j. _1      ) ,   0
  res=. res , fassert  1    0                      -: lartg ( _1        j. -FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg   _1                         0
  res=. res , fassert  1    0                      -: lartg ( _1        j.  FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg ( _1        j.  1      ) ,   0
  res=. res , fassert  1    0                      -: lartg ( _1        j.  FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg ((-FP_UNFL) j. _1      ) ,   0
  res=. res , fassert  1    0                      -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg  (-FP_UNFL)              ,   0
  res=. res , fassert  1    0                      -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg ((-FP_UNFL) j.  1      ) ,   0
  res=. res , fassert  1    0                      -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  0        j. -FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  0        j. -sqr_big) ,   0
  res=. res , fassert  1    0                      -: lartg (  0        j. _1      ) ,   0
  res=. res , fassert  1    0                      -: lartg (  0        j. -FP_UNFL) ,   0
  res=. res , fassert (0 , _1j1  % %: 2          ) -: lartg    0                     , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (0 , _1    j.   FP_UNFL % 4) -: lartg    0                     , (-FP_OVFL) j. _1        NB. LAPACK's DLARTG: (c=0,s=-1,r=FP_OVFL), ZLARTG: (c=0,s=-1+i*FP_UNFL/4,r=FP_OVFL)
  res=. res , fassert  0   _1                      -: lartg    0                     , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert  0   _1                      -: lartg    0                     ,  -FP_OVFL
  res=. res , fassert  0   _1                      -: lartg    0                     , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (0 , _1    j.   FP_UNFL %_4) -: lartg    0                     , (-FP_OVFL) j.  1        NB. LAPACK's DLARTG: (c=0,s=-1,r=FP_OVFL), ZLARTG: (c=0,s=-1-i*FP_UNFL/4,r=FP_OVFL)
  res=. res , fassert (0 , _1j_1 % %: 2          ) -: lartg    0                     , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (0 , _1j1  % %: 2          ) -: lartg    0                     , (-sqr_big) j. -sqr_big
  res=. res , fassert  0   _1                      -: lartg    0                     ,  -sqr_big
  res=. res , fassert (0 , _1j_1 % %: 2          ) -: lartg    0                     , (-sqr_big) j.  sqr_big
  res=. res , fassert (0 ,  1    j.~  FP_UNFL %_4) -: lartg    0                     ,  _1        j. -FP_OVFL
  res=. res , fassert (0 , _1j1  % %: 2          ) -: lartg    0                     ,  _1        j. _1
  res=. res , fassert (0 , _1    j.   FP_UNFL    ) -: lartg    0                     ,  _1        j. -FP_UNFL  NB. LAPACK's DLARTG: (c=0,s=-1,r=1), ZLARTG: (c=0,s=-1+i*FP_UNFL,r=1)
  res=. res , fassert  0   _1                      -: lartg    0                        _1
  res=. res , fassert (0 , _1    j.  -FP_UNFL    ) -: lartg    0                     ,  _1        j.  FP_UNFL  NB. LAPACK's DLARTG: (c=0,s=-1,r=1), ZLARTG: (c=0,s=-1-i*FP_UNFL,r=1)
  res=. res , fassert (0 , _1j_1 % %: 2          ) -: lartg    0                     ,  _1        j.  1
  res=. res , fassert (0 , _1    j.~  FP_UNFL %_4) -: lartg    0                     ,  _1        j.  FP_OVFL
  res=. res , fassert  0    0j1                    -: lartg    0                     , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert (0 ,  1    j.~ -FP_UNFL    ) -: lartg    0                     , (-FP_UNFL) j. _1
  res=. res , fassert (0 , _1j1  % %: 2          ) -: lartg    0                     , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  0   _1                      -: lartg    0                     ,  -FP_UNFL
  res=. res , fassert (0 , _1j_1 % %: 2          ) -: lartg    0                     , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (0 , _1    j.~ -FP_UNFL    ) -: lartg    0                     , (-FP_UNFL) j.  1
  res=. res , fassert  0    0j_1                   -: lartg    0                     , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert  0    0j1                    -: lartg    0                     ,   0        j. -FP_OVFL
  res=. res , fassert  0    0j1                    -: lartg    0                     ,   0        j. -sqr_big
  res=. res , fassert  0    0j1                    -: lartg    0                     ,   0        j. _1
  res=. res , fassert  0    0j1                    -: lartg    0                     ,   0        j. -FP_UNFL
  res=. res , fassert  1    0                      -: lartg    0                         0                     NB. test branching for zero input
  res=. res , fassert  0    0j_1                   -: lartg    0                     ,   0        j.  FP_UNFL
  res=. res , fassert  0    0j_1                   -: lartg    0                     ,   0        j.  1
  res=. res , fassert  0    0j_1                   -: lartg    0                     ,   0        j.  sqr_big
  res=. res , fassert  0    0j_1                   -: lartg    0                     ,   0        j.  FP_OVFL
  res=. res , fassert  0    0j1                    -: lartg    0                     ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert (0 ,  1    j.~  FP_UNFL    ) -: lartg    0                     ,   FP_UNFL  j. _1
  res=. res , fassert (0 ,  1j1   % %: 2         ) -: lartg    0                     ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  0    1                      -: lartg    0                     ,   FP_UNFL
  res=. res , fassert (0 ,  1j_1  % %: 2         ) -: lartg    0                     ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (0 , _1    j.~  FP_UNFL    ) -: lartg    0                     ,   FP_UNFL  j.  1
  res=. res , fassert  0    0j_1                   -: lartg    0                     ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (0 ,  1    j.~  FP_UNFL % 4) -: lartg    0                     ,   1        j. -FP_OVFL
  res=. res , fassert (0 ,  1j1   % %: 2         ) -: lartg    0                     ,   1        j. _1
  res=. res , fassert (0 ,  1    j.   FP_UNFL    ) -: lartg    0                     ,   1        j. -FP_UNFL  NB. LAPACK's DLARTG: (c=0,s= 1,r=1), ZLARTG: (c=0,s= 1+i*FP_UNFL,r=1)
  res=. res , fassert  0    1                      -: lartg    0                         1
  res=. res , fassert (0 ,  1    j.  -FP_UNFL    ) -: lartg    0                     ,   1        j.  FP_UNFL  NB. LAPACK's DLARTG: (c=0,s= 1,r=1), ZLARTG: (c=0,s= 1-i*FP_UNFL,r=1)
  res=. res , fassert (0 ,  1j_1  % %: 2         ) -: lartg    0                     ,   1        j.  1
  res=. res , fassert (0 , _1    j.~  FP_UNFL % 4) -: lartg    0                     ,   1        j.  FP_OVFL
  res=. res , fassert  0    1                      -: lartg    0                     ,   sqr_big
  res=. res , fassert (0 ,  1j1   % %: 2         ) -: lartg    0                     ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (0 ,  1    j.   FP_UNFL % 4) -: lartg    0                     ,   FP_OVFL  j. _1        NB. LAPACK's DLARTG: (c=0,s= 1,r=FP_OVFL), ZLARTG: (c=0,s= 1+i*FP_UNFL/4,r=FP_OVFL)
  res=. res , fassert  0    1                      -: lartg    0                     ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert  0    1                      -: lartg    0                     ,   FP_OVFL
  res=. res , fassert  0    1                      -: lartg    0                     ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (0 ,  1    j.   FP_UNFL %_4) -: lartg    0                     ,   FP_OVFL  j.  1        NB. LAPACK's DLARTG: (c=0,s= 1,r=FP_OVFL), ZLARTG: (c=0,s= 1-i*FP_UNFL/4,r=FP_OVFL)
  res=. res , fassert (0 ,  1j_1  % %: 2         ) -: lartg    0                     ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert  1    0                      -: lartg (  0        j.  FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  0        j.  1      ) ,   0
  res=. res , fassert  1    0                      -: lartg (  0        j.  FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_UNFL  j. -FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_UNFL  j. _1      ) ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_UNFL  j. -FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg    FP_UNFL               ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_UNFL  j.  FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_UNFL  j.  1      ) ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_UNFL  j.  FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  1        j. -FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  1        j. _1      ) ,   0
  res=. res , fassert  1    0                      -: lartg (  1        j. -FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg    1                         0
  res=. res , fassert  1    0                      -: lartg (  1        j.  FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  1        j.  1      ) ,   0
  res=. res , fassert  1    0                      -: lartg (  1        j.  FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  sqr_big  j. -sqr_big) ,   0
  res=. res , fassert  1    0                      -: lartg    sqr_big               ,   0
  res=. res , fassert  1    0                      -: lartg (  sqr_big  j.  sqr_big) ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_OVFL  j. -FP_OVFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_OVFL  j. _1      ) ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_OVFL  j. -FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg    FP_OVFL               ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_OVFL  j.  FP_UNFL) ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_OVFL  j.  1      ) ,   0
  res=. res , fassert  1    0                      -: lartg (  FP_OVFL  j.  FP_OVFL) ,   0
  NB. - edge cases input
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,  -FP_OVFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  FP_UNFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,  -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   0        j.  FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   0        j.  FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg ((-FP_OVFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg  (-FP_OVFL)              , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg  (-FP_OVFL)              , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg  (-FP_OVFL)              ,  -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg  (-FP_OVFL)              , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg  (-FP_OVFL)              , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1 ,   %FP_OVFL) -: lartg  (-FP_OVFL)              ,  _1
  res=. res , fassert (1  0j_1 % %: 2) -: lartg  (-FP_OVFL)              , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg  (-FP_OVFL)              , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg  (-FP_OVFL)              ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg  (-FP_OVFL)              , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg  (-FP_OVFL)              , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg  (-FP_OVFL)              ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg  (-FP_OVFL)              ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg  (-FP_OVFL)              ,   0        j.  FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg  (-FP_OVFL)              ,   0        j.  FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg  (-FP_OVFL)              ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg  (-FP_OVFL)              ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg  (-FP_OVFL)              ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg  (-FP_OVFL)              ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg  (-FP_OVFL)              ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1 , - %FP_OVFL) -: lartg  (-FP_OVFL)              ,   1
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg  (-FP_OVFL)              ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg  (-FP_OVFL)              ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg  (-FP_OVFL)              ,   FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg  (-FP_OVFL)              ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg  (-FP_OVFL)              ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,  -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   0        j.  FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   0        j.  FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg ((-FP_OVFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,  -FP_OVFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  FP_UNFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_OVFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-sqr_big) j. -sqr_big) , (-sqr_big) j. -sqr_big
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-sqr_big) j. -sqr_big) , (-sqr_big) j.  sqr_big
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-sqr_big) j. -sqr_big) ,   sqr_big  j. -sqr_big
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-sqr_big) j. -sqr_big) ,   sqr_big  j.  sqr_big
  res=. res , fassert (1  1    % %: 2) -: lartg  (-sqr_big)              ,  -sqr_big
  res=. res , fassert (1  0j_1 % %: 2) -: lartg  (-sqr_big)              ,   0        j. -sqr_big
  res=. res , fassert (1  0j1  % %: 2) -: lartg  (-sqr_big)              ,   0        j.  sqr_big
  res=. res , fassert (1 _1    % %: 2) -: lartg  (-sqr_big)              ,   sqr_big
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-sqr_big) j.  sqr_big) , (-sqr_big) j. -sqr_big
  res=. res , fassert (1  1    % %: 2) -: lartg ((-sqr_big) j.  sqr_big) , (-sqr_big) j.  sqr_big
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-sqr_big) j.  sqr_big) ,   sqr_big  j. -sqr_big
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-sqr_big) j.  sqr_big) ,   sqr_big  j.  sqr_big
  res=. res , fassert (1 ,~  %FP_OVFL) -: lartg   _1                     ,  -FP_OVFL
  res=. res , fassert (1 ,    FP_UNFL) -: lartg   _1                     ,  -FP_UNFL
  res=. res , fassert (1 , -  FP_UNFL) -: lartg   _1                     ,   FP_UNFL
  res=. res , fassert (_1,~  %FP_OVFL) -: lartg   _1                     ,   FP_OVFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,  -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   0        j.  FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg ((-FP_UNFL) j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert  0  1            -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,  -FP_OVFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert  0  0j1          -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,  -FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   0        j. -FP_OVFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   0        j. -FP_UNFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   0        j.  FP_UNFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   0        j.  FP_OVFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert  0  0j_1         -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert  0 _1            -: lartg ((-FP_UNFL) j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg  (-FP_UNFL)              , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert  0  1            -: lartg  (-FP_UNFL)              , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert  0  1            -: lartg  (-FP_UNFL)              ,  -FP_OVFL
  res=. res , fassert  0  1            -: lartg  (-FP_UNFL)              , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg  (-FP_UNFL)              , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1 ,~   FP_UNFL) -: lartg  (-FP_UNFL)              ,  _1
  res=. res , fassert  0  0j_1         -: lartg  (-FP_UNFL)              , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg  (-FP_UNFL)              , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg  (-FP_UNFL)              ,  -FP_UNFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg  (-FP_UNFL)              , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert  0  0j1          -: lartg  (-FP_UNFL)              , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert  0  0j_1         -: lartg  (-FP_UNFL)              ,   0        j. -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg  (-FP_UNFL)              ,   0        j. -FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg  (-FP_UNFL)              ,   0        j.  FP_UNFL
  res=. res , fassert  0  0j1          -: lartg  (-FP_UNFL)              ,   0        j.  FP_OVFL
  res=. res , fassert  0  0j_1         -: lartg  (-FP_UNFL)              ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg  (-FP_UNFL)              ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg  (-FP_UNFL)              ,   FP_UNFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg  (-FP_UNFL)              ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert  0  0j1          -: lartg  (-FP_UNFL)              ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (_1,~   FP_UNFL) -: lartg  (-FP_UNFL)              ,   1
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg  (-FP_UNFL)              ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert  0 _1            -: lartg  (-FP_UNFL)              ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert  0 _1            -: lartg  (-FP_UNFL)              ,   FP_OVFL
  res=. res , fassert  0 _1            -: lartg  (-FP_UNFL)              ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg  (-FP_UNFL)              ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert  0  0j_1         -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,  -FP_OVFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert  0  1            -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,  -FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   0        j. -FP_OVFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   0        j. -FP_UNFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   0        j.  FP_UNFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   0        j.  FP_OVFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert  0 _1            -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert  0  0j1          -: lartg ((-FP_UNFL) j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,  -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   0        j.  FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg ((-FP_UNFL) j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg (  0        j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  0        j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  0        j. -FP_OVFL) ,  -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  0        j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg (  0        j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  0        j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  0        j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  0        j. -FP_OVFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  0        j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  0        j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  0        j. -FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  0        j. -FP_OVFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  0        j. -FP_OVFL) ,   0        j.  FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  0        j. -FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  0        j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  0        j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  0        j. -FP_OVFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg (  0        j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  0        j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg (  0        j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  0        j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  0        j. -FP_OVFL) ,   FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  0        j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg (  0        j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  0        j. -sqr_big) ,  -sqr_big
  res=. res , fassert (1  1    % %: 2) -: lartg (  0        j. -sqr_big) ,   0        j. -sqr_big
  res=. res , fassert (1 _1    % %: 2) -: lartg (  0        j. -sqr_big) ,   0        j.  sqr_big
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  0        j. -sqr_big) ,   sqr_big
  res=. res , fassert (0  1j1  % %: 2) -: lartg (  0        j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert  0  0j1          -: lartg (  0        j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert  0  0j1          -: lartg (  0        j. -FP_UNFL) ,  -FP_OVFL
  res=. res , fassert  0  0j1          -: lartg (  0        j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg (  0        j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert  0  1            -: lartg (  0        j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg (  0        j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  0        j. -FP_UNFL) ,  -FP_UNFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg (  0        j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert  0 _1            -: lartg (  0        j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert  0  1            -: lartg (  0        j. -FP_UNFL) ,   0        j. -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  0        j. -FP_UNFL) ,   0        j. -FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  0        j. -FP_UNFL) ,   0        j.  FP_UNFL
  res=. res , fassert  0 _1            -: lartg (  0        j. -FP_UNFL) ,   0        j.  FP_OVFL
  res=. res , fassert  0  1            -: lartg (  0        j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg (  0        j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  0        j. -FP_UNFL) ,   FP_UNFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg (  0        j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert  0 _1            -: lartg (  0        j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg (  0        j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert  0  0j_1         -: lartg (  0        j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert  0  0j_1         -: lartg (  0        j. -FP_UNFL) ,   FP_OVFL
  res=. res , fassert  0  0j_1         -: lartg (  0        j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg (  0        j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg (  0        j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert  0  0j_1         -: lartg (  0        j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert  0  0j_1         -: lartg (  0        j.  FP_UNFL) ,  -FP_OVFL
  res=. res , fassert  0  0j_1         -: lartg (  0        j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg (  0        j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert  0 _1            -: lartg (  0        j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg (  0        j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  0        j.  FP_UNFL) ,  -FP_UNFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg (  0        j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert  0  1            -: lartg (  0        j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert  0 _1            -: lartg (  0        j.  FP_UNFL) ,   0        j. -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  0        j.  FP_UNFL) ,   0        j. -FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  0        j.  FP_UNFL) ,   0        j.  FP_UNFL
  res=. res , fassert  0  1            -: lartg (  0        j.  FP_UNFL) ,   0        j.  FP_OVFL
  res=. res , fassert  0 _1            -: lartg (  0        j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg (  0        j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  0        j.  FP_UNFL) ,   FP_UNFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg (  0        j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert  0  1            -: lartg (  0        j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg (  0        j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert  0  0j1          -: lartg (  0        j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert  0  0j1          -: lartg (  0        j.  FP_UNFL) ,   FP_OVFL
  res=. res , fassert  0  0j1          -: lartg (  0        j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg (  0        j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  0        j.  sqr_big) ,  -sqr_big
  res=. res , fassert (1 _1    % %: 2) -: lartg (  0        j.  sqr_big) ,   0        j. -sqr_big
  res=. res , fassert (1  1    % %: 2) -: lartg (  0        j.  sqr_big) ,   0        j.  sqr_big
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  0        j.  sqr_big) ,   sqr_big
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg (  0        j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  0        j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  0        j.  FP_OVFL) ,  -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  0        j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg (  0        j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  0        j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  0        j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  0        j.  FP_OVFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  0        j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  0        j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  0        j.  FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  0        j.  FP_OVFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  0        j.  FP_OVFL) ,   0        j.  FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  0        j.  FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  0        j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  0        j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  0        j.  FP_OVFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg (  0        j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  0        j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg (  0        j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  0        j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  0        j.  FP_OVFL) ,   FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  0        j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg (  0        j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,  -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,   0        j.  FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg (  FP_UNFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert  0  0j1          -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,  -FP_OVFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert  0 _1            -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg (  FP_UNFL  j. -FP_UNFL) ,  -FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   0        j. -FP_OVFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   0        j. -FP_UNFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   0        j.  FP_UNFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   0        j.  FP_OVFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert  0  1            -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert  0  0j_1         -: lartg (  FP_UNFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg    FP_UNFL               , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert  0 _1            -: lartg    FP_UNFL               , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert  0 _1            -: lartg    FP_UNFL               ,  -FP_OVFL
  res=. res , fassert  0 _1            -: lartg    FP_UNFL               , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg    FP_UNFL               , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (_1,~   FP_UNFL) -: lartg    FP_UNFL               ,  _1
  res=. res , fassert  0  0j1          -: lartg    FP_UNFL               , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg    FP_UNFL               , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg    FP_UNFL               ,  -FP_UNFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg    FP_UNFL               , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert  0  0j_1         -: lartg    FP_UNFL               , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert  0  0j1          -: lartg    FP_UNFL               ,   0        j. -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg    FP_UNFL               ,   0        j. -FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg    FP_UNFL               ,   0        j.  FP_UNFL
  res=. res , fassert  0  0j_1         -: lartg    FP_UNFL               ,   0        j.  FP_OVFL
  res=. res , fassert  0  0j1          -: lartg    FP_UNFL               ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg    FP_UNFL               ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg    FP_UNFL               ,   FP_UNFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg    FP_UNFL               ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert  0  0j_1         -: lartg    FP_UNFL               ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1 ,~   FP_UNFL) -: lartg    FP_UNFL               ,   1
  res=. res , fassert (0  1j1  % %: 2) -: lartg    FP_UNFL               ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert  0  1            -: lartg    FP_UNFL               ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert  0  1            -: lartg    FP_UNFL               ,   FP_OVFL
  res=. res , fassert  0  1            -: lartg    FP_UNFL               ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg    FP_UNFL               ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert  0 _1            -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,  -FP_OVFL
  res=. res , fassert (0 _1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert  0  0j_1         -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg (  FP_UNFL  j.  FP_UNFL) ,  -FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   0        j. -FP_OVFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   0        j. -FP_UNFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   0        j.  FP_UNFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   0        j.  FP_OVFL
  res=. res , fassert (0 _1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (0  1j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert  0  0j1          -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL
  res=. res , fassert (0  1j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert  0  1            -: lartg (  FP_UNFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,  -FP_OVFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,   0        j.  FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg (  FP_UNFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (_1,~  %FP_OVFL) -: lartg    1                     ,  -FP_OVFL
  res=. res , fassert ( 1, -  FP_UNFL) -: lartg    1                     ,  -FP_UNFL
  res=. res , fassert ( 1,    FP_UNFL) -: lartg    1                     ,   FP_UNFL
  res=. res , fassert ( 1,~  %FP_OVFL) -: lartg    1                     ,   FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  sqr_big  j. -sqr_big) , (-sqr_big) j. -sqr_big
  res=. res , fassert (1 _1    % %: 2) -: lartg (  sqr_big  j. -sqr_big) , (-sqr_big) j.  sqr_big
  res=. res , fassert (1  1    % %: 2) -: lartg (  sqr_big  j. -sqr_big) ,   sqr_big  j. -sqr_big
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  sqr_big  j. -sqr_big) ,   sqr_big  j.  sqr_big
  res=. res , fassert (1 _1    % %: 2) -: lartg    sqr_big               ,  -sqr_big
  res=. res , fassert (1  0j1  % %: 2) -: lartg    sqr_big               ,   0        j. -sqr_big
  res=. res , fassert (1  0j_1 % %: 2) -: lartg    sqr_big               ,   0        j.  sqr_big
  res=. res , fassert (1  1    % %: 2) -: lartg    sqr_big               ,   sqr_big
  res=. res , fassert (1 _1    % %: 2) -: lartg (  sqr_big  j.  sqr_big) , (-sqr_big) j. -sqr_big
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  sqr_big  j.  sqr_big) , (-sqr_big) j.  sqr_big
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  sqr_big  j.  sqr_big) ,   sqr_big  j. -sqr_big
  res=. res , fassert (1  1    % %: 2) -: lartg (  sqr_big  j.  sqr_big) ,   sqr_big  j.  sqr_big
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,  -FP_OVFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,   0        j.  FP_UNFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_OVFL  j. -FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,  -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,   0        j.  FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   0        j.  FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg (  FP_OVFL  j. -FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg    FP_OVFL               , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg    FP_OVFL               , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg    FP_OVFL               ,  -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg    FP_OVFL               , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg    FP_OVFL               , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1 , - %FP_OVFL) -: lartg    FP_OVFL               ,  _1
  res=. res , fassert (1  0j1  % %: 2) -: lartg    FP_OVFL               , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg    FP_OVFL               , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg    FP_OVFL               ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg    FP_OVFL               , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg    FP_OVFL               , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg    FP_OVFL               ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg    FP_OVFL               ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg    FP_OVFL               ,   0        j.  FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg    FP_OVFL               ,   0        j.  FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg    FP_OVFL               ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg    FP_OVFL               ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg    FP_OVFL               ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg    FP_OVFL               ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg    FP_OVFL               ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1 ,   %FP_OVFL) -: lartg    FP_OVFL               ,   1
  res=. res , fassert (1  1j1  % %: 3) -: lartg    FP_OVFL               ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg    FP_OVFL               ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg    FP_OVFL               ,   FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg    FP_OVFL               ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg    FP_OVFL               ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1 _1j1  % %: 3) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,  -FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,   0        j.  FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   0        j.  FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1  1j1  % %: 3) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1  1j_1 % %: 3) -: lartg (  FP_OVFL  j.  FP_UNFL) ,   FP_OVFL  j.  FP_OVFL
  res=. res , fassert (1 _1    % %: 2) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_OVFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j. -FP_UNFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,  -FP_OVFL
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_UNFL
  res=. res , fassert (1  0j_1 % %: 2) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_OVFL) j.  FP_OVFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,  -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_UNFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) , (-FP_UNFL) j.  FP_OVFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   0        j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,   0        j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,   0        j.  FP_UNFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   0        j.  FP_OVFL
  res=. res , fassert (2 _1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_OVFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j. -FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL
  res=. res , fassert  1  0            -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_UNFL
  res=. res , fassert (2  1j_1 % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_UNFL  j.  FP_OVFL
  res=. res , fassert (1  0j1  % %: 2) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_OVFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j. -FP_UNFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL
  res=. res , fassert (2  1j1  % %: 6) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_UNFL
  res=. res , fassert (1  1    % %: 2) -: lartg (  FP_OVFL  j.  FP_OVFL) ,   FP_OVFL  j.  FP_OVFL
  NB. - input without edge cases
  res=. res , fassert (1  1    % %: 2) -: lartg _1j_1 _1j_1
  res=. res , fassert (2  1j1  % %: 6) -: lartg _1j_1 _1
  res=. res , fassert (1  0j1  % %: 2) -: lartg _1j_1 _1j1
  res=. res , fassert (2  1j_1 % %: 6) -: lartg _1j_1  0j_1
  res=. res , fassert (2 _1j1  % %: 6) -: lartg _1j_1  0j1
  res=. res , fassert (1  0j_1 % %: 2) -: lartg _1j_1  1j_1
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg _1j_1  1
  res=. res , fassert (1 _1    % %: 2) -: lartg _1j_1  1j1
  res=. res , fassert (1  1j_1 % %: 3) -: lartg _1    _1j_1
  res=. res , fassert (1  1    % %: 2) -: lartg _1    _1
  res=. res , fassert (1  1j1  % %: 3) -: lartg _1    _1j1
  res=. res , fassert (1  0j_1 % %: 2) -: lartg _1     0j_1
  res=. res , fassert (1  0j1  % %: 2) -: lartg _1     0j1
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg _1     1j_1
  res=. res , fassert (1 _1    % %: 2) -: lartg _1     1
  res=. res , fassert (1 _1j1  % %: 3) -: lartg _1     1j1
  res=. res , fassert (1  0j_1 % %: 2) -: lartg _1j1  _1j_1
  res=. res , fassert (2  1j_1 % %: 6) -: lartg _1j1  _1
  res=. res , fassert (1  1    % %: 2) -: lartg _1j1  _1j1
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg _1j1   0j_1
  res=. res , fassert (2  1j1  % %: 6) -: lartg _1j1   0j1
  res=. res , fassert (1 _1    % %: 2) -: lartg _1j1   1j_1
  res=. res , fassert (2 _1j1  % %: 6) -: lartg _1j1   1
  res=. res , fassert (1  0j1  % %: 2) -: lartg _1j1   1j1
  res=. res , fassert (1  1j1  % %: 3) -: lartg  0j_1 _1j_1
  res=. res , fassert (1  0j1  % %: 2) -: lartg  0j_1 _1
  res=. res , fassert (1 _1j1  % %: 3) -: lartg  0j_1 _1j1
  res=. res , fassert (1  1    % %: 2) -: lartg  0j_1  0j_1
  res=. res , fassert (1 _1    % %: 2) -: lartg  0j_1  0j1
  res=. res , fassert (1  1j_1 % %: 3) -: lartg  0j_1  1j_1
  res=. res , fassert (1  0j_1 % %: 2) -: lartg  0j_1  1
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg  0j_1  1j1
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg  0j1  _1j_1
  res=. res , fassert (1  0j_1 % %: 2) -: lartg  0j1  _1
  res=. res , fassert (1  1j_1 % %: 3) -: lartg  0j1  _1j1
  res=. res , fassert (1 _1    % %: 2) -: lartg  0j1   0j_1
  res=. res , fassert (1  1    % %: 2) -: lartg  0j1   0j1
  res=. res , fassert (1 _1j1  % %: 3) -: lartg  0j1   1j_1
  res=. res , fassert (1  0j1  % %: 2) -: lartg  0j1   1
  res=. res , fassert (1  1j1  % %: 3) -: lartg  0j1   1j1
  res=. res , fassert (1  0j1  % %: 2) -: lartg  1j_1 _1j_1
  res=. res , fassert (2 _1j1  % %: 6) -: lartg  1j_1 _1
  res=. res , fassert (1 _1    % %: 2) -: lartg  1j_1 _1j1
  res=. res , fassert (2  1j1  % %: 6) -: lartg  1j_1  0j_1
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg  1j_1  0j1
  res=. res , fassert (1  1    % %: 2) -: lartg  1j_1  1j_1
  res=. res , fassert (2  1j_1 % %: 6) -: lartg  1j_1  1
  res=. res , fassert (1  0j_1 % %: 2) -: lartg  1j_1  1j1
  res=. res , fassert (1 _1j1  % %: 3) -: lartg  1    _1j_1
  res=. res , fassert (1 _1    % %: 2) -: lartg  1    _1
  res=. res , fassert (1 _1j_1 % %: 3) -: lartg  1    _1j1
  res=. res , fassert (1  0j1  % %: 2) -: lartg  1     0j_1
  res=. res , fassert (1  0j_1 % %: 2) -: lartg  1     0j1
  res=. res , fassert (1  1j1  % %: 3) -: lartg  1     1j_1
  res=. res , fassert (1  1    % %: 2) -: lartg  1     1
  res=. res , fassert (1  1j_1 % %: 3) -: lartg  1     1j1
  res=. res , fassert (1 _1    % %: 2) -: lartg  1j1  _1j_1
  res=. res , fassert (2 _1j_1 % %: 6) -: lartg  1j1  _1
  res=. res , fassert (1  0j_1 % %: 2) -: lartg  1j1  _1j1
  res=. res , fassert (2 _1j1  % %: 6) -: lartg  1j1   0j_1
  res=. res , fassert (2  1j_1 % %: 6) -: lartg  1j1   0j1
  res=. res , fassert (1  0j1  % %: 2) -: lartg  1j1   1j_1
  res=. res , fassert (2  1j1  % %: 6) -: lartg  1j1   1
  res=. res , fassert (1  1    % %: 2) -: lartg  1j1   1j1
  res=. res , fassert (2 2 $   % %: 2) -: lartg q3 ,: q4
  NB. c is real
  res=. res , fassert 0 -:!.0 qni lartg _1j_1 _1j_1
  res=. res , fassert 0 -:!.0 qni lartg _1j_1 _1j1
  res=. res , fassert 0 -:!.0 qni lartg _1j_1  1j_1
  res=. res , fassert 0 -:!.0 qni lartg _1j_1  1j1
  res=. res , fassert 0 -:!.0 qni lartg _1j1  _1j_1
  res=. res , fassert 0 -:!.0 qni lartg _1j1  _1j1
  res=. res , fassert 0 -:!.0 qni lartg _1j1   1j_1
  res=. res , fassert 0 -:!.0 qni lartg _1j1   1j1
  res=. res , fassert 0 -:!.0 qni lartg  1j_1 _1j_1
  res=. res , fassert 0 -:!.0 qni lartg  1j_1 _1j1
  res=. res , fassert 0 -:!.0 qni lartg  1j_1  1j_1
  res=. res , fassert 0 -:!.0 qni lartg  1j_1  1j1
  res=. res , fassert 0 -:!.0 qni lartg  1j1  _1j_1
  res=. res , fassert 0 -:!.0 qni lartg  1j1  _1j1
  res=. res , fassert 0 -:!.0 qni lartg  1j1   1j_1
  res=. res , fassert 0 -:!.0 qni lartg  1j1   1j1

  NB. rot
  NB. - input contains NaN
  res=. res , fassert 1 1 -: isnan  0     0j_. rot  0     0j_.
  res=. res , fassert 1 1 -: isnan  0     0j_. rot  0    _.
  res=. res , fassert 1 1 -: isnan  0     0j_. rot  0j_.  0
  res=. res , fassert 1 1 -: isnan  0     0j_. rot _.     0
  res=. res , fassert 1 1 -: isnan  0    _.    rot  0     0j_.
  res=. res , fassert 1 1 -: isnan  0    _.    rot  0    _.
  res=. res , fassert 1 1 -: isnan  0    _.    rot  0j_.  0
  res=. res , fassert 1 1 -: isnan  0    _.    rot _.     0
  res=. res , fassert 1 1 -: isnan  0j_.  0    rot  0     0j_.
  res=. res , fassert 1 1 -: isnan  0j_.  0    rot  0    _.
  res=. res , fassert 1 1 -: isnan  0j_.  0    rot  0j_.  0
  res=. res , fassert 1 1 -: isnan  0j_.  0    rot _.     0
  res=. res , fassert 1 1 -: isnan _.     0    rot  0     0j_.
  res=. res , fassert 1 1 -: isnan _.     0    rot  0    _.
  res=. res , fassert 1 1 -: isnan _.     0    rot  0j_.  0
  res=. res , fassert 1 1 -: isnan _.     0    rot _.     0
  res=. res , fassert 1 1 -: isnan  0    _.j_. rot  0    _.j_.
  res=. res , fassert 1 1 -: isnan  0    _.j_. rot  0j_.  0j_.
  res=. res , fassert 1 1 -: isnan  0    _.j_. rot  0j_. _.
  res=. res , fassert 1 1 -: isnan  0    _.j_. rot _.     0j_.
  res=. res , fassert 1 1 -: isnan  0    _.j_. rot _.    _.
  res=. res , fassert 1 1 -: isnan  0    _.j_. rot _.j_.  0
  res=. res , fassert 1 1 -: isnan  0j_.  0j_. rot  0    _.j_.
  res=. res , fassert 1 1 -: isnan  0j_.  0j_. rot  0j_.  0j_.
  res=. res , fassert 1 1 -: isnan  0j_.  0j_. rot  0j_. _.
  res=. res , fassert 1 1 -: isnan  0j_.  0j_. rot _.     0j_.
  res=. res , fassert 1 1 -: isnan  0j_.  0j_. rot _.    _.
  res=. res , fassert 1 1 -: isnan  0j_.  0j_. rot _.j_.  0
  res=. res , fassert 1 1 -: isnan  0j_. _.    rot  0    _.j_.
  res=. res , fassert 1 1 -: isnan  0j_. _.    rot  0j_.  0j_.
  res=. res , fassert 1 1 -: isnan  0j_. _.    rot  0j_. _.
  res=. res , fassert 1 1 -: isnan  0j_. _.    rot _.     0j_.
  res=. res , fassert 1 1 -: isnan  0j_. _.    rot _.    _.
  res=. res , fassert 1 1 -: isnan  0j_. _.    rot _.j_.  0
  res=. res , fassert 1 1 -: isnan _.     0j_. rot  0    _.j_.
  res=. res , fassert 1 1 -: isnan _.     0j_. rot  0j_.  0j_.
  res=. res , fassert 1 1 -: isnan _.     0j_. rot  0j_. _.
  res=. res , fassert 1 1 -: isnan _.     0j_. rot _.     0j_.
  res=. res , fassert 1 1 -: isnan _.     0j_. rot _.    _.
  res=. res , fassert 1 1 -: isnan _.     0j_. rot _.j_.  0
  res=. res , fassert 1 1 -: isnan _.    _.    rot  0    _.j_.
  res=. res , fassert 1 1 -: isnan _.    _.    rot  0j_.  0j_.
  res=. res , fassert 1 1 -: isnan _.    _.    rot  0j_. _.
  res=. res , fassert 1 1 -: isnan _.    _.    rot _.     0j_.
  res=. res , fassert 1 1 -: isnan _.    _.    rot _.    _.
  res=. res , fassert 1 1 -: isnan _.    _.    rot _.j_.  0
  res=. res , fassert 1 1 -: isnan _.j_.  0    rot  0    _.j_.
  res=. res , fassert 1 1 -: isnan _.j_.  0    rot  0j_.  0j_.
  res=. res , fassert 1 1 -: isnan _.j_.  0    rot  0j_. _.
  res=. res , fassert 1 1 -: isnan _.j_.  0    rot _.     0j_.
  res=. res , fassert 1 1 -: isnan _.j_.  0    rot _.    _.
  res=. res , fassert 1 1 -: isnan _.j_.  0    rot _.j_.  0
  res=. res , fassert 1 1 -: isnan _.j_. _.    rot _.j_. _.
  res=. res , fassert 1 1 -: isnan _.j_. _.    rot _.j_.  0j_.
  res=. res , fassert 1 1 -: isnan _.j_. _.    rot _.    _.j_.
  res=. res , fassert 1 1 -: isnan _.j_. _.    rot  0j_. _.j_.
  res=. res , fassert 1 1 -: isnan _.j_.  0j_. rot _.j_. _.
  res=. res , fassert 1 1 -: isnan _.j_.  0j_. rot _.j_.  0j_.
  res=. res , fassert 1 1 -: isnan _.j_.  0j_. rot _.    _.j_.
  res=. res , fassert 1 1 -: isnan _.j_.  0j_. rot  0j_. _.j_.
  res=. res , fassert 1 1 -: isnan _.    _.j_. rot _.j_. _.
  res=. res , fassert 1 1 -: isnan _.    _.j_. rot _.j_.  0j_.
  res=. res , fassert 1 1 -: isnan _.    _.j_. rot _.    _.j_.
  res=. res , fassert 1 1 -: isnan _.    _.j_. rot  0j_. _.j_.
  res=. res , fassert 1 1 -: isnan  0j_. _.j_. rot _.j_. _.
  res=. res , fassert 1 1 -: isnan  0j_. _.j_. rot _.j_.  0j_.
  res=. res , fassert 1 1 -: isnan  0j_. _.j_. rot _.    _.j_.
  res=. res , fassert 1 1 -: isnan  0j_. _.j_. rot  0j_. _.j_.
  res=. res , fassert 1 1 -: isnan _.j_. _.j_. rot _.j_. _.j_.
  res=. res , fassert 1 1 -: isnan         cs  rot  0     0j_.
  res=. res , fassert 1 1 -: isnan         cs  rot  0    _.
  res=. res , fassert 1 1 -: isnan         cs  rot  0j_.  0
  res=. res , fassert 1 1 -: isnan         cs  rot _.     0
  res=. res , fassert 1 1 -: isnan         cs  rot  0    _.j_.
  res=. res , fassert 1 1 -: isnan         cs  rot  0j_.  0j_.
  res=. res , fassert 1 1 -: isnan         cs  rot  0j_. _.
  res=. res , fassert 1 1 -: isnan         cs  rot _.     0j_.
  res=. res , fassert 1 1 -: isnan         cs  rot _.    _.
  res=. res , fassert 1 1 -: isnan         cs  rot _.j_.  0
  res=. res , fassert 1 1 -: isnan         cs  rot _.j_. _.
  res=. res , fassert 1 1 -: isnan         cs  rot _.j_.  0j_.
  res=. res , fassert 1 1 -: isnan         cs  rot _.    _.j_.
  res=. res , fassert 1 1 -: isnan         cs  rot  0j_. _.j_.
  res=. res , fassert 1 1 -: isnan         cs  rot _.j_. _.j_.
  res=. res , fassert 1 1 -: isnan         cs  rot _. qn1 q3
  res=. res , fassert 1 1 -: isnan         cs  rot _. qni q3
  res=. res , fassert 1 1 -: isnan         cs  rot _. qnj q3
  res=. res , fassert 1 1 -: isnan         cs  rot _. qnk q3
  res=. res , fassert 1 1 -: isnan (_. qn1 cs) rot        q3
  res=. res , fassert 1 1 -: isnan (_. qni cs) rot        q3
  res=. res , fassert 1 1 -: isnan (_. qnj cs) rot        q3
  res=. res , fassert 1 1 -: isnan (_. qnk cs) rot        q3
  res=. res , fassert 1 1 -: isnan  0     0j_. rot        q3
  res=. res , fassert 1 1 -: isnan  0    _.    rot        q3
  res=. res , fassert 1 1 -: isnan  0j_.  0    rot        q3
  res=. res , fassert 1 1 -: isnan _.     0    rot        q3
  res=. res , fassert 1 1 -: isnan  0    _.j_. rot        q3
  res=. res , fassert 1 1 -: isnan  0j_.  0j_. rot        q3
  res=. res , fassert 1 1 -: isnan  0j_. _.    rot        q3
  res=. res , fassert 1 1 -: isnan _.     0j_. rot        q3
  res=. res , fassert 1 1 -: isnan _.    _.    rot        q3
  res=. res , fassert 1 1 -: isnan _.j_.  0    rot        q3
  res=. res , fassert 1 1 -: isnan _.j_. _.    rot        q3
  res=. res , fassert 1 1 -: isnan _.j_.  0j_. rot        q3
  res=. res , fassert 1 1 -: isnan _.    _.j_. rot        q3
  res=. res , fassert 1 1 -: isnan  0j_. _.j_. rot        q3
  res=. res , fassert 1 1 -: isnan _.j_. _.j_. rot        q3
  res=. res , fassert (,.~ 0 1) -: isnan  cs               rot q3 ,: _. qn1 q4
  res=. res , fassert (,.~ 0 1) -: isnan  cs               rot q3 ,: _. qni q4
  res=. res , fassert (,.~ 0 1) -: isnan  cs               rot q3 ,: _. qnj q4
  res=. res , fassert (,.~ 0 1) -: isnan  cs               rot q3 ,: _. qnk q4
  res=. res , fassert (,.~ 0 1) -: isnan (cs ,: _. qn1 cs) rot q3
  res=. res , fassert (,.~ 0 1) -: isnan (cs ,: _. qni cs) rot q3
  res=. res , fassert (,.~ 0 1) -: isnan (cs ,: _. qnj cs) rot q3
  res=. res , fassert (,.~ 0 1) -: isnan (cs ,: _. qnk cs) rot q3
  NB. - input without edge cases
  res=. res , fassert 0 0 -: cs rot 0 0
  res=. res , fassert            q3  -:  1     0    rot q3
  res=. res , fassert (-         q3) -: _1     0    rot q3
  res=. res , fassert ( +        q3) -:  0j_1  0    rot q3
  res=. res , fassert (- &.(1&{) q3) -:  0     1    rot q3
  res=. res , fassert (- &.(0&{) q3) -:  0    _1    rot q3
  res=. res , fassert (-+        q3) -:  0     0j1  rot q3
  res=. res , fassert ( +        q3) -:  0     0j_1 rot q3
  res=. res , fassert ((,: -)  q3) -: (1 0 ,: _1 0) rot q3
  res=. res , fassert (q3 ,:   q4) -:  1 0          rot q3 ,: q4
  res=. res , fassert (q3 ,: - q4) -: (1 0 ,: _1 0) rot q3 ,: q4

  'rot' reportv (# ([ , -) +/) res
)
