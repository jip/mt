NB. rot.ijs
NB. Rotations
NB.
NB. lartg   generates a plane rotation of a 2-vector
NB. lartv   applies a plane rotation to a 2-vector
NB. lartvt  applies a transposed plane rotation to a 2-vector
NB.
NB. Version: 1.0.0 2008-07-30
NB. Copyright: Igor Zhuravlov igor at uic.dvgu.ru
NB. License: Version 3 of the GNU GPL or any later version

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. various exponents for effective thresholds
FP_ESFMA=: 1 2 3 4 _1 _2 _3 _4 * >. 1r4 * (FP_BASE ^. FP_SFMIN % (FP_IGUNFL { (FP_EPS , 1)))

NB. effective underflow threshold FP_SFMN
NB.   = 2^_1020 (gradual underflow)
NB.   = 2^_968 (no gradual underflow)
NB. effective overflow threshold FP_SFMX
NB.   = 2^1020 (gradual underflow)
NB.   = 2^968 (no gradual underflow)
'FP_SFMN4 FP_SFMN2 FP_SFMN3 FP_SFMN FP_SFMX4 FP_SFMX2 FP_SFMX3 FP_SFMX'=: FP_BASE ^ FP_ESFMA

NB. miscellaneous
FP_SQRTEPS=: %: FP_EPS
FP_SQRTMX4=: %: FP_SFMX4
FP_MN2MX2=: (FP_SFMN2*(1-FP_EPS)) , FP_SFMX2
FP_MX1MN=: FP_SFMX , 1 , FP_SFMN
FP_MN1MX=: FP_SFMN , 1 , FP_SFMX
FP_MN3MN4MX4MX3=: (FP_SFMN3*(1-FP_EPS)) , (FP_SFMN4*(1-FP_EPS)) , FP_SFMX4 , FP_SFMX3
FP_MXMX21MN2MN=: FP_SFMX , FP_SFMX2 , 1 , FP_SFMN2 , FP_SFMN
FP_1MN2MNMN6MN8=: 1 , FP_SFMN2 , FP_SFMN , 0 0  NB. FP_SFMN4^6 = 0, FP_SFMN4^8 = 0
FP_MNMN21MX2MX=: |. FP_MXMX21MN2MN
FP_MN1324MX8231=: (FP_SFMN*(1-FP_EPS)) , (FP_SFMN3*(1-FP_EPS)) , (FP_SFMN2*(1-FP_EPS)) , (FP_SFMN4*(1-FP_EPS)) , FP_SQRTMX4 , FP_SFMX2 , FP_SFMX3 , FP_SFMX
FP_MX13241MN4231=: FP_SFMX , FP_SFMX3 , FP_SFMX2 , FP_SFMX4 , 1 , FP_SFMN4 , FP_SFMN2 , FP_SFMN3 , FP_SFMN
FP_MN13241MX4231=: FP_SFMN , FP_SFMN3 , FP_SFMN2 , FP_SFMN4 , 1 , FP_SFMX4 , FP_SFMX2 , FP_SFMX3 , FP_SFMX

abs1=: >./ @: | @ +.   NB. max(|Re(y)|,|Im(y)|)
sgnr=: sgn @ (9 & o.)  NB. sgn(Re(y))

NB. ---------------------------------------------------------
NB. lartgc1
NB. case 1: f≠0, g≠0, neither f nor g too big or small,
NB. minimal work

lartgc1=: 3 : 0
  'f2 g2'=. abssq 'f g'=. y
  d1=. (sgnr f) * %: f2 * fg2=. f2 + g2  NB. (% sgnr) = (* sgnr)
  (f2 , (f * + g) , (f * fg2)) % d1
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. lartg                                                   1
NB. Generates a plane rotation of 2-vector:
NB.   [  cs  sn  ]     [ f ]     [ r ]
NB.   [  __      ]  .  [   ]  =  [   ]   where cs^2 + |sn|^2 = 1
NB.   [ -sn  cs  ]     [ g ]     [ 0 ]
NB.
NB. Syntax:
NB.   'cs sn r'=. lartg (f,g)
NB. where
NB.   (f,g) - 2-vector to rotate
NB.   cs sn - representation of 2×2 rotation matrix
NB.   r     - representation of rotated 2-vector (r,0)
NB.
NB. If:
NB.   'cs sn r'=. lartg (f,g)
NB. then
NB.   (r,0) -: G1 mp (f,g)
NB.   (r,0) -: (f,g) mp G2
NB.   (0,r) -: G2 mp (g,f)
NB.   (0,r) -: (g,f) mp G1
NB. where
NB.   G1=. 2 2 $ cs , sn , (- + sn) , cs
NB.   G2=. |: G1
NB.
NB. Notes:
NB. - input NaN or ±∞ leads to inconsistent output or NaN
NB.   error
NB. - other input may lead to ±∞ in output
NB.
NB. References:
NB. [1] D. Bindel, J. Demmel, W. Kahan, O. Marques. (2001) On
NB.     Computing Givens rotations reliably and efficiently.
NB.     LAPACK Working Note 148, University of Tennessee,
NB.     UT-CS-00-449, January 31, 2001.
NB.     http://www.netlib.org/lapack/lawns/downloads/
NB. [2] Anderson, Edward. (2000) Discontinuous Plane
NB.     Rotations and the Symmetric Eigenvalue Problem.
NB.     LAPACK Working Note 150, University of Tennessee,
NB.     UT-CS-00-454, December 4, 2000.
NB.     http://www.netlib.org/lapack/lawns/downloads/

lartg=: (3 : 0) " 1

  'f g'=. y
  'scalef scaleg'=. scalefg=. abs1 y

  if. scaleg = 0 do.
    NB. g=0, f may be 0
    (sgnr f) , 0 , f * cs
    return.

  elseif. scalef = 0 do.
    NB. f=0, g≠0
    NB. optionally scale g by z^(±4) so that z^_2 ≤ ||g|| ≤ z^2
    ios=. FP_MN2MX2 I. scaleg
    gs=. g * ios { FP_MX1MN
    d1=. %: abssq gs
    cs=. 0
    sn=. + gs % d1
    r=. d1 * ios { FP_MN1MX

  elseif. (scalef <: FP_SFMX4) *. (scalef >: FP_SFMN4) *. (scaleg <: FP_SFMX4) do.
    NB. case 1: f≠0, g≠0, neither f nor g too big or small, minimal work
    lartgc1 y
    return.

  elseif. scaleg < FP_SQRTEPS * scalef do.
    NB. case 2: f≠0, g≠0, |f|^2 + |g|^2 rounds to |f|^2
    iog=. {: iofg=. FP_MN2MX2 I. scalefg
    dscalefg=. iofg { FP_MX1MN
    'fs gs'=. y * dscalefg
    cs=. sgnr f
    sn=. cs * (iog { FP_MN1MX) * ({. dscalefg) * (fs * + gs) % abssq fs
    r=. cs * f

  elseif. scalef < FP_SQRTEPS * scaleg do.
    NB. case 3: f≠0, g≠0, |f|^2 + |g|^2 rounds to |g|^2
    iofg=. FP_MN3MN4MX4MX3 I. scalefg               NB. io{f,g} = count{f,g}+2 = {_2,_1,0,1,2}+2 = {0,1,2,3,4}
    dscalefg=. iofg { FP_MXMX21MN2MN
    'f2 g2'=. abssq 'fs gs'=. y * dscalefg
    d1=. (sgnr f) * %: f2 * g2
    cs=. d1 %~ f2 * (-~/ iofg) { FP_1MN2MNMN6MN8    NB. (countf-countg)≤0 => form non-neg. io={0,1,2,3,4}, cs *=  (z^2)^(countf-countg)
    sn=. d1 %~ fs * + gs
    r=. d1 %~ fs * g2 * ({: iofg) { FP_MNMN21MX2MX  NB. r *= (z^2)^countg

  elseif. do.
    NB. case 4: f≠0, g≠0, scale f and g up or down and use formula from case 1
    iof=. FP_MN1324MX8231 I. scalef
    'cs sn r'=. lartgc1 y * iof { FP_MX13241MN4231
    r=. r * iof { FP_MN13241MX4231

  end.

  cs , sn , r
)

NB. ---------------------------------------------------------
NB. lartv                                                 1 1
NB. lartvt                                                1 1
NB. Applies a plane rotation to a 2-vector:
NB. - lartv:
NB.   [ f2 ]    [  cs  sn  ]     [ f ]
NB.   [    ] := [  __      ]  .  [   ]
NB.   [ g2 ]    [ -sn  cs  ]     [ g ]
NB. - lartvt:          __
NB.   [ f2 ]    [ cs  -sn  ]     [ f ]
NB.   [    ] := [          ]  .  [   ]
NB.   [ g2 ]    [ sn   cs  ]     [ g ]
NB.
NB. Syntax:
NB.   'f2 g2'=. (cs,sn) lartv (f,g)
NB.   'f2 g2'=. (cs,sn) lartvt (f,g)
NB. where
NB.   (f,g)   - 2-vector to rotate
NB.   (cs,sn) - representation of 2×2 rotation matrix
NB.   (f2,g2) - 2-vector, rotated version of (f,g)
NB.
NB. If:
NB.   'f2 g2'=. (cs,sn) lartv (f,g)
NB.   'f2h g2h'=. (cs,sn) lartvt (f,g)
NB. then
NB.   (f2,g2)   -: G1 mp (f,g)
NB.   (f2h,g2h) -: G2 mp (f,g)
NB. where
NB.   G1=. 2 2 $ cs , sn , (- + sn) , cs
NB.   G2=. |: G1
NB.
NB. Applications:
NB.   'vf2 vg2'=. |: (vcs ,. vsn) lartv (vf ,. vg)
NB.
NB. Notes:
NB. - input NaN or ±∞ leads to inconsistent output or NaN
NB.   error
NB. - other input may lead to ±∞ in output

lartv=:  ([: : ((+/ @: *) , (+/ @: (* (_1 1 & * @ |. @: +))~))) " 1 1
lartvt=: ([: : ((+/ @: (* (1 _1 & * @: +))~) , (+/ @: (* |.)~))) " 1 1

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB.*tlartg v test lartg

NB.*tlartg v test lartg

tlartg=: (3 : 0) " 1
  'cs sn r'=. lartg y
  G=. 2 2 $ cs , sn , (- + sn) , cs
  (r,0) -: G mp y
)

NB. ---------------------------------------------------------
NB.*tlartg2 v test lartg

tlartg2=: (3 : 0) " 1
  erm=. 4 : '(| x - y) % ((FP_EPS_mt_ * | x) >. (FP_SFMIN_mt_ * FP_PREC_mt_))'
  'cs sn r'=. lartg y
  G=. 2 2 $ cs , sn , (- + sn) , cs
  r erm {. G mp y
)

NB. ---------------------------------------------------------
NB.*testlartg v test lartg

testlartg=: 3 : 0
  *./ tlartg ((9^2^2) , 2) $, (,"0/~) , (j./~) 1e_308 1e_200 1e_100 1e_50 10 1e50 1e100 1e200 1e308
)

NB. =========================================================
Note 'rot testing and timing'
   a=. (0.1 * 1e5 2 $ 100) * (10 ^ ? 1e5 2 $ 70)
   ac=. j./ (0.1 * 2 1e5 2 $ 100) * (10 ^ ? 2 1e5 2 $ 70)
   b=. (2 ^ _1074 + ? 1e5 2 $ (1023 - _1074))
   bc=. j./ (2 ^ _1074 + ? 2 1e5 2 $ (1023 - _1074))
   ts=: 6!:2, 7!:2@]
   load '/home/maya/j602-user/projects/mt/rot.ijs'
   2 ts 'z=. lartg_mt_ a'
13.005 3.10898e7
   2 ts 'z=. lartg_mt_ ac'
13.1008 3.10899e7
   ts 'z=. lartg_mt_ b'
   tlartg_mt_ ((4^2^2) , 2) $, (,"0/~) , (j./~) _ __ 0 _.
   tlartg_mt_ 0 ,.~ ~. , (j./~) 2 ^ _1074 _1020 _765 _510 _255 0 255 510 765 1020 1023
   tlartg_mt_ 0 ,. ~. , (j./~) 2 ^ _1074 _1020 _765 _510 _255 0 255 510 765 1020 1023
   ts 'z=. lartg_mt_ ((13^2^2) , 2) $, (,"0/~) , (j./~) 2 ^ _1074 _1020 _765 _510 _255 _10 0 10 255 510 765 1020 1023'
)
