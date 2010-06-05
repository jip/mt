NB. hrd.ijs
NB. Reduce a general matrix to upper Hessenberg form
NB.
NB. gehrd   reduce a general matrix to upper Hessenberg form
NB.
NB. Version: 1.0.0 2009-06-01
NB. Copyright: Igor Zhuravlov, igor at uic.dvgu.ru
NB. License: GNU GPL

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gehrd                                                   1
NB. Reduce a general matrix A to upper Hessenberg form H by
NB. an unitary similarity transformation:
NB.   Q' * A * Q = H
NB.
NB. Syntax:
NB.   'HQ tau'=. gehrd A ; ss
NB. where
NB.   A   - N×N-matrix with isolated eigenvalues (see gebal)
NB.   ss  - 2-vector, corner start (left and up) and size
NB.         (width and height) of A11 (see gebalp)
NB.   HQ  - N×N-matrix, ???
NB.   tau - (N-1)-vector, ???
NB.
NB. If:
NB.   'HQ tau'=. gehrd A ; ss
NB.   >>>>>>>>>>>>>>>>>D=. diagmat s
NB.   Dinv=. %. D
NB. then
NB.   Dinv -: diagmat % s
NB.   As -: D mp Ap mp Dinv
NB.   As -: Ap (] * (% " 1)) s
NB.
NB. Notes:
NB.
NB. References:
NB. [1] 
NB. [2] 

gehrd=: (3 : 0) " 1

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
  erm=. 4 : '(| x - y) % ((FP_EPS * | x) >. (FP_SFMIN * FP_PREC))'
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
