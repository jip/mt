NB. rot.ijs
NB. Rotations
NB.
NB. lartg  generates a plane rotation of 2-vector
NB.
NB. Version: 1.0.0 2008-07-30
NB. Copyright: Igor Zhuravlov igor@uic.dvgu.ru
NB. License: Version 3 of the GNU GPL or any later version

coclass 'pjlap'

NB. =========================================================
NB. Local verbs

NB. =========================================================
NB. Interface verbs

NB. ---------------------------------------------------------
NB. lartg                                                   1
NB. Generates a plane rotation of 2-vector:
NB.   [  cs  sn  ]     [ f ]     [ r ]
NB.   [  __      ]  .  [   ]  =  [   ]   where cs**2 + |sn|**2 = 1
NB.   [ -sn  cs  ]     [ g ]     [ 0 ]
NB.
NB. Syntax:
NB.   'cs sn r'=. lartg (f,g)
NB. where
NB.   (f,g) - 2-vector to rotate
NB.   cs sn - representation of 2-by-2 rotation matrix
NB.   r     - representation of rotated 2-vector (r 0)
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

  abs1=. >./ @: | @ +.
  abssq=. +/ @: *: @ +.
  sgn=. (0 & <:) @ *

  'f g'=. y
  'scalef scaleg'=. abs1 y
  if. scaleg = 0 do.
    cs=. 1
    sn=. 0
    r=. f
  elseif. scalef = 0 do.
    cs=0
    gs=. g
    'scaleg scale
    r=. f
  elseif. scalef = 0 do.
  elseif. scalef = 0 do.
  elseif. scalef = 0 do.
  elseif. scalef = 0 do.
  elseif. scalef = 0 do.
  cs , sn , r
)

NB. =========================================================
NB.*tlartg v test lartg

tzlartg=: 3 : 0
match=. matchclean_jlapack_;;
smoutput 'cs sn r'=. zlartg zzero_jlapack_ + y
G=. 2 2 $ cs , sn , (- + sn) , cs
smoutput z=. (r,0) match clean G mp y
0 pick z
)

NB. =========================================================
NB. sequence of givens rotations
NB.
NB. Syntax:
NB.   'i vichg Gpacked'=. givenseq v
NB. where
NB.   v -: (v0 v1 ... v(N-1)), complex datatype
NB.   vchg -: (0 ... 0 , v'(i) , v(i+1) ... v(N-1))
NB.   Gpacked -: ((cs01 ... cs(i-1)(i)) ,: (sn01 ... sn(i-1)(i))

givenseq=: (3 : 0) " 1
  vcs=. vsn=. i.0
  ilast=. <: # y
  for_i. i. ilast do.
    'csi sni vi1chg'=. zlartg (i + 1 0) { y
NB. v2=. (i + 0 1) { y
NB. G=. 2 2 $ csi , sni , (+ - sni) , csi
NB. smoutput 2 4 $ 'i' ; 'v2' ; 'G' ; 'v2*G' ; i ; v2 ; G ; (v2 mp G)
    vcs=. vcs , csi
    vsn=. vsn , sni
    y=. vi1chg (i+1) } y
    if. 1 0 -: csi , sni do.
      ilast=. i
      break.
    end.
  end.
  ilast ; vi1chg ; vcs ,: vsn
)

NB. Gseq=. unpackGseq (vcs ,: vsn)
unpackGseq=: 3 : 0
)

NB. triangularize
gtzero3=: 3 : 0
)
