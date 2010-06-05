NB. equ.ijs
NB. Equilibrate a matrix
NB.
NB. geequ  equilibrate a matrix
NB. poequ  equilibrate a Hermitian positive definite matrix

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. geequ
NB. Equilibrate a matrix
NB.
NB. References:
NB. [1] "A parallel matrix scaling algorithm"
NB.     Patrick R. Amestoy, Iain S. Duff, Daniel Ruiz, and Bora Uçar
NB.     Technical report RAL-TR-2008-013, May 6, 2008
NB.
NB. TODO:
NB. - verify

geequ=: 3 : 0
  eps=. 1e_1
  'm n'=. $ y
  d1=. m $ 1
  d2=. n $ 1
  while.
NB.    'vr vc'=. (1 & (I. @ (0 & = @ ])) } &. > @ ((>./ " 1 ; (>./ )) @: |)) y  NB. use this line instead following if row or col may consists of zeros only
   'vr vc'=. ((>./ " 1 ; (>./ )) @: |) y
   (eps < normi <: vr) +. (eps < normi <: vc)
  do.
    dr=. %: vr              NB. Dr(k)   := diag(sqrt(r[i](k)))
    dc=. %: vc              NB. Dc(k)   := diag(sqrt(c[i](k)))
    y=. dr %~ y (% " 1) dc  NB. A(k+1)  := inv(Dr) * A(k) * inv(Dc)
    d1=. d1 % dr            NB. D1(k+1) := D1(k) * inv(Dr(k))
    d2=. d2 % dc            NB. D2(k+1) := D2(k) * inv(Dc(k))
  end.
  y ; d1 ; d2
)

NB. ---------------------------------------------------------
NB. poequ
NB. Equilibrate a Hermitian positive definite matrix

poequ=: [:

NB. =========================================================
Note 'equ testing and timing'
   load '~addons/math/lapack/lapack.ijs'
   load '~addons/math/lapack/geequ.ijs'
   load '/home/jip/j602-user/projects/mtack/mtack.ijs'

   ts=: 6!:2, 7!:2@]
   A=. (* (10 & ^))/ ? 2 10 10 $ 10
   'Aeq d1 d2'=. geequ A
   Aeq -: d1 * A (* " 1) d2
1
   A -: d1 %~ Aeq (% " 1) d2
1







   ] a=. 7 7 $ 6 0 0 0 0 1 0 0 4 0 0.00025 0.0125 0.02 0.125 1 128 64 0 0 _2 16 0 16384 0 1 _400 256 _4000 _2 _256 0 0.0125 2 2 32 0 0 0 0 0 0 0 0 8 0 0.004 0.125 _0.2 3
   a1000f=. 0 ((3 33 333 666) } " 1) 0 (2 25 125 800 999) } dzero_jlapack_ + ? 1000 1000 $ 10
   a1000c=. 0 ((3 33 333 666) } " 1) 0 (2 25 125 800 999) } zzero_jlapack_ + j./ ? 2 1000 1000 $ 10

   (gebalp -: (2b1000 & gebal_jlapack_)) a1000f
   +/ +/ (gebalp ~: (2b1000 & gebal_jlapack_)) a1000f
   tgebal_jlapack_ a1000f
   'ss p'=. gebalpi a1000f
   P=. p2P p
   Pinv=. |: P
   Pinv -: %. P
   Aperm=. gebalp a1000f
   Aperm -: p pt a1000f
   Aperm -: (P mp a1000f) mp Pinv

   5 ts 'gebal_jlapack_ a1000f'
0.239994 2.51897e7
   5 ts 'gebalp a1000f'
0.151756 1.68269e7
   (gebalp a1000f) -: (2b1000 gebal_jlapack_ a1000f)
1
   5 ts 'gebal_jlapack_ a1000c'
0.570236 4.19459e7
   5 ts 'gebalp a1000c'
0.245782 3.36287e7
   (gebalp a1000c) -: (2b1000 gebal_jlapack_ a1000c)
1

   5 ts 'gebal_jlapack_ a1000f'
0.289798 2.51897e7
   5 ts 'gebal_jlapack_ a1000f'
0.289943 2.51897e7
   5 ts 'gebal a1000f'
0.350457 2.9392e7
   5 ts 'gebal a1000f'
0.357612 2.9392e7

)
