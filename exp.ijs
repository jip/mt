NB. exp.ijs
NB. Matrix exponential
NB.
NB. geexp  matrix exponential of a general matrix
NB. diexp  matrix exponential of a diagonalizable matrix
NB. heexp  matrix exponential of a Hermitian matrix

coclass 'pjlap'

NB. =========================================================
NB. Local verbs

NB. =========================================================
NB. Interface verbs

NB. ---------------------------------------------------------
NB. geexp                                                   2
NB. Matrix exponential of a general matrix
NB.
NB. Syntax:
NB.   Aexp=. geexp A
NB. where
NB.
NB. If:
NB.   
NB. then
NB.   
NB.
NB. References:
NB. [1] N. J. Higham, The scaling and squaring method for the
NB.     matrix exponential revisited, SIAM J. Matrix Anal.
NB.     Appl. Vol. 26, No. 4, pp. 1179–1193
NB.
NB. Notes:
NB. - 
NB.
NB. TODO:
NB. - 

geexp=: (3 : 0) " 2

  NB. Padé approximant degrees m
  vm=. 3 5 7 9 13

  NB. coeffcients θ[m]
  theta=. 0.01495585217958292 0.2539398330063230 0.9504178996162932 2.097847961257068 5.371920351148152

  NB. coeffcients b[i] of degree 13 Padé approximant
  b=. 64764752532480000x 32382376266240000x 7771770303897600x 1187353796428800x 129060195264000x 10559470521600x 670442572800x 33522128640x 1323241920x 40840800x 960960x 16380x 182x 1x

  NB. preprocess A

  NB. - shift
  mu=. (tr % #) y
  y=. (- mu) shiftdiag y

  NB. --- NB. balance to reduce 1-norm of A (wait for pjlap project)
  NB. --- 'y s'=. 0 3 {:: gebals (] ; (0 , #) ; (a: " _)) y
>>>>>>>>>>>>
  NB. make report of matrices b[i]*M^i , where i=0..13
  R=. b * 14 powsm M

  for_m. vm do.
    if. (norm1 M) <: (m { theta) do.
      if. m <: 9 do.
      else.
      end.
      X=. ((^ mu) * d) * X (* " 1) di  NB. undo preprocessing
    end.
  end.

  s=. >. 2 ^. (norm1 M) % (13 { theta)
  M=. M % 2x ^ s

  NB. form [13/13] Padé approximant to expm(A)
  VU=. (2 | i. 14) +/. (b (* " 0 2) (14 powsm M))
  r13=. (gesvx @: (-/ ; +/)) VU

  X=. (mp~ ^: s) r13
  X=. ((^ mu) * d) * X (* " 1) di  NB. undo preprocessing
)

NB. ---------------------------------------------------------
NB. diexp                                                   2
NB. Matrix exponential of a diagonalizable matrix
NB.
NB. Syntax:
NB.   Aexp=. diexp (VR ; Λ)
NB. where
NB.
NB. If:
NB.   
NB. then
NB.   
NB.
NB. Notes:
NB. - 
NB.
NB. TODO:
NB. - 

diexp=: ((0 & {::) ([ mp ((* (+ @ |:))~ ^)) (1 & {::)) : [: " 2

NB. ---------------------------------------------------------
NB. heexp                                                   2
NB. Matrix exponential of a Hermitian (symmetric if real)
NB. matrix
NB.
NB. Syntax:
NB.   Aexp=. diexp (VR ; Λ ; VRinv)
NB. where
NB.
NB. If:
NB.   
NB. then
NB.   
NB.
NB. Notes:
NB. - 
NB.
NB. TODO:
NB. - 

heexp=: (0 & {::) mp ((^ @: (1 & {::)) * (2 & {::)) : [: " 2

NB. ---------------------------------------------------------
NB. exp

exp=: [:

NB. =========================================================
Note 'exp testing and timing'
   load '~addons/math/lapack/lapack.ijs'
   load '~addons/math/lapack/gebal.ijs'
   load '/home/jip/j602-user/projects/pure_j_lapack.ijs'

   ts=: 6!:2, 7!:2@]
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