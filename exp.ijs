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
NB.   E=. geexp A
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
NB. - shiftdiag for r[m](A)

geexp=: (3 : 0) " 2

  NB. Padé approximant degrees m
  vm=. 3 5 7 9 13

  NB. coeffcients θ[m]
  theta=. 0.01495585217958292 0.2539398330063230 0.9504178996162932 2.097847961257068 5.371920351148152

  NB. coeffcients b[i] of degree 13 Padé approximant
  b=. 64764752532480000x 32382376266240000x 7771770303897600x 1187353796428800x 129060195264000x 10559470521600x 670442572800x 33522128640x 1323241920x 40840800x 960960x 16380x 182x 1x

  NB. preprocess A

  NB. - shift
  mu=. (trace % #) y
  y=. (- mu) shiftdiag y

  NB. - balance to reduce 1-norm of A
  'y scale'=. 0 3 {:: gebals (] ; (0 , #) ; (a: " _)) y

  NB. find m
  iom=. theta I. norm1 y
  if. iom < # vm do.
    m=. iom { vm
  else.
    m=. {: vm
  end.

  if. m<13 do.
    NB. b[i] coeffs for V (1st row) and U (2nd row)
    bc=. _2 |:@(]\) (>: m) {. b

    NB. A powers (0 2 4 ...)
    pA=. (+: i. <. -: m) gepow y

    NB. U=. A*Σ(b[i+1]*(A^i),i=1,3,..,m  )
    NB. V=.   Σ(b[i  ]*(A^i),i=0,2,..,m-1)
    'V U'=. bc +/@(* " 1) pA
    U=. A mp U

    NB. find r[m](A)
    r=. gesv (V-U) ; (V+U)
  else.
    s=. >. 2 ^. (norm1 y) % {: theta
    y=. y % 2 ^ s
    'V U'=.
    U=. A mp U
    r=. gesv (V-U) ; (V+U)
    r=. mp~ ^: s r
  end.
  E=. (^ mu) * r (] * (% " 1)) scale
)

NB. ---------------------------------------------------------
NB. diexp                                                   1
NB. Matrix exponential of a diagonalizable matrix
NB.
NB. Syntax:
NB.   E=. diexp (RV ; ev ; RVi)
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
NB. heexp                                                   1
NB. Matrix exponential of a Hermitian matrix
NB.
NB. Syntax:
NB.   E=. diexp (RV ; ev)
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
