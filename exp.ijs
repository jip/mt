NB. Matrix exponential
NB.
NB. geexp  Matrix exponential of a general matrix
NB. diexp  Matrix exponential of a diagonalizable matrix
NB. heexp  Matrix exponential of a Hermitian matrix
NB.
NB. XREF: gepow sdiag
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. coeffcients b[i] of degree 13 Padé approximant
PA13_GEEXP=: 64764752532480000x 32382376266240000x 7771770303897600x 1187353796428800x 129060195264000x 10559470521600x 670442572800x 33522128640x 1323241920x 40840800x 960960x 16380x 182x 1x

NB. coeffcients θ[m]
THETA_GEEXP=: 24774315513739r1656496414665010 2287286674603605r9007199254740992 1070075424639547r1125899906842624  1180983412074661r562949953421312 1512061155730925r281474976710656

NB. ---------------------------------------------------------
NB. r=. m geexpm2r A
NB. find r_m(A) = [m/m] Padé approximant to A

geexpm2r=: 4 : 0
  vbyrs=. +/ @ (* " 3 1)        NB. multiply vector by report, then sum
  sdiag02=. sdiag " 0 2
  b0b1=. ((0 0 ; 1 0) & {) @ ]  NB. extract (b[0] , b[1]) from y

  NB. b[i] coeffs for V (1st row) and U (2nd row)
  bc=. _2 (|: @ (]\)) (>: x) {. PA13_GEEXP

  NB. U=. A*Σ(b[i+1]*(A^i),i=0,2,..,m-1)
  NB. V=.   Σ(b[i  ]*(A^i),i=0,2,..,m-1)
  if. x < 13 do.
    NB. A powers (2 [4 [6 [8]]]), shape: p×N×N
    pA=. (+: }. i. >. -: x) gepow y

    NB. - multiply each table in pA by corresp. atom b[i], output: report 2×p×N×N
    NB. - sum multiplied tables, output: report 2×N×N
    NB. - shift 1st table's diagonal by b[0], 2nd table's diagonal by b[1]
    'V U'=. pA (b0b1 sdiag02 (vbyrs (0 1 & }.))) bc
  else.
    NB. report of A powers (2 4 6)
    pA=. 2 4 6 gepow y

    NB. V=. (b[8]*(A^2)+b[10]*(A^4)+b[12]*(A^6)) * (A^6)
    NB. U=. (b[9]*(A^2)+b[11]*(A^4)+b[13]*(A^6)) * (A^6)
    NB. VU=. V ,: U
    VU=. pA (({: @ [) (mp " 2 2) (vbyrs (2 _3 & {.))) bc

    NB. V=.      V + b[6]*A6+b[4]*A4+b[2]*A2+b[0]*I
    NB. U=. A * (U + b[7]*A6+b[5]*A4+b[3]*A2+b[1]*I)
    'V U'=. VU + (pA (b0b1 sdiag02 (vbyrs ((0 1 ,: 2 3) & (] ;. 0)))) bc)
  end.
  U=. y mp U

  NB. find r_m(A) = [m/m] Padé approximant to A
  r=. V (+ %. -) U  NB. r=. gesv V (- ; +) U
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. geexp                                                   2
NB. Matrix exponential of a general matrix
NB.
NB. Syntax:
NB.   E=. geexp A
NB. where
NB.   A - N×N-matrix, a general matrix
NB.   E - N×N-matrix, matrix exponent exp(A)
NB.   N >= 0
NB.
NB. References:
NB. [1] N. J. Higham, The scaling and squaring method for the
NB.     matrix exponential revisited, SIAM J. Matrix Anal.
NB.     Appl. Vol. 26, No. 4, pp. 1179–1193
NB.
NB. TODO:
NB. - FIXME: geexp_mt_ 3 3 $ 5 0 1r10 0 1 0 0 0 4


geexp=: (3 : 0) " 2
  NB. Padé approximant degrees m
  vm=. 3 5 7 9 13

  NB. preprocess A

  NB. - shift
  mu=. (trace % #) y
  y=. (- mu) sdiag y
smoutput 'shifted y' ; y
  NB. - balance to reduce 1-norm
  'y scale'=. 0 3 { gebals (] ; (0 , #) ; (a: " _)) y
smoutput 'balanced y' ; y
  NB. find max m such that norm1(A) <= θ[m] , or let m=13
  m=. ((THETA_GEEXP I. norm1 y) <. <: # vm) { vm

  NB. find r_m(A) = [m/m] Padé approximant to A
  if. m < 13 do.
    r=. m geexpm2r y
  else.
    s=. >. 2 ^. (norm1 y) % {: THETA_GEEXP  NB. find a minimal integer such that norm1(A/2^s)<=θ[13]
    y=. y % 2 ^ s                           NB. scaling
    r=. m geexpm2r y
smoutput 's' ; s ; 'r' ; r ; 'norm1 y' ; (norm1 y)
    r=. mp~ ^: s r                          NB. undo scaling
  end.
  E=. (^ mu) * r (] * (% " 1)) scale        NB. undo preprocessing
)

NB. ---------------------------------------------------------
NB. diexp                                                   1
NB. Matrix exponential of a diagonalizable matrix
NB.
NB. Syntax:
NB.   E=. diexp (rv ; ev ; rvi)
NB. where
NB.   rv  - N×N-matrix, right eigenvectors of A
NB.   ev  - N-vector, eigenvalues of A
NB.   rvi - N×N-matrix, inversion of rv
NB.   E   - N×N-matrix, matrix exponent exp(A)
NB.   N  >= 0

diexp=: (0 & {::) mp ((^ @: (1 & {::)) * (2 & {::)) : [: " 1

NB. ---------------------------------------------------------
NB. heexp                                                   1
NB. Matrix exponential of a Hermitian matrix
NB.
NB. Syntax:
NB.   E=. heexp (rv ; ev)
NB. where
NB.   rv - N×N-matrix, right eigenvectors of A
NB.   ev - N-vector, eigenvalues of A
NB.   E  - N×N-matrix, matrix exponent exp(A)
NB.   N >= 0

heexp=: ((0 & {::) ([ mp ((* (+ @ |:))~ ^)) (1 & {::)) : [: " 1

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
   'fs p'=. gebalpi a1000f
   P=. p2P p
   Pinv=. |: P
   Pinv -: %. P
   Aperm=. gebalp a1000f
   Aperm -: p pp a1000f
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
