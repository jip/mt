NB. pure_j_lapack.ijs
NB. Implement LAPACK routines on pure J

script_z_ '~system/packages/math/matutil.ijs'  NB. diag

NB. =========================================================
NB. Local verbs

prc=: [ (C. " 1) C.  NB. apply permutation x to both rows and columns of table y
p2P=: =/ (i. @ #)    NB. transform permutation vector y to permutation matrix

NB. =========================================================
NB. Interface verbs

NB. ---------------------------------------------------------
NB. gebalp
NB. General square matrix A balancing permutation:
NB.                   ( A00 A01 A02 )
NB. Ap = P*A*inv(P) = ( 0   A11 A12 )
NB.                   ( 0   0   A22 )
NB.
NB. Syntax:
NB.   'Ap ss p'=. gebalpi A
NB. where
NB.   A  - N-by-N matrix
NB.   Ap - N-by-N matrix, balanced by permutation
NB.   ss - 2-vector, corner start (left and up) and size
NB.        (width and height) of A11
NB.   p  - N-vector, standart permutation vector to
NB.        transform A to Aperm
NB.
NB. If:
NB.   'Ap ss p'=. gebalpi A
NB.   P=. p2P p
NB.   Pinv=. %. P
NB. then
NB.   Pinv -: |: P
NB.   Ap -: P mp A mp Pinv
NB.   Ap -: p prc A
NB.   A11 -: (,.~ ss) (] ;. 0) Ap

gebalp=: 3 : 0
  nz=. 0 & ~:
  nzd=. nz @ diag                             NB. non-zeros on diagonal
  mkcp=. < @ (, ` (, @ ]) @. =)               NB. make cycle permutation from x and y

  n=. # y
  f=. 0                                       NB. A11 up left corner start
  s=. n                                       NB. A11 size
  p=. i. n                                    NB. permutations accumulator
  'nzr nzc'=. (((+/ " 1 ,: (+/)) (- " 1) diag) @ nz) y  NB. count non-zeros in rows and cols of A without diagonal

NB.  nzr=. p { nzr                              NB. apply all permutations
  whilst. 1 do.                               NB. traverse rows from bottom to up
    zi=. f + (f ,: s) (i: & 0 ;. 0) nzr       NB. find rightmost zero's index into A11 non-zeros counter
    if. zi >: f + s do. break. end.           NB. no zeros left?
    nza=. (zi { p) (0:`[`(nz @ ({ " 1))) } y  NB. find nzr amendment after cols swapping (indirect zi because col zi may be swapped before)
    cp=. (<: f + s) mkcp zi                   NB. new cycle permutation
    p=. cp C. p                               NB. accumulate permutations
    nzr=. (cp C. nzr) - (p { nza)             NB. in nzr move zero found to end, apply all permutations to nza, then amend nzr by excluded row
    s=. <: s                                  NB. exclude last row
  end.

  nzc=. p { nzc                               NB. apply all permutations
  whilst. 1 do.                               NB. traverse columns from left to right
    zi=. f + (f ,: s) (i. & 0 ;. 0) nzc       NB. find leftmost zero's index into A11 non-zeros counter
    if. zi >: f + s do. break. end.           NB. no zeros left?
    nza=. (zi { p) (0:`[`(nz @ {)) } y        NB. find nzc amendment after rows swapping (indirect zi because row zi may be swapped before)
    cp=. f mkcp zi                            NB. new cycle permutation
    p=. cp C. p                               NB. accumulate permutations
    nzc=. (cp C. nzc) - (p { nza)             NB. in nzc move zero found to start, apply all permutations to nza, then amend nzc by excluded col
    f=. >: f                                  NB. exclude...
    s=. <: s                                  NB. ...leading column
  end.

  (p prc y) ; (f , s) ; p
)

NB. ---------------------------------------------------------
NB. gebals
NB. General square matrix A balancing by scaling:
NB.   A11s = D*A11*inv(D)
NB.
NB. Syntax:
NB.   'A11s s'=. gebalsi A11
NB. where
NB.   A11  - N-by-N matrix from gebalp's partition
NB.   A11s - N-by-N matrix, balanced by scaling
NB.   s    - N-vector, diagonal of scaling matrix D
NB.
NB. If:
NB.   'A11s s'=. gebalsi A11
NB.   D=. diagmat s
NB.   Dinv=. %. D
NB. then
NB.   Dinv -: |: D
NB.   A11s -: D mp A11 mp Dinv
NB.   A11s -: A11 (] * [ *"1 %) s
NB.
NB. Applications:
NB.   'A11s s'=. (,.~ ss) gebals ;. 0 Ap

gebals=: 3 : 0
  n=. # y
  scale=. n $ 1x
  RADIX=. 2x
  whilst. noconv do.
    noconv=. 0
    for_i. i. n do.
      'c r'=. i ([ (+/ @ (| @ ((0:`[`]) }))) ({ " 1 ,. {)) y  NB. 1-norm of i-th col and i-th row without diagonal element 
      if. (r ~: 0) *. (c ~: 0) do.
        m=. >. -: <: RADIX ^. r % c
        f=. RADIX ^ m
        if. (m ~: 0) *. (((r % f) + (c * f)) < (0.95 * (c + r))) do.
          scale=. (f * (i { scale)) i } scale
          noconv=. 1
          y=. ((i  {      y) % f) (i }    ) y
          y=. ((i ({ " 1) y) * f) (<a: ; i) } y
        end.
      end.
    end.
  end.
  y ; scale
)

NB. ---------------------------------------------------------
NB. gebal
NB. Balance square matrix A
NB. 'Ab ss p s'=. gebal A

gebal=: 3 : 0
  'Ap ss p'=. gebalp y
  'A11s s'=. (,.~ ss) gebals ;. 0 Ap
  IOS11=. ({. + (i. @ {:)) ss
  s=. s IOS11 } (# y) $ 1x
  Ab=. A11s (< ;~ IOS11) } Ap
  Ab ; ss ; p ; s
)

NB. ---------------------------------------------------------
NB. geequ
NB. Equilibrate matrix
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
  n=. # y
  d1=. n $ 1
  d2=. n $ 1
  whilst. (eps < vnormi <: dr) +. (eps < vnormi <: dc)  do.
NB.    'dr dc'=. (((1 & (I. @ (0 & = @ ])) }) @: %: &. >) @ ((>./ " 1 ; (>./ )) @: |)) y  NB. use this line instead following if row or column may contain zeros only
    'dr dc'=. (%: &. > @ ((>./ " 1 ; (>./ )) @: |)) y
    d1=. d1 % dr
    d2=. d2 % dc
    y=. d1 * y (* " 1) d2
  end.
  dr ; dc ; y
)


NB. =========================================================
Note 'gebal testing and timing'
   load '~addons/math/lapack/lapack.ijs'
   load '~addons/math/lapack/gebal.ijs'   NB. modified 'permute only' version
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
   Aperm -: p prc a1000f
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
