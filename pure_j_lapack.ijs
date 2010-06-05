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
NB. gebalpi
NB. General square matrix A balancing permutation info:
NB.                      ( A00 A01 A02 )
NB. Aperm = P*A*inv(P) = ( 0   A11 A12 )
NB.                      ( 0   0   A22 )
NB.
NB. Syntax:
NB.   'ss p'=. gebalpi A
NB. where
NB.   A  - N-by-N matrix to balance by permutation
NB.   ss - 2-vector, corner start (left and up) and size
NB.        (width and height) of A11
NB.   p  - N-vector, standart permutation vector to
NB.        transform A to Aperm
NB.
NB. If:
NB.   'ss p'=. gebalpi A
NB.   Aperm=. gebalp A
NB.   P=. p2P p
NB.   Pinv=. %. P
NB. then
NB.   Pinv -: |: P
NB.   Aperm -: P mp A mp Pinv
NB.   Aperm -: p prc A
NB.   A11 -: (,.~ ss) (] ;. 0) Aperm

gebalpi=: 3 : 0
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

  (f , s) ; p
)

NB. ---------------------------------------------------------
NB. gebalp
NB. General square matrix A balancing permutation:
NB.                      ( A00 A01 A02 )
NB. Aperm = P*A*inv(P) = ( 0   A11 A12 )
NB.                      ( 0   0   A22 )
NB.
NB. Syntax:
NB.   Aperm=. gebalp A
NB. where
NB.   A and Aperm - N-by-N matrces

gebalp=: prc~ ((1 & {::) @ gebalpi)

NB. ---------------------------------------------------------
NB. gebalsi
NB. General square matrix A balancing scaling info:
NB.   Ascal = D*A*inv(D)
NB.
NB. Syntax:
NB.   s=. gebalsi A
NB. where
NB.   A - N-by-N matrix to balance by scaling
NB.   s - N-vector, diagonal of scaling matrix D
NB.
NB. If:
NB.   s=. gebalsi A
NB.   Ascal=. gebals A
NB.   D=. diagmat s
NB.   Dinv=. %. D
NB. then
NB.   Dinv -: |: D
NB.   Ascal -: D mp A mp Dinv
NB.   Ascal -: s * A (* " 1) (% s)
NB.
NB. Applications:
NB.   A11scal=. (,.~ ss) gebals ;. 0 Aperm

gebalsi=: 3 : 0
  n=. # y
  s=. n $ 1
  noconv=. 1
  whilst. noconv do.
    'r c'=. ((+/ " 1 ,: (+/)) (- " 1) (+: @ diag)) @ |
    
  end.
)

NB. ############################################################

NB. ---------------------------------------------------------
NB. gebals
NB. General square matrix A balancing scaling:
NB.   Ascal = D*A*inv(D)
NB.
NB. Syntax:
NB.   Ascal=. gebals A
NB. where
NB.   A and Ascal - N-by-N matrces

gebals=: * ((*/ %) @ gebalsi)

NB. ############################################################
NB. ---------------------------------------------------------
NB. 'ss p s'=. gebali mat
gebali=. 3 : 0
  'ss p'=. gebalpi y
  ss ; p ; ss gebalsi y <<<<<<<<<<<<<<<<<<<<<,
  (((0 & {::))~ gebalpi) y
  ((0 & {::)) @ gebalpi
)

NB. ---------------------------------------------------------
NB. Abal=. gebal A
NB. Balance square matrix A

gebal=: 3 : '?????? gebals @ gebalp'

NB. ---------------------------------------------------------
Note 'gebal timing'
   load '~addons/math/lapack/lapack.ijs'
   load '~addons/math/lapack/gebal.ijs'   NB. modified 'permute only' version
   ts=: 6!:2, 7!:2@]
   a1000f=. dzero_jlapack_ + ? 1000 1000 $ 10
   a1000c=. zzero_jlapack_ + j./ ? 2 1000 1000 $ 10
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

NB. ---------------------------------------------------------
Note 'gebal testing'
   ] a=. 7 7 $ 6 0 0 0 0 1 0 0 4 0 0.00025 0.0125 0.02 0.125 1 128 64 0 0 _2 16 0 16384 0 1 _400 256 _4000 _2 _256 0 0.0125 2 2 32 0 0 0 0 0 0 0 0 8 0 0.004 0.125 _0.2 3
   load '/home/jip/j602-user/projects/pure_j_lapack.ijs'
   load '/home/jip/j602-user/projects/lapack/lapack.ijs'
   load '/home/jip/j602-user/projects/lapack/gebal.ijs'
   a1000f=. dzero_jlapack_ + ? 1000 1000 $ 10
   (gebalp -: (2b1000 & gebal_jlapack_)) 0 ((3 33 333 666) } " 1) 0 (2 25 125 800 999) } a1000f
   +/ +/ (gebalp ~: (2b1000 & gebal_jlapack_)) 0 ((3 33 333 666) } " 1) 0 (2 25 125 800 999) } a1000f
   tgebal_jlapack_ (0 ((3 33 333 666) } " 1) 0 (2 25 125 800 999) } a1000f)
   'ss p'=. gebalpi 0 ((3 33 333 666) } " 1) 0 (2 25 125 800 999) } a1000f
   P=. p2P p
   Pinv=. |: P
   Pinv -: %. P
   Aperm=. gebalp 0 ((3 33 333 666) } " 1) 0 (2 25 125 800 999) } a1000f
   Aperm -: p prc 0 ((3 33 333 666) } " 1) 0 (2 25 125 800 999) } a1000f
   Aperm -: (P mp (0 ((3 33 333 666) } " 1) 0 (2 25 125 800 999) } a1000f)) mp Pinv
)
NB. ############################################################
