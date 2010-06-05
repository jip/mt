NB. bal.ijs
NB. Balance a matrix or pair of matrices
NB.
NB. gebalp  Isolate eigenvalues of a general square matrix
NB. gebals  Make the rows and columns of a general square
NB.         matrix as close in 1-norm as possible
NB. gebal   Balance a general square matrix
NB. ggbal   Balance a pair of general square matrices
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

pt=: [ (C. " 1) C.   NB. apply permutation x to both rows and columns of table y
p2P=: =/ (i. @ #)    NB. transform permutation vector y to permutation matrix

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gebalp
NB. Permute a general square matrix A by a similarity
NB. transformation to isolate eigenvalues in diagonals of A00
NB. and A22:
NB.                         ( A00 A01 A02 )
NB.   Ap = P * A * inv(P) = ( 0   A11 A12 )
NB.                         ( 0   0   A22 )
NB. where A00 and A22 are square upper triangular matrices.
NB.
NB. Syntax:
NB.   'Ap fs p'=. gebalp A
NB. where
NB.   A  - N×N-matrix
NB.   Ap - N×N-matrix with isolated eigenvalues, being A
NB.        with permuted rows and columns
NB.   fs - 2-vector, corner start (left and up) and size
NB.        (width and height) of A11
NB.   p  - N-vector, standard permutation vector,
NB.        presentation of permutation matrix P
NB.
NB. If:
NB.   'Ap fs p'=. gebalp A
NB.   P=. p2P p
NB.   Pinv=. %. P
NB. then (with appropriate comparison tolerance)
NB.   Pinv -: |: P
NB.   Ap -: P mp A mp Pinv
NB.   Ap -: p pt A
NB.   A11 -: (,.~ fs) (] ;. 0) Ap
NB.
NB. References:
NB. [1] Kressner, D. 2004. Numerical methods and software for
NB.     general and structured eigenvalue problems. Ph.D.
NB.     thesis, TU Berlin, Institut für Mathematik, Berlin,
NB.     Germany.
NB.
NB. Notes:
NB. - result is identical to LAPACK's dgebal/zgebal with 'P'
NB.   option
NB.
NB. TODO:
NB. - try to apply incremental permut. cp to y

gebalp=: 3 : 0
  nz=. 0 & ~:
  nzd=. nz @ diag                             NB. non-zeros on diagonal

  n=. # y
  f=. 0                                       NB. A11 upper left corner start
  s=. n                                       NB. A11 size
  p=. i. n                                    NB. permutations accumulator
  'nzr nzc'=. (((+/ " 1 ,: (+/)) (- " 1) diag) @ nz) y  NB. count non-zeros in rows and cols of A without diagonal

NB.  nzr=. p { nzr                              NB. apply all permutations
  whilst. 1 do.                               NB. traverse rows from bottom to up
    zi=. f + (f ,: s) (i: & 0 ;. 0) nzr       NB. find rightmost zero's index into A11 non-zeros counter
    if. zi >: f + s do. break. end.           NB. no zeros left?
    nza=. (zi { p) (0:`[`(nz @ ({ " 1))) } y  NB. find nzr amendment after cols swapping (indirect zi because col zi may be swapped before)
    cp=. (<: f + s) ii2cp zi                  NB. new cycle permutation
    p=. cp C. p                               NB. accumulate permutations
    nzr=. (cp C. nzr) - (p { nza)             NB. in nzr move zero found to end, apply all permutations to nza, then amend nzr by excluded row
    s=. <: s                                  NB. exclude last row
  end.

  nzc=. p { nzc                               NB. apply all permutations
  whilst. 1 do.                               NB. traverse columns from left to right
    zi=. f + (f ,: s) (i. & 0 ;. 0) nzc       NB. find leftmost zero's index into A11 non-zeros counter
    if. zi >: f + s do. break. end.           NB. no zeros left?
    nza=. (zi { p) (0:`[`(nz @ {)) } y        NB. find nzc amendment after rows swapping (indirect zi because row zi may be swapped before)
    cp=. f ii2cp zi                           NB. new cycle permutation
    p=. cp C. p                               NB. accumulate permutations
    nzc=. (cp C. nzc) - (p { nza)             NB. in nzc move zero found to start, apply all permutations to nza, then amend nzc by excluded col
    f=. >: f                                  NB. exclude...
    s=. <: s                                  NB. ...leading column
  end.

  (p pt y) ; (f , s) ; p
)

gebalp2=: 3 : 0
  x2b3=. < @ < @ < @ [

  n=. # y
  p=. i. n
  iL=. 0
  iH=. <: n
  whilst. swapped do.
    swapped=. 0
    i=. iL
    while. (i <: iH) *. (-. swapped) do.
      if. 0 = ((iL ,: (iH - iL)) ((+/ @: |) ;. 0) i (x2b3 { {) y) do.
        swapped=. 1
        cp=. i ii2cp iH
        y=. cp pt y
        p=. cp C. p
        iH=. <: iH
      end.
      i=. >: i
    end.
    j=. iL
    while. (j <: iH) *. (-. swapped) do.
      if. 0 = ((iL ,: (iH - iL)) ((+/ @: |) ;. 0) j (x2b3 { ({ " 1)) y) do.
        swapped=. 1
        cp=. j ii2cp iL
        y=. cp pt y
        p=. cp C. p
        iL=. >: iL
      end.
      j=. >: j
    end.
  end.

  y ; (iL , (>: iH - iL)) ; p
)

NB. ---------------------------------------------------------
NB. gebals
NB. Apply a diagonal similarity transformation to rows and
NB. columns of matrix A11 (see gebalp) to make the rows and
NB. columns as close in 1-norm as possible:
NB.   As = D * A * inv(D)
NB. where D is diagonal scaling matrix.
NB.
NB. Syntax:
NB.   'As fs p s'=. gebals Ap ; fs ; p
NB. where
NB.   Ap - N×N-matrix with isolated eigenvalues (see gebalp)
NB.   fs - 2-vector, corner start (left and up) and size
NB.        (width and height) of A11 (see gebalp)
NB.   p  - N-vector, standard permutation vector,
NB.        presentation of permutation matrix P (see gebalp)
NB.   As - N×N-matrix, scaled version of Ap
NB.   s  - N-vector, diagonal of scaling matrix D
NB.
NB. If:
NB.   'As fs p s'=. gebals Ap ; fs ; p
NB.   D=. diagmat s
NB.   Dinv=. %. D
NB. then (with appropriate comparison tolerance)
NB.   Dinv -: diagmat % s
NB.   As -: D mp Ap mp Dinv
NB.   As -: Ap (] * (% " 1)) s
NB.
NB. References:
NB. [1] Kressner, D. 2004. Numerical methods and software for
NB.     general and structured eigenvalue problems. Ph.D.
NB.     thesis, TU Berlin, Institut für Mathematik, Berlin,
NB.     Germany.
NB.
NB. Notes:
NB. - result isn't identical to LAPACK's dgebal/zgebal with 'S'
NB.   option, but is close to
NB.
NB. Applications:
NB.   'As fs p s'=. gebals (] ; (0 , #) ; (a: " _)) A

gebals=: 3 : 0
  'Ap fs p'=. y
  n=. # Ap
  s=. n $ 1x                                        NB. scaling matrix D diagonal
  RADIX=. x: FP_BASE                                NB. floating point base
  whilst. noconv do.
    noconv=. 0                                      NB. don't repeat scaling by default
    for_i. ({. + (i. @ {:)) fs do.                  NB. traverse A11
      'c r'=. (fs ,. 0 2) (+/ @: |) ;. 0 (i ([ ((0:`[`]) }) ({ " 1 ,. {)) Ap)  NB. 1-norm of i-th col and i-th row of A11 without diagonal element
      if. r (*. & (~: & 0)) c do.                   NB. protect against zero row or col
        sum=. c + r
        f=. 1x
        while. c < r % RADIX do.
          c=. c * RADIX
          r=. r % RADIX
          f=. f * RADIX
        end.
        while. c >: r * RADIX do.
          c=. c % RADIX
          r=. r * RADIX
          f=. f % RADIX
        end.
        if. (f ~: 1x) *. ((c + r) < (0.95 * sum)) do.  NB. does correction f meaningful?
          s=. (f * (i { s)) i } s                   NB. correct scale[i]
          Ap=. ((i  {      Ap) % f)        i  } Ap  NB. scale i-th row
          Ap=. ((i ({ " 1) Ap) * f) (<a: ; i) } Ap  NB. scale i-th col
          noconv=. 1                                NB. A11 traversing will be repeated
        end.
      end.
    end.
  end.
  Ap ; fs ; p ; s
)

NB. ---------------------------------------------------------
NB. gebal
NB. Balance a general square matrix A. This involves, first,
NB. isolating eigenvalues (see gebalp); and second, making
NB. the rows and columns of A11 as close in 1-norm as
NB. possible (see gebals):
NB.                     ( A00 A01 A02 )
NB.   Ap = P*A*inv(P) = ( 0   A11 A12 )
NB.                     ( 0   0   A22 )
NB.   Ab = D * Ap * inv(D)
NB. where
NB.   A00 and A22 - square upper triangular matrices
NB.   P           - permutation matrix
NB.   D           - diagonal scaling matrix
NB.
NB. Syntax:
NB.   'Ab fs p s'=. gebal A
NB. where
NB.   A  - N×N-matrix
NB.   Ab - N×N-matrix, balanced version of A
NB.   fs - 2-vector, corner start (left and up) and size
NB.        (width and height) of A11 (see gebalp)
NB.   p  - N-vector, standard permutation vector,
NB.        presentation of permutation matrix P (see gebalp)
NB.   s  - N-vector, diagonal of scaling matrix D (see
NB.        gebals)
NB.
NB. If:
NB.   'Ab fs p s'=. gebal A
NB.   P=. p2P p
NB.   Pinv=. %. P
NB.   D=. diagmat s
NB.   Dinv=. %. D
NB. then (with appropriate comparison tolerance)
NB.   Pinv -: |: P
NB.   Ap -: P mp A mp Pinv
NB.   Ap -: p pt A
NB.   Dinv -: diagmat % s
NB.   Ab -: D mp Ap mp Dinv
NB.   Ab -: Ap (] * (% " 1)) s
NB.   A11 -: (,.~ fs) (] ;. 0) Ap

gebal=: gebals @ gebalp

NB. ---------------------------------------------------------
NB. ggbal
NB. Balance a pair of general square matrices A and B.

ggbal=: [:

NB. =========================================================
Note 'bal testing and timing'
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
