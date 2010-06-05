NB. Balance a matrix or pair of matrices
NB.
NB. gebalxp  Isolate eigenvalues of a general square matrix
NB. gebalxs  Make the rows and columns of a general square
NB.          matrix as close in 1-norm as possible
NB. gebalx   Balance a general square matrix
NB. ggbalx   Balance a pair of general square matrices
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

gebalxp1d=: 1 : 0
:
  '`ioz getv mkt dhs'=. u
  'p hs nz'=. y
  'h s'=. hs
  while.
    zi=. h + (h ,: s) (ioz & 0 ;. 0) nz
    zi < h + s
  do.
    nza=. (zi { p) (0:`[`(0 ~: getv)) } x
    nst=. < (h mkt s) , zi
    p=. nst C. :: ] p
    nz=. (nst C. :: ] nz) - (p { nza)
    'h s'=. dhs (h , s)
  end.
  p ; (h , s)
)

gebalxp2d=: 1 : 0
:
  n=. # x
  'p hs'=. x i:`(({.u)`:6)`(+ <:)`(+&0 _1) gebalxp1d ((i. n) ; (0,n) ;       ({. y) )
  'p hs'=. x i.`(({:u)`:6)`[     `(+&1 _1) gebalxp1d (p      ; hs    ; (p C. ({: y)))
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. geballp
NB.
NB. Description:
NB.   Permute a general square matrix A by a similarity
NB.   transformation to isolate eigenvalues in diagonals of
NB.   square upper triangular matrices B00 and B22:
NB.         ( B00         )
NB.     B = ( B10 B11     ) = P * A * P^_1
NB.         ( B20 B21 B22 )
NB.
NB. Syntax:
NB.   'B hs p'=. geballp A
NB. where
NB.   A  - n×n-matrix
NB.   B  - n×n-matrix with isolated eigenvalues, being A
NB.        with permuted rows and columns
NB.   hs - 2-vector of integers (h,s) 'head' and 'size',
NB.        defines submatrix B11 position in B
NB.   p  - n-vector, full permutation of A
NB.
NB. Algorithm:
NB.   In: A
NB.   Out: B, hs, p
NB.   1) initialize hs, p
NB.   2) count non-zeros in rows (nzr) and columns (nzc) of A
NB.      without diagonal
NB.   3) apply all permutations to nzr (optional)
NB.   4) while there are zeros in nzr fragment (,.hs),
NB.      traverse rows from top to bottom:
NB.      4.1) find nzr amendment (nza) after cols swapping
NB.           (indirect zi because column zi may be swapped
NB.           before)
NB.      4.2) compose non-standard transposition nst
NB.      4.3) try to adjust p by nst
NB.      4.4) adjust nzr: move zero found in (4) to start,
NB.           apply all permutations to nza, then amend nzr
NB.           by excluded row
NB.      4.5) adjust h, s to exclude leading row
NB.   5) apply all permutations to nzc
NB.   6) while there are zeros in nzc fragment (,. hs),
NB.      traverse columns from right to left:
NB.      6.1) find nzc amendment (nza) after rows swapping
NB.           (indirect zi because row zi may be swapped
NB.           before)
NB.      6.2) compose non-standard transposition nst
NB.      6.3) try to adjust p by nst
NB.      6.4) adjust nzc: move zero found in (6) to end,
NB.           apply all permutations to nza, then amend nzc
NB.           by excluded column
NB.      6.5) adjust s to exclude last column
NB.   7) link the following pieces with p to produce output:
NB.      7.1) apply full permutation p to A to produce B
NB.      7.2) assemble hs from h and s
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   Pinv -: |: P
NB.   B -: P mp A mp Pinv
NB.   B -: p sp A
NB.   B11 -: (,.~ hs) (] ;. 0) B
NB. where
NB.   'B hs p'=. geballp A
NB.   P=. p2P p
NB.   Pinv=. %. P

geballp2=: ([ ((sp~ (0&{::)) ; ]) (({`({"1) gebalxp2d) (((+/,:(+/"1)) (-"1) diag)@:(0&~:))))

NB. ---------------------------------------------------------
NB. gebalup
NB.
NB. Description:
NB.   Permute a general square matrix A by a similarity
NB.   transformation to isolate eigenvalues in diagonals of
NB.   square upper triangular matrices B00 and B22:
NB.         ( B00 B01 B02 )
NB.     B = (     B11 B12 ) = P * A * P^_1
NB.         (         B22 )
NB.
NB. Syntax:
NB.   'B hs p'=. gebalup A
NB. where
NB.   A  - n×n-matrix
NB.   B  - n×n-matrix with isolated eigenvalues, being A
NB.        with permuted rows and columns
NB.   hs - 2-vector of integers (h,s) 'head' and 'size',
NB.        defines submatrix B11 position in B
NB.   p  - n-vector, full permutation of A
NB.
NB. Algorithm:
NB.   In: A
NB.   Out: B, hs, p
NB.   1) initialize hs, p
NB.   2) count non-zeros in rows (nzr) and columns (nzc) of A
NB.      without diagonal
NB.   3) apply all permutations to nzr (optional)
NB.   4) while there are zeros in nzr fragment (,.hs),
NB.      traverse rows from bottom to top:
NB.      4.1) find nzr amendment (nza) after cols swapping
NB.           (indirect zi because column zi may be swapped
NB.           before)
NB.      4.2) compose non-standard transposition nst
NB.      4.3) try to adjust p by nst
NB.      4.4) adjust nzr: move zero found in (4) to end,
NB.           apply all permutations to nza, then amend nzr
NB.           by excluded row
NB.      4.5) adjust s to exclude last row
NB.   5) apply all permutations to nzc
NB.   6) while there are zeros in nzc fragment (,.hs),
NB.      traverse columns from left to right:
NB.      6.1) find nzc amendment (nza) after rows swapping
NB.           (indirect zi because row zi may be swapped
NB.           before)
NB.      6.2) compose non-standard transposition nst
NB.      6.3) try to adjust p by nst
NB.      6.4) adjust nzc: move zero found in (6) to start,
NB.           apply all permutations to nza, then amend nzc
NB.           by excluded column
NB.      6.5) adjust h, s to exclude leading column
NB.   7) link the following pieces with p to produce output:
NB.      7.1) apply full permutation p to A to produce B
NB.      7.2) assemble hs from h and s
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   Pinv -: |: P
NB.   B -: P mp A mp Pinv
NB.   B -: p sp A
NB.   B11 -: (,.~ hs) (] ;. 0) B
NB. where
NB.   'B hs p'=. gebalup A
NB.   P=. p2P p
NB.   Pinv=. %. P
NB.
NB. References:
NB. [1] Kressner, D. 2004. Numerical methods and software for
NB.     general and structured eigenvalue problems. Ph.D.
NB.     thesis, TU Berlin, Institut für Mathematik, Berlin,
NB.     Germany.
NB.
NB. Notes:
NB. - models LAPACK's xGEBAL with 'P' option

gebalup2=: ([ ((sp~ (0&{::)) ; ]) ((({"1)`{ gebalxp2d) (((+/"1,:(+/)) (-"1) diag)@:(0&~:))))

geballp=: 3 : 0
  n=. # y
  h=. 0
  s=. n
  p=. i. n
  'nzr nzc'=. ((+/ " 1 ,: (+/)) (- " 1) diag) 0 ~: y

NB.  nzc=. p C. nzc
  while.
    zi=. h + (h ,: s) (i: & 0 ;. 0) nzc
    zi < h + s
  do.
    nza=. (zi { p) (0:`[`(0 ~: {)) } y
    nst=. < (<: h + s) , zi
    p=. nst C. :: ] p
    nzc=. (nst C. :: ] nzc) - (p { nza)
    s=. <: s
  end.

  nzr=. p C. nzr
  while.
    zi=. h + (h ,: s) (i. & 0 ;. 0) nzr
    zi < h + s
  do.
    nza=. (zi { p) (0:`[`(0 ~: ({ " 1))) } y
    nst=. < h , zi
    p=. nst C. :: ] p
    nzr=. (nst C. :: ] nzr) - (p { nza)
    h=. >: h
    s=. <: s
  end.

  (p sp y) ; p ; (h , s)
)

gebalup=: 3 : 0
  n=. # y
  h=. 0
  s=. n
  p=. i. n
  'nzr nzc'=. ((+/ " 1 ,: (+/)) (- " 1) diag) 0 ~: y

NB.  nzr=. p C. nzr
  while.
    zi=. h + (h ,: s) (i: & 0 ;. 0) nzr
    zi < h + s
  do.
    nza=. (zi { p) (0:`[`(0 ~: ({ " 1))) } y
    nst=. < (<: h + s) , zi
    p=. nst C. :: ] p
    nzr=. (nst C. :: ] nzr) - (p { nza)
    s=. <: s
  end.

  nzc=. p C. nzc
  while.
    zi=. h + (h ,: s) (i. & 0 ;. 0) nzc
    zi < h + s
  do.
    nza=. (zi { p) (0:`[`(0 ~: {)) } y
    nst=. < h , zi
    p=. nst C. :: ] p
    nzc=. (nst C. :: ] nzc) - (p { nza)
    h=. >: h
    s=. <: s
  end.

  (p sp y) ; p ; (h , s)
)

NB. ---------------------------------------------------------
NB. gebalus
NB.
NB. Description:
NB.   Apply a diagonal similarity transformation:
NB.     C = D^_1 * B * D
NB.   to make the 1-norms of each row of B11 and its
NB.   corresponding column as close as possible
NB.
NB. Syntax:
NB.   'C hs p s'=. gebalus B ; hs ; p
NB. where
NB.   B  - n×n-matrix with isolated eigenvalues, the output
NB.        of gebalup
NB.   hs - 2-vector of integers (h,s) 'head' and 'size',
NB.        defines submatrix B11 position in B, the output of
NB.        gebalup
NB.   p  - some not changing parameter, the output of gebalup
NB.   C  - n×n-matrix, scaled version of B
NB.   s  - n-vector, diagonal of scaling matrix D
NB.
NB. Algorithm:
NB.   In: B hs p
NB.   Out: C s
NB.   1) initialize s
NB.
NB.
NB. Assertions (with appropriate comparison tolerance):#############
NB.   Dinv -: diagmat % s
NB.   As -: D mp Ap mp Dinv
NB.   As -: Ap (] * (% " 1)) s
NB. where
NB.   'As fs p s'=. gebalus Ap ; fs ; p
NB.   D=. diagmat s
NB.   Dinv=. %. D
NB.
NB. Application:
NB.   'As fs p s'=. gebalus (] ; (0 , #) ; (a: " _)) A
NB.
NB. Notes:
NB. - models LAPACK's xGEBAL with 'S' option, but with
NB.   slightly different result
NB.
NB. References:
NB. [1] Kressner, D. 2004. Numerical methods and software for
NB.     general and structured eigenvalue problems. Ph.D.
NB.     thesis, TU Berlin, Institut für Mathematik, Berlin,
NB.     Germany.

gebalus=: 3 : 0
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
NB. gebalu
NB. Balance a general square matrix A. This involves, first,
NB. isolating eigenvalues (see gebalup); and second, making
NB. the rows and columns of A11 as close in 1-norm as
NB. possible (see gebalus):
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
NB.   'Ab fs p s'=. gebalu A
NB. where
NB.   A  - N×N-matrix
NB.   Ab - N×N-matrix, balanced version of A
NB.   fs - 2-vector, corner start (left and up) and size
NB.        (width and height) of A11 (see gebalup)
NB.   p  - N-vector, standard permutation vector,
NB.        presentation of permutation matrix P (see gebalup)
NB.   s  - N-vector, diagonal of scaling matrix D (see
NB.        gebalus)
NB.
NB. If:
NB.   'Ab fs p s'=. gebalu A
NB.   P=. p2P p
NB.   Pinv=. %. P
NB.   D=. diagmat s
NB.   Dinv=. %. D
NB. then (with appropriate comparison tolerance)
NB.   Pinv -: |: P
NB.   Ap -: P mp A mp Pinv
NB.   Ap -: p sp A
NB.   Dinv -: diagmat % s
NB.   Ab -: D mp Ap mp Dinv
NB.   Ab -: Ap (] * (% " 1)) s
NB.   A11 -: (,.~ fs) (] ;. 0) Ap
NB.
NB. Notes:
NB. - models LAPACK's xGEBAL with 'B' option, but with
NB.   slightly different result

gebalu=: gebalus @ gebalup

NB. ---------------------------------------------------------
NB. ggbalu
NB. Balance a pair of general square matrices A and B.

ggbalu=: [:

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
   Aperm -: p sp a1000f
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
