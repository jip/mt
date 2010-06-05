NB. Balance a matrix or pair of matrices
NB.
NB. gebalxp  Isolate eigenvalues of a general square matrix
NB. gebalxs  Make the rows and columns of a general square
NB.          matrix as close in 1-norm as possible
NB. gebalx   Balance a general square matrix
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

BALSCLFAC=:  2
BALFACTOR=:  0.95
BALESFMIN1=: FP_EMIN + FP_FLEN - 1  NB. _1022 + 53 - 1 = _970
BALESFMAX1=: - BALESFMIN1
BALSFMIN1=: BALSCLFAC ^ BALESFMIN1
BALSFMAX1=: BALSCLFAC ^ BALESFMAX1
BALSFMIN2=: BALSFMIN1 * BALSCLFAC
BALSFMAX2=: BALSFMAX1 % BALSCLFAC

NB. Vector of values:
NB.   BALSFMIN2 * BALSCLFAC^i
NB. where
NB.   i = {0,1,...,z}
NB.   z = ⌈log_{BALSCLFAC}(BALSFMAX2)⌉

BALPOWNEG=: BALSFMIN2 * BALSCLFAC ^ i.  1 + >. BALSCLFAC ^. BALSFMAX2

NB. Vector of values:
NB.   BALSFMAX2 / BALSCLFAC^(z-i)
NB. where
NB.   i = {0,1,...,z}
NB.   z = ⌈log_{BALSCLFAC}(BALSFMAX2)⌉

BALPOWPOS=: BALSFMAX2 % BALSCLFAC ^ i. _1 - >. BALSCLFAC ^. BALSFMAX2

gebalsv=: BALPOWPOS {~ (BALPOWPOS i. (1 { (BALSCLFAC ^ 1 1 _1) & ((*^:(({.<{:)@])^:_) (3&{.)))) (<. dbg 'minL') (BALESFMAX1 - BALPOWPOS I. (1+FP_PREC)*(3{])) (<. dbg 'minR') (BALPOWNEG I. (4{]))

NB. ---------------------------------------------------------
NB. gebalxp1d
NB.
NB. Description:
NB.   Traverse single direction (row-wise or column-wise)
NB.   within gebalxp process
NB.
NB. Syntax:
NB.   vapp=. ioz`getv`mkt`dhs gebalxp1d
NB. where
NB.   ioz  - dyad to scan vector, either (i.) or (i:), is
NB.          evoked as:
NB.            io=. (ioz & 0) vector
NB.   getv - dyad to extract vector from matrix, is either
NB.          ({) or ({"1), is evoked as:
NB.            vector=. iovector getv matrix
NB.   mkt  - dyad to prepare index for transposition, is
NB.          either (+ <:) or ([), is evoked as:
NB.            io=. h mkt s
NB.   dhs  - monad to reduce submatrix B11 by excluding row
NB.          and column which are intersecting in the element
NB.          with IO either (<0 0) or (<_1 _1), is either
NB.          (+&0 _1) or (+&1 _1), is evoked as:
NB.            hs=. dhs hs
NB.   vapp - dyad to traverse single direction, is
NB.          evoked as:
NB.            'p hs'=. A vapp (p ; hs ; nz)
NB.   nz   - n-vector of non-negative integers, count of
NB.          non-zero elements in rows (columns) of A
NB.   A    - n×n-matrix
NB.   p    - n-vector, full permutation of A
NB.   hs   - 2-vector of integers (h,s) 'head' and 'size',
NB.          defines submatrix B11 position in B
NB.
NB. Algorithm:
NB.   In: A, p, hs, nz
NB.   Out: p, hs
NB.   1) extract verbs ioz, getv, mkt, dhs from u
NB.   2) extract p, hs, nz from y
NB.   3) extract h, s from hs
NB.   4) while there are zeros in nz fragment with rIOS
NB.      (,.hs), traverse direction:
NB.      4.1) find nz amendment (nza) after vectors swapping
NB.           (indirect zi because vector zi may be swapped
NB.           before)
NB.      4.2) prepare index for transposition:
NB.             io=. h mkt s
NB.      4.3) compose non-standard transposition nst
NB.      4.4) try to adjust p by nst
NB.           4.4.1) if failed (i.e. if zi=io), then leave p
NB.                  unchanged
NB.      4.5) adjust nz: try to move zero found in (4) to
NB.           edge, apply all permutations to nza, then amend
NB.           nz by excluded vector
NB.           4.5.1) if failed (i.e. if zi=io), then leave nz
NB.                  unchanged
NB.      4.6) adjust hs to exclude leading (tail) row and
NB.           column
NB.   5) assemble output

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

NB. ---------------------------------------------------------
NB. gebalxp2d
NB.
NB. Description:
NB.   Traverse both directions (row-wise and column-wise)
NB.   within gebalxp process
NB.
NB. Syntax:
NB.   vapp=. getv0`getv1 gebalxp2d
NB. where
NB.   getv0 - dyad to extract vector (either row or column)
NB.           from matrix, is either ({) or ({"1), is evoked
NB.           as:
NB.             vector=. iovector getv matrix
NB.   getv1 - dyad to extract vector (either column or row)
NB.           from matrix, of direction opposite to getv0, is
NB.           either ({) or ({"1), is evoked as:
NB.             vector=. iovector getv matrix
NB.   vapp  - dyad to traverse both directions, is evoked as:
NB.             'p hs'=. A vapp (nz0 ,: nz1)
NB.   nz0   - n-vector of non-negative integers, count of
NB.           non-zero elements in either rows or columns
NB.           excluding diagonal
NB.   nz1   - n-vector of non-negative integers, count of
NB.           non-zero elements in either columns or rows
NB.           excluding diagonal, opposite to nz0
NB.   A     - n×n-matrix
NB.   p     - n-vector, full permutation of A
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines submatrix B11 position in B
NB.
NB. Algorithm:
NB.   In: A, nz0, nz1
NB.   Out: p, hs
NB.   1) initialize p, hs so that B = B11 = A
NB.   2) extract nz0, which defines 1st traverse direction
NB.      (row-wise or column-wise), from y
NB.   3) use A, p, hs, nz0 to traverse A through the 1st
NB.      direction to accumulate permutations and to reduce
NB.      B11
NB.   4) extract nz1, which defines 2nd traverse direction
NB.      (row-wise or column-wise), opposite to nz0, from y,
NB.      then apply p to nz1
NB.   5) use A, p, hs, nz1 to traverse A through the 2nd
NB.      direction to accumulate permutations and to further
NB.      reduce B11
NB.   6) return p, hs

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
NB. gebalup
NB.
NB. Description:
NB.   Permute a general square matrix A by a similarity
NB.   transformation to isolate eigenvalues:
NB.     B = P * A * P
NB.
NB. Syntax:
NB.   'B p hs'=. geballp A
NB.   'B p hs'=. gebalup A
NB. where
NB.   A  - n×n-matrix
NB.   B  - n×n-matrix with isolated eigenvalues, being A
NB.        with permuted rows and columns
NB.   p  - n-vector, full permutation of A
NB.   hs - 2-vector of integers (h,s) 'head' and 'size',
NB.        defines submatrix B11 position in B
NB.
NB. Algorithm:
NB.   In: A
NB.   Out: B, p, hs
NB.   1) count non-zeros in rows (nzr) and columns (nzc) of A
NB.      without diagonal, then laminate vectors derived:
NB.        nz=. nzc ,: nzr  NB. geballp
NB.        nz=. nzr ,: nzc  NB. gebalup
NB.   2) traverse A's columns from right to left, then rows
NB.      from top to bottom (geballp), or rows from bottom to
NB.      top, then columns from left to right (gebalup), to
NB.      produce p, hs
NB.   3) apply full permutation p to A to produce B:
NB.        B=. p sp A
NB.   7) link B, p, hs to assemple output
NB.
NB. Storage layout:
NB.   geballp:                    gebalup:
NB.         ( Bl00          )         ( Bu00 B01 B02  )
NB.     B = ( B10  B11      )     B = (      B11 B12  )
NB.         ( B20  B21 Bl22 )         (          Bu22 )
NB. where
NB.   Bl00 Bl22 - square lower triangular matrices with
NB.               isolated eigenvalues in diagonal
NB.   Bu00 Bu22 - square upper triangular matrices with
NB.               isolated eigenvalues in diagonal
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   Pinv -: |: P
NB.   B -: P mp A mp Pinv          NB. apply p to rows and columns of A
NB.   B -: p sp A
NB.   B11 -: (,.~ hs) (] ;. 0) B
NB. where
NB.   'B hs p'=. gebalxp A
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
NB. - gebalup models LAPACK's xGEBAL with 'P' option

geballp=: ([ ((sp~ (0&{::)) ; ]) (({`({"1) gebalxp2d) (((+/,:(+/"1)) (-"1) diag)@:(0&~:))))
gebalup=: ([ ((sp~ (0&{::)) ; ]) ((({"1)`{ gebalxp2d) (((+/"1,:(+/)) (-"1) diag)@:(0&~:))))

NB. ---------------------------------------------------------
NB. gebals
NB.
NB. Description:
NB.   Apply a diagonal similarity transformation:
NB.     C = D^_1 * B * D
NB.   to make the 1-norms of each row of B11 and its
NB.   corresponding column as close as possible
NB.
NB. Syntax:
NB.   'C p hs d'=. gebals B ; p ; hs
NB. where
NB.   B  - n×n-matrix with isolated eigenvalues, the output
NB.        of gebalxp
NB.   p  - some not changing parameter, the output of gebalxp
NB.   hs - 2-vector of integers (h,s) 'head' and 'size',
NB.        defines submatrix B11 position in B, the output of
NB.        gebalxp
NB.   C  - n×n-matrix, scaled version of B
NB.   d  - n-vector, diagonal of scaling matrix D
NB.
NB. Algorithm:
NB.   In: B, p, hs
NB.   Out: C, d
NB.   1) extract B, p and hs from y
NB.   2) initialize d, bt (1st lIO behind tail), rios (rIOS
NB.      of B11 row/column within B row/column)
NB.   3) do...:
NB.      3.1) mark scaling process as converged
NB.      3.2) traverse B11, for each i-th pair of row and
NB.           corresp. column:
NB.           3.2.1) extract them and laminate into rc
NB.           3.2.2) calculate r and c, the row-wise norm1t
NB.                  of rc
NB.           3.2.3) calculate ra and ca, the magnitude of
NB.                  largest element in each row of rc
NB.           3.2.4) save sum of r and c into sum
NB.           3.2.5) calculate:
NB.                    fup := 
NB.      ...while not converged
NB.   4) link d to y to produce output
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   Dinv -: diagmat % d
NB.   C -: Dinv mp B mp D
NB.   C -: B (*"1 % ]) d
NB. where
NB.   'C p hs d'=. gebals B ; p ; hs
NB.   D=. diagmat d
NB.   Dinv=. %. D
NB.
NB. Application:
NB. - scale non-permuted matrix A (default p and hs),
NB.   i.e. balanse without eigenvalues isolating step:
NB.   'C p hs d'=. gebalus (];(a:"_);(0,#)) A
NB.
NB. Notes:
NB. - models LAPACK's xGEBAL with 'S' option, with following
NB.   difference: ra and ca never get value from diagonal
NB.   element
NB.
NB. References:
NB. [1] Kressner, D. 2004. Numerical methods and software for
NB.     general and structured eigenvalue problems. Ph.D.
NB.     thesis, TU Berlin, Institut für Mathematik, Berlin,
NB.     Germany.

gebals=: 3 : 0
  'B p hs'=. y
  'h s'=. hs
  d=. (# B) $ 1
  bt=. h + s
  rios=. ,. hs
  whilst. noconv do.
    noconv=. (] dbg 'noconv') 0
    i=. <: h
    while. bt > i=. >: i do.
      rc=. rios (] ;. 0)"2 1 i ([ 0:`[`]}"1 { ,: {"1) B
      'r c'=. (norm1t"1 dbg 'R C') rc
      'ra ca'=. (|@(liofmax { ])"1 dbg 'RA CA') rc
      sum=. r + c
      g=. r (% dbg 'new g') BALSCLFAC
      fup=. (gebalsv dbg 'fup=gebalsv()') c,1,g,(1>.c>.ca),(ra<.g)
      c=. c (* dbg 'new c') fup
      g=. c (% dbg 'new g') BALSCLFAC
      r=. r (% dbg 'new r') fup
      fdn=. (gebalsv dbg 'fdn=gebalsv()') r,1,g,(r>.ra%fup),(fup<.g<.ca*fup)
      f=. fup (% dbg 'final f') fdn
      c=. c (% dbg 'final c') fdn
      r=. r (* dbg 'final r') fdn
      if. (r+c) < BALFACTOR*sum do.
        smoutput 'Case 4: (C+R)<0.95*S'
        di=. i { d
        if. f (*. & (<&1)) di do.
          smoutput 'Case 4a: F<1'
          if. BALSFMIN1 >: f*di do. continue. end.
        end.
        if. f (*. & (>&1)) di do.
          smoutput 'Case 4b: F>1'
          if. di >: BALSFMAX1%f do. continue. end.
        end.
        smoutput 'Case 5'
        d=. i ((*&f) upd1) dbg 'new d' d       NB. CHECKME!
        noconv=. (] dbg 'noconv') 1
        B=. i ((%&f) upd1    ) dbg 'B after B[I,:]/F' B
        B=. i (((*&f) upd1)"1) dbg 'B after B[:,I]*F' B
      end.
    end.
  end.
  B ; p ; hs ; d
)

NB. ---------------------------------------------------------
NB. geball
NB. gebalu
NB.
NB. Description:
NB.   Balance a general square matrix A. This involves,
NB.   first, isolating eigenvalues (see gebalxp):
NB.     B = P * A * P
NB.   and second, making the rows and columns of B11 as close
NB.   in 1-norm as possible (see gebals):
NB.     C = D^_1 * B * D
NB.
NB. Syntax:
NB.   'C p hs d'=. gebalx A
NB. where
NB.   A  - n×n-matrix
NB.   C  - n×n-matrix, balanced version of A
NB.   p  - n-vector, full permutation of A
NB.   hs - 2-vector of integers (h,s) 'head' and 'size',
NB.        defines submatrix B11 position in B (see gebalxp)
NB.   d  - n-vector, diagonal of scaling matrix D (see
NB.        gebals)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   Pinv -: |: P
NB.   Dinv -: diagmat % d
NB.   B -: P mp A mp Pinv
NB.   C -: Dinv mp B mp D
NB.   C -: B (*"1 % ]) d
NB.   B11 -: (,.~ hs) (] ;. 0) B
NB. where
NB.   'C p hs d'=. gebalx A
NB.   B=. p sp A
NB.   P=. p2P p
NB.   Pinv=. %. P
NB.   D=. diagmat d
NB.   Dinv=. %. D
NB.
NB. Notes:
NB. - gebalu models LAPACK's xGEBAL with 'B' option

geball=: gebals @ geballp
gebalu=: gebals @ (gebalup dbg 'gebalup')

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
