NB. Balance a matrix or pair of matrices
NB.
NB. gebalxp    Isolate eigenvalues of a general square matrix
NB. gebalxs    Make the rows and columns of a general square
NB.            matrix as close in 1-norm as possible
NB. gebalx     Balance a general square matrix
NB.
NB. ggbalxp    Isolate eigenvalues in a pair of general
NB.            square matrices
NB. ggbalxs    Make the rows and columns in a pair of general
NB.            square matrices as close in 1-norm as possible
NB. ggbalx     Balance a pair of general square matrices
NB.
NB. testgebal  Test gebalx by general matrix given
NB. testggbal  Test ggbalx by pair of general square matrices
NB.            given
NB. testbal    Adv. to make verb to test gxbalx by
NB.            matrix(-ces) of generator and shape given
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Scaler constants

BALSCLFAC=:  2
BALFACTOR=:  0.95
BALESFMIN1=: FP_EMIN + FP_FLEN - 1  NB. _1022 + 53 - 1 = _970 , CHECKME: for no gradual underflow case only
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

BALPOWMIN=: BALSFMIN2 * BALSCLFAC ^ i.  1 + >. BALSCLFAC ^. BALSFMAX2

NB. Vector of values:
NB.   BALSFMAX2 / BALSCLFAC^(z-i)
NB. where
NB.   i = {0,1,...,z}
NB.   z = ⌈log_{BALSCLFAC}(BALSFMAX2)⌉

BALPOWMAX=: BALSFMAX2 % BALSCLFAC ^ i. _1 - >. BALSCLFAC ^. BALSFMAX2

NB. ---------------------------------------------------------
NB. gebalxp1d
NB.
NB. Description:
NB.   Adv. to make verb to traverse single direction
NB.   (row-wise or column-wise) within gebalxp process
NB.
NB. Syntax:
NB.   vapp=. ioz`getv`mkt`dhs gebalxp1d
NB. where
NB.   ioz  - dyad to scan vector, either (i.) or (i:), is
NB.          called as:
NB.            io=. (ioz & 0) vector
NB.   getv - dyad to extract vector from matrix, is either
NB.          ({) or ({"1), is called as:
NB.            vector=. iovector getv matrix
NB.   mkt  - dyad to prepare index for transposition, is
NB.          either (+ <:) or ([), is called as:
NB.            io=. h mkt s
NB.   dhs  - monad to reduce submatrix B11 by excluding row
NB.          and column which are intersecting in the element
NB.          with IO either (<0 0) or (<_1 _1), is either
NB.          (+&0 _1) or (+&1 _1), is called as:
NB.            hs=. dhs hs
NB.   vapp - dyad to traverse single direction, is
NB.          called as:
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
NB.   Adv. to make verb to traverse both directions (row-wise
NB.   and column-wise) within gebalxp process
NB.
NB. Syntax:
NB.   vapp=. getv0`getv1 gebalxp2d
NB. where
NB.   getv0 - dyad to extract vector (either row or column)
NB.           from matrix, is either ({) or ({"1), is called
NB.           as:
NB.             vector=. iovector getv matrix
NB.   getv1 - dyad to extract vector (either column or row)
NB.           from matrix, of direction opposite to getv0, is
NB.           either ({) or ({"1), is called as:
NB.             vector=. iovector getv matrix
NB.   vapp  - dyad to traverse both directions, is called as:
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

NB. ---------------------------------------------------------
NB. gebalsf
NB.
NB. Description:
NB.   Calculate scaling factor for gebals process
NB.
NB. Syntax:
NB.   f=. gebalsf a,b,c,d,e
NB. where
NB.   a - scalar to scale up
NB.   b = 1, scalar to scale up
NB.   c - scalar to scale down
NB.   d - scalar to control overflow
NB.   e - scalar to control underflow
NB.   f = BALSCLFAC^min(i,j,k)
NB.   i = ⌈log_{BALSCLFAC}(c/a)⌉, i.e. maximal integer
NB.       safisfying:
NB.         a*f < c/f
NB.   j = ⌈log_{BALSCLFAC}(BALSFMAX2/d)⌉, i.e. maximal
NB.       integer safisfying:
NB.         d*f < BALSFMAX2
NB.   k = ⌈log_{BALSCLFAC}(e/BALSFMIN2)⌉, i.e. maximal
NB.       integer safisfying:
NB.         e/f > BALSFMIN2
NB.
NB. Note:
NB. - conventional (closed) insertion point is calculated by:
NB.     c=. x I. y
NB.   and provides:
NB.     y <: c { x
NB. - alternative (open) insertion point is calculated by:
NB.     o=. x I. (1 + FP_PREC) * y
NB.   and provides:
NB.     y < o { x

gebalsf=: BALPOWMAX {~ (BALPOWMAX i. (1 { (BALSCLFAC ^ 1 1 _1) & ((*^:(({.<{:)@])^:_) (3&{.)))) <. (BALESFMAX1 - BALPOWMAX I. (1+FP_PREC)*(3{])) <. (BALPOWMIN I. (4{]))

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
NB.   'B p hs'=. gebalxp A
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
NB.        B=. p fp A
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
NB.   B -: p fp A
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

geballp=: ([ ((fp~ (0 & {::)) ; ]) (({`({"1) gebalxp2d) (((+/,:(+/"1)) (-"1) diag)@:(0&~:))))
gebalup=: ([ ((fp~ (0 & {::)) ; ]) ((({"1)`{ gebalxp2d) (((+/"1,:(+/)) (-"1) diag)@:(0&~:))))

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
NB.   1) extract B, p, h and s from y
NB.   2) initialize d (D's diagonal), bt (1st lIO behind
NB.      tail), rios (rIOS of B11 row/column within B
NB.      row/column)
NB.   3) do...:
NB.      3.1) mark scaling process as converged
NB.      3.2) traverse B11, for each i-th pair of row and
NB.           corresp. column:
NB.           3.2.1) extract them and laminate into rc, then
NB.                  replace diagonal elements by zeros
NB.           3.2.2) calculate r and c, the norm1t of row and
NB.                  column
NB.           3.2.3) calculate ra and ca, the magnitudes of
NB.                  largest element in row and column
NB.           3.2.4) calculate
NB.                    fup := BALSCLFAC^i
NB.                  where i is mininal integer safisfying:
NB.                    (c*fup) >= (r/(fup*BALSCLFAC))
NB.                  or
NB.                    max(fup,c*fup,ca*fup) >= BALSFMAX2
NB.                  or
NB.                    min(r/(fup*BALSCLFAC),ra/fup) <= BALSFMIN2
NB.           3.2.5) scale up c and ca, scale down r and ra:
NB.                    c  := c  * fup
NB.                    ca := ca * fup
NB.                    r  := r  / fup
NB.                    ra := ra / fup
NB.           3.2.6) calculate
NB.                    fdn := BALSCLFAC^i
NB.                  where i is mininal integer safisfying:
NB.                    (r*fdn) >= (c/(fdn*BALSCLFAC))
NB.                  or
NB.                    max(r*fdn,ra*fdn) >= BALSFMAX2
NB.                  or
NB.                    min(fdn,c/(fdn*BALSCLFAC),ca/fdn) <= BALSFMIN2
NB.           3.2.7) scale down c and scale up r:
NB.                    c := c / fdn
NB.                    r := r * fdn
NB.           3.2.8) calculate scale factor:
NB.                    f := fup / fdn
NB.           3.2.9) if r and c are changed by f
NB.                  considerably, and the following holds:
NB.                    f >= 1 OR d[i] >= 1 OR f*d[i] > BALSFMIN1
NB.                  and
NB.                    f <= 1 OR d[i] <= 1 OR d[i] < BALSFMAX1/BALSCLFAC
NB.                  then:
NB.                  3.2.9.1) scale d[i] by f up
NB.                  3.2.9.2) scale B[i,:] by f down
NB.                  3.2.9.3) scale B[:,i] by f up
NB.                  3.2.9.4) mark scaling process as not
NB.                           yet converged
NB.      ...while not converged
NB.   4) link B, p, hs and d to produce output
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
NB.   i.e. balance without eigenvalues isolating step:
NB.   'C p hs d'=. gebals (];(a:"_);(0,#)) A
NB.
NB. Notes:
NB. - models LAPACK's xGEBAL with 'S' option, with following
NB.   difference: ra and ca never get value from diagonal
NB.   element
NB.
NB. TODO:
NB. - consider parallel approach described in [2]
NB.
NB. References:
NB. [1] Kressner, D. 2004. Numerical methods and software for
NB.     general and structured eigenvalue problems. Ph.D.
NB.     thesis, TU Berlin, Institut für Mathematik, Berlin,
NB.     Germany.
NB. [2] P. R. Amestoy, I. S. Duff, D. Ruiz, and B. Uçar.
NB.     A parallel matrix scaling algorithm. May 6, 2008.
NB.     Technical report RAL-TR-2008-013
NB.     http://www.numerical.rl.ac.uk/reports/reports.html

gebals=: 3 : 0
  'B p hs'=. y
  'h s'=. hs
  d=. (# B) $ 1
  bt=. h + s
  rios=. ,. hs
  whilst. noconv do.
    noconv=. 0
    i=. <: h
    while. bt > i=. >: i do.
      rc=. rios (] ;. 0)"2 1 i ([ 0:`[`]}"1 { ,: {"1) B
      'r c'=. norm1t"1 rc
      if. r (*. & (0&~:)) c do.
        'ra ca'=. |@(liofmax { ])"1 rc
        sum=. r + c
        g=. r % BALSCLFAC
        fup=. gebalsf c,1,g,(1>.c>.ca),(ra<.g)
        c=. c * fup
        g=. c % BALSCLFAC
        r=. r % fup
        fdn=. gebalsf r,1,g,(r>.ra%fup),(fup<.g<.ca*fup)
        f=. fup % fdn
        c=. c % fdn
        r=. r * fdn
        if. (r+c) < BALFACTOR*sum do.
          di=. i { d
          if. f (*. & (<&1)) di do.
            if. BALSFMIN1 >: f*di do. continue. end.
          end.
          if. f (*. & (>&1)) di do.
            if. di >: BALSFMAX1%f do. continue. end.
          end.
          d=. (di*f) i } d
          B=. i  (%&f) upd    B
          B=. i ((*&f) upd)"1 B
          noconv=. 1
        end.
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
NB.   C11 -: (,.~ hs) (] ;. 0) C
NB. where
NB.   'C p hs d'=. gebalx A
NB.   B=. p fp A
NB.   P=. p2P p
NB.   Pinv=. %. P
NB.   D=. diagmat d
NB.   Dinv=. %. D
NB.
NB. Notes:
NB. - gebalu models LAPACK's xGEBAL with 'B' option

geball=: gebals @ geballp
gebalu=: gebals @ gebalup

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgebal
NB.
NB. Description:
NB.   Test:
NB.   - gebal (math/lapack)
NB.   - gebalx (math/mt)
NB.   by general (possibly sparse) matrix given
NB.
NB. Syntax:
NB.   testgebal A
NB. where
NB.   A - n×n-matrix
NB.
NB. TODO:
NB. - consider [1]
NB.
NB. References:
NB. [1] Michael H. Schneider, Stavros A. Zenios. A
NB.     Comparative Study of Algorithms for Matrix Balancing.
NB.     Operations Research. Vol. 38. No. 3. May-June 1990

testgebal=: 3 : 0
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'gebal'

  rcond=. (norm1 con (getrilu1p@getrflu1p)) y

  ('gebal_jlapack_' tmonad (]`]`(rcond"_)`(_."_)`(_."_))) y

  ('geball' tmonad (]`]`(rcond"_)`(_."_)`(_."_))) y
  ('gebalu' tmonad (]`]`(rcond"_)`(_."_)`(_."_))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testbal
NB.
NB. Description:
NB.   Adv. to make verb to test gebalx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testbal
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     (? @ $ 0:) testbal_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testbal_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testbal_mt_ 150 200

testbal=: 1 : 'EMPTY_mt_ [ (testgebal_mt_ @ (u spmat_mt_ 0.25) ^: (=/))'
