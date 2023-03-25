require 'math/lapack2'

NB. Description:
NB.   Reduce a pair of square matrices (A,B) to generalized
NB.   upper Hessenberg form by an unitary (orthogonal)
NB.   transformations, where A is a general matrix and B is
NB.   upper triangular
NB.
NB. Syntax:
NB.   'H T Q Z'=. (compQ ; compZ) xgghrd ilo ; ihi ; A ; B ; Q1 ; Z1
NB. where
NB.   compQ - literal, case-insensitive, in which the head
NB.           specifies what to do with Q:
NB.             'N' - to not compute Q
NB.             'I' - to compute Q
NB.             'V' - to compute Q1*Q
NB.   compZ - literal, case-insensitive, in which the head
NB.           specifies what to do with Z:
NB.             'N' - to not compute Z
NB.             'I' - to compute Z
NB.             'V' - to compute Z1*Z
NB.   ilo   ∈ [1,max(1,ihi)], IO starting row and column,
NB.           1-based
NB.   ihi   ∈ [min(ilo,n),n], IO ending row and column,
NB.           1-based
NB.   A     - n×n-matrix, the upper triangular in rows and
NB.           columns outside ilo:ihi
NB.   B     - n×n-matrix, the upper triangular
NB.   Q1    - n×n-matrix, the unitary (orthogonal) if
NB.           compQ='V' or any noun otherwise
NB.   Z1    - n×n-matrix, the unitary (orthogonal) if
NB.           compZ='V' or any noun otherwise
NB.   H     - n×n-matrix, the upper Hessenberg in rows and
NB.           columns ilo:ihi and upper triangular outside
NB.   T     - n×n-matrix, the upper triangular:
NB.             T = Q^H B Z
NB.   Q     - 0×0-matrix if compQ='N', or n×n-matrix, the
NB.           unitary (orthogonal):
NB.             Q    if compQ='I'
NB.             Q1*Q if compQ='V'
NB.   Z     - 0×0-matrix if compZ='N', or n×n-matrix, the
NB.           unitary (orthogonal):
NB.             Z    if compZ='I'
NB.             Z1*Z if compZ='V'
NB.   n     ≥ 0, the size of A, B, H and T
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dgghrd=: 4 : 0
  'compQ compZ'=. x
  'ilo ihi A B Q1 Z1'=. y
  n=. # A
  assert.                          (ismatrix_jlapack2_ ,  issquare_jlapack2_                                   ) A
  assert.                          (ismatrix_jlapack2_ ,  issquare_jlapack2_ ,  (-: utri_jlapack2_)    ,  n = #) B
  assert. ('nNiI' e.~ {. compQ) +. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isorthogonal_jlapack2_ *. n = #) Q1
  assert. ('nNiI' e.~ {. compZ) +. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isorthogonal_jlapack2_ *. n = #) Z1
  Q1=. n (0 0 $ 0.0)"_`(0.0 $~ 2 # [)`(|:@])@.(1 3 5 I. 'nNiIvV' i. {. compQ) Q1
  Z1=. n (0 0 $ 0.0)"_`(0.0 $~ 2 # [)`(|:@])@.(1 3 5 I. 'nNiIvV' i. {. compZ) Z1
  ldAB=. , 1 >. n
  ldQ=.  , 1 >. 0:^:('nN' e.~ {. compQ) n
  ldZ=.  , 1 >. 0:^:('nN' e.~ {. compZ) n
  cdrc=. dgghrd_jlapack2_ (, compQ) ; (, compZ) ; (, n) ; (, ilo) ; (, ihi) ; (|: A) ; ldAB ; (|: B) ; ldAB ; Q1 ; ldQ ; Z1 ; ldZ ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 6 8 10 12 { cdrc
)

zgghrd=: 4 : 0
  'compQ compZ'=. x
  'ilo ihi A B Q1 Z1'=. y
  n=. # A
  assert.                          (ismatrix_jlapack2_ ,  issquare_jlapack2_                                ) A
  assert.                          (ismatrix_jlapack2_ ,  issquare_jlapack2_ ,  (-: utri_jlapack2_) ,  n = #) B
  assert. ('nNiI' e.~ {. compQ) +. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isunitary_jlapack2_ *. n = #) Q1
  assert. ('nNiI' e.~ {. compZ) +. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isunitary_jlapack2_ *. n = #) Z1
  Q1=. n (0 0 $ 0j0)"_`(0j0 $~ 2 # [)`(|:@])@.(1 3 5 I. 'nNiIvV' i. {. compQ) Q1
  Z1=. n (0 0 $ 0j0)"_`(0j0 $~ 2 # [)`(|:@])@.(1 3 5 I. 'nNiIvV' i. {. compZ) Z1
  ldAB=. , 1 >. n
  ldQ=.  , 1 >. 0:^:('nN' e.~ {. compQ) n
  ldZ=.  , 1 >. 0:^:('nN' e.~ {. compZ) n
  cdrc=. zgghrd_jlapack2_ (, compQ) ; (, compZ) ; (, n) ; (, ilo) ; (, ihi) ; (|: A) ; ldAB ; (|: B) ; ldAB ; Q1 ; ldQ ; Z1 ; ldZ ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 6 8 10 12 { cdrc
)
