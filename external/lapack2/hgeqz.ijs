require 'math/lapack2'

NB. Description:
NB.   Reduce a pair of square matrices (A,B) to generalized
NB.   upper Hessenberg form by an unitary (orthogonal)
NB.   transformations, where A is a general matrix and B is
NB.   upper triangular
NB.
NB. Syntax:
NB.   'S P alpha beta Q Z'=. (job ; compQ ; compZ) xhgeqz ilo ; ihi ; H ; T ; Q1 ; Z1
NB. where
NB.   job   - string, case-insensitive, in which the head
NB.           specifies what to compute:
NB.             'E' - to compute eigenvalues only
NB.             'S' - to compute eigenvalues and the Schur
NB.                   form
NB.   compQ - string, case-insensitive, in which the head
NB.           specifies what to do with Q:
NB.             'N' - to not compute Q
NB.             'I' - to compute Q
NB.             'V' - to compute Q1*Q
NB.   compZ - string, case-insensitive, in which the head
NB.           specifies what to do with Z:
NB.             'N' - to not compute Z
NB.             'I' - to compute Z
NB.             'V' - to compute Z1*Z
NB.   ilo   ∈ [1,max(1,ihi)], IO starting row and column,
NB.           1-based
NB.   ihi   ∈ [min(ilo,n),n], IO ending row and column,
NB.           1-based
NB.   H     - n×n-matrix, upper Hessenberg in rows and
NB.           columns ilo:ihi and upper triangular outside
NB.   T     - n×n-matrix, upper triangular
NB.   Q1    - n×n-matrix, unitary (orthogonal) if compQ='V',
NB.           or any noun otherwise
NB.   Z1    - n×n-matrix, unitary (orthogonal) if compZ='V',
NB.           or any noun otherwise
NB.   S     - n×n-matrix, upper triangular from the
NB.           generalized Schur factorization if job='S', or
NB.           matrix with diagonal matching that of S
NB.   P     - n×n-matrix, upper triangular from the
NB.           generalized Schur factorization if job='S', or
NB.           matrix with diagonal matching that of P
NB.   alpha - n-vector, generalized eigenvalues nominator,
NB.           also a diagonal of S
NB.   beta  - n-vector, generalized eigenvalues denominator,
NB.           also a diagonal of P
NB.   Q     - 0×0-matrix if compQ='N', or n×n-matrix, unitary
NB.           (orthogonal), left Schur vectors:
NB.             Q     for matrix pair (H,T) if compQ='I'
NB.             Q1*Q  for matrix pair (A,B) if compQ='V'
NB.   Z     - 0×0-matrix if compZ='N', or n×n-matrix, unitary
NB.           (orthogonal), right Schur vectors:
NB.             Z     for matrix pair (H,T) if compZ='I'
NB.             Z1*Z  for matrix pair (A,B) if compZ='V'
NB.   n     ≥ 0, the size of H, T, S, P, alpha and beta
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dhgeqz=: 4 : 0
  'job compQ compZ'=. x
  'ilo ihi H T Q1 Z1'=. y
  n=. # H
  assert.                          (ismatrix_jlapack2_ ,  issquare_jlapack2_ ,  (ilo , ihi)&(uhmat_jlapack2_ ishessenberg_jlapack2_)         ) H
  assert.                          (ismatrix_jlapack2_ ,  issquare_jlapack2_ ,  (-: utri_jlapack2_)                                  ,  n = #) T
  assert. ('nNiI' e.~ {. compQ) +. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isorthogonal_jlapack2_                               *. n = #) Q1
  assert. ('nNiI' e.~ {. compZ) +. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isorthogonal_jlapack2_                               *. n = #) Z1
  Q1=. n (0 0 $ 0.0)"_`(0.0 $~ 2 # [)`(|:@])@.(1 3 5 I. 'nNiIvV' i. {. compQ) Q1
  Z1=. n (0 0 $ 0.0)"_`(0.0 $~ 2 # [)`(|:@])@.(1 3 5 I. 'nNiIvV' i. {. compZ) Z1
  ldHT=. , 1 >. n
  ldQ=.  , 1 >. 0:^:('nN' e.~ {. compQ) n
  ldZ=.  , 1 >. 0:^:('nN' e.~ {. compZ) n
  lwork=. , 1 >. n  NB. minimal
  'H T alphar alphai beta Q1 Z1'=. (|:L:0) 7 9 11 12 13 14 16 { dhgeqz_jlapack2_ (, job) ; (, compQ) ; (, compZ) ; (, n) ; (, ilo) ; (, ihi) ; (|: H) ; ldHT ; (|: T) ; ldHT ; (n $ 0.0) ; (n $ 0.0) ; (n $ 0.0) ; Q1 ; ldQ ; Z1 ; ldZ ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to alphar, alphai and beta
  alpha=. alphar j. alphai
  H ; T ; alpha ; beta ; Q1 ; Z1
)

zhgeqz=: 4 : 0
  'job compQ compZ'=. x
  'ilo ihi H T Q1 Z1'=. y
  n=. # H
  assert.                          (ismatrix_jlapack2_ ,  issquare_jlapack2_ ,  (ilo , ihi)&(uhmat_jlapack2_ ishessenberg_jlapack2_)          ) H
  assert.                          (ismatrix_jlapack2_ ,  issquare_jlapack2_ ,  (-: utri_jlapack2_)                                   ,  n = #) T
  assert. ('nNiI' e.~ {. compQ) +. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isunitary_jlapack2_                                   *. n = #) Q1
  assert. ('nNiI' e.~ {. compZ) +. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isunitary_jlapack2_                                   *. n = #) Z1
  Q1=. n (0 0 $ 0j0)"_`(0j0 $~ 2 # [)`(|:@])@.(1 3 5 I. 'nNiIvV' i. {. compQ) Q1
  Z1=. n (0 0 $ 0j0)"_`(0j0 $~ 2 # [)`(|:@])@.(1 3 5 I. 'nNiIvV' i. {. compZ) Z1
  ldHT=. , 1 >. n
  ldQ=.  , 1 >. 0:^:('nN' e.~ {. compQ) n
  ldZ=.  , 1 >. 0:^:('nN' e.~ {. compZ) n
  lwork=. , 1 >. n  NB. minimal
  (|:L:0) 7 9 11 12 13 15 { zhgeqz_jlapack2_ (, job) ; (, compQ) ; (, compZ) ; (, n) ; (, ilo) ; (, ihi) ; (|: H) ; ldHT ; (|: T) ; ldHT ; (n $ 0.0) ; (n $ 0.0) ; Q1 ; ldQ ; Z1 ; ldZ ; (lwork $ 0j0) ; lwork ; (n $ 0.0) ; , _1
    NB. (|:) doesn't affect to alpha and beta
)
