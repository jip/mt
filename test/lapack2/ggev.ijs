require 'math/lapack2'

NB. Description:
NB.   Compute for a pair of square nonsymmetric matrices the
NB.   generalized eigenvalues and, optionally, the left
NB.   and/or right generalized eigenvectors
NB.
NB. Syntax:
NB.   'alpha beta Vl Vr'=. (jobVl ; jobVr) xggev A ;  B
NB.   'alpha beta Vl Vr'=. (jobVl ; jobVr) xggev A ,: B
NB. where
NB.   jobVl - literal, case-insensitive, in which the head
NB.           specifies whether to compute Vl:
NB.             'N' - to not compute
NB.             'V' - to compute
NB.   jobVr - literal, case-insensitive, in which the head
NB.           specifies whether to compute Vr:
NB.             'N' - to not compute
NB.             'V' - to compute
NB.   A,B   - n×n-matrix, a matrix pair
NB.   alpha - n-vector, generalized eigenvalues nominator of
NB.           matrix pair (A,B)
NB.   beta  - n-vector, generalized eigenvalues denominator
NB.           of matrix pair (A,B)
NB.   Vl    - n×n-matrix, left eigenvectors of A if jobVl='V'
NB.           or 0×0-matrix otherwise
NB.   Vr    - n×n-matrix, right eigenvectors of A if
NB.           jobVr='V' or 0×0-matrix otherwise
NB.   n     ≥ 0, the size of A, B, Vl, Vr, alpha and beta
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale
NB. - eigenvectors computed are normalized to have
NB.   Euclidean norm equal to 1 and largest component real
NB. - each eigenvector is scaled so the largest component has
NB.     |Re(V(i))| + |Im(Vi)| = 1

dggev=: 4 : 0
  'jobVl jobVr'=. x
  'A B'=. y
  n=. # A
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_        ) A
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) B
  Vl=. (0 0 [^:('nN' e.~ {. jobVl) }. $ y) $ 0.0
  Vr=. (0 0 [^:('nN' e.~ {. jobVr) }. $ y) $ 0.0
  ldAB=. , 1 >. n
  NB. lwork=. , 1 >. 8 * n  NB. minimal
  lwork=. , 1 >. n * 39  NB. optimal
  'alphar alphai beta Vl Vr'=. (|: L: 0) 8 9 10 11 13 { dggev_jlapack2_ (, jobVl) ; (, jobVr) ; (, n) ; (|: A) ; ldAB ; (|: B) ; ldAB ; (n $ 0.0) ; (n $ 0.0) ; (n $ 0.0) ; Vl ; (, 1 >. # Vl) ; Vr ; (, 1 >. # Vr) ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to alphar, alphai and beta
  alpha=. alphar j. alphai
  if. # cx=. I. alphai ~: 0 do.
    if. 'vV' e.~ {. jobVl do. Vl=. cx cxpair_jlapack2_ Vl end.
    if. 'vV' e.~ {. jobVr do. Vr=. cx cxpair_jlapack2_ Vr end.
  end.
  alpha ; beta ; Vl ; Vr
)

zggev=: 4 : 0
  'jobVl jobVr'=. x
  'A B'=. y
  n=. # A
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_        ) A
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , n = #) B
  Vl=. (0 0 [^:('nN' e.~ {. jobVl) }. $ y) $ 0j0
  Vr=. (0 0 [^:('nN' e.~ {. jobVr) }. $ y) $ 0j0
  ldAB=. , 1 >. n
  NB. lwork=. , 1 >. 2 * n  NB. minimal
  lwork=. , 1 >. n * 33  NB. optimal
  (|: L: 0) 8 9 10 12 { zggev_jlapack2_ (, jobVl) ; (, jobVr) ; (, n) ; (|: A) ; ldAB ; (|: B) ; ldAB ; (n $ 0j0) ; (n $ 0j0) ; Vl ; (, 1 >. # Vl) ; Vr ; (, 1 >. # Vr) ; (lwork $ 0j0) ; lwork ; ((8 * n) $ 0.0) ; , _1
    NB. (|:) doesn't affect to alpha and beta
)
