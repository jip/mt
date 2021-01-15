require 'math/lapack2'

NB. Description:
NB.   Compute for a pair of square nonsymmetric matrices the
NB.   generalized eigenvalues and, optionally, the left
NB.   and/or right generalized eigenvectors
NB.
NB. Syntax:
NB.   'alpha beta Vl Vr'=. (jobVl , jobVr) xggev A ;  B
NB.   'alpha beta Vl Vr'=. (jobVl , jobVr) xggev A ,: B
NB. where
NB.   jobVl - scalar, character, case-insensitive:
NB.             'N' - do not compute Vl
NB.             'V' - to compute Vl
NB.   jobVr - scalar, character, case-insensitive:
NB.             'N' - do not compute Vr
NB.             'V' - to computed Vr
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
NB. - the verbs below are loaded into the current locale
NB. - the computed eigenvectors are normalized to have
NB.   Euclidean norm equal to 1 and largest component real
NB. - each eigenvector is scaled so the largest component has
NB.     |Re(V(i))| + |Im(Vi)| = 1

dggev=: 4 : 0
  assert. *./ (2 = # x) , x e. 'nNvV'
  'jobVl jobVr'=. x
  'A B'=. y
  n=. # A
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_         ) A
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_ *. n = #) B
  select. 3!:0 A
    case. JCMPX do. A=. 9 o. A
    case. JFL   do.
    case.       do. A=. A + 0.0
  end.
  select. 3!:0 B
    case. JCMPX do. B=. 9 o. B
    case. JFL   do.
    case.       do. B=. B + 0.0
  end.
  Vl=. (0 0 [^:(jobVl e. 'nN') }. $ y) $ 0.0
  Vr=. (0 0 [^:(jobVr e. 'nN') }. $ y) $ 0.0
  ldAB=. , 1 >. n
  ldVl=. , 1 >. # Vl
  ldVr=. , 1 >. # Vr
  NB. lwork=. , 1 >. 8 * n  NB. minimal
  lwork=. , 1 >. n * 39  NB. optimal
  cdrc=. dggev_jlapack2_ (, jobVl) ; (, jobVr) ; (, n) ; (|: A) ; ldAB ; (|: B) ; ldAB ; (n $ 0.0) ; (n $ 0.0) ; (n $ 0.0) ; Vl ; ldVl ; Vr ; ldVr ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  'alphar alphai beta Vl Vr'=. (|: L: 0) 8 9 10 11 13 { cdrc  NB. (|:) doesn't affect to alphar, alphai and beta
  alpha=. alphar j. alphai
  if. # cx=. I. alphai ~: 0 do.
    if. jobVl e. 'vV' do. Vl=. cx cxpair_jlapack2_ Vl end.
    if. jobVr e. 'vV' do. Vr=. cx cxpair_jlapack2_ Vr end.
  end.
  alpha ; beta ; Vl ; Vr
)

zggev=: 4 : 0
  assert. *./ (2 = # x) , x e. 'nNvV'
  'jobVl jobVr'=. x
  'A B'=. y
  n=. # A
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_         ) A
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. n = #) B
  if. JCMPX ~: 3!:0 A do. A=. A + 0j0 end.
  if. JCMPX ~: 3!:0 B do. B=. B + 0j0 end.
  Vl=. (0 0 [^:(jobVl e. 'nN') }. $ y) $ 0j0
  Vr=. (0 0 [^:(jobVr e. 'nN') }. $ y) $ 0j0
  ldAB=. , 1 >. n
  ldVl=. , 1 >. # Vl
  ldVr=. , 1 >. # Vr
  NB. lwork=. , 1 >. 2 * n  NB. minimal
  lwork=. , 1 >. n * 33  NB. optimal
  cdrc=. zggev_jlapack2_ (, jobVl) ; (, jobVr) ; (, n) ; (|: A) ; ldAB ; (|: B) ; ldAB ; (n $ 0j0) ; (n $ 0j0) ; Vl ; ldVl ; Vr ; ldVr ; (lwork $ 0j0) ; lwork ; ((8 * n) $ 0.0) ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 8 9 10 12 { cdrc  NB. (|:) doesn't affect to alpha and beta
)
