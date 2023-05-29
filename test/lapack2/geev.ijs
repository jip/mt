require 'math/lapack2'

NB. Description:
NB.   Compute the eigenvalues and, optionally, the left
NB.   and/or right eigenvectors for square nonsymmetric
NB.   matrix
NB.
NB. Syntax:
NB.   'w Vl Vr'=. (jobVl ; jobVr) xgeev A
NB. where
NB.   jobVl - literal, case-insensitive, in which the head
NB.           specifies whether to compute Vl:
NB.             'N' - to not compute
NB.             'V' - to compute
NB.   jobVr - literal, case-insensitive, in which the head
NB.           specifies whether to compute Vr:
NB.             'N' - to not compute
NB.             'V' - to compute
NB.   A     - n×n-matrix
NB.   w     - n-vector, eigenvalues of A
NB.   Vl    - n×n-matrix, left eigenvectors of A or
NB.           0×0-matrix
NB.   Vr    - n×n-matrix, right eigenvectors of A or
NB.           0×0-matrix
NB.   n     ≥ 0, the size of A, Vl, Vr and w
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale
NB. - eigenvectors computed are normalized to have
NB.   Euclidean norm equal to 1 and largest component real

dgeev=: 4 : 0
  'jobVl jobVr'=. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  Vl=. (0 0 [^:('nN' e.~ {. jobVl) $ y) $ 0.0
  Vr=. (0 0 [^:('nN' e.~ {. jobVr) $ y) $ 0.0
  NB. lwork=. , 1 >. n * 4 [^:(x +./@e. 'vV') 3  NB. minimal
  lwork=. , 1 >. n * 130 [^:(x +./@e. 'vV') 34  NB. optimal
  'wr wi Vl Vr'=. (|: L: 0) 6 7 8 10 { dgeev_jlapack2_ (, jobVl) ; (, jobVr) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 0.0) ; (n $ 0.0) ; Vl ; (, 1 >. # Vl) ; Vr ; (, 1 >. # Vr) ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to wr and wi
  w=. wr j. wi
  if. # cx=. I. wi ~: 0 do.
    if. 'vV' e.~ {. jobVl do. Vl=. cx cxpair_jlapack2_ Vl end.
    if. 'vV' e.~ {. jobVr do. Vr=. cx cxpair_jlapack2_ Vr end.
  end.
  w ; Vl ; Vr
)

zgeev=: 4 : 0
  'jobVl jobVr'=. x
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  n=. # y
  Vl=. (0 0 [^:('nN' e.~ {. jobVl) $ y) $ 0j0
  Vr=. (0 0 [^:('nN' e.~ {. jobVr) $ y) $ 0j0
  NB. , lwork=. , 1 >. +: n  NB. minimal
  lwork=. , 1 >. n * 130 [^:(x +./@e. 'vV') 33  NB. optimal
  (|: L: 0) 6 7 9 { zgeev_jlapack2_ (, jobVl) ; (, jobVr) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 0j0) ; Vl ; (, 1 >. # Vl) ; Vr ; (, 1 >. # Vr) ; (lwork $ 0j0) ; lwork ; ((+: n) $ 0.0) ; , _1
    NB. (|:) doesn't affect to w
)
