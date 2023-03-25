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
NB. - the verbs below are loaded into the current locale
NB. - the computed eigenvectors are normalized to have
NB.   Euclidean norm equal to 1 and largest component real

dgeev=: 4 : 0
  'jobVl jobVr'=. x
  assert. jobVl ,&('nNvV' e.~ {.) jobVr
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_ , isreal_jlapack2_) y
  if. JFL ~: 3!:0 y do. y=. 9 o. y end.
  n=. # y
  Vl=. (0 0 [^:('nN' e.~ {. jobVl) $ y) $ 0.0
  Vr=. (0 0 [^:('nN' e.~ {. jobVr) $ y) $ 0.0
  ldVl=. , 1 >. # Vl
  ldVr=. , 1 >. # Vr
  NB. lwork=. , 1 >. n * 4 [^:(x +./@e. 'vV') 3  NB. minimal
  lwork=. , 1 >. n * 130 [^:(x +./@e. 'vV') 34  NB. optimal
  cdrc=. dgeev_jlapack2_ (, jobVl) ; (, jobVr) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 0.0) ; (n $ 0.0) ; Vl ; ldVl ; Vr ; ldVr ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  'wr wi Vl Vr'=. (|: L: 0) 6 7 8 10 { cdrc  NB. (|:) doesn't affect to wr and wi
  w=. wr j. wi
  if. # cx=. I. wi ~: 0 do.
    if. 'vV' e.~ {. jobVl do. Vl=. cx cxpair_jlapack2_ Vl end.
    if. 'vV' e.~ {. jobVr do. Vr=. cx cxpair_jlapack2_ Vr end.
  end.
  w ; Vl ; Vr
)

zgeev=: 4 : 0
  'jobVl jobVr'=. x
  assert. jobVl ,&('nNvV' e.~ {.) jobVr
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  n=. # y
  Vl=. (0 0 [^:('nN' e.~ {. jobVl) $ y) $ 0j0
  Vr=. (0 0 [^:('nN' e.~ {. jobVr) $ y) $ 0j0
  ldVl=. , 1 >. # Vl
  ldVr=. , 1 >. # Vr
  NB. , lwork=. , 1 >. +: n  NB. minimal
  lwork=. , 1 >. n * 130 [^:(x +./@e. 'vV') 33  NB. optimal
  cdrc=. zgeev_jlapack2_ (, jobVl) ; (, jobVr) ; (, n) ; (|: y) ; (, 1 >. n) ; (n $ 0j0) ; Vl ; ldVl ; Vr ; ldVr ; (lwork $ 0j0) ; lwork ; ((+: n) $ 0.0) ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 6 7 9 { cdrc  NB. (|:) doesn't affect to w
)
