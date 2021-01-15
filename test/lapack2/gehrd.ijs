require 'math/lapack2'

NB. Description:
NB.   Reduce a general square matrix to upper Hessenberg form
NB.   by an unitary (orthogonal) similarity transformation
NB.
NB. Syntax:
NB.   'HQf tau'=. xgehrd ilo ; ihi ; A
NB. where
NB.   ilo ∈ [1,max(1,ihi)], IO starting row and column,
NB.         1-based
NB.   ihi ∈ [min(ilo,n),n], IO ending row and column, 1-based
NB.   A   - n×n-matrix, a matrix to reduce, the upper
NB.         triangular in rows and columns outside ilo:ihi
NB.   HQf - n×n-matrix, H and Qf combined
NB.   tau - (n-1)-vector, the scalar factors of elementary
NB.         reflectors applied to A, in elements ilo:ihi
NB.   H   - n×n-matrix, the upper Hessenberg in rows and
NB.         columns ilo:ihi and upper triangular outside
NB.   Qf  - n×n-matrix, columns below the first subdiagonal
NB.         with the tau represent the Q in the factored
NB.         form
NB.   Q   - n×n-matrix, the unitary (orthogonal), which is
NB.         defined as the product of (ihi-ilo) elementary
NB.         reflectors
NB.   n   ≥ 0, the size of A, HQf, H, Qf, Q and tau
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dgehrd=: 3 : 0
  'ilo ihi A'=. y
  n=. # A
  assert. (= <.)                          ilo , ihi
  assert. 1 0&=`((0 , n)&I. , <:/)@.(* n) ilo , ihi
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_) A
  select. 3!:0 A
    case. JCMPX do. A=. 9 o. A
    case. JFL   do.
    case.       do. A=. A + 0.0
  end.
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 4160 32&p.`1:@.(2&>) n  NB. optimal
  cdrc=. dgehrd_jlapack2_ (, n) ; (, ilo) ; (, ihi) ; (|: A) ; (, 1 >. n) ; ((<: n) $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 4 6 { cdrc  NB. (|:) doesn't affect to tau
)

zgehrd=: 3 : 0
  'ilo ihi A'=. y
  n=. # A
  assert. (= <.)                          ilo , ihi
  assert. 1 0&=`((0 , n)&I. , <:/)@.(* n) ilo , ihi
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_) A
  if. JCMPX ~: 3!:0 A do. A=. A + 0j0 end.
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 4160 32&p.`1:@.(2&>) n  NB. optimal
  cdrc=. zgehrd_jlapack2_ (, n) ; (, ilo) ; (, ihi) ; (|: A) ; (, 1 >. n) ; ((<: n) $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 4 6 { cdrc  NB. (|:) doesn't affect to tau
)
