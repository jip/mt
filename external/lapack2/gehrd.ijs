require 'math/lapack2'

NB. Description:
NB.   Reduce a square matrix to upper Hessenberg form by an
NB.   unitary (orthogonal) similarity transformation
NB.
NB. Syntax:
NB.   'HQf tau'=. xgehrd ilo ; ihi ; A
NB. where
NB.   ilo ∈ [1,max(1,ihi)], IO starting row and column,
NB.         1-based
NB.   ihi ∈ [min(ilo,n),n], IO ending row and column, 1-based
NB.   A   - n×n-matrix, a matrix to reduce, upper triangular
NB.         in rows and columns outside ilo:ihi
NB.   HQf - n×n-matrix, H and Qf combined
NB.   tau - (n-1)-vector, scalar factors of elementary
NB.         reflectors applied to A, in elements ilo:ihi
NB.   H   - n×n-matrix, upper Hessenberg in rows and columns
NB.         ilo:ihi and upper triangular outside
NB.   Qf  - n×n-matrix, columns below the first subdiagonal
NB.         with the tau represent the Q in the factored
NB.         form
NB.   Q   - n×n-matrix, unitary (orthogonal), which is
NB.         defined as the product of (ihi-ilo) elementary
NB.         reflectors
NB.   n   ≥ 0, the size of A, HQf, H, Qf, Q and tau
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dgehrd=: 3 : 0
  'ilo ihi A'=. y
  n=. # A
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) A
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 4160 32&p.`1:@.(2&>) n  NB. optimal
  (|: L: 0) 4 6 { dgehrd_jlapack2_ (, n) ; (, ilo) ; (, ihi) ; (|: A) ; (, 1 >. n) ; ((<: n) $ 0.0) ; (lwork $ 0.0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)

zgehrd=: 3 : 0
  'ilo ihi A'=. y
  n=. # A
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) A
  NB. lwork=. , 1 >. n  NB. minimal
  lwork=. , 4160 32&p.`1:@.(2&>) n  NB. optimal
  (|: L: 0) 4 6 { zgehrd_jlapack2_ (, n) ; (, ilo) ; (, ihi) ; (|: A) ; (, 1 >. n) ; ((<: n) $ 0.0) ; (lwork $ 0j0) ; lwork ; , _1
    NB. (|:) doesn't affect to tau
)
