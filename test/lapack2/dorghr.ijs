require 'math/lapack2'

NB. Description:
NB.   Generate a real orthogonal matrix from output of DGEHRD
NB.
NB. Syntax:
NB.   Q=. dorghr ilo ; ihi ; A ; tau
NB. where
NB.   ilo ∈ [1,max(1,ihi)], IO starting row and column,
NB.         1-based
NB.   ihi ∈ [min(ilo,n),n], IO ending row and column, 1-based
NB.   A   - n×n-matrix, real, contains Qf
NB.   Qf  - n×n-matrix, columns ilo:ihi below the first
NB.         subdiagonal contain elementary reflectors as
NB.         returned by DGEHRD, with the tau it represents
NB.         the Q in the factored form
NB.   tau - (n-1)-vector, real, scalar factors of elementary
NB.         reflectors as returned by DGEHRD
NB.   Q   - n×n-matrix, real, orthogonal, which is defined as
NB.         the product of (ihi-ilo) elementary reflectors
NB.         H(i) of order n:
NB.           Q = Π{H(i),i=ilo:ihi-1}
NB.           H(i) = I - v[i] * τ[i] * v[i]'
NB.   n   ≥ 0, the size of A, Q and Qf
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dorghr=: 3 : 0
  'ilo ihi A tau'=. y
  n=. # A
  assert. (ismatrix_jlapack2_ , issquare_jlapack2_) A
  assert. (isvector_jlapack2_ , (<: n) = #        ) tau
  NB. lwork=. , 1 >. ihi-ilo  NB. minimal
  lwork=. , 1 >. 32 * ihi - ilo  NB. optimal
  |: 4 {:: dorghr_jlapack2_ (, n) ; (, ilo) ; (, ihi) ; (|: A) ; (, 1 >. n) ; tau ; (lwork $ 0.0) ; lwork ; , _1
)
