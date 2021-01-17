require 'math/lapack2'

NB. Description:
NB.   Generate a complex unitary matrix from output of ZGEHRD
NB.
NB. Syntax:
NB.   Q=. zunghr ilo ; ihi ; A ; tau
NB. where
NB.   ilo ∈ [1,max(1,ihi)], IO starting row and column,
NB.         1-based
NB.   ihi ∈ [min(ilo,n),n], IO ending row and column, 1-based
NB.   A   - n×n-matrix, contains Qf
NB.   Qf  - n×n-matrix, columns ilo:ihi below the first
NB.         subdiagonal contain elementary reflectors as
NB.         returned by ZGEHRD, with the tau it represents
NB.         the Q in the factored form
NB.   tau - (n-1)-vector, the scalar factors of elementary
NB.         reflectors as returned by ZGEHRD
NB.   Q   - n×n-matrix, complex, unitary, which is defined as
NB.         the product of (ihi-ilo) elementary reflectors of
NB.         order n:
NB.           Q = Π{H(i),i=ilo:ihi-1}
NB.           H(i) = I - v[i] * τ[i] * v[i]'
NB.   n   ≥ 0, the size of A, Q and Qf
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

zunghr=: 3 : 0
  'ilo ihi A tau'=. y
  n=. # A
  assert. (= <.)                          ilo , ihi
  assert. 1 0&=`((0 , n)&I. , <:/)@.(* n) ilo , ihi
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_) A
  assert. (isvector_jlapack2_ *. (<: n) = #        ) tau
  if. JCMPX ~: 3!:0 A do. A=.   A   + 0j0 end.
  if. JCMPX ~: 3!:0 A do. tau=. tau + 0j0 end.
  NB. lwork=. , 1 >. ihi-ilo  NB. minimal
  lwork=. , 1 >. 32 * ihi - ilo  NB. optimal
  cdrc=. zunghr_jlapack2_ (, n) ; (, ilo) ; (, ihi) ; (|: A) ; (, 1 >. n) ; tau ; (lwork $ 0j0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  |: 4 {:: cdrc
)
