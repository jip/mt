require 'math/lapack2'

NB. Description:
NB.   Balance a general square matrix
NB.
NB. Syntax:
NB.   'Abal ilo ihi scale'=. job xgebal A
NB. where
NB.   job   - scalar, character, case-insensitive:
NB.             'N' - to init ilo, ihi and scale
NB.             'P' - to permute only
NB.             'S' - to scale only
NB.             'B' - to do both permute and scale
NB.   A     - n×n-matrix to balance
NB.   Abal  - n×n-matrix, a balanced matrix A
NB.   ilo   ∈ [1,max(1,ihi)], IO starting row and column,
NB.           1-based
NB.   ihi   ∈ [min(ilo,n),n], IO ending row and column,
NB.           1-based
NB.   scale - n-vector, float, details of the permutations
NB.           and scaling factors applied to A
NB.   n     ≥ 0, the size of A, Abal and scale
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dgebal=: 4 : 0
  assert. (1 = # x) *. x e. 'nNpPsSbB'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_) y
  select. 3!:0 y
    case. JCMPX do. y=. 9 o. y
    case. JFL   do.
    case.       do. y=. y + 0.0
  end.
  n=. # y
  cdrc=. dgebal_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (, 0) ; (, 0) ; (n $ 0.0) ; , _1
  assert. 0 = _1 {:: cdrc
  'y ilo ihi scale'=. 3 5 6 7 { cdrc
  (|: y) ; ({. ilo) ; ({. ihi) ; scale
)

zgebal=: 4 : 0
  assert. (1 = # x) *. x e. 'nNpPsSbB'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_) y
  if. JCMPX ~: 3!:0 y do. y=. y + 0j0 end.
  n=. # y
  cdrc=. zgebal_jlapack2_ (, x) ; (, n) ; (|: y) ; (, 1 >. n) ; (, 0) ; (, 0) ; (n $ 0.0) ; , _1
  assert. 0 = _1 {:: cdrc
  'y ilo ihi scale'=. 3 5 6 7 { cdrc
  (|: y) ; ({. ilo) ; ({. ihi) ; scale
)
