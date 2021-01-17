require 'math/lapack2'

NB. Description:
NB.   Balance a pair of general square matrices
NB.
NB. Syntax:
NB.   'Abal Bbal ilo ihi lscale rscale'=. job xggbal A ;  B
NB.   'Abal Bbal ilo ihi lscale rscale'=. job xggbal A ,: B
NB. where
NB.   job       - scalar, character, case-insensitive:
NB.                 'N' - to init ilo, ihi and xscale
NB.                 'P' - to permute only
NB.                 'S' - to scale only
NB.                 'B' - to do both permute and scale
NB.   A,B       - n×n-matrix, a matrix pair
NB.   Abal,Bbal - n×n-matrix, a balanced matrix pair (A,B)
NB.   ilo       ∈ [1,max(1,ihi)], IO starting row and column,
NB.               1-based
NB.   ihi       ∈ [min(ilo,n),n], IO ending row and column,
NB.               1-based
NB.   lscale    - n-vector, float, details of the
NB.               permutations and scaling factors applied
NB.               to the left side of A and B
NB.   rscale    - n-vector, float, details of the
NB.               permutations and scaling factors applied
NB.               to the right side of A and B
NB.   n         ≥ 0, the size of A, B, Abal, Bbal, lscale and
NB.               rscale
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dggbal=: 4 : 0
  assert. (1 = # x) *. x e. 'nNpPsSbB'
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
  ld=. , 1 >. n
  work=. ((6 * n) >.^:(x e. 'sSbB') 1) $ 0.0
  cdrc=. dggbal_jlapack2_ (, x) ; (, n) ; (|: A) ; ld ; (|: B) ; ld ; (, 0) ; (, 0) ; (n $ 0.0) ; (n $ 0.0) ; work ; , _1
  assert. 0 = _1 {:: cdrc
  'A B ilo ihi lscale rscale'=. 3 5 7 8 9 10 { cdrc
  (|: A) ; (|: B) ; ({. ilo) ; ({. ihi) ; lscale ; rscale
)

zggbal=: 4 : 0
  assert. (1 = # x) *. x e. 'nNpPsSbB'
  'A B'=. y
  n=. # A
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_         ) A
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. n = #) B
  if. JCMPX ~: 3!:0 A do. A=. A + 0j0 end.
  if. JCMPX ~: 3!:0 B do. B=. B + 0j0 end.
  ld=. , 1 >. n
  work=. ((6 * n) >.^:(x e. 'sSbB') 1) $ 0.0
  cdrc=. zggbal_jlapack2_ (, x) ; (, n) ; (|: A) ; ld ; (|: B) ; ld ; (, 0) ; (, 0) ; (n $ 0.0) ; (n $ 0.0) ; work ; , _1
  assert. 0 = _1 {:: cdrc
  'A B ilo ihi lscale rscale'=. 3 5 7 8 9 10 { cdrc
  (|: A) ; (|: B) ; ({. ilo) ; ({. ihi) ; lscale ; rscale
)
