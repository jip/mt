require 'math/lapack2'

NB. Description:
NB.   Solve a system of linear equations with a non-complex
NB.   symmetric matrix, and factor the last one
NB.
NB. Syntax:
NB.   'DT1 ipiv X'=. uplo dsysv_aa A ; B
NB. where
NB.   uplo - scalar, character, case-insensitive:
NB.            'L' - use lower triangle of DT1 only, form is:
NB.                    P * L1 * D * L1^T * P^T = A
NB.            'U' - use upper triangle of DT1 only, form is:
NB.                    P * U1 * D * U1^T * P^T = A
NB.   A    - n×n-matrix, real, the symmetric to be factored
NB.          to DT1 and ipiv
NB.   B    - n×nrhs-matrix, real, RHS
NB.   DT1  - n×n-matrix, D and T1 combined
NB.   ipiv - n-vector, integer, the pivot indices that define
NB.          P
NB.   X    - n×nrhs-matrix, real, solutions of equation:
NB.            A * X = B
NB.   D    - n×n-matrix, real, the symmetric tridiagonal
NB.          (opposite diagonal not stored)
NB.   T1   - n×n-matrix, either L1 or U1
NB.   L1   - n×n-matrix, real, the unit lower triangular
NB.          (unit diagonal not stored)
NB.   U1   - n×n-matrix, real, the unit upper triangular
NB.          (unit diagonal not stored)
NB.   P    - n×n-matrix, boolean, the permutation matrix
NB.   n    ≥ 0, the order of system
NB.   nrhs ≥ 0, the number of RHS
NB.
NB. Notes:
NB. - the verbs below are loaded into the current locale

dsysv_aa=: 4 : 0
  'A B'=. y
  'n nrhs'=. $ B
  assert. *./ (1 = # x) , x e. 'lLuU'
  assert. (ismatrix_jlapack2_ *. issquare_jlapack2_ *. isreal_jlapack2_ *. n = #) A
  assert. (ismatrix_jlapack2_ *.                       isreal_jlapack2_         ) B
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
  NB. lwork=. , >./ 1 , (+: n) , (_2 3 p. n)  NB. minimal
  lwork_dsytrf_aa=. 1 >. n * >: 64
  lwork_dsytrs_aa=. 1 >. _2 3 p. n
  lwork=. , lwork_dsytrf_aa >. lwork_dsytrs_aa  NB. optimal
  cdrc=. dsysv_aa_jlapack2_ (, x) ; (, n) ; (, nrhs) ; (|: A) ; ld ; (n $ 00) ; (|: B) ; ld ; (lwork $ 0.0) ; lwork ; , _1
  assert. 0 = _1 {:: cdrc
  (|: L: 0) 4 6 7 { cdrc  NB. (|:) doesn't affect to ipiv
)
