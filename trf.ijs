NB. trf.ijs
NB. Triangular factorization of matrix
NB.
NB. getrfpl1u  Factorization with partial pivoting of a
NB.            general matrix: inv(P) * L1 * U = A, where P
NB.            is permutation matrix, L1 is unit lower
NB.            triangular and U is upper triangular
NB. getrflu1p  Factorization with partial pivoting of a
NB.            general matrix: L * U1 * inv(P) = A, where P
NB.            is permutation matrix, L is lower triangular
NB.            and U1 is unit upper triangular
NB. getrfpu1l  Factorization with partial pivoting of a
NB.            general matrix: inv(P) * U1 * L = A, where P
NB.            is permutation matrix, U1 is unit upper
NB.            triangular and L is lower triangular
NB. getrful1p  Factorization with partial pivoting of a
NB.            general matrix: U * L1 * inv(P) = A, where P
NB.            is permutation matrix, U is upper triangular
NB.            and L1 is unit lower triangular
NB.
NB. hetrfpu    Factorization with full pivoting of a
NB.            Hermitian (symmetric) matrix:
NB.            inv(P) * U1 * T * U1' * P = A, where P is
NB.            permutation matrix, U1 is unit upper
NB.            triangular and T is tridiagonal
NB. hetrfpl    Factorization with full pivoting of a
NB.            Hermitian (symmetric) matrix:
NB.            inv(P) * L1 * T * L1' * P = A, where P is
NB.            permutation matrix, L1 is unit lower
NB.            triangular and T is tridiagonal
NB.
NB. potrfu     Cholesky factorization of a Hermitian
NB.            (symmetric) positive definite matrix:
NB.            U * U' = A, where U is upper triangular
NB. potrfl     Cholesky factorization of a Hermitian
NB.            (symmetric) positive definite matrix:
NB.            L * L' = A, where L is lower triangular
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-10

coclass 'mt'

NB. =========================================================
NB. Local definitions

iofmaxm=: (i.>./) @: |  NB. IO 1st element with max magnitude from list y
iolmaxm=: (i:>./) @: |  NB. IO last element with max magnitude from list y

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. getrfpl1u
NB.
NB. Description:
NB.   Triangular factorization of a general matrix:
NB.     P * L1 * U = A
NB.
NB. Syntax
NB.     'p L1U'=. getrfpl1u A
NB. where
NB.   A  - m*n matrix
NB.
NB. Notes:
NB. - this is a right-looking version of algorithm

NB. ---------------------------------------------------------
NB. Description:
NB.   Single step of non-blocked version of algorithms
NB.
NB. Syntax
NB.   'pi1 pfxi1 sfxLi1 sfxRi1=. getf2xxxxstep (pi ; pfxi ; sfxLi ; sfxRi)
NB. where
NB.   pi     - p(i), m-vector after i-th and before (i+1)-th
NB.            step, permutation
NB.   pfxi   - pfx(i), i*n-matrix, first i rows of matrix
NB.            L1U(i)
NB.   sfxLi  - sfxL(i), the left part of 
NB.   sfxRi  - sfxR(i), the right part of 
NB.   pi1    - p(i+1), m-vector after (i+1)-th step,
NB.            permutation
NB.   pfxi1  - pfx(i+1), (i+1)*n-matrix, first i rows of
NB.            matrix L1U(i+1)
NB.   sfxLi1 - sfxL(i+1), the left part of 
NB.   sfxRi1 - sfxR(i+1), the right part of 

NB. ---------------------------------------------------------
NB. Description:
NB.   Single step of algorithms
NB.
NB. Syntax
NB.   'pi1 pfxi1 sfxLi1 sfxRi1'=. getf2xxxxstep (pi ; pfxi ; sfxLi ; sfxRi)
NB. where
NB.   pfxi   - pfx(i), ##########
NB.
NB. a00 a00 a00 a00      a00 a00 a00 a00
NB. a10 a11 a12 a12      a00 a00 a00 a00
NB. a10 a21 a22 a22      a10 a10 a11 a12
NB. a10 a21 a22 a22      a10 a10 a21 a22
NB.

NB. ---------------------------------------------------------
NB. Verb:          Solves:                  Syntax:
NB. getrfpl1u      inv(P) * L1 * U = A      'pi Af'=. getrfpl1u A
NB. getrflu1p      L * U1 * inv(P) = A      'pi Af'=. getrflu1p A
NB. getrfpu1l      inv(P) * U1 * L = A      'pi Af'=. getrfpu1l A
NB. getrful1p      U * L1 * inv(P) = A      'pi Af'=. getrful1p A
NB.
NB. Description:
NB.   triangular factorization with partial pivoting of a
NB.   general matrix
NB. where:
NB.   A  - m×n-matrix, containing either U, U1, L or L1
NB.   Af - m×n-matrix, contains either U and L1, or U1 and L
NB.   pi - m-vector or n-vector, inversed rows (columns)
NB.        permutation of A
NB.   U  - ?×?-matrix, upper triangular factor
NB.   U1 - ?×?-matrix, unit upper triangular matrix (diagonal is not saved)
NB.   L  - ?×?-matrix, lower triangular matrix
NB.   L1 - ?×?-matrix, unit lower triangular matrix (diagonal is not saved)
NB.   m  ≥ 0
NB.   n  ≥ 0
NB.
NB. If:
NB.   'pi Af'=. getrfpl1u A
NB.   p=. /: pi
NB.   Pi=. ((i.#pi)=/pi)
NB.   P=. ((i.#p)=/p)
NB.   L1=. trl1 Af
NB.   U=. tru Af
NB. then (with appropriate comparison tolerance)
NB.   Pi -: %. P
NB.   Pi -: |: P
NB.   A -: p { L1 mp U
NB.   A -: p C. L1 mp U
NB.   A -: P mp L1 mp U
NB.
NB. Notes:
NB.   - FLOPs: 
NB.   - recursive calls: 3*k*(2^k)-2^(k+1)+3,
NB.     where k = ⌈log_2(n)⌉

NB. (-: (clean@((/: @ (0&({::))) C."1 ((trl mp tru1)@(1&({::))))@getrflu1p))
NB. 'ip LU1'=. getrflu1p A

getrflu1p=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. n) ; y
  elseif. 1=m do.
    dpi=. 0 ios2cp iofmaxm y
    pi=. dpi C. i. n
    y=. ((] 0:} %) (0&({,))) dpi (C."1) y                          NB. permute single column (row), scale by head, keep head unscaled
    pi ; y
  elseif. do.
    k=. n (<. >.@-:) m
    'pia Afa'=. getrflu1p k {. y                                   NB. factorize 1st block recursively
    y=. pia (C."1) k }. y                                          NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trtrsxu1 & (k & ({."1))) y                         NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. getrflu1p y ((- (Afba & mp)) & (k & (}."1))) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. (i. k) , (k + pib)                                      NB. apply 2nd block's permutation to 1st block
    (dpib (C."1) pia) ; ((dpib (C."1) Afa) , (Afba ,. Afbb))       NB. assemble solution
  end.
)

NB. (-: (clean@((/: @ (0&({::))) C. ((trl1 mp tru)@(1&({::))))@getrfpl1u)) A
NB. 'ip L1U'=. getrfpl1u A

getrfpl1u=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. m) ; y
  elseif. 1 = n do.
    dpi=. 0 ios2cp iofmaxm y
    pi=. dpi C. i. m
    y=. ((] 0:} %) (0&({,))) dpi C. y                          NB. permute single column (row), scale by head, keep head unscaled
    pi ; y
  elseif. do.
    k=. m (<. >.@-:) n
    'pia Afa'=. getrfpl1u k {."1 y                             NB. factorize 1st block recursively
    y=. pia C. k }."1 y                                        NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trtrsl1x & (k & {.)) y                         NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. getrfpl1u y ((- (mp & Afba)) & (k & }.)) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. (i. k) , (k + pib)                                  NB. apply 2nd block's permutation to 1st block
    (dpib C. pia) ; ((dpib C. Afa) ,. (Afba , Afbb))           NB. assemble solution
  end.
)

NB. (-: (clean@((/: @ (0&({::))) C. (((tru1~ (-~/@$)) mp (trl~ (-~/@$)))@(1&({::))))@getrfpu1l)) A
NB. 'ip U1L'=. getrfpu1l A

getrfpu1l=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. m) ; y
  elseif. 1 = n do.
    dpi=. (<:m) ios2cp iolmaxm y
    pi=. dpi C. i. m
    y=. ((] _1:} %) (_1&({,))) dpi C. y                           NB. permute single column (row), scale by head, keep head unscaled
    pi ; y
  elseif. do.
    k=. m (<. >.@-:) n
    'pia Afa'=. getrfpu1l (-k) {."1 y                             NB. factorize 1st block recursively
    y=. pia C. (-k) }."1 y                                        NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trtrsu1x & ((-k) & {.)) y                         NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. getrfpu1l y ((- (mp & Afba)) & ((-k) & }.)) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. pib , ((m-k) + i. k)                                   NB. FIXME! remove 2nd; apply 2nd block's permutation to 1st block
    (dpib C. pia) ; ((Afbb , Afba) ,. (dpib C. Afa))              NB. assemble solution
  end.
)

NB. (-: (clean@((/: @ (0&({::))) C."1 (((tru~ (-~/@$)) mp (trl1~ (-~/@$)))@(1&({::))))@getrful1p)) A
NB. 'ip UL1'=. getrful1p A

getrful1p=: 3 : 0
  'm n'=. sh=. $ y
  if. 0 e. sh do.
    (i. n) ; y
  elseif. 1 = m do.
    dpi=. (<:n) ios2cp iolmaxm y
    pi=. dpi C. i. n
    y=. ((] _1:} %) (_1&({,))) dpi (C."1) y                           NB. permute single column (row), scale by head, keep head unscaled
    pi ; y
  elseif. do.
    k=. n (<. >.@-:) m
    'pia Afa'=. getrful1p (-k) {. y                                 NB. factorize 1st block recursively
    y=. pia (C."1) (-k) }. y                                        NB. apply 1st block's permutation to 2nd block, purge original y, reuse name 'y'
    Afba=. Afa (trtrsxl1 & ((-k) & ({."1))) y                         NB. calculate 2nd block's 1st sub-block
    'pib Afbb'=. getrful1p y ((- (Afba & mp)) & ((-k) & (}."1))) Afa  NB. update 2nd block's 2nd sub-block and factorize it recursively
    dpib=. pib , ((n-k) + i. k)                                       NB. FIXME! remove 2nd; apply 2nd block's permutation to 1st block
    (dpib (C."1) pia) ; ((Afbb ,. Afba) , (dpib (C."1) Afa))          NB. assemble solution
  end.
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. tgetrf
NB.
NB. Description:
NB.   Test general matrix TRF algorithms by matrix given
NB.
NB. Syntax:
NB.   tgetrf A
NB. where
NB.   A  - m×n-matrix
NB.
NB. Note:
NB. -

tgetrf=: 3 : 0
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'getrf'

  rcond=. ((_."_)`(norm1 con getri) @. (=/@$)) y  NB. meaninigful for square matrices only

  ('getrf_jlapack_' tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- (((mp & >)/ @ }:) invperm_jlapack_~      (2 & {::) ))) % (FP_EPS*((norm1*c)@[)))))) y NB. berr := ||P*L1*U - A ||/(ε*||A||*n)

  ('getrflu1p'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((/: @ (0 & {::)) C."1 (( trl            mp  tru1          )@(1 & {::))))) % (FP_EPS*((norm1*#)@[)))))) y NB. berr := ||L*U1*P - A ||/(ε*||A||*m)
  ('getrfpl1u'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((/: @ (0 & {::)) C.   (( trl1           mp  tru           )@(1 & {::))))) % (FP_EPS*((norm1*c)@[)))))) y NB. berr := ||P*L1*U - A ||/(ε*||A||*n)
  ('getrfpu1l'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((/: @ (0 & {::)) C.   (((tru1~ (-~/@$)) mp (trl ~ (-~/@$)))@(1 & {::))))) % (FP_EPS*((norm1*c)@[)))))) y NB. berr := ||P*U1*L - A ||/(ε*||A||*n)
  ('getrful1p'      tmonad (]`]`(rcond"_)`(_."_)`(((norm1@(- ((/: @ (0 & {::)) C."1 (((tru ~ (-~/@$)) mp (trl1~ (-~/@$)))@(1 & {::))))) % (FP_EPS*((norm1*#)@[)))))) y NB. berr := ||U*L1*P - A ||/(ε*||A||*m)

  EMPTY
)

NB. ---------------------------------------------------------
NB. thetrf
NB.
NB. Description:
NB.   Test Hermitian matrix TRF algorithms by matrix given
NB.
NB. Syntax:
NB.   thetrf A
NB. where
NB.   A - n×n-matrix, is used to produce Hermitian matrix

thetrf=: 3 : 0
  NB. A=. (+ ct) y                         NB. convert matrix type: ge -> he
  NB. rcond=. (norm1 con getri) A

  EMPTY
)

NB. ---------------------------------------------------------
NB. testtrf
NB.
NB. Description:
NB.   Adv. to make verb to test TRF algorithms by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testtrf
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.            mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testtrf 150 100

testtrf=: 1 : 'EMPTY [ (thetrf ^: (=/@$) [ tgetrf) @ u'
