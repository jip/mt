NB. Raise to an integer power
NB.
NB. gepow      Raise a general matrix to integer power(s)
NB. dipow      Raise a diagonalizable matrix to integer
NB.            power(s)
NB. hepow      Raise a Hermitian (symmetric) matrix to
NB.            integer power(s)
NB.
NB. testgepow  Test gepow by general matrix given
NB. testdipow  Test dipow by diagonalizable matrix given
NB. testhepow  Test hepow by Hermitian (symmetric) matrix
NB.            given
NB. testpow    Adv. to make verb to test xxpow by matrix of
NB.            generator and shape given
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gepow
NB.
NB. Description:
NB.   Raise a general matrix A to integer power(s)
NB.
NB. Syntax:
NB.   P=. p gepow A
NB. where
NB.   A  - n×n-matrix, a general matrix
NB.   p  - sh-array of positive integers, power(s)
NB.   P  - sh×n×n-array if r>0,
NB.        n×n-array    if r=0, a matrix A in power(s) p
NB.   sh - r-vector of non-negative integers, the shape of p
NB.   r  ≥ 0, the rank of p
NB.
NB. References:
NB. [1] http://www.jsoftware.com/jwiki/Essays/Linear_Recurrences

gepow=: 4 : 0

  mpi3=. mp/ ^: (3 = (# @ $))       NB. apply mp/ only to 3-rank arrays (stiff rank)

  pl=. i. >: <. (2 & ^.) (>./) , x  NB. powers list: 2^i
  pc=. mp~ ^: pl y                  NB. powers cache: A^2^i
  pb=. (< @ I. @ (|. " 1) @ #:) x   NB. pl bits boxed array of shape sh

  pc (mpi3 @ ({~ >)) " 3 0 pb       NB. extract and mp A's powers for each pl atom
)

NB. ---------------------------------------------------------
NB. dipow
NB.
NB. Description:
NB.   Raise a diagonalizable matrix to integer power(s)
NB.
NB. Syntax:
NB.   P=. p dipow (iL ; v ; L)
NB.   P=. p dipow (R ; v ; iR)
NB. where
NB.   L  - n×n-matrix, unitary (orthogonal), rows are left
NB.        eigenvectors of A
NB.   R  - n×n-matrix, unitary (orthogonal), columns are
NB.        right eigenvectors of A
NB.   v  - n-vector, eigenvalues of A
NB.   iL = L^_1
NB.   iR = R^_1
NB.   p  - sh-array of positive integers, power(s)
NB.   P  - sh×n×n-array if r>0,
NB.        n×n-array    if r=0, a matrix A in power(s) p
NB.   sh - r-vector of non-negative integers, the shape of p
NB.   r  ≥ 0, the rank of p
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((-: ~.) v) +. ((-: ct) A)  NB. A must be normal (diagonalizable)
NB.   C -: L mp R                 NB. LAPACK doesn't guarantee (C -: idmat # A)
NB.   iL -: R %"1 c               NB. see [1]
NB.   iR -: L % c                 NB. see [1]
NB.   A -: iL mp V mp L
NB.   A -: iL mp v * L
NB.   A -: R mp V mp iR
NB.   A -: R mp v * iR
NB.   F -: diagmat"1 f
NB.   P -: iL mp"2 F mp"2 L
NB.   P -: iL mp"2 f *"1 2 L
NB.   P -: R mp"2 F mp"2 iR
NB.   P -: R mp f *"1 2 iR
NB.   P -: p dipow iL ; v ; L
NB.   P -: p dipow R ; v ; iR
NB. where
NB.   'Lh v R'=. geev A           NB. conjugated left eigenvectors in columns ; eigenvalues ; right eigenvectors in columns
NB.   L=. ct Lh                   NB. left eigenvectors in rows
NB.   iL=. %. L
NB.   iR=. %. R
NB.   V=. diagmat v
NB.   c=. L mp"1 |: R
NB.   C=. diagmat c
NB.   P=. p gepow A
NB.   f=. v ^1 0 p
NB.   F=. p gepow V
NB.
NB. References:
NB. [1] http://icl.cs.utk.edu/lapack-forum/viewtopic.php?p=985#p985
NB.     LAPACK/ScaLAPACK Development ‹ DGEEVX and left eigenvectors
NB.     Julien Langou, Fri Dec 22, 2006 5:15 pm

dipow=: (0 {:: ]) mp"2 ([ ^"1 0~ 1 {:: ]) (*"1 2) 2 {:: ]

NB. ---------------------------------------------------------
NB. hepow
NB.
NB. Description:
NB.   Raise a Hermitian matrix to integer power(s)
NB.
NB. Syntax:
NB.   P=. p dipow (v ; R)
NB. where
NB.   v - n-vector, eigenvalues of A
NB.   R - n×n-matrix, columns are eigenvectors of A
NB.   p  - sh-array of positive integers, power(s)
NB.   P  - sh×n×n-array if r>0,
NB.        n×n-array    if r=0, a matrix A in power(s) p
NB.   sh - r-vector of non-negative integers, the shape of p
NB.   r  ≥ 0, the rank of p
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((-: ct) A)     NB. A must be Hermitian (symmetric)
NB.   iR -: ct R
NB.   A -: R mp V mp iR
NB.   A -: R mp v * iR
NB.   F -: diagmat"1 f
NB.   P -: R mp"2 F mp"2 iR
NB.   P -: R mp"2 f *"1 2 iR
NB.   P -: p hepow v ; R
NB. where
NB.   'v R'=. heev A  NB. eigenvalues ; right eigenvectors in columns
NB.   iR=. %. R
NB.   V=. diagmat v
NB.   P=. p gepow A
NB.   f=. v ^1 0 p
NB.   F=. p gepow V

hepow=: (0 {:: ]) mp"2 ([ ^"1 0~ 1 {:: ]) *"1 2 ct@(0 {:: ])

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgepow
NB.
NB. Description:
NB.   Test gepow by general matrix given
NB.
NB. Syntax:
NB.   testgepow A
NB. where
NB.   A - n×n-matrix
NB.
NB. Notes:
NB. - fixed powers vector (p -: 5 7) is used

testgepow=: 3 : 0
  rcond=. (norm1 con (getrilu1p@getrflu1p)) y

  ('5 7&gepow' tmonad (]`]`(rcond"_)`(_."_)`(_."_))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testdipow
NB.
NB. Description:
NB.   Test dipow by diagonalizable matrix given
NB.
NB. Syntax:
NB.   testdipow A
NB. where
NB.   A - n×n-matrix, diagonalizable
NB.
NB. Notes:
NB. - fixed powers vector (p -: 5 7) is used
NB.
NB. TODO:
NB. - replace geev_jlapack_ by geev

testdipow=: 3 : 0
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'geev'

  rcond=. (norm1 con (getrilu1p@getrflu1p)) y

  'Lh v R'=. geev_jlapack_ y         NB. do eigendecomposition
  v=. j./ (*"1 (-@*@{.)) |: +. v     NB. for each v[i] in v, flip sign of v[i] if Re(v[i])>0, to force
                                     NB. A to be negative definite, this will avoid NaN error in diexp
  assert ((-: ~.) v) +. ((-: ct) A)  NB. A must be normal (diagonalizable)
  L=. ct Lh                          NB. restore L
  iR=. L ([ % (mp"1 |:)) R           NB. reconstruct R^_1 , see [1] in diexp

  ('5 7&dipow' tmonad (]`]`(rcond"_)`(_."_)`(_."_))) (R ; v ; iR)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhepow
NB.
NB. Description:
NB.   Test hepow by Hermitian (symmetric) matrix given
NB.
NB. Syntax:
NB.   testhepow A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric)
NB.
NB. TODO:
NB. - replace heev_jlapack_ by heev

testhepow=: 3 : 0
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'heev'

  rcond=. (norm1 con (hetripl@hetrfpl)) y

  NB. - do eigendecomposition: 'v R'=. heev A
  NB. - for each v[i] in v, flip sign of v[i] if Re(v[i])>0,
  NB.   to force A to be negative definite, this will avoid
  NB.   NaN error in heexp
  NB. - save adjusted boxed duplet back into y
  y=. ((j./@(*"1 (-@*@{.))@:|:@:+.&.>)`] ag) heev_jlapack_ y

  ('5 7&hepow' tmonad (]`]`(rcond"_)`(_."_)`(_."_))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpow
NB.
NB. Description:
NB.   Adv. to make verb to test xxpow by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testpow
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     (? @ $ 0:) testpow_mt_ 200 200
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testpow_mt_ 200 200
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testpow_mt_ 200 200

testpow=: 1 : 'EMPTY_mt_ [ ((testhepow_mt_ @ (u hemat_mt_)) [ (testdipow_mt_ @ (u dimat_mt_ u)) [ (testgepow_mt_ @ u)) ^: (=/)'
