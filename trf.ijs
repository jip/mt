NB. trf.ijs
NB. Triangular factorization for general (LU) and Hermitian
NB. (symmetric) positive definite matrix (Cholesky)
NB.
NB. getrflu   LU factorization of a general matrix, both
NB.           blocked and unblocked versions
NB. getrful   UL factorization of a general matrix, both
NB.           blocked and unblocked versions
NB. rgetrflu  LU factorization of a general matrix, recursive
NB.           version
NB. rgetrful  UL factorization of a general matrix, recursive
NB.           version
NB. potrfl    Cholesky factorization of a Hermitian
NB.           (symmetric) positive definite matrix for lower
NB.           triangular factor, both blocked and unblocked
NB.           versions
NB. potrfu    Cholesky factorization of a Hermitian
NB.           (symmetric) positive definite matrix for upper
NB.           triangular factor, both blocked and unblocked
NB.           versions
NB. rpotrfl   Cholesky factorization of a Hermitian
NB.           (symmetric) positive definite matrix for lower
NB.           triangular factor, recursive version
NB. rpotrfu   Cholesky factorization of a Hermitian
NB.           (symmetric) positive definite matrix for upper
NB.           triangular factor, recursive version
NB.
NB. Version: 1.0.0 2008-11-01
NB. Copyright: Igor Zhuravlov, igor at uic.dvgu.ru
NB. License: GNU GPL

coclass 'mt'

NB. =========================================================
NB. Local constants

NB. =========================================================
NB. Local verbs

iomaxm=: (i.>./) @: |  NB. IO element from list y with max magnitude

NB. ---------------------------------------------------------
NB. kcgetrflu
NB. LU Crout algorithm's kernel for column orientation
NB.
NB. Syntax:
NB.   W=. ios kcgetrflu Z
NB. where
NB.   Z     - m×n general matrix
NB.   ios   - (j , (j+1) , ... , (j+nb-1))
NB.   W     - (A - B * C)
NB.   nb    - block size, nb=1 (non-blocked version) or nb≥1
NB.           (blocked version)
NB.   j     - index of row, along with nb it defines Z
NB.           partitioning, j is multiplier of nb in range:
NB.           0 ≤ j < min(m,n)
NB.   m     ≥ 0
NB.   n     ≥ 0
NB.   A B C - submatrices extracted from Z at step (j/nb):
NB.                  j  nb (n-j-nb)
NB.           j      *  C  *
NB.           (m-j)  B  A  *

kcgetrflu=: ((_1 , ({. @ [)) ,: (((- {.)~ #) , (# @ [)))`(_1 0 ,: (((- , ]) {.)~ #))`(((0 , {.) ,: ({. , #)) @ [) kernel3 (-`mp)

NB. ---------------------------------------------------------
NB. krgetrflu
NB. LU Crout algorithm's kernel for row orientation
NB.
NB. Syntax:
NB.   W=. ios krgetrflu Z
NB. where
NB.   Z     - m×n general matrix
NB.   ios   - (j , (j+1) , ... , (j+nb-1))
NB.   W     - (A - B * C)
NB.   nb    - block size, nb=1 (non-blocked version) or nb≥1
NB.           (blocked version)
NB.   j     - index of row, along with nb it defines Z
NB.           partitioning, j is multiplier of nb in range:
NB.           0 ≤ j < min(m,n)
NB.   m     ≥ 0
NB.   n     ≥ 0
NB.   A B C - submatrices extracted from Z at step (j/nb):
NB.                    j  nb (n-j-nb)
NB.           j        *  *  C
NB.           nb       B  *  A
NB.           (m-j-nb) *  *  *

krgetrflu=: ((({. @ [) , _1:) ,: ((# @ [) , ((- ({. + #))~ ({: @ $))))`((({. , 0:) ,: (# , {.)) @ [)`(0 _1 ,: (({. @ [) , ((- ({. + #))~ ({: @ $)))) kernel3 (-`mp)

NB. ---------------------------------------------------------
NB. kpotrfl
NB. Cholesky algorithm's kernel for lower triangular factor
NB. case: L*L' = Z
NB.
NB. Syntax:
NB.   W=. ios kpotrfl Z
NB. where
NB.   Z       -  n×n Hermitian (symmetric) positive definite
NB.              matrix
NB.   ios     -: (j , (j+1) , ... , (j+nb-1))
NB.   W       -: [B;D] - [A;C] * A'
NB.   nb      -  block size, nb=1 (non-blocked version) or
NB.              nb≥1 (blocked version)
NB.   j       -  index of row, defines Z partitioning, it is
NB.              multiplier of nb in range: 0 ≤ j < min(m,n)
NB.   n       ≥  0
NB.   A B C D -  submatrices extracted from Z at step (j/nb):
NB.                       j  nb (n-j-nb)
NB.              j        *  *  *
NB.              nb       A  B  *
NB.              (n-j-nb) C  D  *

kpotrfl=: ((_1 , ({.@[)) ,: ((- {.)~ #) , (#@[))`(_1 0 ,: ((- {.)~ #) , ({.@[))`((# ((] , 0:) ,: ,) {.)@[) kernel3 (-`(mp ct))

NB. ---------------------------------------------------------
NB. kpotrfu
NB. Cholesky algorithm's kernel for upper triangular factor
NB. case: U*U' = Z
NB.
NB. Syntax:
NB.   W=. ios kpotrfu Z
NB. where
NB.   Z       -  n×n Hermitian (symmetric) positive definite
NB.              matrix
NB.   ios     -: (j , (j+1) , ... , (j+nb-1))
NB.   W       -: [D;B] - [C;A] * A'
NB.   nb      -  block size, nb=1 (non-blocked version) or
NB.              nb≥1 (blocked version)
NB.   j       -  index of row, defines Z partitioning, it is
NB.              multiplier of nb in range: 0 ≤ j < min(m,n)
NB.   n       ≥  0
NB.   A B C D -  submatrices extracted from Z at step (j/nb):
NB.                       j  nb (n-j-nb)
NB.              j        *  D  C
NB.              nb       *  B  A
NB.              (n-j-nb) *  *  *

kpotrfu=: (((((0 , [) ,: (+ , ])) #)~ {.)@[)`(((0 _1 ,: (] , -)) ({. + #))~ #)`(_1 0 ,: ((2 & $) @ # @ [)) kernel3 (-`(mp ct))

NB. =========================================================
NB. Interface verbs

NB. LU Crout algo
NB. 'p LU'=. [nb] getrflu matrix
NB. where nb - block size, nb>0, default is NB_GETRF

getrflu=: (BGETRF&$: : (4 : 0)) " 0 2
  mn=. <./ 'm n'=. $ y
  p=. i. m

  if. n > 1 do.
    for_j. range 0 , (mn - 1) , x do.
      jb=. x <. mn - j
      iosc=. (j + jb) ht2i j
      'p2 rc2'=. 1 getrflu iosc kcgetrflu y   NB. update, then factorize diagonal and subdiagonal blocks by non-blocked version
      dp=. (i. j) , (j + p2)             NB. current iteration's rows permutation (back to y frame from rc2 frame)
      p=. dp C. p                        NB. adjust p by p2: ((p2 { (iosr { p)) iosr } p)
      y=. dp C. y                        NB. apply rows permutation starting from j-th row
      y=. rc2 (< (m ht2i j) ; iosc) } y  NB. write back diag. and subdiag. blocks and undo p2 permutation in it
      NB. extract Ljj
      NB. update diag. and behind diag. blocks
      NB. solve (Ljj*X=C)
      rc2=. (iosc krgetrflu y) getrslx1 (jb {. rc2)
      y=. rc2 (< iosc ; (n ht2i (j + jb))) } y       NB. write back X in place of rc2
    end.
  else.
    dp=. 0 ii2cp iomaxm y                NB. find pivot and form cycle permutation
    p=. dp C. p                          NB. adjust p by p2: ((p2 { (iosr { p)) iosr } p)
    y=. ({. 0 } (%"1 {.)) dp C. y        NB. apply rows permutation starting from j-th row, then scale by head
  end.

  p ; y
)

NB. recursive blocked

NB. References:
NB. [1] J. J. Dongarra, S. Hammarling, D. W. Walker. (1996)
NB.     Key Concepts For Parallel Out-Of-Core LU Factorization.
NB.     LAPACK Working Note 110, UT-CS-96-324, April 1996.
NB. [2] J. J. Dongarra, E. F. D'Azevedo. (1997)
NB.     The Design and Implementation of the Parallel
NB.     Out-of-core ScaLAPACK LU, QR, and Cholesky
NB.     Factorization Routines.
NB.     LAPACK Working Note 118, UT-CS-97-347, January 1997.

rgetrflu=: (3 : 0) " 2
  'm n'=. mn=. $ y
  if. (RGETRF < 7!:5 < 'y') *. (n > 1) do.
    k=. m <. n2=. >. -: n
    'p1 LU1'=. (m , k) rgetrflu ;. 0 y     NB. factorize A's 1st block column recursively
    y=. (0 _1 ,: m , (n - k)) (p1 & C.) ;. 0 y  NB. apply p1 to A's 2nd block column
    A12=. (k , (n - k)) ] ;. 0 y
    LU11=. (2 $ k) ] ;. 0 LU1
    U12=. A12 getrslx1 LU11                     NB. solve L11*U12=A12 for U12
    L21=. (_1 0 ,: (m - k) , k) ] ;. 0 LU1
    A22=. (_1 0 ,: (m - k) , (n - k)) ] ;. 0 y
    'p2 LU22'=. rgetrflu (A22 - L21 mp U12)       NB. factorize updated A22 recursively
    dp=. (i. k) , (k + p2)     NB. (k + p2) C. p1
    p=. dp C. p1
    p ; (dp C. LU1) ,. (U12 , LU22)
  else.
    getrflu y
  end.
)

NB. ---------------------------------------------------------
NB. potrfl                                                  2
NB. Hermitian (symmetric) positive definite matrix Cholesky
NB. blocked factorization by a lower triangular factor
NB.   L * L' = A
NB.
NB. Syntax:
NB.   L=. potrfl A
NB. where
NB.   A - N×N-matrix, Hermitian (symmetric), positive
NB.       definite; strictly upper triangle is not referenced
NB.   L  - N×N-matrix, lower triangular Cholesky factor
NB.
NB. If:
NB.   >>>>>>>>>>>>>>>>>'HQ tau'=. gehrd A ; ss
NB.   D=. diagmat s
NB.   Dinv=. %. D
NB. then
NB.   Dinv -: diagmat % s
NB.   As -: D mp Ap mp Dinv
NB.   As -: Ap (] * (% " 1)) s
NB.
NB. Notes:
NB.   - FLOPs: 
NB.   - strict upper triangle is not modified
NB.
NB. References:
NB. [1] 
NB. [2] 

potrfl=: (BPOTRF&$: : (4 : 0)) " 0 2
  n=. # y
  if. n > 2 do.
    for_j. range 0 , (n - 1) , x do.
      jb=. x <. n - j
      ios=. (j + jb) ht2i j
      bd=. ios kpotrfl y
      b=. 1 potrfl jb {. bd
      y=. b (< ;~ ios) } y
      d=. (jb }. bd) getrsxu ct b
      y=. d (< (n ht2i (j + jb)) ; ios) } y
    end.
NB.---  elseif. n = 2 do.
NB.---    'a bconj b c'=. v=. , y
NB.---    assert. a (*. & (> & 0)) c
NB.---    e=. b % d=. %: a
NB.---    (d , bconj) ,: (e , (%: c - (abssq e)))
NB.---  elseif. do.
  else.
NB.---    assert. y > 0                         NB. check for positive definite
    %: y
  end.
)

NB. U*U'=A
NB.   A -: ((+/ .*) (+ @ |:)) potf2u_mt_ A
NB.   A -: ((+/ .*) (+ @ |:)) potf2_mt_ &. ((] ;. 0) :. (] ;. 0)) A

potrfu=: (BPOTRF&$: : (4 : 0)) " 0 2
  n=. # y
  if. n > 2 do.
    for_j. range (n - 1) , 0 , x do.
      jb=. x <. n - j
      ios=. (j + jb) ht2i j
      db=. ios kpotrfu y
      b=. (,.~ _1 , jb) (1 & potrfu) ;. 0 db
      y=. b (< ;~ ios) } y
      d=. ((j , jb) ] ;. 0 db) getrsxl ct b
      y=. d ((i. j) ; ios) } y
    end.
  elseif. n = 2 do.
    'a b bconj c'=. v=. , y
    assert. a (*. & (> & 0)) c
    e=. b % f=. %: c
    ((%: a - (abssq e)) , e) ,: (bconj , f)
  elseif. do.
    assert. y > 0                         NB. check for positive definite
    %: y
  end.
)

potf2l_old=: (3 : 0) " 2
  n=. # y
  for_j. i. n do.                           NB. traverse diagonal forward, j=0..(n-1)
    ioc=. < (j + i. n - j) ; j              NB. A[i,j], i=j..(n-1)
    c=. ioc { y                             NB. column on and under A[j][j]
    A21=. (_1 0 ,: (n - j) , j) (] ;. 0) y  NB. (n-j)×j-submatrix with upper right corner at A[j,j-1]
    c=. c - (mp (+ @ {.)) A21               NB. GEMV(-1,A21,(conjugated A21 1st row),1,c)
    Ljj=. 9 o. {. c                         NB. L[j][j] is always real
    assert. Ljj > 0                         NB. check for positive definite
    y=. (c % %: Ljj) ioc } y                NB. scale and write back column c
  end.
)

potrfl_old=: (3 : 0) " 2
  n=. # y
  nb=. 32           NB. choosed experimentally
  for_j. range 0 , (n-1) , nb do.
    jb=. nb <. n - j
    bd=. ((2 $ j) ,: ((n - j) , jb)) ] ;. 0 y
    ac=. (_1 0 ,: ((n - j) , j)) ] ;. 0 y
    ah=. (0 0 ,: (jb , j)) ct ;. 0 ac
    bd=. bd - ac mp ah
    b=. (2 $ jb) potf2l_old ;. 0 bd
    y=. b (< ;~ j + i. jb) } y
    ios=. (j + i. jb) ; ((j + jb) ([ + (i. @ (-~))) n)
    y=. 0 (< ios) } y
    d=. (_1 0 ,: ((n - (j + jb)) , jb)) ] ;. 0 bd
    x=. d mp (128!:1) ct b
    y=. x (< |. ios) } y
  end.
)

NB. ---------------------------------------------------------
NB. rpotrfl                                                 2
NB. Hermitian (symmetric) positive definite matrix Cholesky
NB. recursive blocked factorization by a lower triangular
NB. factor
NB.   L * L' = A
NB.
NB. Syntax:
NB.   L=. rpotrfl A
NB. where
NB.   A - N×N-matrix, Hermitian (symmetric), positive
NB.       definite; strictly upper triangle is not referenced
NB.   L  - N×N-matrix, lower triangular Cholesky factor
NB.
NB. If:
NB.   >>>>>>>>>>>>>>>>>'HQ tau'=. gehrd A ; ss
NB.   D=. diagmat s
NB.   Dinv=. %. D
NB. then
NB.   Dinv -: diagmat % s
NB.   As -: D mp Ap mp Dinv
NB.   As -: Ap (] * (% " 1)) s
NB.
NB. Notes:
NB.   - FLOPs: 
NB.   - recursive calls: 3*k*(2^k)-2^(k+1)+3,
NB.     where k  = ⌈log_2(⌈n/nb⌉)⌉
NB.           nb = CPU_CACHE, block size
NB.   - strict upper triangle is not modified
NB.
NB. References:
NB. [1] F. G. Gustavson, I. Jonsson, "Minimal-storage high-
NB.     performance Cholesky factorization via blocking
NB.     recursion", IBM J. Res. Develop., vol. 44, No. 6,
NB.     2000, pp. 823-850.

rpotrfl=: (3 : 0) " 2
  if. RPOTRF < 7!:5 < 'y' do.
    n=. # y
    p=. >. -: n
    q=. n - p
    L0=. (2 $ p) rpotrfl ;. 0 y
    Y=. (0 _1 ,: (p , q)) ] ;. 0 y
    L2H=. Y getrslx L0
    Z=. (,.~ _1 , q) ] ;. 0 y
    (L0 ,. Y) ,: ([ ,. (rpotrfl @ (Z & -) @ mp)~ ct) L2H
  else.
NB. smoutput 'rpotrfl' ; '$ y' ; ($ y)
    potrfl y
  end.
)

rpotrfl_old=: (3 : 0) " 2
  if. RPOTRF < 7!:5 < 'y' do.
    n=. # y
    p=. >. -: n
    q=. n - p
    Yh=. (_1 0 ,: (q , p)) ] ;. 0 y
    Z=. (_1 _1 ,: (2 $ q)) ] ;. 0 y
    L0=. (2 $ p) rpotrfl_old ;. 0 y
    L2=. Yh mp (128!:1) (ct L0)
    L1=. rpotrfl Z - (mp ct) L2
    L0 , L2 ,. L1
  else.
    potf2l_old y
  end.
)

NB. =========================================================
NB. Test suite

NB. 'name type rows cols error time space'=. name vextract ttrf A
ttrf=: 1 : 0
:
  'm n'=. $ y
  't s'=. ts 'LUmix=. ' , x , ' y'
smoutput 'x' ; x ; 'm' ; m ; 'n' ; n ; '($ y)' ; ($ y) ; '$ LUmix' ; ($ LUmix)
  e=. (norm1 y - u LUmix) % (n * FP_EPS * norm1 y)
  x ; (datatype y) ; m ; n ; e ; t ; s
)

NB. result=. testgetrf A
NB. where result is boxed table with columns: name type rows cols error time space
testgetrf=: (3 : 0) " 2
NB.smoutput 'testgetrf ($ y)' ; ($ y)
  r=. 0 7 $ a:
  r=. r , 'getrflu'        (UL2L1 mp UL2U) ttrf y
  r=. r , 'rgetrflu'       (UL2L1 mp UL2U) ttrf y
  r=. r , 'getrf_jlapack_' (((mp & >) @ (2 & {.)) invperm_jlapack_ (2 & {::)) ttrf y
)

NB. result=. testpotrf A
NB. where result is boxed table with columns: name type rows cols error time space
testpotrf=: (3 : 0) " 2
  llh2a=. (mp ct) @ UL2L
  uuh2a=. (mp ct) @ UL2U
NB.smoutput 'testpotrf ($ y)' ; ($ y)
  r=. 0 7 $ a:
  r=. r , 'potrfl'         llh2a ttrf y
  r=. r , 'potrfl_old'     llh2a ttrf y
  r=. r , 'potrf_jlapack_' llh2a ttrf y
  r=. r , 'potrfu'         uuh2a ttrf y
  r=. r , '1 & potrfl'     llh2a ttrf y
  r=. r , 'potf2l_old'     llh2a ttrf y
  r=. r , 'rpotrfl'        llh2a ttrf y
  r=. r , 'rpotrfl_old'    llh2a ttrf y
)

testtrf=: (3 : 0) " 1
  r=. 0 7 $ a:
  r=. r , (testpotrf @ (mkun mkpo) @ {. ^: (=/)) 2 $ y
  r=. r , (testpotrf @ (mkor mkpo) @ {. ^: (=/)) 2 $ y
  r=. r , (testgetrf @ rand01) y
)


NB. test for errors, measure time and space
NB.   tpotrf mat

require '~addons/math/lapack/lapack.ijs'
need_jlapack_ 'potrf'

zztpotrf=: 3 : 0

  n=. # y
  in=. i. n
  I=. =/~ in                           NB. identity rectangular matrix
  U=. <:/~ in                          NB. upper triangle rectangular matrix
  L=. >:/~ in                          NB. lower triangle rectangular matrix

  LU2L=. (n & ({. " 1)) @: (L & *)     NB. L is m×min(m,n) matrix
  LU2U=. (n & {.) @: (U & *)           NB. U is min(m,n)×n matrix
  error=. %: @ (+/) @: , @: | @: -

  smoutput 'routine match error time space'
  timespace=. ts 'L=. rpotrfl y'
  smoutput 'rpotrfl ' , (": y (-: , error) (mp ct) (LU2L L)) , ' ' , (": timespace)
  timespace=. ts 'L=. 1 potrfl y'
  smoutput '1 potrfl ' , (": y (-: , error) (mp ct) (LU2L L)) , ' ' , (": timespace)
  timespace=. ts 'L=. potf2l_old y'
  smoutput 'potf2l_old ' , (": y (-: , error) (mp ct) (LU2L L)) , ' ' , (": timespace)
  timespace=. ts 'L=. potrfl y'
  smoutput 'potrfl ' , (": y (-: , error) (mp ct) (LU2L L)) , ' ' , (": timespace)
  timespace=. ts 'L=. potrfl_old y'
  smoutput 'potrfl_old ' , (": y (-: , error) (mp ct) (LU2L L)) , ' ' , (": timespace)
  timespace=. ts 'L=. potrf_jlapack_ y'
  smoutput 'LAPACK potrf ' , (": y (-: , error) (mp ct) (LU2L L)) , ' ' , (": timespace)
)

NB. test and measure interface verb
NB.   testpotrf n

zztestpotrf=: 3 : 0
   Af=. mkor_mt_ mkpo_mt_ y
NB.   Ac=. mkun_mt_ mkpo_mt_ y
  tpotrf Af
)

NB. =========================================================
Note 'chol testing and timing'
   A =: _4]\ 33 7j_8 7j_10 3j_4, 7j8 28 2j4 _10j_11, 7j10 2j_4 22 3j3, 3j4 _10j11 3j_3 16
   A -: ((+/ .*) (+ @ |:)) (potf2_mt_ A)
1
   B=: 3 3$ 1 4 6  4 0 3  6 3 2
   potf2_mt_ B
error

   A100=. ((+/ .*) (+ @ |:)) j./ ? 2 100 100 $ 100
   A100 -: ((+/ .*) (+ @ |:)) (potf2_mt_ A100)
1
   A500=. ((+/ .*) (+ @ |:)) j./ ? 2 500 500 $ 100
   A500 -: ((+/ .*) (+ @ |:)) (potf2_mt_ A500)
1
   A1000=. ((+/ .*) (+ @ |:)) j./ ? 2 1000 1000 $ 100
   A1000 -: ((+/ .*) (+ @ |:)) (potf2_mt_ A1000)
1
   ts=: 6!:2, 7!:2@]

   A50d=. ((+/ .*) |:) ? 50 50 $ 100
   (potf2_mt_ (+/ @: (+/) @: | @: -) Choleski_mt_) A50d
6.33016e_8
   (A50d (+/ @: (+/) @: | @: -) ((+/ .*) (+ @ |:)) (potf2_mt_ A50d))
1.73168e_9
   (A50d (+/ @: (+/) @: | @: -) ((+/ .*) (+ @ |:)) (Choleski_mt_ A50d))
2.53472e_5

   A1000f=: ((+/ .*) |:) 10.1 + ? (2 $ 1000) $ 10
   A1000c=: ((+/ .*) (+ @ |:)) j./ ? (2 , (2 $ 1000)) $ 10

   ts & > 'z=. potf2l_mt_ A1000f';'z=. potrfl_mt_ A1000f';'z=. rpotrfl_mt_ A1000f'
4.32201 1.78581e7
 1.0309 1.35258e7
1.76633 2.93628e7
   ts & > 'z=. potf2l_mt_ A1000c';'z=. potrfl_mt_ A1000c';'z=. rpotrfl_mt_ A1000c'
7.73323 5.14207e7
2.44583 2.70426e7
3.96861 5.87229e7

)

NB. зависимость времени исполнения potf2l от размера матрицы

require 'plot'

mp=: +/ . *            NB. matrix product
lio=: + ` (* i.)/ " 1  NB. integers grid (2{y) steps from (0{y) by (1{y)

writetable=: 4 : '(toHOST,(": x),"1 LF) 1!:2 y'  NB. table writetable < 'filename'
readtable=: 3 : '>0 ". each cutopen toJ 1!:1 y'  NB. table=. readtable < 'filename'

flops_potf=: 0 0.166667 1.41667 0.166667 & p.

flops_Chol=: 4 : 0
  if. x>:y do.
    flops_potf y
  else.
    x (flops_Chol + (0 0 1 0.5 0.1875 p. ])) (>. -: y)
  end.
)

plot2d=: 3 : 0
  vnb=. lio 1 1 100 NB. , y  NB. NB[i]
  vt=. (# vnb) $ 0      NB. execution time
NB.  A=: ((+/ .*) |:) 10.1 + ? (2 $ y) $ 10
  A=: ((+/ .*) (+ @ |:)) j./ ? (2 , (2 $ y)) $ 10
  for_b. vnb do.
    nb=: b NB. >>>>>> CPU_CACHE=: b
    t=. (6!:2) 'z=. nb potrfl2 A'
NB.    f=. ((flops_Chol y) % t) % 1e6
    if. 0 = (y%10) | b do.
      smoutput 'b_index' ; ((": >: b_index) , '/' , (": y)) ; 'NB' ; b ; 'time' ; t
    end.
    vt=. t b_index } vt
NB.    vf=. f b_index } vf
  end.

  pd 'reset'
  pd 'title Cholesky factorization via potrfl2'
  pd 'type line'
  pd 'xticpos ' , ": lio 10 10 , ((# vnb) % 10)
NB.  pd 'grids 0 1'
NB.  pd 'yrange 0 0.075'
  pd 'xcaption Size, N'
  pd 'ycaption Execution time, sec'
  pd vnb;vt
  pd 'show'
  pd 'save bmp 1280 1000 chol_2d.bmp'

  vt writetable < 'vt.dump'
NB.  vff writetable < 'vff.dump'
)

NB. зависимость времени исполнения от размера квадрата, факторизуемого посредством potf2l

plot3d=: 3 : 0
  vn=. lio 1 1 , y   NB. N[i]
  mtf=. (2 $ y) $ 0  NB. execution time for float
NB.  mff=. (2 $ y) $ 0  NB. FLOPs for float
  for_n. vn do.
    Af=: ((+/ .*) |:) ? (2 $ n) $ 100
NB.    Ac=: ((+/ .*) (+ @ |:)) j./ ? (2 , (2 $ n)) $ 100
    for_b. lio 1 1 , n do. NB. log: ((, (>. @ -: @ {:)) ^: (1 ~: {:) ^: _ ) n
      >>>>>>>>>>>>CPU_CACHE=: b
      tf=. (6!:2) 'z=. rpotrfl Af'
NB.      ff=. ((nb flops_Chol n) % tf) % 1e6
      mtf=. tf (< n_index , b_index) } mtf
NB.      mff=. ff (< n_index , b_index) } mff
    end.
    if. 0 = (y%10) | n do.
      smoutput 'n_index' ; ((": >: n_index) , '/' , (": y)) ; 'n' ; n
    end.
  end.

  pd 'reset'
  pd 'title Cholesky factorization blocked recursion exec time, sec'
NB.  pd 'sub 1 2'
NB.  pd 'new'
  pd 'type surface'
  pd 'viewpoint 2.4 0.25 0.5' NB. _1.6 2.4 0.5'
  pd 'yticpos ' , ": lio 10 10 , (y % 10)
  pd 'grids 0 1'
  pd 'zrange 0 0.075'
  pd 'xcaption N'
  pd 'ycaption NB'
  pd vn;vn;mtf
NB.  pd 'use'
NB.  pd 'type surface'
NB.  pd 'viewpoint _1.6 2.4 0.5'
NB.  pd 'xcaption N'
NB.  pd 'ycaption NB'
NB.  pd vn;vn;mff
NB.  pd 'endsub'
  pd 'show'
  pd 'save bmp 1280 1000 chol_nb.bmp'

  mtf writetable < 'mtf.dump'
NB.  mff writetable < 'mff.dump'
)

plotflops=: 3 : 0
  require 'plot'
  mopcount=. 0 0 0 1r3e6 & p.
  vn=. 100 * i. 10
  vflopspersec=. 10 $ 0
  for_n. vn do.
    m=. ((+/ .*) (+ @ |:)) j./ ? (2 , (2 $ n)) $ 100
    time=. (6!:2) 'potf2 m'
    smoutput 'n' ; n ; 'time' ; time ; 'Mflops/s' ; ((mopcount n) % time) ; 'n_index' ; n_index
    vflopspersec=. ((mopcount n) % time) n_index } vflopspersec
  end.
  'title potf2 Cholesky factorization; ycaption MFLOPs/s; xcaption N' plot vn ; vflopspersec
)


NB. ---------------

NB. test for errors
NB.   tgetrf mat

require '~addons/math/lapack/lapack.ijs'
need_jlapack_ 'getrf'

zzztgetrf=: 3 : 0
  mn=. <./ 'm n'=. $ y
  im=. i. m
  in=. i. n
  I=. im =/ in                           NB. identity rectangular matrix
  U=. im <:/ in                          NB. upper triangle rectangular matrix
  L=. im >/ in                           NB. strict lower triangle rectangular matrix with unit diagonal elements

  LU2L=. (mn & ({. " 1)) @: (I + L & *)  NB. L is m×min(m,n) matrix
  LU2U=. (mn & {.) @: (U & *)            NB. U is min(m,n)×n matrix
  error=. %: @ (+/) @: , @: | @: -

  'p LU'=. rgetrf y
  smoutput 'rgetrf (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'p LU'=. 1 getrf y
  smoutput '1 getrf (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'p LU'=. getrf y
  smoutput 'getrf (match error): ' , ": (p C. y) (-: , error) (LU2L LU) mp (LU2U LU)
  'L U p'=. getrf_jlapack_ y
  smoutput 'LAPACK getrf (match error): ' , ": y (-: , error) p invperm_jlapack_~ L mp U
)

NB. measure time and space
NB.   mgetrf mat

zzzmgetrf=: 3 : 0
  ts=. 6!:2, 7!:2@]
  ts & > 'pLU=. rgetrf y';'pLU=. 1 getrf y';'pLU=. getrf y';'LUp=. getrf_jlapack_ y'
)

NB. test and measure interface verbs
NB.   testgetrf m,n
NB.   testgetrf n

zzztestgetrf=: 3 : 0
  y=. 2 {. y , {: y    NB. n->(n,n) , (m,n)->(m,n)
  (tgetrf ] mgetrf) 0.1 * ? y $ 100
)

NB. =========================================================
NB. Estimation suite

NB. зависимость времени исполнения potf2l от размера матрицы

require 'plot'

lio=: + ` (* i.)/ " 1  NB. integers grid (2{y) steps from (0{y) by (1{y)

writetable=: 4 : '(toHOST,(": x),"1 LF) 1!:2 y'  NB. Syntax: table writetable < 'filename'
readtable=: 3 : '>0 ". each cutopen toJ 1!:1 y'  NB. Syntax: table=. readtable < 'filename'

NB. Syntax: getrfNB (m,n)

zzzgetrfNB=: 3 : 0
  max=. >./ 'm n'=. y=. 2 {. y , {: y    NB. n->(n,n) , (m,n)->(m,n)
  vnb=. max ht2i 1      NB. NB[i]
  vt=. (# vnb) $ 0      NB. execution time
  A=: 0.1 * ? y $ 100
  for_nb. vnb do.
    gnb=. nb
    t=. (6!:2) 'z=. gnb getrf A'
    if. 0 = (max%10) | nb do.
      smoutput 'nb_index' ; ((": >: nb_index) , '/' , (": y)) ; 'NB' ; nb ; 'time' ; t
    end.
    vt=. t nb_index } vt
  end.

  pd 'reset'
  pd 'title LU factorization via getrf'
  pd 'type line'
  pd 'xticpos ' , ": range 10 , (# vnb) , 10
NB.  pd 'grids 0 1'
NB.  pd 'yrange 0 0.075'
  pd 'xcaption Size, NB'
  pd 'ycaption Execution time, sec'
  pd vnb;vt
  pd 'show'
  pd 'save bmp 1280 1000 getrf_2d.bmp'

  vt writetable < 'vt.dump'
)

NB. >>>>>>>>>>>>>>>>>>>>
NB. зависимость времени исполнения от размера квадрата, факторизуемого посредством potf2l

zzzplot3d=: 3 : 0
  vn=. lio 1 1 , y   NB. N[i]
  mtf=. (2 $ y) $ 0  NB. execution time for float
NB.  mff=. (2 $ y) $ 0  NB. FLOPs for float
  for_n. vn do.
    Af=: ((+/ .*) |:) ? (2 $ n) $ 100
NB.    Ac=: ((+/ .*) (+ @ |:)) j./ ? (2 , (2 $ n)) $ 100
    for_b. lio 1 1 , n do. NB. log: ((, (>. @ -: @ {:)) ^: (1 ~: {:) ^: _ ) n
      >>>>>>>>>>>>CPU_CACHE=: b
      tf=. (6!:2) 'z=. rpotrfl Af'
NB.      ff=. ((nb flops_Chol n) % tf) % 1e6
      mtf=. tf (< n_index , b_index) } mtf
NB.      mff=. ff (< n_index , b_index) } mff
    end.
    if. 0 = (y%10) | n do.
      smoutput 'n_index' ; ((": >: n_index) , '/' , (": y)) ; 'n' ; n
    end.
  end.

  pd 'reset'
  pd 'title Cholesky factorization blocked recursion exec time, sec'
NB.  pd 'sub 1 2'
NB.  pd 'new'
  pd 'type surface'
  pd 'viewpoint 2.4 0.25 0.5' NB. _1.6 2.4 0.5'
  pd 'yticpos ' , ": lio 10 10 , (y % 10)
  pd 'grids 0 1'
  pd 'zrange 0 0.075'
  pd 'xcaption N'
  pd 'ycaption NB'
  pd vn;vn;mtf
NB.  pd 'use'
NB.  pd 'type surface'
NB.  pd 'viewpoint _1.6 2.4 0.5'
NB.  pd 'xcaption N'
NB.  pd 'ycaption NB'
NB.  pd vn;vn;mff
NB.  pd 'endsub'
  pd 'show'
  pd 'save bmp 1280 1000 chol_nb.bmp'

  mtf writetable < 'mtf.dump'
NB.  mff writetable < 'mff.dump'
)
