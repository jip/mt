NB. lu.ijs
NB. LU factorization
NB.
NB. getrf  LU factorization of a general matrix

coclass 'mt'

NB. =========================================================
NB. Local verbs

iomaxm=: (i.>./) @: |  NB. IO element from list y with max magnitude

NB. ---------------------------------------------------------
NB. LU Crout algorithm's kernel for row and column orientations
NB.
NB. at step (j%nb), j=0..(min(m,n)-1),
NB. nb=1 for non-blocked version, nb>1 for blocked version,
NB. returns A-B*C, where:
NB.
NB.                 column oriented case     row oriented case
NB.
NB. y partition:           j  nb (n-j-nb)             j  nb (n-j-nb)
NB.                 j      *  C  *           j        *  *  C
NB.                 (m-j)  B  A  *           nb       B  *  A
NB.                                          (m-j-nb) *  *  *

kcgetrf=: ((_1 , ({. @ [)) ,: (((- {.)~ #) , (# @ [)))`(_1 0 ,: (((- , ]) {.)~ #))`(((0 , {.) ,: ({. , #)) @ [) kernel3 (-`mp)

krgetrf=: ((({. @ [) , _1:) ,: ((# @ [) , ((- ({. + #))~ ({: @ $))))`((({. , 0:) ,: (# , {.)) @ [)`(0 _1 ,: (({. @ [) , ((- ({. + #))~ ({: @ $)))) kernel3 (-`mp)

NB. =========================================================
NB. Interface verbs

NB. LU Crout algo
NB. 'p LU'=. [nb] getrf matrix
NB. where nb - block size, nb>0, default is NB_GETRF

getrf=: (NB_GETRF&$: : (4 : 0)) " 0 2
  mn=. <./ 'm n'=. $ y
  p=. i. m

  if. n > 1 do.
    for_j. range 0 , (mn - 1) , x do.
      jb=. x <. mn - j
      iosc=. (j + jb) ht2i j
      'p2 rc2'=. 1 getrf iosc kcgetrf y   NB. update, then factorize diagonal and subdiagonal blocks by non-blocked version
      dp=. (i. j) , (j + p2)             NB. current iteration's rows permutation (back to y frame from rc2 frame)
      p=. dp C. p                        NB. adjust p by p2: ((p2 { (iosr { p)) iosr } p)
      y=. dp C. y                        NB. apply rows permutation starting from j-th row
      y=. rc2 (< (m ht2i j) ; iosc) } y  NB. write back diag. and subdiag. blocks and undo p2 permutation in it
      NB. extract Ljj
      NB. update diag. and behind diag. blocks
      NB. solve (Ljj*X=C)
      rc2=. (iosc krgetrf y) getrslx1 (jb {. rc2)
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

rgetrf=: (3 : 0) " 2
  'm n'=. mn=. $ y
  if. (CPU_CACHE < 7!:5 < 'y') *. (n > 1) do.   NB. (m > 2) *. (n > 2) do. NB. 
    k=. m <. n2=. >. -: n
    'p1 LU1'=. (0 0 ,: m , k) rgetrf ;. 0 y     NB. factorize A's 1st block column recursively
    y=. (0 _1 ,: m , (n - k)) (p1 & C.) ;. 0 y  NB. apply p1 to A's 2nd block column
    A12=. (0 0 ,: k , (n - k)) ] ;. 0 y
    LU11=. (0 0 ,: 2 $ k) ] ;. 0 LU1
    U12=. A12 getrslx1 LU11                     NB. solve L11*U12=A12 for U12
    L21=. (_1 0 ,: (m - k) , k) ] ;. 0 LU1
    A22=. (_1 0 ,: (m - k) , (n - k)) ] ;. 0 y
    'p2 LU22'=. rgetrf (A22 - L21 mp U12)       NB. factorize updated A22 recursively
    dp=. (i. k) , (k + p2)     NB. (k + p2) C. p1
    p=. dp C. p1
    p ; (dp C. LU1) ,. (U12 , LU22)
  else.
    getrf y
  end.
)

NB. =========================================================
NB. Test suite

NB. test for errors
NB.   tgetrf mat

require '~addons/math/lapack/lapack.ijs'
need_jlapack_ 'getrf'

tgetrf=: 3 : 0
  mn=. <./ 'm n'=. $ y
  im=. i. m
  in=. i. n
  I=. im =/ in                           NB. identity rectangular matrix
  U=. im <:/ in                          NB. upper triangle rectangular matrix
  L=. im >/ in                           NB. strict lower triangle rectangular matrix

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

mgetrf=: 3 : 0
  ts=. 6!:2, 7!:2@]
  ts & > 'pLU=. rgetrf y';'pLU=. 1 getrf y';'pLU=. getrf y';'LUp=. getrf_jlapack_ y'
)

NB. test and measure interface verbs
NB.   testgetrf m,n
NB.   testgetrf n

testgetrf=: 3 : 0
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

getrfNB=: 3 : 0
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
