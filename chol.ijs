NB. chol.ijs
NB. Cholesky factorization
NB.
NB. potf2l  Hermitian (symmetric) positive definite matrix
NB.         unblocked Cholesky factorization for a lower
NB.         triangular factor
NB. potf2u  Hermitian (symmetric) positive definite matrix
NB.         unblocked Cholesky factorization for a upper
NB.         triangular factor
NB.
NB. Version: 1.0.0 2008-09-01
NB. Copyright: Igor Zhuravlov, igor@uic.dvgu.ru
NB. License: GNU GPL

coclass 'pjlap'

NB. =========================================================
NB. Local constants

NB. =========================================================
NB. Local verbs

NB. =========================================================
NB. Interface verbs

NB. ---------------------------------------------------------
NB. potf2l                                                  2
NB. Hermitian (symmetric) positive definite matrix unblocked
NB. Cholesky factorization for a lower triangular factor
NB.   L * L' = A
NB.
NB. Syntax:
NB.   L=. potf2l A
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
NB.   - (3 2 1 1r2 p. n) flops
NB.   - 
NB.
NB. References:
NB. [1] 
NB. [2] 

NB. L*L'=A
NB. FLOPs: n/12*(n^2+17*n+2)
NB. занимает
NB. - 6464 для 11×11, 8512 для 12×12
NB. - 289280 для 127×127, 567808 для 128×128
NB. - 2.37722e6 для 362×362, 4.47437e6 для 363×363

potf2l=: (3 : 0) " 2
  n=. # y
  qio=. 2 2 ($ " 1) 2 0 1 2 ({ " 1) 0 ,. (lio n , _1 , n) ,. i. n  NB. quad IOS of A21 at step j
  for_j. i. n do.                      NB. traverse diagonal forward, j=0..(n-1)
    ioc=. < (j + i. n - j) ; j         NB. A[i,j], i=j..(n-1)
    c=. ioc { y                        NB. column on and under A[j][j]
    A21=. (j_index { qio) (] ;. 0) y   NB. (n-j)×j-submatrix with upper right corner at A[j,j-1]
    c=. c - (mp (+ @ {.)) A21          NB. GEMV(-1,A21,(conjugated A21 1st row),1,c)
    Ljj=. 9 o. {. c                    NB. L[j][j] is always real
    assert. Ljj > 0                    NB. check for positive definite
    y=. (c % %: Ljj) ioc } y
  end.

  y * (>:/~ i. n)                      NB. fill strict upper triangular by zeros
)

NB. U*U'=A
NB.   A -: ((+/ .*) (+ @ |:)) potf2u_pjlap_ A
NB.   A -: ((+/ .*) (+ @ |:)) potf2_pjlap_ &. ((] ;. 0) :. (] ;. 0)) A

potf2u=: (3 : 0) " 2
  n=. # y
  qio=. 2 2 ($ " 1) 1 2 1 (# " 1) 0 ,. (lio n , _1 , n) ,. i. n  NB. quad IOS of A12 at step j
  for_j. lio n , _1 , n do.            NB. traverse diagonal backward, j=n..1
    ioc=. < (i. j) ; <: j              NB. A[i,j-1], i=0..(j-1)
    c=. ioc { y                        NB. column on and above A[j][j]
    A12=. (j_index { qio) (] ;. 0) y   NB. j×(n-j)-submatrix with lower left corner at A[j-1,j]
    c=. c - (mp (+ @ {:)) A12          NB. GEMV(-1,A12,(conjugated A12 last row),1,c)
    Ljj=. 9 o. {: c                    NB. L[j][j] is always real
    assert. Ljj > 0                    NB. check for positive definite
    y=. (c % %: Ljj) ioc } y
  end.

  y * (<:/~ i. n)                      NB. fill strict lower triangular by zeros
)

potrfl=: (3 : 0) " 2
  n=. # y
  nb=. 32           NB. choosed experimentally
  for_j. range 0 , (n-1) , nb do.
    jb=. nb <. n-j
    bd=. ((j , j) ,: ((n - j) , jb)) ] ;. 0 y
    ac=. ((j , 0) ,: ((n - j) , j)) ] ;. 0 y
    ah=. ((j , 0) ,: (jb , j)) h ;. 0 y
    bd=. bd - ac mp ah
    b=. (2 $ jb) potf2l ;. 0 bd
    y=. b (< ;~ j + i. jb) } y
    ios=. (j + i. jb) ; ((j + jb) ([ + (i. @ (-~))) n)
    y=. 0 (< ios) } y
    d=. (_1 0 ,: ((n - (j + jb)) , jb)) ] ;. 0 bd
    x=. d mp (128!:1) h b
    y=. x (< |. ios) } y
  end.
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB.*tlartg v test lartg

NB.*tlartg v test lartg

tlartg=: (3 : 0) " 1
  'cs sn r'=. lartg y
  G=. 2 2 $ cs , sn , (- + sn) , cs
  (r,0) -: G mp y
)

NB. ---------------------------------------------------------
NB.*tlartg2 v test lartg

tlartg2=: (3 : 0) " 1
  erm=. 4 : '(| x - y) % ((FP_EPS_pjlap_ * | x) >. (FP_SFMIN_pjlap_ * FP_PREC_pjlap_))'
  'cs sn r'=. lartg y
  G=. 2 2 $ cs , sn , (- + sn) , cs
  r erm {. G mp y
)

NB. ---------------------------------------------------------
NB.*testlartg v test lartg

testlartg=: 3 : 0
  *./ tlartg ((9^2^2) , 2) $, (,"0/~) , (j./~) 1e_308 1e_200 1e_100 1e_50 10 1e50 1e100 1e200 1e308
)

NB. =========================================================
Note 'chol testing and timing'
   A =: _4]\ 33 7j_8 7j_10 3j_4, 7j8 28 2j4 _10j_11, 7j10 2j_4 22 3j3, 3j4 _10j11 3j_3 16
   A -: ((+/ .*) (+ @ |:)) (potf2_pjlap_ A)
1
   B=: 3 3$ 1 4 6  4 0 3  6 3 2
   potf2_pjlap_ B
error

   A100=. ((+/ .*) (+ @ |:)) j./ ? 2 100 100 $ 100
   A100 -: ((+/ .*) (+ @ |:)) (potf2_pjlap_ A100)
1
   A500=. ((+/ .*) (+ @ |:)) j./ ? 2 500 500 $ 100
   A500 -: ((+/ .*) (+ @ |:)) (potf2_pjlap_ A500)
1
   A1000=. ((+/ .*) (+ @ |:)) j./ ? 2 1000 1000 $ 100
   A1000 -: ((+/ .*) (+ @ |:)) (potf2_pjlap_ A1000)
1
   ts=: 6!:2, 7!:2@]

   A50d=. ((+/ .*) |:) ? 50 50 $ 100
   (potf2_pjlap_ (+/ @: (+/) @: | @: -) Choleski_pjlap_) A50d
6.33016e_8
   (A50d (+/ @: (+/) @: | @: -) ((+/ .*) (+ @ |:)) (potf2_pjlap_ A50d))
1.73168e_9
   (A50d (+/ @: (+/) @: | @: -) ((+/ .*) (+ @ |:)) (Choleski_pjlap_ A50d))
2.53472e_5

   A1000f=: ((+/ .*) |:) 10.1 + ? (2 $ 1000) $ 10
   A1000c=: ((+/ .*) (+ @ |:)) j./ ? (2 , (2 $ 1000)) $ 10
   ts & > 'z=. potf2l_pjlap_ A1000f';'z=. potrfl_pjlap_ A1000f';'z=. Cholesky_pjlap_ A1000f'
4.32201 1.78581e7
 1.0309 1.35258e7
1.76633 2.93628e7
   ts & > 'z=. potf2l_pjlap_ A1000c';'z=. potrfl_pjlap_ A1000c';'z=. Cholesky_pjlap_ A1000c'
7.73323 5.14207e7
2.44583 2.70426e7
3.96861 5.87229e7

)

NB. зависимость времени исполнения potf2l от размера матрицы

require 'plot'

mp=: +/ . *            NB. matrix product
h =: +@|:              NB. conjugate transpose
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
      tf=. (6!:2) 'z=. Cholesky Af'
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
