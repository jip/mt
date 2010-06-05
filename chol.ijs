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

mp=: +/ . *  NB. matrix product
h =: +@|:    NB. conjugate transpose

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

CholeskyNB=: 3 : 0
 n=.#y
 if. NBLOCK>:n do.
  potf2l y
 else.
  p=.>.-:n
  q=.n-p
  X=.(p,p){.y
  Y=.(p,-q){.y
  Z=.(-q,q){.y
  L0=.CholeskyNB X
  L2=. h L2h=.Y %. L0
  L1=.CholeskyNB Z - L2 mp L2h
  L0,L2,.L1
 end.
)

Cholesky1=: 3 : 0
 n=.#y
 if. 1>:n do.
NB.  assert. (A=|A)>0=A  NB. check for positive definite
  %:y
 else.
  p=.>.-:n
  q=.n-p
  X=.(p,p){.y
  Y=.(p,-q){.y
  Z=.(-q,q){.y
  L0=. Cholesky1 X
  L2=. h L2h=.Y %. L0
  L1=. Cholesky1 Z - L2 mp L2h
  L0,L2,.L1
 end.
)


Choleski=: 3 : 0
 n=.#A=.y
 if. 1>:n do.
  assert. (A=|A)>0=A  NB. check for positive definite
  %:A
 else.
  p=.>.n%2 [ q=.<.n%2
  X=.(p,p){.A [ Y=.(p,-q){.A [ Z=.(-q,q){.A
  L0=.Choleski X
  L1=.Choleski Z-(T=.(h Y) mp %.X) mp Y
  L0,(T mp L0),.L1
 end.
)

potrf=: (3 : 0) " 2
  n=. # y

NB. *           Compute the Cholesky factorization A = L*L'.
NB. *
NB.             DO 20 J = 1, N, NB
NB. *
NB. *              Update and factorize the current diagonal block and test
NB. *              for non-positive-definiteness.
NB. *
NB.                JB = MIN( NB, N-J+1 )
NB.                CALL ZHERK( 'Lower', 'No transpose', JB, J-1, -ONE,
NB.      $                     A( J, 1 ), LDA, ONE, A( J, J ), LDA )
NB.                CALL ZPOTF2( 'Lower', JB, A( J, J ), LDA, INFO )
NB.                IF( INFO.NE.0 )
NB.      $            GO TO 30
NB.                IF( J+JB.LE.N ) THEN
NB. *
NB. *                 Compute the current block column.
NB. *
NB.                   CALL ZGEMM( 'No transpose', 'Conjugate transpose',
NB.      $                        N-J-JB+1, JB, J-1, -CONE, A( J+JB, 1 ),
NB.      $                        LDA, A( J, 1 ), LDA, CONE, A( J+JB, J ),
NB.      $                        LDA )
NB.                   CALL ZTRSM( 'Right', 'Lower', 'Conjugate transpose',
NB.      $                        'Non-unit', N-J-JB+1, JB, CONE, A( J, J ),
NB.      $                        LDA, A( J+JB, J ), LDA )
NB.                END IF
NB.    20       CONTINUE

NB. CPU L1 cache size
    L1_CACHE=:  8*1024  NB.  8k for Intel Pentium 4
NB. L1_CACHE=: 64*1024  NB. 64k for AMD Athlon

NB. complex square matrix max size to fit entirely in L1 cache
NB. TODO: automate via sp=: 3 : '7!:5 <''y'''
NB. NBLOCK=: 22             NB. or 31 for floating point datatype

nb=: 2 ^ 4
nb2=: +: nb

rpotrf=: (3 : 0) " 2
  n=. #y
  if. n < nb2 do.
    potf2 y
    return.
  else.
    n1=. <. -: n
    j1=. >: n1
    'ios1 ios2'=. n1 ({. ; }.) i. n
    iosA11=. < 2 $ < ios1
    y=. iosA11 ((rpotrf @: {)`[`]) } y
    iosA21=. < ios2 ; ios1
    y=. iosA21 ((rgetrsu @: {)`[`]) } y
    iosA22=. < 2 $ < ios2
    y=. ((iosA21 { y) rsyrk (iosA22 { y)) iosA22 } y
    y=. iosA22 ((rpotrf @: {)`[`]) } y
  end.
)

rgetrsu=: (3 : 0) " 2 2
  [: ''
)

rsyrk=: (3 : 0) " 2 2
  [: ''
)

syrk=: (3 : 0) " 2 2
  [: ''
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
   5 (ts & >) 'z=. potf2_pjlap_ A100';'z=. potf2_pjlap_ A500';'z=. potf2_pjlap_ A1000'
0.0164126    810688
  1.08656 1.28652e7
  7.70693 5.14168e7

   A50d=. ((+/ .*) |:) ? 50 50 $ 100
   (potf2_pjlap_ (+/ @: (+/) @: | @: -) Choleski_pjlap_) A50d
6.33016e_8
   (A50d (+/ @: (+/) @: | @: -) ((+/ .*) (+ @ |:)) (potf2_pjlap_ A50d))
1.73168e_9
   (A50d (+/ @: (+/) @: | @: -) ((+/ .*) (+ @ |:)) (Choleski_pjlap_ A50d))
2.53472e_5


)

NB. зависимость времени исполнения potf2l от размера матрицы

require 'plot'
lio=: + ` (* i.)/ " 1                   NB. integers grid (2{y) steps from (0{y) by (1{y)
writetable=: 4 : '(toHOST,(": x),"1 LF) 1!:2 y'  NB. table writetable < 'filename'
readtable=: 3 : '>0 ". each cutopen toJ 1!:1 y'  NB. table=. readtable < 'filename'
flops_potf2=: 0 0.166667 1.41667 0.166667 & p.
flops_Chol1=: 3 : 0
  if. 2>y do.
    1
  else.
    (flops_Cho1 + (0 0 1 0.5 0.1875 & p.)) >. -: y
  end.
)

plotflopsN=: 3 : 0
  N=. 50
  vn=. lio 400 1 , N
  vFfp=. N $ 0
  vFfr=. N $ 0
  for_n. vn do.
    Af=: ((+/ .*) |:) ? (2 $ n) $ 100
NB.    A500c=: ((+/ .*) (+ @ |:)) j./ ? 2 500 500 $ 100
    'tfp tfr'=. ((6!:2) & >) 'z=. potf2l Af';'z=. Cholesky1 Af'
    if. 0 = (N%10) | n do.
      smoutput 'n_index' ; ((": >: n_index) , '/' , (": N)) ; 'n' ; n ; 'time float plain' ; tfp ; 'time float recursive' ; tfr
    end.
    vFfp=. (tfp %~ flops_potf2 n) n_index } vFfp
    vFfr=. (tfr %~ flops_Chol1 n) n_index } vFfr
  end.
  'type line,marker;keypos ctie;keystyle lbhc;title Cholesky factorization;ycaption FLOPs.;xcaption N;key potf2 Chol/1' plot vn ; vFfp ,: vFfr
  vn ; vFfp ; vFfr
)

NB. зависимость времени исполнения от размера квадрата, факторизуемого посредством potf2l

plotflopsNB=: 3 : 0
  N=. 150
  vn=. lio 1 1 , N
  mFfp=. (2 $ N) $ 0
  mFfr=. (2 $ N) $ 0
  for_n. vn do.
    Af=: ((+/ .*) |:) ? (2 $ n) $ 100
NB.    A500c=: ((+/ .*) (+ @ |:)) j./ ? 2 500 500 $ 100
    for_b. i. n do.
      NBLOCK=: >: b
      'tfp tfr'=. ((6!:2) & >) 'z=. potf2l Af';'z=. CholeskyNB Af'
      mFfp=. tfp (< n_index , b) } mFfp
      mFfr=. tfr (< n_index , b) } mFfr
    end.
    if. 0 = (N%10) | n do.
      smoutput 'n_index' ; ((": >: n_index) , '/' , (": N)) ; 'n' ; n
    end.
  end.
smoutput 'vn' ; ($vn) ; 'mFfp' ; ($mFfp) ; 'mFfr' ; ($mFfp)
  pd 'reset'
  pd 'title Cholesky factorization blocked recursion exec time, sec'
  pd 'xcaption N'
  pd 'ycaption NB'
  pd 'type surface'
  pd vn;vn;mFfp
  pd 'type surface'
  pd vn;vn;mFfr
  pd 'show'
  pd 'save bmp chol_nb.bmp'
  mFfp writetable < 'mFfp.dump'
  mFfr writetable < 'mFfr.dump'
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
