NB. qf.ijs
NB. Orthogonal factorizations LQ QL QR RQ
NB.
NB. gelq2  LQ factorization of a general matrix (non-blocked
NB.        version)
NB. geql2  QL factorization of a general matrix (non-blocked
NB.        version)
NB. geqr2  QR factorization of a general matrix (non-blocked
NB.        version)
NB. gerq2  RQ factorization of a general matrix (non-blocked
NB.        version)
NB.
NB. gelqf  LQ factorization of a general matrix
NB. geqlf  QL factorization of a general matrix
NB. geqrf  QR factorization of a general matrix
NB. gerqf  RQ factorization of a general matrix
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. Nouns, differences between cIOSs at consequent
NB. iterations: cios(i+1)-cios(i)

GELQ2DCIOS=: 3 2 $ 1 0j_1 1 0 0j_1 0j_1    NB. z,τ,L
GEQL2DCIOS=: 3 2 $ 0j_1 _1 0 _1 0j_1 0j_1  NB. z,τ,L
GEQR2DCIOS=: 3 2 $ 0j_1 1 0 1 0j_1 0j_1    NB. z,τ,R
GERQ2DCIOS=: 3 2 $ _1 0j_1 _1 0 0j_1 0j_1  NB. z,τ,R

NB. ---------------------------------------------------------
NB. mkcios0gelq2
NB. mkcios0geql2
NB. mkcios0geqr2
NB. mkcios0gerq2
NB.
NB. Create cIOS at 0-th iteration for corresp. method
NB.
NB. Syntax:
NB.   cios0=. mkcios0gelq2 mn
NB.   cios0=. mkcios0geql2 mn
NB.   cios0=. mkcios0geqr2 mn
NB.   cios0=. mkcios0gerq2 mn
NB. where
NB.   mn    - 2-vector of integers (m,n), shape of matrix to
NB.           factorize
NB.   cios0 - 3×2-table cios(0), cIOSs corresponding to
NB.           iteration 0, see ge*2step verbs

mkcios0gelq2=: 3 : 0
  'm1 n1'=. _1 1 + 'm n'=. y
  3 2 $ 0 _1 0 _1 _1 _2 j. (1 , n1 , 1 1 , m1 , n)
)

mkcios0geql2=: 3 : 0
  'm1 n1'=. 1 _1 + 'm n'=. y
  3 2 $ 0 _1 0 _1 1 0 j. (m1 , 1 1 1 , m , n1)
)

mkcios0geqr2=: 3 : 0
  'm1 n1'=. 1 _1 + 'm n'=. y
  3 2 $ _1 0 _1 0 _2 _1 j. (m1 , 1 1 1 , m , n1)
)

mkcios0gerq2=: 3 : 0
  'm1 n1'=. _1 1 + 'm n'=. y
  3 2 $ _1 0 _1 0 0 1 j. (1 , n1 , 1 1 , m1 , n)
)

NB. ---------------------------------------------------------
NB. geq2step
NB.
NB. Template adv. to form verbs ge*2step
NB.
NB. Syntax:
NB.   vstep=. vref`vapp geq2step
NB. where
NB.   vref  - verb to generate an elementary reflector, see
NB.           larfg* larfp*; is called as:
NB.             z=. vref y
NB.   vapp  - verb to apply an elementary reflector, see
NB.           larfL* larfR*; is called as:
NB.             Aupd=. cios vapp A
NB.   vstep - verb to perform single step of Q-factorization;
NB.           see ge*2step; is called as:
NB.             'Ai1 ciosi1'=. dcios gelq2step (Ai ; ciosi)

geq2step=: 1 : '(< @ (+ (1&{::))) 1} (((1 {:: ]) (< @ (m gerf0 0 1)) (0 {:: ])) 0} ])'

NB. ---------------------------------------------------------
NB. gelq2step
NB. geql2step
NB. geqr2step
NB. gerq2step
NB.
NB. Single step of corresp. method
NB.
NB. Syntax:
NB.   'Ai1 ciosi1'=. dcios gelq2step (Ai ; ciosi)
NB.   'Ai1 ciosi1'=. dcios geql2step (Ai ; ciosi)
NB.   'Ai1 ciosi1'=. dcios geqr2step (Ai ; ciosi)
NB.   'Ai1 ciosi1'=. dcios gerq2step (Ai ; ciosi)
NB. where
NB.   Ai    - (m+1)×n-matrix A(i) to update before i-th
NB.           iteration (i=0:min(m,n)-1)
NB.   ciosi - 3×2-matrix cios(i) of cIOSs (see struct.ijs)
NB.           for i-th iteration; rows (0:2) contains:
NB.             0 - (m-i+1)-vector Y=(α[i],x[i][1:m-(i+1)],0)
NB.                 to reflect, or vector
NB.                 Z=(β[i],v[i][1:m-(i+1)],τ[i])
NB.                 of reflection result, is stored in
NB.                 A[i:m,i], cIOS are:
NB.                   (_1 j. (m-i+1)) , (i j. 1)
NB.             1 - scalar τ[i], is stored in A[m,i], cIOS
NB.                 is:
NB.                   (_1 j. 1) , (i j. 1)
NB.             2 - (m-i)×(n-(i+1))-matrix L to apply an
NB.                 elementary reflector from the left, is
NB.                 stored in A[i:m-1,i+1:n-1], cIOS are:
NB.                   (_2 j. (m-i)) , (_1 j. (n-(i+1)))
NB.   dcios - difference between cIOSs at consequent
NB.           iterations: cios(i+1)-cios(i)
NB.   Ai1   - (m+1)×n-matrix A(i+1) after i-th iteration
NB.   ciosi1 - 3×2-matrix cios(i+1) of cIOSs for (i+1)-th
NB.           iteration

NB.           (), performs actions:
NB.             1) extract A(i) and cios(i) and supply its to
NB.                gerf0; the last is configured to use
NB.                gerund (vref`vapp), to get vectors y(i)
NB.                and z(i)'s cIOS from cios(i)[0], scalar
NB.                τ(i)'s cIOS from cios(i)[1];
NB.             2) box output A(i+1) and write it into 0-th
NB.                item of input;
NB.             3) adjust cios(i) by Δcios;
NB.             4) box output cios(i+1) and write it into
NB.                1-th item of input.

gelq2step=: (larfgfc`(0 2 larfRfcs)) geq2step
geql2step=: (larfgb`(0 2 larfLbsc)) geq2step
geqr2step=: (larfgf`(0 2 larfLfsc)) geq2step
gerq2step=: (larfgbc`(0 2 larfRbcs)) geq2step

NB. ---------------------------------------------------------
NB. Template adv. to form verbs ge*2

geq2=: 1 : '(<./ @ $) ((0&{::)@((2{m)`:6)) ((0 ((0{m)`:6) ]) ; (((1{m)`:6) @ $))'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gelq2
NB. geql2
NB. geqr2
NB. gerq2
NB. emulate xGELQ2
NB. LQf=. gelq2 A

gelq2=: ,.~`mkcios0gelq2`(GELQ2DCIOS & gelq2step) geq2
geql2=: ,  `mkcios0geql2`(GEQL2DCIOS & geql2step) geq2
geqr2=: ,~ `mkcios0geqr2`(GEQR2DCIOS & geqr2step) geq2
gerq2=: ,. `mkcios0gerq2`(GERQ2DCIOS & gerq2step) geq2

NB. ---------------------------------------------------------
NB. gelqf
NB. geqlf
NB. geqrf
NB. gerqf
NB. emulate xGELQF

gelqf=: gelq2  NB. stub for a while
geqlf=: geql2  NB. stub for a while

QFBS=: 3  NB. block size for qf algorithms

GEQRFDCIOS=: QFBS * 3 2 $ 1j_1 1 0 0 1j_1 1j_1  NB. Z,T,R

NB. cios0=. mkcios0geqrf m,n,bs
mkcios0geqrf=: 3 : 0
  'm n bs'=. y
  3 2 $ (0 0,(m+1),0 0 0) j. (m,bs,bs,bs,m,(n-nb))
)

geqrfstep=: [ ([ (larfb map3 0 1 2) (larftfc map 0 1)) (geqr2 upd 0)

geqrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  i=. <. k % QFBS      NB. # of iterations
  y=. (m+1+QFBS) {. y  NB. pre-allocate space for τ and T, filled by zeros
  cios0=. mkcios0geqrf m,n,QFBS
  (m+1) {. 0 {:: i (GEQRFDCIOS & ungl2step) (y ; cios0)
)

gerqf=: gerq2  NB. stub for a while

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB.*tgeqf v test gelq2 geql2 geqr2 gerq2 gelqf geqlf geqrf gerqf

tgeqf=: 3 : 0
)

NB. ---------------------------------------------------------
NB.*testqf v test orthogonal factorizations

testqf=: 3 : 0
)

NB. Af6x6=. 6 6 $ 6 _4 9 6 3 4 8 1 0 9 5 _9 _7 9 6 3 5 _4 8 9 _7 _9 _2 _9 3 _7 _1 7 3 _8 _8 9 8 8 _6 6
NB. Af6x4=. 6 4 $ 4 5 4 _7 _3 6 5 0 9 _6 _4 9 _9 0 3 _3 _1 _1 3 2 _2 3 _6 _5
NB. Af4x6=. 4 6 $ 5 7 _7 6 _6 _3 _5 9 _5 _5 _5 6 0 1 6 8 2 2 4 _7 _5 _3 _1 1
NB. Ac6x6=. 6 6 $ 8j_7 _8j_8 1j_6 _4j1 4j_2 _7j_2 _4j9 _8j_2 _1j_4 1j_7 5j_8 6j_9 5j_9 1j_5 _9j_7 _6j7 _1j_5 1j8 5j5 5j6 _3j_7 2j2 _8j_7 9j8 _2j6 8j_1 _6j7 _5j_3 _5j_6 _8 1j7 _7j8 6j2 3 8j_5 _1
NB. Ac6x4=. 6 4 $ _3j1 _5j_2 1 _9j6 8j_4 _2j2 _8j_1 _6j6 _2j5 _8j5 _5j_1 3j3 _7j_6 _7j_8 _3j_4 4j_2 5j4 8j8 6 3j6 _2j_6 9j5 _9 _1j1
NB. Ac4x6=. 4 6 $ _3j_5 2j2 _1j_4 _5j_8 _4j8 8j_5 8j_2 2j4 _5j_3 _6j_5 8j4 9j_1 7j2 _6j_5 7j_1 _5j6 _1 _1j9 1j_5 _4j_7 _6 _1j_9 _6j1 _7j1
NB. smoutput 'Qf' ; ($ Qf) ; Qf ; 'Tau' ; ($ Tau) ; Tau
