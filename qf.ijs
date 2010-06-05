NB. qf.ijs
NB. Orthogonal factorizations LQ QL QR RQ
NB.
NB. gelq2  LQ factorization of a matrix (non-blocked version)
NB. gelqf  LQ factorization of a matrix
NB. geql2  QL factorization of a matrix (non-blocked version)
NB. geqlf  QL factorization of a matrix
NB. geqr2  QR factorization of a matrix (non-blocked version)
NB. geqrf  QR factorization of a matrix
NB. gerq2  RQ factorization of a matrix (non-blocked version)
NB. gerqf  RQ factorization of a matrix
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. Nouns, differences between cIOSs at consequent
NB. iterations: cios(i+1)-cios(i)

GELQ2DCIOS=: 4 2 $ 1 0j_1 1 0j_1 1 0 0j_1 0j_1
GEQL2DCIOS=: 4 2 $ 0j1 _1 0j1 _1 0 _1 0j_1 0j_1
GEQR2DCIOS=: 4 2 $ 0j_1 1 0j_1 1 0 1 0j_1 0j_1
GERQ2DCIOS=: 4 2 $ _1 0j_1 _1 0j_1 _1 0 0j_1 0j_1

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
NB.   cios0 - 4×2-table cios(0), cIOSs corresponding to
NB.           iteration 0, see ge*2step verbs

mkcios0gelq2=: 3 : 0
  'm1 n1'=. _1 1 + 'm n'=. y
  4 2 $ 0 _2 0 _1 0 _1 _1 _2 j. (1 , n , 1 , n1 , 1 1 , m1 , n)
)

mkcios0geql2=: 3 : 0
  'm1 n1'=. _1 1 * 1 _1 + 'm n'=. y
  nm=. - m
  4 2 $ 1 _1 0 _1 0 _1 1 0 j. (nm , 1 , m1 , 1 1 1 , m , n1)
)

mkcios0geqr2=: 3 : 0
  'm1 n1'=. 1 _1 + 'm n'=. y
  4 2 $ _2 0 _1 0 _1 0 _2 _1 j. (m , 1 , m1 , 1 1 1 , m , n1)
)

mkcios0gerq2=: 3 : 0
  'm1 n1'=. 1 _1 * _1 1 + 'm n'=. y
  nn=. - n
  4 2 $ _1 1 _1 0 _1 0 0 1 j. (1 , nn , 1 , n1 , 1 1 , m1 , n)
)

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
NB.   ciosi - 5×2-matrix cios(i) of cIOSs (see struct.ijs)
NB.           for i-th iteration; rows (0:4) contains:
NB.             0 - (m-i)-vector Y=(α[i],x[i][1:m-(i+1)])
NB.                 to reflect, is stored in A[i:m-1,i],
NB.                 cIOS are:
NB.                   (_2 j. (m-i)) , (i j. 1)
NB.             1 - (m-i+1)-vector Z=(β[i],v[i][1:m-(i+1)],τ[i])
NB.                 of reflection result, is stored in
NB.                 A[i:m,i], cIOS are:
NB.                   (_1 j. (m-i+1)) , (i j. 1)
NB.             2 - scalar τ[i], is stored in A[m,i], cIOS
NB.                 is:
NB.                   (_1 j. 1) , (i j. 1)
NB.             3 - (m-i)×(n-(i+1))-matrix L to apply an
NB.                 elementary reflector from the left, is
NB.                 stored in A[i:m-1,i+1:n-1], cIOS are:
NB.                   (_2 j. (m-i)) , (_1 j. (n-(i+1)))
NB.   dcios - difference between cIOSs at consequent
NB.           iterations: cios(i+1)-cios(i)
NB.   Ai1   - (m+1)×n-matrix A(i+1) after i-th iteration
NB.   ciosi - 5×2-matrix cios(i+1) of cIOSs for (i+1)-th
NB.           iteration

NB. FIXME! pre_conj(y), post_conj(z)
gelq2step=: (< @ (+ (1&{::))) 1} (((1 {:: ]) (< @ ((larfg`(1 3 larfR)) gerf0 0 1 2)) (0 {:: ])) 0} ])

NB. CHECKME!
geql2step=: (< @ (+ (1&{::))) 1} (((1 {:: ]) (< @ (((larfg`(1 3 larfL)) gerf0 0 1 2) dbg 'gerf0')) (0 {:: ])) 0} ])

NB. OK
geqr2step=: (< @ (+ (1&{::))) 1} (((1 {:: ]) (< @ ((larfg`(1 3 larfL)) gerf0 0 1 2)) (0 {:: ])) 0} ])

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. unghr
NB. Generate an unitary matrix Q which is defined as the
NB. product of ({:fs) elementary reflectors of order n, as
NB. returned by gehrd:
NB.    Q = H(f) H(f+1) . . . H(f+s-1)
NB.
NB. Q=. unghr A ; fs
NB. fs - (f,s), defines where Qf is in A

unghr=: 3 : 0
  cfrom=. ({~ cios2ios)~
  'A fs'=. y
  e=. +/ fs
  fjs=. j./ fs
NB. smoutput '2 $ fjs' ; ($ 2 $ fjs) ; (2 $ fjs) ; 'A' ; ($ A) ; A ; '(2 $ fjs) cfrom A' ; ($ (2 $ fjs) cfrom A) ; ((2 $ fjs) cfrom A)
  Qf=. trl1 (2 $ fjs) cfrom A
  Tau=. {. ((e j. 1) , fjs) cfrom A
NB. smoutput 'Qf' ; ($ Qf) ; Qf ; 'Tau' ; ($ Tau) ; Tau
  mp/ (idmat # Qf) -"2 Tau * (* +)"0/~"1 |:Qf
)

NB. ---------------------------------------------------------
NB. gelq2
NB. emulate xGELQ2
NB. LQf=. gelq2 A

gelq2=: 3 : 0
  k=. <./ mn=. $ y
  y=. y ,. 0    NB. append zero column to A to store τ[0:min(m,n)-1]

  NB. link A and cios(0), do iterations, then extract LQf
  0 {:: k (GELQ2DCIOS & gelq2step) (y ; (mkcios0gelq2 mn))
)

gelqf=: gelq2  NB. stub for a while

NB. ---------------------------------------------------------
NB. geql2
NB. emulate xGEQL2
NB. QfL=. geql2 A
NB. FIXME!

geql2=: 3 : 0
  k=. <./ mn=. $ y
  y=. 0 , y    NB. prepend zero row to A to store τ[n-min(m,n):n-1]

  NB. link A and cios(0), do iterations, then extract QfL
  0 {:: k (GEQL2DCIOS & geql2step) (y ; (mkcios0geql2 mn))
)

geqlf=: geql2  NB. stub for a while

NB. ---------------------------------------------------------
NB. geqr2
NB. emulate xGEQR2
NB. RQf=. geqr2 A

geqr2=: 3 : 0
  k=. <./ mn=. $ y
  y=. y , 0    NB. append in-place zero row to A to store τ[0:min(m,n)-1]

  NB. link A and cios(0), do iterations, then extract RQf
  0 {:: k (GEQR2DCIOS & geqr2step) (y ; (mkcios0geqr2 mn))
)

geqrf=: geqr2  NB. stub for a while

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB.*tgeqrf v test geqrf

tgeqrf=: 3 : 0
)

NB. ---------------------------------------------------------
NB.*testorf v test orthogonal factorizations

testorf=: 3 : 0
)

NB. Af6x6=. 6 6 $ 6 _4 9 6 3 4 8 1 0 9 5 _9 _7 9 6 3 5 _4 8 9 _7 _9 _2 _9 3 _7 _1 7 3 _8 _8 9 8 8 _6 6
NB. Af6x4=. 6 4 $ 4 5 4 _7 _3 6 5 0 9 _6 _4 9 _9 0 3 _3 _1 _1 3 2 _2 3 _6 _5
NB. Af4x6=. 4 6 $ 5 7 _7 6 _6 _3 _5 9 _5 _5 _5 6 0 1 6 8 2 2 4 _7 _5 _3 _1 1
NB. Ac6x6=. 6 6 $ 8j_7 _8j_8 1j_6 _4j1 4j_2 _7j_2 _4j9 _8j_2 _1j_4 1j_7 5j_8 6j_9 5j_9 1j_5 _9j_7 _6j7 _1j_5 1j8 5j5 5j6 _3j_7 2j2 _8j_7 9j8 _2j6 8j_1 _6j7 _5j_3 _5j_6 _8 1j7 _7j8 6j2 3 8j_5 _1
NB. Ac6x4=. 6 4 $ _3j1 _5j_2 1 _9j6 8j_4 _2j2 _8j_1 _6j6 _2j5 _8j5 _5j_1 3j3 _7j_6 _7j_8 _3j_4 4j_2 5j4 8j8 6 3j6 _2j_6 9j5 _9 _1j1
NB. Ac4x6=. 4 6 $ _3j_5 2j2 _1j_4 _5j_8 _4j8 8j_5 8j_2 2j4 _5j_3 _6j_5 8j4 9j_1 7j2 _6j_5 7j_1 _5j6 _1 _1j9 1j_5 _4j_7 _6 _1j_9 _6j1 _7j1
NB. ((tru @ (0 & {::)) ; unghr) (geqr2 Af6x6) ; 0 6
