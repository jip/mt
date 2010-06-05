NB. gq.ijs
NB. Generate Q from LQ QL QR RQ output
NB.
NB. unglq  Generate a matrix with orthonormal rows from
NB.        output of gelq2 or gelqf
NB. ungql  Generate a matrix with orthonormal columns from
NB.        output of geql2 or geqlf
NB. ungqr  Generate a matrix with orthonormal columns from
NB.        output of geqr2 or geqrf
NB. ungrq  Generate a matrix with orthonormal rows from
NB.        output of gerq2 or gerqf
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. Nouns, differences between cIOSs at consequent
NB. iterations: cios(i+1)-cios(i)

UNGL2DCIOS=: 6 2 $ _1 0j1 0j1 0j1 _1 0j1 _1 0 _1 _1 _1 0j_1  NB. z,R,vtoscale,τ,diag,rto0
UNG2LDCIOS=: 6 2 $ 0j1 1 0j1 0j1 0j1 1 0 1 1 1 0j_1 1        NB. z,L,vtoscale,τ,diag,cto0
UNG2RDCIOS=: 6 2 $ 0j1 _1 0j1 0j1 0j1 _1 0 _1 _1 _1 0j_1 _1  NB. z,L,vtoscale,τ,diag,cto0
UNGR2DCIOS=: 6 2 $ 1 0j1 0j1 0j1 1 0j1 1 0 1 1 1 0j_1        NB. z,R,vtoscale,τ,diag,rto0

NB. ---------------------------------------------------------
NB. mkcios0ungl2
NB. mkcios0geql2
NB. mkcios0geqr2
NB. mkcios0gerq2
NB.
NB. Create cIOS at 0-th iteration for corresp. method
NB.
NB. Syntax:
NB.   cios0=. mkcios0ungl2 (m,n,k,p)
NB.   cios0=. mkcios0geql2 (m,n,k,p)
NB.   cios0=. mkcios0geqr2 (m,n,k,p)
NB.   cios0=. mkcios0gerq2 (m,n,k,p)
NB. where
NB.   kn    - 2-vector of integers (k,n), shape of matrix to
NB.           generate
NB.   cios0 - 6×2-table cios(0), cIOSs corresponding to
NB.           iteration 0, see ung*step verbs

mkcios0ungl2=: 3 : 0
  p1=. (<: @ {:) 'm n k p'=. y
  6 2 $ (p1,_1,((k-1)-m),_2,p1,_2,p1,_1,p1,p1,p1,0) j. (1,(n+2-p),(k-p),(n-p1),1,(n-p),1 1 1 1 1,p1)
)

mkcios0ung2l=: 3 : 0
  mp1=. >: {. 'mp np kp p1'=. (}: (- , (<: @ ])) {:) 'm n k p'=. y
  6 2 $ (0,np,1,(n-k),1,np,0,np,mp1,np,_1,np) j. ((m+2-p),1,mp1,kp,mp,1 1 1 1 1,p1,1)
)

mkcios0ung2r=: 3 : 0
  p1=. <: {: 'm n k p'=. y
  6 2 $ (_1,p1,_2,((k-1)-n),_2,p1,_1,p1,p1,p1,0,p1) j. ((m+2-p),1,(m-p1),(k-p),(m-p),1 1 1 1 1,p1,1)
)

mkcios0ungr2=: 3 : 0
  'mp np1'=. 0 1 + ((2 & {.) - {:) 'm n k p'=. y
  6 2 $ (mp,0,(m-k),1,mp,1,mp,0,mp,np1,mp,_1) j. (1,(n+2-p),(k-p),np1,1,(n-p),1 1 1 1 1,(<:p))
)

NB. ---------------------------------------------------------
NB. Template adv. to form verbs ung((2[lr])|([lr]2))step

ungq2step=: 1 : '([ ([ ([ (0: amendci 5) ((1&((0{m)`:6)) mapci 3 4)) (((* (((1{m)`:6)@(0&({,))))~) upd2ci 3 2)) ((2{m)`:6)) step'

ungl2step=: (- +)`(+@-)`(0 1 larfRfcc) ungq2step
ung2lstep=:  -   `   - `(0 1 larfLbss) ungq2step
ung2rstep=:  -   `   - `(0 1 larfLfss) ungq2step
ungr2step=: (- +)`(+@-)`(0 1 larfRbcc) ungq2step

NB. ---------------------------------------------------------
NB. ungl2
NB. Generate a matrix with orthonormal rows from output of
NB. gelq2 or gelqf (non-blocked version)
NB.
NB. Syntax:
NB.   Q=. [p] ungl2 LQf
NB. where
NB.   LQf - m×(n+1)-matrix, output of gelq2 or gelqf
NB.   p   - integer in range 0:k, default is k, the number of
NB.         elementary reflectors, whose product defines the
NB.         matrix Q
NB.   Q   - k×n-matrix Q with orthonormal rows, which is
NB.         defined as the first k rows of a product of p
NB.         elementary reflectors of order n:
NB.           Q = H(p)' ... H(2)' H(1)'
NB.   k   = min(m,n)
NB.
NB. If:
NB.   LQf=. gelq2 A
NB.   L=. trl 0 _1 }. LQf
NB.   Q=. ungl2 LQf
NB. then
NB.   I -: (mp ct) Q
NB.   A -: L mp Q
NB.   (-: (((trl @ (0 _1 & }.)) mp ungl2) @ gelq2)) A

ungl2=: ($:~ (<./ @ (0 _1 & +) @ $)) :(4 : 0)
  k=. <./ 'm n'=. 0 _1 + $ y
  cios0=. mkcios0ungl2 (m , n , k , x)
  sizeI=. (k - x) , n
  ciosI=. (x , 0) j. sizeI
  y=. (x idmat sizeI) (cios2ios ciosI) } y
  (k , n) {. 0 {:: x (UNGL2DCIOS & ungl2step) (y ; cios0)
)

NB. ---------------------------------------------------------
NB. ung2l
NB. Generate a matrix with orthonormal columns from output of
NB. geql2 or geqlf (non-blocked version)
NB.
NB. Syntax:
NB.   Q=. [p] ung2l QfL
NB. where
NB.   QfL - (m+1)×n-matrix, output of geql2 or geqlf
NB.   p   - integer in range 0:k, default is k, the number of
NB.         elementary reflectors, whose product defines the
NB.         matrix Q
NB.   Q   - m×k-matrix Q with orthonormal columns, which is
NB.         defined as the last k columns of a product of p
NB.         elementary reflectors of order m:
NB.           Q = H(p) ... H(2) H(1)
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'm n'=. $ A
NB.   QfL=. geql2 A
NB.   Q=. ung2l QfL
NB.   L=. (n - m) trl }. QfL
NB. then
NB.   I -: (mp~ ct) Q
NB.   A -: Q mp L
NB.   (-: ((ung2l mp (((n - m) & trl) @ }.)) @ geql2)) A

ung2l=: ($:~ (<./ @ (_1 0 & +) @ $)) :(4 : 0)
  k=. <./ 'm n'=. _1 0 + $ y
  cios0=. mkcios0ung2l (m , n , k , x)
  sizeI=. m , (k - x)
  ciosI=. (1 , (n - k)) j. sizeI
  y=. ((k - m) idmat sizeI) (cios2ios ciosI) } y
  (- (m , k)) {. 0 {:: x (UNG2LDCIOS & ung2lstep) (y ; cios0)
)

NB. ---------------------------------------------------------
NB. ung2r
NB. Generate a matrix with orthonormal columns from output of
NB. geqr2 or geqrf (non-blocked version)
NB.
NB. Syntax:
NB.   Q=. [p] ung2r QfR
NB. where
NB.   QfR - (m+1)×n-matrix, output of geqr2 or geqrf
NB.   p   - integer in range 0:k, default is k, the number of
NB.         elementary reflectors, whose product defines the
NB.         matrix Q
NB.   Q   - m×k-matrix Q with orthonormal columns, which is
NB.         defined as the first k columns of a product of p
NB.         elementary reflectors of order m:
NB.           Q = H(1) H(2) ... H(p)
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'm n'=. $ A
NB.   QfR=. geqr2 A
NB.   Q=. ung2r QfR
NB.   R=. tru }: QfR
NB. then
NB.   I -: (mp~ ct) Q
NB.   A -: Q mp R
NB.   (-: ((ung2r mp (tru @ }:)) @ geqr2)) A

ung2r=: ($:~ (<./ @ (_1 0 & +) @ $)) :(4 : 0)
  k=. <./ 'm n'=. _1 0 + $ y
  cios0=. mkcios0ung2r (m , n , k , x)
  sizeI=. m , (k - x)
  ciosI=. (0 , x) j. sizeI
  y=. ((- x) idmat sizeI) (cios2ios ciosI) } y
  (m , k) {. 0 {:: x (UNG2RDCIOS & ung2rstep) (y ; cios0)
)

NB. ---------------------------------------------------------
NB. ungr2
NB. Generate a matrix with orthonormal rows from output of
NB. gerq2 or gerqf (non-blocked version)
NB.
NB. Syntax:
NB.   Q=. [p] ungr2 RQf
NB. where
NB.   RQf - m×(n+1)-matrix, output of gerq2 or gerqf
NB.   p   - integer in range 0:k, default is k, the number of
NB.         elementary reflectors, whose product defines the
NB.         matrix Q
NB.   Q   - k×n-matrix Q with orthonormal rows, which is
NB.         defined as the last k rows of a product of p
NB.         elementary reflectors of order n:
NB.           Q = H(1)' H(2)' ... H(p)'
NB.   k   = min(m,n)
NB.
NB. If:
NB.   'm n'=. $ A
NB.   RQf=. gerq2 A
NB.   R=. (n - m) trl 0 1 }. RQf
NB.   Q=. ungr2 RQf
NB. then
NB.   I -: (mp ct) Q
NB.   A -: R mp Q
NB.   (-: (((((n - m) & tru) @ (0 1 & }.)) mp ungr2) @ gerq2)) A

ungr2=: ($:~ (<./ @ (0 _1 & +) @ $)) :(4 : 0)
  k=. <./ 'm n'=. 0 _1 + $ y
  cios0=. mkcios0ungr2 (m , n , k , x)
  sizeI=. (k - x) , n
  ciosI=. ((m - k) , 1) j. sizeI
  y=. (idmat sizeI) (cios2ios ciosI) } y
  (- (k , n)) {. 0 {:: x (UNGR2DCIOS & ungr2step) (y ; cios0)
)

NB. =========================================================
NB. Interface

NB. TODO:
NB. - optimize write-in idmat: zeros, then unit diag
NB. - template adv. ung2

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB.*tgq v test ung*

tgq=: 3 : 0
)

NB. ---------------------------------------------------------
NB.*testgq v test Q genberators

testgq=: 3 : 0
  Ac6x6=. 6 6 $ 8j_7 _8j_8 1j_6 _4j1 4j_2 _7j_2 _4j9 _8j_2 _1j_4 1j_7 5j_8 6j_9 5j_9 1j_5 _9j_7 _6j7 _1j_5 1j8 5j5 5j6 _3j_7 2j2 _8j_7 9j8 _2j6 8j_1 _6j7 _5j_3 _5j_6 _8 1j7 _7j8 6j2 3 8j_5 _1
  Ac6x4=. 6 4 $ _3j1 _5j_2 1 _9j6 8j_4 _2j2 _8j_1 _6j6 _2j5 _8j5 _5j_1 3j3 _7j_6 _7j_8 _3j_4 4j_2 5j4 8j8 6 3j6 _2j_6 9j5 _9 _1j1
  Ac4x6=. 4 6 $ _3j_5 2j2 _1j_4 _5j_8 _4j8 8j_5 8j_2 2j4 _5j_3 _6j_5 8j4 9j_1 7j2 _6j_5 7j_1 _5j6 _1 _1j9 1j_5 _4j_7 _6 _1j_9 _6j1 _7j1

  smoutput 'LQ 6x6' ; (; (((trl @ (0 _1 & }.)) mp ungl2) @ gelq2)) Ac6x6
  smoutput 'LQ 6x4' ; (; (((trl @ (0 _1 & }.)) mp ungl2) @ gelq2)) Ac6x4
  smoutput 'LQ 4x6' ; (; (((trl @ (0 _1 & }.)) mp ungl2) @ gelq2)) Ac4x6

  smoutput 'QL 6x6' ; (; ((ung2lx mp (((6 - 6) & trl) @ }.)) @ geql2)) Ac6x6
  smoutput 'QL 6x4' ; (; ((ung2lx mp (((4 - 6) & trl) @ }.)) @ geql2)) Ac6x4
  smoutput 'QL 4x6' ; (; ((ung2lx mp (((6 - 4) & trl) @ }.)) @ geql2)) Ac4x6

  smoutput 'QR 6x6' ; (; ((ung2rx mp (tru @ }:)) @ geqr2)) Ac6x6
  smoutput 'QR 6x4' ; (; ((ung2rx mp (tru @ }:)) @ geqr2)) Ac6x4
  smoutput 'QR 4x6' ; (; ((ung2rx mp (tru @ }:)) @ geqr2)) Ac4x6

  smoutput 'RQ 6x6' ; (; (((((6 - 6) & tru) @ (0 1 & }.)) mp ungr2) @ gerq2)) Ac6x6
  smoutput 'RQ 6x4' ; (; (((((4 - 6) & tru) @ (0 1 & }.)) mp ungr2) @ gerq2)) Ac6x4
  smoutput 'RQ 4x6' ; (; (((((6 - 4) & tru) @ (0 1 & }.)) mp ungr2) @ gerq2)) Ac4x6
)
