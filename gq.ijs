NB. Generate matrix with orthonormal rows or columns from its
NB. factored form
NB.
NB. unglq      Generate a matrix with orthonormal rows from
NB.            output of gelqf
NB. ungql      Generate a matrix with orthonormal columns
NB.            from output of geqlf
NB. ungqr      Generate a matrix with orthonormal columns
NB.            from output of geqrf
NB. ungrq      Generate a matrix with orthonormal rows from
NB.            output of gerqf
NB. unglz      Generate a matrix with orthonormal rows from
NB.            output of tzlzf
NB. ungzl      Generate a matrix with orthonormal columns
NB.            from output of tzzlf
NB. ungzr      Generate a matrix with orthonormal columns
NB.            from output of tzzrf
NB. ungrz      Generate a matrix with orthonormal rows from
NB.            output of tzrzf
NB. unghrx     Generate an unitary (orthogonal) matrix which
NB.            is defined as the product of elementary
NB.            reflectors as returned by gehrdx
NB.
NB. testungq   Test ungxx by general matrix
NB. testungz   Test ungxx by trapezoidal matrix
NB. testunghr  Test unghrx by square matrix
NB. testgq     Adv. to make verb to test ungxxx by matrix of
NB.            generator and shape given
NB.
NB. Version: 0.11.0 2021-01-17
NB.
NB. Copyright 2010-2021 Igor Zhuravlov
NB.
NB. This file is part of mt
NB.
NB. mt is free software: you can redistribute it and/or
NB. modify it under the terms of the GNU Lesser General
NB. Public License as published by the Free Software
NB. Foundation, either version 3 of the License, or (at your
NB. option) any later version.
NB.
NB. mt is distributed in the hope that it will be useful, but
NB. WITHOUT ANY WARRANTY; without even the implied warranty
NB. of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
NB. See the GNU Lesser General Public License for more
NB. details.
NB.
NB. You should have received a copy of the GNU Lesser General
NB. Public License along with mt. If not, see
NB. <http://www.gnu.org/licenses/>.

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

GQNB=: 32   NB. block size limit
GQNX=: 128  NB. crossover point, GQNX ≥ GQNB

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. unglq
NB.
NB. Description:
NB.   Generate a matrix with orthonormal rows from output of
NB.   gelqf
NB.
NB. Syntax:
NB.   Q=. [k] unglq LQf
NB. where
NB.   LQf - m×(n+1)-matrix, the output of gelqf
NB.   Q   - k×n-matrix with orthonormal rows, which is
NB.         defined as the first k rows of the product of k
NB.         elementary reflectors H(i) of order n:
NB.           Q = Π{H(i)',i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   k   ≤ n, optional, default is min(m,n)
NB.
NB. Storage layout:
NB.   LQf for m=3, n=7:                  LQf for m=7, n=3:
NB.     (  l  v0 v0 v0 v0 v0 v0 τ0  )      (  l  v0 v0 τ0  )
NB.     (  l  l  v1 v1 v1 v1 v1 τ1  )      (  l  l  v1 τ1  )
NB.     (  l  l  l  v2 v2 v2 v2 τ2  )      (  l  l  l  τ2  )
NB.                                        (  l  l  l  *   )
NB.                                        (  l  l  l  *   )
NB.                                        (  l  l  l  *   )
NB.                                        (  l  l  l  *   )
NB. where
NB.   l              - elements of m×k-matrix L
NB.   vi             - vector v(i)
NB.   τi             - scalar value conj(τ(i))
NB.   (0,...,0,1,vi) - n-vector u(i)
NB.   *              - any value, is not used
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. L * Q = L * (Q)
NB.   ((((mp  unglq)~  trl        ) -: (unmlqrn  trlpick        )) }:"1) LQf
NB.   NB. I = Q * (Q^H)
NB.   (( 0         idmat (<. , ])/)@(0 _1 + $) -: clean@(unmlqrc unglq)) LQf
NB.
NB. Notes:
NB. - monadic case implements LAPACK's DORGLQ, ZUNGLQ
NB. - straightforward O(k*n^3) code:
NB.     Q=.   k  {.   mp/ (idmat n) -"2 |. (+ {:"1 Qf) * (* +)~"0/~"1    }:"1   k  {.   Qf

unglq=: ($:~ 0 _1 <./@:+ $) :(}:"1@(([ (({."1~ -@c) larfbrcfr ])  idmat      @(((] ,  +) #)~ c) appendr  ])&:>/)@([ (   (<;.1~      (0 1:`(GQNB dhs2liso  0 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} #@]))@] , <@( idmat      @((- {.) ,  -~/@]) $))  tru1pick        @((   <./ @, 0 _1 + $) {.   ]))`((0  idmat , ) <:@c)@.(0 e. (, 0 _1 + $)))

NB. ---------------------------------------------------------
NB. ungql
NB.
NB. Description:
NB.   Generate a matrix with orthonormal columns from output
NB.   of geqlf
NB.
NB. Syntax:
NB.   Q=. [k] ungql QfL
NB. where
NB.   QfL - (m+1)×n-matrix, the output of geqlf
NB.   Q   - m×k-matrix with orthonormal columns, which is
NB.         defined as the last k columns of the product of k
NB.         elementary reflectors H(i) of order m:
NB.           Q = Π{H(i),i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   k   ≤ m, optional, default is min(m,n)
NB.
NB. Storage layout:
NB.   QfL for m=3, n=7:               QfL for m=7, n=3:
NB.     (  *  *  *  *  τ0 τ1 τ2  )      (  τ0 τ1 τ2  )
NB.     (  l  l  l  l  l  v1 v2  )      (  v0 v1 v2  )
NB.     (  l  l  l  l  l  l  v2  )      (  v0 v1 v2  )
NB.     (  l  l  l  l  l  l  l   )      (  v0 v1 v2  )
NB.                                     (  v0 v1 v2  )
NB.                                     (  l  v1 v2  )
NB.                                     (  l  l  v2  )
NB.                                     (  l  l  l   )
NB. where
NB.   l              - elements of k×n-matrix L
NB.   vi             - vector v(i)
NB.   τi             - scalar value τ(i)
NB.   (vi,1,0,...,0) - m-vector u(i)
NB.   *              - any value, is not used
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. Q * L = (Q) * L
NB.   ((((mp~ ungql)~ (trl~ -~/@$)) -: (unmqlln (trlpick~ -~/@$))) }.  ) QfL
NB.   NB. I = (Q^H) * Q
NB.   (((0 <. -~/) idmat ([ , <.)/)@(_1 0 + $) -: clean@(unmqllc ungql)) QfL
NB.
NB. Notes:
NB. - monadic case implements LAPACK's DORGQL, ZUNGQL
NB. - straightforward O(k*m^3) code:
NB.     Q=. (-k) {."1 mp/ (idmat m) -"2 |. (  {.   Qf) * (* +) "0/~"1 |: }.   (-k) {."1 Qf

ungql=: ($:~ _1 0 <./@:+ $) :(}.  @(([ (({.  ~   #) larfblnbc ]) (idmat~ -~/)@(((] ,~ +) c)~ #) stitcht~ ])&:>/)@([ (|.@(<;.2~ '' ; (0 1:`(GQNB dhs2liso _1 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} c@]))@] , <@((idmat~ -~/)@((- {:) ,~ - /@]) $)) (tru1pick~ -~/@$)@((-@(<./)@, _1 0 + $) {."1 ]))`((-  idmat ,~) <:@#)@.(0 e. (, _1 0 + $)))

NB. ---------------------------------------------------------
NB. ungqr
NB.
NB. Description:
NB.   Generate a matrix with orthonormal columns from output
NB.   of geqrf
NB.
NB. Syntax:
NB.   Q=. [k] ungqr QfR
NB. where
NB.   QfR - (m+1)×n-matrix, the output of geqrf
NB.   Q   - m×k-matrix with orthonormal columns, which is
NB.         defined as the first k columns of the product of
NB.         k elementary reflectors H(i) of order m:
NB.           Q = Π{H(i),i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   k   ≤ m, optional, default is min(m,n)
NB.
NB. Storage layout:
NB.   QfR for m=3, n=7:               QfR for m=7, n=3:
NB.     (  r  r  r  r  r  r  r  )       (  r  r  r   )
NB.     (  v0 r  r  r  r  r  r  )       (  v0 r  r   )
NB.     (  v0 v1 r  r  r  r  r  )       (  v0 v1 r   )
NB.     (  τ0 τ1 τ2 *  *  *  *  )       (  v0 v1 v2  )
NB.                                     (  v0 v1 v2  )
NB.                                     (  v0 v1 v2  )
NB.                                     (  v0 v1 v2  )
NB.                                     (  τ0 τ1 τ2  )
NB. where
NB.   r              - elements of k×n-matrix R
NB.   vi             - vector v(i)
NB.   τi             - scalar value τ(i)
NB.   (0,...,0,1,vi) - m-vector u(i)
NB.   *              - any value, is not used
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. Q * R = (Q) * R
NB.   ((((mp~ ungqr)~  tru        ) -: (unmqrln  trupick        )) }:  ) QfR
NB.   NB. I = (Q^H) * Q
NB.   (( 0         idmat ([ , <.)/)@(_1 0 + $) -: clean@(unmqrlc ungqr)) QfR
NB.
NB. Notes:
NB. - monadic case implements LAPACK's DORGQR, ZUNGQR
NB. - straightforward O(k*m^3) code:
NB.     Q=.   k  {."1 mp/ (idmat m) -"2    (  {:   Qf) * (* +) "0/~"1 |: }:     k  {."1 Qf

ungqr=: ($:~ _1 0 <./@:+ $) :(}:  @(([ (({.  ~ -@#) larfblnfc ])  idmat      @(((] ,~ +) c)~ #) stitchb  ])&:>/)@([ (   (<;.1~ '' ; (0 1:`(GQNB dhs2liso  0 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} c@]))@] , <@( idmat      @((- {:) ,~ - /@]) $))  trl1pick        @((   <./ @, _1 0 + $) {."1 ]))`((0  idmat ,~) <:@#)@.(0 e. (, _1 0 + $)))

NB. ---------------------------------------------------------
NB. ungrq
NB.
NB. Description:
NB.   Generate a matrix with orthonormal rows from output of
NB.   gerqf
NB.
NB. Syntax:
NB.   Q=. [k] ungrq RQf
NB. where
NB.   RQf - m×(n+1)-matrix, the output of gerqf
NB.   Q   - k×n-matrix with orthonormal rows, which is
NB.         defined as the last k rows of the product of k
NB.         elementary reflectors H(i) of order n:
NB.           Q = Π{H(i)',i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   k   ≤ n, optional, default is min(m,n)
NB.
NB. Storage layout:
NB.   RQf for m=3, n=7:                 RQf for m=7, n=3:
NB.     (  τ0 v0 v0 v0 v0 r  r  r  )      (  *  r  r  r  )
NB.     (  τ1 v1 v1 v1 v1 v1 r  r  )      (  *  r  r  r  )
NB.     (  τ2 v2 v2 v2 v2 v2 v2 r  )      (  *  r  r  r  )
NB.                                       (  *  r  r  r  )
NB.                                       (  τ0 r  r  r  )
NB.                                       (  τ1 v1 r  r  )
NB.                                       (  τ2 v2 v2 r  )
NB. where
NB.   r              - elements of m×k-matrix R
NB.   vi             - vector v(i)
NB.   τi             - scalar value τ(i)
NB.   (vi,1,0,...,0) - n-vector u(i)
NB.   *              - any value, is not used
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. R * Q = R * (Q)
NB.   ((((mp  ungrq)~ (tru~ -~/@$)) -: (unmrqrn (trupick~ -~/@$))) }."1) RQf
NB.   NB. I = Q * (Q^H)
NB.   (((0 >. -~/) idmat (<. , ])/)@(0 _1 + $) -: clean@(unmrqrc ungrq)) RQf
NB.
NB. Notes:
NB. - monadic case implements LAPACK's DORGRQ, ZUNGRQ
NB. - straightforward O(k*n^3) code:
NB.     Q=. (-k) {.   mp/ (idmat n) -"2    (+ {."1 Qf) * (* +)~"0/~"1    }."1 (-k) {.   Qf

ungrq=: ($:~ 0 _1 <./@:+ $) :(}."1@(([ (({."1~   c) larfbrcbr ]) (idmat~ -~/)@(((] ,  +) #)~ c) appendl~ ])&:>/)@([ (|.@(<;.2~      (0 1:`(GQNB dhs2liso _1 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} #@]))@] , <@((idmat~ -~/)@((- {.) ,  -~/@]) $)) (trl1pick~ -~/@$)@((-@(<./)@, 0 _1 + $) {.   ]))`((-~ idmat , ) <:@c)@.(0 e. (, 0 _1 + $)))

NB. ---------------------------------------------------------
NB. unglz
NB.
NB. Description:
NB.   Generate a matrix with orthonormal rows from output of
NB.   tzlzf
NB.
NB. Syntax:
NB.   Z=. [k] unglz LZf
NB. where
NB.   LZf - m×(n+1)-matrix, the output of tzlzf
NB.   Z   - k×n-matrix with orthonormal rows, which is
NB.         defined as the last k rows of the product of k
NB.         elementary reflectors H(i) of order n:
NB.           Z = Π{H(i)',i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   k   ≤ n, optional, default is m
NB.   m   ≤ n
NB.
NB. Storage layout:
NB.   LZf for m=3, n=7:
NB.     (  τ0 v0 v0 v0 v0 l  0  0   )
NB.     (  τ1 v1 v1 v1 v1 l  l  0   )
NB.     (  τ2 v2 v2 v2 v2 l  l  l   )
NB. where
NB.   l                      - elements of m×n-matrix L
NB.   vi                     - (n-m)-vector v(i)
NB.   τi                     - scalar value conj(τ(i))
NB.   (vi,0,...,0,1,0,...,0) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. L * Z = L * (Z)
NB.   (((mp  unglz)~ -: [ unmlzrn ({."1~ (1 -  c))~) ({."1~ -@#)) LZf
NB.   NB. I = Z * Z^H
NB.   (idmat@# -: clean@(mp  ct))@unglz LZf
NB.
NB. Notes:
NB. - straightforward O(k*n^3) code:
NB.     Z=. (-k) {.   mp/ (idmat n) -"2 |. (+ {."1 Zf) * (* +)~"0/~"1    }."1 (-k) {.   Zf

unglz=: ($:~ #) :(}."1@(larzbrcfr&:>/)@([ (   (<;.1~      0 1:`( GQNB dhs2liso  0 , >:@(>. GQNB >.@%~ -&GQNX)                        )`($~)} #@])@] , <@(idmat~ -~/)@(,  c)) ((-@(<. #) ,  -~/@$@]) {. ]) ,.  (((0 >. -~) idmat <. ,  ]) #))`((-~ idmat , ) <:@c)@.(0 e. (, 0 _1 + $)))

NB. ---------------------------------------------------------
NB. ungzl
NB.
NB. Description:
NB.   Generate a matrix with orthonormal columns from output
NB.   of tzzlf
NB.
NB. Syntax:
NB.   Z=. [k] ungzl ZfL
NB. where
NB.   ZfL - (m+1)×n-matrix, the output of tzzlf
NB.   Z   - m×k-matrix with orthonormal columns, which is
NB.         defined as the first k columns of the product of
NB.         k elementary reflectors H(i) of order m:
NB.           Z = Π{H(i),i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   k   ≤ m, optional, default is n
NB.   n   ≤ m
NB.
NB. Storage layout:
NB.   ZfL for m=7, n=3:
NB.   (  l  0  0   )
NB.   (  l  l  0   )
NB.   (  l  l  l   )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  τ0 τ1 τ2  )
NB. where
NB.   l                      - elements of m×n-matrix L
NB.   vi                     - (n-m)-vector v(i)
NB.   τi                     - scalar value τ(i)
NB.   (0,...,0,1,0,...,0,vi) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. Z * L = (Z) * L
NB.   (((mp~ ungzl)~ -: [ unmzlln ({.  ~ (1 -~ #))~) ({.  ~   c)) ZfL
NB.   NB. I = Z^H * Z
NB.   (idmat@c -: clean@(mp~ ct))@ungzl ZfL
NB.
NB. Notes:
NB. - straightforward O(k*m^3) code:
NB.     Z=.   k  {."1 mp/ (idmat m) -"2 |. (  {:   Zf) * (* +) "0/~"1 |: }:     k  {."1 Zf

ungzl=: ($:~ c) :(}:  @(larzblnbc&:>/)@([ (|.@(<;.2~ '' ; 0 1:`((GQNB dhs2liso _1 , >:@(>. GQNB >.@%~ -&GQNX))`((1 0 $ 0)"_)@.(0 = ]))`($~)} c@])@] , <@ idmat      @(,~ #)) ((  (<. c) ,~ -~/@$@]) {. ]) , ~ (( 0        idmat <. ,~ ]) c))`((0  idmat ,~) <:@#)@.(0 e. (, _1 0 + $)))

NB. ---------------------------------------------------------
NB. ungzr
NB.
NB. Description:
NB.   Generate a matrix with orthonormal columns from output
NB.   of tzzrf
NB.
NB. Syntax:
NB.   Z=. [k] ungzr ZfR
NB. where
NB.   ZfR - (m+1)×n-matrix, the output of tzzrf
NB.   Z   - m×k-matrix with orthonormal columns, which is
NB.         defined as the last k columns of the product of k
NB.         elementary reflectors H(i) of order m:
NB.           Z = Π{H(i),i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   k   ≤ m, optional, default is n
NB.   n   ≤ m
NB.
NB. Storage layout:
NB.   ZfR for m=7, n=3:
NB.   (  τ0 τ1 τ2  )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  r  r  r   )
NB.   (  0  r  r   )
NB.   (  0  0  r   )
NB. where
NB.   r                      - elements of m×n-matrix R
NB.   vi                     - (n-m)-vector v(i)
NB.   τi                     - scalar value τ(i)
NB.   (vi,0,...,0,1,0,...,0) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. Z * R = (Z) * R
NB.   (((mp~ ungzr)~ -: [ unmzrln ({.  ~ (1 -  #))~) ({.  ~ -@c)) ZfR
NB.   NB. I = Z^H * Z
NB.   (idmat@c -: clean@(mp~ ct))@ungzr ZfR
NB.
NB. Notes:
NB. - straightforward O(k*m^3) code:
NB.     Z=. (-k) {."1 mp/ (idmat m) -"2    (  {.   Zf) * (* +) "0/~"1 |: }.   (-k) {."1 Zf

ungzr=: ($:~ c) :(}.  @(larzblnfc&:>/)@([ (   (<;.1~ '' ; 0 1:`((GQNB dhs2liso  0 , >:@(>. GQNB >.@%~ -&GQNX))`((1 0 $ 0)"_)@.(0 = ]))`($~)} c@])@] , <@(idmat~ -~/)@(,~ #)) ((-@(<. c) ,~ - /@$@]) {. ]) ,   (((0 <. - ) idmat <. ,~ ]) c))`((-  idmat ,~) <:@#)@.(0 e. (, _1 0 + $)))

NB. ---------------------------------------------------------
NB. ungrz
NB.
NB. Description:
NB.   Generate a matrix with orthonormal rows from output of
NB.   tzrzf
NB.
NB. Syntax:
NB.   Z=. [k] ungrz RZf
NB. where
NB.   RZf - m×(n+1)-matrix, the output of tzrzf
NB.   Z   - k×n-matrix with orthonormal rows, which is
NB.         defined as the first k rows of the product of k
NB.         elementary reflectors of order n:
NB.           Z = Π{H(i)',i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   k   ≤ n, optional, default is m
NB.   m   ≤ n
NB.
NB. Storage layout:
NB.   RZf for m=3, n=7:
NB.     (  r  r  r  v0 v0 v0 v0 τ0  )
NB.     (  0  r  r  v1 v1 v1 v1 τ1  )
NB.     (  0  0  r  v2 v2 v2 v2 τ2  )
NB. where
NB.   r                      - elements of m×n-matrix R
NB.   vi                     - (n-m)-vector v(i)
NB.   τi                     - scalar value conj(τ(i))
NB.   (0,...,0,1,0,...,0,vi) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. R * Z = R * (Z)
NB.   (((mp  ungrz)~ -: [ unmrzrn ({."1~ (1 -~ c))~) ({."1~   #)) RZf
NB.   NB. I = Z * Z^H
NB.   (idmat@# -: clean@(mp  ct))@ungrz RZf
NB.
NB. Notes:
NB. - straightforward O(k*n^3) code:
NB.     Z=.   k  {.   mp/ (idmat n) -"2    (+ {:"1 Zf) * (* +)~"0/~"1    }:"1   k  {.   Zf

ungrz=: ($:~ #) :(}:"1@(larzbrcbr&:>/)@([ (|.@(<;.2~      0 1:`((GQNB dhs2liso _1 , >:@(>. GQNB >.@%~ -&GQNX))`((0 1 $ 0)"_)@.(0 = ]))`($~)} #@])@] , <@ idmat      @(,  c)) ((  (<. #) ,  - /@$@]) {. ]) ,.~ (( 0        idmat <. ,  ]) #))`((0  idmat , ) <:@c)@.(0 e. (, 0 _1 + $)))

NB. ---------------------------------------------------------
NB. unghrl
NB.
NB. Description:
NB.   Generate an unitary (orthogonal) matrix from output of
NB.   gehrdl
NB.
NB. Syntax:
NB.   Q=. unghrl HQf
NB. where
NB.   HQf - n×(n+1)-matrix with packed H and Qf (see gehrdl)
NB.   Q   - n×n-matrix, the unitary (orthogonal), which is
NB.         defined as the product of (s-1) elementary
NB.         reflectors H(i) of order n:
NB.           Q = Π{H(i)',i=h+s-2:h}
NB.           H(i) = I - v[i]' * τ[i] * v[i]
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. I = Q * Q^H
NB.   (idmat@# -: clean@(mp  ct))@unghrl HQf
NB.
NB. Notes:
NB. - instead of using h and s parameters, the following
NB.   product is really calculating:
NB.     Q = Π{H(i)',i=n-1:0} ,
NB.   where
NB.     H(0:h-1) = H(h+s-1:n-1) = I .
NB.   This approach delivers excessive calculations in rare
NB.   case ((h>0) OR (h+s<n)).

unghrl=: unglq@(|.!.0)

NB. ---------------------------------------------------------
NB. unghru
NB.
NB. Description:
NB.   Generate an unitary (orthogonal) matrix from output of
NB.   gehrdu
NB.
NB. Syntax:
NB.   Q=. unghru HQf
NB. where
NB.   HQf - (n+1)×n-matrix with packed H and Qf (see gehrdu)
NB.   Q   - n×n-matrix, the unitary (orthogonal), which is
NB.         defined as the product of (s-1) elementary
NB.         reflectors H(i) of order n:
NB.           Q = Π{H(i),i=h:h+s-2}
NB.           H(i) = I - v[i] * τ[i] * v[i]'
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. I = Q^H * Q
NB.   (idmat@# -: clean@(mp~ ct))@unghru HQf
NB.
NB. Notes:
NB. - models LAPACK's DORGHR, ZUNGHR
NB. - instead of using h and s parameters, the following
NB.   product is really calculating:
NB.     Q = Π{H(i),i=0:n-1} ,
NB.   where
NB.     H(0:h-1) = H(h+s-1:n-1) = I .
NB.   This approach delivers excessive calculations in rare
NB.   case ((h>0) OR (h+s<n)).

unghru=: ungqr@(0 _1&(|.!.0))

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testungq
NB.
NB. Description:
NB.   Test ungxx by general matrix
NB.
NB. Syntax:
NB.   testungq A
NB. where
NB.   A - m×n-matrix

testungq=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/dorglq'
  load_mttmp_ :: ] 'math/mt/test/lapack2/dorgql'
  load_mttmp_ :: ] 'math/mt/test/lapack2/dorgqr'
  load_mttmp_ :: ] 'math/mt/test/lapack2/dorgrq'
  load_mttmp_ :: ] 'math/mt/test/lapack2/zunglq'
  load_mttmp_ :: ] 'math/mt/test/lapack2/zungql'
  load_mttmp_ :: ] 'math/mt/test/lapack2/zungqr'
  load_mttmp_ :: ] 'math/mt/test/lapack2/zungrq'

  rcond=. (_."_)`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  ks=. ~. 0 1 , (,~ <.@-:) <./ 'm n'=. $ y

  normw=. norm1 Awide=. |:^:(>/@$) y
  normt=. norm1 Atall=. |:^:(</@$) y

  LQf=. gelqf Awide
  QfL=. geqlf Atall
  QfR=. geqrf Atall
  RQf=. gerqf Awide

  NB. LAPACK, real datatype

  ik=. 0
  while. ik < # ks do.
    ('dorglq_mttmp_' tmonad ((   3&{::  (}:"1@] ; ({. {:"1)) 2&{::)`]`(rcond"_)`(_."_)`lqt02)) Awide ; normw ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    ('dorgql_mttmp_' tmonad ((-@(3&{::) (}.  @] ; ({. {.  )) 2&{::)`]`(rcond"_)`(_."_)`qlt02)) Atall ; normt ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    ('dorgqr_mttmp_' tmonad ((   3&{::  (}:  @] ; ({. {:  )) 2&{::)`]`(rcond"_)`(_."_)`qrt02)) Atall ; normt ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    ('dorgrq_mttmp_' tmonad ((-@(3&{::) (}."1@] ; ({. {."1)) 2&{::)`]`(rcond"_)`(_."_)`rqt02)) Awide ; normw ; RQf ; ik { ks
    ik=. >: ik
  end.

  NB. LAPACK, complex datatype

  ik=. 0
  while. ik < # ks do.
    ('zunglq_mttmp_' tmonad ((   3&{::  (}:"1@] ; ({. {:"1)) 2&{::)`]`(rcond"_)`(_."_)`lqt02)) Awide ; normw ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    ('zungql_mttmp_' tmonad ((-@(3&{::) (}.  @] ; ({. {.  )) 2&{::)`]`(rcond"_)`(_."_)`qlt02)) Atall ; normt ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    ('zungqr_mttmp_' tmonad ((   3&{::  (}:  @] ; ({. {:  )) 2&{::)`]`(rcond"_)`(_."_)`qrt02)) Atall ; normt ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    ('zungrq_mttmp_' tmonad ((-@(3&{::) (}."1@] ; ({. {."1)) 2&{::)`]`(rcond"_)`(_."_)`rqt02)) Awide ; normw ; RQf ; ik { ks
    ik=. >: ik
  end.

  NB. mt, any datatype

  ik=. 0
  while. ik < # ks do.
    ('unglq' tdyad ((#@(0&{::))`(   3&{::  {.    2&{:: )`]`(rcond"_)`(_."_)`lqt02)) Awide ; normw ; LQf ; ik { ks
    ik=. >: ik
  end.
  ('unglq' tmonad ((2&{::)`]`(rcond"_)`(_."_)`lqt02)) Awide ; normw ; LQf ; m

  ik=. 0
  while. ik < # ks do.
    ('ungql' tdyad ((c@(0&{::))`(-@(3&{::) {."1 (2&{::))`]`(rcond"_)`(_."_)`qlt02)) Atall ; normt ; QfL ; ik { ks
    ik=. >: ik
  end.
  ('ungql' tmonad ((2&{::)`]`(rcond"_)`(_."_)`qlt02)) Atall ; normt ; QfL ; n

  ik=. 0
  while. ik < # ks do.
    ('ungqr' tdyad ((c@(0&{::))`(   3&{::  {."1 (2&{::))`]`(rcond"_)`(_."_)`qrt02)) Atall ; normt ; QfR ; ik { ks
    ik=. >: ik
  end.
  ('ungqr' tmonad ((2&{::)`]`(rcond"_)`(_."_)`qrt02)) Atall ; normt ; QfR ; n

  ik=. 0
  while. ik < # ks do.
    ('ungrq' tdyad ((#@(0&{::))`(-@(3&{::) {.    2&{:: )`]`(rcond"_)`(_."_)`rqt02)) Awide ; normw ; RQf ; ik { ks
    ik=. >: ik
  end.
  ('ungrq' tmonad ((2&{::)`]`(rcond"_)`(_."_)`rqt02)) Awide ; normw ; RQf ; m

  coerase < 'mttmp'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testungz
NB.
NB. Description:
NB.   Test ungxx by trapezoidal matrix
NB.
NB. Syntax:
NB.   testungz A
NB. where
NB.   A - m×n-matrix

testungz=: 3 : 0
  rcond=. (_."_)`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  ks=. ~. 0 1 , (,~ <.@-:) <./ 'm n'=. $ y

  normw=. norm1 Awide=. |:^:(>/@$) y
  normt=. norm1 Atall=. |:^:(</@$) y

  LZf=. tzlzf Awide
  ZfL=. tzzlf Atall
  ZfR=. tzzrf Atall
  RZf=. tzrzf Awide

  ik=. 0
  while. ik < # ks do.
    ('unglz' tdyad ((3&{::)`(2&{:: )`]`(rcond"_)`(_."_)`lzt02)) Awide ; normw ; LZf ; ik { ks
    ik=. >: ik
  end.
  ('unglz' tmonad ((2&{::)`]`(rcond"_)`(_."_)`lzt02)) Awide ; normw ; LZf ; m

  ik=. 0
  while. ik < # ks do.
    ('ungzl' tdyad ((3&{::)`(2&{::)`]`(rcond"_)`(_."_)`zlt02)) Atall ; normt ; ZfL ; ik { ks
    ik=. >: ik
  end.
  ('ungzl' tmonad ((2&{::)`]`(rcond"_)`(_."_)`zlt02)) Atall ; normt ; ZfL ; n

  ik=. 0
  while. ik < # ks do.
    ('ungzr' tdyad ((3&{::)`(2&{::)`]`(rcond"_)`(_."_)`zrt02)) Atall ; normt ; ZfR ; ik { ks
    ik=. >: ik
  end.
  ('ungzr' tmonad ((2&{::)`]`(rcond"_)`(_."_)`zrt02)) Atall ; normt ; ZfR ; n

  ik=. 0
  while. ik < # ks do.
    ('ungrz' tdyad ((3&{::)`(2&{:: )`]`(rcond"_)`(_."_)`rzt02)) Awide ; normw ; RZf ; ik { ks
    ik=. >: ik
  end.
  ('ungrz' tmonad ((2&{::)`]`(rcond"_)`(_."_)`rzt02)) Awide ; normw ; RZf ; m

  EMPTY
)

NB. ---------------------------------------------------------
NB. testunghr
NB.
NB. Description:
NB.   Test unghrx by square matrix
NB.
NB. Syntax:
NB.   testunghr A
NB. where
NB.   A - n×n-matrix

testunghr=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/dorghr'
  load_mttmp_ :: ] 'math/mt/test/lapack2/zunghr'

  'rcondl rcondu'=. (geconi , gecon1) y

  'norml normu'=. (normi , norm1) y

  HlQf=. (gehrdl~ 0 , #) y
  HuQf=. (gehrdu~ 0 , #) y

  hst01l=: }:@[ (normi hst01) ((ct@[ mp  mp~)  1 trlpick }:"1@(2&{::))~
  hst01u=: }:@[ (norm1 hst01) ((ct@[ mp~ mp ) _1 trupick }:  @(2&{::))~

  unt01l=: (normi unt01 (mp  ct))@]
  unt01u=: (norm1 unt01 (mp~ ct))@]

  ('dorghr_mttmp_' tmonad (((1 ; c ; }: ; }:@{:)@(2&{::))`]`(rcondu"_)`(_."_)`(hst01u >. unt01u))) y ; normu ; HuQf
  ('zunghr_mttmp_' tmonad (((1 ; c ; }: ; }:@{:)@(2&{::))`]`(rcondu"_)`(_."_)`(hst01u >. unt01u))) y ; normu ; HuQf

  ('unghrl'        tmonad ((                      2&{:: )`]`(rcondl"_)`(_."_)`(hst01l >. unt01l))) y ; norml ; HlQf
  ('unghru'        tmonad ((                      2&{:: )`]`(rcondu"_)`(_."_)`(hst01u >. unt01u))) y ; normu ; HuQf

  coerase < 'mttmp'
  erase 'hst01l hst01u unt01l unt01u'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testgq
NB.
NB. Description:
NB.   Adv. to make verb to test ungxxx by matrix of generator
NB.   and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testgq
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testgq_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testgq_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testgq_mt_ 150 200

testgq=: 1 : 'EMPTY [ (testunghr_mt_^:(=/@$) [ testungz_mt_ [ testungq_mt_)@u'
