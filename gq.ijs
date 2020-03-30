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
NB. testungq   Test ung{lq,ql,qr,rq} by general matrix
NB. testungz   Test ung{lz,zl,zr,rz} by trapezoidal matrix
NB. testunghr  Test unghrx by square matrix
NB. testgq     Adv. to make verb to test ungxxx by matrix of
NB.            generator and shape given
NB.
NB. Version: 0.10.5 2020-03-30
NB.
NB. Copyright 2010-2020 Igor Zhuravlov
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

gqvberr=: 2 : '(norm1_mt_@(<: upddiag_mt_)@(u ct_mt_) % FP_EPS_mt_ * v)@]'  NB. conj. to form verb to calc. berr

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
NB.         defined as the first k rows of a product of k
NB.         elementary reflectors of order n:
NB.           Q = Π{H(i)',i=k-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
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
NB.   NB. L * (Q) = L * Q
NB.   (((unmlqrn  trlpick        ) -: ((mp  unglq)~  trl        )) }:"1) LQf
NB.   NB. I = Q * (Q^H)
NB.   (( 0         idmat (<. , ])/)@(0 _1 + $) (-: clean) (unmlqrc unglq)) LQf
NB.
NB. Notes:
NB. - implements LAPACK's DORGLQ, ZUNGLQ
NB. - straightforward O(k*m^3) code:
NB.   Q=. k {. mp/ (idmat n) -"2 |. (+ {:"1 Qf) * (* +)"0/~"1 + }:"1 Qf

unglq=: ($:~ 0 _1 <./@:+ $) :(}:"1@(([ (({."1~ -@c) larfbrcfr ])  idmat      @(((] ,  +) #)~ c) appendr  ])&:>/)@([ (   (<;.1~      (0 1:`(GQNB dhs2liso  0 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} #@]))@] , <@( idmat      @((- {.) ,  -~/@]) $))  tru1pick        @((   <./ @, 0 _1 + $) {.   ])))

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
NB.         defined as the last k columns of a product of k
NB.         elementary reflectors of order m:
NB.           Q = Π{H(i),i=k-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
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
NB.   NB. (Q) * L = Q * L
NB.   (((unmqlln (trlpick~ -~/@$)) -: ((mp~ ungql)~ (trl~ -~/@$))) }.  ) QfL
NB.   NB. I = (Q^H) * Q
NB.   (((0 <. -~/) idmat ([ , <.)/)@(_1 0 + $) (-: clean) (unmqllc ungql)) QfL
NB.
NB. Notes:
NB. - implements LAPACK's DORGQL, ZUNGQL
NB. - straightforward O(k*m^3) code:
NB.   Q=. (-k) {."1 mp/ (idmat m) -"2 |. ({. Qf) * (* +)"0/~"1 |: }. Qf

ungql=: ($:~ _1 0 <./@:+ $) :(}.  @(([ (({.  ~   #) larfblnbc ]) (idmat~ -~/)@(((] ,~ +) c)~ #) stitcht~ ])&:>/)@([ (|.@(<;.2~ '' ; (0 1:`(GQNB dhs2liso _1 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} c@]))@] , <@((idmat~ -~/)@((- {:) ,~ - /@]) $)) (tru1pick~ -~/@$)@((-@(<./)@, _1 0 + $) {."1 ])))

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
NB.         defined as the first k columns of a product of k
NB.         elementary reflectors of order m:
NB.           Q = Π{H(i),i=0:k-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
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
NB.   NB. (Q) * R = Q * R
NB.   (((unmqrln  trupick        ) -: ((mp~ ungqr)~  tru        )) }:  ) QfR
NB.   NB. I = (Q^H) * Q
NB.   (( 0         idmat ([ , <.)/)@(_1 0 + $) (-: clean) (unmqrlc ungqr)) QfR
NB.
NB. Notes:
NB. - implements LAPACK's DORGQR, ZUNGQR
NB. - straightforward O(k*m^3) code:
NB.   Q=. k {."1 mp/ (idmat m) -"2 ({: Qf) * (* +)"0/~"1 |: }: Qf

ungqr=: ($:~ _1 0 <./@:+ $) :(}:  @(([ (({.  ~ -@#) larfblnfc ])  idmat      @(((] ,~ +) c)~ #) stitchb  ])&:>/)@([ (   (<;.1~ '' ; (0 1:`(GQNB dhs2liso  0 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} c@]))@] , <@( idmat      @((- {:) ,~ - /@]) $))  trl1pick        @((   <./ @, _1 0 + $) {."1 ])))

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
NB.         defined as the last k rows of a product of k
NB.         elementary reflectors of order n:
NB.           Q = Π{H(i)',i=0:k-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
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
NB.   NB. R * (Q) = R * Q
NB.   (((unmrqrn (trupick~ -~/@$)) -: ((mp  ungrq)~ (tru~ -~/@$))) }."1) RQf
NB.   NB. I = Q * (Q^H)
NB.   (((0 >. -~/) idmat (<. , ])/)@(0 _1 + $) (-: clean) (unmrqrc ungrq)) RQf
NB.
NB. Notes:
NB. - implements LAPACK's DORGRQ, ZUNGRQ
NB. - straightforward O(k*m^3) code:
NB.   Q=. (-k) {. mp/ (idmat n) -"2 (+ {."1 Qf) * (* +)"0/~"1 + }."1 Qf

ungrq=: ($:~ 0 _1 <./@:+ $) :(}."1@(([ (({."1~   c) larfbrcbr ]) (idmat~ -~/)@(((] ,  +) #)~ c) appendl~ ])&:>/)@([ (|.@(<;.2~      (0 1:`(GQNB dhs2liso _1 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} #@]))@] , <@((idmat~ -~/)@((- {.) ,  -~/@]) $)) (trl1pick~ -~/@$)@((-@(<./)@, 0 _1 + $) {.   ])))

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
NB.         defined as the last k rows of a product of k
NB.         elementary reflectors of order n:
NB.           Z = Π{H(i)',i=k-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.   k   ≤ n, optional, default is min(m,n)
NB.   m   ≤ n
NB.
NB. Storage layout:
NB.   LZf for m=3, n=7:
NB.     (  τ0 v0 v0 v0 v0 l  0  0   )
NB.     (  τ1 v1 v1 v1 v1 l  l  0   )
NB.     (  τ2 v2 v2 v2 v2 l  l  l   )
NB. where
NB.   l                    - elements of m×n-matrix L
NB.   vi                   - (n-m)-vector v(i)
NB.   τi                   - scalar value conj(τ(i))
NB.   (vi,0,...0,1,0,..,0) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = L * Z
NB.   (-: (trlpick@:({."1~ -@#) mp  unglz)@tzlzf) A
NB.   NB. I = Z * Z^H
NB.   (idmat@# -: (mp  ct))@unglz LZf

unglz=: ($:~ #) :(}."1@(larzbrcfr&:>/)@([ (   (<;.1~      0 1:`(GQNB dhs2liso  0 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} #@])@] , <@(idmat~ -~/)@(,  c)) ((-@(<. #) ,  -~/@$@]) {. ]) ,.  (((0 >. -~) idmat <. ,  ]) #)))

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
NB.         defined as the first k columns of a product of k
NB.         elementary reflectors of order m:
NB.           Z = Π{H(i),i=k-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.   k   ≤ m, optional, default is min(m,n)
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
NB.   l                    - elements of m×n-matrix L
NB.   vi                   - (n-m)-vector v(i)
NB.   τi                   - scalar value τ(i)
NB.   (0,..,0,1,0,...0,vi) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = Z * L
NB.   (-: (trlpick@:({.  ~   c) mp~ ungzl)@tzzlf) A
NB.   NB. I = Z^H * Z
NB.   (idmat@c -: (mp~ ct))@ungzl ZfL

ungzl=: ($:~ c) :(}:  @(larzblnbc&:>/)@([ (|.@(<;.2~ '' ; 0 1:`(GQNB dhs2liso _1 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} c@])@] , <@ idmat      @(,~ #)) ((  (<. c) ,~ -~/@$@]) {. ]) , ~ (( 0        idmat <. ,~ ]) c)))

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
NB.         defined as the last k columns of a product of k
NB.         elementary reflectors of order m:
NB.           Z = Π{H(i),i=0:k-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.   k   ≤ m, optional, default is min(m,n)
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
NB.   r                    - elements of m×n-matrix R
NB.   vi                   - (n-m)-vector v(i)
NB.   τi                   - scalar value τ(i)
NB.   (vi,0,..,0,1,0,...0) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = Z * R
NB.   (-: (trupick@:({.  ~ -@c) mp~ ungzr)@tzzrf) A
NB.   NB. I = Z^H * Z
NB.   (idmat@c -: (mp~ ct))@ungzr ZfR

ungzr=: ($:~ c) :(}.  @(larzblnfc&:>/)@([ (   (<;.1~ '' ; 0 1:`(GQNB dhs2liso  0 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} c@])@] , <@(idmat~ -~/)@(,~ #)) ((-@(<. c) ,~ - /@$@]) {. ]) ,   (((0 <. - ) idmat <. ,~ ]) c)))

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
NB.         defined as the first k rows of a product of k
NB.         elementary reflectors of order n:
NB.           Z = Π{H(i)',i=0:k-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.   k   ≤ n, optional, default is min(m,n)
NB.   m   ≤ n
NB.
NB. Storage layout:
NB.   RZf for m=3, n=7:
NB.     (  r  r  r  v0 v0 v0 v0 τ0  )
NB.     (  0  r  r  v1 v1 v1 v1 τ1  )
NB.     (  0  0  r  v2 v2 v2 v2 τ2  )
NB. where
NB.   r                    - elements of m×n-matrix R
NB.   vi                   - (n-m)-vector v(i)
NB.   τi                   - scalar value conj(τ(i))
NB.   (0,..,0,1,0,...0,vi) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = R * Z
NB.   (-: (trupick@:({."1~   #) mp  ungrz)@tzrzf) A
NB.   NB. I = Z * Z^H
NB.   (idmat@# -: (mp  ct))@ungrz RZf

ungrz=: ($:~ #) :(}:"1@(larzbrcbr&:>/)@([ (|.@(<;.2~      0 1:`(GQNB dhs2liso _1 , >:@(>. GQNB >.@%~ -&GQNX))`($~)} #@])@] , <@ idmat      @(,  c)) ((  (<. #) ,  - /@$@]) {. ]) ,.~ (( 0        idmat <. ,  ]) #)))

NB. ---------------------------------------------------------
NB. unghrl
NB.
NB. Description:
NB.   Generate an unitary (orthogonal) matrix Q which is
NB.   defined as the product of elementary reflectors of
NB.   order n, as returned by gehrdl:
NB.     Q = Π{H(i)',i=f+s-2:f} .
NB.
NB. Syntax:
NB.   Q=. unghrl HQf
NB. where
NB.   HQf - n×(n+1)-matrix with packed H and Qf (see gehrdl)
NB.   Q   - n×n-matrix, the unitary (orthogonal)
NB.
NB. Notes:
NB. - instead of using f and s parameters, the following
NB.   product is really calculating:
NB.     Q = Π{H(i)',i=n-1:0} ,
NB.   where
NB.     H(0:f-1) = H(f+s-1:n-1) = I .
NB.   This approach delivers excessive calculations in rare
NB.   case ((f>0) OR (f+s<n)).

unghrl=: unglq@(|.!.0)

NB. ---------------------------------------------------------
NB. unghru
NB.
NB. Description:
NB.   Generate an unitary (orthogonal) matrix Q which is
NB.   defined as the product of elementary reflectors of
NB.   order n, as returned by gehrdu:
NB.     Q = Π{H(i),i=f:f+s-2} .
NB.
NB. Syntax:
NB.   Q=. unghru HQf
NB. where
NB.   HQf - (n+1)×n-matrix with packed H and Qf (see gehrdu)
NB.   Q   - n×n-matrix, the unitary (orthogonal)
NB.
NB. Notes:
NB. - models LAPACK's DORGHR, ZUNGHR
NB. - instead of using f and s parameters, the following
NB.   product is really calculating:
NB.     Q = Π{H(i),i=0:n-1} ,
NB.   where
NB.     H(0:f-1) = H(f+s-1:n-1) = I .
NB.   This approach delivers excessive calculations in rare
NB.   case ((f>0) OR (f+s<n)).

unghru=: ungqr@(0 _1&(|.!.0))

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testungq
NB.
NB. Description:
NB.   Test Q generation algorithms by general matrix
NB.
NB. Syntax:
NB.   testungq A
NB. where
NB.   A - m×n-matrix
NB.
NB. Formula:
NB. - for unglq, ungrq: berr := ||Q * Q^H - I|| / (FP_EPS * n)
NB. - for ungql, ungqr: berr := ||Q^H * Q - I|| / (FP_EPS * m)

testungq=: 3 : 0
  rcond=. (_."_)`gecon1@.(=/@$) y  NB. meaninigful for square matrices only

  ('unglq' tmonad (gelqf`]`(rcond"_)`(_."_)`(mp  gqvberr c))) y
  ('ungql' tmonad (geqlf`]`(rcond"_)`(_."_)`(mp~ gqvberr #))) y
  ('ungqr' tmonad (geqrf`]`(rcond"_)`(_."_)`(mp~ gqvberr #))) y
  ('ungrq' tmonad (gerqf`]`(rcond"_)`(_."_)`(mp  gqvberr c))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testungz
NB.
NB. Description:
NB.   Test Z generation algorithms by trapezoidal matrix
NB.
NB. Syntax:
NB.   testungq A
NB. where
NB.   A - m×n-matrix
NB.
NB. Formula:
NB. - for unglz, ungrz: berr := ||Z * Z^H - I|| / (FP_EPS * n)
NB. - for ungzl, ungzr: berr := ||Z^H * Z - I|| / (FP_EPS * m)

testungz=: 3 : 0
  rcond=. (_."_)`gecon1@.(=/@$) y  NB. meaninigful for square matrices only
  Awide=. |:^:(>/@$) y
  Atall=. |:^:(</@$) y

  ('unglz' tmonad (tzlzf`]`(rcond"_)`(_."_)`(mp  gqvberr c))) (trl~ -~/@$) Awide
  ('ungzl' tmonad (tzzlf`]`(rcond"_)`(_."_)`(mp~ gqvberr #)))  trl         Atall
  ('ungzr' tmonad (tzzrf`]`(rcond"_)`(_."_)`(mp~ gqvberr #))) (tru~ -~/@$) Atall
  ('ungrz' tmonad (tzrzf`]`(rcond"_)`(_."_)`(mp  gqvberr c)))  tru         Awide

  EMPTY
)

NB. ---------------------------------------------------------
NB. testunghr
NB.
NB. Description:
NB.   Test Q generation algorithms by square matrix
NB.
NB. Syntax:
NB.   testunghr A
NB. where
NB.   A - n×n-matrix
NB.
NB. Formula:
NB. - for unghrl: berr := ||Q * Q^H - I|| / (FP_EPS * n)
NB. - for unghru: berr := ||Q^H * Q - I|| / (FP_EPS * n)

testunghr=: 3 : 0
  rcond=. uncon1 y

  ('unghrl' tmonad ((gehrdl~ 0 , #)`]`(rcond"_)`(_."_)`(mp  gqvberr #))) y
  ('unghru' tmonad ((gehrdu~ 0 , #)`]`(rcond"_)`(_."_)`(mp~ gqvberr #))) y

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
