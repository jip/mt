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
NB. Version: 0.9.0 2012-12-29
NB.
NB. Copyright 2010-2012 Igor Zhuravlov
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

NB. ---------------------------------------------------------
NB. ungl2
NB. ung2l
NB. ung2r
NB. ungr2
NB.
NB. Description:
NB.   Reconstruct matrix Q, augmented by trash vector, from
NB.   its elementary reflectors. Non-blocked version of
NB.   algorithms
NB.
NB. Syntax:
NB.   eQ=. ungxx Qf
NB. where
NB.   Qf - k×(n+1)-matrix for ungx2, or (m+1)×k-matrix for
NB.        ung2x, unit triangular, the Q's factored form
NB.   eQ - Q augmented by trash vector
NB.   Q  - k×n-matrix for ungx2, or m×k-matrix for ung2x,
NB.        with orthonormal rows/columns which is defined as
NB.        the first/last rows/columns of a product of k
NB.        elementary reflectors, see corresp.
NB.        ung{lq,ql,qr,rq}
NB.   k  ≤ n for ungx2, or ≤ m for ung2x
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. by implicit matrix product unmxxxx
NB.   (}:"1@ungl2 -: (unmlqrn  idmat      @(0 _1 + $))) Qf  NB. LQ: Q = I * (Q)
NB.   (}.  @ung2l -: (unmqlln (idmat~ -~/)@(_1 0 + $))) Qf  NB. QL: Q = (Q) * I
NB.   (}:  @ung2r -: (unmqrln  idmat      @(_1 0 + $))) Qf  NB. QR: Q = (Q) * I
NB.   (}."1@ungr2 -: (unmrqrn (idmat~ -~/)@(0 _1 + $))) Qf  NB. RQ: Q = I * (Q)
NB.   NB. orthonormality
NB.   (idmat@# (-: clean) (mp  ct))@:(}:"1)@ungl2 Qf  NB. LQ: I = Q * Q^H
NB.   (idmat@c (-: clean) (mp~ ct))@  }.   @ung2l Qf  NB. QL: I = Q^H * Q
NB.   (idmat@c (-: clean) (mp~ ct))@  }:   @ung2r Qf  NB. QR: I = Q^H * Q
NB.   (idmat@# (-: clean) (mp  ct))@:(}."1)@ungr2 Qf  NB. RQ: I = Q * Q^H
NB.
NB. Application:
NB. - simulate LAPACK's xUNGL2(M,N,K) for M≠K:
NB.     NB. eQ=. k ungl2k Qf
NB.     ungl2k=: ungl2@((1 ; 0     ) setdiag  {.       )
NB. - simulate LAPACK's xUNG2L(M,N,K) for M≠K:
NB.     NB. eQ=. k ung2lk Qf
NB.     ung2lk=: ung2l@((1 ; (-  #)) setdiag ({."1~ -)~)
NB. - simulate LAPACK's xUNG2R(M,N,K) for M≠K:
NB.     NB. eQ=. k ung2rk Qf
NB.     ung2rk=: ung2r@((1 ; 0     ) setdiag  {."1     )
NB. - simulate LAPACK's xUNGR2(M,N,K) for M≠K:
NB.     NB. eQ=. k ungr2k Qf
NB.     ungr2k=: ungr2@((1 ; (-~ c)) setdiag ({.  ~ -)~)
NB.
NB. Notes:
NB. - ungl2 implements LAPACK's DORGL2, ZUNGL2 for M=K
NB. - ung2l implements LAPACK's DORG2L, ZUNG2L for N=K
NB. - ung2r implements LAPACK's DORG2R, ZUNG2R for N=K
NB. - ungr2 implements LAPACK's DORGR2, ZUNGR2 for M=K

ungl2=: (((1 _ ,:~ <:@- &#) ,;.0 [) ([ larfrcfr (, ~ 1 {.~   #)~) 0 ,.  ])^:(#@[) (0 $~ 0 ,  -~/@$)
ung2l=: (((_ 1 ,:~    -~&c) ,;.0 [) ([ larflnbc (,.  1 {.~ -@#)~) 0 , ~ ])^:(c@[) (0 $~ 0 ,~ - /@$)
ung2r=: (((_ 1 ,:~ <:@- &c) ,;.0 [) ([ larflnfc (,.~ 1 {.~   #)~) 0 ,   ])^:(c@[) (0 $~ 0 ,~ - /@$)
ungr2=: (((1 _ ,:~    -~&#) ,;.0 [) ([ larfrcbr (,   1 {.~ -@#)~)       ])^:(#@[) (0 $~ 0 ,  -~/@$)

NB. ---------------------------------------------------------
NB. ungl3
NB. ung3l
NB. ung3r
NB. ungr3
NB.
NB. Description:
NB.   Reconstruct matrix Z, augmented by trash vector, from
NB.   its elementary reflectors. Non-blocked version of
NB.   algorithms
NB.
NB. Syntax:
NB.   eZ=. ungxx Zf
NB. where
NB.   Zf - k×(n+1)-matrix for ungx3, or (m+1)×k-matrix for
NB.        ung3x, the Z represented in factored form
NB.   eZ - Z augmented by trash vector
NB.   Z  - k×n-matrix for ungx3, or m×k-matrix for ung3x,
NB.        with orthonormal rows/columns which is defined as
NB.        the first/last rows/columns of a product of k
NB.        elementary reflectors, see corresp.
NB.        ung{lz,zl,zr,rz}
NB.   k  ≤ n for ungx3, or ≤ m for ung3x
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. by implicit matrix product unmxxxx
NB.   (}."1@ungl3 -: (unmlzrn (idmat~ -~/)@(0 _1 + $))) Zf  NB. LZ: Q = I * (Q)
NB.   (}:  @ung3l -: (unmzlln  idmat      @(_1 0 + $))) Zf  NB. ZL: Q = (q) * I
NB.   (}.  @ung3r -: (unmzrln (idmat~ -~/)@(_1 0 + $))) Zf  NB. ZR: Q = (q) * I
NB.   (}:"1@ungr3 -: (unmrzrn  idmat      @(0 _1 + $))) Zf  NB. RZ: Q = I * (Q)
NB.   NB. orthonormality
NB.   (idmat@# (-: clean) (mp  ct))@:(}."1)@ungl3 Zf  NB. RZ: I = Z * Z^H
NB.   (idmat@c (-: clean) (mp~ ct))@  }:   @ung3l Zf  NB. ZL: I = Z^H * Z
NB.   (idmat@c (-: clean) (mp~ ct))@  }.   @ung3r Zf  NB. ZR: I = Z^H * Z
NB.   (idmat@# (-: clean) (mp  ct))@:(}:"1)@ungr3 Zf  NB. RZ: I = Z * Z^H
NB.
NB. Application:
NB. - change eZ's size k:
NB.     NB. eZ=. k ungl3k Zf
NB.     ungl3k=: ungl3@((1 ; (-~ c)) setdiag ({.  ~ -)~)
NB.     NB. eZ=. k ung3lk Zf
NB.     ung3lk=: ung3l@((1 ; 0     ) setdiag  {."1     )
NB.     NB. eZ=. k ung3rk Zf
NB.     ung3rk=: ung3r@((1 ; (-  #)) setdiag ({."1~ -)~)
NB.     NB. eZ=. k ungr2k Zf
NB.     ungr2k=: ungr2@((1 ; 0     ) setdiag  {.       )

ungl3=: (({  ~ _1 - #) ([ larzrcfr (, ~ 0&(1:`(=i:0:)`($~ #)}))~) ])^:(#@[) (0 (,  $ [) c)
ung3l=: (({"1~      c) ([ larzlnbc (,.  0&(1:`(=i.0:)`($~ #)}))~) ])^:(c@[) (0 (,~ $ [) #)
ung3r=: (({"1~ _1 - c) ([ larzlnfc (,.~ 0&(1:`(=i:0:)`($~ #)}))~) ])^:(c@[) (0 (,~ $ [) #)
ungr3=: (({  ~      #) ([ larzrcbr (,   0&(1:`(=i.0:)`($~ #)}))~) ])^:(#@[) (0 (,  $ [) c)


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
NB.   Q=. unglq LQf
NB. where
NB.   LQf - m×(n+1)-matrix, the output of gelqf
NB.   Q   - k×n-matrix with orthonormal rows, which is
NB.         defined as the first k rows of a product of m
NB.         elementary reflectors of order n:
NB.           Q = Π{H(i)',i=m-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.           H(m-1:k) ≡ H(u(m-1:k),τ(m-1:k)) = H(0,0) = I
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   LQf for m=4, n=5:           LQf for m=5, n=4:
NB.     (  l  v0 v0 v0 v0 τ0  )     (  l  v0 v0 v0 τ0  )
NB.     (  l  l  v1 v1 v1 τ1  )     (  l  l  v1 v1 τ1  )
NB.     (  l  l  l  v2 v2 τ2  )     (  l  l  l  v2 τ2  )
NB.     (  l  l  l  l  v3 τ3  )     (  l  l  l  l  τ3  )
NB.                                 (  l  l  l  l  *   )
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
NB. Application:
NB. - simulate LAPACK's xUNGLQ(M,N,K) for M≠K:
NB.     NB. Q=. k unglqk LQf
NB.     unglqk=:  unglq@  {.
NB.
NB. Notes:
NB. - implements LAPACK's DORGLQ, ZUNGLQ for M=K
NB. - straightforward O(k*m^3) code:
NB.   Q=. k {. mp/ (idmat n) -"2 |. (+ {:"1 Qf) * (* +)"0/~"1 + }:"1 Qf

unglq=: }:"1@((((((GQNB ,  _) ,:~               - &c) ];.0 [) (ungl2@[ ,   larfbrcfr) ]) ({."1~ (- GQNB) - c))^:(GQNB %~ -&#) ungl2@(}.~ 2 #    GQNB  * 0 >. GQNB >.@%~ GQNX -~ #))@ tru1        @({.  ~  0 _1    <./ @:+ $)

NB. ---------------------------------------------------------
NB. ungql
NB.
NB. Description:
NB.   Generate a matrix with orthonormal columns from output
NB.   of geqlf
NB.
NB. Syntax:
NB.   Q=. ungql QfL
NB. where
NB.   QfL - (m+1)×n-matrix, the output of geqlf
NB.   Q   - m×k-matrix with orthonormal columns, which is
NB.         defined as the last k columns of a product of n
NB.         elementary reflectors of order m:
NB.           Q = Π{H(i),i=n-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.           H(n-k-1:0) ≡ H(u(n-k-1:0),τ(0:n-k-1)) = H(0,0) = I
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   QfL for m=4, n=5:           QfL for m=5, n=4:
NB.     (  *  τ1 τ2 τ3 τ4  )        (  τ0 τ1 τ2 τ3  )
NB.     (  l  l  v2 v3 v4  )        (  v0 v1 v2 v3  )
NB.     (  l  l  l  v3 v4  )        (  l  v1 v2 v3  )
NB.     (  l  l  l  l  v4  )        (  l  l  v2 v3  )
NB.     (  l  l  l  l  l   )        (  l  l  l  v3  )
NB.                                 (  l  l  l  l   )
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
NB. Application:
NB. - simulate LAPACK's xUNGQL(M,N,K) for M≠K:
NB.     NB. Q=. k ungqlk QfL
NB.     ungqlk=: (ungql@:({."1)~ -)~
NB.
NB. Notes:
NB. - implements LAPACK's DORGQL, ZUNGQL for N=K
NB. - straightforward O(k*m^3) code:
NB.   Q=. (-k) {."1 mp/ (idmat m) -"2 |. ({. Qf) * (* +)"0/~"1 |: }. Qf

e0=: ([ -~ (negpos"0 $)) {. ]

ungql  =: }.  @((((((GQNB ,~ _) ,:~ (GQNB - 1) +  -~&c) ];.0 [) (ung2l@[ ,.~ larfblnbc) ]) ({.  ~    GQNB  + #)                          )^:(GQNB %~ -&c) ung2l                     @(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ c))@(tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $)
ungql_1=: }.  @((( ((GQNB ,~ _) ,:~ (GQNB - 1) +  -~&c) ];.0 [)              larfblnbc     (- GQNB) ((1 ; (,~ _1 ,~ -~/@$)) setdiag e0) ])^:(GQNB %~ -&c) ung2l                     @(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ c))@(tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $)
ungql_2=: }.  @((( ((GQNB ,~ _) ,:~ (GQNB - 1) +  -~&c) ];.0 [)              larfblnbc     (- GQNB) ((1 ; (,~ _1 ,~ -~/@$)) setdiag e0) ])^:(GQNB %~ -&c) (larfblnbc (idmat~ -~/)@$)@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ c))@(tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $)
ungql_3=: }.  @(                                                                                                                                           larfblnbc (idmat~ -~/)@$                                                )@(tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $)

NB. ---------------------------------------------------------
NB. ungqr
NB.
NB. Description:
NB.   Generate a matrix with orthonormal columns from output
NB.   of geqrf
NB.
NB. Syntax:
NB.   Q=. ungqr QfR
NB. where
NB.   QfR - (m+1)×n-matrix, the output of geqrf
NB.   Q   - m×k-matrix with orthonormal columns, which is
NB.         defined as the first k columns of a product of n
NB.         elementary reflectors of order m:
NB.           Q = Π{H(i),i=0:n-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.           H(k:n-1) ≡ H(u(k:n-1),τ(k:n-1)) = H(0,0) = I
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   QfR for m=4, n=5:           QfR for m=5, n=4:
NB.     (  r  r  r  r  r  )         (  r  r  r  r   )
NB.     (  v0 r  r  r  r  )         (  v0 r  r  r   )
NB.     (  v0 v1 r  r  r  )         (  v0 v1 r  r   )
NB.     (  v0 v1 v2 r  r  )         (  v0 v1 v2 r   )
NB.     (  τ0 τ1 τ2 τ3 *  )         (  v0 v1 v2 v3  )
NB.                                 (  τ0 τ1 τ2 τ3  )
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
NB. Application:
NB. - simulate LAPACK's xUNGQR(M,N,K) for M≠K:
NB.     NB. Q=. k ungqrk QfR
NB.     ungqrk=:  ungqr@:({."1)
NB.
NB. Notes:
NB. - implements LAPACK's DORGQR, ZUNGQR for N=K
NB. - straightforward O(k*m^3) code:
NB.   Q=. k {."1 mp/ (idmat m) -"2 ({: Qf) * (* +)"0/~"1 |: }: Qf

ungqr=: }:  @((((((GQNB ,~ _) ,:~               - &#) ];.0 [) (ung2r@[ ,.  larfblnfc) ]) ({.  ~ (- GQNB) - #))^:(GQNB %~ -&c) ung2r@(}.~ 2 #    GQNB  * 0 >. GQNB >.@%~ GQNX -~ c))@ trl1        @({."1~ _1  0    <./ @:+ $)

NB. ---------------------------------------------------------
NB. ungrq
NB.
NB. Description:
NB.   Generate a matrix with orthonormal rows from output of
NB.   gerqf
NB.
NB. Syntax:
NB.   Q=. ungrq RQf
NB. where
NB.   RQf - m×(n+1)-matrix, the output of gerqf
NB.   Q   - k×n-matrix with orthonormal rows, which is
NB.         defined as the last k rows of a product of m
NB.         elementary reflectors of order n:
NB.           Q = Π{H(i)',i=0:m-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.           H(0:m-k-1) ≡ H(u(0:m-k-1),τ(0:m-k-1)) = H(0,0) = I
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   RQf for m=4, n=5:           RQf for m=5, n=4:
NB.     (  τ0 v0 r  r  r  r  )      (  *  r  r  r  r  )
NB.     (  τ1 v1 v1 r  r  r  )      (  τ1 r  r  r  r  )
NB.     (  τ2 v2 c2 v2 r  r  )      (  τ2 v2 r  r  r  )
NB.     (  τ3 v3 v3 v3 v3 r  )      (  τ3 v3 v3 r  r  )
NB.                                 (  τ4 v4 v4 v4 r  )
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
NB. Application:
NB. - simulate LAPACK's xUNGRQ(M,N,K) for M≠K:
NB.     NB. Q=. k ungrqk RQf
NB.     ungrqk=: (ungrq@  {.   ~ -)~
NB.
NB. Notes:
NB. - implements LAPACK's DORGRQ, ZUNGRQ for M=K
NB. - straightforward O(k*m^3) code:
NB.   Q=. (-k) {. mp/ (idmat n) -"2 (+ {."1 Qf) * (* +)"0/~"1 + }."1 Qf

ungrq=: }."1@((((((GQNB ,  _) ,:~ (GQNB - 1) +  -~&#) ];.0 [) (ungr2@[ , ~ larfbrcbr) ]) ({."1~    GQNB  + c))^:(GQNB %~ -&#) ungr2@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ #))@(trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $)

NB. ---------------------------------------------------------
NB. unglz
NB.
NB. Description:
NB.   Generate a matrix with orthonormal rows from output of
NB.   tzlzf
NB.
NB. Syntax:
NB.   Z=. unglz LZf
NB. where
NB.   LZf - k×(n+1)-matrix, the output of tzlzf
NB.   Z   - k×n-matrix with orthonormal rows, which is
NB.         defined as the last k rows of a product of m
NB.         elementary reflectors of order n:
NB.           Z = Π{H(i)',i=m-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.           H(m-k-1:0) ≡ H(u(m-k-1:0),τ(m-k-1:0)) = H(0,0) = I
NB.   k   ≤ m
NB.   m   ≤ n
NB.
NB. Storage layout:
NB.   LZf for m=3, n=7, l=5:
NB.     (  τ0 v0 v0 v0 v0 l  0  0   )
NB.     (  τ1 v1 v1 v1 v1 l  l  0   )
NB.     (  τ2 v2 v2 v2 v2 l  l  l   )
NB. where
NB.   l                    - elements of L
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value conj(τ(i))
NB.   (vi,0,...0,1,0,..,0) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. L * (Z) = L * Z
NB.   (((unmlzrn (trupick~ -~/@$)) -: ((mp  unglz)~ (tru~ -~/@$))) }."1) LZf
NB.   NB. I = Z * (Z^H)
NB.   (((0 >. -~/) idmat (<. , ])/)@(0 _1 + $) (-: clean) (unmlzrc unglz)) LZf
NB.
NB. Application:
NB. - change Z's size k:
NB.     NB. Z=. k unglzk LZf
NB.     unglzk=: (unglz@  {.   ~ -)~

unglz=: }."1@((((((GQNB ,  _) ,:~               -~&#) ];.0 [) (ungl3@[ , ~ larzbrcfr) ]) ({."1~    GQNB  + c))^:(GQNB %~ -&#) ungl3@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ #))@(idmat@]`(a: <@; dhs2lios@(_1 , ]))`[} #)

NB. ---------------------------------------------------------
NB. ungzl
NB.
NB. Description:
NB.   Generate a matrix with orthonormal columns from output
NB.   of tzzlf
NB.
NB. Syntax:
NB.   Z=. ungzl ZfL
NB. where
NB.   ZfL - (m+1)×k-matrix, the output of tzzlf
NB.   Z   - m×k-matrix with orthonormal columns, which is
NB.         defined as the first k columns of a product of n
NB.         elementary reflectors of order m:
NB.           Q = Π{H(i),i=n-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.           H(n-1:k) ≡ H(u(n-1:k),τ(n-1:k)) = H(0,0) = I
NB.   k   ≤ n
NB.   n   ≤ m
NB.
NB. Storage layout:
NB.   ZfL for m=7, n=3, l=5:
NB.   (  l  0  0   )
NB.   (  l  l  0   )
NB.   (  l  l  l   )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  τ0 τ1 τ2  )
NB. where
NB.   l                    - elements of L
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value τ(i)
NB.   (0,..,0,1,0,...0,vi) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. (Z) * L = Z * L
NB.   (((unmzlln  trupick        ) -: ((mp~ ungzl)~  tru        )) }:  ) ZfL
NB.   NB. I = (Z^H) * Z
NB.   (( 0         idmat ([ , <.)/)@(_1 0 + $) (-: clean) (unmzllc ungzl)) ZfL
NB.
NB. Application:
NB. - change Z's size k:
NB.     NB. Z=. k ungzlk ZfL
NB.     ungzlk=: ungzl@:({."1)

ungzl=: }:  @((((((GQNB ,~ _) ,:~  GQNB      -~ - &c) ];.0 [) (ung3l@[ ,.  larzblnbc) ]) ({.  ~ (- GQNB) - #))^:(GQNB %~ -&c) ung3l@(}.~ 2 #    GQNB  * 0 >. GQNB >.@%~ GQNX -~ c))@(idmat@]`(       dhs2lios@( 0 , ]))`[} c)

NB. ---------------------------------------------------------
NB. ungzr
NB.
NB. Description:
NB.   Generate a matrix with orthonormal columns from output
NB.   of tzzrf
NB.
NB. Syntax:
NB.   Z=. ungzr ZfR
NB. where
NB.   ZfR - (m+1)×k-matrix, the output of tzzrf
NB.   Z   - m×k-matrix with orthonormal columns, which is
NB.         defined as the last k columns of a product of n
NB.         elementary reflectors of order m:
NB.           Q = Π{H(i),i=0:n-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.           H(0:n-k-1) ≡ H(u(0:n-k-1),τ(0:n-k-1)) = H(0,0) = I
NB.   k   ≤ n
NB.   n   ≤ m
NB.
NB. Storage layout:
NB.   ZfR for m=7, n=3, l=5:
NB.   (  τ0 τ1 τ2  )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  v0 v1 v2  )
NB.   (  r  r  r   )
NB.   (  0  r  r   )
NB.   (  0  0  r   )
NB. where
NB.   r                    - elements of R
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value τ(i)
NB.   (vi,0,..,0,1,0,...0) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. (Z) * R = Z * R
NB.   (((unmzrln (trlpick~ -~/@$)) -: ((mp~ ungzr)~ (trl~ -~/@$))) }.  ) ZfR
NB.   NB. I = (Z^H) * Z
NB.   (((0 <. -~/) idmat ([ , <.)/)@(_1 0 + $) (-: clean) (unmzrlc ungzr)) ZfR
NB.
NB. Application:
NB. - change Z's size k:
NB.     NB. Z=. k ungzrk ZfR
NB.     ungzrk=: (ungzr@:({."1)~ -)~

ungzr=: }.  @((((((GQNB ,~ _) ,:~               -~&c) ];.0 [) (ung3r@[ ,.~ larzblnfc) ]) ({.  ~    GQNB  + #))^:(GQNB %~ -&c) ung3r@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ c))@(idmat@]`(       dhs2lios@(_1 , ]))`[} c)

NB. ---------------------------------------------------------
NB. ungrz
NB.
NB. Description:
NB.   Generate a matrix with orthonormal rows from output of
NB.   tzrzf
NB.
NB. Syntax:
NB.   Z=. ungrz RZf
NB. where
NB.   RZf - k×(n+1)-matrix, the output of tzrzf
NB.   Z   - k×n-matrix with orthonormal rows, which is
NB.         defined as the first k rows of a product of m
NB.         elementary reflectors of order n:
NB.           Q = Π{H(i)',i=0:m-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.           H(k:m-1) ≡ H(u(k:m-1),τ(k:m-1)) = H(0,0) = I
NB.   k   ≤ m
NB.   m   ≤ n
NB.
NB. Storage layout:
NB.   RZf for m=3, n=7, l=5:
NB.     (  r  r  r  v0 v0 v0 v0 τ0  )
NB.     (  0  r  r  v1 v1 v1 v1 τ1  )
NB.     (  0  0  r  v2 v2 v2 v2 τ2  )
NB. where
NB.   r                    - elements of R
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value conj(τ(i))
NB.   (0,..,0,1,0,...0,vi) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. R * (Z) = R * Z
NB.   (((unmrzrn  trlpick        ) -: ((mp  ungrz)~  trl        )) }:"1) RZf
NB.   NB. I = Z * (Z^H)
NB.   (( 0         idmat (<. , ])/)@(0 _1 + $) (-: clean) (unmrzrc ungrz)) RZf
NB.
NB. Application:
NB. - change Z's size k:
NB.     NB. Z=. k ungrzk RZf
NB.     ungrzk=: ungrz@  {.

ungrz=: }:"1@((((((GQNB ,  _) ,:~               - &#) ];.0 [) (ungr3@[ ,   larzbrcbr) ]) ({."1~ (- GQNB) - c))^:(GQNB %~ -&#) ungr3@(}.~ 2 #    GQNB  * 0 >. GQNB >.@%~ GQNX -~ #))@(idmat@]`(a: <@; dhs2lios@( 0 , ]))`[} #)

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
NB.   Q   - n×n-matrix, an unitary (orthogonal)
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
NB.   Q   - n×n-matrix, an unitary (orthogonal)
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
NB. - for unglq, ungrq:
NB.     berr := ||Q * Q^H - I|| / (FP_EPS * n)
NB. - for ungql, ungqr:
NB.     berr := ||Q^H * Q - I|| / (FP_EPS * m)

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
NB. - for unglz, ungrz:
NB.     berr := ||Z * (Z^H) - I|| / (FP_EPS * n)
NB. - for ungql, ungqr:
NB.     berr := ||(Z^H) * Z - I|| / (FP_EPS * m)

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
NB. - for unghrl:
NB.     berr := ||Q * Q^H - I|| / (FP_EPS * n)
NB. - for ungql, ungqr :
NB.     berr := ||Q^H * Q - I|| / (FP_EPS * m)

testunghr=: 3 : 0
  ('unghrl' tmonad ((gehrdl~ 0 , #)`]`(uncon1@])`(_."_)`(mp  gqvberr c))) y
  ('unghru' tmonad ((gehrdu~ 0 , #)`]`(uncon1@])`(_."_)`(mp~ gqvberr #))) y

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

testgq=: 1 : 'EMPTY_mt_ [ (testunghr_mt_^:(=/@$) [ testungz_mt_ [ testungq_mt_)@u'

NB. ungl2=: (((1 _ ,:~ <:@- &#) ,;.0 [) ([ larfrcfr (, ~ 1 {.~   #)~) 0 ,.  ])^:(#@[) (0 $~ 0 ,  -~/@$)
NB. ung2l=: (((_ 1 ,:~    -~&c) ,;.0 [) ([ larflnbc (,.  1 {.~ -@#)~) 0 , ~ ])^:(c@[) (0 $~ 0 ,~ - /@$)
NB. ung2r=: (((_ 1 ,:~ <:@- &c) ,;.0 [) ([ larflnfc (,.~ 1 {.~   #)~) 0 ,   ])^:(c@[) (0 $~ 0 ,~ - /@$)
NB. ungr2=: (((1 _ ,:~    -~&#) ,;.0 [) ([ larfrcbr (,   1 {.~ -@#)~)       ])^:(#@[) (0 $~ 0 ,  -~/@$)
NB. 
NB. ungl3=: (({  ~ _1 - #) ([ larzrcfr (, ~ 0&(1:`(=i:0:)`($~ #)}))~) ])^:(#@[) (0 (,  $ [) c)
NB. ung3l=: (({"1~      c) ([ larzlnbc (,.  0&(1:`(=i.0:)`($~ #)}))~) ])^:(c@[) (0 (,~ $ [) #)
NB. ung3r=: (({"1~ _1 - c) ([ larzlnfc (,.~ 0&(1:`(=i:0:)`($~ #)}))~) ])^:(c@[) (0 (,~ $ [) #)
NB. ungr3=: (({  ~      #) ([ larzrcbr (,   0&(1:`(=i.0:)`($~ #)}))~) ])^:(#@[) (0 (,  $ [) c)
NB. 
NB. unglq=: }:"1@((((((GQNB ,  _) ,:~               - &c) ];.0 [) (ungl2@[ ,   larfbrcfr) ]) ({."1~ (- GQNB) - c))^:(GQNB %~ -&#) ungl2@(}.~ 2 #    GQNB  * 0 >. GQNB >.@%~ GQNX -~ #))@ tru1        @({.  ~  0 _1    <./ @:+ $)
NB. ungql=: }.  @((((((GQNB ,~ _) ,:~ (GQNB - 1) +  -~&c) ];.0 [) (ung2l@[ ,.~ larfblnbc) ]) ({.  ~    GQNB  + #))^:(GQNB %~ -&c) ung2l@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ c))@(tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $)
NB. ungqr=: }:  @((((((GQNB ,~ _) ,:~               - &#) ];.0 [) (ung2r@[ ,.  larfblnfc) ]) ({.  ~ (- GQNB) - #))^:(GQNB %~ -&c) ung2r@(}.~ 2 #    GQNB  * 0 >. GQNB >.@%~ GQNX -~ c))@ trl1        @({."1~ _1  0    <./ @:+ $)
NB. ungrq=: }."1@((((((GQNB ,  _) ,:~ (GQNB - 1) +  -~&#) ];.0 [) (ungr2@[ , ~ larfbrcbr) ]) ({."1~    GQNB  + c))^:(GQNB %~ -&#) ungr2@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ #))@(trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $)
NB. 
NB. unglz=: }."1@((((((GQNB ,  _) ,:~               -~&#) ];.0 [) (ungl3@[ , ~ larzbrcfr) ]) ({."1~    GQNB  + c))^:(GQNB %~ -&#) ungl3@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ #))@(idmat@]`(a: <@; dhs2lios@(_1 , ]))`[} #)
NB. ungzl=: }:  @((((((GQNB ,~ _) ,:~  c@]              ) ];.0 [) (ung3l@[ ,.  larzblnbc) ]) ({.  ~ (- GQNB) - #))^:(GQNB %~ -&c) ung3l@(}."1~ (- GQNB)  * 0 >. GQNB >.@%~ GQNX -~ c))@(idmat@]`(       dhs2lios@( 0 , ]))`[} c)
NB. ungzr=: }.  @((((((GQNB ,~ _) ,:~               -~&c) ];.0 [) (ung3r@[ ,.~ larzblnfc) ]) ({.  ~    GQNB  + #))^:(GQNB %~ -&c) ung3r@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ c))@(idmat@]`(       dhs2lios@(_1 , ]))`[} c)
NB. ungrz=: }:"1@((((((GQNB ,  _) ,:~               - &#) ];.0 [) (ungr3@[ ,   larzbrcbr) ]) ({."1~ (- GQNB) - c))^:(GQNB %~ -&#) ungr3@(}.~ 2 #    GQNB  * 0 >. GQNB >.@%~ GQNX -~ #))@(idmat@]`(a: <@; dhs2lios@( 0 , ]))`[} #)
