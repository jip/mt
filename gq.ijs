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

GQNB=: 3 NB. 32   NB. block size limit
GQNX=: 5 NB. 128  NB. crossover point, GQNX ≥ GQNB

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
NB.   (] -: ( trl        @:(}:"1) mp  }:"1@ungl2@ tru1        @({.  ~  0 _1    <./ @:+ $))@gelqf) A
NB.   (] -: ((trl~ -~/@$)@  }.    mp~ }.  @ung2l@(tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))@geqlf) A
NB.   (] -: ( tru        @  }:    mp~ }:  @ung2r@ trl1        @({."1~ _1  0    <./ @:+ $))@geqrf) A
NB.   (] -: ((tru~ -~/@$)@:(}."1) mp  }."1@ungr2@(trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))@gerqf) A
NB.
NB. Application:
NB. - simulate LAPACK's xUNGL2(M,N,K) for M≠K:
NB.     NB. eQf=. k ungl2k Qf
NB.     ungl2k=: ungl2@{.
NB. - simulate LAPACK's xUNG2L(M,N,K) for M≠K:
NB.     NB. eQf=. k ung2lk Qf
NB.     ung2lk=: (ung2l@:({."1)~ -)~
NB. - simulate LAPACK's xUNG2R(M,N,K) for M≠K:
NB.     NB. eQf=. k ung2rk Qf
NB.     ung2rk=: ung2r@:({."1)
NB. - simulate LAPACK's xUNGR2(M,N,K) for M≠K:
NB.     NB. eQf=. k ungr2k Qf
NB.     ungr2k=: (ungr2@{.~ -)~
NB.
NB. Notes:
NB. - ungl2 implements LAPACK's DORGL2, ZUNGL2
NB. - ung2l implements LAPACK's DORG2L, ZUNG2L
NB. - ung2r implements LAPACK's DORG2R, ZUNG2R
NB. - ungr2 implements LAPACK's DORGR2, ZUNGR2

ungl2=: (((1 _ ,:~ <:@- &#) ,;.0 [) (((>:@]  0} *) +@-@{:)@[ ,   larfrcfr) 0 ,.  ])^:(#@[) (0 $~ 0 ,  -~/@$)
ung2l=: (((_ 1 ,:~    -~&c) ,;.0 [) (((>:@] _1} *)   -@{.)@[ ,.~ larflnbc) 0 , ~ ])^:(c@[) (0 $~ 0 ,~ - /@$)
ung2r=: (((_ 1 ,:~ <:@- &c) ,;.0 [) (((>:@]  0} *)   -@{:)@[ ,.  larflnfc) 0 ,   ])^:(c@[) (0 $~ 0 ,~ - /@$)
ungr2=: (((1 _ ,:~    -~&#) ,;.0 [) (((>:@] _1} *) +@-@{.)@[ , ~ larfrcbr) 0 ,.~ ])^:(#@[) (0 $~ 0 ,  -~/@$)

NB. form e(i), then append to eQ(i)
ungr2_2=: ([ larfrcbr ((, (1:`[`(0 $~ #@])}~ 1&(=i:1:)))~ _ _&{. :: (0 $~ 0 , #)))/@(0 , ~ |.)

NB. append row of zeros to eQ(i), then place unit in that last row to turn it into e(i)
ungr2_3=: ([ larfrcbr 1:`(_1 <@; 1 (=i:1:) [)`((, 0:)~  _ _&{. :: (0 $~ 0 , #))})/@(0 , ~ |.)

NB. original ungr2 with idea of ungr2_2 and ungr2_3 of unified vectors processing
ungr2_4=: (((1 _ ,:~    -~&#) ,;.0 [) ([ larfrcbr (, 1 {.~ -@#)~) 0 ,.~ ])^:(#@[) (0 $~ 0 ,  -~/@$)

NB. ungr2_4 with two steps ((eQ(i) ,. 0) , e(i)) merged into one: (eQ(i) , e(i))
ungr2_5=: (((1 _ ,:~    -~&#) ,;.0 [) ([ larfrcbr (, 1 {.~ -@#)~) ])^:(#@[) (0 $~ 0 ,  -~/@$)

NB. ---------------------------------------------------------
NB. ungl3
NB. ung3l
NB. ung3r
NB. ungr3
NB.
NB. Description:
NB.   Reconstruct matrix Q, augmented by trash vector, from
NB.   its elementary reflectors. Non-blocked version of
NB.   algorithms
NB.
NB. Syntax:
NB.   eQ=. ungxx Qf
NB. where
NB.   Qf - k×(n+1)-matrix for ungx3, or (m+1)×k-matrix for
NB.        ung3x, unit triangular, the Q's factored form
NB.   eQ - Q augmented by trash vector
NB.   Q  - k×n-matrix for ungx3, or m×k-matrix for ung3x,
NB.        with orthonormal rows/columns which is defined as
NB.        the first/last rows/columns of a product of k
NB.        elementary reflectors, see corresp.
NB.        ung{lz,zl,zr,rz}
NB.   k  ≤ n for ungx3, or ≤ m for ung3x
NB.
NB. Assertions (with appropriate comparison tolerance):##############
NB.   (] -: ( trl        @:(}:"1) mp  }:"1@ungl2@ tru1        @({.  ~  0 _1    <./ @:+ $))@gelqf) A
NB.   (] -: ((trl~ -~/@$)@  }.    mp~ }.  @ung2l@(tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))@geqlf) A
NB.   (] -: ( tru        @  }:    mp~ }:  @ung2r@ trl1        @({."1~ _1  0    <./ @:+ $))@geqrf) A
NB.   (] -: ((tru~ -~/@$)@:(}."1) mp  }."1@ungr2@(trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))@gerqf) A
NB.
NB. Application:#################
NB. - change eQ's size k:
NB.     NB. eQf=. k ungl2k Qf
NB.     ungl2k=: ungl2@{.
NB.     NB. eQf=. k ung2lk Qf
NB.     ung2lk=: (ung2l@:({."1)~ -)~
NB.     NB. eQf=. k ung2rk Qf
NB.     ung2rk=: ung2r@:({."1)
NB.     NB. eQf=. k ungr2k Qf
NB.     ungr2k=: (ungr2@{.~ -)~

NB. non-tested yet
ungl3=: (((1 _ ,:~    -~&#) ,;.0 [) (((>:@] _1} *) +@-@{.)@[ , ~ larzrcfr) 0 ,.~ ])^:(#@[) (0 $~ 0 ,  -~/@$)
ung3l=: (((_ 1 ,:~ <:@- &c) ,;.0 [) (((>:@]  0} *)   -@{:)@[ ,.  larzlnbc) 0 ,   ])^:(c@[) (0 $~ 0 ,~ - /@$)
ung3r=: (((_ 1 ,:~    -~&c) ,;.0 [) (((>:@] _1} *)   -@{.)@[ ,.~ larzlnfc) 0 , ~ ])^:(c@[) (0 $~ 0 ,~ - /@$)

NB. verified
ungr3=: (((1 _ ,:~ 0 ,~ #@]) ,;.0 [) (((>:@]  (0 (=i.0:) ])} *) +@-@{:)@[ , ~ larzrcbr) ])^:(#@[) (0 (, $ [) c)

NB. simpler way to extract i-th row
ungr3_2=: (({~ #) (((>:@] (0 (=i.0:) ])} *) +@-@{:)@[ , ~ larzrcbr) ])^:(#@[) (0 (, $ [) c)

NB. mimic to ungr2_2
ungr3_3=: ([ larzrcbr ((, (1:`[`(0 $~ #@])}~ 1&(=i.1:)))~ _ _&{. :: (0 $~ 0 , #)))/@(0 , ~ |.)

NB. mimic to ungr2_3
ungr3_4=: ([ larzrcbr 1:`(_1 <@; 1 (=i.1:) [)`((, 0:)~  _ _&{. :: (0 $~ 0 , #))})/@(0 , ~ |.)

NB. ungr3_2 with idea of ungr3_3 and ungr3_4 of unified vectors processing
ungr3_5=: (({~ #) ([ larzrcbr (, (1:`[`(0 $~ #@])}~ 1&(=i.1:)))~) ])^:(#@[) (0 (, $ [) c)

NB. ungr3_5 with shorter e(i) in form: (0 ... 0 1) , where 1's position match with the same in (0,...,0,1,v(i),tau(i))
ungr3_6=: (({~ #) ([ larzrcbr (, 1 {.~ -@>:@(1&(=i.1:)))~) ])^:(#@[) (0 (, $ [) c)

NB. ungr3_6 with fork instead hook in length calc
ungr3_7=: (({~ #) ([ larzrcbr (, 1&([ {.~ _1 - =i.1:))~) ])^:(#@[) (0 (, $ [) c)

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
NB.           H(m-1:k)≡H(v(m-1:k),τ(m-1:k))=H(0,0)=I
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   example for m=4, n=5:       example for m=5, n=4:
NB.   (  l  v0 v0 v0 v0 τ0  )     (  l  v0 v0 v0 τ0  )
NB.   (  l  l  v1 v1 v1 τ1  )     (  l  l  v1 v1 τ1  )
NB.   (  l  l  l  v2 v2 τ2  )     (  l  l  l  v2 τ2  )
NB.   (  l  l  l  l  v3 τ3  )     (  l  l  l  l  τ3  )
NB.                               (  l  l  l  l  *   )
NB. where
NB.   l         - elements of m×k-matrix L, lower triangular
NB.   (1,vi,τi) - vector z(i) which defines an elementary
NB.               reflector H(i) (unit is not stored, v(i)
NB.               may be empty)
NB.   *         - any value, is not used
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. L * (Q) = L * Q
NB.   (((unmlqrn  trlpick        ) -: ((mp  unglq)~  trl        )) }:"1) LQf
NB.   NB. I = Q * (Q^H)
NB.   (( 0         idmat (<. , ])/)@(0 _1 + $) (-: clean) (unmlqrc unglq)) LQf
NB.
NB. Application:
NB. - simulate LAPACK's xUNGLQ(M,N,K) for M≠K:
NB.     NB. eQf=. k unglqk LQf
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
NB.           H(n-1:k)≡H(v(n-1:k),τ(n-1:k))=H(0,0)=I
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   example for m=4, n=5:       example for m=5, n=4:
NB.   (  *  τ1 τ2 τ3 τ4  )        (  τ0 τ1 τ2 τ3  )
NB.   (  l  l  v2 v3 v4  )        (  v0 v1 v2 v3  )
NB.   (  l  l  l  v3 v4  )        (  l  v1 v2 v3  )
NB.   (  l  l  l  l  v4  )        (  l  l  v2 v3  )
NB.   (  l  l  l  l  l   )        (  l  l  l  v3  )
NB.                               (  l  l  l  l   )
NB. where
NB.   l         - elements of k×n-matrix L, lower triangular
NB.   (1,vi,τi) - vector z(i) which defines an elementary
NB.               reflector H(i) (unit is not stored, v(i)
NB.               may be empty)
NB.   *         - any value, is not used
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. (Q) * L = Q * L
NB.   (((unmqlln (trlpick~ -~/@$)) -: ((mp~ ungql)~ (trl~ -~/@$))) }.  ) QfL
NB.   NB. I = (Q^H) * Q
NB.   (((0 <. -~/) idmat ([ , <.)/)@(_1 0 + $) (-: clean) (unmqllc ungql)) QfL
NB.
NB. Application:
NB. - simulate LAPACK's xUNGQL(M,N,K) for M≠K:
NB.     NB. eQf=. k ungqlk QfL
NB.     ungqlk=: (ungql@:({."1)~ -)~
NB.
NB. Notes:
NB. - implements LAPACK's DORGQL, ZUNGQL for N=K
NB. - straightforward O(k*m^3) code:
NB.   Q=. (-k) {."1 mp/ (idmat m) -"2 |. ({. Qf) * (* +)"0/~"1 |: }. Qf

ungql=: }.  @((((((GQNB ,~ _) ,:~ (GQNB - 1) +  -~&c) ];.0 [) (ung2l@[ ,.~ larfblnbc) ]) ({.  ~    GQNB  + #))^:(GQNB %~ -&c) ung2l@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ c))@(tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $)

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
NB.           H(k:n-1)≡H(v(k:n-1),τ(k:n-1))=H(0,0)=I
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   example for m=4, n=5:       example for m=5, n=4:
NB.   (  u  u  u  u  u  )         (  u  u  u  u   )
NB.   (  v0 u  u  u  u  )         (  v0 u  u  u   )
NB.   (  v0 v1 u  u  u  )         (  v0 v1 u  u   )
NB.   (  v0 v1 v2 u  u  )         (  v0 v1 v2 u   )
NB.   (  τ0 τ1 τ2 τ3 *  )         (  v0 v1 v2 v3  )
NB.                               (  τ0 τ1 τ2 τ3  )
NB. where
NB.   u         - elements of k×n-matrix R, upper triangular
NB.   (1,vi,τi) - vector z(i) which defines an elementary
NB.               reflector H(i) (unit is not stored, v(i)
NB.               may be empty)
NB.   *         - any value, is not used
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. (Q) * R = Q * R
NB.   (((unmqrln  trupick        ) -: ((mp~ ungqr)~  tru        )) }:  ) QfR
NB.   NB. I = (Q^H) * Q
NB.   (( 0         idmat ([ , <.)/)@(_1 0 + $) (-: clean) (unmqrlc ungqr)) QfR
NB.
NB. Application:
NB. - simulate LAPACK's xUNGQR(M,N,K) for M≠K:
NB.     NB. eQf=. k ungqrk QfR
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
NB.           H(k:m-1)≡H(v(k:m-1),τ(k:m-1))=H(0,0)=I
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   example for m=4, n=5:       example for m=5, n=4:
NB.   (  τ0 v0 u  u  u  u  )      (  *  u  u  u  u  )
NB.   (  τ1 v1 v1 u  u  u  )      (  τ1 u  u  u  u  )
NB.   (  τ2 v2 c2 v2 u  u  )      (  τ2 v2 u  u  u  )
NB.   (  τ3 v3 v3 v3 v3 u  )      (  τ3 v3 v3 u  u  )
NB.                               (  τ4 v4 v4 v4 u  )
NB. where
NB.   u         - elements of m×k-matrix R, upper triangular
NB.   (1,vi,τi) - vector z(i) which defines an elementary
NB.               reflector H(i) (unit is not stored, v(i)
NB.               may be empty)
NB.   *         - any value, is not used
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. R * (Q) = R * Q
NB.   (((unmrqrn (trupick~ -~/@$)) -: ((mp  ungrq)~ (tru~ -~/@$))) }."1) RQf
NB.   NB. I = Q * (Q^H)
NB.   (((0 >. -~/) idmat (<. , ])/)@(0 _1 + $) (-: clean) (unmrqrc ungrq)) RQf
NB.
NB. Application:
NB. - simulate LAPACK's xUNGRQ(M,N,K) for M≠K:
NB.     NB. eQf=. k ungrqk RQf
NB.     ungrqk=: (ungrq@  {.   ~ -)~
NB.
NB. Notes:
NB. - implements LAPACK's DORGRQ, ZUNGRQ for M=K
NB. - straightforward O(k*m^3) code:
NB.   Q=. (-k) {. mp/ (idmat n) -"2 (+ {."1 Qf) * (* +)"0/~"1 + }."1 Qf

ungrq=: }."1@((((((GQNB ,  _) ,:~ (GQNB - 1) +  -~&#) ];.0 [) (ungr2@[ , ~ larfbrcbr) ]) ({."1~    GQNB  + c))^:(GQNB %~ -&#) ungr2@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ #))@(trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $)

NB. ---------------------------------------------------------
NB. ungrz
NB.
NB. Description:
NB.   Generate a matrix with orthonormal rows from output of
NB.   tzrzf
NB.
NB. Syntax:
NB.   Q=. ungrz RZf
NB. where##########################
NB.   RQf - m×(n+1)-matrix, the output of gerqf
NB.   Q   - k×n-matrix with orthonormal rows, which is
NB.         defined as the last k rows of a product of m
NB.         elementary reflectors of order n:
NB.           Q = Π{H(i)',i=0:m-1}
NB.         where
NB.           H(k:m-1)≡H(v(k:m-1),τ(k:m-1))=H(0,0)=I
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   example for m=4, n=5:       example for m=5, n=4:
NB.   (  τ0 v0 u  u  u  u  )      (  *  u  u  u  u  )
NB.   (  τ1 v1 v1 u  u  u  )      (  τ1 u  u  u  u  )
NB.   (  τ2 v2 c2 v2 u  u  )      (  τ2 v2 u  u  u  )
NB.   (  τ3 v3 v3 v3 v3 u  )      (  τ3 v3 v3 u  u  )
NB.                               (  τ4 v4 v4 v4 u  )
NB. where
NB.   u         - elements of m×k-matrix R, upper triangular
NB.   (1,vi,τi) - vector z(i) which defines an elementary
NB.               reflector H(i) (unit is not stored, v(i)
NB.               may be empty)
NB.   *         - any value, is not used
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. R * (Q) = R * Q
NB.   (((unmrqrn (trupick~ -~/@$)) -: ((mp  ungrq)~ (tru~ -~/@$))) }."1) RQf
NB.   NB. I = Q * (Q^H)
NB.   (((0 >. -~/) idmat (<. , ])/)@(0 _1 + $) (-: clean) (unmrqrc ungrq)) RQf
NB.
NB. Application:
NB. - change eQ's size k:
NB.     NB. eQf=. k ungrzk RZf
NB.     ungrzk=: ungrq@{.

unglz=: }."1@((((((GQNB ,  _) ,:~               -~&#) ];.0 [) (ungl3@[ , ~ larzbrcfr) ]) ({."1~    GQNB  + c))^:(GQNB %~ -&#) ungl3@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ #))@(idmat@[`(a: <@; dhs2lios@(_1 , [))`]}~ #)
ungzl=: }:  @((((((GQNB ,~ _) ,:~  GQNB      -~ - &c) ];.0 [) (ung3l@[ ,.  larzblnbc) ]) ({.  ~ (- GQNB) - #))^:(GQNB %~ -&c) ung3l@(}.~ 2 #    GQNB  * 0 >. GQNB >.@%~ GQNX -~ c))@(idmat@[`(       dhs2lios@( 0 , [))`]}~ c)
ungzr=: }.  @((((((GQNB ,~ _) ,:~               -~&c) ];.0 [) (ung3r@[ ,.~ larzblnfc) ]) ({.  ~    GQNB  + #))^:(GQNB %~ -&c) ung3r@(}.~ 2 # (- GQNB) * 0 >. GQNB >.@%~ GQNX -~ c))@(idmat@[`(       dhs2lios@(_1 , [))`]}~ c)
ungrz=: }:"1@((((((GQNB ,  _) ,:~               - &#) ];.0 [) (ungr3@[ ,   larzbrcbr) ]) ({."1~ (- GQNB) - c))^:(GQNB %~ -&#) ungr3@(}.~ 2 #    GQNB  * 0 >. GQNB >.@%~ GQNX -~ #))@(idmat@[`(a: <@; dhs2lios@( 0 , [))`]}~ #)

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
NB.   Test Q generation algorithms by general matrix given
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
NB. testunghr
NB.
NB. Description:
NB.   Test Q generation algorithms by square matrix given
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

testgq=: 1 : 'EMPTY_mt_ [ (testunghr_mt_^:(=/@$) [ testungq_mt_)@u'
