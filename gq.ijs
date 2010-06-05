NB. gq.ijs
NB. Generate Q from LQ QL QR RQ HRD output
NB.
NB. unglq   Generate a matrix with orthonormal rows from
NB.         output of gelqf
NB. ungql   Generate a matrix with orthonormal columns from
NB.         output of geqlf
NB. ungqr   Generate a matrix with orthonormal columns from
NB.         output of geqrf
NB. ungrq   Generate a matrix with orthonormal rows from
NB.         output of gerqf
NB. unghrl  Generate an unitary (orthogonal) matrix which is
NB.         defined as the product of elementary reflectors as
NB.         returned by gehrdl
NB. unghru  Generate an unitary (orthogonal) matrix which is
NB.         defined as the product of elementary reflectors as
NB.         returned by gehrdu
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

UNGBS=: 32   NB. block size limit
UNGNX=: 128  NB. crossover point

NB. ---------------------------------------------------------
NB. ungi
NB.
NB. Description: Number of iterations
NB. Syntax:      iters=. ungi k
NB. where        k = min(rows,columns)
NB. Formula:     iters = max(0,⌊(k+BS-NX-1)/BS⌋)
NB. Notes:       is memo, since repetitive calls are expected

ungi=: (0 >. <.@(UNGBS %~ (_1+UNGBS-UNGNX)&+))M.

NB. ---------------------------------------------------------
NB. ungb
NB.
NB. Description: Size of submatrix processed by blocked algo
NB. Syntax:      size=. ungb k
NB. where        k = min(rows,columns)
NB. Formula:     size = min(k,BS*iters)
NB. Notes:       is memo, since repetitive calls are expected

ungb=: (<. (UNGBS * ungi))M.

NB. ---------------------------------------------------------
NB. Description:
NB.   Single step of non-blocked version of algorithms
NB.
NB. Syntax:
NB.   eQi1=. (s;Qf) ungxxstep eQi
NB. where
NB.   eQi  - Qi augmented by trash vector
NB.   Qi   - Q(i), the matrix Q after i-th step and before
NB.          (i+1)-th one
NB.   Qf   - k×(n+1)-matrix for ungl2 and ungr2, or
NB.          (m+1)×k-matrix for ung2l and ung2r, unit
NB.          triangular, the Q's factored form
NB.   s    - integer in range either [0,n] for ungl2 and
NB.          ungr2 or [0,m] for ung2l and ung2r
NB.   eQi1 - Qi1 augmented by modified trash vector
NB.   Qi1  - Q(i+1), the matrix Q after (i+1)-th step
NB.   Q    - s×n-matrix for ungl2 and ungr2, or m×s-matrix
NB.          for ung2l and ung2r, with orthonormal
NB.          rows/columns which is defined as the first/last
NB.          rows/columns of a product of s elementary
NB.          reflectors (see corresp. ung{lq,ql,qr,rq})
NB.   k    = min(m,n)
NB.   i    - integer in range [0,k]
NB.
NB. Algorithm:
NB.   1) augment eQ(i) by zero vector
NB.   2) find rIOS of z(i) in Qf which describes a current
NB.      elementary reflector H(i)≡H(v(i),τ(i))
NB.   3) extract vector z(i) from Qf and ravel it
NB.   4) reconstruct an elementary reflector H(i) from z(i)
NB.   5) update augmented eQ(i) by H(i)
NB.   6) reconstruct vector subQ(i) from z(i)
NB.   7) merge subQ(i) with updated eQ(i) to produce eQ(i+1)
NB.
NB. Notes:
NB. - ung{l2,2l,2r,r2} and corresp. ung{lq,ql,qr,rq} are
NB.   topologic equivalents

ungl2step=: ((1 _ ,:~ (0 {:: [) <:@-  #@]) (,;.0) 1 {:: [) (((>:@]  0} *) +@-@{:)@[ ,   larfrcfr) 0 ,.  ]
ung2lstep=: ((_ 1 ,:~ (0 {:: [)    -~ c@]) (,;.0) 1 {:: [) (((>:@] _1} *)   -@{.)@[ ,.~ larflnbc) 0 , ~ ]
ung2rstep=: ((_ 1 ,:~ (0 {:: [) <:@-  c@]) (,;.0) 1 {:: [) (((>:@]  0} *)   -@{:)@[ ,.  larflnfc) 0 ,   ]
ungr2step=: ((1 _ ,:~ (0 {:: [)    -~ #@]) (,;.0) 1 {:: [) (((>:@] _1} *) +@-@{.)@[ , ~ larfrcbr) 0 ,.~ ]

NB. ---------------------------------------------------------
NB. Description:
NB.   Single step of algorithms
NB.
NB. Syntax:
NB.   eQi1=.    Qf  ungxxstep eQi         NB. unglq,ungqr
NB.   eQi1=. (s;Qf) ungxxstep eQi         NB. ungql,ungrq
NB. where
NB.   eQi  - Qi augmented by trash vector
NB.   Qi   - Q(i), the matrix Q at i-th step and before
NB.          (i+1)-th one
NB.   Qf   - k×(n+1)-matrix for unglq and ungrq, or
NB.          (m+1)×k-matrix for ungql and ungqr, unit
NB.          triangular, the Q's factored form
NB.   s    - integer in range either [0,n] for unglq and
NB.          ungrq or [0,m] for ungql and ungqr
NB.   eQi1 - Qi1 augmented by modified trash vector
NB.   Qi1  - Q(i+1), the matrix Q after (i+1)-th step
NB.   Q    - s×n-matrix for unglq and ungrq, or m×s-matrix
NB.          for ungql and ungqr, with orthonormal
NB.          rows/columns which is defined as the first/last
NB.          rows/columns of a product of s elementary
NB.          reflectors (see corresp. ung{lq,ql,qr,rq})
NB.   k    = min(m,n)
NB.   i    - integer in range [0,ungi(k)]
NB.
NB. Algorithm:
NB.   1) augment eQ(i) by zero block
NB.   2) find rIOS of Z(i) in Qf which describes a current
NB.      block reflector H(i)≡H(V(i),Τ(i))
NB.   3) extract matrix Z(i) from Qf
NB.   4) reconstruct a block reflector H(i) from Z(i)
NB.   5) update augmented eQ(i) by H(i)
NB.   6) reconstruct matrix subQ(i) from Z(i) by non-blocked
NB.      version of algorithm
NB.   7) merge subQ(i) with updated eQ(i) to produce eQ(i+1)
NB.
NB. Notes:
NB. - ung{l2,2l,2r,r2} and corresp. ung{lq,ql,qr,rq} are
NB.   topologic equivalents

unglqstep=: ((((UNGBS, _),:~(-&c)) (] ;. 0) [)(((ungl2~ #) @ [) ,  larfbrcfr) ]) ({.~ (_, ((-UNGBS)-c)))
ungqrstep=: ((((UNGBS,~_),:~(-&#)) (] ;. 0) [)(((ung2r~ c) @ [) ,. larfblnfc) ]) ({.~ (_,~((-UNGBS)-#)))

ungqlstep=: ((((UNGBS,~_),:~(0 {:: [) ((UNGBS-1)+-~) c@]) (] ;. 0) 1 {:: [) (((ung2l~ c) @ [) ,.~ larfblnbc) ]) ({.~ (_,~(  UNGBS +#)))
ungrqstep=: ((((UNGBS, _),:~(0 {:: [) ((UNGBS-1)+-~) #@]) (] ;. 0) 1 {:: [) (((ungr2~ #) @ [) , ~ larfbrcbr) ]) ({.~ (_, (  UNGBS +c)))

NB. ---------------------------------------------------------
NB. Description:
NB.   Non-blocked version of algorithms
NB.
NB. Syntax:
NB.   eQ=. s ungxx Qf
NB. where
NB.   Qf - k×(n+1)-matrix for ungl2 and ungr2, or
NB.        (m+1)×k-matrix for ung2l and ung2r, unit
NB.        triangular, the Q's factored form
NB.   s  - integer in range either [0,n] for ungl2 and ungr2
NB.        or [0,m] for ung2l and ung2r
NB.   eQ - Q augmented by trash vector
NB.   Q  - s×n-matrix for ungl2 and ungr2, or m×s-matrix for
NB.        ung2l and ung2r, with orthonormal rows/columns
NB.        which is defined as the first/last rows/columns of
NB.        a product of s elementary reflectors (see corresp.
NB.        ung{lq,ql,qr,rq})
NB.   k  = min(m,n)
NB.
NB. Algorithm:
NB.   1) form eQ(0) as unit matrix of proper size
NB.   2) do iterations: eQ=. (s;Qf) (ungxxstep ^: k) eQ0
NB.
NB. If:
NB.   'm n'=. $ A
NB.   k=. m <. n
NB. then (with appropriate comparison tolerance)
NB.   (] -: clean @ ((         trl   @( 0 _1&}.)) mp  (( 0 _1&}.)@(# ungl2 ])@         tru1   @(      k &{.)))@gelqf) A
NB.   (] -: clean @ ((((-~/@$) trl ])@( 1  0&}.)) mp~ (( 1  0&}.)@(c ung2l ])@((-~/@$) tru1 ])@((_,  -k)&{.)))@geqlf) A
NB.   (] -: clean @ ((         tru   @(_1  0&}.)) mp~ ((_1  0&}.)@(c ung2r ])@         trl1   @((_,   k)&{.)))@geqrf) A
NB.   (] -: clean @ ((((-~/@$) tru ])@( 0  1&}.)) mp  (( 0  1&}.)@(# ungr2 ])@((-~/@$) trl1 ])@((_,~ -k)&{.)))@gerqf) A

ungl2=: ungl2step^:(;`(#@])`( idmat        @((,  c)-#@])))
ung2l=: ung2lstep^:(;`(c@])`((idmat~ (-~/))@((,~ #)-c@])))
ung2r=: ung2rstep^:(;`(c@])`( idmat        @((,~ #)-c@])))
ungr2=: ungr2step^:(;`(#@])`((idmat~ (-~/))@((,  c)-#@])))

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
NB.   Q=. [s] unglq LQf
NB. where
NB.   LQf - m×(n+1)-matrix, the output of gelqf
NB.   Q   - s×n-matrix with orthonormal rows, which is
NB.         defined as the first s rows of a product of m
NB.         elementary reflectors of order n:
NB.           Q = Π{H(i)',i=m-1:0}
NB.         where
NB.           H(m-1:s)≡H(v(m-1:s),τ(m-1:s))=H(0,0)=I
NB.   s   - integer in range [0,n], default is k
NB.   k   = min(m,n)
NB.
NB. Algorithm:
NB.   1) extract Qf from LQf selecting only first min(s,k)
NB.      rows
NB.   2) make Qf the unit upper triangular
NB.   3) call ungl2, the non-blocked version of algorithm,
NB.      with arguments supplied: adjusted s, and Qf with
NB.      part intended for blocked algorithm excluded, to
NB.      produce eQ(0)
NB.   4) find iters, the number of iterations
NB.   5) do iterations:
NB.        eQ=.    Qf  (unglqstep ^: iters) eQ0
NB.   6) cut off last column from eQ to produce Q
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
NB. If:
NB.   Q=. unglq gelqf A
NB. then (with appropriate comparison tolerance)
NB.   (idmat @ ms @ $ -: clean @ (mp  ct)) Q
NB.
NB. Notes:
NB. - emulates LAPACK's xUNGLQ

unglq=: ($:~ ( 0 _1&(ms $))) : ( 0 _1 }. ([ (unglqstep^:(]`(ungi@#@])`((-(ungb@#)) ungl2 ((}.~ (2 # (  ungb@#)))@])))) ( tru1            @((_,~  (<. ( 0 _1&(ms $)))) {. ]))))

NB. ---------------------------------------------------------
NB. ungql
NB.
NB. Description:
NB.   Generate a matrix with orthonormal columns from output
NB.   of geqlf
NB.
NB. Syntax:
NB.   Q=. [s] ungql QfL
NB. where
NB.   QfL - (m+1)×n-matrix, the output of geqlf
NB.   Q   - m×s-matrix with orthonormal columns, which is
NB.         defined as the last s columns of a product of n
NB.         elementary reflectors of order m:
NB.           Q = Π{H(i),i=n-1:0}
NB.         where
NB.           H(n-1:k)≡H(v(n-1:k),τ(n-1:k))=H(0,0)=I
NB.   s   - integer in range [0,m], default is k
NB.   k   = min(m,n)
NB.
NB. Algorithm:
NB.   1) extract Qf from QfL selecting only last min(s,k)
NB.      columns
NB.   2) make Qf the unit upper triangular
NB.   3) call ung2l, the non-blocked version of algorithm,
NB.      with arguments supplied: adjusted s, and Qf with
NB.      part intended for blocked algorithm excluded, to
NB.      produce eQ(0)
NB.   4) find iters, the number of iterations
NB.   5) do iterations:
NB.        eQ=. (s;Qf) (ungqlstep ^: iters) eQ0
NB.   6) cut off first row from eQ to produce Q
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
NB. If:
NB.   Q=. ungql geqlf A
NB. then (with appropriate comparison tolerance)
NB.   (idmat @ ms @ $ -: clean @ (mp~ ct)) Q
NB.
NB. Notes:
NB. - emulates LAPACK's xUNGQL

ungql=: ($:~ (_1  0&(ms $))) : ( 1  0 }. ([ (ungqlstep^:(;`(ungi@c@])`((-(ungb@c)) ung2l ((}.~ (2 # (-@ungb@c)))@])))) ((tru1~ (-~/ @ $))@((_, -@(<. (_1  0&(ms $)))) {. ]))))

NB. ---------------------------------------------------------
NB. ungqr
NB.
NB. Description:
NB.   Generate a matrix with orthonormal columns from output
NB.   of geqrf
NB.
NB. Syntax:
NB.   Q=. [s] ungqr QfR
NB. where
NB.   QfR - (m+1)×n-matrix, the output of geqrf
NB.   Q   - m×s-matrix with orthonormal columns, which is
NB.         defined as the first s columns of a product of n
NB.         elementary reflectors of order m:
NB.           Q = Π{H(i),i=0:n-1}
NB.         where
NB.           H(k:n-1)≡H(v(k:n-1),τ(k:n-1))=H(0,0)=I
NB.   s   - integer in range [0,m], default is k
NB.   k   = min(m,n)
NB.
NB. Algorithm:
NB.   1) extract Qf from QfR selecting only first min(s,k)
NB.      columns
NB.   2) make Qf the unit lower triangular
NB.   3) call ung2r, the non-blocked version of algorithm,
NB.      with arguments supplied: adjusted s, and Qf with
NB.      part intended for blocked algorithm excluded, to
NB.      produce eQ(0)
NB.   4) find iters, the number of iterations
NB.   5) do iterations:
NB.        eQ=.    Qf  (ungqrstep ^: iters) eQ0
NB.   6) cut off last row from eQ to produce Q
NB.
NB. Storage layout:
NB.   example for m=4, n=5:       example for m=5, n=4:
NB.   (  u  u  u  u  u  )         (  u  u  u  u   )
NB.   (  v0 u  u  u  u  )         (  v0 u  u  u   )
NB.   (  v0 v1 u  u  u  )         (  v0 v1 u  u   )
NB.   (  v0 v1 v2 u  u  )         (  v0 v1 v2 u   )
NB.   (  τ0 τ1 τ2 τ3  *  )         (  v0 v1 v2 v3  )
NB.                               (  τ0 τ1 τ2 τ3  )
NB. where
NB.   u         - elements of k×n-matrix R, upper triangular
NB.   (1,vi,τi) - vector z(i) which defines an elementary
NB.               reflector H(i) (unit is not stored, v(i)
NB.               may be empty)
NB.   *         - any value, is not used
NB.
NB. If:
NB.   Q=. ungqr geqrf A
NB. then (with appropriate comparison tolerance)
NB.   (idmat @ ms @ $ -: clean @ (mp~ ct)) Q
NB.
NB. Notes:
NB. - emulates LAPACK's xUNGQR

ungqr=: ($:~ (_1  0&(ms $))) : (_1  0 }. ([ (ungqrstep^:(]`(ungi@c@])`((-(ungb@c)) ung2r ((}.~ (2 # (  ungb@c)))@])))) ( trl1            @((_,   (<. (_1  0&(ms $)))) {. ]))))

NB. ---------------------------------------------------------
NB. ungrq
NB.
NB. Description:
NB.   Generate a matrix with orthonormal rows from output of
NB.   gerqf
NB.
NB. Syntax:
NB.   Q=. [s] ungrq RQf
NB. where
NB.   RQf - m×(n+1)-matrix, the output of gerqf
NB.   Q   - s×n-matrix with orthonormal rows, which is
NB.         defined as the last s rows of a product of m
NB.         elementary reflectors of order n:
NB.           Q = Π{H(i)',i=0:m-1}
NB.         where
NB.           H(k:m-1)≡H(v(k:m-1),τ(k:m-1))=H(0,0)=I
NB.   s   - integer in range [0,n], default is k
NB.   k   = min(m,n)
NB.
NB. Algorithm:
NB.   1) extract Qf from RQf selecting only last min(s,k)
NB.      rows
NB.   2) make Qf the unit lower triangular
NB.   3) call ungr2, the non-blocked version of algorithm,
NB.      with arguments supplied: adjusted s, and Qf with
NB.      part intended for blocked algorithm excluded, to
NB.      produce eQ(0)
NB.   4) find iters, the number of iterations
NB.   5) do iterations:
NB.        eQ=. (s;Qf) (ungrqstep ^: iters) eQ0
NB.   6) cut off first column from eQ to produce Q
NB.
NB. Storage layout:
NB.   example for m=4, n=5:       example for m=5, n=4:
NB.   (  τ0 v0 u  u  u  u   )     (  *  u  u  u  u   )
NB.   (  τ1 v1 v1 u  u  u   )     (  τ1 u  u  u  u   )
NB.   (  τ2 v2 c2 v2 u  u   )     (  τ2 v2 u  u  u   )
NB.   (  τ3 v3 v3 v3 v3 u   )     (  τ3 v3 v3 u  u   )
NB.                               (  τ4 v4 v4 v4 u   )
NB. where
NB.   u         - elements of m×k-matrix R, upper triangular
NB.   (1,vi,τi) - vector z(i) which defines an elementary
NB.               reflector H(i) (unit is not stored, v(i)
NB.               may be empty)
NB.   *         - any value, is not used
NB.
NB. If:
NB.   Q=. ungrq gerqf A
NB. then (with appropriate comparison tolerance)
NB.   (idmat @ ms @ $ -: clean @ (mp ct)) Q
NB.
NB. Notes:
NB. - emulates LAPACK's xUNGRQ

ungrq=: ($:~ ( 0 _1&(ms $))) : ( 0  1 }. ([ (ungrqstep^:(;`(ungi@#@])`((-(ungb@#)) ungr2 ((}.~ (2 # (-@ungb@#)))@])))) ((trl1~ (-~/ @ $))@((_,~-@(<. ( 0 _1&(ms $)))) {. ]))))

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

unghrl=: unglq @ (|. !. 0)

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
NB. - emulates LAPACK's xUNGHR
NB. - instead of using f and s parameters, the following
NB.   product is really calculating:
NB.     Q = Π{H(i),i=0:n-1} ,
NB.   where
NB.     H(0:f-1) = H(f+s-1:n-1) = I .
NB.   This approach delivers excessive calculations in rare
NB.   case ((f>0) OR (f+s<n)).

unghru=: ungqr @ (0 _1 & (|. !. 0))

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. tgq
NB.
NB. Description:
NB.   Test Q generation algorithms by matrix given
NB.
NB. Syntax:
NB.   tgq A
NB. where
NB.   A - m×n-matrix

tgq=: 3 : 0
  cberr=. 2 : '((norm1@(<: upddiag)@(u ct)) % (FP_EPS * v)) @ ]'

  ('unglq' tmonad (gelqf`]`((norm1 con ct)@])`(_."_)`(mp  cberr c))) y  NB. berr := ||Q*Q'-I||/(ε*n)
  ('ungql' tmonad (geqlf`]`((norm1 con ct)@])`(_."_)`(mp~ cberr #))) y  NB. berr := ||Q'*Q-I||/(ε*m)
  ('ungqr' tmonad (geqrf`]`((norm1 con ct)@])`(_."_)`(mp~ cberr #))) y  NB. berr := ||Q'*Q-I||/(ε*m)
  ('ungrq' tmonad (gerqf`]`((norm1 con ct)@])`(_."_)`(mp  cberr c))) y  NB. berr := ||Q*Q'-I||/(ε*n)

  NB. following are tested iif A is square
  ('unghrl' tmonad ((gehrdl~ (0,#))`]`((norm1 con ct)@])`(_."_)`(mp  cberr c))) ^: (=/ @ $) y  NB. berr := ||Q*Q'-I||/(ε*n)
  ('unghru' tmonad ((gehrdu~ (0,#))`]`((norm1 con ct)@])`(_."_)`(mp~ cberr #))) ^: (=/ @ $) y  NB. berr := ||Q'*Q-I||/(ε*n)
  EMPTY
)

NB. ---------------------------------------------------------
NB. testgq
NB.
NB. Description:
NB.   Adv. to make verb to test Q generation algorithms by
NB.   matrix of generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testgq
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.            mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testgq 150 100

testgq=: 1 : 'EMPTY [ tgq @ u'
