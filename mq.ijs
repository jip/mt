NB. mq.ijs
NB. Multiply by Q from LQ QL QR RQ HRD output
NB.
NB. unmlqxx  Multiply a matrix by a matrix with orthonormal
NB.          rows as returned by gelqf
NB. unmqlxx  Multiply a matrix by a matrix with orthonormal
NB.          columns as returned by geqlf
NB. unmqrxx  Multiply a matrix by a matrix with orthonormal
NB.          columns as returned by geqrf
NB. unmrqxx  Multiply a matrix by a matrix with orthonormal
NB.          rows as returned by gerqf
NB. unmhrxx  Multiply a matrix by an unitary (orthogonal)
NB.          matrix with orthonormal columns as returned by
NB.          gehrd
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

UNMQBS=: 32  NB. block size limit

unml2ln=: ((1 1 }. [) (({. @ ]) , ($: }.)) ((larfrcfr {.)~)) ^: (0 < # @ ])    NB. FIXME larfrcfr swap args











NB. ---------------------------------------------------------
NB. Compute geometry for block versions of algorithms
NB.
NB. Symbols:
NB.   shapeQf - 2-vector (heightQf,widthQf), the shape of Qf
NB.   shapeQ  - 2-vector (heightQ,widthQ), the shape of Q
NB.   dshape  - 2-vector to convert shapeQf to shapeQ:
NB.               shapeQ=. dshape + shapeQf
NB.   Qf      - Q's factored form
NB.   Q       - matrix with orthonormal rows or columns which
NB.             is defined as the product of k elementary
NB.             reflectors
NB.
NB. Notes:
NB. - maked as memo, since repetitive calls are expected

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. ungqk
NB.
NB. Description: Find minimum in Q's shape
NB. Syntax:      k=. dshape ungqk shapeQf
NB. Formula:     k = min(heightQ,widthQ)

ungqk=: (<./ @: +) M.

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. ungqiters
NB.
NB. Description: Compute number of iterations
NB. Syntax:      iters=. dshape ungqiters shapeQf
NB. Formula:     iters = ⌊k % UNGQBS⌋

ungqiters=: (<. @ (UNGQBS %~ ungqk)) M.

NB. ---------------------------------------------------------
NB. ungl2step
NB. ung2lstep
NB. ung2rstep
NB. ungr2step
NB.
NB. Description:
NB.   Single step of non-blocked version of algorithms
NB.   ung{lq,ql,qr,rq}
NB.
NB. Syntax:
NB.   eQi1=. Qf ungxxstep eQi
NB. where
NB.   eQi  - Qi with appended or stitched trash vector
NB.   Qi   - Q(i), the matrix Q at i-th step
NB.   Qf   - Q's factored form
NB.   Q    - matrix with orthonormal rows or columns which is
NB.          defined as the product of k elementary
NB.          reflectors
NB.   eQi1 - Qi1 with appended or stitched trash vector
NB.   Qi1  - Q(i+1), the matrix Q after i-th step
NB.
NB. Algorithm:
NB.   1) extract vector z(i) from Qf: zi=. io { Qf
NB.   2) extract scalar τ from z(i), then negate and
NB.      optionally conjugate it
NB.   3) multiply z by new τ
NB.   4) replace io-th element of new z by incremented new τ
NB.   5) update eQi by new z, producing eQiupd
NB.   6) append or stitch new z to eQiupd, producing eQi1

ungl2step=: 4 : 0
  io=. x (<: @ - & #) y
  (io { x) ((((- @ + @ {:) ((>: @ [) (io }) *) ]) @ [) , larfrcfr) y
)

ung2lstep=: 4 : 0
  io=. y (- & ({: @ $)) x
  ((< a: ; io) { x) (larflnbc ,. (((- @ {.) ((>: @ [) (io }) *) ]) @ [)) y
)

ung2rstep=: 4 : 0
  io=. x (<: @ - & ({: @ $)) y
  ((< a: ; io) { x) ((((- @ {:) ((>: @ [) (io }) *) ]) @ [) ,. larflnfc) y
)

ungr2step=: 4 : 0
  io=. y (- & #) x
  (io { x) (larfrcbr , (((- @ + @ {.) ((>: @ [) (io }) *) [) @ ])) y
)

NB. ---------------------------------------------------------
NB. ungl2
NB. ung2l
NB. ung2r
NB. ungr2
NB.
NB. Description:
NB.   Non-blocked version of algorithms ung{lq,ql,qr,rq}
NB.
NB. Syntax:
NB.   eQ=. ungxx Qf
NB. where
NB.   Qf - unit triangular matrix, Q's factored form
NB.   eQ - Q with appended or stitched trash vector
NB.   Q  - matrix with orthonormal rows or columns which is
NB.        defined as the product of k elementary reflectors
NB.
NB. Algorithm:
NB.   1) find iters, the number of iterations for ungxxstep
NB.   2) form eQ0
NB.   3) do iterations: eQ=. Qf (ungxxstep ^: iters) eQ0
NB.
NB. Notes:
NB. - eQ and Qf shapes are the same

ungl2=: (0 $~ (0 (0 }) $)) (ungl2step ^: (0 _1 ungqk ($ @ [)))~ ]  NB. Qf -: tru1 Qf
ung2l=: (0 $~ (0 (1 }) $)) (ung2lstep ^: (_1 0 ungqk ($ @ [)))~ ]  NB. Qf -: ((-~/ @ $) tru1 ]) Qf
ung2r=: (0 $~ (0 (1 }) $)) (ung2rstep ^: (_1 0 ungqk ($ @ [)))~ ]  NB. Qf -: trl1 Qf
ungr2=: (0 $~ (0 (0 }) $)) (ungr2step ^: (0 _1 ungqk ($ @ [)))~ ]  NB. Qf -: ((-~/@$) trl1 ]) Qf

NB. ---------------------------------------------------------
NB. unglqstep
NB. ungqlstep
NB. ungqrstep
NB. ungrqstep
NB.
NB. Description:
NB.   Single step of algorithms ung{lq,ql,qr,rq}
NB.
NB. Syntax:
NB.   eQi1=. Qf ungxxstep eQi
NB. where
NB.   eQi  - Qi with appended or stitched trash vector
NB.   Qi   - Q(i), the matrix Q at i-th step
NB.   Qf   - Q's factored form
NB.   Q    - matrix with orthonormal rows or columns which is
NB.          defined as the product of k elementary
NB.          reflectors
NB.   eQi1 - Qi1 with appended or stitched trash vector
NB.   Qi1  - Q(i+1), the matrix Q after i-th step
NB.
NB. Algorithm:
NB.   1) extend matrix eQi by zeros to form eCi - being
NB.      matrix to update Ci with appended or stitched trash
NB.      vector
NB.   2) extract current block - the matrix Qfi from Qf
NB.   3) supply eQi and Qfi to larfbxxxx, to produce eCiupd
NB.   4) apply non-blocked version of algorithm to Qfi to
NB.      produce matrix eQi1part
NB.   5) assemble eCiupd and eQi1part, to produce eQi1
NB.
NB. TODO:
NB. - tacit versions

unglqstep=: 4 : 0
  'shx shy'=. x (,: & $) y
  riosQfi=. ((shx - shy) - UNGQBS) ,: (UNGQBS + (0 (0 }) shy))
  sizeeC=. - ((0 , UNGQBS) + shy)
  (riosQfi (] ;. 0) x) ((ungl2 @ [) , larfbrcfr) (sizeeC {. y)
)

ungqlstep=: 4 : 0
  'shx shy'=. x (,: & $) y
  riosQfi=. (0 (0 }) shy) ,: (UNGQBS + (0 (1 }) shy))
  sizeeC=. (UNGQBS , 0) + shy
  (riosQfi (] ;. 0) x) (larfblnbc ,. (ung2l @ [)) (sizeeC {. y)
)

ungqrstep=: 4 : 0
  'shx shy'=. x (,: & $) y
  riosQfi=. ((shx - shy) - UNGQBS) ,: (UNGQBS + (0 (1 }) shy))
  sizeeC=. - ((UNGQBS , 0) + shy)
  (riosQfi (] ;. 0) x) ((ung2r @ [) ,. larfblnfc) (sizeeC {. y)
)

ungrqstep=: 4 : 0
  'shx shy'=. x (,: & $) y
  riosQfi=. (0 (1 }) shy) ,: (UNGQBS + (0 (0 }) shy))
  sizeeC=. (0 , UNGQBS) + shy
  (riosQfi (] ;. 0) x) (larfbrcbr , (ungr2 @ [)) (sizeeC {. y)
)

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
NB.         defined as the product of k elementary reflectors
NB.         of order n: Q = H(k-1)' * ... * H(1)' * H(0)'
NB.   k   = min(m,n)
NB.
NB. Algorithm:
NB.   1) find iters, the number of iterations for unglqstep
NB.   2) form eQ0
NB.   3) do iterations: eQ=. Qf (unglqstep ^: iters) eQ0
NB.   4) cut off last column from eQ to produce Q
NB.
NB. If:
NB.   'm n'=. 0 _1 + $ LQf
NB.   k=. m <. n
NB.   Q=. unglq LQf
NB. then
NB.   Q -: unglq (k {. LQf)
NB.   I -: (mp ct) Q
NB.
NB. Notes:
NB. - emulates LAPACK's xUNGLQ with following difference: K
NB.   parameter (amount of leading vectors from LQf to form
NB.   Q) is assumed K=k; to emulate case K<k the last (k-K)
NB.   elements in τ[0:k-1] must be zeroed

unglq=: 0 _1 }. (unglqstep ^: ((0 _1 ungqiters $ @ [)`(ungl2 @ (}.~ (2 $ UNGQBS * 0 _1 ungqiters $)))))~ @ tru1 @ ((0 _1 ungqk $) {. ])

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
NB.         defined as the product of k elementary reflectors
NB.         of order m: Q = H(k-1) * ... * H(1) * H(0)
NB.   k   = min(m,n)
NB.
NB. Algorithm:
NB.   1) find iters, the number of iterations for ungqlstep
NB.   2) form eQ0
NB.   3) do iterations: eQ=. Qf (ungqlstep ^: iters) eQ0
NB.   4) cut off first row from eQ to produce Q
NB.
NB. If:
NB.   'm n'=. _1 0 + $ QfL
NB.   k=. m <. n
NB.   Q=. ungql QfL
NB. then
NB.   Q -: ungql (((m+1),(-k)) {. QfL)
NB.   I -: (mp~ ct) Q
NB.
NB. Notes:
NB. - emulates LAPACK's xUNGQL with following difference: K
NB.   parameter (amount of concluding vectors from QfL to
NB.   form Q) is assumed K=k; to emulate case K<k the first
NB.   (k-K) elements in τ[0:k-1] must be zeroed

ungql=: 1 0 }. (ungqlstep ^: ((_1 0 ungqiters $ @ [)`(ung2l @ (}.~ (2 $ (- UNGQBS) * _1 0 ungqiters $)))))~ @ ((-~/ @ $) tru1 ]) @ ((_ , (- @ (_1 0 ungqk $))) {. ])

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
NB.         defined as the product of k elementary reflectors
NB.         of order m: Q = H(0) * H(1) * ... * H(k-1)
NB.   k   = min(m,n)
NB.
NB. Algorithm:
NB.   1) find iters, the number of iterations for ungqrstep
NB.   2) form eQ0
NB.   3) do iterations: eQ=. Qf (ungqrstep ^: iters) eQ0
NB.   4) cut off last row from eQ to produce Q
NB.
NB. If:
NB.   'm n'=. _1 0 + $ QfR
NB.   k=. m <. n
NB.   Q=. ungqr QfR
NB. then
NB.   Q -: ungql (((m+1),k) {. QfR)
NB.   I -: (mp~ ct) Q
NB.
NB. Notes:
NB. - emulates LAPACK's xUNGQR with following difference: K
NB.   parameter (amount of leading vectors from QfR to form
NB.   Q) is assumed K=k; to emulate case K<k the last (k-K)
NB.   elements in τ[0:k-1] must be zeroed

ungqr=: _1 0 }. (ungqrstep ^: ((_1 0 ungqiters $ @ [)`(ung2r @ (}.~ (2 $ UNGQBS * _1 0 ungqiters $)))))~ @ trl1 @ ((_ , (_1 0 ungqk $)) {. ])

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
NB.         defined as the product of k elementary reflectors
NB.         of order n: Q = H(0)' * H(1)' * ... * H(k-1)'
NB.   k   = min(m,n)
NB.
NB. Algorithm:
NB.   1) find iters, the number of iterations for ungrqstep
NB.   2) form eQ0
NB.   3) do iterations: eQ=. Qf (ungrqstep ^: iters) eQ0
NB.   4) cut off first column from eQ to produce Q
NB.
NB. If:
NB.   'm n'=. 0 _1 + $ RQf
NB.   k=. m <. n
NB.   Q=. ungrq RQf
NB. then
NB.   Q -: ungrq (((-k),(n+1)) {. RQf)
NB.   I -: (mp ct) Q
NB.
NB. Notes:
NB. - emulates LAPACK's xUNGRQ with following difference: K
NB.   parameter (amount of concluding vectors from RQf to
NB.   form Q) is assumed K=k; to emulate case K<k the first
NB.   (k-K) elements in τ[0:k-1] must be zeroed

ungrq=: 0 1 }. (ungrqstep ^: ((0 _1 ungqiters $ @ [)`(ungr2 @ (}.~ (2 $ (- UNGQBS) * 0 _1 ungqiters $)))))~ @ ((-~/ @ $) trl1 ]) @ ((- @ (0 _1 ungqk $)) {. ])

NB. ---------------------------------------------------------
NB. unghr
NB.
NB. Description:
NB.   Generate an unitary (orthogonal) matrix Q which is
NB.   defined as the product of ({:fs) elementary reflectors
NB.   of order n, as returned by gehrd:
NB.     Q = H(f) H(f+1) ... H(f+s-1)
NB.
NB. Syntax:
NB.   Q=. unghr HQf ; fs
NB. where
NB.   HQf - max(n,f+s+1)×n-matrix with packed H and Qf (see
NB.         gehrd)
NB.   fs  - 2-vector of integers (f,s) 'from' and 'size',
NB.         defines submatrix A11 position in matrix A (see
NB.         gebalp)
NB.   Q   -
NB.
NB. TODO: dyad

unghr=: 3 : 0
  riosQf=. (2 2 $ 1 1 0 _1) + ,.~ 1 {:: y
  t=. l=. {. 1 {:: y
  b=. r=. (# 0 {:: y) - +/ 1 {:: y
  (1 0 _1,l) setdiag (1 0 0,t) setdiag ((t,l),:(b,r)) uncut riosQf (ungqr ;. 0) (0 {:: y)
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. tungq
NB. Test Q generation algorithms: unglq unql ungqr ungrq by
NB. general matrix
NB.
NB. Syntax: tungq A
NB. where A - general m×n-matrix

tungq=: 3 : 0
  ('unglq' tmonad (gelqf`]`(%@((*&norm1) ct)@])`(_."_)`(((norm1@(<: upddiag0)@(mp  ct)) % (FP_EPS*{:@$))@]))) y  NB. berr := ||Q*Q'-I||/ε*n
  ('ungql' tmonad (geqlf`]`(%@((*&norm1) ct)@])`(_."_)`(((norm1@(<: upddiag0)@(mp~ ct)) % (FP_EPS*#   ))@]))) y  NB. berr := ||Q'*Q-I||/ε*m
  ('ungqr' tmonad (geqrf`]`(%@((*&norm1) ct)@])`(_."_)`(((norm1@(<: upddiag0)@(mp~ ct)) % (FP_EPS*#   ))@]))) y  NB. berr := ||Q'*Q-I||/ε*m
  ('ungrq' tmonad (gerqf`]`(%@((*&norm1) ct)@])`(_."_)`(((norm1@(<: upddiag0)@(mp  ct)) % (FP_EPS*{:@$))@]))) y  NB. berr := ||Q*Q'-I||/ε*n
  EMPTY
)

NB. ---------------------------------------------------------
NB. testgq
NB. Adv. to test Q generators with random matrices of shape y
NB. Syntax: mkge testgq (m,n)
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testgq 150 100

testgq=: 1 : 'EMPTY [ tungq @ u'
