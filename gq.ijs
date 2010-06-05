NB. gq.ijs
NB. Generate Q from LQ QL QR RQ output
NB.
NB. unglq  Generate a matrix with orthonormal rows from
NB.        output of gelqf
NB. ungql  Generate a matrix with orthonormal columns from
NB.        output of geqlf
NB. ungqr  Generate a matrix with orthonormal columns from
NB.        output of geqrf
NB. ungrq  Generate a matrix with orthonormal rows from
NB.        output of gerqf
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. Block size limits for blocked code
UNGQBSMIN=: 2
UNGQBSMAX=: 32

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
NB. ungqbs
NB.
NB. Description: Compute abjusted block size
NB. Syntax:      bs=. dshape ungqbs shapeQf
NB. Formula:     bs = max(UNGQBSMIN,min(UNGQBSMAX,k))

ungqbs=: (UNGQBSMIN >. UNGQBSMAX <. ungqk) M.

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. ungqiters
NB.
NB. Description: Compute number of iterations
NB. Syntax:      iters=. dshape ungqiters shapeQf
NB. Formula:     iters = ⌊k % bs⌋

ungqiters=: (ungqk (<. @ %) ungqbs) M.

NB. ---------------------------------------------------------
NB. Apply reflector H or H' from either the left or the right
NB. to the matrix C to update it, where
NB.   H := H(v,τ) := I - τ*v*v'
NB.
NB. Algorithm:
NB.   Input:
NB.     eC -: C with appended or stitched trash vector
NB.     z   - vector in form:
NB.             (0[0:any-1],1[any],v[any+1:n-1],τ)
NB.           for forward direction, or:
NB.             (τ,v[1:any-1],1[any],0[any+1:n])
NB.           for backward direction
NB.   Output:
NB.     eCupd - Cupd, being updated C, with appended or
NB.             stitched trash vector
NB.   Actions:
NB.     1) conjugate z (optionally): z=. + z
NB.     2) replace τ (at z tail if forward direction, or at z
NB.        head if backward one) by 0 to get ev, the
NB.        extension of vector v: ev=. 0 i} z
NB.     3) multiple C by ev:
NB.        a) if side is left: deltaC=. ev*ev'*C
NB.        b) if side is right: deltaC=. ev*ev'*C
NB.     4) multiply deltaC by optionally conjugated τ:
NB.          deltaC=. deltaC * 

larfr1fcc=: [ - ((({: @ ]) * (larfr (0 & (_1 })))) +)  NB. forward direction, H := H(conj(v),conj(τ))
larfl1bss=: [ - (({. @ ]) * (larfl (0 & (0 }))))       NB. backward direction, H := H(v,τ)
larfl1fss=: [ - (({: @ ]) * (larfl (0 & (_1 }))))      NB. forward direction, H := H(v,τ)
larfr1bcc=: [ - ((({. @ ]) * (larfr (0 & (0 })))) +)   NB. backward direction, H := H(conj(v),conj(τ))

NB. ---------------------------------------------------------
NB. ungl2step
NB. ung2lstep
NB. ung2rstep
NB. ungr2step
NB.
NB. Description:
NB.   Single step of non-blocked version of algorithms LQ QL
NB.   QR RQ
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
  y ((((- @ + @ {:) ((>: @ [) (io }) *) ]) @ ]) , larfr1fcc) (io { x)
)

ung2lstep=: 4 : 0
  io=. y (- & ({: @ $)) y
  y (larfl1bss ,. (((- @ {.) ((>: @ [) (io }) *) ]) @ ])) ((< a: ; io) { x)
)

ung2rstep=: 4 : 0
  io=. y (<: @ - & ({: @ $)) y
  y ((((- @ {:) ((>: @ [) (io }) *) ]) @ ]) ,. larfl1fss) ((< a: ; io) { x)
)

ungr2step=: 4 : 0
  io=. y (- & #) x
  y (larfr1bcc , (((- @ + @ {.) ((>: @ [) (io }) *) ]) @ ])) (io { x)
)

NB. ---------------------------------------------------------
NB. ungl2
NB. ung2l
NB. ung2r
NB. ungr2
NB.
NB. Description:
NB.   Non-blocked version of algorithms LQ QL QR RQ
NB.
NB. Syntax:
NB.   eQ=. ungxx Qf
NB. where
NB.   Qf  - unit triangular matrix, Q's factored form
NB.   eQ  - Q with appended or stitched trash vector
NB.   Q   - matrix with orthonormal rows or columns which is
NB.         defined as the product of k elementary reflectors
NB.
NB. Algorithm:
NB.   1) find iters, the number of iterations for ungxxstep
NB.   2) form eQ0
NB.   3) do iterations: eQ=. Qf (ungxxstep ^: iters) eQ0
NB.
NB. Notes:
NB. - following identity holds: 
NB.   eQ (-: & $) Qf

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
NB.   Single step of algorithms LQ QL QR RQ
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
NB.   4) aplly non-blocked version of algorithm to Qfi to
NB.      produce matrix eQi1part
NB.   5) assemble eCiupd and eQi1part, to produce eQi1

unglqstep=: 4 : 0
  'shx shy'=. x (,: & $) y
  bs=. 0 _1 ungqbs shx
  riosQfi=. ((shx - shy) - bs) ,: (bs + (0 (0 }) shy))
  sizeeC=. - ((0 , bs) + shy)
  (sizeeC {. y) ((ungl2 @ ]) , larfbrcfr) (riosQfi (] ;. 0) x)
)

ungqlstep=: 4 : 0
  'shx shy'=. x (,: & $) y
  bs=. _1 0 ungqbs shx
  riosQfi=. (0 (0 }) shy) ,: (bs + (0 (1 }) shy))
  sizeeC=. (bs , 0) + shy
  (sizeeC {. y) (larfblcbc ,. (ung2l @ ])) (riosQfi (] ;. 0) x)
)

ungqrstep=: 4 : 0
  'shx shy'=. x (,: & $) y
  bs=. _1 0 ungqbs shx
  riosQfi=. ((shx - shy) - bs) ,: (bs + (0 (1 }) shy))
  sizeeC=. - ((bs , 0) + shy)
  (sizeeC {. y) ((ung2r @ ]) ,. larfblnfc) (riosQfi (] ;. 0) x)
)

ungrqstep=: 4 : 0
  'shx shy'=. x (,: & $) y
  bs=. 0 _1 ungqbs shx
  riosQfi=. (0 (1 }) shy) ,: (bs + (0 (0 }) shy))
  sizeeC=. (0 , bs) + shy
  (sizeeC {. y) (larfbrcbr , (ungr2 @ ])) (riosQfi (] ;. 0) x)
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

unglq=: 3 : 0
  'k n1'=. sizeQf=. ((0 _1 & ungqk) , {:) ($ y)
  y=. tru1 sizeQf {. y
  ibs=. */ 'bs iters'=. 0 _1 (ungqbs , ungqiters) sizeQf
  sizeQf0=. - ((bs|k) , (n1-ibs))
  Q=. 0 _1 }. (y unglqstep ^: iters (ungl2 (sizeQf0 {. y)))
)

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
NB.   Q -: ungql (((m+1),(-n)) {. QfL)
NB.   I -: (mp~ ct) Q

ungql=: 3 : 0
  'm1 k'=. sizeQf=. ({. , (_1 0 & ungqk)) ($ y)
  y=. ((-~/ @ $) tru1 ]) sizeQf {. y
  ibs=. */ 'bs iters'=. _1 0 (ungqbs , ungqiters) sizeQf
  sizeQf0=. (m1-ibs) , (bs|k)
  Q=. 1 0 }. (y ungqlstep ^: iters (ung2l (sizeQf0 {. y)))
)

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
NB.   Q -: ungql (((m+1),n) {. QfR)
NB.   I -: (mp~ ct) Q

ungqr=: 3 : 0
  'm1 k'=. sizeQf=. ({. , (_1 0 & ungqk)) ($ y)
  y=. trl1 sizeQf {. y
  ibs=. */ 'bs iters'=. _1 0 (ungqbs , ungqiters) sizeQf
  sizeQf0=. - ((m1-ibs) , (bs|k))
  Q=. _1 0 }. (y ungqrstep ^: iters (ung2r (sizeQf0 {. y)))
)

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

ungrq=: 3 : 0
  'k n1'=. sizeQf=. ((0 _1 & ungqk) , {:) ($ y)
  y=. ((-~/ @ $) tru1 ]) sizeQf {. y
  ibs=. */ 'bs iters'=. 0 _1 (ungqbs , ungqiters) sizeQf
  sizeQf0=. (bs|k) , (n1-ibs)
  Q=. 0 1 }. (y ungrqstep ^: iters (ungr2 (sizeQf0 {. y)))
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. tungq
NB. Test Q generation algorithms: unglq unql ungqr ungrq by
NB. general matrix y
NB.
NB. Syntax: tungq A
NB. where A - general m×n-matrix

tungq=: 3 : 0
  NB. ### while no qf ###
  'fuzzym fuzzyn'=. fuzzysh=. $ y
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'gelqf geqlf geqrf gerqf'
  gelqf=. ((0&{::) (+ & (fuzzysh & {.)) ((tru0 @ (1&{::)) ,. (2&{::))) @ (2b1110 & gelqf_jlapack_)
  NB. ### /while no qf ###

  ('unglq' tmonad gelqf`]`(((%@mp&norm1) ct)@])`(_."_)`((((<: upddiag0)@(mp ct)) (%&norm1) ]) % (FP_EPS*#@]))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testgq
NB. Test Q generators on random matrices with shape y
NB. Syntax: mkge testgq (m,n)
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testgq 500 500

testgq=: 1 : 'EMPTY [ tungq @ u'
