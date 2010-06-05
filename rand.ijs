NB. Random arrays
NB.
NB. trl1mat  Random unit lower triangular matrix
NB. trlmat   Random lower triangular matrix
NB. tru1mat  Random unit upper triangular matrix
NB. trumat   Random upper triangular matrix
NB. gemat    Random general matrix
NB. hemat    Random Hermitian (symmetric) matrix
NB. unmat    Random unitary (orthogonal) matrix
NB. dimat    Random diagonalizable matrix
NB. pomat    Random Hermitian (symmetric) positive defined
NB.          matrix
NB. ptmat    Random Hermitian (symmetric) positive defined
NB.          tridiagonal matrix
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Misc.

vu2y=: 2 : 'v @ u @ (2 & #)'  NB. v(u(y,y))

NB. ---------------------------------------------------------
NB. randu
NB.
NB. Description:
NB.   Uniform distribution U(a,b) with support (a,b)
NB.
NB. Syntax:
NB.   S=. [supp] randu sh
NB. where
NB.   sh   - n-vector, shape of S
NB.   supp - optional 2-vector (a,b), a<b, support of
NB.          distribution, open interval, default is (0 1)
NB.   S    - sh-array of values s ~ U(a,b)
NB.   n    ≥ 0
NB.
NB. Application:
NB. - verb to generate complex matrix S with elements s
NB.   having Re(s) ~ U(0,1) and Im(s) ~ U(_1,1):
NB.     randuc=: randu j. (_1 1 & randu)

randu=: (? @ $ 0:) :((p.~ (-~/\))~ $:)

NB. ---------------------------------------------------------
NB. rande
NB.
NB. Description:
NB.   Exponential distribution E(μ) with mean μ
NB.
NB. Syntax:
NB.   S=. [μ] rande sh
NB. where
NB.   sh - n-vector, shape of S
NB.   μ  > 0, optional mean, default is 1
NB.   S  - sh-array of values s ~ E(μ)
NB.   n  ≥ 0
NB.
NB. Formula:
NB.   e ← -μ*log(1-u)
NB. where
NB.   u ~ U(0,1)
NB.   e ~ E(μ)
NB.
NB. Notes:
NB. - (randu) is used instead of (1-randu) since they are
NB.   equivalent under following assumptions:
NB.   - randu's support is open interval (0,1)
NB.   - randu's values aren't used outside elsewhere

rande=: (- @ ^. @ randu) : (* $:)

NB. ---------------------------------------------------------
NB. randnf
NB. randnc
NB.
NB. Description:
NB.   Normal distribution N(μ,σ²) of float (complex) numbers
NB.   with mean μ and variance σ²
NB.
NB. Syntax:
NB.   S=. [par] randnf sh
NB.   S=. [par] randnc sh
NB. where
NB.   sh  - n-vector, shape of S
NB.   par - optional 2-vector (μ,σ), σ>0, parameters of
NB.         distribution, default is (0 1)
NB.   S   - sh-array of values s ~ N(μ,σ²)
NB.   n   ≥ 0
NB.
NB. Formula:
NB.   nf1 ← sqrt(e)*cos(u)         NB. Box-Muller
NB.   nf2 ← sqrt(e)*sin(u)         NB.   algorithm
NB.   nc  ← (nf1 + i*nf2)/sqrt(2)
NB.   n3  ← σ*n + μ
NB. where
NB.   u            ~ U(0,2*π)
NB.   e            ~ E(2)
NB.   nf1,nf2,nc,n ~ N(0,1)
NB.   n3           ~ N(μ,σ²)
NB.
NB. TODO:
NB. - currently randnf uses (2*N) random values to make N
NB.   random values

randnf=: ((%: @ (2 & rande)) * (  2 o. (0 2p1 & randu))) : (p. $:)
randnc=: ((%: @      rande ) r. (0 2p1 & randu)) : (p. $:)

NB. ---------------------------------------------------------
NB. randtnf
NB.
NB. Description:
NB.   Truncated normal distribution TN(μ,σ²,a,b) of float
NB.   numbers in range [a,b] with mean μ and variance σ²
NB.
NB. Syntax:
NB.   S=. [par] randtnf sh
NB. where
NB.   sh  - n-vector, shape of S
NB.   par - optional 4-vector (μ,σ,a,b), σ>0, a≤b, (μ,σ) are
NB.         normal distribution parameters, [a,b] is s range,
NB.         default is (0 1 -∞ +∞)
NB.   S   - sh-array of values s ~ N(μ,σ²), a ≤ s ≤ b
NB.   n   ≥ 0
NB.
NB. Application:
NB. - left-side truncated normal distribution TN(1,2²,3,+∞):
NB.   randltnf=: (1 2 3 _) & randtnf

randtnf=: (0 1 __ _ & $:) :(4 : 0)
  'mu sigma a b'=. x
  mrandnf=. (mu,sigma) & randnf
  NB. replace out-bounded elements recursively
  rober=. (I. @ ((a & >) +. (b & <))) ((($: @ mrandnf @ # @ [)`[`] }) ^: (0 < (#@[))) ]
  ($ (rober @ mrandnf @ (*/))) y
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:          Syntax:
NB. trl1mat        L1=. randx trl1mat n
NB. trlmat         L=.  randx trlmat  n
NB. tru1mat        U1=. randx tru1mat n
NB. trumat         U=.  randx trumat  n
NB.
NB. Description:
NB.   Adv. to make verb to generate random square triangular
NB.   matrix
NB. where
NB.   n     ≥ 0, size of matrix
NB.   randx - monadic verb to generate random y-matrix (shape
NB.           is taken from y)
NB.   L1    - n×n-matrix, random unit lower triangular
NB.   L     - n×n-matrix, random lower triangular
NB.   U1    - n×n-matrix, random unit upper triangular
NB.   U     - n×n-matrix, random upper triangular
NB.
NB. Algorithm for basic adverb trlmat:
NB.   In: n randx
NB.   Out: L
NB.   0) count non-zero elements in lower triangle:
NB.        nz := n*(n+1)/2
NB.   1) generate random nz-vector:
NB.        vz=. randx nz
NB.   2) prepare fret:
NB.      2.0) generate zero nz-vector:
NB.             f=. nz $ 0
NB.      2.1) generate lIOS for start of an interval marks:
NB.             lios=. +/\ i. nz
NB.      2.2) write in marks:
NB.             f=. 1 lios } f
NB.   3) cut vz on pieces of length 1 2 3 ..., stacking them
NB.      into lower triangular matrix:
NB.        L=. f ];.1 vz
NB.
NB. Application:
NB. - generate complex lower triangular 4×4-matrix with real
NB.   and imagine parts both having mantissa with uniform
NB.   distribution U(_1,1) and exponent with truncated normal
NB.   distribution TN(0,16,_6,4):
NB.     L=. _1 1 0 16 _6 4 & (gemat j. gemat) trmat 4

trlmat=:  1 : '1&([`(+/\@i.@])`(0 $~ -:@(* >:)@])}) ];.1 u@-:@(* >:)'
trl1mat=: 1 : '(1;a:) & setdiag_mt_ @ (u trlmat_mt_)'
trumat=:  1 : '|: @ (u trlmat_mt_)'
tru1mat=: 1 : '|: @ (u trl1mat_mt_)'

NB. ---------------------------------------------------------
NB. gemat
NB.
NB. Description:
NB.   Make random general matrix
NB.
NB. Syntax:
NB.   G=. [par] gemat sh
NB. where
NB.   sh  - n-vector, shape of G
NB.   par - optional 6-vector (ma,mb,μ,σ,ea,eb), σ>0, where
NB.         (ma,mb) are mantissa's uniform distribution
NB.         parameters, and (μ,σ,ea,eb) are exponent's
NB.         truncated normal distribution parameters, default
NB.         is: (_1 1 0 (FP_FLEN/2) FP_EMIN FP_EMAX)
NB.   G   - sh-matrix, random general
NB.   n   ≥ 0
NB.
NB. Formula:
NB.   g ← mantissa * 2 ^ exponent
NB. where
NB.   mantissa ~ U(ma,mb)
NB.   exponent ~ TN(μ,σ,ea,eb)
NB.
NB. Application:
NB. - generate complex 4×4-matrix with real and imagine parts
NB.   both having mantissa with uniform distribution U(_1,1)
NB.   and exponent with truncated normal distribution
NB.   TN(0,16,_6,4):
NB.     G=. _1 1 0 16 _6 4 (gemat j. gemat) 4 4
NB.
NB. Notes:
NB. - default par provides about 95% of g numbers falls into
NB.   the range [-1/FP_EPS,-FP_EPS]U[FP_EPS,1/FP_EPS]
NB.   ("68-95-99.7 rule")

gemat=: ((_1 1 0 , (-: FP_FLEN) , FP_EMIN , FP_EMAX) & $:) :(((randu~ 2&{.) (* 2&^) (randtnf~ 2&}.))~)

NB. ---------------------------------------------------------
NB. hemat
NB.
NB. Description:
NB.   Adv. to make verb to make random Hermitian (symmetric)
NB.   matrix
NB.
NB. Syntax:
NB.   H=. gemat hemat n
NB. where
NB.   n     ≥ 0, size of matrix H
NB.   gemat - monadic verb to generate random y-matrix (shape
NB.           is taken from y)
NB.   H     - n×n-matrix, random Hermitian (symmetric)
NB.
NB. Application:
NB. - generate symmetric 4×4-matrix with certain distribution
NB.   law
NB.     H=. (_1 1 0 16 _6 4 & gemat) hemat 4
NB. - generate Hermitian 4×4-matrix with certain distribution
NB.   law
NB.     H=. (_1 1 0 16 _6 4 & (gemat j. gemat)) hemat 4

hemat=: vu2y (+ ct_mt_)

NB. ---------------------------------------------------------
NB. unmat
NB.
NB. Description:
NB.   Adv. to make verb to make a random unitary (orthogonal)
NB.   matrix with distribution given by Haar measure
NB.
NB. Syntax:
NB.   U=. randn unmat n
NB. where
NB.   n     ≥ 0, size of matrix U
NB.   randn - monadic verb to generate random y-matrix (shape
NB.           is taken from y) with elements distributed as
NB.           N(0,1), it is either randnf or randnc
NB.   U     - n×n-matrix, random unitary (orthogonal)
NB.
NB. Formula:
NB.   Z ← randn(n,n)
NB.   (Q,R) ← QR(Z)
NB.   d ← diag(R)
NB.   Λ ← diagmat(d / |d|)
NB.   Q ← Q*Λ
NB. where
NB.   n - output matrix Q's size
NB.
NB. References:
NB. [1] Francesco Mezzadri (2007). How to generate random
NB.     matrices from the classical compact groups.
NB.     Notices of the AMS, Vol. 54 (2007), 592-604
NB.     http://arxiv.org/abs/math-ph/0609050v2

unmat=: ((* " 1) ((% |) @ ((<0 1) & |:))) & >/ @ geqrf_mt_ vu2y

NB. ---------------------------------------------------------
NB. dimat
NB.
NB. Description:
NB.   Conj. to make verb to make a random diagonalizable
NB.   square matrix
NB.
NB. Syntax:
NB.   D=. randx dimat randq n
NB. where
NB.   n     ≥ 0, size of matrix D
NB.   randq - monadic verb to generate random unitary
NB.           (orthogonal) (y,y)-matrix (size is taken from
NB.           y)
NB.   randx - monadic verb to generate random y-vector
NB.           (length is taken from y) of distinct or float
NB.           (non-complex) numbers
NB.   D     - n×n-matrix, random diagonalizable square
NB.
NB. Notes:
NB. - D will be Hermitian (symmetric) if randx produces float
NB.   (non-complex), possibly non-distinct numbers

dimat=: 2 : 'u (] mp (* ct_mt_)) v'

NB. ---------------------------------------------------------
NB. pomat
NB.
NB. Description:
NB.   Adv. to make verb to make a random Hermitian
NB.   (symmetric) positive defined matrix
NB.
NB. Syntax:
NB.   P=. randx pomat n
NB. where
NB.   n     ≥ 0, size of matrix P
NB.   randx - monadic verb to generate random non-singular
NB.           n×n-matrix
NB.   P     - n×n-matrix, random Hermitian (symmetric)
NB.           positive defined
NB.
NB. Application:
NB.   P=. (2p1 & rande) pomat 4

pomat=: vu2y (mp ct_mt_)

NB. ---------------------------------------------------------
NB. ptmat
NB.
NB. Description:
NB.   Adv. to make verb to make a random Hermitian
NB.   (symmetric) positive defined tridiagonal matrix
NB.
NB. Syntax:
NB.   T=. randx ptmat n
NB. where
NB.   n     ≥ 0, size of matrix T
NB.   randx - monadic verb to generate random vectors d and
NB.           e; is called as:
NB.             d=. (| @ (9 o. randx)) n
NB.             e=. randx (n-1)
NB.   d     - n-vector of positive numbers, the main diagonal
NB.           of T
NB.   e     - (n-1)-vector, the subdiagonal of T
NB.   T     - n×n-matrix, random Hermitian (symmetric)
NB.           positive defined tridiagonal
NB.
NB. Application:
NB.   T=. (randu r. rande) ptmat 4

ptmat=: 1 : '(u @ <:) (((+@[);1:) setdiag_mt_ (([;_1:) setdiag_mt_ ])) ((a:;~(|@(9 o. u))) setdiag_mt_ idmat_mt_)'
