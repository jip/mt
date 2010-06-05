NB. rand.ijs
NB. Random objects
NB.
NB. randu    Uniform distribution U(a,b)
NB. rande    Exponential distribution E(μ)
NB. randnf   Normal distribution N(μ,σ²) of float numbers
NB. randnc   Normal distribution N(μ,σ²) of complex numbers
NB. randtnf  Truncated normal distribution TN(μ,σ²,a,b) of
NB.          float numbers
NB.
NB. gemat    Random general matrix
NB. hemat    Random Hermitian (symmetric) matrix
NB. unmat    Random unitary (orthogonal) matrix
NB. dimat    Random diagonalizable matrix
NB. pomat    Random Hermitian (symmetric) positive defined
NB.          matrix
NB.
NB. XRef:
NB. - system: mp
NB. - mt: ct qr
NB.
NB. Copyright (C) 2009  Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

qr=: 128!:0                   NB. built-in QR factorization (TODO: implement)
vu2y=: 2 : 'v @ u @ (2 & $)'  NB. v(u(y,y))

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. randu                                               _ _ _
NB. Uniform distribution U(a,b) with support (a,b)
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
NB.   having Re(s) ~ U(0,1) and Im(s) ~ U(_1,1) :
NB.   randuc=: randu j. (_1 1 & randu)

randu=: (? @ $ 0:) :((p.~ (-~/\))~ $:)

NB. ---------------------------------------------------------
NB. rande                                               _ _ _
NB. Exponential distribution E(μ) with mean μ
NB.
NB. Formula:
NB.   e ← -μ*log(1-u)
NB. where
NB.   u ~ U(0,1)
NB.   e ~ E(μ)
NB.
NB. Syntax:
NB.   S=. [μ] rande sh
NB. where
NB.   sh - n-vector, shape of S
NB.   μ  > 0, optional mean, default is 1
NB.   S  - sh-array of values s ~ E(μ)
NB.   n  ≥ 0
NB.
NB. Notes:
NB. - (randu) is used instead of (1-randu) since they are
NB.   equivalent under following assumptions:
NB.   - randu's support is open interval (0,1)
NB.   - randu's values aren't used outside elsewhere

rande=: (- @ ^. @ randu) : (* $:)

NB. ---------------------------------------------------------
NB. randnf                                              _ _ _
NB. randnc                                              _ _ _
NB. Normal distribution N(μ,σ²) of float (complex) numbers
NB. with mean μ and variance σ²
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
NB. Notes:
NB. - randnf uses (2*N) random values to make N random values

randnf=: ((%: @ (2 & rande)) * (  2 o. (0 2p1 & randu))) : (p. $:)
randnc=: ((%: @      rande ) * (_12 o. (0 2p1 & randu))) : (p. $:)

NB. ---------------------------------------------------------
NB. randtnf                                             _ _ _
NB. Truncated normal distribution TN(μ,σ²,a,b) of float
NB. numbers in range [a,b] with mean μ and variance σ²
NB.
NB. Syntax:
NB.   S=. [par] randtnf sh
NB. where
NB.   sh  - n-vector, shape of S
NB.   par - optional 4-vector (μ,σ,a,b), σ>0, a≤b, (μ,σ) are
NB.         normal distribution parameters, [a,b] is s range,
NB.         default is (0 1 -∞ ∞)
NB.   S   - sh-array of values s ~ N(μ,σ²), a ≤ s ≤ b
NB.   n   ≥ 0
NB.
NB. Application:
NB. - left-side truncated normal distribution TN(1,2²,3,∞) :
NB.   randltnf=: (1 2 3 _) & randtnf

randtnf=: (0 1 __ _ & $:) :(4 : 0)
  'mu sigma a b'=. x
  mrandnf=. (mu,sigma) & randnf
  NB. replace out-bounded elements recursively
  rober=. (I. @ ((a & >) +. (b & <))) ((($: @ mrandnf @ # @ [)`[`] }) ^: (0 < (#@[))) ]
  ($ (rober @ mrandnf @ (*/))) y
)

NB. ---------------------------------------------------------
NB. gemat                                               _ _ _
NB. Make random general matrix
NB.
NB. Formula:
NB.   g  ← mantissa * 2 ^ exponent
NB. where
NB.   mantissa ~ U(ma,mb)
NB.   exponent ~ TN(μ,σ,ea,eb)
NB.
NB. Syntax:
NB.   G=. [par] gemat sh
NB. where
NB.   sh  - n-vector, shape of G
NB.   par - optional 6-vector (ma,mb,μ,σ,ea,eb), σ>0, (ma,mb)
NB.         are mantissa's uniform distribution parameters,
NB.         and (μ,σ,ea,eb) are exponent's truncated normal
NB.         distribution parameters, default is:
NB.         (_1 1 0 (FP_FLEN/2) FP_EMIN FP_EMAX)
NB.   G   - random general sh-matrix
NB.   n   ≥ 0
NB.
NB. Notes:
NB. - default par provides about 95% of g numbers falls into
NB.   the range [-1/FP_EPS,-FP_EPS]U[FP_EPS,1/FP_EPS]
NB.   ("68-95-99.7 rule")

gemat=: ((_1 1 0 , (-: FP_FLEN) , FP_EMIN , FP_EMAX) & $:) :(((randu~ 2&{.) (* 2&^) (randtnf~ 2&}.))~)

NB. ---------------------------------------------------------
NB. hemat
NB. Adverb to make random Hermitian (symmetric) matrix
NB.
NB. Syntax:
NB.   H=. gemat hemat n
NB. where
NB.   n     ≥ 0, size of matrix H
NB.   gemat - monadic verb to generate random y-matrix (shape
NB.           is taken from y)
NB.   H     - random Hermitian (symmetric) (n,n)-matrix
NB.
NB. Application:
NB.   P=. (_1 1 2 3 _4 8 & gemat) hemat 4

hemat=: vu2y (+ ct)

NB. ---------------------------------------------------------
NB. unmat
NB. Adverb to make a random unitary (orthogonal) matrix with
NB. distribution given by Haar measure
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
NB. Syntax:
NB.   U=. randn unmat n
NB. where
NB.   n     ≥ 0, size of matrix U
NB.   randn - monadic verb to generate random y-matrix (shape
NB.           is taken from y) with elements distributed as
NB.           N(0,1), it is either randnf or randnc
NB.   U     - random unitary (orthogonal) (n,n)-matrix
NB.
NB. References:
NB. [1] Francesco Mezzadri (2007). How to generate random
NB.     matrices from the classical compact groups.
NB.     Notices of the AMS, Vol. 54 (2007), 592-604
NB.     http://arxiv.org/abs/math-ph/0609050v2

unmat=: ((* " 1) ((% |) @ ((<0 1) & |:))) & >/ @ qr vu2y

NB. ---------------------------------------------------------
NB. dimat
NB. Conjunction to make a random diagonalizable square matrix
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
NB.   D     - random diagonalizable square (n,n)-matrix
NB.
NB. Notes:
NB. - D will be Hermitian (symmetric) if randx produces float
NB.   (non-complex), possibly non-distinct numbers

dimat=: 2 : 'u (] mp (* ct)) v'

NB. ---------------------------------------------------------
NB. pomat
NB. Adverb to make a random Hermitian (symmetric) positive
NB. defined matrix
NB.
NB. Syntax:
NB.   P=. randx pomat n
NB. where
NB.   n     ≥ 0, size of matrix P
NB.   randx - monadic verb to generate random non-singular
NB.           (n,n)-matrix
NB.   P     - random Hermitian (symmetric) positive defined
NB.           (n,n)-matrix
NB.
NB. Application:
NB.   P=. (2p1 & rande) pomat 4

pomat=: vu2y (mp ct)
