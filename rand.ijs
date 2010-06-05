NB. rand.ijs
NB. Random objects
NB.
NB. randu   uniform distribution U(a,b)
NB. rande   exponential distribution E(μ)
NB. randnf  normal distribution N(μ,σ²) of float numbers
NB. randnc  normal distribution N(μ,σ²) of complex numbers
NB.
NB. hemat   random Hermitian (symmetric) matrix
NB. unmat   random unitary (orthogonal) matrix
NB. dimat   random diagonalizable matrix
NB. pomat   random Hermitian (symmetric) positive defined matrix
NB.
NB. XRef:
NB. - system: mp
NB. - mt: ct qr

coclass 'mt'

NB. =========================================================
NB. Local definitions

qr=: 128!:0           NB. built-in QR factorization (TODO: implement)

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
NB.   nf1 ← sqrt(-2*log(1-u1))*cos(2*π*u2)  NB. Box-Muller
NB.   nf2 ← sqrt(-2*log(1-u1))*sin(2*π*u2)  NB.   algorithm
NB.   nc  ← (nf1 + i*nf2)/sqrt(2)
NB. where
NB.   u1,u2      ~ U(0,1)
NB.   nf1,nf2,nc ~ N(0,1)
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

randnf=: ((%: @ (2 & rande)) * (  2 o. (0 2p1 & randu))) : (p. $:)
randnc=: ((%: @      rande ) * (_12 o. (0 2p1 & randu))) : (p. $:)

NB. ---------------------------------------------------------
NB. hemat
NB. Adverb to make random Hermitian (symmetric) matrix with
NB. elements of some probability distribution
NB.
NB. Syntax:
NB.   H=. randx hemat n
NB. where
NB.   n     ≥ 0, size of matrix H
NB.   randx - monadic verb to generate random y-matrix (shape
NB.           is taken from y)
NB.   H     - random Hermitian (symmetric) (n,n)-matrix
NB.
NB. Notes:
NB. - supplying verbs (randu rande randnf randnc) won't give
NB.   values distributed exactly so, but closer

hemat=: 2 : '(+ ct) @ u @ (2 & $)'

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
NB.   U=. randx unmat n
NB. where
NB.   n     ≥ 0, size of matrix U
NB.   randx - monadic verb to generate random y-matrix (shape
NB.           is taken from y)
NB.   U     - random unitary (orthogonal) (n,n)-matrix
NB.
NB. References:
NB. [1] Francesco Mezzadri (2007). How to generate random
NB.     matrices from the classical compact groups.
NB.     Notices of the AMS, Vol. 54 (2007), 592-604
NB.     http://arxiv.org/abs/math-ph/0609050v2

unmat=: 1 : '((* " 1) ((% |) @ ((<0 1) & |:))) & >/ @ qr @ u @ (2 & $)'

NB. ---------------------------------------------------------
NB. dimat
NB. Conjunction to make a random diagonalizable non-Hermitian
NB. (non-symmetric) square matrix
NB.
NB. Syntax:
NB.   D=. randx dimat randq n
NB. where
NB.   n     ≥ 0, size of matrix D
NB.   randq - monadic verb to generate random unitary
NB.           (orthogonal) y-matrix (shape is taken from y)
NB.   randx - monadic verb to generate random y-vector
NB.           (length is taken from y)
NB.   D     - random diagonalizable non-Hermitian (non-
NB.           symmetric) square (n,n)-matrix
NB.
NB. tau/lti:
NB.   NB. reconstruct inv(V) from U:
NB.   NB. α=. diag(U^H * V)
NB.   NB. inv(V)=. diagmat(1/α) * U^H
NB.   invV=. U (((|: @ ]) % (+/ @: *)) +)~ V

dimat=: 2 : ''

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

pomat=: 1 : '(mp ct) @ u @ (2 & $)'
