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
NB. - system: 
NB. - mt: 

coclass 'mt'

NB. =========================================================
NB. Local definitions

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
NB.          distribution, default is open interval (0 1)
NB.   S    - sh-array of values s ~ U(a,b)
NB.   n    ≥ 0

randu=: (? @ $ 0:) : ((p.~ (-~/\))~ $:)

NB. ---------------------------------------------------------
NB. rande                                               _ _ _
NB. Exponential distribution E(μ) with mean μ
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
NB.   - (randu) values aren't used outside elsewhere

rande=: (- @ ^. @ randu) : (* $:)

NB. ---------------------------------------------------------
NB. randnf                                              _ _ _
NB. randnc                                              _ _ _
NB. Normal distribution N(μ,σ²) of float (complex) numbers
NB. with mean μ and variance σ²
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
NB. Formula (Box-Muller algorithm):
NB.   nf ← sqrt(-2*log(1-u1))*cos(2*π*u2)
NB.   nc ← (nf1 + i*nf2)/sqrt(2)
NB. where
NB.   u1,u2         ~ U(0,1)
NB.   nf,nf1,nf2,nc ~ N(0,1)

randnf=: ((%: @ (2 & rande)) * (  2 o. (0 2p1 & randu))) : (p. $:)
randnc=: ((%: @      rande ) * (_12 o. (0 2p1 & randu))) : (p. $:)

NB. ---------------------------------------------------------
NB. hemat
NB. Make random Hermitian (symmetric) matrix with elements
NB. distributed <<<<<<<<<<<<<#################
NB.
NB. Syntax:
NB.   he=. randx hemat n,n
NB. where
NB.   n - 
NB.   randx - monadic verb to generate elements

hemat=: (+ ct) @

NB. ---------------------------------------------------------
NB. runmat                                                  1
NB. rormat                                                  1
NB. Make a random unitary (orthogonal) matrix with
NB. distribution given by Haar measure
NB.
NB. Syntax:
NB.   Q=. runmat n,n
NB.   Q=. rormat n,n
NB. where
NB.   n - size of unitary (orthogonal) matrix Q
NB.   Q - n×n random unitary (orthogonal) matrix
NB.   n ≥ 0
NB.
NB. References:
NB. [1] Francesco Mezzadri (2007). How to generate random
NB.     matrices from the classical compact groups.
NB.     Notices of the AMS, Vol. 54 (2007), 592-604
NB.     http://arxiv.org/abs/math-ph/0609050v2

runmat=: (3 : 0) " 1
  z=. (j./ normalrand 2 , y) % %: 2
  'q r'=. qr z
  d=. (<0 1) |: r
  ph=. (% |) d
  q (* " 1) ph  NB. Q * diag(d/|d|)
)

rormat=: (3 : 0) " 1
  z=. normalrand y
  'q r'=. qr z
  d=. (<0 1) |: r
  ph=. (% |) d
  q (* " 1) ph  NB. Q * diag(d/|d|)
)

NB. ---------------------------------------------------------
NB. rpomat
NB. Adverb to make a random Hermitian (symmetric) positive
NB. defined matrix
NB.
NB. Formula:
NB.   Praw ← Q * diag(rand12(n)) * Q'
NB.   P ← (Praw + Praw') / 2
NB.
NB. Syntax:
NB.   P=. mkq rpomat n,n
NB. where
NB.   n   - size of Hermitian (symmetric) positive defined
NB.         matrix P
NB.   mkq - verb to make random unitary (orthogonal) matrix,
NB.         usually runmat or rormat
NB.   P   - n×n random Hermitian (symmetric) positive defined
NB.         matrix
NB.   n   ≥ 0

rpomat=: (>: @: rand01) ` (-: @ (+ ct) @ (] mp (* ct))) ` (`:6)

