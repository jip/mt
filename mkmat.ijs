NB. mkmat.ijs
NB. Make random matrices with certain properties
NB.
NB. mkun  make random unitary matrix
NB. mkor  make random orthogonal matrix
NB. mkpo  make random Hermitian (symmetric) positive defined
NB.       matrix
NB.
NB. XRef:
NB.   normalrand rand01

coclass 'mt'

NB. =========================================================
NB. Local constants

NB. =========================================================
NB. Local verbs

qr=: 128!:0

NB. =========================================================
NB. Interface verbs

NB. ---------------------------------------------------------
NB. mkun                                                    0
NB. mkor                                                    0
NB. Make a random unitary (orthogonal) matrix with
NB. distribution given by Haar measure
NB.
NB. Syntax:
NB.   Q=. mkun n
NB.   Q=. mkor n
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

mkun=: (3 : 0) " 0
  z=. (j./ normalrand 2 , (2 $ y)) % %: 2
  'q r'=. qr z
  d=. (<0 1) |: r
  ph=. (% |) d
  q (* " 1) ph  NB. Q * diag(d/|d|)
)

mkor=: (3 : 0) " 0
  z=. normalrand (2 $ y)
  'q r'=. qr z
  d=. (<0 1) |: r
  ph=. (% |) d
  q (* " 1) ph  NB. Q * diag(d/|d|)
)

NB. ---------------------------------------------------------
NB. mkpo                                                    0
NB. Adverb to make a random Hermitian (symmetric) positive
NB. defined matrix
NB.
NB. Formula:
NB.   Praw ← Q * diag(rand12(n)) * Q'
NB.   P ← (Praw + Praw') / 2
NB.
NB. Syntax:
NB.   P=. mkq mkpo n
NB. where
NB.   n   - size of Hermitian (symmetric) positive defined
NB.         matrix P
NB.   mkq - verb to make random unitary (orthogonal) matrix,
NB.         usually mkun or mkor
NB.   P   - n×n random Hermitian (symmetric) positive defined
NB.         matrix
NB.   n   ≥ 0

mkpo=: (>: @: rand01) ` (-: @ (+ ct) @ (] mp (* ct))) ` (`:6)
