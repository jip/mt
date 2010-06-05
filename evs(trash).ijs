NB. ---------------------------------------------------------
NB. heevs
NB. Solve system A*X=B, where A is a Hermitian matrix
NB.
NB. Syntax:
NB.   x=. b heevs (rv ; ev)
NB. where
NB.   rv - N×N-matrix, right eigenvectors of A
NB.   ev - N-vector, eigenvalues of A
NB.   x  - N-vector, solution
NB.   N >= 0

heevs=: [: : ((0 {:: ]) % (1 {:: ])) mp (mp~ (+ @ |: @ (0 & {::)))

NB. ---------------------------------------------------------
NB. dievs
NB. Solve system A*x=b, where A is a diagonalizable matrix
NB.
NB. Syntax:
NB.   x=. b dievs (rv ; ev ; rvi)
NB. where
NB.   rv  - N×N-matrix, right eigenvectors of A
NB.   ev  - N-vector, eigenvalues of A
NB.   rvi - N×N-matrix, inversion of rv
NB.   x   - N-vector, solution
NB.   N  >= 0

dievs=: [: : ((0 {:: ]) % (1 {:: ])) mp (mp~ (2 & {::))

NB. ---------------------------------------------------------
NB. poevs
NB. Solve system A*x=b, where A is a Hermitian (symmetric)
NB. positive definite matrix
NB.
NB. Syntax:
NB.   X=. B poevs A
NB. where
NB.   A    - n×n-matrix, right eigenvectors of A
NB.   B    - n-vector or n×nrhs-matrix, RHS
NB.   X    - same shape as B, solution
NB.   n    ≥ 0
NB.   nrhs ≥ 0

poevs=: heevs

