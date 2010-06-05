NB. Random arrays
NB.
NB. gemat    Random general matrix
NB. trl1mat  Random unit lower triangular matrix
NB. trlmat   Random lower triangular matrix
NB. tru1mat  Random unit upper triangular matrix
NB. trumat   Random upper triangular matrix
NB. hemat    Random Hermitian (symmetric) matrix
NB. unmat    Random unitary (orthogonal) matrix
NB. dimat    Random diagonalizable matrix
NB. pomat    Random Hermitian (symmetric) positive definite
NB.          matrix
NB. ptmat    Random Hermitian (symmetric) positive definite
NB.          tridiagonal matrix
NB. spmat    Conj. to make verb to create sparse array
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Misc.

vu2y=: 2 : 'v @ u @ (2 & $)'  NB. v(u(y,y))

NB. ---------------------------------------------------------
NB. randu
NB.
NB. Description:
NB.   Uniform distribution U(a,b) with support (a,b)
NB.
NB. Syntax:
NB.   S=. [supp] randu sh
NB. where
NB.   sh   - r-array, shape of S
NB.   supp - optional 2-vector (a,b), a<b, support of
NB.          distribution, open interval, default is (0 1)
NB.   S    - sh-array of values s ~ U(a,b)
NB.   r    ≥ 0, the rank of S
NB.
NB. Formula:
NB.   s ← a + (b-a)*u
NB. where
NB.   u ~ U(0,1)
NB.   s ~ U(a,b)
NB.
NB. Application:
NB. - verb to generate complex matrix S with elements s
NB.   having Re(s) ~ U(0,1) and Im(s) ~ U(_1,1):
NB.     unirand=: randu j. (_1 1 & randu)

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
NB.   sh - r-array, shape of S
NB.   μ  > 0, optional mean, default is 1
NB.   S  - sh-array of values s ~ E(μ)
NB.   r  ≥ 0, the rank of S
NB.
NB. Formula:
NB.   s ← -μ*log(1-u)
NB. where
NB.   u ~ U(0,1)
NB.   s ~ E(μ)
NB.
NB. Application:
NB. - verb to generate complex matrix S with elements s
NB.   having Re(s) ~ E(1) and Im(s) ~ E(2):
NB.     exprand=: rande j. (2 & rande)
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
NB.   Normal distribution N(μ,σ^2) of real (complex) numbers
NB.   with mean μ and variance σ^2
NB.
NB. Syntax:
NB.   S=. [par] randnf sh
NB.   S=. [par] randnc sh
NB. where
NB.   sh  - r-array, shape of S
NB.   par - optional 2-vector (μ,σ), σ>0, parameters of
NB.         distribution, default is (0 1)
NB.   S   - sh-array of values s ~ N(μ,σ^2)
NB.   r   ≥ 0, the rank of S
NB.
NB. Formula:
NB.   nf1 ← sqrt(e)*cos(u)         NB. Box-Muller
NB.   nf2 ← sqrt(e)*sin(u)         NB.   algorithm
NB.   nc  ← (nf1 + i*nf2)/sqrt(2)  NB. D.R.Brillinger [1]
NB.   n3  ← σ*nxx + μ
NB. where
NB.   u            ~ U(0,2*π)
NB.   e            ~ E(2)
NB.   nf1,nf2,nc,n ~ N(0,1)
NB.   n3           ~ N(μ,σ^2)
NB.
NB. Algorithm for monadic randnf:
NB.   In: sh
NB.   Out: S
NB.   0) calc ceiled half of elements quantity in S:
NB.        n := ⌈Π{sh[i],i=0:r-1}/2⌉
NB.   1) generate random complex n-vector with elements
NB.      having real and imagine parts both independently
NB.      distributed normally:
NB.        Re(v[i]) := nf1[i], i=0:n-1
NB.        Im(v[i]) := nf2[i], i=0:n-1
NB.   2) form n×2-array with separated real and imagine
NB.      parts:
NB.        S[i,0] := Re(v[i])
NB.        S[i,1] := Im(v[i])
NB.   3) re-shape S to shape sh (special code is involved):
NB.        S=. sh ($,) S
NB.
NB. Application:
NB. - verb to generate real matrix S with elements s
NB.   having s ~ N(1,3^2):
NB.     normrand=: 1 3 & randnf
NB. - verb to generate complex matrix S with elements s
NB.   having s ~ N(1,3^2)
NB.     normrand=: 1 3 & randnc
NB. - verb to generate complex matrix S with elements s
NB.   having Re(s) ~ N(0,1) and Im(s) ~ N(1,3^2):
NB.     normrand=: randnf j. (1 3 & randnf)
NB.
NB. References:
NB. [1] D.R.Brillinger "Time series. Data analysis and
NB.     theory", The University of California, Berkeley, 1975
NB.     (Д. Бриллинджер "Временные ряды. Обработка данных и
NB.     теория". Изд-во "Мир". М. 1980, стр. 98)

randnf=: (($,) +.@:((%: @ (2 & rande)) r. (0 2p1 & randu))@>.@-:@(*/)) : (p. $:)
randnc=:           ((%: @      rande ) r. (0 2p1 & randu))             : (p. $:)

NB. ---------------------------------------------------------
NB. randtnf
NB.
NB. Description:
NB.   Truncated normal distribution TN(μ,σ^2,a,b) of real
NB.   numbers in range [a,b] with mean μ and variance σ^2
NB.
NB. Syntax:
NB.   S=. [par] randtnf sh
NB. where
NB.   sh   - r-array, shape of S
NB.   par - optional 4-vector (μ,σ,a,b), σ>0, a≤b, (μ,σ) are
NB.         normal distribution parameters, [a,b] is s range,
NB.         default is (0 1 -∞ +∞)
NB.   S   - sh-array of values s ~ N(μ,σ^2), a ≤ s ≤ b
NB.   r   ≥ 0, the rank of S
NB.
NB. Application:
NB. - verb to generate complex matrix S with elements s
NB.   having left-side truncated normal distribution
NB.   TN(1,3^2,4,+∞):
NB.     tnormrand=: 1 3 4 _ & randtnf

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
NB. gemat
NB.
NB. Description:
NB.   Make random general matrix
NB.
NB. Syntax:
NB.   G=. [par] gemat sh
NB. where
NB.   sh  - r-array, shape of G
NB.   par - optional 6-vector (ma,mb,μ,σ,ea,eb), σ>0, where
NB.         (ma,mb) are mantissa's uniform distribution
NB.         parameters, and (μ,σ,ea,eb) are exponent's
NB.         truncated normal distribution parameters, default
NB.         is: (_1 1 0 (FP_FLEN/2) FP_EMIN FP_EMAX)
NB.   G   - sh-matrix, random general
NB.   r   ≥ 0, the rank of G
NB.
NB. Formula:
NB.   g ← mantissa * 2 ^ exponent
NB. where
NB.   mantissa ~ U(ma,mb)
NB.   exponent ~ TN(μ,σ,ea,eb)
NB.
NB. Application:
NB. - generate real 4×4-matrix G with elements g having:
NB.     mantissa(g) ~ U(_1,1)
NB.     exponent(g) ~ TN(0,3^2,_6,4)
NB.   :
NB.     G=. _1 1 0 3 _6 4 gemat 4 4
NB. - generate complex 4×4-matrix G with elements g having:
NB.     mantissa(Re(g)),mantissa(Im(g)) ~ U(_1,1)
NB.     exponent(Re(g)),exponent(Im(g)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     G=. _1 1 0 3 _6 4 (gemat j. gemat) 4 4
NB. - generate complex 4×4-matrix G with elements g having:
NB.     mantissa(Re(g)) ~ U(0,1)
NB.     exponent(Re(g)) ~ TN(0,1,-∞,+∞)
NB.     mantissa(Im(g)) ~ U(_1,1)
NB.     exponent(Im(g)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     G=. ((0 1 0 1 __ _ & gemat) j. (_1 1 0 3 _6 4 & gemat)) 4 4
NB.
NB. Notes:
NB. - default par provides about 95% of g numbers falls into
NB.   the range [-1/FP_EPS,-FP_EPS]U[FP_EPS,1/FP_EPS]
NB.   ("68-95-99.7 rule")

gemat=: ((_1 1 0 , (-: FP_FLEN) , FP_EMIN , FP_EMAX) & $:) :(((randu~ 2&{.) (* 2&^) (randtnf~ 2&}.))~)

NB. ---------------------------------------------------------
NB. Verb:          Syntax:
NB. trl1mat        L1=. randx trl1mat sh
NB. trlmat         L=.  randx trlmat  sh
NB. tru1mat        U1=. randx tru1mat sh
NB. trumat         U=.  randx trumat  sh
NB.
NB. Description:
NB.   Adv. to make verb to generate random square triangular
NB.   matrix
NB. where
NB.   sh    - size or shape, either n or (n,n)
NB.   randx - monad to generate random y-matrix (shape is
NB.           taken from y)
NB.   L1    - n×n-matrix, random unit lower triangular
NB.   L     - n×n-matrix, random lower triangular
NB.   U1    - n×n-matrix, random unit upper triangular
NB.   U     - n×n-matrix, random upper triangular
NB.
NB. Algorithm for basic adverb trlmat:
NB.   In: n randx
NB.   Out: L
NB.   0) count non-zero elements in lower triangle:
NB.        ne := n*(n+1)/2
NB.   1) generate random ne-vector:
NB.        ve=. randx ne
NB.   2) prepare fret:
NB.      2.0) generate zero ne-vector:
NB.             f=. ne $ 0
NB.      2.1) generate lIOS for start of an interval marks:
NB.             lios=. +/\ i. ne
NB.      2.2) write in marks:
NB.             f=. 1 lios } f
NB.   3) cut ve on pieces of length 1 2 3 ... and stack them
NB.      into lower triangular matrix:
NB.        L=. f ];.1 ve
NB.
NB. Application:
NB. - generate real lower triangular 4×4-matrix L with
NB.   elements l having:
NB.     mantissa(l) ~ U(_1,1)
NB.     exponent(l) ~ TN(0,3^2,_6,4)
NB.   :
NB.     L=. (_1 1 0 3 _6 4 & gemat) trlmat 4 4
NB. - generate complex lower triangular 4×4-matrix L with
NB.   elements l having:
NB.     mantissa(Re(l)),mantissa(Im(l)) ~ U(_1,1)
NB.     exponent(Re(l)),exponent(Im(l)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     L=. (_1 1 0 3 _6 4 & (gemat j. gemat)) trlmat 4
NB. - generate complex lower triangular 4×4-matrix L with
NB.   elements l having:
NB.     mantissa(Re(l)) ~ U(0,1)
NB.     exponent(Re(l)) ~ TN(0,1,-∞,+∞)
NB.     mantissa(Im(l)) ~ U(_1,1)
NB.     exponent(Im(l)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     L=. ((0 1 0 1 __ _ & gemat) j. (_1 1 0 3 _6 4 & gemat)) trlmat 4
NB.
NB. Notes:
NB. - only n*(n+1)/2 numbers from RNG are requested
NB.
NB. TODO:
NB. - fret should be sparse

trlmat=:  1 : '1&([`(+/\@i.@{.@])`(0 $~ -:@(* >:)@{.@])}) ];.1 u@-:@(* >:)@{.'
trl1mat=: 1 : '(1;a:) & setdiag_mt_ @ (u trlmat_mt_)'
trumat=:  1 : '|: @ (u trlmat_mt_)'
tru1mat=: 1 : '|: @ (u trl1mat_mt_)'

NB. ---------------------------------------------------------
NB. hemat
NB.
NB. Description:
NB.   Adv. to make verb to make random Hermitian (symmetric)
NB.   matrix
NB.
NB. Syntax:
NB.   H=. randx hemat sh
NB. where
NB.   sh    - size or shape, either n or (n,n)
NB.   randx - monad to generate random y-matrix (shape is
NB.           taken from y)
NB.   H     - n×n-matrix, random Hermitian (symmetric)
NB.
NB. Assertion:
NB.   (-: ct) H
NB. where
NB.   H=. (gemat j. gemat) hemat n
NB.
NB. Application:
NB. - generate real symmetric 4×4-matrix H with
NB.   elements h having:
NB.     mantissa(h) ~ U(_1,1)
NB.     exponent(h) ~ TN(0,3^2,_6,4)
NB.   :
NB.     H=. (_1 1 0 3 _6 4 & gemat) hemat 4 4
NB. - generate complex Hermitian 4×4-matrix L with elements h
NB.   having:
NB.     mantissa(Re(h)),mantissa(Im(h)) ~ U(_1,1)
NB.     exponent(Re(h)),exponent(Im(h)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     H=. (_1 1 0 3 _6 4 & (gemat j. gemat)) hemat 4
NB. - generate complex Hermitian 4×4-matrix L with elements h
NB.   having:
NB.     mantissa(Re(h)) ~ U(0,1)
NB.     exponent(Re(h)) ~ TN(0,1,-∞,+∞)
NB.     mantissa(Im(h)) ~ U(_1,1)
NB.     exponent(Im(h)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     H=. ((0 1 0 1 __ _ & gemat) j. (_1 1 0 3 _6 4 & gemat)) hemat 4
NB.
NB. Notes:
NB. - only n*(n+1)/2 numbers from RNG are requested

hemat=: 1 : '(-: upddiag @ (u trlmat)) vu2y (+ ct_mt_)'

NB. ---------------------------------------------------------
NB. unmat
NB.
NB. Description:
NB.   Adv. to make verb to make a random unitary (orthogonal)
NB.   matrix with distribution given by Haar measure
NB.
NB. Syntax:
NB.   Q=. randnx unmat sh
NB. where
NB.   sh     - size or shape, either n or (n,n)
NB.   randnx - monad to generate random y-matrix (shape is
NB.            taken from y) with elements distributed as
NB.            N(0,1), it is either randnf or randnc
NB.   Q      - n×n-matrix, random unitary (orthogonal)
NB.
NB. Formula:
NB.   A ← randnx(n,n)
NB.   (Q,R) ← QR(A)
NB.   d ← diag(R)
NB.   Λ ← diagmat(d / |d|)
NB.   Q ← Q*Λ
NB. where
NB.   n - output matrix Q's size
NB.
NB. Assertion:
NB.   I -: clean (mp ct) Q
NB. where
NB.   I=. idmat n
NB.   Q=. randnf unmat n
NB.
NB. Application:
NB. - generate real orthogonal 4×4-matrix Q, where Q is
NB.   derived via QR-factorization from real matrix B having
NB.   elements b ~ N(1,3^2):
NB.     Q=. (1 3 & randnf) unmat 4 4
NB. - generate complex unitary 4×4-matrix U, where U is
NB.   derived via QR-factorization from complex matrix B
NB.   having elements b ~ N(1,3^2):
NB.     U=. (1 3 & randnc) unmat 4
NB.
NB. References:
NB. [1] Francesco Mezzadri (2007). How to generate random
NB.     matrices from the classical compact groups.
NB.     Notices of the AMS, Vol. 54 (2007), 592-604
NB.     http://arxiv.org/abs/math-ph/0609050v2
NB.
NB. Notes:
NB. - current gerqf implementation provides d≥0, so since R
NB.   isn't singular, hence d>0 and Λ is identity matrix

unmat=: vu2y (ungqr_mt_ @ geqrf_mt_)

NB. ---------------------------------------------------------
NB. dimat
NB.
NB. Description:
NB.   Conj. to make verb to make a random diagonalizable
NB.   square matrix
NB.
NB. Syntax:
NB.   A=. randx dimat randq sh
NB. where
NB.   sh    - size or shape, either n or (n,n)
NB.   randq - monad to generate random unitary (orthogonal)
NB.           (y,y)-matrix (size is taken from y)
NB.   randx - monad to generate random y-vector (length is
NB.           taken from y) of distinct or real (non-complex)
NB.           numbers
NB.   A     - n×n-matrix, random diagonalizable square
NB.
NB. Formula:
NB.   A ← U * D * U^H
NB. where
NB.   U=. randq sh
NB.   D=. diadmat d
NB.   d=. randx n
NB.
NB. Application:
NB. - generate real diagonalizable 4×4-matrix A with
NB.   eigenvalues d having:
NB.     mantissa(d) ~ U(1,3)
NB.     exponent(d) ~ TN(0,4^2,_5,6)
NB.   and eigenvectors Q derived via QR-factorization from
NB.   real matrix B having elements b ~ N(1,3^2):
NB.     A=. (1 3 0 4 _5 6 & gemat) dimat ((1 3 & randnf) unmat) 4 4
NB. - generate complex diagonalizable Hermitian 4×4-matrix A
NB.   with eigenvalues d having:
NB.     mantissa(d) ~ U(1,3)
NB.     exponent(d) ~ TN(0,4^2,_5,6)
NB.   and eigenvectors U derived via QR-factorization from
NB.   complex matrix B having elements b ~ N(1,3^2):
NB.     A=. (1 3 0 4 _5 6 & gemat) dimat ((1 3 & randnc) unmat) 4
NB. - generate complex diagonalizable 4×4-matrix A with
NB.   eigenvalues d having:
NB.     mantissa(Re(d)),mantissa(Im(d)) ~ U(1,3)
NB.     exponent(Re(d)),exponent(Im(d)) ~ TN(0,4^2,_5,6)
NB.   and eigenvectors U derived via QR-factorization from
NB.   complex matrix B having elements b ~ N(1,3^2):
NB.     A=. (1 3 0 4 _5 6 & (gemat j. gemat)) dimat ((1 3 & randnc) unmat) 4
NB.
NB. Notes:
NB. - A will be Hermitian (symmetric) if randx produces real
NB.   (non-complex), possibly non-distinct numbers
NB. - A is diagonalizable iif A is normal (A^H * A = A * A^H)
NB.
NB. TODO:
NB. - assertion

dimat=: 2 : 'u@{. (] mp (* ct_mt_)) v'

NB. ---------------------------------------------------------
NB. pomat
NB.
NB. Description:
NB.   Adv. to make verb to make a random Hermitian
NB.   (symmetric) positive definite matrix
NB.
NB. Syntax:
NB.   P=. randx pomat sh
NB. where
NB.   sh    - size or shape, either n or (n,n)
NB.   randx - monad to generate random non-singular
NB.           n×n-matrix
NB.   P     - n×n-matrix, random Hermitian (symmetric)
NB.           positive definite
NB.
NB. Application:
NB. - generate real symmetric positive definite 4×4-matrix P
NB.   with elements p having:
NB.     mantissa(p) ~ U(1,3)
NB.     exponent(p) ~ TN(0,4^2,_5,6)
NB.   :
NB.     P=. (1 3 0 4 _5 6 & gemat) pomat 4 4
NB. - generate complex Hermitian positive definite 4×4-matrix
NB.   P with elements p having:
NB.     mantissa(Re(p)),mantissa(Im(p)) ~ U(1,3)
NB.     exponent(Re(p)),exponent(Im(p)) ~ TN(0,4^2,_5,6)
NB.   :
NB.     P=. (1 3 0 4 _5 6 & (gemat j. gemat)) pomat 4
NB. - generate complex Hermitian positive definite 4×4-matrix
NB.   P with elements p having:
NB.     mantissa(Re(p)) ~ U(1,2)
NB.     exponent(Re(p)) ~ TN(0,1,-∞,+∞)
NB.     mantissa(Im(p)) ~ U(1,3)
NB.     exponent(Im(p)) ~ TN(0,4^2,_5,6)
NB.   :
NB.     P=. ((1 2 0 1 __ _ & gemat) j. (1 3 0 4 _5 6 & gemat)) pomat 4
NB. - generate real symmetric negative definite 4×4-matrix P
NB.   with elements p having:
NB.     mantissa(p) ~ U(1,3)
NB.     exponent(p) ~ TN(0,4^2,_5,6)
NB.   :
NB.     P=. -@((1 3 0 4 _5 6 & gemat) pomat) 4 4
NB.
NB. TODO:
NB. - assertion

pomat=: vu2y (mp ct_mt_)

NB. ---------------------------------------------------------
NB. ptmat
NB.
NB. Description:
NB.   Adv. to make verb to make a random Hermitian
NB.   (symmetric) positive definite tridiagonal matrix
NB.
NB. Syntax:
NB.   T=. randx ptmat sh
NB. where
NB.   sh    - size or shape, either n or (n,n)
NB.   randx - monad to generate random vectors d and e; is
NB.           evoked as:
NB.             d=. (| @ (9 o. randx)) n
NB.             e=. randx (n-1)
NB.   d     - n-vector of positive numbers, the main diagonal
NB.           of D
NB.   e     - (n-1)-vector, the subdiagonal of unit lower
NB.           bidiagonal matrix L1
NB.   T     - n×n-matrix, random Hermitian (symmetric)
NB.           positive definite tridiagonal, defined as:
NB.             T = L1 * D * L1^H
NB.
NB. Application:
NB. - generate real symmetric positive definite tridiagonal
NB.   4×4-matrix T with elements d and e having:
NB.     mantissa(d),mantissa(e) ~ U(1,3)
NB.     mantissa(d),exponent(e) ~ TN(0,4^2,_5,6)
NB.   :
NB.     T=. (1 3 0 4 _5 6 & gemat) ptmat 4 4
NB. - generate complex Hermitian positive definite
NB.   tridiagonal 4×4-matrix T with elements d and e having:
NB.     mantissa(d),mantissa(Re(e)),mantissa(Im(e)) ~ U(1,3)
NB.     exponent(d),exponent(Re(e)),exponent(Im(e)) ~ TN(0,4^2,_5,6)
NB.   :
NB.     T=. (1 3 0 4 _5 6 & (gemat j. gemat)) ptmat 4
NB. - generate complex Hermitian positive definite
NB.   tridiagonal 4×4-matrix T with elements d and e having:
NB.     mantissa(d),mantissa(Re(e)) ~ U(0,1)
NB.     mantissa(d),exponent(Re(e)) ~ TN(0,1,-∞,+∞)
NB.     mantissa(Im(e)) ~ U(1,3)
NB.     exponent(Im(e)) ~ TN(0,4^2,_5,6)
NB.   :
NB.     T=. (gemat j. (1 3 0 4 _5 6 & gemat)) ptmat 4
NB.
NB. TODO:
NB. - T should be sparse

ptmat=: 1 : '(((u@<: ; _1:) setdiag_mt_ idmat_mt_) (mp mp ct_mt_@[) (diagmat_mt_@:|@(9 o. u)))@{.'

NB. ---------------------------------------------------------
NB. spmat
NB.
NB. Description:
NB.   Conj. to make verb to create sparse array
NB.
NB. Syntax:
NB.   vapp=. randx spmat ratio
NB. where
NB.   randx - monad to generate random numbers; is evoked as:
NB.             A=. randx sh
NB.   ratio - scalar in range [0,1], specifies desired
NB.           portion of non-zero elements
NB.   vapp  - monad to create sparse array; is evoked as:
NB.             S=. vapp sh
NB.   sh    - r-vector of non-negative integers, a shape of
NB.           arrays A and S
NB.   r     ≥ 0, rank of A and S
NB.
NB. TODO:
NB. - S should be sparse

spmat=: 2 : '(u@<.@(n*(*/))) ((?@$~ #)~ (*/@$)) } ($&0)'

