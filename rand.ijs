NB. Random arrays
NB.
NB. trxxmat    Adv. to make verb to make random square
NB.            triangular matrix
NB. kmsmat     Adv. to make monad to make random
NB.            Kac-Murdock-Szego (KMS) matrix
NB. gemat      Make random real array
NB. dimat      Conj. to make verb to make random
NB.            diagonalizable matrix
NB. hemat      Adv. to make verb to make random Hermitian
NB.            (symmetric) matrix
NB. pomat      Adv. to make verb to make random Hermitian
NB.            (symmetric) positive definite matrix
NB. ptmatx     Adv. to make verb to make random Hermitian
NB.            (symmetric) positive definite tridiagonal
NB.            matrix
NB. unmat      Adv. to make verb to make random unitary
NB.            (orthogonal) matrix
NB. spmat      Conj. to make verb to make random sparse array
NB.
NB. testtrmat  Test trxxmat by matrix size given
NB. testgemat  Test gemat by matrix shape given
NB. testdimat  Test dimat by matrix size given
NB. testhemat  Test hemat by matrix size given
NB. testpomat  Test pomat by matrix size given
NB. testptmat  Test ptmatx by matrix size given
NB. testspmat  Test spmat by matrix shape given
NB. testrand   Test xxxxmatx by matrix shape given
NB.
NB. Version: 0.10.5 2020-03-30
NB.
NB. Copyright 2010-2020 Igor Zhuravlov
NB.
NB. This file is part of mt
NB.
NB. mt is free software: you can redistribute it and/or
NB. modify it under the terms of the GNU Lesser General
NB. Public License as published by the Free Software
NB. Foundation, either version 3 of the License, or (at your
NB. option) any later version.
NB.
NB. mt is distributed in the hope that it will be useful, but
NB. WITHOUT ANY WARRANTY; without even the implied warranty
NB. of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
NB. See the GNU Lesser General Public License for more
NB. details.
NB.
NB. You should have received a copy of the GNU Lesser General
NB. Public License along with mt. If not, see
NB. <http://www.gnu.org/licenses/>.

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. randu
NB.
NB. Description:
NB.   Uniform distribution U(a,b) with support (a,b)
NB.
NB. Syntax:
NB.   S=. [supp] randu sh
NB. where
NB.   sh   - r-vector of non-negative integers, shape of S
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
NB. - verb to make complex matrix S with elements s having:
NB.     Re(s) ~ U(0,1)
NB.     Im(s) ~ U(_1,1)
NB.   :
NB.     unirand=: randu j. _1 1&randu

randu=: (?@$&0) :((p.~ -~/\)~ $:)

NB. ---------------------------------------------------------
NB. rande
NB.
NB. Description:
NB.   Exponential distribution E(μ) with mean μ
NB.
NB. Syntax:
NB.   S=. [μ] rande sh
NB. where
NB.   sh - r-vector of non-negative integers, shape of S
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
NB. - verb to make complex matrix S with elements s having:
NB.     Re(s) ~ E(1)
NB.     Im(s) ~ E(2)
NB.   :
NB.     exprand=: rande j. 2&rande
NB.
NB. Notes:
NB. - randu is used instead of (1-randu) since they are
NB.   equivalent under following assumptions:
NB.   - randu's support is open interval (0,1)
NB.   - randu's values aren't used outside elsewhere

rande=: -@^.@randu : (* $:)

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
NB.   sh  - r-vector of non-negative integers, shape of S
NB.   par - optional 2-vector (μ,σ), σ>0, parameters of
NB.         distribution, default is (0 1)
NB.   S   - sh-array of values s ~ N(μ,σ^2)
NB.   r   ≥ 0, the rank of S
NB.
NB. Formula:
NB.   nf1 ← sqrt(e)*cos(u)         NB. see...
NB.   nf2 ← sqrt(e)*sin(u)         NB. ...[1]
NB.   nc  ← (nf1 + i*nf2)/sqrt(2)  NB. see [2]
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
NB.   1) make random complex n-vector with elements
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
NB. - verb to make real matrix S with elements s having:
NB.     s ~ N(1,3^2)
NB.   :
NB.     normrand=: 1 3&randnf
NB. - verb to make complex matrix S with elements s having:
NB.     s ~ N(1,3^2)
NB.   :
NB.     normrand=: 1 3&randnc
NB. - verb to make complex matrix S with elements s having:
NB.     Re(s) ~ N(0,1)
NB.     Im(s) ~ N(1,3^2)
NB.   :
NB.     normrand=: randnf j. 1 3&randnf
NB.
NB. Notes:
NB. - randnf requests only ⌈Π{sh[:]}/2⌉ numbers from RNG
NB.
NB. References:
NB. [1] G. E. P. Box, Mervin E. Muller. A Note on the
NB.     Generation of Random Normal Deviates. The Annals of
NB.     Mathematical Statistics, 1958, Vol. 29, No. 2, pp.
NB.     610-611.
NB. [2] D. R. Brillinger. Time series. Data analysis and
NB.     theory. The University of California, Berkeley, 1975
NB.     (Д. Бриллинджер. Временные ряды. Обработка данных и
NB.     теория. Изд-во "Мир". М. 1980, стр. 98).

randnf=: (($,) +.@:(%:@(2&rande) r. 0 2p1&randu)@>.@-:@(*/)) : (p. $:)
randnc=:           (%:@   rande  r. 0 2p1&randu)             : (p. $:)

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
NB.   sh   - r-vector of non-negative integers, shape of S
NB.   par - optional 4-vector (μ,σ,a,b), σ>0, a≤b, (μ,σ) are
NB.         normal distribution parameters, [a,b] is s range,
NB.         default is (0 1 -∞ +∞)
NB.   S   - sh-array of values s ~ N(μ,σ^2), a ≤ s ≤ b
NB.   r   ≥ 0, the rank of S
NB.
NB. Application:
NB. - verb to make complex matrix S with elements s having
NB.   left-side truncated normal distribution:
NB.     s ~ TN(1,3^2,4,+∞)
NB.   :
NB.     tnormrand=: 1 3 4 _&randtnf

randtnf=: 0 1 __ _&$: :(4 : 0)
  'mu sigma a b'=. x
  mrandnf=. (mu,sigma)&randnf
  NB. replace out-bounded elements recursively
  rober=. I.@(a&> +. b&<) $:@mrandnf@#@[`[`]}^:(0 < #@[) ]
  ($ rober@mrandnf@(*/)) y
)

NB. ---------------------------------------------------------
NB. kmsmatrho
NB.
NB. Description:
NB.   Adv. to make nilad to produce rho parameter for
NB.   rho-parametrized KMS matrix K(rho)
NB.
NB. Syntax:
NB.   vapp=. randx kmsmatrho
NB. where
NB.   randx - monad to generate random numbers; is called as:
NB.             A=. randx sh
NB.   vapp  - nilad to generate rho; is called as:
NB.             rho=. vapp any_noun
NB.   sh    - size or shape, the shape of A
NB.
NB. Notes:
NB. - randx is used to acquire datatype only
NB. - datatypes supported: floating, complex

kmsmat121=: 1 : '{.@u@1:'
kmsmatrho=: 1 : '((0.5 1&randu_mt_ r. 0 2p1&randu_mt_) kmsmat121_mt_)`(0.5 1&randu_mt_ kmsmat121_mt_)@.(0 -: 11&o.)@(u kmsmat121_mt_)'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:      Syntax:                 Is called as:
NB. trl1mat    vapp=. randx trl1mat    L1=. vapp sh
NB. trlmat     vapp=. randx trlmat     L=.  vapp sh
NB. tru1mat    vapp=. randx tru1mat    U1=. vapp sh
NB. trumat     vapp=. randx trumat     U=.  vapp sh
NB.
NB. Description:
NB.   Adv. to make verb to make random square triangular
NB.   matrix
NB. where
NB.   randx - monad to make random y-array; is called as:
NB.             A=. randx y
NB.   vapp  - monad to make triangular n×n-matrix T; is
NB.           called as:
NB.             T=. vapp sh
NB.   sh    - size or shape, is either n or (n,any_number)
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
NB.   1) make random ne-vector:
NB.        ve=. randx ne
NB.   2) prepare fret:
NB.      2.0) make zero ne-vector:
NB.             f=. ne $ 0
NB.      2.1) make lISO for start of an interval marks:
NB.             liso=. +/\ i. ne
NB.      2.2) write in marks:
NB.             f=. 1 liso} f
NB.   3) cut ve on pieces of length (1 2 3 ... n) and stack
NB.      them into lower triangular matrix:
NB.        L=. f ];.1 ve
NB.
NB. Application:
NB. - make real lower triangular 4×4-matrix L with elements l
NB.   having:
NB.     mantissa(l) ~ U(_1,1)
NB.     exponent(l) ~ TN(0,3^2,_6,4)
NB.   :
NB.     L=. _1 1 0 3 _6 4&gemat trlmat 4 4
NB. - make complex lower triangular 4×4-matrix L with
NB.   elements l having:
NB.     mantissa(Re(l)),mantissa(Im(l)) ~ U(_1,1)
NB.     exponent(Re(l)),exponent(Im(l)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     L=. _1 1 0 3 _6 4&(gemat j. gemat) trlmat 4
NB. - make complex lower triangular 4×4-matrix L with
NB.   elements l having:
NB.     mantissa(Re(l)) ~ U(0,1)
NB.     exponent(Re(l)) ~ TN(0,1,-∞,+∞)
NB.     mantissa(Im(l)) ~ U(_1,1)
NB.     exponent(Im(l)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     L=. (0 1 0 1 __ _&gemat j. _1 1 0 3 _6 4&gemat) trlmat 4
NB.
NB. Notes:
NB. - only n*(n+1)/2 numbers from RNG are requested
NB.
NB. TODO:
NB. - fret should be sparse

trlmat=:  1 : '1&([`(+/\@i.@{.@])`(0 $~ -:@(* >:)@{.@])}) ];.1 u@-:@(* >:)@{.'
trl1mat=: 1 : '(1;a:)&setdiag_mt_@(u trlmat_mt_)'
trumat=:  1 : '|:@(u trlmat_mt_)'
tru1mat=: 1 : '|:@(u trl1mat_mt_)'

NB. ---------------------------------------------------------
NB. kmsmat
NB.
NB. Description:
NB.   Adv. to make monad to make random Kac-Murdock-Szego
NB.   (KMS) matrix [1].
NB.
NB. Syntax:
NB.   vapp=. randx kmsmat
NB. where
NB.   randx - nilad to make rho; is called as:
NB.             rho=. randx any_noun
NB.   vapp  - monad to make KMS; is called as:
NB.             K=. vapp sh
NB.   sh    - size or shape, is either n or (n,any_number)
NB.   rho   - scalar number, the KMS matrix parameter
NB.   K     - n×n-matrix, random KMS matrix
NB.
NB. Notes:
NB. - KMS matrix is defined as:
NB.     K ≡ K(rho,i,j) := {      rho ^(j-i), i ≤ j
NB.                       { conj(rho)^(i-j), i > j
NB. - KMS matrix is Toeplitz
NB.   KMS matrix is symmetric for real rho and is Hermitian
NB.   for complex rho
NB. - a factorization:
NB.     L1 * D * L1^H = K
NB.   exists for KMS matrix, where
NB.     L1 - unit lower triangular matrix
NB.     D  - diagonal matrix
NB. - KMS matrix is positive definite iif:
NB.     0 < abs(rho) < 1
NB.   and is well conditioned for large n (>400) if:
NB.     0.5 < abs(rho) < 1
NB. - K^_1 is tridiagonal
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   K -: D (] mp (mp ct)) L1  NB. verify factorization: L1 * D * L1^H = K
NB.   K -: (mp ct) L            NB. K is positive definite, i.e. Cholesky factorization exists
NB.   iK -: gtpick iK           NB. K^_1 is tridiagonal
NB. where
NB.   n=. 7                     NB. matrix size, any
NB.   rho=. 0.5j0.5             NB. matrix parameter, any, adhere (0 < abs(rho) < 1) to force K to be positive definite
NB.   mkrho=. rho"_             NB. we must know rho value a priori to define L1 and D
NB.   K=.  mkrho kmsmat n
NB.   L1=. trl1pick %. ((+ - rho) ; _1) setdiag idmat n
NB.   D=. diagmat (1) 0} n $ 1 - *: | rho
NB.   L=. potrfl K              NB. try Cholesky factorization: L * L^H = K
NB.   iK=. %. K
NB.
NB. Application:
NB. - make random real positive definite KMS 4×4-matrix:
NB.     K=. {.@randu@1: kmsmat 4 4
NB. - make random complex positive definite KMS 4×4-matrix:
NB.     K=. {.@(randu r. 0 2p1&randu)@1: kmsmat 4
NB.
NB. References:
NB. [1] W.F. Trench, "Numerical solution of the eigenvalue
NB.     problem for Hermitian Toeplitz matrices", SIAM J.
NB.     Matrix Analysis and Appl., 10 (1989), pp. 135-146.

kmsmat=: 1 : 'he_mt_@(+@u ^ -/~@i.@{.)'

NB. ---------------------------------------------------------
NB. gemat
NB.
NB. Description:
NB.   Make random real array
NB.
NB. Syntax:
NB.   G=. [par] gemat sh
NB. where
NB.   sh  - r-vector of non-negative integers, shape of G
NB.   par - optional 6-vector (ma,mb,μ,σ,ea,eb), σ>0, where
NB.         (ma,mb) are mantissa's uniform distribution
NB.         parameters, and (μ,σ,ea,eb) are exponent's
NB.         truncated normal distribution parameters, default
NB.         is: (_1 1 0 (FP_FLEN/2) FP_EMIN FP_EMAX)
NB.   G   - sh-array, random
NB.   r   ≥ 0, the rank of G
NB.
NB. Formula:
NB.   g ← mantissa * 2 ^ exponent
NB. where
NB.   mantissa ~ U(ma,mb)
NB.   exponent ~ TN(μ,σ,ea,eb)
NB.
NB. Application:
NB. - make real 4×4-matrix G with elements g having:
NB.     mantissa(g) ~ U(_1,1)
NB.     exponent(g) ~ TN(0,3^2,_6,4)
NB.   :
NB.     G=. _1 1 0 3 _6 4 gemat 4 4
NB. - make complex 4×4-matrix G with elements g having:
NB.     mantissa(Re(g)),mantissa(Im(g)) ~ U(_1,1)
NB.     exponent(Re(g)),exponent(Im(g)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     G=. _1 1 0 3 _6 4 (gemat j. gemat) 4 4
NB. - make complex 4×4-matrix G with elements g having:
NB.     mantissa(Re(g)) ~ U(0,1)
NB.     exponent(Re(g)) ~ TN(0,1,-∞,+∞)
NB.     mantissa(Im(g)) ~ U(_1,1)
NB.     exponent(Im(g)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     G=. (0 1 0 1 __ _&gemat j. _1 1 0 3 _6 4&gemat) 4 4
NB.
NB. Notes:
NB. - default par provides about 95% of g numbers falls into
NB.   the range [-1/FP_EPS,-FP_EPS]U[FP_EPS,1/FP_EPS]
NB.   ("68-95-99.7 rule")

gemat=: (_1 1 0 , (-: FP_FLEN) , FP_EMIN , FP_EMAX)&$: :(((randu~ 2&{.) (* 2&^) (randtnf~ 2&}.))~)

NB. ---------------------------------------------------------
NB. dimat
NB.
NB. Description:
NB.   Conj. to make verb to make random diagonalizable
NB.   matrix
NB.
NB. Syntax:
NB.   vapp=. randx dimat randq
NB. where
NB.   randq - monad to make Q; is called as:
NB.             Q=. randq sh
NB.   randx - monad to make d; is called as:
NB.             d=. randx n
NB.   vapp  - monad to make A; is called as:
NB.             A=. vapp sh
NB.   sh    - size or shape, is either n or (n,n)
NB.   Q     - n×n-matrix, random unitary (orthogonal):
NB.             I = Q^H * Q
NB.   d     - n-vector of distinct or real (non-complex)
NB.           numbers
NB.   A     - n×n-matrix, random diagonalizable square, is
NB.           defined as:
NB.             A := Q * D * Q^H
NB.   D     - n×n-matrix, diagonal, is defined as:
NB.             D := diagmat d
NB.
NB. Application:
NB. - make real diagonalizable 4×4-matrix A with eigenvalues
NB.   d having:
NB.     mantissa(d) ~ U(1,3)
NB.     exponent(d) ~ TN(0,4^2,_5,6)
NB.   and eigenvectors Q derived via QR-factorization from
NB.   real matrix B:
NB.     A=. 1 3 0 4 _5 6&gemat dimat (randnf unmat) 4 4
NB. - make complex diagonalizable Hermitian 4×4-matrix A with
NB.   eigenvalues d having:
NB.     mantissa(d) ~ U(1,3)
NB.     exponent(d) ~ TN(0,4^2,_5,6)
NB.   and eigenvectors Q derived via QR-factorization from
NB.   complex matrix B:
NB.     A=. 1 3 0 4 _5 6&gemat dimat (randnc unmat) 4
NB. - make complex diagonalizable 4×4-matrix A with
NB.   eigenvalues d having:
NB.     mantissa(Re(d)),mantissa(Im(d)) ~ U(1,3)
NB.     exponent(Re(d)),exponent(Im(d)) ~ TN(0,4^2,_5,6)
NB.   and eigenvectors Q derived via QR-factorization from
NB.   complex matrix B:
NB.     A=. 1 3 0 4 _5 6&(gemat j. gemat) dimat (randnc unmat) 4
NB.
NB. Notes:
NB. - A will be Hermitian (symmetric) if randx produces real
NB.   (non-complex), possibly non-distinct numbers
NB. - A is diagonalizable iif A is normal (A^H * A = A * A^H)

dimat=: 2 : 'u@{. (] mp_mt_ (* ct_mt_)) v'

NB. ---------------------------------------------------------
NB. hemat
NB.
NB. Description:
NB.   Adv. to make verb to make random Hermitian (symmetric)
NB.   matrix
NB.
NB. Syntax:
NB.   vapp=. randx hemat
NB. where
NB.   randx - monad to make random y-array; is called as:
NB.             A=. randx y
NB.   vapp  - monad to make H; is called as:
NB.             H=. vapp sh
NB.   sh    - size or shape, is either n or (n,any_number)
NB.   H     - n×n-matrix, random Hermitian (symmetric)
NB.
NB. Assertions:
NB.   (-: ct) H
NB. where
NB.   H=. randx hemat n
NB.
NB. Application:
NB. - make real symmetric 4×4-matrix H with  elements h
NB.  having:
NB.     mantissa(h) ~ U(_1,1)
NB.     exponent(h) ~ TN(0,3^2,_6,4)
NB.   :
NB.     H=. _1 1 0 3 _6 4&gemat hemat 4 4
NB. - make complex Hermitian 4×4-matrix L with elements h
NB.   having:
NB.     mantissa(Re(h)),mantissa(Im(h)) ~ U(_1,1)
NB.     exponent(Re(h)),exponent(Im(h)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     H=. _1 1 0 3 _6 4&(gemat j. gemat) hemat 4
NB. - make complex Hermitian 4×4-matrix L with elements h
NB.   having:
NB.     mantissa(Re(h)) ~ U(0,1)
NB.     exponent(Re(h)) ~ TN(0,1,-∞,+∞)
NB.     mantissa(Im(h)) ~ U(_1,1)
NB.     exponent(Im(h)) ~ TN(0,3^2,_6,4)
NB.   :
NB.     H=. (0 1 0 1 __ _&gemat j. _1 1 0 3 _6 4&gemat) hemat 4
NB.
NB. Notes:
NB. - only n*(n+1)/2 numbers from RNG are requested

hemat=: 1 : 'tr2he_mt_@(u trlmat_mt_)@(2&$)'

NB. ---------------------------------------------------------
NB. pomat
NB.
NB. Description:
NB.   Adv. to make verb to make random Hermitian (symmetric)
NB.   positive definite matrix
NB.
NB. Syntax:
NB.   vapp=. randx pomat
NB. where
NB.   randx - monad to make A; is called as:
NB.             A=. randx y
NB.   vapp  - monad to make P; is called as:
NB.             P=. vapp sh
NB.   sh    - size or shape, is either n or (n,any_number)
NB.   A     - n×n-matrix, random general square invertible
NB.   P     - n×n-matrix, random Hermitian (symmetric)
NB.           positive definite; is defined as:
NB.             P := A * A^H
NB.
NB. Application:
NB. - make real symmetric positive definite 4×4-matrix P with
NB.   elements p having:
NB.     mantissa(p) ~ U(1,3)
NB.     exponent(p) ~ TN(0,4^2,_5,6)
NB.   :
NB.     P=. 1 3 0 4 _5 6&gemat pomat 4 4
NB. - make complex Hermitian positive definite 4×4-matrix P
NB.   with elements p having:
NB.     mantissa(Re(p)),mantissa(Im(p)) ~ U(1,3)
NB.     exponent(Re(p)),exponent(Im(p)) ~ TN(0,4^2,_5,6)
NB.   :
NB.     P=. 1 3 0 4 _5 6&(gemat j. gemat) pomat 4
NB. - make complex Hermitian positive definite 4×4-matrix P
NB.   with elements p having:
NB.     mantissa(Re(p)) ~ U(1,2)
NB.     exponent(Re(p)) ~ TN(0,1,-∞,+∞)
NB.     mantissa(Im(p)) ~ U(1,3)
NB.     exponent(Im(p)) ~ TN(0,4^2,_5,6)
NB.   :
NB.     P=. (1 2 0 1 __ _&gemat j. 1 3 0 4 _5 6&gemat) pomat 4
NB. - make real symmetric negative definite 4×4-matrix P
NB.   with elements p having:
NB.     mantissa(p) ~ U(1,3)
NB.     exponent(p) ~ TN(0,4^2,_5,6)
NB.   :
NB.     P=. -@(1 3 0 4 _5 6&gemat pomat) 4 4

pomat=: 1 : 'po_mt_@u@(2&$)'

NB. ---------------------------------------------------------
NB. ptmat
NB. ptmat2
NB.
NB. Description:
NB.   Adv. to make verb to make random Hermitian (symmetric)
NB.   positive definite tridiagonal matrix
NB.
NB. Syntax:
NB.   vapp=.  randx ptmat
NB.   vapp2=. randx ptmat2
NB. where
NB.   randx - monad to make random y-array A; is called as:
NB.             A=. randx y
NB.   vapp  - monad to make T; is called as:
NB.             T=. vapp sh
NB.   vapp2 - monad to make T2; is called as:
NB.             T2=. vapp2 sh
NB.   T     - n×n-matrix, random Hermitian (symmetric)
NB.           positive definite tridiagonal, is defined as:
NB.             T := L * L^H
NB.   T2    - n×n-matrix, random Hermitian (symmetric)
NB.           positive definite tridiagonal, is defined as:
NB.             T := K^_1
NB.   sh    - size or shape, is either n or (n,n)
NB.   L     - n×n-matrix, random lower bidiagonal with
NB.           positive diagonal entries
NB.   K     - n×n-matrix, random K(rho) generated from rho
NB.           parameter derived from randx such that
NB.             0.5 < abs(rho) < 1
NB.
NB. Application:
NB. - make real symmetric positive definite tridiagonal
NB.   4×4-matrix T with elements d and e having:
NB.     mantissa(d),mantissa(e) ~ U(1,3)
NB.     mantissa(d),exponent(e) ~ TN(0,4^2,_5,6)
NB.   :
NB.     T=. 1 3 0 4 _5 6&gemat ptmat 4 4
NB. - make complex Hermitian positive definite tridiagonal
NB.   4×4-matrix T with elements d and e having:
NB.     mantissa(d),mantissa(Re(e)),mantissa(Im(e)) ~ U(1,3)
NB.     exponent(d),exponent(Re(e)),exponent(Im(e)) ~ TN(0,4^2,_5,6)
NB.   :
NB.     T=. 1 3 0 4 _5 6&(gemat j. gemat) ptmat 4
NB. - make complex Hermitian positive definite tridiagonal
NB.   4×4-matrix T with elements d and e having:
NB.     mantissa(d),mantissa(Re(e)) ~ U(0,1)
NB.     mantissa(d),exponent(Re(e)) ~ TN(0,1,-∞,+∞)
NB.     mantissa(Im(e)) ~ U(1,3)
NB.     exponent(Im(e)) ~ TN(0,4^2,_5,6)
NB.   :
NB.     T=. gemat j. 1 3 0 4 _5 6&gemat ptmat 4
NB.
NB. Notes:
NB. - ptmat2 produces well-conditioned matrices for big n
NB.   (>400)
NB.
NB. TODO:
NB. - T should be sparse

ptmat=: 1 : 'po_mt_@(>:@| upddiag_mt_)@(((setdiag_mt_ diagmat_mt_)~ (}: ; _1:))/)@:u@(2,{.)'
ptmat2=: 1 : '(9&o. upddiag_mt_)@gtpick_mt_@%.@(u kmsmatrho_mt_ kmsmat_mt_)'

NB. ---------------------------------------------------------
NB. unmat
NB.
NB. Description:
NB.   Adv. to make verb to make random unitary (orthogonal)
NB.   matrix with distribution given by Haar measure
NB.
NB. Syntax:
NB.   vapp=. randnx unmat
NB. where
NB.   randnx - monad to make A, it is either randnf or
NB.            randnc; is called as:
NB.              A=. randnx y
NB.   vapp   - monad to make Q; is called as:
NB.              Q=. vapp sh
NB.   sh     - size or shape, either n or (n,n)
NB.   A      - n×n-matrix, non-singular, with elements
NB.            distributed as N(0,1)
NB.   Q      - n×n-matrix, random unitary (orthogonal)
NB.
NB. Formula:
NB.   A     := randnx(n,n)
NB.   (Q,R) := QR(A)
NB.   d     := diag(R)
NB.   Λ     := diagmat(d/|d|)
NB.   Q     := Q*Λ
NB. where
NB.   n - size of output matrix Q
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   I -: clean (mp ct) Q
NB.   I -: clean (mp~ct) Q
NB. where
NB.   I=. idmat n
NB.   Q=. randnf unmat n
NB.
NB. Application:
NB. - make real orthogonal 4×4-matrix Q, where Q is derived
NB.   via QR-factorization from real matrix B:
NB.     Q=. randnf unmat 4 4
NB. - make complex unitary 4×4-matrix U, where U is derived
NB.   via QR-factorization from complex matrix B:
NB.     U=. randnc unmat 4
NB.
NB. References:
NB. [1] Francesco Mezzadri. How to generate random matrices
NB.     from the classical compact groups. Notices of the
NB.     AMS, 2007, Vol. 54, pp. 592-604.
NB.     http://arxiv.org/abs/math-ph/0609050v2

unmat=: 1 : '(ungqr_mt_*"1*@diag_mt_)@geqrf_mt_@:u@(2&$)'

NB. ---------------------------------------------------------
NB. spmat
NB.
NB. Description:
NB.   Conj. to make verb to make random sparse array
NB.
NB. Syntax:
NB.   vapp=. randx spmat ratio
NB. where
NB.   randx - monad to make random y-array A; is called as:
NB.             A=. randx y
NB.   ratio - scalar in range [0,1], specifies desired
NB.           portion of non-zero elements
NB.   vapp  - monad to create S; is called as:
NB.             S=. vapp sh
NB.   S     - sh-array, sparsed randomly, being is zero
NB.           sh-array inhabited with values from A
NB.   sh    - r-vector of non-negative integers, a shape of S
NB.   r     ≥ 0, the rank of S
NB.
NB. Application:
NB. - make random real sparse 3×4-matrix S:
NB.     S=. (randnf spmat 0.25) 3 4
NB.
NB. TODO:
NB. - S should be sparse

spmat=: 2 : 'u@<.@(n * */) ((?@$~ #)~ (*/@$))} $&0'

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testtrmat
NB.
NB. Description:
NB.   Test trxxmat by matrix size given
NB.
NB. Syntax:
NB.   testtrmat sz
NB. where
NB.   sz - size or shape, is either n or (n,any_number)
NB.
NB. Notes:
NB. - result is not taken into account by benchmark

testtrmat=: 3 : 0
  (' randu           trl1mat' tmonad (]`]`(trl1con1@])`(_."_)`(_."_))) y
  (' randu           trlmat ' tmonad (]`]`(trlcon1 @])`(_."_)`(_."_))) y
  (' randu           tru1mat' tmonad (]`]`(tru1con1@])`(_."_)`(_."_))) y
  (' randu           trumat ' tmonad (]`]`(trucon1 @])`(_."_)`(_."_))) y

  ('(randu j. randu) trl1mat' tmonad (]`]`(trl1con1@])`(_."_)`(_."_))) y
  ('(randu j. randu) trlmat ' tmonad (]`]`(trlcon1 @])`(_."_)`(_."_))) y
  ('(randu j. randu) tru1mat' tmonad (]`]`(tru1con1@])`(_."_)`(_."_))) y
  ('(randu j. randu) trumat ' tmonad (]`]`(trucon1 @])`(_."_)`(_."_))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testgemat
NB.
NB. Description:
NB.   Test gemat by matrix shape given
NB.
NB. Syntax:
NB.   testgemat sh
NB. where
NB.   sh - r-vector of non-negative integers, the shape of
NB.        matrix
NB.
NB. Notes:
NB. - result is not taken into account by benchmark

testgemat=: 3 : 0
  ('gemat' tmonad (]`]`(gecon1@])`(_."_)`(_."_))) y
  EMPTY
)

NB. ---------------------------------------------------------
NB. testdimat
NB.
NB. Description:
NB.   Test dimat by matrix size given
NB.
NB. Syntax:
NB.   testdimat sz
NB. where
NB.   sz - size or shape, is either n or (n,n)
NB.
NB. Notes:
NB. - result is not taken into account by benchmark

testdimat=: 3 : 0
  (' gemat           dimat (randnf unmat)' tmonad (]`]`(gecon1@])`(_."_)`(_."_))) y
  (' gemat           dimat (randnc unmat)' tmonad (]`]`(gecon1@])`(_."_)`(_."_))) y

  ('(gemat j. gemat) dimat (randnc unmat)' tmonad (]`]`(gecon1@])`(_."_)`(_."_))) y
  ('(gemat j. gemat) dimat (randnc unmat)' tmonad (]`]`(gecon1@])`(_."_)`(_."_))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhemat
NB.
NB. Description:
NB.   Test hemat by matrix size given
NB.
NB. Syntax:
NB.   testhemat sz
NB. where
NB.   sz - size or shape, is either n or (n,any_number)
NB.
NB. Notes:
NB. - result is not taken into account by benchmark

testhemat=: 3 : 0
  (' randu           hemat' tmonad (]`]`(hecon1@])`(_."_)`(_."_))) y
  ('(randu j. randu) hemat' tmonad (]`]`(hecon1@])`(_."_)`(_."_))) y
  EMPTY
)

NB. ---------------------------------------------------------
NB. testpomat
NB.
NB. Description:
NB.   Test pomat by matrix size given
NB.
NB. Syntax:
NB.   testpomat sz
NB. where
NB.   sz - size or shape, is either n or (n,any_number)
NB.
NB. Notes:
NB. - result is not taken into account by benchmark

testpomat=: 3 : 0
  (' randu           pomat' tmonad (]`]`(pocon1@])`(_."_)`(_."_))) y
  ('(randu j. randu) pomat' tmonad (]`]`(pocon1@])`(_."_)`(_."_))) y
  EMPTY
)

NB. ---------------------------------------------------------
NB. testptmat
NB.
NB. Description:
NB.   Test ptmatx by matrix size given
NB.
NB. Syntax:
NB.   testptmat sz
NB. where
NB.   sz - size or shape, is either n or (n,any_number)
NB.
NB. Notes:
NB. - result is not taken into account by benchmark

testptmat=: 3 : 0
  (' randu           ptmat ' tmonad (]`]`(ptcon1@])`(_."_)`(_."_))) y
  ('(randu j. randu) ptmat ' tmonad (]`]`(ptcon1@])`(_."_)`(_."_))) y

  (' randu           ptmat2' tmonad (]`]`(ptcon1@])`(_."_)`(_."_))) y
  ('(randu j. randu) ptmat2' tmonad (]`]`(ptcon1@])`(_."_)`(_."_))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testspmat
NB.
NB. Description:
NB.   Test spmat by matrix shape given
NB.
NB. Syntax:
NB.   testspmat sh
NB. where
NB.   sh - r-vector of non-negative integers, the shape of
NB.        matrix
NB.
NB. Notes:
NB. - result is not taken into account by benchmark

testspmat=: 3 : 0
  (' randu           spmat 0.25' tmonad (]`]`(gecon1@])`(_."_)`(_."_))) y
  ('(randu j. randu) spmat 0.25' tmonad (]`]`(gecon1@])`(_."_)`(_."_))) y
  EMPTY
)

NB. ---------------------------------------------------------
NB. testrand
NB.
NB. Description:
NB.   Test xxxxmatx by matrix shape given
NB.
NB. Syntax:
NB.   testrand sh
NB. where
NB.   sh - r-vector of non-negative integers, the shape of
NB.        matrix

testrand=: EMPTY [ testspmat [ testptmat [ testpomat [ testhemat [ testdimat^:(=/) [ testgemat [ testtrmat
