NB. 'Matrix toolbox' addon's entry point
NB.
NB. DEBUG      Debug level
NB. FP_BASE    Floating point base
NB. FP_ELEN    Exponent field length (bits)
NB. FP_FLEN    Fraction field length (bits)
NB. FP_IGUNFL  Is gradual underflow? (boolean)
NB. FP_EBIAS   Exponent bias for normalized numbers
NB. FP_EPS     Machine epsilon
NB. FP_PREC    Machine precision
NB. FP_EMIN    Min exponent for normalized numbers
NB. FP_UNFL    Min normalized positive number
NB. FP_EMAX    Max exponent for normalized numbers
NB. FP_OVFL    Max normalized positive number
NB. FP_SFMIN   Safe min, such that 1/FP_SFMIN does not
NB.            overflow
NB.
NB. testlow    Adv. to make verb to test low-level
NB.            algorithms by matrix of generator and shape
NB.            given
NB. testmid    Adv. to make verb to test mid-level algorithms
NB.            by matrix of generator and shape given
NB. testhigh   Adv. to make verb to test high-level
NB.            algorithms by matrix of generator and shape
NB.            given
NB. test       Adv. to make verb to test algorithms by matrix
NB.            of generator and shape given
NB.
NB. verify     Ambivalent predicate to verify mt addon
NB.
NB. Copyright 2010,2011,2013,2017,2018,2020,2021,2023,2024
NB.           Igor Zhuravlov
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

NB. =========================================================
NB. Concepts
NB.
NB. Terms:
NB.   modifier  - either adverb or conjunction
NB.   actor     - either verb or modifier
NB.   function  - a type of verb which returns useful result,
NB.               usually is called as:
NB.                 result=. function args
NB.   procedure - a type of verb which returns useless
NB.               result, is opposite to function, usually is
NB.               called as anyone of:
NB.                 procedure args
NB.                 trash=. procedure args
NB.                 EMPTY [ procedure args
NB.   predicate - function returning boolean
NB.   semipredicate
NB.             - function returning either boolean or NULL
NB.   identity  - a verb returning its argument[s]
NB.   arity     - the number of arguments taken by a verb
NB.   nilad, niladic verb
NB.             - 0-ary verb, usually is called as:
NB.                 out=. nilad ''
NB.   monad, monadic verb
NB.             - 1-ary verb
NB.   dyad, dyadic verb
NB.             - 2-ary verb
NB.   ambivalent
NB.             - either monadic or dyadic
NB.   debug     - execute a verb and show debug info obtained
NB.   test      - execute a verb with random arguments
NB.               supplied, it's aimed to:
NB.               - check whether the verb execution was
NB.                 succeed
NB.               - estimate a reciprocal of condition number
NB.                 of input if it was a square non-singular
NB.                 matrix
NB.               - measure a relative forward error if verb
NB.                 is a solver
NB.               - measure a relative backward error
NB.               - benchmark a time and space required for
NB.                 execution
NB.   verify    - execute a verb with thoroughly selected
NB.               arguments to check the result correctness
NB.
NB. Notation:
NB.   ∞          is       ±∞ i.e. either -∞ or +∞
NB.   a,b = ∞    means    a,b ∈ {-∞,+∞}, not a = b = ±∞
NB.
NB. Conventions:
NB. 1) a result returned from a test actor is an inverted
NB.    table whose format is specified in test.ijs

NB. =========================================================
NB. Configuration

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. invertible ambivalent identity to erase names
NB. syntax:
NB.   out=. [x] f&.erasen y
NB. to erase global names created while f was executed

erasen=: ([ 4!:5@1) :. ([ (4!:5@0)@erase@(4!:5@1))

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. User config

NB. Debug level used by dbg conj., the atom:
NB.   0 - execute debuging verb transparently and silently
NB.   1 - show for debuging verb its rank and valency,
NB.       input's and output's shapes
NB.   2 - case (1) plus input's and output's values

DEBUG=: 2

NB. ---------------------------------------------------------
NB. System config

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. IEEE 754-1985 floating point constants for single
NB. (32-bit) and double (64-bit) precision

NB. basic values
FP_BASE=: 2                                              NB. floating point base
FP_ELEN=: IF64 {  8 11                                   NB. exponent field length (bits)
FP_FLEN=: IF64 { 24 53                                   NB. fraction field length (bits)
FP_IGUNFL=: 1                                            NB. is gradual underflow? (boolean)

NB. derived values
FP_EBIAS=: <. (FP_BASE ^ (FP_ELEN - 1)) - 1              NB. exponent bias for normalized numbers (127 or 1023)
FP_EPS=: FP_BASE ^ (- FP_FLEN)                           NB. machine epsilon ε (2^_24 or 2^_53)
FP_PREC=: FP_BASE * FP_EPS                               NB. machine precision β (2^_23 or 2^_52)
FP_EMIN=: <. 1 - FP_EBIAS                                NB. min exponent for normalized numbers (_126 or _1022)
FP_UNFL=: FP_BASE ^ FP_EMIN                              NB. min normalized positive number (2^_126 or 2^_1022)
FP_EMAX=: <. ((FP_BASE ^ FP_ELEN) - FP_BASE) - FP_EBIAS  NB. max exponent for normalized numbers (127 or 1023)
FP_OVFL=: (FP_BASE - FP_PREC) * (FP_BASE ^ FP_EMAX)      NB. max normalized positive number ((1-ε)*2^128 or (1-ε)*2^1024)
FP_SFMIN=: FP_BASE ^ (FP_EMIN >. (- FP_EMAX))            NB. safe min, such that 1/FP_SFMIN does not overflow

NB. =========================================================
NB. Includes

NB. ---------------------------------------------------------
NB. System definitions

require 'math/misc/mathutil'   NB. mp_mt_
require 'general/misc/format'  NB. clipfmt_z_
require 'stats/base'           NB. comb_z_ combrep_z_ perm_z_

NB. ---------------------------------------------------------
NB. Addon definitions

NB. utilities
require 'math/mt/dbg'          NB. Debug
require 'math/mt/fork'         NB. Extended forks
require 'math/mt/util'         NB. Utilities
require 'math/mt/iso'          NB. ISO
require 'math/mt/mm'           NB. Matrix Market exchange formats converter
require 'math/mt/norm'         NB. Norms
require 'math/mt/quatern'      NB. Quaternions
require 'math/mt/struct'       NB. Structure handlers
require 'math/mt/rand'         NB. Random arrays
require 'math/mt/test'         NB. Test

NB. low-level
require 'math/mt/basic'        NB. Basic linear algebra operations
require 'math/mt/bak'          NB. Restore original eigenvectors
require 'math/mt/bal'          NB. Balance
require 'math/mt/cond'         NB. Condition number
require 'math/mt/ref'          NB. Reflection
require 'math/mt/rot'          NB. Rotation
require 'math/mt/gq'           NB. Generate Q from its factored form
require 'math/mt/mq'           NB. Multiply by Q represented in factored form
require 'math/mt/scl'          NB. Scale

NB. mid-level
require 'math/mt/eq'           NB. Eigenvalues and Schur form
require 'math/mt/evc'          NB. Eigenvectors
require 'math/mt/hrd'          NB. Hessenberg reduction
require 'math/mt/pf'           NB. Orthogonal factorization with pivoting
require 'math/mt/qf'           NB. Orthogonal factorization
require 'math/mt/trf'          NB. Triangular factorization
require 'math/mt/tri'          NB. Inverse by trf
require 'math/mt/trs'          NB. Solve linear monomial equation by trf

NB. high-level
require 'math/mt/ev'           NB. Eigenvalues and eigenvectors
require 'math/mt/exp'          NB. Matrix exponential
require 'math/mt/ls'           NB. Solve overdetermined or underdetermined linear monomial equation
require 'math/mt/pow'          NB. Raise matrix to integer power(s)
require 'math/mt/sv'           NB. Solve linear monomial equation

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testlow
NB. testmid
NB. testhigh
NB. test
NB.
NB. Description:
NB.   Adv. to make verb to test algorithms either all or not,
NB.   by matrix of generator and shape given
NB.
NB. Syntax:
NB.   log=. (mkge testxxxx) (m,n)
NB. where
NB.   mkge  - monad to generate random non-singular general
NB.           y-matrix (shape is taken from y)
NB.   (m,n) - 2-vector of integers, shape of random matrices
NB.           to test algorithms; only algorithms which
NB.           accept m and n given will be tested
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test low-level algorithms by random square integer
NB.   matrix with elements distributed uniformly with support
NB.   [0,100):
NB.     log=. ?@$&100 testlow_mt_ 10 10
NB. - test mid-level algorithms by random rectangular real
NB.   matrix with elements distributed uniformly with support
NB.   (0,1):
NB.     log=. ?@$&0 testmid_mt_ 200 150
NB. - test high-level algorithms by random square real matrix
NB.   with elements with limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testhigh_mt_ 200 200
NB. - test all algorithms by random rectangular complex
NB.   matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) test_mt_ 150 200

testlow=: 1 : '(u testmq_mt_) ,&.>~ (u testgq_mt_) ,&.>~ (u testrot_mt_) ,&.>~ (u testref_mt_) ,&.>~ (u testquatern_mt_) ,&.>~ (u testnorm_mt_) ,&.>~ (u testcon_mt_) ,&.>~ (u testbal_mt_) ,&.>~ (u testbak_mt_) ,&.>~ (u testbasic_mt_) ,&.>~ testrand_mt_'

testmid=: 1 : '(u testtrs_mt_) ,&.>~ (u testtri_mt_) ,&.>~ (u testtrf_mt_) ,&.>~ (u testqf_mt_) ,&.>~ (u testpf_mt_) ,&.>~ (u testhrd_mt_) ,&.>~ (u testevc_mt_) ,&.>~ (u testeq_mt_)'

testhigh=: 1 : '(u testmm_mt_) ,&.>~ (u testls_mt_) ,&.>~ (u testsv_mt_) ,&.>~ (u testpow_mt_) ,&.>~ (u testexp_mt_) ,&.>~ (u testev_mt_)'

test=: 1 : '(u testhigh_mt_) ,&.>~ (u testmid_mt_) ,&.>~ (u testlow_mt_)'

NB. =========================================================
NB. Verification suite

NB. ---------------------------------------------------------
NB. verify
NB.
NB. Description:
NB.   Ambivalent predicate to verify mt addon and to stop at
NB.   1st failed
NB.
NB. Syntax:
NB.   isSucceed=. [isVerbose] verify ''
NB. where
NB.   isVerbose - boolean to display assertions probed with
NB.               the result of its execution, optional,
NB.               default is 0 to run silently
NB.   isSucceed - boolean, verification result
NB.
NB. Application:
NB. - if verification is failed then re-run in verbose mode
NB.   to see the 1st assertion failed:
NB.      verify_mt_ ''    NB. run silently
NB.   0                   NB. some assertions are failed
NB.      1 verify_mt_ ''  NB. run verbosely to display each
NB.   ...                 NB.   assertion's result

verify=: 0&$: :(4 : '*./ 0!:(x { 3 2)&.erasen 1 dir ''~addons/math/mt/verify/*.ijs''')
