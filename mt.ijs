NB. 'Matrix toolbox' addon's entry point
NB.
NB. TESTLOGFILE  a: to switch logging off or boxed logfile
NB.              name
NB. TESTLOG      Literal array, being formatted test log
NB. DEBUG        Debug level
NB. FP_BASE      Floating point base
NB. FP_ELEN      Exponent field length (bits)
NB. FP_FLEN      Fraction field length (bits)
NB. FP_IGUNFL    Is gradual underflow? (boolean)
NB. FP_EBIAS     Exponent bias for normalized numbers
NB. FP_EPS       Machine epsilon
NB. FP_PREC      Machine precision
NB. FP_EMIN      Min exponent for normalized numbers
NB. FP_UNFL      Min normalized positive number
NB. FP_EMAX      Max exponent for normalized numbers
NB. FP_OVFL      Max normalized positive number
NB. FP_SFMIN     Safe min, such that 1/FP_SFMIN does not
NB.              overflow
NB.
NB. testlow      Adv. to make verb to test low-level
NB.              algorithms by matrix of generator and shape
NB.              given
NB. testmid      Adv. to make verb to test mid-level
NB.              algorithms by matrix of generator and shape
NB.              given
NB. testhigh     Adv. to make verb to test high-level
NB.              algorithms by matrix of generator and shape
NB.              given
NB. test         Adv. to make verb to test algorithms by
NB.              matrix of generator and shape given
NB. verify       Nilad to verify mt, output result to console
NB.              and return it
NB.
NB. Version: 0.13.3 2021-06-29
NB.
NB. Copyright 2010-2021 Igor Zhuravlov
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
NB. Interface

NB. ---------------------------------------------------------
NB. User config

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Tests logging

TESTLOGFILE=: < jpath '~temp/mt.log'  NB. assign a: to switch off file logging
TESTLOG=: ''                          NB. literal array, being formatted test log

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Debug level used by dbg conj., the atom:
NB.   0 - execute debuging verb transparently and silently
NB.   1 - show for debuging verb its rank and valency,
NB.       input's and output's shapes
NB.   2 - case (1) plus input's and output's values

DEBUG=: 2

NB. ---------------------------------------------------------
NB. System config

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. IEEE 754-1985 double-precision 64 bit floating point
NB. constants

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
require 'math/mt/benchmark'    NB. Benchmark

NB. low-level
require 'math/mt/bak'          NB. Restore original eigenvectors
require 'math/mt/bal'          NB. Balance
require 'math/mt/cond'         NB. Condition number
require 'math/mt/ref'          NB. Reflection
require 'math/mt/rot'          NB. Rotation
require 'math/mt/gq'           NB. Generate Q from its factored form
require 'math/mt/mq'           NB. Multiply by Q represented in factored form
require 'math/mt/scl'          NB. Scale
require 'math/mt/sm'           NB. Solve linear monomial equation with triangular matrix

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
require 'math/mt/pow'          NB. Raise matrix to integer power[s]
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
NB.   vtest=. mkge testxxxx
NB. where
NB.   (m,n) - 2-vector of integers, shape of random matrices
NB.           to test algorithms; only algorithms which
NB.           accept m and n given will be tested
NB.   mkge  - monad to generate random non-singular general
NB.           y-matrix (shape is taken from y)
NB.   vtest - verb to test algorithms; is called as:
NB.             vtest (m,n)
NB.
NB. Application:
NB. - test low-level algorithms by random square integer
NB.   matrix with elements distributed uniformly with support
NB.   [0,100):
NB.    ?@$&100 testlow_mt_ 10 10
NB. - test mid-level algorithms by random rectangular real
NB.   matrix with elements distributed uniformly with support
NB.   (0,1):
NB.     ?@$&0 testmid_mt_ 200 150
NB. - test high-level algorithms by random square real matrix
NB.   with elements with limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testhigh_mt_ 200 200
NB. - test all algorithms by random rectangular complex
NB.   matrix:
NB.     (gemat_mt_ j. gemat_mt_) test_mt_ 150 200

testlow=: 1 : 0
     testrand_mt_  y
  (u testbak_mt_ ) y  NB. square matrices only
  (u testbal_mt_ ) y  NB. square matrices only
  (u testref_mt_ ) y  NB. matrices with min dimention ≤ 200 only
  (u testrot_mt_ ) y  NB. matrix of shape (2 1:} y) is used
  (u testgq_mt_  ) y
  (u testmq_mt_  ) y
  (u testsm_mt_  ) y  NB. square matrices with size ≤ 500 only

  EMPTY
)

testmid=: 1 : 0
  (u testeq_mt_  ) y  NB. square matrices only
  (u testevc_mt_ ) y  NB. square matrices only
  (u testhrd_mt_ ) y  NB. square matrices only
  (u testpf_mt_  ) y
  (u testqf_mt_  ) y
  (u testtrf_mt_ ) y
  (u testtri_mt_ ) y  NB. square matrices only
  (u testtrs_mt_ ) y  NB. square matrices only

  EMPTY
)

testhigh=: 1 : 0
  (u testev_mt_  ) y  NB. square matrices only
  (u testexp_mt_ ) y  NB. square matrices only
  (u testpow_mt_ ) y  NB. square matrices only
  (u testsv_mt_  ) y  NB. square matrices only
  (u testls_mt_  ) y
  (u testmm_mt_  ) y

  EMPTY
)

test=: 1 : 0
  echo fmtlog_mt_ 'sentence';'rcond';'rel fwd err';'rel bwd err';'time, sec.';'space, bytes'

  (u testlow_mt_ ) y  NB. low-level algorithms
  (u testmid_mt_ ) y  NB. mid-level algorithms
  (u testhigh_mt_) y  NB. high-level algorithms

  EMPTY
)

NB. =========================================================
NB. Verification suite

NB. ---------------------------------------------------------
NB. verify
NB.
NB. Description:
NB.   Nilad to verify mt, output result to console and return
NB.   it
NB.
NB. Syntax:
NB.   'probed failed'=. verify_mt_ ''
NB. where
NB.   probed ≥ 0, tests probed counter
NB.   failed ≥ 0, tests failed counter

verify=: verifymm
