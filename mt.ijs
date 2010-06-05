NB. 'Matrix toolbox' addon's entry point
NB.
NB. TESTLOGFILE  Log file name or a: to switch logging off
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
NB. EMPTY        i. 0 0
NB.
NB. test         Adv. to make verb to test algorithms by
NB.              matrix of generator and shape given
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. User config

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Tests logging

TESTLOGFILE=: < jpath '~temp/mt.log'                  NB. assign a: to switch off file logging
TESTLOG=: ''                                          NB. literal array, being formatted test log

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Debug level used by "dbg":
NB.   0 - execute debuging verb transparently and silently
NB.   1 - show for debuging verb its rank and valency,
NB.       input's and output's shapes
NB.   2 - case 1 plus input's and output's values

DEBUG=: 2

NB. ---------------------------------------------------------
NB. System config

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. IEEE 754-1985 double-precision 64 bit floating point
NB. constants

NB. basic values
FP_BASE=: 2                                           NB. floating point base
FP_ELEN=: 11                                          NB. exponent field length (bits)
FP_FLEN=: 53                                          NB. fraction field length (bits)
FP_IGUNFL=: 1                                         NB. is gradual underflow? (boolean)

NB. derived values
FP_EBIAS=: (FP_BASE ^ (FP_ELEN - 1)) - 1              NB. exponent bias for normalized numbers = 1023
FP_EPS=: FP_BASE ^ (- FP_FLEN)                        NB. machine epsilon ε = 2^_53
FP_PREC=: FP_BASE * FP_EPS                            NB. machine precision β = 2^_52
FP_EMIN=: 1 - FP_EBIAS                                NB. min exponent for normalized numbers = _1022
FP_UNFL=: FP_BASE ^ FP_EMIN                           NB. min normalized positive number = 2^_1022
FP_EMAX=: ((FP_BASE ^ FP_ELEN) - FP_BASE) - FP_EBIAS  NB. max exponent for normalized numbers = 1023
FP_OVFL=: (FP_BASE - FP_PREC) * (FP_BASE ^ FP_EMAX)   NB. max normalized positive number = (1-ε)*2^1024
FP_SFMIN=: FP_BASE ^ (FP_EMIN >. (- FP_EMAX))         NB. safe min, such that 1/FP_SFMIN does not overflow

NB. ---------------------------------------------------------
NB. Constants

EMPTY=: i. 0 0

NB. =========================================================
NB. Includes

NB. ---------------------------------------------------------
NB. System verbs

script_z_ '~system/main/printf.ijs'             NB. printf vsprintf
script_z_ '~system/main/myutil.ijs'             NB. timespacex
script_z_ '~system/packages/math/mathutil.ijs'  NB. mp

NB. ---------------------------------------------------------
NB. Addon verbs and nouns
NB.
NB. TODO: s@user/projects@addons/math/mt@g

NB. utilities
require '~user/projects/mt/dbg.ijs'     NB. Debug
require '~user/projects/mt/util.ijs'    NB. Utilities
require '~user/projects/mt/ios.ijs'     NB. IOS
require '~user/projects/mt/norm.ijs'    NB. Norms
require '~user/projects/mt/struct.ijs'  NB. Structure handlers
require '~user/projects/mt/rand.ijs'    NB. Random objects
require '~user/projects/mt/con.ijs'     NB. Condition number

NB. low-level
require '~user/projects/mt/bak.ijs'     NB. Recover eigenvectors after balancing
require '~user/projects/mt/bal.ijs'     NB. Balance
require '~user/projects/mt/equ.ijs'     NB. Equilibrate
require '~user/projects/mt/ref.ijs'     NB. Reflections
require '~user/projects/mt/rot.ijs'     NB. Rotations
require '~user/projects/mt/gq.ijs'      NB. Generate Q from its factored form
require '~user/projects/mt/mq.ijs'      NB. Multiply by Q represented in factored form
require '~user/projects/mt/sm.ijs'      NB. Solve linear monomial equation with triangular matrix

NB. mid-level
require '~user/projects/mt/fri.ijs'     NB. Inverse using the Frobenius formula
require '~user/projects/mt/hrd.ijs'     NB. Hessenberg reduction
require '~user/projects/mt/qf.ijs'      NB. Orthogonal factorization
require '~user/projects/mt/trf.ijs'     NB. Triangular factorization
require '~user/projects/mt/tri.ijs'     NB. Inverse by trf
require '~user/projects/mt/trs.ijs'     NB. Solve linear monomial equation from trf

NB. hi-level
require '~user/projects/mt/ev.ijs'      NB. eigenvalue decomposition
require '~user/projects/mt/exp.ijs'     NB. exponent
require '~user/projects/mt/log.ijs'     NB. logarithm
require '~user/projects/mt/pow.ijs'     NB. power
require '~user/projects/mt/sv.ijs'      NB. solve linear monomial equation

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. test
NB.
NB. Description:
NB.   Adv. to make verb to test algorithms by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkge test
NB. where
NB.   (m,n) - 2-vector of integers, shape of random matrices
NB.           to test algorithms; only algorithms which
NB.           accept m and n given will be tested
NB.   mkge  - monadic verb to generate random non-singular
NB.           general y-matrix (shape is taken from y)
NB.   vtest - verb to test algorithms; is called as:
NB.              vtest (m,n)
NB.
NB. Application:
NB. - test by random square real matrix with limited values'
NB.   amplitudes:
NB.     (_1 1 0 16 _6 4 & gemat_mt_) test_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) test_mt_ 150 200

test=: 1 : 0
  '%-25s %-12s %-12s %-12s %-12s %-12s' printf 'algorithm' ; 'rcond' ; 'rel fwd err' ; 'rel bwd err' ; 'time, sec.' ; 'space, bytes'

  NB. low-level algorithms
  (u testbak_mt_) y
  (u testbal_mt_) y
  (u testequ_mt_) y
  (u testref_mt_) y
  (u testrot_mt_) y
  (u testgq_mt_) y
  (u testmq_mt_) y
  (u testsm_mt_) y

  NB. mid-level algorithms
  (u testfri_mt_) y
  (u testhrd_mt_) y
  (u testqf_mt_) y
  (u testtrf_mt_) y
  (u testtri_mt_) y
  (u testtrs_mt_) y

  NB. hi-level algorithms
  (u testexp_mt_) y
  (u testlog_mt_) y
  (u testpow_mt_) y
  (u testsv_mt_) y

  EMPTY_mt_
)
