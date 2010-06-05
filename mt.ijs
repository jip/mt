NB. Matrix toolbox
NB.
NB. test  Adverb to test algorithms
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

NB. =========================================================
NB. Configuration

coclass 'mt'

NB. ---------------------------------------------------------
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
FP_PREC=: FP_BASE * FP_EPS                            NB. machine precision = 2^_52
FP_EMIN=: 1 - FP_EBIAS                                NB. min exponent for normalized numbers = _1022
FP_UNFL=: FP_BASE ^ FP_EMIN                           NB. min normalized positive number = 2^_1022
FP_EMAX=: ((FP_BASE ^ FP_ELEN) - FP_BASE) - FP_EBIAS  NB. max exponent for normalized numbers = 1023
FP_OVFL=: (FP_BASE - FP_PREC) * (FP_BASE ^ FP_EMAX)   NB. max normalized positive number = (1-ε)*2^1024
FP_SFMIN=: FP_BASE ^ (FP_EMIN >. (- FP_EMAX))         NB. safe min, such that 1/SFMIN does not overflow

NB. ---------------------------------------------------------
NB. Miscellaneous

EMPTY=: i. 0 0

NB. =========================================================
NB. Includes

NB. ---------------------------------------------------------
NB. System verbs

script_z_ '~system/main/numeric.ijs'                  NB. range
script_z_ '~system/main/myutil.ijs'                   NB. timespacex
script_z_ '~system/packages/math/mathutil.ijs'        NB. mp
script_z_ '~system/packages/math/matutil.ijs'         NB. diag

NB. ---------------------------------------------------------
NB. Package verbs and nouns

NB. utility
require '~user/projects/mt/util.ijs'    NB. utilities
require '~user/projects/mt/con.ijs'     NB. condition number
require '~user/projects/mt/rand.ijs'    NB. random objects
require '~user/projects/mt/struct.ijs'  NB. struct handlers

NB. low-level
require '~user/projects/mt/bal.ijs'     NB. balance
require '~user/projects/mt/equ.ijs'     NB. equilibrate
require '~user/projects/mt/ref.ijs'     NB. reflections
require '~user/projects/mt/rot.ijs'     NB. rotations

NB. mid-level
require '~user/projects/mt/hrd.ijs'     NB. Hessenberg reduction
require '~user/projects/mt/orf.ijs'     NB. orthogonal factorization
require '~user/projects/mt/pow.ijs'     NB. integer powers
require '~user/projects/mt/trf.ijs'     NB. triangular factorization
require '~user/projects/mt/tri.ijs'     NB. inverse
require '~user/projects/mt/trs.ijs'     NB. solve linear monomial equation via triang. factorization

NB. hi-level
require '~user/projects/mt/exp.ijs'     NB. exponent
require '~user/projects/mt/log.ijs'     NB. logarithm

NB. =========================================================
NB. Interface

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. test
NB. Adverb to test algorithms
NB.
NB. Syntax:
NB.   r=. mkge test m,n
NB. where
NB.   m,n  - 2-vector of integers, shape of random matrices
NB.          to test algorithms; only algorithms which accept
NB.          m and n given are called
NB.   mkge - monadic verb to generate random non-singular
NB.          general y-matrix (shape is taken from y)
NB.   r    - boxed table with 6 columns:
NB.          - algorithm name
NB.          - the estimated reciprocal of the condition
NB.            number of a matrix in 1-norm
NB.          - relative backward error
NB.          - relative forward error
NB.          - execution time, sec.
NB.          - execution space, bytes
NB.
NB. Application:
NB. - with limited random matrix values' amplitudes
NB.   load '~addons/math/mt/mt.ijs'
NB.   cocurrent 'mt'
NB.   r=. (_1 1 0 16 _6 4 & gemat) test 500 500
NB.   r=. (_1 1 0 16 _6 4 & (gemat j. gemat)) test 500 500

test=: 1 : 0

  require 'printf'
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'gesv getrf potrf'

  '%-25s %-12s %-12s %-12s %-12s %-12s' & printf 'algorithm' ; 'rcond' ; 'rel bw err' ; 'rel fw err' ; 'time, sec.' ; 'space, bytes'
  assert. 2 1 -: (# , #@$) y  NB. y must be 2-vector
  ((u testtrf) [ (u testtri) [ (u testtrs)) y
)
