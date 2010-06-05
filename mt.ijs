NB. mt.ijs
NB. Matrix toolbox

NB. =========================================================
NB. Configuration

coclass 'mt'
NB.coinsert   'base'

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
FP_PREC=: FP_BASE * FP_EPS                            NB. machine precision = 2*ε
FP_EMIN=: 1 - FP_EBIAS                                NB. min exponent for normalized numbers = _1022
FP_UNFL=: FP_BASE ^ FP_EMIN                           NB. min normalized positive number = 2^_1022
FP_EMAX=: ((FP_BASE ^ FP_ELEN) - FP_BASE) - FP_EBIAS  NB. max exponent for normalized numbers = 1023
FP_OVFL=: (FP_BASE - FP_PREC) * (FP_BASE ^ FP_EMAX)   NB. max normalized positive number = (1-ε)*2^1024
FP_SFMIN=: FP_BASE ^ (FP_EMIN >. (- FP_EMAX))         NB. safe min, such that 1/SFMIN does not overflow

NB. ---------------------------------------------------------
NB. Miscellaneous

VERBOSE=: 1                                           NB. boolean 'is verbose?' flag

NB. =========================================================
NB. Includes

NB. ('mt';'z') copath 'base'

NB. ---------------------------------------------------------
NB. System verbs

script_z_ '~system/main/numeric.ijs'                  NB. range
script_z_ '~system/main/myutil.ijs'                   NB. timespacex
script_z_ '~system/packages/math/mathutil.ijs'        NB. mp
script_z_ '~system/packages/math/matutil.ijs'         NB. diag
script_z_ '~system/packages/stats/random.ijs'         NB. rand01
script_z_ '~system/packages/stats/statdist.ijs'       NB. normalrand

NB. ---------------------------------------------------------
NB. Package verbs

require '~user/projects/mt/util.ijs'                  NB. utilities
require '~user/projects/mt/struct.ijs'                NB. struct handlers
require '~user/projects/mt/bal.ijs'                   NB. balance
require '~user/projects/mt/equ.ijs'                   NB. equilibrate
require '~user/projects/mt/exp.ijs'                   NB. exponent
require '~user/projects/mt/hrd.ijs'                   NB. Hessenberg reduction
NB. require '~user/projects/mt/orf.ijs'                   NB. orthogonal factorization (LQ QL QR RQ)
require '~user/projects/mt/pow.ijs'                   NB. integer powers
require '~user/projects/mt/rand.ijs'                  NB. random objects
require '~user/projects/mt/rcond.ijs'                 NB. reciprocal of condition number
require '~user/projects/mt/rot.ijs'                   NB. plane rotations
require '~user/projects/mt/sv.ijs'                    NB. solve linear monomial equations
require '~user/projects/mt/trf.ijs'                   NB. triangular factorization (Cholesky LU)

NB. =========================================================
NB. Interface

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. test                                                    1
NB. Adverb to test algorithms
NB.
NB. Syntax:
NB.   r=. mkge test m,n
NB. where
NB.   m,n  - 2-vector of integers, shape of random matrices
NB.          to test algorithms; if m≠n then algorithms that
NB.          accept square matrices only are skipped
NB.   mkge - monadic verb to generate random non-singular
NB.          general y-matrix (shape is taken from y)
NB.   r    - boxed table with 4 columns: 'algorithm name'
NB.          'error' 'time, sec.' 'space, bytes'
NB.
NB. Application:
NB.   cocurrent 'mt'
NB.   r=. (_1 1 0 16 _6 4 & gemat) test 132 132
NB.   r=. (_1 1 0 16 _6 4 & (gemat j. gemat)) test 132 132

NB.--- test=: 1 : 'u testtrf'
test=: 1 : 0

  require 'printf'
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'gesv getrf potrf'

  '%-25s %-12s %-12s %-12s %-12s' & printf ^: (VERBOSE"_) 'Algorithm' ; 'Backward err' ; 'Forward err' ; 'Time, sec.' ; 'Space, bytes'
  assert. 2 1 -: (# , #@$) y  NB. y must be 2-vector
  ((u testtrf) , (u testtrs)) y
)

