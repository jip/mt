NB. mt.ijs
NB. Matrix toolbox

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
FP_PREC=: FP_BASE * FP_EPS                            NB. machine precision = 2*ε
FP_EMIN=: 1 - FP_EBIAS                                NB. min exponent for normalized numbers = _1022
FP_UNFL=: FP_BASE ^ FP_EMIN                           NB. min normalized positive number = 2^_1022
FP_EMAX=: ((FP_BASE ^ FP_ELEN) - FP_BASE) - FP_EBIAS  NB. max exponent for normalized numbers = 1023
FP_OVFL=: (FP_BASE - FP_PREC) * (FP_BASE ^ FP_EMAX)   NB. max normalized positive number = (1-ε)*2^1024
FP_SFMIN=: FP_BASE ^ (FP_EMIN >. (- FP_EMAX))         NB. safe min, such that 1/SFMIN does not overflow

NB. ---------------------------------------------------------
NB. Optimal block size for blocked versions

BGETRF=: 32                                           NB. getrfl1u getrflu1 getrfu1l getrful1
BPOTRF=: 32                                           NB. potrfl potrfu

NB. ---------------------------------------------------------
NB. Optimal values
NB.   CPU_cache_size * algorithm_RAM_consumption_coefficient
NB. for recursive versions, in (7!:5) units

RGETRF=: 512*1024                                     NB. rgetrfl1u rgetrflu1 rgetrfu1l rgetrful1
RPOTRF=: 512*1024                                     NB. rpotrfl rpotrfu

NB. ---------------------------------------------------------
NB. Miscellaneous

VERBOSE=: 1                                           NB. boolean 'is verbose?' flag

NB. =========================================================
NB. Includes

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
NB.   ('mt';'z') copath 'base'
NB.   r=. gemat test 7 7

test=: 1 : 'u testtrf_mt_'
NB.---  NB. assert. 2 1 -: (# , #@$) y  NB. y must be 2-vector
NB.---  NB. (testbal , testequ , testexp , testhrd , testpow , testrot , testsv , testtrf) y
NB.---  u testtrf_mt_ y
NB.---)

NB. ---------------------------------------------------------
NB. block                                                   1
NB. Estimate optimal block size for each blocked version of
NB. algorithms
NB.
NB. Syntax:
NB.   r=. block m,n
NB. where
NB.   m,n - to estimate algorithms with random m×n matrices;
NB.         if m≠n then algorithms that accept square
NB.         matrices only are skipped
NB.   r   - boxed table with 4 columns: 'algorithm name'
NB.         'datatype' 'const name' 'const value', and with
NB.         rows per each unique pair: 'algorithm name'
NB.         'datatype'

block=: (3 : 0) " 1
  assert. 2 1 -: (# , #@$) y
  NB. (blockbal , blockequ , blockexp , blockhrd , blockpow , blockrot , blocksv , blocktrf) y
  blocktrf y
)

NB. ---------------------------------------------------------
NB. space                                                   1
NB. Estimate optimal value
NB.   CPU_cache_size * algorithm_RAM_consumption_coefficient
NB. in (7!:5) units for each blocked version of algorithms.
NB. Optimal value should guarantee that all input,
NB. intermediate and output data will fit entirely into CPU
NB. cache. This value is used by recursive versions to make
NB. decision whether to stop recursion or not.
NB.
NB. Syntax:
NB.   r=. space m,n
NB. where
NB.   m,n - to estimate algorithms with random m×n matrices;
NB.         if m≠n then algorithms that accept square
NB.         matrices only are skipped
NB.   r   - boxed table with 4 columns: 'algorithm name'
NB.         'datatype' 'const name' 'const value', and with
NB.         rows per each unique pair: 'algorithm name'
NB.         'datatype'
NB.
NB. Examples:
NB. - if there is no intermediate nouns, and operations are
NB.   in-place only, then coeff=1
NB. - if algorithm creates e.g. intermediate nouns of double
NB.   size as input, then coeff=1r3, since following must
NB.   hold: (insize + 2*insize ≤ CPU_cache_size)

space=: (3 : 0) " 1
  assert. 2 1 -: (# , #@$) y
  NB. (spacebal , spaceequ , spaceexp , spacehrd , spacepow , spacerot , spacesv , spacetrf) y
  spacetrf y
)
