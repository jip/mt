NB. mt.ijs
NB. Matrix toolbox

NB. =========================================================
NB. Include system verbs

script_z_ '~system/packages/math/mathutil.ijs'        NB. mp
script_z_ '~system/packages/math/matutil.ijs'         NB. diag
script_z_ '~system/main/numeric.ijs'                  NB. range

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
NB. CPU-dependent constants

CPU_CACHE=: 512*1024                                  NB. CPU cache size in (7!:5) units

NB. ---------------------------------------------------------
NB. Block sizes for blocked algorithm versions

NB_POTRF=: 32                                         NB. see potrfNB in potrf.ijs
NB_GETRF=: 32                                         NB. see getrfNB in getrf.ijs

NB. =========================================================
NB. Includes

require '~user/projects/mt/util.ijs'                  NB. class-wide utility verbs
require '~user/projects/mt/bal.ijs'                   NB. balance
require '~user/projects/mt/chol.ijs'                  NB. Cholesky decomposition
require '~user/projects/mt/equ.ijs'                   NB. equilibrate
require '~user/projects/mt/exp.ijs'                   NB. exponent
require '~user/projects/mt/hess.ijs'                  NB. Hessenberg reduction
require '~user/projects/mt/lu.ijs'                    NB. LU factorization
require '~user/projects/mt/pow.ijs'                   NB. integer powers
require '~user/projects/mt/rot.ijs'                   NB. plane rotations
require '~user/projects/mt/sv.ijs'                    NB. solve systems A*X=B
