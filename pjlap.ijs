NB. pjlap.ijs
NB. Pure J Linear Algebra Package

NB. =========================================================
NB. Include system verbs

script_z_ '~system/packages/math/mathutil.ijs'  NB. mp
script_z_ '~system/packages/math/matutil.ijs'   NB. diag

NB. =========================================================
NB. Configuration

coclass 'pjlap'

NB. ---------------------------------------------------------
NB. Initialize IEEE 754-1985 double-precision 64 bit floating
NB. point constants

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

NB. =========================================================
NB. Includes

require '~user/projects/pjlap/util.ijs'         NB. class-wide utility verbs
require '~user/projects/pjlap/bal.ijs'          NB. balance
require '~user/projects/pjlap/equ.ijs'          NB. equilibrate
require '~user/projects/pjlap/exp.ijs'          NB. exponent
require '~user/projects/pjlap/hess.ijs'         NB. Hessenberg reduction
require '~user/projects/pjlap/lu.ijs'           NB. LU factorization
require '~user/projects/pjlap/pow.ijs'          NB. powers
require '~user/projects/pjlap/rot.ijs'          NB. rotations
require '~user/projects/pjlap/sv.ijs'           NB. solve
