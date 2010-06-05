NB. pjlap.ijs
NB. Pure J Linear Algebra Package

script_z_ '~system/packages/math/mathutil.ijs'  NB. mp
script_z_ '~system/packages/math/matutil.ijs'   NB. diag

require '~user/projects/pjlap/util.ijs'         NB. class-wide utility verbs

require '~user/projects/pjlap/bal.ijs'          NB. balance
require '~user/projects/pjlap/equ.ijs'          NB. equilibrate
require '~user/projects/pjlap/exp.ijs'          NB. exponent
require '~user/projects/pjlap/pow.ijs'          NB. powers
require '~user/projects/pjlap/rot.ijs'          NB. rotations
require '~user/projects/pjlap/sv.ijs'           NB. solve
require '~user/projects/pjlap/trf.ijs'          NB. LU factorization



init=: 3 : 0
  NB. machine epsilon, equivalent to:
  NB.   FP_EPSILON=: {. (-: ^: (1 ~: (1: + {:)) ^: _) 1 0.5
  FP_EPS=: +: 9!:18 ''

  NB. underflow threshold
  FP_UNFL=: {. (-: ^: (0 ~: {:) ^: _) 1 0.5

  NB. minimum exponent
  FP_UNFL_EXP=: 2 ^. FP_UNFL

  NB. overflow threshold
  FP_OVFL=: {: (+: ^: (_ ~: {.) ^: _) 1 0.5

  NB. maximum exponent
  NB. FIXME, indeed:
  NB.    +/ 2 ^ _53 {. i. 1024
  NB. 1.79769e308
   FP_OVFL_EXP=: 2 ^. FP_OVFL

  NB. safe minimum, such that 1/SAFMIN does not overflow
  FP_SAFMIN=: 2 ^ FP_OVFL_EXP <. - FP_UNFL_EXP

  NB. safe minimum, such that 1/SAFMIN does not overflow
  FP_SAFMIN=: 2 ^ FP_OVFL_EXP <. - FP_UNFL_EXP

   2 ^. 1.13687e_13 5.68434e_14
_43 _44
   NB. test
   1 = 1 + 1.13687e_13 5.68434e_14
0 1

   NB. test
   2 ^ _1074 _1075
4.94066e_324 0
   NB. safe minimum
   _ > 0.5 ^ _1023 _1024
1 0

   NB. test
   2 ^ 1023 1024
8.98847e307 _
   NB. safe maximum
   _ > % 0.5 ^ 1023 1024
1 0

  IBETA =  2
  IT =     53
  IRND =   5
  NGRD =   0
  MACHEP = -52
  NEGEP =  -53
  IEXP =   11
  MINEXP = -1022
  MAXEXP = 1024
  EPS =    2.22045e-16
  EPSNEG = 1.11022e-16
  XMIN =   2.22507e-308
  XMAX =   1.79769e+308

 ./testdlamch
  Epsilon                      =   1.110223024625157E-016
  Safe minimum                 =   2.225073858507201E-308
  Base                         =    2.00000000000000
  Precision                    =   2.220446049250313E-016
  Number of digits in mantissa =    53.0000000000000
  Rounding mode                =    1.00000000000000
  Minimum exponent             =   -1021.00000000000
  Underflow threshold          =   2.225073858507201E-308
  Largest exponent             =    1024.00000000000
  Overflow threshold           =   1.797693134862316E+308
  Reciprocal of safe minimum   =   4.494232837155790E+307
