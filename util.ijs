NB. util.ijs
NB. Common mt verbs
NB.
NB. vnormi  inf-norm of vector y
NB. norm1   1-norm of table or vector y
NB. trace   matrix trace
NB. ut2tr   transform upper triangular matrix to packed form
NB. sdiag   add element[s from] x to diagonal of matrix y

coclass 'mt'

NB. =========================================================
NB. Local verbs

NB. =========================================================
NB. Interface verbs

vnormi=: >./ @: |                       NB. inf-norm of vector y
norm1=: >./ @: (+/) @: |                NB. 1-norm of table or vector y
trace=: +/ @ diag                       NB. matrix trace
ut2tr=: ((I. @ , @ (<:/~@i.) @ #) { ,)  NB. transform N×N upper triangular matrix to (N*(N+1)/2)-vector, packed form
lio=: + ` (* i.)/ " 1                   NB. integers grid (2{y) steps from (0{y) by (1{y)
ii2cp=: < @ (, ` (, @ ]) @. =)          NB. make cycle permutation from indices x and y
ht2i=: (] + (i. @ -)) " 0 0             NB. form int list from head y to tail x: y (y+1) ... (x-1)
abssq=: +/ @: *: @ +.                   NB. Re(y)^2 + Im(y)^2
ct=: + @ |:                             NB. conjugate transpose
ts=: 6!:2 , 7!:2 @ ]                    NB. execution time and space
UL2L=: * ((>:/ & i.)/ @ $)              NB. fill strict upper triangular by zeros
UL2SL=: * ((>/ & i.)/ @ $)              NB. fill upper triangular by zeros
UL2U=: * ((<:/ & i.)/ @ $)              NB. fill strict lower triangular by zeros
UL2SU=: * ((</ & i.)/ @ $)              NB. fill lower triangular by zeros
UL2I=: (=/ & i.)/ @ $                   NB. identity matrix with the same shape as y
LU2L1=: UL2I + UL2SL                    NB. extract strict L is m×min(m,n) matrix
>>>>>>>>
LU2L1=: (mn & ({. " 1)) @: (I + L & *)  NB. L is m×min(m,n) matrix
LU2U=. (mn & {.) @: (U & *)            NB. U is min(m,n)×n matrix

NB. ---------------------------------------------------------
NB. kernel3
NB. Conjunction to form top level computing kernel verbs to
NB. calculate:
NB.   submatrixA op0 submatrixB op1 submatrixC
NB. from x and matrix y
NB.
NB. Syntax:
NB.   kernel3verb=. (iosA ` iosB ` iosC) kernel3 (op0 ` op1)
NB. where
NB.   iosA iosB iosC - ambivalent verbs to embed into trains
NB.                    (v ] ;. 0 ]) which extract corresp.
NB.                    submatrix using arguments x and
NB.                    matrix y
NB.   op0 op1        - dyadic verbs to combine submatrices
NB.                    A B C
NB.
NB. If:
NB.   mkiosA=. 0 { [
NB.   mkiosB=. 1 { [
NB.   mkiosC=. 2 { [
NB.   kabc=. mkiosA ` mkiosB ` mkiosC kernel3 (- ` mp)
NB.   iosABC=. 3 2 2 $ 2 2 3 3 2 0 3 2 0 2 2 3
NB.   m=. i. 5 5
NB.   A=. (0 { iosABC) (] ;. 0 ]) m
NB.   B=. (1 { iosABC) (] ;. 0 ]) m
NB.   C=. (2 { iosABC) (] ;. 0 ]) m
NB.   updatedA=. iosABC kabc m
NB. then
NB.   updatedA -: A - B mp C

kernel3=: 2 : '(((0{u)`:6) ] ;. 0 ]) ((0{v)`:6) (((1{u)`:6) ] ;. 0 ]) ((1{v)`:6) (((2{u)`:6) ] ;. 0 ])'

NB. ---------------------------------------------------------
NB. sdiag                                                 1 2
NB. Shift diagonal of y by values from x, i.e. make x*I+y
NB. from scalar or vector x and matrix y
NB.
NB. Syntax:
NB.   s=. x sdiag y
NB. where
NB.   y - N×N-matrix
NB.   x - numeric scalar of N-vector, shift for y's diagonal
NB.   s - N×N-matrix, y+x*idmat(#y)
NB.   N >= 0

NB. xplusdiagy=: + (< 0 1) & |:                 NB. new diagonal: x + diag(y)
NB. linIOSdiagy=: (>: * i.) @ # @ ]             NB. linear IOS of y's diagonal
NB. sdiag=: (xplusdiagy linIOSdiagy } ]) " 1 2  NB. replace diagonal

sdiag=: ((+ (< 0 1) & |:) ((>: * i.) @ # @ ]) } ]) " 1 2

NB. =========================================================
Note 'testing and timing'

   NB. machine epsilon
   (-: ^: (1 ~: (1: + {:)) ^: _) 1 0.5
1.13687e_13 5.68434e_14
   NB. system setting
   9!:18 ''
5.68434e_14
   NB. epsilon exponent
   2 ^. 1.13687e_13 5.68434e_14
_43 _44
   NB. test
   1 = 1 + 1.13687e_13 5.68434e_14
0 1

   NB. underflow threshold
   (-: ^: (0 ~: {:) ^: _) 1 0.5
4.94066e_324 0
   NB. minimum exponent
   2 ^. 4.94066e_324
_1074
   NB. test
   2 ^ _1074 _1075
4.94066e_324 0
   NB. safe minimum
   _ > 0.5 ^ _1023 _1024
1 0

   NB. overflow threshold
   (+: ^: (_ ~: {.) ^: _) 1 0.5
_ 8.98847e307
   NB. maximum exponent
   2 ^. 8.98847e307
1023
   NB. test
   2 ^ 1023 1024
8.98847e307 _
   NB. safe maximum
   _ > % 0.5 ^ 1023 1024
1 0
)
