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

NB. ---------------------------------------------------------
NB. sdiag                                                 1 2
NB. Shift diagonal of y by values from x, i.e. for scalar or
NB. vector x and matrix y make x*I+y
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
