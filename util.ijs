NB. util.ijs
NB. Common pjlap verbs
NB.
NB. vnormi  inf-norm of vector y
NB. norm1   1-norm of table or vector y
NB. trace   matrix trace
NB. ut2tr   transform upper triangular matrix to packed form
NB. sdiag   add element[s from] x to diagonal of matrix y

coclass 'pjlap'

NB. =========================================================
NB. Local verbs

NB. =========================================================
NB. Interface verbs

vnormi=: >./ @: |                       NB. inf-norm of vector y
norm1=: >./ @: (+/) @: |                NB. 1-norm of table or vector y
trace=: +/ @ diag                       NB. matrix trace
ut2tr=: ((I. @ , @ (<:/~@i.) @ #) { ,)  NB. transform N-by-N upper triangular matrix to (N*(N+1)/2)-vector, packed form

NB. ---------------------------------------------------------
NB. sdiag                                                 1 2
NB. Shift diagonal of y by values from x, i.e. for scalar or
NB. vector x and matrix y make x*I+y
NB.
NB. Syntax:
NB.   s=. x sdiag y
NB. where
NB.   y - N-by-N table
NB.   x - numeric scalar of N-vector, shift for y's diagonal
NB.   s - N-by-N table, y+x*idmat(#y)
NB.   N >= 0

NB. xplusdiagy=: + (< 0 1) & |:                 NB. new diagonal: x + diag(y)
NB. linIOSdiagy=: (>: * i.) @ # @ ]             NB. linear IOS of y's diagonal
NB. sdiag=: (xplusdiagy linIOSdiagy } ]) " 1 2  NB. replace diagonal

sdiag=: ((+ (< 0 1) & |:) ((>: * i.) @ # @ ]) } ]) " 1 2

NB. =========================================================
Note 'testing and timing'
)
