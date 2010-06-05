NB. util.ijs
NB. Utilities
NB.
NB. sgn       if y<0 then sgn(y)=-1, else sgn(y)=1
NB. trace     matrix trace
NB. ct        conjugate transpose
NB. lio       integers grid (2{y) steps from (0{y) by (1{y)
NB. ii2cp     make cycle permutation from indices x and y
NB. ht2i      form int list from head y to tail x: y (y+1) ... (x-1)
NB. norm1     magnitude-based 1-norm of matrix or vector
NB. normi     magnitude-based ∞-norm of matrix or vector
NB. norm1t    taxicab-based 1-norm of matrix or vector
NB. normit    taxicab-based ∞-norm of matrix or vector
NB. dbprn     print debug info conditionally
NB. sdiag     add element[s from] x to diagonal of matrix y
NB.
NB. Version: 1.0.0 2009-06-01
NB. Copyright: Igor Zhuravlov, igor at uic.dvgu.ru
NB. License: GNU GPL

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

sgn=: { & 1 1 _1 @ *               NB. if y<0 then sgn(y)=-1, else sgn(y)=1

trace=: +/ @ diag                  NB. matrix trace
ct=: + @ |:                        NB. conjugate transpose

lio=: + ` (* i.)/ " 1              NB. integers grid (2{y) steps from (0{y) by (1{y)
ii2cp=: < @ (, ` (, @ ]) @. =)     NB. make cycle permutation from indices x and y (NB!: if y is empty then output must be a:)
ht2i=: ] + (i. @ -)                NB. form int list from head y to tail x: y (y+1) ... (x-1)

dbprn=: [ ('%-25S %-12g %-12g %-12g %-12g %12d' & printf ^: (VERBOSE"_))

NB. ---------------------------------------------------------
NB. Norms

NB. template adverbs to form norm verbs
mocs=: >./ @ (+/     ) @:          NB. vector: sum of, matrix: max of column sums
mors=: >./ @ (+/ " _1) @:          NB. vector: max of, matrix: max of row sums

NB. magnitude-based norms |y|
norm1=: | mocs                     NB. 1-norm of matrix or vector
normi=: | mors                     NB. ∞-norm of matrix or vector

NB. taxicab-based norms |Re(y)| + |Im(y)|
norm1t=: (+/ " 1 @: | @: +.) mocs  NB. 1-norm of matrix or vector
normit=: (+/ " 1 @: | @: +.) mors  NB. ∞-norm of matrix or vector

NB. ---------------------------------------------------------
NB. sdiag
NB. Shift diagonal of y by values from x, i.e. make x*I+y
NB. from scalar or vector x and matrix y
NB.
NB. Syntax:
NB.   s=. x sdiag y
NB. where
NB.   y - n×n-matrix
NB.   x - numeric scalar of n-vector, shift for y's diagonal
NB.   s - n×n-matrix, equals to (y+x*idmat(#y))
NB.   n >= 0
NB.
NB. TODO:
NB. - implement for non-square matrices

sdiag=: (+ diag) ((>: * i.) @ # @ ]) } ]
