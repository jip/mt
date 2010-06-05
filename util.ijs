NB. util.ijs
NB. Utilities
NB.
NB. ii2cp     Make cycle permutation from indices x and y
NB.
NB. sgn       Simplified signum
NB. condneg   Conditional negate
NB. copysign  Copy sign
NB.
NB. prn       Formatted console output
NB.
NB. norm1     Magnitude-based 1-norm of vector or matrix
NB. normi     Magnitude-based ∞-norm of vector or matrix
NB. norm1t    Taxicab-based 1-norm of vector or matrix
NB. normit    Taxicab-based ∞-norm of vector or matrix
NB. norms     Square-based (Euclidean/Frobenius) norm of vector or matrix
NB.
NB. trace     Matrix trace
NB. ct        Conjugate transpose
NB. sdiag     Add element[s from] x to diagonal of matrix y
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. template adverbs to form norm verbs
mocs=: >./ @ (+/     ) @:          NB. vector: sum of, matrix: max of column sums
mors=: >./ @ (+/ " _1) @:          NB. vector: max of, matrix: max of row sums

NB. =========================================================
NB. Interface

ii2cp=: < @ (, ` (, @ ]) @. =)     NB. make cycle permutation from indices x and y (NB!: if y is empty then output must be a:)

sgn=: 0 & (<: - >)                 NB. if y<0 then -1 else 1 endif
condneg=: (* sgn)~                 NB. if x<0 then -y else y endif
copysign=: condneg |               NB. if x<0 then -|y| else |y| endif

prn=: '%-25S %-12g %-12g %-12g %-12g %12d' & printf

dbg=: 2 : 0
  smoutput 'dbg' ; (n , ' [MONAD] ' , (": u b. 0)) ; 'y' ; (($`($ L: 0) @. (0 < L.)) y) ; < y
  o=. u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (($`($ L: 0) @. (0 < L.)) o) ; < o
  o
:
  smoutput 'dbg' ; 'x' ; (($`($ L: 0) @. (0 < L.)) x) ; x ; (n , ' [DYAD] ' , (": u b. 0)) ; 'y' ; (($`($ L: 0) @. (0 < L.)) y) ; < y
  o=. x u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (($`($ L: 0) @. (0 < L.)) o) ; < o
  o
)

NB. ---------------------------------------------------------
NB. Norms

NB. Magnitude-based norms |y|
norm1=: | mocs                     NB. 1-norm of vector or matrix
normi=: | mors                     NB. ∞-norm of vector or matrix

NB. Taxicab-based norms |Re(y)| + |Im(y)|
norm1t=: (+/ " 1 @: | @: +.) mocs  NB. 1-norm of vector or matrix
normit=: (+/ " 1 @: | @: +.) mors  NB. ∞-norm of vector or matrix

NB. Square-based (Euclidean/Frobenius) norm of vector or matrix
NB. for vector input emulates LAPACK's DZNRM2
norms=: (((((+/^:_) &.: *:) @: %) * ]) (>./^:_)) @: | @: +.

NB. ---------------------------------------------------------
NB. Matrix related

trace=: +/ @ diag                  NB. matrix trace
ct=: + @ |:                        NB. conjugate transpose

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
