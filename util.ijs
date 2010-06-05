NB. util.ijs
NB. Utilities
NB.
NB. sgn       If y<0 then sgn(y)=-1, else sgn(y)=1
NB. trace     matrix trace
NB. ct        Conjugate transpose
NB. lio       Integers grid (2{y) steps from (0{y) by (1{y)
NB. ii2cp     Make cycle permutation from indices x and y
NB. ht2i      Form integers list from head y to tail x
NB. norm1     Magnitude-based 1-norm of vector or matrix
NB. normi     Magnitude-based ∞-norm of vector or matrix
NB. norm1t    Taxicab-based 1-norm of vector or matrix
NB. normit    Taxicab-based ∞-norm of vector or matrix
NB. norms     Square-based (Euclidean/Frobenius) norm of vector or matrix
NB. prn       Formatted conslole output
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

sgn=: 0 & (<: - >)                 NB. if y<0 then sgn(y)=-1, else sgn(y)=1

trace=: +/ @ diag                  NB. matrix trace
ct=: + @ |:                        NB. conjugate transpose

lio=: + ` (* i.)/                  NB. (2{y)-list of integers from (0{y) by (1{y)
ii2cp=: < @ (, ` (, @ ]) @. =)     NB. make cycle permutation from indices x and y (NB!: if y is empty then output must be a:)
ht2i=: ] + (i. @ -)                NB. (x-y)-list of integers from head y to tail (x-1): y (y+1) ... (x-1)
hs2i=: [ + ((] * i. @ *) sgn)~     NB. y-list of integers from head x of size y, models verb's (u;.0) rIOS

condneg=: (- @ ]) ^: (0 > [)       NB. IF x<0 THEN -y ELSE y ENDIF
                                   NB. (* sgn)~ CHECKME!
copysign=: condneg |               NB. IF x<0 THEN -abs(y) ELSE abs(y) ENDIF

prn=: '%-25S %-12g %-12g %-12g %-12g %12d' & printf

NB. ---------------------------------------------------------
NB. Complex IOS (cIOS) / Rectangular IOS (rIOS) / Linear IOS
NB.
NB. Following are equivalents:
NB.    (3 5 _7 ,: 2 _3 4) ] ;. 0 report
NB.    (< 3 4 ; 7 6 5 ; _10 _9 _8 _7) { report
NB.    (cios2ios 3j2 5j_3 _7j4) { report
NB.    (cios2rios 3j2 5j_3 _7j4) ] ;. 0 report

NB. convert cIOS to IOS
NB. TODO: consider M. modifier
cios2ios=: < @ (< @ ({. ^: (1 = #)) @ hs2i/ " 1 @: +.)

NB. model 'from' verb accepting cIOS
cfrom=: ({~ cios2ios)~

NB. model 'amend' adverb accepting cIOS
camend=: 1 : '(cios2ios m) }'

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
