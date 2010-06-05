NB. util.ijs
NB. Utilities
NB.
NB. ii2cp     Make cycle permutation from indices x and y
NB.
NB. sgn       Simplified signum
NB. condneg   Conditional negate
NB. copysign  Copy sign
NB.
NB. hds2ios   IOS from head, delta and size
NB. ht2ios    IOS from head and tail
NB. hs2ios    IOS from head and size
NB. cios2ios  Convert cIOS to IOS
NB.
NB. step      Template adv. to make verbs of single iteration
NB.
NB. prn       Formatted console output
NB. dbg       Conj. to show verb's input and output
NB.
NB. norm1     Magnitude-based 1-norm of vector or matrix
NB. normi     Magnitude-based ∞-norm of vector or matrix
NB. norm1t    Taxicab-based 1-norm of vector or matrix
NB. normit    Taxicab-based ∞-norm of vector or matrix
NB. norms     Square-based (Euclidean/Frobenius) norm of vector or matrix
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
                                   NB. CHECKME: ii2cp=: < @ , @ (, ^: ~:)
                                   NB. CHECKME: ii2cp=: (, ^: (-. @ -:)) & ,

sgn=: 0 & (<: - >)                 NB. if y<0 then -1 else 1 endif
condneg=: (* sgn)~                 NB. if x<0 then -y else y endif
copysign=: condneg |               NB. if x<0 then -|y| else |y| endif

NB. ---------------------------------------------------------
NB. IOS generators

NB. y[2]-vector of integers from head y[0] by delta y[1]
hds2ios=: + ` (* i.)/

NB. (x-y)-vector of integers from head y to tail (x-1):
NB. y (y+1) ... (x-1)
ht2ios=: ] + (i. @ -)

NB. y-vector of integers from head x; models rIOS in (u;.0)
hs2ios=: [ + ((] * i. @ *) sgn)~

NB. ---------------------------------------------------------
NB. IOS converters

NB. Convert cIOS to IOS
NB. Note: IOS with length less than array's rank indexes
NB.       the slice
cios2ios=: < " 1 @ (< @ hs2ios/ " 1 @: +.)

NB. ---------------------------------------------------------
NB. step
NB. Template adv. to make verbs of single iteration
NB.
NB. Syntax:
NB.   vstep=: vchange step
NB. where
NB.   vchange - verb to change matrix, is called as:
NB.               Ai1=. ciosi vchange Ai
NB.   vstep   - verb to do single iteration, is called as:
NB.               'Ai1 ciosi1'=. dcios vstep (Ai ; ciosi)
NB.   Ai      - matrix A(i) to update before i-th
NB.             iteration
NB.   ciosi   - matrix cios(i) of cIOSs for i-th iteration
NB.   dcios   - difference between cIOSs at consequent
NB.             iterations: cios(i+1)-cios(i)
NB.   Ai1     - matrix A(i+1) after i-th iteration
NB.   ciosi1  - matrix cios(i+1) of cIOSs for (i+1)-th
NB.             iteration
NB.   i       ≥ 0
NB.
NB. Algorithm:
NB.   1) Ai1=. ciosi vchange Ai
NB.   2) Tmp=. ((< Ai1) 0} Input
NB.   3) ciosi1=. dcios + ciosi
NB.   4) Output=. (< ciosi1) 1} Tmp
NB.   where
NB.     Input -: (Ai ; ciosi)
NB.     Tmp -: (Ai1 ; ciosi)
NB.     Output -: (Ai1 ; ciosi1)

step=: 1 : '(< @ (+ (1&{::))) 1} (((1 {:: ]) (< @ u) (0 {:: ])) 0} ])'

NB. ---------------------------------------------------------
NB. prn
NB. Formatted console output

prn=: '%-25S %-12g %-12g %-12g %-12g %12d' & printf

NB. ---------------------------------------------------------
NB. dbg
NB. Conj. to show verb's input and output

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
