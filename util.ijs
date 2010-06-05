NB. util.ijs
NB. Common mt staff
NB.
NB. abssq    Re(y)^2 + Im(y)^2
NB. sgn      if y<0 then sgn(y)=-1, else sgn(y)=1
NB. normi    ∞-norm of matrix or vector y
NB. norm1    1-norm of matrix or vector y
NB. trace    matrix trace
NB. ct       conjugate transpose
NB. lio      integers grid (2{y) steps from (0{y) by (1{y)
NB. ii2cp    make cycle permutation from indices x and y
NB. ht2i     form int list from head y to tail x: y (y+1) ... (x-1)
NB. kernel3  form top level computing kernel verbs
NB. sdiag    add element[s from] x to diagonal of matrix y
NB. amend    adverb to emulate (;.)-like IOS for (})

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

abssq=: +/ @: *: @ +.                   NB. Re(y)^2 + Im(y)^2
sgn=: { & 1 1 _1 @ *                    NB. if y<0 then sgn(y)=-1, else sgn(y)=1

normi=: >./ @ (+/ " _1) @: |            NB. ∞-norm of matrix or vector y
norm1=: >./ @ (+/     ) @: |            NB. 1-norm of matrix or vector y

trace=: +/ @ diag                       NB. matrix trace
ct=: + @ |:                             NB. conjugate transpose

lio=: + ` (* i.)/ " 1                   NB. integers grid (2{y) steps from (0{y) by (1{y)
ii2cp=: < @ (, ` (, @ ]) @. =)          NB. make cycle permutation from indices x and y
ht2i=: ] + (i. @ -)                     NB. form int list from head y to tail x: y (y+1) ... (x-1)

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
NB.   updatedA=. iosABC kabc m
NB.   A=. (0 { iosABC) (] ;. 0 ]) m
NB.   B=. (1 { iosABC) (] ;. 0 ]) m
NB.   C=. (2 { iosABC) (] ;. 0 ]) m
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
NB.   y - n×n-matrix
NB.   x - numeric scalar of n-vector, shift for y's diagonal
NB.   s - n×n-matrix, equals to (y+x*idmat(#y))
NB.   n >= 0

NB. xplusdiagy=: + (< 0 1) & |:                 NB. new diagonal: x + diag(y)
NB. linIOSdiagy=: (>: * i.) @ # @ ]             NB. linear IOS of y's diagonal
NB. sdiag=: (xplusdiagy linIOSdiagy } ]) " 1 2  NB. replace diagonal

sdiag=: ((+ (< 0 1) & |:) ((>: * i.) @ # @ ]) } ]) " 1 2

NB. ---------------------------------------------------------
NB. amend
NB. Adverb to emulate (;.)-like IOS for (})
NB.
NB. Syntax:
NB.   amendVerb=. rios amend
NB. where
NB.   rios      - [2×]n-array, rectangular IOS just like for
NB.               "cut" (;.) conj. Vectors are allowed. The
NB.               value (_) as size is allowed
NB.   amendVerb - verb to emulate "amend" (}) adv.
NB.   n         ≥ 0
NB.
NB. TODO:
NB. - amend-in-pace via rios stripping and rank justification
NB.
NB. Applications:
NB. - 

NB. IOS=. IOfrom fs2ios size
NB. if     IOfrom≥0 && size≥0 then IOS are in [IOfrom      ,IOfrom+(size-1)] (grows from IOfrom)
NB. elseif IOfrom<0 && size≥0 then IOS are in [IOfrom-(size-1),IOfrom      ] (grows to IOfrom)
NB. elseif IOfrom≥0 && size<0 then IOS are in [IOfrom-(size+1),IOfrom      ] (falls to IOfrom)
NB. elseif IOfrom<0 && size<0 then IOS are in [IOfrom      ,IOfrom+(size+1)] (falls from IOfrom)
fs2ios=: (4 : '< (x - (sgn x) * 0 <. (sgn x) * (* y) * (<: | y)) + ((* y) * (i. | y))') " 0 0

NB. if m is vector then convert n-vector to 2×n-matrix
NB. transform each (IOfrom,size) pair to boxed IOS list
amend=: 1 : '(< @: (fs2ios/) @: ((0 & ,:) ^: (2 > (# @ $))) m) }'
