NB. struct.ijs
NB. Matrix structure handlers
NB.
NB. idmat    make rectangular identity matrix with shifted diagonal
NB. diagmat  make rectangular diagonal matrix from vector
NB. ltri     extract lower triangular (trapezoidal) matrix
NB. utri     extract upper triangular (trapezoidal) matrix
NB. sltri    extract strictly lower triangular (trapezoidal) matrix
NB. sutri    extract strictly upper triangular (trapezoidal) matrix
NB. ltri1    extract unit lower triangular matrix of shape m×min(m,n) <<<<<<<<
NB. utri1    extract unit upper triangular matrix of shape min(m,n)×n <<<<<<<<
NB.
NB. XRef:
NB. - system: mp normalrand rand01
NB. - mt: ct qr

coclass 'mt'

NB. =========================================================
NB. Local definitions

qr=: 128!:0

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. idmat
NB. Make rectangular identity matrix with shifted diagonal
NB. e.g.
NB.    idmat 3
NB. 1 0 0
NB. 0 1 0
NB. 0 0 1
NB.    idmat 3 4
NB. 1 0 0 0
NB. 0 1 0 0
NB. 0 0 1 0
NB.    1 idmat 3 4
NB. 0 1 0 0
NB. 0 0 1 0
NB. 0 0 0 1

idmat=: (0 $: ]) :(= ({. -~/&i. {:)@])

NB. ---------------------------------------------------------
NB. diagmat
NB. Make rectangular diagonal matrix with y on diagonal
NB. x=rows-columns , x=0 is default
NB. e.g.
NB.    diagmat 3 5 7
NB. 3 0 0
NB. 0 5 0
NB. 0 0 7
NB.    1 diagmat 3 5 7
NB. 3 0 0
NB. 0 5 0
NB. 0 0 7
NB. 0 0 0
NB.    _1 diagmat 3 5 7
NB. 3 0 0 0
NB. 0 5 0 0
NB. 0 0 7 0

diagmat=: (0 $: ]) :(((0 (>. , -@<.) [) + #@]) {. (* =@i.@#)@])

NB. ---------------------------------------------------------
NB. ltri
NB. Extract lower triangular (trapezoidal) matrix
NB. e.g.
NB.   ltri 3 5 $ 2
NB. 2 0 0 0 0
NB. 2 2 0 0 0
NB. 2 2 2 0 0
NB.    1 ltri 3 5 $ 2
NB. 2 2 0 0 0
NB. 2 2 2 0 0
NB. 2 2 2 2 0

ltri=: (0 $: ]) : (] * (>: ({. -~/&i. {:)@$@]))

NB. ---------------------------------------------------------
NB. utri
NB. Extract upper triangular (trapezoidal) matrix
NB. e.g.
NB.    utri 3 5 $ 2
NB. 2 2 2 2 2
NB. 0 2 2 2 2
NB. 0 0 2 2 2
NB.    1 utri 3 5 $ 2
NB. 0 2 2 2 2
NB. 0 0 2 2 2
NB. 0 0 0 2 2

utri=: (0 $: ]) : (] * (<: ({. -~/&i. {:)@$@]))

NB. ---------------------------------------------------------
NB. sltri
NB. Extract strictly lower triangular (trapezoidal) matrix
NB. e.g.
NB.    sltri 3 5 $ 2
NB. 0 0 0 0 0
NB. 2 0 0 0 0
NB. 2 2 0 0 0
NB.    1 sltri 3 5 $ 2
NB. 2 0 0 0 0
NB. 2 2 0 0 0
NB. 2 2 2 0 0

sltri=: (0 $: ]) : (] * (> ({. -~/&i. {:)@$@]))

NB. ---------------------------------------------------------
NB. sutri
NB. Extract strictly upper triangular (trapezoidal) matrix
NB. e.g.
NB.    sutri 3 5 $ 2
NB. 0 2 2 2 2
NB. 0 0 2 2 2
NB. 0 0 0 2 2
NB.    1 sutri 3 5 $ 2
NB. 0 0 2 2 2
NB. 0 0 0 2 2
NB. 0 0 0 0 2

sutri=: (0 $: ]) : (] * (< ({. -~/&i. {:)@$@]))

NB. ---------------------------------------------------------
NB. extract unit lower triangular matrix of shape m×min(m,n)

ltri1=: sltri + ((idmat @ $) : (idmat $))

NB. ---------------------------------------------------------
NB. extract unit upper triangular matrix of shape min(m,n)×n
NB. ambivalent

utri1=: sutri + ((idmat @ $) : (idmat $))

NB. ---------------------------------------------------------
NB. TODO:
NB.   mn=. <./ 'm n'=. $ y
NB.   im=. i. m
NB.   in=. i. n
NB.   I=. im =/ in                           NB. identity rectangular matrix
NB.   U=. im <:/ in                          NB. upper triangle rectangular matrix
NB.   L=. im >/ in                           NB. strict lower triangle rectangular matrix
NB.   LU2L=. (n & ({. " 1)) @: (L & *)       NB. L is m×min(m,n) lower traingular matrix
NB.   LU2L=. (mn & ({. " 1)) @: (I + L & *)  NB. L is m×min(m,n) unit lower traingular matrix
NB.   LU2U=. (n & {.) @: (U & *)             NB. U is min(m,n)×n upper traingular matrix
