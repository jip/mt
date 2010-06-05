NB. orf.ijs
NB. Orthogonal factorizations QR RQ QL LQ RRQR RRRQ RRQL RRLQ
NB.
NB. geqrf  QR factorization of a matrix
NB. geqr2  QR factorization of a matrix (non-blocked version)
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. Difference between consequent cIOSs: cios[i+1]-cios[i]
GEQR2DCIOS=: 4 2 $ 1j_1 1 1j_1 1 0 1 1j_1 1j_1

NB. ---------------------------------------------------------
NB. mkcios0geqr2
NB. Create cIOS for geqr2 at 0-th iteration
NB.
NB. Syntax:
NB.   cios0=. mkcios0geqr2 mn
NB. where
NB.   fs    - 2-vector of integers (m,n), shape of matrix A
NB.   cios0 - 4×2-table cios[0], cIOSs corresponding to
NB.           iteration 0, see geqr2step

mkcios0geqr2=: 3 : 0
  'm1 n1'=. 1 _1 + 'm n'=. y
  4 2 $ (0 0 0 0 , m , 0 0 1) j. (m , 1 , m1 , 1 1 1 , m , n1)
)

geqr2step=: [:

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. geqrf
NB. emulate xGEQR2
NB. RQf=. geqr2 A

geqr2=: 3 : 0
  k=. <./ mn=. $ y
  y=. y , 0    NB. append zero row to A to store τ[0:min(m,n)-1]

  cios0=. mkcios0geqr2 mn  NB. create cios[0]

  NB. link A and cios[0], do iterations, then extract RQf
  0 {:: GEQR2DCIOS (geqr2step ^: (k"_)) (y ; cios0)
)

geqrf=: geqr2  NB. stub for a while

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB.*tgeqrf v test geqrf

tgeqrf=: 3 : 0
)

NB. ---------------------------------------------------------
NB.*testorf v test orthogonal factorizations

testorf=: 3 : 0
)
