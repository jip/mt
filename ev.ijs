NB. ev.ijs
NB. Eigenvalue decomposition (EVD) of a matrix
NB.
NB. geev  #########
NB. heev  #########
NB. ggev  #########
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-10

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Name: laqr1
NB. Description: 1st column of (H-s1*I)*(H-s2*I)
NB. Syntax: vK=. (s1,s2) laqr1 H
NB. where   H - 2×2- or 3×3-matrix

laqr1=: 4 : 0
  's1 s2'=. x
  H1=. (-&s1) upddiag        y
  cH2=. ((-&s2) updl 0) {."1 y  NB. only 1st column of H is taken
  H1 mp cH2
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. geev
NB.
NB. ######Description:
NB.   Solve:
NB.     A * X = B
NB.
NB. Syntax
NB.     X=. A gesv B
NB. where
NB.   A  - m*n matrix

geev=: [:

NB. =========================================================
NB. Test suite
