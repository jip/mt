NB. ev.ijs
NB. Eigenvalue decomposition (EVD)
NB.
NB. geev  #########
NB. heev  #########
NB. ggev  #########
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

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

NB. =========================================================
NB. Test suite

