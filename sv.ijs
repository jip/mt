NB. sv.ijs
NB. Solve linear monomial equation
NB.
NB. gesv  #########
NB. hesv  #########
NB. disv  #########
NB. posv  #########
NB. ddsv  #########
NB. pssv  #########
NB. nssv  #########
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-10

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gesv
NB.
NB. Description:
NB.   Solve:
NB.     A * X = B
NB.
NB. Syntax
NB.     X=. A gesv B
NB. where
NB.   A  - m*n matrix

gesv=: (getrs~ getrf)~

NB. =========================================================
NB. Test suite
