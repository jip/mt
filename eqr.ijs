NB. eqr.ijs
NB. Eigenvalues and eigenvectors of structured matrix
NB.
NB. hseqr  Eigenvalues and, optionally, the Schur decomposition
NB.        of a Hessenberg matrix
NB. pteqr  Eigenvalues and, optionally, eigenvectors of a
NB.        symmetric positive definite tridiagonal matrix
NB. steqr  Eigenvalues and, optionally, eigenvectors of a
NB.        symmetric tridiagonal matrix
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Name: laqr1
NB. Description: 1st column of (H-s1*I)*(H-s2*I)
NB. Syntax: vK=. (s1,s2) laqr1 H
NB. where   H - 2×2- or 3×3-matrix
NB. TODO: tacit

laqr1=: 4 : 0
  's1 s2'=. x
  (((-&s1) upddiag ]) (mp (% norm1t)) (((-&s2) updl 0) {."1 ])) y
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. hseqr
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

