NB. Raise to an integer power
NB.
NB. gepow      Raise a general matrix to integer powers
NB. dipow      Raise a diagonalizable matrix to integer
NB.            powers
NB. hepow      Raise a Hermitian (symmetric) matrix to
NB.            integer powers
NB.
NB. testgepow  Test gepow by general matrix given
NB. testdipow  Test dipow by diagonalizable matrix given
NB. testhepow  Test hepow by Hermitian (symmetric) matrix
NB.            given
NB. testpow    Adv. to make verb to test xxpow by matrix of
NB.            generator and shape given
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
NB. gepow
NB.
NB. Description:
NB.   Raise a general matrix A to in integer powers
NB.
NB. Syntax:
NB.   P=. p gepow A
NB. where
NB.   A - N×N-matrix, a general matrix
NB.   p - non-negative integers array of any shape sh, powers
NB.   P - sh×N×N-array, a matrix A in powers p
NB.   N >= 0
NB.
NB. References:
NB. [1] http://www.jsoftware.com/jwiki/Essays/Linear_Recurrences

gepow=: (4 :0) " _ 2

  mpi3=. mp/ ^: (3 = (# @ $))       NB. apply mp/ only to 3-rank arrays (stiff rank)

  pl=. i. >: <. (2 & ^.) (>./) , x  NB. powers list: 2^i
  pc=. mp~ ^: pl y                  NB. powers cache: A^2^i
  pb=. (< @ I. @ (|. " 1) @ #:) x   NB. pl bits boxed array of shape sh

  pc (mpi3 @ ({~ >)) " 3 0 pb       NB. extract and mp A's powers for each pl atom
)

NB. ---------------------------------------------------------
NB. dipow
NB.
NB. Description:
NB.   Raise a diagonalizable matrix A in integer powers
NB.
NB. Syntax:
NB.   P=. p dipow (RV ; ev ; RVi)
NB. where
NB.   RV  - N×N-matrix, right eigenvectors of input matrix A
NB.   ev  - N-vector, eigenvalues of input matrix A
NB.   RVi - N×N-matrix, inversion of RV
NB.   p   - non-negative integers array of any shape sh, powers
NB.   P   - sh×N×N-array, a matrix A in powers p
NB.   N  >= 0
NB.
NB. If:
NB.   A=. RV mp (diagmat ev) mp RVi
NB.   P=. 2 dipow (RV ; ev ; RVi)
NB. then (with appropriate comparison tolerance)
NB.   P -: A mp A
NB.   P -: RV mp (diagmat (ev ^ 2)) mp RVi

dipow=: ((0 {:: ]) mp"2 ([ ^"1 0~ 1 {:: ]) (*"1 2) 2 {:: ]) " _ 1

NB. ---------------------------------------------------------
NB. hepow
NB.
NB. Description:
NB.   Raise a Hermitian matrix A in integer powers
NB.
NB. Syntax:
NB.   P=. p dipow (RV ; ev)
NB. where
NB.   RV  - N×N-matrix, right eigenvectors of input matrix A
NB.   ev  - N-vector, eigenvalues of input matrix A
NB.   p   - non-negative integers array of any shape sh, powers
NB.   P   - sh×N×N-array, a matrix A in powers p
NB.   N  >= 0
NB.
NB. If:
NB.   A=. RV mp (diagmat ev) mp (+ |: RV)
NB.   P=. 2 dipow (RV ; ev)
NB. then (with appropriate comparison tolerance)
NB.   P -: A mp A
NB.   P -: RV mp (diagmat (ev ^ 2)) mp (+ |: RV)

hepow=: ((0 {:: ]) mp"2 ([ ^"1 0~ 1 {:: ]) *"1 2 +@(|:@(0 {:: ]))) " _ 1

NB. =========================================================
Note 'pow testing and timing'
)
