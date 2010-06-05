NB. brd.ijs
NB. Reduce a general matrix to bidiagonal form by a unitary
NB. transformation
NB.
NB. gebrdl  Reduce a general matrix to lower bidiagonal form by
NB.         a unitary transformation
NB. gebrdu  Reduce a general matrix to lower bidiagonal form by
NB.         a unitary transformation
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gebrdl
NB.
NB. Description:
NB.   
NB.
NB. Syntax:
NB.   
NB. where
NB.
NB. If:
NB.   
NB. then (with appropriate comparison tolerance)
NB.   
NB.
NB. References:
NB. [1] 
NB.
NB. Notes:
NB. - 
NB.
NB. TODO:
NB. - 

gebrdl=: [:

NB. ---------------------------------------------------------
NB. gebrdu
NB.
NB. Description:
NB.   
NB.
NB. Syntax:
NB.   
NB. where
NB.
NB. If:
NB.   
NB. then (with appropriate comparison tolerance)
NB.   
NB.
NB. References:
NB. [1] 
NB.
NB. Notes:
NB. - 
NB.
NB. TODO:
NB. - 

gebrdu=: [:


NB. =========================================================
NB. Test suite

tbrd=: [:

testbrd=: [:
