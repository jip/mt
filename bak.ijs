NB. Restore original eigenvectors by backward transformation
NB. from a balanced matrix or matrix pair
NB.
NB. gebakpx    Undo permutations after gebalpx
NB. gebaksx    Undo scaling after gebalsx
NB. gebakx     Form eigenvectors by backward transformation
NB.            of the matrix balanced by gebalx
NB. ggbakpx    Undo permutations after ggbalpx
NB. ggbaksx    Undo scaling after ggbalsx
NB. ggbakx     Form eigenvectors by backward transformation
NB.            of the pair of matrices balanced by ggbalx
NB.
NB. testgebak  Test gebakxx by general matrix given
NB. testggbak  Test ggbakxx by pair of general matrices given
NB. testbak    Adv. to make verb to test gxbakx by matrix or
NB.            matrix pair of generator and shape given
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
NB. gebakl
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

gebakl=: [:

NB. ---------------------------------------------------------
NB. gebaku
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

gebaku=: [:

NB. =========================================================
NB. Test suite

testbak=: [:
