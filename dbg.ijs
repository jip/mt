NB. Debug
NB.
NB. dbg  Conj. to force verb to show debug info
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

gshapes=: $`($(;<)($ L: 0)) @. (0 < L.)  NB. get shapes

NB. ---------------------------------------------------------
NB. dbg1
NB. dbg2
NB.
NB. Description:
NB.   Conj.s to equip verb by debug output
NB.
NB. Syntax:
NB.   vdbg1=. v dbg1 title
NB.   vdbg2=. v dbg2 title
NB. where
NB.   title - any literal to name v
NB.   v     - verb to switch to debug mode
NB.   vdbg1 - being verb v equipped by output of its rank
NB.           and valency, input's and output's shapes
NB.   vdbg2 - the same output as by vdbg1 plus input's and
NB.           output's values

dbg1=: 2 : 0
  smoutput 'dbg' ; (n , ' [MONAD] ' , (": u b. 0)) ; 'y' ; (gshapes_mt_ y)
  o=. u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (gshapes_mt_ o)
  o
:
  smoutput 'dbg' ; 'x' ; (gshapes_mt_ x) ; (n , ' [DYAD] ' , (": u b. 0)) ; 'y' ; (gshapes_mt_ y)
  o=. x u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (gshapes_mt_ o)
  o
)

dbg2=: 2 : 0
  smoutput 'dbg' ; (n , ' [MONAD] ' , (": u b. 0)) ; 'y' ; (gshapes_mt_ y) ; < y
  o=. u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (gshapes_mt_ o) ; < o
  o
:
  smoutput 'dbg' ; 'x' ; (gshapes_mt_ x) ; x ; (n , ' [DYAD] ' , (": u b. 0)) ; 'y' ; (gshapes_mt_ y) ; < y
  o=. x u y
  smoutput 'dbg' ; (n , ' SUCCESS') ; (gshapes_mt_ o) ; < o
  o
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. dbg
NB.
NB. Description:
NB.   Conj. to force verb to show debug info with verbosity
NB.   defined by function DEBUG which is defined in mt.ijs
NB.
NB. Syntax:
NB.   vdbg=. v dbg title
NB. where
NB.   title - any literal to name v
NB.   v     - verb to switch to debug mode
NB.   vdbg  - being verb v equipped by output of debug info
NB.
NB. Application:
NB. - to debug verb '*' in verb (+/ .*) try:
NB.     C=. A (+/ .(* dbg_mt_ '*')) B

dbg=: 2 : 'u`(u dbg1_mt_ n)`(u dbg2_mt_ n) @. DEBUG_mt_'
