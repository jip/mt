NB. Debug
NB.
NB. dbg  Conj. to force verb to show debug info
NB.
NB. Version: 0.13.0 2021-05-21
NB.
NB. Copyright 2010-2021 Igor Zhuravlov
NB.
NB. This file is part of mt
NB.
NB. mt is free software: you can redistribute it and/or
NB. modify it under the terms of the GNU Lesser General
NB. Public License as published by the Free Software
NB. Foundation, either version 3 of the License, or (at your
NB. option) any later version.
NB.
NB. mt is distributed in the hope that it will be useful, but
NB. WITHOUT ANY WARRANTY; without even the implied warranty
NB. of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
NB. See the GNU Lesser General Public License for more
NB. details.
NB.
NB. You should have received a copy of the GNU Lesser General
NB. Public License along with mt. If not, see
NB. <http://www.gnu.org/licenses/>.

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. get shape
dbgshape=: $`($ (; <) $ L: 0)@.(0 < L.)

NB. failure handler
dbgfailed=: 1 : '(dbsig@dberr [ echo@(m ; ''FAILED'' ; coname))@'''''

NB. success handlers
dbgsucceed1=: 1 : '[ echo@(m ; ''SUCCEED'' ; coname@'''' , ''result'' ; dbgshape_mt_    )'
dbgsucceed2=: 1 : '[ echo@(m ; ''SUCCEED'' ; coname@'''' , ''result'' ; dbgshape_mt_ ; <)'

NB. argument[s] handlers
dbgarg1=: 2 : '] [ echo@(n ; ''MONAD''"_ : (''DYAD''"_) ; m ; coname@'''' , (''y'' ; dbgshape_mt_     ) : ((''x'' ; ''y'') ,@,. ,:& dbgshape_mt_     ))'
dbgarg2=: 2 : '] [ echo@(n ; ''MONAD''"_ : (''DYAD''"_) ; m ; coname@'''' , (''y'' ; dbgshape_mt_ ; < ) : ((''x'' ; ''y'') ,@,. ,:&(dbgshape_mt_ ; <)))'

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

dbg1=: 2 : '(n dbgsucceed1_mt_)@u^:(1:`((u b. 0) dbgarg1_mt_ n)) :: (n dbgfailed_mt_)'

dbg2=: 2 : '(n dbgsucceed2_mt_)@u^:(1:`((u b. 0) dbgarg2_mt_ n)) :: (n dbgfailed_mt_)'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. dbg
NB.
NB. Description:
NB.   Conj. to force verb to show debug info with verbosity
NB.   defined by noun DEBUG which is defined in mt.ijs
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

dbg=: 2 : 'u`(u dbg1_mt_ n)`(u dbg2_mt_ n)@.(''DEBUG_mt_''~)'
