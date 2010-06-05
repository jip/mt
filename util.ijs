NB. Utilities
NB.
NB. sgn       Simplified signum
NB. condneg   Conditional negate
NB. copysign  Copy sign
NB. sorim     Sum of real and imaginary parts' modules
NB. soris     Sum of real and imaginary parts' squares
NB. fmtlog    Format log string
NB. ag        Adv. to apply successive verbs from gerund to
NB.           successive elements of list
NB. ms        Minimum in sum of vectors
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
NB. Miscellaneous

condneg=: -@]^:(0>[)                                    NB. if x<0 then -y else y endif
copysign=: -@]^:((=-)&*)                                NB. if x<0 then -|y| else |y| endif
sorim=: +/"1 @: |  @: +.                                NB. sum of real and imaginary parts' modules, |Re(y)| + |Im(y)|
soris=: +/"1 @: *: @: +.                                NB. sum of real and imaginary parts' squares, Re(y)^2 + Im(y)^2
fmtlog=: '%-25S %-16g %-16g %-16g %-16g %16d' vsprintf  NB. log string format

NB. ---------------------------------------------------------
NB. ag
NB.
NB. Description
NB.   Adv. to apply successive verbs from gerund to
NB.   successive elements of list
NB.
NB. Syntax:
NB.   vapp=: g ag
NB. where
NB.   g    - gerund u0`u1`... ; each monad ui is called as:
NB.            eiupd=. ui ei
NB.   vapp - monad to apply successive ui to successive ei;
NB.          is called as:
NB.             Eupd=. vapp E
NB.   E    = rank-1 array (e0,e1,...)
NB.   Eupd = rank-1 array (e0upd,e1upd,...)
NB.
NB. References:
NB. [0] [Jforum] gerund apply
NB.     Henry Rich, Sat Oct 22 06:37:12 HKT 2005
NB.     http://www.jsoftware.com/pipermail/general/2005-October/025450.html
NB. [1] [Jforum] gerund apply
NB.     Jose Mario Quintana, Sat Oct 22 10:08:38 HKT 2005
NB.     http://www.jsoftware.com/pipermail/general/2005-October/025459.html

ag=: /. (,/@)

NB. ---------------------------------------------------------
NB. ms
NB.
NB. Description: Minimum in [sum of] vector[s]
NB. Syntax:      k=. [(delta0,delta1,...)] ms (value0,value1,...)
NB. where        default deltai is 0
NB. Formula:     k = min(delta0+value0,delta1+value1,...)
NB. Notes:       is memo, since repetitive calls are expected

ms=: <./@:(] :+)M.
