NB. Extended forks
NB.
NB. fork3  Generalized fork with 3-depth execution graph
NB. fork4  Generalized fork with 4-depth execution graph
NB. fork5  Generalized fork with 5-depth execution graph
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
NB. fork2
NB.
NB. Description:
NB.   2-fork, the traditional fork with 2-depth execution
NB.   graph
NB.
NB. Syntax:
NB.   vapp=. a0`b`a1 fork2
NB.   vapp=. a0`b`a1`:6
NB.   vapp=. a0 b a1
NB. where
NB.   ax   - ambivalent verbs to define input nodes of
NB.          execution graph
NB.   b    - dyad to define output node of execution graph
NB.   vapp - 2-fork, is evoked as:
NB.            out=.   vapp y
NB.            out=. x vapp y
NB.
NB. Execution graph:
NB.      b
NB.     / \
NB.    a0  a1
NB.
NB. Notes:
NB. - included just to eludicate generalized forks idea

NB. fork2=: 1 : 0
NB.   '`a0 b a1'=. m
NB.   o0=.    a0 y
NB.   o1=.    a1 y
NB.   o0=. o0 b  o1
NB. :
NB.   '`a0 b a1'=. m
NB.   o0=. x  a0 y
NB.   o1=. x  a1 y
NB.   o0=. o0 b  o1
NB. )

NB. ---------------------------------------------------------
NB. fork3
NB.
NB. Description:
NB.   3-fork, generalized fork with 3-depth execution graph
NB.
NB. Syntax:
NB.   vapp=. a0`b0`a1`c`b1`a2 fork3
NB. where
NB.   ax   - ambivalent verbs to define input nodes of
NB.          execution graph
NB.   bx   - dyads to define intermediate nodes of execution
NB.          graph
NB.   c    - dyad to define output node of execution graph
NB.   vapp - 3-fork, is evoked as:
NB.            out=.   vapp y
NB.            out=. x vapp y
NB.
NB. Execution graph:
NB.       c
NB.      / \
NB.     b0  b1
NB.    / \ / \
NB.   a0  a1  a2
NB.
NB. Notes:
NB. - local nouns are re-used to reduce memory consumption
NB. - another 2-fork compatible traverse order is possible:
NB.     a0 b0 c b1 a2 a1

fork3=: 1 : 0
  '`a0 b0 a1 c b1 a2'=. m
  o0=.    a0 y
  o1=.    a1 y
  o2=.    a2 y
  o0=. o0 b0 o1
  o1=. o1 b1 o2
  o2=. 0
  o0=. o0 c  o1
:
  '`a0 b0 a1 c b1 a2'=. m
  o0=. x  a0 y
  o1=. x  a1 y
  o2=. x  a2 y
  o0=. o0 b0 o1
  o1=. o1 b1 o2
  o2=. 0
  o0=. o0 c  o1
)

NB. ---------------------------------------------------------
NB. fork4
NB.
NB. Description:
NB.   4-fork, generalized fork with 4-depth execution graph
NB.
NB. Syntax:
NB.   vapp=. a0`b0`a1`c0`b1`a2`d`c1`b2`a3 fork4
NB. where
NB.   ax    - ambivalent verbs to define input nodes of
NB.           execution graph
NB.   bx cx - dyads to define intermediate nodes of execution
NB.           graph
NB.   d     - dyad to define output node of execution graph
NB.   vapp  - 4-fork, is evoked as:
NB.            out=.   vapp y
NB.            out=. x vapp y
NB.
NB. Execution graph:
NB.         d
NB.        / \
NB.       c0  c1
NB.      / \ / \
NB.     b0  b1  b2
NB.    / \ / \ / \
NB.   a0  a1  a2  a3
NB.
NB. Notes:
NB. - local nouns are re-used to reduce memory consumption
NB. - another 2-fork compatible traverse order is possible:
NB.     a0 b0 c0 d c1 b2 a3 a1 b1 a2

fork4=: 1 : 0
  '`a0 b0 a1 c0 b1 a2 d c1 b2 a3'=. m
  o0=.    a0 y
  o1=.    a1 y
  o2=.    a2 y
  o3=.    a3 y
  o0=. o0 b0 o1
  o1=. o1 b1 o2
  o2=. o2 b2 o3
  o3=. 0
  o0=. o0 c0 o1
  o1=. o1 c1 o2
  o2=. 0
  o0=. o0 d  o1
:
  '`a0 b0 a1 c0 b1 a2 d c1 b2 a3'=. m
  o0=. x  a0 y
  o1=. x  a1 y
  o2=. x  a2 y
  o3=. x  a3 y
  o0=. o0 b0 o1
  o1=. o1 b1 o2
  o2=. o2 b2 o3
  o3=. 0
  o0=. o0 c0 o1
  o1=. o1 c1 o2
  o2=. 0
  o0=. o0 d  o1
)

NB. ---------------------------------------------------------
NB. fork5
NB.
NB. Description:
NB.   5-fork, generalized fork with 5-depth execution graph
NB.
NB. Syntax:
NB.   vapp=. a0`b0`a1`c0`b1`a2`d0`c1`b2`a3`e`d1`c2`b3`a4 fork5
NB. where
NB.   axx      - ambivalent verbs to define input nodes of
NB.              execution graph
NB.   bx cx dx - dyads to define intermediate nodes of
NB.              execution graph
NB.   e        - dyad to define output node of execution
NB.              graph
NB.   vapp     - 5-fork, is evoked as:
NB.                out=.   vapp y
NB.                out=. x vapp y
NB.
NB. Execution graph:
NB.           e
NB.          / \
NB.         d0  d1
NB.        / \ / \
NB.       c0  c1  c2
NB.      / \ / \ / \
NB.     b0  b1  b2  b3
NB.    / \ / \ / \ / \
NB.   a0  a1  a2  a3  a4
NB.
NB. Notes:
NB. - local nouns are re-used to reduce memory consumption
NB. - another 2-fork compatible traverse order is possible:
NB.     a0 b0 c0 d0 e d1 c2 b3 a4 a1 b1 c1 b2 a3 a2

fork5=: 1 : 0
  '`a0 b0 a1 c0 b1 a2 d0 c1 b2 a3 e d1 c2 b3 a4'=. m
  o0=.    a0 y
  o1=.    a1 y
  o2=.    a2 y
  o3=.    a3 y
  o4=.    a4 y
  o0=. o0 b0 o1
  o1=. o1 b1 o2
  o2=. o2 b2 o3
  o3=. o3 b3 o4
  o4=. 0
  o0=. o0 c0 o1
  o1=. o1 c1 o2
  o2=. o2 c2 o3
  o3=. 0
  o0=. o0 d0 o1
  o1=. o1 d1 o2
  o2=. 0
  o0=. o0 e  o1
:
  '`a0 b0 a1 c0 b1 a2 d0 c1 b2 a3 e d1 c2 b3 a4'=. m
  o0=. x  a0 y
  o1=. x  a1 y
  o2=. x  a2 y
  o3=. x  a3 y
  o4=. x  a4 y
  o0=. o0 b0 o1
  o1=. o1 b1 o2
  o2=. o2 b2 o3
  o3=. o3 b3 o4
  o4=. 0
  o0=. o0 c0 o1
  o1=. o1 c1 o2
  o2=. o2 c2 o3
  o3=. 0
  o0=. o0 d0 o1
  o1=. o1 d1 o2
  o2=. 0
  o0=. o0 e  o1
)
