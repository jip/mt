NB. pow.ijs
NB. Matrix power[s]
NB.
NB. gepow  X-Y-Z matrix of a general matrix y powers 0..(x-1)
NB. dipow  X-Y-Z matrix of a diagonalizable matrix y powers 0..(x-1)
NB. hepow  X-Y-Z matrix of a Hermitian matrix y powers 0..(x-1)

coclass 'pjlap'

NB. =========================================================
NB. Local verbs

NB. =========================================================
NB. Interface verbs

NB. ---------------------------------------------------------
NB. gepow                                               2 0 2
NB. Make X-Y-Z matrix of a general matrix y powers:
NB.   I y y^2 ... y^(x-1)
NB.
NB. Syntax:
NB.   P=. [x] gepow y
NB. where
NB.   y - N-by-N table
NB.   x - integer >= 0, powers count, default is #y
NB.   P - x-by-N-by-N report, matrix y in powers 0..(x-1)
NB.   N >= 0
NB.
NB. References:
NB. [1] http://www.jsoftware.com/jwiki/Essays/Linear_Recurrences

gepow=: (# $: ]) :((4 :0) " 0 2)

  mpi3=. mp/ ^: (3 = (# @ $))           NB. insert mp only between 2-cells of report y (stiff rank)

  pl=. i. >: <. (2 & ^.) x              NB. powers list: 2^i
  pb=. (< @ I. @ (|. " 1) @ #: @ i.) x  NB. pl bits boxed list
  pp=. mp~ ^: pl y                      NB. powers proxy: y^2^i

  pp (mpi3 @ ({~ >)) " 3 0 pb           NB. extract and mp y's powers for each pl atom
)

NB. ---------------------------------------------------------
NB. dipow                                                   2
NB. Matrix powonential of a diagonalizable matrix
NB.
NB. Syntax:
NB.   Apow=. dipow (VR ; Λ)
NB. where
NB.
NB. If:
NB.   
NB. then
NB.   
NB.
NB. Notes:
NB. - 
NB.
NB. TODO:
NB. - 

dipow=: ((0 & {::) ([ mp ((* (+ @ |:))~ ^)) (1 & {::)) : [: " 2

NB. ---------------------------------------------------------
NB. hepow                                                   2
NB. Matrix powonential of a Hermitian (symmetric if real)
NB. matrix
NB.
NB. Syntax:
NB.   Apow=. dipow (VR ; Λ ; VRinv)
NB. where
NB.
NB. If:
NB.   
NB. then
NB.   
NB.
NB. Notes:
NB. - 
NB.
NB. TODO:
NB. - 

hepow=: (0 & {::) mp ((^ @: (1 & {::)) * (2 & {::)) : [: " 2


NB. =========================================================
Note 'pow testing and timing'
)
