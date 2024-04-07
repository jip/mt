require 'math/lapack2'

NB. Description:
NB.   Applies a vector of plane rotations with real cosines
NB.   to elements of vectors
NB.
NB. Syntax:
NB.   'ox oy'=. xlartv ix ; iy ; c ; s
NB. where
NB.   ix - n-vector, 1st components of vectors to be rotated
NB.   iy - n-vector, 2nd components of vectors to be rotated
NB.   c  - n-vector, real, the cosines of the plane rotations
NB.   s  - n-vector, the sines of the plane rotations
NB.   ox - n-vector, 1st components of vectors rotated
NB.   oy - n-vector, 2nd components of vectors rotated
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dlartv=: 3 : 0
  assert. ({. = }.) # S: 0 y
  'x y c s'=. y
  2 4 { dlartv_jlapack2_ (, # x) ; (, x) ; (, 1) ; (, y) ; (, 1) ; (, c) ; (, s) ; (,1)
)

zlartv=: 3 : 0
  assert. ({. = }.) # S: 0 y
  'x y c s'=. y
  2 4 { zlartv_jlapack2_ (, # x) ; (, x) ; (, 1) ; (, y) ; (, 1) ; (, c) ; (, s) ; (,1)
)
