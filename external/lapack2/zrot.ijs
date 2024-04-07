require 'math/lapack2'

NB. Description:
NB.   Applies a plane rotation with real cosines to elements
NB.   of vectors
NB.
NB. Syntax:
NB.   'ox oy'=. zrot ix ; iy ; c ; s
NB. where
NB.   ix - n-vector, 1st components of vectors to be rotated
NB.   iy - n-vector, 2nd components of vectors to be rotated
NB.   c  - scalar, real, the cosine of the plane rotation
NB.   s  - scalar, the sine of the plane rotation
NB.   ox - n-vector, 1st components of vectors rotated
NB.   oy - n-vector, 2nd components of vectors rotated
NB.
NB. Notes:
NB. - verb below is loaded into the current locale

zrot=: 3 : 0
  'x y c s'=. y
  assert. x =/&# y
  assert. c , &# s
  2 4 { zrot_jlapack2_ (, # x) ; (, x) ; (, 1) ; (, y) ; (, 1) ; (, c) ; (, s)
)
