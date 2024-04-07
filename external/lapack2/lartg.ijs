require 'math/lapack2'

NB. Description:
NB.   Generates a plane rotation
NB.
NB. Syntax:
NB.   'c s r'=. xlartg f ; g
NB. where
NB.   f - scalar, the 1st component of vector to be rotated
NB.   g - scalar, the 2nd component of vector to be rotated
NB.   c - scalar, real, the cosine of the rotation
NB.   s - scalar, the sine of the rotation
NB.   r - scalar, the nonzero component of the vector rotated
NB.
NB. Notes:
NB. - verbs below are loaded into the current locale

dlartg=: 3 : 0
  'f g'=. y
  assert. f (,&(0 -: #@$)) g
  ({.L:0) 3 4 5 { dlartg_jlapack2_ (, f) ; (, g) ; (((, 0.0) ; (, 0.0) ; (, 0.0)))
)

zlartg=: 3 : 0
  'f g'=. y
  assert. f (,&(0 -: #@$)) g
  ({.L:0) 3 4 5 { zlartg_jlapack2_ (, f) ; (, g) ; (((, 0.0) ; (, 0j0) ; (, 0j0)))
)
