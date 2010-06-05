NB. util.ijs
NB. Common pjlap verbs

script_z_ '~system/packages/math/mathutil.ijs'  NB. mp
script_z_ '~system/packages/math/matutil.ijs'   NB. diag

coclass 'pjlap'

NB. =========================================================
NB. Local verbs

NB. =========================================================
NB. Interface verbs

vnormi=: >./ @: |            NB. inf-norm of vector y

NB. ---------------------------------------------------------

NB. =========================================================
Note 'testing and timing'
)
