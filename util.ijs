NB. util.ijs
NB. Common pjlap verbs

coclass 'pjlap'

NB. =========================================================
NB. Local verbs

NB. =========================================================
NB. Interface verbs

vnormi=: >./ @: |            NB. inf-norm of vector y
norm1=: >./ @: (+/) @: |     NB. 1-norm of table or vector y
trace=: +/ @ diag            NB. matrix trace

NB. ---------------------------------------------------------

NB. =========================================================
Note 'testing and timing'
)
