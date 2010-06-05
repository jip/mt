NB. mq.ijs
NB. Multiply by Q from LQ QL QR RQ HRD output
NB.
NB. unmlqxx  Multiply a matrix by a matrix with orthonormal
NB.          rows as returned by gelqf
NB. unmqlxx  Multiply a matrix by a matrix with orthonormal
NB.          columns as returned by geqlf
NB. unmqrxx  Multiply a matrix by a matrix with orthonormal
NB.          columns as returned by geqrf
NB. unmrqxx  Multiply a matrix by a matrix with orthonormal
NB.          rows as returned by gerqf
NB. unmhrxx  Multiply a matrix by an unitary (orthogonal)
NB.          matrix with orthonormal columns as returned by
NB.          gehrd
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

UNMQBS=: 32  NB. block size limit

unml2ln=: ((1 1 }. [) (({. @ ]) , ($: }.)) ((larfrcfr {.)~)) ^: (0 < # @ ])    NB. FIXME larfrcfr swap args
