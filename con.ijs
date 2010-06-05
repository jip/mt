NB. Condition number
NB.
NB. con     Conj. to make verb estimating the reciprocal of
NB.         the condition number of a matrix in a given norm
NB. xxconx  Calculate reciprocal of the condition number of a
NB.         matrix in a given norm
NB.
NB. Version: 0.6.0 2010-06-05
NB.
NB. Copyright 2010 Igor Zhuravlov
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

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. con
NB.
NB. Description:
NB.   Conj. to make verb estimating the reciprocal of the
NB.   condition number of a matrix in a given norm
NB.
NB. Syntax:
NB.   vapp=. norm con inv
NB. where
NB.   norm - monadic verb to calculate norm of matrix, is
NB.          called as:
NB.            normA=. norm A
NB.   inv  - monadic verb to inverse square matrix, is called
NB.          as:
NB.            invA=. inv A
NB.   vapp - monadic verb to calculate the reciprocal of the
NB.          condition number of a matrix in a given norm, is
NB.          called as:
NB.            rcondA=. vapp A
NB.   A    - n×n-matrix
NB.
NB. TODO:
NB. - implement more practical norm-estimation approach

con=: 2 : '(* & (% @ u)) v'

NB. ---------------------------------------------------------
NB. gecon1
NB. geconi
NB. hecon1
NB. heconi
NB. pocon1
NB. poconi
NB. ptcon1
NB. ptconi
NB. uncon1
NB.
NB. Description:
NB.   Calculate reciprocal of the condition number of a
NB.   matrix in a given norm
NB.
NB. Syntax:
NB.   normG=. geconx G
NB.   normH=. heconx H
NB.   normP=. poconx P
NB.   normT=. ptconx T
NB.   normQ=. uncon1 Q
NB. where
NB.   G     - n×n-matrix of type: general, band,
NB.           tridiagonal, triangular or triangular band
NB.   normG ≥ 0, reciprocal of the condition number of G
NB.   H     - n×n-matrix of type: Hermitian (symmetric)
NB.   normH ≥ 0, reciprocal of the condition number of H
NB.   P     - n×n-matrix of type: Hermitian (symmetric)
NB.           positive definite
NB.   normP ≥ 0, reciprocal of the condition number of P
NB.   T     - n×n-matrix of type: Hermitian (symmetric)
NB.           positive definite tridiagonal
NB.   normT ≥ 0, reciprocal of the condition number of T
NB.   Q     - n×n-matrix, unitary (orthogonal)
NB.   normQ ≥ 0, reciprocal of the condition number of Q in
NB.           1-norm
NB.
NB. Notes:
NB. - extraneous values in triangular, band matrices must be
NB.   zeroed
NB. - gecon1 simulates LAPACK's xGECON('1'), xGBCON('1'),
NB.   xGTCON('1'), xTBCON('1'), xTRCON('1')
NB. - geconi simulates LAPACK's xGECON('I'), xGBCON('I'),
NB.   xGTCON('I'), xTBCON('I'), xTRCON('I')
NB. - hecon1 simulates LAPACK's xHECON('1'), xSYCON('1')
NB. - heconi simulates LAPACK's xHECON('I'), xSYCON('I')
NB. - pocon1 simulates LAPACK's xPOCON('1'), xPBCON('1')
NB. - poconi simulates LAPACK's xPOCON('I'), xPBCON('I')
NB. - ptcon1 simulates LAPACK's xPTCON('1')
NB. - ptconi simulates LAPACK's xPTCON('I')

gecon1=: norm1 con (getrilu1p@getrflu1p)
geconi=: normI con (getrilu1p@getrflu1p)

hecon1=: norm1 con (hetripl@hetrfpl)
heconi=: normi con (hetripl@hetrfpl)

pocon1=: norm1 con (potril@potrfl)
poconi=: normi con (potril@potrfl)

ptcon1=: norm1 con pttril
ptconi=: normi con pttril

uncon1=: norm1 con ct
