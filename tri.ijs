NB. Inverse by triangular factorization
NB.
NB. trtrixx    Inverse triangular matrix
NB. getrixxxx  Inverse general matrix
NB. hetripx    Inverse Hermitian (symmetric) matrix
NB. potrix     Inverse Hermitian (symmetric) positive
NB.            definite matrix
NB. pttri      Inverse Hermitian (symmetric) positive
NB.            definite tridiagonal matrix
NB.
NB. testtrtri  Test trtrixx by triangular matrix given
NB. testgetri  Test getrixxxx by general matrix given
NB. testhetri  Test hetripx by Hermitian (symmetric) matrix
NB.            given
NB. testpotri  Test potrix by Hermitian (symmetric) positive
NB.            definite matrix given
NB. testpttri  Test pttri by Hermitian (symmetric) positive
NB.            definite tridiagonal matrix given
NB. testtri    Adv. to make verb to test xxtrixx by matrix of
NB.            generator and shape given
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

TRINB=: 64  NB. block size limit

NB. ---------------------------------------------------------
NB. getrilu1pstep
NB.
NB. Description:
NB.   Single step of non-blocked version of getrilu1p
NB.
NB. Syntax:
NB.   'pfxi1 sfxi1'=. getrilu1pstep (pfxi ; sfxi)
NB. where
NB.   pfxi  - (j+TRINB)×n-matrix after i-th step and before
NB.           (i+1)-th one, contains not yet processed part
NB.   sfxi  - (n-j-TRINB)×n-matrix after i-th step and before
NB.           (i+1)-th one, contains already processed part
NB.   pfxi1 - j×n-matrix, being pfxi after (i+1)-th step
NB.   sfxi1 - (n-j)×n-matrix, being sfxi after (i+1)-th step
NB.   j     = (I-(i+1))*TRINB
NB.   i     = 0:I-1
NB.   I     = ⌊n/TRINB⌋
NB.
NB. Algorithm:
NB.   In:
NB.     pfx(i) := A[0:j+TRINB-1,0:n-1]
NB.     sfx(i) := A[j+TRINB:n-1,0:n-1]
NB.   Out:
NB.     pfx(i+1) := A[0:j-1,0:n-1]
NB.     sfx(i+1) := A[j:n-1,0:n-1]
NB.   1) extract current diagonal block of A(i) from pfx(i):
NB.        U0i := A[j:j+TRINB-1,j:j+TRINB-1] ,
NB.      it contains merged current diagonal blocks of L^_1
NB.      and U1
NB.   2) extract current subdiagonal block row of U1 from
NB.      pfx(i):
NB.        U1i := A[j:j+TRINB-1,j+TRINB:n-1]
NB.   3) extract current block row of L^_1 from pfx(i):
NB.      3.1) extract current block row of A(i):
NB.             Ri := A[j:j+TRINB-1,0:n-1] ,
NB.           it contains merged current row blocks of
NB.           L^_1 and U1
NB.      3.2) zeroize in Ri elements behind the main diagonal
NB.   4) update Ri:
NB.      4.1) do:
NB.             Ri := Ri - U1i * sfx(i)
NB.      4.2) solve:
NB.             Riupd * U0i = Ri
NB.           for Riupd, where U0i is unit upper triangular
NB.   5) assemble output:
NB.      5.1) cut off current block row from pfx(i) to
NB.           produce pfx(i+1)
NB.      5.2) append sfx(i) to Riupd to produce sfx(i+1)
NB.      5.3) link sfx(i+1) to pfx(i+1)

getrilu1pstep=: 3 : 0
  'pfx sfx'=. y
  'j n'=. $ pfx
  j=. j - TRINB
  U0i=. (,.~ j , TRINB) (];.0) pfx
  U1i=. (- TRINB,(# sfx)) {. pfx
  Ri=. (-TRINB) {. pfx
  Ri=. ((i. TRINB) </ ((i. n) - j)) } Ri ,: 0  NB. spec code
  Ri=. Ri - U1i mp sfx
  Ri=. U0i trsmu1x Ri
  ((-TRINB) }. pfx) ; (Ri , sfx)
)

NB. ---------------------------------------------------------
NB. getripl1ustep
NB.
NB. Description:
NB.   Single step of non-blocked version of getripl1u
NB.
NB. Syntax:
NB.   'pfxi1 sfxi1'=. getripl1ustep (pfxi ; sfxi)
NB. where
NB.   pfxi  - n×(j+TRINB)-matrix after i-th step and before
NB.           (i+1)-th one, contains not yet processed part
NB.   sfxi  - n×(n-j-TRINB)-matrix after i-th step and before
NB.           (i+1)-th one, contains already processed part
NB.   pfxi1 - n×j-matrix, being pfxi after (i+1)-th step
NB.   sfxi1 - n×(n-j)-matrix, being sfxi after (i+1)-th step
NB.   j     = (I-(i+1))*TRINB
NB.   i     = 0:I-1
NB.   I     = ⌊n/TRINB⌋
NB.
NB. Algorithm:
NB.   In:
NB.     pfx(i) := A[0:n-1,0:j+TRINB-1]
NB.     sfx(i) := A[0:n-1,j+TRINB:n-1]
NB.   Out:
NB.     pfx(i+1) := A[0:n-1,0:j-1]
NB.     sfx(i+1) := A[0:n-1,j:n-1]
NB.   1) extract current diagonal block of A(i) from pfx(i):
NB.        L0i := A[j:j+TRINB-1,j:j+TRINB-1] ,
NB.      it contains merged current diagonal blocks of U^_1
NB.      and L1
NB.   2) extract current subdiagonal block column of L1 from
NB.      pfx(i):
NB.        L1i := A[j+TRINB:n-1,j:j+TRINB-1]
NB.   3) extract current block column of U^_1 from pfx(i):
NB.      3.1) extract current block column of A(i):
NB.             Ci := A[0:n-1,j:j+TRINB-1] ,
NB.           it contains merged current column blocks of
NB.           U^_1 and L1
NB.      3.2) zeroize in Ci elements under main diagonal
NB.   4) update Ci:
NB.      4.1) do:
NB.             Ci := Ci - sfx(i) * L1i
NB.      4.2) solve:
NB.             Ciupd * L0i = Ci
NB.           for Ciupd, where L0i is unit lower triangular
NB.   5) assemble output:
NB.      5.1) cut off current block column from pfx(i) to
NB.           produce pfx(i+1)
NB.      5.2) stitch sfx(i) to Ciupd to produce sfx(i+1)
NB.      5.3) link sfx(i+1) to pfx(i+1)
NB.
NB. References:
NB. [1] J. DuCroz, N. Higham. Stability of Methods for Matrix
NB.     Inversion, UT-CS-90-119, October, 1990.
NB.     LAPACK Working Note 27.
NB.     http://www.netlib.org/lapack/lawns/downloads/

getripl1ustep=: 3 : 0
  'pfx sfx'=. y
  'n j'=. $ pfx
  j=. j - TRINB
  L0i=. (,.~ j , TRINB) (];.0) pfx
  L1i=. (- (c sfx),TRINB) {. pfx
  Ci=. (-TRINB) {."1 pfx
  Ci=. (((i. n) - j) >/ i. TRINB) } Ci ,: 0  NB. spec code
  Ci=. Ci - sfx mp L1i
  Ci=. L0i trsmxl1 Ci
  ((-TRINB) }."1 pfx) ; (Ci ,. sfx)
)

NB. ---------------------------------------------------------
NB. getripu1lstep
NB.
NB. Description:
NB.   Single step of non-blocked version of getripu1l
NB.
NB. Syntax:
NB.   'pfxi1 sfxi1'=. getripu1lstep (pfxi ; sfxi)
NB. where
NB.   pfxi  - n×j-matrix after i-th step and before
NB.           (i+1)-th one, contains already processed part
NB.   sfxi  - n×(n-j)-matrix after i-th step and before
NB.           (i+1)-th one, contains тще нуе processed part
NB.   pfxi1 - n×(j+TRINB)-matrix, being pfxi after (i+1)-th
NB.           step
NB.   sfxi1 - n×(n-j-TRINB)-matrix, being sfxi after (i+1)-th
NB.           step
NB.   j     = n-(I-i)*TRINB
NB.   i     = 0:I-1
NB.   I     = ⌊n/TRINB⌋
NB.
NB. Algorithm:
NB.   In:
NB.     pfx(i) := A[0:n-1,0:j-1]
NB.     sfx(i) := A[0:n-1,j:n-1]
NB.   Out:
NB.     pfx(i+1) := A[0:n-1,0:j+TRINB-1]
NB.     sfx(i+1) := A[0:n-1,j+TRINB:n-1]
NB.   1) extract current diagonal block of A(i) from sfx(i):
NB.        U0i := A[j:j+TRINB-1,j:j+TRINB-1] ,
NB.      it contains merged current diagonal blocks of L^_1
NB.      and U1
NB.   2) extract current superdiagonal block column of U1
NB.      from sfx(i):
NB.        U1i := A[0:j-1,j:j+TRINB-1]
NB.   3) extract current block column of L^_1 from sfx(i):
NB.      3.1) extract current block column of A(i):
NB.             Ri := A[0:n-1,j:j+TRINB-1] ,
NB.           it contains merged current column blocks of
NB.           L^_1 and U1
NB.      3.2) zeroize in Ci elements above the main diagonal
NB.   4) update Ci:
NB.      4.1) do:
NB.             Ci := Ci - sfx(i) * U1i
NB.      4.2) solve:
NB.             U0i * Ciupd = Ci
NB.           for Ciupd, where U0i is unit upper triangular
NB.   5) assemble output:
NB.      5.1) cut off current block column from sfx(i) to
NB.           produce sfx(i+1)
NB.      5.2) stitch Ciupd to pfx(i) to produce pfx(i+1)
NB.      5.3) link sfx(i+1) to pfx(i+1)

getripu1lstep=: 3 : 0
  'pfx sfx'=. y
  'n j'=. $ pfx
  U0i=. ((j , 0) ,: (2 # TRINB)) (];.0) sfx
  U1i=. (j , TRINB) {. sfx
  Ci=. TRINB {."1 sfx
  Ci=. (((i. n) - j) </ (i. TRINB)) } Ci ,: 0  NB. spec code
  Ci=. Ci - pfx mp U1i
  Ci=. U0i trsmxu1 Ci
  (pfx ,. Ci) ; (TRINB }."1 sfx)
)

NB. ---------------------------------------------------------
NB. getriul1pstep
NB.
NB. Description:
NB.   Single step of non-blocked version of getriul1p
NB.
NB. Syntax:
NB.   'pfxi1 sfxi1'=. getristep (pfxi ; sfxi)
NB. where
NB.   pfxi  - j×n-matrix after i-th step and before
NB.           (i+1)-th one, contains already processed part
NB.   sfxi  - (n-j)×n-matrix after i-th step and before
NB.           (i+1)-th one, contains not yet processed part
NB.   pfxi1 - (j+TRINB)×n-matrix, being pfx(i) after (i+1)-th
NB.           step
NB.   sfxi1 - (n-j-TRINB)×n-matrix, being sfx(i) after
NB.           (i+1)-th step
NB.   j     = n-(I-i)*TRINB
NB.   i     = 0:I-1
NB.   I     = ⌊n/TRINB⌋
NB.
NB. Algorithm:
NB.   In:
NB.     pfx(i) := A[0:j-1,0:n-1]
NB.     sfx(i) := A[j:n-1,0:n-1]
NB.   Out:
NB.     pfx(i+1) := A[0:j+TRINB-1,0:n-1]
NB.     sfx(i+1) := A[j+TRINB:n-1,0:n-1]
NB.   1) extract current diagonal block of A(i) from sfx(i):
NB.        L0i := A[j:j+TRINB-1,j:j+TRINB-1] ,
NB.      it contains merged current diagonal blocks of U^_1
NB.      and L1
NB.   2) extract current subdiagonal block row of L1 from
NB.      sfx(i):
NB.        L1i := A[j:j+TRINB-1,0:j-1]
NB.   3) extract current block row of U^_1 from sfx(i):
NB.      3.1) extract current block row of A(i):
NB.             Ri := A[j:j+TRINB-1,0:n-1] ,
NB.           it contains merged current row blocks of U^_1
NB.           and L1
NB.      3.2) zeroize in Ri elements to the left of the main
NB.           diagonal
NB.   4) update Ri:
NB.      4.1) do:
NB.             Ri := Ri - L1i * pfx(i)
NB.      4.2) solve:
NB.             L0i * Riupd = Ri
NB.           for Riupd, where L0i is unit lower triangular
NB.   5) assemble output:
NB.      5.1) cut off current block row from sfx(i) to
NB.           produce sfx(i+1)
NB.      5.2) append Riupd to pfx(i) to produce pfx(i+1)
NB.      5.3) link sfx(i+1) to pfx(i+1)

getriul1pstep=: 3 : 0
  'pfx sfx'=. y
  'j n'=. $ pfx
  L0i=. ((0 , j) ,: (2 # TRINB)) (];.0) sfx
  L1i=. (TRINB , j) {. sfx
  Ri=. TRINB {. sfx
  Ri=. ((i. TRINB) >/ ((i. n) - j)) } Ri ,: 0  NB. spec code
  Ri=. Ri - L1i mp pfx
  Ri=. L0i trsml1x Ri
  (pfx , Ri) ; (TRINB }. sfx)
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:        Syntax:
NB. trtril       iL=.  trtril  L
NB. trtril1      iL1=. trtril1 L1
NB. trtriu       iU=.  trtriu  U
NB. trtriu1      iU1=. trtriu1 U1
NB.
NB. Description:
NB.   Inverse triangular matrix
NB. where:
NB.   L   - n×n-matrix, lower triangular
NB.   iL  - n×n-matrix, lower triangular, an inversion of L
NB.   L1  - n×n-matrix, unit lower triangular (diagonal is
NB.         not saved)
NB.   iL1 - n×n-matrix, unit lower triangular (diagonal is
NB.         not saved), an inversion of L1
NB.   U   - n×n-matrix, upper triangular
NB.   iU  - n×n-matrix, upper triangular, an inversion of U
NB.   U1  - n×n-matrix, unit upper triangular (diagonal is
NB.         not saved)
NB.   iU1 - n×n-matrix, unit upper triangular (diagonal is
NB.         not saved), an inversion of U1
NB.
NB. Algorithm for trtriu:
NB.   In: A - n×n-matrix
NB.   Out: iU
NB.   1) if 1 < # A
NB.      1.1) then
NB.           1.1.1) form (#A)-vector of zeros:
NB.                    fret=. n $ 0
NB.           1.1.2) find splitting edge:
NB.                    k := ⌈n/2⌉
NB.           1.1.3) mark intervals:
NB.                    fret=. 1 (0,k)} fret
NB.           1.1.4) cut A by fret to block matrix bA:
NB.                    bA = (  A00  A01  )  k
NB.                         (  A10  A11  )  n-k
NB.                            k    n-k
NB.           1.1.5) apply trtriu itself to A00 and A11:
NB.                    iA00=. $: A00
NB.                    iA11=. $: A11
NB.           1.1.6) replace A00 and A11 by iA00 and iA11,
NB.                  respectively, in bA
NB.           1.1.7) inverse A01:
NB.                    iA01=. - iA00 mp A01 mp iA11
NB.           1.1.8) replace A01 by iA01 in bA
NB.           1.1.9) assemble iA from block matrix:
NB.                    iA=. icut bA
NB.     1.2) else return reciprocal:
NB.            iA=. % A
NB.
NB. Notes:
NB. - unit diagonal is not referenced
NB. - models LAPACK's xTRTRI with following differences:
NB.   - blocked not partitioned algorithm is used
NB.   - opposite triangle must be zeroed
NB.
NB. References:
NB. [1] JWiki/Essays/Triangular Matrix Inverse
NB.     Roger Hui
NB.     http://www.jsoftware.com/jwiki/Essays/Triangular%20Matrix%20Inverse
NB.
NB. TODO:
NB. - fret would be sparse

trtril=:  %     `(icut@((2:}~ <@:-@((1 1 {:: ]) mp (1 0 {:: ]) mp 0 0 {:: ]))@(<@$:`<`<;.1~ ;~@((0) 1:`(, >.@-:)`(#~)} #))))@.(1 < #)
trtril1=: (1:"0)`(icut@((2:}~ <@:-@((1 1 {:: ]) mp (1 0 {:: ]) mp 0 0 {:: ]))@(<@$:`<`<;.1~ ;~@((0) 1:`(, >.@-:)`(#~)} #))))@.(1 < #)
trtriu=:  %     `(icut@((1:}~ <@:-@((0 0 {:: ]) mp (0 1 {:: ]) mp 1 1 {:: ]))@(<@$:`<`<;.1~ ;~@((0) 1:`(, >.@-:)`(#~)} #))))@.(1 < #)
trtriu1=: (1:"0)`(icut@((1:}~ <@:-@((0 0 {:: ]) mp (0 1 {:: ]) mp 1 1 {:: ]))@(<@$:`<`<;.1~ ;~@((0) 1:`(, >.@-:)`(#~)} #))))@.(1 < #)

NB. ---------------------------------------------------------
NB. Verb:        Factorization used:    Syntax:
NB. getrilu1p    L * U1 * P = A         iA=. getrilu1p LU1p
NB. getripl1u    P * L1 * U = A         iA=. getripl1u pL1U
NB. getripu1l    P * U1 * L = A         iA=. getripu1l pU1L
NB. getriul1p    U * L1 * P = A         iA=. getriul1p UL1p
NB.
NB. Description:
NB.   Inverse a general matrix A, represented in factored
NB.   form [1]
NB. where
NB.   A    - n×n-matrix
NB.   LU1p - 3-vector of boxes, the output of getrflu1p, the
NB.          matrix A represented in factored form
NB.   pL1U - 3-vector of boxes, the output of getrfpl1u, the
NB.          matrix A represented in factored form
NB.   pU1L - 3-vector of boxes, the output of getrfpu1l, the
NB.          matrix A represented in factored form
NB.   LU1p - 3-vector of boxes, the output of getrflu1p, the
NB.          matrix A represented in factored form
NB.   iA   - n×n-matrix, an inversion of A
NB.
NB. Algorithm for getripl1u:#################
NB.   In: ip L1U
NB.   Out: iA
NB.   1) extract ip and L1U from y
NB.   2) extract U from L1U, inverse it and save into y
NB.   3) replace U by y into L1U
NB.   4) find I, the number of iterations of partitioned
NB.      algorithm:
NB.        I := ⌊n/TRINB⌋
NB.      note: partitioned algorithm will be applied to the
NB.            first I*TRINB columns of L1U
NB.   5) invert A:
NB.      5.1) prepare suffix sfx from the last n%TRINB
NB.           columns of y, which won't be touched by
NB.           partitioned algorithm
NB.           5.1.1) extract last n%TRINB columns from iU
NB.           5.1.2) extract lower triangular
NB.           5.1.3) solve
NB.      5.2) I iterations
NB.      5.3) extract suffix
NB.      5.4) apply ip
NB.
NB. Assertions:
NB.   (%. -: (getrilu1p@getrflu1p)) A
NB.   (%. -: (getripl1u@getrfpl1u)) A
NB.   (%. -: (getripu1l@getrfpu1l)) A
NB.   (%. -: (getriul1p@getrful1p)) A
NB.
NB. Notes:
NB. - getripl1u implements LAPACK's xGETRI
NB.
NB. References:
NB. [1] J. DuCroz, N. Higham. Stability of Methods for Matrix
NB.     Inversion, UT-CS-90-119, October, 1990.
NB.     LAPACK Working Note 27.
NB.     http://www.netlib.org/lapack/lawns/downloads/

getrilu1p=: 3 : 0
  'ip LU1'=. y
  n=. c LU1
  y=. trtril trl LU1
  y=. (</~ i. n) } y ,: LU1  NB. spec code
  I=. <. n % TRINB
  ip (C.^:_1) 1 {:: getrilu1pstep ^: I ((TRINB * I) ({. ; ([ (tru1 trsmu1x trl) }.)) y)
)

getripl1u=: 3 : 0
  'ip L1U'=. y
  n=. # L1U
  y=. trtriu tru L1U
  y=. (>/~ i. n) } y ,: L1U  NB. spec code
  I=. <. n % TRINB
  ip (C.^:_1"1) 1 {:: getripl1ustep ^: I ((TRINB * I) (({."1) ; (-@[ (trl1 trsmxl1 tru) (}."1))) y)
)

NB. num_array3=: [num_value2] bool_IOS_table aamend num_value1
NB.
NB. Notes:
NB. - the following definition:
NB.     amend=: 1 : 0
NB.     :
NB.       y=. m } y ,: x
NB.     )
NB.   doesn't work, see http://www.jsoftware.com/pipermail/programming/2010-March/018944.html

amend=: 1 : 0
  y=. m } y ,: 0
:
  y=. m } y ,: x
)
getripl1ustep2=: 4 : 0
  'n jj'=. $ y
  j=. (n - jj) - TRINB
  'L0i L1i'=. TRINB ({. ; }.) j }. x
  (L0i trsmxl1 ((-j) trupick x) - y mp L1i) ,. y
)
getripl1u2=: 3 : 0
  ip=. 0 {:: y
  y=. (1 {:: y) (>/~ i. # ip) amend trtriu tru 1 {:: y
  ip C.^:_1"1 > (getripl1ustep2&:>/) (({. @ ((0 _ ,. TRINB) & (<;._3))) , (<@((-~/@$) (trl1 trsmxl1 tru) ])@((-@(TRINB&|)@#) {."1 ]))) y
)
getripl1ustep3=: 4 : 0
  'n jj'=. $ y
  j=. (n - jj) - TRINB
  'L0i L1i'=. TRINB ({. ; }.) j }. x
  (L0i trsmxl1 (x (((i. n) - j) <:/ i. TRINB) amend 0) - y mp L1i) ,. y
)
getripl1u3=: 3 : 0
  ip=. 0 {:: y
  y=. (1 {:: y) (>/~ i. # ip) amend trtriu tru 1 {:: y
  ip C.^:_1"1 > (getripl1ustep3&:>/) (({. @ ((0 _ ,. TRINB) & (<;._3))) , (<@((-~/@$) (trl1 trsmxl1 tru) ])@((-@(TRINB&|)@#) {."1 ]))) y
)

getripu1l=: 3 : 0
  'ip U1L'=. y
  n=. c U1L
  y=. trtril trl U1L
  y=. (</~ i. n) } y ,: U1L  NB. spec code
  I=. <. n % TRINB
  ip (C.^:_1"1) 0 {:: getripu1lstep ^: I ((TRINB | n) ((((2 # [) {. ]) trsmxu1 (trl@:({."1))) ; (}."1)) y)
)

getriul1p=: 3 : 0
  'ip UL1'=. y
  n=. # UL1
  y=. trtriu tru UL1
  y=. (>/~ i. n) } y ,: UL1  NB. spec code
  I=. <. n % TRINB
  ip (C.^:_1) 0 {:: getriul1pstep ^: I ((TRINB | n) ((((2 # [) {. ]) trsml1x (tru@{.)) ; }.) y)
)

NB. ---------------------------------------------------------
NB. Verb:    Factorization used:            Syntax:
NB. hetripl  P * L1 * T * L1^H * P^_1 = A   iA=. hetripl pL1T
NB. hetripu  P * U1 * T * U1^H * P^_1 = A   iA=. hetripu pU1T
NB.
NB. Description:
NB.   Inverse Hermitian (symmetric) matrix A, represented in
NB.   factored form
NB. where
NB.   A    - n×n-matrix, Hermitian (symmetric)
NB.   pL1T - 3-vector of boxes, the output of hetrfpl, the
NB.          matrix A represented in factored form
NB.   pU1T - 3-vector of boxes, the output of hetrfpu, the
NB.          matrix A represented in factored form
NB.   iA   - n×n-matrix, an inversion of A
NB.
NB. Assertion:
NB.   iA -: %. A
NB. where
NB.   iA=. hetripl hetrfpl A
NB.   or
NB.   iA=. hetripu hetrfpu A
NB.
NB. Notes:
NB. - is similar to LAPACK's xHETRI, but uses another
NB.   factorization, see hetrfx

hetripl=: (/: @ (0 & {::)) sp (pttri @ (2 & {::)) ((ct @ ]) mp mp) (trtril1 @ (1 & {::))
hetripu=: (/: @ (0 & {::)) sp (pttri @ (2 & {::)) ((ct @ ]) mp mp) (trtriu1 @ (1 & {::))

NB. ---------------------------------------------------------
NB. Verb:     Factorization used:          Syntax:
NB. potril    L * L^H = A                  iA=. potril L
NB. potriu    U * U^H = A                  iA=. potriu U
NB.
NB. Description:
NB.   Inverse Hermitian (symmetric) positive definite matrix
NB.    A, represented in factored form
NB. where
NB.   A  - n×n-matrix, Hermitian (symmetric) positive
NB.        definite
NB.   L  - n×n-matrix, the output of potrfl, lower
NB.        triangular Cholesky factor
NB.   U  - n×n-matrix, the output of potrfu, upper
NB.        triangular Cholesky factor
NB.   iA - n×n-matrix, an inversion of A
NB.
NB. Notes:
NB. - potril models LAPACK's xPOTRI('L'), but uses direct,
NB.   not iterative matrix product

potril=: (mp~ ct) @ trtril
potriu=: (mp~ ct) @ trtriu

NB. ---------------------------------------------------------
NB. Verb:     Factorization used:         Syntax:
NB. pttril    L1 * D * L1^H = A           iA=. [L1D] pttril A
NB. pttriu    U1 * D * U1^H = A           iA=. [U1D] pttriu A
NB.
NB. Description:
NB.   Inverse Hermitian (symmetric) positive definite
NB.   tridiagonal matrix A, represented in factored form [1]
NB. where
NB.   A   - n×n-matrix, Hermitian (symmetric) positive definite
NB.         tridiagonal
NB.   L1D - 2-vector of boxes, the output of pttrfl, the
NB.         matrix A represented in factored form, optional
NB.   U1D - 2-vector of boxes, the output of pttrfu, the
NB.         matrix A represented in factored form, optional
NB.   iA  - n×n-matrix, inversion of A
NB.
NB. Algorithm for dyadic pttril:
NB. Reflexive is used to emulate monadic gerund power:
NB. u^:(v0`v1`v2)y ↔ (v0 y)u^:(v1 y)(v2 y)
NB.
NB. Assertions:
NB.
NB. References:
NB. [1] Moawwad El-Mikkawy, El-Desouky Rahmo. A new recursive
NB.     algorithm for inverting general tridiagonal and
NB.     anti-tridiagonal matrices.
NB.     Applied Mathematics and Computation 204 (2008) 368–372
NB.     http://dx.doi.org/10.1016/j.amc.2008.06.053
NB.
NB. TODO:
NB. - A should be sparse


pttri=: ((4 : 0) ^: (((0:`(+@])`(_1&diag)`,.`((-@,. (1&(|.!.0)))~ }.)`diag fork3)@])`(<:@#@])`((,. @ ((]`- ag) @ (*/\)&.|.) @ ((, (%@{:))&>/) @ ((_1&diag&.>)`(diag&.>) ag) @ pttrfl)@])))~
  io=. -c y
  pi=. io { x
  ((io (>: upd1) (+/"1) ((}. pi) (*"1) (2 {."1 y))) % ({. pi)) ,. y
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testtrtri
NB.
NB. Description:
NB.   Test:
NB.   - 128!:1 (built-in)
NB.   - trtrixx (math/mt addon)
NB.   by triangular matrix given
NB.
NB. Syntax:
NB.   testtrtri A
NB. where
NB.   A - n×n-matrix, lower triangular
NB.
NB. Formula:
NB. - berr := ||I - A * A^_1|| / (ε * ||A|| * ||A^_1|| * n)

testtrtri=: 3 : 0
  L1=. |: U1=. tru1 U=. |: y
  rcondU=.  (norm1 con trtriu ) U
  rcondU1=. (norm1 con trtriu1) U1
  rcondL=.  (norm1 con trtril ) y
  rcondL1=. (norm1 con trtril1) L1

  ('(128!:1)'  tmonad (]`]`(rcondU "_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) U

  ('trtril'    tmonad (]`]`(rcondL "_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) y
  ('trtril1'   tmonad (]`]`(rcondL1"_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) L1
  ('trtriu'    tmonad (]`]`(rcondU "_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) U
  ('trtriu1'   tmonad (]`]`(rcondU1"_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) U1

  EMPTY
)

NB. ---------------------------------------------------------
NB. testgetri
NB.
NB. Description:
NB.   Test getrixxxx by general matrix given
NB.
NB. Syntax:
NB.   testgetri A
NB. where
NB.   A - n×n-matrix
NB.
NB. Formula:
NB. - berr := ||I - A * A^_1|| / (ε * ||A|| * ||A^_1|| * n)

testgetri=: 3 : 0
  rcond=. (norm1 con (getrilu1p@getrflu1p)) y

  ('%.'        tmonad (] `]`(rcond"_)`(_."_)`(           (norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@])))   ))               y

  ('getrilu1p' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; getrflu1p) y
  ('getripl1u' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; getrfpl1u) y
  ('getripl1u2' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; getrfpl1u) y
  ('getripl1u3' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; getrfpl1u) y
  ('getripu1l' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; getrfpu1l) y
  ('getriul1p' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; getrful1p) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhetri
NB.
NB. Description:
NB.   Test hetripx by Hermitian (symmetric) matrix given
NB.
NB. Syntax:
NB.   testhetri A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric)
NB.
NB. Formula:
NB. - berr := ||I - A * A^_1|| / (ε * ||A|| * ||A^_1|| * n)

testhetri=: 3 : 0
  rcond=. _. NB. (norm1 con (hetripl@hetrfpl)) y

  ('hetripl' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; hetrfpl) y
  ('hetripu' tmonad (}.`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; hetrfpu) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpotri
NB.
NB. Description:
NB.   Test hetripx by Hermitian (symmetric) positive definite
NB.   matrix given
NB.
NB. Syntax:
NB.   testpotri A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive definite
NB.
NB. Formula:
NB. - berr := ||I - A * A^_1|| / (ε * ||A|| * ||A^_1|| * n)

testpotri=: 3 : 0
  rcond=. (norm1 con (potril@potrfl)) y

  ('potril' tmonad ((1 & {::)`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; potrfl) y
  ('potriu' tmonad ((1 & {::)`]`(rcond"_)`(_."_)`((0 {:: [) ((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))) ]))) (; potrfu) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpttri
NB.
NB. Description:
NB.   Test pttri by Hermitian (symmetric) positive definite
NB.   tridiagonal matrix given
NB.
NB. Syntax:
NB.   testpttri A
NB. where
NB.   A - n×n-matrix, Hermitian (symmetric) positive
NB.       definite tridiagonal
NB.
NB. Formula:
NB. - berr := ||I - A * A^_1|| / (ε * ||A|| * ||A^_1|| * n)

testpttri=: 3 : 0
EMPTY return.
  rcond=. (norm1 con pttri) y

  ('pttril' tdyad ((pttrfl@])`]`]`(rcond"_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) y
  ('pttriu' tdyad ((pttrfu@])`]`]`(rcond"_)`(_."_)`((norm1@(<: upddiag)@mp)%(FP_EPS*(*&norm1)*(#@]))))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testtri
NB.
NB. Description:
NB.   Adv. to make verb to test triangular inversion
NB.   algorithms by matrix of generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testtri
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     (? @ $ 0:) testtri_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testtri_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testtri_mt_ 150 200

testtri=: 1 : 'EMPTY_mt_ [ ((testpttri_mt_ @ (u ptmat_mt_)) [ (testpotri_mt_ @ (u pomat_mt_)) [ (testhetri_mt_ @ (u hemat_mt_)) [ (testgetri_mt_ @ u) [ (testtrtri_mt_ @ (u trlmat_mt_))) ^: (=/)'
