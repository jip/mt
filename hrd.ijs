NB. Reduce to Hessenberg form by an unitary similarity
NB. transformation
NB.
NB. gehrdx     Reduce a general matrix to Hessenberg form
NB. gghrdx     Reduce a pair of general and triangular
NB.            matrices to generalized Hessenberg form
NB.
NB. testgehrd  Test gehrdx by general matrix given
NB. testhrd    Adv. to make verb to test gxhrdx by matrix of
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

HRDNB=: 32   NB. block size limit
HRDNX=: 128  NB. crossover point, HRDNX ≥ HRDNB

NB. ---------------------------------------------------------
NB. lahr2l
NB.
NB. Description:
NB.   Reduce the first HRDNB rows (panel) of a general matrix
NB.   subeA so that elements behind the 1st supdiagonal are
NB.   zero. The reduction is performed by an unitary
NB.   (orthogonal) similarity transformation: Q * subeA * Q^_1
NB.
NB. Syntax:
NB.   'Y V H T'=. lahr2l subeA
NB. where
NB.   subeA - (n-i)×(n-i)-matrix eA[i:n-1,i+1:n]
NB.   Y     - HRDNB×(n-i)-matrix, Y = T * V * subeA, the
NB.           last column contains trash
NB.   V     - HRDNB×(n-i)-matrix, unit upper triangular,
NB.           the last column contains τ[i:i+HRDNB-1]
NB.   T     - HRDNB×HRDNB-matrix, lower triangular
NB.   H     - HRDNB×HRDNB-matrix, lower triangular
NB.   eA    - n×(n+1)-matrix, being A with stitched trash
NB.           column
NB.   A     - n×n-matrix to reduce
NB.   Q     - (n-i)×(n-i)-matrix, block reflector,
NB.             Q = I - V'*T*V
NB.   VH    - HRDNB×(n-i)-matrix, represents reduced rows
NB.           of subeA, lower triangles of VH and H are
NB.           match, strict upper triangles of VH and V are
NB.           match
NB.   i     - integer from set:
NB.             {h+j*HRDNB, j=0:I-1, I=max(0,1+⌊(s-2-HRDNX)/HRDNB⌋)},
NB.           defines subeA position in the eA
NB.   n     ≥ 0, integer, size of matrix A
NB.
NB. Notes:
NB. - v[i]'s direction is forward, but T is lower triangular

lahr2l=: 3 : 0
  V=. 0 {. y
  T=. H=. 0 0 $ 0
  for_j. i. HRDNB do.
    b=. (j { y) - (+ (<: j) {"1 V) mp (j {. y)
    b=. b - (((0 (_1) } b) mp (ct V)) mp (ct T)) mp V  NB. matrix-by-vector ops only
    z1=. 1 j } z=. _1 + upd1 (j , _1) larfg (0 (i. j) } b)
    u=. (* +@{:) z1
    w=. V +@mp (+ - 0 (_1) } u)
    T=. T (0 append) ((w mp T) , (+ {: z1))
    y=. ((w (i. j) } (0 , u)) mp y) j } y
    H=. H (0 append) ((j {. b) , (j { z))
    V=. V (_1 append) z1
  end.
  (HRDNB {. y) ; V ; H ; T
)

NB. ---------------------------------------------------------
NB. lahr2u
NB.
NB. Description:
NB.   Reduce the first HRDNB columns (panel) of a general
NB.   matrix subeA so that elements below the 1st subdiagonal
NB.   are zero. The reduction is performed by an unitary
NB.   (orthogonal) similarity transformation: Q^_1 * subeA * Q
NB.
NB. Syntax:
NB.   'Y V H T'=. lahr2u subeA
NB. where
NB.   subeA - (n-i)×(n-i)-matrix eA[i+1:n,i:n-1]
NB.   Y     - (n-i)×HRDNB-matrix, Y = subeA * V * T, the
NB.           last row contains trash
NB.   V     - (n-i)×HRDNB-matrix, unit lower triangular,
NB.           the last row contains τ[i:i+HRDNB-1]
NB.   T     - HRDNB×HRDNB-matrix, upper triangular
NB.   H     - HRDNB×HRDNB-matrix, upper triangular
NB.   eA    - (n+1)×n-matrix, being A with appended trash row
NB.   A     - n×n-matrix to reduce
NB.   Q     - (n-i)×(n-i)-matrix, block reflector,
NB.             Q = I - V*T*V'
NB.   VH    - (n-i)×HRDNB-matrix, represents reduced columns
NB.           of subeA, upper triangles of VH and H are
NB.           match, strict lower triangles of VH and V are
NB.           match
NB.   i     - integer from set:
NB.             {h+j*HRDNB, j=0:I-1, I=max(0,1+⌊(s-2-HRDNX)/HRDNB⌋)},
NB.           defines subeA position in the eA
NB.   n     ≥ 0, integer, size of matrix A
NB.
NB. Notes:
NB. - implements LAPACK's xLAHR2 with following differences:
NB.   - upper i×HRDNB part of Y is calculated later in gehrdu
NB.   - V and H are returned separately from each other
NB.   - V is formed explicitely

lahr2u=: 3 : 0
  V=. _ 0 {. y
  T=. H=. 0 0 $ 0
  for_j. i. HRDNB do.
    b=. (j {"1 y) - (j {."1 y) mp + (<: j) { V
    b=. b - V mp (ct T) mp (ct V) mp 0 (_1) } b  NB. matrix-by-vector ops only
    z1=. 1 j } z=. (j , _1) larfg (0 (i. j) } b)
    u=. (* {:) z1
    w=. (+ - 0 (_1) } u) +@mp V
    T=. T (0 stitch) ((T mp w) , ({: z1))
    y=. (y mp (w (i. j) } (0 , u))) (< a: ; j) } y
    H=. H (0 stitch) ((j {. b) , (j { z))
    V=. V (_1 stitch) z1
  end.
  (HRDNB {."1 y) ; V ; H ; T
)

NB. ---------------------------------------------------------
NB. gehd2l
NB.
NB. Description:
NB.   Reduce an augmented general matrix eA to lower
NB.   Hessenberg form H by a unitary (orthogonal) similarity
NB.   transformation by non-blocked algorithm:
NB.     Q * A * Q^_1 = H
NB.
NB. Syntax:
NB.   HQf=. hs gehd2l eA
NB. where
NB.   eA  - (n+1)×(n+1)-matrix, being A with appended trash
NB.         row and column of τs, is already reduced in rows
NB.         [0:h-1]
NB.   hs  - 2-vector of integers (h,s) 'head' and 'size',
NB.         defines submatrix A11 to be reduced position in
NB.         matrix A, see gehrdl
NB.   HQf - (n+1)×(n+1)-matrix, combined H and Qf, see gehrdl

gehd2l=: 4 : 0
  A=. ({. x) {. y                                                       NB. skip ...
  y=. ({. x) }. y                                                       NB. ...reduced rows
  for_j. dhs2lios (x + 1 _1) do.                                        NB. (s-1)-vector: h+1,h+2,...,h+s-1
    r=. {. y
    z1=. 1 (0) } z=. larfgfc j }. r
    eL=. z1 larflcfr (}. y)                                             NB. L := H' * L
    eR=. z1 larfrnfr (j }."1 eL)                                        NB. R := R * H
    A=. A , ((j {. r) , z)
    y=. (j {."1 eL) ,. eR
  end.
  0 ((A (+&#) y) ([ dhs2lios ((_1-[),])) ((c y) - (+/ x)))"_ } (A , y)  NB. clear τ[h+s-1:n-1]
)

NB. ---------------------------------------------------------
NB. gehd2u
NB.
NB. Description:
NB.   Reduce an augmented general matrix eA to upper
NB.   Hessenberg form H by a unitary (orthogonal) similarity
NB.   transformation by non-blocked algorithm:
NB.     Q^_1 * A * Q = H
NB.
NB. Syntax:
NB.   HQf=. hs gehd2u eA
NB. where
NB.   eA  - (n+1)×(n+1)-matrix, being A with appended row of
NB.         τs and trash column, is already reduced in columns
NB.         0:h-1
NB.   hs  - 2-vector of integers (h,s) 'head' and 'size',
NB.         defines submatrix A11 to be reduced position in
NB.         matrix A, see gehrdu
NB.   HQf - (n+1)×(n+1)-matrix, combined H and Qf, see gehrdu
NB.
NB. Notes:
NB. - implements LAPACK's xGEHD2 up to storage layout

gehd2u=: 4 : 0
  A=. ({. x) {."1 y                                  NB. skip ...
  y=. ({. x) }."1 y                                  NB. ...reduced columns
  for_j. dhs2lios (x + 1 _1) do.                     NB. (s-1)-vector: h+1,h+2,...,h+s-1
    c=. {."1 y
    z1=. 1 (0) } z=. larfgf j }. c
    eR=. z1 larfrnfc (0 1 }. y)                      NB. R := R * H
    eL=. z1 larflcfc (j }. eR)                       NB. L := H' * L
    A=. A ,. ((j {. c) , z)
    y=. (j {. eR) , eL
  end.
  0 (dhs2lios (_2,((# y) - (+/ x))))"_ } (A ,. y)  NB. clear τ[h+s-1:n-1]
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gehrdl
NB.
NB. Description:
NB.   Reduce a general matrix A to lower Hessenberg form H
NB.   by a unitary (orthogonal) similarity transformation:
NB.     Q * A * Q^_1 = H
NB.
NB. Syntax:
NB.   HQf=. hs gehrdl A
NB. where
NB.   A   - n×n-matrix to reduce
NB.   hs  - 2-vector of integers (h,s) 'head' and 'size',
NB.         defines submatrix A11 to be reduced position in
NB.         matrix A, see gebalpl and storage layout below
NB.   HQf - n×(n+1)-matrix, combined H and Qf
NB.   H   - n×n-matrix, it has zeros behind 0-th diagonal
NB.         elements [0:h-1] and [h+s:n-1], and zeros behind
NB.         1st supdiagonal
NB.   Qf  - (s-1)×(n-h)-matrix, unit upper triangular (unit
NB.         diagonal not stored), represents Q in factored
NB.         form:
NB.           Q = Π{H(i)',i=h+s-2:h} ,
NB.         where each elementary reflector Q(i) is
NB.         represented as:
NB.           Q(i) = I - v[i] * τ[i] * v[i]' ,
NB.         and values v[i][i-h+1:s-2] and τ[i] are
NB.         stored in (i-h)-th row of Qf:
NB.           Qf[i-h,0:s-1] = (0,...,0,1,v[i][i-h+1],...,v[i][s-2],0,...0,τ[i])
NB.         assuming
NB.           v[i][0:i-h-1] = v[i][s-1:n-h-2] = 0 ,
NB.           v[i][i-h] = 1 ,
NB.         v[i][i-h+1:s-2] is stored in A[i,i+2:h+s-1], τ[i]
NB.         is stored in A[i,n]
NB.   Q   - n×n-matrix, being unit matrix with unitary
NB.         (orthogonal) matrix inserted into elements
NB.         Q[h:h+s-1,h:h+s-1]
NB.
NB. Storage layout:
NB.       (  A00          )
NB.   A = (  A10 A11      )
NB.       (  A20 A21 A22  )
NB. where
NB.   A00     - h×h-matrix, lower triangular
NB.   A11     - s×s-matrix to be reduced
NB.   A22     - (n-(h+s))×(n-(h+s))-matrix, lower triangular
NB.   A10,A21 - matrices to be updated
NB. Example for h=1, s=5, n=7:
NB.   input  A                     output HQf
NB.   (  a                    )    (  a                       )
NB.   (  a  a  a  a  a  a     )    (  a  a  β1 v1 v1 v1    τ1 )
NB.   (  a  a  a  a  a  a     )    (  h  h  h  β2 v2 v2    τ2 )
NB.   (  a  a  a  a  a  a     )    (  h  h  h  h  β3 v3    τ3 )
NB.   (  a  a  a  a  a  a     )    (  h  h  h  h  h  β4    τ4 )
NB.   (  a  a  a  a  a  a     )    (  h  h  h  h  h  h        )
NB.   (  a  a  a  a  a  a  a  )    (  a  a  h  h  h  h  a     )
NB.
NB. If:
NB.   n=. # A
NB.   hs=. 0 , n
NB.   HQf=. hs gehrdl A
NB.   H=. 1 trl 0 _1 }. HQf
NB.   Q=. unghrl HQf
NB. then (with appropriate comparison tolerance)
NB.   (idmat n) -: (mp ct) Q
NB.   H -: A (mp~ mp (ct @ ])) Q

gehrdl=: 4 : 0
  'h s'=. x
  n1=. >: n=. # y
  y=. (2 $ n1) {. y
  Atop=. (h , _) {. y
  Aleft=. (h - (n1 , _1)) {. y
  y=. (h + 0 1) }. y
  I=. 0 >. <. (s - (2+HRDNX-HRDNB)) % HRDNB       NB. how many panels will be reduced
  for_i. HRDNB dhs2lios (h , I) do.               NB. reduce i-th panel, i = {h,h+HRDNB,...,h+(I-1)*HRDNB}
    'Y V H T'=. lahr2l y                          NB. use (n-i)×(n-i)-matrix A[i:n-1,i+1:n]
    eV0=. 0 ,. (0 (_1) }"1 V)                     NB. prepend by zero column, replace τs by zeros
    Aleft=. Aleft - (ct eV0) mp (T mp (eV0 mp Aleft))  NB. update (n-i)×(i+1)-matrix A[i:n-1,0:i]
    y=. (HRDNB }. y) - (ct (HRDNB }."1 eV0)) mp Y      NB. apply reflector from the left
    y=. y - (y mp (ct T mp (0 (_1) }"1 V))) mp V       NB. apply reflector from the right
    V=. ((i. HRDNB) </ (i. (n-i))) } H ,: V       NB. write H into V's lower triangle in-place
    Atop=. Atop , (HRDNB {. Aleft) ,. V
    Aleft=. (HRDNB }. Aleft) ,. (HRDNB {."1 y)
    y=. HRDNB }."1 y
  end.
  _1 0 }. (x + 1 _1 * HRDNB * I) gehd2l (Atop , Aleft ,. y)
)

NB. ---------------------------------------------------------
NB. gehrdu
NB.
NB. Description:
NB.   Reduce a general matrix A to upper Hessenberg form H
NB.   by a unitary (orthogonal) similarity transformation:
NB.     Q^_1 * A * Q = H
NB.
NB. Syntax:
NB.   HQf=. hs gehrdu A
NB. where
NB.   A   - n×n-matrix to reduce
NB.   hs  - 2-vector of integers (h,s) 'head' and 'size',
NB.         defines submatrix A11 to be reduced position in
NB.         matrix A, see gebalpu and storage layout below
NB.   HQf - (n+1)×n-matrix, combined H and Qf
NB.   H   - n×n-matrix, it has zeros under 0-th diagonal
NB.         elements [0:h-1] and [h+s:n-1], and zeros below
NB.         1st subdiagonal
NB.   Qf  - (n-h)×(s-1)-matrix, unit lower triangular (unit
NB.         diagonal not stored), represents Q in factored
NB.         form:
NB.           Q = Π{H(i),i=h:h+s-2} ,
NB.         where each elementary reflector Q(i) is
NB.         represented as:
NB.           Q(i) = I - v[i] * τ[i] * v[i]' ,
NB.         and values v[i][i-h+1:s-2] and τ[i] are
NB.         stored in (i-h)-th column of Qf:
NB.           Qf[0:s-1,i-h] = (0,...,0,1,v[i][i-h+1],...,v[i][s-2],0,...0,τ[i])
NB.         assuming
NB.           v[i][0:i-h-1] = v[i][s-1:n-h-2] = 0 ,
NB.           v[i][i-h] = 1 ,
NB.         v[i][i-h+1:s-2] is stored in A[i+2:h+s-1,i], τ[i]
NB.         is stored in A[n,i]
NB.   Q   - n×n-matrix, being unit matrix with unitary
NB.         (orthogonal) matrix inserted into elements
NB.         Q[h:h+s-1,h:h+s-1]
NB.
NB. Storage layout:
NB.       (  A00 A01 A02  )
NB.   A = (      A11 A12  )
NB.       (          A22  )
NB. where
NB.   A00     - h×h-matrix, upper triangular
NB.   A11     - s×s-matrix to be reduced
NB.   A22     - (n-(h+s))×(n-(h+s))-matrix, upper triangular
NB.   A01,A12 - matrices to be updated
NB. Example for h=1, s=5, n=7:
NB.   input  A                     output HQf
NB.   (  a  a  a  a  a  a  a  )    (  a  a  h  h  h  h  a  )
NB.   (     a  a  a  a  a  a  )    (     a  h  h  h  h  a  )
NB.   (     a  a  a  a  a  a  )    (     β1 h  h  h  h  h  )
NB.   (     a  a  a  a  a  a  )    (     v1 β2 h  h  h  h  )
NB.   (     a  a  a  a  a  a  )    (     v1 v2 β3 h  h  h  )
NB.   (     a  a  a  a  a  a  )    (     v1 v2 v3 β4 h  h  )
NB.   (                    a  )    (                    a  )
NB.                                (     τ1 τ2 τ3 τ4       )
NB.
NB. If:
NB.   n=. # A
NB.   hs=. 0 , n
NB.   HQf=. hs gehrdu A
NB.   H=. _1 tru _1 0 }. HQf
NB.   Q=. unghru HQf
NB. then (with appropriate comparison tolerance)
NB.   (idmat n) -: (mp~ ct) Q
NB.   H -: A ((ct @ ]) mp mp) Q
NB.
NB. Notes:
NB. - implements LAPACK's xGEHRD up to storage layout

gehrdu=: 4 : 0
  'h s'=. x
  n1=. >: n=. # y
  y=. (2 $ n1) {. y
  Aleft=. h {."1 y
  Atop=. (h - (_1 , n1)) {. y
  y=. (h + 1 0) }. y
  I=. 0 >. <. (s - (2+HRDNX-HRDNB)) % HRDNB       NB. how many panels will be reduced
  for_i. HRDNB dhs2lios (h , I) do.               NB. reduce i-th panel, i = {h,h+HRDNB,...,h+(I-1)*HRDNB}
    'Y V H T'=. lahr2u y                          NB. use (n-i)×(n-i)-matrix A[i+1:n,i:n-1]
    eV0=. 0 , (0 (_1) } V)                        NB. prepend by zero row, replace τs by zeros
    Atop=. Atop - ((Atop mp eV0) mp T) mp (ct eV0)  NB. update (i+1)×(n-i)-matrix A[0:i,i:n-1]
    y=. (HRDNB }."1 y) - Y mp (ct (HRDNB }. eV0))   NB. apply reflector from the right
    y=. y - V mp (ct (0 (_1) } V) mp T) mp y        NB. apply reflector from the left
    V=. ((i. (n-i)) >/ (i. HRDNB)) } H ,: V       NB. write H into V's upper triangle in-place
    Aleft=. Aleft ,. (HRDNB {."1 Atop) , V
    Atop=. (HRDNB }."1 Atop) , (HRDNB {. y)
    y=. HRDNB }. y
  end.
  0 _1 }. (x + 1 _1 * HRDNB * I) gehd2u (Aleft ,. Atop , y)
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgehrd
NB.
NB. Description:
NB.   Test Hessenberg reduction algorithms:
NB.   - gehrd (math/lapack addon)
NB.   - gehrdl gehrdu (math/mt addon)
NB.   by general matrix given
NB.
NB. Syntax:
NB.   testgehrd A
NB. where
NB.   A - n×n-matrix
NB.
NB. Formula:
NB. - Q * A * Q^_1 = H : berr := ||A - Q^_1 * H * Q|| / (ε * ||A|| * n)
NB. - Q^_1 * A * Q = H : berr := ||A - Q * H * Q^_1|| / (ε * ||A|| * n)
NB.
NB. Syntax: testgehrd A
NB. where A - general n×n-matrix

testgehrd=: 3 : 0
  require '~addons/math/lapack/lapack.ijs'
  need_jlapack_ 'gehrd'

  rcond=. (norm1 con (getrilu1p@getrflu1p)) y

NB. FIXME!  ('gehrd_jlapack_' tmonad (]`({: , (,   &. > / @ }:))`(rcond"_)`(_."_)`((norm1@(- ((mp~ ungqr) & > /)))%((FP_EPS*#*norm1)@[)))) y

  ('gehrdl'  tdyad ((0,#)`]`]`(rcond"_)`(_."_)`((norm1@(- ((( 1 & trl)@( 0 _1&}.)) (] mp~ (mp~ ct)) unghrl)))%((FP_EPS*#*norm1)@[)))) y
  ('gehrdu'  tdyad ((0,#)`]`]`(rcond"_)`(_."_)`((norm1@(- (((_1 & tru)@(_1  0&}.)) (] mp  (mp  ct)) unghru)))%((FP_EPS*#*norm1)@[)))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testhrd
NB.
NB. Description:
NB.   Adv. to make verb to test gxhrdx by matrix of generator
NB.   and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testhrd
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
NB.     (? @ $ 0:) testhrd_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testhrd_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testhrd_mt_ 150 200

testhrd=: 1 : 'EMPTY_mt_ [ (testgehrd_mt_ @ u) ^: (=/)'
