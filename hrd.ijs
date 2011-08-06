NB. Reduce to Hessenberg form by an unitary similarity
NB. transformation
NB.
NB. gehrdx     Reduce a general matrix to Hessenberg form
NB. gghrdx     Reduce a pair of general and triangular
NB.            matrices to generalized Hessenberg form
NB.
NB. testgehrd  Test gehrdx by general matrix given
NB. testgghrd  Test gghrdx by general matrices given
NB. testhrd    Adv. to make verb to test gxhrdx by matrices
NB.            of generator and shape given
NB.
NB. Version: 0.6.8 2010-08-20
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
  j=. 0
  while. j < HRDNB do.
    b=. (j { y) - (+ (<: j) {"1 V) mp (j {. y)
    b=. b - (((0 (_1) } b) mp (ct V)) mp (ct T)) mp V  NB. matrix-by-vector ops only
    z1=. 1 j } z=. _1 + upd (j , _1) larfg (0 (i. j) } b)
    u=. (* +@{:) z1
    w=. V +@mp (+ - 0 (_1) } u)
    T=. T appendl ((w mp T) , (+ {: z1))
    y=. ((w (i. j) } (0 , u)) mp y) j } y
    H=. H appendl ((j {. b) , (j { z))
    V=. V appendr z1
    j=. >: j
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
  j=. 0
  while. j < HRDNB do.
    b=. (j {"1 y) - (j {."1 y) mp + (<: j) { V
    b=. b - V mp (ct T) mp (ct V) mp 0 (_1) } b  NB. matrix-by-vector ops only
    z1=. 1 j } z=. (j , _1) larfg (0 (i. j) } b)
    u=. (* {:) z1
    w=. (+ - 0 (_1) } u) +@mp V
    T=. T stitcht ((T mp w) , ({: z1))
    y=. (y mp (w (i. j) } (0 , u))) (< a: ; j) } y
    H=. H stitcht ((j {. b) , (j { z))
    V=. V stitchb z1
    j=. >: j
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
NB.   eA  - (n+1)×(n+1)-matrix, being A with appended row of
NB.         trash and column of τs, is already reduced in
NB.         rows 0:h-1
NB.   hs  - 2-vector of integers (h,s) 'head' and 'size',
NB.         defines submatrix A11 to be reduced position in
NB.         matrix A, see gehrdl
NB.   HQf - (n+1)×(n+1)-matrix, combined H and Qf, see gehrdl

gehd2l=: 4 : 0
  A=. ({. x) {. y                             NB. skip ...
  y=. ({. x) }. y                             NB. ...reduced rows
  'j jlimit'=. 1 0 + (+/\) x                  NB. 'j jlimit'=. (h+1),(h+s)
  while. j < jlimit do.                       NB. (s-1)-vector: h+1,h+2,...,h+s-1
    r=. {. y
    z1=. 1 (0) } z=. larfgfc j }. r
    eL=. z1 larflcfr (}. y)                   NB. L := H' * L
    eR=. z1 larfrnfr (j }."1 eL)              NB. R := R * H
    A=. A , ((j {. r) , z)
    y=. (j {."1 eL) ,. eR
    j=. >: j
  end.
  0 (< ((c y) th2lios&<: jlimit);_1) } (A,y)  NB. clear τ[h+s-1:n-1]
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
NB.         τs and column of trash, is already reduced in
NB.         columns 0:h-1
NB.   hs  - 2-vector of integers (h,s) 'head' and 'size',
NB.         defines submatrix A11 to be reduced position in
NB.         matrix A, see gehrdu
NB.   HQf - (n+1)×(n+1)-matrix, combined H and Qf, see gehrdu
NB.
NB. Notes:
NB. - implements LAPACK's xGEHD2 up to storage layout

gehd2u=: 4 : 0
  A=. ({. x) {."1 y                            NB. skip ...
  y=. ({. x) }."1 y                            NB. ...reduced columns
  'j jlimit'=. 1 0 + (+/\) x                   NB. 'j jlimit'=. (h+1),(h+s)
  while. j < jlimit do.                        NB. (s-1)-vector: h+1,h+2,...,h+s-1
    c=. {."1 y
    z1=. 1 (0) } z=. larfgf j }. c
    eR=. z1 larfrnfc (0 1 }. y)                NB. R := R * H
    eL=. z1 larflcfc (j }. eR)                 NB. L := H' * L
    A=. A ,. ((j {. c) , z)
    y=. (j {. eR) , eL
    j=. >: j
  end.
  0 (< _1;((# y) th2lios&<: jlimit)) } (A,.y)  NB. clear τ[h+s-1:n-1]
)

NB. ---------------------------------------------------------
NB. gghrdus
NB.
NB. Description:
NB.   Reduce a pair of matrices (A,B) to generalized upper
NB.   Hessenberg form (H,T) by a unitary (orthogonal)
NB.   similarity transformation by non-blocked algorithm [1]:
NB.     Q^_1 * A * Z = H
NB.     Q^_1 * B * Z = T
NB.   in order to reduce the generalized eigenvalue problem:
NB.     A * x = λ * B * x
NB.   to its standard form:
NB.     H * y = λ * T * y
NB.     y = Z^_1 * x
NB.   and to accumulate rotations to (optionally) form Q and
NB.   Z later.
NB.
NB. Syntax:
NB.   'hs HT QZcs'=. gghrdus hs ; AB
NB. where
NB.   hs   - 2-vector of integers (h,s) 'head' and 'size',
NB.          defines submatrices A11 and B11 to be reduced
NB.          position in matrices A and B, respectively, see
NB.          gehrdu
NB.   AB   -: A ,: B
NB.   A    - n×n-matrix, general
NB.   B    - n×n-matrix, upper triangular
NB.   HT   -: H ,: T
NB.   H    - n×n-matrix, upper Hessenberg in
NB.          H[h:h+s-1,h:h+s-1], and upper triangular outside
NB.   T    - n×n-matrix, upper triangular
NB.   QZcs - ((s-2)*(s-1)/2)×2×2-report, where each
NB.          2×2-matrix is 2-vector of laminated rotations:
NB.            ( c(Q)[i,j]  s(Q)[i,j] )
NB.            ( c(Z)[i,j]  s(Z)[i,j] )
NB.         accumulates rotations to form Q and Z later
NB.
NB. Notes:
NB. - models LAPACK's xGGHRD('N','N')
NB.
NB. TODO:
NB. - re-factor to chain of 2 evokings of the same actor

gghrdus=: 3 : 0
  'hs AB'=. y
  'h sz'=. hs
  QZcs=. 0 2 2 $ 0
  j=. h
  while. j < (h+sz-2) do.                          NB. (s-2)-vector: h,h+1,h+2,...,h+s-3
    i=. h+sz-1
    while. i > (j+1) do.                           NB. (h+s-j-2)-vector: h+s-1,h+s-2,...,j+2
      lios=. i - 1 0
      iospair=. < a: ; lios ; (n th2lios j)
      pair=. iospair { AB
      iosrot=. < 0 ; a: ; 0
      'c s r'=. lartg iosrot { pair
      qcs=. c , s
      pair=. (qcs & (rot " 1 2)) &. (|: " 2) pair  NB. rotate columns
      pair=. (r , 0) iosrot } pair                 NB. place (r,0) explicitely
      AB=. pair iospair } AB
      lios=. i - 0 1
      iospair=. < a: ; a: ; lios                   NB. overtake in B since pair must be unirectangular
      pair=. iospair { AB
      iosrot=. < 1 ; i ; < a:
      'c s r'=. lartg iosrot { pair
      zcs=. c , s
      pair=. zcs rot " 1 2 pair                    NB. rotate rows
      pair=. (r , 0) iosrot } pair                 NB. place (r,0) explicitely
      AB=. pair iospair } AB
      QZcs=. QZcs , qcs ,: zcs
      i=. <: i
    end.
    j=. >: j
  end.
  hs ; AB ; QZcs
)

NB. ---------------------------------------------------------
NB. gghrduqz
NB.
NB. Description:
NB.   Apply accumulated rotations from QZcs to Q1 and/or Z1
NB.   to produce Q and/or Z explicitely
NB.
NB. Syntax:
NB.   'HT Q Z'=. (Q1 ; Z1) gghrduqz hs ; HT ; QZcs
NB. where
NB.   hs   - 2-vector of integers (h,s) 'head' and 'size',
NB.          defines Hessenberg submatrix H11 position in
NB.          matrix H, see gehrdu
NB.   HT   - pass-through argument
NB.   QZcs - ((s-2)*(s-1)/2)×2×2-report, where each
NB.          2×2-matrix is 2-vector of laminated rotations:
NB.            ( c(Q)[i,j]  s(Q)[i,j] )
NB.            ( c(Z)[i,j]  s(Z)[i,j] )
NB.         to form Q and Z later
NB.   Q1  - n×n-matrix or (i.0), the unitary (orthogonal),
NB.         typically from the QR factorization of B
NB.   Q   - either (Q1*dQ) or (i.0)
NB.   Z1  - n×n-matrix or (i.0), the unitary (orthogonal)
NB.   Z   - either (Z1*dZ) or (i.0)
NB.
NB. TODO:
NB. - re-factor to chain of 2 evokings of the same actor

gghrduqz=: 4 : 0
  'Q1 Z1'=. x
  'hs HT QZcs'=. y
  'h s'=. hs
  if. Q1 -.@-: (i. 0) do.
    j=. h
    k=. 0
    QZcs=. (< a: ; 0 ; 0) (+ upd) QZcs     NB. conjugate all s of Qcs
    while. j < (h+s-2) do.
      i=. h+s-1
      while. i > (j+1) do.
        iospair=. < a: ; (i - 1 0)
        qcs=. (< k , 0) { QZcs             NB. extract current rotation qcs[k] from Qcs
        Q1=. iospair ((qcs & rot) upd) Q1  NB. update columns by qcs[k]
        k=. >: k
        i=. <: i
      end.
      j=. >: j
    end.
  end.
  if. Z1 -.@-: (i. 0) do.
    j=. h
    k=. 0
    while. j < (h+s-2) do.
      i=. h+s-1
      while. i > (j+1) do.
        iospair=. < a: ; (i - 0 1)
        zcs=. (< k , 1)  { Zcs             NB. extract current rotation zcs[k] from Zcs
        Z1=. iospair ((zcs & rot) upd) Z1  NB. update columns by zcs[k]
        k=. >: k
        i=. <: i
      end.
      j=. >: j
    end.
  end.
  HT ; Q1 ; Z1
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
NB. Assertions (with appropriate comparison tolerance):
NB.   (idmat n) -: po Q
NB.   H -: A (mp~ mp (ct @ ])) Q
NB. where
NB.   n=. # A
NB.   hs=. 0 , n
NB.   HQf=. hs gehrdl A
NB.   H=. 1 trl 0 _1 }. HQf
NB.   Q=. unghrl HQf

gehrdl=: 4 : 0
  'h s'=. x
  n1=. >: n=. # y
  y=. (2 $ n1) {. y
  Atop=. (h , _) {. y
  Aleft=. (h - (n1 , _1)) {. y
  y=. (h + 0 1) }. y
  I=. 0 >. <. (s - (2+HRDNX-HRDNB)) % HRDNB  NB. how many panels will be reduced
  'i ilimit'=. h + (0,HRDNB) * I             NB. 'i ilimit'=. h,(h+HRDNB*I)
  while. i < ilimit do.                      NB. reduce i-th panel, i = {h,h+HRDNB,...,h+(I-1)*HRDNB} or (HRDNB dhs2lios (h,I))
    'Y V H T'=. lahr2l y                     NB. use (n-i)×(n-i)-matrix A[i:n-1,i+1:n]
    eV0=. 0 ,. (0 (_1) }"1 V)                NB. prepend by zero column, replace τs by zeros
    Aleft=. Aleft - (ct eV0) mp (T mp (eV0 mp Aleft))  NB. update (n-i)×(i+1)-matrix A[i:n-1,0:i]
    y=. (HRDNB }. y) - (ct (HRDNB }."1 eV0)) mp Y      NB. apply reflector from the left
    y=. y - (y mp (ct T mp (0 (_1) }"1 V))) mp V       NB. apply reflector from the right
    V=. ((i. HRDNB) </ (i. (n-i))) } H ,: V  NB. write H into V's lower triangle in-place
    Atop=. Atop , (HRDNB {. Aleft) ,. V
    Aleft=. (HRDNB }. Aleft) ,. (HRDNB {."1 y)
    y=. HRDNB }."1 y
    i=. HRDNB + i
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
NB. Assertions (with appropriate comparison tolerance):
NB.   (idmat n) -: (mp~ ct) Q
NB.   H -: A ((ct @ ]) mp mp) Q
NB. where
NB.   n=. # A
NB.   hs=. 0 , n
NB.   HQf=. hs gehrdu A
NB.   H=. _1 tru _1 0 }. HQf
NB.   Q=. unghru HQf
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
  I=. 0 >. <. (s - (2+HRDNX-HRDNB)) % HRDNB  NB. how many panels will be reduced
  'i ilimit'=. h + (0,HRDNB) * I             NB. 'i ilimit'=. h,(h+HRDNB*I)
  while. i < ilimit do.                      NB. reduce i-th panel, i = {h,h+HRDNB,...,h+(I-1)*HRDNB} or (HRDNB dhs2lios (h,I))
    'Y V H T'=. lahr2u y                     NB. use (n-i)×(n-i)-matrix A[i+1:n,i:n-1]
    eV0=. 0 , (0 (_1) } V)                   NB. prepend by zero row, replace τs by zeros
    Atop=. Atop - ((Atop mp eV0) mp T) mp (ct eV0)  NB. update (i+1)×(n-i)-matrix A[0:i,i:n-1]
    y=. (HRDNB }."1 y) - Y mp (ct (HRDNB }. eV0))   NB. apply reflector from the right
    y=. y - V mp (ct (0 (_1) } V) mp T) mp y        NB. apply reflector from the left
    V=. ((i. (n-i)) >/ (i. HRDNB)) } H ,: V  NB. write H into V's upper triangle in-place
    Aleft=. Aleft ,. (HRDNB {."1 Atop) , V
    Atop=. (HRDNB }."1 Atop) , (HRDNB {. y)
    y=. HRDNB }. y
    i=. HRDNB + i
  end.
  0 _1 }. (x + 1 _1 * HRDNB * I) gehd2u (Aleft ,. Atop , y)
)

NB. ---------------------------------------------------------
NB. gghrdu
NB.
NB. Description:
NB.   Reduce a pair of matrices (A,B) to generalized upper
NB.   Hessenberg form (H,T) by a unitary (orthogonal)
NB.   similarity transformation by non-blocked algorithm [1]:
NB.     Q^_1 * A * Z = H
NB.     Q^_1 * B * Z = T
NB.   in order to reduce the generalized eigenvalue problem:
NB.     A * x = λ * B * x                                 (1)
NB.   to its standard form:
NB.     H * y = λ * T * y
NB.     y = Z^_1 * x
NB.   The unitary (orthogonal) matrices Q and Z are
NB.   determined as products of Givens rotations. They may
NB.   either be formed explicitly, or they may be
NB.   postmultiplied into input matrices Q1 and Z1, so that
NB.     Q1 * A * Z1^_1 = (Q1*Q) * H * (Z1*Z)^_1
NB.     Q1 * B * Z1^_1 = (Q1*Q) * T * (Z1*Z)^_1
NB.   If Q1 is the unitary (orthogonal) matrix from the QR
NB.   factorization of B in the original equation (1), then
NB.   gghrdu reduces the original problem to generalized
NB.   Hessenberg form
NB.
NB. Syntax:
NB.   HT=.                 gghrdu hs ; AB
NB.   'HT Q Z'=. (Q1 ; Z1) gghrdu hs ; AB
NB. where
NB.   hs  - 2-vector of integers (h,s) 'head' and 'size',
NB.         defines submatrices A11 and B11 to be reduced
NB.         position in matrices A and B, respectively, see
NB.         gehrdu
NB.   AB  -: A ,: B
NB.   A   - n×n-matrix, general
NB.   B   - n×n-matrix, upper triangular
NB.   HT  -: H ,: T
NB.   H   - n×n-matrix, upper Hessenberg in
NB.         H[h:h+s-1,h:h+s-1], and upper triangular outside
NB.   T   - n×n-matrix, upper triangular
NB.   Q1  - n×n-matrix or (i.0), the unitary (orthogonal),
NB.         typically from the QR factorization of B
NB.   Q   - either (Q1*dQ) or (i.0)
NB.   Z1  - n×n-matrix or (i.0), the unitary (orthogonal)
NB.   Z   - either (Z1*dZ) or (i.0)
NB.
NB. Notes:
NB. - models LAPACK's xGGHRD
NB.
NB. Application:
NB. - model LAPACK's xGGHRD('N','N'):
NB.     HT=. gghrdu hs;AB
NB.     HT=. 0 {:: (2 # a:) gghrdu (hs;AB)
NB. - model LAPACK's xGGHRD('I','I'):
NB.     'HT Q Z'=. (gghrdu~ (;~ @ idmat @ # @ (1 & {::))) (hs;AB)
NB. - model LAPACK's xGGHRD('V','V'):
NB.     'HT Q Z'=. (Q1;Z1) gghrdu (hs;AB)
NB. - model LAPACK's xGGHRD('V','N'):
NB.     'HT Q'=. 2 {. (Q1;a:) gghrdu (hs;AB)
NB.
NB. References:
NB. [1] G. H. Golub and C. F. Van Loan, Matrix Computations,
NB.     Johns Hopkins University Press, Baltimore, Md, USA,
NB.     3rd edition, 1996, p. 378

gghrdu=: (1 {:: gghrdus) : (gghrduqz gghrdus)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgehrd
NB.
NB. Description:
NB.   Test Hessenberg reduction algorithms:
NB.   - gehrd (math/lapack addon)
NB.   - gehrdx (math/mt addon)
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
  require :: ] '~addons/math/lapack/lapack.ijs'
  need_jlapack_ :: ] 'gehrd'

  rcond=. gecon1 y

  ('2b1100 & gehrd_jlapack_' tmonad        ((];1:;#)`(,&>/)`(rcond"_)`(_."_)`((norm1@(- (((_1 & tru)@:(}:  )) (] mp  (mp  ct)) unghru)))%((FP_EPS*#*norm1)@[)))) y

  ('gehrdl'                  tdyad  ((0,#)`]        `]     `(rcond"_)`(_."_)`((norm1@(- ((( 1 & trl)@:(}:"1)) (] mp~ (mp~ ct)) unghrl)))%((FP_EPS*#*norm1)@[)))) y
  ('gehrdu'                  tdyad  ((0,#)`]        `]     `(rcond"_)`(_."_)`((norm1@(- (((_1 & tru)@:(}:  )) (] mp  (mp  ct)) unghru)))%((FP_EPS*#*norm1)@[)))) y

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
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testhrd_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testhrd_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testhrd_mt_ 150 150

testhrd=: 1 : 'EMPTY_mt_ [ (testgehrd_mt_ @ u) ^: (=/)'
