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
NB. Version: 0.6.8 2010-11-30
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
NB.           defines subeA position in eA
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
NB.           defines subeA position in eA
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
NB. gghrdl
NB.
NB. Description:
NB.   Reduce a pair of matrices (A,B) to generalized lower
NB.   Hessenberg form (H,T) by a unitary (orthogonal)
NB.   similarity transformation by non-blocked algorithm:
NB.     Q * A * Z^_1 = H
NB.     Q * B * Z^_1 = T
NB.   in order to reduce the generalized eigenvalue problem:
NB.     x * A = λ * x * B
NB.   to its standard form:
NB.     y * H = λ * y * T
NB.     y = x * Q^_1
NB.   and accumulate rotations to form Q and Z later
NB.
NB. Syntax:
NB.   'HT dQ dZ'=. hs gghrdl A ,: B
NB. where
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines submatrices A11 and B11 to be reduced
NB.           position in matrices A and B, respectively, see
NB.           gehrdl
NB.   A     - n×n-matrix, general
NB.   B     - n×n-matrix, lower triangular
NB.   HT    -: H ,: T
NB.   H     - n×n-matrix, lower Hessenberg inside the
NB.           submatrix H[h:h+s-1,h:h+s-1], and lower
NB.           triangular outside
NB.   T     - n×n-matrix, lower triangular
NB.   dQ,dZ - any×4-matrix, accumulates rotations to form Q
NB.           and Z later, see rotsclx; dQ and dZ may have
NB.           the same shapes

gghrdl=: 4 : 0
  'h s'=. x
  t=. h+s-1
  n=. c y
  dQ=. dZ=. 0 4 $ 0
  i=. h
  liosr1a=. n th2lios i                        NB. (n-h)-vector h:n-1
  liosc2a=. i. h+s                             NB. (h+s)-vector 0:h+s-1
  while. i < <: t do.                          NB. (s-2)-vector: h:h+s-3
    j=. t
    liosr1b=. n th2lios <: j                   NB. (n-h-s+2)-vector h+s-2:n-1
    liosc2b=. i. >: j                          NB. (j+1)-vector 0:h+s-1
    while. j > >: i do.                        NB. (h+s-i-2)-vector (desc) h+s-1:i+2
      lios=. j - 1 0
      NB. step 1: rotate columns lios to kill A[i,j]
      'y cs'=. rot rotga y ; (< 0 ; liosr1a ; lios) ; 0
      y=. (< 1 ; liosr1b ; lios) (cs & rot) upd y
      dZ=. dZ , cs , lios
      lios=. j - 0 1
      NB. step 2: rotate rows lios to kill B[j-1,j]
      'y cs'=. (rot &. |:) rotga y ; (< 1 ; lios ; liosc2b) ; < < a: ; _1
      y=. (< 0 ; lios ; liosc2a) (cs & (rot &. |:)) upd y
      dQ=. dQ , cs , lios
      NB. step 3: update IOS
      liosr1b=. (j-2) , liosr1b
      liosc2b=. }: liosc2b
      j=. <: j
    end.
    liosr1a=. }. liosr1a
    i=. >: i
  end.
  y ; dQ ; + dZ
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
NB.     A * x = λ * B * x
NB.   to its standard form:
NB.     H * y = λ * T * y
NB.     y = Z^_1 * x
NB.   and accumulate rotations to form Q and Z later
NB.
NB. Syntax:
NB.   'HT dQ dZ'=. hs gghrdu A ,: B
NB. where
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines submatrices A11 and B11 to be reduced
NB.           position in matrices A and B, respectively, see
NB.           gehrdu
NB.   A     - n×n-matrix, general
NB.   B     - n×n-matrix, upper triangular
NB.   HT    -: H ,: T
NB.   H     - n×n-matrix, upper Hessenberg inside the
NB.           submatrix H[h:h+s-1,h:h+s-1], and upper
NB.           triangular outside
NB.   T     - n×n-matrix, upper triangular
NB.   dQ,dZ - any×4-matrix, accumulates rotations to form Q
NB.           and Z later, see rotsclx; dQ and dZ may have
NB.           the same shapes
NB.
NB. References:
NB. [1] G. H. Golub and C. F. Van Loan, Matrix Computations,
NB.     Johns Hopkins University Press, Baltimore, Md, USA,
NB.     3rd edition, 1996, p. 378

gghrdu=: 4 : 0
  'h s'=. x
  t=. h+s-1
  n=. c y
  dQ=. dZ=. 0 4 $ 0
  j=. h
  liosc1a=. n th2lios j                        NB. (n-h)-vector h:n-1
  liosr2a=. i. h+s                             NB. (h+s)-vector 0:h+s-1
  while. j < <: t do.                          NB. (s-2)-vector h:h+s-3
    i=. t
    liosc1b=. n th2lios <: i                   NB. (n-h-s+2)-vector h+s-2:n-1
    liosr2b=. i. >: i                          NB. (i+1)-vector 0:h+s-1
    while. i > >: j do.                        NB. (h+s-j-2)-vector (desc) h+s-1:j+2
      lios=. i - 1 0
      NB. step 1: rotate rows lios to kill A[i,j]
      'y cs'=. (rot &. |:) rotga y ; (< 0 ; lios ; liosc1a) ; < < a: ; 0
      y=. (< 1 ; lios ; liosc1b) (cs & (rot &. |:)) upd y
      dQ=. dQ , cs , lios
      lios=. i - 0 1
      NB. step 2: rotate columns lios to kill B[i,i-1]
      'y cs'=. rot rotga y ; (< 1 ; liosr2b ; lios) ; _1
      y=. (< 0 ; liosr2a ; lios) (cs & rot) upd y
      dZ=. dZ , cs , lios
      NB. step 3: update IOS
      liosc1b=. (i-2) , liosc1b
      liosr2b=. }: liosr2b
      i=. <: i
    end.
    liosc1a=. }. liosc1a
    j=. >: j
  end.
  y ; (+ dQ) ; dZ
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
NB.         matrix A, see geballp and storage layout below
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
NB.           Q(i) = I - v[i]' * τ[i] * v[i] ,
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
NB.   ((idmat @ #) -: po) Q
NB.   H -: A (mp~ mp (ct @ ])) Q
NB. where
NB.   HQf=. gehrdl A
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
NB.         matrix A, see gebalup and storage layout below
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
NB.   ((idmat @ #) -: (mp~ ct)) Q
NB.   H -: A ((ct @ ]) mp mp) Q
NB. where
NB.   HQf=. gehrdu A
NB.   H=. _1 tru _1 0 }. HQf
NB.   Q=. unghru HQf
NB.
NB. Notes:
NB. - models LAPACK's xGEHRD

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
NB. Verb:           Syntax:
NB. gghrdlnn        'H T'=.     hs gghrdlnn A ,: B
NB. gghrdlnv        'H T Z'=.   hs gghrdlnv A , B ,: Z1
NB. gghrdlvn        'H T Q'=.   hs gghrdlvn A , B ,: Q1
NB. gghrdlvv        'H T Q Z'=. hs gghrdlvv A , B , Q1 ,: Z1
NB.
NB. Description:
NB.   Reduce a pair of matrices (A,B) to generalized lower
NB.   Hessenberg form (H,T) by a unitary (orthogonal)
NB.   similarity transformation by non-blocked algorithm:
NB.     Q * A * Z^_1 = H
NB.     Q * B * Z^_1 = T
NB.   in order to reduce the generalized eigenvalue problem:
NB.     x * A = λ * x * B                                 (1)
NB.   to its standard form:
NB.     y * H = λ * y * T
NB.     y = x * Q^_1
NB.   The unitary (orthogonal) matrices Q and Z are
NB.   determined as products of Givens rotations. They may
NB.   either be formed explicitly, or they may be derived
NB.   from input matrices Q1 and Z1, so that:
NB.     Q1^_1 * A * Z1 = (ΔQ*Q1)^_1 * H * (ΔZ*Z1)
NB.     Q1^_1 * B * Z1 = (ΔQ*Q1)^_1 * T * (ΔZ*Z1)
NB.   If Q1 is the unitary (orthogonal) matrix from the LQ
NB.   factorization of B in the original equation (1), then
NB.   gghrdl reduces the original problem to generalized
NB.   Hessenberg form
NB. where
NB.   hs   - 2-vector of integers (h,s) 'head' and 'size',
NB.          defines submatrices A11 and B11 to be reduced
NB.          position in matrices A and B, respectively, see
NB.          geballp and gehrdl
NB.   A    - n×n-matrix, general
NB.   B    - n×n-matrix, lower triangular
NB.   H    - n×n-matrix, lower Hessenberg inside the submatrix
NB.          H[h:h+s-1,h:h+s-1], and lower triangular outside
NB.   T    - n×n-matrix, lower triangular
NB.   Q1   - n×n-matrix, the unitary (orthogonal), typically
NB.          from the LQ factorization of B
NB.   Q    - n×n-matrix (ΔQ*Q1)
NB.   Z1   - n×n-matrix, the unitary (orthogonal)
NB.   Z    - n×n-matrix (ΔZ*Z1)
NB.
NB. TODO:
NB. - implement blocked version

gghrdlnn=: [: : (0 {:: gghrdl)
gghrdlnv=: [: : (({: @ ]) ((0 {:: ]) , (rotscll (2 {:: ]))) (gghrdl }:))
gghrdlvn=: [: : (({: @ ]) ((0 {:: ]) , (rotscll (1 {:: ]))) (gghrdl }:))
gghrdlvv=: [: : ((2 }. ]) ((0 {:: ]) , (rotscll &: >"2 0 }.)) (gghrdl (2 & {.)))

NB. ---------------------------------------------------------
NB. Verb:           Syntax:
NB. gghrdunn        'H T'=.     hs gghrdunn A ,: B
NB. gghrdunv        'H T Z'=.   hs gghrdunv A , B ,: Z1
NB. gghrduvn        'H T Q'=.   hs gghrduvn A , B ,: Q1
NB. gghrduvv        'H T Q Z'=. hs gghrduvv A , B , Q1 ,: Z1
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
NB.   either be formed explicitly, or they may be derived
NB.   from input matrices Q1 and Z1, so that:
NB.     Q1 * A * Z1^_1 = (Q1*ΔQ) * H * (Z1*ΔZ)^_1
NB.     Q1 * B * Z1^_1 = (Q1*ΔQ) * T * (Z1*ΔZ)^_1
NB.   If Q1 is the unitary (orthogonal) matrix from the QR
NB.   factorization of B in the original equation (1), then
NB.   gghrdu reduces the original problem to generalized
NB.   Hessenberg form
NB. where
NB.   hs   - 2-vector of integers (h,s) 'head' and 'size',
NB.          defines submatrices A11 and B11 to be reduced
NB.          position in matrices A and B, respectively, see
NB.          gebalup and gehrdu
NB.   A    - n×n-matrix, general
NB.   B    - n×n-matrix, upper triangular
NB.   H    - n×n-matrix, upper Hessenberg inside the submatrix
NB.          H[h:h+s-1,h:h+s-1], and upper triangular outside
NB.   T    - n×n-matrix, upper triangular
NB.   Q1   - n×n-matrix, the unitary (orthogonal), typically
NB.          from the QR factorization of B
NB.   Q    - n×n-matrix (Q1*ΔQ)
NB.   Z1   - n×n-matrix, the unitary (orthogonal)
NB.   Z    - n×n-matrix (Z1*ΔZ)
NB.
NB. Notes:
NB. - gghrdunn models LAPACK's xGGHRD('N','N')
NB. - gghrdunv models LAPACK's xGGHRD('N','V')
NB. - gghrduvn models LAPACK's xGGHRD('V','N')
NB. - gghrduvv models LAPACK's xGGHRD('V','V')
NB.
NB. Application:
NB. - model LAPACK's xGGHRD('N','I'):
NB.     gghrduni=: gghrdunv (, (idmat @ c))
NB. - model LAPACK's xGGHRD('I','N'):
NB.     gghrduin=: gghrduvn (, (idmat @ c))
NB. - model LAPACK's xGGHRD('I','I'):
NB.     gghrduii=: gghrduvv (,~^:2~ (idmat @ c))
NB. - model LAPACK's xGGHRD('I','V'):
NB.     gghrduiv=: gghrduvv ((1 & A.) @ , (idmat @ c))
NB. - model LAPACK's xGGHRD('V','I'):
NB.     gghrduvi=: gghrduvv (, (idmat @ c))
NB.
NB. References:
NB. [1] G. H. Golub and C. F. Van Loan, Matrix Computations,
NB.     Johns Hopkins University Press, Baltimore, Md, USA,
NB.     3rd edition, 1996, p. 378
NB. [2] Bo Kågström, Daniel Kressner, Enrique S.
NB.     Quintana-Ortí, and Gregorio Quintana-Ortí
NB.     Blocked Algorithms for the Reduction to
NB.     Hessenberg-Triangular Form Revisited.
NB.     February 2008.
NB.     LAPACK Working Note 198
NB.     http://www.netlib.org/lapack/lawns/downloads/
NB.
NB. TODO:
NB. - implement blocked version [2]

gghrdunn=: [: : (0 {:: gghrdu)
gghrdunv=: [: : (({: @ ]) ((0 {:: ]) , (rotsclu (2 {:: ]))) (gghrdu }:))
gghrduvn=: [: : (({: @ ]) ((0 {:: ]) , (rotsclu (1 {:: ]))) (gghrdu }:))
gghrduvv=: [: : ((2 }. ]) ((0 {:: ]) , (rotsclu &: >"2 0 }.)) (gghrdu (2 & {.)))

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

testgehrd=: 3 : 0
  require :: ] '~addons/math/lapack/lapack.ijs'
  need_jlapack_ :: ] 'gehrd'

  rcond=. gecon1 y

  ('2b1100 & gehrd_jlapack_' tmonad ((];1:;#)`(,&>/)`(rcond"_)`(_."_)`((norm1@(- (((_1 & tru)@:(}:  )) (] mp  (mp  ct)) unghru)))%((FP_EPS*#*norm1)@[)))) y

  ('gehrdl'                  tdyad ((0,#)`]       `]     `(rcond"_)`(_."_)`((norm1@(- ((( 1 & trl)@:(}:"1)) (] mp~ (mp~ ct)) unghrl)))%((FP_EPS*#*norm1)@[)))) y
  ('gehrdu'                  tdyad ((0,#)`]       `]     `(rcond"_)`(_."_)`((norm1@(- (((_1 & tru)@:(}:  )) (] mp  (mp  ct)) unghru)))%((FP_EPS*#*norm1)@[)))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testgghrd
NB.
NB. Description:
NB.   Test Hessenberg reduction algorithms:
NB.   - gghrdx (math/mt addon)
NB.   by general matrices given
NB.
NB. Syntax:
NB.   testgehrd AB
NB. where
NB.   AB - 2×n×n-report
NB.
NB. Formula:
NB.   berr := max(berr0,berr1,berr2,berr3)
NB. where
NB.   ||M|| := max(||M||_1 , FP_SFMIN)
NB.   β - machine precision
NB.   - gghrdl:
NB.       berr0 := ||Q1^_1 * A * Z1 - Q^_1 * H * Z|| / (β * ||Q1^_1 * A * Z1|| * n)
NB.       berr1 := ||Q1^_1 * L * Z1 - Q^_1 * T * Z|| / (β * ||Q1^_1 * L * Z1|| * n)
NB.       berr2 := ||I - Q^_1 * Q|| / (β * n)
NB.       berr3 := ||I - Z^_1 * Z|| / (β * n)
NB.       (L,Z1) := B
NB.       Q := ΔQ*Q1
NB.       Z := ΔZ*Z1
NB.   - gghrdu:
NB.       berr0 := ||Q1 * A * Z1^_1 - Q * H * Z^_1|| / (β * ||Q1 * A * Z1^_1|| * n)
NB.       berr1 := ||Q1 * R * Z1^_1 - Q * T * Z^_1|| / (β * ||Q1 * R * Z1^_1|| * n)
NB.       berr2 := ||I - Q * Q^_1|| / (β * n)
NB.       berr3 := ||I - Z * Z^_1|| / (β * n)
NB.       (Q1,R) := B
NB.       Q := Q1*ΔQ
NB.       Z := Z1*ΔZ
NB.
NB. Notes:
NB. - B is reconstructed from either (L*Z1) for gghrdl, or
NB.   (Q1*R) for gghrdu, this unifies calculations and
NB.   eliminates gelqf (geqrf) round-off errors

testgghrd=: 3 : 0
  prepl=. ((,~ <) ((3&{) mp~"2 (2&{.)))~ _2&(<\)                                                  NB. L: 'AB HT QZ'=. (A,L,I,:Z1) prepl (H,T,Q,:Z)
  prepu=. ((,~ <) ((2&{) mp "2 (2&{.)))~ _2&(<\)                                                  NB. R: 'AB HT QZ'=. (A,R,Q1,:I) prepu (H,T,Q,:Z)
  safenorm=. FP_SFMIN >. norm1"2                                                                  NB. compute 1-norm safely: ||M|| := max(||M||_1 , FP_SFMIN)
  cdiff1=: 2 : '(0 & {::) safenorm@:- ((((u@{.@]) mp"2 (mp"2 (v@{:)))&>/)@}.)'                    NB. L: (ct cdiff1 ]) : ||updA - Q^_1 * H * Z|| , ||B - Q^_1 * T * Z||
                                                                                                  NB. R: (] cdiff1 ct) : ||updA - Q * H * Z^_1|| , ||B - Q * T * Z^_1||
  adiff2=: 1 : '(safenorm @ (<: upddiag) @ (u ct)"2) @ (2 & {::)'                                 NB. L: (mp~ adiff2) : ||I - Q^_1 * Q|| , ||I - Z^_1 * Z||
                                                                                                  NB. R: (mp  adiff2) : ||I - Q * Q^_1|| , ||I - Z * Z^_1||
  denom1=. safenorm @ (0 & {::)                                                                   NB. ||updA|| , ||B||
  getn=. c @ (0 & {::)                                                                            NB. n
  safediv=. ((({:<.(%/@}:))`((<./@(}:*(1,{:)))%(1&{))@.(1>(1&{)))`(%/@}:)@.(</@}:))%(FP_PREC*{:)  NB. compute u%d safely: u_by_d=. safediv (u,d,n)
  cberr01=. 2 : 'safediv"1 @: ((u cdiff1 v) ,. denom1 ,. getn)'                                   NB. L: (ct cberr01 ]) : (berr0 , berr1) for L
                                                                                                  NB. R: (] cberr01 ct) : (berr0 , berr1) for R
  aberr23=. 1 : '((<. (u adiff2))~ % (FP_PREC * ])) getn'                                         NB. L: (mp~ aberr23) : (berr2 , berr3) for L
                                                                                                  NB. R: (mp  aberr23) : (berr2 , berr3) for R
  berrl=: (>./ @ ((ct cberr01 ]) , (mp~ aberr23)) @ prepl) f.
  berru=: (>./ @ ((] cberr01 ct) , (mp  aberr23)) @ prepu) f.

  rcond=. <./ gecon1"2 y
  'L Z1 Q1 R'=. (((trl ; unglq) @ gelqf) , ((ungqr ; tru) @ geqrf)) {: y
  I=. idmat c y

  ('gghrdlnn' tdyad ((0,c)`]`]`(rcond"_)`(_."_)`(_."_))) y
  ('gghrdlnv' tdyad ((0,c)`]`]`(rcond"_)`(_."_)`(_."_))) (y,I)
  ('gghrdlvn' tdyad ((0,c)`]`]`(rcond"_)`(_."_)`(_."_))) (y,I)
  ('gghrdlvv' tdyad ((0,c)`]`]`(rcond"_)`(_."_)`berrl )) ((L 1} y),I,:Z1)
  ('gghrdunn' tdyad ((0,c)`]`]`(rcond"_)`(_."_)`(_."_))) y
  ('gghrdunv' tdyad ((0,c)`]`]`(rcond"_)`(_."_)`(_."_))) (y,I)
  ('gghrduvn' tdyad ((0,c)`]`]`(rcond"_)`(_."_)`(_."_))) (y,I)
  ('gghrduvv' tdyad ((0,c)`]`]`(rcond"_)`(_."_)`berru )) ((R 1} y),Q1,:I)

  erase 'cdiff1 adiff2 berrl berru'

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

testhrd=: 1 : 'EMPTY_mt_ [ ((testgghrd_mt_ @ u @ (2&,)) [ (testgehrd_mt_ @ u)) ^: (=/)'
