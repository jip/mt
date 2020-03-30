NB. Orthogonal factorization
NB.
NB. gelqf     LQ factorization of a general matrix
NB. geqlf     QL factorization of a general matrix
NB. geqrf     QR factorization of a general matrix
NB. gerqf     RQ factorization of a general matrix
NB.
NB. tzlzf     LZ factorization of a trapezoidal matrix
NB. tzzlf     ZL factorization of a trapezoidal matrix
NB. tzzrf     ZR factorization of a trapezoidal matrix
NB. tzrzf     RZ factorization of a trapezoidal matrix
NB.
NB. testgeqf  Test gexxf by general matrix
NB. testtzqf  Test tzxxf by trapezoidal matrix
NB. testqf    Adv. to make verb to test gexxf and tzxxf by
NB.           matrix of generator and shape given
NB.
NB. Version: 0.10.5 2020-03-30
NB.
NB. Copyright 2010-2020 Igor Zhuravlov
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

QFNB=: 32   NB. block size limit
QFNX=: 128  NB. crossover point, QFNX ≥ QFNB

NB. ---------------------------------------------------------
NB. gelq2
NB.
NB. Description:
NB.   LQ-factorization of the augmented input matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   LQf=. gelq2 eA
NB. where
NB.   eA    - m×(n+1)-matrix, being A augmented by trash
NB.           vector:
NB.             eA -: A ,. Trash
NB.   A     - m×n-matrix, the input to factorize
NB.   Trash - m-vector, will be replaced by Tau
NB.   LQf   - m×(n+1)-matrix, combined L and Qf (unit
NB.           diagonal not stored)
NB.   L     - m×k-matrix, lower triangular
NB.   Qf    - k×(n+1)-matrix, unit upper triangular, the Q
NB.           represented in factored form
NB.   Q     - k×n-matrix with orthonormal rows, which is
NB.           defined as the first k rows of a product of m
NB.           elementary reflectors H(i) of order n:
NB.             Q = Π{H(i)',i=m-1:0}
NB.           where
NB.             H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.             H(m-1:k) ≡ H(u(m-1:k),τ(m-1:k)) = H(0,0) = I
NB.   k     = min(m,n)
NB.
NB. Storage layout for m=3, n=7:
NB.   LQf:
NB.     (  l  v0 v0 v0 v0 v0 v0 τ0  )
NB.     (  l  l  v1 v1 v1 v1 v1 τ1  )
NB.     (  l  l  l  v2 v2 v2 v2 τ2  )
NB.   L:
NB.     (  l  0  0                  )
NB.     (  l  l  0                  )
NB.     (  l  l  l                  )
NB.   Qf:
NB.     (  1  v0 v0 v0 v0 v0 v0 τ0  )
NB.     (  0  1  v1 v1 v1 v1 v1 τ1  )
NB.     (  0  0  1  v2 v2 v2 v2 τ2  )
NB. where
NB.   l              - elements of L
NB.   vi             - vector v(i)
NB.   τi             - scalar value conj(τ(i))
NB.   (0,...,0,1,vi) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = L * (Q)
NB.   (}:"1 -: (unmlqrn  trlpick        @:(}:"1))@gelq2) eA
NB.   NB. I = Q * (Q^H)
NB.   (( 0         idmat (<. , ])/)@(0 _1 + $) (-: clean) (unmlqrc unglq)@gelq2) eA
NB.
NB. Notes:
NB. - models LAPACK's xGELQ2
NB. - gelq2 and gelqf are topologic equivalents
NB. - if L diagonal's non-negativity is required, then
NB.   larfgfc should be replaced by larfpfc

gelq2=: 3 : 0
  pfx=. 0 {."1 y
  sfxT=. 0 {. y
  while. -. 0 e. 0 _1 + $ y do.
    z=. larfgfc {. y
    sfxT=. sfxT , z
    y=. ((1) 0} z) larfrnfr }. y
    pfx=. pfx ,. sfxT ,&:({."1) y
    sfxT=. }."1 sfxT
    y=. }."1 y
  end.
  pfx stitcht sfxT
)

NB. ---------------------------------------------------------
NB. geql2
NB.
NB. Description:
NB.   QL-factorization of the augmented input matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   QfL=. geql2 eA
NB. where
NB.   eA    - (m+1)×n-matrix, being A augmented by trash
NB.           vector:
NB.             eA -: Trash , A
NB.   Trash - n-vector, will be replaced by Tau
NB.   A     - m×n-matrix, the input to factorize
NB.   QfL   - (m+1)×n-matrix, combined Qf (unit diagonal not
NB.           stored) and L
NB.   Qf    - (m+1)×k-matrix, unit upper triangular, the Q
NB.           represented in factored form
NB.   L     - k×n-matrix, lower triangular
NB.   Q     - m×k-matrix with orthonormal columns, which is
NB.           defined as the last k columns of a product of n
NB.           elementary reflectors H(i) of order m:
NB.             Q = Π{H(i),i=n-1:0}
NB.           where
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.             H(n-1:k) ≡ H(u(n-1:k),τ(n-1:k)) = H(0,0) = I
NB.   k     = min(m,n)
NB.
NB. Storage layout for m=7, n=3:
NB.   QfL:                Qf:                 L:
NB.     (  τ0 τ1 τ2  )      (  τ0 τ1 τ2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  l  v1 v2  )      (  1  v1 v2  )      (  l  0  0   )
NB.     (  l  l  v2  )      (  0  1  v2  )      (  l  l  0   )
NB.     (  l  l  l   )      (  0  0  1   )      (  l  l  l   )
NB. where
NB.   l              - elements of L
NB.   vi             - vector v(i)
NB.   τi             - scalar value τ(i)
NB.   (vi,1,0,...,0) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = (Q) * L
NB.   (}.   -: (unmqlln (trlpick~ -~/@$)@  }.   )@geql2) eA
NB.   NB. I = (Q^H) * Q
NB.   (((0 <. -~/) idmat ([ , <.)/)@(_1 0 + $) (-: clean) (unmqllc ungql)@geql2) eA
NB.
NB. Notes:
NB. - models LAPACK's xGEQL2
NB. - geql2 and geqlf are topologic equivalents
NB. - if L diagonal's non-negativity is required, then larfgb
NB.   should be replaced by larfpb

geql2=: 3 : 0
  pfxR=. 0 {."1 y
  sfx=. 0 {. y
  while. -. 0 e. _1 0 + $ y do.
    z=. larfgb {:"1 y
    pfxR=. z ,. pfxR
    y=. ((1) _1} z) larflcbc }:"1 y
    sfx=. sfx ,~ y ,&{: pfxR
    y=. }: y
    pfxR=. }: pfxR
  end.
  pfxR appendr sfx
)

NB. ---------------------------------------------------------
NB. geqr2
NB.
NB. Description:
NB.   QR-factorization of the augmented input matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   QfR=. geqr2 eA
NB. where
NB.   eA    - (m+1)×n-matrix, being A augmented by trash
NB.           vector:
NB.             eA -: A , Trash
NB.   A     - m×n-matrix, the input to factorize
NB.   Trash - n-vector, will be replaced by Tau
NB.   QfR   - (m+1)×n-matrix, combined Qf (unit diagonal not
NB.           stored) and R
NB.   Qf    - (m+1)×k-matrix, unit lower triangular, the Q
NB.           represented in factored form
NB.   R     - k×n-matrix, upper triangular
NB.   Q     - m×k-matrix with orthonormal columns, which is
NB.           defined as the first k columns of a product of n
NB.           elementary reflectors H(i) of order m:
NB.             Q = Π{H(i),i=0:n-1}
NB.           where
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.             H(k:n-1) ≡ H(u(k:n-1),τ(k:n-1)) = H(0,0) = I
NB.   k   = min(m,n)
NB.
NB. Storage layout for m=7, n=3:
NB.   QfR:                Qf:                 R:
NB.     (  r  r  r   )      (  1  0  0   )      (  r  r  r   )
NB.     (  v0 r  r   )      (  v0 1  0   )      (  0  r  r   )
NB.     (  v0 v1 r   )      (  v0 v1 1   )      (  0  0  r   )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  τ0 τ1 τ2  )      (  τ0 τ1 τ2  )
NB. where
NB.   r              - elements of R
NB.   vi             - vector v(i)
NB.   τi             - scalar value τ(i)
NB.   (0,...,0,1,vi) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = (Q) * R
NB.   (}:   -: (unmqrln  trupick        @  }:   )@geqr2) eA
NB.   NB. I = (Q^H) * Q
NB.   (( 0         idmat ([ , <.)/)@(_1 0 + $) (-: clean) (unmqrlc ungqr)@geqr2) eA
NB.
NB. Notes:
NB. - models LAPACK's xGEQR2
NB. - gerq2 and geqrf are topologic equivalents
NB. - if R diagonal's non-negativity is required, then larfgf
NB.   should be replaced by larfpf

geqr2=: 3 : 0
  pfx=. 0 {. y
  sfxL=. 0 {."1 y
  while. -. 0 e. _1 0 + $ y do.
    z=. larfgf {."1 y
    sfxL=. sfxL ,. z
    y=. ((1) 0} z) larflcfc }."1 y
    pfx=. pfx , sfxL ,&{. y
    sfxL=. }. sfxL
    y=. }. y
  end.
  pfx , sfxL
)

NB. ---------------------------------------------------------
NB. gerq2
NB.
NB. Description:
NB.   RQ-factorization of the augmented input matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   RQf=. gerq2 eA
NB. where
NB.   eA    - m×(n+1)-matrix, being A augmented by trash
NB.           vector:
NB.             eA -: Trash ,. A
NB.   Trash - m-vector, will be replaced by Tau
NB.   A     - m×n-matrix, the input to factorize
NB.   RQf   - m×(n+1)-matrix, combined R and  Qf (unit
NB.           diagonal not stored)
NB.   R     - m×k-matrix, upper triangular
NB.   Qf    - k×(n+1)-matrix, unit lower triangular, the Q
NB.           represented in factored form
NB.   Q     - k×n-matrix with orthonormal rows which is
NB.           defined as the last k rows of a product of m
NB.           elementary reflectors H(i) of order n:
NB.             Q = Π{H(i)',i=0:m-1}
NB.           where
NB.             H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.             H(k:m-1) ≡ H(u(k:m-1),τ(k:m-1)) = H(0,0) = I
NB.   k     = min(m,n)
NB.
NB. Storage layout for m=3, n=7:
NB.   RQf:
NB.     (  τ0 v0 v0 v0 v0 r  r  r   )
NB.     (  τ1 v1 v1 v1 v1 v1 r  r   )
NB.     (  τ2 v2 v2 v2 v2 v2 v2 r   )
NB.   R:
NB.     (                 r  r  r   )
NB.     (                 0  r  r   )
NB.     (                 0  0  r   )
NB.   Qf:
NB.     (  τ0 v0 v0 v0 v0 1  0  0   )
NB.     (  τ1 v1 v1 v1 v1 v1 1  0   )
NB.     (  τ2 v2 v2 v2 v2 v2 v2 1   )
NB. where
NB.   r              - elements of R
NB.   vi             - vector v(i)
NB.   τi             - scalar value conj(τ(i))
NB.   (vi,1,0,...,0) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = R * (Q)
NB.   (}."1 -: (unmrqrn (trupick~ -~/@$)@:(}."1))@gerq2) eA
NB.   NB. I = Q * (Q^H)
NB.   (((0 >. -~/) idmat (<. , ])/)@(0 _1 + $) (-: clean) (unmrqrc ungrq)@gerq2) eA
NB.
NB. Notes:
NB. - models LAPACK's xGERQ2
NB. - gerq2 and gerqf are topologic equivalents
NB. - if R diagonal's non-negativity is required, then
NB.   larfgbc should be replaced by larfpbc

gerq2=: 3 : 0
  pfxB=. 0 {. y
  sfx=. 0 {."1 y
  while. -. 0 e. 0 _1 + $ y do.
    z=. larfgbc {: y
    pfxB=. z , pfxB
    y=. ((1) _1} z) larfrnbr }: y
    sfx=. sfx ,.~ y ,&:({:"1) pfxB
    y=. }:"1 y
    pfxB=. }:"1 pfxB
  end.
  pfxB stitchb sfx
)

NB. ---------------------------------------------------------
NB. latlz
NB.
NB. Description:
NB.   LZ-factorization of the augmented trapezoidal matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   LZf=. l latlz eA
NB. where
NB.   eA    - m×(n+1)-matrix, being A augmented by trash
NB.           vector:
NB.             eA -: Trash ,. iA1 ,. A2 ,. iL
NB.   Trash - m-vector, will be replaced by Tau
NB.   iA1   - m×(l-1)-matrix, part to be replaced by Zf
NB.   A2    - m×(n-m-l+1)-matrix, is not changed
NB.   iL    - m×m-matrix, lower triangular
NB.   LZf   - m×(n+1)-matrix, combined unit lower trapezoidal
NB.           Zf and lower trapezoidal L:
NB.             LZf -: Tau ,. oA1 ,. A2 ,. oL
NB.   Tau   - m-vector, scalars τ[0:m-1] for Zf
NB.   oA1   - m×(l-1)-matrix, rows are vectors v[0:m-1] for
NB.           Zf
NB.   oL    - m×m-matrix, lower triangular
NB.   Zf    - m×(n+1)-matrix, the Z represented in factored
NB.           form
NB.   Z     - n×n-matrix with orthonormal rows which is
NB.           defined as the product of m elementary
NB.           reflectors H(i) of order n:
NB.             Q = Π{H(i)',i=m-1:0}
NB.           where
NB.             H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.
NB. Storage layout for m=3, n=9, l=5:
NB.   eA:
NB.     (  *  a1 a1 a1 a1 a2 a2 il 0  0   )
NB.     (  *  a1 a1 a1 a1 a2 a2 il li 0   )
NB.     (  *  a1 a1 a1 a1 a2 a2 il il il  )
NB.   LZf:
NB.     (  τ0 v0 v0 v0 v0 a2 a2 ol 0  0   )
NB.     (  τ1 v1 v1 v1 v1 a2 a2 ol ol 0   )
NB.     (  τ2 v2 v2 v2 v2 a2 a2 ol ol ol  )
NB.   L:
NB.     (     0  0  0  0  a2 a2 ol 0  0   )
NB.     (     0  0  0  0  a2 a2 ol ol 0   )
NB.     (     0  0  0  0  a2 a2 ol ol ol  )
NB.   Zf:
NB.     (  τ0 v0 v0 v0 v0 0  0  1  0  0   )
NB.     (  τ1 v1 v1 v1 v1 0  0  0  1  0   )
NB.     (  τ2 v2 v2 v2 v2 0  0  0  0  1   )
NB. where
NB.   *                    - any scalar value, is not used
NB.   a1                   - elements of iA1
NB.   a2                   - elements of A2
NB.   il                   - elements of iL
NB.   ol                   - elements of oL
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value conj(τ(i))
NB.   (vi,0,...0,1,0,..,0) - n-vector u(i)
NB.
NB. Notes:
NB. - latlz and tzlzf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0

latlz=: 4 : 0
  iso=. < < < (dhs2liso 0 , x) , _1
  m=. # y
  sfxR=. (- m) {."1 y
  y=. (- m) }."1 y
  pfx=. 0 {. y
  while. # y do.
    y=. y ,. {."1 sfxR
    sfxR=. 1 1 }. sfxR
    z=. {. y
    bak=. iso { z
    z=. larfgbc 0 iso} z
    y=. ((1) _1} z) larzrnfr }. y
    pfx=. (pfx ,. 0) , bak iso} z
  end.
  pfx
)

NB. ---------------------------------------------------------
NB. latzl
NB.
NB. Description:
NB.   ZL-factorization of the augmented trapezoidal matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   ZfL=. l latzl eA
NB. where
NB.   eA    - (m+1)×n-matrix, being A augmented by trash
NB.           vector:
NB.             eA -: iL , A1 , iA2 , Trash
NB.   iL    - n×n-matrix, lower triangular
NB.   A1    - (m-n-l+1)×n-matrix, is not changed
NB.   iA2   - (l-1)×n-matrix, part to be replaced by Zf
NB.   Trash - n-vector, will be replaced by Tau
NB.   ZfL   - (m+1)×n-matrix, combined lower trapezoidal L
NB.           and unit lower trapezoidal Zf:
NB.             ZfL -: oL , A1 , oA2 , Tau
NB.   oL    - n×n-matrix, lower triangular
NB.   oA2   - (l-1)×n-matrix, columns are vectors v[0:n-1]
NB.           for Zf
NB.   Tau   - n-vector, scalars τ[0:n-1] for Zf
NB.   Zf    - (m+1)×n-matrix, the Z represented in factored
NB.           form
NB.   Z     - m×m-matrix with orthonormal columns which is
NB.           defined as the product of n elementary
NB.           reflectors H(i) of order m:
NB.             Z = Π{H(i),i=n-1:0}
NB.           where
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.
NB. Storage layout for m=9, n=3, l=5:
NB.   eA:               ZfL:              Zf:               L:
NB.   (  il 0  0   )    (  ol 0  0   )    (  1  0  0   )    (  ol 0  0   )
NB.   (  il il 0   )    (  ol ol 0   )    (  0  1  0   )    (  ol ol 0   )
NB.   (  il il il  )    (  ol ol ol  )    (  0  0  1   )    (  ol ol ol  )
NB.   (  a1 a1 a1  )    (  a1 a1 a1  )    (  0  0  0   )    (  a1 a1 a1  )
NB.   (  a1 a1 a1  )    (  a1 a1 a1  )    (  0  0  0   )    (  a1 a1 a1  )
NB.   (  a2 a2 a2  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a2 a2 a2  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a2 a2 a2  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a2 a2 a2  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  *  *  *   )    (  τ0 τ1 τ2  )    (  τ0 τ1 τ2  )
NB. where
NB.   il                   - elements of iL
NB.   a1                   - elements of A1
NB.   a2                   - elements of iA2
NB.   *                    - any scalar value, is not used
NB.   ol                   - elements of oL
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value τ(i)
NB.   (0,...0,1,0,..,0,vi) - m-vector u(i)
NB.                          elementary re
NB.
NB. Notes:
NB. - latzl and tzzlf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0

latzl=: 4 : 0
  iso=. < < < 0 , dhs2liso _1 , x
  n=. c y
  pfxT=. n {. y
  y=. n }. y
  sfx=. 0 {."1 y
  while. c y do.
    y=. ({: pfxT) , y
    pfxT=. _1 _1 }. pfxT
    z=. {:"1 y
    bak=. iso { z
    z=. larfgf 0 iso} z
    y=. ((1) 0} z) larzlcbc }:"1 y
    sfx=. (bak iso} z) ,. 0 , sfx
  end.
  sfx
)

NB. ---------------------------------------------------------
NB. latzr
NB.
NB. Description:
NB.   ZR-factorization of the augmented trapezoidal matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   ZfR=. l latzr eA
NB. where
NB.   eA    - (m+1)×n-matrix, being A augmented by trash
NB.           vector:
NB.             eA -: Trash , iA1 , A2 , iR
NB.   Trash - n-vector, will be replaced by Tau
NB.   iA1   - (l-1)×n-matrix, part to be replaced by Zf
NB.   A2    - (m-n-l+1)×n-matrix, is not changed
NB.   iR    - n×n-matrix, upper triangular
NB.   ZfR   - (m+1)×n-matrix, combined upper trapezoidal R
NB.           and unit upper trapezoidal Zf:
NB.             ZfR -: Tau , oA1 , A2 , oR
NB.   Tau   - n-vector, scalars τ[0:n-1] for Zf
NB.   oA1   - (l-1)×n-matrix, columns are vectors v[0:n-1]
NB.           for Zf
NB.   oR    - n×n-matrix, upper triangular
NB.   Zf    - (m+1)×n-matrix, the Z represented in factored
NB.           form
NB.   Z     - m×m-matrix with orthonormal columns which is
NB.           defined as the product of n elementary
NB.           reflectors H(i) of order m:
NB.             Z = Π{H(i),i=0:n-1}
NB.           where
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.
NB. Storage layout for m=9, n=3, l=5:
NB.   eA:               ZfR:              Zf:               R:
NB.   (  *  *  *   )    (  τ0 τ1 τ2  )    (  τ0 τ1 τ2  )
NB.   (  a1 a1 a1  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a1 a1 a1  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a1 a1 a1  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a1 a1 a1  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a2 a2 a2  )    (  a2 a2 a2  )    (  0  0  0   )    (  a2 a2 a2  )
NB.   (  a2 a2 a2  )    (  a2 a2 a2  )    (  0  0  0   )    (  a2 a2 a2  )
NB.   (  ir ir ir  )    (  or or or  )    (  1  0  0   )    (  or or or  )
NB.   (  0  ir ir  )    (  0  or or  )    (  0  1  0   )    (  0  or or  )
NB.   (  0  0  ir  )    (  0  0  or  )    (  0  0  1   )    (  0  0  or  )
NB. where
NB.   *                    - any scalar value, is not used
NB.   a1                   - elements of iA1
NB.   a2                   - elements of A2
NB.   ir                   - elements of iR
NB.   or                   - elements of oR
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value τ(i)
NB.   (vi,0,...0,1,0,..,0) - m-vector u(i)
NB.
NB. Notes:
NB. - latzr and tzzrf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0

latzr=: 4 : 0
  iso=. < < < (dhs2liso 0 , x) , _1
  n=. c y
  sfxB=. (- n) {. y
  y=. (- n) }. y
  pfx=. 0 {."1 y
  while. c y do.
    y=. y , {. sfxB
    sfxB=. 1 1 }. sfxB
    z=. {."1 y
    bak=. iso { z
    z=. larfgb 0 iso} z
    y=. ((1) _1} z) larzlcfc }."1 y
    pfx=. (pfx , 0) ,. bak iso} z
  end.
  pfx
)

NB. ---------------------------------------------------------
NB. latrz
NB.
NB. Description:
NB.   RZ-factorization of the augmented trapezoidal matrix by
NB.   non-blocked version of algorithm
NB.
NB. Syntax:
NB.   RZf=. l latrz eA
NB. where
NB.   eA    - m×(n+1)-matrix, being A augmented by trash
NB.           vector:
NB.             eA -: iR ,. A1 ,. iA2 ,. Trash
NB.   iR    - m×m-matrix, upper triangular
NB.   A1    - m×(n-m-l+1)-matrix, is not changed
NB.   iA2   - m×(l-1)-matrix, part to be replaced by Zf
NB.   Trash - m-vector, will be replaced by Tau
NB.   RZf   - m×(n+1)-matrix, combined upper trapezoidal R
NB.           and unit upper trapezoidal Zf:
NB.             RZf -: oR ,. A1 ,. oA2 ,. Tau
NB.   oR    - m×m-matrix, upper triangular
NB.   oA2   - m×(l-1)-matrix, rows are vectors v[0:m-1] for
NB.           Zf
NB.   Tau   - m-vector, scalars τ[0:m-1] for Zf
NB.   Zf    - m×(n+1)-matrix, the Z represented in factored
NB.           form
NB.   Z     - n×n-matrix with orthonormal rows which is
NB.           defined as the product of m elementary
NB.           reflectors H(i) of order n:
NB.             Z = Π{H(i)',i=0:m-1}
NB.           where
NB.             H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.
NB. Storage layout for m=3, n=9, l=5:
NB.   eA:
NB.     (  ir ir ir a1 a1 a2 a2 a2 a2 *   )
NB.     (  0  ir ir a1 a1 a2 a2 a2 a2 *   )
NB.     (  0  0  ir a1 a1 a2 a2 a2 a2 *   )
NB.   RZf:
NB.     (  or or or a1 a1 v0 v0 v0 v0 τ0  )
NB.     (  0  or or a1 a1 v1 v1 v1 v1 τ1  )
NB.     (  0  0  or a1 a1 v2 v2 v2 v2 τ2  )
NB.   R:
NB.     (  or or or a1 a1 0  0  0  0      )
NB.     (  0  or or a1 a1 0  0  0  0      )
NB.     (  0  0  or a1 a1 0  0  0  0      )
NB.   Zf:
NB.     (  1  0  0  0  0  v0 v0 v0 v0 τ0  )
NB.     (  0  1  0  0  0  v1 v1 v1 v1 τ1  )
NB.     (  0  0  1  0  0  v2 v2 v2 v2 τ2  )
NB. where
NB.   ir                   - elements of iR
NB.   a1                   - elements of A1
NB.   a2                   - elements of iA2
NB.   *                    - any scalar value, is not used
NB.   or                   - elements of oR
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value conj(τ(i))
NB.   (0,...0,1,0,..,0,vi) - n-vector u(i)
NB.
NB. Notes:
NB. - models LAPACK's xLATRZ with the following difference:
NB.   - v(i) is saved instead of conj(v(i))
NB.   - conj(τ(i)) is saved instead of τ(i)
NB.   to keep consistence with gexxf
NB. - latrz and tzrzf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0

latrz=: 4 : 0
  iso=. < < < 0 , dhs2liso _1 , x
  m=. # y
  pfxL=. m {."1 y
  y=. m }."1 y
  sfx=. 0 {. y
  while. # y do.
    y=. ({:"1 pfxL) ,. y
    pfxL=. _1 _1 }. pfxL
    z=. {: y
    bak=. iso { z
    z=. larfgfc 0 iso} z
    y=. ((1) 0} z) larzrnbr }: y
    sfx=. (bak iso} z) , 0 ,. sfx
  end.
  sfx
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gelqf
NB.
NB. Description:
NB.   LQ factorization of a general matrix
NB.
NB. Syntax:
NB.   LQf=. gelqf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   LQf - m×(n+1)-matrix, combined L and Qf (unit
NB.         diagonal not stored)
NB.   L   - m×k-matrix, lower triangular
NB.   Qf  - k×(n+1)-matrix, unit upper triangular, the Q
NB.         represented in factored form
NB.   Q   - k×n-matrix with orthonormal rows, which is
NB.         defined as the first k rows of a product of m
NB.         elementary reflectors H(i) of order n:
NB.           Q = Π{H(i)',i=m-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.           H(m-1:k) ≡ H(u(m-1:k),τ(m-1:k)) = H(0,0)=I
NB.   k   = min(m,n)
NB.
NB. Storage layout for m=3, n=7:
NB.   LQf:
NB.     (  l  v0 v0 v0 v0 v0 v0 τ0  )
NB.     (  l  l  v1 v1 v1 v1 v1 τ1  )
NB.     (  l  l  l  v2 v2 v2 v2 τ2  )
NB.   L:
NB.     (  l  0  0                  )
NB.     (  l  l  0                  )
NB.     (  l  l  l                  )
NB.   Qf:
NB.     (  1  v0 v0 v0 v0 v0 v0 τ0  )
NB.     (  0  1  v1 v1 v1 v1 v1 τ1  )
NB.     (  0  0  1  v2 v2 v2 v2 τ2  )
NB. where
NB.   l              - elements of L
NB.   vi             - vector v(i)
NB.   τi             - scalar value conj(τ(i))
NB.   (0,...,0,1,vi) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = L * (Q)
NB.   (     -: (unmlqrn  trlpick        @:(}:"1))@gelqf) A
NB.   NB. I = Q * (Q^H)
NB.   (( 0         idmat (<. , ])/)@        $  (-: clean) (unmlqrc unglq)@gelqf) A
NB.
NB. Notes:
NB. - models LAPACK's xGELQF
NB. - gelq2 and gelqf are topologic equivalents

gelqf=: 3 : 0
  y=. y ,. 0
  pfx=. 0 {."1 y
  sfxT=. 0 {. y
  while. -. 0 e. QFNX < 0 _1 + $ y do.
    nb=. <./ QFNB , 0 _1 + $ y
    Z=. gelq2 nb {. y
    sfxT=. sfxT , Z
    y=. (tru1 Z) larfbrnfr nb }. y
    pfx=. pfx ,. sfxT ,&(nb&({."1)) y
    sfxT=. nb }."1 sfxT
    y=. nb }."1 y
  end.
  pfx ,. sfxT , gelq2 y
)

NB. ---------------------------------------------------------
NB. geqlf
NB.
NB. Description:
NB.   QL factorization of a general matrix
NB.
NB. Syntax:
NB.   QfL=. geqlf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   QfL - (m+1)×n-matrix, combined Qf (unit diagonal not
NB.         stored) and L
NB.   Qf  - (m+1)×k-matrix, unit upper triangular, the Q
NB.         represented in factored form
NB.   L   - k×n-matrix, lower triangular
NB.   Q   - m×k-matrix with orthonormal columns, which is
NB.         defined as the last k columns of a product of n
NB.         elementary reflectors H(i) of order m:
NB.           Q = Π{H(i),i=n-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.           H(n-1:k) ≡ H(u(n-1:k),τ(n-1:k)) = H(0,0) = I
NB.   k   = min(m,n)
NB.
NB. Storage layout for m=7, n=3:
NB.   QfL:                Qf:                 L:
NB.     (  τ0 τ1 τ2  )      (  τ0 τ1 τ2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  l  v1 v2  )      (  1  v1 v2  )      (  l  0  0   )
NB.     (  l  l  v2  )      (  0  1  v2  )      (  l  l  0   )
NB.     (  l  l  l   )      (  0  0  1   )      (  l  l  l   )
NB. where
NB.   l              - elements of L
NB.   vi             - vector v(i)
NB.   τi             - scalar value τ(i)
NB.   (vi,1,0,...,0) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = (Q) * L
NB.   (     -: (unmqlln (trlpick~ -~/@$)@  }.   )@geqlf) A
NB.   NB. I = (Q^H) * Q
NB.   (((0 <. -~/) idmat ([ , <.)/)@        $  (-: clean) (unmqllc ungql)@geqlf) A
NB.
NB. Notes:
NB. - models LAPACK's xGEQLF
NB. - geql2 and geqlf are topologic equivalents

geqlf=: 3 : 0
  y=. 0 , y
  pfxR=. 0 {."1 y
  sfx=. 0 {. y
  while. -. 0 e. QFNX < _1 0 + $ y do.
    nb=. - <./ QFNB , _1 0 + $ y
    Z=. geql2 nb {."1 y
    pfxR=. Z ,. pfxR
    y=. ((tru1~ -~/@$) Z) larfblcbc nb }."1 y
    sfx=. sfx ,~ y ,.&(nb&{.) pfxR
    y=. nb }. y
    pfxR=. nb }. pfxR
  end.
  sfx ,~ pfxR ,.~ geql2 y
)

NB. ---------------------------------------------------------
NB. geqrf
NB.
NB. Description:
NB.   QR factorization of a general matrix
NB.
NB. Syntax:
NB.   QfR=. geqrf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   QfR - (m+1)×n-matrix, combined Qf (unit diagonal not
NB.         stored) and R
NB.   Qf  - (m+1)×k-matrix, unit lower triangular, the Q
NB.         represented in factored form
NB.   R   - k×n-matrix, upper triangular
NB.   Q   - m×k-matrix with orthonormal columns, which is
NB.         defined as the first k columns of a product of n
NB.         elementary reflectors H(i) of order m:
NB.           Q = Π{H(i),i=0:n-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.           H(k:n-1) ≡ H(u(k:n-1),τ(k:n-1)) = H(0,0) = I
NB.   k   = min(m,n)
NB.
NB. Storage layout for m=7, n=3:
NB.   QfR:                Qf:                 R:
NB.     (  r  r  r   )      (  1  0  0   )      (  r  r  r   )
NB.     (  v0 r  r   )      (  v0 1  0   )      (  0  r  r   )
NB.     (  v0 v1 r   )      (  v0 v1 1   )      (  0  0  r   )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  v0 v1 v2  )      (  v0 v1 v2  )
NB.     (  τ0 τ1 τ2  )      (  τ0 τ1 τ2  )
NB. where
NB.   r              - elements of R
NB.   vi             - vector v(i)
NB.   τi             - scalar value τ(i)
NB.   (0,...,0,1,vi) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = (Q) * R
NB.   (     -: (unmqrln  trupick        @  }:   )@geqrf) A
NB.   NB. I = (Q^H) * Q
NB.   (( 0         idmat ([ , <.)/)@        $  (-: clean) (unmqrlc ungqr)@geqrf) A
NB.
NB. Notes:
NB. - models LAPACK's xGEQRF
NB. - gerq2 and geqrf are topologic equivalents

geqrf=: 3 : 0
  y=. y , 0
  pfx=. 0 {. y
  sfxL=. 0 {."1 y
  while. -. 0 e. QFNX < _1 0 + $ y do.
    nb=. <./ QFNB , _1 0 + $ y
    Z=. geqr2 nb {."1 y
    sfxL=. sfxL ,. Z
    y=. (trl1 Z) larfblcfc nb }."1 y
    pfx=. pfx , sfxL ,.&(nb&{.) y
    sfxL=. nb }. sfxL
    y=. nb }. y
  end.
  pfx , sfxL ,. geqr2 y
)

NB. ---------------------------------------------------------
NB. gerqf
NB.
NB. Description:
NB.   RQ factorization of a general matrix
NB.
NB. Syntax:
NB.   RQf=. gerqf A
NB. where
NB.   A   - m×n-matrix, the input to factorize
NB.   RQf - m×(n+1)-matrix, combined R and Qf (unit diagonal
NB.         not stored)
NB.   R   - m×k-matrix, upper triangular
NB.   Qf  - k×(n+1)-matrix, unit lower triangular, the Q
NB.         represented in factored form
NB.   Q   - k×n-matrix with orthonormal rows, which is
NB.         defined as the last k rows of a product of m
NB.         elementary reflectors H(i) of order n:
NB.           Q = Π{H(i)',i=0:m-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.           H(k:m-1) ≡ H(u(k:m-1),τ(k:m-1)) = H(0,0) = I
NB.   k   = min(m,n)
NB.
NB. Storage layout for m=3, n=7:
NB.   RQf:
NB.     (  τ0 v0 v0 v0 v0 r  r  r   )
NB.     (  τ1 v1 v1 v1 v1 v1 r  r   )
NB.     (  τ2 v2 v2 v2 v2 v2 v2 r   )
NB.   R:
NB.     (                 r  r  r   )
NB.     (                 0  r  r   )
NB.     (                 0  0  r   )
NB.   Qf:
NB.     (  τ0 v0 v0 v0 v0 1  0  0   )
NB.     (  τ1 v1 v1 v1 v1 v1 1  0   )
NB.     (  τ2 v2 v2 v2 v2 v2 v2 1   )
NB. where
NB.   r              - elements of R
NB.   vi             - vector v(i)
NB.   τi             - scalar value conj(τ(i))
NB.   (vi,1,0,...,0) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = R * (Q)
NB.   (     -: (unmrqrn (trupick~ -~/@$)@:(}."1))@gerqf) A
NB.   NB. I = Q * (Q^H)
NB.   (((0 >. -~/) idmat (<. , ])/)@        $  (-: clean) (unmrqrc ungrq)@gerqf) A
NB.
NB. Notes:
NB. - models LAPACK's xGERQF
NB. - gerq2 and gerqf are topologic equivalents

gerqf=: 3 : 0
  y=. 0 ,. y
  pfxB=. 0 {. y
  sfx=. 0 {."1 y
  while. -. 0 e. QFNX < 0 _1 + $ y do.
    nb=. - <./ QFNB , 0 _1 + $ y
    Z=. gerq2 nb {. y
    pfxB=. Z , pfxB
    y=. ((trl1~ -~/@$) Z) larfbrnbr nb }. y
    sfx=. sfx ,.~ y ,&(nb&({."1)) pfxB
    y=. nb }."1 y
    pfxB=. nb }."1 pfxB
  end.
  sfx ,.~ pfxB ,~ gerq2 y
)

NB. ---------------------------------------------------------
NB. tzlzf
NB.
NB. Description:
NB.   LZ-factorization of the trapezoidal matrix
NB.
NB. Syntax:
NB.   LZf=. tzlzf A
NB. where
NB.   A   - m×n-matrix:
NB.           A -: iA0 ,. iL
NB.   iA0 - m×(n-m)-matrix, part to be replaced by Zf
NB.   iL  - m×m-matrix, lower triangular
NB.   LZf - m×(n+1)-matrix, combined lower trapezoidal L and
NB.         unit lower trapezoidal Zf:
NB.           LZf -: Tau ,. oA0 ,. oL
NB.   Tau - m-vector, scalars τ[0:m-1] for Zf
NB.   oA0 - m×(n-m)-matrix, rows are vectors v[0:m-1] for Zf
NB.   oL  - m×m-matrix, lower triangular
NB.   Zf  - m×(n+1)-matrix, the Z represented in factored
NB.         form
NB.   Z   - n×n-matrix with orthonormal rows which is defined
NB.         as the product of m elementary reflectors H(i) of
NB.         order n:
NB.           Z = Π{H(i)',i=m-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.
NB. Storage layout for m=3, n=7:
NB.   A:
NB.     (     a0 a0 a0 a0 il 0  0   )
NB.     (     a0 a0 a0 a0 il il 0   )
NB.     (     a0 a0 a0 a0 il il il  )
NB.   LZf:
NB.     (  τ0 v0 v0 v0 v0 ol 0  0   )
NB.     (  τ1 v1 v1 v1 v1 ol ol 0   )
NB.     (  τ2 v2 v2 v2 v2 ol ol ol  )
NB.   L:
NB.     (     0  0  0  0  ol 0  0   )
NB.     (     0  0  0  0  ol ol 0   )
NB.     (     0  0  0  0  ol ol ol  )
NB.   Zf:
NB.     (  τ0 v0 v0 v0 v0 1  0  0   )
NB.     (  τ1 v1 v1 v1 v1 0  1  0   )
NB.     (  τ2 v2 v2 v2 v2 0  0  1   )
NB. where
NB.   a0                   - elements of iA0
NB.   il                   - elements of iL
NB.   ol                   - elements of oL
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value conj(τ(i))
NB.   (vi,0,...0,1,0,..,0) - n-vector u(i)
NB.
NB. Notes:
NB. - latlz and tzlzf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0

tzlzf=: 3 : 0
  y=. 0 ,. y
  l=. -~/ 'm n1'=. $ y
  sfxR=. (- m) {."1 y
  y=. (- m) }."1 y
  pfx=. 0 {. y
  nb=. QFNB <. # y
  I=. idmat nb
  O=. (2 # nb) $ 0
  while. QFNX < # y do.
    y=. y ,. nb {."1 sfxR
    sfxR=. (2 # nb) }. sfxR
    Z=. l latlz nb {. y
    y=. (I ,.~ l {."1 Z) larzbrnfr nb }. y
    I=. O ,. I
    pfx=. pfx appendl Z
  end.
  pfx appendl l latlz y ,. sfxR
)

NB. ---------------------------------------------------------
NB. tzzlf
NB.
NB. Description:
NB.   ZL-factorization of the trapezoidal matrix
NB.
NB. Syntax:
NB.   ZfL=. tzzlf A
NB. where
NB.   A   - m×n-matrix:
NB.           A -: iL , iA0
NB.   iL  - n×n-matrix, lower triangular
NB.   iA0 - (m-n)×n-matrix, part to be replaced by Zf
NB.   ZfL - (m+1)×n-matrix, combined unit upper trapezoidal
NB.         Zf and lower trapezoidal L:
NB.           ZfL -: oL , oA0 , Tau
NB.   oL  - n×n-matrix, lower triangular
NB.   oA0 - (m-n)×n-matrix, columns are vectors v[0:m-1] for
NB.         Zf
NB.   Tau - n-vector, scalars τ[0:n-1] for Zf
NB.   Zf  - (m+1)×n-matrix, the Z represented in factored
NB.         form
NB.   Z   - m×m-matrix with orthonormal columns which is
NB.         defined as the product of n elementary reflectors
NB.         of order m:
NB.           Z = Π{H(i),i=n-1:0}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.
NB. Storage layout for m=7, n=3:
NB.   A:                ZfL:              Zf:               L:
NB.   (  il 0  0   )    (  ol 0  0   )    (  1  0  0   )    (  ol 0  0   )
NB.   (  il il 0   )    (  ol ol 0   )    (  0  1  0   )    (  ol ol 0   )
NB.   (  il il il  )    (  ol ol ol  )    (  0  0  1   )    (  ol ol ol  )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.                     (  τ0 τ1 τ2  )    (  τ0 τ1 τ2  )
NB. where
NB.   il                   - elements of iL
NB.   a0                   - elements of iA0
NB.   ol                   - elements of oL
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value τ(i)
NB.   (0,...0,1,0,..,0,vi) - m-vector u(i)
NB.
NB. Notes:
NB. - latzl and tzzlf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0

tzzlf=: 3 : 0
  y=. y , 0
  l=. -/ 'm1 n'=. $ y
  pfxT=. n {. y
  y=. n }. y
  sfx=. 0 {."1 y
  nb=. - QFNB <. c y
  I=. idmat - nb
  O=. (2 # - nb) $ 0
  while. QFNX < c y do.
    y=. (nb {. pfxT) , y
    pfxT=. (2 # nb) }. pfxT
    Z=. l latzl nb {."1 y
    y=. (I , (- l) {. Z) larzblcbc nb }."1 y
    I=. I , O
    sfx=. Z stitchb sfx
  end.
  (l latzl pfxT , y) stitchb sfx
)

NB. ---------------------------------------------------------
NB. tzzrf
NB.
NB. Description:
NB.   ZR-factorization of the trapezoidal matrix
NB.
NB. Syntax:
NB.   ZfR=. tzzrf A
NB. where
NB.   A   - m×n-matrix:
NB.           A -: iA0 , iR
NB.   iA0 - (m-n)×n-matrix, part to be replaced by Zf
NB.   iR  - n×n-matrix, upper triangular
NB.   ZfR - (m+1)×n-matrix, combined unit upper trapezoidal
NB.         Zf and upper trapezoidal R:
NB.           ZfR -: Tau , oA0 , oR
NB.   Tau - n-vector, scalars τ[0:n-1] for Zf
NB.   oA0 - (m-n)×n-matrix, columns are vectors v[0:m-1] for
NB.         Zf
NB.   oR  - n×n-matrix, upper triangular
NB.   Zf  - (m+1)×n-matrix, the Z represented in factored
NB.         form
NB.   Z   - m×m-matrix with orthonormal columns which is
NB.         defined as the product of n elementary reflectors
NB.         of order m:
NB.           Z = Π{H(i),i=0:n-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * (u(i))^H
NB.
NB. Storage layout for m=7, n=3:
NB.   A:                ZfR:              Zf:               R:
NB.                     (  τ0 τ1 τ2  )    (  τ0 τ1 τ2  )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  ir ir ir  )    (  or or or  )    (  1  0  0   )    (  or or or  )
NB.   (  0  ir ir  )    (  0  or or  )    (  0  1  0   )    (  0  or or  )
NB.   (  0  0  ir  )    (  0  0  or  )    (  0  0  1   )    (  0  0  or  )
NB. where
NB.   a0                   - elements of iA0
NB.   ir                   - elements of iR
NB.   or                   - elements of oR
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value τ(i)
NB.   (vi,0,...0,1,0,..,0) - m-vector u(i)
NB.
NB. Notes:
NB. - latzr and tzzrf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0

tzzrf=: 3 : 0
  y=. 0 , y
  l=. -/ 'm1 n'=. $ y
  sfxB=. (- n) {. y
  y=. (- n) }. y
  pfx=. 0 {."1 y
  nb=. QFNB <. c y
  I=. idmat nb
  O=. (2 # nb) $ 0
  while. QFNX < c y do.
    y=. y , nb {. sfxB
    sfxB=. (2 # nb) }. sfxB
    Z=. l latzr nb {."1 y
    y=. (I ,~ l {. Z) larzblcfc nb }."1 y
    I=. O , I
    pfx=. pfx stitcht Z
  end.
  pfx stitcht l latzr y , sfxB
)

NB. ---------------------------------------------------------
NB. tzrzf
NB.
NB. Description:
NB.   RZ-factorization of the trapezoidal matrix
NB.
NB. Syntax:
NB.   RZf=. tzrzf A
NB. where
NB.   A   - m×n-matrix:
NB.           A -: iR ,. iA1
NB.   iR  - m×m-matrix, upper triangular
NB.   iA1 - m×(n-m)-matrix, part to be replaced by Zf
NB.   RZf - m×(n+1)-matrix, combined upper trapezoidal R and
NB.         unit upper trapezoidal Zf:
NB.           RZf -: oR ,. oA1 ,. Tau
NB.   oR  - m×m-matrix, upper triangular
NB.   oA1 - m×(-m)-matrix, rows are vectors v[0:m-1] for Zf
NB.   Tau - m-vector, scalars τ[0:m-1] for Zf
NB.   Zf  - m×(n+1)-matrix, the Z represented in factored
NB.         form
NB.   Z   - n×n-matrix with orthonormal rows which is defined
NB.         as the product of m elementary reflectors H(i) of
NB.         order n:
NB.           Z = Π{H(i)',i=0:m-1}
NB.         where
NB.           H(i) ≡ H(u(i),τ(i)) := I - (u(i))^H * τ(i) * u(i)
NB.
NB. Storage layout for m=3, n=7:
NB.   A:
NB.     (  ir ir ir a1 a1 a1 a1     )
NB.     (  0  ir ir a1 a1 a1 a1     )
NB.     (  0  0  ir a1 a1 a1 a1     )
NB.   RZf:
NB.     (  or or or v0 v0 v0 v0 τ0  )
NB.     (  0  or or v1 v1 v1 v1 τ1  )
NB.     (  0  0  or v2 v2 v2 v2 τ2  )
NB.   R:
NB.     (  or or or 0  0  0  0      )
NB.     (  0  or or 0  0  0  0      )
NB.     (  0  0  or 0  0  0  0      )
NB.   Zf:
NB.     (  1  0  0  v0 v0 v0 v0 τ0  )
NB.     (  0  1  0  v1 v1 v1 v1 τ1  )
NB.     (  0  0  1  v2 v2 v2 v2 τ2  )
NB. where
NB.   ir                   - elements of iR
NB.   a1                   - elements of iA1
NB.   or                   - elements of oR
NB.   vi                   - l-vector v(i)
NB.   τi                   - scalar value conj(τ(i))
NB.   (0,...0,1,0,..,0,vi) - n-vector u(i)
NB.
NB. Notes:
NB. - models LAPACK's xTZRZF with the following differences:
NB.   - v(i) is saved instead of conj(v(i))
NB.   - conj(τ(i)) is saved instead of τ(i)
NB. - latrz and tzrzf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0

tzrzf=: 3 : 0
NB.QFNB=. QFNX=. 3
  y=. y ,. 0
  l=. -/ 'm n1'=. $ y
  pfxL=. m {."1 y
  y=. m }."1 y
  sfx=. 0 {. y
  nb=. - QFNB <. # y
  I=. idmat - nb
  O=. (2 # - nb) $ 0
  while. QFNX < # y do.
    y=. (nb {."1 pfxL) ,. y
    pfxL=. (2 # nb) }. pfxL
    Z=. (- l) latrz nb {. y
    y=. (I ,. l {."1 Z) larzbrnbr nb }. y
    I=. I ,. O
    sfx=. Z appendr sfx
  end.
  ((- l) latrz pfxL ,. y) appendr sfx
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgeqf
NB.
NB. Description:
NB.   Test orthogonal factorization algorithms:
NB.   - 128!:0 (built-in)
NB.   - gelqf geqlf geqrf gerqf (math/lapack addon)
NB.   - gelqf geqlf geqrf gerqf (math/mt addon)
NB.   by general matrix
NB.
NB. Syntax:
NB.   testgeqf A
NB. where
NB.   A - m×n-matrix
NB.
NB. Formula:
NB. - for LQ: berr := max( ||L - A * Q^H|| / (FP_EPS * ||A|| * n), ||Q * Q^H - I|| / (FP_EPS * n) )
NB. - for QL: berr := max( ||L - Q^H * A|| / (FP_EPS * ||A|| * m), ||Q^H * Q - I|| / (FP_EPS * m) )
NB. - for QR: berr := max( ||R - Q^H * A|| / (FP_EPS * ||A|| * m), ||Q^H * Q - I|| / (FP_EPS * m) )
NB. - for RQ: berr := max( ||R - A * Q^H|| / (FP_EPS * ||A|| * n), ||Q * Q^H - I|| / (FP_EPS * n) )

testgeqf=: 3 : 0
  require :: ] '~addons/math/lapack/lapack.ijs'
  need_jlapack_ :: ] 'gelqf geqlf geqrf gerqf'

  rcond=. (_."_)`gecon1@.(=/@$) y  NB. meaninigful for square matrices only

  ('128!:0'                tmonad (]`]                 `(rcond"_)`(_."_)`((norm1@((1 {:: ])              - (( mp~ ct                   )  0&{::)) % (FP_EPS * (1:`]@.*)@norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@(<: upddiag)@(mp~ ct       )@(0&{::))))) y

  ('2b1110&gelqf_jlapack_' tmonad (]`({. ,  ,. &.>/@}.)`(rcond"_)`(_."_)`((norm1@((0 {:: ])              - (((   <./ @$@]) {."1 unmlqrc)~ 1&{::)) % (FP_EPS * (1:`]@.*)@norm1 * c)@[) >. ((% FP_EPS * c)~ norm1@(<: upddiag)@(unmlqrc unglq)@(1&{::))))) y
  ('2b0111&geqlf_jlapack_' tmonad (]`({: ,~ , ~&.>/@}:)`(rcond"_)`(_."_)`((norm1@((1 {:: ])              - (((-@(<./)@$@]) {.   unmqllc)~ 0&{::)) % (FP_EPS * (1:`]@.*)@norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@(<: upddiag)@(unmqllc ungql)@(0&{::))))) y
  ('2b0111&geqrf_jlapack_' tmonad (]`({: ,~ ,  &.>/@}:)`(rcond"_)`(_."_)`((norm1@((1 {:: ])              - (((   <./ @$@]) {.   unmqrlc)~ 0&{::)) % (FP_EPS * (1:`]@.*)@norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@(<: upddiag)@(unmqrlc ungqr)@(0&{::))))) y
  ('2b1110&gerqf_jlapack_' tmonad (]`({. ,  ,.~&.>/@}.)`(rcond"_)`(_."_)`((norm1@((0 {:: ])              - (((-@(<./)@$@]) {."1 unmrqrc)~ 1&{::)) % (FP_EPS * (1:`]@.*)@norm1 * c)@[) >. ((% FP_EPS * c)~ norm1@(<: upddiag)@(unmrqrc ungrq)@(1&{::))))) y

  ('gelqf'                 tmonad (]`]                 `(rcond"_)`(_."_)`((norm1@( trl        @:(}:"1)@] -  ((   <./ @$@]) {."1 unmlqrc)~       ) % (FP_EPS * (1:`]@.*)@norm1 * c)@[) >. ((% FP_EPS * c)~ norm1@(<: upddiag)@(unmlqrc unglq)        )))) y
  ('geqlf'                 tmonad (]`]                 `(rcond"_)`(_."_)`((norm1@((trl~ -~/@$)@  }.   @] -  ((-@(<./)@$@]) {.   unmqllc)~       ) % (FP_EPS * (1:`]@.*)@norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@(<: upddiag)@(unmqllc ungql)        )))) y
  ('geqrf'                 tmonad (]`]                 `(rcond"_)`(_."_)`((norm1@( tru        @  }:   @] -  ((   <./ @$@]) {.   unmqrlc)~       ) % (FP_EPS * (1:`]@.*)@norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@(<: upddiag)@(unmqrlc ungqr)        )))) y
  ('gerqf'                 tmonad (]`]                 `(rcond"_)`(_."_)`((norm1@((tru~ -~/@$)@:(}."1)@] -  ((-@(<./)@$@]) {."1 unmrqrc)~       ) % (FP_EPS * (1:`]@.*)@norm1 * c)@[) >. ((% FP_EPS * c)~ norm1@(<: upddiag)@(unmrqrc ungrq)        )))) y

  EMPTY
)

NB. ---------------------------------------------------------
NB. testtzqf
NB.
NB. Description:
NB.   Test orthogonal factorization algorithms:
NB.   - tzlzf tzzlf tzzrf tzrzf (math/mt addon)
NB.   by trapezoidal matrix
NB.
NB. Syntax:
NB.   testtzqf A
NB. where
NB.   A - m×n-matrix
NB.
NB. Formula:
NB. - for LZ: berr := max( ||A - L * Z|| / (FP_EPS * ||A|| * m), ||Z * Z^H - I|| / (FP_EPS * n) )
NB. - for ZL: berr := max( ||A - Z * L|| / (FP_EPS * ||A|| * n), ||Z^H * Z - I|| / (FP_EPS * m) )
NB. - for ZR: berr := max( ||A - Z * R|| / (FP_EPS * ||A|| * n), ||Z^H * Z - I|| / (FP_EPS * m) )
NB. - for RZ: berr := max( ||A - R * Z|| / (FP_EPS * ||A|| * m), ||Z * Z^H - I|| / (FP_EPS * n) )

testtzqf=: 3 : 0
  rcond=. (_."_)`gecon1@.(=/@$) y  NB. meaninigful for square matrices only
  Awide=. |:^:(>/@$) y
  Atall=. |:^:(</@$) y

  ('tzlzf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(- (unmlzrn (0:`((a: <@; 0 th2liso~ -~/)@[)`]}~ $)@:(}."1))) % (FP_EPS * (1:`]@.*)@norm1 * c)@[) >. ((% FP_EPS * c)~ norm1@((<: upddiag)~ 0 >. -~/@$)@(unmlzrc unglz))))) (trl~ -~/@$) Awide
  ('tzzlf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(- (unmzlln (0:`(          th2liso~/    @[)`]}~ $)@  }:   )) % (FP_EPS * (1:`]@.*)@norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@( <: upddiag             )@(unmzllc ungzl)))))  trl         Atall
  ('tzzrf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(- (unmzrln (0:`((       0 th2liso~ -~/)@[)`]}~ $)@  }.   )) % (FP_EPS * (1:`]@.*)@norm1 * #)@[) >. ((% FP_EPS * #)~ norm1@((<: upddiag)~ 0 <. -~/@$)@(unmzrlc ungzr))))) (tru~ -~/@$) Atall
  ('tzrzf' tmonad (]`]`(rcond"_)`(_."_)`((norm1@(- (unmrzrn (0:`((a: <@;   th2liso~/   )@[)`]}~ $)@:(}:"1))) % (FP_EPS * (1:`]@.*)@norm1 * c)@[) >. ((% FP_EPS * c)~ norm1@( <: upddiag             )@(unmrzrc ungrz)))))  tru         Awide

  EMPTY
)

NB. ---------------------------------------------------------
NB. testqf
NB.
NB. Description:
NB.   Adv. to make verb to test gexxx by matrix of generator
NB.   and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testqf
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
NB.     ?@$&0 testqf_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testqf_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testqf_mt_ 150 200

testqf=: 1 : 'EMPTY [ (testtzqf_mt_ [ testgeqf_mt_)@u'
