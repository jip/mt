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
NB. Version: 0.12.0 2021-02-01
NB.
NB. Copyright 2010-2021 Igor Zhuravlov
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
NB.   LQf   - m×(n+1)-matrix, L and Qf combined
NB.   L     - m×k-matrix, the lower trapezoidal
NB.   Qf    - k×(n+1)-matrix, the unit upper trapezoidal
NB.           (unit diagonal not stored), represents the Q in
NB.           factored form
NB.   Q     - n×n-matrix, the unitary (orthogonal), which is
NB.           defined as the product of k elementary
NB.           reflectors H(i) of order n:
NB.             Q = Π{H(i)',i=k-1:0}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
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
NB.   NB. A = L * Q
NB.   (}:"1 -: ( trl        @:(}:"1) mp unglq)@gelq2) eA
NB.   NB. A = L * (Q)
NB.   (}:"1 -: (unmlqrn  trlpick        @:(}:"1))@gelq2) eA
NB.
NB. Notes:
NB. - models LAPACK's xGELQ2
NB. - gelq2 and gelqf are topologic equivalents
NB. - if L diagonal's non-negativity is required, then
NB.   larfgfc would be replaced by larfpfc

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
NB.   QfL   - (m+1)×n-matrix, Qf and L combined
NB.   Qf    - (m+1)×k-matrix, the unit upper trapezoidal
NB.           (unit diagonal not stored), represents the Q in
NB.           factored form
NB.   L     - k×n-matrix, the lower trapezoidal
NB.   Q     - m×m-matrix, the unitary (orthogonal), which is
NB.           defined as the product of k elementary
NB.           reflectors H(i) of order m:
NB.             Q = Π{H(i),i=k-1:0}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
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
NB.   NB. A = Q * L
NB.   (}.   -: (ungql mp (trl~ -~/@$)@: }.   )@geql2) eA
NB.   NB. A = (Q) * L
NB.   (}.   -: (unmqlln (trlpick~ -~/@$)@  }.   )@geql2) eA
NB.
NB. Notes:
NB. - models LAPACK's xGEQL2
NB. - geql2 and geqlf are topologic equivalents
NB. - if L diagonal's non-negativity is required, then larfgb
NB.   would be replaced by larfpb

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
NB.   QfR   - (m+1)×n-matrix, Qf and R combined
NB.   Qf    - (m+1)×k-matrix, the unit lower trapezoidal
NB.           (unit diagonal not stored), represents the Q in
NB.           factored form
NB.   R     - k×n-matrix, the upper trapezoidal
NB.   Q     - m×m-matrix, the unitary (orthogonal), which is
NB.           defined as the product of k elementary
NB.           reflectors H(i) of order m:
NB.             Q = Π{H(i),i=0:k-1}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
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
NB.   NB. A = Q * R
NB.   (}:   -: (ungqr mp  tru        @: }:   )@geqr2) eA
NB.   NB. A = (Q) * R
NB.   (}:   -: (unmqrln  trupick        @  }:   )@geqr2) eA
NB.
NB. Notes:
NB. - models LAPACK's xGEQR2
NB. - geqr2 and geqrf are topologic equivalents
NB. - if R diagonal's non-negativity is required, then larfgf
NB.   would be replaced by larfpf

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
NB.   RQf   - m×(n+1)-matrix, R and Qf combined
NB.   R     - m×k-matrix, the upper trapezoidal
NB.   Qf    - k×(n+1)-matrix, the unit lower trapezoidal
NB.           (unit diagonal not stored), represents the Q in
NB.           factored form
NB.   Q     - n×n-matrix, the unitary (orthogonal), which is
NB.           defined as the product of k elementary
NB.           reflectors H(i) of order n:
NB.             Q = Π{H(i)',i=0:k-1}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
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
NB.   NB. A = R * Q
NB.   (}."1 -: ((tru~ -~/@$)@:(}."1) mp ungrq)@gerq2) eA
NB.   NB. A = R * (Q)
NB.   (}."1 -: (unmrqrn (trupick~ -~/@$)@:(}."1))@gerq2) eA
NB.
NB. Notes:
NB. - models LAPACK's xGERQ2
NB. - gerq2 and gerqf are topologic equivalents
NB. - if R diagonal's non-negativity is required, then
NB.   larfgbc would be replaced by larfpbc

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
NB.   iL    - m×m-matrix, the lower triangular
NB.   LZf   - m×(n+1)-matrix, the unit lower trapezoidal Zf
NB.           and the lower trapezoidal L combined:
NB.             LZf -: Tau ,. oA1 ,. A2 ,. oL
NB.   Tau   - m-vector, scalars τ[0:m-1] for Zf
NB.   oA1   - m×(l-1)-matrix, rows are vectors v[0:m-1] for
NB.           Zf
NB.   oL    - m×m-matrix, the lower triangular
NB.   Zf    - m×(n+1)-matrix, the Z represented in factored
NB.           form
NB.   Z     - n×n-matrix, the unitary (orthogonal), which is
NB.           defined as the product of m elementary
NB.           reflectors H(i) of order n:
NB.             Q = Π{H(i)',i=m-1:0}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   m   ≤ n
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
NB.   *                      - any scalar value, is not used
NB.   a1                     - elements of iA1
NB.   a2                     - elements of A2
NB.   il                     - elements of iL
NB.   ol                     - elements of oL
NB.   vi                     - l-vector v(i)
NB.   τi                     - scalar value conj(τ(i))
NB.   (vi,0,...,0,1,0,...,0) - n-vector u(i)
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
NB.   iL    - n×n-matrix, the lower triangular
NB.   A1    - (m-n-l+1)×n-matrix, is not changed
NB.   iA2   - (l-1)×n-matrix, part to be replaced by Zf
NB.   Trash - n-vector, will be replaced by Tau
NB.   ZfL   - (m+1)×n-matrix, the lower trapezoidal L and the
NB.           unit lower trapezoidal Zf combined:
NB.             ZfL -: oL , A1 , oA2 , Tau
NB.   oL    - n×n-matrix, the lower triangular
NB.   oA2   - (l-1)×n-matrix, columns are vectors v[0:n-1]
NB.           for Zf
NB.   Tau   - n-vector, scalars τ[0:n-1] for Zf
NB.   Zf    - (m+1)×n-matrix, the Z represented in factored
NB.           form
NB.   Z     - m×m-matrix, the unitary (orthogonal), which is
NB.           defined as the product of n elementary
NB.           reflectors H(i) of order m:
NB.             Z = Π{H(i),i=n-1:0}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   m   ≥ n
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
NB.   il                     - elements of iL
NB.   a1                     - elements of A1
NB.   a2                     - elements of iA2
NB.   *                      - any scalar value, is not used
NB.   ol                     - elements of oL
NB.   vi                     - l-vector v(i)
NB.   τi                     - scalar value τ(i)
NB.   (0,...,0,1,0,...,0,vi) - m-vector u(i)
NB.                            elementary re
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
NB.   iR    - n×n-matrix, the upper triangular
NB.   ZfR   - (m+1)×n-matrix, the upper trapezoidal R and the
NB.           unit upper trapezoidal Zf combined:
NB.             ZfR -: Tau , oA1 , A2 , oR
NB.   Tau   - n-vector, scalars τ[0:n-1] for Zf
NB.   oA1   - (l-1)×n-matrix, columns are vectors v[0:n-1]
NB.           for Zf
NB.   oR    - n×n-matrix, the upper triangular
NB.   Zf    - (m+1)×n-matrix, the Z represented in factored
NB.           form
NB.   Z     - m×m-matrix, the unitary (orthogonal), which is
NB.           defined as the product of n elementary
NB.           reflectors H(i) of order m:
NB.             Z = Π{H(i),i=0:n-1}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   m   ≥ n
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
NB.   *                      - any scalar value, is not used
NB.   a1                     - elements of iA1
NB.   a2                     - elements of A2
NB.   ir                     - elements of iR
NB.   or                     - elements of oR
NB.   vi                     - l-vector v(i)
NB.   τi                     - scalar value τ(i)
NB.   (vi,0,...,0,1,0,...,0) - m-vector u(i)
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
NB.   iR    - m×m-matrix, the upper triangular
NB.   A1    - m×(n-m-l+1)-matrix, is not changed
NB.   iA2   - m×(l-1)-matrix, part to be replaced by Zf
NB.   Trash - m-vector, will be replaced by Tau
NB.   RZf   - m×(n+1)-matrix, the upper trapezoidal R and the
NB.           unit upper trapezoidal Zf combined:
NB.             RZf -: oR ,. A1 ,. oA2 ,. Tau
NB.   oR    - m×m-matrix, the upper triangular
NB.   oA2   - m×(l-1)-matrix, rows are vectors v[0:m-1] for
NB.           Zf
NB.   Tau   - m-vector, scalars τ[0:m-1] for Zf
NB.   Zf    - m×(n+1)-matrix, the Z represented in factored
NB.           form
NB.   Z     - n×n-matrix, the unitary (orthogonal), which is
NB.           defined as the product of m elementary
NB.           reflectors H(i) of order n:
NB.             Z = Π{H(i)',i=0:m-1}
NB.             H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   m   ≤ n
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
NB.   ir                     - elements of iR
NB.   a1                     - elements of A1
NB.   a2                     - elements of iA2
NB.   *                      - any scalar value, is not used
NB.   or                     - elements of oR
NB.   vi                     - l-vector v(i)
NB.   τi                     - scalar value conj(τ(i))
NB.   (0,...,0,1,0,...,0,vi) - n-vector u(i)
NB.
NB. Notes:
NB. - models LAPACK's xLATRZ with the following difference:
NB.   - v(i) is stored instead of conj(v(i))
NB.   - conj(τ(i)) is stored instead of τ(i)
NB.   to keep consistency with gexxf
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
NB.   LQf - m×(n+1)-matrix, L and Qf combined
NB.   L   - m×k-matrix, the lower trapezoidal
NB.   Qf  - k×(n+1)-matrix, the unit upper trapezoidal (unit
NB.         diagonal not stored), represents the Q in
NB.         factored form
NB.   Q   - n×n-matrix, the unitary (orthogonal), which is
NB.         defined as the product of k elementary
NB.         reflectors H(i) of order n:
NB.           Q = Π{H(i)',i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
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
NB.   NB. A = L * Q
NB.   (     -: ( trl        @:(}:"1) mp unglq)@gelqf) A
NB.   NB. A = L * (Q)
NB.   (     -: (unmlqrn  trlpick        @:(}:"1))@gelqf) A
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
NB.   QfL - (m+1)×n-matrix, Qf and L combined
NB.   Qf  - (m+1)×k-matrix, the unit upper trapezoidal (unit
NB.         diagonal not stored), represents the Q in
NB.         factored form
NB.   L   - k×n-matrix, the lower trapezoidal
NB.   Q   - m×m-matrix, the unitary (orthogonal), which is
NB.         defined as the product of k elementary
NB.         reflectors H(i) of order m:
NB.           Q = Π{H(i),i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
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
NB.   NB. A = Q * L
NB.   (     -: (ungql mp (trl~ -~/@$)@: }.   )@geqlf) A
NB.   NB. A = (Q) * L
NB.   (     -: (unmqlln (trlpick~ -~/@$)@  }.   )@geqlf) A
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
NB.   QfR - (m+1)×n-matrix, Qf and R combined
NB.   Qf  - (m+1)×k-matrix, the unit lower trapezoidal (unit
NB.         diagonal not stored), represents the Q in
NB.         factored form
NB.   R   - k×n-matrix, the upper trapezoidal
NB.   Q   - m×m-matrix, the unitary (orthogonal), which is
NB.         defined as the product of k elementary reflectors
NB.         H(i) of order m:
NB.           Q = Π{H(i),i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
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
NB.   NB. A = Q * R
NB.   (     -: (ungqr mp  tru        @: }:   )@geqrf) A
NB.   NB. A = (Q) * R
NB.   (     -: (unmqrln  trupick        @  }:   )@geqrf) A
NB.
NB. Notes:
NB. - models LAPACK's xGEQRF
NB. - geqr2 and geqrf are topologic equivalents

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
NB.   RQf - m×(n+1)-matrix, R and Qf combined
NB.   R   - m×k-matrix, the upper trapezoidal
NB.   Qf  - k×(n+1)-matrix, the unit lower trapezoidal (unit
NB.         diagonal not stored), represents the Q in
NB.         factored form
NB.   Q   - n×n-matrix, the unitary (orthogonal), which is
NB.         defined as the product of k elementary reflectors
NB.         H(i) of order n:
NB.           Q = Π{H(i)',i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
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
NB.   NB. A = R * Q
NB.   (     -: ((tru~ -~/@$)@:(}."1) mp ungrq)@gerqf) A
NB.   NB. A = R * (Q)
NB.   (     -: (unmrqrn (trupick~ -~/@$)@:(}."1))@gerqf) A
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
NB.   iL  - m×m-matrix
NB.   LZf - m×(n+1)-matrix, the lower trapezoidal L and unit
NB.         lower trapezoidal Zf combined:
NB.           LZf -: Tau ,. oA0 ,. oL
NB.   Tau - m-vector, scalars τ[0:m-1] for Zf
NB.   oA0 - m×(n-m)-matrix, rows are vectors v[0:m-1] for Zf
NB.   oL  - m×m-matrix, the lower triangular
NB.   Zf  - m×(n+1)-matrix, the Z represented in factored
NB.         form
NB.   Z   - n×n-matrix, the unitary (orthogonal), which is
NB.         defined as the product of m elementary reflectors
NB.         H(i) of order n:
NB.           Z = Π{H(i)',i=m-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   m   ≤ n
NB.
NB. Storage layout for m=3, n=7:
NB.   A:
NB.     (     a0 a0 a0 a0 il *  *   )
NB.     (     a0 a0 a0 a0 il il *   )
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
NB.   a0                     - elements of iA0
NB.   il                     - elements of iL
NB.   ol                     - elements of oL
NB.   vi                     - l-vector v(i)
NB.   τi                     - scalar value conj(τ(i))
NB.   (vi,0,...,0,1,0,...,0) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = L * Z
NB.   (-: (({."1~ -@#) mp unglz)@tzlzf) A
NB.   NB. A = L * (Z)
NB.   (-: (unmlzrn ((1 -  c) {."1 ({."1~ -@#)))@tzlzf) A
NB.
NB. Notes:
NB. - latlz and tzlzf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0
NB. - strict upper triangle of iL is ignored

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
NB.   iL  - n×n-matrix
NB.   iA0 - (m-n)×n-matrix, part to be replaced by Zf
NB.   ZfL - (m+1)×n-matrix, the unit upper trapezoidal Zf and
NB.         lower trapezoidal L combined:
NB.           ZfL -: oL , oA0 , Tau
NB.   oL  - n×n-matrix, the lower triangular
NB.   oA0 - (m-n)×n-matrix, columns are vectors v[0:m-1] for
NB.         Zf
NB.   Tau - n-vector, scalars τ[0:n-1] for Zf
NB.   Zf  - (m+1)×n-matrix, the Z represented in factored
NB.         form
NB.   Z   - m×m-matrix, the unitary (orthogonal), which is
NB.         defined as the product of n elementary reflectors
NB.         H(i) of order m:
NB.           Z = Π{H(i),i=n-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   m   ≥ n
NB.
NB. Storage layout for m=7, n=3:
NB.   A:                ZfL:              Zf:               L:
NB.   (  il *  *   )    (  ol 0  0   )    (  1  0  0   )    (  ol 0  0   )
NB.   (  il il *   )    (  ol ol 0   )    (  0  1  0   )    (  ol ol 0   )
NB.   (  il il il  )    (  ol ol ol  )    (  0  0  1   )    (  ol ol ol  )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.                     (  τ0 τ1 τ2  )    (  τ0 τ1 τ2  )
NB. where
NB.   il                     - elements of iL
NB.   a0                     - elements of iA0
NB.   ol                     - elements of oL
NB.   vi                     - l-vector v(i)
NB.   τi                     - scalar value τ(i)
NB.   (0,...,0,1,0,...,0,vi) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = Z * L
NB.   (-: (ungzl mp ({.  ~   c))@tzzlf) A
NB.   NB. A = (Z) * L
NB.   (-: (unmzlln ((1 -~ #) {.   ({.  ~   c)))@tzzlf) A
NB.
NB. Notes:
NB. - latzl and tzzlf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0
NB. - strict upper triangle of iL is ignored

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
NB.   iR  - n×n-matrix
NB.   ZfR - (m+1)×n-matrix, the unit upper trapezoidal Zf and
NB.         upper trapezoidal R combined:
NB.           ZfR -: Tau , oA0 , oR
NB.   Tau - n-vector, scalars τ[0:n-1] for Zf
NB.   oA0 - (m-n)×n-matrix, columns are vectors v[0:m-1] for
NB.         Zf
NB.   oR  - n×n-matrix, the upper triangular
NB.   Zf  - (m+1)×n-matrix, the Z represented in factored
NB.         form
NB.   Z   - m×m-matrix, the unitary (orthogonal), which is
NB.         defined as the product of n elementary reflectors
NB.         H(i) of order m:
NB.           Z = Π{H(i),i=0:n-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   m   ≥ n
NB.
NB. Storage layout for m=7, n=3:
NB.   A:                ZfR:              Zf:               R:
NB.                     (  τ0 τ1 τ2  )    (  τ0 τ1 τ2  )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  a0 a0 a0  )    (  v0 v1 v2  )    (  v0 v1 v2  )    (  0  0  0   )
NB.   (  ir ir ir  )    (  or or or  )    (  1  0  0   )    (  or or or  )
NB.   (  *  ir ir  )    (  0  or or  )    (  0  1  0   )    (  0  or or  )
NB.   (  *  *  ir  )    (  0  0  or  )    (  0  0  1   )    (  0  0  or  )
NB. where
NB.   a0                     - elements of iA0
NB.   ir                     - elements of iR
NB.   or                     - elements of oR
NB.   vi                     - l-vector v(i)
NB.   τi                     - scalar value τ(i)
NB.   (vi,0,...,0,1,0,...,0) - m-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = Z * R
NB.   (-: (ungzr mp ({.  ~ -@c))@tzzrf) A
NB.   NB. A = (Z) * R
NB.   (-: (unmzrln ((1 -  #) {.   ({.  ~ -@c)))@tzzrf) A
NB.
NB. Notes:
NB. - latzr and tzzrf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0
NB. - strict lower triangle of iR is ignored

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
NB.   iR  - m×m-matrix
NB.   iA1 - m×(n-m)-matrix, part to be replaced by Zf
NB.   RZf - m×(n+1)-matrix, the upper trapezoidal R and unit
NB.         upper trapezoidal Zf combined:
NB.           RZf -: oR ,. oA1 ,. Tau
NB.   oR  - m×m-matrix, the upper triangular
NB.   oA1 - m×(n-m)-matrix, rows are vectors v[0:m-1] for Zf
NB.   Tau - m-vector, scalars τ[0:m-1] for Zf
NB.   Zf  - m×(n+1)-matrix, the Z represented in factored
NB.         form
NB.   Z   - n×n-matrix, the unitary (orthogonal), which is
NB.         defined as the product of m elementary reflectors
NB.         H(i) of order n:
NB.           Z = Π{H(i)',i=0:m-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   m   ≤ n
NB.
NB. Storage layout for m=3, n=7:
NB.   A:
NB.     (  ir ir ir a1 a1 a1 a1     )
NB.     (  *  ir ir a1 a1 a1 a1     )
NB.     (  *  *  ir a1 a1 a1 a1     )
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
NB.   ir                     - elements of iR
NB.   a1                     - elements of iA1
NB.   or                     - elements of oR
NB.   vi                     - l-vector v(i)
NB.   τi                     - scalar value conj(τ(i))
NB.   (0,...,0,1,0,...,0,vi) - n-vector u(i)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   NB. A = R * Z
NB.   (-: (({."1~   #) mp ungrz)@tzrzf) A
NB.   NB. A = R * (Z)
NB.   (-: (unmrzrn ((1 -~ c) {."1 ({."1~   #)))@tzrzf) A
NB.
NB. Notes:
NB. - models LAPACK's xTZRZF with the following differences:
NB.   - v(i) is stored instead of conj(v(i))
NB.   - conj(τ(i)) is stored instead of τ(i)
NB. - latrz and tzrzf are topologic equivalents
NB. - in u(i) 0s and 1 are not stored, v(i) is empty for l=0,
NB.   0s and 1 are absent and u(i) is empty when n=0
NB. - strict lower triangle of iR is ignored

tzrzf=: 3 : 0
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
NB.   Test:
NB.   - 128!:0 (built-in)
NB.   - qrd (math/misc addon)
NB.   - xGELQF xGEQLF xGEQRF xGERQF (math/lapack2 addon)
NB.   - gelqf geqlf geqrf gerqf (math/mt addon)
NB.   by general matrix
NB.
NB. Syntax:
NB.   testgeqf A
NB. where
NB.   A - m×n-matrix

testgeqf=: 3 : 0
  load        :: ] 'numeric'
  load_mttmp_ :: ] 'math/misc/mathutil'
  load_mttmp_ :: ] 'math/misc/makemat'
  load_mttmp_ :: ] 'math/misc/matutil'
  load_mttmp_ :: ] 'math/misc/linear'
  load_mttmp_ :: ] 'math/misc/matfacto'
  load_mttmp_ :: ] 'math/mt/test/lapack2/gelqf'
  load_mttmp_ :: ] 'math/mt/test/lapack2/geqlf'
  load_mttmp_ :: ] 'math/mt/test/lapack2/geqrf'
  load_mttmp_ :: ] 'math/mt/test/lapack2/gerqf'

  rcond=. (_."_)`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  norm=. norm1 y

  args=. y ; norm

  ('128!:0'        tmonad ((0&{::)`]                                                     `(rcond"_)`(_."_)`qrt01)) args
  ('qrd_mttmp_'    tmonad ((0&{::)`]                                                     `(rcond"_)`(_."_)`qrt01)) args

  ('dgelqf_mttmp_' tmonad ((0&{::)`( trl        @(0&{::) ;  unglq@(0&{:: stitcht  1&{::))`(rcond"_)`(_."_)`lqt01)) args
  ('dgeqlf_mttmp_' tmonad ((0&{::)`((trl~ -~/@$)@(0&{::) ;~ ungql@(0&{:: appendr~ 1&{::))`(rcond"_)`(_."_)`qlt01)) args
  ('dgeqrf_mttmp_' tmonad ((0&{::)`( tru        @(0&{::) ;~ ungqr@       ;              )`(rcond"_)`(_."_)`qrt01)) args
  ('dgerqf_mttmp_' tmonad ((0&{::)`((tru~ -~/@$)@(0&{::) ;  ungrq@(0&{:: stitchb~ 1&{::))`(rcond"_)`(_."_)`rqt01)) args

  ('zgelqf_mttmp_' tmonad ((0&{::)`( trl        @(0&{::) ;  unglq@(0&{:: stitcht  1&{::))`(rcond"_)`(_."_)`lqt01)) args
  ('zgeqlf_mttmp_' tmonad ((0&{::)`((trl~ -~/@$)@(0&{::) ;~ ungql@(0&{:: appendr~ 1&{::))`(rcond"_)`(_."_)`qlt01)) args
  ('zgeqrf_mttmp_' tmonad ((0&{::)`( tru        @(0&{::) ;~ ungqr@       ;              )`(rcond"_)`(_."_)`qrt01)) args
  ('zgerqf_mttmp_' tmonad ((0&{::)`((tru~ -~/@$)@(0&{::) ;  ungrq@(0&{:: stitchb~ 1&{::))`(rcond"_)`(_."_)`rqt01)) args

  ('gelqf'         tmonad ((0&{::)`( trl        @:(}:"1) ;  unglq                       )`(rcond"_)`(_."_)`lqt01)) args
  ('geqlf'         tmonad ((0&{::)`((trl~ -~/@$)@  }.    ;~ ungql                       )`(rcond"_)`(_."_)`qlt01)) args
  ('geqrf'         tmonad ((0&{::)`( tru        @  }:    ;~ ungqr                       )`(rcond"_)`(_."_)`qrt01)) args
  ('gerqf'         tmonad ((0&{::)`((tru~ -~/@$)@:(}."1) ;  ungrq                       )`(rcond"_)`(_."_)`rqt01)) args

  coerase < 'mttmp'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testtzqf
NB.
NB. Description:
NB.   Test:
NB.   - xTZRZF (math/lapack2 addon)
NB.   - tzlzf tzzlf tzzrf tzrzf (math/mt addon)
NB.   by trapezoidal matrix
NB.
NB. Syntax:
NB.   testtzqf A
NB. where
NB.   A - m×n-matrix
NB.
NB. TODO:
NB. - add xQRT12 test

testtzqf=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/tzrzf'

  rcond=. (_."_)`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  normw=. norm1 Awide=. |:^:(>/@$) y
  normt=. normi Atall=. |:^:(</@$) y

  NB. LAPACK doesn't clean strict lower triangle in R, so we need a rzt01 variant
  rzt01a=: ((1 {:: [) %~^:(0 < [) (norm1 % FP_EPS * 1 >. c)@((- trupick@(0&{::))~ (unmrzrn ((1 -~ c) {."1 trupick@({."1~ #)))))`0:@.(0 e. $@]) >. (norm1 % FP_EPS * 1 >. c)@(<: upddiag)@(unmrzrc ungrz)`0:@.(0 e. $)@]

  ('dtzrzf_mttmp_' tmonad ((0&{::)`(0&{:: ,.  1&{::)`(rcond"_)`(_."_)`rzt01a )) Awide ; normw
  ('ztzrzf_mttmp_' tmonad ((0&{::)`(0&{:: ,.  1&{::)`(rcond"_)`(_."_)`rzt01a )) Awide ; normw

  ('tzlzf'         tmonad ((0&{::)`]                `(rcond"_)`(_."_)`lzt01  )) Awide ; normw
  ('tzzlf'         tmonad ((0&{::)`]                `(rcond"_)`(_."_)`zlt01  )) Atall ; normt
  ('tzzrf'         tmonad ((0&{::)`]                `(rcond"_)`(_."_)`zrt01  )) Atall ; normt
  ('tzrzf'         tmonad ((0&{::)`]                `(rcond"_)`(_."_)`rzt01  )) Awide ; normw

  coerase < 'mttmp'
  erase 'rzt01a'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testqf
NB.
NB. Description:
NB.   Adv. to make verb to test gexxf by matrix of generator
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
