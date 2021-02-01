NB. Orthogonal factorization with pivoting
NB.
NB. gelpf     LQ factorization with row pivoting of a general
NB.           matrix
NB. geplf     QL factorization with column pivoting of a
NB.           general matrix
NB. geprf     QR factorization with column pivoting of a
NB.           general matrix
NB. gerpf     RQ factorization with row pivoting of a general
NB.           matrix
NB.
NB. testgepf  Test gexxf by general matrix
NB. testpf    Adv. to make verb to test gexxf by matrix of
NB.           generator and shape given
NB.
NB. Version: 0.12.0 2021-02-01
NB.
NB. Copyright 2011-2021 Igor Zhuravlov
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
NB. References:
NB. [1] C. H. Bischof, G. Quintana-Ortí. Computing Rank-
NB.     Revealing QR Factorizations of Dense Matrices. ACM
NB.     Transactions on Mathematical Software, June 1998,
NB.     Vol. 24, No. 2, pp. 226-253.
NB. [2] C. H. Bischof, G. Quintana-Ortí. Algorithm 782: Codes
NB.     for Rank-Revealing QR Factorizations of Dense
NB.     Matrices. ACM Transactions on Mathematical Software,
NB.     June 1998, Vol. 24, No. 2, pp. 254-257.
NB.
NB. TODO:
NB. - geplf: Q * L * P = A
NB. - gerpf: P * R * Q = A

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Constants

PFFP=: 0.9

PFSF=: 100  NB. Safe Factor, a security factor to avoid
            NB. arithmetic exceptions

NB. ---------------------------------------------------------
NB. lauc1f
NB.
NB. Description:
NB.   Apply incremental condition estimation to determine
NB.   whether the k-th column of A, stored in vector w, would
NB.   be acceptable as a pivot column with respect to the
NB.   threshold thresh
NB.
NB. Syntax:
NB.   out=. lauc1f ismin;iv;w;gamma;thresh
NB. where
NB.   ismin  - an estimate for the smallest singular value of
NB.            R(0:k-1,0:k-1)
NB.   iv     - k-vector, an approximate smallest left
NB.            singular vector of R(0:k-1,0:k-1)
NB.   w      - k-vector of elements R(0:k-1,k)
NB.   gamma  - k-th diagonal element R(k,k)
NB.   thresh - if osmin is smaller than thresh, the k-th
NB.            column is rejected
NB.   out    - is either 0 if k-th column of R is rejected,
NB.            or
NB.              out -: osmin;ov
NB.            if the k-th column of R is found acceptable
NB.   osmin  - an estimate for the smallest singular value of
NB.            R(0:k,0:k)
NB.   ov     - (k+1)-vector, an approximate smallest left
NB.            singular vector of R(0:k,0:k)
NB.   R      - m×n-matrix, the upper trapezoidal
NB.
NB. Notes:
NB. - models RRQR's xLAUC1 [1, 2]

lauc1f=: 3 : 0
  'smin v w gamma thresh'=. y
  if. thresh > | gamma do.
    0
  else.
    'smin cs'=. laic12 smin ; gamma , w mp + v
    if. thresh > smin do.
      0
    else.
      smin ; v ((* {:) , 0 { ]) cs
    end.
  end.
)

NB. ---------------------------------------------------------
NB. gelpc
NB.
NB. Description:
NB.   Continue a partial LQ factorization on
NB.   A(offset:offset+dsrd-1,:). If A(0:offset-1,:) has been
NB.   reduced to lower trapezoidal form, then gelpc applies
NB.   the traditional row pivoting strategy to identify dsrd
NB.   more independent rows of A with the restriction that
NB.   the condition number of the leading lower triangle of A
NB.   would not be larger than 1/rcond.
NB.   If lacptd (≤ dsrd) such rows are found, then the
NB.   condition number of lower triangle in
NB.   A(0:offset+lacptd-1,0:offset+lacptd-1) is less than
NB.   1/rcond. If lacptd < dsrd, then the LQ factorization
NB.   of A is completed, otherwise only dsrd new steps were
NB.   performed.
NB.
NB. Syntax:
NB.   'dLQf oA10 oA11 op ov osvlues omxnm'=. gelpc iA10;iA11;dsrd;rcond;ip;iv;isvlues;imxnm
NB. where
NB.   iA10    - (m-offset)×offset-matrix, left part of not yet
NB.             factorized rows
NB.   iA11    - (m-offset)×(n+1-offset)-matrix, right part of
NB.             not yet factorized rows, augmented by
NB.             trash vector
NB.   dsrd    > 0, number of independent rows to extract
NB.   rcond   > 0, 1/rcond specifies an upper bound of the
NB.             condition number
NB.   ip      - (m-offset)-vector, rows inversed permutation
NB.             of (iA10 ,. iA11)
NB.   iv      - offset-vector, an approximate smallest right
NB.             singular vector of iL
NB.   isvlues - 2-vector, estimates of the singular values:
NB.               sigma_max(A) , sigma_min(iL)
NB.   imxnm   - ≥ 0, the norm of largest row in iL
NB.   dLQf    - lacptd×(n+1)-matrix, identified more
NB.             independent rows of A,
NB.             (offset+lacptd)×(n+1)-matrix (LQf , dLQf)
NB.             contains oQ in factorized form and oL:
NB.               oL -: trl }:"1 LQf , dLQf
NB.               oQ -: unglq LQf , dLQf
NB.   oA10    - (m-offset-lacptd)×(offset+lacptd)-matrix,
NB.             left part of still not factorized rows
NB.   oA11    - (m-offset-lacptd)×(n+1-offset-lacptd)-matrix,
NB.             right part of still not factorized rows
NB.   op      - (m-offset)-vector, rows inversed
NB.             permutation of (dLQf , oA10 ,. oA11)
NB.   ov      - (offset+lacptd)-vector, an approximate
NB.             smallest right singular vector of oL
NB.   osvlues - 4-vector, estimates of the singular values:
NB.               sigma_max(A) , sigma_r(oL) , sigma_q(oL) , sigma_min(A)
NB.             where
NB.               r -: offset+lacptd
NB.               q -: </ m , n , r+1
NB.   omxnm   - ≥ 0, the norm of largest row in oL
NB.   LQf     - offset×(n+1)-matrix, already factorized
NB.             rows of A, contains iQ in factorized form
NB.             and iL:
NB.               iL -: trl }:"1 LQf
NB.               iQ -: unglq LQf
NB.   iL      -:trl (2 # offset) {. LQf
NB.   oL      -:trl (2 # offset+lacptd) {. LQf , dLQf
NB.   iQ      - offset×(n+1)-matrix with orthonormal rows,
NB.             which is defined as the product of (offset)
NB.             elementary reflectors
NB.   oQ      - (offset+lacptd)×(n+1)-matrix with orthonormal
NB.             rows, which is defined as the product of
NB.             (offset+lacptd) elementary reflectors
NB.   L       - m×k-matrix, the lower trapezoidal
NB.   A       - m×n-matrix to factorize
NB.   k       = min(m,n)
NB.
NB. Storage layout:
NB.   input:                         output:
NB.   (  LQf     LQf   ) offset      (  LQf            LQf   ) offset
NB.   (  iA10    iA11  ) m-offset    (  dLQf           dLQf  ) lacptd
NB.      offset  n+1-offset          (  oA10           oA11  ) m-offset-lacptd
NB.                                     offset+lacptd  n+1-offset-lacptd
NB.
NB. Notes:
NB. - mxnm is updated only when offset = 0

gelpc=: 3 : 0
  'A10 A11 dsrd rcond p v svlues mxnm'=. y
  n=. <: n1=. A10 +&c A11
  'm offset'=. +/\. $ A10
  k=. m <. n
  dLQf=. (0,n1)$0
  if. offset do.
    'smax smin'=. 0 1{svlues
  end.
  NB. Initialize partial row norms. The first item of
  NB. rnorms store the exact row norms
  rnorms=. ,:~ normsr }:"1 A11
  NB. Compute factorization
  i=. >: offset
  lasti=. k <. offset + dsrd
  while. i <: lasti do.
    NB. Determine ith pivot row and swap if necessary
    io=. liofmax {. rnorms
    if. io do.
      dip=. < 0 , io
      A10=. dip C. A10
      A11=. dip C. A11
      rnorms=. dip C."1 rnorms
      p=. ((# dLQf)&+&.:> dip) C. p
    end.
    rnorms=. }."1 rnorms
    w=. {. A10
    NB. Generate elementary reflector H(i)
    gamma=. {. z=. larfgfc {. A11
    NB. Apply elementary reflector H(I) to the corresponding
    NB. block of matrix A
    y=. (< (<0) ; 0) { A11=. (1 (0)} z) larfrnfr A11
    NB. Update partial row norms
    if. i < lasti do.
      'temp temp2'=. 2 %/\ (| y) , rnorms
      temp=. 0 >. 1 - *: temp
      temp2=. >: 20 %~ temp * *: temp2
      rnorms2=. }. normsr }:"1 A11
      rnorms3=. temp (((* %:)~ {.) 0} ]) rnorms
      rnorms2=. (,:~ 1 = temp2)} rnorms3 ,: rnorms2
      rnorms=. (,:~ 0 = {. rnorms)} rnorms2 ,: rnorms
    end.
    NB. Check new row for independence
    if. 1 = i do.
      mxnm=. smin=. smax=. | gamma
      v=. 1
      if. 0 = mxnm do.
        svlues=. smin 2} svlues
        A11=. z (< 0 ; <<0)} A11
        break.
      end.
    else.
      smaxpr=. mxnm * 3 %: i
      out=. lauc1f smin ; v ; w ; gamma ; smaxpr * rcond
      if. 0 = L. out do.
        NB. Row rejected
        A11=. z 0} A11
        break.
      end.
      NB. Row accepted
      'smin v'=. out
      smax=. smaxpr
    end.
    dLQf=. dLQf , w , z
    A10=. (}. A10) ,. y
    A11=. 1 1 }. A11
    i=. >: i
  end.
  svlues=. (smax,smin) 0 1} svlues
  if. dsrd = # dLQf do.
    NB. DSRD independent rows have been found
    svlues=. smin 2 3} svlues
  else.
    NB. All remaining rows rejected
    if. <./ $ A11 do.
      NB. Factor remaining columns
      A11=. gelqf }:"1 A11   NB. FIXME TWICE excessive re-shaping
    end.
    eAsfx=. A10 ,. A11
    NB. Use incremental condition estimation to get an
    NB. estimate of the smallest singular value
    j=. 0 = c A10
    jsupr=. k - c A10
    while. j < jsupr do.
      'w gamma'=. (}: ; {:) (< j ; (i. >: # v)) { eAsfx
      'smin cs'=. laic12 smin ; gamma , w mp + v
      v=. v ((* {:) , 0 { ]) cs
      if. 0 = j do.
        svlues=. smin 2} svlues
      end.
      j=. >: j
    end.
    svlues=. smin 3} svlues
  end.
  dLQf ; A10 ; A11 ; p ; v ; svlues ; mxnm
)

NB. ---------------------------------------------------------
NB. geprc
NB.
NB. Description:
NB.   Continue a partial QR factorization on
NB.   A(:,offset:offset+dsrd-1). If A(:,0:offset-1) has been
NB.   reduced to upper trapezoidal form, then geprc applies
NB.   the traditional column pivoting strategy to identify
NB.   dsrd more independent columns of A with the restriction
NB.   that the condition number of the leading upper triangle
NB.   of A would not be larger than 1/rcond.
NB.   If lacptd (≤ dsrd) such columns are found, then the
NB.   condition number of upper triangle in
NB.   A(0:offset+lacptd-1,0:offset+lacptd-1) is less than
NB.   1/rcond. If lacptd < dsrd, then the QR factorization
NB.   of A is completed, otherwise only dsrd new steps were
NB.   performed.
NB.
NB. Syntax:
NB.   'dQfR oA01 oA11 op ov osvlues omxnm'=. geprc iA01;iA11;dsrd;rcond;ip;iv;isvlues;imxnm
NB. where
NB.   iA01    - offset×(n-offset)-matrix, top part of not yet
NB.             factorized columns
NB.   iA11    - (m+1-offset)×(n-offset)-matrix, bottom part
NB.             of not yet factorized columns, augmented by
NB.             trash vector
NB.   dsrd    > 0, number of independent columns to extract
NB.   rcond   > 0, 1/rcond specifies an upper bound of the
NB.             condition number
NB.   ip      - (n-offset)-vector, columns inversed
NB.             permutation of (iA01 , iA11)
NB.   iv      - offset-vector, an approximate smallest left
NB.             singular vector of iR
NB.   isvlues - 2-vector, estimates of the singular values:
NB.               sigma_max(A) , sigma_min(iR)
NB.   imxnm   - ≥ 0, the norm of largest column in iR
NB.   dQfR    - (m+1)×lacptd-matrix, identified more
NB.             independent columns of A,
NB.             (m+1)×(offset+lacptd)-matrix (QfR ,. dQfR)
NB.             contains oQ in factorized form and oR:
NB.               oR -: tru }: QfR ,. dQfR
NB.               oQ -: ungqr QfR ,. dQfR
NB.   oA01    - (offset+lacptd)×(n-offset-lacptd)-matrix, top
NB.             part of still not factorized columns
NB.   oA11    - (m+1-offset-lacptd)×(n-offset-lacptd)-matrix,
NB.             bottom part of still not factorized columns
NB.   op      - (n-offset)-vector, columns inversed
NB.             permutation of (dQfR ,. oA01 , oA11)
NB.   ov      - (offset+lacptd)-vector, an approximate
NB.             smallest left singular vector of oR
NB.   osvlues - 4-vector, estimates of the singular values:
NB.               sigma_max(A) , sigma_r(oR) , sigma_q(oR) , sigma_min(A)
NB.             where
NB.               r -: offset+lacptd
NB.               q -: </ m , n , r+1
NB.   omxnm   - ≥ 0, the norm of largest column in oR
NB.   QfR     - (m+1)×offset-matrix, already factorized
NB.             columns of A, contains iQ in factorized form
NB.             and iR:
NB.               iR -: tru }: QfR
NB.               iQ -: ungqr QfR
NB.   iQ      - (m+1)×offset-matrix with orthonormal columns,
NB.             which is defined as the product of (offset)
NB.             elementary reflectors
NB.   oQ      - (m+1)×(offset+lacptd)-matrix with orthonormal
NB.             columns, which is defined as the product of
NB.             (offset+lacptd) elementary reflectors
NB.   iR      -:tru (2 # offset) {. QfR
NB.   oR      -:tru (2 # offset+lacptd) {. QfR ,. dQfR
NB.   R       - k×n-matrix, the upper trapezoidal
NB.   A       - m×n-matrix to factorize
NB.   k       = min(m,n)
NB.
NB. Storage layout:
NB.   input:                           output:
NB.   (  QfR     iA01  ) offset        (  QfR     dQfR    oA01  ) offset+lacptd
NB.   (  QfR     iA11  ) m+1-offset    (  QfR     dQfR    oA11  ) m+1-offset-lacptd
NB.      offset  n-offset                 offset  lacptd  n-offset-lacptd
NB.
NB. Notes:
NB. - models RRQR's xGEQPC [1, 2] with following differences:
NB.   - Q returned is represented as the product of
NB.     elementary reflectors, and is merged with R
NB. - mxnm is updated only when offset = 0

geprc=: 3 : 0
  'A01 A11 dsrd rcond p v svlues mxnm'=. y
  m=. <: m1=. A01 +&# A11
  'offset n'=. +/\ $ A01
  k=. m <. n
  dQfR=. (m1,0)$0
  if. offset do.
    'smax smin'=. 0 1{svlues
  end.
  NB. Initialize partial column norms. The first item of
  NB. cnorms store the exact column norms
  cnorms=. ,:~ normsc }: A11
  NB. Compute factorization
  i=. >: offset
  lasti=. k <. offset + dsrd
  while. i <: lasti do.
    NB. Determine ith pivot column and swap if necessary
    io=. liofmax {. cnorms
    if. io do.
      dip=. < 0 , io
      A01=. dip C."1 A01
      A11=. dip C."1 A11
      cnorms=. dip C."1 cnorms
      p=. ((c dQfR)&+&.:> dip) C. p
    end.
    cnorms=. }."1 cnorms
    w=. {."1 A01
    NB. Generate elementary reflector H(i)
    gamma=. {. z=. larfgf {."1 A11
    NB. Apply elementary reflector H(I) to the corresponding
    NB. block of matrix A
    y=. (< 0 ; <<0) { A11=. (1 (0)} z) larflcfc A11
    NB. Update partial column norms
    if. i < lasti do.
      'temp temp2'=. 2 %/\ (| y) , cnorms
      temp=. 0 >. 1 - *: temp
      temp2=. >: 20 %~ temp * *: temp2
      cnorms2=. }. normsc }: A11
      cnorms3=. temp (((* %:)~ {.) 0} ]) cnorms
      cnorms2=. (,:~ 1 = temp2)} cnorms3 ,: cnorms2
      cnorms=. (,:~ 0 = {. cnorms)} cnorms2 ,: cnorms
    end.
    NB. Check new column for independence
    if. 1 = i do.
      mxnm=. smin=. smax=. | gamma
      v=. 1
      if. 0 = mxnm do.
        svlues=. smin 2} svlues
        A11=. z (< (<0) ; 0)} A11
        break.
      end.
    else.
      smaxpr=. mxnm * 3 %: i
      out=. lauc1f smin ; v ; w ; gamma ; smaxpr * rcond
      if. 0 = L. out do.
        NB. Column rejected
        A11=. z (< a: ; 0)} A11
        break.
      end.
      NB. Column accepted
      'smin v'=. out
      smax=. smaxpr
    end.
    dQfR=. dQfR ,. w , z
    A01=. (}."1 A01) , y
    A11=. 1 1 }. A11
    i=. >: i
  end.
  svlues=. (smax,smin) 0 1} svlues
  if. dsrd = c dQfR do.
    NB. DSRD independent columns have been found
    svlues=. smin 2 3} svlues
  else.
    NB. All remaining columns rejected
    if. <./ $ A11 do.
      NB. Factor remaining columns
      A11=. geqrf }: A11   NB. FIXME TWICE excessive re-shaping
    end.
    eAsfx=. A01 , A11
    NB. Use incremental condition estimation to get an
    NB. estimate of the smallest singular value
    j=. 0 = # A01
    jsupr=. k - # A01
    while. j < jsupr do.
      'w gamma'=. (}: ; {:) (< (i. >: # v) ; j) { eAsfx
      'smin cs'=. laic12 smin ; gamma , w mp + v
      v=. v ((* {:) , 0 { ]) cs
      if. 0 = j do.
        svlues=. smin 2} svlues
      end.
      j=. >: j
    end.
    svlues=. smin 3} svlues
  end.
  dQfR ; A01 ; A11 ; p ; v ; svlues ; mxnm
)

NB. ---------------------------------------------------------
NB. gelpw
NB.
NB. Description:
NB.   Apply one block step of the Householder LQ
NB.   factorization algorithm with restricted pivoting. It is
NB.   called by gelpb to factorize a window of the matrix.
NB.   Let A be the partial LQ factorization of
NB.   (offset+lwsize)×n-matrix C, i.e. we have computed an
NB.   orthogonal matrix iQ and a permutation matrix iP such
NB.   that
NB.     iP * iL * iQ = C
NB.   Then in addition let iv be an approximate smallest
NB.   right singular vector of iL in the sense that
NB.     sigma_min(iL) ~ twonorm(iL * iv) = smin
NB.   and
NB.     sigma_max(iL) ~ (offset ^ 1r3) * mxnm = smax
NB.   with
NB.     cond_no(iL) ~ smax/smin ≤ 1/rcond.
NB.   Then gelpw tries to identify nb rows in (iA10,.iA11)
NB.   such that
NB.     cond_no(oL) < 1/rcond.
NB.   On exit,
NB.     oP * oL * oQ = C
NB.   is again a partial LQ factorization of C, but last
NB.   lacptd rows of oL have been reduced via a series of
NB.   elementary reflectors. Further
NB.     sigma_min(oL) ~ twonorm(oL * ov) = smin
NB.   and
NB.     sigma_max(oL) ~ ((offset+lacptd) ^ 1r3) * mxnm = smax
NB.   with
NB.     cond_no(oL) ~ smax/smin ≤ 1/rcond.
NB.   In the ideal case, lacptd = nb, that is, we found nb
NB.   independent rows in the window consisting of the
NB.   first lwsize rows of A.
NB.
NB. Syntax:
NB.   'dLQf oA10 oA11 op osmin ov'=. gelpw iA10;iA11;nb;rcond;ip;mxnm;ismin;iv
NB. where
NB.   iA10   - lwsize×offset-matrix, left part of not yet
NB.            factorized part of eA
NB.   iA11   - lwsize×(n+1-offset)-matrix, right part of not
NB.            yet factorized part of eA
NB.   nb     > 0, the number of independent rows to identify.
NB.            This equals the desired blocksize in gelpb
NB.   rcond  > 0, 1/rcond specifies an upper bound on the
NB.            condition number
NB.   ip     - lwsize-vector, rows inversed permutation of
NB.            A, represents permutation m×m-matrix
NB.   mxnm   ≥ 0, the norm of the largest rown in iL
NB.   ismin  ≥ 0, an estimate of smallest singular value of
NB.            iL
NB.   iv     - offset-vector, an approximate right
NB.            null-vector of iL
NB.   dLQf   - lacptd×(n+1)-matrix, factorized rows,
NB.            (offset+lacptd)×(n+1)-matrix (LQf , dLQf)
NB.            contains oQ in factorized form and oL:
NB.              oL -: trl }:"1 LQf , dLQf
NB.              oQ -: unglq LQf , dLQf
NB.   oA10   - (lwsize-lacptd)×(offset+lacptd)-matrix, left
NB.            part of not yet factorized or rejected rows
NB.   oA11   - (lwsize-lacptd)×(n+1-offset-lacptd)-matrix,
NB.            right part of not yet factorized or rejected
NB.            rows, augmented by trash vector
NB.   op     - lwsize-vector, rows inversed permutation of
NB.            A, represents permutation m×m-matrix
NB.   osmin  ≥ 0, an estimate of smallest singular value of
NB.            oL
NB.   ov     - (offset+lacptd)-vector, an approximate right
NB.            null-vector of oL
NB.   LQf    - offset×(n+1)-matrix, already factorized
NB.            rows of A, contains iQ in factorized form
NB.            and iL:
NB.              iL -: trl }:"1 LQf
NB.              iQ -: unglq LQf
NB.   iL     -:trl (2 # offset) {. LQf
NB.   oL     -:trl (2 # offset+lacptd) {. LQf , dLQf
NB.   iQ     - offset×(n+1)-matrix with orthonormal rows,
NB.            which is defined as the product of (offset)
NB.            elementary reflectors
NB.   oQ     - (offset+lacptd)×(n+1)-matrix with orthonormal
NB.            rows, which is defined as the product of
NB.            (offset+lacptd) elementary reflectors
NB.   lacptd ≤ nb
NB.   lwsize ≤ n-offset, the size of pivot window
NB.   eA     - m×(n+1)-matrix, being A, augmented by zero
NB.            vector
NB.              eA -: A ,. 0
NB.
NB. Storage layout:
NB.   input:                       output:
NB.   (  LQf     LQf   ) offset    (  LQf            LQf   ) offset
NB.   (  iA10    iA11  ) lwsize    (  dLQf           dLQf  ) lacptd
NB.      offset  n+1-offset        (  oA10           0A11  ) lwsize-lacptd
NB.                                   offset+lacptd  n+1-offset-lacptd

gelpw=: 3 : 0
  'A10 A11 nb rcond p mxnm smin v'=. y
  offset=. c A10
  n=. <: n1=. offset + c A11
  dLQf=. (0,n1)$0
  NB. Initialize partial row norms (stored in the first item
  NB. of rnorms) and exact row norms (stored in the second
  NB. item of rnorms) for the first batch of rows
  rnorms=. ,:~ normsr }:"1 A11
  NB. Main loop
  lastk=. <: (n - offset) <. # p
  while. nb > # dLQf do.
    NB. Determine pivot candidate
    io=. liofmax {. rnorms
    if. io do.
      dip=. < 0 , io
      A10=. dip C. A10
      A11=. dip C. A11
      rnorms=. dip C."1 rnorms
      p=. ((# dLQf)&+&.:> dip) C. p
    end.
    NB. Determine (offset+lacptd)st diagonal element gamma
    NB. of matrix A would elementary reflector be applied
    gamma=. rnorms (negpos~ 9&o.) & (0&({,)) A11
    NB. Update estimate for largest singular value
    smax=. mxnm * 3 %: >: c A10
    u=. {."1 A11
    w=. {. A10
    NB. Is candidate pivot row acceptable ?
    out=. lauc1f smin ; v ; w ; gamma ; smax * rcond
    if. 0 = L. out do. break. end.
    NB. Pivot candidate was accepted
    'smin v'=. out
    NB. Generate Householder vector
    z=. larfgfc {. A11
    dLQf=. dLQf , w , z
    NB. Apply Householder reflection to A11
    y=. (< (<0) ; 0) { A11=. (1 (0)} z) larfrnfr A11
    A10=. (}. A10) ,. y
    A11=. 1 1 }. A11
    rnorms=. }."1 rnorms
    NB. Update partial column norms
    if. lastk > # dLQf do.
      'temp temp2'=. 2 %/\ (| y) , rnorms
      temp=. 0 >. 1 - *: temp
      temp2=. >: 20 %~ temp * *: temp2
      rnorms2=. }. normsr }:"1 A11
      rnorms3=. temp (((* %:)~ {.) 0} ]) rnorms
      rnorms2=. (,:~ 1 = temp2)} rnorms3 ,: rnorms2
      rnorms=. (,:~ 0 = {. rnorms)} rnorms2 ,: rnorms
    end.
  end.
  NB. Reject all remaining rows in pivot window
  dLQf ; A10 ; A11 ; p ; smin ; v
)

NB. ---------------------------------------------------------
NB. geprw
NB.
NB. Description:
NB.   Apply one block step of the Householder QR
NB.   factorization algorithm with restricted pivoting. It is
NB.   called by geprb to factorize a window of the matrix.
NB.   Let A be the partial QR factorization of
NB.   m×(offset+lwsize)-matrix C, i.e. we have computed an
NB.   orthogonal matrix iQ and a permutation matrix iP such
NB.   that
NB.     iQ * iR * iP = C
NB.   Then in addition let iv be an approximate smallest left
NB.   singular vector of iR in the sense that
NB.     sigma_min(iR) ~ twonorm(iR' * iv) = smin
NB.   and
NB.     sigma_max(iR) ~ (offset ^ 1r3) * mxnm = smax
NB.   with
NB.     cond_no(iR) ~ smax/smin ≤ 1/rcond.
NB.   Then geprw tries to identify nb columns in (iA01,iA11)
NB.   such that
NB.     cond_no(oR) < 1/rcond.
NB.   On exit,
NB.     oQ * oR * oP = C
NB.   is again a partial QR factorization of C, but last
NB.   lacptd columns of oR have been reduced via a series of
NB.   elementary reflectors. Further
NB.     sigma_min(oR) ~ twonorm(oR' * ov) = smin
NB.   and
NB.     sigma_max(oR) ~ ((offset+lacptd) ^ 1r3) * mxnm = smax
NB.   with
NB.     cond_no(oR) ~ smax/smin ≤ 1/rcond.
NB.   In the ideal case, lacptd = nb, that is, we found nb
NB.   independent columns in the window consisting of the
NB.   first lwsize columns of A.
NB.
NB. Syntax:
NB.   'dQfR oA01 oA11 op osmin ov'=. geprw iA01;iA11;nb;rcond;ip;mxnm;ismin;iv
NB. where
NB.   iA01   - offset×lwsize-matrix, top part of not yet
NB.            factorized part of eA
NB.   iA11   - (m+1-offset)×lwsize-matrix, bottom part of not
NB.            yet factorized part of eA
NB.   nb     > 0, the number of independent columns to
NB.            identify. This equals the desired blocksize in
NB.            geprb
NB.   rcond  > 0, 1/rcond specifies an upper bound on the
NB.            condition number
NB.   ip     - lwsize-vector, columns inversed permutation of
NB.            A, represents permutation n×n-matrix
NB.   mxnm   ≥ 0, the norm of the largest column in iR
NB.   ismin  ≥ 0, an estimate of smallest singular value of
NB.            iR
NB.   iv     - offset-vector, an approximate left null-vector
NB.            of iR
NB.   dQfR   - (m+1)×lacptd-matrix, factorized columns,
NB.            (m+1)×(offset+lacptd)-matrix (QfR ,. dQfR)
NB.            contains oQ in factorized form and oR:
NB.              oR -: tru }: QfR ,. dQfR
NB.              oQ -: ungqr QfR ,. dQfR
NB.   oA01   - (offset+lacptd)×(lwsize-lacptd)-matrix, top
NB.            part of not yet factorized or rejected columns
NB.   oA11   - (m+1-offset-lacptd)×(lwsize-lacptd)-matrix,
NB.            bottom part of not yet factorized or rejected
NB.            columns, augmented by trash vector
NB.   op     - lwsize-vector, columns inversed permutation of
NB.            A, represents permutation n×n-matrix
NB.   osmin  ≥ 0, an estimate of smallest singular value of
NB.            oR
NB.   ov     - (offset+lacptd)-vector, an approximate left
NB.            null-vector of oR
NB.   QfR    - (m+1)×offset-matrix, already factorized
NB.            columns of A, contains iQ in factorized form
NB.            and iR:
NB.              iR -: tru }: QfR
NB.              iQ -: ungqr QfR
NB.   iQ     - (m+1)×offset-matrix with orthonormal columns,
NB.            which is defined as the product of (offset)
NB.            elementary reflectors
NB.   oQ     - (m+1)×(offset+lacptd)-matrix with orthonormal
NB.            columns, which is defined as the product of
NB.            (offset+lacptd) elementary reflectors
NB.   iR      -:tru (2 # offset) {. QfR
NB.   oR      -:tru (2 # offset+lacptd) {. QfR ,. dQfR
NB.   lacptd ≤ nb
NB.   lwsize ≤ n-offset, the size of pivot window
NB.   eA     - (m+1)×n-matrix, being A, augmented by zero
NB.            vector
NB.              eA -: A , 0
NB.
NB. Storage layout:
NB.   input:                           output:
NB.   (  QfR     iA01  ) offset        (  QfR     dQfR    oA01  ) offset+lacptd
NB.   (  QfR     iA11  ) m+1-offset    (  QfR     dQfR    oA11  ) m+1-offset-lacptd
NB.      offset  lwsize                   offset  lacptd  lwsize-lacptd
NB.
NB. Notes:
NB. - models RRQR's xGEQPW [1, 2]

geprw=: 3 : 0
  'A01 A11 nb rcond p mxnm smin v'=. y
  offset=. # A01
  m=. <: m1=. offset + # A11
  dQfR=. (m1,0)$0
  NB. Initialize partial column norms (stored in the first
  NB. item of cnorms) and exact column norms (stored in the
  NB. second item of cnorms) for the first batch of columns
  cnorms=. ,:~ normsc }: A11
  NB. Main loop
  lastk=. <: (m - offset) <. # p
  while. nb > c dQfR do.
    NB. Determine pivot candidate
    io=. liofmax {. cnorms
    if. io do.
      dip=. < 0 , io
      A01=. dip C."1 A01
      A11=. dip C."1 A11
      cnorms=. dip C."1 cnorms
      p=. ((c dQfR)&+&.:> dip) C. p
    end.
    NB. Determine (offset+lacptd)st diagonal element gamma
    NB. of matrix A would elementary reflector be applied
    gamma=. cnorms (negpos~ 9&o.) & (0&({,)) A11
    NB. Update estimate for largest singular value
    smax=. mxnm * 3 %: >: # A01
    u=. {. A11
    w=. {."1 A01
    NB. Is candidate pivot column acceptable ?
    out=. lauc1f smin ; v ; w ; gamma ; smax * rcond
    if. 0 = L. out do. break. end.
    NB. Pivot candidate was accepted
    'smin v'=. out
    NB. Generate Householder vector
    z=. larfgf {."1 A11
    dQfR=. dQfR ,. w , z
    NB. Apply Householder reflection to A11
    y=. (< 0 ; <<0) { A11=. (1 (0)} z) larflcfc A11
    A01=. (}."1 A01) , y
    A11=. 1 1 }. A11
    cnorms=. }."1 cnorms
    NB. Update partial column norms
    if. lastk > c dQfR do.
      'temp temp2'=. 2 %/\ (| y) , cnorms
      temp=. 0 >. 1 - *: temp
      temp2=. >: 20 %~ temp * *: temp2
      cnorms2=. }. normsc }: A11
      cnorms3=. temp (((* %:)~ {.) 0} ]) cnorms
      cnorms2=. (,:~ 1 = temp2)} cnorms3 ,: cnorms2
      cnorms=. (,:~ 0 = {. cnorms)} cnorms2 ,: cnorms
    end.
  end.
  NB. Reject all remaining columns in pivot window
  dQfR ; A01 ; A11 ; p ; smin ; v
)

NB. ---------------------------------------------------------
NB. gelpb
NB.
NB. Description:
NB.   Compute a LQ factorization:
NB.     P * (  L00   0    ) * Q = A
NB.         (  L10   L11  )
NB.   of matrix A. The permutation P is chosen with the goal
NB.   to reveal the rank of A by a suitably dimensioned
NB.   trailing submatrix L11 with norm(L11) being small.
NB.
NB. Syntax:
NB.   'ip L Qn orcond rank svlues'=. ircond gelpb A
NB. where
NB.   A      - m×n-matrix, the input to factorize
NB.   ircond > 0, 1/ircond specifies an upper bound on the
NB.            condition number of L00
NB.   ip     - m-vector, rows inversed permutation of A,
NB.            represents permutation m×m-matrix P
NB.   L      - m×k-matrix, the lower trapezoidal
NB.   Qn     - n×n-matrix, the unitary (orthogonal), embeds
NB.            Q:
NB.              Q -: k {. Qn
NB.   orcond > 0, 1/orcond is an estimate for the condition
NB.            number of L00
NB.   rank   ≥ 0, an estimate for the numerical rank of A
NB.            with respect to the threshold 1/ircond in the
NB.            sense that
NB.              rank = arg_max(cond_no(L(0:r-1,0:r-1))<1/ircond)
NB.            This may be an underestimate of the rank if
NB.            the leading rows were not well-conditioned
NB.   svlues - 4-vector, estimates of the singular values,
NB.            see gelpf
NB.   Q      - k×n-matrix with orthonormal rows, which is
NB.            defined as the first k rows of the product of
NB.            k elementary reflectors
NB.   k      = min(m,n)
NB.
NB. Storage layout:
NB.   L  =  (  L00   0    ) rank
NB.         (  L10   L11  ) m-rank
NB.            rank  n-rank

gelpb=: 4 : 0
  k=. <./ 'm n'=. $ y
  p=. i. m
  if. 0 = k do.
    p ; ((m,0)$0) ; ((2#n)$0) ; 0 ; 0 ; 4$0
    return.
  end.
  y=. y ,. 0
  NB. Move row with largest residual norm left front
  'y A10 A11 p v svlues mxnm'=. gelpc ((m,0)$0);y;1;x;p;(i.0);(4$_.);_.
  if. # y do.
    if. 1 = k do.
      p ; (trl }:"1 y) ; (n unglq y) ; (%/ 1 0 { svlues) ; 1 ; svlues
      return.
    else.
      smin=. 1 { svlues
    end.
  else.
    p ; (trl }:"1 A11) ; (n unglq A11) ; 0 ; 0 ; 4$0
    return.
  end.
  NB. Factor remaining rows using blocked code with
  NB. restricted pivoting strategy
  nb=. QFNB <. k
  NB. The size of the pivot window is chosen to be nb+nllity
  nllity=. k <. 10 >. <. 0.5 0.05 mp nb , m
  norej=. m
  while. (# y) < k <. norej do.
    NB. Invariant: the 1st row of A1x is the 1st row
    NB. in currently considered block row
    kb=. <./ nb , (k , norej) - # y
    NB. The goal now is to find kb independent rows
    NB. among the remaining norej not yet rejected rows
    lwsize=. (norej - # y) <. kb + nllity
    'A00 A10'=. lwsize ({. ; }.) A10
    'A01 A11'=. lwsize ({. ; }.) A11
    pwi=. (i. lwsize) + # y
    'dLQf A00 A01 pw smin v'=. gelpw A00;A01;kb;x;(pwi{p);mxnm;smin;v
    p=. pw pwi} p
    if. # dLQf do.
      NB. Accumulate Householder vectors in a block reflector
      if. # A11 do.
        NB. Apply block reflector to A11
        A11=. ((# y) tru1 dLQf) larfbrnfr A11
      end.
      y=. y , dLQf
    end.
    NB. Move rejected rows to the end if there is space
    if. kb > # dLQf do.
      kklwsz=. y +&# A01
      if. norej > kklwsz do.
        'ppfx psfx'=. (# y) ({. ; }.) p
        p=. ppfx , (# A01) |. psfx
        A10=. (A10 ,. (# dLQf) {."1 A11) , A00
        A11=. ((# dLQf) }."1 A11) , A01
        norej=. norej - # A01
      else.
        A10=. A00 , A10 ,. (# dLQf) {."1 A11
        A11=. A01 , (# dLQf) }."1 A11
        break.
      end.
    else.
      A10=. A00 , A10 ,. (# dLQf) {."1 A11
      A11=. A01 , (# dLQf) }."1 A11
    end.
  end.
  svlues=. (smin , mxnm * 3 %: # y) 1 0} svlues
  if. k > # y do.
    NB. Process rejected rows
    'ppfx psfx'=. (# y) ({. ; }.) p
    'eApfx A10 A11 psfx v svlues mxnm'=. gelpc A10;A11;(# A11);x;psfx;v;svlues;mxnm
    p=. ppfx , psfx
    y=. y , eApfx
  else.
    svlues=. smin 2 3} svlues
  end.
  rank=. # y
  y=. y , A10 ,. A11
  p ; (trl }:"1 y) ; (n unglq y) ; (%/ 1 0 { svlues) ; rank ; svlues
)

NB. ---------------------------------------------------------
NB. geprb
NB.
NB. Description:
NB.   Compute a QR factorization:
NB.     Q * (  R00  R01  ) * P = A
NB.         (  0    R11  )
NB.   of matrix A. The permutation P is chosen with the goal
NB.   to reveal the rank of A by a suitably dimensioned
NB.   trailing submatrix R11 with norm(R11) being small.
NB.
NB. Syntax:
NB.   'Qm R ip orcond rank svlues'=. ircond geprb A
NB. where
NB.   A      - m×n-matrix, the input to factorize
NB.   ircond > 0, 1/ircond specifies an upper bound on the
NB.            condition number of R00
NB.   Qm     - m×m-matrix, the unitary (orthogonal), embeds
NB.            Q:
NB.              Q -: k {."1 Qm
NB.   R      - k×n-matrix, the upper trapezoidal
NB.   ip     - n-vector, columns inversed permutation of A,
NB.            represents permutation n×n-matrix P
NB.   orcond > 0, 1/orcond is an estimate for the condition
NB.            number of R00
NB.   rank   ≥ 0, an estimate for the numerical rank of A
NB.            with respect to the threshold 1/ircond in the
NB.            sense that
NB.              rank = arg_max(cond_no(R(0:r-1,0:r-1))<1/ircond)
NB.            This may be an underestimate of the rank if
NB.            the leading columns were not well-conditioned
NB.   svlues - 4-vector, estimates of the singular values,
NB.            see geprf
NB.   Q      - m×k-matrix with orthonormal columns, which is
NB.            defined as the first k columns of the product
NB.            of k elementary reflectors
NB.   k      = min(m,n)
NB.
NB. Storage layout:
NB.   R  =  (  R00   R01  ) rank
NB.         (  0     R11  ) m-rank
NB.            rank  n-rank
NB.
NB. Notes:
NB. - models RRQR's xGEQPB(3) [1, 2] with following
NB.   differences:
NB.   - matrix C is an identity matrix
NB.   - rejected columns are moved to R11 not from the left
NB.     side, but from the right
NB.   - rejected columns are moved to R11 in non-reversed
NB.     order
NB.   - ircond get its default value in caller geprf

geprb=: 4 : 0
  k=. <./ 'm n'=. $ y
  p=. i. n
  if. 0 = k do.
    ((2#m)$0) ; ((0,n)$0) ; p ; 0 ; 0 ; 4$0
    return.
  end.
  y=. y , 0
  NB. Move column with largest residual norm up front
  'y A01 A11 p v svlues mxnm'=. geprc ((0,n)$0);y;1;x;p;(i.0);(4$_.);_.
  if. c y do.
    if. 1 = k do.
      (m ungqr y) ; (tru }: y) ; p ; (%/ 1 0 { svlues) ; 1 ; svlues
      return.
    else.
      smin=. 1 { svlues
    end.
  else.
    (m ungqr A11) ; (tru }: A11) ; p ; 0 ; 0 ; 4$0
    return.
  end.
  NB. Factor remaining columns using blocked code with
  NB. restricted pivoting strategy
  nb=. QFNB <. k
  NB. The size of the pivot window is chosen to be nb+nllity
  nllity=. k <. 10 >. <. 0.5 0.05 mp nb , n
  norej=. n
  while. (c y) < k <. norej do.
    NB. Invariant: the 1st column of Ax1 is the 1st column
    NB. in currently considered block column
    kb=. <./ nb , (k , norej) - c y
    NB. The goal now is to find kb independent columns
    NB. among the remaining norej not yet rejected columns
    lwsize=. (norej - c y) <. kb + nllity
    'A00 A01'=. lwsize ({."1 ; }."1) A01
    'A10 A11'=. lwsize ({."1 ; }."1) A11
    pwi=. (i. lwsize) + c y
    'dQfR A00 A10 pw smin v'=. geprw A00;A10;kb;x;(pwi{p);mxnm;smin;v
    p=. pw pwi} p
    if. c dQfR do.
      NB. Accumulate Householder vectors in a block reflector
      if. c A11 do.
        NB. Apply block reflector to A11
        A11=. ((-c y) trl1 dQfR) larfblcfc A11
      end.
      y=. y ,. dQfR
    end.
    NB. Move rejected columns to the end if there is space
    if. kb > c dQfR do.
      kklwsz=. y +&c A10
      if. norej > kklwsz do.
        'ppfx psfx'=. (c y) ({. ; }.) p
        p=. ppfx , (c A10) |. psfx
        A01=. (A01 , (c dQfR) {. A11) ,. A00
        A11=. ((c dQfR) }. A11) ,. A10
        norej=. norej - c A10
      else.
        A01=. A00 ,. A01 , (c dQfR) {. A11
        A11=. A10 ,. (c dQfR) }. A11
        break.
      end.
    else.
      A01=. A00 ,. A01 , (c dQfR) {. A11
      A11=. A10 ,. (c dQfR) }. A11
    end.
  end.
  svlues=. (smin , mxnm * 3 %: c y) 1 0} svlues
  if. k > c y do.
    NB. Process rejected columns
    'ppfx psfx'=. (c y) ({. ; }.) p
    'eApfx A01 A11 psfx v svlues mxnm'=. geprc A01;A11;(c A11);x;psfx;v;svlues;mxnm
    p=. ppfx , psfx
    y=. y ,. eApfx
  else.
    svlues=. smin 2 3} svlues
  end.
  rank=. c y
  y=. y ,. A01 , A11
  (m ungqr y) ; (tru }: y) ; p ; (%/ 1 0 { svlues) ; rank ; svlues
)

NB. ---------------------------------------------------------
NB. trlpr
NB.
NB. Description:
NB.   Compute an estimate for the numerical rank of an
NB.   lower trapezoidal matrix
NB.
NB. Syntax:
NB.   rank=. rcond trlpr L
NB. where
NB.   L     - m×k-matrix, the lower trapezoidal, k≤m
NB.   rcond > 0, 1/rcond specifies an upper bound on the
NB.           condition number of L
NB.   rank  ≥ 0, an estimate for the numerical rank of L

trlpr=: 4 : 0
  smin=. smax=. | 0 ({,) y
  rank=. 0 < smin
  if. rank do.
    k=. c y
    x1=. x2=. 1
    while. rank < k do.
      'w gamma'=. (}: ; {:) (< rank ; i. >: rank) { y
      'smin cs1'=. laic12 smin ; gamma , x1 mp + w
      'smax cs2'=. laic11 smax ; gamma , x2 mp + w
      if. smin < smax*x do. break. end.
      x1=. x1 ((* {:) , 0 { ]) cs1
      x2=. x2 ((* {:) , 0 { ]) cs2
      rank=. >: rank
    end.
  end.
  rank
)

NB. ---------------------------------------------------------
NB. trprr
NB.
NB. Description:
NB.   Compute an estimate for the numerical rank of an
NB.   upper trapezoidal matrix
NB.
NB. Syntax:
NB.   rank=. rcond trprr R
NB. where
NB.   R     - k×n-matrix, the upper trapezoidal, k≤n
NB.   rcond > 0, 1/rcond specifies an upper bound on the
NB.           condition number of R
NB.   rank  ≥ 0, an estimate for the numerical rank of R
NB.
NB. Notes:
NB. - models RRQR's xTRRNK [1, 2] with following
NB.   difference: rectangular R is allowed

trprr=: 4 : 0
  smin=. smax=. | 0 ({,) y
  rank=. 0 < smin
  if. rank do.
    k=. # y
    x1=. x2=. 1
    while. rank < k do.
      'w gamma'=. (}: ; {:) (< (i. >: rank) ; rank) { y
      'smin cs1'=. laic12 smin ; gamma , w mp + x1
      'smax cs2'=. laic11 smax ; gamma , w mp + x2
      if. smin < smax*x do. break. end.
      x1=. x1 ((* {:) , 0 { ]) cs1
      x2=. x2 ((* {:) , 0 { ]) cs2
      rank=. >: rank
    end.
  end.
  rank
)

NB. ---------------------------------------------------------
NB. gelpg3
NB.
NB. Description:
NB.   Retriangularize a special matrix, and accumulate
NB.   rotations applied
NB.
NB. Syntax:
NB.   'dQ L'=. gelpg3 A
NB. where
NB.   A  - m×k-matrix, has zeros above 1st subdiagonal,
NB.        excepting the first row
NB.   L  - m×k-matrix, the lower trapezoidal
NB.   dQ - r×4-matrix, rotations accumulated, where each row
NB.        defines one rotation and is 4-vector of values:
NB.          c , s , iof , iog
NB.
NB. Storage layout:
NB. - example for input 6×4-matrix A:
NB.     x x x x
NB.     x 0 0 0
NB.     x x 0 0
NB.     x x x 0
NB.     x x x x
NB.     x x x x

gelpg3=: 3 : 0
  'm k'=. $ y
  dQ=. 0 4$0
  if. (1 < k) *. 0 < m do.
    NB. Compute Givens rotations needed to nullify the first
    NB. row of matrix A, apply to A, accumulate rotations
    NB. applied
    i=. <: k
    liso=. i. m
    while. i do.
      'y cs'=. rot&.|: rotga y ; (< liso ; i - 1 0) ; 0
      dQ=. dQ , (+ cs) , i - 1 0
      i=. <: i
    end.
  end.
  dQ ; y
)

NB. ---------------------------------------------------------
NB. geprg3
NB.
NB. Description:
NB.   Retriangularize a special matrix, and accumulate
NB.   rotations applied
NB.
NB. Syntax:
NB.   'dQ R'=. geprg3 A
NB. where
NB.   A  - k×n-matrix, has zeros below 1st superdiagonal,
NB.        excepting the first column
NB.   R  - k×n-matrix, the upper trapezoidal
NB.   dQ - r×4-matrix, rotations accumulated, where each row
NB.        defines one rotation and is 4-vector of values:
NB.          c , s , iof , iog
NB.
NB. Storage layout:
NB. - example for input 4×6-matrix A:
NB.     x x x x x x
NB.     x 0 x x x x
NB.     x 0 0 x x x
NB.     x 0 0 0 x x
NB.
NB. Notes:
NB. - models RRQR's xGRET(3) [1, 2] with following
NB.   differences:
NB.   - matrix C is an identity matrix
NB.   - Givens rotations:
NB.     - are applied to A in rows-oriented manner
NB.     - are accumulated in dQ instead of applying to Q
NB.       immediately
NB.   - block algorithm only

geprg3=: 3 : 0
  'k n'=. $ y
  dQ=. 0 4$0
  if. (1 < k) *. 0 < n do.
    NB. Compute Givens rotations needed to nullify the first
    NB. column of matrix A, apply to A, accumulate rotations
    NB. applied
    i=. <: k
    liso=. i. n
    while. i do.
      'y cs'=. rot rotga y ; (< (i - 1 0) ; liso) ; < < a: ; 0
      dQ=. dQ , (+ cs) , i - 1 0
      i=. <: i
    end.
  end.
  dQ ; y
)

NB. ---------------------------------------------------------
NB. hslph3
NB.
NB. Description:
NB.   Reduce the lower Hessenberg matrix A to lower
NB.   triangular form, and accumulate rotations applied
NB.
NB. Syntax:
NB.   'dQ L'=. hslph3 H
NB. where
NB.   H  - m×k-matrix, the lower Hessenberg, k≤m
NB.   L  - m×k-matrix, the lower trapezoidal
NB.   dQ - r×4-matrix, rotations accumulated, where each row
NB.        defines one rotation and is 4-vector of values:
NB.          c , s , iof , iog

hslph3=: 3 : 0
  'm k'=. $ y
  dQ=. 0 4$0
  if. (1 < k) *. 0 < m do.
    NB. Compute Givens rotations needed to reduce lower
    NB. Hessenberg matrix H to lower triangular form L,
    NB. apply to H, accumulate rotations applied
    i=. 1
    liso=. i. m
    while. i < k do.
      'y cs'=. rot&.|: rotga y ; (< liso ; i - 1 0) ; 0
      liso=. }. liso
      dQ=. dQ , (+ cs) , i - 1 0
      i=. >: i
    end.
  end.
  dQ ; y
)

NB. ---------------------------------------------------------
NB. hsprh3
NB.
NB. Description:
NB.   Reduce the upper Hessenberg matrix A to upper
NB.   triangular form, and accumulate rotations applied
NB.
NB. Syntax:
NB.   'dQ R'=. hsprh3 H
NB. where
NB.   H  - k×n-matrix, the upper Hessenberg, k≤n
NB.   R  - k×n-matrix, the upper trapezoidal
NB.   dQ - r×4-matrix, rotations accumulated, where each row
NB.        defines one rotation and is 4-vector of values:
NB.          c , s , iof , iog
NB.
NB. Notes:
NB. - models RRQR's xHESS(3) [1, 2] with following
NB.   differences:
NB.   - matrix C is an identity matrix
NB.   - Givens rotations:
NB.     - are applied to A in rows-oriented manner
NB.     - are accumulated in dQ instead of applying to Q
NB.       immediately
NB.   - block algorithm only

hsprh3=: 3 : 0
  'k n'=. $ y
  dQ=. 0 4$0
  if. (1 < k) *. 0 < n do.
    NB. Compute Givens rotations needed to reduce upper
    NB. Hessenberg matrix H to upper triangular form R,
    NB. apply to H, accumulate rotations applied
    i=. 1
    liso=. i. n
    while. i < k do.
      'y cs'=. rot rotga y ; (< (i - 1 0) ; liso) ; < < a: ; 0
      liso=. }. liso
      dQ=. dQ , (+ cs) , i - 1 0
      i=. >: i
    end.
  end.
  dQ ; y
)

NB. ---------------------------------------------------------
NB. trlpc
NB.
NB. Description:
NB.   Carry out Pan-Tang Algorithm 3 for the stage rank. This
NB.   is a mofified version of the original algorithm. The
NB.   improved features are the following:
NB.   - use of Bischof's ICE to reduce the computational cost
NB.   - reorganization of the main loop to save computations
NB.   - no permutation is carried out if not strictly needed
NB.
NB. Syntax:
NB.   'op oL dQ isovrn rcnr rcnrp1 svlues'=. trprc ip;iL;rank
NB. where
NB.   ip     - m-vector, rows inversed permutation of iL,
NB.            represents permutation m×m-matrix
NB.   iL     - m×k-matrix, input lower trapezoidal
NB.   rank   - integer in range [0,k], the estimate for the
NB.            rank of iL
NB.   op     - m-vector, rows inversed permutation of oL,
NB.            represents permutation m×m-matrix
NB.   oL     - m×k-matrix, output lower trapezoidal
NB.   dQ     - r×4-matrix, rotations accumulated, where each
NB.            row defines one rotation and is 4-vector of
NB.            values:
NB.              c , s , iof , iog
NB.   isovrn - boolean, signals about allowed maximum number
NB.            of steps has beed exceeded, that is, the
NB.            matrix presents a slow convergence
NB.   rcnr   - > 0, 1/rcnr specifies an estimate of the
NB.            condition number of L00
NB.   rcnrp1 - > 0, 1/rcnrp1 specifies an estimate of the
NB.            condition number of L00p1
NB.   svlues - 4-vector, estimates of the singular values,
NB.            see gelpf
NB.   A      - m×n-matrix, the input to factorize
NB.   k      = min(m,n)
NB.
NB. Storage layout:
NB.   L  =  (  L00   0    ) rank    =  (  L00p1   0  ) rank+1
NB.         (  L10   R11  ) m-rank     (  *       *  ) m-rank-1
NB.            rank  n-rank               rank+1  n-rank-1
NB.
NB. Notes:
NB. - if the leading block of L is singular or near singular,
NB.   there will be no permutation because in that case the
NB.   right (and left) singular vectors are the canonical
NB.   ones ((0,0,...,0,1)^T). That is, there will not be
NB.   permutation if rcond <: PFSF * FP_SFMIN

trlpc=: 3 : 0
  'p L rank'=. y
  'm k'=. $ L
  k1=. <: k
  dQ=. 0 4 $ 0
  if. 0 = k do.
    svlues=. 4 $ isovrn=. rcnr=. rcnrp1=. rank=. 0
  elseif. rank = k do.
    NB. Apply Chan algorithm
    f=. PFFP
    NB. Move the best row of A(k1:m-1,0:n-1) to (k-1)-th
    NB. position
    io=. liofmax (< m (th2liso ; ]) k1) { L
    if. io do.
      dip=. < k1 ([ , +) io
      L=. dip C. L
      p=. dip C. p
    end.
    NB. Estimate the largest singular value, the smallest
    NB. singular value, and its corresponding right singular
    NB. vector
    smin=. smax=. | 0 ({,) L
    x1=. x2=. j=. 1
    while. j < rank do.
      'w gamma'=. (}: ; {:) (< j ; i. >: j) { L
      'smax cs'=. laic11 smax ; gamma , w mp + x1
      x1=. x1 ((* {:) , 0 { ]) cs
      'smin cs'=. laic12 smin ; gamma , w mp + x2
      x2=. x2 ((* {:) , 0 { ]) cs
      j=. >: j
    end.
    NB. Determine if matrix A is singular or nearly singular
    if. smin > smax * FP_SFMIN * PFSF do.
      NB. Matrix is not singular or not nearly singular.
      NB. Follow usual method: Estimate the left singular
      NB. vector corresponding to the smallest singular value
      NB. of lower triangular block A(0:rank-1,0:rank-1)
      x2=. ((2 # rank) {. L) trsmrlnn x2
      NB. Find the index with largest absolute value in
      NB. vector x2
      io=. liofmax x2
      NB. Permut if necessary
      if. ((f * | io { x2) > | {: x2) *. (io < <: rank) do.
        NB. Exchange cyclically to the left the
        NB. rows/elements between io and rank-1, that is:
        NB. io->rank-1, io+1->io, io+2->io+1,...,
        NB. rank-1->rank-2
        dip=. < 1 |. rank th2liso io
        L=. dip C. L
        p=. dip C. p
        NB. Retriangularize matrix A after the permutation
        'dQi L'=. hslph3 L
        dQ=. dQ , dQi
      end.
    end.
    isovrn=. 0
    svlues=. 1 3#smax,smin
    rcnr=. rcnrp1=. smin%smax
  else.
    NB. Apply modified Pan&Tang algorithm 3
    ns=. 0
    mxstps=. 25 + m
    f=. PFFP % %: >: rank
    NB. Compute the norms of rows of matrix A
    rnorms=. normsr L
    NB. Estimate the smallest singular value of
    NB. A(0:rank-1,0:rank-1) and its corresponding right
    NB. singular vector. smin will contain the smallest
    NB. singular value and x1 will contain the right singular
    NB. vector
    smin=. | 0 ({,) L
    x1=. j=. 1
    while. j < rank do.
      'w gamma'=. (}: ; {:) (< j ; i. >: j) { L
      'smin cs'=. laic12 smin ; gamma , w mp + x1
      x1=. x1 ((* {:) , 0 { ]) cs
      j=. >: j
    end.
    NB. Initialize loop variables
    nca=. 0
    nctba=. m-rank
    ii=. rank
    while. (ns < mxstps) *. (nca < nctba) do.
      NB. Estimate the smallest singular value of
      NB. A(0:rank,0:rank) and its corresponding right
      NB. singular vector as if row ii of matrix A were on
      NB. row of number rank
      ij=. ii <. <: k
      diag=. (< ii , ij) { L
      i=. <: ij
      while. i >: rank do.
        diag=. (mp~ lartg) ((< ii , i ) { L) , diag
        i=. <: i
      end.
      'mnrp1 cs'=. laic12 smin ; diag , ((< ii ; i. rank) { L) mp + x1
      if. mnrp1 >: f * | diag do.
        NB. Row ii accepted on the bottom part of matrix A
        nca=. >: nca
        if. ii = <: m do.
          ii=. >: rank
        else.
          ii=. >: ii
        end.
      else.
        NB. Row ii not accepted on the left part of matrix A
        NB.
        NB. Permut row ii to position number rank
        NB.
        NB. Exchange cyclically to the bottom the
        NB. rows/elements between rank and ii, that is,
        NB. rank->rank+1, rank+1->rank+2,...,ii-1->ii,
        NB. ii->rank
        dip=. < |. (>: ii) th2liso rank
        L=. dip C. L
        p=. dip C. p
        rnorms=. dip C. rnorms
        NB. Retriangularize matrix A after the permutation,
        NB. adjust Q accordingly
        'dQi L'=. gelpg3 L
        dQ=. dQ , dQi
        NB. Estimate the largest singular value
        mxrp1=. (3 %: >: rank) * normitr (>: rank) {. rnorms
        NB. Estimate the left singular vector
        if. mnrp1 > mxrp1 * PFSF * FP_SFMIN do.
          NB. Matrix is not singular or not nearly singular
          NB.
          NB. First, end the estimation of the right singular
          NB. vector
          x2=. x1 ((* {:) , 0 { ]) cs
          NB. Obtain the left singular vector from the right
          NB. one
          x2=. ((2 # >: rank) {. L) trsmrlnn x2
          NB. Find the index with largest absolute value
          io=. liofmax x2
          NB. Permut row io to position rank
          if. io < rank do.
            NB. Exchange cyclically to the top the
            NB. rows/elements between io and rank, that
            NB. is, io->rank,io+1->io,io+2->io+1,...,
            NB. rank->rank-1
            dip=. < 1 |. (>: rank) th2liso io
            L=. dip C. L
            p=. dip C. p
            rnorms=. dip C. rnorms
            NB. Retriangularize matrix A after the
            NB. permutation
            'dQi L'=. hslph3 L
            dQ=. dQ , dQi
            NB. Estimate the smallest singular value of
            NB. A(0:rank-1,0:rank-1) and its corresponding
            NB. right singular vector. smin will contain the
            NB. smallest singular value and x1 will contain
            NB. the right singular vector
            smin=. | 0 ({,) L
            x1=. j=. 1
            while. j < rank do.
              'w gamma'=. (}: ; {:) (< j ; i. >: j) { L
              'smin cs'=. laic12 smin ; gamma , x1 mp + w
              x1=. x1 ((* {:) , 0 { ]) cs
              j=. >: j
            end.
          end.
        end.
        NB. Update loop variables
        nca=. 0
        ns=. >: ns
        if. ii = <: m do.
          ii=. >: rank
        else.
          ii=. >: ii
        end.
      end.
    end.
    NB. Final pivoting
    io=. rank + liofmax rank }. rnorms
    if. ((f * | io { rnorms) > | rank { rnorms) *. (io > rank) do.
      NB. Exchange cyclically to the bottom the
      NB. rows/elements between rank and ii, that is,
      NB. rank->rank+1, rank+1->rank+2,...,io-1->io,
      NB. io->rank
      dip=. < |. (>: io) th2liso rank
      L=. dip C. L
      p=. dip C. p
      NB. Retriangularize matrix A after permutation
      'dQi L'=. gelpg3 L
      dQ=. dQ , dQi
    end.
    NB. Compute vector svlues and variables rcnr and rcnrp1
    NB.
    NB. Compute the largest singular value of
    NB. A(0:rank-1,0:rank-1)
    smax=. | 0 ({,) L
    x2=. j=. 1
    while. j < rank do.
      'w gamma'=. (}: ; {:) (< j ; i. >: j) { L
      'smax cs'=. laic11 smax ; gamma , x2 mp + w
      x2=. x2 ((* {:) , 0 { ]) cs
      j=. >: j
    end.
    svlues=. smax , smin
    NB. Compute the largest singular value and the smallest
    NB. singular value of A(0:rank,0:rank)
    'w gamma'=. (}: ; {:) (< rank ; i. >: rank) { L
    'mxrp1 cs'=. laic11 smax ; gamma , x2 mp + w
    'smin cs'=. laic12 smin ; gamma , x1 mp + w
    x1=. x1 ((* {:) , 0 { ]) cs
    svlues=. svlues , smin
    NB. Compute the smallest singular value of
    NB. A(0:k-1,0:k-1)
    j=. >: rank
    while. j < k do.
      'w gamma'=. (}: ; {:) (< j ; i. >: j) { L
      'smin cs'=. laic12 smin ; gamma , x1 mp + w
      x1=. x1 ((* {:) , 0 { ]) cs
      j=. >: j
    end.
    isovrn=. ns >: mxstps
    svlues=. svlues , smin
    NB. Compute rcnr and rcnrp1
    rcnr=. %/ 1 0 { svlues
    rcnrp1=. (2 { svlues) % mxrp1
  end.
  p ; L ; dQ ; isovrn ; rcnr ; rcnrp1 ; svlues
)

NB. ---------------------------------------------------------
NB. trprc
NB.
NB. Description:
NB.   Carry out Pan-Tang Algorithm 3 for the stage rank. This
NB.   is a mofified version of the original algorithm. The
NB.   improved features are the following:
NB.   - use of Bischof's ICE to reduce the computational cost
NB.   - reorganization of the main loop to save computations
NB.   - no permutation is carried out if not strictly needed
NB.
NB. Syntax:
NB.   'dQ oR op isovrn rcnr rcnrp1 svlues'=. trprc iR;ip;rank
NB. where
NB.   iR     - k×n-matrix, input upper trapezoidal
NB.   ip     - n-vector, columns inversed permutation of iR,
NB.            represents permutation n×n-matrix
NB.   rank   - integer in range [0,k], the estimate for the
NB.            rank of iR
NB.   dQ     - r×4-matrix, rotations accumulated, where each
NB.            row defines one rotation and is 4-vector of
NB.            values:
NB.              c , s , iof , iog
NB.   oR     - k×n-matrix, output upper trapezoidal
NB.   op     - n-vector, columns inversed permutation of oR,
NB.            represents permutation n×n-matrix
NB.   isovrn - boolean, signals about allowed maximum number
NB.            of steps has beed exceeded, that is, the
NB.            matrix presents a slow convergence
NB.   rcnr   - > 0, 1/rcnr specifies an estimate of the
NB.            condition number of R00
NB.   rcnrp1 - > 0, 1/rcnrp1 specifies an estimate of the
NB.            condition number of R00p1
NB.   svlues - 4-vector, estimates of the singular values,
NB.            see geprf
NB.   A      - m×n-matrix, the input to factorize
NB.   k      = min(m,n)
NB.
NB. Storage layout:
NB.   R  =  (  R00   R01  ) rank    =  (  R00p1   *  ) rank+1
NB.         (  0     R11  ) m-rank     (  0       *  ) m-rank-1
NB.            rank  n-rank               rank+1  n-rank-1
NB.
NB. Notes:
NB. - models RRQR's xTRQYC(3) [1, 2] with following
NB.   difference: matrix C is an identity matrix
NB. - if the leading block of R is singular or near singular,
NB.   there will be no permutation because in that case the
NB.   right (and left) singular vectors are the canonical
NB.   ones ((0,0,...,0,1)^T). That is, there will not be
NB.   permutation if rcond <: PFSF * FP_SFMIN

trprc=: 3 : 0
  'R p rank'=. y
  'k n'=. $ R
  k1=. <: k
  dQ=. 0 4 $ 0
  if. 0 = k do.
    svlues=. 4 $ isovrn=. rcnr=. rcnrp1=. rank=. 0
  elseif. rank = k do.
    NB. Apply Chan algorithm
    f=. PFFP
    NB. Move the best column of A(0:m-1,k1:n-1) to (k-1)-th
    NB. position
    io=. liofmax (< n (] ; th2liso) k1) { R
    if. io do.
      dip=. < k1 ([ , +) io
      R=. dip C."1 R
      p=. dip C. p
    end.
    NB. Estimate the largest singular value, the smallest
    NB. singular value, and its corresponding left singular
    NB. vector
    smin=. smax=. | 0 ({,) R
    x1=. x2=. j=. 1
    while. j < rank do.
      'w gamma'=. (}: ; {:) (< (i. >: j) ; j) { R
      'smax cs'=. laic11 smax ; gamma , w mp + x1
      x1=. x1 ((* {:) , 0 { ]) cs
      'smin cs'=. laic12 smin ; gamma , w mp + x2
      x2=. x2 ((* {:) , 0 { ]) cs
      j=. >: j
    end.
    NB. Determine if matrix A is singular or nearly singular
    if. smin > smax * FP_SFMIN * PFSF do.
      NB. Matrix is not singular or not nearly singular.
      NB. Follow usual method: Estimate the right singular
      NB. vector corresponding to the smallest singular value
      NB. of upper triangular block A(0:rank-1,0:rank-1)
      x2=. ((2 # rank) {. R) trsmlunn x2
      NB. Find the index with largest absolute value in
      NB. vector x2
      io=. liofmax x2
      NB. Permut if necessary
      if. ((f * | io { x2) > | {: x2) *. (io < <: rank) do.
        NB. Exchange cyclically to the left the
        NB. columns/elements between io and rank-1, that is:
        NB. io->rank-1, io+1->io, io+2->io+1,...,
        NB. rank-1->rank-2
        dip=. < 1 |. rank th2liso io
        R=. dip C."1 R
        p=. dip C. p
        NB. Retriangularize matrix A after the permutation
        'dQi R'=. hsprh3 R
        dQ=. dQ , dQi
      end.
    end.
    isovrn=. 0
    svlues=. 1 3#smax,smin
    rcnr=. rcnrp1=. smin%smax
  else.
    NB. Apply modified Pan&Tang algorithm 3
    ns=. 0
    mxstps=. 25 + n
    f=. PFFP % %: >: rank
    NB. Compute the norms of columns of matrix A
    cnorms=. normsc R
    NB. Estimate the smallest singular value of
    NB. A(0:rank-1,0:rank-1) and its corresponding left
    NB. singular vector. smin will contain the smallest
    NB. singular value and x1 will contain the left singular
    NB. vector
    smin=. | 0 ({,) R
    x1=. j=. 1
    while. j < rank do.
      'w gamma'=. (}: ; {:) (< (i. >: j) ; j) { R
      'smin cs'=. laic12 smin ; gamma , w mp + x1
      x1=. x1 ((* {:) , 0 { ]) cs
      j=. >: j
    end.
    NB. Initialize loop variables
    nca=. 0
    nctba=. n-rank
    ii=. rank
    while. (ns < mxstps) *. (nca < nctba) do.
      NB. Estimate the smallest singular value of
      NB. A(0:rank,0:rank) and its corresponding left
      NB. singular vector as if column ii of matrix A were on
      NB. column of number rank
      ij=. ii <. <: k
      diag=. (< ij , ii) { R
      i=. <: ij
      while. i >: rank do.
        diag=. (mp~ lartg) ((< i , ii ) { R) , diag
        i=. <: i
      end.
      'mnrp1 cs'=. laic12 smin ; diag , ((< (i. rank) ; ii) { R) mp + x1
      if. mnrp1 >: f * | diag do.
        NB. Column ii accepted on the right part of matrix A
        nca=. >: nca
        if. ii = <: n do.
          ii=. >: rank
        else.
          ii=. >: ii
        end.
      else.
        NB. Column ii not accepted on the right part of
        NB. matrix A
        NB.
        NB. Permut column ii to position number rank
        NB.
        NB. Exchange cyclically to the right the
        NB. columns/elements between rank and ii, that is,
        NB. rank->rank+1, rank+1->rank+2,...,ii-1->ii,
        NB. ii->rank
        dip=. < |. (>: ii) th2liso rank
        R=. dip C."1 R
        p=. dip C. p
        cnorms=. dip C. cnorms
        NB. Retriangularize matrix A after the permutation,
        NB. adjust Q accordingly
        'dQi R'=. geprg3 R
        dQ=. dQ , dQi
        NB. Estimate the largest singular value
        mxrp1=. (3 %: >: rank) * normitr (>: rank) {. cnorms
        NB. Estimate the right singular vector
        if. mnrp1 > mxrp1 * PFSF * FP_SFMIN do.
          NB. Matrix is not singular or not nearly singular
          NB.
          NB. First, end the estimation of the left singular
          NB. vector
          x2=. x1 ((* {:) , 0 { ]) cs
          NB. Obtain the right singular vector from the left
          NB. one
          x2=. ((2 # >: rank) {. R) trsmlunn x2
          NB. Find the index with largest absolute value
          io=. liofmax x2
          NB. Permut column io to position rank
          if. io < rank do.
            NB. Exchange cyclically to the left the
            NB. columns/elements between io and rank, that
            NB. is, io->rank,io+1->io,io+2->io+1,...,
            NB. rank->rank-1
            dip=. < 1 |. (>: rank) th2liso io
            R=. dip C."1 R
            p=. dip C. p
            cnorms=. dip C. cnorms
            NB. Retriangularize matrix A after the
            NB. permutation
            'dQi R'=. hsprh3 R
            dQ=. dQ , dQi
            NB. Estimate the smallest singular value of
            NB. A(0:rank-1,0:rank-1) and its corresponding
            NB. left singular vector. smin will contain the
            NB. smallest singular value and x1 will contain
            NB. the left singular vector
            smin=. | 0 ({,) R
            x1=. j=. 1
            while. j < rank do.
              'w gamma'=. (}: ; {:) (< (i. >: j) ; j) { R
              'smin cs'=. laic12 smin ; gamma , w mp + x1
              x1=. x1 ((* {:) , 0 { ]) cs
              j=. >: j
            end.
          end.
        end.
        NB. Update loop variables
        nca=. 0
        ns=. >: ns
        if. ii = <: n do.
          ii=. >: rank
        else.
          ii=. >: ii
        end.
      end.
    end.
    NB. Final pivoting
    io=. rank + liofmax rank }. cnorms
    if. ((f * | io { cnorms) > | rank { cnorms) *. (io > rank) do.
      NB. Exchange cyclically to the right the
      NB. columns/elements between rank and ii, that is,
      NB. rank->rank+1, rank+1->rank+2,...,io-1->io,
      NB. io->rank
      dip=. < |. (>: io) th2liso rank
      R=. dip C."1 R
      p=. dip C. p
      NB. Retriangularize matrix A after permutation
      'dQi R'=. geprg3 R
      dQ=. dQ , dQi
    end.
    NB. Compute vector svlues and variables rcnr and rcnrp1
    NB.
    NB. Compute the largest singular value of
    NB. A(0:rank-1,0:rank-1)
    smax=. | 0 ({,) R
    x2=. j=. 1
    while. j < rank do.
      'w gamma'=. (}: ; {:) (< (i. >: j) ; j) { R
      'smax cs'=. laic11 smax ; gamma , w mp + x2
      x2=. x2 ((* {:) , 0 { ]) cs
      j=. >: j
    end.
    svlues=. smax , smin
    NB. Compute the largest singular value and the smallest
    NB. singular value of A(0:rank,0:rank)
    'w gamma'=. (}: ; {:) (< (i. >: rank) ; rank) { R
    'mxrp1 cs'=. laic11 smax ; gamma , w mp + x2
    'smin cs'=. laic12 smin ; gamma , w mp + x1
    x1=. x1 ((* {:) , 0 { ]) cs
    svlues=. svlues , smin
    NB. Compute the smallest singular value of
    NB. A(0:k-1,0:k-1)
    j=. >: rank
    while. j < k do.
      'w gamma'=. (}: ; {:) (< (i. >: j) ; j) { R
      'smin cs'=. laic12 smin ; gamma , w mp + x1
      x1=. x1 ((* {:) , 0 { ]) cs
      j=. >: j
    end.
    isovrn=. ns >: mxstps
    svlues=. svlues , smin
    NB. Compute rcnr and rcnrp1
    rcnr=. %/ 1 0 { svlues
    rcnrp1=. (2 { svlues) % mxrp1
  end.
  dQ ; R ; p ; isovrn ; rcnr ; rcnrp1 ; svlues
)

NB. ---------------------------------------------------------
NB. trlpy
NB.
NB. Description:
NB.   Detect the right rank for lower trapezoidal matrix L.
NB.   The algorithm used here is an version of Pan and Tang's
NB.   RRQR algorithm number 3. This algorithm is applied to
NB.   matrix L until the right rank is obtained. If the input
NB.   ordering of matrix L is not accepted, the matrix will
NB.   be permuted and retriangularized until the rank is
NB.   revealed.
NB.
NB. Syntax:
NB.   'op oL oQn orcond rank svlues'=. ircond trlpy ip;iL;iQn;trash;trash;trash
NB. where
NB.   ircond > 0, 1/ircond is an upper bound on the condition
NB.            number of iL00
NB.   ip     - m-vector, rows inversed permutation of iL,
NB.            represents permutation m×m-matrix
NB.   iL     - m×k-matrix, the lower trapezoidal
NB.            before rank detecting
NB.   iQn    - n×n-matrix, the unitary (orthogonal), embeds
NB.            iQ:
NB.              iQ -: k {. iQn
NB.   op     - m-vector, rows inversed permutation of oL,
NB.            represents permutation m×m-matrix
NB.   oL     - m×k-matrix, the lower trapezoidal
NB.            after rank detecting
NB.   oQn    - n×n-matrix, the unitary (orthogonal), embeds
NB.            oQ:
NB.              oQ -: k {. oQn
NB.   orcond > 0, 1/orcond is an estimate for the condition
NB.            number of oL00
NB.   rank   - integer in range [0,k], an estimate of the
NB.            rank offered by this algorithm
NB.   svlues - 4-vector, estimates of the singular values,
NB.            see gelpf
NB.   iQ     - k×n-matrix with orthonormal rows, which is
NB.            defined as the first k rows of the product of
NB.            k elementary reflectors, before rank detecting
NB.   oQ     - k×n-matrix with orthonormal rows, which is
NB.            defined as the first k rows of the product of
NB.            k elementary reflectors, after rank detecting
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   L  =  (  L00   0    ) rank    =  (  L00p1   0  ) rank+1
NB.         (  L10   L11  ) m-rank     (  *       *  ) m-rank-1
NB.            rank  n-rank               rank+1  n-rank-1

trlpy=: 4 : 0
  'p L Q'=. 3 {. y
  k=. c L
  dQ=. 0 4 $ 0
  if. k do.
    NB. Compute the initial estimate for the rank
    rank=. x trlpr L
    NB. Loop for the detection of the actual rank. The
    NB. variable rank is updated until the rank is found. To
    NB. avoid infinite loops, the variable rank either
    NB. increases or decreases.
    nornkdtd=. is1st=. godown=. 1
    whilst. nornkdtd do.
      NB. Get tighter bounds for the value rank
      'p L dQi isovrn rcnr rcnrp1 svlues'=. trlpc p;L;rank
      dQ=. dQ , dQi
      NB. Check if the numerical rank is larger, equal or
      NB. smaller than the contents of rank
      if. ((rcnrp1 < x) +. rank = k) *. rcnr >: x do.
        nornkdtd=. 0
      elseif. x <: rcnr <. rcnrp1 do.
        if. is1st +. godown do.
          rank=. >: rank
        else.
          nornkdtd=. 0
          rank=. _.
        end.
      elseif. x > rcnr >. rcnrp1 do.
        if. 1 = rank do.
          nornkdtd=. 0
          if. 0 = | 0 ({,) L do.
            svlues=. 4 $ rcnr=. rank=. 0
          end.
        else.
          godown=. 0
          rank=. <: rank
        end.
      else.
        nornkdtd=. 0
        rank=. _.
      end.
      is1st=. 0
    end.
    if. isovrn do. rank=. _. end.
  else.
    svlues=. 4 $ rcnr=. rank=. 0
  end.
  p ; L ; (Q rotscll dQ) ; rcnr ; rank ; svlues
)

NB. ---------------------------------------------------------
NB. trpry
NB.
NB. Description:
NB.   Detect the right rank for upper trapezoidal matrix R.
NB.   The algorithm used here is an version of Pan and Tang's
NB.   RRQR algorithm number 3. This algorithm is applied to
NB.   matrix R until the right rank is obtained. If the input
NB.   ordering of matrix R is not accepted, the matrix will
NB.   be permuted and retriangularized until the rank is
NB.   revealed.
NB.
NB. Syntax:
NB.   'oQm oR op orcond rank svlues'=. ircond trpry iQm;iR;ip;trash;trash;trash
NB. where
NB.   ircond > 0, 1/ircond is an upper bound on the condition
NB.            number of iR00
NB.   iQm    - m×m-matrix, the unitary (orthogonal), embeds
NB.            iQ:
NB.              iQ -: k {."1 iQm
NB.   iR     - k×n-matrix, the upper trapezoidal
NB.            before rank detecting
NB.   ip     - n-vector, columns inversed permutation of iR,
NB.            represents permutation n×n-matrix
NB.   oQm    - m×m-matrix, the unitary (orthogonal), embeds
NB.            oQ:
NB.              oQ -: k {."1 oQm
NB.   oR     - k×n-matrix, the upper trapezoidal
NB.            after rank detecting
NB.   op     - n-vector, columns inversed permutation of oR,
NB.            represents permutation n×n-matrix
NB.   orcond > 0, 1/orcond is an estimate for the condition
NB.            number of oR00
NB.   rank   - integer in range [0,k], an estimate of the
NB.            rank offered by this algorithm
NB.   svlues - 4-vector, estimates of the singular values,
NB.            see geprf
NB.   iQ     - m×k-matrix with orthonormal columns, which is
NB.            defined as the first k columns of the product
NB.            of k elementary reflectors, before rank
NB.            detecting
NB.   oQ     - m×k-matrix with orthonormal columns, which is
NB.            defined as the first k columns of the product
NB.            of k elementary reflectors, after rank
NB.            detecting
NB.   k   = min(m,n)
NB.
NB. Storage layout:
NB.   R  =  (  R00   R01  ) rank    =  (  R00p1   *  ) rank+1
NB.         (  0     R11  ) m-rank     (  0       *  ) m-rank-1
NB.            rank  n-rank               rank+1  n-rank-1
NB.
NB. Notes:
NB. - models RRQR's xTRQPY [1, 2] with following differences:
NB.   - rcond get its default value in caller geprf
NB.   - in case of problems in the rank computation the rank
NB.     variable gets value NaN (_.)

trpry=: 4 : 0
  'Q R p'=. 3 {. y
  k=. # R
  dQ=. 0 4 $ 0
  if. k do.
    NB. Compute the initial estimate for the rank
    rank=. x trprr R
    NB. Loop for the detection of the actual rank. The
    NB. variable rank is updated until the rank is found. To
    NB. avoid infinite loops, the variable rank either
    NB. increases or decreases.
    nornkdtd=. is1st=. goright=. 1
    whilst. nornkdtd do.
      NB. Get tighter bounds for the value rank
      'dQi R p isovrn rcnr rcnrp1 svlues'=. trprc R;p;rank
      dQ=. dQ , dQi
      NB. Check if the numerical rank is larger, equal or
      NB. smaller than the contents of rank
      if. ((rcnrp1 < x) +. rank = k) *. rcnr >: x do.
        nornkdtd=. 0
      elseif. x <: rcnr <. rcnrp1 do.
        if. is1st +. goright do.
          rank=. >: rank
        else.
          nornkdtd=. 0
          rank=. _.
        end.
      elseif. x > rcnr >. rcnrp1 do.
        if. 1 = rank do.
          nornkdtd=. 0
          if. 0 = | 0 ({,) R do.
            svlues=. 4 $ rcnr=. rank=. 0
          end.
        else.
          goright=. 0
          rank=. <: rank
        end.
      else.
        nornkdtd=. 0
        rank=. _.
      end.
      is1st=. 0
    end.
    if. isovrn do. rank=. _. end.
  else.
    svlues=. 4 $ rcnr=. rank=. 0
  end.
  (Q rotsclu dQ) ; R ; p ; rcnr ; rank ; svlues
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gelpf
NB.
NB. Description:
NB.   LQ factorization with row pivoting of a general
NB.   matrix A:
NB.     P * (  L00   0    ) * Q = A
NB.         (  L10   L11  )
NB.   The permutation P is chosen with the goal to reveal the
NB.   rank of A by a suitably dimensioned trailing submatrix
NB.   L11 with norm(L11) being small.
NB.
NB. Syntax:
NB.   'ip L Qn orcond rank svlues'=. [ircond] gelpf A
NB. where
NB.   A      - m×n-matrix, the input to factorize
NB.   ircond > 0, optional, default is FP_EPS. 1/ircond
NB.            specifies an upper bound on the condition
NB.            number of L00
NB.   ip     - m-vector, rows inversed permutation of L,
NB.            represents permutation m×m-matrix
NB.   L      - m×k-matrix, the lower trapezoidal
NB.   Qn     - n×n-matrix, the unitary (orthogonal), embeds
NB.            Q:
NB.              Q -: k {. Qn
NB.   orcond > 0, 1/orcond is an estimate for the condition
NB.            number of L00
NB.   rank   - integer in range [0,k], an estimate of the
NB.            rank offered by this algorithm, or NaN if
NB.            there was a problem to compute a rank
NB.   svlues - 4-vector, estimates of the singular values:
NB.              sigma_max(L00) , sigma_min(L00) , sigma_min(L00p1) , sigma_min(L)
NB.            If the triangular factorization is a
NB.            rank-revealing one (which will be the case if
NB.            the leading columns were well-conditioned),
NB.            then svlues indeed contains estimates of
NB.            singular values:
NB.              sigma_max(A) , sigma_r(A) , sigma_(r+1)(A) , sigma_min(A)
NB.            By examining these values, one can confirm
NB.            that the rank is well defined with respect to
NB.            the threshold chosen
NB.   Q      - k×n-matrix with orthonormal rows, which is
NB.            defined as the first k rows of the product of
NB.            k elementary reflectors
NB.   k      = min(m,n)
NB.
NB. Storage layout:
NB.   L  =  (  L00   0    ) rank    =  (  L00p1   0  ) rank+1
NB.         (  L10   L11  ) m-rank     (  *       *  ) m-rank-1
NB.            rank  n-rank               rank+1  n-rank-1
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   _. ~: rank
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A -: p { L mp Q
NB.   A -: p C. L mp Q
NB.   A -: ip C.^:_1 L mp Q
NB.   A -: P mp L mp Q
NB.   (idmat c L) -: clean (mp ct) Q
NB. where
NB.   'ip L Qn rcond rank svlues'=. gelpf A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.   Q=. (c L) {. Qn

gelpf=: FP_EPS&$: :([ ]`trlpy@.(0 < 3 {:: ]) gelpb)

NB. NB. 'Qm L ip orcond rank svlues'=. [ircond] geplf A
NB. geplf=: FP_EPS&$: :([ ]`trply@.(0 < 3 {:: ]) geplb)

NB. ---------------------------------------------------------
NB. geprf
NB.
NB. Description:
NB.   QR factorization with column pivoting of a general
NB.   matrix A:
NB.     Q * (  R00  R01  ) * P = A
NB.         (   0   R11  )
NB.   The permutation P is chosen with the goal to reveal the
NB.   rank of A by a suitably dimensioned trailing submatrix
NB.   R11 with norm(R11) being small.
NB.
NB. Syntax:
NB.   'Qm R ip orcond rank svlues'=. [ircond] geprf A
NB. where
NB.   A      - m×n-matrix, the input to factorize
NB.   ircond > 0, optional, default is FP_EPS. 1/ircond
NB.            specifies an upper bound on the condition
NB.            number of R00
NB.   Qm     - m×m-matrix, the unitary (orthogonal), embeds
NB.            Q:
NB.              Q -: k {."1 Qm
NB.   R      - k×n-matrix, the upper trapezoidal
NB.   ip     - n-vector, columns inversed permutation of R,
NB.            represents permutation n×n-matrix
NB.   orcond > 0, 1/orcond is an estimate for the condition
NB.            number of R00
NB.   rank   - integer in range [0,k], an estimate of the
NB.            rank offered by this algorithm, or NaN if
NB.            there was a problem to compute a rank
NB.   svlues - 4-vector, estimates of the singular values:
NB.              sigma_max(R00) , sigma_min(R00) , sigma_min(R00p1) , sigma_min(R)
NB.            If the triangular factorization is a
NB.            rank-revealing one (which will be the case if
NB.            the leading columns were well-conditioned),
NB.            then svlues indeed contains estimates of
NB.            singular values:
NB.              sigma_max(A) , sigma_r(A) , sigma_(r+1)(A) , sigma_min(A)
NB.            By examining these values, one can confirm
NB.            that the rank is well defined with respect to
NB.            the threshold chosen
NB.   Q      - m×k-matrix with orthonormal columns, which is
NB.            defined as the first k columns of the product
NB.            of k elementary reflectors
NB.   k      = min(m,n)
NB.
NB. Storage layout:
NB.   R  =  (  R00   R01  ) rank    =  (  R00p1   *  ) rank+1
NB.         (  0     R11  ) m-rank     (  0       *  ) m-rank-1
NB.            rank  n-rank               rank+1  n-rank-1
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   _. ~: rank
NB.   P -: %. iP
NB.   P -: |: iP
NB.   P -: ip2P ip
NB.   A -: p {"1 Q mp R
NB.   A -: p C."1 Q mp R
NB.   A -: ip C.^:_1"1 Q mp R
NB.   A -: Q mp R mp iP         NB. apply ip to columns
NB.   (idmat # R) -: clean (mp~ ct) Q
NB. where
NB.   'Qm R ip rcond rank svlues'=. geprf A
NB.   p=. /: ip
NB.   iP=. p2P ip
NB.   P=. p2P p
NB.   Q=. (# R) {."1 Qm
NB.
NB. Notes:
NB. - simulates LAPACK's xGEQP3
NB. - implements RRQR's xGEQPY(3) [1, 2] with following
NB.   differences:
NB.   - matrix C is an identity matrix
NB.   - if there was a problem to compute a rank, then rank
NB.     gets value NaN

geprf=: FP_EPS&$: :([ ]`trpry@.(0 < 3 {:: ]) geprb)

NB. NB. 'ip R Qn orcond rank svlues'=. [ircond] gerpf A
NB. gerpf=: FP_EPS&$: :([ ]`trrpy@.(0 < 3 {:: ]) gerpb)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testgepf
NB.
NB. Description:
NB.   Test gexxf by general matrix
NB.
NB. Syntax:
NB.   testgepf A
NB. where
NB.   A - m×n-matrix
NB.
NB. TODO:
NB. - add xQRT12 test

testgepf=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/geqp3'

  rcond=. (_."_)`geconi@.(=/@$) y  NB. meaninigful for square matrices only

  norm=. norm1 y

  args=. y ; norm

  NB. current pf verbs return Q in non-factorized form,
  NB. this differs from xGEQP3 so specialized adapters
  NB. are needed
  NB. berrA=. (A ; normA) xpt01a (ip ; x ; Qn ; trash)
  NB. berrA=. (A ; normA) pxt01a (Qm ; x ; ip ; trash)
  lpt01a=: ((1 {:: [) %~^:(0 < [)   C.  ~&(0&{::)          (normi % FP_EPS * >./@$)@:- (1 {:: ]) ([ mp ({.~   c)~ ) 2 {:: ])`0:@.(0 e. $@(0 {:: [))
  plt01a=: ((1 {:: [) %~^:(0 < [) ((C."1  (0&{::))~ 2&{::) (norm1 % FP_EPS * >./@$)@:- (0 {:: ]) (({."1~ -@#) mp ]) 1 {:: ])`0:@.(0 e. $@(0 {:: [))
  prt01a=: ((1 {:: [) %~^:(0 < [) ((C."1  (0&{::))~ 2&{::) (norm1 % FP_EPS * >./@$)@:- (0 {:: ]) (({."1~   #) mp ]) 1 {:: ])`0:@.(0 e. $@(0 {:: [))
  rpt01a=: ((1 {:: [) %~^:(0 < [)   C.  ~&(0&{::)          (normi % FP_EPS * >./@$)@:- (1 {:: ]) ([ mp ({.~ -@c)~ ) 2 {:: ])`0:@.(0 e. $@(0 {:: [))
  NB. berrQn=. trash0 xqt11a (trash1 ; trash2 ; Qn ; trash3)
  NB. berrQm=. trash0 qxt11a (Qm ; trash1 ; trash2 ; trash3)
  lqt11a=: (normi@(<: upddiag)@(mp  ct) % FP_EPS * #)@(2 {:: ])
  qlt11a=: (norm1@(<: upddiag)@(mp~ ct) % FP_EPS * #)@(0 {:: ])
  qrt11a=: (norm1@(<: upddiag)@(mp~ ct) % FP_EPS * #)@(0 {:: ])
  rqt11a=: (normi@(<: upddiag)@(mp  ct) % FP_EPS * #)@(2 {:: ])

  ('dgeqp3_mttmp_' tmonad (((; 0 #~ c)@(0&{::))`(<:@(1&{::) ; 0&{:: , 2&{::)`(rcond"_)`(_."_)`(prt01  >. qrt11 ))) args
  ('dgeqp3_mttmp_' tmonad (((; 1 #~ c)@(0&{::))`(<:@(1&{::) ; 0&{:: , 2&{::)`(rcond"_)`(_."_)`(prt01  >. qrt11 ))) args
  ('zgeqp3_mttmp_' tmonad (((; 0 #~ c)@(0&{::))`(<:@(1&{::) ; 0&{:: , 2&{::)`(rcond"_)`(_."_)`(prt01  >. qrt11 ))) args
  ('zgeqp3_mttmp_' tmonad (((; 1 #~ c)@(0&{::))`(<:@(1&{::) ; 0&{:: , 2&{::)`(rcond"_)`(_."_)`(prt01  >. qrt11 ))) args

  ('gelpf'         tmonad ((            0&{:: )`]                      `(rcond"_)`(_."_)`(lpt01a >. lqt11a))) args
  ('geplf'         tmonad ((            0&{:: )`]                      `(rcond"_)`(_."_)`(plt01a >. qlt11a))) args
  ('geprf'         tmonad ((            0&{:: )`]                      `(rcond"_)`(_."_)`(prt01a >. qrt11a))) args
  ('gerpf'         tmonad ((            0&{:: )`]                      `(rcond"_)`(_."_)`(rpt01a >. rqt11a))) args

  coerase < 'mttmp'
  erase 'lpt01a plt01a prt01a rpt01a lqt11a qlt11a qrt11a'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testpf
NB.
NB. Description:
NB.   Adv. to make verb to test gexxx by matrix of generator
NB.   and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testpf
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
NB.     ?@$&0 testpf_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testpf_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testpf_mt_ 150 200

testpf=: 1 : 'EMPTY [ testgepf_mt_@u'
