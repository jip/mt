NB. hrd.ijs
NB. Reduce a general matrix to upper Hessenberg form
NB.
NB. gehrdu  Reduce a general matrix to upper Hessenberg form
NB. gehd2u  Reduce a general matrix to upper Hessenberg form
NB.         (non-blocked version)
NB.
NB. hehrd   Reduce a general matrix to Hessenberg form
NB.         (3-diagonal Hermitian matrix)
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. Differences between cIOSs at consequent iterations:
NB. cios(i+1)-cios(i)
GEHD2UDCIOS=: 5 2 $ 1j_1 1 1j_1 1 0 1 0 1j_1 1j_1 1j_1

NB. ---------------------------------------------------------
NB. mkcios0gehd2u
NB. Create cIOS for gehd2u at 0-th iteration
NB.
NB. Syntax:
NB.   cios0=. mkcios0gehd2u (fs , n)
NB. where
NB.   fs    - 2-vector of integers (f,s) 'from' and 'size',
NB.           defines submatrix to be reduced position in
NB.           matrix A, see gebalp
NB.   n     ≥ 0, integer, size of matrix A
NB.   cios0 - 5×2-table cios[0], cIOSs corresponding to
NB.           iteration 0, see gehd2ustep

mkcios0gehd2u=: 3 : 0
  'f1 s1 e'=. (1 _1&+ , (+/)) @ }: 'f s n'=. y
  5 2 $ (f1 , f , f1 , f , e , f , 0 , 3 $ f1) j. (s1 , 1 , s , 1 1 1 , e , s1 , s1 , (n-f1))
)

NB. ---------------------------------------------------------
NB. gehd2ustep
NB. Single step of Hessenberg reduction
NB.
NB. Syntax:
NB.   'Ai1 ciosi1'=. dcios gehd2ustep (Ai ; ciosi)
NB. where
NB.   Ai    - n×n-matrix A(i) with isolated eigenvalues and
NB.           submatrix A11 to update (see gebal) before i-th
NB.           iteration (i=f:f+s-2)
NB.   ciosi - 5×2-matrix cios(i) of cIOSs (see struct.ijs)
NB.           for i-th iteration; rows (0:4) contains (let
NB.           e=f+s):
NB.             0 - (e-(i+1))-vector Y=(α[i],x[i][1:e-(i+2)])
NB.                 to reflect, is stored in A[i+1:e-1,i],
NB.                 cIOS are:
NB.                   ((i+1) j. (e-(i+1))) , (i j. 1)
NB.             1 - (e-i)-vector Z=(β[i],v[i][1:e-(i+2)],τ[i])
NB.                 of reflection result, is stored in
NB.                 A[i+1:e,i], cIOS are:
NB.                   ((i+1) j. (e-i)) , (i j. 1)
NB.             2 - scalar τ[i], is stored in A[e,i], cIOS
NB.                 is:
NB.                   (e j. 1) , (i j. 1)
NB.             3 - e×(e-(i+1))-matrix R to apply an
NB.                 elementary reflector from the right, is
NB.                 stored in A[0:e-1,i+1:e-1], cIOS are:
NB.                   (0 j. e) , ((i+1) j. (e-(i+1)))
NB.             4 - (e-(i+1))×(n-(i+1))-matrix L to apply an
NB.                 elementary reflector from the left, is
NB.                 stored in A[i+1:e-1,i+1:n-1], cIOS are:
NB.                   ((i+1) j. (e-(i+1))) , ((i+1) j. (n-(i+1)))
NB.   dcios - difference between cIOSs at consequent
NB.           iterations: cios(i+1)-cios(i)
NB.   Ai1   - n×n-matrix A(i+1) after i-th iteration
NB.   ciosi - 5×2-matrix cios(i+1) of cIOSs for (i+1)-th
NB.           iteration

gehd2ustep=: (< @ (+ (1&{::))) 1} (((1 {:: ]) (< @ ((larfg`(1 3 4 larfRL)) gerf0 0 1 2)) (0 {:: ])) 0} ])

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. gehd2u
NB. Reduce a general matrix A to upper Hessenberg form H by
NB. an unitary similarity transformation:
NB.   Q' * A * Q = H
NB. This is a non-blocked version of algorithm.
NB.
NB. Syntax:
NB.   B=. gehd2u (A ; fs)
NB. where
NB.   A   - n×n-matrix with isolated eigenvalues (see gebalp)
NB.   fs  - 2-vector of integers (f,s) 'from' and 'size',
NB.         defines submatrix A11 position in matrix A (see
NB.         gebalp)
NB.   B   - max(n,f+s+1)×n-matrix with packed H and Qf (see
NB.         storage layout below)
NB.   A11 - s×s-matrix to be reduced
NB.   H   - s×s-matrix, upper Hessenberg form of A11
NB.   Qf  - s×(s-1)-matrix, unit lower triangular, represents
NB.         unitary (orthogonal) matrix Q
NB.
NB. Storage layout:
NB.   A11 - stored in A[f:f+s-1,f:f+s-1]
NB.   H   - stored in place of A11 in upper triangular part
NB.         and in 1st subdiagonal
NB.   Qf  - stored below 1st subdiagonal in A[f+2:f+s,f:f+s-2],
NB.         represents unitary (orthogonal) matrix Q in
NB.         factored form:
NB.           Q = Q[f] * Q[f+1] * ... * Q[f+s-2] ,
NB.         where each elementary reflector Q[i] is
NB.         represented as:
NB.           Q[i]=I-τ[i]*v[i]*v[i]' ,
NB.         and values v[i][i-f+1:s-2] and τ[i] are
NB.         stored in (i-f)-th column of Qf:
NB.           Qf[0:s-1,i-f]=(0,...,0,1,v[i][i-f+1],...,v[i][s-2],τ[i])
NB.         assuming v[i][0:i-f-1]=0, v[i][i-f]=1,
NB.         v[i][i-f+1:s-2] is stored in A[i+2:f+s-1,i], τ[i]
NB.         is stored in A[f+s,i]
NB.   Example for f=1, s=5, n=7:
NB.     input                        output
NB.     (  a  a  a  a  a  a  a  )    (  a  a  h  h  h  h  a  )
NB.     (     a  a  a  a  a  a  )    (     a  h  h  h  h  a  )
NB.     (     a  a  a  a  a  a  )    (     β1 h  h  h  h  h  )
NB.     (     a  a  a  a  a  a  )    (     v1 β2 h  h  h  h  )
NB.     (     a  a  a  a  a  a  )    (     v1 v2 β3 h  h  h  )
NB.     (     a  a  a  a  a  a  )    (     v1 v2 v3 β4 h  h  )
NB.     (                    a  )    (     τ1 τ2 τ3 τ4    a  )
NB.     if (f+s=n) then zero row will be appended to matrix to store τ[f:f+s-2]
NB.
NB. If:
NB.   B=. gehd2u (A ; fs)
NB.   H=. _1 utri B
NB.   Q=. unghr (B ; fs)
NB. then
NB.   (idmat n) -: ((ct Q) mp Q)
NB.   H -: (ct Q) mp A mp Q
NB.
NB. Notes:
NB. - emulates xGEHD2 up to storage layout scheme
NB.
NB. References:
NB. [1] 
NB. [2] 

gehd2u=: 3 : 0
  n=. # 0 {:: y                 NB. A's size
  'f s'=. fs=. 1 {:: y          NB. extract 'from' and 'size' values of A11 within A

  if. n = (f+s) do.
    y=. ((,&0 &. >) {. y) 0} y  NB. append zero row to A under box to store τ[f:f+s-2]
  end.

  NB. Create cios[0] and replace in-place fs by it in y under box
  y=. (< mkcios0gehd2u (fs,n)) 1} y

  NB. do iterations, then extract B
  0 {:: GEHD2UDCIOS (gehd2ustep ^: (<: @ (s"_))) y
)

gehrdu=: gehd2u  NB. stub for a while

NB. =========================================================
NB. Test suite

tgehrd=: 3 : 0
)

thehrd=: 3 : 0
)

testhrd=: 3 : 0
)
