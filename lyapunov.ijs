NB. lyapunov.ijs
NB. Solve Luapunov equation
NB.
NB. lyapchol  - solve continuous-time Lyapunov equation
NB.             directly for Cholesky factor
NB. dlyapchol - solve discrete-time Lyapunov equation
NB.             directly for Cholesky factor
NB. lyap      - solve continuous-time Lyapunov equation
NB. dlyap     - solve discrete-time Lyapunov equation
NB.
NB. TODO
NB. - replace some of assert. by throw.
NB. - handle non-squared case
NB. - implement enhanced LAPACK interfaces
NB. - CO2DI: Cayley transform
NB.
NB. 2007-12-25 1.0.0 Igor Zhuravlov |.'ur.ugvd.ciu@rogi'

require '~system/packages/math/matutil.ijs'     NB. diag
require '~user/projects/lapack/lapack.ijs'      NB. -> require '~addons/math/lapack/lapack.ijs'
NB. need_jlapack_ 'gees geev gerqf potrf trtrs'
require '~user/projects/lapack/gees.ijs'        NB. -> /dev/null
require '~user/projects/lapack/geev.ijs'        NB. -> /dev/null
require '~user/projects/lapack/gerqf.ijs'       NB. -> /dev/null
require '~user/projects/lapack/potrf.ijs'       NB. -> /dev/null
require '~user/projects/lapack/trtrs.ijs'       NB. -> /dev/null

NB. ===========================================================================
NB. lyapchol
NB. Solve stable non-negative definite continuous-time Lyapunov equation
NB. A*X + X*A' + B*B' = 0 directly for Cholesky factor U, X = U*U'
NB.
NB. Usage:
NB.   U=. A lyapchol B
NB.   U=. (R;Q) lyapchol B
NB. where:
NB.   A   - N-by-N stable matrix, i.e. all eigenvalues of A must
NB.         have negative real parts
NB.   R,Q - N-by-N matrices from Schur factorization: Q*R*Q' = A
NB.   B   - N-by-M matrix
NB.   U   - N-by-N upper triangular matrix
NB.   N   > 0
NB.   M   > 0
NB.
NB. Solution must exist and be unique, i.e. vectors eigenvalue(A) and
NB.   eigenvalue(B) have not to have common values
NB.
NB. Reference:
NB.   Solution of stable continuous- or discrete-time Lyapunov equations
NB.   (Cholesky factor), URL: http://www.slicot.org/shared/doc/SB03OD.html

lyapchol=: 4 : 0
  n=. # y
  if. 0 = n do. i. 0 return. end.           NB. early termination
  vmatrixorvector_jlapack_ y
  if. L. x do.
    'Q R'=. x                               NB. Schur factorization Q*R*Q' = A
    vsquare_jlapack_ Q
    assert. n = # Q
    assert. Q (-: &: $) R                   NB. R and Q shapes are match
  else.
    vsquare_jlapack_ x
    assert. n = # x
    'Q R'=. 2 {. zgees_jlapack_ x           NB. Schur factorization Q*R*Q' = A
  end.
  evA=. diag R                              NB. eigenvalues(A)
  evB=. 1 {:: zgeev_jlapack_ y              NB. eigenvalues(B)
  assert. -. 0 e. 0 > 9 o. evA              NB. A is stable
  assert. evA (*:/ @: (e. , e.~)) evB       NB. solution exists and is unique
  R1=. zgerqfR y                            NB. RQ factorization R1*P1=B
  R2=. zgerqfR (h Q) mp R1                  NB. RQ factorization R2*P2=Q'*R1

  NB. solve triangular Luapunov equation directly for Cholesky factor V
  V=. 2 {:: R (sorzhouiter ^: n) R2 ; (i. n) ; ((2 $ n) $ 0)

  R3=. zgerqfR Q mp V                       NB. RQ factorization R3*P3=Q*V
  R3=. (*"1 * @: diag) R3                   NB. negate columns having negative diagonal elements
  smoutput 2 6 $ 'Q' ; 'R' ; 'R1' ; 'R2' ; 'V' ; 'R3' ; Q ; R ; R1 ; R2 ; V ; R3
  R3
)

NB. ===========================================================================
NB. sorzhouiter
NB. Execute single iteration # (n-j) of Sorensen-Zhou algorithm
NB.
NB. Usage:
NB.   'B1 ijupd Uupd'=. R sorzhouiter B ; ij ; U
NB. where:
NB.   R     - N-by-N upper triangular stable matrix, i.e. all eigenvalues of R
NB.           must have negative real parts
NB.   B     - N-by-M matrix, updated at step #j
NB.   ij    - i. j
NB.   U     - N-by-N upper triangular matrix with all but first j columns
NB.           updated
NB.   B1    - (N-1)-by-M matrix B without last row and updated after step #j
NB.   ijupd - }: ij
NB.   Uupd  - N-by-N matrix U with updated column #j
NB.   N     > 0
NB.   M     > 0
NB.
NB. Reference:
NB.   Danny C. Sorensen and Yunkai Zhou, "Direct methods for matrix
NB.   Sylvester and Lyapunov equations", J. Appl. Math, vol. 2003, no. 6,
NB.   pp. 277-303, 2003. doi:10.1155/S1110757X03212055

sorzhouiter=: 4 : 0
  'B1 bh ij j'=. ; split }: y   NB. split B on B1;bh at last row and ij on ijupd;j at last atom
  jj=. 2 $ j                    NB. IO: lambda in R, tau in U, R1 in R
  ijj=. <ij;j                   NB. IO: r in R, u in U
  lambda=. (< jj) { x
  tau=. (N2 bh) % (%: _2 * 9 o. lambda)
  bh=. bh % tau
  chk=. (0 < # B1) *. (0 < tau)
  if. chk do.
    atmp=. (jj {. x) shiftR1 } (lambda ; ij)     NB. "R1+conj(lambda)*I" where R1 is "jj {. x"
    btmp=. (B1 mp (+ - bh)) + (ijj { x) * - tau  NB. -tau*r-(1/tau)*B1*b
    u=. ztrtrs_jlapack_ atmp ; btmp              NB. solve for u complex upper triangular system atmp*u=btmp
  else.
    u=. j $ 0
  end.

  NB. update B1 if needed
  NB. replace in U leading j elements in j-th column if needed
  NB. replace (j+1)-th diagonal element in U
  ((- & (u */ bh))^:chk B1) ; ij ; (tau (< jj) } u (ijj })^:chk 2 {:: y)
)

NB. ===========================================================================
NB. dlyapchol - solve convergent non-negative definite discrete-time Lyapunov
NB.             equation A*X*A' + X + B*B' = 0 directly for Cholesky factor

NB. ===========================================================================
NB. lyap      - solve continuous-time Lyapunov equation A*X + X*A' + B*B' = 0

NB. ===========================================================================
NB. dlyap     - solve discrete-time Lyapunov equation A*X*A' + X + B*B' = 0

NB. ===========================================================================
NB. utilities

split=: (}:;{:) &. >                    NB. split under box at last item
N2=: [: %: [: +/ [: *: |                NB. 2-norm of vector
h=: +@|:                                NB. conjugate transpose of table
zgerqfR=: 0 {:: zgerqf_jlapack_         NB. take R from RQ factorization
shiftdiag=: (diag @: [) + (+ @ ])       NB. diag(x)+conjugate(y)
iosdiag=: <"1 @ ,.~                     NB. indices of y-th diagonal elements
shiftR1=: (shiftdiag (0 & {::))`(iosdiag @: (1 & {::) @: ])`[

NB. ===========================================================================
NB. test matrices

testlyapchol=: 3 : 0
  A=. 4 4 $ _1 _1 2 2 37 _10 _4 _2 _12 0 7 7 _12 4 _6 _9
  B=. 4 4 $ 2 2 4 5 0 8 7 7 0 0 2 6 0 0 0 2
  X=. 4 4 $ 18.3023846590474 43.8883259794941 75.7212100367117 _56.8758547174409 43.8883259794941 139.292151950139 128.468952844055 _100.964634818162 75.7212100367117 128.468952844055 492.482132895331 _365.531487118111 _56.8758547174409 _100.964634818162 _365.531487118111 274.871182227256
  U=. (] ;.0) dpotrf_jlapack_ (] ;.0) X        NB. upper triangular Cholesky factor U, X=U*U'
  u=. A lyapchol B
  x=. u mp h u
  smoutput 2 10 $ 'A' ; 'B' ; 'X' ; 'U' ; 'X-UU''' ; 'AX+XA''+BB''' ; 'solution: u' ; 'x=u*u''' ; 'Ax+xA''+BB''' ; 'U-u' ; A ; B ; X ; U ; (X - U mp h U) ; ((A mp X) + (X mp h A) + (B mp h B)) ; u ; x ; ((A mp x) + (x mp h A) + (B mp h B)) ; (U-u)
  U (-: !. 1e_11) u
)
