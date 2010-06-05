NB. lyapunov.ijs
NB. Solve Luapunov equation
NB.
NB. lyapchol   solve continuous-time Lyapunov equation
NB.            directly for Cholesky factor
NB. dlyapchol  solve discrete-time Lyapunov equation
NB.            directly for Cholesky factor
NB. lyap       solve continuous-time Lyapunov equation
NB. dlyap      solve discrete-time Lyapunov equation
NB.
NB. TODO
NB. - replace some of assert. by throw.
NB. - gebal before geev or geevx instead of
NB.   gehrd before gees
NB.   see:
NB.   http://www.netlib.org/lapack/lug/node50.html
NB.   http://www.netlib.org/lapack/lug/node51.html
NB.   http://www.netlib.org/lapack/lug/node70.html
NB.   http://www.netlib.org/lapack/lug/node94.html
NB.
NB. 2008-03-30 1.0.0 Igor Zhuravlov |.'ur.ugvd.ciu@rogi'

script_z_ '~system/packages/math/matutil.ijs'  NB. diag
require '~user/projects/lapack/lapack.ijs'      NB. -> require '~addons/math/lapack/lapack.ijs'
NB. need_jlapack_ 'gees geev gerqf potrf trtrs'
require '~user/projects/lapack/gees.ijs'        NB. -> /dev/null
require '~user/projects/lapack/geev.ijs'        NB. -> /dev/null
require '~user/projects/lapack/gerqf.ijs'       NB. -> /dev/null
require '~user/projects/lapack/potrf.ijs'       NB. -> /dev/null
require '~user/projects/lapack/trtrs.ijs'       NB. -> /dev/null
require '~user/projects/tau/util.ijs'           NB. shiftdiag rndmat rndmat_neig

coclass 'tau'

NB. =========================================================
NB. Utilities

split=: (}:;{:) &. >                    NB. split under box at last item
N2=: [: %: [: +/ [: *: |                NB. 2-norm of vector

NB. ---------------------------------------------------------
NB. sorzhouiter
NB. Execute single iteration # (n-j) of Sorensen-Zhou algorithm
NB.
NB. Syntax:
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
    atmp=. (+ lambda) shiftdiag (jj {. x)        NB. conj(lambda)*I + R1
    btmp=. (B1 mp (+ - bh)) + (ijj { x) * - tau  NB. -tau*r-(1/tau)*B1*b
    u=. trtrs_jlapack_ atmp ; btmp               NB. solve for u complex upper triangular system atmp*u=btmp
  else.
    u=. j $ 0
  end.

  NB. update B1 if needed
  NB. replace in U leading j elements in j-th column if needed
  NB. replace in U (j+1)-th diagonal element
  ((- & (u */ bh))^:chk B1) ; ij ; (tau (< jj) } u (ijj })^:chk 2 {:: y)
)

NB. =========================================================
NB. lyapchol
NB. Solve stable non-negative definite continuous-time
NB. Lyapunov equation:
NB.   A*X + X*A' + B*B' = 0
NB. directly for Cholesky factor U, X = U*U'. The matrix A
NB. must be stable, that is, all the eigenvalues of A must
NB. have negative real parts.
NB.
NB. Usage:
NB.   U=. A lyapchol B
NB.   U=. (R;Q) lyapchol B
NB. where:
NB.   A   - N-by-N stable matrix
NB.   R,Q - N-by-N matrices from non-real Schur factorization:
NB.         Q*R*Q' = A
NB.   B   - N-by-M matrix
NB.   U   - N-by-N upper triangular matrix
NB.   N  >= 0
NB.   M  >= 0
NB.
NB. Reference:
NB.   Solution of stable continuous- or discrete-time Lyapunov equations
NB.   (Cholesky factor), URL: http://www.slicot.org/shared/doc/SB03OD.html

lyapchol=: 4 : 0
  vmatrixorvector_jlapack_ y
  n=. # y
  if. L. x do.
    'Q R'=. x                               NB. Schur factorization Q*R*Q' = A
    vsquare_jlapack_ Q
    assert. n = # Q
    assert. Q (-: &: $) R                   NB. R and Q shapes are match
  else.
    vsquare_jlapack_ x
    assert. n = # x
    'Q R'=. 13 gees_jlapack_ x              NB. Schur factorization Q*R*Q' = A
  end.
  assert. -. 0 e. 0 > 9 o. diag R           NB. A is stable
  R1=. 8 gerqf_jlapack_ y                   NB. RQ factorization R1*P1=B
  R2=. 8 gerqf_jlapack_ (h Q) mp R1         NB. RQ factorization R2*P2=Q'*R1

  NB. solve triangular Luapunov equation directly for Cholesky factor V
  V=. 2 {:: R (sorzhouiter ^: n) R2 ; (i. n) ; ((2 $ n) $ 0)

  R3=. 8 gerqf_jlapack_ Q mp V              NB. RQ factorization R3*P3=Q*V
  R3=. (*"1 * @: diag) R3                   NB. negate columns having negative diagonal elements
NB. ---  smoutput 2 6 $ 'Q' ; 'R' ; 'R1' ; 'R2' ; 'V' ; 'R3' ; Q ; R ; R1 ; R2 ; V ; R3
  R3
)

NB. ---------------------------------------------------------
NB. dlyapchol
NB. Solve convergent non-negative definite discrete-time
NB. Lyapunov equation
NB.   A*X*A' + X + B*B' = 0
NB. directly for Cholesky factor U, X = U*U'. The matrix A
NB. must be d-stable, that is, all the eigenvalues of A must
NB. lie inside the unit circle.

NB. ---------------------------------------------------------
NB. lyap
NB. Solve continuous-time Lyapunov equation
NB.   A*X + X*A' + B*B' = 0

NB. ---------------------------------------------------------
NB. dlyap
NB. Solve discrete-time Lyapunov equation
NB.   A*X*A' + X + B*B' = 0

NB. =========================================================
NB. Test suite

NB. Syntax: is_passed=. tlyapchol A;B

tlyapchol=: 3 : 0
'A B'=. y
U=. A lyapchol B
X=. (mp h) U
err=. clean (A mp X) + (X (mp h) A) + ((mp h) B)
(*./ ^: 2) 0 = err
)

NB. Syntax: testlyapchol ''

testlyapchol=: 3 : 0
ma0=. 0 0 $ 0
mb0=. 0 0 $ 0
ma1=. rndmat_neig 4
mb1=. ? 4 4 $ 10
ma2=. rndmat_neig 10
mb2=. ? 10 5 $ 10
ma3=. rndmat_neig 100
mb3=. ? 100 50 $ 10
tlyapchol &> (<ma0;mb0) , (<ma1;mb1) , (<ma2;mb2) , (<ma3;mb3)
)
