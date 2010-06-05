NB. expm.ijs
NB. Matrix exponent and Cauchy integral from state-space
NB. representation of LTI system
NB.
NB. prexpm  prepare time-invariant parts for expm
NB. expm    calculate matrix exponent and Cauchy integral
NB.
NB. References:
NB. [1] Podchukaev V.A. Theory of informational processes
NB.     and systems. - M., 2006. (Подчукаев В. А. Теория
NB.     информационных процессов и систем. - М.: Гардарики,
NB.     2006 - 209 с.)
NB.     URL: http://www.sgau.ru/uit/Book3.htm
NB. [2] Andrievskiy B.R., Fradkov A.L. Selected chapters of
NB.     automatic control theory with MATLAB examples. -
NB.     SPb., 2000 (Андриевский Б.Р., Фрадков А. Л.
NB.     Избранные главы теории автоматического управления
NB.     с примерами на языке MATLAB. - СПб.: Наука, 2000. -
NB.     475 с., ил. 86)
NB.
NB. Resources:
NB. - http://www.jsoftware.com/jwiki/...
NB. - http://www.dvgu.ru/forum/...
NB.
NB. TODO:
NB. - gebal before geev or geevx instead of
NB.   gehrd before gees
NB.   see:
NB.   http://www.netlib.org/lapack/lug/node50.html
NB.   http://www.netlib.org/lapack/lug/node51.html
NB.   http://www.netlib.org/lapack/lug/node70.html
NB.   http://www.netlib.org/lapack/lug/node94.html
NB.
NB. 2008-02-29 1.0.0 Igor Zhuravlov |.'ur.ugvd.ciu@rogi'

script_z_ '~system/packages/math/mathutil.ijs'  NB. mp
script_z_ '~system/main/numeric.ijs'            NB. clean
script_z_ '~system/packages/math/makemat.ijs'   NB. idmat
require '~user/projects/lapack/lapack.ijs'      NB. '~addons/math/lapack/lapack.ijs'
require '~user/projects/lapack/geev.ijs'        NB. need_jlapack_ 'geev gesvd gesvx'
require '~user/projects/lapack/gesvd.ijs'       NB. (line above makes it excessive)
require '~user/projects/lapack/gesvx.ijs'       NB. -//-
require '~user/projects/tau/util.ijs'           NB. powsm rndmat rndmat_neig

coclass 'tau'

NB. =========================================================
NB. Local utilities

getCols=: {: @ $     NB. get columns count of table y
shybyx=: |.!.0 "0 1  NB. shift y with step x

NB. ---------------------------------------------------------
NB. makeGi
NB. Calculate G[i] for M[i]
NB.
NB. Syntax:
NB.   Gi=. makeGi mi , lambdai , Ng
NB. where:
NB.   mi      - lambdai multiplicity
NB.   lambdai - i-th eigenvalue of G
NB.   Ng      = #G , matrix G's minimal polynom's order
NB.   Gi      - mi-by-Ng table, G[i]
NB.
NB. Test:
NB.    makeGi 4 1j1 7
NB. 1 1 1 1  1  1   1
NB. 0 1 2 3  4  5   6
NB. 0 0 2 6 12 20  30
NB. 0 0 0 6 24 60 120

phi=: {: $ 1:                          NB. phi(lambdai) polynome coeffs
makeDphi=: p.. @ ] ^: [                NB. mi-by-Ng table of derivatives
shiftDphi=: (-@[) shybyx makeDphi      NB. shift each row by #row to the right
makeGi=: (i. @ {.) shiftDphi phi       NB. make coeffs table then shift

NB. ---------------------------------------------------------
NB. makeTi
NB. Calculate T[i] for M[i]
NB.
NB. Syntax:
NB.   Ti=. makeTi mi , lambdai , Ng
NB. where:
NB.   mi      - lambdai multiplicity
NB.   lambdai - i-th eigenvalue of G
NB.   Ng      = #G , matrix G's minimal polynom's order
NB.   Ti      - mi-by-Ng table, T[i]
NB.
NB. Test:
NB.    makeTi 4 1j1 7
NB. 1 1j1 0j2 _2j2   _4 _4j_4  0j_8
NB. 0   1 1j1  0j2 _2j2    _4 _4j_4
NB. 0   0   1  1j1  0j2  _2j2    _4
NB. 0   0   0    1  1j1   0j2  _2j2

lipows=: (^ i.)/ @ (_2 & {.)  NB. lambdai ^ i. Ng
IOnmi=: - @ i. @ (0 { ])      NB. 0 _1 _2 ... -(mi-1)
makeTi=: IOnmi shybyx lipows  NB. mi-by-Ng table from vector's shifts

NB. ---------------------------------------------------------
NB. prepV
NB. Prepare eigenvalues of matrix y: remove dups, then count.
NB. Outputs 3 columns:
NB. - mi, lambdai multiplicity
NB. - lambdai, i-th eigenvalue of G
NB. - Ng = #y , table y's minimal polynom's order
NB.
NB. If:
NB.   'vm vlambda vNg' =. |: prepV G
NB. then
NB.   Ng  -: (+/ vm)
NB.   vNg -: (# vNg) $ Ng
NB.
NB. Test:
NB.    prepV diagmat 4 4 3 2j2 2j_2 1j1 1j1
NB. 2    4 7
NB. 1    3 7
NB. 1  2j2 7
NB. 1 2j_2 7
NB. 2  1j1 7

prepV=: ((1: #. =) ,. ~. ,. #) @ (2 & geev_jlapack_)

NB. ---------------------------------------------------------
NB. makeG
NB. Build square matrix G of augmented LTI system: G = ( A  B )
NB.                                                    ( 0  0 )
NB. Syntax:
NB.   G=. makeG A;B[;trash]
NB. where
NB.   A   - Nx-by-Nx table, state matrix of LTI system
NB.   B   - Nx-by-Nu table, control input matrix of LTI system
NB.   G   - Ng-by-Ng table, augmented LTI system
NB.   Ng  = Nx+Nu
NB.   Nx >= 0
NB.   Nu >= 0

makeG=: (0 {:: ]) ((+ & getCols) {. ,.) (1 {:: ])

NB. ---------------------------------------------------------
NB. makeLtM
NB. Calculate matrix M or vector L(t) for equation M*A(t)=L(t)
NB.
NB. Syntax:
NB.   Lt=. ts makeLtM V
NB.   M=. 0 makeLtM V
NB. where:
NB.   V  - (#vm)-by-3 table, prepared eigenvalues of G, output of prepV
NB.   ts > 0, sample time
NB.   M  - Ng-by-Ng table for equation M*A(t)=L(t)
NB.   Lt - Ng-vector, RHS L(t) for equation M*A(t)=L(t)
NB.   Ng = #G , table G's minimal polynom's order
NB.
NB. Tests:
NB.    1.0 makeLtM prepV diagmat 4 4 3 2j2 2j_2 1j1 1j1
NB. 54.5982 54.5982 20.0855 _3.07493j6.71885 _3.07493j_6.71885 1.46869j2.28736 1.46869j2.28736
NB.    0 makeLtM prepV diagmat 4 4 3 2j2 2j_2 1j1 1j1
NB. 1    4   16      64  256      1024    4096
NB. 0    1    8      48  256      1280    6144
NB. 1    3    9      27   81       243     729
NB. 1  2j2  0j8  _16j16  _64 _128j_128  0j_512
NB. 1 2j_2 0j_8 _16j_16  _64  _128j128   0j512
NB. 1  1j1  0j2    _2j2   _4     _4j_4    0j_8
NB. 0    1  2j2     0j6 _8j8       _20 _24j_24

tspows=: (^ i.) {.                    NB. ts ^ i. mi
prepMi=: (makeGi * makeTi) @ ]        NB. prepare Mi=. Gi*Ti
prepLti=: (^ @ (* (1 { ]))) * tspows  NB. prepare L[i](t)
prepLtMi=: prepLti`prepMi @. (0 = [)  NB. choose L[i](t) or M[i] to prepare
makeLtM=: ; @: (< @ prepLtMi " 1)     NB. make and merge all L[i](t) or M[i]

NB. ---------------------------------------------------------
NB. makeAt
NB. Solve equation M*A(t)=L(t) for A(t)
NB.
NB. Syntax:
NB.   At=. NxPMV makeAt ts
NB. where:
NB.   ts    > 0, sample time
NB.   NxPMV - output of prexpm, being (Nx;P;M;V)
NB.   At    - Ng-vector, solution A(t) of equation M*A(t)=L(t)
NB.   Ng    = #G , table G's minimal polynom's order

makeLt=: (makeLtM (3 & {::))~
makeAt=: gesvx_jlapack_ @ ((2 {:: [) ; makeLt)

NB. =========================================================
NB. prexpm
NB. Prepare time-invariant parts for expm
NB.
NB. Syntax:
NB.   'Nx P M V'=. prexpm A;B[;trash]
NB. where:
NB.   A - Nx-by-Nx table, state matrix of LTI system, should
NB.       be stable, i.e. all eigenvalues of A must have
NB.       negative real parts
NB.   B - Nx-by-Nu table, control input matrix of LTI system
NB.   V - (#vm)-by-3 table, prepared eigenvalues of G, output
NB.       of prepV
NB.   M - Ng-by-Ng table, matrix for equation M*A(t)=L(t),
NB.       output of makeLtM
NB.   P - Ng-by-Ng-by-Ng table, powers 0..(Ng-1) of G,
NB.       output of powsm
NB.   G - Ng-by-Ng table, augmented LTI system, output of
NB.       makeG
NB.   Ng = Nx + Nu
NB.   Nx >= 0
NB.   Nu >= 0

prexpm=: getCols @ (0 {:: ]) ; (powsm ; (0 & makeLtM ; ]) @ prepV) @ makeG

NB. ---------------------------------------------------------
NB. expm
NB. Calculate matrix exponent
NB.   Phi(ts) = exp(A*ts)
NB. and Cauchy integral
NB.   Gamma(ts) = Integral(exp(A*(ts-t))*dt,t=0..ts)*B
NB. for sampe time ts via Lagrange-Sylvester interpolation
NB. polynome [1, p. 88]
NB.
NB. Syntax:
NB.   'Phi Gamma'=. NxPMV expm ts
NB. where:
NB.   NxPMV - output of prexpm, being (Nx;P;M;V)
NB.   ts    > 0, sample time
NB.   Phi   - Nx-by-Nx table, matrix exponent
NB.   Gamma - Nx-by-Nu table, Cauchy intergal
NB.   Nu    = (#P)-Nx
NB.
NB. Applications:
NB.   'Phi Gamma'=. (A;B) ((prexpm @ [) expm ]) ts

splitbyx=: {."1 ; }."1              NB. split table y at column x
extract=: [ splitbyx {.             NB. extract 1st x rows from table y
augsys=: +/ @ ((1 {:: [) * makeAt)  NB. find exponent of augmented system
expm=: (0 {:: [) extract augsys     NB. extract Phi and Gamma from exp(G*ts) [2, p. 145]

NB. =========================================================
NB. Test suite

NB. Syntax: is_passed=. texpm A;B;ts

texpm=: 3 : 0
'A B ts'=. y
P=. 0 {:: (prexpm A;B) expm ts
eigA=. /:~ ^ ts * (2 geev_jlapack_ A)
eigP=. /:~ (2 geev_jlapack_ P)
err=. clean %: +/ *: | eigA - eigP
0 = err
)

NB. Syntax: testexpm ''

testexpm=: 3 : 0
'A0 A1 A2 A3'=. (rndmat_neig &. >) 4;6;8;10
'B0 B1 B2 B3'=. (rndmat &. >) 4 3;6 5;8 7;10 9
'ts0 ts1 ts2 ts3'=. 0.01 * >: ? 4 $ 100
texpm &> (< A0;B0;ts0) , (< A1;B1;ts1) , (< A2;B2;ts2) , (< A3;B3;ts3)
)
