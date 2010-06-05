NB. expm.ijs
NB. Calculate matrix exponent and Cauchy integral from state and
NB. control input matrices of LTI system
NB.
NB. prexpm  prepare time-invariant parts for expm
NB. expm    calculate matrix exponent and Cauchy integral
NB.
NB. References:
NB. [1] Podchukaev V.A. Theory of informational processes and systems. - M.,
NB.     2006. (Подчукаев В. А. Теория информационных процессов и систем. -
NB.     М.: Гардарики, 2006 - 209 с.)
NB.     URL: http://www.sgau.ru/uit/Book3.htm
NB. [2] Andrievskiy B.R., Fradkov A.L. Selected chapters of automatic
NB.     control theory with MATLAB examples. - SPb., 2000 (Андриевский Б.
NB.     Р., Фрадков А. Л. Избранные главы теории автоматического управления
NB.     с примерами на языке MATLAB. - СПб.: Наука, 2000. - 475 с., ил. 86)
NB.
NB. Resources:
NB. - http://www.jsoftware.com/jwiki/...
NB. - http://www.dvgu.ru/forum/...
NB.
NB. 2008-02-15 1.0.0 Igor Zhuravlov |.'ur.ugvd.ciu@rogi'

require '~system/packages/math/mathutil.ijs'  NB. mp
require '~system/main/numeric.ijs'            NB. clean
require '~system/packages/math/makemat.ijs'   NB. idmat
require '~user/projects/lapack/lapack.ijs'    NB. '~addons/math/lapack/lapack.ijs'
require '~user/projects/lapack/geev.ijs'      NB. need_jlapack_ 'geev gesvx'
require '~user/projects/lapack/gesvx.ijs'     NB. (line above makes it excessive)
require '~user/projects/tau/util.ijs'         NB. makeP

coclass 'tau'

NB. =========================================================
NB. Utilities

getA=: 0 {:: ]       NB. extract A
getB=: 1 {:: ]       NB. extract B

getNx=: 0 {:: [      NB. extract Nx
getP=: 1 {:: [       NB. extract P
getM=: 2 {:: [       NB. extract M
getV=: 3 {:: [       NB. extract V

getIO=: 0 { ]        NB. extract IO[i]
getmi=: 1 { ]        NB. extract m[i]
getli=: 2 { ]        NB. extract lambda[i]
getNg=: 3 { ]        NB. extract Ng

getCols=: {: @ $     NB. get columns count of table y
shybyx=: |.!.0 "0 1  NB. shift y with step x

NB. ---------------------------------------------------------
NB. makeGi
NB. Calculate G[i] for M[i]
NB.
NB. Syntax:
NB.   Gi=. makeGi lambdai , mi , IO , Ng
NB. where:
NB.   Ng      = #G , matrix G minimal polynom's order
NB.   IO      - IO 1st row (atom) of corresp. M[i] (L[i](t)) in M (L(t))
NB.   mi      - lambdai multiplicity
NB.   lambdai - i-th eigenvalue of G
NB.   Gi      - mi-by-Ng matrix, G[i]
NB.
NB. Test:
NB.    makeGi 3 4 1j1 7
NB. 1 1 1 1  1  1   1
NB. 0 1 2 3  4  5   6
NB. 0 0 2 6 12 20  30
NB. 0 0 0 6 24 60 120

phi=: getNg $ 1:                   NB. phi(lambdai) polynome coeffs
makeDphi=: p.. @ ] ^: [            NB. mi-by-Ng table of derivatives
shiftDphi=: (-@[) shybyx makeDphi  NB. shift each row by #row to the right
makeGi=: (i.@getmi) shiftDphi phi  NB. make coeffs table then shift

NB. ---------------------------------------------------------
NB. makeTi
NB. Calculate T[i] for M[i]
NB.
NB. Syntax:
NB.   Ti=. makeTi lambdai , mi , IO , Ng
NB. where:
NB.   Ng      = #G , matrix G minimal polynom's order
NB.   IOs     - IO 1st row (atom) of corresp. M[i] (L[i](t)) in M (L(t))
NB.   mi      - lambdai multiplicity
NB.   lambdai - i-th eigenvalue of G
NB.   Ti      - mi-by-Ng matrix, T[i]
NB.
NB. Test:
NB.    makeTi 3 4 1j1 7
NB. 1 1j1 0j2 _2j2   _4 _4j_4  0j_8
NB. 0   1 1j1  0j2 _2j2    _4 _4j_4
NB. 0   0   1  1j1  0j2  _2j2    _4
NB. 0   0   0    1  1j1   0j2  _2j2

lipows=: getli (^ i.) getNg   NB. lambdai ^ i. Ng
IOnmi=: - @ i. @ getmi        NB. 0 _1 _2 ... -(mi-1)
makeTi=: IOnmi shybyx lipows  NB. mi-by-Ng table from vector's shifts

NB. ---------------------------------------------------------
NB. prepV
NB. Prepare eigenvalues from y: remove dups, then count.
NB. Outputs 4 columns:
NB. - IO, IO 1st row (atom) of corresp. M[i] (L[i](t)) in M (L(t))
NB. - mi, lambdai multiplicity
NB. - lambdai, i-th eigenvalue of G
NB. - Ng = #G , matrix G minimal polynom's order
NB.
NB. If:
NB.   'vIO vm vlambda vNg' =. |: prepV eigenvalues_of_G
NB. then
NB.   Ng  -: (+/ vm)
NB.   Ng  -: vm (+ & {:) vIO
NB.   vNg -: (# vNg) $ Ng
NB.
NB. Test:
NB.    prepV 4 4 3 2j2 2j_2 1j1 1j1
NB. 0 2    4 7
NB. 2 1    3 7
NB. 3 1  2j2 7
NB. 4 1 2j_2 7
NB. 5 2  1j1 7

cnt=: 1: #. =              NB. count atoms
IOs=: _1 & shybyx @ (+/\)  NB. IO 1st row (atom) of corresponding M[i] (L[i](t)) in M (L(t))
prepV=: (IOs ,. ]) @ cnt ,. ~. ,. #

NB. ---------------------------------------------------------
NB. makeG
NB. Build square matrix G of augmented LTI system: G = ( A  B )
NB.                                                    ( 0  0 )
NB. Syntax:
NB.   G=. makeG A;B
NB. where
NB.   A - Nx-by-Nx state matrix of LTI system
NB.   B - Nx-by-Nu control input matrix of LTI system
NB.   G - Ng-by-Ng matrix, augmented LTI system
NB.   Ng = Nx+Nu
NB.   Nx >= 0
NB.   Nu >= 0

makeG=: getA ((+ & getCols) {. ,.) getB

NB. ---------------------------------------------------------
NB. makeLtM
NB. Calculate matrix M or vector L(t) for equation M*A(t)=L(t)
NB.
NB. Syntax:
NB.   Lt=. T makeLtM V
NB.   M=. 0 makeLtM V
NB. where:
NB.   V  - (#vm)-by-4 matrix, prepared eigenvalues of G, output of prepV
NB.   T  - sampling period, T>0
NB.   M  - Ng-by-Ng matrix for equation M*A(t)=L(t)
NB.   Lt - Ng-vector, RHS L(t) for equation M*A(t)=L(t)
NB.   Ng = #G , matrix G minimal polynom's order
NB.
NB. Tests:
NB.    1.0 makeLtM prepV 4 4 3 2j2 2j_2 1j1 1j1
NB. 54.5982 54.5982 20.0855 _3.07493j6.71885 _3.07493j_6.71885 1.46869j2.28736 1.46869j2.28736
NB.
NB.    0 makeLtM prepV 4 4 3 2j2 2j_2 1j1 1j1
NB. 1    4   16      64  256      1024    4096
NB. 0    1    8      48  256      1280    6144
NB. 1    3    9      27   81       243     729
NB. 1  2j2  0j8  _16j16  _64 _128j_128  0j_512
NB. 1 2j_2 0j_8 _16j_16  _64  _128j128   0j512
NB. 1  1j1  0j2    _2j2   _4     _4j_4    0j_8
NB. 0    1  2j2     0j6 _8j8       _20 _24j_24

Tpows=: [ (^ i.) getmi                NB. T ^ i. mi
prepMi=: makeGi * makeTi              NB. prepare Mi=. Gi*Ti
prepLti=: (^ @ ([ * getli)) * Tpows   NB. prepare L[i](t)
prepLtMi=: prepLti`prepMi @. (0 = [)  NB. choose L[i](t) or M[i] to prepare
makeLtM=: ; @: (< @ prepLtMi " 1)     NB. make and merge all L[i](t) or M[i]

NB. ---------------------------------------------------------
NB. makeAt
NB. Solve equation M*A(t)=L(t) for A(t)
NB.
NB. Syntax:
NB.   At=. NxPMV makeAt T
NB. where:
NB.   T     - sampling period, T>0
NB.   NxPMV - output of prexpm, being (Nx;P;M;V)
NB.   At    - Ng-vector, solution A(t) of equation M*A(t)=L(t)
NB.   Ng    = #G , matrix G minimal polynom's order

makeAt=: gesvx_jlapack_ @ (getM ; ] makeLtM getV)

NB. =========================================================
NB. prexpm
NB. Prepare time-invariant parts for expm
NB.
NB. Syntax:
NB.   'Nx P M V'=. prexpm A;B
NB. where:
NB.   A - Nx-by-Nx state matrix of LTI system
NB.   B - Nx-by-Nu control input matrix of LTI system
NB.   V - (#vm)-by-4 matrix, prepared eigenvalues of G, output of prepV
NB.   M - Ng-by-Ng matrix for equation M*A(t)=L(t), output of makeLtM
NB.   P - Ng-by-Ng-by-Ng matrix, powers 0..(Ng-1) of G, output of makeP
NB.   G - Ng-by-Ng matrix, augmented LTI system, output of makeG
NB.   Ng = Nx + Nu
NB.   Nx >= 0
NB.   Nu >= 0

prexpm=: getCols @ getA ; (makeP ; (0 & makeLtM ; ]) @ prepV @ (2 & geev_jlapack_)) @ makeG

NB. ---------------------------------------------------------
NB. expm
NB. Calculate matrix exponent
NB.   Phi(T) = exp(A*T)
NB. and Cauchy integral
NB.   Gamma(T) = Integral(exp(A*(T-t))*dt,t=0..T)*B
NB. for sampling period T via Lagrange-Sylvester interpolation
NB. polynome [1, p. 88]
NB.
NB. Syntax:
NB.   'Phi Gamma'=. NxPMV expm T
NB. where:
NB.   NxPMV - output of prexpm, being (Nx;P;M;V)
NB.   T     - sampling period, T>0
NB.   Phi   - Nx-by-Nx matrix, matrix exponent
NB.   Gamma - Nx-by-Nu matrix, Cauchy intergal
NB.   Nu    = (#P)-Nx
NB.
NB. Applications:
NB.   'Phi Gamma'=. (A;B) ((prexpm @ [) expm ]) T

splitbyx=: {."1 ; }."1         NB. split table y at column x
extract=: [ splitbyx {.        NB. extract 1st x rows from table y
augsys=: +/ @ (getP * makeAt)  NB. find exponent of augmented system
expm=: getNx extract augsys    NB. extract Phi and Gamma from exp(G*T) [2, p. 145]

NB. =========================================================
NB. Test suite

NB. Syntax: testexpm ''

testexpm=: 3 : 0
A=. 4 4 $ 9 7 6 8 _5 5 _1 5 2 _1 2 _6 0 _9 6 _2
B=. 4 3 $ 8 _8 _7 8 _3 _6 _2 6 _6 7 _9 2
T0=. 1.0
Phi0=. 4 4 $ 4583.5214041683 3163.7993176507 5728.9541875919 1799.6339167756 _3146.5211228348 _2110.7848405469 _3873.3554848169 _1230.7932384536 _441.7654701830 _225.1134698769 _478.0047389178 _167.2890633915 2539.3940584840 1703.7459778865 3126.7462652618 992.0945826301
Gamma0=. 4 3 $ 7707.3287490973 _3474.2601844433 _9957.9575055456 _5156.1280004707 2378.5356749258 6632.0542895394 _575.5472194441 330.9423568862 700.3399369133 4161.2894243602 _1917.7014396919 _5354.9631406784
T1=. 0.35
Phi1=. 4 4 $ 20.9001489724 9.8494136392 22.5809734768 8.1019006385 _13.1348881505 _4.8340120112 _11.3683336856 _4.9325702778 0.3194420921 3.2026954449 2.1278681154 0.5226543224 10.7428054074 3.2534025545 9.1133995995 2.9041972626
Gamma1=. 4 3 $ 27.9959185456 _16.6659816764 _32.3553687499 _11.4394990436 8.1965846460 14.2018691109 3.1376033523 0.9624659161 _6.5665088754 7.3544330863 _4.9172188143 _10.8119095000

smoutput (Phi0;Gamma0) (- clean) &. > ((prexpm A;B) expm T0)
smoutput (Phi1;Gamma1) (- clean) &. > ((prexpm A;B) expm T1)
)
