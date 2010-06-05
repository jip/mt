NB. expm2.ijs
NB. catch error

script_z_ '~system/packages/math/mathutil.ijs'  NB. mp
script_z_ '~system/main/numeric.ijs'            NB. clean
script_z_ '~system/packages/math/makemat.ijs'   NB. idmat
require '~addons/math/lapack/lapack.ijs'
need_jlapack_ 'dgeev zgeev dgesv zgesv'
require '~user/projects/tau/util.ijs'           NB. makeP

coclass 'tau'

NB. =========================================================
NB. Utilities

getCols=: {: @ $     NB. get columns count of table y
shybyx=: |.!.0 "0 1  NB. shift y with step x
dumper=: 0&$: :(4 : 0)
  if. x=1 do.
    'M Lt'=: y
  elseif. x=2 do.
    At=: y
  end.
  smoutput 'dumper(): y='
  smoutput y
  y
)

NB. ---------------------------------------------------------
NB. makeGi
NB. Calculate G[i] for M[i]
NB.
NB. Syntax:
NB.   Gi=. makeGi mi , lambdai , Ng
NB. where:
NB.   mi      - lambdai multiplicity
NB.   lambdai - i-th eigenvalue of G
NB.   Ng      = #G , matrix G minimal polynom's order
NB.   Gi      - mi-by-Ng matrix, G[i]
NB.
NB. Test:
NB.    makeGi 4 1j1 7
NB. 1 1 1 1  1  1   1
NB. 0 1 2 3  4  5   6
NB. 0 0 2 6 12 20  30
NB. 0 0 0 6 24 60 120

phi=: (2 { ]) $ 1:                     NB. phi(lambdai) polynome coeffs
makeDphi=: p.. @ ] ^: [                NB. mi-by-Ng table of derivatives
shiftDphi=: (-@[) shybyx makeDphi      NB. shift each row by #row to the right
makeGi=: (i. @ (0 { ])) shiftDphi phi  NB. make coeffs table then shift

NB. ---------------------------------------------------------
NB. makeTi
NB. Calculate T[i] for M[i]
NB.
NB. Syntax:
NB.   Ti=. makeTi mi , lambdai , Ng
NB. where:
NB.   mi      - lambdai multiplicity
NB.   lambdai - i-th eigenvalue of G
NB.   Ng      = #G , matrix G minimal polynom's order
NB.   Ti      - mi-by-Ng matrix, T[i]
NB.
NB. Test:
NB.    makeTi 4 1j1 7
NB. 1 1j1 0j2 _2j2   _4 _4j_4  0j_8
NB. 0   1 1j1  0j2 _2j2    _4 _4j_4
NB. 0   0   1  1j1  0j2  _2j2    _4
NB. 0   0   0    1  1j1   0j2  _2j2

lipows=: (1 { ]) (^ i.) (2 { ])   NB. lambdai ^ i. Ng
IOnmi=: - @ i. @ (0 { ])          NB. 0 _1 _2 ... -(mi-1)
makeTi=: IOnmi shybyx lipows      NB. mi-by-Ng table from vector's shifts

NB. ---------------------------------------------------------
NB. prepV
NB. Prepare eigenvalues of matrix y: remove dups, then count.
NB. Outputs 3 columns:
NB. - mi, lambdai multiplicity
NB. - lambdai, i-th eigenvalue of G
NB. - Ng = #y , matrix y minimal polynom's order
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

prepV=: ((1: #. =) ,. ~. ,. #) @ (1 {:: dgeev_jlapack_)

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

makeG=: (0 {:: ]) ((+ & getCols) {. ,.) (1 {:: ])

NB. ---------------------------------------------------------
NB. makeLtM
NB. Calculate matrix M or vector L(t) for equation M*A(t)=L(t)
NB.
NB. Syntax:
NB.   Lt=. ts makeLtM V
NB.   M=. 0 makeLtM V
NB. where:
NB.   V  - (#vm)-by-3 matrix, prepared eigenvalues of G, output of prepV
NB.   ts - sampling period, ts>0
NB.   M  - Ng-by-Ng matrix for equation M*A(t)=L(t)
NB.   Lt - Ng-vector, RHS L(t) for equation M*A(t)=L(t)
NB.   Ng = #G , matrix G minimal polynom's order
NB.
NB. Tests:
NB.    1.0 makeLtM prepV diagmat 4 4 3 2j2 2j_2 1j1 1j1
NB. 54.5982 54.5982 20.0855 _3.07493j6.71885 _3.07493j_6.71885 1.46869j2.28736 1.46869j2.28736
NB.
NB.    0 makeLtM prepV diagmat 4 4 3 2j2 2j_2 1j1 1j1
NB. 1    4   16      64  256      1024    4096
NB. 0    1    8      48  256      1280    6144
NB. 1    3    9      27   81       243     729
NB. 1  2j2  0j8  _16j16  _64 _128j_128  0j_512
NB. 1 2j_2 0j_8 _16j_16  _64  _128j128   0j512
NB. 1  1j1  0j2    _2j2   _4     _4j_4    0j_8
NB. 0    1  2j2     0j6 _8j8       _20 _24j_24

tspows=: [ (^ i.) (0 { ])               NB. ts ^ i. mi
prepMi=: makeGi * makeTi                NB. prepare Mi=. Gi*Ti
prepLti=: (^ @ ([ * (1 { ]))) * tspows  NB. prepare L[i](t)
prepLtMi=: prepLti`prepMi @. (0 = [)    NB. choose L[i](t) or M[i] to prepare
makeLtM=: ; @: (< @ prepLtMi " 1)       NB. make and merge all L[i](t) or M[i]

NB. ---------------------------------------------------------
NB. makeAt
NB. Solve equation M*A(t)=L(t) for A(t)
NB.
NB. Syntax:
NB.   At=. NxPMV makeAt ts
NB. where:
NB.   ts    - sampling period, ts>0
NB.   NxPMV - output of prexpm, being (Nx;P;M;V)
NB.   At    - Ng-vector, solution A(t) of equation M*A(t)=L(t)
NB.   Ng    = #G , matrix G minimal polynom's order

M=: (7 7$1 1.69783466514846393j4.84003649065663843 _20.5433106807080286j16.4351634688407167 _114.4259359297168j_71.5261830687437055 151.913015476576021j_675.265338457235771 3526.23206277198688j_411.22436149359072 7977.29994900821202j16368.9008822425258 1 1.69783466514846393j_4.84003649065663843 _20.5433106807080321j_16.4351634688407131 _114.425935929716786j71.5261830687437339 151.913015476576192j675.265338457235771 3526.23206277198688j411.224361493589868 7977.29994900820748j_16368.9008822425294 1 2.66145784252566076 7.08335784754134412 18.8520582947545918 50.1739583963255456 133.535875064456832 355.400101948825466 1 _0.0571271728225874048 0.00326351387470176892 _0.000186435321128999897 1.06505228103709516e_5 _6.08434257238971751e_7 3.4758128964473318e_8 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 2 0 0 0 0) " _
Lt=: (_2.44158297029179217j0.983055384353716422 _2.44158297029179261j_0.983055384353715755 4.55867025493447731 0.967961961441793028 1 0.569999999999999951 0.324899999999999967) " _
NB. At=: (1.00000000000000067j_4.44089209850062616e_16 0.570000000000000173j_3.72414106402510651e_28 0.162450000000000011j_4.7284015923396702e_30 0.0308871141985214044j_2.48294851946182879e_14 0.00476727342798542817j1.25320521616969433e_14 0.0003378506917448109j_2.14772796224925727e_15 6.85254686141710935e_5j3.54507642750158538e_16) " _
At=: (1 0.57 0.16245 0.0308871141985214044 0.00476727342798542817 0.0003378506917448109 6.85254686141710935e_5) " _

NB. makeAt=: , @: (2 & dumper) @ dgesv_jlapack_ @ (1 & dumper) @ (M ; (,. @: Lt))
makeAt=: dumper @: At

NB. =========================================================
NB. prexpm
NB. Prepare time-invariant parts for expm
NB.
NB. Syntax:
NB.   'Nx P M V'=. prexpm A;B
NB. where:
NB.   A - Nx-by-Nx state matrix of LTI system
NB.   B - Nx-by-Nu control input matrix of LTI system
NB.   V - (#vm)-by-3 matrix, prepared eigenvalues of G, output of prepV
NB.   M - Ng-by-Ng matrix for equation M*A(t)=L(t), output of makeLtM
NB.   P - Ng-by-Ng-by-Ng matrix, powers 0..(Ng-1) of G, output of makeP
NB.   G - Ng-by-Ng matrix, augmented LTI system, output of makeG
NB.   Ng = Nx + Nu
NB.   Nx >= 0
NB.   Nu >= 0

prexpm=: getCols @ (0 {:: ]) ; (makeP ; (0 & makeLtM ; ]) @ prepV) @ makeG
NB. Nx=: 4
NB. P=: 7 7 7$1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 2 1 _3 _1 1 _1 5 _1 0 2 1 5 3 _2 2 _4 1 _2 0 1 _4 5 4 4 3 _1 2 _2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _8 10 _11 2 8 _4 22 7 _5 9 0 _2 5 _15 0 _10 _21 _14 _16 _17 18 29 1 9 0 22 17 _5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _38 44 41 46 40 31 _20 37 _29 _22 _30 _18 _13 9 _102 28 _97 _10 _36 _79 132 75 _7 _76 _46 34 _17 107 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 192 _18 427 138 136 303 _534 _91 5 _311 _112 _78 _206 391 _476 246 225 294 48 69 _158 _225 195 _499 _68 86 _264 785 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1946 _964 367 _650 _36 457 _988 _1369 705 _476 382 46 _429 1003 722 _200 3321 1154 460 2027 _4360 _1983 1499 294 1214 818 175 617 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2340 _2122 _9999 _5594 _2224 _5771 11490 _2485 2063 6569 4172 1774 3772 _7115 14056 _7946 5371 _4102 _1432 4307 _11582 1193 1697 14097 6536 4298 9202 _16517 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
NB. M=: 7 7$1 1.69783466514846393j4.84003649065663843 _20.5433106807080286j16.4351634688407167 _114.4259359297168j_71.5261830687437055 151.913015476576021j_675.265338457235771 3526.23206277198688j_411.22436149359072 7977.29994900821202j16368.9008822425258 1 1.69783466514846393j_4.84003649065663843 _20.5433106807080321j_16.4351634688407131 _114.425935929716786j71.5261830687437339 151.913015476576192j675.265338457235771 3526.23206277198688j411.224361493589868 7977.29994900820748j_16368.9008822425294 1 2.66145784252566076 7.08335784754134412 18.8520582947545918 50.1739583963255456 133.535875064456832 355.400101948825466 1 _0.0571271728225874048 0.00326351387470176892 _0.000186435321128999897 1.06505228103709516e_5 _6.08434257238971751e_7 3.4758128964473318e_8 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 2 0 0 0 0
NB. V=: 5 3$1 1.69783466514846326j4.84003649065663755 7 1 1.69783466514846326j_4.84003649065663755 7 1 2.66145784252566076 7 1 _0.0571271728225874187 7 3 0 7
NB. prexpm=: (Nx ; P ; M ; V) " _

NB. ---------------------------------------------------------
NB. expm
NB. Calculate matrix exponent
NB.   Phi(ts) = exp(A*ts)
NB. and Cauchy integral
NB.   Gamma(ts) = Integral(exp(A*(ts-t))*dt,t=0..ts)*B
NB. for sampling period ts via Lagrange-Sylvester interpolation
NB. polynome [1, p. 88]
NB.
NB. Syntax:
NB.   'Phi Gamma'=. NxPMV expm ts
NB. where:
NB.   NxPMV - output of prexpm, being (Nx;P;M;V)
NB.   ts    - sampling period, ts>0
NB.   Phi   - Nx-by-Nx matrix, matrix exponent
NB.   Gamma - Nx-by-Nu matrix, Cauchy intergal
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

NB. Syntax: texpm A;B;ts

texpm=: 3 : 0
'A B ts'=. y
P =. 0 {:: (prexpm A;B) expm ts
eigA=. /:~ ^ ts * (1 {:: dgeev_jlapack_ A)
eigP=. /:~ (1 {:: dgeev_jlapack_ P)      NB. FIXME: this sort gives wrong result
err=. %: +/ *: | eigA - eigP
smoutput (eigA ,. eigP) ; (,. | eigA - eigP) ; err
0 = err
)

NB. Syntax: testexpm ''

testexpm=: 3 : 0
NB.'A0 B0'=.  4 splitbyx _4 + ?  4  7 $ 10
NB.'A1 B1'=.  6 splitbyx _4 + ?  6 11 $ 10
NB.'A2 B2'=.  8 splitbyx _4 + ?  8 13 $ 10
NB.'A3 B3'=. 10 splitbyx _4 + ? 10 17 $ 10

NB.'ts0 ts1 ts2 ts3'=. 0.01 * >: ? 4 $ 100

NB.texpm &> (< A0;B0;ts0) , (< A1;B1;ts1) , (< A2;B2;ts2) , (< A3;B3;ts3)

'A0 B0'=: 4 splitbyx 4 7 $ 2  1 _3 _1 1 _1  5 _1  0  2  1 5  3 _2 2 _4  1 _2 0  1 _4 5  4  4  3 _1  2 _2
ts0=: 0.57
smoutput 2 3 $ 'A0' ; 'B0' ; 'ts0' ; A0 ; B0 ; ts0

texpm A0;B0;ts0

)
