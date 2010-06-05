NB. expm.ijs
NB. Calculate matrix exponent and Cauchy integral from state and
NB. control input matrices of LTI system
NB.
NB. prexpm  prepare time-invariant parts for expm
NB. expm    calculate matrix exponent and Cauchy integral
NB.
NB. References:
NB. - Podchukaev V.A. Theory of informational processes and systems. - M.,
NB.   2006. (Подчукаев В. А. Теория информационных процессов и систем. -
NB.   М.: Гардарики, 2006 - 209 с.)
NB.   URL: http://www.sgau.ru/uit/Book3.htm
NB. - Andrievskiy B.R., Fradkov A.L. Selected chapters of automatic
NB.   control theory with MATLAB examples. - SPb., 2000 (Андриевский Б.
NB.   Р., Фрадков А. Л. Избранные главы теории автоматического управления
NB.   с примерами на языке MATLAB. - СПб.: Наука, 2000. - 475 с., ил. 86)
NB.
NB. Resources:
NB. - http://www.jsoftware.com/jwiki/...
NB. - http://www.dvgu.ru/forum/...
NB.
NB. TODO:
NB. - consider s/@:/@/g when possible
NB. - consider complex A or B => non-self-adjoined eigenvalues
NB. - consider B is vector
NB.
NB. 2008-02-04 0.0.0 Igor Zhuravlov |.'ur.ugvd.ciu@rogi'

require '~user/projects/lapack/lapack.ijs'     NB. '~addons/math/lapack/lapack.ijs'
NB. need_jlapack_ 'geev gels'
require '~user/projects/lapack/geev.ijs'       NB. '~addons/math/lapack/geev.ijs'
require '~user/projects/lapack/gels.ijs'       NB. '~addons/math/lapack/gels.ijs'

coclass 'tau'

NB. ===========================================================================
NB. Utilities

getA=: 0 {:: ]        NB. extract A
getB=: 1 {:: ]        NB. extract B

getM=: 0 {:: [        NB. extract M
getV=: 1 {:: [        NB. extract V
getP=: 2 {:: [        NB. extract P
getNx=: 3 {:: [       NB. extract Nx

getli=: 0 { ]         NB. extract lambdai
getic=: 1 { ]         NB. extract ic
getmi=: 2 { ]         NB. extract mi
getIOs=: 3 { ]        NB. extract IOs
getNm=: 4 { ]         NB. extract Nm
getk=: 5 { ]          NB. extract k

getCols=: {: @ $      NB. get columns count of table y
reim=: 9 11 o."0 _ ]  NB. extract Re(y) and Im(y) from complex y
c2r=: ,/ @ reim       NB. realificate complex y
idmat=: =@i.          NB. identity matrix of size y
mp=: +/ .*            NB. matrix product of x and y
shybyx=: (|.!.0)"0 1  NB. shift y with step x

NB. ---------------------------------------------------------
NB. makeGi
NB. Calculate G[i] for M[i]
NB.
NB. Syntax:
NB.   Gi=. makeGi lambdai , ic , mi , IOs , Nm , k
NB. where:
NB.   k       - matrix G minimal polynom's order, k = Ng
NB.   Nm      = # M
NB.   IOs     - IO 1st row (atom) of corresp. M[i] (L[i](t)) in M (L(t))
NB.   mi      - multiplicity, taking self-adjoiners into account
NB.   ic      - datatype flag: 0 for real, 1 for complex
NB.   lambdai - i-th eigenvalue of G
NB.   Gi      - mi-by-k matrix, G[i]
NB.
NB. Test:
NB.    makeGi 1j1 1 4 7 15 9
NB. 1 1 1 1  1  1   1   1   1
NB. 0 1 2 3  4  5   6   7   8
NB. 0 0 2 6 12 20  30  42  56
NB. 0 0 0 6 24 60 120 210 336

xIO=: i. @ x:
xIOk=: xIO @ getk               NB. i. x: k
xIOmi=: xIO @ getmi             NB. i. x: mi
makeRji=: ! * ! @ [             NB. (!y)%(!y-x)
makeGi=: xIOmi makeRji"0/ xIOk  NB. call (mi makeRji k) for each pair (mi,k)

NB. ---------------------------------------------------------
NB. makeTi
NB. Calculate T1[i] or T2[i] and T3[i] for M[i]
NB.
NB. Syntax:
NB.   T1i=. makeTi lambdai , ic , mi , IOs , Nm , k
NB.   'T2i T3i'=. reim makeTi lambdai , ic , mi , IOs , Nm , k
NB. where:
NB.   k           - matrix G minimal polynom's order, k = Ng
NB.   Nm          = # M
NB.   IOs         - IO 1st row (atom) of corresp. M[i] (L[i](t)) in M (L(t))
NB.   mi          - multiplicity, taking self-adjoiners into account
NB.   ic          - datatype flag: 0 for real, 1 for complex
NB.   lambdai     - i-th eigenvalue of G
NB.   T1i,T2i,T3i - mi-by-k matrix, T1[i], T2[i] or T3[i] respectively
NB.
NB. Tests:
NB.    makeTi 4 0 2 0 15 9
NB. 1 4 16 64 256 1024 4096 16384 65536
NB. 0 1  4 16  64  256 1024  4096 16384
NB.
NB.    x: reim makeTi 1j1 1 4 7 15 9
NB. 1 1 0 _2 _4 _4  0  8 16
NB. 0 1 1  0 _2 _4 _4  0  8
NB. 0 0 1  1  0 _2 _4 _4  0
NB. 0 0 0  1  1  0 _2 _4 _4
NB.
NB. 0 1 2  2  0 _4 _8 _8  0
NB. 0 0 1  2  2  0 _4 _8 _8
NB. 0 0 0  1  2  2  0 _4 _8
NB. 0 0 0  0  1  2  2  0 _4

lipows=: getli (^ i.) getk    NB. lambdai ^ i. k
IOnmi=: - @ i. @ getmi        NB. 0 _1 _2 ... -(mi-1)
makeTi=: IOnmi shybyx lipows  NB. form mi-by-k table from vector's shifts

NB. ---------------------------------------------------------
NB. prepV
NB. Prepare eigenvalues from y: remove dups and self-adjoiners,
NB. classify and count. Outputs 5 columns:
NB. - lambdai, i-th eigenvalue
NB. - ic, datatype flag: 0 for real, 1 for complex
NB. - mi, multiplicity, taking self-adjoiners into account
NB. - IOs, IO 1st row (atom) of corresp. M[i] (L[i](t)) in M (L(t))
NB. - Nm = # M
NB. - k, matrix G minimal polynom's order, k = Ng
NB.
NB. If:
NB.   'vlambda vic vm vIOs vNm vk' =. |: prepV eigenvalues_of_G
NB. then
NB.   vk -: (# vm) $ k
NB.   vNm -: (# vm) $ Nm
NB.   (+/ vm) = k
NB.   (# vm) = (# vlambda) = (r - C)
NB.   ic -: ((0 ~: im) vlambda)
NB. where
NB.   r - quantity of unique eigenvalues
NB.   C - quantity of unique complex eigenvalues without
NB.       self-adjoiners
NB.
NB. Test:
NB.    prepV 4 4 3 2j2 2j_2 1j1 1j_1 1j1 1j_1   NB. r=6, C=2
NB.   4 0 2 0 15 9
NB.   3 0 1 2 15 9
NB. 2j2 1 2 3 15 9
NB. 1j1 1 4 7 15 9

cnt=: 1: #. =                   NB. count atoms
IOs=: _1 & shybyx @ (+/\) @: *  NB. IO 1st row (atom) of corresponding M[i] (L[i](t)) in M (L(t))
ic=: (0 < im) @ ~.              NB. datatype flag
im=: 11 o. ]                    NB. take imaginary part
nnegim=: ] #~ 0 <: im           NB. filter out atoms with negative imaginary part
prepV=: (~. ,. (cnt (] ,. ([ (* ([ ,. (IOs ,. (+/ @: *))) ]) (>: @ ]))) ic)) @ nnegim ,. #

NB. ---------------------------------------------------------
NB. makeG
NB. Build G matrix of augmented LTI system: G = ( A  B ) Nx
NB.                                             ( 0  0 ) Nu
NB.                                               Nx Nu
NB. Syntax:
NB.   G=. makeG A;B

makeG=: getA ((+ & getCols) {. ,.) getB

NB. ---------------------------------------------------------
NB. makeP
NB. Make report of G powers
NB.
NB. Syntax:
NB.   P=. makeP G
NB. where
NB.   G - Ng-by-Ng matrix, augmented LTI system, output of makeG
NB.   P - Ng-by-Ng-by-Ng matrix, powers 0..(Ng-1) of G
NB.
NB. Test:
NB.    ;/ makeP_tau_ 3 3 $ 5 5 0 9 0 5 1 6 7
NB. ┌─────┬─────┬────────┐
NB. │1 0 0│5 5 0│70 25 25│
NB. │0 1 0│9 0 5│50 75 35│
NB. │0 0 1│1 6 7│66 47 79│
NB. └─────┴─────┴────────┘
NB.
NB. Notes:
NB. - powers are calculated via repeated squaring to
NB.   reduce operations from O(Ng*Ng) to O(Ng*log Ng), see
NB.   http://www.jsoftware.com/jwiki/Essays/Linear_Recurrences
NB. - 0-th power (identity matrix) is substituted directly
NB.   without calculation

p2b=: < @ I. @ |.               NB. cvt bits of y to powers, then box it
pows=: p2b"1 @ #: @ i.          NB. call p2b for each power y represented binary
topow=: mp/ @ (mp~ @ ] ^: [)    NB. produce table y to powers from list x, then product all
topows=: > @ (topow &. >)       NB. apply dyad topow under boxes, then open
prepP=: (pows @ [) topows ]     NB. form boxed list of y ^ i. x
repl0=: (idmat @ [)`0:`prepP }  NB. replace 0-th element in y by identity matrix of size x
makeP=: # repl0 <               NB. call: (#G) repl0 (<G)

NB. ---------------------------------------------------------
NB. makeLtM
NB. Calculate matrix M or vector L(t) for equation M*A(t)=L(t)
NB.
NB. Syntax:
NB.   Lt=. T makeLtM V
NB.   M=. 0 makeLtM V
NB. where:
NB.   V  - (#vm)-by-5 matrix, prepared eigenvalues of G, output of prepV
NB.   T  - sampling period, T>0
NB.   M  - Nm-by-Ng matrix for equation M*A(t)=L(t)
NB.   Lt - Nm-vector, RHS L(t) for equation M*A(t)=L(t)
NB.   Nm = +/ vm * (ic+1), see prepV
NB.   Ng = #G
NB.
NB. Test:
NB.    1.1 makeLtM prepV 4 4 3 2j2 2j_2 1j1 1j_1 1j1 1j_1
NB. 81.4509 89.596 27.1126 _5.31123 _5.84235 7.29669 8.02636 1.36268 1.49895 1.64884 1.81372 2.67733 2.94507 3.23958 3.56353
NB.    0 makeLtM prepV 4 4 3 2j2 2j_2 1j1 1j_1 1j1 1j_1
NB. 1 4 16  64 256 1024 4096 16384  65536
NB. 0 1  8  48 256 1280 6144 28672 131072
NB. 1 3  9  27  81  243  729  2187   6561
NB. 1 2  0 _16 _64 _128    0  1024   4096
NB. 0 1  4   0 _64 _320 _768     0   8192
NB. 0 2  8  16   0 _128 _512 _1024      0
NB. 0 0  4  24  64    0 _768 _3584  _8192
NB. 1 1  0  _2  _4   _4    0     8     16
NB. 0 1  2   0  _8  _20  _24     0     64
NB. 0 0  2   6   0  _40 _120  _168      0
NB. 0 0  0   6  24    0 _240  _840  _1344
NB. 0 1  2   2   0   _4   _8    _8      0
NB. 0 0  2   6   8    0  _24   _56    _64
NB. 0 0  0   6  24   40    0  _168   _448
NB. 0 0  0   0  24  120  240     0  _1344

Tpows=: [ (^ i.) getmi                NB. T ^ i. mi
prepMi=: makeGi * makeTi              NB. prepare Mi=. Gi*Ti
prepLti=: (^ @ ([ * getli)) * Tpows   NB. prepare L[i](t)
prepLtMi=: prepLti`prepMi @. (0 = [)  NB. choose L[i](t) or M[i] to prepare
makeLtMi=: getic (c2r ^: [) prepLtMi  NB. complete L[i](t) or M[i], realificate if required
makeLtM=: ; @: (< @: makeLtMi " 1)    NB. make and merge all L[i](t) or M[i]

NB. ---------------------------------------------------------
NB. makeAt
NB. Solve equation M*A(t)=L(t) for A(t)
NB.
NB. Syntax:
NB.   At=. MVPN makeAt T
NB. where:
NB.   T    - sampling period, T>0
NB.   MVPN - output of prexpm, being (M;V;P;Nx)
NB.   At   - Nm-vector, solution A(t) of equation M*A(t)=L(t)
NB.   Nm = +/ vm * (ic+1), see prepV

makeAt=: gels_jlapack_ @: (getM ; ] makeLtM getV)

NB. ===========================================================================
NB. prexpm
NB. Prepare time-invariant parts for expm
NB.
NB. Syntax:
NB.   'M V P Nx'=. prexpm A;B
NB. where:
NB.   A - Nx-by-Nx state matrix of LTI system
NB.   B - Nx-by-Nu control input matrix of LTI system
NB.   M - Nm-by-Ng matrix for equation M*A(t)=L(t), output of makeLtM
NB.   V - (#vm)-by-5 matrix, prepared eigenvalues of G, output of prepV
NB.   P - Ng-by-Ng-by-Ng matrix, powers 0..(Ng-1) of G, output of makeP
NB.   G - Ng-by-Ng matrix, augmented LTI system, output of makeG
NB.   Nm = +/ vm * (ic+1), see prepV
NB.   Ng = Nx + Nu
NB.   Nx >= 0
NB.   Nu >= 0

prexpm=: (((0 & makeLtM ; ]) @ prepV @ (2b010 & geev_jlapack_)) , (< @ makeP)) @ makeG , (< @ getCols @ getA)

NB. ---------------------------------------------------------
NB. expm
NB. Calculate matrix exponent
NB.   exp(A*T)
NB. and Cauchy integral
NB.   Integral(exp(A*(T-t))*dt,t=0..T)*B
NB. for sampling period T via Lagrange-Sylvester interpolation
NB. polynome
NB.
NB. Syntax:
NB.   'EA IE'=. MVPN expm T
NB. where:
NB.   MVPN - output of prexpm, being (M;V;P;Nx)
NB.   T    - sampling period, T>0
NB.   EA   - Nx-by-Nx matrix, matrix exponent
NB.   IE   - Nx-by-Nu matrix, Cauchy intergal
NB.   Nu   = (#P)-Nx

splitbyx=: {."1 ; }."1         NB. split table y at column x
extract=: [ splitbyx {.        NB. extract 1st x rows from table y
augsys=: +/ @ (getP * makeAt)  NB. find exponent of augmented system
expm=: getNx extract augsys

NB. ===========================================================================
NB. Test suite

NB. pass_result=: texpm A;B;T;EA;IE

texpm=: 3 : 0
match=. matchclean_jlapack_;;
smoutput 'EA IE'=. ((prexpm @ (0 1 & {)) expm (2 & {::)) y
smoutput a=. EA match 3 {:: y
smoutput b=. IE match 4 {:: y
(0 pick a) *. 0 pick b
)

NB. testexpm ''

testexpm=: 3 : 0
A=. 4 4 $ 0.8462 0.8381 0.8318 0.3046 0.5252 0.0196 0.5028 0.1897 0.2026 0.6813 0.7095 0.1934 0.6721 0.3795 0.4289 0.6822
B=. 4 3 $ 0.3028 0.3784 0.4966 0.5417 0.8600 0.8998 0.1509 0.8537 0.8216 0.6979 0.5936 0.6449
T0=. 1.0
EA0=. 4 4 $ 3.4818 2.3874 2.8282 1.1815 1.2996 1.8528 1.4749 0.6245 1.0820 1.5227 2.8538 0.6937 2.0741 1.6433 1.9914 2.5156
IE0=. 4 3 $ 1.4393 2.4832 2.7270 1.0422 1.8137 1.9133 0.8258 2.1815 2.2029 1.6991 2.2054 2.3720
T1=. 0.35
EA1=. 4 4 $ 1.415 0.4041 0.4353 0.1687 0.2383 1.071 0.2455 0.09722 0.1339 0.2975 1.333 0.1037 0.3373 0.2218 0.2581 1.302
IE1=. 4 3 $ 0.1878 0.2837 0.3344 0.2219 0.3642 0.3825 0.1037 0.4019 0.394 0.3177 0.3173 0.3441

Aocw=: 4 4 $ _0.0069 0.0139 0 _9.81 _0.0905 _0.3149 235.8928 0 0.0004 _0.0034 _0.4282 0 0 0 1 0
Bocw=: 4 1 $ 0
Tocw=: 0.1486
EAocw=: 4 4 $ 0.9990 0.0020 _0.0709 _1.4566 _0.0121 0.9460 33.0649 0.0091 0.0001 _0.0005 0.9301 0 0 0 0.1435 1
IEocw=: 4 1 $ 0

NB. texpm &> (< A;B;T0;EA0;IE0) , (< A;B;T1;EA1;IE1) , (< Aocw;Bocw;Tocw;EAocw;IEocw)
texpm &> (< Aocw;Bocw;Tocw;EAocw;IEocw)
)
