NB. util.ijs
NB. Linear time-invariant (LTI) system's utilities
NB.
NB. expm       matrix exponent and Cauchy integral
NB. logm       matrix logarithm
NB. powm       matrix y to power x
NB. powsm      X-Y-Z matrix (report) of matrix y powers 0..(x-1)
NB.
NB. shiftdiag  add element[s from] x to diagonal of matrix y
NB. logspace   create vector of logarithmically spaced numbers
NB. plot       plot data as multiplot or subplots grid
NB.
NB. rndmat     generate random matrix
NB. rndmatpe   generate random matrix with positive eigenvalues
NB. rndmatne   generate random matrix with negative eigenvalues
NB.
NB. Resources:
NB. - http://www.jsoftware.com/jwiki/...
NB. - http://www.dvgu.ru/forum/...
NB.
NB. TODO:
NB. - rework pow[s]m to base on eigen-decomposition:
NB.   if   A   = L * D * R, where D = diagmat(d)
NB.   then A^n = L * diagmat(d^n) * R
NB.
NB. Version: 1.0.0 2008-03-30
NB. Copyright: Igor Zhuravlov igor@uic.dvgu.ru
NB. License: Version 3 of the GNU GPL or any later version

script_z_ '~system/packages/math/mathutil.ijs'  NB. mp
script_z_ '~system/packages/math/makemat.ijs'   NB. idmat diagmat
require '~user/projects/lapack/lapack.ijs'      NB. '~addons/math/lapack/lapack.ijs'
require '~user/projects/lapack/gesvd.ijs'       NB. need_jlapack_ 'gesvd'

require 'plot'

coclass 'tau'

NB. =========================================================
NB. Local utilities

NB. =========================================================
NB. Interface verbs

NB. ---------------------------------------------------------
NB. expm                                                  1 0
NB. Calculate matrix exponent
NB.   Φ(ts) = exp(A*ts)
NB. and Cauchy integral
NB.   Γ(ts) = Integral(exp(A*(ts-t))*dt,t=0..ts)*B
NB. for sampe time ts.
NB.
NB. Syntax:
NB.   'Phi Gamma'=. VPIVB expm ts
NB. where:
NB.   VPIVB - boxed quad (V;Λ;inv(V);B), where
NB.           V      - Nx-by-Nx matrix, where columns are
NB.                    right eigenvector of A
NB.           Λ      - Nx-vector of eigenvalues {λ[i]} of A,
NB.                    also known as poles of system
NB.           inv(V) - Nx-by-Nx matrix, inversion of V
NB.           B      - Nx-by-Nu control input matrix
NB.           A      - Nx-by-Nx state matrix
NB.   ts    > 0, sample time
NB.   Phi   - Nx-by-Nx table, matrix exponent
NB.   Gamma - Nx-by-Nu table, Cauchy intergal
NB.
NB. Applications:
NB.   'Phi Gamma'=. (A;B) ((prexpm @ [) expm ]) ts

expm=: (4 : 0) " 1 0
  'V Lambda invV B'=. x

  NB. A*V=V*diagmat(Λ) <=> A=V*diagmat(Λ)*inv(V) =>
  NB. exp(A*ts) = V*exp(diagmat(Λ*ts))*inv(V) = V*diagmat(exp(Λ*ts))*inv(V)
  Phi=. V mp ((^ Lambda * y) * invV)

  NB. Integral(exp(λ*(ts-t))*dt,t=0..ts) = γ = (exp(λ*ts)-1)/λ, if λ≠0
  NB.                                        = ts,              if λ=0, =>
  NB. Integral(exp(A*(ts-t))*dt,t=0..ts) = V*diagmat({γ[i]})*inv(V)
  Gamma=. V mp ((y ([ ((I. @: (0 = ])) }) (((- & 1) @: (^ @: *)) % ])) Lambda) * invV) mp B

  Phi;Gamma
)

NB. ---------------------------------------------------------
NB. powm
NB. Raise matrix y to power x
NB.
NB. Syntax:
NB.   P=. [x] powm y
NB. where
NB.   y - N-by-N table
NB.   x - integer >= 0, power, default is #y
NB.   P - N-by-N table, matrix y in power x
NB.   N >= 0


NB. ---------------------------------------------------------
NB. powsm
NB. Make report of table y powers: I y y^2 ... y^(x-1)
NB.
NB. Syntax:
NB.   P=. [x] powsm y
NB. where
NB.   y - N-by-N table
NB.   x - integer >= 0, powers count, default is #y
NB.   P - x-by-N-by-N report, matrix y in powers 0..(x-1)
NB.   N >= 0
NB.
NB. Test:
NB.    ;/ powsm 3 3 $ 5 5 0 9 0 5 1 6 7
NB. ┌─────┬─────┬────────┐
NB. │1 0 0│5 5 0│70 25 25│
NB. │0 1 0│9 0 5│50 75 35│
NB. │0 0 1│1 6 7│66 47 79│
NB. └─────┴─────┴────────┘
NB.
NB. Notes:
NB. - powers are calculated via repeated squaring, see
NB.   http://www.jsoftware.com/jwiki/Essays/Linear_Recurrences
NB. - 0-th power (identity matrix) is substituted directly
NB.   without calculation

p2b=: < @ I. @ |.                 NB. cvt bits of y to powers, then box it
pows=: p2b"1 @ #: @ i. @ [        NB. call p2b for each power x represented binary
topow=: mp/ @ (mp~ @ ] ^: [)      NB. produce table y to powers from list x, then product all
topows=: topow &. >               NB. apply dyad topow under boxes
make1=: idmat @ # &. > @ ]        NB. make boxed identity matrix of size N
repl1=: (make1) 0: } topows       NB. replace by identity matrix in 1st box
prepP=: pows > @ repl1 < @ ]      NB. form report of y ^ i. x
make0=: (0 $~ 0 , $ @ ]) " _      NB. make zero report: 0 N N $ 0
check0=: make0`prepP @. (0 ~: [)  NB. choose report type depending on x=0
powsm=: (# $: ]) :check0          NB. force dyadic call: (#@]) check0 ]

NB. ---------------------------------------------------------
NB. shiftdiag                                             1 2
NB. For scalar/vector x and matrix y make x*I+y
NB.
NB. Syntax:
NB.   s=. x shiftdiag y
NB. where
NB.   y - N-by-N table
NB.   x - numeric scalar of N-vector, shift for y's diagonal
NB.   s - N-by-N table, y+x*idmat(#y)
NB.   N >= 0

xplusdiagy=: + (< 0 1) & |:                     NB. new diagonal: x + diag(y)
linIOSdiagy=: (>: * i.) @ # @ ]                 NB. linear IOS of y's diagonal
shiftdiag=: (xplusdiagy linIOSdiagy } ]) " 1 2  NB. replace diagonal

NB. ---------------------------------------------------------
NB. logspace
NB. Create vector of logarithmically spaced numbers
NB.
NB. Syntax:
NB.   grid=. logspace (a , b , n)
NB. where:
NB.   a    - set start value to 10^a
NB.   b    - set final value to 10^b
NB.   n   >= 0, integer, atoms in grid
NB.   grid - n-vector of logarithmically spaced numbers in
NB.          range [10^a, 10^b]

logspace=: (10 & ^) @ steps @ (- & 0 0 1)

NB. ---------------------------------------------------------
NB. plot                                                1 1 1
NB. Plot data as multiplot or as subplots. When multiplot,
NB. plots grid may have row and column titles, and some
NB. axes linked across. When subplots, each plot may have its
NB. own options.
NB.
NB. Synax:
NB.   p=. [p] plot ps;pf;so;rt;ct;kt;xx;<dat
NB.   p=. [p] plot ps;pf;so;kt;xx;<dat
NB. where
NB.   ps  - string, plot start commands and options. For
NB.         multiplot: 'multi' command, axes linking, plot
NB.         title, plot-wide key options etc. For subplots:
NB.         plot title, 'sub' command, plot-wide key options
NB.         etc.
NB.   pf  - string, plot finish commands: 'show', 'save' etc.
NB.   so  - subplot options, any one of:
NB.         - string, options for all [sub]subplot
NB.         - boxed [R-by-]C array, each element #i contains
NB.           under box any one of:
NB.           si1          - string, single subplot options
NB.           si1;...;siSi - boxed strings, options for each
NB.                          subsubplot
NB.         Not any option takes effect for multiplot, see
NB.         help.
NB.         Only last subsubplot's key option takes effect.
NB.   rt  - row titles for multiplot, any one of:
NB.         s1            - string, prefix for row titles,
NB.                         empty means no row titles
NB.         s1;s2;...;sR  - boxed strings, row titles
NB.   ct  - column titles for multiplot, any one of:
NB.         s1            - string, prefix for column titles,
NB.                         empty means no column titles
NB.         s1;s2;...;sC  - boxed strings, column titles
NB.   kt  - subplot key titles, any one of:
NB.         - string, prefix for all [sub]subplot's key
NB.           titles, empty means no key titles
NB.         - boxed [R-by-]C array, each element #i contains
NB.           under box any one of:
NB.           kti1           - key titles for subplot
NB.           kti1;...;ktiSi - key titles for subsubplots
NB.           where each ktij is any one of:
NB.           s1            - string, prefix for #i subplot's
NB.                           or all #ij subsubplots' key
NB.                           titles, empty means no titles
NB.           s1;s2;...;sKi - boxed strings, key titles
NB.         Only last subsubplot's key titles take effect.
NB.   xx  - abscissa grid[s], any one of:
NB.         - [K-by-]N array of non-complex numeric atoms,
NB.           plot-wide x axis (axes), N=0 means default axis
NB.         - boxed [R-by-]C array, x axes, each element #i
NB.           contains under box any one of:
NB.           - xxi1           - [Ki-by-]Ni array of non-
NB.                              complex numeric atoms, x axis
NB.                              (axes) for subplot #i, or all
NB.                              subsubplots #ij, Ni=0 means
NB.                              default axis
NB.           - xxi1;...;xxiSi - boxed x axes for subsubplots
NB.   dat - boxed [R-by-]C array, data to plot; each element
NB.         #i contains under box any one of:
NB.         - dati1                  - data for subplot
NB.         - dati1;dati2;...;datiSi - data for subsubplots
NB.         where each datij is [Kij-by-]Nij numeric array.
NB.         If Si>1, then subplot consists of subsubplots.
NB.         If some of subsubplots have keys then the only
NB.         last keyed subsubplot's keys are displayed.
NB.         If Kij>1, then correspondent subsubplot #ij has
NB.         Kij graphics, and may have Kij keys (non-empty
NB.         corresp. entry in kt).
NB.   p   - plot object
NB.   R  >= 0, integer, rows in multiplot grid
NB.   C  >= 0, integer, columns in multiplot grid
NB.   Si >= 0, integer, subsubplots in subplot #i
NB.   Nij>= 0, integer, abscissa grid length for subplot #ij
NB.   Kij>= 0, integer, keys in subplot #ij
NB.
NB. Applications:
NB.   rc=. $ {. LTI1_Yimp
NB.   ps=. 'reset;multi ' , (": rc) , ';xgroup 0;ygroup 0;title Impulse response'
NB.   pf=. 'show'
NB.   so=. 'type line,marker;keypos ctie;keystyle lbhc;penstyle 0 1 2'
NB.   plot ps;pf;so;'Out #';'In #';'LTI #';tgrid; (< " 2) 0 3 |: 9 o. LTI1_Yimp , LTI2_Yimp ,: LTI2_Yimp1
NB.
NB. TODO:
NB. - check Nx|Ny|Nu == 0

plot=: (3 : 0) " 0 1
  (conew 'jzplot') plot y                    NB. supply new plot object
:
  ib=. * @ L.                                NB. is a box?
  ip=. -. @ (ib +. ('' & -:)) @ ]            NB. is a prefix?
  fsh=. (_2 & ({.!.1)) @ $                   NB. fill shape: (c)->(1 c) or (r c)->(r c)
  rsh=. {. @ fsh                             NB. rows in shape (1 for vector)
  p2t=. ((<"1 @: ((, ":) " 1 0)) i.)~ ^: ip  NB. convert prefix y to boxed key titles vector of length x
  p2k=. ; @ ((<'"') ,.~ (<' "') ,. p2t)      NB. wrap titles by double quotes, then raze to make keys list string

  if. 8 -: # y do.
    NB. mplot special case: additional po preprocessing

    'ps pf so rt ct kt xx dat'=. y
    'r c'=. fsh dat

    NB. expand row and column prefixes to titles list strings; wrapping each title by double quotes
    ps=. ps , ';ycaption' , (r p2k rt) , ';xcaption' , (c p2k ct)

  else.
    'ps pf so kt xx dat'=. y
  end.

  NB. expand key prefixes to key titles list strings; wrapping each key by double quotes
  kt=. dat ((p2k~ rsh)~ ` ($: L: _1) @. (ib @ [)) kt
gkt=: kt
smoutput 'kt' ; <kt

  NB. pair-wise conditional append non-empty key titles (kt), prepended by string ';key', to subplot options (so)
  so=. so ((, ((';key' & ,) ^: (* @ #))) L: 0) kt
gso=: so
smoutput 'so' ; <so

  NB. for each xx and dat pairs (xxi,dati): dati=. if empty(xxi) then (dati) else (xxi;dati)
  xdat=. xx ((; ^: (* @ # @ [)) ` ($: L: _1) @. (ib @ ])) dat
gxdat=: xdat
smoutput 'xdat' ; <xdat

  NB. couple each leaf from so with corresp. rank-1 box from xdat
NB. +++  sdat=. so ((, & <) L: 0 _1) xdat
  sdat=. so ((, & <) L: 0 _2) xdat
gsdat=: sdat
smoutput 'sdat' ; <sdat

  NB. under each sdat box: if it is subsubplots data then convert
  NB. it from boxed Si-vector to boxed Si-by-2 table,
  NB. finally ravel whole [R-by-]C array to form multiplot pd data
  pdat=. , ((> ^: (ib @: (0 & {::))) &. >) sdat
gpdat=: pdat
smoutput 'plot pdat=' ; <pdat

  pd__x ps

  if. 8 -: # y do.
    pd__x pdat
  else.
    NB. subplots special case: send each [sub]subplot's data
    NB. to pd consequently to enshure proper pd initialization

    NB. 1) prepare pdat to form list (pdat1;pdat2;...;pdatRC)
    NB.    where: RC -: R*C
    NB.           pdati is either (soi;xdati) or (soi1;xdati1;soi2;xdati2;...;soiSi;xdatiSi)
    NB. 2) until prepared pdat is empty:
    NB.    2.1) apply pd to unboxed head from list
    NB.    2.2) call itself with tail supplied
    ((($: @ }. @ [) (pd__x @ (0 & {::))) ^: (0 < #)) (, ; pdat)

  end.

  pd__x pf
  x
)

NB. ---------------------------------------------------------
NB. rndmat
NB. Generate rectangular random matrix of values in range [0,10)
NB.
NB. Syntax:
NB.   mat=. rndmat size
NB.   mat=. rndmat rows cols

rndmat=: 0.1 * (? @ (100 $~ ({. , {:)))

NB. ---------------------------------------------------------
NB. rndmatpe
NB. Generate square random matrix with positive eigenvalues
NB. in range (0,10]
NB.
NB. Syntax:
NB.   mat=. rndmatpe size

rndmatpe=: 3 : 0
d=. diagmat 0.1 + 0.1 * ? y $ 100
o=. 2b001 geev_jlapack_ rndmat y  NB. o=. 2b100 gesvd_jlapack_ rndmat y
m=. o mp d mp + |: o
)

NB. ---------------------------------------------------------
NB. rndmatne
NB. Generate square random matrix with negative eigenvalues
NB. in range [-10,0)
NB.
NB. Syntax:
NB.   mat=. rndmatne size

rndmatne=: 3 : 0
d=. diagmat _0.1 + _0.1 * ? y $ 100
o=. 2b100 gesvd_jlapack_ rndmat y
m=. o mp d mp + |: o
)

NB. =========================================================
NB. Test suite

NB. Syntax: testplot ''
NB. TODO: rewrite

testplot=: 3 : 0
  mba=. < @: ? @ ($ & 10)  NB. make random boxed array

  po=. 'multi 2 3;xgroup 0;ygroup 0;title Test mplot'
  so=. 2 3 $ (<'type line') , (<'type line;keypos ctie;keystyle lbhc;penstyle 0 1 2 3') , (<'type line';'type stick') , (<'type line';'type dot;penstyle 5;keypos ctie;keystyle lbhc;penstyle 0 1 2') , (<'type line;keypos ctie;keystyle lbhc';'type bar;keypos ctie;keystyle lbhc') , (<'type line;keypos ctie;keystyle lbhc';'type stick;keypos ctie;keystyle lbhc')
  rt=. '1st row';'2nd row'
  ct=. 'Col. #'
  kt=. 2 3 $ '' ; 'a' ; ('';'') ; ('' ; 'b') ; (('c1';'c2';'c3';'c4') ; (< 'd1';'d2';'d3')) ; (<'e' ; < 'f1';'f2')
  x=. 2 3 $ (< i.0) , (< i. 5) , (< i. 5) , (< (i.0) ; (i. 4)) , (< (i.0) ; (i. 3 6)) , (< (i. 5) ; (i. 6))
  dat=. 2 3 $ (mba 6) , (mba 4 5) , (< (mba 5) , (mba 5)) , (< (mba 6) , (mba 3 4)) , (< (mba 4 5) , (mba 3 6)) , (< (mba 3 5) , (mba 2 6))
  p=. mplot po;so;rt;ct;kt;x;<dat

  po=. 'title Test splot;sub 0 0 _1x _36x;sub 2 3'
  so=. 2 3 $ (<'new;type line;title Subplot #00;ycaption ycap00;xcaption xcap 00') , (<'new;type line;title Subplot #01;keypos ctie;keystyle lbhc;penstyle 0 1 2 3') , (<'new;type line;title Subplot #02';'type stick') , (<'new;type area;title Subplot #10';'type line;keypos ctie;keystyle lbhc;penstyle 0 1 2') , (<'new;type line;keypos ctie;keystyle lbhc';'type bar;title Subplot #11;keypos ctie;keystyle lbhc') , (<'new;type line;title Subplot #12;keypos ctie;keystyle lbhc';'type stick;keypos ctie;keystyle lbhc')
  kt=. 2 3 $ '' ; 'a' ; ('';'') ; ('' ; 'b') ; (('z1';'z2';'z3';'z4') ; (< 'c1';'c2';'c3')) ; (<'d' ; < 'e1';'e2')
  x=. 2 3 $ (< i.0) , (< i. 5) , (< i. 5) , (< (i.0) ; (i. 4)) , (< (i.0) ; (i. 3 6)) , (< (i. 5) ; (i. 6))
  dat=. 2 3 $ (mba 6) , (mba 4 5) , (< (mba 5) , (mba 5)) , (< (mba 6) , (mba 3 4)) , (< (mba 4 5) , (mba 3 6)) , (< (mba 3 5) , (mba 2 6))
  p=. splot po;so;kt;x;<dat

po=. 'type line'
'R C'=. RC=. 2 3

mpt=. 'mplot test'
rt=. 'rt 0' ; 'rt 1'
ct=. 'ct 0' ; 'ct 1' ; 'ct 2'
mkt=. 'kt 0' ; 'kt 1' ; 'kt 2' ; 'kt 3'
mx=. 0.1 * i. 10
mdata=. ? (4 10 , RC) $ 10
mplot_tau_ po;mpt;rt;ct;mkt;mx;<mdata

spt=. 'splot test'
so=. RC $ 'type line';'type line,marker';'type line,stick';'type bar';'type hist';'type sbar'
st=. RC $ 't 00' ; 't 01' ; 't 02' ; 't 10' ; 't 11' ; 't 12'
yt=. RC $ 'yt 00' ; 'yt 01' ; 'yt 02' ; 'yt 10' ; 'yt 11' ; 'yt 12'
xt=. RC $ 'xt 00' ; 'xt 01' ; 'xt 02' ; 'xt 10' ; 'xt 11' ; 'xt 12'
skt=. ('k #';(<'kt 0' ; 'kt 1' ; 'kt 2' ; 'kt 3')) (1 2;0 1) } RC $ a:
sx=. RC $ < _1 _0.5 0 0.5 1
sdata=. ((? 4 5 $ 100);(? 4 5 $ 100)) (0 1;1 2) } <"1 i. (RC, 5)
splot_tau_ po;spt;so;st;yt;xt;skt;sx;<sdata


)
