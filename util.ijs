NB. util.ijs
NB. Linear time-invariant (LTI) system's utilities
NB.
NB. makeP        make report of table y powers
NB. timeplot     plot the report in time-domain as multi-plot
NB. freqplot     plot the report in frequency-domain as multi-plot
NB.
NB. rndmat       generate random matrix
NB. rndmat_neig  generate random matrix with negative eigenvalues
NB.
NB. Resources:
NB. - http://www.jsoftware.com/jwiki/...
NB. - http://www.dvgu.ru/forum/...
NB.
NB. 2008-02-28 1.0.0 Igor Zhuravlov |.'ur.ugvd.ciu@rogi'

script_z_ '~system/packages/math/mathutil.ijs'  NB. mp
script_z_ '~system/packages/math/makemat.ijs'   NB. idmat diagmat
require '~user/projects/lapack/lapack.ijs'      NB. '~addons/math/lapack/lapack.ijs'
require '~user/projects/lapack/gesvd.ijs'       NB. need_jlapack_ 'gesvd'

coclass 'tau'

NB. =========================================================
NB. makeP
NB. Make report of table y powers: I y y^2 ... y^(x-1)
NB.
NB. Syntax:
NB.   P=. [x] makeP y
NB. where
NB.   y - N-by-N matrix
NB.   x >= 0, powers count, #y is default
NB.   P - x-by-N-by-N report, powers 0..(x-1) of y
NB.   N >= 0
NB.
NB. Test:
NB.    ;/ makeP 3 3 $ 5 5 0 9 0 5 1 6 7
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
NB. - memoization in use

p2b=: < @ I. @ |.                 NB. cvt bits of y to powers, then box it
pows=: p2b"1 @ #: @ i. @ [        NB. call p2b for each power x represented binary
topow=: mp/ @ (mp~ @ ] ^: [)      NB. produce table y to powers from list x, then product all
topows=: topow &. >               NB. apply dyad topow under boxes
make1=: idmat @ # &. > @ ]        NB. make boxed identity matrix of size N
repl1=: (make1) 0: } topows       NB. replace by identity matrix in 1st box
prepP=: pows > @ repl1 < @ ]      NB. form report of y ^ i. x
make0=: (0 $~ 0 , $ @ ]) " _      NB. make zero report: 0 N N $ 0
check0=: make0`prepP @. (0 ~: [)  NB. choose report type depending on x=0
makeP=: (# $: ]) :check0 M.       NB. force dyadic call: (#@]) check0 ]

NB. ---------------------------------------------------------
NB. timeplot
NB. Plot the report in time-domain as R-by-C multi-plot
NB.
NB. Synax:
NB.   [tplot[;trows[;tcols]]] timeplot [x;]data1[;data2[;...]]
NB. where
NB.   tplot - string, optional title for entire plot
NB.   trows - optional title for rows, any one of:
NB.           s            - string, prefix for row titles
NB.           s1;s2;...;sR - boxed strings, row titles
NB.   tcols - optional title for cols, any one of:
NB.           s            - string, prefix for column titles
NB.           s1;s2;...;sC - boxed strings, column titles
NB.   x     - N-vector, optional x tics, default is
NB.             (i. # data1)
NB.   datai - N-by-R-by-C report, data to plot in R-by-C
NB.           multi-plot. If more than one report is
NB.           supplied, then all vectors ((<r,c) {"2 datai)
NB.           are plotted in the same sub-plot with
NB.           coordinate (<r,c)
NB.   R >= 0
NB.   C >= 0
NB.   N >= 0
NB.
NB. Applications:
NB.   ('LTI resp.';(<'1st output ch.';'2nd output ch.');'Input ch. #') timeplot t;Y1;Y2;Y3
NB.   (<'LTI #1';'LTI #2') (timeplot & >) (< t1;X11;X12;X13) , (< t2;Y21;Y22)
NB.
NB. TODO:
NB. - check Nx|Ny|Nu == 0

timeplot=: (('Plot';'Out #';'In #')&$: :(4 : 0)) " 1 _
0$0
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
NB. rndmat_neig
NB. Generate square random matrix with negative eigenvalues
NB. in range [-10,0)
NB.
NB. Syntax:
NB.   mat=. rndmat_neig size

rndmat_neig=: 3 : 0
d=. diagmat _0.1 + _0.1 * ? y $ 100
o=. 2b100 gesvd_jlapack_ rndmat y
m=. o mp d mp +|:o
)
