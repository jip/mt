NB. util.ijs
NB. Linear time-invariant (LTI) system's utilities
NB.
NB. ht         Hermitian (conjugate) transpose of table
NB. powm       raise table y to power x
NB. powsm      make report of table y powers
NB. shiftdiag  add element[s from] x to diagonal of table y
NB.
NB. rndmat     generate random matrix
NB. rndmatpe   generate random matrix with positive eigenvalues
NB. rndmatne   generate random matrix with negative eigenvalues
NB.
NB. Resources:
NB. - http://www.jsoftware.com/jwiki/...
NB. - http://www.dvgu.ru/forum/...
NB.
NB. Version: 1.0.0 2008-03-30
NB. Copyright: Igor Zhuravlov |.'ur.ugvd.ciu@rogi'
NB. License: GPL v3 or later

script_z_ '~system/packages/math/mathutil.ijs'  NB. mp
script_z_ '~system/packages/math/makemat.ijs'   NB. idmat diagmat
require '~user/projects/lapack/lapack.ijs'      NB. '~addons/math/lapack/lapack.ijs'
require '~user/projects/lapack/gesvd.ijs'       NB. need_jlapack_ 'gesvd'

coclass 'tau'

NB. =========================================================

NB. Hermitian (conjugate) transpose of table
ht=: +@|:

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
o=. 2b100 gesvd_jlapack_ rndmat y
m=. o mp d mp +|:o
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
m=. o mp d mp +|:o
)
