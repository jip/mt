NB. tau_utils.ijs
NB. Utilities for Control Theory toolbox
NB.
NB. makeP  make report of table y powers
NB.
NB. Resources:
NB. - http://www.jsoftware.com/jwiki/...
NB. - http://www.dvgu.ru/forum/...
NB.
NB. TODO:
NB. - consider s/@:/@/g when possible
NB. - consider B is vector
NB.
NB. 2008-02-11 1.0.0 Igor Zhuravlov |.'ur.ugvd.ciu@rogi'

require '~system/packages/math/mathutil.ijs'  NB. mp
require '~system/packages/math/makemat.ijs'   NB. idmat

coclass 'tau'

NB. ===========================================================================
NB. makeP
NB. Make report of table y powers
NB.
NB. Syntax:
NB.   P=. [x] makeP y
NB. where
NB.   y - N-by-N matrix, N >= 0
NB.   x - powers count, #y is default, x >= 0
NB.   P - x-by-N-by-N report, powers 0..(x-1) of y
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

p2b=: < @ I. @ |.                 NB. cvt bits of y to powers, then box it
pows=: p2b"1 @ #: @ i. @ [        NB. call p2b for each power x represented binary
topow=: mp/ @ (mp~ @ ] ^: [)      NB. produce table y to powers from list x, then product all
topows=: topow &. >               NB. apply dyad topow under boxes
make1=: idmat @ # &. > @ ]        NB. make boxed identity matrix of size N
repl1=: (make1) 0: } topows       NB. replace by identity matrix in 1st box
prepP=: pows > @ repl1 < @ ]      NB. form report of y ^ i. x
make0=: (0 $~ 0 , $ @ ]) " _      NB. make zero report: 0 N N $ 0
check0=: make0`prepP @. (0 ~: [)  NB. choose report type depending on x=0
makeP=: (# $: ]) :check0          NB. force dyadic call: (#@]) check0 ]

