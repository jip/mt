NB. orf.ijs
NB. Orthogonal factorizations QR RQ QL LQ RRQR RRRQ RRQL RRLQ
NB.
NB. geqrf   Fact...
NB.
NB. Copyright (C) 2009 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2009-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. geqrf
NB. emulate xGEQR2
NB. 'Tau Qf R'=. geqrf A

geqrf=: 3 : 0
  k=. <./ 'm n'=. $ y
  Tau=. k $ 0
  for_i. i. k do.
    i1=. >: i
    iosV=. < (m ht2i i) ; i              NB. ios of rotating vector (α,x[1],...,x[n-1]) and result (β,v[1],...,v[n-1])
    betavtau=. larfp iosV { y            NB. (β[i],v[i][1],...,v[i][n-1],τ[i])
    tau=. {: betavtau                    NB. τ[i]
    Tau=. tau i } Tau                    NB. write τ[i] into τ in-place
    if. i1 < n do.                       NB. if (i ≤ (n-1)) ...
      y=. (}: betavtau) iosV } y         NB. write back (β[i],v[i][1],...,v[i][n-1]) into A
      iosL=. < (m ht2i i) ; (n ht2i i1)  NB. IOS updating submatrix H'*A
      y=. (((+ tau) _1 } betavtau) larfl (iosL { y)) iosL } y
    end.
  end.
  Tau ; (tru ; trl1) y
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB.*tgeqrf v test geqrf

tgeqrf=: 3 : 0
)

NB. ---------------------------------------------------------
NB.*testorf v test orthogonal factorizations

testorf=: 3 : 0
)
