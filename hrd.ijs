NB. hrd.ijs
NB. Reduce a general matrix to upper Hessenberg form
NB.
NB. gehrdu  Reduce a general matrix to upper Hessenberg form
NB. gehrdl  Reduce a general matrix to upper Hessenberg form
NB.
NB. hehrd   Reduce a general matrix to Hessenberg form
NB.         (3-diagonal Hermitian matrix)
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
NB. gehrd u                                                  1
NB. Reduce a general matrix A to upper Hessenberg form H by
NB. an unitary similarity transformation:
NB.   Q' * A * Q = H
NB.
NB. Syntax:
NB.   'tau Qf H'=. gehrd A ; ss
NB. where
NB.   A   - N×N-matrix with isolated eigenvalues (see gebal)
NB.   ss  - 2-vector, corner start (left and up) and size
NB.         (width and height) of A11 (see gebalp)
NB.   HQ  - N×N-matrix, ???
NB.   tau - (N-1)-vector, ???
NB.
NB. If:
NB.   'HQ tau'=. gehrd A ; ss
NB.   >>>>>>>>>>>>>>>>>D=. diagmat s
NB.   Dinv=. %. D
NB. then
NB.   Dinv -: diagmat % s
NB.   As -: D mp Ap mp Dinv
NB.   As -: Ap (] * (% " 1)) s
NB.
NB. Notes:
NB. - gehrdu emulates xGEHD2
NB.
NB. References:
NB. [1] 
NB. [2] 

gehrdu=: 3 : 0
  'A ss'=. y
  n=. # A
  e=. +/ ss
  Tau=. (<: n) $ 0
  for_i. (<: e) ht2i ({. ss) do.
    i1=. >: i
    iosV=. < (e ht2i i1) ; i      NB. ios of rotating vector (α,x[1],...,x[n-1]) and result (β,v[1],...,v[n-1])
    betavtau=. larfg iosV { A     NB. (β[i],v[i][1],...,v[i][n-1],τ[i])
    tau=. {: betavtau             NB. τ[i]
    Tau=. tau i } Tau             NB. write τ[i] into τ in-place
    if. 0 ~: tau do.
      A=. (}: betavtau) iosV } A  NB. write back (β[i],v[i][1],...,v[i][n-1]) into A
      i1e=. e ht2i i1
      iosR=. < (i. e) ; i1e       NB. IOS updating submatrix A*H
      A=. (betavtau larfr (iosR { A)) iosR } A
      iosL=. < i1e ; (n ht2i i1)  NB. IOS updating submatrix H'*A
      A=. (betavtau larfl (iosL { A)) iosL } A
NB. smoutput 'iV' ; iosV ; 'iV { A' ; ($ iosV { A) ; (iosV { A) ; 'iR' ; iosR ; 'iR { A' ; ($ iosR { A) ; (iosR { A) ; 'iL' ; iosL ; 'iL { A' ; ($ iosL { A) ; (iosL { A)
    end.
  end.
  Tau ; (_1 & tru ; (trl1 @ }.)) A
)

gehrdu2=: 3 : 0
  'A ss'=. y
  n=. # A
  e=. +/ 'f s'=. ss                                 NB. 'from' and 'size' values define submatrix to be reduced
  Tau=. (<: n) $ 0
  ciosV=. ((>: f) j. (<: s)) , (f j. 1)             NB. cIOS of rotating vector (α,x[1],...,x[n-1]) and result (β,v[1],...,v[n-1])
  ciosR=. (0 j. e) , (e (] j. -) (>: f))            NB. cIOS of submatrix updated from the right
  ciosL=. ((>: f) j. (<: s)) , (n (] j. -) (>: f))  NB. cIOS of submatrix updated from the left
  for_i. (<: e) ht2i f do.                          NB. i=s:(f+s-1)
NB. smoutput 'ciosV' ; ciosV ; 'cios2ios ciosV' ; (cios2ios ciosV) ; 'ciosV cfrom A' ; ($ ciosV cfrom A) ; (ciosV cfrom A)
    bvt=. larfg ciosV cfrom A                       NB. btv=(β[i],v[i][1],...,v[i][n-1],τ[i])
    tau=. {: bvt                                    NB. τ[i]
    Tau=. tau i } Tau                               NB. write τ[i] into τ in-place
    if. 0 ~: tau do.
      NB. 1) write back (β[i],v[i][1],...,v[i][n-1]) into A
      NB. 2) apply rotation from the right
      NB. 3) apply rotation from the left
NB. smoutput 'cV' ; (cios2ios ciosV) ; 'cV { A' ; ($ ciosV cfrom A) ; (ciosV cfrom A) ; 'cR' ; (cios2ios ciosR) ; 'cR { A' ; ($ ciosR cfrom A) ; (ciosR cfrom A) ; 'cL' ; (cios2ios ciosL) ; 'cL { A' ; ($ ciosL cfrom A) ; (ciosL cfrom A)
      A=. bvt ([ (larfl larf1 ciosL) [ (larfr larf1 ciosR) (((ciosV camend)~ }:)~)) A
      ciosV=. ciosV + 1j_1 1
      ciosR=. ciosR + 0 1j_1
      ciosL=. ciosL + 1j_1 1j_1
    end.
  end.
  Tau ; (_1 & tru ; (trl1 @ }.)) A
)

NB. =========================================================
NB. Test suite

tgehrd=: 3 : 0
)

testhrd=: 3 : 0
)
