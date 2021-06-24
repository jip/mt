NB. Test
NB.
NB. tmonad  Conj. to make monad to test computational monad
NB. tdyad   Conj. to make monad to test computational dyad
NB. drvevx  Dyads to compute the normalization error of
NB.         eigenvectors produced by nonsymmetric eigenvalue
NB.         problem solvers
NB. drgev   Adv. to make dyad to compute the relative
NB.         backward error of eigenvectors produced by
NB.         generalized nonsymmetric eigenvalue problem
NB.         solvers
NB. xxt01   Actors to compute the relative backward error for
NB.         the matrix reconstructed from gexxf,tzxxf output
NB. t02x    Modifiers to make dyad to compute the relative
NB.         backward error for the solution[s] computed
NB. xxt02   Dyads to compute the relative backward error for
NB.         the matrix reconstructed partially from ungxx
NB.         output
NB. t03     Dyad to compute the relative backward error for
NB.         the matrix times its inverse reconstructed
NB. xxt03   Actors to compute the relative backward error for
NB.         the matrix [partial] multiplication by unmxxxx
NB. t04x    Actors to compute the relative forward error for
NB.         the solution[s] computed
NB. xxt11   Dyads to compute the relative backward error for
NB.         the unitary (orthogonal) matrix reconstructed
NB.         from gexpf gepxf output
NB. qrt14   Checks whether X is in the row space of op(A)
NB. qrt16x  Adv. to make dyad to compute the residual for a
NB.         solution[s] computed of an overdetermined or
NB.         underdetermined system involving a matrix of full
NB.         rank, or its [conjugate-]transpose
NB. qrt171  Adv. to make dyad to compute the ratio for
NB.         zero-residual problem
NB. t211    Dyad to compute the relative backward error of
NB.         eigenvectors produced by symmetric eigenvalue
NB.         problem solvers
NB. t22x    Dyads to compute the error of eigenvectors
NB.         produced by nonsymmetric eigenvalue problem
NB.         solvers
NB. t511x   Dyads to compute the error of generalized Schur
NB.         form produced by hgexxsxx
NB. t513x   Monads to compute the error of Schur vectors
NB.         produced by hgexxsxx
NB. t52xx   Dyads to compute the error of Schur vectors
NB.         produced by tgevcxxx
NB.
NB. Version: 0.13.2 2021-06-24
NB.
NB. Copyright 2010-2021 Igor Zhuravlov
NB.
NB. This file is part of mt
NB.
NB. mt is free software: you can redistribute it and/or
NB. modify it under the terms of the GNU Lesser General
NB. Public License as published by the Free Software
NB. Foundation, either version 3 of the License, or (at your
NB. option) any later version.
NB.
NB. mt is distributed in the hope that it will be useful, but
NB. WITHOUT ANY WARRANTY; without even the implied warranty
NB. of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
NB. See the GNU Lesser General Public License for more
NB. details.
NB.
NB. You should have received a copy of the GNU Lesser General
NB. Public License along with mt. If not, see
NB. <http://www.gnu.org/licenses/>.

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. drvev
NB.
NB. Description:
NB.   Monad to compute the normalization error of
NB.   eigenvectors produced by nonsymmetric eigenvalue
NB.   problem solvers
NB.
NB. Syntax:
NB.   nerrV=. drvev V
NB. where
NB.   V     - n×n-matrix, either left or right eigenvectors
NB.   nerrV ≥ 0, the normalization error
NB.   n     ≥ 0, the size of V
NB.
NB. Formula:
NB.   err     := max(errL,errR)
NB.   vrmax[j] := max(|Re(V[i,j])|)
NB.                i | Im(V[i,j]) == 0
NB.
NB.   vmax[j]  := max(|V[:,j]|)
NB.                i
NB.
NB.   if ∃ j | vrmax[j] / vmax[j] < 1 - 2 * FP_PREC then
NB.     errV := 1 / FP_PREC
NB.   else
NB.     errV := max(min(| ||V[:,j]||_E - 1 | , 1) / FP_PREC)
NB.                j
NB.   endif
NB. where
NB.   V    - either L or R
NB.   errV - either errL or errR
NB.
NB. Notes:
NB. - models LAPACK's xDRVEV

drvev=: (% FP_PREC)"_`(>./@(FP_PREC %~ 1 ([ <. |@:-) normsc))@.((1 _2 p. FP_PREC) *./@:<: (%~&(>./@:|) ({."1 (* 0 = *) {:"1)@:+.))

NB. ---------------------------------------------------------
NB. cberrAfact
NB.
NB. Description:
NB.   Conj. to make dyad to compute the relative backward
NB.   error of matrix factorization
NB.
NB. Syntax:
NB.   berr=. (A ; normA) (normx cberrAfact mul) factors
NB. where
NB.   normx   - monad to compute matrix norm; is called as:
NB.               normM=. normx M
NB.   mul     - monad to compute Aapprox; is called as:
NB.               Aapprox=. mul factors
NB.   A       - n×n-matrix to decompose
NB.   normA   ≥ 0, the norm of A
NB.   factors - any noun, boxed factors of Aapprox
NB.   Aapprox - same shape as A, approximate A
NB.   berr    ≥ 0, the relative backward error
NB.   n       ≥ 0, the size of A and Aapprox
NB.
NB. Formula:
NB.   n := size(A)
NB.   ||A|| := max(normx(A) , FP_SFMIN)
NB.   ||F|| := normx(A - Aapprox)
NB.   if ||A|| > ||F|| then
NB.     berr := (||F|| / ||A||) / (FP_PREC * n)
NB.   elseif 1 > ||A|| then
NB.     berr := (min(||F|| , n * ||A||) / ||A||) / (FP_PREC * n)
NB.   else
NB.     berr := min(||F|| / ||A|| , n) / (FP_PREC * n)
NB.   endif
NB.
NB. Notes:
NB. - models LAPACK's xGET51(1), xGET51(2) and 1st check in
NB.   DSYT21(1) and ZHET21(1) when normx is norm1

cberrAfact=: 2 : '(u@((- 0&{::)~ v) ((% {:) <. 0 { ])`((<. */) % 1 { ])@.(1 ([ > {) ])`(% {:)@.(< {:) #@(0 {:: [) , FP_SFMIN >. 1 {:: [) % (FP_PREC * #@(0 {:: [))'

NB. ---------------------------------------------------------
NB. aberrU
NB.
NB. Description:
NB.   Adv. to make monad to compute the relative backward
NB.   error of unitary (orthogonal) matrix
NB.
NB. Syntax:
NB.   berrU=. (compI aberrU) Uapprox
NB. where
NB.   compI   - monad to compute Iapprox; is called as:
NB.               Iapprox=. compI Uapprox
NB.   Uapprox - n×n-matrix, approximate unitary (orthogonal)
NB.   Iapprox - n×n-matrix, approximate identity matrix
NB.   berrU   ≥ 0, the relative backward error for Uapprox
NB.   n       ≥ 0, the size of A, Uapprox and Iapprox
NB.
NB. Formula:
NB.   n := size(Uapprox)
NB.   berrU := min(||Iapprox - I||_1 , n) / (FP_PREC * n)
NB.
NB. Notes:
NB. - models LAPACK's 2nd check in DSYT21(1) and ZHET21(1)

aberrU=: 1 : 'norm1@(<: upddiag)@u (<. % FP_PREC * ]) #'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. tmonad
NB. tdyad
NB.
NB. Description:
NB.   Conj. to make monad to test computational verb
NB.
NB. Syntax:
NB.   vtestm=. mname tmonad        vgety`vgeto`vrcond`vferr`vberr
NB.   vtestd=. dname tdyad   vgetx`vgety`vgeto`vrcond`vferr`vberr
NB. where
NB.   vgetx  - monad to extract left argument for vd; is
NB.            called as:
NB.              argx=. vgetx y
NB.   vgety  - monad to extract right argument for vm or vd;
NB.            is called as:
NB.              argy=. vgety y
NB.   vgeto  - monad to extract output from ret;
NB.            is called as:
NB.              out=. vgeto ret
NB.   vrcond - monad to find rcond; is called as:
NB.              rcond=. vrcond y
NB.   vferr  - dyad to find ferr; is called as:
NB.              ferr=. y vferr out
NB.   vberr  - dyad to find berr; is called as:
NB.              berr=. y vberr out
NB.   mname  - literal, the name of monad vm to test
NB.   dname  - literal, the name of dyad vd to test
NB.   vtestm - monad to test monad vm and to log result:
NB.              mname rcond ferr berr time space
NB.            on the screen, in the global var TESTLOG and,
NB.            optionally, in the log file; is called as:
NB.              vtestm y
NB.   vtestd - monad to test dyad vd and to log result:
NB.              dname rcond ferr berr time space
NB.            on the screen, in the global var TESTLOG and,
NB.            optionally, in the log file; is called as:
NB.              vtestd y
NB.   vm     - monad to test; is called as:
NB.              ret=. vm argy
NB.   vd     - dyad to test; is called as:
NB.              ret=. argx vd argy
NB.   y      - some input for vtestm or vtestd
NB.   argx   - some left argument for vd
NB.   argy   - some right argument for vm or vd
NB.   ret    - some output from vm or vd
NB.   out    - rectified ret, i.e. filtered output
NB.   ferr   ≥ 0 or NaN, the relative forward error
NB.   berr   ≥ 0 or NaN, the relative backward error
NB.   rcond  ≥ 0, the estimated reciprocal of the condition
NB.            number of the input matrix; +∞ if matrix is
NB.            singular; NaN if matrix is non-square
NB.
NB. Application:
NB. - to test geqrf:
NB.     NB. to estimate rcond in 1-norm
NB.     vrcond=. (_."_)`gecon1@.(=/@$)
NB.     NB. to calc. berr, assuming:
NB.     NB.   berr := ||A - realA||_1 / (FP_EPS * ||A||_1 * m)
NB.     vberr=. ((- %&norm1 [) % FP_EPS * (norm1 * #)@[) unmqr
NB.     NB. do the job
NB.     ('geqrf' tmonad ]`]`vrcond`(_."_)`vberr) A
NB. - to test getrs:
NB.     NB. to estimate rcond in ∞-norm
NB.     vrcond=. (_."_)`geconi@.(=/@$)@(0&{::)
NB.     NB. to calc. ferr, assuming:
NB.     NB.   ferr := ||x - realx||_inf / ||realx||_inf
NB.     vferr=. ((- %&normi [) 1&{::)~
NB.     NB. to calc. componentwise berr [LUG 75], assuming:
NB.     NB.   berr := max_i(|b - A * realx|_i / (|A| * |realx| + |b|)_i)
NB.     vberr=. (mp&>/@[ |@- (0 {:: [) mp ]) >./@% (((0 {:: [) mp&| ]) + |@mp&>/@[)
NB.     NB. do the job
NB.     ('getrs' tdyad (0&{::)`(mp&>/)`]`vrcond`vferr`vberr) (A;x)

tmonad=: 2 : 0
  '`vgety vgeto vrcond vferr vberr'=. n
  try. rcond=. vrcond y catch. rcond=. _ end.
  try.
    argy=. vgety y
    try.
      't s'=. timespacex 'ret=. ' , m , ' argy'
      try.
        out=. vgeto ret
        try. ferr=. y vferr out catch. ferr=. _. end.
        try. berr=. y vberr out catch. berr=. _. end.
      catch.
        'ferr berr'=. 2 # _.
      end.
    catch.
      dbsig 3  NB. jump to upper catch block
    end.
  catch.
    'ferr berr t s'=. 4 # _.
  end.
  logline=. fmtlog_mt_ m ; rcond ; ferr ; berr ; t ; s
  logline (1!:2) 2
  wd^:IFQT 'msgs'
  logline (1!:3~ ,&LF)~^:(0 < #@]) TESTLOGFILE_mt_
  TESTLOG_mt_=: TESTLOG_mt_ , logline
  EMPTY
)

tdyad=: 2 : 0
  '`vgetx vgety vgeto vrcond vferr vberr'=. n
  try. rcond=. vrcond y catch. rcond=. _ end.
  try.
    argx=. vgetx y
    argy=. vgety y
    try.
      't s'=. timespacex 'ret=. argx ' , m , ' argy'
      try.
        out=. vgeto ret
        try. ferr=. y vferr out catch. ferr=. _. end.
        try. berr=. y vberr out catch. berr=. _. end.
      catch.
        'ferr berr'=. 2 # _.
      end.
    catch.
      dbsig 3  NB. jump to upper catch block
    end.
  catch.
    'ferr berr t s'=. 4 # _.
  end.
  logline=. fmtlog_mt_ m ; rcond ; ferr ; berr ; t ; s
  logline (1!:2) 2
  wd^:IFQT 'msgs'
  logline (1!:3~ ,&LF)~^:(0 < #@]) TESTLOGFILE_mt_
  TESTLOG_mt_=: TESTLOG_mt_ , logline
  EMPTY
)

NB. ---------------------------------------------------------
NB. drvevl
NB. drvevr
NB.
NB. Description:
NB.   Dyads to compute the normalization error of
NB.   eigenvectors produced by nonsymmetric eigenvalue
NB.   problem solvers
NB.
NB. Syntax:
NB.   errL=. trash drgevl (trash ; L     ; trash)
NB.   errR=. trash drgevr (trash ; trash ; R    )
NB. where
NB.   L    - n×n-matrix, columns with left eigenvectors
NB.   R    - n×n-matrix, columns with right eigenvectors
NB.   errL ≥ 0, float, scalar, error for L
NB.   errR ≥ 0, float, scalar, error for R

drvevl=: drvev@(1 {:: ])
drvevr=: drvev@(2 {:: ])

NB. ---------------------------------------------------------
NB. drgev
NB.
NB. Description:
NB.   Adv. to make dyad to compute the relative backward
NB.   error of eigenvectors produced by generalized
NB.   nonsymmetric eigenvalue problem solvers
NB.
NB. Syntax:
NB.   vberr=. mmul`vmul`norma`normb drgev
NB. where
NB.   mmul     - dyad to multiply matrices; is called as:
NB.                M3=. M1 mmul M2
NB.   vmul     - dyad to multiply matrix by diagonal matrix
NB.              represented as vector; is called as:
NB.                M2=. v1 vmul M1
NB.   norma    - monad to compute norm of matrix; is called
NB.              as:
NB.                normM=. norma M1
NB.   normb    - monad to compute norms of matrix rows
NB.              (columns); is called as:
NB.                v2=. normb M1
NB.   vberr    - dyad to compute berr; is called as:
NB.                berr=. (AB ; normAB) vberr (e1e2 ; V)
NB.   AB       - 2×n×n-brick, matrix pair (A,B) to
NB.              eigen-decompose
NB.   e1e2     - 2×n-matrix, laminated vectors of
NB.              generalized eigenvalues α and β
NB.   V        - either L or R
NB.   M1,M2,M3 - n×n-matrix
NB.   v1,v2    - n-vector
NB.   L        - n×n-matrix, columns with left eigenvectors
NB.   R        - n×n-matrix, columns with right eigenvectors
NB.   berr     ≥ 0, float, scalar, backward error
NB.   normAB   -:normA , normB
NB.   normA    ≥ 0, float, scalar, the norm of matrix A
NB.   normB    ≥ 0, float, scalar, the norm of matrix B
NB.
NB. Formula:
NB.   n := size(A)
NB.   berr := max(berr0,berr1)
NB.   - for ggevlxx:
NB.     ||L|| := max(||L||_inf , FP_PREC)
NB.     ||R|| := max(||R||_inf , FP_PREC)
NB.     - for L:
NB.         berr0 := (||(C2 * (C1 * E2)) * L   * A - (C2 * (C1 * E2)) * L   * B||_1 / ||L||) / FP_PREC
NB.         berr1 := normit(normitc(L) - 1) / (FP_PREC * n)
NB.     - for R:
NB.         berr0 := (||A * R^H * ((E2 * С1) * С2) - B * R^H * ((E1 * С1) * С2)||_1 / ||R||) / FP_PREC
NB.         berr1 := normit(normitc(R) - 1) / (FP_PREC * n)
NB.     - for L and R:
NB.         berr0 := max(berr0(ggevlvx),berr0(ggevlxv))
NB.         berr1 := max(berr1(ggevlvx),berr1(ggevlxv))
NB.   - for ggevuxx:
NB.     ||L|| := max(||L||_1 , FP_PREC)
NB.     ||R|| := max(||R||_1 , FP_PREC)
NB.     - for L:
NB.         berr0 := (||(C2 * (C1 * E2)) * L^H * A - (C2 * (C1 * E2)) * L^H * B||_1 / ||L||) / FP_PREC
NB.         berr1 := normit(normitc(L) - 1) / (FP_PREC * n)
NB.     - for R:
NB.         berr0 := (||A * R   * ((E2 * С1) * С2) - B * R   * ((E1 * С1) * С2)||_1 / ||R||) / FP_PREC
NB.         berr1 := normit(normitc(R) - 1) / (FP_PREC * n)
NB.     - for L and R:
NB.         berr0 := max(berr0(ggevuvx),berr0(ggevuxv))
NB.         berr1 := max(berr1(ggevuvx),berr1(ggevuxv))
NB.   C1 := (diag(coeff1))^_1
NB.   C2 := (diag(coeff2))^_1
NB.   coeff1(i) := if(sorim(α(i)) > (1/FP_SFMIN) / ||B||
NB.                or sorim(β(i)) > (1/FP_SFMIN) / ||A||
NB.                or max(sorim(α(i)),sorim(β(i))) < 1)
NB.                then max(max(sorim(α(i)),sorim(β(i))),FP_SFMIN)
NB.                else 1
NB.   coeff2(i) := max(sorim(α(i))*||B||,sorim(β(i))*||A||, FP_SFMIN)
NB.   ||A|| := max(||A||_1 , FP_SFMIN)
NB.   ||B|| := max(||B||_1 , FP_SFMIN)
NB. where
NB.   v(i) - i-th element of vector v
NB.
NB. Application:
NB. - to compute berr for L from ggevlvx:
NB.     NB. berr=. AB vberrlL (e1e2 ; L)
NB.     vberrlL=:  mp_mt_~        "2` *   `normi_mt_`normitr_mt_ drgev_mt_
NB. - to compute berr for R from ggevlxv:
NB.     NB. berr=. AB vberrlR (e1e2 ; R)
NB.     vberrlR=: (mp_mt_  ct_mt_)"2`(*"1)`normi_mt_`normitr_mt_ drgev_mt_
NB. - to compute berr for L from ggevuvx:
NB.     NB. berr=. AB vberruL (e1e2 ; L)
NB.     vberruL=: (mp_mt_~ ct_mt_)"2` *   `norm1_mt_`normitc_mt_ drgev_mt_
NB. - to compute berr for R from ggevuxv:
NB.     NB. berr=. AB vberruR (e1e2 ; R)
NB.     vberruR=: mp_mt_          "2`(*"1)`norm1_mt_`normitc_mt_ drgev_mt_
NB.
NB. Notes:
NB. - models LAPACK's xDRGEV and xGET52

drgev=: 1 : 0
  '`mmul vmul norma normb'=. m
  'e1e2 V'=. y
  n=. # V
  if. 0 = n do. 0 return. end.
  banorm=. |. FP_SFMIN >. 1 {:: x                         NB. 2-vector, float
  alfbetmax=. (% FP_SFMIN) % 1 >. banorm                  NB. 2-vector, float
  abs1ab=. sorim"1 e1e2                                   NB. 2×n-matrix, float
  abmax=. >./ abs1ab                                      NB. n-vector, float
  cond=. +./ (abs1ab > alfbetmax) , 1 > abmax             NB. n-vector, boolean
  e1e2=. e1e2 %"1 cond} 1 ,: FP_SFMIN >. abmax            NB. 2×n-matrix
  abcoeff=. (|. e1e2) %"1 >./ FP_SFMIN , abs1ab * banorm  NB. 2×n-matrix
  Err=. -/ abcoeff vmul (0 {:: x) mmul V                  NB. n×n-matrix
  errnrm=. (norma Err) % FP_PREC >. norma V               NB. scalar, float
  result1=. errnrm % FP_PREC                              NB. scalar, float
  enrmer=. normir <: normb V                              NB. scalar, float
  result2=. enrmer % FP_PREC * n                          NB. scalar, float
  result1 >. result2                                      NB. scalar, float
)

NB. ---------------------------------------------------------
NB. get01
NB.
NB. Description:
NB.   Conj. to make dyad to compute the relative backward
NB.   error for the general matrix reconstructed from
NB.   gexxxxxxx output
NB.
NB. Syntax:
NB.   berrG=. (G ; normG) (normx get01 getSize) Gapprox
NB. where
NB.   normx   - monad to compute matrix norm; is called as:
NB.               normM=. normx M
NB.   getSize - monad to get size of G; is called as:
NB.               size=. getSize G
NB.   G       - m×n-matrix, general
NB.   Gapprox - same shape as G, the approximate G
NB.   normG   ≥ 0, the norm of G
NB.   size    ∈ {m,n}, the size of G
NB.   berrG   ≥ 0, the relative backward error for Gapprox
NB.   m       ≥ 0, rows in G
NB.   n       ≥ 0, columns in G
NB.
NB. Formula:
NB.   (m,n) := shape(G)
NB.   if 0 = m or 0 = n then
NB.     berrG := 0
NB.   elseif 0 = ||G|| and 0 < ||G - Gapprox|| then
NB.     berrG := 1 / FP_EPS
NB.   else
NB.     berrG := ((||G - Gapprox|| / size) / ||G||) / FP_EPS
NB.   endif
NB. where
NB.   for getrflu1p                       : Gapprox := L * U1 * P, size := m, ||matrix|| := ||matrix||_inf
NB.   for getrfpl1u and gesvxxx and xGETRF: Gapprox := P * L1 * U, size := n, ||matrix|| := ||matrix||_1
NB.   for getrfpu1l                       : Gapprox := P * U1 * L, size := n, ||matrix|| := ||matrix||_1
NB.   for getrful1p                       : Gapprox := U * L1 * P, size := m, ||matrix|| := ||matrix||_inf
NB.
NB. Notes:
NB. - models LAPACK's xGET01

get01=: 2 : '((FP_EPS , v@]) %~/@,@,.`(%@FP_EPS)@.(</@:*@]) (1 {:: [) , u@(- 0&{::)~)`0:@.(0 e. $@])'

NB. ---------------------------------------------------------
NB. het01
NB.
NB. Description:
NB.   Dyad to compute the relative backward error for the
NB.   matrix reconstructed from hexxxxx, poxxxxx and ptxxxxx
NB.   output
NB.
NB. Syntax:
NB.   berrH=. (H ; normH) het01 Happrox
NB. where
NB.   H       - n×n-matrix, the Hermitian (symmetric),
NB.             possibly positive definite, and possibly
NB.             tridiagonal
NB.   normH   ≥ 0, the norm of H
NB.   Happrox - same shape as H, the approximate H
NB.   berrH   ≥ 0, the relative backward error for Happrox
NB.   n       ≥ 0, the size of H
NB.
NB. Formula:
NB.   n := size(H)
NB.   if 0 = n then
NB.     berrH := 0
NB.   elseif 0 = ||H||_1 and 0 < ||H - Happrox||_1 or ∃ i | 0 ≠ Im(Happrox(i,i)) then
NB.     berrH := 1 / FP_EPS
NB.   else
NB.     berrH := ((||H - Happrox||_1 / n) / ||H||_1) / FP_EPS
NB.   endif
NB. where
NB.   for hetrfpl and hesvxxx: Happrox := P * L1 * T * L1^H * P^H
NB.   for hetrfpu and        : Happrox := P * U1 * T * U1^H * P^H
NB.   for potrfl  and posvxxx: Happrox := L * L^H
NB.   for potrfu             : Happrox := U * U^H
NB.   for pttrfl  and ptsvxxx: Happrox := L1 * D * L1^H
NB.   for pttrfu             : Happrox := U1 * D * U1^H
NB.
NB. Notes:
NB. - models:
NB.   - LAPACK's DSYT01_AA, ZHET01_AA and xPTT01 with the
NB.     following differences:
NB.     - ∞-norm is used instead of 1-norm
NB.     - since diag(Happrox) may contain complex values,
NB.       then an additional check for (0 = Im(Happrox(i,i)))
NB.       is needed as same as in ZHET01 and ZPOT01
NB.   - LAPACK's xPOT01 with the following differences:
NB.     - ∞-norm is used instead of 1-norm
NB.     - a condition (0 < ||H - Happrox||_1) is checked, too
NB.       as in other xxxT01

het01=: ((FP_EPS , #@]) %~/@,@,.`(%@FP_EPS)@.(</@:*@]) (1 {:: [) , normi@(- 0&{::)~)`(%@FP_EPS)@.(0 +./@:~: 11 o. diag@])`0:@.(0 = #@])

NB. ---------------------------------------------------------
NB. hst01
NB.
NB. Description:
NB.   Adv. to make dyad to compute the relative backward
NB.   error for the Hessenberg matrix reconstructed from
NB.   gehrdx output
NB.
NB. Syntax:
NB.   berrS=. (S ; normS) (normx hst01) Sapprox
NB. where
NB.   normx   - monad to compute matrix norm; is called as:
NB.               normM=. normx M
NB.   S       - n×n-matrix, the lower or upper Hessenberg
NB.   Sapprox - same shape as S, the approximate S
NB.   normS   ≥ 0, the norm of S
NB.   berrS   ≥ 0, the relative backward error for Sapprox
NB.   n       ≥ 0, the size of S
NB.
NB. Formula:
NB.   n := size(S)
NB.   ||S|| := max(normx(S) , FP_SFMIN)
NB.   if 0 = n then
NB.     berrS := 0
NB.   else
NB.     berrS := (min(||S - Sapprox|| , ||S||) / max((FP_SFMIN * n) / FP_PREC , ||S|| * FP_PREC)) / n
NB.   endif
NB. where
NB.   for gehrdl: Sapprox := Q^H * H * Q  , ||matrix|| := ||matrix||_inf
NB.   for gehrdu: Sapprox := Q   * H * Q^H, ||matrix|| := ||matrix||_1
NB.
NB. Notes:
NB. - models LAPACK's xHST01 when normx is norm1

hst01=: 1 : '%~/@(#@] ([ , (((>. FP_PREC %~ FP_SFMIN&*)~ FP_PREC&*) {.) , 1 { ]) u@(- 0&{::)~ (] , <.) FP_SFMIN >. 1 {:: [)`0:@.(0 = #@])'

NB. ---------------------------------------------------------
NB. lqt01
NB. qlt01
NB. qrt01
NB. rqt01
NB.
NB. Description:
NB.   Dyads to compute the relative backward error for the
NB.   matrix reconstructed from gexqf geqxf output
NB.
NB. Syntax:
NB.   berrA=. (A ; normA) lqt01 (Lapprox ; Qapprox)
NB.   berrA=. (A ; normA) qlt01 (Qapprox ; Lapprox)
NB.   berrA=. (A ; normA) qrt01 (Qapprox ; Rapprox)
NB.   berrA=. (A ; normA) rqt01 (Rapprox ; Qapprox)
NB. where
NB.   A       - m×n-matrix
NB.   normA   ≥ 0, the norm of A
NB.   xapprox - approximate factors of A
NB.   berrA   ≥ 0, the relative backward error of
NB.             Q-factorization
NB.   m       ≥ 0, rows in A
NB.   n       ≥ 0, columns in A
NB.
NB. Formula:
NB.   (m,n) := shape(A)
NB.   berrA := max(berr0,berr1)
NB.   if 0 < ||A||_1 then
NB.     berr0 := ((||F1||_1 / max(1, size)) / ||A||_1) / FP_EPS
NB.   else
NB.     berr0 := 0
NB.   endif
NB.   berr1 := (||F2||_1 / max(1, size)) / FP_EPS
NB. where
NB.   for gelqf: F1 := L - A * Q^H, F2 := Q * Q^H - I, size := n
NB.   for geqlf: F1 := L - Q^H * A, F2 := Q^H * Q - I, size := m
NB.   for geqrf: F1 := R - Q^H * A, F2 := Q^H * Q - I, size := m
NB.   for gerqf: F1 := R - A * Q^H, F2 := Q * Q^H - I, size := n
NB.
NB. Notes:
NB. - lqt01 models LAPACK's xLQT01
NB. - qlt01 models LAPACK's xQLT01
NB. - qrt01 models LAPACK's xQRT01
NB. - rqt01 models LAPACK's xRQT01

lqt01=: (FP_EPS %~ (1 {:: [) %~ (norm1 % 1 >. c)@(((mp~ 0&{::)~ ct@(1&{::)) - 0 {:: ]))`0:@.(0 = 1 {:: [) >. (FP_EPS %~ (1 >. c) %~ norm1@(<: upddiag)@(mp  ct))@(1 {:: ])
qlt01=: (FP_EPS %~ (1 {:: [) %~ (norm1 % 1 >. #)@(((mp  0&{::)~ ct@(0&{::)) - 1 {:: ]))`0:@.(0 = 1 {:: [) >. (FP_EPS %~ (1 >. #) %~ norm1@(<: upddiag)@(mp~ ct))@(0 {:: ])
qrt01=: qlt01
rqt01=: lqt01

NB. ---------------------------------------------------------
NB. lpt01
NB. plt01
NB. prt01
NB. rpt01
NB.
NB. Description:
NB.   Dyads to compute the relative backward error for the
NB.   matrix reconstructed from gexpf gepxf output
NB.
NB. Syntax:
NB.   berrA=. (A ; normA) lpt01 (ip ; LQf)
NB.   berrA=. (A ; normA) plt01 (ip ; QfL)
NB.   berrA=. (A ; normA) prt01 (ip ; QfR)
NB.   berrA=. (A ; normA) rpt01 (ip ; RQf)
NB. where
NB.   A       - m×n-matrix
NB.   normA   ≥ 0, the norm of A
NB.   ip      - inversed permutation of rows or columns of A
NB.   xQf,Qfx - Q-factorization of A
NB.   berrA   ≥ 0, the relative backward error for matrix A
NB.             permuted
NB.   m       ≥ 0, rows in A
NB.   n       ≥ 0, columns in A
NB.
NB. Formula:
NB.   (m,n) := shape(A)
NB.   if 0 = m or 0 = n then
NB.     berrA := 0
NB.   else
NB.     berrA := ||F|| / (FP_EPS * max(m,n))
NB.     if 0 < ||A|| then
NB.       berrA := berrA / ||A||
NB.     endif
NB.   endif
NB. where
NB.   for gelpf: F := P * A - L * Q, ||matrix|| := ||matrix||_inf
NB.   for geplf: F := A * P - Q * L, ||matrix|| := ||matrix||_1
NB.   for geprf: F := A * P - Q * R, ||matrix|| := ||matrix||_1
NB.   for gerpf: F := P * A - R * Q, ||matrix|| := ||matrix||_inf
NB.
NB. Notes:
NB. - prt01 models LAPACK's xQPT01

lpt01=: ((1 {:: [) %~^:(0 < [) C.  ~&(0&{::) (normi % FP_EPS * >./@$)@:- (unmlqrn  trlpick        @:(}:"1))@(1 {:: ]))`0:@.(0 e. $@(0 {:: [))
plt01=: ((1 {:: [) %~^:(0 < [) C."1~&(0&{::) (norm1 % FP_EPS * >./@$)@:- (unmqlln (trlpick~ -~/@$)@  }.   )@(1 {:: ]))`0:@.(0 e. $@(0 {:: [))
prt01=: ((1 {:: [) %~^:(0 < [) C."1~&(0&{::) (norm1 % FP_EPS * >./@$)@:- (unmqrln  trupick        @  }:   )@(1 {:: ]))`0:@.(0 e. $@(0 {:: [))
rpt01=: ((1 {:: [) %~^:(0 < [) C.  ~&(0&{::) (normi % FP_EPS * >./@$)@:- (unmrqrn (trupick~ -~/@$)@:(}."1))@(1 {:: ]))`0:@.(0 e. $@(0 {:: [))

NB. ---------------------------------------------------------
NB. lzt01
NB. zlt01
NB. zrt01
NB. rzt01
NB.
NB. Description:
NB.   Dyads to compute the relative backward error for the
NB.   matrix reconstructed from tzxxf output
NB.
NB. Syntax:
NB.   berrA=. (A ; normA) lzt01 LZf
NB.   berrA=. (A ; normA) zlt01 ZfL
NB.   berrA=. (A ; normA) zrt01 ZfR
NB.   berrA=. (A ; normA) rzt01 RZf
NB. where
NB.   A       - m×n-matrix
NB.   normA   ≥ 0, the norm of A
NB.   xZf,Zfx - Z-factorization of trapezoidal part of A
NB.   berrA   ≥ 0, the relative backward error for xZf, Zfx
NB.   m       ≥ 0, rows in A
NB.   n       ≥ 0, columns in A
NB.
NB. Formula:
NB.   (m,n) := shape(A)
NB.   berrA := max(berr0,berr1)
NB.   if 0 = m or 0 = n then
NB.     berrA := 0
NB.   else
NB.     if 0 < ||A|| then
NB.       berr0 := (||F1|| / (FP_EPS * size)) / ||A||
NB.     else
NB.       berr0 := ||F1|| / (FP_EPS * size)
NB.     endif
NB.     berr1 := ||F2|| / (FP_EPS * size)
NB.   endif
NB. where
NB.   for tzlzf: F1 := A - L * Z, F2 := Z * Z^H - I, size := n, ||matrix|| := ||matrix||_1
NB.   for tzzlf: F1 := A - Z * L, F2 := Z^H * Z - I, size := m, ||matrix|| := ||matrix||_inf
NB.   for tzzrf: F1 := A - Z * R, F2 := Z^H * Z - I, size := m, ||matrix|| := ||matrix||_inf
NB.   for tzrzf: F1 := A - R * Z, F2 := Z * Z^H - I, size := n, ||matrix|| := ||matrix||_1
NB.
NB. Notes:
NB. - rzt01 models LAPACK's xRZT01 and xRZT02
NB. - shortened geometry is used:
NB.   - L (R) is square triangular min(m,n)×min(m,n)-matrix
NB.   - Z is m×n-matrix

lzt01=: ((1 {:: [) %~^:(0 < [) (norm1 % FP_EPS * 1 >. c)@((- (trlpick~ -~/@$)@(0&{::))~ (unmlzrn ((1 -  c) {."1 ({."1~ -@#)))))`0:@.(0 e. $@]) >. (norm1 % FP_EPS * 1 >. c)@((<: upddiag)~ 0 >. -~/@$)@(unmlzrc unglz)`0:@.(0 e. $)@]
zlt01=: ((1 {:: [) %~^:(0 < [) (normi % FP_EPS * 1 >. #)@((-  trlpick        @(0&{::))~ (unmzlln ((1 -~ #) {.   ({.  ~   c)))))`0:@.(0 e. $@]) >. (normi % FP_EPS * 1 >. #)@( <: upddiag             )@(unmzllc ungzl)`0:@.(0 e. $)@]
zrt01=: ((1 {:: [) %~^:(0 < [) (normi % FP_EPS * 1 >. #)@((- (trupick~ -~/@$)@(0&{::))~ (unmzrln ((1 -  #) {.   ({.  ~ -@c)))))`0:@.(0 e. $@]) >. (normi % FP_EPS * 1 >. #)@((<: upddiag)~ 0 <. -~/@$)@(unmzrlc ungzr)`0:@.(0 e. $)@]
rzt01=: ((1 {:: [) %~^:(0 < [) (norm1 % FP_EPS * 1 >. c)@((-  trupick        @(0&{::))~ (unmrzrn ((1 -~ c) {."1 ({."1~   #)))))`0:@.(0 e. $@]) >. (norm1 % FP_EPS * 1 >. c)@( <: upddiag             )@(unmrzrc ungrz)`0:@.(0 e. $)@]

NB. ---------------------------------------------------------
NB. unt01
NB.
NB. Description:
NB.   Conj. to make monad to compute the relative backward
NB.   error for the matrix reconstructed from gehrdx output
NB.
NB. Syntax:
NB.   berrU=. (normx unt01 compI) Uapprox
NB. where
NB.   normx   - monad to compute matrix norm; is called as:
NB.               normM=. normx M
NB.   compI   - monad to compute Iapprox; is called as:
NB.               Iapprox=. compI Uapprox
NB.   Uapprox - m×n-matrix, approximate unitary (orthogonal)
NB.   Iapprox - m×m- or n×n-matrix, approximate identity
NB.             matrix
NB.   normU   ≥ 0, the norm of U
NB.   berrU   ≥ 0, the relative backward error for Uapprox
NB.   m       ≥ 0, rows in Uapprox
NB.   n       ≥ 0, columns in Uapprox
NB.
NB. Formula:
NB.   (m,n) := shape(Uapprox)
NB.   if 0 = m or 0 = n then
NB.     berrU := 0
NB.   else
NB.     berrU := (||Iapprox - I|| / max(m , n)) / FP_PREC
NB.   endif
NB. where
NB.   for gehrdl: Iapprox := Q   * Q^H, ||matrix|| := ||matrix||_inf
NB.   for gehrdu: Iapprox := Q^H * Q  , ||matrix|| := ||matrix||_1
NB.
NB. Notes:
NB. - models LAPACK's DORT01, ZUNT01 when normx is norm1

unt01=: 2 : '%~/@(FP_PREC , >./@$ , u@(<: upddiag)@v)`0:@.(0 e. $@])'

NB. ---------------------------------------------------------
NB. t02m  (dyadic conj.)
NB. t02v  (dyadic adv.)
NB.
NB. Description:
NB.   Modifiers to make dyad to compute the relative backward
NB.   error for the solution[s] computed
NB.
NB. Syntax:
NB.   vberrX=. calcB t02m norm1tx
NB.   vberrx=. calcb t02v
NB. where
NB.   calcB   - dyad to compute Bapprox; is called as:
NB.               Bapprox=. Xapprox calcB A
NB.   calcb   - dyad to compute bapprox; is called as:
NB.               bapprox=. xapprox calcb A
NB.   norm1tx - monad to compute column-wise or row-wise
NB.             vector 1-taxicab-norm for list of vectors; is
NB.             called as:
NB.               normVectors=. norm1tx vectors
NB.   vberrX  - dyad to compute the relative backward error
NB.             for solutions computed; is called as:
NB.               berrX=. (A ; B ; X ; trash ; normA) vberrX Xapprox
NB.   vberrx  - dyad to compute the relative backward error
NB.             for the solution computed; is called as:
NB.               berrX=. (A ; b ; x ; trash ; normA) vberrx xapprox
NB.   A       - n×n-matrix of linear system to solve
NB.   B       - n×nrhs-matrix or nrhs×n-matrix, exact RHS
NB.   b       - n-vector, the exact RHS
NB.   Bapprox - same shape as B, approximate RHS:
NB.               Bapprox := op(A) * Xapprox  or
NB.               Bapprox := Xapprox * op(A)
NB.   bapprox - n-vector, the approximate RHS:
NB.               bapprox := op(A) * xapprox  or
NB.               bapprox := xapprox * op(A)
NB.   X       - same shape as B, exact solutions of equation:
NB.               op(A) * X = B  or
NB.               X * op(A) = B
NB.   x       - n-vector, the exact solution of equation:
NB.               op(A) * x = b  or
NB.               x * op(A) = b
NB.   Xapprox - same shape as B, approximate solutions
NB.   xapprox - n-vector, the approximate solution
NB.   normA   ≥ 0, the norm of op(A)
NB.   berrX   ≥ 0, the relative backward error
NB.   n       ≥ 0, the order of system
NB.   nrhs    ≥ 0, the number of RHS
NB.
NB. Formula:
NB.   (n,nrhs) := shape(B)
NB.   foreach i-th computed solution Xapprox from nrhs solutions do
NB.     if 0 = n or 0 = nrhs then
NB.       berrX[i] := 0
NB.     elseif 0 = ||op(A)|| or 0 = ||Xapprox|| then
NB.       berrX[i] := 1 / FP_EPS
NB.     else
NB.       berrX[i] := ((||B - Bapprox|| / ||op(A)||) / ||Xapprox||) / FP_EPS
NB.     endif
NB.   endfor
NB.   berrX := max(berrX[i])
NB. where
NB.   ||vector|| := norm1t(vector)
NB.   ||matrix|| := ||matrix||_1     when A is at left of X or either A^T or A^H is from right of X
NB.              := ||matrix||_inf   when A is at right of X or either A^T or A^H is from left of X
NB.   Bapprox    := op(A) * Xapprox  for (op(A) * X = B) equations
NB.              := Xapprox * op(A)  for (X * op(A) = B) equations
NB.
NB. Notes:
NB. - models LAPACK's xTRT02, xGET02, xGTT02, xPOT02, xPTT02

t02m=: 2 : 'max@((FP_EPS , 4 {:: [) %~/@,@,.`(%@FP_EPS)@.(0 ([ = {) ])"1 ] ,.&v       (u 0&{::)~ - 1 {:: [)`(%@FP_EPS)@.(0 = 4 {:: [)`0:@.(0 e. $@])'
t02v=: 1 : '    ((FP_EPS , 4 {:: [) %~/@,@,.`(%@FP_EPS)@.(0 ([ = {) ])   ] , &norm1tc (u 0&{::)~ - 1 {:: [)`(%@FP_EPS)@.(0 = 4 {:: [)`0:@.(0 =  #@])'

NB. ---------------------------------------------------------
NB. lqt02
NB. qlt02
NB. qrt02
NB. rqt02
NB.
NB. Description:
NB.   Dyads to compute the relative backward error for the
NB.   matrix reconstructed partially from ungxq,ungqx output
NB.
NB. Syntax:
NB.   berrQ=. (A ; normA ; LQf ; k) lqt02 Qapprox
NB.   berrQ=. (A ; normA ; QfL ; k) qlt02 Qapprox
NB.   berrQ=. (A ; normA ; QfR ; k) qrt02 Qapprox
NB.   berrQ=. (A ; normA ; RQf ; k) rqt02 Qapprox
NB. where
NB.   A       - m×n-matrix
NB.   normA   ≥ 0, the norm of A
NB.   xQf,Qfx - Q-factorization of A
NB.   k       ∈ [0,min(m,n)], an amount of elementary
NB.             reflectors taken to form Qapprox
NB.   Qapprox - m×n-matrix, approximate unitary (orthogonal),
NB.             which is defined as the product of k
NB.             elementary reflectors
NB.   berrQ   ≥ 0, the relative backward error for Qapprox
NB.   m       ≥ 0, rows in A
NB.   n       ≥ 0, columns in A
NB.
NB. Formula:
NB.   (m,n) := shape(A)
NB.   berrQ := max(berr0,berr1)
NB.   if 0 < ||A||_1 then
NB.     berr0 := ((||F1||_1 / max(1, size)) / ||A||_1) / FP_EPS
NB.   else
NB.     berr0 := 0
NB.   endif
NB.   berr1 := (||F2||_1 / max(1, size)) / FP_EPS
NB. where
NB.   k ∊ {0, 1, min(m,n)/2, min(m,n)}
NB.   for unglq: F1 := L(0  :k-1,0  :m-1) - A(0  :k-1,0:n-1)   * Qapprox         ^H, F2 := Qapprox   * Qapprox^H - I, size := n, Qapprox := H(k-1)' * ... * H(0  )'
NB.   for ungql: F1 := L(0  :n-1,n-k:n-1) - Qapprox         ^H * A(0:m-1,n-k:n-1)  , F2 := Qapprox^H * Qapprox   - I, size := m, Qapprox := H(k-1)  * ... * H(0  )
NB.   for ungqr: F1 := R(0  :n-1,0  :k-1) - Qapprox         ^H * A(0:m-1,0  :k-1)  , F2 := Qapprox^H * Qapprox   - I, size := m, Qapprox := H(0  )  * ... * H(k-1)
NB.   for ungrq: F1 := R(m-k:m-1,0  :m-1) - A(m-k:m-1,0:n-1)   * Qapprox         ^H, F2 := Qapprox   * Qapprox^H - I, size := n, Qapprox := H(0  )' * ... * H(k-1)'
NB.
NB. Notes:
NB. - m≤n for LQ and RQ, m≥n for QL and QR
NB. - lqt02 models LAPACK's xLQT02
NB. - qlt02 models LAPACK's xQLT02
NB. - qrt02 models LAPACK's xQRT02
NB. - rqt02 models LAPACK's xRQT02

lqt02=: %~/@(FP_EPS , (1 {:: [) , (1 >. c@(0 {:: [)) , norm1@((   3&{::  {.    trl        @            (2&{::))@[ - ((mp~    3&{::  {.    0&{:: )~ ct)))`0:@.(0 = 1 {:: [) >. FP_EPS %~ ((% 1 >. c@(0&{::))~ norm1@(<: upddiag)@(mp  ct))
qlt02=: %~/@(FP_EPS , (1 {:: [) , (1 >. #@(0 {:: [)) , norm1@((-@(3&{::) {."1 (trl~ -~/@$)@            (2&{::))@[ - ((mp  -@(3&{::) {."1 (0&{::))~ ct)))`0:@.(0 = 1 {:: [) >. FP_EPS %~ ((% 1 >. #@(0&{::))~ norm1@(<: upddiag)@(mp~ ct))
qrt02=: %~/@(FP_EPS , (1 {:: [) , (1 >. #@(0 {:: [)) , norm1@((   3&{::  {."1  tru        @            (2&{::))@[ - ((mp     3&{::  {."1 (0&{::))~ ct)))`0:@.(0 = 1 {:: [) >. FP_EPS %~ ((% 1 >. #@(0&{::))~ norm1@(<: upddiag)@(mp~ ct))
rqt02=: %~/@(FP_EPS , (1 {:: [) , (1 >. c@(0 {:: [)) , norm1@((-@(3&{::) {.   (tru~ -~/@$)@            (2&{::))@[ - ((mp~ -@(3&{::) {.    0&{:: )~ ct)))`0:@.(0 = 1 {:: [) >. FP_EPS %~ ((% 1 >. c@(0&{::))~ norm1@(<: upddiag)@(mp  ct))

NB. ---------------------------------------------------------
NB. lzt02
NB. zlt02
NB. zrt02
NB. rzt02
NB.
NB. Description:
NB.   Dyads to compute the relative backward error for the
NB.   matrix reconstructed partially from ungxz,ungzx output
NB.
NB. Syntax:
NB.   berrZ=. (A ; normA ; LZf ; k) lzt02 Zapprox
NB.   berrZ=. (A ; normA ; ZfL ; k) zlt02 Zapprox
NB.   berrZ=. (A ; normA ; ZfR ; k) zrt02 Zapprox
NB.   berrZ=. (A ; normA ; RZf ; k) rzt02 Zapprox
NB. where
NB.   A       - m×n-matrix
NB.   normA   ≥ 0, the norm of A
NB.   xZf,Zfx - Z-factorization of trapezoidal part of A
NB.   k       ∈ [0,min(m,n)], an amount of elementary
NB.             reflectors taken to form Zapprox
NB.   Zapprox - m×n-matrix, approximate unitary (orthogonal),
NB.             which is defined as the product of k
NB.             elementary reflectors
NB.   berrZ   ≥ 0, the relative backward error for Zapprox
NB.   m       ≥ 0, rows in A
NB.   n       ≥ 0, columns in A
NB.
NB. Formula:
NB.   (m,n) := shape(A)
NB.   berrZ := max(berr0,berr1)
NB.   if 0 < ||A||_1 then
NB.     berr0 := ((||F1||_1 / max(1, size)) / ||A||_1) / FP_EPS
NB.   else
NB.     berr0 := 0
NB.   endif
NB.   berr1 := (||F2||_1 / max(1, size)) / FP_EPS
NB. where
NB.   k ∊ {0, 1, min(m,n)/2, min(m,n)}
NB.   for unglz: F1 := A(m-k-1:m-1,n-k-1:n-1) - L(m-k-1:m-1,m-k-1:m-1) *  Zapprox(m-k-1:m-1,n-k-1:n-1), F2 := Zapprox   * Zapprox^H - I, size := n, Zapprox := H(k-1)' * ... * H(0  )'
NB.   for ungzl: F1 := A(0    :k-1,0    :k-1) - L(0    :k-1,0    :k-1) *~ Zapprox(0    :k-1,0    :k-1), F2 := Zapprox^H * Zapprox   - I, size := m, Zapprox := H(k-1)  * ... * H(0  )
NB.   for ungzr: F1 := A(m-k-1:m-1,n-k-1:n-1) - R(n-k-1:n-1,n-k-1:n-1) *~ Zapprox(m-k-1:m-1,n-k-1:n-1), F2 := Zapprox^H * Zapprox   - I, size := m, Zapprox := H(0  )  * ... * H(k-1)
NB.   for ungrz: F1 := A(0    :k-1,0    :k-1) - R(0    :k-1,0    :k-1) *  Zapprox(0    :k-1,0    :k-1), F2 := Zapprox   * Zapprox^H - I, size := n, Zapprox := H(0  )' * ... * H(k-1)'
NB.
NB. Notes:
NB. - m≤n for LZ and RZ, m≥n for ZL and ZR
NB. - shortened geometry is used:
NB.   - L (R) is square triangular min(m,n)×min(m,n)-matrix
NB.   - Z is m×n-matrix

lzt02=: %~/@(FP_EPS , (1 {:: [) , (1 >. c@(0 {:: [)) , norm1@(((0 {:: [) trlpick@({.~ 2 # -@#) ]) - (2 {:: [) (({.~ 2 # -@#) mp  ({."1~ -@#)@]) ]))`0:@.(0 = 1 {:: [) >. FP_EPS %~ ((% 1 >. c@(0&{::))~ norm1@(<: upddiag)@(mp  ct))
zlt02=: %~/@(FP_EPS , (1 {:: [) , (1 >. #@(0 {:: [)) , norm1@(((0 {:: [) trlpick@({.~ 2 #   c) ]) - (2 {:: [) (({.~ 2 #   c) mp~ ({.  ~   c)@]) ]))`0:@.(0 = 1 {:: [) >. FP_EPS %~ ((% 1 >. #@(0&{::))~ norm1@(<: upddiag)@(mp~ ct))
zrt02=: %~/@(FP_EPS , (1 {:: [) , (1 >. #@(0 {:: [)) , norm1@(((0 {:: [) trupick@({.~ 2 # -@c) ]) - (2 {:: [) (({.~ 2 # -@c) mp~ ({.  ~ -@c)@]) ]))`0:@.(0 = 1 {:: [) >. FP_EPS %~ ((% 1 >. #@(0&{::))~ norm1@(<: upddiag)@(mp~ ct))
rzt02=: %~/@(FP_EPS , (1 {:: [) , (1 >. c@(0 {:: [)) , norm1@(((0 {:: [) trupick@({.~ 2 #   #) ]) - (2 {:: [) (({.~ 2 #   #) mp  ({."1~   #)@]) ]))`0:@.(0 = 1 {:: [) >. FP_EPS %~ ((% 1 >. c@(0&{::))~ norm1@(<: upddiag)@(mp  ct))

NB. ---------------------------------------------------------
NB. t03
NB.
NB. Description:
NB.   Dyad to compute the relative backward error for the
NB.   matrix times its inverse reconstructed
NB.
NB. Syntax:
NB.   berriA=. (A ; rcondA ; normA) t03 iAapprox
NB. where
NB.   A        - n×n-matrix, general or triangular or
NB.              Hermitian (symmetric), possibly positive
NB.              definite, possibly tridiagonal
NB.   rcondA   ≥ 0, the reciprocal of the condition number of
NB.              A
NB.   normA    ≥ 0, the norm of A
NB.   iAapprox - n×n-matrix, the approximate of A inverse
NB.   berriA   ≥ 0, the relative backward error for iAapprox
NB.   n        ≥ 0, the size of A and iAapprox
NB.
NB. Formula:
NB.   n := size(A)
NB.   if 0 = n then
NB.     berriA := 0
NB.   elseif 0 = ||A||_1 or 0 = ||A^_1||_1 then
NB.     berriA := 1 / FP_EPS
NB.   else
NB.     berriA := ((||A * A^_1 - I||_1 * rcond(A)) / n) / FP_EPS
NB.   endif
NB. where
NB.   for trtril   : A := L                      , rcond(A) := trlcon1 (L)
NB.   for trtril1  : A := L1                     , rcond(A) := trl1con1(L1)
NB.   for trtriu   : A := U                      , rcond(A) := trucon1 (U)
NB.   for trtriu1  : A := U1                     , rcond(A) := tru1con1(U1)
NB.   for getrilu1p: A := L * U1 * P             , rcond(A) := gecon1  (A)
NB.   for getripl1u: A := P * L1 * U             , rcond(A) := gecon1  (A)
NB.   for getripu1l: A := P * U1 * L             , rcond(A) := gecon1  (A)
NB.   for getriul1p: A := U * L1 * P             , rcond(A) := gecon1  (A)
NB.   for hetripl  : A := P * L1 * T * L1^H * P^H, rcond(A) := hecon1  (A)
NB.   for hetripu  : A := P * U1 * T * U1^H * P^H, rcond(A) := hecon1  (A)
NB.   for potril   : A := L * L^H                , rcond(A) := pocon1  (A)
NB.   for potriu   : A := U * U^H                , rcond(A) := pocon1  (A)
NB.   for pttril   : A := L1 * D * L1^H          , rcond(A) := ptcon1  (A)
NB.   for pttriu   : A := U1 * D * U1^H          , rcond(A) := ptcon1  (A)
NB.
NB. Notes:
NB. - models LAPACK's xTRT01, xGET03, xPOT03

t03=: %~/@(FP_EPS , #@] , norm1@(<: upddiag)@(mp~ 0&{::)~ * 1 {:: [)`(%@FP_EPS)@.(0 e. ((, 2&{::)~ norm1))`0:@.(0 = #@])

NB. ---------------------------------------------------------
NB. lqt03
NB. qlt03
NB. qrt03
NB. rqt03
NB.
NB. Description:
NB.   Adv. to make dyad to compute the relative backward
NB.   error for the matrix partial multiplication by unmxxxx
NB.
NB. Syntax:
NB.   berrP=. (C ; normC ; LQf ; k) (compP lqt03) Papprox
NB.   berrP=. (C ; normC ; QfL ; k) (compP qlt03) Papprox
NB.   berrP=. (C ; normC ; QfR ; k) (compP qrt03) Papprox
NB.   berrP=. (C ; normC ; RQf ; k) (compP rqt03) Papprox
NB. where
NB.   compP   - dyad to compute P; is called as:
NB.               P=. C compP Q
NB.   C       - m×n-matrix when (C * op(Q)), or n×m-matrix
NB.             when (op(Q) * C)
NB.   normC   ≥ 0, the norm of C
NB.   xQf,Qfx - matrix to extract Q in factorized form
NB.   k       ∈ [0,min(m,n)], an amount of elementary
NB.             reflectors taken to form Q
NB.   Papprox - m×n-matrix or n×m-matrix, approximate P
NB.   berrP   ≥ 0, the relative backward error for P
NB.   Q       - n×n-matrix, unitary (orthogonal), which is
NB.             defined as the product of k elementary
NB.             reflectors of order n
NB.   op(Q)   = Q for unmxxxn or op(Q) = Q^H for unmxxxc
NB.   P       - m×n-matrix (C * op(Q)), or n×m-matrix
NB.             (op(Q) * C), the matrix inner product
NB.   m       ≥ 0, rows or columns in C
NB.   n       ≥ 0, columns or rows in C
NB.
NB. Formula:
NB.   n := size(Q)
NB.   m := columns(C) if side='L' or rows(C) if side='R'
NB.   if 0 < ||C||_1 then
NB.     berrP := ||F1||_1 / (max(1, n) * ||C||_1 * FP_EPS)
NB.   else
NB.     berrP := ||F1||_1 / (max(1, n) *           FP_EPS)
NB.   endif
NB. where
NB.   k ∊ {0, 1, min(m,n)/2, min(m,n)}
NB.   F1 := Papprox - C compP Q(k)
NB.   for unmlqxx: Q(k) := H(k-1)' * ... * H(0  )'
NB.   for unmqlxx: Q(k) := H(k-1)  * ... * H(0  )
NB.   for unmqrxx: Q(k) := H(0  )  * ... * H(k-1)
NB.   for unmrqxx: Q(k) := H(0  )' * ... * H(k-1)'
NB.
NB. Notes:
NB. - m≤n for LQ and RQ, m≥n for QL and QR
NB. - xxt03 models corresp. LAPACK's xxxT03 with he following
NB.   difference: only a sole test is performed, either
NB.   (Q * C), (Q^H * C), (C * Q) or (C * Q^H)
NB.   - lqt03 models LAPACK's xLQT03
NB.   - qlt03 models LAPACK's xQLT03
NB.   - qrt03 models LAPACK's xQRT03
NB.   - rqt03 models LAPACK's xRQT03
NB.
NB. Application:
NB. - implement LAPACK's xLQT03:
NB.     result1=. (    C  ; (norm1 C) ; LQf ; k) ((mp~ ] ) lqt03) Papprox
NB.     result2=. ((|: C) ; (normi C) ; LQf ; k) ((mp  ] ) lqt03) Papprox
NB.     result3=. (    C  ; (norm1 C) ; LQf ; k) ((mp~ ct) lqt03) Papprox
NB.     result4=. ((|: C) ; (normi C) ; LQf ; k) ((mp  ct) lqt03) Papprox
NB.     result=. result1 , result2 , result3 , result4
NB. - implement LAPACK's xQLT03:
NB.     result1=. (    C  ; (norm1 C) ; QfL ; k) ((mp~ ] ) qlt03) Papprox
NB.     result2=. ((|: C) ; (normi C) ; QfL ; k) ((mp  ] ) qlt03) Papprox
NB.     result3=. (    C  ; (norm1 C) ; QfL ; k) ((mp~ ct) qlt03) Papprox
NB.     result4=. ((|: C) ; (normi C) ; QfL ; k) ((mp  ct) qlt03) Papprox
NB.     result=. result1 , result2 , result3 , result4
NB. - implement LAPACK's xQRT03:
NB.     result1=. (    C  ; (norm1 C) ; QfR ; k) ((mp~ ] ) qrt03) Papprox
NB.     result2=. ((|: C) ; (normi C) ; QfR ; k) ((mp  ] ) qrt03) Papprox
NB.     result3=. (    C  ; (norm1 C) ; QfR ; k) ((mp~ ct) qrt03) Papprox
NB.     result4=. ((|: C) ; (normi C) ; QfR ; k) ((mp  ct) qrt03) Papprox
NB.     result=. result1 , result2 , result3 , result4
NB. - implement LAPACK's xRQT03:
NB.     result1=. (    C  ; (norm1 C) ; RQf ; k) ((mp~ ] ) rqt03) Papprox
NB.     result2=. ((|: C) ; (normi C) ; RQf ; k) ((mp  ] ) rqt03) Papprox
NB.     result3=. (    C  ; (norm1 C) ; RQf ; k) ((mp~ ct) rqt03) Papprox
NB.     result4=. ((|: C) ; (normi C) ; RQf ; k) ((mp  ct) rqt03) Papprox
NB.     result=. result1 , result2 , result3 , result4

lqt03=: 1 : 'norm1@(- 0&{:: u    3&{::  (<:@c@] unglq {.  )                              2&{:: )~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. <:@c@(2 {:: [)'
qlt03=: 1 : 'norm1@(- 0&{:: u -@(3&{::) (<:@#@] ungql {."1)                              2&{:: )~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. <:@#@(2 {:: [)'
qrt03=: 1 : 'norm1@(- 0&{:: u    3&{::  (<:@#@] ungqr {."1)                              2&{:: )~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. <:@#@(2 {:: [)'
rqt03=: 1 : 'norm1@(- 0&{:: u -@(3&{::) (<:@c@] ungrq {.  )                              2&{:: )~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. <:@c@(2 {:: [)'

NB. ---------------------------------------------------------
NB. lzt03
NB. zlt03
NB. zrt03
NB. rzt03
NB.
NB. Description:
NB.   Adv. to make dyad to compute the relative backward
NB.   error for the matrix partial multiplication by unmxxxx
NB.
NB. Syntax:
NB.   berrP=. (C ; normC ; LZf ; k) (compP lzt03) Papprox
NB.   berrP=. (C ; normC ; ZfL ; k) (compP zlt03) Papprox
NB.   berrP=. (C ; normC ; ZfR ; k) (compP zrt03) Papprox
NB.   berrP=. (C ; normC ; RZf ; k) (compP rzt03) Papprox
NB. where
NB.   compP   - dyad to compute P; is called as:
NB.               P=. C compP Z
NB.   C       - m×n-matrix when (C * op(Z)), or n×m-matrix
NB.             when (op(Z) * C)
NB.   normC   ≥ 0, the norm of C
NB.   xZf,Zfx - matrix to extract Z in factorized form
NB.   k       ∈ [0,min(m,n)], an amount of elementary
NB.             reflectors taken to form Z
NB.   Papprox - m×n-matrix or n×m-matrix, approximate P
NB.   berrP   ≥ 0, the relative backward error for P
NB.   Z       - n×n-matrix, unitary (orthogonal), which is
NB.             defined as the product of k elementary
NB.             reflectors of order n
NB.   op(Z)   = Z for unmxxxn or op(Z) = Z^H for unmxxxc
NB.   P       - m×n-matrix (C * op(Z)), or n×m-matrix
NB.             (op(Z) * C), the matrix inner product
NB.   m       ≥ 0, rows or columns in C
NB.   n       ≥ 0, columns or rows in C
NB.
NB. Formula:
NB.   n := size(Z)
NB.   m := columns(C) if side='L' or rows(C) if side='R'
NB.   if 0 < ||C||_1 then
NB.     berrP := ||F1||_1 / (max(1, n) * ||C||_1 * FP_EPS)
NB.   else
NB.     berrP := ||F1||_1 / (max(1, n) *           FP_EPS)
NB.   endif
NB. where
NB.   k ∊ {0, 1, min(m,n)/2, min(m,n)}
NB.   F1 := Papprox - C compP Z(k)
NB.   for unmlzxx: Z(k) := H(k-1)' * ... * H(0  )'
NB.   for unmzlxx: Z(k) := H(k-1)  * ... * H(0  )
NB.   for unmzrxx: Z(k) := H(0  )  * ... * H(k-1)
NB.   for unmrzxx: Z(k) := H(0  )' * ... * H(k-1)'
NB.
NB. Notes:
NB. - m≤n for LZ and RZ, m≥n for ZL and ZR
NB. - Zf is feeded to ungxx instead of xZf/Zfx

lzt03=: 1 : 'norm1@(- 0&{:: u -@(3&{::) (<:@c@] unglz {.  ) (((}."1~ -) ,.  idmat@]) #)@(2&{::))~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. <:@c@(2 {:: [)'
zlt03=: 1 : 'norm1@(- 0&{:: u    3&{::  (<:@#@] ungzl {."1) (( }.  ~    , ~ idmat@]) c)@(2&{::))~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. <:@#@(2 {:: [)'
zrt03=: 1 : 'norm1@(- 0&{:: u -@(3&{::) (<:@#@] ungzr {."1) (((}.  ~ -) ,   idmat@]) c)@(2&{::))~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. <:@#@(2 {:: [)'
rzt03=: 1 : 'norm1@(- 0&{:: u    3&{::  (<:@c@] ungrz {.  ) (( }."1~    ,.~ idmat@]) #)@(2&{::))~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. <:@c@(2 {:: [)'

NB. ---------------------------------------------------------
NB. hst03l
NB. hst03u
NB.
NB. Description:
NB.   Dyads to compute the relative backward error for the
NB.   matrix multiplication by unmhrxxx
NB.
NB. Syntax:
NB.   berrP=. (C ; normC ; HlQf ; P) hst03l Papprox
NB.   berrP=. (C ; normC ; HuQf ; P) hst03u Papprox
NB. where
NB.   C       - m×n-matrix when (C * op(Q)), or n×m-matrix
NB.             when (op(Q) * C)
NB.   normC   ≥ 0, the norm of C
NB.   HlQf    - n×(n+1)-matrix with packed Qf (see gehrdl)
NB.   HuQf    - (n+1)×n-matrix with packed Qf (see gehrdu)
NB.   P       - m×n-matrix (C * op(Q)), or n×m-matrix
NB.             (op(Q) * C), the matrix inner product
NB.   Papprox - n×n-matrix, approximate P
NB.   berrP   ≥ 0, the relative backward error for P
NB.   Q       - n×n-matrix, the unitary (orthogonal), which
NB.             is defined as the product of n elementary
NB.             reflectors of order n
NB.   op(Q)   = Q for unmhrxxn or op(Q) = Q^H for unmhrxxc
NB.   m       ≥ 0, rows or columns in C
NB.   n       ≥ 0, columns or rows in C
NB.
NB. Formula:
NB.   n := size(Q)
NB.   if 0 < ||C||_1 then
NB.     berrP := ||Papprox - P||_1 / (max(1, n) * ||C||_1 * FP_EPS)
NB.   else
NB.     berrP := ||Papprox - P||_1 / (max(1, n) *           FP_EPS)
NB.   endif

hst03l=: norm1@(- 3&{::)~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. #@(2 {:: [)
hst03u=: norm1@(- 3&{::)~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. c@(2 {:: [)

NB. ---------------------------------------------------------
NB. t04m  (dyadic adv.)
NB. t04v  (dyad)
NB.
NB. Description:
NB.   Actors to compute the relative forward error for
NB.   the solution[s] computed
NB.
NB. Syntax:
NB.   ferrX=. (A ; B ; X ; rcondA ; trash) (normitx t04m) Xapprox
NB.   ferrX=. (A ; b ; x ; rcondA ; trash)          t04v  xapprox
NB. where
NB.   normitx - monad to compute column-wise or row-wise
NB.             vector inf-taxicab-norm for list of vectors;
NB.             is called as:
NB.               normVectors=. normitx vectors
NB.   A       - n×n-matrix of linear system to solve
NB.   B       - n×nrhs-matrix or nrhs×n-matrix, exact RHS
NB.   b       - n-vector, the exact RHS
NB.   X       - same shape as B, exact solutions of equation:
NB.               op(A) * X = B  or
NB.               X * op(A) = B
NB.   x       - n-vector, the exact solution of equation:
NB.               op(A) * x = b  or
NB.               x * op(A) = b
NB.   Xapprox - same shape as B, approximate solutions
NB.   xapprox - n-vector, the approximate solution
NB.   rcondA  ≥ 0, the reciprocal of the condition number of
NB.             op(A)
NB.   ferrX   ≥ 0, the relative forward error
NB.   n       ≥ 0, the order of system
NB.   nrhs    ≥ 0, the number of RHS
NB.
NB. Formula:
NB.   (n,nrhs) := shape(B)
NB.   foreach i-th pair (X,Xapprox) from nrhs solutions do
NB.     if 0 = n or 0 = nrhs then
NB.       ferrX[i] := 0
NB.     elseif 0 = rcond(op(A)) or 0 = ||X|| and 0 < ||X-Xapprox|| then
NB.       ferrX[i] := 1 / FP_EPS
NB.     else
NB.       ferrX[i] := (||X - Xapprox|| / ||X||) * rcond(op(A))
NB.     endif
NB.   endfor
NB.   ferrX := max(ferrX[i])
NB.   if 1 > ferrX * FP_EPS then
NB.     ferrX := ferrX / FP_EPS
NB.   endif
NB. where
NB.   ||vector|| := normit(vector)
NB.
NB. Notes:
NB. - models LAPACK's xGET04

t04m=: 1 : '%&FP_EPS^:(1 > FP_EPS&*)@max@((3 {:: [) >/@:*@}.`(*`%/ ,: %@FP_EPS)}@, ] (- ,:&u       ]) 2 {:: [)`(%@FP_EPS)@.(0 = 3 {:: [)`0:@.(0 e. $@])'
t04v=:      %&FP_EPS^:(1 > FP_EPS&*)@    ((3 {:: [) *`%/@,`(%@FP_EPS)@.(>/@:*@])   ] (- , &normitr ]) 2 {:: [)`(%@FP_EPS)@.(0 = 3 {:: [)`0:@.(0 =  #@])

NB. ---------------------------------------------------------
NB. lqt11
NB. qlt11
NB. qrt11
NB. rqt11
NB.
NB. Description:
NB.   Dyad to compute the relative backward error for the
NB.   unitary (orthogonal) matrix reconstructed from gexpf
NB.   gepxf output
NB.
NB. Syntax:
NB.   berrQ=. trash0 lqt11 (trash1 ; LQf)
NB.   berrQ=. trash0 qlt11 (trash1 ; QfL)
NB.   berrQ=. trash0 qrt11 (trash1 ; QfR)
NB.   berrQ=. trash0 rqt11 (trash1 ; RQf)
NB. where
NB.   xQf,Qfx - output of gexpf or gepxf, contains Qf
NB.   berrQ   ≥ 0, the relative backward error for Qf
NB.
NB. Formula:
NB.   (m,n) := shape(LQf or QfL or QfR or RQf)
NB.   berrQ := ||F|| / (FP_EPS * size)
NB. where
NB.   for gelpf: F := Q * Q^H - I, size := n, ||matrix|| := ||matrix||_inf
NB.   for geplf: F := Q^H * Q - I, size := m, ||matrix|| := ||matrix||_1
NB.   for geprf: F := Q^H * Q - I, size := m, ||matrix|| := ||matrix||_1
NB.   for gerpf: F := Q * Q^H - I, size := n, ||matrix|| := ||matrix||_inf
NB.
NB. Notes:
NB. - qrt11 models LAPACK's xQRT11 with the following
NB.   difference in Q generation method:
NB.   - LAPACK uses two consequtive calls to DORMQR/ZUNMQR
NB.   - mt uses single dyadic call to ungxx

lqt11=: ((normi@(<: upddiag)@(mp  ct)@unglq~ (% FP_EPS&*) ]) <:@c)@(1 {:: ])
qlt11=: ((norm1@(<: upddiag)@(mp~ ct)@ungql~ (% FP_EPS&*) ]) <:@#)@(1 {:: ])
qrt11=: ((norm1@(<: upddiag)@(mp~ ct)@ungqr~ (% FP_EPS&*) ]) <:@#)@(1 {:: ])
rqt11=: ((normi@(<: upddiag)@(mp  ct)@ungrq~ (% FP_EPS&*) ]) <:@c)@(1 {:: ])

NB. ---------------------------------------------------------
NB. qrt14
NB.
NB. Description:
NB.   Checks whether X is in the row space of op(A)
NB. where
NB.   op(A) is A or A^T or A^H
NB.
NB. Syntax:
NB.   errX=. (A ; B ; trash ; normA) qrt14 Xapprox
NB. where
NB.   A       - m×n-matrix of full rank
NB.   B       - m×nrhs-matrix or n×nrhs-matrix or m-vector or
NB.             n-vector, exact RHS
NB.   normA   ≥ 0, the norm of op(A)
NB.   Xapprox - n×nrhs-matrix or m×nrhs-matrix or n-vector or
NB.             m-vector, approximate solutions of equation:
NB.               op(A) * X = B
NB.   errX    ≥ 0, the error
NB.   m       ≥ 0, the number of rows in A
NB.   n       ≥ 0, the number of columns in A
NB.   nrhs    ≥ 0, the number of RHS and the number of
NB.             columns in B, Bapprox, X and Xapprox

NB. Formula (for op(A) = A):
NB.   In: A, Xapprox
NB.   Out: a measure of distance from X to row space of op(A)
NB.   1) scale both X and A such that their norms are in the
NB.      range [sqrt(FP_EPS), 1/sqrt(FP_EPS)]
NB.   2) compute a Q-factorization:
NB.      2.1) if op(A) = A^T or A^H (which implies m≥n) then:
NB.             Q * R = [A,X]
NB.      2.2) else (which implies op(A)=A, m<n):
NB.             L * Q = op([A,X^H])
NB.   3) T := trailing triangle
NB.   4) return ||T|| / (FP_EPS * max(m,n,nrhs))
NB.
NB. Notes:
NB. - is called for overdetermined system, one of:
NB.   - op(A) = A, m < n
NB.   - op(A) = A^T or A^H, m ≥ n
NB. - models LAPACK's xQRT14

qrt14=: ((((FP_EPS * >./@(, $))~ c) * normm@((tru@}:;.0~ (_ ,:~ 2 # c))@geqr2@(,. , 0:)`((trl@:(}:"1);.0~ (_ ,:~ 2 # #))@gelq2@((, ct) ,. 0:))@.(</@$@[))&((scl~ ,&1)~^:(0 < [)~ normm))~ 0&{::)~`0:@.((0 e. (, $@(0&{::  )))~ c)

NB. ---------------------------------------------------------
NB. qrt16m
NB. qrt16v
NB.
NB. Description:
NB.   Adv. to make dyad to compute the residual for a
NB.   solution[s] computed of an overdetermined or
NB.   underdetermined system involving a matrix of full rank,
NB.   or its [conjugate-]transpose
NB.
NB. Syntax:
NB.   resX=. (A ; B ; trash ; normA) (calcB qrt16m) Xapprox
NB.   resx=. (A ; b ; trash ; normA) (calcb qrt16v) xapprox
NB. where
NB.   calcB   - dyad to compute Bapprox; is called as:
NB.               Bapprox=. Xapprox calcB A
NB.   calcb   - dyad to compute bapprox; is called as:
NB.               bapprox=. xapprox calcb A
NB.   A       - m×n-matrix of full rank
NB.   B       - nrhs×n-matrix or nrhs×m-matrix, exact RHS
NB.   b       - n-vector or m-vector, the exact RHS
NB.   Bapprox - same shape as B, approximate RHS
NB.   bapprox - same shape as b, the approximate RHS
NB.   Xapprox - n×nrhs-matrix or m×nrhs-matrix, approximate
NB.             solutions of equation:
NB.               op(A) * X = B
NB.   xapprox - n-vector or m-vector, the approximate
NB.             solution of equation:
NB.               op(A) * x = b
NB.   normA   ≥ 0, the norm of op(A)
NB.   resX    ≥ 0, the max of residuals of solutions
NB.   resx    ≥ 0, the residual of the solution
NB.   m       ≥ 0, the number of rows in A
NB.   n       ≥ 0, the number of columns in A
NB.   nrhs    ≥ 0, the number of RHS
NB.
NB. Formula:
NB.   m := rows(A)
NB.   (n,nrhs) := shape(X)
NB.   if 0 = m or 0 = n or 0 = nrhs then
NB.     resX := 0
NB.   else
NB.     foreach i-th computed solution Xapprox from nrhs solutions do
NB.       if 0 = ||op(A)|| and 0 = ||B - Bapprox|| then
NB.         resX[i] := 0
NB.       elseif 0 = ||op(A)|| or 0 = ||Xapprox|| then
NB.         resX[i] := 1 / FP_EPS
NB.       else
NB.         resX[i] := ((||B - Bapprox|| / ||op(A)||) / ||Xapprox||) / FP_EPS
NB.       endif
NB.     endfor
NB.     resX := max(resX[i])
NB.   endif
NB. where
NB.   ||vector|| := norm1t(vector)
NB.   ||matrix|| := ||matrix||_1     when op(A) is A
NB.              := ||matrix||_inf   when op(A) is either A^T or A^H
NB.   Bapprox    := op(A) * Xapprox
NB.
NB. Notes:
NB. - qrt16m models LAPACK's xQRT16

qrt16m=: 1 : 'max@(((FP_EPS * >./@$@(0 {:: [)) , 3 {:: [) (0:`(%@FP_EPS`(%~/@,@,.)@.((*.&(*@{:)) |.))@.(+.&(*@{:)))"1 ] ,.&norm1tc (u 0&{::)~ -~ 1 {:: [)`0:@.((0 e. (, $@(0&{::  )))~ c)'
qrt16v=: 1 : '    (((FP_EPS * >./@$@(0 {:: [)) , 3 {:: [) (0:`(%@FP_EPS`(%~/@,@,.)@.((*.&(*@{:)) |.))@.(+.&(*@{:)))   ] , &norm1t  (u 0&{::)~ -~ 1 {:: [)`0:@.((0 e.    $@(0 {:: [) )   )'

NB. ---------------------------------------------------------
NB. qrt171
NB.
NB. Description:
NB.   Adv. to make dyad to compute the ratio:
NB.     || R' * op(A) || / ( ||A|| * ||B|| * max(m,n,nrhs) * FP_EPS )
NB.   for zero-residual problem
NB. where
NB.   R = B - op(A) * Xapprox
NB.   op(A) is A or A^T or A^H
NB.
NB. Syntax:
NB.   ratio=. (A ; B ; trash ; normA) (calcB qrt171) Xapprox
NB. where
NB.   calcB   - dyad to compute Bapprox; is called as:
NB.               Bapprox=. Xapprox calcB A
NB.   A       - m×n-matrix of full rank
NB.   B       - m×nrhs-matrix or n×nrhs-matrix or m-vector or
NB.             n-vector, exact RHS
NB.   Bapprox - same shape as B, approximate RHS
NB.   Xapprox - n×nrhs-matrix or m×nrhs-matrix or n-vector or
NB.             m-vector, approximate solutions of equation:
NB.               op(A) * X = B
NB.   normA   ≥ 0, the norm of op(A)
NB.   ratio   ≥ 0, the ratio
NB.   m       ≥ 0, the number of rows in A
NB.   n       ≥ 0, the number of columns in A
NB.   nrhs    ≥ 0, the number of RHS and the number of
NB.             columns in B, Bapprox, X and Xapprox
NB.
NB. Formula:
NB.   (m,n) := shape(A)
NB.   nrhs := columns(X)
NB.   if 0 = m or 0 = n or 0 = nrhs then
NB.     err := 0
NB.   else
NB.     C := B - op(A) * Xapprox
NB.     normRS := || C ||_max
NB.     if normRS > (FP_SFMIN / FP_PREC) then
NB.       C := (normRS , 1) scl C
NB.     endif
NB.     err := || C^H * A ||_1
NB.     if || A ||_1 > 0 then
NB.       err := err / || A ||_1
NB.     endif
NB.     if normRS > (FP_SFMIN / FP_PREC) then
NB.       err := err * normRS
NB.     endif
NB.     if || B ||_1 > 0 then
NB.       err := err / || B ||_1
NB.     endif
NB.     err := err / (FP_EPS * max(m,n,nrhs))
NB.   endif
NB.
NB. Notes:
NB. - is called for LS problem, one of:
NB.   - op(A) = A, m ≥ n
NB.   - op(A) = A^T or A^H, m < n
NB. - models LAPACK's xQRT17(1) with the following
NB.   difference:
NB.   - normA is ||op(A)||_1 not ||A||_1
NB. - executing:
NB.     err := err * normRS
NB.   is undoing a scaling, indeed:
NB.     err := (normRS , 1) scl^-1 err

qrt171=: 1 : '(0:`0:`0:`]`(normm@])`((u 0&{::)~ -~ 1 {:: [)`[`((scl~ ,&1)~^:((FP_SFMIN % FP_PREC) < [))`[`0:`(%~`(%~^:(0 < [))/@((2 {. ]) , (*^:((FP_SFMIN % FP_PREC) < [) 2&{)))`((2 {:: ]) , (3 {:: ]) %~^:(0 < [) norm1@((mp~ ct)~ ct^:(</@$)@(0&{::)))`]`]`((2}~ ((((FP_EPS * >./@(, {.)) <@, 1 { ])~ $@(0&{::))~ (c , norm1)@(1&{::)))@[) fork5)`0:@.((0 e. (, $@(0&{::  )))~ c)'

NB. ---------------------------------------------------------
NB. t211
NB.
NB. Description:
NB.   Dyad to compute the relative backward error of
NB.   eigenvectors produced by symmetric eigenvalue problem
NB.   solvers
NB.
NB. Syntax:
NB.   berr=. (H ; normH) t211 (w ; V)
NB. where
NB.   H     - n×n-matrix, the Hermitian (symmetric)
NB.   normH ≥ 0, the norm of H
NB.   w     - n-vector, eigenvalues of H
NB.   V     - n×n-matrix, eigenvectors of H
NB.   berr  ≥ 0, the relative backward error
NB.   n     ≥ 0, the size of H and V
NB.
NB. Formula:
NB.   berr := max(berrH,berrV)
NB. where
NB.   berrH ≥ 0, the relative backward error of H
NB.           factorization as:
NB.             Happrox := V * diagmat(w) * V^H
NB.   berrV ≥ 0, the relative backward error of V unitarity
NB.           (orthogonality) as:
NB.             Iapprox := V * V^H
NB.
NB. Notes:
NB. - models LAPACK's DSYT21(1) and ZHET21(1)

t211=: normi cberrAfact (0&{:: (] mp (* ct)) 1&{::) >. (mp ct) aberrU@(1 {:: ])

NB. ---------------------------------------------------------
NB. t22l
NB. t22r
NB.
NB. Description:
NB.   Dyads to compute the error of eigenvectors produced by
NB.   nonsymmetric eigenvalue problem solvers
NB.
NB. Syntax:
NB.   errL=. (A ; normA ; trash) t22l (w ; L     ; trash)
NB.   errR=. (A ; trash ; normA) t22r (w ; trash ; R    )
NB. where
NB.   A     - n×n-matrix to eigen-decompose
NB.   normA ≥ 0, norm of A
NB.   w     - n-vector, eigenvalues of A
NB.   L     - n×n-matrix, columns with left eigenvectors
NB.   R     - n×n-matrix, columns with right eigenvectors
NB.   errL  ≥ 0, float, scalar, error for L
NB.   errR  ≥ 0, float, scalar, error for R
NB.
NB. Formula:
NB.   err  := max(errL,errR)
NB.   errL := min((||A^H * L - L * W^H|| / max(||L|| , FP_PREC)) / max(||A|| , FP_SFMIN) , 1) / FP_PREC
NB.   errR := min((||A   * R - R * W  || / max(||R|| , FP_PREC)) / max(||A|| , FP_SFMIN) , 1) / FP_PREC
NB. where
NB.   ||matrix|| := ||matrix||_inf  for errL
NB.              := ||matrix||_1    for errR
NB.   W          := diagmat(w)
NB.
NB. Notes:
NB. - models LAPACK's xGET22

t22l=: FP_PREC %~ 1 <. ((((mp~ ct@(0&{::))~ 1&{::) - (1 {:: ]) *"1 +@(0 {:: ])) (% FP_PREC&>.)&normi 1 {:: ]) % FP_SFMIN >. 1 {:: [
t22r=: FP_PREC %~ 1 <. ((((mp~     0&{:: )~ 2&{::) - (2 {:: ]) *"1   (0 {:: ])) (% FP_PREC&>.)&norm1 2 {:: ]) % FP_SFMIN >. 2 {:: [

NB. ---------------------------------------------------------
NB. t511l
NB. t511u
NB.
NB. Description:
NB.   Dyads to compute the error of generalized Schur form
NB.   produced by hgexxsxx
NB.
NB. Syntax:
NB.   errl=. (Al ; normA) t511l (Bl ; Ql ; Zl)
NB.   erru=. (Au ; normA) t511u (Bu ; Qu ; Zu)
NB. where
NB.   Al    - n×n-matrix from Hessenberg-triangular pair to
NB.           reduce to generalized Schur form by hgezqsxx:
NB.             Ql^H * Bl * Zl = Al
NB.   Au    - n×n-matrix from Hessenberg-triangular pair to
NB.           reduce to generalized Schur form by hgeqzsxx:
NB.             Qu * Bu * Zu^H = Au
NB.   normA ≥ 0, the norm of A
NB.   Bl    - n×n-matrix from generalized Schur form's matrix
NB.           pair (S,P) produced by hgezqsxx
NB.   Bu    - n×n-matrix from generalized Schur form's matrix
NB.           pair (S,P) produced by hgeqzsxx
NB.   Ql    - n×n-matrix, the unitary (orthogonal) matrix of
NB.           left Schur vectors produced by hgezqsxx
NB.   Qu    - n×n-matrix, the unitary (orthogonal) matrix of
NB.           left Schur vectors produced by hgeqzsxx
NB.   Zl    - n×n-matrix, the unitary (orthogonal) matrix of
NB.           right Schur vectors produced by hgezqsxx
NB.   Zu    - n×n-matrix, the unitary (orthogonal) matrix of
NB.           right Schur vectors produced by hgeqzsxx
NB.   errl  ≥ 0, float, scalar, error of generalized Schur
NB.           form produced by hgezqsxx
NB.   erru  ≥ 0, float, scalar, error of generalized Schur
NB.           form produced by hgeqzsxx
NB.
NB. Formula:
NB.   errX := max(errH,errT)
NB. where
NB.   S,P,Q,Z - output from hgexxsvv
NB.   errX is calculated by cberrAfact
NB.   - hgezqxxx:
NB.       Happrox := Q^H * S * Z
NB.       Tapprox := Q^H * P * Z
NB.   - hgeqzxxx:
NB.       Happrox := Q * S * Z^H
NB.       Tapprox := Q * P * Z^H
NB.
NB. Notes:
NB. - models LAPACK's xGET51(1)

t511l=: normi cberrAfact (1&{:: (mp~ ct)~ 0&{::  mp     2&{::)
t511u=: norm1 cberrAfact (1&{::  mp       0&{:: (mp ct) 2&{::)

NB. ---------------------------------------------------------
NB. t513l
NB. t513u
NB.
NB. Description:
NB.   Monads to compute the error of Schur vectors produced
NB.   by hgexxsxx
NB.
NB. Syntax:
NB.   errl=. t513l Vl
NB.   erru=. t513u Vu
NB. where
NB.   Vl   - n×n-matrix, the unitary (orthogonal) matrix of
NB.          Schur vectors produced by hgezqsxx
NB.   errl ≥ 0, float, scalar, error of Schur vectors
NB.          produced by hgezqsxx
NB.   Vu   - n×n-matrix, the unitary (orthogonal) matrix of
NB.          Schur vectors produced by hgeqzsxx
NB.   erru ≥ 0, float, scalar, error of left Schur vectors
NB.          produced by hgeqzsxx
NB.
NB. Formula:
NB.   n := size(Vl)
NB.   errX := max(errQ,errZ)
NB. where
NB.   S,P,Q,Z - output from hgexxsvv
NB.   errX is calculated by aberrU
NB.   - hgezqxxx:
NB.       errQ := ||Q^H * Q - I|| / (FP_PREC * n)
NB.       errZ := ||Z^H * Z - I|| / (FP_PREC * n)
NB.   - hgeqzxxx:
NB.       errQ := ||Q * Q^H - I|| / (FP_PREC * n)
NB.       errZ := ||Z * Z^H - I|| / (FP_PREC * n)
NB.
NB. Notes:
NB. - models LAPACK's xGET51(3)

t513l=: (mp~ ct) aberrU
t513u=: (mp  ct) aberrU

NB. ---------------------------------------------------------
NB. t52ll
NB. t52lr
NB. t52lb
NB. t52ul
NB. t52ur
NB. t52ub
NB.
NB. Description:
NB.   Dyads to compute the error of Schur vectors produced
NB.   by tgevcxxx
NB.
NB. Syntax:
NB.   errl=. SPl     t52ll Yapprox       NB. to test tgevcll
NB.   errl=. SPQZHTl t52ll YQapprox      NB. to test tgevcllb
NB.   errl=. SPl     t52lr Xapprox       NB. to test tgevclr
NB.   errl=. SPQZHTl t52lr XZapprox      NB. to test tgevclrb
NB.   errl=. SPl     t52lb YXapprox      NB. to test tgevclb
NB.   errl=. SPQZHTl t52lb YQXZapprox    NB. to test tgevclbb
NB.   erru=. SPu     t52ul Yapprox       NB. to test tgevcul
NB.   erru=. SPQZHTu t52ul QYapprox      NB. to test tgevculb
NB.   erru=. SPu     t52ur Xapprox       NB. to test tgevcur
NB.   erru=. SPQZHTu t52ur ZXapprox      NB. to test tgevcurb
NB.   erru=. SPu     t52ub YXapprox      NB. to test tgevcub
NB.   erru=. SPQZHTu t52ub QYZXapprox    NB. to test tgevcubb
NB. where
NB.   SPl        - 2×n×n-matrix (S,:P), generalized Schur
NB.                form, produced by hgezqsxx
NB.   SPu        - 2×n×n-matrix (S,:P), generalized Schur
NB.                form, produced by hgeqzsxx
NB.   SPQZHTl    - 6×n×n-matrix (S,P,Q,Z,H,:T), composed from
NB.                output of hgezqsvv and gghrdlnn
NB.   SPQZHTu    - 6×n×n-matrix (S,P,Q,Z,H,:T), composed from
NB.                output of hgeqzsvv and gghrdunn
NB.   Yapprox    - n×n-matrix, left eigenvectors produced by
NB.                tgevcxl
NB.   Xapprox    - n×n-matrix, right eigenvectors produced by
NB.                tgevcxr
NB.   YXapprox   - 2×n×n-matrix, left and right eigenvectors
NB.                produced by tgevcxb,
NB.                  YXapprox -: Yapprox ,: Xapprox
NB.   YQapprox   - n×n-matrix, left eigenvectors Y*Q produced
NB.                by tgevcllb
NB.   XZapprox   - n×n-matrix, right eigenvectors X*Z
NB.                produced by tgevclrb
NB.   QYapprox   - n×n-matrix, left eigenvectors Q*Y produced
NB.                by tgevculb
NB.   ZXapprox   - n×n-matrix, right eigenvectors Z*X
NB.                produced by tgevcurb
NB.   YQXZapprox - 2×n×n-matrix, left and right eigenvectors
NB.                produced by tgevclbb,
NB.                  YQXZapprox -: YQapprox ,: XZapprox
NB.   QYZXapprox - 2×n×n-matrix, left and right eigenvectors
NB.                produced by tgevcubb,
NB.                  QYZXapprox -: QYapprox ,: ZXapprox
NB.   errl       ≥ 0, float, scalar, error of eigenvectors
NB.                produced by tgevclxx
NB.   erru       ≥ 0, float, scalar, error of eigenvectors
NB.                produced by tgevcuxx
NB.
NB. Formula:
NB.   berr := max(berr0,berr1)
NB. where
NB.   ||M|| := max(||M||_1 , FP_SFMIN)
NB.   ||v|| := max(|Re(v(i))|+|Im(v(i))|)
NB.   e1(i) - i-th eigenvalue, also i-th element on S
NB.           diagonal
NB.   e2(i) - i-th eigenvalue, also i-th element on P
NB.           diagonal
NB.   l(i)  - i-th left eigenvector
NB.   lb(i) - i-th back transformed left eigenvector
NB.   r(i)  - i-th right eigenvector
NB.   rb(i) - i-th back transformed right eigenvector
NB.   - tgevcll:
NB.       berr0 := max(||l(i) * (e2(i)*S - e1(i)*P)  || / (FP_PREC * max(|| e2(i)*S   ||,|| e1(i)*P   ||)))
NB.       berr1 := max(| ||l(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevclr:
NB.       berr0 := max(||r(i) * (e2(i)*S - e1(i)*P)^H|| / (FP_PREC * max(||(e2(i)*S)^H||,||(e1(i)*P)^H||)))
NB.       berr1 := max(| ||r(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevclb:
NB.       berr0 := berr(tgevcll)
NB.       berr1 := berr(tgevclr)
NB.   - tgevcllb:
NB.       berr0 := max(||l(i) * (e2(i)*H - e1(i)*T)  || / (FP_PREC * max(|| e2(i)*H   ||,|| e1(i)*T   ||)))
NB.       berr1 := max(| ||lb(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevclrb:
NB.       berr0 := max(||r(i) * (e2(i)*H - e1(i)*T)^H|| / (FP_PREC * max(||(e2(i)*H)^H||,||(e1(i)*T)^H||)))
NB.       berr1 := max(| ||rb(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevclbb:
NB.       berr0 := berr(tgevcllb)
NB.       berr1 := berr(tgevclrb)
NB.   - tgevcul:
NB.       berr0 := max(||(e2(i)*S - e1(i)*P)^H * l(i)|| / (FP_PREC * max(||(e2(i)*S)^H||,||(e1(i)*P)^H||)))
NB.       berr1 := max(| ||l(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevcur:
NB.       berr0 := max(||(e2(i)*S - e1(i)*P)   * r(i)|| / (FP_PREC * max(|| e2(i)*S   ||,|| e1(i)*P   ||)))
NB.       berr1 := max(| ||r(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevcub:
NB.       berr0 := berr(tgevcul)
NB.       berr1 := berr(tgevcur)
NB.   - tgevculb:
NB.       berr0 := max(||(e2(i)*H - e1(i)*T)^H * l(i)|| / (FP_PREC * max(||(e2(i)*H)^H||,||(e1(i)*T)^H||)))
NB.       berr1 := max(| ||lb(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevcurb:
NB.       berr0 := max(||(e2(i)*H - e1(i)*T)   * r(i)|| / (FP_PREC * max(|| e2(i)*H   ||,|| e1(i)*T   ||)))
NB.       berr1 := max(| ||rb(i)|| - 1 |) / (FP_PREC * n)
NB.   - tgevcubb:
NB.       berr0 := berr(tgevculb)
NB.       berr1 := berr(tgevcurb)
NB.
NB. Notes:
NB. - models LAPACK's xGET52
NB. - t52xx are non-iterative and are require O(N^3) RAM

t52ll=: (     normir@:((((  norm1r@:((mp"1 2 -/"3)~      )                                            )       %  (FP_PREC * FP_SFMIN >.      >./"1 @:( norm1       "2)@[))~ (_2&{. *"_ 1|:@|.@:(diag"2)@(2&{.)))~)) >. (     normir@:<:@:normitr % FP_PREC * c)@]
t52lr=: (     normir@:((((                                    norm1r@:((mp"2 1~ -/"3)~ +    )         )       %  (FP_PREC * FP_SFMIN >.      >./"1 @:(       normi "2)@[))~ (_2&{. *"_ 1|:@|.@:(diag"2)@(2&{.)))~)) >. (     normir@:<:@:normitr % FP_PREC * c)@]
t52lb=: (>./@:normir@:((((((norm1r@:( mp"1 2      ~    {.) ,: norm1r@:((mp"2 1      )  + @{:))~ -/"3)~) (>./@:%) (FP_PREC * FP_SFMIN >. |:@:(>./"2)@:((norm1,normi)"2)@[))~ (_2&{. *"_ 1|:@|.@:(diag"2)@(2&{.)))~)) >. (>./@:normir@:<:@:normitr % FP_PREC * c)@]

t52ul=: (     normir@:((((  norm1r@:((mp"1 2 -/"3)~ ct   )                                            )       %  (FP_PREC * FP_SFMIN >.      >./"1 @:( normi       "2)@[))~ (_2&{. *"_ 1|:@|.@:(diag"2)@(2&{.)))~)) >. (     normir@:<:@:normitc % FP_PREC * c)@]
t52ur=: (     normir@:((((                                    norm1r@:((mp"2 1~ -/"3)~ |:   )         )       %  (FP_PREC * FP_SFMIN >.      >./"1 @:(       norm1 "2)@[))~ (_2&{. *"_ 1|:@|.@:(diag"2)@(2&{.)))~)) >. (     normir@:<:@:normitc % FP_PREC * c)@]
t52ub=: (>./@:normir@:((((((norm1r@:( mp"1 2      ~ ct@{.) ,: norm1r@:((mp"2 1      )  |:@{:))~ -/"3)~) (>./@:%) (FP_PREC * FP_SFMIN >. |:@:(>./"2)@:((normi,norm1)"2)@[))~ (_2&{. *"_ 1|:@|.@:(diag"2)@(2&{.)))~)) >. (>./@:normir@:<:@:normitc % FP_PREC * c)@]
