NB. Eigenvalues and Schur form
NB.
NB. hgexxe     Eigenvalues of pair of structured matrices
NB. hgexxs     Eigenvalues and the Schur form of pair of
NB.            structured matrices
NB.
NB. testhgeqe  Test hgexxe by general matrices given
NB. testhgeqs  Test hgexxs by general matrices given
NB. testeqz    Adv. to make verb to test hgexxx by matrices
NB.            of generator and shape given
NB.
NB. Version: 0.6.8 2010-11-30
NB.
NB. Copyright 2010 Igor Zhuravlov
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
NB. Name: laqr1
NB. Description: 1st column of (H-s1*I)*(H-s2*I)
NB. Syntax: vK=. (s1,s2) laqr1 H
NB. where   H - 2×2-matrix or 3×3-matrix
NB. TODO: tacit
NB. Notes: implements LAPACK's xLAQR1

laqr1=: 4 : 0
  's1 s2'=. x
  (((-&s1) upddiag ]) (mp (% norm1t)) (((-&s2) updl 0) {."1 ])) y
)

NB. ---------------------------------------------------------
NB. hgezqeo
NB. hgeqzeo
NB.
NB. Description:
NB.   Calculate generalized eigenvalues of hs-segment
NB.
NB. Syntax:
NB.   'HTupd signbc'=. hs hgexxeo H ,: T
NB. where
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines eigenvalues range
NB.   H     - n×n-matrix, either lower (hgezqeo) or upper
NB.           (hgeqzeo) Hessenberg inside the submatrix
NB.           H[h:h+s-1,h:h+s-1], and lower (hgezqeo) or
NB.           upper (hgeqzeo) triangular outside
NB.   T     - n×n-matrix, either lower (hgezqeo) or upper
NB.           (hgeqzeo) triangular
NB.   HTupd -: Hupd ,: Tupd
NB.   Hupd  - n×n-matrix, being H with hs-segment of diagonal
NB.           replaced by alpha (see hgexx)
NB.   Tupd  - n×n-matrix, being T with hs-segment of diagonal
NB.           replaced by beta (see hgexx)
NB.   signbc - s-vector, scaling factors to form Q,Z later

hgezqeo=: 4 : 0
  'Hd Td'=. (0 , x) diag"2 y
  absb=. | Td
  'signbc Td'=. (,:~ FP_SFMIN < absb) } 1 0 2 |: (1 ,: * Td) ,: (0 ,: absb)
  ((((Hd * signbc) ,: Td) (;"1) 0 , x) setdiag"1 2 y) ; signbc
)

hgeqzeo=: 4 : 0
  'Hd Td'=. (0 , x) diag"2 y
  absb=. | Td
  'signbc Td'=. (,:~ FP_SFMIN < absb) } 1 0 2 |: (1 ,: + * Td) ,: (0 ,: absb)
  ((((Hd * signbc) ,: Td) (;"1) 0 , x) setdiag"1 2 y) ; signbc
)

NB. ---------------------------------------------------------
NB. hgezqso
NB. hgeqzso
NB.
NB. Description:
NB.   Calculate generalized eigenvalues of hs-segment and
NB.   reduce corresponding rows (hgezqso) or columns
NB.   (hgeqzso) to generalized Schur form
NB.
NB. Syntax:
NB.   'HTupd signbc'=. hs hgexxso H ,: T
NB. where
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines eigenvalues range
NB.   H     - n×n-matrix, either lower (hgezqso) or upper
NB.           (hgeqzso) Hessenberg inside the submatrix
NB.           H[h:h+s-1,h:h+s-1], and lower (hgezqso) or
NB.           upper (hgeqzso) triangular outside
NB.   T     - n×n-matrix, either lower (hgezqeo) or upper
NB.           (hgeqzso) triangular
NB.   HTupd -: Hupd ,: Tupd
NB.   Hupd  - n×n-matrix, being H with rows (hgezqso) or
NB.           columns (hgeqzso) from hs-segment transformed
NB.           to Shur form
NB.   Tupd  - n×n-matrix, being T with rows (hgezqso) or
NB.           columns (hgeqzso) from hs-segment transformed
NB.           to Shur form
NB.   signbc - s-vector, scaling factors to form Q,Z later

hgezqso=: 4 : 0
  lios=. dhs2lios x
  'y signbc'=. x hgezqeo y
  subHT=. lios {"2 y
  (((,:~ lios >/ (i. c y)) } subHT ,: subHT *"2 signbc) lios }"2 y) ; signbc
)

hgeqzso=: 4 : 0
  lios=. dhs2lios x
  'y signbc'=. x hgeqzeo y
  subHT=. lios {"1 y
  (((,:~ (i. c y) </ lios) } subHT ,: subHT *"1 signbc) lios }"1 y) ; signbc
)

NB. ---------------------------------------------------------
NB. hgezq
NB. hgeqz
NB.
NB. Description:
NB.   Adv. to make verbs to find eigenvalues of either lower
NB.   (hgezq) or upper (hgeqz) Hessenberg-triangular pair
NB.   (H,T) and, optionally, to reduce this pair to
NB.   generalized Schur form
NB.
NB. Syntax:
NB.   vapp=. hgexxxo`init`reset`step hgexx
NB. where
NB.   hgexxxo - monad to calculate generalized eigenvalues
NB.             of hs-segment and, optionally, to reduce
NB.             corresponding columns to generalized Schur
NB.             form, is either hgezqeo or hgezqso for
NB.             hgezq, or hgeqzeo or hgeqzso for hgeqz, is
NB.             called as:
NB.               'HTupd signbc'=. hgexxxo hs ; HT
NB.   init    - dyad to initialize counters, is called as:
NB.               'ifrstm ilastm'=. (h , ilast) init (0 , n-1)
NB.   reset   - dyad to reset counters optionally, is called
NB.             as:
NB.               'ifrstm ilastm'=. (h , ilast) reset (ifrstm , ilastm)
NB.   step    - monad to change counter optionally, is
NB.             called as:
NB.               ifrstm=. ifirst step ifrstm
NB.   vapp    - dyad to find eigenvalues of lower (upper)
NB.             Hessenberg-triangular pair (H,T) and,
NB.             optionally, to reduce this pair to
NB.             generalized Schur form, is called as:
NB.               'HTupd dQ dZ'=. hs vapp HT
NB.             see hgeqzxxxx
NB.
NB. Notes:
NB. - non-converged eigenvalues will be set to NaN

hgezq=: 1 : 0
:
  '`hgezqxo init reset step'=. m
  e=. +/ 'h s'=. x
  dQ=. dZ=. 0 4 $ 0
  abnorm=. (0 2 ,. ,.~ x) norms"2 ;. 0 y
smoutput ": 'abnorm' ; abnorm
  'atol btol'=. abtol=. FP_SFMIN >. FP_PREC * abnorm
  'ascale bscale'=. abscale=. % FP_SFMIN >. abnorm
  'y signbc'=. ((c y) (] , -) e) hgezqxo y  NB. process eigenvalues (columns) h+s:n-1
  dZ=. dZ , 4 {."1 signbc ,. (c y) th2lios e

  NB. Eigenvalues h+s:n-1 have been found.

  NB. Initialize dynamic indices
  ilast=. <: e
  'ifrstm ilastm'=. (h , ilast) init (0 , <: c y)
                  NB. ifrstm - the column of the last
                  NB.          splitting column to the left
                  NB.          of the column ilast, this is
                  NB.          always at least h
  iiter=. 0       NB. counts iterations since the last
                  NB. eigenvalue was found, to tell when to
                  NB. use an extraordinary shift
  eshift=. 0
  maxit=. 30 * s  NB. the maximum number of ZQ sweep allowed
  jiter=. 0

  NB. Main ZQ iteration loop
  NB. Row operations modify columns ifrstm:*
  NB. Column operations modify rows *:ilastm
  while. jiter < maxit do.
smoutput 'loop 170 JITER = ',":jiter
    goto60=. 1  NB. set default branching
    NB. split the matrix if possible, by two tests:
    NB.   1. H[j-1,j]=0 OR j=h
    NB.   2. T[j,j]=0
    if. ilast ~: h do.
      if. atol < sorim (< 0 , ilast - 1 0) { y do.
        if. btol < | (< 1 , ,~ ilast) { y do.
          NB. general case: j < ilast
          j=. <: ilast
          while. j >: h do.
smoutput 'loop 40 J = ',":j
            NB. test 1: H[j-1,j]=0 OR j=h
            if. j = h do.
              ilazro=. 1
smoutput 'ILAZRO = .T.'
            elseif. atol >: sorim (< 0 , j - 1 0) { y do.
              y=. 0 (< 0 , j - 1 0) } y
              ilazro=. 1
smoutput 'ILAZRO = .T.'
            elseif. do.
              ilazro=. 0
smoutput 'ILAZRO = .F.'
            end.
            NB. test 2: T[j,j]=0
            if. btol > | (< 1 , ,~ j) { y do.
smoutput 'Test 2 = .T.'
              y=. 0 (< 1 , ,~ j) } y
              NB. test 2a: check for 2 consecutive small
              NB, superdiagonals in H
              ilazr2=. 0
              if. -. ilazro do.
                'Hj1j Hjj1 Hjj'=. sorim ((<"1) 0 ,. (_1 0,0 1,:0 0) + j) { y
                if. 0 >: (Hj1j , Hjj) mp (ascale * (Hjj1 , -atol)) do.
                  ilazr2=. 1
smoutput 'ILAZR2 = .T.'
                end.
              end.
              NB. If both tests (1 & 2) pass, i.e., the
              NB. leading diagonal element of T in the block
              NB. is zero, then split a 1x1 block off at the
              NB. left, i.e. at the j-th row/column. The
              NB. leading diagonal element of the remainder
              NB. can also be zero, so this may have to be
              NB. done repeatedly.
              if. ilazro +. ilazr2 do.
smoutput 'ILAZRO OR ILAZR2 = .T.'
                jch=. j
                lios=. (>: ilastm) th2lios j
                while. jch < ilast do.
smoutput 'loop 20 JCH = ',":jch
smoutput ": 'before rotation subHT' ; (< a: ; lios ; (jch + 0 1)) { y
                  'y cs'=. rot rotga y ; (< 0 ; lios ; (jch + 0 1)) ; 0
                  lios=. }. lios
                  y=. (< 1 ; lios ; (jch + 0 1)) (cs & rot) upd y
smoutput ": 'after rotation subHT' ; (< a: ; lios ; (jch + 0 1)) { y
                  dZ=. dZ , (+ cs) , jch + 0 1
                  if. ilazr2 do.
                    y=. (< 0 , jch - 1 0) (* & ({. cs)) upd y
                    ilazr2=. 0
                  end.
                  if. btol <: sorim (< 1 , jch + 1 1) { y do.
smoutput ": 'after loop 20 HT' ; y
                    if. ilast > >: jch do.
                      ifirst=. >: jch
                      goto60=. 0
smoutput '(1) GOTO 70'
                    else.
smoutput '(3) GOTO 60'
                    end.
                    goto_l60.
                  end.
                  y=. 0 (< 1 , jch + 1 1) } y
                  jch=. >: jch
                end.
smoutput ": 'after loop 20 HT' ; y
smoutput '(2) GOTO 50'
              else.
smoutput 'ILAZRO OR ILAZR2 = .F.'
                NB. Only test 2 passed - chase the zero to
                NB. T[ilast,ilast], then process as in the
                NB. case T[ilast,ilast]=0
                jch=. j
                liosr=. (>: ilastm) th2lios <: j
                liosc=. (2 + j) th2lios ifrstm
                while. jch < ilast do.
smoutput 'loop 30 JCH = ',":jch
smoutput ": 'before rotation subHT' ; (< a: ; liosr ; (jch + 0 1)) { y
                  'y cs'=. rot rotga y ; (< 1 ; (2 }. liosr) ; (jch + 0 1)) ; 0
                  y=. (< 0 ; liosr ; (jch + 0 1)) (cs & rot) upd y
smoutput ": 'after rotation subHT' ; (< a: ; liosr ; (jch + 0 1)) { y
                  dZ=. dZ , (+ cs) , jch + 0 1
smoutput ": 'before rotation subHT' ; (< a: ; (jch - 0 1) ; liosc) { y
                  'y cs'=. (rot &. |:) rotga y ; (< 0 ; (jch - 0 1) ; liosc) ; < < a: ; _1
                  y=. (< 1 ; (jch - 0 1) ; (_2 }. liosc)) (cs & (rot &. |:)) upd y
smoutput ": 'after rotation subHT' ; (< a: ; (jch - 0 1) ; liosc) { y
                  dQ=. dQ , cs , jch - 0 1
                  liosr=. }. liosr
                  liosc=. liosc , 2 + jch
                  jch=. >: jch
                end.
smoutput ": 'after loop 30 HT' ; y
              end.
smoutput '(3) GOTO 50'
              goto_l50.
            elseif. ilazro do.
smoutput 'ILAZRO = .T.'
              ifirst=. j
              goto60=. 0
smoutput '(2) GOTO 70'
              goto_l60.
            end.
            j=. <: j
          end.
          NB. drop-through is impossible
          ((< _.) setdiag"2 y) ; ,~ a:  NB. set all eigenvalues to NaN
          return.
        else.
          y=. 0 (< 1 , ,~ ilast) } y
smoutput '(1) GOTO 50'
        end.
        label_l50.
smoutput 'label 50'
        NB. T[ilast,ilast]=0 - clear H[ilast-1,ilast] to
        NB. split off a 1x1 block
        lios=. (>: ilast) th2lios ifrstm
smoutput ": 'before rotation subHT' ; (< a: ; (ilast - 0 1) ; lios) { y
        'y cs'=. (rot &. |:) rotga y ; (< 0 ; (ilast - 0 1) ; lios) ; < a: ; _1
        y=. (< 1 ; (ilast - 0 1) ; (}: lios)) (cs & (rot &. |:)) upd y
smoutput ": 'after rotation subHT' ; (< a: ; (ilast - 0 1) ; lios) { y
        dQ=. dQ , cs , ilast - 0 1
      else.
        y=. 0 (< 0 , ilast - 1 0) } y
smoutput '(2) GOTO 60'
      end.
    else.
smoutput '(1) GOTO 60'
    end.
    label_l60.
    if. goto60 do.
smoutput 'label 60'
      NB. H[ilast-1,ilast]=0 - standartize B, set alpha and
      NB. beta
      'y signbc'=. (ilast , 1) hgezqxo y  NB. process ilast-th eigenvalue (column)
      dQ=. dQ , 4 {."1 signbc ,. ilast
      NB. goto next block - exit if finished
      ilast=. <: ilast
      if. ilast < h do.
        NB. normal exit
        'y signbc'=. (0 , 0 >. <: h) hgezqxo y  NB. process eigenvalues (columns) 0:h-1
        dQ=. dQ , 4 {."1 signbc ,. i. 0 >. <: h
        y ; dQ ; dZ
        return.
      end.
      NB. reset counters
      iiter=. 0
      eshift=. 0
      'ifrstm ilastm'=. (h , ilast) reset (ifrstm , ilastm)
    else.
smoutput 'label 70'
      NB. ZQ step
      NB. This iteration only involves rows/columns
      NB. ifirst:ilast. We assume ifirst<ilast, and that the
      NB. diagonal of B is non-zero
      iiter=. >: iiter
      ifrstm=. ifirst step ifrstm
      NB. Compute the shift.
      NB. At this point, ifirst<ilast, and the diagonal
      NB. elements of T[ifirst:ilast,ifirst:ilast] are larger
      NB. than btol in magnitude
      if. 10 | iiter do.
smoutput 'Wilkinson shift'
        NB. The Wilkinson shift, i.e., the eigenvalues of the
        NB. bottom-right 2x2 block of T^_1*H which is nearest
        NB. to the bottom-right element.
        NB. We factor T as D*L, where L is unit lower
        NB. triangular, and compute L^_1*(D^_1*H)
smoutput ": 'abscale' ; abscale
smoutput ": 'HT[ilast-1,ilast-1]' ; (< a: ; ;~ ilast - 1) { y
        'L21 DA11 DA12 DA21 DA22'=. %/ (6 0 1 2 3 ,: 7 4 4 7 7) ({,) abscale * ((< a: ; ;~ ilast - 1 0) { y)
        IBA22=. DA22 - L21 * DA12
smoutput ": 'L21' ; L21 ; 'DA11' ; DA11 ; 'DA12' ; DA12 ; 'DA21' ; DA21 ; 'DA22' ; DA22 ; 'IBA22' ; IBA22
        t1=. -: DA11 + IBA22
        rtdisc=. %: (t1 , DA21 , -DA11) mp (t1 , DA12 , DA22)
        temp=. +/ (*/) +. rtdisc , t1 - IBA22
        shift=. t1 - temp condneg rtdisc
smoutput ": 't1' ; t1 ; 'rtdisc' ; rtdisc ; 'temp' ; temp ; 'shift' ; shift
      else.
        NB. Exceptional shift. Chosen for no paticularly good
        NB. reason
        eshift=. eshift + + %/ abscale * (;/ 0 1 ,. (0 _1 ,: _1 _1) + ilast) { y
smoutput ": 'eshift' ; eshift
        shift=. eshift
      end.
      NB. now check for two consecutive small subdiagonals
      HTd=. (0 1 ,"0 1 ilast (] , -) ifirst) diag"1 2/ y
      ctemp=. (- (shift&*))/ abscale * {. HTd
      temp=. (sorim }."1 ctemp) ,: ascale * sorim (< 1 ; 0 ; <<0) { HTd
      tempr=. >./ temp
smoutput ": 'HTd' ; HTd ; 'ctemp' ; ctemp ; 'temp' ; temp ; 'tempr' ; tempr
      temp=. temp %"1 ((0 , 1 - FP_EPS) I. tempr) } 1 , tempr ,: 1
      'istart ctemp'=. (+&ifirst , {&ctemp) (ilast - ifirst) | >: (>:/ temp * atol ,: sorim (< 1 ; 0 ; <<_1) { HTd) i: 1
smoutput ": 'temp' ; temp ; 'istart' ; istart ; 'ctemp' ; ctemp
smoutput 'label 90'
      NB. do an implicit-shift ZQ sweep
      NB. initial Z
smoutput ": 'fg' ; ctemp , ascale * (< 0 , istart + 0 1) { y
      cs=. }: lartg ctemp , ascale * (< 0 , istart + 0 1) { y
      NB. sweep
      j=. istart
      liosr=. (>: ilastm) th2lios j
      liosc=. (j + 2) th2lios ifrstm
      while. j < ilast do.
smoutput 'loop 150 J = ',":j
        lios=. j + 0 1
        NB. is a first iteration?
smoutput ": 'before rotation subHT' ; (< a: ; liosr ; lios) { y
        if. j = istart do.
          y=. (< a: ; liosr ; lios) (cs & rot"2) upd y
        else.
smoutput ": 'fg' ; (< 0 ; ({. liosr) ; lios) { y
          'y cs'=. rot rotga y ; (< 0 ; liosr ; lios) ; 0
          liosr=. }. liosr
          y=. (< 1 ; liosr ; lios) (cs & rot) upd y
        end.
smoutput ": 'cs' ; cs
smoutput ": 'after rotation subHT' ; (< a: ; liosr ; lios) { y
        dZ=. dZ , (+ cs) , lios
        lios=. j + 1 0
smoutput ": 'before rotation subHT' ; (< a: ; lios ; liosc) { y
        'y cs'=. (rot &. |:) rotga y ; (< 1 ; lios ; liosc) ; < < a: ; _1
        NB. isn't a last iteration?
        if. j < <: ilast do.
          liosc=. liosc , j + 2
        end.
        y=. (< 0 ; lios ; liosc) (cs & (rot &. |:)) upd y
smoutput ": 'after rotation subHT' ; (< a: ; lios ; liosc) { y
        dQ=. dQ , cs , lios
        j=. >: j
      end.
smoutput ": 'after loop 150 HT' ; y
    end.
    jiter=. >: jiter
  end.
  NB. drop-through means non-convergence, set incorrect eigenvalues 0:ilast to NaN
  ((_. ; 0 0 , >: ilast) setdiag"2 y) ; ,~ a:
)

hgeqz=: 1 : 0
:
  '`hgeqzxo init reset step'=. m
  e=. +/ 'h s'=. x
  dQ=. dZ=. 0 4 $ 0
  abnorm=. (0 2 ,. ,.~ x) norms"2 ;. 0 y
smoutput ": 'abnorm' ; abnorm
  'atol btol'=. abtol=. FP_SFMIN >. FP_PREC * abnorm
  'ascale bscale'=. abscale=. % FP_SFMIN >. abnorm
  'y signbc'=. ((c y) (] , -) e) hgeqzxo y  NB. process eigenvalues (columns) h+s:n-1
  dZ=. dZ , 4 {."1 signbc ,. (c y) th2lios e

  NB. Eigenvalues h+s:n-1 have been found.

  NB. Initialize dynamic indices
  ilast=. <: e
  'ifrstm ilastm'=. (h , ilast) init (0 , <: c y)
                  NB. ifrstm - the row of the last splitting
                  NB.          row above row ilast, this is
                  NB.          always at least h
  iiter=. 0       NB. counts iterations since the last
                  NB. eigenvalue was found, to tell when to
                  NB. use an extraordinary shift
  eshift=. 0
  maxit=. 30 * s  NB. the maximum number of QZ sweep allowed
  jiter=. 0

  NB. Main QZ iteration loop
  NB. Column operations modify rows ifrstm:*
  NB. Row operations modify columns *:ilastm
  while. jiter < maxit do.
smoutput 'loop 170 JITER = ',":jiter
    goto60=. 1  NB. set default branching
    NB. split the matrix if possible, by to tests:
    NB.   1. H[j,j-1]=0 OR j=h
    NB.   2. T[j,j]=0
    if. ilast ~: h do.
      if. atol < sorim (< 0 , ilast - 0 1) { y do.
        if. btol < | (< 1 , ,~ ilast) { y do.
          NB. general case: j < ilast
          j=. <: ilast
          while. j >: h do.
smoutput 'loop 40 J = ',":j
            NB. test 1: H[j,j-1]=0 OR j=h
            if. j = h do.
              ilazro=. 1
smoutput 'ILAZRO = .T.'
            elseif. atol >: sorim (< 0 , j - 0 1) { y do.
              y=. 0 (< 0 , j - 0 1) } y
              ilazro=. 1
smoutput 'ILAZRO = .T.'
            elseif. do.
              ilazro=. 0
smoutput 'ILAZRO = .F.'
            end.
            NB. test 2: T[j,j]=0
            if. btol > | (< 1 , ,~ j) { y do.
smoutput 'Test 2 = .T.'
              y=. 0 (< 1 , ,~ j) } y
              NB. test 2a: check for 2 consecutive small
              NB. subdiagonals in H
              ilazr2=. 0
              if. -. ilazro do.
                'Hjj1 Hj1j Hjj'=. sorim ((<"1) 0 ,. (0 _1,1 0,:0 0) + j) { y
                if. 0 >: (Hjj1 , Hjj) mp (ascale * (Hj1j , -atol)) do.
                  ilazr2=. 1
smoutput 'ILAZR2 = .T.'
                end.
              end.
              NB. If both tests (1 & 2) pass, i.e., the
              NB. leading diagonal element of T in the block
              NB. is zero, then split a 1x1 block off at the
              NB. top, i.e. at the j-th row/column. The
              NB. leading diagonal element of the remainder
              NB. can also be zero, so this may have to be
              NB. done repeatedly.
              if. ilazro +. ilazr2 do.
smoutput 'ILAZRO OR ILAZR2 = .T.'
                jch=. j
                lios=. (>: ilastm) th2lios j
                while. jch < ilast do.
smoutput 'loop 20 JCH = ',":jch
smoutput ": 'before rotation subHT' ; (< a: ; (jch + 0 1) ; lios) { y
                  'y cs'=. (rot &. |:) rotga y ; (< 0 ; (jch + 0 1) ; lios) ; < < a: ; 0
                  lios=. }. lios
                  y=. (< 1 ; (jch + 0 1) ; lios) (cs & (rot &. |:)) upd y
smoutput ": 'after rotation subHT' ; (< a: ; (jch + 0 1) ; lios) { y
                  dQ=. dQ , (+ cs) , jch + 0 1
                  if. ilazr2 do.
                    y=. (< 0 , jch - 0 1) (* & ({. cs)) upd y
                    ilazr2=. 0
                  end.
                  if. btol <: sorim (< 1 , jch + 1 1) { y do.
smoutput ": 'after loop 20 HT' ; y
                    if. ilast > >: jch do.
                      ifirst=. >: jch
                      goto60=. 0
smoutput '(1) GOTO 70'
                    else.
smoutput '(3) GOTO 60'
                    end.
                    goto_u60.
                  end.
                  y=. 0 (< 1 , jch + 1 1) } y
                  jch=. >: jch
                end.
smoutput ": 'after loop 20 HT' ; y
smoutput '(2) GOTO 50'
              else.
smoutput 'ILAZRO OR ILAZR2 = .F.'
                NB. Only test 2 passed - chase the zero to
                NB. T[ilast,ilast], then process as in the
                NB. case T[ilast,ilast]=0
                jch=. j
                liosc=. (>: ilastm) th2lios <: j
                liosr=. (2 + j) th2lios ifrstm
                while. jch < ilast do.
smoutput 'loop 30 JCH = ',":jch
smoutput ": 'before rotation subHT' ; (< a: ; (jch + 0 1) ; liosc) { y
                  'y cs'=. (rot &. |:) rotga y ; (< 1 ; (jch + 0 1) ; (2 }. liosc)) ; < < a: ; 0
                  y=. (< 0 ; (jch + 0 1) ; liosc) (cs & (rot &. |:)) upd y
smoutput ": 'after rotation subHT' ; (< a: ; (jch + 0 1) ; liosc) { y
                  dQ=. dQ , (+ cs) , jch + 0 1
smoutput ": 'before rotation subHT' ; (< a: ; liosr ; (jch - 0 1)) { y
                  'y cs'=. rot rotga y ; (< 0 ; liosr ; (jch - 0 1)) ; _1
                  y=. (< 1 ; (_2 }. liosr) ; (jch - 0 1)) (cs & rot) upd y
smoutput ": 'after rotation subHT' ; (< a: ; liosr ; (jch - 0 1)) { y
                  dZ=. dZ , cs , jch - 0 1
                  liosc=. }. liosc
                  liosr=. liosr , 2 + jch
                  jch=. >: jch
                end.
smoutput ": 'after loop 30 HT' ; y
              end.
smoutput '(3) GOTO 50'
              goto_u50.
            elseif. ilazro do.
smoutput 'ILAZRO = .T.'
              ifirst=. j
              goto60=. 0
smoutput '(2) GOTO 70'
              goto_u60.
            end.
            j=. <: j
          end.
          NB. drop-through is impossible
          ((< _.) setdiag"2 y) ; ,~ a:  NB. set all eigenvalues to NaN
          return.
        else.
          y=. 0 (< 1 , ,~ ilast) } y
smoutput '(1) GOTO 50'
        end.
        label_u50.
smoutput 'label 50'
        NB. T[ilast,ilast]=0 - clear H[ilast,ilast-1] to
        NB. split off a 1x1 block
        lios=. (>: ilast) th2lios ifrstm
smoutput ": 'before rotation subHT' ; (< a: ; lios ; (ilast - 0 1)) { y
        'y cs'=. rot rotga y ; (< 0 ; lios ; (ilast - 0 1)) ; _1
        y=. (< 1 ; (}: lios) ; (ilast - 0 1)) (cs & rot) upd y
smoutput ": 'after rotation subHT' ; (< a: ; lios ; (ilast - 0 1)) { y
        dZ=. dZ , cs , ilast - 0 1
      else.
        y=. 0 (< 0 , ilast - 0 1) } y
smoutput '(2) GOTO 60'
      end.
    else.
smoutput '(1) GOTO 60'
    end.
    label_u60.
    if. goto60 do.
smoutput 'label 60'
      NB. H[ilast,ilast-1]=0 - standartize B, set alpha and
      NB. beta
      'y signbc'=. (ilast , 1) hgeqzxo y  NB. process ilast-th eigenvalue (column)
      dZ=. dZ , 4 {."1 signbc ,. ilast
      NB. goto next block - exit if finished
      ilast=. <: ilast
      if. ilast < h do.
        NB. normal exit
        'y signbc'=. (0 , 0 >. <: h) hgeqzxo y  NB. process eigenvalues (columns) 0:h-1
        dZ=. dZ , 4 {."1 signbc ,. i. 0 >. <: h
        y ; dQ ; dZ
        return.
      end.
      NB. reset counters
      iiter=. 0
      eshift=. 0
      'ifrstm ilastm'=. (h , ilast) reset (ifrstm , ilastm)
    else.
smoutput 'label 70'
      NB. QZ step
      NB. This iteration only involves rows/columns
      NB. ifirst:ilast. We assume ifirst<ilast, and that the
      NB. diagonal of B is non-zero
      iiter=. >: iiter
      ifrstm=. ifirst step ifrstm
      NB. Compute the shift.
      NB. At this point, ifirst<ilast, and the diagonal
      NB. elements of T[ifirst:ilast,ifirst:ilast] are larger
      NB. than btol in magnitude
      if. 10 | iiter do.
smoutput 'Wilkinson shift'
        NB. The Wilkinson shift, i.e., the eigenvalues of the
        NB. bottom-right 2x2 block of H*T^_1 which is nearest
        NB. to the bottom-right element.
        NB. We factor T as U*D, where U is unit upper
        NB. triangular, and compute (H*D^_1)*U^_1
smoutput ": 'abscale' ; abscale
smoutput ": 'HT[ilast-1,ilast-1]' ; (< a: ; ;~ ilast - 1) { y
        'U12 AD11 AD21 AD12 AD22'=. %/ (5 0 2 1 3 ,: 7 4 4 7 7) ({,) abscale * ((< a: ; ;~ ilast - 1 0) { y)
        ABI22=. AD22 - U12 * AD21
smoutput ": 'U12' ; U12 ; 'AD11' ; AD11 ; 'AD21' ; AD21 ; 'AD12' ; AD12 ; 'AD22' ; AD22 ; 'ABI22' ; ABI22
        t1=. -: AD11 + ABI22
        rtdisc=. %: (t1 , AD12 , -AD11) mp (t1 , AD21 , AD22)
        temp=. +/ (*/) +. rtdisc , t1 - ABI22
        shift=. t1 - temp condneg rtdisc
smoutput ": 't1' ; t1 ; 'rtdisc' ; rtdisc ; 'temp' ; temp ; 'shift' ; shift
      else.
        NB. Exceptional shift. Chosen for no paticularly good
        NB. reason
        eshift=. eshift + + %/ abscale * (;/ 0 1 ,. (_1 0 ,: _1 _1) + ilast) { y
smoutput ": 'eshift' ; eshift
        shift=. eshift
      end.
      NB. now check for two consecutive small subdiagonals
      HTd=. (0 _1 ,"0 1 ilast (] , -) ifirst) diag"1 2/ y
      ctemp=. (- (shift&*))/ abscale * {. HTd
      temp=. (sorim }."1 ctemp) ,: ascale * sorim (< 1 ; 0 ; <<0) { HTd
      tempr=. >./ temp
smoutput ": 'HTd' ; HTd ; 'ctemp' ; ctemp ; 'temp' ; temp ; 'tempr' ; tempr
      temp=. temp %"1 ((0 , 1 - FP_EPS) I. tempr) } 1 , tempr ,: 1
      'istart ctemp'=. (+&ifirst , {&ctemp) (ilast - ifirst) | >: (>:/ temp * atol ,: sorim (< 1 ; 0 ; <<_1) { HTd) i: 1
smoutput ": 'temp' ; temp ; 'istart' ; istart ; 'ctemp' ; ctemp
smoutput 'label 90'
      NB. do an implicit-shift QZ sweep
      NB. initial Q
smoutput ": 'fg' ; ctemp , ascale * (< 0 , istart + 1 0) { y
      cs=. }: lartg ctemp , ascale * (< 0 , istart + 1 0) { y
      NB. sweep
      j=. istart
      liosc=. (>: ilastm) th2lios j
      liosr=. (j + 2) th2lios ifrstm
      while. j < ilast do.
smoutput 'loop 150 J = ',":j
        lios=. j + 0 1
        NB. is a first iteration?
smoutput ": 'before rotation subHT' ; (< a: ; lios ; liosc) { y
        if. j = istart do.
          y=. (< a: ; lios ; liosc) (cs & (rot &. |:)"2) upd y
        else.
smoutput ": 'fg' ; (< 0 ; lios ; ({. liosc)) { y
          'y cs'=. (rot &. |:) rotga y ; (< 0 ; lios ; liosc) ; < < a: ; 0
          liosc=. }. liosc
          y=. (< 1 ; lios ; liosc) (cs & (rot &. |:)) upd y
        end.
smoutput ": 'cs' ; cs
smoutput ": 'after rotation subHT' ; (< a: ; lios ; liosc) { y
        dQ=. dQ , (+ cs) , lios
        lios=. j + 1 0
smoutput ": 'before rotation subHT' ; (< a: ; liosr ; lios) { y
        'y cs'=. rot rotga y ; (< 1 ; liosr ; lios) ; _1
        NB. isn't a last iteration?
        if. j < <: ilast do.
          liosr=. liosr , j + 2
        end.
        y=. (< 0 ; liosr ; lios) (cs & rot) upd y
smoutput ": 'after rotation subHT' ; (< a: ; liosr ; lios) { y
        dZ=. dZ , cs , lios
        j=. >: j
      end.
smoutput ": 'after loop 150 HT' ; y
    end.
    jiter=. >: jiter
  end.
  NB. drop-through means non-convergence, set incorrect eigenvalues 0:ilast to NaN
  ((_. ; 0 0 , >: ilast) setdiag"2 y) ; ,~ a:
)

NB. ---------------------------------------------------------
NB. hgezqe
NB. hgezqs
NB. hgeqze
NB. hgeqzs
NB.
NB. Description:
NB.   Shortcuts, see hgeqzxxxx
NB.
NB. Syntax:
NB.   'HTupd dQ dZ'=. hs hgexxx HT

hgezqe=: hgezqeo`[`(2 1&{@,`[@.((<{:)~{.))`[ hgezq
hgezqs=: hgezqso`]`]                      `] hgezq

hgeqze=: hgeqzeo`[`(2 1&{@,`[@.((<{:)~{.))`[ hgeqz
hgeqzs=: hgeqzso`]`]                      `] hgeqz

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. hgezqenn
NB. hgezqenv
NB. hgezqevn
NB. hgezqevv
NB. hgezqsnn
NB. hgezqsnv
NB. hgezqsvn
NB. hgezqsvv
NB.
NB. Description:
NB.   Eigenvalues and, optionally, generalized Schur form of
NB.   pair of structured matrices
NB.
NB. Syntax:
NB.   ab=.      hs hgezqenn H ,: T
NB.   'ab Z'=.  hs hgezqenv H , T ,: Z1
NB.   'ab Q'=.  hs hgezqevn H , T ,: Q1
NB.   'ab QZ'=. hs hgezqevv H , T , Q1 ,: Z1
NB.   SP=.      hs hgezqsnn H ,: T
NB.   SPZ=.     hs hgezqsnv H , T ,: Z1
NB.   SPQ=.     hs hgezqsvn H , T ,: Q1
NB.   SPQZ=.    hs hgezqsvv H , T , Q1 ,: Z1
NB. where
NB.   hs   - 2-vector of integers (h,s) 'head' and 'size',
NB.          defines submatrices H11 and T11 position in H
NB.          and T, respectively, see ggballp and gehrdl
NB.   H    - n×n-matrix, lower Hessenberg inside the
NB.          submatrix H[h:h+s-1,h:h+s-1], and lower
NB.          triangular outside
NB.   T    - n×n-matrix, lower triangular
NB.   ab   -: a ,: b
NB.   a    - n-vector, defines eigenvalues
NB.   b    - n-vector, defines eigenvalues
NB.   Q1   - n×n-matrix, the unitary (orthogonal)
NB.   Q    - n×n-matrix (ΔQ*Q1)
NB.   Z1   - n×n-matrix, the unitary (orthogonal)
NB.   Z    - n×n-matrix (ΔZ*Z1)
NB.   QZ   -: Q ,: Z
NB.   SP   -: S ,: P
NB.   SPQ  -: SP , Q
NB.   SPZ  -: SP , Z
NB.   SPQZ -: SPQ , Z
NB.   S    - n×n-matrix, lower triangular from the
NB.          generalized Schur factorization
NB.   P    - n×n-matrix, lower triangular from the
NB.          generalized Schur factorization
NB.
NB. Assertions:
NB.   ab -: diag"2 SP  NB. generalized eigenvalues alpha,beta
NB.                    NB. are defined by diagonals of S,P ,
NB.                    NB. respectively
NB.
NB. Notes:
NB. - non-converged eigenvalues are set to NaN

hgezqenn=:            diag"2 @ (0 {::                               hgezqe         )
hgezqenv=: (2 {  ]) ((diag"2 @ (0 {:: ])) ;  (rotscll  (2 & {::))) (hgezqe (2 & {.))
hgezqevn=: (2 {  ]) ((diag"2 @ (0 {:: ])) ;  (rotscll  (1 & {::))) (hgezqe (2 & {.))
hgezqevv=: (2 }. ]) ((diag"2 @ (0 {:: ])) ;  (rotscll"2 &.: > }.)) (hgezqe (2 & {.))

hgezqsnn=:                      0 {::                               hgezqs
hgezqsnv=: (2 {  ]) ((         (0 {   ])) ,: (rotscll  (2 & {::))) (hgezqs (2 & {.))
hgezqsvn=: (2 {  ]) ((         (0 {   ])) ,: (rotscll  (1 & {::))) (hgezqs (2 & {.))
hgezqsvv=: (2 }. ]) ((         (0 {:: ])) ,  (rotscll"2 &: > }. )) (hgezqs (2 & {.))

NB. ---------------------------------------------------------
NB. hgeqzenn
NB. hgeqzenv
NB. hgeqzevn
NB. hgeqzevv
NB. hgeqzsnn
NB. hgeqzsnv
NB. hgeqzsvn
NB. hgeqzsvv
NB.
NB. Description:
NB.   Eigenvalues and, optionally, Schur form:
NB.     H = Q * S * Z**H
NB.     T = Q * P * Z**H
NB.   of pair of structured matrices
NB.
NB. Syntax:
NB.   ab=.      hs hgeqzenn H ,: T
NB.   'ab Z'=.  hs hgeqzenv H , T ,: Z1
NB.   'ab Q'=.  hs hgeqzevn H , T ,: Q1
NB.   'ab QZ'=. hs hgeqzevv H , T , Q1 ,: Z1
NB.   SP=.      hs hgeqzsnn H ,: T
NB.   SPZ=.     hs hgeqzsnv H , T ,: Z1
NB.   SPQ=.     hs hgeqzsvn H , T ,: Q1
NB.   SPQZ=.    hs hgeqzsvv H , T , Q1 ,: Z1
NB. where
NB.   hs   - 2-vector of integers (h,s) 'head' and 'size',
NB.          defines submatrices H11 and T11 position in H
NB.          and T, respectively, see ggbalup and gehrdu
NB.   H    - n×n-matrix, upper Hessenberg inside the
NB.          submatrix H[h:h+s-1,h:h+s-1], and upper
NB.          triangular outside
NB.   T    - n×n-matrix, upper triangular
NB.   ab   -: a ,: b
NB.   a    - n-vector, defines eigenvalues
NB.   b    - n-vector, defines eigenvalues
NB.   Q1   - n×n-matrix, the unitary (orthogonal)
NB.   Q    - n×n-matrix (Q1*ΔQ)
NB.   Z1   - n×n-matrix, the unitary (orthogonal)
NB.   Z    - n×n-matrix (Z1*ΔZ)
NB.   QZ   -: Q ,: Z
NB.   SP   -: S ,: P
NB.   SPQ  -: SP , Q
NB.   SPZ  -: SP , Z
NB.   SPQZ -: SPQ , Z
NB.   S    - n×n-matrix, upper triangular from the
NB.          generalized Schur factorization
NB.   P    - n×n-matrix, upper triangular from the
NB.          generalized Schur factorization
NB.
NB. Notes:
NB. - non-converged eigenvalues are set to NaN
NB. - hgeqzenn models LAPACK's xHGEQZ('E','N','N')
NB. - hgeqzenv models LAPACK's xHGEQZ('E','N','V')
NB. - hgeqzevn models LAPACK's xHGEQZ('E','V','N')
NB. - hgeqzevv models LAPACK's xHGEQZ('E','V','V')
NB. - hgeqzsnn models LAPACK's xHGEQZ('S','N','N')
NB. - hgeqzsnv models LAPACK's xHGEQZ('S','N','V')
NB. - hgeqzsvn models LAPACK's xHGEQZ('S','V','N')
NB. - hgeqzsvv models LAPACK's xHGEQZ('S','V','V')
NB.
NB. Assertions:
NB.   ab -: diag"2 SP  NB. generalized eigenvalues alpha,beta
NB.                    NB. are defined by diagonals of S,P ,
NB.                    NB. respectively
NB.
NB. Application:
NB. - model LAPACK's xHGEQZ('N','I'):
NB.     gghrduni=: hgeqzenv (, (idmat @ c))
NB. - model LAPACK's xHGEQZ('I','N'):
NB.     gghrduin=: hgeqzevn (, (idmat @ c))
NB. - model LAPACK's xHGEQZ('I','I'):
NB.     gghrduii=: hgeqzevv (,~^:2~ (idmat @ c))
NB. - model LAPACK's xHGEQZ('I','V'):
NB.     gghrduiv=: hgeqzevv ((1 & A.) @ , (idmat @ c))
NB. - model LAPACK's xHGEQZ('V','I'):
NB.     gghrduvi=: hgeqzevv (, (idmat @ c))
NB. - detect case of non-convergence (0=converged,
NB.   1=non-converged), any of:
NB.     128!:5 < ab
NB.     128!:5 < SP
NB.     128!:5 < diag"2 SP
NB.
NB.   hgeqzsnn_mt_ 2 3 ; 0 {:: ggbalu_mt_ AB
NB.   ] 'Q1 R'=. ((tru_mt_@}:) ;~ ungqr_mt_) geqrf_mt_ (0;1) {:: ggbalu_mt_ AB
NB.   'HT Q Z'=. (Q1 ; idmat_mt_ 7) gghrdu_mt_ 2 3 ; R 1} 0 {:: ggbalu_mt_ AB
NB.   'SP dQ dZ'=. hgexxsnn hs ; HT

hgeqzenn=:            diag"2 @ (0 {::                               hgeqze         )
hgeqzenv=: (2 {  ]) ((diag"2 @ (0 {:: ])) ;  (rotsclu  (2 & {::))) (hgeqze (2 & {.))
hgeqzevn=: (2 {  ]) ((diag"2 @ (0 {:: ])) ;  (rotsclu  (1 & {::))) (hgeqze (2 & {.))
hgeqzevv=: (2 }. ]) ((diag"2 @ (0 {:: ])) ;  (rotsclu"2 &.: > }.)) (hgeqze (2 & {.))

hgeqzsnn=:                      0 {::                               hgeqzs
hgeqzsnv=: (2 {  ]) ((         (0 {   ])) ,: (rotsclu  (2 & {::))) (hgeqzs (2 & {.))
hgeqzsvn=: (2 {  ]) ((         (0 {   ])) ,: (rotsclu  (1 & {::))) (hgeqzs (2 & {.))
hgeqzsvv=: (2 }. ]) ((         (0 {:: ])) ,  (rotsclu"2 &: > }. )) (hgeqzs (2 & {.))

NB. =========================================================
NB. Test suite
