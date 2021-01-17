NB. Eigenvalues and Schur form
NB.
NB. hgexxe    Eigenvalues of pair of structured matrices
NB. hgexxs    Eigenvalues and the Schur form of pair of
NB.           structured matrices
NB.
NB. testhgeq  Test hgexxxxx by square matrices
NB. testeq    Adv. to make verb to test hgexxxxx by matrices
NB.           of generator and shape given
NB.
NB. Version: 0.11.0 2021-01-17
NB.
NB. Copyright 2011-2021 Igor Zhuravlov
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
  ((-&s1 upddiag ]) (mp (% norm1t)) ((-&s2 updl 0) {."1 ])) y
)

NB. ---------------------------------------------------------
NB. hgexxeo
NB.
NB. Description:
NB.   Compute generalized eigenvalues of hs-segment
NB.
NB. Syntax:
NB.   'HTupd signbc'=. hs hgexxeo H ,: T
NB. where
NB.   hs     - 2-vector of integers (h,s) 'head' and 'size',
NB.            defines eigenvalues range
NB.   H      - n×n-matrix, either lower or upper Hessenberg
NB.            inside the submatrix H[h:h+s-1,h:h+s-1], and
NB.            lower or upper triangular outside
NB.   T      - n×n-matrix, either lower or upper triangular
NB.   HTupd  -:Hupd ,: Tupd
NB.   Hupd   - n×n-matrix, being H with hs-segment of
NB.            diagonal replaced by alpha (see hgexx)
NB.   Tupd   - n×n-matrix, being T with hs-segment of
NB.            diagonal replaced by beta (see hgexx)
NB.   signbc - s-vector, scaling factors to form Q,Z later

hgexxeo=: 4 : 0
  'Hd Td'=. (0 , x) diag"2 y
  absb=. | Td
  'signbc Td'=. (,:~ FP_SFMIN < absb)} 1 0 2 |: (1 ,: + * Td) ,: 0 ,: absb
  ((((Hd * signbc) ,: Td) (;"1) 0 , x) setdiag"1 2 y) ; signbc
)

NB. ---------------------------------------------------------
NB. hgezqso
NB. hgeqzso
NB.
NB. Description:
NB.   Compute generalized eigenvalues of hs-segment and
NB.   reduce corresponding rows (hgezqso) or columns
NB.   (hgeqzso) to generalized Schur form
NB.
NB. Syntax:
NB.   'HTupd signbc'=. hs hgexxso H ,: T
NB. where
NB.   hs     - 2-vector of integers (h,s) 'head' and 'size',
NB.            defines eigenvalues range
NB.   H      - n×n-matrix, either lower (hgezqso) or upper
NB.            (hgeqzso) Hessenberg inside the submatrix
NB.            H[h:h+s-1,h:h+s-1], and lower (hgezqso) or
NB.            upper (hgeqzso) triangular outside
NB.   T      - n×n-matrix, either lower (hgezqso) or upper
NB.            (hgeqzso) triangular
NB.   HTupd  -:Hupd ,: Tupd
NB.   Hupd   - n×n-matrix, being H with rows (hgezqso) or
NB.            columns (hgeqzso) from hs-segment transformed
NB.            to Shur form
NB.   Tupd   - n×n-matrix, being T with rows (hgezqso) or
NB.            columns (hgeqzso) from hs-segment transformed
NB.            to Shur form
NB.   signbc - s-vector, scaling factors to form Q,Z later

hgezqso=: 4 : 0
  liso=. dhs2liso x
  'y signbc'=. x hgexxeo y
  subHT=. liso {"2 y
  (((,:~ liso >/ (i. c y))} subHT ,: subHT *"2 signbc) liso}"2 y) ; signbc
)

hgeqzso=: 4 : 0
  liso=. dhs2liso x
  'y signbc'=. x hgexxeo y
  subHT=. liso {"1 y
  (((,:~ (i. c y) </ liso)} subHT ,: subHT *"1 signbc) liso}"1 y) ; signbc
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
NB.   hgexxxo - monad to compute generalized eigenvalues of
NB.             hs-segment and, optionally, to reduce
NB.             corresponding columns to generalized Schur
NB.             form, is either hgezqso (hgezq only), hgeqzso
NB.             (hgeqz only) or hgexxeo, is called as:
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
NB.               'HTupd dQ1 dZ1'=. hs vapp HT
NB.             see hgeqzxxxx
NB.
NB. Notes:
NB. - non-converged eigenvalues are set to NaN

hgezq=: 1 : 0
:
  '`hgezqxo init reset step'=. m
  e=. +/ 'h s'=. x
  dQ1=. dZ1=. 0 4 $ 0
  abnorm=. (0 2 ,. ,.~ x) norms"2;.0 y
  'atol btol'=. abtol=. FP_SFMIN >. FP_PREC * abnorm
  'ascale bscale'=. abscale=. % FP_SFMIN >. abnorm
  'y signbc'=. ((c y) (] , -) e) hgezqxo y  NB. process eigenvalues (columns) h+s:n-1
  dZ1=. dZ1 , 4 {."1 signbc ,. (c y) th2liso e

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
            NB. test 1: H[j-1,j]=0 OR j=h
            if. j = h do.
              ilazro=. 1
            elseif. atol >: sorim (< 0 , j - 1 0) { y do.
              y=. 0 (< 0 , j - 1 0)} y
              ilazro=. 1
            else.
              ilazro=. 0
            end.
            NB. test 2: T[j,j]=0
            if. btol > | (< 1 , ,~ j) { y do.
              y=. 0 (< 1 , ,~ j)} y
              NB. test 2a: check for 2 consecutive small
              NB. superdiagonals in H
              ilazr2=. 0
              if. -. ilazro do.
                'Hj1j Hjj1 Hjj'=. sorim ((<"1) 0 ,. (_1 0,0 1,:0 0) + j) { y
                if. 0 >: (Hj1j , Hjj) mp (ascale * (Hjj1 , -atol)) do.
                  ilazr2=. 1
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
                jch=. j
                liso=. (>: ilastm) th2liso j
                while. jch < ilast do.
                  'y cs'=. rot&.|: rotga y ; (< 0 ; liso ; (jch + 0 1)) ; 0
                  liso=. }. liso
                  y=. (< 1 ; liso ; (jch + 0 1)) cs&(rot&.|:) upd y
                  dZ1=. dZ1 , (+ cs) , jch + 0 1
                  if. ilazr2 do.
                    y=. (< 0 , jch - 1 0) *&({. cs) upd y
                    ilazr2=. 0
                  end.
                  if. btol <: sorim (< 1 , jch + 1 1) { y do.
                    if. ilast > >: jch do.
                      ifirst=. >: jch
                      goto60=. 0
                    end.
                    goto_l60.
                  end.
                  y=. 0 (< 1 , jch + 1 1)} y
                  jch=. >: jch
                end.
              else.
                NB. Only test 2 passed - chase the zero to
                NB. T[ilast,ilast], then process as in the
                NB. case T[ilast,ilast]=0
                jch=. j
                lisor=. (>: ilastm) th2liso <: j
                lisoc=. (2 + j) th2liso ifrstm
                while. jch < ilast do.
                  'y cs'=. rot&.|: rotga y ; (< 1 ; (2 }. lisor) ; (jch + 0 1)) ; 0
                  y=. (< 0 ; lisor ; (jch + 0 1)) cs&(rot&.|:) upd y
                  dZ1=. dZ1 , (+ cs) , jch + 0 1
                  'y cs'=. rot rotga y ; (< 0 ; (jch - 0 1) ; lisoc) ; < < a: ; _1
                  y=. (< 1 ; (jch - 0 1) ; (_2 }. lisoc)) cs&rot upd y
                  dQ1=. dQ1 , cs , jch - 0 1
                  lisor=. }. lisor
                  lisoc=. lisoc , 2 + jch
                  jch=. >: jch
                end.
              end.
              goto_l50.
            elseif. ilazro do.
              ifirst=. j
              goto60=. 0
              goto_l60.
            end.
            j=. <: j
          end.
          NB. drop-through is impossible
          ((< _.) setdiag"2 y) ; ,~ a:  NB. set all eigenvalues to NaN
          return.
        else.
          y=. 0 (< 1 , ,~ ilast)} y
        end.
        label_l50.
        NB. T[ilast,ilast]=0 - clear H[ilast-1,ilast] to
        NB. split off a 1x1 block
        liso=. (>: ilast) th2liso ifrstm
        'y cs'=. rot rotga y ; (< 0 ; (ilast - 0 1) ; liso) ; < < a: ; _1
        y=. (< 1 ; (ilast - 0 1) ; (}: liso)) cs&rot upd y
        dQ1=. dQ1 , cs , ilast - 0 1
      else.
        y=. 0 (< 0 , ilast - 1 0)} y
      end.
    end.
    label_l60.
    if. goto60 do.
      NB. H[ilast-1,ilast]=0 - standartize B, set alpha and
      NB. beta
      'y signbc'=. (ilast , 1) hgezqxo y  NB. process ilast-th eigenvalue (column)
      dQ1=. dQ1 , 4 {."1 signbc ,. ilast
      NB. goto next block - exit if finished
      ilast=. <: ilast
      if. ilast < h do.
        NB. normal exit
        'y signbc'=. (0 , 0 >. <: h) hgezqxo y  NB. process eigenvalues (columns) 0:h-1
        dQ1=. dQ1 , 4 {."1 signbc ,. i. 0 >. <: h
        y ; dQ1 ; dZ1
        return.
      end.
      NB. reset counters
      iiter=. 0
      eshift=. 0
      'ifrstm ilastm'=. (h , ilast) reset (ifrstm , ilastm)
    else.
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
        NB. The Wilkinson shift, i.e., the eigenvalues of the
        NB. bottom-right 2x2 block of T^_1*H which is nearest
        NB. to the bottom-right element.
        NB. We factor T as D*L, where L is unit lower
        NB. triangular, and compute L^_1*(D^_1*H)
        'L21 DA11 DA12 DA21 DA22'=. %/ (6 0 1 2 3 ,: 7 4 4 7 7) ({,) abscale * ((< a: ; ;~ ilast - 1 0) { y)
        IBA22=. DA22 - L21 * DA12
        t1=. -: DA11 + IBA22
        rtdisc=. %: (t1 , DA21 , -DA11) mp (t1 , DA12 , DA22)
        temp=. +/ (*/) +. rtdisc , t1 - IBA22
        shift=. t1 - temp negneg rtdisc
      else.
        NB. Exceptional shift. Chosen for no paticularly good
        NB. reason
        eshift=. eshift + + %/ abscale * (;/ 0 1 ,. (0 _1 ,: _1 _1) + ilast) { y
        shift=. eshift
      end.
      NB. now check for two consecutive small subdiagonals
      HTd=. (0 1 ,"0 1 ilast (] , -) ifirst) diag"1 2/ y
      ctemp=. (- (shift&*))/ abscale * {. HTd
      temp=. (sorim }."1 ctemp) ,: ascale * sorim (< 1 ; 0 ; <<0) { HTd
      tempr=. >./ temp
      temp=. temp %"1 ((0 , 1 - FP_EPS) I. tempr)} 1 , tempr ,: 1
      'istart ctemp'=. (+&ifirst , {&ctemp) (ilast - ifirst) | >: (>:/ temp * atol ,: sorim (< 1 ; 0 ; <<_1) { HTd) i: 1
      NB. do an implicit-shift ZQ sweep
      NB. initial Z
      cs=. lartg ctemp , ascale * (< 0 , istart + 0 1) { y
      NB. sweep
      j=. istart
      lisor=. (>: ilastm) th2liso j
      lisoc=. (j + 2) th2liso ifrstm
      while. j < ilast do.
        liso=. j + 0 1
        NB. is a first iteration?
        if. j = istart do.
          y=. (< a: ; lisor ; liso) cs&(rot&.|:)"2 upd y
        else.
          'y cs'=. rot&.|: rotga y ; (< 0 ; lisor ; liso) ; 0
          lisor=. }. lisor
          y=. (< 1 ; lisor ; liso) cs&(rot&.|:) upd y
        end.
        dZ1=. dZ1 , (+ cs) , liso
        liso=. j + 1 0
        'y cs'=. rot rotga y ; (< 1 ; liso ; lisoc) ; < < a: ; _1
        NB. isn't a last iteration?
        if. j < <: ilast do.
          lisoc=. lisoc , j + 2
        end.
        y=. (< 0 ; liso ; lisoc) cs&rot upd y
        dQ1=. dQ1 , cs , liso
        j=. >: j
      end.
    end.
    jiter=. >: jiter
  end.
  NB. drop-through means non-convergence, set incorrect eigenvalues 0:ilast to NaN
  ((_. ; 0 0 , >: ilast) setdiag"2 y) ; dQ1 ; dZ1
)

hgeqz=: 1 : 0
:
  '`hgeqzxo init reset step'=. m
  e=. +/ 'h s'=. x
  dQ1=. dZ1=. 0 4 $ 0
  abnorm=. (0 2 ,. ,.~ x) norms"2;.0 y
  'atol btol'=. abtol=. FP_SFMIN >. FP_PREC * abnorm
  'ascale bscale'=. abscale=. % FP_SFMIN >. abnorm
  'y signbc'=. ((c y) (] , -) e) hgeqzxo y  NB. process eigenvalues (columns) h+s:n-1
  dZ1=. dZ1 , 4 {."1 signbc ,. (c y) th2liso e

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
            NB. test 1: H[j,j-1]=0 OR j=h
            if. j = h do.
              ilazro=. 1
            elseif. atol >: sorim (< 0 , j - 0 1) { y do.
              y=. 0 (< 0 , j - 0 1)} y
              ilazro=. 1
            else.
              ilazro=. 0
            end.
            NB. test 2: T[j,j]=0
            if. btol > | (< 1 , ,~ j) { y do.
              y=. 0 (< 1 , ,~ j)} y
              NB. test 2a: check for 2 consecutive small
              NB. subdiagonals in H
              ilazr2=. 0
              if. -. ilazro do.
                'Hjj1 Hj1j Hjj'=. sorim ((<"1) 0 ,. (0 _1,1 0,:0 0) + j) { y
                if. 0 >: (Hjj1 , Hjj) mp (ascale * (Hj1j , -atol)) do.
                  ilazr2=. 1
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
                jch=. j
                liso=. (>: ilastm) th2liso j
                while. jch < ilast do.
                  'y cs'=. rot rotga y ; (< 0 ; (jch + 0 1) ; liso) ; < < a: ; 0
                  liso=. }. liso
                  y=. (< 1 ; (jch + 0 1) ; liso) cs&rot upd y
                  dQ1=. dQ1 , (+ cs) , jch + 0 1
                  if. ilazr2 do.
                    y=. (< 0 , jch - 0 1) *&({. cs) upd y
                    ilazr2=. 0
                  end.
                  if. btol <: sorim (< 1 , jch + 1 1) { y do.
                    if. ilast > >: jch do.
                      ifirst=. >: jch
                      goto60=. 0
                    end.
                    goto_u60.
                  end.
                  y=. 0 (< 1 , jch + 1 1)} y
                  jch=. >: jch
                end.
              else.
                NB. Only test 2 passed - chase the zero to
                NB. T[ilast,ilast], then process as in the
                NB. case T[ilast,ilast]=0
                jch=. j
                lisoc=. (>: ilastm) th2liso <: j
                lisor=. (2 + j) th2liso ifrstm
                while. jch < ilast do.
                  'y cs'=. rot rotga y ; (< 1 ; (jch + 0 1) ; (2 }. lisoc)) ; < < a: ; 0
                  y=. (< 0 ; (jch + 0 1) ; lisoc) cs&rot upd y
                  dQ1=. dQ1 , (+ cs) , jch + 0 1
                  'y cs'=. rot&.|: rotga y ; (< 0 ; lisor ; (jch - 0 1)) ; _1
                  y=. (< 1 ; (_2 }. lisor) ; (jch - 0 1)) cs&(rot&.|:) upd y
                  dZ1=. dZ1 , cs , jch - 0 1
                  lisoc=. }. lisoc
                  lisor=. lisor , 2 + jch
                  jch=. >: jch
                end.
              end.
              goto_u50.
            elseif. ilazro do.
              ifirst=. j
              goto60=. 0
              goto_u60.
            end.
            j=. <: j
          end.
          NB. drop-through is impossible
          ((< _.) setdiag"2 y) ; ,~ a:  NB. set all eigenvalues to NaN
          return.
        else.
          y=. 0 (< 1 , ,~ ilast)} y
        end.
        label_u50.
        NB. T[ilast,ilast]=0 - clear H[ilast,ilast-1] to
        NB. split off a 1x1 block
        liso=. (>: ilast) th2liso ifrstm
        'y cs'=. rot&.|: rotga y ; (< 0 ; liso ; (ilast - 0 1)) ; _1
        y=. (< 1 ; (}: liso) ; (ilast - 0 1)) cs&(rot&.|:) upd y
        dZ1=. dZ1 , cs , ilast - 0 1
      else.
        y=. 0 (< 0 , ilast - 0 1)} y
      end.
    end.
    label_u60.
    if. goto60 do.
      NB. H[ilast,ilast-1]=0 - standartize B, set alpha and
      NB. beta
      'y signbc'=. (ilast , 1) hgeqzxo y  NB. process ilast-th eigenvalue (column)
      dZ1=. dZ1 , 4 {."1 signbc ,. ilast
      NB. goto next block - exit if finished
      ilast=. <: ilast
      if. ilast < h do.
        NB. normal exit
        'y signbc'=. (0 , 0 >. <: h) hgeqzxo y  NB. process eigenvalues (columns) 0:h-1
        dZ1=. dZ1 , 4 {."1 signbc ,. i. 0 >. <: h
        y ; dQ1 ; dZ1
        return.
      end.
      NB. reset counters
      iiter=. 0
      eshift=. 0
      'ifrstm ilastm'=. (h , ilast) reset (ifrstm , ilastm)
    else.
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
        NB. The Wilkinson shift, i.e., the eigenvalues of the
        NB. bottom-right 2x2 block of H*T^_1 which is nearest
        NB. to the bottom-right element.
        NB. We factor T as U*D, where U is unit upper
        NB. triangular, and compute (H*D^_1)*U^_1
        'U12 AD11 AD21 AD12 AD22'=. %/ (5 0 2 1 3 ,: 7 4 4 7 7) ({,) abscale * ((< a: ; ;~ ilast - 1 0) { y)
        ABI22=. AD22 - U12 * AD21
        t1=. -: AD11 + ABI22
        rtdisc=. %: (t1 , AD12 , -AD11) mp (t1 , AD21 , AD22)
        temp=. +/ (*/) +. rtdisc , t1 - ABI22
        shift=. t1 - temp negneg rtdisc
      else.
        NB. Exceptional shift. Chosen for no paticularly good
        NB. reason
        eshift=. eshift + + %/ abscale * (;/ 0 1 ,. (_1 0 ,: _1 _1) + ilast) { y
        shift=. eshift
      end.
      NB. now check for two consecutive small subdiagonals
      HTd=. (0 _1 ,"0 1 ilast (] , -) ifirst) diag"1 2/ y
      ctemp=. (- (shift&*))/ abscale * {. HTd
      temp=. (sorim }."1 ctemp) ,: ascale * sorim (< 1 ; 0 ; <<0) { HTd
      tempr=. >./ temp
      temp=. temp %"1 ((0 , 1 - FP_EPS) I. tempr)} 1 , tempr ,: 1
      'istart ctemp'=. (+&ifirst , {&ctemp) (ilast - ifirst) | >: (>:/ temp * atol ,: sorim (< 1 ; 0 ; <<_1) { HTd) i: 1
      NB. do an implicit-shift QZ sweep
      NB. initial Q
      cs=. lartg ctemp , ascale * (< 0 , istart + 1 0) { y
      NB. sweep
      j=. istart
      lisoc=. (>: ilastm) th2liso j
      lisor=. (j + 2) th2liso ifrstm
      while. j < ilast do.
        liso=. j + 0 1
        NB. is a first iteration?
        if. j = istart do.
          y=. (< a: ; liso ; lisoc) cs&rot"2 upd y
        else.
          'y cs'=. rot rotga y ; (< 0 ; liso ; lisoc) ; < < a: ; 0
          lisoc=. }. lisoc
          y=. (< 1 ; liso ; lisoc) cs&rot upd y
        end.
        dQ1=. dQ1 , (+ cs) , liso
        liso=. j + 1 0
        'y cs'=. rot&.|: rotga y ; (< 1 ; lisor ; liso) ; _1
        NB. isn't a last iteration?
        if. j < <: ilast do.
          lisor=. lisor , j + 2
        end.
        y=. (< 0 ; lisor ; liso) cs&(rot&.|:) upd y
        dZ1=. dZ1 , cs , liso
        j=. >: j
      end.
    end.
    jiter=. >: jiter
  end.
  NB. drop-through means non-convergence, set incorrect eigenvalues 0:ilast to NaN
  ((_. ; 0 0 , >: ilast) setdiag"2 y) ; dQ1 ; dZ1
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
NB.   'HTupd dQ1 dZ1'=. hs hgexxx HT

hgezqe=: hgexxeo`[`(2 1&{@,`[@.((<{:)~{.))`[ hgezq
hgezqs=: hgezqso`]`]                      `] hgezq

hgeqze=: hgexxeo`[`(2 1&{@,`[@.((<{:)~{.))`[ hgeqz
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
NB.   Compute eigenvalues (hgezqxxx) and reduce to
NB.   generalized Schur form (hgezqsxx):
NB.     dQ1^H * S * dZ1 = H
NB.     dQ1^H * P * dZ1 = T
NB.   the generalized lower Hessenberg form (H,T) using
NB.   single-shift ZQ method. Matrix pairs of this type are
NB.   produced by the reduction to generalized lower
NB.   Hessenberg form of a matrix pair (A,B):
NB.     Q1^H * H * Z1 = A
NB.     Q1^H * T * Z1 = B
NB.   as computed by gghrdlxx. The unitary (orthogonal)
NB.   matrices dQ1 and dZ1 may either be formed explicitly,
NB.   or they may be postmultiplied by input matrices Q1 and
NB.   Z1, so that:
NB.     (dQ1*Q1)^H * S * (dZ1*Z1) = Q1^H * H * Z1
NB.     (dQ1*Q1)^H * P * (dZ1*Z1) = Q1^H * T * Z1
NB.   To avoid overflow, eigenvalues of the matrix pair (H,T)
NB.   (equivalently, of (A,B)) are computed as a pair of
NB.   values. Each i-th eigenvector (row) from L and R has a
NB.   corresponding eigenvalue represented as a pair of i-th
NB.   elements from vectors e1 and e2:
NB.     E1=. diagmat e1=. diag S
NB.     E2=. diagmat e2=. diag P
NB.   If E2 is nonsingular then:
NB.     E=. diagmat e1%e2
NB.   is a diagonal matrix of eigenvalues, and GNEP can be
NB.   expressed as:
NB.     L * A = E * L * B
NB.     A * R^H = B * R^H * E
NB.   and if E1 is nonsingular then:
NB.     E=. diagmat e2%e1
NB.   is a diagonal matrix of eigenvalues, and GNEP can be
NB.   expressed as:
NB.     E * L * A = L * B
NB.     A * R^H * E = B * R^H * E
NB.
NB. Syntax:
NB.   e1e2=.        hs hgezqenn H ,: T
NB.   'e1e2 Z2'=.   hs hgezqenv H ,  T ,: Z1
NB.   'e1e2 Q2'=.   hs hgezqevn H ,  T ,: Q1
NB.   'e1e2 Q2Z2'=. hs hgezqevv H ,  T ,  Q1 ,: Z1
NB.   'S P'=.       hs hgezqsnn H ,: T
NB.   'S P Z2'=.    hs hgezqsnv H ,  T ,: Z1
NB.   'S P Q2'=.    hs hgezqsvn H ,  T ,: Q1
NB.   'S P Q2 Z2'=. hs hgezqsvv H ,  T ,  Q1 ,: Z1
NB. where
NB.   hs   - 2-vector of integers (h,s) 'head' and 'size',
NB.          defines submatrices H11 and T11 position in H
NB.          and T, respectively, see ggballp and gehrdl
NB.   H    - n×n-matrix, the lower Hessenberg inside the
NB.          submatrix H[h:h+s-1,h:h+s-1], and lower
NB.          triangular outside
NB.   T    - n×n-matrix, the lower triangular
NB.   e1e2 - 2×n-matrix of eigenvalues e1 and e2:
NB.            e1e2 -: e1 ,: e2
NB.   Q1   - n×n-matrix, the unitary (orthogonal)
NB.   Q2   - n×n-matrix, the unitary (orthogonal), the left
NB.          Schur vectors of (H,T) pair if Q1=I, and of
NB.          (A,B) pair otherwise
NB.   Z1   - n×n-matrix, the unitary (orthogonal)
NB.   Z2   - n×n-matrix, the unitary (orthogonal), the right
NB.          Schur vectors of (H,T) pair if Z1=I, and of
NB.          (A,B) pair otherwise
NB.   Q2Z2 -:Q2 ,: Z2
NB.   S    - n×n-matrix, the lower triangular
NB.   P    - n×n-matrix, the lower triangular
NB.
NB. Notes:
NB. - non-converged eigenvalues are set to NaN
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   e1e2     -: diag"2 (S ,: P)
NB.   Q2       -: dQ1 mp Q1
NB.   Z2       -: dZ1 mp Z1
NB.   (H ,: T) -: dQ1 (mp~ ct)~"2 (S ,: P) mp"2 dZ1
NB.   (C ,: D) -: Q2 (mp~ ct)~"2 (S ,: P) mp"2 Z2
NB.   I        -: Q2^H * Q2
NB.   I        -: Z2^H * Z2
NB. where
NB.   C - n×n-matrix, general
NB.   D - n×n-matrix, general
NB.   n=. # C
NB.   hs=. 0 , n
NB.   I=. idmat n
NB.   'B Z0'=. (trl@:(}:"1) ,: unglq)@gelqf D
NB.   A=. C (mp ct) Z0
NB.   'H T Q1 Z1'=. hs gghrdlvv A , B , I ,: Z0
NB.   e1e2=. hs hgezqenn H ,: T
NB.   'S P Q2 Z2'=. hs hgezqsvv H , T , Q1 ,: Z1
NB.   'S P dQ1 dZ1'=. hs hgezqsvv H , T , ,:~ I
NB.
NB. Application:
NB. - detect case of non-convergence (0=converged,
NB.   1=non-converged), any of:
NB.     128!:5 < e1e2
NB.     128!:5 < S,:P         NB. too expensive, use the next
NB.     128!:5 < diag"2 S,:P

hgezqenn=:            diag"2@(0 {::                            hgezqe     )
hgezqenv=: (2 {  ]) ((diag"2@(0 {:: ])) ; (rotscll  2&{::  )) (hgezqe 2&{.)
hgezqevn=: (2 {  ]) ((diag"2@(0 {:: ])) ; (rotscll  1&{::  )) (hgezqe 2&{.)
hgezqevv=: (2 }. ]) ((diag"2@(0 {:: ])) ; (rotscll"2&.:> }.)) (hgezqe 2&{.)

hgezqsnn=:                    0 {::                            hgezqs
hgezqsnv=: (2 {  ]) ((        0 {:: ] ) , (rotscll  2&{::  )) (hgezqs 2&{.)
hgezqsvn=: (2 {  ]) ((        0 {:: ] ) , (rotscll  1&{::  )) (hgezqs 2&{.)
hgezqsvv=: (2 }. ]) ((        0 {:: ] ) , (rotscll"2&: > }.)) (hgezqs 2&{.)

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
NB.   Compute eigenvalues (hgeqzxxx) and reduce to
NB.   generalized Schur form (hgeqzsxx):
NB.     dQ1 * S * dZ1^H = H
NB.     dQ1 * P * dZ1^H = T
NB.   the generalized upper Hessenberg form (H,T) using
NB.   single-shift QZ method. Matrix pairs of this type are
NB.   produced by the reduction to generalized upper
NB.   Hessenberg form of a matrix pair (A,B):
NB.     Q1 * H * Z1^H = A
NB.     Q1 * T * Z1^H = B
NB.   as computed by gghrduxx. The unitary (orthogonal)
NB.   matrices dQ1 and dZ1 may either be formed explicitly,
NB.   or they may be premultiplied by input matrices Q1 and
NB.   Z1, so that:
NB.     (Q1*dQ1) * S * (Z1*dZ1)^H = Q1 * H * Z1^H
NB.     (Q1*dQ1) * P * (Z1*dZ1)^H = Q1 * T * Z1^H
NB.   To avoid overflow, eigenvalues of the matrix pair (H,T)
NB.   (equivalently, of (A,B)) are computed as a pair of
NB.   values. Each i-th eigenvector (column) from L and R has
NB.   a corresponding eigenvalue represented as a pair of
NB.   i-th elements from vectors e1 and e2:
NB.     E1=. diagmat e1=. diag S
NB.     E2=. diagmat e2=. diag P
NB.   If E2 is nonsingular then:
NB.     E=. diagmat e1%e2
NB.   is a diagonal matrix of eigenvalues, and GNEP can be
NB.   expressed as:
NB.     L^H * A = E * L^H * B
NB.     A * R = B * R * E
NB.   and if E1 is nonsingular then:
NB.     E=. diagmat e2%e1
NB.   is a diagonal matrix of eigenvalues, and GNEP can be
NB.   expressed as:
NB.     E * L^H * A = L^H * B
NB.     A * R * E = B * R * E
NB.
NB. Syntax:
NB.   e1e2=.        hs hgeqzenn H ,: T
NB.   'e1e2 Z2'=.   hs hgeqzenv H ,  T ,: Z1
NB.   'e1e2 Q2'=.   hs hgeqzevn H ,  T ,: Q1
NB.   'e1e2 Q2Z2'=. hs hgeqzevv H ,  T ,  Q1 ,: Z1
NB.   'S P'=.       hs hgeqzsnn H ,: T
NB.   'S P Z2'=.    hs hgeqzsnv H ,  T ,: Z1
NB.   'S P Q2'=.    hs hgeqzsvn H ,  T ,: Q1
NB.   'S P Q2 Z2'=. hs hgeqzsvv H ,  T ,  Q1 ,: Z1
NB. where
NB.   hs   - 2-vector of integers (h,s) 'head' and 'size',
NB.          defines submatrices H11 and T11 position in H
NB.          and T, respectively, see ggbalup and gehrdu
NB.   H    - n×n-matrix, the upper Hessenberg inside the
NB.          submatrix H[h:h+s-1,h:h+s-1], and upper
NB.          triangular outside
NB.   T    - n×n-matrix, the upper triangular
NB.   e1e2 - 2×n-matrix of eigenvalues e1 and e2:
NB.            e1e2 -: e1 ,: e2
NB.   Q1   - n×n-matrix, the unitary (orthogonal)
NB.   Q2   - n×n-matrix, the unitary (orthogonal), the left
NB.          Schur vectors of (H,T) pair if Q1=I, and of
NB.          (A,B) pair otherwise
NB.   Z1   - n×n-matrix, the unitary (orthogonal)
NB.   Z2   - n×n-matrix, the unitary (orthogonal), the right
NB.          Schur vectors of (H,T) pair if Z1=I, and of
NB.          (A,B) pair otherwise
NB.   Q2Z2 -:Q2 ,: Z2
NB.   S    - n×n-matrix, the upper triangular
NB.   P    - n×n-matrix, the upper triangular
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
NB. Assertions (with appropriate comparison tolerance):
NB.   e1e2     -: diag"2 (S ,: P)
NB.   Q2       -: Q1 mp dQ1
NB.   Z2       -: Z1 mp dZ1
NB.   (H ,: T) -: dQ1 mp"2 (S ,: P) (mp ct)"2 dZ1
NB.   (C ,: D) -: Q2 mp"2 (S ,: P) (mp ct)"2 Z2
NB.   I        -: Q2 * Q2^H
NB.   I        -: Z2 * Z2^H
NB. where
NB.   C - n×n-matrix, general
NB.   D - n×n-matrix, general
NB.   n=. # C
NB.   hs=. 0 , n
NB.   I=. idmat n
NB.   'Q0 B'=. (ungqr ,: tru@}:)@geqrf D
NB.   A=. Q0 (mp~ ct)~ C
NB.   'H T Q1 Z1'=. hs gghrduvv A , B , Q0 ,: I
NB.   e1e2=. hs hgeqzenn H ,: T
NB.   'S P Q2 Z2'=. hs hgeqzsvv H , T , Q1 ,: Z1
NB.   'S P dQ1 dZ1'=. hs hgeqzsvv H , T , ,:~ I
NB.
NB. Application:
NB. - models LAPACK's xHGEQZ('E','N','I'):
NB.     NB. 'e1e2 dZ1'=. hs hgeqzeni H ,: T
NB.     hgeqzeni=: hgeqzenv (, idmat@c)
NB. - models LAPACK's xHGEQZ('E','I','N'):
NB.     NB. 'e1e2 dQ1'=. hs hgeqzein H ,: T
NB.     hgeqzein=: hgeqzevn (, idmat@c)
NB. - models LAPACK's xHGEQZ('E','I','I'):
NB.     NB. 'e1e2 dQ1dZ1'=. hs hgeqzeii H ,: T
NB.     hgeqzeii=: hgeqzevv (,~^:2~ idmat@c)
NB. - models LAPACK's xHGEQZ('E','I','V'):
NB.     NB. 'e1e2 dQ1Z2'=. hs hgeqzeiv H , T ,: Z1
NB.     hgeqzeiv=: hgeqzevv (1&A.@, idmat@c)
NB. - models LAPACK's xHGEQZ('E','V','I'):
NB.     NB. 'e1e2 Q2dZ1'=. hs hgeqzevi H , T ,: Q1
NB.     hgeqzevi=: hgeqzevv (, idmat@c)
NB. - models LAPACK's xHGEQZ('S','N','I'):
NB.     NB. 'S P dZ1'=. hs hgeqzsni H ,: T
NB.     hgeqzsni=: hgeqzsnv (, idmat@c)
NB. - models LAPACK's xHGEQZ('S','I','N'):
NB.     NB. 'S P dQ1'=. hs hgeqzsin H ,: T
NB.     hgeqzsin=: hgeqzsvn (, idmat@c)
NB. - models LAPACK's xHGEQZ('S','I','I'):
NB.     NB. 'S P dQ1 dZ1'=. hs hgeqzsii H ,: T
NB.     hgeqzsii=: hgeqzsvv (,~^:2~ idmat@c)
NB. - models LAPACK's xHGEQZ('S','I','V'):
NB.     NB. 'S P dQ1 Z2'=. hs hgeqzsiv H , T ,: Z1
NB.     hgeqzsiv=: hgeqzsvv (1&A.@, idmat@c)
NB. - models LAPACK's xHGEQZ('S','V','I'):
NB.     NB. 'S P Q2 dZ1'=. hs hgeqzsvi H , T ,: Q1
NB.     hgeqzsvi=: hgeqzsvv (, idmat@c)
NB. - detect case of non-convergence (0=converged,
NB.   1=non-converged), any of:
NB.     128!:5 < e1e2
NB.     128!:5 < S,:P         NB. too expensive, use the next
NB.     128!:5 < diag"2 S,:P
NB.
NB. References:
NB. [1] C. B. Moler, G. W. Stewart. An Algorithm for
NB.     Generalized Matrix Eigenvalue Problems. SIAM J.
NB.     Numer. Anal., 10(1973), pp. 241-256.

hgeqzenn=:            diag"2@(0 {::                            hgeqze     )
hgeqzenv=: (2 {  ]) ((diag"2@(0 {:: ])) ; (rotsclu  2&{::  )) (hgeqze 2&{.)
hgeqzevn=: (2 {  ]) ((diag"2@(0 {:: ])) ; (rotsclu  1&{::  )) (hgeqze 2&{.)
hgeqzevv=: (2 }. ]) ((diag"2@(0 {:: ])) ; (rotsclu"2&.:> }.)) (hgeqze 2&{.)

hgeqzsnn=:                    0 {::                            hgeqzs
hgeqzsnv=: (2 {  ]) ((        0 {:: ] ) , (rotsclu  2&{::  )) (hgeqzs 2&{.)
hgeqzsvn=: (2 {  ]) ((        0 {:: ] ) , (rotsclu  1&{::  )) (hgeqzs 2&{.)
hgeqzsvv=: (2 }. ]) ((        0 {:: ] ) , (rotsclu"2&: > }.)) (hgeqzs 2&{.)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testhgeq
NB.
NB. Description:
NB.   Test:
NB.   - xHGEQZ (math/lapack2 addon)
NB.   - hgexxxxx (math/mt addon)
NB.   by pair of square matrices
NB.
NB. Syntax:
NB.   testhgeq AB
NB. where
NB.   AB - 2×n×n-brick

testhgeq=: 3 : 0
  load_mttmp_ :: ] 'math/mt/test/lapack2/hgeqz'

  n=. c y
  hs=. 0 , n
  I=. idmat n

  'Hl Tl'=. HTl=. hs gghrdlnn (((mp  ct@unglq) ,: trlpick@(_1 }."1 ])) gelqf)/ y
  'Hu Tu'=. HTu=. hs gghrdunn (((mp~ ct@ungqr) ,: trupick@(_1 }.   ])) geqrf)/ y

  rcondl=. (geconi Hl) <. trlconi Tl
  rcondu=. (gecon1 Hu) <. trucon1 Tu

  normsl=. ;/ normi"2 HTl
  normsu=. ;/ norm1"2 HTu

  argslapack=. normsu , ;/ HTu , ,:~ I  NB. arguments for xHGEQZ            t511u t513u
  argsmtl=.    normsl ,  < HTl          NB. arguments for hgezqxnn          t511l t513l
  argsmtvl=.   normsl ,  < HTl ,     I  NB. arguments for hgezqxnv hgezqxvn t511l t513l
  argsmtvvl=.  normsl ,  < HTl , ,:~ I  NB. arguments for hgezqxvv          t511l t513l
  argsmtu=.    normsu ,  < HTu          NB. arguments for hgeqzxnn          t511u t513u
  argsmtvu=.   normsu ,  < HTu ,     I  NB. arguments for hgeqzxnv hgeqzxvn t511u t513u
  argsmtvvu=.  normsu ,  < HTu , ,:~ I  NB. arguments for hgeqzxvv          t511u t513u

  t511u1=: (t511u"1~ (  2             0      ,:  3            1)&{  )~      (0 4 5 ,: 1 4 5)&{
  t511l2=: (t511l"1~ (((2 ; 0)&{::) ; 0&{::) ,: (2 ; 1)&{:: ; 1 &{::)~ <"2@((0 2 3 ,: 1 2 3)&{)
  t511u2=: (t511u"1~ (((2 ; 0)&{::) ; 0&{::) ,: (2 ; 1)&{:: ; 1 &{::)~ <"2@((0 2 3 ,: 1 2 3)&{)

  t513u4=:               t513u@(4 {:: ])
  t513u5=:               t513u@(5 {:: ])
  t513u45=: (4 {:: ]) >.&t513u  5 {:: ]

  t513l1=:                   t513l@(1     {:: ])
  t513u1=:                   t513u@(1     {:: ])
  t513l2=:                   t513l@(2     {   ])
  t513u2=:                   t513u@(2     {   ])
  t513l01=: ((1;0) {:: ]) >.&t513l  (1;1) {:: ]
  t513u01=: ((1;0) {:: ]) >.&t513u  (1;1) {:: ]
  t513l23=: (2     {   ]) >.&t513l  3     {   ]
  t513u23=: (2     {   ]) >.&t513u  3     {   ]

  ('''enn''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(_."_                ))) argslapack
  ('''eni''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u5 ))) argslapack
  ('''env''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u5 ))) argslapack
  ('''ein''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u4 ))) argslapack
  ('''eii''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u45))) argslapack
  ('''eiv''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u45))) argslapack
  ('''evn''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u4 ))) argslapack
  ('''evi''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u45))) argslapack
  ('''evv''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u45))) argslapack

  ('''snn''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(_."_                ))) argslapack
  ('''sni''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u5 ))) argslapack
  ('''snv''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u5 ))) argslapack
  ('''sin''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u4 ))) argslapack
  ('''sii''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(t511u1 >./@, t513u45))) argslapack
  ('''siv''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(t511u1 >./@, t513u45))) argslapack
  ('''svn''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u4 ))) argslapack
  ('''svi''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(t511u1 >./@, t513u45))) argslapack
  ('''svv''&dhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(t511u1 >./@, t513u45))) argslapack

  ('''enn''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(_."_                ))) argslapack
  ('''eni''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u5 ))) argslapack
  ('''env''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u5 ))) argslapack
  ('''ein''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u4 ))) argslapack
  ('''eii''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u45))) argslapack
  ('''eiv''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u45))) argslapack
  ('''evn''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u4 ))) argslapack
  ('''evi''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u45))) argslapack
  ('''evv''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u45))) argslapack

  ('''snn''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(_."_                ))) argslapack
  ('''sni''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u5 ))) argslapack
  ('''snv''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u5 ))) argslapack
  ('''sin''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u4 ))) argslapack
  ('''sii''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(t511u1 >./@, t513u45))) argslapack
  ('''siv''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(t511u1 >./@, t513u45))) argslapack
  ('''svn''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(             t513u4 ))) argslapack
  ('''svi''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(t511u1 >./@, t513u45))) argslapack
  ('''svv''&zhgeqz_mttmp_' tmonad ((0 1}~ 1 ; #@(2&{::))  `]`(rcondu"_)`(_."_)`(t511u1 >./@, t513u45))) argslapack

  ('hgezqenn'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondl"_)`(_."_)`(_."_                ))) argsmtl
  ('hgezqenv'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondl"_)`(_."_)`(             t513l1 ))) argsmtvl
  ('hgezqevn'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondl"_)`(_."_)`(             t513l1 ))) argsmtvl
  ('hgezqevv'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondl"_)`(_."_)`(             t513l01))) argsmtvvl
  ('hgezqsnn'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondl"_)`(_."_)`(_."_                ))) argsmtl
  ('hgezqsnv'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondl"_)`(_."_)`(             t513l2 ))) argsmtvl
  ('hgezqsvn'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondl"_)`(_."_)`(             t513l2 ))) argsmtvl
  ('hgezqsvv'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondl"_)`(_."_)`(t511l2 >./@, t513l23))) argsmtvvl

  ('hgeqzenn'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondu"_)`(_."_)`(_."_                ))) argsmtu
  ('hgeqzenv'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondu"_)`(_."_)`(             t513u1 ))) argsmtvu
  ('hgeqzevn'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondu"_)`(_."_)`(             t513u1 ))) argsmtvu
  ('hgeqzevv'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondu"_)`(_."_)`(             t513u01))) argsmtvvu
  ('hgeqzsnn'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondu"_)`(_."_)`(_."_                ))) argsmtu
  ('hgeqzsnv'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondu"_)`(_."_)`(             t513u2 ))) argsmtvu
  ('hgeqzsvn'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondu"_)`(_."_)`(             t513u2 ))) argsmtvu
  ('hgeqzsvv'              tdyad  ((0 , c@(2&{::))`(2&{::)`]`(rcondu"_)`(_."_)`(t511u2 >./@, t513u23))) argsmtvvu

  coerase < 'mttmp'
  erase 't511u1 t511l2 t511u2 t513u4 t513u5 t513u45 t513l1 t513u1 t513l2 t513u2 t513l01 t513u01 t513l23 t513u23'

  EMPTY
)

NB. ---------------------------------------------------------
NB. testeq
NB.
NB. Description:
NB.   Adv. to make verb to test hgexxxxx by matrices of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testeq
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random square real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     ?@$&0 testeq_mt_ 150 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testeq_mt_ 150 150
NB. - test by random square complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testeq_mt_ 150 150

testeq=: 1 : 'EMPTY [ testhgeq_mt_@u@(2&,)^:(=/)'
