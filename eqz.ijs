NB. Eigenvalues and Schur form of pair of structured matrices
NB.
NB. hgeqzxe     Eigenvalues of pair of structured matrices
NB. hgeqzxs     Eigenvalues and the Schur form of pair of
NB.             structured matrices
NB.
NB. testhgeqze  Test hgeqzxe by general matrices given
NB. testhgeqzs  Test hgeqzxs by general matrices given
NB. testeqz     Adv. to make verb to test hgeqzxx by matrices
NB.             of generator and shape given
NB.
NB. Version: 0.6.8 2010-10-14
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
NB. hgeqzleo
NB. hgeqzueo
NB.
NB. Description:
NB.   Calculate generalized eigenvalues of hs-segment
NB.
NB. Syntax:
NB.   'HTupd signbc'=. hgeqzxeo hs ; HT
NB. where
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines eigenvalues range
NB.   HT    -: H ,: T
NB.   H     - n×n-matrix, either lower (hgeqzleo) or upper
NB.           (hgeqzueo) Hessenberg inside the submatrix
NB.           H[h:h+s-1,h:h+s-1], and lower (hgeqzlenn) or
NB.           upper (hgeqzuenn) triangular outside
NB.   T     - n×n-matrix, either lower (hgeqzleo) or upper
NB.           (hgeqzueo) triangular
NB.   HTupd -: Hupd ,: Tupd
NB.   Hupd  - n×n-matrix, being H with hs-segment of diagonal
NB.           replaced by alpha (see hgeqzx)
NB.   Tupd  - n×n-matrix, being T with hs-segment of diagonal
NB.           replaced by beta (see hgeqzx)
NB.   signbc - s-vector, scaling factors to form Q,Z later

hgeqzueo=: 3 : 0
  'hs HT'=. y
  'Hd Td'=. (0 , hs) diag"2 HT
  absb=. | Td
  'signbc Td'=. (,:~ FP_SFMIN < absb) } (1 ,: + * Td) ,: (0 ,: absb)
  ((((Hd * signbc) ,: Td) (;"1) 0 , hs) setdiag"1 2 HT) ; signbc
)

NB. ---------------------------------------------------------
NB. hgeqzlso
NB. hgeqzuso
NB.
NB. Description:
NB.   Calculate generalized eigenvalues of hs-segment and
NB.   reduce corresponding columns to generalized Schur form
NB.
NB. Syntax:
NB.   'HTupd signbc'=. hgeqzxso hs ; HT
NB. where
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines eigenvalues range
NB.   HT    -: H ,: T
NB.   H     - n×n-matrix, either lower (hgeqzleo) or upper
NB.           (hgeqzueo) Hessenberg inside the submatrix
NB.           H[h:h+s-1,h:h+s-1], and lower (hgeqzlenn) or
NB.           upper (hgeqzuenn) triangular outside
NB.   T     - n×n-matrix, either lower (hgeqzleo) or upper
NB.           (hgeqzueo) triangular
NB.   HTupd -: Hupd ,: Tupd
NB.   Hupd  - n×n-matrix, being H with either rows (hgeqzlso)
NB.           or columns (hgeqzuso) from hs-segment
NB.           transformed to Shur form
NB.   Tupd  - n×n-matrix, being T with either rows (hgeqzlso)
NB.           or columns (hgeqzuso) from hs-segment
NB.           transformed to Shur form
NB.   signbc - s-vector, scaling factors to form Q,Z later

hgeqzuso=: 3 : 0
  'hs HT'=. y
  lios=. dhs2lios hs
  'HT signbc'=. hgeqzueo y
  subHT=. lios {"1 HT
  (((,:~ (i. c HT) </ lios) } subHT ,: subHT *"1 signbc) ((lios }"1) dbg 'alpha&beta subHT') HT) ; signbc
)

NB. ---------------------------------------------------------
NB. hgeqzlnn
NB. hgeqzunn
NB.
NB. Description:
NB.   Adv. to make verbs to find eigenvalues of lower (upper)
NB.   Hessenberg-triangular pair (H,T) and, optionally, to
NB.   reduce it to generalized Schur form
NB.
NB. Syntax:
NB.   vapp=. hgeqzxxo`init`reset`step hgeqzxnn
NB. where
NB.   hgeqzxxo - monad to calculate generalized eigenvalues
NB.              of hs-segment and, optionally, to reduce
NB.              corresponding columns to generalized Schur
NB.              form, is either hgeqzxeo or hgeqzxso, is
NB.              called as:
NB.                'HTupd signbc'=. hgeqzxxo hs ; HT
NB.   init     - dyad to initialize counters, is called as:
NB.                'ifrstm ilastm'=. (h , ilast) init (0 , n-1)
NB.   reset    - dyad to reset counters optionally, is called
NB.              as:
NB.                'ifrstm ilastm'=. (h , ilast) reset (ifrstm , ilastm)
NB.   step     - monad to change counter optionally, is
NB.              called as:
NB.                ifrstm=. ifirst step ifrstm
NB.   vapp     - monad to find eigenvalues of lower (upper)
NB.              Hessenberg-triangular pair (H,T) and,
NB.              optionally, to reduce it to generalized
NB.              Schur form, is called as:
NB.                'HTupd dQ dZ'=. vapp hs ; HT
NB.              see hgeqzxxnn
NB.
NB. Notes:
NB. - non-converged eigenvalues will be set to NaN

hgeqzunn=: 1 : 0
  '`hgeqzuxo init reset step'=. m
  'hs HT'=. y
  e=. +/ 'h s'=. hs
  dQ=. dZ=. 4 0 $ 0
  abnorm=. (0 2 ,. ,.~ hs) norms"2 ;. 0 HT
  'atol btol'=. abtol=. FP_SFMIN >. FP_PREC * abnorm
  'ascale bscale'=. abscale=. % FP_SFMIN >. abnorm
  'HT signbc'=. hgeqzuxo ((c HT) (] , -) e) ; HT  NB. process eigenvalues (columns) h+s:n-1
  dZ=. dZ , 4 {."1 signbc ,. (c HT) th2lios e

  NB. Eigenvalues h+s:n-1 have been found.

  NB. Initialize dynamic indices
  ilast=. <: e
  'ifrstm ilastm'=. (h , ilast) init (0 , <: c HT)
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
      if. atol < sorim (< 0 , ilast - 0 1) { HT do.
        if. btol < | (< 1 , ,~ ilast) { HT do.
          NB. general case: j < ilast
          j=. <: ilast
          while. j >: h do.
smoutput 'loop 40 J = ',":j
            NB. test 1: H[j,j-1]=0 OR j=h
            if. j = h do.
              ilazro=. 1
smoutput 'ILAZRO = .T.'
            elseif. atol >: sorim (< 0 , j - 0 1) { HT do.
              HT=. 0 (< 0 , j - 0 1) } HT
              ilazro=. 1
smoutput 'ILAZRO = .T.'
            elseif. do.
              ilazro=. 0
smoutput 'ILAZRO = .F.'
            end.
            NB. test 2: T[j,j]=0
            if. btol > | (< 1 , ,~ j) { HT do.
smoutput 'Test 2 = .T.'
              HT=. 0 (< 1 , ,~ j) } HT
              NB. test 2a: check for 2 consecutive small subdiagonals in H
              ilazr2=. 0
              if. -. ilazro do.
                'Hjj1 Hj1j Hjj'=. sorim ((<"1) 0 ,. (0 _1,1 0,:0 0) + j) { HT
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
                  'HT cs'=. (((rot &. |:) rotga) dbg 'rotga(H)') HT ; (< 0 ; (jch + 0 1) ; lios) ; < < a: ; 0
                  lios=. }. lios
                  HT=. (< 1 ; (jch + 0 1) ; lios) (((cs & (rot &. |:)) upd) dbg 'rot(T)') HT
                  dQ=. dQ , (+ cs) , jch + 0 1
                  if. ilazr2 do.
                    HT=. (< 0 , jch - 0 1) (* & ({. cs)) upd HT
                    ilazr2=. 0
                  end.
                  if. btol <: sorim (< 1 , jch + 1 1) { HT do.
                    if. ilast > >: jch do.
                      ifirst=. >: jch
                      goto60=. 0
smoutput '(1) GOTO 70'
                    else.
smoutput '(3) GOTO 60'
                    end.
                    goto_60.
                  end.
                  HT=. 0 (< 1 , jch + 1 1) } HT
                  jch=. >: jch
                end.
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
                  'HT cs'=. (((rot &. |:) rotga) dbg 'rotga(T)') HT ; (< 1 ; (jch + 0 1) ; (2 }. liosc)) ; < < a: ; 0
                  HT=. (< 0 ; (jch + 0 1) ; liosc) (((cs & (rot &. |:)) upd) dbg 'rot(H)') HT
                  dQ=. dQ , (+ cs) , jch + 0 1
                  'HT cs'=. ((rot rotga) dbg 'rotga(H)') HT ; (< 0 ; liosr ; (jch - 0 1)) ; _1
                  HT=. (< 1 ; (_2 }. liosr) ; (jch - 0 1)) (((cs & rot) upd) dbg 'rot(T)') HT
                  dZ=. dZ , cs , jch - 0 1
                  liosc=. }. liosc
                  liosr=. liosr , 2 + jch
                  jch=. >: jch
                end.
              end.
smoutput '(3) GOTO 50'
              goto_50.
            elseif. ilazro do.
smoutput 'ILAZRO = .T.'
              ifirst=. j
              goto60=. 0
smoutput '(2) GOTO 70'
              goto_60.
            end.
            j=. <: j
          end.
          NB. drop-through is impossible
          ((< _.) setdiag"2 HT) ; ,~ a:  NB. set all eigenvalues to NaN
          return.
        else.
          HT=. 0 (< 1 , ,~ ilast) } HT
smoutput '(1) GOTO 50'
        end.
        label_50.
smoutput 'label 50'
        NB. T[ilast,ilast]=0 - clear H[ilast,ilast-1] to
        NB. split off a 1x1 block
        lios=. (>: ilast) th2lios ifrstm
        'HT cs'=. ((rot rotga) dbg 'rotga(H)') HT ; (< 0 ; lios ; (ilast - 0 1)) ; _1
        HT=. (< 1 ; (}: lios) ; (ilast - 0 1)) (((cs & rot) upd) dbg 'rot(T)') HT
        dZ=. dZ , cs , ilast - 0 1
      else.
        HT=. 0 (< 0 , ilast - 0 1) } HT
smoutput '(2) GOTO 60'
      end.
    else.
smoutput '(1) GOTO 60'
    end.
    label_60.
smoutput 'label 60'
    if. goto60 do.
      NB. H[ilast,ilast-1]=0 - standartize B, set alpha and
      NB. beta
      'HT signbc'=. hgeqzuxo (ilast , 1) ; HT  NB. process ilast-th eigenvalue (column)
      dZ=. dZ , 4 {."1 signbc , ilast
      NB. goto next block - exit if finished
      ilast=. <: ilast
      if. ilast < h do.
        NB. normal exit
        'HT signbc'=. hgeqzuxo (0 , 0 >. <: h) ; HT  NB. process eigenvalues (columns) 0:h-1
        dZ=. dZ , 4 {."1 signbc ,. i. 0 >. <: h
        HT ; dQ ; dZ
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
        NB. The Wilkinson shift, i.e., the eigenvalues of the
        NB. bottom-right 2x2 block of H*T^_1 which is nearest
        NB. to the bottom-right element.
        NB. We factor T as U*D, where U is unit upper
        NB. triangular, and compute (H*D^_1)*U^_1
        'U12 AD11 AD21 AD12 AD22'=. %/ (5 0 2 1 3 ,: 7 4 4 7 7) ({,) abscale * ((< a: ; ;~ ilast - 1 0) { HT)
        ABI22=. AD22 - U12 * AD21
        t1=. -: AD11 + ABI22
        rtdisc=. %: (t1 , AD12 , -AD11) mp (t1 , AD21 , AD22)
        temp=. +/ (*/) +. rtdisc , t1 - ABI22
        shift=. t1 - temp condneg rtdisc
      else.
        NB. Exceptional shift. Chosen for no paticularly good
        NB. reason
        eshift=. eshift + + %/ abscale * (;/ 0 1 ,. (_1 0 ,: _1 _1)+ ilast) { HT
        shift=. eshift
      end.
      NB. now check for two consecutive small subdiagonals
      HTd=. (0 _1 ,"0 1 ilast (] , -) ifirst) diag"1 2/ HT
      ctemp=. (- (shift&*))/ abscale * {. HTd
      temp=. (sorim }."1 ctemp) ,: ascale * sorim (< 1 ; 0 ; <<0) { HTd
      tempr=. >./ temp
      temp=. temp %"1 ((0 , 1 - FP_EPS) I. tempr) } 1 , tempr ,: 1
      'istart ctemp'=. (+&ifirst , {&ctemp) (ilast - ifirst) | >: (>:/ temp * atol ,: sorim (< 1 ; 0 ; <<_1) { HTd) i: 1
smoutput 'label 90'
      NB. do an implicit-shift QZ sweep
      NB. initial Q
      cs=. }: lartg ctemp , ascale * (< 0 , istart + 1 0) { HT
      NB. sweep
      j=. istart
      liosc=. (>: ilastm) th2lios <: j
      liosr=. (j + 2) th2lios ifrstm
      while. j < ilast do.
smoutput 'loop 150 J = ',":j
        lios=. j + 0 1
        NB. is a first iteration?
        if. j = istart do.
          HT=. (< a: ; lios ; (}. liosc)) (((cs & (rot &. |:)"2) upd) dbg 'rot(HT)') HT
        else.
          'HT cs'=. (((rot &. |:) rotga) dbg 'rotga(H)') HT ; (< 0 ; lios ; liosc) ; < < a: ; 0
          liosc=. }. liosc
          HT=. (< 1 ; lios ; liosc) (((cs & (rot &. |:)) upd) dbg 'rot(T)') HT
        end.
        dQ=. dQ , (+ cs) , lios
        lios=. j + 1 0
        'HT cs'=. ((rot rotga) dbg 'rotga(T)') HT ; (< 1 ; liosr ; lios) ; _1
        NB. isn't a last iteration?
        if. j < <: ilast do.
          liosr=. liosr , j + 2
        end.
        HT=. (< 0 ; liosr ; lios) (((cs & rot) upd) dbg 'rot(H)') HT
        dZ=. dZ , cs , lios
        j=. >: j
      end.
    end.
    jiter=. >: jiter
  end.
  NB. drop-through means non-convergence, set incorrect eigenvalues 0:ilast to NaN
  ((_. ; 0 0 , >: ilast) setdiag"2 HT) ; ,~ a:
)

NB. ---------------------------------------------------------
NB. hgeqzlenn
NB. hgeqzuenn
NB.
NB. Description:
NB.   Eigenvalues of pair of structured matrices
NB.
NB. Syntax:
NB.   'ab dQ dZ'=. hgeqzxenn hs ; HT
NB. where
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines submatrices H11 and T11 position in H
NB.           and T, respectively (see ggbalxp)
NB.   HT    -: H ,: T
NB.   ab    -: alpha ,: beta
NB.   dQ,dZ - either (i.0) if QZ iteration did not converge,
NB.           or any×4-matrix, accumulates scalings and
NB.           rotations to form Q and Z later, see rotsclx;
NB.           dQ and dZ may have the same shapes
NB.   alpha - n-vector, defines eigenvalues
NB.   beta  - n-vector, defines eigenvalues
NB.   H     - n×n-matrix, either lower (hgeqzlenn) or upper
NB.           (hgeqzuenn) Hessenberg inside the submatrix
NB.           H[h:h+s-1,h:h+s-1], and lower (hgeqzlenn) or
NB.           upper (hgeqzuenn) triangular outside
NB.   T     - n×n-matrix, either lower (hgeqzlenn) or upper
NB.           (hgeqzuenn) triangular
NB.
NB. Notes:
NB. - hgeqzuenn implements LAPACK's xHGEQZ('E','N')
NB. - non-converged eigenvalues are set to NaN

hgeqzuenn=: 0 ((diag"2&.>) upd) (hgeqzueo`[`(2 1&{@,`[@.((<{:)~{.))`[ hgeqzunn)

NB. ---------------------------------------------------------
NB. hgeqzlsnn
NB. hgeqzusnn
NB.
NB. Description:
NB.   Reduce Hessenberg-triangular pair (H,T) to generalized
NB.   Schur form:
NB.     H = Q*S*Z**H
NB.     T = Q*P*Z**H
NB.
NB. Syntax:
NB.   'SP dQ dZ'=. hgeqzxsnn hs ; HT
NB. where
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines submatrices H11 and T11 position in H
NB.           and T, respectively (see ggbalxp)
NB.   HT    -: H ,: T
NB.   SP    -: S ,: P
NB.   dQ,dZ - any×4-matrix, accumulates scalings and
NB.           rotations to form Q and Z later, see rotsclx;
NB.           dQ and dZ may have the same shapes
NB.   H     - n×n-matrix, either lower (hgeqzlsnn) or upper
NB.           (hgeqzusnn) Hessenberg inside the submatrix
NB.           H[h:h+s-1,h:h+s-1], and lower (hgeqzlenn) or
NB.           upper (hgeqzuenn) triangular outside
NB.   T     - n×n-matrix, either lower (hgeqzlsnn) or upper
NB.           (hgeqzusnn) triangular
NB.   S     - n×n-matrix, , ...
NB.   P     - n×n-matrix, , ...
NB.
NB. Notes:
NB. - hgeqzusnn implements LAPACK's xHGEQZ('S','N')
NB. - non-converged eigenvalues are set to NaN
NB. - generalized eigenvalues are defined by diagonals:
NB.     alpha -: diag S
NB.     beta -: diag P

hgeqzusnn=:                      hgeqzuso`]`]                      `] hgeqzunn

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB.   ab=.                 hgeqzxe hs ; HT
NB.   'ab Q Z'=. (Q1 ; Z1) hgeqzxe hs ; HT

hgeqzle=: (0 {:: hgeqzlenn) : (({.@] , (rotscll &. > }.)) hgeqzlenn)
hgeqzue=: (0 {:: hgeqzuenn) : (({.@] , (rotsclu &. > }.)) hgeqzuenn)

NB. ---------------------------------------------------------
NB.   SP=.                 hgeqzxs hs ; HT
NB.   'SP Q Z'=. (Q1 ; Z1) hgeqzxs hs ; HT
NB.
NB. Assertions (with appropriate comparison tolerance):#############
NB.   Q -: ungql QfR
NB.   I -: (mp~ ct) Q
NB.   A -: Q mp R
NB.   (] -: ((         tru   @  }:   ) mp~ ungqr)@geqrf) A
NB. where
NB.   QfR=. geqrf A
NB.   Q=. ungqr QfR
NB.   R=. tru }: QfR
NB.   hgeqzusnn_mt_ 2 3 ; 0 {:: ggbalu_mt_ AB
NB.   ] 'Q1 R'=. ((tru_mt_@}:) ;~ ungqr_mt_) geqrf_mt_ (0;1) {:: ggbalu_mt_ AB
NB.   'HT Q Z'=. (Q1 ; idmat_mt_ 7) gghrdu_mt_ 2 3 ; R 1} 0 {:: ggbalu_mt_ AB
NB.   'SP dQ dZ'=. hgeqzxsnn hs ; HT

hgeqzls=: (0 {:: hgeqzlsnn) : (({.@] , (rotscll &. > }.)) hgeqzlsnn)
hgeqzus=: (0 {:: hgeqzusnn) : (({.@] , (rotsclu &. > }.)) hgeqzusnn)

NB. =========================================================
NB. Test suite
