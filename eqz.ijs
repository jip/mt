NB. Eigenvalues and eigenvectors of pair of structured
NB. matrices
NB.
NB. hseqzex  Eigenvalues and, optionally, eigenvectors of pair
NB.          of structured matrices
NB. hseqzsx  Eigenvalues, the Schur form and, optionally,
NB.          eigenvectors of pair of structured matrices
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
NB. hseqzelo
NB. hseqzeuo
NB.
NB. Description:
NB.   Eigenvalues of pair of structured matrices outside of
NB.   hs-segment, is called from hseqzexnn
NB.
NB. Syntax:
NB.   'HTdupd signbc'=. hseqzexo HTd
NB. where
NB.
NB. Notes:
NB. - unfortunately, the following produces rank error:
NB.     'signbc Td'=. mask } &. |: (1 ,. + * Td) ,: (0 ,. absb)

hseqzeuo=: 3 : 0
  'Hd Td'=. y
  absb=. | Td
  mask=. FP_SFMIN < absb
  signbc=. mask } 1 ,. + * Td
  ((Hd * signbc) ,: (mask } 0 ,. absb)) ; signbc
)

NB. ---------------------------------------------------------
NB. hseqzslo
NB. hseqzsuo
NB.
NB. Description:
NB.   Eigenvalues of pair of structured matrices outside of
NB.   hs-segment, is called from hseqzsxnn
NB.
NB. Syntax:
NB.   'HTupd signbc'=. hseqzsxo hs ; HT
NB. where

hseqzsuo=: 3 : 0
  'hs HT'=. y
  lios=. dhs2lios hs
  'HTd signbc'=. hseqzeuo (0 , hs) diag"2 HT
  subHT=. (HTd ;"1 0 a:) setdiag"1 2 lios {"1 HT
  (((i. c HT) </ lios } subHT ,: subHT *"1 signbc) lios }"1 HT) ; signbc
)

NB. ---------------------------------------------------------
NB. rotga
NB.
NB. Description:
NB.   Adv. to make verb to generate and apply rotation
NB.
NB. Syntax:
NB.   vapp=. vrota rotga
NB. where
NB.   vrota   - dyad to apply rotation; is called as:
NB.               subAupd=. cs vrota subA
NB.             and is any of:
NB.               rot        NB. apply rotation to rows
NB.               rot &. |:  NB. apply rotation to columns
NB.   vapp    - monad to generate and apply rotation; is
NB.             called as:
NB.               'Aupd cs'=. vapp A ; iossubA ; iosfg
NB.   cs      - 2-vector (c,s), curtailed output of lartg,
NB.             defines rotation matrix
NB.   A       - n×n-matrix to update
NB.   Aupd    - n×n-matrix, updated A, being A with subA
NB.             replaced by subAupd
NB.   subA    - 2×m-matrix or m×2-matrix, array of 2-vectors
NB.             to apply rotation
NB.   subAupd - matrix of the same shape as subA, the rotated
NB.             subA
NB.   iossubA - ios for subA (subAupd) within A (Aupd)
NB.   iosfg   - ios within subA of 2-vector (f,g) to generate
NB.             rotation
NB.
NB. Notes:
NB. - rotated 2-vector (r,0) is written into A explicitely to
NB.   avoid rotation roundoff errors

rotga=: 1 : 0
  'A iossubA iosfg'=. y
  subA=. iossubA { A
  csr=. lartg iosfg { subA
  (((({: csr) , 0) iosfg } (}: csr) u subA) iossubA } A) ; (}: csr)
)

NB. ---------------------------------------------------------
NB. hseqzelnn
NB. hseqzeunn
NB.
NB. Description:
NB.   Eigenvalues of pair of structured matrices
NB.
NB. Syntax:
NB.   'hs ab dQZ'=. hseqzexnn hs ; HT
NB. where
NB.   HT    -: H ,: T
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines submatrices H11 and T11 position in H
NB.           and T, respectively (see ggbalxp)
NB.   ab    -: alpha ,: beta
NB.   alpha - n-vector, ...
NB.   beta  - n-vector, ...
NB.   H     - n×n-matrix, either lower (hseqzelnn) or upper
NB.           (hseqzeunn) Hessenberg inside the submatrix
NB.           H[h:h+s-1,h:h+s-1], and lower (hseqzelnn) or
NB.           upper (hseqzeunn)  triangular outside
NB.   T     - n×n-matrix, either lower (hseqzelnn) or upper
NB.           triangular
NB.   dQZ   - m×2×3-report, where each 2×3-matrix
NB.           is 2-vector of laminated rotations:
NB.             ( cs(Q)[i,j]  sn(Q)[i,j]  k0)
NB.             ( cs(Z)[i,j]  sn(Z)[i,j]  k1)
NB.           accumulates rotations to form Q and Z later
NB.
NB. Notes:
NB. - hseqzeunn implements LAPACK's xHGEQZ('E','N')
NB. - hseqzsunn implements LAPACK's xHGEQZ('S','N')
NB.
NB. TODO:
NB. - init ab by NaN

hseqzeunn=: 3 : 0
  'hs HT'=. y
  'h s'=. hs
  n=. c HT
  abnorm=. (,.~ hs) norms"2 ;. 0 HT
  'atol btol'=. abtol=. FP_SFMIN >. FP_PREC * abnorm
  'ascale bscale'=. abscale=. % FP_SFMIN >. abnorm
  'HTd signbc'=. hseqzeuo (0 , (+/ hs)) diag"2 HT  NB. eigenvalues[h+s:n-1]
  HT=. (HTd (;"1) 0 , +/ hs) setdiag"1 2 HT

  NB. main QZ iteration loop

  NB. initialize dynamic indices
  ilast=. h+s-1
  ifirstm=. h
  ilastm=. ilast
  iiter=. 0
  eshift=. 0
  maxit=. 30 * s
  jiter=. 0

  while. jiter < maxit do.
    NB. split the matrix if possible, by to tests:
    NB. 1. H[j,j-1]=0 OR j=h
    NB. 2. T[j,j]=0
    if. ilast ~: h do.
      if. atol < sorim (< 0 , ilast - 0 1) { HT do.
        if. btol < | (< 1 , 2 # ilast) { HT do.
          NB. general case: j < ilast
          j=. <: ilast
          while. j >: h do.
            NB. test 1: H[j,j-1]=0 OR j=h
            if. j = h do.
              ilazro=. 1
            elseif. atol >: sorim (< 0 , j - 0 1) { HT do.
              HT=. 0 (< 0 , j - 0 1) } HT
              ilazro=. 1
            elseif. do.
              ilazro=. 0
            end.
            NB. test 2: T[j,j]=0
            if. btol > | (< 1 , 2 # j) { HT do.
              HT=. 0 (< 1 , 2 # j) } HT
              NB. test 2a: check for 2 consecutive smallsubdiagonals in A
              ilazr2=. 0
              if. -. ilazro do.
                'Hjj1 Hj1j Hjj'=. (<"1 (0 ,. (0 _1,1 0,:0 0) + j)) { HT
                if. 0 >: (Hjj1 , Hjj) mp (ascale * (Hj1j , -atol)) do.
                  ilazr2=. 1
                end.
              end.
              NB. if both tests (1 & 2) pass, then split a 1x1 block off at the top.
              NB. ...
              if. ilazro +. ilazr2 do.
                jch=. j
                while. jch < ilast do.
                  'HT cs'=. (rot &. |:) rotga HT ; (< 0 ; (jch + 0 1) ; (ilastm ht2lios jch)) ; (< a: ; 0)
                  HT=. (< 1 ; (jch + 0 1) ; (ilastm ht2lios >: jch)) (cs & (rot &. |:)) upd HT
                  dQ=. dQ , cs , jch
                  if. ilazr2 do.
                    HT=. (< 0 , jch - 0 1) (* & ({. cs)) upd HT
                  end.
                  ilazr2=. 0
                  if. btol <: sorim (< 1 , jch + 1 1) { HT do.
                    if. ilast <: >: jch do.
                      goto_60.
                    else.
                      goto_70.
                    end.
                  end.
                  HT=. 0 (< 1 , jch + 1 1) } HT
                  jch=. >: jch
                end.
              else.
                NB. only test 2 passed - chase the zero to T[ilast,ilast],
                NB. then process as in the case T[ilast,ilast]=0
                jch=. j
                while. jch < ilast do.
                  'HT cs'=. (rot &. |:) rotga HT ; (< 1 ; (jch + 0 1) ; (ilastm ht2lios >: jch)) ; (< a: ; 0)
                  HT=. (< 0 ; (jch + 0 1) ; (ilastm ht2lios <: jch)) (cs & (rot &. |:)) upd HT
                  dQ=. dQ , cs , jch
                  'HT cs'=. rot rotga HT ; (< 0 ; ((>: jch) ht2lios ifirstm) ; (jch - 0 1)) ; (< _1 ; _1 0)
                  HT=. (< 1 ; ((<: jch) ht2lios ifirstm) ; (jch - 0 1)) (cs & rot) upd HT
                  dZ=. dZ , cs , jch
                  jch=. >: jch
                end.
              end.
              goto_50.
            elseif. ilazro do.
              ifirst=. j
              goto_70.
            end.
            j=. <: j
          end.
          NB. drop-through is impossible
          hs ; (_. (< a: ; i. >: ilast) } diag"2 HT) ; dQZ  NB. set incorrect eigenvalues[0:ilast] to NaN
          return.
        else.
          HT=. 0 (< 1 , 2 # ilast) } HT
        end.
        label_50.
        NB. T[ilast,ilast]=0 - clear H[ilast,ilast-1] to split off a 1x1 block
        'HT cs'=. rot rotga HT ; (< 0 ; (ilast ht2lios ifirstm) ; (ilast - 0 1)) ; (< _1 ; _1 0)
        HT=. (< 1 ; ((<: ilast) ht2lios ifirstm) ; (ilast - 0 1)) (cs & rot) upd HT
        dZ=. dZ , cs , ilast
      else.
        HT=. 0 (< 0 , ilast - 0 1) } HT
      end.
    end.
    label_60.
    NB. H[ilast,ilast-1]=0 - standartize B, set alpha and beta
    'HTd signbc'=. hseqzeuo (0 , ilast , 1) diag"2 HT  NB. eigenvalues[ilast]
    HT=. (, HTd) (< a: ; ;~ ilast) } HT
    NB. goto next block - exit if finished
    ilast=. <: ilast
    if. ilast < h do.
      NB. normal exit
      'HTd signbc'=. hseqzeuo (0 0 , 0 >. <: h) diag"2 HT  NB. calc eigenvalues[0:h-1]
      HTd=. HTd ,. (0 , h) diag"1 2 HT                     NB. read eigenvalues[h:]
      hs ; HTd ; dQZ
      return.
    end.
    NB. reset counters
    iiter=. 0
    eshift=. 0
    ilastm=. ilast
    if. ifirstm > ilast do.
      ifirstm=. h
    end.
    goto_160.
    NB. QZ step
    NB. This iteration only involves rows/columns
    NB. ifirst:ilast. We assume ifirst<ilast, and that the
    NB. diagonal of B is non-zero
    label_70.
    iiter=. >: iiter
    ifirstm=. ifirst
    NB. compute the shift
    NB. at this point, ifirst<ilast, and the diagonal
    NB. elements of T[ifirst:ilast,ifirst:ilast] are larger
    NB. than btol in magnitude
    if. 0 -. 10 | iiter do.
      NB. The Wilkinson shift, i.e., the eigenvalues of the
      NB. bottom-right 2x2 block of A*B^_1 which is nearest
      NB. to the bottom-right element.
      NB. We factor B as U*D, where U is unit upper
      NB. triangular, and compute (A*D^_1)*U^_1
      'U12 AD11 AD21 AD12 AD22'=. %/ (5 0 2 1 3 ,: 7 4 4 7 7) ({,) abscale * ((< a: ; ;~ ilast - 1 0) { HT)
      ABI22=. AD22 - U12 * AD21
      t1=. -: AD11 + ABI22
      rtdisc=. %: (t1 , AD12 , -AD11) mp (t1 , AD21 , AD22)
      temp=. +/ (*/) +. rtdisc , t1 - ABI22
      shift=. t1 - temp condneg rtdisc
    else.
      NB. Exceptional shift. Chosen for no paticularly good reason
      eshift=. eshift + + %/ abscale * (< a: ; (;/) ilast - 1 0) { HT
      shift=. eshift
    end.
    NB. now check for two consecutive small subdiagonals
    NB. ????????????????????
    HTd=. (0 , ilast (] , -) >: ishift) diag"2 HT
    Hd=. (_1 , ifirst , (ilast - (ifirst + 2))) diag {. HT
    temp=. sorim (- (shift&*))/ abscale * HTd
    temp2=. ascale * sorim }. Hd
    tempr=. temp > temp2
    'temp temp2'=. (temp ,: temp2) % ((0 , 1 - FP_EPS) I. tempr) } 1 , tempr , 1
    flags=. >:/ (temp ,: temp2) * (atol ,: (sorim }: Hd))
    io=. flags i: 1
    if. io < # flags do.
      istart=. j
      goto_90.
    end.
    NB. j=. <: ilast
    NB. while. j > ifirst do.
    NB.   ctemp=. (1 , -shift) mp abscale * (< a: ; ;~ j) { HT
    NB.   temp=. (sorim ctemp) , (ascale * sorim (< 0 , j + 1 0) { HT)
    NB.   tempr=. >./ temp
    NB.   'temp=. temp % 
    NB.   j=. >: j
    NB. end.
    istart=. ifirst
    ctemp=. (1 , -shift) mp abscale * (< a: ; ;~ ifirst) { HT
    label_90.
    NB. do an implicit-shift QZ sweep
    if. istart < ilast do.
      NB. initial Q
      cs=. }: clartg ctemp , ascale * (< 0 , istart + 1 0) { HT
      NB. sweep
      j=. istart
      whilest.
        ios=. < 0 ; 0 1 (+ ; ]) j
        csr=. lartg ios { HT
        HT=. (({: csr) , 0) ios } HT
        cs=. }: csr
        j < ilast
      do.
        ios=. < a: ; (j + 0 1) ; (ilastm ht2lios j)
        HT=. ios (cs & (rot &. |:)"2) upd HT
        dQ=. dQ , (cs , j + 0 1)  NB. FIXME!
        ios=. < 1 ; j ; (j + 1 0)
        csr=. lartg ios { HT
        HT=. (({: csr) , 0) ios } HT
        ios=. < 0 ; ((ilastm <. j + 2) ht2lios ifirstm) ; (j + 1 0)
        HT=. ios (cs & rot) upd HT
        ios=. < 1 ; (j ht2lios ifirstm) ; (j + 1 0)
        HT=. ios (cs & rot) upd HT
        dZ=. dZ , (cs , j + 1 0)  NB. FIXME!
        j=. >: j
      end.
    end.
    label_160.
    jiter=. >: jiter
  end.
  NB. drop-through means non-convergence, set incorrect eigenvalues[0:ilast] to NaN
  hs ; (_. (< a: ; i. >: ilast) } diag"2 HT) ; dQZ
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB.   ab=.                 hseqzex hs ; HT
NB.   'ab Q Z'=. (Q1 ; Z1) hseqzex hs ; HT

hseqzel=: (1 {:: hseqzelnn) : (hseqzelqz hseqzelnn)
hseqzeu=: (1 {:: hseqzeunn) : (hseqzeuqz hseqzeunn)

NB. ---------------------------------------------------------
NB.   SP=.                 hseqzsx hs ; HT
NB.   'SP Q Z'=. (Q1 ; Z1) hseqzsx hs ; HT

hseqzsl=: (1 {:: hseqzslnn) : (hseqzslqz hseqzslnn)
hseqzsu=: (1 {:: hseqzsunn) : (hseqzsuqz hseqzsunn)

NB. =========================================================
NB. Test suite
