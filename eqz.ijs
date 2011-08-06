NB. Eigenvalues and eigenvectors of pair of structured
NB. matrices
NB.
NB. hgeqzex  Eigenvalues and, optionally, eigenvectors of pair
NB.          of structured matrices
NB. hgeqzsx  Eigenvalues, the Schur form and, optionally,
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
NB. hgeqzelo
NB. hgeqzeuo
NB.
NB. Description:
NB.   Calculate generalized eigenvalues, which are correspond
NB.   to diagonal elements of H, T given
NB.
NB. Syntax:
NB.   'HTdupd signbc'=. hgeqzexo HTd
NB. where
NB.
NB. Notes:
NB. - unfortunately, the following produces rank error:
NB.     'signbc Td'=. mask } &. |: (1 ,. + * Td) ,: (0 ,. absb)

hgeqzeuo=: 3 : 0
  'Hd Td'=. y
  absb=. | Td
  mask=. FP_SFMIN < absb
  signbc=. mask } 1 ,. + * Td
  ((Hd * signbc) ,: (mask } 0 ,. absb)) ; signbc
)

NB. ---------------------------------------------------------
NB. hgeqzslo
NB. hgeqzsuo
NB.
NB. Description:
NB.   Calculate generalized eigenvalues of hs-segment and
NB.   reduce corresponding columns to generalized Schur form
NB.
NB. Syntax:
NB.   'HTupd signbc'=. hgeqzsxo hs ; HT
NB. where
NB.
NB. Notes:
NB. - hs applied here marks segment either before or after
NB.   orginal hs-segment, supplied to hgeqzxxxx

hgeqzsuo=: 3 : 0
  'hs HT'=. y
  lios=. dhs2lios hs
  'HTd signbc'=. hgeqzeuo (0 , hs) diag"2 HT
  subHT=. (HTd ;"1 a:) setdiag"1 2 lios {"1 HT
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
NB.   avoid lartg roundoff errors

rotga=: 1 : 0
  'A iossubA iosfg'=. y
  subA=. iossubA { A
  csr=. lartg iosfg { subA
  (((({: csr) , 0) iosfg } (}: csr) u subA) iossubA } A) ; (}: csr)
)

NB. ---------------------------------------------------------
NB. hgeqzelnn
NB. hgeqzeunn
NB.
NB. Description:
NB.   Eigenvalues of pair of structured matrices
NB.
NB. Syntax:
NB.   'ab dQ dZ'=. hgeqzexnn hs ; HT
NB. where
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines submatrices H11 and T11 position in H
NB.           and T, respectively (see ggbalxp)
NB.   HT    -: H ,: T
NB.   ab    -: alpha ,: beta
NB.   dQ,dZ - either (i.0) if QZ iteration did not converge,
NB.           or any×4-matrix, accumulates scalings and
NB.           rotations to form Q and Z later, see hgeqzxqz;
NB.           dQ and dZ may have the same shapes
NB.   alpha - n-vector, ...
NB.   beta  - n-vector, ...
NB.   H     - n×n-matrix, either lower (hgeqzelnn) or upper
NB.           (hgeqzeunn) Hessenberg inside the submatrix
NB.           H[h:h+s-1,h:h+s-1], and lower (hgeqzelnn) or
NB.           upper (hgeqzeunn) triangular outside
NB.   T     - n×n-matrix, either lower (hgeqzelnn) or upper
NB.           triangular
NB.
NB. Notes:
NB. - hgeqzeunn implements LAPACK's xHGEQZ('E','N')
NB. - non-converged eigenvalues are set to NaN

hgeqzeunn=: 3 : 0
  'hs HT'=. y
  e=. +/ 'h s'=. hs
  dQ=. dZ=. 4 0 $ 0
  abnorm=. (,.~ hs) norms"2 ;. 0 HT
  'atol btol'=. abtol=. FP_SFMIN >. FP_PREC * abnorm
  'ascale bscale'=. abscale=. % FP_SFMIN >. abnorm
  'HTd signbc'=. hgeqzeuo (0 , e) diag"2 HT  NB. calculate eigenvalues h+s:n-1
  HT=. (HTd (;"1) 0 , e) setdiag"1 2 HT
  dZ=. dZ , 4 {."1 signbc ,. (c HT) ht2lios e

  NB. main QZ iteration loop

  NB. initialize dynamic indices
  NB.
  NB. Eigenvalues h+s:n-1 have been found.
  NB. Column operations modify rows ifirstm:*
  NB. Column operations modify columns *:ilastm
  NB. ifirstm - the row of the last splitting row above row
  NB.           ilast, this is always at least h
  NB. iiter   - counts iterations since the last eigenvalue
  NB.           was found, to tell when to use an
  NB.           extraordinary shift
  NB. maxit   - the maximum number of QZ sweep allowed

  ilast=. <: e
  ifirstm=. h
  ilastm=. ilast
  iiter=. 0
  eshift=. 0
  maxit=. 30 * s
  jiter=. 0

  while. jiter < maxit do.
    goto60=. 1  NB. mark default branching
    NB. split the matrix if possible, by to tests:
    NB.   1. H[j,j-1]=0 OR j=h
    NB.   2. T[j,j]=0
    if. ilast ~: h do.
      if. atol < sorim (< 0 , ilast - 0 1) { HT do.
        if. btol < | (< 1 , ,~ ilast) { HT do.
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
            if. btol > | (< 1 , ,~ j) { HT do.
              HT=. 0 (< 1 , ,~ j) } HT
              NB. test 2a: check for 2 consecutive small subdiagonals in H
              ilazr2=. 0
              if. -. ilazro do.
                'Hjj1 Hj1j Hjj'=. ((<"1) 0 ,. (0 _1,1 0,:0 0) + j) { HT
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
                while. jch < ilast do.
                  'HT cs'=. (rot &. |:) rotga HT ; (< 0 ; (jch + 0 1) ; ((>: ilastm) ht2lios jch)) ; (< a: ; 0)
                  HT=. (< 1 ; (jch + 0 1) ; (ilastm ht2lios & >: jch)) (cs & (rot &. |:)) upd HT
                  dQ=. dQ , (+ cs) , jch + 0 1
                  if. ilazr2 do.
                    HT=. (< 0 , jch - 0 1) (* & ({. cs)) upd HT
                    ilazr2=. 0
                  end.
                  if. btol <: sorim (< 1 , jch + 1 1) { HT do.
                    if. ilast > >: jch do.
                      ifirst=. >: jch
                      goto60=. 0
                    end.
                    goto_60.
                  end.
                  HT=. 0 (< 1 , jch + 1 1) } HT
                  jch=. >: jch
                end.
              else.
                NB. Only test 2 passed - chase the zero to
                NB. T[ilast,ilast], then process as in the
                NB. case T[ilast,ilast]=0
                jch=. j
                while. jch < ilast do.
                  'HT cs'=. (rot &. |:) rotga HT ; (< 1 ; (jch + 0 1) ; (ilastm ht2lios & >: jch)) ; (< a: ; 0)
                  HT=. (< 0 ; (jch + 0 1) ; ((>: ilastm) ht2lios <: jch)) (cs & (rot &. |:)) upd HT
                  dQ=. dQ , (+ cs) , jch + 0 1
                  'HT cs'=. rot rotga HT ; (< 0 ; ((2 + jch) ht2lios ifirstm) ; (jch - 0 1)) ; (< _1 ; _1 0)
                  HT=. (< 1 ; (jch ht2lios ifirstm) ; (jch - 0 1)) (cs & rot) upd HT
                  dZ=. dZ , cs , jch - 0 1
                  jch=. >: jch
                end.
              end.
              goto_50.
            elseif. ilazro do.
              ifirst=. j
              goto60=. 0
              goto_60.
            end.
            j=. <: j
          end.
          NB. drop-through is impossible
          ((2 , c HT) $ _.) ; ,~ a:  NB. set all eigenvalues to NaN
          return.
        else.
          HT=. 0 (< 1 , ,~ ilast) } HT
        end.
        label_50.
        NB. T[ilast,ilast]=0 - clear H[ilast,ilast-1] to
        NB. split off a 1x1 block
        'HT cs'=. rot rotga HT ; (< 0 ; ((>: ilast) ht2lios ifirstm) ; (ilast - 0 1)) ; (< _1 ; 1 0)
        HT=. (< 1 ; (ilast ht2lios ifirstm) ; (ilast - 0 1)) (cs & rot) upd HT
        dZ=. dZ , cs , ilast - 0 1
      else.
        HT=. 0 (< 0 , ilast - 0 1) } HT
      end.
    end.
    label_60.
    if. goto60 do.
      NB. H[ilast,ilast-1]=0 - standartize B, set alpha and
      NB. beta
      'HTd signbc'=. hgeqzeuo (0 , ilast , 1) diag"2 HT  NB. ilast-th eigenvalue
      HT=. (, HTd) (< a: ; ;~ ilast) } HT
      dZ=. dZ , 4 {."1 signbc , ilast
      NB. goto next block - exit if finished
      ilast=. <: ilast
      if. ilast < h do.
        NB. normal exit
        'HTd signbc'=. hgeqzeuo (0 0 , 0 >. <: h) diag"2 HT  NB. calc eigenvalues 0:h-1
        HTd=. HTd ,. (0 , h) diag"1 2 HT                     NB. augment by eigenvalues h:n-1
        dZ=. dZ , 4 {."1 signbc ,. i. 0 >. <: h
        HTd ; dQ ; dZ
        return.
      end.
      NB. reset counters
      iiter=. 0
      eshift=. 0
      ilastm=. ilast
      if. ifirstm > ilast do.
        ifirstm=. h
      end.
    else.
      NB. QZ step
      NB. This iteration only involves rows/columns
      NB. ifirst:ilast. We assume ifirst<ilast, and that the
      NB. diagonal of B is non-zero
      iiter=. >: iiter
      ifirstm=. ifirst
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
      HTd=. (0 _1 ,"0 1 ilast (] , -) ifirst) diag"1 2 HT
      ctemp=. (- (shift&*))/ abscale * {. HTd
      temp=. (sorim }."1 ctemp) ,: ascale * sorim (< 1 ; 0 ; <<0) { HTd
      tempr=. >./ temp
      temp=. temp % ((0 , 1 - FP_EPS) I. tempr) } 1 , tempr ,: 1
      flags=. >:/ temp * atol ,: sorim (< 1 ; 0 ; <<_1) { HTd
      io=. flags i: 1
      if. io < # flags do.
        istart=. >: io + ifirst
        ctemp=. (>: io) { ctemp
      else.
        istart=. ifirst
        ctemp=. {. ctemp
      end.
      NB. do an implicit-shift QZ sweep
      if. istart < ilast do.
        NB. initial Q
        cs=. }: clartg ctemp , ascale * (< 0 , istart + 1 0) { HT
        NB. sweep
        j=. istart
        whilest.
          ios=. < 0 ; 0 1 (+ ; (<:@])) j
          csr=. lartg ios { HT
          HT=. (({: csr) , 0) ios } HT
          cs=. }: csr
          j < ilast
        do.
          ios=. < a: ; (j + 0 1) ; ((>: ilastm) ht2lios j)
          HT=. ios (cs & (rot &. |:)"2) upd HT
          dQ=. dQ , cs , j + 0 1
          ios=. < 1 ([ ; + ; + , ]) j
          cs=. }: csr=. lartg ios { HT
          HT=. (({: csr) , 0) ios } HT
          ios=. < 0 ; (((>: ilastm) <. j + 2) ht2lios ifirstm) ; (j + 1 0)
          HT=. ios (cs & rot) upd HT
          ios=. < 1 ; ((>: j) ht2lios ifirstm) ; (j + 1 0)
          HT=. ios (cs & rot) upd HT
          dZ=. dZ , cs , j + 1 0
          j=. >: j
        end.
      end.
    end.
    jiter=. >: jiter
  end.
  NB. drop-through means non-convergence, set incorrect eigenvalues[0:ilast] to NaN
  (_. (< a: ; i. >: ilast) } diag"2 HT) ; ,~ a:
)

NB. ---------------------------------------------------------
NB. hgeqzslnn
NB. hgeqzsunn
NB.
NB. Description:
NB.   Reduce Hessenberg-triangular pair (H,T) to generalized
NB.   Schur form:
NB.     H = Q*S*Z**H
NB.     T = Q*P*Z**H
NB.
NB. Syntax:
NB.   'SP dQ dZ'=. hgeqzexnn hs ; HT
NB. where
NB.   hs    - 2-vector of integers (h,s) 'head' and 'size',
NB.           defines submatrices H11 and T11 position in H
NB.           and T, respectively (see ggbalxp)
NB.   HT    -: H ,: T
NB.   SP    -: S ,: P
NB.   dQ,dZ - any×4-matrix, accumulates scalings and
NB.           rotations to form Q and Z later, see hgeqzxqz;
NB.           dQ and dZ may have the same shapes
NB.   H     - n×n-matrix, either lower (hgeqzelnn) or upper
NB.           (hgeqzeunn) Hessenberg inside the submatrix
NB.           H[h:h+s-1,h:h+s-1], and lower (hgeqzelnn) or
NB.           upper (hgeqzeunn) triangular outside
NB.   T     - n×n-matrix, either lower (hgeqzelnn) or upper
NB.           triangular
NB.   S     - n×n-matrix, , ...
NB.   P     - n×n-matrix, , ...
NB.
NB. Notes:
NB. - hgeqzsunn implements LAPACK's xHGEQZ('S','N')
NB.
NB. TODO:
NB. - 

hgeqzsunn=: [:

NB. ---------------------------------------------------------
NB. hgeqzuqz
NB.
NB. Description:
NB.   Apply scalings and rotations accumulated in dQ (dZ) to
NB.   Q1 (Z1) to produce Q (Z) explicitely
NB.
NB. Syntax:
NB.   'a Q Z'=. (Q1 ; Z1) gghrduqz a ; dQ ; dZ
NB. where
NB.   a       - pass-through argument
NB.   dQ,dZ   - any×4-matrix, where each row is 4-vector of
NB.             values, either:
NB.               m , io , 0 , 0
NB.             or:
NB.               c , s , iof , iog
NB.             accumulates scalings and rotations to form Q
NB.             and Z; dQ and dZ may have the same shapes
NB.   Q1,Z1   - n×n-matrix or (i.0), the unitary
NB.             (orthogonal), typically from the gghrdx
NB.             reduction of matrix pair (A,B)
NB.   Q       - either (i.0) when Q1 -: (i.0) , or n×n-matrix
NB.             (Q1*ΔQ) otherwise
NB.   Z       - either (i.0) when Z1 -: (i.0) , or n×n-matrix
NB.             (Z1*ΔZ) otherwise
NB.   m       - scalar to scale column
NB.   io      - lIO column to scale
NB.   (c,s)   - 2-vector, defines rotation for 2-vectors
NB.             (f,g), being a curtailed output of larfg
NB.   iof,iog - lIOS columns to form 2-vectors (f,g) to
NB.             rotate

hgeqzuqz=: 4 : 0
  'Q1 Z1'=. x
  'a dQ dZ'=. y
  if. Q1 -.@-: (i. 0) do.
    #############
    rQZ=. (< a: ; 0 ; 1) (+ upd) rQZ     NB. conjugate all s of Q's part of rQZ
    j=. h
    k=. 0
    while. j < e do.                     NB. (s-1)-vector: h,h+1,h+2,...,h+s-2
      i=. e
      while. i > (j+1) do.               NB. (h+s-j-2)-vector: h+s-1,h+s-2,...,j+2
        iospair=. < a: ; (i - 1 0)
        q=. (< k , 0) { rQZ              NB. extract current rotation q[k]
        Q1=. iospair ((q & rot) upd) Q1  NB. update columns by q[k]
        k=. >: k
        i=. <: i
      end.
      j=. >: j
    end.
  end.
  if. Z1 -.@-: (i. 0) do.
    j=. h
    k=. 0
    while. j < e do.                     NB. (s-1)-vector: h,h+1,h+2,...,h+s-2
      i=. e
      while. i > (j+1) do.               NB. (h+s-j-2)-vector: h+s-1,h+s-2,...,j+2
        iospair=. < a: ; (i - 0 1)
        z=. (< k , 1) { rQZ              NB. extract current rotation z[k]
        Z1=. iospair ((z & rot) upd) Z1  NB. update columns by z[k]
        k=. >: k
        i=. <: i
      end.
      j=. >: j
    end.
  end.
  a ; Q1 ; Z1
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB.   ab=.                 hgeqzex hs ; HT
NB.   'ab Q Z'=. (Q1 ; Z1) hgeqzex hs ; HT

hgeqzel=: (0 {:: hgeqzelnn) : (hgeqzelqz hgeqzelnn)
hgeqzeu=: (0 {:: hgeqzeunn) : (hgeqzeuqz hgeqzeunn)

NB. ---------------------------------------------------------
NB.   SP=.                 hgeqzsx hs ; HT
NB.   'SP Q Z'=. (Q1 ; Z1) hgeqzsx hs ; HT

hgeqzsl=: (0 {:: hgeqzslnn) : (hgeqzslqz hgeqzslnn)
hgeqzsu=: (0 {:: hgeqzsunn) : (hgeqzsuqz hgeqzsunn)

NB. =========================================================
NB. Test suite
