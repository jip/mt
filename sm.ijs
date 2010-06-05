NB. sm.ijs
NB. Solve matrix linear monomial equation with triangular
NB. matrix
NB.
NB. trsmlx    Solve equation L * X = B, where L is a lower
NB.           triangular matrix
NB. trsml1x   Solve equation L1 * X = B, where L1 is an unit
NB.           lower triangular matrix
NB. trsml1hx  Solve equation L1^H * X = B, where L1 is an
NB.           unit lower triangular matrix
NB. trsml1tx  Solve equation L1^T * X = B, where L1 is an
NB.           unit lower triangular matrix
NB. trsmlhx   Solve equation L^H * X = B, where L is a lower
NB.           triangular matrix
NB. trsmltx   Solve equation L^T * X = B, where L is a lower
NB.           triangular matrix
NB. trsmux    Solve equation U * X = B, where U is an upper
NB.           triangular matrix
NB. trsmu1hx  Solve equation U1^H * X = B, where U1 is an
NB.           unit upper triangular matrix
NB. trsmu1tx  Solve equation U1^T * X = B, where U1 is an
NB.           unit upper triangular matrix
NB. trsmuhx   Solve equation U^H * X = B, where U is an
NB.           upper triangular matrix
NB. trsmutx   Solve equation U^T * X = B, where U is an
NB.           upper triangular matrix
NB. trsmu1x   Solve equation U1 * X = B, where U1 is an unit
NB.           upper triangular matrix
NB. trsmxl    Solve equation X * L = B, where L is a lower
NB.           triangular matrix
NB. trsmxl1   Solve equation X * L1 = B, where L1 is an unit
NB.           lower triangular matrix
NB. trsmxl1h  Solve equation X * L1^H = B, where L1 is an
NB.           unit lower triangular matrix
NB. trsmxl1t  Solve equation X * L1^T = B, where L1 is an
NB.           unit lower triangular matrix
NB. trsmxlh   Solve equation X * L^H = B, where L is a lower
NB.           triangular matrix
NB. trsmxlt   Solve equation X * L^T = B, where L is a lower
NB.           triangular matrix
NB. trsmxu    Solve equation X * U = B, where U is an upper
NB.           triangular matrix
NB. trsmxu1   Solve equation X * U1 = B, where U1 is an unit
NB.           upper triangular matrix
NB. trsmxu1h  Solve equation X * U1^H = B, where U1 is an
NB.           unit upper triangular matrix
NB. trsmxu1t  Solve equation X * U1^T = B, where U1 is an
NB.           unit upper triangular matrix
NB. trsmxuh   Solve equation X * U^H = B, where U is an
NB.           upper triangular matrix
NB. trsmxut   Solve equation X * U^T = B, where U is an
NB.           upper triangular matrix
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb:          Solves:           Syntax:
NB. trsmlx         L    * X = B      Xv=. A trsmlx   Bv
NB. trsml1x        L1   * X = B      Xv=. A trsml1x  Bv
NB. trsml1hx       L1^H * X = B      Xv=. A trsml1hx Bv
NB. trsml1tx       L1^T * X = B      Xv=. A trsml1tx Bv
NB. trsmlhx        L^H  * X = B      Xv=. A trsmlhx  Bv
NB. trsmltx        L^T  * X = B      Xv=. A trsmltx  Bv
NB. trsmux         U    * X = B      Xv=. A trsmux   Bv
NB. trsmu1x        U1   * X = B      Xv=. A trsmu1x  Bv
NB. trsmu1hx       U1^H * X = B      Xv=. A trsmu1hx Bv
NB. trsmu1tx       U1^T * X = B      Xv=. A trsmu1tx Bv
NB. trsmuhx        U^H  * X = B      Xv=. A trsmuhx  Bv
NB. trsmutx        U^T  * X = B      Xv=. A trsmutx  Bv
NB. trsmxl         X * L    = B      Xh=. A trsmxl   Bh
NB. trsmxl1        X * L1   = B      Xh=. A trsmxl1  Bh
NB. trsmxl1h       X * L1^H = B      Xh=. A trsmxl1h Bh
NB. trsmxl1t       X * L1^T = B      Xh=. A trsmxl1t Bh
NB. trsmxlh        X * L^H  = B      Xh=. A trsmxlh  Bh
NB. trsmxlt        X * L^T  = B      Xh=. A trsmxlt  Bh
NB. trsmxu         X * U    = B      Xh=. A trsmxu   Bh
NB. trsmxu1        X * U1   = B      Xh=. A trsmxu1  Bh
NB. trsmxu1h       X * U1^H = B      Xh=. A trsmxu1h Bh
NB. trsmxu1t       X * U1^T = B      Xh=. A trsmxu1t Bh
NB. trsmxuh        X * U^H  = B      Xh=. A trsmxuh  Bh
NB. trsmxut        X * U^T  = B      Xh=. A trsmxut  Bh
NB.
NB. Description:
NB.   Solve linear monomial equation with triangular matrix
NB. where:
NB.   A    - n×n-matrix, containing either U, U1, L or L1
NB.   U    - n×n-matrix, upper triangular
NB.   U1   - n×n-matrix, unit upper triangular (diagonal is
NB.          not saved)
NB.   L    - n×n-matrix, lower triangular
NB.   L1   - n×n-matrix, unit lower triangular (diagonal is
NB.          not saved)
NB.   Bv   - n-vector or n×nrhs-matrix, the RHS
NB.   Bh   - n-vector or nrhs×n-matrix, the RHS
NB.   Xv   - same shape as Bv, the solution
NB.   Xh   - same shape as Bh, the solution
NB.   nrhs ≥ 0
NB.
NB. Notes:
NB. - opposite triangle is not referenced
NB. - unit diagonal is not referenced
NB. - emulate BLAS's xTRSV for vector RHS
NB. - emulate BLAS's xTRSM for matrix RHS with following
NB.   discrepancy: alpha parameter is not supported
NB.   (alpha≡1); to involve alpha, use the following pattern:
NB.     X=. A trsmxxxx alpha*B

trsmlx=:   ((((    (#@])  {   (1 {:: [)) - (] mp~ (((#@]) (1 dhs2lios (    * ,[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))) % (((#@])     (* >:)  (# @ (0 {:: [)))   ({,) (0 {:: [))) , ~ ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(li)=(1 dhs2lios (   i*n,i)), lio(lii)=   i*(n+1)
trsml1x=:  ( ((    (#@])  {   (1 {:: [)) - (] mp~ (((#@]) (1 dhs2lios (    * ,[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))                                                          ) , ~ ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(li)=(1 dhs2lios (   i*n,i))
trsml1hx=: ( (((_1-(#@])) {   (1 {:: [)) - (] mp~ (((#@]) (] dhs2lios ((_1-[),[)) (# @ (0 {:: [))) +@({,) (0 {:: [)))                                                          ) ,   ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i))
trsml1tx=: ( (((_1-(#@])) {   (1 {:: [)) - (] mp~ (((#@]) (] dhs2lios ((_1-[),[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))                                                          ) ,   ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i))
trsmlhx=:  (((((_1-(#@])) {   (1 {:: [)) - (] mp~ (((#@]) (] dhs2lios ((_1-[),[)) (# @ (0 {:: [))) +@({,) (0 {:: [)))) % (((#@]) (_1-(* >:)) (# @ (0 {:: [))) +@({,) (0 {:: [))) ,   ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i)), lio(lii)=-1-i*(n+1)
trsmltx=:  (((((_1-(#@])) {   (1 {:: [)) - (] mp~ (((#@]) (] dhs2lios ((_1-[),[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))) % (((#@]) (_1-(* >:)) (# @ (0 {:: [)))   ({,) (0 {:: [))) ,   ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i)), lio(lii)=-1-i*(n+1)
trsmux=:   (((((_1-(#@])) {   (1 {:: [)) - (] mp~ (((#@]) (1 dhs2lios ((_1-*),[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))) % (((#@]) (_1-(* >:)) (# @ (0 {:: [)))   ({,) (0 {:: [))) ,   ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(ui)=(1 dhs2lios (-1-i*n,i)), lio(lii)=-1-i*(n+1)
trsmu1x=:  ( (((_1-(#@])) {   (1 {:: [)) - (] mp~ (((#@]) (1 dhs2lios ((_1-*),[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))                                                          ) ,   ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(ui)=(1 dhs2lios (-1-i*n,i))
trsmu1hx=: ( ((    (#@])  {   (1 {:: [)) - (] mp~ (((#@]) (] dhs2lios (  2   #[)) (# @ (0 {:: [))) +@({,) (0 {:: [)))                                                          ) , ~ ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(ui)=(n dhs2lios (   i  ,i))
trsmu1tx=: ( ((    (#@])  {   (1 {:: [)) - (] mp~ (((#@]) (] dhs2lios (  2   #[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))                                                          ) , ~ ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(ui)=(n dhs2lios (   i  ,i))
trsmuhx=:  ((((    (#@])  {   (1 {:: [)) - (] mp~ (((#@]) (] dhs2lios (  2   #[)) (# @ (0 {:: [))) +@({,) (0 {:: [)))) % (((#@])     (* >:)  (# @ (0 {:: [))) +@({,) (0 {:: [))) , ~ ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(ui)=(n dhs2lios (   i  ,i)), lio(lii)=   i*(n+1)
trsmutx=:  ((((    (#@])  {   (1 {:: [)) - (] mp~ (((#@]) (] dhs2lios (  2   #[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))) % (((#@])     (* >:)  (# @ (0 {:: [)))   ({,) (0 {:: [))) , ~ ]) ^: (;`(#@])`(0 {.   ]))  NB. lios(ui)=(n dhs2lios (   i  ,i)), lio(lii)=   i*(n+1)
trsmxl=:   (((((_1-(c@])) {"1 (1 {:: [)) - (] mp  (((c@]) (] dhs2lios ((_1-[),[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))) % (((c@]) (_1-(* >:)) (# @ (0 {:: [)))   ({,) (0 {:: [))) ,.  ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i)), lio(lii)=-1-i*(n+1)
trsmxl1=:  ( (((_1-(c@])) {"1 (1 {:: [)) - (] mp  (((c@]) (] dhs2lios ((_1-[),[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))                                                          ) ,.  ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(n dhs2lios (-1-i  ,i))
trsmxl1t=: ( ((    (c@])  {"1 (1 {:: [)) - (] mp  (((c@]) (1 dhs2lios (    * ,[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))                                                          ) ,.~ ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (   i*n,i))
trsmxl1h=: ( ((    (c@])  {"1 (1 {:: [)) - (] mp  (((c@]) (1 dhs2lios (    * ,[)) (# @ (0 {:: [))) +@({,) (0 {:: [)))                                                          ) ,.~ ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (   i*n,i))
trsmxlt=:  ((((    (c@])  {"1 (1 {:: [)) - (] mp  (((c@]) (1 dhs2lios (    * ,[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))) % (((c@])     (* >:)  (# @ (0 {:: [)))   ({,) (0 {:: [))) ,.~ ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (   i*n,i)), lio(lii)=   i*(n+1)
trsmxlh=:  ((((    (c@])  {"1 (1 {:: [)) - (] mp  (((c@]) (1 dhs2lios (    * ,[)) (# @ (0 {:: [))) +@({,) (0 {:: [)))) % (((c@])     (* >:)  (# @ (0 {:: [))) +@({,) (0 {:: [))) ,.~ ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (   i*n,i)), lio(lii)=   i*(n+1)
trsmxu=:   ((((    (c@])  {"1 (1 {:: [)) - (] mp  (((c@]) (] dhs2lios (  2   #[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))) % (((c@])     (* >:)  (# @ (0 {:: [)))   ({,) (0 {:: [))) ,.~ ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(n dhs2lios (   i  ,i)), lio(lii)=   i*(n+1)
trsmxu1=:  ( ((    (c@])  {"1 (1 {:: [)) - (] mp  (((c@]) (] dhs2lios (  2   #[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))                                                          ) ,.~ ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(n dhs2lios (   i  ,i))
trsmxu1h=: ( ((    (c@])  {"1 (1 {:: [)) - (] mp  (((c@]) (1 dhs2lios ((_1-*),[)) (# @ (0 {:: [))) +@({,) (0 {:: [)))                                                          ) ,.  ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (-1-i*n,i))
trsmxu1t=: ( ((    (c@])  {"1 (1 {:: [)) - (] mp  (((c@]) (1 dhs2lios ((_1-*),[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))                                                          ) ,.  ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (-1-i*n,i))
trsmxuh=:  ((((    (c@])  {"1 (1 {:: [)) - (] mp  (((c@]) (1 dhs2lios ((_1-*),[)) (# @ (0 {:: [))) +@({,) (0 {:: [)))) % (((c@]) (_1-(* >:)) (# @ (0 {:: [))) +@({,) (0 {:: [))) ,.  ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (-1-i*n,i)), lio(lii)=-1-i*(n+1)
trsmxut=:  ((((    (c@])  {"1 (1 {:: [)) - (] mp  (((c@]) (1 dhs2lios ((_1-*),[)) (# @ (0 {:: [)))   ({,) (0 {:: [)))) % (((c@]) (_1-(* >:)) (# @ (0 {:: [)))   ({,) (0 {:: [))) ,.  ]) ^: (;`(c@])`(0 {."1 ]))  NB. lios(li)=(1 dhs2lios (-1-i*n,i)), lio(lii)=-1-i*(n+1)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testtrsm
NB.
NB. Description:
NB.   Test linear monomial equation solving algorithms:
NB.   - trtrs (math/lapack addon)
NB.   - trsmxxxx (math/mt addon)
NB.   by triangular matrix given
NB.
NB. Syntax:
NB.   testtrsm (A;X)
NB. where
NB.   A - n×n-matrix, triangular
NB.   X - n×n-matrix, exact solution
NB.
NB. Formula:
NB. - ferr := max(||X - exactX|| / ||X||)
NB. - berr := max(||B - op(A) * X|| / (||op(A)|| * ||X||*eps))

testtrsm=: 3 : 0
  'A X'=. y
  'L L1 U U1'=. bT=. (trl ; trl1 ; tru ; tru1) A
  'conL conL1 conU conU1'=. ((((norm1 con trtril)&.>)`((norm1 con trtril1)&.>)`((norm1 con trtriu)&.>)`((norm1 con trtriu1)&.>)) ag) bT

  ('trtrslx'   tdyad ((0&{::)`( mp       & >/)`]`(conL "_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((( mp       & >/)@[) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (L ;X)
  ('trtrsl1x'  tdyad ((0&{::)`( mp       & >/)`]`(conL1"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((( mp       & >/)@[) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (L1;X)
  ('trtrsl1hx' tdyad ((0&{::)`((mp~ ct)~ & >/)`]`(conL1"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ ct)~ & >/)@[) - ((mp~ ct) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1@|:@]))))))) (L1;X)
  ('trtrsl1tx' tdyad ((0&{::)`((mp~ |:)~ & >/)`]`(conL1"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ |:)~ & >/)@[) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1@|:@]))))))) (L1;X)
  ('trtrslhx'  tdyad ((0&{::)`((mp~ ct)~ & >/)`]`(conL "_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ ct)~ & >/)@[) - ((mp~ ct) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1@|:@]))))))) (L ;X)
  ('trtrsltx'  tdyad ((0&{::)`((mp~ |:)~ & >/)`]`(conL "_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ |:)~ & >/)@[) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1@|:@]))))))) (L ;X)
  ('trtrsux'   tdyad ((0&{::)`( mp       & >/)`]`(conU "_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((( mp       & >/)@[) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (U ;X)
  ('trtrsu1x'  tdyad ((0&{::)`( mp       & >/)`]`(conU1"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((( mp       & >/)@[) - ( mp~     (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1@|:@]))))))) (U1;X)
  ('trtrsu1hx' tdyad ((0&{::)`((mp~ ct)~ & >/)`]`(conU1"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ ct)~ & >/)@[) - ((mp~ ct) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1@|:@]))))))) (U1;X)
  ('trtrsu1tx' tdyad ((0&{::)`((mp~ |:)~ & >/)`]`(conU1"_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ |:)~ & >/)@[) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1@|:@]))))))) (U1;X)
  ('trtrsuhx'  tdyad ((0&{::)`((mp~ ct)~ & >/)`]`(conU "_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ ct)~ & >/)@[) - ((mp~ ct) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1@|:@]))))))) (U ;X)
  ('trtrsutx'  tdyad ((0&{::)`((mp~ |:)~ & >/)`]`(conU "_)`(normi@(((- (% & (normi"1@|:)) [) (1&{::))~))`(normi@((norm1t"1@|:@((((mp~ |:)~ & >/)@[) - ((mp~ |:) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1@|:@]))))))) (U ;X)
  ('trtrsxl'   tdyad ((0&{::)`( mp~      & >/)`]`(conL "_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((( mp     ~ & >/)@[) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (L ;X)
  ('trtrsxl1'  tdyad ((0&{::)`( mp~      & >/)`]`(conL1"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((( mp     ~ & >/)@[) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (L1;X)
  ('trtrsxl1h' tdyad ((0&{::)`((mp  ct)~ & >/)`]`(conL1"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  ct)~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1   @]))))))) (L1;X)
  ('trtrsxl1t' tdyad ((0&{::)`((mp  |:)~ & >/)`]`(conL1"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  |:)~ & >/)@[) - ((mp  |:) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1   @]))))))) (L1;X)
  ('trtrsxlh'  tdyad ((0&{::)`((mp  ct)~ & >/)`]`(conL "_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  ct)~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1   @]))))))) (L ;X)
  ('trtrsxlt'  tdyad ((0&{::)`((mp  |:)~ & >/)`]`(conL "_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  |:)~ & >/)@[) - ((mp  |:) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1   @]))))))) (L ;X)
  ('trtrsxu'   tdyad ((0&{::)`( mp~      & >/)`]`(conU "_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((( mp     ~ & >/)@[) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (U ;X)
  ('trtrsxu1'  tdyad ((0&{::)`( mp~      & >/)`]`(conU1"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((( mp     ~ & >/)@[) - ( mp      (0 & {::))~)) % (((FP_EPS*norm1@(0 {:: [))*(norm1t"1   @]))))))) (U1;X)
  ('trtrsxu1h' tdyad ((0&{::)`((mp  ct)~ & >/)`]`(conU1"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  ct)~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1   @]))))))) (U1;X)
  ('trtrsxu1t' tdyad ((0&{::)`((mp  |:)~ & >/)`]`(conU1"_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  |:)~ & >/)@[) - ((mp  |:) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1   @]))))))) (U1;X)
  ('trtrsxuh'  tdyad ((0&{::)`((mp  ct)~ & >/)`]`(conU "_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  ct)~ & >/)@[) - ((mp  ct) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1   @]))))))) (U ;X)
  ('trtrsxut'  tdyad ((0&{::)`((mp  |:)~ & >/)`]`(conU "_)`(normi@(((- (% & (normi"1   )) [) (1&{::))~))`(normi@((norm1t"1   @((((mp  |:)~ & >/)@[) - ((mp  |:) (0 & {::))~)) % (((FP_EPS*normi@(0 {:: [))*(norm1t"1   @]))))))) (U ;X)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testsm
NB.
NB. Description:
NB.   Adv. to make verb to test linear monomial equation
NB.   solving algorithms trsmxxxx by matrix of generator and
NB.   shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testsm
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.            mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB.   NB. with limited random matrix values' amplitudes
NB.   cocurrent 'mt'
NB.   (_1 1 0 16 _6 4 & gemat) testsm 500 500
NB.   (_1 1 0 16 _6 4 & (gemat j. gemat)) testsm 500 500

testsm=: 1 : 'EMPTY_mt_ [ (testtrsm_mt_ @ (u ; u)) ^: (=/)'
