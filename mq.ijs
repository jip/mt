NB. mq.ijs
NB. Multiply a general matrix by a matrix with orthonormal
NB. rows or columns, which is represented in factored form
NB.
NB. unmlqxx   Multiply a general matrix by a matrix with
NB.           orthonormal rows, which is represented in
NB.           factored form, as returned by gelqf
NB. unmqlxx   Multiply a general matrix by a matrix with
NB.           orthonormal columns, which is represented in
NB.           factored form, as returned by geqlf
NB. unmqrxx   Multiply a general matrix by a matrix with
NB.           orthonormal columns, which is represented in
NB.           factored form, as returned by geqrf
NB. unmrqxx   Multiply a general matrix by a matrix with
NB.           orthonormal rows, which is represented in
NB.           factored form, as returned by gerqf
NB. unmhrxxx  Multiply a general matrix by an unitary
NB.           (orthogonal) matrix, which is represented in
NB.           factored form, as returned by gehrdl and gehrdu
NB. unmbrxx   Multiply a general matrix by an unitary
NB.           (orthogonal) matrix, which is represented in
NB.           factored form, as returned by gebrd
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-01-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

UNMQBS=: 4 NB. 32   NB. block size limit

NB. ---------------------------------------------------------
NB. Restrained version of Take verb ({.)
NB. Quick-and simple bonded version. FIXME!

rt=: 2 : 'm&{. ^: ((|m)<v)'

NB. ---------------------------------------------------------
NB. Description:
NB.   Do a partial step of unm{l2,2l,2r,r2}{l,r}{n,c}
NB.
NB. Syntax
NB.   bCi1=. Qf unmxxxx (pfxCi ; sfxCi)
NB. where
NB.   pfxCi - pfxC(i), the prefix of eC(i), either already
NB.           processed or not yet processed part
NB.   sfxCi - sfxC(i), the suffix of eC(i), either not yet
NB.           processed or already processed part, inversely
NB.           to pfxC(i)
NB.   eC(i) - matrix C augmented by trash vector, after i-th
NB.           and before (i+1)-th step, it may be restored by
NB.           merging pfxC(i) and sfxC(i)
NB.   C     - m×n-matrix to multiply
NB.   Qf    - unit triangular matrix, it represents Q in
NB.           factored form, and contains vectors vtau[0:k-1]
NB.   Q     - matrix with orthonormal rows or columns, which
NB.           is defined as the product of elementary
NB.           reflectors
NB.   bCi1  - bC(i+1), either pfxC(i+1) or sfxC(i+1), being
NB.           pfxC(i) or sfxC(i) updated by i-th elementary
NB.           reflector from Qf
NB.
NB. Algorithm:
NB.   In:  Qf , (pfxCi;sfxCi) as (aCi;bCi) or (bCi;aCi)
NB.   Out: bCi1
NB.   0) form rios, rIOS of vtau[i], depending on aC(i)'s
NB.      size
NB.   1) extract vtau[i] from Qf and ravel it:
NB.        vtaui=. rios (, ;. 0) Qf
NB.   2) apply an elementary reflector defined by vtau[i] to
NB.      bC(i) to produce bC(i+1):
NB.        bCi1=. vtaui larfxxxx bCi
NB.
NB. Notes:
NB. - this is only a 1st part of step indeed; finalizing
NB.   action for last and non-last iterations vary, so it is
NB.   factored out to caller

unml2lnstep=: (,;.0~ (1 _ ,:~ 2 #      #@(0&({::)))) larflcfr 1 {:: ]
unml2lcstep=: (,;.0~ (1 _ ,:~ 2 #      #@(0&({::)))) larflnfr 1 {:: ]
unml2rnstep=: (,;.0~ (1 _ ,:~ 2 #      c@(0&({::)))) larfrcfr 1 {:: ]
unml2rcstep=: (,;.0~ (1 _ ,:~ 2 #      c@(0&({::)))) larfrnfr 1 {:: ]

unm2llnstep=: (,;.0~ (_ 1 ,:~ 2 # _1 - #@(1&({::)))) larflnbc 0 {:: ]
unm2llcstep=: (,;.0~ (_ 1 ,:~ 2 # _1 - #@(1&({::)))) larflcbc 0 {:: ]
unm2lrnstep=: (,;.0~ (_ 1 ,:~ 2 # _1 - c@(1&({::)))) larfrnbc 0 {:: ]
unm2lrcstep=: (,;.0~ (_ 1 ,:~ 2 # _1 - c@(1&({::)))) larfrcbc 0 {:: ]

unm2rlnstep=: (,;.0~ (_ 1 ,:~ 2 #      #@(0&({::)))) larflnfc 1 {:: ]
unm2rlcstep=: (,;.0~ (_ 1 ,:~ 2 #      #@(0&({::)))) larflcfc 1 {:: ]
unm2rrnstep=: (,;.0~ (_ 1 ,:~ 2 #      c@(0&({::)))) larfrnfc 1 {:: ]
unm2rrcstep=: (,;.0~ (_ 1 ,:~ 2 #      c@(0&({::)))) larfrcfc 1 {:: ]

unmr2lnstep=: (,;.0~ (1 _ ,:~ 2 # _1 - #@(1&({::)))) larflcbr 0 {:: ]
unmr2lcstep=: (,;.0~ (1 _ ,:~ 2 # _1 - #@(1&({::)))) larflnbr 0 {:: ]
unmr2rnstep=: (,;.0~ (1 _ ,:~ 2 # _1 - c@(1&({::)))) larfrcbr 0 {:: ]
unmr2rcstep=: (,;.0~ (1 _ ,:~ 2 # _1 - c@(1&({::)))) larfrnbr 0 {:: ]

NB. ---------------------------------------------------------
NB. Description:
NB.   Do a step of unm{lq,ql,qr,rq}{l,r}{n,c}
NB.
NB. Syntax
NB.   bCi1=. Qf unmxxxx (pfxCi ; sfxCi)
NB. where
NB.   pfxCi - pfxC(i), the prefix of eC(i), either already
NB.           processed or not yet processed part
NB.   sfxCi - sfxC(i), the suffix of eC(i), either not yet
NB.           processed or already processed part, inversely
NB.           to pfxC(i)
NB.   eC(i) - matrix C augmented by trash vector, after i-th
NB.           and before (i+1)-th step, it may be restored by
NB.           merging pfxC(i) and sfxC(i)
NB.   C     - m×n-matrix to multiply
NB.   Qf    - unit triangular matrix, it represents Q in
NB.           factored form, and contains vectors vtau[0:k-1]
NB.   Q     - matrix with orthonormal rows or columns, which
NB.           is defined as the product of elementary
NB.           reflectors
NB.   bCi1  - bC(i+1), either pfxC(i+1) or sfxC(i+1), being
NB.           pfxC(i) or sfxC(i) updated by i-th elementary
NB.           reflector from Qf
NB.
NB. Algorithm:
NB.   In:  Qf , (pfxCi;sfxCi) as (aCi;bCi) or (bCi;aCi)
NB.   Out: bCi1
NB.   0) form rios, rIOS of matrix VTau[i] which defines the
NB.      block reflector, it depends on aC(i)'s size
NB.   1) try to extract VTau[i] from Qf:
NB.        VTaui=. rios (] ;. 0) Qf
NB.      1.0) if failure occures (i.e. execution is at last
NB.           step and 0<UNMQBS|k ), then replace in rios
NB.           wrong explicit length (UNMQBS) by implicit one
NB.           (∞) and retry
NB.   2) apply a block reflector defined by VTaui to bCi:
NB.        bCi1=. VTaui larfbxxxx bCi

unmlqlnstep=: (];.0~ ::(];.0~ _&(2:})) ((UNMQBS, _) ,:~ 2 #             #@(0&({::)))) larfblcfr 1 {:: ]
unmlqlcstep=: (];.0~ ::(];.0~ _&(2:})) ((UNMQBS, _) ,:~ 2 #             #@(0&({::)))) larfblnfr 1 {:: ]
unmlqrnstep=: (];.0~ ::(];.0~ _&(2:})) ((UNMQBS, _) ,:~ 2 #             c@(0&({::)))) larfbrcfr 1 {:: ]
unmlqrcstep=: (];.0~ ::(];.0~ _&(2:})) ((UNMQBS, _) ,:~ 2 #             c@(0&({::)))) larfbrnfr 1 {:: ]

unmqllnstep=: (];.0~ ::(];.0~ _&(3:})) ((UNMQBS,~_) ,:~ 2 # (-UNMQBS) - #@(1&({::)))) larfblnbc 0 {:: ]
unmqllcstep=: (];.0~ ::(];.0~ _&(3:})) ((UNMQBS,~_) ,:~ 2 # (-UNMQBS) - #@(1&({::)))) larfblcbc 0 {:: ]
unmqlrnstep=: (];.0~ ::(];.0~ _&(3:})) ((UNMQBS,~_) ,:~ 2 # (-UNMQBS) - c@(1&({::)))) larfbrnbc 0 {:: ]
unmqlrcstep=: (];.0~ ::(];.0~ _&(3:})) ((UNMQBS,~_) ,:~ 2 # (-UNMQBS) - c@(1&({::)))) larfbrcbc 0 {:: ]

unmqrlnstep=: (];.0~ ::(];.0~ _&(3:})) ((UNMQBS,~_) ,:~ 2 #             #@(0&({::)))) larfblnfc 1 {:: ]
unmqrlcstep=: (];.0~ ::(];.0~ _&(3:})) ((UNMQBS,~_) ,:~ 2 #             #@(0&({::)))) larfblcfc 1 {:: ]
unmqrrnstep=: (];.0~ ::(];.0~ _&(3:})) ((UNMQBS,~_) ,:~ 2 #             c@(0&({::)))) larfbrnfc 1 {:: ]
unmqrrcstep=: (];.0~ ::(];.0~ _&(3:})) ((UNMQBS,~_) ,:~ 2 #             c@(0&({::)))) larfbrcfc 1 {:: ]

unmrqlnstep=: (];.0~ ::(];.0~ _&(2:})) ((UNMQBS, _) ,:~ 2 # (-UNMQBS) - #@(1&({::)))) larfblcbr 0 {:: ]
unmrqlcstep=: (];.0~ ::(];.0~ _&(2:})) ((UNMQBS, _) ,:~ 2 # (-UNMQBS) - #@(1&({::)))) larfblnbr 0 {:: ]
unmrqrnstep=: (];.0~ ::(];.0~ _&(2:})) ((UNMQBS, _) ,:~ 2 # (-UNMQBS) - c@(1&({::)))) larfbrcbr 0 {:: ]
unmrqrcstep=: (];.0~ ::(];.0~ _&(2:})) ((UNMQBS, _) ,:~ 2 # (-UNMQBS) - c@(1&({::)))) larfbrnbr 0 {:: ]

NB. ---------------------------------------------------------
NB. Verb     Action  Side   Tran  Syntax
NB. unml2ln  Q * C   left   none  eCprod=. Qf unml2ln (C, 0)
NB. unml2lc  Q'* C   left   ct    eCprod=. Qf unml2lc (C, 0)
NB. unml2rn  C * Q   right  none  eCprod=. Qf unml2rn (C,.0)
NB. unml2rc  C * Q'  right  ct    eCprod=. Qf unml2rc (C,.0)
NB. unm2lln  Q * C   left   none  eCprod=. Qf unm2lln (0, C)
NB. unm2llc  Q'* C   left   ct    eCprod=. Qf unm2llc (0, C)
NB. unm2lrn  C * Q   right  none  eCprod=. Qf unm2lrn (0,.C)
NB. unm2lrc  C * Q'  right  ct    eCprod=. Qf unm2lrc (0,.C)
NB. unm2rln  Q * C   left   none  eCprod=. Qf unm2rln (C, 0)
NB. unm2rlc  Q'* C   left   ct    eCprod=. Qf unm2rlc (C, 0)
NB. unm2rrn  C * Q   right  none  eCprod=. Qf unm2rrn (C,.0)
NB. unm2rrc  C * Q'  right  ct    eCprod=. Qf unm2rrc (C,.0)
NB. unmr2ln  Q * C   left   none  eCprod=. Qf unmr2ln (0, C)
NB. unmr2lc  Q'* C   left   ct    eCprod=. Qf unmr2lc (0, C)
NB. unmr2rn  C * Q   right  none  eCprod=. Qf unmr2rn (0,.C)
NB. unmr2rc  C * Q'  right  ct    eCprod=. Qf unmr2rc (0,.C)
NB.
NB. Description:
NB.   Multiply a general matrix C, augmented by trash vector,
NB.   by matrix Q. This is non-blocked version of algorithm
NB. where
NB.   C      - m×n-matrix to multiply
NB.   Qf     - unit triangular matrix, it represents Q in
NB.            factored form as returned by ge{lq,ql,qr,rq}2,
NB.            and contains vectors vtau[0:k-1]
NB.   Q      - matrix with orthonormal rows or columns, which
NB.            is defined as the product of elementary
NB.            reflectors
NB.   eCprod - being product of matrix Q and augmented matrix
NB.            C, trash vector is modified on exit
NB.
NB. Algorithm:
NB.   In: Qf , eC
NB.   Out: eCprod
NB.   0) form (pfxC;sfxC) as (aC;bC) or (bC;aC)
NB.   1) find I, the decremented number of iterations
NB.   2) start iterations via power (^:) on (pfxCi;sfxCi) as
NB.      (aCi;bCi) or (bCi;aCi)
NB.      2.0) apply unmxxxxstep:
NB.             tmp0=. Qf unmxxxxstep (pfxCi;sfxCi)
NB.      2.1) combine aCi and tmp0 to produce (pfxCi1;sfxCi1)
NB.           as (aCi1;bCi1) or (bCi1;aCi1)
NB.   3) apply unmxxxxstep to (pfxCi1;sfxCi1) to produce bCi1
NB.   4) combine aCi1 and bCi1 to produce eCprod
NB.
NB. If:
NB.   2 -: # $ A        NB. A is a 2-rank array (i.e. matrix)
NB.   -:/ @ $ A         NB. A is a square matrix (it's not
NB.                     NB.   necessary and is assumed just
NB.                     NB.   to simplify assertions)
NB.   n=. # A           NB. size of matrix A
NB. then
NB.   (idmat n) (-: (clean@(_1  0&}.))) ((   tru1 gelqf A) unml2ln ((ct unglq gelqf A) ,   0))
NB.   (idmat n) (-: (clean@(_1  0&}.))) ((   tru1 gelqf A) unml2lc ((   unglq gelqf A) ,   0))
NB.   (idmat n) (-: (clean@( 0 _1&}.))) ((   tru1 gelqf A) unml2rn ((ct unglq gelqf A) ,.  0))
NB.   (idmat n) (-: (clean@( 0 _1&}.))) ((   tru1 gelqf A) unml2rc ((   unglq gelqf A) ,.  0))
NB.   (idmat n) (-: (clean@( 1  0&}.))) ((_1 tru1 geqlf A) unm2lln ((ct ungql geqlf A) , ~ 0))
NB.   (idmat n) (-: (clean@( 1  0&}.))) ((_1 tru1 geqlf A) unm2llc ((   ungql geqlf A) , ~ 0))
NB.   (idmat n) (-: (clean@( 0  1&}.))) ((_1 tru1 geqlf A) unm2lrn ((ct ungql geqlf A) ,.~ 0))
NB.   (idmat n) (-: (clean@( 0  1&}.))) ((_1 tru1 geqlf A) unm2lrc ((   ungql geqlf A) ,.~ 0))
NB.   (idmat n) (-: (clean@(_1  0&}.))) ((   trl1 geqrf A) unm2rln ((ct ungqr geqrf A) ,   0))
NB.   (idmat n) (-: (clean@(_1  0&}.))) ((   trl1 geqrf A) unm2rlc ((   ungqr geqrf A) ,   0))
NB.   (idmat n) (-: (clean@( 0 _1&}.))) ((   trl1 geqrf A) unm2rrn ((ct ungqr geqrf A) ,.  0))
NB.   (idmat n) (-: (clean@( 0 _1&}.))) ((   trl1 geqrf A) unm2rrc ((   ungqr geqrf A) ,.  0))
NB.   (idmat n) (-: (clean@( 1  0&}.))) (( 1 trl1 gerqf A) unmr2ln ((ct ungrq gerqf A) , ~ 0))
NB.   (idmat n) (-: (clean@( 1  0&}.))) (( 1 trl1 gerqf A) unmr2lc ((   ungrq gerqf A) , ~ 0))
NB.   (idmat n) (-: (clean@( 0  1&}.))) (( 1 trl1 gerqf A) unmr2rn ((ct ungrq gerqf A) ,.~ 0))
NB.   (idmat n) (-: (clean@( 0  1&}.))) (( 1 trl1 gerqf A) unmr2rc ((   ungrq gerqf A) ,.~ 0))
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - emulates LAPACK's xUNM{L2,2L,2R,R2}
NB. - unml2{ln,lc,rn,rc} and unmlq{ln,lc,rn,rc} respectively
NB.   are topologic equivalents

unml2ln=: [((0{::]), unml2lnstep)( (((0{::])(( ,  {.  ); ( 1  0&}.@]))unml2lnstep)^:(<:@#@[))(;~(0&{.)))
unml2lc=: [((0{::]), unml2lcstep)([(((0{::])(((, ~{:)~);~(_1  0&}.@[))unml2lcstep)^:(<:@#@[))(({.; }.)~(_1+#))~)
unml2rn=: [((0{::]),.unml2rnstep)([(((0{::])(((,.~tc)~);~( 0 _1&}.@[))unml2rnstep)^:(<:@#@[))((({.~(_&,)); (}.~(0&,)))(_1+#))~)
unml2rc=: [((0{::]),.unml2rcstep)( (((0{::])(( ,. hc  ); ( 0  1&}.@]))unml2rcstep)^:(<:@#@[))(;~(_ 0&{.)))

unm2lln=: [(unm2llnstep, (1{::]))([((unm2llnstep(( ,  {.  ); ( 1  0&}.@]))(1{::]))^:(<:@c@[))((({.;~}.)~( 1-c))~))
unm2llc=: [(unm2llcstep, (1{::]))( ((unm2llcstep(((, ~{:)~);~(_1  0&}.@[))(1{::]))^:(<:@c@[))(; (0&{.)))
unm2lrn=: [(unm2lrnstep,.(1{::]))( ((unm2lrnstep(((,.~tc)~);~( 0 _1&}.@[))(1{::]))^:(<:@c@[))(; (_ 0&{.)))
unm2lrc=: [(unm2lrcstep,.(1{::]))([((unm2lrcstep(( ,. hc  ); ( 0  1&}.@]))(1{::]))^:(<:@c@[))((({.~(_&,));~(}.~(0&,)))( 1-c))~)

unm2rln=: [((0{::]), unm2rlnstep)([(((0{::])(((, ~{:)~);~(_1  0&}.@[))unm2rlnstep)^:(<:@c@[))((({.; }.)~(_1+c))~))
unm2rlc=: [((0{::]), unm2rlcstep)( (((0{::])(( ,  {.  ); ( 1  0&}.@]))unm2rlcstep)^:(<:@c@[))(;~(0&{.)))
unm2rrn=: [((0{::]),.unm2rrnstep)( (((0{::])(( ,. hc  ); ( 0  1&}.@]))unm2rrnstep)^:(<:@c@[))(;~(_ 0&{.)))
unm2rrc=: [((0{::]),.unm2rrcstep)([(((0{::])(((,.~tc)~);~( 0 _1&}.@[))unm2rrcstep)^:(<:@c@[))((({.~(_&,)); (}.~(0&,)))(_1+c))~)

unmr2ln=: [(unmr2lnstep, (1{::]))( ((unmr2lnstep(((, ~{:)~);~(_1  0&}.@[))(1{::]))^:(<:@#@[))(; (0&{.)))
unmr2lc=: [(unmr2lcstep, (1{::]))([((unmr2lcstep(( ,  {.  ); ( 1  0&}.@]))(1{::]))^:(<:@#@[))((({.;~}.)~( 1-#))~))
unmr2rn=: [(unmr2rnstep,.(1{::]))([((unmr2rnstep(( ,. hc  ); ( 0  1&}.@]))(1{::]))^:(<:@#@[))((({.~(_&,));~(}.~(0&,)))( 1-#))~)
unmr2rc=: [(unmr2rcstep,.(1{::]))( ((unmr2rcstep(((,.~tc)~);~( 0 _1&}.@[))(1{::]))^:(<:@#@[))(; (_ 0&{.)))

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb     Action  Side   Tran  Syntax
NB. unmlqln  Q * C   left   none  QeC=.  Qf unmlqln C
NB. unmlqlc  Q'* C   left   ct    cQeC=. Qf unmlqlc C
NB. unmlqrn  C * Q   right  none  eCQ=.  Qf unmlqlc C
NB. unmlqrc  C * Q'  right  ct    eCcQ=. Qf unmlqlc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form Qf as returned by gelqf
NB. where
NB.
NB. Algorithm:
NB.   iters==⌊(#+BS-1)/BS⌋
NB.
NB. If:
NB.   LQf=. gelqf A
NB.   L=. trl (0 _1 }. LQf)
NB.   Qf=. tru1 LQf
NB.   Q=. unglq LQf
NB.   k=. <./ $ A
NB. then
NB.   ((idmat @ c) (-: clean) ((tru1         ) unmlqln (ct @ unglq)) @ gelqf) A
NB.   ((idmat @ c) (-: clean) ((tru1 @ (k&{.)) unmlqlc (     unglq)) @ gelqf) A
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - emulates LAPACK's xUNMLQ
NB. - unml2{lc,ln,rc,rn} and unmlq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmlqln=: _1 0 }. ((unml2ln`(0 {:: (  ((0 {:: ]) ((,  (  UNMQBS  rt #)) ; (UNMQBS* 1  0)&}.@]) unmlqlnstep)^:(<.@(UNMQBS %~ (<:UNMQBS) + #@[)) (    ;~     0&{.                      )))@.(UNMQBS < #@[)~ (tru1@({.~(<./@(0 _1 + $)))))~ ,&0)
unmlqlc=: _1 0 }. ((unml2lc`(1 {:: ([ ((0 {:: ]) ((, ~((-UNMQBS) rt #))~;~(UNMQBS*_1  0)&}.@[) unmlqlcstep)^:(<.@(UNMQBS %~ (<:UNMQBS) + #@[)) ((({.; }.)~(UNMQBS rounddown (_1+#)))~)))@.(UNMQBS < #@[)~ (tru1@({.~(<./@(0 _1 + $)))))~ ,&0)




NB.    timespacex & > ('HCln=. Qf unml2ln ecA';'HClc=. Qf unml2lc ecA';'HCrn=. Qf unml2rn ecA2';'HCrno=. Qf unml2rno ecA2';'HCrc=. Qf unml2rc ecA2';'HCrco=. Qf unml2rco ecA2')
NB. 1.84525 7.35002e6
NB. 1.43646 9.46349e6
NB. 1.49872 9.46349e6
NB. 1.99539 8.39872e6
NB. 1.43921 7.35002e6
NB. 2.01253 8.39878e6
NB.    load '/home/jip/j602-user/projects/mt/mt.ijs'
NB.    cocurrent 'mt'
NB.    A=. j./ _9 + ? 2 10 10 $ 19
NB.    B=. j./ _9 + ? 2 7 10 $ 19
NB.    C=. j./ _9 + ? 2 10 7 $ 19
