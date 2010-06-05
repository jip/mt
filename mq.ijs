NB. Multiply a general matrix by a matrix with orthonormal
NB. rows or columns, which is represented in factored form
NB.
NB. unmlqxx    Multiply a general matrix by a matrix with
NB.            orthonormal rows, which is represented in
NB.            factored form, as returned by gelqf
NB. unmqlxx    Multiply a general matrix by a matrix with
NB.            orthonormal columns, which is represented in
NB.            factored form, as returned by geqlf
NB. unmqrxx    Multiply a general matrix by a matrix with
NB.            orthonormal columns, which is represented in
NB.            factored form, as returned by geqrf
NB. unmrqxx    Multiply a general matrix by a matrix with
NB.            orthonormal rows, which is represented in
NB.            factored form, as returned by gerqf
NB. unmhrxxx   Multiply a general matrix by an unitary
NB.            (orthogonal) matrix, which is represented in
NB.            factored form, as returned by gehrdx
NB. unmbrxx    Multiply a general matrix by an unitary
NB.            (orthogonal) matrix, which is represented in
NB.            factored form, as returned by gebrdx
NB.
NB. testunmq   Test Q multiplication qf-algorithms by general
NB.            matrix given
NB. testunmhr  Test Q multiplication hrd-algorithms by square
NB.            matrix given
NB. testmq     Adv. to make verb to test unmxxxxx by matrix of
NB.            generator and shape given
NB.
NB. Copyright (C) 2010 Igor Zhuravlov
NB. For license terms, see the file COPYING in this distribution
NB. Version: 1.0.0 2010-06-01

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

MQBS=: 32   NB. block size limit

NB. ---------------------------------------------------------

arounddown=: 1 : '- m&|'  NB. adverb to round down by an integer constant

NB. ---------------------------------------------------------
NB. Description:
NB.   Single step of non-blocked version of algorithms
NB.
NB. Syntax
NB.   'pfxCi1 sfxCi1'=. Qf unmxxxx (pfxCi ; sfxCi)
NB. where
NB.   pfxCi   - pfxC(i), the prefix of eC(i), either already
NB.             processed or not yet processed part
NB.   sfxCi   - sfxC(i), the suffix of eC(i), either not yet
NB.             processed or already processed part,
NB.             inversely to pfxC(i)
NB.   eC(i)   - matrix C augmented by trash vector, after
NB.             i-th and before (i+1)-th step, it may be
NB.             restored by merging pfxC(i) and sfxC(i)
NB.   C       - m×n-matrix to multiply
NB.   Qf      - unit triangular matrix, it represents Q in
NB.             factored form, and contains vectors
NB.             vtau[0:k-1]
NB.   Q       - matrix with orthonormal rows or columns,
NB.             which is defined as the product of
NB.             elementary reflectors
NB.   pfxCi1  - pfxC(i+1), the prefix of eC(i+1), either
NB.             already processed or not yet processed part
NB.   sfxCi1  - sfxC(i+1), the suffix of eC(i+1), either not
NB.             yet processed or already processed part,
NB.             inversely to pfxC(i)
NB.   eC(i+1) - matrix C augmented by modified trash vector,
NB.             after (i+1)-th step, it may be restored by
NB.             merging pfxC(i+1) and sfxC(i+1)
NB.
NB. Algorithm:
NB.   In:  Qf pfxC(i) sfxC(i)
NB.   Out: pfxC(i+1) sfxC(i+1)
NB.   0) form rios, rIOS of vtau[i] which defines an
NB.      elementary reflector, depending on aC(i)'s size
NB.   1) extract vtau[i] from Qf and ravel it:
NB.        vtaui=. rios (, ;. 0) Qf
NB.   2) apply an elementary reflector defined by vtau[i] to
NB.      either pfxC(i) or sfxC(i) to produce tmp:
NB.        tmp=. vtaui larfxxxx xfxCi
NB.   3) combine tmp and either pfxC(i) or sfxC(i) to
NB.      produce pfxC(i+1) and sfxC(i+1)

unml2lnstep=: (0 {:: ]) ((,    1  _&rt)  ;  1  0&}.@])  (,;.0~ (1 _ ,:~ 2 #      #@(0&({::)))) larflcfr 1 {:: ]
unml2lcstep=: (0 {:: ]) ((, ~ _1  _&rt)~ ;~_1  0&}.@[)  (,;.0~ (1 _ ,:~ 2 #      #@(0&({::)))) larflnfr 1 {:: ]
unml2rnstep=: (0 {:: ]) ((,.~  _ _1&rt)~ ;~ 0 _1&}.@[)  (,;.0~ (1 _ ,:~ 2 #      c@(0&({::)))) larfrcfr 1 {:: ]
unml2rcstep=: (0 {:: ]) ((,.   _  1&rt)  ;  0  1&}.@])  (,;.0~ (1 _ ,:~ 2 #      c@(0&({::)))) larfrnfr 1 {:: ]

unm2llnstep=: (1 {:: ]) ((,    1  _&rt)  ;  1  0&}.@])~ (,;.0~ (_ 1 ,:~ 2 # _1 - #@(1&({::)))) larflnbc 0 {:: ]
unm2llcstep=: (1 {:: ]) ((, ~ _1  _&rt)~ ;~_1  0&}.@[)~ (,;.0~ (_ 1 ,:~ 2 # _1 - #@(1&({::)))) larflcbc 0 {:: ]
unm2lrnstep=: (1 {:: ]) ((,.~  _ _1&rt)~ ;~ 0 _1&}.@[)~ (,;.0~ (_ 1 ,:~ 2 # _1 - c@(1&({::)))) larfrnbc 0 {:: ]
unm2lrcstep=: (1 {:: ]) ((,.   _  1&rt)  ;  0  1&}.@])~ (,;.0~ (_ 1 ,:~ 2 # _1 - c@(1&({::)))) larfrcbc 0 {:: ]

unm2rlnstep=: (0 {:: ]) ((, ~ _1  _&rt)~ ;~_1  0&}.@[)  (,;.0~ (_ 1 ,:~ 2 #      #@(0&({::)))) larflnfc 1 {:: ]
unm2rlcstep=: (0 {:: ]) ((,    1  _&rt)  ;  1  0&}.@])  (,;.0~ (_ 1 ,:~ 2 #      #@(0&({::)))) larflcfc 1 {:: ]
unm2rrnstep=: (0 {:: ]) ((,.   _  1&rt)  ;  0  1&}.@])  (,;.0~ (_ 1 ,:~ 2 #      c@(0&({::)))) larfrnfc 1 {:: ]
unm2rrcstep=: (0 {:: ]) ((,.~  _ _1&rt)~ ;~ 0 _1&}.@[)  (,;.0~ (_ 1 ,:~ 2 #      c@(0&({::)))) larfrcfc 1 {:: ]

unmr2lnstep=: (1 {:: ]) ((, ~ _1  _&rt)~ ;~_1  0&}.@[)~ (,;.0~ (1 _ ,:~ 2 # _1 - #@(1&({::)))) larflcbr 0 {:: ]
unmr2lcstep=: (1 {:: ]) ((,    1  _&rt)  ;  1  0&}.@])~ (,;.0~ (1 _ ,:~ 2 # _1 - #@(1&({::)))) larflnbr 0 {:: ]
unmr2rnstep=: (1 {:: ]) ((,.   _  1&rt)  ;  0  1&}.@])~ (,;.0~ (1 _ ,:~ 2 # _1 - c@(1&({::)))) larfrcbr 0 {:: ]
unmr2rcstep=: (1 {:: ]) ((,.~  _ _1&rt)~ ;~ 0 _1&}.@[)~ (,;.0~ (1 _ ,:~ 2 # _1 - c@(1&({::)))) larfrnbr 0 {:: ]

NB. ---------------------------------------------------------
NB. Verb     Action   Side   Tran  Syntax
NB. unml2ln  Q * C    left   none  eCprod=. Qf unml2ln (C, 0)
NB. unml2lc  Q^H * C  left   ct    eCprod=. Qf unml2lc (C, 0)
NB. unml2rn  C * Q    right  none  eCprod=. Qf unml2rn (C,.0)
NB. unml2rc  C * Q^H  right  ct    eCprod=. Qf unml2rc (C,.0)
NB. unm2lln  Q * C    left   none  eCprod=. Qf unm2lln (0, C)
NB. unm2llc  Q^H * C  left   ct    eCprod=. Qf unm2llc (0, C)
NB. unm2lrn  C * Q    right  none  eCprod=. Qf unm2lrn (0,.C)
NB. unm2lrc  C * Q^H  right  ct    eCprod=. Qf unm2lrc (0,.C)
NB. unm2rln  Q * C    left   none  eCprod=. Qf unm2rln (C, 0)
NB. unm2rlc  Q^H * C  left   ct    eCprod=. Qf unm2rlc (C, 0)
NB. unm2rrn  C * Q    right  none  eCprod=. Qf unm2rrn (C,.0)
NB. unm2rrc  C * Q^H  right  ct    eCprod=. Qf unm2rrc (C,.0)
NB. unmr2ln  Q * C    left   none  eCprod=. Qf unmr2ln (0, C)
NB. unmr2lc  Q^H * C  left   ct    eCprod=. Qf unmr2lc (0, C)
NB. unmr2rn  C * Q    right  none  eCprod=. Qf unmr2rn (0,.C)
NB. unmr2rc  C * Q^H  right  ct    eCprod=. Qf unmr2rc (0,.C)
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
NB.   2) start iterations by Power (^:) on (pfxCi;sfxCi) as
NB.      (aCi;bCi) or (bCi;aCi)
NB.      2.0) apply unmxxxxstep:
NB.             tmp=. Qf unmxxxxstep (pfxCi;sfxCi)
NB.      2.1) combine aCi and tmp to produce (pfxCi1;sfxCi1)
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
NB. - implements LAPACK's xUNM{L2,2L,2R,R2}
NB. - unml2{ln,lc,rn,rc} and unmlq{ln,lc,rn,rc} respectively
NB.   are topologic equivalents

unml2ln=: , &>/@(  unml2lnstep^:(#@[) (;~ 0  &{.                         ) )
unml2lc=: , &>/@([ unml2lcstep^:(#@[) (( {.       ;   }.      )~ (_1 + #))~)
unml2rn=: ,.&>/@([ unml2rnstep^:(#@[) ((({.~ _&,) ;  (}.~ 0&,))  (_1 + #))~)
unml2rc=: ,.&>/@(  unml2rcstep^:(#@[) (;~ _ 0&{.                         ) )

unm2lln=: , &>/@([ unm2llnstep^:(c@[) (( {.       ;~  }.      )~ ( 1 - c))~)
unm2llc=: , &>/@(  unm2llcstep^:(c@[) (;  0  &{.                         ) )
unm2lrn=: ,.&>/@(  unm2lrnstep^:(c@[) (;  _ 0&{.                         ) )
unm2lrc=: ,.&>/@([ unm2lrcstep^:(c@[) ((({.~ _&,) ;~ (}.~ 0&,))  ( 1 - c))~)

unm2rln=: , &>/@([ unm2rlnstep^:(c@[) (( {.       ;   }.      )~ (_1 + c))~)
unm2rlc=: , &>/@(  unm2rlcstep^:(c@[) (;~ 0  &{.                         ) )
unm2rrn=: ,.&>/@(  unm2rrnstep^:(c@[) (;~ _ 0&{.                         ) )
unm2rrc=: ,.&>/@([ unm2rrcstep^:(c@[) ((({.~ _&,) ;  (}.~ 0&,))  (_1 + c))~)

unmr2ln=: , &>/@(  unmr2lnstep^:(#@[) (;  0&{.                           ) )
unmr2lc=: , &>/@([ unmr2lcstep^:(#@[) (( {.       ;~  }.      )~ ( 1 - #))~)
unmr2rn=: ,.&>/@([ unmr2rnstep^:(#@[) ((({.~ _&,) ;~ (}.~ 0&,))  ( 1 - #))~)
unmr2rc=: ,.&>/@(  unmr2rcstep^:(#@[) (;  _ 0&{.                         ) )

NB. ---------------------------------------------------------
NB. Description:
NB.   Single step of algorithms
NB.
NB. Syntax
NB.   'pfxCi1 sfxCi1'=. Qf unmxxxx (pfxCi ; sfxCi)
NB. where
NB.   pfxCi   - pfxC(i), the prefix of eC(i), either already
NB.             processed or not yet processed part
NB.   sfxCi   - sfxC(i), the suffix of eC(i), either not yet
NB.             processed or already processed part,
NB.             inversely to pfxC(i)
NB.   eC(i)   - matrix C augmented by trash vector, after
NB.             i-th and before (i+1)-th step, it may be
NB.             restored by merging pfxC(i) and sfxC(i)
NB.   C       - m×n-matrix to multiply
NB.   Qf      - unit triangular matrix, it represents Q in
NB.             factored form, and contains vectors
NB.             vtau[0:k-1]
NB.   Q       - matrix with orthonormal rows or columns,
NB.             which is defined as the product of
NB.             elementary reflectors
NB.   pfxCi1  - pfxC(i+1), the prefix of eC(i+1), either
NB.             already processed or not yet processed part
NB.   sfxCi1  - sfxC(i+1), the suffix of eC(i+1), either not
NB.             yet processed or already processed part,
NB.             inversely to pfxC(i)
NB.   eC(i+1) - matrix C augmented by modified trash vector,
NB.             after (i+1)-th step, it may be restored by
NB.             merging pfxC(i+1) and sfxC(i+1)
NB.
NB. Algorithm:
NB.   In:  Qf pfxC(i) sfxC(i)
NB.   Out: pfxC(i+1) sfxC(i+1)
NB.   0) form rios, rIOS of matrix VTau[i] which defines the
NB.      block reflector, it depends on aC(i)'s size
NB.   1) try to extract VTau[i] from Qf:
NB.        VTaui=. rios (] ;. 0) Qf
NB.      1.0) if failure occures (i.e. execution is at first
NB.           or last step and 0<MQBS|k ), then replace in
NB.           rios wrong explicit length (MQBS) by implicit
NB.           one (∞) and retry
NB.   2) apply a block reflector defined by VTau[i] to either
NB.      pfxC(i) or sfxC(i) to produce tmp:
NB.        tmp=. VTaui larfbxxxx xfxCi
NB.   3) combine tmp and either pfxC(i) or sfxC(i) to
NB.      produce pfxC(i+1) and sfxC(i+1)
NB.
NB. Notes:
NB. - assigning MQBS=:1 transforms this verbs to their
NB.   non-blocked counterparts

unmlqlnstep=:                          (0 {:: ])              ((,   (MQBS* 1  _)&rt)             ;  (MQBS* 1  0)&}.@])  (  ];.0~ ::(];.0~ _&(2:})) ((MQBS*1 _) ,:~ 2 # #@(0&({::))                                                      )) larfblcfr 1 {:: ]
unmlqlcstep=:                          (0 {:: ])              ((, ~ (MQBS*_1  _)&rt)~            ;~ (MQBS*_1  0)&}.@[)  (  ];.0~ ::(];.0~ _&(2:})) ((MQBS*1 _) ,:~ 2 # #@(0&({::))                                                      )) larfblnfr 1 {:: ]
unmlqrnstep=:                          (0 {:: ])              ((,.~ (MQBS* _ _1)&rt)~            ;~ (MQBS* 0 _1)&}.@[)  (  ];.0~ ::(];.0~ _&(2:})) ((MQBS*1 _) ,:~ 2 # c@(0&({::))                                                      )) larfbrcfr 1 {:: ]
unmlqrcstep=:                          (0 {:: ])              ((,.  (MQBS* _  1)&rt)             ;  (MQBS* 0  1)&}.@])  (  ];.0~ ::(];.0~ _&(2:})) ((MQBS*1 _) ,:~ 2 # c@(0&({::))                                                      )) larfbrnfr 1 {:: ]

unmqllnstep=:                          (1 {:: ])              ((,   (MQBS* 1  _)&rt)             ;  (MQBS* 1  0)&}.@])~ ([ ];.0~ ::(];.0~ _&(3:})) ((MQBS*_ 1) ,:~ ((((((<:@-) {.)~ {.) ,  ((MQBS arounddown @ -~ {:)~ >:@{:))~ $)~ #&>))) larfblnbc 0 {:: ]
unmqllcstep=: (((<@(_1-(MQBS|<:@(c@[-#@(1 {:: ])))))`0:`]) }) (((1 {:: ]) , ~ ({.  ~ (0&({::)))) ;~ (}.  ~ (0&({::))))~ ([ ];.0~ ::(];.0~ _&(3:})) ((MQBS*_ 1) ,:~ ((((((<:@-) {.)~ {.) ,  ((MQBS arounddown @ -~ {:)~ >:@{:))~ $)~ #&>))) larfblcbc 0 {:: ]
unmqlrnstep=: (((<@(_1-(MQBS|<:@(c@[-c@(1 {:: ])))))`0:`]) }) (((1 {:: ]) ,.~ ({."1~ (0&({::)))) ;~ (}."1~ (0&({::))))~ ([ ];.0~ ::(];.0~ _&(3:})) ((MQBS*_ 1) ,:~ ((((((<:@-) {.)~ {.) ,  ((MQBS arounddown @ -~ {:)~ >:@{:))~ $)~ c&>))) larfbrnbc 0 {:: ]
unmqlrcstep=:                          (1 {:: ])              ((,.  (MQBS* _  1)&rt)             ;  (MQBS* 0  1)&}.@])~ ([ ];.0~ ::(];.0~ _&(3:})) ((MQBS*_ 1) ,:~ ((((((<:@-) {.)~ {.) ,  ((MQBS arounddown @ -~ {:)~ >:@{:))~ $)~ c&>))) larfbrcbc 0 {:: ]

unmqrlnstep=:                          (0 {:: ])              ((, ~ (MQBS*_1  _)&rt)~            ;~ (MQBS*_1  0)&}.@[)  (  ];.0~ ::(];.0~ _&(3:})) ((MQBS*_ 1) ,:~ 2 # #@(0&({::))                                                      )) larfblnfc 1 {:: ]
unmqrlcstep=:                          (0 {:: ])              ((,   (MQBS* 1  _)&rt)             ;  (MQBS* 1  0)&}.@])  (  ];.0~ ::(];.0~ _&(3:})) ((MQBS*_ 1) ,:~ 2 # #@(0&({::))                                                      )) larfblcfc 1 {:: ]
unmqrrnstep=:                          (0 {:: ])              ((,.  (MQBS* _  1)&rt)             ;  (MQBS* 0  1)&}.@])  (  ];.0~ ::(];.0~ _&(3:})) ((MQBS*_ 1) ,:~ 2 # c@(0&({::))                                                      )) larfbrnfc 1 {:: ]
unmqrrcstep=:                          (0 {:: ])              ((,.~ (MQBS* _ _1)&rt)~            ;~ (MQBS* 0 _1)&}.@[)  (  ];.0~ ::(];.0~ _&(3:})) ((MQBS*_ 1) ,:~ 2 # c@(0&({::))                                                      )) larfbrcfc 1 {:: ]

unmrqlnstep=: (((<@(_1-(MQBS|<:@(#@[-#@(1 {:: ])))))`0:`]) }) (((1 {:: ]) , ~ ({.  ~ (0&({::)))) ;~ (}.  ~ (0&({::))))~ ([ ];.0~ ::(];.0~ _&(2:})) ((MQBS*1 _) ,:~ ((((((<:@-) {:)~ {.) ,~ ((MQBS arounddown @ -~ {.)~ >:@{:))~ $)~ #&>))) larfblcbr 0 {:: ]
unmrqlcstep=:                          (1 {:: ])              ((,   (MQBS* 1  _)&rt)             ;  (MQBS* 1  0)&}.@])~ ([ ];.0~ ::(];.0~ _&(2:})) ((MQBS*1 _) ,:~ ((((((<:@-) {:)~ {.) ,~ ((MQBS arounddown @ -~ {.)~ >:@{:))~ $)~ #&>))) larfblnbr 0 {:: ]
unmrqrnstep=:                          (1 {:: ])              ((,.  (MQBS* _  1)&rt)             ;  (MQBS* 0  1)&}.@])~ ([ ];.0~ ::(];.0~ _&(2:})) ((MQBS*1 _) ,:~ ((((((<:@-) {:)~ {.) ,~ ((MQBS arounddown @ -~ {.)~ >:@{:))~ $)~ c&>))) larfbrcbr 0 {:: ]
unmrqrcstep=: (((<@(_1-(MQBS|<:@(#@[-c@(1 {:: ])))))`0:`]) }) (((1 {:: ]) ,.~ ({."1~ (0&({::)))) ;~ (}."1~ (0&({::))))~ ([ ];.0~ ::(];.0~ _&(2:})) ((MQBS*1 _) ,:~ ((((((<:@-) {:)~ {.) ,~ ((MQBS arounddown @ -~ {.)~ >:@{:))~ $)~ c&>))) larfbrnbr 0 {:: ]

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb     Action   Side   Tran  Syntax
NB. unmlqln  Q * C    left   none  QeC=.  Qf unmlqln C
NB. unmlqlc  Q^H * C  left   ct    cQeC=. Qf unmlqlc C
NB. unmlqrn  C * Q    right  none  eCQ=.  Qf unmlqlc C
NB. unmlqrc  C * Q^H  right  ct    eCcQ=. Qf unmlqlc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form Qf as returned by gelqf
NB.
NB. ############Algorithm:
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
NB. - implements LAPACK's xUNMLQ
NB. - unml2{lc,ln,rc,rn} and unmlq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmlqln=: _1  0 }. ((unml2ln`(, &>/@(  unmlqlnstep^:(>.@(MQBS %~ #@[)) (;~ 0  &{.                                             ) ))@.(MQBS < #@[)~  tru1            @({.  ~     0 _1&(ms $) ))~ ,   &0)
unmlqlc=: _1  0 }. ((unml2lc`(, &>/@([ unmlqlcstep^:(>.@(MQBS %~ #@[)) (( {.       ;   }.      )~ (MQBS arounddown @ (_1 + #)))~))@.(MQBS < #@[)~  tru1            @({.  ~     0 _1&(ms $) ))~ ,   &0)
unmlqrn=:  0 _1 }. ((unml2rn`(,.&>/@([ unmlqrnstep^:(>.@(MQBS %~ #@[)) ((({.~ _&,) ;  (}.~ 0&,))  (MQBS arounddown @ (_1 + #)))~))@.(MQBS < #@[)~  tru1            @({.  ~     0 _1&(ms $) ))~ ,.  &0)
unmlqrc=:  0 _1 }. ((unml2rc`(,.&>/@(  unmlqrcstep^:(>.@(MQBS %~ #@[)) (;~ _ 0&{.                                             ) ))@.(MQBS < #@[)~  tru1            @({.  ~     0 _1&(ms $) ))~ ,.  &0)

NB. ---------------------------------------------------------
NB. Verb     Action   Side   Tran  Syntax
NB. unmqlln  Q * C    left   none  QeC=.  Qf unmqlln C
NB. unmqllc  Q^H * C  left   ct    cQeC=. Qf unmqllc C
NB. unmqlrn  C * Q    right  none  eCQ=.  Qf unmqllc C
NB. unmqlrc  C * Q^H  right  ct    eCcQ=. Qf unmqllc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form Qf as returned by geqlf
NB.
NB. ############Algorithm:
NB.   iters==⌊(#+BS-1)/BS⌋
NB.
NB. If:
NB.   LQf=. gelqf A
NB.   L=. trl (0 _1 }. LQf)
NB.   Qf=. tru1 LQf
NB.   Q=. unglq LQf
NB.   k=. <./ $ A
NB. then
NB.   ((idmat @ c) (-: clean) ((_1 & tru1 @ ((_ _10)&{.)) unmqllc (     ungql)) @ geqlf) A
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - implements LAPACK's xUNMLQ
NB. - unml2{lc,ln,rc,rn} and unmlq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmqlln=:  1  0 }. ((unm2lln`(, &>/@([ unmqllnstep^:(>.@(MQBS %~ c@[)) (( {.       ;~  }.      )~ (MQBS                  - c ))~))@.(MQBS < c@[)~ (tru1~ (-~/ @ $))@({."1~ -@(_1  0&(ms $))))~ , ~ &0)
unmqllc=:  1  0 }. ((unm2llc`(, &>/@(  unmqllcstep^:(>.@(MQBS %~ c@[)) (;  0  &{.                                             ) ))@.(MQBS < c@[)~ (tru1~ (-~/ @ $))@({."1~ -@(_1  0&(ms $))))~ , ~ &0)
unmqlrn=:  0  1 }. ((unm2lrn`(,.&>/@(  unmqlrnstep^:(>.@(MQBS %~ c@[)) (;  _ 0&{.                                             ) ))@.(MQBS < c@[)~ (tru1~ (-~/ @ $))@({."1~ -@(_1  0&(ms $))))~ ,.~ &0)
unmqlrc=:  0  1 }. ((unm2lrc`(,.&>/@([ unmqlrcstep^:(>.@(MQBS %~ c@[)) ((({.~ _&,) ;~ (}.~ 0&,))  (MQBS                  - c ))~))@.(MQBS < c@[)~ (tru1~ (-~/ @ $))@({."1~ -@(_1  0&(ms $))))~ ,.~ &0)

NB. ---------------------------------------------------------
NB. Verb     Action   Side   Tran  Syntax
NB. unmqrln  Q * C    left   none  QeC=.  Qf unmqrln C
NB. unmqrlc  Q^H * C  left   ct    cQeC=. Qf unmqrlc C
NB. unmqrrn  C * Q    right  none  eCQ=.  Qf unmqrlc C
NB. unmqrrc  C * Q^H  right  ct    eCcQ=. Qf unmqrlc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form Qf as returned by geqrf
NB.
NB. ############Algorithm:
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
NB. - implements LAPACK's xUNMLQ
NB. - unml2{lc,ln,rc,rn} and unmlq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmqrln=: _1  0 }. ((unm2rln`(, &>/@([ unmqrlnstep^:(>.@(MQBS %~ c@[)) (( {.       ;   }.      )~ (MQBS arounddown @ (_1 + c)))~))@.(MQBS < c@[)~  trl1            @({."1~    _1  0&(ms $) ))~ ,   &0)
unmqrlc=: _1  0 }. ((unm2rlc`(, &>/@(  unmqrlcstep^:(>.@(MQBS %~ c@[)) (;~ 0  &{.                                             ) ))@.(MQBS < c@[)~  trl1            @({."1~    _1  0&(ms $) ))~ ,   &0)
unmqrrn=:  0 _1 }. ((unm2rrn`(,.&>/@(  unmqrrnstep^:(>.@(MQBS %~ c@[)) (;~ _ 0&{.                                             ) ))@.(MQBS < c@[)~  trl1            @({."1~    _1  0&(ms $) ))~ ,.  &0)
unmqrrc=:  0 _1 }. ((unm2rrc`(,.&>/@([ unmqrrcstep^:(>.@(MQBS %~ c@[)) ((({.~ _&,) ;  (}.~ 0&,))  (MQBS arounddown @ (_1 + c)))~))@.(MQBS < c@[)~  trl1            @({."1~    _1  0&(ms $) ))~ ,.  &0)

NB. ---------------------------------------------------------
NB. Verb     Action   Side   Tran  Syntax
NB. unmrqln  Q * C    left   none  QeC=.  Qf unmrqln C
NB. unmrqlc  Q^H * C  left   ct    cQeC=. Qf unmrqlc C
NB. unmrqrn  C * Q    right  none  eCQ=.  Qf unmrqlc C
NB. unmrqrc  C * Q^H  right  ct    eCcQ=. Qf unmrqlc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form Qf as returned by gerqf
NB.
NB. ############Algorithm:
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
NB. - implements LAPACK's xUNMLQ
NB. - unml2{lc,ln,rc,rn} and unmlq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmrqln=:  1  0 }. ((unmr2ln`(, &>/@(  unmrqlnstep^:(>.@(MQBS %~ #@[)) (;  0  &{.                                             ) ))@.(MQBS < #@[)~ (trl1~ (-~/ @ $))@({.  ~ -@( 0 _1&(ms $))))~ , ~ &0)
unmrqlc=:  1  0 }. ((unmr2lc`(, &>/@([ unmrqlcstep^:(>.@(MQBS %~ #@[)) (( {.       ;~  }.      )~ (MQBS                  - # ))~))@.(MQBS < #@[)~ (trl1~ (-~/ @ $))@({.  ~ -@( 0 _1&(ms $))))~ , ~ &0)
unmrqrn=:  0  1 }. ((unmr2rn`(,.&>/@([ unmrqrnstep^:(>.@(MQBS %~ #@[)) ((({.~ _&,) ;~ (}.~ 0&,))  (MQBS                  - # ))~))@.(MQBS < #@[)~ (trl1~ (-~/ @ $))@({.  ~ -@( 0 _1&(ms $))))~ ,.~ &0)
unmrqrc=:  0  1 }. ((unmr2rc`(,.&>/@(  unmrqrcstep^:(>.@(MQBS %~ #@[)) (;  _ 0&{.                                             ) ))@.(MQBS < #@[)~ (trl1~ (-~/ @ $))@({.  ~ -@( 0 _1&(ms $))))~ ,.~ &0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmhrlln  Q * C    left   none  QeC=.  Qf unmhrlln C
NB. unmhrllc  Q^H * C  left   ct    cQeC=. Qf unmhrllc C
NB. unmhrlrn  C * Q    right  none  eCQ=.  Qf unmhrlrn C
NB. unmhrlrc  C * Q^H  right  ct    eCcQ=. Qf unmhrlrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by unitary (orthogonal)
NB.   matrix Q, which is represented in factored form Qf as
NB.   returned by gehrdl
NB.
NB. ###########Algorithm:
NB.   iters==⌊(#+BS-1)/BS⌋
NB.
NB. If:
NB.   LQf=. gelqf A
NB.   L=. trl (0 _1 }. LQf)
NB.   Qf=. tru1 LQf
NB.   Q=. unglq LQf
NB.   k=. <./ $ A
NB. then
NB.   
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - implements LAPACK's xUNMLQ
NB. - unml2{lc,ln,rc,rn} and unmlq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmhrlln=: (unmlqln~ (|. !. 0))~
unmhrllc=: (unmlqlc~ (|. !. 0))~
unmhrlrn=: (unmlqrn~ (|. !. 0))~
unmhrlrc=: (unmlqrc~ (|. !. 0))~

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmhruln  Q * C    left   none  QeC=.  Qf unmhruln C
NB. unmhrulc  Q^H * C  left   ct    cQeC=. Qf unmhrulc C
NB. unmhrurn  C * Q    right  none  eCQ=.  Qf unmhrurn C
NB. unmhrurc  C * Q^H  right  ct    eCcQ=. Qf unmhrurc C
NB.
NB. Description:
NB.   Multiply a general matrix C by unitary (orthogonal)
NB.   matrix Q, which is represented in factored form Qf as
NB.   returned by gehrdu
NB.
NB. ###########Algorithm:
NB.   iters==⌊(#+BS-1)/BS⌋
NB.
NB. If:
NB.   LQf=. gelqf A
NB.   L=. trl (0 _1 }. LQf)
NB.   Qf=. tru1 LQf
NB.   Q=. unglq LQf
NB.   k=. <./ $ A
NB. then
NB.   
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - implements LAPACK's xUNMLQ
NB. - unml2{lc,ln,rc,rn} and unmlq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmhruln=: (unmqrln~ (0 _1 & (|. !. 0)))~
unmhrulc=: (unmqrlc~ (0 _1 & (|. !. 0)))~
unmhrurn=: (unmqrrn~ (0 _1 & (|. !. 0)))~
unmhrurc=: (unmqrrc~ (0 _1 & (|. !. 0)))~

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testunmq
NB.
NB. Description:
NB.   Test Q multiplication qf-algorithms by general matrix
NB.   given
NB.
NB. Syntax:
NB.   testunmq (A;C)
NB. where
NB.   A - m×n-matrix, is used to get Qf
NB.   C - m×n-matrix, is used as multiplier

testunmq=: 3 : 0
  'A C'=. y
  rcond=. ((_."_)`(norm1 con (getri@getrf)) @. (=/@$)) C  NB. meaninigful for square matrices only
  'LQf QfL QfR RQf'=. xQf=. (gelqf ; geqlf ; geqrf ; gerqf) A
  'Qlq Qql Qqr Qrq'=. (((unglq~ (<:@c))&.>)`((ungql~ (<:@#))&.>)`((ungqr~ (<:@#))&.>)`((ungrq~ (<:@c))&.>)) ag xQf

  ('unmlqln' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (LQf;(ct C);    Qlq )  NB. berr := ||Q * C   - Q * C  || / (ε * ||C|| * n)
  ('unmlqlc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (LQf;(ct C);(ct Qlq))  NB. berr := ||Q^H * C - Q^H * C|| / (ε * ||C|| * n)
  ('unmlqrn' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (LQf;    C ;    Qlq )  NB. berr := ||C * Q   - C * Q  || / (ε * ||C|| * n)
  ('unmlqrc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (LQf;    C ;(ct Qlq))  NB. berr := ||C * Q^H - C * Q^H|| / (ε * ||C|| * n)

  ('unmqlln' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (QfL;    C ;    Qql )  NB. berr := ||Q * C   - Q * C  || / (ε * ||C|| * m)
  ('unmqllc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (QfL;    C ;(ct Qql))  NB. berr := ||Q^H * C - Q^H * C|| / (ε * ||C|| * m)
  ('unmqlrn' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (QfL;(ct C);    Qql )  NB. berr := ||C * Q - C * Q    || / (ε * ||C|| * m)
  ('unmqlrc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (QfL;(ct C);(ct Qql))  NB. berr := ||C * Q^H - C * Q^H|| / (ε * ||C|| * m)

  ('unmqrln' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (QfR;    C ;    Qqr )  NB. berr := ||Q * C   - Q * C  || / (ε * ||C|| * m)
  ('unmqrlc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (QfR;    C ;(ct Qqr))  NB. berr := ||Q^H * C - Q^H * C|| / (ε * ||C|| * m)
  ('unmqrrn' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (QfR;(ct C);    Qqr )  NB. berr := ||C * Q   - C * Q  || / (ε * ||C|| * m)
  ('unmqrrc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (QfR;(ct C);(ct Qqr))  NB. berr := ||C * Q^H - C * Q^H|| / (ε * ||C|| * m)

  ('unmrqln' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (RQf;(ct C);    Qrq )  NB. berr := ||Q * C   - Q * C  || / (ε * ||C|| * n)
  ('unmrqlc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (RQf;(ct C);(ct Qrq))  NB. berr := ||Q^H * C - Q^H * C|| / (ε * ||C|| * n)
  ('unmrqrn' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (RQf;    C ;    Qrq )  NB. berr := ||C * Q   - C * Q  || / (ε * ||C|| * n)
  ('unmrqrc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (RQf;    C ;(ct Qrq))  NB. berr := ||C * Q^H - C * Q^H|| / (ε * ||C|| * n)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testunmhr
NB.
NB. Description:
NB.   Test Q multiplication hrd-algorithms by square general
NB.   matrix given
NB.
NB. Syntax:
NB.   testunmhr (A;C)
NB. where
NB.   A - n×n-matrix, is used to get Qf
NB.   C - n×n-matrix, is used as multiplier

testunmhr=: 3 : 0
  'A C'=. y
  rcond=. (norm1 con (getri@getrf)) C
  'HlQf HuQf'=. xQf=. ((gehrdl~ (0,#)) ; (gehrdu~ (0,c))) A
  'Qhrl Qhru'=. ((unghrl&.>)`(unghru&.>)) ag xQf

  ('unmhrlln' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (HlQf;(ct C);    Qhrl )  NB. berr := ||Q * C   - Q * C  || / (ε * ||C|| * n)
  ('unmhrllc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (HlQf;(ct C);(ct Qhrl))  NB. berr := ||Q^H * C - Q^H * C|| / (ε * ||C|| * n)
  ('unmhrlrn' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (HlQf;    C ;    Qhrl )  NB. berr := ||C * Q   - C * Q  || / (ε * ||C|| * n)
  ('unmhrlrc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*c)@(1 {:: [))))))) (HlQf;    C ;(ct Qhrl))  NB. berr := ||C * Q^H - C * Q^H|| / (ε * ||C|| * n)

  ('unmhruln' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (HuQf;(ct C);    Qhru )  NB. berr := ||Q * C   - Q * C  || / (ε * ||C|| * m)
  ('unmhrulc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp~ & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (HuQf;(ct C);(ct Qhru))  NB. berr := ||Q^H * C - Q^H * C|| / (ε * ||C|| * m)
  ('unmhrurn' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (HuQf;    C ;    Qhru )  NB. berr := ||C * Q   - C * Q  || / (ε * ||C|| * m)
  ('unmhrurc' tdyad ((0 & {::)`(1 & {::)`]`(rcond"_)`(_."_)`(((norm1@((- ((mp  & >/)@}.))~)) % (FP_EPS*((norm1*#)@(1 {:: [))))))) (HuQf;    C ;(ct Qhru))  NB. berr := ||C * Q^H - C * Q^H|| / (ε * ||C|| * m)

  EMPTY
)

NB. ---------------------------------------------------------
NB. testmq
NB.
NB. Description:
NB.   Adv. to make verb to test unmxxxxx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   vtest=. mkmat testmq
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   vtest - monad to test algorithms by matrix mat; is
NB.           called as:
NB.             vtest (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     (? @ $ 0:) testmq_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     (_1 1 0 4 _6 4 & gemat_mt_) testmq_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testmq_mt_ 150 200

testmq=: 1 : 'EMPTY_mt_ [ (testunmhr_mt_ ^: (=/@$@(0&({::))) [ testunmq_mt_) @ (u ; u)'
