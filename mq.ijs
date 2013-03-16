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
NB. unmlzxx    Multiply a general matrix by a matrix with
NB.            orthonormal rows, which is represented in
NB.            factored form, as returned by tzlzf
NB. unmzlxx    Multiply a general matrix by a matrix with
NB.            orthonormal columns, which is represented in
NB.            factored form, as returned by tzzlf
NB. unmzrxx    Multiply a general matrix by a matrix with
NB.            orthonormal columns, which is represented in
NB.            factored form, as returned by tzzrf
NB. unmrzxx    Multiply a general matrix by a matrix with
NB.            orthonormal rows, which is represented in
NB.            factored form, as returned by tzrzf
NB.
NB. testunmq   Test unmxxxx by general matrix given
NB. testunmhr  Test unmhrxxx by square matrix given
NB. testunmz   Test unmxxxx by trapezoidal matrix given
NB. testmq     Adv. to make verb to test unmxxxxx by matrix
NB.            of generator and shape given
NB.
NB. Version: 0.9.0 2012-12-29
NB.
NB. Copyright 2010-2012 Igor Zhuravlov
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
NB. Blocked code constants

MQNB=: 32   NB. block size limit

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
NB.   1) form rios, rIOS of vtau[i] which defines an
NB.      elementary reflector, depending on aC(i)'s size
NB.   2) extract vtau[i] from Qf and ravel it:
NB.        vtaui=. rios (, ;. 0) Qf
NB.   3) apply an elementary reflector defined by vtau[i] to
NB.      either pfxC(i) or sfxC(i) to produce tmp:
NB.        tmp=. vtaui larfxxxx xfxCi
NB.   4) combine tmp and either pfxC(i) or sfxC(i) to
NB.      produce pfxC(i+1) and sfxC(i+1)

unml2lnstep=: (0 {:: ]) ((,    1  _&rt)  ;   }.   @])  (,;.0~ (1 _ ,:~ 2 #      #@(0&{::))) larflcfr 1 {:: ]
unml2lcstep=: (0 {:: ]) ((, ~ _1  _&rt)~ ;~  }:   @[)  (,;.0~ (1 _ ,:~ 2 #      #@(0&{::))) larflnfr 1 {:: ]
unml2rnstep=: (0 {:: ]) ((,.~  _ _1&rt)~ ;~ (}:"1)@[)  (,;.0~ (1 _ ,:~ 2 #      c@(0&{::))) larfrcfr 1 {:: ]
unml2rcstep=: (0 {:: ]) ((,.   _  1&rt)  ;  (}."1)@])  (,;.0~ (1 _ ,:~ 2 #      c@(0&{::))) larfrnfr 1 {:: ]

unm2llnstep=: (1 {:: ]) ((,    1  _&rt)  ;   }.   @])~ (,;.0~ (_ 1 ,:~ 2 # _1 - #@(1&{::))) larflnbc 0 {:: ]
unm2llcstep=: (1 {:: ]) ((, ~ _1  _&rt)~ ;~  }:   @[)~ (,;.0~ (_ 1 ,:~ 2 # _1 - #@(1&{::))) larflcbc 0 {:: ]
unm2lrnstep=: (1 {:: ]) ((,.~  _ _1&rt)~ ;~ (}:"1)@[)~ (,;.0~ (_ 1 ,:~ 2 # _1 - c@(1&{::))) larfrnbc 0 {:: ]
unm2lrcstep=: (1 {:: ]) ((,.   _  1&rt)  ;  (}."1)@])~ (,;.0~ (_ 1 ,:~ 2 # _1 - c@(1&{::))) larfrcbc 0 {:: ]

unm2rlnstep=: (0 {:: ]) ((, ~ _1  _&rt)~ ;~  }:   @[)  (,;.0~ (_ 1 ,:~ 2 #      #@(0&{::))) larflnfc 1 {:: ]
unm2rlcstep=: (0 {:: ]) ((,    1  _&rt)  ;   }.   @])  (,;.0~ (_ 1 ,:~ 2 #      #@(0&{::))) larflcfc 1 {:: ]
unm2rrnstep=: (0 {:: ]) ((,.   _  1&rt)  ;  (}."1)@])  (,;.0~ (_ 1 ,:~ 2 #      c@(0&{::))) larfrnfc 1 {:: ]
unm2rrcstep=: (0 {:: ]) ((,.~  _ _1&rt)~ ;~ (}:"1)@[)  (,;.0~ (_ 1 ,:~ 2 #      c@(0&{::))) larfrcfc 1 {:: ]

unmr2lnstep=: (1 {:: ]) ((, ~ _1  _&rt)~ ;~  }:   @[)~ (,;.0~ (1 _ ,:~ 2 # _1 - #@(1&{::))) larflcbr 0 {:: ]
unmr2lcstep=: (1 {:: ]) ((,    1  _&rt)  ;   }.   @])~ (,;.0~ (1 _ ,:~ 2 # _1 - #@(1&{::))) larflnbr 0 {:: ]
unmr2rnstep=: (1 {:: ]) ((,.   _  1&rt)  ;  (}."1)@])~ (,;.0~ (1 _ ,:~ 2 # _1 - c@(1&{::))) larfrcbr 0 {:: ]
unmr2rcstep=: (1 {:: ]) ((,.~  _ _1&rt)~ ;~ (}:"1)@[)~ (,;.0~ (1 _ ,:~ 2 # _1 - c@(1&{::))) larfrnbr 0 {:: ]

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
NB. Assertions:
NB.   (idmat n) (-: clean@: }:   ) ((   tru1 gelqf A) unml2ln ((ct unglq gelqf A) ,   0))
NB.   (idmat n) (-: clean@: }:   ) ((   tru1 gelqf A) unml2lc ((   unglq gelqf A) ,   0))
NB.   (idmat n) (-: clean@:(}:"1)) ((   tru1 gelqf A) unml2rn ((ct unglq gelqf A) ,.  0))
NB.   (idmat n) (-: clean@:(}:"1)) ((   tru1 gelqf A) unml2rc ((   unglq gelqf A) ,.  0))
NB.   (idmat n) (-: clean@: }.   ) ((_1 tru1 geqlf A) unm2lln ((ct ungql geqlf A) , ~ 0))
NB.   (idmat n) (-: clean@: }.   ) ((_1 tru1 geqlf A) unm2llc ((   ungql geqlf A) , ~ 0))
NB.   (idmat n) (-: clean@:(}."1)) ((_1 tru1 geqlf A) unm2lrn ((ct ungql geqlf A) ,.~ 0))
NB.   (idmat n) (-: clean@:(}."1)) ((_1 tru1 geqlf A) unm2lrc ((   ungql geqlf A) ,.~ 0))
NB.   (idmat n) (-: clean@: }:   ) ((   trl1 geqrf A) unm2rln ((ct ungqr geqrf A) ,   0))
NB.   (idmat n) (-: clean@: }:   ) ((   trl1 geqrf A) unm2rlc ((   ungqr geqrf A) ,   0))
NB.   (idmat n) (-: clean@:(}:"1)) ((   trl1 geqrf A) unm2rrn ((ct ungqr geqrf A) ,.  0))
NB.   (idmat n) (-: clean@:(}:"1)) ((   trl1 geqrf A) unm2rrc ((   ungqr geqrf A) ,.  0))
NB.   (idmat n) (-: clean@: }.   ) (( 1 trl1 gerqf A) unmr2ln ((ct ungrq gerqf A) , ~ 0))
NB.   (idmat n) (-: clean@: }.   ) (( 1 trl1 gerqf A) unmr2lc ((   ungrq gerqf A) , ~ 0))
NB.   (idmat n) (-: clean@:(}."1)) (( 1 trl1 gerqf A) unmr2rn ((ct ungrq gerqf A) ,.~ 0))
NB.   (idmat n) (-: clean@:(}."1)) (( 1 trl1 gerqf A) unmr2rc ((   ungrq gerqf A) ,.~ 0))
NB. where
NB.   2 -: # $ A  NB. A is a 2-rank array (i.e. matrix)
NB.   -:/@$ A     NB. A is a square matrix (it's not
NB.               NB.   necessary and is assumed just
NB.               NB.   to simplify assertions)
NB.   n=. # A     NB. size of matrix A
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - implements LAPACK's {DOR,ZUN}M{L2,2L,2R,R2}
NB. - unml2{ln,lc,rn,rc} and unmlq{ln,lc,rn,rc} respectively
NB.   are topologic equivalents

unml2ln_8=: (larflcfr&:>/@,~ (-@c <\ ,)@|.   )~ <
unml2lc_8=: (larflnfr&:>/@,~ (-@c <\ ,)      )~ <
unml2rn_8=: (larfrcfr&:>/@,~ (-@c <\ ,)      )~ <
unml2rc_8=: (larfrnfr&:>/@,~ (-@c <\ ,)@|.   )~ <

unm2lln=: , &>/@([ unm2llnstep^:(c@[) (( {.       ;~  }.      )~ ( 1 - c))~)
unm2llc=: , &>/@(  unm2llcstep^:(c@[) (;  0& {.                          ) )
unm2lrn=: ,.&>/@(  unm2lrnstep^:(c@[) (;  0&({."1)                       ) )
unm2lrc=: ,.&>/@([ unm2lrcstep^:(c@[) ((({.~ _&,) ;~ (}.~ 0&,))  ( 1 - c))~)

unm2lln_8=: (larflnbc&:>/@,~ (-@c <\ ,)@|.@|:)~ <
unm2llc_8=: (larflcbc&:>/@,~ (-@c <\ ,)   @|:)~ <
unm2lrn_8=: (larfrnbc&:>/@,~ (-@c <\ ,)   @|:)~ <
unm2lrc_8=: (larfrcbc&:>/@,~ (-@c <\ ,)@|.@|:)~ <

unm2rln=: , &>/@([ unm2rlnstep^:(c@[) (( {.       ;   }.      )~ (_1 + c))~)
unm2rlc=: , &>/@(  unm2rlcstep^:(c@[) (;~ 0& {.                          ) )
unm2rrn=: ,.&>/@(  unm2rrnstep^:(c@[) (;~ 0&({."1)                       ) )
unm2rrc=: ,.&>/@([ unm2rrcstep^:(c@[) ((({.~ _&,) ;  (}.~ 0&,))  (_1 + c))~)

unm2rln_8=: (larflnfc&:>/@,~ (-@c <\ ,)   @|:)~ <
unm2rlc_8=: (larflcfc&:>/@,~ (-@c <\ ,)@|.@|:)~ <
unm2rrn_8=: (larfrnfc&:>/@,~ (-@c <\ ,)@|.@|:)~ <
unm2rrc_8=: (larfrcfc&:>/@,~ (-@c <\ ,)   @|:)~ <

unmr2ln=: , &>/@(  unmr2lnstep^:(#@[) (;  0& {.                          ) )
unmr2lc=: , &>/@([ unmr2lcstep^:(#@[) (( {.       ;~  }.      )~ ( 1 - #))~)
unmr2rn=: ,.&>/@([ unmr2rnstep^:(#@[) ((({.~ _&,) ;~ (}.~ 0&,))  ( 1 - #))~)
unmr2rc=: ,.&>/@(  unmr2rcstep^:(#@[) (;  0&({."1)                       ) )

unmr2ln_8=: (larflcbr&:>/@,~ (-@c <\ ,)      )~ <
unmr2lc_8=: (larflnbr&:>/@,~ (-@c <\ ,)@|.   )~ <
unmr2rn_8=: (larfrcbr&:>/@,~ (-@c <\ ,)@|.   )~ <
unmr2rc_8=: (larfrnbr&:>/@,~ (-@c <\ ,)      )~ <

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmlqln   Q * C    left   none  B=. LQf unmlqln C
NB. unmlqlc   Q^H * C  left   ct    B=. LQf unmlqlc C
NB. unmlqrn   C * Q    right  none  B=. LQf unmlqrn C
NB. unmlqrc   C * Q^H  right  ct    B=. LQf unmlqrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form LQf as returned by gelqf
NB. where
NB.   B,C - m×n-matrices
NB.   LQf - n×(m+1)-matrix (ln,lc cases) or m×(n+1)-matrix
NB.         (rn,rc), contains Qf (unit diagonal not stored),
NB.         the output of gelqf
NB.   Qf  - k×(m+1)-matrix (ln,lc) or k×(n+1)-matrix (rn,rc),
NB.         unit upper triangular, the Q represented in
NB.         factored form
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as a product of k
NB.         elementary reflectors:
NB.           Q = Π{H(i)',i=k-1:0}
NB.   k   = min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (idmat@c (-: clean) (unmlqln ct@(<:@c unglq ]))@gelqf) A
NB.   (idmat@c (-: clean) (unmlqlc    (<:@c unglq ]))@gelqf) A
NB.   (idmat@c (-: clean) (unmlqrn ct@(<:@c unglq ]))@gelqf) A
NB.   (idmat@c (-: clean) (unmlqrc    (<:@c unglq ]))@gelqf) A
NB.
NB. Notes:
NB. - implements LAPACK's DORMLQ, ZUNMLQ
NB. - unml2{lc,ln,rc,rn} and unmlq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmlqln=: }:  @(((unml2ln`((larfblcfr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,  &0)
unmlqlc=: }:  @(((unml2lc`((larfblnfr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,  &0)
unmlqrn=: }:"1@(((unml2rn`((larfbrcfr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,. &0)
unmlqrc=: }:"1@(((unml2rc`((larfbrnfr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,. &0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmqlln   Q * C    left   none  B=. QfL unmqlln C
NB. unmqllc   Q^H * C  left   ct    B=. QfL unmqllc C
NB. unmqlrn   C * Q    right  none  B=. QfL unmqlrn C
NB. unmqlrc   C * Q^H  right  ct    B=. QfL unmqlrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form QfL as returned by geqlf
NB. where
NB.   B,C - m×n-matrices
NB.   QfL - (m+1)×n-matrix (ln,lc cases), (n+1)×m-matrix
NB.         (rn,rc), contains Qf (unit diagonal not stored),
NB.         the output of geqlf
NB.   Qf  - (m+1)×k-matrix (ln,lc) or (n+1)×k-matrix (rn,rc),
NB.         unit upper triangular, the Q represented in
NB.         factored form
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as a product of k
NB.         elementary reflectors:
NB.           Q = Π{H(i),i=k-1:0}
NB.   k   = min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (idmat@# (-: clean) (unmqlln ct@(<:@# ungql ]))@geqlf) A
NB.   (idmat@# (-: clean) (unmqllc    (<:@# ungql ]))@geqlf) A
NB.   (idmat@# (-: clean) (unmqlrn ct@(<:@# ungql ]))@geqlf) A
NB.   (idmat@# (-: clean) (unmqlrc    (<:@# ungql ]))@geqlf) A
NB.
NB. Notes:
NB. - implements LAPACK's DORMQL, ZUNMQL
NB. - unm2l{lc,ln,rc,rn} and unmql{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmqlln=: }.  @(((unm2lln`((larfblnbc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ , ~&0)
unmqllc=: }.  @(((unm2llc`((larfblcbc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ , ~&0)
unmqlrn=: }."1@(((unm2lrn`((larfbrnbc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ ,.~&0)
unmqlrc=: }."1@(((unm2lrc`((larfbrcbc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ ,.~&0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmqrln   Q * C    left   none  B=. QfR unmqrln C
NB. unmqrlc   Q^H * C  left   ct    B=. QfR unmqrlc C
NB. unmqrrn   C * Q    right  none  B=. QfR unmqrrn C
NB. unmqrrc   C * Q^H  right  ct    B=. QfR unmqrrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form QfR as returned by geqrf
NB. where
NB.   B,C - m×n-matrices
NB.   QfR - (m+1)×n-matrix (ln,lc cases), (n+1)×m-matrix
NB.         (rn,rc), contains Qf (unit diagonal not stored),
NB.         the output of geqrf
NB.   Qf  - (m+1)×k-matrix (ln,lc) or (n+1)×k-matrix (rn,rc),
NB.         unit lower triangular, the Q represented in
NB.         factored form
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as a product of k
NB.         elementary reflectors:
NB.           Q = Π{H(i),i=0:k-1}
NB.   k   = min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (idmat@# (-: clean) (unmqrln ct@(<:@# ungqr ]))@geqrf) A
NB.   (idmat@# (-: clean) (unmqrlc    (<:@# ungqr ]))@geqrf) A
NB.   (idmat@# (-: clean) (unmqrrn ct@(<:@# ungqr ]))@geqrf) A
NB.   (idmat@# (-: clean) (unmqrrc    (<:@# ungqr ]))@geqrf) A
NB.
NB. Notes:
NB. - implements LAPACK's DORMQR, ZUNMQR
NB. - unm2r{lc,ln,rc,rn} and unmqr{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmqrln=: }:  @(((unm2rln`((larfblnfc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,  &0)
unmqrlc=: }:  @(((unm2rlc`((larfblcfc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,  &0)
unmqrrn=: }:"1@(((unm2rrn`((larfbrnfc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,. &0)
unmqrrc=: }:"1@(((unm2rrc`((larfbrcfc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,. &0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmrqln   Q * C    left   none  B=. RQf unmrqln C
NB. unmrqlc   Q^H * C  left   ct    B=. RQf unmrqlc C
NB. unmrqrn   C * Q    right  none  B=. RQf unmrqrn C
NB. unmrqrc   C * Q^H  right  ct    B=. RQf unmrqrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form RQf as returned by gerqf
NB.
NB. where
NB.   B,C - m×n-matrices
NB.   LQf - n×(m+1)-matrix (ln,lc cases) or m×(n+1)-matrix
NB.         (rn,rc), contains Qf (unit diagonal not stored),
NB.         the output of gerqf
NB.   Qf  - k×(m+1)-matrix (ln,lc) or k×(n+1)-matrix (rn,rc),
NB.         unit lower triangular, the Q represented in
NB.         factored form
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as a product of k
NB.         elementary reflectors:
NB.           Q = Π{H(i)',i=0:k-1}
NB.   k   = min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (idmat@c (-: clean) (unmrqln ct@(<:@c ungrq ]))@gerqf) A
NB.   (idmat@c (-: clean) (unmrqlc    (<:@c ungrq ]))@gerqf) A
NB.   (idmat@c (-: clean) (unmrqrn ct@(<:@c ungrq ]))@gerqf) A
NB.   (idmat@c (-: clean) (unmrqrc    (<:@c ungrq ]))@gerqf) A
NB.
NB. Notes:
NB. - implements LAPACK's DORMRQ, ZUNMRQ
NB. - unmr2{lc,ln,rc,rn} and unmrq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmrqln=: }.  @(((unmr2ln`((larfblcbr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ , ~&0)
unmrqlc=: }.  @(((unmr2lc`((larfblnbr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ , ~&0)
unmrqrn=: }."1@(((unmr2rn`((larfbrcbr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ ,.~&0)
unmrqrc=: }."1@(((unmr2rc`((larfbrnbr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ ,.~&0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmhrlln  Q * C    left   none  B=. HQf unmhrlln C
NB. unmhrllc  Q^H * C  left   ct    B=. HQf unmhrllc C
NB. unmhrlrn  C * Q    right  none  B=. HQf unmhrlrn C
NB. unmhrlrc  C * Q^H  right  ct    B=. HQf unmhrlrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by unitary (orthogonal)
NB.   matrix Q, which is represented in factored form HQf as
NB.   returned by gehrdl
NB. where
NB.   B,C - m×n-matrices
NB.   HQf - m×(m+1)-matrix (ln,lc cases) or n×(n+1)-matrix
NB.         (rn,rc), contains Qf (unit diagonal not stored),
NB.         the output of gehrdl
NB.   Qf  - (s-1)×(m-h)-matrix (ln,lc) or (s-1)×(n-h)-matrix
NB.         (rn,rc), unit upper triangular, the Q represented
NB.         in factored form, located in HQf[h:h+s-2,h+1:end]
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), being
NB.         unit matrix with unitary (orthogonal) matrix
NB.         inserted into elements Q[h:h+s-1,h:h+s-1] :
NB.           Q = Π{H(i)',i=h+s-2:h}
NB.   hs  - 2-vector of integers (h,s) 'head' and 'size',
NB.         defines submatrix Qf position in matrix HQf, see
NB.         see gehrdl
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (idmat@c (-: clean) (unmhrlln ct@unghrl)@gehrdl) A
NB.   (idmat@c (-: clean) (unmhrllc ct@unghrl)@gehrdl) A
NB.   (idmat@c (-: clean) (unmhrlrn ct@unghrl)@gehrdl) A
NB.   (idmat@c (-: clean) (unmhrlrc ct@unghrl)@gehrdl) A
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - instead of using f and s parameters, the following
NB.   product is really calculating:
NB.     Q = Π{H(i)',i=n-1:0} ,
NB.   where
NB.     H(0:f-1) = H(f+s-1:n-1) = I .
NB.   This approach delivers excessive calculations in rare
NB.   case ((f>0) OR (f+s<n)).

unmhrlln=: (unmlqln~ |.!.0)~
unmhrllc=: (unmlqlc~ |.!.0)~
unmhrlrn=: (unmlqrn~ |.!.0)~
unmhrlrc=: (unmlqrc~ |.!.0)~

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmhruln  Q * C    left   none  B=. HQf unmhruln C
NB. unmhrulc  Q^H * C  left   ct    B=. HQf unmhrulc C
NB. unmhrurn  C * Q    right  none  B=. HQf unmhrurn C
NB. unmhrurc  C * Q^H  right  ct    B=. HQf unmhrurc C
NB.
NB. Description:
NB.   Multiply a general matrix C by unitary (orthogonal)
NB.   matrix Q, which is represented in factored form HQf as
NB.   returned by gehrdu
NB. where
NB.   B,C - m×n-matrices
NB.   HQf - (m+1)×m-matrix (ln,lc cases) or (n+1)×n-matrix
NB.         (rn,rc), contains Qf (unit diagonal not stored),
NB.         the output of gehrdu
NB.   Qf  - (m-h)×(s-1)-matrix (ln,lc) or (n-h)×(s-1)-matrix
NB.         (rn,rc), unit lower triangular, the Q represented
NB.         in factored form, located in HQf[h+1:end,h:h+s-2]
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), being
NB.         unit matrix with unitary (orthogonal) matrix
NB.         inserted into elements Q[h:h+s-1,h:h+s-1] :
NB.           Q = Π{H(i),i=h:h+s-2}
NB.   hs  - 2-vector of integers (h,s) 'head' and 'size',
NB.         defines submatrix Qf position in matrix HQf, see
NB.         see gehrdu
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (idmat@c (-: clean) (unmhrlln ct@unghrl)@gehrdl) A
NB.   (idmat@c (-: clean) (unmhrllc ct@unghrl)@gehrdl) A
NB.   (idmat@c (-: clean) (unmhrlrn ct@unghrl)@gehrdl) A
NB.   (idmat@c (-: clean) (unmhrlrc ct@unghrl)@gehrdl) A
NB.
NB. Notes:
NB. - implements LAPACK's DORMHR, ZUNMHR
NB. - instead of using f and s parameters, the following
NB.   product is really calculating:
NB.     Q = Π{H(i),i=0:n-1} ,
NB.   where
NB.     H(0:f-1) = H(f+s-1:n-1) = I .
NB.   This approach delivers excessive calculations in rare
NB.   case ((f>0) OR (f+s<n)).

unmhruln=: (unmqrln~ 0 _1&(|.!.0))~
unmhrulc=: (unmqrlc~ 0 _1&(|.!.0))~
unmhrurn=: (unmqrrn~ 0 _1&(|.!.0))~
unmhrurc=: (unmqrrc~ 0 _1&(|.!.0))~

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmrzln   Q * C    left   none  B=. RQf unmrzln C
NB. unmrzlc   Q^H * C  left   ct    B=. RQf unmrzlc C
NB. unmrzrn   C * Q    right  none  B=. RQf unmrzrn C
NB. unmrzrc   C * Q^H  right  ct    B=. RQf unmrzrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form RQf as returned by gerqf
NB.
NB. where
NB.   B,C - m×n-matrices
NB.   LQf - n×(m+1)-matrix (ln,lc cases) or m×(n+1)-matrix
NB.         (rn,rc), contains Qf (unit diagonal not stored),
NB.         the output of gerqf
NB.   Qf  - k×(m+1)-matrix (ln,lc) or k×(n+1)-matrix (rn,rc),
NB.         unit lower triangular, the Q represented in
NB.         factored form
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as a product of k
NB.         elementary reflectors:
NB.           Q = Π{H(i)',i=0:k-1}
NB.   k   = min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (idmat@c (-: clean) (unmrzln ct@(<:@c ungrq ]))@gerqf) A
NB.   (idmat@c (-: clean) (unmrzlc    (<:@c ungrq ]))@gerqf) A
NB.   (idmat@c (-: clean) (unmrzrn ct@(<:@c ungrq ]))@gerqf) A
NB.   (idmat@c (-: clean) (unmrzrc    (<:@c ungrq ]))@gerqf) A
NB.
NB. Notes:
NB. - implements LAPACK's xORMRZ, xUNMRZ
NB. - unmr2{lc,ln,rc,rn} and unmrz{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmrzln_7=: }:  @(((unmr3ln`((larzblcbr&:>/@,~     (- MQNB)&(<\) )~ <)@.(MQNB < #@[))~  tru1                                )~ ,  &0)
unmrzlc_7=: }:  @(((unmr3lc`((larzblnbr&:>/@,~ |.@((- MQNB)&(<\)))~ <)@.(MQNB < #@[))~  tru1                                )~ ,  &0)
unmrzrn_7=: }:"1@(((unmr3rn`((larzbrcbr&:>/@,~ |.@((- MQNB)&(<\)))~ <)@.(MQNB < #@[))~  tru1                                )~ ,. &0)
unmrzrc_7=: }:"1@(((unmr3rc`((larzbrnbr&:>/@,~     (- MQNB)&(<\) )~ <)@.(MQNB < #@[))~  tru1                                )~ ,. &0)

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
NB.
NB. Formula:
NB. - for LQ:
NB.   - for Q * C  : berr := ||Q * C   - Q * C  || / (FP_EPS * ||C|| * s)
NB.   - for Q^H * C: berr := ||Q^H * C - Q^H * C|| / (FP_EPS * ||C|| * s)
NB.   - for C * Q  : berr := ||C * Q   - C * Q  || / (FP_EPS * ||C|| * s)
NB.   - for C * Q^H: berr := ||C * Q^H - C * Q^H|| / (FP_EPS * ||C|| * s)
NB. - for QL:
NB.   - for Q * C  : berr := ||Q * C   - Q * C  || / (FP_EPS * ||C|| * m)
NB.   - for Q^H * C: berr := ||Q^H * C - Q^H * C|| / (FP_EPS * ||C|| * m)
NB.   - for C * Q  : berr := ||C * Q - C * Q    || / (FP_EPS * ||C|| * m)
NB.   - for C * Q^H: berr := ||C * Q^H - C * Q^H|| / (FP_EPS * ||C|| * m)
NB. - for QR:
NB.   - for Q * C  : berr := ||Q * C   - Q * C  || / (FP_EPS * ||C|| * m)
NB.   - for Q^H * C: berr := ||Q^H * C - Q^H * C|| / (FP_EPS * ||C|| * m)
NB.   - for C * Q  : berr := ||C * Q - C * Q    || / (FP_EPS * ||C|| * m)
NB.   - for C * Q^H: berr := ||C * Q^H - C * Q^H|| / (FP_EPS * ||C|| * m)
NB. - for RQ:
NB.   - for Q * C  : berr := ||Q * C   - Q * C  || / (FP_EPS * ||C|| * n)
NB.   - for Q^H * C: berr := ||Q^H * C - Q^H * C|| / (FP_EPS * ||C|| * n)
NB.   - for C * Q  : berr := ||C * Q   - C * Q  || / (FP_EPS * ||C|| * n)
NB.   - for C * Q^H: berr := ||C * Q^H - C * Q^H|| / (FP_EPS * ||C|| * n)

testunmq=: 3 : 0
  'A C'=. y
  rcond=. (_."_)`gecon1@.(=/@$) C  NB. meaninigful for square matrices only
  'LQf QfL QfR RQf'=. xQf=. (gelqf ; geqlf ; geqrf ; gerqf) A
  'Qlq Qql Qqr Qrq'=. (((unglq~ <:@c)&.>)`((ungql~ <:@#)&.>)`((ungqr~ <:@#)&.>)`((ungrq~ <:@c)&.>)) ag xQf

  ('unmlqln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [)))) (LQf;(ct C);    Qlq )
  ('unmlqlc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [)))) (LQf;(ct C);(ct Qlq))
  ('unmlqrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [)))) (LQf;    C ;    Qlq )
  ('unmlqrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [)))) (LQf;    C ;(ct Qlq))

  ('unmqlln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [)))) (QfL;    C ;    Qql )
  ('unmqllc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [)))) (QfL;    C ;(ct Qql))
  ('unmqlrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [)))) (QfL;(ct C);    Qql )
  ('unmqlrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [)))) (QfL;(ct C);(ct Qql))

  ('unmqrln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [)))) (QfR;    C ;    Qqr )
  ('unmqrlc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [)))) (QfR;    C ;(ct Qqr))
  ('unmqrrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [)))) (QfR;(ct C);    Qqr )
  ('unmqrrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [)))) (QfR;(ct C);(ct Qqr))

  ('unmrqln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [)))) (RQf;(ct C);    Qrq )
  ('unmrqlc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [)))) (RQf;(ct C);(ct Qrq))
  ('unmrqrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [)))) (RQf;    C ;    Qrq )
  ('unmrqrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [)))) (RQf;    C ;(ct Qrq))

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
NB.
NB. Formula:
NB. - for lower HRD:
NB.   - for Q * C  : berr := ||Q * C   - Q * C  || / (FP_EPS * ||C|| * n)
NB.   - for Q^H * C: berr := ||Q^H * C - Q^H * C|| / (FP_EPS * ||C|| * n)
NB.   - for C * Q  : berr := ||C * Q   - C * Q  || / (FP_EPS * ||C|| * n)
NB.   - for C * Q^H: berr := ||C * Q^H - C * Q^H|| / (FP_EPS * ||C|| * n)
NB. - for upper HRD:
NB.   - for Q * C  : berr := ||Q * C   - Q * C  || / (FP_EPS * ||C|| * m)
NB.   - for Q^H * C: berr := ||Q^H * C - Q^H * C|| / (FP_EPS * ||C|| * m)
NB.   - for C * Q  : berr := ||C * Q   - C * Q  || / (FP_EPS * ||C|| * m)
NB.   - for C * Q^H: berr := ||C * Q^H - C * Q^H|| / (FP_EPS * ||C|| * m)

testunmhr=: 3 : 0
  'A C'=. y
  rcond=. gecon1 C
  'HlQf HuQf'=. xQf=. ((gehrdl~ 0 , #) ; (gehrdu~ 0 , c)) A
  'Qhrl Qhru'=. ((unghrl&.>)`(unghru&.>)) ag xQf

  ('unmhrlln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`((norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [))))) (HlQf;(ct C);    Qhrl )
  ('unmhrllc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`((norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [))))) (HlQf;(ct C);(ct Qhrl))
  ('unmhrlrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`((norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [))))) (HlQf;    C ;    Qhrl )
  ('unmhrlrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`((norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * c)@(1 {:: [))))) (HlQf;    C ;(ct Qhrl))

  ('unmhruln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`((norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [))))) (HuQf;(ct C);    Qhru )
  ('unmhrulc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`((norm1@(- mp~&>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [))))) (HuQf;(ct C);(ct Qhru))
  ('unmhrurn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`((norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [))))) (HuQf;    C ;    Qhru )
  ('unmhrurc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`((norm1@(- mp &>/@}.)~ % (FP_EPS * norm1 * #)@(1 {:: [))))) (HuQf;    C ;(ct Qhru))

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
NB.     ?@$&0 testmq_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     _1 1 0 4 _6 4&gemat_mt_ testmq_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     (gemat_mt_ j. gemat_mt_) testmq_mt_ 150 200

testmq=: 1 : 'EMPTY_mt_ [ (testunmhr_mt_^:(=/@$@(0&{::)) [ testunmq_mt_)@(u ; u)'
