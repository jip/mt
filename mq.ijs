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
NB. unmhrxxx   Multiply a general matrix by an unitary
NB.            (orthogonal) matrix, which is represented in
NB.            factored form, as returned by gehrdx
NB. unmbrxx    Multiply a general matrix by an unitary
NB.            (orthogonal) matrix, which is represented in
NB.            factored form, as returned by gebrdx
NB.
NB. testunmq   Test unmxxxx by general matrix
NB. testunmhr  Test unmhrxxx by square matrix
NB. testunmz   Test unmxxxx by trapezoidal matrix
NB. testmq     Adv. to make verb to test unmxxxxx by matrix
NB.            of generator and shape given
NB.
NB. Version: 0.10.5 2020-03-30
NB.
NB. Copyright 2010-2020 Igor Zhuravlov
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

mqvberr=: 2 : 'norm1_mt_@(- u&>/@}.)~ % (FP_EPS_mt_ * (1:`]@.*)@norm1_mt_ * v)@(1 {:: [)'  NB. conj. to form verb to calc. berr

NB. ---------------------------------------------------------
NB. Blocked code constants

MQNB=: 32   NB. block size limit

NB. ---------------------------------------------------------
NB. Verb     Action   Side   Tran  Syntax
NB. unml2ln  Q   * C  left   none  eCprod=. Qf unml2ln (C, 0)
NB. unml2lc  Q^H * C  left   ct    eCprod=. Qf unml2lc (C, 0)
NB. unml2rn  C * Q    right  none  eCprod=. Qf unml2rn (C,.0)
NB. unml2rc  C * Q^H  right  ct    eCprod=. Qf unml2rc (C,.0)
NB. unm2lln  Q   * C  left   none  eCprod=. Qf unm2lln (0, C)
NB. unm2llc  Q^H * C  left   ct    eCprod=. Qf unm2llc (0, C)
NB. unm2lrn  C * Q    right  none  eCprod=. Qf unm2lrn (0,.C)
NB. unm2lrc  C * Q^H  right  ct    eCprod=. Qf unm2lrc (0,.C)
NB. unm2rln  Q   * C  left   none  eCprod=. Qf unm2rln (C, 0)
NB. unm2rlc  Q^H * C  left   ct    eCprod=. Qf unm2rlc (C, 0)
NB. unm2rrn  C * Q    right  none  eCprod=. Qf unm2rrn (C,.0)
NB. unm2rrc  C * Q^H  right  ct    eCprod=. Qf unm2rrc (C,.0)
NB. unmr2ln  Q   * C  left   none  eCprod=. Qf unmr2ln (0, C)
NB. unmr2lc  Q^H * C  left   ct    eCprod=. Qf unmr2lc (0, C)
NB. unmr2rn  C * Q    right  none  eCprod=. Qf unmr2rn (0,.C)
NB. unmr2rc  C * Q^H  right  ct    eCprod=. Qf unmr2rc (0,.C)
NB.
NB. Description:
NB.   Multiply a general matrix C, augmented by trash vector,
NB.   by matrix Q. This is non-blocked version of algorithm
NB. where
NB.   C      - m×(n+1)-matrix or (m+1)×n-matrix to multiply
NB.   Qf     - unit triangular matrix, it represents Q in
NB.            factored form as returned by ge{lq,ql,qr,rq}2,
NB.            and contains vectors Vtau[0:k-1]
NB.   Q      - matrix with orthonormal rows or columns, which
NB.            is defined as the product of elementary
NB.            reflectors
NB.   eCprod - being product of matrix Q and augmented matrix
NB.            C, trash vector is modified on exit
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((((}:  @unml2ln ,  &0)~  tru1        @({.  ~  0 _1    <./ @:+ $)) -: (mp~    (unglq~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) LQf
NB.   ((((}:  @unml2lc ,  &0)~  tru1        @({.  ~  0 _1    <./ @:+ $)) -: (mp~ ct@(unglq~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) LQf
NB.   ((((}:"1@unml2rn ,. &0)~  tru1        @({.  ~  0 _1    <./ @:+ $)) -: (mp     (unglq~ <:@c)))~ 0 ?@$~ 0 _1     + $) LQf
NB.   ((((}:"1@unml2rc ,. &0)~  tru1        @({.  ~  0 _1    <./ @:+ $)) -: (mp  ct@(unglq~ <:@c)))~ 0 ?@$~ 0 _1     + $) LQf
NB.   ((((}.  @unm2lln , ~&0)~ (tru1~ -~/@$)@({."1~  _1 0 -@(<./)@:+ $)) -: (mp~    (ungql~ <:@#)))~ 0 ?@$~ _1 0     + $) QfL
NB.   ((((}.  @unm2llc , ~&0)~ (tru1~ -~/@$)@({."1~  _1 0 -@(<./)@:+ $)) -: (mp~ ct@(ungql~ <:@#)))~ 0 ?@$~ _1 0     + $) QfL
NB.   ((((}."1@unm2lrn ,.~&0)~ (tru1~ -~/@$)@({."1~  _1 0 -@(<./)@:+ $)) -: (mp     (ungql~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) QfL
NB.   ((((}."1@unm2lrc ,.~&0)~ (tru1~ -~/@$)@({."1~  _1 0 -@(<./)@:+ $)) -: (mp  ct@(ungql~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) QfL
NB.   ((((}:  @unm2rln ,  &0)~  trl1        @({."1~  _1 0    <./ @:+ $)) -: (mp~    (ungqr~ <:@#)))~ 0 ?@$~ _1 0     + $) QfR
NB.   ((((}:  @unm2rlc ,  &0)~  trl1        @({."1~  _1 0    <./ @:+ $)) -: (mp~ ct@(ungqr~ <:@#)))~ 0 ?@$~ _1 0     + $) QfR
NB.   ((((}:"1@unm2rrn ,. &0)~  trl1        @({."1~  _1 0    <./ @:+ $)) -: (mp     (ungqr~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) QfR
NB.   ((((}:"1@unm2rrc ,. &0)~  trl1        @({."1~  _1 0    <./ @:+ $)) -: (mp  ct@(ungqr~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) QfR
NB.   ((((}.  @unmr2ln , ~&0)~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $)) -: (mp~    (ungrq~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) RQf
NB.   ((((}.  @unmr2lc , ~&0)~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $)) -: (mp~ ct@(ungrq~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) RQf
NB.   ((((}."1@unmr2rn ,.~&0)~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $)) -: (mp     (ungrq~ <:@c)))~ 0 ?@$~ 0 _1     + $) RQf
NB.   ((((}."1@unmr2rc ,.~&0)~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $)) -: (mp  ct@(ungrq~ <:@c)))~ 0 ?@$~ 0 _1     + $) RQf
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - unml2xx implement LAPACK's DORML2, ZUNML2
NB. - unm2lxx implement LAPACK's DORM2L, ZUNM2L
NB. - unm2rxx implement LAPACK's DORM2R, ZUNM2R
NB. - unmr2xx implement LAPACK's DORMR2, ZUNMR2
NB. - unm{l2,2l,2r,r2}{ln,lc,rn,rc} and
NB.   unm{lq,ql,qr,rq}{ln,lc,rn,rc} respectively are
NB.   topologic equivalents

unml2ln=: (larflcfr&:>/@,~ (-@c <\ ,)@|.   )~ <
unml2lc=: (larflnfr&:>/@,~ (-@c <\ ,)      )~ <
unml2rn=: (larfrcfr&:>/@,~ (-@c <\ ,)      )~ <
unml2rc=: (larfrnfr&:>/@,~ (-@c <\ ,)@|.   )~ <

unm2lln=: (larflnbc&:>/@,~ (-@c <\ ,)@|.@|:)~ <
unm2llc=: (larflcbc&:>/@,~ (-@c <\ ,)   @|:)~ <
unm2lrn=: (larfrnbc&:>/@,~ (-@c <\ ,)   @|:)~ <
unm2lrc=: (larfrcbc&:>/@,~ (-@c <\ ,)@|.@|:)~ <

unm2rln=: (larflnfc&:>/@,~ (-@c <\ ,)   @|:)~ <
unm2rlc=: (larflcfc&:>/@,~ (-@c <\ ,)@|.@|:)~ <
unm2rrn=: (larfrnfc&:>/@,~ (-@c <\ ,)@|.@|:)~ <
unm2rrc=: (larfrcfc&:>/@,~ (-@c <\ ,)   @|:)~ <

unmr2ln=: (larflcbr&:>/@,~ (-@c <\ ,)      )~ <
unmr2lc=: (larflnbr&:>/@,~ (-@c <\ ,)@|.   )~ <
unmr2rn=: (larfrcbr&:>/@,~ (-@c <\ ,)@|.   )~ <
unmr2rc=: (larfrnbr&:>/@,~ (-@c <\ ,)      )~ <

NB. ---------------------------------------------------------
NB. Verb     Action   Side   Tran  Syntax
NB. unml3ln  Z   * C  left   none  eCprod=. Zf unml3ln (0, C)
NB. unml3lc  Z^H * C  left   ct    eCprod=. Zf unml3lc (0, C)
NB. unml3rn  C * Z    right  none  eCprod=. Zf unml3rn (0,.C)
NB. unml3rc  C * Z^H  right  ct    eCprod=. Zf unml3rc (0,.C)
NB. unm3lln  Z   * C  left   none  eCprod=. Zf unm3lln (C, 0)
NB. unm3llc  Z^H * C  left   ct    eCprod=. Zf unm3llc (C, 0)
NB. unm3lrn  C * Z    right  none  eCprod=. Zf unm3lrn (C,.0)
NB. unm3lrc  C * Z^H  right  ct    eCprod=. Zf unm3lrc (C,.0)
NB. unm3rln  Z   * C  left   none  eCprod=. Zf unm3rln (0, C)
NB. unm3rlc  Z^H * C  left   ct    eCprod=. Zf unm3rlc (0, C)
NB. unm3rrn  C * Z    right  none  eCprod=. Zf unm3rrn (0,.C)
NB. unm3rrc  C * Z^H  right  ct    eCprod=. Zf unm3rrc (0,.C)
NB. unmr3ln  Z   * C  left   none  eCprod=. Zf unmr3ln (C, 0)
NB. unmr3lc  Z^H * C  left   ct    eCprod=. Zf unmr3lc (C, 0)
NB. unmr3rn  C * Z    right  none  eCprod=. Zf unmr3rn (C,.0)
NB. unmr3rc  C * Z^H  right  ct    eCprod=. Zf unmr3rc (C,.0)
NB.
NB. Description:
NB.   Multiply a general matrix C, augmented by trash vector,
NB.   by matrix Z. This is non-blocked version of algorithm
NB. where
NB.   C      - m×(n+1)-matrix or (m+1)×n-matrix to multiply
NB.   Zf     - unit triangular matrix, it represents Z in
NB.            factored form as returned by tz{lz,zl,zr,rz}3,
NB.            and contains vectors Vtau[0:k-1]
NB.   Z      - matrix with orthonormal rows or columns, which
NB.            is defined as the product of elementary
NB.            reflectors
NB.   eCprod - being product of matrix Z and augmented matrix
NB.            C, trash vector is modified on exit
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((((}.  @unml3ln , ~&0)~ (idmat@[`(a: <@; dhs2liso@(_1 , [))`]}~ #)) -: (mp~    (unglz~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) LZf
NB.   ((((}.  @unml3lc , ~&0)~ (idmat@[`(a: <@; dhs2liso@(_1 , [))`]}~ #)) -: (mp~ ct@(unglz~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) LZf
NB.   ((((}."1@unml3rn ,.~&0)~ (idmat@[`(a: <@; dhs2liso@(_1 , [))`]}~ #)) -: (mp     (unglz~ <:@c)))~ 0 ?@$~ 0 _1     + $) LZf
NB.   ((((}."1@unml3rc ,.~&0)~ (idmat@[`(a: <@; dhs2liso@(_1 , [))`]}~ #)) -: (mp  ct@(unglz~ <:@c)))~ 0 ?@$~ 0 _1     + $) LZf
NB.   ((((}:  @unm3lln ,  &0)~ (idmat@[`(       dhs2liso@( 0 , [))`]}~ c)) -: (mp~    (ungzl~ <:@#)))~ 0 ?@$~ _1 0     + $) ZfL
NB.   ((((}:  @unm3llc ,  &0)~ (idmat@[`(       dhs2liso@( 0 , [))`]}~ c)) -: (mp~ ct@(ungzl~ <:@#)))~ 0 ?@$~ _1 0     + $) ZfL
NB.   ((((}:"1@unm3lrn ,. &0)~ (idmat@[`(       dhs2liso@( 0 , [))`]}~ c)) -: (mp     (ungzl~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) ZfL
NB.   ((((}:"1@unm3lrc ,. &0)~ (idmat@[`(       dhs2liso@( 0 , [))`]}~ c)) -: (mp  ct@(ungzl~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) ZfL
NB.   ((((}.  @unm3rln , ~&0)~ (idmat@[`(       dhs2liso@(_1 , [))`]}~ c)) -: (mp~    (ungzr~ <:@#)))~ 0 ?@$~ _1 0     + $) ZfR
NB.   ((((}.  @unm3rlc , ~&0)~ (idmat@[`(       dhs2liso@(_1 , [))`]}~ c)) -: (mp~ ct@(ungzr~ <:@#)))~ 0 ?@$~ _1 0     + $) ZfR
NB.   ((((}."1@unm3rrn ,.~&0)~ (idmat@[`(       dhs2liso@(_1 , [))`]}~ c)) -: (mp     (ungzr~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) ZfR
NB.   ((((}."1@unm3rrc ,.~&0)~ (idmat@[`(       dhs2liso@(_1 , [))`]}~ c)) -: (mp  ct@(ungzr~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) ZfR
NB.   ((((}:  @unmr3ln ,  &0)~ (idmat@[`(a: <@; dhs2liso@( 0 , [))`]}~ #)) -: (mp~    (ungrz~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) RZf
NB.   ((((}:  @unmr3lc ,  &0)~ (idmat@[`(a: <@; dhs2liso@( 0 , [))`]}~ #)) -: (mp~ ct@(ungrz~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) RZf
NB.   ((((}:"1@unmr3rn ,. &0)~ (idmat@[`(a: <@; dhs2liso@( 0 , [))`]}~ #)) -: (mp     (ungrz~ <:@c)))~ 0 ?@$~ 0 _1     + $) RZf
NB.   ((((}:"1@unmr3rc ,. &0)~ (idmat@[`(a: <@; dhs2liso@( 0 , [))`]}~ #)) -: (mp  ct@(ungrz~ <:@c)))~ 0 ?@$~ 0 _1     + $) RZf
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - unmr3xx implement LAPACK's DORMR3, ZUNMR3
NB. - unm{l3,3l,3r,r3}{ln,lc,rn,rc} and
NB.   unm{lz,zl,zr,rz}{ln,lc,rn,rc} respectively are
NB.   topologic equivalents

unml3ln=: (larzlcfr&:>/@,~ (-@c <\ ,)@|.   )~ <
unml3lc=: (larzlnfr&:>/@,~ (-@c <\ ,)      )~ <
unml3rn=: (larzrcfr&:>/@,~ (-@c <\ ,)      )~ <
unml3rc=: (larzrnfr&:>/@,~ (-@c <\ ,)@|.   )~ <

unm3lln=: (larzlnbc&:>/@,~ (-@c <\ ,)@|.@|:)~ <
unm3llc=: (larzlcbc&:>/@,~ (-@c <\ ,)   @|:)~ <
unm3lrn=: (larzrnbc&:>/@,~ (-@c <\ ,)   @|:)~ <
unm3lrc=: (larzrcbc&:>/@,~ (-@c <\ ,)@|.@|:)~ <

unm3rln=: (larzlnfc&:>/@,~ (-@c <\ ,)   @|:)~ <
unm3rlc=: (larzlcfc&:>/@,~ (-@c <\ ,)@|.@|:)~ <
unm3rrn=: (larzrnfc&:>/@,~ (-@c <\ ,)@|.@|:)~ <
unm3rrc=: (larzrcfc&:>/@,~ (-@c <\ ,)   @|:)~ <

unmr3ln=: (larzlcbr&:>/@,~ (-@c <\ ,)      )~ <
unmr3lc=: (larzlnbr&:>/@,~ (-@c <\ ,)@|.   )~ <
unmr3rn=: (larzrcbr&:>/@,~ (-@c <\ ,)@|.   )~ <
unmr3rc=: (larzrnbr&:>/@,~ (-@c <\ ,)      )~ <

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmlqln   Q   * C  left   none  B=. LQf unmlqln C
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
NB.   ((unmlqln -: (mp~    (unglq~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) LQf
NB.   ((unmlqlc -: (mp~ ct@(unglq~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) LQf
NB.   ((unmlqrn -: (mp     (unglq~ <:@c))~) 0 ?@$~ 0 _1     + $) LQf
NB.   ((unmlqrc -: (mp  ct@(unglq~ <:@c))~) 0 ?@$~ 0 _1     + $) LQf
NB.
NB. Notes:
NB. - implement LAPACK's DORMLQ, ZUNMLQ
NB. - unml2{lc,ln,rc,rn} and unmlq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmlqln=: }:  @(((unml2ln`((larfblcfr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,  &0)
unmlqlc=: }:  @(((unml2lc`((larfblnfr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,  &0)
unmlqrn=: }:"1@(((unml2rn`((larfbrcfr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,. &0)
unmlqrc=: }:"1@(((unml2rc`((larfbrnfr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,. &0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmqlln   Q   * C  left   none  B=. QfL unmqlln C
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
NB.   ((unmqlln -: (mp~    (ungql~ <:@#))~) 0 ?@$~ _1 0     + $) QfL
NB.   ((unmqllc -: (mp~ ct@(ungql~ <:@#))~) 0 ?@$~ _1 0     + $) QfL
NB.   ((unmqlrn -: (mp     (ungql~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) QfL
NB.   ((unmqlrc -: (mp  ct@(ungql~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) QfL
NB.
NB. Notes:
NB. - implement LAPACK's DORMQL, ZUNMQL
NB. - unm2l{lc,ln,rc,rn} and unmql{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmqlln=: }.  @(((unm2lln`((larfblnbc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ , ~&0)
unmqllc=: }.  @(((unm2llc`((larfblcbc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ , ~&0)
unmqlrn=: }."1@(((unm2lrn`((larfbrnbc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ ,.~&0)
unmqlrc=: }."1@(((unm2lrc`((larfbrcbc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ ,.~&0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmqrln   Q   * C  left   none  B=. QfR unmqrln C
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
NB.   ((unmqrln -: (mp~    (ungqr~ <:@#))~) 0 ?@$~ _1 0     + $) QfR
NB.   ((unmqrlc -: (mp~ ct@(ungqr~ <:@#))~) 0 ?@$~ _1 0     + $) QfR
NB.   ((unmqrrn -: (mp     (ungqr~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) QfR
NB.   ((unmqrrc -: (mp  ct@(ungqr~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) QfR
NB.
NB. Notes:
NB. - implement LAPACK's DORMQR, ZUNMQR
NB. - unm2r{lc,ln,rc,rn} and unmqr{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmqrln=: }:  @(((unm2rln`((larfblnfc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,  &0)
unmqrlc=: }:  @(((unm2rlc`((larfblcfc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,  &0)
unmqrrn=: }:"1@(((unm2rrn`((larfbrnfc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,. &0)
unmqrrc=: }:"1@(((unm2rrc`((larfbrcfc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,. &0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmrqln   Q   * C  left   none  B=. RQf unmrqln C
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
NB.   RQf - l×(m+1)-matrix (ln,lc cases) or l×(n+1)-matrix
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
NB.   ((unmrqln -: (mp~    (ungrq~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) RQf
NB.   ((unmrqlc -: (mp~ ct@(ungrq~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) RQf
NB.   ((unmrqrn -: (mp     (ungrq~ <:@c))~) 0 ?@$~ 0 _1     + $) RQf
NB.   ((unmrqrc -: (mp  ct@(ungrq~ <:@c))~) 0 ?@$~ 0 _1     + $) RQf
NB.
NB. Notes:
NB. - implement LAPACK's DORMRQ, ZUNMRQ
NB. - unmr2{lc,ln,rc,rn} and unmrq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmrqln=: }.  @(((unmr2ln`((larfblcbr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ , ~&0)
unmrqlc=: }.  @(((unmr2lc`((larfblnbr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ , ~&0)
unmrqrn=: }."1@(((unmr2rn`((larfbrcbr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ ,.~&0)
unmrqrc=: }."1@(((unmr2rc`((larfbrnbr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ ,.~&0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmlzln   Z   * C  left   none  B=. LZf unmlzln C
NB. unmlzlc   Z^H * C  left   ct    B=. LZf unmlzlc C
NB. unmlzrn   C * Z    right  none  B=. LZf unmlzrn C
NB. unmlzrc   C * Z^H  right  ct    B=. LZf unmlzrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Z, which is
NB.   represented in factored form LZf as returned by tzlzf
NB.
NB. where
NB.   B,C - m×n-matrices
NB.   LZf - k×(m+1)-matrix (ln,lc cases) or k×(n+1)-matrix
NB.         (rn,rc), contains Zf (identity submatrix not
NB.         stored), the output of tzlzf
NB.   Zf  - k×(m+1)-matrix (ln,lc) or k×(n+1)-matrix (rn,rc),
NB.         trailing k×k-submatrix is identity, the Z
NB.         represented in factored form
NB.   Z   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as a product of k
NB.         elementary reflectors:
NB.           Z = Π{H(i)',i=k-1:0}
NB.   k   ≤ min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((unmlzln -: (mp~    (unglz~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) LZf
NB.   ((unmlzlc -: (mp~ ct@(unglz~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) LZf
NB.   ((unmlzrn -: (mp     (unglz~ <:@c))~) 0 ?@$~ 0 _1     + $) LZf
NB.   ((unmlzrc -: (mp  ct@(unglz~ <:@c))~) 0 ?@$~ 0 _1     + $) LZf
NB.
NB. Notes:
NB. - unml3{lc,ln,rc,rn} and unmlz{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmlzln=: }.  @(((unml3ln`((larzblcfr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; dhs2liso@(_1 , [))`]}~ #))~ , ~&0)
unmlzlc=: }.  @(((unml3lc`((larzblnfr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; dhs2liso@(_1 , [))`]}~ #))~ , ~&0)
unmlzrn=: }."1@(((unml3rn`((larzbrcfr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; dhs2liso@(_1 , [))`]}~ #))~ ,.~&0)
unmlzrc=: }."1@(((unml3rc`((larzbrnfr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; dhs2liso@(_1 , [))`]}~ #))~ ,.~&0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmzlln   Z   * C  left   none  B=. ZfL unmzlln C
NB. unmzllc   Z^H * C  left   ct    B=. ZfL unmzllc C
NB. unmzlrn   C * Z    right  none  B=. ZfL unmzlrn C
NB. unmzlrc   C * Z^H  right  ct    B=. ZfL unmzlrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Z, which is
NB.   represented in factored form ZfL as returned by tzzlf
NB.
NB. where
NB.   B,C - m×n-matrices
NB.   ZfL - (m+1)×k-matrix (ln,lc cases) or (n+1)×k-matrix
NB.         (rn,rc), contains Zf (identity submatrix not
NB.         stored), the output of tzzlf
NB.   Zf  - (m+1)×k-matrix (ln,lc) or (n+1)×k-matrix (rn,rc),
NB.         leading k×k-submatrix is identity, the Z
NB.         represented in factored form
NB.   Z   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as a product of k
NB.         elementary reflectors:
NB.           Z = Π{H(i),i=k-1:0}
NB.   k   ≤ min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((unmzlln -: (mp~    (ungzl~ <:@#))~) 0 ?@$~ _1 0     + $) ZfL
NB.   ((unmzllc -: (mp~ ct@(ungzl~ <:@#))~) 0 ?@$~ _1 0     + $) ZfL
NB.   ((unmzlrn -: (mp     (ungzl~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) ZfL
NB.   ((unmzlrc -: (mp  ct@(ungzl~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) ZfL
NB.
NB. Notes:
NB. - unm3l{lc,ln,rc,rn} and unmzl{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmzlln=: }:  @(((unm3lln`((larzblnbc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       dhs2liso@( 0 , [))`]}~ c))~ ,  &0)
unmzllc=: }:  @(((unm3llc`((larzblcbc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       dhs2liso@( 0 , [))`]}~ c))~ ,  &0)
unmzlrn=: }:"1@(((unm3lrn`((larzbrnbc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       dhs2liso@( 0 , [))`]}~ c))~ ,. &0)
unmzlrc=: }:"1@(((unm3lrc`((larzbrcbc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       dhs2liso@( 0 , [))`]}~ c))~ ,. &0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmzrln   Z   * C  left   none  B=. ZfR unmzrln C
NB. unmzrlc   Z^H * C  left   ct    B=. ZfR unmzrlc C
NB. unmzrrn   C * Z    right  none  B=. ZfR unmzrrn C
NB. unmzrrc   C * Z^H  right  ct    B=. ZfR unmzrrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Z, which is
NB.   represented in factored form ZfR as returned by tzzrf
NB.
NB. where
NB.   B,C - m×n-matrices
NB.   ZfR - (m+1)×k-matrix (ln,lc cases) or (n+1)×k-matrix
NB.         (rn,rc), contains Zf (identity submatrix not
NB.         stored), the output of tzzrf
NB.   Zf  - (m+1)×k-matrix (ln,lc) or (n+1)×k-matrix (rn,rc),
NB.         trailing k×k-submatrix is identity, the Z
NB.         represented in factored form
NB.   Z   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as a product of k
NB.         elementary reflectors:
NB.           Z = Π{H(i),i=0:k-1}
NB.   k   ≤ min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((unmzrln -: (mp~    (ungzr~ <:@#))~) 0 ?@$~ _1 0     + $) ZfR
NB.   ((unmzrlc -: (mp~ ct@(ungzr~ <:@#))~) 0 ?@$~ _1 0     + $) ZfR
NB.   ((unmzrrn -: (mp     (ungzr~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) ZfR
NB.   ((unmzrrc -: (mp  ct@(ungzr~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) ZfR
NB.
NB. Notes:
NB. - unm3r{lc,ln,rc,rn} and unmzr{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmzrln=: }.  @(((unm3rln`((larzblnfc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       dhs2liso@(_1 , [))`]}~ c))~ , ~&0)
unmzrlc=: }.  @(((unm3rlc`((larzblcfc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       dhs2liso@(_1 , [))`]}~ c))~ , ~&0)
unmzrrn=: }."1@(((unm3rrn`((larzbrnfc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       dhs2liso@(_1 , [))`]}~ c))~ ,.~&0)
unmzrrc=: }."1@(((unm3rrc`((larzbrcfc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       dhs2liso@(_1 , [))`]}~ c))~ ,.~&0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmrzln   Z   * C  left   none  B=. RZf unmrzln C
NB. unmrzlc   Z^H * C  left   ct    B=. RZf unmrzlc C
NB. unmrzrn   C * Z    right  none  B=. RZf unmrzrn C
NB. unmrzrc   C * Z^H  right  ct    B=. RZf unmrzrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Z, which is
NB.   represented in factored form RZf as returned by tzrzf
NB.
NB. where
NB.   B,C - m×n-matrices
NB.   RZf - k×(m+1)-matrix (ln,lc cases) or k×(n+1)-matrix
NB.         (rn,rc), contains Zf (identity submatrix not
NB.         stored), the output of tzrzf
NB.   Zf  - k×(m+1)-matrix (ln,lc) or k×(n+1)-matrix (rn,rc),
NB.         leading k×k-submatrix is identity, the Z
NB.         represented in factored form
NB.   Z   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as a product of k
NB.         elementary reflectors:
NB.           Z = Π{H(i)',i=0:k-1}
NB.   k   ≤ min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((unmrzln -: (mp~    (ungrz~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) RZf
NB.   ((unmrzlc -: (mp~ ct@(ungrz~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) RZf
NB.   ((unmrzrn -: (mp     (ungrz~ <:@c))~) 0 ?@$~ 0 _1     + $) RZf
NB.   ((unmrzrc -: (mp  ct@(ungrz~ <:@c))~) 0 ?@$~ 0 _1     + $) RZf
NB.
NB. Notes:
NB. - implement LAPACK's xORMRZ, xUNMRZ
NB. - unmr3{lc,ln,rc,rn} and unmrz{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unmrzln=: }:  @(((unmr3ln`((larzblcbr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; dhs2liso@( 0 , [))`]}~ #))~ ,  &0)
unmrzlc=: }:  @(((unmr3lc`((larzblnbr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; dhs2liso@( 0 , [))`]}~ #))~ ,  &0)
unmrzrn=: }:"1@(((unmr3rn`((larzbrcbr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; dhs2liso@( 0 , [))`]}~ #))~ ,. &0)
unmrzrc=: }:"1@(((unmr3rc`((larzbrnbr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; dhs2liso@( 0 , [))`]}~ #))~ ,. &0)

NB. ---------------------------------------------------------
NB. Verb      Action   Side   Tran  Syntax
NB. unmhrlln  Q   * C  left   none  B=. HQf unmhrlln C
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
NB.   (idmat@c (-: clean) (unmhrlln ct@unghrl)@(gehrdl~ 0 , c)) A
NB.   (idmat@c (-: clean) (unmhrllc    unghrl)@(gehrdl~ 0 , c)) A
NB.   (idmat@c (-: clean) (unmhrlrn ct@unghrl)@(gehrdl~ 0 , c)) A
NB.   (idmat@c (-: clean) (unmhrlrc    unghrl)@(gehrdl~ 0 , c)) A
NB.
NB. Notes:
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
NB. unmhruln  Q   * C  left   none  B=. HQf unmhruln C
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
NB.   (idmat@# (-: clean) (unmhruln ct@unghru)@(gehrdu~ 0 , #)) A
NB.   (idmat@# (-: clean) (unmhrulc    unghru)@(gehrdu~ 0 , #)) A
NB.   (idmat@# (-: clean) (unmhrurn ct@unghru)@(gehrdu~ 0 , #)) A
NB.   (idmat@# (-: clean) (unmhrurc    unghru)@(gehrdu~ 0 , #)) A
NB.
NB. Notes:
NB. - models LAPACK's DORMHR, ZUNMHR
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

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testunmq
NB.
NB. Description:
NB.   Test Q multiplication qf-algorithms by general matrix
NB.
NB. Syntax:
NB.   testunmq (A;C)
NB. where
NB.   A - m×n-matrix, is used to produce Qf
NB.   C - m×n-matrix, is used as multiplier
NB.
NB. Formula:
NB. - for LQ:
NB.   - for Q   * C: berr := ||(Q  ) * C - Q   * C|| / (FP_EPS * ||C|| * n)
NB.   - for Q^H * C: berr := ||(Q^H) * C - Q^H * C|| / (FP_EPS * ||C|| * n)
NB.   - for C * Q  : berr := ||C * (Q  ) - C * Q  || / (FP_EPS * ||C|| * n)
NB.   - for C * Q^H: berr := ||C * (Q^H) - C * Q^H|| / (FP_EPS * ||C|| * n)
NB. - for QL:
NB.   - for Q   * C: berr := ||(Q  ) * C - Q   * C|| / (FP_EPS * ||C|| * m)
NB.   - for Q^H * C: berr := ||(Q^H) * C - Q^H * C|| / (FP_EPS * ||C|| * m)
NB.   - for C * Q  : berr := ||C * (Q  ) - C * Q  || / (FP_EPS * ||C|| * m)
NB.   - for C * Q^H: berr := ||C * (Q^H) - C * Q^H|| / (FP_EPS * ||C|| * m)
NB. - for QR:
NB.   - for Q   * C: berr := ||(Q  ) * C - Q   * C|| / (FP_EPS * ||C|| * m)
NB.   - for Q^H * C: berr := ||(Q^H) * C - Q^H * C|| / (FP_EPS * ||C|| * m)
NB.   - for C * Q  : berr := ||C * (Q  ) - C * Q  || / (FP_EPS * ||C|| * m)
NB.   - for C * Q^H: berr := ||C * (Q^H) - C * Q^H|| / (FP_EPS * ||C|| * m)
NB. - for RQ:
NB.   - for Q   * C: berr := ||(Q  ) * C - Q   * C|| / (FP_EPS * ||C|| * n)
NB.   - for Q^H * C: berr := ||(Q^H) * C - Q^H * C|| / (FP_EPS * ||C|| * n)
NB.   - for C * Q  : berr := ||C * (Q  ) - C * Q  || / (FP_EPS * ||C|| * n)
NB.   - for C * Q^H: berr := ||C * (Q^H) - C * Q^H|| / (FP_EPS * ||C|| * n)

testunmq=: 3 : 0
  'A C'=. y
  rcond=. (_."_)`gecon1@.(=/@$) C  NB. meaninigful for square matrices only

  Qlq=. (unglq~ <:@c) LQf=. gelqf A
  Qql=. (ungql~ <:@#) QfL=. geqlf A
  Qqr=. (ungqr~ <:@#) QfR=. geqrf A
  Qrq=. (ungrq~ <:@c) RQf=. gerqf A

  ('unmlqln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr c))) LQf ; (ct C) ;    Qlq
  ('unmlqlc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr c))) LQf ; (ct C) ; ct Qlq
  ('unmlqrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr c))) LQf ;     C  ;    Qlq
  ('unmlqrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr c))) LQf ;     C  ; ct Qlq

  ('unmqlln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr #))) QfL ;     C  ;    Qql
  ('unmqllc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr #))) QfL ;     C  ; ct Qql
  ('unmqlrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr #))) QfL ; (ct C) ;    Qql
  ('unmqlrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr #))) QfL ; (ct C) ; ct Qql

  ('unmqrln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr #))) QfR ;     C  ;    Qqr
  ('unmqrlc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr #))) QfR ;     C  ; ct Qqr
  ('unmqrrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr #))) QfR ; (ct C) ;    Qqr
  ('unmqrrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr #))) QfR ; (ct C) ; ct Qqr

  ('unmrqln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr c))) RQf ; (ct C) ;    Qrq
  ('unmrqlc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr c))) RQf ; (ct C) ; ct Qrq
  ('unmrqrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr c))) RQf ;     C  ;    Qrq
  ('unmrqrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr c))) RQf ;     C  ; ct Qrq

  EMPTY
)

NB. ---------------------------------------------------------
NB. testunmz
NB.
NB. Description:
NB.   Test Z multiplication zf-algorithms by trapezoidal
NB.   matrix
NB.
NB. Syntax:
NB.   testunmz (A;C)
NB. where
NB.   A - m×n-matrix, is used to produce Zf
NB.   C - m×n-matrix, is used as multiplier
NB.
NB. Formula:
NB. - for LZ:
NB.   - for Z   * C: berr := ||(Z  ) * C - Z   * C|| / (FP_EPS * ||C|| * n)
NB.   - for Z^H * C: berr := ||(Z^H) * C - Z^H * C|| / (FP_EPS * ||C|| * n)
NB.   - for C * Z  : berr := ||C * (Z  ) - C * Z  || / (FP_EPS * ||C|| * n)
NB.   - for C * Z^H: berr := ||C * (Z^H) - C * Z^H|| / (FP_EPS * ||C|| * n)
NB. - for ZL:
NB.   - for Z   * C: berr := ||(Z  ) * C - Z   * C|| / (FP_EPS * ||C|| * m)
NB.   - for Z^H * C: berr := ||(Z^H) * C - Z^H * C|| / (FP_EPS * ||C|| * m)
NB.   - for C * Z  : berr := ||C * (Z  ) - C * Z  || / (FP_EPS * ||C|| * m)
NB.   - for C * Z^H: berr := ||C * (Z^H) - C * Z^H|| / (FP_EPS * ||C|| * m)
NB. - for ZR:
NB.   - for Z   * C: berr := ||(Z  ) * C - Z   * C|| / (FP_EPS * ||C|| * m)
NB.   - for Z^H * C: berr := ||(Z^H) * C - Z^H * C|| / (FP_EPS * ||C|| * m)
NB.   - for C * Z  : berr := ||C * (Z  ) - C * Z  || / (FP_EPS * ||C|| * m)
NB.   - for C * Z^H: berr := ||C * (Z^H) - C * Z^H|| / (FP_EPS * ||C|| * m)
NB. - for RZ:
NB.   - for Z   * C: berr := ||(Z  ) * C - Z   * C|| / (FP_EPS * ||C|| * n)
NB.   - for Z^H * C: berr := ||(Z^H) * C - Z^H * C|| / (FP_EPS * ||C|| * n)
NB.   - for C * Z  : berr := ||C * (Z  ) - C * Z  || / (FP_EPS * ||C|| * n)
NB.   - for C * Z^H: berr := ||C * (Z^H) - C * Z^H|| / (FP_EPS * ||C|| * n)

testunmz=: 3 : 0
  'A C'=. y
  rcond=. (_."_)`gecon1@.(=/@$) C  NB. meaninigful for square matrices only

  Awide=. |:^:(>/@$) A
  Atall=. |:^:(</@$) A

  Cwide=. |:^:(>/@$) C
  Ctall=. |:^:(</@$) C

  Zlz=. (unglz~ <:@c) LZf=. tzlzf (trl~ -~/@$) Awide
  Zzl=. (ungzl~ <:@#) ZfL=. tzzlf  trl         Atall
  Zzr=. (ungzr~ <:@#) ZfR=. tzzrf (tru~ -~/@$) Atall
  Zrz=. (ungrz~ <:@c) RZf=. tzrzf  tru         Awide

  ('unmlzln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr c))) LZf ; Ctall ;    Zlz
  ('unmlzlc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr c))) LZf ; Ctall ; ct Zlz
  ('unmlzrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr c))) LZf ; Cwide ;    Zlz
  ('unmlzrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr c))) LZf ; Cwide ; ct Zlz

  ('unmzlln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr #))) ZfL ; Ctall ;    Zzl
  ('unmzllc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr #))) ZfL ; Ctall ; ct Zzl
  ('unmzlrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr #))) ZfL ; Cwide ;    Zzl
  ('unmzlrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr #))) ZfL ; Cwide ; ct Zzl

  ('unmzrln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr #))) ZfR ; Ctall ;    Zzr
  ('unmzrlc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr #))) ZfR ; Ctall ; ct Zzr
  ('unmzrrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr #))) ZfR ; Cwide ;    Zzr
  ('unmzrrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr #))) ZfR ; Cwide ; ct Zzr

  ('unmrzln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr c))) RZf ; Ctall ;    Zrz
  ('unmrzlc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr c))) RZf ; Ctall ; ct Zrz
  ('unmrzrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr c))) RZf ; Cwide ;    Zrz
  ('unmrzrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr c))) RZf ; Cwide ; ct Zrz

  EMPTY
)

NB. ---------------------------------------------------------
NB. testunmhr
NB.
NB. Description:
NB.   Test Q multiplication hrd-algorithms by square matrix
NB.
NB. Syntax:
NB.   testunmhr (A;C)
NB. where
NB.   A - n×n-matrix, is used to produce Qf
NB.   C - n×n-matrix, is used as multiplier
NB.
NB. Formula:
NB. - for lower HRD:
NB.   - for Q   * C: berr := ||(Q  ) * C - Q   * C|| / (FP_EPS * ||C|| * n)
NB.   - for Q^H * C: berr := ||(Q^H) * C - Q^H * C|| / (FP_EPS * ||C|| * n)
NB.   - for C * Q  : berr := ||C * (Q  ) - C * Q  || / (FP_EPS * ||C|| * n)
NB.   - for C * Q^H: berr := ||C * (Q^H) - C * Q^H|| / (FP_EPS * ||C|| * n)
NB. - for upper HRD:
NB.   - for Q   * C: berr := ||(Q  ) * C - Q   * C|| / (FP_EPS * ||C|| * m)
NB.   - for Q^H * C: berr := ||(Q^H) * C - Q^H * C|| / (FP_EPS * ||C|| * m)
NB.   - for C * Q  : berr := ||C * (Q  ) - C * Q  || / (FP_EPS * ||C|| * m)
NB.   - for C * Q^H: berr := ||C * (Q^H) - C * Q^H|| / (FP_EPS * ||C|| * m)

testunmhr=: 3 : 0
  'A C'=. y
  rcond=. gecon1 C

  Qhrl=. unghrl HlQf=. (gehrdl~ 0 , c) A
  Qhru=. unghru HuQf=. (gehrdu~ 0 , #) A

  ('unmhrlln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr c))) HlQf ; C ;    Qhrl
  ('unmhrllc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr c))) HlQf ; C ; ct Qhrl
  ('unmhrlrn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr c))) HlQf ; C ;    Qhrl
  ('unmhrlrc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr c))) HlQf ; C ; ct Qhrl

  ('unmhruln' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr #))) HuQf ; C ;    Qhru
  ('unmhrulc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp~ mqvberr #))) HuQf ; C ; ct Qhru
  ('unmhrurn' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr #))) HuQf ; C ;    Qhru
  ('unmhrurc' tdyad ((0&{::)`(1&{::)`]`(rcond"_)`(_."_)`(mp  mqvberr #))) HuQf ; C ; ct Qhru

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

testmq=: 1 : 'EMPTY [ (testunmhr_mt_^:(=/@$@(0&{::)) [ testunmz_mt_ [ testunmq_mt_)@(u ; u)'
