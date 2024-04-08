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
NB.
NB. testunmq   Test unmxxxx by general matrix
NB. testunmz   Test unmxxxx by trapezoidal matrix
NB. testunmhr  Test unmhrxxx by square matrix
NB. testmq     Adv. to make verb to test unmxxxxx by matrix
NB.            of generator and shape given
NB.
NB. Copyright 2010,2011,2013,2017,2018,2020,2021,2023,2024
NB.           Igor Zhuravlov
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

NB. =========================================================
NB. Configuration

coclass 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. Blocked code constants

MQNB=: 32  NB. block size limit

NB. ---------------------------------------------------------
NB. Verb       Action     Side     Tran    Syntax
NB. unml2ln    Q   * C    left     none    eCprod=. Qf unml2ln (C, 0)
NB. unml2lc    Q^H * C    left     ct      eCprod=. Qf unml2lc (C, 0)
NB. unml2rn    C * Q      right    none    eCprod=. Qf unml2rn (C,.0)
NB. unml2rc    C * Q^H    right    ct      eCprod=. Qf unml2rc (C,.0)
NB. unm2lln    Q   * C    left     none    eCprod=. Qf unm2lln (0, C)
NB. unm2llc    Q^H * C    left     ct      eCprod=. Qf unm2llc (0, C)
NB. unm2lrn    C * Q      right    none    eCprod=. Qf unm2lrn (0,.C)
NB. unm2lrc    C * Q^H    right    ct      eCprod=. Qf unm2lrc (0,.C)
NB. unm2rln    Q   * C    left     none    eCprod=. Qf unm2rln (C, 0)
NB. unm2rlc    Q^H * C    left     ct      eCprod=. Qf unm2rlc (C, 0)
NB. unm2rrn    C * Q      right    none    eCprod=. Qf unm2rrn (C,.0)
NB. unm2rrc    C * Q^H    right    ct      eCprod=. Qf unm2rrc (C,.0)
NB. unmr2ln    Q   * C    left     none    eCprod=. Qf unmr2ln (0, C)
NB. unmr2lc    Q^H * C    left     ct      eCprod=. Qf unmr2lc (0, C)
NB. unmr2rn    C * Q      right    none    eCprod=. Qf unmr2rn (0,.C)
NB. unmr2rc    C * Q^H    right    ct      eCprod=. Qf unmr2rc (0,.C)
NB.
NB. Description:
NB.   Multiply a general matrix C, augmented by trash vector,
NB.   by matrix Q. This is non-blocked version of algorithm
NB. where
NB.   C      - m×(n+1)-matrix or (m+1)×n-matrix to multiply
NB.   Qf     - unit trapezoidal matrix, it represents Q in
NB.            factored form as returned by ge{lq,ql,qr,rq}2,
NB.            and contains vectors Vtau[0:k-1]
NB.   Q      - unitary (orthogonal) matrix, which is defined
NB.            as the product of elementary reflectors
NB.   eCprod - the product of matrix Q and the augmented
NB.            matrix C, trash vector is modified on exit
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

unml2ln=: ((larflcfr&:>/@,~ (-@c <\ ,)@|.   )~ <)^:(0 < #@[)
unml2lc=: ((larflnfr&:>/@,~ (-@c <\ ,)      )~ <)^:(0 < #@[)
unml2rn=: ((larfrcfr&:>/@,~ (-@c <\ ,)      )~ <)^:(0 < #@[)
unml2rc=: ((larfrnfr&:>/@,~ (-@c <\ ,)@|.   )~ <)^:(0 < #@[)

unm2lln=: ((larflnbc&:>/@,~ (-@c <\ ,)@|.@|:)~ <)^:(0 < c@[)
unm2llc=: ((larflcbc&:>/@,~ (-@c <\ ,)   @|:)~ <)^:(0 < c@[)
unm2lrn=: ((larfrnbc&:>/@,~ (-@c <\ ,)   @|:)~ <)^:(0 < c@[)
unm2lrc=: ((larfrcbc&:>/@,~ (-@c <\ ,)@|.@|:)~ <)^:(0 < c@[)

unm2rln=: ((larflnfc&:>/@,~ (-@c <\ ,)   @|:)~ <)^:(0 < c@[)
unm2rlc=: ((larflcfc&:>/@,~ (-@c <\ ,)@|.@|:)~ <)^:(0 < c@[)
unm2rrn=: ((larfrnfc&:>/@,~ (-@c <\ ,)@|.@|:)~ <)^:(0 < c@[)
unm2rrc=: ((larfrcfc&:>/@,~ (-@c <\ ,)   @|:)~ <)^:(0 < c@[)

unmr2ln=: ((larflcbr&:>/@,~ (-@c <\ ,)      )~ <)^:(0 < #@[)
unmr2lc=: ((larflnbr&:>/@,~ (-@c <\ ,)@|.   )~ <)^:(0 < #@[)
unmr2rn=: ((larfrcbr&:>/@,~ (-@c <\ ,)@|.   )~ <)^:(0 < #@[)
unmr2rc=: ((larfrnbr&:>/@,~ (-@c <\ ,)      )~ <)^:(0 < #@[)

NB. ---------------------------------------------------------
NB. Verb       Action     Side     Tran    Syntax
NB. unml3ln    Z   * C    left     none    eCprod=. Zf unml3ln (0, C)
NB. unml3lc    Z^H * C    left     ct      eCprod=. Zf unml3lc (0, C)
NB. unml3rn    C * Z      right    none    eCprod=. Zf unml3rn (0,.C)
NB. unml3rc    C * Z^H    right    ct      eCprod=. Zf unml3rc (0,.C)
NB. unm3lln    Z   * C    left     none    eCprod=. Zf unm3lln (C, 0)
NB. unm3llc    Z^H * C    left     ct      eCprod=. Zf unm3llc (C, 0)
NB. unm3lrn    C * Z      right    none    eCprod=. Zf unm3lrn (C,.0)
NB. unm3lrc    C * Z^H    right    ct      eCprod=. Zf unm3lrc (C,.0)
NB. unm3rln    Z   * C    left     none    eCprod=. Zf unm3rln (0, C)
NB. unm3rlc    Z^H * C    left     ct      eCprod=. Zf unm3rlc (0, C)
NB. unm3rrn    C * Z      right    none    eCprod=. Zf unm3rrn (0,.C)
NB. unm3rrc    C * Z^H    right    ct      eCprod=. Zf unm3rrc (0,.C)
NB. unmr3ln    Z   * C    left     none    eCprod=. Zf unmr3ln (C, 0)
NB. unmr3lc    Z^H * C    left     ct      eCprod=. Zf unmr3lc (C, 0)
NB. unmr3rn    C * Z      right    none    eCprod=. Zf unmr3rn (C,.0)
NB. unmr3rc    C * Z^H    right    ct      eCprod=. Zf unmr3rc (C,.0)
NB.
NB. Description:
NB.   Multiply a general matrix C, augmented by trash vector,
NB.   by matrix Z. This is non-blocked version of algorithm
NB. where
NB.   C      - m×(n+1)-matrix or (m+1)×n-matrix to multiply
NB.   Zf     - unit trapezoidal matrix, it represents Z in
NB.            factored form as returned by tz{lz,zl,zr,rz}3,
NB.            and contains vectors Vtau[0:k-1]
NB.   Z      - unitary (orthonormal) matrix which is defined
NB.            as the product of elementary reflectors
NB.   eCprod - the product of matrix Z and the augmented
NB.            matrix C, trash vector is modified on exit
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((((}.  @unml3ln , ~&0)~ (idmat@[`(a: <@; liso4dhs@(_1 , [))`]}~ #)) -: (mp~    (unglz~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) LZf
NB.   ((((}.  @unml3lc , ~&0)~ (idmat@[`(a: <@; liso4dhs@(_1 , [))`]}~ #)) -: (mp~ ct@(unglz~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) LZf
NB.   ((((}."1@unml3rn ,.~&0)~ (idmat@[`(a: <@; liso4dhs@(_1 , [))`]}~ #)) -: (mp     (unglz~ <:@c)))~ 0 ?@$~ 0 _1     + $) LZf
NB.   ((((}."1@unml3rc ,.~&0)~ (idmat@[`(a: <@; liso4dhs@(_1 , [))`]}~ #)) -: (mp  ct@(unglz~ <:@c)))~ 0 ?@$~ 0 _1     + $) LZf
NB.   ((((}:  @unm3lln ,  &0)~ (idmat@[`(       liso4dhs@( 0 , [))`]}~ c)) -: (mp~    (ungzl~ <:@#)))~ 0 ?@$~ _1 0     + $) ZfL
NB.   ((((}:  @unm3llc ,  &0)~ (idmat@[`(       liso4dhs@( 0 , [))`]}~ c)) -: (mp~ ct@(ungzl~ <:@#)))~ 0 ?@$~ _1 0     + $) ZfL
NB.   ((((}:"1@unm3lrn ,. &0)~ (idmat@[`(       liso4dhs@( 0 , [))`]}~ c)) -: (mp     (ungzl~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) ZfL
NB.   ((((}:"1@unm3lrc ,. &0)~ (idmat@[`(       liso4dhs@( 0 , [))`]}~ c)) -: (mp  ct@(ungzl~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) ZfL
NB.   ((((}.  @unm3rln , ~&0)~ (idmat@[`(       liso4dhs@(_1 , [))`]}~ c)) -: (mp~    (ungzr~ <:@#)))~ 0 ?@$~ _1 0     + $) ZfR
NB.   ((((}.  @unm3rlc , ~&0)~ (idmat@[`(       liso4dhs@(_1 , [))`]}~ c)) -: (mp~ ct@(ungzr~ <:@#)))~ 0 ?@$~ _1 0     + $) ZfR
NB.   ((((}."1@unm3rrn ,.~&0)~ (idmat@[`(       liso4dhs@(_1 , [))`]}~ c)) -: (mp     (ungzr~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) ZfR
NB.   ((((}."1@unm3rrc ,.~&0)~ (idmat@[`(       liso4dhs@(_1 , [))`]}~ c)) -: (mp  ct@(ungzr~ <:@#)))~ 0 ?@$~ _1 0 |.@:+ $) ZfR
NB.   ((((}:  @unmr3ln ,  &0)~ (idmat@[`(a: <@; liso4dhs@( 0 , [))`]}~ #)) -: (mp~    (ungrz~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) RZf
NB.   ((((}:  @unmr3lc ,  &0)~ (idmat@[`(a: <@; liso4dhs@( 0 , [))`]}~ #)) -: (mp~ ct@(ungrz~ <:@c)))~ 0 ?@$~ 0 _1 |.@:+ $) RZf
NB.   ((((}:"1@unmr3rn ,. &0)~ (idmat@[`(a: <@; liso4dhs@( 0 , [))`]}~ #)) -: (mp     (ungrz~ <:@c)))~ 0 ?@$~ 0 _1     + $) RZf
NB.   ((((}:"1@unmr3rc ,. &0)~ (idmat@[`(a: <@; liso4dhs@( 0 , [))`]}~ #)) -: (mp  ct@(ungrz~ <:@c)))~ 0 ?@$~ 0 _1     + $) RZf
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - unmr3xx implement LAPACK's DORMR3, ZUNMR3
NB. - unm{l3,3l,3r,r3}{ln,lc,rn,rc} and
NB.   unm{lz,zl,zr,rz}{ln,lc,rn,rc} respectively are
NB.   topologic equivalents

unml3ln=: ((larzlcfr&:>/@,~ (-@c <\ ,)@|.   )~ <)^:(0 < #@[)
unml3lc=: ((larzlnfr&:>/@,~ (-@c <\ ,)      )~ <)^:(0 < #@[)
unml3rn=: ((larzrcfr&:>/@,~ (-@c <\ ,)      )~ <)^:(0 < #@[)
unml3rc=: ((larzrnfr&:>/@,~ (-@c <\ ,)@|.   )~ <)^:(0 < #@[)

unm3lln=: ((larzlnbc&:>/@,~ (-@c <\ ,)@|.@|:)~ <)^:(0 < c@[)
unm3llc=: ((larzlcbc&:>/@,~ (-@c <\ ,)   @|:)~ <)^:(0 < c@[)
unm3lrn=: ((larzrnbc&:>/@,~ (-@c <\ ,)   @|:)~ <)^:(0 < c@[)
unm3lrc=: ((larzrcbc&:>/@,~ (-@c <\ ,)@|.@|:)~ <)^:(0 < c@[)

unm3rln=: ((larzlnfc&:>/@,~ (-@c <\ ,)   @|:)~ <)^:(0 < c@[)
unm3rlc=: ((larzlcfc&:>/@,~ (-@c <\ ,)@|.@|:)~ <)^:(0 < c@[)
unm3rrn=: ((larzrnfc&:>/@,~ (-@c <\ ,)@|.@|:)~ <)^:(0 < c@[)
unm3rrc=: ((larzrcfc&:>/@,~ (-@c <\ ,)   @|:)~ <)^:(0 < c@[)

unmr3ln=: ((larzlcbr&:>/@,~ (-@c <\ ,)      )~ <)^:(0 < #@[)
unmr3lc=: ((larzlnbr&:>/@,~ (-@c <\ ,)@|.   )~ <)^:(0 < #@[)
unmr3rn=: ((larzrcbr&:>/@,~ (-@c <\ ,)@|.   )~ <)^:(0 < #@[)
unmr3rc=: ((larzrnbr&:>/@,~ (-@c <\ ,)      )~ <)^:(0 < #@[)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Verb       Action     Side     Tran    Syntax
NB. unmlqln    Q   * C    left     none    B=. LQf unmlqln C
NB. unmlqlc    Q^H * C    left     ct      B=. LQf unmlqlc C
NB. unmlqrn    C * Q      right    none    B=. LQf unmlqrn C
NB. unmlqrc    C * Q^H    right    ct      B=. LQf unmlqrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form LQf as returned by gelqf
NB. where
NB.   B,C - m×n-matrices
NB.   LQf - n×(m+1)-matrix (ln,lc cases) or m×(n+1)-matrix
NB.         (rn,rc), contains Qf, the output of gelqf
NB.   Qf  - k×(m+1)-matrix (ln,lc) or k×(n+1)-matrix (rn,rc),
NB.         unit upper trapezoidal (unit diagonal is not
NB.         stored), represents the Q in factored form
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as the product of
NB.         k elementary reflectors H(i):
NB.           Q = Π{H(i)',i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
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
NB. - unml2{ln,lc,rn,rc} and unmlq{ln,lc,rn,rc} respectively
NB.   are topologic equivalents

unmlqln=: }:  @(((unml2ln`((larfblcfr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,  &0)`(i.@$@])@.(0 e. $@])
unmlqlc=: }:  @(((unml2lc`((larfblnfr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,  &0)`(i.@$@])@.(0 e. $@])
unmlqrn=: }:"1@(((unml2rn`((larfbrcfr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,. &0)`(i.@$@])@.(0 e. $@])
unmlqrc=: }:"1@(((unml2rc`((larfbrnfr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~  tru1        @({.  ~  0 _1    <./ @:+ $))~ ,. &0)`(i.@$@])@.(0 e. $@])

NB. ---------------------------------------------------------
NB. Verb       Action     Side     Tran    Syntax
NB. unmqlln    Q   * C    left     none    B=. QfL unmqlln C
NB. unmqllc    Q^H * C    left     ct      B=. QfL unmqllc C
NB. unmqlrn    C * Q      right    none    B=. QfL unmqlrn C
NB. unmqlrc    C * Q^H    right    ct      B=. QfL unmqlrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form QfL as returned by geqlf
NB. where
NB.   B,C - m×n-matrices
NB.   QfL - (m+1)×n-matrix (ln,lc cases), (n+1)×m-matrix
NB.         (rn,rc), contains Qf, the output of geqlf
NB.   Qf  - (m+1)×k-matrix (ln,lc) or (n+1)×k-matrix (rn,rc),
NB.         unit upper trapezoidal (unit diagonal is not
NB.         stored), represents the Q in factored form
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as the product of
NB.         k elementary reflectors H(i):
NB.           Q = Π{H(i),i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
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
NB. - unm2l{ln,lc,rn,rc} and unmql{ln,lc,rn,rc} respectively
NB.   are topologic equivalents

unmqlln=: }.  @(((unm2lln`((larfblnbc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ , ~&0)`(i.@$@])@.(0 e. $@])
unmqllc=: }.  @(((unm2llc`((larfblcbc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ , ~&0)`(i.@$@])@.(0 e. $@])
unmqlrn=: }."1@(((unm2lrn`((larfbrnbc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ ,.~&0)`(i.@$@])@.(0 e. $@])
unmqlrc=: }."1@(((unm2lrc`((larfbrcbc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (tru1~ -~/@$)@({."1~ _1  0 -@(<./)@:+ $))~ ,.~&0)`(i.@$@])@.(0 e. $@])

NB. ---------------------------------------------------------
NB. Verb       Action     Side     Tran    Syntax
NB. unmqrln    Q   * C    left     none    B=. QfR unmqrln C
NB. unmqrlc    Q^H * C    left     ct      B=. QfR unmqrlc C
NB. unmqrrn    C * Q      right    none    B=. QfR unmqrrn C
NB. unmqrrc    C * Q^H    right    ct      B=. QfR unmqrrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form QfR as returned by geqrf
NB. where
NB.   B,C - m×n-matrices
NB.   QfR - (m+1)×n-matrix (ln,lc cases), (n+1)×m-matrix
NB.         (rn,rc), contains Qf, the output of geqrf
NB.   Qf  - (m+1)×k-matrix (ln,lc) or (n+1)×k-matrix (rn,rc),
NB.         unit lower trapezoidal (unit diagonal is not
NB.         stored), represents the Q in factored form
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as the product of
NB.         k elementary reflectors H(i):
NB.           Q = Π{H(i),i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
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
NB. - unm2r{ln,lc,rn,rc} and unmqr{ln,lc,rn,rc} respectively
NB.   are topologic equivalents

unmqrln=: }:  @(((unm2rln`((larfblnfc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,  &0)`(i.@$@])@.(0 e. $@])
unmqrlc=: }:  @(((unm2rlc`((larfblcfc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,  &0)`(i.@$@])@.(0 e. $@])
unmqrrn=: }:"1@(((unm2rrn`((larfbrnfc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,. &0)`(i.@$@])@.(0 e. $@])
unmqrrc=: }:"1@(((unm2rrc`((larfbrcfc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~  trl1        @({."1~ _1  0    <./ @:+ $))~ ,. &0)`(i.@$@])@.(0 e. $@])

NB. ---------------------------------------------------------
NB. Verb       Action     Side     Tran    Syntax
NB. unmrqln    Q   * C    left     none    B=. RQf unmrqln C
NB. unmrqlc    Q^H * C    left     ct      B=. RQf unmrqlc C
NB. unmrqrn    C * Q      right    none    B=. RQf unmrqrn C
NB. unmrqrc    C * Q^H    right    ct      B=. RQf unmrqrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form RQf as returned by gerqf
NB. where
NB.   B,C - m×n-matrices
NB.   RQf - l×(m+1)-matrix (ln,lc cases) or l×(n+1)-matrix
NB.         (rn,rc), contains Qf, the output of gerqf
NB.   Qf  - k×(m+1)-matrix (ln,lc) or k×(n+1)-matrix (rn,rc),
NB.         unit lower trapezoidal (unit diagonal is not
NB.         stored), represents the Q in factored form
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as the product of
NB.         k elementary reflectors H(i):
NB.           Q = Π{H(i)',i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
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
NB. - unmr2{ln,lc,rn,rc} and unmrq{ln,lc,rn,rc} respectively
NB.   are topologic equivalents

unmrqln=: }.  @(((unmr2ln`((larfblcbr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ , ~&0)`(i.@$@])@.(0 e. $@])
unmrqlc=: }.  @(((unmr2lc`((larfblnbr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ , ~&0)`(i.@$@])@.(0 e. $@])
unmrqrn=: }."1@(((unmr2rn`((larfbrcbr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ ,.~&0)`(i.@$@])@.(0 e. $@])
unmrqrc=: }."1@(((unmr2rc`((larfbrnbr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (trl1~ -~/@$)@({.  ~  0 _1 -@(<./)@:+ $))~ ,.~&0)`(i.@$@])@.(0 e. $@])

NB. ---------------------------------------------------------
NB. Verb       Action     Side     Tran    Syntax
NB. unmlzln    Z   * C    left     none    B=. LZf unmlzln C
NB. unmlzlc    Z^H * C    left     ct      B=. LZf unmlzlc C
NB. unmlzrn    C * Z      right    none    B=. LZf unmlzrn C
NB. unmlzrc    C * Z^H    right    ct      B=. LZf unmlzrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Z, which is
NB.   represented in factored form LZf as returned by tzlzf
NB. where
NB.   B,C - m×n-matrices
NB.   LZf - k×(m+1)-matrix (ln,lc cases) or k×(n+1)-matrix
NB.         (rn,rc), contains Zf (identity submatrix is not
NB.         stored), the output of tzlzf
NB.   Zf  - k×(m+1)-matrix (ln,lc) or k×(n+1)-matrix (rn,rc),
NB.         trailing k×k-submatrix is identity, the Z
NB.         represented in factored form
NB.   Z   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as the product of
NB.         k elementary reflectors H(i):
NB.           Z = Π{H(i)',i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   k   ≤ min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((unmlzln -: (mp~    (unglz~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) LZf
NB.   ((unmlzlc -: (mp~ ct@(unglz~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) LZf
NB.   ((unmlzrn -: (mp     (unglz~ <:@c))~) 0 ?@$~ 0 _1     + $) LZf
NB.   ((unmlzrc -: (mp  ct@(unglz~ <:@c))~) 0 ?@$~ 0 _1     + $) LZf
NB.
NB. Notes:
NB. - unml3{ln,lc,rn,rc} and unmlz{ln,lc,rn,rc} respectively
NB.   are topologic equivalents
NB. - if not all reflectors are needed then part of LZf
NB.   would be zeroed, e.g.:
NB.     original LZf with m=5, n=9:
NB.       (  τ0 v0 v0 v0 v0 l  0  0  0  0  )
NB.       (  τ1 v1 v1 v1 v1 l  l  0  0  0  )
NB.       (  τ2 v2 v2 v2 v2 l  l  l  0  0  )
NB.       (  τ3 v3 v3 v3 v3 l  l  l  l  0  )
NB.       (  τ4 v4 v4 v4 v4 l  l  l  l  l  )
NB.     LZf used as x argument when k=2:
NB.       (  τ3 v3 v3 v3 v3 0  0  0  l  0  )
NB.       (  τ4 v4 v4 v4 v4 0  0  0  l  l  )
NB.     note zeroed elements in columns 5,6,7

unmlzln=: }.  @(((unml3ln`((larzblcfr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; liso4dhs@(_1 , [))`]}~ #))~ , ~&0)`(i.@$@])@.(0 e. $@])
unmlzlc=: }.  @(((unml3lc`((larzblnfr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; liso4dhs@(_1 , [))`]}~ #))~ , ~&0)`(i.@$@])@.(0 e. $@])
unmlzrn=: }."1@(((unml3rn`((larzbrcfr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; liso4dhs@(_1 , [))`]}~ #))~ ,.~&0)`(i.@$@])@.(0 e. $@])
unmlzrc=: }."1@(((unml3rc`((larzbrnfr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; liso4dhs@(_1 , [))`]}~ #))~ ,.~&0)`(i.@$@])@.(0 e. $@])

NB. ---------------------------------------------------------
NB. Verb       Action     Side     Tran    Syntax
NB. unmzlln    Z   * C    left     none    B=. ZfL unmzlln C
NB. unmzllc    Z^H * C    left     ct      B=. ZfL unmzllc C
NB. unmzlrn    C * Z      right    none    B=. ZfL unmzlrn C
NB. unmzlrc    C * Z^H    right    ct      B=. ZfL unmzlrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Z, which is
NB.   represented in factored form ZfL as returned by tzzlf
NB. where
NB.   B,C - m×n-matrices
NB.   ZfL - (m+1)×k-matrix (ln,lc cases) or (n+1)×k-matrix
NB.         (rn,rc), contains Zf (identity submatrix is not
NB.         stored), the output of tzzlf
NB.   Zf  - (m+1)×k-matrix (ln,lc) or (n+1)×k-matrix (rn,rc),
NB.         leading k×k-submatrix is identity, the Z
NB.         represented in factored form
NB.   Z   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as the product of
NB.         k elementary reflectors H(i):
NB.           Z = Π{H(i),i=k-1:0}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   k   ≤ min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((unmzlln -: (mp~    (ungzl~ <:@#))~) 0 ?@$~ _1 0     + $) ZfL
NB.   ((unmzllc -: (mp~ ct@(ungzl~ <:@#))~) 0 ?@$~ _1 0     + $) ZfL
NB.   ((unmzlrn -: (mp     (ungzl~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) ZfL
NB.   ((unmzlrc -: (mp  ct@(ungzl~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) ZfL
NB.
NB. Notes:
NB. - unm3l{ln,lc,rn,rc} and unmzl{ln,lc,rn,rc} respectively
NB.   are topologic equivalents
NB. - if not all reflectors are needed then part of ZfL
NB.   would be zeroed, e.g.:
NB.     original ZfL with m=9, n=5:    ZfL used as x argument when k=2:
NB.       (  l  0  0  0  0   )           (  l  0   )
NB.       (  l  l  0  0  0   )           (  l  l   )
NB.       (  l  l  l  0  0   )           (  0  0   )
NB.       (  l  l  l  l  0   )           (  0  0   )
NB.       (  l  l  l  l  l   )           (  0  0   )
NB.       (  v0 v1 v2 v3 v4  )           (  v0 v1  )
NB.       (  v0 v1 v2 v3 v4  )           (  v0 v1  )
NB.       (  v0 v1 v2 v3 v4  )           (  v0 v1  )
NB.       (  v0 v1 v2 v3 v4  )           (  v0 v1  )
NB.       (  τ0 τ1 τ2 τ3 τ4  )           (  τ0 τ1  )
NB.                                    note zeroed elements in rows 2,3,4

unmzlln=: }:  @(((unm3lln`((larzblnbc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       liso4dhs@( 0 , [))`]}~ c))~ ,  &0)`(i.@$@])@.(0 e. $@])
unmzllc=: }:  @(((unm3llc`((larzblcbc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       liso4dhs@( 0 , [))`]}~ c))~ ,  &0)`(i.@$@])@.(0 e. $@])
unmzlrn=: }:"1@(((unm3lrn`((larzbrnbc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       liso4dhs@( 0 , [))`]}~ c))~ ,. &0)`(i.@$@])@.(0 e. $@])
unmzlrc=: }:"1@(((unm3lrc`((larzbrcbc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       liso4dhs@( 0 , [))`]}~ c))~ ,. &0)`(i.@$@])@.(0 e. $@])

NB. ---------------------------------------------------------
NB. Verb       Action     Side     Tran    Syntax
NB. unmzrln    Z   * C    left     none    B=. ZfR unmzrln C
NB. unmzrlc    Z^H * C    left     ct      B=. ZfR unmzrlc C
NB. unmzrrn    C * Z      right    none    B=. ZfR unmzrrn C
NB. unmzrrc    C * Z^H    right    ct      B=. ZfR unmzrrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Z, which is
NB.   represented in factored form ZfR as returned by tzzrf
NB. where
NB.   B,C - m×n-matrices
NB.   ZfR - (m+1)×k-matrix (ln,lc cases) or (n+1)×k-matrix
NB.         (rn,rc), contains Zf (identity submatrix is not
NB.         stored), the output of tzzrf
NB.   Zf  - (m+1)×k-matrix (ln,lc) or (n+1)×k-matrix (rn,rc),
NB.         trailing k×k-submatrix is identity, the Z
NB.         represented in factored form
NB.   Z   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as the product of
NB.         k elementary reflectors H(i):
NB.           Z = Π{H(i),i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i) * τ(i) * u(i)'
NB.   k   ≤ min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((unmzrln -: (mp~    (ungzr~ <:@#))~) 0 ?@$~ _1 0     + $) ZfR
NB.   ((unmzrlc -: (mp~ ct@(ungzr~ <:@#))~) 0 ?@$~ _1 0     + $) ZfR
NB.   ((unmzrrn -: (mp     (ungzr~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) ZfR
NB.   ((unmzrrc -: (mp  ct@(ungzr~ <:@#))~) 0 ?@$~ _1 0 |.@:+ $) ZfR
NB.
NB. Notes:
NB. - unm3r{ln,lc,rn,rc} and unmzr{ln,lc,rn,rc} respectively
NB.   are topologic equivalents
NB. - if not all reflectors are needed then part of ZfR
NB.   would be zeroed, e.g.:
NB.     original ZfR with m=9, n=5:    ZfR used as x argument when k=2:
NB.       (  τ0 τ1 τ2 τ3 τ4  )           (  τ3 τ4  )
NB.       (  v0 v1 v2 v3 v4  )           (  v3 v4  )
NB.       (  v0 v1 v2 v3 v4  )           (  v3 v4  )
NB.       (  v0 v1 v2 v3 v4  )           (  v3 v4  )
NB.       (  v0 v1 v2 v3 v4  )           (  v3 v4  )
NB.       (  r  r  r  r  r   )           (  0  0   )
NB.       (  0  r  r  r  r   )           (  0  0   )
NB.       (  0  0  r  r  r   )           (  0  0   )
NB.       (  0  0  0  r  r   )           (  r  r   )
NB.       (  0  0  0  0  r   )           (  0  r   )
NB.                                    note zeroed elements in rows 5,6,7

unmzrln=: }.  @(((unm3rln`((larzblnfc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       liso4dhs@(_1 , [))`]}~ c))~ , ~&0)`(i.@$@])@.(0 e. $@])
unmzrlc=: }.  @(((unm3rlc`((larzblcfc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       liso4dhs@(_1 , [))`]}~ c))~ , ~&0)`(i.@$@])@.(0 e. $@])
unmzrrn=: }."1@(((unm3rrn`((larzbrnfc&:>/@,~ |.@  {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       liso4dhs@(_1 , [))`]}~ c))~ ,.~&0)`(i.@$@])@.(0 e. $@])
unmzrrc=: }."1@(((unm3rrc`((larzbrcfc&:>/@,~      {.   @(<;.3~ ,:~@(MQNB ,~ #)))~ <)@.(MQNB < c@[))~ (idmat@[`(       liso4dhs@(_1 , [))`]}~ c))~ ,.~&0)`(i.@$@])@.(0 e. $@])

NB. ---------------------------------------------------------
NB. Verb       Action     Side     Tran    Syntax
NB. unmrzln    Z   * C    left     none    B=. RZf unmrzln C
NB. unmrzlc    Z^H * C    left     ct      B=. RZf unmrzlc C
NB. unmrzrn    C * Z      right    none    B=. RZf unmrzrn C
NB. unmrzrc    C * Z^H    right    ct      B=. RZf unmrzrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Z, which is
NB.   represented in factored form RZf as returned by tzrzf
NB. where
NB.   B,C - m×n-matrices
NB.   RZf - k×(m+1)-matrix (ln,lc cases) or k×(n+1)-matrix
NB.         (rn,rc), contains Zf (identity submatrix is not
NB.         stored), the output of tzrzf
NB.   Zf  - k×(m+1)-matrix (ln,lc) or k×(n+1)-matrix (rn,rc),
NB.         leading k×k-submatrix is identity, the Z
NB.         represented in factored form
NB.   Z   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), unitary
NB.         (orthogonal), which is defined as the product of
NB.         k elementary reflectors H(i):
NB.           Z = Π{H(i)',i=0:k-1}
NB.           H(i) ≡ H(u(i),τ(i)) := I - u(i)' * τ(i) * u(i)
NB.   k   ≤ min(m,n)
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   ((unmrzln -: (mp~    (ungrz~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) RZf
NB.   ((unmrzlc -: (mp~ ct@(ungrz~ <:@c))~) 0 ?@$~ 0 _1 |.@:+ $) RZf
NB.   ((unmrzrn -: (mp     (ungrz~ <:@c))~) 0 ?@$~ 0 _1     + $) RZf
NB.   ((unmrzrc -: (mp  ct@(ungrz~ <:@c))~) 0 ?@$~ 0 _1     + $) RZf
NB.
NB. Notes:
NB. - unmr3{ln,lc,rn,rc} and unmrz{ln,lc,rn,rc} respectively
NB.   are topologic equivalents
NB. - implement LAPACK's DORMRZ, ZUNMRZ with the following
NB.   difference: if not all reflectors are needed then
NB.   part of RZf would be zeroed, e.g.:
NB.     original RZf with m=5, n=9:
NB.       (  r  r  r  r  r v0 v0 v0 v0 τ0  )
NB.       (  0  r  r  r  r v1 v1 v1 v1 τ1  )
NB.       (  0  0  r  r  r v2 v2 v2 v2 τ2  )
NB.       (  0  0  0  r  r v3 v3 v3 v3 τ3  )
NB.       (  0  0  0  0  r v4 v4 v4 v4 τ4  )
NB.     RZf used as x argument when k=2:
NB.       (  r  r  0  0  0 v0 v0 v0 v0 τ0  )
NB.       (  0  r  0  0  0 v1 v1 v1 v1 τ1  )
NB.     note zeroed elements in columns 2,3,4

unmrzln=: }:  @(((unmr3ln`((larzblcbr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; liso4dhs@( 0 , [))`]}~ #))~ ,  &0)`(i.@$@])@.(0 e. $@])
unmrzlc=: }:  @(((unmr3lc`((larzblnbr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; liso4dhs@( 0 , [))`]}~ #))~ ,  &0)`(i.@$@])@.(0 e. $@])
unmrzrn=: }:"1@(((unmr3rn`((larzbrcbr&:>/@,~ |.@:({."1)@(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; liso4dhs@( 0 , [))`]}~ #))~ ,. &0)`(i.@$@])@.(0 e. $@])
unmrzrc=: }:"1@(((unmr3rc`((larzbrnbr&:>/@,~      {."1 @(<;.3~ ,:~@(MQNB ,  c)))~ <)@.(MQNB < #@[))~ (idmat@[`(a: <@; liso4dhs@( 0 , [))`]}~ #))~ ,. &0)`(i.@$@])@.(0 e. $@])

NB. ---------------------------------------------------------
NB. Verb        Action     Side     Tran    Syntax
NB. unmhrlln    Q   * C    left     none    B=. HQf unmhrlln C
NB. unmhrllc    Q^H * C    left     ct      B=. HQf unmhrllc C
NB. unmhrlrn    C * Q      right    none    B=. HQf unmhrlrn C
NB. unmhrlrc    C * Q^H    right    ct      B=. HQf unmhrlrc C
NB.
NB. Description:
NB.   Multiply a general matrix C by unitary (orthogonal)
NB.   matrix Q, which is represented in factored form HQf as
NB.   returned by gehrdl
NB. where
NB.   B,C - m×n-matrices
NB.   HQf - m×(m+1)-matrix (ln,lc cases) or n×(n+1)-matrix
NB.         (rn,rc), contains Qf, the output of gehrdl
NB.   Qf  - (s-1)×(m-h)-matrix (ln,lc) or (s-1)×(n-h)-matrix
NB.         (rn,rc), the unit upper trapezoidal (unit
NB.         diagonal is not stored), represents the Q in
NB.         factored form, located in HQf[h:h+s-2,h+1:end]
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), is the
NB.         unit matrix with unitary (orthogonal) matrix
NB.         inserted into elements Q[h:h+s-1,h:h+s-1] :
NB.           Q = Π{H(i)',i=h+s-2:h}
NB.           H(i) = I - v[i]' * τ[i] * v[i]
NB.   hs  - 2-vector of integers (h,s) 'head' and 'size',
NB.         defines submatrix Qf position in matrix HQf, see
NB.         see gehrdl
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (idmat@c -: clean@(unmhrlln ct@unghrl)@(gehrdl~ 0 , c)) A
NB.   (idmat@c -: clean@(unmhrllc    unghrl)@(gehrdl~ 0 , c)) A
NB.   (idmat@c -: clean@(unmhrlrn ct@unghrl)@(gehrdl~ 0 , c)) A
NB.   (idmat@c -: clean@(unmhrlrc    unghrl)@(gehrdl~ 0 , c)) A
NB.
NB. Notes:
NB. - instead of using f and s arguments, the following
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
NB. Verb        Action     Side     Tran    Syntax
NB. unmhruln    Q   * C    left     none    B=. HQf unmhruln C
NB. unmhrulc    Q^H * C    left     ct      B=. HQf unmhrulc C
NB. unmhrurn    C * Q      right    none    B=. HQf unmhrurn C
NB. unmhrurc    C * Q^H    right    ct      B=. HQf unmhrurc C
NB.
NB. Description:
NB.   Multiply a general matrix C by unitary (orthogonal)
NB.   matrix Q, which is represented in factored form HQf as
NB.   returned by gehrdu
NB. where
NB.   B,C - m×n-matrices
NB.   HQf - (m+1)×m-matrix (ln,lc cases) or (n+1)×n-matrix
NB.         (rn,rc), contains Qf, the output of gehrdu
NB.   Qf  - (m-h)×(s-1)-matrix (ln,lc) or (n-h)×(s-1)-matrix
NB.         (rn,rc), the unit lower trapezoidal (unit
NB.         diagonal is not stored), represents the Q in
NB.         factored form, located in HQf[h+1:end,h:h+s-2]
NB.   Q   - m×m-matrix (ln,lc) or n×n-matrix (rn,rc), is the
NB.         unit matrix with unitary (orthogonal) matrix
NB.         inserted into elements Q[h:h+s-1,h:h+s-1] :
NB.           Q = Π{H(i),i=h:h+s-2}
NB.           H(i) = I - v[i] * τ[i] * v[i]'
NB.   hs  - 2-vector of integers (h,s) 'head' and 'size',
NB.         defines submatrix Qf position in matrix HQf, see
NB.         see gehrdu
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   (idmat@# -: clean@(unmhruln ct@unghru)@(gehrdu~ 0 , #)) A
NB.   (idmat@# -: clean@(unmhrulc    unghru)@(gehrdu~ 0 , #)) A
NB.   (idmat@# -: clean@(unmhrurn ct@unghru)@(gehrdu~ 0 , #)) A
NB.   (idmat@# -: clean@(unmhrurc    unghru)@(gehrdu~ 0 , #)) A
NB.
NB. Notes:
NB. - models LAPACK's DORMHR, ZUNMHR
NB. - instead of using f and s arguments, the following
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
NB.   Test:
NB.   - DORMxx ZUNMxx (math/lapack2 addon)
NB.   - unmxxxx (math/mt addon)
NB.   by general matrices
NB.
NB. Syntax:
NB.   log=. testunmq (A ; C)
NB. where
NB.   A   - m×n-matrix, is used to produce Qf
NB.   C   - m×n-matrix, is used as multiplier
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.
NB. Notes:
NB. - LAPACK's DORMxQ and ZUNMxQ requires A to have at least
NB.   1 row, so Head ({.) is used when (# A) equals 0

testunmq=: 3 : 0
  load_mttmp_ 'math/mt/external/lapack2/dormlq'
  load_mttmp_ 'math/mt/external/lapack2/dormql'
  load_mttmp_ 'math/mt/external/lapack2/dormqr'
  load_mttmp_ 'math/mt/external/lapack2/dormrq'
  load_mttmp_ 'math/mt/external/lapack2/zunmlq'
  load_mttmp_ 'math/mt/external/lapack2/zunmql'
  load_mttmp_ 'math/mt/external/lapack2/zunmqr'
  load_mttmp_ 'math/mt/external/lapack2/zunmrq'

  'A C'=. y

  rcond=. nan`geconi@.(=/@$) C  NB. meaninigful for square matrices only

  ks=. ~. 0 1 , (,~ <.@-:) <./ 'm n'=. $ A

  Awide=. |:^:(>/@$) A
  Atall=. |:^:(</@$) A

  normw=. norm1 Cwide=. |:^:(>/@$) C
  normt=. norm1 Ctall=. |:^:(</@$) C

  LQf=. gelqf Awide
  QfL=. geqlf Atall
  QfR=. geqrf Atall
  RQf=. gerqf Awide

  log=. nolog ''

  NB. LAPACK, real datatype

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''ln''&dormlq_mttmp_' tmonad (((   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ ] ) lqt03))) Ctall ; normt ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''lt''&dormlq_mttmp_' tmonad (((   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ ct) lqt03))) Ctall ; normt ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rn''&dormlq_mttmp_' tmonad (((   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  ] ) lqt03))) Cwide ; normw ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rt''&dormlq_mttmp_' tmonad (((   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  ct) lqt03))) Cwide ; normw ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''ln''&dormql_mttmp_' tmonad (((-@(3&{::) (}.   ; {.  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp~ ] ) qlt03))) Ctall ; normt ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''lt''&dormql_mttmp_' tmonad (((-@(3&{::) (}.   ; {.  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp~ ct) qlt03))) Ctall ; normt ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rn''&dormql_mttmp_' tmonad (((-@(3&{::) (}.   ; {.  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp  ] ) qlt03))) Cwide ; normw ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rt''&dormql_mttmp_' tmonad (((-@(3&{::) (}.   ; {.  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp  ct) qlt03))) Cwide ; normw ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''ln''&dormqr_mttmp_' tmonad (((   3&{::  (}:   ; {:  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp~ ] ) qrt03))) Ctall ; normt ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''lt''&dormqr_mttmp_' tmonad (((   3&{::  (}:   ; {:  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp~ ct) qrt03))) Ctall ; normt ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rn''&dormqr_mttmp_' tmonad (((   3&{::  (}:   ; {:  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp  ] ) qrt03))) Cwide ; normw ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rt''&dormqr_mttmp_' tmonad (((   3&{::  (}:   ; {:  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp  ct) qrt03))) Cwide ; normw ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''ln''&dormrq_mttmp_' tmonad (((-@(3&{::) (}."1 ; {."1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ ] ) rqt03))) Ctall ; normt ; RQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''lt''&dormrq_mttmp_' tmonad (((-@(3&{::) (}."1 ; {."1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ ct) rqt03))) Ctall ; normt ; RQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rn''&dormrq_mttmp_' tmonad (((-@(3&{::) (}."1 ; {."1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  ] ) rqt03))) Cwide ; normw ; RQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rt''&dormrq_mttmp_' tmonad (((-@(3&{::) (}."1 ; {."1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  ct) rqt03))) Cwide ; normw ; RQf ; ik { ks
    ik=. >: ik
  end.

  NB. LAPACK, complex datatype

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''ln''&zunmlq_mttmp_' tmonad (((   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ ] ) lqt03))) Ctall ; normt ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''lc''&zunmlq_mttmp_' tmonad (((   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ ct) lqt03))) Ctall ; normt ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rn''&zunmlq_mttmp_' tmonad (((   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  ] ) lqt03))) Cwide ; normw ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rc''&zunmlq_mttmp_' tmonad (((   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  ct) lqt03))) Cwide ; normw ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''ln''&zunmql_mttmp_' tmonad (((-@(3&{::) (}.   ; {.  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp~ ] ) qlt03))) Ctall ; normt ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''lc''&zunmql_mttmp_' tmonad (((-@(3&{::) (}.   ; {.  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp~ ct) qlt03))) Ctall ; normt ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rn''&zunmql_mttmp_' tmonad (((-@(3&{::) (}.   ; {.  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp  ] ) qlt03))) Cwide ; normw ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rc''&zunmql_mttmp_' tmonad (((-@(3&{::) (}.   ; {.  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp  ct) qlt03))) Cwide ; normw ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''ln''&zunmqr_mttmp_' tmonad (((   3&{::  (}:   ; {:  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp~ ] ) qrt03))) Ctall ; normt ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''lc''&zunmqr_mttmp_' tmonad (((   3&{::  (}:   ; {:  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp~ ct) qrt03))) Ctall ; normt ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rn''&zunmqr_mttmp_' tmonad (((   3&{::  (}:   ; {:  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp  ] ) qrt03))) Cwide ; normw ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rc''&zunmqr_mttmp_' tmonad (((   3&{::  (}:   ; {:  )@:({."1) 2&{::) , {.)`]`(rcond"_)`nan`((mp  ct) qrt03))) Cwide ; normw ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''ln''&zunmrq_mttmp_' tmonad (((-@(3&{::) (}."1 ; {."1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ ] ) rqt03))) Ctall ; normt ; RQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''lc''&zunmrq_mttmp_' tmonad (((-@(3&{::) (}."1 ; {."1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ ct) rqt03))) Ctall ; normt ; RQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rn''&zunmrq_mttmp_' tmonad (((-@(3&{::) (}."1 ; {."1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  ] ) rqt03))) Cwide ; normw ; RQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rc''&zunmrq_mttmp_' tmonad (((-@(3&{::) (}."1 ; {."1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  ct) rqt03))) Cwide ; normw ; RQf ; ik { ks
    ik=. >: ik
  end.

  NB. mt, any datatype

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmlqln' tdyad ((   3&{::   {.    2&{::)`(0&{::)`]`(rcond"_)`nan`((mp~ ] ) lqt03))) Ctall ; normt ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmlqlc' tdyad ((   3&{::   {.    2&{::)`(0&{::)`]`(rcond"_)`nan`((mp~ ct) lqt03))) Ctall ; normt ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmlqrn' tdyad ((   3&{::   {.    2&{::)`(0&{::)`]`(rcond"_)`nan`((mp  ] ) lqt03))) Cwide ; normw ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmlqrc' tdyad ((   3&{::   {.    2&{::)`(0&{::)`]`(rcond"_)`nan`((mp  ct) lqt03))) Cwide ; normw ; LQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmqlln' tdyad ((-@(3&{::) ({."1) 2&{::)`(0&{::)`]`(rcond"_)`nan`((mp~ ] ) qlt03))) Ctall ; normt ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmqllc' tdyad ((-@(3&{::) ({."1) 2&{::)`(0&{::)`]`(rcond"_)`nan`((mp~ ct) qlt03))) Ctall ; normt ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmqlrn' tdyad ((-@(3&{::) ({."1) 2&{::)`(0&{::)`]`(rcond"_)`nan`((mp  ] ) qlt03))) Cwide ; normw ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmqlrc' tdyad ((-@(3&{::) ({."1) 2&{::)`(0&{::)`]`(rcond"_)`nan`((mp  ct) qlt03))) Cwide ; normw ; QfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmqrln' tdyad ((   3&{::  ({."1) 2&{::)`(0&{::)`]`(rcond"_)`nan`((mp~ ] ) qrt03))) Ctall ; normt ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmqrlc' tdyad ((   3&{::  ({."1) 2&{::)`(0&{::)`]`(rcond"_)`nan`((mp~ ct) qrt03))) Ctall ; normt ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmqrrn' tdyad ((   3&{::  ({."1) 2&{::)`(0&{::)`]`(rcond"_)`nan`((mp  ] ) qrt03))) Cwide ; normw ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmqrrc' tdyad ((   3&{::  ({."1) 2&{::)`(0&{::)`]`(rcond"_)`nan`((mp  ct) qrt03))) Cwide ; normw ; QfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmrqln' tdyad ((-@(3&{::)  {.    2&{::)`(0&{::)`]`(rcond"_)`nan`((mp~ ] ) rqt03))) Ctall ; normt ; RQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmrqlc' tdyad ((-@(3&{::)  {.    2&{::)`(0&{::)`]`(rcond"_)`nan`((mp~ ct) rqt03))) Ctall ; normt ; RQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmrqrn' tdyad ((-@(3&{::)  {.    2&{::)`(0&{::)`]`(rcond"_)`nan`((mp  ] ) rqt03))) Cwide ; normw ; RQf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmrqrc' tdyad ((-@(3&{::)  {.    2&{::)`(0&{::)`]`(rcond"_)`nan`((mp  ct) rqt03))) Cwide ; normw ; RQf ; ik { ks
    ik=. >: ik
  end.

  coerase < 'mttmp'

  log
)

NB. ---------------------------------------------------------
NB. testunmz
NB.
NB. Description:
NB.   Test:
NB.   - DORMRZ ZUNMRZ (math/lapack2 addon)
NB.   - unmxxxx (math/mt addon)
NB.   by trapezoidal and general matrices
NB.
NB. Syntax:
NB.   log=. testunmz (A ; C)
NB. where
NB.   A   - m×n-matrix, is used to produce Zf
NB.   C   - m×n-matrix, is used as multiplier
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.
NB. Notes:
NB. - LAPACK's DORMRZ and ZUNMRZ requires A to have at least
NB.   1 row, so Head ({.) is used when (# A) equals 0

testunmz=: 3 : 0
  load_mttmp_ 'math/mt/external/lapack2/tzrzf'
  load_mttmp_ 'math/mt/external/lapack2/dormrz'
  load_mttmp_ 'math/mt/external/lapack2/zunmrz'

  'A C'=. y

  rcond=. nan`geconi@.(=/@$) C  NB. meaninigful for square matrices only

  ks=. ~. 0 1 , (,~ <.@-:) <./ 'm n'=. $ A

  Awide=. |:^:(>/@$) A
  Atall=. |:^:(</@$) A

  normw=. norm1 Cwide=. |:^:(>/@$) C
  normt=. norm1 Ctall=. |:^:(</@$) C

  LZf=. tzlzf Awide
  ZfL=. tzzlf Atall
  ZfR=. tzzrf Atall

  log=. nolog ''

  NB. LAPACK, real datatype

  NB. LAPACK stores RZf differently, so we need variants for
  NB. rzt03 (check DORMRZ by DORMRZ, heh) and RZf itself

  rzt03a=: 1 : 'norm1@(- 0&{:: u ''ln'' dormrz_mttmp_ 3&{:: (<:@(-~/)@$@] ; (}:"1 ; {:"1)@{. , <@idmat@<:@c@]) 2&{::)~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. <:@c@(2 {:: [)'

  try.
    RZf=. tru (0&{:: ,. 1&{::) dtzrzf_mttmp_ Awide
  catch.
    RZf=. _.
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''ln''&dormrz_mttmp_' tmonad ((<:@(-~/)@$@(2&{::) ; (   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ ] ) rzt03a))) Ctall ; normt ; RZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''lt''&dormrz_mttmp_' tmonad ((<:@(-~/)@$@(2&{::) ; (   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ |:) rzt03a))) Ctall ; normt ; RZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rn''&dormrz_mttmp_' tmonad ((<:@(-~/)@$@(2&{::) ; (   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  ] ) rzt03a))) Cwide ; normw ; RZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rt''&dormrz_mttmp_' tmonad ((<:@(-~/)@$@(2&{::) ; (   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  |:) rzt03a))) Cwide ; normw ; RZf ; ik { ks
    ik=. >: ik
  end.

  NB. LAPACK, complex datatype

  NB. LAPACK stores RZf differently, so we need variants for
  NB. rzt03 (check ZUNMRZ by ZUNMRZ, heh) and RZf itself

  rzt03b=: 1 : 'norm1@(- 0&{:: u ''ln'' zunmrz_mttmp_ 3&{:: (<:@(-~/)@$@] ; (}:"1 ; {:"1)@{. , <@idmat@<:@c@]) 2&{::)~ % FP_EPS * 1:^:(0&=)@(1 {:: [) * 1 >. <:@c@(2 {:: [)'

  try.
    RZf=. tru (0&{:: ,. 1&{::) ztzrzf_mttmp_ Awide
  catch.
    RZf=. _.
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''ln''&zunmrz_mttmp_' tmonad ((<:@(-~/)@$@(2&{::) ; (   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ ] ) rzt03b))) Ctall ; normt ; RZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''lc''&zunmrz_mttmp_' tmonad ((<:@(-~/)@$@(2&{::) ; (   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp~ ct) rzt03b))) Ctall ; normt ; RZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rn''&zunmrz_mttmp_' tmonad ((<:@(-~/)@$@(2&{::) ; (   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  ] ) rzt03b))) Cwide ; normw ; RZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('''rc''&zunmrz_mttmp_' tmonad ((<:@(-~/)@$@(2&{::) ; (   3&{::  (}:"1 ; {:"1)@  {.    2&{::) , {.)`]`(rcond"_)`nan`((mp  ct) rzt03b))) Cwide ; normw ; RZf ; ik { ks
    ik=. >: ik
  end.

  NB. mt, any datatype

  RZf=. tzrzf Awide

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmlzln' tdyad ((-@(3&{::)  {.    (((}."1~ -) ,.  idmat@]) #)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp~ ] ) lzt03))) Ctall ; normt ; LZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmlzlc' tdyad ((-@(3&{::)  {.    (((}."1~ -) ,.  idmat@]) #)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp~ ct) lzt03))) Ctall ; normt ; LZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmlzrn' tdyad ((-@(3&{::)  {.    (((}."1~ -) ,.  idmat@]) #)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp  ] ) lzt03))) Cwide ; normw ; LZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmlzrc' tdyad ((-@(3&{::)  {.    (((}."1~ -) ,.  idmat@]) #)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp  ct) lzt03))) Cwide ; normw ; LZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmzlln' tdyad ((   3&{::  ({."1) (( }.  ~    , ~ idmat@]) c)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp~ ] ) zlt03))) Ctall ; normt ; ZfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmzllc' tdyad ((   3&{::  ({."1) (( }.  ~    , ~ idmat@]) c)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp~ ct) zlt03))) Ctall ; normt ; ZfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmzlrn' tdyad ((   3&{::  ({."1) (( }.  ~    , ~ idmat@]) c)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp  ] ) zlt03))) Cwide ; normw ; ZfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmzlrc' tdyad ((   3&{::  ({."1) (( }.  ~    , ~ idmat@]) c)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp  ct) zlt03))) Cwide ; normw ; ZfL ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmzrln' tdyad ((-@(3&{::) ({."1) (((}.  ~ -) ,   idmat@]) c)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp~ ] ) zrt03))) Ctall ; normt ; ZfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmzrlc' tdyad ((-@(3&{::) ({."1) (((}.  ~ -) ,   idmat@]) c)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp~ ct) zrt03))) Ctall ; normt ; ZfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmzrrn' tdyad ((-@(3&{::) ({."1) (((}.  ~ -) ,   idmat@]) c)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp  ] ) zrt03))) Cwide ; normw ; ZfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmzrrc' tdyad ((-@(3&{::) ({."1) (((}.  ~ -) ,   idmat@]) c)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp  ct) zrt03))) Cwide ; normw ; ZfR ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmrzln' tdyad ((   3&{::   {.    (( }."1~    ,.~ idmat@]) #)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp~ ] ) rzt03))) Ctall ; normt ; RZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmrzlc' tdyad ((   3&{::   {.    (( }."1~    ,.~ idmat@]) #)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp~ ct) rzt03))) Ctall ; normt ; RZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmrzrn' tdyad ((   3&{::   {.    (( }."1~    ,.~ idmat@]) #)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp  ] ) rzt03))) Cwide ; normw ; RZf ; ik { ks
    ik=. >: ik
  end.

  ik=. 0
  while. ik < # ks do.
    log=. log lcat ('unmrzrc' tdyad ((   3&{::   {.    (( }."1~    ,.~ idmat@]) #)@(2&{::))`(0&{::)`]`(rcond"_)`nan`((mp  ct) rzt03))) Cwide ; normw ; RZf ; ik { ks
    ik=. >: ik
  end.

  coerase < 'mttmp'
  erase 'rzt03a rzt03b'

  log
)

NB. ---------------------------------------------------------
NB. testunmhr
NB.
NB. Description:
NB.   Test:
NB.   - DORMHR ZUNMHR (math/lapack2 addon)
NB.   - unmhrxxx (math/mt addon)
NB.   by square matrices
NB.
NB. Syntax:
NB.   log=. testunmhr (A ; C)
NB. where
NB.   A   - n×n-matrix, is used to produce Qf
NB.   C   - n×n-matrix, is used as multiplier
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.
NB. Notes:
NB. - LAPACK's DORMHR and ZUNMHR requires A to have at least
NB.   1 row, so Head ({.) is used when (# A) equals 0

testunmhr=: 3 : 0
  load_mttmp_ 'math/mt/external/lapack2/gehrd'
  load_mttmp_ 'math/mt/external/lapack2/dorghr'
  load_mttmp_ 'math/mt/external/lapack2/zunghr'
  load_mttmp_ 'math/mt/external/lapack2/dormhr'
  load_mttmp_ 'math/mt/external/lapack2/zunmhr'

  'A C'=. y

  rcond=. nan`geconi@.(=/@$) C  NB. meaninigful for square matrices only

  normC=. norm1 C

  NB. LAPACK, real datatype

  NB. LAPACK stores HQf differently, so we need a HQf variant
  try.
    HuQf=. (0&{:: , 1&{::) 'HuQf2 tau2'=. dgehrd_mttmp_ (1 ; # ; ]) A
    Qu=. dorghr_mttmp_ 1 ; ((# ; 1&{.^:(0 = #)) HuQf2) , < tau2
  catch.
    HuQf=. _.
    Qu=. _.
  end.

  log=.          ('''ln''&dormhr_mttmp_' tmonad (((1 ; c ; 1&{.^:(0 = #)@:}: ; }:@{:)@(2&{::) , {.)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp~    Qu
  log=. log lcat ('''lt''&dormhr_mttmp_' tmonad (((1 ; c ; 1&{.^:(0 = #)@:}: ; }:@{:)@(2&{::) , {.)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp~ |: Qu
  log=. log lcat ('''rn''&dormhr_mttmp_' tmonad (((1 ; c ; 1&{.^:(0 = #)@:}: ; }:@{:)@(2&{::) , {.)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp     Qu
  log=. log lcat ('''rt''&dormhr_mttmp_' tmonad (((1 ; c ; 1&{.^:(0 = #)@:}: ; }:@{:)@(2&{::) , {.)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp  |: Qu

  NB. LAPACK, complex datatype

  NB. LAPACK stores HQf differently, so we need a HQf variant
  try.
    HuQf=. (0&{:: , 1&{::) 'HuQf2 tau2'=. zgehrd_mttmp_ (1 ; # ; ]) A
    Qu=. zunghr_mttmp_ 1 ; ((# ; 1&{.^:(0 = #)) HuQf2) , < tau2
  catch.
    HuQf=. _.
    Qu=. _.
  end.

  log=. log lcat ('''ln''&zunmhr_mttmp_' tmonad (((1 ; c ; 1&{.^:(0 = #)@:}: ; }:@{:)@(2&{::) , {.)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp~    Qu
  log=. log lcat ('''lc''&zunmhr_mttmp_' tmonad (((1 ; c ; 1&{.^:(0 = #)@:}: ; }:@{:)@(2&{::) , {.)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp~ ct Qu
  log=. log lcat ('''rn''&zunmhr_mttmp_' tmonad (((1 ; c ; 1&{.^:(0 = #)@:}: ; }:@{:)@(2&{::) , {.)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp     Qu
  log=. log lcat ('''rc''&zunmhr_mttmp_' tmonad (((1 ; c ; 1&{.^:(0 = #)@:}: ; }:@{:)@(2&{::) , {.)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp  ct Qu

  NB. mt, any datatype

  Ql=. unghrl HlQf=. (gehrdl~ 0 , c) A
  Qu=. unghru HuQf=. (gehrdu~ 0 , #) A

  log=. log lcat ('unmhrlln' tdyad ((2&{::)`(0&{::)`]`(rcond"_)`nan`hst03l)) C ; normC ; HlQf ; C mp~    Ql
  log=. log lcat ('unmhrllc' tdyad ((2&{::)`(0&{::)`]`(rcond"_)`nan`hst03l)) C ; normC ; HlQf ; C mp~ ct Ql
  log=. log lcat ('unmhrlrn' tdyad ((2&{::)`(0&{::)`]`(rcond"_)`nan`hst03l)) C ; normC ; HlQf ; C mp     Ql
  log=. log lcat ('unmhrlrc' tdyad ((2&{::)`(0&{::)`]`(rcond"_)`nan`hst03l)) C ; normC ; HlQf ; C mp  ct Ql

  log=. log lcat ('unmhruln' tdyad ((2&{::)`(0&{::)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp~    Qu
  log=. log lcat ('unmhrulc' tdyad ((2&{::)`(0&{::)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp~ ct Qu
  log=. log lcat ('unmhrurn' tdyad ((2&{::)`(0&{::)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp     Qu
  log=. log lcat ('unmhrurc' tdyad ((2&{::)`(0&{::)`]`(rcond"_)`nan`hst03u)) C ; normC ; HuQf ; C mp  ct Qu

  coerase < 'mttmp'

  log
)

NB. ---------------------------------------------------------
NB. testmq
NB.
NB. Description:
NB.   Adv. to make verb to test unmxxxxx by matrix of
NB.   generator and shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testmq) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testmq_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testmq_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testmq_mt_ 150 200

testmq=: 1 : '(nolog_mt_`testunmhr_mt_@.(=/@$@(0&{::)) ,&.>~ testunmz_mt_ ,&.>~ testunmq_mt_)@(u ; u)'
