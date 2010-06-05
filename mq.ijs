NB. mq.ijs
NB. Multiply by matrix with orthonormal rows or columns,
NB. which is represented in factored form
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

UNMQBS=: 32   NB. block size limit
UNMQNX=: 128  NB. crossover point, UNMQNX ≥ UNMQBS

NB. ---------------------------------------------------------
unml2lnstep=: ((,;.0)~(1 _,:~2$(0{$@(0&{::))))larflcfr(1{::])
unml2lcstep=: ((,;.0)~(1 _,:~2$(0{$@(0&{::))))larflnfr(1{::])
unml2rnstep=: ((,;.0)~(1 _,:~2$(1{$@(0&{::))))larfrcfr(1{::])
unml2rcstep=: ((,;.0)~(1 _,:~2$(1{$@(0&{::))))larfrnfr(1{::])

NB. ---------------------------------------------------------
unm2llnstep=: ((,;.0)~(_ 1,:~2$(-@>:@(0{$@(1&{::)))))larflnbc(0{::])
unm2llcstep=: ((,;.0)~(_ 1,:~2$(-@>:@(0{$@(1&{::)))))larflcbc(0{::])
unm2lrnstep=: ((,;.0)~(_ 1,:~2$(-@>:@(1{$@(1&{::)))))larfrnbc(0{::])
unm2lrcstep=: ((,;.0)~(_ 1,:~2$(-@>:@(1{$@(1&{::)))))larfrcbc(0{::])

NB. ---------------------------------------------------------
NB. Verb     Action  Side   Tran  Syntax
NB. unml2ln  Q * C   left   none  QeC=.  Qf unml2ln (C, 0)
NB. unml2lc  Q'* C   left   ct    cQeC=. Qf unml2lc (C, 0)
NB. unml2rn  C * Q   right  none  eCQ=.  Qf unml2lc (C,.0)
NB. unml2rc  C * Q'  right  ct    eCcQ=. Qf unml2lc (C,.0)
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form Qf as returned by gelqf.
NB.   This is non-blocked version of algorithm
NB. where
NB.   C      - m×n-matrix to multiply
NB.   Qf     - unit upper triangular k×(n+1)-matrix Qf,
NB.            it represents Q in factored form################################
NB.   Q      - k×n-matrix with orthonormal rows, which is
NB.            defined as the product of k elementary
NB.            reflectors of order n: Q = Π{H(i)',i=k-1:0}
NB.   eCprod - being product of matrices C and Q with
NB.            appended or stitched modified trash vector
NB.
NB. Algorithm:
NB.   1) form y-argument for (^:) verb:
NB.        y=. i0 ; eC
NB.      where i0 is IO 1st reflecting vector v in Qf, it is
NB.      either 0 or ((# Qf) - 1)
NB.   2) find number of iterations:
NB.        iters=. # Qf
NB.   3) start iterations:
NB.        y=. Qf step ^: iters y
NB.      3.1) extract reflecting vector:
NB.             vi=. rIOS (, ;. 0) eCi
NB.      3.2) reflect submatrix subeCi of eCi:
NB.             subeCupdi=. Qf larfxxfr subeCi
NB.      3.3) merge not reflected and reflected parts of eCi
NB.           to form eCi1
NB.      3.4) link next IO and eCi1 to form y-argument of the
NB.           following iteration
NB.   4) extract eCi from result of (^:) verb:
NB.        eCprod=. 1 {:: y
NB.
NB. If:
NB.   LQf=. gelqf A
NB.   L=. trl (0 _1 }. LQf)
NB.   Qf=. tru1 LQf
NB.   Q=. unglq LQf
NB. then
NB.   (Q mp (ct A)) ((-: & clean ) }:) (Qf unml2ln ((ct A) , 0))                           CLARIFYME!
NB.   (idmat n) ((-: & clean ) }:) ((tru1 gelqf A) unml2lc ((unglq tru1 gelqf A) , 0))     CLARIFYME!
NB.   clean (tru1 gelqf B) unml2rn ((ct unglq tru1 gelqf B) ,. 0)                          CLARIFYME!
NB.   clean (tru1 gelqf B) unml2rc (B ,. 0)                                                CLARIFYME!
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - emulates LAPACK's xUNML2
NB. - unml2{lc,ln,rc,rn} and unmlq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unml2ln=: [((0{::]), unml2lnstep)( (((0{::])(( ,  {.  ); ( 1  0&}.@]))unml2lnstep)^:(<:@#@[))(;~(0&{.)))
unml2lc=: [((0{::]), unml2lcstep)([(((0{::])(((, ~{:)~);~(_1  0&}.@[))unml2lcstep)^:(<:@#@[))((({.;}.)~(<:@#))~))
unml2rn=: [((0{::]),.unml2rnstep)([(((0{::])(((,.~tc)~);~( 0 _1&}.@[))unml2rnstep)^:(<:@#@[))((({.~(_&,));(}.~(0&,)))(<:@#))~)
unml2rc=: [((0{::]),.unml2rcstep)( (((0{::])(( ,. hc  ); ( 0  1&}.@]))unml2rcstep)^:(<:@#@[))(;~(_ 0&{.)))

NB. ---------------------------------------------------------
NB. Verb     Action  Side   Tran  Syntax
NB. unm2lln  Q * C   left   none  QeC=.  Qf unm2lln (0, C)
NB. unm2llc  Q'* C   left   ct    cQeC=. Qf unm2llc (0, C)
NB. unm2lrn  C * Q   right  none  eCQ=.  Qf unm2llc (0,.C)
NB. unm2lrc  C * Q'  right  ct    eCcQ=. Qf unm2llc (0,.C)
NB.
NB. Description:
NB.   Multiply a general matrix C by matrix Q, which is
NB.   represented in factored form Qf as returned by geqlf.
NB.   This is non-blocked version of algorithm
NB. where
NB.   C      - m×n-matrix to multiply###############################
NB.   Qf     - unit upper triangular k×(n+1)-matrix Qf,
NB.            it represents Q in factored form
NB.   Q      - k×n-matrix with orthonormal rows, which is
NB.            defined as the product of k elementary
NB.            reflectors of order n: Q = Π{H(i)',i=k-1:0}
NB.   eCprod - being product of matrices C and Q with
NB.            appended or stitched modified trash vector
NB.
NB. Algorithm:
NB.   1) form y-argument for (^:) verb:
NB.        y=. i0 ; eC
NB.      where i0 is IO 1st reflecting vector v in Qf, it is
NB.      either 0 or ((# Qf) - 1)
NB.   2) find number of iterations:
NB.        iters=. # Qf
NB.   3) start iterations:
NB.        y=. Qf step ^: iters y
NB.      3.1) extract reflecting vector:
NB.             vi=. rIOS (, ;. 0) eCi
NB.      3.2) reflect submatrix subeCi of eCi:
NB.             subeCupdi=. Qf larfxxfr subeCi
NB.      3.3) merge not reflected and reflected parts of eCi
NB.           to form eCi1
NB.      3.4) link next IO and eCi1 to form y-argument of the
NB.           following iteration
NB.   4) extract eCi from result of (^:) verb:
NB.        eCprod=. 1 {:: y
NB.
NB. If:
NB.   LQf=. gelqf A
NB.   L=. trl (0 _1 }. LQf)
NB.   Qf=. tru1 LQf
NB.   Q=. unglq LQf
NB. then
NB.   clean (_1 tru1 geqlf A) unm2lln ((ct ungql (11 _10 {. geqlf A)) ,~ 0)                                      CLARIFYME!
NB.   ((ct ungql (11 _10 {. geqlf A)) mp A) ((,.@; & clean ) }.) ((_1 tru1 11 _10 {. geqlf A) unm2llc (A ,~ 0))  CLARIFYME!
NB.   (idmat n) ((-: & clean ) (0 1&}.)) ((_1 tru1 geqlf A) unm2lrn ((ct ungql _1 tru1 geqlf A) ,.~ 0))          CLARIFYME!
NB.   (idmat n) ((-: & clean ) (0 1&}.)) ((_1 tru1 geqlf A) unm2lrc ((ungql _1 tru1 geqlf A) ,.~ 0))             CLARIFYME!
NB.
NB. Notes:
NB. - input's and output's shapes are the same
NB. - emulates LAPACK's xUNML2
NB. - unm2l{lc,ln,rc,rn} and unmlq{lc,ln,rc,rn} respectively
NB.   are topologic equivalents

unm2lln=: [(unm2llnstep, (1{::]))([((unm2llnstep(( ,  {.  ); ( 1  0&}.@]))(1{::]))^:(<:@(1{$)@[))((({.;~}.)~(-@<:@(1{$)))~))
unm2llc=: [(unm2llcstep, (1{::]))( ((unm2llcstep(((, ~{:)~);~(_1  0&}.@[))(1{::]))^:(<:@(1{$)@[))(; (0&{.)))
unm2lrn=: [(unm2lrnstep,.(1{::]))( ((unm2lrnstep(((,.~tc)~);~( 0 _1&}.@[))(1{::]))^:(<:@(1{$)@[))(; (_ 0&{.)))
unm2lrc=: [(unm2lrcstep,.(1{::]))([((unm2lrcstep(( ,. hc  ); ( 0  1&}.@]))(1{::]))^:(<:@(1{$)@[))((((}.~(0&,));({.~(_&,)))(-@<:@(1{$)))~))


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
