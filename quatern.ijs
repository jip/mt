NB. Quaternions
NB.
NB. qnxx       Get/set component(s)
NB. qnmarkxxx  Mark component(s)
NB. qnconxx    Conjugate component(s)
NB. qnlen      Length
NB. qnsign     Signum
NB. qninv      Inverse
NB. qnmul      Multiply
NB. qndivl     Divide (left quotient)
NB. qndivr     Divide (right quotient)
NB. qnf        Adv. to quaternificate verb
NB.
NB. Version: 0.11.0 2021-01-17
NB.
NB. Copyright 2011-2021 Igor Zhuravlov
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
NB. Concepts
NB.
NB. Notation:
NB.   ℝ  ≡ [-∞,+∞]
NB.   a,b,c,d ∊ ℝ       Basis elements multiplication table:
NB.   ℂi ≡ ℝ  + ℝ *i          i   j   k
NB.   ℂj ≡ ℝ  + ℝ *j      i  _1   k  -j
NB.   ℂk ≡ ℝ  + ℝ *k      j  -k  _1   i
NB.   ℍ  ≡ ℂi + ℂi*j      k   j  -i  _1
NB.   q  = a + b*i + c*j + d*k ∊ ℍ
NB.   q∞ - quaternion infinity
NB.
NB. Terms:
NB.   real infinity
NB.     - has a = ∞ and b,c,d = 0
NB.   imaginary infinity
NB.     - has only one from b,c,d = ∞ and "others together
NB.       with a" = 0
NB.   directed infinity
NB.     - is (q * ∞) where all q components a,b,c,d ≠ ∞
NB.     - or it has only one component = ∞ and others ≠ ∞
NB.   quaternion infinity
NB.     - has multiple components which are = ∞
NB.
NB. Conventions:
NB.   Math:                     J:
NB.     x = a + b*i ∊ ℂi          x -: a j. b
NB.     y = c + d*i ∊ ℂi          y -: c j. d
NB.     z = a + c*j ∊ ℂj          z -: a j. c
NB.     w = a + d*k ∊ ℂk          w -: a j. d
NB.     u = b + c*k ∊ ℂk          u -: b j. c
NB.     v = b + d*j ∊ ℂj          v -: b j. d
NB.     s = c + b*k ∊ ℂk          s -: c j. b
NB.     q = x + y*j               q -: x , y
NB.       = w + j*s
NB.       = w + u*i
NB.       = z + i*v
NB.     sgn(q * ∞) = sgn(q)       (qnsign q)
NB.     sgn(q∞) is undefined      throws NaN error
NB.
NB. Notes:
NB. - elements from ℝ, ℂi, ℂj, ℂk must not be mixed with each
NB.   other by math operators in J

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Get/set component(s)
NB.
NB. Verb    Syntax (monad)    Syntax (dyad)
NB. qn1     a=. qn1  q        qa=. a qn1  q
NB. qni     b=. qni  q        qb=. b qni  q
NB. qnj     c=. qnj  q        qc=. c qnj  q
NB. qnk     d=. qnk  q        qd=. d qnk  q
NB. qn1i    x=. qn1i q        qx=. x qn1i q
NB. qn1j    z=. qn1j q        qz=. z qn1j q
NB. qn1k    w=. qn1k q        qw=. w qn1k q
NB. qnij    u=. qnij q        qu=. u qnij q
NB. qnik    v=. qnik q        qv=. v qnik q
NB. qnjk    y=. qnjk q        qy=. y qnjk q

qn1=: ( 9 o. {.) : ((j.  qni) 0} ])
qni=: (11 o. {.) : ((j.~ qn1) 0} ])
qnj=: ( 9 o. {:) : ((j.  qnk) 1} ])
qnk=: (11 o. {:) : ((j.~ qnj) 1} ])

qn1i=: {.            : (0})
qnjk=: {:            : (1})
qn1j=: j./@(9   &o.) : ((j.~ +.)~ 11&o.)
qnik=: j./@(  11&o.) : ((j.  +.)~  9&o.)
qn1k=: j./@(9 11&o.) : ((< 2 2 2 $ 0 0 2 0 1 1 0 1) j./@:{ ,&:+.)
qnij=: j./@(11 9&o.) : ((< 2 2 2 $ 1 0 0 1 0 0 2 1) j./@:{ ,&:+.)

NB. ---------------------------------------------------------
NB. Mark component(s)
NB.
NB. Verb         Action                 Syntax
NB. qnmark1      a + 0*i + 0*j + 0*k    qa=.   qnmark1   q
NB. qnmarki      0 + b*i + 0*j + 0*k    qb=.   qnmarki   q
NB. qnmarkj      0 + 0*i + c*j + 0*k    qc=.   qnmarkj   q
NB. qnmarkk      0 + 0*i + 0*j + d*k    qd=.   qnmarkk   q
NB. qnmark1i     a + b*i + 0*j + 0*k    qab=.  qnmark1i  q
NB. qnmark1j     a + 0*i + c*j + 0*k    qac=.  qnmark1j  q
NB. qnmark1k     a + 0*i + 0*j + d*k    qad=.  qnmark1k  q
NB. qnmarkij     0 + b*i + c*j + 0*k    qbc=.  qnmarkij  q
NB. qnmarkik     0 + b*i + 0*j + d*k    qbd=.  qnmarkik  q
NB. qnmarkjk     0 + 0*i + c*j + d*k    qcd=.  qnmarkjk  q
NB. qnmark1ij    a + b*i + c*j + 0*k    q1ij=. qnmark1ij q
NB. qnmark1ik    a + b*i + 0*j + d*k    q1ik=. qnmark1ik q
NB. qnmark1jk    a + 0*i + c*j + d*k    q1jk=. qnmark1jk q
NB. qnmarkijk    0 + b*i + c*j + d*k    qijk=. qnmarkijk q

qnmark1=: (0 ,~    qn1) : [:
qnmarki=: (0 ,~ j.@qni) : [:
qnmarkj=: (0 ,     qnj) : [:
qnmarkk=: (0 ,  j.@qnk) : [:

qnmark1i=: 0&(1})          : [:
qnmarkjk=: 0&(0})          : [:
qnmark1j=:      9&o.       : [:
qnmarkik=: j.@(11&o.)      : [:
qnmark1k=: (qn1 ,  j.@qnk) : [:
qnmarkij=: (qnj ,~ j.@qni) : [:

qnmark1ij=: (1}~    qnj) : [:
qnmark1ik=: (1}~ j.@qnk) : [:
qnmark1jk=: (0}~    qn1) : [:
qnmarkijk=: (0}~ j.@qni) : [:

NB. ---------------------------------------------------------
NB. Conjugate component(s)
NB.
NB. Verb       Action                  Syntax
NB. qncon1     -a + b*i + c*j + d*k    qc=. qncon1  q
NB. qnconi      a - b*i + c*j + d*k    qc=. qnconi  q
NB. qnconj      a + b*i - c*j + d*k    qc=. qnconj  q
NB. qnconk      a + b*i + c*j - d*k    qc=. qnconk  q
NB. qnconij     a - b*i - c*j + d*k    qc=. qnconij q
NB. qnconjk     a + b*i - c*j - d*k    qc=. qnconjk q
NB. qnconik     a - b*i + c*j - d*k    qc=. qnconik q
NB. qnconv      a - b*i - c*j - d*k    qc=. qnconv  q
NB.
NB. References:
NB. [1] E. A. Karataev. Inner conjugation of quaternions.
NB.     Volzhskiy, 2002 (Е. А. Каратаев. Внутреннее
NB.     сопряжение кватернионов. Волжский, 2002).

NB. inner single conjugation
qncon1=: -@+`]    "0"1 : [:  NB. by 1
qnconi=:   +`]    "0"1 : [:  NB. by i
qnconj=:   ]`(-@+)"0"1 : [:  NB. by j
qnconk=:   ]`   + "0"1 : [:  NB. by k

NB. inner double conjugation
qnconik=: +                : [:  NB. by i,k
qnconjk=: ]`-"0"1          : [:  NB. by j,k
qnconij=: qnconik@:qnconjk : [:  NB. by i,j

NB. vector conjugation
qnconv=:  +`-"0"1          : [:  NB. by i,j,k

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Other operators

NB. ---------------------------------------------------------
NB. qnlen
NB.
NB. Description:
NB.   Length
NB.
NB. Syntax:
NB.   len=. qnlen q
NB. where
NB.   q   - a single quaternion or laminated quatertions
NB.   len - a scalar or (# q)-vector of non-negative numbers
NB.
NB. Assertions:
NB.   0 *./@:<: qnlen q

qnlen=: norms"1

NB. ---------------------------------------------------------
NB. qnsign
NB.
NB. Description:
NB.   Signum
NB.
NB. Syntax:
NB.   qsgn=. qnsign q
NB. where
NB.   q,qsgn - a single quaternion or laminated quatertions
NB.
NB. Assertions:
NB.   1 *./@:= qnlen qnsign q
NB.
NB. Notes:
NB. - throws NaN error for q∞

qnsign=: (      (;  (% +/&.:*:)&>/@(% :: (*@[)"0 L: 0)      >./@(dbsig@33^:(1 < _ +/@:= ]))) |@(,@:+.^:(JCMPX = 3!:0)))@(dbsig@33^:(0 0&(-:!.0)))`(_.j_. _.j_."_)@.(isnan@<)"1 : [:

NB. ---------------------------------------------------------
NB. qninv
NB.
NB. Description:
NB.   Inverse
NB.
NB. Syntax:
NB.   qi=. qninv q
NB. where
NB.   q,qi - a single quaternion or laminated quatertions
NB.
NB. Assertions (with appropriate comparison tolerance):
NB.   q -: qninv qninv q
NB.   1 0 *./@:(-:"1) q qnmul qninv q
NB.
NB. Notes:
NB. - throws NaN error for q∞

qninv=: (qnconv (; ((% +/@: *:)&>/@(% :: (*@[)"0 L: 0) % ]) >./@(dbsig@33^:(1 < _ +/@:= ]))) |@,@(+.^:(JCMPX = 3!:0)))@(dbsig@33^:(0 0&(-:!.0)))`(_.j_. _.j_."_)@.(isnan@<)"1 : [:

NB. ---------------------------------------------------------
NB. qnmul
NB.
NB. Description:
NB.   Product
NB.
NB. Syntax:
NB.   qp=. qx qnmul qy
NB. where
NB.   qx,qy,qp - a single quaternion or laminated quatertions

qnmul=: [: : (,.~@:(+./"1)@isnan`(,:&_.j_.)}@(mp"1 2 (,:"1 (-@:({:"1) ,. {."1)@:+)))

NB. ---------------------------------------------------------
NB. qndivl
NB. qndivr
NB.
NB. Description:
NB.   Divide
NB.
NB. Syntax:
NB.   qq=. qn qndivx qd
NB. where
NB.   qn,qd,qq - a single quaternion or laminated quatertions
NB.
NB. Notes:
NB. - throws NaN error for qn = q∞ or qd = q∞ or qd = 0

qndivl=: [: : (qnmul  qninv)  NB. via left quotient
qndivr=: [: : (qnmul~ qninv)  NB. via right quotient

NB. ---------------------------------------------------------
NB. qnf
NB.
NB. Description:
NB.   Adv. to quaternificate verb
NB.
NB. Formula [1]:
NB.   f(q) = Re(f(λ)) + sgn(q-a) * Im(f(λ))
NB. where
NB.   q = a + b*i + c*j + d*k ∊ ℍ
NB.   λ = a + |q-a|*i         ∊ ℂ
NB.   sgn(q) = q/|q|,     if q ≠ 0
NB.          = undefined, if q = 0
NB.
NB. Syntax:
NB.   o=. (f qnf) y
NB. where
NB.   f - monad to compute real or complex function value of
NB.       real or complex argument, is called as:
NB.         o=. f y
NB.   y - quaternion argument
NB.   o - quaternion function value of y
NB.
NB. Examples:
NB.    NB. e^1
NB.    ^ 1
NB. 2.71828
NB.    NB. quaternificated e^y of real y has b=c=d=0
NB.    ^ qnf 1 0
NB. 2.71828 0
NB.    ^ 1j2
NB. _1.1312j2.47173
NB.    NB. quaternificated e^y of complex y has c=d=0
NB.    ^ qnf 1j2 0
NB. _1.1312j2.47173 0
NB.    NB. quaternificated e^y of quaternion y
NB.    ^ qnf 1j2 3j4
NB. 1.69392j_0.78956 _1.18434j_1.57912
NB.
NB. Notes:
NB. - 0-rank approach
NB. - the following intermediate quaternion is used:
NB.     |q-a| + b*i + c*j + d*k
NB.   to aviod |q-a| value computing twice
NB.
NB. References:
NB. [1] L. G. Bairak. Integral Formula of Cauchy for
NB.     Quaternions. 2010 (Л. Г. Байрак. Интегральная формула
NB.     Коши для кватернионов. 2010).

qnf=: 1 : 'qn1_mt_ (u@(j. qn1_mt_) ((9 o. [) qn1_mt_ (* 11&o.)~) (% qn1_mt_)@]) qnlen_mt_@qnmarkijk_mt_ qn1_mt_ ]'
