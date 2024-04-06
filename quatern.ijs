NB. Quaternions
NB.
NB. qnxx         Get/set component(s)
NB. qnmarkxxx    Mark component(s)
NB. qnconxx      Conjugate component(s)
NB. qnlen        Length
NB. qnsign       Signum
NB. qninv        Inverse
NB. qnmul        Multiply
NB. qndivl       Divide (left quotient)
NB. qndivr       Divide (right quotient)
NB. qnf          Adv. to quaternificate verb
NB.
NB. testqn1      Test verbs receiving the only quaternion
NB. testqn2      Test verbs receiving a couple of quaternions
NB. testqnf      Test qnf adverb
NB. testquatern  Adv. to make verb to test quaternion
NB.              algorithms by 2-vectors of generator given
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

NB. =========================================================
NB. Concepts
NB.
NB. Notation:
NB.   ℝ  ≡ [-∞,+∞]
NB.   a,b,c,d ∈ ℝ       Basis elements multiplication table:
NB.   ℂi ≡ ℝ  + ℝ *i          i   j   k
NB.   ℂj ≡ ℝ  + ℝ *j      i  _1   k  -j
NB.   ℂk ≡ ℝ  + ℝ *k      j  -k  _1   i
NB.   ℍ  ≡ ℂi + ℂi*j      k   j  -i  _1
NB.   q  = a + b*i + c*j + d*k ∈ ℍ
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
NB.     x = a + b*i ∈ ℂi          x -: a j. b
NB.     y = c + d*i ∈ ℂi          y -: c j. d
NB.     z = a + c*j ∈ ℂj          z -: a j. c
NB.     w = a + d*k ∈ ℂk          w -: a j. d
NB.     u = b + c*k ∈ ℂk          u -: b j. c
NB.     v = b + d*j ∈ ℂj          v -: b j. d
NB.     s = c + b*k ∈ ℂk          s -: c j. b
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
NB. Configuration

coclass 'mt'

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

qnsign=: (      (;  (% +/!.0&.:*:)&>/@(% :: (*@[)"0 L: 0)      >./@(dbsig@33^:(1 < _ +/@:= ]))) |@(,@:+.^:(JCMPX = 3!:0)))@(dbsig@33^:(0 0&(-:!.0)))`(2 # nan)@.(isnan@<)"1 : [:

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

qninv=: (qnconv (; ((% +/!.0@: *:)&>/@(% :: (*@[)"0 L: 0) % ]) >./@(dbsig@33^:(1 < _ +/@:= ]))) |@,@(+.^:(JCMPX = 3!:0)))@(dbsig@33^:(0 0&(-:!.0)))`(2 # nan)@.(isnan@<)"1 : [:

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

qnmul=: [: : (,.~@:(+./"1)@isnan`(,: nan)}@(mp"1 2 (,:"1 (-@:({:"1) ,. {."1)@:+)))

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
NB.   q = a + b*i + c*j + d*k ∈ ℍ
NB.   λ = a + |q-a|*i         ∈ ℂ
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

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testqn1
NB.
NB. Description:
NB.   Test verbs receiving the only quaternion:
NB.   - qnxx (math/mt addon)
NB.   - qnmarkxxx (math/mt addon)
NB.   - qnconxx (math/mt addon)
NB.   - qnlen (math/mt addon)
NB.   - qnsign (math/mt addon)
NB.   - qninv (math/mt addon)
NB.   by 2-vectors
NB.
NB. Syntax:
NB.   log=. testqn1 Qn
NB. where
NB.   Qn  - m×4-matrix of m 3-tuples (q,a,b)
NB.   log - 6-vector of boxes, test log, see test.ijs
NB.   q   - 2-vector, quaternion
NB.   a,b - any scalars

testqn1=: 3 : 0
  'Q A B'=. 2 ({."1 ; ({."1 ; {:"1)@(9 o. }."1)) y

  log=.          ('qn1"1'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qni"1'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnj"1'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnk"1'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qn1i"1' tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qn1j"1' tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qn1k"1' tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnij"1' tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnik"1' tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnjk"1' tmonad (]`]`nan`0:`0:)) Q

  log=. log lcat ('qn1"0 1'  tdyad (1&{::`(0&{::)`]`nan`0:`0:)) Q ; A
  log=. log lcat ('qni"0 1'  tdyad (1&{::`(0&{::)`]`nan`0:`0:)) Q ; A
  log=. log lcat ('qnj"0 1'  tdyad (1&{::`(0&{::)`]`nan`0:`0:)) Q ; A
  log=. log lcat ('qnk"0 1'  tdyad (1&{::`(0&{::)`]`nan`0:`0:)) Q ; A
  log=. log lcat ('qn1i"0 1' tdyad (1&{::`(0&{::)`]`nan`0:`0:)) Q ; A j. B
  log=. log lcat ('qn1j"0 1' tdyad (1&{::`(0&{::)`]`nan`0:`0:)) Q ; A j. B
  log=. log lcat ('qn1k"0 1' tdyad (1&{::`(0&{::)`]`nan`0:`0:)) Q ; A j. B
  log=. log lcat ('qnij"0 1' tdyad (1&{::`(0&{::)`]`nan`0:`0:)) Q ; A j. B
  log=. log lcat ('qnik"0 1' tdyad (1&{::`(0&{::)`]`nan`0:`0:)) Q ; A j. B
  log=. log lcat ('qnjk"0 1' tdyad (1&{::`(0&{::)`]`nan`0:`0:)) Q ; A j. B

  log=. log lcat ('qnmark1"1'   tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmarki"1'   tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmarkj"1'   tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmarkk"1'   tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmark1i"1'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmark1j"1'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmark1k"1'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmarkij"1'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmarkik"1'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmarkjk"1'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmark1ij"1' tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmark1ik"1' tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmark1jk"1' tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnmarkijk"1' tmonad (]`]`nan`0:`0:)) Q

  log=. log lcat ('qncon1'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnconi'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnconj'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnconk'  tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnconij' tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnconjk' tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnconik' tmonad (]`]`nan`0:`0:)) Q
  log=. log lcat ('qnconv'  tmonad (]`]`nan`0:`0:)) Q

  log=. log lcat ('qnlen' tmonad (]`]`nan`nan`nan)) Q

  log=. log lcat ('qnsign' tmonad (]`]`nan`nan`nan)) Q

  log=. log lcat ('qninv' tmonad (]`]`nan`nan`nan)) Q
)

NB. ---------------------------------------------------------
NB. testqn2
NB.
NB. Description:
NB.   Test verbs receiving a couple of quaternions:
NB.   - qnmul (math/mt addon)
NB.   - qndinl qndivr (math/mt addon)
NB.   by 2-vectors
NB.
NB. Syntax:
NB.   log=. testqn2 Q1Q2
NB. where
NB.   Q1Q2  - m×4-matrix of m quaternion pairs (q1,q2)
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.   q1,q2 - 2-vector, any quaternions

testqn2=: 3 : 0
  log=.          ('qnmul'  tdyad (2&({."1)`(2&(}."1))`]`nan`nan`nan)) y
  log=. log lcat ('qndivl' tdyad (2&({."1)`(2&(}."1))`]`nan`nan`nan)) y
  log=. log lcat ('qndivr' tdyad (2&({."1)`(2&(}."1))`]`nan`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testqnf
NB.
NB. Description:
NB.   Test quaternificate adverb:
NB.   - qnf (math/mt addon)
NB.   by 2-vectors
NB.
NB. Syntax:
NB.   log=. testqnf Q
NB. where
NB.   Q   - m×2-matrix of m quaternions
NB.   log - 6-vector of boxes, test log, see test.ijs

testqnf=: 3 : 0
  log=.          ('*: qnf"1' tmonad (]`]`nan`nan`nan)) y
  log=. log lcat ('%: qnf"1' tmonad (]`]`nan`nan`nan)) y
  log=. log lcat ('^  qnf"1' tmonad (]`]`nan`nan`nan)) y
  log=. log lcat ('^. qnf"1' tmonad (]`]`nan`nan`nan)) y
)

NB. ---------------------------------------------------------
NB. testquatern
NB.
NB. Description:
NB.   Adv. to make verb to test quaternion algorithms by
NB.   matrix of generator and shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testquatern) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Application:
NB. - test by random rectangular real matrix with elements
NB.   distributed uniformly with support (0,1):
NB.     log=. ?@$&0 testquatern_mt_ 200 150
NB. - test by random square real matrix with elements with
NB.   limited value's amplitude:
NB.     log=. _1 1 0 4 _6 4&gemat_mt_ testquatern_mt_ 200 200
NB. - test by random rectangular complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testquatern_mt_ 150 200

testquatern=: 1 : 'testqnf_mt_@u@(2&(1})) ,&.>~ (testqn2_mt_ ,&.>~ testqn1_mt_)@u@(4&(1}))'
