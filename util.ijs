NB. Utilities
NB.
NB. isnan       Mark NaN values
NB. max         Max-of, 0 for empty list
NB. negneg      Conditional negate
NB. negpos      Conditional negate
NB. copysign    Copy sign
NB. sorim       Sum of real and imaginary parts' modules
NB. soris       Sum of real and imaginary parts' squares
NB. lcat        Concatenate logs
NB. nolog       Nilad to generate neutral for test actors
NB. tmonad      Conj. to make monad to test computational monad
NB. tdyad       Conj. to make monad to test computational dyad
NB. assert      Advanced version of the (assert.) control
NB. cut3        Split list by delimiter taken from its tail
NB. cut2        Split list by delimiter
NB. cut         Split list by delimiter
NB. cutl2       Split list by any delimiter
NB. cutl        Split list by any delimiter
NB.
NB. verifyutil  Verify util verbs
NB.
NB. Version: 0.13.0 2021-05-21
NB.
NB. Copyright 2010-2021 Igor Zhuravlov
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

NB. Format log string
fmtlog=: ;@:(40 17 17 17 17 _16&(({.{.@('d<n/a>'&(8!:2)))&.>))

NB. mark...
ispos0=:  _ =!.0 %  NB. ... +0  values in y
isneg0=: __ =!.0 %  NB. ... -0  values in y

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Miscellaneous

isnan=: 128!:5  NB. mark NaN values in y

max=: >./`0:@.(0 = #)  NB. max-of, 0 for empty list

negneg=:   ($@] $    (isneg0 +. 0&>)@{.`((,:  -)@{:)}@,:)`(_."0@])@.(+.&(isnan@<))  NB. if x<0 then - y  else  y  endif
negpos=:   ($@] $    (isneg0 +: 0&>)@{.`((,:  -)@{:)}@,:)`(_."0@])@.(+.&(isnan@<))  NB. if x≥0 then - y  else  y  endif
copysign=: ($@] $ =/&(isneg0 +. 0&>)   `((,:~ -)@{:)}@,:)`(_."0@])@.(+.&(isnan@<))  NB. if x<0 then -|y| else |y| endif

sorim=: | `(+/!.0"1@:| @:+.)@.(JCMPX = 3!:0)  NB. sum of real and imaginary parts' modules, |Re(y)| + |Im(y)|
soris=: *:`(+/!.0"1@:*:@:+.)@.(JCMPX = 3!:0)  NB. sum of real and imaginary parts' squares, Re(y)^2 + Im(y)^2

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. test suite utilities

lcat=: ,&.>&:(,&.>/"2^:(<:@#@$))

NB. ---------------------------------------------------------
NB. nolog
NB.
NB. Description:
NB.   Nilad to generate neutral for test actors
NB.
NB. Syntax:
NB.   emptylog=. nolog ''
NB. where
NB.   emptylog - an empty log which could be joined with
NB.              another log
NB.
NB. Assertions:
NB.   log -: log      ,&.> emptylog
NB.   log -: emptylog ,&.> log
NB. where
NB.   log - some another log
NB.
NB. Notes:
NB. - test actor may finish with the neutral result if input
NB.   is not applicable (say, non-square matrix for a method
NB.   requiring square input)

nolog=: 1 5 # EMPTY ; $@0

NB. ---------------------------------------------------------
NB. tmonad
NB. tdyad
NB.
NB. Description:
NB.   Conj. to make monad to test computational verb
NB.
NB. Syntax:
NB.   'omsent rcond ferr berr time space'=. (imsent tmonad        vgety`vgeto`vrcond`vferr`vberr) y
NB.   'odsent rcond ferr berr time space'=. (idsent tdyad   vgetx`vgety`vgeto`vrcond`vferr`vberr) y
NB. where
NB.   imsent - string, J sentence for monadic execution; is
NB.            called as:
NB.              ret=.      sentence argy
NB.   idsent - string, J sentence for dyadic execution; is
NB.            called as:
NB.              ret=. argx sentence argy
NB.   vgetx  - monad to extract left argument for vd; is
NB.            called as:
NB.              argx=. vgetx y
NB.   vgety  - monad to extract right argument for vm or vd;
NB.            is called as:
NB.              argy=. vgety y
NB.   vgeto  - monad to extract output from ret;
NB.            is called as:
NB.              out=. vgeto ret
NB.   vrcond - monad to find rcond; is called as:
NB.              rcond=. vrcond y
NB.   vferr  - dyad to find ferr; is called as:
NB.              ferr=. y vferr out
NB.   vberr  - dyad to find berr; is called as:
NB.              berr=. y vberr out
NB.   y      - some input for monad made from conj.
NB.   omsent -: ,: imsent
NB.   odsent -: ,: idsent
NB.   rcond  ≥ 0, the estimated reciprocal of the condition
NB.            number of the input matrix; +∞ if matrix is
NB.            singular; NaN if matrix is non-square
NB.   ferr   ≥ 0 or NaN, the relative forward error
NB.   berr   ≥ 0 or NaN, the relative backward error
NB.   time   ≥ 0, estimation [1] of sentence execution time
NB.   space  ≥ 0, the number of bytes used to execute the
NB.            sentence
NB.   argx   - some left argument for sentence
NB.   argy   - some right argument for sentence
NB.   ret    - some result of sentence execution
NB.   out    - rectified ret, i.e. filtered output
NB.
NB. Application:
NB. - to test geqrf:
NB.     NB. to estimate rcond in 1-norm
NB.     vrcond=. (_."_)`gecon1@.(=/@$)
NB.     NB. to calc. berr, assuming:
NB.     NB.   berr := ||A - realA||_1 / (FP_EPS * ||A||_1 * m)
NB.     vberr=. ((- %&norm1 [) % FP_EPS * (norm1 * #)@[) unmqr
NB.     NB. do the job
NB.     'sent rcond ferr berr time space'=. ('geqrf' tmonad ]`]`vrcond`(_."_)`vberr) A
NB. - to test getrs:
NB.     NB. to estimate rcond in ∞-norm
NB.     vrcond=. (_."_)`geconi@.(=/@$)@(0&{::)
NB.     NB. to calc. ferr, assuming:
NB.     NB.   ferr := ||x - realx||_inf / ||realx||_inf
NB.     vferr=. ((- %&normi [) 1&{::)~
NB.     NB. to calc. componentwise berr [LUG 75], assuming:
NB.     NB.   berr := max_i(|b - A * realx|_i / (|A| * |realx| + |b|)_i)
NB.     vberr=. (mp&>/@[ |@- (0 {:: [) mp ]) >./@% (((0 {:: [) mp&| ]) + |@mp&>/@[)
NB.     NB. do the job
NB.     'sent rcond ferr berr time space'=. ('getrs' tdyad (0&{::)`(mp&>/)`]`vrcond`vferr`vberr) (A;x)
NB.
NB. References:
NB. [1] Magne Haveraaen, Hogne Hundvebakke. Some Statistical
NB.     Performance Estimation Techniques for Dynamic
NB.     Machines. Appeared in Weihai Yu & al. (eds.): Norsk
NB.     Informatikk-konferanse 2001, Tapir, Trondheim Norway
NB.     2001, pp. 176-185.
NB.     https://www.ii.uib.no/saga/papers/perfor-5d.pdf
NB.
NB. Notes:
NB. 1) recommended observations count to provide standard
NB.    deviation <= 1% for CPU run-time estimator minimum on
NB.    systems with load ≤ 80 is equal to 5 [1]
NB. 2) side effect: the result is sent to the console

tmonad=: 2 : 0
  '`vgety vgeto vrcond vferr vberr'=. n
  try. rcond=. vrcond y catch. rcond=. _ end.
  try.
    ybak=. memu argy=. vgety y
    try.
      't s'=. , (5 1 # i. 2) <./`]/. (5 1 # timex`(7!:2))`:0 'ret=. ' , m , ' argy'
      if. -. argy -: ybak do. m=. m , ' NB. error: y changed' end.
      try.
        out=. vgeto ret
        try. ferr=. y vferr out catch. ferr=. _. end.
        try. berr=. y vberr out catch. berr=. _. end.
      catch.
        'ferr berr'=. 2 # _.
      end.
    catch.
      dbsig 3  NB. jump to upper catch block
    end.
  catch.
    'ferr berr t s'=. 4 # _.
  end.
  erase 'argy ybak'
  logline=. (,: m) ; rcond ; ferr ; berr ; t ; s
  (fmtlog_mt_ logline) (1!:2) 2
  wd^:IFQT 'msgs'
  logline
)

tdyad=: 2 : 0
  '`vgetx vgety vgeto vrcond vferr vberr'=. n
  try. rcond=. vrcond y catch. rcond=. _ end.
  try.
    xbak=. memu argx=. vgetx y
    ybak=. memu argy=. vgety y
    try.
      't s'=. , (5 1 # i. 2) <./`]/. (5 1 # timex`(7!:2))`:0 'ret=. argx ' , m , ' argy'
      select. #. (argx -: xbak) , argy -: ybak
        case. 2 do. m=. m , ' NB. error: y was changed'
        case. 1 do. m=. m , ' NB. error: x was changed'
        case. 0 do. m=. m , ' NB. error: x and y were changed'
      end.
      try.
        out=. vgeto ret
        try. ferr=. y vferr out catch. ferr=. _. end.
        try. berr=. y vberr out catch. berr=. _. end.
      catch.
        'ferr berr'=. 2 # _.
      end.
    catch.
      dbsig 3  NB. jump to upper catch block
    end.
  catch.
    'ferr berr t s'=. 4 # _.
  end.
  erase 'argy ybak argx xbak'
  logline=. (,: m) ; rcond ; ferr ; berr ; t ; s
  (fmtlog_mt_ logline) (1!:2) 2
  wd^:IFQT 'msgs'
  logline
)

NB. end of test suite utilities
NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. verification suite utilities

reportv=: 1 : '] [ echo@(m , '': assertions probed: '' , ":@{. , '', failed: '' , ":@{:)'

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. flt staff

NB. ---------------------------------------------------------
NB. assert
NB.
NB. Description:
NB.   Advanced version of the (assert.) control
NB.
NB. Syntax:
NB.   trash=. [msg] assert chk
NB. where
NB.   msg - string, optional, will be shown if assertion is
NB.         failed
NB.   chk - numeric vector
NB.
NB. Notes:
NB. - fixes system's (assert) to match (assert.) control
NB. - is equipped with error message
NB.
NB. References:
NB. [1] Igor Zhuravlov. [Jprogramming] assert verb from
NB.     stdlib mismatches assert. control
NB.     2019-12-30 00:43:46 UTC.
NB.     http://www.jsoftware.com/pipermail/programming/2019-December/054693.html

assert=: 0 0 $ dbsig^:((1 +./@:~: ])`(12"_))^:(9!:34@'')

NB. ---------------------------------------------------------
NB. Notes:
NB. - the following definitions are assumed in the examples
NB.   below:
NB.     delimiters=. LF , ' '
NB.     list=. 'foo bar  baz' , LF , 'qux' , LF2 , 'quux' , LF , ' corge ' , LF , 'flob'

NB. ---------------------------------------------------------
NB. cut3
NB.
NB. Description:
NB.   Split list by delimiter taken from its tail
NB.
NB. Syntax:
NB.   splitted_list=. cut3 list
NB.
NB. Examples:
NB.    cut3 list
NB. +----+----+------------------------+
NB. |foo |ar  |az qux  quux  corge  flo|
NB. +----+----+------------------------+

cut3=: <;._2

NB. ---------------------------------------------------------
NB. cut2
NB.
NB. Description:
NB.   Split list by delimiter
NB.
NB. Syntax:
NB.   splitted_list=. [delimiter] cut2 list
NB.
NB. Notes:
NB. - default delimiter is ' ' character
NB. - doesn't drop repeating delimiters
NB. - monad (cut2) is an inverse of (' '&joinstring)
NB. - dyad (cut2) is an inverse of (joinstring)
NB.
NB. Examples:
NB.    cut2 list
NB. +---+---++--------------+-----+-----+
NB. |foo|bar||baz qux  quux |corge| flob|
NB. +---+---++--------------+-----+-----+
NB.    LF cut2 list
NB. +------------+---++----+-------+----+
NB. |foo bar  baz|qux||quux| corge |flob|
NB. +------------+---++----+-------+----+

cut2=: ' '&$: : (cut3@,~)

NB. ---------------------------------------------------------
NB. cut
NB.
NB. Description:
NB.   Split list by delimiter
NB.
NB. Syntax:
NB.   splitted_list=. [delimiter] cut list
NB.
NB. Notes:
NB. - default delimiter is ' ' character
NB. - drops repeating delimiters
NB. - identic to (cut) verb from the Standard Library
NB.
NB. Examples:
NB.    cut list
NB. +---+---+--------------+-----+-----+
NB. |foo|bar|baz qux  quux |corge| flob|
NB. +---+---+--------------+-----+-----+
NB.    LF cut list
NB. +------------+---+----+-------+----+
NB. |foo bar  baz|qux|quux| corge |flob|
NB. +------------+---+----+-------+----+

cut=: -.&a:@cut2

NB. ---------------------------------------------------------
NB. cutl2
NB.
NB. Description:
NB.   Split list by any delimiter
NB.
NB. Syntax:
NB.   splitted_list=. delimiters cutl2 list
NB.
NB. Notes:
NB. - like (cut2) verb, but delimiter can be any element
NB.   from (delimiters) argument
NB. - like (cutl) verb, but doesn't drop repeating delimiters
NB.
NB. Examples:
NB.    delimiters cutl2 list
NB. +---+---++---+---++----++-----++----+
NB. |foo|bar||baz|qux||quux||corge||flob|
NB. +---+---++---+---++----++-----++----+

cutl2=: ((, {.) (e. cut3 [) ])~

NB. ---------------------------------------------------------
NB. cutl
NB.
NB. Description:
NB.   Split list by any delimiter
NB.
NB. Syntax:
NB.   splitted_list=. delimiters cutl list
NB.
NB. Notes:
NB. - like (cut) verb, but delimiter can be any element
NB.   from (delimiters) argument
NB. - drops repeating delimiters
NB.
NB. Examples:
NB.    delimiters cutl list
NB. +---+---+---+---+----+-----+----+
NB. |foo|bar|baz|qux|quux|corge|flob|
NB. +---+---+---+---+----+-----+----+

cutl=: -.&a:@cutl2

NB. end of flt staff
NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NB. =========================================================
NB. Verification suite

NB. ---------------------------------------------------------
NB. verifyutil
NB.
NB. Description:
NB.   Nilad to verify util actors, output result to console
NB.   and return it
NB.
NB. Syntax:
NB.   'probed failed'=. verifyutil ''
NB. where
NB.   probed ≥ 0, assertions probed counter
NB.   failed ≥ 0, assertions failed counter

verifyutil=: 3 : 0
  delimiters=. LF , ' '
  string=. 'foo bar  baz' , LF , 'qux' , LF2 , 'quux' , LF , ' corge ' , LF , 'flob'

  NB. verify itself
  res=.       0 0 1 0 0 -: isnan  _1 _0.0 _. 0.0 1
  res=. res , 0 0 0 1 0 -: ispos0 _1 _0.0 _. 0.0 1
  res=. res , 0 1 0 0 0 -: isneg0 _1 _0.0 _. 0.0 1

  NB. max
  res=. res ,  0 -:  max ''
  res=. res , __ -:  max _.   _.
  res=. res ,  _ -:  max _.    _ __
  res=. res ,  1 -:  max  0    1 _.
  res=. res ,  1 -:  max  0    1 _. __
  res=. res , ispos0 max  0.0 _1
  res=. res , isneg0 max _0.0 _1

  NB. negneg

  res=. res ,        isnan _.   negneg _.
  res=. res ,        isnan _.   negneg __
  res=. res ,        isnan _.   negneg _1
  res=. res ,        isnan _.   negneg _0.0
  res=. res ,        isnan _.   negneg  0
  res=. res ,        isnan _.   negneg  1
  res=. res ,        isnan _.   negneg  _
  res=. res , 1 1 -: isnan _.   negneg _. _.
  res=. res , 1 1 -: isnan _.   negneg _.  2
  res=. res , 1 1 -: isnan _.   negneg __  2
  res=. res ,        isnan __   negneg _.
  res=. res ,        isnan _1   negneg _.
  res=. res ,        isnan _0.0 negneg _.
  res=. res ,        isnan  0   negneg _.
  res=. res ,        isnan  1   negneg _.
  res=. res ,        isnan  _   negneg _.
  res=. res , 1 1 -: isnan __   negneg _. _.
  res=. res , 1 1 -: isnan _1   negneg _. _.
  res=. res , 1 1 -: isnan _0.0 negneg _. _.
  res=. res , 1 1 -: isnan  0   negneg _. _.
  res=. res , 1 1 -: isnan  1   negneg _. _.
  res=. res , 1 1 -: isnan  _   negneg _. _.
  res=. res , 1 1 -: isnan __   negneg _.  2
  res=. res , 1 1 -: isnan _1   negneg _.  2
  res=. res , 1 1 -: isnan _0.0 negneg _.  2
  res=. res , 1 1 -: isnan  0   negneg _.  2
  res=. res , 1 1 -: isnan  1   negneg _.  2
  res=. res , 1 1 -: isnan  _   negneg _.  2

  res=. res , '' -: __   negneg ''
  res=. res , '' -: _1   negneg ''
  res=. res , '' -: _0.0 negneg ''
  res=. res , '' -:  0   negneg ''
  res=. res , '' -:  1   negneg ''
  res=. res , '' -:  _   negneg ''

  res=. res ,  _ -: __   negneg __
  res=. res ,  _ -: _1   negneg __
  res=. res ,  _ -: _0.0 negneg __
  res=. res , __ -:  0   negneg __
  res=. res , __ -:  1   negneg __
  res=. res , __ -:  _   negneg __

  res=. res ,  1 -: __   negneg _1
  res=. res ,  1 -: _1   negneg _1
  res=. res ,  1 -: _0.0 negneg _1
  res=. res , _1 -:  0   negneg _1
  res=. res , _1 -:  1   negneg _1
  res=. res , _1 -:  _   negneg _1

  res=. res , ispos0 __   negneg _0.0
  res=. res , ispos0 _1   negneg _0.0
  res=. res , ispos0 _0.0 negneg _0.0
  res=. res , isneg0  0   negneg _0.0
  res=. res , isneg0  1   negneg _0.0
  res=. res , isneg0  _   negneg _0.0

  res=. res , isneg0 __   negneg 0.0
  res=. res , isneg0 _1   negneg 0.0
  res=. res , isneg0 _0.0 negneg 0.0
  res=. res , ispos0  0   negneg 0.0
  res=. res , ispos0  1   negneg 0.0
  res=. res , ispos0  _   negneg 0.0

  res=. res , _1 -: __   negneg 1
  res=. res , _1 -: _1   negneg 1
  res=. res , _1 -: _0.0 negneg 1
  res=. res ,  1 -:  0   negneg 1
  res=. res ,  1 -:  1   negneg 1
  res=. res ,  1 -:  _   negneg 1

  res=. res , __ -: __   negneg _
  res=. res , __ -: _1   negneg _
  res=. res , __ -: _0.0 negneg _
  res=. res ,  _ -:  0   negneg _
  res=. res ,  _ -:  1   negneg _
  res=. res ,  _ -:  _   negneg _

  res=. res , ( _  1 0 0 _1 __&-: *. (0 0 0 1 0 0 -: isneg0) *. 0 0 1 0 0 0 -: ispos0) _1 negneg __ _1 _0.0 0 1 _
  res=. res , (__ _1 0 0  1  _&-: *. (0 0 1 0 0 0 -: isneg0) *. 0 0 0 1 0 0 -: ispos0)  1 negneg __ _1 _0.0 0 1 _

  NB. negpos

  res=. res ,        isnan _.   negpos _.
  res=. res ,        isnan _.   negpos __
  res=. res ,        isnan _.   negpos _1
  res=. res ,        isnan _.   negpos _0.0
  res=. res ,        isnan _.   negpos  0
  res=. res ,        isnan _.   negpos  1
  res=. res ,        isnan _.   negpos  _
  res=. res , 1 1 -: isnan _.   negpos _. _.
  res=. res , 1 1 -: isnan _.   negpos _.  2
  res=. res , 1 1 -: isnan _.   negpos __  2
  res=. res ,        isnan __   negpos _.
  res=. res ,        isnan _1   negpos _.
  res=. res ,        isnan _0.0 negpos _.
  res=. res ,        isnan  0   negpos _.
  res=. res ,        isnan  1   negpos _.
  res=. res ,        isnan  _   negpos _.
  res=. res , 1 1 -: isnan __   negpos _. _.
  res=. res , 1 1 -: isnan _1   negpos _. _.
  res=. res , 1 1 -: isnan _0.0 negpos _. _.
  res=. res , 1 1 -: isnan  0   negpos _. _.
  res=. res , 1 1 -: isnan  1   negpos _. _.
  res=. res , 1 1 -: isnan  _   negpos _. _.
  res=. res , 1 1 -: isnan __   negpos _.  2
  res=. res , 1 1 -: isnan _1   negpos _.  2
  res=. res , 1 1 -: isnan _0.0 negpos _.  2
  res=. res , 1 1 -: isnan  0   negpos _.  2
  res=. res , 1 1 -: isnan  1   negpos _.  2
  res=. res , 1 1 -: isnan  _   negpos _.  2

  res=. res , '' -: __   negpos ''
  res=. res , '' -: _1   negpos ''
  res=. res , '' -: _0.0 negpos ''
  res=. res , '' -:  0   negpos ''
  res=. res , '' -:  1   negpos ''
  res=. res , '' -:  _   negpos ''

  res=. res , __ -: __   negpos __
  res=. res , __ -: _1   negpos __
  res=. res , __ -: _0.0 negpos __
  res=. res ,  _ -:  0   negpos __
  res=. res ,  _ -:  1   negpos __
  res=. res ,  _ -:  _   negpos __

  res=. res , _1 -: __   negpos _1
  res=. res , _1 -: _1   negpos _1
  res=. res , _1 -: _0.0 negpos _1
  res=. res ,  1 -:  0   negpos _1
  res=. res ,  1 -:  1   negpos _1
  res=. res ,  1 -:  _   negpos _1

  res=. res , isneg0 __   negpos _0.0
  res=. res , isneg0 _1   negpos _0.0
  res=. res , isneg0 _0.0 negpos _0.0
  res=. res , ispos0  0   negpos _0.0
  res=. res , ispos0  1   negpos _0.0
  res=. res , ispos0  _   negpos _0.0

  res=. res , ispos0 __   negpos 0.0
  res=. res , ispos0 _1   negpos 0.0
  res=. res , ispos0 _0.0 negpos 0.0
  res=. res , isneg0  0   negpos 0.0
  res=. res , isneg0  1   negpos 0.0
  res=. res , isneg0  _   negpos 0.0

  res=. res ,  1 -: __   negpos 1
  res=. res ,  1 -: _1   negpos 1
  res=. res ,  1 -: _0.0 negpos 1
  res=. res , _1 -:  0   negpos 1
  res=. res , _1 -:  1   negpos 1
  res=. res , _1 -:  _   negpos 1

  res=. res ,  _ -: __   negpos _
  res=. res ,  _ -: _1   negpos _
  res=. res ,  _ -: _0.0 negpos _
  res=. res , __ -:  0   negpos _
  res=. res , __ -:  1   negpos _
  res=. res , __ -:  _   negpos _

  res=. res , (__ _1 0 0  1  _&-: *. (0 0 1 0 0 0 -: isneg0) *. 0 0 0 1 0 0 -: ispos0) _1 negpos __ _1 _0.0 0 1 _
  res=. res , ( _  1 0 0 _1 __&-: *. (0 0 0 1 0 0 -: isneg0) *. 0 0 1 0 0 0 -: ispos0)  1 negpos __ _1 _0.0 0 1 _

  NB. copysign

  res=. res ,        isnan _.   copysign _.
  res=. res ,        isnan _.   copysign __
  res=. res ,        isnan _.   copysign _1
  res=. res ,        isnan _.   copysign _0.0
  res=. res ,        isnan _.   copysign  0
  res=. res ,        isnan _.   copysign  1
  res=. res ,        isnan _.   copysign  _
  res=. res , 1 1 -: isnan _.   copysign _. _.
  res=. res , 1 1 -: isnan _.   copysign _.  2
  res=. res , 1 1 -: isnan _.   copysign __  2
  res=. res ,        isnan __   copysign _.
  res=. res ,        isnan _1   copysign _.
  res=. res ,        isnan _0.0 copysign _.
  res=. res ,        isnan  0   copysign _.
  res=. res ,        isnan  1   copysign _.
  res=. res ,        isnan  _   copysign _.
  res=. res , 1 1 -: isnan __   copysign _. _.
  res=. res , 1 1 -: isnan _1   copysign _. _.
  res=. res , 1 1 -: isnan _0.0 copysign _. _.
  res=. res , 1 1 -: isnan  0   copysign _. _.
  res=. res , 1 1 -: isnan  1   copysign _. _.
  res=. res , 1 1 -: isnan  _   copysign _. _.
  res=. res , 1 1 -: isnan __   copysign _.  2
  res=. res , 1 1 -: isnan _1   copysign _.  2
  res=. res , 1 1 -: isnan _0.0 copysign _.  2
  res=. res , 1 1 -: isnan  0   copysign _.  2
  res=. res , 1 1 -: isnan  1   copysign _.  2
  res=. res , 1 1 -: isnan  _   copysign _.  2

  res=. res , '' -: __   copysign ''
  res=. res , '' -: _1   copysign ''
  res=. res , '' -: _0.0 copysign ''
  res=. res , '' -:  0   copysign ''
  res=. res , '' -:  1   copysign ''
  res=. res , '' -:  _   copysign ''

  res=. res , __ -: __   copysign __
  res=. res , __ -: _1   copysign __
  res=. res , __ -: _0.0 copysign __
  res=. res ,  _ -:  0   copysign __
  res=. res ,  _ -:  1   copysign __
  res=. res ,  _ -:  _   copysign __

  res=. res , _1 -: __   copysign _1
  res=. res , _1 -: _1   copysign _1
  res=. res , _1 -: _0.0 copysign _1
  res=. res ,  1 -:  0   copysign _1
  res=. res ,  1 -:  1   copysign _1
  res=. res ,  1 -:  _   copysign _1

  res=. res , isneg0 __   copysign _0.0
  res=. res , isneg0 _1   copysign _0.0
  res=. res , isneg0 _0.0 copysign _0.0
  res=. res , ispos0  0   copysign _0.0
  res=. res , ispos0  1   copysign _0.0
  res=. res , ispos0  _   copysign _0.0

  res=. res , isneg0 __   copysign  0.0
  res=. res , isneg0 _1   copysign  0.0
  res=. res , isneg0 _0.0 copysign  0.0
  res=. res , ispos0  0   copysign  0.0
  res=. res , ispos0  1   copysign  0.0
  res=. res , ispos0  _   copysign  0.0

  res=. res , _1 -: __   copysign  1
  res=. res , _1 -: _1   copysign  1
  res=. res , _1 -: _0.0 copysign  1
  res=. res ,  1 -:  0   copysign  1
  res=. res ,  1 -:  1   copysign  1
  res=. res ,  1 -:  _   copysign  1

  res=. res , __ -: __   copysign  _
  res=. res , __ -: _1   copysign  _
  res=. res , __ -: _0.0 copysign  _
  res=. res ,  _ -:  0   copysign  _
  res=. res ,  _ -:  1   copysign  _
  res=. res ,  _ -:  _   copysign  _

  res=. res , (__ _1 0 0 _1 __&-: *. (0 0 1 1 0 0 -: isneg0) *. 0 0 0 0 0 0 -: ispos0) _1 copysign __ _1 _0.0 0 1 _
  res=. res , ( _  1 0 0  1  _&-: *. (0 0 0 0 0 0 -: isneg0) *. 0 0 1 1 0 0 -: ispos0)  1 copysign __ _1 _0.0 0 1 _

  NB. sorim
  res=. res , '' -: sorim ''
  res=. res , 0 2 3 -: sorim 0 2 _3
  res=. res , (_. _ _ 0 2 2 7 7 7 7&-: *. 1 0 0 0 0 0 0 0 0 0 -: isnan) sorim _. __ _ 0 _2 2 _3j_4 _3j4 3j_4 3j4

  NB. soris
  res=. res , '' -: soris ''
  res=. res , 0 4 9 -: soris 0 2 _3
  res=. res , (_. _ _ 0 4 4 25 25 25 25&-: *. 1 0 0 0 0 0 0 0 0 0 -: isnan) soris _. __ _ 0 _2 2 _3j_4 _3j4 3j_4 3j4

  NB. cut3
  res=. res , '' -: cut3 ''
  res=. res , ('foo ' ; 'ar  ' ; 'az' , LF , 'qux' , LF2 , 'quux' , LF , ' corge ' , LF , 'flo') -: cut3 string

  NB. cut2
  res=. res , (, a:) -: cut2 ''
  res=. res , ('foo' ; 'bar' ; '' ; ('baz' , LF , 'qux' , LF2 , 'quux' , LF) ; 'corge' ; LF , 'flob') -: cut2 string
  res=. res , ('foo bar  baz' ; 'qux' ; '' ; 'quux' ; ' corge ' ; 'flob') -: LF cut2 string

  NB. cut
  res=. res , '' -: cut ''
  res=. res , ('foo' ; 'bar' ; ('baz' , LF , 'qux' , LF2 , 'quux' , LF) ; 'corge' ; LF , 'flob') -: cut string
  res=. res , ('foo bar  baz' ; 'qux' ; 'quux' ; ' corge ' ; 'flob') -: LF cut string

  NB. cutl2
  res=. res , (, a:) -: delimiters cutl2 ''
  res=. res , ('foo' ; 'bar' ; '' ; 'baz' ; 'qux' ; '' ; 'quux' ; '' ; 'corge' ; '' ; 'flob') -: delimiters cutl2 string

  NB. cutl
  res=. res , '' -: delimiters cutl ''
  res=. res , ('foo' ; 'bar' ; 'baz' ; 'qux' ; 'quux' ; 'corge' ; 'flob') -: delimiters cutl string

  'util' reportv (# ([ , -) +/) res
)
