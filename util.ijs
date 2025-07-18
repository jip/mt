NB. Utilities
NB.
NB. ispos0    Mark +0 values
NB. isneg0    Mark -0 values
NB. isnan     Mark NaN values
NB. nan       Produce NaN of input datatype
NB. max       Max-of, 0 for empty list
NB. negneg    Conditional negate
NB. negpos    Conditional negate
NB. copysign  Copy sign
NB. sorim     Sum of real and imaginary parts' modules
NB. soris     Sum of real and imaginary parts' squares
NB. lcat      Concatenate logs
NB. nolog     Nilad to generate neutral for test actors
NB. tmonad    Conj. to make monad to test computational monad
NB. tdyad     Conj. to make monad to test computational dyad
NB. assert    Advanced version of the (assert.) control
NB. mrg       Adv. to merge arguments interleaving their
NB.           items as determined by a list
NB. cut3      Split list by delimiter taken from its tail
NB. cut2      Split list by delimiter
NB. cut       Split list by delimiter
NB. cutl2     Split list by any delimiter
NB. cutl      Split list by any delimiter
NB. info      Get information about mt interace name
NB. env       Get environment
NB. erasen    Erase global names created between invocations
NB.
NB. Copyright 2010,2011,2013,2017,2018,2020,2021,2023,2024,
NB.           2025 Igor Zhuravlov
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

NB. Format log string
NB. note: fix (8!:2) for complex [NaN] input
fmtlog=: ;@:(40 17 17 17 17 _16&(({.{.@('d<n/a>'&(8!:2 :: (,: 'n/a'))))&.>))

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Miscellaneous

NB. mark...
ispos0=:  _ =!.0 %  NB. ... +0  values in y
isneg0=: __ =!.0 %  NB. ... -0  values in y
isnan=: 128!:5      NB. ... NaN values in y

nan=:   _."_`(_.j_."_)@.(JCMPX = 3!:0)  NB. produce NaN of input datatype

max=: >./`0:@.(0 = #)  NB. max-of, 0 for empty list

negneg=:   ($@] $    (isneg0 +. 0&>)@{.`((,:  -)@{:)}@,:)`(($ $ nan)@])@.(+.&(isnan@<))  NB. if x<0 then - y  else  y  endif
negpos=:   ($@] $    (isneg0 +: 0&>)@{.`((,:  -)@{:)}@,:)`(($ $ nan)@])@.(+.&(isnan@<))  NB. if x≥0 then - y  else  y  endif
copysign=: ($@] $ =/&(isneg0 +. 0&>)   `((,:~ -)@{:)}@,:)`(($ $ nan)@])@.(+.&(isnan@<))  NB. if x<0 then -|y| else |y| endif

sorim=: | `(+/!.0"1@:| @:+.)@.(JCMPX = 3!:0)  NB. sum of real and imaginary parts' modules, |Re(y)| + |Im(y)|
soris=: *:`(+/!.0"1@:*:@:+.)@.(JCMPX = 3!:0)  NB. sum of real and imaginary parts' squares, Re(y)^2 + Im(y)^2

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. test suite utilities
NB.
NB. Notes:
NB. - see test log format in mt.ijs

NB. ---------------------------------------------------------
NB. lcat
NB.
NB. Description:
NB.   Concatenate logs
NB.
NB. Syntax:
NB.   logz=.      lcat logs
NB.   logz=. logx lcat logy
NB. where
NB.   logs      - matrix of boxes, each row is a test log
NB.   logx,logy - 6-vector of boxes, test log to concatenate
NB.   logz      - 6-vector of boxes, concatenated test log

lcat=: $:/ : (,&.>&:(,&.>/"2^:(<:@#@$)))

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
NB.   log -: log      lcat emptylog
NB.   log -: emptylog lcat log
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
NB.     vrcond=. nan`gecon1@.(=/@$)
NB.     NB. to calc. berr, assuming:
NB.     NB.   berr := ||A - realA||_1 / (FP_EPS * ||A||_1 * m)
NB.     vberr=. ((- %&norm1 [) % FP_EPS * (norm1 * #)@[) unmqr
NB.     NB. do the job
NB.     'sent rcond ferr berr time space'=. ('geqrf' tmonad ]`]`vrcond`nan`vberr) A
NB. - to test getrs:
NB.     NB. to estimate rcond in ∞-norm
NB.     vrcond=. nan`geconi@.(=/@$)@(0&{::)
NB.     NB. to calc. ferr, assuming:
NB.     NB.   ferr := ||x - realx||_inf / ||realx||_inf
NB.     vferr=. ((- %&normi [) 1&{::)~
NB.     NB. to calc. componentwise berr [LUG 75], assuming:
NB.     NB.   berr := max_i(|b - A * realx|_i / (|A| * |realx| + |b|)_i)
NB.     vberr=. (mp&>/@[ |@- (0 {:: [) mp ]) >./@% (((0 {:: [) mp&| ]) + |@mp&>/@[)
NB.     NB. do the job
NB.     'sent rcond ferr berr time space'=. ('getrs' tdyad (0&{::)`(mp&>/)`]`vrcond`vferr`vberr) (A;x)
NB.
NB. Notes:
NB. 1) recommended observations count to provide standard
NB.    deviation <= 1% for CPU run-time estimator minimum on
NB.    systems with load ≤ 80 is equal to 5 [1]
NB. 2) side effect: the result is sent to the console
NB.
NB. References:
NB. [1] Magne Haveraaen, Hogne Hundvebakke. Some Statistical
NB.     Performance Estimation Techniques for Dynamic
NB.     Machines. Appeared in Weihai Yu & al. (eds.): Norsk
NB.     Informatikk-konferanse 2001, Tapir, Trondheim Norway
NB.     2001, pp. 176-185.
NB.     https://www.ii.uib.no/saga/papers/perfor-5d.pdf

tmonad=: 2 : 0
  '`vgety vgeto vrcond vferr vberr'=. n
  try. rcond=. vrcond y catch. rcond=. _ end.
  try.
    ybak=. 1 memu argy=. vgety y
    try.
      't s'=. , ((5 1 # i. 2)) (<./`]/.) (5 1 # timex`(7!:2))`:0 'ret=. ' , m , ' argy'
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
    xbak=. 1 memu argx=. vgetx y
    ybak=. 1 memu argy=. vgety y
    try.
      't s'=. , ((5 1 # i. 2)) (<./`]/.) (5 1 # timex`(7!:2))`:0 'ret=. argx ' , m , ' argy'
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
NB. Examples:
NB. - when asserts are enabled:
NB.      9!:34 ''  NB. check asserts are enabled
NB.   1
NB.      NB. mt                               NB. stdlib
NB.      assert_mt_ 1 1 0                     assert_z_ 1 1 0
NB.   |assertion failure: assert_mt_       |assertion failure: assert_z_
NB.   |       assert_mt_ 1 1 0             |       assert_z_ 1 1 0
NB.      assert_mt_ 1 1 1                     assert_z_ 1 1 1
NB.      assert_mt_ 1 1 2                     assert_z_ 1 1 2          NB. no failure occured - it's a bug #1
NB.   |assertion failure: assert_mt_
NB.   |       assert_mt_ 1 1 2
NB.
NB.      'Oops!' assert_mt_ 1 1 0             'Oops!' assert_z_ 1 1 0
NB.   |Oops!: assert_mt_                   |Oops!: assert_z_
NB.   |   'Oops!'    assert_mt_ 1 1 0      |   'Oops!'    assert_z_ 1 1 0
NB.      'Oops!' assert_mt_ 1 1 1             'Oops!' assert_z_ 1 1 1
NB.      'Oops!' assert_mt_ 1 1 2             'Oops!' assert_z_ 1 1 2  NB. no failure occured - it's a bug #1
NB.   |Oops!: assert_mt_
NB.   |   'Oops!'    assert_mt_ 1 1 2
NB.
NB. - when asserts are disabled:
NB.      9!:35 [ 0  NB. disable asserts
NB.      9!:34 ''   NB. check asserts are disabled
NB.   0
NB.      NB. mt                               NB. stdlib
NB.      assert_mt_ 1 1 1                     assert_z_ 1 1 1
NB.      assert_mt_ 1 1 2                     assert_z_ 1 1 2
NB.      assert_mt_ 1 1 0                     assert_z_ 1 1 0          NB. failure occured - it's a bug #2
NB.                                        |assertion failure: assert_z_
NB.                                        |       assert_z_ 1 1 0
NB.
NB.      'Oops!' assert_mt_ 1 1 1             'Oops!' assert_z_ 1 1 1
NB.      'Oops!' assert_mt_ 1 1 2             'Oops!' assert_z_ 1 1 2
NB.      'Oops!' assert_mt_ 1 1 0             'Oops!' assert_z_ 1 1 0  NB. failure occured - it's a bug #2
NB.                                        |Oops!: assert_z_
NB.                                        |   'Oops!'    assert_z_ 1 1 0
NB.      9!:35 [ 1  NB. restore default setting
NB.      9!:34 ''   NB. check asserts are enabled
NB.   1
NB.
NB. Notes:
NB. - ambivalent procedure
NB. - fixes system's (assert_z_) to match (assert.) control
NB. - depends on 9!:34 (Enable assert.) setting
NB. - values of rank>1 are supported accidentally, too:
NB.      NB. mt                               NB. stdlib
NB.      assert_mt_ 1 1 ,: 1 0                assert_z_ 1 1 ,: 1 0  NB. no failure occured
NB.   |assertion failure: assert_mt_
NB.   |       assert_mt_ 1 1,:1 0
NB.
NB. References:
NB. [1] Igor Zhuravlov. [Jprogramming] assert verb from
NB.     stdlib mismatches assert. control
NB.     2019-12-30 00:43:46 UTC
NB.     https://www.jsoftware.com/pipermail/programming/2019-December/054693.html

assert=: 0 0 $ dbsig^:((1 +./@:~: ])`(12"_))^:(9!:34@'')

NB. ---------------------------------------------------------
NB. mrg
NB.
NB. Description:
NB.   Adv. to merge x and y by interleaving their items as
NB.   determined by a list m of shape x +&# y
NB.
NB. Examples:
NB.    'ABC' 0 1 1 0 1 0 mrg '012'
NB. A01B2C
NB.    2 0 1 2 2 1 0 1 0 mrg 'ABC' , '012' , 'xyz'
NB. xA0yz1B2C
NB.
NB. References:
NB. [1] JPhrases 3B. Merge & Amend
NB.     https://code.jsoftware.com/wiki/JPhrases/MergeAmend

mrg=: 1 : '/:@/:@(m"_) { ,'

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

NB. ---------------------------------------------------------
NB. info
NB.
NB. Description:
NB.   Get information about mt interace name
NB.
NB. Syntax:
NB.   i=. [s] info n
NB. where
NB.   s - string, optional, where every character specifies
NB.       which section must be included:
NB.         'd' - description
NB.         's' - syntax
NB.         'f' - formula
NB.         'l' - storage layout
NB.         'g' - algorithm
NB.         'a' - assertions
NB.         'e' - examples
NB.         'p' - application
NB.         'n' - notes
NB.         't' - TODO
NB.         'r' - references
NB.       default is '' (all sections)
NB.   n - string, an interface name
NB.   i - string, an information about name
NB.
NB. Examples:
NB.   NB. get info about 'trsmllnu_mt_' name
NB.   info_mt_ 'trsmllnu_mt_'
NB.   NB. the same as above
NB.   info_mt_ 'trsmllnu'
NB.   NB. get syntax and notes info about 'trsmllnu_mtbli_'
NB.   'sn' info_mt_ 'trsmllnu_mtbli_'
NB.
NB. Notes:
NB. - if (n) argument is a simple name not explicit locative
NB.   then it's searched in mt locale

info=: ''&$: : (4 : 0)
  ('section ''' , x , ''' is not recognized') assert '' -: x -. 'dsflgaepntr'
  i=. 4!:4 < y
  ('name ''' , y , ''' isn''t from script') assert * i
  content=. 1!:1 i { 4!:3 ''
  el=. '((_[[:alpha:]][[:alnum:]]*)+|_)_'       NB. pattern for explicit [chained] locale[s]
  nm=. '[[:alpha:]][[:alnum:]_]*(' , el , ')?'  NB. pattern for name
  ga=. '(' , nm , '=: .+\n)'                    NB. pattern for global assignment sentence
  y=. ((el , '$') ; '(' , el , ')?') rxrplc y   NB. replace locale name[s] by its pattern
  specs=. ('^(NB\..*\n)+(?=((\n' , ga , '+)*\n)?' , ga , '*' , y , '=: )') rxfirst content
  if. # x do.
    specs=. ; (<"0 x) (('^NB\. ' , '.*\n(NB\. .+\n)+?(NB\.\n|$)' ,~ (('Description' ; 'Syntax' ; 'Formula' ; 'Storage layout' ; 'Algorithm' ; 'Assertions' ; 'Examples' ; 'Application' ; 'Notes' ; 'TODO' ; 'References') {::~ 'dsflgaepntr' i. [)) rxfirst ])L:0 < specs
  end.
  }:^:(LF2 -: _2&{.) ('(?<=^|\n)NB\. ?';'') rxrplc specs
)

NB. ---------------------------------------------------------
NB. env
NB.
NB. Description:
NB.   Nilad to get environment
NB.
NB. Syntax:
NB.   e=. [t] env ''
NB. where
NB.   t - boolean, e type, optional, default is 1:
NB.         0 = return e as a vector of boxes with data
NB.         1 = return e as a string
NB.   e - string or 22-vector of boxes:
NB.          0 {:: e  - string, architecture
NB.          1 {:: e  > 0, integer, cores
NB.          2 {:: e  > 0, integer, max threads
NB.          3 {:: e  > 0, integer, worker threads
NB.          4..11    - 3-vector of non-negative integers,
NB.                     threads in threadpools 0..7:
NB.                       (#idle,#unfinished,#threads)
NB.         12 {:: e  ≥ 0, integer, executing thread# (0
NB.                     means master thread)
NB.         13..15    ≥ -1, integer, thresholds to switch
NB.                     (+/ .*) to GEMM (built-in BLIS
NB.                     implementation) for datatypes:
NB.                       13 {:: e  - integer
NB.                       14 {:: e  - floating
NB.                       15 {:: e  - complex
NB.         16..17    - strings about LAPACK interface:
NB.                       16 {:: e  - library [path]file name
NB.                       17 {:: e  - version string
NB.         18..19    - strings about BLAS interface:
NB.                       18 {:: e  - library [path]file name
NB.                       19 {:: e  - version string
NB.         20..21    - strings about BLIS interface:
NB.                       20 {:: e  - library [path]file name
NB.                       21 {:: e  - version string
NB.
NB. Application:
NB. - create 2 threads:
NB.     {{0 T.''}}^:2 ''
NB. - create (#cores - 1) threads in threadpool 0:
NB.     {{0 T.0}}^:] <: {. 8 T. ''
NB. - destroy 2 threads:
NB.     {{55 T.''}}^:2 ''
NB. - destroy (#cores - 1) threads in threadpool 0:
NB.     {{55 T.0}}^:] <: {. 8 T. ''
NB. - set threshold for floating matrices of size 1024×1024
NB.   or larger:
NB.     (<. 1024^3) (9!:58) 1
NB. - always use BLIS for any complex matrix:
NB.     0 (9!:58) 2
NB. - never use BLIS for any integer matrix:
NB.     _1 (9!:58) 0
NB. - try to load LAPACK interfaces if presented in system:
NB.     load 'math/lapack2'
NB. - try to load BLAS interfaces if presented in system:
NB.     load 'math/mt/external/blas/blas'
NB. - try to load BLIS interfaces if presented in system:
NB.     load 'math/mt/external/blis/blis'
NB.
NB. Notes:
NB. - see [1,2] about 9!:58 foreign
NB.
NB. References:
NB. [1] Bill Lam. [Jprogramming] J902 beta program started
NB.     2020-05-14 23:23:13 UTC
NB.     https://www.jsoftware.com/pipermail/programming/2020-May/055784.html
NB. [2] Bill Lam. [Jprogramming] J902 beta program started
NB.     2020-05-15 00:36:53 UTC
NB.     https://www.jsoftware.com/pipermail/programming/2020-May/055787.html

env=: 1&$: :(4 : 0)
  e=.     < 9!:56 'cpu'
  e=. e , < 9!:56 'cores'
  e=. e , < {: 8 T. ''
  e=. e , < 1 T. ''
  e=. e , < 2 T. 0
  e=. e , < 2 T. 1
  e=. e , < 2 T. 2
  e=. e , < 2 T. 3
  e=. e , < 2 T. 4
  e=. e , < 2 T. 5
  e=. e , < 2 T. 6
  e=. e , < 2 T. 7
  e=. e , < 3 T. ''
  e=. e , < (9!:58) 0
  e=. e , < (9!:58) 1
  e=. e , < (9!:58) 2
  e=. e , < (3 : 'liblapack_jlapack2_')`('n/a'"_)@.(nc@<@'liblapack_jlapack2_') ''
  e=. e , < ver_jlapack2_ :: 'n/a' ''
  e=. e , < (3 : 'LIB_mtbla_'         )`('n/a'"_)@.(nc@<@'LIB_mtbla_'         ) ''
  e=. e , < ver_mtbla_    :: 'n/a' ''
  e=. e , < (3 : 'LIB_mtbli_'         )`('n/a'"_)@.(nc@<@'LIB_mtbli_'         ) ''
  e=. e , < ver_mtbli_    :: 'n/a' ''
  if. x do.
    tpl1=. cut3 {{)n
Hardware
  threads in pool# (#idle,#unfinished,#threads):
Thresholds to switch (+/ .*) to GEMM (built-in BLAS implementation)
(switches if threshold <= m*n*p, _1 = switch never)
  for datatypes:
Interfaces to external libraries (optional)
  LAPACK
  BLAS
  BLIS
}}
    tpl2=. cut3 {{)n
  architecture:
  cores:
  maxthreads:
  worker threads:
    #0:
    #1:
    #2:
    #3:
    #4:
    #5:
    #6:
    #7:
  executing thread# (0 means master thread):
    integer :
    floating:
    complex :
    library:
    version:
    library:
    version:
    library:
    version:
}}
    e=. tpl2 ((, ' ' , ":)L:0) e             NB. merge pairs (title,value)
    e=. tpl1 ((_31 {. #: 1040129243)) mrg e  NB. interleave with value-less titles
    e=. LF joinstring e                      NB. raze with LF interleaved
  end.
)

NB. ---------------------------------------------------------
NB. erasen
NB.
NB. Description:
NB.   Erase global names created between invocations
NB.
NB. Syntax:
NB.   out=. [x] f&.erasen y
NB. where
NB.   out -: [x] f y
NB.
NB. Notes:
NB. - is an invertible ambivalent identity
NB. - erases names created while f was executed

erasen=: ([ 4!:5@1) :. ([ (4!:5@0)@erase@(4!:5@1))
