NB. Utilities
NB.
NB. max        Max-of, 0 for empty list
NB. maxc       Max-of, '' for empty list
NB. negneg     Conditional negate
NB. negpos     Conditional negate
NB. copysign   Copy sign
NB. sorim      Sum of real and imaginary parts' modules
NB. soris      Sum of real and imaginary parts' squares
NB. fmtlog     Format log string
NB. assert     Advanced version of the (assert.) control
NB. cut3       Split list by delimiter taken from its tail
NB. cut2       Split list by delimiter
NB. cut        Split list by delimiter
NB. cutl2      Split list by any delimiter
NB. cutl       Split list by any delimiter
NB. benchmark  Adv. to make ambivalent verb to benchmark
NB.            sentences using matrices of generator and
NB.            shape given
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

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Miscellaneous

max=:  >./`      0: @.(0 = #)                                   NB. max-of, 0 for empty list
maxc=: >./`(c {. 0:)@.(0 = #)                                   NB. max-of, '' for empty list

negneg=: -@]^:(0>[)                                             NB. if x<0 then -y else y endif
negpos=: -@]^:(0<:[)                                            NB. if x≥0 then -y else y endif

copysign=: (=/&:*`((,:~ -)@{:))}@,:                             NB. if x<0 then -|y| else |y| endif

sorim=: | `(+/"1@:| @:+.)@.(JCMPX = 3!:0)                       NB. sum of real and imaginary parts' modules, |Re(y)| + |Im(y)|
soris=: *:`(+/"1@:*:@:+.)@.(JCMPX = 3!:0)                       NB. sum of real and imaginary parts' squares, Re(y)^2 + Im(y)^2

fmtlog=: ;@:(40 17 17 17 17 _16&(({.{.@('d<n/a>'&(8!:2)))&.>))  NB. log string format

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

NB. ---------------------------------------------------------
NB. benchmark
NB.
NB. Description:
NB.   Adv. to make ambivalent verb to benchmark sentences
NB.   using matrices of generator and shape given
NB.
NB. Syntax:
NB.   d=. [rx] (mkmat atest) benchmark (m,n)
NB. where
NB.   mkmat - monad to generate a material for test matrices;
NB.           is called as:
NB.             mat=. mkmat (m,n)
NB.   atest - adv. to make monadic procedure to run tests; is
NB.           called as:
NB.             trash=. (mkmat atest) (m,n)
NB.   rx    - string, optional, a regular expression to
NB.           filter out sentences, default is:
NB.             ^((?!_mt(mm|tmp|lap|fla|bl[ai])_\b|\b128!:([01]|10)\b|%\.|\+\/ \. ?\*).)+
NB.           meaning: benchmark all sentences except ones
NB.           listed below:
NB.             128!:0
NB.             128!:1
NB.             128!:10
NB.             +/ .*
NB.             +/ . *
NB.             %.
NB.             with name which contains any of the next:
NB.               _mtmm_
NB.               _mttmp_
NB.               _mtlap_
NB.               _mtfla_
NB.               _mtbla_
NB.               _mtbli_
NB.   (m,n) - a shape of test matrices to be used by tests
NB.   d     > 0, a total duration (in seconds) of tested
NB.           sentences execution time
NB.
NB. Notes:
NB. - side effects:
NB.    - augments the TESTLOG_mt_ global noun
NB.    - outputs to the console
NB.
NB. Application:
NB.      load 'math/mt'
NB.      ] sizes=. 100 * #\ i. 5
NB.   100 200 300 400 500
NB.      NB. benchmark solvers by real matrices with
NB.      NB. elements distributed uniformly with support
NB.      NB. (0,1) and by 3 RHS
NB.      ] times=. ?@$&0 testsv_mt_ benchmark"1 sizes ,. 3
NB.   12345 123456 1234567 12345678 123456789
NB.      NB. benchmark triangular factorizators by square
NB.      NB. complex matrices
NB.      mkmat=. gemat_mt_ j. gemat_mt_
NB.      ] times=. mkmat testtrf_mt_ benchmark"1 ,.~ sizes
NB.   54321 654321 7654321 87654321 987654321
NB.      exit ''

benchmark=: 1 : 0
  '^((?!_mt(mm|tmp|lap|fla|bl[ai])_\b|\b128!:([01]|10)\b|%\.|\+\/ \. ?\*).)+' u benchmark_mt_ y
:
  bkp=. TESTLOG_mt_
  TESTLOG_mt_=: 1 5 # EMPTY ; ''
  u y
  ds=. +/!.0 ". S: 0 a: -.~ '(\d+\.\d+)(?=\s+(n/a|\d+)$)'&rxfirst L: 0 (#~ x&rxeq S: 0) cut3_mt_ TESTLOG_mt_  NB. FIXME
  TESTLOG_mt_=: bkp ,&.> TESTLOG_mt_
  ds
)

