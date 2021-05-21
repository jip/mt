NB. Benchmark
NB.
NB. benchmark  Adv. to make dyad to benchmark mt addon
NB.
NB. Version: 0.13.0 2021-05-21
NB.
NB. Copyright 2017-2021 Igor Zhuravlov
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
NB. Conventions:
NB.   1) mt's benchmarking process uses matrices of size(s)
NB.      and generator given
NB.   2) benchmarking process produces scalar result for each
NB.      matrix size given
NB.   3) benchmarking result is a sum of estimated execution
NB.      durations (in seconds) for some mt's verbs set
NB.   4) verbs set is defined by mt's self-test verb test_mt_
NB.   5) an execution duration for each verb is estimated as
NB.      proposed in [1]
NB.
NB. Notes:
NB. - size of n×n-matrix:      100   150   200  500  600  800
NB.   space used, Mb: float:   0.125 0.25  0.5  2    4     8
NB.                   complex: 0.25  0.5   1    4    8    16
NB.
NB. References:
NB. [1] Magne Haveraaen, Hogne Hundvebakke. Some Statistical
NB.     Performance Estimation Techniques for Dynamic
NB.     Machines. Appeared in Weihai Yu & al. (eds.): Norsk
NB.     Informatikk-konferanse 2001, Tapir, Trondheim Norway
NB.     2001, pp. 176-185.
NB.     https://www.ii.uib.no/saga/papers/perfor-5d.pdf

NB. =========================================================
NB. Local definitions

NB. recommended observations count to provide standard
NB. deviation <= 1% for CPU run-time estimator minimum on
NB. systems with load ≤ 80 [1]
BMKRUNS=: 5

NB. remove trailing blanks
NB. reference: m144, JPhrases 4B. Locating & Selecting
rtb=: #~ ([: +./\. ' '&~:)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. benchmark
NB.
NB. Description:
NB.   Adv. to make dyad to benchmark mt addon using a matrix
NB.   generator given
NB.
NB. Syntax:
NB.   vapp=. mkmat benchmark
NB. where
NB.   mkmat - monad to generate a square matrix; is called
NB.           as:
NB.             mat=. mkmat 2 # n
NB.   vapp  - monad to benchmark mt addon algorithms with
NB.           matrix generator given; is called as:
NB.             out=. descr vapp sizes
NB.           formatted result (sizes ,. out) will also be
NB.           stored into file:
NB.             jpath '~temp/mt.benchmark.',descr,'.result'
NB.           also, a raw data derived from measurement and
NB.           used for calculation will be stored into file:
NB.             jpath '~temp/mt.benchmark.',descr,'.raw'
NB.           also, a raw measurement log used to form raw
NB.           data will be stored into file:
NB.             jpath '~temp/mt.benchmark.',descr,'.log'
NB.   n     > 0, the size of square matrices used for tests
NB.   descr - literal vector, benchmark description
NB.   sizes - s-vector, integer, test matrices sizes
NB.   out   - s-vector, float, benchmark results for corresp.
NB.           sizes
NB.   s     > 0, the number of sizes
NB.
NB. Application:
NB. - typical workflow to compare different hardware
NB.   environments performance for the same J interpreter
NB.   version:
NB.     # on one host:
NB.     user@i5-3230M:~/j64-805> ./jconsole.sh
NB.        load 'math/mt'
NB.        mkmat=. _1 1 0 3 _6 4&gemat_mt_
NB.        bmkrun=. mkmat benchmark_mt_
NB.        sizes=. 100 150 200 300 400 500
NB.        ] out1=. 'i5-3230M_j64-805' bmkrun sizes
NB.        NB. a lot of output, very long execution...
NB.     75.305 283.243 1040.378 9666.966 26840.712 66151.747
NB.        exit ''
NB.     # on another host:
NB.     user@i7-4790:~/j64-805> ./jconsole.sh
NB.        load 'math/mt'
NB.        mkmat=. _1 1 0 3 _6 4&gemat_mt_
NB.        bmkrun=. mkmat benchmark_mt_
NB.        sizes=. 100 150 200 300 400 500
NB.        ] out2=. 'i7-4790_j64-805' bmkrun sizes
NB.        NB. a lot of output, very long execution...
NB.     42.323 145.205 395.336 1782.877 8592.176 21163.120
NB.        NB. now we can compare a speed-up of 2nd host against 1st one:
NB.        sizes ,. 75.305 283.243 1040.378 9666.966 26840.712 66151.747 % out2
NB.     100 1.77929
NB.     150 1.95064
NB.     200 2.63163
NB.     300 5.42212
NB.     400 3.12386
NB.     500  3.1258
NB.        exit ''
NB. - typical workflow to compare different J interpreter
NB.   versions performance for the same hardware environment:
NB.     # in one interpreter:
NB.     user@i7-4790:~/j64-805> ./jconsole.sh
NB.        load 'math/mt'
NB.        mkmat=. _1 1 0 3 _6 4&gemat_mt_
NB.        bmkrun=. mkmat benchmark_mt_
NB.        sizes=. 100 150 200 300 400 500
NB.        ] out1=. 'i7-4790_j64-805' bmkrun sizes
NB.        NB. a lot of output, very long execution...
NB.     42.323 145.205 395.336 1782.877 8592.176 21163.120
NB.        exit ''
NB.     # in another interpreter:
NB.     user@i7-4790:~/j64-806> ./jconsole.sh
NB.        load 'math/mt'
NB.        mkmat=. _1 1 0 3 _6 4&gemat_mt_
NB.        bmkrun=. mkmat benchmark_mt_
NB.        sizes=. 100 150 200 300 400 500
NB.        ] out2=. 'i7-4790_j64-806' bmkrun sizes
NB.        NB. a lot of output, very long execution...
NB.     43.179 150.100 400.637 1814.688 8729.023 20626.823
NB.        NB. now we can compare a speed-up of 2nd interpreter against 1st one:
NB.        sizes ,. 42.323 145.205 395.336 1782.877 8592.176 21163.120 % out2
NB.     100 0.980176
NB.     150 0.967388
NB.     200 0.986769
NB.     300  0.98247
NB.     400 0.984323
NB.     500    1.026
NB.        exit ''

benchmark=: 1 : 0
:
  datatypename=. datatype u 10

  timespacex_bak=.  timespacex_z_
  TESTLOGFILE_bak=. TESTLOGFILE_mt_
  TESTLOG_bak=.     TESTLOG_mt_

  timespacex_z_=: timex , (_."_)
  TESTLOGFILE_mt_=: < jpath '~temp/mt.benchmark.' , x , '.log'

  measurements=. i. 0 0 0  NB. (# y)×(# verbs)×BMKRUNS-brick of boxed strings 'n/a' or with formatted float
  titles=. ''  NB. (# y)-vector of boxed titles
  iosize=. 0
  while. iosize < # y do.
    size=. iosize { y
    measurement=. ''  NB. (# verbs)×BMKRUNS-matrix of boxed strings 'n/a' or with formatted float
    titles=. titles , < (6!:0 'YYYY-MM-DD hh:mm:ss') , ', datatype: ' , datatypename , ', size: ' , (": size) , LF
    run=. 1
    while. run <: BMKRUNS_mt_ do.
      echo (6!:0 'YYYY-MM-DD hh:mm:ss') , ', datatype: ' , datatypename , ', size: ' , (": size) , ', run#: ' , (": run) , ' of ' , (": BMKRUNS_mt_)
      TESTLOG_mt_=: ''
      (u test_mt_) 2 # size
      measurement=. measurement ,. :: (,.@]) '(\d+\.\d+|n/a)(?=\s+(n/a)$)' rxall TESTLOG_mt_
      run=. >: run
    end.
    measurements=. measurements , measurement
    iosize=. >: iosize
  end.

  NB. extract verb names from last run
  measuredverbs=. rtb_mt_ L: 0 '^.{40}'&rxall TESTLOG_mt_

  NB. not all measured verbs belongs to mt addon, so filter aliens out
  iofellows=. measuredverbs e. namelist_mt_ 3
  fellowsnames=. iofellows # measuredverbs
  measurements=. iofellows&(#"2) measurements

  NB. save pretty-printed raw results; optional
  (titles ;@,. <@clipfmt"2 fellowsnames ,."1 2 measurements) fwrite '~temp/mt.benchmark.' , x , '.raw'

  NB. replace 'n/a' by _. , convert to float, unbox
  measurements=. (measurements = < 'n/a')} measurements ,: < '_.'
  measurements=. > ". L: 0 measurements  NB. (# y)×(# verbs)×BMKRUNS-brick of floats and possibly NaNs

  NB. for each verb: exclude NaNs, take minimum execution time [1]
  rounds=. <./@(#~ -.@(128!:5))"1 measurements  NB. (# y)×(# verbs)-matrix of floats

  NB. for each size: estimate benchmark value
  out=. +/"1 rounds  NB. (# y)-vector of floats

  NB. prepare brick and save
  (clipfmt y ,. out) fwrite '~temp/mt.benchmark.' , x , '.result'

  timespacex_z_=: timespacex_bak f.
  TESTLOGFILE_mt_=: TESTLOGFILE_bak
  TESTLOG_mt_=: TESTLOG_bak

  out
)
