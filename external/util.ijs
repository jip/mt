NB. Tests' utilities
NB.
NB. issquare  Same as in the (math/lapack2) addon
NB. basicxxx  Utilities to either check or modify argument
NB. initnoun  Define global noun and return its value
NB. dlsym     Obtain address of a symbol in a shared object
NB.           or executable
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

coclass 'mttst'
coinsert 'mt'

NB. =========================================================
NB. Interface

issquare=: =/@$  NB. same as in the (math/lapack2) addon

NB. check
NB. - ranks
basiccr0=: 0 2 2         -: #@$S:0
basiccr1=: 2 1 0         -: #@$S:0
basiccr2=: 0 1 0 2       -: #@$S:0
basiccr3=: 0 2 0 2       -: #@$S:0
basiccr4=: 0 2 2 0 2     -: #@$S:0
basiccr5=: 0 1 0 1 0 2   -: #@$S:0
basiccr6=: 0 2 1 0 0 1 0 -: #@$S:0
NB. - shape
basiccs0=: issquare@(0&{::)
basiccs1=: issquare@(1&{::)
basiccs3=: issquare@(3&{::)
basiccs4=: issquare@(4&{::)
basiccs5=: issquare@(5&{::)
NB. - compare shapes
basiccmp=: -:/@($L:0)@(1 2&{)

NB. modify
NB. - conjugate under ISO specified
basiccj0=: +&.>&.(1      &{)
basiccj1=: +&.>&.(0 1 3  &{)
basiccj2=: +&.>&.(0 2 4 5&{)
NB. - swap elements
basicswp=: (< 1 2)&C.

NB. ---------------------------------------------------------
NB. initnoun
NB.
NB. Description:
NB.   Define global noun and return its value.
NB.
NB. Syntax:
NB.   val=. name initnoun val
NB. where
NB.   name - string, global noun's name
NB.   val  - noun's value to initialize

initnoun=: 4 : '(x)=: y'

NB. ---------------------------------------------------------
NB. dlsym
NB.
NB. Description:
NB.   Obtain address of a symbol in a shared object or
NB.   executable
NB.
NB. Syntax:
NB.   'addr errmsg'=. dlsym libpath ; symname
NB. where
NB.   libpath - string, library's FQFN
NB.   symname - string, symbol name
NB.   addr    - integer, an address or 0 if symbol is not
NB.             found
NB.   errmsg  - string, an error message if addr=0 or ''

dlsym=: 3 : 0
  'lib sym'=. y
  select. UNAME
    case. 'Win' do.
      dl=. 'kernel32.dll'
      sig=. ' GetProcAddress x x *c'
      h=. 0 {:: (dl , ' LoadLibrary * *c') cd < lib
    case. 'Linux' ; 'OpenBSD' ; 'FreeBSD' do.
      dl=. '/usr/lib' , (IF64 {:: '' ; '64') , '/libdl.so'
      sig=. ' dlsym x x *c'
      h=. 0 {:: (dl , ' dlopen * *c i') cd lib ; 1  NB. lazy binding
    case. 'Android' ; 'Darwin' ; 'Unknown' ; 'Wasm' do.
      ('UNAME ' , UNAME , ' isn''t supported yet') dbsig 11
    case. do.
      ('UNAME ' , UNAME , ' isn''t recognized'   ) dbsig 11
  end.
  if. h do.
    p=. 0 {:: (dl , sig) cd h ; sym
  else.
    p=. 0  NB. cannot open lib
  end.
  (<p) 0} cderx ''
)
