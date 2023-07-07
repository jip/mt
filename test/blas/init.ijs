NB. Core utilities
NB.
NB. LIB  Noun, string, full path to BLAS library
NB.
NB. Version: 0.14.0 2023-03-21
NB.
NB. Copyright 2010-2023 Igor Zhuravlov
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

coclass 'mtbla'
coinsert 'mt'

NB. =========================================================
NB. Local definitions

NB. ---------------------------------------------------------
NB. findlib
NB.
NB. Description:
NB.   Nilad to find BLAS dynamic library and return the FQFN
NB.   of the first library found.
NB.
NB. Syntax:
NB.   fqfn=. findlib ''
NB. where
NB.   fqfn - string, FQFN of the first library found
NB.
NB. Notes:
NB. - will throw exception if no library was found
NB. - lib name has priority over path while lookup i.e.
NB.   search pairs (path[i],name[j]) are ordered as:
NB.     foreach_j(foreach_i(lookup(path[i],name[j])))

findlib=: 3 : 0
  select. UNAME
    case. 'Win' do.
      paths=. < jpath '~'
      names=. 'libblas.dll' ; 'librefblas.dll' ; 'liblapack3.dll' ; 'libopenblas.dll' ; 'libopenblas_serial.dll' ; 'libopenblas_pthreads.dll' ; 'libblis.dll' ; 'libflame.dll'
    case. 'Linux' ; 'OpenBSD' ; 'FreeBSD' do.
      h=. getenv 'HOME'
      paths=. (h ,L:0 '/lib' ; '/blas/lib' ; '/openblas/lib') , '/opt/OpenBLAS' ; '/opt/OpenBLAS/lib' ; '/usr/local/lib' ; '/usr/lib'
      arch=. '' ;~ IF64 { '32' ,: '64'
      paths=. , ;L:_1 { paths ,&< arch
      names=. 'libblas.so' ; 'librefblas.so' ; 'liblapack3.so' ; 'libopenblas.so' ; 'libopenblas_serial.so' ; 'libopenblas_pthreads.so' ; 'libblas.so.3' ; 'librefblas.so.3' ; 'liblapack3.so.3' ; 'libopenblas.so.0' ; 'libopenblas_serial.so.0' ; 'libopenblas_pthreads.so.0' ; 'libblis.so' ; 'libflame.so'
    case. 'Android' do.
      arch=. LF -.~ 2!:0 'getprop ro.product.cpu.abi'
      arch=. ((IF64 { ('arm64-v8a' ; 'x86_64')&i. , 2:) {:: ('armeabi-v7a' ; 'x86')&,) < arch  NB. fix it
      paths=. (jpath '~bin/../libexec/' , arch) ; ({.~ i:&'/') LIBFILE                         NB. ({.~ i:&'/') is the same as fpath_j_
      names=. 'libblas.so' ; 'librefblas.so' ; 'liblapack3.so' ; 'libopenblas.so' ; 'libopenblas_serial.so' ; 'libopenblas_pthreads.so' ; 'libblas.so.3' ; 'librefblas.so.3' ; 'liblapack3.so.3' ; 'libopenblas.so.0' ; 'libopenblas_serial.so.0' ; 'libopenblas_pthreads.so.0' ; 'libblis.so' ; 'libflame.so'
    case. 'Darwin' ; 'Unknown' ; 'Wasm' do.
      ('UNAME ' , UNAME , ' isn''t supported yet') dbsig 11
    case. do.
      ('UNAME ' , UNAME , ' isn''t recognized'   ) dbsig 11
  end.
  paths=. termsep_j_ L: 0 (#~ fexist S: 0) paths
  fqfns=. (#~ (0 ~: 0 {:: dlsym@(;&'caxpy_')) S: 0) , |: ;L:_1 { paths ,&< names
  ('BLAS library file(s) ' , (; ,&' 'L:0 names) , 'was(-ere) not found') assert * # fqfns
  0 {:: fqfns
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Init LIB global noun by full path of first library found.
NB. Throw error if no library was found.
NB.
NB. Notes:
NB. - OS-specific

('LIB_mtbla_' initnoun findlib)^:(0 ~: nc@<@'LIB_mtbla_') ''
