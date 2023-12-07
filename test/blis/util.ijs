NB. Utilities
NB.
NB. obja               Allocate BLIS object for noun
NB. objf               Free BLIS object
NB. int_type_size_str  Query a string with the size of gint_t
NB.                    signed integer
NB. arch               Query a string that contains the name
NB.                    of the configuration
NB. ver                Get version string
NB.
NB. Version: 0.14.0 2023-07-04
NB.
NB. Copyright 2023 Igor Zhuravlov
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

coclass 'mtbli'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. OAPI functions

NB. allocate BLIS object for noun
NB. obj=. obja noun
NB. TODO:
NB.   change interface to accept multiple nouns
NB.   allocate objects for noun
NB.   'obj0 obj1 ...'=. o4n noun0 ; noun1 ; ...

obja=: 3 : 0
  ydat=. symdat < 'y'
  yobj=. mema SIZEOF_OBJ_T
  ytyp=. (8 16 i. 3!:0 y) { :: (('obja: datatype ' , (datatype y) , ' isn''t supported yet')&dbsig@11) _2 ic DOUBLE , DCOMPLEX
  select. # $ y
    case. 2 do.
      'm n'=. $ y
      obj_create_with_attached_buffer_cd     ytyp ; m ; n ; (< ydat) ; n ; 1 ; < < yobj
    case. 1 do.
      n=. # y
      obj_create_with_attached_buffer_cd     ytyp ; 1 ; n ; (< ydat) ; 1 ; 1 ; < < yobj
    case. 0 do.
      obj_create_1x1_with_attached_buffer_cd ytyp ;         (< ydat) ;         < < yobj
    case.   do.
      ('obja: rank ' , (": # $ y) , ' isn''t supported') dbsig 11
  end.
  yobj
)

NB. free BLIS object by address given
objf=: memf

NB. ---------------------------------------------------------
NB. Service functions

NB. query a string with the size of gint_t signed integer
NB. strSize=. int_type_size_str_mtbli_ ''
int_type_size_str=: (0 , JSTR) memr@,~ info_get_int_type_size_str_cd

NB. query a string that contains the name of the
NB.   configuration
NB. strArch=. arch_mtbli_ ''
arch=: (0 , JSTR) memr@,~ arch_string_cd@arch_query_id_cd

NB. Get version string
NB. strVer=. ver_mtbli_ ''
ver=: (0 , JSTR) memr@,~ info_get_version_str_cd
