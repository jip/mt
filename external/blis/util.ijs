NB. Utilities
NB.
NB. symhdr             Get address of header for named value
NB. updmem             Adv. to make bivalent verb to update
NB.                    memory
NB. obja               Allocate BLIS object for noun
NB. objf               Free BLIS object related with noun
NB. int_type_size_str  Query a string with the size of gint_t
NB.                    signed integer
NB. arch               Query a string that contains the name
NB.                    of the configuration
NB. ver                Get version string
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

coclass 'mtbli'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

symhdr=: 15!:12  NB. get address of header for named value

NB. ---------------------------------------------------------
NB. NB. conj. to make bivalent verb to read bytes from memory y,
NB. NB.   decode by v, process by u [with left argument x],
NB. NB.   encode back by inv(v), write back to y
NB. NB. note: a (bivalent) conj. from
NB. NB.   addons/misc/miscutils/langexten.ijs is used here
NB. NB.   inlined
NB. updmem=: 2 : 'u&.v^:(1:`(] memr)) memw ]'

NB. adv. to make bivalent verb to read bytes from memory y,
NB.   process by u [with left argument x], write back to y
NB. note: a (bivalent) conj. from
NB.   addons/misc/miscutils/langexten.ijs is used here
NB.   inlined
updmem=: 1 : 'u^:(1:`(] memr)) memw ]'

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. OAPI functions

NB. ---------------------------------------------------------
NB. obja
NB.
NB. Description:
NB.   Allocate BLIS object for noun
NB.
NB. Syntax:
NB.   'hdr obj'=. obja noun
NB. where
NB.   noun - any noun
NB.   hdr  - integer, pointer to an enveloping noun's header
NB.   obj  - integer, pointer to an allocated memory
NB.
NB. Storage layout:
NB.          J symbol table
NB.   sym -> (name,hdr)                  NB. result of symget
NB.
NB.          J symbol headers
NB.   hdr -> (offset,...,len,...)        NB. result of symhdr
NB.
NB.          J symbol data
NB.   dat -> data_of_length_len          NB. result of symdat
NB.   dat = (hdr+offset)
NB.
NB.          J heap
NB.   obj -> BLISobject(...,buffer,...)  NB. result of obja
NB.   buffer := dat
NB.
NB. Application:
NB. - allocate BLIS objects in a batch:
NB.     'obj0 obj1'=. obja_mtbli_ S: 0 (noun0 ; noun1)
NB.
NB. Notes:
NB. - side effect: (refcount++) in enveloping noun

obja=: 3 : 0
  yhdr=. symhdr < 'y'
  ydat=. symdat < 'y'
  yobj=. mema SIZEOF_OBJ_T
  ytyp=. (4 8 16 i. (3!:0) y) { :: (('obja: datatype ' , (datatype y) , ' isn''t supported yet')&dbsig@11) _2 ic ((INT , DOUBLE , DCOMPLEX))
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
  >: updmem yhdr , (((4 * SZI) , 1 , JINT))  NB. refcount++
  yhdr , yobj
)

NB. ---------------------------------------------------------
NB. objf
NB.
NB. Description:
NB.   Procedure to free BLIS object related with noun
NB.
NB. Syntax:
NB.   trash=. objf (hdr , obj)
NB. where
NB.   hdr - integer, pointer to an enveloping noun's header
NB.   obj - integer, pointer to an allocated memory
NB.
NB.
NB. Application:
NB. - free BLIS objects in a batch:
NB.     trash=. objf_mtbli_"1 (obj0 ,: obj1)
NB.
NB. Notes:
NB. - side effect: (refcount--) in enveloping noun

objf=: EMPTY [ (<: updmem)@(,&(((4 * SZI) , 1 , JINT)))`memf"0

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
