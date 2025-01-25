NB. Utilities
NB.
NB. JINT2        integer2 datatype ID
NB. JINT4        integer4 datatype ID
NB.
NB. symhdr       Get address of header for named value
NB. memr         Advanced version of the (memr_z_) verb
NB. memw         Advanced version of the (memw_z_) verb
NB. meme         Adv. to make bivalent verb to edit memory
NB. obja         Allocate BLIS object for noun
NB. objf         Free BLIS object related with noun
NB. ver          Get version string
NB. arch         Get a string that contains the name of the
NB.              configuration
NB. *_get_*_str  Get a general configuration string
NB. info_get_*_impl_string
NB.              Get a string with implementation type for an
NB.              operation given
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

coclass 'mtbli'

NB. =========================================================
NB. Local definitions

NB. =========================================================
NB. Interface

NB. Datatype IDs
JINT2=: 6  NB. integer2
JINT4=: 7  NB. integer4

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Memory management

symhdr=: 15!:12  NB. get address of header for named value

NB. ---------------------------------------------------------
NB. memr
NB.
NB. Description:
NB.   Advanced version of the (memr_z_) verb
NB.
NB. Syntax:
NB.   data=. memr address , byte_offset , count [, type]
NB. where
NB.   type - scalar, optional, any of the following:
NB.                 1 or JB01          NB. boolean
NB.                 2 or JCHAR         NB. literal, default
NB.                 4 or JINT or JPTR  NB. integer
NB.                 5                  NB. integer1
NB.                 6 or JINT2         NB. integer2
NB.                 7 or JINT4         NB. integer4
NB.                 8 or JFL           NB. floating
NB.                 9                  NB. floating2
NB.                10                  NB. floating4
NB.                11                  NB. floating16
NB.                16 or JCMPX         NB. complex
NB.            131072 or JCHAR2        NB. unicode
NB.            262144 or JCHAR4        NB. unicode4
NB.
NB. Notes:
NB. - extends system's (memr_z_) to match (datatype_z_) verb

memr=: 3 : 0
  typ=. 3&{ :: (JCHAR"_) y
  select. typ
    case. JINT2 do.
      dat=. JINT2 c. _1 ic memr_z_ (JCHAR ,~ +:@ {.)&.(2 3&{) y  NB. sizeof(JINT2) == 2*sizeof(JCHAR)
    case. JINT4 do.
      dat=. JINT4 c. _2 ic memr_z_ (JCHAR ,~ 4 * {.)&.(2 3&{) y  NB. sizeof(JINT4) == 4*sizeof(JCHAR)
    NB. case. 5 ; 9 ; 10 ; 11 do.
    NB.   ('memr_mtbli_: datatype ' , (":typ) , ' isn''t supported yet') dbsig 11
    case. do.
      dat=. memr_z_ y
  end.
  dat
)

NB. ---------------------------------------------------------
NB. memw
NB.
NB. Description:
NB.   Advanced version of the (memw_z_) verb
NB.
NB. Syntax:
NB.   trash=. data memw address , byte_offset , count [, type]
NB. where
NB.   type - scalar, optional, any of the following:
NB.                 1 or JB01          NB. boolean
NB.                 2 or JCHAR         NB. literal, default
NB.                 4 or JINT or JPTR  NB. integer
NB.                 5                  NB. integer1
NB.                 6 or JINT2         NB. integer2
NB.                 7 or JINT4         NB. integer4
NB.                 8 or JFL           NB. floating
NB.                 9                  NB. floating2
NB.                10                  NB. floating4
NB.                11                  NB. floating16
NB.                16 or JCMPX         NB. complex
NB.            131072 or JCHAR2        NB. unicode
NB.            262144 or JCHAR4        NB. unicode4
NB.
NB. Notes:
NB. - extends system's (memw_z_) to match (datatype_z_) verb

memw=: 4 : 0
  typ=. 3&{ :: (JCHAR"_) y
  select. typ
    case. JINT2 do.
      (1 ic x) memw_z_ (JCHAR ,~ +:@ {.)&.(2 3&{) y  NB. sizeof(JINT2) == 2*sizeof(JCHAR)
    case. JINT4 do.
      (2 ic x) memw_z_ (JCHAR ,~ 4 * {.)&.(2 3&{) y  NB. sizeof(JINT4) == 4*sizeof(JCHAR)
    NB. case. 5 ; 9 ; 10 ; 11 do.
    NB.   ('memw_mtbli_: datatype ' , (":typ) , ' isn''t supported yet') dbsig 11
    case. do.
      x memw_z_ y
  end.
  EMPTY
)

NB. ---------------------------------------------------------
NB. NB. conj. to make bivalent verb to read bytes from memory y,
NB. NB.   decode by v, process by u [with left argument x],
NB. NB.   encode back by inv(v), write back to y
NB. NB. note: a (bivalent) conj. from
NB. NB.   addons/misc/miscutils/langexten.ijs is used here
NB. NB.   inlined
NB. meme=: 2 : 'u&.v^:(1:`(] memr_mtbli_)) memw_mtbli_ ]'

NB. adv. to make bivalent verb to read bytes from memory y,
NB.   process by u [with left argument x], write back to y
NB. note: a (bivalent) conj. from
NB.   addons/misc/miscutils/langexten.ijs is used here
NB.   inlined
meme=: 1 : 'u^:(1:`(] memr_mtbli_)) memw_mtbli_ ]'

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
  yobj=. mema OBJ_T_SIZE
  ytyp=. (4 5 6 7 8 16 i. (3!:0) y) { :: (('obja: datatype ' , (datatype y) , ' isn''t supported yet')&dbsig@11) INT , INT , INT , INT , DOUBLE , DCOMPLEX
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
  >: meme yhdr , (4 * SZI) , 1 , JINT  NB. refcount++
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

objf=: EMPTY [ (<: meme)@(,&((4 * SZI) , 1 , JINT))`memf"0

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Service functions

NB. ---------------------------------------------------------
NB. General library information

NB. Get version string
NB. str=. ver ''
ver=: (0 , JSTR) memr@,~ info_get_version_str_cd

NB. ---------------------------------------------------------
NB. Specific configuration

NB. Get a string that contains the name of the configuration
NB. str=. arch ''
arch=: (0 , JSTR) memr@,~ arch_string_cd@arch_query_id_cd

NB. ---------------------------------------------------------
NB. General configuration

NB. Get a string with the size of gint_t signed integer
NB. str=. info_get_int_type_size_str ''
info_get_int_type_size_str=: (0 , JSTR) memr@,~ info_get_int_type_size_str_cd

NB. Get a string with the thread implementation method
NB. str=. thread_get_thread_impl_str ''
thread_get_thread_impl_str=: (0 , JSTR) memr@,~ thread_get_thread_impl_str_cd@thread_get_thread_impl_cd

NB. ---------------------------------------------------------
NB. Kernel information

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Micro-kernel implementation type query

NB. str=. info_get_gemm_ukr_impl_string method ; dt
info_get_gemm_ukr_impl_string=:       (0 , JSTR) memr@,~ info_get_gemm_ukr_impl_string_cd

NB. str=. info_get_gemmtrsm_l_ukr_impl_string method ; dt
info_get_gemmtrsm_l_ukr_impl_string=: (0 , JSTR) memr@,~ info_get_gemmtrsm_l_ukr_impl_string_cd

NB. str=. info_get_gemmtrsm_u_ukr_impl_string method ; dt
info_get_gemmtrsm_u_ukr_impl_string=: (0 , JSTR) memr@,~ info_get_gemmtrsm_u_ukr_impl_string_cd

NB. str=. info_get_trsm_l_ukr_impl_string method ; dt
info_get_trsm_l_ukr_impl_string=:     (0 , JSTR) memr@,~ info_get_trsm_l_ukr_impl_string_cd

NB. str=. info_get_trsm_u_ukr_impl_string method ; dt
info_get_trsm_u_ukr_impl_string=:     (0 , JSTR) memr@,~ info_get_trsm_u_ukr_impl_string_cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Operation implementation type query

NB. str=. info_get_gemm_impl_string dt
info_get_gemm_impl_string=:  (0 , JSTR) memr@,~ info_get_gemm_impl_string_cd

NB. str=. info_get_gemmt_impl_string dt
info_get_gemmt_impl_string=: (0 , JSTR) memr@,~ info_get_gemmt_impl_string_cd

NB. str=. info_get_hemm_impl_string dt
info_get_hemm_impl_string=:  (0 , JSTR) memr@,~ info_get_hemm_impl_string_cd

NB. str=. info_get_herk_impl_string dt
info_get_herk_impl_string=:  (0 , JSTR) memr@,~ info_get_herk_impl_string_cd

NB. str=. info_get_her2k_impl_string dt
info_get_her2k_impl_string=: (0 , JSTR) memr@,~ info_get_her2k_impl_string_cd

NB. str=. info_get_symm_impl_string dt
info_get_symm_impl_string=:  (0 , JSTR) memr@,~ info_get_symm_impl_string_cd

NB. str=. info_get_syrk_impl_string dt
info_get_syrk_impl_string=:  (0 , JSTR) memr@,~ info_get_syrk_impl_string_cd

NB. str=. info_get_syr2k_impl_string dt
info_get_syr2k_impl_string=: (0 , JSTR) memr@,~ info_get_syr2k_impl_string_cd

NB. str=. info_get_trmm_impl_string dt
info_get_trmm_impl_string=:  (0 , JSTR) memr@,~ info_get_trmm_impl_string_cd

NB. str=. info_get_trmm3_impl_string dt
info_get_trmm3_impl_string=: (0 , JSTR) memr@,~ info_get_trmm3_impl_string_cd

NB. str=. info_get_trsm_impl_string dt
info_get_trsm_impl_string=:  (0 , JSTR) memr@,~ info_get_trsm_impl_string_cd
