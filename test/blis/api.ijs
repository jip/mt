NB. API definitions
NB.
NB. xxxxxcd  Cover verbs to call BLIS subroutine or function
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
NB. Concepts
NB.
NB. Conventions:
NB. 1) LIB_mtbli_ global noun must exist
NB. 2) bit fields are presented as byte vector

NB. =========================================================
NB. Configuration

coclass 'mtbli'
coinsert 'mt'

NB. =========================================================
NB. Local definitions

lib=. dquote LIB
ifw=. IFWIN # '+'

NB. bitwise operations on bytes
NB. note: avoids conversion to integer
and=: (2b010001 b.)&.(a.&i.)
xor=: (2b010110 b.)&.(a.&i.)
or=:  (2b010111 b.)&.(a.&i.)
not=: (2b011010 b.)&.(a.&i.)
sft=: (2b100001 b.)&.(a.&i.)  NB. unsigned bitwise shift

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

NB. =========================================================
NB. Interface

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. BLIS types

NB. ---------------------------------------------------------
NB. Enumerated argument types

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. num_t
NB. Semantic meaning: Matrix/vector operand...

FLOAT=:             2 ic 2b00000000000000000000000000000000  NB. ...contains single-precision real elements.
SCOMPLEX=:          2 ic 2b00000000000000000000000000000001  NB. ...contains single-precision complex elements.
DOUBLE=:            2 ic 2b00000000000000000000000000000010  NB. ...contains double-precision real elements.
DCOMPLEX=:          2 ic 2b00000000000000000000000000000011  NB. ...contains double-precision complex elements.
INT=:               2 ic 2b00000000000000000000000000000100  NB. ...contains integer elements of type gint_t.
CONSTANT=:          2 ic 2b00000000000000000000000000000101  NB. ...contains polymorphic representation of a constant value.

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. dom_t
NB. Semantic meaning: Matrix/vector operand...

REAL=:              2 ic 2b00000000000000000000000000000000  NB. ...contains real domain elements.
COMPLEX=:           2 ic 2b00000000000000000000000000000001  NB. ...contains complex domain elements.

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. prec_t
NB. Semantic meaning: Matrix/vector operand...

SINGLE_PREC=:       2 ic 2b00000000000000000000000000000000  NB. ...contains single-precision elements.
DOUBLE_PREC=:       2 ic 2b00000000000000000000000000000010  NB. ...contains double-precision elements.

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. trans_t
NB. Semantic meaning: Corresponding matrix operand...

NO_TRANSPOSE=:      2 ic 2b00000000000000000000000000000000  NB. ...will be used as given.
TRANSPOSE=:         2 ic 2b00000000000000000000000000001000  NB. ...will be implicitly transposed.
CONJ_NO_TRANSPOSE=: 2 ic 2b00000000000000000000000000010000  NB. ...will be implicitly conjugated.
CONJ_TRANSPOSE=:    2 ic 2b00000000000000000000000000011000  NB. ...will be implicitly transposed and conjugated.

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. conj_t
NB. Semantic meaning: Corresponding matrix/vector operand...

NO_CONJUGATE=:      2 ic 2b00000000000000000000000000000000  NB. ...will be used as given.
CONJUGATE=:         2 ic 2b00000000000000000000000000010000  NB. ...will be implicitly conjugated.

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. side_t
NB. Semantic meaning: Corresponding matrix operand...

LEFT=:              2 ic 2b00000000000000000000000000000000  NB. ...appears on the left.
RIGHT=:             2 ic 2b00000000000000000000000000000001  NB. ...appears on the right.

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. struc_t
NB. Semantic meaning: Matrix operand...

GENERAL=:           2 ic 2b00000000000000000000000000000000  NB. ...has no structure.
HERMITIAN=:         2 ic 2b00001000000000000000000000000000  NB. ...has Hermitian structure.
SYMMETRIC=:         2 ic 2b00010000000000000000000000000000  NB. ...has symmetric structure.
TRIANGULAR=:        2 ic 2b00011000000000000000000000000000  NB. ...has triangular structure.

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. uplo_t
NB. Semantic meaning: Corresponding matrix operand...

ZEROS=:             2 ic 2b00000000000000000000000000000000  NB. ...is filled by zeros.
UPPER=:             2 ic 2b00000000000000000000000001100000  NB. ...is stored in (and will be accessed only from) the upper triangle.
LOWER=:             2 ic 2b00000000000000000000000011000000  NB. ...is stored in (and will be accessed only from) the lower triangle.
DENSE=:             2 ic 2b00000000000000000000000011100000  NB. ...is stored as a full matrix (ie: in both triangles).

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. diag_t
NB. Semantic meaning: Corresponding matrix operand...

NONUNIT_DIAG=:      2 ic 2b00000000000000000000000000000000  NB. ...has a non-unit diagonal that should be explicitly read from.
UNIT_DIAG=:         2 ic 2b00000000000000000000000100000000  NB. ...has a unit diagonal that should be implicitly assumed (and not read from).

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. obj_t->objbits_t
NB. info bit field masks

TARGET_DT_SHIFT=: 10
EXEC_DT_SHIFT=:   13
COMP_DT_SHIFT=:   29

DOMAIN_BIT=:        2 ic 2b00000000000000000000000000000001  NB. 1 <<  0 ==          1
PRECISION_BIT=:     2 ic 2b00000000000000000000000000000010  NB. 1 <<  1 ==          2
DATATYPE_BITS=:     2 ic 2b00000000000000000000000000000111  NB. 7 <<  0 ==          7
TRANS_BIT=:         2 ic 2b00000000000000000000000000001000  NB. 1 <<  3 ==          8
CONJ_BIT=:          2 ic 2b00000000000000000000000000010000  NB. 1 <<  4 ==         16
CONJTRANS_BITS=:    2 ic 2b00000000000000000000000000011000  NB. 3 <<  3 ==         24
UPPER_BIT=:         2 ic 2b00000000000000000000000000100000  NB. 1 <<  5 ==         32
DIAG_BIT=:          2 ic 2b00000000000000000000000001000000  NB. 1 <<  6 ==         64
LOWER_BIT=:         2 ic 2b00000000000000000000000010000000  NB. 1 <<  7 ==        128
UPLO_BITS=:         2 ic 2b00000000000000000000000011100000  NB. 7 <<  5 ==        224
UNIT_DIAG_BIT=:     2 ic 2b00000000000000000000000100000000  NB. 1 <<  8 ==        256
TARGET_DT_BITS=:    2 ic 2b00000000000000000001110000000000  NB. 7 << 10 ==       7168
EXEC_DT_BITS=:      2 ic 2b00000000000000001110000000000000  NB. 7 << 13 ==      57344
STRUC_BITS=:        2 ic 2b00011000000000000000000000000000  NB. 3 << 27 ==  402653184
COMP_DT_BITS=:      2 ic 2b11100000000000000000000000000000  NB. 7 << 29 == 3758096384

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. enumerated type value definitions

BITVAL_UPPER=:      UPPER_BIT or DIAG_BIT
BITVAL_LOWER=:      LOWER_BIT or DIAG_BIT
BITVAL_COMPLEX=:    DOMAIN_BIT
BITVAL_CONST_TYPE=: CONSTANT

NB. ---------------------------------------------------------
NB. Global scalar constants

MINUS_TWO=: 0 {:: dlsym LIB ; 'BLIS_MINUS_TWO'
MINUS_ONE=: 0 {:: dlsym LIB ; 'BLIS_MINUS_ONE'
ZERO=:      0 {:: dlsym LIB ; 'BLIS_ZERO'
ONE=:       0 {:: dlsym LIB ; 'BLIS_ONE'
TWO=:       0 {:: dlsym LIB ; 'BLIS_TWO'

NB. ---------------------------------------------------------
NB. Structures

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Object type

NB. typedef struct obj_s
NB. {                                       32bit            64bit
NB.         // Basic fields                 offset length    offset length
NB.         struct obj_s* root;               0      4        0       8
NB.
NB.         dim_t         off[2];             4    2*4        8     2*8
NB.         dim_t         dim[2];            12    2*4       24     2*8
NB.         doff_t        diag_off;          20      4       40       8
NB.
NB.         objbits_t     info;              24     =4       48      =4
NB.         objbits_t     info2;             28     =4       52      =4
NB.         siz_t         elem_size;         32      4       56       8
NB.
NB.         void*         buffer;            36      4       64       8
NB.         inc_t         rs;                40      4       72       8
NB.         inc_t         cs;                44      4       80       8
NB.         inc_t         is;                48      4       88       8
NB.
NB.         // Bufferless scalar storage
NB.         atom_t        scalar;            52    2*4       96     2*8
NB.
NB.         // Pack-related fields
NB.         dim_t         m_padded;          60      4      112       8
NB.         dim_t         n_padded;          64      4      120       8
NB.         inc_t         ps;                68      4      128       8
NB.         inc_t         pd;                72      4      136       8
NB.         dim_t         m_panel;           76      4      144       8
NB.         dim_t         n_panel;           80      4      152       8
NB.
NB.         // User-customizable fields
NB.         obj_pack_fn_t pack_fn;           84      4      160       8
NB.         void*         pack_params;       88      4      168       8
NB.         obj_ker_fn_t  ker_fn;            92      4      176       8
NB.         void*         ker_params;        96      4      184       8
NB.                                 sizeof= 100             192
NB. } obj_t;

SIZEOF_OBJ_T=: IF64 { 100 192

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Context type

NB. typedef struct cntx_s                                        32bit            64bit
NB. {                                                            offset length    offset length
NB.         blksz_t   blkszs[ BLIS_NUM_BLKSZS ];                                     0   22*(4+4) = 176 ints   = 1408 bytes
NB.         bszid_t   bmults[ BLIS_NUM_BLKSZS ];                                  1408               22 shorts =   88 bytes
NB.
NB.         func_t    ukrs[ BLIS_NUM_UKRS ];                                      1496       49*4 = 196 ptrs   = 1568 bytes
NB.         mbool_t   ukr_prefs[ BLIS_NUM_UKR_PREFS ];                            3064       14*4 =  56 bools  =   56 bytes
NB.
NB.         void_fp   l3_sup_handlers[ BLIS_NUM_LEVEL3_OPS ];                     3120               11 ptrs   =   88 bytes
NB.
NB.         ind_t     method;                                                     3208                1 short  =    4 bytes
NB.                                                      align                    3212                              4 bytes
NB.                                                      sizeof=                  3216
NB. } cntx_t;

SIZEOF_CNTX_T=: IF64 { 1680 3216

NB. cntx_t* bli_gks_query_cntx( void );
NB.
NB. Application:
NB. - JE: j9.5.0-beta3/j64avx2/linux
NB.       cntx_addr=. gksquerycntxcd ''
NB.       blkszs=. 22 2 4 $ memr cntx_addr , 0 , 176 , JINT
NB.       bmults_raw=. memr cntx_addr , 1408 , 88 , JCHAR
NB.       $ bmults=. _2 (3!:4) bmults_raw
NB.    22
NB.       ukrs=. 49 4 $ memr cntx_addr , 1496 , 196 , JPTR
NB.       ukr_prefs_raw=. memr cntx_addr , 3064 , 56 , JCHAR
NB.       ukr_prefs=. 14 4 $ 1 = a. i. ukr_prefs_raw
NB.       l3_sup_handlers=. memr cntx_addr , 3120 , 11 , JPTR
NB.       method_raw=. memr cntx_addr , 3208 , 4 , JCHAR
NB.       ] method=. _2 (3!:4) method_raw
NB.    1
NB.       a. i. memr cntx_addr , 3212 , 4 , JCHAR
NB.    0 0 0 0

gks_query_cntx_cd=: (lib,' bli_gks_query_cntx > ',ifw,' x')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Runtime type

NB. typedef struct rntm_s
NB. {                                                                                    32bit            64bit
NB.         // "External" fields: these may be queried by the end-user.                  offset length    offset length
NB.         timpl_t   thread_impl;                                                                         0     1 short =  4 bytes
NB.
NB.         bool      auto_factor;                                                                         4     1 bool  =  1 byte
NB.                                                                                                        5                3 bytes
NB.         dim_t     num_threads;                                                                         8     1 int   =  8 bytes
NB.         dim_t     thrloop[ BLIS_NUM_LOOPS ];                                                          16     6 ints  = 48 bytes
NB.         bool      pack_a; // enable/disable packing of left-hand matrix A.                            64     1 bool  =  1 byte
NB.         bool      pack_b; // enable/disable packing of right-hand matrix B.                           65     1 bool  =  1 byte
NB.         bool      l3_sup; // enable/disable small matrix handling in level-3 ops.                     66     1 bool  =  1 byte
NB.                                                                              align                    67                5 bytes
NB.                                                                              sizeof=                  72               72 bytes
NB. } rntm_t;

SIZEOF_RNTM_T=: IF64 { _. 72

NB. BLIS_EXPORT_BLIS void bli_rntm_init_from_global( rntm_t* rntm );
NB.
NB. Application:
NB. - JE: j9.5.0-beta3/j64avx2/linux
NB.       rntm_l=. mema SIZEOF_RNTM_T_mtbli_
NB.       EMPTY [ rntm_init_from_global_cd < < rntm_l
NB.       a. i. memr rntm_l , 0 , 72 , JCHAR
NB.    1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0
NB.       _8 ]\ a. i. memr rntm_l , 0 , 72 , JCHAR
NB.    1 0 0 0 0 0 0 0
NB.    1 0 0 0 0 0 0 0
NB.    1 0 0 0 0 0 0 0
NB.    1 0 0 0 0 0 0 0
NB.    1 0 0 0 0 0 0 0
NB.    1 0 0 0 0 0 0 0
NB.    1 0 0 0 0 0 0 0
NB.    1 0 0 0 0 0 0 0
NB.    0 0 1 0 0 0 0 0
NB.       thread_impl_raw=. memr rntm_l , 0 , 4 , JCHAR
NB.       ] thread_impl=. _2 (3!:4) thread_impl_raw
NB.    1
NB.       auto_factor_raw=. memr rntm_l , 4 , 1 , JCHAR
NB.       ] auto_factor=. 1 = a. i. auto_factor_raw
NB.    0
NB.       ] num_threads=. memr rntm_l , 8 , 1 , JINT
NB.    1
NB.       ] thrloop=. memr rntm_l , 16 , 6 , JINT
NB.    1 1 1 1 1 1
NB.       pack_a_raw=. memr rntm_l , 64 , 1 , JCHAR
NB.       ] pack_a=. 1 = a. i. pack_a_raw
NB.    0
NB.       pack_b_raw=. memr rntm_l , 65 , 1 , JCHAR
NB.       ] pack_b=. 1 = a. i. pack_b_raw
NB.    0
NB.       l3_sup_raw=. memr rntm_l , 66 , 1 , JCHAR
NB.       ] l3_sup=. 1 = a. i. l3_sup_raw
NB.    1
NB.       a. i. memr rntm_l , 67 , 5 , JCHAR
NB.    0 0 0 0 0

rntm_init_from_global_cd=: (lib,' bli_rntm_init_from_global ',ifw,' n *')&cd

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Initialization and Cleanup

NB. void bli_init( void );

NB. void bli_finalize( void );

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. datatype

NB. num_t bli_dt_proj_to_real( num_t dt );
dt_proj_to_real=: (not BITVAL_COMPLEX)&and

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Object management

NB. ---------------------------------------------------------
NB. Object creation function reference

NB. void bli_obj_create
NB.      (
NB.        num_t   dt,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        inc_t   rs,
NB.        inc_t   cs,
NB.        obj_t*  obj
NB.      );

obj_create_cd=: (lib,' bli_obj_create ',ifw,' n x x x x x *')&cd

NB. void bli_obj_free
NB.      (
NB.        obj_t*  obj
NB.      );

obj_free_cd=: (lib,' bli_obj_free ',ifw,' n *')&cd

NB. void bli_obj_create_without_buffer
NB.      (
NB.        num_t   dt,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        obj_t*  obj
NB.      );

obj_create_without_buffer_cd=: (lib,'bli_obj_create_without_buffer ',ifw,' n x x x *')&cd

NB. void bli_obj_attach_buffer
NB.      (
NB.        void*   p,
NB.        inc_t   rs,
NB.        inc_t   cs,
NB.        inc_t   is,
NB.        obj_t*  obj
NB.      );

obj_attach_buffer_cd=: (lib,'bli_obj_attach_buffer ',ifw,' n & x x x *')&cd

NB. void bli_obj_create_with_attached_buffer
NB.      (
NB.        num_t   dt,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        void*   p,
NB.        inc_t   rs,
NB.        inc_t   cs,
NB.        obj_t*  obj
NB.      );

obj_create_with_attached_buffer_cd=: (lib,' bli_obj_create_with_attached_buffer ',ifw,' n x x x & x x *')&cd

NB. void bli_obj_alloc_buffer
NB.      (
NB.        inc_t   rs,
NB.        inc_t   cs,
NB.        inc_t   is,
NB.        obj_t*  obj
NB.      );

obj_alloc_buffer_cd=: (lib,'bli_obj_alloc_buffer ',ifw,' n x x x *')&cd

NB. void bli_obj_create_1x1
NB.      (
NB.        num_t   dt,
NB.        obj_t*  obj
NB.      );

obj_create_1x1_cd=: (lib,' bli_obj_create_1x1 ',ifw,' n x *')&cd

NB. void bli_obj_create_1x1_with_attached_buffer
NB.      (
NB.        num_t   dt,
NB.        void*   p,
NB.        obj_t*  obj
NB.      );

obj_create_1x1_with_attached_buffer_cd=: (lib,' bli_obj_create_1x1_with_attached_buffer ',ifw,' n x & *')&cd

NB. void bli_obj_create_conf_to
NB.      (
NB.        obj_t*  s,
NB.        obj_t*  d
NB.      );

obj_create_conf_to_cd=: (lib,' bli_obj_create_conf_to ',ifw,' n & *')&cd

NB. void bli_obj_scalar_init_detached
NB.      (
NB.        num_t   dt,
NB.        obj_t*  obj
NB.      );

obj_scalar_init_detached_cd=: (lib,' bli_obj_scalar_init_detached ',ifw,' n x *')&cd

NB. ---------------------------------------------------------
NB. Object initialization function reference

NB. C: void bli_obj_init_full_shallow_copy_of( obj_t* a, obj_t* b );
NB. J: trash=. a obj_init_full_shallow_copy_of b
NB. note: copy from a to b
obj_init_full_shallow_copy_of=: (memw~ memr)~&(,&(0 , SIZEOF_OBJ_T , JCHAR))

NB. ---------------------------------------------------------
NB. Object accessor function reference

NB. dim_t bli_obj_row_off( obj_t* obj );
obj_row_off=:                                          memr@(,&((IF64 {  4   8) , 1 , JINT ))

NB. dim_t bli_obj_col_off( obj_t* obj );
obj_col_off=:                                          memr@(,&((IF64 {  8  16) , 1 , JINT ))

NB. dim_t bli_obj_length( obj_t* obj );
obj_length=:                                           memr@(,&((IF64 { 12  24) , 1 , JINT ))

NB. dim_t bli_obj_width( obj_t* obj );
obj_width=:                                            memr@(,&((IF64 { 16  32) , 1 , JINT ))

NB. doff_t bli_obj_diag_offset( obj_t* obj );
obj_diag_offset=:                                      memr@(,&((IF64 { 20  40) , 1 , JINT ))

NB. num_t bli_obj_dt( obj_t* obj );
obj_dt=:                            DATATYPE_BITS  and memr@(,&((IF64 { 24  48) , 4 , JCHAR))

NB. dom_t bli_obj_domain( obj_t* obj );
obj_domain=:                        DOMAIN_BIT     and memr@(,&((IF64 { 24  48) , 4 , JCHAR))

NB. prec_t bli_obj_prec( obj_t* obj );
obj_prec=:                          PRECISION_BIT  and memr@(,&((IF64 { 24  48) , 4 , JCHAR))

NB. trans_t bli_obj_conjtrans_status( obj_t* obj );
obj_conjtrans_status=:              CONJTRANS_BITS and memr@(,&((IF64 { 24  48) , 4 , JCHAR))

NB. trans_t bli_obj_onlytrans_status( obj_t* obj );
obj_onlytrans_status=:              TRANS_BIT      and memr@(,&((IF64 { 24  48) , 4 , JCHAR))

NB. conj_t bli_obj_conj_status( obj_t* obj );
obj_conj_status=:                   CONJ_BIT       and memr@(,&((IF64 { 24  48) , 4 , JCHAR))

NB. struc_t bli_obj_struc( obj_t* obj );
obj_struc=:                         STRUC_BITS     and memr@(,&((IF64 { 24  48) , 4 , JCHAR))

NB. uplo_t bli_obj_uplo( obj_t* obj );
obj_uplo=:                          UPLO_BITS      and memr@(,&((IF64 { 24  48) , 4 , JCHAR))

NB. diag_t bli_obj_diag( obj_t* obj );
obj_diag=:                          UNIT_DIAG_BIT  and memr@(,&((IF64 { 24  48) , 4 , JCHAR))

NB. num_t bli_obj_target_dt( obj_t* obj );
obj_target_dt=: TARGET_DT_SHIFT sft TARGET_DT_BITS and memr@(,&((IF64 { 24  48) , 4 , JCHAR))

NB. siz_t bli_obj_elem_size( obj_t* obj );
obj_elem_size=:                                        memr@(,&((IF64 { 32  56) , 1 , JINT ))

NB. void* bli_obj_buffer( obj_t* obj );
obj_buffer=:                                           memr@(,&((IF64 { 36  64) , 1 , JPTR ))

NB. inc_t bli_obj_row_stride( obj_t* obj );
obj_row_stride=:                                       memr@(,&((IF64 { 40  72) , 1 , JINT ))

NB. inc_t bli_obj_col_stride( obj_t* obj );
obj_col_stride=:                                       memr@(,&((IF64 { 44  80) , 1 , JINT ))

NB. inc_t bli_obj_imag_stride( obj_t* obj );
obj_imag_stride=:                                      memr@(,&((IF64 { 48  88) , 1 , JINT ))

NB. dim_t bli_obj_padded_length( obj_t* obj );
obj_padded_length=:                                    memr@(,&((IF64 { 60 112) , 1 , JINT ))

NB. dim_t bli_obj_padded_width( obj_t* obj );
obj_padded_width=:                                     memr@(,&((IF64 { 64 120) , 1 , JINT ))

NB. dim_t bli_obj_panel_length( obj_t* obj );
obj_panel_length=:                                     memr@(,&((IF64 { 76 144) , 1 , JINT ))

NB. dim_t bli_obj_panel_width( obj_t* obj );
obj_panel_width=:                                      memr@(,&((IF64 { 80 152) , 1 , JINT ))

NB. bool bli_obj_has_trans( obj_t* obj );
obj_has_trans=: TRANS_BIT -: obj_onlytrans_status

NB. dim_t bli_obj_length_after_trans( obj_t* obj );
obj_length_after_trans=: obj_length`obj_width@.obj_has_trans

NB. dim_t bli_obj_width_after_trans( obj_t* obj );
obj_width_after_trans=: obj_width`obj_length@.obj_has_trans

NB. dim_t bli_obj_vector_dim( obj_t* obj );
obj_vector_dim=: obj_length`obj_width@.(1 -: obj_length)

NB. bool bli_obj_is_1x1( obj_t* x );
obj_is_1x1=: obj_length *.&(1&=) obj_width

NB. inc_t bli_obj_vector_inc( obj_t* obj );
obj_vector_inc=: obj_row_stride`obj_col_stride@.(1 -: obj_length)`1:@.obj_is_1x1

NB. bool bli_obj_is_upper( obj_t* obj );
obj_is_upper=: BITVAL_UPPER -: obj_uplo

NB. bool bli_obj_is_lower( obj_t* obj );
obj_is_lower=: BITVAL_LOWER -: obj_uplo

NB. bool bli_obj_is_upper_or_lower( obj_t* obj );
obj_is_upper_or_lower=: obj_is_upper +. obj_is_lower

NB. bool bli_obj_is_const( obj_t* obj );
obj_is_const=: BITVAL_CONST_TYPE -: obj_dt

NB. bool bli_obj_is_complex( obj_t* obj );
obj_is_complex=: (BITVAL_COMPLEX -: obj_domain) *. -.@obj_is_const

NB. void* bli_obj_buffer_at_off( obj_t* obj );
obj_buffer_at_off=: obj_buffer + obj_elem_size * (obj_col_off , obj_row_off) mp obj_col_stride , obj_row_stride

NB. ---------------------------------------------------------
NB. Object mutator function reference

NB. C: void bli_obj_set_dt( num_t dt, obj_t* obj );
NB. J: trash=. dt obj_set_dt obj
obj_set_dt=:        ( or                       (not DATATYPE_BITS )&and) updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. C: void bli_obj_set_target_dt( num_t dt, obj_t* obj );
NB. J: trash=. dt obj_set_target_dt obj
obj_set_target_dt=: ((or TARGET_DT_SHIFT&sft)~ (not TARGET_DT_BITS)&and) updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. C: void bli_obj_set_exec_dt( num_t dt, obj_t* obj );
NB. J: trash=. dt obj_set_exec_dt obj
obj_set_exec_dt=:   ((or EXEC_DT_SHIFT  &sft)~ (not EXEC_DT_BITS  )&and) updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. C: void bli_obj_set_comp_dt( num_t dt, obj_t* obj );
NB. J: trash=. dt obj_set_comp_dt obj
obj_set_comp_dt=:   ((or COMP_DT_SHIFT  &sft)~ (not COMP_DT_BITS  )&and) updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. C: void bli_obj_set_conjtrans( trans_t trans, obj_t* obj );
NB. J: trash=. trans obj_set_conjtrans obj
obj_set_conjtrans=: ( or                       (not CONJTRANS_BITS)&and) updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. C: void bli_obj_set_onlytrans( trans_t trans, obj_t* obj );
NB. J: trash=. trans obj_set_onlytrans obj
obj_set_onlytrans=: ( or                       (not TRANS_BIT     )&and) updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. C: void bli_obj_set_conj( conj_t conj, obj_t* obj );
NB. J: trash=. conj obj_set_conj obj
obj_set_conj=:      ( or                       (not CONJ_BIT      )&and) updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. C: void bli_obj_set_struc( struc_t struc, obj_t* obj );
NB. J: trash=. struc obj_set_struc obj
obj_set_struc=:     ( or                       (not STRUC_BITS    )&and) updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. C: void bli_obj_set_uplo( uplo_t uplo, obj_t* obj );
NB. J: trash=. uplo obj_set_uplo obj
obj_set_uplo=:      ( or                       (not UPLO_BITS     )&and) updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. C: void bli_obj_set_diag( diag_t diag, obj_t* obj );
NB. J: trash=. diag obj_set_diag obj
obj_set_diag=:      ( or                       (not UNIT_DIAG_BIT )&and) updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. void bli_obj_toggle_uplo( obj_t* obj );
NB. note: invert both lower and upper bits
obj_toggle_uplo=:   (LOWER_BIT or UPPER_BIT)&xor                         updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. C: void bli_obj_apply_trans( trans_t trans, obj_t* obj );
NB. J: trash=. trans obj_apply_trans obj
obj_apply_trans=:   xor                                                  updmem ,&((IF64 { 24 48) , 4 , JCHAR)

NB. C: void bli_obj_apply_conj( conj_t conj, obj_t* obj );
NB. J: trash=. conj obj_apply_conj obj
obj_apply_conj=: obj_apply_trans

NB. C: void bli_obj_set_length( dim_t m, obj_t* obj );
NB. J: trash=. m obj_set_length obj
obj_set_length=:        memw ,&((IF64 { 12  24) , 1 , JINT)

NB. C: void bli_obj_set_width( dim_t n, obj_t* obj );
NB. J: trash=. n obj_set_width obj
obj_set_width=:         memw ,&((IF64 { 16  32) , 1 , JINT)

NB. C: void bli_obj_set_diag_offset( doff_t doff, obj_t* obj );
NB. J: trash=. doff obj_set_diag_offset obj
obj_set_diag_offset=:   memw ,&((IF64 { 20  40) , 1 , JINT)

NB. C: void bli_obj_set_buffer( void* p, obj_t* obj );
NB. J: trash=. p obj_set_buffer obj
obj_set_buffer=:        memw ,&((IF64 { 36  64) , 1 , JINT)

NB. C: void bli_obj_set_row_stride( inc_t rs, obj_t* obj );
NB. J: trash=. rs obj_set_row_stride obj
obj_set_row_stride=:    memw ,&((IF64 { 40  72) , 1 , JINT)

NB. C: void bli_obj_set_col_stride( inc_t cs, obj_t* obj );
NB. J: trash=. cs obj_set_col_stride obj
obj_set_col_stride=:    memw ,&((IF64 { 44  80) , 1 , JINT)

NB. C: void bli_obj_set_padded_length( dim_t m, obj_t* obj );
NB. J: trash=. m obj_set_padded_length obj
obj_set_padded_length=: memw ,&((IF64 { 60 112) , 1 , JINT)

NB. C: void bli_obj_set_padded_width( dim_t n, obj_t* obj );
NB. J: trash=. n obj_set_padded_width obj
obj_set_padded_width=:  memw ,&((IF64 { 64 120) , 1 , JINT)

NB. C: void bli_obj_set_panel_length( dim_t m, obj_t* obj );
NB. J: trash=. m obj_set_panel_length obj
obj_set_panel_length=:  memw ,&((IF64 { 76 144) , 1 , JINT)

NB. C: void bli_obj_set_panel_width( dim_t n, obj_t* obj );
NB. J: trash=. n obj_set_panel_width obj
obj_set_panel_width=:   memw ,&((IF64 { 80 152) , 1 , JINT)

NB. C: void bli_obj_set_off( mdim_t mdim, dim_t offset, obj_t* obj );
NB. J: trash=. (mdim , offset) obj_set_off obj
obj_set_off=: (1 { [) memw (1 , JINT) ,~ (, ((IF64 {  4  8 ,:  8 16) {~ {.))~

NB. C: void bli_obj_set_dim( mdim_t mdim, dim_t dim_val, obj_t* obj );
NB. J: trash=. (mdim , dim_val) obj_set_dim obj
obj_set_dim=: (1 { [) memw (1 , JINT) ,~ (, ((IF64 { 12 16 ,: 24 32) {~ {.))~

NB. C: void bli_obj_set_offs( dim_t offm, dim_t offn, obj_t* obj );
NB. J: trash=. (offm , offn) obj_set_offs obj
obj_set_offs=: (bli_obj_set_off"1 0~ 0 1&,.)~"1 0

NB. C: void bli_obj_set_dims( dim_t m, dim_t n, obj_t* obj );
NB. J: trash=. (m , n) obj_set_dims obj
obj_set_dims=:        obj_set_length       `obj_set_width       "0

NB. C: void bli_obj_set_padded_dims( dim_t m, dim_t n, obj_t* obj );
NB. J: trash=. (m , n) obj_set_padded_dims obj
obj_set_padded_dims=: obj_set_padded_length`obj_set_padded_width"0

NB. C: void bli_obj_set_strides( inc_t rs, inc_t cs, obj_t* obj );
NB. J: trash=. (rs , cs) obj_set_strides obj
obj_set_strides=:     obj_set_row_stride   `obj_set_col_stride  "0

NB. C: void bli_obj_set_panel_dims( dim_t m, dim_t n, obj_t* obj );
NB. J: trash=. (m , n) obj_set_panel_dims obj
obj_set_panel_dims=:  obj_set_panel_length `obj_set_panel_width "0

NB. ---------------------------------------------------------
NB. Other object function reference

NB. C: void bli_obj_induce_trans( obj_t* obj );
NB. J: trash=. obj_induce_trans obj
obj_induce_trans=: 3 : 0
  (obj_set_dims       ~ obj_width      , obj_length    ) y
  (obj_set_strides    ~ obj_col_stride , obj_row_stride) y
  (obj_set_offs       ~ obj_col_off    , obj_row_off   ) y
  (obj_set_diag_offset~ -@obj_diag_offset              ) y

  obj_toggle_uplo^:obj_is_upper_or_lower y

  NB. Induce transposition among packed fields.
  (obj_set_padded_dims~ obj_padded_width , obj_padded_length) y
  (obj_set_panel_dims ~ obj_panel_width  , obj_panel_length ) y

  NB. Note that this verb DOES NOT touch the transposition
  NB. bit! If the calling code is using this function to
  NB. handle an object whose transposition bit is set prior
  NB. to computation, that code needs to manually clear or
  NB. toggle the bit, via bli_obj_set_onlytrans() or
  NB. bli_obj_toggle_trans(), respectively.
)

NB. C: void bli_obj_alias_to( obj_t* a, obj_t* b );
NB. J: trash=. a obj_alias_to b
obj_alias_to=: obj_init_full_shallow_copy_of

NB. C: void bli_obj_real_part( obj_t* c, obj_t* r );
NB. J: trash=. c obj_real_part r
obj_real_part=: 4 : 0
  x obj_alias_to y

  if. obj_is_complex x do.
    NB. Change the datatypes.
    y obj_set_dt       ~ dt_proj_to_real obj_dt        x
    y obj_set_target_dt~ dt_proj_to_real obj_target_dt x
    y obj_set_exec_dt  ~ dt_proj_to_real obj_exec_dt   x
    y obj_set_comp_dt  ~ dt_proj_to_real obj_comp_dt   x

    NB. Don't touch the attached scalar datatype.

    NB. Update the element size.
    y obj_set_elem_size~ -: obj_elem_size x

    NB. Update the strides.
    y obj_set_strides~ +: (obj_row_stride , obj_col_stride) x

    NB. Buffer is left unchanged.
  end.
)

NB. C: void bli_obj_imag_part( obj_t* c, obj_t* i );
NB. J: trash=. c obj_imag_part i
obj_imag_part=: 4 : 0
  if. obj_is_complex x do.
    x obj_alias_to y

    NB. Change the datatype.
    y obj_set_dt       ~ dt_proj_to_real obj_dt        x
    y obj_set_target_dt~ dt_proj_to_real obj_target_dt x
    y obj_set_exec_dt  ~ dt_proj_to_real obj_exec_dt   x
    y obj_set_comp_dt  ~ dt_proj_to_real obj_comp_dt   x

    NB. Don't touch the attached scalar datatype.

    NB. Update the element size.
    es_c=. obj_elem_size x
    (-: es_c) obj_set_elem_size y

    NB. Update the strides.
    y obj_set_strides~ +: (obj_row_stride , obj_col_stride) x

    NB. Update the buffer.
    is_c=. obj_imag_stride   x
    p=.    obj_buffer_at_off x
    (p + is_c * -: es_c) obj_set_buffer y
  end.
)

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Computational function reference

NB. ---------------------------------------------------------
NB. Object API

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1v operations

NB. y := y + conj?(x)
NB.
NB. void bli_addv
NB.      (
NB.        obj_t*  x,
NB.        obj_t*  y,
NB.      );

addv_cd=: (lib,' bli_addv ',ifw,' n & *')&cd

NB. void bli_amaxv
NB.      (
NB.        obj_t*  x,
NB.        obj_t*  index
NB.      );

amaxv_cd=: (lib,' bli_amaxv ',ifw,' n & *')&cd

NB. y := y + conj?(alpha) * conj?(x)
NB.
NB. void bli_axpyv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x,
NB.        obj_t*  y
NB.      );

axpyv_cd=: (lib,' bli_axpyv ',ifw,' n & & *')&cd

NB. y := conj?(beta) * y + conj?(alpha) * conj?(x)
NB.
NB. void bli_axpbyv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x,
NB.        obj_t*  beta,
NB.        obj_t*  y
NB.      );

axpbyv_cd=: (lib,' bli_axpbyv ',ifw,' n & & & *')&cd

NB. y := conj?(x)
NB.
NB. void bli_copyv
NB.      (
NB.        obj_t*  x,
NB.        obj_t*  y
NB.      );

copyv_cd=: (lib,' bli_copyv ',ifw,' n & *')&cd

NB. rho := conj?(x)^T * conj?(y)
NB.
NB. void bli_dotv
NB.      (
NB.        obj_t*  x,
NB.        obj_t*  y,
NB.        obj_t*  rho
NB.      );

dotv_cd=: (lib,' bli_dotv ',ifw,' n & & *')&cd

NB. rho := conj?(beta) * rho + conj?(alpha) * conj?(x)^T * conj?(y)
NB.
NB. void bli_dotxv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x,
NB.        obj_t*  y,
NB.        obj_t*  beta,
NB.        obj_t*  rho
NB.      );

dotxv_cd=: (lib,' bli_dotxv ',ifw,' n & & & & *')&cd

NB. void bli_invertv
NB.      (
NB.        obj_t*  x
NB.      );

invertv_cd=: (lib,' bli_invertv ',ifw,' n *')&cd

NB. x := ( 1.0 / conj?(alpha) ) * x
NB.
NB. void bli_invscalv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x
NB.      );

invscalv_cd=: (lib,' bli_invscalv ',ifw,' n & *')&cd

NB. x := conj?(alpha) * x
NB.
NB. void bli_scalv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x
NB.      );

scalv_cd=: (lib,' bli_scalv ',ifw,' n & *')&cd

NB. y := conj?(alpha) * conj?(x)
NB.
NB. void bli_scal2v
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x,
NB.        obj_t*  y
NB.      );

scal2v_cd=: (lib,' bli_scal2v ',ifw,' n & & *')&cd

NB. x := conj?(alpha)
NB.
NB. void bli_setv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x
NB.      );

setv_cd=: (lib,' bli_setv ',ifw,' n & *')&cd

NB. real(x) := real(alpha)
NB.
NB. void bli_setrv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x
NB.      );

setrv_cd=: (lib,' bli_setrv ',ifw,' n & *')&cd

NB. imag(x) := real(alpha)
NB.
NB. void bli_setiv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x
NB.      );

setiv_cd=: (lib,' bli_setiv ',ifw,' n & *')&cd

NB. y := y - conj?(x)
NB.
NB. void bli_subv
NB.      (
NB.        obj_t*  x,
NB.        obj_t*  y
NB.      );

subv_cd=: (lib,' bli_subv ',ifw,' n & *')&cd

NB. void bli_swapv
NB.      (
NB.        obj_t*  x,
NB.        obj_t*  y
NB.      );

swapv_cd=: (lib,' bli_swapv ',ifw,' n * *')&cd

NB. y := conj?(beta) * y + conj?(x)
NB.
NB. void bli_xpbyv
NB.      (
NB.        obj_t*  x,
NB.        obj_t*  beta,
NB.        obj_t*  y
NB.      );

xpbyv_cd=: (lib,' bli_xpbyv ',ifw,' n & & *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1d operations

NB. void bli_addd
NB.      (
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

addd_cd=: (lib,' bli_addd ',ifw,' n & *')&cd

NB. void bli_axpyd
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

axpyd_cd=: (lib,' bli_axpyd ',ifw,' n & & *')&cd

NB. void bli_copyd
NB.      (
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

copyd_cd=: (lib,' bli_copyd ',ifw,' n & *')&cd

NB. void bli_invertd
NB.      (
NB.        obj_t*  a
NB.      );

invertd_cd=: (lib,' bli_invertd ',ifw,' n *')&cd

NB. void bli_invscald
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a
NB.      );

invscald_cd=: (lib,' bli_invscald ',ifw,' n & *')&cd

NB. void bli_scald
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a
NB.      );

scald_cd=: (lib,' bli_scald ',ifw,' n & *')&cd

NB. void bli_scal2d
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

scal2d_cd=: (lib,' bli_scal2d ',ifw,' n & & *')&cd

NB. void bli_setd
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a
NB.      );

setd_cd=: (lib,' bli_setd ',ifw,' n & *')&cd

NB. void bli_setid
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a
NB.      );

setid_cd=: (lib,' bli_setid ',ifw,' n & *')&cd

NB. void bli_shiftd
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a
NB.      );

shiftd_cd=: (lib,' bli_shiftd ',ifw,' n & *')&cd

NB. void bli_subd
NB.      (
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

subd_cd=: (lib,' bli_subd ',ifw,' n & *')&cd

NB. void bli_xpbyd
NB.      (
NB.        obj_t*  a,
NB.        obj_t*  beta,
NB.        obj_t*  b
NB.      );

xpbyd_cd=: (lib,' bli_xpbyd ',ifw,' n & & *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1m operations

NB. B := B + trans?(A)
NB.
NB. void bli_addm
NB.      (
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

addm_cd=: (lib,' bli_addm ',ifw,' n & *')&cd

NB. B := B + conj?(alpha) * trans?(A)
NB.
NB. void bli_axpym
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

axpym_cd=: (lib,' bli_axpym ',ifw,' n & & *')&cd

NB. B := trans?(A)
NB.
NB. void bli_copym
NB.      (
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

copym_cd=: (lib,' bli_copym ',ifw,' n & *')&cd

NB. A := ( 1.0 / conj?(alpha) ) * A
NB.
NB. void bli_invscalm
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a
NB.      );

invscalm_cd=: (lib,' bli_invscalm ',ifw,' n & *')&cd

NB. A := conj?(alpha) * A
NB.
NB. void bli_scalm
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a
NB.      );

scalm_cd=: (lib,' bli_scalm ',ifw,' n & *')&cd

NB. B := conj?(alpha) * trans?(A)
NB.
NB. void bli_scal2m
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

scal2m_cd=: (lib,' bli_scal2m ',ifw,' n & & *')&cd

NB. A := conj?(alpha)
NB.
NB. void bli_setm
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a
NB.      );

setm_cd=: (lib,' bli_setm ',ifw,' n & *')&cd

NB. real(A) := real(alpha)
NB.
NB. void bli_setrm
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a
NB.      );

setrm_cd=: (lib,' bli_setrm ',ifw,' n & *')&cd

NB. imag(A) := real(alpha)
NB.
NB. void bli_setim
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a
NB.      );

setim_cd=: (lib,' bli_setim ',ifw,' n & *')&cd

NB. B := B - trans?(A)
NB.
NB. void bli_subm
NB.      (
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

subm_cd=: (lib,' bli_subm ',ifw,' n & *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1f operations

NB. z := z + conj?(alphax) * conj?(x) + conj?(alphay) * conj?(y)
NB.
NB. void bli_axpy2v
NB.      (
NB.        obj_t*  alphax,
NB.        obj_t*  alphay,
NB.        obj_t*  x,
NB.        obj_t*  y,
NB.        obj_t*  z
NB.      );

axpy2v_cd=: (lib,' bli_axpy2v ',ifw,' n & & & & *')&cd

NB. rho := conj?(x)^T * conj?(y)
NB. z   := z + conj?(alpha) * conj?(x)
NB.
NB. void bli_dotaxpyv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x,
NB.        obj_t*  y,
NB.        obj_t*  rho,
NB.        obj_t*  z
NB.      );

dotaxpyv_cd=: (lib,' bli_dotaxpyv ',ifw,' n & & & * *')&cd

NB. y := y + alpha * conja(A) * conjx(x)
NB.
NB. void bli_axpyf
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  x,
NB.        obj_t*  y
NB.      );

axpyf_cd=: (lib,' bli_axpyf ',ifw,' n & & & *')&cd

NB. y := conj?(beta) * y + conj?(alpha) * conj?(A)^T * conj?(x)
NB.
NB. void bli_dotxf
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  x,
NB.        obj_t*  beta,
NB.        obj_t*  y
NB.      );

dotxf_cd=: (lib,' bli_dotxf ',ifw,' n & & & & *')&cd

NB. y := conj?(beta) * y + conj?(alpha) * conj?(A)^T * conj?(w)
NB. z :=               z + conj?(alpha) * conj?(A)   * conj?(x)
NB.
NB. void bli_dotxaxpyf
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  w,
NB.        obj_t*  x,
NB.        obj_t*  beta,
NB.        obj_t*  y,
NB.        obj_t*  z
NB.      );

dotxaxpyf_cd=: (lib,' bli_dotxaxpyf ',ifw,' n & & & & & * *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-2 operations

NB. y := conj?(beta) * y + conj?(alpha) * trans?(A) * conj?(x)
NB.
NB. void bli_gemv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  x,
NB.        obj_t*  beta,
NB.        obj_t*  y
NB.      );

gemv_cd=: (lib,' bli_gemv ',ifw,' n & & & & *')&cd

NB. A := A + conj?(alpha) * conj?(x) * conj?(y)^T
NB.
NB. void bli_ger
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x,
NB.        obj_t*  y,
NB.        obj_t*  a
NB.      );

ger_cd=: (lib,' bli_ger ',ifw,' n & & & *')&cd

NB. y := conj?(beta) * y + conj?(alpha) * conj?(A) * conj?(x)
NB.
NB. void bli_hemv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  x,
NB.        obj_t*  beta,
NB.        obj_t*  y
NB.      );

hemv_cd=: (lib,' bli_hemv ',ifw,' n & & & & *')&cd

NB. A := A + conj?(alpha) * conj?(x) * conj?(x)^H
NB.
NB. void bli_her
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x,
NB.        obj_t*  a
NB.      );

her_cd=: (lib,' bli_her ',ifw,' n & & *')&cd

NB. A := A + alpha * conj?(x) * conj?(y)^H + conj(alpha) * conj?(y) * conj?(x)^H
NB.
NB. void bli_her2
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x,
NB.        obj_t*  y,
NB.        obj_t*  a
NB.      );

her2_cd=: (lib,' bli_her2 ',ifw,' n & & & *')&cd

NB. y := conj?(beta) * y + conj?(alpha) * conj?(A) * conj?(x)
NB.
NB. void bli_symv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  x,
NB.        obj_t*  beta,
NB.        obj_t*  y
NB.      );

symv_cd=: (lib,' bli_symv ',ifw,' n & & & & *')&cd

NB. A := A + conj?(alpha) * conj?(x) * conj?(x)^T
NB.
NB. void bli_syr
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x,
NB.        obj_t*  a
NB.      );

syr_cd=: (lib,' bli_syr ',ifw,' n & & *')&cd

NB. A := A + alpha * conj?(x) * conj?(y)^T + conj(alpha) * conj?(y) * conj?(x)^T
NB.
NB. void bli_syr2
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  x,
NB.        obj_t*  y,
NB.        obj_t*  a
NB.      );

syr2_cd=: (lib,' bli_syr2 ',ifw,' n & & & *')&cd

NB. x := conj?(alpha) * transa(A) * x
NB.
NB. void bli_trmv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  x
NB.      );

trmv_cd=: (lib,' bli_trmv ',ifw,' n & & *')&cd

NB. Solve the linear system
NB. transa(A) * x = alpha * y
NB.
NB. void bli_trsv
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  y
NB.      );

trsv_cd=: (lib,' bli_trsv ',ifw,' n & & *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-3 operations

NB. C := beta * C + alpha * trans?(A) * trans?(B)
NB.
NB. void bli_gemm
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b,
NB.        obj_t*  beta,
NB.        obj_t*  c,
NB.      );

gemm_cd=: (lib,' bli_gemm ',ifw,' n & & & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(B)
NB.
NB. void bli_gemmt
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b,
NB.        obj_t*  beta,
NB.        obj_t*  c
NB.      );

gemmt_cd=: (lib,' bli_gemmt ',ifw,' n & & & & *')&cd

NB. C := beta * C + alpha * conj?(A) * trans?(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * trans?(B) * conj?(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_hemm
NB.      (
NB.        side_t  sidea,
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b,
NB.        obj_t*  beta,
NB.        obj_t*  c
NB.      );

hemm_cd=: (lib,' bli_hemm ',ifw,' n x & & & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(A)^H
NB.
NB. void bli_herk
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  beta,
NB.        obj_t*  c
NB.      );

herk_cd=: (lib,' bli_herk ',ifw,' n & & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(B)^H + conj(alpha) * trans?(B) * trans?(A)^H
NB.
NB. void bli_her2k
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b,
NB.        obj_t*  beta,
NB.        obj_t*  c
NB.      );

her2k_cd=: (lib,' bli_her2k ',ifw,' n & & & & *')&cd

NB. C := beta * C + alpha * conj?(A) * trans?(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * trans?(B) * conj?(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_symm
NB.      (
NB.        side_t  sidea,
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b,
NB.        obj_t*  beta,
NB.        obj_t*  c
NB.      );

symm_cd=: (lib,' bli_symm ',ifw,' n x & & & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(A)^T
NB.
NB. void bli_syrk
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  beta,
NB.        obj_t*  c
NB.      );

syrk_cd=: (lib,' bli_syrk ',ifw,' n & & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(B)^T + alpha * trans?(B) * trans?(A)^T
NB.
NB. void bli_syr2k
NB.      (
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b,
NB.        obj_t*  beta,
NB.        obj_t*  c
NB.      );

syr2k_cd=: (lib,' bli_syr2k ',ifw,' n & & & & *')&cd

NB. B := alpha * transa(A) * B  if sidea is BLIS_LEFT, or
NB. B := alpha * B * transa(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_trmm
NB.      (
NB.        side_t  sidea,
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

trmm_cd=: (lib,' bli_trmm ',ifw,' n x & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * trans?(B) * trans?(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_trmm3
NB.      (
NB.        side_t  sidea,
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b,
NB.        obj_t*  beta,
NB.        obj_t*  c
NB.      );

trmm3_cd=: (lib,' bli_trmm3 ',ifw,' n x & & & & *')&cd

NB. Solve the linear system with multiple right-hand sides
NB.   transa(A) * X = alpha * B  if sidea is BLIS_LEFT, or
NB.   X * transa(A) = alpha * B  if sidea is BLIS_RIGHT
NB.
NB. void bli_trsm
NB.      (
NB.        side_t  sidea,
NB.        obj_t*  alpha,
NB.        obj_t*  a,
NB.        obj_t*  b
NB.      );

trsm_cd=: (lib,' bli_trsm ',ifw,' n x & & *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Utility operations

NB. void bli_asumv
NB.      (
NB.        obj_t*  x,
NB.        obj_t*  asum
NB.      );

asumv_cd=: (lib,' bli_asumv ',ifw,' n & *')&cd

NB. void bli_norm[1fi]m
NB.      (
NB.        obj_t*  a,
NB.        obj_t*  norm
NB.      );

norm1m_cd=: (lib,' bli_norm1m ',ifw,' n & *')&cd
normfm_cd=: (lib,' bli_normfm ',ifw,' n & *')&cd
normim_cd=: (lib,' bli_normim ',ifw,' n & *')&cd

NB. void bli_norm[1fi]v
NB.      (
NB.        obj_t*  x,
NB.        obj_t*  norm
NB.      );

norm1v_cd=: (lib,' bli_norm1v ',ifw,' n & *')&cd
normfv_cd=: (lib,' bli_normfv ',ifw,' n & *')&cd
normiv_cd=: (lib,' bli_normiv ',ifw,' n & *')&cd

NB. void bli_mkherm
NB.      (
NB.        obj_t*  a
NB.      );

mkherm_cd=: (lib,' bli_mkherm ',ifw,' n *')&cd

NB. void bli_mksymm
NB.      (
NB.        obj_t*  a
NB.      );

mksymm_cd=: (lib,' bli_mksymm ',ifw,' n *')&cd

NB. void bli_mktrim
NB.      (
NB.        obj_t*  a
NB.      );

mktrim_cd=: (lib,' bli_mktrim ',ifw,' n *')&cd

NB. void bli_fprintv
NB.      (
NB.        FILE*   file,
NB.        char*   s1,
NB.        obj_t*  x,
NB.        char*   format,
NB.        char*   s2
NB.      );

fprintv_cd=: (lib,' bli_fprintv ',ifw,' n & &c & &c &c')&cd

NB. void bli_fprintm
NB.      (
NB.        FILE*   file,
NB.        char*   s1,
NB.        obj_t*  a,
NB.        char*   format,
NB.        char*   s2
NB.      );

fprintm_cd=: (lib,' bli_fprintm ',ifw,' n & &c & &c &c')&cd

NB. void bli_printv
NB.      (
NB.        char*   s1,
NB.        obj_t*  x,
NB.        char*   format,
NB.        char*   s2
NB.      );

printv_cd=: (lib,' bli_printv ',ifw,' n &c & &c &c')&cd

NB. void bli_printm
NB.      (
NB.        char*   s1,
NB.        obj_t*  a,
NB.        char*   format,
NB.        char*   s2
NB.      );

printm_cd=: (lib,' bli_printm ',ifw,' n &c & &c &c')&cd

NB. void bli_randv
NB.      (
NB.        obj_t*  x
NB.      );

randv_cd=: (lib,' bli_randv ',ifw,' n *')&cd

NB. void bli_randm
NB.      (
NB.        obj_t*  a
NB.      );

randm_cd=: (lib,' bli_randm ',ifw,' n *')&cd

NB. scale_new^2 * sumsq_new = x[0]^2 + x[1]^2 + ... x[m-1]^2 + scale_old^2 * sumsq_old
NB.
NB. void bli_sumsqv
NB.      (
NB.        obj_t*  x,
NB.        obj_t*  scale,
NB.        obj_t*  sumsq
NB.      );

sumsqv_cd=: (lib,' bli_sumsqv ',ifw,' n & * *')&cd

NB. void bli_getsc
NB.      (
NB.        obj_t*   chi,
NB.        double*  zeta_r,
NB.        double*  zeta_i
NB.      );

getsc_cd=: (lib,' bli_getsc ',ifw,' n & * *')&cd

NB. err_t bli_getijv
NB.       (
NB.         dim_t    i,
NB.         obj_t*   x,
NB.         double*  ar,
NB.         double*  ai
NB.       );

getijv_cd=: (lib,' bli_getijv ',ifw,' x x & *d *d')&cd

NB. err_t bli_getijm
NB.       (
NB.         dim_t    i,
NB.         dim_t    j,
NB.         obj_t*   b,
NB.         double*  ar,
NB.         double*  ai
NB.       );

getijm_cd=: (lib,' bli_getijm ',ifw,' x x x & *d *d')&cd

NB. void bli_setsc
NB.      (
NB.        double  zeta_r,
NB.        double  zeta_i,
NB.        obj_t*  chi
NB.      );

setsc_cd=: (lib,' bli_setsc ',ifw,' n d d *')&cd

NB. err_t bli_setijv
NB.      (
NB.        double  ar,
NB.        double  ai,
NB.        dim_t   i,
NB.        obj_t*  x
NB.      );

setijv_cd=: (lib,' bli_setijv ',ifw,' x d d x *')&cd

NB. err_t bli_setijm
NB.      (
NB.        double  ar,
NB.        double  ai,
NB.        dim_t   i,
NB.        dim_t   j,
NB.        obj_t*  b
NB.      );

setijm_cd=: (lib,' bli_setijm ',ifw,' x d d x x *')&cd

NB. void bli_eqsc
NB.      (
NB.        obj_t*  chi,
NB.        obj_t*  psi,
NB.        bool*   is_eq
NB.      );
NB.
NB. Application:
NB.      NB. let chi and psi be some J numeric nouns
NB.      objs=. obja L: _1 chi ; psi  NB. create BLIS objects
NB.      bsc=. {: a.                  NB. J byte scalar to receive the result of comparison
NB.        NB. it would be sufficient to initialize by boolean (bsc=. 0),
NB.        NB. but we want to ensure bsc will change its value
NB.        NB. so use out-of-boolean-domain value ASCII 255 here
NB.        NB. exploiting the fact JB01 and JINT datatypes use
NB.        NB. the same amount of bytes for scalar data
NB.      bscdat=. symdat < 'bsc'      NB. pointer to scalar's value
NB.      EMPTY [ eqsc_cd_mtbli_ <"0 objs , < bscdat
NB.        NB. call bli_eqsc with arguments - pointers to data
NB.      a. i. bsc                    NB. cast the result to integer
NB.   1                               NB. value is changed to some value within boolean domain
NB.      objf L: 0 objs               NB. destroy allocated objects
NB.      ({: a.) memw bscdat , 0 1    NB. reset to prepare to re-use in another comparison, optionally

eqsc_cd=: (lib,' bli_eqsc ',ifw,' n & & *b')&cd

NB. void bli_eqv
NB.      (
NB.        obj_t*  x,
NB.        obj_t*  y,
NB.        bool*   is_eq
NB.      );

eqv_cd=: (lib,' bli_eqv ',ifw,' n & & *b')&cd

NB. void bli_eqm
NB.      (
NB.        obj_t*  a,
NB.        obj_t*  b,
NB.        bool*   is_eq
NB.      );

eqm_cd=: (lib,' bli_eqm ',ifw,' n & & *b')&cd

NB. ---------------------------------------------------------
NB. Typed API

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1v operations

NB. y := y + conjx(x)
NB.
NB. void bli_?addv
NB.      (
NB.        conj_t  conjx,
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy
NB.      );

daddv_cd=: (lib,' bli_daddv ',ifw,' n x x &d x *d')&cd
zaddv_cd=: (lib,' bli_zaddv ',ifw,' n x x &j x *j')&cd

NB. mimic BLAS routines i?amax()
NB.
NB. void bli_?amaxv
NB.      (
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        dim_t*  index
NB.      );

damaxv_cd=: (lib,' bli_damaxv ',ifw,' n x &d x *x')&cd
zamaxv_cd=: (lib,' bli_zamaxv ',ifw,' n x &j x *x')&cd

NB. y := y + alpha * conjx(x)
NB.
NB. void bli_?axpyv
NB.      (
NB.        conj_t  conjx,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy
NB.      );

daxpyv_cd=: (lib,' bli_daxpyv ',ifw,' n x x &d &d x *d x')&cd
zaxpyv_cd=: (lib,' bli_zaxpyv ',ifw,' n x x &j &j x *j x')&cd

NB. y := beta * y + alpha * conjx(x)
NB.
NB. void bli_?axpbyv
NB.      (
NB.        conj_t  conjx,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  beta,
NB.        ctype*  y, inc_t incy
NB.      );

daxpbyv_cd=: (lib,' bli_daxpbyv ',ifw,' n x x &d &d x &d *d x')&cd
zaxpbyv_cd=: (lib,' bli_zaxpbyv ',ifw,' n x x &j &j x &j *j x')&cd

NB. y := conjx(x)
NB.
NB. void bli_?copyv
NB.      (
NB.        conj_t  conjx,
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy
NB.      );

dcopyv_cd=: (lib,' bli_dcopyv ',ifw,' n x x &d x *d x')&cd
zcopyv_cd=: (lib,' bli_zcopyv ',ifw,' n x x &j x *j x')&cd

NB. rho := conjx(x)^T * conjy(y)
NB.
NB. void bli_?dotv
NB.      (
NB.        conj_t  conjx,
NB.        conj_t  conjy,
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy,
NB.        ctype*  rho
NB.      );

ddotv_cd=: (lib,' bli_ddotv ',ifw,' n x x x &d x &d x *d')&cd
zdotv_cd=: (lib,' bli_zdotv ',ifw,' n x x x &j x &j x *j')&cd

NB. rho := beta * rho + alpha * conjx(x)^T * conjy(y)
NB.
NB. void bli_?dotxv
NB.      (
NB.        conj_t  conjx,
NB.        conj_t  conjy,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy,
NB.        ctype*  beta,
NB.        ctype*  rho
NB.      );

ddotxv_cd=: (lib,' bli_ddotxv ',ifw,' n x x x &d &d x &d x &d *d')&cd
zdotxv_cd=: (lib,' bli_zdotxv ',ifw,' n x x x &j &j x &j x &j *j')&cd

NB. Invert all elements
NB.
NB. void bli_?invertv
NB.      (
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx
NB.      );

dinvertv_cd=: (lib,' bli_dinvertv ',ifw,' n x *d x')&cd
zinvertv_cd=: (lib,' bli_zinvertv ',ifw,' n x *j x')&cd

NB. x := ( 1.0 / conjalpha(alpha) ) * x
NB.
NB. void bli_?invscalv
NB.      (
NB.        conj_t  conjalpha,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx
NB.      );

dinvscalv_cd=: (lib,' bli_dinvscalv ',ifw,' n x x &d *d x')&cd
zinvscalv_cd=: (lib,' bli_zinvscalv ',ifw,' n x x &j *j x')&cd

NB. x := conjalpha(alpha) * x
NB.
NB. void bli_?scalv
NB.      (
NB.        conj_t  conjalpha,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx
NB.      );

dscalv_cd=: (lib,' bli_dscalv ',ifw,' n x x &d *d x')&cd
zscalv_cd=: (lib,' bli_zscalv ',ifw,' n x x &j *j x')&cd

NB. y := alpha * conjx(x)
NB.
NB. void bli_?scal2v
NB.      (
NB.        conj_t  conjx,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy
NB.      );

dscal2v_cd=: (lib,' bli_dscal2v ',ifw,' n x x &d &d x *d x')&cd
zscal2v_cd=: (lib,' bli_zscal2v ',ifw,' n x x &j &j x *j x')&cd

NB. x := conjalpha(alpha)
NB.
NB. void bli_?setv
NB.      (
NB.        conj_t  conjalpha,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx
NB.      );

dsetv_cd=: (lib,' bli_dsetv ',ifw,' n x x &d *d x')&cd
zsetv_cd=: (lib,' bli_zsetv ',ifw,' n x x &d *d x')&cd

NB. y := y - conjx(x)
NB.
NB. void bli_?subv
NB.      (
NB.        conj_t  conjx,
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy
NB.      );

dsubv_cd=: (lib,' bli_dsubv ',ifw,' n x x &d x *d x')&cd
zsubv_cd=: (lib,' bli_zsubv ',ifw,' n x x &j x *j x')&cd

NB. Swap corresponding elements of two n-length vectors x and y
NB.
NB. void bli_?swapv
NB.      (
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy
NB.      );

dswapv_cd=: (lib,' bli_dswapv ',ifw,' n x *d x *d x')&cd
zswapv_cd=: (lib,' bli_zswapv ',ifw,' n x *j x *j x')&cd

NB. y := beta * y + conjx(x)
NB.
NB. void bli_?xpbyv
NB.      (
NB.        conj_t  conjx,
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  beta,
NB.        ctype*  y, inc_t incy
NB.      );

dxpbyv_cd=: (lib,' bli_dxpbyv ',ifw,' n x x &d x &d *d x')&cd
zxpbyv_cd=: (lib,' bli_zxpbyv ',ifw,' n x x &j x &j *j x')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1d operations

NB. void bli_?addd
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

daddd_cd=: (lib,' bli_daddd ',ifw,' n x x x x x &d x x *d x x')&cd
zaddd_cd=: (lib,' bli_zaddd ',ifw,' n x x x x x &j x x *j x x')&cd

NB. void bli_?axpyd
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

daxpyd_cd=: (lib,' bli_daxpyd ',ifw,' n x x x x x &d &d x x *d x x')&cd
zaxpyd_cd=: (lib,' bli_zaxpyd ',ifw,' n x x x x x &j &j x x *j x x')&cd

NB. void bli_?copyd
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

dcopyd_cd=: (lib,' bli_dcopyd ',ifw,' n x x x x x &d x x *d x x')&cd
zcopyd_cd=: (lib,' bli_zcopyd ',ifw,' n x x x x x &j x x *j x x')&cd

NB. void bli_?invertd
NB.      (
NB.        doff_t  diagoffa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dinvertd_cd=: (lib,' bli_dinvertd ',ifw,' n x x x *d x x')&cd
zinvertd_cd=: (lib,' bli_zinvertd ',ifw,' n x x x *j x x')&cd

NB. void bli_?invscald
NB.      (
NB.        conj_t  conjalpha,
NB.        doff_t  diagoffa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dinvscald_cd=: (lib,' bli_dinvscald ',ifw,' n x x x x &d *d x x')&cd
zinvscald_cd=: (lib,' bli_zinvscald ',ifw,' n x x x x &j *j x x')&cd

NB. void bli_?scald
NB.      (
NB.        conj_t  conjalpha,
NB.        doff_t  diagoffa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dscald_cd=: (lib,' bli_dscald ',ifw,' n x x x x &d *d x x')&cd
zscald_cd=: (lib,' bli_zscald ',ifw,' n x x x x &j *j x x')&cd

NB. void bli_?scal2d
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

dscal2d_cd=: (lib,' bli_dscal2d ',ifw,' n x x x x x &d &d x x *d x x')&cd
zscal2d_cd=: (lib,' bli_zscal2d ',ifw,' n x x x x x &j &j x x *j x x')&cd

NB. void bli_?setd
NB.      (
NB.        conj_t  conjalpha,
NB.        doff_t  diagoffa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dsetd_cd=: (lib,' bli_dsetd ',ifw,' n x x x x &d *d x x')&cd
zsetd_cd=: (lib,' bli_zsetd ',ifw,' n x x x x &j *j x x')&cd

NB. void bli_?setid
NB.      (
NB.        doff_t    diagoffa,
NB.        dim_t     m,
NB.        dim_t     n,
NB.        ctype_r*  alpha,
NB.        ctype*    a, inc_t rsa, inc_t csa
NB.      );

dsetid_cd=: (lib,' bli_dsetid ',ifw,' n x x x &d *d x x')&cd
zsetid_cd=: (lib,' bli_zsetid ',ifw,' n x x x &j *j x x')&cd

NB. void bli_?shiftd
NB.      (
NB.        doff_t  diagoffa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dshiftd_cd=: (lib,' bli_dshiftd ',ifw,' n x x x &d *d x x')&cd
zshiftd_cd=: (lib,' bli_zshiftd ',ifw,' n x x x &j *j x x')&cd

NB. void bli_?subd
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

dsubd_cd=: (lib,' bli_dsubd ',ifw,' n x x x x x &d x x *d x x')&cd
zsubd_cd=: (lib,' bli_zsubd ',ifw,' n x x x x x &j x x *j x x')&cd

NB. void bli_?xpbyd
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   beta,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

dxpbyd_cd=: (lib,' bli_dxpbyd ',ifw,' n x x x x x &d x x &d *d x x')&cd
zxpbyd_cd=: (lib,' bli_zxpbyd ',ifw,' n x x x x x &j x x &j *j x x')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1m operations

NB. B := B + transa(A)
NB.
NB. void bli_?addm
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        uplo_t   uploa,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

daddm_cd=: (lib,' bli_daddm ',ifw,' n x x x x x x &d x x *d x x')&cd
zaddm_cd=: (lib,' bli_zaddm ',ifw,' n x x x x x x &j x x *j x x')&cd

NB. B := B + alpha * transa(A)
NB.
NB. void bli_?axpym
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        uplo_t   uploa,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

daxpym_cd=: (lib,' bli_daxpym ',ifw,' n x x x x x x &d &d x x *d x x')&cd
zaxpym_cd=: (lib,' bli_zaxpym ',ifw,' n x x x x x x &j &j x x *j x x')&cd

NB. B := transa(A)
NB.
NB. void bli_?copym
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        uplo_t   uploa,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

dcopym_cd=: (lib,' bli_dcopym ',ifw,' n x x x x x x &d x x *d x x')&cd
zcopym_cd=: (lib,' bli_zcopym ',ifw,' n x x x x x x &j x x *j x x')&cd

NB. A := ( 1.0 / conjalpha(alpha) ) * A
NB.
NB. void bli_?invscalm
NB.      (
NB.        conj_t  conjalpha,
NB.        doff_t  diagoffa,
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dinvscalm_cd=: (lib,' bli_dinvscalm ',ifw,' n x x x x x &d *d x x')&cd
zinvscalm_cd=: (lib,' bli_zinvscalm ',ifw,' n x x x x x &j *j x x')&cd

NB. A := conjalpha(alpha) * A
NB.
NB. void bli_?scalm
NB.      (
NB.        conj_t  conjalpha,
NB.        doff_t  diagoffa,
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dscalm_cd=: (lib,' bli_dscalm ',ifw,' n x x x x x &d *d x x')&cd
zscalm_cd=: (lib,' bli_zscalm ',ifw,' n x x x x x &j *j x x')&cd

NB. B := alpha * transa(A)
NB.
NB. void bli_?scal2m
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        uplo_t   uploa,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

dscal2m_cd=: (lib,' bli_dscal2m ',ifw,' n x x x x x x &d &d x x *d x x')&cd
zscal2m_cd=: (lib,' bli_zscal2m ',ifw,' n x x x x x x &j &j x x *j x x')&cd

NB. Set all elements of an m x n matrix A to conjalpha(alpha)
NB.
NB. void bli_?setm
NB.      (
NB.        conj_t  conjalpha,
NB.        doff_t  diagoffa,
NB.        diag_t  diaga,
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dsetm_cd=: (lib,' bli_dsetm ',ifw,' n x x x x x x &d *d x x')&cd
zsetm_cd=: (lib,' bli_zsetm ',ifw,' n x x x x x x &j *j x x')&cd

NB. B := B - transa(A)
NB.
NB. void bli_?subm
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        uplo_t   uploa,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

dsubm_cd=: (lib,' bli_dsubm ',ifw,' n x x x x x x &d x x *d x x')&cd
zsubm_cd=: (lib,' bli_zsubm ',ifw,' n x x x x x x &j x x *j x x')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1f operations

NB. z := z + alphax * conjx(x) + alphay * conjy(y)
NB.
NB. void bli_?axpy2v
NB.      (
NB.        conj_t  conjx,
NB.        conj_t  conjy,
NB.        dim_t   m,
NB.        ctype*  alphax,
NB.        ctype*  alphay,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy,
NB.        ctype*  z, inc_t incz
NB.      );

daxpy2v_cd=: (lib,' bli_daxpy2v ',ifw,' n x x x &d &d &d x &d x *d x')&cd
zaxpy2v_cd=: (lib,' bli_zaxpy2v ',ifw,' n x x x &j &d &j x &j x *j x')&cd

NB. rho := conjxt(x^T) * conjy(y)
NB. z   := z + alpha * conjx(x)
NB.
NB. void bli_?dotaxpyv
NB.      (
NB.        conj_t  conjxt,
NB.        conj_t  conjx,
NB.        conj_t  conjy,
NB.        dim_t   m,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy,
NB.        ctype*  rho,
NB.        ctype*  z, inc_t incz
NB.      );

ddotaxpyv_cd=: (lib,' bli_ddotaxpyv ',ifw,' n x x x x &d &d x &d x *d *d x')&cd
zdotaxpyv_cd=: (lib,' bli_zdotaxpyv ',ifw,' n x x x x &j &j x &j x *j *j x')&cd

NB. y := y + alpha * conja(A) * conjx(x)
NB.
NB. void bli_?axpyf
NB.      (
NB.        conj_t  conja,
NB.        conj_t  conjx,
NB.        dim_t   m,
NB.        dim_t   b,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t inca, inc_t lda,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy
NB.      );

daxpyf_cd=: (lib,' bli_daxpyf ',ifw,' n x x x x &d &d x x &d x *d x')&cd
zaxpyf_cd=: (lib,' bli_zaxpyf ',ifw,' n x x x x &j &j x x &j x *j x')&cd

NB. y := beta * y + alpha * conjat(A^T) * conjx(x)
NB.
NB. void bli_?dotxf
NB.      (
NB.        conj_t  conjat,
NB.        conj_t  conjx,
NB.        dim_t   m,
NB.        dim_t   b,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t inca, inc_t lda,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  beta,
NB.        ctype*  y, inc_t incy
NB.      );

ddotxf_cd=: (lib,' bli_ddotxf ',ifw,' n x x x x &d &d x x &d x &d *d x')&cd
zdotxf_cd=: (lib,' bli_zdotxf ',ifw,' n x x x x &j &j x x &j x &j *d x')&cd

NB. y := beta * y + alpha * conjat(A^T) * conjw(w)
NB. z :=        z + alpha * conja(A)    * conjx(x)
NB.
NB. void bli_?dotxaxpyf
NB.      (
NB.        conj_t  conjat,
NB.        conj_t  conja,
NB.        conj_t  conjw,
NB.        conj_t  conjx,
NB.        dim_t   m,
NB.        dim_t   b,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t inca, inc_t lda,
NB.        ctype*  w, inc_t incw,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  beta,
NB.        ctype*  y, inc_t incy,
NB.        ctype*  z, inc_t incz
NB.      );

ddotxaxpyf_cd=: (lib,' bli_ddotxaxpyf ',ifw,' n x x x x x x &d &d x x &d x &d x &d *d x *d x')&cd
zdotxaxpyf_cd=: (lib,' bli_zdotxaxpyf ',ifw,' n x x x x x x &j &j x x &j x &j x &j *j x *j x')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-2 operations

NB. y := beta * y + alpha * transa(A) * conjx(x)
NB.
NB. void bli_?gemv
NB.      (
NB.        trans_t  transa,
NB.        conj_t   conjx,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   x, inc_t incx,
NB.        ctype*   beta,
NB.        ctype*   y, inc_t incy
NB.      );

dgemv_cd=: (lib,' bli_dgemv ',ifw,' n x x x x &d &d x x &d x &d *d x')&cd
zgemv_cd=: (lib,' bli_zgemv ',ifw,' n x x x x &j &j x x &j x &j *j x')&cd

NB. A := A + alpha * conjx(x) * conjy(y)^T
NB.
NB. void bli_?ger
NB.      (
NB.        conj_t  conjx,
NB.        conj_t  conjy,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dger_cd=: (lib,' bli_dger ',ifw,' n x x x x &d &d x &d x *d x x')&cd
zger_cd=: (lib,' bli_zger ',ifw,' n x x x x &j &j x &j x *j x x')&cd

NB. y := beta * y + alpha * conja(A) * conjx(x)
NB.
NB. void bli_?hemv
NB.      (
NB.        uplo_t  uploa,
NB.        conj_t  conja,
NB.        conj_t  conjx,
NB.        dim_t   m,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t rsa, inc_t csa,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  beta,
NB.        ctype*  y, inc_t incy
NB.      );

dhemv_cd=: (lib,' bli_dhemv ',ifw,' n x x x x &d &d x x &d x &d *d x')&cd
zhemv_cd=: (lib,' bli_zhemv ',ifw,' n x x x x &j &j x x &j x &j *j x')&cd

NB. A := A + alpha * conjx(x) * conjx(x)^H
NB.
NB. void bli_?her
NB.      (
NB.        uplo_t  uploa,
NB.        conj_t  conjx,
NB.        dim_t   m,
NB.        rtype*  alpha,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dher_cd=: (lib,' bli_dher ',ifw,' n x x x &d &d x *d x x')&cd
zher_cd=: (lib,' bli_zher ',ifw,' n x x x &j &j x *j x x')&cd

NB. A := A + alpha * conjx(x) * conjy(y)^H + conj(alpha) * conjy(y) * conjx(x)^H
NB.
NB. void bli_?her2
NB.      (
NB.        uplo_t  uploa,
NB.        conj_t  conjx,
NB.        conj_t  conjy,
NB.        dim_t   m,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dher2_cd=: (lib,' bli_dher2 ',ifw,' n x x x x &d &d x &d x *d x x')&cd
zher2_cd=: (lib,' bli_zher2 ',ifw,' n x x x x &j &j x &j x *j x x')&cd

NB. y := beta * y + alpha * conja(A) * conjx(x)
NB.
NB. void bli_?symv
NB.      (
NB.        uplo_t  uploa,
NB.        conj_t  conja,
NB.        conj_t  conjx,
NB.        dim_t   m,
NB.        ctype*  alpha,
NB.        ctype*  a, inc_t rsa, inc_t csa,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  beta,
NB.        ctype*  y, inc_t incy
NB.      );

dsymv_cd=: (lib,' bli_dsymv ',ifw,' n x x x x &d &d x x &d x &d *d x')&cd
zsymv_cd=: (lib,' bli_zsymv ',ifw,' n x x x x &j &j x x &j x &j *j x')&cd

NB. A := A + alpha * conjx(x) * conjx(x)^T
NB.
NB. void bli_?syr
NB.      (
NB.        uplo_t  uploa,
NB.        conj_t  conjx,
NB.        dim_t   m,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dsyr_cd=: (lib,' bli_dsyr ',ifw,' n x x x &d &d x *d x x')&cd
zsyr_cd=: (lib,' bli_zsyr ',ifw,' n x x x &j &j x *j x x')&cd

NB. A := A + alpha * conjx(x) * conjy(y)^T + conj(alpha) * conjy(y) * conjx(x)^T
NB.
NB. void bli_?syr2
NB.      (
NB.        uplo_t  uploa,
NB.        conj_t  conjx,
NB.        conj_t  conjy,
NB.        dim_t   m,
NB.        ctype*  alpha,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dsyr2_cd=: (lib,' bli_dsyr2 ',ifw,' n x x x x &d &d x &d x *d x x')&cd
zsyr2_cd=: (lib,' bli_zsyr2 ',ifw,' n x x x x &j &j x &j x *j x x')&cd

NB. x := alpha * transa(A) * x
NB.
NB. void bli_?trmv
NB.      (
NB.        uplo_t   uploa,
NB.        trans_t  transa,
NB.        diag_t   diaga,
NB.        dim_t    m,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   x, inc_t incx
NB.      );

dtrmv_cd=: (lib,' bli_dtrmv ',ifw,' n x x x x &d &d x x *d x')&cd
ztrmv_cd=: (lib,' bli_ztrmv ',ifw,' n x x x x &j &j x x *j x')&cd

NB. transa(A) * x = alpha * y
NB.
NB. void bli_?trsv
NB.      (
NB.        uplo_t   uploa,
NB.        trans_t  transa,
NB.        diag_t   diaga,
NB.        dim_t    m,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   y, inc_t incy
NB.      );

dtrsv_cd=: (lib,' bli_dtrsv ',ifw,' n x x x x &d &d x x *d x')&cd
ztrsv_cd=: (lib,' bli_ztrsv ',ifw,' n x x x x &j &j x x *j x')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-3 operations

NB. C := beta * C + alpha * transa(A) * transb(B)
NB.
NB. void bli_?gemm
NB.      (
NB.        trans_t  transa,
NB.        trans_t  transb,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        dim_t    k,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb,
NB.        ctype*   beta,
NB.        ctype*   c, inc_t rsc, inc_t csc
NB.      );

dgemm_cd=: (lib,' bli_dgemm ',ifw,' n x x x x x &d &d x x &d x x &d *d x x')&cd
zgemm_cd=: (lib,' bli_zgemm ',ifw,' n x x x x x &j &j x x &j x x &j *j x x')&cd

NB. void bli_?gemm_ex
NB.      (
NB.        trans_t  transa,
NB.        trans_t  transb,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        dim_t    k,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb,
NB.        ctype*   beta,
NB.        ctype*   c, inc_t rsc, inc_t csc,
NB.        cntx_t*  cntx,
NB.        rntm_t*  rntm
NB.      );

dgemm_ex_cd=: (lib,' bli_dgemm_ex ',ifw,' n x x x x x &d &d x x &d x x &d *d x x & &')&cd
zgemm_ex_cd=: (lib,' bli_zgemm_ex ',ifw,' n x x x x x &j &j x x &j x x &j *j x x & &')&cd

NB. C := beta * C + alpha * transa(A) * transb(B)
NB.
NB. void bli_?gemmt
NB.      (
NB.        uplo_t   uploc,
NB.        trans_t  transa,
NB.        trans_t  transb,
NB.        dim_t    m,
NB.        dim_t    k,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb,
NB.        ctype*   beta,
NB.        ctype*   c, inc_t rsc, inc_t csc
NB.      );

dgemmt_cd=: (lib,' bli_dgemmt ',ifw,' n x x x x x &d &d x x &d x x &d *d x x')&cd
zgemmt_cd=: (lib,' bli_zgemmt ',ifw,' n x x x x x &j &j x x &j x x &j *j x x')&cd

NB. C := beta * C + alpha * conja(A) * transb(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * transb(B) * conja(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_?hemm
NB.      (
NB.        side_t   sidea,
NB.        uplo_t   uploa,
NB.        conj_t   conja,
NB.        trans_t  transb,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb,
NB.        ctype*   beta,
NB.        ctype*   c, inc_t rsc, inc_t csc
NB.      );

zhemm_cd=: (lib,' bli_zhemm ',ifw,' n x x x x x x &j &j x x &j x x &j *j x x')&cd

NB. C := beta * C + alpha * transa(A) * transa(A)^H
NB.
NB. void bli_?herk
NB.      (
NB.        uplo_t   uploc,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    k,
NB.        rtype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        rtype*   beta,
NB.        ctype*   c, inc_t rsc, inc_t csc
NB.      );

zherk_cd=: (lib,' bli_zherk ',ifw,' n x x x x &j &j x x &j *j x x')&cd

NB. C := beta * C + alpha * transa(A) * transb(B)^H + conj(alpha) * transb(B) * transa(A)^H
NB.
NB. void bli_?her2k
NB.      (
NB.        uplo_t   uploc,
NB.        trans_t  transa,
NB.        trans_t  transb,
NB.        dim_t    m,
NB.        dim_t    k,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb,
NB.        rtype*   beta,
NB.        ctype*   c, inc_t rsc, inc_t csc
NB.      );

zher2k_cd=: (lib,' bli_zher2k ',ifw,' n x x x x x &j &j x x &j x x &d *j x x')&cd

NB. C := beta * C + alpha * conja(A) * transb(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * transb(B) * conja(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_?symm
NB.      (
NB.        side_t   sidea,
NB.        uplo_t   uploa,
NB.        conj_t   conja,
NB.        trans_t  transb,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb,
NB.        ctype*   beta,
NB.        ctype*   c, inc_t rsc, inc_t csc
NB.      );

dsymm_cd=: (lib,' bli_dsymm ',ifw,' n x x x x x x &d &d x x &d x x &d *d x x')&cd
zsymm_cd=: (lib,' bli_zsymm ',ifw,' n x x x x x x &j &j x x &j x x &j *j x x')&cd

NB. C := beta * C + alpha * transa(A) * transa(A)^T
NB.
NB. void bli_?syrk
NB.      (
NB.        uplo_t   uploc,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    k,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   beta,
NB.        ctype*   c, inc_t rsc, inc_t csc
NB.      );

dsyrk_cd=: (lib,' bli_dsyrk ',ifw,' n x x x x &d &d x x &d *d x x')&cd
zsyrk_cd=: (lib,' bli_zsyrk ',ifw,' n x x x x &j &j x x &j *j x x')&cd

NB. C := beta * C + alpha * transa(A) * transb(B)^T + alpha * transb(B) * transa(A)^T
NB.
NB. void bli_?syr2k
NB.      (
NB.        uplo_t   uploc,
NB.        trans_t  transa,
NB.        trans_t  transb,
NB.        dim_t    m,
NB.        dim_t    k,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb,
NB.        ctype*   beta,
NB.        ctype*   c, inc_t rsc, inc_t csc
NB.      );

dsyr2k_cd=: (lib,' bli_dsyr2k ',ifw,' n x x x x x &d &d x x &d x x &d *d x x')&cd
zsyr2k_cd=: (lib,' bli_zsyr2k ',ifw,' n x x x x x &j &j x x &j x x &j *j x x')&cd

NB. B := alpha * transa(A) * B  if sidea is BLIS_LEFT, or
NB. B := alpha * B * transa(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_?trmm
NB.      (
NB.        side_t   sidea,
NB.        uplo_t   uploa,
NB.        trans_t  transa,
NB.        diag_t   diaga,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

dtrmm_cd=: (lib,' bli_dtrmm ',ifw,' n x x x x x x &d &d x x *d x x')&cd
ztrmm_cd=: (lib,' bli_ztrmm ',ifw,' n x x x x x x &j &j x x *j x x')&cd

NB. C := beta * C + alpha * transa(A) * transb(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * transb(B) * transa(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_?trmm3
NB.      (
NB.        side_t   sidea,
NB.        uplo_t   uploa,
NB.        trans_t  transa,
NB.        diag_t   diaga,
NB.        trans_t  transb,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb,
NB.        ctype*   beta,
NB.        ctype*   c, inc_t rsc, inc_t csc
NB.      );

dtrmm3_cd=: (lib,' bli_dtrmm3 ',ifw,' n x x x x x x x &d &d x x &d x x &d *d x x')&cd
ztrmm3_cd=: (lib,' bli_ztrmm3 ',ifw,' n x x x x x x x &j &j x x &j x x &j *j x x')&cd

NB. Solve the linear system with multiple right-hand sides
NB.   transa(A) * X = alpha * B  if sidea is BLIS_LEFT, or
NB.   X * transa(A) = alpha * B  if sidea is BLIS_RIGHT
NB.
NB. void bli_?trsm
NB.      (
NB.        side_t   sidea,
NB.        uplo_t   uploa,
NB.        trans_t  transa,
NB.        diag_t   diaga,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   alpha,
NB.        ctype*   a, inc_t rsa, inc_t csa,
NB.        ctype*   b, inc_t rsb, inc_t csb
NB.      );

dtrsm_cd=: (lib,' bli_dtrsm ',ifw,' n x x x x x x &d &d x x *d x x')&cd
ztrsm_cd=: (lib,' bli_ztrsm ',ifw,' n x x x x x x &j &j x x *j x x')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Utility operations

NB. void bli_?asumv
NB.      (
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        rtype*  asum
NB.      );

dasumv_cd=: (lib,' bli_dasumv ',ifw,' n x &d x *d')&cd
zasumv_cd=: (lib,' bli_zasumv ',ifw,' n x &j x *d')&cd

NB. void bli_?norm[1fi]m
NB.      (
NB.        doff_t  diagoffa,
NB.        doff_t  diaga,
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  a, inc_t rs_a, inc_t cs_a,
NB.        rtype*  norm
NB.      );

dnorm1m_cd=: (lib,' bli_dnorm1m ',ifw,' n x x x x x &d x x *d')&cd
dnormfm_cd=: (lib,' bli_dnormfm ',ifw,' n x x x x x &d x x *d')&cd
dnormim_cd=: (lib,' bli_dnormim ',ifw,' n x x x x x &d x x *d')&cd

znorm1m_cd=: (lib,' bli_znorm1m ',ifw,' n x x x x x &j x x *d')&cd
znormfm_cd=: (lib,' bli_znormfm ',ifw,' n x x x x x &j x x *d')&cd
znormim_cd=: (lib,' bli_znormim ',ifw,' n x x x x x &j x x *d')&cd

NB. void bli_?norm[1fi]v
NB.      (
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        rtype*  norm
NB.      );

dnorm1v_cd=: (lib,' bli_dnorm1v ',ifw,' n x &d x *d')&cd
dnormfv_cd=: (lib,' bli_dnormfv ',ifw,' n x &d x *d')&cd
dnormiv_cd=: (lib,' bli_dnormiv ',ifw,' n x &d x *d')&cd

znorm1v_cd=: (lib,' bli_znorm1v ',ifw,' n x &j x *d')&cd
znormfv_cd=: (lib,' bli_znormfv ',ifw,' n x &j x *d')&cd
znormiv_cd=: (lib,' bli_znormiv ',ifw,' n x &j x *d')&cd

NB. void bli_?mkherm
NB.      (
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        ctype*  a, inc_t rs_a, inc_t cs_a
NB.      );

zmkherm_cd=: (lib,' bli_zmkherm ',ifw,' n x x *j x x')&cd

NB. void bli_?mksymm
NB.      (
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        ctype*  a, inc_t rs_a, inc_t cs_a
NB.      );

dmksymm_cd=: (lib,' bli_dmksymm ',ifw,' n x x *d x x')&cd
zmksymm_cd=: (lib,' bli_zmksymm ',ifw,' n x x *j x x')&cd

NB. void bli_?mktrim
NB.      (
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        ctype*  a, inc_t rs_a, inc_t cs_a
NB.      );

dmktrim_cd=: (lib,' bli_dmktrim ',ifw,' n x x *d x x')&cd
zmktrim_cd=: (lib,' bli_zmktrim ',ifw,' n x x *j x x')&cd

NB. void bli_?fprintv
NB.      (
NB.        FILE*   file,
NB.        char*   s1,
NB.        dim_t   m,
NB.        ctype*  x, inc_t incx,
NB.        char*   format,
NB.        char*   s2
NB.      );

dfprintv_cd=: (lib,' bli_dfprintv ',ifw,' n & &c x &d x &c &c')&cd
zfprintv_cd=: (lib,' bli_zfprintv ',ifw,' n & &c x &j x &c &c')&cd

NB. void bli_?fprintm
NB.      (
NB.        FILE*   file,
NB.        char*   s1,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  a, inc_t rs_a, inc_t cs_a,
NB.        char*   format,
NB.        char*   s2
NB.      );

dfprintm_cd=: (lib,' bli_dfprintm ',ifw,' n & &c x x &d x x &c &c')&cd
zfprintm_cd=: (lib,' bli_zfprintm ',ifw,' n & &c x x &j x x &c &c')&cd

NB. void bli_?printv
NB.      (
NB.        char*   s1,
NB.        dim_t   m,
NB.        ctype*  x, inc_t incx,
NB.        char*   format,
NB.        char*   s2
NB.      );

dprintv_cd=: (lib,' bli_dprintv ',ifw,' n &c x &d x &c &c')&cd
zprintv_cd=: (lib,' bli_zprintv ',ifw,' n &c x &j x &c &c')&cd

NB. void bli_?printm
NB.      (
NB.        char*   s1,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  a, inc_t rs_a, inc_t cs_a,
NB.        char*   format,
NB.        char*   s2
NB.      );

dprintm_cd=: (lib,' bli_dprintm ',ifw,' n &c x x &d x x &c &c')&cd
zprintm_cd=: (lib,' bli_zprintm ',ifw,' n &c x x &j x x &c &c')&cd

NB. void bli_?randv
NB.      (
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx
NB.      );

drandv_cd=: (lib,' bli_drandv ',ifw,' n x *d x')&cd
zrandv_cd=: (lib,' bli_zrandv ',ifw,' n x *j x')&cd

NB. void bli_?randm
NB.      (
NB.        doff_t  diagoffa,
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  a, inc_t rs_a, inc_t cs_a
NB.      );

drandm_cd=: (lib,' bli_drandm ',ifw,' n x x x x *d x x')&cd
zrandm_cd=: (lib,' bli_zrandm ',ifw,' n x x x x *j x x')&cd

NB. void bli_?sumsqv
NB.      (
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        rtype*  scale,
NB.        rtype*  sumsq
NB.      );

dsumsqv_cd=: (lib,' bli_dsumsqv ',ifw,' n x &d x *d *d')&cd
zsumsqv_cd=: (lib,' bli_zsumsqv ',ifw,' n x &j x *d *d')&cd

NB. void bli_?getsc
NB.      (
NB.        ctype*   chi,
NB.        double*  zeta_r,
NB.        double*  zeta_i
NB.      );

dgetsc_cd=: (lib,' bli_dgetsc ',ifw,' n &d *d *d')&cd
zgetsc_cd=: (lib,' bli_zgetsc ',ifw,' n &j *d *d')&cd

NB. err_t bli_?getijv
NB.      (
NB.        dim_t    i,
NB.        ctype*   x, inc_t incx,
NB.        double*  ar,
NB.        double*  ai
NB.      );

dgetijv_cd=: (lib,' bli_dgetijv ',ifw,' x x &d x *d *d')&cd
zgetijv_cd=: (lib,' bli_zgetijv ',ifw,' x x &j x *d *d')&cd

NB. err_t bli_?getijm
NB.      (
NB.        dim_t    i,
NB.        dim_t    j,
NB.        ctype*   b, inc_t rs_b, inc_t cs_b,
NB.        double*  ar,
NB.        double*  ai
NB.      );

dgetijm_cd=: (lib,' bli_dgetijm ',ifw,' x x x &d x x *d *d')&cd
zgetijm_cd=: (lib,' bli_zgetijm ',ifw,' x x x &j x x *d *d')&cd

NB. void bli_?setsc
NB.      (
NB.        double  zeta_r,
NB.        double  zeta_i,
NB.        ctype*  chi
NB.      );

dsetsc_cd=: (lib,' bli_dsetsc ',ifw,' n d d *d')&cd
zsetsc_cd=: (lib,' bli_zsetsc ',ifw,' n d d *j')&cd

NB. err_t bli_?setijv
NB.      (
NB.        double  ar,
NB.        double  ai,
NB.        dim_t   i,
NB.        ctype*  x, inc_t incx
NB.      );

dsetijv_cd=: (lib,' bli_dsetijv ',ifw,' x d d x *d x')&cd
zsetijv_cd=: (lib,' bli_zsetijv ',ifw,' x d d x *j x')&cd

NB. err_t bli_?setijm
NB.      (
NB.        double  ar,
NB.        double  ai,
NB.        dim_t   i,
NB.        dim_t   j,
NB.        ctype*  b, inc_t rs_b, inc_t cs_b
NB.      );

dsetijm_cd=: (lib,' bli_dsetijm ',ifw,' x d d x x *d x x')&cd
zsetijm_cd=: (lib,' bli_zsetijm ',ifw,' x d d x x *j x x')&cd

NB. void bli_?eqsc
NB.      (
NB.        conj_t  conjchi,
NB.        ctype*  chi,
NB.        ctype*  psi,
NB.        bool*   is_eq
NB.      );

deqsc_cd=: (lib,' bli_deqsc ',ifw,' n x &d &d *b')&cd
zeqsc_cd=: (lib,' bli_zeqsc ',ifw,' n x &j &j *b')&cd

NB. void bli_?eqv
NB.      (
NB.        conj_t  conjx,
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy,
NB.        bool*   is_eq
NB.      );

deqv_cd=: (lib,' bli_deqv ',ifw,' n x x &d x &d x *b')&cd
zeqv_cd=: (lib,' bli_zeqv ',ifw,' n x x &j x &j x *b')&cd

NB. void bli_?eqm
NB.      (
NB.        doff_t   diagoffa,
NB.        diag_t   diaga,
NB.        uplo_t   uploa,
NB.        trans_t  transa,
NB.        dim_t    m,
NB.        dim_t    n,
NB.        ctype*   a, inc_t rs_a, inc_t cs_a,
NB.        ctype*   b, inc_t rs_b, inc_t cs_b,
NB.        bool*    is_eq
NB.      );

deqm_cd=: (lib,' bli_deqm ',ifw,' n x x x x x x &d x x &d x x *b')&cd
zeqm_cd=: (lib,' bli_zeqm ',ifw,' n x x x x x x &j x x &j x x *b')&cd

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Query function reference

NB. char* bli_info_get_int_type_size_str( void );
info_get_int_type_size_str_cd=: (lib,' bli_info_get_int_type_size_str > ',ifw,' x')&cd

NB. ---------------------------------------------------------
NB. General library information

NB. char* bli_info_get_version_str( void );
info_get_version_str_cd=: (lib,' bli_info_get_version_str > ',ifw,' x')&cd

NB. ---------------------------------------------------------
NB. Specific configuration

NB. arch_t bli_arch_query_id( void );
arch_query_id_cd=: (lib,' bli_arch_query_id > ',ifw,' x'  )&cd

NB. char* bli_arch_string( arch_t id );
NB. application:
NB.      0 _1 2 memr@,~ arch_string_cd_mtbli_ arch_query_id_cd_mtbli_ ''
NB.   haswell
arch_string_cd=:   (lib,' bli_arch_string > '  ,ifw,' x x')&cd

NB. ---------------------------------------------------------
NB. General configuration

NB. gint_t bli_info_get_int_type_size( void );
info_get_int_type_size_cd=:          (lib,' bli_info_get_int_type_size > '         ,ifw,' x')&cd

NB. gint_t bli_info_get_num_fp_types( void );
info_get_num_fp_types_cd=:           (lib,' bli_info_get_num_fp_types > '          ,ifw,' x')&cd

NB. gint_t bli_info_get_max_type_size( void );
info_get_max_type_size_cd=:          (lib,' bli_info_get_max_type_size > '         ,ifw,' x')&cd

NB. gint_t bli_info_get_page_size( void );
info_get_page_size_cd=:              (lib,' bli_info_get_page_size > '             ,ifw,' x')&cd

NB. gint_t bli_info_get_simd_size( void );
info_get_simd_num_registers_cd=:     (lib,' bli_info_get_simd_num_registers > '    ,ifw,' x')&cd

NB. gint_t bli_info_get_simd_num_registers( void );
info_get_simd_size_cd=:              (lib,' bli_info_get_simd_size > '             ,ifw,' x')&cd

NB. gint_t bli_info_get_simd_align_size( void );
info_get_simd_align_size_cd=:        (lib,' bli_info_get_simd_align_size > '       ,ifw,' x')&cd

NB. gint_t bli_info_get_stack_buf_max_size( void );
info_get_stack_buf_max_size_cd=:     (lib,' bli_info_get_stack_buf_max_size > '    ,ifw,' x')&cd

NB. gint_t bli_info_get_stack_buf_align_size( void );
info_get_stack_buf_align_size_cd=:   (lib,' bli_info_get_stack_buf_align_size > '  ,ifw,' x')&cd

NB. gint_t bli_info_get_heap_addr_align_size( void );
info_get_heap_addr_align_size_cd=:   (lib,' bli_info_get_heap_addr_align_size > '  ,ifw,' x')&cd

NB. gint_t bli_info_get_heap_stride_align_size( void );
info_get_heap_stride_align_size_cd=: (lib,' bli_info_get_heap_stride_align_size > ',ifw,' x')&cd

NB. gint_t bli_info_get_pool_addr_align_size( void );
info_get_pool_addr_align_size_cd=:   (lib,' bli_info_get_pool_addr_align_size > '  ,ifw,' x')&cd

NB. gint_t bli_info_get_enable_stay_auto_init( void );
info_get_enable_stay_auto_init_cd=:  (lib,' bli_info_get_enable_stay_auto_init > ' ,ifw,' x')&cd

NB. gint_t bli_info_get_enable_blas( void );
info_get_enable_blas_cd=:            (lib,' bli_info_get_enable_blas > '           ,ifw,' x')&cd

NB. gint_t bli_info_get_blas_int_type_size( void );
info_get_blas_int_type_size_cd=:     (lib,' bli_info_get_blas_int_type_size > '    ,ifw,' x')&cd

NB. ---------------------------------------------------------
NB. Kernel information

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Micro-kernel implementation type query

NB. Possible implementation (ie: the ind_t method argument) types are:
BLIS_1M=:  0  NB. Implementation based on the 1m method. (This is the default induced method when real domain kernels are present but complex kernels are missing.)
BLIS_NAT=: 1  NB. Implementation based on "native" execution (ie: NOT an induced method).

NB. Possible microkernel types (ie: the return values for bli_info_get_*_ukr_impl_string()) are:
REFERENCE_UKERNEL=: 0  NB. ("refrnce"): This value is returned when the queried microkernel is provided by the reference implementation.
VIRTUAL_UKERNEL=:   1  NB. ("virtual"): This value is returned when the queried microkernel is driven by a the "virtual" microkernel provided by an induced method. This happens for any method value that is not BLIS_NAT (ie: native), but only applies to the complex domain.
OPTIMIZED_UKERNEL=: 2  NB. ("optimzd"): This value is returned when the queried microkernel is provided by an implementation that is neither reference nor virtual, and thus we assume the kernel author would deem it to be "optimized". Such a microkernel may not be optimal in the literal sense of the word, but nonetheless is intended to be optimized, at least relative to the reference microkernels.
NOTAPPLIC_UKERNEL=: 3  NB. ("notappl"): This value is returned usually when performing a gemmtrsm or trsm microkernel type query for any method value that is not BLIS_NAT (ie: native). That is, induced methods cannot be (purely) used on trsm-based microkernels because these microkernels perform more a triangular inversion, which is not matrix multiplication.

NB. char* bli_info_get_gemm_ukr_impl_string( ind_t method, num_t dt )
info_get_gemm_ukr_impl_string_cd=:       (lib,' bli_info_get_gemm_ukr_impl_string > '      ,ifw,' x x x')&cd

NB. char* bli_info_get_gemmtrsm_l_ukr_impl_string( ind_t method, num_t dt )
info_get_gemmtrsm_l_ukr_impl_string_cd=: (lib,' bli_info_get_gemmtrsm_l_ukr_impl_string > ',ifw,' x x x')&cd

NB. char* bli_info_get_gemmtrsm_u_ukr_impl_string( ind_t method, num_t dt )
info_get_gemmtrsm_u_ukr_impl_string_cd=: (lib,' bli_info_get_gemmtrsm_u_ukr_impl_string > ',ifw,' x x x')&cd

NB. char* bli_info_get_trsm_l_ukr_impl_string( ind_t method, num_t dt )
info_get_trsm_l_ukr_impl_string_cd=:     (lib,' bli_info_get_trsm_l_ukr_impl_string > '    ,ifw,' x x x')&cd

NB. char* bli_info_get_trsm_u_ukr_impl_string( ind_t method, num_t dt )
info_get_trsm_u_ukr_impl_string_cd=:     (lib,' bli_info_get_trsm_u_ukr_impl_string > '    ,ifw,' x x x')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Operation implementation type query

NB. char* bli_info_get_gemm_impl_string( num_t dt );
info_get_gemm_impl_string_cd=:  (lib,' bli_info_get_gemm_impl_string > ' ,ifw,' x x')&cd

NB. char* bli_info_get_hemm_impl_string( num_t dt );
info_get_hemm_impl_string_cd=:  (lib,' bli_info_get_hemm_impl_string > ' ,ifw,' x x')&cd

NB. char* bli_info_get_herk_impl_string( num_t dt );
info_get_herk_impl_string_cd=:  (lib,' bli_info_get_herk_impl_string > ' ,ifw,' x x')&cd

NB. char* bli_info_get_her2k_impl_string( num_t dt );
info_get_her2k_impl_string_cd=: (lib,' bli_info_get_her2k_impl_string > ',ifw,' x x')&cd

NB. char* bli_info_get_symm_impl_string( num_t dt );
info_get_symm_impl_string_cd=:  (lib,' bli_info_get_symm_impl_string > ' ,ifw,' x x')&cd

NB. char* bli_info_get_syrk_impl_string( num_t dt );
info_get_syrk_impl_string_cd=:  (lib,' bli_info_get_syrk_impl_string > ' ,ifw,' x x')&cd

NB. char* bli_info_get_syr2k_impl_string( num_t dt );
info_get_syr2k_impl_string_cd=: (lib,' bli_info_get_syr2k_impl_string > ',ifw,' x x')&cd

NB. char* bli_info_get_trmm_impl_string( num_t dt );
info_get_trmm_impl_string_cd=:  (lib,' bli_info_get_trmm_impl_string > ' ,ifw,' x x')&cd

NB. char* bli_info_get_trmm3_impl_string( num_t dt );
info_get_trmm3_impl_string_cd=: (lib,' bli_info_get_trmm3_impl_string > ',ifw,' x x')&cd

NB. char* bli_info_get_trsm_impl_string( num_t dt );
info_get_trsm_impl_string_cd=:  (lib,' bli_info_get_trsm_impl_string > ' ,ifw,' x x')&cd

NB. =========================================================
NB. Clean-up

erase 'lib ifw'
