NB. API definitions
NB.
NB. SZB       BLIS integer datatype size, in bytes
NB. BINT      J datatype ID for BLIS integer datatype
NB. BIC       15!:0 parameter code for BLIS integer datatype
NB. *_OFFS_*  Field offsets in BLIS structures
NB. *_SIZE    sizeof() BLIS structures
NB.
NB. xxxxxcd  Cover verbs to call BLIS subroutine or function
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
NB. Concepts
NB.
NB. Conventions:
NB. 1) LIB_mtbli_ global noun must exist
NB. 2) BLIS datatypes sizeof():
NB.      arch_t               4
NB.      atom_t               16
NB.      bli_pthread_mutex_t  SZI*4+8
NB.      conj_t               4
NB.      diag_t               4
NB.      dim_t                SZB
NB.      doff_t               SZB
NB.      dom_t                4
NB.      err_t                4
NB.      gint_t               SZB
NB.      inc_t                SZB
NB.      ind_t                4
NB.      invdiag_t            4
NB.      kerid_t              4
NB.      kimpl_t              4
NB.      mdim_t               4
NB.      num_t                4
NB.      objbits_t            4
NB.      pack_t               4
NB.      packbuf_t            4
NB.      packord_t            4
NB.      prec_t               4
NB.      side_t               4
NB.      siz_t                SZB
NB.      struc_t              4
NB.      timpl_t              4
NB.      trans_t              4
NB.      uplo_t               4

NB. =========================================================
NB. Configuration

coclass 'mtbli'

NB. =========================================================
NB. Local definitions

lib=. dquote LIB
ifw=. IFWIN # '+'

NB. bitwise operations with integers
NOT=: 2b011010 b. : [:    NB. Not
SFT=: [: : (2b100001 b.)  NB. Shift

NB. *********************************************************
NB. BLIS types

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. info bit field offsets

COMP_PREC_SHIFT=: 22
SCALAR_DT_SHIFT=: 24

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. info bit field masks

DATATYPE_BITS=:         JINT4 c. 2b00000000000000000000000000000111
  DOMAIN_BIT=:          JINT4 c. 2b00000000000000000000000000000001
  PRECISION_BIT=:       JINT4 c. 2b00000000000000000000000000000110
CONJTRANS_BITS=:        JINT4 c. 2b00000000000000000000000000011000
  TRANS_BIT=:           JINT4 c. 2b00000000000000000000000000001000
  CONJ_BIT=:            JINT4 c. 2b00000000000000000000000000010000
UPLO_BITS=:             JINT4 c. 2b00000000000000000000000011100000
  UPPER_BIT=:           JINT4 c. 2b00000000000000000000000000100000
  DIAG_BIT=:            JINT4 c. 2b00000000000000000000000001000000
  LOWER_BIT=:           JINT4 c. 2b00000000000000000000000010000000
UNIT_DIAG_BIT=:         JINT4 c. 2b00000000000000000000000100000000
INVERT_DIAG_BIT=:       JINT4 c. 2b00000000000000000000001000000000
PACK_SCHEMA_BITS=:      JINT4 c. 2b00000000000000001111110000000000
  PACK_PANEL_BIT=:      JINT4 c. 2b00000000000000000000010000000000
  PACK_FORMAT_BITS=:    JINT4 c. 2b00000000000000000111100000000000
  PACK_BIT=:            JINT4 c. 2b00000000000000001000000000000000
PACK_REV_IF_UPPER_BIT=: JINT4 c. 2b00000000000000010000000000000000
PACK_REV_IF_LOWER_BIT=: JINT4 c. 2b00000000000000100000000000000000
PACK_BUFFER_BITS=:      JINT4 c. 2b00000000000011000000000000000000
STRUC_BITS=:            JINT4 c. 2b00000000001100000000000000000000
COMP_PREC_BIT=:         JINT4 c. 2b00000000110000000000000000000000
SCALAR_DT_BITS=:        JINT4 c. 2b00000111000000000000000000000000
  SCALAR_DOMAIN_BIT=:   JINT4 c. 2b00000001000000000000000000000000
  SCALAR_PREC_BIT=:     JINT4 c. 2b00000110000000000000000000000000

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Enumerated type value definitions

BITVAL_REAL=:               JINT4 c. 2b00000000000000000000000000000000
BITVAL_COMPLEX=:            JINT4 c. 2b00000000000000000000000000000001
BITVAL_SINGLE_PREC=:        JINT4 c. 2b00000000000000000000000000000000
BITVAL_DOUBLE_PREC=:        JINT4 c. 2b00000000000000000000000000000010
  BITVAL_FLOAT_TYPE=:       JINT4 c. 2b00000000000000000000000000000000
  BITVAL_SCOMPLEX_TYPE=:    JINT4 c. 2b00000000000000000000000000000001
  BITVAL_DOUBLE_TYPE=:      JINT4 c. 2b00000000000000000000000000000010
  BITVAL_DCOMPLEX_TYPE=:    JINT4 c. 2b00000000000000000000000000000011
  BITVAL_INT_TYPE=:         JINT4 c. 2b00000000000000000000000000000100
  BITVAL_CONST_TYPE=:       JINT4 c. 2b00000000000000000000000000000101
BITVAL_NO_TRANS=:           JINT4 c. 2b00000000000000000000000000000000
BITVAL_TRANS=:              JINT4 c. 2b00000000000000000000000000001000
BITVAL_NO_CONJ=:            JINT4 c. 2b00000000000000000000000000000000
BITVAL_CONJ=:               JINT4 c. 2b00000000000000000000000000010000
BITVAL_CONJ_TRANS=:         JINT4 c. 2b00000000000000000000000000011000
BITVAL_ZEROS=:              JINT4 c. 2b00000000000000000000000000000000
BITVAL_UPPER=:              JINT4 c. 2b00000000000000000000000001100000
BITVAL_LOWER=:              JINT4 c. 2b00000000000000000000000011000000
BITVAL_DENSE=:              JINT4 c. 2b00000000000000000000000011100000
BITVAL_NONUNIT_DIAG=:       JINT4 c. 2b00000000000000000000000000000000
BITVAL_UNIT_DIAG=:          JINT4 c. 2b00000000000000000000000100000000
BITVAL_INVERT_DIAG=:        JINT4 c. 2b00000000000000000000001000000000
BITVAL_NOT_PACKED=:         JINT4 c. 2b00000000000000000000000000000000
  BITVAL_1E=:               JINT4 c. 2b00000000000000000000100000000000
  BITVAL_1R=:               JINT4 c. 2b00000000000000000001000000000000
  BITVAL_RO=:               JINT4 c. 2b00000000000000000001100000000000
  BITVAL_PACKED_UNSPEC=:    JINT4 c. 2b00000000000000001000000000000000
  BITVAL_PACKED_PANELS=:    JINT4 c. 2b00000000000000001000010000000000
  BITVAL_PACKED_PANELS_1E=: JINT4 c. 2b00000000000000001000110000000000
  BITVAL_PACKED_PANELS_1R=: JINT4 c. 2b00000000000000001001010000000000
  BITVAL_PACKED_PANELS_RO=: JINT4 c. 2b00000000000000001001110000000000
BITVAL_PACK_FWD_IF_UPPER=:  JINT4 c. 2b00000000000000000000000000000000
BITVAL_PACK_REV_IF_UPPER=:  JINT4 c. 2b00000000000000010000000000000000
BITVAL_PACK_FWD_IF_LOWER=:  JINT4 c. 2b00000000000000000000000000000000
BITVAL_PACK_REV_IF_LOWER=:  JINT4 c. 2b00000000000000100000000000000000
BITVAL_BUFFER_FOR_A_BLOCK=: JINT4 c. 2b00000000000000000000000000000000
BITVAL_BUFFER_FOR_B_PANEL=: JINT4 c. 2b00000000000001000000000000000000
BITVAL_BUFFER_FOR_C_PANEL=: JINT4 c. 2b00000000000010000000000000000000
BITVAL_BUFFER_FOR_GEN_USE=: JINT4 c. 2b00000000000011000000000000000000
BITVAL_GENERAL=:            JINT4 c. 2b00000000000000000000000000000000
BITVAL_HERMITIAN=:          JINT4 c. 2b00000000000100000000000000000000
BITVAL_SYMMETRIC=:          JINT4 c. 2b00000000001000000000000000000000
BITVAL_TRIANGULAR=:         JINT4 c. 2b00000000001100000000000000000000

NB. =========================================================
NB. Interface

NB. *********************************************************
NB. BLIS types

NB. BLIS integer datatype
NB.
NB.   BLIS_INT_TYPE_SIZE  SZB  BINT   BIC
NB.   ------------------  ---  -----  ----
NB.   32                  4    JINT4  'i'
NB.   64                  8    JINT   'l'

NB. const char* bli_info_get_int_type_size_str( void );
NB. note: see info_get_int_type_size_str in external/blis/util.ijs
info_get_int_type_size_str_cd=: (lib,' bli_info_get_int_type_size_str > ',ifw,' x')&cd

SZB=: _3 SFT ". info_get_int_type_size_str ''  NB. sizeof()
('external/blis/api: BLIS_INT_TYPE_SIZE = ' , (": 3 SFT SZB) , ' isn''t supported') assert SZB e. 4 8

BINT=: JINT4 [^:(4 = SZB) JINT  NB. J datatype ID

BIC=: 'i'    [^:(4 = SZB) 'l'   NB. 15!:0 parameter code

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Enumerated type definitions

NB. ---------------------------------------------------------
NB. Operational parameter types

NB. trans_t
NB. Semantic meaning: Corresponding matrix operand...

NO_TRANSPOSE=:       BITVAL_NO_TRANS       NB. ...will be used as given.
TRANSPOSE=:          BITVAL_TRANS          NB. ...will be implicitly transposed.
CONJ_NO_TRANSPOSE=:  BITVAL_CONJ           NB. ...will be implicitly conjugated.
CONJ_TRANSPOSE=:     BITVAL_CONJ_TRANS     NB. ...will be implicitly transposed and conjugated.

NB. conj_t
NB. Semantic meaning: Corresponding matrix/vector operand...

NO_CONJUGATE=:       BITVAL_NO_CONJ        NB. ...will be used as given.
CONJUGATE=:          BITVAL_CONJ           NB. ...will be implicitly conjugated.

NB. uplo_t
NB. Semantic meaning: Corresponding matrix operand...

ZEROS=:              BITVAL_ZEROS          NB. ...is filled by zeros.
LOWER=:              BITVAL_LOWER          NB. ...is stored in (and will be accessed only from) the lower triangle.
UPPER=:              BITVAL_UPPER          NB. ...is stored in (and will be accessed only from) the upper triangle.
DENSE=:              BITVAL_DENSE          NB. ...is stored as a full matrix (ie: in both triangles).

NB. side_t
NB. Semantic meaning: Corresponding matrix operand...

LEFT=:               JINT4 c. 0            NB. ...appears on the left.
RIGHT=:              JINT4 c. 1            NB. ...appears on the right.

NB. diag_t
NB. Semantic meaning: Corresponding matrix operand...

NONUNIT_DIAG=:       BITVAL_NONUNIT_DIAG   NB. ...has a non-unit diagonal that should be explicitly read from.
UNIT_DIAG=:          BITVAL_UNIT_DIAG      NB. ...has a unit diagonal that should be implicitly assumed (and not read from).

NB. invdiag_t
NB. Semantic meaning: Corresponding matrix operand...

NO_INVERT_DIAG=:     JINT4 c. 0            NB. ...has a diagonal that should be explicitly read from.
INVERT_DIAG=:        BITVAL_INVERT_DIAG    NB. ...has a diagonal that should be inverted.

NB. struc_t
NB. Semantic meaning: Matrix operand...

GENERAL=:            BITVAL_GENERAL        NB. ...has no structure.
HERMITIAN=:          BITVAL_HERMITIAN      NB. ...has Hermitian structure.
SYMMETRIC=:          BITVAL_SYMMETRIC      NB. ...has symmetric structure.
TRIANGULAR=:         BITVAL_TRIANGULAR     NB. ...has triangular structure.

NB. ---------------------------------------------------------
NB. Data types

NB. num_t
NB. Semantic meaning: Matrix/vector operand...

FLOAT=:              BITVAL_FLOAT_TYPE     NB. ...contains single-precision real elements.
DOUBLE=:             BITVAL_DOUBLE_TYPE    NB. ...contains double-precision real elements.
SCOMPLEX=:           BITVAL_SCOMPLEX_TYPE  NB. ...contains single-precision complex elements.
DCOMPLEX=:           BITVAL_DCOMPLEX_TYPE  NB. ...contains double-precision complex elements.
INT=:                BITVAL_INT_TYPE       NB. ...contains integer elements of type gint_t.
CONSTANT=:           BITVAL_CONST_TYPE     NB. ...contains polymorphic representation of a constant value.

NB. dom_t
NB. Semantic meaning: Matrix/vector operand...

REAL=:               BITVAL_REAL           NB. ...contains real domain elements.
COMPLEX=:            BITVAL_COMPLEX        NB. ...contains complex domain elements.

NB. prec_t
NB. Semantic meaning: Matrix/vector operand...

SINGLE_PREC=:        BITVAL_SINGLE_PREC    NB. ...contains single-precision elements.
DOUBLE_PREC=:        BITVAL_DOUBLE_PREC    NB. ...contains double-precision elements.

NB. pack_t
NB. Semantic meaning: Pack schema type

NOT_PACKED=:         BITVAL_NOT_PACKED
PACKED_UNSPEC=:      BITVAL_PACKED_UNSPEC
PACKED_VECTOR=:      BITVAL_PACKED_UNSPEC
PACKED_MATRIX=:      BITVAL_PACKED_UNSPEC
PACKED_PANELS=:      BITVAL_PACKED_PANELS
PACKED_PANELS_1E=:   BITVAL_PACKED_PANELS_1E
PACKED_PANELS_1R=:   BITVAL_PACKED_PANELS_1R
PACKED_PANELS_RO=:   BITVAL_PACKED_PANELS_RO

NB. packord_t
NB. Semantic meaning: Pack order type

PACK_FWD_IF_UPPER=:  BITVAL_PACK_FWD_IF_UPPER
PACK_REV_IF_UPPER=:  BITVAL_PACK_REV_IF_UPPER
PACK_FWD_IF_LOWER=:  BITVAL_PACK_FWD_IF_LOWER
PACK_REV_IF_LOWER=:  BITVAL_PACK_REV_IF_LOWER

NB. packbuf_t
NB. Semantic meaning: Pack buffer type

BUFFER_FOR_A_BLOCK=: BITVAL_BUFFER_FOR_A_BLOCK
BUFFER_FOR_B_PANEL=: BITVAL_BUFFER_FOR_B_PANEL
BUFFER_FOR_C_PANEL=: BITVAL_BUFFER_FOR_C_PANEL
BUFFER_FOR_GEN_USE=: BITVAL_BUFFER_FOR_GEN_USE

NB. mdim_t
NB. Semantic meaning: Matrix dimension type

M=:                  JINT4 c. 0
N=:                  JINT4 c. 1

NB. ---------------------------------------------------------
NB. Implementation types

NB. ind_t
NB. Semantic meaning: Induced method type

M1=:                 JINT4 c. 0            NB. Implementation based on the 1m method. (This is the default induced method when real domain kernels are present but complex kernels are missing.)
NAT=:                JINT4 c. 1            NB. Implementation based on "native" execution (ie: NOT an induced method).

NB. timpl_t
NB. Semantic meaning: Threading implementation type

SINGLE=:             JINT4 c. 0            NB. disable multithreading
OPENMP=:             JINT4 c. 1            NB. enable threading via OpenMP
POSIX=:              JINT4 c. 2            NB. enable threading via pthreads
HPX=:                JINT4 c. 3            NB. enable threading via HPX

NB. kimpl_t
NB. Semantic meaning: Possible microkernel types (ie: the
NB. return values for bli_info_get_*_ukr_impl_string())

REFERENCE_UKERNEL=:  JINT4 c. 0            NB. ("refrnce"): This value is returned when the queried microkernel is provided by the reference implementation.
VIRTUAL_UKERNEL=:    JINT4 c. 1            NB. ("virtual"): This value is returned when the queried microkernel is driven by a the "virtual" microkernel provided by an induced method. This happens for any method value that is not BLIS_NAT (ie: native), but only applies to the complex domain.
OPTIMIZED_UKERNEL=:  JINT4 c. 2            NB. ("optimzd"): This value is returned when the queried microkernel is provided by an implementation that is neither reference nor virtual, and thus we assume the kernel author would deem it to be "optimized". Such a microkernel may not be optimal in the literal sense of the word, but nonetheless is intended to be optimized, at least relative to the reference microkernels.
NOTAPPLIC_UKERNEL=:  JINT4 c. 3            NB. ("notappl"): This value is returned usually when performing a gemmtrsm or trsm microkernel type query for any method value that is not BLIS_NAT (ie: native). That is, induced methods cannot be (purely) used on trsm-based microkernels because these microkernels perform more a triangular inversion, which is not matrix multiplication.

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Structures

NB. ---------------------------------------------------------
NB. Stack type
NB.
NB. typedef struct
NB. {                              offset                        length              notes
NB.   siz_t elem_size;             0                             SZB
NB.   siz_t block_len;             SZB                           SZB
NB.   siz_t max_blocks;            SZB*2                         SZB
NB.   siz_t size;                  SZB*3                         SZB
NB.   siz_t capacity;              SZB*4                         SZB
NB.                                SZB*5                         SZI-min(SZI,SZB)    align
NB.   void** blocks;               SZI+SZB*5-min(SZI,SZB)        SZI
NB.
NB.   bli_pthread_mutex_t lock;    SZI*2+SZB*5-min(SZI,SZB)      SZI*4+8
NB.                                SZI*6+SZB*5-min(SZI,SZB)+8                        sizeof(stck_t)
NB. } stck_t;

'STCK_T_OFFS_BL STCK_T_OFFS_MB STCK_T_OFFS_SZ STCK_T_OFFS_CP trash STCK_T_OFFS_BS STCK_T_OFFS_LK STCK_T_SIZE'=: +/\ SZI ((5 # ]) , ([ - <.) , [ , 8 + 4 * [) SZB

NB. ---------------------------------------------------------
NB. Context type
NB.
NB. typedef struct cntx_s
NB. {                             offset                             sizeof()                      notes
NB.   stck_t  blkszs;             0                                  SZI*6+SZB*5-min(SZI,SZB)+8
NB.   stck_t  bmults;             SZI* 6+SZB* 5-min(SZI,SZB)  + 8    SZI*6+SZB*5-min(SZI,SZB)+8
NB.
NB.   stck_t  ukrs;               SZI*12+SZB*10-min(SZI,SZB)*2+16    SZI*6+SZB*5-min(SZI,SZB)+8
NB.   stck_t  ukr2s;              SZI*18+SZB*15-min(SZI,SZB)*3+24    SZI*6+SZB*5-min(SZI,SZB)+8
NB.   stck_t  ukr_prefs;          SZI*24+SZB*20-min(SZI,SZB)*4+32    SZI*6+SZB*5-min(SZI,SZB)+8
NB.
NB.   stck_t  l3_sup_handlers;    SZI*30+SZB*25-min(SZI,SZB)*5+40    SZI*6+SZB*5-min(SZI,SZB)+8
NB.                               SZI*36+SZB*30-min(SZI,SZB)*6+48                                  sizeof(cntx_t)
NB. } cntx_t;

'CNTX_T_OFFS_BM CNTX_T_OFFS_US CNTX_T_OFFS_U2 CNTX_T_OFFS_UP CNTX_T_OFFS_L3 CNTX_T_SIZE'=: +/\ 6 # SZI (8 + ((+ 6&*)~ 5&*) - <.) SZB

NB. ---------------------------------------------------------
NB. Runtime type
NB.
NB. typedef struct rntm_s
NB. {
NB.   // "External" fields: these may be queried by the end-user.                  offset                   sizeof()          notes
NB.   timpl_t   thread_impl;                                                       0                        4
NB.
NB.   bool      auto_factor;                                                       4                        1
NB.                                                                                5                        3                 align
NB.   dim_t     num_threads;                                                       8                        SZB
NB.   dim_t     thrloop[ BLIS_NUM_LOOPS ];                                         SZB               + 8    SZB*6
NB.   bool      pack_a; // enable/disable packing of left-hand matrix A.           SZB*7             + 8    1
NB.   bool      pack_b; // enable/disable packing of right-hand matrix B.          SZB*7             + 9    1
NB.   bool      l3_sup; // enable/disable small matrix handling in level-3 ops.    SZB*7             +10    1
NB.                                                                                SZB*7             +11    min(SZI,SZB)-3    align
NB. } rntm_t;                                                                      SZB*7+min(SZI,SZB)+ 8                      sizeof(rntm_t)
NB.
NB. #define BLIS_NUM_LOOPS 6

'RNTM_T_OFFS_AF trash RNTM_T_OFFS_NT RNTM_T_OFFS_TL RNTM_T_OFFS_PA RNTM_T_OFFS_PB RNTM_T_OFFS_PC trash RNTM_T_SIZE'=: +/\ SZI (4 1 3 , (1 6 * ]) , 1 1 1 , _3 + <.) SZB

NB. ---------------------------------------------------------
NB. Object type
NB.
NB. typedef struct obj_s
NB. {
NB.   // Basic fields                 offset                          sizeof()            notes
NB.   struct obj_s* root;             0                               SZI
NB.
NB.   dim_t         off[2];           SZI                             SZB*2
NB.   dim_t         dim[2];           SZI  +SZB* 2                    SZB*2
NB.   doff_t        diag_off;         SZI  +SZB* 4                    SZB
NB.
NB.   objbits_t     info;             SZI  +SZB* 5                    4
NB.   objbits_t     info2;            SZI  +SZB* 5+ 4                 4
NB.   siz_t         elem_size;        SZI  +SZB* 5+ 8                 SZB
NB.
NB.   void*         buffer;           SZI  +SZB* 6+ 8                 SZI
NB.   inc_t         rs;               SZI*2+SZB* 6+ 8                 SZB
NB.   inc_t         cs;               SZI*2+SZB* 7+ 8                 SZB
NB.   inc_t         is;               SZI*2+SZB* 8+ 8                 SZB
NB.                                   SZI*2+SZB* 9+ 8                 SZI-min(SZI,SZB)    align
NB.   // Bufferless scalar storage
NB.   atom_t        scalar;           SZI*3+SZB* 9+ 8-min(SZI,SZB)    16
NB.
NB.   // Pack-related fields
NB.   dim_t         m_padded;         SZI*3+SZB* 9+24-min(SZI,SZB)    SZB
NB.   dim_t         n_padded;         SZI*3+SZB*10+24-min(SZI,SZB)    SZB
NB.   inc_t         ps;               SZI*3+SZB*11+24-min(SZI,SZB)    SZB
NB.   inc_t         pd;               SZI*3+SZB*12+24-min(SZI,SZB)    SZB
NB.   dim_t         m_panel;          SZI*3+SZB*13+24-min(SZI,SZB)    SZB
NB.   dim_t         n_panel;          SZI*3+SZB*14+24-min(SZI,SZB)    SZB
NB.                                   SZI*3+SZB*15+24-min(SZI,SZB)                        sizeof(obj_t)
NB. } obj_t;

'OBJ_T_OFFS_OF OBJ_T_OFFS_DM OBJ_T_OFFS_DO OBJ_T_OFFS_IN OBJ_T_OFFS_I2 OBJ_T_OFFS_ES OBJ_T_OFFS_BF OBJ_T_OFFS_RS OBJ_T_OFFS_CS OBJ_T_OFFS_IS trash OBJ_T_OFFS_SC OBJ_T_OFFS_MD OBJ_T_OFFS_ND OBJ_T_OFFS_PS OBJ_T_OFFS_PD OBJ_T_OFFS_ML OBJ_T_OFFS_NL OBJ_T_SIZE'=: +/\ SZI ([ , (2 2 1 * ]) , (4 4 , ,~) , (16 ,~ ([ - <.) ,~ 3 # ]) , 6 # ]) SZB

NB. *********************************************************
NB. Global scalar constants

TWO=:         0 {:: dlsym LIB ; 'BLIS_TWO'
ONE=:         0 {:: dlsym LIB ; 'BLIS_ONE'
ZERO=:        0 {:: dlsym LIB ; 'BLIS_ZERO'
MINUS_ONE=:   0 {:: dlsym LIB ; 'BLIS_MINUS_ONE'
MINUS_TWO=:   0 {:: dlsym LIB ; 'BLIS_MINUS_TWO'
ONE_I=:       0 {:: dlsym LIB ; 'BLIS_ONE_I'
MINUS_ONE_I=: 0 {:: dlsym LIB ; 'BLIS_MINUS_ONE_I'
NAN=:         0 {:: dlsym LIB ; 'BLIS_NAN'

NB. *********************************************************
NB. Functions

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Context

NB. C: cntx_t* bli_gks_query_cntx( void );
NB. J: cntx_addr=. gks_query_cntx_cd_mtbli_ ''
gks_query_cntx_cd=: (lib,' bli_gks_query_cntx > ',ifw,' ',BIC)&cd

NB. C: void bli_cntx_print( const cntx_t* cntx );
NB. J: trash=. cntx_print_cd < < cntx_addr
NB.
NB. Application:
NB.      cntx_addr=. gks_query_cntx_cd_mtbli_ ''
NB.      EMPTY [ cntx_print_cd < < cntx_addr  NB. outputs to stdout

cntx_print_cd=: (lib,' bli_cntx_print ',ifw,' n &')&cd

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Runtime

NB. C: void bli_rntm_init_from_global( rntm_t* rntm );
NB. J: trash=. rntm_init_from_global_cd < < rntm_addr
NB.
NB. Application:
NB.      rntm_addr=. mema SIZEOF_RNTM_T_mtbli_
NB.      EMPTY [ rntm_init_from_global_cd_mtbli_ < < rntm_addr
NB.      a. i. memr rntm_addr , 0 , SIZEOF_RNTM_T_mtbli_ , JCHAR
NB.   1 0 0 0 1 0 0 0 7 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0

rntm_init_from_global_cd=: (lib,' bli_rntm_init_from_global ',ifw,' n *')&cd

NB. C: void bli_rntm_print( const rntm_t* rntm );
NB. J: trash=. rntm_print_cd < < rntm_addr
NB.
NB. Application:
NB.      rntm_addr=. mema SIZEOF_RNTM_T_mtbli_
NB.      EMPTY [ rntm_init_from_global_cd_mtbli_ < < rntm_addr
NB.      rntm_print_cd_mtbli_ < < rntm_addr  NB. outputs to stdout
NB.   thread impl: 1
NB.   rntm contents    nt  jc  pc  ic  jr  ir
NB.   autofac? 1 |    7   1   1   1   1   1
NB.
NB. Notes:
NB. - BLIS library doesn't export bli_rntm_print name:
NB.        dlsym_mttst_ LIB_mtbli_ ; 'bli_rntm_print'
NB.     +-+-----------------------------------------------------------+
NB.     |0|/home/user/lib/libblis.so: undefined symbol: bli_rntm_print|
NB.     +-+-----------------------------------------------------------+
NB.   so rntm_print_cd implements it

NB. rntm_print_cd=: (lib,' bli_rntm_print ',ifw,' n &')&cd

rntm_print_cd=: 3 : 0
  assert ((2 = L.) , 1 = #) y
  y=. > y
  assert ((1 = L.) , 1 = #) y
  y=. > y
  assert ((-: <.) , 1 = #) y

  ti=.                     memr y , 0                1 , JINT4
  af=.                     memr y , RNTM_T_OFFS_AF , 1 , JB01
  nt=.                     memr y , RNTM_T_OFFS_NT , 1 , BINT
  'pr ir jr ic pc jc'=.    memr y , RNTM_T_OFFS_TL , 6 , BINT
  'pack_a pack_b l3_sup'=. memr y , RNTM_T_OFFS_PA , 3 , JB01
  echo 'thread impl: ' , (": ti) , LF , 'rntm contents    nt  jc  pc  ic  jr  ir' , LF , 'autofac? ' , (1 ": af) , ' | ' ,(4 ": nt , jc , pc , ic , jr , ir)
)

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Multi-threading

NB. ---------------------------------------------------------
NB. Globally at runtime

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Globally at runtime: the automatic way

NB. dim_t bli_thread_get_num_threads( void );
thread_get_num_threads_cd=: (lib,' bli_thread_get_num_threads > ',ifw,' ',BIC)&cd

NB. void bli_thread_set_num_threads( dim_t n_threads );
thread_set_num_threads_cd=: (lib,' bli_thread_set_num_threads ',ifw,' n ',BIC)&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Globally at runtime: the manual way

NB. void bli_thread_set_ways( dim_t jc, dim_t pc, dim_t ic, dim_t jr, dim_t ir );
thread_set_ways_cd=:  (lib,' bli_thread_set_ways ',ifw,' n ',BIC,' ',BIC,' ',BIC,' ',BIC,' ',BIC  )&cd

NB. dim_t bli_thread_get_jc_nt( void );
thread_get_jc_nt_cd=: (lib,' bli_thread_get_jc_nt > ',ifw,' ',BIC)&cd

NB. dim_t bli_thread_get_pc_nt( void );
thread_get_pc_nt_cd=: (lib,' bli_thread_get_pc_nt > ',ifw,' ',BIC)&cd

NB. dim_t bli_thread_get_ic_nt( void );
thread_get_ic_nt_cd=: (lib,' bli_thread_get_ic_nt > ',ifw,' ',BIC)&cd

NB. dim_t bli_thread_get_jr_nt( void );
thread_get_jr_nt_cd=: (lib,' bli_thread_get_jr_nt > ',ifw,' ',BIC)&cd

NB. dim_t bli_thread_get_ir_nt( void );
thread_get_ir_nt_cd=: (lib,' bli_thread_get_ir_nt > ',ifw,' ',BIC)&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Globally at runtime: overriding the default threading implementation

NB. void bli_thread_set_thread_impl( timpl_t ti );
thread_set_thread_impl_cd=: (lib,' bli_thread_set_thread_impl ',ifw,' n i')&cd

NB. timpl_t bli_thread_get_thread_impl( void );
thread_get_thread_impl_cd=: (lib,' bli_thread_get_thread_impl > ',ifw,' i')&cd

NB. const char* bli_thread_get_thread_impl_str( timpl_t ti );
NB. note: see thread_get_thread_impl_str in external/blis/util.ijs
thread_get_thread_impl_str_cd=: (lib,' bli_thread_get_thread_impl_str > ',ifw,' x i')&cd

NB. ---------------------------------------------------------
NB. Locally at runtime

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Initializing a rntm_t

NB. 1) via BLIS_RNTM_INITIALIZER macro

NB. rntm=. BLIS_RNTM_INITIALIZER ''
NB. Application:
NB.      rntm=. BLIS_RNTM_INITIALIZER_mtbli_ ''
NB.      NB. set threads number
NB.      rntm_set_num_threads_cd_mtbli_ 3 ; < < rntm
NB.      NB. destroy rntm
NB.      memf rntm

BLIS_RNTM_INITIALIZER=: 3 : 0
  rntm=. mema RNTM_T_SIZE
  SINGLE  memw rntm , 0                1 , JINT4
  0       memw rntm , RNTM_T_OFFS_AF , 1 , JB01
  1       memw rntm , RNTM_T_OFFS_NT , 1 , BINT
  (6 # 1) memw rntm , RNTM_T_OFFS_TL , 6 , BINT
  0 0 1   memw rntm , RNTM_T_OFFS_PA , 3 , JB01
  rntm
)

NB. 2) via bli_rntm_init_from_global() function: see above

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Locally at runtime: the automatic way

NB. void bli_rntm_set_num_threads( dim_t n_threads, rntm_t* rntm );
rntm_set_num_threads_cd=: (lib,' bli_rntm_set_num_threads ',ifw,' n ',BIC,' *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Locally at runtime: the manual way

NB. void bli_rntm_set_ways( dim_t jc, dim_t pc, dim_t ic, dim_t jr, dim_t ir, rntm_t* rntm );
rntm_set_ways_cd=: (lib,' bli_rntm_set_ways ',ifw,' n ',BIC,' ',BIC,' ',BIC,' ',BIC,' ',BIC,' *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Locally at runtime: overriding the default threading implementation

NB. void bli_rntm_set_thread_impl( timpl_t ti, rntm_t* rntm );
rntm_set_thread_impl_cd=: (lib,' bli_rntm_set_thread_impl ',ifw,' n i *')&cd

NB. ---------------------------------------------------------
NB. Other thread function reference

NB. void bli_thread_reset( void );
bli_thread_reset_cd=: (lib,' bli_thread_reset ',ifw,' n')&cd

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

obj_create_cd=: (lib,' bli_obj_create ',ifw,' n i ',BIC,' ',BIC,' ',BIC,' ',BIC,' *')&cd

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

obj_create_without_buffer_cd=: (lib,'bli_obj_create_without_buffer ',ifw,' n i ',BIC,' ',BIC,' *')&cd

NB. void bli_obj_attach_buffer
NB.      (
NB.        const void*   p,
NB.              inc_t   rs,
NB.              inc_t   cs,
NB.              inc_t   is,
NB.        const obj_t*  obj
NB.      );

obj_attach_buffer_cd=: (lib,'bli_obj_attach_buffer ',ifw,' n & ',BIC,' ',BIC,' ',BIC,' *')&cd

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

obj_create_with_attached_buffer_cd=: (lib,' bli_obj_create_with_attached_buffer ',ifw,' n i ',BIC,' ',BIC,' & ',BIC,' ',BIC,' *')&cd

NB. void bli_obj_alloc_buffer
NB.      (
NB.        inc_t   rs,
NB.        inc_t   cs,
NB.        inc_t   is,
NB.        obj_t*  obj
NB.      );

obj_alloc_buffer_cd=: (lib,'bli_obj_alloc_buffer ',ifw,' n ',BIC,' ',BIC,' ',BIC,' *')&cd

NB. void bli_obj_create_1x1
NB.      (
NB.        num_t   dt,
NB.        obj_t*  obj
NB.      );

obj_create_1x1_cd=: (lib,' bli_obj_create_1x1 ',ifw,' n i *')&cd

NB. void bli_obj_create_1x1_with_attached_buffer
NB.      (
NB.        num_t   dt,
NB.        void*   p,
NB.        obj_t*  obj
NB.      );

obj_create_1x1_with_attached_buffer_cd=: (lib,' bli_obj_create_1x1_with_attached_buffer ',ifw,' n i & *')&cd

NB. void bli_obj_create_conf_to
NB.      (
NB.        const obj_t*  s,
NB.              obj_t*  d
NB.      );

obj_create_conf_to_cd=: (lib,' bli_obj_create_conf_to ',ifw,' n & *')&cd

NB. void bli_obj_scalar_init_detached
NB.      (
NB.        num_t   dt,
NB.        obj_t*  obj
NB.      );

obj_scalar_init_detached_cd=: (lib,' bli_obj_scalar_init_detached ',ifw,' n i *')&cd

NB. ---------------------------------------------------------
NB. Object initialization function reference

NB. C: void bli_obj_init_full_shallow_copy_of( const obj_t* a, obj_t* b );
NB. J: trash=. a obj_init_full_shallow_copy_of b
NB. note: copy from a to b
obj_init_full_shallow_copy_of=: (memw~ memr)~&(,&(0 , OBJ_T_SIZE , JCHAR))

NB. Initialize object with default properties (info field)
NB. C: void bli_obj_set_defaults( obj_t* obj );
NB. J: trash=. obj_set_defaults obj
obj_set_defaults=: (BITVAL_DENSE OR BITVAL_GENERAL) memw ,&( OBJ_T_OFFS_IN            , 1 , JINT4)

NB. ---------------------------------------------------------
NB. Object accessor function reference

NB. info field

NB. num_t bli_obj_dt( const obj_t* obj );
obj_dt=:                                               DATATYPE_BITS         AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. dom_t bli_obj_domain( const obj_t* obj );
obj_domain=:                                           DOMAIN_BIT            AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. prec_t bli_obj_prec( const obj_t* obj );
obj_prec=:                                             PRECISION_BIT         AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. trans_t bli_obj_conjtrans_status( const obj_t* obj );
obj_conjtrans_status=:                                 CONJTRANS_BITS        AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. trans_t bli_obj_onlytrans_status( const obj_t* obj );
obj_onlytrans_status=:                                 TRANS_BIT             AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. conj_t bli_obj_conj_status( const obj_t* obj );
obj_conj_status=:                                      CONJ_BIT              AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. struc_t bli_obj_struc( const obj_t* obj );
obj_struc=:                                            STRUC_BITS            AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. uplo_t bli_obj_uplo( const obj_t* obj );
obj_uplo=:                                             UPLO_BITS             AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. diag_t bli_obj_diag( const obj_t* obj );
obj_diag=:                                             UNIT_DIAG_BIT         AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. bool bli_obj_has_inverted_diag( const obj_t* obj );
obj_has_inverted_diag=:    BITVAL_INVERT_DIAG       -: INVERT_DIAG_BIT       AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. bool bli_obj_is_pack_rev_if_upper( const obj_t* obj );
obj_is_pack_rev_if_upper=: BITVAL_PACK_REV_IF_UPPER -: PACK_REV_IF_UPPER_BIT AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. bool bli_obj_is_pack_rev_if_lower( const obj_t* obj );
obj_is_pack_rev_if_lower=: BITVAL_PACK_REV_IF_LOWER -: PACK_REV_IF_LOWER_BIT AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. pack_t bli_obj_pack_schema( const obj_t* obj );
obj_pack_schema=:                                      PACK_SCHEMA_BITS      AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. bool bli_obj_is_packed( const obj_t* obj );
obj_is_packed=:            (JINT4 c. 0)             ~: PACK_BIT              AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. bool bli_obj_is_panel_packed( const obj_t* obj );
obj_is_panel_packed=:      (JINT4 c. 0)             ~: PACK_PANEL_BIT        AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. packbuf_t bli_obj_pack_buffer_type( const obj_t* obj );
obj_pack_buffer_type=:                                 PACK_BUFFER_BITS      AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. prec_t bli_obj_comp_prec( const obj_t* obj );
obj_comp_prec=:            (- COMP_PREC_SHIFT)     SFT COMP_PREC_BIT         AND memr@(,&( OBJ_T_OFFS_IN          , 1 , JINT4))

NB. info2 field

NB. num_t bli_obj_scalar_dt( const obj_t* obj );
obj_scalar_dt=:            (- SCALAR_DT_SHIFT)     SFT SCALAR_DT_BITS        AND memr@(,&( OBJ_T_OFFS_I2          , 1 , JINT4))

NB. dom_t bli_obj_scalar_domain( const obj_t* obj );
obj_scalar_domain=:        (- SCALAR_DT_SHIFT)     SFT SCALAR_DOMAIN_BIT     AND memr@(,&( OBJ_T_OFFS_I2          , 1 , JINT4))

NB. prec_t bli_obj_scalar_prec( const obj_t* obj );
obj_scalar_prec=:          (- SCALAR_DT_SHIFT)     SFT SCALAR_PREC_BIT       AND memr@(,&( OBJ_T_OFFS_I2          , 1 , JINT4))

NB. other fields

NB. obj_t* bli_obj_root( const obj_t* obj );
obj_root=:                                                                       memr@(,&( 0                      , 1 , JPTR ))

NB. dim_t bli_obj_row_off( const obj_t* obj );
obj_row_off=:                                                                    memr@(,&( OBJ_T_OFFS_OF          , 1 , BINT ))

NB. dim_t bli_obj_col_off( const obj_t* obj );
obj_col_off=:                                                                    memr@(,&((OBJ_T_OFFS_OF + SZB  ) , 1 , BINT ))

NB. C: dim_t bli_obj_off( mdim_t mdim, const obj_t* obj );
NB. J: off=. mdim obj_off obj
obj_off=:                                                                        memr@(, (1 , BINT) ,~ OBJ_T_OFFS_OF + SZB&*)~

NB. dim_t bli_obj_length( const obj_t* obj );
obj_length=:                                                                     memr@(,&((OBJ_T_OFFS_DM + SZB*M) , 1 , BINT ))

NB. dim_t bli_obj_width( const obj_t* obj );
obj_width=:                                                                      memr@(,&((OBJ_T_OFFS_DM + SZB*N) , 1 , BINT ))

NB. C: dim_t bli_obj_dim( mdim_t mdim, const obj_t* obj );
NB. J: dim=. mdim obj_dim obj
obj_dim=:                                                                        memr@(, (1 , BINT) ,~ OBJ_T_OFFS_DM + SZB&*)~

NB. dim_t bli_obj_min_dim( const obj_t* obj );
NB. dim=. obj_min_dim obj
obj_min_dim=: obj_length <. obj_width

NB. dim_t bli_obj_max_dim( const obj_t* obj );
NB. dim=. obj_max_dim obj
obj_max_dim=: obj_length >. obj_width

NB. dim_t bli_obj_length_after_trans( const obj_t* obj );
obj_length_after_trans=: obj_length`obj_width@.obj_has_trans

NB. dim_t bli_obj_width_after_trans( const obj_t* obj );
obj_width_after_trans=: obj_width`obj_length@.obj_has_trans

NB. doff_t bli_obj_diag_offset( const obj_t* obj );
obj_diag_offset=:                                                                memr@(,&( OBJ_T_OFFS_DO          , 1 , BINT ))

NB. doff_t bli_obj_diag_offset_after_trans( const obj_t* obj );
obj_diag_offset_after_trans=: -^:(obj_has_trans`obj_diag_offset)

NB. Element size query
NB. siz_t bli_obj_elem_size( const obj_t* obj );
obj_elem_size=:                                                                  memr@(,&( OBJ_T_OFFS_ES          , 1 , BINT ))

NB. void* bli_obj_buffer( const obj_t* obj );
obj_buffer=:                                                                     memr@(,&( OBJ_T_OFFS_BF          , 1 , JPTR ))

NB. inc_t bli_obj_row_stride( const obj_t* obj );
obj_row_stride=:                                                                 memr@(,&( OBJ_T_OFFS_RS          , 1 , BINT ))

NB. inc_t bli_obj_col_stride( const obj_t* obj );
obj_col_stride=:                                                                 memr@(,&( OBJ_T_OFFS_CS          , 1 , BINT ))

NB. inc_t bli_obj_imag_stride( const obj_t* obj );
obj_imag_stride=:                                                                memr@(,&( OBJ_T_OFFS_IS          , 1 , BINT ))

NB. void* bli_obj_internal_scalar_buffer( const obj_t* obj );
obj_internal_scalar_buffer=: OBJ_T_OFFS_SC&+

NB. dim_t bli_obj_padded_length( const obj_t* obj );
obj_padded_length=:                                                              memr@(,&( OBJ_T_OFFS_MD          , 1 , BINT ))

NB. dim_t bli_obj_padded_width( const obj_t* obj );
obj_padded_width=:                                                               memr@(,&( OBJ_T_OFFS_ND          , 1 , BINT ))

NB. inc_t bli_obj_panel_stride( const obj_t* obj );
obj_panel_stride=:                                                               memr@(,&( OBJ_T_OFFS_PS          , 1 , BINT ))

NB. inc_t bli_obj_panel_dim( const obj_t* obj );
obj_panel_dim=:                                                                  memr@(,&( OBJ_T_OFFS_PD          , 1 , BINT ))

NB. dim_t bli_obj_panel_length( const obj_t* obj );
obj_panel_length=:                                                               memr@(,&( OBJ_T_OFFS_ML          , 1 , BINT ))

NB. dim_t bli_obj_panel_width( const obj_t* obj );
obj_panel_width=:                                                                memr@(,&( OBJ_T_OFFS_NL          , 1 , BINT ))

NB. bool bli_obj_has_trans( const obj_t* obj );
obj_has_trans=: BITVAL_TRANS -: obj_onlytrans_status

NB. bool bli_obj_has_notrans( const obj_t* obj );
obj_has_notrans=: BITVAL_NO_TRANS -: obj_onlytrans_status

NB. bool bli_obj_has_conj( const obj_t* obj );
obj_has_conj=: BITVAL_CONJ -: obj_conj_status

NB. bool bli_obj_has_noconj( const obj_t* obj );
obj_has_noconj=: BITVAL_NO_CONJ -: obj_conj_status

NB. bool bli_obj_has_nonunit_diag( const obj_t* obj );
obj_has_nonunit_diag=: BITVAL_NONUNIT_DIAG -: obj_diag

NB. bool bli_obj_has_unit_diag( const obj_t* obj );
obj_has_unit_diag=: BITVAL_UNIT_DIAG -: bli_obj_diag

NB. dim_t bli_obj_length_after_trans( const obj_t* obj );
obj_length_after_trans=: obj_length`obj_width@.obj_has_trans

NB. dim_t bli_obj_width_after_trans( const obj_t* obj );
obj_width_after_trans=: obj_width`obj_length@.obj_has_trans

NB. dim_t bli_obj_vector_dim( const obj_t* obj );
obj_vector_dim=: obj_length`obj_width@.(1 -: obj_length)

NB. dim_t bli_obj_length_stored( const obj_t* obj );
obj_length_stored=: (bli_obj_length + 0 <.   obj_diag_offset)`(obj_length <. obj_width  - obj_diag_offset)@.obj_is_upper

NB. dim_t bli_obj_width_stored( const obj_t* obj );
obj_width_stored=:  (bli_obj_width  + 0 <. -@obj_diag_offset)`(obj_width  <. obj_length + obj_diag_offset)@.obj_is_lower

NB. dim_t bli_obj_length_stored_after_trans( const obj_t* obj );
obj_length_stored_after_trans=: obj_length_stored`obj_width_stored @.obj_has_trans

NB. dim_t bli_obj_width_stored_after_trans( const obj_t* obj );
obj_width_stored_after_trans=:  obj_width_stored `obj_length_stored@.obj_has_trans

NB. dim_t bli_obj_vector_dim( const obj_t* x );
obj_vector_dim=: obj_length`obj_width@.(1 = bli_obj_length)

NB. inc_t bli_obj_vector_inc( const obj_t* x );
obj_vector_inc=: obj_row_stride`obj_col_stride@.(1 = obj_length)`1:@.obj_is_1x1

NB. bool bli_obj_is_vector( const obj_t* x );
obj_is_vector=: obj_length +.&(1&=) obj_width

NB. bool bli_obj_is_row_vector( const obj_t* x );
obj_is_row_vector=: 1 = obj_length

NB. bool bli_obj_is_col_vector( const obj_t* x );
obj_is_col_vector=: 1 = obj_width

NB. bool bli_obj_has_zero_dim( const obj_t* x );
obj_has_zero_dim=: obj_length +.&(0&=) obj_width

NB. bool bli_obj_is_general( const obj_t* obj );
obj_is_general=: BITVAL_GENERAL -: obj_struc

NB. bool bli_obj_is_hermitian( const obj_t* obj );
obj_is_hermitian=: BITVAL_HERMITIAN -: obj_struc

NB. bool bli_obj_is_symmetric( const obj_t* obj );
obj_is_symmetric=: BITVAL_SYMMETRIC -: obj_struc

NB. bool bli_obj_is_triangular( const obj_t* obj );
obj_is_triangular=: BITVAL_TRIANGULAR -: obj_struc

NB. bool bli_obj_is_1x1( const obj_t* x );
obj_is_1x1=: obj_length *.&(1&=) obj_width

NB. inc_t bli_obj_vector_inc( const obj_t* obj );
obj_vector_inc=: obj_row_stride`obj_col_stride@.(1 -: obj_length)`1:@.obj_is_1x1

NB. bool bli_obj_is_upper( const obj_t* obj );
obj_is_upper=: BITVAL_UPPER -: obj_uplo

NB. bool bli_obj_is_lower( const obj_t* obj );
obj_is_lower=: BITVAL_LOWER -: obj_uplo

NB. bool bli_obj_is_upper_or_lower( const obj_t* obj );
obj_is_upper_or_lower=: obj_is_upper +. obj_is_lower

NB. bool bli_obj_is_dense( const obj_t* obj );
obj_is_dense=: BITVAL_DENSE -: obj_uplo

NB. bool bli_obj_is_zeros( const obj_t* obj );
obj_is_zeros=: BITVAL_ZEROS -: obj_uplo

NB. bool bli_obj_is_float( const obj_t* obj );
obj_is_float=: BITVAL_FLOAT_TYPE -: obj_dt

NB. bool bli_obj_is_double( const obj_t* obj );
obj_is_double=: BITVAL_DOUBLE_TYPE -: obj_dt

NB. bool bli_obj_is_scomplex( const obj_t* obj );
obj_is_scomplex=: BITVAL_SCOMPLEX_TYPE -: obj_dt

NB. bool bli_obj_is_dcomplex( const obj_t* obj );
obj_is_dcomplex=: BITVAL_DCOMPLEX_TYPE -: obj_dt

NB. bool bli_obj_is_int( const obj_t* obj );
obj_is_int=: BITVAL_INT_TYPE -: obj_dt

NB. bool bli_obj_is_const( const obj_t* obj );
obj_is_const=: BITVAL_CONST_TYPE -: obj_dt

NB. bool bli_obj_is_real( const obj_t* obj );
obj_is_real=: (BITVAL_REAL -: obj_domain) *. -.@obj_is_const

NB. bool bli_obj_is_complex( const obj_t* obj );
obj_is_complex=: (BITVAL_COMPLEX -: obj_domain) *. -.@obj_is_const

NB. bool bli_obj_is_single_prec( const obj_t* obj );
obj_is_single_prec=: BITVAL_SINGLE_PREC -: obj_prec

NB. bool bli_obj_is_double_prec( const obj_t* obj );
obj_is_double_prec=: BITVAL_DOUBLE_PREC -: obj_prec

NB. num_t bli_obj_dt_proj_to_single_prec( const obj_t* obj );
obj_dt_proj_to_single_prec=: (NOT BITVAL_SINGLE_PREC) AND obj_dt

NB. num_t bli_obj_dt_proj_to_double_prec( const obj_t* obj );
obj_dt_proj_to_double_prec=: (NOT BITVAL_DOUBLE_PREC) OR  obj_dt

NB. num_t bli_obj_dt_proj_to_real( const obj_t* obj );
obj_dt_proj_to_real=: (NOT BITVAL_COMPLEX) AND obj_dt

NB. num_t bli_obj_dt_proj_to_complex( const obj_t* obj );
obj_dt_proj_to_complex=: BITVAL_COMPLEX OR obj_dt

NB. Check if two objects are aliases of one another
NB. C: bool bli_obj_is_alias_of( const obj_t* a, const obj_t* b );
NB. J: bool=. a obj_is_alias_of b
obj_is_alias_of=: =&obj_buffer

NB. void* bli_obj_buffer_at_off( const obj_t* obj );
obj_buffer_at_off=: obj_buffer + obj_elem_size * (obj_col_off , obj_row_off) mp obj_col_stride , obj_row_stride

NB. inc_t bli_obj_row_stride_mag( const obj_t* obj );
obj_row_stride_mag=:  |@obj_row_stride

NB. inc_t bli_obj_col_stride_mag( const obj_t* obj );
obj_col_stride_mag=:  |@obj_col_stride

NB. inc_t bli_obj_imag_stride_mag( const obj_t* obj );
obj_imag_stride_mag=: |@obj_imag_stride

NB. bool bli_obj_root_is_general( const obj_t* obj );
obj_root_is_general=: obj_is_general@obj_root

NB. bool bli_obj_root_is_hermitian( const obj_t* obj );
obj_root_is_hermitian=: obj_is_hermitian@obj_root

NB. bool bli_obj_root_is_symmetric( const obj_t* obj );
obj_root_is_symmetric=: obj_is_symmetric@obj_root

NB. bool bli_obj_root_is_triangular( const obj_t* obj );
obj_root_is_triangular=: obj_is_triangular@obj_root

NB. bool bli_obj_root_is_herm_or_symm( const obj_t* obj );
obj_root_is_herm_or_symm=: (obj_is_hermitian OR obj_is_symmetric)@obj_root

NB. bool bli_obj_root_is_upper( const obj_t* obj );
obj_root_is_upper=: obj_is_upper@obj_root

NB. bool bli_obj_root_is_lower( const obj_t* obj );
obj_root_is_lower=: obj_is_lower@obj_root

NB. ---------------------------------------------------------
NB. Object mutator function reference

NB. info field

NB. C: void bli_obj_apply_trans( trans_t trans, obj_t* obj );
NB. J: trash=. trans obj_apply_trans obj
obj_apply_trans=:             XOR                                                          meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_apply_conj( conj_t conj, obj_t* obj );
NB. J: trash=. conj obj_apply_conj obj
obj_apply_conj=: obj_apply_trans

NB. C: void bli_obj_set_conjtrans( trans_t trans, obj_t* obj );
NB. J: trash=. trans obj_set_conjtrans obj
obj_set_conjtrans=:           ( OR                       (NOT CONJTRANS_BITS        )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_onlytrans( trans_t trans, obj_t* obj );
NB. J: trash=. trans obj_set_onlytrans obj
obj_set_onlytrans=:           ( OR                       (NOT TRANS_BIT             )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_conj( conj_t conj, obj_t* obj );
NB. J: trash=. conj obj_set_conj obj
obj_set_conj=:                ( OR                       (NOT CONJ_BIT              )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_struc( struc_t struc, obj_t* obj );
NB. J: trash=. struc obj_set_struc obj
obj_set_struc=:               ( OR                       (NOT STRUC_BITS            )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_uplo( uplo_t uplo, obj_t* obj );
NB. J: trash=. uplo obj_set_uplo obj
obj_set_uplo=:                ( OR                       (NOT UPLO_BITS             )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_diag( diag_t diag, obj_t* obj );
NB. J: trash=. diag obj_set_diag obj
obj_set_diag=:                ( OR                       (NOT UNIT_DIAG_BIT         )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_invert_diag( invdiag_t invdiag, obj_t* obj );
NB. J: trash=. invdiag obj_set_invert_diag obj
obj_set_invert_diag=:         ( OR                       (NOT INVERT_DIAG_BIT       )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_dt( num_t dt, obj_t* obj );
NB. J: trash=. dt obj_set_dt obj
obj_set_dt=:                  ( OR                       (NOT DATATYPE_BITS         )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_comp_prec( prec_t dt, obj_t* obj );
NB. J: trash=. dt obj_set_comp_prec obj
obj_set_comp_prec=:           ((OR COMP_PREC_SHIFT&SFT)~ (NOT COMP_PREC_BIT         )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. void bli_obj_toggle_trans( obj_t* obj );
obj_toggle_trans=: TRANSPOSE&obj_apply_trans

NB. void bli_obj_toggle_conj( obj_t* obj );
obj_toggle_conj=: CONJUGATE&obj_apply_conj

NB. C: void bli_obj_toggle_uplo( obj_t* obj );
NB. J: trash=. obj_toggle_uplo obj
NB. note: invert both lower and upper bits
obj_toggle_uplo=:                                        (LOWER_BIT OR UPPER_BIT    )&XOR  meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_pack_schema( pack_t schema, obj_t* obj );
NB. J: trash=. schema obj_set_pack_schema obj
obj_set_pack_schema=:         ( OR                       (NOT PACK_SCHEMA_BITS      )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_pack_order_if_upper( packord_t ordif, obj_t* obj );
NB. J: trash=. ordif obj_set_pack_order_if_upper obj
obj_set_pack_order_if_upper=: ( OR                       (NOT PACK_REV_IF_UPPER_BIT )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_pack_order_if_lower( packord_t ordif, obj_t* obj );
NB. J: trash=. ordif obj_set_pack_order_if_lower obj
obj_set_pack_order_if_lower=: ( OR                       (NOT PACK_REV_IF_LOWER_BIT )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. C: void bli_obj_set_pack_buffer_type( packbuf_t buf_type, obj_t* obj );
NB. J: trash=. buf_type obj_set_pack_buffer_type obj
obj_set_pack_buffer_type=:    ( OR                       (NOT PACK_BUFFER_BITS      )&AND) meme ,&( OBJ_T_OFFS_IN , 1 , JINT4)

NB. info2 field

NB. C: void bli_obj_set_scalar_dt( num_t dt, obj_t* obj );
NB. J: trash=. dt obj_set_scalar_dt obj
obj_set_scalar_dt=:           ((OR SCALAR_DT_SHIFT&SFT)~ (NOT SCALAR_DT_BITS        )&AND) meme ,&( OBJ_T_OFFS_I2 , 1 , JINT4)

NB. C: void bli_obj_set_scalar_domain( dom_t dt, obj_t* obj );
NB. J: trash=. dt obj_set_scalar_domain obj
obj_set_scalar_domain=:       ((OR SCALAR_DT_SHIFT&SFT)~ (NOT SCALAR_DOMAIN_BIT     )&AND) meme ,&( OBJ_T_OFFS_I2 , 1 , JINT4)

NB. C: void bli_obj_set_scalar_prec( prec_t dt, obj_t* obj );
NB. J: trash=. dt obj_set_scalar_prec obj
obj_set_scalar_prec=:         ((OR SCALAR_DT_SHIFT&SFT)~ (NOT SCALAR_PREC_BIT       )&AND) meme ,&( OBJ_T_OFFS_I2 , 1 , JINT4)

NB. other fields

NB. void bli_obj_set_as_root( obj_t* obj );
obj_set_as_root=:                                   memw ,&( 0                        , 1 , JPTR )

NB. C: void bli_obj_set_off( mdim_t mdim, dim_t offset, obj_t* obj );
NB. J: trash=. (mdim , offset) obj_set_off obj
obj_set_off=:                               (1 { [) memw (1 , BINT) ,~ (, OBJ_T_OFFS_OF + SZB * {.)~

NB. C: void bli_obj_set_offs( dim_t offm, dim_t offn, obj_t* obj );
NB. J: trash=. (offm , offn) obj_set_offs obj
obj_set_offs=:                                      memw ,&( OBJ_T_OFFS_OF            , 2 , BINT )

NB. C: void bli_obj_inc_off( mdim_t mdim, dim_t offset, obj_t* obj );
NB. J: trash=. (mdim , offset) obj_inc_off obj
obj_inc_off=:            (1 { [) + meme (1 , BINT) ,~ (, OBJ_T_OFFS_OF + SZB * {.)~

NB. C: void bli_obj_inc_offs( dim_t offm, dim_t offn, obj_t* obj );
NB. J: trash=. (offm , offn) obj_inc_offs obj
obj_inc_offs=:                   + meme ,&(OBJ_T_OFFS_OF , 2 , BINT )

NB. C: void bli_obj_set_length( dim_t m, obj_t* obj );
NB. J: trash=. m obj_set_length obj
obj_set_length=:                                    memw ,&((OBJ_T_OFFS_DM + SZB * M) , 1 , BINT )

NB. C: void bli_obj_set_width( dim_t n, obj_t* obj );
NB. J: trash=. n obj_set_width obj
obj_set_width=:                                     memw ,&((OBJ_T_OFFS_DM + SZB * N) , 1 , BINT )

NB. C: void bli_obj_set_dim( mdim_t mdim, dim_t dim_val, obj_t* obj );
NB. J: trash=. (mdim , dim_val) obj_set_dim obj
obj_set_dim=:                               (1 { [) memw (1 , BINT) ,~ (, OBJ_T_OFFS_DM + SZB * {.)~

NB. C: void bli_obj_set_dims( dim_t m, dim_t n, obj_t* obj );
NB. J: trash=. (m , n) obj_set_dims obj
obj_set_dims=:                                      memw ,&( OBJ_T_OFFS_DM            , 2 , BINT )

NB. C: void bli_obj_set_dims_with_trans( trans_t trans, dim_t m, dim_t n, obj_t* obj );
NB. J: trash=. (trans , m , n) obj_set_dims_with_trans obj
obj_set_dims_with_trans=: (obj_set_dims~ |.^:((BITVAL_TRANS ~: TRANS_BIT AND NOT@{.)`}.))~

NB. C: void bli_obj_set_diag_offset( doff_t doff, obj_t* obj );
NB. J: trash=. doff obj_set_diag_offset obj
obj_set_diag_offset=:                               memw ,&( OBJ_T_OFFS_DO            , 1 , BINT )

NB. void bli_obj_negate_diag_offset( obj_t* obj );
obj_negate_diag_offset=:         - meme ,&(OBJ_T_OFFS_DO , 1 , BINT )

NB. void bli_obj_inc_diag_offset( doff_t offset, obj_t* obj );
obj_inc_diag_offset=:            + meme ,&(OBJ_T_OFFS_DO , 1 , BINT )

NB. Element size modification
NB. C: void bli_obj_set_elem_size( siz_t size, obj_t* obj );
NB. J: trash=. size obj_set_elem_size obj
obj_set_elem_size=:                                 memw ,&( OBJ_T_OFFS_ES            , 1 , BINT )

NB. C: void bli_obj_set_buffer( void* p, obj_t* obj );
NB. J: trash=. p obj_set_buffer obj
obj_set_buffer=:                                    memw ,&( OBJ_T_OFFS_BF            , 1 , JPTR )

NB. C: void bli_obj_set_row_stride( inc_t rs, obj_t* obj );
NB. J: trash=. rs obj_set_row_stride obj
obj_set_row_stride=:                                memw ,&( OBJ_T_OFFS_RS            , 1 , BINT )

NB. C: void bli_obj_set_col_stride( inc_t cs, obj_t* obj );
NB. J: trash=. cs obj_set_col_stride obj
obj_set_col_stride=:                                memw ,&( OBJ_T_OFFS_CS            , 1 , BINT )

NB. C: void bli_obj_set_imag_stride( inc_t is, obj_t* obj );
NB. J: trash=. is obj_set_imag_stride obj
obj_set_imag_stride=:                               memw ,&( OBJ_T_OFFS_IS            , 1 , BINT )

NB. Bufferless scalar field modification
NB. C: void bli_obj_copy_internal_scalar( const obj_t* a, obj_t* b );
NB. J: trash=. a obj_copy_internal_scalar b
obj_copy_internal_scalar=: (memw~ memr)~&(,(OBJ_T_OFFS_SC , 16 , JCHAR))

NB. C: void bli_obj_set_padded_length( dim_t m, obj_t* obj );
NB. J: trash=. m obj_set_padded_length obj
obj_set_padded_length=:                             memw ,&( OBJ_T_OFFS_MD            , 1 , BINT )

NB. C: void bli_obj_set_padded_width( dim_t n, obj_t* obj );
NB. J: trash=. n obj_set_padded_width obj
obj_set_padded_width=:                              memw ,&( OBJ_T_OFFS_ND            , 1 , BINT )

NB. C: void bli_obj_set_panel_stride( inc_t ps, obj_t* obj );
NB. J: trash=. ps obj_set_panel_stride obj
obj_set_panel_stride=:                              memw ,&( OBJ_T_OFFS_PS            , 1 , BINT )

NB. C: void bli_obj_set_panel_dim( inc_t pd, obj_t* obj );
NB. J: trash=. pd obj_set_panel_dim obj
obj_set_panel_dim=:                                 memw ,&( OBJ_T_OFFS_PD            , 1 , BINT )

NB. C: void bli_obj_set_panel_length( dim_t m, obj_t* obj );
NB. J: trash=. m obj_set_panel_length obj
obj_set_panel_length=:                              memw ,&( OBJ_T_OFFS_ML            , 1 , BINT )

NB. C: void bli_obj_set_panel_width( dim_t n, obj_t* obj );
NB. J: trash=. n obj_set_panel_width obj
obj_set_panel_width=:                               memw ,&( OBJ_T_OFFS_NL            , 1 , BINT )

NB. C: void bli_obj_set_strides( inc_t rs, inc_t cs, obj_t* obj );
NB. J: trash=. (rs , cs) obj_set_strides obj
obj_set_strides=:     obj_set_row_stride   `obj_set_col_stride  "0

NB. C: void bli_obj_set_padded_dims( dim_t m, dim_t n, obj_t* obj );
NB. J: trash=. (m , n) obj_set_padded_dims obj
obj_set_padded_dims=: obj_set_padded_length`obj_set_padded_width"0

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

NB. Make a full alias (shallow copy)
NB. C: void bli_obj_alias_to( const obj_t* a, obj_t* b );
NB. J: trash=. a obj_alias_to b
obj_alias_to=: obj_init_full_shallow_copy_of

NB. Create an alias with a trans value applied
NB. C: void bli_obj_alias_with_trans( trans_t trans, const obj_t* a, obj_t* b );
NB. J: trash=. trans obj_alias_with_trans a ; b
NB. note: trans may include a conj component
obj_alias_with_trans=: (obj_apply_trans 1&{::) obj_alias_to&>/

NB. Create an alias with a conj value applied
NB. C: void bli_obj_alias_with_conj( conj_t conja, const obj_t* a, obj_t* b );
NB. J: trash=. conja obj_alias_with_conj a ; b
obj_alias_with_conj=:  (obj_apply_conj  1&{::) obj_alias_to&>/

NB. Alias only the real part
NB. C: void bli_obj_real_part( const obj_t* c, obj_t* r );
NB. J: trash=. c obj_real_part r
obj_real_part=: 4 : 0
  x obj_alias_to y

  if. obj_is_complex x do.
    NB. Change the datatypes.
    y obj_set_dt~ obj_dt_proj_to_real x

    NB. Don't touch the attached scalar datatype.

    NB. Update the element size.
    y obj_set_elem_size~ -: obj_elem_size x

    NB. Update the strides.
    y obj_set_strides~ +: (obj_row_stride , obj_col_stride) x

    NB. Buffer is left unchanged.
  end.
)

NB. Alias only the imaginary part
NB. C: void bli_obj_imag_part( const obj_t* c, obj_t* i );
NB. J: trash=. c obj_imag_part i
obj_imag_part=: 4 : 0
  if. obj_is_complex x do.
    x obj_alias_to y

    NB. Change the datatype.
    y obj_set_dt~ obj_dt_proj_to_real x

    NB. Don't touch the attached scalar datatype.

    NB. Update the element size.
    es_c=. obj_elem_size x
    (-: es_c) obj_set_elem_size y

    NB. Update the strides.
    y obj_set_strides~ +: (obj_row_stride , obj_col_stride) x

    NB. Update the buffer.
    is_c=. obj_imag_stride x
    p=.    obj_buffer      x
    (p + is_c * -: es_c) obj_set_buffer y
  end.
)

NB. void bli_obj_print
NB.      (
NB.        const char*  label,
NB.        const obj_t* obj
NB.      );

obj_print=: (lib,' bli_obj_print ',ifw,' n &c &')&cd


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
NB.        const obj_t*  x,
NB.        const obj_t*  y,
NB.      );

addv_cd=: (lib,' bli_addv ',ifw,' n & *')&cd

NB. void bli_amaxv
NB.      (
NB.        const obj_t*  x,
NB.        const obj_t*  index
NB.      );

amaxv_cd=: (lib,' bli_amaxv ',ifw,' n & *')&cd

NB. y := y + conj?(alpha) * conj?(x)
NB.
NB. void bli_axpyv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x,
NB.        const obj_t*  y
NB.      );

axpyv_cd=: (lib,' bli_axpyv ',ifw,' n & & *')&cd

NB. y := conj?(beta) * y + conj?(alpha) * conj?(x)
NB.
NB. void bli_axpbyv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x,
NB.        const obj_t*  beta,
NB.        const obj_t*  y
NB.      );

axpbyv_cd=: (lib,' bli_axpbyv ',ifw,' n & & & *')&cd

NB. y := conj?(x)
NB.
NB. void bli_copyv
NB.      (
NB.        const obj_t*  x,
NB.        const obj_t*  y
NB.      );

copyv_cd=: (lib,' bli_copyv ',ifw,' n & *')&cd

NB. rho := conj?(x)^T * conj?(y)
NB.
NB. void bli_dotv
NB.      (
NB.        const obj_t*  x,
NB.        const obj_t*  y,
NB.        const obj_t*  rho
NB.      );

dotv_cd=: (lib,' bli_dotv ',ifw,' n & & *')&cd

NB. rho := conj?(beta) * rho + conj?(alpha) * conj?(x)^T * conj?(y)
NB.
NB. void bli_dotxv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x,
NB.        const obj_t*  y,
NB.        const obj_t*  beta,
NB.        const obj_t*  rho
NB.      );

dotxv_cd=: (lib,' bli_dotxv ',ifw,' n & & & & *')&cd

NB. void bli_invertv
NB.      (
NB.        const obj_t*  x
NB.      );

invertv_cd=: (lib,' bli_invertv ',ifw,' n *')&cd

NB. x := ( 1.0 / conj?(alpha) ) * x
NB.
NB. void bli_invscalv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x
NB.      );

invscalv_cd=: (lib,' bli_invscalv ',ifw,' n & *')&cd

NB. x := conj?(alpha) * x
NB.
NB. void bli_scalv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x
NB.      );

scalv_cd=: (lib,' bli_scalv ',ifw,' n & *')&cd

NB. y := conj?(alpha) * conj?(x)
NB.
NB. void bli_scal2v
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x,
NB.        const obj_t*  y
NB.      );

scal2v_cd=: (lib,' bli_scal2v ',ifw,' n & & *')&cd

NB. x := conj?(alpha)
NB.
NB. void bli_setv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x
NB.      );

setv_cd=: (lib,' bli_setv ',ifw,' n & *')&cd

NB. real(x) := real(alpha)
NB.
NB. void bli_setrv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x
NB.      );

setrv_cd=: (lib,' bli_setrv ',ifw,' n & *')&cd

NB. imag(x) := real(alpha)
NB.
NB. void bli_setiv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x
NB.      );

setiv_cd=: (lib,' bli_setiv ',ifw,' n & *')&cd

NB. y := y - conj?(x)
NB.
NB. void bli_subv
NB.      (
NB.        const obj_t*  x,
NB.        const obj_t*  y
NB.      );

subv_cd=: (lib,' bli_subv ',ifw,' n & *')&cd

NB. void bli_swapv
NB.      (
NB.        const obj_t*  x,
NB.        const obj_t*  y
NB.      );

swapv_cd=: (lib,' bli_swapv ',ifw,' n * *')&cd

NB. y := conj?(beta) * y + conj?(x)
NB.
NB. void bli_xpbyv
NB.      (
NB.        const obj_t*  x,
NB.        const obj_t*  beta,
NB.        const obj_t*  y
NB.      );

xpbyv_cd=: (lib,' bli_xpbyv ',ifw,' n & & *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1d operations

NB. void bli_addd
NB.      (
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

addd_cd=: (lib,' bli_addd ',ifw,' n & *')&cd

NB. void bli_axpyd
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

axpyd_cd=: (lib,' bli_axpyd ',ifw,' n & & *')&cd

NB. void bli_copyd
NB.      (
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

copyd_cd=: (lib,' bli_copyd ',ifw,' n & *')&cd

NB. void bli_invertd
NB.      (
NB.        const obj_t*  a
NB.      );

invertd_cd=: (lib,' bli_invertd ',ifw,' n *')&cd

NB. void bli_invscald
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a
NB.      );

invscald_cd=: (lib,' bli_invscald ',ifw,' n & *')&cd

NB. void bli_scald
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a
NB.      );

scald_cd=: (lib,' bli_scald ',ifw,' n & *')&cd

NB. void bli_scal2d
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

scal2d_cd=: (lib,' bli_scal2d ',ifw,' n & & *')&cd

NB. void bli_setd
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a
NB.      );

setd_cd=: (lib,' bli_setd ',ifw,' n & *')&cd

NB. void bli_setid
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a
NB.      );

setid_cd=: (lib,' bli_setid ',ifw,' n & *')&cd

NB. void bli_shiftd
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a
NB.      );

shiftd_cd=: (lib,' bli_shiftd ',ifw,' n & *')&cd

NB. void bli_subd
NB.      (
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

subd_cd=: (lib,' bli_subd ',ifw,' n & *')&cd

NB. void bli_xpbyd
NB.      (
NB.        const obj_t*  a,
NB.        const obj_t*  beta,
NB.        const obj_t*  b
NB.      );

xpbyd_cd=: (lib,' bli_xpbyd ',ifw,' n & & *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1m operations

NB. B := B + trans?(A)
NB.
NB. void bli_addm
NB.      (
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

addm_cd=: (lib,' bli_addm ',ifw,' n & *')&cd

NB. B := B + conj?(alpha) * trans?(A)
NB.
NB. void bli_axpym
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

axpym_cd=: (lib,' bli_axpym ',ifw,' n & & *')&cd

NB. B := trans?(A)
NB.
NB. void bli_copym
NB.      (
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

copym_cd=: (lib,' bli_copym ',ifw,' n & *')&cd

NB. A := ( 1.0 / conj?(alpha) ) * A
NB.
NB. void bli_invscalm
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a
NB.      );

invscalm_cd=: (lib,' bli_invscalm ',ifw,' n & *')&cd

NB. A := conj?(alpha) * A
NB.
NB. void bli_scalm
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a
NB.      );

scalm_cd=: (lib,' bli_scalm ',ifw,' n & *')&cd

NB. B := conj?(alpha) * trans?(A)
NB.
NB. void bli_scal2m
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

scal2m_cd=: (lib,' bli_scal2m ',ifw,' n & & *')&cd

NB. A := conj?(alpha)
NB.
NB. void bli_setm
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a
NB.      );

setm_cd=: (lib,' bli_setm ',ifw,' n & *')&cd

NB. real(A) := real(alpha)
NB.
NB. void bli_setrm
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a
NB.      );

setrm_cd=: (lib,' bli_setrm ',ifw,' n & *')&cd

NB. imag(A) := real(alpha)
NB.
NB. void bli_setim
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a
NB.      );

setim_cd=: (lib,' bli_setim ',ifw,' n & *')&cd

NB. B := B - trans?(A)
NB.
NB. void bli_subm
NB.      (
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

subm_cd=: (lib,' bli_subm ',ifw,' n & *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1f operations

NB. z := z + conj?(alphax) * conj?(x) + conj?(alphay) * conj?(y)
NB.
NB. void bli_axpy2v
NB.      (
NB.        const obj_t*  alphax,
NB.        const obj_t*  alphay,
NB.        const obj_t*  x,
NB.        const obj_t*  y,
NB.        const obj_t*  z
NB.      );

axpy2v_cd=: (lib,' bli_axpy2v ',ifw,' n & & & & *')&cd

NB. rho := conj?(x)^T * conj?(y)
NB. z   := z + conj?(alpha) * conj?(x)
NB.
NB. void bli_dotaxpyv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  xt,
NB.        const obj_t*  x,
NB.        const obj_t*  y,
NB.        const obj_t*  rho,
NB.        const obj_t*  z
NB.      );

dotaxpyv_cd=: (lib,' bli_dotaxpyv ',ifw,' n & & & & * *')&cd

NB. y := y + alpha * conja(A) * conjx(x)
NB.
NB. void bli_axpyf
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  x,
NB.        const obj_t*  y
NB.      );

axpyf_cd=: (lib,' bli_axpyf ',ifw,' n & & & *')&cd

NB. y := conj?(beta) * y + conj?(alpha) * conj?(A)^T * conj?(x)
NB.
NB. void bli_dotxf
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  x,
NB.        const obj_t*  beta,
NB.        const obj_t*  y
NB.      );

dotxf_cd=: (lib,' bli_dotxf ',ifw,' n & & & & *')&cd

NB. y := conj?(beta) * y + conj?(alpha) * conj?(A)^T * conj?(w)
NB. z :=               z + conj?(alpha) * conj?(A)   * conj?(x)
NB.
NB. void bli_dotxaxpyf
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  at,
NB.        const obj_t*  a,
NB.        const obj_t*  w,
NB.        const obj_t*  x,
NB.        const obj_t*  beta,
NB.        const obj_t*  y,
NB.        const obj_t*  z
NB.      );

dotxaxpyf_cd=: (lib,' bli_dotxaxpyf ',ifw,' n & & & & & & * *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-2 operations

NB. y := conj?(beta) * y + conj?(alpha) * trans?(A) * conj?(x)
NB.
NB. void bli_gemv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  x,
NB.        const obj_t*  beta,
NB.        const obj_t*  y
NB.      );

gemv_cd=: (lib,' bli_gemv ',ifw,' n & & & & *')&cd

NB. A := A + conj?(alpha) * conj?(x) * conj?(y)^T
NB.
NB. void bli_ger
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x,
NB.        const obj_t*  y,
NB.        const obj_t*  a
NB.      );

ger_cd=: (lib,' bli_ger ',ifw,' n & & & *')&cd

NB. y := conj?(beta) * y + conj?(alpha) * conj?(A) * conj?(x)
NB.
NB. void bli_hemv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  x,
NB.        const obj_t*  beta,
NB.        const obj_t*  y
NB.      );

hemv_cd=: (lib,' bli_hemv ',ifw,' n & & & & *')&cd

NB. A := A + conj?(alpha) * conj?(x) * conj?(x)^H
NB.
NB. void bli_her
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x,
NB.        const obj_t*  a
NB.      );

her_cd=: (lib,' bli_her ',ifw,' n & & *')&cd

NB. A := A + alpha * conj?(x) * conj?(y)^H + conj(alpha) * conj?(y) * conj?(x)^H
NB.
NB. void bli_her2
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x,
NB.        const obj_t*  y,
NB.        const obj_t*  a
NB.      );

her2_cd=: (lib,' bli_her2 ',ifw,' n & & & *')&cd

NB. y := conj?(beta) * y + conj?(alpha) * conj?(A) * conj?(x)
NB.
NB. void bli_symv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  x,
NB.        const obj_t*  beta,
NB.        const obj_t*  y
NB.      );

symv_cd=: (lib,' bli_symv ',ifw,' n & & & & *')&cd

NB. A := A + conj?(alpha) * conj?(x) * conj?(x)^T
NB.
NB. void bli_syr
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x,
NB.        const obj_t*  a
NB.      );

syr_cd=: (lib,' bli_syr ',ifw,' n & & *')&cd

NB. A := A + alpha * conj?(x) * conj?(y)^T + conj(alpha) * conj?(y) * conj?(x)^T
NB.
NB. void bli_syr2
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  x,
NB.        const obj_t*  y,
NB.        const obj_t*  a
NB.      );

syr2_cd=: (lib,' bli_syr2 ',ifw,' n & & & *')&cd

NB. x := conj?(alpha) * transa(A) * x
NB.
NB. void bli_trmv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  x
NB.      );

trmv_cd=: (lib,' bli_trmv ',ifw,' n & & *')&cd

NB. Solve the linear system
NB. transa(A) * x = alpha * y
NB.
NB. void bli_trsv
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  y
NB.      );

trsv_cd=: (lib,' bli_trsv ',ifw,' n & & *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-3 operations

NB. C := beta * C + alpha * trans?(A) * trans?(B)
NB.
NB. void bli_gemm
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b,
NB.        const obj_t*  beta,
NB.        const obj_t*  c,
NB.      );

gemm_cd=: (lib,' bli_gemm ',ifw,' n & & & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(B)
NB.
NB. void bli_gemmt
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b,
NB.        const obj_t*  beta,
NB.        const obj_t*  c
NB.      );

gemmt_cd=: (lib,' bli_gemmt ',ifw,' n & & & & *')&cd

NB. C := beta * C + alpha * conj?(A) * trans?(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * trans?(B) * conj?(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_hemm
NB.      (
NB.              side_t  sidea,
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b,
NB.        const obj_t*  beta,
NB.        const obj_t*  c
NB.      );

hemm_cd=: (lib,' bli_hemm ',ifw,' n i & & & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(A)^H
NB.
NB. void bli_herk
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  beta,
NB.        const obj_t*  c
NB.      );

herk_cd=: (lib,' bli_herk ',ifw,' n & & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(B)^H + conj(alpha) * trans?(B) * trans?(A)^H
NB.
NB. void bli_her2k
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b,
NB.        const obj_t*  beta,
NB.        const obj_t*  c
NB.      );

her2k_cd=: (lib,' bli_her2k ',ifw,' n & & & & *')&cd

NB. C := beta * C + alpha * conj?(A) * trans?(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * trans?(B) * conj?(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_symm
NB.      (
NB.              side_t  sidea,
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b,
NB.        const obj_t*  beta,
NB.        const obj_t*  c
NB.      );

symm_cd=: (lib,' bli_symm ',ifw,' n i & & & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(A)^T
NB.
NB. void bli_syrk
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  beta,
NB.        const obj_t*  c
NB.      );

syrk_cd=: (lib,' bli_syrk ',ifw,' n & & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(B)^T + alpha * trans?(B) * trans?(A)^T
NB.
NB. void bli_syr2k
NB.      (
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b,
NB.        const obj_t*  beta,
NB.        const obj_t*  c
NB.      );

syr2k_cd=: (lib,' bli_syr2k ',ifw,' n & & & & *')&cd

NB. B := alpha * transa(A) * B  if sidea is BLIS_LEFT, or
NB. B := alpha * B * transa(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_trmm
NB.      (
NB.              side_t  sidea,
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

trmm_cd=: (lib,' bli_trmm ',ifw,' n i & & *')&cd

NB. C := beta * C + alpha * trans?(A) * trans?(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * trans?(B) * trans?(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_trmm3
NB.      (
NB.              side_t  sidea,
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b,
NB.        const obj_t*  beta,
NB.        const obj_t*  c
NB.      );

trmm3_cd=: (lib,' bli_trmm3 ',ifw,' n i & & & & *')&cd

NB. Solve the linear system with multiple right-hand sides
NB.   transa(A) * X = alpha * B  if sidea is BLIS_LEFT, or
NB.   X * transa(A) = alpha * B  if sidea is BLIS_RIGHT
NB.
NB. void bli_trsm
NB.      (
NB.              side_t  sidea,
NB.        const obj_t*  alpha,
NB.        const obj_t*  a,
NB.        const obj_t*  b
NB.      );

trsm_cd=: (lib,' bli_trsm ',ifw,' n i & & *')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Utility operations

NB. void bli_asumv
NB.      (
NB.        const obj_t*  x,
NB.        const obj_t*  asum
NB.      );

asumv_cd=: (lib,' bli_asumv ',ifw,' n & *')&cd

NB. void bli_norm[1fi]m
NB.      (
NB.        const obj_t*  a,
NB.        const obj_t*  norm
NB.      );

norm1m_cd=: (lib,' bli_norm1m ',ifw,' n & *')&cd
normfm_cd=: (lib,' bli_normfm ',ifw,' n & *')&cd
normim_cd=: (lib,' bli_normim ',ifw,' n & *')&cd

NB. void bli_norm[1fi]v
NB.      (
NB.        const obj_t*  x,
NB.        const obj_t*  norm
NB.      );

norm1v_cd=: (lib,' bli_norm1v ',ifw,' n & *')&cd
normfv_cd=: (lib,' bli_normfv ',ifw,' n & *')&cd
normiv_cd=: (lib,' bli_normiv ',ifw,' n & *')&cd

NB. void bli_mkherm
NB.      (
NB.        const obj_t*  a
NB.      );

mkherm_cd=: (lib,' bli_mkherm ',ifw,' n *')&cd

NB. void bli_mksymm
NB.      (
NB.        const obj_t*  a
NB.      );

mksymm_cd=: (lib,' bli_mksymm ',ifw,' n *')&cd

NB. void bli_mktrim
NB.      (
NB.        const obj_t*  a
NB.      );

mktrim_cd=: (lib,' bli_mktrim ',ifw,' n *')&cd

NB. void bli_fprintv
NB.      (
NB.              FILE*   file,
NB.        const char*   s1,
NB.        const obj_t*  x,
NB.        const char*   format,
NB.        const char*   s2
NB.      );

fprintv_cd=: (lib,' bli_fprintv ',ifw,' n & &c & &c &c')&cd

NB. void bli_fprintm
NB.      (
NB.              FILE*   file,
NB.        const char*   s1,
NB.        const obj_t*  a,
NB.        const char*   format,
NB.        const char*   s2
NB.      );

fprintm_cd=: (lib,' bli_fprintm ',ifw,' n & &c & &c &c')&cd

NB. void bli_printv
NB.      (
NB.        const char*   s1,
NB.        const obj_t*  x,
NB.        const char*   format,
NB.        const char*   s2
NB.      );

printv_cd=: (lib,' bli_printv ',ifw,' n &c & &c &c')&cd

NB. void bli_printm
NB.      (
NB.        const char*   s1,
NB.        const obj_t*  a,
NB.        const char*   format,
NB.        const char*   s2
NB.      );

printm_cd=: (lib,' bli_printm ',ifw,' n &c & &c &c')&cd

NB. void bli_randv
NB.      (
NB.        const obj_t*  x
NB.      );

randv_cd=: (lib,' bli_randv ',ifw,' n *')&cd

NB. void bli_randm
NB.      (
NB.        const obj_t*  a
NB.      );

randm_cd=: (lib,' bli_randm ',ifw,' n *')&cd

NB. scale_new^2 * sumsq_new = x[0]^2 + x[1]^2 + ... x[m-1]^2 + scale_old^2 * sumsq_old
NB.
NB. void bli_sumsqv
NB.      (
NB.        const obj_t*  x,
NB.        const obj_t*  scale,
NB.        const obj_t*  sumsq
NB.      );

sumsqv_cd=: (lib,' bli_sumsqv ',ifw,' n & * *')&cd

NB. void bli_getsc
NB.      (
NB.        const obj_t*   chi,
NB.              double*  zeta_r,
NB.              double*  zeta_i
NB.      );

getsc_cd=: (lib,' bli_getsc ',ifw,' n & * *')&cd

NB. err_t bli_getijv
NB.      (
NB.              dim_t    i,
NB.        const obj_t*   x,
NB.              double*  ar,
NB.              double*  ai
NB.      );

getijv_cd=: (lib,' bli_getijv ',ifw,' i ',BIC,' & *d *d')&cd

NB. err_t bli_getijm
NB.      (
NB.              dim_t    i,
NB.              dim_t    j,
NB.        const obj_t*   b,
NB.              double*  ar,
NB.              double*  ai
NB.      );

getijm_cd=: (lib,' bli_getijm ',ifw,' i ',BIC,' ',BIC,' & *d *d')&cd

NB. void bli_setsc
NB.      (
NB.              double  zeta_r,
NB.              double  zeta_i,
NB.        const obj_t*  chi
NB.      );

setsc_cd=: (lib,' bli_setsc ',ifw,' n d d *')&cd

NB. err_t bli_setijv
NB.      (
NB.              double  ar,
NB.              double  ai,
NB.              dim_t   i,
NB.        const obj_t*  x
NB.      );

setijv_cd=: (lib,' bli_setijv ',ifw,' i d d ',BIC,' *')&cd

NB. err_t bli_setijm
NB.      (
NB.              double  ar,
NB.              double  ai,
NB.              dim_t   i,
NB.              dim_t   j,
NB.        const obj_t*  b
NB.      );

setijm_cd=: (lib,' bli_setijm ',ifw,' i d d ',BIC,' ',BIC,' *')&cd

NB. void bli_eqsc
NB.      (
NB.        const obj_t*  chi,
NB.        const obj_t*  psi,
NB.              bool*   is_eq
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
NB.      EMPTY [ eqsc_cd_mtbli_ <"0 ({:L:0 objs) , < bscdat
NB.        NB. call bli_eqsc with arguments - pointers to data
NB.      a. i. bsc                    NB. cast the result to integer
NB.   1                               NB. value is changed to some value within boolean domain
NB.      objf_mtbli_ L: 0 objs        NB. destroy allocated objects
NB.      ({: a.) memw bscdat , 0 1    NB. reset to prepare to re-use in another comparison, optionally

eqsc_cd=: (lib,' bli_eqsc ',ifw,' n & & *b')&cd

NB. void bli_eqv
NB.      (
NB.        const obj_t*  x,
NB.        const obj_t*  y,
NB.              bool*   is_eq
NB.      );

eqv_cd=: (lib,' bli_eqv ',ifw,' n & & *b')&cd

NB. void bli_eqm
NB.      (
NB.        const obj_t*  a,
NB.        const obj_t*  b,
NB.              bool*   is_eq
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
NB.              conj_t  conjx,
NB.              dim_t   n,
NB.        const ctype*  x, inc_t incx,
NB.              ctype*  y, inc_t incy
NB.      );

daddv_cd=: (lib,' bli_daddv ',ifw,' n i ',BIC,' &d ',BIC,' *d ',BIC)&cd
zaddv_cd=: (lib,' bli_zaddv ',ifw,' n i ',BIC,' &j ',BIC,' *j ',BIC)&cd

NB. mimic BLAS routines i?amax()
NB.
NB. void bli_?amaxv
NB.      (
NB.              dim_t   n,
NB.        const ctype*  x, inc_t incx,
NB.              dim_t*  index
NB.      );

damaxv_cd=: (lib,' bli_damaxv ',ifw,' n ',BIC,' &d ',BIC,' *',BIC)&cd
zamaxv_cd=: (lib,' bli_zamaxv ',ifw,' n ',BIC,' &j ',BIC,' *',BIC)&cd

NB. y := y + alpha * conjx(x)
NB.
NB. void bli_?axpyv
NB.      (
NB.              conj_t  conjx,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.        const ctype*  x, inc_t incx,
NB.              ctype*  y, inc_t incy
NB.      );

daxpyv_cd=: (lib,' bli_daxpyv ',ifw,' n i ',BIC,' &d &d ',BIC,' *d ',BIC)&cd
zaxpyv_cd=: (lib,' bli_zaxpyv ',ifw,' n i ',BIC,' &j &j ',BIC,' *j ',BIC)&cd

NB. y := beta * y + alpha * conjx(x)
NB.
NB. void bli_?axpbyv
NB.      (
NB.              conj_t  conjx,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  beta,
NB.              ctype*  y, inc_t incy
NB.      );

daxpbyv_cd=: (lib,' bli_daxpbyv ',ifw,' n i ',BIC,' &d &d ',BIC,' &d *d ',BIC)&cd
zaxpbyv_cd=: (lib,' bli_zaxpbyv ',ifw,' n i ',BIC,' &j &j ',BIC,' &j *j ',BIC)&cd

NB. y := conjx(x)
NB.
NB. void bli_?copyv
NB.      (
NB.              conj_t  conjx,
NB.              dim_t   n,
NB.        const ctype*  x, inc_t incx,
NB.              ctype*  y, inc_t incy
NB.      );

dcopyv_cd=: (lib,' bli_dcopyv ',ifw,' n i ',BIC,' &d ',BIC,' *d ',BIC)&cd
zcopyv_cd=: (lib,' bli_zcopyv ',ifw,' n i ',BIC,' &j ',BIC,' *j ',BIC)&cd

NB. rho := conjx(x)^T * conjy(y)
NB.
NB. void bli_?dotv
NB.      (
NB.              conj_t  conjx,
NB.              conj_t  conjy,
NB.              dim_t   n,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  y, inc_t incy,
NB.              ctype*  rho
NB.      );

ddotv_cd=: (lib,' bli_ddotv ',ifw,' n i i ',BIC,' &d ',BIC,' &d ',BIC,' *d')&cd
zdotv_cd=: (lib,' bli_zdotv ',ifw,' n i i ',BIC,' &j ',BIC,' &j ',BIC,' *j')&cd

NB. rho := beta * rho + alpha * conjx(x)^T * conjy(y)
NB.
NB. void bli_?dotxv
NB.      (
NB.              conj_t  conjx,
NB.              conj_t  conjy,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  y, inc_t incy,
NB.        const ctype*  beta,
NB.              ctype*  rho
NB.      );

ddotxv_cd=: (lib,' bli_ddotxv ',ifw,' n i i ',BIC,' &d &d ',BIC,' &d ',BIC,' &d *d')&cd
zdotxv_cd=: (lib,' bli_zdotxv ',ifw,' n i i ',BIC,' &j &j ',BIC,' &j ',BIC,' &j *j')&cd

NB. Invert all elements
NB.
NB. void bli_?invertv
NB.      (
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx
NB.      );

dinvertv_cd=: (lib,' bli_dinvertv ',ifw,' n ',BIC,' *d ',BIC)&cd
zinvertv_cd=: (lib,' bli_zinvertv ',ifw,' n ',BIC,' *j ',BIC)&cd

NB. x := ( 1.0 / conjalpha(alpha) ) * x
NB.
NB. void bli_?invscalv
NB.      (
NB.              conj_t  conjalpha,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.              ctype*  x, inc_t incx
NB.      );

dinvscalv_cd=: (lib,' bli_dinvscalv ',ifw,' n i ',BIC,' &d *d ',BIC)&cd
zinvscalv_cd=: (lib,' bli_zinvscalv ',ifw,' n i ',BIC,' &j *j ',BIC)&cd

NB. x := conjalpha(alpha) * x
NB.
NB. void bli_?scalv
NB.      (
NB.              conj_t  conjalpha,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.              ctype*  x, inc_t incx
NB.      );

dscalv_cd=: (lib,' bli_dscalv ',ifw,' n i ',BIC,' &d *d ',BIC)&cd
zscalv_cd=: (lib,' bli_zscalv ',ifw,' n i ',BIC,' &j *j ',BIC)&cd

NB. y := alpha * conjx(x)
NB.
NB. void bli_?scal2v
NB.      (
NB.              conj_t  conjx,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.        const ctype*  x, inc_t incx,
NB.              ctype*  y, inc_t incy
NB.      );

dscal2v_cd=: (lib,' bli_dscal2v ',ifw,' n i ',BIC,' &d &d ',BIC,' *d ',BIC)&cd
zscal2v_cd=: (lib,' bli_zscal2v ',ifw,' n i ',BIC,' &j &j ',BIC,' *j ',BIC)&cd

NB. x := conjalpha(alpha)
NB.
NB. void bli_?setv
NB.      (
NB.              conj_t  conjalpha,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.              ctype*  x, inc_t incx
NB.      );

dsetv_cd=: (lib,' bli_dsetv ',ifw,' n i ',BIC,' &d *d ',BIC)&cd
zsetv_cd=: (lib,' bli_zsetv ',ifw,' n i ',BIC,' &d *d ',BIC)&cd

NB. y := y - conjx(x)
NB.
NB. void bli_?subv
NB.      (
NB.              conj_t  conjx,
NB.              dim_t   n,
NB.        const ctype*  x, inc_t incx,
NB.              ctype*  y, inc_t incy
NB.      );

dsubv_cd=: (lib,' bli_dsubv ',ifw,' n i ',BIC,' &d ',BIC,' *d ',BIC)&cd
zsubv_cd=: (lib,' bli_zsubv ',ifw,' n i ',BIC,' &j ',BIC,' *j ',BIC)&cd

NB. Swap corresponding elements of two n-length vectors x and y
NB.
NB. void bli_?swapv
NB.      (
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx,
NB.        ctype*  y, inc_t incy
NB.      );

dswapv_cd=: (lib,' bli_dswapv ',ifw,' n ',BIC,' *d ',BIC,' *d ',BIC)&cd
zswapv_cd=: (lib,' bli_zswapv ',ifw,' n ',BIC,' *j ',BIC,' *j ',BIC)&cd

NB. y := beta * y + conjx(x)
NB.
NB. void bli_?xpbyv
NB.      (
NB.              conj_t  conjx,
NB.              dim_t   n,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  beta,
NB.              ctype*  y, inc_t incy
NB.      );

dxpbyv_cd=: (lib,' bli_dxpbyv ',ifw,' n i ',BIC,' &d ',BIC,' &d *d ',BIC)&cd
zxpbyv_cd=: (lib,' bli_zxpbyv ',ifw,' n i ',BIC,' &j ',BIC,' &j *j ',BIC)&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1d operations

NB. void bli_?addd
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

daddd_cd=: (lib,' bli_daddd ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zaddd_cd=: (lib,' bli_zaddd ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. void bli_?axpyd
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

daxpyd_cd=: (lib,' bli_daxpyd ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zaxpyd_cd=: (lib,' bli_zaxpyd ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. void bli_?copyd
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

dcopyd_cd=: (lib,' bli_dcopyd ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zcopyd_cd=: (lib,' bli_zcopyd ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. void bli_?invertd
NB.      (
NB.        doff_t  diagoffa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  a, inc_t rsa, inc_t csa
NB.      );

dinvertd_cd=: (lib,' bli_dinvertd ',ifw,' n ',BIC,' ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zinvertd_cd=: (lib,' bli_zinvertd ',ifw,' n ',BIC,' ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. void bli_?invscald
NB.      (
NB.              conj_t  conjalpha,
NB.              doff_t  diagoffa,
NB.              dim_t   m,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dinvscald_cd=: (lib,' bli_dinvscald ',ifw,' n i ',BIC,' ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zinvscald_cd=: (lib,' bli_zinvscald ',ifw,' n i ',BIC,' ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. void bli_?scald
NB.      (
NB.              conj_t  conjalpha,
NB.              doff_t  diagoffa,
NB.              dim_t   m,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dscald_cd=: (lib,' bli_dscald ',ifw,' n i ',BIC,' ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zscald_cd=: (lib,' bli_zscald ',ifw,' n i ',BIC,' ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. void bli_?scal2d
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

dscal2d_cd=: (lib,' bli_dscal2d ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zscal2d_cd=: (lib,' bli_zscal2d ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. void bli_?setd
NB.      (
NB.              conj_t  conjalpha,
NB.              doff_t  diagoffa,
NB.              dim_t   m,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dsetd_cd=: (lib,' bli_dsetd ',ifw,' n i ',BIC,' ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zsetd_cd=: (lib,' bli_zsetd ',ifw,' n i ',BIC,' ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. void bli_?setid
NB.      (
NB.              doff_t    diagoffa,
NB.              dim_t     m,
NB.              dim_t     n,
NB.        const ctype_r*  alpha,
NB.              ctype*    a, inc_t rsa, inc_t csa
NB.      );

dsetid_cd=: (lib,' bli_dsetid ',ifw,' n ',BIC,' ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zsetid_cd=: (lib,' bli_zsetid ',ifw,' n ',BIC,' ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. void bli_?shiftd
NB.      (
NB.              doff_t  diagoffa,
NB.              dim_t   m,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dshiftd_cd=: (lib,' bli_dshiftd ',ifw,' n ',BIC,' ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zshiftd_cd=: (lib,' bli_zshiftd ',ifw,' n ',BIC,' ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. void bli_?subd
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

dsubd_cd=: (lib,' bli_dsubd ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zsubd_cd=: (lib,' bli_zsubd ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. void bli_?xpbyd
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const ctype*   beta,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

dxpbyd_cd=: (lib,' bli_dxpbyd ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zxpbyd_cd=: (lib,' bli_zxpbyd ',ifw,' n ',BIC,' i i ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1m operations

NB. B := B + transa(A)
NB.
NB. void bli_?addm
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              uplo_t   uploa,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

daddm_cd=: (lib,' bli_daddm ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zaddm_cd=: (lib,' bli_zaddm ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. B := B + alpha * transa(A)
NB.
NB. void bli_?axpym
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              uplo_t   uploa,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

daxpym_cd=: (lib,' bli_daxpym ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zaxpym_cd=: (lib,' bli_zaxpym ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. B := transa(A)
NB.
NB. void bli_?copym
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              uplo_t   uploa,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

dcopym_cd=: (lib,' bli_dcopym ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zcopym_cd=: (lib,' bli_zcopym ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. A := ( 1.0 / conjalpha(alpha) ) * A
NB.
NB. void bli_?invscalm
NB.      (
NB.              conj_t  conjalpha,
NB.              doff_t  diagoffa,
NB.              uplo_t  uploa,
NB.              dim_t   m,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dinvscalm_cd=: (lib,' bli_dinvscalm ',ifw,' n i ',BIC,' i ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zinvscalm_cd=: (lib,' bli_zinvscalm ',ifw,' n i ',BIC,' i ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. A := conjalpha(alpha) * A
NB.
NB. void bli_?scalm
NB.      (
NB.              conj_t  conjalpha,
NB.              doff_t  diagoffa,
NB.              uplo_t  uploa,
NB.              dim_t   m,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dscalm_cd=: (lib,' bli_dscalm ',ifw,' n i ',BIC,' i ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zscalm_cd=: (lib,' bli_zscalm ',ifw,' n i ',BIC,' i ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. B := alpha * transa(A)
NB.
NB. void bli_?scal2m
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              uplo_t   uploa,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

dscal2m_cd=: (lib,' bli_dscal2m ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zscal2m_cd=: (lib,' bli_zscal2m ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. Set all elements of an m x n matrix A to conjalpha(alpha)
NB.
NB. void bli_?setm
NB.      (
NB.              conj_t  conjalpha,
NB.              doff_t  diagoffa,
NB.              diag_t  diaga,
NB.              uplo_t  uploa,
NB.              dim_t   m,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dsetm_cd=: (lib,' bli_dsetm ',ifw,' n i ',BIC,' i i ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zsetm_cd=: (lib,' bli_zsetm ',ifw,' n i ',BIC,' i i ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. B := B - transa(A)
NB.
NB. void bli_?subm
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              uplo_t   uploa,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

dsubm_cd=: (lib,' bli_dsubm ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zsubm_cd=: (lib,' bli_zsubm ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-1f operations

NB. z := z + alphax * conjx(x) + alphay * conjy(y)
NB.
NB. void bli_?axpy2v
NB.      (
NB.              conj_t  conjx,
NB.              conj_t  conjy,
NB.              dim_t   m,
NB.        const ctype*  alphax,
NB.        const ctype*  alphay,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  y, inc_t incy,
NB.              ctype*  z, inc_t incz
NB.      );

daxpy2v_cd=: (lib,' bli_daxpy2v ',ifw,' n i i ',BIC,' &d &d &d ',BIC,' &d ',BIC,' *d ',BIC)&cd
zaxpy2v_cd=: (lib,' bli_zaxpy2v ',ifw,' n i i ',BIC,' &j &d &j ',BIC,' &j ',BIC,' *j ',BIC)&cd

NB. rho := conjxt(x^T) * conjy(y)
NB. z   := z + alpha * conjx(x)
NB.
NB. void bli_?dotaxpyv
NB.      (
NB.              conj_t  conjxt,
NB.              conj_t  conjx,
NB.              conj_t  conjy,
NB.              dim_t   m,
NB.        const ctype*  alpha,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  y, inc_t incy,
NB.              ctype*  rho,
NB.              ctype*  z, inc_t incz
NB.      );

ddotaxpyv_cd=: (lib,' bli_ddotaxpyv ',ifw,' n i i i ',BIC,' &d &d ',BIC,' &d ',BIC,' *d *d ',BIC)&cd
zdotaxpyv_cd=: (lib,' bli_zdotaxpyv ',ifw,' n i i i ',BIC,' &j &j ',BIC,' &j ',BIC,' *j *j ',BIC)&cd

NB. y := y + alpha * conja(A) * conjx(x)
NB.
NB. void bli_?axpyf
NB.      (
NB.              conj_t  conja,
NB.              conj_t  conjx,
NB.              dim_t   m,
NB.              dim_t   b,
NB.        const ctype*  alpha,
NB.        const ctype*  a, inc_t inca, inc_t lda,
NB.        const ctype*  x, inc_t incx,
NB.              ctype*  y, inc_t incy
NB.      );

daxpyf_cd=: (lib,' bli_daxpyf ',ifw,' n i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' *d ',BIC)&cd
zaxpyf_cd=: (lib,' bli_zaxpyf ',ifw,' n i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' *j ',BIC)&cd

NB. y := beta * y + alpha * conjat(A^T) * conjx(x)
NB.
NB. void bli_?dotxf
NB.      (
NB.              conj_t  conjat,
NB.              conj_t  conjx,
NB.              dim_t   m,
NB.              dim_t   b,
NB.        const ctype*  alpha,
NB.        const ctype*  a, inc_t inca, inc_t lda,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  beta,
NB.              ctype*  y, inc_t incy
NB.      );

ddotxf_cd=: (lib,' bli_ddotxf ',ifw,' n i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' &d *d ',BIC)&cd
zdotxf_cd=: (lib,' bli_zdotxf ',ifw,' n i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' &j *d ',BIC)&cd

NB. y := beta * y + alpha * conjat(A^T) * conjw(w)
NB. z :=        z + alpha * conja(A)    * conjx(x)
NB.
NB. void bli_?dotxaxpyf
NB.      (
NB.              conj_t  conjat,
NB.              conj_t  conja,
NB.              conj_t  conjw,
NB.              conj_t  conjx,
NB.              dim_t   m,
NB.              dim_t   b,
NB.        const ctype*  alpha,
NB.        const ctype*  a, inc_t inca, inc_t lda,
NB.        const ctype*  w, inc_t incw,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  beta,
NB.              ctype*  y, inc_t incy,
NB.              ctype*  z, inc_t incz
NB.      );

ddotxaxpyf_cd=: (lib,' bli_ddotxaxpyf ',ifw,' n i i i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' &d ',BIC,' &d *d ',BIC,' *d ',BIC)&cd
zdotxaxpyf_cd=: (lib,' bli_zdotxaxpyf ',ifw,' n i i i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' &j ',BIC,' &j *j ',BIC,' *j ',BIC)&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-2 operations

NB. y := beta * y + alpha * transa(A) * conjx(x)
NB.
NB. void bli_?gemv
NB.      (
NB.              trans_t  transa,
NB.              conj_t   conjx,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const ctype*   x, inc_t incx,
NB.        const ctype*   beta,
NB.              ctype*   y, inc_t incy
NB.      );

dgemv_cd=: (lib,' bli_dgemv ',ifw,' n i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' &d *d ',BIC)&cd
zgemv_cd=: (lib,' bli_zgemv ',ifw,' n i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' &j *j ',BIC)&cd

NB. A := A + alpha * conjx(x) * conjy(y)^T
NB.
NB. void bli_?ger
NB.      (
NB.              conj_t  conjx,
NB.              conj_t  conjy,
NB.              dim_t   m,
NB.              dim_t   n,
NB.        const ctype*  alpha,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  y, inc_t incy,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dger_cd=: (lib,' bli_dger ',ifw,' n i i ',BIC,' ',BIC,' &d &d ',BIC,' &d ',BIC,' *d ',BIC,' ',BIC)&cd
zger_cd=: (lib,' bli_zger ',ifw,' n i i ',BIC,' ',BIC,' &j &j ',BIC,' &j ',BIC,' *j ',BIC,' ',BIC)&cd

NB. y := beta * y + alpha * conja(A) * conjx(x)
NB.
NB. void bli_?hemv
NB.      (
NB.              uplo_t  uploa,
NB.              conj_t  conja,
NB.              conj_t  conjx,
NB.              dim_t   m,
NB.        const ctype*  alpha,
NB.        const ctype*  a, inc_t rsa, inc_t csa,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  beta,
NB.              ctype*  y, inc_t incy
NB.      );

dhemv_cd=: (lib,' bli_dhemv ',ifw,' n i i i ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' &d *d ',BIC)&cd
zhemv_cd=: (lib,' bli_zhemv ',ifw,' n i i i ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' &j *j ',BIC)&cd

NB. A := A + alpha * conjx(x) * conjx(x)^H
NB.
NB. void bli_?her
NB.      (
NB.              uplo_t  uploa,
NB.              conj_t  conjx,
NB.              dim_t   m,
NB.        const rtype*  alpha,
NB.        const ctype*  x, inc_t incx,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dher_cd=: (lib,' bli_dher ',ifw,' n i i ',BIC,' &d &d ',BIC,' *d ',BIC,' ',BIC)&cd
zher_cd=: (lib,' bli_zher ',ifw,' n i i ',BIC,' &d &j ',BIC,' *j ',BIC,' ',BIC)&cd

NB. A := A + alpha * conjx(x) * conjy(y)^H + conj(alpha) * conjy(y) * conjx(x)^H
NB.
NB. void bli_?her2
NB.      (
NB.              uplo_t  uploa,
NB.              conj_t  conjx,
NB.              conj_t  conjy,
NB.              dim_t   m,
NB.        const ctype*  alpha,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  y, inc_t incy,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dher2_cd=: (lib,' bli_dher2 ',ifw,' n i i i ',BIC,' &d &d ',BIC,' &d ',BIC,' *d ',BIC,' ',BIC)&cd
zher2_cd=: (lib,' bli_zher2 ',ifw,' n i i i ',BIC,' &j &j ',BIC,' &j ',BIC,' *j ',BIC,' ',BIC)&cd

NB. y := beta * y + alpha * conja(A) * conjx(x)
NB.
NB. void bli_?symv
NB.      (
NB.              uplo_t  uploa,
NB.              conj_t  conja,
NB.              conj_t  conjx,
NB.              dim_t   m,
NB.        const ctype*  alpha,
NB.        const ctype*  a, inc_t rsa, inc_t csa,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  beta,
NB.              ctype*  y, inc_t incy
NB.      );

dsymv_cd=: (lib,' bli_dsymv ',ifw,' n i i i ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' &d *d ',BIC)&cd
zsymv_cd=: (lib,' bli_zsymv ',ifw,' n i i i ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' &j *j ',BIC)&cd

NB. A := A + alpha * conjx(x) * conjx(x)^T
NB.
NB. void bli_?syr
NB.      (
NB.              uplo_t  uploa,
NB.              conj_t  conjx,
NB.              dim_t   m,
NB.        const ctype*  alpha,
NB.        const ctype*  x, inc_t incx,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dsyr_cd=: (lib,' bli_dsyr ',ifw,' n i i ',BIC,' &d &d ',BIC,' *d ',BIC,' ',BIC)&cd
zsyr_cd=: (lib,' bli_zsyr ',ifw,' n i i ',BIC,' &j &j ',BIC,' *j ',BIC,' ',BIC)&cd

NB. A := A + alpha * conjx(x) * conjy(y)^T + conj(alpha) * conjy(y) * conjx(x)^T
NB.
NB. void bli_?syr2
NB.      (
NB.              uplo_t  uploa,
NB.              conj_t  conjx,
NB.              conj_t  conjy,
NB.              dim_t   m,
NB.        const ctype*  alpha,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  y, inc_t incy,
NB.              ctype*  a, inc_t rsa, inc_t csa
NB.      );

dsyr2_cd=: (lib,' bli_dsyr2 ',ifw,' n i i i ',BIC,' &d &d ',BIC,' &d ',BIC,' *d ',BIC,' ',BIC)&cd
zsyr2_cd=: (lib,' bli_zsyr2 ',ifw,' n i i i ',BIC,' &j &j ',BIC,' &j ',BIC,' *j ',BIC,' ',BIC)&cd

NB. x := alpha * transa(A) * x
NB.
NB. void bli_?trmv
NB.      (
NB.              uplo_t   uploa,
NB.              trans_t  transa,
NB.              diag_t   diaga,
NB.              dim_t    m,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   x, inc_t incx
NB.      );

dtrmv_cd=: (lib,' bli_dtrmv ',ifw,' n i i i ',BIC,' &d &d ',BIC,' ',BIC,' *d ',BIC)&cd
ztrmv_cd=: (lib,' bli_ztrmv ',ifw,' n i i i ',BIC,' &j &j ',BIC,' ',BIC,' *j ',BIC)&cd

NB. transa(A) * x = alpha * y
NB.
NB. void bli_?trsv
NB.      (
NB.              uplo_t   uploa,
NB.              trans_t  transa,
NB.              diag_t   diaga,
NB.              dim_t    m,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   y, inc_t incy
NB.      );

dtrsv_cd=: (lib,' bli_dtrsv ',ifw,' n i i i ',BIC,' &d &d ',BIC,' ',BIC,' *d ',BIC)&cd
ztrsv_cd=: (lib,' bli_ztrsv ',ifw,' n i i i ',BIC,' &j &j ',BIC,' ',BIC,' *j ',BIC)&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Level-3 operations

NB. C := beta * C + alpha * transa(A) * transb(B)
NB.
NB. void bli_?gemm
NB.      (
NB.              trans_t  transa,
NB.              trans_t  transb,
NB.              dim_t    m,
NB.              dim_t    n,
NB.              dim_t    k,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const ctype*   b, inc_t rsb, inc_t csb,
NB.        const ctype*   beta,
NB.              ctype*   c, inc_t rsc, inc_t csc
NB.      );

dgemm_cd=: (lib,' bli_dgemm ',ifw,' n i i ',BIC,' ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zgemm_cd=: (lib,' bli_zgemm ',ifw,' n i i ',BIC,' ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. void bli_?gemm_ex
NB.      (
NB.              trans_t  transa,
NB.              trans_t  transb,
NB.              dim_t    m,
NB.              dim_t    n,
NB.              dim_t    k,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const ctype*   b, inc_t rsb, inc_t csb,
NB.        const ctype*   beta,
NB.              ctype*   c, inc_t rsc, inc_t csc,
NB.        const cntx_t*  cntx,
NB.        const rntm_t*  rntm
NB.      );

dgemm_ex_cd=: (lib,' bli_dgemm_ex ',ifw,' n i i ',BIC,' ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC,' & &')&cd
zgemm_ex_cd=: (lib,' bli_zgemm_ex ',ifw,' n i i ',BIC,' ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC,' & &')&cd

NB. C := beta * C + alpha * transa(A) * transb(B)
NB.
NB. void bli_?gemmt
NB.      (
NB.              uplo_t   uploc,
NB.              trans_t  transa,
NB.              trans_t  transb,
NB.              dim_t    m,
NB.              dim_t    k,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const ctype*   b, inc_t rsb, inc_t csb,
NB.        const ctype*   beta,
NB.              ctype*   c, inc_t rsc, inc_t csc
NB.      );

dgemmt_cd=: (lib,' bli_dgemmt ',ifw,' n i i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zgemmt_cd=: (lib,' bli_zgemmt ',ifw,' n i i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. C := beta * C + alpha * conja(A) * transb(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * transb(B) * conja(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_?hemm
NB.      (
NB.              side_t   sidea,
NB.              uplo_t   uploa,
NB.              conj_t   conja,
NB.              trans_t  transb,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const ctype*   b, inc_t rsb, inc_t csb,
NB.        const ctype*   beta,
NB.              ctype*   c, inc_t rsc, inc_t csc
NB.      );

dhemm_cd=: (lib,' bli_dhemm ',ifw,' n i i i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zhemm_cd=: (lib,' bli_zhemm ',ifw,' n i i i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. C := beta * C + alpha * transa(A) * transa(A)^H
NB.
NB. void bli_?herk
NB.      (
NB.              uplo_t   uploc,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    k,
NB.        const rtype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const rtype*   beta,
NB.              ctype*   c, inc_t rsc, inc_t csc
NB.      );

dherk_cd=: (lib,' bli_dherk ',ifw,' n i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zherk_cd=: (lib,' bli_zherk ',ifw,' n i i ',BIC,' ',BIC,' &d &j ',BIC,' ',BIC,' &d *j ',BIC,' ',BIC)&cd

NB. C := beta * C + alpha * transa(A) * transb(B)^H + conj(alpha) * transb(B) * transa(A)^H
NB.
NB. void bli_?her2k
NB.      (
NB.              uplo_t   uploc,
NB.              trans_t  transa,
NB.              trans_t  transb,
NB.              dim_t    m,
NB.              dim_t    k,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const ctype*   b, inc_t rsb, inc_t csb,
NB.        const rtype*   beta,
NB.              ctype*   c, inc_t rsc, inc_t csc
NB.      );

dher2k_cd=: (lib,' bli_dher2k ',ifw,' n i i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zher2k_cd=: (lib,' bli_zher2k ',ifw,' n i i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &d *j ',BIC,' ',BIC)&cd

NB. C := beta * C + alpha * conja(A) * transb(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * transb(B) * conja(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_?symm
NB.      (
NB.              side_t   sidea,
NB.              uplo_t   uploa,
NB.              conj_t   conja,
NB.              trans_t  transb,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const ctype*   b, inc_t rsb, inc_t csb,
NB.        const ctype*   beta,
NB.              ctype*   c, inc_t rsc, inc_t csc
NB.      );

dsymm_cd=: (lib,' bli_dsymm ',ifw,' n i i i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zsymm_cd=: (lib,' bli_zsymm ',ifw,' n i i i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. C := beta * C + alpha * transa(A) * transa(A)^T
NB.
NB. void bli_?syrk
NB.      (
NB.              uplo_t   uploc,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    k,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const ctype*   beta,
NB.              ctype*   c, inc_t rsc, inc_t csc
NB.      );

dsyrk_cd=: (lib,' bli_dsyrk ',ifw,' n i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zsyrk_cd=: (lib,' bli_zsyrk ',ifw,' n i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. C := beta * C + alpha * transa(A) * transb(B)^T + alpha * transb(B) * transa(A)^T
NB.
NB. void bli_?syr2k
NB.      (
NB.              uplo_t   uploc,
NB.              trans_t  transa,
NB.              trans_t  transb,
NB.              dim_t    m,
NB.              dim_t    k,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const ctype*   b, inc_t rsb, inc_t csb,
NB.        const ctype*   beta,
NB.              ctype*   c, inc_t rsc, inc_t csc
NB.      );

dsyr2k_cd=: (lib,' bli_dsyr2k ',ifw,' n i i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
zsyr2k_cd=: (lib,' bli_zsyr2k ',ifw,' n i i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. B := alpha * transa(A) * B  if sidea is BLIS_LEFT, or
NB. B := alpha * B * transa(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_?trmm
NB.      (
NB.              side_t   sidea,
NB.              uplo_t   uploa,
NB.              trans_t  transa,
NB.              diag_t   diaga,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

dtrmm_cd=: (lib,' bli_dtrmm ',ifw,' n i i i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
ztrmm_cd=: (lib,' bli_ztrmm ',ifw,' n i i i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. C := beta * C + alpha * transa(A) * transb(B)  if sidea is BLIS_LEFT, or
NB. C := beta * C + alpha * transb(B) * transa(A)  if sidea is BLIS_RIGHT
NB.
NB. void bli_?trmm3
NB.      (
NB.              side_t   sidea,
NB.              uplo_t   uploa,
NB.              trans_t  transa,
NB.              diag_t   diaga,
NB.              trans_t  transb,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.        const ctype*   b, inc_t rsb, inc_t csb,
NB.        const ctype*   beta,
NB.              ctype*   c, inc_t rsc, inc_t csc
NB.      );

dtrmm3_cd=: (lib,' bli_dtrmm3 ',ifw,' n i i i i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &d *d ',BIC,' ',BIC)&cd
ztrmm3_cd=: (lib,' bli_ztrmm3 ',ifw,' n i i i i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &j *j ',BIC,' ',BIC)&cd

NB. Solve the linear system with multiple right-hand sides
NB.   transa(A) * X = alpha * B  if sidea is BLIS_LEFT, or
NB.   X * transa(A) = alpha * B  if sidea is BLIS_RIGHT
NB.
NB. void bli_?trsm
NB.      (
NB.              side_t   sidea,
NB.              uplo_t   uploa,
NB.              trans_t  transa,
NB.              diag_t   diaga,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   alpha,
NB.        const ctype*   a, inc_t rsa, inc_t csa,
NB.              ctype*   b, inc_t rsb, inc_t csb
NB.      );

dtrsm_cd=: (lib,' bli_dtrsm ',ifw,' n i i i i ',BIC,' ',BIC,' &d &d ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
ztrsm_cd=: (lib,' bli_ztrsm ',ifw,' n i i i i ',BIC,' ',BIC,' &j &j ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Utility operations

NB. void bli_?asumv
NB.      (
NB.              dim_t   n,
NB.        const ctype*  x, inc_t incx,
NB.              rtype*  asum
NB.      );

dasumv_cd=: (lib,' bli_dasumv ',ifw,' n ',BIC,' &d ',BIC,' *d')&cd
zasumv_cd=: (lib,' bli_zasumv ',ifw,' n ',BIC,' &j ',BIC,' *d')&cd

NB. void bli_?norm[1fi]m
NB.      (
NB.              doff_t  diagoffa,
NB.              doff_t  diaga,
NB.              uplo_t  uploa,
NB.              dim_t   m,
NB.              dim_t   n,
NB.        const ctype*  a, inc_t rs_a, inc_t cs_a,
NB.              rtype*  norm
NB.      );

dnorm1m_cd=: (lib,' bli_dnorm1m ',ifw,' n ',BIC,' ',BIC,' i ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' *d')&cd
dnormfm_cd=: (lib,' bli_dnormfm ',ifw,' n ',BIC,' ',BIC,' i ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' *d')&cd
dnormim_cd=: (lib,' bli_dnormim ',ifw,' n ',BIC,' ',BIC,' i ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' *d')&cd

znorm1m_cd=: (lib,' bli_znorm1m ',ifw,' n ',BIC,' ',BIC,' i ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' *d')&cd
znormfm_cd=: (lib,' bli_znormfm ',ifw,' n ',BIC,' ',BIC,' i ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' *d')&cd
znormim_cd=: (lib,' bli_znormim ',ifw,' n ',BIC,' ',BIC,' i ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' *d')&cd

NB. void bli_?norm[1fi]v
NB.      (
NB.              dim_t   n,
NB.        const ctype*  x, inc_t incx,
NB.              rtype*  norm
NB.      );

dnorm1v_cd=: (lib,' bli_dnorm1v ',ifw,' n ',BIC,' &d ',BIC,' *d')&cd
dnormfv_cd=: (lib,' bli_dnormfv ',ifw,' n ',BIC,' &d ',BIC,' *d')&cd
dnormiv_cd=: (lib,' bli_dnormiv ',ifw,' n ',BIC,' &d ',BIC,' *d')&cd

znorm1v_cd=: (lib,' bli_znorm1v ',ifw,' n ',BIC,' &j ',BIC,' *d')&cd
znormfv_cd=: (lib,' bli_znormfv ',ifw,' n ',BIC,' &j ',BIC,' *d')&cd
znormiv_cd=: (lib,' bli_znormiv ',ifw,' n ',BIC,' &j ',BIC,' *d')&cd

NB. void bli_?mkherm
NB.      (
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        ctype*  a, inc_t rs_a, inc_t cs_a
NB.      );

dmkherm_cd=: (lib,' bli_dmkherm ',ifw,' n i ',BIC,' *d ',BIC,' ',BIC)&cd
zmkherm_cd=: (lib,' bli_zmkherm ',ifw,' n i ',BIC,' *j ',BIC,' ',BIC)&cd

NB. void bli_?mksymm
NB.      (
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        ctype*  a, inc_t rs_a, inc_t cs_a
NB.      );

dmksymm_cd=: (lib,' bli_dmksymm ',ifw,' n i ',BIC,' *d ',BIC,' ',BIC)&cd
zmksymm_cd=: (lib,' bli_zmksymm ',ifw,' n i ',BIC,' *j ',BIC,' ',BIC)&cd

NB. void bli_?mktrim
NB.      (
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        ctype*  a, inc_t rs_a, inc_t cs_a
NB.      );

dmktrim_cd=: (lib,' bli_dmktrim ',ifw,' n i ',BIC,' *d ',BIC,' ',BIC)&cd
zmktrim_cd=: (lib,' bli_zmktrim ',ifw,' n i ',BIC,' *j ',BIC,' ',BIC)&cd

NB. void bli_?fprintv
NB.      (
NB.              FILE*   file,
NB.        const char*   s1,
NB.              dim_t   m,
NB.        const ctype*  x, inc_t incx,
NB.        const char*   format,
NB.        const char*   s2
NB.      );

dfprintv_cd=: (lib,' bli_dfprintv ',ifw,' n & &c ',BIC,' &d ',BIC,' &c &c')&cd
zfprintv_cd=: (lib,' bli_zfprintv ',ifw,' n & &c ',BIC,' &j ',BIC,' &c &c')&cd

NB. void bli_?fprintm
NB.      (
NB.              FILE*   file,
NB.        const char*   s1,
NB.              dim_t   m,
NB.              dim_t   n,
NB.        const ctype*  a, inc_t rs_a, inc_t cs_a,
NB.        const char*   format,
NB.        const char*   s2
NB.      );

dfprintm_cd=: (lib,' bli_dfprintm ',ifw,' n & &c ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &c &c')&cd
zfprintm_cd=: (lib,' bli_zfprintm ',ifw,' n & &c ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &c &c')&cd

NB. void bli_?printv
NB.      (
NB.        const char*   s1,
NB.              dim_t   m,
NB.        const ctype*  x, inc_t incx,
NB.        const char*   format,
NB.        const char*   s2
NB.      );

dprintv_cd=: (lib,' bli_dprintv ',ifw,' n &c ',BIC,' &d ',BIC,' &c &c')&cd
zprintv_cd=: (lib,' bli_zprintv ',ifw,' n &c ',BIC,' &j ',BIC,' &c &c')&cd

NB. void bli_?printm
NB.      (
NB.        const char*   s1,
NB.              dim_t   m,
NB.              dim_t   n,
NB.        const ctype*  a, inc_t rs_a, inc_t cs_a,
NB.        const char*   format,
NB.        const char*   s2
NB.      );

dprintm_cd=: (lib,' bli_dprintm ',ifw,' n &c ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &c &c')&cd
zprintm_cd=: (lib,' bli_zprintm ',ifw,' n &c ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &c &c')&cd

NB. void bli_?randv
NB.      (
NB.        dim_t   n,
NB.        ctype*  x, inc_t incx
NB.      );

drandv_cd=: (lib,' bli_drandv ',ifw,' n ',BIC,' *d ',BIC)&cd
zrandv_cd=: (lib,' bli_zrandv ',ifw,' n ',BIC,' *j ',BIC)&cd

NB. void bli_?randm
NB.      (
NB.        doff_t  diagoffa,
NB.        uplo_t  uploa,
NB.        dim_t   m,
NB.        dim_t   n,
NB.        ctype*  a, inc_t rs_a, inc_t cs_a
NB.      );

drandm_cd=: (lib,' bli_drandm ',ifw,' n ',BIC,' i ',BIC,' ',BIC,' *d ',BIC,' ',BIC)&cd
zrandm_cd=: (lib,' bli_zrandm ',ifw,' n ',BIC,' i ',BIC,' ',BIC,' *j ',BIC,' ',BIC)&cd

NB. void bli_?sumsqv
NB.      (
NB.              dim_t   n,
NB.        const ctype*  x, inc_t incx,
NB.              rtype*  scale,
NB.              rtype*  sumsq
NB.      );

dsumsqv_cd=: (lib,' bli_dsumsqv ',ifw,' n ',BIC,' &d ',BIC,' *d *d')&cd
zsumsqv_cd=: (lib,' bli_zsumsqv ',ifw,' n ',BIC,' &j ',BIC,' *d *d')&cd

NB. void bli_?getsc
NB.      (
NB.        const ctype*   chi,
NB.              double*  zeta_r,
NB.              double*  zeta_i
NB.      );

dgetsc_cd=: (lib,' bli_dgetsc ',ifw,' n &d *d *d')&cd
zgetsc_cd=: (lib,' bli_zgetsc ',ifw,' n &j *d *d')&cd

NB. err_t bli_?getijv
NB.      (
NB.              dim_t    i,
NB.        const void*    x, inc_t incx,
NB.              double*  ar,
NB.              double*  ai
NB.      );

dgetijv_cd=: (lib,' bli_dgetijv ',ifw,' i ',BIC,' & ',BIC,' *d *d')&cd
zgetijv_cd=: (lib,' bli_zgetijv ',ifw,' i ',BIC,' & ',BIC,' *d *d')&cd

NB. err_t bli_?getijm
NB.      (
NB.              dim_t    i,
NB.              dim_t    j,
NB.        const void*    b, inc_t rs_b, inc_t cs_b,
NB.              double*  ar,
NB.              double*  ai
NB.      );

dgetijm_cd=: (lib,' bli_dgetijm ',ifw,' i ',BIC,' ',BIC,' & ',BIC,' ',BIC,' *d *d')&cd
zgetijm_cd=: (lib,' bli_zgetijm ',ifw,' i ',BIC,' ',BIC,' & ',BIC,' ',BIC,' *d *d')&cd

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
NB.        void*   x, inc_t incx
NB.      );

dsetijv_cd=: (lib,' bli_dsetijv ',ifw,' i d d ',BIC,' * ',BIC)&cd
zsetijv_cd=: (lib,' bli_zsetijv ',ifw,' i d d ',BIC,' * ',BIC)&cd

NB. err_t bli_?setijm
NB.      (
NB.        double  ar,
NB.        double  ai,
NB.        dim_t   i,
NB.        dim_t   j,
NB.        void*   b, inc_t rs_b, inc_t cs_b
NB.      );

dsetijm_cd=: (lib,' bli_dsetijm ',ifw,' i d d ',BIC,' ',BIC,' * ',BIC,' ',BIC)&cd
zsetijm_cd=: (lib,' bli_zsetijm ',ifw,' i d d ',BIC,' ',BIC,' * ',BIC,' ',BIC)&cd

NB. void bli_?eqsc
NB.      (
NB.              conj_t  conjchi,
NB.        const ctype*  chi,
NB.        const ctype*  psi,
NB.              bool*   is_eq
NB.      );

deqsc_cd=: (lib,' bli_deqsc ',ifw,' n i &d &d *b')&cd
zeqsc_cd=: (lib,' bli_zeqsc ',ifw,' n i &j &j *b')&cd

NB. void bli_?eqv
NB.      (
NB.              conj_t  conjx,
NB.              dim_t   n,
NB.        const ctype*  x, inc_t incx,
NB.        const ctype*  y, inc_t incy,
NB.              bool*   is_eq
NB.      );

deqv_cd=: (lib,' bli_deqv ',ifw,' n i ',BIC,' &d ',BIC,' &d ',BIC,' *b')&cd
zeqv_cd=: (lib,' bli_zeqv ',ifw,' n i ',BIC,' &j ',BIC,' &j ',BIC,' *b')&cd

NB. void bli_?eqm
NB.      (
NB.              doff_t   diagoffa,
NB.              diag_t   diaga,
NB.              uplo_t   uploa,
NB.              trans_t  transa,
NB.              dim_t    m,
NB.              dim_t    n,
NB.        const ctype*   a, inc_t rs_a, inc_t cs_a,
NB.        const ctype*   b, inc_t rs_b, inc_t cs_b,
NB.              bool*    is_eq
NB.      );

deqm_cd=: (lib,' bli_deqm ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' &d ',BIC,' ',BIC,' *b')&cd
zeqm_cd=: (lib,' bli_zeqm ',ifw,' n ',BIC,' i i i ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' &j ',BIC,' ',BIC,' *b')&cd

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. Query function reference

NB. ---------------------------------------------------------
NB. General library information

NB. const char* bli_info_get_version_str( void );
NB. note: see ver in external/blis/util.ijs
info_get_version_str_cd=: (lib,' bli_info_get_version_str > ',ifw,' x')&cd

NB. ---------------------------------------------------------
NB. Specific configuration

NB. arch_t bli_arch_query_id( void );
arch_query_id_cd=: (lib,' bli_arch_query_id > ',ifw,' i')&cd

NB. const char* bli_arch_string( arch_t id );
NB. note: see arch in external/blis/util.ijs
arch_string_cd=: (lib,' bli_arch_string > ',ifw,' x i')&cd

NB. ---------------------------------------------------------
NB. General configuration

NB. gint_t bli_info_get_int_type_size( void );
info_get_int_type_size_cd=:              (lib,' bli_info_get_int_type_size > '             ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_num_fp_types( void );
info_get_num_fp_types_cd=:               (lib,' bli_info_get_num_fp_types > '              ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_max_type_size( void );
info_get_max_type_size_cd=:              (lib,' bli_info_get_max_type_size > '             ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_page_size( void );
info_get_page_size_cd=:                  (lib,' bli_info_get_page_size > '                 ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_simd_num_registers( void );
info_get_simd_num_registers_cd=:         (lib,' bli_info_get_simd_num_registers > '        ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_simd_size( void );
info_get_simd_size_cd=:                  (lib,' bli_info_get_simd_size > '                 ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_simd_align_size( void );
info_get_simd_align_size_cd=:            (lib,' bli_info_get_simd_align_size > '           ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_stack_buf_max_size( void );
info_get_stack_buf_max_size_cd=:         (lib,' bli_info_get_stack_buf_max_size > '        ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_stack_buf_align_size( void );
info_get_stack_buf_align_size_cd=:       (lib,' bli_info_get_stack_buf_align_size > '      ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_heap_addr_align_size( void );
info_get_heap_addr_align_size_cd=:       (lib,' bli_info_get_heap_addr_align_size > '      ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_heap_stride_align_size( void );
info_get_heap_stride_align_size_cd=:     (lib,' bli_info_get_heap_stride_align_size > '    ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_pool_addr_align_size_a( void );
info_get_pool_addr_align_size_a_cd=:     (lib,' bli_info_get_pool_addr_align_size_a > '    ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_pool_addr_align_size_b( void );
info_get_pool_addr_align_size_b_cd=:     (lib,' bli_info_get_pool_addr_align_size_b > '    ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_pool_addr_align_size_c( void );
info_get_pool_addr_align_size_c_cd=:     (lib,' bli_info_get_pool_addr_align_size_c > '    ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_pool_addr_align_size_gen( void );
info_get_pool_addr_align_size_gen_cd=:   (lib,' bli_info_get_pool_addr_align_size_gen > '  ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_pool_addr_offset_size_a( void );
info_get_pool_addr_offset_size_a_cd=:    (lib,' bli_info_get_pool_addr_offset_size_a > '   ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_pool_addr_offset_size_b( void );
info_get_pool_addr_offset_size_b_cd=:    (lib,' bli_info_get_pool_addr_offset_size_b > '   ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_pool_addr_offset_size_c( void );
info_get_pool_addr_offset_size_c_cd=:    (lib,' bli_info_get_pool_addr_offset_size_c > '   ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_pool_addr_offset_size_gen( void );
info_get_pool_addr_offset_size_gen_cd=:  (lib,' bli_info_get_pool_addr_offset_size_gen > ' ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_blas( void );
info_get_enable_blas_cd=:                (lib,' bli_info_get_enable_blas > '               ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_cblas( void );
info_get_enable_cblas_cd=:               (lib,' bli_info_get_enable_cblas > '              ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_blas_int_type_size( void );
info_get_blas_int_type_size_cd=:         (lib,' bli_info_get_blas_int_type_size > '        ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_pba_pools( void );
info_get_enable_pba_pools_cd=:           (lib,' bli_info_get_enable_pba_pools > '          ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_sba_pools( void );
info_get_enable_sba_pools_cd=:           (lib,' bli_info_get_enable_sba_pools > '          ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_threading( void );
info_get_enable_threading_cd=:           (lib,' bli_info_get_enable_threading > '          ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_openmp( void );
info_get_enable_openmp_cd=:              (lib,' bli_info_get_enable_openmp > '             ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_pthreads( void );
info_get_enable_pthreads_cd=:            (lib,' bli_info_get_enable_pthreads > '           ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_hpx( void );
info_get_enable_hpx_cd=:                 (lib,' bli_info_get_enable_hpx > '                ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_openmp_as_default( void );
info_get_enable_openmp_as_default_cd=:   (lib,' bli_info_get_enable_openmp_as_default > '  ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_pthreads_as_default( void );
info_get_enable_pthreads_as_default_cd=: (lib,' bli_info_get_enable_pthreads_as_default > ',ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_hpx_as_default( void );
info_get_enable_hpx_as_default_cd=:      (lib,' bli_info_get_enable_hpx_as_default > '     ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_thread_jrir_slab( void );
info_get_thread_jrir_slab_cd=:           (lib,' bli_info_get_thread_jrir_slab > '          ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_thread_jrir_rr( void );
info_get_thread_jrir_rr_cd=:             (lib,' bli_info_get_thread_jrir_rr > '            ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_thread_jrir_tlb( void );
info_get_thread_jrir_tlb_cd=:            (lib,' bli_info_get_thread_jrir_tlb > '           ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_tls( void );
info_get_enable_tls_cd=:                 (lib,' bli_info_get_enable_tls > '                ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_memkind( void );
info_get_enable_memkind_cd=:             (lib,' bli_info_get_enable_memkind > '            ,ifw,' ',BIC)&cd

NB. gint_t bli_info_get_enable_sandbox( void );
info_get_enable_sandbox_cd=:             (lib,' bli_info_get_enable_sandbox > '            ,ifw,' ',BIC)&cd

NB. ---------------------------------------------------------
NB. Kernel information

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Micro-kernel implementation type query

NB. const char* bli_info_get_gemm_ukr_impl_string( ind_t method, num_t dt );
NB. note: see info_get_gemm_ukr_impl_string in external/blis/util.ijs
info_get_gemm_ukr_impl_string_cd=:       (lib,' bli_info_get_gemm_ukr_impl_string > '      ,ifw,' x i i')&cd

NB. const char* bli_info_get_gemmtrsm_l_ukr_impl_string( ind_t method, num_t dt );
NB. note: see info_get_gemmtrsm_l_ukr_impl_string in external/blis/util.ijs
info_get_gemmtrsm_l_ukr_impl_string_cd=: (lib,' bli_info_get_gemmtrsm_l_ukr_impl_string > ',ifw,' x i i')&cd

NB. const char* bli_info_get_gemmtrsm_u_ukr_impl_string( ind_t method, num_t dt );
NB. note: see info_get_gemmtrsm_u_ukr_impl_string in external/blis/util.ijs
info_get_gemmtrsm_u_ukr_impl_string_cd=: (lib,' bli_info_get_gemmtrsm_u_ukr_impl_string > ',ifw,' x i i')&cd

NB. const char* bli_info_get_trsm_l_ukr_impl_string( ind_t method, num_t dt );
NB. note: see info_get_trsm_l_ukr_impl_string in external/blis/util.ijs
info_get_trsm_l_ukr_impl_string_cd=:     (lib,' bli_info_get_trsm_l_ukr_impl_string > '    ,ifw,' x i i')&cd

NB. const char* bli_info_get_trsm_u_ukr_impl_string( ind_t method, num_t dt );
NB. note: see info_get_trsm_u_ukr_impl_string in external/blis/util.ijs
info_get_trsm_u_ukr_impl_string_cd=:     (lib,' bli_info_get_trsm_u_ukr_impl_string > '    ,ifw,' x i i')&cd

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. Operation implementation type query

NB. const char* bli_info_get_gemm_impl_string( num_t dt );
NB. note: see info_get_gemm_impl_string in external/blis/util.ijs
info_get_gemm_impl_string_cd=:  (lib,' bli_info_get_gemm_impl_string > ' ,ifw,' x i')&cd

NB. const char* bli_info_get_gemmt_impl_string( num_t dt );
NB. note: see info_get_gemmt_impl_string in external/blis/util.ijs
info_get_gemmt_impl_string_cd=: (lib,' bli_info_get_gemmt_impl_string > ',ifw,' x i')&cd

NB. const char* bli_info_get_hemm_impl_string( num_t dt );
NB. note: see info_get_hemm_impl_string in external/blis/util.ijs
info_get_hemm_impl_string_cd=:  (lib,' bli_info_get_hemm_impl_string > ' ,ifw,' x i')&cd

NB. const char* bli_info_get_herk_impl_string( num_t dt );
NB. note: see info_get_herk_impl_string in external/blis/util.ijs
info_get_herk_impl_string_cd=:  (lib,' bli_info_get_herk_impl_string > ' ,ifw,' x i')&cd

NB. const char* bli_info_get_her2k_impl_string( num_t dt );
NB. note: see info_get_her2k_impl_string in external/blis/util.ijs
info_get_her2k_impl_string_cd=: (lib,' bli_info_get_her2k_impl_string > ',ifw,' x i')&cd

NB. const char* bli_info_get_symm_impl_string( num_t dt );
NB. note: see info_get_symm_impl_string in external/blis/util.ijs
info_get_symm_impl_string_cd=:  (lib,' bli_info_get_symm_impl_string > ' ,ifw,' x i')&cd

NB. const char* bli_info_get_syrk_impl_string( num_t dt );
NB. note: see info_get_syrk_impl_string in external/blis/util.ijs
info_get_syrk_impl_string_cd=:  (lib,' bli_info_get_syrk_impl_string > ' ,ifw,' x i')&cd

NB. const char* bli_info_get_syr2k_impl_string( num_t dt );
NB. note: see info_get_syr2k_impl_string in external/blis/util.ijs
info_get_syr2k_impl_string_cd=: (lib,' bli_info_get_syr2k_impl_string > ',ifw,' x i')&cd

NB. const char* bli_info_get_trmm_impl_string( num_t dt );
NB. note: see info_get_trmm_impl_string in external/blis/util.ijs
info_get_trmm_impl_string_cd=:  (lib,' bli_info_get_trmm_impl_string > ' ,ifw,' x i')&cd

NB. const char* bli_info_get_trmm3_impl_string( num_t dt );
NB. note: see info_get_trmm3_impl_string in external/blis/util.ijs
info_get_trmm3_impl_string_cd=: (lib,' bli_info_get_trmm3_impl_string > ',ifw,' x i')&cd

NB. const char* bli_info_get_trsm_impl_string( num_t dt );
NB. note: see info_get_trsm_impl_string in external/blis/util.ijs
info_get_trsm_impl_string_cd=:  (lib,' bli_info_get_trsm_impl_string > ' ,ifw,' x i')&cd

NB. =========================================================
NB. Clean-up

erase 'lib ifw trash'
