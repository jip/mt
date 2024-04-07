NB. API definitions
NB.
NB. xxxxx_cd  Cover verbs to call BLAS subroutine or function
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
NB. Concepts
NB.
NB. Conventions:
NB. 1) LIB_mtbla_ global noun must exist

NB. =========================================================
NB. Configuration

coclass 'mtbla'

NB. =========================================================
NB. Local definitions

lib=. dquote LIB
ifw=. IFWIN # '+'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Level 1 BLAS

NB. SUBROUTINE _ROTG ( A, B, C, S )
drotg_cd=: (lib,' drotg_ ',ifw,' n &d &d *d *d')&cd
zrotg_cd=: (lib,' zrotg_ ',ifw,' n &j &j *d *j')&cd

NB. SUBROUTINE DROTMG( D1, D2, A, B, PARAM )
drotmg_cd=: (lib,' drotmg_ ',ifw,' n *d *d *d &d *d')&cd

NB. SUBROUTINE DROT ( N, X, INCX, Y, INCY, C, S )
drot_cd=: (lib,' drot_ ',ifw,' n &i *d &i *d &i &d &d')&cd

NB. SUBROUTINE ZDROT( N, CX, INCX, CY, INCY, C, S )
zdrot_cd=: (lib,' zdrot_ ',ifw,' n &i *j &i *j &i &d &d')&cd

NB. SUBROUTINE DROTM ( N, X, INCX, Y, INCY, PARAM )
drotm_cd=: (lib,' drotm_ ',ifw,' n &i *d &i *d &i &d')&cd

NB. SUBROUTINE _SWAP ( N, X, INCX, Y, INCY )
dswap_cd=: (lib,' dswap_ ',ifw,' n &i *d &i *d &i')&cd
zswap_cd=: (lib,' zswap_ ',ifw,' n &i *j &i *j &i')&cd

NB. SUBROUTINE _SCAL ( N, ALPHA, X, INCX )
dscal_cd=: (lib,' dscal_ ',ifw,' n &i &d *d &i')&cd
zscal_cd=: (lib,' zscal_ ',ifw,' n &i &j *j &i')&cd
zdscal_cd=: (lib,' zdscal_ ',ifw,' n &i &d *j &i')&cd

NB. SUBROUTINE _COPY ( N, X, INCX, Y, INCY )
dcopy_cd=: (lib,' dcopy_ ',ifw,' n &i &d &i *d &i')&cd
zcopy_cd=: (lib,' zcopy_ ',ifw,' n &i &j &i *j &i')&cd

NB. SUBROUTINE _AXPY ( N, ALPHA, X, INCX, Y, INCY )
daxpy_cd=: (lib,' daxpy_ ',ifw,' n &i &d &d &i *d &i')&cd
zaxpy_cd=: (lib,' zaxpy_ ',ifw,' n &i &j &j &i *j &i')&cd

NB. FUNCTION DDOT ( N, X, INCX, Y, INCY )
ddot_cd=: (lib,' ddot_ ',ifw,' d &i &d &i &d &i')&cd

NB. FUNCTION ZDOTU ( N, X, INCX, Y, INCY )
zdotu_cd=: (lib,' zdotu_ ',ifw,' j &i &f &i &f &i')&cd

NB. FUNCTION ZDOTC ( N, X, INCX, Y, INCY )
zdotc_cd=: (lib,' zdotc_ ',ifw,' j &i &f &i &f &i')&cd

NB. FUNCTION _NRM2 ( N, X, INCX )
dnrm2_cd=: (lib,' dnrm2_ ',ifw,' d &i &d &i')&cd
dznrm2_cd=: (lib,' dznrm2_ ',ifw,' d &i &j &i')&cd

NB. FUNCTION _ASUM ( N, X, INCX )
dasum_cd=: (lib,' dasum_ ',ifw,' d &i &d &i')&cd
dzasum_cd=: (lib,' dzasum_ ',ifw,' d &i *j &i')&cd

NB. FUNCTION I_AMAX( N, X, INCX )
idamax_cd=: (lib,' idamax_ ',ifw,' i &i &d &i')&cd
izamax_cd=: (lib,' izamax_ ',ifw,' i &i &j &i')&cd

NB. FUNCTION DCABS1(Z)
dcabs1_cd=: (lib,' dcabs1_ ',ifw,' d &j')&cd

NB. ---------------------------------------------------------
NB. Level 2 BLAS

NB. SUBROUTINE _GEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
dgemv_cd=: (lib,' dgemv_ ',ifw,' n &c &i &i &d &d &i &d &i &d *d &i')&cd
zgemv_cd=: (lib,' zgemv_ ',ifw,' n &c &i &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
dsymv_cd=: (lib,' dsymv_ ',ifw,' n &c &i &d &d &i &d &i &d *d &i')&cd

NB. SUBROUTINE ZHEMV ( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
zhemv_cd=: (lib,' zhemv_ ',ifw,' n &c &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE _TRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
dtrmv_cd=: (lib,' dtrmv_ ',ifw,' n &c &c &c &i &d &i *d &i')&cd
ztrmv_cd=: (lib,' ztrmv_ ',ifw,' n &c &c &c &i &j &i *j &i')&cd

NB. SUBROUTINE _TRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
dtrsv_cd=: (lib,' dtrsv_ ',ifw,' n &c &c &c &i &d &i *d &i')&cd
ztrsv_cd=: (lib,' ztrsv_ ',ifw,' n &c &c &c &i &j &i *j &i')&cd

NB. SUBROUTINE DGER ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
dger_cd=: (lib,' dger_ ',ifw,' n &i &i &d &d &i &d &i *d &i')&cd

NB. SUBROUTINE ZGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
zgerc_cd=: (lib,' zgerc_ ',ifw,' n &i &i &j &j &i &j &i *j &i')&cd

NB. SUBROUTINE ZGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
zgeru_cd=: (lib,' zgeru_ ',ifw,' n &i &i &j &j &i &j &i *j &i')&cd

NB. SUBROUTINE DSYR ( UPLO, N, ALPHA, X, INCX, A, LDA )
dsyr_cd=: (lib,' dsyr_ ',ifw,' n &c &i &d &d &i *d &i')&cd

NB. SUBROUTINE DSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
dsyr2_cd=: (lib,' dsyr2_ ',ifw,' n &c &i &d &d &i &d &i *d &i')&cd

NB. SUBROUTINE ZHER ( UPLO, N, ALPHA, X, INCX, A, LDA )
zher_cd=: (lib,' zher_ ',ifw,' n &c &i &d &j &i *j &i')&cd

NB. SUBROUTINE ZHER2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
zher2_cd=: (lib,' zher2_ ',ifw,' n &c &i &j &j &i &j &i *j &i')&cd

NB. ---------------------------------------------------------
NB. Level 3 BLAS

NB. SUBROUTINE _GEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
dgemm_cd=: (lib,' dgemm_ ',ifw,' n &c &c &i &i &i &d &d &i &d &i &d *d &i')&cd
zgemm_cd=: (lib,' zgemm_ ',ifw,' n &c &c &i &i &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE _SYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
dsymm_cd=: (lib,' dsymm_ ',ifw,' n &c &c &i &i &d &d &i &d &i &d *d &i')&cd
zsymm_cd=: (lib,' zsymm_ ',ifw,' n &c &c &i &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE ZHEMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
zhemm_cd=: (lib,' zhemm_ ',ifw,' n &c &c &i &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE _TRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )
dtrmm_cd=: (lib,' dtrmm_ ',ifw,' n &c &c &c &c &i &i &d &d &i *d &i')&cd
ztrmm_cd=: (lib,' ztrmm_ ',ifw,' n &c &c &c &c &i &i &j &j &i *j &i')&cd

NB. SUBROUTINE _SYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
dsyrk_cd=: (lib,' dsyrk_ ',ifw,' n &c &c &i &i &d &d &i &d *d &i')&cd
zsyrk_cd=: (lib,' zsyrk_ ',ifw,' n &c &c &i &i &j &j &i &j *j &i')&cd

NB. SUBROUTINE ZHERK ( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
zherk_cd=: (lib,' zherk_ ',ifw,' n &c &c &i &i &d &j &i &d *j &i')&cd

NB. SUBROUTINE _SYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
dsyr2k_cd=: (lib,' dsyr2k_ ',ifw,' n &c &c &i &i &d &d &i &d &i &d *d &i')&cd
zsyr2k_cd=: (lib,' zsyr2k_ ',ifw,' n &c &c &i &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE ZHER2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
zher2k_cd=: (lib,' zher2k_ ',ifw,' n &c &c &i &i &j &j &i &j &i &d *j &i')&cd

NB. SUBROUTINE _TRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )
dtrsm_cd=: (lib,' dtrsm_ ',ifw,' n &c &c &c &c &i &i &d &d &i *d &i')&cd
ztrsm_cd=: (lib,' ztrsm_ ',ifw,' n &c &c &c &c &i &i &j &j &i *j &i')&cd

NB. =========================================================
NB. Clean-up

erase 'lib ifw'
