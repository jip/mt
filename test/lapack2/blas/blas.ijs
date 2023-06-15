NB. Interface to BLAS
NB.
NB. xxxxxcd   Cover verbs to call BLAS subroutine or function
NB. basicxxx  Utilities to either check or modify argument
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

coclass 'mttst'
coinsert 'mt'

NB. =========================================================
NB. Includes

require 'math/lapack2'

NB. =========================================================
NB. Local definitions

lib=. dquote liblapack_jlapack2_
ifw=. IFWIN # '+'

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. Level 1 BLAS

NB. SUBROUTINE _ROTG ( A, B, C, S )
drotgcd=: (lib,' drotg_ ',ifw,' n &d &d *d *d')&cd
zrotgcd=: (lib,' zrotg_ ',ifw,' n &j &j *d *j')&cd

NB. SUBROUTINE DROTMG( D1, D2, A, B, PARAM )
drotmgcd=: (lib,' drotmg_ ',ifw,' n *d *d *d &d *d')&cd

NB. SUBROUTINE DROT ( N, X, INCX, Y, INCY, C, S )
drotcd=: (lib,' drot_ ',ifw,' n &i *d &i *d &i &d &d')&cd

NB. SUBROUTINE ZDROT( N, CX, INCX, CY, INCY, C, S )
zdrotcd=: (lib,' zdrot_ ',ifw,' n &i *j &i *j &i &d &d')&cd

NB. SUBROUTINE DROTM ( N, X, INCX, Y, INCY, PARAM )
drotmcd=: (lib,' drotm_ ',ifw,' n &i *d &i *d &i &d')&cd

NB. SUBROUTINE _SWAP ( N, X, INCX, Y, INCY )
dswapcd=: (lib,' dswap_ ',ifw,' n &i *d &i *d &i')&cd
zswapcd=: (lib,' zswap_ ',ifw,' n &i *j &i *j &i')&cd

NB. SUBROUTINE _SCAL ( N, ALPHA, X, INCX )
dscalcd=: (lib,' dscal_ ',ifw,' n &i &d *d &i')&cd
zscalcd=: (lib,' zscal_ ',ifw,' n &i &j *j &i')&cd
zdscalcd=: (lib,' zdscal_ ',ifw,' n &i &d *j &i')&cd

NB. SUBROUTINE _COPY ( N, X, INCX, Y, INCY )
dcopycd=: (lib,' dcopy_ ',ifw,' n &i &d &i *d &i')&cd
zcopycd=: (lib,' zcopy_ ',ifw,' n &i &j &i *j &i')&cd

NB. SUBROUTINE _AXPY ( N, ALPHA, X, INCX, Y, INCY )
daxpycd=: (lib,' daxpy_ ',ifw,' n &i &d &d &i *d &i')&cd
zaxpycd=: (lib,' zaxpy_ ',ifw,' n &i &j &j &i *j &i')&cd

NB. FUNCTION DDOT ( N, X, INCX, Y, INCY )
ddotcd=: (lib,' ddot_ ',ifw,' d &i &d &i &d &i')&cd

NB. FUNCTION ZDOTU ( N, X, INCX, Y, INCY )
zdotucd=: (lib,' zdotu_ ',ifw,' j &i &f &i &f &i')&cd

NB. FUNCTION ZDOTC ( N, X, INCX, Y, INCY )
zdotccd=: (lib,' zdotc_ ',ifw,' j &i &f &i &f &i')&cd

NB. FUNCTION _NRM2 ( N, X, INCX )
dnrm2cd=: (lib,' dnrm2_ ',ifw,' d &i &d &i')&cd
dznrm2cd=: (lib,' dznrm2_ ',ifw,' d &i &j &i')&cd

NB. FUNCTION _ASUM ( N, X, INCX )
dasumcd=: (lib,' dasum_ ',ifw,' d &i &d &i')&cd
dzasumcd=: (lib,' dzasum_ ',ifw,' d &i *j &i')&cd

NB. FUNCTION I_AMAX( N, X, INCX )
idamaxcd=: (lib,' idamax_ ',ifw,' i &i &d &i')&cd
izamaxcd=: (lib,' izamax_ ',ifw,' i &i &j &i')&cd

NB. FUNCTION DCABS1(Z)
dcabs1cd=: (lib,' dcabs1_ ',ifw,' d &j')&cd

NB. ---------------------------------------------------------
NB. Level 2 BLAS

NB. SUBROUTINE _GEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
dgemvcd=: (lib,' dgemv_ ',ifw,' n &c &i &i &d &d &i &d &i &d *d &i')&cd
zgemvcd=: (lib,' zgemv_ ',ifw,' n &c &i &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
dsymvcd=: (lib,' dsymv_ ',ifw,' n &c &i &d &d &i &d &i &d *d &i')&cd

NB. SUBROUTINE ZHEMV ( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
zhemvcd=: (lib,' zhemv_ ',ifw,' n &c &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE _TRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
dtrmvcd=: (lib,' dtrmv_ ',ifw,' n &c &c &c &i &d &i *d &i')&cd
ztrmvcd=: (lib,' ztrmv_ ',ifw,' n &c &c &c &i &j &i *j &i')&cd

NB. SUBROUTINE _TRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
dtrsvcd=: (lib,' dtrsv_ ',ifw,' n &c &c &c &i &d &i *d &i')&cd
ztrsvcd=: (lib,' ztrsv_ ',ifw,' n &c &c &c &i &j &i *j &i')&cd

NB. SUBROUTINE DGER ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
dgercd=: (lib,' dger_ ',ifw,' n &i &i &d &d &i &d &i *d &i')&cd

NB. SUBROUTINE ZGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
zgerccd=: (lib,' zgerc_ ',ifw,' n &i &i &j &j &i &j &i *j &i')&cd

NB. SUBROUTINE ZGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
zgerucd=: (lib,' zgeru_ ',ifw,' n &i &i &j &j &i &j &i *j &i')&cd

NB. SUBROUTINE DSYR ( UPLO, N, ALPHA, X, INCX, A, LDA )
dsyrcd=: (lib,' dsyr_ ',ifw,' n &c &i &d &d &i *d &i')&cd

NB. SUBROUTINE DSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
dsyr2cd=: (lib,' dsyr2_ ',ifw,' n &c &i &d &d &i &d &i *d &i')&cd

NB. SUBROUTINE ZHER ( UPLO, N, ALPHA, X, INCX, A, LDA )
zhercd=: (lib,' zher_ ',ifw,' n &c &i &d &j &i *j &i')&cd

NB. SUBROUTINE ZHER2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
zher2cd=: (lib,' zher2_ ',ifw,' n &c &i &j &j &i &j &i *j &i')&cd

NB. ---------------------------------------------------------
NB. Level 3 BLAS

NB. SUBROUTINE _GEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
dgemmcd=: (lib,' dgemm_ ',ifw,' n &c &c &i &i &i &d &d &i &d &i &d *d &i')&cd
zgemmcd=: (lib,' zgemm_ ',ifw,' n &c &c &i &i &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE _SYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
dsymmcd=: (lib,' dsymm_ ',ifw,' n &c &c &i &i &d &d &i &d &i &d *d &i')&cd
zsymmcd=: (lib,' zsymm_ ',ifw,' n &c &c &i &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE ZHEMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
zhemmcd=: (lib,' zhemm_ ',ifw,' n &c &c &i &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE _TRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )
dtrmmcd=: (lib,' dtrmm_ ',ifw,' n &c &c &c &c &i &i &d &d &i *d &i')&cd
ztrmmcd=: (lib,' ztrmm_ ',ifw,' n &c &c &c &c &i &i &j &j &i *j &i')&cd

NB. SUBROUTINE _SYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
dsyrkcd=: (lib,' dsyrk_ ',ifw,' n &c &c &i &i &d &d &i &d *d &i')&cd
zsyrkcd=: (lib,' zsyrk_ ',ifw,' n &c &c &i &i &j &j &i &j *j &i')&cd

NB. SUBROUTINE ZHERK ( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
zherkcd=: (lib,' zherk_ ',ifw,' n &c &c &i &i &d &j &i &d *j &i')&cd

NB. SUBROUTINE _SYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
dsyr2kcd=: (lib,' dsyr2k_ ',ifw,' n &c &c &i &i &d &d &i &d &i &d *d &i')&cd
zsyr2kcd=: (lib,' zsyr2k_ ',ifw,' n &c &c &i &i &j &j &i &j &i &j *j &i')&cd

NB. SUBROUTINE ZHER2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
zher2kcd=: (lib,' zher2k_ ',ifw,' n &c &c &i &i &j &j &i &j &i &d *j &i')&cd

NB. SUBROUTINE _TRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )
dtrsmcd=: (lib,' dtrsm_ ',ifw,' n &c &c &c &c &i &i &d &d &i *d &i')&cd
ztrsmcd=: (lib,' ztrsm_ ',ifw,' n &c &c &c &c &i &i &j &j &i *j &i')&cd

erase 'lib ifw'

NB. ---------------------------------------------------------
NB. Utilities

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
basiccs0=: issquare_jlapack2_@(0&{::)
basiccs1=: issquare_jlapack2_@(1&{::)
basiccs3=: issquare_jlapack2_@(3&{::)
basiccs4=: issquare_jlapack2_@(4&{::)
basiccs5=: issquare_jlapack2_@(5&{::)
NB. - compare shapes
basiccmp=: -:/@($L:0)@(1 2&{)

NB. modify
NB. - conjugate under ISO specified
basiccj0=: 1      &(+&.> upd)
basiccj1=: 0 1 3  &(+&.> upd)
basiccj2=: 0 2 4 5&(+&.> upd)
NB. - swap elements
basicswp=: (< 1 2)&C.
