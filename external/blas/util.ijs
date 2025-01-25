NB. Utilities
NB.
NB. ver  Get version string
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

coclass 'mtbla'

NB. =========================================================
NB. Interface

NB. Get version string
NB. str=. ver ''

ver=: 3 : 0
  ifw=. IFWIN # '+'
  if. 0 {:: dlsym LIB ; 'openblas_get_config' do.
    NB. OpenBLAS
    memr 0 _1 2 ,~ ('"',LIB,'" openblas_get_config >',ifw,' x') cd ''
  elseif. 0 {:: dlsym LIB ; 'ilaver_' do.
    NB. LAPACK
    3 }. ; (('.' , ":)L:0) ('"',LIB,'" ilaver_ ',ifw,' n *i *i *i') cd ((3 # < , 0))
  else.
    NB. the reference BLAS has no version identifier
    'unknown'
  end.
)
