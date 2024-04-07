NB. Interface to BLIS
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
coinsert 'mttst mt'

NB. =========================================================
NB. Includes

require                 'math/mt/test/util'       NB. Tests' utilities
require                 'math/mt/test/blis/init'  NB. Init LIB_mtbli_
require^:(*@#@".@'LIB') 'math/mt/test/blis/util'  NB. BLIS Utilities
require^:(*@#@".@'LIB') 'math/mt/test/blis/api'   NB. API definitions
