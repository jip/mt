NB. Verify mm verb
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
NB. Local definitions

NB. y is an array, direct  conversion (str <- arr) must fail
dir0=. 0:@ mm_mtmm_      :: 1

NB. y is a string, inverse conversion (str -> arr) must fail
inv0=. 0:@(mm_mtmm_^:_1) :: 1

NB. y is an array, cyclic conversion (arr <- str <- arr) must
NB. succeed: result must match with y by datatype, rank,
NB. shape and values
cyc1=. (-: *. -:&(3!:0)) ]&.mm_mtmm_

NB. =========================================================
NB. Verification suite

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. direct conversions which must fail

NB. error: 0-rank arrays aren't supported
NB. reason: J sentence creates a scalar which has rank<2
dir0 42

NB. error: 1-rank arrays aren't supported
NB. reason: J sentence creates a vector which has rank<2
dir0 1 2

NB. error: Matrix Market exchange formats support the 0 only as a sparse element
NB. reason: J sentence creates a sparse array with sparse element not equal to 0
dir0 1 (3 2 $ 0 1 1 0 1 2)} 1 $. 2 3 ; 0 1 ; _

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. inverse conversions which must fail

NB. error: line longer than 1024 bytes was detected
inv0 {{)n
%%matrixmarket matrix array integer general
0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9
0
}}

NB. error: not a Matrix Market exchange format
NB. reason: header row contains not 5 qualifiers
inv0 {{)n
%%matrixmarket matrix array integer general trash
1 1
11
}}

NB. error: not a Matrix Market exchange format
NB. reason: header row contains not 5 qualifiers
inv0 {{)n
%%matrixmarket matrix array integer 
1 1
11
}}

NB. error: banner '%%MatriksMarket' is not recognized
NB. reason: 1st qualifier in header row is not '%%MatrixMarket'
inv0 {{)n
%%MatriksMarket matrix array integer general
1 1
11
}}

NB. error: object 'vector' is not recognized
NB. reason: 2nd qualifier in header row is not 'matrix'
inv0 {{)n
%%matrixmarket vector array integer general
1
11
}}

NB. error: format 'raveled' is not recognized
NB. reason: 3rd qualifier (format) in header row is invalid
inv0 {{)n
%%matrixmarket matrix raveled integer general
1
11
}}

NB. error: 'boolean' field is not recognized
NB. reason: 4th qualifier (field) in header row is invalid
inv0 {{)n
%%matrixmarket matrix array boolean general
1
11
}}

NB. error: 'hermitian' field is not recognized
NB. reason: header row is invalid:
NB. - 4th and 5th qualifiers both specify symmetry
NB. - no qualifier specifying field
inv0 {{)n
%%matrixmarket matrix coordinate hermitian general
2 2 3
11
22
33
}}

NB. error: 'triangular' symmetry is not recognized
NB. reason: 5th qualifier (symmetry) in header row is invalid
inv0 {{)n
%%matrixmarket matrix array integer triangular
1
11
}}

NB. error: 'pattern' and 'skew-symmetric' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix coordinate pattern skew-symmetric
3 3 3
2 1
3 1
3 2
}}

NB. error: 'pattern' and 'hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix coordinate pattern hermitian
2 2 3
1 1
2 1
2 2
}}

NB. error: 'pattern' and 'skew-hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix coordinate pattern skew-hermitian
3 3 3
2 1
3 1
3 2
}}

NB. error: 'integer' and 'hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix coordinate integer hermitian
2 2 3
1 1 11
2 1 21
2 2 22
}}

NB. error: 'integer' and 'skew-hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix coordinate integer skew-hermitian
2 2 1
2 1 21
}}

NB. error: 'real' and 'hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix coordinate real hermitian
2 2 3
1 1 11.1
2 1 21.2
2 2 22.3
}}

NB. error: 'real' and 'skew-hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix coordinate real skew-hermitian
2 2 1
2 1 21.1
}}

NB. error: 'pattern' and 'skew-symmetric' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix array pattern skew-symmetric
3 3
1
0
1
}}

NB. error: 'pattern' and 'hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix array pattern hermitian
2 2
1
1
0
}}

NB. error: 'pattern' and 'skew-hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix array pattern skew-hermitian
3 3
1
0
1
}}

NB. error: 'integer' and 'hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix array integer hermitian
2 2
11
12
21
22
}}

NB. error: 'integer' and 'skew-hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix array integer skew-hermitian
2 2
21
}}

NB. error: 'real' and 'hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix array real hermitian
2 2
11.1
12.2
21.3
22.4
}}

NB. error: 'real' and 'skew-hermitian' qualifiers are incompatible
NB. reason: header row contains incompatible qualifiers
inv0 {{)n
%%matrixmarket matrix array real skew-hermitian
2 2
21.1
}}

NB. error: size values '1 2.3' must be integer'
NB. reason: size row is invalid (contains non-integer)
inv0 {{)n
%%matrixmarket matrix array integer general
1 2.3
11
22
}}

NB. error: size format is Dim1 ... DimN Len but '2 2' found
NB. reason: size row is invalid (contains less than 3 numbers)
inv0 {{)n
%%matrixmarket matrix coordinate pattern general
2 2
1
2
}}

NB. error: size format is Dim1 ... DimN but '2' found
NB. reason: size row is invalid (contains less than 2 numbers)
inv0 {{)n
%%matrixmarket matrix array integer general
2
11
22
}}

NB. error: 2 elements was expected, but 3 data rows found (1)
NB. reason: either the size would be '2 2 3' or there would be 2 data rows
inv0 {{)n
%%matrixmarket matrix coordinate pattern general
2 2 2
1 1
1 2
2 2
}}

NB. error: each data row must contain 2 ISO (1)
NB. reason: last data row contains insufficient quantity of numbers
inv0 {{)n
%%matrixmarket matrix coordinate pattern general
2 2 3
1 1
1 2
2
}}

NB. error: each data row must contain 2 ISO (1)
NB. reason: last data row contains excessive quantity of numbers
inv0 {{)n
%%matrixmarket matrix coordinate pattern general
2 2 3
1 1
1 2
2 2 2
}}

NB. error: each data row must contain 2 ISO and a single matrix entry (1)
NB. reason: last data row contains insufficient quantity of numbers
inv0 {{)n
%%matrixmarket matrix coordinate integer general
2 2 3
1 1 11
1 2 12
2 22
}}

NB. error: each data row must contain 2 ISO and a single matrix entry (1)
NB. reason: last data row contains excessive quantity of numbers
inv0 {{)n
%%matrixmarket matrix coordinate integer general
2 2 3
1 1 11
1 2 12
2 2 22 222
}}

NB. error: each data row must contain 2 ISO and a single matrix entry (1)
NB. reason: last data row contains insufficient quantity of numbers
inv0 {{)n
%%matrixmarket matrix coordinate real general
2 2 3
1 1 11.1
1 2 12.2
2 22.3
}}

NB. error: each data row must contain 2 ISO and a single matrix entry (1)
NB. reason: last data row contains excessive quantity of numbers
inv0 {{)n
%%matrixmarket matrix coordinate real general
2 2 3
1 1 11.1
1 2 12.2
2 2 22.3 22.3
}}

NB. error: each data row must contain 2 ISO and a real and imaginary part of single matrix entry (1)
NB. reason: last data row contains insufficient quantity of numbers
inv0 {{)n
%%matrixmarket matrix coordinate complex general
2 2 3
1 1 11.1 0.11
1 2 12.2 0.12
2 2 22.3
}}

NB. error: each data row must contain 2 ISO and a real and imaginary part of single matrix entry (1)
NB. reason: last data row contains excessive quantity of numbers
inv0 {{)n
%%matrixmarket matrix coordinate complex general
2 2 3
1 1 11.1 0.11
1 2 12.2 0.12
2 2 22.3 0.22 0.222
}}

NB. error: 'symmetric' matrix must have all dimensions the same
NB. reason: size row is invalid (matrix dimensions must be the same)
inv0 {{)n
%%matrixmarket matrix coordinate integer symmetric
2 3 3
1 1 11
2 1 21
2 2 22
}}

NB. error: 'skew-symmetric' matrix must have all dimensions the same
NB. reason: size row is invalid (matrix dimensions must be the same)
inv0 {{)n
%%matrixmarket matrix coordinate integer skew-symmetric
2 3 3
1 1 11
2 1 21
2 2 22
}}

NB. error: 'hermitian' matrix must have all dimensions the same
NB. reason: size row is invalid (matrix dimensions must be the same)
inv0 {{)n
%%matrixmarket matrix coordinate complex hermitian
2 3 3
1 1 11 0
2 1 21 0.21
2 2 22 0
}}

NB. error: 'skew-hermitian' matrix must have all dimensions the same
NB. reason: size row is invalid (matrix dimensions must be the same)
inv0 {{)n
%%matrixmarket matrix coordinate complex skew-hermitian
2 3 3
1 1 0 0.11
2 1 0 0.21
2 2 0 0.22
}}

NB. error: not more than 3 elements was expected, but 4 data rows found
NB. reason: excessive data rows found
inv0 {{)n
%%matrixmarket matrix coordinate integer symmetric
2 2 4
1 1 11
1 2 12
2 1 21
2 2 22
}}

NB. error: there are non-recognized values
NB. reason: last data row contains not a number
inv0 {{)n
%%matrixmarket matrix coordinate integer general
2 2 4
1 1 11
1 2 12
2 1 21
2 2 ab
}}

NB. error: each data row must contain 2 ISO (2)
NB. reason: last data row contains excessive number delimited by TABs span, which is considered as not a delimiter in Matrix Market exchange formats, but a delimiter in J
inv0 {{)n
%%matrixmarket matrix coordinate pattern general
2 2 4
1 1
1 2
2 1
2 2	3
}}

NB. error: each data row must contain 2 ISO and a single matrix entry (2)
NB. reason: last data row contains excessive number delimited by TABs span, which is considered as not a delimiter in Matrix Market exchange formats, but a delimiter in J
inv0 {{)n
%%matrixmarket matrix coordinate integer general
2 2 4
1 1 11
1 2 12
2 1 21
2 2 22	222
}}

NB. error: each data row must contain 2 ISO and a single matrix entry (2)
NB. reason: last data row contains excessive number delimited by TABs span, which is considered as not a delimiter in Matrix Market exchange formats, but a delimiter in J
inv0 {{)n
%%matrixmarket matrix coordinate real general
2 2 4
1 1 11.1
1 2 12.2
2 1 21.3
2 2 22.4	0.22
}}

NB. error: each data row must contain 2 ISO and a real and imaginary part of single matrix entry (2)
NB. reason: last data row contains excessive number delimited by TABs span, which is considered as not a delimiter in Matrix Market exchange formats, but a delimiter in J
inv0 {{)n
%%matrixmarket matrix coordinate complex general
2 2 4
1 1 11.1 0.11
1 2 12.2 0.12
2 1 21.3 0.21
2 2 22.4 0.22	0.222
}}

NB. error: indices must be 1-based
NB. reason: indices are 0-based
inv0 {{)n
%%matrixmarket matrix coordinate pattern general
2 2 4
0 0
0 1
1 0
1 1
}}

NB. error: index exceeding dimension is detected
NB. reason: last data row contains 2nd index (4) excessing 2nd dimension (3)
inv0 {{)n
%%matrixmarket matrix coordinate pattern general
2 3 4
1 1
1 2
2 1
2 4
}}

NB. error: index exceeding dimension is detected
NB. reason: last data row contains 2nd index (4) excessing 2nd dimension (3)
inv0 {{)n
%%matrixmarket matrix coordinate pattern general
2 3 4
1 1
1 2
2 1
2 4
}}

NB. error: integer data type was expected but real data type is detected
NB. reason: last data row contains real not integer value
inv0 {{)n
%%matrixmarket matrix coordinate integer general
2 2 4
1 1 11
1 2 12
2 1 21
2 2 22.1
}}

NB. error: diagonal values must be real
NB. reason: last data row contains value for diagonal element of Hermitian matrix which has non-zero imaginary part
inv0 {{)n
%%matrixmarket matrix coordinate complex hermitian
2 2 3
1 1 11.1   0.0
2 1 21.2 -0.21
2 2 22.3  0.22
}}

NB. error: diagonal values must be imaginary
NB. reason: last data row contains value for diagonal element of skew-Hermitian matrix which has non-zero real part
inv0 {{)n
%%matrixmarket matrix coordinate complex skew-hermitian
2 2 3
1 1   0.0 0.11
2 1 -21.1 0.21
2 2  22.2 0.22
}}

NB. error: entries above the main diagonal are detected
NB. reason: 2nd data row contains ISO for strict upper triangular part of matrix
inv0 {{)n
%%matrixmarket matrix coordinate integer symmetric
2 2 3
1 1 11
1 2 12
2 2 22
}}

NB. error: entries either on or above the main diagonal are detected
NB. reason: last data row contains ISO for upper triangular part of matrix
inv0 {{)n
%%matrixmarket matrix coordinate integer symmetric
2 2 1
1 2 12
}}

NB. error: entries either on or above the main diagonal are detected
NB. reason: last data row contains ISO for diagonal element of matrix
inv0 {{)n
%%matrixmarket matrix coordinate integer skew-symmetric
2 2 1
2 2 22
}}

NB. error: each data row must contain a real and imaginary part of single matrix entry (1)
NB. reason: last data row contains insufficient data (imagimary part omitted)
inv0 {{)n
%%matrixmarket matrix array complex general
2 2
11.1 0.11
12.2 0.12
22.3
}}

NB. error: 'symmetric' matrix must have all dimensions the same
NB. reason: size row is invalid (matrix dimensions must be the same)
inv0 {{)n
%%matrixmarket matrix array integer symmetric
2 3
11
21
22
}}

NB. error: 'skew-symmetric' matrix must have all dimensions the same
NB. reason: size row is invalid (matrix dimensions must be the same)
inv0 {{)n
%%matrixmarket matrix array integer skew-symmetric
2 3
-21
0
32
}}

NB. error: 'hermitian' matrix must have all dimensions the same
NB. reason: size row is invalid (matrix dimensions must be the same)
inv0 {{)n
%%matrixmarket matrix array complex hermitian
2 3
11.1    0
21.2 0.21
22.3 0
}}

NB. error: 'skew-hermitian' matrix must have all dimensions the same
NB. reason: size row is invalid (matrix dimensions must be the same)
inv0 {{)n
%%matrixmarket matrix coordinate complex skew-hermitian
2 3
   0 0.11
21.1 0.21
   0 0.22
}}

NB. error: 4 elements was expected, but 3 data rows found
NB. reason: either the size would be '1 3' or '3 1' instead of '2 2', or there would be 4 data rows
inv0 {{)n
%%matrixmarket matrix array integer general
2 2
11
12
21
}}

NB. error: there are non-recognized values
NB. reason: last data row contains not a number
inv0 {{)n
%%matrixmarket matrix array integer general
2 2
11
12
21
ab
}}

NB. error: each data row must contain a real and imaginary part of single matrix entry (2)
NB. reason: last data row contains excessive number delimited by TABs span, which is considered as not a delimiter in Matrix Market exchange formats, but a delimiter in J
inv0 {{)n
%%matrixmarket matrix array complex general
2 2
11.1 0.11
12.2 0.12
21.3 0.21
22.4 0.22	0.222
}}

NB. error: each data row must contain just a single matrix entry
NB. reason: last data row contains excessive number delimited by TABs span, which is considered as not a delimiter in Matrix Market exchange formats, but a delimiter in J
inv0 {{)n
%%matrixmarket matrix array real general
2 2
11.1
12.2
21.3
22.4	0.22
}}

NB. error: integer data type was expected but real data type is detected
NB. reason: last data row contains real not integer value
inv0 {{)n
%%matrixmarket matrix array integer general
2 2
11
12
21
22.1
}}

NB. error: diagonal values must be real
NB. reason: last data row contains value for diagonal element of Hermitian matrix which has non-zero imaginary part
inv0 {{)n
%%matrixmarket matrix array complex hermitian
2 2
11.1    0
21.2 0.21
22.3 0.22
}}

NB. error: diagonal values must be imaginary
NB. reason: last data row contains value for diagonal element of skew-Hermitian matrix which has non-zero real part
inv0 {{)n
%%matrixmarket matrix array complex skew-hermitian
2 2
   0 0.11
21.1 0.21
22.2 0.22
}}

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. cyclic conversions which must succeed

NB. ---------------------------------------------------------
NB. input without edge cases

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. coordinate

NB. type: coordinate pattern general
NB. rank: 2 (matrix)
NB. shape: 2 3
NB. lenght: 5
cyc1 1 (5 2 $ 0 0 0 1 0 2 1 0 1 2)} 1 $. 2 3 ; 0 1 ; 0

NB. type: coordinate pattern general
NB. rank: 3 (brick)
NB. shape: 2 3 4
NB. lenght: 9
cyc1 1 (9 3 $ 0 0 0 0 0 1 0 0 2 0 1 0 0 1 2 1 0 0 1 0 2 1 1 0 1 1 1)} 1 $. 2 3 4 ; 0 1 2 ; 0

NB. type: coordinate pattern symmetric
NB. rank: 2 (matrix)
NB. shape: 2 2
NB. lenght: 2
cyc1 1 (3 2 $ 0 0 0 1 1 0)} 1 $. 2 2 ; 0 1 ; 0

NB. type: coordinate pattern symmetric
NB. rank: 3 (brick)
NB. shape: 2 2 2
NB. lenght: 3
cyc1 1 (5 3 $ 0 0 0 0 0 1 0 1 0 1 0 0 1 1 1)} 1 $. 2 2 2 ; 0 1 2 ; 0

NB. type: coordinate integer general
NB. rank: 2 (matrix)
NB. shape: 2 3
NB. lenght: 5
cyc1 11 12 13 21 23 (5 2 $ 0 0 0 1 0 2 1 0 1 2)} 1 $. 2 3 ; 0 1 ; 00

NB. type: coordinate integer symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
NB. lenght: 5
cyc1 11 21 31 21 22 31 33 (7 2 $ 0 0 0 1 0 2 1 0 1 1 2 0 2 2)} 1 $. 3 3 ; 0 1 ; 00

NB. type: coordinate integer symmetric
NB. rank: 3 (brick)
NB. shape: 3 3 3
NB. lenght: 9
cyc1 111 121 131 121 122 131 133 121 122 122 222 232 232 233 131 133 232 233 133 233 333 (21 3 $ 0 0 0 0 0 1 0 0 2 0 1 0 0 1 1 0 2 0 0 2 2 1 0 0 1 0 1 1 1 0 1 1 1 1 1 2 1 2 1 1 2 2 2 0 0 2 0 2 2 1 1 2 1 2 2 2 0 2 2 1 2 2 2)} 1 $. 3 3 3 ; 0 1 2 ; 00

NB. type: coordinate integer skew-symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
NB. lenght: 2
cyc1 21 _21 _32 32 (4 2 $ 0 1 1 0 1 2 2 1)} 1 $. 3 3 ; 0 1 ; 00

NB. type: coordinate integer skew-symmetric
NB. rank: 3 (brick)
NB. shape: 4 4 4
NB. lenght: 4
cyc1 _132 142 132 _143 _142 143 132 _142 _132 _243 142 243 _132 143 132 243 _143 _243 142 _143 _142 _243 143 243 (24 3 $ 0 1 2 0 1 3 0 2 1 0 2 3 0 3 1 0 3 2 1 0 2 1 0 3 1 2 0 1 2 3 1 3 0 1 3 2 2 0 1 2 0 3 2 1 0 2 1 3 2 3 0 2 3 1 3 0 1 3 0 2 3 1 0 3 1 2 3 2 0 3 2 1)} 1 $. 4 4 4 ; 0 1 2 ; 00

NB. type: coordinate real general
NB. rank: 2 (matrix)
NB. shape: 2 3
NB. lenght: 5
cyc1 11.1 12.2 13.3 21.4 _23.5 (5 2 $ 0 0 0 1 0 2 1 0 1 2)} 1 $. 2 3

NB. type: coordinate real general
NB. rank: 3 (brick)
NB. shape: 2 3 4
NB. lenght: 9
cyc1 111.1 112.2 113.3 121.4 123.5 211.6 213.7 221.8 222.9 (9 3 $ 0 0 0 0 0 1 0 0 2 0 1 0 0 1 2 1 0 0 1 0 2 1 1 0 1 1 1)} 1 $. 2 3 4

NB. type: coordinate real symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
NB. lenght: 5
cyc1 11.1 21.2 31.4 21.2 22.3 31.4 33.5 (7 2 $ 0 0 0 1 0 2 1 0 1 1 2 0 2 2)} 1 $. 3 3

NB. type: coordinate real symmetric
NB. rank: 3 (brick)
NB. shape: 3 3 3
NB. lenght: 9
cyc1 111.1 121.2 131.4 121.2 122.3 131.4 133.5 121.2 122.3 122.3 222.6 232.7 232.7 233.8 131.4 133.5 232.7 233.8 133.5 233.8 333.9 (21 3 $ 0 0 0 0 0 1 0 0 2 0 1 0 0 1 1 0 2 0 0 2 2 1 0 0 1 0 1 1 1 0 1 1 1 1 1 2 1 2 1 1 2 2 2 0 0 2 0 2 2 1 1 2 1 2 2 2 0 2 2 1 2 2 2)} 1 $. 3 3 3

NB. type: coordinate real skew-symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
NB. lenght: 2
cyc1 _21.1 21.1 32.2 _32.2 (4 2 $ 0 1 1 0 1 2 2 1)} 1 $. 3 3

NB. type: coordinate real skew-symmetric
NB. rank: 3 (brick)
NB. shape: 4 4 4
NB. lenght: 4
cyc1 _132.1 _142.2 132.1 _143.3 142.2 143.3 132.1 142.2 _132.1 243.4 _142.2 _243.4 _132.1 143.3 132.1 _243.4 _143.3 243.4 _142.2 _143.3 142.2 243.4 143.3 _243.4 (24 3 $ 0 1 2 0 1 3 0 2 1 0 2 3 0 3 1 0 3 2 1 0 2 1 0 3 1 2 0 1 2 3 1 3 0 1 3 2 2 0 1 2 0 3 2 1 0 2 1 3 2 3 0 2 3 1 3 0 1 3 0 2 3 1 0 3 1 2 3 2 0 3 2 1)} 1 $. 4 4 4

NB. type: coordinate complex general
NB. rank: 2 (matrix)
NB. shape: 2 3
NB. lenght: 5
cyc1 11.1j0.11 12.2j0.12 13.3j0.13 21.4j0.21 23.5j0.23 (5 2 $ 0 0 0 1 0 2 1 0 1 2)} 1 $. 2 3 ; 0 1 ; 0j0

NB. type: coordinate complex general
NB. rank: 3 (brick)
NB. shape: 2 3 4
NB. lenght: 9
cyc1 111.1j0.111 112.2j0.112 113.3j0.113 121.4j0.121 123.5j0.123 211.6j0.211 213.7j0.213 221.8j0.221 222.9j0.222 (9 3 $ 0 0 0 0 0 1 0 0 2 0 1 0 0 1 2 1 0 0 1 0 2 1 1 0 1 1 1)} 1 $. 2 3 4 ; 0 1 2 ; 0j0

NB. type: coordinate complex symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
NB. lenght: 5
cyc1 11.1j0.11 21.2j0.21 31.4j0.31 21.2j0.21 22.3j0.22 31.4j0.31 33.5j0.33 (7 2 $ 0 0 0 1 0 2 1 0 1 1 2 0 2 2)} 1 $. 3 3 ; 0 1 ; 0j0

NB. type: coordinate complex symmetric
NB. rank: 3 (brick)
NB. shape: 3 3 3
NB. lenght: 9
cyc1 111.1j0.111 121.2j0.121 131.4j0.131 121.2j0.121 122.3j0.122 131.4j0.131 133.5j0.133 121.2j0.121 122.3j0.122 122.3j0.122 222.6j0.222 232.7j0.232 232.7j0.232 233.8j0.233 131.4j0.131 133.5j0.133 232.7j0.232 233.8j0.233 133.5j0.133 233.8j0.233 333.9j0.333 (21 3$0 0 0 0 0 1 0 0 2 0 1 0 0 1 1 0 2 0 0 2 2 1 0 0 1 0 1 1 1 0 1 1 1 1 1 2 1 2 1 1 2 2 2 0 0 2 0 2 2 1 1 2 1 2 2 2 0 2 2 1 2 2 2)} 1 $. 3 3 3 ; 0 1 2 ; 0j0

NB. type: coordinate complex skew-symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
NB. lenght: 2
cyc1 _21.1j_0.21 21.1j0.21 _32.2j0.22 32.2j_0.22 (4 2 $ 0 1 1 0 1 2 2 1)} 1 $. 3 3 ; 0 1 ; 0j0

NB. type: coordinate complex skew-symmetric
NB. rank: 3 (brick)
NB. shape: 4 4 4
NB. lenght: 4
cyc1 _132.1j_0.132 _142.2j_0.142 132.1j0.132 _143.3j_0.143 142.2j0.142 143.3j0.143 132.1j0.132 142.2j0.142 _132.1j_0.132 _243.4j_0.243 _142.2j_0.142 243.4j0.243 _132.1j_0.132 143.3j0.143 132.1j0.132 243.4j0.243 _143.3j_0.143 _243.4j_0.243 _142.2j_0.142 _143.3j_0.143 142.2j0.142 _243.4j_0.243 143.3j0.143 243.4j0.243 (24 3 $ 0 1 2 0 1 3 0 2 1 0 2 3 0 3 1 0 3 2 1 0 2 1 0 3 1 2 0 1 2 3 1 3 0 1 3 2 2 0 1 2 0 3 2 1 0 2 1 3 2 3 0 2 3 1 3 0 1 3 0 2 3 1 0 3 1 2 3 2 0 3 2 1)} 1 $. 4 4 4 ; 0 1 2 ; 0j0

NB. type: coordinate complex hermitian
NB. rank: 2 (matrix)
NB. shape: 3 3
NB. lenght: 5
cyc1 11.1 21.2j_0.21 31.4j_0.31 21.2j0.21 22.3 31.4j0.31 33.5 (7 2 $ 0 0 0 1 0 2 1 0 1 1 2 0 2 2)} 1 $. 3 3 ; 0 1 ; 0j0

NB. type: coordinate complex hermitian
NB. rank: 3 (brick)
NB. shape: 3 3 3
NB. lenght: 9
cyc1 111 121 131 121 122 132j_0.132 131 132j0.132 121 122 132j0.132 122 222 232 132j_0.132 232 233 131 132j_0.132 132j0.132 232 233 233 333 (24 3 $ 0 0 0 0 0 1 0 0 2 0 1 0 0 1 1 0 1 2 0 2 0 0 2 1 1 0 0 1 0 1 1 0 2 1 1 0 1 1 1 1 1 2 1 2 0 1 2 1 1 2 2 2 0 0 2 0 1 2 1 0 2 1 1 2 1 2 2 2 1 2 2 2)} 1 $. 3 3 3 ; 0 1 2 ; 0j0

NB. type: coordinate complex skew-hermitian
NB. rank: 2 (matrix)
NB. shape: 3 3
NB. lenght: 5
cyc1 0j0.11 _21j0.21 _31j0.31 21j0.21 0j0.22 _32j0.32 31j0.31 32j0.32 (8 2 $ 0 0 0 1 0 2 1 0 1 1 1 2 2 0 2 1)} 1 $. 3 3 ; 0 1 ; 0j0

NB. type: coordinate complex skew-hermitian
NB. rank: 3 (brick)
NB. shape: 3 3 3
NB. lenght: 9
cyc1 0j0.111 0j0.121 0j0.131 0j0.121 0j0.122 _132j0.132 0j0.131 132j0.132 0j0.121 0j0.122 132j0.132 0j0.122 0j0.222 0j0.232 _132j0.132 0j0.232 0j0.233 0j0.131 _132j0.132 132j0.132 0j0.232 0j0.233 0j0.233 0j0.333 (24 3 $ 0 0 0 0 0 1 0 0 2 0 1 0 0 1 1 0 1 2 0 2 0 0 2 1 1 0 0 1 0 1 1 0 2 1 1 0 1 1 1 1 1 2 1 2 0 1 2 1 1 2 2 2 0 0 2 0 1 2 1 0 2 1 1 2 1 2 2 2 1 2 2 2)} 1 $. 3 3 3 ; 0 1 2 ; 0j0

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. array

NB. type: array pattern general
NB. rank: 2 (matrix)
NB. shape: 2 3
cyc1 1 (5 2 $ 0 0 0 1 0 2 1 0 1 2)} 2 3 $ 0

NB. type: array pattern general
NB. rank: 3 (brick)
NB. shape: 2 3 4
cyc1 1 (9 3 $ 0 0 0 0 0 1 0 0 2 0 1 0 0 1 2 1 0 0 1 0 2 1 1 0 1 1 1)} 2 3 4 $ 0

NB. type: array pattern symmetric
NB. rank: 2 (matrix)
NB. shape: 2 2
cyc1 1 (3 2 $ 0 0 0 1 1 0)} 2 2 $ 0

NB. type: array pattern symmetric
NB. rank: 3 (brick)
NB. shape: 2 2 2
cyc1 1 (5 3 $ 0 0 0 0 0 1 0 1 0 1 0 0 1 1 1)} 2 2 2 $ 0

NB. type: array integer general
NB. rank: 2 (matrix)
NB. shape: 2 3
cyc1 10 20 +/ 1 2 3

NB. type: array integer general
NB. rank: 3 (brick)
NB. shape: 2 3 4
cyc1 100 200 +/ 10 20 30 +/ 1 2 3 4

NB. type: array integer symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
cyc1 3 3 $ 11 21 31 21 22 32 31 32 33

NB. type: array integer symmetric
NB. rank: 3 (brick)
NB. shape: 4 4 4
cyc1 4 4 4 $ 111 121 131 141 121 122 132 142 131 132 133 143 141 142 143 144 121 122 132 142 122 222 232 242 132 232 233 243 142 242 243 244 131 132 133 143 132 232 233 243 133 233 333 343 143 243 343 344 141 142 143 144 142 242 243 244 143 243 343 344 144 244 344 444

NB. type: array integer skew-symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
cyc1 3 3 $ 0 21 31 _21 0 32 _31 _32 0

NB. type: array integer skew-symmetric
NB. rank: 3 (brick)
NB. shape: 4 4 4
cyc1 _132 _142 132 _143 142 143 132 142 _132 _243 _142 243 _132 143 132 243 _143 _243 _142 _143 142 _243 143 243 (24 3 $ 0 1 2 0 1 3 0 2 1 0 2 3 0 3 1 0 3 2 1 0 2 1 0 3 1 2 0 1 2 3 1 3 0 1 3 2 2 0 1 2 0 3 2 1 0 2 1 3 2 3 0 2 3 1 3 0 1 3 0 2 3 1 0 3 1 2 3 2 0 3 2 1)} 4 4 4 $ 00

NB. type: array real general
NB. rank: 2 (matrix)
NB. shape: 2 3
cyc1 2 3 $ 11.1 12.3 13.5 21.2 22.4 23.6

NB. type: array real general
NB. rank: 3 (brick)
NB. shape: 2 3 4
cyc1 2 3 4 $ 111.01 112.04 113.07 114.1 121.02 122.05 123.08 124.11 131.03 132.06 133.09 134.12 211.13 212.16 213.19 214.22 221.14 222.17 223.2 224.23 231.15 232.18 233.21 234.24

NB. type: array real symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
cyc1 3 3 $ 11.1 21.2 31.3 21.2 22.4 32.5 31.3 32.5 33.6

NB. type: array real symmetric
NB. rank: 3 (brick)
NB. shape: 4 4 4
cyc1 4 4 4 $ 111.01 121.02 131.03 141.04 121.02 122.05 132.06 142.07 131.03 132.06 133.08 143.09 141.04 142.07 143.09 144.1 121.02 122.05 132.06 142.07 122.05 222.11 232.12 242.13 132.06 232.12 233.14 243.15 142.07 242.13 243.15 244.16 131.03 132.06 133.08 143.09 132.06 232.12 233.14 243.15 133.08 233.14 333.17 343.18 143.09 243.15 343.18 344.19 141.04 142.07 143.09 144.1 142.07 242.13 243.15 244.16 143.09 243.15 343.18 344.19 144.1 244.16 344.19 444.2

NB. type: array real skew-symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
cyc1 _21.1 31.2 21.1 _32.3 _31.2 32.3 (6 2 $ 0 1 0 2 1 0 1 2 2 0 2 1)} 3 3 $ 0.0

NB. type: array real skew-symmetric
NB. rank: 3 (brick)
NB. shape: 4 4 4
cyc1 _132.1 _142.2 132.1 _143.3 142.2 143.3 132.1 142.2 _132.1 _243.4 _142.2 243.4 _132.1 143.3 132.1 243.4 _143.3 _243.4 _142.2 _143.3 142.2 _243.4 143.3 243.4 (24 3 $ 0 1 2 0 1 3 0 2 1 0 2 3 0 3 1 0 3 2 1 0 2 1 0 3 1 2 0 1 2 3 1 3 0 1 3 2 2 0 1 2 0 3 2 1 0 2 1 3 2 3 0 2 3 1 3 0 1 3 0 2 3 1 0 3 1 2 3 2 0 3 2 1)} 4 4 4 $ 0.0

NB. type: array complex general
NB. rank: 2 (matrix)
NB. shape: 2 3
cyc1 2 3 $ 11.1j0.11 12.3j0.12 13.5j0.13 21.2j0.21 22.4j0.22 23.6j0.23

NB. type: array complex general
NB. rank: 3 (brick)
NB. shape: 2 3 4
cyc1 2 3 4 $ 111.01j0.111 112.04j0.112 113.07j0.113 114.1j0.114 121.02j0.121 122.05j0.122 123.08j0.123 124.11j0.124 131.03j0.131 132.06j0.132 133.09j0.133 134.12j0.134 211.13j0.211 212.16j0.212 213.19j0.213 214.22j0.214 221.14j0.221 222.17j0.222 223.2j0.223 224.23j0.224 231.15j0.231 232.18j0.232 233.21j0.233 234.24j0.234

NB. type: array complex symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
cyc1 3 3 $ 11.1j0.11 21.2j0.21 31.3j0.31 21.2j0.21 22.4j0.22 32.5j0.32 31.3j0.31 32.5j0.32 33.6j0.33

NB. type: array complex symmetric
NB. rank: 3 (brick)
NB. shape: 4 4 4
cyc1 4 4 4 $ 111.01j0.111 121.02j0.121 131.03j0.131 141.04j0.141 121.02j0.121 122.05j0.122 132.06j0.132 142.07j0.142 131.03j0.131 132.06j0.132 133.08j0.133 143.09j0.143 141.04j0.141 142.07j0.142 143.09j0.143 144.1j0.144 121.02j0.121 122.05j0.122 132.06j0.132 142.07j0.142 122.05j0.122 222.11j0.222 232.12j0.232 242.13j0.242 132.06j0.132 232.12j0.232 233.14j0.233 243.15j0.243 142.07j0.142 242.13j0.242 243.15j0.243 244.16j0.244 131.03j0.131 132.06j0.132 133.08j0.133 143.09j0.143 132.06j0.132 232.12j0.232 233.14j0.233 243.15j0.243 133.08j0.133 233.14j0.233 333.17j0.333 343.18j0.343 143.09j0.143 243.15j0.243 343.18j0.343 344.19j0.344 141.04j0.141 142.07j0.142 143.09j0.143 144.1j0.144 142.07j0.142 242.13j0.242 243.15j0.243 244.16j0.244 143.09j0.143 243.15j0.243 343.18j0.343 344.19j0.344 144.1j0.144 244.16j0.244 344.19j0.344 444.2j0.444

NB. type: array complex skew-symmetric
NB. rank: 2 (matrix)
NB. shape: 3 3
cyc1 _21.1j_0.21 _31.2j_0.31 21.1j0.21 _32.3j_0.32 31.2j0.31 32.3j0.32 (6 2 $ 0 1 0 2 1 0 1 2 2 0 2 1)} 3 3 $ 0j0

NB. type: array complex skew-symmetric
NB. rank: 3 (brick)
NB. shape: 4 4 4
cyc1 _132.1j_0.132 _142.2j_0.142 132.1j0.132 _143.3j_0.143 142.2j0.142 143.3j0.143 132.1j0.132 142.2j0.142 _132.1j_0.132 _243.4j_0.243 _142.2j_0.142 243.4j0.243 _132.1j_0.132 143.3j0.143 132.1j0.132 243.4j0.243 _143.3j_0.143 _243.4j_0.243 _142.2j_0.142 _143.3j_0.143 142.2j0.142 _243.4j_0.243 143.3j0.143 243.4j0.243 (24 3 $ 0 1 2 0 1 3 0 2 1 0 2 3 0 3 1 0 3 2 1 0 2 1 0 3 1 2 0 1 2 3 1 3 0 1 3 2 2 0 1 2 0 3 2 1 0 2 1 3 2 3 0 2 3 1 3 0 1 3 0 2 3 1 0 3 1 2 3 2 0 3 2 1)} 4 4 4 $ 0j0

NB. type: array complex hermitian
NB. rank: 2 (matrix)
NB. shape: 3 3
cyc1 3 3 $ 11.1 21.2j_0.21 31.3j_0.31 21.2j0.21 22.4 32.5j_0.32 31.3j0.31 32.5j0.32 33.6

NB. type: array complex hermitian
NB. rank: 3 (brick)
NB. shape: 4 4 4
cyc1 4 4 4 $ 111.01 121.02 131.03 141.04 121.02 122.05 132.06j_0.132 142.07j_0.142 131.03 132.06j0.132 133.08 143.09j_0.143 141.04 142.07j0.142 143.09j0.143 144.1 121.02 122.05 132.06j0.132 142.07j0.142 122.05 222.11 232.12 242.13 132.06j_0.132 232.12 233.14 243.15j_0.243 142.07j_0.142 242.13 243.15j0.243 244.16 131.03 132.06j_0.132 133.08 143.09j0.143 132.06j0.132 232.12 233.14 243.15j0.243 133.08 233.14 333.17 343.18 143.09j_0.143 243.15j_0.243 343.18 344.19 141.04 142.07j_0.142 143.09j_0.143 144.1 142.07j0.142 242.13 243.15j_0.243 244.16 143.09j0.143 243.15j0.243 343.18 344.19 144.1 244.16 344.19 444.2

NB. type: array complex skew-hermitian
NB. rank: 2 (matrix)
NB. shape: 3 3
cyc1 3 3 $ 0j0.11 _21.2j0.21 _31.3j0.31 21.2j0.21 0j0.22 _32.5j0.32 31.3j0.31 32.5j0.32 0j0.33

NB. type: array complex skew-hermitian
NB. rank: 3 (brick)
NB. shape: 4 4 4
cyc1 4 4 4 $ 0j0.111 0j0.121 0j0.131 0j0.141 0j0.121 0j0.122 _132.06j0.132 _142.07j0.142 0j0.131 132.06j0.132 0j0.133 _143.09j0.143 0j0.141 142.07j0.142 143.09j0.143 0j0.144 0j0.121 0j0.122 132.06j0.132 142.07j0.142 0j0.122 0j0.222 0j0.232 0j0.242 _132.06j0.132 0j0.232 0j0.233 _243.15j0.243 _142.07j0.142 0j0.242 243.15j0.243 0j0.244 0j0.131 _132.06j0.132 0j0.133 143.09j0.143 132.06j0.132 0j0.232 0j0.233 243.15j0.243 0j0.133 0j0.233 0j0.333 0j0.343 _143.09j0.143 _243.15j0.243 0j0.343 0j0.344 0j0.141 _142.07j0.142 _143.09j0.143 0j0.144 142.07j0.142 0j0.242 _243.15j0.243 0j0.244 143.09j0.143 243.15j0.243 0j0.343 0j0.344 0j0.144 0j0.244 0j0.344 0j0.444

NB. ---------------------------------------------------------
NB. input with edge cases

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. coordinate

cyc1 $. 0 0   $ 0
cyc1 $. 0 0   $ 00
cyc1 $. 0 0   $ 0.0
cyc1 $. 0 0   $ 0j0
cyc1 $. 0 0 0 $ 0
cyc1 $. 0 0 0 $ 00
cyc1 $. 0 0 0 $ 0.0
cyc1 $. 0 0 0 $ 0j0
cyc1 $. 1 1   $ 0
cyc1 $. 1 1   $ 00
cyc1 $. 1 1   $ 0.0
cyc1 $. 1 1   $ 0j0
cyc1 $. 1 1 1 $ 0
cyc1 $. 1 1 1 $ 00
cyc1 $. 1 1 1 $ 0.0
cyc1 $. 1 1 1 $ 0j0
cyc1 $. 1 1   $ 1
cyc1 $. 1 1   $ 01
cyc1 $. 1 1   $ 1.0
cyc1 $. 1 1   $ 1j0
cyc1 $. 1 1 1 $ 1
cyc1 $. 1 1 1 $ 01
cyc1 $. 1 1 1 $ 1.0
cyc1 $. 1 1 1 $ 1j0
cyc1 $. 1 1   $ 12
cyc1 $. 1 1   $ 1.2
cyc1 $. 1 1   $ 1j2
cyc1 $. 1 1 1 $ 12
cyc1 $. 1 1 1 $ 1.2
cyc1 $. 1 1 1 $ 1j2
cyc1 $. 2 2   $ 0
cyc1 $. 2 2   $ 00
cyc1 $. 2 2   $ 0.0
cyc1 $. 2 2   $ 0j0
cyc1 $. 2 2 2 $ 0
cyc1 $. 2 2 2 $ 00
cyc1 $. 2 2 2 $ 0.0
cyc1 $. 2 2 2 $ 0j0
cyc1 $. 2 2   $ 1
cyc1 $. 2 2   $ 01
cyc1 $. 2 2   $ 1.0
cyc1 $. 2 2   $ 1j0
cyc1 $. 2 2 2 $ 1
cyc1 $. 2 2 2 $ 01
cyc1 $. 2 2 2 $ 1.0
cyc1 $. 2 2 2 $ 1j0
cyc1 $. 2 2   $ 12
cyc1 $. 2 2   $ 1.2
cyc1 $. 2 2   $ 1j2
cyc1 $. 2 2 2 $ 12
cyc1 $. 2 2 2 $ 1.2
cyc1 $. 2 2 2 $ 1j2
cyc1 $. 2 3   $ 0
cyc1 $. 2 3   $ 00
cyc1 $. 2 3   $ 0.0
cyc1 $. 2 3   $ 0j0
cyc1 $. 2 3 4 $ 0
cyc1 $. 2 3 4 $ 00
cyc1 $. 2 3 4 $ 0.0
cyc1 $. 2 3 4 $ 0j0
cyc1 $. 2 3   $ 1
cyc1 $. 2 3   $ 01
cyc1 $. 2 3   $ 1.0
cyc1 $. 2 3   $ 1j0
cyc1 $. 2 3 4 $ 1
cyc1 $. 2 3 4 $ 01
cyc1 $. 2 3 4 $ 1.0
cyc1 $. 2 3 4 $ 1j0
cyc1 $. 2 3   $ 12
cyc1 $. 2 3   $ 1.2
cyc1 $. 2 3   $ 1j2
cyc1 $. 2 3 4 $ 12
cyc1 $. 2 3 4 $ 1.2
cyc1 $. 2 3 4 $ 1j2

NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NB. array

cyc1    0 0   $ 0
cyc1    0 0   $ 00
cyc1    0 0   $ 0.0
cyc1    0 0   $ 0j0
cyc1    0 0 0 $ 0
cyc1    0 0 0 $ 00
cyc1    0 0 0 $ 0.0
cyc1    0 0 0 $ 0j0
cyc1    1 1   $ 0
cyc1    1 1   $ 00
cyc1    1 1   $ 0.0
cyc1    1 1   $ 0j0
cyc1    1 1 1 $ 0
cyc1    1 1 1 $ 00
cyc1    1 1 1 $ 0.0
cyc1    1 1 1 $ 0j0
cyc1    1 1   $ 1
cyc1    1 1   $ 01
cyc1    1 1   $ 1.0
cyc1    1 1   $ 1j0
cyc1    1 1 1 $ 1
cyc1    1 1 1 $ 01
cyc1    1 1 1 $ 1.0
cyc1    1 1 1 $ 1j0
cyc1    1 1   $ 12
cyc1    1 1   $ 1.2
cyc1    1 1   $ 1j2
cyc1    1 1 1 $ 12
cyc1    1 1 1 $ 1.2
cyc1    1 1 1 $ 1j2
cyc1    2 2   $ 0
cyc1    2 2   $ 00
cyc1    2 2   $ 0.0
cyc1    2 2   $ 0j0
cyc1    2 2 2 $ 0
cyc1    2 2 2 $ 00
cyc1    2 2 2 $ 0.0
cyc1    2 2 2 $ 0j0
cyc1    2 2   $ 1
cyc1    2 2   $ 01
cyc1    2 2   $ 1.0
cyc1    2 2   $ 1j0
cyc1    2 2 2 $ 1
cyc1    2 2 2 $ 01
cyc1    2 2 2 $ 1.0
cyc1    2 2 2 $ 1j0
cyc1    2 2   $ 12
cyc1    2 2   $ 1.2
cyc1    2 2   $ 1j2
cyc1    2 2 2 $ 12
cyc1    2 2 2 $ 1.2
cyc1    2 2 2 $ 1j2
cyc1    2 3   $ 0
cyc1    2 3   $ 00
cyc1    2 3   $ 0.0
cyc1    2 3   $ 0j0
cyc1    2 3 4 $ 0
cyc1    2 3 4 $ 00
cyc1    2 3 4 $ 0.0
cyc1    2 3 4 $ 0j0
cyc1    2 3   $ 1
cyc1    2 3   $ 01
cyc1    2 3   $ 1.0
cyc1    2 3   $ 1j0
cyc1    2 3 4 $ 1
cyc1    2 3 4 $ 01
cyc1    2 3 4 $ 1.0
cyc1    2 3 4 $ 1j0
cyc1    2 3   $ 12
cyc1    2 3   $ 1.2
cyc1    2 3   $ 1j2
cyc1    2 3 4 $ 12
cyc1    2 3 4 $ 1.2
cyc1    2 3 4 $ 1j2
