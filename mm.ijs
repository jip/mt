NB. Matrix Market exchange formats converter
NB.
NB. mm      Convert J numeric array to/from suitable Matrix
NB.         Market exchange format string
NB.
NB. testmm  Test mm verb
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
NB. Notation:
NB.   MM - Matrix Market exchange formats
NB.   mm - a MM facility implemented by this module
NB.
NB. Conventions:
NB. 1) MM allows the 0 only as a sparse element in sparse
NB.    matrices, this reduces an amount of exportable J arrays
NB. 2) mm supports arrays of any rank>1, this extends MM
NB.    limited by rank=2 only
NB. 3) mm supports arrays with skew-Hermitian symmetry, this
NB.    extends MM
NB. 4) mm supports 'array pattern' qualifiers combination,
NB.    this extends MM
NB. 5) arrays of shape (rank # n) where n<2 are always
NB.    considered as general i.e. without any symmetry
NB. 6) ±inf and nan aren't still supported in exporting J
NB.    array to MM format, this reduces MM
NB.
NB. TODO:
NB. - replace Format (":) by Format (8!:n) in arr->str
NB.   conversion when and if (8!:n) will be extended to
NB.   support ±inf and nan
NB. - add datatypes: extended, rational, quaternion,
NB.   octonion
NB.
NB. References:
NB. [1] Text File Formats
NB.     https://math.nist.gov/MatrixMarket/formats.html

NB. =========================================================
NB. Configuration

coclass 'mtmm'
coinsert 'mt'

NB. =========================================================
NB. Local definitions

NB. header parts allowed for...
BANNER=:     '%%matrixmarket'                                                             NB. ...banner
OBJECT=:     'matrix'                                                                     NB. ...object
FORMATS=:    'array'   ; 'coordinate'                                                     NB. ...format
FIELDS=:     'pattern' ; 'integer'   ; 'real'           ; 'complex'                       NB. ...field
SYMMETRIES=: 'general' ; 'symmetric' ; 'skew-symmetric' ; 'hermitian' ; 'skew-hermitian'  NB. ...symmetry

NB. count SPACE spans
NB. counter=. cspans string
NB. notes: based on optimized m69 from JPhrases 4B
cspans=: +/@(> |.!.0)@(' '&E.)

NB. mask=. trlxmask4iso iso
NB. for symmetric array (which has all dimensions the same),
NB. mask to select ISO pointing to...
trl0mask=: *./@(2 </\ _2&C.)@|:  NB. ...strict lower triangle
trlmask=: *./@(2 <:/\ _2&C.)@|:  NB. ...lower triangle

NB. bvect=. trlxchk4iso iso
NB. where bvect - rank-vector, boolean, contains 1s only if check is succeed
NB. check whether ISO are pointing to...
trl0chk=: *./@((2 </\ _2&C.)&.|:)  NB. ...strict lower triangle
trlchk=: *./@((2 <:/\ _2&C.)&.|:)  NB. ...lower triangle

NB. generate all #s in radix y
NB. reference: JPhrases 8A
odometer=: #: i.@(*/)

NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NB. flt staff

NB. ---------------------------------------------------------
NB. Boolean short circuit
NB. (evaluated from right to left)
NB. u and v are predicates
NB. let operators priority be: not > and > or

and=:    2 : '   0:`    u @.v'  NB.      u  and      v                uandv=.    u and    v
notAnd=: 2 : '   0:`(-.@u)@.v'  NB. (not u) and      v                notuandv=. u notAnd v
orNot=:  2 : '   1:`    u @.v'  NB.      u              or (not v)    uornotv=.  u orNot  v
andNot=: 2 : '   u `    0:@.v'  NB.      u  and (not v)               uandnotv=. u andNot v
or=:     2 : '   u `    1:@.v'  NB.      u              or      v     uorv=.     u or     v
notOr=:  2 : '-.@u `    1:@.v'  NB. (not u)             or      v     notuorv=.  u notOr  v

NB. like fread, but may throw 'file name' error
fread2=: (1!:1)@fboxname@boxopen@jpath

NB. predicate to check is array symmetric
NB. based on:
NB.   https://code.jsoftware.com/wiki/Essays/Symmetric_Array
NB.   Author: Roger Hui
issym=:    (-: 0&|:) or (3 > #@$) and (     -: _2&|:)

NB. predicate to check is array skew-symmetric
isskw=:    (-: 0&|:) or (3 > #@$) and (-    -: _2&|:)

NB. predicate to check is array Hermitian
ishmt=:    (-: 0&|:) or (3 > #@$) and (   + -: _2&|:)

NB. predicate to check is array skew-Hermitian
isskwhmt=: (-: 0&|:) or (3 > #@$) and (-@:+ -: _2&|:)

NB. end of flt staff
NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NB. ---------------------------------------------------------
NB. isosym
NB.
NB. Description:
NB.   Compute ISO nub to/from symmetric array for shape given
NB.
NB. Syntax:
NB.   arriso=. isosym     shape
NB.   iso=.    isosym^:_1 shape
NB. where
NB.   shape  -:(rank # length)  NB. an array's shape
NB.   arriso - shape-array, ISO, symmetric
NB.   iso    - (rank ([ ! <:@+) length)-vector, ISO nub
NB.            from ravelled symmetric array
NB.   rank   > 1, integer, an array's rank
NB.
NB. Assertions:
NB.      nub -: iso ({ ,) arr
NB. where
NB.      'rank length'=. 2 4
NB.      shape=. rank # length
NB.      nub=. 11 21 31 41 22 32 42 33 43 44
NB.      ] arriso=. isosym shape
NB.   0 1 2 3
NB.   1 4 5 6
NB.   2 5 7 8
NB.   3 6 8 9
NB.      ] arr=. arriso { nub
NB.   11 21 31 41
NB.   21 22 32 42
NB.   31 32 33 43
NB.   41 42 43 44
NB.      ] iso=. isosym^:_1 shape
NB.   0 4 8 12 5 9 13 10 14 15

isosym=: (3 : 0) :. (#. _2&C.&.|:@(# combrep {.))
  iso=. (</.~ /:~"1) odometer y            NB. ISO for all elements, grouped by sorted ISO
  vals=. (# S: 0 # i.@#) iso               NB. replicate value for each ISO
  iso=. ; iso
  vals iso} y $ 0
)

NB. ---------------------------------------------------------
NB. isoskw
NB.
NB. Description:
NB.   Compute ISO nub to/from skew-symmetric array for shape
NB.   given
NB.
NB. Syntax:
NB.   arriso=. isoskw     shape
NB.   iso=.    isoskw^:_1 shape
NB. where
NB.   shape  -:(rank # length)  NB. an array's shape
NB.   arriso - shape-array, ISO, skew-symmetric
NB.   iso    - (rank ! length)-vector, ISO nub
NB.            from ravelled skew-symmetric array
NB.   rank   ∈ [2, length], integer, an array's rank
NB.
NB. Assertions:
NB.      nub -: iso ({ ,) arr
NB. where
NB.      'rank length'=. 2 4
NB.      shape=. rank # length
NB.      nub=. 21 31 41 32 42 43
NB.      ] nubx=. 0 , (, -@|.) nub
NB.   0 21 31 41 32 42 43 _43 _42 _32 _41 _31 _21
NB.      ] arriso=. isoskw shape
NB.   0 _1 _2 _3
NB.   1  0 _4 _5
NB.   2  4  0 _6
NB.   3  5  6  0
NB.      ] arr=. arriso { nubx
NB.    0 _21 _31 _41
NB.   21   0 _32 _42
NB.   31  32   0 _43
NB.   41  42  43   0
NB.      ] iso=. isoskw^:_1 shape
NB.   4 8 12 9 13 14

isoskw=: (3 : 0) :. (#. _2&C.&.|:@(# comb {.))
  iso=. (# perm {.) y                      NB. ISO for all but diagonal elements, sorted by
                                           NB. ordered ISO, there will be (rank ! length)
                                           NB. groups of equal size (! rank)
  vals=. (# (!@[ # #\@i.@!) {.) y          NB. assign a natural number to each ISO group,
                                           NB. each ISO within group will have the same value
  par=. (C.!.2) (/:"1) iso                 NB. parity (is valid since there is no diagonals here)
  vals=. (          1 = par)} (,: -) vals  NB. negate evenly permuted values
  vals iso} y $ 00
)

NB. ---------------------------------------------------------
NB. isohmt
NB.
NB. Description:
NB.   Compute ISO nub from Hermitian or skew-Hermitian array
NB.   for shape given
NB.
NB. Syntax:
NB.   arriso=. isohmt     shape
NB.   iso=.    isohmt^:_1 shape
NB. where
NB.   shape  -:(rank # length)  NB. an array's shape
NB.   arriso - shape-array, ISO, skew-symmetric
NB.   iso    - (rank ([ ! <:@+) length)-vector, ISO nub
NB.            from ravelled Hermitian or skew-Hermitian
NB.            array
NB.   rank   > 1, integer, an array's rank
NB.
NB. Assertions:
NB.      nub -: iso ({ ,) arr
NB. where
NB.      'rank length'=. 2 4
NB.      shape=. rank # length
NB.      nub=. 11 2j1 3j1 4j1 22 3j2 4j2 33 4j3 44
NB.      ] nubx=. (, +@|.@}.) nub
NB.   11 2j1 3j1 4j1 22 3j2 4j2 33 4j3 44 44 4j_3 33 4j_2 3j_2 22 4j_1 3j_1 2j_1
NB.      ] arriso=. isohmt shape
NB.   0 _1 _2 _3
NB.   1  4 _5 _6
NB.   2  5  7 _8
NB.   3  6  8  9
NB.      ] arr=. arriso { nubx
NB.    11 2j_1 3j_1 4j_1
NB.   2j1   22 3j_2 4j_2
NB.   3j1  3j2   33 4j_3
NB.   4j1  4j2  4j3   44
NB.      ] iso=. isohmt^:_1 shape
NB.   0 4 8 12 5 9 13 10 14 15

isohmt=: (3 : 0) :. (isosym^:_1)
  iso=. (</.~ /:~"1) odometer y            NB. ISO for all elements, grouped by sorted ISO
  vals=. (# S: 0 # i.@#) iso               NB. replicate value for each ISO
  iso=. ; iso
  ndmask=. (# y) = (#@~."1) iso            NB. ISO for non-diagonals mask
  par=. (C.!.2) (/:"1) iso                 NB. parity (is valid only for non-diagonals)
  vals=. (ndmask *. 1 = par)} (,: -) vals  NB. negate evenly permuted non-diagonal values
  vals iso} y $ 00
)

NB. ---------------------------------------------------------
NB. mmic
NB.
NB. Description:
NB.   Import MM coordinate object
NB.
NB. Syntax:
NB.   'riso rdat'=. (le ; shape ; ioField ; ioSymmetry) mmic bdat
NB. where
NB.   le         ≥ 0, quantity of elements expected
NB.   shape      - r-vector, object's shape
NB.   ioField    ≥ 0, IO(field) in FIELDS
NB.   ioSymmetry ≥ 0, IO(symmetry) in SYMMETRIES
NB.   bdat       - boxed strings, data lines
NB.   riso       - raw ISO
NB.   rdat       - raw data for object
NB.   r          > 1, object's rank
NB.
NB. References:
NB. [1] Henry Rich. [Jbeta] S: results seems inconsistent,
NB.     also is L:
NB.     2022-12-28 03:58:23 UTC
NB.     https://www.jsoftware.com/pipermail/beta/2022-December/010685.html
NB. [2] Github / Jsoftware / jsource / issues / 194
NB.     Amend gives length error for sparse array
NB.     Igor Zhuravlov, 2024-04-22 07:02 UTC
NB.     https://github.com/jsoftware/jsource/issues/194

mmic=: 4 : 0
  'le shape ioField ioSymmetry'=. x
  rank=. # shape
  NB. check elements quantity presented
  lp=. # y
  ((": le) , ' elements was expected, but ' , (": lp) , ' data rows found (1)') assert le = lp
  NB. check max columns quantity presented for non-empty objects
  if. lp do.
    select. ioField
      case. 0 do.  NB. pattern
        ('each data row must contain just ' , (": rank) , ' ISO (1)') assert ((<: rank) = cspans) S: 0 y  NB. count SPACE spans
      case. 1 ; 2 do.  NB. integer or real
        ('each data row must contain just ' , (": rank) , ' ISO and a single matrix entry (1)') assert (rank = cspans) S: 0 y  NB. count SPACE spans
      case. 3 do.  NB. complex
        ('each data row must contain just ' , (": rank) , ' ISO and a real and imaginary part of single matrix entry (1)') assert ((>: rank) = cspans) S: 0 y  NB. count SPACE spans
    end.
  end.
  NB. check elements quantity depending on symmetry
  if. ioSymmetry do.
    NB. symmetry of any kind
    ('''' , (ioSymmetry {:: SYMMETRIES) , ''' matrix must have all dimensions the same') assert ({. = }.) shape
    if. 2 = ioSymmetry do.  NB. skew-symmetric
      lemax=. rank    !       {. shape
    else.  NB. symmetric, Hermitian or skew-Hermitian
      lemax=. rank ([ ! <:@+) {. shape
    end.
  else.
    NB. no symmetry
    lemax=. */ shape
  end.
  ('not more than ' , (": lemax) , ' elements was expected, but ' , (": lp) , ' data rows found') assert lemax >: lp  NB. some elements may be omitted
  NB. convert strings with data lines to J array
  'rp cp'=. 2 ({.!.1) $ y=. _. ". > y  NB. rows and columns presented, fill is required when y is empty
  'there are non-recognized values' assert -. isnan < y
  NB. ((": le) , ' elements was expected, but ' , (": rp) , ' data rows found (2)') assert le = rp  NB. how is this possible to violate?
  fret=. ''  NB. makes sense for complex field only
  se=. ioField {:: ((0 ; 00 ; 0.0 ; 0j0))  NB. sparse element
  NB. check columns quantity for non-empty objects
  if. rp do.
    select. ioField
      case. 0 do.  NB. pattern
        ('each data row must contain ' , (": rank) , ' ISO (2)') assert cp = rank
      case. 1 ; 2 do.  NB. integer or real
        ('each data row must contain ' , (": rank) , ' ISO and a single matrix entry (2)') assert cp = >: rank
      case. 3 do.  NB. complex
        ('each data row must contain ' , (": rank) , ' ISO and a real and imaginary part of single matrix entry (2)') assert cp = 2 + rank
        fret=. '' ; ((1 j. <: rank) , 1j1) # 1 1  NB. to separate ISO from values
    end.
    NB. extract ISO and values:
    NB. - pattern: iso_vector
    NB. - integer, real: iso_vector value_scalar
    NB. - complex: iso_vector re_of_value_scalar im_of_value_scalar
    'iso dat'=. (;&1)`(}:"1 ; {:"1)`(}:"1 ; {:"1)`(fret&(<`(<@:(j./"1));.1))@.ioField y
    iso=. JINT c. iso  NB. convert to integer type
    'indices must be 1-based' assert 0 (< ,) iso
    iso=. <: iso  NB. translate MM's 1-based ISO to J's 0-based ones
    'index exceeding dimension is detected' assert iso <"1 shape
  else.
    iso=. i. 0 , rank
    dat=. 0 {. se
  end.
  NB. simplify datatype if mismatches with ioField
  select. ioField
    case. 1 do.  NB. integer
      if. JINT ~: 3!:0 dat do.
        'non-integer values detected' assert (-: <.) dat
        dat=. JINT c. dat
      end.
    case. 2 do.  NB. real
      if. JFL ~: 3!:0 dat do.
        'non-real values detected' assert (-: +) dat
        dat=. JFL c. dat
      end.
    case. do.
  end.
  NB. restore elements known due to symmetry
  NB. restore ISO omitted
  count=. ! rank
  select. ioSymmetry
    case. 1 ; 3 ; 4 do.  NB. symmetric, Hermitian or skew-Hermitian
      if. 3 = ioSymmetry do.  NB. Hermitian
        'diagonal values must be real' assert 0 = 11 o. (rank ~: (#@~."1) iso) # dat
      elseif. 4 = ioSymmetry do.  NB. skew-Hermitian
        'diagonal values must be imaginary' assert 0 = 9 o. (rank ~: (#@~."1) iso) # dat
      end.
      NB. elements presented must have ISO satisfying:
      NB.   ISO[0] <= ISO[1] <= ... <= ISO[R-4] <= ISO[R-3] <= ISO[R-1] <= ISO[R-2]
      NB. - the last two elements are swapped since data
      NB.   provided is on and below the main diagonal
      NB. - ISO for non-diagonal elements are permuted from
      NB.   (/:~ ISO) by odd permutation
      NB. - ISO for diagonal elements are not a permutation
      NB.   so has parity=0
      'entries above the main diagonal are detected' assert trlchk iso
      NB. restore ISO omitted
      NB. - generate the set of permutations for each ISO
      NB.   - get nub since ISO may repeat e.g. 0 2 2
      NB.   - its cardinality are different
      NB. - box sets to avoid reshaping later
      iso=. (i. count) <@~.@:A."_ 1 iso
      counts=. # S: 0^:(*@#) iso  NB. see [1]
      iso=. ; iso
    case. 2 do.  NB. skew-symmetric
      NB. elements presented must have ISO satisfying:
      NB.   ISO[0] < ISO[1] < ... < ISO[R-4] < ISO[R-3] < ISO[R-1] < ISO[R-2]
      NB. - the last two elements are swapped since data
      NB.   provided is below the main diagonal
      NB. - every such ISO is permuted from (/:~ ISO) by odd
      NB.   permutation
      'entries either on or above the main diagonal are detected' assert trl0chk iso
      NB. restore ISO omitted
      NB. - generate the set of permutations for each ISO
      NB.   - its cardinality are the same
      NB. - merge sets (no reshaping)
      iso=. ,/ (i. count) (A."_ 1) iso
  end.
  NB. restore values omitted
  NB. - replicate for each set
  select. ioSymmetry
    case. 1 do.  NB. symmetric
      if. ioField do. NB. not pattern
        dat=. counts # dat
      end.
    case. 2 do.  NB. skew-symmetric
      NB. - negate values from strict upper triangle i.e. if
      NB.   the corresponding ISO is derived from (/:~ ISO)
      NB.   by even permutation
      dat=. (1 = (C.!.2) (/:"1) iso)} (,:    -) count  # dat
    case. 3 do.  NB. Hermitian
      NB. - conjugate values from strict upper triangle i.e.
      NB.   if the corresponding ISO is derived from
      NB.   (/:~ ISO) by even permutation
      dat=. (1 = (C.!.2) (/:"1) iso)} (,: +   ) counts # dat
    case. 4 do.  NB. skew-Hermitian
      NB. - conjugate and negate values from strict upper
      NB.   triangle i.e. if the corresponding ISO is derived
      NB.   from (/:~ ISO) by even permutation
      dat=. (1 = (C.!.2) (/:"1) iso)} (,: +@:-) counts # dat
  end.
  NB. place values at ISO positions in sparse array
  y=. 1 $. shape ; (i. rank) ; se
  if. # dat do.  NB. work-around [2]
    y=. dat iso} y
  end.
)

NB. ---------------------------------------------------------
NB. mmia
NB.
NB. Description:
NB.   Import MM array object
NB.
NB. Syntax:
NB.   rdat=. (shape ; ioField ; ioSymmetry) mmia bdat
NB. where
NB.   shape      - r-vector, object's shape
NB.   ioField    ≥ 0, IO(field) in FIELDS
NB.   ioSymmetry ≥ 0, IO(symmetry) in SYMMETRIES
NB.   bdat       - boxed strings, data lines
NB.   rdat       - raw data for object
NB.   r          > 1, object's rank

mmia=: 4 : 0
  'shape ioField ioSymmetry'=. x
  rank=. # shape
  lp=. # y
  NB. check max columns quantity presented
  if. (3 = ioField) *. 0 < lp do.  NB. complex
    'each data row must contain a real and imaginary part of single matrix entry (1)' assert (' '&e.S:0) y
  end.
  NB. check elements quantity depending on symmetry
  if. ioSymmetry do.
    NB. symmetry of any kind
    ('''' , (ioSymmetry {:: SYMMETRIES) , ''' matrix must have all dimensions the same') assert ({. = }.) shape
    if. 2 = ioSymmetry do.  NB. skew-symmetric
      le=. rank    !       {. shape
    else.  NB. symmetric, Hermitian or skew-Hermitian
      le=. rank ([ ! <:@+) {. shape
    end.
  else.
    NB. no symmetry
    le=. */ shape
  end.
  ((": le) , ' elements was expected, but ' , (": lp) , ' data rows found') assert le = lp  NB. all elements must be presented
  NB. convert strings with data lines to J array
  'rp cp'=. 2 ({.!.1) $ y=. _. ". > y  NB. rows and columns presented, cp is 2 for complex field and 1 otherwise
  'there are non-recognized values' assert -. isnan < y
  NB. ((": le) , ' elements was expected, but ' , (": rp) , ' elements found') assert le = rp  NB. how is this possible to violate?
  NB. check columns quantity
  if. 3 = ioField do.  NB. complex
    'each data row must contain a real and imaginary part of single matrix entry (2)' assert (2 = cp) +. (, 0) -: $ y
    y=. JCMPX c. j./"1 y  NB. if y was a 0×0-matrix then (j./"1 y) will have datatype boolean so explicit conversion may be needed
  else.  NB. pattern or integer or real
    'each data row must contain just a single matrix entry' assert 1 = cp
    NB. simplify datatype if mismatches with ioField
    select. ioField
      case. 0 do.  NB. pattern
        if. JB01 ~: 3!:0 y do.
          'non-pattern values detected' assert (-: 1&=) y
          y=. JB01 c. y
        end.
      case. 1 do.  NB. integer
        if. JINT ~: 3!:0 y do.
          'non-integer values detected' assert (-: <.) y
          y=. JINT c. y
        end.
      case. 2 do.  NB. real
        if. JFL ~: 3!:0 y do.
          'non-real values detected' assert (-: +) y
          y=. JFL c. y
        end.
    end.
  end.
  NB. restore elements known due to symmetry
  select. ioSymmetry
    case. 0 do.  NB. general
      y=. (|:"2) (_2 C. shape) ($ ,) y
    case. 1 do.  NB. symmetric
      iso=. isosym shape
      y=. iso { y
    case. 2 do.  NB. skew-symmetric
      iso=. isoskw shape
      y=. iso { 0 , (, -@|.) y
    case. 3 do.  NB. Hermitian
      mask=. rank > (#@~."1) (_2&C.&.|:) (# combrep {.) shape
      'diagonal values must be real' assert 0 = 11 o. mask # y
      iso=. isohmt shape
      y=. iso { (, +@|.@}.) y
    case. 4 do.  NB. skew-Hermitian
      mask=. rank > (#@~."1) (_2&C.&.|:) (# combrep {.) shape
      'diagonal values must be imaginary' assert 0 = 9 o. mask # y
      iso=. isohmt shape
      y=. iso { (, -@:+@|.@}.) y
  end.
  y
)

NB. =========================================================
NB. Interface

NB. ---------------------------------------------------------
NB. mm
NB.
NB. Description:
NB.   Convert J numeric array to/from suitable MM string
NB.
NB. Syntax:
NB.   str=. mm     arr
NB.   arr=. mm^:_1 str
NB. where
NB.   str - string, arr in MM
NB.   arr - r-rank array, numeric
NB.   r   > 1, array's rank
NB.
NB. Assertions:
NB.   ((-: *. -:&(3!:0)) ]&.mm) arr
NB.
NB. Application:
NB.   str=. mm_mtmm_ arr
NB.   str fwrite '~temp/AM.mm'
NB.   arr=. mm_mtmm_^:_1 fread '~temp/AM.mm'
NB.
NB. Notes:
NB. - there is a J bug with locale in obverse [1], so
NB.   inversion uses locatives
NB.
NB. References:
NB. [1] https://code.jsoftware.com/wiki/System/Interpreter/Bugs/Errors#Obverse_in_Locale

mm=: (3 : 0) :. (3 : 0)
  NB. str->arr
  y=. CRLF cutl_mt_ y  NB. cut by spans of CR and LF
  'line longer than 1024 bytes was detected' assert_mt_ ((1024 >: #) S: 0) y
  header=. cut_mt_ tolower 0 {:: y  NB. to lower case, then cut by SPACE spans
  y=. (#~ ('%' ~: {.) S: 0) y  NB. remove header and comments
  y=. (#~ a:&~:) dltb L: 0 y   NB. remove empty lines
  'not a Matrix Market exchange format' assert_mt_ 5 = # header
  ('banner '''   , (0 {:: header) , ''' is not recognized') assert_mt_ BANNER_mtmm_ -: 0 {:: header
  ('object '''   , (1 {:: header) , ''' is not recognized') assert_mt_ OBJECT_mtmm_ -: 1 {:: header
  ioFormat=.   FORMATS_mtmm_    i. 2 { header
  ioField=.    FIELDS_mtmm_     i. 3 { header
  ioSymmetry=. SYMMETRIES_mtmm_ i. 4 { header
  ('format '''   , (2 {:: header) , ''' is not recognized') assert_mt_ ioFormat   < # FORMATS_mtmm_
  ('field '''    , (3 {:: header) , ''' is not recognized') assert_mt_ ioField    < # FIELDS_mtmm_
  ('symmetry ''' , (4 {:: header) , ''' is not recognized') assert_mt_ ioSymmetry < # SYMMETRIES_mtmm_
  ('''' , (ioField {:: FIELDS_mtmm_) , ''' and ''' , (ioSymmetry {:: SYMMETRIES_mtmm_) , ''' qualifiers are incompatible') assert_mt_ ioSymmetry <: 1 2 4 {~ 0 2 I. ioField
  size=. ". 0 {:: y
  ('size values ''' , (": size) , ''' must be integer') assert_mt_ (3!:0 size) e. JB01 , JINT
  y=. }. y
  if. ioFormat do.  NB. coordinate
    ('size format is Dim1 ... DimN Len but ''' , (": size) , ''' found') assert_mt_ 2 < # size
    'shape le'=. (}: ; {:) size  NB. shape, quantity expected
    y=. (le ; shape ; ioField ; ioSymmetry) mmic_mtmm_ y
  else.  NB. array
    ('size format is Dim1 ... DimN but ''' , (": size) , ''' found') assert_mt_ 1 < # size
    y=. (size ; ioField ; ioSymmetry) mmia_mtmm_ y
  end.
  y
)
  NB. arr->str
  rank=. # shape=. $ y  NB. rank and shape presented
  ((": rank) , '-rank arrays aren''t supported') assert 1 < rank
  'ioFormat ioField'=. 2 4 #: (JB01 , JINT , JFL , JCMPX , 1024 4096 8192 16384) i. 3!:0 y
  ioSymmetry=. shape 0:`((0:`((0 4 {~ isskwhmt)`3:@.ishmt)@.(3 = ioField))`2:@.isskw)@.(*@ioField)`1:@.issym@]`0:@.(({. +./@:~: or (2 > [) }.)@[) y
  if. ioFormat do.
    NB. coordinate
    'Matrix Market exchange formats support the 0 only as a sparse element' assert 0 = 3 $. y
    y=. 8 $. y    NB. remove sparse elements
    iso=. 4 $. y  NB. generate all ISO
    NB. compose data
    y=. (>: iso) [`(,. 5&$.)`(,. 5&$.)`(,. +.@(5&$.))@.ioField y
    NB. filter out repeating elements known due to symmetry
    y=. y [`(#~ trlmask)`(#~ trl0mask)`(#~ trlmask)`(#~ trlmask)@.ioSymmetry iso
  else.
    NB. array
    NB. compose data
    y=. +.^:(3 = ioField) , (|:"2)^:(0 = ioSymmetry) y
    NB. filter out repeating elements known due to symmetry
    y=. y [`({~ isosym^:_1)`({~ isoskw^:_1)`({~ isohmt^:_1)`({~ isohmt^:_1)@.ioSymmetry shape
  end.
  size=. ": (# y) ,~^:ioFormat shape  NB. append a quantity presented if coordinate format
  format=.   ioFormat   {:: FORMATS
  field=.    ioField    {:: FIELDS
  symmetry=. ioSymmetry {:: SYMMETRIES
  str=. BANNER , ' ' , OBJECT , ' ' , format , ' ' , field , ' ' , symmetry , LF , size , LF , , (((=&'_')`(,:&'-')}) ":!.(IF64 { 9 17) ,. y) ,. LF
)

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. testmm
NB.
NB. Description:
NB.   Test mm by matrix with the following structure:
NB.   - general (not structured)
NB.   - symmetric
NB.   - skew-symmetric
NB.   - Hermitian
NB.   - skew-Hermitian
NB.   and with the following sparsity:
NB.   - dense (non sparse)
NB.   - sparse
NB.
NB. Syntax:
NB.   log=. testmm A
NB. where
NB.   A   - array of any rank > 1
NB.   log - 6-vector of boxes, test log
NB.
NB. Notes:
NB. - berr shows boolean 'is matched exactly'

testmm=: 3 : 0
  rcondA=. gecon1 y

  log=.          ('mm_mtmm_'     tmonad (]       `(mm_mtmm_^:_1)`(rcondA"_)`nan`-:)) y
  log=. log lcat ('mm_mtmm_^:_1' tmonad (mm_mtmm_`]             `(rcondA"_)`nan`-:)) y
)

NB. ---------------------------------------------------------
NB. testmm_mt_
NB.
NB. Description:
NB.   Adv. to make verb to test mm by matrices of generator
NB.   and shape given
NB.
NB. Syntax:
NB.   log=. (mkmat testmm_mt_) (m,n)
NB. where
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log
NB.
NB. Application:
NB. - test by random square boolean matrix:
NB.     log=. ?@$&2 testmm_mt_ 15 15
NB. - test by random integer matrix:
NB.     log=. ?@$&100 testmm_mt_ 10 15
NB. - test by random real matrix:
NB.     log=. ?@$&0 testmm_mt_ 15 10
NB. - test by random square complex matrix:
NB.     log=. (gemat_mt_ j. gemat_mt_) testmm_mt_ 10 10

testmm_mt_=: 1 : 'lcat_mt_@((nolog_mt_`(lcat_mt_@(testmm_mtmm_@sy4gel_mt_`(testmm_mtmm_@ss4gel_mt_)`(testmm_mtmm_@he4gel_mt_)`(testmm_mtmm_@sh4gel_mt_)`:0))@.(=/@$) ,:~ testmm_mtmm_)@(u spmat_mt_ 0.25) ,~ (nolog_mt_`(lcat_mt_@(testmm_mtmm_@sy4gel_mt_`(testmm_mtmm_@ss4gel_mt_)`(testmm_mtmm_@he4gel_mt_)`(testmm_mtmm_@sh4gel_mt_)`:0))@.(=/@$) ,:~ testmm_mtmm_)@u)'
