NB. <this file purpose; brief description>
NB.
NB. <name>        <brief description>
NB. ...
NB. <name>        <brief description>
NB.
NB. <testxxname>  Test <name> by <matrix of type xx>
NB. ...
NB. <testyyname>  Test <name> by <matrix of type yy>
NB. <testname>    Adv. to make verb to test <names> by
NB.               matrix of generator and shape given
NB.
NB. Version: <n.n.n> <yyyy-mm-dd>
NB.
NB. Copyright <yyyy> <full name>
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

<possible section delimiters list>

NB. delimiter                                                  levels

NB. #########################################################   0  0  0
NB. #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#   1
NB. # # # # # # # # # # # # # # # # # # # # # # # # # # # # #   2  1
NB. =========================================================   3  2  1  0  0  0
NB. =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=   4
NB. = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   5  3
NB. *********************************************************   6  4  2  1
NB. *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*   7
NB. * * * * * * * * * * * * * * * * * * * * * * * * * * * * *   8  5
NB. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++   9  6  3  2  1
NB. +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+  10
NB. + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  11  7
NB. ---------------------------------------------------------  12  8  4  3  2  1
NB. - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  13  9

<end of possible section delimiters list>

NB. =========================================================
NB. Concepts <optional section>
NB.
NB. Terms: <optional section>
NB.   <Any text>
NB.
NB. Notation: <optional section>
NB.   <Any text>
NB.
NB. Conventions: <optional section>
NB.   <Any text>
NB.
NB. Assertions: <optional section>
NB.   <Any text>
NB.
NB. Examples: <optional section>
NB.   <Any text>
NB.
NB. Notes: <optional section>
NB.   <Any text>
NB.
NB. TODO: <optional section>
NB.   <Any text>
NB.
NB. References: <optional section>
NB. [1] <reference>
NB. ...
NB. [x] <reference>

NB. =========================================================
NB. Configuration

coclass 'mt'
<code (optional)>

NB. =========================================================
NB. Local definitions <optional section>

<isolated definitions, optional section>

<localname>=: <definition>  NB. <individual description>
...
<localname>=: <definition>  NB. <individual description>

<OR>

NB. <multi-line...
NB. ...individual description>
<localname>=: <definition>

...

NB. <multi-line...
NB. ...individual description>
<localname>=: <definition>

<OR>

NB. <common description>
<localname>=: <definition>
...
<localname>=: <definition>

<end of optional isolated definitions section>

NB. ---------------------------------------------------------
NB. <grouped definitions description, optional section>

<localname>=: <definition>  NB. <individual description>
...
<localname>=: <definition>  NB. <individual description>

<OR>

NB. <multi-line...
NB. ...individual description>
<localname>=: <definition>

...

NB. <multi-line...
NB. ...individual description>
<localname>=: <definition>

<OR>

NB. <common description>
<localname>=: <definition>
...
<localname>=: <definition>

NB. ---------------------------------------------------------
NB. <LocalAdverb>
NB.
NB. Description: <optional section>
NB.   Adv. to make verb to <description>
NB. where <optional section>
NB.   <item 1> - <description>
NB.   <item 2> - <description>
NB.   ...
NB.
NB. Syntax: <optional section>
NB.   <o>=. [<x>] (<u> <LocalAdverb>) <y>
NB. where <optional section>
NB.   <u> - monad/dyad to <description>, is called as:
NB.           <o>=. [<x>] <u> <y>
NB.   <y> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <x> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <o> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.
NB. Formula: <optional section>
NB. - <description>:
NB.     <formula>
NB. ...
NB. - <description>:
NB.     <formula>
NB.
NB. Storage layout: <optional section>
NB.   <layout>
NB. where <optional section>
NB.   <description>
NB.
NB. Algorithm: <optional section>
NB.   In: <input nouns; comma-separated list>
NB.   Out: <output nouns; comma-separated list>
NB.   1) <1st step description>
NB.   ...
NB.   x) <last step description>
NB.
NB. Assertions (with appropriate comparison tolerance): <optional section>
NB.   <assertion>
NB.   ...
NB.   <assertion>
NB. where <optional section>
NB.   <copula>
NB.   ...
NB.   <copula>
NB.
NB. Examples: <optional section>
NB.    <command to execute>
NB. <output>
NB. ...
NB.    <command to execute>
NB. <output>
NB.
NB. Application: <optional section>
NB. - <description>
NB.     <command to execute>
NB. ...
NB. - <description>
NB.     <command to execute>
NB.
NB. Notes: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. TODO: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. References: <optional section>
NB. [1] <reference>
NB. ...
NB. [x] <reference>

<LocalAdverb>=: <definition>

NB. ---------------------------------------------------------
NB. <LocalConjunction>
NB.
NB. Description: <optional section>
NB.   Conj. to make verb to <description>
NB. where <optional section>
NB.   <item 1> - <description>
NB.   <item 2> - <description>
NB.   ...
NB.
NB. Syntax: <optional section>
NB.   <o>=. [<x>] (<u> <LocalConjunction> <v>) <y>
NB. where <optional section>
NB.   <u> - monad/dyad to <description>, is called as:
NB.           <o>=. [<x>] <u> <y>
NB.   <v> - monad/dyad to <description>, is called as:
NB.           <o>=. [<x>] <v> <y>
NB.   <y> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <x> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <o> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.
NB. Formula: <optional section>
NB. - <description>:
NB.     <formula>
NB. ...
NB. - <description>:
NB.     <formula>
NB.
NB. Storage layout: <optional section>
NB.   <layout>
NB. where <optional section>
NB.   <description>
NB.
NB. Algorithm: <optional section>
NB.   In: <input nouns; comma-separated list>
NB.   Out: <output nouns; comma-separated list>
NB.   1) <1st step description>
NB.   ...
NB.   x) <last step description>
NB.
NB. Assertions (with appropriate comparison tolerance): <optional section>
NB.   <assertion>
NB.   ...
NB.   <assertion>
NB. where <optional section>
NB.   <copula>
NB.   ...
NB.   <copula>
NB.
NB. Examples: <optional section>
NB.    <command to execute>
NB. <output>
NB. ...
NB.    <command to execute>
NB. <output>
NB.
NB. Application: <optional section>
NB. - <description>
NB.     <command to execute>
NB. ...
NB. - <description>
NB.     <command to execute>
NB.
NB. Notes: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. TODO: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. References: <optional section>
NB. [1] <reference>
NB. ...
NB. [x] <reference>

<LocalConjunction>=: <definition>

NB. ---------------------------------------------------------
NB. <LocalVerb>
NB.
NB. Description:
NB.   <description>
NB. where <optional section>
NB.   <item 1> - <description>
NB.   <item 2> - <description>
NB.   ...
NB.
NB. Syntax:
NB.   <o>=. [<x>] <LocalVerb> <y>
NB. where <optional section>
NB.   <y> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <x> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <o> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.
NB. Formula: <optional section>
NB. - <description>:
NB.     <formula>
NB. ...
NB. - <description>:
NB.     <formula>
NB.
NB. Storage layout: <optional section>
NB.   <layout>
NB. where <optional section>
NB.   <description>
NB.
NB. Algorithm: <optional section>
NB.   In: <input nouns; comma-separated list>
NB.   Out: <output nouns; comma-separated list>
NB.   1) <1st step description>
NB.   ...
NB.   x) <last step description>
NB.
NB. Assertions (with appropriate comparison tolerance): <optional section>
NB.   <assertion>
NB.   ...
NB.   <assertion>
NB. where <optional section>
NB.   <copula>
NB.   ...
NB.   <copula>
NB.
NB. Examples: <optional section>
NB.    <command to execute>
NB. <output>
NB. ...
NB.    <command to execute>
NB. <output>
NB.
NB. Application: <optional section>
NB. - <description>
NB.     <command to execute>
NB. ...
NB. - <description>
NB.     <command to execute>
NB.
NB. Notes: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. TODO: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. References: <optional section>
NB. [1] <reference>
NB. ...
NB. [x] <reference>

<LocalVerb>=: <definition>

NB. =========================================================
NB. Interface

<isolated definitions, optional section>

<name>=: <definition>  NB. <individual description>
...
<name>=: <definition>  NB. <individual description>

<OR>

NB. <multi-line...
NB. ...individual description>
<name>=: <definition>

...

NB. <multi-line...
NB. ...individual description>
<name>=: <definition>

<OR>

NB. <common description>
<name>=: <definition>
...
<name>=: <definition>

<end of optional isolated definitions section>

NB. ---------------------------------------------------------
NB. <grouped definitions description, optional section>

<name>=: <definition>  NB. <individual description>
...
<name>=: <definition>  NB. <individual description>

<OR>

NB. <multi-line...
NB. ...individual description>
<name>=: <definition>

...

NB. <multi-line...
NB. ...individual description>
<name>=: <definition>

<OR>

NB. <common description>
<name>=: <definition>
...
<name>=: <definition>

NB. ---------------------------------------------------------
NB. <Adverb>
NB.
NB. Description: <optional section>
NB.   Adv. to make verb to <description>
NB. where <optional section>
NB.   <item 1> - <description>
NB.   <item 2> - <description>
NB.   ...
NB.
NB. Syntax: <optional section>
NB.   <o>=. [<x>] (<u> <Adverb>) <y>
NB. where <optional section>
NB.   <u> - monad/dyad to <description>, is called as:
NB.           <o>=. [<x>] <u> <y>
NB.   <y> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <x> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <o> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.
NB. Formula: <optional section>
NB. - <description>:
NB.     <formula>
NB. ...
NB. - <description>:
NB.     <formula>
NB.
NB. Storage layout: <optional section>
NB.   <layout>
NB. where <optional section>
NB.   <description>
NB.
NB. Algorithm: <optional section>
NB.   In: <input nouns; comma-separated list>
NB.   Out: <output nouns; comma-separated list>
NB.   1) <1st step description>
NB.   ...
NB.   x) <last step description>
NB.
NB. Assertions (with appropriate comparison tolerance): <optional section>
NB.   <assertion>
NB.   ...
NB.   <assertion>
NB. where <optional section>
NB.   <copula>
NB.   ...
NB.   <copula>
NB.
NB. Examples: <optional section>
NB.    <command to execute>
NB. <output>
NB. ...
NB.    <command to execute>
NB. <output>
NB.
NB. Application: <optional section>
NB. - <description>
NB.     <command to execute>
NB. ...
NB. - <description>
NB.     <command to execute>
NB.
NB. Notes: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. TODO: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. References: <optional section>
NB. [1] <reference>
NB. ...
NB. [x] <reference>

<Adverb>=: <definition>

NB. ---------------------------------------------------------
NB. <Conjunction>
NB.
NB. Description: <optional section>
NB.   Conj. to make verb to <description>
NB. where <optional section>
NB.   <item 1> - <description>
NB.   <item 2> - <description>
NB.   ...
NB.
NB. Syntax: <optional section>
NB.   <o>=. [<x>] (<u> <Conjunction> <v>) <y>
NB. where <optional section>
NB.   <u> - monad/dyad to <description>, is called as:
NB.           <o>=. [<x>] <u> <y>
NB.   <v> - monad/dyad to <description>, is called as:
NB.           <o>=. [<x>] <v> <y>
NB.   <y> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <x> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <o> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.
NB. Formula: <optional section>
NB. - <description>:
NB.     <formula>
NB. ...
NB. - <description>:
NB.     <formula>
NB.
NB. Storage layout: <optional section>
NB.   <layout>
NB. where <optional section>
NB.   <description>
NB.
NB. Algorithm: <optional section>
NB.   In: <input nouns; comma-separated list>
NB.   Out: <output nouns; comma-separated list>
NB.   1) <1st step description>
NB.   ...
NB.   x) <last step description>
NB.
NB. Assertions (with appropriate comparison tolerance): <optional section>
NB.   <assertion>
NB.   ...
NB.   <assertion>
NB. where <optional section>
NB.   <copula>
NB.   ...
NB.   <copula>
NB.
NB. Examples: <optional section>
NB.    <command to execute>
NB. <output>
NB. ...
NB.    <command to execute>
NB. <output>
NB.
NB. Application: <optional section>
NB. - <description>
NB.     <command to execute>
NB. ...
NB. - <description>
NB.     <command to execute>
NB.
NB. Notes: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. TODO: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. References: <optional section>
NB. [1] <reference>
NB. ...
NB. [x] <reference>

<Conjunction>=: <definition>

NB. ---------------------------------------------------------
NB. <Verb>
NB.
NB. Description:
NB.   <description>
NB. where <optional section>
NB.   <item 1> - <description>
NB.   <item 2> - <description>
NB.   ...
NB.
NB. Syntax:
NB.   <o>=. [<x>] <Verb> <y>
NB. where <optional section>
NB.   <y> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <x> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.   <o> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.
NB. Formula: <optional section>
NB. - <description>:
NB.     <formula>
NB. ...
NB. - <description>:
NB.     <formula>
NB.
NB. Storage layout: <optional section>
NB.   <layout>
NB. where <optional section>
NB.   <description>
NB.
NB. Algorithm: <optional section>
NB.   In: <input nouns; comma-separated list>
NB.   Out: <output nouns; comma-separated list>
NB.   1) <1st step description>
NB.   ...
NB.   x) <last step description>
NB.
NB. Assertions (with appropriate comparison tolerance): <optional section>
NB.   <assertion>
NB.   ...
NB.   <assertion>
NB. where <optional section>
NB.   <copula>
NB.   ...
NB.   <copula>
NB.
NB. Examples: <optional section>
NB.    <command to execute>
NB. <output>
NB. ...
NB.    <command to execute>
NB. <output>
NB.
NB. Application: <optional section>
NB. - <description>
NB.     <command to execute>
NB. ...
NB. - <description>
NB.     <command to execute>
NB.
NB. Notes: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. TODO: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. References: <optional section>
NB. [1] <reference>
NB. ...
NB. [x] <reference>

<Verb>=: <definition>

NB. =========================================================
NB. Test suite

NB. ---------------------------------------------------------
NB. <testxxname>
NB.
NB. Description:
NB.   Test <description> by <matrix of type xx>
NB. where <optional section>
NB.   <item 1> - <description>
NB.   <item 2> - <description>
NB.   ...
NB.
NB. Syntax:
NB.   <testxxname> <y>
NB. where <optional section>
NB.   <y> - scalar/n-vector/m×n-matrix/sh-array <description>
NB.
NB. Formula: <optional section>
NB. - <description>:
NB.     <formula>
NB. ...
NB. - <description>:
NB.     <formula>
NB.
NB. Storage layout: <optional section>
NB.   <layout>
NB. where <optional section>
NB.   <description>
NB.
NB. Algorithm: <optional section>
NB.   In: <input nouns; comma-separated list>
NB.   Out: <output nouns; comma-separated list>
NB.   1) <1st step description>
NB.   ...
NB.   x) <last step description>
NB.
NB. Assertions (with appropriate comparison tolerance): <optional section>
NB.   <assertion>
NB.   ...
NB.   <assertion>
NB. where <optional section>
NB.   <copula>
NB.   ...
NB.   <copula>
NB.
NB. Examples: <optional section>
NB.    <command to execute>
NB. <output>
NB. ...
NB.    <command to execute>
NB. <output>
NB.
NB. Application: <optional section>
NB. - <description>
NB.     <command to execute>
NB. ...
NB. - <description>
NB.     <command to execute>
NB.
NB. Notes: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. TODO: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. References: <optional section>
NB. [1] <reference>
NB. ...
NB. [x] <reference>

<testxxname>=: <definition>

NB. ---------------------------------------------------------
NB. <testname>
NB.
NB. Description:
NB.   Adv. to make verb to test <names> by matrix of
NB.   generator and shape given
NB. where <optional section>
NB.   <item 1> - <description>
NB.   <item 2> - <description>
NB.   ...
NB.
NB. Syntax:
NB.   log=. (mkmat <testname>) (m,n)
NB. where <optional section>
NB.   mkmat - monad to generate a matrix; is called as:
NB.             mat=. mkmat (m,n)
NB.   (m,n) - 2-vector of integers, the shape of matrix mat
NB.   log   - 6-vector of boxes, test log, see test.ijs
NB.
NB. Formula: <optional section>
NB. - <description>:
NB.     <formula>
NB. ...
NB. - <description>:
NB.     <formula>
NB.
NB. Storage layout: <optional section>
NB.   <layout>
NB. where <optional section>
NB.   <description>
NB.
NB. Algorithm: <optional section>
NB.   In: <input nouns; comma-separated list>
NB.   Out: <output nouns; comma-separated list>
NB.   1) <1st step description>
NB.   ...
NB.   x) <last step description>
NB.
NB. Assertions (with appropriate comparison tolerance): <optional section>
NB.   <assertion>
NB.   ...
NB.   <assertion>
NB. where <optional section>
NB.   <copula>
NB.   ...
NB.   <copula>
NB.
NB. Examples: <optional section>
NB.    <command to execute>
NB. <output>
NB. ...
NB.    <command to execute>
NB. <output>
NB.
NB. Application: <optional section>
NB. - <description>
NB.     <command to execute>
NB. ...
NB. - <description>
NB.     <command to execute>
NB.
NB. Notes: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. TODO: <optional section>
NB. - <text>
NB. ...
NB. - <text>
NB.
NB. References: <optional section>
NB. [1] <reference>
NB. ...
NB. [x] <reference>

<testname>=: 1 : 'EMPTY [ <definition>'

NB. =========================================================
NB. Clean-up <optional section>

<code>
