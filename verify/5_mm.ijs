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

NB. verification files directory
DATA_DIR=. '~addons/math/mt/verify/mm_data/'

NB. =========================================================
NB. Verification suite

NB. direct conversion which must fail
0:@ mm_mtmm_      :: 1 fread2_mtmm_ DATA_DIR , '103.ijs'
0:@ mm_mtmm_      :: 1 fread2_mtmm_ DATA_DIR , '104.ijs'
0:@ mm_mtmm_      :: 1 fread2_mtmm_ DATA_DIR , '105.ijs'

NB. inverse conversion which must fail
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '000.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '001.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '002.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '003.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '004.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '005.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '006.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '007.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '008.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '009.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '010.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '011.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '013.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '014.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '015.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '016.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '017.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '018.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '019.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '020.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '021.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '022.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '023.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '024.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '025.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '026.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '027.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '028.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '029.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '030.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '031.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '032.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '033.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '034.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '035.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '036.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '037.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '038.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '039.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '040.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '041.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '068.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '069.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '070.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '071.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '072.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '073.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '074.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '093.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '094.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '095.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '098.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '099.mm'
0:@(mm_mtmm_^:_1) :: 1 fread2_mtmm_ DATA_DIR , '100.mm'

NB. direct conversion which must succeed
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '012_apg_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '042_cpg_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '043_cpg_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '044_cps_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '045_cps_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '046_cig_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '047_aps_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '048_cis_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '049_cis_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '050_cik_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '051_cik_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '052_crg_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '053_crg_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '054_crs_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '055_crs_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '056_crk_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '057_crk_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '058_ccg_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '059_ccg_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '060_ccs_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '061_ccs_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '062_cck_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '063_cck_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '064_cch_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '065_cch_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '066_ccw_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '067_ccw_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '075_aig_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '076_aig_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '077_ais_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '078_ais_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '079_aik_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '080_aik_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '081_arg_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '082_arg_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '083_ars_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '084_ars_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '085_ark_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '086_ark_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '087_acg_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '088_acg_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '089_acs_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '090_acs_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '091_ack_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '092_ack_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '096_ach_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '097_ach_3'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '101_acw_2'
(,&'.mm' ((-: mm_mtmm_    )  ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '102_acw_3'

NB. inverse conversion which must succeed
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '012_apg_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '042_cpg_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '043_cpg_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '044_cps_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '045_cps_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '046_cig_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '047_aps_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '048_cis_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '049_cis_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '050_cik_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '051_cik_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '052_crg_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '053_crg_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '054_crs_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '055_crs_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '056_crk_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '057_crk_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '058_ccg_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '059_ccg_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '060_ccs_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '061_ccs_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '062_cck_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '063_cck_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '064_cch_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '065_cch_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '066_ccw_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '067_ccw_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '075_aig_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '076_aig_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '077_ais_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '078_ais_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '079_aik_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '080_aik_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '081_arg_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '082_arg_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '083_ars_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '084_ars_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '085_ark_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '086_ark_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '087_acg_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '088_acg_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '089_acs_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '090_acs_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '091_ack_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '092_ack_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '096_ach_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '097_ach_3'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '101_acw_2'
(,&'.mm' ((-: mm_mtmm_^:_1)~ ".)&fread2_mtmm_ ,&'.ijs') DATA_DIR , '102_acw_3'

NB. inverse-of-direct conversion which must succeed
NB. - dense
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    0 0   $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    0 0   $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    0 0   $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    0 0   $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    0 0 0 $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    0 0 0 $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    0 0 0 $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    0 0 0 $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1   $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1   $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1   $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1   $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1 1 $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1 1 $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1 1 $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1 1 $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1   $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1   $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1   $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1   $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1 1 $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1 1 $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1 1 $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1 1 $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1   $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1   $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1   $ 1j2
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1 1 $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1 1 $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    1 1 1 $ 1j2
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2   $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2   $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2   $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2   $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2 2 $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2 2 $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2 2 $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2 2 $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2   $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2   $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2   $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2   $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2 2 $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2 2 $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2 2 $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2 2 $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2   $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2   $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2   $ 1j2
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2 2 $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2 2 $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 2 2 $ 1j2
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3   $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3   $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3   $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3   $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3 4 $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3 4 $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3 4 $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3 4 $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3   $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3   $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3   $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3   $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3 4 $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3 4 $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3 4 $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3 4 $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3   $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3   $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3   $ 1j2
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3 4 $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3 4 $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_)    2 3 4 $ 1j2
NB. - sparse
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 0 0   $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 0 0   $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 0 0   $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 0 0   $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 0 0 0 $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 0 0 0 $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 0 0 0 $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 0 0 0 $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1   $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1   $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1   $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1   $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1 1 $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1 1 $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1 1 $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1 1 $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1   $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1   $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1   $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1   $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1 1 $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1 1 $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1 1 $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1 1 $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1   $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1   $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1   $ 1j2
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1 1 $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1 1 $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 1 1 1 $ 1j2
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2   $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2   $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2   $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2   $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2 2 $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2 2 $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2 2 $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2 2 $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2   $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2   $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2   $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2   $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2 2 $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2 2 $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2 2 $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2 2 $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2   $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2   $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2   $ 1j2
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2 2 $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2 2 $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 2 2 $ 1j2
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3   $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3   $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3   $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3   $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3 4 $ 0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3 4 $ 00
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3 4 $ 0.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3 4 $ 0j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3   $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3   $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3   $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3   $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3 4 $ 1
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3 4 $ 01
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3 4 $ 1.0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3 4 $ 1j0
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3   $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3   $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3   $ 1j2
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3 4 $ 12
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3 4 $ 1.2
((-: *. -:&(3!:0)) ]&.mm_mtmm_) $. 2 3 4 $ 1j2
