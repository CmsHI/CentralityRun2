#!/bin/bash

#5 TeV 285090 Prompt Reco and EPOS
root -l -b -q 'compare_data_mc_1D.C+("root:://eoscms//eos/cms/store/group/phys_heavyions/dhangal/pr_forests/PAMinimumBias/v3/000/285/090/HiForest_", "/afs/cern.ch/work/q/qixu/public/centrality/Run2Prep/5TeV/mergeskim_pPb_EposMinBias_5TeV_8022_Forest_corrCentrTable_161108_122954.root")'

#5 TeV 285216 Prompt Reco and EPOS
#root -l -b -q 'compare_data_mc_1D.C+("", "", "(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1)")'

#5 TeV 285244

#5 TeV 285368
#5 TeV 285369
#5 TeV 285371
#5 TeV 285374

#8 TeV 285479
#8 TeV 285480

#8 TeV 285517 Express Stream and HIJING
#root -l -b -q 'compare_data_mc_1D.C+("", "")'

#8 TeV 285517 Express Stream and EPOS 
#root -l -b -q 'compare_data_mc_1D.C+("", "")'
