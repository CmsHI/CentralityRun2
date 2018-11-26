#ifndef VAR_HEADER_H
#define VAR_HEADER_H

//Headers{{{
#include <TROOT.h>
#include <TString.h>
//}}}

const Int_t NVar = 13;
TString Vars[NVar] = {"hiHF", "hiHFplus", "hiHFplusEta4", "hiHFminus",
							"hiHFminusEta4", "hiHFECut", "hiHFECutPlus",
							"hiHFECutMinus", "hiEB", "hiEE", "hiNtracks",
							"hiHFhit", "hiNpix"};
TString VarName[NVar] = {"HF", "HFplus", "HFplusEta4", "HFminus",
								"HFminusEta4", "HFECut", "HFECutplus",
								"HFECutminus", "EB", "EE", "Ntracks", "HFhit",
								"Npix"};
//number of bins and max value of variable for collision events
const Int_t nBinC[NVar] = {600, 300, 150, 300, 150, 600, 300, 300, 300,
									300, 400, 3600, 1000};
const Double_t VarMaxC[NVar] = {6000, 3000, 1500, 3000, 1500, 6000, 3000,
										3000, 3000, 3000, 4000, 180000, 100000};
const Double_t NormRangeCut[NVar] = {500, 150, 100, 150, 100, 500, 150, 150,
												300, 300, 300, 18000, 8000};
const Double_t VarW[NVar] = {4, 2, 2, 2, 2, 4, 2, 2, 2, 2, 4, 5, 10};

//number of bins and max value of variable for noise events
const Int_t nBinN[NVar] = {50, 50, 25, 50, 25, 50, 50, 50, 50, 50, 20,
									100, 50};
const Double_t VarMaxN[NVar] = {50, 50, 25, 50, 25, 50, 50, 50, 50, 50,
										20, 1000, 500};

const Int_t NSel = 44;//2018PbPb
//const Int_t NSel = 4;//tmp
//const Int_t NSel = 16;//2015PbPb

//Name of selection 2018 PbPb{{{
TString SelName[NSel] = {"", "2Trk", "PV", "CC", "Coin2th3", "Coin3th3",
								"Coin4th3", "Coin5th3", "Coin2th4", "Coin3th4",
								"Coin4th4", "Coin5th4", "2Trk_CC",
								"2Trk_Coin3th3", "2Trk_Coin4th3", "2Trk_Coin5th3",
								"2Trk_Coin2th4", "2Trk_Coin3th4", "2Trk_Coin4th4",
								"PV_CC", "PV_Coin3th3", "PV_Coin4th3",
								"PV_Coin5th3", "PV_Coin2th4", "PV_Coin3th4",
								"PV_Coin4th4", "CC_Coin3th3", "CC_Coin4th3",
								"CC_Coin5th3", "CC_Coin2th4", "CC_Coin3th4",
								"CC_Coin4th4", "2Trk_CC_Coin3th3",
								"2Trk_CC_Coin4th3", "2Trk_CC_Coin5th3",
								"2Trk_CC_Coin2th4", "2Trk_CC_Coin3th4",
								"2Trk_CC_Coin4th4", "PV_CC_Coin3th3",
								"PV_CC_Coin4th3", "PV_CC_Coin5th3",
								"PV_CC_Coin2th4", "PV_CC_Coin3th4",
								"PV_CC_Coin4th4"};
//}}}

//Selection cut 2018 PbPb{{{
TString Selection[NSel] = {
									"1",
									"(nTrk>=2)",
									"(pprimaryVertexFilter==1)",
									"(pclusterCompatibilityFilter==1)",
									"(phfCoincFilter2Th3==1)",
									"(phfCoincFilter3Th3==1)",
									"(phfCoincFilter4Th3==1)",
									"(phfCoincFilter5Th3==1)",
									"(phfCoincFilter2Th4==1)",
									"(phfCoincFilter3Th4==1)",
									"(phfCoincFilter4Th4==1)",
									"(phfCoincFilter5Th4==1)",
									"(nTrk>=2&&pclusterCompatibilityFilter==1)",
									"(nTrk>=2&&phfCoincFilter3Th3==1)",
									"(nTrk>=2&&phfCoincFilter4Th3==1)",
									"(nTrk>=2&&phfCoincFilter5Th3==1)",
									"(nTrk>=2&&phfCoincFilter2Th4==1)",
									"(nTrk>=2&&phfCoincFilter3Th4==1)",
									"(nTrk>=2&&phfCoincFilter4Th4==1)",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1)",
									"(pprimaryVertexFilter==1&&phfCoincFilter3Th3==1)",
									"(pprimaryVertexFilter==1&&phfCoincFilter4Th3==1)",
									"(pprimaryVertexFilter==1&&phfCoincFilter5Th3==1)",
									"(pprimaryVertexFilter==1&&phfCoincFilter2Th4==1)",
									"(pprimaryVertexFilter==1&&phfCoincFilter3Th4==1)",
									"(pprimaryVertexFilter==1&&phfCoincFilter4Th4==1)",
									"(pclusterCompatibilityFilter==1&&phfCoincFilter3Th3==1)",
									"(pclusterCompatibilityFilter==1&&phfCoincFilter4Th3==1)",
									"(pclusterCompatibilityFilter==1&&phfCoincFilter5Th3==1)",
									"(pclusterCompatibilityFilter==1&&phfCoincFilter2Th4==1)",
									"(pclusterCompatibilityFilter==1&&phfCoincFilter3Th4==1)",
									"(pclusterCompatibilityFilter==1&&phfCoincFilter4Th4==1)",
									"(nTrk>=2&&pclusterCompatibilityFilter==1&&phfCoincFilter3Th3==1)",
									"(nTrk>=2&&pclusterCompatibilityFilter==1&&phfCoincFilter4Th3==1)",
									"(nTrk>=2&&pclusterCompatibilityFilter==1&&phfCoincFilter5Th3==1)",
									"(nTrk>=2&&pclusterCompatibilityFilter==1&&phfCoincFilter2Th4==1)",
									"(nTrk>=2&&pclusterCompatibilityFilter==1&&phfCoincFilter3Th4==1)",
									"(nTrk>=2&&pclusterCompatibilityFilter==1&&phfCoincFilter4Th4==1)",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1&&phfCoincFilter3Th3==1)",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1&&phfCoincFilter4Th3==1)",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1&&phfCoincFilter5Th3==1)",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1&&phfCoincFilter2Th4==1)",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1&&phfCoincFilter3Th4==1)",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1&&phfCoincFilter4Th4==1)"};
//}}}

/*
//Name of selection 2018 PbPb{{{
TString SelName[NSel] = {
								"2Trk_CC_Coin2th4", "2Trk_CC_Coin3th4",
								"2Trk_CC_Coin4th4", "PV_CC_Coin3th3"
								};
//}}}

//Selection cut 2018 PbPb{{{
TString Selection[NSel] = {
									"(nTrk>=2&&pclusterCompatibilityFilter==1&&phfCoincFilter2Th4==1)",
									"(nTrk>=2&&pclusterCompatibilityFilter==1&&phfCoincFilter3Th4==1)",
									"(nTrk>=2&&pclusterCompatibilityFilter==1&&phfCoincFilter4Th4==1)",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1&&phfCoincFilter3Th3==1)"
									};
//}}}
*/
/*
//Name of selection 2015 PbPb{{{
TString SelName[NSel] = {"", "MB", "PV", "CC", "Coin", "MB_PV", "MB_CC",
								"MB_Coin", "PV_CC", "PV_Coin", "CC_Coin",
								"MB_PV_CC", "MB_PV_Coin", "MB_CC_Coin", 
								"PV_CC_Coin", "MB_PV_CC_Coin"};
//}}}

//Selection cut 2015 PbPb{{{
TString Selection[NSel] = {"1",
									"(HLT_HIL1MinimumBiasHF2AND_part1_v1==1||HLT_HIL1MinimumBiasHF2AND_part2_v1==1||HLT_HIL1MinimumBiasHF2AND_part3_v1==1||HLT_HIL1MinimumBiasHF2AND_part4_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part6_v1==1||HLT_HIL1MinimumBiasHF2AND_part7_v1==1||HLT_HIL1MinimumBiasHF2AND_part8_v1==1||HLT_HIL1MinimumBiasHF2AND_part9_v1==1)",
									"(pprimaryVertexFilter==1)",
									//"pBeamScrapingFilter==1",
									"(pclusterCompatibilityFilter==1)",
									"(phfCoincFilter3==1)",
									"( (HLT_HIL1MinimumBiasHF2AND_part1_v1==1||HLT_HIL1MinimumBiasHF2AND_part2_v1==1||HLT_HIL1MinimumBiasHF2AND_part3_v1==1||HLT_HIL1MinimumBiasHF2AND_part4_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part6_v1==1||HLT_HIL1MinimumBiasHF2AND_part7_v1==1||HLT_HIL1MinimumBiasHF2AND_part8_v1==1||HLT_HIL1MinimumBiasHF2AND_part9_v1==1)&&(pprimaryVertexFilter==1) )",
									"( (HLT_HIL1MinimumBiasHF2AND_part1_v1==1||HLT_HIL1MinimumBiasHF2AND_part2_v1==1||HLT_HIL1MinimumBiasHF2AND_part3_v1==1||HLT_HIL1MinimumBiasHF2AND_part4_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part6_v1==1||HLT_HIL1MinimumBiasHF2AND_part7_v1==1||HLT_HIL1MinimumBiasHF2AND_part8_v1==1||HLT_HIL1MinimumBiasHF2AND_part9_v1==1)&&(pclusterCompatibilityFilter==1) )",
									"( (HLT_HIL1MinimumBiasHF2AND_part1_v1==1||HLT_HIL1MinimumBiasHF2AND_part2_v1==1||HLT_HIL1MinimumBiasHF2AND_part3_v1==1||HLT_HIL1MinimumBiasHF2AND_part4_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part6_v1==1||HLT_HIL1MinimumBiasHF2AND_part7_v1==1||HLT_HIL1MinimumBiasHF2AND_part8_v1==1||HLT_HIL1MinimumBiasHF2AND_part9_v1==1)&&(phfCoincFilter3==1) )",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1)",
									"(pprimaryVertexFilter==1&&phfCoincFilter3==1)",
									"(pclusterCompatibilityFilter==1&&phfCoincFilter3==1)",
									"( (HLT_HIL1MinimumBiasHF2AND_part1_v1==1||HLT_HIL1MinimumBiasHF2AND_part2_v1==1||HLT_HIL1MinimumBiasHF2AND_part3_v1==1||HLT_HIL1MinimumBiasHF2AND_part4_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part6_v1==1||HLT_HIL1MinimumBiasHF2AND_part7_v1==1||HLT_HIL1MinimumBiasHF2AND_part8_v1==1||HLT_HIL1MinimumBiasHF2AND_part9_v1==1)&&(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1) )",
									"( (HLT_HIL1MinimumBiasHF2AND_part1_v1==1||HLT_HIL1MinimumBiasHF2AND_part2_v1==1||HLT_HIL1MinimumBiasHF2AND_part3_v1==1||HLT_HIL1MinimumBiasHF2AND_part4_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part6_v1==1||HLT_HIL1MinimumBiasHF2AND_part7_v1==1||HLT_HIL1MinimumBiasHF2AND_part8_v1==1||HLT_HIL1MinimumBiasHF2AND_part9_v1==1)&&(pprimaryVertexFilter==1&&phfCoincFilter3==1) )",
									"( (HLT_HIL1MinimumBiasHF2AND_part1_v1==1||HLT_HIL1MinimumBiasHF2AND_part2_v1==1||HLT_HIL1MinimumBiasHF2AND_part3_v1==1||HLT_HIL1MinimumBiasHF2AND_part4_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part6_v1==1||HLT_HIL1MinimumBiasHF2AND_part7_v1==1||HLT_HIL1MinimumBiasHF2AND_part8_v1==1||HLT_HIL1MinimumBiasHF2AND_part9_v1==1)&&(pclusterCompatibilityFilter==1&&phfCoincFilter3==1) )",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1&&phfCoincFilter3==1)",
									"( (HLT_HIL1MinimumBiasHF2AND_part1_v1==1||HLT_HIL1MinimumBiasHF2AND_part2_v1==1||HLT_HIL1MinimumBiasHF2AND_part3_v1==1||HLT_HIL1MinimumBiasHF2AND_part4_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part5_v1==1||HLT_HIL1MinimumBiasHF2AND_part6_v1==1||HLT_HIL1MinimumBiasHF2AND_part7_v1==1||HLT_HIL1MinimumBiasHF2AND_part8_v1==1||HLT_HIL1MinimumBiasHF2AND_part9_v1==1)&&(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1&&phfCoincFilter3==1) )"};
//}}}
*/
/*
//Name of selection 2015 Hydjet{{{
TString SelName[NSel] = {"", "MB", "PV", "CC", "Coin", "MB_PV", "MB_CC",
								"MB_Coin", "PV_CC", "PV_Coin", "CC_Coin",
								"MB_PV_CC", "MB_PV_Coin", "MB_CC_Coin", 
								"PV_CC_Coin", "MB_PV_CC_Coin"};
//}}}

//Selection cut 2015 Hydjet{{{
TString Selection[NSel] = {"1",
									"(HLT_HIL1MinimumBiasHF2AND_v1==1)",
									"(pprimaryVertexFilter==1)",
									//"pBeamScrapingFilter==1",
									"(pclusterCompatibilityFilter==1)",
									"(phfCoincFilter3==1)",
									"( (HLT_HIL1MinimumBiasHF2AND_v1==1)&&(pprimaryVertexFilter==1) )",
									"( (HLT_HIL1MinimumBiasHF2AND_v1==1)&&(pclusterCompatibilityFilter==1) )",
									"( (HLT_HIL1MinimumBiasHF2AND_v1==1)&&(phfCoincFilter3==1) )",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1)",
									"(pprimaryVertexFilter==1&&phfCoincFilter3==1)",
									"(pclusterCompatibilityFilter==1&&phfCoincFilter3==1)",
									"( (HLT_HIL1MinimumBiasHF2AND_v1==1)&&(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1) )",
									"( (HLT_HIL1MinimumBiasHF2AND_v1==1)&&(pprimaryVertexFilter==1&&phfCoincFilter3==1) )",
									"( (HLT_HIL1MinimumBiasHF2AND_v1==1)&&(pclusterCompatibilityFilter==1&&phfCoincFilter3==1) )",
									"(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1&&phfCoincFilter3==1)",
									"( (HLT_HIL1MinimumBiasHF2AND_v1==1)&&(pprimaryVertexFilter==1&&pclusterCompatibilityFilter==1&&phfCoincFilter3==1) )"};
//}}}
*/
#endif