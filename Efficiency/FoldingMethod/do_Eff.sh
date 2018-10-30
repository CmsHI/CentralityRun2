#!/bin/bash
#sel_EffMC(bool isPbPb, bool PV, bool BS, const Int_t GF, const Int_t Coin, Int_t MCver)
#PV = pPAprimaryVertexFilter
#BS = pBeamScrapingFilter
#GF = 0: no filter, 1: pVertexFilterCutG, 2: pVertexFilterCutGplus
#Cont = phfCoincFilter level
#MCver(!isPbPb) = 0: EPOS, 1: HIJING, 2: AMPT
#MCver(isPbPb) = 0: Hydjet, 1: CymbalEv8_new, 2: CymbalEv8, 3: CymbalEv7
#						4: CymbalEv6, 5: CymbalEv5, 6: CymbalEv4, 7: CymbalEv3,
#						8: CymbalEv2, 9: CymbalEv1
Var1="hiNtracks"
Var2="hiHF"
Var3="hiHFhit"
Var4="hiET"
Var5="hiEB"
Var6="hiEE"
Var7="hiNpix"
Var8="vz"

#root -l -b -q 'MCCutCheck.C+('\"$Var8\"', 0)'
#root -l -b -q 'MCCutCheck.C+('\"$Var8\"', 1)'
#root -l -b -q 'MCCutCheck.C+('\"$Var8\"', 2)'

#root -l -b -q 'MCeff.C+(false, '\"$Var7\"', true, true, true, 0, 1, 0)'
#root -l -b -q 'MCeff.C+(false, '\"$Var7\"', true, true, true, 0, 1, 1)'
#root -l -b -q 'MCeff.C+(false, '\"$Var7\"', true, true, true, 0, 1, 2)'
#root -l -b -q 'MCeff.C+(false, '\"$Var7\"', true, true, true, 2, 1, 0)'
#root -l -b -q 'MCeff.C+(false, '\"$Var7\"', true, true, true, 2, 1, 1)'
#root -l -b -q 'MCeff.C+(false, '\"$Var7\"', true, true, true, 2, 1, 2)'

#root -l -b -q 'foldEfftoData.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var3\"', true, true, true, 0, 1, 1, 0)'
#root -l -b -q 'foldEfftoData.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var3\"', true, true, true, 0, 1, 1, 1)'
#root -l -b -q 'foldEfftoData.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var3\"', true, true, true, 0, 1, 1, 2)'
#root -l -b -q 'foldEfftoData.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var7\"', true, true, true, 2, 1, 0, 0)'
#root -l -b -q 'foldEfftoData.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var7\"', true, true, true, 2, 1, 0, 1)'
#root -l -b -q 'foldEfftoData.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var7\"', true, true, true, 2, 1, 0, 2)'

#root -l -b -q 'MCweight.C+(false, '\"$Var7\"', '\"$Var3\"', '\"$Var7\"', '\"$Var2\"', true, true, true, 0, 1, 1, 0)'
#root -l -b -q 'MCweight.C+(false, '\"$Var7\"', '\"$Var3\"', '\"$Var7\"', '\"$Var2\"', true, true, true, 0, 1, 1, 1)'
#root -l -b -q 'MCweight.C+(false, '\"$Var7\"', '\"$Var3\"', '\"$Var7\"', '\"$Var2\"', true, true, true, 0, 1, 1, 2)'
root -l -b -q 'MCweight.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var7\"', true, true, true, 2, 1, 1, 0)'
root -l -b -q 'MCweight.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var7\"', true, true, true, 2, 1, 1, 1)'
root -l -b -q 'MCweight.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var7\"', true, true, true, 2, 1, 1, 2)'

#root -l -b -q 'WeightComp.C+(false, true, true, true, 0, 1, 0)'
#root -l -b -q 'WeightComp.C+(false, true, true, true, 0, 1, 1)'
#root -l -b -q 'WeightComp.C+(false, true, true, true, 0, 1, 2)'
#root -l -b -q 'WeightComp.C+(false, true, true, true, 2, 1, 0)'
#root -l -b -q 'WeightComp.C+(false, true, true, true, 2, 1, 1)'
#root -l -b -q 'WeightComp.C+(false, true, true, true, 2, 1, 2)'

#root -l -b -q 'MCCutCheck.C+(0)'
#root -l -b -q 'MCCutCheck.C+(1)'
#root -l -b -q 'MCCutCheck.C+(2)'
