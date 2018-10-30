#!/bin/bash
#CompDist(bool isPbPb, TString var, bool Tr, bool PV, bool BS, const Int_t GF, const Int_t Coin)
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

#root -l -b -q 'CompDist.C+(false, '\"$Var1\"', true, true, true, 2, 1)'
#root -l -b -q 'CompDist.C+(false, '\"$Var2\"', true, true, true, 2, 1)'
#root -l -b -q 'CompDist.C+(false, '\"$Var3\"', true, true, true, 0, 1)'
#root -l -b -q 'CompDist.C+(false, '\"$Var4\"', true, true, true, 2, 1)'
#root -l -b -q 'CompDist.C+(false, '\"$Var5\"', true, true, true, 0, 1)'
#root -l -b -q 'CompDist.C+(false, '\"$Var6\"', true, true, true, 2, 1)'
#root -l -b -q 'CompDist.C+(false, '\"$Var7\"', true, true, true, 2, 1)'

#root -l -b -q 'MCratio.C+(false, '\"$Var1\"', true, true, true, 2, 1)'
#root -l -b -q 'MCratio.C+(false, '\"$Var2\"', true, true, true, 2, 1)'
#root -l -b -q 'MCratio.C+(false, '\"$Var3\"', true, true, true, 2, 1)'
#root -l -b -q 'MCratio.C+(false, '\"$Var4\"', true, true, true, 2, 1)'
#root -l -b -q 'MCratio.C+(false, '\"$Var5\"', true, true, true, 2, 1)'
#root -l -b -q 'MCratio.C+(false, '\"$Var6\"', true, true, true, 2, 1)'
#root -l -b -q 'MCratio.C+(false, '\"$Var7\"', true, true, true, 2, 1)'

#root -l -b -q 'CompDataMC.C+(false, '\"$Var7\"', '\"$Var1\"', true, true, true, 0, 1, 1, 0)'
#root -l -b -q 'CompDataMC.C+(false, '\"$Var7\"', '\"$Var1\"', true, true, true, 0, 1, 1, 1)'
#root -l -b -q 'CompDataMC.C+(false, '\"$Var7\"', '\"$Var1\"', true, true, true, 0, 1, 1, 2)'
root -l -b -q 'CompDataMC.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var7\"', true, true, true, 2, 1, 1, 0)'
root -l -b -q 'CompDataMC.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var7\"', true, true, true, 2, 1, 1, 1)'
root -l -b -q 'CompDataMC.C+(false, '\"$Var2\"', '\"$Var3\"', '\"$Var7\"', '\"$Var7\"', true, true, true, 2, 1, 1, 2)'
