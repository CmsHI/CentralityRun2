#!/bin/bash
#do_fit.C(TString var, bool PV, bool BS, const Int_t GF, const Int_t Coin, Int_t MCver)
#PV = pPAprimaryVertexFilter
#BS = pBeamScrapingFilter
#GF = 0: no filter, 1: pVertexFilterCutG, 2: pVertexFilterCutGplus
#Cont = phfCoincFilter level
#MCver = 0: EPOS, 1: HIJING

Var1="hiNtracks"
Var2="hiHFplusEta4"
Var3="hiHFplus"
Var4="hiHFminusEta4"
Var5="hiHFminus"
Var6="hiHF"
Var7="hiHFhit"

root -l -b -q 'fit.C+('\"$Var1\"', true, true, 0, 1, 0)'
