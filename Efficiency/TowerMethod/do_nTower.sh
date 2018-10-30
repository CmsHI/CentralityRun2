#!/bin/bash
#make_nTower.C(bool PV, bool BS, const Int_t GF, const Int_t Coin, Int_t MCver)
#PV = pPAprimaryVertexFilter
#BS = pBeamScrapingFilter
#GF = 0: no filter, 1: pVertexFilterCutG, 2: pVertexFilterCutGplus
#Cont = phfCoincFilter level
#MCver = 0: EPOS, 1: HIJING

#root -l -b -q 'make_nTower.C+(true, true, 0, 1, 0)'

root -l -b -q 'draw_nTower.C+(true, true, 0, 1, 0)'
