//Headers{{{
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TLegend.h>
#include <TString.h>
#include <TCut.h>
#include <TMath.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <TEfficiency.h>
#include "../../analysisUtils.h"
#include "../../utils.h"
#include <TSystem.h>
#include <TParameter.h>
#include "bin.h"
//}}}

void MCweight(bool isPbPb = false, TString wvar1 = "vz", TString wvar2 = "vz", TString wvar3 = "vz", TString var = "hiHF", bool Tr = true, bool PV = true, bool BS = true, const Int_t GF = 0, const Int_t Coin = 1, const Int_t iter = 1, const Int_t MCver = 0)
{
	gStyle->SetOptStat(0000);

//get MC file{{{
	TFile* infMC = 0;
	TString MCname;
  if (isPbPb)
	{
		if(MCver==0)
		{
			infMC = TFile::Open("root://cmsxrootd.fnal.gov///store/user/rbi/merged/Hydjet_Quenched_MinBias_5020GeV_750-HINPbPbWinter16DR-NoPU_75X_mcRun2_HeavyIon_forest_v2/0.root");
			MCname = "HYDJET";
		}
	  else if (MCver==1)
		{
			infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/mverweij/MC/PbPb/Hydjet_Quenched_Cymbal5Ev8_PbPbMinBias_5020GeV/crab_HiForestCymbal5Ev8offGENSIM/170113_161526/merge/HiForest_0.root");
			MCname = "Cymbal5Ev8_new";
		}
	  else if (MCver==2)
		{
			infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/mverweij/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev8/crab_HiForestCymbal5Ev8/160912_154217/merge/HiForest_0.root");
			MCname = "Cymbal5Ev8";
		}
	  else if (MCver==3)
		{
			infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/mverweij/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev7/crab_HiForestCymbal5Ev7/160912_154053/merge/HiForest_0.root");
			MCname = "Cymbal5Ev7";
		}
	  else if (MCver==4)
		{
			infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/mverweij/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev6/crab_HiForestCymbal5Ev6/160912_153731/merge/HiForest_0.root");
			MCname = "Cymbal5Ev6";
		}
	  else if (MCver==5)
		{
			infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/cmst3/user/mverweij/jetsPbPb/Run2/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev5/crab_HiForestCymbal5Ev5/merge/HiForest_0.root");
			MCname = "Cymbal5Ev5";
		}
	  else if (MCver==6)
		{
			infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/cmst3/user/mverweij/jetsPbPb/Run2/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev4/crab_HiForestCymbal5Ev4/merge/HiForest_0.root");
			MCname = "Cymbal5Ev4";
		}
	  else if (MCver==7)
		{
			infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/mverweij/MC/PbPb/FixedHydjet_Quenched_MinBias_Cymbal5Ev3/crab_HiForestCymbal5Ev3Fixed/160912_093619/merge/HiForest_0.root");
			MCname = "Cymbal5Ev3";
		}
	  else if (MCver==8)
		{
			infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/cmst3/group/hintt/mverweij/CS/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev2/merge/HiForest_0.root");
			MCname = "Cymbal5Ev2";
		}
	  else if (MCver==9)
		{
			infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/cmst3/group/hintt/mverweij/CS/MC/PbPb/Hydjet_Cymbal5Ev1_MinBias_5020GeV/merge/HiForest_0.root");
			MCname = "Cymbal5Ev1";
		}
	}
	else
	{
		if (MCver==0)
		{
			infMC = TFile::Open("root://eoscms//eos/cms/store/group/phys_heavyions/jmartinb/pPb2016/MC/EPOS8TeV/HiForestAOD_EPOSLHC_pPb_8TeV.root");
			MCname = "EPOS";
		}
	  else if (MCver==1)
		{
			infMC = TFile::Open("root://eoscms//eos/cms/store/group/phys_heavyions/jmartinb/pPb2016/MC/HIJING8TeV/HiForestAOD_HIJING_pPb_8TeV.root");
			MCname = "HIJING";
		}
		else if (MCver==2)
		{
			infMC = TFile::Open("root://eoscms//eos/cms/store/group/phys_heavyions/jmartinb/pPb2016/MC/AMPT8TeV/HiForestAOD_AMPT_pPb_8TeV.root");
			MCname = "AMPT";
		}
	}
	TTree *tin = (TTree*) infMC->Get("hiEvtAnalyzer/HiTree");
	tin->AddFriend("hltanalysis/HltTree");
	tin->AddFriend("skimanalysis/HltTree");
//}}}  

//get variables{{{
	Float_t hiHF,hiHFhit,hiEB,hiEE,hiET;
	Float_t vz;
	int hiNtracks,hiNpix;
	Int_t HLT_HIL1MinimumBiasHF1AND_v1,pcollisionEventSelection;
	Int_t HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2;
	Int_t HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2;
	Int_t HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2;
	Int_t HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2;
	Int_t HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2;
	Int_t HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2;
	Int_t HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2;
	Int_t HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2;
	Int_t pPAprimaryVertexFilter;
	Int_t pBeamScrapingFilter;
	//Int_t pVertexFilterCutG;
	Int_t pVertexFilterCutGplus;
	Int_t phfCoincFilter1;
	Int_t phfCoincFilter2;
	Int_t phfCoincFilter3;
	Int_t phfCoincFilter4;
	Int_t ProcessID;
	tin->SetBranchAddress("hiHF",&hiHF);
	tin->SetBranchAddress("hiHFhit",&hiHFhit);
	tin->SetBranchAddress("hiEB",&hiEB);
	tin->SetBranchAddress("hiEE",&hiEE);
	tin->SetBranchAddress("hiET",&hiET);
	tin->SetBranchAddress("hiNtracks",&hiNtracks);
	tin->SetBranchAddress("hiNpix",&hiNpix);
	tin->SetBranchAddress("vz",&vz);
	if(isPbPb)
	{
	  tin->SetBranchAddress("HLT_HIL1MinimumBiasHF1AND_v1",&HLT_HIL1MinimumBiasHF1AND_v1);
	  tin->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
	}
	else
	{
	  tin->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2);
	  tin->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2);
	  tin->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2);
	  tin->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2);
	  tin->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2);
	  tin->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2);
	  tin->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2);
	  tin->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2);
	  tin->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
	  tin->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
	  //tin->SetBranchAddress("pVertexFilterCutG",&pVertexFilterCutG);
	  tin->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
	  tin->SetBranchAddress("phfCoincFilter1",&phfCoincFilter1);
	  tin->SetBranchAddress("phfCoincFilter2",&phfCoincFilter2);
	  tin->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3);
	  tin->SetBranchAddress("phfCoincFilter4",&phfCoincFilter4);
	  tin->SetBranchAddress("ProcessID",&ProcessID);
	}
//}}}

//set color{{{
	Int_t color_tot = 0;
	Int_t color_pass = 0;
	if(MCname == "EPOS")
	{
		color_tot = 94;
		color_pass = 2;
	}
	else if(MCname == "HIJING")
	{
		color_tot = 64;
		color_pass = 4;
	}
	else if(MCname == "AMPT")
	{
		color_tot = 827;
		color_pass = 843;
	}
//}}}

	TString Pname;
	if(isPbPb) Pname = "PbPb";
	else Pname = "pPb";

//get weighting{{{
	TFile* fweight[iter];
	TH1D* hweight[iter];
	for(Int_t iiter = 0; iiter < iter; iiter++)
	{
		if(iiter == 0) fweight[iiter] = new TFile(Form("DATAefficiency/eff_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iiter, MCname.Data()), "READ");
		else if(iiter == 1) fweight[iiter] = new TFile(Form("DATAefficiency/eff_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", wvar2.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iiter, MCname.Data()), "READ");
		else if(iiter == 2) fweight[iiter] = new TFile(Form("DATAefficiency/eff_%s_by_%s_and_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", wvar3.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iiter, MCname.Data()), "READ");
		if(iiter == 0) hweight[iiter] = (TH1D*) fweight[iiter]->Get(Form("h_%s_ratio", wvar1.Data()));
		else if(iiter == 1) hweight[iiter] = (TH1D*) fweight[iiter]->Get(Form("h_%s_ratio", wvar2.Data()));
		else if(iiter == 2) hweight[iiter] = (TH1D*) fweight[iiter]->Get(Form("h_%s_ratio", wvar3.Data()));
	}
//}}}

//Define output file{{{
	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	string saveDIR = mainDIR + "/MCefficiency/";
	void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
	if (dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.c_str(), kTRUE);
	TFile* fout = 0;
	if(iter == 1) fout = new TFile(Form("MCefficiency/MC_eff_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()), "RECREATE");
	else if(iter == 2) fout = new TFile(Form("MCefficiency/MC_eff_%s_by_%s_and_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()), "RECREATE");
	else if(iter == 3) fout = new TFile(Form("MCefficiency/MC_eff_%s_by_%s_and_%s_by_%s_and_%s_by_%s_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), wvar3.Data(), wvar1.Data(), wvar2.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()), "RECREATE");
//}}}


//Ntrack{{{
	const Int_t nBins = 54;
	Double_t xMin = 0;
	Int_t nBinsCut1 = 21;
	Int_t nBinsCut2 = 46;
	Double_t Binwidth1 = 1.0;
	Double_t Binwidth2 = 2.0;
	Double_t Binwidth3 = 20.0;
//}}}

/*
//HF{{{
	const Int_t nBins = 75;
	Double_t xMin = 0;
	Int_t nBinsCut1 = 31;
	Int_t nBinsCut2 = 66;
	Double_t Binwidth1 = 1.0;
	Double_t Binwidth2 = 2.0;
	Double_t Binwidth3 = 20.0;
//}}}
*/
/*
//HFhit{{{
	const Int_t nBins = 565;
	Double_t xMin = -100;
	Int_t nBinsCut1 = 551;
	Int_t nBinsCut2 = 561;
	Double_t Binwidth1 = 2.0;
	Double_t Binwidth2 = 200.0;
	Double_t Binwidth3 = 1000.0;
//}}}
*/
/*
//ET{{{
	const Int_t nBins = 120;
	Double_t xMin = 0;
	Int_t nBinsCut1 = 101;
	Int_t nBinsCut2 = 116;
	Double_t Binwidth1 = 0.1;
	Double_t Binwidth2 = 1.0;
	Double_t Binwidth3 = 5.0;
//}}}
*/
/*
//EB{{{
	const Int_t nBins = 29;
	Double_t xMin = 0;
	Int_t nBinsCut1 = 21;
	Int_t nBinsCut2 = 50;
	Double_t Binwidth1 = 1.0;
	Double_t Binwidth2 = 20.0;
	Double_t Binwidth3 = 5.0;
//}}}
*/
/*
//EE{{{
	const Int_t nBins = 120;
	Double_t xMin = 0;
	Int_t nBinsCut1 = 101;
	Int_t nBinsCut2 = 116;
	Double_t Binwidth1 = 0.1;
	Double_t Binwidth2 = 1.0;
	Double_t Binwidth3 = 5.0;
//}}}
*/
/*
//Npix{{{
	const Int_t nBins = 207;
	Double_t xMin = 0;
	Int_t nBinsCut1 = 101;
	Int_t nBinsCut2 = 201;
	Double_t Binwidth1 = 1.0;
	Double_t Binwidth2 = 5.0;
	Double_t Binwidth3 = 200.0;
//}}}
*/
/*
//vz bin{{{
	const Int_t nBins = 106;
	Double_t xMin = -25;
	Int_t nBinsCut1 = 4;
	Int_t nBinsCut2 = 104;
	Double_t Binwidth1 = 5.0;
	Double_t Binwidth2 = 0.2;
	Double_t Binwidth3 = 5.0;
//}}}
*/

	Double_t BinsArr[nBins+1] = {0.};
	for(Int_t i = 0; i < nBins+1; i++)
	{
		if(i == 0) BinsArr[i] = xMin;
		else if(i< nBinsCut1) BinsArr[i] = BinsArr[i-1]+Binwidth1;
		else if(i< nBinsCut2) BinsArr[i] = BinsArr[i-1]+Binwidth2;
		else BinsArr[i] = BinsArr[i-1]+Binwidth3;
	}
	TH1D* htot = new TH1D(Form("htot_%s", var.Data()), Form(";%s;", var.Data()), nBins, BinsArr);
	htot->Sumw2();
	TH1D* hpass = (TH1D*) htot->Clone(Form("hpass_%s", var.Data()));
	TH1D* dummyH = (TH1D*) htot->Clone("hDummy");
  
//Fill histo{{{
	Long64_t nEntries = tin->GetEntries();
	for(Long64_t ievt = 0; ievt < nEntries; ievt++)
	{
		tin->GetEntry(ievt);
		if(ievt % 10000 == 0) cout << ievt << " / " << nEntries << endl;

//set variable name{{{
		Double_t wVar1 = 0.;
		Double_t wVar2 = 0.;
		Double_t wVar3 = 0.;
		Double_t xVar = 0.;
		if(wvar1 == "hiET") wVar1 = hiET;
		if(wvar1 == "hiEB") wVar1 = hiEB;
		if(wvar1 == "hiEE") wVar1 = hiEE;
		if(wvar1 == "hiHF") wVar1 = hiHF;
		if(wvar1 == "hiHFhit") wVar1 = hiHFhit;
		if(wvar1 == "hiNpix") wVar1 = hiNpix;
		if(wvar1 == "hiNtracks") wVar1 = hiNtracks;
		if(wvar1 == "vz") wVar1 = vz;
		if(wvar2 == "hiET") wVar2 = hiET;
		if(wvar2 == "hiEB") wVar2 = hiEB;
		if(wvar2 == "hiEE") wVar2 = hiEE;
		if(wvar2 == "hiHF") wVar2 = hiHF;
		if(wvar2 == "hiHFhit") wVar2 = hiHFhit;
		if(wvar2 == "hiNpix") wVar2 = hiNpix;
		if(wvar2 == "hiNtracks") wVar2 = hiNtracks;
		if(wvar3 == "hiET") wVar3 = hiET;
		if(wvar3 == "hiEB") wVar3 = hiEB;
		if(wvar3 == "hiEE") wVar3 = hiEE;
		if(wvar3 == "hiHF") wVar3 = hiHF;
		if(wvar3 == "hiHFhit") wVar3 = hiHFhit;
		if(wvar3 == "hiNpix") wVar3 = hiNpix;
		if(wvar3 == "hiNtracks") wVar3 = hiNtracks;
		if(var == "hiET") xVar = hiET;
		if(var == "hiEB") xVar = hiEB;
		if(var == "hiEE") xVar = hiEE;
		if(var == "hiHF") xVar = hiHF;
		if(var == "hiHFhit") xVar = hiHFhit;
		if(var == "hiNpix") xVar = hiNpix;
		if(var == "hiNtracks") xVar = hiNtracks;
//}}}

		Double_t wval = 1.;
		for(Int_t iiter = 0; iiter < iter; iiter++)
		{
			Int_t wbin = 0;
			if(iiter == 0) wbin = hweight[iiter]->FindFixBin(wVar1);
			else if(iiter == 1) wbin = hweight[iiter]->FindFixBin(wVar2);
			else if(iiter == 2) wbin = hweight[iiter]->FindFixBin(wVar3);
			Double_t wvaltmp = hweight[iiter]->GetBinContent(wbin);
			if(wvaltmp <= 0) wvaltmp = 1.;
			wval = wval*wvaltmp;
		}
		htot->Fill(xVar, wval);

//apply cut{{{
		if(isPbPb)
		{
			if (!HLT_HIL1MinimumBiasHF1AND_v1 || !pcollisionEventSelection) continue;
		}
		else
		{
			if(!(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2)) continue;
			if(PV && !pPAprimaryVertexFilter) continue;
			if(BS && !pBeamScrapingFilter) continue;
			//if(GF == 1 && !pVertexFilterCutG) continue;
			if(GF == 2 && !pVertexFilterCutGplus) continue;
			if(Coin == 1 && !phfCoincFilter1) continue;
			if(Coin == 2 && !phfCoincFilter2) continue;
			if(Coin == 3 && !phfCoincFilter3) continue;
			if(Coin == 4 && !phfCoincFilter4) continue;
			if(MCname == "EPOS")
			{
				 if( !(ProcessID==101||ProcessID==105) ) continue;
			}
		}
//}}}

		hpass->Fill(xVar, wval);
	}
//}}}

	for(Int_t i = 0; i < nBins; i++)
	{
		htot->SetBinContent(i+1, htot->GetBinContent(i+1)/htot->GetBinWidth(i+1));
		hpass->SetBinContent(i+1, hpass->GetBinContent(i+1)/hpass->GetBinWidth(i+1));
	}

	TParameter<double>* effVal = new TParameter<double>(Form("v_eff_%s", var.Data()),(double)(hpass->Integral("width")/htot->Integral("width")));

	TEfficiency* effH = new TEfficiency(*hpass,*htot);
	effH->SetName(Form("h_eff_%s", var.Data()));

	TCanvas* c_dist = new TCanvas(Form("c_dist_%s", var.Data()), "", 0, 0, 600, 600);
	htot->SetMarkerStyle(20);
	htot->SetMarkerColor(color_tot);
	htot->SetLineColor(color_tot);
	hpass->SetMarkerStyle(20);
	hpass->SetMarkerColor(color_pass);
	hpass->SetLineColor(color_pass);
	htot->Draw("pe");
	hpass->Draw("samepe");
	TLegend* l1 = new TLegend(0.7, 0.7, 0.9, 0.9);
	l1->AddEntry(htot, "total", "pe");
	l1->AddEntry(hpass, "selected", "pe");
	l1->Draw();
	if(iter == 1) c_dist->SaveAs(Form("MCefficiency/MCdist_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
	else if(iter == 2) c_dist->SaveAs(Form("MCefficiency/MCdist_%s_by_%s_and_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
	else if(iter == 3) c_dist->SaveAs(Form("MCefficiency/MCdist_%s_by_%s_and_%s_by_%s_and_%s_by_%s_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), wvar3.Data(), wvar1.Data(), wvar2.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));

	TH1D* hratio = (TH1D*) hpass->Clone(Form("hratio_%s",var.Data()));

	TCanvas* c_eff = new TCanvas(Form("c_eff_%s", var.Data()), "", 0, 0, 600, 600);
	TGraphAsymmErrors* geff = effH->CreateGraph();
	geff->SetName(Form("g_eff_%s", var.Data()));
	dummyH->Draw();
	geff->Draw("same");
//	hratio->Divide(htot);
//	hratio->Draw("samepe");
	drawText(Form("eff = %0.3f", effVal->GetVal()), 0.2, 0.2+0.3, 1, 30);
	drawText(MCname.Data(), 0.2, 0.2+0.1, 1, 30);
	TLegend* l2 = new TLegend(0.7, 0.5, 0.9, 0.7);
//	l2->AddEntry(geff, "Teff", "pe");
//	l2->AddEntry(hratio, "ratio", "pe");
//	l2->Draw();
	if(iter == 1) c_eff->SaveAs(Form("MCefficiency/MC_eff_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
	else if(iter == 2) c_eff->SaveAs(Form("MCefficiency/MC_eff_%s_by_%s_and_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
	else if(iter == 3) c_eff->SaveAs(Form("MCefficiency/MC_eff_%s_by_%s_and_%s_by_%s_and_%s_by_%s_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), wvar3.Data(), wvar1.Data(), wvar2.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));

	fout->cd();
	geff->Write();
	effH->Write();
	effVal->Write();
	htot->Write();
	hpass->Write();
	//hratio->Write();
}
