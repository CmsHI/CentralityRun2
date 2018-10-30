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

void MCeff(bool isPbPb = false, TString var = "vz", bool Tr = true, bool PV = true, bool BS = true, const Int_t GF = 0, const Int_t Coin = 1, const Int_t MCver = 0)
{
	gStyle->SetOptStat(0000);

//make directory{{{
	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	string saveDIR = mainDIR + "/MCefficiency/";
	void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
	if (dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.c_str(), kTRUE);
//}}}

	TString Pname;
	if(isPbPb) Pname = "PbPb";
	else Pname = "pPb";
	TCanvas* c_dist = new TCanvas("c_dist", "", 0, 0, 600, 600);
	TCanvas* c_eff = new TCanvas("c_eff", "", 0, 0, 600, 600);

//SetColor{{{
	Int_t color_tot = 0;
	Int_t color_pass = 0;
	if(MCver == 0)
	{
		color_tot = 94;
		color_pass = 2;
	}
	else if(MCver == 1)
	{
		color_tot = 64;
		color_pass = 4;
	}
	else if(MCver == 2)
	{
		color_tot = 827;
		color_pass = 843;
	}
//}}}

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
//}}}  

	TTree* t = (TTree*) infMC->Get("hiEvtAnalyzer/HiTree");
	t->AddFriend("hltanalysis/HltTree");
	t->AddFriend("skimanalysis/HltTree");

/*
//Ntrack{{{
	const Int_t nBins = 54;
	Double_t xMin = 0;
	Int_t nBinsCut1 = 21;
	Int_t nBinsCut2 = 46;
	Double_t Binwidth1 = 1.0;
	Double_t Binwidth2 = 2.0;
	Double_t Binwidth3 = 20.0;
//}}}
*/
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
	const Int_t nBins = 201;
	Double_t xMin = -100;
	Int_t nBinsCut1 = 5;
	Int_t nBinsCut2 = 14;
	Double_t Binwidth1 = 50.0;
	Double_t Binwidth2 = 100.0;
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

//Npix{{{
	const Int_t nBins = 107;
	Double_t xMin = 0;
	Int_t nBinsCut1 = 51;
	Int_t nBinsCut2 = 101;
	Double_t Binwidth1 = 2.0;
	Double_t Binwidth2 = 10.0;
	Double_t Binwidth3 = 200.0;
//}}}

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
	//	cout << "input: " << i << "th: " << BinsArr[i] << endl;
	}

/*
	for(Int_t j = 0; j < nBins+1; j++)
	{
		cout << "output: " << j << "th: " << BinsArr[j] << endl;
	}
*/
	TH1D* htot = new TH1D(Form("htot_%s",var.Data()), Form(";%s;",var.Data()), nBins, BinsArr);
	htot->Sumw2();
	TH1D* hpass = (TH1D*) htot->Clone(Form("hpass_%s",var.Data()));
	TH1D* dummyH = (TH1D*) htot->Clone("hDummy");

//Event selection cut{{{
	TCut PVfilter, BSfilter, Gfilter, Coinfilter, Trigger, PrID;
	if(isPbPb)
	{
		if(PV) PVfilter = "";
		else PVfilter = "";
		if(BS) BSfilter = "";
		else BSfilter = "";
		Gfilter = "";
		if(Coin == 1) Coinfilter = "phfCoincFilter1==1";
		else if(Coin == 2) Coinfilter = "phfCoincFilter2==1";
		else if(Coin == 3) Coinfilter = "phfCoincFilter3==1";
		else if(Coin == 4) Coinfilter = "phfCoincFilter4==1";
		else Coinfilter = "";
		if(Tr) Trigger = "HLT_HIL1MinimumBiasHF1AND_v1==1";
		else Trigger = "";
	}
	else
	{
		if(PV) PVfilter = "pPAprimaryVertexFilter==1";
		else PVfilter = "";
		if(BS) BSfilter = "pBeamScrapingFilter==1";
		else BSfilter = "";
		if(GF == 1) Gfilter = "pVertexFilterCutG==1";
		else if(GF == 2) Gfilter = "pVertexFilterCutGplus==1";
		else Gfilter = "";
		if(Coin == 1) Coinfilter = "phfCoincFilter1==1";
		else if(Coin == 2) Coinfilter = "phfCoincFilter2==1";
		else if(Coin == 3) Coinfilter = "phfCoincFilter3==1";
		else if(Coin == 4) Coinfilter = "phfCoincFilter4==1";
		else Coinfilter = "";
		Trigger = "(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2==1)";
		if(MCname == "EPOS") PrID = "ProcessID==101||ProcessID==105";
		else PrID = "";
	}
//}}}

	t->Draw(Form("%s>>%s", var.Data(), htot->GetName()), PrID);
	t->Draw(Form("%s>>%s", var.Data(), hpass->GetName()), Trigger&&PVfilter&&BSfilter&&Gfilter&&Coinfilter&&PrID);

	for(Int_t i = 0; i < nBins; i++)
	{
		htot->SetBinContent(i+1, htot->GetBinContent(i+1)/htot->GetBinWidth(i+1));
		hpass->SetBinContent(i+1, hpass->GetBinContent(i+1)/hpass->GetBinWidth(i+1));
	}

//Draw dist{{{
	c_dist->cd();
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
	c_dist->SaveAs(Form("MCefficiency/MC_dist_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter0_%s.pdf", var.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, MCname.Data()));
//}}}

//Efficiency{{{
  TParameter<double>* effVal = new TParameter<double>(Form("v_eff_%s", var.Data()),((double)hpass->Integral("width")/(double)htot->Integral("width")));
  TEfficiency* effH = new TEfficiency(*hpass,*htot);
  effH->SetName(Form("h_eff_%s",var.Data()));
//}}}

//Draw efficiency{{{
	Double_t effout = ((double)hpass->Integral()/(double)htot->Integral());
	cout << "efficiency by ratio: " << effout << endl;

	TH1D* hratio = (TH1D*) hpass->Clone(Form("hratio_%s",var.Data()));

	c_eff->cd();
	TGraphAsymmErrors* geff = effH->CreateGraph();
	geff->SetName(Form("g_eff_%s",var.Data()));
	dummyH->Draw();
	geff->Draw("same");
//	hratio->Divide(htot);
//	hratio->Draw("samepe");
	drawText(Form("eff. = %0.3f",effVal->GetVal()),0.2,0.2+0.3, 1, 30);
	drawText(MCname.Data(),0.2,0.2+0.1, 1, 30);
//	TLegend* l2 = new TLegend(0.7, 0.5, 0.9, 0.7);
//	l2->AddEntry(geff, "Teff", "pe");
//	l2->AddEntry(hratio, "ratio", "pe");
//	l2->Draw();
	c_eff->SaveAs(Form("MCefficiency/MC_eff_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter0_%s.pdf",var.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, MCname.Data()));
//}}}

//Save{{{
	TFile* fout = new TFile(Form("MCefficiency/MC_eff_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter0_%s.root", var.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, MCname.Data()), "RECREATE");
	fout->cd();
	geff->Write();
	effH->Write();
	effVal->Write();
	htot->Write();
	hpass->Write();
	hratio->Write();
//}}}
}
