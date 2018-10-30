//Headers{{{
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "stdio.h"
#include "string.h"
#include <iostream>
#include "TEfficiency.h"
#include "../../analysisUtils.h"
#include "../../utils.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"
#include "TParameter.h"
#include "TObjArray.h"
#include "bin.h"
//}}}

void foldEfftoData(bool isPbPb = false, TString wvar1 = "hiHF", TString wvar2 = "hiHF", TString wvar3 = "hiHF", TString var = "hiNtracks", bool Tr = true, bool PV = true, bool BS = true, const Int_t GF = 0, const Int_t Coin = 1, const Int_t iter = 0, const Int_t MCver = 0)
{
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0000);

//set names{{{
	TString Pname;
	if(isPbPb) Pname = "PbPb";
	else Pname = "pPb";
	TString MCname;
	if (isPbPb)
	{
		if(MCver==0) MCname = "HYDJET";
		else if (MCver==1) MCname = "Cymbal5Ev8_new";
		else if (MCver==2) MCname = "Cymbal5Ev8";
		else if (MCver==3) MCname = "Cymbal5Ev7";
		else if (MCver==4) MCname = "Cymbal5Ev6";
		else if (MCver==5) MCname = "Cymbal5Ev5";
		else if (MCver==6) MCname = "Cymbal5Ev4";
		else if (MCver==7) MCname = "Cymbal5Ev3";
		else if (MCver==8) MCname = "Cymbal5Ev2";
		else if (MCver==9) MCname = "Cymbal5Ev1";
	}
	else
	{
		if (MCver==0) MCname = "EPOS";
		else if (MCver==1) MCname = "HIJING";
		else if (MCver==2) MCname = "AMPT";
	}
//}}}

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

//get tree{{{
	TChain* tref = new TChain("hiEvtAnalyzer/HiTree","");
	TChain* tskimanalysis = new TChain("skimanalysis/HltTree","");
	TChain* thltanalysis = new TChain("hltanalysis/HltTree","");
	if(isPbPb)
	{
	    const char* fileData = Form("root://cms-xrd-global.cern.ch//eos/cms/store/group/phys_heavyions/jmartinb/HiForestAODPbPb5TeV/HIMinimumBias2/EventTree_PbPb_data_HIMinimumBias2_run_15Feb.root");
	    tref->Add(fileData);
	    tskimanalysis->Add(fileData);
	    thltanalysis->Add(fileData);
	}
	else
	{
		const char* fileData = Form("root://eoscms//eos/cms/store/group/phys_heavyions/jmartinb/pPb2016/DATA/HiForestAOD_MinimumBiasDATA_pPb8TeV_run286471.root");
		//for(Int_t i = 0; i < 5; i++)
		//{
	   // tref->Add(Form("%s%d.root", fileData, i));
	   // tskimanalysis->Add(Form("%s%d.root", fileData, i));
	   // thltanalysis->Add(Form("%s%d.root", fileData, i));
		//}
	    tref->Add(Form("%s", fileData));
	    tskimanalysis->Add(Form("%s", fileData));
	    thltanalysis->Add(Form("%s", fileData));
	}
	tref->AddFriend(tskimanalysis);
	tref->AddFriend(thltanalysis);
//}}}

	Long64_t nEntries = tref->GetEntries();
	nEntries = 500000;
	//nEntries = 1000;

//get variables{{{
	Float_t hiHF,hiHFhit,hiEB,hiEE,hiET;
	int hiNtracks,hiNpix;
	Float_t vz;
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
	tref->SetBranchAddress("hiHF",&hiHF);
	tref->SetBranchAddress("hiHFhit",&hiHFhit);
	tref->SetBranchAddress("hiEB",&hiEB);
	tref->SetBranchAddress("hiEE",&hiEE);
	tref->SetBranchAddress("hiET",&hiET);
	tref->SetBranchAddress("hiNtracks",&hiNtracks);
	tref->SetBranchAddress("hiNpix",&hiNpix);
	tref->SetBranchAddress("vz",&vz);
	if(isPbPb)
	{
	  tref->SetBranchAddress("HLT_HIL1MinimumBiasHF1AND_v1",&HLT_HIL1MinimumBiasHF1AND_v1);
	  tref->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
	}
	else
	{
	  tref->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2);
	  tref->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2);
	  tref->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2);
	  tref->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2);
	  tref->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2);
	  tref->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2);
	  tref->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2);
	  tref->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2);
	  tref->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
	  tref->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
	  //tref->SetBranchAddress("pVertexFilterCutG",&pVertexFilterCutG);
	  tref->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
	  tref->SetBranchAddress("phfCoincFilter1",&phfCoincFilter1);
	  tref->SetBranchAddress("phfCoincFilter2",&phfCoincFilter2);
	  tref->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3);
	  tref->SetBranchAddress("phfCoincFilter4",&phfCoincFilter4);
	}
//}}}

	Double_t BinsArr[nBins+1] = {0.};
	for(Int_t i = 0; i < nBins+1; i++)
	{
		if(i == 0) BinsArr[i] = xMin;
		else if(i< nBinsCut1) BinsArr[i] = BinsArr[i-1]+Binwidth1;
		else if(i< nBinsCut2) BinsArr[i] = BinsArr[i-1]+Binwidth2;
		else BinsArr[i] = BinsArr[i-1]+Binwidth3;
	}

//get histogram{{{
	const char* hName = Form("h_%s",var.Data());
	TH1D* h = new TH1D(hName,Form("%s;%s;entries",hName,var.Data()),nBins,BinsArr);
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	TH1D* hW = (TH1D*)(h->Clone(Form("%s_weighted",hName)));
	TString Finname;
	if(iter == 0) Finname = Form("MC_eff_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data());
	else if(iter == 1) Finname = Form("MC_eff_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data());
	else if(iter == 2) Finname = Form("MC_eff_%s_by_%s_and_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data());
	TFile* fEff = TFile::Open(Form("MCefficiency/%s", Finname.Data()));
	if (!fEff)
	{
	  cout << "[ERROR] No efficiency file found" << endl;
	  return;
	}
	TEfficiency* eff = static_cast<TEfficiency*>(fEff->FindObjectAny(Form("h_eff_%s",var.Data())));
	if (!eff)
	{
		cout << "[ERROR] No Tefficiency found for " << var.Data() << endl;
		return;
	}
	TH1D* htot = static_cast<TH1D*>(fEff->FindObjectAny(Form("htot_%s",var.Data())));
	if(!htot)
	{
		cout << "[ERROR] No TH1D found for " << var.Data() << endl;
		return;
	}
//}}}

	for (Long64_t jentry=0; jentry<nEntries;jentry++)
	{
		tref->GetEntry(jentry);
    
		if( jentry % 10000 == 0 ) std::cout << jentry << "/" << nEntries << std::endl;

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
		}
//}}}

//Fill hist{{{

//set filling variable{{{
		double currVar(0.);
		Int_t cutoff = 0;
		if(var == "hiET")
		{
			currVar = hiET;
			cutoff = ETCutoff;
		}
		if(var == "hiEB")
		{
			currVar = hiEB;
			cutoff = EBCutoff;
		}
		if(var == "hiEE")
		{
			currVar = hiEE;
			cutoff = EECutoff;
		}
		if(var == "hiHF")
		{
			currVar = hiHF;
			cutoff = HFCutoff;
		}
		if(var == "hiHFhit")
		{
			currVar = hiHFhit;
			cutoff = HFhitCutoff;
		}
		if(var == "hiNpix")
		{
			currVar = hiNpix;
			cutoff = NpixCutoff;
		}
		if(var == "hiNtracks")
		{
			currVar = hiNtracks;
			cutoff = NTrkCutoff;
		}
		if(var == "vz")
		{
			currVar = vz;
			cutoff = vzCutoff;
		}
//}}}

		h->Fill(currVar);

		int bin = eff->FindFixBin(currVar);
		double effVal = eff->GetEfficiency(bin);
		Double_t correction = 0;
		if(bin > cutoff || effVal <= 0) correction = 1.;
		//if(effVal <= 0) correction = 1.;
		else correction = effVal;
		hW->Fill(currVar, 1./correction);
//}}}
	}

	for(Int_t i = 0; i < nBins+1; i++)
	{
		h->SetBinContent(i+1, h->GetBinContent(i+1)/h->GetBinWidth(i+1));
		hW->SetBinContent(i+1, hW->GetBinContent(i+1)/hW->GetBinWidth(i+1));
	}

	TString Foutname;
	if(iter == 0) Foutname = Form("eff_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data());
	else if(iter == 1) Foutname = Form("eff_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data());
	else if(iter == 2) Foutname = Form("eff_%s_by_%s_and_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data());
	TFile* fsave = new TFile(Form("DATAefficiency/%s", Foutname.Data()),"RECREATE");

	double effVal = ((double)h->Integral("width")/(double)hW->Integral("width"));
	cout << "The efficiency from " << var << " is " << effVal << endl;

//Draw corrected data{{{
	TCanvas *c_eff = new TCanvas(Form("c_eff_%s",var.Data()), "c_eff", 600, 600);
	c_eff->cd();
	c_eff->SetLogy();
	SetHistColor(h,2);
	SetHistColor(hW,1);
	h->SetMarkerStyle(20);
	h->SetMarkerSize(0.8);
	hW->SetMarkerStyle(21);
	hW->SetMarkerSize(0.8);

	TLegend* l1 = new TLegend(0.5,0.6,0.90,0.8);
	legStyle(l1);
	l1->AddEntry(h, "Data", "p");
	l1->AddEntry(hW, "Data corrected", "p");

	hW->Draw();
	h->Draw("same");
	l1->Draw("same");

	double offSet(0);
	if(var == "hiEB") offSet = 0.2;
	drawText(Form("eff. = %0.3f",effVal),0.2,0.2+0.3, 1, 30);
	drawText(Form("%s eff. weights", MCname.Data()),0.2+offSet,0.2+0.1, 1, 30);
// drawText(Form("%s",trigcap.Data()),0.2+offSet,0.2+0.05);
// drawText(Form("%s",evSelCutCap.Data()),0.2+offSet,0.2);

	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	string saveDIR = mainDIR + "/DATAefficiency/";
	void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
	if (dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.c_str(), kTRUE);
	if(iter == 0) c_eff->SaveAs(Form("DATAefficiency/Data_eff_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
	else if(iter == 1) c_eff->SaveAs(Form("DATAefficiency/Data_eff_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
	else if(iter == 2) c_eff->SaveAs(Form("DATAefficiency/Data_eff_%s_by_%s_and_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
//}}}

	TCanvas *c_weight = new TCanvas(Form("c_weight_%s",var.Data()), "c_weight", 600, 600);
	c_weight->cd();
	c_weight->SetLogy();
	TH1D* hratio = (TH1D*) hW->Clone(Form("%s_ratio", hName));
	hratio->Divide(htot);
	hratio->Draw();
	if(iter == 0) c_weight->SaveAs(Form("DATAefficiency/weight_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
	else if(iter == 1) c_weight->SaveAs(Form("DATAefficiency/weight_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
	else if(iter == 2) c_weight->SaveAs(Form("DATAefficiency/weight_%s_by_%s_and_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));

	fsave->cd();
	hW->Write();
	h->Write();
	hratio->Write();
	fsave->Close();
}
