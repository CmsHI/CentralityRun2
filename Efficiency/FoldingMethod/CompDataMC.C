//Headers{{{
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TPad.h>
#include <TLegend.h>
#include <TString.h>
#include <TLine.h>
#include <iostream>
#include "bin.h"
//}}}

void CompDataMC(bool isPbPb = false, TString wvar1 = "hiHF", TString wvar2 = "hiHF", TString wvar3 = "hiHF", TString var = "hiHF", bool Tr = true, bool PV = true, bool BS = true, const Int_t GF = 0, const Int_t Coin = 1, const Int_t iter = 0, const Int_t MCver = 1)
{
	gStyle->SetOptStat(0000);
	TString RDname = "root://eoscms//eos/cms/store/group/phys_heavyions/jmartinb/pPb2016/DATA/HiForestAOD_MinimumBiasDATA_pPb8TeV_run286471.root";

//set names{{{
	TString MCname;
	if(MCver == 0) MCname = "EPOS";
	else if(MCver == 1) MCname = "HIJING";
	else if(MCver == 2) MCname = "AMPT";
	TString Pname;
	if(isPbPb) Pname = "PbPb";
	else Pname = "pPb";
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

//Get tree from file{{{
	TChain *t_RD_evt;
	TChain *t_RD_skim;
	TChain *t_RD_hlt;
	t_RD_evt = new TChain("hiEvtAnalyzer/HiTree");
	t_RD_skim = new TChain("skimanalysis/HltTree");
	t_RD_hlt = new TChain("hltanalysis/HltTree");

/*
	for(Int_t i = 0; i < 5; i++)
	{
		t_RD_evt->Add(Form("%s%d.root", RDname.Data(), i));
		t_RD_skim->Add(Form("%s%d.root", RDname.Data(), i));
		t_RD_hlt->Add(Form("%s%d.root", RDname.Data(), i));
	}
*/
	t_RD_evt->Add(Form("%s", RDname.Data()));
	t_RD_evt->AddFriend(t_RD_skim);
	t_RD_evt->AddFriend(t_RD_hlt);
//}}}

//Cuts{{{
	TCut trigger;
	if(Tr) trigger = "(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1==1)";
	else trigger = "";
	TCut PVfilter;
	if(PV) PVfilter = "pPAprimaryVertexFilter==1";
	else PVfilter = "";
	TCut BSfilter;
	if(BS) BSfilter = "pBeamScrapingFilter==1";
	else BSfilter = "";
	TCut Gfilter;
	if(GF == 1) Gfilter = "pVertexFilterCutG==1";
	else if(GF == 2) Gfilter = "pVertexFilterCutGplus==1";
	else Gfilter = "";
	TCut Coinfilter;
	if(Coin == 1) Coinfilter = "phfCoincFilter1==1";
	else if(Coin == 2) Coinfilter = "phfCoincFilter2==1";
	else if(Coin == 3) Coinfilter = "phfCoincFilter3==1";
	else if(Coin == 4) Coinfilter = "phfCoincFilter4==1";
	else Coinfilter = "";
	//TCut Lumi = "(67<lumi&&lumi<1540)";
	TCut Lumi = "";
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
	const Int_t nBins = 21;
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
	}

	TH1D* hRD_pass = new TH1D("hRD_pass", Form(";%s;", var.Data()), nBins, BinsArr);
	hRD_pass->Sumw2();
	t_RD_evt->Draw(Form("%s>>hRD_pass", var.Data()), "", trigger&&PVfilter&&BSfilter&&Gfilter&&Coinfilter&&Lumi);

	for(Int_t i = 0; i < nBins; i++)
	{
		hRD_pass->SetBinContent(i+1, hRD_pass->GetBinContent(i+1)/hRD_pass->GetBinWidth(i+1));
	}

//Get MC hist{{{
	TString Finname;
	if(iter == 0) Finname = Form("MC_eff_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data());
	else if(iter == 1) Finname = Form("MC_eff_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data());
	else if(iter == 2) Finname = Form("MC_eff_%s_by_%s_and_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data());
	else if(iter == 3) Finname = Form("MC_eff_%s_by_%s_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.root", var.Data(), wvar1.Data(), Pname.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data());
	TFile* fin = new TFile(Form("MCefficiency/%s", Finname.Data()), "READ");
	TH1D* hMC_tot = (TH1D*) fin->Get(Form("htot_%s", var.Data()));;
	TH1D* hMC_pass = (TH1D*) fin->Get(Form("hpass_%s", var.Data()));
	TH1D* hMC_pass_1 = (TH1D*) hMC_pass->Clone(Form("hpass_%s_1", var.Data()));
	hMC_tot->Sumw2();
	hMC_pass->Sumw2();
	hMC_pass_1->Sumw2();
//}}}

	TCanvas* c1 = new TCanvas("c1", "", 0, 0, 600, 600);
	c1->Divide(2, 1);
	c1->cd(1)->SetPad(0.01, 0.25, 0.99, 0.99);
	c1->cd(1)->SetBottomMargin(0.0);

	hRD_pass->Scale(1./hRD_pass->Integral("width"));
	hMC_pass->Scale(1./hMC_pass->Integral("width"));
	hMC_pass_1->Scale(1./hMC_pass_1->Integral("width"));
//	hRD_pass->Scale(1./hRD_pass->Integral());
//	hMC_pass->Scale(1./hMC_pass->Integral());
//	hMC_pass_1->Scale(1./hMC_pass_1->Integral());
//	hMC_pass->Scale(1./hMC_tot->Integral("width"));
//	hMC_tot->Scale(1./hMC_tot->Integral("width"));

//Draw total{{{
	hRD_pass->SetMarkerStyle(20);
	hRD_pass->SetMarkerColor(1);
	hRD_pass->SetLineColor(1);
//	hMC_tot->SetMarkerStyle(20);
//	hMC_tot->SetMarkerColor(color_tot);
//	hMC_tot->SetLineColor(color_tot);
	hMC_pass->SetMarkerStyle(20);
	hMC_pass->SetMarkerColor(color_pass);
	hMC_pass->SetLineColor(color_pass);
	hMC_pass_1->SetMarkerStyle(20);
	hMC_pass_1->SetMarkerColor(color_pass);
	hMC_pass_1->SetLineColor(color_pass);
	Double_t max1 = hRD_pass->GetBinContent(hRD_pass->GetMaximumBin());
//	Double_t max2 = hMC_tot->GetBinContent(hMC_tot->GetMaximumBin());
	Double_t max3 = hMC_pass_1->GetBinContent(hMC_pass_1->GetMaximumBin());
//	Double_t maxi = 0;
	Double_t maxf = 0;
//	maxi = max(max1, max2);
//	maxf = max(maxi, max3);
	maxf = max(max1, max3);
	hRD_pass->GetYaxis()->SetRangeUser(0, maxf*1.2);
	hRD_pass->Draw("pe");
//	hMC_tot->Draw("samep");
	//hMC_tot->Draw("");
	hMC_pass_1->Draw("samep");
	TLegend* leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
	leg1->AddEntry(hRD_pass, "Data", "pl");
//	leg1->AddEntry(hMC_tot, Form("%s_total", MCname.Data()), "pl");
	leg1->AddEntry(hMC_pass_1, Form("%s_pass", MCname.Data()), "pl");
	leg1->Draw();

	c1->cd(2)->SetPad(0.01, 0.01, 0.99, 0.25);
	c1->cd(2)->SetTopMargin(0.0);
	c1->cd(2)->SetBottomMargin(0.1);
	c1->cd(2)->SetLogy();
	//TPad* p1 = (TPad*) (c1->cd(2));
	//p1->SetLogy();
	hMC_pass->Divide(hRD_pass);
	hMC_pass->Draw();
	TLine* line1 = new TLine(xMin, 1, BinsArr[nBins], 1);
	line1->SetLineWidth(1);
	line1->SetLineColor(1);
	line1->Draw();
	if(iter == 0) c1->SaveAs(Form("Comp_data_MC/Comp_data_MC_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
	else if(iter == 1) c1->SaveAs(Form("Comp_data_MC/Comp_data_MC_%s_by_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
	else if(iter == 2) c1->SaveAs(Form("Comp_data_MC/Comp_data_MC_%s_by_%s_and_%s_by_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), wvar2.Data(), wvar1.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
	else if(iter == 3) c1->SaveAs(Form("Comp_data_MC/Comp_data_MC_%s_by_%s_Tr%d_PV%d_BS%d_GF%d_Coin%d_iter%d_%s.pdf", var.Data(), wvar1.Data(), Tr, PV, BS, GF, Coin, iter, MCname.Data()));
//}}}
}
