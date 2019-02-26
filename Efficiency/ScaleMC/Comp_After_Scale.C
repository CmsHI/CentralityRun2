//Headers{{{
#include "../../Utilities/Style_Header.h"
#include "../../Utilities/Var_Header.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TString.h>
#include <TCut.h>
#include <TMath.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <TEfficiency.h>
#include <TSystem.h>
#include <TParameter.h>
#include <TGraphAsymmErrors.h>
//}}}

void Comp_After_Scale(const Int_t ivar = 5, const Int_t Coin = 2, const Int_t Th = 4)
{
	SetStyle();

//make directory{{{
	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	string saveDIR = mainDIR + "/AfterScale/";
	void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
	if (dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.c_str(), kTRUE);
//}}}

	const Int_t NOF = 4;//Number Of Files
	TString akaD = "PR326478";
	TString akaM[NOF] = {"Hydjet_Cymbal5F", "AMPT_String", "AMPT_NoString",
								"EPOS"};
	const Int_t Range[NOF] = {200, 200, 200, 200};

	TFile* fref = new TFile(Form("../../HistFiles/PbPb2018_%s_histo_newCNT.root", akaD.Data()), "READ");
	TFile* fin[NOF];
	TCanvas* c1 = new TCanvas("c1", "", 0, 0, 600, 600);
	TCanvas* c2 = new TCanvas("c2", "", 0, 0, 600, 600);
	TCanvas* c3 = new TCanvas("c3", "", 0, 0, 600, 600);
	TCanvas* c4 = new TCanvas("c4", "", 0, 0, 600, 600);
	TCanvas* c5 = new TCanvas("c5", "", 0, 0, 600, 600);
	TCanvas* c6 = new TCanvas("c6", "", 0, 0, 600, 600);

//get reference hist{{{
	TH1D* href = (TH1D*) fref->Get(Form("h%s_PV_CC_Coin%dth%d", VarName[ivar].Data(), Coin, Th));
	FormTH1Marker(href, 0, 0, 1.4);
	href->Scale(1./( href->Integral(href->GetXaxis()->FindBin(NormRangeCut[ivar]), href->GetXaxis()->FindBin(VarMaxC[ivar])) ));
	href->Scale(1., "width");
//}}}

//get hist{{{
	TH1D* h1[NOF];
	TH1D* h2[NOF];
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		fin[ifile] = new TFile(Form("MCefficiency/MC_eff_2018_%s_Range%d_coin%dth%d_%s_by_%s.root", VarName[ivar].Data(), Range[ifile], Coin, Th, akaM[ifile].Data(), akaD.Data()), "READ");
		h1[ifile] = (TH1D*) fin[ifile]->Get(Form("hsel_%s", VarName[ivar].Data()));
		FormTH1Marker(h1[ifile], ifile+1, ifile+1, 1.4);
		h1[ifile]->Scale(1./( h1[ifile]->Integral(h1[ifile]->GetXaxis()->FindBin(NormRangeCut[ivar]), h1[ifile]->GetXaxis()->FindBin(VarMaxC[ivar])) ));
		h1[ifile]->Scale(1., "width");
		h2[ifile] = (TH1D*) fin[ifile]->Get(Form("heff_%s", VarName[ivar].Data()));
		FormTH1Marker(h2[ifile], ifile+1, ifile+1, 1.4);
	}
//}}}

//make new bin{{{
	Double_t tmpBinArr[nBinC[ivar]];
	tmpBinArr[0] = 0;
	Int_t NewnBin = 0;
	for(Int_t ibin = 1; ibin < nBinC[ivar]+2; ibin++)
	{
		if(tmpBinArr[ibin-1] < NormRangeCut[ivar]) tmpBinArr[ibin] = tmpBinArr[ibin-1]+( VarMaxC[ivar]/(double)nBinC[ivar] );
		else tmpBinArr[ibin] = tmpBinArr[ibin-1]+( VarW[ivar]*VarMaxC[ivar]/(double)nBinC[ivar] );
		if(tmpBinArr[ibin] >= VarMaxC[ivar])
		{
			NewnBin = ibin+1;
			break;
		}
	}
	Double_t NewBinArr[NewnBin];
	NewBinArr[0] = 0;
	for(Int_t ibin = 1; ibin < NewnBin; ibin++)
	{
		if(NewBinArr[ibin-1] < NormRangeCut[ivar]) NewBinArr[ibin] = NewBinArr[ibin-1]+( VarMaxC[ivar]/(double)nBinC[ivar] );
		else NewBinArr[ibin] = NewBinArr[ibin-1]+( VarW[ivar]*VarMaxC[ivar]/(double)nBinC[ivar] );
	}
//}}}

//total dist{{{
	c1->cd();
	c1->SetLogy();
	TLegend* leg1 = new TLegend(0.4, 0.7, 0.9, 0.9);
	FormLegend(leg1, 0.04);
	TH1D* htmp1 = new TH1D("htmp1", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
	htmp1->SetAxisRange(1e-9, 0.8, "Y");
	htmp1->Draw();
	href->Draw("samep");
	leg1->AddEntry(href, Form("%s", akaD.Data()), "p");
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		h1[ifile]->Draw("samep");
		leg1->AddEntry(h1[ifile], Form("%s %d~", akaM[ifile].Data(), Range[ifile]), "p");
	}
	leg1->Draw();
	c1->SaveAs(Form("AfterScale/%s_dist_scaled_by_%s_Range%d_Coin%dth%d.pdf", VarName[ivar].Data(), akaD.Data(), Range[0], Coin, Th));
//}}}

//total dist zoom{{{
	c2->cd();
	c2->SetLogy();
	TLegend* leg2 = new TLegend(0.4, 0.7, 0.9, 0.9);
	FormLegend(leg2, 0.04);
	TH1D* htmp2 = new TH1D("htmp2", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
	htmp2->SetAxisRange(0, NormRangeCut[ivar], "X");
	htmp2->SetAxisRange(1e-9, 0.8, "Y");
	htmp2->Draw();
	href->Draw("samepe");
	leg2->AddEntry(href, Form("%s", akaD.Data()), "p");
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		h1[ifile]->Draw("samepe");
		leg2->AddEntry(h1[ifile], Form("%s %d~", akaM[ifile].Data(), Range[ifile]), "p");
	}
	leg2->Draw();
	c2->SaveAs(Form("AfterScale/%s_dist_scaled_by_%s_Range%d_Coin%dth%d_zoom.pdf", VarName[ivar].Data(), akaD.Data(), Range[0], Coin, Th));
//}}}

//total ratio{{{
	c3->cd();
	TLegend* leg3 = new TLegend(0.4, 0.2, 0.9, 0.45);
	FormLegend(leg3, 0.04);
	TH1D* htmp3 = new TH1D("htmp3", Form("%s;%s;ratio", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
	htmp3->SetAxisRange(0, 1.2, "Y");
	htmp3->Draw();
	href->Rebin(15);
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		h1[ifile]->Rebin(15);
		h1[ifile]->Divide(h1[ifile], href, 1, 1, "b");
/*
		for(Int_t i = 1; i < h1[ifile]->GetNbinsX(); i++)
		{
			h1[ifile]->SetBinError(i, 0.);
		}
*/
		h1[ifile]->Draw("samep");
		leg3->AddEntry(h1[ifile], Form("%s %d~", akaM[ifile].Data(), Range[ifile]), "p");
	}
	leg3->Draw();
	SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
	c3->SaveAs(Form("AfterScale/%s_ratio_scaled_by_%s_Range%d_Coin%dth%d.pdf", VarName[ivar].Data(), akaD.Data(), Range[0], Coin, Th));
//}}}

//total ratio zoom{{{
	c4->cd();
	TLegend* leg4 = new TLegend(0.4, 0.2, 0.9, 0.45);
	FormLegend(leg4, 0.04);
	TH1D* htmp4 = new TH1D("htmp4", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
	htmp4->SetAxisRange(0, NormRangeCut[ivar], "X");
	htmp4->SetAxisRange(0, 1.2, "Y");
	htmp4->Draw();
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		h1[ifile]->Draw("samep");
		leg4->AddEntry(h1[ifile], Form("%s %d~", akaM[ifile].Data(), Range[ifile]), "p");
	}
	leg4->Draw();
	SetLine(1, 0, 1, NormRangeCut[ivar], 1, 0, 2);
	c4->SaveAs(Form("AfterScale/%s_ratio_scaled_by_%s_Range%d_Coin%dth%d_zoom.pdf", VarName[ivar].Data(), akaD.Data(), Range[0], Coin, Th));
//}}}

//total eff{{{
	c5->cd();
	TLegend* leg5 = new TLegend(0.4, 0.2, 0.9, 0.45);
	FormLegend(leg5, 0.04);
	TH1D* htmp5 = new TH1D("htmp5", Form("%s;%s;ratio", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
	htmp5->SetAxisRange(0, 1.2, "Y");
	htmp5->Draw();
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		h2[ifile]->Draw("samepe");
		leg5->AddEntry(h2[ifile], Form("%s %d~", akaM[ifile].Data(), Range[ifile]), "p");
	}
	leg5->Draw();
	SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
	c5->SaveAs(Form("AfterScale/%s_eff_scaled_by_%s_Range%d_Coin%dth%d.pdf", VarName[ivar].Data(), akaD.Data(), Range[0], Coin, Th));
//}}}

//total eff zoom{{{
	c6->cd();
	TLegend* leg6 = new TLegend(0.4, 0.2, 0.9, 0.45);
	FormLegend(leg6, 0.04);
	TH1D* htmp6 = new TH1D("htmp6", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
	htmp6->SetAxisRange(0, VarMaxZ[ivar], "X");
	htmp6->SetAxisRange(0, 1.2, "Y");
	htmp6->Draw();
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		h2[ifile]->Draw("samepe");
		leg6->AddEntry(h2[ifile], Form("%s %d~", akaM[ifile].Data(), Range[ifile]), "p");
	}
	leg6->Draw();
	SetLine(1, 0, 1, VarMaxZ[ivar], 1, 0, 2);
	c6->SaveAs(Form("AfterScale/%s_eff_scaled_by_%s_Range%d_Coin%dth%d_zoom.pdf", VarName[ivar].Data(), akaD.Data(), Range[0], Coin, Th));
//}}}
}
