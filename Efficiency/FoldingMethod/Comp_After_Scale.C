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

void Comp_After_Scale()
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
	TString akaD = "PR326617";
	TString akaM[NOF] = {"Hydjet_Cymbal5F", "AMPT_String", "AMPT_NoString",
								"EPOS"};

	TFile* fref;
	TFile* fin[NOF][2];
	TCanvas* c1[NVar][2];
	TCanvas* c2[NVar][2];
	TCanvas* c3[NVar][2];
	TCanvas* c4[NVar][2];
	TH1D* href[NVar][2];
	TH1D* h1[NOF][NVar][2];
	TString Coin1[2] = {"Coin3th3", "Coin2th4"};
	TString Coin2[2] = {"coin3th3", "coin2th4"};

//Define canvas{{{
	for(Int_t ivar = 0; ivar < NVar; ivar++)
	{
		for(Int_t i = 0; i < 2; i++)
		{
			c1[ivar][i] = new TCanvas(Form("c1_%s_%d", VarName[ivar].Data(), i), "", 0, 0, 600, 600);
			c2[ivar][i] = new TCanvas(Form("c2_%s_%d", VarName[ivar].Data(), i), "", 0, 0, 600, 600);
			c3[ivar][i] = new TCanvas(Form("c3_%s_%d", VarName[ivar].Data(), i), "", 0, 0, 600, 600);
			c4[ivar][i] = new TCanvas(Form("c4_%s_%d", VarName[ivar].Data(), i), "", 0, 0, 600, 600);
		}
	}
//}}}

	for(Int_t ivar = 0; ivar < NVar; ivar++)
	{
		fref = new TFile(Form("../../HistFiles/PbPb2018_%s_histo.root", akaD.Data()), "READ");

//get reference hist{{{
		for(Int_t i = 0; i < 2; i++)
		{
			href[ivar][i] = (TH1D*) fref->Get(Form("h%s_PV_CC_%s", VarName[ivar].Data(), Coin1[i].Data()));
			FormTH1Marker(href[ivar][i], 0, 0, 1.4);
			href[ivar][i]->Scale(1./( href[ivar][i]->Integral(href[ivar][i]->GetXaxis()->FindBin(NormRangeCut[ivar]), href[ivar][i]->GetXaxis()->FindBin(VarMaxC[ivar])) ));
			href[ivar][i]->Scale(1., "width");
		}
//}}}

//get hist{{{
		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			for(Int_t i = 0; i < 2; i++)
			{
				fin[ifile][i] = new TFile(Form("MCefficiency/MC_eff_2018_%s_%s_%s_by_%s.root", VarName[ivar].Data(), Coin2[i].Data(), akaM[ifile].Data(), akaD.Data()), "READ");
				h1[ifile][ivar][i] = (TH1D*) fin[ifile][i]->Get(Form("hsel_%s", VarName[ivar].Data()));
				FormTH1Marker(h1[ifile][ivar][i], ifile+1, ifile+1, 1.4);
				h1[ifile][ivar][i]->Scale(1./( h1[ifile][ivar][i]->Integral(h1[ifile][ivar][i]->GetXaxis()->FindBin(NormRangeCut[ivar]), h1[ifile][ivar][i]->GetXaxis()->FindBin(VarMaxC[ivar])) ));
				h1[ifile][ivar][i]->Scale(1., "width");
			}
		}
//}}}
	}

	for(Int_t ivar = 0; ivar < NVar; ivar++)
	{
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

		for(Int_t i = 0; i < 2; i++)
		{
//total dist{{{
			c1[ivar][i]->cd();
			c1[ivar][i]->SetLogy();
			TLegend* leg1 = new TLegend(0.55, 0.7, 0.9, 0.9);
			FormLegend(leg1, 0.04);
			TH1D* htmp1 = new TH1D("htmp1", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
			htmp1->SetAxisRange(1e-9, 0.8, "Y");
			htmp1->Draw();
			href[ivar][i]->Draw("samepe");
			leg1->AddEntry(href[ivar][i], Form("%s", akaD.Data()), "p");
			for(Int_t ifile = 0; ifile < NOF; ifile++)
			{
				h1[ifile][ivar][i]->Draw("samepe");
				leg1->AddEntry(h1[ifile][ivar][i], Form("%s", akaM[ifile].Data()), "p");
			}
			leg1->Draw();
			c1[ivar][i]->SaveAs(Form("AfterScale/%s_dist_scaled_by_%s_%s.pdf", VarName[ivar].Data(), akaD.Data(), Coin1[i].Data()));
//}}}

//total dist zoom{{{
			c2[ivar][i]->cd();
			c2[ivar][i]->SetLogy();
			TLegend* leg2 = new TLegend(0.55, 0.7, 0.9, 0.9);
			FormLegend(leg2, 0.04);
			TH1D* htmp2 = new TH1D("htmp2", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
			htmp2->SetAxisRange(0, NormRangeCut[ivar], "X");
			htmp2->SetAxisRange(1e-9, 0.8, "Y");
			htmp2->Draw();
			href[ivar][i]->Draw("samepe");
			leg2->AddEntry(href[ivar][i], Form("%s", akaD.Data()), "p");
			for(Int_t ifile = 0; ifile < NOF; ifile++)
			{
				h1[ifile][ivar][i]->Draw("samepe");
				leg2->AddEntry(h1[ifile][ivar][i], Form("%s", akaM[ifile].Data()), "p");
			}
			leg2->Draw();
			c2[ivar][i]->SaveAs(Form("AfterScale/%s_dist_scaled_by_%s_%s_zoom.pdf", VarName[ivar].Data(), akaD.Data(), Coin1[i].Data()));
//}}}

//total ratio{{{
		c3[ivar][i]->cd();
		TLegend* leg3 = new TLegend(0.55, 0.7, 0.9, 0.9);
		FormLegend(leg3, 0.04);
		TH1D* htmp3 = new TH1D("htmp3", Form("%s;%s;ratio", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp3->SetAxisRange(0, 1.2, "Y");
		htmp3->Draw();
		for(Int_t ifile = 1; ifile < NOF; ifile++)
		{
			h1[ifile][ivar][i]->Divide(h1[ifile][ivar][i], href[ivar][i], 1, 1, "b");
			h1[ifile][ivar][i]->Draw("samepe");
			leg3->AddEntry(h1[ifile][ivar][i], Form("%s", akaM[ifile].Data()), "p");
		}
		leg3->Draw();
		SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
		c3[ivar][i]->SaveAs(Form("AfterScale/%s_ratio_scaled_by_%s_%s.pdf", VarName[ivar].Data(), akaD.Data(), Coin1[i].Data()));
//}}}

//total ratio zoom{{{
			c4[ivar][i]->cd();
			TLegend* leg4 = new TLegend(0.55, 0.7, 0.9, 0.9);
			FormLegend(leg4, 0.04);
			TH1D* htmp4 = new TH1D("htmp4", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
			htmp4->SetAxisRange(0, NormRangeCut[ivar], "X");
			htmp4->SetAxisRange(0, 1.2, "Y");
			htmp4->Draw();
			for(Int_t ifile = 1; ifile < NOF; ifile++)
			{
				h1[ifile][ivar][i]->Draw("samepe");
				leg4->AddEntry(h1[ifile][ivar][i], Form("%s", akaM[ifile].Data()), "p");
			}
			leg4->Draw();
			SetLine(1, 0, 1, NormRangeCut[ivar], 1, 0, 2);
			c4[ivar][i]->SaveAs(Form("AfterScale/%s_ratio_scaled_by_%s_%s_zoom.pdf", VarName[ivar].Data(), akaD.Data(), Coin1[i].Data()));
//}}}
		}
	}
}
