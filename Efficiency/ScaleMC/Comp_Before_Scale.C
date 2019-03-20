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

void Comp_Before_Scale()
{
	SetStyle();

//make directory{{{
	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	string saveDIR = mainDIR + "/BeforeScale/";
	void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
	if (dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.c_str(), kTRUE);
//}}}

	const Int_t NOF = 5;//Number Of Files
	TString aka[NOF] = {"PR326483", "Hydjet_Cymbal5F", "AMPT_String",
							"AMPT_NoString", "EPOS"};

	TFile* fin[NOF];
	const Int_t NCut = 3;
	TString Coin[NCut] = {"total", "Coin2th4", "Coin3th3"};
	TString CutName[NCut] = {"", "PV_CC_Coin2th4", "PV_CC_Coin3th3"};
	TString Suffix = "newCNT";
	const Int_t rebin[NVar] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

//Define canvas{{{
	TCanvas* cdist[NVar][NCut];
	TCanvas* cdist_z[NVar][NCut];
	TCanvas* cratio[NVar][NCut];
	TCanvas* cratio_z[NVar][NCut];

	for(Int_t ivar = 0; ivar < NVar; ivar++)
	{
		for(Int_t icut = 0; icut < NCut; icut++)
		{
			cdist[ivar][icut] = new TCanvas(Form("cdist_%s_%d", VarName[ivar].Data(), icut), "", 0, 0, 600, 600);
			cdist_z[ivar][icut] = new TCanvas(Form("cdist_z_%s_%d", VarName[ivar].Data(), icut), "", 0, 0, 600, 600);
			cratio[ivar][icut] = new TCanvas(Form("cratio_%s_%d", VarName[ivar].Data(), icut), "", 0, 0, 600, 600);
			cratio_z[ivar][icut] = new TCanvas(Form("cratio_z_%s_%d", VarName[ivar].Data(), icut), "", 0, 0, 600, 600);
		}
	}
//}}}

//Get histogram{{{
	TH1D* h1[NOF][NVar][NCut];
	TH1D* h1copy[NOF][NVar][NCut];

	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		fin[ifile] = new TFile(Form("../../HistFiles/PbPb2018_%s_histo_%s.root", aka[ifile].Data(), Suffix.Data()), "READ");
		for(Int_t ivar = 0; ivar < NVar; ivar++)
		{
			Double_t NormFactor;
			for(Int_t icut = 0; icut < NCut; icut++)
			{
				h1[ifile][ivar][icut] = (TH1D*) fin[ifile]->Get(Form("h%s_%s", VarName[ivar].Data(), CutName[icut].Data()));
				h1copy[ifile][ivar][icut] = (TH1D*) h1[ifile][ivar][icut]->Clone();
				FormTH1Marker(h1copy[ifile][ivar][icut], ifile, ifile, 1.4);
				if(icut == 0) NormFactor = h1copy[ifile][ivar][icut]->Integral( h1copy[ifile][ivar][icut]->GetXaxis()->FindBin(NormRangeMin[ivar]), h1copy[ifile][ivar][icut]->GetXaxis()->FindBin(NormRangeMax[ivar]) );
				h1copy[ifile][ivar][icut]->Scale(1./NormFactor);
				h1copy[ifile][ivar][icut]->Scale(1., "width");
			}
		}
	}
//}}}

	for(Int_t ivar = 0; ivar < NVar; ivar++)
	{
//make new bin{{{
		Double_t tmpBinArr[nBinC[ivar]];
		tmpBinArr[0] = 0;
		Int_t NewnBin = 0;
		for(Int_t ibin = 1; ibin < nBinC[ivar]+2; ibin++)
		{
			if(tmpBinArr[ibin-1] < BinBoundary[ivar]) tmpBinArr[ibin] = tmpBinArr[ibin-1]+( VarMaxC[ivar]/(double)nBinC[ivar] );
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
			if(NewBinArr[ibin-1] < BinBoundary[ivar]) NewBinArr[ibin] = NewBinArr[ibin-1]+( VarMaxC[ivar]/(double)nBinC[ivar] );
			else NewBinArr[ibin] = NewBinArr[ibin-1]+( VarW[ivar]*VarMaxC[ivar]/(double)nBinC[ivar] );
		}
//}}}

/*
//Set error 0{{{
		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			for(Int_t ibin = 0; ibin < nBinC[ivar]; ibin++)
			{
				h1copy[ifile][ivar]->SetBinError(ibin, 0.);
				h2copy[ifile][ivar][0][0]->SetBinError(ibin, 0.);
				h2copy[ifile][ivar][0][1]->SetBinError(ibin, 0.);
				h2copy[ifile][ivar][1][0]->SetBinError(ibin, 0.);
				h2copy[ifile][ivar][1][1]->SetBinError(ibin, 0.);
			}
		}
//}}}
*/

//draw dist{{{
		TLegend* leg_dist[NCut];
		TH1D* htmp1 = new TH1D("htmp1", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp1->SetAxisRange(1e-9, 0.8, "Y");
		for(Int_t icut = 0; icut < NCut; icut++)
		{
			cdist[ivar][icut]->cd();
			cdist[ivar][icut]->SetLogy();
			leg_dist[icut] = new TLegend(0.5, 0.65, 0.9, 0.9);
			FormLegend(leg_dist[icut], 0.04);
			htmp1->Draw();
			TH1D* h1rebin[NOF];
			for(Int_t ifile = 0; ifile < NOF; ifile++)
			{
				h1rebin[ifile] = (TH1D*) h1copy[ifile][ivar][icut]->Clone();
				h1rebin[ifile]->Rebin(rebin[ivar]);
				h1rebin[ifile]->Draw("pesame");
				leg_dist[icut]->AddEntry(h1rebin[ifile], Form("%s", aka[ifile].Data()), "p");
			}
			leg_dist[icut]->Draw();
			cdist[ivar][icut]->SaveAs(Form("BeforeScale/%s_dist_%s_%s_%s.pdf", VarName[ivar].Data(), Coin[icut].Data(), aka[0].Data(), Suffix.Data()));
		}
//}}}

//draw dist zoom{{{
		TLegend* leg_dist_z[NCut];
		TH1D* htmp2 = new TH1D("htmp2", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp2->SetAxisRange(0, VarMaxZ[ivar], "X");
		htmp2->SetAxisRange(1e-9, 0.8, "Y");
		for(Int_t icut = 0; icut < NCut; icut++)
		{
			leg_dist_z[icut] = new TLegend(0.5, 0.65, 0.9, 0.9);
			FormLegend(leg_dist_z[icut], 0.04);
			cdist_z[ivar][icut]->cd();
			cdist_z[ivar][icut]->SetLogy();
			htmp2->Draw();
			for(Int_t ifile = 0; ifile < NOF; ifile++)
			{
				h1copy[ifile][ivar][icut]->Draw("pesame");
				leg_dist_z[icut]->AddEntry(h1copy[ifile][ivar][icut], Form("%s", aka[ifile].Data()), "p");
				leg_dist_z[icut]->Draw();
			}
			cdist_z[ivar][icut]->SaveAs(Form("BeforeScale/%s_dist_%s_zoom_%s_%s.pdf", VarName[ivar].Data(), Coin[icut].Data(), aka[0].Data(), Suffix.Data()));
		}
//}}}

//draw ratio{{{
		TLegend* leg_ratio[NCut];
		TH1D* htmp3 = new TH1D("htmp3", Form("%s;%s;ratio", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp3->SetAxisRange(0, 1.2, "Y");
		for(Int_t icut = 0; icut < NCut; icut++)
		{
			leg_ratio[icut] = new TLegend(0.5, 0.65, 0.9, 0.9);
			FormLegend(leg_ratio[icut], 0.04);
			cratio[ivar][icut]->cd();
			htmp3->Draw();
			for(Int_t ifile = 1; ifile < NOF; ifile++)
			{
				h1copy[ifile][ivar][icut]->Divide(h1copy[ifile][ivar][icut], h1copy[0][ivar][icut], 1, 1, "b");
				h1copy[ifile][ivar][icut]->Draw("pesame");
				leg_ratio[icut]->AddEntry(h1copy[ifile][ivar][icut], Form("%s", aka[ifile].Data()), "p");
			}
			leg_ratio[icut]->Draw();
			SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
			cratio[ivar][icut]->SaveAs(Form("BeforeScale/%s_ratio_%s_%s_%s.pdf", VarName[ivar].Data(), Coin[icut].Data(), aka[0].Data(), Suffix.Data()));
		}
//}}}

//total ratio zoom{{{
		TLegend* leg_ratio_z[NCut];
		TH1D* htmp4 = new TH1D("htmp4", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp4->SetAxisRange(0, VarMaxZ[ivar], "X");
		htmp4->SetAxisRange(0, 1.2, "Y");
		for(Int_t icut = 0; icut < NCut; icut++)
		{
			leg_ratio_z[icut] = new TLegend(0.5, 0.65, 0.9, 0.9);
			FormLegend(leg_ratio_z[icut], 0.04);
			cratio_z[ivar][icut]->cd();
			htmp4->Draw();
			for(Int_t ifile = 1; ifile < NOF; ifile++)
			{
				h1copy[ifile][ivar][icut]->Draw("pesame");
				leg_ratio_z[icut]->AddEntry(h1copy[ifile][ivar][icut], Form("%s", aka[ifile].Data()), "p");
			}
			leg_ratio_z[icut]->Draw();
			SetLine(1, 0, 1, BinBoundary[ivar], 1, 0, 2);
			cratio_z[ivar][icut]->SaveAs(Form("BeforeScale/%s_ratio_%s_zoom_%s_%s.pdf", VarName[ivar].Data(), Coin[icut].Data(), aka[0].Data(), Suffix.Data()));
		}
//}}}
	}
}
