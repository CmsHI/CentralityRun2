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
	TString aka[NOF] = {"PR326478", "Hydjet_Cymbal5F", "AMPT_String",
							"AMPT_NoString", "EPOS"};

	TFile* fin[NOF];
	TString PV[2] = {"PV", "2Trk"};
	TString Coin[2] = {"Coin3Th3", "Coin2Th4"};

//Define canvas{{{
	TCanvas* c1[NVar];
	TCanvas* c2[NVar];
	TCanvas* c3[NVar];
	TCanvas* c4[NVar];
	TCanvas* c5[NVar][2][2];
	TCanvas* c6[NVar][2][2];
	TCanvas* c7[NVar][2][2];
	TCanvas* c8[NVar][2][2];

	for(Int_t ivar = 0; ivar < NVar; ivar++)
	{
		c1[ivar] = new TCanvas(Form("c1_%s", VarName[ivar].Data()), "", 0, 0, 600, 600);
		c2[ivar] = new TCanvas(Form("c2_%s", VarName[ivar].Data()), "", 0, 0, 600, 600);
		c3[ivar] = new TCanvas(Form("c3_%s", VarName[ivar].Data()), "", 0, 0, 600, 600);
		c4[ivar] = new TCanvas(Form("c4_%s", VarName[ivar].Data()), "", 0, 0, 600, 600);
		for(Int_t i = 0; i < 2; i++)
		{
			for(Int_t j = 0; j < 2; j++)
			{
				c5[ivar][i][j] = new TCanvas(Form("c5_%s_%d_%d", VarName[ivar].Data(), i, j), "", 0, 0, 600, 600);
				c6[ivar][i][j] = new TCanvas(Form("c6_%s_%d_%d", VarName[ivar].Data(), i, j), "", 0, 0, 600, 600);
				c7[ivar][i][j] = new TCanvas(Form("c7_%s_%d_%d", VarName[ivar].Data(), i, j), "", 0, 0, 600, 600);
				c8[ivar][i][j] = new TCanvas(Form("c8_%s_%d_%d", VarName[ivar].Data(), i, j), "", 0, 0, 600, 600);
			}
		}
	}
//}}}

//Get histogram{{{
	TH1D* h1[NOF][NVar];
	TH1D* h2[NOF][NVar][2][2];
	TH1D* h1copy[NOF][NVar];
	TH1D* h2copy[NOF][NVar][2][2];

	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		fin[ifile] = new TFile(Form("../../HistFiles/PbPb2018_%s_histo_newCNT.root", aka[ifile].Data()), "READ");
		for(Int_t ivar = 0; ivar < NVar; ivar++)
		{
			h1[ifile][ivar] = (TH1D*) fin[ifile]->Get(Form("h%s_", VarName[ivar].Data()));
			h2[ifile][ivar][0][0] = (TH1D*) fin[ifile]->Get(Form("h%s_PV_CC_Coin3th3", VarName[ivar].Data()));
			h2[ifile][ivar][0][1] = (TH1D*) fin[ifile]->Get(Form("h%s_PV_CC_Coin2th4", VarName[ivar].Data()));
			h2[ifile][ivar][1][0] = (TH1D*) fin[ifile]->Get(Form("h%s_2Trk_CC_Coin3th3", VarName[ivar].Data()));
			h2[ifile][ivar][1][1] = (TH1D*) fin[ifile]->Get(Form("h%s_2Trk_CC_Coin2th4", VarName[ivar].Data()));
			h1copy[ifile][ivar] = (TH1D*) h1[ifile][ivar]->Clone();
			h2copy[ifile][ivar][0][0] = (TH1D*) h2[ifile][ivar][0][0]->Clone();
			h2copy[ifile][ivar][0][1] = (TH1D*) h2[ifile][ivar][0][1]->Clone();
			h2copy[ifile][ivar][1][0] = (TH1D*) h2[ifile][ivar][1][0]->Clone();
			h2copy[ifile][ivar][1][1] = (TH1D*) h2[ifile][ivar][1][1]->Clone();
			FormTH1Marker(h1copy[ifile][ivar], ifile, ifile, 1.4);
			Double_t ScaleFactor = h1copy[ifile][ivar]->Integral();
			h1copy[ifile][ivar]->Scale(1./ScaleFactor);
			h1copy[ifile][ivar]->Scale(1., "width");
			for(Int_t i = 0; i < 2; i++)
			{
				for(Int_t j = 0; j < 2; j++)
				{
					FormTH1Marker(h2copy[ifile][ivar][i][j], ifile, ifile, 1.4);
					h2copy[ifile][ivar][i][j]->Scale(1./ScaleFactor);
					h2copy[ifile][ivar][i][j]->Scale(1., "width");
				}
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

/*
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
*/

//total dist{{{
		c1[ivar]->cd();
		c1[ivar]->SetLogy();
		TLegend* leg1 = new TLegend(0.5, 0.65, 0.9, 0.9);
		FormLegend(leg1, 0.04);
		TH1D* htmp1 = new TH1D("htmp1", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp1->SetAxisRange(1e-9, 0.8, "Y");
		htmp1->Draw();
		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			h1copy[ifile][ivar]->Draw("psame");
			leg1->AddEntry(h1copy[ifile][ivar], Form("%s", aka[ifile].Data()), "p");
		}
		leg1->Draw();
		c1[ivar]->SaveAs(Form("BeforeScale/%s_dist_noselection_%s_newCNT.pdf", VarName[ivar].Data(), aka[0].Data()));
//}}}

//total dist zoom{{{
		c2[ivar]->cd();
		c2[ivar]->SetLogy();
		TLegend* leg2 = new TLegend(0.5, 0.65, 0.9, 0.9);
		FormLegend(leg2, 0.04);
		TH1D* htmp2 = new TH1D("htmp2", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp2->SetAxisRange(0, VarMaxZ[ivar], "X");
		htmp2->SetAxisRange(1e-9, 0.8, "Y");
		htmp2->Draw();
		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			h1copy[ifile][ivar]->Draw("psame");
			leg2->AddEntry(h1copy[ifile][ivar], Form("%s", aka[ifile].Data()), "p");
		}
		leg2->Draw();
		c2[ivar]->SaveAs(Form("BeforeScale/%s_dist_noselection_zoom_%s_newCNT.pdf", VarName[ivar].Data(), aka[0].Data()));
//}}}

//total ratio{{{
		c3[ivar]->cd();
		TLegend* leg3 = new TLegend(0.5, 0.65, 0.9, 0.9);
		FormLegend(leg3, 0.04);
		TH1D* htmp3 = new TH1D("htmp3", Form("%s;%s;ratio", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp3->SetAxisRange(0, 1.2, "Y");
		htmp3->Draw();
		for(Int_t ifile = 1; ifile < NOF; ifile++)
		{
			h1copy[ifile][ivar]->Divide(h1copy[ifile][ivar], h1copy[0][ivar], 1, 1, "b");
			h1copy[ifile][ivar]->Draw("psame");
			leg3->AddEntry(h1copy[ifile][ivar], Form("%s", aka[ifile].Data()), "p");
		}
		leg3->Draw();
		SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
		c3[ivar]->SaveAs(Form("BeforeScale/%s_ratio_noselection_%s_newCNT.pdf", VarName[ivar].Data(), aka[0].Data()));
//}}}

//total ratio zoom{{{
		c4[ivar]->cd();
		TLegend* leg4 = new TLegend(0.5, 0.65, 0.9, 0.9);
		FormLegend(leg4, 0.04);
		TH1D* htmp4 = new TH1D("htmp4", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp4->SetAxisRange(0, VarMaxZ[ivar], "X");
		htmp4->SetAxisRange(0, 1.2, "Y");
		htmp4->Draw();
		for(Int_t ifile = 1; ifile < NOF; ifile++)
		{
			h1copy[ifile][ivar]->Draw("psame");
			leg4->AddEntry(h1copy[ifile][ivar], Form("%s", aka[ifile].Data()), "p");
		}
		leg4->Draw();
		SetLine(1, 0, 1, NormRangeCut[ivar], 1, 0, 2);
		c4[ivar]->SaveAs(Form("BeforeScale/%s_ratio_noselectioni_zoom_%s_newCNT.pdf", VarName[ivar].Data(), aka[0].Data()));
//}}}

//selected dist{{{
		for(Int_t i = 0; i < 2; i++)
		{
			for(Int_t j = 0; j < 2; j++)
			{
				c5[ivar][i][j]->cd();
				c5[ivar][i][j]->SetLogy();
				TLegend* leg5 = new TLegend(0.5, 0.65, 0.9, 0.9);
				FormLegend(leg5, 0.04);
				TH1D* htmp5 = new TH1D("htmp5", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
				htmp5->SetAxisRange(1e-9, 0.8, "Y");
				htmp5->Draw();
				for(Int_t ifile = 0; ifile < NOF; ifile++)
				{
					h2copy[ifile][ivar][i][j]->Draw("psame");
					leg5->AddEntry(h2copy[ifile][ivar][i][j], Form("%s", aka[ifile].Data()), "p");
				}
				leg5->Draw();
				c5[ivar][i][j]->SaveAs(Form("BeforeScale/%s_dist_%s_%s_%s_newCNT.pdf", VarName[ivar].Data(), PV[i].Data(), Coin[j].Data(), aka[0].Data()));
			}
		}
//}}}

//selected dist zoom{{{
		for(Int_t i = 0; i < 2; i++)
		{
			for(Int_t j = 0; j < 2; j++)
			{
				c6[ivar][i][j]->cd();
				c6[ivar][i][j]->SetLogy();
				TH1D* htmp6 = new TH1D("htmp6", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
				htmp6->SetAxisRange(0, VarMaxZ[ivar], "X");
				htmp6->SetAxisRange(1e-9, 0.8, "Y");
				htmp6->Draw();
				TLegend* leg6 = new TLegend(0.5, 0.65, 0.9, 0.9);
				FormLegend(leg6, 0.04);
				for(Int_t ifile = 0; ifile < NOF; ifile++)
				{
					h2copy[ifile][ivar][i][j]->Draw("psame");
					leg6->AddEntry(h2copy[ifile][ivar][i][j], Form("%s", aka[ifile].Data()), "p");
				}
				leg6->Draw();
				c6[ivar][i][j]->SaveAs(Form("BeforeScale/%s_dist_%s_%s_zoom_%s_newCNT.pdf", VarName[ivar].Data(), PV[i].Data(), Coin[j].Data(), aka[0].Data()));
			}
		}
//}}}

//selected ratio{{{
		for(Int_t i = 0; i < 2; i++)
		{
			for(Int_t j = 0; j < 2; j++)
			{
				c7[ivar][i][j]->cd();
				TH1D* htmp7 = new TH1D("htmp7", Form("%s;%s;ratio", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
				htmp7->SetAxisRange(0, 1.2, "Y");
				htmp7->Draw();
				TLegend* leg7 = new TLegend(0.5, 0.65, 0.9, 0.9);
				FormLegend(leg7, 0.04);
				for(Int_t ifile = 1; ifile < NOF; ifile++)
				{
					h2copy[ifile][ivar][i][j]->Divide(h2copy[ifile][ivar][i][j], h2[0][ivar][i][j], 1, 1, "b");
					h2copy[ifile][ivar][i][j]->Draw("psame");
					leg7->AddEntry(h2copy[ifile][ivar][i][j], Form("%s", aka[ifile].Data()), "p");
				}
				leg7->Draw();
				SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
				c7[ivar][i][j]->SaveAs(Form("BeforeScale/%s_ratio_%s_%s_%s_newCNT.pdf", VarName[ivar].Data(), PV[i].Data(), Coin[j].Data(), aka[0].Data()));
			}
		}
//}}}

//selected ratio zoom{{{
		for(Int_t i = 0; i < 2; i++)
		{
			for(Int_t j = 0; j < 2; j++)
			{
				c8[ivar][i][j]->cd();
				TH1D* htmp8 = new TH1D("htmp8", Form("%s;%s;ratio", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
				htmp8->SetAxisRange(0, VarMaxZ[ivar], "X");
				htmp8->SetAxisRange(0, 1.2, "Y");
				htmp8->Draw();
				TLegend* leg8 = new TLegend(0.5, 0.65, 0.9, 0.9);
				FormLegend(leg8, 0.04);
				for(Int_t ifile = 1; ifile < NOF; ifile++)
				{
					h2copy[ifile][ivar][i][j]->Draw("psame");
					leg8->AddEntry(h2copy[ifile][ivar][i][j], Form("%s", aka[ifile].Data()), "p");
				}
				leg8->Draw();
				SetLine(1, 0, 1, NormRangeCut[ivar], 1, 0, 2);
				c8[ivar][i][j]->SaveAs(Form("BeforeScale/%s_ratio_%s_%s_zoom_%s_newCNT.pdf", VarName[ivar].Data(), PV[i].Data(), Coin[j].Data(), aka[0].Data()));
			}
		}
//}}}
	}
}
