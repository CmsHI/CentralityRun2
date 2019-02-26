//Headers{{{
#include "../../Utilities/Style_Header.h"
#include "../../Utilities/Var_Header.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TGraphAsymmErrors.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <TEfficiency.h>
#include <TSystem.h>
#include <TParameter.h>
//}}}

void MCeffComp()
{
	SetStyle();

//make directory{{{
	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	string saveDIR = mainDIR + "/MCefficiency/";
	void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
	if (dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.c_str(), kTRUE);
//}}}

	const Int_t NOF = 4;//Number of Files
	TFile* fin[NOF];
	TString aka[NOF] = {"Hydjet_Cymbal5F", "AMPT_String", "AMPT_NoString",
							"EPOS"};

//get hist{{{
	TH1D* htot[NOF][NVar];
	TH1D* h1[NOF][NVar][2];

	Double_t ntot[NOF][NVar];
	Double_t n1[NOF][NVar][2];

	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		fin[ifile] = new TFile(Form("../../HistFiles/PbPb2018_%s_histo_newCNT.root", aka[ifile].Data()));

		for(Int_t ivar = 0; ivar < NVar; ivar++)
		{
			htot[ifile][ivar] = (TH1D*) fin[ifile]->Get(Form("h%s_", VarName[ivar].Data()));
			h1[ifile][ivar][0] = (TH1D*) fin[ifile]->Get(Form("h%s_PV_CC_Coin3th3", VarName[ivar].Data()));
			h1[ifile][ivar][1] = (TH1D*) fin[ifile]->Get(Form("h%s_PV_CC_Coin2th4", VarName[ivar].Data()));
			FormTH1Marker(htot[ifile][ivar], 0, 0, 1.4);
			FormTH1Marker(h1[ifile][ivar][0], ifile, ifile, 1.4);
			FormTH1Marker(h1[ifile][ivar][1], ifile, ifile, 1.4);
			ntot[ifile][ivar] = htot[ifile][ivar]->GetEntries();
			n1[ifile][ivar][0] = h1[ifile][ivar][0]->GetEntries();
			n1[ifile][ivar][1] = h1[ifile][ivar][1]->GetEntries();
			h1[ifile][ivar][0]->Divide(h1[ifile][ivar][0], htot[ifile][ivar], 1, 1, "b");
			h1[ifile][ivar][1]->Divide(h1[ifile][ivar][1], htot[ifile][ivar], 1, 1, "b");

		}
	}
//}}}

	TLatex* lt1 = new TLatex();
	FormLatex(lt1, 12, 0.04);

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

//Coin3th3 cut{{{
		TCanvas* c1 = new TCanvas("c1", "", 0, 0, 600, 600);
		c1->cd();
		TLegend* leg1 = new TLegend(0.5, 0.6, 0.9, 0.9);
		FormLegend(leg1, 0.04);

		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			if(ifile == 0) h1[ifile][ivar][0]->Draw("pe");
			else h1[ifile][ivar][0]->Draw("samepe");
			leg1->AddEntry(h1[ifile][ivar][0], Form("%s", aka[ifile].Data()), "p");
			lt1->DrawLatex(0.15, 0.55 - ( (double) ifile)*0.05, Form("%s eff. = %.3f %%", aka[ifile].Data(), 100*( n1[ifile][ivar][0]/ntot[ifile][ivar] )));
		}
		leg1->Draw();
		SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
		c1->SaveAs(Form("MCefficiency/MCCompeff_dist_%s_Coin3th3_newCNT_fine.pdf", VarName[ivar].Data()));

		TCanvas* c2 = new TCanvas("c2", "", 0, 0, 600, 600);
		c2->cd();
		TLegend* leg2 = new TLegend(0.5, 0.6, 0.9, 0.9);
		FormLegend(leg2, 0.04);
		TH1D* htmp1 = new TH1D("htmp1", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp1->SetAxisRange(0, VarMaxZ[ivar], "X");
		htmp1->SetAxisRange(0, 1.2, "Y");
		htmp1->Draw();

		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			h1[ifile][ivar][0]->Draw("samepe");
			leg2->AddEntry(h1[ifile][ivar][0], Form("%s", aka[ifile].Data()), "p");
			lt1->DrawLatex(0.15, 0.55 - ( (double) ifile)*0.05, Form("%s eff. = %.3f %%", aka[ifile].Data(), 100*( n1[ifile][ivar][0]/ntot[ifile][ivar] )));
		}
		leg2->Draw();
		SetLine(1, 0, 1, VarMaxZ[ivar], 1, 0, 2);
		c2->SaveAs(Form("MCefficiency/MCCompeff_dist_%s_Coin3th3_newCNT_zoom_fine.pdf", VarName[ivar].Data()));
//}}}

//Coin2th4 cut{{{
		TCanvas* c3 = new TCanvas("c3", "", 0, 0, 600, 600);
		c3->cd();
		TLegend* leg3 = new TLegend(0.5, 0.6, 0.9, 0.9);
		FormLegend(leg3, 0.04);

		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			if(ifile == 0) h1[ifile][ivar][1]->Draw("pe");
			else h1[ifile][ivar][1]->Draw("samepe");
			leg3->AddEntry(h1[ifile][ivar][1], Form("%s", aka[ifile].Data()), "p");
			lt1->DrawLatex(0.15, 0.55 - ( (double) ifile)*0.05, Form("%s eff. = %.3f %%", aka[ifile].Data(), 100*( n1[ifile][ivar][1]/ntot[ifile][ivar] )));
		}
		leg3->Draw();
		SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
		c3->SaveAs(Form("MCefficiency/MCCompeff_dist_%s_Coin2th4_newCNT_fine.pdf", VarName[ivar].Data()));

		TCanvas* c4 = new TCanvas("c4", "", 0, 0, 600, 600);
		c4->cd();
		TLegend* leg4 = new TLegend(0.5, 0.6, 0.9, 0.9);
		FormLegend(leg4, 0.04);
		TH1D* htmp2 = new TH1D("htmp2", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp2->SetAxisRange(0, VarMaxZ[ivar], "X");
		htmp2->SetAxisRange(0, 1.2, "Y");
		htmp2->Draw();

		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			h1[ifile][ivar][1]->Draw("samepe");
			leg4->AddEntry(h1[ifile][ivar][1], Form("%s", aka[ifile].Data()), "p");
			lt1->DrawLatex(0.15, 0.55 - ( (double) ifile)*0.05, Form("%s eff. = %.3f %%", aka[ifile].Data(), 100*( n1[ifile][ivar][1]/ntot[ifile][ivar] )));
		}
		leg4->Draw();
		SetLine(1, 0, 1, VarMaxZ[ivar], 1, 0, 2);
		c4->SaveAs(Form("MCefficiency/MCCompeff_dist_%s_Coin2th4_newCNT_zoom_fine.pdf", VarName[ivar].Data()));
//}}}
	}
}
