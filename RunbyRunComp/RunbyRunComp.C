//Headers{{{
#include "../Utilities/Style_Header.h"
#include "../Utilities/Var_Header.h"
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
#include <TSystem.h>
#include <TParameter.h>
//}}}

void RunbyRunComp()
{
	SetStyle();

//make directory{{{
	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());

	string CompDIR = mainDIR + "/Comp/";
	void * dirc = gSystem->OpenDirectory(CompDIR.c_str());
	if (dirc) gSystem->FreeDirectory(dirc);
	else gSystem->mkdir(CompDIR.c_str(), kTRUE);

	string RatioDIR = mainDIR + "/Ratio/";
	void * dirr = gSystem->OpenDirectory(RatioDIR.c_str());
	if (dirr) gSystem->FreeDirectory(dirr);
	else gSystem->mkdir(RatioDIR.c_str(), kTRUE);
//}}}

	TString version = "v8";
	TString aka1 = "PR326617";
	TFile* fref = new TFile(Form("../HistFiles/PbPb2018_%s_histo.root", aka1.Data()), "READ");
	const Int_t NOF = 2;//Number Of Files to be compared
	TString aka2[NOF] = {"PR326622", "PR326790"};

	TFile* fcomp[NOF];
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		fcomp[ifile] = new TFile(Form("../HistFiles/PbPb2018_%s_histo.root", aka2[ifile].Data()), "READ");
	}

	TH1D* hRef[NVar][NSel];
	TH1D* hComp[NOF][NVar][NSel];

//Get histos from file and draw{{{
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

		for(Int_t isel = 0; isel < NSel; isel++)
		{
			hRef[ivar][isel] = (TH1D*) fref->Get(Form("h%s_%s", VarName[ivar].Data(), SelName[isel].Data()));
			FormTH1Marker(hRef[ivar][isel], 0, 0, 1.4);
			hRef[ivar][isel]->Scale(1./hRef[ivar][isel]->Integral());
			hRef[ivar][isel]->Scale(1., "width");
			for(Int_t ifile = 0; ifile < NOF; ifile++)
			{
				hComp[ifile][ivar][isel] = (TH1D*) fcomp[ifile]->Get(Form("h%s_%s", VarName[ivar].Data(), SelName[isel].Data()));
				FormTH1Marker(hComp[ifile][ivar][isel], ifile+1, ifile+1, 1.4);
				hComp[ifile][ivar][isel]->Scale(1./hComp[ifile][ivar][isel]->Integral());
				hComp[ifile][ivar][isel]->Scale(1., "width");
			}

//dist{{{
			TCanvas* c1 = new TCanvas("c1", "", 0, 0, 600, 600);
			c1->cd();
			c1->SetLogy();
			TLegend* leg1 = new TLegend(0.57, 0.75, 0.9, 0.93);
			FormLegend(leg1, 0.04);
			TH1D* htmp1 = new TH1D("htmp1", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
			htmp1->SetAxisRange(1e-9, 0.8, "Y");
			htmp1->Draw();
			hRef[ivar][isel]->Draw("samepe");
			leg1->AddEntry(hRef[ivar][isel], Form("%s", aka1.Data()), "p");
			for(Int_t ifile = 0; ifile < NOF; ifile++)
			{
				hComp[ifile][ivar][isel]->Draw("samepe");
				leg1->AddEntry(hComp[ifile][ivar][isel], Form("%s", aka2[ifile].Data()), "p");
			}
			leg1->Draw();
			c1->SaveAs(Form("Comp/RunbyRunComp_%s_%s_%s.pdf", VarName[ivar].Data(),  SelName[isel].Data(), version.Data()));
//}}}

//dist zoom{{{
			TCanvas* c2 = new TCanvas("c2", "", 0, 0, 600, 600);
			c2->cd();
			c2->SetLogy();
			TLegend* leg2 = new TLegend(0.57, 0.75, 0.9, 0.93);
			FormLegend(leg2, 0.04);
			TH1D* htmp2 = new TH1D("htmp2", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
			htmp2->SetAxisRange(0, NormRangeCut[ivar], "X");
			htmp2->SetAxisRange(1e-9, 0.8, "Y");
			htmp2->Draw();
			hRef[ivar][isel]->Draw("samepe");
			leg2->AddEntry(hRef[ivar][isel], Form("%s", aka1.Data()), "p");
			for(Int_t ifile = 0; ifile < NOF; ifile++)
			{
				hComp[ifile][ivar][isel]->Draw("samepe");
				leg2->AddEntry(hComp[ifile][ivar][isel], Form("%s", aka2[ifile].Data()), "p");
			}
			leg2->Draw();
			c2->SaveAs(Form("Comp/RunbyRunComp_%s_%s_%s_zoom.pdf", VarName[ivar].Data(),  SelName[isel].Data(), version.Data()));
//}}}

//ratio{{{
			TCanvas* c3 = new TCanvas("c3", "", 0, 0, 600, 600);
			c3->cd();
			TLegend* leg3 = new TLegend(0.4, 0.25, 0.9, 0.43);
			FormLegend(leg3, 0.04);
			TH1D* htmp3 = new TH1D("htmp3", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
			htmp3->SetAxisRange(0, 1.5, "Y");
			htmp3->Draw();
			for(Int_t ifile = 0; ifile < NOF; ifile++)
			{
				hComp[ifile][ivar][isel]->Divide(hComp[ifile][ivar][isel], hRef[ivar][isel], 1, 1, "b");
				hComp[ifile][ivar][isel]->Draw("samepe");
				leg3->AddEntry(hComp[ifile][ivar][isel], Form("%s/%s", aka2[ifile].Data(), aka1.Data()));
			}
			leg3->Draw();
			SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
			c3->SaveAs(Form("Ratio/RunbyRunRatio_%s_%s_%s.pdf", VarName[ivar].Data(), SelName[isel].Data(), version.Data()));
//}}}

//ratio zoom{{{
			TCanvas* c4 = new TCanvas("c4", "", 0, 0, 600, 600);
			c4->cd();
			TLegend* leg4 = new TLegend(0.4, 0.25, 0.9, 0.43);
			FormLegend(leg4, 0.04);
			TH1D* htmp4 = new TH1D("htmp4", Form("%s;%s;", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
			htmp4->SetAxisRange(0, NormRangeCut[ivar], "X");
			htmp4->SetAxisRange(0, 1.5, "Y");
			htmp4->Draw();
			for(Int_t ifile = 0; ifile < NOF; ifile++)
			{
				hComp[ifile][ivar][isel]->Draw("samepe");
				leg4->AddEntry(hComp[ifile][ivar][isel], Form("%s/%s", aka2[ifile].Data(), aka1.Data()));
			}
			leg4->Draw();
			SetLine(1, 0, 1, NormRangeCut[ivar], 1, 0, 2);
			c4->SaveAs(Form("Ratio/RunbyRunRatio_%s_%s_%s_zoom.pdf", VarName[ivar].Data(), SelName[isel].Data(), version.Data()));
//}}}

		}
	}
//}}}
}
