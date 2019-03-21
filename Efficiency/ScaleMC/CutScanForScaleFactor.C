//Headers{{{
#include "../../Utilities/Style_Header.h"
#include "../../Utilities/Var_Header.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TCut.h>
#include <TMath.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <TEfficiency.h>
#include <TSystem.h>
#include <TParameter.h>
//}}}

//External functions{{{
TH1D* ScaleTH1(TTree* t1, Int_t ivar, Int_t iiter, Double_t ScaleX, Double_t NormMin, Double_t NormMax, Int_t rebin)
{
//scale histogram{{{

//make new bin{{{
	Double_t tmpBinArr[nBinC[ivar]];
	tmpBinArr[ivar] = 0;
	Int_t NewnBin = 0;
	for(Int_t ibin = 1; ibin < nBinC[0]+2; ibin++)
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

	TH1D* h1 = new TH1D(Form("h1_%d", iiter), "", NewnBin-1, NewBinArr);
	t1->Draw(Form("%f*%s>>h1_%d", ScaleX, Vars[ivar].Data(), iiter), "");
	FormTH1(h1, 1);
	h1->Rebin(rebin);
	h1->Scale(1./( h1->Integral(h1->GetXaxis()->FindBin(NormMin), h1->GetXaxis()->FindBin(NormMax)) ));
	h1->Scale(1., "width");
	return h1;
//}}}
}

Double_t Getchi2(TH1D* href, TH1D* hcomp, Double_t xmin = -1, Double_t xmax = -1)
{
//calculate chi2{{{
	Double_t dev = 0;
	Int_t Ndof = 0;

	Double_t minbin = 1;
	Double_t maxbin = href->GetNbinsX()+1;

	if(xmin != -1) minbin = href->GetXaxis()->FindBin(xmin);
	if(xmax != -1) maxbin = href->GetXaxis()->FindBin(xmax);

	for(Int_t i = minbin; i < maxbin; i++)
	{
		Double_t y1 = href->GetBinContent(i);
		Double_t y2 = hcomp->GetBinContent(i);
		Double_t e1 = href->GetBinError(i);
		Double_t e2 = hcomp->GetBinError(i);

		Double_t dy = y1-y2;
		Double_t de = e1*e1+e2*e2;
		if(de > 0)
		{
			dev += dy*dy/de;
			Ndof += 1;
		}
	}
	return dev/Ndof;
//}}}
}
//}}}

void CutScanForScaleFactor(const Int_t ivar = 5, const Int_t MCN = 0)
{
	SetStyle();

//make directory{{{
	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	string saveDIR = mainDIR + "/Scaled";
	void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
	if (dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.c_str(), kTRUE);
//}}}

	TString akaD = "PR326483";
	TString akaM;
	TString Suffix = "newCNT";

//Set Names{{{
	TFile* fdata = new TFile(Form("../../HistFiles/PbPb2018_%s_histo_%s.root", akaD.Data(), Suffix.Data()), "READ");
	TString fMC;
	if(MCN == 0)
	{
		fMC = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/HydjetCymbal5F_Forest/HYDJET_CYMBAL5F_PbPb_5020GeV/HydjetCymbal5F_5020GeV_PbPb_Forest/181128_165108/HydjetCymbal5F_5020GeV_PbPb_Forest.root";
		akaM = "Hydjet_Cymbal5F";
	}
	else if(MCN == 1)
	{
		fMC = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/AMPT_StringMelting_Forest/AMPT_StringMelting_Forest.root";
		akaM = "AMPT_String";
	}
	else if(MCN == 2)
	{
		fMC = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/AMPT_noStringMelting_Forest/AMPT_No_StringMelting_Forest.root";
		akaM = "AMPT_NoString";
	}
	else if(MCN == 3)
	{
		fMC = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/EPOS_Forest/EPOS_Forest.root";
		akaM = "EPOS";
	}
	else
	{
		cout << "Out of MC sample" << endl;
		return;
	}
//}}}

//Rebin set{{{
	Int_t rebin = 1;
	if(ivar == 0 || ivar == 1 || ivar == 2 || ivar == 3 || ivar == 4 ||
		ivar == 5 || ivar == 6 || ivar == 7)
	{
		if(MCN == 1) rebin = 15;
		else if(MCN == 2) rebin = 25;
		else if(MCN == 3) rebin = 5;
	}
//}}}

//Get trees{{{
	TChain* t_evt = new TChain("hiEvtAnalyzer/HiTree");
	TChain* t_skim = new TChain("skimanalysis/HltTree");
	//TChain* t_hlt = new TChain("hltanalysis/HltTree");
	//TChain* t_trk = new TChain("ppTrack/trackTree");
	t_evt->Add(Form("%s", fMC.Data()));
	t_skim->Add(Form("%s", fMC.Data()));
	//t_hlt->Add(Form("%s", fMC.Data()));
	//t_trk->Add(Form("%s", fMC.Data()));
	t_evt->AddFriend(t_skim);
	//t_evt->AddFriend(t_hlt);
	//t_evt->AddFriend(t_trk);
//}}}

	const Int_t Niter = 60;
	FILE* ftxt;
	ftxt = fopen(Form("Scaled/ScaleFactor_scan_%s_%s_%s_%s.txt", VarName[ivar].Data(), akaD.Data(), akaM.Data(), Suffix.Data()), "w");
	TFile* fout = new TFile(Form("Scaled/Scaled_MC_eff_scan_%s_%s_%s_%s.root", VarName[ivar].Data(), akaD.Data(), akaM.Data(), Suffix.Data()), "RECREATE");
	fout->cd();

//get data{{{
	TH1D* href = (TH1D*) fdata->Get(Form("h%s_PV_CC_Coin2th4", VarName[ivar].Data()));
	FormTH1(href, 0);
	href->Rebin(rebin);
	href->Scale(1./href->Integral( href->GetXaxis()->FindBin(NormRangeMin[ivar]), href->GetXaxis()->FindBin(NormRangeMax[ivar]) ));
	href->Scale(1., "width");
//}}}

	for(Int_t ichi = 0; ichi < Nchi; ichi++)
	{
		TH1D* hMC[Niter];
		TGraphErrors* g1 = new TGraphErrors(Niter);
		Double_t chi2val[Niter];
		Double_t ScaleX[Niter];
		Double_t ScaleMin = 0.8;
		Double_t variation = 0.6;
		Int_t ibest = 0;

//scale iteration{{{
		for(Int_t iiter = 0; iiter < Niter; iiter++)
		{
			ScaleX[iiter] = ScaleMin+(variation/Niter)*iiter;
			hMC[iiter] = ScaleTH1(t_evt, ivar, iiter, ScaleX[iiter], NormRangeMin[ivar], NormRangeMax[ivar], rebin);
			chi2val[iiter] = Getchi2(hMC[iiter], href, Chi2Range[ivar][ichi]);
			g1->SetPoint(iiter, ScaleX[iiter], chi2val[iiter]);
			if(chi2val[iiter] < chi2val[ibest]) ibest = iiter;
		}

		for(Int_t iiter = 0; iiter < Niter; iiter++)
		{
			//delete unusing histograms in order to prevent resource consuming
			if(iiter != ibest) delete hMC[iiter];
		}
//}}}

		fprintf(ftxt, "%d: %f \t %f \n", (int)Chi2Range[ivar][ichi], ScaleX[ibest], chi2val[ibest]);

//Draw comp{{{
		TCanvas* c1 = new TCanvas("c1", "", 0, 0, 600, 600);
		c1->cd();
		c1->SetLogy();
		FormTH1(href, 0);
		href->Draw("hist");
		hMC[ibest]->Draw("samepe");
		c1->SaveAs(Form("Scaled/Scaled_MC_Data_Comp_scan_%s_cut%d_%s_%s_%s.pdf", VarName[ivar].Data(), (int)Chi2Range[ivar][ichi], akaD.Data(), akaM.Data(), Suffix.Data()));
//}}}

//Draw chi2{{{
		TCanvas* c2 = new TCanvas("c2", "", 0, 0, 600, 600);
		c2->cd();
		g1->SetTitle(";Scale factor;Reduced #chi^2 (normalized by ndf)");
		g1->GetXaxis()->CenterTitle();
		g1->GetYaxis()->CenterTitle();
		g1->SetMarkerStyle(20);
		g1->Draw("ap");
		c2->SaveAs(Form("Scaled/Scale_factor_chi2_scan_%s_cut%d_%s_%s_%s.pdf", VarName[ivar].Data(), (int)Chi2Range[ivar][ichi], akaD.Data(), akaM.Data(), Suffix.Data()));
//}}}

//Draw ratio{{{
		TCanvas* c3 = new TCanvas("c3", "", 0, 0, 600, 600);
		c3->cd();
		TH1D* hratio = (TH1D*) hMC[ibest]->Clone(Form("hratio_%s", VarName[ivar].Data()));
		hratio->SetAxisRange(0, 1.2, "Y");
		hratio->Divide(href);
		hratio->Draw("pe");
		SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
		c3->SaveAs(Form("Scaled/Scaled_ratio_scan_%s_cut%d_%s_%s_%s.pdf", VarName[ivar].Data(), (int)Chi2Range[ivar][ichi], akaD.Data(), akaM.Data(), Suffix.Data()));
//}}}

//Save plots{{{
		g1->SetName(Form("Scalefactor_%s_%d", VarName[ivar].Data(), (int)Chi2Range[ivar][ichi]));
		g1->Write();
		href->SetName(Form("h%s_ref%d", VarName[ivar].Data(), (int)Chi2Range[ivar][ichi]));
		href->Write();
		hMC[ibest]->SetName(Form("h%s_scaled%d", VarName[ivar].Data(), (int)Chi2Range[ivar][ichi]));
		hMC[ibest]->Write();
		hratio->Write();
//}}}
	}
	fclose(ftxt);
}
