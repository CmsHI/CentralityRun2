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
TH1D* ScaleTH1(TTree* t1, Int_t iiter, Double_t ScaleX, Int_t ivar)
{
//scale histogram{{{

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

	TH1D* h1 = new TH1D(Form("h1_%d", iiter), "", NewnBin-1, NewBinArr);
	t1->Draw(Form("%f*%s>>h1_%d", ScaleX, Vars[ivar].Data(), iiter), "");
	FormTH1(h1, 1);
	h1->Scale(1./( h1->Integral(h1->GetXaxis()->FindBin(NormRangeCut[ivar]), h1->GetXaxis()->FindBin(VarMaxC[ivar])) ));
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

void GetScaleFactor(const Int_t ivar = 5)
{
	SetStyle();

//make directory{{{
	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	string saveDIR = mainDIR + "/Scaled";
	void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
	if (dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.c_str(), kTRUE);
//}}}

	TString akaD = "PR326478";
	TString akaM = "EPOS";
	TFile* fdata = new TFile(Form("../../HistFiles/PbPb2018_%s_histo_newCNT.root", akaD.Data()), "READ");
	//TString fMC = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/HydjetCymbal5F_Forest/HYDJET_CYMBAL5F_PbPb_5020GeV/HydjetCymbal5F_5020GeV_PbPb_Forest/181128_165108/HydjetCymbal5F_5020GeV_PbPb_Forest.root";//Hydjet Cymbal5F
	//TString fMC = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/AMPT_StringMelting_Forest/AMPT_StringMelting_Forest.root";//AMPT String Melting
	//TString fMC = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/AMPT_noStringMelting_Forest/AMPT_No_StringMelting_Forest.root";//AMPT No String Melting
	TString fMC = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/HydjetCymbal5F_Forest/HYDJET_CYMBAL5F_PbPb_5020GeV_Forest.root";//EPOS

//Get trees{{{
	TChain* t_evt = new TChain("hiEvtAnalyzer/HiTree");
	TChain* t_skim = new TChain("skimanalysis/HltTree");
	TChain* t_hlt = new TChain("hltanalysis/HltTree");
	TChain* t_trk = new TChain("ppTrack/trackTree");
	t_evt->Add(Form("%s", fMC.Data()));
	t_skim->Add(Form("%s", fMC.Data()));
	t_hlt->Add(Form("%s", fMC.Data()));
	t_trk->Add(Form("%s", fMC.Data()));
	t_evt->AddFriend(t_skim);
	t_evt->AddFriend(t_hlt);
	t_evt->AddFriend(t_trk);
//}}}

	TH1D* href = (TH1D*) fdata->Get(Form("h%s_", VarName[ivar].Data()));
	FormTH1(href, 0);
	href->Scale(1./( href->Integral(href->GetXaxis()->FindBin(NormRangeCut[ivar]), href->GetXaxis()->FindBin(VarMaxC[ivar])) ));
	href->Scale(1., "width");

	const Int_t Niter = 60;
	TH1D* hMC[Niter];
	TGraphErrors* g1 = new TGraphErrors(Niter);
	Double_t chi2val[Niter];
	Double_t ScaleX[Niter];
	Double_t ScaleMin = 0.8;
	Double_t variation = 0.6;
	if(ivar == 9) ScaleMin = 1.5;
	Int_t ibest = 0;
	const Double_t RangeCut = 500;

//scale iteration{{{
	for(Int_t iiter = 0; iiter < Niter; iiter++)
	{
		ScaleX[iiter] = ScaleMin+(variation/Niter)*iiter;
		hMC[iiter] = ScaleTH1(t_evt, iiter, ScaleX[iiter], ivar);
		chi2val[iiter] = Getchi2(hMC[iiter], href, RangeCut);
		g1->SetPoint(iiter, ScaleX[iiter], chi2val[iiter]);
		if(chi2val[iiter] < chi2val[ibest]) ibest = iiter;
	}
	for(Int_t iiter = 0; iiter < Niter; iiter++)
	{
		//delete unusing histograms in order to prevent resource consuming
		if(iiter != ibest) delete hMC[iiter];
	}
//}}}

	FILE* ftxt;
	ftxt = fopen(Form("Scaled/ScaleFactor_%s_Range%d_%s_by_%s.txt", VarName[ivar].Data(), (int) RangeCut, akaM.Data(), akaD.Data()), "w");
	fprintf(ftxt, "%f", ScaleX[ibest]);
	fclose(ftxt);

//Draw comp{{{
	TCanvas* c1 = new TCanvas("c1", "", 0, 0, 600, 600);
	c1->cd();
	c1->SetLogy();
	FormTH1(href, 0);
	href->Draw("hist");
	hMC[ibest]->Draw("samepe");
	c1->SaveAs(Form("Scaled/Scaled_MC_Data_Comp_%s_%s_by_%s.pdf", VarName[ivar].Data(), akaM.Data(), akaD.Data()));
//}}}

//Draw chi2{{{
	TCanvas* c2 = new TCanvas("c2", "", 0, 0, 600, 600);
	c2->cd();
	g1->SetTitle(";Scale factor;Reduced #chi^2 (normalized by ndf)");
	g1->GetXaxis()->CenterTitle();
	g1->GetYaxis()->CenterTitle();
	g1->SetMarkerStyle(20);
	g1->Draw("ap");
	c2->SaveAs(Form("Scaled/Scale_factor_chi2_scan_%s_%s_by_%s.pdf", VarName[ivar].Data(), akaM.Data(), akaD.Data()));
//}}}

//Draw ratio{{{
	TCanvas* c3 = new TCanvas("c3", "", 0, 0, 600, 600);
	c3->cd();
	TH1D* hratio = (TH1D*) hMC[ibest]->Clone(Form("hratio_%s", VarName[ivar].Data()));
	hratio->SetAxisRange(0, 1.2, "Y");
	hratio->Divide(href);
	hratio->Draw("pe");
	SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
	c3->SaveAs(Form("Scaled/Scaled_ratio_%s_%s_by_%s.pdf", VarName[ivar].Data(), akaM.Data(), akaD.Data()));
//}}}

	hMC[ibest]->SetName(Form("h%s", VarName[ivar].Data()));

//Save plots{{{
	TFile* fout = new TFile(Form("Scaled/Scaled_MC_eff_%s_%s_by_%s.root", VarName[ivar].Data(), akaM.Data(), akaD.Data()), "RECREATE");
	fout->cd();
	g1->SetName(Form("Scalefactor_%s", VarName[ivar].Data()));
	g1->Write();
	href->SetName(Form("h%s_ref", VarName[ivar].Data()));
	href->Write();
	hMC[ibest]->SetName(Form("h%s_scaled", VarName[ivar].Data()));
	hMC[ibest]->Write();
	hratio->Write();
//}}}
}
