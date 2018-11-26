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
TH1D* ScaleTH1(TTree* t1, Int_t iiter, Double_t ScaleX, Double_t RangeCut)
{
//scale histogram{{{

//make new bin{{{
	Double_t tmpBinArr[nBinC[0]];
	tmpBinArr[0] = 0;
	Int_t NewnBin = 0;
	for(Int_t ibin = 1; ibin < nBinC[0]+2; ibin++)
	{
		if(tmpBinArr[ibin-1] < NormRangeCut[0]) tmpBinArr[ibin] = tmpBinArr[ibin-1]+( VarMaxC[0]/(double)nBinC[0] );
		else tmpBinArr[ibin] = tmpBinArr[ibin-1]+( VarW[0]*VarMaxC[0]/(double)nBinC[0] );
		if(tmpBinArr[ibin] >= VarMaxC[0])
		{
			NewnBin = ibin+1;
			break;
		}
	}
	Double_t NewBinArr[NewnBin];
	NewBinArr[0] = 0;
	for(Int_t ibin = 1; ibin < NewnBin; ibin++)
	{
		if(NewBinArr[ibin-1] < NormRangeCut[0]) NewBinArr[ibin] = NewBinArr[ibin-1]+( VarMaxC[0]/(double)nBinC[0] );
		else NewBinArr[ibin] = NewBinArr[ibin-1]+( VarW[0]*VarMaxC[0]/(double)nBinC[0] );
	}
//}}}

	TH1D* h1 = new TH1D(Form("h1_%d", iiter), "", NewnBin-1, NewBinArr);
	t1->Draw(Form("%f*%s>>h1_%d", ScaleX, Vars[0].Data(), iiter), "");
	FormTH1(h1, 1);
	h1->Scale(1./( h1->Integral(h1->GetXaxis()->FindBin(RangeCut), h1->GetXaxis()->FindBin(VarMaxC[0])) ));
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

void CutScanForScaleFactor()
{
	SetStyle();

//make directory{{{
	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	string saveDIR = mainDIR + "/Scaled";
	void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
	if (dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.c_str(), kTRUE);
//}}}

	TString akaD = "PR326587";
	TString akaM = "Hydjet_Cymbal5F";
	TFile* fdata = new TFile(Form("../../HistFiles/PbPb2018_%s_histo.root", akaD.Data()), "READ");
	TString fMC = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/HiForestAOD_HydjetCymbal5F_PbPb_MC_5TeV_CMSSW10_3_1_event_71000.root";//Hydjet Cymbal5F
	//TString fMC = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/HiForest_AMPT_PbPb_5020GeV_stringMelting_v10.root";//AMPT String Melting
	//TString fMC = "/afs/cern.ch/work/h/hckim/public/Run2018/CentMC_HiForest/add_AMPT_nostringmelting_HiForestAOD_v10KNU.root";//AMPT No String Melting
	//TString fMC = "/eos/cms/store/group/phys_heavyions/ygo/PbPb2018/MC/EPOS/EposLHC_PbPb_5TeV_v10_1031.root";//EPOS

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

	Double_t RangeCut[3] = {80, 300, 500};
	const Int_t Niter = 60;
	FILE* ftxt;
	ftxt = fopen(Form("Scaled/ScaleFactor_sel_%s_%s_%s.txt", VarName[0].Data(), akaD.Data(), akaM.Data()), "w");
	TFile* fout = new TFile(Form("Scaled/Scaled_MC_eff_sel_%s_%s_%s.root", VarName[0].Data(), akaD.Data(), akaM.Data()), "RECREATE");
	fout->cd();

	for(Int_t icut = 0; icut < 3; icut++)
	{
		TH1D* href = (TH1D*) fdata->Get(Form("h%s_PV_CC_Coin2th4", VarName[0].Data()));
		FormTH1(href, 0);
		href->Scale(1./href->Integral());
		href->Scale(1., "width");

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
			hMC[iiter] = ScaleTH1(t_evt, iiter, ScaleX[iiter], RangeCut[icut]);
			chi2val[iiter] = Getchi2(hMC[iiter], href, RangeCut[icut]);
			g1->SetPoint(iiter, ScaleX[iiter], chi2val[iiter]);
			if(chi2val[iiter] < chi2val[ibest]) ibest = iiter;
		}

		for(Int_t iiter = 0; iiter < Niter; iiter++)
		{
			//delete unusing histograms in order to prevent resource consuming
			if(iiter != ibest) delete hMC[iiter];
		}
//}}}

		fprintf(ftxt, "%f \n", ScaleX[ibest]);

//Draw comp{{{
		TCanvas* c1 = new TCanvas("c1", "", 0, 0, 600, 600);
		c1->cd();
		c1->SetLogy();
		FormTH1(href, 0);
		href->Draw("hist");
		hMC[ibest]->Draw("samepe");
		c1->SaveAs(Form("Scaled/Scaled_MC_Data_Comp_sel_%s_cut%d_%s_%s.pdf", VarName[0].Data(), (int)RangeCut[icut], akaD.Data(), akaM.Data()));
//}}}

//Draw chi2{{{
		TCanvas* c2 = new TCanvas("c2", "", 0, 0, 600, 600);
		c2->cd();
		g1->SetTitle(";Scale factor;Reduced #chi^2 (normalized by ndf)");
		g1->GetXaxis()->CenterTitle();
		g1->GetYaxis()->CenterTitle();
		g1->SetMarkerStyle(20);
		g1->Draw("ap");
		c2->SaveAs(Form("Scaled/Scale_factor_chi2_scan_sel_%s_cut%d_%s_%s.pdf", VarName[0].Data(), (int)RangeCut[icut], akaD.Data(), akaM.Data()));
//}}}

//Draw ratio{{{
		TCanvas* c3 = new TCanvas("c3", "", 0, 0, 600, 600);
		c3->cd();
		TH1D* hratio = (TH1D*) hMC[ibest]->Clone(Form("hratio_%s", VarName[0].Data()));
		hratio->SetAxisRange(0, 1.2, "Y");
		hratio->Divide(href);
		hratio->Draw("pe");
		SetLine(1, 0, 1, VarMaxC[0], 1, 0, 2);
		c3->SaveAs(Form("Scaled/Scaled_ratio_sel_%s_cut%d_%s_%s.pdf", VarName[0].Data(), (int)RangeCut[icut], akaD.Data(), akaM.Data()));
//}}}

//Save plots{{{
		g1->SetName(Form("Scalefactor_%s_%d", VarName[0].Data(), (int)RangeCut[icut]));
		g1->Write();
		href->SetName(Form("h%s_ref%d", VarName[0].Data(), (int)RangeCut[icut]));
		href->Write();
		hMC[ibest]->SetName(Form("h%s_scaled%d", VarName[0].Data(), (int)RangeCut[icut]));
		hMC[ibest]->Write();
		hratio->Write();
//}}}
	}
	fclose(ftxt);
}
