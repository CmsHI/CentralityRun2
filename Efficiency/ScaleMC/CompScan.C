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

void CompScan(const Int_t ivar = 0)
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
	TString akaD = "PR326483";
	TString akaM[NOF] = {"Hydjet_Cymbal5F", "AMPT_String", "AMPT_NoString",
								"EPOS"};

	const Int_t NOC = 6;//Number Of Cut
	const Double_t RangeCut[NOC] = {80, 100, 200, 300, 500, 1000};
	TFile* fin[NOF];
	TCanvas* c1[NOF];
	TGraphErrors* g1[NOF][NOC];

//Get chi2 dist{{{
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		fin[ifile] = new TFile(Form("Scaled/Scaled_MC_eff_scan_%s_%s_%s.root", VarName[ivar].Data(), akaD.Data(), akaM[ifile].Data()), "READ");
		for(Int_t icut = 0; icut < NOC; icut++)
		{
			g1[ifile][icut] = (TGraphErrors*) fin[ifile]->Get(Form("Scalefactor_%s_%d", VarName[ivar].Data(), (int)RangeCut[icut]));
			g1[ifile][icut]->SetMarkerColor(colorArr[icut]);
		}
	}
//}}}

//Draw chi2 dist{{{
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		c1[ifile] = new TCanvas(Form("c1_%d", ifile), "", 0, 0, 600, 600);
		c1[ifile]->cd();
		TLegend* leg = new TLegend(0.5, 0.7, 0.9, 0.9);
		FormLegend(leg, 0.04);
		TString word[NOC];
		Double_t minchi2[NOC];
		Double_t chi2val[NOC];
		ifstream in;
		in.open(Form("Scaled/ScaleFactor_scan_%s_%s_%s.txt", VarName[ivar].Data(), akaD.Data(), akaM[ifile].Data()));
		if(in.is_open())
		{
			for(Int_t icut = 0; icut < NOC; icut++)
			{
				if(icut == 0) g1[ifile][icut]->Draw("ap");
				else g1[ifile][icut]->Draw("samep");
				leg->AddEntry(g1[ifile][icut], Form("Range: %d ~", (int)RangeCut[icut]), "p");
				in >> word[icut] >> minchi2[icut] >> chi2val[icut];
				SetLine(1, minchi2[icut], 0, minchi2[icut], 100, icut, 2);
			}
			leg->Draw();
			TLatex* ltx1 = new TLatex();
			FormLatex(ltx1, 12, 0.06);
			ltx1->DrawLatex(0.5, 0.6, Form("%s", akaM[ifile].Data()));
			c1[ifile]->SaveAs(Form("Scaled/chi2_cut_comp_scan_%s_%s_%s.pdf", VarName[ivar].Data(), akaD.Data(), akaM[ifile].Data()));
		}
	}
//}}}
}
