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

void CompScan()
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
	TString akaD = "PR326587";
	TString akaM[NOF] = {"Hydjet_Cymbal5F", "AMPT_String", "AMPT_NoString",
								"EPOS"};

	Double_t RangeCut[3] = {80, 300, 500};
	TFile* fin[NOF];
	TCanvas* c1[NOF];
	TGraphErrors* g1[NOF][3];

//Get chi2 dist{{{
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		fin[ifile] = new TFile(Form("Scaled/Scaled_MC_eff_sel_%s_%s_%s.root", VarName[0].Data(), akaD.Data(), akaM[ifile].Data()), "READ");
		for(Int_t icut = 0; icut < 3; icut++)
		{
			g1[ifile][icut] = (TGraphErrors*) fin[ifile]->Get(Form("Scalefactor_%s_%d", VarName[0].Data(), (int)RangeCut[icut]));
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
		Double_t minchi2[3];
		ifstream in;
		in.open(Form("Scaled/ScaleFactor_sel_%s_%s_%s.txt", VarName[0].Data(), akaD.Data(), akaM[ifile].Data()));
		if(in.is_open())
		{
			for(Int_t icut = 0; icut < 3; icut++)
			{
				if(icut == 0) g1[ifile][icut]->Draw("ap");
				else g1[ifile][icut]->Draw("samep");
				leg->AddEntry(g1[ifile][icut], Form("Range: %d ~ #infty", (int)RangeCut[icut]), "p");
				in >> minchi2[icut];
				SetLine(1, minchi2[icut], 0, minchi2[icut], 100, icut, 2);
			}
			leg->Draw();
			TLatex* ltx1 = new TLatex();
			FormLatex(ltx1, 12, 0.06);
			ltx1->DrawLatex(0.5, 0.6, Form("%s", akaM[ifile].Data()));
			c1[ifile]->SaveAs(Form("Scaled/chi2_cut_comp_sel_%s_%s.pdf", akaD.Data(), akaM[ifile].Data()));
		}
	}
//}}}
}
