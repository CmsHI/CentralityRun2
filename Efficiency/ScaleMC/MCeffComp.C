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
	const Int_t NCut = 7;
	TString CutName[NCut-1] = {"PV", "CC", "Coin2Th4", "Coin3Th3",
									"PV+CC+Coin2Th4", "PV+CC+Coin3Th3"};
	TString CutLeg[NCut-1] = {"PV", "CC", "Coin2th4", "Coin3th3",
									"PV_CC_Coin2th4", "PV_CC_Coin3th3"};

//get hist{{{
	TH1D* h1[NOF][NVar][NCut];
	Double_t n1[NOF][NVar][NCut];
	TEfficiency* effH[NOF][NVar][NCut-1];

	TFile* fout = new TFile("MCefficiency/PureMCeff.root", "RECREATE");

	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		fin[ifile] = new TFile(Form("../../HistFiles/PbPb2018_%s_histo_newCNT.root", aka[ifile].Data()));

		for(Int_t ivar = 0; ivar < NVar; ivar++)
		{
			for(Int_t icut = 0; icut < NCut; icut++)
			{
				if(icut == 0) h1[ifile][ivar][icut] = (TH1D*) fin[ifile]->Get(Form("h%s_", VarName[ivar].Data()));
				else h1[ifile][ivar][icut] = (TH1D*) fin[ifile]->Get(Form("h%s_%s", VarName[ivar].Data(), CutLeg[icut-1].Data()));
				FormTH1Marker(h1[ifile][ivar][icut], icut, icut, 1.4);
				n1[ifile][ivar][icut] = h1[ifile][ivar][icut]->GetEntries();
				if(icut > 0)
				{
					effH[ifile][ivar][icut-1] = new TEfficiency(*h1[ifile][ivar][icut], *h1[ifile][ivar][0]);
					if(icut == 5) effH[ifile][ivar][icut-1]->SetName(Form("heff_%s_%s_coin2th4", aka[ifile].Data(), VarName[ivar].Data()));
					else if(icut == 6) effH[ifile][ivar][icut-1]->SetName(Form("heff_%s_%s_coin3th3", aka[ifile].Data(), VarName[ivar].Data()));
					h1[ifile][ivar][icut]->Divide(h1[ifile][ivar][icut], h1[ifile][ivar][0], 1, 1, "b");
				}
			}
			fout->cd();
			effH[ifile][ivar][4]->Write();
			effH[ifile][ivar][5]->Write();
		}
	}
	fout->Close();
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

//Cut Comp{{{

//total{{{
		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			TCanvas* cCutComp = new TCanvas("cCutComp", "", 0, 0, 600, 600);
			cCutComp->cd();
			TLegend* legCutComp = new TLegend(0.5, 0.5, 0.9, 0.78);
			FormLegend(legCutComp, 0.04);
			TH1D* hCutComptmp = new TH1D("hCutComptmp", Form("%s;%s;", aka[ifile].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
			hCutComptmp->SetAxisRange(0, 1.2, "Y");
			hCutComptmp->Draw();
			for(Int_t icut = 1; icut < NCut; icut++)
			{
				h1[ifile][ivar][icut]->Draw("samepe");
				legCutComp->AddEntry(h1[ifile][ivar][icut], Form("%s", CutLeg[icut-1].Data()), "p");
				lt1->DrawLatex(0.25, 0.45 - ( (double) icut)*0.05, Form("%s eff. = %.3f %%", CutName[icut-1].Data(), 100*( n1[ifile][ivar][icut]/n1[ifile][ivar][0] )));
			}
			legCutComp->Draw();
			cCutComp->SaveAs(Form("MCefficiency/MCCompCut_%s_%s.pdf", VarName[ivar].Data(), aka[ifile].Data()));
		}
//}}}

//zoom{{{
		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			TCanvas* cCutCompZoom = new TCanvas("cCutCompZoom", "", 0, 0, 600, 600);
			cCutCompZoom->cd();
			TLegend* legCutCompZoom = new TLegend(0.5, 0.5, 0.9, 0.78);
			FormLegend(legCutCompZoom, 0.04);
			TH1D* hCutComptmpZoom = new TH1D("hCutComptmpZoom", Form("%s;%s;", aka[ifile].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
			hCutComptmpZoom->SetAxisRange(0, 1.2, "Y");
			hCutComptmpZoom->SetAxisRange(0, VarMaxZ[ivar], "X");
			hCutComptmpZoom->Draw();
			for(Int_t icut = 1; icut < NCut; icut++)
			{
				h1[ifile][ivar][icut]->Draw("samepe");
				legCutCompZoom->AddEntry(h1[ifile][ivar][icut], Form("%s", CutLeg[icut-1].Data()), "p");
				lt1->DrawLatex(0.25, 0.45 - ( (double) icut)*0.05, Form("%s eff. = %.3f %%", CutName[icut-1].Data(), 100*( n1[ifile][ivar][icut]/n1[ifile][ivar][0] )));
			}
			legCutCompZoom->Draw();
			cCutCompZoom->SaveAs(Form("MCefficiency/MCCompCut_%s_%s_zoom.pdf", VarName[ivar].Data(), aka[ifile].Data()));
		}
//}}}

//}}}

//MC Comp{{{

//Coin2Th4{{{

//total{{{
		TCanvas* c1 = new TCanvas("c1", "", 0, 0, 600, 600);
		c1->cd();
		TH1D* htmp1 = new TH1D("htmp1", Form(";%s;", VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp1->SetAxisRange(0, 1.2, "Y");
		TLegend* leg1 = new TLegend(0.5, 0.55, 0.9, 0.78);
		FormLegend(leg1, 0.04);
		htmp1->Draw();

		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			FormTH1Marker(h1[ifile][ivar][5], ifile, ifile, 1.4);
			h1[ifile][ivar][5]->Draw("samepe");
			leg1->AddEntry(h1[ifile][ivar][5], Form("%s", aka[ifile].Data()), "p");
			lt1->DrawLatex(0.25, 0.45 - ( (double) ifile)*0.05, Form("%s eff. = %.3f %%", aka[ifile].Data(), 100*( n1[ifile][ivar][5]/n1[ifile][ivar][0] )));
		}
		leg1->Draw();
		SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
		c1->SaveAs(Form("MCefficiency/MCCompEff_dist_%s_Coin2th4_newCNT.pdf", VarName[ivar].Data()));
//}}}

//zoom{{{
		TCanvas* c1Zoom = new TCanvas("c1Zoom", "", 0, 0, 600, 600);
		c1Zoom->cd();
		TLegend* leg1Zoom = new TLegend(0.5, 0.55, 0.9, 0.78);
		FormLegend(leg1Zoom, 0.04);
		TH1D* htmp1Zoom = new TH1D("htmp1Zoom", Form(";%s;", VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp1Zoom->SetAxisRange(0, VarMaxZ[ivar], "X");
		htmp1Zoom->SetAxisRange(0, 1.2, "Y");
		htmp1Zoom->Draw();

		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			h1[ifile][ivar][5]->Draw("samepe");
			leg1Zoom->AddEntry(h1[ifile][ivar][5], Form("%s", aka[ifile].Data()), "p");
			lt1->DrawLatex(0.25, 0.45 - ( (double) ifile)*0.05, Form("%s eff. = %.3f %%", aka[ifile].Data(), 100*( n1[ifile][ivar][5]/n1[ifile][ivar][0] )));
		}
		leg1Zoom->Draw();
		SetLine(1, 0, 1, VarMaxZ[ivar], 1, 0, 2);
		c1Zoom->SaveAs(Form("MCefficiency/MCCompEff_dist_%s_Coin2th4_newCNT_zoom.pdf", VarName[ivar].Data()));
//}}}

//}}}

//Coin3Th3{{{

//total{{{
		TCanvas* c2 = new TCanvas("c2", "", 0, 0, 600, 600);
		c2->cd();
		TH1D* htmp2 = new TH1D("htmp2", Form(";%s;", VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp2->SetAxisRange(0, 1.2, "Y");
		TLegend* leg2 = new TLegend(0.5, 0.55, 0.9, 0.78);
		FormLegend(leg2, 0.04);
		htmp2->Draw();

		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			FormTH1Marker(h1[ifile][ivar][6], ifile, ifile, 1.4);
			h1[ifile][ivar][6]->Draw("samepe");
			leg2->AddEntry(h1[ifile][ivar][6], Form("%s", aka[ifile].Data()), "p");
			lt1->DrawLatex(0.25, 0.45 - ( (double) ifile)*0.05, Form("%s eff. = %.3f %%", aka[ifile].Data(), 100*( n1[ifile][ivar][6]/n1[ifile][ivar][0] )));
		}
		leg2->Draw();
		SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
		c2->SaveAs(Form("MCefficiency/MCCompEff_dist_%s_Coin3th3_newCNT.pdf", VarName[ivar].Data()));
//}}}

//zoom{{{
		TCanvas* c2Zoom = new TCanvas("c2Zoom", "", 0, 0, 600, 600);
		c2Zoom->cd();
		TLegend* leg2Zoom = new TLegend(0.5, 0.55, 0.9, 0.78);
		FormLegend(leg2Zoom, 0.04);
		TH1D* htmp2Zoom = new TH1D("htmp2Zoom", Form(";%s;", VarName[ivar].Data()), NewnBin-1, NewBinArr);
		htmp2Zoom->SetAxisRange(0, VarMaxZ[ivar], "X");
		htmp2Zoom->SetAxisRange(0, 1.2, "Y");
		htmp2Zoom->Draw();

		for(Int_t ifile = 0; ifile < NOF; ifile++)
		{
			h1[ifile][ivar][6]->Draw("samepe");
			leg2Zoom->AddEntry(h1[ifile][ivar][6], Form("%s", aka[ifile].Data()), "p");
			lt1->DrawLatex(0.25, 0.45 - ( (double) ifile)*0.05, Form("%s eff. = %.3f %%", aka[ifile].Data(), 100*( n1[ifile][ivar][6]/n1[ifile][ivar][0] )));
		}
		leg2Zoom->Draw();
		SetLine(1, 0, 1, VarMaxZ[ivar], 1, 0, 2);
		c2Zoom->SaveAs(Form("MCefficiency/MCCompEff_dist_%s_Coin3th3_newCNT_zoom.pdf", VarName[ivar].Data()));
//}}}

//}}}

//}}}
	}
}
