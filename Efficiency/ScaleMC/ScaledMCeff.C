//Headers{{{
#include "../../Utilities/Style_Header.h"
#include "../../Utilities/Var_Header.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TCut.h>
#include <TMath.h>
#include <iostream>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>
#include <TParameter.h>
//}}}

void ScaledMCeff(const Int_t ivar = 5, const Int_t Coin = 2, const Int_t Th = 4, const Double_t RangeCut = 200)
{
	SetStyle();
	TString akaD = "PR326478";
	TString akaM = "EPOS";

//make directory{{{
	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	string saveDIR = mainDIR + "/MCefficiency/";
	void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
	if (dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.c_str(), kTRUE);
//}}}

	//TString filename = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/HydjetCymbal5F_Forest/HYDJET_CYMBAL5F_PbPb_5020GeV/HydjetCymbal5F_5020GeV_PbPb_Forest/181128_165108/HydjetCymbal5F_5020GeV_PbPb_Forest.root";//Hydjet Cymbal5F
	//TString filename = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/AMPT_StringMelting_Forest/AMPT_StringMelting_Forest.root";//AMPT String Melting
	//TString filename = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/AMPT_noStringMelting_Forest/AMPT_No_StringMelting_Forest.root";//AMPT No String Melting
	TString filename = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/EPOS_Forest/EPOS_Forest.root";//EPOS

//get tree{{{
	TChain* t_evt = new TChain("hiEvtAnalyzer/HiTree");
	TChain* t_skim = new TChain("skimanalysis/HltTree");
	TChain* t_hlt = new TChain("hltanalysis/HltTree");
	TChain* t_trk = new TChain("ppTrack/trackTree");
	TChain* t_rec = new TChain("rechitanalyzerpp/tower");
	t_evt->Add(Form("%s", filename.Data()));
	t_skim->Add(Form("%s", filename.Data()));
	t_hlt->Add(Form("%s", filename.Data()));
	t_trk->Add(Form("%s", filename.Data()));
	t_rec->Add(Form("%s", filename.Data()));
	t_evt->AddFriend(t_skim);
	t_evt->AddFriend(t_hlt);
	t_evt->AddFriend(t_trk);
	t_evt->AddFriend(t_rec);
//}}}

//get variables{{{
	Float_t hiHF, hiHFplus, hiHFplusEta4, hiHFminus, hiHFminusEta4,
			hiHFECut, hiHFECutPlus, hiHFECutMinus, hiHFhit, hiEB, hiEE;
	Int_t MAXTOWSIZE = 10000;
	Float_t e[MAXTOWSIZE];
	Float_t et[MAXTOWSIZE];
	Float_t eta[MAXTOWSIZE];
	Int_t hiNtracks, hiNpix, n;
	Int_t pprimaryVertexFilter;
	Int_t pclusterCompatibilityFilter;
	t_evt->SetBranchAddress("hiHF",&hiHF);
	t_evt->SetBranchAddress("hiHFplus",&hiHFplus);
	t_evt->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
	t_evt->SetBranchAddress("hiHFminus",&hiHFminus);
	t_evt->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
	t_evt->SetBranchAddress("hiHFECut",&hiHFECut);
	t_evt->SetBranchAddress("hiHFECutPlus",&hiHFECutPlus);
	t_evt->SetBranchAddress("hiHFECutMinus",&hiHFECutMinus);
	t_evt->SetBranchAddress("hiEB",&hiEB);
	t_evt->SetBranchAddress("hiEE",&hiEE);
	t_evt->SetBranchAddress("hiNtracks",&hiNtracks);
	t_evt->SetBranchAddress("hiHFhit",&hiHFhit);
	t_evt->SetBranchAddress("hiNpix",&hiNpix);
	t_evt->SetBranchAddress("eta",&eta);
	t_evt->SetBranchAddress("e",&e);
	t_evt->SetBranchAddress("et",&et);
	t_evt->SetBranchAddress("n",&n);
	t_evt->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
	t_evt->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
//}}}

//define histogram{{{

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

	TH1D* htot = new TH1D(Form("htot_%s", VarName[ivar].Data()), Form("%s;%s;entries", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
	TH1D* hsel = new TH1D(Form("hsel_%s", VarName[ivar].Data()), Form("%s;%s;entries", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);

	FormTH1Marker(htot, 0, 0, 1.4);
	FormTH1Marker(hsel, 1, 0, 1.4);
//}}}

//Get scale factor{{{
	Double_t towSF = 1.;//scale factor for tower energy
	Double_t SF = 1.;//scale factor for variable
/*
	ifstream in1;
	in1.open(Form("Scaled/ScaleFactor_HF_%s_by_%s.txt", akaM.Data(), akaD.Data()));
	if(in1.is_open())
	{
		in1 >> towSF;
	}
	else
	{
		cout << "There is No such file!!! Please confirm the name!!!" << endl;
		return;
	}
	in1.close();
*/
	ifstream in2;
	in2.open(Form("Scaled/ScaleFactor_%s_Range%d_%s_by_%s.txt", VarName[ivar].Data(), (int) RangeCut, akaM.Data(), akaD.Data()));
	if(in2.is_open())
	{
		in2 >> SF;
	}
	else
	{
		cout << "There is No such file!!! Please confirm the name!!!" << endl;
		return;
	}
	in2.close();
//}}}

	Long64_t nEntries = t_evt->GetEntries();
	//t_evt->IncrementTotalBuffers(1937220430);
	//nEntries = 40000;
	//nEntries = 1000;
	for(Long64_t jentry=0; jentry<nEntries;jentry++)
	{
		t_evt->GetEntry(jentry);

		if( jentry % 10000 == 0 ) std::cout << jentry << "/" << nEntries << std::endl;

//variable to fill{{{
		Double_t currVar = 0.;
		if(ivar == 0) currVar = hiHF;
		if(ivar == 1) currVar = hiHFplus;
		if(ivar == 2) currVar = hiHFplusEta4;
		if(ivar == 3) currVar = hiHFminus;
		if(ivar == 4) currVar = hiHFminusEta4;
		if(ivar == 5) currVar = hiHFECut;
		if(ivar == 6) currVar = hiHFECutPlus;
		if(ivar == 7) currVar = hiHFECutMinus;
		if(ivar == 8) currVar = hiEB;
		if(ivar == 9) currVar = hiEE;
		if(ivar == 10) currVar = hiNtracks;
		if(ivar == 11) currVar = hiHFhit;
		if(ivar == 12) currVar = hiNpix;
//}}}

//Fill histogram{{{
		htot->Fill(SF*currVar);

		if(!pprimaryVertexFilter) continue;
		if(!pclusterCompatibilityFilter) continue;

		Int_t nplus = 0;
		Int_t nminus = 0;
		for(Int_t in = 0; in < n; in++)
		{
			if(abs(eta[in]) < 3.0 || abs(eta[in]) > 6.0) continue;
			if(eta[in] >= 3.0 && eta[in] <= 6.0 && towSF*e[in] >= Th) nplus++;
			if(eta[in] <= -3.0 && eta[in] >= -6.0 && towSF*e[in] >= Th) nminus++;
		}
		if(nplus < Coin && nminus < Coin) continue;
		hsel->Fill(SF*currVar);
//}}}
	}

	TH1D* heff = (TH1D*) hsel->Clone(Form("heff_%s", VarName[ivar].Data()));
	heff->Divide(heff, htot, 1, 1, "b");

	TParameter<double>* effVal = new TParameter<double>(Form("v_eff_%s", VarName[ivar].Data()), ((double) hsel->Integral("width")/(double) htot->Integral("width")));
	TEfficiency* effH = new TEfficiency(*hsel, *htot);
	effH->SetName(Form("heff_%s_coin%dth%d", VarName[ivar].Data(), Coin, Th));
	TGraphAsymmErrors* geff = effH->CreateGraph();
	geff->SetName(Form("geff_%s_coin%dth%d", VarName[ivar].Data(), Coin, Th));
	TH1D* htmp = new TH1D("htmp", Form(";%s;efficiency", VarName[ivar].Data()), NewnBin-1, NewBinArr);
	htmp->GetYaxis()->SetRangeUser(0, 1.2);

	TLatex* lt1 = new TLatex();
	FormLatex(lt1, 12, 0.04);
	TCanvas* c1 = new TCanvas("c1", "", 0, 0, 600, 600);
	c1->cd();
	htmp->Draw();
	geff->Draw("same");
	SetLine(1, 0, 1, VarMaxC[ivar], 1, 0, 2);
	lt1->DrawLatex(0.5, 0.6, Form("eff. = %.3f %%", 100*effVal->GetVal()));
	lt1->DrawLatex(0.5, 0.5, Form("%s", akaM.Data()));
	c1->SaveAs(Form("MCefficiency/eff_dist_%s_Range%d_coin%dth%d_%s_by_%s.pdf", VarName[ivar].Data(), (int)RangeCut, Coin, Th, akaM.Data(), akaD.Data()));

	TFile* fout = new TFile(Form("MCefficiency/MC_eff_2018_%s_Range%d_coin%dth%d_%s_by_%s.root", VarName[ivar].Data(), (int)RangeCut, Coin, Th, akaM.Data(), akaD.Data()), "RECREATE");
	fout->cd();
	hsel->Write();
	heff->Write();
	effVal->Write();
	effH->Write();
	geff->Write();
	fout->Close();
}
