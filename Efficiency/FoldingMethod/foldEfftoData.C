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
#include <TSystem.h>
#include <TParameter.h>
//}}}

void foldEfftoData(const Int_t ivar = 5, const Int_t Coin = 2, const Int_t Th = 4,  const Double_t RangeCut = 200)
{
	SetStyle();
	TString akaM = "EPOS";
	TString akaD = "PR326478";

//make directory{{{
	string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	string saveDIR = mainDIR + "/DATAefficiency/";
	void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
	if (dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.c_str(), kTRUE);
//}}}

//get tree{{{
	const Int_t NOF = 1;//Number Of Files
	TString filename1[NOF] = {"/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/PromptRecoForests/HIMinimumBias0/v5/run326478/181129_165457/0000/HiForestAOD_326478_part1.root"};
	TChain* t_evt = new TChain("hiEvtAnalyzer/HiTree");
	TChain* t_skim = new TChain("skimanalysis/HltTree");
	TChain* t_hlt = new TChain("hltanalysis/HltTree");
	TChain* t_trk = new TChain("ppTrack/trackTree");
	for(Int_t ifile = 0; ifile < NOF; ifile++)
	{
		t_evt->Add(Form("%s", filename1[ifile].Data()));
		t_skim->Add(Form("%s", filename1[ifile].Data()));
		t_hlt->Add(Form("%s", filename1[ifile].Data()));
		t_trk->Add(Form("%s", filename1[ifile].Data()));
	}
	t_evt->AddFriend(t_skim);
	t_evt->AddFriend(t_hlt);
	t_evt->AddFriend(t_trk);
//}}}

	Long64_t nEntries = t_evt->GetEntries();
	//nEntries = 500000;
	//nEntries = 1000;

//get variables{{{
	Float_t hiHF, hiHFplus, hiHFplusEta4, hiHFminus, hiHFminusEta4,
			hiHFECut, hiHFECutPlus, hiHFECutMinus, hiHFhit, hiEE, hiEB;
	Int_t hiNtracks, hiNpix;
	Int_t pprimaryVertexFilter;
	Int_t pclusterCompatibilityFilter;
	Int_t phfCoincFilter3Th3;
	Int_t phfCoincFilter4Th3;
	Int_t phfCoincFilter5Th3;
	Int_t phfCoincFilter2Th4;
	Int_t phfCoincFilter3Th4;
	Int_t phfCoincFilter4Th4;
	t_evt->SetBranchAddress("hiHF",&hiHF);
	t_evt->SetBranchAddress("hiHFplus",&hiHFplus);
	t_evt->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
	t_evt->SetBranchAddress("hiHFminus",&hiHFminus);
	t_evt->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
	t_evt->SetBranchAddress("hiHFECut",&hiHFECut);
	t_evt->SetBranchAddress("hiHFECutPlus",&hiHFECutPlus);
	t_evt->SetBranchAddress("hiHFECutMinus",&hiHFECutMinus);
	t_evt->SetBranchAddress("hiHFhit",&hiHFhit);
	t_evt->SetBranchAddress("hiEB",&hiEB);
	t_evt->SetBranchAddress("hiEE",&hiEE);
	t_evt->SetBranchAddress("hiNtracks",&hiNtracks);
	t_evt->SetBranchAddress("hiNpix",&hiNpix);
	t_evt->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
	t_evt->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
	t_evt->SetBranchAddress("phfCoincFilter3Th3",&phfCoincFilter3Th3);
	t_evt->SetBranchAddress("phfCoincFilter4Th3",&phfCoincFilter4Th3);
	t_evt->SetBranchAddress("phfCoincFilter5Th3",&phfCoincFilter5Th3);
	t_evt->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4);
	t_evt->SetBranchAddress("phfCoincFilter3Th4",&phfCoincFilter3Th4);
	t_evt->SetBranchAddress("phfCoincFilter4Th4",&phfCoincFilter4Th4);
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

	TH1D* h1 = new TH1D(Form("h1_%s", VarName[ivar].Data()), Form("%s;%s;entries", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
	TH1D* h2 = new TH1D(Form("h2_%s", VarName[ivar].Data()), Form("%s;%s;entries", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
	TH1D* h1W = new TH1D(Form("h1W_%s", VarName[ivar].Data()), Form("%s;%s;entries", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
	TH1D* h2W = new TH1D(Form("h2W_%s", VarName[ivar].Data()), Form("%s;%s;entries", VarName[ivar].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);

	FormTH1Marker(h1, 0, 0, 1.4);
	FormTH1Marker(h2, 0, 0, 1.4);
	FormTH1Marker(h1W, 1, 0, 1.4);
	FormTH1Marker(h2W, 1, 0, 1.4);
//}}}

//Get efficiency{{{
	TFile* fEff = new TFile(Form("../ScaleMC/MCefficiency/MC_eff_2018_%s_Range%d_coin%dth%d_%s_by_%s.root", VarName[ivar].Data(), (int)RangeCut, Coin, Th, akaM.Data(), akaD.Data()), "READ");
	TEfficiency* effH = (TEfficiency*)(fEff->Get(Form("heff_%s_coin%dth%d", VarName[ivar].Data(), Coin, Th)));
//}}}

	for(Long64_t jentry=0; jentry<nEntries;jentry++)
	//for(Long64_t jentry=0; jentry<100000;jentry++)
	{
		t_evt->GetEntry(jentry);

		if( jentry % 10000 == 0 ) std::cout << jentry << "/" << nEntries << std::endl;

//apply cut{{{
		if(!pprimaryVertexFilter) continue;
		if(!pclusterCompatibilityFilter) continue;
		if( (Coin == 3 && Th == 3) && !phfCoincFilter3Th3) continue;
		if( (Coin == 4 && Th == 3) && !phfCoincFilter4Th3) continue;
		if( (Coin == 5 && Th == 3) && !phfCoincFilter5Th3) continue;
		if( (Coin == 2 && Th == 4) && !phfCoincFilter2Th4) continue;
		if( (Coin == 3 && Th == 4) && !phfCoincFilter3Th4) continue;
		if( (Coin == 4 && Th == 4) && !phfCoincFilter4Th4) continue;
//}}}

//Fill hist{{{
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
		h1->Fill(currVar);

		Int_t bin = effH->FindFixBin(currVar);
		Double_t effVal = effH->GetEfficiency(bin);
		Int_t cutoff = effH->FindFixBin(NormRangeCut[ivar]);
		Double_t correction = 0.;
		if(bin > cutoff || effVal <= 0) correction = 1.;
		else correction = effVal;
		//if(effVal <= 0) correction = 1.;
		h1W->Fill(currVar, 1./correction);
//}}}
	}

	TFile* fout = new TFile(Form("DATAefficiency/Data_eff_2018_%s_Range%d_coin%dth%d_%s_by_%s.root", VarName[ivar].Data(), (int)RangeCut, Coin, Th, akaM.Data(), akaD.Data()), "RECREATE");

	Double_t DeffVal = ((double) h1->Integral("width")/(double) h1W->Integral("width"));
	cout << "The efficiency from " << VarName[ivar].Data() << " with PV is " << DeffVal << endl;

//Draw corrected data{{{
	TCanvas *cdist1 = new TCanvas(Form("cdist1_%s", VarName[ivar].Data()), "", 0, 0, 600, 600);
	cdist1->cd();
	cdist1->SetLogy();
	h1W->Scale(1., "width");
	h1->Scale(1., "width");
	h1W->Draw("pe");
	h1->Draw("samepe");

	TLegend* leg1 = new TLegend(0.5,0.6,0.90,0.8);
	FormLegend(leg1, 0.04);
	leg1->AddEntry(h1, "Data", "p");
	leg1->AddEntry(h1W, "Data corrected", "p");
	leg1->Draw("same");

	TLatex *lt1 = new TLatex();
	FormLatex(lt1, 12, 0.04);
	lt1->DrawLatex(0.2, 0.5, Form("eff. = %0.3f %%", 100*DeffVal));
	lt1->DrawLatex(0.2, 0.3, Form("%s eff. weights", akaM.Data()));
	cdist1->SaveAs(Form("DATAefficiency/Data_dist_%s_Range%d_PV_coin%d_th%d_%s_by_%s.pdf", VarName[ivar].Data(), (int)RangeCut, Coin, Th, akaM.Data(), akaD.Data()));
//}}}

	TCanvas* ceff1 = new TCanvas(Form("ceff1_%s", VarName[ivar].Data()), "", 0, 0, 600, 600);
	ceff1->cd();
	TH1D* hratio1 = (TH1D*) h1W->Clone(Form("hratio1_%s", VarName[ivar].Data()));
	hratio1->Divide(hratio1, h1, 1, 1, "b");
	hratio1->SetAxisRange(0, 1.2, "Y");
	hratio1->Draw("pe");
	lt1->DrawLatex(0.2, 0.5, Form("eff. = %0.3f %%", 100*DeffVal));
	lt1->DrawLatex(0.2, 0.3, Form("%s eff. weights", akaM.Data()));
	ceff1->SaveAs(Form("DATAefficiency/Data_eff_2018_%s_Range%d_PV_coin%dth%d_%s_by_%s.pdf", VarName[ivar].Data(), (int)RangeCut, Coin, Th, akaM.Data(), akaD.Data()));

	fout->cd();
	h1->Write();
	h1W->Write();
	hratio1->Write();
	fout->Close();
}
