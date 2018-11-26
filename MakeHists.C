//Headers{{{
#include "Utilities/Style_Header.h"
#include "Utilities/Var_Header.h"
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
//}}}

void MakeHists()
{
//make directory{{{
	SetStyle();
	TString mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	TString saveDIR = mainDIR + "/HistFiles";
	void * dirp = gSystem->OpenDirectory(saveDIR.Data());
	if(dirp) gSystem->FreeDirectory(dirp);
	else gSystem->mkdir(saveDIR.Data(), kTRUE);
//}}}

	TString fname1 = "/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/PromptRecoForests/HIMinimumBias0/v2/000/326/622/HiForest_326622.root";
	//TString fname = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/HiForest_AMPT_PbPb_5020GeV_stringMelting_v10.root";//AMPT String melt
	//TString fname = "/afs/cern.ch/work/h/hckim/public/Run2018/CentMC_HiForest/add_AMPT_nostringmelting_HiForestAOD_v10KNU.root";//AMPT NoString melt
	//TString fname = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/HiForestAOD_HydjetCymbal5F_PbPb_MC_5TeV_CMSSW10_3_1_event_71000.root";//Hydjet Cymbal5F
	//TString fname = "/eos/cms/store/group/phys_heavyions/ygo/PbPb2018/MC/EPOS/EposLHC_PbPb_5TeV_v10_1031.root";//EPOS
	TString aka = "PR326622";
	//TString aka = "EPOS";
	Bool_t isColl = true;//collision: true, noise: false

//get data{{{
	TChain* t_evt = new TChain("hiEvtAnalyzer/HiTree");
	TChain* t_skim = new TChain("skimanalysis/HltTree");
	TChain* t_hlt = new TChain("hltanalysis/HltTree");
	TChain* t_trk = new TChain("ppTrack/trackTree");

	t_evt->Add(Form("%s", fname1.Data()));
	t_skim->Add(Form("%s", fname1.Data()));
	t_hlt->Add(Form("%s", fname1.Data()));
	t_trk->Add(Form("%s", fname1.Data()));
	t_evt->AddFriend(t_skim);
	t_evt->AddFriend(t_hlt);
	t_evt->AddFriend(t_trk);
//}}}

	TFile* fout = new TFile(Form("HistFiles/PbPb2018_%s_histo.root", aka.Data()), "RECREATE");
	fout->cd();

	if(isColl)
	{
//collision case{{{
		for(Int_t ivar = 0; ivar < NVar; ivar++)
		{
			cout << VarName[ivar].Data() << " loop" << endl;

//make new bin{{{
			Double_t tmpBinArr[nBinC[ivar]];
			tmpBinArr[0] = 0;
			Int_t NewnBin = 0;
			for(Int_t ibin = 1; ibin < nBinC[ivar]+2; ibin++)
			{
				if(tmpBinArr[ibin-1] < NormRangeCut[ivar]) tmpBinArr[ibin] = tmpBinArr[ibin-1]+( VarMaxC[ivar]/(double)nBinC[ivar] );
				else tmpBinArr[ibin] = tmpBinArr[ibin-1]+( VarW[ivar]*VarMaxC[ivar]/(double)nBinC[ivar] );
//cout << ibin << " th: " << tmpBinArr[ibin] << endl;
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
				cout << "variable: " << ivar << ", selection: " << isel << endl;
				TH1D* h = new TH1D(Form("h%s_%s", VarName[ivar].Data(), SelName[isel].Data()), Form("%s;%s;", SelName[isel].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
				t_evt->Draw(Form("%s>>h%s_%s", Vars[ivar].Data(), VarName[ivar].Data(), SelName[isel].Data()), Form("%s", Selection[isel].Data()));
				h->Write();
			}
		}
//}}}
	}
	else
	{
//noise case{{{
		for(Int_t ivar = 0; ivar < NVar; ivar++)
		{
			for(Int_t isel = 0; isel < NSel; isel++)
			{
				cout << "variable: " << ivar << ", selection: " << isel << endl;
				TH1D* h = new TH1D(Form("h%s_%s", VarName[ivar].Data(), SelName[isel].Data()), Form("%s;%s;", SelName[isel].Data(), VarName[ivar].Data()), nBinN[ivar], 0, VarMaxN[ivar]);
				t_evt->Draw(Form("%s>>h%s_%s", Vars[ivar].Data(), VarName[ivar].Data(), SelName[isel].Data()), Form("%s", Selection[isel].Data()));
				h->Write();
			}
		}
//}}}
	}
	fout->Close();
}
