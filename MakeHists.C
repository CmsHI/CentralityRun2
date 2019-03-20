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

/*
	TString fname[4] = {"/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/PromptRecoForests/HIMinimumBias0/v5/run326501/HIMinimumBias0/crab_HIMB0_Forest_run326501_das/181129_165420/0000/HiForestAOD_326501_part1.root",
							"/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/PromptRecoForests/HIMinimumBias0/v5/run326501/HIMinimumBias0/crab_HIMB0_Forest_run326501_das/181129_165420/0000/HiForestAOD_326501_part2.root",
							"/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/PromptRecoForests/HIMinimumBias0/v5/run326501/HIMinimumBias0/crab_HIMB0_Forest_run326501_das/181129_165420/0000/HiForestAOD_326501_part3.root",
							"/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/PromptRecoForests/HIMinimumBias0/v5/run326501/HIMinimumBias0/crab_HIMB0_Forest_run326501_das/181129_165420/0000/HiForestAOD_326501_part4.root",
							};//Data
*/
	//TString fname = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/HydjetCymbal5F_Forest/HYDJET_CYMBAL5F_PbPb_5020GeV/HydjetCymbal5F_5020GeV_PbPb_Forest/181128_165108/HydjetCymbal5F_5020GeV_PbPb_Forest.root";//Hydjet Cymbal5F
	TString fname = "/eos/cms/store/group/phys_heavyions/bdiab/Centrality/HIForest_MC_Hydjet_Drum5F/HiForest_MC_Hydjet_Drum5F_merged.root";//Hydjet Cymbal5F official
	//TString fname = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/AMPT_StringMelting_Forest/AMPT_StringMelting_Forest.root";//AMPT String melt
	//TString fname = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/AMPT_noStringMelting_Forest/AMPT_No_StringMelting_Forest.root";//AMPT NoString melt
	//TString fname = "/eos/cms/store/group/phys_heavyions/dileptons/hanseul/FOREST/EPOS_Forest/EPOS_Forest.root";//EPOS
	//TString fname = "/afs/cern.ch/work/h/hanseul/public/starlight/HiForestAOD_Starlight_prodmode5_eve9999.root";//Starlight Single
	//TString fname = "/afs/cern.ch/work/h/hanseul/public/starlight/HiForestAOD_Starlight_prodmode6_eve9999.root";//Starlight Double
	//TString aka = "PR326501";
	TString aka = "Hydjet_Cymbal5F";
	Bool_t isColl = true;//collision: true, noise: false
	TString Suffix = "official";

//get data{{{
	TChain* t_evt = new TChain("hiEvtAnalyzer/HiTree");
	TChain* t_skim = new TChain("skimanalysis/HltTree");
	//TChain* t_hlt = new TChain("hltanalysis/HltTree");
	//TChain* t_trk = new TChain("ppTrack/trackTree");


	t_evt->Add(Form("%s", fname.Data()));
	t_skim->Add(Form("%s", fname.Data()));
	//t_hlt->Add(Form("%s", fname.Data()));
	//t_trk->Add(Form("%s", fname.Data()));


/*
	for(Int_t i = 0; i < 4; i++)
	{
		t_evt->Add(Form("%s", fname[i].Data()));
		t_skim->Add(Form("%s", fname[i].Data()));
		//t_hlt->Add(Form("%s", fname[i].Data()));
		//t_trk->Add(Form("%s", fname[i].Data()));
	}
*/
	t_evt->AddFriend(t_skim);
	//t_evt->AddFriend(t_hlt);
	//t_evt->AddFriend(t_trk);
//}}}

	TFile* fout = new TFile(Form("HistFiles/PbPb2018_%s_histo_%s.root", aka.Data(), Suffix.Data()), "RECREATE");
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
				if(tmpBinArr[ibin-1] < BinBoundary[ivar]) tmpBinArr[ibin] = tmpBinArr[ibin-1]+( VarMaxC[ivar]/(double)nBinC[ivar] );
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
				if(NewBinArr[ibin-1] < BinBoundary[ivar]) NewBinArr[ibin] = NewBinArr[ibin-1]+( VarMaxC[ivar]/(double)nBinC[ivar] );
				else NewBinArr[ibin] = NewBinArr[ibin-1]+( VarW[ivar]*VarMaxC[ivar]/(double)nBinC[ivar] );
			}
//}}}

			for(Int_t isel = 0; isel < NSel; isel++)
			{
				cout << "variable: " << ivar << ", selection: " << isel << endl;
				TH1D* h = new TH1D(Form("h%s_%s", VarName[ivar].Data(), SelName[isel].Data()), Form("%s;%s;", SelName[isel].Data(), VarName[ivar].Data()), NewnBin-1, NewBinArr);
				TH1D* hCB = new TH1D(Form("h%s_%s_CB", VarName[ivar].Data(), SelName[isel].Data()), Form("%s;%s;", SelName[isel].Data(), VarName[ivar].Data()), nBinC[ivar], 0, VarMaxC[ivar]);
				t_evt->Draw(Form("%s>>h%s_%s", Vars[ivar].Data(), VarName[ivar].Data(), SelName[isel].Data()), Form("%s", Selection[isel].Data()));
				t_evt->Draw(Form("%s>>h%s_%s_CB", Vars[ivar].Data(), VarName[ivar].Data(), SelName[isel].Data()), Form("%s", Selection[isel].Data()));
				h->Write();
				hCB->Write();
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
