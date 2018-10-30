// Author Yeonju Go
// last modification : 2017/03/07 
// 
//Headers{{{
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TPad.h"
#include "TLatex.h"
#include "stdio.h"
#include "../../utils.h"
#include "nTower_Header.h"
//}}}

void makeHist_nTower(float thr=0, const Int_t var = 0, bool PV = true, bool BS = true, const Int_t GF = 0, const Int_t Coin = 1, const Int_t MCver = 0);

void make_nTower(bool PV = true, bool BS = true, const Int_t GF = 0, const Int_t Coin = 1, const Int_t MCver = 0){

    for(int i=0;i<nbin_e;i++){
        cout << "START Energy threshold : " << ethr[i] << endl;
        makeHist_nTower(ethr[i],0, PV, BS, GF, Coin, MCver); 
    }

    for(int i=0;i<nbin_et;i++){
        cout << "START transverse energy threshold : " << etthr[i] << endl;
        makeHist_nTower(etthr[i],1, PV, BS, GF, Coin, MCver); 
    }
}


void makeHist_nTower(float thr, const Int_t var, bool PV, bool BS, const Int_t GF, const Int_t Coin, const Int_t MCver)
{
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle -> SetOptStat(0);
    TH1::SetDefaultSumw2();

//Cuts{{{
	TCut PVfilter;
	if(PV) PVfilter = "pPAprimaryVertexFilter==1";
	else PVfilter = "";
	TCut BSfilter;
	if(BS) BSfilter = "pBeamScrapingFilter==1";
	else BSfilter = "";
	TCut Gfilter;
	if(GF == 1) Gfilter = "pVertexFilterCutG==1";
	else if(GF == 2) Gfilter = "pVertexFilterCutGplus==1";
	else Gfilter = "";
	TCut Coinfilter;
	if(Coin == 1) Coinfilter = "phfCoincFilter1==1";
	else if(Coin == 2) Coinfilter = "phfCoincFilter2==1";
	else if(Coin == 3) Coinfilter = "phfCoincFilter3==1";
	else if(Coin == 4) Coinfilter = "phfCoincFilter4==1";
	else Coinfilter = "";
	TCut Trigger = "(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2==1)";
	TCut Lumi = "(67<lumi&&lumi<1540)";
//}}}

// Get Trees from data & mc files.{{{
    TFile *f[2]; // 0:data, 1:MC

    TChain* t_skim[2];
    TChain* t_hlt[2];
    TChain* t_rec[2];
	TChain* t_evt;
    const char* fname_data_pre = "root://eoscms//eos/cms/store/group/phys_heavyions/dhangal/pr_forests/PAMinimumBias/v5/000/285/530/HiForest_";
    const char* fname_mc_EPOS_pre = "root://eoscms//eos/cms/store/group/phys_heavyions/jmartinb/pPb2016/MC/EPOS8TeV/HiForest_pPbEPOS8TeV_RECODEBUG";
    const char* fname_mc_HIJING_pre = "root://eoscms//eos/cms/store/group/phys_heavyions/jmartinb/pPb2016/MC/EPOS8TeV/HiForest_pPbEPOS8TeV_RECODEBUG";//not available yet

    t_evt = new TChain("hiEvtAnalyzer/HiTree");
    for(int i=0;i<2;i++){ 
        t_skim[i] = new TChain("skimanalysis/HltTree");
        t_hlt[i] = new TChain("hltanalysis/HltTree");
        t_rec[i] = new TChain("rechitanalyzer/tower");
        if(i==0)
			{
				for(Int_t j = 0; j < 5; j++)
				{
					t_skim[0]->Add(Form("%s%d.root", fname_data_pre, j));
					t_hlt[0]->Add(Form("%s%d.root", fname_data_pre, j));
					t_rec[0]->Add(Form("%s%d.root", fname_data_pre, j));
					t_evt->Add(Form("%s%d.root", fname_data_pre, j));
					t_rec[0] -> AddFriend(t_evt);
				}
			}
			else if(MCver == 0)
			{
				t_skim[1]->Add(Form("%s.root", fname_mc_EPOS_pre));
				t_hlt[1]->Add(Form("%s.root", fname_mc_EPOS_pre));
				t_rec[1]->Add(Form("%s.root", fname_mc_EPOS_pre));
			}
			else if(MCver == 1)
			{
				t_skim[1]->Add(Form("%s.root", fname_mc_HIJING_pre));
				t_hlt[1]->Add(Form("%s.root", fname_mc_HIJING_pre));
				t_rec[1]->Add(Form("%s.root", fname_mc_HIJING_pre));
			}
        t_rec[i] -> AddFriend(t_hlt[i]);
        t_rec[i] -> AddFriend(t_skim[i]);
    }
//}}}

// [et] Define et histograms (data/mc , pf tree/rechit tree){{{
    cout << "LET'S DEFINE HISTOGRAMS" << endl;
    TH1F* h1F[2]; //[0:data, 1:mc]
    for(int i=0; i<2; i++){
        h1F[i] = new TH1F(Form("h1F_sample%d",i), ";# of HF towers;Events",100,0,1000);
        h1F[i] -> SetMarkerStyle(marker[i]);
        h1F[i] -> SetMarkerSize(0.9);
        if(i==0) SetHistColor(h1F[i],2);
        else SetHistColor(h1F[i],1);
    }
//}}}

// [rechit] Fill histogram by using <Draw> function in TTree.
    cout << "LET'S FILL HISTOGRAMS FROM TREE" << endl;
    //********************************************************
    TString towerCut; 
    if(var==0) towerCut = Form("abs(tower.eta) > 3.0 && abs(tower.eta)<5.0 && tower.e > %.1f",thr);
    else if(var==1) towerCut = Form("abs(tower.eta) > 3.0 && abs(tower.eta)<5.0 && tower.et > %.1f",thr);
    
    std::clock_t start_loop, end_loop;
    start_loop = std::clock();
    for(int i=0; i<2; i++){
        std::clock_t startT_draw, endT_draw;
        startT_draw = std::clock();
//fill histogram{{{
        if(i==0)
			{
				t_rec[i]->Draw(Form("Sum$(%s)>>+%s",towerCut.Data(),h1F[i]->GetName()), trigger&&PVfilter&&BSfilter&&Gfilter&&Coinfilter);
			}
        else if(i==1) t_rec[i]->Draw(Form("Sum$(%s)>>+%s",towerCut.Data(),h1F[i]->GetName()), "");
//}}}
        h1F[i] = (TH1F*)gDirectory->Get(h1F[i]->GetName());
        endT_draw = std::clock();
        std::cout.precision(6);      // get back to default precision
        std::cout << "DRAW finished in             : " << (endT_draw - startT_draw) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    }

    end_loop = std::clock();
    std::cout.precision(6);      // get back to default precision
    std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "exited loop" << std::endl;

	TString MCname;
	if(MCver == 0) MCname = "EPOS";
	else if(MCver == 1) MCname = "HIJING";
	TString Varname;
	if(var == 0) Varname = "e";
	else if (var == 1) Varname = "et";

    TFile* outFile = new TFile(Form("histfiles/nTower_%sThr%.1f_PV%d_BS%d_GF%d_Coin%d_%s.root",Varname.Data(),thr, PV, BS, GF, Coin, MCname.Data()),"recreate");
    outFile->cd();
    h1F[0]->Write();
    h1F[1]->Write();
    outFile->Close();
} 
