// to calculate the passing fraction of HF tower filter for the starlight and hydjet sample. 
// tower energy threshold could be changed : 1,2,3,4,5 GeV
// the number of tower could be changed : 1,2,3,4,5
// Three kind of options. 
// 1. for generated particles (using HiGenParticleAna/hi tree)
// 2. for reco particles (could use rechitanalyzer/hf or simply skimanalysis/HltTree) 
// 3. for L1 bits (now, 24 Nov, starlight and hydjet samples are both high threshold which is wrong. so need new sample with low L1 threshold)
//
// This macro gives you a set of numbers for three options. 
// We hope option 1,2 and 3 have the same results!
//
// Author Yeonju Go
// 24 Nov 2015

#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1D.h"
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
#include "TLorentzVector.h"
#include "stdio.h"
#include <iostream>
#include "../utils.h"
#include "mbemulationanalysis/L1UpgradeTree.h"
#include "mbemulationanalysis/L1UpgradeTree.C"

const double HFin = 2.87;
const double HFout = 5.02;

void GetHFfraction(const char* fname="root://eoscms//eos/cms/store/group/phys_heavyions/chflores/Foresting_RunPrep2015/STARLIGHTProd/HiForest_Starlight_Merge19112015.root",
        //const char* fname="/mnt/hadoop/cms/store/user/dgulhan/mergedForest/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_merged_forest_1.root",
        float eThr = 3,
        int nTower= 1,
        bool doL1 = 0,
        TString cap="")
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    SetHistTitleStyle();
    SetyjPadStyle();
    TFile *fin = TFile::Open(fname);
    cout << "open file : " << fname << endl;
    TTree *t = (TTree*) fin -> Get("skimanalysis/HltTree");
    TTree *t_gen = (TTree*) fin -> Get("HiGenParticleAna/hi");
    t-> AddFriend(t_gen);
    //TTree *t_l1 = (TTree*) fin -> Get("EmulatorResults/L1UpgradeTree");
    //TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
    //TTree *t = (TTree*) fin -> Get("hltanalysis/HltTree");
    //t-> AddFriend(t_evt);
    //t-> AddFriend(t_skim);
    //t-> AddFriend(t_l1);
    cout << "Energh Threshold : " << eThr << " GeV, # of towers for firing HF : " << nTower << endl;
    //////////////////////////////////////////////////////////////// 
    //////////////////////////////////////////////////////////////// 
    // GEN FRACTION
    cout << " :::::::::::::::::: GEN FRACTION :::::::::::::::::: " << endl;
    Int_t mult;
    std::vector<float>* pt=0;
    std::vector<float>* eta=0;
    std::vector<float>* phi=0;
    std::vector<int>*   pdg=0;
    t->SetBranchAddress("mult",&mult);
    t->SetBranchAddress("pt",&pt);
    t->SetBranchAddress("eta",&eta);
    t->SetBranchAddress("phi",&phi);
    t->SetBranchAddress("pdg",&pdg);

    int totEvt = t->GetEntries();
    int hfANDevt = 0;
    int hfORevt = 0;
    int hfXORevt = 0;
    for(int jentry=0; jentry<totEvt;jentry++){
        if(jentry%1000==0) cout << jentry << " / " << totEvt << endl;
        t->GetEntry(jentry);
        int nPosTower = 0;
        int nNegTower = 0;
        for(int i=0; i<mult;i++){
            float m = -1.;
            if(fabs(pdg->at(i)) == 111) m = 0.13498; //pion 0
            else if(fabs(pdg->at(i)) == 211) m = 0.13957; //pion +-
            else if(fabs(pdg->at(i)) == 321) m = 0.49368; //kaon +-
            else if(fabs(pdg->at(i)) == 130) m = 0.49761; //kaon 0 long
            else if(fabs(pdg->at(i)) == 2212) m = 0.93827; //proton
            else if(fabs(pdg->at(i)) == 2112) m = 0.93957; //neutron
            else if(fabs(pdg->at(i)) == 22) m = 0.; //photon
            else if(fabs(pdg->at(i)) == 11) m = 0.000511; //electron
            else if(fabs(pdg->at(i)) == 13) m = 0.10566; //muon
            else continue;
            TLorentzVector * v = new TLorentzVector();  
            v->SetPtEtaPhiM(pt->at(i),eta->at(i),phi->at(i),m);
            double e = v->Energy();
            delete v;
            if(e<eThr) continue;
            if( (eta->at(i)>HFin ) && ( (eta->at(i))<HFout ) ) nPosTower++;
            else if( (eta->at(i)>-HFout ) && ( (eta->at(i))<-HFin ) ) nNegTower++;
            else continue;
        }
        if(nPosTower>=nTower && nNegTower>=nTower) { hfANDevt++; hfORevt++; }
        else if( ( nPosTower>=nTower && nNegTower<nTower ) || ( nPosTower<nTower && nNegTower>=nTower )) { hfORevt++; hfXORevt++; }
    }

    cout << " HF AND Efficiency : \t" << hfANDevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfANDevt/totEvt*100 << "\t%"<< endl;
    cout << " HF OR Efficiency : \t" << hfORevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfORevt/totEvt*100 << "\t%"<< endl;
    cout << " HF XOR Efficiency : \t" << hfXORevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfXORevt/totEvt*100 << "\t%"<< endl;

    //////////////////////////////////////////////////////////////// 
    //////////////////////////////////////////////////////////////// 
    // RECO FRACTION
    cout << " :::::::::::::::::: RECO FRACTION ::::::::::::::::: " << endl;
    hfANDevt = 0;
    hfORevt = 0;
    hfXORevt = 0;
    if(nTower==1){
        hfANDevt = t->GetEntries(Form("phfCoincFilter==1"));
        hfORevt = t->GetEntries(Form("phfPosFilter==1 || phfNegFilter==1"));
        hfXORevt = t->GetEntries(Form("(phfPosFilter==1 || phfNegFilter==1) && phfCoincFilter==0"));
    } else {
        hfANDevt = t->GetEntries(Form("phfCoincFilter%d==1",nTower));
        hfORevt = t->GetEntries(Form("phfPosFilter%d==1 || phfNegFilter%d==1",nTower,nTower));
        hfXORevt = t->GetEntries(Form("(phfPosFilter%d==1 || phfNegFilter%d==1) && phfCoincFilter%d==0",nTower,nTower,nTower));
    } 

    cout << " HF AND Efficiency : \t" << hfANDevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfANDevt/totEvt*100 << "\t%"<< endl;
    cout << " HF OR Efficiency : \t" << hfORevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfORevt/totEvt*100 << "\t%"<< endl;
    cout << " HF XOR Efficiency : \t" << hfXORevt << "\t/\t" << totEvt << "\t=\t"<< (double)hfXORevt/totEvt*100 << "\t%"<< endl;

    //////////////////////////////////////////////////////////////// 
    //////////////////////////////////////////////////////////////// 
    // L1 FRACTION
    hfANDevt = 0;
    hfORevt = 0;
    hfXORevt = 0;
    if(doL1){
        cout << " :::::::::::::::::: L1 FRACTION :::::::::::::::::: " << endl;
        L1UpgradeTree * tl1 = new L1UpgradeTree(fname);
        tl1->Loop(99999999);
        tl1->GetEntries();
    }
}
