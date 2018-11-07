#include "fit.C"
void fitAll(){

    //////////////////////////////////////////
    // Retrieve input files and trees from them 
    const char* dir = "/eos/cms/store/group/phys_heavyions/ygo/XeXe2017/"; //lxplus XeXe
    TFile *infData, *infhydjet, *infepos, *infampt, *infUPC;
    TTree *tref, *thydjet, *tepos, *tampt, *tupc;
    infData = TFile::Open("/eos/cms/store/group/phys_heavyions/ygo/XeXe2017/DATA_MB/crab_HIMinimumBias1/HiForestAOD_MB1.root");
    tref = (TTree*)infData->Get("hiEvtAnalyzer/HiTree");
    tref->AddFriend("hltanalysisReco/HltTree");
    tref->AddFriend("skimanalysis/HltTree");
    //tref->AddFriend("anaTrack/trackTree");

    infhydjet = TFile::Open("/eos/cms/store/group/phys_heavyions/ygo/XeXe2017/MC_EPOSLHC/HiForestAOD_EPOSLHC.root");
    thydjet = (TTree*)infhydjet->Get("hiEvtAnalyzer/HiTree");
    thydjet->AddFriend("skimanalysis/HltTree");
    thydjet->AddFriend("hltanalysisReco/HltTree");
/*
    infepos = TFile::Open("/eos/cms/store/group/phys_heavyions/ygo/XeXe2017/MC_HYDJET/HiForestAOD_HYDJET.root");
    tepos = (TTree*)infepos->Get("hiEvtAnalyzer/HiTree");
    tepos->AddFriend("skimanalysis/HltTree");
    tepos->AddFriend("hltanalysisReco/HltTree");
    
    infampt = TFile::Open("/eos/cms/store/group/phys_heavyions/jmartinb/XeXe2016/MC/HiForestAOD_AMPT_XeXe544TeV.root");
    tampt = (TTree*)infampt->Get("hiEvtAnalyzer/HiTree");
    tampt->AddFriend("skimanalysis/HltTree");
    tampt->AddFriend("hltanalysisReco/HltTree");
*/  
    
    infUPC = TFile::Open("/eos/cms/store/group/phys_heavyions/kilee/HiForestAOD_Starlight_Single_th4.root");
    tupc = (TTree*)infUPC->Get("hiEvtAnalyzer/HiTree");
    //tupc->AddFriend("hltanalysisReco/HltTree");
    tupc->AddFriend("skimanalysis/HltTree");

    //////////////////////////////////////////
    // Trigger condition 
    TString st_trig = "HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_v1";
    const int nTrig = 20;
    for(int i=1;i<=nTrig;++i){
        st_trig += Form(" || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part%d_v1",i);
    }
    TCut trigger = Form("%s",st_trig.Data());
    TCut PVcut = "pPAprimaryVertexFilter";
    TCut totEvtSel = PVcut && "pBeamScrapingFilter && phfCoincFilter3";
    TCut totCut = trigger && totEvtSel; 


    //////////////////////////////////////////
    // Run! 
    //void fit(TTree* tref=0, TTree* t=0, TTree* tupc=0, string var = "hiHF", bool doDiffBins=0, int nbins = 2000, double xmin=0 double xmax=4000, double normRangeMin=200, double normRangeMax=3500, TCut dataCut = "pPAprimaryVertexFilter && pBeamScrapingFilter && phfCoincFilter4", string cutname = "PV_BS_HF4", string mc="epos", double scaleMin=0.8, bool doMCevtCut=1, bool doNSDselection=0, string dataCap="XeXe 5.44 TeV")
    static const int Nvar = 1;
    string vars[Nvar] = {"hiHF"};//,"hiHFECut"};/*,"hiET","hiEB","hiEE","hiNtracks","hiNpix","hiNpixelTracks"};*/
    double scaleMin_hydjet[Nvar] = {0.95};
    
    //static const int Nvar = 2;
    //string vars[Nvar] = {"hiHF","hiHFECut"};/*,"hiET","hiEB","hiEE","hiNtracks","hiNpix","hiNpixelTracks"};*/
    //double scaleMin_hydjet[Nvar] = {0.8,0.8};

    for(int i = 0; i < Nvar; ++i){
        fit(tref,thydjet,tupc,vars[i].data(),0,2000,0,4000,200,3500, totCut, "HFORtirg_PV_BS_HF3","HYDJET",scaleMin_hydjet[i],0,0,"XeXe 5.44 TeV MB1");
        //fit(tref,thydjet,tupc,vars[i].data(),0,1000,0,4000,200,3500, totCut, "HFORtirg_PV_BS_HF3","HYDJET",scaleMin_hydjet[i],0,0,"XeXe 5.44 TeV MB1");
        //fit(tref,thydjet,tupc,vars[i].data(),1,1000,0,4000,200,3500, totCut, "HFORtirg_PV_BS_HF3","HYDJET",scaleMin_hydjet[i],0,0,"XeXe 5.44 TeV MB1");
        //fit(tref,thydjet,tupc,vars[i].data(),0,1000,0,4000,200,3500, totCut, "HFORtirg_PV_BS_HF3","HYDJET",scaleMin_hydjet[i],1,0,"XeXe 5.44 TeV MB1");
    }

}
