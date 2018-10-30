#include "par.h"

    const char *fname = "root://eoscms//store/group/phys_heavyions/kjung/pPb_EposMinBias_5TeV_8022_Forest/EPOS5TeV_GEN_SIM/crab_pPb_EposMinBias_5TeV_8022_Forest_corrCentrTable/161108_122954/0000"; //new
    const char *histoname = "HiForestAOD_";
void mergedselectdirectory(){

    TChain *t_evt = new TChain("hiEvtAnalyzer/HiTree");
    TChain *t_skim = new TChain("skimanalysis/HltTree");
    TChain *t_hlt = new TChain("hltanalysis/HltTree");
    TChain *t_tower = new TChain("rechitanalyzer/tower");
    TChain *t_castor = new TChain("rechitanalyzer/castor");
    TChain *t_zdcrechit = new TChain("rechitanalyzer/zdcrechit");
    TChain *t_zdcdigi = new TChain("rechitanalyzer/zdcdigi");
    TChain *t_nt = new TChain("rechitanalyzer/ntEvent");

    for(int i=1;i<192;i++){
        cout << "open file : " << Form("%s/%s%d.root",fname,histoname,i) << endl;
        t_evt->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_skim->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_hlt->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_tower->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_castor->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_zdcrechit->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_zdcdigi->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_nt->AddFile(Form("%s/%s%d.root",fname,histoname,i));
    }

    TFile *fout = new TFile("/afs/cern.ch/work/q/qixu/public/centrality/Run2Prep/5TeV/mergeskim_pPb_EposMinBias_5TeV_8022_Forest_corrCentrTable_161108_122954.root","recreate");
        fout->cd();
        
        TDirectory *dir = fout->mkdir("hiEvtAnalyzer");
        dir->cd();
        t_evt->CloneTree(-1,"fast");
        
        TDirectory *dir = fout->mkdir("skimanalysis");
        dir->cd();
        t_skim->CloneTree(-1,"fast");
        
        TDirectory *dir = fout->mkdir("hltanalysis");
        dir->cd();
        t_hlt->CloneTree(-1,"fast");
        
        TDirectory *dir = fout->mkdir("rechitanalyzer");
        dir->cd();
        t_tower->CloneTree(-1,"fast");
        t_castor->CloneTree(-1,"fast");
        t_zdcrechit->CloneTree(-1,"fast");
        t_zdcdigi->CloneTree(-1,"fast");
        t_nt->CloneTree(-1,"fast");

        fout->Write();
        fout->Close();
}
