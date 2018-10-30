// Get event selection filter efficiency as a function of run number. 
// 
//
// Author Yeonju Go
// 23 Nov 2015
//
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
#include "stdio.h"
#include "../../utils.h"

const int nfilter = 5;
const char* evtfilter[] = {"","pprimaryVertexFilter", "phfCoincFilter3", "pclusterCompatibilityFilter", "pcollisionEventSelection"};
void GetfilterEff( int run=262921,
        const char* trig = "HLT_HIL1MinimumBiasHF1AND_v1")
{
    TH1::SetDefaultSumw2();

    const char* fname = Form("/home/goyeonju/CMS/Files/centrality/EventTree_PbPb_data_HIMinimumBias2_run%d_15Feb.root",run);
    TFile *fin = TFile::Open(fname);
    TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
    TTree *t_hlt = (TTree*) fin -> Get("hltanalysis/HltTree");
    t_skim -> AddFriend(t_hlt);

    double Nevt_total = t_skim -> GetEntries(Form("(%s==1)",trig));
    double Nevt_evtSelected[nfilter];
    double Nevt_eff[nfilter];
    for(int i=0;i<nfilter;i++){  
        Nevt_evtSelected[i] = t_skim -> GetEntries(Form("(%s==1)&&(%s)",trig,evtfilter[i]));
        Nevt_eff[i] = Nevt_evtSelected[i]/Nevt_total;
    }

    ofstream outputf(Form("filterEff_%s.txt", trig), ios::app);
    outputf << run << "\t" << Nevt_total << "\t"  
        << setprecision(4) <<Nevt_eff[0]*100 << "\t" 
        << setprecision(4) <<Nevt_eff[1]*100 << "\t" 
        << setprecision(4) <<Nevt_eff[2]*100 << "\t" 
        << setprecision(4) <<Nevt_eff[3]*100 << "\t" 
        << setprecision(4) <<Nevt_eff[4]*100 << "\t" 
        << endl;
    outputf.close();
}

void filterEff_runDep(){

    int run[] =
       {262620, 262640, 262656, 262694, 262695, 262703,
       262726, 262735, 262768, 262784, 262811, 262816,
       262818, 262834, 262837, 262893, 262921, 262988,
       263005, 263022, 263035, 263261, 263286, 263322,
       263333, 263349, 263362, 263379, 263400, 263410,
       263502, 263584, 263604, 263614};

    const int Nrun = sizeof(run)/sizeof(int);

    for(int i=0; i<Nrun; i++){
        //GetfilterEff(run[i], "HLT_HIL1MinimumBiasHF1AND_v1");
        //GetfilterEff(run[i], "HLT_HIL1MinimumBiasHF2AND_v1");
        GetfilterEff(run[i], "HLT_HIL1MinimumBiasHF2AND_part1_v1");
    }
}
