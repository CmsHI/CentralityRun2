// to get the efficiency numbers for each event selection filters.
// Author Yeonju Go 
// 20 Nov 2015
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
#include <iostream>
#include "yjUtility.h"
#include "par.h"

void filterEffCalculator()
{
    //trig : HLT_L1MinimumBiasHF1OR_part1_v1
    //const char* evtfilter[] = {"","HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1","pBeamScrapingFilter","pPAprimaryVertexFilter","pHBHENoiseFilterResultProducer","pVertexFilterCutG","phfCoincFilter3"};
    //const char* evtfilter[] = {"1","pBeamScrapingFilter","pPAprimaryVertexFilter","pHBHENoiseFilterResultProducer","pVertexFilterCutG","phfCoincFilter3"};
    //const char* trig = " L1_MinimumBiasHF0_OR_BptxAND_Final";
        /*
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 ||\
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 ||\
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || \
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || \
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || \
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 ||\
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1 \
        */
    ofstream fout("Eff.log");

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    SetHistTitleStyle();
    SetyjPadStyle();
    TChain *t = new TChain("hiEvtAnalyzer/HiTree");
    TChain *t_skim = new TChain("skimanalysis/HltTree");
    TChain *t_hlt = new TChain("hltanalysis/HltTree");
    for(int i=startpoint;i<nfile;i++){
        cout << "open file : " << Form("%s/%s%d.root",fname,histoname,i) << endl;
        t->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_skim->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_hlt->AddFile(Form("%s/%s%d.root",fname,histoname,i));
    }
    t-> AddFriend(t_skim);
    t-> AddFriend(t_hlt);

    int totalN = t->GetEntries();
    int minBiasN = t->GetEntries(Form("%s==1",trig));
    fout << "Number of "<< trig << " = " << minBiasN << endl;
    fout << "trig efficiency " << " = " << (double)minBiasN/(double)totalN*100 << "\t%" << endl;
    for(int i=1;i<nfilter;i++){
        const char* tmpTrig = Form("(%s==1)",evtfilter[i]);
        int filterN = t->GetEntries(tmpTrig);
        //fout << evtfilter[i] << "::: efficiency = \t" << filterN << "\t/\t" << minBiasN << " =\t" << (double)filterN/(double)minBiasN*100 <<"\t%" << endl;
        fout << evtfilter[i] << "::: efficiency = \t" << filterN << "\t/\t" << totalN << " =\t" << (double)filterN/(double)totalN*100 <<"\t%" << endl;
    }
    fout << endl;
    for(int i=1;i<nfilter;i++){
        TString tmpTrig1 = Form("(%s==1",trig);
        for(int j=1;j<=i;j++){
            tmpTrig1 += Form(" && %s==1",evtfilter[j]);
        }
        tmpTrig1 += ")";
        int filterN = t->GetEntries(tmpTrig1);
        //fout << evtfilter[i] << "::: accu efficiency = \t" << filterN << "\t/\t" << minBiasN << " =\t" << (double)filterN/(double)minBiasN*100 <<"\t%" << endl;
        fout << evtfilter[i] << "::: accu efficiency = \t" << filterN << "\t/\t" << totalN << " =\t" << (double)filterN/(double)totalN*100 <<"\t%" << endl;
    }
}

