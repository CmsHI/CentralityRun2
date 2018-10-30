// to draw the distribution of centrality variables
// with the different event selection filters in one canvas.
// should modify evtfilter definition part
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
#include <iostream>
#include "yjUtility.h"
#include "par.h"

//////////// modify here as you want!! ////////////////
const double dy= 0.5;
const int colhere[] = {1,1,2,4,6,9};
void Get1DEffPlots(TTree* t_evt=0, TString v1="hiHF",int xbin=200, double xmin=0, double xmax=4500, TCut cut="",TString cap="", bool isPassed=0,bool isAOD=0);

void filter_1D_TrigDiff(
     //   const char* fname="/afs/cern.ch/work/q/qixu/public/centrality",
        /*
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 ||\
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 ||\
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || \
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || \
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || \
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 ||\
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1 \
        */ bool isAOD=0)
{
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);
    SetHistTitleStyle(0.06,0.04);
    SetyjPadStyle();

    TChain *t_evt = new TChain("hiEvtAnalyzer/HiTree");
    TChain *t_skim = new TChain("skimanalysis/HltTree");
    TChain *t_hlt = new TChain("hltanalysis/HltTree");
    for(int i=startpoint;i<nfile;i++){
        cout << "open file : " << Form("%s/%s%d.root",fname,histoname,i) << endl;
        t_evt->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_skim->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_hlt->AddFile(Form("%s/%s%d.root",fname,histoname,i));
    }
    t_evt-> AddFriend(t_skim);
    t_evt-> AddFriend(t_hlt);

    double Nevt_t = t_evt -> GetEntries(cut);
    cout << "cut is '" << cut.GetTitle() << "' and # of events = " << Nevt_t << endl;

/*
    Get1DEffPlots(t_evt, "hiHF",100,0,5000,cut,Form("%s fullRange",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiHFplusEta4",100,0,5000,cut,Form("%s fullRange",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiHFhit",100,0,4000,cut,Form("%s fullRange",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiNpix",100,0,10000,cut,Form("%s fullRange",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiBin",100,0,10000,cut,Form("%s fullRange",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiNtracks",100,0,10000,cut,Form("%s fullRange",type.Data()),1,isAOD);
   // Get1DEffPlots(t_evt, "hiZDC",100,0,50000,cut,Form("%s fullRange",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiET",100,0,2000,cut,Form("%s fullRange",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiEE",100,0,4000,cut,Form("%s fullRange",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiEB",100,0,5000,cut,Form("%s fullRange",type.Data()),1,isAOD);
*/
    Get1DEffPlots(t_evt, "hiHF",100,0,100,cut,Form("%s",type.Data()),1,isAOD);
 /*
    Get1DEffPlots(t_evt, "hiHFplus",100,0,100,cut,Form("%s",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiHFminus",100,0,100,cut,Form("%s",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiHFplusEta4",70,0,70,cut,Form("%s",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiHFminusEta4",100,0,100,cut,Form("%s",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiHFhit",100,0,1000,cut,Form("%s",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiNpix",100,0,300,cut,Form("%s",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiBin",200,0,200,cut,Form("%s",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiNtracks",100,0,200,cut,Form("%s",type.Data()),1,isAOD);
   // Get1DEffPlots(t_evt, "hiZDC",100,0,5000,cut,Form("%s",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiET",20,0,20,cut,Form("%s",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiEE",20,0,20,cut,Form("%s",type.Data()),1,isAOD);
    Get1DEffPlots(t_evt, "hiEB",80,0,80,cut,Form("%s",type.Data()),1,isAOD);
    */

}

void Get1DEffPlots(TTree* t_evt, TString v1, int xbin, double xmin, double xmax, TCut cut, TString cap, bool isPassed, bool isAOD)
{
    TCanvas *c_tot = new TCanvas(Form("c_tot_%s_%s",v1.Data(),cap.Data()), "c_tot", 300,600);
    c_tot->Divide(1,2);
    TCut totcut[nfilter];
    for(int i=0; i<nfilter;i++){
        totcut[i] = "";
        for(int j=1;j<=i;j++){
        totcut[i] = totcut[i] && Form("%s==%d",evtfilter[j],(int)isPassed);
        }
        totcut[2] = totcut[0] && Form("HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1");
        cout<<totcut[i].GetTitle()<<endl;
    }

    TCanvas *c_temp= new TCanvas(Form("c_temp_%s",cap.Data()), "", 300,300);
    TH1D *h1D[nfilter];
    TH1D *h1D_eff[nfilter];
    for(int i=0; i<nfilter; i++){
        h1D[i] = new TH1D(Form("h_%s_%s",v1.Data(),evtfiltershort[i]), Form(";%s;Events", v1.Data()), xbin, xmin,xmax );
        h1D[i]->SetMarkerStyle(24+i); 
        h1D[i]->SetMarkerSize(0.8);
        h1D_eff[i] = (TH1D*)h1D[i]->Clone(Form("hEff_%s_%s",v1.Data(),evtfiltershort[i]));
        h1D_eff[i] -> SetTitle(Form(";%s;Filter Rate",v1.Data()));
        SetHistColor(h1D[i],colhere[i]);
        SetHistColor(h1D_eff[i],colhere[i]);
    }
    for(int i=0; i<nfilter; i++){
        t_evt->Draw(Form("%s>>+%s",v1.Data(), h1D[i]->GetName()), totcut[i]);
        h1D[i]=(TH1D*)gDirectory->Get(h1D[i]->GetName());
        if(i!=0) h1D_eff[i] -> Divide(h1D[i],h1D[0],1,1,"B");
        cout<<h1D[i]->GetEntries()<<endl;
    }
    if(v1=="hiBin")
    TLegend* l1 = new TLegend(0.30, 0.19, 0.9, 0.58, Form("%s",cap.Data()));
    else
    TLegend* l1 = new TLegend(0.30, 0.54, 0.9, 0.94, Form("%s",cap.Data()));
    legStyle(l1);
    l1->SetTextSize(0.04);
    for(int i=0; i<1;i++){
 //       if(i==0) l1->AddEntry(h1D[i], "NO FILTER", "l");  else
            l1->AddEntry(h1D[i], evtfiltershort[i]); 
    }
    l1->AddEntry(h1D[1], "HLT_PAL1MinimumBiasHF");
    l1->AddEntry(h1D[1], "OR_SinglePixelTrack_part1_v1","");
    l1->AddEntry(h1D[2], "HLT_PAL1MinimumBiasHF");
    l1->AddEntry(h1D[2], "AND_SinglePixelTrack_v1","");

    if(v1=="hiBin")h1D[0]->SetMinimum(0);
    for(int i=0;i<3;i++){
        c_tot->cd(1);
        if(i==0) h1D[i]->DrawCopy("hist");
        else h1D[i] -> DrawCopy("ep same");
        l1 -> Draw();
        if(i!=0){
            c_tot->cd(2);
            //h1D_eff[i]->GetYaxis()-> SetRangeUser(h1D_eff[nfilter-1]->GetMinimum()*0.9,1.0);
            h1D_eff[i]->GetYaxis()-> SetRangeUser(0.0,1.0);
            if(i==1) {
                h1D_eff[i] -> DrawCopy("ep"); 
                jumSun(xmin,1,xmax,1);
            }
            else h1D_eff[i] -> DrawCopy("ep same"); 
        }
    }
    c_tot->SaveAs(Form("figure/filter_1D_TrigDff_%s_%s.gif",v1.Data(),cap.Data()));

    TFile* outf = new TFile(Form("histfiles/filterEff_1D_%s_%s.root",v1.Data(),cap.Data()), "RECREATE");
    outf->cd();
    for(int i=0; i<nfilter; i++){
        h1D[i]->Write();
        h1D_eff[i]->Write();
        TParameter<float>* eff = new TParameter<float>(Form("IntEff_%s",evtfiltershort[i]),0);
        eff->SetVal(h1D[i]->Integral()/h1D[0]->Integral());
        cout<<eff->GetVal()<<endl;
        eff->Write();
    }
//    c_tot->Write();
//    outf->Close();

    for(int i=0; i<nfilter; i++){
        delete h1D[i];
        delete h1D_eff[i];
    }
    delete c_temp;
}
