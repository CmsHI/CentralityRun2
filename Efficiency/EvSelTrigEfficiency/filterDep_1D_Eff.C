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
#include "analysisUtils.h"
#include "par.h"

//////////// modify here as you want!! ////////////////
const double dy= 0.5;
const int colhere[] = {1,1,2,4,6,9};
void Get1DEffPlots(TTree* t_evt=0, TString v1="hiHF",int xbin=200, double xmin=0, double xmax=4500, double *bin=NULL, TCut cut="",TString cap="", bool isPassed=0,bool isUni=0);
TF1* FitEff(TH1* h, double xmin, double xmax, TString form);
TFile* outf = new TFile(Form("histfiles/filterEff_1D_%s.root",type.Data()), "RECREATE");

void filterDep_1D_Eff(
     //   const char* fname="/afs/cern.ch/work/q/qixu/public/centrality",
        /*
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 ||\
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 ||\
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || \
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || \
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || \
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 ||\
        HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1 \
        */ bool isUni=0)
{
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);
    SetHistTitleStyle(0.06,0.04);
    SetyjPadStyle();

    TChain *t_evt = new TChain("hiEvtAnalyzer/HiTree");
    TChain *t_skim = new TChain("skimanalysis/HltTree");
    TChain *t_hlt = new TChain("hltanalysis/HltTree");
//    TChain *t_gen = new TChain("HiGenParticleAna/hi");
    TChain *t_pp = new TChain("ppTrack/trackTree");
    for(int i=startpoint;i<nfile;i++){
        cout << "open file : " << Form("%s/%s%d.root",fname,histoname,i) << endl;
        t_evt->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_skim->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_hlt->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_pp->AddFile(Form("%s/%s%d.root",fname,histoname,i));
//        t_gen->AddFile(Form("%s/%s%d.root",fname,histoname,i));
    }
    t_evt-> AddFriend(t_skim);
    t_evt-> AddFriend(t_hlt);
    t_evt-> AddFriend(t_pp);
//    t_evt-> AddFriend(t_gen);

    double Nevt_t = t_evt -> GetEntries(cut);
    cout << "cut is '" << cut.GetTitle() << "' and # of events = " << Nevt_t << endl;

/*
    Get1DEffPlots(t_evt, "hiHF",100,0,5000,cut,Form("%s fullRange",type.Data()),1,isUni);
    Get1DEffPlots(t_evt, "hiHFplusEta4",100,0,5000,cut,Form("%s fullRange",type.Data()),1,isUni);
    Get1DEffPlots(t_evt, "hiHFhit",100,0,4000,cut,Form("%s fullRange",type.Data()),1,isUni);
    Get1DEffPlots(t_evt, "hiNpix",100,0,10000,cut,Form("%s fullRange",type.Data()),1,isUni);
    Get1DEffPlots(t_evt, "hiBin",100,0,10000,cut,Form("%s fullRange",type.Data()),1,isUni);
    Get1DEffPlots(t_evt, "hiNtracks",100,0,10000,cut,Form("%s fullRange",type.Data()),1,isUni);
   // Get1DEffPlots(t_evt, "hiZDC",100,0,50000,cut,Form("%s fullRange",type.Data()),1,isUni);
    Get1DEffPlots(t_evt, "hiET",100,0,2000,cut,Form("%s fullRange",type.Data()),1,isUni);
    Get1DEffPlots(t_evt, "hiEE",100,0,4000,cut,Form("%s fullRange",type.Data()),1,isUni);
    Get1DEffPlots(t_evt, "hiEB",100,0,5000,cut,Form("%s fullRange",type.Data()),1,isUni);
*/
    const int nHFEta4 = 140;
    double binHFEta4[nHFEta4+1]={};
        binHFEta4[0] = 0;
        binHFEta4[1] = 0.4;
        binHFEta4[2] = 1.0;
    for(int i = 3; i<nHFEta4+1;i++){
        binHFEta4[i] = i-1;
    }
    double *bin=NULL;

    Get1DEffPlots(t_evt, "hiHF",200,0,200,bin,cut,Form("%s",type.Data()),1,1);
    Get1DEffPlots(t_evt, "hiHFplus",100,0,100,bin,cut,Form("%s",type.Data()),1,1);
    Get1DEffPlots(t_evt, "hiHFminus",100,0,100,bin,cut,Form("%s",type.Data()),1,1);
    Get1DEffPlots(t_evt, "hiHFplusEta4",nHFEta4,0,0,binHFEta4,cut,Form("%s",type.Data()),1,0);
    Get1DEffPlots(t_evt, "hiHFminusEta4",nHFEta4,0,0,binHFEta4,cut,Form("%s",type.Data()),1,0);
    Get1DEffPlots(t_evt, "hiHFhit",100,0,1000,bin,cut,Form("%s",type.Data()),1,1);
    Get1DEffPlots(t_evt, "hiNpix",100,0,300,bin,cut,Form("%s",type.Data()),1,1);
    Get1DEffPlots(t_evt, "hiBin",200,0,200,bin,cut,Form("%s",type.Data()),1,1);
    Get1DEffPlots(t_evt, "hiNtracks",200,0,200,bin,cut,Form("%s",type.Data()),1,1);
   // Get1DEffPlots(t_evt, "hiZDC",100,0,5000,cut,Form("%s",type.Data()),1,1);
    Get1DEffPlots(t_evt, "hiET",20,0,20,bin,cut,Form("%s",type.Data()),1,1);
    Get1DEffPlots(t_evt, "hiEE",20,0,20,bin,cut,Form("%s",type.Data()),1,1);
    Get1DEffPlots(t_evt, "hiEB",80,0,80,bin,cut,Form("%s",type.Data()),1,1);
    Get1DEffPlots(t_evt, "nTrk",320,0,320,bin,cut,Form("%s",type.Data()),1,1);
    Get1DEffPlots(t_evt, "nTrkVtx",300,0,300,bin,cut,Form("%s",type.Data()),1,1);
    outf->Close();
    
}

void Get1DEffPlots(TTree* t_evt, TString v1, int xbin, double xmin, double xmax, double *bin, TCut cut, TString cap, bool isPassed, bool isUni)
{
    TCanvas *c_tot = new TCanvas(Form("c_tot_%s_%s",v1.Data(),cap.Data()), "c_tot", 300,600);
    c_tot->Divide(1,2);
    TCut totcut[nfilter];
    for(int i=0; i<nfilter;i++){
        totcut[i] = "";
        for(int j=1;j<=i;j++){
        totcut[i] = totcut[i] && Form("%s==%d",evtfilter[j],(int)isPassed);
        }
    }

    TCanvas *c_temp= new TCanvas(Form("c_temp_%s",cap.Data()), "", 300,300);
    TH1D *h1D[nfilter];
    TH1D *h1D_eff[nfilter];
    TF1 *fun_eff[nfilter];
    TParameter<float>* eff[nfilter];
    for(int i=0; i<nfilter; i++){
        if(isUni)
        h1D[i] = new TH1D(Form("h_%s_%s",v1.Data(),evtfiltershort[i]), Form(";%s;Events", v1.Data()), xbin, xmin, xmax );
        else{
        h1D[i] = new TH1D(Form("h_%s_%s",v1.Data(),evtfiltershort[i]), Form(";%s;Events", v1.Data()), xbin, bin );
        xmin = bin[0];
        xmax = bin[xbin];
        }
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
        if(i!=0) h1D_eff[i] -> Divide(h1D[i],h1D[1],1,1,"B");
        else h1D_eff[i] -> Divide(h1D[i],h1D[0],1,1,"B");
        if(i==1) h1D_eff[i] -> Divide(h1D[i],h1D[0],1,1,"B");
        eff[i] = new TParameter<float>(Form("IntEff_%s",evtfiltershort[i]),0);
        if(i>1) eff[i]->SetVal(h1D[i]->GetEntries()/h1D[1]->GetEntries());
        else eff[i]->SetVal(h1D[i]->GetEntries()/h1D[0]->GetEntries());
    }
    if(v1=="hiBin")
    TLegend* l1 = new TLegend(0.30, 0.19, 0.9, 0.58, Form("%s",cap.Data()));
    else
    TLegend* l1 = new TLegend(0.40, 0.64, 0.9, 0.94, Form("%s",cap.Data()));
    legStyle(l1);
    l1->SetTextSize(0.04);
    for(int i=0; i<nfilter;i++){
 //       if(i==0) l1->AddEntry(h1D[i], "NO FILTER", "l");  else
            l1->AddEntry(h1D[i], evtfiltershort[i]); 
    }

    if(v1=="hiBin")h1D[0]->SetMinimum(0);
    for(int i=0;i<nfilter;i++){
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
                jumSun(xmin,1,xmax/5.,1);
            }
            else h1D_eff[i] -> DrawCopy("ep same"); 
        }
//        fun_eff[i] = FitEff(h1D_eff[i],xmin,xmax,"exp");
//        fun_eff[i] -> SetName(Form("funEff_%s",evtfiltershort[i]));
    }
    c_tot->SaveAs(Form("figure/filterDep_1D_Eff_%s_%s.gif",v1.Data(),cap.Data()));

    outf->cd();
    outf->mkdir(Form("%s",v1.Data()));
    TDirectory *dir = (TDirectory*)outf->Get(Form("%s",v1.Data()));
    dir->cd();
    for(int i=0; i<nfilter; i++){
        h1D[i]->Write();
        h1D_eff[i]->Write();
        eff[i]->Write();
//        fun_eff[i]->Write();
    }
    c_tot->Write();

    for(int i=0; i<nfilter; i++){
        delete h1D[i];
        delete h1D_eff[i];
    }
    delete c_temp;
}
/*
void Get1DEffPlots_loop(TTree* t_evt, TString v1, int xbin, double xmin, double xmax, TCut cut, TString cap, bool isPassed, bool isUni)
{
    TCanvas *c_tot = new TCanvas(Form("c_tot_%s_%s",v1.Data(),cap.Data()), "c_tot", 300,600);
    c_tot->Divide(1,2);
   
    TCut totcut[nfilter];
    for(int i=0; i<nfilter;i++){
        totcut[i] = "";
        for(int j=1;j<=i;j++){
        totcut[i] = totcut[i] && Form("%s==%d",evtfilter[j],(int)isPassed);
        }
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

    // add DS selection, use loops
    for(int ievt = 0;ievt<t_evt->GetEntries();ievt++){
        if(!isDSEvent(t_evt,ievt)) continue;
        t_evt->GetEntry(ievt);
    for(int i=0; i<nfilter; i++){
        if(!passEvent(t_evt,totcut[i],ievt)) continue;
        //t_evt->Draw(Form("%s>>+%s",v1.Data(), h1D[i]->GetName()), totcut[i]);
        h1D[i]->Fill(v1);
    }
    }
    for(int i=0; i<nfilter; i++){
        if(i!=0) h1D_eff[i] -> Divide(h1D[i],h1D[1],1,1,"B");
    }
    if(v1=="hiBin")
    TLegend* l1 = new TLegend(0.30, 0.19, 0.9, 0.58, Form("%s",cap.Data()));
    else
    TLegend* l1 = new TLegend(0.40, 0.64, 0.9, 0.94, Form("%s",cap.Data()));
    legStyle(l1);
    l1->SetTextSize(0.04);
    for(int i=0; i<nfilter;i++){
 //       if(i==0) l1->AddEntry(h1D[i], "NO FILTER", "l");  else
            l1->AddEntry(h1D[i], evtfiltershort[i]); 
    }

    if(v1=="hiBin")h1D[0]->SetMinimum(0);
    for(int i=0;i<nfilter;i++){
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
    c_tot->SaveAs(Form("figure/filterDep_1DDS_Eff_%s_%s.gif",v1.Data(),cap.Data()));

    TFile* outf = new TFile(Form("histfiles/filterEff_1DDS_%s_%s.root",v1.Data(),cap.Data()), "RECREATE");
    outf->cd();
    for(int i=0; i<nfilter; i++){
        h1D[i]->Write();
        h1D_eff[i]->Write();
        TParameter<float>* eff = new TParameter<float>(Form("IntEff_%s",evtfiltershort[i]),0);
        eff->SetVal(h1D[i]->GetEntries()/h1D[0]->GetEntries());
        eff->Write();
    }
    c_tot->Write();
    outf->Close();

    for(int i=0; i<nfilter; i++){
        delete h1D[i];
        delete h1D_eff[i];
    }
    delete c_temp;
}
*/

TF1* FitEff(TH1* hEff, double xmin, double xmax, TString form){
    TF1 *fun = NULL;
    if(form == "exp") //exponential fitting
    fun = new TF1("fun","[0] + [1] * exp(-[2]*x)",xmin,xmax);
    hEff -> Fit("fun","R");
    return fun;
}
