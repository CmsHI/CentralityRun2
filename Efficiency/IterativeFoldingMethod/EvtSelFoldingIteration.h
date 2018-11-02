#ifndef EvtSelFoldingIteration_h
#define EvtSelFoldingIteration_h
// This macro 
// Created : 12 June 2017
// Author : Yeonju Go
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
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
#include "stdio.h"
#include "string.h"
#include <iostream>
#include "TEfficiency.h"
#include "../../analysisUtils.h"
#include "../../yjUtility.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"
#include "TParameter.h"
#include "TObjArray.h"

class EffHists {
    public:
        TH1D *h1D_mc_tot;
        TH1D *h1D_mc_evtSel;
        TH1D *h1D_mc_tot_weighted;
        TH1D *h1D_mc_evtSel_weighted;
        TH1D *h1D_mc_filterRate;
        TH1D *h1D_data_tot;
        TH1D *h1D_data_evtSel;
        TH1D *h1D_data_corrected;
        TH1D *h1D_weightFactor;
        double MC_integFilterRate=0;
        double DATA_integEvtSelEff=0;
};

EffHists getEvtSelEff(bool doNoTrigEff, const string mcName, const TString var, const int num, TH1D* weightHist[], const TString weightVar[]){
    EffHists res;
    TString noTrigEffCap = "";
    if(doNoTrigEff) noTrigEffCap = "_noTrigEff";
    int run = 286471;
    
    TString fout_name = Form("output/evtSelEff_foldingMethod%s_run%d_%s_%s.root",noTrigEffCap.Data(),run,mcName.data(),var.Data());
    if(num!=0) {
        string temp = "WeightedBy";
        for(int i=0;i<num;++i){
            temp = temp + "_" + weightVar[i];
        }    
        fout_name=Form("output/evtSelEff_foldingMethod%s_run%d_%s_%s_%s.root",noTrigEffCap.Data(),run,mcName.data(),var.Data(),temp.data());
        //cout << fout_name << endl;
    }
   
    TFile* temp_f = new TFile(fout_name,"READ");
    if ( !(temp_f->IsZombie()) ){
        cout << "This file already exists : " << fout_name << endl;
        return res; 
    } else {
        delete temp_f;
    }
    
    // STRAT!!
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);
    cout << "::::::: getEvtSelEff() Starts for "<< var << endl;
    cout << "::::::: Iteration  : "<<num<< endl;
    if(num!=0) {
        for(int i=0;i<num;++i)  { 
            cout << "::::::: WeightHist : "<< weightHist[i]->GetName() << ", WeightVariable = " << weightVar[i] << endl;
        }
    }
    /// LOAD TREES AND CREATE HISTOS
    TChain * tdata = new TChain("hiEvtAnalyzer/HiTree","");
    TChain * tskimanalysis = new TChain("skimanalysis/HltTree","");
    TChain * thltanalysis = new TChain("hltanalysis/HltTree","");
    const char* fileData = Form("/home/goyeonju/CMS/Files/centrality/pPb2016/HiForestAOD_MinimumBiasDATA_pPb8TeV_run%d.root",run);
    tdata->Add(fileData);
    tskimanalysis->Add(fileData);
    thltanalysis->Add(fileData);
    tdata->AddFriend(tskimanalysis);
    tdata->AddFriend(thltanalysis);

    // MC File and Tree
    TChain* tmc = new TChain("hiEvtAnalyzer/HiTree","");
    TChain * tskimanalysismc = new TChain("skimanalysis/HltTree","");
    TChain * thltanalysismc = new TChain("hltanalysis/HltTree","");
    const char* fileMC = Form("/home/goyeonju/CMS/Files/centrality/pPb2016/withGen/HiForestAOD_%s_pPb_8TeV.root",mcName.data());
    tmc->Add(fileMC);
    tskimanalysismc->Add(fileMC);
    thltanalysismc->Add(fileMC);
    tmc->AddFriend(tskimanalysismc);
    tmc->AddFriend(thltanalysismc);

    Float_t hiHF,hiHFhit,hiEB,hiEE,hiET,hiHFminus,hiHFplus,hiHFminusEta4,hiHFplusEta4;
    Float_t hiHF_mc,hiHFhit_mc,hiEB_mc,hiEE_mc,hiET_mc,hiHFminus_mc,hiHFplus_mc,hiHFminusEta4_mc,hiHFplusEta4_mc;
    int hiNtracks,hiNpix;
    int hiNtracks_mc,hiNpix_mc, ProcessID_mc;
    Int_t HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8;
    Int_t HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_mc,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_mc,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_mc,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_mc,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_mc,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_mc,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_mc,
    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_mc;
    Int_t pPAprimaryVertexFilter,pBeamScrapingFilter,phfCoincFilter1,pVertexFilterCutGplus;
    Int_t pPAprimaryVertexFilter_mc,pBeamScrapingFilter_mc,phfCoincFilter1_mc,pVertexFilterCutGplus_mc;
    tdata->SetBranchAddress("hiHF",&hiHF);
    tdata->SetBranchAddress("hiHFhit",&hiHFhit);
    tdata->SetBranchAddress("hiEB",&hiEB);
    tdata->SetBranchAddress("hiEE",&hiEE);
    tdata->SetBranchAddress("hiET",&hiET);
    tdata->SetBranchAddress("hiNtracks",&hiNtracks);
    tdata->SetBranchAddress("hiNpix",&hiNpix);
    tdata->SetBranchAddress("hiHFminus",&hiHFminus);
    tdata->SetBranchAddress("hiHFplus",&hiHFplus);
    tdata->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
    tdata->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
    tdata->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1);
    tdata->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2);
    tdata->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3);
    tdata->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4);
    tdata->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5);
    tdata->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6);
    tdata->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7);
    tdata->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8);
    tdata->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
    tdata->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
    tdata->SetBranchAddress("phfCoincFilter1",&phfCoincFilter1);
    tdata->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);

    tmc->SetBranchAddress("hiHF",&hiHF_mc);
    tmc->SetBranchAddress("hiHFhit",&hiHFhit_mc);
    tmc->SetBranchAddress("hiEB",&hiEB_mc);
    tmc->SetBranchAddress("hiEE",&hiEE_mc);
    tmc->SetBranchAddress("hiET",&hiET_mc);
    tmc->SetBranchAddress("hiNtracks",&hiNtracks_mc);
    tmc->SetBranchAddress("hiNpix",&hiNpix_mc);
    tmc->SetBranchAddress("hiHFminus",&hiHFminus_mc);
    tmc->SetBranchAddress("hiHFplus",&hiHFplus_mc);
    tmc->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_mc);
    tmc->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_mc);
    tmc->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_mc);
    tmc->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_mc);
    tmc->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_mc);
    tmc->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_mc);
    tmc->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_mc);
    tmc->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_mc);
    tmc->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_mc);
    tmc->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2",&HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_mc);
    tmc->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter_mc);
    tmc->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_mc);
    tmc->SetBranchAddress("phfCoincFilter1",&phfCoincFilter1_mc);
    tmc->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_mc);
    if(mcName=="EPOSLHC") tmc->SetBranchAddress("ProcessID",&ProcessID_mc);

    int nbins = 700;
    double xmin = 0;
    double xmax = 350; 
    if(var == "hiNpix") xmax = 1800;
    if(var == "hiNtracks") xmax = 300;
    if(var == "hiHFhit"){xmin=-120; xmax = 7500;}
    if(var == "hiEE"||var=="hiET") xmax = 150;
    if(var == "hiEB") xmax = 250;
    if(var == "hiHFplus") xmax = 100;
    if(var == "hiHFminus") xmax = 200;
    if(var == "hiHFplusEta4") xmax = 60;
    if(var == "hiHFminusEta4") xmax = 90;
   // int binWidth3_int = ( xmax*2/3 ) /5/2;//the number of bins of this range = 5
   // int binWidth2_int = ( xmax*4/15 ) /10;//the number of bins of this range = 10
   // int binWidth1_int = ( xmax*1/15 ) /10;//the number of bins of this range = 10 
   // double binWidth3_double = ( xmax*2./3. ) /5./4.;//the number of bins of this range = 5
   // double binWidth2_double = ( xmax*4./15. ) /10.;//the number of bins of this range = 10
   // double binWidth1_double = ( xmax*1./15. ) /15.;//the number of bins of this range = 15 
    double binWidth0(0),binWidth1(0),binWidth2(0),binWidth3(0),binWidth4(0); 
    if(var == "hiNpix"){
        binWidth0 = 2; 
        binWidth1 = 4; 
        binWidth2 = 10;
        binWidth3 = 48;
        binWidth4 = 300;
    } else if(var == "hiNtracks"){
        binWidth0 = 1;
        binWidth1 = 1;
        binWidth2 = 2;
        binWidth3 = 5;
        binWidth4 = 20;
    } else if(var=="hiHF") {
        binWidth0 = 0.25; 
        binWidth1 = 0.5; 
        binWidth2 = 2.0; 
        binWidth3 = 10.0; 
        binWidth4 = 100.0; 
    } else if(var=="hiHFhit"){
        binWidth0 = 5.0; 
        binWidth1 = 15.0; 
        binWidth2 = 70.0; 
        binWidth3 = 150.0; 
        binWidth4 = 500.0; 
    } else if(var=="hiHFplus" || var=="hiHFminus"){
        binWidth0 = 0.15; 
        binWidth1 = 0.25; 
        binWidth2 = 1.0; 
        binWidth3 = 5.0; 
        binWidth4 = 50.0; 
    } else if(var=="hiHFplusEta4" || var=="hiHFminusEta4"){
        binWidth0 = 0.05; 
        binWidth1 = 0.12; 
        binWidth2 = 0.8; 
        binWidth3 = 3.0; 
        binWidth4 = 30.0; 
    }
    cout << "binWidth 1,2,3,4 = " << binWidth1 << ", "<<  binWidth2 << ", "<< binWidth3  << ", "<< binWidth4 << endl; 
    int nDiffBins=0;
    float tempBins[300];
    tempBins[0]=0;
    if(var=="hiHFhit") tempBins[0]=-120;
    for(int i=0;tempBins[i]<xmax;++i){
        if(tempBins[i]<xmax*1./80.) tempBins[i+1] = tempBins[i]+binWidth0;
        else if(tempBins[i]>=xmax*1./80. && tempBins[i]<xmax*1./40.) tempBins[i+1] = tempBins[i]+binWidth1;
        else if(tempBins[i]>=xmax*1./40. && tempBins[i]<xmax*1./10.) tempBins[i+1] = tempBins[i]+binWidth2;
        else if(tempBins[i]>=xmax*1./10. && tempBins[i]<xmax*1./2. ) tempBins[i+1] = tempBins[i]+binWidth3;
        else if(tempBins[i]>=xmax*1./2. ) tempBins[i+1] = tempBins[i]+binWidth4;
        else cout << i << "th out of "<<nDiffBins<<"! somethings wrong";
        nDiffBins++; 
    }
    cout << "Number of bins = " << nDiffBins << endl;
    double bins[nDiffBins+1];
    for(int i=0;i<nDiffBins+1;++i){
        bins[i] = tempBins[i];
        if(i!=nDiffBins) cout << bins[i]<< ",   ";
        else cout << bins[i]<< endl;
    } 

    TH1D* h1D_mc_tot = new TH1D("mc_tot", Form(";%s;Entries",var.Data()),nDiffBins,bins);
    TH1D* h1D_mc_evtSel = (TH1D*) h1D_mc_tot->Clone("mc_evtSel");
    TH1D* h1D_mc_tot_weighted= (TH1D*) h1D_mc_tot->Clone("mc_tot_weighted");
    TH1D* h1D_mc_evtSel_weighted= (TH1D*) h1D_mc_tot->Clone("mc_evtSel_weighted");
    TH1D* h1D_mc_filterRate= (TH1D*) h1D_mc_tot->Clone("mc_filterRate");
    TH1D* h1D_data_tot= (TH1D*) h1D_mc_tot->Clone("data_tot");
    TH1D* h1D_data_evtSel= (TH1D*) h1D_mc_tot->Clone("data_evtSel");
    TH1D* h1D_data_corrected= (TH1D*) h1D_mc_tot->Clone("data_corrected");
    TH1D* h1D_weightFactor= (TH1D*) h1D_mc_tot->Clone("weightFactor");
    
   // TH1D* h1D_mc_tot = new TH1D(Form("mc_tot_%s",var.Data()), Form(";%s;Entries",var.Data()),nDiffBins,bins);
   // TH1D* h1D_mc_evtSel = (TH1D*) h1D_mc_tot->Clone(Form("mc_evtSel_%s",var.Data()));
   // TH1D* h1D_mc_tot_weighted= (TH1D*) h1D_mc_tot->Clone(Form("mc_tot_weighted_%s",var.Data()));
   // TH1D* h1D_mc_evtSel_weighted= (TH1D*) h1D_mc_tot->Clone(Form("mc_evtSel_weighted_%s",var.Data()));
   // TH1D* h1D_mc_filterRate= (TH1D*) h1D_mc_tot->Clone(Form("mc_filterRate_%s",var.Data()));
   // TH1D* h1D_data_tot= (TH1D*) h1D_mc_tot->Clone(Form("data_tot_%s",var.Data()));
   // TH1D* h1D_data_evtSel= (TH1D*) h1D_mc_tot->Clone(Form("data_evtSel_%s",var.Data()));
   // TH1D* h1D_data_corrected= (TH1D*) h1D_mc_tot->Clone(Form("data_corrected_%s",var.Data()));
   // TH1D* h1D_weightFactor= (TH1D*) h1D_mc_tot->Clone(Form("weightFactor_%s",var.Data()));

    // MC Event Loop
    double Nden_mc = 0; 
    double Nnum_mc = 0; 
    int nEntries_mc = tmc->GetEntries();
    for(int iev=0;iev<nEntries_mc;++iev){
        tmc->GetEntry(iev);
        //if ( iev % 50000 == 0 ) std::cout << "MC :: " <<  iev << "/" << nEntries_mc << " : " << iev*100.0/nEntries_mc << std::endl;
        if(mcName=="EPOSLHC" && !(ProcessID_mc==101||ProcessID_mc==105)) continue;
        
        // Trigger selecttion
        if(doNoTrigEff && !(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_mc))
            continue;
        
        double currVar(0.);
        int nWeight = num;
        if(num==0) nWeight=1;
        double currWeightVar[nWeight];
        if(var == "hiET") currVar = hiET_mc;
        if(var == "hiEB")  currVar = hiEB_mc;
        if(var == "hiEE")  currVar = hiEE_mc;
        if(var == "hiHF")  currVar = hiHF_mc;
        if(var == "hiHFhit")  currVar = hiHFhit_mc;
        if(var == "hiNpix")  currVar = hiNpix_mc;
        if(var == "hiNtracks")  currVar = hiNtracks_mc;
        if(var == "hiHFminusEta4")  currVar = hiHFminusEta4_mc;
        if(var == "hiHFplusEta4")  currVar = hiHFplusEta4_mc;
        if(var == "hiHFplus")  currVar = hiHFminus_mc;
        if(var == "hiHFminus")  currVar = hiHFplus_mc;
        if(num!=0){
            for(int i=0;i<num;++i){
                if(weightVar[i] == "hiET") currWeightVar[i] = hiET_mc;
                if(weightVar[i] == "hiEB")  currWeightVar[i] = hiEB_mc;
                if(weightVar[i] == "hiEE")  currWeightVar[i] = hiEE_mc;
                if(weightVar[i] == "hiHF")  currWeightVar[i] = hiHF_mc;
                if(weightVar[i] == "hiHFhit")  currWeightVar[i] = hiHFhit_mc;
                if(weightVar[i] == "hiNpix")  currWeightVar[i] = hiNpix_mc;
                if(weightVar[i] == "hiNtracks")  currWeightVar[i] = hiNtracks_mc;
                if(weightVar[i] == "hiHFminusEta4")  currWeightVar[i] = hiHFminusEta4_mc;
                if(weightVar[i] == "hiHFplusEta4")  currWeightVar[i] = hiHFplusEta4_mc;
                if(weightVar[i] == "hiHFminus")  currWeightVar[i] = hiHFminus_mc;
                if(weightVar[i] == "hiHFplus")  currWeightVar[i] = hiHFplus_mc;
            }
        }
        //denominator 
        h1D_mc_tot->Fill(currVar);

        double weight = 1.; 
        if(num!=0){
            double tempWeight[num];
            for(int i=0;i<num;++i){
                int bin = -1;
                double binVal = -1; 
                bin = weightHist[i]->FindBin(currWeightVar[i]);
                binVal = weightHist[i]->GetBinContent(bin);
                tempWeight[i] = ((binVal>0) ? (binVal) : 1.);
                weight = weight*tempWeight[i];
            }
           // if(var=="hiNtracks" && weightVar[0]=="hiHF" && weightVar[1]=="hiNpix") {
           //     if(iev==0) cout << " Variable = " << var << ", Weight Varible = " << weightVar[0] << endl; 
           //     cout << "hiNtracks = " << currVar 
           //         << "\thiHF = " << currWeightVar[0] << "\t weight_hiHF = " << tempWeight[0] 
           //         << "\thiNpix = "<< currWeightVar[1] << "\tweight_hiNpix = " << tempWeight[1]
           //         << "\ttotalWeight = " << weight << endl;
           // }
        }
        h1D_mc_tot_weighted->Fill(currVar,weight);
        Nden_mc = Nden_mc+weight;

        //Selected
        if(doNoTrigEff==0){
        if(!(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_mc ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_mc))
            continue;
        }

        if(!(pPAprimaryVertexFilter_mc==1 && pBeamScrapingFilter_mc==1 && phfCoincFilter1_mc==1 && pVertexFilterCutGplus_mc==1))
            continue;
        //numerator
        h1D_mc_evtSel->Fill(currVar);
        h1D_mc_evtSel_weighted->Fill(currVar,weight);
        Nnum_mc = Nnum_mc+weight;
    }// MC Event Loop for Efficiency 
    h1D_mc_filterRate->Divide(h1D_mc_evtSel_weighted, h1D_mc_tot_weighted, 1.,1., "B");

    // DATA Event Loop
    double Nden_data = 0; 
    double Nnum_data = 0; 
    int nEntries_data = tdata->GetEntries();
    for(int iev=0;iev<nEntries_data;++iev){
        tdata->GetEntry(iev);
        if ( iev % 100000 == 0 ) std::cout << "DATA :: " <<  iev << "/" << nEntries_data << " : " << iev*100.0/nEntries_data << std::endl;

        double currVar(0.);
        if(var == "hiET") currVar = hiET;
        if(var == "hiEB")  currVar = hiEB;
        if(var == "hiEE")  currVar = hiEE;
        if(var == "hiHF")  currVar = hiHF;
        if(var == "hiHFhit")  currVar = hiHFhit;
        if(var == "hiNpix")  currVar = hiNpix;
        if(var == "hiNtracks")  currVar = hiNtracks;
        if(var == "hiHFminus")  currVar = hiHFminus;
        if(var == "hiHFplus")  currVar = hiHFplus;
        if(var == "hiHFminusEta4")  currVar = hiHFminusEta4;
        if(var == "hiHFplusEta4")  currVar = hiHFplusEta4;

        h1D_data_tot->Fill(currVar);

        int bin = -1;
        double binVal = -1; 
        double weight = -1; 
        bin = h1D_mc_filterRate->FindBin(currVar);
        binVal = h1D_mc_filterRate->GetBinContent(bin);
        weight = ((binVal>0) ? (1./binVal) : 1.);

        Nden_data = Nden_data+weight;

        //Selected
        if(!(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1 ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2 ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3 ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4 ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5 ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6 ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7 ||
                    HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8))
            continue;
        if(!(pPAprimaryVertexFilter==1 && pBeamScrapingFilter==1 && phfCoincFilter1==1 && pVertexFilterCutGplus==1))
            continue;

        h1D_data_evtSel->Fill(currVar);
        h1D_data_corrected->Fill(currVar,weight);
        Nnum_data = Nnum_data+weight;
    }// DATA Event Loop for Efficiency 
    h1D_weightFactor->Divide(h1D_data_corrected, h1D_mc_tot_weighted, 1.,1., "B");
    //h1D_weightFactor->Divide(h1D_data_corrected, h1D_mc_tot, 1.,1., "B");


    res.h1D_mc_tot = (TH1D*) h1D_mc_tot->Clone(Form("h1D_%s",h1D_mc_tot->GetName()));
    res.h1D_mc_evtSel = (TH1D*) h1D_mc_evtSel->Clone(Form("h1D_%s",h1D_mc_evtSel->GetName()));
    res.h1D_mc_tot_weighted = (TH1D*) h1D_mc_tot_weighted->Clone(Form("h1D_%s",h1D_mc_tot_weighted->GetName()));
    res.h1D_mc_evtSel_weighted = (TH1D*) h1D_mc_evtSel_weighted->Clone(Form("h1D_%s",h1D_mc_evtSel_weighted->GetName()));
    res.h1D_mc_filterRate = (TH1D*) h1D_mc_filterRate->Clone(Form("h1D_%s",h1D_mc_filterRate->GetName()));
    res.h1D_data_tot = (TH1D*) h1D_data_tot->Clone(Form("h1D_%s",h1D_data_tot->GetName()));
    res.h1D_data_evtSel = (TH1D*) h1D_data_evtSel->Clone(Form("h1D_%s",h1D_data_evtSel->GetName()));
    res.h1D_data_corrected = (TH1D*) h1D_data_corrected->Clone(Form("h1D_%s",h1D_data_corrected->GetName()));
    res.h1D_weightFactor= (TH1D*) h1D_weightFactor->Clone(Form("h1D_%s",h1D_weightFactor->GetName()));
    res.MC_integFilterRate = (double)Nnum_mc/(double)Nden_mc; 
    res.DATA_integEvtSelEff = (double)(res.h1D_data_evtSel->Integral()/res.h1D_data_corrected->Integral());

    //cout << "res.DATA_integEvtSelEff = (double)(res.h1D_data_evtSel->Integral()/res.h1D_data_corrected->Integral())"
    //    << res.DATA_integEvtSelEff << " = " << res.h1D_data_evtSel->Integral() <<"/" << res.h1D_data_corrected->Integral();

    delete h1D_mc_tot;
    delete h1D_mc_evtSel;
    delete h1D_mc_tot_weighted;
    delete h1D_mc_evtSel_weighted;
    delete h1D_mc_filterRate;
    delete h1D_data_tot;
    delete h1D_data_evtSel;
    delete h1D_data_corrected;
    delete h1D_weightFactor;

    TFile* fout = new TFile(fout_name,"RECREATE");
    fout->cd();
    res.h1D_mc_tot->Write();
    res.h1D_mc_evtSel->Write();
    res.h1D_mc_tot_weighted->Write();
    res.h1D_mc_evtSel_weighted->Write();
    res.h1D_mc_filterRate->Write();
    res.h1D_data_tot->Write();
    res.h1D_data_evtSel->Write();
    res.h1D_data_corrected->Write();
    res.h1D_weightFactor->Write();
    fout->Close();
    
    delete fout;
    delete tdata;
    delete tskimanalysis;
    delete thltanalysis;
    delete tmc;
    delete tskimanalysismc;
    delete thltanalysismc;
    
    return res;
}
#endif
