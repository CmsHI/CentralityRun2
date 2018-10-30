// to get the trigger turn-on curves.
// baseTrig goes to denominator.
// have to be careful of prescale part!!!
// 
// Author Yeonju Go 
// 19 Nov. 2015
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
#include "../../utils.h"

void triggerTurnOn(const char* fname="root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/ExpressPhysics/Merged/ExpressHiForest_run262163-262172_1.4M.root",
        ///store/group/phys_heavyions/velicanu/forest/Run2015E/ExpressPhysics/Merged/ExpressHiForest_run262163_277k.root", 
        const char* baseTrig = "L1_MinimumBiasHF1_OR", 
        const char* xvar="hiHF",
        double xmin = 0,
        double xmax = 200,
        int nBin = 50,
        TString cap="")
{

    //baseTrig : HLT_ZeroBias_v2 53450 , L1Tech_BPTX_plus_AND_minus.v0 224871, L1_ZeroBias 54037,  
    //trig : HLT_L1MinimumBiasHF1OR_part1_v1, HLT_L1MinimumBiasHF2OR_v1, HLT_L1MinimumBiasHF1AND_v1, HLT_L1MinimumBiasHF2AND_v1
    const int nTrig = 3;
//    TString trig[] = {"HLT_L1MinimumBiasHF1AND_v1","HLT_L1MinimumBiasHF2OR_v1","HLT_L1MinimumBiasHF2AND_v1"};
    TString trig[] = {"L1_MinimumBiasHF1_AND","L1_MinimumBiasHF2_OR","L1_MinimumBiasHF2_AND"};

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    SetHistTitleStyle();
    SetyjPadStyle();
    TFile *fin = TFile::Open(fname);
    cout << "open file : " << fname << endl;
    TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
    TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
    TTree *t = (TTree*) fin -> Get("hltanalysis/HltTree");
    t-> AddFriend(t_evt);
    t-> AddFriend(t_skim);

    TH1D* h[nTrig];
    TH1D* hbase = new TH1D(Form("hbase"),Form("%s;%s;Entries",baseTrig,xvar),nBin,xmin,xmax);
    const char* lumiCut = "((1==1))";
    //const char* lumiCut = "((Run>262163))";
    //const char* lumiCut = "((Run==262163 && LumiBlock>5) || (Run>262163))";
    //t->Draw( Form("%s>>hbase",xvar), Form("%s_Prescl*(%s==1 && %s)",baseTrig,baseTrig,lumiCut) );
    //t->Draw( Form("%s>>hbase",xvar), Form("(%s==1 && %s)",baseTrig,lumiCut) );
    //t->Draw( Form("%s>>hbase",xvar), Form("(%s==1 && %s)",baseTrig,lumiCut) );
    int baseTrigEntry = t->GetEntries(Form("%s==1 && %s",baseTrig,lumiCut));
    int trigEntry[nTrig];
    cout << "BASE TRIGGER ::: " << baseTrig << " ::: " << endl; 
    for(int i=0;i<nTrig;i++){
        //const char* tmpTrig = Form("%s_Prescl*(%s==1 && %s)",trig[i].Data(),trig[i].Data(), lumiCut);
        //const char* tmpTrig = Form("%s_Prescl*(%s==1 && %s==1 && %s)",trig[i].Data(),trig[i].Data(), baseTrig, lumiCut);
        const char* tmpTrig = Form("(%s==1 && %s==1 && %s)",trig[i].Data(), baseTrig, lumiCut);
        h[i] = new TH1D(Form("h%d",i), Form("%s;%s;Entries",trig[i].Data(),xvar),nBin,xmin,xmax);
        t->Draw( Form("%s>>%s",xvar,h[i]->GetName()), tmpTrig );
        trigEntry[i]=t->GetEntries(tmpTrig);
        cout << trig[i] << "::: efficiency = " << trigEntry[i] << "/" << baseTrigEntry << " = " << (double)trigEntry[i]/(double)baseTrigEntry*100 <<" %" << endl;
    }

    TLegend* l1 = new TLegend(0.45, 0.7, 0.9, 0.95);
    legStyle(l1);
   
    TCanvas* c_tot = new TCanvas(Form("c_tot"), "c_tot", 500,1000);
    c_tot->Divide(1,2);
    c_tot->cd(1);
    gPad->SetLogy();
    //hbase->GetYaxis()->SetRangeUser(0.8,hbase->GetMaximum()*10);
    hbase->DrawCopy("hist");
    l1->AddEntry(hbase, Form("%s",baseTrig),"l");
    for(int i=0;i<nTrig;i++){
        c_tot->cd(1);
        SetHistColor(h[i],col[i+2]);
        l1->AddEntry(h[i], Form("%s",trig[i].Data()),"l");
       // h[i]->GetYaxis()->SetRangeUser(0.8,hbase->GetMaximum()*10);
        h[i]->DrawCopy("hist same");
    }
    l1->Draw("same");

    //Efficiency

    TH1D* hEff[nTrig];
    for(int i=0;i<nTrig;i++){
        hEff[i] = (TH1D*) h[i]->Clone(Form("hEff%d",i));
        hEff[i]->SetTitle(Form(";%s;Firing Rate",xvar));
        c_tot->cd(2);
        hEff[i]->Divide(h[i], hbase, 1,1,"B");  
        hEff[i]->GetYaxis()->SetRangeUser(0,1.2);
        if(i==0) hEff[i]->Draw("pe");
        else hEff[i]->Draw("pe same");
    }
    jumSun(xmin,1,xmax,1);
    c_tot->SaveAs(Form("pdf/trigTurnOn_%s_%s%s.png",xvar,baseTrig,cap.Data()));
}

int main(){
//    const char* fin = "root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/ExpressPhysics/Merged/ExpressHiForest_run262163_277k.root";// 5 < LumiBlock < 242 
//    const char* fin = "root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/ExpressPhysics/Merged/ExpressHiForest_run262167_473k.root";// 1 <LumiBlock < 12
    const char* fin = "root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/ExpressPhysics/Merged/ExpressHiForest_run262163-262172_1.4M.root";//all forest

    //baseTrig : HLT_ZeroBias_v2 53450 , L1Tech_BPTX_plus_AND_minus.v0 224871, L1_ZeroBias 54037,  
/*
    triggerTurnOn(fin,"HLT_ZeroBias_v2" ,"hiHF", 0, 200, 50);
    triggerTurnOn(fin,"HLT_ZeroBias_v2" ,"hiNpix", 0, 1400, 50);
    //triggerTurnOn(fin,"HLT_ZeroBias_v2" ,"hiBin", 0, 200, 50);

    triggerTurnOn(fin,"L1Tech_BPTX_plus_AND_minus.v0" ,"hiHF", 0, 200, 50);
    triggerTurnOn(fin,"L1Tech_BPTX_plus_AND_minus.v0" ,"hiNpix", 0, 1400, 50);
    //triggerTurnOn(fin,"L1Tech_BPTX_plus_AND_minus.v0" ,"hiBin", 0, 200, 50);

    triggerTurnOn(fin,"L1_ZeroBias" ,"hiHF", 0, 200, 50);
    triggerTurnOn(fin,"L1_ZeroBias" ,"hiNpix", 0, 1400, 50);
    //triggerTurnOn(fin,"L1_ZeroBias" ,"hiBin", 0, 200, 50);
*/     
    triggerTurnOn(fin,"L1_MinimumBiasHF1_OR" ,"hiHF", 0, 200, 50);
   // triggerTurnOn(fin,"L1_MinimumBiasHF1_OR" ,"hiNpix", 0, 1500, 50);
//    triggerTurnOn(fin,"L1_MinimumBiasHF1_OR" ,"hiNtracksPtCut", 0, 100, 50);
   
    //triggerTurnOn(fin,"L1_MinimumBiasHF1_OR" ,"hiHF", 0, 30, 10, "_zoom");
   // triggerTurnOn(fin,"L1_MinimumBiasHF1_OR" ,"hiNpix", 0, 200, 50, "_zoom");
//    triggerTurnOn(fin,"L1_MinimumBiasHF1_OR" ,"hiNtracksPtCut", 0, 20, 20, "_zoom");

    /*
    triggerTurnOn(fin,"HLT_L1MinimumBiasHF1OR_part1_v1","hiHF", 0, 200, 50);
    triggerTurnOn(fin,"HLT_L1MinimumBiasHF1OR_part1_v1" ,"hiNpix", 0, 1500, 50);
    triggerTurnOn(fin,"HLT_L1MinimumBiasHF1OR_part1_v1" ,"hiNtracksPtCut", 0, 100, 50);
   
    triggerTurnOn(fin,"HLT_L1MinimumBiasHF1OR_part1_v1","hiHF", 0, 20, 50, "_zoom");
    triggerTurnOn(fin,"HLT_L1MinimumBiasHF1OR_part1_v1" ,"hiNpix", 0, 100, 50, "_zoom");
    triggerTurnOn(fin,"HLT_L1MinimumBiasHF1OR_part1_v1" ,"hiNtracksPtCut", 0, 50, 50, "_zoom");
*/
    return 0;
}
