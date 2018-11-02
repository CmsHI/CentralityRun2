// Created : 20 June 2017
// Modified : 23 June 2017
// Author : Yeonju Go // 
#include "../../yjUtility.h"
const int colHere[]={1,4,2,kYellow+2,8,41}; 
const int lineHere[]={1,1,1,1,1,1}; 
const int markerHere[]={24,25,26,27,28,1,1,1,1}; 
double getEffVal(TString fname);
void drawRatioBtwWeightedDist_data(bool doNoTrigEff=1, string mcName = "HIJING", TString var="hiHF", string wVar1="hiNpix", string wVar2="hiNtracks", string wVar3="hiHFplusEta4", string wVar4="hiHFminusEta4", bool doNorm=1) { 
    TH1::SetDefaultSumw2(); 
    gStyle -> SetOptStat(0); 
    
    const int Nwei = 5;
    TString weight[Nwei];
    TString leg[Nwei];
   // for(int i=0; i<Nvar; ++i){
   //     if(i==0) leg[i] = "0th";
   //     else if(i==1) leg[i] = wVar[i-1];
   //     else leg[i] = leg[i-1]+"_"+wVar[i-1];

   //     if(i==0) weight[i] = "";
   //     else if(i==1) weight[i] = wVar[i-1];
   //     else weight[i] = weight[i-1]+"_"+wVar[i-1];
   // }
    leg[0] = "0th";
    leg[1] = wVar1;
    leg[2] = wVar1+"_"+wVar2;
    leg[3] = wVar1+"_"+wVar2+"_"+wVar3;
    leg[4] = wVar1+"_"+wVar2+"_"+wVar3+"_"+wVar4;
    weight[0] = "";
    weight[1] = wVar1;
    weight[2] = wVar1+"_"+wVar2;
    weight[3] = wVar1+"_"+wVar2+"_"+wVar3;
    weight[4] = wVar1+"_"+wVar2+"_"+wVar3+"_"+wVar4;

    TString noTrigEffCap = "";
    if(doNoTrigEff) noTrigEffCap = "_noTrigEff";   

    double eff[Nwei]={0};
    TFile* fin[Nwei];
    for(int i=0; i<Nwei; ++i){
        TString fname = "";
        if(i==0) fname = Form("output/weight_v6_17June30/evtSelEff_foldingMethod%s_run286471_%s_%s.root",noTrigEffCap.Data(),mcName.data(),var.Data());
        else fname = Form("output/weight_v6_17June30/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s.root",noTrigEffCap.Data(),mcName.data(),var.Data(),weight[i].Data());
       // if(i==0) fname = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s.root",noTrigEffCap.Data(),mcName.data(),var.Data());
       // else fname = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s.root",noTrigEffCap.Data(),mcName.data(),var.Data(),weight[i].Data());
        fin[i] = new TFile(fname);
        eff[i] = getEffVal(fname);
    }
    TH1D* h1D_weighted[Nwei];
    TH1D* h1D_data_evtSel = (TH1D*) fin[Nwei-1]->Get("h1D_data_evtSel");
    TH1D* h1D_data_corrected = (TH1D*) fin[Nwei-1]->Get("h1D_data_corrected");
    for(int i=0; i<Nwei; ++i){
        h1D_weighted[i] = (TH1D*) fin[i]->Get("h1D_mc_tot_weighted");
        h1D_weighted[i]->SetName(Form("h1D_weighted_%d",i));
    }

    // DRAW
    SetyjPadStyle(); 
    TCanvas* c1 = new TCanvas("c1","",400,600); 
    ratioPanelCanvas(c1);
    c1->cd(1);
    SetHistTextSize(h1D_weighted[0],0.9);
    SetHistTextSize(h1D_data_corrected,0.9);
    h1D_data_corrected->Scale(1.,"width");
    if(doNorm) h1D_data_corrected->Scale(1./h1D_data_corrected->Integral());
    h1D_data_corrected->SetMarkerStyle(27);
    h1D_data_corrected->SetMarkerColor(1);
    for(int i=0; i<Nwei; ++i){
        h1D_weighted[i]->Scale(1.,"width");
        if(doNorm) h1D_weighted[i]->Scale(1./h1D_weighted[i]->Integral());
        h1D_weighted[i]->SetLineStyle(lineHere[i]);
        h1D_weighted[i]->SetLineColor(colHere[i]);
    }
    //calculate the maximum range
    double maxY = cleverRange(h1D_data_corrected);
    for(int i=0; i<Nwei; ++i){
        double tempMaxY = cleverRange(h1D_weighted[i]);
        if(tempMaxY > tempMaxY) maxY = tempMaxY;
    }
    h1D_data_corrected->GetYaxis()->SetRangeUser(0,maxY);
    h1D_data_corrected->DrawCopy("p");
    for(int i=0; i<Nwei; ++i) {
        h1D_weighted[i]->GetYaxis()->SetRangeUser(0,maxY);
        h1D_weighted[i]->DrawCopy("hist same");
    }

    TLegend* l1 = new TLegend(0.2,0.7,0.95,0.95);
    legStyle(l1);
    l1->AddEntry(h1D_data_corrected,"Data corrected");
    for(int i=0; i<Nwei; ++i){
        if(i==0) l1->AddEntry(h1D_weighted[i],Form("MC;Eff=%.3f",eff[i]));
        else l1->AddEntry(h1D_weighted[i],Form("MC weightBy_%s;Eff=%.3f",leg[i].Data(),eff[i]));
    }
    l1->Draw("same");
    c1->SetLogx();


    c1->cd(2);
    TH1D* ratio[Nwei];
    for(int i=0; i<Nwei; ++i){
        ratio[i] = (TH1D*) h1D_weighted[i]->Clone(Form("ratio_%d",i));
        ratio[i]->Divide(h1D_weighted[i], h1D_data_corrected, 1., 1.);
        SetHistTextSize(ratio[i],0.9);
        ratio[i]->SetMarkerSize(0.8);
        ratio[i]->SetMarkerStyle(20);
        ratio[i]->SetMarkerColor(colHere[i]);
        ratio[i]->SetMarkerStyle(markerHere[i]);
        ratio[i]->SetLineColor(colHere[i]);
        ratio[i]->SetLineStyle(1);
        ratio[i]->SetTitle(Form(";%s;Ratio",var.Data()));
        ratio[i]->GetYaxis()->SetRangeUser(0,2);
        if(i==0) ratio[i]->Draw("p");
        else ratio[i]->Draw("p same");
    }
    
    double xmax = 350;
    if(var=="hiNtracks") xmax =300; 
    if(var=="hiNpix") xmax =1800; 
    jumSun(0,1,xmax,1);
   
    TLegend* l2 = new TLegend(0.6,0.7,0.95,0.95);
    legStyle(l2);
    for(int i=0; i<Nwei; ++i){
        l2->AddEntry(ratio[i],Form("%s/data",leg[i].Data()),"pl");
    }
    l2->Draw("same");
    c1->SetLogx();
    if(doNorm) c1->SaveAs(Form("figures/Ratio_MCiteratedOverData%s_%s_%s_weightedBy_%s_%s.pdf",noTrigEffCap.Data(),mcName.data(),var.Data(),wVar1.data(),wVar2.data()));
    else c1->SaveAs(Form("figures/Ratio_MCiteratedOverData%s_noNorm_%s_%s_weightedBy_%s_%s.pdf",noTrigEffCap.Data(),mcName.data(),var.Data(),wVar1.data(),wVar2.data()));
}

double getEffVal(TString fname){
    TFile* fin = new TFile(fname);
    TH1D* h1D_evtSel = (TH1D*) fin->Get("h1D_data_evtSel");
    TH1D* h1D_corrected = (TH1D*) fin->Get("h1D_data_corrected");
    double Eff = h1D_evtSel->Integral()/h1D_corrected->Integral();
    delete fin;
    return Eff;
}
