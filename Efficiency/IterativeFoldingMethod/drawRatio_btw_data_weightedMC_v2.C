// Created : 05 July 2017
// Modified : 05 July 2017
// Author : Yeonju Go  
// This macro is the same with drawRatio_btw_data_weightedMC.C 
// except that each iterated MC is divided by each data corrected distribution in terms of ratio 

#include "../../yjUtility.h"
const int colHere[]={1,4,2,kYellow+2,8,41}; 
const int markerHere[]={24,25,26,27,28,1,1,1,1}; 
double getEffVal(TString fname);
void drawRatio(const int Nvar, string var, string wVar[], string mcName = "HIJING", bool isLogx=1, bool doNoTrigEff=1, bool doNorm=1); 

void drawRatio_btw_data_weightedMC_v2(bool doNoTrigEff=0,bool isLogx=1) {
    const int Nvar=5;
    string var[] = {"hiNtracks","hiHF","hiNpix","hiHFplusEta4","hiHFminusEta4"};
   
    for(int i=0; i<Nvar;++i){
    for(int j=0; j<Nvar;++j){
    for(int k=0; k<Nvar;++k){
    for(int l=0; l<Nvar;++l){
    for(int m=0; m<Nvar;++m){
        if(i==m ||  i==l || i==k || i==j) continue;
        if(j==m ||  j==l || j==k) continue;
        if(k==m ||  k==l) continue;
        if(l==m) continue; 
        string wVarArr[Nvar-1];
        wVarArr[0]=var[j];
        wVarArr[1]=var[k];
        wVarArr[2]=var[l];
        wVarArr[3]=var[m];
        drawRatio(Nvar, var[i], wVarArr,"HIJING", isLogx, doNoTrigEff);
        drawRatio(Nvar, var[i], wVarArr,"EPOSLHC", isLogx, doNoTrigEff);
        drawRatio(Nvar, var[i], wVarArr,"AMPT", isLogx, doNoTrigEff);
    }}}}}
}

void drawRatio(const int Nvar, string var, string wVar[], string mcName, bool isLogx, bool doNoTrigEff, bool doNorm) 
{ 
    TH1::SetDefaultSumw2(); 
    gStyle -> SetOptStat(0); 
    
    const int Nwei = 5;
    TString weight[Nwei];
    TString leg[Nwei];
    for(int i=0; i<Nvar; ++i){
        if(i==0) leg[i] = "0th";
        else if(i==1) leg[i] = wVar[i-1];
        else leg[i] = leg[i-1]+"_"+wVar[i-1];

        if(i==0) weight[i] = "";
        else if(i==1) weight[i] = wVar[i-1];
        else weight[i] = weight[i-1]+"_"+wVar[i-1];
    }

    TString noTrigEffCap = "";
    if(doNoTrigEff) noTrigEffCap = "_noTrigEff";   

    double eff[Nwei]={0};
    TFile* fin[Nwei];
    for(int i=0; i<Nwei; ++i){
        TString fname = "";
        // Get weighted MC and data corrected distributrions in root files produced from evtSelEff_folding_onlyMixing.C 
       // if(i==0) fname = Form("output/weight_v6_17June30/evtSelEff_foldingMethod%s_run286471_%s_%s.root",noTrigEffCap.Data(),mcName.data(),var.data());
       // else fname = Form("output/weight_v6_17June30/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s.root",noTrigEffCap.Data(),mcName.data(),var.data(),weight[i].Data());
        if(i==0) fname = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s.root",noTrigEffCap.Data(),mcName.data(),var.data());
        else fname = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s.root",noTrigEffCap.Data(),mcName.data(),var.data(),weight[i].Data());
        fin[i] = new TFile(fname);
        eff[i] = getEffVal(fname);
    }
    TH1D* h1D_weighted[Nwei];
    TH1D* h1D_data_evtSel[Nwei];
    TH1D* h1D_data_corrected[Nwei];
    for(int i=0; i<Nwei; ++i){
        h1D_weighted[i] = (TH1D*) fin[i]->Get("h1D_mc_tot_weighted");
        h1D_weighted[i]->SetName(Form("h1D_weighted_%d",i));
        h1D_data_evtSel[i] = (TH1D*) fin[i]->Get("h1D_data_evtSel");
        h1D_data_evtSel[i]->SetName(Form("h1D_data_evtSel_%d",i));
        h1D_data_corrected[i] = (TH1D*) fin[i]->Get("h1D_data_corrected");
        h1D_data_corrected[i]->SetName(Form("h1D_data_corrected_%d",i));
    }

    // DRAW
    SetyjPadStyle(); 
    TCanvas* c1 = new TCanvas("c1","",400,600); 
    ratioPanelCanvas(c1);
    c1->cd(1);
    if(isLogx) gPad->SetLogx();
    SetHistTextSize(h1D_weighted[0],0.9);
    SetHistTextSize(h1D_data_corrected[0],0.9);
    for(int i=0; i<Nwei; ++i){
        h1D_weighted[i]->Scale(1.,"width");
        h1D_data_corrected[i]->Scale(1.,"width");
        if(doNorm) h1D_weighted[i]->Scale(1./h1D_weighted[i]->Integral());
        if(doNorm) h1D_data_corrected[i]->Scale(1./h1D_data_corrected[i]->Integral());
        h1D_weighted[i]->SetLineColor(colHere[i]);
        h1D_data_corrected[i]->SetMarkerStyle(markerHere[i]);
        //h1D_data_corrected[i]->SetMarkerStyle(27);
        h1D_data_corrected[i]->SetMarkerColor(colHere[i]);
    }
    
    //calculate the maximum range
    double maxY = 0; 
    for(int i=0; i<Nwei; ++i){
        double tempMaxY = cleverRange(h1D_weighted[i]);
        double tempMaxY_data = cleverRange(h1D_data_corrected[i]);
        if(tempMaxY_data > tempMaxY) tempMaxY = tempMaxY_data;
        if(tempMaxY > maxY) maxY = tempMaxY;
    }
    maxY = 1.2*maxY;
    for(int i=0; i<Nwei; ++i) {
        h1D_data_corrected[i]->GetYaxis()->SetRangeUser(0,maxY);
        if(i==0) h1D_data_corrected[i]->DrawCopy("p");
        else h1D_data_corrected[i]->DrawCopy("p same");
        h1D_weighted[i]->GetYaxis()->SetRangeUser(0,maxY);
        h1D_weighted[i]->DrawCopy("hist same");
    }

    TLegend* l1 = new TLegend(0.3,0.45,0.95,0.95);
    //TLegend* l1 = new TLegend(0.2,0.7,0.95,0.95);
    legStyle(l1);
    for(int i=0; i<Nwei; ++i){
        l1->AddEntry(h1D_data_corrected[i],Form("Data corrected;iter%d",i));
    }
    for(int i=0; i<Nwei; ++i){
        l1->AddEntry(h1D_weighted[i],Form("MC weighted;iter%d;Eff=%.3f",i,eff[i]));
       // if(i==0) l1->AddEntry(h1D_weighted[i],Form("MC;Eff=%.3f",eff[i]));
       // else l1->AddEntry(h1D_weighted[i],Form("MC weightBy_%s;Eff=%.3f",leg[i].Data(),eff[i]));
    }
    l1->Draw("same");


    c1->cd(2);
    if(isLogx) gPad->SetLogx();
    TH1D* ratio[Nwei];
    for(int i=0; i<Nwei; ++i){
        ratio[i] = (TH1D*) h1D_weighted[i]->Clone(Form("ratio_%d",i));
        ratio[i]->Divide(h1D_weighted[i], h1D_data_corrected[i], 1., 1.);
        SetHistTextSize(ratio[i],0.9);
        ratio[i]->SetMarkerSize(0.8);
        ratio[i]->SetMarkerStyle(20);
        ratio[i]->SetMarkerColor(colHere[i]);
        ratio[i]->SetMarkerStyle(markerHere[i]);
        ratio[i]->SetLineColor(colHere[i]);
        ratio[i]->SetLineStyle(1);
        ratio[i]->SetTitle(Form(";%s;MC/DATA",var.data()));
        ratio[i]->GetYaxis()->SetRangeUser(0,3);
        if(i==0) ratio[i]->Draw("p");
        else ratio[i]->Draw("p same");
    }
    drawText(Form("weightedBy_%s",weight[Nvar-1].Data()),0.2,0.9,0,kBlack,10); 
    double xmax = 350;
    if(var=="hiNtracks") xmax =300; 
    if(var=="hiNpix") xmax =1800; 
    jumSun(0,1,xmax,1);
   
   // TLegend* l2 = new TLegend(0.6,0.7,0.95,0.95);
   // legStyle(l2);
   // for(int i=0; i<Nwei; ++i){
   //     l2->AddEntry(ratio[i],Form("%s/data",leg[i].Data()),"pl");
   // }
   // l2->Draw("same");

    if(doNorm) c1->SaveAs(Form("figures/Ratio_MCiteratedOverData%s_%s_%s_weightedBy_%s_isLogx%d.pdf",noTrigEffCap.Data(),mcName.data(),var.data(),weight[Nvar-1].Data(),(int)isLogx));
    else c1->SaveAs(Form("figures/Ratio_MCiteratedOverData%s_noNorm_%s_%s_weightedBy_%s_isLogx%d.pdf",noTrigEffCap.Data(),mcName.data(),var.data(),weight[Nvar-1].Data(),(int)isLogx));

    delete c1;
    delete l1;
    for(int i=0; i<Nwei; ++i){
        delete ratio[i]; 
        delete h1D_weighted[i]; 
        delete h1D_data_evtSel[i]; 
        delete h1D_data_corrected[i]; 
        delete fin[i]; 
    }
}

double getEffVal(TString fname){
    TFile* fin = new TFile(fname);
    TH1D* h1D_evtSel = (TH1D*) fin->Get("h1D_data_evtSel");
    TH1D* h1D_corrected = (TH1D*) fin->Get("h1D_data_corrected");
    double Eff = h1D_evtSel->Integral()/h1D_corrected->Integral();
    delete fin;
    return Eff;
}
