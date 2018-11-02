// Created : 20 June 2017
// Author : Yeonju Go // 
#include "../../yjUtility.h"
const int colHere[]={8,4,2,kYellow+2,}; 
double getEffVal(TString fname);
void summary_EvtSelEff_all() { 
    TH1::SetDefaultSumw2(); 
    gStyle -> SetOptStat(0); 
   
    string mcName[]={"HIJING","EPOSLHC","AMPT"}; 
    const int Nwei = 6;
    TString weight[Nwei];
    //TString leg[Nwei] = {"0th","Ntracks_Npix","Ntracks_Npix_Ntracks_Npix","Ntracks_Npix_Ntracks_Npix_Ntracks_Npix"};
    TString var = "hiNpix";
    TString weight_var = "hiNtracks";
    TString temp_even = var+"_"+weight_var;
    TString temp_odd = weight_var;
    TString temp_pair = "_"+var+"_"+weight_var;
    for(int i=0; i<Nwei; ++i){
        if(i==0) { 
            weight[i]=""; 
        } else if(i%2==1) {
            if(i==1) { weight[i] = temp_odd;}
            else {
                temp_odd = temp_odd  + temp_pair;
                weight[i] = temp_odd; 
            }
        } else {
            if(i==2) { 
                weight[i] = temp_even;
            } else {
                temp_even = temp_even + temp_pair;
                weight[i] = temp_even; 
            }
        }
        cout << "weight " << i << " = " << weight[i] << endl;
    } 

    const int Nmc = 3;
    TH1D* h1D_ESE[Nmc];
    for(int imc=0; imc<Nmc; ++imc){
        h1D_ESE[imc] = new TH1D(Form("h1D_ESE_%s",mcName[imc].data()), Form("%s (iteration with %s);The number of iteration;Event Selection Efficiency",var.Data(),weight_var.Data()), Nwei, 0, Nwei);
        for(int i=0; i<Nwei; ++i){
            TString fname =""; 
            if(i==0) fname = Form("output/evtSelEff_foldingMethod_run286471_%s_%s.root",mcName[imc].data(),var.Data());
            else fname = Form("output/evtSelEff_foldingMethod_run286471_%s_%s_WeightedBy_%s.root",mcName[imc].data(),var.Data(),weight[i].Data());
            double effVal = getEffVal(fname);
            h1D_ESE[imc]->SetBinContent(i+1,effVal); 
        } 
    }

    // DRAW
    SetyjPadStyle(); 
    TCanvas* c1 = new TCanvas("c1","",500,500); 
    for(int imc=0; imc<Nmc; ++imc){
        SetHistTextSize(h1D_ESE[imc]);
        h1D_ESE[imc]->SetMarkerStyle(20);
        h1D_ESE[imc]->SetMarkerColor(colHere[imc]);
        h1D_ESE[imc]->GetYaxis()->SetRangeUser(0.90,1.03);
        if(imc==0) h1D_ESE[imc]->DrawCopy("p");
        else h1D_ESE[imc]->DrawCopy("p same");
    }
    jumSun(0,1,Nwei,1);
    
    TLegend* l1 = new TLegend(0.6,0.78,0.95,0.95);
    legStyle(l1);
    for(int imc=0; imc<Nmc; ++imc){
        l1->AddEntry(h1D_ESE[imc],Form("%s",mcName[imc].data()));
    }
    l1->Draw("same");
    
    c1->SaveAs(Form("figures/summary_ESE_%s_iteration_with_%s.pdf",var.Data(),weight_var.Data()));
}
double getEffVal(TString fname){
    TFile* fin = new TFile(fname);
    TH1D* h1D_evtSel = (TH1D*) fin->Get("h1D_data_evtSel");
    TH1D* h1D_corrected = (TH1D*) fin->Get("h1D_data_corrected");
    double Eff = h1D_evtSel->Integral()/h1D_corrected->Integral();
    delete fin;
    return Eff;
}
