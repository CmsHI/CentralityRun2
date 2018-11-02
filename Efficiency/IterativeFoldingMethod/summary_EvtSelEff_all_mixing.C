// Created : 20 June 2017
// Author : Yeonju Go // 
#include "../../yjUtility.h"
#include  <THStack.h>
const int Nmc=3;//HIJING EPOS AMPT
string mcName[]={"HIJING","EPOSLHC","AMPT"}; 

class EFFSET{
    public:
       float eff[Nmc];
       float eff_withTrig[Nmc];
};

const int colHere[]={1,4,2,8,4,2,kYellow+2}; 
const int fstHere[]={3004,3005,3012,3006}; 
double getEffVal(TString fname);
EFFSET summary_EvtSelEff(bool doNoTrigEff, const int nth, const int Nvar, string var, string wVar[]);

void summary_EvtSelEff_all_mixing(bool doNoTrigEff=0, const int nth=5) { 
    const int Nvar=5;
    string var[] = {"hiNtracks","hiHF","hiNpix","hiHFplusEta4", "hiHFminusEta4"}; 
    //string var[] = {"hiNtracks","hiHF","hiNpix","hiHFplus", "hiHFminus"}; 
    string var_tot = "";
    string var_tot_forFileName = "";
    for(int i=0; i<Nvar;++i){
        if(i==0) var_tot = var[i];
        else var_tot = var_tot + ", "+var[i];
        if(i==0) var_tot_forFileName = var[i];
        else var_tot_forFileName = var_tot_forFileName + "_"+var[i];
    }
    TString noTrigEffCap = "";
    if(doNoTrigEff) noTrigEffCap = "_noTrigEff"; 

    // HISTOGRAM DEFINE
    TH1D* h1D_avg = new TH1D("h1D_avg",";Event Selection Efficiency;Entries",20,0.90,1.00);
    TH1D* h1D_avg_trig = new TH1D("h1D_avg_trig",";(Trigger + Event Selection) Efficiency;Entries",20,0.90,1.00);
    TH1D* h1D_avg_eachVar[Nvar];
    TH1D* h1D_avg_trig_eachVar[Nvar];
    TH1D* h1D_avg_eachMC[Nmc];
    TH1D* h1D_avg_trig_eachMC[Nmc];
    TH1D* h1D_avg_eachMC_eachVar[Nmc][Nvar];
    TH1D* h1D_avg_trig_eachMC_eachVar[Nmc][Nvar];

    THStack *hs = new THStack("hs","");
    THStack *hs_trig = new THStack("hs_trig","");
    THStack *hs_eachVar[Nvar];
    THStack *hs_trig_eachVar[Nvar];
    for(int imc=0; imc<Nmc; ++imc){
        h1D_avg_eachMC[imc] = (TH1D*) h1D_avg->Clone(Form("%s_%s",h1D_avg->GetName(),mcName[imc].data()));
        h1D_avg_trig_eachMC[imc] = (TH1D*) h1D_avg_trig->Clone(Form("%s_%s",h1D_avg_trig->GetName(),mcName[imc].data()));
        for(int i=0; i<Nvar; ++i){
            double xmin=0.90;
            int xbin=30;
            if(var[i]=="hiHFhit"){ xmin=0.80; xbin=40; }
            if(imc==0) {
                hs_eachVar[i] = new THStack(Form("hs_%s",var[i].data()),"");
                hs_trig_eachVar[i] = new THStack(Form("hs_trig_%s",var[i].data()),"");
                h1D_avg_eachVar[i] = new TH1D(Form("h1D_avg_%s",var[i].data()),";Event Selection Efficiency;Entries",xbin,xmin,1.00);
                h1D_avg_trig_eachVar[i] = new TH1D(Form("h1D_avg_trig_%s",var[i].data()),";(Trigger + Event Selection) Efficiency;Entries",xbin,xmin,1.00);
            }
            h1D_avg_eachMC_eachVar[imc][i] = (TH1D*) h1D_avg_eachVar[i]->Clone(Form("%s_%s_%s",h1D_avg->GetName(),mcName[imc].data(),var[i].data()));
            h1D_avg_trig_eachMC_eachVar[imc][i] = (TH1D*) h1D_avg_trig_eachVar[i]->Clone(Form("%s_%s_%s",h1D_avg_trig->GetName(),mcName[imc].data(),var[i].data()));
        }
    }
    double eff_mean_trig[Nmc]={0.0};  
    double eff_mean[Nmc]={0.0};  
    int nCount = 0;
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
        EFFSET a = summary_EvtSelEff(doNoTrigEff, 3, Nvar, var[i], wVarArr);
        for(int imc=0; imc<Nmc; ++imc){
            h1D_avg->Fill(a.eff[imc]); 
            h1D_avg_trig->Fill(a.eff_withTrig[imc]); 
            h1D_avg_eachMC[imc]->Fill(a.eff[imc]);
            h1D_avg_trig_eachMC[imc]->Fill(a.eff_withTrig[imc]);
            h1D_avg_eachVar[i]->Fill(a.eff[imc]); 
            h1D_avg_trig_eachVar[i]->Fill(a.eff_withTrig[imc]); 
            h1D_avg_eachMC_eachVar[imc][i]->Fill(a.eff[imc]);
            h1D_avg_trig_eachMC_eachVar[imc][i]->Fill(a.eff_withTrig[imc]);
            nCount++;
            eff_mean_trig[imc] += a.eff_withTrig[imc];
            eff_mean[imc] += a.eff[imc];
        }

    }}}}}
    double tot_eff_mean_trig=0.0; 
    double tot_eff_mean=0.0; 
    for(int imc=0; imc<Nmc; ++imc){
        tot_eff_mean_trig += eff_mean_trig[imc]; 
        tot_eff_mean += eff_mean[imc]; 
        eff_mean_trig[imc] = eff_mean_trig[imc]/(double)nCount*3.0;
        eff_mean[imc] = eff_mean[imc]/(double)nCount*3.0;
        //eff_mean[imc] = eff_mean[imc]/(double)nCount/3.0;
    }
    tot_eff_mean_trig = tot_eff_mean_trig/(double)nCount;
    tot_eff_mean = tot_eff_mean/(double)nCount;
    float xpos = 0.2;
    float ypos = 0.91;
    float dy = 0.035;
    //DRAW EVT. SEL. EFFICIENCY
    TLegend* ll = new TLegend(0.2,0.67,0.55,0.84); 
    legStyle(ll);
    TLegend* ll_trig = new TLegend(0.2,0.67,0.55,0.84); 
    legStyle(ll_trig);
    SetyjPadStyle(); 
    TCanvas* cc = new TCanvas("cc", "", 500,500); 
    SetHistTextSize(h1D_avg);
    h1D_avg->Draw();
    for(int imc=0; imc<Nmc; ++imc){
        h1D_avg_eachMC[imc]->SetFillColor(colHere[imc+1]);
        h1D_avg_eachMC[imc]->SetFillStyle(fstHere[imc]);
        ll->AddEntry(h1D_avg_eachMC[imc],Form("%s; %.2f %s",mcName[imc].data(),eff_mean[imc]*100,"%"));
        //hs->Add(h1D_avg_eachMC[imc]);
        if(imc==0) h1D_avg_eachMC[imc]->Draw();
        else h1D_avg_eachMC[imc]->Draw("same");
    }
    //hs->Draw("same"); 
    ll->Draw("same");
    drawText(Form("%dth iteration with : ",Nvar-1),xpos,ypos);
    drawText(Form("%s",var_tot.data()),xpos,ypos-dy);
    float mean = h1D_avg->GetMean();
    float rms = h1D_avg->GetRMS();
    drawText(Form("Mean : %.2f %s",tot_eff_mean*100, "%"),xpos,ypos-0.5,0,kBlack,20);
    //drawText(Form("(%.2f #pm %.2f) %s",mean*100, rms*100, "%"),xpos,ypos-0.5,0,kBlack,20);
     
    //DRAW (TRIGGER+EVT.SEL.) EFFICIENCY
    TCanvas* cc2 = new TCanvas("cc2", "", 500,500); 
    SetHistTextSize(h1D_avg_trig);
    h1D_avg_trig->Draw();
    for(int imc=0; imc<Nmc; ++imc){
        h1D_avg_trig_eachMC[imc]->SetFillColor(colHere[imc+1]);
        h1D_avg_trig_eachMC[imc]->SetFillStyle(fstHere[imc]);
        ll_trig->AddEntry(h1D_avg_trig_eachMC[imc],Form("%s; %.2f %s",mcName[imc].data(),eff_mean_trig[imc]*100,"%"));
        if(imc==0) h1D_avg_trig_eachMC[imc]->Draw();
        else h1D_avg_trig_eachMC[imc]->Draw("same");
        //hs_trig->Add(h1D_avg_trig_eachMC[imc]);
    }
    //hs_trig->Draw("same"); 
    ll_trig->Draw("same"); 
    drawText(Form("%dth iteration with : ",Nvar-1),xpos,ypos);
    drawText(Form("%s",var_tot.data()),xpos,ypos-dy);
    mean = h1D_avg_trig->GetMean();
    rms = h1D_avg_trig->GetRMS();
    drawText(Form("Mean : %.2f %s",tot_eff_mean_trig*100, "%"),xpos,ypos-0.5,0,kBlack,20);
    //drawText(Form("Mean : %.2f %s",tot_eff_mean_trig*100, "%"),xpos,ypos-0.5,0,kBlack,20);
    //drawText(Form("(%.2f #pm %.2f) %s",mean*100, rms*100, "%"),xpos,ypos-0.5,0,kBlack,20);
   
    //DRAW (EVT.SEL.) EFFICIENCY FOR EACH VARIABLE
    TCanvas* c_ese[Nvar];
    for(int i=0; i<Nvar; ++i){
        c_ese[i] = new TCanvas(Form("c_ese_%s",var[i].data()), "", 500,500); 
        for(int imc=0; imc<Nmc; ++imc){
            SetHistTextSize(h1D_avg_eachVar[i]);
            h1D_avg_eachMC_eachVar[imc][i]->SetFillColor(colHere[imc+1]);
            hs_eachVar[i]->Add(h1D_avg_eachMC_eachVar[imc][i]);
        }  
        h1D_avg_eachVar[i]->Draw();
        hs_eachVar[i]->Draw("same"); 
        ll->Draw("same"); 
        drawText(Form("%dth iteration with : ",Nvar-1),xpos,ypos);
        drawText(Form("%s",var_tot.data()),xpos,ypos-dy);
        mean = h1D_avg_eachVar[i]->GetMean();
        rms = h1D_avg_eachVar[i]->GetRMS();
        drawText(Form("Only from %s",var[i].data()),xpos,ypos-0.5+1.5*dy,0,kBlack,17);
        drawText(Form("(%.2f #pm %.2f) %s",mean*100, rms*100, "%"),xpos,ypos-0.5,0,kBlack,20);
        c_ese[i]->SaveAs(Form("figures/summary_ESE%s_onlyByVar_%s__%s.pdf",noTrigEffCap.Data(),var[i].data(),var_tot_forFileName.data()));
    }
    
    //DRAW (TRIGGER+EVT.SEL.) EFFICIENCY FOR EACH VARIABLE
    TCanvas* c_ese_trig[Nvar];
    for(int i=0; i<Nvar; ++i){
        c_ese_trig[i] = new TCanvas(Form("c_ese_trig_%s",var[i].data()), "", 500,500); 
        for(int imc=0; imc<Nmc; ++imc){
            SetHistTextSize(h1D_avg_trig_eachVar[i]);
            h1D_avg_trig_eachMC_eachVar[imc][i]->SetFillColor(colHere[imc+1]);
            hs_trig_eachVar[i]->Add(h1D_avg_trig_eachMC_eachVar[imc][i]);
        }  
        h1D_avg_trig_eachVar[i]->Draw();
        hs_trig_eachVar[i]->Draw("same"); 
        ll->Draw("same"); 
        drawText(Form("%dth iteration with : ",Nvar-1),xpos,ypos);
        drawText(Form("%s",var_tot.data()),xpos,ypos-dy);
        mean = h1D_avg_trig_eachVar[i]->GetMean();
        rms = h1D_avg_trig_eachVar[i]->GetRMS();
        drawText(Form("Only from %s",var[i].data()),xpos,ypos-0.5+1.5*dy,0,kBlack,17);
        drawText(Form("(%.2f #pm %.2f) %s",mean*100, rms*100, "%"),xpos,ypos-0.5,0,kBlack,20);
        c_ese_trig[i]->SaveAs(Form("figures/summary_ESE%s_withTrigger_onlyByVar_%s__%s.pdf",noTrigEffCap.Data(),var[i].data(),var_tot_forFileName.data()));
    }
    
    //SAVE CANVAS  
    cc->SaveAs(Form("figures/summary_ESE%s_total_%s.pdf",noTrigEffCap.Data(),var_tot_forFileName.data()));
    cc2->SaveAs(Form("figures/summary_ESE%s_total_withTrigger_%s.pdf",noTrigEffCap.Data(),var_tot_forFileName.data()));
}    

EFFSET summary_EvtSelEff(bool doNoTrigEff, const int nth, const int Nvar, string var, string wVar[]) { 
    TH1::SetDefaultSumw2(); 
    gStyle -> SetOptStat(0); 
   
    TString weight[Nvar];
    TString leg[Nvar];
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
    
    TH1D* h1D_ESE[Nmc];
    for(int imc=0; imc<Nmc; ++imc){
        h1D_ESE[imc] = new TH1D(Form("h1D_ESE_%s",mcName[imc].data()), Form("%s (weighted by %s);The number of iteration;Event Selection Efficiency",var.data(),weight[Nvar-1].Data()), Nvar, 0, Nvar);
        for(int i=0; i<Nvar; ++i){
            TString fname =""; 
            if(i==0) fname = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s.root",noTrigEffCap.Data(),mcName[imc].data(),var.data());
            else fname = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s.root",noTrigEffCap.Data(),mcName[imc].data(),var.data(),weight[i].Data());
           // if(i==0) fname = Form("output/weight_v6_17June30/evtSelEff_foldingMethod%s_run286471_%s_%s.root",noTrigEffCap.Data(),mcName[imc].data(),var.data());
           // else fname = Form("output/weight_v6_17June30/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s.root",noTrigEffCap.Data(),mcName[imc].data(),var.data(),weight[i].Data());
            double effVal = getEffVal(fname);
            h1D_ESE[imc]->SetBinContent(i+1,effVal); 
        } 
    }
    
    float temp_eff[Nmc]={0.0};
    float temp_eff_withTrig[Nmc]={0.0};

    // DRAW
    SetyjPadStyle(); 
    TCanvas* c1 = new TCanvas("c1","",500,500); 
    for(int imc=0; imc<Nmc; ++imc){
        SetHistTextSize(h1D_ESE[imc]);
        h1D_ESE[imc]->SetMarkerStyle(20);
        h1D_ESE[imc]->SetMarkerColor(colHere[imc]);
        h1D_ESE[imc]->SetLineColor(colHere[imc]);
        h1D_ESE[imc]->GetYaxis()->SetRangeUser(0.90,1.03);
        h1D_ESE[imc]->SetNdivisions(510);
        if(imc==0) h1D_ESE[imc]->DrawCopy("pl hist");
        else h1D_ESE[imc]->DrawCopy("pl hist same");
       
        temp_eff[imc] = h1D_ESE[imc]->GetBinContent(nth); 
    }
    jumSun(0,1,Nvar,1);
    
    TLegend* l1 = new TLegend(0.2,0.78,0.55,0.95);
    legStyle(l1);
    for(int imc=0; imc<Nmc; ++imc)
        l1->AddEntry(h1D_ESE[imc],Form("%s",mcName[imc].data()));
    l1->Draw("same");
    for(int imc=1; imc<Nvar; ++imc)
        onSun(imc,0.9,imc,1);    

    // DRAW applying Trigger Efficiency
    TH1D* h1D_ESE_trigEff[Nmc];
    float trigEff[] = {0.97,0.97,0.97};//by data
    //float trigEff[] = {0.987,0.971,0.991};//by each MC 
    TCanvas* c2 = new TCanvas("c2","",500,500); 
    for(int imc=0; imc<Nmc; ++imc){
        h1D_ESE_trigEff[imc]=(TH1D*) h1D_ESE[imc]->Clone(Form("%s_trigEff",h1D_ESE[imc]->GetName()));
        SetHistTextSize(h1D_ESE_trigEff[imc]);
        h1D_ESE_trigEff[imc]->Scale(trigEff[imc]);
        h1D_ESE_trigEff[imc]->SetMarkerStyle(26);
        h1D_ESE_trigEff[imc]->SetMarkerColor(colHere[imc]);
        h1D_ESE_trigEff[imc]->SetLineColor(colHere[imc]);
        h1D_ESE_trigEff[imc]->SetLineStyle(2);
        h1D_ESE_trigEff[imc]->GetYaxis()->SetRangeUser(0.90,1.03);
        h1D_ESE_trigEff[imc]->SetNdivisions(510);
        if(doNoTrigEff){
            if(imc==0) h1D_ESE_trigEff[imc]->DrawCopy("pl hist");
            else h1D_ESE_trigEff[imc]->DrawCopy("pl hist same");
        }

        temp_eff_withTrig[imc] = h1D_ESE_trigEff[imc]->GetBinContent(nth); 
    }
    jumSun(0,1,Nvar,1);
    
    TLegend* l2 = new TLegend(0.5,0.78,0.95,0.95);
    legStyle(l2);
    for(int imc=0; imc<Nmc; ++imc)
        l2->AddEntry(h1D_ESE_trigEff[imc],Form("%s (*trigEff)",mcName[imc].data()));
    if(doNoTrigEff) l2->Draw("same");

    for(int imc=1; imc<Nvar; ++imc)
        onSun(imc,0.9,imc,1);    
    
    c1->cd(); 
    if(doNoTrigEff){
        for(int imc=0; imc<Nmc; ++imc)
            h1D_ESE_trigEff[imc]->DrawCopy("pl hist same");
        l2->Draw("same");
    }
    c1->SaveAs(Form("figures/summary_ESE%s_%s_iteration_with_%s.pdf",noTrigEffCap.Data(),var.data(),weight[Nvar-1].Data()));
    //c2->SaveAs(Form("figures/summary_TrigEff_and_ESE%s_%s_iteration_with_%s.pdf",noTrigEffCap.Data(),var.data(),weight[Nvar-1].Data()));
    EFFSET a;
    for(int imc=0; imc<Nmc; ++imc){
        a.eff[imc] = temp_eff[imc];
        a.eff_withTrig[imc] = temp_eff_withTrig[imc];
    }
    delete c1;
    delete c2;
    delete l1;
    delete l2;
    for(int imc=0; imc<Nmc; ++imc){
        delete h1D_ESE[imc]; 
        delete h1D_ESE_trigEff[imc]; 
    }
    return a;
}
double getEffVal(TString fname){
    TFile* fin = new TFile(fname);
    TH1D* h1D_evtSel = (TH1D*) fin->Get("h1D_data_evtSel");
    TH1D* h1D_corrected = (TH1D*) fin->Get("h1D_data_corrected");
    double Eff = h1D_evtSel->Integral()/h1D_corrected->Integral();
    delete fin;
    return Eff;
}
