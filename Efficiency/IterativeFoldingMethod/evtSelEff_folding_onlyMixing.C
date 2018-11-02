// This macro folds the MC efficiency into data by weighting data events with MC efficiency
// Created : 12 June 2017
// Modified : 16 June 2017
// Author : Yeonju Go // 

#include "EvtSelFoldingIteration.h" 
const int colHere[]={1,kAzure+2,8,4,2,kYellow+2}; 
void drawEvtSelEff(EffHists a, TString var, string mcName, TString cap="", const int num=0, bool doNoTrigEff=1); 
void evtSelEff_folding_onlyMixing(bool doNoTrigEff=0, string mcName = "HIJING", bool do0=1, bool do1=0, bool do2=0, bool do3=0, bool do4=0, bool do5=0) { 
    TH1::SetDefaultSumw2(); 
    gStyle -> SetOptStat(0); 
    TString var[] = {"hiHF","hiNpix","hiNtracks","hiHFplusEta4","hiHFminusEta4"};
    //TString var[] = {"hiHF","hiNpix","hiNtracks","hiHFplus","hiHFminus","hiHFhit"};
    int Nvar=5;
    //EffHists getEvtSelEff(bool doNoTrigEff, string mcName, TString var, int num, TH1D* weightHist[], TString weightVar[])

    //No correction
    if(do0){
        EffHists iter0[Nvar];
        TH1D* weightHist[1];
        weightHist[0] = new TH1D("weightHist","",1,0,40);
        TString weightVar[1];
        weightVar[0] ="";
        for(int iv=0;iv<Nvar;++iv){
            iter0[iv] = getEvtSelEff(doNoTrigEff,mcName,var[iv],0,weightHist,weightVar); 
            if(iter0[iv].MC_integFilterRate==0){
               //cout << "okok" << endl;
               continue;
            }
            drawEvtSelEff(iter0[iv],var[iv],mcName,"",0,doNoTrigEff);
        }
        delete weightHist[0];
    }
    TString noTrigEffCap = ""; 
    if(doNoTrigEff) noTrigEffCap = "_noTrigEff"; 
    //1st iteration weighing factor
    if(do1){
        EffHists iter1[Nvar][Nvar];
        for(int iv_w=0;iv_w<Nvar;++iv_w){ // 1st weighting by this variable
            for(int iv=0;iv<Nvar;++iv){ // 
                if(iv_w==iv) continue;
                
                TString weightVar[1];
                weightVar[0] =var[iv_w];
                TString fname = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[0].Data()); 
                TFile* fin = new TFile(fname);
                TH1D* weightHist[1];
                weightHist[0] = (TH1D*) fin->Get("h1D_weightFactor");
                iter1[iv_w][iv] = getEvtSelEff(doNoTrigEff,mcName,var[iv],1,weightHist,weightVar); 
                if(iter1[iv_w][iv].MC_integFilterRate==0) continue;
                TString cap = Form("1stWeightedBy_%s",var[iv_w].Data());
                drawEvtSelEff(iter1[iv_w][iv],var[iv],mcName,cap,1,doNoTrigEff);
                delete fin;
            }
        }
    }
   
    //2nd iteration weighing factor
    if(do2){
        const int Nweight =2;
        EffHists iter2[Nvar][Nvar][Nvar];
        for(int iv_2w=0;iv_2w<Nvar;++iv_2w){ // 2nd weighting by this variable
        for(int iv_w=0;iv_w<Nvar;++iv_w){ // 1st weighting by this variable
            for(int iv=0;iv<Nvar;++iv){ //drawing variable
                if(iv==iv_2w) continue;
                if(iv_w==iv_2w) continue;
                if(iv==iv_w) continue;

                TString weightVar[Nweight];
                weightVar[0] =var[iv_w];
                weightVar[1] =var[iv_2w];
                // Get weightHist from the files
                TString fname[Nweight];
                fname[0] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[0].Data());
                fname[1] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[1].Data(),weightVar[0].Data());

                TH1D* weightHist[Nweight];
                TFile* fin[Nweight];
                for(int i=0;i<Nweight;++i){
                    fin[i] = new TFile(fname[i]);
                    weightHist[i] = (TH1D*) fin[i]->Get("h1D_weightFactor");
                    weightHist[i] ->SetName(Form("h1D_weightFactor_%d",i));
                }
                iter2[iv_2w][iv_w][iv] = getEvtSelEff(doNoTrigEff,mcName,var[iv],2,weightHist,weightVar); 
                if(iter2[iv_2w][iv_w][iv].MC_integFilterRate==0) continue;
                TString cap = Form("1stWeightedBy_%s_2ndWeightedBy_%s",var[iv_w].Data(),var[iv_2w].Data());
                drawEvtSelEff(iter2[iv_2w][iv_w][iv],var[iv],mcName,cap,Nweight,doNoTrigEff);
                for(int i=0;i<Nweight;++i){
                    delete fin[i]; 
                }
            }
        }
        }
    }

    //3rd iteration weighing factor
    if(do3){
         const int Nweight =3;
         EffHists iter3[Nvar][Nvar][Nvar][Nvar];
         for(int iv_3w=0;iv_3w<Nvar;++iv_3w){ // 3rd weighting by this variable
         for(int iv_2w=0;iv_2w<Nvar;++iv_2w){ // 2nd weighting by this variable
         for(int iv_w=0;iv_w<Nvar;++iv_w){ // 1st weighting by this variable
             for(int iv=0;iv<Nvar;++iv){ //
                 if(iv==iv_3w || iv==iv_2w || iv==iv_w) continue;
                 if(iv_3w==iv_w || iv_2w==iv_w) continue;
                 if(iv_3w==iv_2w) continue;

                 TString weightVar[Nweight];
                 weightVar[0] =var[iv_w];
                 weightVar[1] =var[iv_2w];
                 weightVar[2] =var[iv_3w];
                 // Get weightHist from the files
                 TString fname[Nweight];
                 fname[0] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[0].Data());
                 fname[1] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[1].Data(),weightVar[0].Data());
                 fname[2] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[2].Data(),weightVar[0].Data(),weightVar[1].Data());
                 TH1D* weightHist[Nweight];
                 TFile* fin[Nweight];
                 for(int i=0;i<Nweight;++i){
                     fin[i] = new TFile(fname[i]);
                     weightHist[i] = (TH1D*) fin[i]->Get("h1D_weightFactor"); 
                     weightHist[i]->SetName(Form("h1D_weightFactor_%d",i));
                 }
                 //Get Efficiency Hist
                 iter3[iv_3w][iv_2w][iv_w][iv] = getEvtSelEff(doNoTrigEff,mcName,var[iv],Nweight,weightHist,weightVar); 
                 if(iter3[iv_3w][iv_2w][iv_w][iv].MC_integFilterRate==0) continue;
                 TString cap = Form("1stWeightedBy_%s_2ndWeightedBy_%s_3rdWeightedBy_%s",var[iv_w].Data(),var[iv_2w].Data(),var[iv_3w].Data());
                 drawEvtSelEff(iter3[iv_3w][iv_2w][iv_w][iv],var[iv],mcName,cap,Nweight,doNoTrigEff);
                 for(int i=0;i<Nweight;++i){
                    delete fin[i]; 
                 }
             }            
         }}}
    }

    //4th iteration weighing factor
    if(do4){
         const int Nweight =4;
         EffHists iter4[Nvar][Nvar][Nvar][Nvar][Nvar];
         for(int iv_4w=0;iv_4w<Nvar;++iv_4w){ // 4th weighting by this variable
         for(int iv_3w=0;iv_3w<Nvar;++iv_3w){ // 3rd weighting by this variable
         for(int iv_2w=0;iv_2w<Nvar;++iv_2w){ // 2nd weighting by this variable
         for(int iv_w=0;iv_w<Nvar;++iv_w){ // 1st weighting by this variable
             for(int iv=0;iv<Nvar;++iv){ //
                 if(iv==iv_4w || iv==iv_3w || iv==iv_2w || iv==iv_w) continue;
                 if(iv_w==iv_4w || iv_w==iv_3w || iv_w==iv_2w) continue;
                 if(iv_2w==iv_4w || iv_2w==iv_3w) continue;
                 if(iv_3w==iv_4w) continue;

                 TString weightVar[Nweight];
                 weightVar[0] =var[iv_w];
                 weightVar[1] =var[iv_2w];
                 weightVar[2] =var[iv_3w];
                 weightVar[3] =var[iv_4w];
                 // Get weightHist from the files
                 TString fname[Nweight];
                 fname[0] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[0].Data());
                 fname[1] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[1].Data(),weightVar[0].Data());
                 fname[2] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[2].Data(),weightVar[0].Data(),weightVar[1].Data());
                 fname[3] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s_%s_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[3].Data(),weightVar[0].Data(),weightVar[1].Data(),weightVar[2].Data());
                 TH1D* weightHist[Nweight];
                 TFile* fin[Nweight];
                 for(int i=0;i<Nweight;++i){
                     fin[i] = new TFile(fname[i]);
                     weightHist[i] = (TH1D*) fin[i]->Get("h1D_weightFactor");
                     weightHist[i]->SetName(Form("h1D_weightFactor_%d",i));
                 }
                 //Get Efficiency Hist
                 iter4[iv_4w][iv_3w][iv_2w][iv_w][iv] = getEvtSelEff(doNoTrigEff,mcName,var[iv],Nweight,weightHist,weightVar);
                 if(iter4[iv_4w][iv_3w][iv_2w][iv_w][iv].MC_integFilterRate==0) continue;
                 TString cap = Form("1stWeightedBy_%s_2ndWeightedBy_%s_3rdWeightedBy_%s_4thWeightedBy_%s",var[iv_w].Data(),var[iv_2w].Data(),var[iv_3w].Data(),var[iv_4w].Data());
                 drawEvtSelEff(iter4[iv_4w][iv_3w][iv_2w][iv_w][iv],var[iv],mcName,cap,Nweight,doNoTrigEff);
                 for(int i=0;i<Nweight;++i){
                    delete fin[i];
                 }
             }
         }}}}
    }

    //5th iteration weighing factor
    if(do5){
         const int Nweight =5;
         EffHists iter[Nvar][Nvar][Nvar][Nvar][Nvar][Nvar];
         for(int iv_5w=0;iv_5w<Nvar;++iv_5w){ // 5th weighting by this variable
         for(int iv_4w=0;iv_4w<Nvar;++iv_4w){ // 4th weighting by this variable
         for(int iv_3w=0;iv_3w<Nvar;++iv_3w){ // 3rd weighting by this variable
         for(int iv_2w=0;iv_2w<Nvar;++iv_2w){ // 2nd weighting by this variable
         for(int iv_w=0;iv_w<Nvar;++iv_w){ // 1st weighting by this variable
             for(int iv=0;iv<Nvar;++iv){ //
                 if(iv==iv_5w || iv==iv_4w || iv==iv_3w || iv==iv_2w || iv==iv_w) continue;
                 if(iv_w==iv_5w || iv_w==iv_4w || iv_w==iv_3w || iv_w==iv_2w) continue;
                 if(iv_2w==iv_5w || iv_2w==iv_4w || iv_2w==iv_3w) continue;
                 if(iv_3w==iv_5w || iv_3w==iv_4w) continue;
                 if(iv_4w==iv_5w) continue;
                 
                 TString weightVar[Nweight];
                 weightVar[0] =var[iv_w];
                 weightVar[1] =var[iv_2w];
                 weightVar[2] =var[iv_3w];
                 weightVar[3] =var[iv_4w];
                 weightVar[4] =var[iv_5w];
                 // Get weightHist from the files
                 TString fname[Nweight];
                 fname[0] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[0].Data());
                 fname[1] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[1].Data(),weightVar[0].Data());
                 fname[2] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[2].Data(),weightVar[0].Data(),weightVar[1].Data());
                 fname[3] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s_%s_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[3].Data(),weightVar[0].Data(),weightVar[1].Data(),weightVar[2].Data());
                 fname[4] = Form("output/evtSelEff_foldingMethod%s_run286471_%s_%s_WeightedBy_%s_%s_%s_%s.root",noTrigEffCap.Data(),mcName.data(),weightVar[4].Data(),weightVar[0].Data(),weightVar[1].Data(),weightVar[2].Data(),weightVar[3].Data());
                 TH1D* weightHist[Nweight];
                 TFile* fin[Nweight];
                 for(int i=0;i<Nweight;++i){
                     fin[i] = new TFile(fname[i]);
                     weightHist[i] = (TH1D*) fin[i]->Get("h1D_weightFactor");
                     weightHist[i]->SetName(Form("h1D_weightRactor_%d",i));
                 }
                 //Get Efficiency Hist
                 iter[iv_5w][iv_4w][iv_3w][iv_2w][iv_w][iv] = getEvtSelEff(doNoTrigEff,mcName,var[iv],Nweight,weightHist,weightVar);
                 if(iter[iv_5w][iv_4w][iv_3w][iv_2w][iv_w][iv].MC_integFilterRate==0) continue;
                 TString cap = Form("1stWeightedBy_%s_2ndWeightedBy_%s_3rdWeightedBy_%s_4thWeightedBy_%s_5thWeightedBy_%s",var[iv_w].Data(),var[iv_2w].Data(),var[iv_3w].Data(),var[iv_4w].Data(),var[iv_5w].Data());
                 drawEvtSelEff(iter[iv_5w][iv_4w][iv_3w][iv_2w][iv_w][iv],var[iv],mcName,cap,Nweight,doNoTrigEff);
                 for(int i=0;i<Nweight;++i){
                    delete fin[i];
                 }

             }
         }}}}}
    }
}


void drawEvtSelEff(EffHists a, TString var, string mcName, TString cap, const int num, bool doNoTrigEff){
    SetyjPadStyle();
    
    TString noTrigEffCap = ""; 
    if(doNoTrigEff) noTrigEffCap = "noTrigEff_"; 
    
    //filterRate
    TCanvas* c_filterRate = new TCanvas(Form("c_filterRate_%s_%s",var.Data(),cap.Data()), "",500,500);
    TEfficiency* pEff=0;
    pEff = new TEfficiency(*a.h1D_mc_evtSel_weighted,*a.h1D_mc_tot_weighted);
    pEff->SetTitle(Form(";%s;%s Filter Rate",var.Data(),mcName.data()));
    pEff->SetMarkerStyle(20);
    pEff->Draw("AP");
    drawText(Form("%s Filter Rate = %.2f %s",mcName.data(),a.MC_integFilterRate*100,"%"),0.2,0.4,0,kBlack,26);
    double tempFilterRate = (a.h1D_mc_evtSel_weighted->Integral()/a.h1D_mc_tot_weighted->Integral()); 
    //drawText(Form("%s Filter Rate temp = %.2f %s",mcName.data(),tempFilterRate*100,"%"),0.4,0.5,0,kBlack,17);
    if(num>0) drawText(cap,0.93,0.3,1,kBlack,12);
    double xmax = a.h1D_mc_tot->GetBinLowEdge(a.h1D_mc_tot->GetNbinsX())+a.h1D_mc_tot->GetBinWidth(a.h1D_mc_tot->GetNbinsX()); 
    jumSun(0,1,xmax,1);
    c_filterRate->SaveAs(Form("figures/%sfilterRate_%s_%s_%s.pdf",noTrigEffCap.Data(),mcName.data(),var.Data(),cap.Data()));

    //Correction Factor
    TCanvas* c_corrFactor = new TCanvas(Form("c_corrFactor_%s_%s",var.Data(),cap.Data()), "",500,500);
    SetHistTextSize(a.h1D_weightFactor);
    a.h1D_weightFactor->SetMarkerStyle(20);
    a.h1D_weightFactor->SetMarkerColor(1);
    a.h1D_weightFactor->SetTitle(Form(";%s;Weight Factor = Data_corrected/%s",var.Data(),mcName.data()));
    a.h1D_weightFactor->DrawCopy("p");
    if(num>0) drawText(cap,0.15,0.8,0,kBlack,12); 
    jumSun(0,1,xmax,1);
    c_corrFactor->SaveAs(Form("figures/%scorrectionFactor_%s_%s_%s.pdf",noTrigEffCap.Data(),mcName.data(),var.Data(),cap.Data()));
            
    //Distribution
    TCanvas* c_dist = new TCanvas(Form("c_dist_%s_%s",var.Data(),cap.Data()), "",500,500);
    a.h1D_mc_tot->Scale(1.,"width");
    a.h1D_mc_tot_weighted->Scale(1.,"width");
    a.h1D_data_evtSel->Scale(1.,"width");
    a.h1D_mc_tot->Scale(1./a.h1D_mc_tot->Integral());
    a.h1D_mc_tot_weighted->Scale(1./a.h1D_mc_tot_weighted->Integral());
    double data_evtSel_area = a.h1D_data_evtSel->Integral();
    a.h1D_data_evtSel->Scale(1./data_evtSel_area);
   // a.h1D_mc_tot->Scale(1./a.h1D_mc_tot->Integral(),"width");
   // a.h1D_mc_tot_weighted->Scale(1./a.h1D_mc_tot_weighted->Integral(),"width");
   // a.h1D_data_evtSel->Scale(1./a.h1D_data_evtSel->Integral(),"width");
    
    a.h1D_mc_tot->SetMarkerSize(0.8);
    a.h1D_mc_tot_weighted->SetMarkerSize(0.8);
    a.h1D_data_evtSel->SetMarkerSize(0.8);
    a.h1D_mc_tot->SetMarkerStyle(21);
    a.h1D_mc_tot_weighted->SetMarkerStyle(20);
    a.h1D_data_evtSel->SetMarkerStyle(20);
    a.h1D_mc_tot->SetMarkerColor(kAzure+2);
    a.h1D_mc_tot_weighted->SetMarkerColor(8);
    a.h1D_data_evtSel->SetMarkerColor(1);

    //calculate the maximum range
    double maxY = cleverRange(a.h1D_mc_tot);
    double tempMaxY1 = cleverRange(a.h1D_mc_tot_weighted);
    double tempMaxY2 = cleverRange(a.h1D_data_evtSel);
    if(tempMaxY1 >= maxY) maxY = tempMaxY1;
    if(tempMaxY2 >= maxY) maxY = tempMaxY2;
    a.h1D_data_corrected->GetYaxis()->SetRangeUser(0,maxY);
    a.h1D_mc_tot->GetYaxis()->SetRangeUser(0,maxY);
    a.h1D_mc_tot_weighted->GetYaxis()->SetRangeUser(0,maxY);

    
    SetHistTextSize(a.h1D_mc_tot); 
    a.h1D_mc_tot->SetTitle(Form(";%s;Arbitrary Normalization",var.Data()));
    a.h1D_mc_tot->DrawCopy();
    if(num>0) a.h1D_mc_tot_weighted->DrawCopy("same");
    a.h1D_data_evtSel->DrawCopy("same");
    TLegend* l1 = new TLegend(0.4,0.7,0.9,0.9);
    legStyle(l1);
    l1->AddEntry(a.h1D_data_evtSel,"DATA selected","pl");
    l1->AddEntry(a.h1D_mc_tot,Form("%s (triggered)",mcName.data()),"pl");
    if(num>0) l1->AddEntry(a.h1D_mc_tot_weighted,Form("%s (triggered) weighted",mcName.data()),"pl");
    l1->Draw("same");
    if(num>0) drawText(cap,0.93,0.6,1,kBlack,12); 
    
    c_dist->SaveAs(Form("figures/%sdist_%s_%s_%s.pdf",noTrigEffCap.Data(),mcName.data(),var.Data(),cap.Data()));
    
    //EventSelectionEfficiency
    TCanvas* c_evtSelEff = new TCanvas(Form("c_evtSelEff_%s_%s",var.Data(),cap.Data()), "",500,500);
    c_evtSelEff->SetLogy();
    a.h1D_data_corrected->Scale(1.,"width");
    a.h1D_data_evtSel->Scale(data_evtSel_area);
    SetHistTextSize(a.h1D_data_corrected);
    a.h1D_data_corrected->SetMarkerSize(0.8);
    a.h1D_data_corrected->SetMarkerColor(1);
    a.h1D_data_evtSel->SetMarkerColor(2);
    a.h1D_data_corrected->SetMarkerStyle(21);
    a.h1D_data_evtSel->SetMarkerStyle(20);
    a.h1D_data_corrected->DrawCopy(); 
    a.h1D_data_evtSel->DrawCopy("same");
    TLegend* l2 = new TLegend(0.6,0.7,0.9,0.9);
    legStyle(l2);
    l2->AddEntry(a.h1D_data_evtSel,"Data","pl");
    l2->AddEntry(a.h1D_data_corrected,"Data corrected","pl");
    l2->Draw("same");
    drawText(Form("Evt. Sel. Eff. = %.2f %s",a.DATA_integEvtSelEff*100,"%"),0.2,0.3,0,kBlack,26);
    if(num>0) drawText(cap,0.15,0.2,0,kBlack,12); 
    //double tempEff = a.h1D_data_evtSel->Integral("width")/a.h1D_data_corrected->Integral("width");
    //drawText(Form("Evt. Sel. Eff. (temp) = %.2f %s",tempEff*100,"%"),0.3,0.4,0,kBlack,17);
    c_evtSelEff->SaveAs(Form("figures/%sevtSelEff_%s_%s_%s.pdf",noTrigEffCap.Data(),mcName.data(),var.Data(),cap.Data()));

    delete c_filterRate;    
    delete c_corrFactor;    
    delete c_dist;    
    delete c_evtSelEff;    
    delete l1;    
    delete l2;    
    delete pEff;    

}
