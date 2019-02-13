// Author : Yeonju Go

//#include "../yjUtility.h"
#include "/afs/cern.ch/work/y/ygo/private/PbPb2018/CentralityRun2/Utilities/utils.h"
#include <TObject.h>
#include <TH3D.h>

void efficiencyWithZDC(TString sample = "MB", bool doPV=true, bool doHF=true, bool doNOTPV=false, bool doNOTHF=false){
    // WARNING !! doPV/HF and doNOTPV/HF can NOT be 'true' at the same time, and vice versa.
    // WARNING !! If both(e.g. doPV and doNOTPV) are false you just take everything! 
    cout << " :::::: efficiencyWithZDC.C :::::: " << endl;
    
    TString cap = sample;
    if(doPV) cap += "_PV"; 
    if(doHF) cap += "_HF"; 
    if(doNOTPV) cap += "_notPV"; 
    if(doNOTHF) cap += "_notHF"; 

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    //SetyjPadStyle();
    SetHistTitleStyle(0.05,0.03);
    //gStyle->SetLabelOffset(1.2);

    ///////////////////////////////////////////////////////////////////////////
    // Import input file and trees 

    int nFiles = 1000;
    TString fname[nFiles];
    TString capRun = "";

    if(sample=="MB"){ 
        nFiles = 2;
        fname[0] = "/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/PromptRecoForests/v1/000/326/501/HiForest_326501_part1.root";
        fname[1] = "/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/PromptRecoForests/v1/000/326/501/HiForest_326501_part2.root";
        capRun = "Run326501_PromptReco";
    } else if(sample=="ZB"){
        nFiles = 1;
        fname[0]="/eos/cms/store/group/phys_heavyions/dileptons/hanseul/data/HIForward/326942/HIForward_AOD_run326942.root";
        capRun = "Run326942_PromptReco";
        //nFiles = 1;
        //fname[0]="/eos/cms/store/group/phys_heavyions/dileptons/hanseul/data/HIForward/791/AOD/326791/181120_173444/0000/HIForward_PromptReco_run326791.root";
        //capRun = "Run326791_PromptReco";
        //nFiles = 48;
        //for(int i=0;i<nFiles;++i)
        //    fname[i] = Form("/eos/cms/store/group/phys_heavyions/dileptons/hanseul/data/HIForward/790/AOD/326790/181120_173216/0000/HiForestAOD_%d.root",i+1);
        //capRun = "Run326790_PromptReco";
        //capRun = "Run326784-326812_PromptReco";
        //nFiles = 1;
        //fname[0] = "/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIPhysicsForward/v2_reduced/000/326/580/HiForest_326580.root";
        //capRun = "Run326580_Streamer";
    } else if(sample=="EB"){
        nFiles = 32;
        for(int i=0;i<nFiles;++i){
            fname[i] = Form("/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIEmptyBX_RAW/v3/000/326/501/HiForest_%d_numEvent1000.root",i);
        }
        capRun = "Run326501_PromptRaw";
    } else{ 
        cout << " [ERROR] No this kind of sample! " << endl;
    }
    TString treeName = "rechitanalyzerpp/zdcrechit";
    TChain* t = new TChain(Form("%s",treeName.Data()));
    TChain* t_evt = new TChain("hiEvtAnalyzer/HiTree");
    TChain* t_sk = new TChain("skimanalysis/HltTree");
    TChain* t_hlt = new TChain("hltanalysis/HltTree");
    for(int i=0;i<nFiles;++i){
        t->Add(fname[i]);
        t_evt->Add(fname[i]);
        t_sk->Add(fname[i]);
        t_hlt->Add(fname[i]);
    }
    //t->AddFriend(t_evt);
    //t->AddFriend(t_sk);
    //t->AddFriend(t_hlt);

    t->SetBranchStatus("*",0);     // disable all branches
    t->SetBranchStatus("n",1);     // 
    t->SetBranchStatus("e",1);     // 
    t->SetBranchStatus("zside",1);     // 
    t_evt->SetBranchStatus("*",0);     // 
    t_evt->SetBranchStatus("hi*",1);     // 
    t_hlt->SetBranchStatus("*",0);     // 
    t_hlt->SetBranchStatus("*Zero*",1);     // 
    t_hlt->SetBranchStatus("*Unpaired*",1);     // 
    t_hlt->SetBranchStatus("*HIL1NotBptxOR*",1);     // 

    ///////////////////////////////////////////////////////////////////////////
    // SetBranchAddress
    Int_t MAXHITSIZE = 18;
    Int_t n;
    Float_t e[MAXHITSIZE];
    Int_t zside[MAXHITSIZE];
    t->SetBranchAddress("n", &n);
    t->SetBranchAddress("e", &e);
    t->SetBranchAddress("zside", &zside);

    float hiHF, hiHFplus, hiHFminus, hiHFplusEta4, hiHFminusEta4;
    Int_t hiNpix=0;
    Int_t hiNtracks=0;
    Int_t hiNpixPlus=0;
    Int_t hiNpixMinus=0;
    Int_t hiNpixelTracks=0;
    Int_t hiNpixelTracksPlus=0;
    Int_t hiNpixelTracksMinus=0;
    t_evt->SetBranchAddress("hiHF",&hiHF);
    t_evt->SetBranchAddress("hiHFplus",&hiHFplus);
    t_evt->SetBranchAddress("hiHFminus",&hiHFminus);
    t_evt->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
    t_evt->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
    t_evt->SetBranchAddress("hiNtracks",&hiNtracks);
    t_evt->SetBranchAddress("hiNpix",&hiNpix);
    t_evt->SetBranchAddress("hiNpixPlus",&hiNpixPlus);
    t_evt->SetBranchAddress("hiNpixMinus",&hiNpixMinus);
    t_evt->SetBranchAddress("hiNpixelTracks",&hiNpixelTracks);
    t_evt->SetBranchAddress("hiNpixelTracksPlus",&hiNpixelTracksPlus);
    t_evt->SetBranchAddress("hiNpixelTracksMinus",&hiNpixelTracksMinus);

    Int_t phfCoincFilter2Th4 = -1;
    Int_t pprimaryVertexFilter = -1;
    Int_t pclusterCompatibilityFilter = -1;
    t_sk->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4);
    t_sk->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
    t_sk->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
    
    Int_t HLT_HIZeroBias_v1;
    Int_t HLT_HIL1NotBptxOR_v1, HLT_HIL1UnpairedBunchBptxMinus_v1, HLT_HIL1UnpairedBunchBptxPlus_v1;
    t_hlt->SetBranchAddress("HLT_HIZeroBias_v1",&HLT_HIZeroBias_v1);
    t_hlt->SetBranchAddress("HLT_HIL1UnpairedBunchBptxMinus_v1",&HLT_HIL1UnpairedBunchBptxMinus_v1);//when there's one beam only for minus side
    t_hlt->SetBranchAddress("HLT_HIL1UnpairedBunchBptxPlus_v1",&HLT_HIL1UnpairedBunchBptxPlus_v1);//when there's one beam only for plus side
    t_hlt->SetBranchAddress("HLT_HIL1NotBptxOR_v1",&HLT_HIL1NotBptxOR_v1);//when there's no beam
    //t_hlt->SetBranchAddress("",&);


    ///////////////////////////////////////////////////////////////////////////
    // Define histograms 
    int eBin = 100;
    double emax = 200000;
    double emax_z = 200000/4.;
    int nPixelMax = 4000; 
    int nPixelBin = 100; 
    if(sample=="EB"){
        emax = 10000;
        emax_z = 10000/4.;
        nPixelMax = 30;
        nPixelBin = 30;
    }
    if(sample=="ZB"){
        emax = 200000;
        emax_z = 20000;
        nPixelMax = 30;
        nPixelBin = 30;
    }
    if(doNOTPV || doNOTHF){
        nPixelMax = 30;
    }
    TH2D* h2D_e = new TH2D("h2D_e", ";#Sigma E_{ZDC Plus};#Sigma E_{ZDC Minus}", eBin,0,emax, eBin,0,emax);
    TH2D* h2D_e_z = new TH2D("h2D_e_z", ";#Sigma E_{ZDC Plus};#Sigma E_{ZDC Minus}", eBin*2,0,emax_z, eBin*2,0,emax/4.);
    TH2D* h2D_e_pixel = new TH2D("h2D_e_pixel", ";#Sigma E_{ZDC Plus};N_{pixel tracks}", eBin*2,0,emax_z, nPixelMax,0,nPixelMax);
    TH3D* h3D_e_p = new TH3D("h3D_e_p", ";#Sigma E_{ZDC Plus};#Sigma E_{ZDC Minus};N_{pixel tracks} plus", eBin,0,emax, eBin,0,emax, nPixelBin, 0, nPixelMax);
    TH3D* h3D_e_m = new TH3D("h3D_e_m", ";#Sigma E_{ZDC Plus};#Sigma E_{ZDC Minus};N_{pixel tracks} minus", eBin,0,emax, eBin,0,emax, nPixelBin, 0, nPixelMax);
    TH3D* h3D_e = new TH3D("h3D_e", ";#Sigma E_{ZDC Plus};#Sigma E_{ZDC Minus};N_{pixel tracks}", eBin,0,emax, eBin,0,emax, nPixelBin, 0, nPixelMax*2);
    float titleOffset = 1.5;
    h3D_e_p->GetXaxis()->SetTitleOffset(titleOffset);
    h3D_e_p->GetYaxis()->SetTitleOffset(titleOffset);
    h3D_e_p->GetZaxis()->SetTitleOffset(titleOffset);
    h3D_e_m->GetXaxis()->SetTitleOffset(titleOffset);
    h3D_e_m->GetYaxis()->SetTitleOffset(titleOffset);
    h3D_e_m->GetZaxis()->SetTitleOffset(titleOffset);
    h3D_e->GetXaxis()->SetTitleOffset(titleOffset);
    h3D_e->GetYaxis()->SetTitleOffset(titleOffset);
    h3D_e->GetZaxis()->SetTitleOffset(titleOffset);


    ///////////////////////////////////////////////////////////////////////////
    // EVENT LOOP
    int nEvt= 0;
    int nTotEvt= 0;
    Long64_t nentries = t->GetEntries();
    //for(int ientries = 0; ientries < 20000; ++ientries){
    for(int ientries = 0; ientries < nentries; ++ientries){
        t->GetEntry(ientries);
        t_hlt->GetEntry(ientries);
        t_evt->GetEntry(ientries);
        t_sk->GetEntry(ientries);
        if (ientries % 50000 == 0)
            printf("current entry = %d out of %lld : %.3f %%\n", ientries, nentries, (double)ientries / nentries * 100);
       
        if (sample=="ZB" && HLT_HIZeroBias_v1==0) continue; //only Zerobias events
        if (sample=="EB" && !(HLT_HIL1UnpairedBunchBptxMinus_v1 || HLT_HIL1UnpairedBunchBptxPlus_v1|| HLT_HIL1NotBptxOR_v1)) continue; //only Empty Bunch events
        nTotEvt++;

        if (doPV && pprimaryVertexFilter!=1) continue; 
        if (doHF && phfCoincFilter2Th4!=1) continue; 
        if (doNOTPV && pprimaryVertexFilter!=0) continue; 
        if (doNOTHF && phfCoincFilter2Th4!=0) continue; 
        nEvt++;   

        float zdcplus = 0; 
        float zdcminus = 0; 
        float zdctot= 0; 
        for (int iz = 0; iz<n; ++iz){
            if(zside[iz]==-1){ // ZDC minus
               zdcminus+= e[iz]; 
            }else{ // ZDC plus 
               zdcplus += e[iz]; 
            }
        }
        zdctot = zdcplus+zdcminus;
        if(ientries<5)
        cout << hiHF << ", " << hiNpix << ", "<< hiNpixPlus << ", " << hiNpixelTracksPlus << ", " << zdcplus << ", "  << zdcminus << ", "  << zdctot << endl; 

        h2D_e->Fill(zdcplus,zdcminus);
        h2D_e_z->Fill(zdcplus,zdcminus);
        h2D_e_pixel->Fill(zdcplus,hiNpixelTracks);
        h3D_e->Fill(zdcplus,zdcminus,hiNpixelTracks);
        h3D_e_p->Fill(zdcplus,zdcminus,hiNpixelTracksPlus);
        h3D_e_m->Fill(zdcplus,zdcminus,hiNpixelTracksMinus);
    }//event loop
    
    TCanvas* c = new TCanvas("c","",600,600); 
    canvasStyle(c,0.17,0.12,0.12);
    h2D_e->DrawCopy("colz");
    gPad->SetLogz();
    float xpos = 0.2;
    float ypos = 0.85;
    float dy = 0.05;
    drawText(Form("%s",capRun.Data()),xpos,ypos);
    drawText(Form("%s",cap.Data()),xpos,ypos-dy);
    drawText(Form("%d out of %d total events",nEvt,nTotEvt),xpos,ypos-2*dy);
    c->SaveAs(Form("figures/%s_h2D_ZDCplus_ZDCminus_%s.png",cap.Data(),capRun.Data()));
    
    TCanvas* c_z = new TCanvas("c_z","",600,600); 
    canvasStyle(c_z,0.17,0.12,0.12);
    h2D_e_z->DrawCopy("colz");
    gPad->SetLogz();
    drawText(Form("%s",capRun.Data()),xpos,ypos);
    drawText(Form("%s",cap.Data()),xpos,ypos-dy);
    drawText(Form("%d out of %d total events",nEvt,nTotEvt),xpos,ypos-2*dy);
    c_z->SaveAs(Form("figures/%s_h2D_ZDCplus_ZDCminus_zoomIn_%s.png",cap.Data(),capRun.Data()));
    
    TCanvas* c_e_pixel = new TCanvas("c_e_pixel","",600,600); 
    canvasStyle(c_e_pixel,0.17,0.12,0.12);
    h2D_e_pixel->DrawCopy("colz");
    gPad->SetLogz();
    drawText(Form("%s",capRun.Data()),xpos,ypos);
    drawText(Form("%s",cap.Data()),xpos,ypos-dy);
    drawText(Form("%d out of %d total events",nEvt,nTotEvt),xpos,ypos-2*dy);
    c_e_pixel->SaveAs(Form("figures/%s_h2D_ZDCplus_NpixelTracks_%s.png",cap.Data(),capRun.Data()));
    
    TCanvas* c_3d_p = new TCanvas("c_3d_p","",600,600); 
    canvasStyle(c_3d_p);
    //c_3d_p->SetTheta(14.13841);
    //c_3d_p->SetPhi(1.767446);
    c_3d_p->SetTheta(29.24675);
    c_3d_p->SetPhi(14.89516);
    //gStyle->SetCanvasPreferGL(1);
    //gStyle->SetPalette(1);
    h3D_e_p->Draw("lego2");
    drawText(Form("%s",capRun.Data()),xpos,ypos);
    drawText(Form("%s",cap.Data()),xpos,ypos-dy);
    drawText(Form("%d out of %d total events",nEvt,nTotEvt),xpos,ypos-2*dy);
    c_3d_p->SaveAs(Form("figures/%s_h3D_ZDCplus_ZDCminus_NpixelTracksPlus_%s.png",cap.Data(),capRun.Data()));
    
    TCanvas* c_3d_m = new TCanvas("c_3d_m","",600,600); 
    canvasStyle(c_3d_m);
    c_3d_m->SetTheta(10.26551);
    c_3d_m->SetPhi(5.078483);
    h3D_e_m->Draw("lego2");
    drawText(Form("%s",capRun.Data()),xpos,ypos);
    drawText(Form("%s",cap.Data()),xpos,ypos-dy);
    drawText(Form("%d out of %d total events",nEvt,nTotEvt),xpos,ypos-2*dy);
    c_3d_m->SaveAs(Form("figures/%s_h3D_ZDCplus_ZDCminus_NpixelTracksMinus_%s.png",cap.Data(),capRun.Data()));
    
    TCanvas* c_3d= new TCanvas("c_3d","",600,600); 
    canvasStyle(c_3d);
    c_3d->SetTheta(29.24675);
    c_3d->SetPhi(14.89516);
    h3D_e->Draw("lego2");
    drawText(Form("%s",capRun.Data()),xpos,ypos);
    drawText(Form("%s",cap.Data()),xpos,ypos-dy);
    drawText(Form("%d out of %d total events",nEvt,nTotEvt),xpos,ypos-2*dy);
    c_3d->SaveAs(Form("figures/%s_h3D_ZDCplus_ZDCminus_NpixelTracks_%s.png",cap.Data(),capRun.Data()));
    
    //gPad->SetLogz();
    
    TFile* fout = new TFile(Form("output/EfficiencyUsingZDC_%s_%s.root",cap.Data(),capRun.Data()),"RECREATE");
    fout->cd();
    h2D_e_pixel->Write();
    h2D_e_z->Write();
    h2D_e->Write();
    h3D_e->Write();
    h3D_e_m->Write();
    h3D_e_p->Write();
    fout->Close();

}

