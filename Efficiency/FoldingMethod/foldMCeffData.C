// This macro folds the MC efficiency into data by weighting data events with MC efficiency
//
// Javier Martin Blanco
// 25 Jan 2017
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
#include "stdio.h"
#include "string.h"
#include <iostream>
#include "TEfficiency.h"
#include "../../analysisUtils.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"
#include "TParameter.h"
#include "TObjArray.h"

//////////// modify here as you want!! ////////////////
double getXmax(string var)
{
  if(var == "hiET") return 30.;
  if(var == "hiEB") return 100;
  if(var == "hiEE") return 50;
  if(var == "hiHF") return 50;
  if(var == "hiHFhit") return 1000;
  if(var == "hiNpix") return 500;
  if(var == "hiNtracks") return 30;
  
  return 50;
}

TH1* createHisto(const char* hName, string var)
{
  int nBins = 500;
  
  double xMin = 0;
  double xMax = 8000;
  
  if(var == "hiET"){
    xMax = 2000;
  }
  
  if(var == "hiEB"){
    xMax = 5000;
  }
  
  if(var == "hiEE"){
    xMax = 4000;
  }
  
  if(var == "hiHF"){
    xMax = 6000;
  }
  
  if(var == "hiHFhit"){
    xMax = 200000;
  }
  
  if(var == "hiNpix"){
    xMax = 50000;
  }
  
  if(var == "hiNtracks"){
    xMax = 4000;
    nBins = (int)(xMax/2);
  }

  TH1D* h = new TH1D(hName,Form("%s;%s;entries",hName,var.data()),nBins,xMin,xMax);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  
  return h;
}

void foldMCeffData(string mc = "Hydjet", int run = 262620, TString trigcap = "HLT_HIL1MinimumBiasHF1AND_v1", TString evSelCutCap = "pcollisionEventSelection")
{
  TH1::SetDefaultSumw2();
  gStyle -> SetOptStat(0);
  
  /// LOAD TREES AND CREATE HISTOS
  TChain *tref;
  if ( run > 0 )
  {
    tref = new TChain("hiEvtAnalyzer/HiTree","");
    TChain * tskimanalysis = new TChain("skimanalysis/HltTree","");
    TChain * thltanalysis = new TChain("hltanalysis/HltTree","");
    
    const char* fileData = Form("root://cms-xrd-global.cern.ch//eos/cms/store/group/phys_heavyions/jmartinb/HiForestAODPbPb5TeV/HIMinimumBias2/EventTree_PbPb_data_HIMinimumBias2_run%d_15Feb.root",run);
    tref->Add(fileData);
    tskimanalysis->Add(fileData);
    thltanalysis->Add(fileData);
    
    tref->AddFriend(tskimanalysis);
    tref->AddFriend(thltanalysis);
  }
  else
  {
    char line[1024];
    ifstream in("files_HIMinimumBias2.txt");
    tref = new TChain("hiEvtAnalyzer/HiTree","");
    TChain * tskimanalysis = new TChain("skimanalysis/HltTree","");
    TChain * thltanalysis = new TChain("hltanalysis/HltTree","");
    while (in.getline(line,1024,'\n'))
    {
      tref->Add(line);
      tskimanalysis->Add(line);
      thltanalysis->Add(line);
    }
    
    tref->AddFriend(tskimanalysis);
    tref->AddFriend(thltanalysis);
  }

  Long64_t nEntries = tref->GetEntries();
  nEntries = 500000;
  Float_t hiHF,hiHFhit,hiEB,hiEE,hiET;
  int hiNtracks,hiNpix;
  Int_t HLT_HIL1MinimumBiasHF1AND_v1,pcollisionEventSelection;
  tref->SetBranchAddress("hiHF",&hiHF);
  tref->SetBranchAddress("hiHFhit",&hiHFhit);
  tref->SetBranchAddress("hiEB",&hiEB);
  tref->SetBranchAddress("hiEE",&hiEE);
  tref->SetBranchAddress("hiET",&hiET);
  tref->SetBranchAddress("hiNtracks",&hiNtracks);
  tref->SetBranchAddress("hiNpix",&hiNpix);
  tref->SetBranchAddress("hiNpix",&hiNpix);
  tref->SetBranchAddress("HLT_HIL1MinimumBiasHF1AND_v1",&HLT_HIL1MinimumBiasHF1AND_v1);
  tref->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);

  static const int Nvar = 7;
  string vars[Nvar] = {"hiHF","hiHFhit","hiNtracks","hiNpix","hiEB","hiEE","hiET"};
  
  TObjArray* aHistos = new TObjArray();
  aHistos->SetOwner(1);
  TObjArray* aHistosW = new TObjArray();
  aHistosW->SetOwner(1);
  TObjArray* aEff = new TObjArray();
  TH1* h(0x0);
  TH1* hW(0x0);
  TFile* fEff(0x0);
  TEfficiency* eff(0x0);
  for (int i = 0 ; i < Nvar ; i++)
  {
    const char* hName = Form("h_%s",vars[i].c_str());
    h = createHisto(hName,vars[i].c_str());
    hW = static_cast<TH1*>(h->Clone(Form("%s_weighted",hName)));
    
    aHistos->Add(h);
    aHistosW->Add(hW);
    
    fEff = TFile::Open(Form("MCefficiency/eff_%s_%s_%s_%s.root",vars[i].c_str(),trigcap.Data(),evSelCutCap.Data(),mc.data()));
    if (!fEff)
    {
      cout << "[ERROR] No efficiency file found for " << vars[i].c_str() << endl;
      return;
    }
    eff = static_cast<TEfficiency*>(fEff->FindObjectAny(Form("h_eff_%s_%s_%s_%s",vars[i].c_str(),trigcap.Data(),evSelCutCap.Data(),mc.data())));
    if (!eff)
    {
      cout << "[ERROR] No Tefficiency found for " << vars[i].c_str() << endl;
      return;
    }
    aEff->Add(eff);
  }
  
  /// LOOP IN EVENTS TO FILL HISTOS
  cout << "Event loop:" << endl;
  for (Long64_t jentry=0; jentry<nEntries;jentry++)
  {
    tref->GetEntry(jentry);
    
    if ( jentry % 10000 == 0 ) std::cout << jentry << "/" << nEntries << std::endl;
    
    if (!HLT_HIL1MinimumBiasHF1AND_v1 || !pcollisionEventSelection) continue;
      
    for(int i = 0; i < Nvar; ++i)
    {
      double currVar(0.);
      if(vars[i] == "hiET") currVar = hiET;
      if(vars[i] == "hiEB")  currVar = hiEB;
      if(vars[i] == "hiEE")  currVar = hiEE;
      if(vars[i] == "hiHF")  currVar = hiHF;
      if(vars[i] == "hiHFhit")  currVar = hiHFhit;
      if(vars[i] == "hiNpix")  currVar = hiNpix;
      if(vars[i] == "hiNtracks")  currVar = hiNtracks;
      
      const char* hName = Form("h_%s",vars[i].c_str());
      h = static_cast<TH1*>(aHistos->FindObject(hName));
      hW = static_cast<TH1*>(aHistosW->FindObject(Form("%s_weighted",hName)));
      eff = static_cast<TEfficiency*>(aEff->FindObject(Form("h_eff_%s_%s_%s_%s",vars[i].c_str(),trigcap.Data(),evSelCutCap.Data(),mc.data())));
      
      h->Fill(currVar);
      
      int bin = eff->FindFixBin(currVar);
      double effVal = eff->GetEfficiency(bin);
      hW->Fill(currVar,((effVal > 0) ? (1./effVal) : 1.));
    }
    
  }

  /// COMPUTE EFFICIENCIES AND MAKE PLOTS
  
  for(int i = 0; i < Nvar; ++i)
  {
    const char* hName = Form("h_%s",vars[i].c_str());
    h = static_cast<TH1*>(aHistos->FindObject(hName));
    hW = static_cast<TH1*>(aHistosW->FindObject(Form("%s_weighted",hName)));
    
    double effVal = (double)(h->Integral()/hW->Integral());
    cout << "The efficiency from " << vars[i] << " is " << effVal << endl;
    
    TCanvas *c_eff = new TCanvas(Form("c_eff_%s_%s_%s_%s",vars[i].c_str(),trigcap.Data(),evSelCutCap.Data(),mc.data()), "c_eff", 500,500);
    c_eff->SetLogy();
    
    SetHistColor(h,2);
    SetHistColor(hW,1);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
    hW->SetMarkerStyle(21);
    hW->SetMarkerSize(0.8);
    
    TLegend* l1 = new TLegend(0.5,0.6,0.90,0.8);
    legStyle(l1);
    l1->AddEntry(h, "Data", "p");
    l1->AddEntry(hW, "Data corrected", "p");
    
    hW->GetXaxis()->SetRangeUser(0.,getXmax(vars[i]));
    hW->Draw();
    h->Draw("same");
    l1->Draw("same");
    
    double offSet(0);
    if(vars[i] == "hiEB") offSet = 0.2;
    drawText(Form("eff. = %0.3f",effVal),0.2,0.2+0.5);
    drawText(Form("%s eff. weights",mc.data()),0.2+offSet,0.2+0.1);
    drawText(Form("%s",trigcap.Data()),0.2+offSet,0.2+0.05);
    drawText(Form("%s",evSelCutCap.Data()),0.2+offSet,0.2);
    
    string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
    string saveDIR = mainDIR + "/DATAefficiency/";
    void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
    if (dirp) gSystem->FreeDirectory(dirp);
    else gSystem->mkdir(saveDIR.c_str(), kTRUE);
    c_eff->SaveAs(Form("DATAefficiency/eff_%s_%s_%s_%s.pdf",vars[i].c_str(),trigcap.Data(),evSelCutCap.Data(),mc.data()));
    
    TFile* fsave = new TFile(Form("DATAefficiency/eff_%s_%s_%s_%s.root",vars[i].c_str(),trigcap.Data(),evSelCutCap.Data(),mc.data()),"RECREATE");
    c_eff->Write();
    fsave->Close();
  }

}
