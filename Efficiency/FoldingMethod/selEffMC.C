// to draw the distribution of centrality variables
// with the different event selection filters in one canvas.
// should modify evtfilter definition part
//
// Javier Martin Blanco
// 24 Jan 2017
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

void Get1DEffPlots(TTree* t=0, TString var="hiHF", TCut cut="", TString trigcap="", TString evSelCutCap="", string mc="");

void selEffMC(string mc = "Hydjet", TString trigcap = "HLT_HIL1MinimumBiasHF1AND_v1", TString evSelCutCap = "pcollisionEventSelection")
{
  TCut cut = Form("%s && %s",trigcap.Data(),evSelCutCap.Data());
  
  gStyle -> SetOptStat(0);
  
  TFile* infMC;
  if (mc=="Hydjet") infMC = TFile::Open("root://cmsxrootd.fnal.gov///store/user/rbi/merged/Hydjet_Quenched_MinBias_5020GeV_750-HINPbPbWinter16DR-NoPU_75X_mcRun2_HeavyIon_forest_v2/0.root");
  else if (mc=="CymbalEv8_new") infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/mverweij/MC/PbPb/Hydjet_Quenched_Cymbal5Ev8_PbPbMinBias_5020GeV/crab_HiForestCymbal5Ev8offGENSIM/170113_161526/merge/HiForest_0.root");
  else if (mc=="CymbalEv8") infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/mverweij/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev8/crab_HiForestCymbal5Ev8/160912_154217/merge/HiForest_0.root");
  else if (mc=="CymbalEv7") infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/mverweij/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev7/crab_HiForestCymbal5Ev7/160912_154053/merge/HiForest_0.root");
  else if (mc=="CymbalEv6") infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/mverweij/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev6/crab_HiForestCymbal5Ev6/160912_153731/merge/HiForest_0.root");
  else if (mc=="CymbalEv5") infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/cmst3/user/mverweij/jetsPbPb/Run2/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev5/crab_HiForestCymbal5Ev5/merge/HiForest_0.root");
  else if (mc=="CymbalEv4") infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/cmst3/user/mverweij/jetsPbPb/Run2/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev4/crab_HiForestCymbal5Ev4/merge/HiForest_0.root");
  else if (mc=="CymbalEv3") infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/mverweij/MC/PbPb/FixedHydjet_Quenched_MinBias_Cymbal5Ev3/crab_HiForestCymbal5Ev3Fixed/160912_093619/merge/HiForest_0.root");
  else if (mc=="CymbalEv2") infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/cmst3/group/hintt/mverweij/CS/MC/PbPb/Hydjet_Quenched_MinBias_Cymbal5Ev2/merge/HiForest_0.root");
  else if (mc=="CymbalEv1") infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/cmst3/group/hintt/mverweij/CS/MC/PbPb/Hydjet_Cymbal5Ev1_MinBias_5020GeV/merge/HiForest_0.root");
  
  TTree *t = (TTree*)infMC->Get("hiEvtAnalyzer/HiTree");
  t->AddFriend("hltanalysis/HltTree");
  t->AddFriend("skimanalysis/HltTree");
  
  double Nevt_t = t -> GetEntries(cut);
  cout << "cut is '" << cut.GetTitle() << "' and # of events = " << Nevt_t << endl;

  Get1DEffPlots(t, "hiHF",cut,trigcap,evSelCutCap,mc);
  Get1DEffPlots(t, "hiNtracks",cut,trigcap,evSelCutCap,mc);
  Get1DEffPlots(t, "hiHFhit",cut,trigcap,evSelCutCap,mc);
  Get1DEffPlots(t, "hiNpix",cut,trigcap,evSelCutCap,mc);
  Get1DEffPlots(t, "hiET",cut,trigcap,evSelCutCap,mc);
  Get1DEffPlots(t, "hiEE",cut,trigcap,evSelCutCap,mc);
  Get1DEffPlots(t, "hiEB",cut,trigcap,evSelCutCap,mc);
  
}

void Get1DEffPlots(TTree* t, TString var, TCut cut, TString trigcap, TString evSelCutCap, string mc)
{
  TString cap(Form("%s_%s_%s",trigcap.Data(),evSelCutCap.Data(),mc.data()));
  
  int nBins = 100;
  
  double xMin = 0;
  double xMax = 8000;
  double rangePlot = 100;
  
  if(var == "hiET"){
    xMax = 2000;
    nBins = (int)(xMax/2);
    rangePlot = 30;
  }
  
  if(var == "hiEB"){
    xMax = 5000;
    nBins = (int)(xMax/2.5);
    rangePlot = 100;
  }
  
  if(var == "hiEE"){
    xMax = 4000;
    nBins = (int)(xMax/2);
    rangePlot = 50;
  }
  
  if(var == "hiHF"){
    xMax = 6000;
    nBins = (int)(xMax/3);
    rangePlot = 50;
  }
  
  if(var == "hiHFhit"){
    xMax = 200000;
    nBins = (int)(xMax/50);
    rangePlot = 1000;
  }
  
  if(var == "hiNpix"){
    xMax = 50000;
    nBins = (int)(xMax/25);
    rangePlot = 500;
  }
  
  if(var == "hiNtracks"){
    xMax = 4000;
    nBins = (int)(xMax);
    rangePlot = 30;
  }

  
  TH1D* htot = new TH1D(Form("htot_%s",var.Data()), Form("efficiency;%s;",var.Data()), nBins,xMin,xMax);
  TH1D* hpass = (TH1D*)htot->Clone(Form("hpass_%sd",var.Data()));
  TH1D* dummyH = (TH1D*)htot->Clone("hDummy");
  
  t->Draw(Form("%s>>%s",var.Data(),htot->GetName()), "");
  t->Draw(Form("%s>>%s",var.Data(),hpass->GetName()), cut);
  
  TParameter<double>* effVal = new TParameter<double>(Form("v_eff_%s_%s",var.Data(),cap.Data()),(double)(hpass->GetEntries()/htot->GetEntries()));
  
  TEfficiency* effH = new TEfficiency(*hpass,*htot);
  effH->SetName(Form("h_eff_%s_%s",var.Data(),cap.Data()));
  
  TCanvas *c_eff = new TCanvas(Form("c_eff_%s_%s",var.Data(),cap.Data()), "c_eff", 500,500);
  TGraphAsymmErrors* geff = effH->CreateGraph();
  geff->SetName(Form("g_eff_%s_%s",var.Data(),cap.Data()));

  dummyH->GetXaxis()->SetRangeUser(0.,rangePlot);
  dummyH->Draw();
  geff->Draw("same");
  
  drawText(Form("eff. = %0.3f",effVal->GetVal()),0.2,0.2+0.5);
  drawText(mc.data(),0.2,0.2+0.1);
  drawText(Form("%s",trigcap.Data()),0.2,0.2+0.05);
  drawText(Form("%s",evSelCutCap.Data()),0.2,0.2);
  
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/MCefficiency/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  c_eff->SaveAs(Form("MCefficiency/eff_%s_%s.pdf",var.Data(),cap.Data()));
  
  TFile* fsave = new TFile(Form("MCefficiency/eff_%s_%s.root",var.Data(),cap.Data()),"RECREATE");
  c_eff->Write();
  geff->Write();
  effH->Write();
  effVal->Write();
  fsave->Close();

}
