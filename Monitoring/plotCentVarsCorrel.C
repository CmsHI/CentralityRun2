// to draw the correlation plots of centrality variables
//
// Author Javier Martin
// 9 Nov 2016
//
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
//basic c++ header, string ...
#include <iostream>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <math.h>

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
#include "TSystem.h"
#include "../utils.h"

//const char* trigcut = "(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1)";
const char* trigcut = "HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1";
//const char* trigcap = "HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part*";
const char* trigcap = "HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1";
const char* evSelCut = "pBeamScrapingFilter && pPAprimaryVertexFilter && phfCoincFilter1 && pVertexFilterCutG";
//const char* evSelCut = "!pBeamScrapingFilter";
//const char* evSelCut = "!pPAprimaryVertexFilter";
//const char* evSelCut = "!phfCoincFilter1";
//const char* evSelCut = "!pVertexFilterCutG";
//const char* evSelCutCap = evSelCut;
const char* evSelCutCap = "BS+PV+HFC+PVG";
const string cap = "pAExpress_5TeV_run285090";
//const char* cap = "MBFullEvtSel";
TCut selCut = Form("%s && %s",trigcut,evSelCut);
double maxZ = 10000;

void Get2DPlots(TTree* t_evt=0, TString v1="hiHF", TString v2="hiNpix", int xbin=200, double xmin=0, double xmax=4500, int ybin=200, double ymin=0, double ymax=10000, TCut cut="");

//void plotCentVarsCorrel(const char* fname="root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pPb_EposMinBias_5TeV_8022_Forest/EPOS5TeV_GEN_SIM/crab_pPb_EposMinBias_5TeV_8022_Forest_corrCentrTable/161108_122954/0000/HiForestAOD_92.root")
void plotCentVarsCorrel(const char* fname="root://eoscms//eos/cms/store/group/phys_heavyions/kjung/ExpressForests/v1/Merged/HiForest_run285090_Express.root")
{
  cout << "Aplying the following event selection: " << selCut.GetTitle() << endl;
  
  TH1::SetDefaultSumw2();
  gStyle -> SetOptStat(0);
  SetHistTitleStyle(0.06,0.04);
  SetyjPadStyle();
  
  TFile *fin = TFile::Open(fname);
  TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
  TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
  TTree *t_hlt = (TTree*) fin -> Get("hltanalysis/HltTree");
  t_evt -> AddFriend(t_hlt);
  t_evt -> AddFriend(t_skim);
  
  cout << t_evt->GetEntries() << endl;
  
  double hiHFMax = 300;
  double hiHFSSMax = 200;
  double hiHFSSTruncMax = 80;
  double hiBinMax = 200;
  double hiHFhitMax = 8000;
  double hiETMax = 150;
  double hiEEMax = 120;
  double hiEBMax = 1200;
  double hiNpixMax = 2000;
  double hiNpixelTracksMax = 500;
  double hiNtracksMax = 400;
  double hiNtracksCutEtaMax = 350;
  double hiNtracksCutMax = 150;
  double hiZDCMax = 20000;

  int nbin = 100;
  
  Get2DPlots(t_evt, "hiHFplus","hiHFminus",nbin,0,hiHFSSMax,nbin,0,hiHFSSMax,selCut);
  Get2DPlots(t_evt, "hiHFhitPlus","hiHFhitMinus",nbin,0,hiHFhitMax,nbin,0,hiHFhitMax,selCut);
  Get2DPlots(t_evt, "hiHF","hiNpix",nbin,0,hiHFMax,nbin,0,hiNpixMax,selCut);
  Get2DPlots(t_evt, "hiHF","hiZDC",nbin,0,hiHFMax,nbin,0,hiZDCMax,selCut);
  Get2DPlots(t_evt, "hiHF","hiNtracks",nbin,0,hiHFMax,nbin,0,hiNtracksMax,selCut);
  Get2DPlots(t_evt, "hiHFplusEta4","hiHFminusEta4",nbin,0,hiHFSSTruncMax,nbin,0,hiHFSSTruncMax,selCut);
  Get2DPlots(t_evt, "hiHFplusEta4","hiNtracks",nbin,0,hiHFSSTruncMax,nbin,0,hiNtracksMax,selCut);
  Get2DPlots(t_evt, "hiHFplusEta4","hiZDC",nbin,0,hiHFSSTruncMax,nbin,0,hiZDCMax,selCut);
  Get2DPlots(t_evt, "hiHFplusEta4","hiZDCplus",nbin,0,hiHFSSTruncMax,nbin,0,hiZDCMax,selCut);
  Get2DPlots(t_evt, "hiHFplusEta4","hiZDCminus",nbin,0,hiHFSSTruncMax,nbin,0,hiZDCMax,selCut);
  Get2DPlots(t_evt, "hiNpix","hiZDC",nbin,0,hiNpixMax,nbin,0,hiZDCMax,selCut);
  Get2DPlots(t_evt, "hiZDC","hiNtracks",nbin,0,hiZDCMax,nbin,0,hiNtracksMax,selCut);
  
}

void Get2DPlots(TTree* t_evt, TString v1, TString v2, int xbin, double xmin, double xmax, int ybin, double ymin, double ymax, TCut cut)
{
  TCanvas *c1 = new TCanvas(Form("%sVS%s",v2.Data(), v1.Data()), Form("%sVS%s",v2.Data(), v1.Data()), 500,500);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetLeftMargin(0.1677852);
  c1->SetRightMargin(0.1342282);
  c1->SetTopMargin(0.05839416);
  c1->SetBottomMargin(0.1715328);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  c1->cd();
  
  TH2D *h2D = new TH2D(Form("h2D_%s_%s_filter_%s",v1.Data(),v2.Data(),cap.data()), Form(";%s;%s", v1.Data(), v2.Data()), xbin, xmin, xmax, ybin, ymin, ymax);
  t_evt->Draw(Form("%s:%s>>+%s",v2.Data(), v1.Data(), h2D->GetName() ), selCut);
  h2D = (TH2D*)gDirectory->Get(h2D->GetName());
  h2D->GetZaxis()->SetRangeUser(0.,maxZ);

  h2D->Draw("colz");
  gPad->SetLogz();
  
  drawText(cap.data(),0.2,0.2+0.7);
  //drawText(cut1.GetTitle(),0.2,0.2);
  drawText(Form("%s",trigcap),0.2,0.2+0.65);
  drawText(Form("%s",evSelCutCap),0.2,0.2+0.6);
  
  // save canvas
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/2DcentralityDistributions/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  c1->SaveAs(Form("2DcentralityDistributions/centVarsCorrel_%s_%s_%s.pdf", v1.Data(),v2.Data(),cap.data()));
}

