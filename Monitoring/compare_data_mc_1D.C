// comparison centrality objects between DATA and MC
// Author : Yeonju Go
// 24 Nov 2015

//basic c++ header, string ...
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <math.h>
//tree, hist, vector ...
#include <TROOT.h>
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TMath.h>
#include <math.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>
#include <TCut.h>
#include <TString.h>
//canvas, legend, latex ... //cosmetic
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TAxis.h>
//random
#include <TRandom.h>
#include <TStopwatch.h>
#include <ctime>
//private setup
#include "../utils.h"

// global vars definition
const char* trigcap("");
const char* evSelCut("");
const char* evSelCutCap("");
const char* mcCut("");
TString lumiCut;

void compareTwo(TTree* t1=0 ,TTree* t2=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, TCut cut1="(1)", TCut cut2="(1)", const string cap = "");
void compare_data_mc_1D(const char* fname_data="root://eoscms//eos/cms/store/group/phys_heavyions/kjung/ExpressForests/v1/Merged/HiForest_run285090_Express.root", const char* fname_mc="/afs/cern.ch/work/q/qixu/public/centrality/Run2Prep/5TeV/mergeskim_pPb_EposMinBias_5TeV_8022_Forest_corrCentrTable_161108_122954.root")
{
/*
  TFile* f1 = TFile::Open(fname_data);
  TTree* t1 = (TTree*) f1 -> Get("hiEvtAnalyzer/HiTree");
  TTree* t1_skim = (TTree*) f1 -> Get("skimanalysis/HltTree");
  TTree* t1_hlt = (TTree*) f1 -> Get("hltanalysis/HltTree");
*/
	TChain* t1 = new TChain("hiEvtAnalyzer/HiTree");
	TChain* t1_skim = new TChain("skimanalysis/HltTree");
	TChain* t1_hlt = new TChain("hltanalysis/HltTree");

	for(Int_t i = 0; i < 19; i++)
	{
		t1->Add(Form("%s%d.root", fname_data, i));
		t1_skim->Add(Form("%s%d.root", fname_data, i));
		t1_hlt->Add(Form("%s%d.root", fname_data, i));
	}

  t1->AddFriend(t1_skim);
  t1->AddFriend(t1_hlt);

  TFile* f2 = TFile::Open(fname_mc);
  TTree* t2 = (TTree*) f2 -> Get("hiEvtAnalyzer/HiTree");
  TTree* t2_skim = (TTree*) f2 -> Get("skimanalysis/HltTree");
  TTree* t2_hlt = (TTree*) f2 -> Get("hltanalysis/HltTree");
  t2->AddFriend(t2_skim);
  t2->AddFriend(t2_hlt);
  
  //total
  double hiHFMax = 300;
  double hiHFSSMax = 200;
  double hiHFSSTruncMax = 80;
  double hiBinMax = 100;
  double hiHFhitMax = 8000;
  double hiETMax = 50;
  double hiEEMax = 50;
  double hiEBMax = 100;
  double hiNpixMax = 2000;
  double hiNpixelTracksMax = 500;
  double hiNtracksMax = 400;
  double hiNtracksCutEtaMax = 350;
  double hiNtracksCutMax = 150;
  double hiZDCMax = 40000;
  
  int nBin = 100;
   const char* trigcut = "(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1)";
  //const char* trigcut = "HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_ForExpress_v1";
  trigcap = "HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part*";
  //trigcap = "HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_ForExpress_v1";
  evSelCut = "pBeamScrapingFilter && pPAprimaryVertexFilter && phfCoincFilter1 && pVertexFilterCutG";
	lumiCut = "(lumi > 60 && lumi < 200) || (lumi > 270 && lumi < 360)";
  //const string cap = "pAExpress_5TeV_run285090";
  const string cap = "pAPromptReco_5TeV_run285090";
>>>>>>> cf422d9785422a5895839bf37f5302c36ce1e606
  evSelCutCap = "BS+PV+HFC+PVG";
  TCut selCut = Form("%s && %s && %s",trigcut,evSelCut, lumiCut.Data());
  cout << "Aplying the following event selection: " << selCut.GetTitle() << endl;
  
  mcCut = ""; //"pcollisionEventSelection";
  compareTwo(t1, t2, "hiHF",nBin,0,hiHFMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiHFplus",nBin,0,hiHFSSMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiHFminus",nBin,0,hiHFSSMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiHFplusEta4",nBin,0,hiHFSSTruncMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiHFminusEta4",nBin,0,hiHFSSTruncMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiBin",nBin,0,hiBinMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiHFhit",nBin,0,hiHFhitMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiET",nBin,0,hiETMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiEB",nBin,0,hiEBMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiEE",nBin,0,hiEEMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiNpix",hiNpixMax/10,0,hiNpixMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiNtracks",hiNtracksMax/4,0,hiNtracksMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiNtracksPtCut",hiNtracksCutMax/4,0,hiNtracksCutMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiNtracksEtaCut",hiNtracksCutEtaMax/4,0,hiNtracksCutEtaMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiNtracksEtaPtCut",hiNtracksCutMax/4,0,hiNtracksCutMax,selCut,mcCut, cap);
  //  compareTwo(t1, t2, "hiNpixelTracks",nBin,0,hiNpixelTracksMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiZDC",nBin,0,hiZDCMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiZDCplus",nBin,0,hiZDCMax,selCut,mcCut, cap);
  compareTwo(t1, t2, "hiZDCminus",nBin,0,hiZDCMax,selCut,mcCut, cap);
  
} // main function

void compareTwo(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, TCut cut1, TCut cut2, const string cap)  {
  SetHistTitleStyle();
  SetyjPadStyle();
  gStyle->SetOptStat(0);
  static int j = 0;
  TCanvas* c=  new TCanvas(Form("c_%s_%d",var.Data(),j),"", 500,900);
  c->Divide(1,2);
  c->cd(1);
  gPad->SetLogy();
  TH1D* h1 = new TH1D(Form("h1_%s_%d",var.Data(),j), Form(";%s;",var.Data()), nBins,xMin,xMax);
  TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%d",var.Data(),j));
  h1->Sumw2();
  h2->Sumw2();
  t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), cut1);
  t2->Draw(Form("%s>>%s",var.Data(),h2->GetName()), cut2);
  h1->Scale( 1. / t1->GetEntries(cut1));
  h2->Scale( 1. / t2->GetEntries(cut2));
  SetHistColor(h1,2);
  SetHistColor(h2,1);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.8);
  h1->SetNdivisions(505);
  h2->SetNdivisions(505);
  TLegend* l1 = new TLegend(0.5,0.6,0.90,0.8);
  legStyle(l1);
  l1->AddEntry(h1, "Data", "p");
  l1->AddEntry(h2, "EPOS", "l");
  
  double range = cleverRange(h1,h2,1.5,1.e-4);
  h1->DrawCopy("L");
  h2->DrawCopy("hist same");
  l1->Draw("same");
  
  c->cd(2);
  h1->Divide(h2);
  h1->SetYTitle("DATA / EPOS");
  h1->GetYaxis()->SetRangeUser(0.0,2.0);
  SetHistColor(h1,2);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.8);
  h1->DrawCopy("");
  //h2->DrawCopy("le1");
  jumSun(xMin,1,xMax,1);
  drawText(cap.data(),0.2,0.2+0.7);
  //drawText(cut1.GetTitle(),0.2,0.2);
  drawText(Form("%s",trigcap),0.2,0.2+0.12);
  drawText(Form("%s",evSelCutCap),0.2,0.2+0.06);
  drawText(Form("MC sel: %s",mcCut),0.2,0.2);
  
  // save canvas
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/comparisonDataMCDistributions/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  c-> SaveAs(Form("comparisonDataMCDistributions/compare_%s_%s.pdf",var.Data(),cap.data()));
  j++;
}
