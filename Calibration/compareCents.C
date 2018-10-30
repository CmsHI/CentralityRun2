// comparison centrality bin values using two different calibrations
// IMPORTANT NOTE: the two forest must contain the same number of events, and they must have been produced locally using the same data files so the events and the ordering in the tree are the same
// Author : Javier Martin Blanco 15/02/2017

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
#include <TH1.h>
#include <TH1D.h>

//private setup
#include "../Utilities/treeUtils.h"

//////////////////////////////////////////////
// Definition of global variables
//////////////////////////////////////////////

int hiBin1,hiBin2;
UInt_t run, lumi;
int HLT_Trigger_part1, HLT_Trigger_part2, HLT_Trigger_part3, HLT_Trigger_part4, HLT_Trigger_part5, HLT_Trigger_part6, HLT_Trigger_part7, HLT_Trigger_part8, HLT_Trigger_part9, HLT_Trigger_part10, HLT_Trigger_part11, HLT_Trigger_part12, HLT_Trigger_part13, HLT_Trigger_part14, HLT_Trigger_part15, HLT_Trigger_part16, HLT_Trigger_part17, HLT_Trigger_part18, HLT_Trigger_part19, HLT_Trigger_part20;
int pBeamScrapingFilter, pprimaryVertexFilter, phfCoincFilter3;

bool loadSelBranches(const char* triggerName, TChain* t);

void compareCents(const char* inputFile1 = "inputFiles1.txt",
                  const char* inputFile2 = "inputFiles2.txt",
                  const char* triggerName = "HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1",
                    )
{
  //
  
  // Read files and chain TTrees
  TChain* tree1 = readFiles(inputFile1);
  TChain* tree2 = readFiles(inputFile2);
  
  tree1->SetBranchAddress("hiBin",	&hiBin1);
  tree2->SetBranchAddress("hiBin",	&hiBin2);
  
  loadSelBranches(triggerName,tree1);
  //
  
  unsigned int nEvents1 = tree1->GetEntries();
  unsigned int nEvents2 = tree2->GetEntries();
  
  if (nEvents1 != nEvents2) cout << "[ERROR] Trees do not contain the same number of events << endl;
  else {cout << "[ERROR] No events to be processed" << endl; return;}
  
  // Load needed branches
  if (!loadBranches(triggerName,tree)) {cout << "[ERROR] Error loading branches" << endl; return;}
  //
  
  TH1::SetDefaultSumw2(true);
  TH1D* h = new TH1D("hiBinComparison",";#Delta hiBin; N_{events}", 200,-100.,100.);
  
  // Event loop
  for(unsigned int iev = 0; iev < nEvents; iev++)
  {
    if(iev%5000 == 0) cout<<"Processing event: " << iev << " / " << nEvents << endl;
    tree->GetEntry(iev);
    
    if ((HLT_Trigger_part1 || HLT_Trigger_part2 || HLT_Trigger_part3 || HLT_Trigger_part4 || HLT_Trigger_part5 ||
         HLT_Trigger_part6 || HLT_Trigger_part7 || HLT_Trigger_part8 || HLT_Trigger_part9 || HLT_Trigger_part10 ||
         HLT_Trigger_part11 || HLT_Trigger_part12 || HLT_Trigger_part13 || HLT_Trigger_part14 || HLT_Trigger_part15 ||
         HLT_Trigger_part16 || HLT_Trigger_part17 || HLT_Trigger_part18 || HLT_Trigger_part19 || HLT_Trigger_part20) && pBeamScrapingFilter && pprimaryVertexFilter && phfCoincFilter3)
    {
      double DhiBin = hiBin1 - hiBin2;
      h->Fill(DhiBin);
    }
  }
  
  gStyle -> SetOptStat(0);
  TCanvas* c =  new TCanvas(Form("c_%s",varN.Data()),"", 500,500);
  h->Draw();
  
  c->SaveAs("hiBinComparison.pdf");
}

bool loadSelBranches(const char* triggerName, TChain* t)
{
  TString tName(triggerName);
  
  TObjArray* tArray(0x0);
  if (tName.Contains(";")) tArray = tName.Tokenize(";");
  else
  {
    tArray = new TObjArray();
    tArray->Add(new TObjString(triggerName));
  }
  tArray->SetOwner(kTRUE);
  
  TIter it(tArray);
  TObjString* obj(0x0);
  int cnt(1);
  while ((obj = static_cast<TObjString*>(it.Next())))
  {
    TString tChar = obj->GetString();
    if (!t->GetBranch(tChar.Data()))
    {
      cout << "[ERROR] Specified trigger path: " << tChar.Data() << " does not exist" << endl;
      return false;
    }
    else if (cnt==1) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part1);
    else if (cnt==2) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part2);
    else if (cnt==3) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part3);
    else if (cnt==4) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part4);
    else if (cnt==5) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part5);
    else if (cnt==6) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part6);
    else if (cnt==7) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part7);
    else if (cnt==8) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part8);
    else if (cnt==9) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part9);
    else if (cnt==10) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part10);
    else if (cnt==11) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part11);
    else if (cnt==12) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part12);
    else if (cnt==13) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part13);
    else if (cnt==14) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part14);
    else if (cnt==15) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part15);
    else if (cnt==16) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part16);
    else if (cnt==17) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part17);
    else if (cnt==18) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part18);
    else if (cnt==19) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part19);
    else if (cnt==20) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part20);
    else { cout << "[ERROR] Code does not support such a long trigger list, please update" << endl; return false;}
    
    cnt++;
  }
  
  if (t->GetBranch("pBeamScrapingFilter")) t->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter);
  if (t->GetBranch("pprimaryVertexFilter")) t->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
  if (t->GetBranch("phfCoincFilter3")) t->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3);
  
  tArray->Delete();
  return true;
}

