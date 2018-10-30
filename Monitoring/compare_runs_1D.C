// Author : Yeonju Go
#include "../utils.h"
#include <string>
#include <iostream>
#include "TSystem.h"

// global vars definition
const char* trigcap("");
const char* evSelCut("");
const char* evSelCutCap("");

const int colhere[] = {1,2,4,kYellow+2,6,kGreen+2,46,kOrange,kViolet,kOrange+10};
void draw_compare_runs_1D(string var, TCut trigCut, const string cap, bool presclOn = true);

void compare_runs_1D()
{
//  const char* trigcut = "(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1)";
  const char* trigcut = "HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_ForExpress_v1";
  //  trigcap = "HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part*";
  trigcap = "HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_ForExpress_v1";
  evSelCut = "pBeamScrapingFilter && pPAprimaryVertexFilter && phfCoincFilter1 && pVertexFilterCutG";
  const string cap = "pAExpress_5TeV";
  evSelCutCap = "BS+PV+HFC+PVG";
  TCut selCut = Form("%s && %s",trigcut,evSelCut);
  cout << "Aplying the following event selection: " << selCut.GetTitle() << endl;
  
  TStopwatch timer;
  timer.Start();
  
  draw_compare_runs_1D("hiBin",selCut,cap);
  draw_compare_runs_1D("hiHF",selCut,cap);
  draw_compare_runs_1D("hiHFplus",selCut,cap);
  draw_compare_runs_1D("hiHFminus",selCut,cap);
  draw_compare_runs_1D("hiHFplusEta4",selCut,cap);
  draw_compare_runs_1D("hiHFminusEta4",selCut,cap);
  draw_compare_runs_1D("hiNpix",selCut,cap);
  draw_compare_runs_1D("hiNtracks",selCut,cap);
  draw_compare_runs_1D("hiNtracksPtCut",selCut,cap);
  draw_compare_runs_1D("hiNtracksEtaCut",selCut,cap);
  draw_compare_runs_1D("hiNtracksEtaPtCut",selCut,cap);
  draw_compare_runs_1D("sumPtVtx",selCut,cap);
  draw_compare_runs_1D("hiET",selCut,cap);;
  draw_compare_runs_1D("hiEE",selCut,cap);
  draw_compare_runs_1D("hiEB",selCut,cap);
  draw_compare_runs_1D("hiZDC",selCut,cap);
  draw_compare_runs_1D("hiZDCplus",selCut,cap);
  draw_compare_runs_1D("hiZDCminus",selCut,cap);
  draw_compare_runs_1D("zVtx",selCut,cap);//change many things !!
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
}

void draw_compare_runs_1D(string var, TCut trigCut, const string cap, bool presclOn){
  SetHistTitleStyle();
  SetyjPadStyle();
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  
  ////////////////////// for HIMinimumBias1
  const char* runPrefix = "284";
//  const char* dir = Form("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/ExpressForests/v1/000/%s/",runPrefix);
  //const char* dir = "root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForestAODPbPb5TeV/HIMinimumBias1/";//EOS

  
//  string run[] =
//  {"755","756","757"};

  string run[] =
  {"090","216"};
  
  const int Nrun = sizeof(run)/sizeof(string);
  const int Ncomp = 6;
  int Nset =0;
  if(Nrun%Ncomp!=0) Nset = (int)Nrun/Ncomp + 1;
  else Nset = (int)Nrun/Ncomp;
  
  const char* dir[Nrun] = {"root://eoscms//eos/cms/store/group/phys_heavyions/kjung/ExpressForests/v1/Merged/HiForest_run285090_Express.root","root://eoscms//eos/cms/store/group/phys_heavyions/mverweij/forest/ExpressForests/v1/000/285/216/chunks/HiForest_0.root"};
  
  int nBin = 50;
//  if(var=="hiBin") nBin=200;
  double xmax = 6000;
  if(var=="hiBin") xmax=200;
  if(var == "hiNpix") xmax = 2000;
  if(var == "hiHFhit") xmax = 8000;
  if(var == "hiNtracks") xmax = 400;
  if(var == "hiNtracksPtCut" || var == "hiNtracksEtaCut" || var == "hiNtracksEtaPtCut") xmax = 350;
  if(var == "sumPtVtx") xmax = 500;
  if(var == "hiEE") xmax = 50;
  if(var == "hiET") xmax = 50;
  if(var == "hiEB") xmax = 100;
  if(var == "hiZDC" || var == "hiZDCplus" || var == "hiZDCminus") xmax = 100;
  if(var == "hiHF") xmax = 300;
  if(var == "hiHFplus" || var == "hiHFminus") xmax = 200;
  if(var == "hiHFplusEta4" || var == "hiHFminusEta4") xmax = 80;
  
  // Create dir to save canvases
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/runEvolution/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  
  int nEvents[Nrun];
  
  TH1D* h[Nrun];
  for(int i=0; i<Nrun; ++i){
    TFile* f1;
//    f1 = TFile::Open(Form("%sEventTree_PbPb_data_%s_run%s_15Feb.root",dir,dataset,run[i].data()));
//    f1 = TFile::Open(Form("%s%s/HiForest_0.root",dir,run[i].data()));
    f1 = TFile::Open(Form("%s",dir[i]));
    if(f1->IsZombie()) { cout << run[i].data() << " doesn't exist!! " << endl; continue;}
    cout << "Open file : " << f1->GetName() << endl;
    TTree* t1 = (TTree*) f1 -> Get("hiEvtAnalyzer/HiTree");
    TTree* t1_skim = (TTree*) f1 -> Get("skimanalysis/HltTree");
    TTree* t1_hlt = (TTree*) f1 -> Get("hltanalysis/HltTree");
    TTree* t1_track= (TTree*) f1 -> Get("ppTrack/trackTree");
    t1->AddFriend(t1_skim);
    t1->AddFriend(t1_hlt);
    t1->AddFriend(t1_track);
    nEvents[i] = t1->GetEntries(trigCut);
    cout << run[i].data() << " # of events : " << nEvents[i] << endl;
    if(var == "zVtx") h[i] = new TH1D(Form("h%d",i), Form(";%s (cm);Event Fraction",var.data()), nBin, -25, 25); // for zVtx
    else h[i] = new TH1D(Form("h%d",i), Form(";%s;Event Fraction",var.data()), nBin, 0, xmax);
    //draw and considering Prescl
    if(presclOn) t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),trigCut);
    else t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),trigCut);
  }
  
  TCanvas* can[Nset];
  TLegend* leg[Nset];
  for(int j=0;j<Nset;j++){
    leg[j] = new TLegend(0.7,0.65,0.95,0.95);
    legStyle(leg[j]);
    can[j] = new TCanvas(Form("can%d_%s",j,var.data()),"",500,500);
    if(var != "hiBin") can[j]->SetLogy();
    //compare between 'Ncomp' runs
    for(int i=0;i<Ncomp;i++){
      if((Ncomp*j+i)>=Nrun) continue;
      h[Ncomp*j+i]->Scale(1./nEvents[Ncomp*j+i],"width"); // this is for all!!
      SetHistColor(h[Ncomp*j+i],colhere[i]);
      h[Ncomp*j+i]->SetMarkerStyle(20);
      h[Ncomp*j+i]->SetMarkerSize(0.8);
      if(i==0) h[Ncomp*j+i]->DrawCopy("L");
      else h[Ncomp*j+i]->DrawCopy("L same");
      leg[j]->AddEntry(h[Ncomp*j+i], Form("Run %s",run[Ncomp*j+i].data()),"el");
    }
    leg[j]->Draw();
    
    drawText(cap.data(),0.2,0.2+0.7);
    drawText(trigcap,0.2,0.2+0.06);
    drawText(evSelCutCap,0.2,0.2);
    can[j]->SaveAs(Form("runEvolution/compareBtwRuns_%s_run%s%s_%s.pdf",cap.data(),runPrefix,run[Ncomp*j].data(),var.data()));
  }
  cout << "c"<<endl;
  TCanvas* cratio[Nset];
  for(int j=0;j<Nset;j++){
    cratio[j] = new TCanvas(Form("cratio%d_%s",j,var.data()),"",500,500);
    if(var != "hiBin") {
      can[j]->SetLogy();
      h[Ncomp*j]->Rebin(2);
    }
    //compare between 'Ncomp' runs
    for(int i=1;i<Ncomp;i++){
      if((Ncomp*j+i)>=Nrun) continue;
      if(var != "hiBin") h[Ncomp*j+i]->Rebin(2);
      h[Ncomp*j+i]->Divide(h[Ncomp*j]);
      h[Ncomp*j+i]->SetMarkerStyle(20);
      h[Ncomp*j+i]->SetMarkerSize(0.85);
      h[Ncomp*j+i]->GetYaxis()->SetRangeUser(0.5,1.5);
      h[Ncomp*j+i]->SetTitle(Form(";%s;Run XXX / Run %s%s",var.data(),runPrefix,run[Ncomp*j].data()));
      if(i==1) h[Ncomp*j+i]->DrawCopy("L");
      else h[Ncomp*j+i]->DrawCopy("L same");
    }
    leg[j]->Draw();
    drawText(cap.data(),0.2,0.2+0.7);
    drawText(trigcap,0.2,0.2+0.06);
    drawText(evSelCutCap,0.2,0.2);
    if(var == "zVtx") jumSun(-25,1,25,1);//for zVtx
    else jumSun(0,1,xmax,1);
    cratio[j]->SaveAs(Form("runEvolution/compareBtwRuns_%s_ratio_run%s%s_%s.pdf",cap.data(),runPrefix,run[Ncomp*j].data(),var.data()));
  }
}

