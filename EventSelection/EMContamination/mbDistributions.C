// comparison centrality objects between DATA 2016 and DATA 2015
// Author : Javier Martin Blanco 10/11/2016

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
#include "../../Utilities/treeUtils.h"

void mbDistributions(const char* inputFile = "inputFiles.txt",
                    const char* triggerName = "HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1",
                    const char* sampleLabel = "STARlight_single"
                    )
{
  // Read files and chain TTrees
  TChain* tree = readFiles(inputFile);
  //
  
  unsigned int nEvents = tree->GetEntries();
  if (nEvents>0) cout << "[INFO] Number of events to be processed = " << nEvents << endl;
  else {cout << "[ERROR] No events to be processed" << endl; return;}
  
  // Load needed branches
  if (!loadBranches(triggerName,tree)) {cout << "[ERROR] Error loading branches" << endl; return;}
  //
  
  
  // Definition of variable ranges and number of bins for plots
  double hiHFMax = 4000;
  double hiHFSSMax = 2000;
  double hiHFSSTruncMax = 80;
  double hiBinMax = 100;
  double hiHFhitMax = 10000;
  double hiETMax = 50;
  double hiEEMax = 50;
  double hiMBMax = 100;
  double hiNpixMax = 40000;
  double hiNpixelTracksMax = 500;
  double hiNtracksMax = 4000;
  double hiNtracksCutEtaMax = 350;
  double hiNtracksCutMax = 150;
  double hiZDCMax = 40000;
  
  int nBin = 1000;
  
  const int nvar = 19;
  const char* varname[nvar] = {"hiHF", "hiHFplus", "hiHFplusEta4", "hiHFminus",
    "hiHFminusEta4","hiNtracks","hiHFhit","hiNpix","nVtx","nTrk","trkEta","trkPt","zVtx","trkDxy1","trkDz1","nTrkVtx","highPurity","normChi2Vtx","vtxDist2D"};
  
    // Definition of histos to store distributions
  TObjArray* avarName = new TObjArray();
  for (int i = 0 ; i < nvar ; i++)
  {
    avarName->Add(new TObjString(varname[i]));
  }
  
  TObjArray* aSave = new TObjArray();
  
  TH1::SetDefaultSumw2(true);
  for (int i = 0 ; i < nvar ; i++)
  {
    TObjArray* arr = new TObjArray();
    arr->SetName(Form("histos_%s",varname[i]));
    
    double xMin = 0.;
    double xMax = 300;
    
    if( !strcmp(varname[i],"hiNtracks") ) {xMax = hiNtracksMax; nBin = (int)hiNtracksMax;}
    if( !strcmp(varname[i],"hiHF")) {xMax = hiHFMax; nBin = (int)(hiHFMax*4*10);}
    if( !strcmp(varname[i],"hiHFplus") || !strcmp(varname[i],"hiHFminus")) {xMax = hiHFSSMax; nBin = (int)(hiHFSSMax/2.);}
    if( !strcmp(varname[i],"hiHFplusEta4") || !strcmp(varname[i],"hiHFminusEta4")) {xMax = hiHFSSTruncMax; nBin = (int)hiHFSSTruncMax*8;}
    if( !strcmp(varname[i],"hiHFhit") ){xMax = hiHFhitMax; nBin = (int)(hiHFhitMax/10.);}
    if( !strcmp(varname[i],"hiNpix") ){xMax = hiNpixMax; nBin = (int)(hiNpixMax/10.);}
    if( !strcmp(varname[i],"zVtx") ){xMin = -30.; xMax = 30.; nBin = 1200;}
    if( !strcmp(varname[i],"trkDxy1") ){xMin = -150; xMax = 150.; nBin = 1200;}
    if( !strcmp(varname[i],"trkDz1") ){xMin = -600; xMax = 600.; nBin = 4800;}
    if( !strcmp(varname[i],"nVtx") ){xMax = 20; nBin = 20;}
    if( !strcmp(varname[i],"nTrkVtx") ){xMax = 500; nBin = 500;}
    if( !strcmp(varname[i],"highPurity") ){xMax = 2; nBin = 2;}
    if( !strcmp(varname[i],"nTrk") ){xMax = 600; nBin = 600;}
    if( !strcmp(varname[i],"trkEta") ){xMin = -4. ; xMax = 4.; nBin = 40;}
    if( !strcmp(varname[i],"trkPt") ){xMin = 0. ; xMax = 1000.; nBin = 500;}
    if( !strcmp(varname[i],"normChi2Vtx") ){xMax = 200; nBin = 2000;}
    if( !strcmp(varname[i],"vtxDist2D") ){xMin = -3;xMax = 3; nBin = 1920;}

    
    TH1D* hMB = new TH1D(Form("hMB_%s",varname[i]), Form(";%s;",varname[i]), nBin,xMin,xMax);
    TH1D* h = (TH1D*)hMB->Clone(Form("h_%s",varname[i]));
    TH1D* h_BS = (TH1D*)hMB->Clone(Form("h_BS_%s",varname[i]));
    TH1D* h_PV = (TH1D*)hMB->Clone(Form("h_PV_%s",varname[i]));
    TH1D* h_HF1 = (TH1D*)hMB->Clone(Form("h_HF1_%s",varname[i]));
    TH1D* h_HF2 = (TH1D*)hMB->Clone(Form("h_HF2_%s",varname[i]));
    TH1D* h_HF3 = (TH1D*)hMB->Clone(Form("h_HF3_%s",varname[i]));
    TH1D* h_HF4 = (TH1D*)hMB->Clone(Form("h_HF4_%s",varname[i]));
    TH1D* h_BS_PV = (TH1D*)hMB->Clone(Form("h_BS_PV_%s",varname[i]));
    TH1D* h_BS_PV_HF1 = (TH1D*)hMB->Clone(Form("h_BS_PV_HF1_%s",varname[i]));
    TH1D* h_BS_PV_HF2 = (TH1D*)hMB->Clone(Form("h_BS_PV_HF2_%s",varname[i]));
    TH1D* h_BS_PV_HF3 = (TH1D*)hMB->Clone(Form("h_BS_PV_HF3_%s",varname[i]));
    TH1D* h_BS_PV_HF4 = (TH1D*)hMB->Clone(Form("h_BS_PV_HF4_%s",varname[i]));
    
    TH1D* hMB_BS = (TH1D*)hMB->Clone(Form("hMB_BS_%s",varname[i]));
    TH1D* hMB_PV = (TH1D*)hMB->Clone(Form("hMB_PV_%s",varname[i]));
    TH1D* hMB_HF1 = (TH1D*)hMB->Clone(Form("hMB_HF1_%s",varname[i]));
    TH1D* hMB_HF2 = (TH1D*)hMB->Clone(Form("hMB_HF2_%s",varname[i]));
    TH1D* hMB_HF3 = (TH1D*)hMB->Clone(Form("hMB_HF3_%s",varname[i]));
    TH1D* hMB_HF4 = (TH1D*)hMB->Clone(Form("hMB_HF4_%s",varname[i]));
    TH1D* hMB_BS_PV = (TH1D*)hMB->Clone(Form("hMB_BS_PV_%s",varname[i]));
    TH1D* hMB_BS_PV_HF1 = (TH1D*)hMB->Clone(Form("hMB_BS_PV_HF1_%s",varname[i]));
    TH1D* hMB_BS_PV_HF2 = (TH1D*)hMB->Clone(Form("hMB_BS_PV_HF2_%s",varname[i]));
    TH1D* hMB_BS_PV_HF3 = (TH1D*)hMB->Clone(Form("hMB_BS_PV_HF3_%s",varname[i]));
    TH1D* hMB_BS_PV_HF4 = (TH1D*)hMB->Clone(Form("hMB_BS_PV_HF4_%s",varname[i]));
    
    TH1D* hMB_NotBS = (TH1D*)hMB->Clone(Form("hMB_NotBS_%s",varname[i]));
    TH1D* hMB_NotPV = (TH1D*)hMB->Clone(Form("hMB_NotPV_%s",varname[i]));
    TH1D* hMB_NotHF1 = (TH1D*)hMB->Clone(Form("hMB_NotHF1_%s",varname[i]));
    TH1D* hMB_NotHF2 = (TH1D*)hMB->Clone(Form("hMB_NotHF2_%s",varname[i]));
    TH1D* hMB_NotHF3 = (TH1D*)hMB->Clone(Form("hMB_NotHF3_%s",varname[i]));
    TH1D* hMB_NotHF4 = (TH1D*)hMB->Clone(Form("hMB_NotHF4_%s",varname[i]));
    
    arr->Add(h);
    arr->Add(h_BS);
    arr->Add(h_PV);
    arr->Add(h_HF1);
    arr->Add(h_HF2);
    arr->Add(h_HF3);
    arr->Add(h_HF4);
    arr->Add(h_BS_PV);
    arr->Add(h_BS_PV_HF1);
    arr->Add(h_BS_PV_HF2);
    arr->Add(h_BS_PV_HF3);
    arr->Add(h_BS_PV_HF4);
    arr->Add(hMB);
    arr->Add(hMB_BS);
    arr->Add(hMB_PV);
    arr->Add(hMB_HF1);
    arr->Add(hMB_HF2);
    arr->Add(hMB_HF3);
    arr->Add(hMB_HF4);
    arr->Add(hMB_BS_PV);
    arr->Add(hMB_BS_PV_HF1);
    arr->Add(hMB_BS_PV_HF2);
    arr->Add(hMB_BS_PV_HF3);
    arr->Add(hMB_BS_PV_HF4);
    
    arr->Add(hMB_NotBS);
    arr->Add(hMB_NotPV);
    arr->Add(hMB_NotHF1);
    arr->Add(hMB_NotHF2);
    arr->Add(hMB_NotHF3);
    arr->Add(hMB_NotHF4);
    
    aSave->Add(arr);
  }
  
  TObjArray* adummy(0x0);
  // Event loop
  for(unsigned int iev = 0; iev < nEvents; iev++)
  {
    if(iev%5000 == 0) cout<<"Processing event: " << iev << " / " << nEvents << endl;
    tree->GetEntry(iev);
    if(run!=304906) continue;  
    for (int i = 0 ; i < nvar ; i++)
    {
      TString varN = static_cast<TObjString*>(avarName->At(i))->GetString();
      adummy = static_cast<TObjArray*>(aSave->FindObject(Form("histos_%s",varN.Data())));
      
      float parameter = -1;
      int loopPar = 1;
      if(!strcmp(varN.Data(),"hiHF")) parameter = hf;
      if(!strcmp(varN.Data(),"hiHFplus")) parameter = hfplus;
      if(!strcmp(varN.Data(),"hiHFminus")) parameter = hfminus;
      if(!strcmp(varN.Data(),"hiHFplusEta4")) parameter = hfpluseta4;
      if(!strcmp(varN.Data(),"hiHFminusEta4")) parameter = hfminuseta4;
      if(!strcmp(varN.Data(),"hiNtracks")) parameter = hiNtrks;
      if(!strcmp(varN.Data(),"hiHFhit")) parameter = hfhit;
      if(!strcmp(varN.Data(),"hiNpix")) parameter = npix;
      if(!strcmp(varN.Data(),"nVtx")) parameter = nvtx;
      if(!strcmp(varN.Data(),"nTrk")) parameter = ntrk;
      if(!strcmp(varN.Data(),"zVtx")) parameter = zvtx[maxptvtx];
      if(!strcmp(varN.Data(),"trkDxy1")) parameter = trkdxy1[maxptvtx];
      if(!strcmp(varN.Data(),"trkDz1")) parameter = trkdz1[maxptvtx];
      if(!strcmp(varN.Data(),"nTrkVtx")) parameter = ntrkvtx[maxptvtx];
      if(!strcmp(varN.Data(),"highPurity")) parameter = (int)highpurity[maxptvtx];
      if(!strcmp(varN.Data(),"normChi2Vtx") )parameter = normchi2vtx[maxptvtx];
      if(!strcmp(varN.Data(),"vtxDist2D") )parameter = vtxdist2d[maxptvtx];
      if(!strcmp(varN.Data(),"trkEta") ) {loopPar = ntrk;}
      if(!strcmp(varN.Data(),"trkPt") ) {loopPar = ntrk;}
      
      static_cast<TH1*>(adummy->FindObject(Form("h_%s",varN.Data())))->Fill(parameter); // Histogram without any selection
      
      for (int n = 0 ; n < loopPar ; n++ )
      {
          if(!strcmp(varN.Data(),"trkEta") ) parameter = trkEta[n];
          if(!strcmp(varN.Data(),"trkPt") )  parameter = trkPt[n];

          if (pBeamScrapingFilter)
          {
              static_cast<TH1*>(adummy->FindObject(Form("h_BS_%s",varN.Data())))->Fill(parameter);

              if (pprimaryVertexFilter)
              {
                  static_cast<TH1*>(adummy->FindObject(Form("h_BS_PV_%s",varN.Data())))->Fill(parameter);

                  if (phfCoincFilter1)
                  {
                      static_cast<TH1*>(adummy->FindObject(Form("h_BS_PV_HF1_%s",varN.Data())))->Fill(parameter);
                  }

                  if (phfCoincFilter2)
                  {
                      static_cast<TH1*>(adummy->FindObject(Form("h_BS_PV_HF2_%s",varN.Data())))->Fill(parameter);
                  }

                  if (phfCoincFilter3)
                  {
                      static_cast<TH1*>(adummy->FindObject(Form("h_BS_PV_HF3_%s",varN.Data())))->Fill(parameter);
                  }

                  if (phfCoincFilter4)
                  {
                      static_cast<TH1*>(adummy->FindObject(Form("h_BS_PV_HF4_%s",varN.Data())))->Fill(parameter);
                  }
              }
          }
      }
      
      if ((HLT_Trigger_part1 || HLT_Trigger_part2 || HLT_Trigger_part3 || HLT_Trigger_part4 || HLT_Trigger_part5 ||
           HLT_Trigger_part6 || HLT_Trigger_part7 || HLT_Trigger_part8 || HLT_Trigger_part9 || HLT_Trigger_part10 ||
           HLT_Trigger_part11 || HLT_Trigger_part12 || HLT_Trigger_part13 || HLT_Trigger_part14 || HLT_Trigger_part15 ||
           HLT_Trigger_part16 || HLT_Trigger_part17 || HLT_Trigger_part18 || HLT_Trigger_part19 || HLT_Trigger_part20))
      {
        
        for (int n = 0 ; n < loopPar ; n++ )
        {
          if(!strcmp(varN.Data(),"trkEta") ) parameter = trkEta[n];
          if(!strcmp(varN.Data(),"trkPt") )  parameter = trkPt[n];
          
          static_cast<TH1*>(adummy->FindObject(Form("hMB_%s",varN.Data())))->Fill(parameter);
          
          if (pBeamScrapingFilter)
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_BS_%s",varN.Data())))->Fill(parameter);
            
            if (pprimaryVertexFilter)
            {
              static_cast<TH1*>(adummy->FindObject(Form("hMB_BS_PV_%s",varN.Data())))->Fill(parameter);
              
              if (phfCoincFilter1)
              {
                static_cast<TH1*>(adummy->FindObject(Form("hMB_BS_PV_HF1_%s",varN.Data())))->Fill(parameter);
              }
              
              if (phfCoincFilter2)
              {
                static_cast<TH1*>(adummy->FindObject(Form("hMB_BS_PV_HF2_%s",varN.Data())))->Fill(parameter);
              }
              
              if (phfCoincFilter3)
              {
                static_cast<TH1*>(adummy->FindObject(Form("hMB_BS_PV_HF3_%s",varN.Data())))->Fill(parameter);
              }
              
              if (phfCoincFilter4)
              {
                static_cast<TH1*>(adummy->FindObject(Form("hMB_BS_PV_HF4_%s",varN.Data())))->Fill(parameter);
              }
            }
          }
          else
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_NotBS_%s",varN.Data())))->Fill(parameter);
          }
          
          if (pprimaryVertexFilter)
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_PV_%s",varN.Data())))->Fill(parameter);
          }
          else
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_NotPV_%s",varN.Data())))->Fill(parameter);
          }
          
          
          if (phfCoincFilter1)
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_HF1_%s",varN.Data())))->Fill(parameter);
          }
          else
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_NotHF1_%s",varN.Data())))->Fill(parameter);
          }
          
          if (phfCoincFilter2)
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_HF2_%s",varN.Data())))->Fill(parameter);
          }
          else
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_NotHF2_%s",varN.Data())))->Fill(parameter);
          }
          
          if (phfCoincFilter3)
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_HF3_%s",varN.Data())))->Fill(parameter);
          }
          else
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_NotHF3_%s",varN.Data())))->Fill(parameter);
          }
          
          if (phfCoincFilter4)
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_HF4_%s",varN.Data())))->Fill(parameter);
          }
          else
          {
            static_cast<TH1*>(adummy->FindObject(Form("hMB_NotHF4_%s",varN.Data())))->Fill(parameter);
          }
        }
      }
    }
  }
  
  // save results
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/Results";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  TFile *f = new TFile(Form("%s/mbDistributions_run304906_%s.root",saveDIR.c_str(),sampleLabel),"RECREATE");
  //TFile *f = new TFile(Form("%s/mbDistributions_%s.root",saveDIR.c_str(),sampleLabel),"RECREATE");
  aSave->Write("mbDistributions", TObject::kOverwrite | TObject::kSingleKey);
  f->Close();
}
