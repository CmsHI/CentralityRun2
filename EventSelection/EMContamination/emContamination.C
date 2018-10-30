// compute fraction of EM events after event selection
// Author : Javier Martin Blanco 25/01/2017

//basic c++ header, string ...
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <cstring>
#include <vector>
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
#include <TF1.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>
#include <TCut.h>
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
#include "../../Utilities/drawUtils.h"
#include "./mbDistributions.C"

void emContamination(const char* hfFilter = "2", // Number of required towers in HF coincidence filter
                     const std::vector<bool> recreateDistributions = {false,false,false}, // {data,MC,EM}
                     const std::vector<const char*> inputFile = {"dataInputFiles.txt","emInputFiles_single.txt","emInputFiles_double.txt"},
                     const std::vector<const char*> sampleLabel = {"DATA","STARlight_single","STARlight_double"},
                     const std::vector<double> xSection = {1,0.2,0.2},
                     const char* colLabel = "XeXe 5.44 TeV",
                     const char* triggerName = "HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1"
                 )
{
  
  // Check the input parameters
  const int nSamples = 3;
  if (strcmp(hfFilter,"1") && strcmp(hfFilter,"2") && strcmp(hfFilter,"3") && strcmp(hfFilter,"4")) {cout << "[ERROR] hfFilter must be the number of towers for the filter. Options are: 1 , 2 , 3 or 4 " << endl; return;}
  if (recreateDistributions.size()!=nSamples) {cout << "[ERROR] recreateDistributions must have " << nSamples << " entries. Usage: {data,starlight_single,starlight_double}" << endl; return;}
  if (inputFile.size()!=nSamples) {cout << "[ERROR] inputfile must have " << nSamples << " entries. Usage: {data,starlight_single,starlight_double}" << endl; return;}
  if (sampleLabel.size()!=nSamples) {cout << "[ERROR] sampleLabel must have " << nSamples << " entries. Usage: {data,starlight_single,starlight_double}" << endl; return;}
  if (xSection.size()!=nSamples) {cout << "[ERROR] xSection must have 3 entries. Usage: {data,starlight_single,starlight_double}" << endl; return;}
  
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string resultsDIR = mainDIR + "/Results";
  
  TObjArray* arrHistosarr = new TObjArray();
  for (int n = 0 ; n < nSamples ; n++)
  {
    // Run macro to create the mb distributions
    if (recreateDistributions[n])
    {
      cout << "[INFO] Creating distributions for " << sampleLabel[n]  << endl;
      mbDistributions(inputFile[n],triggerName,sampleLabel[n]);
    }
    
    // Read file with mb distributions
    TFile *f = TFile::Open(Form("%s/mbDistributions_%s.root",resultsDIR.c_str(),sampleLabel[n]));
    if (!f) {cout << "[ERROR] No file mbDistributions_" << sampleLabel[n] << ".root was found" << endl; return;}
    //
    // Get array with histos
    TObjArray* arrHistos = static_cast<TObjArray*>(f->FindObjectAny("mbDistributions"));
    if (!arrHistos) {cout << "ERROR: No histos array found in file " << f->GetName() << endl; return;}
    arrHistosarr->Add(arrHistos);
    //
  }
  
  TObjArray* aSave = new TObjArray();
  
  const int nvar = 10;
  const char* varname[nvar] = {"hiHF", "hiHFplus", "hiHFplusEta4", "hiHFminus",
    "hiHFminusEta4","hiNtracks","hiHFhit","hiNpix","nTrk","nTrkVtx"};
  
  // Create directory to save the plots
  string saveDIR = mainDIR + "/Results/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  
  // Define cross section weigthing functions
  TF1* fEMS = new TF1("wEMS","[0]",0.,1000.);
  TF1* fEMD = new TF1("wEMD","[0]",0.,1000.);
  
  // Make the plots
  TObjArray* arrVar(0x0);
  TH1* hMB(0x0);
  TH1* hMB_evtSel(0x0);
  for (int i = 0 ; i < nvar ; i++)
  {
    TString varN(varname[i]);
    
    gStyle -> SetOptStat(0);
    TCanvas* c=  new TCanvas(Form("c_%s",varN.Data()),"", 900,500);
    c->Divide(2,1);
    
    TLegend* l1 = new TLegend(0.40,0.64,1.00,0.76);
    legStyle(l1);
    
    double nMB(0.),nEMS(0.),nEMD(0.);
    double nMB_evtSel(0.),nEMS_evtSel(0.),nEMD_evtSel(0.);
    double wEMS(0.),wEMD(0.);
    
    for (int n = 0 ; n < nSamples ; n++)
    {
      TObjArray* arrHistos = static_cast<TObjArray*>(arrHistosarr->At(n));
      if (arrHistos) arrVar = static_cast<TObjArray*>(arrHistos->FindObject(Form("histos_%s",varN.Data())));
      else {cout << "[ERROR] No array of histos found for " << sampleLabel[n] << endl; return;}
      
      hMB = static_cast<TH1*>(arrVar->FindObject(Form("hMB_%s",varN.Data()))->Clone());
      hMB_evtSel = static_cast<TH1*>(arrVar->FindObject(Form("hMB_BS_PV_HF%s_%s",hfFilter,varN.Data()))->Clone());
      if (!hMB || !hMB_evtSel)
      {
        cout << "[ERROR] No MB histos found in array " << endl;
        return;
      }
      
      if (n==0)
      {
        nMB = hMB->GetEntries();
        nMB_evtSel = hMB_evtSel->GetEntries();
      }
      else if (n==1)
      {
        nEMS = hMB->GetEntries();
        nEMS_evtSel = hMB_evtSel->GetEntries();
        wEMS = (nMB/xSection[0])*(xSection[1]/nEMS);
        fEMS->SetParameter(0,wEMS);
        hMB->Multiply(fEMS);
        hMB_evtSel->Multiply(fEMS);
      }
      else if (n==2)
      {
        nEMD = hMB->GetEntries();
        nEMD_evtSel = hMB_evtSel->GetEntries();
        wEMD = (nMB/xSection[0])*(xSection[2]/nEMD);
        fEMD->SetParameter(0,wEMD);
        hMB->Multiply(fEMD);
        hMB_evtSel->Multiply(fEMD);
      }
      else {cout << "[ERROR] There are more samples than expected" << endl; return;}
      
      c->cd(1);
      gPad->SetLogy();
      hMB->GetXaxis()->SetRangeUser(0,50);
      hMB->SetLineColor(n+1);
      hMB->SetMarkerColor(n+1);
      hMB->SetMarkerStyle(20);
      hMB->DrawCopy(n==0?"":"same");
      
      l1->AddEntry(hMB, sampleLabel[n], "p");
      
      c->cd(2);
      hMB_evtSel->GetXaxis()->SetRangeUser(0,50);
      hMB_evtSel->SetLineColor(n+1);
      hMB_evtSel->SetMarkerColor(n+1);
      hMB_evtSel->SetMarkerStyle(20);
      hMB_evtSel->Draw(n==0?"":"same");
    }
    // Compute EM contamination
    double emCont = ((wEMS*nEMS)+(wEMD*nEMD))*100./(nMB);
    double emCont_evtSel = ((wEMS*nEMS_evtSel)+(wEMD*nEMD_evtSel))*100./(nMB_evtSel);
    
    c->cd(1);
    l1->Draw("same");
    drawText(colLabel,0.54,0.2+0.65);
    drawText("MB triggered",0.54,0.2+0.60);
    drawText(Form("EM contamination in HF%s = %0.3f %%",hfFilter,emCont),0.15,0.15);
    
    c->cd(2);
    drawText(Form("MB + vtx +BS + HF%s filters",hfFilter),0.15,0.2+0.65);
    drawText(Form("EM contamination after evt. selection (HF%s) = %0.3f %%",hfFilter,emCont_evtSel),0.15,0.2+0.60);
    
    c->SaveAs(Form("%s/emContamination_%s_HF%s.pdf",saveDIR.c_str(),varN.Data(),hfFilter));
    aSave->Add(c);
  }

  // save results
  TFile *fSave = new TFile(Form("%s/emContamination_HF%s.root",saveDIR.c_str(),hfFilter),"RECREATE");
  aSave->Write("emContamination", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close();
}
