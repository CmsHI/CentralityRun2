// computes level of noise rejection by different event filters
// Author : Javier Martin Blanco 24/01/2018

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
#include "./ebDistributions.C"

void noiseRejection(bool recreateDistributions = false,
                    const char* inputFile = "inputFiles.txt",
                    const char* triggerName = "HLT_HIL1NotBptxOR_v1",
                    const char* colLabel = "XeXe 5.44 TeV"
                    )
{
  
  // Run macro to create the noise distributions
  if (recreateDistributions) ebDistributions(inputFile,triggerName);
  //
  
  // Read file with noise distributions
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string resultsDIR = mainDIR + "/Results/";
  TFile *f = TFile::Open(Form("%s/ebDistributions.root",resultsDIR.c_str()));
  if (!f) {cout << "[ERROR] No file ebDistributions.root was found" << endl; return;}
  //

  // Get array with histos
  TObjArray* arrHistosEB = static_cast<TObjArray*>(f->FindObjectAny("ebDistributions"));
  if (!arrHistosEB) {cout << "ERROR: No histos array found in file " << f << endl; return;}
  //
  
  TObjArray* aSave = new TObjArray();
  
  const int nvar = 8;
  const char* varname[8] = {"hiHF", "hiHFplus", "hiHFplusEta4", "hiHFminus", "hiHFminusEta4","hiNtracks","hiHFhit","hiNpix"};
  
  // Create directory to save results
  string saveDIR = mainDIR + "/Results/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  
  // Make the plots
  TObjArray* arrVarEB(0x0);
  TH1* hEB(0x0);
  TH1* hEB_PV(0x0);
  TH1* hEB_BS(0x0);
  TH1* hEB_CC(0x0);
  for (int i = 0 ; i < nvar ; i++)
  {
    TString varN(varname[i]);
    
    arrVarEB = static_cast<TObjArray*>(arrHistosEB->FindObject(Form("histos_%s",varN.Data())));
    
    hEB = static_cast<TH1*>(arrVarEB->FindObject(Form("hEB_%s",varN.Data()))->Clone());
    hEB_PV = static_cast<TH1*>(arrVarEB->FindObject(Form("hEB_PV_%s",varN.Data()))->Clone());
    hEB_BS = static_cast<TH1*>(arrVarEB->FindObject(Form("hEB_BS_%s",varN.Data()))->Clone());
    hEB_CC = static_cast<TH1*>(arrVarEB->FindObject(Form("hEB_CC_%s",varN.Data()))->Clone());
    
    if (!hEB || !hEB_PV || !hEB_BS || !hEB_CC) {cout << "[ERROR]: No histos found in array for variable " << varN.Data() << endl; return;}
    
    double nEB = hEB->GetEntries();
    double nEB_PV = hEB_PV->GetEntries();
    double nEB_BS = hEB_BS->GetEntries();
    double nEB_CC = hEB_CC->GetEntries();
    
    double rej_PV = ((nEB-nEB_PV)/nEB)*100.;
    double rej_BS = ((nEB-nEB_BS)/nEB)*100.;
    double rej_CC = ((nEB-nEB_CC)/nEB)*100.;
    
    gStyle -> SetOptStat(0);
    TCanvas* c=  new TCanvas(Form("c_%s",varN.Data()),"", 900,500);
    c->Divide(2,1);
    c->cd(1);
    gPad->SetLogy();
    hEB->GetXaxis()->SetRangeUser(0,20);
    hEB->SetLineColor(1);
    hEB->SetMarkerColor(1);
    hEB->SetMarkerStyle(20);
    hEB->DrawCopy();

    TLegend* l1 = new TLegend(0.40,0.64,1.00,0.76);
    legStyle(l1);
    
    l1->AddEntry(hEB, triggerName, "p");
    l1->Draw("same");
    drawText(colLabel,0.2,0.2+0.65);
    
    c->cd(2);
    hEB_PV->GetXaxis()->SetRangeUser(0,20);
    hEB_PV->SetLineColor(1);
    hEB_PV->SetMarkerColor(1);
    hEB_PV->SetMarkerStyle(20);
    hEB_PV->Draw();
    TLegend* l2 = new TLegend(0.098,0.7,0.7,0.73);
    legStyle(l2);
    l2->AddEntry(hEB_PV, Form("%s + vtx filter",triggerName), "p");
    drawText(Form("vtx filter rejected %s = %0.3f %%",triggerName,rej_PV),0.15,0.3);
    drawText(Form("bs filter rejected %s = %0.3f %%",triggerName,rej_BS),0.15,0.25);
    drawText(Form("cc filter rejected %s = %0.3f %%",triggerName,rej_CC),0.15,0.2);
    
    l2->Draw("same");
    
    c->SaveAs(Form("%s/noiseRejection_%s.pdf",saveDIR.c_str(),varN.Data()));
    aSave->Add(c);
  }
  
  // save results
  TFile *fSave = new TFile(Form("%s/noiseRejection.root",saveDIR.c_str()),"RECREATE");
  aSave->Write("noiseRejection", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close();
  f->Close();
  
  return;
}
