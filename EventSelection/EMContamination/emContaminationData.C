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
#include <TEfficiency.h>
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

TEfficiency* getMCeff(const char* varname,
                      const char* mcResults,
                      const char* hfFilter = "2");

bool computeMCeff(const char* mcResults,
                  const char* hfFilter = "2");

void emContaminationData(const char* hfFilter = "2", // Number of required towers in HF coincidence filter
                         const std::vector<bool> recreateDistributions = {false,false}, // {data,MC,EM}
                         const std::vector<const char*> inputFile = {"dataInputFiles.txt","hydjetInputFiles.txt"},
                         const std::vector<const char*> sampleLabel = {"DATA","HYDJET"},
                         const char* colLabel = "XeXe 5.44 TeV",
                         const char* triggerName = "HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1")
{
  const int nSamples = 2;
  if (strcmp(hfFilter,"1") && strcmp(hfFilter,"2") && strcmp(hfFilter,"3") && strcmp(hfFilter,"4")) {cout << "[ERROR] hfFilter must be the number of towers for the filter. Options are: 1 , 2 , 3 or 4 " << endl; return;}
  if (recreateDistributions.size()!=nSamples) {cout << "[ERROR] recreateDistributions must have " << nSamples << " entries. Usage: {data,MC}" << endl; return;}
  if (inputFile.size()!=nSamples) {cout << "[ERROR] inputfile must have " << nSamples << " entries. Usage: {data,MC}" << endl; return;}
  if (sampleLabel.size()!=nSamples) {cout << "[ERROR] sampleLabel must have " << nSamples << " entries. Usage: {data,MC}" << endl; return;}
  
  for (int n = 0 ; n < nSamples ; n++)
  {
    // Run macro to create the mb distributions
    if (recreateDistributions[n])
    {
      cout << "[INFO] Creating distributions for " << sampleLabel[n]  << endl;
      mbDistributions(inputFile[n],triggerName,sampleLabel[n]);
    }
    
  }
  
  // Create directory to save the plots
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/Results";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  
  TString shfFilter(hfFilter);
  TString hfFilterPlus("3");
  if (!strcmp(hfFilter,"1")) hfFilterPlus = "2";
  else if (!strcmp(hfFilter,"2")) hfFilterPlus = "3";
  else if (!strcmp(hfFilter,"3")) hfFilterPlus = "4";
  else {cout << "[ERROR] Cannot compute EM contamination for hfFilter > 3" << endl; return;}
  
  const int nvar = 6;
  const char* varname[nvar] = {"hiHF", "hiHFplus", "hiHFminus","hiNtracks","hiHFhit","hiNpix"};
  
  // Read file with mb data distributions
  TFile *f = TFile::Open(Form("%s/mbDistributions_%s.root",saveDIR.c_str(),sampleLabel[0]));
  if (!f) {cout << "[ERROR] No file mbDistributions_" << sampleLabel[0] << ".root was found" << endl; return;}
  
  // Get array with histos
  TObjArray* arrHistos = static_cast<TObjArray*>(f->FindObjectAny("mbDistributions"));
  if (!arrHistos) {cout << "[ERROR] No histos array found in file " << f->GetName() << endl; return;}
  //
  
  TObjArray* aSave = new TObjArray();
  TH1* hMB_evtSel(0x0);
  TH1* hMB_evtSel_corr(0x0);
  TH1* hMB_evtSel_corrPlus(0x0);
  for (int i = 0 ; i < nvar ; i++)
  {
    TString varN(varname[i]);
    TObjArray* arrVar = static_cast<TObjArray*>(arrHistos->FindObject(Form("histos_%s",varN.Data())));
    if (!arrVar) {cout << "[ERROR] No histos array found for var " << varN.Data() << " found in file mbDistributions_" << sampleLabel[0] << ".root" << endl; return;}
    
    TEfficiency* eff = getMCeff(varN.Data(),sampleLabel[1],shfFilter.Data());
    TEfficiency* effPlus = getMCeff(varN.Data(),sampleLabel[1],hfFilterPlus.Data());
    
    if (!eff || !effPlus) {cout << "[ERROR] No efficiencies found" << endl; return;}
    
    hMB_evtSel = static_cast<TH1*>(arrVar->FindObject(Form("hMB_BS_PV_HF%s_%s",shfFilter.Data(),varN.Data())));
    if (!hMB_evtSel) {cout << "[ERROR] No histo " << Form("hMB_BS_PV_HF%s_%s",shfFilter.Data(),varN.Data()) << " found" << endl; return;}
    hMB_evtSel_corr = static_cast<TH1*>(hMB_evtSel->Clone(Form("hMB_BS_PV_HF%s_%s_corr",shfFilter.Data(),varN.Data())));
    hMB_evtSel = static_cast<TH1*>(arrVar->FindObject(Form("hMB_BS_PV_HF%s_%s",hfFilterPlus.Data(),varN.Data())));
    if (!hMB_evtSel) {cout << "[ERROR] No histo " << Form("hMB_BS_PV_HF%s_%s",hfFilterPlus.Data(),varN.Data()) << " found" << endl; return;}
    hMB_evtSel_corrPlus = static_cast<TH1*>(hMB_evtSel->Clone(Form("hMB_BS_PV_HF%s_%s_corr",hfFilterPlus.Data(),varN.Data())));
    
    int nBins = hMB_evtSel_corr->GetNbinsX();
    for (int i = 1 ; i <= nBins ; i++)
    {
      double effVal = eff->GetEfficiency(i);
      double val = hMB_evtSel_corr->GetBinContent(i);
      
      if (effVal>0.) hMB_evtSel_corr->SetBinContent(i,val/effVal);
      else {cout << "[WARNING] efficiency < 0 for bin " << i << endl;}
      
      double effValPlus = effPlus->GetEfficiency(i);
      double valPlus = hMB_evtSel_corrPlus->GetBinContent(i);
      
      if (effValPlus>0.) hMB_evtSel_corrPlus->SetBinContent(i,valPlus/effValPlus);
      else {cout << "[WARNING] efficiency < 0 for bin " << i << endl;}
    }
    
    double emCont = 100. - (hMB_evtSel_corrPlus->Integral(1,nBins)*100./hMB_evtSel_corr->Integral(1,nBins));
    
    gStyle -> SetOptStat(0);
    TCanvas* c =  new TCanvas(Form("c_%s",varN.Data()),"", 500,500);
    
    hMB_evtSel_corr->GetXaxis()->SetRangeUser(0,100);
    hMB_evtSel_corr->SetLineColor(1);
    hMB_evtSel_corr->SetMarkerColor(1);
    hMB_evtSel_corr->SetMarkerStyle(20);
    hMB_evtSel_corr->Draw();
    
    hMB_evtSel_corrPlus->GetXaxis()->SetRangeUser(0,100);
    hMB_evtSel_corrPlus->SetLineColor(2);
    hMB_evtSel_corrPlus->SetMarkerColor(2);
    hMB_evtSel_corrPlus->SetMarkerStyle(20);
    hMB_evtSel_corrPlus->Draw("same");
    
    TLegend* l1 = new TLegend(0.40,0.64,1.00,0.76);
    legStyle(l1);
    l1->AddEntry(hMB_evtSel_corr, Form("MB+BS+PV+HF%s / eff",shfFilter.Data()), "p");
    l1->AddEntry(hMB_evtSel_corrPlus, Form("MB+BS+PV+HF%s / eff",hfFilterPlus.Data()), "p");
    
    l1->Draw("same");
    drawText(colLabel,0.54,0.2+0.65);
    drawText(Form("eff from %s",sampleLabel[1]),0.40,0.15);
    drawText(Form("EM contamination in HF%s = %0.3f %%",shfFilter.Data(),emCont),0.40,0.35);

    c->SaveAs(Form("%s/emContaminationData_eff%s_%s_HF%s.pdf",saveDIR.c_str(),sampleLabel[1],varN.Data(),shfFilter.Data()));
    aSave->Add(c);
    
  }

  // save results
  TFile *fSave = new TFile(Form("%s/emContaminationData_eff%s_HF%s.root",saveDIR.c_str(),sampleLabel[1],hfFilter),"RECREATE");
  aSave->Write("emContaminationData", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close();
}

TEfficiency* getMCeff(const char* varname,
                      const char* mcResults,
                      const char* hfFilter)
{
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string resultsDIR = mainDIR + "/Results";
  
  // Read file with eff distributions
  TFile *f = TFile::Open(Form("%s/filterEff_%s_HF%s.root",resultsDIR.c_str(),mcResults,hfFilter));
  if (!f) {
    cout << "[INFO] No file filterEff_" << mcResults << "_HF" << hfFilter << ".root was found" << endl;
    cout << "[INFO] Computing efficiencies for HF" << hfFilter << " selection" << endl;
    if (!computeMCeff(mcResults,hfFilter)) {cout << "[ERROR] Efficiency could not be computed" << endl; return 0x0;}
    f = TFile::Open(Form("%s/filterEff_%s_HF%s.root",resultsDIR.c_str(),mcResults,hfFilter));
  }
  
  TObjArray* arrHistos = static_cast<TObjArray*>(f->FindObjectAny("filterEff"));
  if (!arrHistos) {cout << "[ERROR] No efficiency array found in file " << f->GetName() << endl; return 0x0;}
  
  TEfficiency* eff = static_cast<TEfficiency*>(arrHistos->FindObject(Form("effMB_BS_PV_HF%s_%s",hfFilter,varname))->Clone());
  f->Close();
  if (!eff) {cout << "[ERROR] No efficiency " << Form("effMB_BS_PV_HF%s_%s",hfFilter,varname) << " was found" << endl; return 0x0;}
  else {return eff;}
}

bool computeMCeff(const char* mcResults,
                  const char* hfFilter)
{
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string resultsDIR = mainDIR + "/Results";
  
  // Read file with mb distributions
  TFile *f = TFile::Open(Form("%s/mbDistributions_%s.root",resultsDIR.c_str(),mcResults));
  if (!f) {cout << "[ERROR] No file mbDistributions_" << mcResults << ".root was found" << endl; return false;}
  
  // Get array with histos
  TObjArray* arrHistos = static_cast<TObjArray*>(f->FindObjectAny("mbDistributions"));
  if (!arrHistos) {cout << "[ERROR] No histos array found in file " << f->GetName() << endl; return false;}
  //
  
  TObjArray* aSave = new TObjArray();
  
  TH1* hMB(0x0);
  TH1* hMB_evtSel(0x0);
  
  const int nvar = 6;
  const char* varname[nvar] = {"hiHF", "hiHFplus", "hiHFminus","hiNtracks","hiHFhit","hiNpix"};
  
  for (int i = 0 ; i < nvar ; i++)
  {
    TString varN(varname[i]);
    TObjArray* arrVar = static_cast<TObjArray*>(arrHistos->FindObject(Form("histos_%s",varN.Data())));
    if (!arrVar) {cout << "[ERROR] No histos array found for var " << varN.Data() << " found in file mbDistributions_" << mcResults << ".root" << endl; return false;}
    
    hMB = static_cast<TH1*>(arrVar->FindObject(Form("hMB_%s",varN.Data()))->Clone());
    hMB_evtSel = static_cast<TH1*>(arrVar->FindObject(Form("hMB_BS_PV_HF%s_%s",hfFilter,varN.Data()))->Clone());
    
    if (!hMB || !hMB_evtSel) {cout << "[ERROR] No histos found to compute efficiency " << endl; return false;}
    
    TEfficiency* eff = new TEfficiency(*hMB_evtSel,*hMB);
    eff->SetName(Form("effMB_BS_PV_HF%s_%s",hfFilter,varN.Data()));
    aSave->Add(eff);
  }
  
  // save results
  TFile *fSave = new TFile(Form("%s/filterEff_%s_HF%s.root",resultsDIR.c_str(),mcResults,hfFilter),"RECREATE");
  aSave->Write("filterEff", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close();
  
  return true;
}
