#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"

using namespace std;

bool descend(float i,float j) { return (i<j); }

void fitSlices(TH2 *hCorr, TF1 *func) {
  
  int nBins = hCorr->GetNbinsX();
  
  TH1D* hMean = new TH1D(Form("%s_1",hCorr->GetName()),"",nBins,hCorr->GetXaxis()->GetXmin(),hCorr->GetXaxis()->GetXmax());
  TH1D* hSigma = new TH1D(Form("%s_2",hCorr->GetName()),"",nBins,hCorr->GetXaxis()->GetXmin(),hCorr->GetXaxis()->GetXmax());
  
  for(int i = 1; i < nBins+1; i++){
    
    int bin = nBins - i;
    TH1D* h = hCorr->ProjectionY(Form("%s_bin%d",hCorr->GetName(),bin),i,i);
    
    func->SetParameter(0,h->GetMaximum());
    func->SetParameter(1,h->GetMean());
    func->SetParameter(2,h->GetRMS());
    
    h->Fit(func,"Q");
    
    hMean->SetBinContent(i,func->GetParameter(1));
    hMean->SetBinError(i,func->GetParError(1));
    hSigma->SetBinContent(i,func->GetParameter(2));
    hSigma->SetBinError(i,func->GetParError(2));
    
  }
  
}

void makeMCCentralityTable(int nbins = 200, const string label = "HFtowersPlusTrunc", const char * tag = "CentralityTable_HFtowersPlusTrunc200_EPOS5TeV_v80x01_mc") {
  
  TH1::SetDefaultSumw2();
  
  TString inFileName = "/afs/cern.ch/user/j/jmartinb/private/Analysis/Centrality/Forest/CMSSW_8_0_22/test/ForestpPbEPOS5TeV/crab_projects/crab_HIForest_pPbEPOS502TeV_20161106/results/HiForestAODSIM_pPbEPOS5TeV.root";
  TFile *inFile = TFile::Open(inFileName);
  TTree *t = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");
  
  TString outFileName = "CentralityTable_HFtowersPlusTrunc200_EPOS5TeV_d20161105_v1.root";
  TFile * outFile = new TFile(outFileName,"recreate");
  TDirectory* dir = outFile->mkdir(tag);
  dir->cd();
  TNtuple * nt = new TNtuple("nt","","value:bin:b:npart:ncoll:nhard");
  
  const int runNum = 1;
  CentralityBins * bins = new CentralityBins(Form("run%d",runNum), tag, nbins);
  bins->table_.reserve(nbins);
  
  ofstream txtfile("output_EPOS5TeV_HFtowersPlusTrunc_d20161105.txt");
  txtfile << "Input tree: " << inFileName << endl;
  
  double binboundaries[nbins+1];
  vector<float> values;
  
  float b, npart, ncoll, nhard, hf, hfplus, hfpluseta4, hfminuseta4, hfminus, hfhit, ee, eb, zdc, zdcplus, zdcminus;
  int run, lumi, npix, npixtrks, ntrks;
  t->SetBranchAddress("run",	&run);
  t->SetBranchAddress("lumi",	&lumi);
  t->SetBranchAddress("b",	&b);
  t->SetBranchAddress("Npart",	&npart);
  t->SetBranchAddress("Ncoll",	&ncoll);
  t->SetBranchAddress("Nhard",	&nhard);
  t->SetBranchAddress("hiHF",		&hf);
  t->SetBranchAddress("hiHFplus",	&hfplus);
  t->SetBranchAddress("hiHFplusEta4",	&hfpluseta4);
  t->SetBranchAddress("hiHFminus",	&hfminus);
  t->SetBranchAddress("hiHFminusEta4",	&hfminuseta4);
  t->SetBranchAddress("hiNpix",		&npix);
  t->SetBranchAddress("hiNpixelTracks",	&npixtrks);
  t->SetBranchAddress("hiNtracks",	&ntrks);
  
  bool binB = label.compare("b") == 0;
  bool binNpart = label.compare("Npart") == 0;
  bool binNcoll = label.compare("Ncoll") == 0;
  bool binNhard = label.compare("Nhard") == 0;
  bool binHF = label.compare("HFtowers") == 0;
  bool binHFplus = label.compare("HFtowersPlus") == 0;
  bool binHFminus = label.compare("HFtowersMinus") == 0;
  bool binHFplusTrunc = label.compare("HFtowersPlusTrunc") == 0;
  bool binHFminusTrunc = label.compare("HFtowersMinusTrunc") == 0;
  bool binNpix = label.compare("PixelHits") == 0;
  bool binNpixTrks = label.compare("PixelTracks") == 0;
  bool binNtrks = label.compare("Tracks") == 0;
  
  unsigned int Nevents = t->GetEntries();
  txtfile << "Number of events = " << Nevents << endl << endl;
  
  for(unsigned int iev = 0; iev < Nevents; iev++) {
    
    if(iev%50000 == 0) cout<<"Processing event: " << iev << " / " << Nevents << endl;
    t->GetEntry(iev);
    
    float parameter = -1;
    if(binB) parameter = b;
    if(binNpart) parameter = npart;
    if(binNcoll) parameter = ncoll;
    if(binNhard) parameter = nhard;
    if(binHF) parameter = hf;
    if(binHFplus) parameter = hfplus;
    if(binHFminus) parameter = hfminus;
    if(binHFplusTrunc) parameter = hfpluseta4;
    if(binHFminusTrunc) parameter = hfminuseta4;
    if(binNpix) parameter = npix;
    if(binNpixTrks) parameter = npixtrks;
    if(binNtrks) parameter = ntrks;
    
    values.push_back(parameter);
    
  }
  
  if(binB) sort(values.begin(),values.end(),descend);
  else sort(values.begin(),values.end());
  
  txtfile << "-------------------------------------" << endl;
  txtfile << label.data() << " based cuts are: " << endl;
  txtfile << "(";
  
  int size = values.size();
  binboundaries[nbins] = values[size-1];
  
  for(int i = 0; i < nbins; i++) {
    int entry = (int)(i*(size/nbins));
    if(entry < 0 || i == 0) binboundaries[i] = 0;
    else binboundaries[i] = values[entry];
  }
  
  for(int i = 0; i < nbins; i++) {
    if(binboundaries[i] < 0) binboundaries[i] = 0;
    txtfile << binboundaries[i] << ", ";
  }
  txtfile << binboundaries[nbins] << ")" << endl << "-------------------------------------" << endl;
  
  dir->cd();
  TH2D* hNpart = new TH2D("hNpart","",nbins,binboundaries,500,0,500);
  TH2D* hNcoll = new TH2D("hNcoll","",nbins,binboundaries,2000,0,2000);
  TH2D* hNhard = new TH2D("hNhard","",nbins,binboundaries,250,0,250);
  TH2D* hb = new TH2D("hb","",nbins,binboundaries,600,0,30);
  
  for(unsigned int iev = 0; iev < Nevents; iev++) {
    
    if( iev % 50000 == 0 ) cout<<"Processing event : " << iev << endl;
    t->GetEntry(iev);
    
    float parameter = -1;
    if(binB) parameter = b;
    if(binNpart) parameter = npart;
    if(binNcoll) parameter = ncoll;
    if(binNhard) parameter = nhard;
    if(binHF) parameter = hf;
    if(binHFplus) parameter = hfplus;
    if(binHFminus) parameter = hfminus;
    if(binHFplusTrunc) parameter = hfpluseta4;
    if(binHFminusTrunc) parameter = hfminuseta4;
    if(binNpix) parameter = npix;
    if(binNpixTrks) parameter = npixtrks;
    if(binNtrks) parameter = ntrks;
    
    hNpart->Fill(parameter,npart);
    hNcoll->Fill(parameter,ncoll);
    hNhard->Fill(parameter,nhard);
    hb->Fill(parameter,b);
    
    int bin = hNpart->GetXaxis()->FindBin(parameter) - 1;
    if(bin < 0) bin = 0;
    if(bin >= nbins) bin = nbins - 1;
    nt->Fill(parameter,bin,b,npart,ncoll,nhard);
    
  }
  
  TF1* fGaus = new TF1("fb","gaus(0)",0,2);
  fitSlices(hNpart,fGaus);
  fitSlices(hNcoll,fGaus);
  fitSlices(hNhard,fGaus);
  fitSlices(hb,fGaus);
  
  TH1D* hNpartMean = (TH1D*)gDirectory->Get("hNpart_1");
  TH1D* hNpartSigma = (TH1D*)gDirectory->Get("hNpart_2");
  TH1D* hNcollMean = (TH1D*)gDirectory->Get("hNcoll_1");
  TH1D* hNcollSigma = (TH1D*)gDirectory->Get("hNcoll_2");
  TH1D* hNhardMean = (TH1D*)gDirectory->Get("hNhard_1");
  TH1D* hNhardSigma = (TH1D*)gDirectory->Get("hNhard_2");
  TH1D* hbMean = (TH1D*)gDirectory->Get("hb_1");
  TH1D* hbSigma = (TH1D*)gDirectory->Get("hb_2");
  
  txtfile<<"-------------------------------------"<<endl;
  txtfile<<"# Bin NpartMean NpartSigma NcollMean NcollSigma bMean bSigma BinEdge"<<endl;
  for(int i = 0; i < nbins; i++){
    
    int ii = nbins-i;
    bins->table_[i].n_part_mean = hNpartMean->GetBinContent(ii);
    bins->table_[i].n_part_var = hNpartSigma->GetBinContent(ii);
    bins->table_[i].n_coll_mean = hNcollMean->GetBinContent(ii);
    bins->table_[i].n_coll_var = hNcollSigma->GetBinContent(ii);
    bins->table_[i].b_mean = hbMean->GetBinContent(ii);
    bins->table_[i].b_var = hbSigma->GetBinContent(ii);
    bins->table_[i].n_hard_mean = hNhardMean->GetBinContent(ii);
    bins->table_[i].n_hard_var = hNhardSigma->GetBinContent(ii);
    bins->table_[i].bin_edge = binboundaries[ii-1];
    
    txtfile << i << " " << hNpartMean->GetBinContent(ii) << " " << hNpartSigma->GetBinContent(ii) << " " << hNcollMean->GetBinContent(ii) << " " << hNcollSigma->GetBinContent(ii) << " " << hbMean->GetBinContent(ii) << " " <<hbSigma->GetBinContent(ii) << " " << binboundaries[ii-1] << " " <<endl;
  }
  txtfile<<"-------------------------------------"<<endl;
  
  outFile->cd();
  dir->cd();
  bins->Write();
  nt->Write();  
  bins->Delete();
  outFile->Write();
  txtfile.close();
  
}
