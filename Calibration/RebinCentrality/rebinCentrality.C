#include "TTree.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <iostream>

using namespace std;


int getBin(float v, float* array);
void getBins(float* array);



void rebinCentrality(){
  
  TTree* t[20];
  int Ntree = 3;
  
  string dirs[] = {"hiEvtAnalyzer","hltanalysis","skimanalysis"};
  string trees[] = {"HiTree","HltTree","HltTree"};
  
  TFile* inf = TFile::Open("root://eoscms//eos/cms/store/group/phys_heavyions/dhangal/pr_forests/PAMinimumBias/v3/000/285/216/HiForest_1.root");
  
  t[0] = (TTree*)inf->Get(Form("%s/%s",dirs[0].data(),trees[0].data()));
  
  float hiHF,hiHFplus,hiHFminus,hiHFminusEta4,hiHFplusEta4;
  int hiNtracks;
  UInt_t lumi;
  
  int oldbin, newbin;
  float hf;
  float hfbins[100];
  getBins(hfbins);
  
  t[0]->SetBranchAddress("lumi",&lumi);
  t[0]->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
  t[0]->SetBranchAddress("hiBin",&oldbin);
  t[0]->SetBranchAddress("hiHF",&hiHF);
  t[0]->SetBranchAddress("hiHFplus",&hiHFplus);
  t[0]->SetBranchAddress("hiHFminus",&hiHFminus);
  t[0]->SetBranchAddress("hiNtracks",&hiNtracks);
  t[0]->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
  
  TFile* outf = new TFile("output_rebinned.root","recreate");
 
  outf->mkdir(dirs[0].data());
  outf->cd(dirs[0].data());
  TTree* tnew = new TTree(trees[0].data(),"");
  
  tnew->Branch("lumi",&lumi,"lumi/F");
  tnew->Branch("hiHF",&hiHF,"hiHF/F");
  tnew->Branch("hiHFplus",&hiHFplus,"hiHFplus/F");
  tnew->Branch("hiHFminus",&hiHFminus,"hiHFminus/F");
  tnew->Branch("hiNtracks",&hiNtracks,"hiNtracks/I");
  tnew->Branch("hiHFminusEta4",&hiHFminusEta4,"hiHFminusEta4/F");
  tnew->Branch("hiHFplusEta4",&hiHFplusEta4,"hiHFplusEta4/F");
  tnew->Branch("hiBin",&newbin,"hiBin/I");
  tnew->Branch("oldBin",&oldbin,"oldBin/I");
  
  int Nevents = t[0]->GetEntries();
  //   Nevents = 1000;
  
  for(int ie = 0; ie < Nevents; ++ie){
    t[0]->GetEntry(ie);
    
    if ( ie % 1000 == 0 ) std::cout << ie << "/" << Nevents << std::endl;
    
    newbin = getBin(hiHFplusEta4,hfbins);
    //      cout<<"HF : "<<hf<<" || old bin : "<<oldbin<<"  || new bin : "<<newbin<<"   "<<endl;
    
    tnew->Fill();
  }
  
  outf->cd("");
  outf->mkdir(dirs[1].data());
  outf->cd(dirs[1].data());
  t[1] = ((TTree*)inf->Get(Form("%s/%s",dirs[1].data(),trees[1].data())))->CloneTree();
  
  outf->cd("");
  outf->mkdir(dirs[2].data());
  outf->cd(dirs[2].data());
  t[2] = ((TTree*)inf->Get(Form("%s/%s",dirs[2].data(),trees[2].data())))->CloneTree();
  
  outf->Write();
  
}



int getBin(float v, float* array){
  int bin = 0;
  
  while(v < array[bin] && bin < 99){
    bin++;
  }
  
  return bin;
}


void getBins(float* array){
  
  //These are the 2016 boundaries (2013 ones scaled)
  
  float bins[] = {
    3.420414e+01,
    3.006925e+01,
    2.757431e+01,
    2.573436e+01,
    2.426939e+01,
    2.303442e+01,
    2.199445e+01,
    2.106447e+01,
    2.023449e+01,
    1.946451e+01,
    1.877453e+01,
    1.812455e+01,
    1.752956e+01,
    1.696958e+01,
    1.644459e+01,
    1.594460e+01,
    1.546961e+01,
    1.501962e+01,
    1.458963e+01,
    1.417465e+01,
    1.376966e+01,
    1.338967e+01,
    1.302467e+01,
    1.266468e+01,
    1.231469e+01,
    1.198470e+01,
    1.165971e+01,
    1.134472e+01,
    1.103472e+01,
    1.073973e+01,
    1.044974e+01,
    1.016475e+01,
    9.884753e+00,
    9.614759e+00,
    9.354766e+00,
    9.094772e+00,
    8.844779e+00,
    8.594785e+00,
    8.354791e+00,
    8.119797e+00,
    7.889803e+00,
    7.664808e+00,
    7.444814e+00,
    7.224819e+00,
    7.004825e+00,
    6.794830e+00,
    6.584836e+00,
    6.384840e+00,
    6.184845e+00,
    5.984850e+00,
    5.794855e+00,
    5.604860e+00,
    5.419865e+00,
    5.234869e+00,
    5.054873e+00,
    4.879878e+00,
    4.704882e+00,
    4.534887e+00,
    4.369891e+00,
    4.204895e+00,
    4.039899e+00,
    3.879903e+00,
    3.724907e+00,
    3.569911e+00,
    3.419914e+00,
    3.274918e+00,
    3.124922e+00,
    2.979925e+00,
    2.839929e+00,
    2.704932e+00,
    2.569936e+00,
    2.439939e+00,
    2.309942e+00,
    2.184945e+00,
    2.059948e+00,
    1.934952e+00,
    1.819955e+00,
    1.704957e+00,
    1.594960e+00,
    1.484963e+00,
    1.379966e+00,
    1.274968e+00,
    1.179971e+00,
    1.079973e+00,
    9.849754e-01,
    8.949776e-01,
    7.999800e-01,
    7.099823e-01,
    6.199845e-01,
    5.349866e-01,
    4.449889e-01,
    3.599910e-01,
    2.699932e-01,
    1.849954e-01,
    1.049974e-01,
    5.499863e-02,
    3.499912e-02,
    1.999950e-02,
    4.999875e-03,
    0.000000e+00};
  
  for(int i = 0; i < 100; ++i){
    array[i] = bins[i];
  }
  
}

//void getBins(float* array){
// These are the 2013boundaries

//   float bins[] = {
//3.420414e+01,
//3.006925e+01,
//2.757431e+01,
//2.573436e+01,
//2.426939e+01,
//2.303442e+01,
//2.199445e+01,
//2.106447e+01,
//2.023449e+01,
//1.946451e+01,
//1.877453e+01,
//1.812455e+01,
//1.752956e+01,
//1.696958e+01,
//1.644459e+01,
//1.594460e+01,
//1.546961e+01,
//1.501962e+01,
//1.458963e+01,
//1.417465e+01,
//1.376966e+01,
//1.338967e+01,
//1.302467e+01,
//1.266468e+01,
//1.231469e+01,
//1.198470e+01,
//1.165971e+01,
//1.134472e+01,
//1.103472e+01,
//1.073973e+01,
//1.044974e+01,
//1.016475e+01,
//9.884753e+00,
//9.614759e+00,
//9.354766e+00,
//9.094772e+00,
//8.844779e+00,
//8.594785e+00,
//8.354791e+00,
//8.119797e+00,
//7.889803e+00,
//7.664808e+00,
//7.444814e+00,
//7.224819e+00,
//7.004825e+00,
//6.794830e+00,
//6.584836e+00,
//6.384840e+00,
//6.184845e+00,
//5.984850e+00,
//5.794855e+00,
//5.604860e+00,
//5.419865e+00,
//5.234869e+00,
//5.054873e+00,
//4.879878e+00,
//4.704882e+00,
//4.534887e+00,
//4.369891e+00,
//4.204895e+00,
//4.039899e+00,
//3.879903e+00,
//3.724907e+00,
//3.569911e+00,
//3.419914e+00,
//3.274918e+00,
//3.124922e+00,
//2.979925e+00,
//2.839929e+00,
//2.704932e+00,
//2.569936e+00,
//2.439939e+00,
//2.309942e+00,
//2.184945e+00,
//2.059948e+00,
//1.934952e+00,
//1.819955e+00,
//1.704957e+00,
//1.594960e+00,
//1.484963e+00,
//1.379966e+00,
//1.274968e+00,
//1.179971e+00,
//1.079973e+00,
//9.849754e-01,
//8.949776e-01,
//7.999800e-01,
//7.099823e-01,
//6.199845e-01,
//5.349866e-01,
//4.449889e-01,
//3.599910e-01,
//2.699932e-01,
//1.849954e-01,
//1.049974e-01,
//5.499863e-02,
//3.499912e-02,
//1.999950e-02,
//4.999875e-03,
//0.000000e+00};
//
//   for(int i = 0; i < 100; ++i){
//      array[i] = bins[i];
//   }
//
//}
