// to draw the ''''correlation plot'''' of centrality variables
// with the different event selection filters in one canvas.
// should modify evtfilter definition part
//
// Author Yeonju Go
// 23 Nov 2015
//
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
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
#include "yjUtility.h"
#include "par.h"

//////////// modify here as you want!! ////////////////
const double dy= 0.5;

void Get2DEffPlots(TTree* t_evt=0, TString v1="hiHF", TString v2="hiNpix", int xbin=200, double xmin=0, double xmax=4500, int ybin=200, double ymin=0, double ymax=10000, TCut cut="", const char* cap="",int isPassed=1);

void filterDep_2D()
{
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(10);
    SetHistTitleStyle(0.06,0.04);
    SetyjPadStyle();

    TChain *t_evt = new TChain("hiEvtAnalyzer/HiTree");
    TChain *t_skim = new TChain("skimanalysis/HltTree");
    TChain *t_hlt = new TChain("hltanalysis/HltTree");
    for(int i=startpoint;i<nfile;i++){
        cout << "open file : " << Form("%s/%s%d.root",fname,histoname,i) << endl;
        t_evt->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_skim->AddFile(Form("%s/%s%d.root",fname,histoname,i));
        t_hlt->AddFile(Form("%s/%s%d.root",fname,histoname,i));
    }
    t_evt-> AddFriend(t_skim);
    t_evt-> AddFriend(t_hlt);

    double hiNpixMax = 2000;
    double hiZDCMax = 80000;
    double hiHFMax = 250;
    double hiHFTrucMax = 80;
    double hiHFMinusMax = 150;
    double hiHFhitMax = 3000;
    int nbin = 100;
    for(int i=0;i<3;i++){
        Get2DEffPlots(t_evt, "hiHFplus","hiHFminus",nbin,0,hiHFMax,nbin,0,hiHFMinusMax,cut,cut.GetTitle(),i);
        Get2DEffPlots(t_evt, "hiHFhitPlus","hiHFhitMinus",nbin,0,hiHFhitMax,nbin,0,hiHFhitMax,cut,cut.GetTitle(),i);
        Get2DEffPlots(t_evt, "hiHF","hiNpix",nbin,0,hiHFMax,nbin,0,hiNpixMax,cut,cut.GetTitle(),i);
        Get2DEffPlots(t_evt, "hiHFplusEta4","hiNtracks",nbin,0,hiHFTrucMax,nbin,0,hiHFMax,cut,cut.GetTitle(),i);
        //Get2DEffPlots(t_evt, "hiHF","hiZDC",nbin,0,hiHFMax,nbin,0,hiZDCMax,cut,cut.GetTitle(),cut,cut.GetTitle(),i);
        //Get2DEffPlots(t_evt, "hiNpix","hiZDC",nbin,0,hiNpixMax,nbin,0,hiZDCMax,cut,cut.GetTitle(),cut,cut.GetTitle(),i);
    }
}
void Get2DEffPlots(TTree* t_evt, TString v1, TString v2, int xbin, double xmin, double xmax, int ybin, double ymin, double ymax, TCut cut, const char* cap,int isPassed){
    TCut totcut[nfilter];
        totcut[0] = evtfilter[0];
        totcut[1] = (totcut[0] && trig);
        if(isPassed==0)
        totcut[1] = (totcut[0] && Form("!(%s)",trig));
    for(int i=0; i<nfilter;i++){
        for(int j=2;j<=i;j++){
        if(isPassed==1)
        totcut[i] =  (totcut[i] && totcut[1] && Form("%s",evtfilter[j]));
        if(isPassed==0){
        totcut[i] = (totcut[0] && trig && Form("!(%s)",evtfilter[j]));
        }
        if(isPassed==2)
        totcut[i] = (totcut[1] && Form("%s",evtfilter[j]));
        }
        cout<<totcut[i].GetTitle()<<endl;
    } 

    TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);
    c_temp->cd();
    TH2D *h2D[nfilter];
    for(int i=0; i<nfilter; i++){
        h2D[i] = new TH2D(Form("h2D_%s_%s_filter%d_%s_isPassed%d",v1.Data(),v2.Data(),i,cap,(int)isPassed), Form("%s;%s;%s",totcut[i].GetTitle(), v1.Data(), v2.Data()), xbin, xmin, xmax, ybin, ymin, ymax);
        t_evt->Draw(Form("%s:%s>>+%s",v2.Data(), v1.Data(), h2D[i]->GetName() ), totcut[i]);
      //  h2D[i]=(TH2D*)gDirectory->Get(h2D[i]->GetName());
    }
    TCanvas *c_tot = new TCanvas("c_temp", "c_temp", 300*(ceil((float)nfilter/2)), 300*2);
    TCanvas *c[nfilter];
    c_tot->Divide((int)ceil((float)nfilter/2),2);
        TLatex t;
        t.SetTextSize(0.06);
        t.SetNDC();
        float min = h2D[0]->GetMinimum()*(1.+gStyle->GetHistTopMargin());
        float max = h2D[0]->GetMaximum()*(1.+gStyle->GetHistTopMargin());
    for(int i=0;i<nfilter;i++){
        c_tot->cd(i+1);
        h2D[i]->SetTitle("");
        h2D[i]->Draw("colz");
        h2D[i]->GetZaxis()->SetRangeUser(min,max);
        if(isPassed==1)
        t.DrawLatex(0.2,0.85,evtfilteraccushort[i]);
        if(isPassed==0 && i>=2)
        t.DrawLatex(0.2,0.85,Form("%s && !(%s)",evtfiltershort[1],evtfiltershort[i]));
        if(isPassed==2 && i>=2)
        t.DrawLatex(0.2,0.85,Form("%s && %s",evtfiltershort[1],evtfiltershort[i]));
        if(isPassed==0 && i==0)
        t.DrawLatex(0.2,0.85,Form("%s",evtfiltershort[i]));
        if(isPassed==0 && i==1)
        t.DrawLatex(0.2,0.85,Form("!(%s)",evtfiltershort[i]));
        if(isPassed==2 && i<2)
        t.DrawLatex(0.2,0.85,Form("%s",evtfiltershort[i]));
        gPad->SetLogz();

        c[i] = new TCanvas(Form("c%d",i), "", 300, 300);
        c[i]->cd();
        h2D[i]->SetTitle("");
        h2D[i]->Draw("colz");
        if(isPassed==1)
        t.DrawLatex(0.2,0.85,evtfilteraccushort[i]);
        if(isPassed==0 && i>=2)
        t.DrawLatex(0.2,0.85,Form("%s && !(%s)",evtfiltershort[1],evtfiltershort[i]));
        if(isPassed==2 && i>=2)
        t.DrawLatex(0.2,0.85,Form("%s && %s",evtfiltershort[1],evtfiltershort[i]));
        if(isPassed==0 && i<2)
        t.DrawLatex(0.2,0.85,Form("%s",evtfiltershort[i]));
        if(isPassed==2 && i<2)
        t.DrawLatex(0.2,0.85,Form("%s",evtfiltershort[i]));
        gPad->SetLogz();
        //c[i]->SaveAs(Form("figure/filterDep_2D_%s_%s_%s_isPassed%d_%s.gif", v1.Data(),v2.Data(),cap,(int)isPassed, evtfilter[i]));
    }
    c_tot->SaveAs(Form("figure/filterDep_2D_%s_%s_%s_isPassed%d.gif", v1.Data(),v2.Data(),cap,(int)isPassed));
}

