// Author Yeonju Go
// last modification : 2017/03/07
// 
//Headers{{{
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1F.h"
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
#include "TLatex.h"
#include "stdio.h"
#include "../../utils.h"
#include "nTower_Header.h"
//}}}

const int rebinN = 5;
float normHistHere(TH1* hNom=0, TH1* hDen=0, TH1* hRatio=0, double cut_i=700, double cut_f=900);
void nTower_plots(double *thr, const Int_t nbin, const Int_t var = 0, bool PV = true, bool BS = true, const Int_t GF = 0, const Int_t Coin = 1, const Int_t MCver = 0, float norm_i = 150, float norm_f = 700);

void draw_nTower(bool PV = true, bool BS = true, const Int_t GF = 0, const Int_t Coin = 1, const Int_t MCver = 0){
    nTower_plots(ethr, nbin_e, 0, PV, BS, GF, Coin, MCver, 150, 700);
    nTower_plots(etthr, nbin_et, 1, PV, BS, GF, Coin, MCver, 150, 700);
}

void nTower_plots(double *thr, const Int_t nbin, const Int_t var, bool PV, bool BS, const Int_t GF, const Int_t Coin, const Int_t MCver, float norm_i, float norm_f){
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle -> SetOptStat(0);

	TString Varname;
	if(var == 0) Varname = "e";
	else if(var == 1) Varname = "et";
	TString MCname;
	if(MCver == 0) MCname = "EPOS";
	else if(var == 1) MCname = "HIJING";

//make binning{{{
    double thrArr[nbin+1];
    for(int i=0;i<nbin+1;i++){
        double binwidth = (thr[1]-thr[0])/2.;
        if(i!=nbin) thrArr[i] = thr[i]-binwidth;
        else thrArr[i] = thr[nbin-1]+binwidth;
    }
//}}}

    double ratio[nbin];
    TH1D* h1D_e = new TH1D(Form("h1D_%s_ratio",Varname.Data()),Form(";(# of HF tower) %s threshold (GeV);DATA/MC",Varname.Data()), nbin,thrArr);
    h1D_e->SetTitleOffset(1.5, "Y");
    for(int i=0;i<nbin;i++){
        cout << "START " << Varname.Data() << " threshold : " << thr[i] << endl;

//Get plots{{{
        TFile* f = new TFile(Form("histfiles/nTower_%sThr%.1f_PV%d_BS%d_GF%d_Coin%d_%s.root", Varname.Data(), thr[i], PV, BS, GF, Coin, MCname.Data()), "READ");
        TH1F* hdata = (TH1F*) f->Get("h1F_sample0");
        TH1F* hmc = (TH1F*) f->Get("h1F_sample1");
        TH1F* hratio = (TH1F*) hdata->Clone("h1F_ratio");
        ratio[i] = normHistHere(hdata, hmc, hratio, norm_i, norm_f);
        cout << i << "th ratio = " << ratio[i] << endl;
        h1D_e->SetBinContent(i+1,ratio[i]);
//}}}

        TCanvas *c = new TCanvas("c","c", 500,1000);
        c->Divide(1,2);
        c->cd(1);
//MC/data comparison{{{
        gPad->SetLogy();
        hdata->SetNdivisions(505);
        hmc->SetNdivisions(505);
        hratio->SetNdivisions(505);
        Double_t range = cleverRange(hdata,hmc);
        hdata->GetYaxis()->SetRangeUser(0.1,range);
        hmc->GetYaxis()->SetRangeUser(0.1,range);
        hdata->DrawCopy();
        hmc->DrawCopy("hist same");
        drawText(Form("tower %s > %.1f",Varname.Data(),thr[i]), 0.26,0.88);
        jumSun(norm_i, 0.1, norm_i, range);
        jumSun(norm_f, 0.1, norm_f, range);
//}}}
        c->cd(2);    
//ratio{{{
        hratio->Rebin(rebinN);
        hratio->Scale(1./rebinN); 
        hratio->GetYaxis()->SetTitle("DATA/MC");
        hratio->GetYaxis()->SetRangeUser(0,2);
        hratio->DrawCopy();
        jumSun(0, 1, 1000, 1);
        drawText(Form("DATA/MC = %.3f", ratio[i]), 0.26,0.88);
//}}}
        c->SaveAs(Form("pdf/nTower_%sThr%.1f_normRange%dto%d_PV%d_BS%d_GF%d_Coin%d_%s.pdf", Varname.Data(), thr[i], (int)norm_i, (int)norm_f, PV, BS, GF, Coin, MCname.Data()));
    }
    TCanvas* c_e = new TCanvas("c_e", "",300,300);
    h1D_e->Draw("hist");
    c_e->SaveAs(Form("pdf/%sThr_vs_ratio_distribution_PV%d_BS%d_GF%d_Coin%d_%s.pdf", Varname.Data(), PV, BS, GF, Coin, MCname.Data()));
}


float normHistHere(TH1* hNom, TH1* hDen, TH1* hRatio, double cut_i, double cut_f){
    int cutBinFrom = hDen->FindBin(cut_i);
    int cutBinTo = hDen->FindBin(cut_f);
    int lastBin = hDen->FindBin(990);

    TH1::SetDefaultSumw2();
    hDen->Scale(hNom->Integral(cutBinFrom,cutBinTo)/hDen->Integral(cutBinFrom,cutBinTo));
    hRatio->Divide(hNom,hDen);

    //cout << "total integral of DATA in the whole range: " << hData->Integral()<< endl;
    //cout << "total integral of MC in the whole range: " << hMC->Integral()<< endl;
    //cout << "full efficiency in the whole range: " <<  hData->Integral()/hMC->Integral()<< endl;

    double nomInt = hNom->Integral(0,cutBinTo)+hNom->Integral(cutBinTo+1,lastBin);
    double denInt = hDen->Integral(0,cutBinTo)+hNom->Integral(cutBinTo+1,lastBin);
    return nomInt/denInt;
} 
