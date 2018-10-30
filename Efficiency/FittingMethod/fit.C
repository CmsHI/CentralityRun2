//Headers{{{
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "../../utils.h"
//}}}

using namespace std;

static double rangeMin = 200;
static bool doUPC = 0;
static bool doLargeRangeMin = 0;

double chi2(TH1* h1, TH1* h2){
//Get chi2{{{
    double c = 0;
    int ndof = 0;

    for(int i = h1->GetXaxis()->FindBin(rangeMin); i < h1->GetNbinsX()+1; ++i){
        double y1 = h1->GetBinContent(i);
        double y2 = h2->GetBinContent(i);
        double e1 = h1->GetBinError(i);
        double e2 = h2->GetBinError(i);

        double dy = y1-y2;
        double de2 = e1*e1+e2*e2;
        if(de2 > 0) { 
            c += dy*dy/de2;
            ndof += 1;
        }
    }
    return c/ndof;
//}}}
}

TH1* scale(TString var, TTree* nt, double sf = 1, TCut cut = ""){

//set the minimum and maximum{{{
    int nbins = 100;

    double xmin = 0;
    double xmax = 8000;

    if(var == "hiET"){
        xmax = 2000;
        rangeMin = 100;
        if(doLargeRangeMin) rangeMin = 200;
    }

    if(var == "hiEB"){
        xmax = 5000;
        rangeMin = 200;
        if(doLargeRangeMin) rangeMin = 400;
    }

    if(var == "hiEE"){
        xmax = 4000;
        rangeMin = 200;
        if(doLargeRangeMin) rangeMin = 400;
    }

    if(var == "hiHF"){
        xmax = 1000;
        rangeMin = 200;
        if(doLargeRangeMin) rangeMin = 1000;
    }

    if(var == "hiHFhit"){
        xmax = 10000;
        rangeMin = 1000;
        if(doLargeRangeMin) rangeMin = 50000;
    }

    if(var == "hiNpix"){
        xmax = 50000;
        rangeMin = 1500;
        if(doLargeRangeMin) rangeMin = 3000;
    }

    if(var == "hiNtracks"){
        xmax = 500;
        rangeMin = 200;
        if(doLargeRangeMin) rangeMin = 400;
    }
	if(var == "hiHFminusEta4"){
		xmax = 100;
		rangeMin = 10;
	}
	if(var == "hiHFminus"){
		xmax = 300;
		rangeMin = 10;
	}
	if(var == "hiHFplusEta4"){
		xmax = 100;
		rangeMin = 10;
	}
	if(var == "hiHFplus"){
		xmax = 300;
		rangeMin = 10;
	}
//}}}

    TH1D* h = new TH1D("h",Form(";%s [GeV];event fraction",var.Data()),nbins,xmin,xmax);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();

    nt->Draw(Form("%f*%s>>h",sf,var.Data()),cut);

    double rangeIntegral = h->Integral(h->GetXaxis()->FindBin(rangeMin),h->GetNbinsX()+1);
    if(rangeIntegral > 0){
        h->Scale(1./rangeIntegral);
    }
    return h;
}

void fit(TString var = "hiHF", bool PV = true, bool BS = true, const Int_t GF = 0, const Int_t Coin = 1, const Int_t MCver = 0){

//Get tree{{{
	TString RDname = "root://eoscms//eos/cms/store/group/phys_heavyions/dhangal/pr_forests/PAMinimumBias/v5/000/285/530/HiForest_";
	TString MCEname = "root://eoscms//eos/cms/store/group/phys_heavyions/jmartinb/pPb2016/MC/EPOS8TeV/HiForest_pPbEPOS8TeV_RECODEBUG";
	TString MCHname = "root://eoscms//eos/cms/store/group/phys_heavyions/jmartinb/pPb2016/MC/EPOS8TeV/HiForest_pPbEPOS8TeV_RECODEBUG";//not available yet
	TString UPCname = "root://eoscms//eos/cms/store/group/phys_heavyions/jmartinb/pPb2016/MC/EPOS8TeV/HiForest_pPbEPOS8TeV_RECODEBUG";//not available yet
	TString MCname;
	TChain* t_RD_evt, *t_RD_skim, *t_RD_hlt;
	TChain* t_MC_evt, *t_MC_skim, *t_MC_hlt;
	TChain* t_UPC_evt, *t_UPC_skim, *t_UPC_hlt;
	t_RD_evt = new TChain("hiEvtAnalyzer/HiTree");
	t_RD_skim = new TChain("skimanalysis/HltTree");
	t_RD_hlt = new TChain("hltanalysis/HltTree");
	t_MC_evt = new TChain("hiEvtAnalyzer/HiTree");
	t_MC_skim = new TChain("skimanalysis/HltTree");
	t_MC_hlt = new TChain("hltanalysis/HltTree");
	t_UPC_evt = new TChain("hiEvtAnalyzer/HiTree");
	t_UPC_skim = new TChain("skimanalysis/HltTree");
	t_UPC_hlt = new TChain("hltanalysis/HltTree");
	for(Int_t i = 0; i < 5; i++)
	{
		t_RD_evt->Add(Form("%s%d.root", RDname.Data(), i));
		t_RD_skim->Add(Form("%s%d.root", RDname.Data(), i));
		t_RD_hlt->Add(Form("%s%d.root", RDname.Data(), i));
	}
	if(MCver == 0)
	{
		t_MC_evt->Add(Form("%s.root", MCEname.Data()));
		t_MC_skim->Add(Form("%s.root", MCEname.Data()));
		t_MC_hlt->Add(Form("%s.root", MCEname.Data()));
		MCname = "EPOS";
	}
	else if(MCver == 1)
	{
		t_MC_evt->Add(Form("%s.root", MCHname.Data()));
		t_MC_skim->Add(Form("%s.root", MCHname.Data()));
		t_MC_hlt->Add(Form("%s.root", MCHname.Data()));
		MCname = "HIJING";
	}
	t_RD_evt->AddFriend(t_RD_skim);
	t_RD_evt->AddFriend(t_RD_hlt);
	t_MC_evt->AddFriend(t_MC_skim);
	t_MC_evt->AddFriend(t_MC_hlt);
	if(doUPC)
	{
		t_UPC_evt->Add(Form("%s.root", UPCname.Data()));
		t_UPC_skim->Add(Form("%s.root", UPCname.Data()));
		t_UPC_hlt->Add(Form("%s.root", UPCname.Data()));
		t_UPC_evt->AddFriend(t_UPC_skim);
		t_UPC_evt->AddFriend(t_UPC_hlt);
	}
//}}}

//Cuts{{{
	TCut trigger = "(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1==1||HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1==1)";
	TCut PVfilter;
	if(PV) PVfilter = "pPAprimaryVertexFilter==1";
	else PVfilter = "";
	TCut BSfilter;
	if(BS) BSfilter = "pBeamScrapingFilter==1";
	else BSfilter = "";
	TCut Gfilter;
	if(GF == 1) Gfilter = "pVertexFilterCutG==1";
	else if(GF == 2) Gfilter = "pVertexFilterCutGplus==1";
	else Gfilter = "";
	TCut Coinfilter;
	if(Coin == 1) Coinfilter = "phfCoincFilter1==1";
	else if(Coin == 2) Coinfilter = "phfCoincFilter2==1";
	else if(Coin == 3) Coinfilter = "phfCoincFilter3==1";
	else if(Coin == 4) Coinfilter = "phfCoincFilter4==1";
	else Coinfilter = "";
	TCut Lumi = "(67<lumi&&lumi<1540)";
//}}}

    TH1::SetDefaultSumw2();
    TH1* hupc, *href;

    if(doUPC){
        //WARNING : No trigger for UPC
        hupc = scale(var.Data(),t_UPC_evt,1., PVfilter&&BSfilter&&Gfilter&&Coinfilter&&Lumi);
        hupc->SetLineColor(4);
        hupc->SetName("hupc");
    }
    href = scale(var.Data(),t_RD_evt,1.,trigger&&PVfilter&&BSfilter&&Gfilter&&Coinfilter&&Lumi);

    cout<<"DATA ENTRIES : "<<t_RD_evt->GetEntries(trigger&&PVfilter&&BSfilter&&Gfilter&&Coinfilter)<<endl;

//scale MC and get best chi2{{{
    static const int nsf = 40;
    TH1* hsf[nsf];
    double chi2_sf[nsf];
    double sf[nsf];
    int ibest = 0;
    TGraphErrors* gsf = new TGraphErrors(nsf);

    double variation = 0.2;
    double scaleMin = 0.75;

    for(int i = 0; i < nsf; ++i){  
        //h[i]->SetName(Form("h%d",i));
        sf[i] = scaleMin+(variation/nsf)*i;
        hsf[i] = scale(var.Data(), t_MC_evt, sf[i]);
        hsf[i]->SetLineColor(2);
        chi2_sf[i] = chi2(hsf[i], href);
        gsf->SetPoint(i, sf[i], chi2_sf[i]);
        if(chi2_sf[i] < chi2_sf[ibest]) ibest = i;
    }
//}}}

    double eff = href->Integral()/(href->Integral(href->GetXaxis()->FindBin(rangeMin)+1,href->GetNbinsX()+1) + hsf[ibest]->Integral(0,hsf[ibest]->GetXaxis()->FindBin(rangeMin)));

    cout << "the best scale is : " << sf[ibest] << " i = " << ibest << endl;   
    cout << "efficiency+contamination is : " << eff << endl;

//Draw eff distribution{{{
    TCanvas* c1 = new TCanvas("c1","",600,600);
    gsf->GetXaxis()->SetTitle("scale factor");
    gsf->GetYaxis()->SetTitle("reduced #chi^2 (normalized by ndf)");
    gsf->GetXaxis()->CenterTitle();
    gsf->GetYaxis()->CenterTitle();
    gsf->SetMarkerStyle(20);
    gsf->Draw("Ap");
    gsf->GetYaxis()->SetRangeUser(0.78,3.62);
    c1->SaveAs(Form("figures/figure_%s_%s_PV%d_BS%d_GF%d_Coin%d_%s.pdf", "chi2", var.Data(), PV, BS, GF, Coin, MCname.Data()));
//}}}

//Compare data and best scale MC{{{
    TCanvas* c2 = new TCanvas("c2","",600,600);
    c2->SetLogy();
    hsf[ibest]->SetMaximum(10.*hsf[ibest]->GetMaximum());
    hsf[ibest]->Draw("hist");
    href->SetMarkerStyle(20);
    href->SetMarkerSize(0.8);
    href->Draw("same");
    TLine* vertical = new TLine(rangeMin, hsf[ibest]->GetMinimum(), rangeMin, hsf[ibest]->GetMaximum());
    vertical->SetLineStyle(2);
    vertical->Draw("same");

    TLegend * leg1 = new TLegend(0.4,0.6,0.90,0.85);
    //TLegend * leg1 = new TLegend(0.4,0.6,0.95,0.9);
    leg1->SetFillStyle(1);
    leg1->SetFillColor(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.025);
    leg1->AddEntry(href,"data","p");
    leg1->AddEntry(hsf[ibest],Form("%s based fit", MCname.Data()), "l");
    leg1->AddEntry(hsf[ibest],Form("%s scaling : %f", MCname.Data(), sf[ibest]),"");
    leg1->AddEntry(hsf[ibest],Form("Eff+Contam : %f", eff),"");
    leg1->AddEntry(hsf[ibest],Form("#chi^{2}/ndf : %f",chi2_sf[ibest]),"");
    leg1->Draw();

    c2->SaveAs(Form("figures/figure_%s_%s_PV%d_BS%d_GF%d_Coin%d_%s.pdf", "fit", var.Data(), PV, BS, GF, Coin, MCname.Data()));
//}}}

//Draw data best scale difference{{{
    TCanvas* c3 = new TCanvas("c3","c3",600,600);
    c3->SetLogx();
    TH1D* hdiff = (TH1D*)href->Clone("hdiff");
    hdiff->SetTitle(Form(";%s [GeV];data - fit",var.Data()));
    hdiff->Add(hsf[ibest],-1);
    hdiff->Draw();
    TLine* horizontal = new TLine(0,0,8000,0);
    horizontal->SetLineStyle(2);
    horizontal->Draw("same");

    double supc = 0;
    if(doUPC){
        supc = hdiff->Integral(1,5)/hupc->Integral(1,5);
        hupc->Scale(supc);
        hupc->Draw("hist same");
    }

    TLegend * leg2 = new TLegend(0.4,0.6,0.95,0.9);
    leg2->SetFillStyle(1);
    leg2->SetFillColor(0);
    leg2->SetBorderSize(0);

    leg2->SetTextSize(0.025);
    leg2->AddEntry(hdiff,"data - fit","p");

    if(doUPC){
        leg2->AddEntry(hupc,"Starlight","l");
        leg2->AddEntry(hupc,Form("Vertical scaling : %f",supc),"");
    }
    leg2->Draw();

    c3->Print(Form("figures/figure_%s_%s_PV%d_BS%d_GF%d_Coin%d_%s.pdf","diff", var.Data(), PV, BS, GF, Coin, MCname.Data()));
//}}}
}
