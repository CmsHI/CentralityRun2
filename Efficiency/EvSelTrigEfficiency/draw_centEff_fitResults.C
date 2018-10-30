#include "Riostream.h"
#include "../utils.h"

void draw_centEff_fitResults(string cutname = "hiHF_hfCoinc3_pVtx") {

    ifstream in;
    //in.open("temp.txt");
    in.open(Form("centEff_fitResults_%s.txt",cutname.data()));

    string run;
    int runtemp;
    Float_t eff, chi2;
    Int_t nlines = 0;
    //TFile *f = new TFile("hist_effcon_runDep.root","RECREATE");
    TCanvas *c1 = new TCanvas("c1","",10,10,900,500);
    c1->SetGrid();
    TH1F *h = new TH1F("h1","Eff+Contamination vs. run",1,0,1);
    h->SetStats(0);
    h->SetTitle("Eff+Contamination vs. run; Run number;Efficiency+Contamination (%)");
    h->SetAxisRange(95,105,"Y");
    h->SetCanExtend(TH1::kAllAxes);
    h->SetLineColor(2);
    TH1F *hgaus = new TH1F("hgaus","Eff+Contamination",20,99,101);

    while (1) {
        in >> run>> eff >> chi2;
        //run = itoa(runtemp);
        if (!in.good()) break;
        if (nlines < 5) printf("run=%s, eff+contamination=%3f, reduced chi2=%3f\n",run.data(),eff,chi2);
        h->Fill(run.data(),eff);
        hgaus->Fill(eff);
        nlines++;
    }
    in.close();
    printf(" found %d points\n",nlines);
    c1->cd();
    h->LabelsDeflate();
    h->Draw("hist");
    c1->SaveAs(Form("figures/centEff_runDep_%s.png",cutname.data()));
 
    TCanvas *c2 = new TCanvas("c2","",500,500);
    double gausMean, gausResol;
    double histMean, histResol;
    TF1* fgaus = cleverGaus(hgaus);
    gausMean = fgaus->GetParameter(1);
    gausResol = fgaus->GetParameter(2);
    hgaus->SetTitle(";Efficiency+Contamination (%);");
    hgaus->SetNdivisions(504);
    hgaus->Draw("hist");
    fgaus->Draw("same");
    histMean = hgaus->GetMean();
    histResol= hgaus->GetRMS();// GetRMS()=GetStdDev() only in ROOT!
    float xpos(0.62), ypos(0.76);
    float dy=0.05;
    drawText(Form("from historgam"),xpos,ypos,kBlack,18);
    drawText(Form("mean = %.3f",histMean),xpos,ypos-dy,kBlack,18);
    drawText(Form("sigma = %.3f",histResol),xpos,ypos-2*dy,kBlack,18);
    ypos=0.55;
    drawText(Form("from gaus fit"),xpos,ypos,kBlack,18);
    drawText(Form("mean = %.3f",gausMean),xpos,ypos-dy,kBlack,18);
    drawText(Form("sigma = %.3f",gausResol),xpos,ypos-2*dy,kBlack,18);
    xpos=0.62-0.21;
    drawText("Dataset : HIMinimumBias2", xpos, 0.18+2*dy);
    drawText("Trig : HLT_HIL1MinimumBiasHF1AND_v1", xpos, 0.18+dy);
    drawText("Evt. sel. filter : hfCoinc3 && pVtx", xpos, 0.18);
    c2->SaveAs(Form("figures/centEff_runDep_gaus_%s.png",cutname.data()));

}
