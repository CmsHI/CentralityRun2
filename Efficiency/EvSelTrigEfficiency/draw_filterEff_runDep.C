#include "Riostream.h"
#include "../utils.h"
// have to be the same as filterEff_runDep.C's !! 
const int nfilter = 5;
const char* evtfilter[] = {"","pprimaryVertexFilter", "phfCoincFilter3", "pclusterCompatibilityFilter", "pcollisionEventSelection"};
const int colhere[] = {1,28,8,46,9};

void draw_filterEff_runDep(const char* trig = "HLT_HIL1MinimumBiasHF2AND"){ 
    ifstream in;
    in.open(Form("filterEff_%s.txt", trig));

    Float_t eff[nfilter];
    string run;
    Float_t entries;
    Int_t nlines = 0;

    TCanvas *c1 = new TCanvas("c1","",10,10,900,500);
    c1->SetGrid();
    TH1F *h[nfilter];
    TH1F *hgaus[nfilter];
    for(int i=0; i<nfilter; i++){
        h[i] = new TH1F(Form("h%d",i),"Evt selection efficiency vs. run",1,0,1);
        h[i]->SetStats(0);
        h[i]->SetTitle("Event selection efficiency vs. run; Run number; Event selection efficiency (%)");
        //h[i]->SetAxisRange(40,100,"Y");//for HF1AND trigger
        h[i]->SetAxisRange(85,100,"Y");
        h[i]->SetCanExtend(TH1::kAllAxes);
        h[i]->SetLineColor(colhere[i]);
        h[i]->SetLineWidth(2);
    
        hgaus[i] = new TH1F(Form("hgaus%d",i),"Evt selection efficiency",40,90,100);
        //hgaus[i] = new TH1F(Form("hgaus%d",i),"Evt selection efficiency",15,40,100);//for HF1AND trigger
        hgaus[i]->SetTitle(";Event selection effieiency (%);");
        //hgaus[i]->SetNdivisions(504);
    }

    while (1) {
        in >> run >> entries >> eff[0] >> eff[1] >> eff[2] >> eff[3] >> eff[4];
        if (!in.good()) break;
        for(int i=0; i<nfilter; i++){
            h[i]->Fill(run.data(),eff[i]);
            hgaus[i]->Fill(eff[i]);
        }
        nlines++;
    }
    in.close();
    c1->cd();

    TLegend* l1 = new TLegend(0.68, 0.18, 0.9, 0.35);
    //TLegend* l1 = new TLegend(0.68, 0.68, 0.9, 0.85);
    legStyle(l1);
    for(int i=1; i<nfilter; i++){
        h[i]->LabelsDeflate();
        l1->AddEntry(h[i], evtfilter[i]);
        if(i==1) h[i]->Draw("hist");
        else  h[i]->Draw("hist same");
    }
    l1->Draw("same");
    drawText(trig, 0.25, 0.18);
    c1->SaveAs(Form("pdf/evtfilterEff_runDep_%s.png",trig));

    double gausMean[nfilter], gausResol[nfilter];
    double histMean[nfilter], histResol[nfilter];
    TF1* fgaus[nfilter];
    TCanvas *c2[nfilter];
    
    for(int i=1; i<nfilter; i++){
        c2[i] = new TCanvas(Form("c2_%d",i),"",500,500);
        fgaus[i]= cleverGaus(hgaus[i]);
        gausMean[i] = fgaus[i]->GetParameter(1);
        gausResol[i] = fgaus[i]->GetParameter(2);
        histMean[i] = hgaus[i]->GetMean();
        histResol[i]= hgaus[i]->GetRMS();// GetRMS()=GetStdDev() only in ROOT!
        hgaus[i]->Draw("hist");
        fgaus[i]->Draw("same");

        float xpos(0.20), ypos(0.76);
        float dy=0.05;
        drawText(Form("from historgam"),xpos,ypos,kBlack,18);
        drawText(Form("mean = %.3f",histMean[i]),xpos,ypos-dy,kBlack,18);
        drawText(Form("sigma = %.3f",histResol[i]),xpos,ypos-2*dy,kBlack,18);
        ypos=0.55;
        drawText(Form("from gaus fit"),xpos,ypos,kBlack,18);
        drawText(Form("mean = %.3f",gausMean[i]),xpos,ypos-dy,kBlack,18);
        drawText(Form("sigma = %.3f",gausResol[i]),xpos,ypos-2*dy,kBlack,18);
        drawText(trig, 0.2, 0.18+dy);
        drawText(evtfilter[i], 0.2, 0.18);
        c2[i]->SaveAs(Form("pdf/evtfilterEff_runDep_gaus_%s_%s.png",trig,evtfilter[i]));

    }
}
