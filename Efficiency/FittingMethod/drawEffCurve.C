// Created : 13 Feb 2018
// Modified : 13 Feb 2018
// Author : Yeonju Go
#include "../yjUtility.h"

//Double_t myfit(Double_t *x, Double_t *par) {
//    Double_t xx = x[0]+par[1];
//    Int_t binNum=histSig->FindBin(xx);
//    return par[0]*(TMath::Erf(xx));
//}
Double_t effFunction(Double_t *x, Double_t *par)
{
    Double_t xx = x[0];
    Double_t fitVal;
    if (xx <= par[0]) fitVal = 0.5*(1+TMath::Erf((xx-par[1])/(par[2])));
    else fitVal = 0.5*(1+TMath::Erf((xx-par[1])/(xx*par[2]/par[0])));
    //else fitVal = 0.5*(1+TMath::Erf((xx-par[1])/(xx*par[2]/par[0])));
    //if (xx <= par[0]) fitVal = 0.5*(1+TMath::Erf((xx-par[1])*xx/(xx*par[2])));
    //else fitVal = 0.5*(1+TMath::Erf((xx-par[3])/(xx*par[4])));
    //if (xx <= par[0]) fitVal = 0.5*(1+TMath::Erf((xx-par[1])*xx/(xx*par[2])));
    //else fitVal = 0.5*(1+TMath::Erf((xx-par[1])/(xx*par[2])));
    
   // if (xx <= par[0]) return 0.5*(1+TMath::Erf((xx-par[1])*xx/(xx*par[2])));
   // else return 0.5*(1+TMath::Erf((xx-par[1])/(xx*par[2])));
    return fitVal; 
}
void drawEffCurve(int rangeMin = 300, bool doEMsub = 0){
    gStyle->SetOptStat(0);
    SetHistTitleStyle();
    //TString var = "hiNpix";
    //TString var = "hiNtracks";
    TString var = "hiHF";
    TString cap = "";
    int rangeMax = 3000;
    int nBins = 2000;
    const char* evtCut = "PV_BS_HF3_HFORtrig";
   // const char* evtCut = "PV_BS_HF4_HFORtrig";
    const char* mcCap = "noNSDselection_epos";
    int functionN = 4;
    //const char* emSubCap = "";
    //output/fitting_hiHF_PV_BS_HF4_HFORtrig_noNSDselection_epos_20180404_MCcut1_nBins2000_rangeMin300_rangeMax2500.root
    TString totCap = Form("%s_%s_%s_20180309_MCcut0_nBins%d_rangeMin%d_rangeMax%d_doEMsub%d_doMCWeight0",var.Data(),evtCut,mcCap,nBins,rangeMin,rangeMax,(int)doEMsub);
    //TString totCap = Form("%s_%s_%s_20180309_MCcut0_nBins%d_rangeMin%d_rangeMax%d_doEMsub%d_doMCWeight0",var.Data(),evtCut,mcCap,nBins,rangeMin,rangeMax,(int)doEMsub);
    //TString totCap = Form("%s_%s_%s_20180308_MCcut0_nBins%d_rangeMin%d_rangeMax%d_doEMsub%d_doMCWeight0",var.Data(),evtCut,mcCap,nBins,rangeMin,rangeMax,(int)doEMsub);
    //TString totCap = Form("%s_%s_%s_20180308_MCcut0_nBins%d_rangeMin%d_rangeMax%d_doEMsub0_doMCWeight0",var.Data(),evtCut,mcCap,nBins,rangeMin,rangeMax);
    //TString totCap = Form("%s_%s_%s_20180308_MCcut0_nBins%d_rangeMin%d_rangeMax%d",var.Data(),evtCut,mcCap,nBins,rangeMin,rangeMax); 
    //TString totCap = Form("%s_%s_%s_20180216_MCcut0_nBins%d_rangeMin%d_rangeMax%d",var.Data(),evtCut,mcCap,nBins,rangeMin,rangeMax); 
    //TString totCap = Form("%s_%s_%s_20180215_MCcut0_nBins%d_rangeMin%d_rangeMax%d",var.Data(),evtCut,mcCap,nBins,rangeMin,rangeMax); 
    //TString totCap = Form("%s_%s_%s_20180213_MCcut0_nBins%d_rangeMin%d_rangeMax%d",var.Data(),evtCut,mcCap,nBins,rangeMin,rangeMax); 
    TFile* f = new TFile(Form("output/fitting_%s.root",totCap.Data()));
    TH1D* hdata = (TH1D*) f->Get("hdata");
    TH1D* hmc = (TH1D*) f->Get("hmc");
   
    double xmax = 50; 

    //EventSelectionEfficiency
    TCanvas* c_evtSelEff = new TCanvas(Form("c_evtSelEff_%s_%s",var.Data(),cap.Data()), "",300,600);
    c_evtSelEff->Divide(1,2);
    c_evtSelEff->cd(1);
    //gPad->SetLogy();
   // hmc->Scale(1.,"width");
   // hdata->Scale(1.,"width");
    int rangeMinOld = rangeMin;
    rangeMin=50;
    double data_area = hdata->Integral("width");
    double mc_area = hmc->Integral("width");
    double eff = hdata->Integral("width")/(hdata->Integral(hdata->GetXaxis()->FindBin(rangeMin)+1,hdata->GetNbinsX()+1,"width") + hmc->Integral(0,hmc->GetXaxis()->FindBin(rangeMin),"width"));
    cout << "data total = " << hdata->Integral("width") << endl;
    cout << "mc until rangeMin = " << hmc->Integral(0,hmc->GetXaxis()->FindBin(rangeMin),"width") << endl;
    cout << "data from rangeMin = " << (hdata->Integral(hdata->GetXaxis()->FindBin(rangeMin)+1,hdata->GetNbinsX()+1,"width")) << endl;
    //hdata->Scale(data_evtSel_area);
    SetHistTextSize(hmc);
    hmc->SetMarkerSize(0.8);
    hmc->SetMarkerColor(1);
    hdata->SetMarkerColor(2);
    hmc->SetMarkerStyle(21);
    hdata->SetMarkerStyle(20);
    hmc->GetXaxis()->SetRangeUser(0,xmax);
    hdata->GetXaxis()->SetRangeUser(0,xmax);
    hmc->GetYaxis()->SetRangeUser(0,0.04);
    hdata->GetYaxis()->SetRangeUser(0,0.04);
    hmc->SetTitle(";#SigmaE_{HF} (GeV);Arbitrary Normalization");
    hmc->DrawCopy(); 
    hdata->DrawCopy("same");
    double DATA_integEvtSelEff = data_area/mc_area; 
    drawText(Form("Evt. Sel. Eff. = %.2f %s",eff*100,"%"),0.35,0.6);
   // drawText(Form("Evt. Sel. Eff. = %.2f %s",DATA_integEvtSelEff*100,"%"),0.40,0.6);
    //drawText(Form("Evt. Sel. Eff. = %.2f %s",DATA_integEvtSelEff*100,"%"),0.2,0.3,0,kBlack,26);
    
    drawText(Form("Norm. range > %d",(int)rangeMinOld),0.12,0.85);
    
    c_evtSelEff->cd(2);
    //TH1D* ratio = new TH1D("ratio","",2000,0,4000);
    TH1D* ratio = (TH1D*) hdata->Clone("ratio");
    SetHistTextSize(ratio);
    //ih->SetTitleOffset(titleoffset, "X");
    //SetHistTitleOffsetStyle(ratio);
    ratio->Divide(hdata,hmc);
    ratio->GetYaxis()->SetTitle("Eff.+Contam. (=DATA/EPOS)");
    ratio->GetXaxis()->SetTitle("#Sigma E_{HF} (GeV)");
    //ratio->GetXaxis()->SetTitleCenter();
    //if(doEMsub)ratio->GetYaxis()->SetTitle("Eff+Contamination (=DATA_EMsub/EPOS)");
    //ratio->GetXaxis()->SetRangeUser(5,20);
    ratio->GetXaxis()->SetRangeUser(0,xmax);
    ratio->GetYaxis()->SetRangeUser(0,1.2);
    ratio->Draw();
    jumSun(0,1,xmax,1);
   
    //////////////////////////////////////////////////
    // FITTING
    //
    //TF1 f_myfit = 
    // func[ieta] = new TF1(Form("func_%d", ieta), "[0]/(TMath::Exp((-x+[1])/[2]+1)+[3])+[4]", 3.3, 20); 
    TF1 *f_temp;
    if(functionN==0) f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])/(TMath::Sqrt(x*x)*[2])))", 0, xmax); //nominal
    else if(functionN==1) f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])*x/(TMath::Sqrt(x*x)*[2])))", 0, xmax);//well-matched at high efficiency region
    else if(functionN==2) f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])/(TMath::Sqrt(x)*[2])))", 0, xmax);
    else if(functionN==3) f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])*x/(TMath::Sqrt(x)*[2])))", 0, xmax);
    else if(functionN==4) f_temp = new TF1("func", effFunction, 0, xmax,3);
    //TF1 *f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])/(TMath::Sqrt(x*x)*[2])))", 0, xmax); 
    //TF1 *f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])/(TMath::Sqrt(x*x)*[2]))*([3]*x+[4]))", 0, xmax); 
    //TF1 *f_temp = new TF1("func", " [0]+[1]*(1+TMath::Erf((x-[2])/(TMath::Sqrt(x)*[3])))", 0, xmax); 
    f_temp->SetNpx(1000); 
    if(functionN==3) f_temp->SetParameters(0.50,1.34114e+01,7.5804); 
    else if(functionN==4) f_temp->SetParameters(12,1.34114e+01,7.5804,13.3955,0.2563); 
    else f_temp->SetParameters(0.5,1.34114e+01,7.45804e-01); 
    //f_temp->SetParameters(1.31443e-02,0.50,1.34114e+01,7.45804e-01); 
    if(functionN!=4) f_temp->FixParameter(0, 0.5); 
    if(functionN==4) f_temp->SetParLimits(0, 11,25); 
    //TF1 *f_temp = new TF1("func", "[0]+[1]*(1+TMath::Erf((x-[2])/TMath::Sqrt(x)))", 0, xmax); 
    //TF1 *f_temp = new TF1("func", "[0]*TMath::ATan([1]*x) +[2]", 0, xmax); 
    //TF1 *f_temp = new TF1("func", " (exp([0]-x)*TMath::Erf((x-[1])/[2]))+[3]", 0, xmax); 
    //TF1 *f_temp = new TF1("func", " (exp(-[0]*x)*TMath::Erf((x-[1])/[2]))+[3]", 0, xmax); 
    //TF1 *f_temp = new TF1("func", "([0]*TMath::Erf((x-[1])/[2]))+[3]+[4]*x*x", 0, xmax); 
    //TF1 *f_temp = new TF1("func", "(TMath::Erf([0]+[1]*x+[2]*x*x))", 0, xmax); 
    //TF1 *f_temp = new TF1("func", "([0]*TMath::Erf((x-[1]/[2])))+[3]", 0, xmax); 
    //TF1 *f_temp = new TF1("func", "([0]*TMath::Erf((x-[1])/[2]))+[3]", 0, xmax); 
    //TF1 *f_temp = new TF1("func", "[0]/(TMath::Exp((-x+[1])/[2]+1)+[3])+[4]", 0, xmax); 
    //TF1 *f_erf = new TF1("erf(x)","[0]TMath::erf(x-[1])+[2]",0,xmax);
    //TF1 *f_erf = new TF1("erf(x)","ROOT::Math::erf(x+[1])+[2]",0,xmax);
    TF1 *f_gaussian_cdf = new TF1("gaussian_cdf(x)","ROOT::Math::gaussian_cdf(x)",0,xmax);
//    f_temp->SetParameters(15, 13.18, 2.54,0.50); 
 //   f_temp->SetParLimits(2, 0, 50); 
    //f_temp->SetParameters(0.59, 3.18, 0.97,0.40); 
//    f_temp->SetParameters(-4.795, 1.653, 0.923, 3.018, 2.369); 
    //f_temp->SetParameters(1,0,1,0,10);
    ratio->Fit(f_temp);
    ratio->Fit(f_temp);
    ratio->Fit(f_temp);
    ratio->Fit(f_temp);
    ratio->Fit(f_temp);

    f_temp->Draw("same");
    cout << "x="<< xmax << ": " << f_temp->Eval(xmax) << endl;
    cout << "x="<< 100 << ": " << f_temp->Eval(100) << endl;
   // drawText(Form("bin width at %d GeV = %.1f",(int)xmax,ratio->GetBinWidth(ratio->FindBin(xmax))), 0.4,0.3+0.05);
   // drawText(Form("Max. Eff. = %.2f %%",100*f_temp->Eval(50)), 0.4,0.3);
    //f_erf->Draw("same");
    //f_gaussian_cdf->Draw("same");
    //f_erf->Draw("same");
    if(functionN==0) drawText(Form("0.5*(1+Erf( #frac{(x-%.4f)}{%.4f*x} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)),0.40,0.3); 
    else if(functionN==1) drawText(Form("0.5*(1+Erf( #frac{x*(x-%.4f)}{%.4f*x} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)),0.40,0.3); 
    else if(functionN==2)drawText(Form("0.5*(1+Erf( #frac{x-%.4f}{%.4f*#sqrt{x}} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)),0.40,0.3); 
    else if(functionN==3)drawText(Form("0.5*(1+Erf( #frac{x*(x-%.4f)}{%.4f*#sqrt{x}} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)),0.40,0.3); 
    else if(functionN==4) {
        double xpos = 0.34;
        double ypos = 0.5;
        drawText(Form("if x#leq%.2f",f_temp->GetParameter(0)),xpos,ypos); 
        drawText(Form("0.5*(1+Erf( #frac{x-%.4f}{%.4f} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)),xpos+0.02,ypos-0.06); 
        drawText(Form("if x>%.2f",f_temp->GetParameter(0)),xpos, ypos-0.15); 
        drawText(Form("0.5*(1+Erf( #frac{x-%.4f}{%.4f*x} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)/f_temp->GetParameter(0)),xpos+0.02,ypos-0.21); 
        //drawText(Form("(x#leq%.1f); 0.5*(1+Erf( #frac{x*(x-%.4f)}{%.4f*x} ))",f_temp->GetParameter(0),f_temp->GetParameter(1),f_temp->GetParameter(2)),0.35,0.3); 
        //drawText(Form("(x>%.1f); 0.5*(1+Erf( #frac{x-%.4f}{%.4f*x} ))",f_temp->GetParameter(0),f_temp->GetParameter(3),f_temp->GetParameter(4)),0.35,0.3-0.1); 
        
    }
    //drawText(Form("#frac{1}{2}(1+Erf(#frac{x-%.2f}{%.2f*#sqrt{x}}))",f_temp->GetParameter(1),f_temp->GetParameter(2)),0.35,0.3); 
    //drawText(Form("0.5*(1+Erf((x-%.2f)/(%.2f*#sqrt{x})))",f_temp->GetParameter(1),f_temp->GetParameter(2)),0.35,0.3); 

    ////////////////////////////////
    // Cross check
    int NBINS = hdata->GetNbinsX();
    TH1D* hinv = new TH1D("hinv", "",NBINS,0,4000); 
    //TH1D* hinv = (TH1D*) hdata->Clone("hinv");
    for(int i=1;i<hdata->GetNbinsX()+1;++i){
        double xVal = hdata->GetBinContent(i);
        double xErr = hdata->GetBinError(i);
        double x = hdata->GetBinCenter(i);
        double scale = 1./f_temp->Eval(x);
        if(xVal==0){
            hinv->SetBinContent(i,0);
            hinv->SetBinError(i,0);
        } else{
            hinv->SetBinContent(i,xVal*scale);
            hinv->SetBinError(i,xErr);
        }
        //if(i<20) cout << "x = " << x << ", xVal = " << xVal << ", xcale = " << scale << ", xVal*scale = " << xVal*scale << endl;
        //if(i==hdata->GetNbinsX()) cout << "x = " << x <<endl;
    }
    c_evtSelEff->cd(1);
    hinv->SetMarkerColor(4);
    hinv->SetMarkerStyle(32);
    hinv->Draw("p same"); 
    
    // LEGEND 
    TLegend* l2 = new TLegend(0.5,0.7,0.9,0.9);
    legStyle(l2);
    if(doEMsub){
        l2->AddEntry(hdata,"Data (EM subtracted)","pl");
        l2->AddEntry(hmc,"EPOS","pl");
        l2->AddEntry(hinv,"Data/Eff."); 
    } else{
        l2->AddEntry(hdata,"Data","pl");
        l2->AddEntry(hmc,"EPOS","pl");
        l2->AddEntry(hinv,"Data/Eff."); 
    } 
    l2->Draw("same");
    //cout << "inverted hinv integral = " << hinv->Integral()<< endl;

    //double effinv = hdata->Integral(0,hdata->GetNbinsX()+1,"width")/(hdata->Integral(hdata->GetXaxis()->FindBin(rangeMin)+1,hdata->GetNbinsX()+1,"width") + hmc->Integral(0,hmc->GetXaxis()->FindBin(rangeMin),"width"));
    double effinv = (hdata->Integral(0,hdata->GetNbinsX()+1))/(hdata->Integral(hdata->GetXaxis()->FindBin(rangeMin)+1,hdata->GetNbinsX()+1) + hinv->Integral(0,hinv->GetXaxis()->FindBin(rangeMin)));
    double effinv_2 = hdata->Integral(0,hdata->GetNbinsX()+1)/hinv->Integral();
    cout << hmc->Integral(0,hmc->GetXaxis()->FindBin(rangeMin)) << endl;
    cout << hinv->Integral(0,hinv->GetXaxis()->FindBin(rangeMin)) << endl;
    cout << hdata->Integral() << endl;
    cout << "nbins = " << hdata->GetNbinsX()+1 << endl;
    //cout << "data (hdata->Integral(hdata->GetXaxis()->FindBin(rangeMin)+1,hdata->GetNbinsX()+1)<<endl;
    cout << hdata->Integral(hdata->GetXaxis()->FindBin(rangeMin)+1,hdata->GetNbinsX()+1)<<endl;
    cout << "eff  = " << eff << endl;
    cout << "eff inverted = " << effinv << endl;
    cout << "eff inverted 2 = " << effinv_2 << endl;
    cout << "data total = " << hdata->Integral("width") << endl;
    cout << "mc until rangeMin = " << hmc->Integral(0,hmc->GetXaxis()->FindBin(rangeMin),"width") << endl;
    cout << "data from rangeMin = " << (hdata->Integral(hdata->GetXaxis()->FindBin(rangeMin)+1,hdata->GetNbinsX()+1,"width")) << endl;

    drawText("Evt. Sel. Eff.(from Data/Eff.)",0.35,0.6-0.05);
    drawText(Form("= %.2f %s",effinv*100,"%"),0.35+0.25,0.6-2*0.05);

    c_evtSelEff->SaveAs(Form("figures/effCurve_%s_%s_func%d.pdf",var.Data(),totCap.Data(),functionN));
    TString fout_name = Form("output/effCurve_%s_func%d.root",totCap.Data(),functionN);
    TFile* fout = new TFile(fout_name,"recreate");
    fout->cd();
    f_temp->Write();
    ratio->Write();
    hmc->Write();
    hdata->Write();
}
