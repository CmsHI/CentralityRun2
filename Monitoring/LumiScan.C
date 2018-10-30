#include <TROOT.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>

void LumiScan()
{
	TChain* chain1 = new TChain("hiEvtAnalyzer/HiTree");
	for(Int_t i = 0; i < 19; i++)
	{
		chain1->Add(Form("root:://eoscms//eos/cms/store/group/phys_heavyions/dhangal/pr_forests/PAMinimumBias/v3/000/285/090/HiForest_%d.root", i));
	}

	TCanvas* c1 = new TCanvas("c1", "", 0, 0, 600, 600);
	TH1D* h1 = new TH1D();
	chain1->Draw("lumi>>h1", "");
}
