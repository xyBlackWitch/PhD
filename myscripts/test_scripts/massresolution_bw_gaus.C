class RhoTuple;

#include "TString.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/setPandaStyle.C"
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/PandaSmartLabel.C"
#include "../common_jenny.cpp"

void massresolution_bw_gaus(TString inputFile=""){
	
	double param[5];
	double paramf[5];



	TFile * input = new TFile(inputFile, "READ");
	TCanvas * c = (TCanvas*) input->Get("c_h_xi_m");
	TH1D *  h = (TH1D*) c->GetPrimitive("h_xi_m");
	
	double maxbin = h->GetMaximumBin();
	double max = h->GetBinLowEdge(maxbin);
	cout << max << endl;

	double mean = h->GetMean();
	cout << mean << endl;

	TCanvas * c2 = new TCanvas("c2", "c", 0,0,1000,700);

	TString func = TString::Format("gaus");
	TString func2 = TString::Format("TMath::BreitWigner(x,%.4f,[0])", mean);
	TString func3 = TString::Format("gaus(1) + TMath::BreitWigner(x, %.4f,[0])", mean);

	TF1 * bw = new TF1("bw", func2, 1.7,1.95);
	h->Fit(bw, "Q0R");
	bw->GetParameters(&param[0]);

	TF1 * gausFit = new TF1("gausFit",func, 1.78,1.86 );
	h->Fit(gausFit,"Q0R");
	gausFit->GetParameters(&param[1]);

	TF1 * total = new TF1("total",func3 , 1.6,1.95);
	total->SetParameters(param);
	total->SetLineColor(kBlack);
	h->Fit(total, "Q0R");

	total->GetParameters(&paramf[0]);

	TF1 * totalf = new TF1("totalf",func3+"*[4]" , 1.6,1.95);
	totalf->SetParameters(paramf);

	totalf->SetParName(0, "#Gamma");
	totalf->SetParName(1, "Const.");
	totalf->SetParName(2, "Mean");
	totalf->SetParName(3, "#sigma");
	totalf->SetParName(4, "A");

	h->Fit(totalf, "0R");




	h->Draw();
//	gausFit->Draw("SAME");
//	bw->Draw("SAME");
	totalf->Draw("SAME");

	
}
