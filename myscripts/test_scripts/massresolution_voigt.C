class RhoTuple;

#include "TString.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/setPandaStyle.C"
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/PandaSmartLabel.C"
#include "../common_jenny.cpp"

void massresolution_voigt(TString inputFile=""){
	
	double param[5];



	TFile * input = new TFile(inputFile, "READ");
	TCanvas * c = (TCanvas*) input->Get("c_h_xi_m");
	TH1D *  h = (TH1D*) c->GetPrimitive("h_xi_m");
	
	double maxbin = h->GetMaximumBin();
	double max = h->GetBinCenter(maxbin);
	cout << max << endl;

	double mean = h->GetMean();
	cout << mean << endl;

	TCanvas * c2 = new TCanvas("c2", "c", 0,0,1000,700);

	TString func = TString::Format("TMath::Voigt(x-%.4f, [0], [1])", max);


	TF1 * f = new TF1("f",func, 1.6,1.93 );

	h->Fit(f,"0R");

	f->GetParameters(&param[0]);


	TF1 * voigt = new TF1("voigt", func+"*[2]", 1.6, 1.93);

	voigt->SetParameter(0, param[0]);
	voigt->SetParameter(1, param[1]);


	voigt->SetParName(0, "#sigma");
	voigt->SetParName(1, "#Gamma");
	voigt->SetParName(2, "A");

	h->Fit(voigt, "0R");


	h->Draw();
	voigt->Draw("SAME");

	
}
