class RhoTuple;

#include "TString.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/setPandaStyle.C"
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/PandaSmartLabel.C"
#include "../common_jenny.cpp"


void massresolution_breitwigner(TString inputFile=""){


	TFile * input = new TFile(inputFile, "READ");
	TTree * mc =(TTree*) input->Get("ntp");
//	TCanvas * c = (TCanvas*) input->Get("Canvas_1_n2");
//	TH1D *  h = (TH1D*) c->GetPrimitive("htemp");
	
	TCanvas * c2 = new TCanvas("c2", "c", 0,0,1000,700);
	TH1D * h = new TH1D("h", "Mass distribution for #Xi(1280)^{-}; m [GeV/c^{2}]; counts", 5000, 1.6, 2);
	mc->Project("h", "m", "N==1");

	TF1 * voigt = jenny::GetBreitWignerFit(h, 1.65, 1.93);




	h->Draw();
	voigt->Draw("SAME");




	
}
