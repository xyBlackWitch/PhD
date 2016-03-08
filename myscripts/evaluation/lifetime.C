/**
* @file common.cpp
* @mainpage common.cpp Helper Functions
* @author Jennifer Puetz (jennifer.puetz@gmx.de)
* @date 2015
* @brief calculate decay vertex of particle
*~~~
*/

#include "TFile.h"
#include "TTree.h"
#include "../common_jenny.cpp"

void lifetime(TString pre="", bool save=true, bool close=false){

	TString inputFile = pre + "/output_ana.root";


	// Get input
	TFile * inFile = new TFile(inputFile, "READ");
	TTree * lambda0 = (TTree*) inFile->Get("ntpLambda0");
	TTree * antiLambda0 = (TTree*) inFile->Get("ntpAntiLambda0");

	//define cuts

	TString cut = "&& VtxFit_HowGood==1 && MassFit_HowGood>0";
	TString Vtxcut = "&& VtxFit_HowGood==1";


	TH1D * h_lambda0_vtx = new TH1D("h_lambda0_vtx", "c#tau for #Lambda^{0}; c#tau/cm; counts", 200, 0,200);
	lambda0->Project("h_lambda0_vtx", "VtxFit_ctau", "McTruthMatch && Lambda0_HitTag " + cut);
	jenny::CreateDrawAndSaveHistogram(h_lambda0_vtx, pre+"/plots/", "lambda0_lifetime", save, close, true);

	TH1D * h_antiLambda0_vtx = new TH1D("h_antiLambda0_vtx", "c#tau #bar{#Lambda}^{0}; c#tau/cm; counts", 200, 0,200);
	antiLambda0->Project("h_antiLambda0_vtx", "VtxFit_ctaud", "McTruthMatch && antiLambda0_HitTag " + cut);
	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_vtx, pre+"/plots/", "antiLambda0_lifetime", save, close, true);


}
