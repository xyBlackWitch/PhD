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

void decay_vertex(TString pre="", bool save=true, bool close=false){

	TString inputFile = pre + "/output_ana.root";


	// Get input
	TFile * inFile = new TFile(inputFile, "READ");
	TTree * lambda0 = (TTree*) inFile->Get("ntpLambda0");
	TTree * antiLambda0 = (TTree*) inFile->Get("ntpAntiLambda0");

	//define cuts

	TString cut = "&& VtxFit_HowGood==1 && MassFit_HowGood>0";


	//Create save and Draw histograms

	TH2D * h_lambda0_vtx = new TH2D("h_lambda0_vtx", "decay vertex of #Lambda^{0}; z/cm; R/cm", 200, 0,30, 200, 0,20);
	lambda0->Project("h_lambda0_vtx", "sqrt(VtxFit_vx**2+ VtxFit_vy**2):VtxFit_vz", "McTruthMatch" + cut);
	jenny::CreateDrawAndSaveHistogram2D(h_lambda0_vtx, pre+"/plots/", "lambda0_decay_vtx", save, close, 0, 800);

	TH2D * h_antiLambda0_vtx = new TH2D("h_antiLambda0_vtx", "decay vertex of #bar{#Lambda}^{0}; z/cm; R/cm", 200, 0,30, 200, 0,20);
	antiLambda0->Project("h_antiLambda0_vtx", "sqrt(VtxFit_vx**2+ VtxFit_vy**2):VtxFit_vz", "McTruthMatch" + cut);
	jenny::CreateDrawAndSaveHistogram2D(h_antiLambda0_vtx, pre+"/plots/", "antiLambda0_decay_vtx", save, close, 0 , 800);

	TH1D * h_lambda0_vtxz = new TH1D("h_lambda0_vtxz", "decay vertex of #Lambda^{0}; c \cdot t_{prop.}; counts", 200, 0,120);
	lambda0->Project("h_lambda0_vtxz", "VtxFit_ctau", "McTruthMatch" + cut);
	jenny::CreateDrawAndSaveHistogram(h_lambda0_vtxz, pre+"/plots/", "lambda0_decay_length", save, close, true);


}
