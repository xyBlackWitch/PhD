class RhoTuple;

#include "TString.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/setPandaStyle.C"
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/PandaSmartLabel.C"
#include "../common_jenny.cpp"

void massresolution(){
	
	setPandaStyle();
	//input
	
	TString inPath = "/home/ikp1/puetz/panda/mysimulations/analysis/pbarp_lambda0_antilambda0/1000_events/idealtracking/";
	TString inFile = inPath + "output_ana.root";
	
	TFile * data = new TFile(inFile, "READ");
	
	//getting data
	
	TTree * ntpCC = (TTree*) data->Get("ntpCrossCheck"); 
	
	//Create Histogram
	
	TH1D * beforefit = new TH1D("beforefit", "Massresolution of Lambda0 before 4C-fit; m_{reco}-m_{mc}; counts", 100,-0.1,0.1);
	ntpCC->Project("beforefit", "d0m-truth_d0m","McTruthMatch==1");
	
	TH1D * afterfit = new TH1D("afterfit", "Massresolution of Lambda0 after 4C-fit; m_{reco}-m_{mc}; counts", 100,-0.1,0.1);
	ntpCC->Project("afterfit", "f4c_d0m-truth_d0m", "McTruthMatch==1");
	
	jenny::CreateDrawAndSaveNHistograms(beforefit, afterfit, "beforefit", "afterfit", inPath, "mass", false, false, true);
}
