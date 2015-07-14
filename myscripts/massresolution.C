class RhoTuple;

#include "TString.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"

void massresolution(){
	
	//input
	
	TString inPath = "/private/puetz/fairsoft_mar15/pandaroot/mysimulations/analysis/pbarp_lambda0_antilambda0/1000_events/idealtracking/";
	TString inFile = inPath + "output_ana.root";
	
	TFile * data = new TFile(inFile, "READ");
	
	//getting data
	
	TTree * ntpCC = (TTree*) data->Get("ntpCrossCheck"); 
	
	//Create Histogram
	
	TH1D * beforefit = new TH1D("beforefit", "Massresolution of Lambda0 before 4C-fit; m_{reco}-m_{mc}; counts", 100,-0.1,0.1);
	ntpCC->Project("beforefit", "d0m-truth_d0m","McTruthMatch==1");
	
	TH1D * afterfit = new TH1D("afterfit", "Massresolution of Lambda0 after 4C-fit; m_{reco}-m_{mc}; counts", 100,-0.1,0.1);
	ntpCC->Project("afterfit", "f4c_d0m-truth_d0m", "McTruthMatch==1");
	
	
	TCanvas * c1 = new TCanvas("c1", "c1", 0,0,800,500);
	beforefit->Draw();
	
	TCanvas * c2 = new TCanvas("c2", "c2", 0,0,800,500);
	afterfit->Draw();
}