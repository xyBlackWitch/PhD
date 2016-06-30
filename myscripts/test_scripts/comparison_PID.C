/**@file comparison_PID.C
 * @mainpage none
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2016
 * @brief script to compare different PIDs
 * @details This script compares the reconstruction efficiency for the full decay tree for different PIDs.
 *
 */

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "../common_jenny.cpp"
#include "TObjString.h"


void comparison_PID(TString folder="", TString pre="output_ana", float nevts=10000.){


	TString ideal = TString::Format("%s_idealPID.root", pre.Data());
	TString all = TString::Format("./%s/%s_allPID.root", folder.Data(), pre.Data());
	TString loose = TString::Format("./%s/%s_loosePID.root", folder.Data(), pre.Data());
	TString tight = TString::Format("./%s/%s_tightPID.root", folder.Data(), pre.Data());
	TString verytight = TString::Format("./%s/%s_verytightPID.root", folder.Data(), pre.Data());
	TString best = TString::Format("./%s/%s_bestPID.root", folder.Data(), pre.Data());

	TFile * idealFile = new TFile(ideal, "READ");
	TFile * allFile = new TFile(all, "READ");
	TFile * looseFile = new TFile(loose, "READ");
	TFile * tightFile = new TFile(tight, "READ");
	TFile * verytightFile = new TFile(verytight, "READ");
	TFile * bestFile = new TFile(best, "READ");

	TTree * idealTree = (TTree*) idealFile->Get("ntpXiSys");
	TTree * allTree = (TTree*) allFile->Get("ntpXiSys");
	TTree * looseTree = (TTree*) looseFile->Get("ntpXiSys");
	TTree * tightTree = (TTree*) tightFile->Get("ntpXiSys");
	TTree * verytightTree = (TTree*) verytightFile->Get("ntpXiSys");
	TTree * bestTree = (TTree*) bestFile->Get("ntpXiSys");

	TString cut = "McTruthMatch & 4CFit_prob>0.01";

	TCanvas * c = new TCanvas("c", "", 0,0,800,500);
	TH1D * h_eff = new TH1D("h_eff", "#bar{#Xi}^{+}#Xi(1820)^{-}; PID selection; efficiency [%]", 6,0,6);
	jenny::comparisonStyle(h_eff);

	TH1D * h = new TH1D("h", "tht", 100, 0, 10);

	idealTree->Project("h", "XiSys_tht", cut);
	float nideal = h->GetEntries();
	h_eff->Fill("ideal", nideal/nevts*100);

	allTree->Project("h", "XiSys_tht", cut);
	float nall = h->GetEntries();
	h_eff->Fill("all", nall/nevts*100);

	looseTree->Project("h", "XiSys_tht", cut);
	float nloose = h->GetEntries();
	h_eff->Fill("loose", nloose/nevts*100);


	TH1D * h_tight = new TH1D("h_tight", "tht", 100, 0, 10);
	tightTree->Project("h_tight", "XiSys_tht", cut);
	float ntight = h_tight->GetEntries();
	h_eff->Fill("tight", ntight/nevts*100);


	TH1D * h_verytight = new TH1D("h_verytight", "tht", 100, 0, 10);
	verytightTree->Project("h_verytight", "XiSys_tht", cut);
	float nverytight = h_verytight->GetEntries();
	h_eff->Fill("verytight", nverytight/nevts*100);


	bestTree->Project("h", "XiSys_tht", cut);
	float nbest = h->GetEntries();
	h_eff->Fill("best", nbest/nevts*100);

	h_eff->Draw("PL");
	TFile * out = new TFile("comparison_PID_"+ folder + ".root", "RECREATE");
	h_eff->Write();

	c->Print("comparison_PID_"+ folder + ".pdf");

}
