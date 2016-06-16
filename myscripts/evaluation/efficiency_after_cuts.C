/**
* @file efficiency_after_cuts.C
* @mainpage efficiency_after_cuts.C Method to evaluate the efficiencies after each cut
* @author Jennifer Puetz (jennifer.puetz@gmx.de)
* @date 2016
* @brief evaluate effieciens after each cut
* @details This file holds a Method to evaluate the efficiencies after each cut for a given particle.
* It is intended to be used from ROOT macros.
*
*/

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"

void efficiency_after_cuts(TString oldtrunk="", TString newtrunk="", TString ntuple="ntpLambda0", TString value="Lambda0", bool hittag=true){

	TString histvalue= TString::Format("%s_tht", value.Data());

	if (hittag){
		TString cut = TString::Format("McTruthMatch & %s_HitTag==1", value.Data());
		TString Vtxcut = TString::Format("McTruthMatch & %s_HitTag==1 & VtxFit_HowGood==1", value.Data());
		TString masscut = TString::Format("McTruthMatch & %s_HitTag==1 & VtxFit_HowGood==1 & MassFit_prob>0.01", value.Data());
	}
	else{
		TString cut = TString::Format("McTruthMatch", value.Data());
		TString Vtxcut = TString::Format("McTruthMatch & VtxFit_HowGood==1", value.Data());
		TString masscut = TString::Format("McTruthMatch & VtxFit_HowGood==1 & MassFit_prob>0.01", value.Data());
	}


	if(oldtrunk=="" || newtrunk==""){
		cout << "Enter two files for the comparision!" << endl;
	}

	//*** load ntuple for particle from files
	TFile * dataOld = new TFile(oldtrunk, "READ");
	TFile * dataNew = new TFile(newtrunk, "READ");

	TTree * oldTree = (TTree*) dataOld->Get(ntuple);
	TTree * newTree = (TTree*) dataNew->Get(ntuple);

	float nevts = 1.5e6;

	//*** Create Canvas
	TCanvas * c = new TCanvas("c", "", 0,0,800,500);
	c->SetGridy();

	gStyle->SetOptStat(0);
	//*** Create Histogram
	TH1D * h_old = new TH1D("h_old", "Efficiency after each cut ; ; efficieny [%]", 3,0,3);
	h_old->SetLineColor(kBlue);
	h_old->SetMarkerStyle(21);
	h_old->SetMarkerColor(kBlue);
	h_old->GetYaxis()->SetRangeUser(0,100);
	h_old->GetYaxis()->SetTitleSize(0.05);
	h_old->GetYaxis()->SetLabelSize(0.05);
	h_old->GetXaxis()->SetLabelSize(0.06);

	TH1D * h_new = new TH1D("h_new", "Efficiency after each cut; ; efficieny [%]", 3,0,3);
	h_new->SetLineColor(kRed);
	h_new->SetMarkerStyle(21);
	h_new->SetMarkerColor(kRed);
	h_new->GetYaxis()->SetRangeUser(0,100);
	h_new->GetYaxis()->SetTitleSize(0.05);
	h_new->GetYaxis()->SetLabelSize(0.05);
	h_new->GetXaxis()->SetLabelSize(0.06);


	//*** Get Efficiencies
	//* Mass window
	TH1D * h_masswindow = new TH1D("h_masswindow", "", 100,0,10);

	oldTree->Project("h_masswindow", histvalue, cut);
	float old_masswindow = h_masswindow->GetEntries();
	h_old->Fill("massWindow", old_masswindow/nevts*100);


	newTree->Project("h_masswindow", histvalue, cut);
	float new_masswindow = h_masswindow->GetEntries();
	h_new->Fill("massWindow", new_masswindow/nevts*100);


	//* VtxFit probability cut
	TH1D * h_vtxcut = new TH1D("h_vtxcut", "", 100,0,10);

	oldTree->Project("h_vtxcut", histvalue, Vtxcut);
	float old_vtxcut = h_vtxcut->GetEntries();
	h_old->Fill("VtxFit", old_vtxcut/nevts*100);


	newTree->Project("h_vtxcut", histvalue, Vtxcut);
	float new_vtxcut = h_vtxcut->GetEntries();
	h_new->Fill("VtxFit", new_vtxcut/nevts*100);


	//* MassFit probability cut
	TH1D * h_masscut = new TH1D("h_masscut", "", 100,0,10);

	oldTree->Project("h_masscut", histvalue, masscut);
	float old_masscut = h_masscut->GetEntries();
	h_old->Fill("MassFit", old_masscut/nevts*100);


	newTree->Project("h_masscut", histvalue, masscut);
	float new_masscut = h_masscut->GetEntries();
	h_new->Fill("MassFit", new_masscut/nevts*100);


	//*** Create a legend for the lines
  	TLegend * legend = new TLegend(0.63,0.69,0.88,0.86, "");
  	legend->AddEntry(h_old, "old trunk", "lp");
  	legend->AddEntry(h_new, "new trunk", "lp");


  	//*** Draw everything to canvas
	h_old->Draw("PL");
	h_new->Draw("SAME PL");
	legend->Draw();

	c->Print(value+"_efficiencies.pdf");


}
