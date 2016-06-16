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

void efficiency_after_cuts_Full_Tree(TString oldtrunk="", TString newtrunk="", TString ntuple="ntpXiSys", TString value="XiSys"){

	TString histvalue= TString::Format("%s_tht", value.Data());
	TString cut = TString::Format("McTruthMatch ", value.Data());
	TString cut4C = TString::Format("McTruthMatch & 4CFit_prob>0.01", value.Data());



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
	TH1D * h_old = new TH1D("h_old", "Efficiency after each cut ; ; efficieny [%]", 2,0,2);
	h_old->SetLineColor(kBlue);
	h_old->SetMarkerStyle(21);
	h_old->SetMarkerColor(kBlue);
	h_old->GetYaxis()->SetRangeUser(0,100);
	h_old->GetYaxis()->SetTitleSize(0.05);
	h_old->GetYaxis()->SetLabelSize(0.05);
	h_old->GetXaxis()->SetLabelSize(0.06);

	TH1D * h_new = new TH1D("h_new", "Efficiency after each cut; ; efficieny [%]", 2,0,2);
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


	//* 4CFit probability cut
	TH1D * h_cut4C = new TH1D("h_cut4C", "", 100,0,10);

	oldTree->Project("h_cut4C", histvalue, cut4C);
	float old_cut4C = h_cut4C->GetEntries();
	h_old->Fill("4CFit", old_cut4C/nevts*100);


	newTree->Project("h_cut4C", histvalue, cut4C);
	float new_cut4C = h_cut4C->GetEntries();
	h_new->Fill("4CFit", new_cut4C/nevts*100);



	//*** Create a legend for the lines
  	TLegend * legend = new TLegend(0.63,0.69,0.88,0.86, "");
  	legend->AddEntry(h_old, "old trunk", "lp");
  	legend->AddEntry(h_new, "new trunk", "lp");


  	//*** Draw everything to canvas
	h_old->Draw("PL");
	h_new->Draw("SAME PL");
	legend->Draw();

	c->Print(value + "_efficiencies.pdf");
}
