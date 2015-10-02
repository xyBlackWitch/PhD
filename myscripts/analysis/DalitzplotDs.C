/**
* @file DalitzplotPbarSystem.C
* @mainpage DalitzplotPbarSystem.C Analysis macro to create Dalitzplots from simple EvtGen
*
* @author Jennifer Puetz (jennifer.puetz@fz-juelich.de)
* @date 2015
* @brief create dalitzplots
*
*/


#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

enum id{
	pbarp = 0, Ds1p = 1, Ds1m=2,
	Ds = 3, pi0 = 4
};



void DalitzplotDs(){


	//*** Data input
	TString inputFile ="/home/ikp1/puetz/panda/mysimulations/test/EvtGen/Partwave/evtOutput_Ds.root";
	TFile * data = new TFile(inputFile, "READ");

	TString outPath = "/home/ikp1/puetz/panda/mysimulations/test/EvtGen/Partwave/";
	TFile * out = new TFile(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Partwave.root", "RECREATE");


	TTree * sim = (TTree*) data->Get("ntp");
	int nevents = sim->GetEntriesFast();


	TH2D * dalitz_Ds = new TH2D("dalitz_Ds", "Dalitz plot for Partwave model; m^{2}(D_{s1}^{-}, #pi^{0})/GeV^{2}/c^{4}; m^{2}(D_{s}^{*}, #pi^{0})/GeV^{2}/c^{4}", 200,0,50,200,0,50);

	gStyle->SetOptStat(0);

	TLorentzVector lDs, lDs1m, lDs, lpi0;
	TLorentzVector PDs1mPi0, Ppi0Ds;

	for (int n=0; n<nevents; n++){

		sim->GetEntry(n);

		double Erho = sim->GetLeaf("E")->GetValue(Ds1p);
		double EDs1m = sim->GetLeaf("E")->GetValue(Ds1m);
		double EDs = sim->GetLeaf("E")->GetValue(Ds);
		double Epi0 = sim->GetLeaf("E")->GetValue(pi0);

		double Pxrho = sim->GetLeaf("px")->GetValue(Ds1p);
		double PxDs1m = sim->GetLeaf("px")->GetValue(Ds1m);
		double PxDs = sim->GetLeaf("px")->GetValue(Ds);
		double Pxpi0 = sim->GetLeaf("px")->GetValue(pi0);

		double Pyrho = sim->GetLeaf("py")->GetValue(Ds1p);
		double PyDs1m = sim->GetLeaf("py")->GetValue(Ds1m);
		double PyDs = sim->GetLeaf("py")->GetValue(Ds);
		double Pypi0 = sim->GetLeaf("py")->GetValue(pi0);

		double Pzrho = sim->GetLeaf("pz")->GetValue(Ds1p);
		double PzDs1m = sim->GetLeaf("pz")->GetValue(Ds1m);
		double PzDs = sim->GetLeaf("pz")->GetValue(Ds);
		double Pzpi0 = sim->GetLeaf("pz")->GetValue(pi0);

		lDs.SetPxPyPzE(Pxrho, Pyrho, Pzrho, Erho);
		lDs1m.SetPxPyPzE(PxDs1m, PyDs1m, PzDs1m, EDs1m);
		lDs.SetPxPyPzE(PxDs, PyDs, PzDs, EDs);
		lpi0.SetPxPyPzE(Pxpi0, Pypi0, Pzpi0, Epi0);


		PDs1mPi0 = lDs + lDs1m;
		Ppi0Ds = lDs + lpi0;


		dalitz_Ds->Fill(Ppi0Ds.M2(),PDs1mPi0.M2());
	}

	TH1D * profiley = dalitz_Ds->ProjectionY("profiley", 1, -1, "o");

//	out->cd();
//
//	dalitz_Ds->Write();
//
//	out->Save();


	TCanvas * c = new TCanvas("c", "Dalitz plot Partwave model", 0,0,1500,1000);
	c->Divide(2,1);

	c->cd(1);
	profiley->SetTitle("Projection of y axis");
	profiley->GetYaxis()->SetTitle("counts");
	profiley->Draw();

	c->cd(2);
//	dalitz_Ds->GetZaxis()->SetRangeUser(0,325);
	dalitz_Ds->Draw("COLZ");

	//****write histograms

//	c->Print(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Partwave.pdf");
//	c->Print(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Partwave.png");


}
