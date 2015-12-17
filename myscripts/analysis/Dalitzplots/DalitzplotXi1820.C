/**
* @file DalitzplotXi1820.C
* @mainpage DalitzplotXi1820.C Analysis macro to create Dalitzplots from simple EvtGen
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
	pbarp = 0, xi1820 = 1, xi=2,
	lambda0 = 3, kaon = 4
};


void DalitzplotXi1820(){


	//*** Data input
	TString inputFile = "/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/branching/spin_3half/output_ana.root";
	TFile * data = new TFile(inputFile, "READ");

	TString outPath = "/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/branching/spin_3half/plots";
	TFile * out = new TFile(outPath+"/root-files/Dalitzplot_MC.root", "RECREATE");


	TTree * sim = (TTree*) data->Get("ntpMC");
	int nevents = sim->GetEntriesFast();


	TH2D * dalitz_Xilk = new TH2D("dalitz_Xilk", "Dalitz plot for MC; m^{2}(#Lambda^{0},K^{-})/GeV^{2}/c^{4}; m^{2}(#bar{#Xi}, K^{-})/GeV^{2}/c^{4}", 150,2.5,3.8,150,3.2,4.6);

	gStyle->SetOptStat(0);

	TLorentzVector lXi, lk, lla, lXi1820;
	TLorentzVector PXiK, PlaK, PXil;

	for (int n=0; n<nevents; n++){

		sim->GetEntry(n);

		double Eaxi = sim->GetLeaf("e")->GetValue(xi);
		double Ek = sim->GetLeaf("e")->GetValue(kaon);
		double Ela = sim->GetLeaf("e")->GetValue(lambda0);
		double EXi1820 = sim->GetLeaf("e")->GetValue(xi1820);

		double Pxaxi = sim->GetLeaf("px")->GetValue(xi);
		double Pxk = sim->GetLeaf("px")->GetValue(kaon);
		double Pxla = sim->GetLeaf("px")->GetValue(lambda0);
		double PxXi1820 = sim->GetLeaf("px")->GetValue(xi1820);

		double Pyaxi = sim->GetLeaf("py")->GetValue(xi);
		double Pyk = sim->GetLeaf("py")->GetValue(kaon);
		double Pyla = sim->GetLeaf("py")->GetValue(lambda0);
		double PyXi1820 = sim->GetLeaf("py")->GetValue(xi1820);

		double Pzaxi = sim->GetLeaf("pz")->GetValue(xi);
		double Pzk = sim->GetLeaf("pz")->GetValue(kaon);
		double Pzla = sim->GetLeaf("pz")->GetValue(lambda0);
		double PzXi1820 = sim->GetLeaf("pz")->GetValue(xi1820);

		lXi.SetPxPyPzE(Pxaxi, Pyaxi, Pzaxi, Eaxi);
		lk.SetPxPyPzE(Pxk, Pyk, Pzk, Ek);
		lla.SetPxPyPzE(Pxla, Pyla, Pzla, Ela);
		lXi1820.SetPxPyPzE(PxXi1820, PyXi1820, PzXi1820, EXi1820);


		PXiK = lXi + lk;
		PlaK = lla + lk;
		PXil = lXi + lla;


		dalitz_Xilk->Fill(PlaK.M2(),PXiK.M2());

	}

	out->cd();

	dalitz_Xilk->Write();

	out->Save();

//
	TCanvas * c = new TCanvas("c", "Dalitz plot PHSP model", 0,0,1500,1000);
//	dalitz_Xilk->GetZaxis()->SetRangeUser(0,40);
	dalitz_Xilk->Draw("COLZ");

	//****write histograms
	c->Print(outPath+"/png-files/Dalitzplots_MC.png");


}
