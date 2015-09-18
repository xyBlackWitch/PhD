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
	pbarp = 0, xi1820 = 1, xi=2,
	lambda0 = 3, kaon = 4
};


void DalitzplotPbarSystem(){


	//*** Data input
	TString inputFile ="/private/puetz/mysimulations/analysis/cascade_1820_lambda0_K/evtOutput_Partwave01.root";
	TFile * data = new TFile(inputFile, "READ");

	TString outPath = "/private/puetz/mysimulations/analysis/cascade_1820_lambda0_K/";
	TFile * out = new TFile(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Xi1820_Partwave01.root", "RECREATE");


	TTree * sim = (TTree*) data->Get("ntp");
	int nevents = sim->GetEntriesFast();


	TH2D * dalitz_aXilk = new TH2D("dalitz_aXilk", "Dalitz plot for Partwave model; m^{2}(#Lambda^{0},K^{-})/GeV^{2}/c^{4}; m^{2}(#bar{#Xi}, K^{-})/GeV^{2}/c^{4}", 200,2.5,3.8,200,3.2,4.6);
//	TH2D * dalitz_aXilk = new TH2D("dalitz_aXilk", "Dalitz plot; m^{2}(#bar{#Xi}, #Lambda^{0})/GeV^{2}/c^{4}; m^{2}(#Lambda^{0},K^{-})/GeV^{2}/c^{4}",200,5.8,7.3, 200,3.2,3.8);
	TH1D * mass = new TH1D("mass", "Mass distribution of #Xi(1820)^{-}; m/GeV/c^{2}; counts", 200, 1.61, 1.94);
	TH1D * massXibar = new TH1D("massXibar", "Mass distribution for #bar{#Xi}; m/GeV/c^{2}; counts", 200, 1.2, 1.4);

	gStyle->SetOptStat(0);

	TLorentzVector laXi, lk, lla, lXi1820;
	TLorentzVector PaXiK, PlaK, PaXil;
	TLorentzVector test;

	for (int n=0; n<nevents; n++){

		sim->GetEntry(n);

		double Eaxi = sim->GetLeaf("E")->GetValue(xi);
		double Ek = sim->GetLeaf("E")->GetValue(kaon);
		double Ela = sim->GetLeaf("E")->GetValue(lambda0);
		double EXi1820 = sim->GetLeaf("E")->GetValue(xi1820);

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

		laXi.SetPxPyPzE(Pxaxi, Pyaxi, Pzaxi, Eaxi);
		lk.SetPxPyPzE(Pxk, Pyk, Pzk, Ek);
		lla.SetPxPyPzE(Pxla, Pyla, Pzla, Ela);
		lXi1820.SetPxPyPzE(PxXi1820, PyXi1820, PzXi1820, EXi1820);


		PaXiK = laXi + lk;
		PlaK = lla + lk;
		PaXil = laXi + lla;

//		cout << "laXi: " << endl;
//		laXi.Print();
//		cout<< " LXi1820: " << endl;
//		lXi1820.Print();
//
//		test = lXi1820 + laXi;
//
//		cout << "sum:" << endl;
//
//		test.Print();

		dalitz_aXilk->Fill(PlaK.M2(),PaXiK.M2());
		mass->Fill(lXi1820.M());
		massXibar->Fill(laXi.M());

	}

	out->cd();

	dalitz_aXilk->Write();

	out->Save();


	TCanvas * c_mass = new TCanvas("c_mass", "Mass distribution Partwave01 model", 0,0,1000,900);
//	c->Divide(2,1);

//	c->cd(1);
	mass->Draw();

	c_mass->Print(outPath+"Massdist_pbarsys_simpleEvtGen_Xi1820_Partwave01.pdf");
	c_mass->Print(outPath+"Massdist_pbarsys_simpleEvtGen_Xi1820_Partwave01.png");

	TCanvas * c = new TCanvas("c", "Dalitz plot Partwave01 model", 0,0,1000,900);
//	c->cd(2);
	dalitz_aXilk->Draw("COLZ");

	//****write histograms

	c->Print(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Xi1820_Partwave01.pdf");
	c->Print(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Xi1820_Partwave01.png");


}
