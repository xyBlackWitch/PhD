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


void DalitzplotXi1820_EvtGen(){


	//*** Data input
	TString inputFile ="/home/ikp1/puetz/panda/myscripts/simChain/EvtGen/XiMinus_1690_1500000_events.root";
	TFile * data = new TFile(inputFile, "READ");

	TString outPath = "/home/ikp1/puetz/panda/myscripts/simChain/EvtGen/";
	TFile * out = new TFile(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Xi1820_Xi1690.root", "RECREATE");


	TTree * sim = (TTree*) data->Get("ntp");
	int nevents = sim->GetEntriesFast();


	TH2D * dalitz_Xilk = new TH2D("dalitz_Xilk", "Dalitz plot; m^{2}(#Lambda^{0},K^{-})/GeV^{2}/c^{4}; m^{2}(#bar{#Xi}, K^{-})/GeV^{2}/c^{4}", 150,2.7,3.3,150,3.2,4.1);
//	TH2D * dalitz_Xilk = new TH2D("dalitz_Xilk", "Dalitz plot; m^{2}(#bar{#Xi}, #Lambda^{0})/GeV^{2}/c^{4}; m^{2}(#Lambda^{0},K^{-})/GeV^{2}/c^{4}",200,5.8,7.3, 200,3.2,3.8);
	TH1D * massXi = new TH1D("mass", "Mass distribution of #Xi(1820)^{-}; m/GeV/c^{2}; counts", 150, 1.61, 1.94);

	gStyle->SetOptStat(0);

	TLorentzVector lXi, lk, lla, lXi1820;
	TLorentzVector PXiK, PlaK, PXil;

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

		lXi.SetPxPyPzE(Pxaxi, Pyaxi, Pzaxi, Eaxi);
		lk.SetPxPyPzE(Pxk, Pyk, Pzk, Ek);
		lla.SetPxPyPzE(Pxla, Pyla, Pzla, Ela);
		lXi1820.SetPxPyPzE(PxXi1820, PyXi1820, PzXi1820, EXi1820);


		PXiK = lXi + lk;
		PlaK = lla + lk;
		PXil = lXi + lla;


		dalitz_Xilk->Fill(PlaK.M2(),PXiK.M2());
//		mass->Fill(lXi1820.M());
//		massXi->Fill(lXi.M());

	}

	out->cd();

	dalitz_Xilk->Write();

	out->Save();


//	TCanvas * c_mass = new TCanvas("c_mass", "Mass distribution Xi1690 model", 0,0,1000,900);
//	c->Divide(2,1);

//	c->cd(1);
//	mass->Draw();

//	c_mass->Print(outPath+"Massdist_pbarsys_simpleEvtGen_Xi1820_Xi1690.pdf");
//	c_mass->Print(outPath+"Massdist_pbarsys_simpleEvtGen_Xi1820_Xi1690.png");

	TCanvas * c = new TCanvas("c", "Dalitz plot", 0,0,1500,800);
//	c->cd(2);
//	dalitz_Xilk->GetZaxis()->SetRangeUser(0,40);
	dalitz_Xilk->Draw("COLZ");

	//****write histograms

	c->Print(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Xi1820_Xi1690.pdf");
	c->Print(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Xi1820_Xi1690.png");


}
