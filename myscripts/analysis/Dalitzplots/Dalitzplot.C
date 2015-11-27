/**
* @file Dalitzplot.C
* @mainpage Dalittplot.C Analysis macro to create Dalitzplots from simple EvtGen
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


void Dalitzplot(){


	//*** Data input
	TString inputFile = "/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/evtOutput_Partwave10.root";
	TFile * data = new TFile(inputFile, "READ");

	//*** out path
	TString outPath = "/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/";
	TFile * out = new TFile(outPath+"Dalitzplots_Xi1820_simpleEvtGen_Partwave10.root", "RECREATE");

	//** get Data from File
	TTree * sim = (TTree*) data->Get("ntp");
	int nevents = sim->GetEntriesFast();

	TCanvas * c = new TCanvas("c", "Dalitz plot for Xi(1820)-", 0,0,800,500);
//	c->Divide(2,2);

	TH2D * dalitz_kpkpi = new TH2D("dalitz_kpkpi", "Dalitz plot for #Xi(1820)^{-}; m^{2}_{12}/GeV^{2}/c^{4}; m^{2}_{13}/GeV^{2}/c^{4}",200,2.5,3.2, 200,0.35,0.8);
	TH2D * dalitz_kpppi = new TH2D("dalitz_kpppi", "Dalitz plot for #Xi(1820)^{-}; m^{2}_{12}/GeV^{2}/c^{4}; m^{2}_{23}/GeV^{2}/c^{4}",200,0,5, 200,0,5);
	TH2D * dalitz_kpippi = new TH2D("dalitz_kpippi", "Dalitz plot for #Xi(1820)^{-}; m^{2}_{13}/GeV^{2}/c^{4}; m^{2}_{23}/GeV^{2}/c^{4}",200,0,5, 200,0,5);

	gStyle->SetOptStat(0);

	TLorentzVector lk, lp, lpi;
	TLorentzVector Pkp, Pkpi, Pppi;

	for (int n=0; n<nevents; n++){

		sim->GetEntry(n);

		double Ek = sim->GetLeaf("E")->GetValue(4);
		double Ep = sim->GetLeaf("E")->GetValue(5);
		double Epi = sim->GetLeaf("E")->GetValue(6);

		double Pxk = sim->GetLeaf("px")->GetValue(4);
		double Pxp = sim->GetLeaf("px")->GetValue(5);
		double Pxpi = sim->GetLeaf("px")->GetValue(6);

		double Pyk = sim->GetLeaf("py")->GetValue(4);
		double Pyp = sim->GetLeaf("py")->GetValue(5);
		double Pypi = sim->GetLeaf("py")->GetValue(6);

		double Pzk = sim->GetLeaf("pz")->GetValue(4);
		double Pzp = sim->GetLeaf("pz")->GetValue(5);
		double Pzpi = sim->GetLeaf("pz")->GetValue(6);

		lk.SetPxPyPzE(Pxk, Pyk, Pzk, Ek);
		lp.SetPxPyPzE(Pxp, Pyp, Pzp, Ep);
		lpi.SetPxPyPzE(Pxpi, Pypi, Pzpi, Epi);

		Pkp = lk + lp;
		Pppi = lp + lpi;
		Pkpi = lk + lpi;


		dalitz_kpkpi->Fill(Pkp.M2(), Pkpi.M2());
		dalitz_kpppi->Fill(Pkp.M2(), Pppi.M2());
		dalitz_kpippi->Fill(Pkpi.M2(), Pppi.M2());
	}


//	c->cd(1);
	dalitz_kpkpi->GetZaxis()->SetRangeUser(0,16);
	dalitz_kpkpi->Draw("COLZ");

//	c->cd(2);
//	dalitz_kpppi->Draw("COLZ");
//
//	c->cd(3);
//	dalitz_kpippi->Draw("COLZ");


	//****write histograms
	out->cd();

	dalitz_kpkpi->Write();

	out->Save();

	c->Print(outPath+"Dalitzplots_Xi1820_simpleEvtGen_Partwave10.pdf");
	c->Print(outPath+"Dalitzplots_Xi1820_simpleEvtGen_Partwave10.png");


}
