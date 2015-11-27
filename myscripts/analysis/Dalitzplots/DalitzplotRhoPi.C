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
	pbarp = 0, rho0 = 1, pi0=2,
	pip = 3, pim = 4
};



void DalitzplotRhoPi(){


	//*** Data input
	TString inputFile ="/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/EvtGen/rho0pi/evtOutput_PHSP.root";
	TFile * data = new TFile(inputFile, "READ");

	TString outPath = "/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/EvtGen/rho0pi/";
	TFile * out = new TFile(outPath+"Dalitzplots_pbarsys_simpleEvtGen_PHSP.root", "RECREATE");


	TTree * sim = (TTree*) data->Get("ntp");
	int nevents = sim->GetEntriesFast();


	TH2D * dalitz_p0pipi = new TH2D("dalitz_p0pipi", "Dalitz plot for PHSP model; m^{2}(#pi^{-},#pi^{+})/GeV^{2}/c^{4}; m^{2}(#pi^{+}, #pi^{0})/GeV^{2}/c^{4}", 200,0,6,200,0,6);

	gStyle->SetOptStat(0);

	TLorentzVector lrho, lpi0, lpip, lpim;
	TLorentzVector Ppi0pip, Ppimpip;

	for (int n=0; n<nevents; n++){

		sim->GetEntry(n);

		double Erho = sim->GetLeaf("E")->GetValue(rho0);
		double Epi0 = sim->GetLeaf("E")->GetValue(pi0);
		double Epip = sim->GetLeaf("E")->GetValue(pip);
		double Epim = sim->GetLeaf("E")->GetValue(pim);

		double Pxrho = sim->GetLeaf("px")->GetValue(rho0);
		double Pxpi0 = sim->GetLeaf("px")->GetValue(pi0);
		double Pxpip = sim->GetLeaf("px")->GetValue(pip);
		double Pxpim = sim->GetLeaf("px")->GetValue(pim);

		double Pyrho = sim->GetLeaf("py")->GetValue(rho0);
		double Pypi0 = sim->GetLeaf("py")->GetValue(pi0);
		double Pypip = sim->GetLeaf("py")->GetValue(pip);
		double Pypim = sim->GetLeaf("py")->GetValue(pim);

		double Pzrho = sim->GetLeaf("pz")->GetValue(rho0);
		double Pzpi0 = sim->GetLeaf("pz")->GetValue(pi0);
		double Pzpip = sim->GetLeaf("pz")->GetValue(pip);
		double Pzpim = sim->GetLeaf("pz")->GetValue(pim);

		lrho.SetPxPyPzE(Pxrho, Pyrho, Pzrho, Erho);
		lpi0.SetPxPyPzE(Pxpi0, Pypi0, Pzpi0, Epi0);
		lpip.SetPxPyPzE(Pxpip, Pypip, Pzpip, Epip);
		lpim.SetPxPyPzE(Pxpim, Pypim, Pzpim, Epim);


		Ppi0pip = lpip + lpi0;
		Ppimpip = lpip + lpim;


		dalitz_p0pipi->Fill(Ppimpip.M2(),Ppi0pip.M2());
	}

	TH1D * profiley = dalitz_p0pipi->ProjectionY("profiley", 1, -1, "o");

	out->cd();

	dalitz_p0pipi->Write();

	out->Save();


	TCanvas * c = new TCanvas("c", "Dalitz plot PHSP model", 0,0,1500,1000);
	c->Divide(2,1);

	c->cd(1);
	profiley->SetTitle("Projection of y axis");
	profiley->GetYaxis()->SetTitle("counts");
	profiley->Draw();

	c->cd(2);
//	dalitz_p0pipi->GetZaxis()->SetRangeUser(0,325);
	dalitz_p0pipi->Draw("COLZ");

	//****write histograms

	c->Print(outPath+"Dalitzplots_pbarsys_simpleEvtGen_PHSP.pdf");
	c->Print(outPath+"Dalitzplots_pbarsys_simpleEvtGen_PHSP.png");


}
