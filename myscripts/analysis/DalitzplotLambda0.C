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
	pbarp = 0, lambda0 = 1, antilambda0=2,
	p = 5, pim = 6
};



void DalitzplotLambda0(){


	//*** Data input
	TString inputFile ="/private/puetz/mysimulations/analysis/cascade_1820_lambda0_K/EvtGen/lambda0/evtOuput_Partwave.root";
	TFile * data = new TFile(inputFile, "READ");

	TString outPath = "/private/puetz/mysimulations/analysis/cascade_1820_lambda0_K/EvtGen/lambda0/";
	TFile * out = new TFile(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Partwave.root", "RECREATE");


	TTree * sim = (TTree*) data->Get("ntp");
	int nevents = sim->GetEntriesFast();


	TH2D * dalitz_l0ppi = new TH2D("dalitz_l0ppi", "Dalitz plot for Partwave model; m^{2}(#Lambda^{0},#pi^{-})/GeV^{2}/c^{4}; m^{2}(#pi^{-}, p)/GeV^{2}/c^{4}", 200,0,6,200,0,6);

	gStyle->SetOptStat(0);

	TLorentzVector llambda0, lantilambda0, lp, lpim;
	TLorentzVector Plambda0pim, Ppimp;

	for (int n=0; n<nevents; n++){

		sim->GetEntry(n);

		double Elambda0 = sim->GetLeaf("E")->GetValue(lambda0);
		double Eantilambda0 = sim->GetLeaf("E")->GetValue(antilambda0);
		double Ep = sim->GetLeaf("E")->GetValue(p);
		double Epim = sim->GetLeaf("E")->GetValue(pim);

		double Pxlambda0 = sim->GetLeaf("px")->GetValue(lambda0);
		double Pxantilambda0 = sim->GetLeaf("px")->GetValue(antilambda0);
		double Pxp = sim->GetLeaf("px")->GetValue(p);
		double Pxpim = sim->GetLeaf("px")->GetValue(pim);

		double Pylambda0 = sim->GetLeaf("py")->GetValue(lambda0);
		double Pyantilambda0 = sim->GetLeaf("py")->GetValue(antilambda0);
		double Pyp = sim->GetLeaf("py")->GetValue(p);
		double Pypim = sim->GetLeaf("py")->GetValue(pim);

		double Pzlambda0 = sim->GetLeaf("pz")->GetValue(lambda0);
		double Pzantilambda0 = sim->GetLeaf("pz")->GetValue(antilambda0);
		double Pzp = sim->GetLeaf("pz")->GetValue(p);
		double Pzpim = sim->GetLeaf("pz")->GetValue(pim);

		llambda0.SetPxPyPzE(Pxlambda0, Pylambda0, Pzlambda0, Elambda0);
		lantilambda0.SetPxPyPzE(Pxantilambda0, Pyantilambda0, Pzantilambda0, Eantilambda0);
		lp.SetPxPyPzE(Pxp, Pyp, Pzp, Ep);
		lpim.SetPxPyPzE(Pxpim, Pypim, Pzpim, Epim);


		Plambda0pim = lpim + llambda0;
		Ppimp = lp + lpim;


		dalitz_l0ppi->Fill(Ppimp.M2(),Plambda0pim.M2());
	}

	out->cd();

	dalitz_l0ppi->Write();

	out->Save();


	TCanvas * c = new TCanvas("c", "Dalitz plot Partwave model", 0,0,1000,900);
	dalitz_l0ppi->GetZaxis()->SetRangeUser(0,325);
	dalitz_l0ppi->Draw("COLZ");

	//****write histograms

	c->Print(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Partwave.pdf");
	c->Print(outPath+"Dalitzplots_pbarsys_simpleEvtGen_Partwave.png");


}
