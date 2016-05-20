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
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/PandaSmartLabel.C"
#include "/home/ikp1/puetz/panda/PandaSoftware/pandaroot/trunk/source/macro/setPandaStyle.C"
#include "../../common_jenny.cpp"

enum id{
	pbarp = 0, xi=1,
	lambda0 = 2, kaon = 3
};


void DalitzplotXi1820_XiPlus_lambda_k(){


	//*** Data input
	TString inputFile = "/home/ikp1/puetz/panda/mysimulations/analysis/XiPlus_Lambda_K/output_ana.root";
	TFile * data = new TFile(inputFile, "READ");

	TString outPath = "/home/ikp1/puetz/panda/mysimulations/analysis/XiPlus_Lambda_K/plots";
	TFile * out = new TFile(outPath+"/root-files/Dalitzplot_MC.root", "RECREATE");


	TTree * sim = (TTree*) data->Get("ntpMC");
	int nevents = sim->GetEntriesFast();


	TH2D * dalitz_Xilk = new TH2D("dalitz_Xilk", "Dalitz plot for MC; m^{2}(#Lambda^{0},K^{-})/GeV^{2}/c^{4}; m^{2}(#bar{#Xi}, K^{-})/GeV^{2}/c^{4}", 150,2.5,3.8,150,3.2,4.8);
	gStyle->SetOptStat(0);

	TLorentzVector lXi, lk, lla, lXi1820;
	TLorentzVector PXiK, PlaK, PXil;

	for (int n=0; n<nevents; n++){

		sim->GetEntry(n);

		double Eaxi = sim->GetLeaf("e")->GetValue(xi);
		double Ek = sim->GetLeaf("e")->GetValue(kaon);
		double Ela = sim->GetLeaf("e")->GetValue(lambda0);

		double Pxaxi = sim->GetLeaf("px")->GetValue(xi);
		double Pxk = sim->GetLeaf("px")->GetValue(kaon);
		double Pxla = sim->GetLeaf("px")->GetValue(lambda0);

		double Pyaxi = sim->GetLeaf("py")->GetValue(xi);
		double Pyk = sim->GetLeaf("py")->GetValue(kaon);
		double Pyla = sim->GetLeaf("py")->GetValue(lambda0);

		double Pzaxi = sim->GetLeaf("pz")->GetValue(xi);
		double Pzk = sim->GetLeaf("pz")->GetValue(kaon);
		double Pzla = sim->GetLeaf("pz")->GetValue(lambda0);

		lXi.SetPxPyPzE(Pxaxi, Pyaxi, Pzaxi, Eaxi);
		lk.SetPxPyPzE(Pxk, Pyk, Pzk, Ek);
		lla.SetPxPyPzE(Pxla, Pyla, Pzla, Ela);


		PXiK = lXi + lk;
		PlaK = lla + lk;
		PXil = lXi + lla;


		dalitz_Xilk->Fill(PlaK.M2(),PXiK.M2());

	}

	setPandaStyle();

	out->cd();

	dalitz_Xilk->Write();

	out->Save();


	TCanvas * c = new TCanvas("c", "Dalitz plot PHSP model", 0,0,1500,1000);
	dalitz_Xilk->Draw("COLZ");


	PandaSmartLabel("L");

	//****write histograms
	c->Print(outPath+"/png-files/Dalitzplots_MC.png");
	c->Print(outPath+"/pdf-files/Dalitzplots_MC.pdf");


}
