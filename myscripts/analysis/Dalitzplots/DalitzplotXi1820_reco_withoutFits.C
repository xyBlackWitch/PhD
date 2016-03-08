/**
* @file Dalitzplotxi1820.C
* @mainpage Dalitzplotxi1820.C Analysis macro to create Dalitzplots from simple EvtGen
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



void DalitzplotXi1820_reco_withoutFits(){


	//*** Data input
	TString inputFile = "/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/branching/spin_3half/output_ana.root";
	TFile * data = new TFile(inputFile, "READ");

	TString outPath = "/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/best/plots";
//	TFile * out = new TFile(outPath+"/root-files/Dalitzplot_reco.root", "RECREATE");


	TTree * xi = (TTree*) data->Get("ntpXiPlus");
	TTree * xi1820 = (TTree*) data->Get("ntpXiMinus1820");

	int nevents = xi1820->GetEntriesFast();


	TH2D * dalitz_Xilk = new TH2D("dalitz_Xilk", "Dalitz plot for reco; m^{2}(#Lambda^{0},K^{-})/GeV^{2}/c^{4}; m^{2}(#bar{#Xi}, K^{-})/GeV^{2}/c^{4}", 150,2.5,3.8,150,3.2,4.6);

	gStyle->SetOptStat(0);

	TLorentzVector lXi, lk, lla ;
	TLorentzVector PXiK, PlaK, PXil;

	for (int n=0; n<10000; n++){

		xi1820->GetEntry(n);

		double Ek = xi1820->GetLeaf("XiMinus_d1e")->GetValue(0);
		double Ela = xi1820->GetLeaf("XiMinus_d0e")->GetValue(0);

		double Pyk = xi1820->GetLeaf("XiMinus_d1py")->GetValue(0);
		double Pyla = xi1820->GetLeaf("XiMinus_d0py")->GetValue(0);

		double Pzk = xi1820->GetLeaf("XiMinus_d1pz")->GetValue(0);
		double Pzla = xi1820->GetLeaf("XiMinus_d0pz")->GetValue(0);

		double Pxk = xi1820->GetLeaf("XiMinus_d1px")->GetValue(0);
		double Pxla = xi1820->GetLeaf("XiMinus_d0px")->GetValue(0);


		xi->GetEntry(n);

		double Eaxi = xi ->GetLeaf("xiplus_e")->GetValue(0);

		double Pxaxi = xi->GetLeaf("xiplus_px")->GetValue(0);
		double Pyaxi = xi->GetLeaf("xiplus_py")->GetValue(0);
		double Pzaxi = xi->GetLeaf("xiplus_pz")->GetValue(0);


		lXi.SetPxPyPzE(Pxaxi, Pyaxi, Pzaxi, Eaxi);
		lk.SetPxPyPzE(Pxk, Pyk, Pzk, Ek);
		lla.SetPxPyPzE(Pxla, Pyla, Pzla, Ela);



		PXiK = lXi + lk;
		PlaK = lla + lk;
		PXil = lXi + lla;


		dalitz_Xilk->Fill(PlaK.M2(),PXiK.M2());

	}

//	out->cd();
//
//	dalitz_Xilk->Write();
//
//	out->Save();
//
//
	TCanvas * c = new TCanvas("c", "Dalitz plot PHSP model", 0,0,1500,1000);
//	dalitz_Xilk->GetZaxis()->SetRangeUser(0,40);
	dalitz_Xilk->Draw("COLZ");
//
//	//****write histograms
//	c->Print(outPath+"/png-files/Dalitzplots.png");


}
