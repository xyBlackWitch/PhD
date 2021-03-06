/**
* @file Dalitzplotxisys.C
* @mainpage Dalitzplotxisys.C Analysis macro to create Dalitzplots from simple EvtGen
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



void DalitzplotXi1820_reco_XiPlus_lambda_k(TString inputFile=""){

	if(inputFile==""){
	//*** Data input
		inputFile = "/home/ikp1/puetz/panda/mysimulations/analysis/XiPlus_Lambda_K/output_ana.root";
	}

	TFile * data = new TFile(inputFile, "READ");

	TString outPath = "/home/ikp1/puetz/panda/mysimulations/analysis/XiPlus_Lambda_K/plots";
	TFile * out = new TFile(outPath+"/root-files/Dalitzplot_reco.root", "RECREATE");


	TTree * xisys = (TTree*) data->Get("ntpXiSys");

	int nevents = xisys->GetEntriesFast();


	TH2D * dalitz_Xilk = new TH2D("dalitz_Xilk", "Dalitz plot for reco; M^{2}(#Lambda^{0},K^{-})[GeV^{2}/c^{4}]; M^{2}(#bar{#Xi}, K^{-})[GeV^{2}/c^{4}]", 150,2.5,3.8,150,3.2,4.8);

	gStyle->SetOptStat(0);

	TLorentzVector lXi, lk, lla ;
	TLorentzVector PXiK, PlaK, PXil;
	TLorentzVector ini;
	ini.SetXYZT(0,0,4.6,sqrt(4.6*4.6+0.938272046*0.938272046)+0.938272046);


	for (int n=0; n<nevents; n++){

		xisys->GetEntry(n);

		int truth = xisys->GetLeaf("McTruthMatch")->GetValue(0);
		double prob = xisys->GetLeaf("4CFit_prob")->GetValue(0);


		if (truth==1 && prob>0.01){

			double Pyk = xisys->GetLeaf("4cFit_d2py")->GetValue(0);
			double Pyla = xisys->GetLeaf("4cFit_d1py")->GetValue(0);

			double Pzk = xisys->GetLeaf("4cFit_d2pz")->GetValue(0);
			double Pzla = xisys->GetLeaf("4cFit_d1pz")->GetValue(0);

			double Pxk = xisys->GetLeaf("4cFit_d2px")->GetValue(0);
			double Pxla = xisys->GetLeaf("4cFit_d1px")->GetValue(0);

			double Ek = sqrt(Pxk**2+Pyk**2+Pzk**2+0.493677*0.493677);//xisys->GetLeaf("4cFit_d1d1e")->GetValue(0);
			double Ela = sqrt(Pxla**2+Pyla**2+Pzla**2+1.115683*1.115683); //xisys->GetLeaf("4cFit_d1d0e")->GetValue(0);


			double Maxi = xisys ->GetLeaf("XiSys_d0m")->GetValue(0);
			double Pxaxi = xisys->GetLeaf("4cFit_d0px")->GetValue(0);
			double Pyaxi = xisys->GetLeaf("4cFit_d0py")->GetValue(0);
			double Pzaxi = xisys->GetLeaf("4cFit_d0pz")->GetValue(0);
			double Eaxi = sqrt(Pxaxi**2+Pyaxi**2+Pzaxi**2+1.32171*1.32171);//xisys ->GetLeaf("4cFit_d0e")->GetValue(0);


			lXi.SetPxPyPzE(Pxaxi, Pyaxi, Pzaxi, Eaxi);
			lk.SetPxPyPzE(Pxk, Pyk, Pzk, Ek);
			lla.SetPxPyPzE(Pxla, Pyla, Pzla, Ela);

			PXiK = lXi + lk;
			PlaK = lla + lk;

			dalitz_Xilk->Fill(PlaK.M2(),PXiK.M2());

		}

	}



	out->cd();

	dalitz_Xilk->Write();

	out->Save();

	setPandaStyle();

	TCanvas * c = new TCanvas("c", "Dalitz plot", 0,0,1500,1000);

	dalitz_Xilk->GetZaxis()->SetRangeUser(0,50);
	dalitz_Xilk->Draw("COLZ");


	PandaSmartLabel("L");

	//****write histograms
	c->Print(outPath+"/png-files/Dalitzplots_reco.png");
	c->Print(outPath+"/pdf-files/Dalitzplots_reco.pdf");


}
