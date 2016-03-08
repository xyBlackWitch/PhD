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



void DalitzplotXi1820_reco(){


	//*** Data input
	TString inputFile = "/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/branching/spin_3half/output_ana.root";
	TFile * data = new TFile(inputFile, "READ");

	TString outPath = "/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/branching/spin_3half/plots";
//	TFile * out = new TFile(outPath+"/root-files/Dalitzplot_reco.root", "RECREATE");


	TTree * xisys = (TTree*) data->Get("ntpXiSys");

	int nevents = xisys->GetEntriesFast();


	TH2D * dalitz_Xilk = new TH2D("dalitz_Xilk", "Dalitz plot for reco; m^{2}(#Lambda^{0},K^{-})/GeV^{2}/c^{4}; m^{2}(#bar{#Xi}, K^{-})/GeV^{2}/c^{4}", 150,2.5,4,150,3.2,4.9);
	TH1D * mass_Xi_before_cut = new TH1D("mass_Xi_before_cut", "mass disitribution of m(#bar{#Xi^{+}}); m(#bar{#Xi})/GeV/c^{2}; counts", 150,-2,3);
	TH1D * mass_Xi_after_cut = new TH1D("mass_Xi_after_cut", "mass disitribution of m(#bar{#Xi^{+}}); m(#bar{#Xi})/GeV/c^{2}; counts", 150,-2,3);

	gStyle->SetOptStat(0);

	TLorentzVector lXi, lk, lla ;
	TLorentzVector PXiK, PlaK, PXil;

	for (int n=0; n<nevents; n++){

		xisys->GetEntry(n);

		int truth = xisys->GetLeaf("McTruthMatch")->GetValue(0);
		double prob = xisys->GetLeaf("4CFit_prob")->GetValue(0);


		if (truth==1 && prob>0.01){
			double Ek = xisys->GetLeaf("4cFit_d1d1e")->GetValue(0);
			double Ela = xisys->GetLeaf("4cFit_d1d0e")->GetValue(0);

			double Pyk = xisys->GetLeaf("4cFit_d1d1py")->GetValue(0);
			double Pyla = xisys->GetLeaf("4cFit_d1d0py")->GetValue(0);

			double Pzk = xisys->GetLeaf("4cFit_d1d1pz")->GetValue(0);
			double Pzla = xisys->GetLeaf("4cFit_d1d0pz")->GetValue(0);

			double Pxk = xisys->GetLeaf("4cFit_d1d1px")->GetValue(0);
			double Pxla = xisys->GetLeaf("4cFit_d1d0px")->GetValue(0);




			double Eaxi = xisys ->GetLeaf("4cFit_d0e")->GetValue(0);
			double Maxi = xisys ->GetLeaf("XiSys_d0m")->GetValue(0);
			double Pxaxi = xisys->GetLeaf("4cFit_d0px")->GetValue(0);
			double Pyaxi = xisys->GetLeaf("4cFit_d0py")->GetValue(0);
			double Pzaxi = xisys->GetLeaf("4cFit_d0pz")->GetValue(0);


			lXi.SetPxPyPzE(Pxaxi, Pyaxi, Pzaxi, Eaxi);
			lk.SetPxPyPzE(Pxk, Pyk, Pzk, Ek);
			lla.SetPxPyPzE(Pxla, Pyla, Pzla, Ela);



			PXiK = lXi + lk;
			PlaK = lla + lk;
			PXil = lXi + lla;

			double mx = PlaK.M2();
			double my = PXiK.M2();

			double cut1 = 3.9285+0.6345*sqrt(1-(mx-3.199)*(mx-3.199)/(0.559*0.559));
			double cut2 = 3.9285-0.6345*sqrt(1-(mx-3.199)*(mx-3.199)/(0.559*0.559));

			double Xim = lXi.M();



//			if (my>cut1 || my<cut2){
//				cout << "m(Xi+) after fit: " << Xim << endl;
//				cout << "m(Xi+) before fit: " << Maxi << endl;
//				cout << "Lambda: "<< endl;
//				lla.Print();
//				cout << "m(K-): "<< lk.M() <<  endl;

				mass_Xi_after_cut->Fill(lXi.M());
				mass_Xi_before_cut->Fill(Maxi);

				dalitz_Xilk->Fill(PlaK.M2(),PXiK.M2());
//			}
		}

	}



//	out->cd();
//
//	dalitz_Xilk->Write();
//
//	out->Save();

	setPandaStyle();
//	TCanvas * c = new TCanvas("c", "Dalitz plot PHSP model", 0,0,1500,1000);
//	TF1 * ellipse1 = new TF1("ellipse1", "3.9285+0.6345*sqrt(1-(x-3.199)*(x-3.199)/(0.559*0.559))",2.628,3.758);
//	TF1 * ellipse2 = new TF1("ellipse2", "3.9285-0.6345*sqrt(1-(x-3.199)*(x-3.199)/(0.559*0.559))",2.628,3.758);
//	dalitz_Xilk->GetZaxis()->SetRangeUser(0,40);
//	dalitz_Xilk->Draw("COLZ");

//	TCanvas * c1 = new TCanvas("c1", "c1", 0,0,1500,1000);
//	mass_Xi_before_cut->SetLineColor("kRed");
//	mass_Xi_before_cut->Draw();

	jenny::CreateDrawAndSaveNHistograms(mass_Xi_before_cut, mass_Xi_after_cut, "before fit", "after fit", "", "", false, false);

	PandaSmartLabel("L");

//	//****write histograms
//	c->Print(outPath+"/png-files/Dalitzplots_reco.png");
//	c->Print(outPath+"/pdf-files/Dalitzplots_reco.pdf");


}
