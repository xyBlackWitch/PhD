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
#include "../common_jenny.cpp"



void lorentzvector(TString inputFile=""){

	if(inputFile==""){
	//*** Data input
		inputFile = "/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/branching/spin_3half/output_ana_new.root";
	}

	TFile * data = new TFile(inputFile, "READ");

	TString outPath = "/home/ikp1/puetz/panda/mysimulations/analysis/cascade_1820_lambda0_K/branching/spin_3half/plots";


	TTree * xisys = (TTree*) data->Get("ntpXiSys");

	int nevents = xisys->GetEntriesFast();
	double mom = 4.6;
	double mp=0.938272046;
//
//
//	TH2D * dalitz_Xilk = new TH2D("dalitz_Xilk", "Dalitz plot for reco; M^{2}(#Lambda^{0},K^{-})[GeV^{2}/c^{4}]; M^{2}(#bar{#Xi}, K^{-})[GeV^{2}/c^{4}]", 150,2.5,3.8,150,3.2,4.8);
//	TH1D * mass_Xi_before_cut = new TH1D("mass_Xi_before_cut", "mass disitribution of m(#bar{#Xi^{+}}); m(#bar{#Xi})/GeV/c^{2}; counts", 150,-2,3);
//	TH1D * mass_Xi_after_cut = new TH1D("mass_Xi_after_cut", "mass disitribution of m(#bar{#Xi^{+}}); m(#bar{#Xi})/GeV/c^{2}; counts", 150,-2,3);
//
//	gStyle->SetOptStat(0);

	TLorentzVector lXi, lk, lla ;
	TLorentzVector Pf, diff;
	TLorentzVector PXiK, PlaK, PXil;
	TLorentzVector ini;

	ini.SetXYZT(0,0, mom, sqrt(mp*mp+ mom*mom)+mp);
	cout << "Initial Vector: "<< ini.Print() << endl;

	for (int n=0; n<nevents; n++){
		if ((n%10000)==0) cout << "evt "<< n <<endl;

		xisys->GetEntry(n);

		int truth = xisys->GetLeaf("McTruthMatch")->GetValue(0);
		double prob = xisys->GetLeaf("4CFit_prob")->GetValue(0);

//
		if (truth==1 && prob>0.01){
//
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


			double x = Pxk + Pxla + Pxaxi;
			double y = Pyk + Pyla + Pyaxi;
			double z = Pzk + Pzla + Pzaxi - ini.Z();
			double E = Ek + Ela + Eaxi - ini.E();

			lXi.SetPxPyPzE(Pxaxi, Pyaxi, Pzaxi, Eaxi);
			lk.SetPxPyPzE(Pxk, Pyk, Pzk, Ek);
			lla.SetPxPyPzE(Pxla, Pyla, Pzla, Ela);

			double px = (x<1e-6)? 0 : x;
			double py = (y<1e-6)? 0 : y;
			double pz = (z<1e-6)? 0 : z;

			Pf = lXi + lk + lla;
			diff = ini - Pf ;

			double diffmag = diff.M();
			double diffmom = diff.P();

			if (diffmom<1e-5 & diffmag<5e-3){
				continue;
			}
			else{
				cout << "event number: " << n << endl;
				cout << "difference to initial vector: (" << px <<", " << py << ", " << pz << ", " << E << ")" << endl;
				cout << "Probability: " << prob << endl;

				PXiK = lXi + lk;
				PlaK = lla + lk;


				double mx = PlaK.M2();
				double my = PXiK.M2();

				cout << "M2(Xi+,K-)= " << my << " M2(Lambda, K-)= " << mx << endl;
			}


////
////			double cut1 = 3.9285+0.6345*sqrt(1-(mx-3.199)*(mx-3.199)/(0.559*0.559));
////			double cut2 = 3.9285-0.6345*sqrt(1-(mx-3.199)*(mx-3.199)/(0.559*0.559));
//
//			double Xim = lXi.M();
//
//
//
////			if (my>cut1 || my<cut2){
////				cout << "m(Xi+) after fit: " << Xim << endl;
////				cout << "m(Xi+) before fit: " << Maxi << endl;
////				cout << "Lambda: "<< endl;
////				lla.Print();
////				cout << "m(K-): "<< lk.M() <<  endl;
//
//				mass_Xi_after_cut->Fill(lXi.M());
//				mass_Xi_before_cut->Fill(Maxi);
//
////				dalitz_Xilk->Fill(PlaK.M2(),PXiK.M2());
////			}
		}
//
	}



//	out->cd();
//
//	dalitz_Xilk->Write();
//
//	out->Save();

//	setPandaStyle();
////	TCanvas * c = new TCanvas("c", "Dalitz plot PHSP model", 0,0,1500,1000);
////	TF1 * ellipse1 = new TF1("ellipse1", "3.9285+0.6345*sqrt(1-(x-3.199)*(x-3.199)/(0.559*0.559))",2.628,3.758);
////	TF1 * ellipse2 = new TF1("ellipse2", "3.9285-0.6345*sqrt(1-(x-3.199)*(x-3.199)/(0.559*0.559))",2.628,3.758);
////	dalitz_Xilk->GetZaxis()->SetRangeUser(0,40);
//	dalitz_Xilk->Draw("COLZ");
//
////	TCanvas * c1 = new TCanvas("c1", "c1", 0,0,1500,1000);
////	mass_Xi_before_cut->SetLineColor("kRed");
////	mass_Xi_before_cut->Draw();
//
////	jenny::CreateDrawAndSaveNHistograms(mass_Xi_before_cut, mass_Xi_after_cut, "before fit", "after fit", "", "", false, false);
//
//	PandaSmartLabel("L");

//	//****write histograms
//	c->Print(outPath+"/png-files/Dalitzplots_reco.png");
//	c->Print(outPath+"/pdf-files/Dalitzplots_reco.pdf");


}
