
#include "TCanvas.h"
#include "../common_jenny.cpp"


void mass_diff_4C_cut(TString infile=""){

	if(infile==""){
		cout << "Infile is missing!" << endl;
		exit(0);
	}

	TFile * data = new TFile(infile, "READ");

	TTree * xisys = (TTree*) data->Get("ntpXiSys");

	int nevents = xisys->GetEntries();

	TLorentzVector lambda, kaon, XiPlus, Xi, lk, XiXip;


	TH1D * diff_mass = new TH1D("diff_mass", "differenz of invariant mass; m^{2}(#Xi(1820)^{-}) - m^{2}(Lambda, K^{-}); counts", 100, -1,1);
	TH1D * diff_mom_sys = new TH1D("diff_mom_sys", "differenz p of sys; p_{sys}-p(#Xi(1820)^{-},#bar{#Xi}^{+}); counts", 100,-1,1);
	TH1D * diff_mom = new TH1D("diff_mom", "differenz of momentum; p(#Xi(1820)^{-}) - p(#Lambda, K^{-}); counts", 100, -1,1);
//	TH1D * diff_momx = new TH1D("diff_momx", "differenz of momentum for x; px(#Xi(1820)^{-}) - px(#Lambda, K^{-}); counts", 100, -1,1);
//	TH1D * diff_momy = new TH1D("diff_momy", "differenz of momentum for y; py(#Xi(1820)^{-}) - py(#Lambda, K^{-}); counts", 100, -1,1);
//	TH1D * diff_momz = new TH1D("diff_momz", "differenz of momentum for z; pz(#Xi(1820)^{-}) - pz(#Lambda, K^{-}); counts", 100, -1,1);
	TH1D * diff_energy = new TH1D("diff_energy", "differenz of energy; E(#Xi(1820)^{-}) - E(#Lambda, K^{-}); counts", 100, -1,1);
	TH1D * diff_mass_sys = new TH1D("diff_mass_sys", "invariant mass; m(#Xi(1820)^{-},#bar{#Xi}^{+}); counts", 100, 2,5);


	for(int i=0; i<nevents; i++){



		xisys->GetEntry(i);

		int mct = xisys->GetLeaf("McTruthMatch")->GetValue(0);
		double prob = xisys->GetLeaf("4CFit_prob")->GetValue(0);

		if (mct==1 && prob>0.01){
			//*************XiPlus************************
			double Xi_m = xisys->GetLeaf("4cFit_d1m")->GetValue(0);
			double Xi_p = xisys->GetLeaf("4cFit_d1p")->GetValue(0);
			double Xi_px = xisys->GetLeaf("4cFit_d1px")->GetValue(0);
			double Xi_py = xisys->GetLeaf("4cFit_d1py")->GetValue(0);
			double Xi_pz = xisys->GetLeaf("4cFit_d1pz")->GetValue(0);
			double Xi_e = xisys->GetLeaf("4cFit_d1e")->GetValue(0);

			Xi.SetPxPyPzE(Xi_px, Xi_py, Xi_pz, Xi_e);

			//*************lambda*************************

			double lambda_px = xisys->GetLeaf("4cFit_d1d0px")->GetValue(0);
			double lambda_py = xisys->GetLeaf("4cFit_d1d0py")->GetValue(0);
			double lambda_pz = xisys->GetLeaf("4cFit_d1d0pz")->GetValue(0);
			double lambda_e = xisys->GetLeaf("4cFit_d1d0e")->GetValue(0);

			lambda.SetPxPyPzE(lambda_px, lambda_py, lambda_pz, lambda_e);

			//**************kaon**************************

			double kaon_px = xisys->GetLeaf("4cFit_d1d1px")->GetValue(0);
			double kaon_py = xisys->GetLeaf("4cFit_d1d1py")->GetValue(0);
			double kaon_pz = xisys->GetLeaf("4cFit_d1d1pz")->GetValue(0);
			double kaon_e = xisys->GetLeaf("4cFit_d1d1e")->GetValue(0);

			kaon.SetPxPyPzE(kaon_px, kaon_py, kaon_pz, kaon_e);

			//*************XiPlus*************************

			double xip_px = xisys->GetLeaf("4cFit_d0px")->GetValue(0);
			double xip_py = xisys->GetLeaf("4cFit_d0py")->GetValue(0);
			double xip_pz = xisys->GetLeaf("4cFit_d0pz")->GetValue(0);
			double xip_e = xisys->GetLeaf("4cFit_d0e")->GetValue(0);

			XiPlus.SetPxPyPzE(xip_px, xip_py, xip_pz, xip_e);

			double xisys_p = 4.6;

			lk = lambda + kaon;
			XiXip = XiPlus + Xi;

			double diff_m = Xi_m**2 - lk.M2();
			double diff_p_sys =xisys_p - XiXip.P();
			double diff_p = Xi_p - lk.P();
			double diff_px = Xi_px - lk.Px();
			double diff_py = Xi_py - lk.Py();
			double diff_pz = Xi_pz - lk.Pz();
			double diff_e = Xi_e - lk.E();

			diff_mass->Fill(diff_m);
			diff_mom_sys->Fill(diff_p_sys);
			diff_mom->Fill(diff_p);
//			diff_momx->Fill(diff_px);
//			diff_momy->Fill(diff_py);
//			diff_momz->Fill(diff_pz);
			diff_energy->Fill(diff_e);
			diff_mass_sys->Fill(XiXip.M2());
		}

	}

//	jenny::CreateDrawAndSaveHistogram(diff_mass, "", "", false, false);
	jenny::CreateDrawAndSaveHistogram(diff_mom_sys, "", "", false, false);
//	jenny::CreateDrawAndSaveHistogram(diff_mom, "", "", false, false);
//	jenny::CreateDrawAndSaveHistogram(diff_momx, "", "", false, false);
//	jenny::CreateDrawAndSaveHistogram(diff_momy, "", "", false, false);
//	jenny::CreateDrawAndSaveHistogram(diff_momz, "", "", false, false);
//	jenny::CreateDrawAndSaveHistogram(diff_energy, "", "", false, false);
	jenny::CreateDrawAndSaveHistogram(diff_mass_sys, "", "", false, false);

}
