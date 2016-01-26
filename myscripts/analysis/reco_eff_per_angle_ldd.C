/**
 * @file reco_eff_per_angle_ldd.C
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Creates a plot with the reco efficiency as a function of the angle.
 * @details Methode Creates a plot with the reco efficiency as a function of the angle.. data_ldd is comming from the analysis which is stored in a file called output_ana.root
 * Basically ROOT.
 */

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "../common_jenny.cpp"

reco_eff_per_angle_ldd(TString pre1 ="", TString pre2=""){

	TFile * data_ldd = new TFile(pre1, "READ");
	TFile * data = new TFile(pre2, "READ");

	TTree * ntpMC_ldd = (TTree*) data_ldd->Get("ntpMC");
	TTree * ntpAntiProton_ldd = (TTree*) data_ldd->Get("ntpAntiProton");

	TTree * ntpMC = (TTree*) data->Get("ntpMC");
	TTree * ntpAntiProton = (TTree*) data->Get("ntpAntiProton");

	TString tht_cut_aprot;
	TString tht_cut_mc;

	TH1D * h_mc_ldd = new TH1D("h_mc_ldd", "h_mc_ldd", 25, 0,25);
	TH1D * h_aprot_ldd = new TH1D("h_aprot_ldd", "h_aprot_ldd", 25,0,25);

	TH1D * h_mc = new TH1D("h_mc", "h_mc", 25, 0,25);
	TH1D * h_aprot = new TH1D("h_aprot", "h_aprot", 25,0,25);

	TH1D * h_eff_ldd = new TH1D("h_eff_ldd", " reco efficiency ; #Theta/Deg; efficiency",25,0,25);
	TH1D * h_eff = new TH1D("h_eff", "reco efficiency; #Theta/Deg; efficiency",25,0,25);

	for(double theta = 0; theta<25; theta++){


		double tht_min =theta/180*TMath::Pi();
		double tht_max = (theta+1)/180*TMath::Pi();

		tht_cut_aprot.Form("AntiProton_MC_tht>= %f && AntiProton_MC_tht< %f", tht_min, tht_max);
		tht_cut_mc.Form("tht>=%f && tht< %f", tht_min, tht_max);

		ntpAntiProton_ldd->Project("h_aprot_ldd", "AntiProton_MC_tht*180/pi", "McTruthMatch && Mother==-3122 && AntiProton_HitTag && " + tht_cut_aprot);
		ntpAntiProton->Project("h_aprot", "AntiProton_MC_tht*180/pi", "McTruthMatch && Mother==-3122 && AntiProton_HitTag && " + tht_cut_aprot);

		double nevts_reco_ldd = h_aprot_ldd->GetEntries();
		double nevts_reco = h_aprot->GetEntries();

		ntpMC_ldd->Project("h_mc_ldd", "tht", "moth==7 && pdg==-2212 && " + tht_cut_mc);
		ntpMC->Project("h_mc", "tht", "moth==7 && pdg==-2212 && " + tht_cut_mc);

		double nevts_mc_ldd = h_mc_ldd->GetEntries();
		double nevts_mc = h_mc->GetEntries();

		double eff_ldd = (nevts_mc_ldd==0)? 0: nevts_reco_ldd/nevts_mc_ldd;
		double eff = (nevts_mc==0)? 0: nevts_reco/nevts_mc;

		h_eff_ldd->Fill(theta, eff_ldd);
		h_eff->Fill(theta, eff);

	}
	h_eff->GetYaxis()->SetRangeUser(0,1.1);
	h_eff_ldd->GetYaxis()->SetRangeUser(0,1.1);
	jenny::CreateDrawAndSaveNHistograms(h_eff_ldd, h_eff, "with Lambda Discs", "without Lambda Discs", "/home/ikp1/puetz/panda/myscripts/simChain/test/plots/", "antiproton_new_lambda_disk_minus20percent_z", true, false);

}



