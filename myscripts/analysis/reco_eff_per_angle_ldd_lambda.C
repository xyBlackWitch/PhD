/**
 * @file reco_eff_per_angle_ldd_lambda.C
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

reco_eff_per_angle_ldd_lambda(TString pre1 ="", TString pre2=""){

	TFile * data_ldd = new TFile(pre1, "READ");
	TFile * data = new TFile(pre2, "READ");

	TTree * ntpMC_ldd = (TTree*) data_ldd->Get("ntpMC");
	TTree * ntpAntiLambda0_ldd = (TTree*) data_ldd->Get("ntpAntiLambda0");

	TTree * ntpMC = (TTree*) data->Get("ntpMC");
	TTree * ntpAntiLambda0 = (TTree*) data->Get("ntpAntiLambda0");

	TString tht_cut_alam;
	TString tht_cut_mc;

	TH1D * h_mc_ldd = new TH1D("h_mc_ldd", "h_mc_ldd", 25, 0,25);
	TH1D * h_alam_ldd = new TH1D("h_alam_ldd", "h_alam_ldd", 25,0,25);

	TH1D * h_mc = new TH1D("h_mc", "h_mc", 25, 0,25);
	TH1D * h_alam = new TH1D("h_alam", "h_alam", 25,0,25);

	TH1D * h_eff_ldd = new TH1D("h_eff_ldd", " reco efficiency for #bar{#Lambda}; #Theta/Deg; efficiency",25,0,25);
	TH1D * h_eff = new TH1D("h_eff", "reco efficiency; #Theta/Deg; efficiency",25,0,25);

	for(double theta = 0; theta<25; theta++){


		double tht_min =theta/180*TMath::Pi();
		double tht_max = (theta+1)/180*TMath::Pi();

		tht_cut_alam.Form("McTruth_tht>= %f && McTruth_tht< %f", tht_min, tht_max);
		tht_cut_mc.Form("tht>=%f && tht< %f", tht_min, tht_max);

		ntpAntiLambda0_ldd->Project("h_alam_ldd", "McTruth_tht*180/pi", "McTruthMatch  && antiLambda0_HitTag && VtxFit_HowGood==1 && MassFit_prob>0.01 && " + tht_cut_alam);
		ntpAntiLambda0->Project("h_alam", "McTruth_tht*180/pi", "McTruthMatch  && antiLambda0_HitTag && VtxFit_HowGood==1 && MassFit_prob>0.01 && " + tht_cut_alam);

		double nevts_reco_ldd = h_alam_ldd->GetEntries();
		double nevts_reco = h_alam->GetEntries();

		ntpMC_ldd->Project("h_mc_ldd", "tht", "moth==2 && pdg==-3122 && " + tht_cut_mc);
		ntpMC->Project("h_mc", "tht", "moth==2 && pdg==-3122 && " + tht_cut_mc);

		double nevts_mc_ldd = h_mc_ldd->GetEntries();
		double nevts_mc = h_mc->GetEntries();

		double eff_ldd = (nevts_mc_ldd==0)? 0: nevts_reco_ldd/nevts_mc_ldd;
		double eff = (nevts_mc==0)? 0: nevts_reco/nevts_mc;

		h_eff_ldd->Fill(theta, eff_ldd);
		h_eff->Fill(theta, eff);

	}
	h_eff->GetYaxis()->SetRangeUser(0,1.1);
	h_eff_ldd->GetYaxis()->SetRangeUser(0,1.1);
	jenny::CreateDrawAndSaveNHistograms(h_eff_ldd, h_eff, "with Lambda Discs", "without Lambda Discs", "/home/ikp1/puetz/panda/myscripts/simChain/test/plots/", "lambda_old_lambda_disk", false, false);

}



