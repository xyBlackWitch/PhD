/**
 * @file eval_final_states.C
 * @mainpage eval_final_states.C Get Data from analysis file and creates and saves different histograms
 *
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Get Data from analysis file and create and save different histograms
 * @details This files get the data for final state particle from the analysis of
 * pbar p -> Xi+ Xi(1820)-
 * 			 |	  |-> Lambda0 + K-
 * 			 |			|-> p + Pi-
 * 			 |
 * 			 |->  AntiLambda0 + Pi+
 * 			 		|-> pbar + Pi+
 * 	 Different useful histogramms are created and saved.
 */

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "/home/ikp1/puetz/panda/myscripts/common_jenny.cpp"


void eval_final_states_antiparticle(TString prefix="",  TString path="", bool save=kTRUE, bool close=kFALSE){

//	TString prefix;
//	if (path != "") prefix = path+"/"+pre;
//	else prefix = pre;



	//*** Input file
	TString inFile = TString::Format("%s/output_ana.root", prefix.Data());


	//*** get Data from Tree
	TFile * data = new TFile(inFile, "READ");

	TTree * ntpPiMinus = (TTree*)data->Get("ntpPiMinus");
	TTree * ntpPiPlus = (TTree*) data->Get("ntpPiPlus");
	TTree * ntpKaonMinus = (TTree*) data->Get("ntpKaonMinus");
	TTree * ntpKaonPlus = (TTree*) data->Get("ntpkaonplus");
	TTree * ntpProton = (TTree*) data->Get("ntpProton");
	TTree * ntpAntiProton = (TTree*) data->Get("ntpAntiProton");


	//***************************************************************************************************************
	//**** Get information about pi+ ********************************************************************************

	gStyle->SetOptStat(1111);
	gStyle->SetStatX(0.85);
	gStyle->SetStatY(0.85);

	TH2D * h_pimius_pt_vs_pz = new TH2D("h_piplus_pt_vs_pz", "Transverse vs. longitudinal momentum for #pi^{+}; pz/GeV/c; pt/GeV/c", 200,-1,2.5, 200,0,1);
	ntpPiPlus->Project("h_piplus_pt_vs_pz", "piplus_MC_pt: piplus_MC_pz", "McTruthMatch && piplus_HitTag");
	jenny::CreateDrawAndSaveHistogram(h_piplus_pt_vs_pz, prefix+"/plots/", "piplus_pt_vs_pz",save, close);

	TH1D * h_piplus_costht = new TH1D("h_piplus_costht", "cos(#Theta) distribution for #pi^{+}; cos(#Theta); counts", 200, -1,1.05);
	ntpPiPlus->Project("h_piplus_costht", "cos(piplus_tht)", "McTruthMatch && piplus_HitTag");
	jenny::CreateDrawAndSaveHistogram(h_piplus_costht, prefix+"/plots/", "piplus_costht",save, close);

	TH1D * h_piplus_MC_costht = new TH1D("h_piplus_MC_costht", "generated cos(#Theta) distribution for #pi^{+}; cos(#Theta); counts", 200, -1,1.05);
	ntpPiPlus->Project("h_piplus_MC_costht", "piplus_MC_CosTheta", "McTruthMatch && piplus_HitTag");
	jenny::CreateDrawAndSaveHistogram(h_piplus_MC_costht, prefix+"/plots/", "piplus_MC_costht",save, close);

	TH1D * h_piplus_number_per_evt = new TH1D("h_piplus_number_per_evt", "number of #pi^{+} per Event; # #pi^{+}/Event; counts", 7,0,7);
	ntpPiPlus->Project("h_piplus_number_per_evt", "cand", "McTruthMatch && piplus_HitTag");
	jenny::CreateDrawAndSaveHistogram(h_piplus_number_per_evt, prefix+"/plots/", "piplus_number_per_evt",save, close);


	//**** Get information about pi-(Lambda0) ********************************************************************


	TH2D * h_piminus_Al0_pt_vs_pz = new TH2D("h_piminus_Al0_pt_vs_pz", "Transverse vs. longitudinal momentum for #pi^{-}(#Lambda^{0}); pz/GeV/c; pt/GeV/c", 200,-0.1,0.7, 200,0,0.2);
	ntpPiMinus->Project("h_piminus_Al0_pt_vs_pz", "piminus_MC_pt: piminus_MC_pz", "McTruthMatch && piminus_HitTag && Mother==3122");
	jenny::CreateDrawAndSaveHistogram(h_piminus_Al0_pt_vs_pz, prefix+"/plots/", "piminus_Al0_pt_vs_pz",save, close);


	TH1D * h_piminus_Al0_costht = new TH1D("h_piminus_Al0_costht", "cos(#Theta) distribution for #pi^{-}(#Lambda^{0}); cos(#Theta); counts", 200, -1,1.05);
	ntpPiMinus->Project("h_piminus_Al0_costht", "cos(piminus_tht)", "McTruthMatch && piminus_HitTag && Mother==3122");
	jenny::CreateDrawAndSaveHistogram(h_piminus_Al0_costht, prefix+"/plots/", "piminus_Al0_costht",save, close);

	TH1D * h_piminus_Al0_MC_costht = new TH1D("h_piminus_Al0_MC_costht", "generated cos(#Theta) distribution for #pi^{-}(#Lambda^{0}); cos(#Theta); counts", 200, -1,1.05);
	ntpPiMinus->Project("h_piminus_Al0_MC_costht", "piminus_MC_CosTheta", "McTruthMatch && piminus_HitTag && Mother==3122");
	jenny::CreateDrawAndSaveHistogram(h_piminus_Al0_MC_costht, prefix+"/plots/", "piminus_Al0_MC_costht",save, close);

	TH1D * h_piminus_Al0_number_per_evt = new TH1D("h_piminus_Al0_number_per_evt", "number of #pi^{-}(#Lambda^{0}) per Event; # #pi^{-}(#Lambda^{0})/Event; counts", 7,0,7);
	ntpPiMinus->Project("h_piminus_Al0_number_per_evt", "cand", "McTruthMatch && piminus_HitTag && Mother==3122");
	jenny::CreateDrawAndSaveHistogram(h_piminus_Al0_number_per_evt, prefix+"/plots/", "piminus_Al0_number_per_evt",save, close);



	//**** Get information about pi-(Xi-) ********************************************************************


	TH2D * h_piminus_XiPlus_pt_vs_pz = new TH2D("h_piminus_XiPlus_pt_vs_pz", "Transverse vs. longitudinal momentum for #pi^{-}(#Xi); pz/GeV/c; pt/GeV/c", 200,-0.1,1, 200,0,0.3);
	ntpPiMinus->Project("h_piminus_XiPlus_pt_vs_pz", "piminus_MC_pt: piminus_MC_pz", "McTruthMatch && piminus_HitTag && Mother==3312");
	jenny::CreateDrawAndSaveHistogram(h_piminus_XiPlus_pt_vs_pz, prefix+"/plots/", "piminus_XiPlus_pt_vs_pz",save, close);


	TH1D * h_piminus_XiPlus_costht = new TH1D("h_piminus_XiPlus_costht", "cos(#Theta) distribution for #pi^{-}(#Xi); cos(#Theta); counts", 200, -1,1.05);
	ntpPiMinus->Project("h_piminus_XiPlus_costht", "cos(piminus_tht)", "McTruthMatch && piminus_HitTag && Mother==3312");
	jenny::CreateDrawAndSaveHistogram(h_piminus_XiPlus_costht, prefix+"/plots/", "piminus_XiPlus_costht",save, close);

	TH1D * h_piminus_XiPlus_MC_costht = new TH1D("h_piminus_XiPlus_MC_costht", "generated cos(#Theta) distribution for #pi^{-}(#Xi); cos(#Theta); counts", 200, -1,1.05);
	ntpPiMinus->Project("h_piminus_XiPlus_MC_costht", "piminus_MC_CosTheta", "McTruthMatch && piminus_HitTag && Mother==3312");
	jenny::CreateDrawAndSaveHistogram(h_piminus_XiPlus_MC_costht, prefix+"/plots/", "piminus_XiPlus_MC_costht",save, close);

	TH1D * h_piminus_XiPlus_number_per_evt = new TH1D("h_piminus_XiPlus_number_per_evt", "number of #pi^{-}(#Xi) per Event; # #pi^{-}(#Xi)/Event; counts", 7,0,7);
	ntpPiMinus->Project("h_piminus_XiPlus_number_per_evt", "cand", "McTruthMatch && piminus_HitTag && Mother==3312");
	jenny::CreateDrawAndSaveHistogram(h_piminus_XiPlus_number_per_evt, prefix+"/plots/", "piminus_XiPlus_number_per_evt",save, close);



	//**** Get information about K- or K+ *************************************************************************

	if (ntpKaonMinus!=0x0){
		TH2D * h_KaonMinus_pt_vs_pz = new TH2D("h_KaonMinus_pt_vs_pz", "Transverse vs. longitudinal momentum for K^{-}; pz/GeV/c; pt/GeV/c", 200,-0.5,2.1, 200,0,0.7);
		ntpKaonMinus->Project("h_KaonMinus_pt_vs_pz", "kaonminus_MC_pt: kaonminus_MC_pz", "McTruthMatch && kaonminus_HitTag");
		jenny::CreateDrawAndSaveHistogram(h_KaonMinus_pt_vs_pz, prefix+"/plots/", "kaonminus_pt_vs_pz",save, close);


		TH1D * h_KaonMinus_costht = new TH1D("h_KaonMinus_costht", "cos(#Theta) distribution for K^{-}; cos(#Theta); counts", 200, -1,1.05);
		ntpKaonMinus->Project("h_KaonMinus_costht", "cos(kaonminus_tht)", "McTruthMatch && kaonminus_HitTag");
		jenny::CreateDrawAndSaveHistogram(h_KaonMinus_costht, prefix+"/plots/", "kaonminus_costht",save, close);

		TH1D * h_KaonMinus_MC_costht = new TH1D("h_KaonMinus_MC_costht", "generated cos(#Theta) distribution for K^{-}; cos(#Theta); counts", 200, -1,1.05);
		ntpKaonMinus->Project("h_KaonMinus_MC_costht", "kaonminus_MC_CosTheta", "McTruthMatch && kaonminus_HitTag");
		jenny::CreateDrawAndSaveHistogram(h_KaonMinus_MC_costht, prefix+"/plots/", "kaonminus_MC_costht",save, close);

		TH1D * h_KaonMinus_number_per_evt = new TH1D("h_KaonMinus_number_per_evt", "number of K^{-} per Event; # K^{-}/Event; counts", 7,0,7);
		ntpKaonMinus->Project("h_KaonMinus_number_per_evt", "cand", "McTruthMatch && kaonminus_HitTag");
		jenny::CreateDrawAndSaveHistogram(h_KaonMinus_number_per_evt, prefix+"/plots/", "kaonminus_number_per_evt",save, close);
	}
	else if(ntpKaonPlus!=0x0){
		TH2D * h_KaonPlus_pt_vs_pz = new TH2D("h_KaonPlus_pt_vs_pz", "Transverse vs. longitudinal momentum for K^{+}; pz/GeV/c; pt/GeV/c", 200,0,2.1, 200,0,0.7);
		ntpKaonPlus->Project("h_KaonPlus_pt_vs_pz", "kaonplus_MC_pt: kaonplus_MC_pz", "McTruthMatch && kaonplus_HitTag");
		jenny::CreateDrawAndSaveHistogram(h_KaonPlus_pt_vs_pz, prefix+"/plots/", "KaonPlus_pt_vs_pz",save, close);


		TH1D * h_KaonPlus_costht = new TH1D("h_KaonPlus_costht", "cos(#Theta) distribution for K^{+}; cos(#Theta); counts", 200, -1,1.05);
		ntpKaonPlus->Project("h_KaonPlus_costht", "cos(kaonplus_tht)", "McTruthMatch && kaonplus_HitTag");
		jenny::CreateDrawAndSaveHistogram(h_KaonPlus_costht, prefix+"/plots/", "KaonPlus_costht",save, close);

		TH1D * h_KaonPlus_MC_costht = new TH1D("h_KaonPlus_MC_costht", "generated cos(#Theta) distribution for K^{+}; cos(#Theta); counts", 200, -1,1.05);
		ntpKaonPlus->Project("h_KaonPlus_MC_costht", "kaonplus_MC_CosTheta", "McTruthMatch && kaonplus_HitTag");
		jenny::CreateDrawAndSaveHistogram(h_KaonPlus_MC_costht, prefix+"/plots/", "KaonPlus_MC_costht",save, close);

		TH1D * h_KaonPlus_number_per_evt = new TH1D("h_KaonPlus_number_per_evt", "number of K^{+} per Event; # K^{+}/Event; counts", 7,0,7);
		ntpKaonPlus->Project("h_KaonPlus_number_per_evt", "cand", "McTruthMatch && kaonplus_HitTag");
		jenny::CreateDrawAndSaveHistogram(h_KaonPlus_number_per_evt, prefix+"/plots/", "KaonPlus_number_per_evt",save, close);
	}
	else{
		cout << "No particle of kind K- or K+!" << endl;
	}


	//**** Get information about p   *************************************************************************
	TH2D * h_Proton_pt_vs_pz = new TH2D("h_Proton_pt_vs_pz", "Transverse vs. longitudinal momentum for p; pz/GeV/c; pt/GeV/c", 200,-0.1,2.7, 200,0,0.8);
	ntpProton->Project("h_Proton_pt_vs_pz", "proton_MC_pt: proton_MC_pz", "McTruthMatch && proton_HitTag && Mother==3122");
	jenny::CreateDrawAndSaveHistogram(h_Proton_pt_vs_pz, prefix+"/plots/", "proton_pt_vs_pz",save, close);


	TH1D * h_Proton_costht = new TH1D("h_Proton_costht", "cos(#Theta) distribution for p; cos(#Theta); counts", 200, -1,1.05);
	ntpProton->Project("h_Proton_costht", "cos(proton_tht)", "McTruthMatch && proton_HitTag && Mother==3122");
	jenny::CreateDrawAndSaveHistogram(h_Proton_costht, prefix+"/plots/", "proton_costht",save, close);

//	TH1D * h_Proton_MC_costht = new TH1D("h_Proton_MC_costht", "generated cos(#Theta) distribution for p; cos(#Theta); counts", 200, -1,1.05);
//	ntpProton->Project("h_Proton_MC_costht", "proton_MC_CosTheta", "McTruthMatch && proton_HitTag && Mother==3122");
//	jenny::CreateDrawAndSaveHistogram(h_Proton_MC_costht, prefix+"/plots/", "proton_MC_costht",save, close);

	TH1D * h_Proton_number_per_evt = new TH1D("h_Proton_number_per_evt", "number of p per Event; # p/Event; counts", 10,0,10);
	ntpProton->Project("h_Proton_number_per_evt", "cand", "McTruthMatch && proton_HitTag && Mother==3122");
	jenny::CreateDrawAndSaveHistogram(h_Proton_number_per_evt, prefix+"/plots/", "proton_number_per_evt",save, close);

	//**** Get information about p+  *************************************************************************

	TH2D * h_antiProton_pt_vs_pz = new TH2D("h_antiProton_pt_vs_pz", "Transverse vs. longitudinal momentum for #bar{p}; pz/GeV/c; pt/GeV/c", 200,0.4,2.8, 200,0,0.7);
	ntpAntiProton->Project("h_antiProton_pt_vs_pz", "AntiProton_MC_pt: AntiProton_MC_pz", "McTruthMatch && AntiProton_HitTag && Mother==-3122");
	jenny::CreateDrawAndSaveHistogram(h_antiProton_pt_vs_pz, prefix+"/plots/", "AntiProton_pt_vs_pz",save, close);


	TH1D * h_antiProton_costht = new TH1D("h_antiProton_costht", "cos(#Theta) distribution for #bar{p}; cos(#Theta); counts", 200, -1,1.05);
	ntpAntiProton->Project("h_antiProton_costht", "cos(AntiProton_tht)", "McTruthMatch && AntiProton_HitTag && Mother==-3122");
	jenny::CreateDrawAndSaveHistogram(h_antiProton_costht, prefix+"/plots/", "AntiProton_costht",save, close);

	TH1D * h_antiProton_MC_costht = new TH1D("h_antiProton_MC_costht", "generated cos(#Theta) distribution for #bar{p}; cos(#Theta); counts", 200, -1,1.05);
	ntpAntiProton->Project("h_antiProton_MC_costht", "AntiProton_MC_CosTheta", "McTruthMatch && AntiProton_HitTag && Mother==-3122");
	jenny::CreateDrawAndSaveHistogram(h_antiProton_MC_costht, prefix+"/plots/", "AntiProton_MC_costht",save, close);

	TH1D * h_antiProton_number_per_evt = new TH1D("h_antiProton_number_per_evt", "number of #bar{p} per Event; # #bar{p}/Event; counts", 10,0,10);
	ntpAntiProton->Project("h_antiProton_number_per_evt", "cand", "McTruthMatch && AntiProton_HitTag && Mother==-3122");
	jenny::CreateDrawAndSaveHistogram(h_antiProton_number_per_evt, prefix+"/plots/", "AntiProton_number_per_evt",save, close);
}
