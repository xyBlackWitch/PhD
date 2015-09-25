/**
 * @file eval_cascade.C
 * @mainpage eval_cascade.C Get Data from analysis file and creates and saves different histograms
 *
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Get Data from analysis file and create and save different histograms
 * @details This files get the data from the analysis of
 * pbar p -> Xi+ Xi-
 * 			 |	  |-> Lambda0 + Pi-
 * 			 |			|-> p + Pi-
 * 			 |
 * 			 |->  AntiLambda0 + Pi+
 * 			 		|-> pbar + Pi+
 * 	 Different useful histogramms are created and saved.
 */

class RhoTuple;
class RhoCandidate;

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "../common_jenny.cpp"



void eval_cascade(bool save=true, bool close=false){

	//Input file
	TString inPath = "/private/puetz/mysimulations/analysis/pbarp_Xiplus_Ximinus/idealtracking/1000_events/momentum_2.7/run2/";
	TFile * inFile = new TFile(inPath+"output_ana.root", "READ");

	//Output file
	TString outPath = "";

	//Get NTuple
	TTree * ntpMC = (TTree*) inFile->Get("ntpMC");
	TTree * ntpPiMinus = (TTree*) inFile->Get("ntpPiMinus");
	TTree * ntpPiPlus = (TTree*) inFile->Get("ntpPiPlus");
	TTree * ntpProton = (TTree*) inFile->Get("ntpProton");
	TTree * ntpAntiProton= (TTree*) inFile->Get("ntpAntiProton");
	TTree * ntpLambda0 = (TTree*) inFile->Get("ntpLambda0");
	TTree * ntpAntiLambda0= (TTree*) inFile->Get("ntpAntiLambda0");
	TTree * ntpXiPlus = (TTree*) inFile->Get("ntpXiPlus");
	TTree * ntpXiMinus = (TTree*) inFile->Get("ntpXiMinus");
	TTree * ntpXiSys = (TTree*) inFile->Get("ntpXiSys");


	TString cut1 = "VtxFit_prob>0.01";
	TString cut2 = cut1 + "&& fMass_Prob>0.01";



	//************************** Projections for PiMinus from Xi- decay ********************************************

	//Create histograms
	TH1D * h_piminus_tht = new TH1D("h_piminus_tht", "#Theta distribution for #pi^{-}(#Xi^{-}); #theta/rad; counts" ,100,-1,1);
	ntpPiMinus->Project("h_piminus_tht", "piminus_tht", "McTruthMatch==1 && Mother==3312");

	TH1D * h_piminus_pz = new TH1D("h_piminus_pz", "distribution of longitudinal momentum for #pi^{-}(#Xi^{-}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiMinus->Project("h_piminus_pz", "piminus_pz", "McTruthMatch==1 && Mother==3312");

	TH1D * h_piminus_npart = new TH1D("h_piminus_npart", "number of reconstructed #pi^{-}(#Xi^{-}) per Event; reco. #pi^{-} per Event; counts", 7,1,7);
	ntpPiMinus->Project("h_piminus_npart", "ncand", "McTruthMatch==1 && Mother==3312");

	TH2D * h_piminus_pt_vs_pz = new TH2D("h_piminus_pt_vs_pz", "p_{t} vs p_{z} for #pi^{-}(#Xi^{-}); p_{z}/Gev/c; p_{t}/Gev/c", 100,-0.1,0.7,100,0,0.2);
	ntpPiMinus->Project("h_piminus_pt_vs_pz", "piminus_pt:piminus_pz", "McTruthMatch==1 && Mother==3312");

	TH2D * h_piminus_costht_vs_pz = new TH2D("h_piminus_costht_vs_pz", "cos(#theta) vs p_{z} for #pi^{-}(#Xi^{-}); p_{z}/Gev/c; cos(#theta)", 100,-0.1,0.7,100,-1.,1.);
	ntpPiMinus->Project("h_piminus_costht_vs_pz", "PiMinus_CosTheta:piminus_pz", "McTruthMatch==1 && Mother==3312");


	//historgrams for efficiency
	TH1D * h_piminus_MC_tht = new TH1D("h_piminus_MC_tht", "cos(#Theta) distribution for #pi^{-}_{MC}(#Xi^{-}); cos(#theta); counts" ,100,-1,1);
	ntpPiMinus->Project("h_piminus_MC_tht", "cos(piminus_MC_tht)", "McTruthMatch==1 && Mother==3312");

	TH1D * h_piminus_MC_pz = new TH1D("h_piminus_MC_pz", "distribution of longitudinal momentum for #pi^{-}_{MC}(#Xi^{-}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiMinus->Project("h_piminus_MC_pz", "piminus_MC_pz", "McTruthMatch==1 && Mother==3312");

	TH1D * h_MC_tht = new TH1D("h_MC_tht", "#Theta_{MC} distribution for #pi^{-}(#Xi^{-}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_MC_tht", "cos(tht)", "moth==1 && pdg==-211");

	TH1D * h_MC_pz = new TH1D("h_MC_pz", "distribution of longitudinal momentum (MC) for #pi^{-}(#Xi^{-}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_MC_pz", "pz", "moth==1 && pdg==-211");

	TH1D * h_eff_pz = new TH1D("h_eff_pz", "reconstruction efficiency for #pi^{-} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_eff_pz->Divide(h_piminus_MC_pz, h_MC_pz);
	h_eff_pz->SetStats(0);
	h_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_eff_tht = new TH1D("h_eff_tht", "reconstruction efficiency for #pi^{-} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_eff_tht->Divide(h_piminus_MC_tht, h_MC_tht);
	h_eff_tht->SetStats(0);
	h_eff_tht->GetYaxis()->SetRangeUser(0,1.1);




	//	save histograms
//	jenny::CreateDrawAndSaveHistogram(h_piminus_MC_tht, outPath, "h_piminus_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminus_pz, outPath, "h_piminus_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminus_npart, outPath, "h_piminus_npart", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminus_pt_vs_pz, outPath, "h_piminus_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminus_costht_vs_pz, outPath, "h_piminus_costht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_eff_tht, outPath, "h_eff_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_eff_pz, outPath, "h_eff_pz", save, close);



	//************************** Projections for piplus from Xi- decay ********************************************

	//Create histograms
	TH1D * h_piplus_tht = new TH1D("h_piplus_tht", "cos(#Theta) distribution for #pi^{+}(#Xi^{+}); cos(#theta); counts" ,100,-1,1);
	ntpPiPlus->Project("h_piplus_tht", "PiPlus_tht", "McTruthMatch==1 && Mother==-3312");

	TH1D * h_piplus_pz = new TH1D("h_piplus_pz", "distribution of longitudinal momentum for #pi^{+}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiPlus->Project("h_piplus_pz", "PiPlus_pz", "McTruthMatch==1 && Mother==-3312");

	TH1D * h_piplus_npart = new TH1D("h_piplus_npart", "number of reconstructed #pi^{+}(#Xi^{+}) per Event; reco. #pi^{+} per Event; counts", 7,1,7);
	ntpPiPlus->Project("h_piplus_npart", "ncand", "McTruthMatch==1 && Mother==-3312");

	TH2D * h_piplus_pt_vs_pz = new TH2D("h_piplus_pt_vs_pz", "p_{t} vs p_{z} for #pi^{+}(#Xi^{+}); p_{z}/Gev/c; p_{t}/Gev/c", 100,-0.1,0.7,100,0,0.2);
	ntpPiPlus->Project("h_piplus_pt_vs_pz", "PiPlus_pt:PiPlus_pz", "McTruthMatch==1 && Mother==-3312");

	TH2D * h_piplus_costht_vs_pz = new TH2D("h_piplus_costht_vs_pz", "cos(#theta) vs p_{z} for #pi^{+}(#Xi^{+}); p_{z}/Gev/c; cos(#theta)", 100,-0.1,0.7,100,-1.,1.);
	ntpPiPlus->Project("h_piplus_costht_vs_pz", "PiPlus_CosTheta:PiPlus_pz", "McTruthMatch==1 && Mother==-3312");


	//historgrams for efficiency
	TH1D * h_piplus_McTruth_tht = new TH1D("h_piplus_McTruth_tht", "#Theta distribution for #pi^{+}_{MC}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpPiPlus->Project("h_piplus_McTruth_tht", "cos(PiPlus_MC_tht)", "McTruthMatch==1 && Mother==-3312");

	TH1D * h_piplus_McTruth_pz = new TH1D("h_piplus_McTruth_pz", "distribution of longitudinal momentum for #pi^{+}_{MC}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiPlus->Project("h_piplus_McTruth_pz", "PiPlus_MC_pz", "McTruthMatch==1 && Mother==-3312");

	TH1D * h_piplus_MC_tht = new TH1D("h_piplus_MC_tht", "#Theta_{MC} distribution for #pi^{+}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_piplus_MC_tht", "cos(tht)", "moth==2 && pdg==211");

	TH1D * h_piplus_MC_pz = new TH1D("h_piplus_MC_pz", "distribution of longitudinal momentum (MC) for #pi^{+}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_piplus_MC_pz", "pz", "moth==2 && pdg==211");

	TH1D * h_piplus_eff_pz = new TH1D("h_piplus_eff_pz", "reconstruction efficiency for #pi^{+} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_piplus_eff_pz->Divide(h_piplus_McTruth_pz, h_piplus_MC_pz);
	h_piplus_eff_pz->SetStats(0);
	h_piplus_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_piplus_eff_tht = new TH1D("h_piplus_eff_tht", "reconstruction efficiency for #pi^{+} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_piplus_eff_tht->Divide(h_piplus_McTruth_tht, h_piplus_MC_tht);
	h_piplus_eff_tht->SetStats(0);
	h_piplus_eff_tht->GetYaxis()->SetRangeUser(0,1.1);




	//save histograms
//	jenny::CreateDrawAndSaveHistogram(h_piplus_MC_tht, outPath, "h_piplus_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_pz, outPath, "h_piplus_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_npart, outPath, "h_piplus_npart", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_pt_vs_pz, outPath, "h_piplus_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_costht_vs_pz, outPath, "h_piplus_costht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_eff_tht, outPath, "h_piplus_eff_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piplus_eff_pz, outPath, "h_piplus_eff_pz", save, close);


	//*********************** plots for comparison of pi+ and pi- *********************************************

//	jenny::CreateDrawAndSaveNHistograms(h_piminus_npart, h_piplus_npart, "#pi^{-}","#pi^{+}",  outPath, "h_piplus_piminus_npart", save, close);




	//************************** Projections for proton from Lambda0 decay ********************************************

	//Create histograms
	TH1D * h_proton_tht = new TH1D("h_proton_tht", "#Theta distribution for #pi^{-}(#Lambda^{0}); #theta/rad; counts" ,100,-1,1);
	ntpProton->Project("h_proton_tht", "proton_tht", "McTruthMatch==1 && Mother==3122");

	TH1D * h_proton_pz = new TH1D("h_proton_pz", "distribution of longitudinal momentum for #pi^{-}(#Lambda^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpProton->Project("h_proton_pz", "proton_pz", "McTruthMatch==1 && Mother==3122");

	TH1D * h_proton_npart = new TH1D("h_proton_npart", "number of reconstructed #pi^{-}(#Lambda^{0}) per Event; reco. #pi^{-} per Event; counts", 7,1,7);
	ntpProton->Project("h_proton_npart", "ncand", "McTruthMatch==1 && Mother==3122");

	TH2D * h_proton_pt_vs_pz = new TH2D("h_proton_pt_vs_pz", "p_{t} vs p_{z} for #pi^{-}(#Lambda^{0}); p_{z}/Gev/c; p_{t}/Gev/c", 100,-0.1,0.7,100,0,0.2);
	ntpProton->Project("h_proton_pt_vs_pz", "proton_pt:proton_pz", "McTruthMatch==1 && Mother==3122");

	TH2D * h_proton_costht_vs_pz = new TH2D("h_proton_costht_vs_pz", "cos(#theta) vs p_{z} for #pi^{-}(#Lambda^{0}); p_{z}/Gev/c; cos(#theta)", 100,-0.1,0.7,100,-1.,1.);
	ntpProton->Project("h_proton_costht_vs_pz", "proton_CosTheta:proton_pz", "McTruthMatch==1 && Mother==3122");


	//historgrams for efficiency
	TH1D * h_proton_MC_tht = new TH1D("h_proton_MC_tht", "cos(#Theta) distribution for #pi^{-}_{MC}(#Lambda^{0}); cos(#theta); counts" ,100,-1,1);
	ntpProton->Project("h_proton_MC_tht", "cos(proton_MC_tht)", "McTruthMatch==1 && Mother==3122");

	TH1D * h_proton_MC_pz = new TH1D("h_proton_MC_pz", "distribution of longitudinal momentum for #pi^{-}_{MC}(#Lambda^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpProton->Project("h_proton_MC_pz", "proton_MC_pz", "McTruthMatch==1 && Mother==3122");

	TH1D * h_MC_tht = new TH1D("h_MC_tht", "#Theta_{MC} distribution for #pi^{-}(#Lambda^{0}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_MC_tht", "cos(tht)", "moth==3 && pdg==2212");

	TH1D * h_MC_pz = new TH1D("h_MC_pz", "distribution of longitudinal momentum (MC) for #pi^{-}(#Lambda^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_MC_pz", "pz", "moth==3 && pdg==2212");

	TH1D * h_eff_pz = new TH1D("h_eff_pz", "reconstruction efficiency for #pi^{-} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_eff_pz->Divide(h_proton_MC_pz, h_MC_pz);
	h_eff_pz->SetStats(0);
	h_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_eff_tht = new TH1D("h_eff_tht", "reconstruction efficiency for #pi^{-} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_eff_tht->Divide(h_proton_MC_tht, h_MC_tht);
	h_eff_tht->SetStats(0);
	h_eff_tht->GetYaxis()->SetRangeUser(0,1.1);




	//	save histograms
//	jenny::CreateDrawAndSaveHistogram(h_proton_MC_tht, outPath, "h_proton_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_proton_pz, outPath, "h_proton_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_proton_npart, outPath, "h_proton_npart", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_proton_pt_vs_pz, outPath, "h_proton_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_proton_costht_vs_pz, outPath, "h_proton_costht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_eff_tht, outPath, "h_eff_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_eff_pz, outPath, "h_eff_pz", save, close);



	//************************** Projections for antiProton from AntiLambda0 decay ********************************************

	//Create histograms
	TH1D * h_antiProton_tht = new TH1D("h_antiProton_tht", "cos(#Theta) distribution for #pi^{+}(#bar{#Lambda}^{0}); cos(#theta); counts" ,100,-1,1);
	ntpAntiProton->Project("h_antiProton_tht", "antiProton_tht", "McTruthMatch==1 && Mother==-3122");

	TH1D * h_antiProton_pz = new TH1D("h_antiProton_pz", "distribution of longitudinal momentum for #pi^{+}(#bar{#Lambda}^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpAntiProton->Project("h_antiProton_pz", "antiProton_pz", "McTruthMatch==1 && Mother==-3122");

	TH1D * h_antiProton_npart = new TH1D("h_antiProton_npart", "number of reconstructed #pi^{+}(#bar{#Lambda}^{0}) per Event; reco. #pi^{+} per Event; counts", 7,1,7);
	ntpAntiProton->Project("h_antiProton_npart", "ncand", "McTruthMatch==1 && Mother==-3122");

	TH2D * h_antiProton_pt_vs_pz = new TH2D("h_antiProton_pt_vs_pz", "p_{t} vs p_{z} for #pi^{+}(#bar{#Lambda}^{0}); p_{z}/Gev/c; p_{t}/Gev/c", 100,-0.1,0.7,100,0,0.2);
	ntpAntiProton->Project("h_antiProton_pt_vs_pz", "antiProton_pt:antiProton_pz", "McTruthMatch==1 && Mother==-3122");

	TH2D * h_antiProton_costht_vs_pz = new TH2D("h_antiProton_costht_vs_pz", "cos(#theta) vs p_{z} for #pi^{+}(#bar{#Lambda}^{0}); p_{z}/Gev/c; cos(#theta)", 100,-0.1,0.7,100,-1.,1.);
	ntpAntiProton->Project("h_antiProton_costht_vs_pz", "antiProton_CosTheta:antiProton_pz", "McTruthMatch==1 && Mother==-3122");


	//historgrams for efficiency
	TH1D * h_antiProton_McTruth_tht = new TH1D("h_antiProton_McTruth_tht", "#Theta distribution for #pi^{+}_{MC}(#bar{#Lambda}^{0}); #theta/rad; counts" ,100,-1,1);
	ntpAntiProton->Project("h_antiProton_McTruth_tht", "cos(antiProton_MC_tht)", "McTruthMatch==1 && Mother==-3122");

	TH1D * h_antiProton_McTruth_pz = new TH1D("h_antiProton_McTruth_pz", "distribution of longitudinal momentum for #pi^{+}_{MC}(#bar{#Lambda}^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpAntiProton->Project("h_antiProton_McTruth_pz", "antiProton_MC_pz", "McTruthMatch==1 && Mother==-3122");

	TH1D * h_antiProton_MC_tht = new TH1D("h_antiProton_MC_tht", "#Theta_{MC} distribution for #pi^{+}(#bar{#Lambda}^{0}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_antiProton_MC_tht", "cos(tht)", "moth==5 && pdg==-2212");

	TH1D * h_antiProton_MC_pz = new TH1D("h_antiProton_MC_pz", "distribution of longitudinal momentum (MC) for #pi^{+}(#bar{#Lambda}^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_antiProton_MC_pz", "pz", "moth==5 && pdg==-2212");

	TH1D * h_antiProton_eff_pz = new TH1D("h_antiProton_eff_pz", "reconstruction efficiency for #pi^{+} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_antiProton_eff_pz->Divide(h_antiProton_McTruth_pz, h_antiProton_MC_pz);
	h_antiProton_eff_pz->SetStats(0);
	h_antiProton_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_antiProton_eff_tht = new TH1D("h_antiProton_eff_tht", "reconstruction efficiency for #pi^{+} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_antiProton_eff_tht->Divide(h_antiProton_McTruth_tht, h_antiProton_MC_tht);
	h_antiProton_eff_tht->SetStats(0);
	h_antiProton_eff_tht->GetYaxis()->SetRangeUser(0,1.1);




	//save histograms
//	jenny::CreateDrawAndSaveHistogram(h_antiProton_MC_tht, outPath, "h_antiProton_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiProton_pz, outPath, "h_antiProton_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiProton_npart, outPath, "h_antiProton_npart", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiProton_pt_vs_pz, outPath, "h_antiProton_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiProton_costht_vs_pz, outPath, "h_antiProton_costht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiProton_eff_tht, outPath, "h_antiProton_eff_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiProton_eff_pz, outPath, "h_antiProton_eff_pz", save, close);


	//*********************** plots for comparison of p and pbar *********************************************

//	jenny::CreateDrawAndSaveNHistograms(h_proton_npart, h_antiProton_npart, "p","#bar{p}",  outPath, "h_antiProton_proton_npart", save, close);



	//************************** Projections for PiMinus from Xi- decay ********************************************

	//Create histograms
	TH1D * h_piminuslam_tht = new TH1D("h_piminuslam_tht", "#Theta distribution for #pi^{-}(#Lambda^{0}); #theta/rad; counts" ,100,-1,1);
	ntpPiMinus->Project("h_piminuslam_tht", "piminus_tht", "McTruthMatch==1 && Mother==3122");

	TH1D * h_piminuslam_pz = new TH1D("h_piminuslam_pz", "distribution of longitudinal momentum for #pi^{-}(#Lambda^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiMinus->Project("h_piminuslam_pz", "piminus_pz", "McTruthMatch==1 && Mother==3122");

	TH1D * h_piminuslam_npart = new TH1D("h_piminuslam_npart", "number of reconstructed #pi^{-}(#Lambda^{0}) per Event; reco. #pi^{-} per Event; counts", 7,1,7);
	ntpPiMinus->Project("h_piminuslam_npart", "ncand", "McTruthMatch==1 && Mother==3122");

	TH2D * h_piminuslam_pt_vs_pz = new TH2D("h_piminuslam_pt_vs_pz", "p_{t} vs p_{z} for #pi^{-}(#Lambda^{0}); p_{z}/Gev/c; p_{t}/Gev/c", 100,-0.1,0.7,100,0,0.2);
	ntpPiMinus->Project("h_piminuslam_pt_vs_pz", "piminus_pt:piminus_pz", "McTruthMatch==1 && Mother==3122");

	TH2D * h_piminuslam_costht_vs_pz = new TH2D("h_piminuslam_costht_vs_pz", "cos(#theta) vs p_{z} for #pi^{-}(#Lambda^{0}); p_{z}/Gev/c; cos(#theta)", 100,-0.1,0.7,100,-1.,1.);
	ntpPiMinus->Project("h_piminuslam_costht_vs_pz", "PiMinus_CosTheta:piminus_pz", "McTruthMatch==1 && Mother==3122");


	//historgrams for efficiency
	TH1D * h_piminuslam_MC_tht = new TH1D("h_piminuslam_MC_tht", "cos(#Theta) distribution for #pi^{-}_{MC}(#Lambda^{0}); cos(#theta); counts" ,100,-1,1);
	ntpPiMinus->Project("h_piminuslam_MC_tht", "cos(piminus_MC_tht)", "McTruthMatch==1 && Mother==3122");

	TH1D * h_piminuslam_MC_pz = new TH1D("h_piminuslam_MC_pz", "distribution of longitudinal momentum for #pi^{-}_{MC}(#Lambda^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiMinus->Project("h_piminuslam_MC_pz", "piminus_MC_pz", "McTruthMatch==1 && Mother==3122");

	TH1D * h_MC_tht = new TH1D("h_MC_tht", "#Theta_{MC} distribution for #pi^{-}(#Lambda^{0}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_MC_tht", "cos(tht)", "moth==3 && pdg==-211");

	TH1D * h_MC_pz = new TH1D("h_MC_pz", "distribution of longitudinal momentum (MC) for #pi^{-}(#Lambda^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_MC_pz", "pz", "moth==3 && pdg==-211");

	TH1D * h_eff_pz = new TH1D("h_eff_pz", "reconstruction efficiency for #pi^{-} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_eff_pz->Divide(h_piminuslam_MC_pz, h_MC_pz);
	h_eff_pz->SetStats(0);
	h_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_eff_tht = new TH1D("h_eff_tht", "reconstruction efficiency for #pi^{-} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_eff_tht->Divide(h_piminuslam_MC_tht, h_MC_tht);
	h_eff_tht->SetStats(0);
	h_eff_tht->GetYaxis()->SetRangeUser(0,1.1);




	//	save histograms
//	jenny::CreateDrawAndSaveHistogram(h_piminuslam_MC_tht, outPath, "h_piminuslam_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminuslam_pz, outPath, "h_piminuslam_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminuslam_npart, outPath, "h_piminuslam_npart", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminuslam_pt_vs_pz, outPath, "h_piminuslam_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_piminuslam_costht_vs_pz, outPath, "h_piminuslam_costht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_eff_tht, outPath, "h_eff_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_eff_pz, outPath, "h_eff_pz", save, close);



	//************************** Projections for piplus from Xi- decay ********************************************

	//Create histograms
	TH1D * h_pipluslam_tht = new TH1D("h_pipluslam_tht", "cos(#Theta) distribution for #pi^{+}(#bar{#Lambda}^{0}); cos(#theta); counts" ,100,-1,1);
	ntpPiPlus->Project("h_pipluslam_tht", "PiPlus_tht", "McTruthMatch==1 && Mother==-3122");

	TH1D * h_pipluslam_pz = new TH1D("h_pipluslam_pz", "distribution of longitudinal momentum for #pi^{+}(#bar{#Lambda}^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiPlus->Project("h_pipluslam_pz", "PiPlus_pz", "McTruthMatch==1 && Mother==-3122");

	TH1D * h_pipluslam_npart = new TH1D("h_pipluslam_npart", "number of reconstructed #pi^{+}(#bar{#Lambda}^{0}) per Event; reco. #pi^{+} per Event; counts", 7,1,7);
	ntpPiPlus->Project("h_pipluslam_npart", "ncand", "McTruthMatch==1 && Mother==-3122");

	TH2D * h_pipluslam_pt_vs_pz = new TH2D("h_pipluslam_pt_vs_pz", "p_{t} vs p_{z} for #pi^{+}(#bar{#Lambda}^{0}); p_{z}/Gev/c; p_{t}/Gev/c", 100,-0.1,0.7,100,0,0.2);
	ntpPiPlus->Project("h_pipluslam_pt_vs_pz", "PiPlus_pt:PiPlus_pz", "McTruthMatch==1 && Mother==-3122");

	TH2D * h_pipluslam_costht_vs_pz = new TH2D("h_pipluslam_costht_vs_pz", "cos(#theta) vs p_{z} for #pi^{+}(#bar{#Lambda}^{0}); p_{z}/Gev/c; cos(#theta)", 100,-0.1,0.7,100,-1.,1.);
	ntpPiPlus->Project("h_pipluslam_costht_vs_pz", "PiPlus_CosTheta:PiPlus_pz", "McTruthMatch==1 && Mother==-3122");


	//historgrams for efficiency
	TH1D * h_pipluslam_McTruth_tht = new TH1D("h_pipluslam_McTruth_tht", "#Theta distribution for #pi^{+}_{MC}(#bar{#Lambda}^{0}); #theta/rad; counts" ,100,-1,1);
	ntpPiPlus->Project("h_pipluslam_McTruth_tht", "cos(PiPlus_MC_tht)", "McTruthMatch==1 && Mother==-3122");

	TH1D * h_pipluslam_McTruth_pz = new TH1D("h_pipluslam_McTruth_pz", "distribution of longitudinal momentum for #pi^{+}_{MC}(#bar{#Lambda}^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpPiPlus->Project("h_pipluslam_McTruth_pz", "PiPlus_MC_pz", "McTruthMatch==1 && Mother==-3122");

	TH1D * h_pipluslam_MC_tht = new TH1D("h_pipluslam_MC_tht", "#Theta_{MC} distribution for #pi^{+}(#bar{#Lambda}^{0}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_pipluslam_MC_tht", "cos(tht)", "moth==5 && pdg==211");

	TH1D * h_pipluslam_MC_pz = new TH1D("h_pipluslam_MC_pz", "distribution of longitudinal momentum (MC) for #pi^{+}(#bar{#Lambda}^{0}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_pipluslam_MC_pz", "pz", "moth==5 && pdg==211");

	TH1D * h_pipluslam_eff_pz = new TH1D("h_pipluslam_eff_pz", "reconstruction efficiency for #pi^{+} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_pipluslam_eff_pz->Divide(h_pipluslam_McTruth_pz, h_pipluslam_MC_pz);
	h_pipluslam_eff_pz->SetStats(0);
	h_pipluslam_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_pipluslam_eff_tht = new TH1D("h_pipluslam_eff_tht", "reconstruction efficiency for #pi^{+} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_pipluslam_eff_tht->Divide(h_pipluslam_McTruth_tht, h_pipluslam_MC_tht);
	h_pipluslam_eff_tht->SetStats(0);
	h_pipluslam_eff_tht->GetYaxis()->SetRangeUser(0,1.1);




	//save histograms
//	jenny::CreateDrawAndSaveHistogram(h_pipluslam_MC_tht, outPath, "h_pipluslam_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_pipluslam_pz, outPath, "h_pipluslam_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_pipluslam_npart, outPath, "h_pipluslam_npart", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_pipluslam_pt_vs_pz, outPath, "h_pipluslam_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_pipluslam_costht_vs_pz, outPath, "h_pipluslam_costht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_pipluslam_eff_tht, outPath, "h_pipluslam_eff_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_pipluslam_eff_pz, outPath, "h_pipluslam_eff_pz", save, close);


	//*********************** plots for comparison of pi+ and pi- *********************************************

//	jenny::CreateDrawAndSaveNHistograms(h_piminuslam_npart, h_pipluslam_npart, "#pi^{-}","#pi^{+}",  outPath, "h_pipluslam_piminuslam_npart", save, close);





	//************************ projections for lambda0 *******************************************************

	TH1D * h_lambda0_tht = new TH1D("h_lambda0_tht", "#theta distribution for #Lambda^{0}; #theta/rad; counts", 100,0,1);
	ntpLambda0->Project("h_lambda0_tht", "Lambda0_tht", "McTruthMatch==1");

	TH1D * h_lambda0_pz = new TH1D("h_lambda0_pz", "distribution of longitudinal momentum for #Lambda^{0}; p_{z}/GeV/c; counts", 100,0,2.7);
	ntpLambda0->Project("h_lambda0_pz", "Lambda0_pz", "McTruthMatch==1");


	TString cut1 = "VtxFit_prob>0.01";
	TString cut2 = cut1 + "&& fMass_Prob>0.01";


	//mass resolution with different cuts
	TH1D * h_lambda0_mass = new TH1D("h_lambda0_mass", "mass resolution of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpLambda0->Project("h_lambda0_mass", "(Lambda0_m-MCTruth_m)/1.115683", "McTruthMatch==1");

	TH1D * h_lambda0_mass_cut1 = new TH1D("h_lambda0_mass_cut1", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpLambda0->Project("h_lambda0_mass_cut1", "(Lambda0_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut1);

	TH1D * h_lambda0_mass_cut2 = new TH1D("h_lambda0_mass_cut2", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpLambda0->Project("h_lambda0_mass_cut2", "(Lambda0_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut2);


	//momentum resolution for different cuts

	TH1D * h_lambda0_p_diff = new TH1D("h_lambda0_p_diff", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpLambda0->Project("h_lambda0_p_diff", "VtxFit_momres", "McTruthMatch==1");

	TH1D * h_lambda0_p_cut1 = new TH1D("h_lambda0_p_cut1", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpLambda0->Project("h_lambda0_p_cut1", "VtxFit_momres", "McTruthMatch==1&&"+cut1);

	TH1D * h_lambda0_p_cut2 = new TH1D("h_lambda0_p_cut2", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpLambda0->Project("h_lambda0_p_cut2", "VtxFit_momres", "McTruthMatch==1&&"+cut2);




	//Vertex resolution

	TH1D * h_lambda0_vtx_resx = new TH1D("h_lambda0_vtx_resx", "vertex resolution for #Lambda^{0} (x coordinate); #Delta x[cm]; counts", 30,-3,3);
	ntpLambda0->Project("h_lambda0_vtx_resx", "VtxFit_diffvx", "McTruthMatch==1");

	TH1D * h_lambda0_vtx_resy = new TH1D("h_lambda0_vtx_resy", "vertex resolution for #Lambda^{0} (y coordinate); #Delta y[cm]; counts", 30,-3,3);
	ntpLambda0->Project("h_lambda0_vtx_resy", "VtxFit_diffvy", "McTruthMatch==1");

	TH1D * h_lambda0_vtx_resz = new TH1D("h_lambda0_vtx_resz", "vertex resolution for #Lambda^{0} (z coordinate); #Delta z[cm]; counts", 30,-3,3);
	ntpLambda0->Project("h_lambda0_vtx_resz", "VtxFit_diffvz", "McTruthMatch==1");

	TH1D * h_lambda0_vtx_pullx = new TH1D("h_lambda0_vtx_pullx", "pull distribution for decay vertex of #Lambda^{0} (x coordinate); #Delta x / #sigma_{x}; counts", 30, -3,3);
	ntpLambda0->Project("h_lambda0_vtx_pullx", "VtxFit_pullvx", "McTruthMatch==1");

	TH1D * h_lambda0_vtx_pully = new TH1D("h_lambda0_vtx_pully", "pull distribution for decay vertex of #Lambda^{0} (y coordinate); #Delta y / #sigma_{y}; counts", 30, -3,3);
	ntpLambda0->Project("h_lambda0_vtx_pully", "VtxFit_pullvy", "McTruthMatch==1");

	TH1D * h_lambda0_vtx_pullz = new TH1D("h_lambda0_vtx_pullz", "pull distribution for decay vertex of #Lambda^{0} (z coordinate); #Delta z / #sigma_{z}; counts", 30, -3,3);
	ntpLambda0->Project("h_lambda0_vtx_pullz", "VtxFit_pullvz", "McTruthMatch==1");

	//2d projections

	TH2D * h_lambda0_pt_vs_pz = new TH2D("h_lambda0_pt_vs_pz", "p_{t} vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; p_{t}/GeV/c", 100,0,2.7, 100,0,0.4);
	ntpLambda0->Project("h_lambda0_pt_vs_pz", "Lambda0_pt:Lambda0_pz", "McTruthMatch==1");

	TH2D * h_lambda0_tht_vs_pz = new TH2D("h_lambda0_tht_vs_pz", "cos(#theta) vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; cos(#theta)", 100,0,2.7, 100,-1,1);
	ntpLambda0->Project("h_lambda0_tht_vs_pz", "cos(Lambda0_tht):Lambda0_pz", "McTruthMatch==1");

	//reconstruction efficiency

	TH1D * h_lambda0_McTruth_tht = new TH1D("h_lambda0_McTruth_tht", "#Theta distribution for #Lambda^{0}_{MC}(#Xi^{-}); #theta/rad; counts" ,100,-1,1);
	ntpLambda0->Project("h_lambda0_McTruth_tht", "cos(lambda0_MC_tht)", "McTruthMatch==1 && Mother==3312");

	TH1D * h_lambda0_McTruth_pz = new TH1D("h_lambda0_McTruth_pz", "distribution of longitudinal momentum for #Lambda^{0}_{MC}(#Xi^{-}); pz/GeV/c; counts" ,100,0,2.7);
	ntpLambda0->Project("h_lambda0_McTruth_pz", "lambda0_MC_pz", "McTruthMatch==1 && Mother==3312");

	TH1D * h_lambda0_MC_tht = new TH1D("h_lambda0_MC_tht", "#Theta_{MC} distribution for #Lambda^{0}(#Xi^{-}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_lambda0_MC_tht", "cos(tht)", "moth==1 && pdg==3122");

	TH1D * h_lambda0_MC_pz = new TH1D("h_lambda0_MC_pz", "distribution of longitudinal momentum (MC) for #Lambda^{0}(#Xi^{-}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_lambda0_MC_pz", "pz", "moth==1 && pdg==3122");

	TH1D * h_lambda0_eff_pz = new TH1D("h_lambda0_eff_pz", "reconstruction efficiency for #Lambda^{0} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_lambda0_eff_pz->Divide(h_lambda0_McTruth_pz, h_lambda0_MC_pz);
	h_lambda0_eff_pz->SetStats(0);
	h_lambda0_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_lambda0_eff_tht = new TH1D("h_lambda0_eff_tht", "reconstruction efficiency for #Lambda^{0} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_lambda0_eff_tht->Divide(h_lambda0_McTruth_tht, h_lambda0_MC_tht);
	h_lambda0_eff_tht->SetStats(0);
	h_lambda0_eff_tht->GetYaxis()->SetRangeUser(0,1.1);


	// save histograms
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_tht, outPath, "h_lambda0_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_pz, outPath, "h_lambda0_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_mass, outPath, "h_lambda0_massres", save, close);
//	jenny::CreateDrawAndSaveNHistograms(h_lambda0_mass, h_lambda0_mass_cut1,h_lambda0_mass_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_lambda0_mass_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_p_diff2, outPath, "h_lambda0_p_diff", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_p_mc, outPath, "h_lambda0", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_pt_vs_pz, outPath, "h_lambda0_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_tht_vs_pz, outPath, "h_lambda0_tht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveNHistograms(h_lambda0_p_diff, h_lambda0_p_cut1,h_lambda0_p_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_lambda0_pres_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_vtx_resx, outPath, "h_lambda0_vtx_resx", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_vtx_resy, outPath, "h_lambda0_vtx_resy", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_vtx_resz, outPath, "h_lambda0_vtx_resz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_vtx_pullx, outPath, "h_lambda0_vtx_pullx", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_vtx_pully, outPath, "h_lambda0_vtx_pully", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_vtx_pullz, outPath, "h_lambda0_vtx_pullz", save, close);

//	jenny::CreateDrawAndSaveHistogram(h_lambda0_eff_pz, outPath, "h_lambda0_eff_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_lambda0_eff_tht , outPath, "h_lamda0_eff_tht", save, close);



	//************************ projections for antiLambda0 *******************************************************

	TH1D * h_antiLambda0_tht = new TH1D("h_antiLambda0_tht", "#theta distribution for #Lambda^{0}; #theta/rad; counts", 100,0,1);
	ntpAntiLambda0->Project("h_antiLambda0_tht", "antiLambda0_tht", "McTruthMatch==1");

	TH1D * h_antiLambda0_pz = new TH1D("h_antiLambda0_pz", "distribution of longitudinal momentum for #Lambda^{0}; p_{z}/GeV/c; counts", 100,0,2.7);
	ntpAntiLambda0->Project("h_antiLambda0_pz", "antiLambda0_pz", "McTruthMatch==1");


	TString cut1 = "VtxFit_prob>0.01";
	TString cut2 = cut1 + "&& fMass_Prob>0.01";


	//mass resolution with different cuts
	TH1D * h_antiLambda0_mass = new TH1D("h_antiLambda0_mass", "mass resolution of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpAntiLambda0->Project("h_antiLambda0_mass", "(antiLambda0_m-MCTruth_m)/1.115683", "McTruthMatch==1");

	TH1D * h_antiLambda0_mass_cut1 = new TH1D("h_antiLambda0_mass_cut1", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpAntiLambda0->Project("h_antiLambda0_mass_cut1", "(antiLambda0_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut1);

	TH1D * h_antiLambda0_mass_cut2 = new TH1D("h_antiLambda0_mass_cut2", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpAntiLambda0->Project("h_antiLambda0_mass_cut2", "(antiLambda0_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut2);


	//momentum resolution for different cuts

	TH1D * h_antiLambda0_p_diff = new TH1D("h_antiLambda0_p_diff", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpAntiLambda0->Project("h_antiLambda0_p_diff", "VtxFit_momres", "McTruthMatch==1");

	TH1D * h_antiLambda0_p_cut1 = new TH1D("h_antiLambda0_p_cut1", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpAntiLambda0->Project("h_antiLambda0_p_cut1", "VtxFit_momres", "McTruthMatch==1&&"+cut1);

	TH1D * h_antiLambda0_p_cut2 = new TH1D("h_antiLambda0_p_cut2", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpAntiLambda0->Project("h_antiLambda0_p_cut2", "VtxFit_momres", "McTruthMatch==1&&"+cut2);




	//Vertex resolution

	TH1D * h_antiLambda0_vtx_resx = new TH1D("h_antiLambda0_vtx_resx", "vertex resolution for #Lambda^{0} (x coordinate); #Delta x[cm]; counts", 30,-3,3);
	ntpAntiLambda0->Project("h_antiLambda0_vtx_resx", "VtxFit_diffvx", "McTruthMatch==1");

	TH1D * h_antiLambda0_vtx_resy = new TH1D("h_antiLambda0_vtx_resy", "vertex resolution for #Lambda^{0} (y coordinate); #Delta y[cm]; counts", 30,-3,3);
	ntpAntiLambda0->Project("h_antiLambda0_vtx_resy", "VtxFit_diffvy", "McTruthMatch==1");

	TH1D * h_antiLambda0_vtx_resz = new TH1D("h_antiLambda0_vtx_resz", "vertex resolution for #Lambda^{0} (z coordinate); #Delta z[cm]; counts", 30,-3,3);
	ntpAntiLambda0->Project("h_antiLambda0_vtx_resz", "VtxFit_diffvz", "McTruthMatch==1");

	TH1D * h_antiLambda0_vtx_pullx = new TH1D("h_antiLambda0_vtx_pullx", "pull distribution for decay vertex of #Lambda^{0} (x coordinate); #Delta x / #sigma_{x}; counts", 30, -3,3);
	ntpAntiLambda0->Project("h_antiLambda0_vtx_pullx", "VtxFit_pullvx", "McTruthMatch==1");

	TH1D * h_antiLambda0_vtx_pully = new TH1D("h_antiLambda0_vtx_pully", "pull distribution for decay vertex of #Lambda^{0} (y coordinate); #Delta y / #sigma_{y}; counts", 30, -3,3);
	ntpAntiLambda0->Project("h_antiLambda0_vtx_pully", "VtxFit_pullvy", "McTruthMatch==1");

	TH1D * h_antiLambda0_vtx_pullz = new TH1D("h_antiLambda0_vtx_pullz", "pull distribution for decay vertex of #Lambda^{0} (z coordinate); #Delta z / #sigma_{z}; counts", 30, -3,3);
	ntpAntiLambda0->Project("h_antiLambda0_vtx_pullz", "VtxFit_pullvz", "McTruthMatch==1");

	//2d projections

	TH2D * h_antiLambda0_pt_vs_pz = new TH2D("h_antiLambda0_pt_vs_pz", "p_{t} vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; p_{t}/GeV/c", 100,0,2.7, 100,0,0.4);
	ntpAntiLambda0->Project("h_antiLambda0_pt_vs_pz", "antiLambda0_pt:antiLambda0_pz", "McTruthMatch==1");

	TH2D * h_antiLambda0_tht_vs_pz = new TH2D("h_antiLambda0_tht_vs_pz", "cos(#theta) vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; cos(#theta)", 100,0,2.7, 100,-1,1);
	ntpAntiLambda0->Project("h_antiLambda0_tht_vs_pz", "cos(antiLambda0_tht):antiLambda0_pz", "McTruthMatch==1");

	//reconstruction efficiency

	TH1D * h_antiLambda0_McTruth_tht = new TH1D("h_antiLambda0_McTruth_tht", "#Theta distribution for #Lambda^{0}_{MC}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpAntiLambda0->Project("h_antiLambda0_McTruth_tht", "cos(antiLambda0_MC_tht)", "McTruthMatch==1 && Mother==3312");

	TH1D * h_antiLambda0_McTruth_pz = new TH1D("h_antiLambda0_McTruth_pz", "distribution of longitudinal momentum for #Lambda^{0}_{MC}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpAntiLambda0->Project("h_antiLambda0_McTruth_pz", "antiLambda0_MC_pz", "McTruthMatch==1 && Mother==3312");

	TH1D * h_antiLambda0_MC_tht = new TH1D("h_antiLambda0_MC_tht", "#Theta_{MC} distribution for #Lambda^{0}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_antiLambda0_MC_tht", "cos(tht)", "moth==2 && pdg==-3122");

	TH1D * h_antiLambda0_MC_pz = new TH1D("h_antiLambda0_MC_pz", "distribution of longitudinal momentum (MC) for #Lambda^{0}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_antiLambda0_MC_pz", "pz", "moth==2 && pdg==-3122");

	TH1D * h_antiLambda0_eff_pz = new TH1D("h_antiLambda0_eff_pz", "reconstruction efficiency for #Lambda^{0} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_antiLambda0_eff_pz->Divide(h_antiLambda0_McTruth_pz, h_antiLambda0_MC_pz);
	h_antiLambda0_eff_pz->SetStats(0);
	h_antiLambda0_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_antiLambda0_eff_tht = new TH1D("h_antiLambda0_eff_tht", "reconstruction efficiency for #Lambda^{0} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_antiLambda0_eff_tht->Divide(h_antiLambda0_McTruth_tht, h_antiLambda0_MC_tht);
	h_antiLambda0_eff_tht->SetStats(0);
	h_antiLambda0_eff_tht->GetYaxis()->SetRangeUser(0,1.1);


	// save histograms
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_tht, outPath, "h_antiLambda0_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_pz, outPath, "h_antiLambda0_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_mass, outPath, "h_antiLambda0_massres", save, close);
//	jenny::CreateDrawAndSaveNHistograms(h_antiLambda0_mass, h_antiLambda0_mass_cut1,h_antiLambda0_mass_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_antiLambda0_mass_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_p_diff2, outPath, "h_antiLambda0_p_diff", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_p_mc, outPath, "h_antiLambda0", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_pt_vs_pz, outPath, "h_antiLambda0_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_tht_vs_pz, outPath, "h_antiLambda0_tht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveNHistograms(h_antiLambda0_p_diff, h_antiLambda0_p_cut1,h_antiLambda0_p_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_antiLambda0_pres_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_vtx_resx, outPath, "h_antiLambda0_vtx_resx", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_vtx_resy, outPath, "h_antiLambda0_vtx_resy", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_vtx_resz, outPath, "h_antiLambda0_vtx_resz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_vtx_pullx, outPath, "h_antiLambda0_vtx_pullx", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_vtx_pully, outPath, "h_antiLambda0_vtx_pully", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_vtx_pullz, outPath, "h_antiLambda0_vtx_pullz", save, close);

//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_eff_pz, outPath, "h_antiLambda0_eff_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_antiLambda0_eff_tht , outPath, "h_lamda0_eff_tht", save, close);




	//************************ projections for Xi- *******************************************************

	TH1D * h_ximinus_tht = new TH1D("h_ximinus_tht", "#theta distribution for #Lambda^{0}; #theta/rad; counts", 100,0,1);
	ntpXiMinus->Project("h_ximinus_tht", "ximinus_tht", "McTruthMatch==1");

	TH1D * h_ximinus_pz = new TH1D("h_ximinus_pz", "distribution of longitudinal momentum for #Lambda^{0}; p_{z}/GeV/c; counts", 100,0,2.7);
	ntpXiMinus->Project("h_ximinus_pz", "ximinus_pz", "McTruthMatch==1");


	TString cut1 = "VtxFit_prob>0.01";
	TString cut2 = cut1 + "&& fMass_Prob>0.01";


	//mass resolution with different cuts
	TH1D * h_ximinus_mass = new TH1D("h_ximinus_mass", "mass resolution of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpXiMinus->Project("h_ximinus_mass", "(ximinus_m-MCTruth_m)/1.115683", "McTruthMatch==1");

	TH1D * h_ximinus_mass_cut1 = new TH1D("h_ximinus_mass_cut1", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpXiMinus->Project("h_ximinus_mass_cut1", "(ximinus_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut1);

	TH1D * h_ximinus_mass_cut2 = new TH1D("h_ximinus_mass_cut2", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpXiMinus->Project("h_ximinus_mass_cut2", "(ximinus_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut2);


	//momentum resolution for different cuts

	TH1D * h_ximinus_p_diff = new TH1D("h_ximinus_p_diff", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpXiMinus->Project("h_ximinus_p_diff", "VtxFit_momres", "McTruthMatch==1");

	TH1D * h_ximinus_p_cut1 = new TH1D("h_ximinus_p_cut1", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpXiMinus->Project("h_ximinus_p_cut1", "VtxFit_momres", "McTruthMatch==1&&"+cut1);

	TH1D * h_ximinus_p_cut2 = new TH1D("h_ximinus_p_cut2", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpXiMinus->Project("h_ximinus_p_cut2", "VtxFit_momres", "McTruthMatch==1&&"+cut2);




	//Vertex resolution

	TH1D * h_ximinus_vtx_resx = new TH1D("h_ximinus_vtx_resx", "vertex resolution for #Lambda^{0} (x coordinate); #Delta x[cm]; counts", 30,-3,3);
	ntpXiMinus->Project("h_ximinus_vtx_resx", "VtxFit_diffvx", "McTruthMatch==1");

	TH1D * h_ximinus_vtx_resy = new TH1D("h_ximinus_vtx_resy", "vertex resolution for #Lambda^{0} (y coordinate); #Delta y[cm]; counts", 30,-3,3);
	ntpXiMinus->Project("h_ximinus_vtx_resy", "VtxFit_diffvy", "McTruthMatch==1");

	TH1D * h_ximinus_vtx_resz = new TH1D("h_ximinus_vtx_resz", "vertex resolution for #Lambda^{0} (z coordinate); #Delta z[cm]; counts", 30,-3,3);
	ntpXiMinus->Project("h_ximinus_vtx_resz", "VtxFit_diffvz", "McTruthMatch==1");

	TH1D * h_ximinus_vtx_pullx = new TH1D("h_ximinus_vtx_pullx", "pull distribution for decay vertex of #Lambda^{0} (x coordinate); #Delta x / #sigma_{x}; counts", 30, -3,3);
	ntpXiMinus->Project("h_ximinus_vtx_pullx", "VtxFit_pullvx", "McTruthMatch==1");

	TH1D * h_ximinus_vtx_pully = new TH1D("h_ximinus_vtx_pully", "pull distribution for decay vertex of #Lambda^{0} (y coordinate); #Delta y / #sigma_{y}; counts", 30, -3,3);
	ntpXiMinus->Project("h_ximinus_vtx_pully", "VtxFit_pullvy", "McTruthMatch==1");

	TH1D * h_ximinus_vtx_pullz = new TH1D("h_ximinus_vtx_pullz", "pull distribution for decay vertex of #Lambda^{0} (z coordinate); #Delta z / #sigma_{z}; counts", 30, -3,3);
	ntpXiMinus->Project("h_ximinus_vtx_pullz", "VtxFit_pullvz", "McTruthMatch==1");

	//2d projections

	TH2D * h_ximinus_pt_vs_pz = new TH2D("h_ximinus_pt_vs_pz", "p_{t} vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; p_{t}/GeV/c", 100,0,2.7, 100,0,0.4);
	ntpXiMinus->Project("h_ximinus_pt_vs_pz", "ximinus_pt:ximinus_pz", "McTruthMatch==1");

	TH2D * h_ximinus_tht_vs_pz = new TH2D("h_ximinus_tht_vs_pz", "cos(#theta) vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; cos(#theta)", 100,0,2.7, 100,-1,1);
	ntpXiMinus->Project("h_ximinus_tht_vs_pz", "cos(ximinus_tht):ximinus_pz", "McTruthMatch==1");

	//reconstruction efficiency

	TH1D * h_ximinus_McTruth_tht = new TH1D("h_ximinus_McTruth_tht", "#Theta distribution for #Lambda^{0}_{MC}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpXiMinus->Project("h_ximinus_McTruth_tht", "cos(ximinus_MC_tht)", "McTruthMatch==1 && Mother==88888");

	TH1D * h_ximinus_McTruth_pz = new TH1D("h_ximinus_McTruth_pz", "distribution of longitudinal momentum for #Lambda^{0}_{MC}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpXiMinus->Project("h_ximinus_McTruth_pz", "ximinus_MC_pz", "McTruthMatch==1 && Mother==88888");

	TH1D * h_ximinus_MC_tht = new TH1D("h_ximinus_MC_tht", "#Theta_{MC} distribution for #Lambda^{0}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_ximinus_MC_tht", "cos(tht)", "moth==0 && pdg==3312");

	TH1D * h_ximinus_MC_pz = new TH1D("h_ximinus_MC_pz", "distribution of longitudinal momentum (MC) for #Lambda^{0}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_ximinus_MC_pz", "pz", "moth==0 && pdg==3312");

	TH1D * h_ximinus_eff_pz = new TH1D("h_ximinus_eff_pz", "reconstruction efficiency for #Lambda^{0} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_ximinus_eff_pz->Divide(h_ximinus_McTruth_pz, h_ximinus_MC_pz);
	h_ximinus_eff_pz->SetStats(0);
	h_ximinus_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_ximinus_eff_tht = new TH1D("h_ximinus_eff_tht", "reconstruction efficiency for #Lambda^{0} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_ximinus_eff_tht->Divide(h_ximinus_McTruth_tht, h_ximinus_MC_tht);
	h_ximinus_eff_tht->SetStats(0);
	h_ximinus_eff_tht->GetYaxis()->SetRangeUser(0,1.1);


	// save histograms
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_tht, outPath, "h_ximinus_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_pz, outPath, "h_ximinus_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_mass, outPath, "h_ximinus_massres", save, close);
//	jenny::CreateDrawAndSaveNHistograms(h_ximinus_mass, h_ximinus_mass_cut1,h_ximinus_mass_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_ximinus_mass_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_p_diff2, outPath, "h_ximinus_p_diff", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_p_mc, outPath, "h_ximinus", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_pt_vs_pz, outPath, "h_ximinus_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_tht_vs_pz, outPath, "h_ximinus_tht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveNHistograms(h_ximinus_p_diff, h_ximinus_p_cut1,h_ximinus_p_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_ximinus_pres_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_vtx_resx, outPath, "h_ximinus_vtx_resx", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_vtx_resy, outPath, "h_ximinus_vtx_resy", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_vtx_resz, outPath, "h_ximinus_vtx_resz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_vtx_pullx, outPath, "h_ximinus_vtx_pullx", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_vtx_pully, outPath, "h_ximinus_vtx_pully", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_vtx_pullz, outPath, "h_ximinus_vtx_pullz", save, close);

//	jenny::CreateDrawAndSaveHistogram(h_ximinus_eff_pz, outPath, "h_ximinus_eff_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_ximinus_eff_tht , outPath, "h_lamda0_eff_tht", save, close);





	//************************ projections for Xi+ *******************************************************

	TH1D * h_xiplus_tht = new TH1D("h_xiplus_tht", "#theta distribution for #Lambda^{0}; #theta/rad; counts", 100,0,1);
	ntpXiPlus->Project("h_xiplus_tht", "xiplus_tht", "McTruthMatch==1");

	TH1D * h_xiplus_pz = new TH1D("h_xiplus_pz", "distribution of longitudinal momentum for #Lambda^{0}; p_{z}/GeV/c; counts", 100,0,2.7);
	ntpXiPlus->Project("h_xiplus_pz", "xiplus_pz", "McTruthMatch==1");



	//mass resolution with different cuts
	TH1D * h_xiplus_mass = new TH1D("h_xiplus_mass", "mass resolution of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpXiPlus->Project("h_xiplus_mass", "(xiplus_m-MCTruth_m)/1.115683", "McTruthMatch==1");

	TH1D * h_xiplus_mass_cut1 = new TH1D("h_xiplus_mass_cut1", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpXiPlus->Project("h_xiplus_mass_cut1", "(xiplus_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut1);

	TH1D * h_xiplus_mass_cut2 = new TH1D("h_xiplus_mass_cut2", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpXiPlus->Project("h_xiplus_mass_cut2", "(xiplus_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut2);


	//momentum resolution for different cuts

	TH1D * h_xiplus_p_diff = new TH1D("h_xiplus_p_diff", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpXiPlus->Project("h_xiplus_p_diff", "VtxFit_momres", "McTruthMatch==1");

	TH1D * h_xiplus_p_cut1 = new TH1D("h_xiplus_p_cut1", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpXiPlus->Project("h_xiplus_p_cut1", "VtxFit_momres", "McTruthMatch==1&&"+cut1);

	TH1D * h_xiplus_p_cut2 = new TH1D("h_xiplus_p_cut2", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpXiPlus->Project("h_xiplus_p_cut2", "VtxFit_momres", "McTruthMatch==1&&"+cut2);




	//Vertex resolution

	TH1D * h_xiplus_vtx_resx = new TH1D("h_xiplus_vtx_resx", "vertex resolution for #Lambda^{0} (x coordinate); #Delta x[cm]; counts", 30,-3,3);
	ntpXiPlus->Project("h_xiplus_vtx_resx", "VtxFit_diffvx", "McTruthMatch==1");

	TH1D * h_xiplus_vtx_resy = new TH1D("h_xiplus_vtx_resy", "vertex resolution for #Lambda^{0} (y coordinate); #Delta y[cm]; counts", 30,-3,3);
	ntpXiPlus->Project("h_xiplus_vtx_resy", "VtxFit_diffvy", "McTruthMatch==1");

	TH1D * h_xiplus_vtx_resz = new TH1D("h_xiplus_vtx_resz", "vertex resolution for #Lambda^{0} (z coordinate); #Delta z[cm]; counts", 30,-3,3);
	ntpXiPlus->Project("h_xiplus_vtx_resz", "VtxFit_diffvz", "McTruthMatch==1");

	TH1D * h_xiplus_vtx_pullx = new TH1D("h_xiplus_vtx_pullx", "pull distribution for decay vertex of #Lambda^{0} (x coordinate); #Delta x / #sigma_{x}; counts", 30, -3,3);
	ntpXiPlus->Project("h_xiplus_vtx_pullx", "VtxFit_pullvx", "McTruthMatch==1");

	TH1D * h_xiplus_vtx_pully = new TH1D("h_xiplus_vtx_pully", "pull distribution for decay vertex of #Lambda^{0} (y coordinate); #Delta y / #sigma_{y}; counts", 30, -3,3);
	ntpXiPlus->Project("h_xiplus_vtx_pully", "VtxFit_pullvy", "McTruthMatch==1");

	TH1D * h_xiplus_vtx_pullz = new TH1D("h_xiplus_vtx_pullz", "pull distribution for decay vertex of #Lambda^{0} (z coordinate); #Delta z / #sigma_{z}; counts", 30, -3,3);
	ntpXiPlus->Project("h_xiplus_vtx_pullz", "VtxFit_pullvz", "McTruthMatch==1");

	//2d projections

	TH2D * h_xiplus_pt_vs_pz = new TH2D("h_xiplus_pt_vs_pz", "p_{t} vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; p_{t}/GeV/c", 100,0,2.7, 100,0,0.4);
	ntpXiPlus->Project("h_xiplus_pt_vs_pz", "xiplus_pt:xiplus_pz", "McTruthMatch==1");

	TH2D * h_xiplus_tht_vs_pz = new TH2D("h_xiplus_tht_vs_pz", "cos(#theta) vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; cos(#theta)", 100,0,2.7, 100,-1,1);
	ntpXiPlus->Project("h_xiplus_tht_vs_pz", "cos(xiplus_tht):xiplus_pz", "McTruthMatch==1");

	//reconstruction efficiency

	TH1D * h_xiplus_McTruth_tht = new TH1D("h_xiplus_McTruth_tht", "#Theta distribution for #Lambda^{0}_{MC}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpXiPlus->Project("h_xiplus_McTruth_tht", "cos(xiplus_MC_tht)", "McTruthMatch==1 && Mother==88888");

	TH1D * h_xiplus_McTruth_pz = new TH1D("h_xiplus_McTruth_pz", "distribution of longitudinal momentum for #Lambda^{0}_{MC}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpXiPlus->Project("h_xiplus_McTruth_pz", "xiplus_MC_pz", "McTruthMatch==1 && Mother==88888");

	TH1D * h_xiplus_MC_tht = new TH1D("h_xiplus_MC_tht", "#Theta_{MC} distribution for #Lambda^{0}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_xiplus_MC_tht", "cos(tht)", "moth==0 && pdg==-3312");

	TH1D * h_xiplus_MC_pz = new TH1D("h_xiplus_MC_pz", "distribution of longitudinal momentum (MC) for #Lambda^{0}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_xiplus_MC_pz", "pz", "moth==0 && pdg==-3312");

	TH1D * h_xiplus_eff_pz = new TH1D("h_xiplus_eff_pz", "reconstruction efficiency for #Lambda^{0} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_xiplus_eff_pz->Divide(h_xiplus_McTruth_pz, h_xiplus_MC_pz);
	h_xiplus_eff_pz->SetStats(0);
	h_xiplus_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_xiplus_eff_tht = new TH1D("h_xiplus_eff_tht", "reconstruction efficiency for #Lambda^{0} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_xiplus_eff_tht->Divide(h_xiplus_McTruth_tht, h_xiplus_MC_tht);
	h_xiplus_eff_tht->SetStats(0);
	h_xiplus_eff_tht->GetYaxis()->SetRangeUser(0,1.1);


	// save histograms
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_tht, outPath, "h_xiplus_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_pz, outPath, "h_xiplus_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_mass, outPath, "h_xiplus_massres", save, close);
//	jenny::CreateDrawAndSaveNHistograms(h_xiplus_mass, h_xiplus_mass_cut1,h_xiplus_mass_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_xiplus_mass_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_p_diff2, outPath, "h_xiplus_p_diff", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_p_mc, outPath, "h_xiplus", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_pt_vs_pz, outPath, "h_xiplus_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_tht_vs_pz, outPath, "h_xiplus_tht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveNHistograms(h_xiplus_p_diff, h_xiplus_p_cut1,h_xiplus_p_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_xiplus_pres_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_vtx_resx, outPath, "h_xiplus_vtx_resx", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_vtx_resy, outPath, "h_xiplus_vtx_resy", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_vtx_resz, outPath, "h_xiplus_vtx_resz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_vtx_pullx, outPath, "h_xiplus_vtx_pullx", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_vtx_pully, outPath, "h_xiplus_vtx_pully", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_vtx_pullz, outPath, "h_xiplus_vtx_pullz", save, close);

//	jenny::CreateDrawAndSaveHistogram(h_xiplus_eff_pz, outPath, "h_xiplus_eff_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xiplus_eff_tht , outPath, "h_lamda0_eff_tht", save, close);



	//************************ projections for Xi-Xi+ system *******************************************************

	TH1D * h_xisys_tht = new TH1D("h_xisys_tht", "#theta distribution for #Lambda^{0}; #theta/rad; counts", 100,0,1);
	ntpXiSys->Project("h_xisys_tht", "xisys_tht", "McTruthMatch==1");

	TH1D * h_xisys_pz = new TH1D("h_xisys_pz", "distribution of longitudinal momentum for #Lambda^{0}; p_{z}/GeV/c; counts", 100,0,2.7);
	ntpXiSys->Project("h_xisys_pz", "xisys_pz", "McTruthMatch==1");



	//mass resolution with different cuts
	TH1D * h_xisys_mass = new TH1D("h_xisys_mass", "mass resolution of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpXiSys->Project("h_xisys_mass", "(xisys_m-MCTruth_m)/1.115683", "McTruthMatch==1");

	TH1D * h_xisys_mass_cut1 = new TH1D("h_xisys_mass_cut1", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpXiSys->Project("h_xisys_mass_cut1", "(xisys_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut1);

	TH1D * h_xisys_mass_cut2 = new TH1D("h_xisys_mass_cut2", "mass of #Lambda^{0}; m/Gev/c^{2}; counts", 100,-0.15,0.15);
	ntpXiSys->Project("h_xisys_mass_cut2", "(xisys_m-MCTruth_m)/1.115683", "McTruthMatch==1 &&"+cut2);


	//momentum resolution for different cuts

	TH1D * h_xisys_p_diff = new TH1D("h_xisys_p_diff", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpXiSys->Project("h_xisys_p_diff", "4CFit_momres", "McTruthMatch==1");

	TH1D * h_xisys_p_cut1 = new TH1D("h_xisys_p_cut1", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpXiSys->Project("h_xisys_p_cut1", "4CFit_momres", "McTruthMatch==1&&"+cut1);

	TH1D * h_xisys_p_cut2 = new TH1D("h_xisys_p_cut2", "momentum of #Lambda^{0}; p/GeV/c; counts", 100,-1,1);
	ntpXiSys->Project("h_xisys_p_cut2", "4CFit_momres", "McTruthMatch==1&&"+cut2);



	//Vertex resolution 4CFit

	TH1D * h_xisys_vtx_resx = new TH1D("h_xisys_vtx_resx", "vertex resolution for #Lambda^{0} (x coordinate); #Delta x[cm]; counts", 30,-3,3);
	ntpXiSys->Project("h_xisys_vtx_resx", "4CFit_diffvx", "McTruthMatch==1");

	TH1D * h_xisys_vtx_resy = new TH1D("h_xisys_vtx_resy", "vertex resolution for #Lambda^{0} (y coordinate); #Delta y[cm]; counts", 30,-3,3);
	ntpXiSys->Project("h_xisys_vtx_resy", "4CFit_diffvy", "McTruthMatch==1");

	TH1D * h_xisys_vtx_resz = new TH1D("h_xisys_vtx_resz", "vertex resolution for #Lambda^{0} (z coordinate); #Delta z[cm]; counts", 30,-3,3);
	ntpXiSys->Project("h_xisys_vtx_resz", "4CFit_diffvz", "McTruthMatch==1");

	TH1D * h_xisys_vtx_pullx = new TH1D("h_xisys_vtx_pullx", "pull distribution for decay vertex of #Lambda^{0} (x coordinate); #Delta x / #sigma_{x}; counts", 30, -3,3);
	ntpXiSys->Project("h_xisys_vtx_pullx", "4CFit_pullvx", "McTruthMatch==1");

	TH1D * h_xisys_vtx_pully = new TH1D("h_xisys_vtx_pully", "pull distribution for decay vertex of #Lambda^{0} (y coordinate); #Delta y / #sigma_{y}; counts", 30, -3,3);
	ntpXiSys->Project("h_xisys_vtx_pully", "4CFit_pullvy", "McTruthMatch==1");

	TH1D * h_xisys_vtx_pullz = new TH1D("h_xisys_vtx_pullz", "pull distribution for decay vertex of #Lambda^{0} (z coordinate); #Delta z / #sigma_{z}; counts", 30, -3,3);
	ntpXiSys->Project("h_xisys_vtx_pullz", "4CFit_pullvz", "McTruthMatch==1");

	//2d projections

	TH2D * h_xisys_pt_vs_pz = new TH2D("h_xisys_pt_vs_pz", "p_{t} vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; p_{t}/GeV/c", 100,0,2.7, 100,0,0.4);
	ntpXiSys->Project("h_xisys_pt_vs_pz", "xisys_pt:xisys_pz", "McTruthMatch==1");

	TH2D * h_xisys_tht_vs_pz = new TH2D("h_xisys_tht_vs_pz", "cos(#theta) vs. p_{z} for #Lambda^{0}; p_{z}/GeV/c; cos(#theta)", 100,0,2.7, 100,-1,1);
	ntpXiSys->Project("h_xisys_tht_vs_pz", "cos(xisys_tht):xisys_pz", "McTruthMatch==1");

	//reconstruction efficiency

	TH1D * h_xisys_McTruth_tht = new TH1D("h_xisys_McTruth_tht", "#Theta distribution for #Lambda^{0}_{MC}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpXiSys->Project("h_xisys_McTruth_tht", "cos(xisys_MC_tht)", "McTruthMatch==1 && Mother==88888");

	TH1D * h_xisys_McTruth_pz = new TH1D("h_xisys_McTruth_pz", "distribution of longitudinal momentum for #Lambda^{0}_{MC}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpXiSys->Project("h_xisys_McTruth_pz", "xisys_MC_pz", "McTruthMatch==1 && Mother==88888");

	TH1D * h_xisys_MC_tht = new TH1D("h_xisys_MC_tht", "#Theta_{MC} distribution for #Lambda^{0}(#Xi^{+}); #theta/rad; counts" ,100,-1,1);
	ntpMC->Project("h_xisys_MC_tht", "cos(tht)", "moth==-1 && pdg==88888");

	TH1D * h_xisys_MC_pz = new TH1D("h_xisys_MC_pz", "distribution of longitudinal momentum (MC) for #Lambda^{0}(#Xi^{+}); pz/GeV/c; counts" ,100,0,2.7);
	ntpMC->Project("h_xisys_MC_pz", "pz", "moth==-1 && pdg==88888");

	TH1D * h_xisys_eff_pz = new TH1D("h_xisys_eff_pz", "reconstruction efficiency for #Lambda^{0} depending on p_{z}; p_{z}/GeV/c; efficiency", 100,0,2.7);
	h_xisys_eff_pz->Divide(h_xisys_McTruth_pz, h_xisys_MC_pz);
	h_xisys_eff_pz->SetStats(0);
	h_xisys_eff_pz->GetYaxis()->SetRangeUser(0,1.01);

	TH1D * h_xisys_eff_tht = new TH1D("h_xisys_eff_tht", "reconstruction efficiency for #Lambda^{0} depending on cos(#theta); cos(#theta); efficiency", 100,-1,1);
	h_xisys_eff_tht->Divide(h_xisys_McTruth_tht, h_xisys_MC_tht);
	h_xisys_eff_tht->SetStats(0);
	h_xisys_eff_tht->GetYaxis()->SetRangeUser(0,1.1);


	// save histograms
//	jenny::CreateDrawAndSaveHistogram(h_xisys_tht, outPath, "h_xisys_tht", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_pz, outPath, "h_xisys_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_mass, outPath, "h_xisys_massres", save, close);
//	jenny::CreateDrawAndSaveNHistograms(h_xisys_mass, h_xisys_mass_cut1,h_xisys_mass_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_xisys_mass_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_p_diff2, outPath, "h_xisys_p_diff", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_p_mc, outPath, "h_xisys", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_pt_vs_pz, outPath, "h_xisys_pt_vs_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_tht_vs_pz, outPath, "h_xisys_tht_vs_pz", save, close);
//	jenny::CreateDrawAndSaveNHistograms(h_xisys_p_diff, h_xisys_p_cut1,h_xisys_p_cut2, "no cut", "VtxFit_prob>0.01", "cut on MassFit", outPath, "h_xisys_pres_cuts", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_vtx_resx, outPath, "h_xisys_vtx_resx", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_vtx_resy, outPath, "h_xisys_vtx_resy", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_vtx_resz, outPath, "h_xisys_vtx_resz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_vtx_pullx, outPath, "h_xisys_vtx_pullx", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_vtx_pully, outPath, "h_xisys_vtx_pully", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_vtx_pullz, outPath, "h_xisys_vtx_pullz", save, close);

//	jenny::CreateDrawAndSaveHistogram(h_xisys_eff_pz, outPath, "h_xisys_eff_pz", save, close);
//	jenny::CreateDrawAndSaveHistogram(h_xisys_eff_tht , outPath, "h_lamda0_eff_tht", save, close);


}
