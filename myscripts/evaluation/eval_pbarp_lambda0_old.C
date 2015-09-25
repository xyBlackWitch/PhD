class RhoTuple;

#include <vector>
#include <iostream>

#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TString.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"


void eval_pbarp_lambda0(TString inputDir="/private/puetz/fairsoft_mar15/pandaroot/mysimulations/analysis/pbarp_lambda0_antilambda0/10000_events/", bool saveOutput = true){

  TString filename_ana = inputDir + "output_ana.root";
  
  //***open input
  TFile * data = new TFile(filename_ana, "READ");
  TTree * ntpPiPlus= (TTree*)data->Get("ntpPiPlus");
  TTree * ntpPiMinus = (TTree*)data->Get("ntpPiMinus");
  TTree * ntpMc = (TTree*)data->Get("ntpMC");
  TTree * ntpProton = (TTree*)data->Get("ntpProton");
  TTree * ntpAntiProton = (TTree*)data->Get("ntpAntiProton");
	TTree * ntpLambda0 = (TTree*)data->Get("ntpLambda0");



  //projection for MC

  TH2F * h_theta_vs_p_mc = new TH2F("h_theta_vs_p_mc", "#Theta vs. P for MC; P / GeV/c ; #Theta", 100,0.3,1.4,100,0,3);
  ntpMC->Project("h_theta_vs_p_mc", "tht:p");//, "moth==0&&pdg==3122" );

  TH2F * h_costheta_vs_p_mc = new TH2F("h_cosTheta_vs_p_mc", "cos(#Theta) vs. P_{L} for MC; P_{L} / GeV/c ; cos(#Theta)", 100,-1,1.8,100,-1,1);
  ntpMC->Project("h_cosTheta_vs_p_mc", "cos(tht):pz");//, "moth==0&&pdg==3122");
  
  TH2F * h_pt_vs_pl_mc = new TH2F("h_pt_vs_pl_mc", "P_{t} vs.P_{l} for MC; P_{l}/GeV/c; P_{t}/GeV/c", 100,-1,1.8,100,0,0.4);
  ntpMC->Project("h_pt_vs_pl_mc", "sqrt(px*px+py*py):pz");//, "moth==0&&pdg==3122");

  

	TH2F * h_theta_vs_p_mc_antiProton = new TH2F("h_theta_vs_p_mc_antiProton", "#Theta vs. P for AntiProton  MC; P / GeV/c ; #Theta", 100,0.3,1.4,100,0,3);
  ntpMC->Project("h_theta_vs_p_mc_antiProton", "tht:p", "moth==2&&pdg==-2212" );

  TH2F * h_costheta_vs_p_mc_antiProton = new TH2F("h_cosTheta_vs_p_mc_antiProton", "cos(#Theta) vs. P_{L} for AntiProtonMC; P_{L} / GeV/c ; cos(#Theta)", 100,-1,1.8,100,-1,1);
  ntpMC->Project("h_cosTheta_vs_p_mc_antiProton", "cos(tht):pz", "moth==2&&pdg==-2212");
  
  TH2F * h_pt_vs_pl_mc_antiProton = new TH2F("h_pt_vs_pl_mc_antiProton", "P_{t} vs.P_{l} for  MC AntiProton; P_{l}/GeV/c; P_{t}/GeV/c", 100,-1,1.8,100,0,0.4);
  ntpMC->Project("h_pt_vs_pl_mc_antiProton", "sqrt(px*px+py*py):pz", "moth==2&&pdg==-2212");

	
	
	TH2F * h_theta_vs_p_mc_proton = new TH2F("h_theta_vs_p_mc_proton", "#Theta vs. P for  MC Proton; P / GeV/c ; #Theta", 100,0.3,1.4,100,0,3);
  ntpMC->Project("h_theta_vs_p_mc_proton", "tht:p", "moth==1&&pdg==2212" );

  TH2F * h_costheta_vs_p_mc_proton = new TH2F("h_cosTheta_vs_p_mc_proton", "cos(#Theta) vs. P_{L} for  MC Proton; P_{L} / GeV/c ; cos(#Theta)", 100,-1,1.8,100,-1,1);
  ntpMC->Project("h_cosTheta_vs_p_mc_proton", "cos(tht):pz", "moth==1&&pdg==2212");
  
  TH2F * h_pt_vs_pl_mc_proton = new TH2F("h_pt_vs_pl_mc_proton", "P_{t} vs.P_{l} for MC Proton; P_{l}/GeV/c; P_{t}/GeV/c", 100,-1,1.8,100,0,0.4);
  ntpMC->Project("h_pt_vs_pl_mc_proton", "sqrt(px*px+py*py):pz", "moth==1&&pdg==2212");

	

	TH2F * h_theta_vs_p_mc_piminus = new TH2F("h_theta_vs_p_mc_piminus", "#Theta vs. P for MC #pi^{-}; P / GeV/c ; #Theta", 100,0.3,1.4,100,0,3);
  ntpMC->Project("h_theta_vs_p_mc_piminus", "tht:p", "moth==1&&pdg==-211" );

  TH2F * h_costheta_vs_p_mc_piminus = new TH2F("h_cosTheta_vs_p_mc_piminus", "cos(#Theta) vs. P_{L} for MC #pi^{-}; P_{L} / GeV/c ; cos(#Theta)", 100,-1,1.8,100,-1,1);
  ntpMC->Project("h_cosTheta_vs_p_mc_piminus", "cos(tht):pz", "moth==1&&pdg==-211");
  
  TH2F * h_pt_vs_pl_mc_piminus = new TH2F("h_pt_vs_pl_mc_piminus", "P_{t} vs.P_{l} for MC #pi^{-}; P_{l}/GeV/c; P_{t}/GeV/c", 100,-1,1.8,100,0,0.4);
  ntpMC->Project("h_pt_vs_pl_mc_piminus", "sqrt(px*px+py*py):pz", "moth==1&&pdg==-211");

	

	TH2F * h_theta_vs_p_mc_piplus = new TH2F("h_theta_vs_p_mc_piplus", "#Theta vs. P for MC #pi^{+}; P / GeV/c ; #Theta", 100,0.3,1.4,100,0,3);
  ntpMC->Project("h_theta_vs_p_mc_piplus", "tht:p", "moth==2&&pdg==211" );

  TH2F * h_costheta_vs_p_mc_piplus = new TH2F("h_cosTheta_vs_p_mc_piplus", "cos(#Theta) vs. P_{L} for MC #pi^{+}; P_{L} / GeV/c ; cos(#Theta)", 100,-1,1.8,100,-1,1);
  ntpMC->Project("h_cosTheta_vs_p_mc_piplus", "cos(tht):pz", "moth==2&&pdg==211");
  
  TH2F * h_pt_vs_pl_mc_piplus = new TH2F("h_pt_vs_pl_mc_piplus", "P_{t} vs.P_{l} for MC #pi^{+}; P_{l}/GeV/c; P_{t}/GeV/c", 100,-1,1.8,100,0,0.4);
  ntpMC->Project("h_pt_vs_pl_mc_piplus", "sqrt(px*px+py*py):pz", "moth==2&&pdg==211");



	TH2F * h_theta_vs_p_mc_lambda0 = new TH2F("h_theta_vs_p_mc_lambda0", "#Theta vs. P for MC #Lambda^{0}; P / GeV/c ; #Theta", 100,0.3,1.4,100,0,3);
  ntpMC->Project("h_theta_vs_p_mc_lambda0", "tht:p", "moth==0&&pdg==3122" );

  TH2F * h_costheta_vs_p_mc_lambda0 = new TH2F("h_cosTheta_vs_p_mc_lambda0", "cos(#Theta) vs. P_{L} for MC #Lambda^{0}; P_{L} / GeV/c ; cos(#Theta)", 100,-1,1.8,100,-1,1);
  ntpMC->Project("h_cosTheta_vs_p_mc_lambda0", "cos(tht):pz", "moth==0&&pdg==3122");
  
  TH2F * h_pt_vs_pl_mc_lambda0 = new TH2F("h_pt_vs_pl_mc_lambda0", "P_{t} vs.P_{l} for MC #Lambda^{0}; P_{l}/GeV/c; P_{t}/GeV/c", 100,0,1.8,100,0,0.4);
  ntpMC->Project("h_pt_vs_pl_mc_lambda0", "sqrt(px*px+py*py):pz", "moth==0&&pdg==3122");



	TH2F * h_theta_vs_p_mc_antiLambda0 = new TH2F("h_theta_vs_p_mc_antiLambda0", "#Theta vs. P for MC #bar{#Lambda}^{0}; P / GeV/c ; #Theta", 100,0.3,1.4,100,0,3);
  ntpMC->Project("h_theta_vs_p_mc_antiLambda0", "tht:p", "moth==0&&pdg==-3122" );

  TH2F * h_costheta_vs_p_mc_antiLambda0 = new TH2F("h_cosTheta_vs_p_mc_antiLambda0", "cos(#Theta) vs. P_{L} for MC #bar{#Lambda}^{0}; P_{L} / GeV/c ; cos(#Theta)", 100,-1,1.8,100,-1,1);
  ntpMC->Project("h_cosTheta_vs_p_mc_antiLambda0", "cos(tht):pz", "moth==0&&pdg==-3122");
  
  TH2F * h_pt_vs_pl_mc_antiLambda0 = new TH2F("h_pt_vs_pl_mc_antiLambda0", "P_{t} vs.P_{l} for MC #bar{#Lambda}^{0}; P_{l}/GeV/c; P_{t}/GeV/c", 100,0,1.8,100,0,0.4);
  ntpMC->Project("h_pt_vs_pl_mc_antiLambda0", "sqrt(px*px+py*py):pz", "moth==0&&pdg==-3122");

 

	 //projections for PiPlus
	
  TH2F * h_pt_vs_pl_piplus = new TH2F("h_pt_vs_pl_piplus", "P_{t} vs. P_{l} for #pi^{+}; P_{l}; P_{t}", 100,-0.05,0.4,100,0,0.2);
  ntpPiPlus->Project("h_pt_vs_pl_piplus", "PiPlus_pt:PiPlus_pz", "MCTruthMatch==1");

  TH2F * h_cosTheta_vs_p_piplus = new TH2F("h_cosTheta_vs_p_piplus", "cos(#Theta) vs. P_{L} for #pi^{+}; P_{L} / GeV/c ; cos(#Theta)", 100,-0.05,0.4,100,-1.1,1.1);
  ntpPiPlus->Project("h_cosTheta_vs_p_piplus",  "PiPlus_CosTheta:PiPlus_pz", "MCTruthMatch==1");

  TH2F * h_theta_vs_p_piplus = new TH2F("h_theta_vs_p_piplus", "#Theta vs. P for #pi^{+}; P / GeV/c ; #Theta", 100,0,0.4,100,0,2);
	ntpPiPlus->Project("h_theta_vs_p_piplus", "PiPlus_tht:PiPlus_p", "MCTruthMatch==1");

  TH2F * h_pt_vs_pl_piplus_McTruth = new TH2F("h_pt_vs_pl_piplus_McTruth", "P_{t} vs. P_{l} for #pi^{+} (McTruth); P_{l}; P_{t}", 100,-0.05,0.4,100,0,0.2);
  ntpPiPlus->Project("h_pt_vs_pl_piplus_McTruth", "PiPlus_MC_pt:PiPlus_MC_pz", "MCTruthMatch==1&&Mother==-3122");

  TH2F * h_cosTheta_vs_p_piplus_McTruth = new TH2F("h_cosTheta_vs_p_piplus_McTruth", "cos(#Theta) vs. P_{L} for #pi^{+}; P_{L} / GeV/c ; cos(#Theta)", 100,-0.05,0.4,100,-1.1,1.1);
  ntpPiPlus->Project("h_cosTheta_vs_p_piplus_McTruth",  "PiPlus_MC_CosTheta:PiPlus_MC_pz", "MCTruthMatch==1&&Mother==-3122");

  TH2F * h_theta_vs_p_piplus_McTruth = new TH2F("h_theta_vs_p_piplus_McTruth", "#Theta vs. P for #pi^{+}; P / GeV/c ; #Theta", 100,0,0.4,100,0,2);
  ntpPiPlus->Project("h_theta_vs_p_piplus_McTruth", "PiPlus_MC_tht:PiPlus_MC_p", "MCTruthMatch==1&&Mother==-3122");

	//reco probability for PiPlus
  TH1F * h_reco_p_piplus = new TH1F("h_reco_p_piplus", "reconstructed momentum for #pi^{+}; P/GeV/c; counts", 100, 0,0.4);
	ntpPiPlus->Project("h_reco_p_piplus", "PiPlus_MC_p", "MCTruthMatch==1&&Mother==-3122");

	TH1F * h_mc_p_piplus = new TH1F("h_mc_p_piplus", "MC momentum for #pi^{+}; P/GeV/c; counts",100,0,0.4);
	ntpMC->Project("h_mc_p_piplus", "p", "moth==2&&pdg==211");

	TH1F * h_prob_vs_p_piplus = new TH1F("h_prob_vs_p_piplus","reconstruction probability for #pi^{+};P/GeV/c;reco. probability",100,0,0.4);
	h_prob_vs_p_piplus->Divide(h_reco_p_piplus,h_mc_p_piplus);
	h_prob_vs_p_piplus->SetStats(0);
	h_prob_vs_p_piplus->GetYaxis()->SetRangeUser(0,1.1);

  TH1F * h_reco_theta_piplus = new TH1F("h_reco_theta_piplus", "reconstructed angle for #pi^{+}; #theta/rad; counts", 100,0,2);
	ntpPiPlus->Project("h_reco_theta_piplus", "PiPlus_MC_tht", "MCTruthMatch==1&&Mother==-3122");

	TH1F * h_mc_theta_piplus = new TH1F("h_mc_theta_piplus", "MC angle for #pi^{+};#theta/rad; counts",100,0,2);
	ntpMC->Project("h_mc_theta_piplus", "tht", "moth==2&&pdg==211");

	TH1F * h_prob_vs_theta_piplus = new TH1F("h_prob_vs_theta_piplus","reconstruction probability for #pi^{+};#theta/rad;reco. probability",100,0,2);
	h_prob_vs_theta_piplus->Divide(h_reco_theta_piplus,h_mc_theta_piplus);
	h_prob_vs_theta_piplus->SetStats(0);
	h_prob_vs_theta_piplus->GetYaxis()->SetRangeUser(0,1.1);
	
  //projections for PiMinus

  TH2F * h_pt_vs_pl_piminus = new TH2F("h_pt_vs_pl_piminus", "P_{t} vs. P_{l} for #pi^{-}; P_{l}; P_{t}", 100,-0.05,0.4,100,0,0.2);
  ntpPiMinus->Project("h_pt_vs_pl_piminus", "PiMinus_pt:PiMinus_pz", "MCTruthMatch==1");

  TH2F * h_cosTheta_vs_p_piminus = new TH2F("h_cosTheta_vs_p_piminus", "cos(#Theta) vs. P_{L} for #pi^{-}; P_{L} / GeV/c ; cos(#Theta)", 100,-0.05,0.4,100,-1.1,1.1);
  ntpPiMinus->Project("h_cosTheta_vs_p_piminus",  "PiMinus_CosTheta:PiMinus_pz", "MCTruthMatch==1");

  TH2F * h_theta_vs_p_piminus = new TH2F("h_theta_vs_p_piminus", "#Theta vs. P for #pi^{-}; P / GeV/c ; #Theta", 100,0,0.4,100,0,2);
	ntpPiMinus->Project("h_theta_vs_p_piminus", "PiMinus_tht:PiMinus_p", "MCTruthMatch==1");

  TH2F * h_pt_vs_pl_piminus_McTruth = new TH2F("h_pt_vs_pl_piminus_McTruth", "P_{t} vs. P_{l} for #pi^{-} (McTruth); P_{l}; P_{t}", 100,-0.05,0.4,100,0,0.2);
  ntpPiMinus->Project("h_pt_vs_pl_piminus_McTruth", "PiMinus_MC_pt:PiMinus_MC_pz", "MCTruthMatch==1&&Mother==3122");

  TH2F * h_cosTheta_vs_p_piminus_McTruth = new TH2F("h_cosTheta_vs_p_piminus_McTruth", "cos(#Theta) vs. P_{L} for #pi^{-}; P_{L} / GeV/c ; cos(#Theta)", 100,-0.05,0.4,100,-1.1,1.1);
  ntpPiMinus->Project("h_cosTheta_vs_p_piminus_McTruth",  "PiMinus_MC_CosTheta:PiMinus_MC_pz", "MCTruthMatch==1&&Mother==3122");

  TH2F * h_theta_vs_p_piminus_McTruth = new TH2F("h_theta_vs_p_piminus_McTruth", "#Theta vs. P for #pi^{-}; P / GeV/c ; #Theta", 100,0,0.4,100,0,2);
  ntpPiMinus->Project("h_theta_vs_p_piminus_McTruth", "PiMinus_MC_tht:PiMinus_MC_p", "MCTruthMatch==1&&Mother==3122");

	//reco probability for PiMinus
  TH1F * h_reco_p_piminus = new TH1F("h_reco_p_piminus", "reconstructed momentum for #pi^{-}; P/GeV/c; counts", 100, 0,0.4);
	ntpPiMinus->Project("h_reco_p_piminus", "PiMinus_MC_p", "MCTruthMatch==1&&Mother==3122");

	TH1F * h_mc_p_piminus = new TH1F("h_mc_p_piminus", "MC momentum for #pi^{-}; #theta/rad; counts",100,0,0.4);
	ntpMC->Project("h_mc_p_piminus", "p", "moth==1&&pdg==-211");

	TH1F * h_prob_vs_p_piminus = new TH1F("h_prob_vs_p_piminus","reconstruction probability for #pi^{-};P/GeV/c;reco. probability",100,0,0.4);
	h_prob_vs_p_piminus->Divide(h_reco_p_piminus,h_mc_p_piminus);
	h_prob_vs_p_piminus->SetStats(0);
	h_prob_vs_p_piminus->GetYaxis()->SetRangeUser(0,1.1);

  TH1F * h_reco_theta_piminus = new TH1F("h_reco_theta_piminus", "reconstructed angle for #pi^{-}; #theta/rad;counts", 100,0,2);
	ntpPiMinus->Project("h_reco_theta_piminus", "PiMinus_MC_tht", "MCTruthMatch==1&&Mother==3122");

	TH1F * h_mc_theta_piminus = new TH1F("h_mc_theta_piminus", "MC angle for #pi^{-}; #theta/rad; counts",100,0,2);
	ntpMC->Project("h_mc_theta_piminus", "tht", "moth==1&&pdg==-211");

	TH1F * h_prob_vs_theta_piminus = new TH1F("h_prob_vs_theta_piminus","reconstruction probability for #pi^{-};#theta/rad;reco. probability",100,0,2);
	h_prob_vs_theta_piminus->Divide(h_reco_theta_piminus,h_mc_theta_piminus);
	h_prob_vs_theta_piminus->SetStats(0);
	h_prob_vs_theta_piminus->GetYaxis()->SetRangeUser(0,1.1);

	//projections for Proton
  TH2F * h_pt_vs_pl_proton = new TH2F("h_pt_vs_pl_proton", "P_{t} vs. P_{l} for proton; P_{l}; P_{t}", 100,-1,1.8,100,0,0.4);
  ntpProton->Project("h_pt_vs_pl_proton", "Proton_pt:Proton_pz", "MCTruthMatch==1");
 
  TH2F * h_cosTheta_vs_p_proton = new TH2F("h_cosTheta_vs_p_proton", "cos(#Theta) vs. P_{L} for proton; P_{L} / GeV/c ; cos(#Theta)", 100,-1,1.8,100,-1.1,1.1);
  ntpProton->Project("h_cosTheta_vs_p_proton",  "Proton_CosTheta:Proton_pz", "MCTruthMatch==1");

  TH2F * h_theta_vs_p_proton = new TH2F("h_theta_vs_p_proton", "#Theta vs. P for proton; P / GeV/c ; #Theta", 100,-1,1.8,100,0,1);
  ntpProton->Project("h_theta_vs_p_proton", "Proton_tht:Proton_p", "MCTruthMatch==1");

  TH2F * h_pt_vs_pl_proton_McTruth = new TH2F("h_pt_vs_pl_proton_McTruth", "P_{t} vs. P_{l} for proton (McTruth); P_{l}; P_{t}", 100,0,1.5,100,0,0.4);
  ntpProton->Project("h_pt_vs_pl_proton_McTruth", "Proton_MC_pt:Proton_MC_pz", "MCTruthMatch==1&&Mother==3122");
 
  TH2F * h_cosTheta_vs_p_proton_McTruth = new TH2F("h_cosTheta_vs_p_proton_McTruth", "cos(#Theta) vs. P_{L} for proton; P_{L} / GeV/c ; cos(#Theta)", 100,0,1.5,100,-1.1,1.1);
  ntpProton->Project("h_cosTheta_vs_p_proton_McTruth",  "Proton_MC_CosTheta:Proton_MC_pz", "MCTruthMatch==1&&Mother==3122");

  TH2F * h_theta_vs_p_proton_McTruth = new TH2F("h_theta_vs_p_proton_McTruth", "#Theta vs. P for proton; P / GeV/c ; #Theta", 100,0,1.8,100,0,1);
  ntpProton->Project("h_theta_vs_p_proton_McTruth", "Proton_MC_tht:Proton_MC_p", "MCTruthMatch==1&&Mother==3122");

	//reco probability for Proton
  TH1F * h_reco_p_proton = new TH1F("h_reco_p_proton", "reconstructed momentum for proton; P/GeV/c; counts", 100, 0,1.8);
	ntpProton->Project("h_reco_p_proton", "Proton_MC_p", "MCTruthMatch==1&&Mother==3122");

	TH1F * h_mc_p_proton = new TH1F("h_mc_p_proton", "MC momentum for #pi^{-}; P/GeV/c; counts",100,0,1.8);
	ntpMC->Project("h_mc_p_proton", "p", "moth==1&&pdg==2212");

	TH1F * h_prob_vs_p_proton = new TH1F("h_prob_vs_p_proton","reconstruction probability for proton;P/GeV/c;reco. probability",100,0,1.8);
	h_prob_vs_p_proton->Divide(h_reco_p_proton,h_mc_p_proton);
	h_prob_vs_p_proton->SetStats(0);
	h_prob_vs_p_proton->GetYaxis()->SetRangeUser(0,1.1);

  TH1F * h_reco_theta_proton = new TH1F("h_reco_theta_proton", "reconstructed angle for proton; #theta/rad; counts", 100,0,1);
	ntpProton->Project("h_reco_theta_proton", "Proton_MC_tht", "MCTruthMatch==1&&Mother==3122");

	TH1F * h_mc_theta_proton = new TH1F("h_mc_theta_proton", "MC angle for #pi^{-}; #theta/rad; counts",100,0,1);
	ntpMC->Project("h_mc_theta_proton", "tht", "moth==1&&pdg==2212");

	TH1F * h_prob_vs_theta_proton = new TH1F("h_prob_vs_theta_proton","reconstruction probability for proton;#theta/rad;reco. probability",100,0,1);
	h_prob_vs_theta_proton->Divide(h_reco_theta_proton,h_mc_theta_proton);
	h_prob_vs_theta_proton->SetStats(0);
	h_prob_vs_theta_proton->GetYaxis()->SetRangeUser(0,1.1);

	
	//projections for AntiProton
  TH2F * h_pt_vs_pl_antiProton = new TH2F("h_pt_vs_pl_antiProton", "P_{t} vs. P_{l} for antiproton; P_{l}; P_{t}", 100,-1,1.8,100,0,0.4);
  ntpAntiProton->Project("h_pt_vs_pl_antiProton", "AntiProton_pt:AntiProton_pz", "MCTruthMatch==1");
 
  TH2F * h_cosTheta_vs_p_antiProton = new TH2F("h_cosTheta_vs_p_antiProton", "cos(#Theta) vs. P_{L} for antiproton; P_{L} / GeV/c ; cos(#Theta)", 100,-1,1.8,100,-1.1,1.1);
  ntpAntiProton->Project("h_cosTheta_vs_p_antiProton",  "AntiProton_CosTheta:AntiProton_pz", "MCTruthMatch==1");

  TH2F * h_theta_vs_p_antiProton = new TH2F("h_theta_vs_p_antiProton", "#Theta vs. P for antiproton; P / GeV/c ; #Theta", 100,0,1.8,100,0,1);
  ntpAntiProton->Project("h_theta_vs_p_antiProton", "AntiProton_tht:AntiProton_p", "MCTruthMatch==1");

  TH2F * h_pt_vs_pl_antiProton_McTruth = new TH2F("h_pt_vs_pl_antiProton_McTruth", "P_{t} vs. P_{l} for antiproton (McTruth); P_{l}; P_{t}", 100,0,1.5,100,0,0.4);
  ntpAntiProton->Project("h_pt_vs_pl_antiProton_McTruth", "AntiProton_MC_pt:AntiProton_MC_pz", "MCTruthMatch==1&&Mother==-3122");
 
  TH2F * h_cosTheta_vs_p_antiProton_McTruth = new TH2F("h_cosTheta_vs_p_antiProton_McTruth", "cos(#Theta) vs. P_{L} for antiproton; P_{L} / GeV/c ; cos(#Theta)", 100,0,1.5,100,-1.1,1.1);
  ntpAntiProton->Project("h_cosTheta_vs_p_antiProton_McTruth",  "AntiProton_MC_CosTheta:AntiProton_MC_pz", "MCTruthMatch==1&&Mother==-3122");

  TH2F * h_theta_vs_p_antiProton_McTruth = new TH2F("h_theta_vs_p_antiProton_McTruth", "#Theta vs. P for antiproton; P / GeV/c ; #Theta", 100,0,1.8,100,0,1);
  ntpAntiProton->Project("h_theta_vs_p_antiProton_McTruth", "AntiProton_MC_tht:AntiProton_MC_p", "MCTruthMatch==1&&Mother==-3122");

	//reco probability for Anitproton
  TH1F * h_reco_p_antiProton = new TH1F("h_reco_p_antiProton", "reconstructed momentum for #pi^{-}; P/GeV/c; counts", 100, 0,1.8);
	ntpAntiProton->Project("h_reco_p_antiProton", "AntiProton_MC_p", "MCTruthMatch==1&&Mother==-3122");

	TH1F * h_mc_p_antiProton = new TH1F("h_mc_p_antiProton", "MC momentum for antiproton; P/GeV/c; counts",100,0,1.8);
	ntpMC->Project("h_mc_p_antiProton", "p", "moth==2&&pdg==-2212");

	TH1F * h_prob_vs_p_antiProton = new TH1F("h_prob_vs_p_antiProton","reconstruction probability for antiproton;P/GeV/c;reco. probability",100,0,1.8);
	h_prob_vs_p_antiProton->Divide(h_reco_p_antiProton,h_mc_p_antiProton);
	h_prob_vs_p_antiProton->SetStats(0);
	h_prob_vs_p_antiProton->GetYaxis()->SetRangeUser(0,1.1);

  TH1F * h_reco_theta_antiProton = new TH1F("h_reco_theta_antiProton", "reconstructed angle for #pi^{-}; #theta/rad; counts", 100,0,1);
	ntpAntiProton->Project("h_reco_theta_antiProton", "AntiProton_MC_tht", "MCTruthMatch==1&&Mother==-3122");

	TH1F * h_mc_theta_antiProton = new TH1F("h_mc_theta_antiProton", "MC angle for antiproton; #theta/rad; counts",100,0,1);
	ntpMC->Project("h_mc_theta_antiProton", "tht", "moth==2&&pdg==-2212");

	TH1F * h_prob_vs_theta_antiProton = new TH1F("h_prob_vs_theta_antiProton","reconstruction probability for antiproton;#theta/rad;reco. probability",100,0,1);
	h_prob_vs_theta_antiProton->Divide(h_reco_theta_antiProton,h_mc_theta_antiProton);
	h_prob_vs_theta_antiProton->SetStats(0);
	h_prob_vs_theta_antiProton->GetYaxis()->SetRangeUser(0,1.1);
	
  //projections for lambda0
  TH1F * h_fvtx_Chi2 = new TH1F("h_fvtx_Chi2", "#chi^{2} distribution for vertex fit; #Chi^{2}; Counts", 100,0,10);
  ntpLambda0->Project("h_fvtx_Chi2", "fvtx_Chi2");

  TH1F * h_fvtx_Prob = new TH1F("h_fvtx_Prob", "#chi^{2} probability for vertex fit; Prob(#Chi^{2}); Counts", 100,0,1);
  ntpLambda0->Project("h_fvtx_Prob", "fvtx_Prob");

  TH1F * h_lambda0_mass = new TH1F("h_lambda0_mass", "Invariant Mass for  #Lambda^{0}; mass GeV/c^{2}; Counts", 100,0.8,2);
  ntpLambda0->Project("h_lambda0_mass", "Lambda0_m", "McTruthMatch==1");

  TH1F * h_lambda0_mom = new TH1F("h_lambda0_mom", "Momentum for  #Lambda^{0}; P/GeV/c; Counts", 100,0,1.8);
  ntpLambda0->Project("h_lambda0_mom", "Lambda0_p", "McTruthMatch==1");

  TH1F * h_lambda0_energy = new TH1F("h_lambda0_energy", "Energy for  #Lambda^{0}; E/GeV; Counts", 100,0,10);
  ntpLambda0->Project("h_lambda0_energy", "Lambda0_e", "McTruthMatch==1");

  TH1F * h_lambda0_pl = new TH1F("h_lambda0_pl", "P_{L}; P_{L}/GeV/c; Counts", 100,0.3,1.4);
  ntpLambda0->Project("h_lambda0_pl", "Lambda0_pz", "McTruthMatch==1");

  TH2F * h_lambda0_pt_vs_pl = new TH2F("h_lambda0_pt_vs_pl", "P_{t} vs. P_{L}; P_{L}/GeV/c for #Lambda^{0}; P_{t}/GeV/c", 100,0.3,1.4,100,0,0.4);
  ntpLambda0->Project("h_lambda0_pt_vs_pl", "Lambda0_pt:Lambda0_pz", "McTruthMatch==1");

  TH2F * h_lambda0_pt_vs_pl_McTruth = new TH2F("h_lambda0_pt_vs_pl_McTruth", "P_{t} vs. P_{L} for #Lambda^{0} (McTruth); P_{L}/GeV/c; P_{t}/GeV/c", 100,0.3,1.4,100,0,0.4);
  ntpLambda0->Project("h_lambda0_pt_vs_pl_McTruth", "truth_pt:truth_pz", "McTruthMatch==1&&Mother==88888");

  TH1F * h_lambda0_vertex_x = new TH1F("h_lambda0_vertex_x", "x-Position of vertex for  #Lambda^{0}; x/cm; Counts", 100,-5,5);
  ntpLambda0->Project("h_lambda0_vertex_x", "Lambda0_vx", "McTruthMatch==1");

  TH1F * h_lambda0_vertex_y = new TH1F("h_lambda0_vertex_y", "y-Position of vertex for  #Lambda^{0}; y/cm; Counts", 100,-5,5);
  ntpLambda0->Project("h_lambda0_vertex_y", "Lambda0_vy", "McTruthMatch==1");

  TH1F * h_lambda0_vertex_z = new TH1F("h_lambda0_vertex_z", "z-Position of vertex for  #Lambda^{0}; z/cm; Counts", 100,-5,5);
  ntpLambda0->Project("h_lambda0_vertex_z", "Lambda0_vz", "McTruthMatch==1");
  
	//reco probability for Lambda0
  TH1F * h_reco_p_lambda0 = new TH1F("h_reco_p_lambda0", "reconstructed momentum for #pi^{-}; P/GeV/c; counts", 100, 0,1.8);
	ntpLambda0->Project("h_reco_p_lambda0", "truth_p", "McTruthMatch==1&&Mother==88888");

	TH1F * h_mc_p_lambda0 = new TH1F("h_mc_p_lambda0", "MC momentum for #Lambda^{0}; P/GeV/c; counts",100,0,1.8);
	ntpMC->Project("h_mc_p_lambda0", "p", "moth==0&&pdg==3122");

	TH1F * h_prob_vs_p_lambda0 = new TH1F("h_prob_vs_p_lambda0","reconstruction probability for #Lambda^{0};P/GeV/c;reco. probability",100,0,1.8);
	h_prob_vs_p_lambda0->Divide(h_reco_p_lambda0,h_mc_p_lambda0);
	h_prob_vs_p_lambda0->SetStats(0);
	h_prob_vs_p_lambda0->GetYaxis()->SetRangeUser(0,1.1);

  TH1F * h_reco_theta_lambda0 = new TH1F("h_reco_theta_lambda0", "reconstructed angle for #pi^{-}; #theta/rad; counts", 100,0,0.5);
	ntpLambda0->Project("h_reco_theta_lambda0", "truth_tht", "McTruthMatch==1&&Mother==88888");

	TH1F * h_mc_theta_lambda0 = new TH1F("h_mc_theta_lambda0", "MC angle for #Lambda^{0}; #theta/rad; counts",100,0,0.5);
	ntpMC->Project("h_mc_theta_lambda0", "tht", "moth==0&&pdg==3122");

	TH1F * h_prob_vs_theta_lambda0 = new TH1F("h_prob_vs_theta_lambda0","reconstruction probability for #Lambda^{0};#theta/rad;reco. probability",100,0,0.5);
	h_prob_vs_theta_lambda0->Divide(h_reco_theta_lambda0,h_mc_theta_lambda0);
	h_prob_vs_theta_lambda0->SetStats(0);
	h_prob_vs_theta_lambda0->GetYaxis()->SetRangeUser(0,1.1);


  //projections for AntiLambda0
  TH1F * h_fvtx_Chi2_antiLambda0 = new TH1F("h_fvtx_Chi2_antiLambda0", "#chi^{2} distribution for vertex fit; #Chi^{2}; Counts", 100,0,10);
  ntpAntiLambda0->Project("h_fvtx_Chi2_antiLambda0", "fvtx_Chi2");

  TH1F * h_fvtx_Prob_antiLambda0 = new TH1F("h_fvtx_Prob_antiLambda0", "#chi^{2} probability for vertex fit; Prob(#Chi^{2}); Counts", 100,0,1);
  ntpAntiLambda0->Project("h_fvtx_Prob_antiLambda0", "fvtx_Prob");

  TH1F * h_antiLambda0_mass = new TH1F("h_antiLambda0_mass", "Invariant Mass for  #bar{#Lambda}^{0}; mass GeV/c^{2}; Counts", 100,0.8,2);
  ntpAntiLambda0->Project("h_antiLambda0_mass", "AntiLambda0_m", "McTruthMatch==1");

  TH1F * h_antiLambda0_mom = new TH1F("h_antiLambda0_mom", "Momentum for  #bar{#Lambda}^{0}; P/GeV/c; Counts", 100,0,1.8);
  ntpAntiLambda0->Project("h_antiLambda0_mom", "AntiLambda0_p", "McTruthMatch==1");

  TH1F * h_antiLambda0_energy = new TH1F("h_antiLambda0_energy", "Energy for  #bar{#Lambda}^{0}; E/GeV; Counts", 100,0,10);
  ntpAntiLambda0->Project("h_antiLambda0_energy", "AntiLambda0_e", "McTruthMatch==1");

  TH1F * h_antiLambda0_pl = new TH1F("h_antiLambda0_pl", "P_{L}; P_{L}/GeV/c; Counts", 100,0.3,1.4);
  ntpAntiLambda0->Project("h_antiLambda0_pl", "AntiLambda0_pz", "McTruthMatch==1");

  TH2F * h_antiLambda0_pt_vs_pl = new TH2F("h_antiLambda0_pt_vs_pl", "P_{t} vs. P_{L} for #bar{#Lambda}^{0}; P_{L}/GeV/c; P_{t}/GeV/c", 100,0.3,1.4,100,0,0.4);
  ntpAntiLambda0->Project("h_antiLambda0_pt_vs_pl", "AntiLambda0_pt:AntiLambda0_pz", "McTruthMatch==1");

  TH2F * h_antiLambda0_pt_vs_pl_McTruth = new TH2F("h_antiLambda0_pt_vs_pl_McTruth", "P_{t} vs. P_{L} for #bar{#Lambda}^{0} (McTruth); P_{L}/GeV/c; P_{t}/GeV/c", 100,0.3,1.4,100,0,0.4);
  ntpAntiLambda0->Project("h_antiLambda0_pt_vs_pl_McTruth", "truth_pt:truth_pz", "McTruthMatch==1&&Mother==88888");

  TH1F * h_antiLambda0_vertex_x = new TH1F("h_antiLambda0_vertex_x", "x-Position of vertex for  #bar{#Lambda}^{0}; x/cm; Counts", 100,-5,5);
  ntpAntiLambda0->Project("h_antiLambda0_vertex_x", "AntiLambda0_vx", "McTruthMatch==1");

  TH1F * h_antiLambda0_vertex_y = new TH1F("h_antiLambda0_vertex_y", "y-Position of vertex for  #bar{#Lambda}^{0}; y/cm; Counts", 100,-5,5);
  ntpAntiLambda0->Project("h_antiLambda0_vertex_y", "AntiLambda0_vy", "McTruthMatch==1");

  TH1F * h_antiLambda0_vertex_z = new TH1F("h_antiLambda0_vertex_z", "z-Position of vertex for  #bar{#Lambda}^{0}; z/cm; Counts", 100,-5,5);
  ntpAntiLambda0->Project("h_antiLambda0_vertex_z", "AntiLambda0_vz", "McTruthMatch==1");


	//reco probability for AntiLambda0
  TH1F * h_reco_p_antiLambda0 = new TH1F("h_reco_p_antiLambda0", "reconstructed momentum for #pi^{-}; P/GeV/c; counts", 100, 0,1.8);
	ntpAntiLambda0->Project("h_reco_p_antiLambda0", "truth_p", "McTruthMatch==1&&Mother==88888");

	TH1F * h_mc_p_antiLambda0 = new TH1F("h_mc_p_antiLambda0", "MC momentum for #bar{#Lambda}^{0}; P/GeV/c; counts",100,0,1.8);
	ntpMC->Project("h_mc_p_antiLambda0", "p", "moth==0&&pdg==-3122");

	TH1F * h_prob_vs_p_antiLambda0 = new TH1F("h_prob_vs_p_antiLambda0","reconstruction probability for #bar{#Lambda}^{0};P/GeV/c;reco. probability",100,0,1.8);
	h_prob_vs_p_antiLambda0->Divide(h_reco_p_antiLambda0,h_mc_p_antiLambda0);
	h_prob_vs_p_antiLambda0->SetStats(0);
	h_prob_vs_p_antiLambda0->GetYaxis()->SetRangeUser(0,1.1);

  TH1F * h_reco_theta_antiLambda0 = new TH1F("h_reco_theta_antiLambda0", "reconstructed angle for #pi^{-};#theta/rad ; counts", 100,0,0.5);
	ntpAntiLambda0->Project("h_reco_theta_antiLambda0", "truth_tht", "McTruthMatch==1&&Mother==88888");

	TH1F * h_mc_theta_antiLambda0 = new TH1F("h_mc_theta_antiLambda0", "MC angle for #bar{#Lambda}^{0}; #theta/rad ; counts",100,0,0.5);
	ntpMC->Project("h_mc_theta_antiLambda0", "tht", "moth==0&&pdg==-3122");

	TH1F * h_prob_vs_theta_antiLambda0 = new TH1F("h_prob_vs_theta_antiLambda0","reconstruction probability for #bar{#Lambda}^{0};#theta/rad ;reco. probability",100,0,0.5);
	h_prob_vs_theta_antiLambda0->Divide(h_reco_theta_antiLambda0,h_mc_theta_antiLambda0);
	h_prob_vs_theta_antiLambda0->SetStats(0);
	h_prob_vs_theta_antiLambda0->GetYaxis()->SetRangeUser(0,1.1);
 

	//***Create Histograms

	//MC
  TCanvas * c_theta_vs_p_mc = new TCanvas("c_theta_vs_p_mc", "theta_vs_p_mc", 0,0,800,500);
  h_theta_vs_p_mc->Draw("colz");

  TCanvas * c_costheta_vs_p_mc = new TCanvas("c_costheta_vs_p_mc", "costheta_vs_p_mc", 0,0,800,500);
  h_cosTheta_vs_p_mc->Draw("colz");

  TCanvas * c_pt_vs_pl_mc = new TCanvas("c_pt_vs_pl_mc", "P_{t} vs P_{l} mc", 0,0,800,500);
  h_pt_vs_pl_mc->Draw("colz");



  TCanvas * c_theta_vs_p_mc_piplus = new TCanvas("c_theta_vs_p_mc_piplus", "theta_vs_p_mc_piplus", 0,0,800,500);
  h_theta_vs_p_mc_piplus->Draw("colz");

  TCanvas * c_costheta_vs_p_mc_piplus = new TCanvas("c_costheta_vs_p_mc_piplus", "costheta_vs_p_mc_piplus", 0,0,800,500);
  h_cosTheta_vs_p_mc_piplus->Draw("colz");

  TCanvas * c_pt_vs_pl_mc_piplus = new TCanvas("c_pt_vs_pl_mc_piplus", "P_{t} vs P_{l} mc", 0,0,800,500);
  h_pt_vs_pl_mc_piplus->Draw("colz");



	TCanvas * c_theta_vs_p_mc_piminus = new TCanvas("c_theta_vs_p_mc_piminus", "theta_vs_p_mc_piminus", 0,0,800,500);
  h_theta_vs_p_mc_piminus->Draw("colz");

  TCanvas * c_costheta_vs_p_mc_piminus = new TCanvas("c_costheta_vs_p_mc_piminus", "costheta_vs_p_mc_piminus", 0,0,800,500);
  h_cosTheta_vs_p_mc_piminus->Draw("colz");

  TCanvas * c_pt_vs_pl_mc_piminus = new TCanvas("c_pt_vs_pl_mc_piminus", "P_{t} vs P_{l} mc", 0,0,800,500);
  h_pt_vs_pl_mc_piminus->Draw("colz");

 


	TCanvas * c_theta_vs_p_mc_proton = new TCanvas("c_theta_vs_p_mc_proton", "theta_vs_p_mc_proton", 0,0,800,500);
  h_theta_vs_p_mc_proton->Draw("colz");

  TCanvas * c_costheta_vs_p_mc_proton = new TCanvas("c_costheta_vs_p_mc_proton", "costheta_vs_p_mc_proton", 0,0,800,500);
  h_cosTheta_vs_p_mc_proton->Draw("colz");

  TCanvas * c_pt_vs_pl_mc_proton = new TCanvas("c_pt_vs_pl_mc_proton", "P_{t} vs P_{l} mc", 0,0,800,500);
  h_pt_vs_pl_mc_proton->Draw("colz");

     

	TCanvas * c_theta_vs_p_mc_antiProton = new TCanvas("c_theta_vs_p_mc_antiProton", "theta_vs_p_mc_antiProton", 0,0,800,500);
  h_theta_vs_p_mc_antiProton->Draw("colz");

  TCanvas * c_costheta_vs_p_mc_antiProton = new TCanvas("c_costheta_vs_p_mc_antiProton", "costheta_vs_p_mc_antiProton", 0,0,800,500);
  h_cosTheta_vs_p_mc_antiProton->Draw("colz");

  TCanvas * c_pt_vs_pl_mc_antiProton = new TCanvas("c_pt_vs_pl_mc_antiProton", "P_{t} vs P_{l} mc", 0,0,800,500);
  h_pt_vs_pl_mc_antiProton->Draw("colz");

	

	TCanvas * c_theta_vs_p_mc_lambda0 = new TCanvas("c_theta_vs_p_mc_lambda0", "theta_vs_p_mc_lambda0", 0,0,800,500);
  h_theta_vs_p_mc_lambda0->Draw("colz");

  TCanvas * c_costheta_vs_p_mc_lambda0 = new TCanvas("c_costheta_vs_p_mc_lambda0", "costheta_vs_p_mc_lambda0", 0,0,800,500);
  h_cosTheta_vs_p_mc_lambda0->Draw("colz");

  TCanvas * c_pt_vs_pl_mc_lambda0 = new TCanvas("c_pt_vs_pl_mc_lambda0", "P_{t} vs P_{l} mc", 0,0,800,500);
  h_pt_vs_pl_mc_lambda0->Draw("colz");
 
	

	TCanvas * c_theta_vs_p_mc_antiLambda0 = new TCanvas("c_theta_vs_p_mc_antiLambda0", "theta_vs_p_mc_antiLambda0", 0,0,800,500);
  h_theta_vs_p_mc_antiLambda0->Draw("colz");

  TCanvas * c_costheta_vs_p_mc_antiLambda0 = new TCanvas("c_costheta_vs_p_mc_antiLambda0", "costheta_vs_p_mc_antiLambda0", 0,0,800,500);
  h_cosTheta_vs_p_mc_antiLambda0->Draw("colz");

  TCanvas * c_pt_vs_pl_mc_antiLambda0 = new TCanvas("c_pt_vs_pl_mc_antiLambda0", "P_{t} vs P_{l} mc", 0,0,800,500);
  h_pt_vs_pl_mc_antiLambda0->Draw("colz");

  //PiPlus 
  TCanvas * c_pt_vs_pl_piplus = new TCanvas("c_pt_vs_pl_piplus", "P_{t} vs P_{l} PiPlus", 0,0,800,500);
  h_pt_vs_pl_piplus->Draw("colz");

  TCanvas * c_costheta_vs_p_piplus = new TCanvas("c_costheta_vs_p_piplus", "costheta_vs_p_mc", 0,0,800,500);
  h_cosTheta_vs_p_piplus->Draw("colz");

  TCanvas * c_theta_vs_p_piplus = new TCanvas("c_theta_vs_p_piplus", "theta_vs_p_mc", 0,0,800,500);
  h_theta_vs_p_piplus->Draw("colz");

  TCanvas * c_pt_vs_pl_piplus_McTruth = new TCanvas("c_pt_vs_pl_piplus_McTruth", "P_{t} vs P_{l} PiPlus", 0,0,800,500);
  h_pt_vs_pl_piplus_McTruth->Draw("colz");

  TCanvas * c_costheta_vs_p_piplus_McTruth = new TCanvas("c_costheta_vs_p_piplus_McTruth", "costheta_vs_p_mc", 0,0,800,500);
  h_cosTheta_vs_p_piplus_McTruth->Draw("colz");

  TCanvas * c_theta_vs_p_piplus_McTruth = new TCanvas("c_theta_vs_p_piplus_McTruth", "theta_vs_p_mc", 0,0,800,500);
  h_theta_vs_p_piplus_McTruth->Draw("colz");

	TCanvas * c_prob_vs_p_piplus = new TCanvas("c_prob_vs_p_piplus", "reco. probability for #pi^{+}", 0,0,800,500);
	h_prob_vs_p_piplus->Draw();

	TCanvas * c_prob_vs_theta_piplus = new TCanvas("c_prob_vs_theta_piplus", "reco. probability for #pi^{+}", 0,0,800,500);
	h_prob_vs_theta_piplus->Draw();

  //PiMinus 
  TCanvas * c_pt_vs_pl_piminus = new TCanvas("c_pt_vs_pl_piminus", "P_{t} vs P_{l} PiMinus", 0,0,800,500);
  h_pt_vs_pl_piminus->Draw("colz");

  TCanvas * c_costheta_vs_p_piminus = new TCanvas("c_costheta_vs_p_piminus", "costheta_vs_p_mc", 0,0,800,500);
  h_cosTheta_vs_p_piminus->Draw("colz");

  TCanvas * c_theta_vs_p_piminus = new TCanvas("c_theta_vs_p_piminus", "theta_vs_p_mc", 0,0,800,500);
  h_theta_vs_p_piminus->Draw("colz");

  TCanvas * c_pt_vs_pl_piminus_McTruth = new TCanvas("c_pt_vs_pl_piminus_McTruth", "P_{t} vs P_{l} PiMinus", 0,0,800,500);
  h_pt_vs_pl_piminus_McTruth->Draw("colz");

  TCanvas * c_costheta_vs_p_piminus_McTruth = new TCanvas("c_costheta_vs_p_piminus_McTruth", "costheta_vs_p_mc", 0,0,800,500);
  h_cosTheta_vs_p_piminus_McTruth->Draw("colz");

  TCanvas * c_theta_vs_p_piminus_McTruth = new TCanvas("c_theta_vs_p_piminus_McTruth", "theta_vs_p_mc", 0,0,800,500);
  h_theta_vs_p_piminus_McTruth->Draw("colz");

	TCanvas * c_prob_vs_p_piminus = new TCanvas("c_prob_vs_p_piminus", "reco. probability for #pi^{-}", 0,0,800,500);
	h_prob_vs_p_piminus->Draw();

	TCanvas * c_prob_vs_theta_piminus = new TCanvas("c_prob_vs_theta_piminus", "reco. probability for #pi^{-}", 0,0,800,500);
	h_prob_vs_theta_piminus->Draw();

	//Proton	
  TCanvas * c_pt_vs_pl_proton = new TCanvas("c_pt_vs_pl_proton", "P_{t} vs P_{l} PiMinus", 0,0,800,500);
  h_pt_vs_pl_proton->Draw("colz");

  TCanvas * c_costheta_vs_p_proton = new TCanvas("c_costheta_vs_p_proton", "costheta_vs_p_mc", 0,0,800,500);
  h_cosTheta_vs_p_proton->Draw("colz");

  TCanvas * c_theta_vs_p_proton = new TCanvas("c_theta_vs_p_proton", "theta_vs_p_mc", 0,0,800,500);
  h_theta_vs_p_proton->Draw("colz");

  TCanvas * c_pt_vs_pl_proton_McTruth = new TCanvas("c_pt_vs_pl_proton_McTruth", "P_{t} vs P_{l} PiMinus", 0,0,800,500);
  h_pt_vs_pl_proton_McTruth->Draw("colz");

  TCanvas * c_costheta_vs_p_proton_McTruth = new TCanvas("c_costheta_vs_p_proton_McTruth", "costheta_vs_p_mc", 0,0,800,500);
  h_cosTheta_vs_p_proton_McTruth->Draw("colz");

  TCanvas * c_theta_vs_p_proton_McTruth = new TCanvas("c_theta_vs_p_proton_McTruth", "theta_vs_p_mc", 0,0,800,500);
  h_theta_vs_p_proton_McTruth->Draw("colz");

	TCanvas * c_prob_vs_p_proton = new TCanvas("c_prob_vs_p_proton", "reco. probability for proton", 0,0,800,500);
	h_prob_vs_p_proton->Draw();

	TCanvas * c_prob_vs_theta_proton = new TCanvas("c_prob_vs_tht_proton", "reco. probability for proton", 0,0,800,500);
	h_prob_vs_theta_proton->Draw();

	//AntiProton
  TCanvas * c_pt_vs_pl_antiProton = new TCanvas("c_pt_vs_pl_antiProton", "P_{t} vs P_{l} PiMinus", 0,0,800,500);
  h_pt_vs_pl_antiProton->Draw("colz");

  TCanvas * c_costheta_vs_p_antiProton = new TCanvas("c_costheta_vs_p_antiProton", "costheta_vs_p_mc", 0,0,800,500);
  h_cosTheta_vs_p_antiProton->Draw("colz");

  TCanvas * c_theta_vs_p_antiProton = new TCanvas("c_theta_vs_p_antiProton", "theta_vs_p_mc", 0,0,800,500);
  h_theta_vs_p_antiProton->Draw("colz");

  TCanvas * c_pt_vs_pl_antiProton_McTruth = new TCanvas("c_pt_vs_pl_antiProton_McTruth", "P_{t} vs P_{l} PiMinus", 0,0,800,500);
  h_pt_vs_pl_antiProton_McTruth->Draw("colz");

  TCanvas * c_costheta_vs_p_antiProton_McTruth = new TCanvas("c_costheta_vs_p_antiProton_McTruth", "costheta_vs_p_mc", 0,0,800,500);
  h_cosTheta_vs_p_antiProton_McTruth->Draw("colz");

  TCanvas * c_theta_vs_p_antiProton_McTruth = new TCanvas("c_theta_vs_p_antiProton_McTruth", "theta_vs_p_mc", 0,0,800,500);
  h_theta_vs_p_antiProton_McTruth->Draw("colz");
 
	TCanvas * c_prob_vs_p_antiProton = new TCanvas("c_prob_vs_p_antiProton", "reco. probability for antiproton", 0,0,800,500);
	h_prob_vs_p_antiProton->Draw();

	TCanvas * c_prob_vs_theta_antiProton = new TCanvas("c_prob_vs_tht_antiProton", "reco. probability for antiproton", 0,0,800,500);
	h_prob_vs_theta_antiProton->Draw();
 
	//Lambda0
  TCanvas * c_fvtx_Chi2 = new TCanvas("c_fvxt_Chi2", "fvxt_Chi2", 0,0,800,500);
  h_fvtx_Chi2->Draw();
  
  TCanvas * c_fvtx_Prob = new TCanvas("c_fvxt_Prob", "fvxt_Prob", 0,0,800,500);
  h_fvtx_Prob->Draw();
  
  TCanvas * c_lambda0_mass = new TCanvas("c_lambda0_mass", "lambda0_mass", 0,0,800,500);
  h_lambda0_mass->Draw();

  TCanvas * c_lambda0_mom = new TCanvas("c_lambda0_mom", "lambda0_mom", 0,0,800,500);
  h_lambda0_mom->Draw();

  TCanvas * c_lambda0_energy = new TCanvas("c_lambda0_energy", "lambda0_energy", 0,0,800,500);
  h_lambda0_energy->Draw();

  TCanvas * c_lambda0_pl = new TCanvas("c_lambda0_pl", "lambda0_pl", 0,0,800,500);
  h_lambda0_pl->Draw();

  TCanvas * c_lambda0_pt_vs_pl = new TCanvas("c_lambda0_pt_vs_pl", "lambda0_pt_vs_pl", 0,0,800,500);
  h_lambda0_pt_vs_pl->Draw("colz");

  TCanvas * c_lambda0_pt_vs_pl_McTruth = new TCanvas("c_lambda0_pt_vs_pl_McTruth", "lambda0_pt_vs_pl_McTruth", 0,0,800,500);
  h_lambda0_pt_vs_pl_McTruth->Draw("colz");

  TCanvas * c_lambda0_vertex_x = new TCanvas("c_lambda0_vertex_x", "lambda0_vertex_x", 0,0,800,500);
  h_lambda0_vertex_x->Draw();
 
  TCanvas * c_lambda0_vertex_y = new TCanvas("c_lambda0_vertex_y", "lambda0_vertex_y", 0,0,800,500);
  h_lambda0_vertex_y->Draw();

  TCanvas * c_lambda0_vertex_z = new TCanvas("c_lambda0_vertex_z", "lambda0_vertex_z", 0,0,800,500);
  h_lambda0_vertex_z->Draw();

	TCanvas * c_prob_vs_p_lambda0 = new TCanvas("c_prob_vs_p_lambda0", "reco. probability for #Lambda^{0}", 0,0,800,500);
	h_prob_vs_p_lambda0->Draw();

	TCanvas * c_prob_vs_theta_lambda0 = new TCanvas("c_prob_vs_theta_lambda0", "reco. probability for #Lambda^{0}", 0,0,800,500);
	h_prob_vs_theta_lambda0->Draw();

	//AntiLambda0
  TCanvas * c_fvtx_Chi2_antiLambda0 = new TCanvas("c_fvxt_Chi2_antiLambda0", "fvxt_Chi2", 0,0,800,500);
  h_fvtx_Chi2_antiLambda0->Draw();
  
  TCanvas * c_fvtx_Prob_antiLambda0 = new TCanvas("c_fvxt_Prob_antiLambda0", "fvxt_Prob", 0,0,800,500);
  h_fvtx_Prob_antiLambda0->Draw();
  
  TCanvas * c_antiLambda0_mass = new TCanvas("c_antiLambda0_mass", "antiLambda0_mass", 0,0,800,500);
  h_antiLambda0_mass->Draw();

  TCanvas * c_antiLambda0_mom = new TCanvas("c_antiLambda0_mom", "antiLambda0_mom", 0,0,800,500);
  h_antiLambda0_mom->Draw();

  TCanvas * c_antiLambda0_energy = new TCanvas("c_antiLambda0_energy", "antiLambda0_energy", 0,0,800,500);
  h_antiLambda0_energy->Draw();

  TCanvas * c_antiLambda0_pl = new TCanvas("c_antiLambda0_pl", "antiLambda0_pl", 0,0,800,500);
  h_antiLambda0_pl->Draw();

  TCanvas * c_antiLambda0_pt_vs_pl = new TCanvas("c_antiLambda0_pt_vs_pl", "antiLambda0_pt_vs_pl", 0,0,800,500);
  h_antiLambda0_pt_vs_pl->Draw("colz");

  TCanvas * c_antiLambda0_pt_vs_pl_McTruth = new TCanvas("c_antiLambda0_pt_vs_pl_McTruth", "antiLambda0_pt_vs_pl_McTruth", 0,0,800,500);
  h_antiLambda0_pt_vs_pl_McTruth->Draw("colz");

  TCanvas * c_antiLambda0_vertex_x = new TCanvas("c_antiLambda0_vertex_x", "antiLambda0_vertex_x", 0,0,800,500);
  h_antiLambda0_vertex_x->Draw();
 
  TCanvas * c_antiLambda0_vertex_y = new TCanvas("c_antiLambda0_vertex_y", "antiLambda0_vertex_y", 0,0,800,500);
  h_antiLambda0_vertex_y->Draw();

  TCanvas * c_antiLambda0_vertex_z = new TCanvas("c_antiLambda0_vertex_z", "antiLambda0_vertex_z", 0,0,800,500);
  h_antiLambda0_vertex_z->Draw();

	TCanvas * c_prob_vs_p_antiLambda0 = new TCanvas("c_prob_vs_p_antiLambda0", "reco. probability for #bar{#Lambda}^{0}", 0,0,800,500);
	h_prob_vs_p_antiLambda0->Draw();

	TCanvas * c_prob_vs_theta_antiLambda0 = new TCanvas("c_prob_vs_theta_antiLambda0", "reco. probability for #bar{#Lambda}^{0}", 0,0,800,500);
	h_prob_vs_theta_antiLambda0->Draw();

	
  //***Save output if saveOutput is true	
  if (saveOutput==true) {
		//MC	
    c_theta_vs_p_mc->Print(inputDir + "plots/theta_vs_p_MC.root");
    c_costheta_vs_p_mc->Print(inputDir + "plots/costheta_vs_p_MC.root");
    c_pt_vs_pl_mc->Print(inputDir + "plots/pt_vs_pl_MC.root");

    c_theta_vs_p_mc_piplus->Print(inputDir + "plots/theta_vs_p_MC_piplus.root");
    c_costheta_vs_p_mc_piplus->Print(inputDir + "plots/costheta_vs_p_MC_piplus.root");
    c_pt_vs_pl_mc_piplus->Print(inputDir + "plots/pt_vs_pl_MC_piplus.root");
 

		c_theta_vs_p_mc_piminus->Print(inputDir + "plots/theta_vs_p_MC_piminus.root");
    c_costheta_vs_p_mc_piminus->Print(inputDir + "plots/costheta_vs_p_MC_piminus.root");
    c_pt_vs_pl_mc_piminus->Print(inputDir + "plots/pt_vs_pl_MC_piminus.root");

    c_theta_vs_p_mc_proton->Print(inputDir + "plots/theta_vs_p_MC_proton.root");
    c_costheta_vs_p_mc_proton->Print(inputDir + "plots/costheta_vs_p_mc_proton.root");
    c_pt_vs_pl_mc_proton->Print(inputDir + "plots/pt_vs_pl_mc_proton.root");
 

		c_theta_vs_p_mc_antiProton->Print(inputDir + "plots/theta_vs_p_MC_antiProton.root");
    c_costheta_vs_p_mc_antiProton->Print(inputDir + "plots/costheta_vs_p_MC_antiProton.root");
    c_pt_vs_pl_mc_antiProton->Print(inputDir + "plots/pt_vs_pl_MC_antiProton.root");
		
		c_theta_vs_p_mc_antiLambda0->Print(inputDir + "plots/theta_vs_p_MC_antiLambda0.root");
    c_costheta_vs_p_mc_antiLambda0->Print(inputDir + "plots/costheta_vs_p_MC_antiLambda0.root");
    c_pt_vs_pl_mc_antiLambda0->Print(inputDir + "plots/pt_vs_pl_MC_antiLambda0.root");

		c_theta_vs_p_mc_lambda0->Print(inputDir + "plots/theta_vs_p_MC_lambda0.root");
    c_costheta_vs_p_mc_lambda0->Print(inputDir + "plots/costheta_vs_p_MC_lambda0.root");
    c_pt_vs_pl_mc_lambda0->Print(inputDir + "plots/pt_vs_pl_MC_lambda0.root");

		//PiPlus
    c_theta_vs_p_piplus->Print(inputDir + "plots/theta_vs_p_piplus.root");
    c_costheta_vs_p_piplus->Print(inputDir + "plots/costheta_vs_p_piplus.root");
    c_pt_vs_pl_piplus->Print(inputDir + "plots/pt_vs_pl_piplus.root");
		
    c_theta_vs_p_piplus_McTruth->Print(inputDir + "plots/theta_vs_p_piplus_McTruth.root");
    c_costheta_vs_p_piplus_McTruth->Print(inputDir + "plots/costheta_vs_p_piplus_McTruth.root");
    c_pt_vs_pl_piplus_McTruth->Print(inputDir + "plots/pt_vs_pl_piplus_McTruth.root");

		c_prob_vs_p_piplus->Print(inputDir + "plots/prob_vs_p_piplus.root");
		c_prob_vs_theta_piplus->Print(inputDir + "plots/prob_vs_theta_piplus.root");

		//PiMinus
    c_theta_vs_p_piminus->Print(inputDir + "plots/theta_vs_p_piminus.root");
    c_costheta_vs_p_piminus->Print(inputDir + "plots/costheta_vs_p_piminus.root");
    c_pt_vs_pl_piminus->Print(inputDir + "plots/pt_vs_pl_piminus.root");

    c_theta_vs_p_piminus_McTruth->Print(inputDir + "plots/theta_vs_p_piminus_McTruth.root");
    c_costheta_vs_p_piminus_McTruth->Print(inputDir + "plots/costheta_vs_p_piminus_McTruth.root");
    c_pt_vs_pl_piminus_McTruth->Print(inputDir + "plots/pt_vs_pl_piminus_McTruth.root");

		c_prob_vs_p_piminus->Print(inputDir + "plots/prob_vs_p_piminus.root");
		c_prob_vs_theta_piminus->Print(inputDir + "plots/prob_vs_theta_piminus.root");


		//Proton
    c_theta_vs_p_proton->Print(inputDir + "plots/theta_vs_p_proton.root");
    c_costheta_vs_p_proton->Print(inputDir + "plots/costheta_vs_p_proton.root");
    c_pt_vs_pl_proton->Print(inputDir + "plots/pt_vs_pl_proton.root");

    c_theta_vs_p_proton_McTruth->Print(inputDir + "plots/theta_vs_p_proton_McTruth.root");
    c_costheta_vs_p_proton_McTruth->Print(inputDir + "plots/costheta_vs_p_proton_McTruth.root");
    c_pt_vs_pl_proton_McTruth->Print(inputDir + "plots/pt_vs_pl_proton_McTruth.root");

		c_prob_vs_p_proton->Print(inputDir + "plots/prob_vs_p_proton.root");
		c_prob_vs_tht_proton->Print(inputDir + "plots/prob_vs_tht_proton.root");

		//AntiProton
    c_theta_vs_p_antiProton->Print(inputDir + "plots/theta_vs_p_antiProton.root");
    c_costheta_vs_p_antiProton->Print(inputDir + "plots/costheta_vs_p_antiProton.root");
    c_pt_vs_pl_antiProton->Print(inputDir + "plots/pt_vs_pl_antiProton.root");
	
    c_theta_vs_p_antiProton_McTruth->Print(inputDir + "plots/theta_vs_p_antiProton_McTruth.root");
    c_costheta_vs_p_antiProton_McTruth->Print(inputDir + "plots/costheta_vs_p_antiProton_McTruth.root");
    c_pt_vs_pl_antiProton_McTruth->Print(inputDir + "plots/pt_vs_pl_antiProton_McTruth.root");

		c_prob_vs_p_antiProton->Print(inputDir + "plots/prob_vs_p_antiProton.root");
		c_prob_vs_tht_antiProton->Print(inputDir + "plots/prob_vs_tht_antiProton.root");

		//Lambda0
		c_fvtx_Chi2->Print(inputDir + "plots/fvtx_Chi2_Lambda0.root");	
		c_fvtx_Prob->Print(inputDir + "plots/fvtx_Prob_Lambda0.root");	
		c_lambda0_mass->Print(inputDir + "plots/lambda0_mass.root");	
		c_lambda0_mom->Print(inputDir + "plots/lambda0_mom.root");	
		c_lambda0_energy->Print(inputDir + "plots/lambda0_energy.root");	
		c_lambda0_pl->Print(inputDir + "plots/lambda0_pl.root");	
		c_lambda0_pt_vs_pl->Print(inputDir + "plots/lambda0_pt_vs_pl.root");
		c_lambda0_pt_vs_pl_McTruth->Print(inputDir + "plots/lambda0_pt_vs_pl_McTruth.root");
		c_lambda0_vertex_x->Print(inputDir + "plots/lambda0_vertex_x.root");	
		c_lambda0_vertex_y->Print(inputDir + "plots/lambda0_vertex_y.root");
		c_lambda0_vertex_z->Print(inputDir + "plots/lambda0_vertex_z.root");	
	
		c_prob_vs_p_lambda0->Print(inputDir + "plots/prob_vs_p_lambda0.root");
		c_prob_vs_theta_lambda0->Print(inputDir + "plots/prob_vs_theta_lambda0.root");

		//AntiLambda0	
		c_fvtx_Chi2_antiLambda0->Print(inputDir + "plots/fvtx_Chi2_antiLambda0.root");	
		c_fvtx_Prob_antiLambda0->Print(inputDir + "plots/fvtx_Prob_antiLambda0.root");	
		c_antiLambda0_mass->Print(inputDir + "plots/antiLambda0_mass.root");	
		c_antiLambda0_mom->Print(inputDir + "plots/antiLambda0_mom.root");	
		c_antiLambda0_energy->Print(inputDir + "plots/antiLambda0_energy.root");	
		c_antiLambda0_pl->Print(inputDir + "plots/antiLambda0_pl.root");	
		c_antiLambda0_pt_vs_pl->Print(inputDir + "plots/antiLambda0_pt_vs_pl.root");
		c_antiLambda0_pt_vs_pl_McTruth->Print(inputDir + "plots/antiLambda0_pt_vs_pl_McTruth.root");
		c_antiLambda0_vertex_x->Print(inputDir + "plots/antiLambda0_vertex_x.root");	
		c_antiLambda0_vertex_y->Print(inputDir + "plots/antiLambda0_vertex_y.root");
		c_antiLambda0_vertex_z->Print(inputDir + "plots/antiLambda0_vertex_z.root");	
		
		c_prob_vs_p_antiLambda0->Print(inputDir + "plots/prob_vs_p_antiLambda0.root");
		c_prob_vs_theta_antiLambda0->Print(inputDir + "plots/prob_vs_theta_antiLambda0.root");
	
  }


  //***close histograms
 
	//MC 
  c_theta_vs_p_mc->Close();
  c_costheta_vs_p_mc->Close();
  c_pt_vs_pl_mc->Close();

  c_theta_vs_p_mc_piplus->Close();
  c_costheta_vs_p_mc_piplus->Close();
  c_pt_vs_pl_mc_piplus->Close();

  c_theta_vs_p_mc_piminus->Close();
  c_costheta_vs_p_mc_piminus->Close();
  c_pt_vs_pl_mc_piminus->Close();

  c_theta_vs_p_mc_proton->Close();
  c_costheta_vs_p_mc_proton->Close();
  c_pt_vs_pl_mc_proton->Close();
  
	c_theta_vs_p_mc_antiProton->Close();
  c_costheta_vs_p_mc_antiProton->Close();
  c_pt_vs_pl_mc_antiProton->Close();

	c_theta_vs_p_mc_lambda0->Close();
  c_costheta_vs_p_mc_lambda0->Close();
  c_pt_vs_pl_mc_lambda0->Close();

	c_theta_vs_p_mc_antiLambda0->Close();
  c_costheta_vs_p_mc_antiLambda0->Close();
  c_pt_vs_pl_mc_antiLambda0->Close();
		
	//PiPlus   
  c_pt_vs_pl_piplus->Close();
  c_costheta_vs_p_piplus->Close();
  c_theta_vs_p_piplus->Close();
  
  c_pt_vs_pl_piplus_McTruth->Close();
  c_costheta_vs_p_piplus_McTruth->Close();
  c_theta_vs_p_piplus_McTruth->Close();
	
	c_prob_vs_p_piplus ->Close();
	c_prob_vs_theta_piplus ->Close();
	
	//PiMinus
  c_pt_vs_pl_piminus->Close();
  c_costheta_vs_p_piminus->Close();
  c_theta_vs_p_piminus->Close();
  
  c_pt_vs_pl_piminus_McTruth->Close();
  c_costheta_vs_p_piminus_McTruth->Close();
  c_theta_vs_p_piminus_McTruth->Close();
	
	c_prob_vs_p_piminus ->Close();
	c_prob_vs_theta_piminus ->Close();
		
	//Proton
  c_pt_vs_pl_proton->Close();
  c_costheta_vs_p_proton->Close();
  c_theta_vs_p_proton->Close();
  
  c_pt_vs_pl_proton_McTruth->Close();
  c_costheta_vs_p_proton_McTruth->Close();
  c_theta_vs_p_proton_McTruth->Close();
	
	c_prob_vs_p_proton ->Close();
	c_prob_vs_tht_proton ->Close();
	
	//AntiProton
  c_pt_vs_pl_antiProton->Close();
  c_costheta_vs_p_antiProton->Close();
  c_theta_vs_p_antiProton->Close();
	
  c_pt_vs_pl_antiProton_McTruth->Close();
  c_costheta_vs_p_antiProton_McTruth->Close();
  c_theta_vs_p_antiProton_McTruth->Close();
	
	c_prob_vs_p_antiProton ->Close();
	c_prob_vs_tht_antiProton ->Close();
	
 	//Lambda0 
  c_fvtx_Chi2->Close();
  c_fvtx_Prob->Close();
  c_lambda0_pl->Close();
  c_lambda0_pt_vs_pl->Close();
  c_lambda0_pt_vs_pl_McTruth->Close();
  c_lambda0_vertex_x->Close();
  c_lambda0_vertex_y->Close();
  c_lambda0_vertex_z->Close();
  c_lambda0_mass->Close();
  c_lambda0_mom->Close();
  c_lambda0_energy->Close();
	
	c_prob_vs_p_lambda0 ->Close();
	c_prob_vs_theta_lambda0 ->Close();
	
	//AntiLambda0	
  c_fvtx_Chi2_antiLambda0->Close();
  c_fvtx_Prob_antiLambda0->Close();
  c_antiLambda0_pl->Close();
  c_antiLambda0_pt_vs_pl->Close();
  c_antiLambda0_pt_vs_pl_McTruth->Close();
  c_antiLambda0_vertex_x->Close();
  c_antiLambda0_vertex_y->Close();
  c_antiLambda0_vertex_z->Close();
  c_antiLambda0_mass->Close();
  c_antiLambda0_mom->Close();
  c_antiLambda0_energy->Close();
	
	c_prob_vs_p_antiLambda0 ->Close();
 	c_prob_vs_theta_antiLambda0 ->Close();
	
}

int main(){
  evaluation();
  return 0;
}

