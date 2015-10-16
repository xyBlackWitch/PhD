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


void eval_pbarp_2pi(TString inputDir="/home/ikp1/puetz/panda/fairsoft_mar15/pandaroot/mysimulations/test/2pi_10000_events/", bool saveOutput = true){

  TString filename_ana = inputDir + "output_ana.root";
  
  //***open input
  TFile * data = new TFile(filename_ana, "READ");
  TTree * ntpBeam = (TTree*)data->Get("ntpBeam");
  TTree * ntpPiPlus= (TTree*)data->Get("ntpPiPlus");
  TTree * ntpPiMinus = (TTree*)data->Get("ntpPiMinus");
  TTree * ntpMc = (TTree*)data->Get("ntpMC");

  //projection for MC

  TH2F * h_theta_vs_p_mc = new TH2F("h_theta_vs_p_mc", "#Theta vs. P for MC; P / GeV/c ; #Theta", 100,0,10,100,0,3);
  ntpMC->Project("h_theta_vs_p_mc", "tht:p", "moth==0" );

  TH2F * h_costheta_vs_p_mc = new TH2F("h_cosTheta_vs_p_mc", "cos(#Theta) vs. P_{L} for MC; P_{L} / GeV/c ; cos(#Theta)", 100,-2,7,100,-1,1);
  ntpMC->Project("h_cosTheta_vs_p_mc", "cos(tht):pz", "moth==0");
  
  TH2F * h_pt_vs_pl_mc = new TH2F("h_pt_vs_pl_mc", "P_{t} vs.P_{l} for MC; P_{l}/GeV/c; P_{t}/GeV/c", 100,-1,7,100,0,2);
  ntpMC->Project("h_pt_vs_pl_mc", "sqrt(px*px+py*py):pz", "moth==0");



  //projections for PiPlus

  TH2F * h_pt_vs_pl_piplus = new TH2F("h_pt_vs_pl_piplus", "P_{t} vs. P_{l} for #pi^{+}; P_{l}; P_{t}", 100,-1,7,100,0,2);
  ntpPiPlus->Project("h_pt_vs_pl_piplus", "PiPlus_pt:PiPlus_pz", "MCTruthMatch==1&&Mother==88888");

  TH2F * h_cosTheta_vs_p_piplus = new TH2F("h_cosTheta_vs_p_piplus", "cos(#Theta) vs. P_{L} for #pi^{+}; P_{L} / GeV/c ; cos(#Theta)", 100,-2,7,100,-1.1,1.1);
  ntpPiPlus->Project("h_cosTheta_vs_p_piplus",  "PiPlus_CosTheta:PiPlus_pz", "MCTruthMatch==1&&Mother==88888");

  TH2F * h_theta_vs_p_piplus = new TH2F("h_theta_vs_p_piplus", "#Theta vs. P for #pi^{+}; P / GeV/c ; #Theta", 100,0,10,100,0,3);
  ntpPiPlus->Project("h_theta_vs_p_piplus", "PiPlus_tht:PiPlus_p", "MCTruthMatch==1&&Mother==88888");

  TH2F * h_pt_vs_pl_piplus_McTruth = new TH2F("h_pt_vs_pl_piplus_McTruth", "P_{t} vs. P_{l} for #pi^{+} (McTruth); P_{l}; P_{t}", 100,-1,7,100,0,2);
  ntpPiPlus->Project("h_pt_vs_pl_piplus_McTruth", "PiPlus_MC_pt:PiPlus_MC_pz", "MCTruthMatch==1&&Mother==88888");

  TH2F * h_cosTheta_vs_p_piplus_McTruth = new TH2F("h_cosTheta_vs_p_piplus_McTruth", "cos(#Theta) vs. P_{L} for #pi^{+} (McTruth); P_{L} / GeV/c ; cos(#Theta)", 100,-2,7,100,-1.1,1.1);
  ntpPiPlus->Project("h_cosTheta_vs_p_piplus_McTruth",  "PiPlus_MC_CosTheta:PiPlus_MC_pz", "MCTruthMatch==1&&Mother==88888");

  TH2F * h_theta_vs_p_piplus_McTruth = new TH2F("h_theta_vs_p_piplus_McTruth", "#Theta vs. P for #pi^{+} (McTruth); P / GeV/c ; #Theta", 100,0,10,100,0,3);
  ntpPiPlus->Project("h_theta_vs_p_piplus_McTruth", "PiPlus_MC_tht:PiPlus_MC_p", "MCTruthMatch==1&&Mother==88888");

  TH1F * h_reco_PiPlus = new TH1F("h_reco_PiPlus", "reconstructed; P(MC)/GeV/c; Counts", 100,0,7);
  ntpPiPlus->Project("h_reco_PiPlus", "PiPlus_MC_p", "MCTruthMatch==1&&Mother==88888");

  TH1F * h_MC_PiPlus = new TH1F("h_MC_PiPlus", "MC; P(MC)/GeV/c; Counts", 100,0,7);  
  ntpMC->Project("h_MC_PiPlus", "p", "moth==0&&pdg==211");

  TH1F * h_prob_PiPlus = (TH1F*) h_reco_PiPlus->Clone("h_prob_PiPlus");
  h_prob_PiPlus->Divide(h_MC_PiPlus);
  h_prob_PiPlus->SetTitle("reconstruction Probability vs. Momentum for #pi^{+} ");	
  h_prob_PiPlus->GetXaxis()->SetTitle("P/GeV/c");
  h_prob_PiPlus->GetYaxis()->SetTitle("reconst. Probability");
  h_prob_PiPlus->SetStats(kFALSE);

  TH1F * h_reco_theta_PiPlus = new TH1F("h_reco_theta_PiPlus", "reconstructed; #theta(MC); Counts", 100,0,3);
  ntpPiPlus->Project("h_reco_theta_PiPlus", "PiPlus_MC_tht", "MCTruthMatch==1&&Mother==88888");

  TH1F * h_MC_theta_PiPlus = new TH1F("h_MC_theta_PiPlus", "MC; #theta(MC); Counts", 100,0,3);  
  ntpMC->Project("h_MC_theta_PiPlus", "tht", "moth==0&&pdg==211");

  TH1F * h_prob_theta_PiPlus = (TH1F*) h_reco_theta_PiPlus->Clone("h_prob_theta_PiPlus");
  h_prob_theta_PiPlus->Divide(h_MC_theta_PiPlus);
  h_prob_theta_PiPlus->SetTitle("reconstruction Probability vs. #theta for #pi^{+} ");
  h_prob_theta_PiPlus->GetXaxis()->SetTitle("#theta");
  h_prob_theta_PiPlus->GetYaxis()->SetTitle("reconst. Probability");
  h_prob_theta_PiPlus->SetStats(kFALSE);
 
	
  //projections for PiMinus
	
  TH2F * h_pt_vs_pl_piminus = new TH2F("h_pt_vs_pl_piminus", "P_{t} vs. P_{l} for #pi^{-}; P_{l}; P_{t}", 100,-1,7,100,0,2);
  ntpPiMinus->Project("h_pt_vs_pl_piminus", "PiMinus_pt:PiMinus_pz", "MCTruthMatch==1&&Mother==88888");
 
  TH2F * h_cosTheta_vs_p_piminus = new TH2F("h_cosTheta_vs_p_piminus", "cos(#Theta) vs. P_{L} for #pi^{-}; P_{L} / GeV/c ; cos(#Theta)", 100,-2,7,100,-1.1,1.1);
  ntpPiMinus->Project("h_cosTheta_vs_p_piminus",  "PiMinus_CosTheta:PiMinus_pz", "MCTruthMatch==1&&Mother==88888");

  TH2F * h_theta_vs_p_piminus = new TH2F("h_theta_vs_p_piminus", "#Theta vs. P for #pi^{-}; P / GeV/c ; #Theta", 100,0,10,100,0,3);
  ntpPiMinus->Project("h_theta_vs_p_piminus", "PiMinus_tht:PiMinus_p", "MCTruthMatch==1&&Mother==88888");

  TH2F * h_pt_vs_pl_piminus_McTruth = new TH2F("h_pt_vs_pl_piminus_McTruth", "P_{t} vs. P_{l} for #pi^{-} (McTruth); P_{l}; P_{t}", 100,-1,7,100,0,2);
  ntpPiMinus->Project("h_pt_vs_pl_piminus_McTruth", "PiMinus_MC_pt:PiMinus_MC_pz", "MCTruthMatch==1&&Mother==88888");
 
  TH2F * h_cosTheta_vs_p_piminus_McTruth = new TH2F("h_cosTheta_vs_p_piminus_McTruth", "cos(#Theta) vs. P_{L} for #pi^{-} (McTruth); P_{L} / GeV/c ; cos(#Theta)", 100,-2,7,100,-1.1,1.1);
  ntpPiMinus->Project("h_cosTheta_vs_p_piminus_McTruth",  "PiMinus_MC_CosTheta:PiMinus_MC_pz", "MCTruthMatch==1&&Mother==88888");

  TH2F * h_theta_vs_p_piminus_McTruth = new TH2F("h_theta_vs_p_piminus_McTruth", "#Theta vs. P for #pi^{-} (McTruth); P / GeV/c ; #Theta", 100,0,10,100,0,3);
  ntpPiMinus->Project("h_theta_vs_p_piminus_McTruth", "PiMinus_MC_tht:PiMinus_MC_p", "MCTruthMatch==1&&Mother==88888");

  TH1F * h_reco_PiMinus = new TH1F("h_reco_PiMinus", "reconstructed; P(MC)/GeV/c; Counts", 100,0,7);
  ntpPiMinus->Project("h_reco_PiMinus", "PiMinus_MC_p", "MCTruthMatch==1&&Mother==88888");

  TH1F * h_MC_PiMinus = new TH1F("h_MC_PiMinus", "MC; P(MC)/GeV/c; Counts", 100,0,7);  
  ntpMC->Project("h_MC_PiMinus", "p", "moth==-1&&pdg==-211");
	
  TH1F * h_prob_PiMinus = (TH1F*) h_reco_PiMinus->Clone("h_prob_PiMinus");
  h_prob_PiMinus->Divide(h_MC_PiMinus);
  h_prob_PiMinus->SetTitle("reconstruction Probability vs. Momentum for #pi^{-} "); 
  h_prob_PiMinus->GetXaxis()->SetTitle("P/GeV/c");
  h_prob_PiMinus->GetYaxis()->SetTitle("reconst. Probability");
  h_prob_PiMinus->SetStats(kFALSE);

  TH1F * h_reco_theta_PiMinus = new TH1F("h_reco_theta_PiMinus", "reconstructed; #theta(MC); Counts", 100,0,3);
  ntpPiMinus->Project("h_reco_theta_PiMinus", "PiMinus_MC_tht", "MCTruthMatch==1&&Mother==88888");

  TH1F * h_MC_theta_PiMinus = new TH1F("h_MC_theta_PiMinus", "MC; #theta(MC); Counts", 100,0,3);  
  ntpMC->Project("h_MC_theta_PiMinus", "tht", "moth==-1&&pdg==-211");

  TH1F * h_prob_theta_PiMinus = (TH1F*) h_reco_theta_PiMinus->Clone("h_prob_theta_PiMinus");
  h_prob_theta_PiMinus->Divide(h_MC_theta_PiMinus);
  h_prob_theta_PiMinus->SetTitle("reconstruction Probability vs. #theta for #pi^{-} ");
  h_prob_theta_PiMinus->GetXaxis()->SetTitle("#theta");
  h_prob_theta_PiMinus->GetYaxis()->SetTitle("reconst. Probability");
  h_prob_theta_PiMinus->SetStats(kFALSE);

	
  //projections for beam
  TH1F * h_fvtx_Chi2 = new TH1F("h_fvtx_Chi2", "#chi^{2} distribution for vertex fit; #Chi^{2}; Counts", 100,0,10);
  ntpBeam->Project("h_fvtx_Chi2", "fvtx_Chi2");

  TH1F * h_fvtx_Prob = new TH1F("h_fvtx_Prob", "#chi^{2} probability for vertex fit; Prob(#Chi^{2}); Counts", 100,0,1);
  ntpBeam->Project("h_fvtx_Prob", "fvtx_Prob");

  TH1F * h_f4c_Chi2 = new TH1F("h_f4c_Chi2", "#chi^{2} distribution for the 4C fit; #Chi^{2}; Counts", 100,0,10);
  ntpBeam->Project("h_f4c_Chi2", "f4c_Chi2");

  TH1F * h_f4c_Prob = new TH1F("h_f4c_Prob", "#chi^{2} probability for the 4C fit; Prob(#Chi^{2}); Counts", 100,0,1);
  ntpBeam->Project("h_f4c_Prob", "f4c_prob");

  TH1F * h_beam_mass = new TH1F("h_beam_mass", "Invariant Mass for #bar{p} and p; mass GeV/c^{2}; Counts", 100,0,6);
  ntpBeam->Project("h_beam_mass", "Beam_m", "McTruth==1");

  TH1F * h_beam_mom = new TH1F("h_beam_mom", "Momentum for #bar{p} and p; P/GeV/c; Counts", 100,0,8);
  ntpBeam->Project("h_beam_mom", "Beam_p", "McTruth==1");

  TH1F * h_beam_energy = new TH1F("h_beam_energy", "Energy for #bar{p} and p; E/GeV; Counts", 100,0,10);
  ntpBeam->Project("h_beam_energy", "Beam_e", "McTruth==1");

  TH1F * h_beam_pl = new TH1F("h_beam_pl", "P_{L} and p; P_{L}/GeV/c; Counts", 100,0,8);
  ntpBeam->Project("h_beam_pl", "Beam_pz", "McTruth==1");

  TH1F * h_beam_vertex_x = new TH1F("h_beam_vertex_x", "x-Position of vertex for #bar{p} and p; x/cm; Counts", 100,-5,5);
  ntpBeam->Project("h_beam_vertex_x", "Beam_vx", "McTruth==1");

  TH1F * h_beam_vertex_y = new TH1F("h_beam_vertex_y", "y-Position of vertex for #bar{p} and p; y/cm; Counts", 100,-5,5);
  ntpBeam->Project("h_beam_vertex_y", "Beam_vy", "McTruth==1");

  TH1F * h_beam_vertex_z = new TH1F("h_beam_vertex_z", "z-Position of vertex for #bar{p} and p; z/cm; Counts", 100,-5,5);
  ntpBeam->Project("h_beam_vertex_z", "Beam_vz", "McTruth==1");

  //Create Histograms
  TCanvas * c_theta_vs_p_mc = new TCanvas("c_theta_vs_p_mc", "theta_vs_p_mc", 0,0,800,500);
  h_theta_vs_p_mc->Draw("colz");

  TCanvas * c_costheta_vs_p_mc = new TCanvas("c_costheta_vs_p_mc", "costheta_vs_p_mc", 0,0,800,500);
  h_cosTheta_vs_p_mc->Draw("colz");

  TCanvas * c_pt_vs_pl_mc = new TCanvas("c_pt_vs_pl_mc", "P_{t} vs P_{l} mc", 0,0,800,500);
  h_pt_vs_pl_mc->Draw("colz");

  
  
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

  TCanvas * c_prob_PiPlus = new TCanvas ("c_prob_PiPlus", "prob_vs_p_PiPlus", 0,0,800,500);
  h_prob_PiPlus->Draw(); 

  TCanvas * c_prob_theta_PiPlus = new TCanvas ("c_prob_theta_PiPlus", "prob_vs_costheta_PiPlus", 0,0,800,500);
  h_prob_theta_PiPlus->Draw(); 


  
  
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
	*/
  TCanvas * c_prob_PiMinus = new TCanvas ("c_prob_PiMinus", "prob_vs_p_PiMinus", 0,0,800,500);
  h_prob_PiMinus->Draw(); 

  TCanvas * c_prob_theta_PiMinus = new TCanvas ("c_prob_theta_PiMinus", "prob_vs_costheta_PiMinus", 0,0,800,500);
  h_prob_theta_PiMinus->Draw(); 
	


  
  TCanvas * c_fvtx_Chi2 = new TCanvas("c_fvxt_Chi2", "fvxt_Chi2", 0,0,800,500);
  h_fvtx_Chi2->Draw();
  
  TCanvas * c_fvtx_Prob = new TCanvas("c_fvxt_Prob", "fvxt_Prob", 0,0,800,500);
  h_fvtx_Prob->Draw();
  
  TCanvas * c_f4c_Chi2 = new TCanvas("c_f4c_Chi2", "f4c_Chi2", 0,0,800,500);
  h_f4c_Chi2->Draw();

  TCanvas * c_f4c_Prob = new TCanvas("c_f4c_Prob", "f4c_Prob", 0,0,800,500);
  h_f4c_Prob->Draw();

  TCanvas * c_beam_mass = new TCanvas("c_beam_mass", "beam_mass", 0,0,800,500);
  h_beam_mass->Draw();

  TCanvas * c_beam_mom = new TCanvas("c_beam_mom", "beam_mom", 0,0,800,500);
  h_beam_mom->Draw();

  TCanvas * c_beam_energy = new TCanvas("c_beam_energy", "beam_energy", 0,0,800,500);
  h_beam_energy->Draw();

  TCanvas * c_beam_pl = new TCanvas("c_beam_pl", "beam_pl", 0,0,800,500);
  h_beam_pl->Draw();

  TCanvas * c_beam_vertex_x = new TCanvas("c_beam_vertex_x", "beam_vertex_x", 0,0,800,500);
  h_beam_vertex_x->Draw();
 
  TCanvas * c_beam_vertex_y = new TCanvas("c_beam_vertex_y", "beam_vertex_y", 0,0,800,500);
  h_beam_vertex_y->Draw();

  TCanvas * c_beam_vertex_z = new TCanvas("c_beam_vertex_z", "beam_vertex_z", 0,0,800,500);
  h_beam_vertex_z->Draw();


	
  //***Save output if saveOutput is true	
  if (saveOutput==true) {
		
    c_theta_vs_p_mc->Print(inputDir + "plots/theta_vs_p_mc.root");
    c_costheta_vs_p_mc->Print(inputDir + "plots/costheta_vs_p_mc.root");
    c_pt_vs_pl_mc->Print(inputDir + "plots/pt_vs_pl_mc.root");


    c_theta_vs_p_piplus->Print(inputDir + "plots/theta_vs_p_piplus.root");
    c_costheta_vs_p_piplus->Print(inputDir + "plots/costheta_vs_p_piplus.root");
    c_pt_vs_pl_piplus->Print(inputDir + "plots/pt_vs_pl_piplus.root");
    c_theta_vs_p_piplus_McTruth->Print(inputDir + "plots/theta_vs_p_piplus_McTruth.root");
    c_costheta_vs_p_piplus_McTruth->Print(inputDir + "plots/costheta_vs_p_piplus_McTruth.root");
    c_pt_vs_pl_piplus_McTruth->Print(inputDir + "plots/pt_vs_pl_piplus_McTruth.root");
		c_prob_PiPlus->Print(inputDir+ "plots/PiPlus_reco_Prob_vs_P.root");
		c_prob_theta_PiPlus->Print(inputDir+ "plots/PiPlus_reco_Prob_vs_costheta.root");


    c_theta_vs_p_piminus->Print(inputDir + "plots/theta_vs_p_piminus.root");
    c_costheta_vs_p_piminus->Print(inputDir + "plots/costheta_vs_p_piminus.root");
    c_pt_vs_pl_piminus->Print(inputDir + "plots/pt_vs_pl_piminus.root");
    c_theta_vs_p_piminus_McTruth->Print(inputDir + "plots/theta_vs_p_piminus_McTruth.root");
    c_costheta_vs_p_piminus_McTruth->Print(inputDir + "plots/costheta_vs_p_piminus_McTruth.root");
    c_pt_vs_pl_piminus_McTruth->Print(inputDir + "plots/pt_vs_pl_piminus_McTruth.root");
		c_prob_PiMinus->Print(inputDir+ "plots/PiMinus_reco_Prob_vs_P.root");
		c_prob_theta_PiMinus->Print(inputDir+ "plots/PiMinus_reco_Prob_vs_costheta.root");


		c_fvtx_Chi2->Print(inputDir + "plots/fvtx_Chi2.root");	
		c_fvtx_Prob->Print(inputDir + "plots/fvtx_Prob.root");	
		c_f4c_Chi2->Print(inputDir + "plots/f4c_Chi2.root");	
		c_f4c_Prob->Print(inputDir + "plots/f4c_Prob.root");	
		c_beam_mass->Print(inputDir + "plots/beam_mass.root");	
		c_beam_mom->Print(inputDir + "plots/beam_mom.root");	
		c_beam_energy->Print(inputDir + "plots/beam_energy.root");	
		c_beam_pl->Print(inputDir + "plots/beam_pl.root");	
		c_beam_vertex_x->Print(inputDir + "plots/beam_vertex_x.root");	
		c_beam_vertex_y->Print(inputDir + "plots/beam_vertex_y.root");
		c_beam_vertex_z->Print(inputDir + "plots/beam_vertex_z.root");	

	
  }


//  ***close histograms
  
		c_theta_vs_p_mc->Close();
		c_costheta_vs_p_mc->Close();
		c_pt_vs_pl_mc->Close();
		
		c_pt_vs_pl_piplus->Close();
		c_costheta_vs_p_piplus->Close();
		c_theta_vs_p_piplus->Close();
 		c_pt_vs_pl_piplus_McTruth->Close();
		c_costheta_vs_p_piplus_McTruth->Close();
		c_theta_vs_p_piplus_McTruth->Close();
 		c_prob_PiPlus->Close();
 		c_prob_theta_PiPlus->Close();
	
			
		c_pt_vs_pl_piminus->Close();
		c_costheta_vs_p_piminus->Close();
		c_theta_vs_p_piminus->Close();
		c_pt_vs_pl_piminus_McTruth->Close();
		c_costheta_vs_p_piminus_McTruth->Close();
		c_theta_vs_p_piminus_McTruth->Close();
		c_prob_PiMinus->Close();
		c_prob_theta_PiMinus->Close();

		 
		c_fvtx_Chi2->Close();
		c_fvtx_Prob->Close();
		c_f4c_Chi2->Close();
		c_f4c_Prob->Close();
		c_beam_pl->Close();
		c_beam_vertex_x->Close();
		c_beam_vertex_y->Close();
		c_beam_vertex_z->Close();
		c_beam_mass->Close();
		c_beam_mom->Close();
		c_beam_energy->Close();

  
}

int main(){
  evaluation();
  return 0;
}

