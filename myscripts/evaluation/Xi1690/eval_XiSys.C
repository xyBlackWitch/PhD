/**
 * @file eval_final_states.C
 * @mainpage eval_final_states.C Get Data from analysis file and creates and saves different histograms
 *
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Get Data from analysis file and create and save different histograms
 * @details This files get the data for Xi(1690) particle from the analysis of
 * pbar p -> Xi+ Xi(1690)-
 * 			 |	  |-> XiMinus + K-
 * 			 |			|-> p + Pi-
 * 			 |
 * 			 |->  XiPlus + Pi+
 * 			 		|-> pbar + Pi+
 * 	 and their c.c
 * 	 Different useful histogramms are created and saved.
 */

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "/home/ikp1/puetz/panda/myscripts/common_jenny.cpp"


void eval_XiSys(TString prefix="", bool save=kTRUE, bool close=kFALSE){

	//*** Input file
	if(prefix==""){
		TString inFile ="output_ana.root";
	}
	else{
		TString inFile = TString::Format("%s/output_ana.root", prefix.Data());
	}


	//*** get Data from Tree
	TFile * data = new TFile(inFile, "READ");

	TTree * ntpXiSys = (TTree*) data->Get("ntpXiSys");

	//*******************************************************************************************************************


	TString cut = "&& 4CFit_prob>0.01 ";

	//**** Get information about XiMinus ********************************************************************************
	TH2D * h_xi_pt_vs_pz = new TH2D("h_xi_pt_vs_pz", "Transverse vs. longitudinal momentum for #Xi^{-}(1690)#bar{#Xi}^{+}-System; p_{z}/GeV/c; p_{t}/GeV/c", 200,-5,5,200,-5,5);
	ntpXiSys->Project("h_xi_pt_vs_pz", "McTruth_pt: McTruth_pz", "McTruthMatch");
	jenny::CreateDrawAndSaveHistogram(h_xi_pt_vs_pz, prefix+"/plots/", "XiSys_pt_vs_pz", save, close);


	TH1D * h_xi_m_nocut = new TH1D("h_xi_m_nocut", "Mass distribution for #Xi^{-}(1690)#bar{#Xi}^{+}-System; m/GeV/c^{2}; counts", 500,3.2,3.3);
	ntpXiSys->Project("h_xi_m_nocut", "4cFit_m", "McTruthMatch");

	TH1D * h_xi_m_cut = new TH1D("h_xi_m_cut", "Mass distribution for #Xi^{-}(1690)#bar{#Xi}^{+}-System with cut; m/GeV/c^{2}; counts", 500,3.2,3.3);
	ntpXiSys->Project("h_xi_m_cut", "4cFit_m", "McTruthMatch "+cut);

	jenny::CreateDrawAndSaveNHistograms(h_xi_m_nocut, h_xi_m_cut, "no cut", "4cFit_prob>0.01", prefix+"/plots/", "XiSys_m_diffcuts", save, close);


//	TH1D * h_xi_m_masscut2 = new TH1D("h_xi_m_masscut2", "Mass distribution for #Xi^{-}(1690)#bar{#Xi}^{+}-System with 4C-Fit cut; m/GeV/c^{2}; counts", 500,3.2,3.3);
//	ntpXiSys->Project("h_xi_m_masscut2", "4cFit_m", "McTruthMatch "+cut);
//	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_xi_m_masscut2, prefix+"/plots/", "XiSys_m_masscut2", save, close, false, 0.005, 0.05, true);
//

	TH1D * h_xi_vtxres_x = new TH1D("h_xi_vtxres_x", "resolution for x coordinate of vertex for #Xi^{-}(1690)#bar{#Xi}^{+}-System; x-x_{MC}; counts", 500,-0.02,0.02);
	ntpXiSys->Project("h_xi_vtxres_x", "4CFit_diffvx", "McTruthMatch "+cut);
	jenny::CreateDrawAndSaveHistogramFWHM(h_xi_vtxres_x, prefix+"/plots/", "XiSys_vtxres_x", save, close);

	TH1D * h_xi_vtxres_y = new TH1D("h_xi_vtxres_y", "resolution for y coordinate of vertex for #Xi^{-}(1690)#bar{#Xi}^{+}-System; y-y_{MC}; counts", 500,-0.02,0.02);
	ntpXiSys->Project("h_xi_vtxres_y", "4CFit_diffvy", "McTruthMatch "+cut);
	jenny::CreateDrawAndSaveHistogramFWHM(h_xi_vtxres_y, prefix+"/plots/", "XiSys_vtxres_y", save, close);

	TH1D * h_xi_vtxres_z = new TH1D("h_xi_vtxres_z", "resolution for z coordinate of vertex for #Xi^{-}(1690)#bar{#Xi}^{+}-System; z-z_{MC}; counts", 500,-0.02,0.02);
	ntpXiSys->Project("h_xi_vtxres_z", "4CFit_diffvz", "McTruthMatch "+cut);
	jenny::CreateDrawAndSaveHistogramFWHM(h_xi_vtxres_z, prefix+"/plots/", "XiSys_vtxres_z", save, close);

	TH1D * h_xi_costht = new TH1D("h_xi_costht", "cos(#Theta) distribution for #Xi^{-}(1690)#bar{#Xi}^{+}-System; cos(#Theta); counts", 500,-1,1.01);
	ntpXiSys->Project("h_xi_costht", "cos(4cFit_tht)","McTruthMatch "+cut);
	jenny::CreateDrawAndSaveHistogram(h_xi_costht, prefix+"/plots/", "XiSys_costht", save, close);

	TH1D * h_xi_tht = new TH1D("h_xi_tht", "#Theta distribution for #Xi^{-}(1690)#bar{#Xi}^{+}-System; #Theta/rad; counts", 500,-0.5,3);
	ntpXiSys->Project("h_xi_tht", "4cFit_tht","McTruthMatch "+cut);
	jenny::CreateDrawAndSaveHistogram(h_xi_tht, prefix+"/plots/", "XiSys_tht", save, close);

	TH1D * h_xi_chisq = new TH1D("h_xi_chisq", "#chi^{2} distribution for #Xi^{-}(1690)#bar{#Xi}^{+}-System; #chi^{2}; counts", 500,0,100);
	ntpXiSys->Project("h_xi_chisq", "4CFit_chisq","McTruthMatch ");
	jenny::CreateDrawAndSaveHistogram(h_xi_chisq, prefix+"/plots/", "XiSys_chisq", save, close,true);

	TH1D * h_xi_prob = new TH1D("h_xi_prob", "probability distribution for #Xi^{-}(1690)#bar{#Xi}^{+}-System; prob; counts", 500,0,1);
	ntpXiSys->Project("h_xi_prob", "4CFit_prob","McTruthMatch ");
	jenny::CreateDrawAndSaveHistogram(h_xi_prob, prefix+"/plots/", "XiSys_prob", save, close, true);

}
