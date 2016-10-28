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


void eval_Xi1690(TString prefix="", bool save=kTRUE, bool close=kFALSE){

	//*** Input file
	if(prefix==""){
		TString inFile ="output_ana.root";
	}
	else{
		TString inFile = TString::Format("%s/output_ana.root", prefix.Data());
	}


	//*** get Data from Tree
	TFile * data = new TFile(inFile, "READ");

	TTree * ntpXiMinus1690 = (TTree*) data->Get("ntpXiMinus1690");
	TTree * ntpXiPlus1690 = (TTree*) data->Get("ntpXiPlus1690");

	//*******************************************************************************************************************


	TString vtxcut = "&& HitTag==1 && VtxFit_HowGood==1 ";


	//**** Get information about XiMinus ********************************************************************************
	if(ntpXiMinus1690!=0X0){
		TH2D * h_xi_pt_vs_pz = new TH2D("h_xi_pt_vs_pz", "Transverse vs. longitudinal momentum for #Xi^{-}(1690); p_{z} [GeV/c]; p_{t} [GeV/c]", 200,1,3.9,200,0,0.8);
		ntpXiMinus1690->Project("h_xi_pt_vs_pz", "MCTruth_pt: MCTruth_pz", "McTruthMatch");
		jenny::CreateDrawAndSaveHistogram(h_xi_pt_vs_pz, prefix+"/plots/", "XiMinus1690_pt_vs_pz", save, close);

		TH2D * h_xi_pt_vs_pz_cut = new TH2D("h_xi_pt_vs_pz_cut", "Transverse vs. longitudinal momentum for #Xi^{-}(1690); p_{z} [GeV/c]; p_{t} [GeV/c]", 200,1,3.9,200,0,0.8);
		ntpXiMinus1690->Project("h_xi_pt_vs_pz_cut", "MCTruth_pt: MCTruth_pz", "McTruthMatch"+vtxcut);
		jenny::CreateDrawAndSaveHistogram(h_xi_pt_vs_pz_cut, prefix+"/plots/", "XiMinus1690_pt_vs_pz_cut", save, close);


		TH1D * h_xi_m_nocut = new TH1D("h_xi_m_nocut", "Mass distribution for #Xi^{-}(1690); M [GeV/c^{2}]; counts",  500,1.49,1.89);
		ntpXiMinus1690->Project("h_xi_m_nocut", "VtxFit_m", "McTruthMatch");

		TH1D * h_xi_m_vtxcut = new TH1D("h_xi_m_vtxcut", "Mass distribution for #Xi^{-}(1690) with vtxcut; M [GeV/c^{2}]; counts", 500,1.49,1.89);
		ntpXiMinus1690->Project("h_xi_m_vtxcut", "VtxFit_m", "McTruthMatch "+vtxcut);

		jenny::CreateDrawAndSaveNHistograms(h_xi_m_nocut, h_xi_m_vtxcut, "mass window", "VtxFit_prob>0.01", prefix+"/plots/", "XiMinus1690_m_diffcuts", save, close);


		TH1D * h_xi_m = new TH1D("h_xi_m", "Mass distribution for #Xi^{-}(1690) with vertex cut and mass cut; M [GeV/c^{2}]; counts", 500,1.49,1.89);
		ntpXiMinus1690->Project("h_xi_m", "VtxFit_m", "McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogramVoigtFit(h_xi_m, prefix+"/plots/", "XiMinus1690_m_masscut", save, close, 1.6311 ,1.718, 1.690,0.029);

		TH1D * h_xi_m_mc = new TH1D("h_xi_m_mc", "Mass distribution for #Xi^{-}(1690) to MC truth; (M_{#Xi*}-M_{PDG}) [GeV/c^{2}]; counts", 500,-50,50);
		ntpXiMinus1690->Project("h_xi_m_mc", "(VtxFit_m-MCTruth_m)*1000", "McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_xi_m_mc, prefix+"/plots/", "XiMinus1690_m_masscut_diff_mc", save, close);


		TH1D * h_xi_vtxres_x = new TH1D("h_xi_vtxres_x", "resolution for x coordinate of vertex for #Xi^{-}(1690); x-x_{MC} [mm]; counts", 500,-2,2);
		ntpXiMinus1690->Project("h_xi_vtxres_x", "VtxFit_diffvx*10", "McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_xi_vtxres_x, prefix+"/plots/", "XiMinus1690_vtxres_x", save, close);

		TH1D * h_xi_vtxres_y = new TH1D("h_xi_vtxres_y", "resolution for y coordinate of vertex for #Xi^{-}(1690); y-y_{MC} [mm]; counts", 500,-2,2);
		ntpXiMinus1690->Project("h_xi_vtxres_y", "VtxFit_diffvy*10", "McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_xi_vtxres_y, prefix+"/plots/", "XiMinus1690_vtxres_y", save, close);

		TH1D * h_xi_vtxres_z = new TH1D("h_xi_vtxres_z", "resolution for z coordinate of vertex for #Xi^{-}(1690); z-z_{MC} [mm]; counts", 500,-5,5);
		ntpXiMinus1690->Project("h_xi_vtxres_z", "VtxFit_diffvz*10", "McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_xi_vtxres_z, prefix+"/plots/", "XiMinus1690_vtxres_z", save, close);

		TH1D * h_xi_costht = new TH1D("h_xi_costht", "cos(#Theta) distribution for #Xi^{-}(1690); cos(#Theta); counts", 500,-1,1.01);
		ntpXiMinus1690->Project("h_xi_costht", "cos(VtxFit_tht)","McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogram(h_xi_costht, prefix+"/plots/", "XiMinus1690_costht", save, close);

		TH1D * h_xi_tht = new TH1D("h_xi_tht", "#Theta distribution for #Xi^{-}(1690); #Theta [rad]; counts", 500,0,0.5);
		ntpXiMinus1690->Project("h_xi_tht", "VtxFit_tht","McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogram(h_xi_tht, prefix+"/plots/", "XiMinus1690_tht", save, close);

		TH1D * h_xi_chisq = new TH1D("h_xi_chisq", "#chi^{2} distribution for #Xi^{-}(1690); #chi^{2}; counts", 500,0,100);
		ntpXiMinus1690->Project("h_xi_chisq", "VtxFit_chisq","McTruthMatch ");
		jenny::CreateDrawAndSaveHistogram(h_xi_chisq, prefix+"/plots/", "XiMinus1690_chisq", save, close, true);

		TH1D * h_xi_prob = new TH1D("h_xi_prob", "probability distribution for #Xi^{-}(1690); prob; counts", 500,0,1);
		ntpXiMinus1690->Project("h_xi_prob", "VtxFit_prob","McTruthMatch ");
		jenny::CreateDrawAndSaveHistogram(h_xi_prob, prefix+"/plots/", "XiMinus1690_prob", save, close,true);
	}




	else if (ntpXiPlus1690 !=0x0){

		//**** Get information about XiPlus ****************************************************************************

		TH2D * h_axi_pt_vs_pz = new TH2D("h_axi_pt_vs_pz", "Transverse vs. longitudinal momentum for #bar{#Xi}^{+}(1690); p_{z} [GeV/c]; p_{t} [GeV/c]", 200,1,3.9,200,0,0.7);
		ntpXiPlus1690->Project("h_axi_pt_vs_pz", "MCTruth_pt: MCTruth_pz", "McTruthMatch");
		jenny::CreateDrawAndSaveHistogram(h_axi_pt_vs_pz, prefix+"/plots/", "XiPlus1690_pt_vs_pz", save, close);

		TH2D * h_axi_pt_vs_pz_cut = new TH2D("h_axi_pt_vs_pz_cut", "Transverse vs. longitudinal momentum for #bar{#Xi}^{+}(1690); p_{z} [GeV/c]; p_{t} [GeV/c]", 200,1,3.9,200,0,0.7);
		ntpXiPlus1690->Project("h_axi_pt_vs_pz_cut", "MCTruth_pt: MCTruth_pz", "McTruthMatch"+vtxcut);
		jenny::CreateDrawAndSaveHistogram(h_axi_pt_vs_pz_cut, prefix+"/plots/", "XiPlus1690_pt_vs_pz_cut", save, close);

		TH1D * h_axi_m_nocut = new TH1D("h_axi_m_nocut", "Mass distribution for #bar{#Xi}^{+}(1690); M [GeV/c^{2}]; counts", 500,1.49,1.89);
		ntpXiPlus1690->Project("h_axi_m_nocut", "VtxFit_m", "McTruthMatch");

		TH1D * h_axi_m_vtxcut = new TH1D("h_axi_m_vtxcut", "Mass distribution for #bar{#Xi}^{+}(1690) with vtxcut; M [GeV/c^{2}]; counts", 500,1.49,1.89);
		ntpXiPlus1690->Project("h_axi_m_vtxcut", "VtxFit_m", "McTruthMatch "+vtxcut);

		jenny::CreateDrawAndSaveNHistograms(h_axi_m_nocut, h_axi_m_vtxcut, "mass window", "VtxFit_prob>0.01", prefix+"/plots/", "XiPlus1690_m_diffcuts", save, close);


		TH1D * h_axi_m = new TH1D("h_axi_m", "Mass distribution for #bar{#Xi}^{+}(1690) with vertex cut and mass cut; M [GeV/c^{2}]; counts", 500,1.49,1.89);
		ntpXiPlus1690->Project("h_axi_m", "VtxFit_m", "McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogramVoigtFit(h_axi_m, prefix+"/plots/", "XiPlus1690_m_masscut", save, close, 1.6311 ,1.719, 1.690,0.029);

		TH1D * h_axi_m_mc = new TH1D("h_axi_m_mc", "Mass distribution for #bar{#Xi}^{+}(1690) to MC truth; (M_{#bar{#Xi}*}-M_{PDG}) [GeV/c^{2}]; counts", 500,-50,50);
		ntpXiPlus1690->Project("h_axi_m_mc", "(VtxFit_m-MCTruth_m)*1000", "McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_axi_m_mc, prefix+"/plots/", "XiPlus1690_m_masscut", save, close);

		TH1D * h_axi_vtxres_x = new TH1D("h_axi_vtxres_x", "resolution for x coordinate of vertex for #bar{#Xi}^{+}(1690); x-x_{MC} [mm]; counts", 500,-1,1);
		ntpXiPlus1690->Project("h_axi_vtxres_x", "VtxFit_diffvx*10", "McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_axi_vtxres_x, prefix+"/plots/", "XiPlus1690_vtxres_x", save, close);

		TH1D * h_axi_vtxres_y = new TH1D("h_axi_vtxres_y", "resolution for y coordinate of vertex for #bar{#Xi}^{+}(1690); y-y_{MC} [mm]; counts", 500,-1,1);
		ntpXiPlus1690->Project("h_axi_vtxres_y", "VtxFit_diffvy*10", "McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_axi_vtxres_y, prefix+"/plots/", "XiPlus1690_vtxres_y", save, close);

		TH1D * h_axi_vtxres_z = new TH1D("h_axi_vtxres_z", "resolution for z coordinate of vertex for #bar{#Xi}^{+}(1690); z-z_{MC} [mm]; counts", 500,-5,5);
		ntpXiPlus1690->Project("h_axi_vtxres_z", "VtxFit_diffvz*10", "McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_axi_vtxres_z, prefix+"/plots/", "XiPlus1690_vtxres_z", save, close);


		TH1D * h_axi_costht = new TH1D("h_axi_costht", "cos(#Theta) distribution for #bar{#Xi}^{+}(1690); cos(#Theta); counts", 500,-1,1.01);
		ntpXiPlus1690->Project("h_axi_costht", "cos(VtxFit_tht)","McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogram(h_axi_costht, prefix+"/plots/", "XiPlus1690_costht", save, close);

		TH1D * h_axi_tht = new TH1D("h_axi_tht", "#Theta distribution for #bar{#Xi}^{+}(1690); #Theta [rad]; counts", 500,0,0.5);
		ntpXiPlus1690->Project("h_axi_tht", "VtxFit_tht","McTruthMatch "+vtxcut);
		jenny::CreateDrawAndSaveHistogram(h_axi_tht, prefix+"/plots/", "XiPlus1690_tht", save, close);

		TH1D * h_xi_chisq = new TH1D("h_xi_chisq", "#chi^{2} distribution for #bar{#Xi}^{+}(1690); #chi^{2}; counts", 500,0,100);
		ntpXiPlus1690->Project("h_xi_chisq", "VtxFit_chisq","McTruthMatch ");
		jenny::CreateDrawAndSaveHistogram(h_xi_chisq, prefix+"/plots/", "XiPlus1690_chisq", save, close, true);

		TH1D * h_xi_prob = new TH1D("h_xi_prob", "probability distribution for #bar{#Xi}^{+}(1690); prob; counts", 500,0,1);
		ntpXiPlus1690->Project("h_xi_prob", "VtxFit_prob","McTruthMatch ");
		jenny::CreateDrawAndSaveHistogram(h_xi_prob, prefix+"/plots/", "XiPlus1690_prob", save, close, true);
	}
	else {
		cout << "No particle of kind Xi(1690)- or Xi(1690)+!" << endl;
	}
}
