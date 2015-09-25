/**
 * @file eval_final_states.C
 * @mainpage eval_final_states.C Get Data from analysis file and creates and saves different histograms
 *
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Get Data from analysis file and create and save different histograms
 * @details This files get the data for Xi(1820) particle from the analysis of
 * pbar p -> Xi+ Xi(1820)-
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
#include "../common_jenny.cpp"


void eval_Xi1820(TString path="", bool save=kTRUE, bool close=kFALSE){

	//*** Input file
	TString inFile = TString::Format("%s/output_ana.root", path.Data());


	//*** get Data from Tree
	TFile * data = new TFile(inFile, "READ");

	TTree * ntpXiMinus1820 = (TTree*) data->Get("ntpXiMinus1820");
	TTree * ntpXiPlus1820 = (TTree*) data->Get("ntpXiPlus1820");

	//*******************************************************************************************************************

	gStyle->SetOptStat(1111);
	gStyle->SetStatX(0.85);
	gStyle->SetStatY(0.85);

	TString vtxcut = "&& VtxFit_HowGood==1 ";
	TString masscut = "&& MassFit_HowGood>0 ";

	//**** Get information about XiMinus ********************************************************************************
	if(ntpXiMinus1820!=0X0){
		TH2D * h_xi_pt_vs_pz = new TH2D("h_xi_pt_vs_pz", "Transversal vs. longitudinal momentum for #Xi(1820)^{-}; p_{z}/GeV/c; p_{t}/GeV/c", 40,-0.5,3.9,40,0,0.8);
		ntpXiMinus1820->Project("h_xi_pt_vs_pz", "MCTruth_pt: MCTruth_pz", "McTruthMatch");
		jenny::CreateDrawAndSaveHistogram(h_xi_pt_vs_pz, path+"/plots/", "XiMinus_pt_vs_pz", save, close);

		gStyle->SetOptStat(0);

		TH1D * h_xi_m_nocut = new TH1D("h_xi_m_nocut", "Mass distribution for #Xi(1820)^{-}; m/GeV/c^{2}; counts", 100,1.5,2);
		ntpXiMinus1820->Project("h_xi_m_nocut", "XiMinus_m", "McTruthMatch");

		TH1D * h_xi_m_vtxcut = new TH1D("h_xi_m_vtxcut", "Mass distribution for #Xi(1820)^{-} with vtxcut; m/GeV/c^{2}; counts", 100,1.5,2);
		ntpXiMinus1820->Project("h_xi_m_vtxcut", "XiMinus_m", "McTruthMatch "+vtxcut);

		TH1D * h_xi_m_masscut = new TH1D("h_xi_m_masscut", "Mass distribution for #Xi(1820)^{-} with vertex cut and mass cut; m/GeV/c^{2}; counts", 100,1.5,2);
		ntpXiMinus1820->Project("h_xi_m_masscut", "XiMinus_m", "McTruthMatch "+vtxcut+masscut);

		jenny::CreateDrawAndSaveNHistograms(h_xi_m_nocut, h_xi_m_vtxcut, h_xi_m_masscut, "no cut", "VtxFit_prob>0.01", "VtxFit_prob>0.01 && MassFit_prob>0.01", path+"/plots/", "XiMinus_m_diffcuts", save, close);

		gStyle->SetOptStat(1111);

		TH1D * h_xi_m_masscut2 = new TH1D("h_xi_m_masscut2", "Mass distribution for #Xi(1820)^{-} with vertex cut and mass cut; m/GeV/c^{2}; counts", 100,1.5,2);
		ntpXiMinus1820->Project("h_xi_m_masscut2", "XiMinus_m", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramFit(h_xi_m_masscut2, path+"/plots/", "XiMinus_m_masscut2", save, close);

		TH1D * h_xi_vtxres_x = new TH1D("h_xi_vtxres_x", "resolution for x coordinate of vertex for #Xi(1820)^{-}; x-x_{MC}; counts", 100,-1,1);
		ntpXiMinus1820->Project("h_xi_vtxres_x", "VtxFit_diffvx", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_xi_vtxres_x, path+"/plots/", "XiMinus_vtxres_x", save, close, false, 0.1, 1, true);

		TH1D * h_xi_vtxres_y = new TH1D("h_xi_vtxres_y", "resolution for y coordinate of vertex for #Xi(1820)^{-}; y-y_{MC}; counts", 100,-1,1);
		ntpXiMinus1820->Project("h_xi_vtxres_y", "VtxFit_diffvy", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_xi_vtxres_y, path+"/plots/", "XiMinus_vtxres_y", save, close, false, 0.1, 1, true);

		TH1D * h_xi_vtxres_z = new TH1D("h_xi_vtxres_z", "resolution for z coordinate of vertex for #Xi(1820)^{-}; z-z_{MC}; counts", 100,-1,1);
		ntpXiMinus1820->Project("h_xi_vtxres_z", "VtxFit_diffvz", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_xi_vtxres_z, path+"/plots/", "XiMinus_vtxres_z", save, close, false, 0.03, 1);//, true);

		TH1D * h_xi_costht = new TH1D("h_xi_costht", "cos(#Theta) distribution for #Xi(1820)^{-}; cos(#Theta); counts", 100,-1,1);
		ntpXiMinus1820->Project("h_xi_costht", "cos(VtxFit_tht)","McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogram(h_xi_costht, path+"/plots/", "XiMinus_costht", save, close);

		TH1D * h_xi_tht = new TH1D("h_xi_tht", "#Theta distribution for #Xi(1820)^{-}; #Theta/rad; counts", 100,-1,1);
		ntpXiMinus1820->Project("h_xi_tht", "VtxFit_tht","McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogram(h_xi_tht, path+"/plots/", "XiMinus_tht", save, close);
	}




	else if (ntpXiPlus1820 !=0x0){

		//**** Get information about XiPlus ****************************************************************************

		TH2D * h_axi_pt_vs_pz = new TH2D("h_axi_pt_vs_pz", "Transversal vs. longitudinal momentum for #bar{#Xi(1820)}; p_{z}/GeV/c; p_{t}/GeV/c", 40,0.5,3.1,40,0,0.7);
		ntpXiPlus1820->Project("h_axi_pt_vs_pz", "MCTruth_pt: MCTruth_pz", "McTruthMatch");
		jenny::CreateDrawAndSaveHistogram(h_axi_pt_vs_pz, path+"/plots/", "XiPlus_pt_vs_pz", save, close);

		gStyle->SetOptStat(0);
		TH1D * h_axi_m_nocut = new TH1D("h_axi_m_nocut", "Mass distribution for #bar{#Xi(1820)}; m/GeV/c^{2}; counts", 100,1.5,2);
		ntpXiPlus1820->Project("h_axi_m_nocut", "Xiplus_m", "McTruthMatch");

		TH1D * h_axi_m_vtxcut = new TH1D("h_axi_m_vtxcut", "Mass distribution for #bar{#Xi(1820)} with vtxcut; m/GeV/c^{2}; counts", 100,1.5,2);
		ntpXiPlus1820->Project("h_axi_m_vtxcut", "Xiplus_m", "McTruthMatch "+vtxcut);

		TH1D * h_axi_m_masscut = new TH1D("h_axi_m_masscut", "Mass distribution for #bar{#Xi(1820)} with vertex cut and mass cut; m/GeV/c^{2}; counts", 100,1.5,2);
		ntpXiPlus1820->Project("h_axi_m_masscut", "Xiplus_m", "McTruthMatch "+vtxcut+masscut);

		jenny::CreateDrawAndSaveNHistograms(h_axi_m_nocut, h_axi_m_vtxcut, h_axi_m_masscut, "no cut", "VtxFit_prob>0.01", "VtxFit_prob>0.01 && MassFit_prob>0.01", path+"/plots/", "XiPlus_m_diffcuts", save, close);

		gStyle->SetOptStat(1111);

		TH1D * h_axi_m_masscut2 = new TH1D("h_axi_m_masscut2", "Mass distribution for #bar{#Xi(1820)} with vertex cut and mass cut; m/GeV/c^{2}; counts", 100,1.5,2);
		ntpXiPlus1820->Project("h_axi_m_masscut2", "Xiplus_m", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramFit(h_axi_m_masscut2, path+"/plots/", "XiPlus_m_masscut", save, close);

		TH1D * h_axi_vtxres_x = new TH1D("h_axi_vtxres_x", "resolution for x coordinate of vertex for #bar{#Xi(1820)}; x-x_{MC}; counts", 100,-1,1);
		ntpXiPlus1820->Project("h_axi_vtxres_x", "VtxFit_diffvx", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_axi_vtxres_x, path+"/plots/", "XiPlus_vtxres_x", save, close, false, 0.1, 1, true);

		TH1D * h_axi_vtxres_y = new TH1D("h_axi_vtxres_y", "resolution for y coordinate of vertex for #bar{#Xi(1820)}; y-y_{MC}; counts", 100,-1,1);
		ntpXiPlus1820->Project("h_axi_vtxres_y", "VtxFit_diffvy", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_axi_vtxres_y, path+"/plots/", "XiPlus_vtxres_y", save, close, false, 0.1, 1, true);

		TH1D * h_axi_vtxres_z = new TH1D("h_axi_vtxres_z", "resolution for z coordinate of vertex for #bar{#Xi(1820)}; z-z_{MC}; counts", 100,-1,1);
		ntpXiPlus1820->Project("h_axi_vtxres_z", "VtxFit_diffvz", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_axi_vtxres_z, path+"/plots/", "XiPlus_vtxres_z", save, close, false, 0.05,1);//, true);


		TH1D * h_axi_costht = new TH1D("h_axi_costht", "cos(#Theta) distribution for #bar{#Xi(1820)}; cos(#Theta); counts", 100,-1,1);
		ntpXiPlus1820->Project("h_axi_costht", "cos(VtxFit_tht)","McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogram(h_axi_costht, path+"/plots/", "XiPlus_costht", save, close);

		TH1D * h_axi_tht = new TH1D("h_axi_tht", "#Theta distribution for #bar{#Xi(1820)}; #Theta/rad; counts", 100,-1,1);
		ntpXiPlus1820->Project("h_axi_tht", "VtxFit_tht","McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogram(h_axi_tht, path+"/plots/", "XiPlus_tht", save, close);
	}
	else {
		cout << "No particle of kind Xi(1820)- or Xi(1820)+!" << endl;
	}
}
