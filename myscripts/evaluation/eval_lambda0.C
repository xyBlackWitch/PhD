/**
 * @file eval_final_states.C
 * @mainpage eval_final_states.C Get Data from analysis file and creates and saves different histograms
 *
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Get Data from analysis file and create and save different histograms
 * @details This files get the data for lambda0 particle from the analysis of
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
#include "../common_jenny.cpp"


void eval_lambda0(TString path="", bool save=kTRUE, bool close=kFALSE){

	//*** Input file
	TString inFile = TString::Format("%s/output_ana.root", path.Data());


	//*** get Data from Tree
	TFile * data = new TFile(inFile, "READ");

	TTree * ntpLambda0 = (TTree*) data->Get("ntpLambda0");
	TTree * ntpAntiLambda0 = (TTree*) data->Get("ntpAntiLambda0");

	//*******************************************************************************************************************

	gStyle->SetOptStat(1111);
	gStyle->SetStatX(0.85);
	gStyle->SetStatY(0.85);

	TString vtxcut = "&& VtxFit_HowGood==1 ";
	TString masscut = "&& MassFit_HowGood>0 ";

	//**** Get information about Lambda0 ********************************************************************************

	TH2D * h_l0_pt_vs_pz = new TH2D("h_l0_pt_vs_pz", "Transversal vs. longitudinal momentum for #Lambda^{0}; p_{z}/GeV/c; p_{t}/GeV/c", 40,-0.5,3.9,40,0,0.8);
	ntpLambda0->Project("h_l0_pt_vs_pz", "McTruth_pt: McTruth_pz", "McTruthMatch && Lambda0_HitTag");
	jenny::CreateDrawAndSaveHistogram(h_l0_pt_vs_pz, path+"/plots/", "lambda0_pt_vs_pz", save, close);

	gStyle->SetOptStat(0);

	TH1D * h_l0_m_nocut = new TH1D("h_l0_m_nocut", "Mass distribution for #Lambda^{0}; m/GeV/c^{2}; counts", 100,1,1.3);
	ntpLambda0->Project("h_l0_m_nocut", "Lambda0_m", "McTruthMatch && Lambda0_HitTag");

	TH1D * h_l0_m_vtxcut = new TH1D("h_l0_m_vtxcut", "Mass distribution for #Lambda^{0} with vtxcut; m/GeV/c^{2}; counts", 100,1,1.3);
	ntpLambda0->Project("h_l0_m_vtxcut", "Lambda0_m", "McTruthMatch && Lambda0_HitTag "+vtxcut);

	TH1D * h_l0_m_masscut = new TH1D("h_l0_m_masscut", "Mass distribution for #Lambda^{0} with vertex cut and mass cut; m/GeV/c^{2}; counts", 100,1,1.3);
	ntpLambda0->Project("h_l0_m_masscut", "Lambda0_m", "McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);

	jenny::CreateDrawAndSaveNHistograms(h_l0_m_nocut, h_l0_m_vtxcut, h_l0_m_masscut, "no cut", "VtxFit_prob>0.01", "VtxFit_prob>0.01 && MassFit_prob>0.01", path+"/plots/", "lambda0_m_diffcuts", save, close);

	gStyle->SetOptStat(1111);

	TH1D * h_l0_m_masscut2 = new TH1D("h_l0_m_masscut2", "Mass distribution for #Lambda^{0} with vertex cut and mass cut; m/GeV/c^{2}; counts", 100,1,1.3);
	ntpLambda0->Project("h_l0_m_masscut2", "Lambda0_m", "McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramFit(h_l0_m_masscut2, path+"/plots/", "lambda0_m_masscut2", save, close);

	TH1D * h_l0_vtxres_x = new TH1D("h_l0_vtxres_x", "resolution for x coordinate of vertex for #Lambda^{0}; x-x_{MC}; counts", 100,-1,1);
	ntpLambda0->Project("h_l0_vtxres_x", "VtxFit_diffvx", "McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_l0_vtxres_x, path+"/plots/", "lambda0_vtxres_x", save, close, false, 0.1, 1, true);

	TH1D * h_l0_vtxres_y = new TH1D("h_l0_vtxres_y", "resolution for y coordinate of vertex for #Lambda^{0}; y-y_{MC}; counts", 100,-1,1);
	ntpLambda0->Project("h_l0_vtxres_y", "VtxFit_diffvy", "McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_l0_vtxres_y, path+"/plots/", "lambda0_vtxres_y", save, close, false, 0.1, 1, true);

	TH1D * h_l0_vtxres_z = new TH1D("h_l0_vtxres_z", "resolution for z coordinate of vertex for #Lambda^{0}; z-z_{MC}; counts", 100,-1,1);
	ntpLambda0->Project("h_l0_vtxres_z", "VtxFit_diffvz", "McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_l0_vtxres_z, path+"/plots/", "lambda0_vtxres_z", save, close, false, 0.03, 1);//, true);

	TH1D * h_l0_costht = new TH1D("h_l0_costht", "cos(#Theta) distribution for #Lambda^{0}; cos(#Theta); counts", 100,-1,1);
	ntpLambda0->Project("h_l0_costht", "cos(VtxFit_tht)","McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_l0_costht, path+"/plots/", "lambda0_costht", save, close);

	TH1D * h_l0_tht = new TH1D("h_l0_tht", "#Theta distribution for #Lambda^{0}; #Theta/rad; counts", 100,-1,1);
	ntpLambda0->Project("h_l0_tht", "VtxFit_tht","McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_l0_tht, path+"/plots/", "lambda0_tht", save, close);






	//**** Get information about AntiLambda0 ****************************************************************************

	TH2D * h_al0_pt_vs_pz = new TH2D("h_al0_pt_vs_pz", "Transversal vs. longitudinal momentum for #bar{#Lambda}^{0}; p_{z}/GeV/c; p_{t}/GeV/c", 40,0.5,3.1,40,0,0.7);
	ntpAntiLambda0->Project("h_al0_pt_vs_pz", "MCTruth_pt: MCTruth_pz", "McTruthMatch && AntiLambda0_HitTag");
	jenny::CreateDrawAndSaveHistogram(h_al0_pt_vs_pz, path+"/plots/", "AntiLambda0_pt_vs_pz", save, close);

	gStyle->SetOptStat(0);
	TH1D * h_al0_m_nocut = new TH1D("h_al0_m_nocut", "Mass distribution for #bar{#Lambda}^{0}; m/GeV/c^{2}; counts", 100,1,1.3);
	ntpAntiLambda0->Project("h_al0_m_nocut", "AntiLambda0_m", "McTruthMatch && AntiLambda0_HitTag");

	TH1D * h_al0_m_vtxcut = new TH1D("h_al0_m_vtxcut", "Mass distribution for #bar{#Lambda}^{0} with vtxcut; m/GeV/c^{2}; counts", 100,1,1.3);
	ntpAntiLambda0->Project("h_al0_m_vtxcut", "AntiLambda0_m", "McTruthMatch && AntiLambda0_HitTag "+vtxcut);

	TH1D * h_al0_m_masscut = new TH1D("h_al0_m_masscut", "Mass distribution for #bar{#Lambda}^{0} with vertex cut and mass cut; m/GeV/c^{2}; counts", 100,1,1.3);
	ntpAntiLambda0->Project("h_al0_m_masscut", "AntiLambda0_m", "McTruthMatch && AntiLambda0_HitTag "+vtxcut+masscut);

	jenny::CreateDrawAndSaveNHistograms(h_al0_m_nocut, h_al0_m_vtxcut, h_al0_m_masscut, "no cut", "VtxFit_prob>0.01", "VtxFit_prob>0.01 && MassFit_prob>0.01", path+"/plots/", "AntiLambda0_m_diffcuts", save, close);

	gStyle->SetOptStat(1111);

	TH1D * h_al0_m_masscut2 = new TH1D("h_al0_m_masscut2", "Mass distribution for #bar{#Lambda}^{0} with vertex cut and mass cut; m/GeV/c^{2}; counts", 100,1,1.3);
	ntpAntiLambda0->Project("h_al0_m_masscut2", "AntiLambda0_m", "McTruthMatch && AntiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramFit(h_al0_m_masscut2, path+"/plots/", "AntiLambda0_m_masscut", save, close);

	TH1D * h_al0_vtxres_x = new TH1D("h_al0_vtxres_x", "resolution for x coordinate of vertex for #bar{#Lambda}^{0}; x-x_{MC}; counts", 100,-1,1);
	ntpAntiLambda0->Project("h_al0_vtxres_x", "VtxFit_diffvx", "McTruthMatch && AntiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_al0_vtxres_x, path+"/plots/", "AntiLambda0_vtxres_x", save, close, false, 0.1, 1, true);

	TH1D * h_al0_vtxres_y = new TH1D("h_al0_vtxres_y", "resolution for y coordinate of vertex for #bar{#Lambda}^{0}; y-y_{MC}; counts", 100,-1,1);
	ntpAntiLambda0->Project("h_al0_vtxres_y", "VtxFit_diffvy", "McTruthMatch && AntiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_al0_vtxres_y, path+"/plots/", "AntiLambda0_vtxres_y", save, close, false, 0.1, 1, true);

	TH1D * h_al0_vtxres_z = new TH1D("h_al0_vtxres_z", "resolution for z coordinate of vertex for #bar{#Lambda}^{0}; z-z_{MC}; counts", 100,-1,1);
	ntpAntiLambda0->Project("h_al0_vtxres_z", "VtxFit_diffvz", "McTruthMatch && AntiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_al0_vtxres_z, path+"/plots/", "AntiLambda0_vtxres_z", save, close, false, 0.05,1);//, true);


	TH1D * h_al0_costht = new TH1D("h_al0_costht", "cos(#Theta) distribution for #bar{#Lambda}^{0}; cos(#Theta); counts", 100,-1,1);
	ntpAntiLambda0->Project("h_al0_costht", "cos(VtxFit_tht)","McTruthMatch && AntiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_al0_costht, path+"/plots/", "AntiLambda0_costht", save, close);

	TH1D * h_al0_tht = new TH1D("h_al0_tht", "#Theta distribution for #bar{#Lambda}^{0}; #Theta/rad; counts", 100,-1,1);
	ntpAntiLambda0->Project("h_al0_tht", "VtxFit_tht","McTruthMatch && AntiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_al0_tht, path+"/plots/", "AntiLambda0_tht", save, close);
}
