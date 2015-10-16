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

	//*** pt vs pz
	TH2D * h_l0_pt_vs_pz = new TH2D("h_l0_pt_vs_pz", "Transversal vs. longitudinal momentum for #Lambda^{0}; p_{z}/GeV/c; p_{t}/GeV/c", 200,-0.5,3.9,200,0,0.8);
	ntpLambda0->Project("h_l0_pt_vs_pz", "McTruth_pt: McTruth_pz", "McTruthMatch && Lambda0_HitTag");
	jenny::CreateDrawAndSaveHistogram(h_l0_pt_vs_pz, path+"/plots/", "lambda0_pt_vs_pz", save, close);

	TH2D * h_l0_pt_vs_pz_cut = new TH2D("h_l0_pt_vs_pz_cut", "Transversal vs. longitudinal momentum for #Lambda^{0} with cut; p_{z}/GeV/c; p_{t}/GeV/c", 200,-0.5,3.9,200,0,0.8);
	ntpLambda0->Project("h_l0_pt_vs_pz_cut", "McTruth_pt: McTruth_pz", "McTruthMatch && Lambda0_HitTag" + vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_l0_pt_vs_pz_cut, path+"/plots/", "lambda0_pt_vs_pz_cut", save, close);

	gStyle->SetOptStat(0);

	//*** mass distribution
	TH1D * h_l0_m_nocut = new TH1D("h_l0_m_nocut", "Mass distribution for #Lambda^{0}; m/GeV/c^{2}; counts", 500,1,1.2);
	ntpLambda0->Project("h_l0_m_nocut", "Lambda0_m", "McTruthMatch && Lambda0_HitTag");

	TH1D * h_l0_m_vtxcut = new TH1D("h_l0_m_vtxcut", "Mass distribution for #Lambda^{0} with vtxcut; m/GeV/c^{2}; counts", 500,1,1.2);
	ntpLambda0->Project("h_l0_m_vtxcut", "Lambda0_m", "McTruthMatch && Lambda0_HitTag "+vtxcut);

	TH1D * h_l0_m_masscut = new TH1D("h_l0_m_masscut", "Mass distribution for #Lambda^{0} with vertex cut and mass cut; m/GeV/c^{2}; counts", 500,1,1.2);
	ntpLambda0->Project("h_l0_m_masscut", "Lambda0_m", "McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);

	jenny::CreateDrawAndSaveNHistograms(h_l0_m_nocut, h_l0_m_vtxcut, h_l0_m_masscut, "no cut", "VtxFit_prob>0.01", "VtxFit_prob>0.01 && MassFit_prob>0.01", path+"/plots/", "lambda0_m_diffcuts", save, close);

	gStyle->SetOptStat(1111);

	TH1D * h_l0_m_masscut2 = new TH1D("h_l0_m_masscut2", "Mass distribution for #Lambda^{0} with vertex cut and mass cut; m/GeV/c^{2}; counts", 500,1,1.2);
	ntpLambda0->Project("h_l0_m_masscut2", "Lambda0_m", "McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_l0_m_masscut2, path+"/plots/", "lambda0_m_masscut2", save, close, false, 0.01,0.1, true);

	//*** vertex resolution
	TH1D * h_l0_vtxres_x = new TH1D("h_l0_vtxres_x", "resolution for x coordinate of vertex for #Lambda^{0}; x-x_{MC}; counts", 500,-1,1.05);
	ntpLambda0->Project("h_l0_vtxres_x", "VtxFit_diffvx", "McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_l0_vtxres_x, path+"/plots/", "lambda0_vtxres_x", save, close, false, 0.1, 1, true);

	TH1D * h_l0_vtxres_y = new TH1D("h_l0_vtxres_y", "resolution for y coordinate of vertex for #Lambda^{0}; y-y_{MC}; counts", 500,-1,1.05);
	ntpLambda0->Project("h_l0_vtxres_y", "VtxFit_diffvy", "McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_l0_vtxres_y, path+"/plots/", "lambda0_vtxres_y", save, close, false, 0.1, 1, true);

	TH1D * h_l0_vtxres_z = new TH1D("h_l0_vtxres_z", "resolution for z coordinate of vertex for #Lambda^{0}; z-z_{MC}; counts", 500,-1,1.05);
	ntpLambda0->Project("h_l0_vtxres_z", "VtxFit_diffvz", "McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_l0_vtxres_z, path+"/plots/", "lambda0_vtxres_z", save, close, false, 0.03, 1, true);

	//*** angular distribution
	TH1D * h_l0_costht = new TH1D("h_l0_costht", "cos(#Theta) distribution for #Lambda^{0}; cos(#Theta); counts", 500,-1,1.05);
	ntpLambda0->Project("h_l0_costht", "cos(VtxFit_tht)","McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_l0_costht, path+"/plots/", "lambda0_costht", save, close);

	TH1D * h_l0_tht = new TH1D("h_l0_tht", "#Theta distribution for #Lambda^{0}; #Theta/rad; counts", 500,-1,1.05);
	ntpLambda0->Project("h_l0_tht", "VtxFit_tht","McTruthMatch && Lambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_l0_tht, path+"/plots/", "lambda0_tht", save, close);

	//*** fit values

	TH1D * h_l0_chi2 = new TH1D("h_l0_chi2", "#chi^{2} distribution for #Lambda^{0}; #chi^2; counts", 500,0,100);
	ntpLambda0->Project("h_l0_chi2", "VtxFit_chisq");
	jenny::CreateDrawAndSaveHistogram(h_l0_chi2, path+"/plots/", "h_l0_chi2", save, close);


	TH1D * h_l0_prob = new TH1D("h_l0_prob", "propability distribution for #Lambda^{0}; prob; counts", 500,0,1);
	ntpLambda0->Project("h_l0_prob", "VtxFit_prob");
	jenny::CreateDrawAndSaveHistogram(h_l0_prob, path+"/plots/", "h_l0_prob", save, close);



	//**** Get information about AntiLambda0 ****************************************************************************

	//*** pt vs pz
	TH2D * h_al0_pt_vs_pz = new TH2D("h_al0_pt_vs_pz", "Transversal vs. longitudinal momentum for #bar{#Lambda}^{0}; p_{z}/GeV/c; p_{t}/GeV/c", 200,0.5,3.1,200,0,0.7);
	ntpAntiLambda0->Project("h_al0_pt_vs_pz", "McTruth_pt: McTruth_pz", "McTruthMatch && antiLambda0_HitTag");
	jenny::CreateDrawAndSaveHistogram(h_al0_pt_vs_pz, path+"/plots/", "antiLambda0_pt_vs_pz", save, close);

	TH2D * h_al0_pt_vs_pz_cut = new TH2D("h_al0_pt_vs_pz_cut", "Transversal vs. longitudinal momentum for #bar{#Lambda}^{0} with cut; p_{z}/GeV/c; p_{t}/GeV/c", 200,0.5,3.1,200,0,0.7);
	ntpAntiLambda0->Project("h_al0_pt_vs_pz_cut", "McTruth_pt: McTruth_pz", "McTruthMatch && antiLambda0_HitTag"+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_al0_pt_vs_pz_cut, path+"/plots/", "antiLambda0_pt_vs_pz_cut", save, close);

	gStyle->SetOptStat(0);

	//*** mass distribution

	TH1D * h_al0_m_nocut = new TH1D("h_al0_m_nocut", "Mass distribution for #bar{#Lambda}^{0}; m/GeV/c^{2}; counts", 500,1,1.2);
	ntpAntiLambda0->Project("h_al0_m_nocut", "antiLambda0_m", "McTruthMatch && antiLambda0_HitTag");

	TH1D * h_al0_m_vtxcut = new TH1D("h_al0_m_vtxcut", "Mass distribution for #bar{#Lambda}^{0} with vtxcut; m/GeV/c^{2}; counts", 500,1,1.2);
	ntpAntiLambda0->Project("h_al0_m_vtxcut", "antiLambda0_m", "McTruthMatch && antiLambda0_HitTag "+vtxcut);

	TH1D * h_al0_m_masscut = new TH1D("h_al0_m_masscut", "Mass distribution for #bar{#Lambda}^{0} with vertex cut and mass cut; m/GeV/c^{2}; counts", 500,1,1.2);
	ntpAntiLambda0->Project("h_al0_m_masscut", "antiLambda0_m", "McTruthMatch && antiLambda0_HitTag "+vtxcut+masscut);

	jenny::CreateDrawAndSaveNHistograms(h_al0_m_nocut, h_al0_m_vtxcut, h_al0_m_masscut, "no cut", "VtxFit_prob>0.01", "VtxFit_prob>0.01 && MassFit_prob>0.01", path+"/plots/", "antiLambda0_m_diffcuts", save, close);

	gStyle->SetOptStat(1111);

	TH1D * h_al0_m_masscut2 = new TH1D("h_al0_m_masscut2", "Mass distribution for #bar{#Lambda}^{0} with vertex cut and mass cut; m/GeV/c^{2}; counts", 500,1,1.2);
	ntpAntiLambda0->Project("h_al0_m_masscut2", "antiLambda0_m", "McTruthMatch && antiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_al0_m_masscut2, path+"/plots/", "antiLambda0_m_masscut", save, close, false, 0.01, 0.1, true);


	//*** vertex resolution

	TH1D * h_al0_vtxres_x = new TH1D("h_al0_vtxres_x", "resolution for x coordinate of vertex for #bar{#Lambda}^{0}; x-x_{MC}; counts", 500,-1,1.05);
	ntpAntiLambda0->Project("h_al0_vtxres_x", "VtxFit_diffvx", "McTruthMatch && antiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_al0_vtxres_x, path+"/plots/", "antiLambda0_vtxres_x", save, close, false, 0.1, 1, true);

	TH1D * h_al0_vtxres_y = new TH1D("h_al0_vtxres_y", "resolution for y coordinate of vertex for #bar{#Lambda}^{0}; y-y_{MC}; counts", 500,-1,1.05);
	ntpAntiLambda0->Project("h_al0_vtxres_y", "VtxFit_diffvy", "McTruthMatch && antiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_al0_vtxres_y, path+"/plots/", "antiLambda0_vtxres_y", save, close, false, 0.08, 1, true);

	TH1D * h_al0_vtxres_z = new TH1D("h_al0_vtxres_z", "resolution for z coordinate of vertex for #bar{#Lambda}^{0}; z-z_{MC}; counts", 500,-1,1.05);
	ntpAntiLambda0->Project("h_al0_vtxres_z", "VtxFit_diffvz", "McTruthMatch && antiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_al0_vtxres_z, path+"/plots/", "antiLambda0_vtxres_z", save, close, false, 0.1,1, true);

	//*** angular distribution
	TH1D * h_al0_costht = new TH1D("h_al0_costht", "cos(#Theta) distribution for #bar{#Lambda}^{0}; cos(#Theta); counts", 500,-1,1.05);
	ntpAntiLambda0->Project("h_al0_costht", "cos(VtxFit_tht)","McTruthMatch && antiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_al0_costht, path+"/plots/", "antiLambda0_costht", save, close);

	TH1D * h_al0_tht = new TH1D("h_al0_tht", "#Theta distribution for #bar{#Lambda}^{0}; #Theta/rad; counts", 500,-1,1.05);
	ntpAntiLambda0->Project("h_al0_tht", "VtxFit_tht","McTruthMatch && antiLambda0_HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_al0_tht, path+"/plots/", "antiLambda0_tht", save, close);

	//*** fit values

	TH1D * h_al0_chi2 = new TH1D("h_al0_chi2", "#chi^{2} distribution for #bar{#Lambda}^{0}; #chi^2; counts", 500,0,100);
	ntpAntiLambda0->Project("h_al0_chi2", "VtxFit_chisq");
	jenny::CreateDrawAndSaveHistogram(h_al0_chi2, path+"/plots/", "h_al0_chi2", save, close);

	TH1D * h_al0_prob = new TH1D("h_al0_prob", "propability distribution for #bar{#Lambda}^{0}; prob; counts", 500,0,1);
	ntpAntiLambda0->Project("h_al0_prob", "VtxFit_prob");
	jenny::CreateDrawAndSaveHistogram(h_al0_prob, path+"/plots/", "h_al0_prob", save, close);
}
