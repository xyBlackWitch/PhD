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
#include "../../common_jenny.cpp"


void eval_lambda0(TString prefix="", bool save=kTRUE, bool close=kFALSE){

	//*** Input file
	if(prefix==""){
		TString inFile ="output_ana.root";
	}
	else{
		TString inFile = TString::Format("%s/output_ana.root", prefix.Data());
	}


	//*** get Data from Tree
	TFile * data = new TFile(inFile, "READ");

	TTree * ntpLambda0 = (TTree*) data->Get("ntpLambda0");
	TTree * ntpAntiLambda0 = (TTree*) data->Get("ntpAntiLambda0");

	//*******************************************************************************************************************

	gStyle->SetOptStat(1111);
	gStyle->SetStatX(0.85);
	gStyle->SetStatY(0.85);

	TString vtxcut = "&& VtxFit_HowGood==1 ";
	TString masscut = "&& MassFit_prob>0.01 ";

	//**** Get information about Lambda0 ********************************************************************************

	//*** pt vs pz
	TH2D * h_l0_pt_vs_pz = new TH2D("h_l0_pt_vs_pz", "Transverse vs. longitudinal momentum for #Lambda; p_{z} [GeV/c]; p_{t} [GeV/c]", 200,-0.5,3.9,200,0,0.8);
	ntpLambda0->Project("h_l0_pt_vs_pz", "McTruth_pt: McTruth_pz", "McTruthMatch && HitTag");
	jenny::CreateDrawAndSaveHistogram(h_l0_pt_vs_pz, prefix+"/plots/", "lambda0_pt_vs_pz", save, close);

	TH2D * h_l0_pt_vs_pz_cut = new TH2D("h_l0_pt_vs_pz_cut", "Transverse vs. longitudinal momentum for #Lambda with cut; p_{z} [GeV/c]; p_{t} [GeV/c]", 200,0,4,200,0,0.8);
	ntpLambda0->Project("h_l0_pt_vs_pz_cut", "McTruth_pt: McTruth_pz", "McTruthMatch && HitTag" + vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_l0_pt_vs_pz_cut, prefix+"/plots/", "lambda0_pt_vs_pz_cut", save, close);

//	gStyle->SetOptStat(0);

	//*** mass distribution
	TH1D * h_l0_m_nocut = new TH1D("h_l0_m_nocut", "Mass distribution for #Lambda; M [GeV/c^{2}]; counts", 500,1.07,1.16);
	ntpLambda0->Project("h_l0_m_nocut", "VtxFit_m", "McTruthMatch && HitTag");

	TH1D * h_l0_m_vtxcut = new TH1D("h_l0_m_vtxcut", "Mass distribution for #Lambda with vtxcut; M [GeV/c^{2}]; counts", 500,1.07,1.16);
	ntpLambda0->Project("h_l0_m_vtxcut", "VtxFit_m", "McTruthMatch && HitTag "+vtxcut);

	TH1D * h_l0_m_masscut = new TH1D("h_l0_m_masscut", "Mass distribution for #Lambda with vertex cut and mass cut; M [GeV/c^{2}]; counts", 500,1.07,1.16);
	ntpLambda0->Project("h_l0_m_masscut", "VtxFit_m", "McTruthMatch && HitTag "+vtxcut+masscut);

	jenny::CreateDrawAndSaveNHistograms(h_l0_m_nocut, h_l0_m_vtxcut, h_l0_m_masscut, "mass window", "VtxFit_prob>0.01", "VtxFit_prob>0.01 && MassFit_prob>0.01", prefix+"/plots/", "lambda0_m_diffcuts", save, close);
//	gStyle->SetOptStat(1111);

	TH1D * h_l0_m_masscut2\ = new TH1D("h_l0_m_masscut2", "Mass distribution for #Lambda with vertex cut and mass cut; M [GeV/c^{2}]; counts", 500,1.07,1.16);
	ntpLambda0->Project("h_l0_m_masscut2", "VtxFit_m", "McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_l0_m_masscut2, prefix+"/plots/", "lambda0_m_masscut2", save, close);//, false, 0.004,0.5);

	TH1D * h_l0_m_mc = new TH1D("h_l0_m_mc", "Mass distribution for #Lambda to MC truth; (M_{#Lambda}-M_{PDG}) [MeV/c^{2}]; counts", 500,-10,10);
	ntpLambda0->Project("h_l0_m_mc", "(VtxFit_m-McTruth_m)*1000", "McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_l0_m_mc, prefix+"/plots/", "lambda0_m_diff_to_mc", save, close);

	//*** vertex resolution
	TH1D * h_l0_vtxres_x = new TH1D("h_l0_vtxres_x", "resolution for x coordinate of vertex for #Lambda; x-x_{MC} [mm]; counts", 500,-11,11);
	ntpLambda0->Project("h_l0_vtxres_x", "VtxFit_diffvx*10", "McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramFWHM(h_l0_vtxres_x, prefix+"/plots/", "lambda0_vtxres_x", save, close);

	TH1D * h_l0_vtxres_y = new TH1D("h_l0_vtxres_y", "resolution for y coordinate of vertex for #Lambda; y-y_{MC} [mm]; counts", 500,-11,11);
	ntpLambda0->Project("h_l0_vtxres_y", "VtxFit_diffvy*10", "McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramFWHM(h_l0_vtxres_y, prefix+"/plots/", "lambda0_vtxres_y", save, close);

	TH1D * h_l0_vtxres_z = new TH1D("h_l0_vtxres_z", "resolution for z coordinate of vertex for #Lambda; z-z_{MC} [mm]; counts", 500,-11,11);
	ntpLambda0->Project("h_l0_vtxres_z", "VtxFit_diffvz*10", "McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramFWHM(h_l0_vtxres_z, prefix+"/plots/", "lambda0_vtxres_z", save, close);

	//*** angular distribution
	TH1D * h_l0_costht = new TH1D("h_l0_costht", "cos(#Theta) distribution for #Lambda; cos(#Theta); counts", 500,-1,1.05);
	ntpLambda0->Project("h_l0_costht", "cos(VtxFit_tht)","McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_l0_costht, prefix+"/plots/", "lambda0_costht", save, close);

	TH1D * h_l0_tht = new TH1D("h_l0_tht", "#Theta distribution for #Lambda; #Theta [rad]; counts", 500,-1,1.05);
	ntpLambda0->Project("h_l0_tht", "VtxFit_tht","McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_l0_tht, prefix+"/plots/", "lambda0_tht", save, close);

	//*** fit values

	TH1D * h_l0_chi2 = new TH1D("h_l0_chi2", "#chi^{2} distribution for #Lambda; #chi^{2}; counts", 500,0,100);
	ntpLambda0->Project("h_l0_chi2", "VtxFit_chisq");
	jenny::CreateDrawAndSaveHistogram(h_l0_chi2, prefix+"/plots/", "h_l0_chi2", save, close,true);


	TH1D * h_l0_prob = new TH1D("h_l0_prob", "propability distribution for #Lambda; prob; counts", 500,0,1);
	ntpLambda0->Project("h_l0_prob", "VtxFit_prob");
	jenny::CreateDrawAndSaveHistogram(h_l0_prob, prefix+"/plots/", "h_l0_prob", save, close, true);



	//**** Get information about AntiLambda0 ****************************************************************************

	//*** pt vs pz
	TH2D * h_al0_pt_vs_pz = new TH2D("h_al0_pt_vs_pz", "Transverse vs. longitudinal momentum for #bar{#Lambda}; p_{z} [GeV/c]; p_{t} [GeV/c]", 200,0.5,3.1,200,0,0.7);
	ntpAntiLambda0->Project("h_al0_pt_vs_pz", "McTruth_pt: McTruth_pz", "McTruthMatch && HitTag");
	jenny::CreateDrawAndSaveHistogram(h_al0_pt_vs_pz, prefix+"/plots/", "antiLambda0_pt_vs_pz", save, close);

	TH2D * h_al0_pt_vs_pz_cut = new TH2D("h_al0_pt_vs_pz_cut", "Transverse vs. longitudinal momentum for #bar{#Lambda} with cut; p_{z} [GeV/c]; p_{t} [GeV/c]", 200,0,4,200,0,0.7);
	ntpAntiLambda0->Project("h_al0_pt_vs_pz_cut", "McTruth_pt: McTruth_pz", "McTruthMatch && HitTag"+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_al0_pt_vs_pz_cut, prefix+"/plots/", "antiLambda0_pt_vs_pz_cut", save, close);

//	gStyle->SetOptStat(0);

	//*** mass distribution

	TH1D * h_al0_m_nocut = new TH1D("h_al0_m_nocut", "Mass distribution for #bar{#Lambda}; M [GeV/c^{2}]; counts", 500,1.07,1.16);
	ntpAntiLambda0->Project("h_al0_m_nocut", "VtxFit_m", "McTruthMatch && HitTag");

	TH1D * h_al0_m_vtxcut = new TH1D("h_al0_m_vtxcut", "Mass distribution for #bar{#Lambda} with vtxcut; M [GeV/c^{2}]; counts", 500,1.07,1.16);
	ntpAntiLambda0->Project("h_al0_m_vtxcut", "VtxFit_m", "McTruthMatch && HitTag "+vtxcut);

	TH1D * h_al0_m_masscut = new TH1D("h_al0_m_masscut", "Mass distribution for #bar{#Lambda} with vertex cut and mass cut; M [GeV/c^{2}]; counts", 500,1.07,1.16);
	ntpAntiLambda0->Project("h_al0_m_masscut", "VtxFit_m", "McTruthMatch && HitTag "+vtxcut+masscut);

	jenny::CreateDrawAndSaveNHistograms(h_al0_m_nocut, h_al0_m_vtxcut, h_al0_m_masscut, "mass window", "VtxFit_prob>0.01", "VtxFit_prob>0.01 && MassFit_prob>0.01", prefix+"/plots/", "antiLambda0_m_diffcuts", save, close);

	//gStyle->SetOptStat(1111);

	TH1D * h_al0_m_masscut2 = new TH1D("h_al0_m_masscut2", "Mass distribution for #bar{#Lambda} with vertex cut and mass cut; M [GeV/c^{2}]; counts", 500,1.07,1.16);
	ntpAntiLambda0->Project("h_al0_m_masscut2", "VtxFit_m", "McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_al0_m_masscut2, prefix+"/plots/", "antiLambda0_m_masscut", save, close);//, false, 0.004,0.5);

	TH1D * h_al0_m_mc = new TH1D("h_al0_m_mc", "Mass distribution for #bar{#Lambda} to MC truth; (M_{#Lambda}-M_{PDG}) [MeV/c^{2}]; counts", 500,-10,10);
	ntpAntiLambda0->Project("h_al0_m_mc", "(VtxFit_m-McTruth_m)*1000", "McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramDoulbeFit(h_al0_m_mc, prefix+"/plots/", "antiLambda0_m_diff_to_mc", save, close);


	//*** vertex resolution

	TH1D * h_al0_vtxres_x = new TH1D("h_al0_vtxres_x", "resolution for x coordinate of vertex for #bar{#Lambda}; x-x_{MC} [mm]; counts", 500,-11,11);
	ntpAntiLambda0->Project("h_al0_vtxres_x", "VtxFit_diffvx*10", "McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramFWHM(h_al0_vtxres_x, prefix+"/plots/", "antiLambda0_vtxres_x", save, close);

	TH1D * h_al0_vtxres_y = new TH1D("h_al0_vtxres_y", "resolution for y coordinate of vertex for #bar{#Lambda}; y-y_{MC} [mm]; counts", 500,-11,11);
	ntpAntiLambda0->Project("h_al0_vtxres_y", "VtxFit_diffvy*10", "McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramFWHM(h_al0_vtxres_y, prefix+"/plots/", "antiLambda0_vtxres_y", save, close);

	TH1D * h_al0_vtxres_z = new TH1D("h_al0_vtxres_z", "resolution for z coordinate of vertex for #bar{#Lambda}; z-z_{MC} [mm]; counts", 500,-11,11);
	ntpAntiLambda0->Project("h_al0_vtxres_z", "VtxFit_diffvz*10", "McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogramFWHM(h_al0_vtxres_z, prefix+"/plots/", "antiLambda0_vtxres_z", save, close);

	//*** angular distribution
	TH1D * h_al0_costht = new TH1D("h_al0_costht", "cos(#Theta) distribution for #bar{#Lambda}; cos(#Theta); counts", 500,-1,1.05);
	ntpAntiLambda0->Project("h_al0_costht", "cos(VtxFit_tht)","McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_al0_costht, prefix+"/plots/", "antiLambda0_costht", save, close);

	TH1D * h_al0_tht = new TH1D("h_al0_tht", "#Theta distribution for #bar{#Lambda}; #Theta [rad]; counts", 500,-1,1.05);
	ntpAntiLambda0->Project("h_al0_tht", "VtxFit_tht","McTruthMatch && HitTag "+vtxcut+masscut);
	jenny::CreateDrawAndSaveHistogram(h_al0_tht, prefix+"/plots/", "antiLambda0_tht", save, close);

	//*** fit values

	TH1D * h_al0_chi2 = new TH1D("h_al0_chi2", "#chi^{2} distribution for #bar{#Lambda}; #chi^{2}; counts", 500,0,100);
	ntpAntiLambda0->Project("h_al0_chi2", "VtxFit_chisq");
	jenny::CreateDrawAndSaveHistogram(h_al0_chi2, prefix+"/plots/", "h_al0_chi2", save, close,true);

	TH1D * h_al0_prob = new TH1D("h_al0_prob", "propability distribution for #bar{#Lambda}; prob; counts", 500,0,1);
	ntpAntiLambda0->Project("h_al0_prob", "VtxFit_prob");
	jenny::CreateDrawAndSaveHistogram(h_al0_prob, prefix+"/plots/", "h_al0_prob", save, close, true);
}
