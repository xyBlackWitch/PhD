/**
 * @file eval_final_states.C
 * @mainpage eval_final_states.C Get Data from analysis file and creates and saves different histograms
 *
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Get Data from analysis file and create and save different histograms
 * @details This files get the data for Xi particle from the analysis of
 * pbar p -> Xi+ Xi(1820)-
 * 			 |	  |-> XiMinus + K-
 * 			 |			|-> p + Pi-
 * 			 |
 * 			 |->  XiPlus + Pi+
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


void eval_MC(TString path="", bool save=kTRUE, bool close=kFALSE){


	//*** Input file
	if(prefix==""){
		TString inFile ="output_ana.root";
	}
	else{
		TString inFile = TString::Format("%s/output_ana.root", prefix.Data());
	}


	//*** get Data from Tree
	TFile * data = new TFile(inFile, "READ");

	TTree * ntpMC = (TTree*) data->Get("ntpMC");

	//**********************AntiXiPlus******************************
//	cout << "Creating histograms for XiPlus" << endl;
//	TH2D * h_aXi_pt_vs_pz = new TH2D("h_aXi_pt_vs_pz", "transversal vs longitudinal momentum for #bar{#Xi}^{+}; pz/GeV/c; pt/GeV/c", 200,0,4,200,0,0.7);
//	ntpMC->Project("h_aXi_pt_vs_pz", "sqrt(px**2+py**2):pz", "pdg==-3312 && moth==0");
//	jenny::CreateDrawAndSaveHistogram(h_aXi_pt_vs_pz, path+"/plots/", "XiPlus_MC_pt_vs_pz", save, close);
//
//	TH1D * h_aXi_tht = new TH1D("h_aXi_tht", "angle distribution for #bar{#Xi}; #theta/rad; coutns", 500,0,3);
//	ntpMC->Project("h_aXi_tht", "tht", "pdg==-3312 && moth==0");
//	jenny::CreateDrawAndSaveHistogram(h_aXi_tht, path+"/plots/", "XiPlus_MC_tht", save, close);

//	//**********************XiMinus1820*****************************
//	cout << "Creating histograms for XiMinus1820" << endl;
//
//	TH2D * h_Xi1820_pt_vs_pz = new TH2D("h_Xi1820_pt_vs_pz", "transverse vs longitudinal momentum for #Xi(1820); pz/GeV/c; pt/GeV/c", 200,0,4,200,0,0.7);
//	ntpMC->Project("h_Xi1820_pt_vs_pz", "sqrt(px**2+py**2):pz", "pdg==23314 && moth==0");
//	jenny::CreateDrawAndSaveHistogram(h_Xi1820_pt_vs_pz, path+"/plots/", "XiMinus1820_MC_pt_vs_pz_new", save, close);

//	TH1D * h_Xi1820_tht = new TH1D("h_Xi1820_tht", "#theta distribution for #Xi(1820); #theta/rad; counts", 500,0,1);
//	ntpMC->Project("h_Xi1820_tht", "tht", "pdg==23314 && moth==0");
//	jenny::CreateDrawAndSaveHistogram(h_Xi1820_tht, path+"/plots/", "XiMinus1820_MC_tht", save, close);


//	//**********************Lambda0*********************************
//	cout << "Creating histograms for Lambda0" << endl;

//	TH2D * h_l0_pt_vs_pz = new TH2D("h_l0_pt_vs_pz", "transversal vs longitudinal momentum for #Lambda^{0}; pz/GeV/c; pt/GeV/c", 200,0,4,200,0,0.8);
//	ntpMC->Project("h_l0_pt_vs_pz", "sqrt(px**2+py**2):pz", "pdg==3122 && moth==1");
//	jenny::CreateDrawAndSaveHistogram(h_l0_pt_vs_pz, path+"/plots/", "Lambda0_MC_pt_vs_pz", save, close);
//
//	TH1D * h_l0_tht = new TH1D("h_l0_tht", "transversal vs longitudinal momentum for #Lambda^{0}; pz/GeV/c; pt/GeV/c", 500,0,3);
//	ntpMC->Project("h_l0_tht", "tht", "pdg==3122 && moth==1");
//	jenny::CreateDrawAndSaveHistogram(h_l0_tht, path+"/plots/", "Lambda0_MC_tht", save, close);
//
//
//	//**********************AntiLambda0*****************************
	cout << "Creating histograms for AntiLambda0" << endl;
	TH2D * h_al0_pt_vs_pz = new TH2D("h_al0_pt_vs_pz", "transversal vs longitudinal momentum for #bar{#Lambda}^{0}; pz/GeV/c; pt/GeV/c", 200,-1,4,200,0,0.8);
	ntpMC->Project("h_al0_pt_vs_pz", "sqrt(px**2+py**2):pz", "pdg==-3122 && moth==2");
	jenny::CreateDrawAndSaveHistogram(h_al0_pt_vs_pz, path+"/plots/", "AntiLambda0_MC_pt_vs_pz", save, close);
//
//	TH1D * h_al0_tht = new TH1D("h_al0_tht", "transversal vs longitudinal momentum for #bar{#Lambda}^{0}; pz/GeV/c; pt/GeV/c", 500,0,3);
//	ntpMC->Project("h_al0_tht", "tht", "pdg==-3122 && moth==2");
//	jenny::CreateDrawAndSaveHistogram(h_al0_tht, path+"/plots/", "AntiLambda0_MC_pt_vs_pz", save, close);





}
