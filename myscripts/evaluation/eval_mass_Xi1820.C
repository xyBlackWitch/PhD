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


void eval_mass_Xi1820(TString path="", bool save=kTRUE, bool close=kFALSE){

	//*** Input file
	TString inFile = TString::Format("%s/output_ana.root", path.Data());


	//*** get Data from Tree
	TFile * data = new TFile(inFile, "READ");

	TTree * ntpXiMinus1820 = (TTree*) data->Get("ntpXiMinus1820");
	TTree * ntpMC = (TTree*) data->Get("ntpMC");

	//*******************************************************************************************************************

	gStyle->SetOptStat(1111);
	gStyle->SetStatX(0.85);
	gStyle->SetStatY(0.85);

	TString vtxcut = "&& VtxFit_HowGood==1 ";


	//**** Get information about XiMinus ********************************************************************************
//		gStyle->SetOptStat(0);

//		TH1D * h_xi_m_nocut = new TH1D("h_xi_m_nocut", "Mass distribution for #Xi(1820)^{-}; m/GeV/c^{2}; counts", 500,1.7,2);
//		ntpXiMinus1820->Project("h_xi_m_nocut", "VtxFit_m", "McTruthMatch");

		TH1D * h_xi_m_vtxcut = new TH1D("h_xi_m_vtxcut", "Mass distribution for #Xi(1820)^{-} with vtxcut; m/GeV/c^{2}; counts", 500,1.6,2);
		ntpXiMinus1820->Project("h_xi_m_vtxcut", "VtxFit_m", "McTruthMatch "+vtxcut);

		TF1 * fitfunc = jenny::tripleGaussFit(h_xi_m_vtxcut, false, 0.02, 0.1);

//		TH1D * h_xi_m_mc = new TH1D("h_xi_m_mc", "MC mass distribution for  #Xi(1820)^{-} ;m/GeV/c^{2}; counts(MC)", 500,1.7,2);
//		ntpMC->Project("h_xi_m_vtxcut", "m", "pdg==23314 && moth==0");

	  	TCanvas * c = new TCanvas("c", "c", 0,0,1500,1000);

//		h_xi_m_nocut->SetLineColor(kBlue);
//		h_xi_m_vtxcut->SetLineColor(kRed);
//		h_xi_m_mc->SetLineColor(kBlack);
//
//
//	  	TLegend * legend = new TLegend(0.7,0.62,0.86,0.795, "");
//	  	legend->AddEntry(h_xi_m_nocut, "mass window", "l");
//	  	legend->AddEntry(h_xi_m_vtxcut, "VtxFit_prob>0.01", "l");
//	  	legend->AddEntry(h_xi_m_mc, "MC", "l");
//
//	  	h_xi_m_nocut->Draw();
	  	h_xi_m_vtxcut->Draw();
	  	fitfunc->Draw("SAME");
//	  	h_xi_m_mc->Draw("SAME Y+");
//
//	  	legend->Draw();





}
