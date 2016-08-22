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
#include "../../common_jenny.cpp"


void eval_mass_Xi1690(TString prefix="", bool save=kTRUE, bool close=kFALSE, TString path=""){

	if (path != "") prefix=path+prefix;

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


	//*******************************************************************************************************************

	TString vtxcut = "&& HitTag==1 && VtxFit_HowGood==1 ";


	//**** Get information about XiMinus ********************************************************************************

		TH1D * h_xi_m_nocut = new TH1D("h_xi_m_nocut", "Mass distribution for #Xi(1820)^{-}; m/GeV/c^{2}; counts", 500,1.6,2);
		ntpXiMinus1690->Project("h_xi_m_nocut", "VtxFit_m", "McTruthMatch");

		TH1D * h_xi_m_vtxcut = new TH1D("h_xi_m_vtxcut", "Mass distribution for #Xi(1820)^{-} with vtxcut; m/GeV/c^{2}; counts", 500,1.6,2);
		ntpXiMinus1690->Project("h_xi_m_vtxcut", "VtxFit_m", "McTruthMatch "+vtxcut);


		TF1 * fitfunc = jenny::GetVoigtFitMC(h_xi_m_vtxcut, 1.67, 1.93, 1.690, 0.03);



	  	TCanvas * c = new TCanvas("c", "c", 0,0,1500,1000);

		h_xi_m_nocut->SetLineColor(kBlue);
		h_xi_m_vtxcut->SetLineColor(kBlack);



	  	TLegend * legend = new TLegend(0.7,0.62,0.86,0.795, "");
	  	legend->AddEntry(h_xi_m_nocut, "mass window", "l");
	  	legend->AddEntry(h_xi_m_vtxcut, "VtxFit_prob>0.01 (Best Cand.)", "l");


	  	h_xi_m_nocut->Draw();
	  	h_xi_m_vtxcut->Draw("SAME");
	  	fitfunc->Draw("SAME");


	  	legend->Draw();





}
