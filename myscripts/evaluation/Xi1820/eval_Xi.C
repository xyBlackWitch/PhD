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
#include "../../common_jenny.cpp"


void eval_Xi(TString prefix="", bool save=kTRUE, bool close=kFALSE){

//	TString prefix;
//	if (path != "") prefix = path+"/"+pre;
//	else prefix = pre;

	//*** Input file
	TString inFile = TString::Format("%s/output_ana.root", prefix.Data());


	//*** get Data from Tree
	TFile * data = new TFile(inFile, "READ");

	TTree * ntpXiMinus = (TTree*) data->Get("ntpXiMinus");
	TTree * ntpXiPlus = (TTree*) data->Get("ntpXiPlus");

	//*******************************************************************************************************************

	//gStyle->SetOptStat(1111);
	//gStyle->SetStatX(0.85);
	//gStyle->SetStatY(0.85);

	TString vtxcut = "& HitTag==1 & VtxFit_HowGood==1 ";
	TString masscut = "&& MassFit_HowGood>0 ";

	//**** Get information about XiMinus ********************************************************************************
	if(ntpXiMinus!=0X0){
		TH2D * h_xi_pt_vs_pz = new TH2D("h_xi_pt_vs_pz", "Transverse vs. longitudinal momentum for #Xi^{-}; p_{z} [GeV/c]; p_{t} [GeV/c]", 200,0.5,3.5,200,0,0.8);
		ntpXiMinus->Project("h_xi_pt_vs_pz", "MCTruth_pt: MCTruth_pz", "McTruthMatch");
		jenny::CreateDrawAndSaveHistogram(h_xi_pt_vs_pz, prefix+"/plots/", "XiMinus_pt_vs_pz", save, close);

		TH2D * h_xi_pt_vs_pz_cut = new TH2D("h_xi_pt_vs_pz_cut", "Transverse vs. longitudinal momentum for #Xi^{-} after cut; p_{z} [GeV/c]; p_{t} [GeV/c]", 200,0.5,3.5,200,0,0.8);
		ntpXiMinus->Project("h_xi_pt_vs_pz_cut", "MCTruth_pt: MCTruth_pz", "McTruthMatch" +vtxcut + masscut);
		jenny::CreateDrawAndSaveHistogram(h_xi_pt_vs_pz_cut, prefix+"/plots/", "XiMinus_pt_vs_pz_cut", save, close);

		//gStyle->SetOptStat(0);

		TH1D * h_xi_m_nocut = new TH1D("h_xi_m_nocut", "Mass distribution for #Xi^{-}; M [GeV/c^{2}]; counts", 500,1.26,1.38);
		ntpXiMinus->Project("h_xi_m_nocut", "VtxFit_m", "McTruthMatch");

		TH1D * h_xi_m_vtxcut = new TH1D("h_xi_m_vtxcut", "Mass distribution for #Xi^{-} with vtxcut; M [GeV/c^{2}]; counts", 500,1.26,1.38);
		ntpXiMinus->Project("h_xi_m_vtxcut", "VtxFit_m", "McTruthMatch "+vtxcut);

		TH1D * h_xi_m_masscut = new TH1D("h_xi_m_masscut", "Mass distribution for #Xi^{-} with vertex cut and mass cut; M [GeV/c^{2}]; counts", 500,1.26,1.38);
		ntpXiMinus->Project("h_xi_m_masscut", "VtxFit_m", "McTruthMatch "+vtxcut+masscut);

		jenny::CreateDrawAndSaveNHistograms(h_xi_m_nocut, h_xi_m_vtxcut, h_xi_m_masscut, "mass window", "VtxFit_prob>0.01", "VtxFit_prob>0.01 && MassFit_prob>0.01", prefix+"/plots/", "XiMinus_m_diffcuts", save, close);

		//gStyle->SetOptStat(1111);

		TH1D * h_xi_m_masscut2 = new TH1D("h_xi_m_masscut2", "Mass distribution for #Xi^{-} with vertex cut and mass cut; M [GeV/c^{2}]; counts", 500,1.26,1.38);
		ntpXiMinus->Project("h_xi_m_masscut2", "VtxFit_m", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_xi_m_masscut2, prefix+"/plots/", "XiMinus_m_masscut", save, close);//, false, 0.01,0.1);

		TH1D * h_xi_m_mc = new TH1D("h_xi_m_mc", "Mass distribution for #Xi^{-} to MC truth; (M_{#Xi}-M_{PDG}) [MeV/c^{2}]; counts", 500,-20, 20);
		ntpXiMinus->Project("h_xi_m_mc", "(VtxFit_m-MCTruth_m)*1000", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_xi_m_mc, prefix+"/plots/", "XiMinus_m_diff_to_mc", save, close);//, false, 0.01,0.1);

		TH1D * h_xi_vtxres_x = new TH1D("h_xi_vtxres_x", " ; x-x_{MC} [mm]; counts", 500,-10,10);
		ntpXiMinus->Project("h_xi_vtxres_x", "VtxFit_diffvx*10", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_xi_vtxres_x, prefix+"/plots/", "XiMinus_vtxres_x", save, close);

		TH1D * h_xi_vtxres_y = new TH1D("h_xi_vtxres_y", " ; y-y_{MC} [mm]; counts", 500,-10,10);
		ntpXiMinus->Project("h_xi_vtxres_y", "VtxFit_diffvy*10", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_xi_vtxres_y, prefix+"/plots/", "XiMinus_vtxres_y", save, close);

		TH1D * h_xi_vtxres_z = new TH1D("h_xi_vtxres_z", " ; z-z_{MC} [mm]; counts", 500,-10,10);
		ntpXiMinus->Project("h_xi_vtxres_z", "VtxFit_diffvz*10", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_xi_vtxres_z, prefix+"/plots/", "XiMinus_vtxres_z", save, close);

		TH1D * h_xi_costht = new TH1D("h_xi_costht", "cos(#Theta) distribution for #Xi^{-}; cos(#Theta); counts", 500,-1,1.05);
		ntpXiMinus->Project("h_xi_costht", "cos(VtxFit_tht)","McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogram(h_xi_costht, prefix+"/plots/", "XiMinus_costht", save, close);

		TH1D * h_xi_tht = new TH1D("h_xi_tht", "#Theta distribution for #Xi^{-}; #Theta [rad]; counts", 500,0,0.5);
		ntpXiMinus->Project("h_xi_tht", "VtxFit_tht","McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogram(h_xi_tht, prefix+"/plots/", "XiMinus_tht", save, close);

		TH1D * h_xi_chisq = new TH1D("h_xi_chisq", "#chi^{2} distribution for #Xi^{-}; #chi^{2}; counts", 500,0,10);
		ntpXiMinus->Project("h_xi_chisq", "VtxFit_chisq","McTruthMatch ");
		jenny::CreateDrawAndSaveHistogram(h_xi_chisq, prefix+"/plots/", "XiMinus_chisq", save, close);

		TH1D * h_xi_prob = new TH1D("h_xi_prob", "probapility distribution for #Xi^{-}; prob; counts", 500,0,1);
		ntpXiMinus->Project("h_xi_prob", "VtxFit_prob","McTruthMatch ");
		jenny::CreateDrawAndSaveHistogram(h_xi_prob, prefix+"/plots/", "XiMinus_prob", save, close,true);
	}




	else if (ntpXiPlus !=0x0){

		//**** Get information about XiPlus ****************************************************************************

		TH2D * h_axi_pt_vs_pz = new TH2D("h_axi_pt_vs_pz", "Transverse vs. longitudinal momentum for #bar{#Xi}^{+}; p_{z} [GeV/c]; p_{t} [GeV/c]", 200,0.5,4,200,0,0.7);
		ntpXiPlus->Project("h_axi_pt_vs_pz", "MCTruth_pt: MCTruth_pz", "McTruthMatch");
		jenny::CreateDrawAndSaveHistogram(h_axi_pt_vs_pz, prefix+"/plots/", "XiPlus_pt_vs_pz", save, close);

		TH2D * h_axi_pt_vs_pz_cut = new TH2D("h_axi_pt_vs_pz_cut", "Transverse vs. longitudinal momentum for #bar{#Xi}^{+} after cut; p_{z} [GeV/c]; p_{t} [GeV/c]", 200,0.5,4,200,0,0.7);
		ntpXiPlus->Project("h_axi_pt_vs_pz_cut", "MCTruth_pt: MCTruth_pz", "McTruthMatch"+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogram(h_axi_pt_vs_pz_cut, prefix+"/plots/", "XiPlus_pt_vs_pz_cut", save, close);

		//gStyle->SetOptStat(0);
		TH1D * h_axi_m_nocut = new TH1D("h_axi_m_nocut", "Mass distribution for #bar{#Xi}^{+}; M [GeV/c^{2}]; counts", 500,1.26,1.38);
		ntpXiPlus->Project("h_axi_m_nocut", "VtxFit_m", "McTruthMatch");

		TH1D * h_axi_m_vtxcut = new TH1D("h_axi_m_vtxcut", "Mass distribution for #bar{#Xi}^{+} with vtxcut; M [GeV/c^{2}]; counts", 500,1.26,1.38);
		ntpXiPlus->Project("h_axi_m_vtxcut", "VtxFit_m", "McTruthMatch "+vtxcut);

		TH1D * h_axi_m_masscut = new TH1D("h_axi_m_masscut", "Mass distribution for #bar{#Xi}^{+} with vertex cut and mass cut; M [GeV/c^{2}]; counts", 500,1.26,1.38);
		ntpXiPlus->Project("h_axi_m_masscut", "VtxFit_m", "McTruthMatch "+vtxcut+masscut);

		jenny::CreateDrawAndSaveNHistograms(h_axi_m_nocut, h_axi_m_vtxcut, h_axi_m_masscut, "mass window", "VtxFit_prob>0.01", "VtxFit_prob>0.01 && MassFit_prob>0.01", prefix+"/plots/", "XiPlus_m_diffcuts", save, close);

		//gStyle->SetOptStat(1111);

		TH1D * h_axi_m_masscut2 = new TH1D("h_axi_m_masscut2", "Mass distribution for #bar{#Xi}^{+} with vertex cut and mass cut; M [GeV/c^{2}]; counts", 500,1.26,1.38);
		ntpXiPlus->Project("h_axi_m_masscut2", "VtxFit_m", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_axi_m_masscut2, prefix+"/plots/", "XiPlus_m_masscut", save, close);

		TH1D * h_axi_m_mc = new TH1D("h_axi_m_mc", "Mass distribution for #bar{#Xi}^{+} to MC truth; (M_{#bar{#Xi}}-M_{PDG}) [MeV/c^{2}]; counts", 500,-20, 20);
		ntpXiPlus->Project("h_axi_m_mc", "(VtxFit_m-MCTruth_m)*1000", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramDoulbeFit(h_axi_m_mc, prefix+"/plots/", "XiPlus_m_diff_to_mc", save, close);

		TH1D * h_axi_vtxres_x = new TH1D("h_axi_vtxres_x", "resolution for x coordinate of vertex for #bar{#Xi}^{+}; x-x_{MC} [mm]; counts", 500,-10,10);
		ntpXiPlus->Project("h_axi_vtxres_x", "VtxFit_diffvx*10", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_axi_vtxres_x, prefix+"/plots/", "XiPlus_vtxres_x", save, close);

		TH1D * h_axi_vtxres_y = new TH1D("h_axi_vtxres_y", "resolution for y coordinate of vertex for #bar{#Xi}^{+}; y-y_{MC} [mm]; counts", 500,-10,10);
		ntpXiPlus->Project("h_axi_vtxres_y", "VtxFit_diffvy*10", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_axi_vtxres_y, prefix+"/plots/", "XiPlus_vtxres_y", save, close);

		TH1D * h_axi_vtxres_z = new TH1D("h_axi_vtxres_z", "resolution for z coordinate of vertex for #bar{#Xi}^{+}; z-z_{MC} [mm]; counts", 500,-20,20);
		ntpXiPlus->Project("h_axi_vtxres_z", "VtxFit_diffvz*10", "McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogramFWHM(h_axi_vtxres_z, prefix+"/plots/", "XiPlus_vtxres_z", save, close);


		TH1D * h_axi_costht = new TH1D("h_axi_costht", "cos(#Theta) distribution for #bar{#Xi}^{+}; cos(#Theta); counts", 500,-1,1.05);
		ntpXiPlus->Project("h_axi_costht", "cos(VtxFit_tht)","McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogram(h_axi_costht, prefix+"/plots/", "XiPlus_costht", save, close);

		TH1D * h_axi_tht = new TH1D("h_axi_tht", "#Theta distribution for #bar{#Xi}^{+}; #Theta [rad]; counts", 500,0,0.5);
		ntpXiPlus->Project("h_axi_tht", "VtxFit_tht","McTruthMatch "+vtxcut+masscut);
		jenny::CreateDrawAndSaveHistogram(h_axi_tht, prefix+"/plots/", "XiPlus_tht", save, close);


		TH1D * h_xi_chisq = new TH1D("h_xi_chisq", "#chi^{2} distribution for #bar{#Xi}^{+}; #chi^{2}; counts", 500,0,3);
		ntpXiPlus->Project("h_xi_chisq", "VtxFit_chisq","McTruthMatch ");
		jenny::CreateDrawAndSaveHistogram(h_xi_chisq, prefix+"/plots/", "XiPlus_chisq", save, close);

		TH1D * h_xi_prob = new TH1D("h_xi_prob", "probability distribution for #bar{#Xi}^{+}; prob; counts", 500,0,1);
		ntpXiPlus->Project("h_xi_prob", "VtxFit_prob","McTruthMatch ");
		jenny::CreateDrawAndSaveHistogram(h_xi_prob, prefix+"/plots/", "XiPlus_prob", save, close, true);
	}
	else {
		cout << "No particle of kind Xi- or Xi+!" << endl;
	}
}
