/**
 * @file vertex_studies_lambda0.C
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Creates histograms for vertex studies for Lambda0 - AntiLambda production
 * @details Methode creates Histogramms for vertex studies for Lambda0 - AntiLambda production. Data is comming from the analysis which is stored in a file called output_ana.root
 * Basically ROOT.
 */


class RhoTuple;

#include "../common_jenny.cpp"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"


void vertex_studies_pbarp_to_2Lamda0_kalman(int nevts=0, bool saveoutput=true, bool close=false){

	//Paths
	TString inPath = "/private/puetz/mysimulations/analysis/pbarp_lambda0_antilambda0/500000_events/";
	TString outPath = inPath + "plots/";


	//Input file
	TFile * file = new TFile(inPath + "ana_total.root", "READ");

	//Reading data from file
//	TTree * ntpMC = (TTree*) file->Get("ntpMC");
	TTree * ntpBeam = (TTree*) file->Get("fntpCrossCheck");

	//***Projection for the Beam

	//simulated vertex position for Beam
//
//	TH1D * h_vx_mc = new TH1D("h_vx_mc", "x coordinate of vertex position; x[cm]; counts", 85,-5,5);
//	ntpMC->Project("h_vx_mc", "x", "moth==0 && (pdg==3122 || pdg==-3122) ");
//
//	TH1D * h_vy_mc = new TH1D("h_vy_mc", "y coordinate of vertex position; y[cm]; counts", 85,-5,5);
//	ntpMC->Project("h_vy_mc", "y", "moth==0 && (pdg==3122 || pdg==-3122)");
//
//	TH1D * h_vz_mc = new TH1D("h_vz_mc", "z coordinate of vertex position; z[cm]; counts", 85,0,5);
//	ntpMC->Project("h_vz_mc", "z", "moth==0 && (pdg==3122 || pdg==-3122)");


	//reconstructed vertex position for Beam
	TH1D * h_vx_reco = new TH1D("h_vx_reco", "x coordinate of vertex position; x[cm]; counts", 58,-5,5);
	ntpBeam->Project("h_vx_reco", "KalmanFit_vx", "McTruthMatch==1");

	TH1D * h_vy_reco = new TH1D("h_vy_reco", "y coordinate of vertex position; y[cm]; counts", 58,-5,5);
	ntpBeam->Project("h_vy_reco", "KalmanFit_vy", "McTruthMatch==1");

	TH1D * h_vz_reco = new TH1D("h_vz_reco", "z coordinate of vertex position; z[cm]; counts", 58,-5,5);
	ntpBeam->Project("h_vz_reco", "KalmanFit_vz", "McTruthMatch==1");


	//resolution vertex position for Beam
	TH1D * h_vx_res = new TH1D("h_vx_res", "x resolution of vertex position; #Delta x[cm]; counts", 100,-1,1);
	ntpBeam->Project("h_vx_res", "KalmanFit_vx-McTruth_vx", "McTruthMatch==1");

	TH1D * h_vy_res = new TH1D("h_vy_res", "y resolution of vertex position; #Delta y[cm]; counts", 100,-1,1);
	ntpBeam->Project("h_vy_res", "KalmanFit_vy-McTruth_vy", "McTruthMatch==1");

	TH1D * h_vz_res = new TH1D("h_vz_res", "z resolution of vertex position; #Delta z[cm]; counts", 100,-1.5,1.5);
	ntpBeam->Project("h_vz_res", "KalmanFit_vz-McTruth_vz", "McTruthMatch==1");



	//resolution vertex position and poca
	TH1D * h_vx_poca = new TH1D("h_vx_poca", "vertex position - poca for x coordinate; x_{reco} - x_{poca}[cm]; counts", 100,-1,1);
	ntpBeam->Project("h_vx_poca", "KalmanFit_vx-pocvx", "McTruthMatch==1");

	TH1D * h_vy_poca = new TH1D("h_vy_poca", "vertex position - poca for y coordinate; y_{reco} - y_{poca}[cm]; counts", 100,-1,1);
	ntpBeam->Project("h_vy_poca", "KalmanFit_vy-pocvy", "McTruthMatch==1");

	TH1D * h_vz_poca = new TH1D("h_vz_poca", "vertex position - poca for z coordinate; z_{reco} - z_{poca}[cm]; counts", 100,-1.5,1.5);
	ntpBeam->Project("h_vz_poca", "KalmanFit_vz-pocvz", "McTruthMatch==1");


	//2-dimensional plots for poca and vertex

	TH2D * h_vxy_vz_reco = new TH2D("h_vxy_vz_reco", "Vertex ; z[cm]; r[cm]", 58,0,5,58,0,5);
	ntpBeam->Project("h_vxy_vz_reco", "sqrt(KalmanFit_vx**2+KalmanFit_vy**2):KalmanFit_vz", "McTruthMatch==1");


	//Goodness of fit

	TH1D * h_chi2 = new TH1D("h_chi2", "#chi^{2} distribution for vertex Fit; #chi^{2}; counts", 58,0,10);
	ntpBeam->Project("h_chi2", "KalmanFit_chisq", "McTruthMatch==1");

	TH1D * h_Prob = new TH1D("h_Prob", "Probability for vertex Fit; Prob; counts", 58,0,1.1);
	ntpBeam->Project("h_Prob", "KalmanFit_prob", "McTruthMatch==1");


	//**********Performance plots******************
	//pulls for vertex

	TH1D * h_pull_vx = new TH1D("h_pull_vx", "pull for x coordinate of vertex position; (x_{reco}-x_{MC})/#sigma_{x}; counts", 1000, -3,3);
	ntpBeam->Project("h_pull_vx", "KalmanFit_pullvx", "McTruthMatch==1");

	TH1D * h_pull_vy = new TH1D("h_pull_vy", "pull for y coordinate of vertex position; (y_{reco}-y_{MC})/#sigma_{y}; counts", 1000, -3,3);
	ntpBeam->Project("h_pull_vy", "KalmanFit_pullvy", "McTruthMatch==1");

	TH1D * h_pull_vz = new TH1D("h_pull_vz", "pull for z coordinate of vertex position; (z_{reco}-z_{MC})/#sigma_{z}; counts", 1000, -3,3);
	ntpBeam->Project("h_pull_vz", "KalmanFit_pullvz", "McTruthMatch==1");


	//vtx res vs p_t

	TH2D * h_vx_vs_pt = new TH2D("h_vx_vs_pt", " transversal momentum vs. vertex resolution(x coordinate); #Delta x/cm; pt/GeV", 200, -1,1, 200, 0,0.4);
	ntpBeam->Project("h_vx_vs_pt", "KalmanFit_pt:(KalmanFit_vx-McTruth_vx)", "McTruthMatch==1");
	TH1D * h_pt_vx = h_vx_vs_pt->ProfileY("h_pt_vx", 1, -1, "o");
	h_pt_vx->SetTitle("vertex resolution(x coordinate) vs. transversam momentum");
	h_pt_vx->GetYaxis()->SetTitle("#Delta x/cm");
	h_pt_vx->GetYaxis()->SetRangeUser(-0.8,0.8);

	TH2D * h_vy_vs_pt = new TH2D("h_vy_vs_pt", "transversal momentum vs. vertex resolution (y coordinate); #Delta y/cm; p_{t}/GeV", 200, -1,1, 200, 0,0.4);
	ntpBeam->Project("h_vy_vs_pt", "KalmanFit_pt:(KalmanFit_vy-McTruth_vy)", "McTruthMatch==1");
	TH1D * h_pt_vy = h_vy_vs_pt->ProfileY("h_pt_vy", 1, -1, "o");
	h_pt_vy->SetTitle("vertex resolution(y coordinate) vs. transversam momentum");
	h_pt_vy->GetYaxis()->SetTitle("#Delta y/cm");
	h_pt_vy->GetYaxis()->SetRangeUser(-0.8,0.8);

	TH2D * h_vz_vs_pt = new TH2D("h_vz_vs_pt", "transversal momentum vs. vertex resolution (z coordinate); #Delta z/cm; p_{t}/GeV", 200, -1,1, 200, 0,0.4);
	ntpBeam->Project("h_vz_vs_pt", "KalmanFit_pt:(KalmanFit_vz-McTruth_vz)", "McTruthMatch==1");
	TH1D * h_pt_vz = h_vz_vs_pt->ProfileY("h_pt_vz", 1, -1, "o");
	h_pt_vz->SetTitle("vertex resolution(z coordinate) vs. transversam momentum");
	h_pt_vz->GetYaxis()->SetTitle("#Delta z/cm");
	h_pt_vz->GetYaxis()->SetRangeUser(-0.8,0.8);

	//vtx res vs p_z

	TH2D * h_vx_vs_pz = new TH2D("h_vx_vs_pz", "longitudinal momentum vs. vertex resolution(x coordinate); #Delta x/cm; p_{z}/GeV", 200, -1,1, 200, 0.5,2.5);
	ntpBeam->Project("h_vx_vs_pz", "KalmanFit_pz:(KalmanFit_vx-McTruth_vx)", "McTruthMatch==1");
	TH2D * h_vy_vs_pz = new TH2D("h_vy_vs_pz", "longitudinal momentum vs. vertex resolution(y coordinate); #Delta y/cm; p_{z}/GeV", 200, -1,1, 200, 0.5,2.5);
	ntpBeam->Project("h_vy_vs_pz", "KalmanFit_pz:(KalmanFit_vy-McTruth_vy)", "McTruthMatch==1");
	TH2D * h_vz_vs_pz = new TH2D("h_vz_vs_pz", "longitudinal momentum vs. vertex resolution(z coordinate); #Delta z/cm; p_{z}/GeV", 200, -1,1, 200, 0.5,2.5);
	ntpBeam->Project("h_vz_vs_pz", "KalmanFit_pz:(KalmanFit_vz-McTruth_vz)", "McTruthMatch==1");

	//vtx res vs p

	TH2D * h_vx_vs_p = new TH2D("h_vx_vs_p", " momentum vs. vertex resolution(x coordinate); #Delta x/cm; p/GeV", 200, -1,1, 200, 0.5,2.5);
	ntpBeam->Project("h_vx_vs_p", "KalmanFit_p:(KalmanFit_vx-McTruth_vx)", "McTruthMatch==1");
	TH2D * h_vy_vs_p = new TH2D("h_vy_vs_p", " momentum vs. vertex resolution(y coordinate); #Delta y/cm; p/GeV", 200, -1,1, 200, 0.5,2.5);
	ntpBeam->Project("h_vy_vs_p", "KalmanFit_p:(KalmanFit_vy-McTruth_vy)", "McTruthMatch==1");
	TH2D * h_vz_vs_p = new TH2D("h_vz_vs_p", " momentum vs. vertex resolution(z coordinate); #Delta z/cm; p/GeV", 200, -1,1, 200, 0.5,2.5);
	ntpBeam->Project("h_vz_vs_p", "KalmanFit_p:(KalmanFit_vz-McTruth_vz)", "McTruthMatch==1");


	//vtx res vs tht

	TH2D * h_vx_vs_tht = new TH2D("h_vx_vs_tht", " Theta vs. vertex resolution(x coordinate); #Delta x/cm; tht/rad", 200, -1,1, 200, 0,0.15);
	ntpBeam->Project("h_vx_vs_tht", "KalmanFit_tht:(KalmanFit_vx-McTruth_vx)", "McTruthMatch==1");
	TH1D * h_tht_vx = h_vx_vs_tht->ProfileY("h_tht_vx", 1, -1, "o");
	h_tht_vx->SetTitle("vertex resolution(x coordinate) vs. Theta");
	h_tht_vx->GetYaxis()->SetTitle("#Delta x/cm");

	TH2D * h_vy_vs_tht = new TH2D("h_vy_vs_tht", " Theta vs. vertex resolution(y coordinate); #Delta y/cm; tht/rad", 200, -1,1, 200, 0,0.15);
	ntpBeam->Project("h_vy_vs_tht", "KalmanFit_tht:(KalmanFit_vy-McTruth_vy)", "McTruthMatch==1");
	TH1D * h_tht_vy = h_vy_vs_tht->ProfileY("h_tht_vy", 1, -1, "o");
	h_tht_vy->SetTitle("vertex resolution(y coordinate) vs. Theta");
	h_tht_vy->GetYaxis()->SetTitle("#Delta y/cm");

	TH2D * h_vz_vs_tht = new TH2D("h_vz_vs_tht", " Theta vs. vertex resolution(z coordinate); #Delta z/cm; tht/rad", 200, -1,1, 200, 0,0.15);
	ntpBeam->Project("h_vz_vs_tht", "KalmanFit_tht:(KalmanFit_vz-McTruth_vz)", "McTruthMatch==1");
	TH1D * h_tht_vz = h_vz_vs_tht->ProfileY("h_tht_vz", 1, -1, "o");
	h_tht_vz->SetTitle("vertex resolution(z coordinate) vs. Theta");
	h_tht_vz->GetYaxis()->SetTitle("#Delta z/cm");


	//*** create, draw and save histogram

//	jenny::CreateDrawAndSaveHistogram(h_vx_mc, outPath, "h_vx_mc_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_mc, outPath, "h_vy_mc_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_mc, outPath, "h_vz_mc_kalman", saveoutput, close);
//
//	jenny::CreateDrawAndSaveHistogram(h_vx_reco, outPath, "h_vx_reco_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_reco, outPath, "h_vy_reco_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_reco, outPath, "h_vz_reco_kalman", saveoutput, close);

	bool autoRange = false;

	jenny::CreateDrawAndSaveHistogramWithFit(h_vx_res, outPath, "h_vx_res_kalman", saveoutput, close, autoRange, 0.1, 1, true);
	jenny::CreateDrawAndSaveHistogramWithFit(h_vy_res, outPath, "h_vy_res_kalman", saveoutput, close, autoRange, 0.1, 1, true);
	jenny::CreateDrawAndSaveHistogramWithFit(h_vz_res, outPath, "h_vz_res_kalman", saveoutput, close, autoRange, 0.05, 1.5, true);

//	jenny::CreateDrawAndSaveHistogramDoubleFit(h_vx_poca, outPath, "h_vx_poca_kalman", saveoutput, close);//, autoRange, 0.1, 1);
//	jenny::CreateDrawAndSaveHistogramDoubleFit(h_vy_poca, outPath, "h_vy_poca_kalman", saveoutput, close);//, autoRange, 0.6, 1);
//	jenny::CreateDrawAndSaveHistogramDoubleFit(h_vz_poca, outPath, "h_vz_poca_kalman", saveoutput, close);//, autoRange, 0.5, 15);
//
//
//	jenny::CreateDrawAndSaveHistogram(h_vxy_vz_reco, outPath, "h_vxy_vz_reco_kalman", saveoutput, close);
//
//
	jenny::CreateDrawAndSaveHistogram(h_chi2, outPath, "h_chi2_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_Prob, outPath, "h_Prob_kalman", saveoutput, close);

//	jenny::CreateDrawAndSaveHistogramWithFit(h_pull_vx, outPath, "h_pull_vx_kalman", saveoutput, close, autoRange, 0.15, 3);
//	jenny::CreateDrawAndSaveHistogramWithFit(h_pull_vy, outPath, "h_pull_vy_kalman", saveoutput, close, autoRange, 0.2, 3, true);
//	jenny::CreateDrawAndSaveHistogramWithFit(h_pull_vz, outPath, "h_pull_vz_kalman", saveoutput, close, autoRange, 0.4, 3, true);

	jenny::CreateDrawAndSaveHistogram(h_vx_vs_pt, outPath, "h_vx_vs_pt_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_pt_vx, outPath, "h_pt_vx_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vy_vs_pt, outPath, "h_vy_vs_pt_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_pt_vy, outPath, "h_pt_vy_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vz_vs_pt, outPath, "h_vz_vs_pt_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_pt_vz, outPath, "h_pt_vz_kalman", saveoutput, close);

	jenny::CreateDrawAndSaveHistogram(h_vx_vs_pz, outPath, "h_vx_vs_pz_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vy_vs_pz, outPath, "h_vy_vs_pz_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vz_vs_pz, outPath, "h_vz_vs_pz_kalman", saveoutput, close);

	jenny::CreateDrawAndSaveHistogram(h_vx_vs_p, outPath, "h_vx_vs_p_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vy_vs_p, outPath, "h_vy_vs_p_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vz_vs_p, outPath, "h_vz_vs_p_kalman", saveoutput, close);

	jenny::CreateDrawAndSaveHistogram(h_vx_vs_tht, outPath, "h_vx_vs_tht_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_tht_vx, outPath, "h_tht_vx_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vy_vs_tht, outPath, "h_vy_vs_tht_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_tht_vy, outPath, "h_tht_vy_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_vz_vs_tht, outPath, "h_vz_vs_tht_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_tht_vz, outPath, "h_tht_vz_kalman", saveoutput, close);



	if(close) exit(0);

}
