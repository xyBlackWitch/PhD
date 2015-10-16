/**
 * @file vertex_studies_lambda0.C
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Creates histograms for vertex studies of Lambda0 decay.
 * @details Methode creates Histogramms for vertex studies of a Lambda0 decay. Data is comming from the analysis which is stored in a file called output_ana.root
 * Basically ROOT.
 */


class RhoTuple;

#include <vector>
#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "../common_jenny.cpp"





void vertex_studies_lambda0_kalman(int nevts=0, bool saveoutput=true, bool close=false){


	//Paths

	TString inPath = "/home/ikp1/puetz/panda/mysimulations/test/boxgenerator/lambda0/500000_events/";
	TString outPath = inPath +"plots/";


	//Input File

	TString inputFile = inPath + "ana_total.root";

	//open input
	TFile * file = new TFile(inputFile, "READ");

	TTree * fntpLambda = (TTree*) file->Get("fntpLambda0");
//	TTree * fntpMc = (TTree*) file->Get("fntpMc");


	//****Projection for Lambda0


//	//simulated vertex position for lambda0
//
//	TH1D * h_vx_Mc = new TH1D("h_vx_Mc", "x coodinate of vertex position; x[cm]; counts", 100,-10,10);
//	fntpMc->Project("h_vx_Mc", "x", "moth==-1");
//
//	TH1D * h_vy_Mc = new TH1D("h_vy_Mc", "y coodinate of vertex position; y[cm]; counts", 100,-10,10);
//	fntpMc->Project("h_vy_Mc", "y", "moth==-1");
//
//	TH1D * h_vz_Mc = new TH1D("h_vz_Mc", "z coodinate of vertex position; z[cm]; counts", 100,-10,10);
//	fntpMc->Project("h_vz_Mc", "z", "moth==-1");

	//Truth information

//	TH1D * h_vx_truth = new TH1D("h_vx_truth", "x coodinate of vertex position; x[cm]; counts", 1000,-100,100);
//	fntpLambda0->Project("h_vx_truth", "McTruth_vx", "McTruthMatch==1");
//
//	TH1D * h_vy_truth = new TH1D("h_vy_truth", "y coodinate of vertex position; y[cm]; counts", 100,-100,100);
//	fntpLambda0->Project("h_vy_truth", "McTruth_vy", "McTruthMatch==1");
//
//	TH1D * h_vz_truth = new TH1D("h_vz_truth", "z coodinate of vertex position; z[cm]; counts", 100,-100,100);
//	fntpLambda0->Project("h_vz_truth", "McTruth_vz", "McTruthMatch==1");


	//reconstructed vertex position for lambda0
	TH1D * h_vx_reco = new TH1D("h_vx_reco", "x coordinate of reconstructed vertex position; x[cm]; counts", 1000,-100,100);
	fntpLambda->Project("h_vx_reco", "KalmanFit_vx", "McTruthMatch==1");

	TH1D * h_vy_reco = new TH1D("h_vy_reco", "y coordinate of reconstructed vertex position; y[cm]; counts", 1000,-100,100);
	fntpLambda->Project("h_vy_reco", "KalmanFit_vy", "McTruthMatch==1");

	TH1D * h_vz_reco = new TH1D("h_vz_reco", "z coordinate of reconstructed vertex position; z[cm]; counts", 1000,-100,100);
	fntpLambda->Project("h_vz_reco", "KalmanFit_vz", "McTruthMatch==1");



	//resolution for reconstructed vertex position for lambda0

	TH1D * h_vx_res = new TH1D("h_vx_res", "x resolution of reconstructed vertex position; #Delta x[cm]; counts", 1000,-1,1);
	fntpLambda->Project("h_vx_res", "KalmanFit_vx-McTruth_vx", "McTruthMatch==1");

	TH1D * h_vy_res = new TH1D("h_vy_res", "y resolution of reconstructed vertex position; #Delta y[cm]; counts", 1000,-1,1);
	fntpLambda->Project("h_vy_res", "KalmanFit_vy-McTruth_vy", "McTruthMatch==1");

	TH1D * h_vz_res = new TH1D("h_vz_res", "z resolution of reconstructed vertex position; #Delta z[cm]; counts", 1000,-1,1);
	fntpLambda->Project("h_vz_res", "KalmanFit_vz-McTruth_vz", "McTruthMatch==1");


	//resolution for reconstructed vertex position and poca for lambda0

	TH1D * h_vx_poc_res = new TH1D("h_vx_poc_res", "differenz reconstructed vertex position and poca position; x_{reco} - x_{poca}[cm]; counts", 100,-1,1);
	fntpLambda->Project("h_vx_poc_res", "KalmanFit_vx-Lambda0_pocvx", "McTruthMatch==1");

	TH1D * h_vy_poc_res = new TH1D("h_vy_poc_res", "differenz reconstructed vertex position and poca position; y_{reco} - y_{poca}[cm]; counts", 100,-1,1);
	fntpLambda->Project("h_vy_poc_res", "KalmanFit_vy-Lambda0_pocvy", "McTruthMatch==1");

	TH1D * h_vz_poc_res = new TH1D("h_vz_poc_res", "differenz reconstructed vertex position and poca position; z_{reco} - z_{poca}[cm]; counts", 100,-1,1);
	fntpLambda->Project("h_vz_poc_res", "KalmanFit_vz-Lambda0_pocvz", "McTruthMatch==1");


//	//Combined histograms
//
//	TH1D * h_reco_div_Mc_vx = new TH1D("h_reco_div_Mc_vx", "Ratio Reco/Mc for the x coordinate of vertex position; x_{reco}/x_{Mc}", 100,-100,100);
//	h_reco_div_Mc_vx->Divide(h_vx_reco,h_vx_truth);
//	h_reco_div_Mc_vx->SetStats(0);
////	h_reco_div_Mc_vx->GetYaxis()->SetRangeUser(0,1.1);
//
//	TH1D * h_reco_div_Mc_vy = new TH1D("h_reco_div_Mc_vy", "Ratio Reco/Mc for the y coordinate of vertex position; y_{reco}/y_{Mc}", 100,-100,100);
//	h_reco_div_Mc_vy->Divide(h_vy_reco,h_vy_truth);
//	h_reco_div_Mc_vy->SetStats(0);
////	h_reco_div_Mc_vy->GetYaxis()->SetRangeUser(0,1.1);
//
//	TH1D * h_reco_div_Mc_vz = new TH1D("h_reco_div_Mc_vz", "Ratio Reco/Mc for the z coordinate of vertex position; z_{reco}/z_{Mc}", 100,0,100);
//	h_reco_div_Mc_vz->Divide(h_vz_reco,h_vz_truth);
//	h_reco_div_Mc_vz->SetStats(0);
////	h_reco_div_Mc_vz->GetYaxis()->SetRangeUser(0,1.1);


	//2-dimensional plots for poca and vertex

	TH2D * h_vxy_vz_reco = new TH2D("h_vxy_vz_reco", "Vertex ; z[cm]; r[cm]", 100,0,5,100,0,5);
	fntpLambda->Project("h_vxy_vz_reco", "sqrt(KalmanFit_vx**2+KalmanFit_vy**2):KalmanFit_vz", "McTruthMatch==1");


	//Goodness of fit

	TH1D * h_chi2 = new TH1D("h_chi2", "#chi^{2} distribution for vertex Fit; #chi^{2}; counts", 1000,0,10);
	fntpLambda->Project("h_chi2", "KalmanFit_chisq", "McTruthMatch==1");

	TH1D * h_Prob = new TH1D("h_Prob", "Probability for vertex Fit; Prob; counts", 1000,0,1);
	fntpLambda->Project("h_Prob", "KalmanFit_prob", "McTruthMatch==1");



	//**********Performance plots******************
	//pulls for vertex

	TH1D * h_pull_vx = new TH1D("h_pull_vx", "pull for x coordinate of momentum; (px_{reco}-px_{Mc})/#sigma_{px}; counts", 500, -5,5);
	fntpLambda0->Project("h_pull_vx", "KalmanFit_mcpullpx", "McTruthMatch==1");

	TH1D * h_pull_vy = new TH1D("h_pull_vy", "pull for y coordinate of momentum; (py_{reco}-py_{Mc})/#sigma_{py}; counts", 500, -5,5);
	fntpLambda0->Project("h_pull_vy", "KalmanFit_mcpullpy", "McTruthMatch==1");

	TH1D * h_pull_vz = new TH1D("h_pull_vz", "pull for z coordinate of momentum; (pz_{reco}-pz_{Mc})/#sigma_{pz}; counts", 500, -5,5);
	fntpLambda0->Project("h_pull_vz", "KalmanFit_mcpullpz", "McTruthMatch==1");


	//vtx res vs p_t

	TH2D * h_vx_vs_pt = new TH2D("h_vx_vs_pt", " transversal momentum vs. vertex resolution(x coordinate); #Delta x/cm; pt/GeV", 200, -1,1, 200, 0,2);
	fntpLambda0->Project("h_vx_vs_pt", "KalmanFit_pt:(KalmanFit_vx-McTruth_vx)", "McTruthMatch==1");
	TH1D * h_pt_vx = h_vx_vs_pt->ProfileY("h_pt_vx", 1, -1, "o");
	h_pt_vx->SetTitle("vertex resolution(x coordinate) vs. transversam momentum");
	h_pt_vx->GetYaxis()->SetTitle("#Delta x/cm");
	h_pt_vx->GetYaxis()->SetRangeUser(-0.8,0.8);

	TH2D * h_vy_vs_pt = new TH2D("h_vy_vs_pt", "transversal momentum vs. vertex resolution (y coordinate); #Delta y/cm; p_{t}/GeV", 200, -1,1, 200, 0,2);
	fntpLambda0->Project("h_vy_vs_pt", "KalmanFit_pt:(KalmanFit_vy-McTruth_vy)", "McTruthMatch==1");
	TH1D * h_pt_vy = h_vy_vs_pt->ProfileY("h_pt_vy", 1, -1, "o");
	h_pt_vy->SetTitle("vertex resolution(y coordinate) vs. transversam momentum");
	h_pt_vy->GetYaxis()->SetTitle("#Delta y/cm");
	h_pt_vy->GetYaxis()->SetRangeUser(-0.8,0.8);

	TH2D * h_vz_vs_pt = new TH2D("h_vz_vs_pt", "transversal momentum vs. vertex resolution (z coordinate); #Delta z/cm; p_{t}/GeV", 200, -1,1, 200, 0,2);
	fntpLambda0->Project("h_vz_vs_pt", "KalmanFit_pt:(KalmanFit_vz-McTruth_vz)", "McTruthMatch==1");
	TH1D * h_pt_vz = h_vz_vs_pt->ProfileY("h_pt_vz", 1, -1, "o");
	h_pt_vz->SetTitle("vertex resolution(z coordinate) vs. transversam momentum");
	h_pt_vz->GetYaxis()->SetTitle("#Delta z/cm");
	h_pt_vz->GetYaxis()->SetRangeUser(-0.8,0.8);

	//vtx res vs p_z

	TH2D * h_vx_vs_pz = new TH2D("h_vx_vs_pz", "longitudinal momentum vs. vertex resolution(x coordinate); #Delta x/cm; p_{z}/GeV", 200, -1,1, 200, 0,3);
	fntpLambda0->Project("h_vx_vs_pz", "KalmanFit_pz:(KalmanFit_vx-McTruth_vx)", "McTruthMatch==1");
	TH2D * h_vy_vs_pz = new TH2D("h_vy_vs_pz", "longitudinal momentum vs. vertex resolution(y coordinate); #Delta y/cm; p_{z}/GeV", 200, -1,1, 200, 0,3);
	fntpLambda0->Project("h_vy_vs_pz", "KalmanFit_pz:(KalmanFit_vy-McTruth_vy)", "McTruthMatch==1");
	TH2D * h_vz_vs_pz = new TH2D("h_vz_vs_pz", "longitudinal momentum vs. vertex resolution(z coordinate); #Delta z/cm; p_{z}/GeV", 200, -1,1, 200, 0,3);
	fntpLambda0->Project("h_vz_vs_pz", "KalmanFit_pz:(KalmanFit_vz-McTruth_vz)", "McTruthMatch==1");

	//vtx res vs p_z

	TH2D * h_vx_vs_p = new TH2D("h_vx_vs_p", " momentum vs. vertex resolution(x coordinate); #Delta x/cm; p/GeV", 200, -1,1, 200, 0,5);
	fntpLambda0->Project("h_vx_vs_p", "KalmanFit_p:(KalmanFit_vx-McTruth_vx)", "McTruthMatch==1");
	TH2D * h_vy_vs_p = new TH2D("h_vy_vs_p", " momentum vs. vertex resolution(y coordinate); #Delta y/cm; p/GeV", 200, -1,1, 200, 0,5);
	fntpLambda0->Project("h_vy_vs_p", "KalmanFit_p:(KalmanFit_vy-McTruth_vy)", "McTruthMatch==1");
	TH2D * h_vz_vs_p = new TH2D("h_vz_vs_p", " momentum vs. vertex resolution(z coordinate); #Delta z/cm; p/GeV", 200, -1,1, 200, 0,5);
	fntpLambda0->Project("h_vz_vs_p", "KalmanFit_p:(KalmanFit_vz-McTruth_vz)", "McTruthMatch==1");


	//vtx res vs tht

	TH2D * h_vx_vs_tht = new TH2D("h_vx_vs_tht", " Theta vs. vertex resolution(x coordinate); #Delta x/cm; tht/rad", 200, -1,1, 200, 0,0.6);
	fntpLambda0->Project("h_vx_vs_tht", "KalmanFit_tht:(KalmanFit_vx-McTruth_vx)", "McTruthMatch==1");
	TH1D * h_tht_vx = h_vx_vs_tht->ProfileY("h_tht_vx", 1, -1, "o");
	h_tht_vx->SetTitle("vertex resolution(x coordinate) vs. Theta");
	h_tht_vx->GetYaxis()->SetTitle("#Delta x/cm");

	TH2D * h_vy_vs_tht = new TH2D("h_vy_vs_tht", " Theta vs. vertex resolution(y coordinate); #Delta y/cm; tht/rad", 200, -1,1, 200, 0,0.6);
	fntpLambda0->Project("h_vy_vs_tht", "KalmanFit_tht:(KalmanFit_vy-McTruth_vy)", "McTruthMatch==1");
	TH1D * h_tht_vy = h_vy_vs_tht->ProfileY("h_tht_vy", 1, -1, "o");
	h_tht_vy->SetTitle("vertex resolution(y coordinate) vs. Theta");
	h_tht_vy->GetYaxis()->SetTitle("#Delta y/cm");

	TH2D * h_vz_vs_tht = new TH2D("h_vz_vs_tht", " Theta vs. vertex resolution(z coordinate); #Delta z/cm; tht/rad", 200, -1,1, 200, 0,0.6);
	fntpLambda0->Project("h_vz_vs_tht", "KalmanFit_tht:(KalmanFit_vz-McTruth_vz)", "McTruthMatch==1");
	TH1D * h_tht_vz = h_vz_vs_tht->ProfileY("h_tht_vz", 1, -1, "o");
	h_tht_vz->SetTitle("vertex resolution(z coordinate) vs. Theta");
	h_tht_vz->GetYaxis()->SetTitle("#Delta z/cm");



	//***Create Canvas and draw histogramm
//
//	jenny::CreateDrawAndSaveHistogram(h_vx_Mc, outPath, "h_vx_Mc_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_Mc, outPath, "h_vy_Mc_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_Mc, outPath, "h_vz_Mc_kalman", saveoutput, close);

	bool autoRange = false;
//
//	jenny::CreateDrawAndSaveHistogramWithFit(h_vx_res, outPath, "h_vx_res_kalman", saveoutput, close, autoRange, 0.1,1, true);
//	jenny::CreateDrawAndSaveHistogramWithFit(h_vy_res, outPath, "h_vy_res_kalman", saveoutput, close, autoRange, 0.1,1, true);
//	jenny::CreateDrawAndSaveHistogramWithFit(h_vz_res, outPath, "h_vz_res_kalman", saveoutput, close, autoRange, 0.1, 1, true);

//	jenny::CreateDrawAndSaveHistogramDoubleFit(h_vx_poc_res, outPath, "h_vx_poc_res_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogramDoubleFit(h_vy_poc_res, outPath, "h_vy_poc_res_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogramDoubleFit(h_vz_poc_res, outPath, "h_vz_poc_res_kalman", saveoutput, close);
//
//	jenny::CreateDrawAndSaveHistogram(h_vx_reco, outPath, "h_vx_reco_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_reco, outPath, "h_vy_reco_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_reco, outPath, "h_vz_reco_kalman", saveoutput, close);
//
//	jenny::CreateDrawAndSaveHistogram(h_reco_div_Mc_vx, outPath, "h_reco_div_Mc_vx_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_reco_div_Mc_vy, outPath, "h_reco_div_Mc_vy_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_reco_div_Mc_vz, outPath, "h_reco_div_Mc_vz_kalman", saveoutput, close);
//
//
//	jenny::CreateDrawAndSaveHistogram(h_vxy_vz_reco, outPath, "h_vxy_vz_reco_kalman", saveoutput, close);
//
//
//	jenny::CreateDrawAndSaveHistogram(h_chi2, outPath, "h_chi2_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_Prob, outPath, "h_Prob_kalman", saveoutput, close);

	jenny::CreateDrawAndSaveHistogram(h_pull_vx, outPath, "h_pull_vx_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_pull_vy, outPath, "h_pull_vy_kalman", saveoutput, close);
	jenny::CreateDrawAndSaveHistogram(h_pull_vz, outPath, "h_pull_vz_kalman", saveoutput, close);

//	jenny::CreateDrawAndSaveHistogram(h_vx_vs_pt, outPath, "h_vx_vs_pt_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_pt_vx, outPath, "h_pt_vx_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_vs_pt, outPath, "h_vy_vs_pt_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_pt_vy, outPath, "h_pt_vy_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_vs_pt, outPath, "h_vz_vs_pt_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_pt_vz, outPath, "h_pt_vz_kalman", saveoutput, close);
//
//	jenny::CreateDrawAndSaveHistogram(h_vx_vs_pz, outPath, "h_vx_vs_pz_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_vs_pz, outPath, "h_vy_vs_pz_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_vs_pz, outPath, "h_vz_vs_pz_kalman", saveoutput, close);
//
//	jenny::CreateDrawAndSaveHistogram(h_vx_vs_p, outPath, "h_vx_vs_p_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_vs_p, outPath, "h_vy_vs_p_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_vs_p, outPath, "h_vz_vs_p_kalman", saveoutput, close);
//
//	jenny::CreateDrawAndSaveHistogram(h_vx_vs_tht, outPath, "h_vx_vs_tht_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_tht_vx, outPath, "h_tht_vx_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_vs_tht, outPath, "h_vy_vs_tht_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_tht_vy, outPath, "h_tht_vy_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_vs_tht, outPath, "h_vz_vs_tht_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_tht_vz, outPath, "h_tht_vz_kalman", saveoutput, close);

	if(close) exit(0);

}
