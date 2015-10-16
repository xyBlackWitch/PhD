/**
 * @file vertex_studies_lambda0.C
 * @author Jennifer Puetz (j.puetz@fz-juelich.de)
 * @date 2015
 * @brief Creates histograms for vertex studies of Xi- decay.
 * @details Methode creates Histogramms for vertex studies of a Xi- decay. Data is comming from the analysis which is stored in a file called output_ana.root
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
#include "TCanvas.h"
#include "TH2.h"
#include "../common_jenny.cpp"


void vertex_studies_XiMinus_kalman(int nevts=0, bool saveoutput=true, bool close=false){

	//Paths

	TString inPath = "/home/ikp1/puetz/panda/mysimulations/test/boxgenerator/Xi/500000_events_beammom_3/";
	TString outPath = inPath + "plots/";


	//Input File

	TString inputFile = inPath + "ana_total.root";

	//open input
	TFile * file = new TFile(inputFile, "READ");

	TTree * ntpXiMinus = (TTree*) file->Get("ntpXiMinus");
//	TTree * ntpMC = (TTree*) file->Get("ntpMC");

	//****Projection for XiMinus

//	// simulated vertex position for XiMinus
//
//	TH1D* h_vx_mc = new TH1D("h_vx_mc", "x coordinate of the vertex position; x[cm]; counts", 1000,-5,5);
//	ntpMC->Project("h_vx_mc", "x", "moth==0 && pdg==3312");
//
//	TH1D* h_vy_mc = new TH1D("h_vy_mc", "y coordinate of the vertex position; y[cm]; counts", 1000,-5,5);
//	ntpMC->Project("h_vy_mc", "y", "moth==0 && pdg==3312");
//
//	TH1D* h_vz_mc = new TH1D("h_vz_mc", "z coordinate of the vertex position; z[cm]; counts", 1000,-5,5);
//	ntpMC->Project("h_vz_mc", "z", "moth==0 && pdg==3312");



	//reconstructed vertex position for XiMinus
	TH1D* h_vx_reco = new TH1D("h_vx_reco", "x coordinate of the vertex position; x[cm]; counts", 1000,-5,5);
	ntpXiMinus->Project("h_vx_reco", "KalmanFit_vx", "McTruthMatch==1");

	TH1D* h_vy_reco = new TH1D("h_vy_reco", "y coordinate of the vertex position; y[cm]; counts", 1000,-5,5);
	ntpXiMinus->Project("h_vy_reco", "KalmanFit_vy", "McTruthMatch==1");

	TH1D* h_vz_reco = new TH1D("h_vz_reco", "z coordinate of the vertex position; z[cm]; counts", 1000,-5,5);
	ntpXiMinus->Project("h_vz_reco", "KalmanFit_vz", "McTruthMatch==1 ");



	//resolution of reconstructed vertex position

	TH1D* h_vx_res = new TH1D("h_vx_res", "x resolution of the vertex position; #Delta_{x}[cm]; counts", 1000,-1,1);
	ntpXiMinus->Project("h_vx_res", "KalmanFit_vx-MCTruth_vx", "McTruthMatch==1 ");

	TH1D* h_vy_res = new TH1D("h_vy_res", "y resolution of the vertex position; #Delta_{y}[cm]; counts", 1000,-1,1);
	ntpXiMinus->Project("h_vy_res", "KalmanFit_vy-MCTruth_vy", "McTruthMatch==1 ");

	TH1D* h_vz_res = new TH1D("h_vz_res", "z resolution of the vertex position; #Delta_{z}[cm]; counts", 1000,-1,1);
	ntpXiMinus->Project("h_vz_res", "KalmanFit_vz-MCTruth_vz", "McTruthMatch==1");


	//resolution of reconstructed vertex position and the poca

	TH1D* h_vx_poc_res = new TH1D("h_vx_poc_res", "differenz reconstructed vertex position and poca position; x_{reco} - x_{poca}[cm]; counts", 1000,-1,1);
	ntpXiMinus->Project("h_vx_poc_res", "KalmanFit_vx - XiMinus_pocvx", "McTruthMatch==1 ");

	TH1D* h_vy_poc_res = new TH1D("h_vy_poc_res", "differenz reconstructed vertex position and poca position; y_{reco} - y_{poca}[cm]; counts", 1000,-1,1);
	ntpXiMinus->Project("h_vy_poc_res", "KalmanFit_vy - XiMinus_pocvy", "McTruthMatch==1 ");

	TH1D* h_vz_poc_res = new TH1D("h_vz_poc_res", "differenz reconstructed vertex position and poca position; z_{reco} - z_{poca}[cm]; counts", 1000,-1.5,1.5);
	ntpXiMinus->Project("h_vz_poc_res", "KalmanFit_vz - XiMinus_pocvz", "McTruthMatch==1");




//	//ratio plot of vertex postion for XiMinus
//	TH1D* h_reco_div_mc_x = new TH1D("h_reco_div_mc_x", "ratio reco/MC for x position of vertex; x[cm]; ratio reco/MC", 1000,-5,5);
//	h_reco_div_mc_x->Divide(h_vx_reco, h_vx_mc);
//	h_reco_div_mc_x->SetStats(0);
//
//	TH1D* h_reco_div_mc_y = new TH1D("h_reco_div_mc_y", "ratio reco/MC for y position of vertex; y[cm]; ratio reco/MC", 1000,-5,5);
//	h_reco_div_mc_y->Divide(h_vy_reco, h_vy_mc);
//	h_reco_div_mc_y->SetStats(0);
//
//	TH1D* h_reco_div_mc_z = new TH1D("h_reco_div_mc_z", "ratio reco/MC for z position of vertex; z[cm]; ratio reco/MC", 1000,-5,5);
//	h_reco_div_mc_z->Divide(h_vz_reco, h_vz_mc);
//	h_reco_div_mc_z->SetStats(0);



	//2-dimensional plots vertex

	TH2D * h_vxy_vz_reco = new TH2D("h_vxy_vz_reco", "Vertex ; z[cm]; xy[cm]", 200,0,5,200,0,5);
	ntpXiMinus->Project("h_vxy_vz_reco", "sqrt(KalmanFit_vx**2+KalmanFit_vy**2):KalmanFit_vz", "McTruthMatch==1 ");



	//Goodness of fit

	TH1D * h_chi2 = new TH1D("h_chi2", "#chi^{2} distribution for vertex Fit; #chi^{2}; counts", 500,0,10);
	ntpXiMinus->Project("h_chi2", "KalmanFit_chisq", "McTruthMatch==1");

	TH1D * h_Prob = new TH1D("h_Prob", "Probability for vertex Fit; Prob; counts", 500,0,1);
	ntpXiMinus->Project("h_Prob", "KalmanFit_prob", "McTruthMatch==1");


	//**********Performance plots******************
	//pulls for vertex

	TH1D * h_pull_vx = new TH1D("h_pull_vx", "pull for x coordinate of momentum; (px_{reco}-px_{MC})/#sigma_{px}; counts", 500, -5,5);
	ntpXiMinus->Project("h_pull_vx", "KalmanFit_mcpullpx", "McTruthMatch==1");

	TH1D * h_pull_vy = new TH1D("h_pull_vy", "pull for y coordinate of momentum; (py_{reco}-py_{MC})/#sigma_{py}; counts", 500, -5,5);
	ntpXiMinus->Project("h_pull_vy", "KalmanFit_mcpullpy", "McTruthMatch==1");

	TH1D * h_pull_vz = new TH1D("h_pull_vz", "pull for z coordinate of momentum; (pz_{reco}-pz_{MC})/#sigma_{pz}; counts", 500, -5,5);
	ntpXiMinus->Project("h_pull_vz", "KalmanFit_mcpullpz", "McTruthMatch==1");


	//vtx res vs p_t

	TH2D * h_vx_vs_pt = new TH2D("h_vx_vs_pt", " transversal momentum vs. vertex resolution(x coordinate); #Delta x/cm; pt/GeV", 200, -1,1, 200, 0,2);
	ntpXiMinus->Project("h_vx_vs_pt", "KalmanFit_pt:(KalmanFit_vx-MCTruth_vx)", "McTruthMatch==1");
	TH1D * h_pt_vx = h_vx_vs_pt->ProfileY("h_pt_vx", 1, -1, "o");
	h_pt_vx->SetTitle("vertex resolution(x coordinate) vs. transversam momentum");
	h_pt_vx->GetYaxis()->SetTitle("#Delta x/cm");
	h_pt_vx->GetYaxis()->SetRangeUser(-0.8,0.8);

	TH2D * h_vy_vs_pt = new TH2D("h_vy_vs_pt", "transversal momentum vs. vertex resolution (y coordinate); #Delta y/cm; p_{t}/GeV", 200, -1,1, 200, 0,2);
	ntpXiMinus->Project("h_vy_vs_pt", "KalmanFit_pt:(KalmanFit_vy-MCTruth_vy)", "McTruthMatch==1");
	TH1D * h_pt_vy = h_vy_vs_pt->ProfileY("h_pt_vy", 1, -1, "o");
	h_pt_vy->SetTitle("vertex resolution(y coordinate) vs. transversam momentum");
	h_pt_vy->GetYaxis()->SetTitle("#Delta y/cm");
	h_pt_vy->GetYaxis()->SetRangeUser(-0.8,0.8);

	TH2D * h_vz_vs_pt = new TH2D("h_vz_vs_pt", "transversal momentum vs. vertex resolution (z coordinate); #Delta z/cm; p_{t}/GeV", 200, -1,1, 200, 0,2);
	ntpXiMinus->Project("h_vz_vs_pt", "KalmanFit_pt:(KalmanFit_vz-MCTruth_vz)", "McTruthMatch==1");
	TH1D * h_pt_vz = h_vz_vs_pt->ProfileY("h_pt_vz", 1, -1, "o");
	h_pt_vz->SetTitle("vertex resolution(z coordinate) vs. transversam momentum");
	h_pt_vz->GetYaxis()->SetTitle("#Delta z/cm");
	h_pt_vz->GetYaxis()->SetRangeUser(-0.8,0.8);

	//vtx res vs p_z

	TH2D * h_vx_vs_pz = new TH2D("h_vx_vs_pz", "longitudinal momentum vs. vertex resolution(x coordinate); #Delta x/cm; p_{z}/GeV", 200, -1,1, 200, 2,3);
	ntpXiMinus->Project("h_vx_vs_pz", "KalmanFit_pz:(KalmanFit_vx-MCTruth_vx)", "McTruthMatch==1");
	TH2D * h_vy_vs_pz = new TH2D("h_vy_vs_pz", "longitudinal momentum vs. vertex resolution(y coordinate); #Delta y/cm; p_{z}/GeV", 200, -1,1, 200, 2,3);
	ntpXiMinus->Project("h_vy_vs_pz", "KalmanFit_pz:(KalmanFit_vy-MCTruth_vy)", "McTruthMatch==1");
	TH2D * h_vz_vs_pz = new TH2D("h_vz_vs_pz", "longitudinal momentum vs. vertex resolution(z coordinate); #Delta z/cm; p_{z}/GeV", 200, -1,1, 200, 2,3);
	ntpXiMinus->Project("h_vz_vs_pz", "KalmanFit_pz:(KalmanFit_vz-MCTruth_vz)", "McTruthMatch==1");

	//vtx res vs p_z

	TH2D * h_vx_vs_p = new TH2D("h_vx_vs_p", " momentum vs. vertex resolution(x coordinate); #Delta x/cm; p/GeV", 200, -1,1, 200, 0,5);
	ntpXiMinus->Project("h_vx_vs_p", "KalmanFit_p:(KalmanFit_vx-MCTruth_vx)", "McTruthMatch==1");
	TH2D * h_vy_vs_p = new TH2D("h_vy_vs_p", " momentum vs. vertex resolution(y coordinate); #Delta y/cm; p/GeV", 200, -1,1, 200, 0,5);
	ntpXiMinus->Project("h_vy_vs_p", "KalmanFit_p:(KalmanFit_vy-MCTruth_vy)", "McTruthMatch==1");
	TH2D * h_vz_vs_p = new TH2D("h_vz_vs_p", " momentum vs. vertex resolution(z coordinate); #Delta z/cm; p/GeV", 200, -1,1, 200, 0,5);
	ntpXiMinus->Project("h_vz_vs_p", "KalmanFit_p:(KalmanFit_vz-MCTruth_vz)", "McTruthMatch==1");


	//vtx res vs tht

	TH2D * h_vx_vs_tht = new TH2D("h_vx_vs_tht", " Theta vs. vertex resolution(x coordinate); #Delta x/cm; tht/rad", 200, -1,1, 200, 0,0.6);
	h_vx_vs_tht->GetZaxis()->SetRangeUser(0,100);
	ntpXiMinus->Project("h_vx_vs_tht", "KalmanFit_tht:(KalmanFit_vx-MCTruth_vx)", "McTruthMatch==1");
	TH1D * h_tht_vx = h_vx_vs_tht->ProfileY("h_tht_vx", 1, -1, "o");
	h_tht_vx->SetTitle("vertex resolution(x coordinate) vs. Theta");
	h_tht_vx->GetYaxis()->SetTitle("#Delta x/cm");

	TH2D * h_vy_vs_tht = new TH2D("h_vy_vs_tht", " Theta vs. vertex resolution(y coordinate); #Delta y/cm; tht/rad", 200, -1,1, 200, 0,0.6);
	ntpXiMinus->Project("h_vy_vs_tht", "KalmanFit_tht:(KalmanFit_vy-MCTruth_vy)", "McTruthMatch==1");
	TH1D * h_tht_vy = h_vy_vs_tht->ProfileY("h_tht_vy", 1, -1, "o");
	h_tht_vy->SetTitle("vertex resolution(y coordinate) vs. Theta");
	h_tht_vy->GetYaxis()->SetTitle("#Delta y/cm");

	TH2D * h_vz_vs_tht = new TH2D("h_vz_vs_tht", " Theta vs. vertex resolution(z coordinate); #Delta z/cm; tht/rad", 200, -1,1, 200, 0,0.6);
	ntpXiMinus->Project("h_vz_vs_tht", "KalmanFit_tht:(KalmanFit_vz-MCTruth_vz)", "McTruthMatch==1");
	TH1D * h_tht_vz = h_vz_vs_tht->ProfileY("h_tht_vz", 1, -1, "o");
	h_tht_vz->SetTitle("vertex resolution(z coordinate) vs. Theta");
	h_tht_vz->GetYaxis()->SetTitle("#Delta z/cm");



	//****Create Canvas, draw histogram and save it
//
//	jenny::CreateDrawAndSaveHistogram(h_vx_mc, outPath, "h_vx_mc_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_mc, outPath, "h_vy_mc_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_mc, outPath, "h_vz_mc_kalman", saveoutput, close);

//	jenny::CreateDrawAndSaveHistogram(h_vx_reco, outPath, "h_vx_reco_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_reco, outPath, "h_vy_reco_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_reco, outPath, "h_vz_reco_kalman", saveoutput, close);

	bool noAutorange =false;

//	jenny::CreateDrawAndSaveHistogramWithFit(h_vx_res, outPath, "h_vx_res_fit_kalman", saveoutput, close, noAutorange, 0.06, 0.5,true);
//	jenny::CreateDrawAndSaveHistogramWithFit(h_vy_res, outPath, "h_vy_res_fit_kalman", saveoutput, close, noAutorange, 0.07, 1, true);
//	jenny::CreateDrawAndSaveHistogramWithFit(h_vz_res, outPath, "h_vz_res_fit_kalman", saveoutput, close, noAutorange, 0.06, 1,true);
//
//
//	jenny::CreateDrawAndSaveHistogram(h_vx_res, outPath, "h_vx_res_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_res, outPath, "h_vy_res_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_res, outPath, "h_vz_res_kalman", saveoutput, close);
//
//	jenny::CreateDrawAndSaveHistogramWithFit(h_vx_poc_res, outPath, "h_vx_poc_res_kalman", saveoutput, close, noAutorange, 0.04, 1.5);
//	jenny::CreateDrawAndSaveHistogramWithFit(h_vy_poc_res, outPath, "h_vy_poc_res_kalman", saveoutput, close, noAutorange, 0.04, 1.5);
//	jenny::CreateDrawAndSaveHistogramWithFit(h_vz_poc_res, outPath, "h_vz_poc_res_kalman", saveoutput, close, noAutorange, 0.04, 1.5);
//
//	jenny::CreateDrawAndSaveHistogram(h_reco_div_mc_x, outPath, "h_reco_div_mc_x_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_reco_div_mc_y, outPath, "h_reco_div_mc_y_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_reco_div_mc_z, outPath, "h_reco_div_mc_z_kalman", saveoutput, close);
//
//
//	jenny::CreateDrawAndSaveHistogram(h_vxy_vz_reco, outPath, "h_vxy_vz_reco_kalman", saveoutput, close);
//
//	jenny::CreateDrawAndSaveHistogram(h_chi2, outPath, "h_chi2_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_Prob, outPath, "h_Prob_kalman", saveoutput, close);
//
	jenny::CreateDrawAndSaveHistogram(h_pull_vx, outPath, "h_pull_vx_kalman", saveoutput, close);//, noAutorange, 0.4, 3);
	jenny::CreateDrawAndSaveHistogram(h_pull_vy, outPath, "h_pull_vy_kalman", saveoutput, close);//, noAutorange, 0.25, 3);
	jenny::CreateDrawAndSaveHistogram(h_pull_vz, outPath, "h_pull_vz_kalman", saveoutput, close);//, noAutorange, 0.4, 3);
//
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

//	jenny::CreateDrawAndSaveHistogram(h_vx_vs_tht, outPath, "h_vx_vs_tht_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_tht_vx, outPath, "h_tht_vx_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_vs_tht, outPath, "h_vy_vs_tht_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_tht_vy, outPath, "h_tht_vy_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_vs_tht, outPath, "h_vz_vs_tht_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_tht_vz, outPath, "h_tht_vz_kalman", saveoutput, close);

	if(close) exit(0);

}
