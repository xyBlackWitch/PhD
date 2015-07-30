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
#include "common_jenny.cpp"





void vertex_studies_lambda0_kalman(int nevts=0, bool saveoutput=true, bool close=false){


	//Paths

	TString inPath = "/private/puetz/mysimulations/test/boxgenerator/lambda0/10000_events/";
	TString outPath = inPath +"plots/";


	//Input File

	TString inputFile = inPath + "output_ana_kalman.root";

	//open input
	TFile * file = new TFile(inputFile, "READ");

	TTree * ntpLambda = (TTree*) file->Get("ntpLambda0");
	TTree * ntpMC = (TTree*) file->Get("ntpMC");


	//****Projection for Lambda0


	//simulated vertex position for lambda0

	TH1D * h_vx_mc = new TH1D("h_vx_mc", "x coodinate of vertex position; x[cm]; counts", 100,-10,10);
	ntpMC->Project("h_vx_mc", "x", "moth==-1");

	TH1D * h_vy_mc = new TH1D("h_vy_mc", "y coodinate of vertex position; y[cm]; counts", 100,-10,10);
	ntpMC->Project("h_vy_mc", "y", "moth==-1");

	TH1D * h_vz_mc = new TH1D("h_vz_mc", "z coodinate of vertex position; z[cm]; counts", 100,-10,10);
	ntpMC->Project("h_vz_mc", "z", "moth==-1");

	//Truth information

	TH1D * h_vx_truth = new TH1D("h_vx_truth", "x coodinate of vertex position; x[cm]; counts", 100,-100,100);
	ntpLambda0->Project("h_vx_truth", "McTruth_vx", "McTruthMatch==1");

	TH1D * h_vy_truth = new TH1D("h_vy_truth", "y coodinate of vertex position; y[cm]; counts", 100,-100,100);
	ntpLambda0->Project("h_vy_truth", "McTruth_vy", "McTruthMatch==1");

	TH1D * h_vz_truth = new TH1D("h_vz_truth", "z coodinate of vertex position; z[cm]; counts", 100,-100,100);
	ntpLambda0->Project("h_vz_truth", "McTruth_vz", "McTruthMatch==1");


	//reconstructed vertex position for lambda0
	TH1D * h_vx_reco = new TH1D("h_vx_reco", "x coordinate of reconstructed vertex position; x[cm]; counts", 100,-100,100);
	ntpLambda->Project("h_vx_reco", "KalmanFit_vx", "McTruthMatch==1");

	TH1D * h_vy_reco = new TH1D("h_vy_reco", "y coordinate of reconstructed vertex position; y[cm]; counts", 100,-100,100);
	ntpLambda->Project("h_vy_reco", "KalmanFit_vy", "McTruthMatch==1");

	TH1D * h_vz_reco = new TH1D("h_vz_reco", "z coordinate of reconstructed vertex position; z[cm]; counts", 100,-100,100);
	ntpLambda->Project("h_vz_reco", "KalmanFit_vz", "McTruthMatch==1");



	//resolution for reconstructed vertex position for lambda0

	TH1D * h_vx_res = new TH1D("h_vx_res", "x resolution of reconstructed vertex position; #Delta x[cm]; counts", 100,-1,1);
	ntpLambda->Project("h_vx_res", "KalmanFit_vx-McTruth_vx", "McTruthMatch==1");

	TH1D * h_vy_res = new TH1D("h_vy_res", "y resolution of reconstructed vertex position; #Delta y[cm]; counts", 100,-1,1);
	ntpLambda->Project("h_vy_res", "KalmanFit_vy-McTruth_vy", "McTruthMatch==1");

	TH1D * h_vz_res = new TH1D("h_vz_res", "z resolution of reconstructed vertex position; #Delta z[cm]; counts", 100,-1,1);
	ntpLambda->Project("h_vz_res", "KalmanFit_vz-McTruth_vz", "McTruthMatch==1");


	//resolution for reconstructed vertex position and poca for lambda0

	TH1D * h_vx_poc_res = new TH1D("h_vx_poc_res", "differenz reconstructed vertex position and poca position; x_{reco} - x_{poca}[cm]; counts", 100,-1,1);
	ntpLambda->Project("h_vx_poc_res", "KalmanFit_vx-Lambda0_pocvx", "McTruthMatch==1");

	TH1D * h_vy_poc_res = new TH1D("h_vy_poc_res", "differenz reconstructed vertex position and poca position; y_{reco} - y_{poca}[cm]; counts", 100,-1,1);
	ntpLambda->Project("h_vy_poc_res", "KalmanFit_vy-Lambda0_pocvy", "McTruthMatch==1");

	TH1D * h_vz_poc_res = new TH1D("h_vz_poc_res", "differenz reconstructed vertex position and poca position; z_{reco} - z_{poca}[cm]; counts", 100,-1,1);
	ntpLambda->Project("h_vz_poc_res", "KalmanFit_vz-Lambda0_pocvz", "McTruthMatch==1");


	//Combined histograms

	TH1D * h_reco_div_mc_vx = new TH1D("h_reco_div_mc_vx", "Ratio Reco/MC for the x coordinate of vertex position; x_{reco}/x_{MC}", 100,-100,100);
	h_reco_div_mc_vx->Divide(h_vx_reco,h_vx_truth);
	h_reco_div_mc_vx->SetStats(0);
//	h_reco_div_mc_vx->GetYaxis()->SetRangeUser(0,1.1);

	TH1D * h_reco_div_mc_vy = new TH1D("h_reco_div_mc_vy", "Ratio Reco/MC for the y coordinate of vertex position; y_{reco}/y_{MC}", 100,-100,100);
	h_reco_div_mc_vy->Divide(h_vy_reco,h_vy_truth);
	h_reco_div_mc_vy->SetStats(0);
//	h_reco_div_mc_vy->GetYaxis()->SetRangeUser(0,1.1);

	TH1D * h_reco_div_mc_vz = new TH1D("h_reco_div_mc_vz", "Ratio Reco/MC for the z coordinate of vertex position; z_{reco}/z_{MC}", 100,0,100);
	h_reco_div_mc_vz->Divide(h_vz_reco,h_vz_truth);
	h_reco_div_mc_vz->SetStats(0);
//	h_reco_div_mc_vz->GetYaxis()->SetRangeUser(0,1.1);


	//2-dimensional plots for poca and vertex

	TH2D * h_vxy_vz_reco = new TH2D("h_vxy_vz_reco", "Vertex ; z[cm]; r[cm]", 100,0,5,100,0,5);
	ntpLambda->Project("h_vxy_vz_reco", "sqrt(KalmanFit_vx**2+KalmanFit_vy**2):KalmanFit_vz", "McTruthMatch==1");


	//Goodness of fit

	TH1D * h_chi2 = new TH1D("h_chi2", "#chi^{2} distribution for vertex Fit; #chi^{2}; counts", 100,0,10);
	ntpLambda->Project("h_chi2", "KalmanFit_chisq", "McTruthMatch==1");

	TH1D * h_Prob = new TH1D("h_Prob", "Probability for vertex Fit; Prob; counts", 100,0,1);
	ntpLambda->Project("h_Prob", "KalmanFit_prob", "McTruthMatch==1");


	//***Create Canvas and draw histogramm
//
//	jenny::CreateDrawAndSaveHistogram(h_vx_mc, outPath, "h_vx_mc_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_mc, outPath, "h_vy_mc_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_mc, outPath, "h_vz_mc_kalman", saveoutput, close);

	bool autoRange = false;

	jenny::CreateDrawAndSaveHistogramWithFit(h_vx_res, outPath, "h_vx_res_kalman", saveoutput, close, autoRange, 0.1,1);
	jenny::CreateDrawAndSaveHistogramWithFit(h_vy_res, outPath, "h_vy_res_kalman", saveoutput, close, autoRange, 0.03,1);
	jenny::CreateDrawAndSaveHistogramWithFit(h_vz_res, outPath, "h_vz_res_kalman", saveoutput, close, autoRange, 0.06, 1);

//	jenny::CreateDrawAndSaveHistogram(h_vx_poc_res, outPath, "h_vx_poc_res_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_poc_res, outPath, "h_vy_poc_res_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_poc_res, outPath, "h_vz_poc_res_kalman", saveoutput, close);
//
//	jenny::CreateDrawAndSaveHistogram(h_vx_reco, outPath, "h_vx_reco_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vy_reco, outPath, "h_vy_reco_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_vz_reco, outPath, "h_vz_reco_kalman", saveoutput, close);
//
//	jenny::CreateDrawAndSaveHistogram(h_reco_div_mc_vx, outPath, "h_reco_div_mc_vx_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_reco_div_mc_vy, outPath, "h_reco_div_mc_vy_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_reco_div_mc_vz, outPath, "h_reco_div_mc_vz_kalman", saveoutput, close);
//
//
//	jenny::CreateDrawAndSaveHistogram(h_vxy_vz_reco, outPath, "h_vxy_vz_reco_kalman", saveoutput, close);
//
//
//	jenny::CreateDrawAndSaveHistogram(h_chi2, outPath, "h_chi2_kalman", saveoutput, close);
//	jenny::CreateDrawAndSaveHistogram(h_Prob, outPath, "h_Prob_kalman", saveoutput, close);


	if(close) exit(0);

}
